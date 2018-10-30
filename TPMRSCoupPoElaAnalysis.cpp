//
//  TPMRSCoupPoElaAnalysis.cpp
//  PZ
//
//  Created by Omar and Manouchehr on 8/28/16.
//
//

#include "TPMRSCoupPoElaAnalysis.h"
#include "pzpostprocanalysis.h"
#include "pzfstrmatrix.h"
#include "TPZElasticCriterion.h"

/// Brief default costructor
TPMRSCoupPoElaAnalysis::TPMRSCoupPoElaAnalysis() : TPZAnalysis()
{
    
    /// Brief define the simulation data
    m_SimulationData = NULL;
    
    /// Brief Vector of compmesh pointers. m_meshvec[0] = flowHdiv, m_meshvec[1] = PressureL2
    m_meshvec.Resize(2);
    
    /// Brief Part of residue at n+1 (current) state
    m_R_n.Resize(0,0);
    
    /// Brief Part of residue at n (last) state
    m_R.Resize(0,0);
    
    /// Brief Solution at n+1 (current) state
    m_X_n.Resize(0,0);
    
    /// Brief memory at n+1 state
    m_memory_n.Resize(0);
    
    /// Brief Solution  at n (last) state
    m_X.Resize(0,0);
    
    /// Brief memory at n (past) state
    m_memory.Resize(0);
    
    /// Brief Strain-Stress solution data
    m_strain_stress_duplets.Resize(0);
    
    /// Brief Residue error
    m_error    = 1.0;
    
    /// Brief Correction variation
    m_dx_norm  = 1.0;
    
    /// Brief number of newton corrections
    m_k_iterations = 0;
    
}

/// Brief default destructor
TPMRSCoupPoElaAnalysis::~TPMRSCoupPoElaAnalysis(){
    
}

/// Brief copy constructor
TPMRSCoupPoElaAnalysis::TPMRSCoupPoElaAnalysis(const TPMRSCoupPoElaAnalysis &copy)
{
    m_SimulationData = copy.m_SimulationData;
    m_meshvec        = copy.m_meshvec;
    m_R_n            = copy.m_R_n;
    m_R              = copy.m_R;
    m_X_n            = copy.m_X_n;
    m_X              = copy.m_X;
    m_error          = copy.m_error;
    m_dx_norm        = copy.m_dx_norm;
    
}

/// Brief Copy assignemnt operator
TPMRSCoupPoElaAnalysis & TPMRSCoupPoElaAnalysis::operator=(const TPMRSCoupPoElaAnalysis &other)
{
    if (this != & other) {  /// prevent self-assignment
        
        m_SimulationData = other.m_SimulationData;
        m_meshvec        = other.m_meshvec;
        m_R_n            = other.m_R_n;
        m_R              = other.m_R;
        m_X_n            = other.m_X_n;
        m_X              = other.m_X;
        m_error          = other.m_error;
        m_dx_norm        = other.m_dx_norm;
    }
    return *this;
}

/// Brief Resize and fill residue and solution vectors
void TPMRSCoupPoElaAnalysis::AdjustVectors()
{
    
    if(fSolution.Rows() == 0 /* || fRhs.Rows() == 0 */)
    {
        DebugStop();
    }
    
    TPZBuildMultiphysicsMesh::AddElements(m_meshvec, this->Mesh());
    TPZBuildMultiphysicsMesh::AddConnects(m_meshvec, this->Mesh());
    TPZBuildMultiphysicsMesh::TransferFromMeshes(m_meshvec, this->Mesh());
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(m_meshvec, this->Mesh());
    
    m_X.Resize(fSolution.Rows(),1);
    m_X.Zero();
    m_X_n.Resize(fSolution.Rows(),1);
    m_X_n.Zero();
    m_R_n.Resize(fSolution.Rows(),1);
    m_R_n.Zero();
    m_R.Resize(fSolution.Rows(),1);
    m_R.Zero();
}

void TPMRSCoupPoElaAnalysis::QuasiNewtonIteration()
{
    
    if(m_k_iterations == 1)
    {
        this->Assemble();
    }
    else
    {
        this->AssembleResidual();
    }
    this->Rhs() += m_R; // total residue
    this->Rhs() *= -1.0;

#ifdef PZDEBUG
//        this->Solver().Matrix()->Print("J = ", std::cout,EMathematicaInput);
//        this->Rhs().Print("R = ", std::cout,EMathematicaInput);
#endif
    
    this->Solve(); // update correction
    m_dx_norm = Norm(this->Solution()); // correction variation
    
#ifdef PZDEBUG
//    this->Solution().Print("dx = ", std::cout,EMathematicaInput);
#endif
    
    m_X_n += this->Solution(); // update state
    
    
    /// Check the update state at current state (n+1) for PMRS_PoroElastic and PMRS_PoroPlastic
  
    this->Update_at_n_State();

    
    this->AssembleResidual();
    m_R_n = this->Rhs();
    
#ifdef PZDEBUG
//    m_X.Print("X = ", std::cout,EMathematicaInput);
//    m_X_n.Print("Xn = ", std::cout,EMathematicaInput);
//    m_R.Print("R = ", std::cout,EMathematicaInput);
//    m_R_n.Print("Rn = ", std::cout,EMathematicaInput);
#endif
    
    m_R_n += m_R; /// total residue
    
#ifdef PZDEBUG
//    m_R_n.Print("Rt = ", std::cout,EMathematicaInput);
#endif
    
    m_error =  Norm(m_R_n); /// residue error
    
}

void TPMRSCoupPoElaAnalysis::ExcecuteOneStep(){
    
    this->SimulationData()->SetCurrentStateQ(false);
    
    /// Check the update state at last state (n) for PMRS_PoroElastic and PMRS_PoroPlastic
    this->UpdateState();

    this->AssembleResidual();
    m_R = this->Rhs();
    
    
    this->SimulationData()->SetCurrentStateQ(true);
    
    /// Check the update state at current state (n+1) for PMRS_PoroElastic and PMRS_PoroPlastic

    this->Update_at_n_State();
    
    m_error = 1.0;
    
    STATE epsilon_res = this->SimulationData()->epsilon_res();
    STATE epsilon_cor = this->SimulationData()->epsilon_cor();
    int n_it  =   this->SimulationData()->n_iterations();
    
    m_SimulationData->Set_must_accept_solution_Q(true); /// For now acceting any solution in the party
    for (int k = 1; k <= n_it; k++)
    {
        this->Set_k_ietrarions(k);
        this->QuasiNewtonIteration();
        
        if(m_error < epsilon_res || (m_dx_norm < epsilon_cor && k > 3 ) )
        {
            std::cout << "PMRS:: Converged with iterations:  " << k << "; error: " << m_error <<  "; dx: " << m_dx_norm << std::endl;
            m_X = m_X_n;
            return;
        }
    }
    std::cout << "PMRS:: Exit max iterations with min dt:  " << m_SimulationData->dt() << "; (secs) " << "; error: " << m_error <<  "; dx: " << m_dx_norm << std::endl;
    
}


/// Brief update last state (at n state) solution for PMRS_PoroElastic
void TPMRSCoupPoElaAnalysis::UpdateState()
{
    this->LoadSolution(m_X);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(m_meshvec, this->Mesh());
}

/// Brief update current state (at n+1 state) solution for PMRS_PoroPlastic
void TPMRSCoupPoElaAnalysis::Update_at_n_State()
{
    this->LoadSolution(m_X_n);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(m_meshvec, this->Mesh());
}


/// Brief the Standard Post process function
void TPMRSCoupPoElaAnalysis::PostProcessStep()
{
    /// Post Process when you want to use datavec in your solution or you don't use memory and integration point
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(m_meshvec, this->Mesh());
    const int dim = this->Mesh()->Dimension();
    int div = m_SimulationData->n_div();
    
    TPZStack<std::string> outputcontrolsres;
    TPZStack<std::string> outputcontrolsgeo;
    
    /// Reservoir outpouts
    outputcontrolsres.Push("p");
    outputcontrolsres.Push("phi");
    outputcontrolsres.Push("kappa");
    
    /// Geomechanics outpouts
    outputcontrolsgeo.Push("ux");
    outputcontrolsgeo.Push("uy");
    outputcontrolsgeo.Push("sxx");
    outputcontrolsgeo.Push("syy");
    outputcontrolsgeo.Push("szz");
    outputcontrolsgeo.Push("exx");
    outputcontrolsgeo.Push("eyy");
    outputcontrolsgeo.Push("ezz");
    outputcontrolsgeo.Push("epxx");
    outputcontrolsgeo.Push("epyy");
    outputcontrolsgeo.Push("epzz");
    
    
    std::string plotfile = m_SimulationData->name_vtk_file();
    
    this->DefineGraphMesh(dim,outputcontrolsres,outputcontrolsgeo,plotfile);
    this->PostProcess(div,dim);
    
    std::cout << "Standard post-processing finished." << std::endl;
    
}


/// Brief execute the evolutionary problem
void TPMRSCoupPoElaAnalysis::Run_Evolution(TPZVec<REAL> &x)
{
    int n = m_SimulationData->n_steps();
    REAL time = 0.0;
    REAL dt = this->SimulationData()->dt();

    for (int i = 0; i < n; i++)
    {
        this->ExcecuteOneStep();
        this->PostProcessStep();

//        this->AppendStrain_Stress(x);
//        this->AppendStrain_Pososity(x);
//        this->AppendStrain_Permeability(x);
//        this->AppendStrain_Pressure(x);
        
        time = (i+1)* dt;
        std::cout<< "PMRS:: Current time (s) = " << time << std::endl;
        this->SimulationData()->SetTime(time);
        
    }
}

/// Brief Compute the strain and the stress at x euclidean point for each time
void TPMRSCoupPoElaAnalysis::AppendStrain_Stress(TPZVec<REAL> & x)
{
    // Finding the geometic element that x bleongs to.
    REAL Tol = 1.0e-4;
    TPZGeoMesh * geometry = this->Mesh()->Reference();
    this->Mesh()->LoadReferences();
    
    int dim = geometry->Dimension();
    bool IsTargetElementQ = false;
    long n_elemenst = geometry->NElements();
    
    int sx_var = 10;
    int sy_var = 13;
    int ey_var = 19;

    TPZVec<STATE> sx;
    TPZVec<STATE> sy;
    TPZVec<STATE> ey;
    
    std::pair<REAL,REAL> duplet;
    
    TPZVec<REAL> parametric_space(dim,0.0);
    
    
    for (long iel = 0; iel < n_elemenst; iel++)
    {
        TPZGeoEl * gel = geometry->Element(iel);
        
#ifdef PZDEBUG
        if (!gel)
        {
            DebugStop();
        }
#endif

        if (gel->Dimension() != dim) {
            continue;
        }
        
        IsTargetElementQ = gel->ComputeXInverse(x, parametric_space, Tol);
        
        if(IsTargetElementQ){
            TPZCompEl * cel = gel->Reference();
            cel->Solution(parametric_space, sx_var, sx);
            cel->Solution(parametric_space, sy_var, sy);
            cel->Solution(parametric_space, ey_var, ey);

            duplet.first = ey[0];
            duplet.second = sy[0] - sx[0];
            
            duplet.first = -duplet.first;
            duplet.second = fabs(duplet.second);
        }
    }
    m_strain_stress_duplets.Push(duplet);
}

/// Brief Compute the strain and the Pososity at x euclidean point for each time
void TPMRSCoupPoElaAnalysis::AppendStrain_Pososity(TPZVec<REAL> & x)
{
    /// Finding the geometic element that x bleongs to
    REAL Tol = 1.0e-4;
    TPZGeoMesh * geometry = this->Mesh()->Reference();
    this->Mesh()->LoadReferences();
    
    int dim = geometry->Dimension();
    bool IsTargetElementQ = false;
    long n_elemenst = geometry->NElements();
    
    int phi_var =  4;
    int ex_var  = 16;
    int ey_var  = 19;
    TPZVec<STATE> phi;
    TPZVec<STATE> ex;
    TPZVec<STATE> ey;
    
    std::pair<REAL,REAL> duplet;
    
    TPZVec<REAL> parametric_space(dim,0.0);
    for (long iel = 0; iel < n_elemenst; iel++)
    {
        TPZGeoEl * gel = geometry->Element(iel);
        
#ifdef PZDEBUG
        if (!gel)
        {
            DebugStop();
        }
#endif
        
        if (gel->Dimension() != dim)
        {
            continue;
        }
        
        
        IsTargetElementQ = gel->ComputeXInverse(x, parametric_space, Tol);
        
        if(IsTargetElementQ){
            TPZCompEl * cel = gel->Reference();
            cel->Solution(parametric_space, phi_var, phi);
            cel->Solution(parametric_space, ex_var, ex);
            cel->Solution(parametric_space, ey_var, ey);
            duplet.first = ex[0]+ ey[0];
            duplet.second = phi[0];
            
            duplet.first = fabs(duplet.first);
            duplet.second = fabs(duplet.second);
        }
    }
    
    m_strain_porosity_duplets.Push(duplet);
}

/// Brief Compute the strain and the Permeability at x euclidean point for each time
void TPMRSCoupPoElaAnalysis::AppendStrain_Permeability(TPZVec<REAL> & x)
{
    
    /// Finding the geometic element that x bleongs to
    REAL Tol = 1.0e-4;
    TPZGeoMesh * geometry = this->Mesh()->Reference();
    this->Mesh()->LoadReferences();
    
    int dim = geometry->Dimension();
    bool IsTargetElementQ = false;
    long n_elemenst = geometry->NElements();
    
    int ky_var =  3;
    int ex_var = 16;
    int ey_var = 19;
    TPZVec<STATE> k;
    TPZVec<STATE> ex;
    TPZVec<STATE> ey;
    
    std::pair<REAL,REAL> duplet;
    
    TPZVec<REAL> parametric_space(dim,0.0);
    for (long iel = 0; iel < n_elemenst; iel++)
    {
        TPZGeoEl * gel = geometry->Element(iel);
        
#ifdef PZDEBUG
        if (!gel)
        {
            DebugStop();
        }
#endif
        
        if (gel->Dimension() != dim)
        {
            continue;
        }
        
        
        IsTargetElementQ = gel->ComputeXInverse(x, parametric_space, Tol);
        
        if(IsTargetElementQ){
            TPZCompEl * cel = gel->Reference();
            cel->Solution(parametric_space, ky_var, k);
            cel->Solution(parametric_space, ex_var, ex);
            cel->Solution(parametric_space, ey_var, ey);
            duplet.first = ex[0] + ey[0];
            duplet.second = k[0];
            
            duplet.first = fabs(duplet.first);
            duplet.second = fabs(duplet.second);
        }
    }
    
   m_strain_permeability_duplets.Push(duplet);
}

/// Brief Compute the strain and the Permeability at x euclidean point for each time
void TPMRSCoupPoElaAnalysis::AppendStrain_Pressure(TPZVec<REAL> & x)
{
    
    /// Finding the geometic element that x bleongs to
    REAL Tol = 1.0e-4;
    TPZGeoMesh * geometry = this->Mesh()->Reference();
    this->Mesh()->LoadReferences();
    
    int dim = geometry->Dimension();
    bool IsTargetElementQ = false;
    long n_elemenst = geometry->NElements();
    
    int p_var  =  0;
    int ex_var = 16;
    int ey_var = 19;
    TPZVec<STATE> p;
    TPZVec<STATE> ex;
    TPZVec<STATE> ey;
    
    std::pair<REAL,REAL> duplet;
    
    TPZVec<REAL> parametric_space(dim,0.0);
    for (long iel = 0; iel < n_elemenst; iel++)
    {
        TPZGeoEl * gel = geometry->Element(iel);
        
#ifdef PZDEBUG
        if (!gel)
        {
            DebugStop();
        }
#endif
        
        if (gel->Dimension() != dim)
        {
            continue;
        }
        
        IsTargetElementQ = gel->ComputeXInverse(x, parametric_space, Tol);
        
        if(IsTargetElementQ){
            TPZCompEl * cel = gel->Reference();
            cel->Solution(parametric_space, p_var, p);
            cel->Solution(parametric_space, ex_var, ex);
            cel->Solution(parametric_space, ey_var, ey);
            duplet.first = ex[0] + ey[0];
            duplet.second = p[0];
            
            duplet.first = fabs(duplet.first);
            duplet.second = fabs(duplet.second);
        }
    }
    
    m_strain_pressure_duplets.Push(duplet);
}

/// Brief Compute the strain and the stress at x euclidean point for each time
void TPMRSCoupPoElaAnalysis::PlotStrainStress(std::string file_name)
{
    
#ifdef PZDEBUG
    if (m_strain_stress_duplets.size() == 0)
    {
        DebugStop();
    }
#endif
    
    int n_data = m_strain_stress_duplets.size();
    TPZFMatrix<REAL> points(n_data,2,0.0);
    for(int i = 0; i < n_data; i++){
        points(i,0) = m_strain_stress_duplets[i].first;
        points(i,1) = m_strain_stress_duplets[i].second;
    }
    
    {
        std::ofstream out(file_name.c_str());
        points.Print("data = ",out,EMathematicaInput);
    }
}

/// Brief Compute the strain and the Porosity at x euclidean point for each time
void TPMRSCoupPoElaAnalysis::PlotStrainPorosity(std::string file_name)
{
    
#ifdef PZDEBUG
    if (m_strain_porosity_duplets.size() == 0)
    {
        DebugStop();
    }
#endif
    
    int n_data = m_strain_porosity_duplets.size();
    TPZFMatrix<REAL> points(n_data,2,0.0);
    for(int i = 0; i < n_data; i++)
    {
        points(i,0) = m_strain_porosity_duplets[i].first;
        points(i,1) = m_strain_porosity_duplets[i].second;
    }
    
    {
        std::ofstream out(file_name.c_str());
        points.Print("phi = ",out,EMathematicaInput);
    }
}

/// Brief Compute the strain and the Porosity at x euclidean point for each time
void TPMRSCoupPoElaAnalysis::PlotStrainPermeability(std::string file_name)
{
    
#ifdef PZDEBUG
    if (m_strain_permeability_duplets.size() == 0)
    {
        DebugStop();
    }
#endif
    
    int n_data = m_strain_permeability_duplets.size();
    TPZFMatrix<REAL> points(n_data,2,0.0);
    for(int i = 0; i < n_data; i++){
        points(i,0) = m_strain_permeability_duplets[i].first;
        points(i,1) = m_strain_permeability_duplets[i].second;
    }
    
    {
        std::ofstream out(file_name.c_str());
        points.Print("k = ",out,EMathematicaInput);
    }
    
}

/// Brief Compute the strain and the Porosity at x euclidean point for each time
void TPMRSCoupPoElaAnalysis::PlotStrainPressure(std::string file_name)
{
    
#ifdef PZDEBUG
    if (m_strain_pressure_duplets.size() == 0)
    {
        DebugStop();
    }
#endif
    
    int n_data = m_strain_pressure_duplets.size();
    TPZFMatrix<REAL> points(n_data,2,0.0);
    for(int i = 0; i < n_data; i++)
    {
        points(i,0) = m_strain_pressure_duplets[i].first;
        points(i,1) = m_strain_pressure_duplets[i].second;
    }
    
    {
        std::ofstream out(file_name.c_str());
        points.Print("p = ",out,EMathematicaInput);
    }
    
}
