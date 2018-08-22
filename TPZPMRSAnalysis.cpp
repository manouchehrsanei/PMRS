//
//  TPZPMRSAnalysis.cpp
//  PZ
//
//  Created by Omar and Manouchehr on 8/28/16.
//
//

#include "TPZPMRSAnalysis.h"
#include "pzpostprocanalysis.h"
#include "pzfstrmatrix.h"
#include "TPZPMRSCouplPoroPlast.h"
#include "TPZElasticCriterion.h"

/** @brief default costructor */
TPZPMRSAnalysis::TPZPMRSAnalysis() : TPZAnalysis()
{
    
    /** @brief define the simulation data */
    m_SimulationData = NULL;
    
    /** @brief Vector of compmesh pointers. m_meshvec[0] = flowHdiv, m_meshvec[1] = PressureL2 */
    m_meshvec.Resize(2);
    
    /** @brief Part of residue at n+1 (current) state */
    m_R_n.Resize(0,0);
    
    /** @brief Part of residue at n (last) state  */
    m_R.Resize(0,0);
    
    /** @brief Solution at n+1 (current) state */
    m_X_n.Resize(0,0);
    
    /** @brief memory at n+1 state */
    m_memory_n.Resize(0);
    
    /** @brief Solution  at n (last) state */
    m_X.Resize(0,0);
    
    /** @brief memory at n (past) state */
    m_memory.Resize(0);
    
    /** @brief Strain-Stress solution data */
    m_strain_stress_duplets.Resize(0);
    
    /** @brief Residue error */
    m_error    = 1.0;
    
    /** @brief Correction variation */
    m_dx_norm  = 1.0;
    
    /** @brief number of newton corrections */
    m_k_iterations = 0;
    
}

/** @brief default destructor $ */
TPZPMRSAnalysis::~TPZPMRSAnalysis(){
    
}

/** @brief copy constructor $ */
TPZPMRSAnalysis::TPZPMRSAnalysis(const TPZPMRSAnalysis &copy)
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

/** @brief Copy assignemnt operator $ */
TPZPMRSAnalysis & TPZPMRSAnalysis::operator=(const TPZPMRSAnalysis &other)
{
    if (this != & other) {  // prevent self-assignment
        
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

/** @brief Resize and fill residue and solution vectors */
void TPZPMRSAnalysis::AdjustVectors()
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

void TPZPMRSAnalysis::QuasiNewtonIteration()
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
    
    
    // Check the update state at current state (n+1) for PMRS_PoroElastic and PMRS_PoroPlastic
    if(IsPoroElastic)
    {
        this->Standard_Update_at_n_State();
    } else
    {
        this->Update_at_n_State();
    }
    
    this->AssembleResidual();
    m_R_n = this->Rhs();
    
#ifdef PZDEBUG
//    m_X.Print("X = ", std::cout,EMathematicaInput);
//    m_X_n.Print("Xn = ", std::cout,EMathematicaInput);
//    m_R.Print("R = ", std::cout,EMathematicaInput);
//    m_R_n.Print("Rn = ", std::cout,EMathematicaInput);
#endif
    
    m_R_n += m_R; // total residue
    
#ifdef PZDEBUG
//    m_R_n.Print("Rt = ", std::cout,EMathematicaInput);
#endif
    
    m_error =  Norm(m_R_n); // residue error
    
}

void TPZPMRSAnalysis::ExcecuteOneStep(){
    
    this->SimulationData()->SetCurrentStateQ(false);
    
    // Check the update state at last state (n) for PMRS_PoroElastic and PMRS_PoroPlastic
    if(IsPoroElastic)
    {
        this->Standard_UpdateState();
    } else
    {
        this->UpdateState();
    }

    this->AssembleResidual();
    m_R = this->Rhs();
    
    this->SimulationData()->SetCurrentStateQ(true);
    
    
    // Check the update state at current state (n+1) for PMRS_PoroElastic and PMRS_PoroPlastic
    if(IsPoroElastic)
    {
        this->Standard_Update_at_n_State();
    } else
    {
        this->Update_at_n_State();
    }
    
    m_error = 1.0;
    
    STATE epsilon_res = this->SimulationData()->epsilon_res();
    STATE epsilon_cor = this->SimulationData()->epsilon_cor();
    int n  =   this->SimulationData()->n_iterations();
    
    m_SimulationData->Set_must_accept_solution_Q(true); // For now acceting any solution in the party
    for (int k = 1; k <= n; k++)
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


/** @brief update last state (at n state) solution for PMRS_PoroElastic */
void TPZPMRSAnalysis::Standard_UpdateState()
{
    this->LoadSolution(m_X);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(m_meshvec, this->Mesh());
}

/** @brief update current state (at n+1 state) solution for PMRS_PoroPlastic */
void TPZPMRSAnalysis::Standard_Update_at_n_State()
{
    this->LoadSolution(m_X_n);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(m_meshvec, this->Mesh());
}

/** @brief update last state (at n state) solution for PMRS_PoroPlastic */
void TPZPMRSAnalysis::UpdateState()
{
    this->LoadSolution(m_X);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(m_meshvec, this->Mesh());
    
    int n_material = m_SimulationData->MaterialIds().size();
    if (n_material == 1) {
        int material_id =  m_SimulationData->MaterialIds()[0].first;
        TPZMaterial * material = this->Mesh()->FindMaterial(material_id);
        TPZPMRSCouplPoroPlast <TPZElasticCriterion, TPZPoroElastoPlasticMem> * rock_material = dynamic_cast<TPZPMRSCouplPoroPlast <TPZElasticCriterion, TPZPoroElastoPlasticMem> *>(material);
        
        if (!rock_material) { // There is no volumetric material of type TPZPMRSCouplPoroPlast <TPZElasticCriterion, TPZElastoPlasticMem>
            DebugStop();
        }
        
        SetMemory(rock_material->GetMemory());
    }
    else{
        //        Implement for several material ids
        DebugStop();
    }
    
}

/** @brief update current state (at n+1 state) solution for PMRS_PoroPlastic */
void TPZPMRSAnalysis::Update_at_n_State()
{
    this->LoadSolution(m_X_n);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(m_meshvec, this->Mesh());
    
    int n_material = m_SimulationData->MaterialIds().size();
    if (n_material == 1) {
        int material_id =  m_SimulationData->MaterialIds()[0].first;
        TPZMaterial * material = this->Mesh()->FindMaterial(material_id);
        TPZPMRSCouplPoroPlast <TPZElasticCriterion, TPZPoroElastoPlasticMem> * rock_material = dynamic_cast<TPZPMRSCouplPoroPlast <TPZElasticCriterion, TPZPoroElastoPlasticMem> *>(material);
        
        if (!rock_material) { // There is no volumetric material of type TPZPMRSCouplPoroPlast <TPZElasticCriterion, TPZElastoPlasticMem>
            DebugStop();
        }
        
        SetMemory_n(rock_material->GetMemory());
    }
    else{
        //        Implement for several material ids
        DebugStop();
    }
    
}


/** @brief the Standard Post process function */
void TPZPMRSAnalysis::PostProcessStepStandard()
{
    // * Post Process when you want to use datavec in your solution or you don't use memory and integration point *
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(m_meshvec, this->Mesh());
    const int dim = this->Mesh()->Dimension();
    int div = m_SimulationData->n_div();
    
    TPZManVector<std::string,50> scalnames = m_SimulationData->scalar_names();
    TPZManVector<std::string,50> vecnames = m_SimulationData->vector_names();
    
    std::string plotfile = m_SimulationData->name_vtk_file();
    
    this->DefineGraphMesh(dim,scalnames,vecnames,plotfile);
    this->PostProcess(div,dim);
    
    std::cout << "Standard post-processing finished." << std::endl;
    
}


/** @brief the Post process function */
void TPZPMRSAnalysis::PostProcessStep()
{
    // * Post Process when you want to use memory and integration point *
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(m_meshvec, this->Mesh());
    const int dim = this->Mesh()->Dimension();
    int div = m_SimulationData->n_div();
 
    TPZManVector<std::string,50> scalnames = m_SimulationData->scalar_names();
    TPZManVector<std::string,50> vecnames = m_SimulationData->vector_names();
    
    std::string plotfile = m_SimulationData->name_vtk_file();
    
    // Scalars and vectors in the integration point
    TPZStack<std::string> scalnames_intPoints;
    TPZStack<std::string> vecnames_intPoints;
    
    TPZStack<std::string> vars;
    for (auto varname : scalnames)
    {
            scalnames_intPoints.push_back(varname);
    }
    for (auto varname : vecnames)
    {
            vecnames_intPoints.push_back(varname);
    }
    
    TPZVec<int> PostProcMatIds(1);
    {
        PostProcMatIds[0] = 1;
        TPZPostProcAnalysis postProcessAnalysis;
        postProcessAnalysis.SetCompMesh(this->Mesh());
        TPZVec<std::string> vars(scalnames_intPoints.size()+vecnames_intPoints.size());
        for (unsigned int i = 0; i < scalnames_intPoints.size(); ++i)
        {
            vars[i] = scalnames_intPoints[i];
        }
        for (unsigned int i = 0; i < vecnames_intPoints.size(); ++i)
        {
            vars[scalnames_intPoints.size()+i] = vecnames_intPoints[i];
        }
        postProcessAnalysis.SetPostProcessVariables(PostProcMatIds, vars);
        TPZFStructMatrix structMatrix(postProcessAnalysis.Mesh());
        structMatrix.SetNumThreads(0);
        postProcessAnalysis.SetStructuralMatrix(structMatrix);
        postProcessAnalysis.TransferSolution();
        postProcessAnalysis.DefineGraphMesh(dim,scalnames_intPoints,vecnames_intPoints,plotfile);
        postProcessAnalysis.PostProcess(div, dim);
    }
    
    std::cout << "Post-processing finished." << std::endl;
    
}

/** @brief execute the evolutionary problem */
void TPZPMRSAnalysis::Run_Evolution(TPZVec<REAL> &x)
{
    int n = m_SimulationData->n_steps();
    REAL time = 0.0;
    REAL dt = this->SimulationData()->dt();

    for (int i = 0; i < n; i++)
    {
        this->ExcecuteOneStep();
        
        if(IsPoroElastic)
        {
            this->PostProcessStepStandard();
        } else
        {
        this->PostProcessStep();
        }
        
//        this->AppendStrain_Stress(x);
//        this->AppendStrain_Pososity(x);
//        this->AppendStrain_Permeability(x);
//        this->AppendStrain_Pressure(x);
        
        time = (i+1)* dt;
        std::cout<< "PMRS:: Current time (s) = " << time << std::endl;
        this->SimulationData()->SetTime(time);
        
    }
}

/** @brief Compute the strain and the stress at x euclidean point for each time */
void TPZPMRSAnalysis::AppendStrain_Stress(TPZVec<REAL> & x)
{
    // Finding the geometic element that x bleongs to.
    REAL Tol = 1.0e-4;
    TPZGeoMesh * geometry = this->Mesh()->Reference();
    this->Mesh()->LoadReferences();
    
    int dim = geometry->Dimension();
    bool IsTargetElementQ = false;
    long n_elemenst = geometry->NElements();
    
    int sx_var = 46;
    int sy_var = 47;
//    int eex_var = 26;
//    int epx_var = 33;
    int eey_var = 27;
    int epy_var = 34;
    TPZVec<STATE> sx;
    TPZVec<STATE> sy;
    TPZVec<STATE> eey;
    TPZVec<STATE> epy;
    
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
            cel->Solution(parametric_space, eey_var, eey);
            cel->Solution(parametric_space, epy_var, epy);
            duplet.first = eey[0] + epy[0];
            duplet.second = sy[0] - sx[0];
            
            duplet.first = -duplet.first;
            duplet.second = fabs(duplet.second);
        }
    }
    m_strain_stress_duplets.Push(duplet);
}

/** @brief Compute the strain and the Pososity at x euclidean point for each time */
void TPZPMRSAnalysis::AppendStrain_Pososity(TPZVec<REAL> & x)
{
    // Finding the geometic element that x bleongs to.
    REAL Tol = 1.0e-4;
    TPZGeoMesh * geometry = this->Mesh()->Reference();
    this->Mesh()->LoadReferences();
    
    int dim = geometry->Dimension();
    bool IsTargetElementQ = false;
    long n_elemenst = geometry->NElements();
    
    int phi_var = 42;
    int eex_var = 26;
    int epx_var = 33;
    int eey_var = 27;
    int epy_var = 34;
    TPZVec<STATE> phi;
    TPZVec<STATE> eex;
    TPZVec<STATE> epx;
    TPZVec<STATE> eey;
    TPZVec<STATE> epy;
    
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
            cel->Solution(parametric_space, eex_var, eex);
            cel->Solution(parametric_space, epx_var, epx);
            cel->Solution(parametric_space, eey_var, eey);
            cel->Solution(parametric_space, epy_var, epy);
            duplet.first = eex[0] + epx[0] + eey[0] + epy[0];
            duplet.second = phi[0];
            
            duplet.first = fabs(duplet.first);
            duplet.second = fabs(duplet.second);
        }
    }
    
    m_strain_porosity_duplets.Push(duplet);
}

/** @brief Compute the strain and the Permeability at x euclidean point for each time */
void TPZPMRSAnalysis::AppendStrain_Permeability(TPZVec<REAL> & x)
{
    
    // Finding the geometic element that x bleongs to.
    REAL Tol = 1.0e-4;
    TPZGeoMesh * geometry = this->Mesh()->Reference();
    this->Mesh()->LoadReferences();
    
    int dim = geometry->Dimension();
    bool IsTargetElementQ = false;
    long n_elemenst = geometry->NElements();
    
    int ky_var  = 44;
    int eex_var = 26;
    int epx_var = 33;
    int eey_var = 27;
    int epy_var = 34;
    TPZVec<STATE> k;
    TPZVec<STATE> eex;
    TPZVec<STATE> epx;
    TPZVec<STATE> eey;
    TPZVec<STATE> epy;
    
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
            cel->Solution(parametric_space, eex_var, eex);
            cel->Solution(parametric_space, epx_var, epx);
            cel->Solution(parametric_space, eey_var, eey);
            cel->Solution(parametric_space, epy_var, epy);
            duplet.first = eex[0] + epx[0] + eey[0] + epy[0];
            duplet.second = k[0];
            
            duplet.first = fabs(duplet.first);
            duplet.second = fabs(duplet.second);
        }
    }
    
   m_strain_permeability_duplets.Push(duplet);
}

/** @brief Compute the strain and the Permeability at x euclidean point for each time */
void TPZPMRSAnalysis::AppendStrain_Pressure(TPZVec<REAL> & x)
{
    
    // Finding the geometic element that x bleongs to.
    REAL Tol = 1.0e-4;
    TPZGeoMesh * geometry = this->Mesh()->Reference();
    this->Mesh()->LoadReferences();
    
    int dim = geometry->Dimension();
    bool IsTargetElementQ = false;
    long n_elemenst = geometry->NElements();
    
    int p_var   = 40;
    int eex_var = 26;
    int epx_var = 33;
    int eey_var = 27;
    int epy_var = 34;
    TPZVec<STATE> p;
    TPZVec<STATE> eex;
    TPZVec<STATE> epx;
    TPZVec<STATE> eey;
    TPZVec<STATE> epy;
    
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
            cel->Solution(parametric_space, eex_var, eex);
            cel->Solution(parametric_space, epx_var, epx);
            cel->Solution(parametric_space, eey_var, eey);
            cel->Solution(parametric_space, epy_var, epy);
            duplet.first = eex[0] + epx[0] + eey[0] + epy[0];
            duplet.second = p[0];
            
            duplet.first = fabs(duplet.first);
            duplet.second = fabs(duplet.second);
        }
    }
    
    m_strain_pressure_duplets.Push(duplet);
}

/** @brief Compute the strain and the stress at x euclidean point for each time */
void TPZPMRSAnalysis::PlotStrainStress(std::string file_name)
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

/** @brief Compute the strain and the Porosity at x euclidean point for each time */
void TPZPMRSAnalysis::PlotStrainPorosity(std::string file_name)
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

/** @brief Compute the strain and the Porosity at x euclidean point for each time */
void TPZPMRSAnalysis::PlotStrainPermeability(std::string file_name)
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

/** @brief Compute the strain and the Porosity at x euclidean point for each time */
void TPZPMRSAnalysis::PlotStrainPressure(std::string file_name)
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
