//
//  TPMRSGeomechanicAnalysis.cpp
//  PMRS
//
//  Created by Omar and Manouchehr on 9/11/18.
//

#include "TPMRSGeomechanicAnalysis.h"

TPMRSGeomechanicAnalysis::TPMRSGeomechanicAnalysis() : TPZAnalysis(){
    
    m_simulation_data = NULL;
    m_X_n.Resize(0, 0);
    m_X.Resize(0, 0);
    m_error = 0;
    m_dx_norm = 0;
    m_k_iterations = 0;
    m_post_processor = NULL;
    m_var_names.resize(0);
    
}

TPMRSGeomechanicAnalysis::~TPMRSGeomechanicAnalysis(){
    
}

TPMRSGeomechanicAnalysis::TPMRSGeomechanicAnalysis(const TPMRSGeomechanicAnalysis & other) :  TPZAnalysis(other){
    
    m_simulation_data   = other.m_simulation_data;
    m_X_n               = other.m_X_n;
    m_X                 = other.m_X;
    m_error             = other.m_error;
    m_dx_norm           = other.m_dx_norm;
    m_k_iterations      = other.m_k_iterations;
    m_post_processor    = other.m_post_processor;
    m_var_names         = other.m_var_names;
    
}

void TPMRSGeomechanicAnalysis::ConfigurateAnalysis(DecomposeType decomposition, TPZSimulationData * simulation_data){
    SetSimulationData(simulation_data);
    TPZStepSolver<STATE> step;
    unsigned int number_threads = m_simulation_data->n_threads();
    
    if(!Mesh()){
        std::cout << "Call SetCompMesh method." << std::endl;
        DebugStop();
    }
    
    switch (decomposition) {
        case ECholesky:
        {
            TPZSkylineStructMatrix struct_mat(Mesh());
            struct_mat.SetNumThreads(number_threads);
            this->SetStructuralMatrix(struct_mat);
        }
            break;
        case ELDLt:
        {
            TPZSymetricSpStructMatrix struct_mat(Mesh());
            struct_mat.SetNumThreads(number_threads);
            this->SetStructuralMatrix(struct_mat);
        }
            break;
        default:
        {
            DebugStop();
        }
            break;
    }
    step.SetDirect(decomposition);
    this->SetSolver(step);
    this->Solution().Resize(Mesh()->Solution().Rows(), 1);
    m_X.Resize(Mesh()->Solution().Rows(), 1);
    m_X_n.Resize(Mesh()->Solution().Rows(), 1);
    
    int n_threads = m_simulation_data->n_threads();
    m_post_processor = new TPZPostProcAnalysis;
    m_post_processor->SetCompMesh(Mesh());
    
    int n_regions = m_simulation_data->NumberOfRegions();
    TPZManVector<std::pair<int, std::pair<TPZManVector<int,12>,TPZManVector<int,12>> >,12>  material_ids = m_simulation_data->MaterialIds();
    TPZManVector<int,10> post_mat_id(n_regions);
    for (int iregion = 0; iregion < n_regions; iregion++)
    {
        int matid = material_ids[iregion].first;
        post_mat_id[iregion] = matid;
    }
    
    // @TODO:: MS, please transfer from xml file
    m_var_names.Push("ux");
    m_var_names.Push("uy");
    m_var_names.Push("sxx");
    m_var_names.Push("syy");
    m_var_names.Push("szz");
    m_var_names.Push("exx");
    m_var_names.Push("eyy");
    m_var_names.Push("ezz");
    m_var_names.Push("epxx");
    m_var_names.Push("epyy");
    m_var_names.Push("epzz");
    
    if (m_simulation_data->Dimension() == 3) {
        m_var_names.Push("uz");
    }
    
    m_post_processor->SetPostProcessVariables(post_mat_id, m_var_names);
    
    TPZFStructMatrix structmatrix(m_post_processor->Mesh());
    structmatrix.SetNumThreads(n_threads);
    m_post_processor->SetStructuralMatrix(structmatrix);
    
}

void TPMRSGeomechanicAnalysis::ExecuteNewtonInteration(){
    Assemble();
    Rhs() *= -1.0;
    Solve();
}

void TPMRSGeomechanicAnalysis::ExecuteOneTimeStep(){
    
    if (m_simulation_data->IsInitialStateQ()) {
        m_X = Solution();
    }

    // The process will update just the current state
    m_simulation_data->SetCurrentStateQ(true);
    TPZFMatrix<STATE> dx(Solution());
    bool residual_stop_criterion_Q = false;
    bool correction_stop_criterion_Q = false;
    REAL norm_res, norm_dx;
    REAL r_norm = m_simulation_data->epsilon_res();
    REAL dx_norm = m_simulation_data->epsilon_cor();
    int n_it = m_simulation_data->n_iterations();
    
    for (int i = 1; i <= n_it; i++) {
        this->ExecuteNewtonInteration();
        dx += Solution();
        norm_dx  = Norm(Solution());
        LoadSolution(dx);
        LoadMemorySolution();
        norm_res = Norm(this->Rhs());
        residual_stop_criterion_Q   = norm_res < r_norm;
        correction_stop_criterion_Q = norm_dx  < dx_norm;
        
        m_k_iterations = i;
        m_error = norm_res;
        m_dx_norm = norm_dx;
        
        if (residual_stop_criterion_Q ||  correction_stop_criterion_Q) {
#ifdef PZDEBUG
            std::cout << "TPMRSGeomechanicAnalysis:: Nonlinear process converged with residue norm = " << norm_res << std::endl;
            std::cout << "TPMRSGeomechanicAnalysis:: Correction norm = " << norm_dx << std::endl;
            std::cout << "TPMRSGeomechanicAnalysis:: Number of iterations = " << i << std::endl;
#endif

            break;
        }
    }
    
    if (residual_stop_criterion_Q == false) {
        std::cout << "TPMRSGeomechanicAnalysis:: Nonlinear process not converged with residue norm = " << norm_res << std::endl;
    }
}

void TPMRSGeomechanicAnalysis::UpdateState(){
    m_simulation_data->SetTransferCurrentToLastQ(true);
    LoadMemorySolution();
    m_simulation_data->SetTransferCurrentToLastQ(false);
}

void TPMRSGeomechanicAnalysis::PostProcessTimeStep(std::string & file){
    
    int dim = Mesh()->Dimension();
    int div = m_simulation_data->n_div();
    TPZStack< std::string> vecnames;
    m_post_processor->TransferSolution();
    m_post_processor->DefineGraphMesh(dim,m_var_names,vecnames,file);
    m_post_processor->PostProcess(div,dim);
}

void TPMRSGeomechanicAnalysis::LoadMemorySolution(){
    
    bool state = m_simulation_data->IsCurrentStateQ();
    if (state) {
        m_simulation_data->Set_must_accept_solution_Q(true);
        AssembleResidual();
        m_simulation_data->Set_must_accept_solution_Q(false);
    }else{
        m_simulation_data->Set_must_accept_solution_Q(true);
        AssembleResidual();
        m_simulation_data->Set_must_accept_solution_Q(false);
    }
}


void TPMRSGeomechanicAnalysis::LoadCurrentState(){
    LoadSolution(m_X_n);
    DebugStop();
}

void TPMRSGeomechanicAnalysis::LoadLastState(){
    LoadSolution(m_X);
    DebugStop();
}


