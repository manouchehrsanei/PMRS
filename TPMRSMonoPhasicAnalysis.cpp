//
//  TPMRSMonoPhasicAnalysis.cpp
//  PMRS
//
//  Created by Omar and Manouchehr on 9/11/18.
//

#include "TPMRSMonoPhasicAnalysis.h"

TPMRSMonoPhasicAnalysis::TPMRSMonoPhasicAnalysis() : TPZAnalysis(){
    
    m_simulation_data = NULL;
    m_X_n.Resize(0, 0);
    m_X.Resize(0, 0);
    m_mesh_vec.Resize(0);
    m_error = 0;
    m_dx_norm = 0;
    m_k_iterations = 0;
    m_post_processor = NULL;
    m_var_names.resize(0);
    
}

TPMRSMonoPhasicAnalysis::~TPMRSMonoPhasicAnalysis(){

}

TPMRSMonoPhasicAnalysis::TPMRSMonoPhasicAnalysis(const TPMRSMonoPhasicAnalysis & other){
    
    m_simulation_data   = other.m_simulation_data;
    m_X_n               = other.m_X_n;
    m_X                 = other.m_X;
    m_mesh_vec          = other.m_mesh_vec;
    m_error             = other.m_error;
    m_dx_norm           = other.m_dx_norm;
    m_k_iterations      = other.m_k_iterations;
    m_post_processor    = other.m_post_processor;
    m_var_names         = other.m_var_names;
    
}

void TPMRSMonoPhasicAnalysis::ConfigurateAnalysis(DecomposeType decomposition, TPZManVector<TPZCompMesh * , 2> & mesh_vec,TPZSimulationData * simulation_data){
    SetSimulationData(simulation_data);
    TPZStepSolver<STATE> step;
    unsigned int number_threads = m_simulation_data->n_threads();
    
    if(!Mesh()){
        std::cout << "Call SetCompMesh method." << std::endl;
        DebugStop();
    }
    m_mesh_vec = mesh_vec;

    
    switch (decomposition) {
        case ELU:
        {
            
#ifdef USING_MKL2
            TPZSpStructMatrix struct_mat(Mesh());
            struct_mat.SetNumThreads(number_threads);
            this->SetStructuralMatrix(struct_mat);
#else
            TPZSkylineNSymStructMatrix struct_mat(Mesh());
            struct_mat.SetNumThreads(number_threads);
            this->SetStructuralMatrix(struct_mat);
#endif
            
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
    TPZManVector<std::pair<int, TPZManVector<int,12>>,12>  material_ids = m_simulation_data->MaterialIds();
    TPZManVector<int,10> post_mat_id(n_regions);
    for (int iregion = 0; iregion < n_regions; iregion++)
    {
        int matid = material_ids[iregion].first;
        post_mat_id[iregion] = matid;
    }
    
    // @TODO:: MS, please transfer from xml file
    m_var_names.Push("p");
    m_var_names.Push("phi");
    m_post_processor->SetPostProcessVariables(post_mat_id, m_var_names);
    
    TPZFStructMatrix structmatrix(m_post_processor->Mesh());
    structmatrix.SetNumThreads(n_threads);
    m_post_processor->SetStructuralMatrix(structmatrix);
    
}

void TPMRSMonoPhasicAnalysis::ExecuteNewtonInteration(){
    Assemble();
//    Solver().Matrix()->Print("j = ",std::cout, EMathematicaInput);
    Rhs() *= -1.0;
//    Rhs().Print("r = ",std::cout, EMathematicaInput);
    Solve();
//    Solution().Print("dp = ",std::cout,EMathematicaInput);
}

void TPMRSMonoPhasicAnalysis::ExecuteOneTimeStep(){
    
    if (m_simulation_data->IsInitialStateQ()) {
        m_X = Solution();
    }
    
    m_simulation_data->SetCurrentStateQ(false);
    AcceptTimeStepSolution();
    
//    // Initial guess
//    m_X_n = m_X;
    m_simulation_data->SetCurrentStateQ(true);
    this->AcceptTimeStepSolution();
    
    TPZFMatrix<STATE> dx;
    bool residual_stop_criterion_Q = false;
    bool correction_stop_criterion_Q = false;
    REAL norm_res, norm_dx;
    REAL r_norm = m_simulation_data->epsilon_res();
    REAL dx_norm = m_simulation_data->epsilon_cor();
    int n_it = m_simulation_data->n_iterations();
    
    for (int i = 1; i <= n_it; i++) {
        this->ExecuteNewtonInteration();
        dx = Solution();
        norm_dx  = Norm(dx);
        m_X_n += dx;

        this->AcceptTimeStepSolution();
        norm_res = Norm(Rhs());
        residual_stop_criterion_Q   = norm_res < r_norm;
        correction_stop_criterion_Q = norm_dx  < dx_norm;
        
        m_k_iterations = i;
        m_error = norm_res;
        m_dx_norm = norm_dx;
        
        if (residual_stop_criterion_Q ||  correction_stop_criterion_Q) {
#ifdef PZDEBUG
            std::cout << "TPMRSMonoPhasicAnalysis:: Nonlinear process converged with residue norm = " << norm_res << std::endl;
            std::cout << "TPMRSMonoPhasicAnalysis:: Number of iterations = " << i << std::endl;
            std::cout << "TPMRSMonoPhasicAnalysis:: Correction norm = " << norm_dx << std::endl;
#endif
            break;
        }
    }
    
    if (residual_stop_criterion_Q == false) {
        std::cout << "TPMRSMonoPhasicAnalysis:: Nonlinear process not converged with residue norm = " << norm_res << std::endl;
    }
}

void TPMRSMonoPhasicAnalysis::UpdateState(){
    m_X = m_X_n;
}

void TPMRSMonoPhasicAnalysis::PostProcessTimeStep(std::string & file){
    
    int dim = Mesh()->Dimension();
    int div = m_simulation_data->n_div();
    TPZStack< std::string> vecnames;
    m_post_processor->TransferSolution();
    m_post_processor->DefineGraphMesh(dim,m_var_names,vecnames,file);
    m_post_processor->PostProcess(div,dim);
}

void TPMRSMonoPhasicAnalysis::AcceptTimeStepSolution(){

    bool state = m_simulation_data->IsCurrentStateQ();
    if (state) {
        m_simulation_data->Set_must_accept_solution_Q(true);
        LoadCurrentState();
        AssembleResidual();
        m_simulation_data->Set_must_accept_solution_Q(false);
    }else{
        m_simulation_data->Set_must_accept_solution_Q(true);
        LoadLastState();
        AssembleResidual();
        m_simulation_data->Set_must_accept_solution_Q(false);
    }
}


void TPMRSMonoPhasicAnalysis::LoadCurrentState(){
    LoadSolution(m_X_n);
    if(m_simulation_data->Get_is_dual_formulation_Q()){
        TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(m_mesh_vec, Mesh());
    }
}

void TPMRSMonoPhasicAnalysis::LoadLastState(){
    LoadSolution(m_X);
    if(m_simulation_data->Get_is_dual_formulation_Q()){
        TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(m_mesh_vec, Mesh());
    }
}


