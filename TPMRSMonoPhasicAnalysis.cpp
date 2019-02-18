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
    m_error          = 0;
    m_dx_norm        = 0;
    m_k_iterations   = 0;
    m_post_processor = NULL;
    m_var_names.resize(0);
    m_vec_var_names.resize(0);
    
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
    m_vec_var_names     = other.m_vec_var_names;
    
}

void TPMRSMonoPhasicAnalysis::ConfigurateAnalysis(DecomposeType decomposition, TPZManVector<TPZCompMesh * , 2> & mesh_vec,TPMRSSimulationData * simulation_data){
    SetSimulationData(simulation_data);
    TPZStepSolver<STATE> step;
    unsigned int number_threads = m_simulation_data->n_threads();
    
    if(!Mesh()){
        std::cout << "Call SetCompMesh method." << std::endl;
        DebugStop();
    }
    m_mesh_vec = mesh_vec;

    
    switch (decomposition) {
        case ECholesky:
        {
            TPZSkylineStructMatrix struct_mat(Mesh());
            struct_mat.SetNumThreads(number_threads);
            this->SetStructuralMatrix(struct_mat);
        }
            break;
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
#ifdef USING_MKL2
            TPZSymetricSpStructMatrix struct_mat(Mesh());
            struct_mat.SetNumThreads(number_threads);
            this->SetStructuralMatrix(struct_mat);
#else
            TPZParFrontStructMatrix<TPZFrontSym<STATE> > struct_mat(Mesh());
            struct_mat.SetNumThreads(number_threads);
            this->SetStructuralMatrix(struct_mat);
#endif
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
    
    TPZManVector<std::string,50> scalnames,vecnames;
    TPZStack<std::string> names;
    scalnames = m_simulation_data->s_names_res();
    vecnames = m_simulation_data->v_names_res();
    for (auto i : scalnames) {
        names.push_back(i);
    }
    for (auto i : vecnames) {
        names.push_back(i);
    }

    m_post_processor->SetPostProcessVariables(post_mat_id, names);
    TPZFStructMatrix structmatrix(m_post_processor->Mesh());
    structmatrix.SetNumThreads(n_threads);
    m_post_processor->SetStructuralMatrix(structmatrix);
}

void TPMRSMonoPhasicAnalysis::ExecuteNewtonInteration(){

    Assemble();
    Solver().Matrix()->SetIsDecomposed(0);// Force numerical factorization
//    Solver().Matrix()->Print("j = ",std::cout, EMathematicaInput);
    Rhs() *= -1.0;
//    Rhs().Print("r = ",std::cout, EMathematicaInput);
    Solve();
//    Solution().Print("dp = ",std::cout,EMathematicaInput);
}

#define CheapNONM_Q

void TPMRSMonoPhasicAnalysis::ExecuteNinthOrderNewtonInteration(REAL & norm_dx){
    
    Assemble();
    TPZMatrix<REAL> * j_x = Solver().Matrix()->Clone();
    Rhs() *= -1.0;
    TPZFMatrix<REAL> r_x = Rhs();
    Solve();

#ifndef CheapNONM_Q
    TPZAutoPointer<TPZMatrix<REAL>> inv_j_x = Solver().Matrix()->Clone();
#endif
    
    TPZFMatrix<STATE> x_k,x,y,z,x_k_new;
    x_k = m_X_n;
    TPZFMatrix<STATE> dx,dyx;
    dx = Solution();
    y = x_k + dx;
    m_X_n = y;
    LoadMemorySolution();
    
#ifndef CheapNONM_Q
    Assemble();
#else
    AssembleResidual();
#endif
    
    TPZFMatrix<REAL> r_y = Rhs();
    Rhs() = r_x;
    Solve();
    dyx = Solution();
    TPZAutoPointer<TPZMatrix<REAL>> inv_j_y = Solver().Matrix()->Clone();
    
    z = x_k + 0.5*(dyx + dx);
    m_X_n = z;
    LoadMemorySolution();
    Rhs() *= -1.0;
    TPZFMatrix<REAL> r_z = Rhs();
    Solve();
    
    TPZFMatrix<REAL> temp_y = Solution();
    TPZFMatrix<REAL> temp_yy;
    j_x->Multiply(temp_y, temp_yy);
    Rhs() = temp_yy;
    Solve();
    TPZFMatrix<REAL> dz_1 = Solution();
    
#ifndef CheapNONM_Q
    Solver().UpdateFrom(inv_j_x);
#endif
    Rhs() = r_z;
    Solve();
    
    TPZFMatrix<REAL> dz_2 = Solution();
    TPZFMatrix<REAL> z_k_new = z + 0.5*(dz_1 + dz_2);

#ifndef CheapNONM_Q
    Solver().UpdateFrom(inv_j_y);
#endif
    m_X_n = z_k_new;
    LoadMemorySolution();
    Rhs() *= -1.0;
    r_z = Rhs();
    Solve();
    
    temp_y = Solution();
    j_x->Multiply(temp_y, temp_yy);
    Rhs() = temp_yy;
    Solve();
    dz_1 = Solution();
    
#ifndef CheapNONM_Q
    Solver().UpdateFrom(inv_j_x);
#endif
    Rhs() = r_z;
    Solve();
    
    dz_2 = Solution();
    
    x_k_new = z_k_new + 0.5*(dz_1 + dz_2);
    norm_dx = Norm(x_k_new - x_k);
    m_X_n = x_k_new;

}

#define NMO9_Q

void TPMRSMonoPhasicAnalysis::ExecuteOneTimeStep(){
    
    /// The process will update just the current state
    m_simulation_data->SetCurrentStateQ(true);
    LoadCurrentState();
    
    TPZFMatrix<STATE> dx;
    bool residual_stop_criterion_Q = false;
    bool correction_stop_criterion_Q = false;
    REAL norm_res, norm_dx;
    REAL r_norm = m_simulation_data->epsilon_res();
    REAL dx_norm = m_simulation_data->epsilon_cor();
    int n_it = m_simulation_data->n_iterations();
    
    for (int i = 1; i <= n_it; i++) {

#ifdef NMO9_Q
        /// https://www.sciencedirect.com/science/article/abs/pii/S0096300318302893
        if (i <= 2) {
            this->ExecuteNewtonInteration();
            dx = Solution();
            norm_dx  = Norm(dx);
            m_X_n += dx;
        }
        else{
            this->ExecuteNinthOrderNewtonInteration(norm_dx);
        }
#else
        this->ExecuteNewtonInteration();
        dx = Solution();
        norm_dx  = Norm(dx);
        m_X_n += dx;

#endif

        
        LoadMemorySolution();
        norm_res = Norm(Rhs());
        
        residual_stop_criterion_Q   = norm_res < r_norm;
        correction_stop_criterion_Q = norm_dx  < dx_norm;
        
        m_k_iterations = i;
        m_error = norm_res;
        m_dx_norm = norm_dx;
        

        if (residual_stop_criterion_Q /*&& correction_stop_criterion_Q */) {
#ifdef PZDEBUG
            std::cout << "TPMRSMonoPhasicAnalysis:: Nonlinear process converged with residue norm = " << norm_res << std::endl;
            std::cout << "TPMRSMonoPhasicAnalysis:: Correction norm = " << norm_dx << std::endl;
            std::cout << "TPMRSMonoPhasicAnalysis:: Number of iterations = " << i << std::endl;
#endif
            break;
        }
    }
    
    if (residual_stop_criterion_Q == false) {
        std::cout << "TPMRSMonoPhasicAnalysis:: Nonlinear process not converged with residue norm = " << norm_res << std::endl;
    }
}

void TPMRSMonoPhasicAnalysis::UpdateState(){
    m_simulation_data->SetTransferCurrentToLastQ(true);
    LoadMemorySolution();
    m_simulation_data->SetTransferCurrentToLastQ(false);
    m_X = m_X_n;
}

void TPMRSMonoPhasicAnalysis::PostProcessTimeStep(std::string & file){
    
    int dim = Mesh()->Dimension();
    int div = m_simulation_data->n_div();
    
    TPZManVector<std::string,50> scalnames;
    TPZManVector<std::string,50> vecnames;

    scalnames = m_simulation_data->s_names_res();
    vecnames = m_simulation_data->v_names_res();
    
    m_post_processor->TransferSolution();
    m_post_processor->DefineGraphMesh(dim, scalnames, vecnames,file);
    m_post_processor->PostProcess(div,dim);
}

void TPMRSMonoPhasicAnalysis::LoadMemorySolution(){

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


void TPMRSMonoPhasicAnalysis::ExecuteUndrainedResponseStep(){
    
    m_X_n.Zero();
    m_simulation_data->SetInitialStateQ(true);
    ExecuteNewtonInteration();
    m_X_n += Solution();
    LoadCurrentState();
    AssembleResidual();
    REAL norm = Norm(Rhs());
    std::cout << "TPMRSMonoPhasicAnalysis:: Undrained initial pressure projected with residual norm = " << norm << std::endl;
    m_X = m_X_n;
    m_simulation_data->SetInitialStateQ(false);
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


