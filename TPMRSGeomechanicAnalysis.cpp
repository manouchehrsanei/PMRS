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
    m_n_update_jac = 1;
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
    m_n_update_jac      = other.m_n_update_jac;
    m_post_processor    = other.m_post_processor;
    m_var_names         = other.m_var_names;
    
}

void TPMRSGeomechanicAnalysis::ConfigurateAnalysis(DecomposeType decomposition, TPMRSSimulationData * simulation_data){
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
#ifdef USING_MKL
            TPZSymetricSpStructMatrix struct_mat(Mesh());
            struct_mat.SetNumThreads(number_threads);
            this->SetStructuralMatrix(struct_mat);
#else
            TPZSkylineStructMatrix struct_mat(Mesh());
            struct_mat.SetNumThreads(number_threads);
            this->SetStructuralMatrix(struct_mat);
#endif
        }
            break;
        case ELU:
        {
            
#ifdef USING_MKL
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
    
    TPZManVector<std::string,50> scalnames,vecnames,tensnames;
    TPZStack<std::string> names;
    scalnames = m_simulation_data->s_names_geo();
    vecnames  = m_simulation_data->v_names_geo();
    tensnames = m_simulation_data->t_names_geo();
    
    for (auto i : scalnames) {
        names.push_back(i);
    }
    for (auto i : vecnames) {
        names.push_back(i);
    }
    for (auto i : tensnames) {
        names.push_back(i);
    }
    
    m_post_processor->SetPostProcessVariables(post_mat_id, names);
    TPZFStructMatrix structmatrix(m_post_processor->Mesh());
    structmatrix.SetNumThreads(n_threads);
    m_post_processor->SetStructuralMatrix(structmatrix);
}

void TPMRSGeomechanicAnalysis::ExecuteM1Interation(REAL & norm_dx){
    
    TPZFMatrix<STATE> x_k;
    x_k = Solution();
    if ((m_k_iterations)%m_n_update_jac) {
        AssembleResidual();
    }else{
        Assemble();
        Solver().Matrix()->SetIsDecomposed(0);// Force numerical factorization
        std::cout << "Jacobian updated at iteration = " << m_k_iterations << endl;
    }
    Rhs() *= -1.0;
    Solve();
    
    m_X_n = x_k + Solution();
    norm_dx = Norm(m_X_n - x_k);
    LoadSolution(m_X_n);
    
}

void TPMRSGeomechanicAnalysis::ExecuteM3Interation(REAL & norm_dx){
    
    TPZFMatrix<STATE> d_eps_x_x, d_eps_y_x ,x_k, y, dx;
    x_k = Solution();
    
    if ((m_k_iterations)%m_n_update_jac) {
        AssembleResidual();
    }else{
        Assemble();
        Solver().Matrix()->SetIsDecomposed(0);// Force numerical factorization
        std::cout << "First Jacobian updated at iteration = " << m_k_iterations << endl;
    }
    
    Rhs() *= -1.0;
    TPZFMatrix<REAL> r_x = Rhs();
    Solve();
    
    // Newton method (SecantQ = false)
    // Quase-Newton method (SecantQ = true)
    bool SecantQ = m_simulation_data->Get_is_secant_reservoir_Q();
    
    d_eps_x_x = Solution();
    y = x_k + d_eps_x_x;
    
    
    m_X_n = y;
    LoadSolution(y);
    LoadMemorySolution();
    
    if (SecantQ) {
        DebugStop();
        AssembleResidual();
    }else{
        if ((m_k_iterations)%m_n_update_jac) {
            AssembleResidual();
        }else{
            Assemble();
            Solver().Matrix()->SetIsDecomposed(0);// Force numerical factorization
            std::cout << "Second Jacobian updated at iteration = " << m_k_iterations << endl;
        }
    }
    
    Rhs() = r_x;
    Solve();
    d_eps_y_x = Solution();
    dx = 0.5*(d_eps_y_x + d_eps_x_x);
    
    m_X_n = x_k + dx;
    norm_dx = Norm(m_X_n - x_k);
    LoadSolution(m_X_n);
    
}

void TPMRSGeomechanicAnalysis::ExecuteM6Interation(REAL & norm_dx){
    
    TPZFMatrix<STATE> d_eps_x_x, d_eps_y_x, d_eps_y_z ,x_k, y, dx, z;
    x_k = Solution();
    
    if ((m_k_iterations)%m_n_update_jac) {
        AssembleResidual();
    }else{
        Assemble();
        Solver().Matrix()->SetIsDecomposed(0);// Force numerical factorization
        std::cout << "First Jacobian updated at iteration = " << m_k_iterations << endl;
    }
    
    TPZMatrix<REAL> * j_x = Solver().Matrix()->Clone();
    
    Rhs() *= -1.0;
    TPZFMatrix<REAL> r_x = Rhs();
    Solve();
    
    // Newton method (SecantQ = false)
    // Quase-Newton method (SecantQ = true)
    bool SecantQ = m_simulation_data->Get_is_secant_reservoir_Q();
    
    d_eps_x_x = Solution();
    y = x_k + d_eps_x_x;
    
    m_X_n = y;
    LoadSolution(y);
    LoadMemorySolution();
    if (SecantQ) {
        DebugStop();
        AssembleResidual();
    }else{
        if ((m_k_iterations)%m_n_update_jac) {
            AssembleResidual();
        }else{
            Assemble();
            Solver().Matrix()->SetIsDecomposed(0);// Force numerical factorization
            std::cout << "Second Jacobian updated at iteration = " << m_k_iterations << endl;
        }
    }
    
    Rhs() *= -1.0;
    TPZFMatrix<REAL> r_y = Rhs();
    
    Rhs() = r_x;
    Solve();
    
    d_eps_y_x = Solution();
    dx = 0.5*(d_eps_y_x + d_eps_x_x);
    z = x_k + dx;
    
    m_X_n = z;
    LoadSolution(z);
    LoadMemorySolution();
    
    Rhs() *= -1.0;
    TPZFMatrix<REAL> r_z = Rhs();
    Solve();
    d_eps_y_z = Solution();
    
    TPZFMatrix<REAL> j_x_d_eps_y_z;
    j_x->Multiply(d_eps_y_z, j_x_d_eps_y_z);
    Rhs() = j_x_d_eps_y_z;
    Solve();
    TPZFMatrix<REAL> dz_1 = Solution();
    
//    m_X_n = x_k;
//    LoadSolution(x_k);
//    LoadMemorySolution();
//    if (SecantQ) {
//        DebugStop();
//        AssembleResidual();
//    }else{
//        if ((m_k_iterations)%m_n_update_jac) {
//            AssembleResidual();
//        }else{
//            Assemble();
//            Solver().Matrix()->SetIsDecomposed(0);// Force numerical factorization
//            std::cout << "Third Jacobian updated at iteration = " << m_k_iterations << endl;
//        }
//    }
    
    Rhs() = r_z;
    Solve();
    TPZFMatrix<REAL> dz_2 = Solution();
    
    m_X_n = z + 0.5*(dz_1 + dz_2);
    norm_dx = Norm(m_X_n - x_k);
    LoadSolution(m_X_n);
}

void TPMRSGeomechanicAnalysis::ExecuteInteration(REAL & norm_dx){
    
    std::string method = m_simulation_data->name_nonlinear_Newton_method();
    
    if (method.compare("M1") == 0) {
        ExecuteM1Interation(norm_dx);
        return;
    }
    
    if (method.compare("M3") == 0) {
        if (m_k_iterations == 1) {
            ExecuteM1Interation(norm_dx);
        }else{
            ExecuteM3Interation(norm_dx);
        }
        return;
    }
    
    if (method.compare("M6") == 0) {
        if (m_k_iterations == 1) {
            ExecuteM1Interation(norm_dx);
        }else{
            ExecuteM6Interation(norm_dx);
        }
        return;
    }
    
    std::cout << "Nonlinear method not implemented = " << method.c_str() << std::endl;
    DebugStop();
    return;
}

#define CheapNONM_Q

void TPMRSGeomechanicAnalysis::ExecuteNinthOrderNewtonInteration(REAL & norm_dx){
    
    TPZFMatrix<STATE> x_k,x,y,z,x_k_new;
    x_k = Solution();
    
    if ((m_k_iterations-1)%m_n_update_jac) {
        AssembleResidual();
    }else{
        Assemble();
        Solver().Matrix()->SetIsDecomposed(0);// Force numerical factorization
        std::cout << "Jacobian updated at iteration = " << m_k_iterations << endl;
    }

#ifndef    CheapNONM_Q
    TPZMatrix<REAL> * j_x = Solver().Matrix()->Clone();
#else
    TPZAutoPointer<TPZMatrix<REAL>> j_x = Solver().Matrix();
#endif
    
    Rhs() *= -1.0;
    TPZFMatrix<REAL> r_x = Rhs();
    Solve();
    
#ifndef    CheapNONM_Q
    TPZAutoPointer<TPZMatrix<REAL>> inv_j_x = Solver().Matrix()->Clone();
#endif

    TPZFMatrix<STATE> dx,dyx;
    dx = Solution();
    y = x_k + dx;
    m_X_n = y;
    LoadSolution(y);
    LoadMemorySolution();
    
#ifndef    CheapNONM_Q
    Assemble();
#else
    AssembleResidual();
#endif
    
    TPZFMatrix<REAL> r_y = Rhs();
    Rhs() = r_x;
    Solve();
    dyx = Solution();
    
#ifndef    CheapNONM_Q
    TPZAutoPointer<TPZMatrix<REAL>> inv_j_y = Solver().Matrix()->Clone();
#endif
    
    z = x_k + 0.5*(dyx + dx);
    m_X_n = z;
    LoadSolution(z);
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
    
#ifndef    CheapNONM_Q
    Solver().UpdateFrom(inv_j_x);
#endif
    Rhs() = r_z;
    Solve();
    
    TPZFMatrix<REAL> dz_2 = Solution();
    TPZFMatrix<REAL> z_k_new = z + 0.5*(dz_1 + dz_2);
    
    
#ifndef    CheapNONM_Q
    Solver().UpdateFrom(inv_j_y);
#endif
    m_X_n = z_k_new;
    LoadSolution(z_k_new);
    LoadMemorySolution();
    Rhs() *= -1.0;
    r_z = Rhs();
    Solve();
    
    temp_y = Solution();
    j_x->Multiply(temp_y, temp_yy);
    Rhs() = temp_yy;
    Solve();
    dz_1 = Solution();
    
#ifndef    CheapNONM_Q
    Solver().UpdateFrom(inv_j_x);
#endif
    Rhs() = r_z;
    Solve();
    
    dz_2 = Solution();
    
    x_k_new = z_k_new + 0.5*(dz_1 + dz_2);
    norm_dx = Norm(x_k_new - x_k);
    m_X_n = x_k_new;
    LoadSolution(x_k_new);
}


#define NMO9_Q

//#define InnerLoopPerformance_Q

bool TPMRSGeomechanicAnalysis::ExecuteOneTimeStep(bool enforced_execution_Q){

    
    
    
    /// The nonlinear process will update just the current state
    m_simulation_data->SetCurrentStateQ(true);
    TPZFMatrix<STATE> dx(Solution());
    LoadSolution(); /// Load the current analysis solution on the cmesh
    bool residual_stop_criterion_Q   = false;
    bool correction_stop_criterion_Q = false;
    REAL norm_res, norm_dx;
    REAL r_norm  = m_simulation_data->epsilon_res();
    REAL dx_norm = m_simulation_data->epsilon_cor();
    int n_it     = m_simulation_data->n_iterations();
    
#ifdef InnerLoopPerformance_Q
    TPZFMatrix<REAL> residuals(n_it,2);
    residuals.Zero();
#endif
    
    for (int i = 1; i <= n_it; i++) {
        m_k_iterations = i;
        
        ExecuteInteration(dx_norm);
//        norm_dx  = Norm(Solution());
//        dx += Solution();
//
//        LoadSolution(dx);
//        m_X_n = dx;
        
        LoadMemorySolution();
        norm_res = Norm(this->Rhs());
        if (m_simulation_data->Get_must_use_sub_stepping_Q() && !enforced_execution_Q) {
            break;
        }
#ifdef InnerLoopPerformance_Q
        residuals(i-1,0) = norm_dx;
        residuals(i-1,1) = norm_res;
#endif
        residual_stop_criterion_Q   = norm_res < r_norm;
        correction_stop_criterion_Q = norm_dx  < dx_norm;
        
        m_error   = norm_res;
        m_dx_norm = norm_dx;
        std::cout << "TPMRSGeomechanicAnalysis:: residue norm = " << norm_res << std::endl;
        if (residual_stop_criterion_Q /*&&  correction_stop_criterion_Q*/) {
            std::cout << "TPMRSGeomechanicAnalysis:: Nonlinear process converged with residue norm = " << norm_res << std::endl;
            std::cout << "TPMRSGeomechanicAnalysis:: Correction norm = " << norm_dx << std::endl;
            std::cout << "TPMRSGeomechanicAnalysis:: Number of iterations = " << i << std::endl;
#ifdef InnerLoopPerformance_Q
            residuals.Print("rgeo = ",std::cout,EMathematicaInput);
#endif
            break;
        }
    }
    
    if (residual_stop_criterion_Q == false) {
        if (m_simulation_data->Get_must_use_sub_stepping_Q()) {
            std::cout << "TPMRSGeomechanicAnalysis:: Nonlinear need supstepping,  residue norm = " << norm_res << std::endl;
        }else{
            std::cout << "TPMRSGeomechanicAnalysis:: Nonlinear process not converged with residue norm = " << norm_res << std::endl;
        }
        
    }
    
    return m_simulation_data->Get_must_use_sub_stepping_Q();
}

void TPMRSGeomechanicAnalysis::ExecuteUndrainedResponseStep(){
    
    /// The process will update just the current state
    m_simulation_data->SetCurrentStateQ(true);
    TPZFMatrix<STATE> dx(Solution());
    bool residual_stop_criterion_Q   = false;
    bool correction_stop_criterion_Q = false;
    REAL norm_res, norm_dx;
    REAL r_norm  = m_simulation_data->epsilon_res();
    REAL dx_norm = m_simulation_data->epsilon_cor();
    int n_it = m_simulation_data->n_iterations();
    
    for (int i = 1; i <= n_it; i++) {
        this->ExecuteInteration(dx_norm);
//        dx += Solution();
//        norm_dx  = Norm(Solution());
//        LoadSolution(dx);
        LoadMemorySolution();
        norm_res = Norm(this->Rhs());
        residual_stop_criterion_Q   = norm_res < r_norm;
        correction_stop_criterion_Q = norm_dx  < dx_norm;
        
        m_k_iterations = i;
        m_error = norm_res;
        m_dx_norm = norm_dx;
        
        if (residual_stop_criterion_Q &&  correction_stop_criterion_Q) {
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
    
    TPZManVector<std::string,50> scalnames,vecnames,tensnames;
    scalnames = m_simulation_data->s_names_geo();
    vecnames  = m_simulation_data->v_names_geo();
    tensnames = m_simulation_data->t_names_geo();
    
    m_post_processor->TransferSolution();
    m_post_processor->DefineGraphMesh(dim,scalnames,vecnames,tensnames,file);
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


