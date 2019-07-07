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
    m_n_update_jac   = 1;
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
    m_n_update_jac      = other.m_n_update_jac;
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
        case ELDLt:
        {
#ifdef USING_MKL
            TPZSymetricSpStructMatrix struct_mat(Mesh());
            struct_mat.SetNumThreads(number_threads);
            this->SetStructuralMatrix(struct_mat);
#else
            TPZParFrontStructMatrix<TPZFrontSym<STATE> > struct_mat(Mesh());
            struct_mat.SetNumThreads(number_threads);
            this->SetStructuralMatrix(struct_mat);
//
//            TPZSkylineStructMatrix struct_mat(Mesh());
//            struct_mat.SetNumThreads(number_threads);
//            this->SetStructuralMatrix(struct_mat);
            
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

void TPMRSMonoPhasicAnalysis::ExecuteM1Interation(REAL & norm_dx){
    
    // Newton method (SecantQ = false)
    // Quase-Newton method (SecantQ = true)
    bool SecantQ = m_simulation_data->Get_is_secant_reservoir_Q();
    
    if ((m_k_iterations)%m_n_update_jac) {
        AssembleResidual();
    }else{
        if (SecantQ) {
            if (m_k_iterations <= 2) {
                Assemble();
                Solver().Matrix()->SetIsDecomposed(0);// Force numerical factorization
                std::cout << "First Jacobian updated at iteration = " << m_k_iterations << endl;
            }else{
                AssembleResidual();
            }
        }else{
            Assemble();
            Solver().Matrix()->SetIsDecomposed(0);// Force numerical factorization
            std::cout << "First Jacobian updated at iteration = " << m_k_iterations << endl;
        }
    }
    
    Rhs() *= -1.0;
    Solve();
    
    TPZFMatrix<STATE> dx = Solution();
    norm_dx = Norm(dx);
    m_X_n += dx; // update.
}


void TPMRSMonoPhasicAnalysis::ExecuteMTInteration(REAL & norm_dx){

    //    https://www.newcastle.edu.au/__data/assets/pdf_file/0004/22576/44_Accelerated-initial-stiffness-schemes-for-elastoplasticity.pdf
    
    TPZFMatrix<STATE> x_k;
    x_k = m_X_n;
    
    // Newton method (SecantQ = false)
    // Quase-Newton method (SecantQ = true)
    bool SecantQ = m_simulation_data->Get_is_secant_reservoir_Q();
    
    if ((m_k_iterations)%m_n_update_jac) {
        AssembleResidual();
    }else{
        if (SecantQ) {
            if (m_k_iterations <= 2) {
                Assemble();
                Solver().Matrix()->SetIsDecomposed(0);// Force numerical factorization
                std::cout << "First Jacobian updated at iteration = " << m_k_iterations << endl;
            }else{
                AssembleResidual();
            }
        }else{
            Assemble();
            Solver().Matrix()->SetIsDecomposed(0);// Force numerical factorization
            std::cout << "First Jacobian updated at iteration = " << m_k_iterations << endl;
        }
    }
    
    Rhs() *= -1.0;
    Solve();
    
    m_X_n = x_k +  m_alpha*Solution();
    norm_dx = Norm(m_X_n - x_k);
    
    LoadMemorySolution();
    if (SecantQ) {
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
    Solve();
    
    m_X_tilde = Solution();
    /// Compute the new alpha Equation 17.
    
    int n_equ = m_X_tilde.Rows();
    REAL num = 0, dem = 0;
    for (int i = 0; i < n_equ; i++) {
        num += (m_X_n(i,0) - x_k(i,0))*(m_X_tilde(i,0));
        dem += (m_X_n(i,0) - x_k(i,0))*(m_X_n(i,0) - x_k(i,0));
    }
    REAL s = num/dem;
    m_alpha += s;
    
    /// Perform equation 15.
    m_X_n += m_X_tilde;
    norm_dx = Norm(m_X_n - x_k);
    LoadSolution(m_X_n);
    
}

void TPMRSMonoPhasicAnalysis::ExecuteM3Interation(REAL & norm_dx){
    

    TPZFMatrix<STATE> d_eps_x_x, d_eps_y_x ,x_k, y;
    x_k = m_X_n;
    
    // Newton method (SecantQ = false)
    // Quase-Newton method (SecantQ = true)
    bool SecantQ = m_simulation_data->Get_is_secant_reservoir_Q();
    
    if ((m_k_iterations)%m_n_update_jac) {
        AssembleResidual();
    }else{
        if (SecantQ) {
            if (m_k_iterations <= 2) {
                Assemble();
                Solver().Matrix()->SetIsDecomposed(0);// Force numerical factorization
                std::cout << "First Jacobian updated at iteration = " << m_k_iterations << endl;
            }else{
                AssembleResidual();
            }
        }else{
            Assemble();
            Solver().Matrix()->SetIsDecomposed(0);// Force numerical factorization
            std::cout << "First Jacobian updated at iteration = " << m_k_iterations << endl;
        }
    }
    
    Rhs() *= -1.0;
    TPZFMatrix<REAL> r_x = Rhs();
    Solve();

    d_eps_x_x = Solution();
    y = x_k + d_eps_x_x;
    
    m_X_n = y;
    LoadMemorySolution();
    
//    if (SecantQ) {
//        AssembleResidual();
//    }else{
//        if ((m_k_iterations)%m_n_update_jac) {
//            AssembleResidual();
//        }else{
//            Assemble();
//            Solver().Matrix()->SetIsDecomposed(0);// Force numerical factorization
//            std::cout << "Second Jacobian updated at iteration = " << m_k_iterations << endl;
//        }
//    }
    
    Rhs() = r_x;
    Solve();
    
    d_eps_y_x = Solution();
    
    m_X_n = x_k + 0.5*(d_eps_y_x + d_eps_x_x);// update.
    norm_dx = Norm(m_X_n-x_k);
    
    
}

void TPMRSMonoPhasicAnalysis::ExecuteM6Interation(REAL & norm_dx){
    
    TPZFMatrix<STATE> d_eps_x_x, d_eps_y_x, d_eps_y_z ,x_k, y, dzx, z;
    x_k = m_X_n;
    
    // Newton method (SecantQ = false)
    // Quase-Newton method (SecantQ = true)
    bool SecantQ = m_simulation_data->Get_is_secant_reservoir_Q();
    
    if ((m_k_iterations)%m_n_update_jac) {
        AssembleResidual();
    }else{
        if (SecantQ) {
            if (m_k_iterations <= 2) {
                Assemble();
                Solver().Matrix()->SetIsDecomposed(0);// Force numerical factorization
                std::cout << "First Jacobian updated at iteration = " << m_k_iterations << endl;
            }else{
                AssembleResidual();
            }
        }else{
            Assemble();
            Solver().Matrix()->SetIsDecomposed(0);// Force numerical factorization
            std::cout << "First Jacobian updated at iteration = " << m_k_iterations << endl;
        }
    }
    
    TPZMatrix<REAL> * j_x = Solver().Matrix()->Clone();
    
    Rhs() *= -1.0;
    TPZFMatrix<REAL> r_x = Rhs();
    Solve();
    
    
    d_eps_x_x = Solution();
    y = x_k + d_eps_x_x;
    
    m_X_n = y;
    LoadMemorySolution();
    
    if (SecantQ) {
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
    dzx = 0.5*(d_eps_y_x + d_eps_x_x);
    z = x_k + dzx;
    
    m_X_n = z;
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
    
    m_X_n = x_k;
    LoadMemorySolution();
    
//    if (SecantQ) {
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
    
    m_X_n = z + 0.5*(dz_1 + dz_2);// update.
    norm_dx = Norm(m_X_n-x_k);
    
}

#define CheapNONM_Q

void TPMRSMonoPhasicAnalysis::ExecuteInteration(REAL & norm_dx){
    
    std::string method = m_simulation_data->name_nonlinear_Newton_method();
    
    if (method.compare("M1") == 0) {
        ExecuteM1Interation(norm_dx);
        return;
    }
    
    if (method.compare("MT") == 0) {
        ExecuteMTInteration(norm_dx);
        return;
    }

    if (method.compare("M3") == 0) {
        if (m_k_iterations <= 2) {
            ExecuteM1Interation(norm_dx);
        }else{
            ExecuteM3Interation(norm_dx);
        }
        return;
    }
    
    if (method.compare("M6") == 0) {
        if (m_k_iterations <= 2) {
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

void TPMRSMonoPhasicAnalysis::ExecuteNinthOrderNewtonInteration(REAL & norm_dx){
    
    if ((m_k_iterations)%m_n_update_jac) {
        AssembleResidual();
    }else{
        Assemble();
        Solver().Matrix()->SetIsDecomposed(0);// Force numerical factorization
        std::cout << "Jacobian updated at iteration = " << m_k_iterations << endl;
    }
    
    // Newton method (SecantQ = false)
    // Quase-Newton method (SecantQ = true)
    bool SecantQ = m_simulation_data->Get_is_secant_reservoir_Q();
    
    
//    TPZFMatrix<STATE> dx,x_k;
//    x_k = m_X_n;
//    dx = Solution();
//
//    TPZFMatrix<STATE> y = x_k + dx;
    
    
    
//    TPZFMatrix<STATE> x_k,x,y,z,x_k_new;
//    TPZAutoPointer<TPZMatrix<REAL>> j_x = Solver().Matrix();
//    if (SecantQ) {
//        <#statements#>
//    }
    
    

#ifndef    CheapNONM_Q
    TPZMatrix<REAL> * j_x = Solver().Matrix()->Clone();
#else
    TPZAutoPointer<TPZMatrix<REAL>> j_x = Solver().Matrix();
#endif
    
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
    
#ifndef    CheapNONM_Q
    TPZAutoPointer<TPZMatrix<REAL>> inv_j_y = Solver().Matrix()->Clone();
#endif
    
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
//#define InnerLoopPerformance_Q

void TPMRSMonoPhasicAnalysis::ExecuteOneTimeStep(){
    
    /// The nonlinear process will update just the current state
    m_simulation_data->SetCurrentStateQ(true);
    m_k_iterations = 0;
    
    if (m_simulation_data->Get_use_internal_accel_ResQ()) {
        ApplyAcceleration();
    }
    
    LoadMemorySolution();
    
    TPZFMatrix<STATE> dx;
    bool residual_stop_criterion_Q = false;
    bool correction_stop_criterion_Q = false;
    REAL norm_res, norm_dx;
    REAL error_tol  = m_simulation_data->epsilon_res();
    REAL dx_tol     = m_simulation_data->epsilon_cor();
    int n_it        = m_simulation_data->n_iterations();
    
#ifdef InnerLoopPerformance_Q
    TPZFMatrix<REAL> residuals(n_it,2);
    residuals.Zero();
#endif
    
    for (int i = 1; i <= n_it; i++) {
        m_k_iterations = i;
        
        ExecuteInteration(norm_dx);
        
        if (m_simulation_data->Get_use_internal_accel_ResQ()) {
            ApplyAcceleration();
        }
        
        LoadMemorySolution();
        norm_res = Norm(Rhs());
        
#ifdef InnerLoopPerformance_Q
        residuals(i-1,0) = norm_dx;
        residuals(i-1,1) = norm_res;
#endif
        
        residual_stop_criterion_Q   = norm_res < error_tol;
        correction_stop_criterion_Q = norm_dx  < dx_tol;
    
        m_error = norm_res;
        m_dx_norm = norm_dx;
        
        std::cout << "TPMRSMonoPhasicAnalysis:: residue norm = " << norm_res << std::endl;
        if (residual_stop_criterion_Q /*&& correction_stop_criterion_Q */) {
            std::cout << "TPMRSMonoPhasicAnalysis:: Nonlinear process converged with residue norm = " << norm_res << std::endl;
            std::cout << "TPMRSMonoPhasicAnalysis:: Correction norm = " << norm_dx << std::endl;
            std::cout << "TPMRSMonoPhasicAnalysis:: Number of iterations = " << i << std::endl;
#ifdef InnerLoopPerformance_Q
            residuals.Print("rres = ",std::cout,EMathematicaInput);
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
    REAL norm_dx;
    ExecuteInteration(norm_dx);
    m_X_n += Solution();
    LoadCurrentState();
    AssembleResidual();
    REAL norm = Norm(Rhs());
    std::cout << "TPMRSMonoPhasicAnalysis:: Undrained initial pressure projected with correction norm = " << norm_dx << std::endl;
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



/// Define an acceleration method for the internal loop k iteration for reservoir module
void TPMRSMonoPhasicAnalysis::AccelerationRes(int k, int n){
    
//    k--;
    int n_terms;
    {
        if (k < n) {
            n_terms = k;
        }
        
        if (k >= n) {
            n_terms = n;
        }
    }
    
    switch (n_terms) {
        case 0:
        {
            m_x_p.Resize(1);
            m_x_p[0] = m_X_n;
        }
            break;
        case 1:
        {
            m_x_p.Resize(2);
            m_x_p[1] = m_X_n;
            m_X_n = ApplyTransformation(m_x_p[1], m_x_p[1], m_x_p[0]);
            
        }
            break;
        case 2:  /// S(A_n)
        {
            m_x_p.Resize(3);
            m_x_p[2] = m_X_n;
            m_X_n = ApplyTransformation(m_x_p[2], m_x_p[1], m_x_p[0]);
            
        }
            break;
        case 3:  /// S(A_n)
        {
            
            m_x_p.Resize(4);
            m_x_p[3] = m_X_n;
            m_X_n = ApplyTransformation(m_x_p[3], m_x_p[2], m_x_p[1]);
            
        }
            break;
        case 4:  /// S2(A_n)
        {
            m_x_p.Resize(5);
            m_x_p[4] = m_X_n;
            
            TPZFMatrix<REAL> Sk1,Sk2,Sk3;
            Sk1 = ApplyTransformation(m_x_p[2], m_x_p[1], m_x_p[0]);
            Sk2 = ApplyTransformation(m_x_p[3], m_x_p[2], m_x_p[1]);
            Sk3 = ApplyTransformation(m_x_p[4], m_x_p[3], m_x_p[2]);
            m_X_n = ApplyTransformation(Sk3,Sk2,Sk1);
            
        }
            break;
        case 5:  /// S2andhalf(A_n)
        {
            
            m_x_p.Resize(6);
            m_x_p[5] = m_X_n;
            
            TPZFMatrix<REAL> Sk1,Sk2,Sk3;
            Sk1 = ApplyTransformation(m_x_p[3], m_x_p[2], m_x_p[1]);
            Sk2 = ApplyTransformation(m_x_p[4], m_x_p[3], m_x_p[2]);
            Sk3 = ApplyTransformation(m_x_p[5], m_x_p[4], m_x_p[3]);
            m_X_n = ApplyTransformation(Sk3,Sk2,Sk1);
            
        }
            break;
        case 6:  /// S3(A_n)
        {
            m_x_p.Resize(7);
            m_x_p[6] = m_X_n;
            
            TPZFMatrix<REAL>Sk1,Sk2,Sk3,S2k1,S2k2,S2k3;
            Sk1 = ApplyTransformation(m_x_p[2], m_x_p[1], m_x_p[0]);
            Sk2 = ApplyTransformation(m_x_p[3], m_x_p[2], m_x_p[1]);
            Sk3 = ApplyTransformation(m_x_p[4], m_x_p[3], m_x_p[2]);
            S2k1 = ApplyTransformation(Sk3,Sk2,Sk1);
            
            Sk1 = ApplyTransformation(m_x_p[3], m_x_p[2], m_x_p[1]);
            Sk2 = ApplyTransformation(m_x_p[4], m_x_p[3], m_x_p[2]);
            Sk3 = ApplyTransformation(m_x_p[5], m_x_p[4], m_x_p[3]);
            S2k2 = ApplyTransformation(Sk3,Sk2,Sk1);
            
            Sk1 = ApplyTransformation(m_x_p[4], m_x_p[3], m_x_p[2]);
            Sk2 = ApplyTransformation(m_x_p[5], m_x_p[4], m_x_p[3]);
            Sk3 = ApplyTransformation(m_x_p[6], m_x_p[5], m_x_p[4]);
            S2k3 = ApplyTransformation(Sk3,Sk2,Sk1);
            
            m_X_n = ApplyTransformation(S2k3,S2k2,S2k1);
            
        }
            break;
        default:
            break;
    }
    
}

TPZFMatrix<REAL> TPMRSMonoPhasicAnalysis::ApplyTransformation(TPZFMatrix<REAL> & An_p_1, TPZFMatrix<REAL> & An, TPZFMatrix<REAL> & An_m_1){
    std::string nonlinear_acceleration = m_simulation_data->name_nonlinear_acceleration();
    TPZFMatrix<REAL> S(An_p_1);
    if (nonlinear_acceleration == "FDM"){
        S = FDMTransformation(An_p_1, An, An_m_1);
    }
    else if (nonlinear_acceleration == "SDM"){
        S = SDMTransformation(An_p_1, An, An_m_1);
    }
    
#ifdef Adaptive_Acceleration_Q
    REAL theta_ratio = 1.0;
    REAL max_theta   = m_simulation_data->Get_max_theta_value();
    bool apply_transformation_Q = true;
    {
        int n_dof = S.Rows();
        REAL num = 0.0;
        REAL den = 0.0;
        for (int i = 0; i < n_dof ; i++) {
            num += (S(i,0)-An_p_1(i,0))*(An(i,0)-An_m_1(i,0));
            den += (An_p_1(i,0)-An(i,0))*(An_p_1(i,0)-An(i,0));
        }
        if (!IsZero(den)) {
            theta_ratio = num / den;
        }
        apply_transformation_Q = fabs(theta_ratio-1.0) < max_theta;
    }
    
    if (apply_transformation_Q) {
        return S;
    }else{
        std::cout << " TPMRSSegregatedAnalysis:: Transformation is avoided swithcing to SFI." << std::endl;
        return An_p_1;
    }
#else
    return S;
#endif
    
    
}

TPZFMatrix<REAL> TPMRSMonoPhasicAnalysis::FDMTransformation(TPZFMatrix<REAL> & An_p_1, TPZFMatrix<REAL> & An, TPZFMatrix<REAL> & An_m_1){
    TPZFMatrix<REAL> S(An_p_1);
    int n_dof = S.Rows();
    
    REAL num = 0.0;
    REAL den = 0.0;
    REAL w;
    for (int i = 0; i < n_dof ; i++) {
        w    = An_m_1(i,0)-An(i,0);
        num += w*(An(i,0)-An_p_1(i,0));
        den += w*(An_m_1(i,0) - 2*An(i,0) + An_p_1(i,0));
    }
    REAL s;
    if (IsZero(den)) {
        s = num / den;
    }else{
        s = num / den;
    }
    S = An_p_1-An;
    S *= s;
    S += An_p_1;
    return S;
}

TPZFMatrix<REAL> TPMRSMonoPhasicAnalysis::SDMTransformation(TPZFMatrix<REAL> & An_p_1, TPZFMatrix<REAL> & An, TPZFMatrix<REAL> & An_m_1){
    
    TPZFMatrix<REAL> S(An_p_1);
    int n_dof = S.Rows();
    
    REAL num = 0.0;
    REAL den = 0.0;
    REAL w;
    for (int i = 0; i < n_dof ; i++) {
        w    = An_m_1(i,0) - 2*An(i,0) + An_p_1(i,0);
        num += w*(An(i,0)-An_p_1(i,0));
        den += w*w;
    }
    REAL s;
    if (IsZero(den)) {
        s = num / den;
    }else{
        s = num / den;
    }
    S = An_p_1-An;
    S *= s;
    S += An_p_1;
    return S;
}

void TPMRSMonoPhasicAnalysis::ApplyAcceleration(){
    
    // Applying the selected nonlinear acceleration
    std::string nonlinear_acceleration = m_simulation_data->name_nonlinear_acceleration();
    bool non_linear_acceleration_Q = (nonlinear_acceleration == "FDM") || (nonlinear_acceleration == "SDM");
    if (non_linear_acceleration_Q) {
        
        int n_terms = m_simulation_data->n_state_acceleration(); /// n=2->S, n=4->S2, and n=6->S3
        
        int n_vec = m_x_p.size();
        if (m_k_iterations > n_terms) {
            for (int i = 0; i < n_vec - 1; i++) {
                m_x_p[i] = m_x_p[i+1];
            }
            if(n_vec!=0){
                m_x_p[n_vec-1] = m_X_n;
            }
        }
        
        AccelerationRes(m_k_iterations,n_terms);
        
    }
    
}
