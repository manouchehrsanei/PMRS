//
//  TPMRSFullyCoupledAnalysis.cpp
//  PZ
//
//  Created by Omar and Manouchehr on 8/28/16.
//
//

#include "TPMRSFullyCoupledAnalysis.h"
#include "pzpostprocanalysis.h"
#include "pzfstrmatrix.h"
#include "TPZElasticCriterion.h"

TPMRSFullyCoupledAnalysis::TPMRSFullyCoupledAnalysis() : TPZAnalysis()
{
    m_simulation_data = NULL;
    m_meshvec.Resize(2);
    m_X_n.Resize(0,0);
    m_X.Resize(0,0);
    m_error        = 1.0;
    m_dx_norm      = 1.0;
    m_k_iterations = 0;
    
}

TPMRSFullyCoupledAnalysis::~TPMRSFullyCoupledAnalysis(){
    
}

TPMRSFullyCoupledAnalysis::TPMRSFullyCoupledAnalysis(const TPMRSFullyCoupledAnalysis &other)
{
    m_simulation_data = other.m_simulation_data;
    m_meshvec         = other.m_meshvec;
    m_X_n             = other.m_X_n;
    m_X               = other.m_X;
    m_error           = other.m_error;
    m_dx_norm         = other.m_dx_norm;
    
}

TPMRSFullyCoupledAnalysis & TPMRSFullyCoupledAnalysis::operator=(const TPMRSFullyCoupledAnalysis &other)
{
    /// check for self-assignment
    if(&other == this){
        return *this;
    }
    
    m_simulation_data = other.m_simulation_data;
    m_meshvec         = other.m_meshvec;
    m_X_n             = other.m_X_n;
    m_X               = other.m_X;
    m_error           = other.m_error;
    m_dx_norm         = other.m_dx_norm;
    
    return *this;
}

void TPMRSFullyCoupledAnalysis::ConfiguratePostProcessor(){
    
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
    
    scalnames = m_simulation_data->s_names_res();
    vecnames  = m_simulation_data->v_names_res();
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

void TPMRSFullyCoupledAnalysis::AdjustVectors()
{
    
    if(fSolution.Rows() == 0 /* || fRhs.Rows() == 0 */)
    {
        DebugStop();
    }
    
    /// @TODO:: check the need for this!
    TPZBuildMultiphysicsMesh::AddElements(m_meshvec, this->Mesh());
    TPZBuildMultiphysicsMesh::AddConnects(m_meshvec, this->Mesh());
    TPZBuildMultiphysicsMesh::TransferFromMeshes(m_meshvec, this->Mesh());
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(m_meshvec, this->Mesh());
    
    m_X.Resize(fSolution.Rows(),1);
    m_X.Zero();
    m_X_n.Resize(fSolution.Rows(),1);
    m_X_n.Zero();
}

void TPMRSFullyCoupledAnalysis::ExecuteNewtonInteration(){
    Assemble();
    Solver().Matrix()->SetIsDecomposed(0);// Force numerical factorization
    Rhs() *= -1.0;
    Solve();
}

#define CheapNONM_Q

void TPMRSFullyCoupledAnalysis::ExecuteNinthOrderNewtonInteration(REAL & norm_dx){
    
    TPZFMatrix<STATE> x_k,x,y,z,x_k_new;
    x_k = Solution();
    
    Assemble();
    TPZMatrix<REAL> * j_x = Solver().Matrix()->Clone();
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
//    LoadSolution(y);
    LoadMemorySolution(y);
    
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
//    LoadSolution(z);
    LoadMemorySolution(z);
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
//    LoadSolution(z_k_new);
    LoadMemorySolution(z_k_new);
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
    LoadMemorySolution(x_k_new);
}

void TPMRSFullyCoupledAnalysis::ExecuteOneTimeStep(int i_time_step)
{

    TPZFMatrix<REAL> last_solution = this->Solution();
    
    REAL psudo_time;
    REAL dt = m_simulation_data->dt();
    int n_level = 0;
    bool enforced_execution_Q = false;
    int n_sub_steps = power(2,m_simulation_data->Get_n_sub_step_level());
    for (int i = 1; i <= n_sub_steps; i++) {
        REAL delta = 1.0/REAL(n_sub_steps);
        REAL alpha = delta*i;
        psudo_time = (i_time_step+1)*dt - (1.0-alpha)*dt;
        {   /// Interpolate loads
            /// The Fully coupled boundary conditions
            ConfigureFCBC(psudo_time);
        }
        bool check_for_sub_stepping_Q = this->ExcecuteOneStepApproximation(enforced_execution_Q);
        if (check_for_sub_stepping_Q && !enforced_execution_Q) {
            n_level++;
            m_simulation_data->Set_n_sub_step_level(n_level);
            n_sub_steps = power(2,n_level);
            if (n_level > 5) {
                n_sub_steps = 64;
                std::cout << "TPMRSSegregatedAnalysis:: The level for substepping is not enough = " << n_level << std::endl;
                std::cout << "TPMRSSegregatedAnalysis:: The number of substeps is fixed at = " << n_sub_steps << std::endl;
                std::cout << "--------------------- Reached the plasticity change tolerance -------------- " << std::endl;
                enforced_execution_Q = true;
            }
            /// It is required to restart the simulation
            i = 0;
            LoadState(last_solution);
            m_simulation_data->Set_must_use_sub_stepping_Q(false);
            std::cout << "TPMRSSegregatedAnalysis:: Increase the level for substepping = " << n_level << std::endl;
            std::cout << "TPMRSSegregatedAnalysis:: Current number of substeps = " << n_sub_steps << std::endl;
            std::cout << "--------------------- Restarting step simulation -------------- " << std::endl;
            std::cout << std::endl;
            std::cout << std::endl;
        }
        else{
            m_simulation_data->Set_must_use_sub_stepping_Q(true);
            UpdateState();
            m_simulation_data->Set_must_use_sub_stepping_Q(false);
        }
    }
    
    if(n_sub_steps > 1){
        std::cout << "TPMRSSegregatedAnalysis:: Geomechanics solved with level of substepping = " << n_level << std::endl;
        std::cout << "TPMRSSegregatedAnalysis:: Current number of substeps = " << n_sub_steps << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;
    }
    
    m_simulation_data->Set_must_use_sub_stepping_Q(false);
    m_simulation_data->Set_n_sub_step_level(0);
}

#define NMO9_Q

bool TPMRSFullyCoupledAnalysis::ExcecuteOneStepApproximation(bool enforced_execution_Q){
    
    /// The nonlinear process will update just the current state
    m_simulation_data->SetCurrentStateQ(true);
    TPZFMatrix<STATE> dx(Solution());
    LoadMemorySolution(dx);
    
    bool residual_stop_criterion_Q = false;
    bool correction_stop_criterion_Q = false;
    REAL norm_res, norm_dx;
    REAL r_norm = m_simulation_data->epsilon_res();
    REAL dx_norm = m_simulation_data->epsilon_cor();
    int n_it = m_simulation_data->n_iterations();
//    dx.Print("x = ");
    for (int i = 1; i <= n_it; i++) {
        
#ifdef NMO9_Q
        /// https://www.sciencedirect.com/science/article/abs/pii/S0096300318302893
        if (i <= 2) {
            this->ExecuteNewtonInteration();
            dx += Solution();
            norm_dx  = Norm(Solution());
            LoadMemorySolution(dx);
//            m_X_n = dx;
            
        }
        else{
            this->ExecuteNinthOrderNewtonInteration(norm_dx);
        }
#else
        
        this->ExecuteNewtonInteration();
        dx += Solution();
//        Solution().Print("dx = ");
//        dx.Print("xnew = ");
        norm_dx  = Norm(Solution());
        LoadMemorySolution(dx);
//        m_X_n = dx;
        
#endif
        
//        this->Rhs().Print("r(xnew) = ");
        norm_res = Norm(this->Rhs());
        if (m_simulation_data->Get_must_use_sub_stepping_Q() && !enforced_execution_Q) {
            break;
        }
        
        residual_stop_criterion_Q   = norm_res < r_norm;
        correction_stop_criterion_Q = norm_dx  < dx_norm;
        
        m_k_iterations = i;
        m_error = norm_res;
        m_dx_norm = norm_dx;
        
        if (residual_stop_criterion_Q || correction_stop_criterion_Q) {
#ifdef PZDEBUG
            std::cout << "TPMRSFullyCoupledAnalysis:: Nonlinear process converged with residue norm = " << norm_res << std::endl;
            std::cout << "TPMRSFullyCoupledAnalysis:: Correction norm = " << norm_dx << std::endl;
            std::cout << "TPMRSFullyCoupledAnalysis:: Number of iterations = " << i << std::endl;
#endif
            
            break;
        }
    }
    
    return m_simulation_data->Get_must_use_sub_stepping_Q();
}


/// update last state (at n state) solution for PMRS_PoroElastic
void TPMRSFullyCoupledAnalysis::UpdateState()
{
    m_simulation_data->SetTransferCurrentToLastQ(true);
    LoadMemorySolution(Solution());
    m_simulation_data->SetTransferCurrentToLastQ(false);
}


void TPMRSFullyCoupledAnalysis::ConfigureFCBC(REAL t){
    
    TPZCompMesh * cmesh = this->Mesh();
    if (!cmesh)
    {
        DebugStop();
    }
    
    int n_regions = m_simulation_data->NumberOfRegions();
    TPZManVector<std::pair<int, std::pair<TPZManVector<int,12>,TPZManVector<int,12>> >,12>  material_ids = m_simulation_data->MaterialIds();
    
    std::map<int, std::string>::iterator it_bc_id_to_type;
    std::map< std::string,std::pair<int,std::vector<std::string> > >::iterator  it_condition_type_to_index_value_names;
    std::map<int , TPMRSInterpolator >::iterator it_bc_id_to_values;
    
    for (int iregion = 0; iregion < n_regions; iregion++)
    {
        /// Inserting boundary conditions
        int n_bc = material_ids[iregion].second.first.size();
        for (int ibc = 0; ibc < n_bc; ibc++)
        {
            int bc_id = material_ids[iregion].second.first [ibc];
            
            it_bc_id_to_type = m_simulation_data->BCIdToConditionTypeFullyCoupled().find(bc_id);
            it_bc_id_to_values = m_simulation_data->BCIdToBCValuesFullyCoupled().find(bc_id);
            
            it_condition_type_to_index_value_names = m_simulation_data->ConditionTypeToBCIndexFullyCoupled().find(it_bc_id_to_type->second);
            
            int bc_index = it_condition_type_to_index_value_names->second.first;
            int n_bc_values = it_bc_id_to_values->second.n_functions();
            
            TPZMaterial * bc_mat = cmesh->FindMaterial(bc_id);
            if (!bc_mat) {
                DebugStop();
            }
            TPZBndCond * bc = dynamic_cast<TPZBndCond *>(bc_mat);
            if (!bc) {
                DebugStop();
            }
            
            bc->SetType(bc_index);
            bc->Val2().Resize(n_bc_values, 1);
            std::vector<REAL> f_values = it_bc_id_to_values->second.f(t);
            for (int i = 0; i < n_bc_values; i++) {
                REAL value = f_values[i];
                bc->Val2()(i,0) = value;
            }
        }
    }
    
}

void TPMRSFullyCoupledAnalysis::PostProcessTimeStep(std::string file)
{
    int dim = Mesh()->Dimension();
    int div = m_simulation_data->n_div();
    
    TPZManVector<std::string,50> scalnames_g,vecnames_g,tensnames_g;
    TPZManVector<std::string,50> scalnames_r,vecnames_r;
    TPZStack<std::string,50> scalnames,vecnames,tensnames;
    
    scalnames_g = m_simulation_data->s_names_geo();
    vecnames_g  = m_simulation_data->v_names_geo();
    tensnames_g = m_simulation_data->t_names_geo();
    
    scalnames_r = m_simulation_data->s_names_res();
    vecnames_r  = m_simulation_data->v_names_res();
    
    for (auto i : scalnames_g) {
        scalnames.push_back(i);
    }
    for (auto i : scalnames_r) {
        scalnames.push_back(i);
    }
    
    for (auto i : vecnames_g) {
        vecnames.push_back(i);
    }
    for (auto i : vecnames_r) {
        vecnames.push_back(i);
    }
    
    for (auto i : tensnames_g) {
        tensnames.push_back(i);
    }
    
    m_post_processor->TransferSolution();
    m_post_processor->DefineGraphMesh(dim,scalnames,vecnames,tensnames,file);
    m_post_processor->PostProcess(div,dim);
    
}

void TPMRSFullyCoupledAnalysis::ExecuteTimeEvolution(){
    
    /// vtk file
    std::string name = m_simulation_data->name_vtk_file();
    std::string file = name + "_fc.vtk";
    
    int n_time_steps = m_simulation_data->n_steps();
    REAL dt          = m_simulation_data->dt();
    REAL time_value  = 0.0;
    
    
    std::cout << std::endl;
    std::cout << "-------------------------------------------------------------" <<std::endl;
    std::cout << "-------------------------------------------------------------" <<std::endl;
    std::cout << "TPMRSFullyCoupledAnalysis:: Begining for recurrent simulation process." <<std::endl;
    std::cout << std::endl << std::endl;
    
    for (int it = 0; it < n_time_steps; it++) {
        /// Interpolate BC data
        time_value = dt * (it+1);
        ConfigureFCBC(time_value);
        
#ifdef USING_BOOST
            boost::posix_time::ptime fss_t1 = boost::posix_time::microsec_clock::local_time();
#endif
            this->ExecuteOneTimeStep(it);
#ifdef USING_BOOST
            boost::posix_time::ptime fss_t2 = boost::posix_time::microsec_clock::local_time();
#endif
            
#ifdef USING_BOOST
            REAL fc_solving_time = boost::numeric_cast<double>((fss_t2-fss_t1).total_milliseconds());
            std::cout << "TPMRSFullyCoupledAnalysis:: Fully coupled process closed in :" << setw(10) <<  fc_solving_time/1000.0 << setw(5)   << " seconds." << std::endl;
            std::cout << std::endl;
            
#endif
        
        bool postprocess_time_Q = ShouldPostprocessQ(time_value);
        if (postprocess_time_Q) {
            PostProcessTimeStep(file);
        }
        UpdateState();
    }
    
    std::cout << std::endl;
    std::cout << "-------------------------------------------------------------" <<std::endl;
    std::cout << "-------------------------------------------------------------" <<std::endl;
    std::cout << "TPMRSFullyCoupledAnalysis:: Ending for recurrent simulation process." <<std::endl;
    std::cout << std::endl;
    std::cout << std::endl << std::endl;
    
}

bool TPMRSFullyCoupledAnalysis::ShouldPostprocessQ(REAL time){
    TPZStack<REAL,500> & reporting_times = m_simulation_data->ReportingTimes();
    bool postprocess_time_Q = false;
    int n_times = reporting_times.size();
    for (int it = 0; it < n_times; it++) {
        bool check = IsZero(reporting_times[it] - time);
        if (check) {
            postprocess_time_Q = true;
        }
    }
    return postprocess_time_Q;
}

void TPMRSFullyCoupledAnalysis::LoadMemorySolution(TPZFMatrix<REAL> & x){
    
    bool state = m_simulation_data->IsCurrentStateQ();
    if (state) {
        m_simulation_data->Set_must_accept_solution_Q(true);
        LoadState(x);
        AssembleResidual();
        m_simulation_data->Set_must_accept_solution_Q(false);
    }else{
        m_simulation_data->Set_must_accept_solution_Q(true);
        LoadState(x);
        AssembleResidual();
        m_simulation_data->Set_must_accept_solution_Q(false);
    }
}

void TPMRSFullyCoupledAnalysis::LoadState(TPZFMatrix<REAL> & x){
    Solution() = x;
    LoadSolution(x);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(m_meshvec, Mesh());
}
