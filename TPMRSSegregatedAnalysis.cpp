//
//  TPMRSSegregatedAnalysis.cpp
//  PMRS
//
//  Created by Omar and Manouchehr on 9/13/18.
//

#include "TPMRSSegregatedAnalysis.h"

TPMRSSegregatedAnalysis::TPMRSSegregatedAnalysis(){
    m_simulation_data       = NULL;
    m_geomechanic_analysis  = NULL;
    m_reservoir_analysis    = NULL;
    m_iterations_summary.Resize(0, 0);
    m_cpu_time_summary.Resize(0, 0);
    m_residuals_summary.Resize(0, 0);
    m_x_p.Resize(0);
    m_x_u.Resize(0);
}

TPMRSSegregatedAnalysis::~TPMRSSegregatedAnalysis(){
    
}

TPMRSSegregatedAnalysis::TPMRSSegregatedAnalysis(const TPMRSSegregatedAnalysis & other){
    m_simulation_data       = other.m_simulation_data;
    m_geomechanic_analysis  = other.m_geomechanic_analysis;
    m_reservoir_analysis    = other.m_reservoir_analysis;
    m_iterations_summary    = other.m_iterations_summary;
    m_cpu_time_summary      = other.m_cpu_time_summary;
    m_residuals_summary     = other.m_residuals_summary;
    m_x_p                   = other.m_x_p;
    m_x_u                   = other.m_x_u;
}

TPMRSGeomechanicAnalysis * TPMRSSegregatedAnalysis::GetGeomechanicsSolver(){
    return m_geomechanic_analysis;
}

TPMRSMonoPhasicAnalysis * TPMRSSegregatedAnalysis::GetReservoirSolver(){
    return m_reservoir_analysis;
}

void TPMRSSegregatedAnalysis::ApplyMemoryLink(TPZCompMesh * cmesh_o, TPZCompMesh * cmesh_d){
    
    if (!m_simulation_data) {
        DebugStop();
    }
    
    int n_regions = m_simulation_data->NumberOfRegions();
    TPZManVector<std::pair<int, std::pair<TPZManVector<int,12>,TPZManVector<int,12>> >,12>  material_ids = m_simulation_data->MaterialIds();
    TPZManVector<int,10> volumetric_mat_id(n_regions);
    for (int iregion = 0; iregion < n_regions; iregion++)
    {
        int matid = material_ids[iregion].first;
        volumetric_mat_id[iregion] = matid;
        
        TPZMaterial * material_o = cmesh_o->FindMaterial(matid);
        TPZMaterial * material_d = cmesh_d->FindMaterial(matid);
        if (!material_o || !material_d) {
            DebugStop();
        }
        
        TPZMatWithMem<TPMRSMemory> * mat_with_memory_o = dynamic_cast<TPZMatWithMem<TPMRSMemory> * >(material_o);
        TPZMatWithMem<TPMRSMemory> * mat_with_memory_d = dynamic_cast<TPZMatWithMem<TPMRSMemory> * >(material_d);
        if (!mat_with_memory_o || !mat_with_memory_d) {
            DebugStop();
        }
        
        mat_with_memory_d->SetMemory(mat_with_memory_o->GetMemory());
    }
}

void TPMRSSegregatedAnalysis::ConfigurateAnalysis(DecomposeType decompose_geo, DecomposeType decompose_res, TPMRSSimulationData * simulation_data, TPZCompMesh * cmesh_geomechanics, TPZCompMesh * cmesh_reservoir, TPZManVector<TPZCompMesh * , 2> & mesh_vec){

    if (!simulation_data || !cmesh_geomechanics || !cmesh_reservoir) {
        DebugStop();
    }
    
    this->SetSimulationData(simulation_data);
    bool mustOptimizeBandwidth = true;

    if (simulation_data->ElasticityOrder()==simulation_data->DiffusionOrder() && simulation_data->Get_is_dual_formulation_Q()) {
        this->AdjustIntegrationOrder(cmesh_reservoir,cmesh_geomechanics);
        this->ApplyMemoryLink(cmesh_reservoir,cmesh_geomechanics);
    }else{
        this->AdjustIntegrationOrder(cmesh_geomechanics,cmesh_reservoir);
        this->ApplyMemoryLink(cmesh_geomechanics,cmesh_reservoir);
    }

    /// The Geomechanics Simulator
    m_geomechanic_analysis = new TPMRSGeomechanicAnalysis;
    m_geomechanic_analysis->SetCompMesh(cmesh_geomechanics,mustOptimizeBandwidth);
    m_geomechanic_analysis->ConfigurateAnalysis(decompose_geo, m_simulation_data);
    
    /// The Reservoir Simulator
    m_reservoir_analysis = new TPMRSMonoPhasicAnalysis;
    m_reservoir_analysis->SetCompMesh(cmesh_reservoir,mustOptimizeBandwidth);
    m_reservoir_analysis->ConfigurateAnalysis(decompose_res, mesh_vec, m_simulation_data);
    
    /// Loading spatial properties
    FillMemory(m_reservoir_analysis->Mesh());

}

void TPMRSSegregatedAnalysis::FillMemory(TPZCompMesh * cmesh){
    
    if (!m_simulation_data || !cmesh) {
        DebugStop();
    }
    
    
    int n_regions = m_simulation_data->NumberOfRegions();
    TPZManVector<std::pair<int, std::pair<TPZManVector<int,12>,TPZManVector<int,12>> >,12>  material_ids = m_simulation_data->MaterialIds();
    TPZManVector<int,10> volumetric_mat_id(n_regions);
    for (int iregion = 0; iregion < n_regions; iregion++)
    {
        int matid = material_ids[iregion].first;
        volumetric_mat_id[iregion] = matid;
        
        TPZMaterial * material = cmesh->FindMaterial(matid);
        if (!material) {
            DebugStop();
        }
        
        TPZMatWithMem<TPMRSMemory> * mat_with_memory = dynamic_cast<TPZMatWithMem<TPMRSMemory> * >(material);
        if (!material) {
            DebugStop();
        }
        
        std::tuple<TPMRSUndrainedParameters, TPMRSPoroMechParameters, TPMRSPhiParameters,TPMRSKappaParameters,TPMRSPlasticityParameters> chunk =    m_simulation_data->MaterialProps()[iregion];
        
        TPMRSUndrainedParameters udrained_parameters(std::get<0>(chunk));
        TPMRSPoroMechParameters poro_parameters(std::get<1>(chunk));
        std::vector<REAL> undrained_pars = udrained_parameters.GetParameters();
        std::vector<REAL> poroperm_pars = poro_parameters.GetParameters();

        REAL phi_0   = undrained_pars[2];
        REAL kappa_0 = undrained_pars[3];
        REAL Eyoung  = poroperm_pars[0];
        REAL nu      = poroperm_pars[1];
        REAL Ks      = poroperm_pars[2];
        REAL Kdr     = Eyoung/(3.0*(1.0-2.0*nu));
        REAL alpha   = 1.0 - (Kdr/Ks);
        std::shared_ptr<TPZAdmChunkVector<TPMRSMemory>> & memory_vector = mat_with_memory->GetMemory();
        
        int ndata = memory_vector->NElements();
        for (int i = 0; i < ndata; i++) {
            memory_vector.get()->operator [](i).Setphi_0(phi_0);
            memory_vector.get()->operator [](i).Setkappa_0(kappa_0);
            memory_vector.get()->operator [](i).SetKs(Ks);
            memory_vector.get()->operator [](i).SetKdr(Kdr);
            memory_vector.get()->operator [](i).SetAlpha(alpha);
        }
        
    }
}


void TPMRSSegregatedAnalysis::AdjustIntegrationOrder(TPZCompMesh * cmesh_o, TPZCompMesh * cmesh_d){
    
    int dim   = cmesh_o->Dimension();
    int nel_o = cmesh_o->NElements();
    int nel_d = cmesh_d->NElements();

    
    if (nel_o != nel_d) {
        std::cout << "The geometrical partitions are not the same." << std::endl;
        DebugStop();
    }
    
    /// Assuming the cmesh_o as directive.
    cmesh_d->LoadReferences();
    for (long el = 0; el < nel_o; el++) {
        TPZCompEl *cel_o = cmesh_o->Element(el);
        if (!cel_o) {
            continue;
        }
        
        TPZGeoEl * gel = cel_o->Reference();
        if (!gel) {
            continue;
        }
        
        /// Finding the other computational element
        TPZCompEl * cel_d = gel->Reference();
        if (!cel_d) {
            continue;
        }
        
        if (gel->Dimension()!=dim) { /// boundary memory is no shared
//            cel_o->SetFreeIntPtIndices();
//            cel_o->ForcePrepareIntPtIndices();
//            cel_d->SetFreeIntPtIndices();
//            cel_d->ForcePrepareIntPtIndices();
            continue;
        }
        
//        cel_o->SetFreeIntPtIndices();
//        cel_o->ForcePrepareIntPtIndices();
        const TPZIntPoints & rule = cel_o->GetIntegrationRule();
        TPZIntPoints * cloned_rule = rule.Clone();
        
        TPZManVector<int64_t,20> indices;
        cel_o->GetMemoryIndices(indices);
        cel_d->SetFreeIntPtIndices();
        cel_d->SetMemoryIndices(indices);
        cel_d->SetIntegrationRule(cloned_rule);
    }
    
#ifdef PZDEBUG
    std::ofstream out_geo("Cmesh_origin_adjusted.txt");
    cmesh_o->Print(out_geo);
#endif
    

#ifdef PZDEBUG
    std::ofstream out_res("Cmesh_destination_adjusted.txt");
    cmesh_d->Print(out_res);
#endif

}

//#define Acceleration_P_output_Q

void TPMRSSegregatedAnalysis::ExecuteOneTimeStep(int i_time_step, int k){
    
#ifdef USING_BOOST
    boost::posix_time::ptime res_t1 = boost::posix_time::microsec_clock::local_time();
#endif
    
    m_reservoir_analysis->ExecuteOneTimeStep();
    
#ifdef USING_BOOST
    boost::posix_time::ptime res_t2 = boost::posix_time::microsec_clock::local_time();
#endif
    
#ifdef USING_BOOST
    REAL res_solving_time = boost::numeric_cast<double>((res_t2-res_t1).total_milliseconds());
    std::cout << "TPMRSMonoPhasicAnalysis:: Newton process closed in :" << setw(10) <<  res_solving_time/1000.0 << setw(5)   << " seconds." << std::endl;
    std::cout << std::endl;
    
    m_cpu_time_summary(i_time_step,1) += res_solving_time;

#endif
    
#ifdef USING_BOOST
    boost::posix_time::ptime geo_t1 = boost::posix_time::microsec_clock::local_time();
#endif
    
    ExecuteTheGeomechanicalApproximation(i_time_step);
    
#ifdef USING_BOOST
    boost::posix_time::ptime geo_t2 = boost::posix_time::microsec_clock::local_time();
#endif
    
#ifdef USING_BOOST
    REAL geo_solving_time = boost::numeric_cast<double>((geo_t2-geo_t1).total_milliseconds());
    std::cout << "TPMRSGeomechanicAnalysis:: Newton process closed in :" << setw(10) <<  geo_solving_time/1000.0 << setw(5)   << " seconds." << std::endl;
    std::cout << std::endl;
    
    m_cpu_time_summary(i_time_step,2) += geo_solving_time;
#endif
    

    // Applying the selected nonlinear acceleration
    std::string nonlinear_acceleration = m_simulation_data->name_nonlinear_acceleration();
    bool non_linear_acceleration_Q = (nonlinear_acceleration == "Shank") || (nonlinear_acceleration == "FDM") || (nonlinear_acceleration == "SDM");
    if (non_linear_acceleration_Q) {
        
        int n_terms = m_simulation_data->n_state_acceleration(); /// n=2->S, n=4->S2, and n=6->S3
        /// Acceleration for the whole thing geo + res
        AccelerationRes(k,n_terms);
        AccelerationGeo(k,n_terms);

        int n_vec = m_x_p.size();
        if (k-1 >= n_terms) {
            for (int i = 0; i < n_vec - 1; i++) {
                m_x_p[i] = m_x_p[i+1];
                m_x_u[i] = m_x_u[i+1];
            }
            if(n_vec!=0){
                m_x_p[n_vec-1] = m_reservoir_analysis->X_n();
                m_x_u[n_vec-1] = m_geomechanic_analysis->Solution();
            }
        }
        
    }
    
}

void TPMRSSegregatedAnalysis::ExecuteTheGeomechanicalApproximation(int i_time_step){
    
    TPZFMatrix<REAL> res_dx  = m_reservoir_analysis->X_n()-m_reservoir_analysis->X();
    TPZFMatrix<REAL> last_du = m_geomechanic_analysis->Solution();
    
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
            
            /// The biot pressure term
            m_reservoir_analysis->X_n() = m_reservoir_analysis->X() + alpha*res_dx;
            m_reservoir_analysis->LoadMemorySolution();
            
            /// The Geomechanic boundary conditions
            ConfigureGeomechanicsBC(psudo_time);
        }
        bool check_for_sub_stepping_Q = m_geomechanic_analysis->ExecuteOneTimeStep(enforced_execution_Q);
        if (check_for_sub_stepping_Q && !enforced_execution_Q) {
            n_level++;
            m_simulation_data->Set_n_sub_step_level(n_level);
            n_sub_steps = power(2,n_level);
            if (n_level > 6) {
                n_sub_steps = 128;
                std::cout << "TPMRSSegregatedAnalysis:: The level for substepping is not enough = " << n_level << std::endl;
                std::cout << "TPMRSSegregatedAnalysis:: The number of substeps is fixed at = " << n_sub_steps << std::endl;
                std::cout << "--------------------- Reached the plasticity change tolerance -------------- " << std::endl;
                enforced_execution_Q = true;
            }
            /// It is required to restart the simulation
            i = 0;
            m_geomechanic_analysis->LoadSolution(last_du);
            m_simulation_data->Set_must_use_sub_stepping_Q(false);
            std::cout << "TPMRSSegregatedAnalysis:: Increase the level for substepping = " << n_level << std::endl;
            std::cout << "TPMRSSegregatedAnalysis:: Current number of substeps = " << n_sub_steps << std::endl;
            std::cout << "--------------------- Restarting step simulation -------------- " << std::endl;
            std::cout << std::endl;
            std::cout << std::endl;
        }
        else{
            m_simulation_data->Set_must_use_sub_stepping_Q(true);
            m_geomechanic_analysis->UpdateState();
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


void TPMRSSegregatedAnalysis::AccelerationGeo(int k, int n){
    
    k--;
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
            m_x_u.Resize(1);
            m_x_u[0] = m_geomechanic_analysis->Solution();
        }
            break;
        case 1:
        {
            m_x_u.Resize(2);
            m_x_u[1] = m_geomechanic_analysis->Solution();
            m_geomechanic_analysis->Solution() = ApplyTransformation(m_x_u[1], m_x_u[1], m_x_u[0]);
            
        }
            break;
        case 2:  /// S(A_n)
        {
            m_x_u.Resize(3);
            m_x_u[2] = m_geomechanic_analysis->Solution();
            m_geomechanic_analysis->Solution() = ApplyTransformation(m_x_u[2], m_x_u[1], m_x_u[0]);
            
        }
            break;
        case 3:  /// S(A_n)
        {
            
            m_x_u.Resize(4);
            m_x_u[3] = m_geomechanic_analysis->Solution();
            m_geomechanic_analysis->Solution() = ApplyTransformation(m_x_u[3], m_x_u[2], m_x_u[1]);
            
        }
            break;
        case 4:  /// S2(A_n)
        {
            m_x_u.Resize(5);
            m_x_u[4] = m_geomechanic_analysis->Solution();
            
            TPZFMatrix<REAL> Sk1,Sk2,Sk3;
            Sk1 = ApplyTransformation(m_x_u[2], m_x_u[1], m_x_u[0]);
            Sk2 = ApplyTransformation(m_x_u[3], m_x_u[2], m_x_u[1]);
            Sk3 = ApplyTransformation(m_x_u[4], m_x_u[3], m_x_u[2]);
            m_geomechanic_analysis->Solution() = ApplyTransformation(Sk3,Sk2,Sk1);
            
        }
            break;
        case 5:  /// S2andhalf(A_n)
        {
            
            m_x_u.Resize(6);
            m_x_u[5] = m_geomechanic_analysis->Solution();
            
            TPZFMatrix<REAL> Sk1,Sk2,Sk3;
            Sk1 = ApplyTransformation(m_x_u[3], m_x_u[2], m_x_u[1]);
            Sk2 = ApplyTransformation(m_x_u[4], m_x_u[3], m_x_u[2]);
            Sk3 = ApplyTransformation(m_x_u[5], m_x_u[4], m_x_u[3]);
            m_geomechanic_analysis->Solution() = ApplyTransformation(Sk3,Sk2,Sk1);
            
        }
            break;
        case 6:  /// S3(A_n)
        {
            m_x_u.Resize(7);
            m_x_u[6] = m_geomechanic_analysis->Solution();
            
            TPZFMatrix<REAL>Sk1,Sk2,Sk3,S2k1,S2k2,S2k3;
            Sk1 = ApplyTransformation(m_x_u[2], m_x_u[1], m_x_u[0]);
            Sk2 = ApplyTransformation(m_x_u[3], m_x_u[2], m_x_u[1]);
            Sk3 = ApplyTransformation(m_x_u[4], m_x_u[3], m_x_u[2]);
            S2k1 = ApplyTransformation(Sk3,Sk2,Sk1);
            
            Sk1 = ApplyTransformation(m_x_u[3], m_x_u[2], m_x_u[1]);
            Sk2 = ApplyTransformation(m_x_u[4], m_x_u[3], m_x_u[2]);
            Sk3 = ApplyTransformation(m_x_u[5], m_x_u[4], m_x_u[3]);
            S2k2 = ApplyTransformation(Sk3,Sk2,Sk1);
            
            Sk1 = ApplyTransformation(m_x_u[4], m_x_u[3], m_x_u[2]);
            Sk2 = ApplyTransformation(m_x_u[5], m_x_u[4], m_x_u[3]);
            Sk3 = ApplyTransformation(m_x_u[6], m_x_u[5], m_x_u[4]);
            S2k3 = ApplyTransformation(Sk3,Sk2,Sk1);
            
            m_geomechanic_analysis->Solution() = ApplyTransformation(S2k3,S2k2,S2k1);
            
        }
            break;
        default:
            break;
    }
    
}

void TPMRSSegregatedAnalysis::AccelerationRes(int k, int n){
    
    k--;
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
                m_x_p[0] = m_reservoir_analysis->X_n();
            }
            break;
        case 1:
        {
            m_x_p.Resize(2);
            m_x_p[1] = m_reservoir_analysis->X_n();
            m_reservoir_analysis->X_n() = ApplyTransformation(m_x_p[1], m_x_p[1], m_x_p[0]);
            
        }
            break;
        case 2:  /// S(A_n)
        {
            m_x_p.Resize(3);
            m_x_p[2] = m_reservoir_analysis->X_n();
            m_reservoir_analysis->X_n() = ApplyTransformation(m_x_p[2], m_x_p[1], m_x_p[0]);

        }
            break;
        case 3:  /// S(A_n)
        {
            
            m_x_p.Resize(4);
            m_x_p[3] = m_reservoir_analysis->X_n();
            m_reservoir_analysis->X_n() = ApplyTransformation(m_x_p[3], m_x_p[2], m_x_p[1]);

        }
            break;
        case 4:  /// S2(A_n)
        {
            m_x_p.Resize(5);
            m_x_p[4] = m_reservoir_analysis->X_n();
            
            TPZFMatrix<REAL> Sk1,Sk2,Sk3;
            Sk1 = ApplyTransformation(m_x_p[2], m_x_p[1], m_x_p[0]);
            Sk2 = ApplyTransformation(m_x_p[3], m_x_p[2], m_x_p[1]);
            Sk3 = ApplyTransformation(m_x_p[4], m_x_p[3], m_x_p[2]);
            m_reservoir_analysis->X_n() = ApplyTransformation(Sk3,Sk2,Sk1);
            
        }
            break;
        case 5:  /// S2andhalf(A_n)
        {
            
            m_x_p.Resize(6);
            m_x_p[5] = m_reservoir_analysis->X_n();
            
            TPZFMatrix<REAL> Sk1,Sk2,Sk3;
            Sk1 = ApplyTransformation(m_x_p[3], m_x_p[2], m_x_p[1]);
            Sk2 = ApplyTransformation(m_x_p[4], m_x_p[3], m_x_p[2]);
            Sk3 = ApplyTransformation(m_x_p[5], m_x_p[4], m_x_p[3]);
            m_reservoir_analysis->X_n() = ApplyTransformation(Sk3,Sk2,Sk1);
            
        }
            break;
        case 6:  /// S3(A_n)
        {
            m_x_p.Resize(7);
            m_x_p[6] = m_reservoir_analysis->X_n();
            
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
            
            m_reservoir_analysis->X_n() = ApplyTransformation(S2k3,S2k2,S2k1);
            
        }
            break;
        default:
            break;
    }
    
#ifdef Acceleration_P_output_Q
    for (auto i : m_x_p) {
        i.Print("pi = ",std::cout,EMathematicaInput);
    }
    m_reservoir_analysis->X_n().Print("pstar = ",std::cout,EMathematicaInput);
#endif
    
}

//#define Adaptive_Acceleration_Q

TPZFMatrix<REAL> TPMRSSegregatedAnalysis::ApplyTransformation(TPZFMatrix<REAL> & An_p_1, TPZFMatrix<REAL> & An, TPZFMatrix<REAL> & An_m_1){
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

TPZFMatrix<REAL> TPMRSSegregatedAnalysis::FDMTransformation(TPZFMatrix<REAL> & An_p_1, TPZFMatrix<REAL> & An, TPZFMatrix<REAL> & An_m_1){
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

TPZFMatrix<REAL> TPMRSSegregatedAnalysis::SDMTransformation(TPZFMatrix<REAL> & An_p_1, TPZFMatrix<REAL> & An, TPZFMatrix<REAL> & An_m_1){
    
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

void TPMRSSegregatedAnalysis::PostProcessTimeStep(std::string & geo_file, std::string & res_file){
    m_reservoir_analysis->PostProcessTimeStep(res_file);
    m_geomechanic_analysis->PostProcessTimeStep(geo_file);
}


//#define Animated_Convergence_Q
//#define Noisy_Q


void TPMRSSegregatedAnalysis::ExecuteTimeEvolution(){
    
    // Resize the summaries matrices
    ConfigurateHistorySummaries();
    
    /// vtk files
    std::string name = m_simulation_data->name_vtk_file();
    std::string file_geo = name + "_geo.vtk";
    std::string file_res = name + "_res.vtk";
    
    int n_max_fss_iterations = m_simulation_data->n_fss_iterations();
    int n_enforced_fss_iterations = m_simulation_data->n_enf_fss_iterations();
    int n_time_steps = m_simulation_data->n_steps();
    REAL r_norm = m_simulation_data->epsilon_res();
    REAL dx_norm = m_simulation_data->epsilon_cor();
    REAL dt = m_simulation_data->dt();
    REAL time_value = 0.0;
    
#ifdef Noisy_Q
    m_reservoir_analysis->Solution().Print("psol = ",std::cout,EMathematicaInput);
    m_reservoir_analysis->X().Print("p = ",std::cout,EMathematicaInput);
    m_reservoir_analysis->X_n().Print("pn = ",std::cout,EMathematicaInput);
    m_geomechanic_analysis->Solution().Print("du = ",std::cout,EMathematicaInput);
#endif
    
    std::cout << std::endl;
    std::cout << "-------------------------------------------------------------" <<std::endl;
    std::cout << "-------------------------------------------------------------" <<std::endl;
    std::cout << "TPMRSSegregatedAnalysis:: Begining for recurrent simulation process." <<std::endl;
    std::cout << std::endl << std::endl;
    
    int64_t n_dof_geo = m_geomechanic_analysis->Solution().Rows();
    int64_t n_dof_res = m_reservoir_analysis->X_n().Rows();
    std::cout << "TPMRSSegregatedAnalysis:: Number of DoF for geomechanics module   = " << n_dof_geo <<std::endl;
    std::cout << "TPMRSSegregatedAnalysis:: Number of DoF for reservoir module      = " << n_dof_res <<std::endl;
    
    m_p_m = m_reservoir_analysis->X_n();
    m_u_m = m_geomechanic_analysis->Solution();
    
    bool error_stop_criterion_Q = false;
    bool dx_stop_criterion_Q = false;
    for (int it = 0; it < n_time_steps; it++) {
        /// Interpolate BC data
        time_value = dt * (it+1);
        ConfigureGeomechanicsBC(time_value);
        ConfigureReservoirBC(time_value);
        
        /// Reseting initial du == 0
        m_geomechanic_analysis->Solution().Zero();
        
        for (int k = 1; k <= n_max_fss_iterations; k++) {

#ifdef USING_BOOST
            boost::posix_time::ptime fss_t1 = boost::posix_time::microsec_clock::local_time();
#endif
            
            this->ExecuteOneTimeStep(it,k);
            
#ifdef USING_BOOST
            boost::posix_time::ptime fss_t2 = boost::posix_time::microsec_clock::local_time();
#endif
            
#ifdef USING_BOOST
            
            REAL fss_solving_time = boost::numeric_cast<double>((fss_t2-fss_t1).total_milliseconds());
            std::cout << "TPMRSSegregatedAnalysis:: Fixed stress process closed in :" << setw(10) <<  fss_solving_time/1000.0 << setw(5)   << " seconds." << std::endl;
            std::cout << std::endl;
            
            m_cpu_time_summary(it,0) = time_value;
            m_cpu_time_summary(it,3) += fss_solving_time;
            
#endif
            m_iterations_summary(it,0) = time_value;
            m_iterations_summary(it,1) += m_reservoir_analysis->Get_k_iterations();
            m_iterations_summary(it,2) += m_geomechanic_analysis->Get_k_iterations();
            m_iterations_summary(it,3) = k;
            
            m_residuals_summary(it,0) = time_value;
            m_residuals_summary(it,1) += m_reservoir_analysis->Get_error();
            m_residuals_summary(it,2) += m_geomechanic_analysis->Get_error();
            
            /// Relative error http://mathworld.wolfram.com/RelativeError.html
            REAL p_error = Norm(m_reservoir_analysis->X_n() - m_p_m);
            REAL du_error = Norm(m_geomechanic_analysis->Solution() - m_u_m);
            REAL last_p_norm = Norm(m_p_m);
            REAL last_du_norm = Norm(m_u_m);
            m_fss_p_norm  = p_error/last_p_norm;
            m_fss_du_norm = du_error/last_du_norm;
    
            if (IsZero(p_error) && IsZero(last_p_norm)) {
                m_fss_p_norm = p_error;
            }
            
            if (IsZero(du_error) && IsZero(last_du_norm)) {
                m_fss_du_norm = du_error;
            }
            
            m_p_m = m_reservoir_analysis->X_n();
            m_u_m = m_geomechanic_analysis->Solution();
            
            
            m_residuals_summary(it,3) = m_fss_p_norm;
            m_residuals_summary(it,4) = m_fss_du_norm;
            
            std::cout << "fss_dp_norm = " << m_fss_p_norm << std::endl;
            std::cout << "fss_du_norm = " << m_fss_du_norm << std::endl;
            
            error_stop_criterion_Q = (m_reservoir_analysis->Get_error() < r_norm) && (m_geomechanic_analysis->Get_error() < r_norm);
            dx_stop_criterion_Q = (m_fss_p_norm < dx_norm) && (m_fss_du_norm < dx_norm);
            
#ifdef Animated_Convergence_Q
            this->PostProcessTimeStep(file_geo, file_res);
#endif
            if ((error_stop_criterion_Q && (k > n_enforced_fss_iterations)) && dx_stop_criterion_Q) {
                std::cout << "TPMRSSegregatedAnalysis:: Iterative process converged with fss iterations = " << k << std::endl;
                std::cout << "TPMRSSegregatedAnalysis:: Iterative process converged with residue norm for res = " << m_reservoir_analysis->Get_error() << std::endl;
                std::cout << "TPMRSSegregatedAnalysis:: Iterative process converged with residue norm for geo = " << m_geomechanic_analysis->Get_error() << std::endl;
                std::cout << "-----------------------------------------------------------------------------------------" << std::endl;
                std::cout << std::endl;
                std::cout << std::endl;
#ifdef Noisy_Q
                m_reservoir_analysis->Solution().Print("psol = ",std::cout,EMathematicaInput);
                m_reservoir_analysis->X().Print("p = ",std::cout,EMathematicaInput);
                m_reservoir_analysis->X_n().Print("pn = ",std::cout,EMathematicaInput);
                m_geomechanic_analysis->Solution().Print("du = ",std::cout,EMathematicaInput);
#endif
                break;
            }
        }
        
        bool postprocess_time_Q = ShouldPostprocessQ(time_value);
        if (postprocess_time_Q) {
            this->PostProcessTimeStep(file_geo, file_res);
        }
        UpdateState();
        
#ifdef Noisy_Q
        m_reservoir_analysis->Solution().Print("psol = ",std::cout,EMathematicaInput);
        m_reservoir_analysis->X().Print("p = ",std::cout,EMathematicaInput);
        m_reservoir_analysis->X_n().Print("pn = ",std::cout,EMathematicaInput);
        m_geomechanic_analysis->Solution().Print("du = ",std::cout,EMathematicaInput);
#endif

    }
    
    std::cout << std::endl;
    std::cout << "-------------------------------------------------------------" <<std::endl;
    std::cout << "-------------------------------------------------------------" <<std::endl;
    std::cout << "TPMRSSegregatedAnalysis:: Ending for recurrent simulation process." <<std::endl;
    std::cout << std::endl;
    std::cout << std::endl << std::endl;
    
}

void TPMRSSegregatedAnalysis::UpdateState(){
    m_reservoir_analysis->UpdateState();
    m_geomechanic_analysis->UpdateState();
    
}

bool TPMRSSegregatedAnalysis::ShouldPostprocessQ(REAL time){
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

void TPMRSSegregatedAnalysis::ConfigureGeomechanicsBC(REAL t, bool IsInitialConditionsQ){
    
    TPZCompMesh * cmesh = m_geomechanic_analysis->Mesh();
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
        int matid = material_ids[iregion].first;
        
        TPZMaterial  * material = cmesh->FindMaterial(matid);
        
        /// Update the elastic response
        {
            std::tuple<TPMRSUndrainedParameters, TPMRSPoroMechParameters, TPMRSPhiParameters,TPMRSKappaParameters,TPMRSPlasticityParameters> chunk =    m_simulation_data->MaterialProps()[iregion];
            
            REAL E,nu;
            
            /// Elastic predictor
            if (IsInitialConditionsQ)
            {
                TPMRSUndrainedParameters undrained_parameters(std::get<0>(chunk));
                std::vector<REAL> u_pars = undrained_parameters.GetParameters();
                E = u_pars[0];
                nu = u_pars[1];
            }else
            {
                TPMRSPoroMechParameters poroperm_parameters(std::get<1>(chunk));
                std::vector<REAL> poroperm_pars = poroperm_parameters.GetParameters();
                E = poroperm_pars[0];
                nu = poroperm_pars[1];
            }

            
            TPZElasticResponse ER;
            ER.SetEngineeringData(E, nu);
            
            /// Plastic corrector
            TPMRSPlasticityParameters plasticity_parameters(std::get<4>(chunk));
            std::vector<REAL> p_pars = plasticity_parameters.GetParameters();
            
            if (p_pars.size() == 0)
            {
                /// Elastic material
                TPZElasticCriterion Elastic;
                Elastic.SetElasticResponse(ER);
                
                TPMRSElastoPlastic<TPZElasticCriterion, TPMRSMemory> * vol_material = dynamic_cast< TPMRSElastoPlastic<TPZElasticCriterion, TPMRSMemory> * >(material);
                vol_material->SetPlasticIntegrator(Elastic);
                
            }else
            {
                /// Elastoplastic material
                switch (plasticity_parameters.GetModel())
                {
                    case plasticity_parameters.ep_mc:
                    {
                        
                        /// Mohr Coulomb data
                        REAL cohesion    = p_pars[0];
                        REAL phi         = (p_pars[1]*M_PI/180);
                        REAL psi         = phi;
                        
                        TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC;
                        LEMC.SetElasticResponse(ER);
                        LEMC.fYC.SetUp(phi, psi, cohesion, ER);
                        
                        TPMRSElastoPlastic <TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPMRSMemory> * vol_material = dynamic_cast< TPMRSElastoPlastic <TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPMRSMemory> * >(material);
                        vol_material->SetPlasticIntegrator(LEMC);
                    }
                        break;
                    case plasticity_parameters.ep_ds: {
                        // Dimaggio Sandler data
                        
                        STATE G   = E / (2.0 * (1.0 + nu));
                        STATE K   = E / (3.0 * (1.0 - 2 * nu));
                        
                        REAL A    = p_pars[0];
                        REAL B    = p_pars[1];
                        REAL C    = p_pars[2];
                        REAL D    = p_pars[3];
                        REAL R    = p_pars[4];
                        REAL W    = p_pars[5];
                        REAL X0   = p_pars[6];
                        REAL phi = 0, psi = 1.0, N = 0;
                        
                        REAL Pc = X0/3.0;
                        TPZTensor<REAL> sigma;
                        sigma.Zero();
        
                        sigma.XX() = Pc;
                        sigma.YY() = Pc;
                        sigma.ZZ() = Pc;
                        
                        
                        TPZPlasticStepPV<TPZSandlerExtended, TPZElasticResponse> LEDS;
                        LEDS.SetElasticResponse(ER);
                        LEDS.fYC.SetUp(A, B, C, D, K, G, W, R, phi, N, psi);
                        
                        
                        // Initial damage data
                        REAL k_0;
                        LEDS.InitialDamage(sigma, k_0);
                        LEDS.fN.m_hardening = k_0;
                        LEDS.fYC.SetInitialDamage(k_0);
                        LEDS.fN.m_eps_t.Zero();
                        LEDS.fN.m_eps_p.Zero();
        
                        
                        TPMRSElastoPlastic <TPZPlasticStepPV<TPZSandlerExtended, TPZElasticResponse>, TPMRSMemory> * vol_material = new TPMRSElastoPlastic <TPZPlasticStepPV<TPZSandlerExtended, TPZElasticResponse>, TPMRSMemory>(matid);
                        
                        vol_material->SetPlasticIntegrator(LEDS);
                    
                    }
                        break;

                    default:
                    {
                        std::cout << "Material not implemented. " << std::endl;
                        DebugStop();
                    }
                        break;
                }
            }
        }
        
        /// Inserting boundary conditions
        int n_bc = material_ids[iregion].second.first.size();
        for (int ibc = 0; ibc < n_bc; ibc++)
        {
            int bc_id = material_ids[iregion].second.first [ibc];
            
            if (IsInitialConditionsQ) {
                it_bc_id_to_type = m_simulation_data->BCIdToConditionTypeGeomechanicsInitial().find(bc_id);
                it_bc_id_to_values = m_simulation_data->BCIdToBCValuesGeomechanicsInitial().find(bc_id);
            }else{
                it_bc_id_to_type = m_simulation_data->BCIdToConditionTypeGeomechanics().find(bc_id);
                it_bc_id_to_values = m_simulation_data->BCIdToBCValuesGeomechanics().find(bc_id);
            }
            
            
            it_condition_type_to_index_value_names = m_simulation_data->ConditionTypeToBCIndexGeomechanics().find(it_bc_id_to_type->second);
            
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

void TPMRSSegregatedAnalysis::ConfigureReservoirBC(REAL t, bool IsInitialConditionsQ){
 
    TPZCompMesh * cmesh = m_reservoir_analysis->Mesh();
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
            
            if (IsInitialConditionsQ) {
                it_bc_id_to_type = m_simulation_data->BCIdToConditionTypeReservoirInitial().find(bc_id);
                it_bc_id_to_values = m_simulation_data->BCIdToBCValuesReservoirInitial().find(bc_id);
            }else{
                it_bc_id_to_type = m_simulation_data->BCIdToConditionTypeReservoir().find(bc_id);
                it_bc_id_to_values = m_simulation_data->BCIdToBCValuesReservoir().find(bc_id);
            }
            
            
            it_condition_type_to_index_value_names = m_simulation_data->ConditionTypeToBCIndexReservoir().find(it_bc_id_to_type->second);
            
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
            std::vector<REAL> f_values = it_bc_id_to_values->second.f(t);
            for (int i = 0; i < n_bc_values; i++) {
                REAL value = f_values[i];
                bc->Val2()(i,0) = value;
            }
        }
    }
    
}


void TPMRSSegregatedAnalysis::ExecuteStaticSolution(){
    std::cout << std::endl;
    std::cout << "-------------------------------------------------------------" <<std::endl;
    std::cout << "-------------------------------------------------------------" <<std::endl;
    std::cout << "TPMRSSegregatedAnalysis:: Opening for initialization process." <<std::endl;
    
    int64_t n_dof_geo = m_geomechanic_analysis->Solution().Rows();
    int64_t n_dof_res = m_reservoir_analysis->X_n().Rows();
    std::cout << "TPMRSSegregatedAnalysis:: Number of DoF for geomechanics module   = " << n_dof_geo <<std::endl;
    std::cout << "TPMRSSegregatedAnalysis:: Number of DoF for reservoir module      = " << n_dof_res <<std::endl;
    
    std::string name = m_simulation_data->name_vtk_file();
    std::string file_geo = name + "_geo.vtk";
    std::string file_res = name + "_res.vtk";
    
    REAL dt = m_simulation_data->dt();
    REAL dt_large = 1.0e10;
    m_simulation_data->Setdt(dt_large);
    
    /// Initial reservoir state
    m_reservoir_analysis->ExecuteOneTimeStep();
    m_reservoir_analysis->UpdateState();
    
    /// Compute stress state corresponding to reservoir state
    m_geomechanic_analysis->ExecuteOneTimeStep(true);
    m_geomechanic_analysis->UpdateState();
    m_simulation_data->SetTransferCurrentToLastQ(true);
    m_geomechanic_analysis->UpdateState();
    m_simulation_data->SetTransferCurrentToLastQ(false);
    
    
    UpdateInitialSigmaAndPressure();
    if (m_simulation_data->Get_is_draw_initial_data_Q()) {
        m_reservoir_analysis->PostProcessTimeStep(file_res);
        m_geomechanic_analysis->PostProcessTimeStep(file_geo);
    }
    m_simulation_data->Setdt(dt);
    
    /// Clean state variables
    m_geomechanic_analysis->Solution().Zero();
    m_geomechanic_analysis->X_n().Zero();
    std::cout << "TPMRSSegregatedAnalysis:: Ending for initialization process." <<std::endl;
    std::cout << "-------------------------------------------------------------" <<std::endl;
    std::cout << "-------------------------------------------------------------" <<std::endl;
    std::cout << std::endl;
    std::cout << std::endl << std::endl;
}

void TPMRSSegregatedAnalysis::ExecuteUndrainedStaticSolution()
{
    std::string name = m_simulation_data->name_vtk_file();
    std::string file_geo = name + "_geo.vtk";
    std::string file_res = name + "_res.vtk";
    
    /// Interpolate BC data
    REAL time_value = 0.0;
    ConfigureGeomechanicsBC(time_value);
    ConfigureReservoirBC(time_value);
    std::cout << std::endl;
    std::cout << "-------------------------------------------------------------" <<std::endl;
    std::cout << "-------------------------------------------------------------" <<std::endl;
    std::cout << "TPMRSSegregatedAnalysis:: Begining for undarined response process." <<std::endl;
    std::cout << std::endl << std::endl;
    { /// Computing the undrained response t_0+
        m_reservoir_analysis->ExecuteOneTimeStep();
        m_reservoir_analysis->UpdateState();
        m_simulation_data->SetTransferCurrentToLastQ(true);
        m_reservoir_analysis->UpdateState();
        m_simulation_data->SetTransferCurrentToLastQ(false);
        
        m_geomechanic_analysis->ExecuteOneTimeStep(true);
        m_geomechanic_analysis->UpdateState();
        m_simulation_data->SetTransferCurrentToLastQ(true);
        m_geomechanic_analysis->UpdateState();
        m_simulation_data->SetTransferCurrentToLastQ(false);
        bool reset_undrained_data_Q = m_simulation_data->Get_reset_undarined_respose_data_Q();
        UpdateUndrainedSigmaAndPressureData(reset_undrained_data_Q,false);
        if (m_simulation_data->Get_is_draw_initial_data_Q()) {
            PostProcessTimeStep(file_geo, file_res);
        }
    }
    std::cout << "TPMRSSegregatedAnalysis:: Ending for undarined response process." <<std::endl;
    std::cout << "-------------------------------------------------------------" <<std::endl;
    std::cout << "-------------------------------------------------------------" <<std::endl;
    std::cout << std::endl;
    std::cout << std::endl << std::endl;
}

void TPMRSSegregatedAnalysis::UpdateInitialSigmaAndPressure(bool reset_u_Q, bool keep_phi_Q) {
    
    TPZCompMesh * cmesh = m_geomechanic_analysis->Mesh();
    
    if (!m_simulation_data || !cmesh) {
        DebugStop();
    }
    
    TPZManVector<REAL,3> u_null(3,0.0);
    int n_regions = m_simulation_data->NumberOfRegions();
    TPZManVector<std::pair<int, std::pair<TPZManVector<int,12>,TPZManVector<int,12>> >,12>  material_ids = m_simulation_data->MaterialIds();
    TPZManVector<int,10> volumetric_mat_id(n_regions);
    for (int iregion = 0; iregion < n_regions; iregion++)
    {
        int matid = material_ids[iregion].first;
        volumetric_mat_id[iregion] = matid;
        
        TPZMaterial * material = cmesh->FindMaterial(matid);
        if (!material) {
            DebugStop();
        }
        
        TPZMatWithMem<TPMRSMemory> * mat_with_memory = dynamic_cast<TPZMatWithMem<TPMRSMemory> * >(material);
        if (!material) {
            DebugStop();
        }
        
        std::tuple<TPMRSUndrainedParameters, TPMRSPoroMechParameters, TPMRSPhiParameters,TPMRSKappaParameters,TPMRSPlasticityParameters> chunk =    m_simulation_data->MaterialProps()[iregion];
        
        TPMRSPoroMechParameters poro_parameters(std::get<1>(chunk));
        std::vector<REAL> poroperm_pars = poro_parameters.GetParameters();
        
        std::shared_ptr<TPZAdmChunkVector<TPMRSMemory>> & memory_vector = mat_with_memory->GetMemory();
        
        int ndata = memory_vector->NElements();
        for (int i = 0; i < ndata; i++) {
            

            REAL p_0 = memory_vector.get()->operator [](i).p_n();
            memory_vector.get()->operator [](i).Setp_0(p_0);
            memory_vector.get()->operator [](i).Setp(p_0);
            
            REAL f_0 = 0.0;
            TPZManVector<REAL,3> f_vec_0(3);
            f_vec_0[0] = f_0;
            f_vec_0[1] = f_0;
            f_vec_0[2] = f_0;
            memory_vector.get()->operator [](i).Setf(f_0);
            memory_vector.get()->operator [](i).Setf_vec(f_vec_0);
            
            REAL geo_delta_phi_0 = 0.0;
            memory_vector.get()->operator [](i).Setdelta_phi(geo_delta_phi_0);
            
            TPZPlasticState<REAL> plas_state_0 = memory_vector.get()->operator [](i).GetPlasticState_n();
            memory_vector.get()->operator [](i).SetPlasticState_0(plas_state_0);
            memory_vector.get()->operator [](i).SetPlasticState(plas_state_0);
            memory_vector.get()->operator [](i).SetPlasticStateSubStep(plas_state_0);
            
            TPZTensor<REAL> sigma_0 = memory_vector.get()->operator [](i).GetSigma_n();
            memory_vector.get()->operator [](i).SetSigma_0(sigma_0);
            memory_vector.get()->operator [](i).SetSigma(sigma_0);
            
            if(keep_phi_Q){
                REAL phi_0 = memory_vector.get()->operator [](i).phi_0();
                memory_vector.get()->operator [](i).Setphi_n(phi_0);
                
                REAL kappa_0 = memory_vector.get()->operator [](i).kappa_0();
                memory_vector.get()->operator [](i).Setkappa_n(kappa_0);
                
            }else{
                REAL phi_0 = memory_vector.get()->operator [](i).phi_n();
                memory_vector.get()->operator [](i).Setphi_0(phi_0);
                
                REAL kappa_0 = memory_vector.get()->operator [](i).kappa_n();
                memory_vector.get()->operator [](i).Setkappa_0(kappa_0);
                
            }
        
            
            /// Cleaning u and plastic state
            if (reset_u_Q) {
                
                memory_vector.get()->operator [](i).Setu_0(u_null);
                memory_vector.get()->operator [](i).Setu(u_null);
                memory_vector.get()->operator [](i).Setu_n(u_null);
                memory_vector.get()->operator [](i).Setu_sub_step(u_null);
                
                memory_vector.get()->operator [](i).GetPlasticState_0().m_eps_p.Zero();
                memory_vector.get()->operator [](i).GetPlasticState_n().m_eps_p.Zero();
                memory_vector.get()->operator [](i).GetPlasticState().m_eps_p.Zero();
                memory_vector.get()->operator [](i).GetPlasticStateSubStep().m_eps_p.Zero();
                
            }else{
                
                TPZManVector<REAL,3> u_0 = memory_vector.get()->operator [](i).Getu_n();
                memory_vector.get()->operator [](i).Setu_0(u_0);
                memory_vector.get()->operator [](i).Setu(u_0);
                memory_vector.get()->operator [](i).Setu_n(u_0);
                memory_vector.get()->operator [](i).Setu_sub_step(u_0);
                
            }
            
        }
        
    }
    
}

void TPMRSSegregatedAnalysis::UpdateUndrainedSigmaAndPressureData(bool reset_u_Q, bool keep_plastic_state_Q) {
    
    TPZCompMesh * cmesh = m_geomechanic_analysis->Mesh();
    
    if (!m_simulation_data || !cmesh) {
        DebugStop();
    }
    
    TPZManVector<REAL,3> u_null(3,0.0);
    int n_regions = m_simulation_data->NumberOfRegions();
    TPZManVector<std::pair<int, std::pair<TPZManVector<int,12>,TPZManVector<int,12>> >,12>  material_ids = m_simulation_data->MaterialIds();
    TPZManVector<int,10> volumetric_mat_id(n_regions);
    for (int iregion = 0; iregion < n_regions; iregion++)
    {
        int matid = material_ids[iregion].first;
        volumetric_mat_id[iregion] = matid;
        
        TPZMaterial * material = cmesh->FindMaterial(matid);
        if (!material) {
            DebugStop();
        }
        
        TPZMatWithMem<TPMRSMemory> * mat_with_memory = dynamic_cast<TPZMatWithMem<TPMRSMemory> * >(material);
        if (!material) {
            DebugStop();
        }
        
        std::tuple<TPMRSUndrainedParameters, TPMRSPoroMechParameters, TPMRSPhiParameters,TPMRSKappaParameters,TPMRSPlasticityParameters> chunk =    m_simulation_data->MaterialProps()[iregion];
        
        TPMRSPoroMechParameters poro_parameters(std::get<1>(chunk));
        std::vector<REAL> poroperm_pars = poro_parameters.GetParameters();
        
        std::shared_ptr<TPZAdmChunkVector<TPMRSMemory>> & memory_vector = mat_with_memory->GetMemory();
        
        int ndata = memory_vector->NElements();
        for (int i = 0; i < ndata; i++) {
            
            
            REAL p = memory_vector.get()->operator [](i).p_n();
            memory_vector.get()->operator [](i).Setp(p);
            
            REAL f_0 = 0.0;
            TPZManVector<REAL,3> f_vec_0(3);
            f_vec_0[0] = f_0;
            f_vec_0[1] = f_0;
            f_vec_0[2] = f_0;
            memory_vector.get()->operator [](i).Setf(f_0);
            memory_vector.get()->operator [](i).Setf_vec(f_vec_0);
            
            REAL geo_delta_phi_0 = 0.0;
            memory_vector.get()->operator [](i).Setdelta_phi(geo_delta_phi_0);
            
            TPZPlasticState<REAL> plas_state = memory_vector.get()->operator [](i).GetPlasticState_n();
            memory_vector.get()->operator [](i).SetPlasticState(plas_state);
            memory_vector.get()->operator [](i).SetPlasticStateSubStep(plas_state);
            
            TPZTensor<REAL> sigma = memory_vector.get()->operator [](i).GetSigma_n();
            memory_vector.get()->operator [](i).SetSigma(sigma);
            
            /// porosity update
            {
                
                REAL alpha = memory_vector.get()->operator [](i).Alpha();
                REAL Kdr   = memory_vector.get()->operator [](i).Kdr();
                REAL phi_0 = memory_vector.get()->operator [](i).phi_0();
                
                REAL p_0   = memory_vector.get()->operator [](i).p_0();
                REAL p_n   = memory_vector.get()->operator [](i).p_n();
                
                REAL S = (1.0-alpha)*(alpha-phi_0)/Kdr;
                
                REAL eps_v_t_0   = memory_vector.get()->operator [](i).GetPlasticState_0().m_eps_t.I1();
                REAL eps_v_t_n   = memory_vector.get()->operator [](i).GetPlasticState_n().m_eps_t.I1();
                
                REAL phi_n = phi_0 + (alpha) * (eps_v_t_n-eps_v_t_0) + S * (p_n - p_0);
                memory_vector.get()->operator [](i).Setphi_n(phi_n);
            }
            
//            if(keep_plastic_state_Q){
////                REAL phi_0 = memory_vector.get()->operator [](i).phi_0();
////                memory_vector.get()->operator [](i).Setphi_n(phi_0);
////
////                REAL kappa_0 = memory_vector.get()->operator [](i).kappa_0();
////                memory_vector.get()->operator [](i).Setkappa_n(kappa_0);
//
//            }else{
//                REAL phi_0 = memory_vector.get()->operator [](i).phi_n();
//                memory_vector.get()->operator [](i).Setphi_0(phi_0);
//
//                REAL kappa_0 = memory_vector.get()->operator [](i).kappa_n();
//                memory_vector.get()->operator [](i).Setkappa_0(kappa_0);
//
//            }
            
            
            /// Cleaning u and plastic state
            if (reset_u_Q) {
                
                memory_vector.get()->operator [](i).Setu_0(u_null);
                memory_vector.get()->operator [](i).Setu(u_null);
                memory_vector.get()->operator [](i).Setu_n(u_null);
                memory_vector.get()->operator [](i).Setu_sub_step(u_null);
                
                memory_vector.get()->operator [](i).GetPlasticState_0().m_eps_p.Zero();
                memory_vector.get()->operator [](i).GetPlasticState_n().m_eps_p.Zero();
                memory_vector.get()->operator [](i).GetPlasticState().m_eps_p.Zero();
                memory_vector.get()->operator [](i).GetPlasticStateSubStep().m_eps_p.Zero();
                
            }else{
                
                TPZManVector<REAL,3> u_0 = memory_vector.get()->operator [](i).Getu_n();
                memory_vector.get()->operator [](i).Setu_0(u_0);
                memory_vector.get()->operator [](i).Setu(u_0);
                memory_vector.get()->operator [](i).Setu_n(u_0);
                memory_vector.get()->operator [](i).Setu_sub_step(u_0);
                
            }
            
        }
        
    }
    
}

void TPMRSSegregatedAnalysis::SetSimulationData(TPMRSSimulationData * simulation_data)
{
    m_simulation_data = simulation_data;
}

void TPMRSSegregatedAnalysis::ConfigurateHistorySummaries(){
    int n_time_steps = m_simulation_data->n_steps();
    m_iterations_summary.Resize(n_time_steps, 4); // (time,res_iteraions,geo_iterations,fss_iteraions)
    m_cpu_time_summary.Resize(n_time_steps, 4); // (time,res_cpu_time,geo_cpu_time,fss_cpu_time)
    m_residuals_summary.Resize(n_time_steps, 5); // (time,res_resdials,geo_resdials,fss_corrections)
    
    m_iterations_summary.Zero();
    m_cpu_time_summary.Zero();
    m_residuals_summary.Zero();
}

TPZFMatrix<REAL> & TPMRSSegregatedAnalysis::IterationsSummary(){
    return m_iterations_summary;
}

TPZFMatrix<REAL> & TPMRSSegregatedAnalysis::TimeSummary(){
    return m_cpu_time_summary;
}

TPZFMatrix<REAL> & TPMRSSegregatedAnalysis::ResidualsSummary(){
    return m_residuals_summary;
}

