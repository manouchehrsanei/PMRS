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
    bool mustOptimizeBandwidth = false;

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
        REAL alpha   = poroperm_pars[2];
        REAL Kdr     = Eyoung/(3.0*(1.0-2.0*nu));
        
        std::shared_ptr<TPZAdmChunkVector<TPMRSMemory>> & memory_vector = mat_with_memory->GetMemory();
        
        int ndata = memory_vector->NElements();
        for (int i = 0; i < ndata; i++) {
            memory_vector.get()->operator [](i).Setphi_0(phi_0);
            memory_vector.get()->operator [](i).Setkappa_0(kappa_0);
            memory_vector.get()->operator [](i).SetAlpha(alpha);
            memory_vector.get()->operator [](i).SetKdr(Kdr);
        }
        
    }
}


void TPMRSSegregatedAnalysis::AdjustIntegrationOrder(TPZCompMesh * cmesh_o, TPZCompMesh * cmesh_d){
    

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

//#define QNAcceleration_Q
//#define AitkenAcceleration_Q

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
    
#ifdef AitkenAcceleration_Q
    AitkenAccelerationRes(k);
#endif
    
#ifdef QNAcceleration_Q
    QNAccelerationRes(k);
#endif

    
#ifdef USING_BOOST
    boost::posix_time::ptime geo_t1 = boost::posix_time::microsec_clock::local_time();
#endif
    
    m_geomechanic_analysis->ExecuteOneTimeStep();
    
#ifdef USING_BOOST
    boost::posix_time::ptime geo_t2 = boost::posix_time::microsec_clock::local_time();
#endif
    
#ifdef USING_BOOST
    REAL geo_solving_time = boost::numeric_cast<double>((geo_t2-geo_t1).total_milliseconds());
    std::cout << "TPMRSGeomechanicAnalysis:: Newton process closed in :" << setw(10) <<  geo_solving_time/1000.0 << setw(5)   << " seconds." << std::endl;
    std::cout << std::endl;
    
    m_cpu_time_summary(i_time_step,2) += geo_solving_time;
#endif
    
#ifdef AitkenAcceleration_Q
    AitkenAccelerationGeo(k);
#endif
    
#ifdef QNAcceleration_Q
    QNAccelerationGeo(k);
#endif
    
}

void TPMRSSegregatedAnalysis::QNAccelerationRes(int k){
    if (k>1) {
        m_xp_m = m_reservoir_analysis->Rhs();
        m_reservoir_analysis->Rhs() = m_xp_m - m_xp_m_1;
        m_reservoir_analysis->Solve();
        TPZFMatrix<REAL> dp = m_reservoir_analysis->Solution();
//                m_reservoir_analysis->X_n().Print("p = ", std::cout);
//                dp.Print("dp = ", std::cout);
        m_reservoir_analysis->Solution() = m_reservoir_analysis->X_n() + dp;
//                m_reservoir_analysis->Solution().Print("p = ", std::cout);
        m_reservoir_analysis->LoadMemorySolution();
        m_xp_m_1 = m_reservoir_analysis->Rhs();
        
    }else{
        m_xp_m_1 = m_reservoir_analysis->Rhs();
    }
}

void TPMRSSegregatedAnalysis::QNAccelerationGeo(int k){
    if (k>1) {
        TPZFMatrix<REAL> du = m_geomechanic_analysis->Solution();
        m_xu_m = m_geomechanic_analysis->Rhs();
        m_geomechanic_analysis->Rhs() = m_xu_m - m_xu_m_1;
        m_geomechanic_analysis->Solve();
        TPZFMatrix<REAL> delta_du = m_geomechanic_analysis->Solution();
        m_geomechanic_analysis->Solution() = du + delta_du;
        m_geomechanic_analysis->LoadSolution(m_geomechanic_analysis->Solution());
        m_geomechanic_analysis->LoadMemorySolution();
        m_xu_m_1 = m_geomechanic_analysis->Rhs();
        
    }else{
        m_xu_m_1 = m_geomechanic_analysis->Rhs();
    }
}

void TPMRSSegregatedAnalysis::AitkenAccelerationRes(int k){
    
    if (k>2) {
        m_xp_m = m_reservoir_analysis->X_n();
        int n_dof = m_xp_m.Rows();
        REAL denom = 0.0;
        REAL numer = 0.0;
        for (int i = 0; i < n_dof; i++) {
            numer += (m_xp_m_2(i,0)-m_xp_m_1(i,0))*(m_xp_m_1(i,0)-m_xp_m(i,0));
            denom += (m_xp_m_2(i,0)-m_xp_m_1(i,0))*(m_xp_m_2(i,0)-2.0*m_xp_m_1(i,0)+m_xp_m(i,0));
        }
        //        m_xp_m.Print("p = ", std::cout);
        REAL s;
        if(IsZero(denom)){
            s = numer / denom;
        }else{
            s = numer / denom;
        }
        m_reservoir_analysis->X_n() = m_xp_m + s*(m_xp_m-m_xp_m_1);
        m_reservoir_analysis->LoadMemorySolution();
        m_xp_m_2 = m_xp_m_1;
        m_xp_m_1 = m_xp_m;
        //        m_reservoir_analysis->X_n().Print("p new = ", std::cout);
    }else if(k>1){
        m_xp_m_1 = m_reservoir_analysis->X_n();
    }else{
        m_xp_m_2 = m_reservoir_analysis->X_n();
    }
    
}

void TPMRSSegregatedAnalysis::AitkenAccelerationGeo(int k){
    
    if (k>2) {
        m_xu_m = m_geomechanic_analysis->Solution();
        int n_dof = m_xp_m.Rows();
        REAL denom = 0.0;
        REAL numer = 0.0;
        for (int i = 0; i < n_dof; i++) {
            numer += (m_xu_m_2(i,0)-m_xu_m_1(i,0))*(m_xu_m_1(i,0)-m_xu_m(i,0));
            denom += (m_xu_m_2(i,0)-m_xu_m_1(i,0))*(m_xu_m_2(i,0)-2.0*m_xu_m_1(i,0)+m_xu_m(i,0));
        }
//        m_xu_m.Print("u = ", std::cout);
        REAL s;
        if(IsZero(denom)){
            s = numer / denom;
        }else{
            s = numer / denom;
        }
        m_geomechanic_analysis->Solution() = m_xu_m + s*(m_xu_m-m_xu_m_1);
        m_geomechanic_analysis->LoadMemorySolution();
        m_xu_m_2 = m_xu_m_1;
        m_xu_m_1 = m_xu_m;
//        m_geomechanic_analysis->Solution().Print("u new = ", std::cout);
    }else if(k>1){
        m_xu_m_1 = m_geomechanic_analysis->Solution();
    }else{
        m_xu_m_2 = m_geomechanic_analysis->Solution();
    }
}

void TPMRSSegregatedAnalysis::PostProcessTimeStep(std::string & geo_file, std::string & res_file){
    m_reservoir_analysis->PostProcessTimeStep(res_file);
    m_geomechanic_analysis->PostProcessTimeStep(geo_file);
}


//#define EC_Q
//#define Animated_Convergence_Q

void TPMRSSegregatedAnalysis::ExecuteTimeEvolution(){
    
    // Resize the summaries matrices
    ConfigurateHistorySummaries();
    
    /// vtk files
    std::string name = m_simulation_data->name_vtk_file();
    std::string file_geo = name + "_geo.vtk";
    std::string file_res = name + "_res.vtk";
    
    int n_max_fss_iterations = m_simulation_data->n_fss_iterations();
    int n_enforced_fss_iterations = m_simulation_data->n_enf_fss_iterations();
    int n_time_steps = m_simulation_data->ReportingTimes().size();
    REAL r_norm = m_simulation_data->epsilon_res();
    REAL dx_norm = m_simulation_data->epsilon_cor();
    REAL dt = m_simulation_data->dt();
    REAL time_value;
    
#ifdef QNAcceleration_Q
    /// Loading initial data
    m_xp_m = m_xp_m_1 = m_xp_m_2 = m_reservoir_analysis->X_n();
    m_xu_m = m_xu_m_1 = m_xu_m_2 = m_geomechanic_analysis->X_n();
#endif
    
#ifdef EC_Q
    /// Loading initial data
    m_xp_m = m_xp_m_1 = m_xp_m_2 = m_reservoir_analysis->X_n();
    m_xu_m = m_xu_m_1 = m_xu_m_2 = m_geomechanic_analysis->X_n();
#endif
    
    m_p_m = m_reservoir_analysis->X_n();
    m_u_m = m_geomechanic_analysis->Solution();
    
    bool error_stop_criterion_Q = false;
    bool dx_stop_criterion_Q = false;
    for (int it = 0; it < n_time_steps; it++) {
        time_value = dt * (it+1);
        
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
            REAL fss_dp_norm = Norm(m_reservoir_analysis->X_n() - m_p_m);
            REAL fss_du_norm = Norm(m_geomechanic_analysis->Solution() - m_u_m);
            m_p_m = m_reservoir_analysis->X_n();
            m_u_m = m_geomechanic_analysis->Solution();
            m_residuals_summary(it,3) = fss_dp_norm + fss_du_norm;
            
            std::cout << "fss_dp_norm = " << fss_dp_norm << std::endl;
            std::cout << "fss_du_norm = " << fss_du_norm << std::endl;
            
            error_stop_criterion_Q = (m_reservoir_analysis->Get_error() < r_norm) && (m_geomechanic_analysis->Get_error() < r_norm);
            dx_stop_criterion_Q = (fss_dp_norm < dx_norm) && (fss_du_norm < dx_norm);
            
#ifdef Animated_Convergence_Q
            this->PostProcessTimeStep(file_geo, file_res);
#endif
            if ((error_stop_criterion_Q && (k > n_enforced_fss_iterations)) && dx_stop_criterion_Q) {
                std::cout << "TPMRSSegregatedAnalysis:: Iterative process converged with residue norm for res = " << m_reservoir_analysis->Get_error() << std::endl;
                std::cout << "TPMRSSegregatedAnalysis:: Iterative process converged with residue norm for geo = " << m_geomechanic_analysis->Get_error() << std::endl;
                std::cout << "-----------------------------------------------------------------------------------------" << std::endl;
                std::cout << std::endl;
                std::cout << std::endl;
//                m_geomechanic_analysis->AssembleResidual();
                UpdateState();
                this->PostProcessTimeStep(file_geo, file_res);
                
                
#ifdef EC_Q
                /// Enhance pressure
                {
                    REAL s = 0.5;
                    m_xp_m = m_reservoir_analysis->X_n();
                    int n_dof = m_xp_m.Rows();
                    REAL dt = m_simulation_data->dt();
                    REAL f1,f2,f3,t1,t2,t3, t;
                    t = (it+1) * dt + s*dt;
                    t3 = (it+1) * dt;
                    t2 = t3 - dt;
                    t1 = t2 - dt;
                    
                    //        m_xp_m_2.Print("xp_2 = ",std::cout);
                    //        m_xp_m_1.Print("xp_1 = ",std::cout);
                    //        m_xp_m.Print("xp = ",std::cout);
                    for (int i = 0; i < n_dof; i++) {
                        f1 = m_xp_m_2(i,0);
                        f2 = m_xp_m_1(i,0);
                        f3 = m_xp_m(i,0);
                        m_reservoir_analysis->X_n()(i,0) = quadratic_extrapolation(f1, f2, f3, t1, t2, t3, t);
                    }
                    //        m_reservoir_analysis->X_n().Print("p = ",std::cout);
                    m_reservoir_analysis->LoadMemorySolution();
                    m_xp_m_2 = m_xp_m_1;
                    m_xp_m_1 = m_xp_m;
                }
#endif
                
#ifdef EC_Q
                /// Enhance displacement
                {
                    REAL s = 0.5;
                    m_xu_m = m_geomechanic_analysis->X_n();
                    int n_dof = m_xp_m.Rows();
                    REAL dt = m_simulation_data->dt();
                    REAL f1,f2,f3,t1,t2,t3, t;
                    t = (it+1) * dt + s*dt;
                    t3 = (it+1) * dt;
                    t2 = t3 - dt;
                    t1 = t2 - dt;
                    //        m_xu_m_2.Print("xu_2 = ",std::cout);
                    //        m_xu_m_1.Print("xu_1 = ",std::cout);
                    //        m_xu_m.Print("xu = ",std::cout);
                    for (int i = 0; i < n_dof; i++) {
                        f1 = m_xu_m_2(i,0);
                        f2 = m_xu_m_1(i,0);
                        f3 = m_xu_m(i,0);
                        m_geomechanic_analysis->X_n()(i,0) = quadratic_extrapolation(f1, f2, f3, t1, t2, t3, t);
                    }
                    //        m_geomechanic_analysis->X_n().Print("u = ",std::cout);
                    m_geomechanic_analysis->LoadMemorySolution();
                    m_xu_m_2 = m_xu_m_1;
                    m_xu_m_1 = m_xu_m;
                }
#endif
                break;
            }
        }
    }
    
}

void TPMRSSegregatedAnalysis::UpdateState(){
    m_reservoir_analysis->UpdateState();
    m_geomechanic_analysis->UpdateState();
    
}

void TPMRSSegregatedAnalysis::ConfigurateBConditions(bool IsInitialConditionsQ){
    
    TPZCompMesh * cmesh = m_geomechanic_analysis->Mesh();
    if (!cmesh)
    {
        DebugStop();
    }
    
    int n_regions = m_simulation_data->NumberOfRegions();
    TPZManVector<std::pair<int, std::pair<TPZManVector<int,12>,TPZManVector<int,12>> >,12>  material_ids = m_simulation_data->MaterialIds();
    
    std::map<int, std::string>::iterator it_bc_id_to_type;
    std::map< std::string,std::pair<int,std::vector<std::string> > >::iterator  it_condition_type_to_index_value_names;
    std::map<int , std::vector<REAL> >::iterator it_bc_id_to_values;
    
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
            
//            /// Updating bulk modulus for porosity model
//            std::get<2>(m_simulation_data->MaterialProps()[iregion]).SetBulkModulus(E, nu);
            
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
                        
                        REAL Pc = -137;
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
                it_bc_id_to_type = m_simulation_data->BCIdToConditionTypeGeomechanicsUndrained().find(bc_id);
                it_bc_id_to_values = m_simulation_data->BCIdToBCValuesGeomechanicsUndrained().find(bc_id);
            }else{
                it_bc_id_to_type = m_simulation_data->BCIdToConditionTypeGeomechanics().find(bc_id);
                it_bc_id_to_values = m_simulation_data->BCIdToBCValuesGeomechanics().find(bc_id);
            }
            
            
            it_condition_type_to_index_value_names = m_simulation_data->ConditionTypeToBCIndexGeomechanics().find(it_bc_id_to_type->second);
            
            int bc_index = it_condition_type_to_index_value_names->second.first;
            int n_bc_values = it_bc_id_to_values->second.size();
        
            TPZMaterial * bc_mat = cmesh->FindMaterial(bc_id);
            if (!bc_mat) {
                DebugStop();
            }
            TPZBndCond * bc = dynamic_cast<TPZBndCond *>(bc_mat);
            if (!bc) {
                DebugStop();
            }
            
            bc->SetType(bc_index);
            for (int i = 0; i < n_bc_values; i++) {
                REAL value = it_bc_id_to_values->second[i];
                bc->Val2()(i,0) = value;
            }
        }
    }
}

void TPMRSSegregatedAnalysis::ExecuteStaticSolution(){
    std::cout << std::endl;
    std::cout << "TPMRSSegregatedAnalysis:: Opening for initialization process." <<std::endl;
    std::string name = m_simulation_data->name_vtk_file();
    std::string file_geo = name + "_geo.vtk";
    std::string file_res = name + "_res.vtk";
    m_geomechanic_analysis->ExecuteUndrainedResponseStep();
    m_geomechanic_analysis->UpdateState();
    this->UpdateInitialSigmaAndPressure();
    m_reservoir_analysis->ExecuteUndrainedResponseStep();
    m_reservoir_analysis->PostProcessTimeStep(file_res);
    m_geomechanic_analysis->PostProcessTimeStep(file_geo);
    std::cout << "TPMRSSegregatedAnalysis:: Ending for initialization process." <<std::endl;
    std::cout << std::endl << std::endl;
}

void TPMRSSegregatedAnalysis::UpdateInitialSigmaAndPressure() {
    
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
        REAL cf      = poroperm_pars[3];
        
        std::shared_ptr<TPZAdmChunkVector<TPMRSMemory>> & memory_vector = mat_with_memory->GetMemory();
        
        int ndata = memory_vector->NElements();
        for (int i = 0; i < ndata; i++) {
            
            /// Because we reused the same memory items
            TPZTensor<REAL> sigma_total_0 = memory_vector.get()->operator [](i).GetSigma_n();
            REAL alpha = memory_vector.get()->operator [](i).Alpha();
            REAL phi_0 = memory_vector.get()->operator [](i).phi_0();
            REAL Kdr = memory_vector.get()->operator [](i).Kdr();
            REAL Se = ((1.0-alpha)*(alpha - phi_0)/(Kdr)) + phi_0*cf;
            REAL p_0 = -alpha*(sigma_total_0.I1()/3)/(Se * Kdr+alpha*alpha);
            
            memory_vector.get()->operator [](i).Setp_0(p_0);
            memory_vector.get()->operator [](i).Setp(p_0);
            memory_vector.get()->operator [](i).Setp_n(p_0);
            sigma_total_0.Zero();/// Converted to effecttive because initial deformation is Zero.
            memory_vector.get()->operator [](i).SetSigma_0(sigma_total_0);
            memory_vector.get()->operator [](i).SetSigma(sigma_total_0);
            memory_vector.get()->operator [](i).SetSigma_n(sigma_total_0);
            memory_vector.get()->operator [](i).Setdelta_phi(0.0); //  Initial porosity correction is zero.
            /// Cleaning u
            memory_vector.get()->operator [](i).Setu_0(u_null);
            memory_vector.get()->operator [](i).Setu(u_null);
            memory_vector.get()->operator [](i).Setu_n(u_null);
            /// Cleaning elasto-plastic states
            memory_vector.get()->operator [](i).GetPlasticState_0().CleanUp();
            memory_vector.get()->operator [](i).GetPlasticState().CleanUp();
            memory_vector.get()->operator [](i).GetPlasticState_n().CleanUp();
        }
        
    }
    
}

void TPMRSSegregatedAnalysis::SetSimulationData(TPMRSSimulationData * simulation_data)
{
    m_simulation_data = simulation_data;
}

void TPMRSSegregatedAnalysis::ConfigurateHistorySummaries(){
    int n_time_steps = m_simulation_data->ReportingTimes().size();
    m_iterations_summary.Resize(n_time_steps, 4); // (time,res_iteraions,geo_iterations,fss_iteraions)
    m_cpu_time_summary.Resize(n_time_steps, 4); // (time,res_cpu_time,geo_cpu_time,fss_cpu_time)
    m_residuals_summary.Resize(n_time_steps, 4); // (time,res_resdials,geo_resdials,fss_corrections)
    
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


REAL TPMRSSegregatedAnalysis::linear_extrapolation(REAL & f_1, REAL & f_2, REAL & t_1, REAL & t_2, REAL & t){
    REAL f = f_1 + ((f_1-f_2)*(t-t_1))/(t_1 - t_2);
    return f;
}

REAL TPMRSSegregatedAnalysis::quadratic_extrapolation(REAL & f1, REAL & f2, REAL & f3, REAL & t1, REAL & t2,  REAL & t3, REAL & t){
    REAL f = f1 + (t - t1)*((-f1 + f2)/(-t1 + t2) + ((t - t2)*(-((-f1 + f2)/(-t1 + t2)) + (-f2 + f3)/(-t2 + t3)))/(-t1 + t3));
    return f;
}
