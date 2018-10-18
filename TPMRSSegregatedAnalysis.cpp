//
//  TPMRSSegregatedAnalysis.cpp
//  PMRS
//
//  Created by Omar Durán on 9/13/18.
//

#include "TPMRSSegregatedAnalysis.h"


TPMRSSegregatedAnalysis::TPMRSSegregatedAnalysis(){
    m_simulation_data       = NULL;
    m_geomechanic_analysis  = NULL;
    m_reservoir_analysis    = NULL;
}

TPMRSSegregatedAnalysis::~TPMRSSegregatedAnalysis(){
    
}

TPMRSSegregatedAnalysis::TPMRSSegregatedAnalysis(const TPMRSSegregatedAnalysis & other){
    m_simulation_data       = other.m_simulation_data;
    m_geomechanic_analysis  = other.m_geomechanic_analysis;
    m_reservoir_analysis    = other.m_reservoir_analysis;
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

void TPMRSSegregatedAnalysis::ConfigurateAnalysis(DecomposeType decompose_geo, DecomposeType decompose_res, TPZSimulationData * simulation_data, TPZCompMesh * cmesh_geomechanics, TPZCompMesh * cmesh_reservoir, TPZManVector<TPZCompMesh * , 2> & mesh_vec){

    if (!simulation_data || !cmesh_geomechanics || !cmesh_reservoir) {
        DebugStop();
    }
    
    this->SetSimulationData(simulation_data);
    bool mustOptimizeBandwidth = false;

    if (simulation_data->ElasticityOrder()==simulation_data->DiffusionOrder() && simulation_data->Get_is_dual_formulation_Q()) {
        this->ApplyMemoryLink(cmesh_reservoir,cmesh_geomechanics);
        this->AdjustIntegrationOrder(cmesh_reservoir,cmesh_geomechanics);
    }else{
        this->ApplyMemoryLink(cmesh_geomechanics,cmesh_reservoir);
        this->AdjustIntegrationOrder(cmesh_geomechanics,cmesh_reservoir);
    }

    // The Geomechanics Simulator
    m_geomechanic_analysis = new TPMRSGeomechanicAnalysis;
    m_geomechanic_analysis->SetCompMesh(cmesh_geomechanics,mustOptimizeBandwidth);
    m_geomechanic_analysis->ConfigurateAnalysis(decompose_geo, m_simulation_data);
    
    // The Reservoir Simulator
    m_reservoir_analysis = new TPMRSMonoPhasicAnalysis;
    m_reservoir_analysis->SetCompMesh(cmesh_reservoir,mustOptimizeBandwidth);
    m_reservoir_analysis->ConfigurateAnalysis(decompose_res, mesh_vec, m_simulation_data);
    
    // Loading spatial properties
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

        REAL phi_0 = undrained_pars[2];
        REAL kappa_0 = undrained_pars[3];
        REAL alpha = poroperm_pars[2];
        REAL Se = poroperm_pars[3];
        
        std::shared_ptr<TPZAdmChunkVector<TPMRSMemory>> & memory_vector = mat_with_memory->GetMemory();
        
        int ndata = memory_vector->NElements();
        for (int i = 0; i < ndata; i++) {
            memory_vector.get()->operator [](i).Setphi_0(phi_0);
            memory_vector.get()->operator [](i).Setkappa_0(kappa_0);
            memory_vector.get()->operator [](i).SetAlpha(alpha);
            memory_vector.get()->operator [](i).SetSe(Se);
        }
        
    }
}


void TPMRSSegregatedAnalysis::AdjustIntegrationOrder(TPZCompMesh * cmesh_o, TPZCompMesh * cmesh_d){
    
    // Assuming the cmesh_o as directive.
    
    cmesh_d->LoadReferences();
    int nel_o = cmesh_o->NElements();
    int nel_d = cmesh_d->NElements();
    
    if (nel_o != nel_d) {
        std::cout << "The geometrical partitions are not the same." << std::endl;
        DebugStop();
    }
    
    for (long el = 0; el < nel_o; el++) {
        TPZCompEl *cel_o = cmesh_o->Element(el);
        if (!cel_o) {
            continue;
        }
        
        TPZGeoEl * gel = cel_o->Reference();
        if (!gel) {
            continue;
        }
        
        // Finding the other computational element
        TPZCompEl * cel_d = gel->Reference();
        if (!cel_d) {
            continue;
        }
        cel_o->SetFreeIntPtIndices();
        cel_o->ForcePrepareIntPtIndices();
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

void TPMRSSegregatedAnalysis::ExecuteOneTimeStep(){
    m_reservoir_analysis->ExecuteOneTimeStep();
    m_geomechanic_analysis->ExecuteOneTimeStep();
}

void TPMRSSegregatedAnalysis::PostProcessTimeStep(std::string & geo_file, std::string & res_file){
    m_reservoir_analysis->PostProcessTimeStep(res_file);
    m_geomechanic_analysis->PostProcessTimeStep(geo_file);
}

void TPMRSSegregatedAnalysis::ExecuteTimeEvolution(){
    
    // vtk files
    std::string name = m_simulation_data->name_vtk_file();
    std::string file_geo = name + "_geo.vtk";
    std::string file_res = name + "_res.vtk";
    
    int n_max_fss_iterations = 20; // @TODO:: MS, please to xml file structure
    int n_enforced_fss_iterations = 5; // @TODO:: MS, please to xml file structure
    int n_time_steps = m_simulation_data->ReportingTimes().size();
    REAL r_norm = m_simulation_data->epsilon_res();
    REAL dx_norm = m_simulation_data->epsilon_cor();

    bool error_stop_criterion_Q = false;
    bool dx_stop_criterion_Q = false;
    for (int it = 0; it < n_time_steps; it++) {
        for (int k = 1; k <= n_max_fss_iterations; k++) {
            this->ExecuteOneTimeStep();
            error_stop_criterion_Q = (m_reservoir_analysis->Get_error() < r_norm) && (m_geomechanic_analysis->Get_error() < r_norm);
            dx_stop_criterion_Q = (m_reservoir_analysis->Get_dx_norm() < dx_norm) && (m_geomechanic_analysis->Get_dx_norm() < dx_norm);
//            this->PostProcessTimeStep(file_geo, file_res);
            if ((error_stop_criterion_Q && (k > n_enforced_fss_iterations)) || dx_stop_criterion_Q) {
                std::cout << "TPMRSSegregatedAnalysis:: Iterative process converged with residue norm for res = " << m_reservoir_analysis->Get_error() << std::endl;
                std::cout << "TPMRSSegregatedAnalysis:: Iterative process converged with residue norm for geo = " << m_geomechanic_analysis->Get_error() << std::endl;
//                m_geomechanic_analysis->AssembleResidual();
                UpdateState();
                this->PostProcessTimeStep(file_geo, file_res);
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
    if (!cmesh) {
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
        
        // Update the elastic response
        {
            std::tuple<TPMRSUndrainedParameters, TPMRSPoroMechParameters, TPMRSPhiParameters,TPMRSKappaParameters,TPMRSPlasticityParameters> chunk =    m_simulation_data->MaterialProps()[iregion];
            
            REAL E,nu;
            
            // Elastic predictor
            if (IsInitialConditionsQ) {
                TPMRSUndrainedParameters undrained_parameters(std::get<0>(chunk));
                std::vector<REAL> u_pars = undrained_parameters.GetParameters();
                E = u_pars[0];
                nu = u_pars[1];
            }else{
                TPMRSPoroMechParameters poroperm_parameters(std::get<1>(chunk));
                std::vector<REAL> poroperm_pars = poroperm_parameters.GetParameters();
                E = poroperm_pars[0];
                nu = poroperm_pars[1];
            }
            
            // Updating bulk modulus for porosity model
            std::get<2>(m_simulation_data->MaterialProps()[iregion]).SetBulkModulus(E, nu);
            
            TPZElasticResponse ER;
            ER.SetUp(E, nu);
            
            // Plastic corrector
            TPMRSPlasticityParameters plasticity_parameters(std::get<4>(chunk));
            std::vector<REAL> p_pars = plasticity_parameters.GetParameters();
            
            if (p_pars.size() == 0) {
                // Elastic material
                TPZElasticCriterion Elastic;
                Elastic.SetElasticResponse(ER);
                
                TPMRSElastoPlastic<TPZElasticCriterion, TPMRSMemory> * vol_material = dynamic_cast< TPMRSElastoPlastic<TPZElasticCriterion, TPMRSMemory> * >(material);
                vol_material->SetPlasticIntegrator(Elastic);
                
            }else{
                // Elastoplastic material
                switch (plasticity_parameters.GetModel()) {
                    case plasticity_parameters.ep_mc: {
                        
                        // Mohr Coulomb data
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
                    default:{
                        std::cout << "Material not implemented. " << std::endl;
                        DebugStop();
                    }
                        break;
                }
            }
        }
        
        // Inserting boundary conditions
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
        
        std::shared_ptr<TPZAdmChunkVector<TPMRSMemory>> & memory_vector = mat_with_memory->GetMemory();
        
        int ndata = memory_vector->NElements();
        for (int i = 0; i < ndata; i++) {
            // Because we reused the same memory items
            TPZTensor<REAL> sigma_total_0 = memory_vector.get()->operator [](i).GetSigma_n();
            REAL p_0 = -(sigma_total_0.I1()/3);
            memory_vector.get()->operator [](i).Setp_0(p_0);
            memory_vector.get()->operator [](i).Setp(p_0);
            memory_vector.get()->operator [](i).Setp_n(p_0);
            sigma_total_0.Zero();// Converted to effecttive because initial deformation is Zero.
            memory_vector.get()->operator [](i).SetSigma_0(sigma_total_0);
            memory_vector.get()->operator [](i).SetSigma(sigma_total_0);
            memory_vector.get()->operator [](i).SetSigma_n(sigma_total_0);
            // Cleaning u
            memory_vector.get()->operator [](i).Setu_0(u_null);
            memory_vector.get()->operator [](i).Setu(u_null);
            memory_vector.get()->operator [](i).Setu_n(u_null);
            // Cleaning elasto-plastic states
            memory_vector.get()->operator [](i).GetPlasticState_0().CleanUp();
            memory_vector.get()->operator [](i).GetPlasticState().CleanUp();
            memory_vector.get()->operator [](i).GetPlasticState_n().CleanUp();

        }
        
    }
    
}
