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
    TPZManVector<std::pair<int, TPZManVector<int,12>>,12>  material_ids = m_simulation_data->MaterialIds();
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
    bool mustOptimizeBandwidth = true;

    if (simulation_data->ElasticityOrder()>simulation_data->DiffusionOrder()) {
        this->ApplyMemoryLink(cmesh_geomechanics,cmesh_reservoir);
        this->AdjustIntegrationOrder(cmesh_geomechanics,cmesh_reservoir);
    }else{
        this->ApplyMemoryLink(cmesh_reservoir,cmesh_geomechanics);
        this->AdjustIntegrationOrder(cmesh_reservoir,cmesh_geomechanics);
    }

    
    // The Geomechanics Simulator
    m_geomechanic_analysis = new TPMRSGeomechanicAnalysis;
    m_geomechanic_analysis->SetCompMesh(cmesh_geomechanics,mustOptimizeBandwidth);
    m_geomechanic_analysis->ConfigurateAnalysis(decompose_geo, m_simulation_data);
    
    // The Reservoir Simulator
    m_reservoir_analysis = new TPMRSMonoPhasicAnalysis;
    m_reservoir_analysis->SetCompMesh(cmesh_reservoir,mustOptimizeBandwidth);
    m_reservoir_analysis->ConfigurateAnalysis(decompose_res, mesh_vec, m_simulation_data);

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
    
    std::string file_res("ReservoirFlow.vtk");
    std::string file_geo("Geomechanic.vtk");
    
    int n_max_fss_iterations = 10; // @TODO:: MS, please to xml file structure
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
            
            if ((error_stop_criterion_Q && (k > n_enforced_fss_iterations)) || dx_stop_criterion_Q) {
                std::cout << "TPMRSSegregatedAnalysis:: Iterative process converged with residue norm for res = " << m_reservoir_analysis->Get_error() << std::endl;
                std::cout << "TPMRSSegregatedAnalysis:: Iterative process converged with residue norm for geo = " << m_geomechanic_analysis->Get_error() << std::endl;
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
