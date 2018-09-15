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

void TPMRSSegregatedAnalysis::ApplyMemoryLink(){
    
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
        
        TPZMaterial * material_geo = m_geomechanic_analysis->Mesh()->FindMaterial(matid);
        TPZMaterial * material_res = m_reservoir_analysis->Mesh()->FindMaterial(matid);
        if (!material_geo || !material_res) {
            DebugStop();
        }
        
        TPZMatWithMem<TPMRSMemory> * mat_with_memory_geo = dynamic_cast<TPZMatWithMem<TPMRSMemory> * >(material_geo);
        TPZMatWithMem<TPMRSMemory> * mat_with_memory_res = dynamic_cast<TPZMatWithMem<TPMRSMemory> * >(material_res);
        if (!mat_with_memory_geo || !mat_with_memory_res) {
            DebugStop();
        }
        mat_with_memory_res->SetMemory(mat_with_memory_geo->GetMemory());
    }
}

void TPMRSSegregatedAnalysis::ConfigurateAnalysis(DecomposeType decompose_geo, DecomposeType decompose_res, TPZSimulationData * simulation_data, TPZCompMesh * cmesh_geomechanics, TPZCompMesh * cmesh_reservoir, TPZManVector<TPZCompMesh * , 2> & mesh_vec){

    if (!simulation_data || !cmesh_geomechanics || !cmesh_reservoir) {
        DebugStop();
    }
    
    this->SetSimulationData(simulation_data);
    bool mustOptimizeBandwidth = true;

    // The Geomechanics Simulator
    m_geomechanic_analysis = new TPMRSGeomechanicAnalysis;
    m_geomechanic_analysis->SetCompMesh(cmesh_geomechanics,mustOptimizeBandwidth);
    m_geomechanic_analysis->ConfigurateAnalysis(decompose_geo, m_simulation_data);
    
    // The Reservoir Simulator
    m_reservoir_analysis = new TPMRSMonoPhasicAnalysis;
    m_reservoir_analysis->SetCompMesh(cmesh_reservoir,mustOptimizeBandwidth);
    m_reservoir_analysis->ConfigurateAnalysis(decompose_res, mesh_vec, m_simulation_data);
    
    this->ApplyMemoryLink();
    
}

void TPMRSSegregatedAnalysis::ExecuteOneTimeStep(bool must_accept_solution_Q){
    m_reservoir_analysis->ExecuteOneTimeStep(must_accept_solution_Q);
    m_geomechanic_analysis->ExecuteOneTimeStep(must_accept_solution_Q);
}

void TPMRSSegregatedAnalysis::PostProcessTimeStep(std::string & geo_file, std::string & res_file){
    m_reservoir_analysis->PostProcessTimeStep(res_file);
    m_geomechanic_analysis->PostProcessTimeStep(geo_file);
}

void TPMRSSegregatedAnalysis::ExecuteTimeEvolution(){
    
    std::string file_res("ReservoirFlow.vtk");
    std::string file_geo("Geomechanic.vtk");
    
    int n_max_fss_iterations = 20; // @TODO:: MS, please to xml file structure
    int n_enforced_fss_iterations = 10; // @TODO:: MS, please to xml file structure
    int n_time_steps = m_simulation_data->ReportingTimes().size();
    REAL r_norm = m_simulation_data->epsilon_res();
    REAL dx_norm = m_simulation_data->epsilon_cor();

    bool error_stop_criterion_Q = false;
    bool dx_stop_criterion_Q = false;
    for (int it = 0; it < n_time_steps; it++) {
        for (int k = 1; k <= n_max_fss_iterations; k++) {
            this->ExecuteOneTimeStep(false);
            error_stop_criterion_Q = (m_reservoir_analysis->Get_error() < r_norm) && (m_geomechanic_analysis->Get_error() < r_norm);
            dx_stop_criterion_Q = (m_reservoir_analysis->Get_dx_norm() < dx_norm) && (m_geomechanic_analysis->Get_dx_norm() < dx_norm);
            
            if ((error_stop_criterion_Q && (k > n_enforced_fss_iterations)) || dx_stop_criterion_Q) {
                this->ExecuteOneTimeStep(true);
                this->PostProcessTimeStep(file_geo, file_res);
                std::cout << "TPMRSSegregatedAnalysis:: Iterative process converged with residue norm for res = " << m_reservoir_analysis->Get_error() << std::endl;
                std::cout << "TPMRSSegregatedAnalysis:: Iterative process converged with residue norm for geo = " << m_geomechanic_analysis->Get_error() << std::endl;
                break;
            }
        }
    }
    
}
