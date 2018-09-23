//
//  TPMRSSegregatedAnalysis.h
//  PMRS
//
//  Created by Omar Durán on 9/13/18.
//

#ifndef TPMRSSegregatedAnalysis_h
#define TPMRSSegregatedAnalysis_h

#include <stdio.h>
#include "TPZMatWithMem.h"
#include "TPMRSMemory.h"
#include "TPMRSGeomechanicAnalysis.h"
#include "TPMRSMonoPhasicAnalysis.h"

class TPMRSSegregatedAnalysis {
    
private:
    
    /// Pointer to Simulation data object
    TPZSimulationData * m_simulation_data;
    
    /// Pointer to geomechanic analysis object
    TPMRSGeomechanicAnalysis * m_geomechanic_analysis;
    
    /// Pointer to reservoir analysis object
    TPMRSMonoPhasicAnalysis * m_reservoir_analysis;
    
public:
    
    /// Default constructor
    TPMRSSegregatedAnalysis();
    
    /// Destructor
    ~TPMRSSegregatedAnalysis();
    
    /// Copy constructor
    TPMRSSegregatedAnalysis(const TPMRSSegregatedAnalysis & other);
    
    /// Set the pointer of Simulation data object
    void SetSimulationData(TPZSimulationData * simulation_data){
        m_simulation_data = simulation_data;
    }

    /// Attach materials with the same share pointer (For now just volumeric linking)
    void ApplyMemoryLink();
    
    /// Configurate internal analysis objects and linked them through the memory shared pointer
    void ConfigurateAnalysis(DecomposeType decompose_geo, DecomposeType decompose_res, TPZSimulationData * simulation_data, TPZCompMesh * cmesh_geomechanics, TPZCompMesh * cmesh_reservoir, TPZManVector<TPZCompMesh * , 2> & mesh_vec);
    
    /// Execute the evolution for a single time step
    void ExecuteOneTimeStep(bool must_accept_solution_Q);
    
    /// Post-processing the variables for a single time step
    void PostProcessTimeStep(std::string & geo_file, std::string & res_file);
    
    /// Execute the transient evolution using Fixed Stress Split Iteration
    void ExecuteTimeEvolution();
    
    /// Update solution state x = x_n
    void UpdateState();
    
};

#endif /* TPMRSSegregatedAnalysis_h */
