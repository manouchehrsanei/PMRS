//
//  TPMRSSegregatedAnalysis.h
//  PMRS
//
//  Created by Omar Durán on 9/13/18.
//

#ifndef TPMRSSegregatedAnalysis_h
#define TPMRSSegregatedAnalysis_h

#include <stdio.h>
#include "pzbndcond.h"
#include "TPZMatWithMem.h"
#include "TPMRSMemory.h"
#include "TPMRSGeomechanicAnalysis.h"
#include "TPMRSMonoPhasicAnalysis.h"
#include "TPZElasticResponse.h"
#include "TPZElasticCriterion.h"
#include "TPZPlasticStepPV.h"
#include "TPZSandlerExtended.h"
#include "TPZYCMohrCoulombPV.h"
#include "TPMRSElastoPlastic.h"

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
    void ApplyMemoryLink(TPZCompMesh * cmesh_o, TPZCompMesh * cmesh_d);
    
    /// Adjust integration orders
    void AdjustIntegrationOrder(TPZCompMesh * cmesh_o, TPZCompMesh * cmesh_d);
    
    /// Configurate internal analysis objects and linked them through the memory shared pointer
    void ConfigurateAnalysis(DecomposeType decompose_geo, DecomposeType decompose_res, TPZSimulationData * simulation_data, TPZCompMesh * cmesh_geomechanics, TPZCompMesh * cmesh_reservoir, TPZManVector<TPZCompMesh * , 2> & mesh_vec);
    
    /// Fill Memory with alpha, phi_0, k_0 and Se
    void FillMemory(TPZCompMesh * cmesh);
    
    /// Execute the evolution for a single time step
    void ExecuteOneTimeStep();
    
    /// Post-processing the variables for a single time step
    void PostProcessTimeStep(std::string & geo_file, std::string & res_file);
    
    /// Execute the transient evolution using Fixed Stress Split Iteration
    void ExecuteTimeEvolution();
    
    /// Update solution state x = x_n
    void UpdateState();
    
    /// Update the initial solution for sigma and pressure
    void UpdateInitialSigmaAndPressure();
    
    /// Configurate boundary conditions (IsInitialConditionsQ is false set recurrent BC's)
    void ConfigurateBConditions(bool IsInitialConditionsQ);
    
    /// Execute initial problem
    void ExecuteStaticSolution();
    

    
};

#endif /* TPMRSSegregatedAnalysis_h */
