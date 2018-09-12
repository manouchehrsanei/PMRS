//
//  TPMRSGeomechanicAnalysis.h
//  PMRS
//
//  Created by Omar Durán on 9/11/18.
//

#ifndef TPMRSGeomechanicAnalysis_h
#define TPMRSGeomechanicAnalysis_h

#include <stdio.h>
#include "pzanalysis.h"
#include "TPZSimulationData.h"
#include "pzpostprocanalysis.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "TPZSSpStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzfstrmatrix.h"
#include "pzstepsolver.h"
#include "pzbuildmultiphysicsmesh.h"

class TPMRSGeomechanicAnalysis : public TPZAnalysis {
    
private:
    
    /// Pointer of Simulation data object
    TPZSimulationData * m_simulation_data;
    
    /// Solution at n+1 state
    TPZFMatrix<STATE> m_X_n;
    
    /// Solution at n (past) state
    TPZFMatrix<STATE> m_X;
    
    /// Residue error
    STATE m_error;
    
    /// Correction variation
    STATE m_dx_norm;
    
    /// number of Newton iterations
    int m_k_iterations;
    
    /// Post-processor object
    TPZPostProcAnalysis * m_post_processor;
    
    /// Variables being postprocessed
    TPZStack<std::string> m_var_names;
    
public:
    
    /// Default constructor
    TPMRSGeomechanicAnalysis();
    
    /// Destructor
    ~TPMRSGeomechanicAnalysis();
    
    /// Copy constructor
    TPMRSGeomechanicAnalysis(const TPMRSGeomechanicAnalysis & other);
    
    /// Set the pointer of Simulation data object
    void SetSimulationData(TPZSimulationData * simulation_data){
        m_simulation_data = simulation_data;
    }
    
    /// Configurate the solver being used to compute the approximation
    void ConfigurateAnalysis(DecomposeType decomposition, TPZSimulationData * simulation_data);
    
    /// Execute a single newton iteration
    void ExecuteNewtonInteration();
    
    /// Execute the evolution for a single pseudo time step
    void ExecuteOneTimeStep();
    
    /// Post-processing the variables for a single pseudo time step
    void PostProcessTimeStep(std::string & file);
    
    /// Update the memory with the converged pseudo time step solution
    void AcceptPseudoTimeStepSolution();
    
    /// Load the current state for the hdiv and 2 meshes
    void LoadCurrentState();
    
    /// Load the last state for the hdiv and 2 meshes
    void LoadLastState();
    
};

#endif /* TPMRSGeomechanicAnalysis_h */
