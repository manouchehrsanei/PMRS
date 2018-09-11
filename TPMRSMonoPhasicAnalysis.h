//
//  TPMRSMonoPhasicAnalysis.h
//  PMRS
//
//  Created by Omar Durán on 9/9/18.
//

#ifndef TPMRSMonoPhasicAnalysis_h
#define TPMRSMonoPhasicAnalysis_h

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

class TPMRSMonoPhasicAnalysis : public TPZAnalysis {
    
private:
    
    /// Pointer of Simulation data object
    TPZSimulationData * m_simulation_data;
    
    /// Solution at n+1 state
    TPZFMatrix<STATE> m_X_n;
    
    /// Solution at n (past) state
    TPZFMatrix<STATE> m_X;
    
    /// Vector of compmesh pointers. fmeshvec[0] = flowHdiv, fmeshvec[1] = PressureL2 */
    TPZManVector<TPZCompMesh * , 2> m_mesh_vec;
    
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
    TPMRSMonoPhasicAnalysis();
    
    /// Destructor
    ~TPMRSMonoPhasicAnalysis();
    
    /// Copy constructor
    TPMRSMonoPhasicAnalysis(const TPMRSMonoPhasicAnalysis & other);
    
    /// Set the pointer of Simulation data object
    void SetSimulationData(TPZSimulationData * simulation_data){
        m_simulation_data = simulation_data;
    }
    
    /// Configurate the solver being used to compute the approximation
    void ConfigurateAnalysis(DecomposeType decomposition, TPZManVector<TPZCompMesh * , 2> & mesh_vec,TPZSimulationData * simulation_data);
    
    /// Execute a single newton iteration
    void ExecuteNewtonInteration();
    
    /// Execute the evolution for a single time step
    void ExecuteOneTimeStep();
    
    /// Post-processing the variables for a single time step
    void PostProcessTimeStep(std::string & file);
    
    /// Update the memory with the converged time step solution
    void AcceptTimeStepSolution();
    
    /// Load the current state for the hdiv and 2 meshes
    void LoadCurrentState();
    
    /// Load the last state for the hdiv and 2 meshes
    void LoadLastState();
    
};

#endif /* TPMRSMonoPhasicAnalysis_h */
