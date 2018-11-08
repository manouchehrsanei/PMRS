//
//  TPMRSMonoPhasicAnalysis.h
//  PMRS
//
//  Created by Omar and Manouchehr on 9/11/18.
//

#ifndef TPMRSMonoPhasicAnalysis_h
#define TPMRSMonoPhasicAnalysis_h

#include <stdio.h>
#include "pzanalysis.h"
#include "TPMRSSimulationData.h"
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
    TPMRSSimulationData * m_simulation_data;
    
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
    
    /// Variables being postprocessed
    TPZStack<std::string> m_vec_var_names;
    
public:
    
    /// Default constructor
    TPMRSMonoPhasicAnalysis();
    
    /// Destructor
    ~TPMRSMonoPhasicAnalysis();
    
    /// Copy constructor
    TPMRSMonoPhasicAnalysis(const TPMRSMonoPhasicAnalysis & other);
    
    /// Set the pointer of Simulation data object
    void SetSimulationData(TPMRSSimulationData * simulation_data){
        m_simulation_data = simulation_data;
    }
    
    /// Configurate the solver being used to compute the approximation
    void ConfigurateAnalysis(DecomposeType decomposition, TPZManVector<TPZCompMesh * , 2> & mesh_vec,TPMRSSimulationData * simulation_data);
    
    /// Execute a single newton iteration
    void ExecuteNewtonInteration();
    
    /// Execute the evolution for a single time step
    void ExecuteOneTimeStep();
    
    /// Project the undrained initial pressure over the pressure space 
    void ExecuteUndrainedResponseStep();
    
    /// Post-processing the variables for a single time step
    void PostProcessTimeStep(std::string & file);
    
    /// Update the memory with the converged time step solution
    void LoadMemorySolution();
    
    /// Load the current state for the hdiv and 2 meshes
    void LoadCurrentState();
    
    /// Load the last state for the hdiv and 2 meshes
    void LoadLastState();
    
    /// Update solution state x = x_n
    void UpdateState();
    
    /// Brief Set Residue error
    void Set_error(STATE error)
    {
        m_error = error;
    }
    
    /// Brief Get Residue error
    STATE Get_error()
    {
        return m_error;
    }
    
    
    /// Brief Set Correction variation
    void Set_dx_norm(STATE dxnorm)
    {
        m_dx_norm = dxnorm;
    }
    
    /// Brief Get Correction variation
    STATE Get_dx_norm()
    {
        return m_dx_norm;
    }
    
    
    /// Brief Set number of Newton iterations
    void Set_k_iterations(int kiterations)
    {
        m_k_iterations = kiterations;
    }
    
    /// Brief Get number of Newton iterations
    int Get_k_iterations()
    {
        return m_k_iterations;
    }
    
    
};

#endif /* TPMRSMonoPhasicAnalysis_h */
