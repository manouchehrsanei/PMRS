//
//  TPMRSGeomechanicAnalysis.h
//  PMRS
//
//  Created by Omar and Manouchehr on 9/11/18.
//

#ifndef TPMRSGeomechanicAnalysis_h
#define TPMRSGeomechanicAnalysis_h

#include <stdio.h>
#include "pzanalysis.h"
#include "TPMRSSimulationData.h"
#include "pzpostprocanalysis.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "TPZSSpStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "pzfstrmatrix.h"
#include "pzstepsolver.h"
#include "pzbuildmultiphysicsmesh.h"

class TPMRSGeomechanicAnalysis : public TPZAnalysis {
    
private:
    
    /// Pointer of Simulation data object
    TPMRSSimulationData * m_simulation_data;
    
    /// Solution at n+1 state
    TPZFMatrix<STATE> m_X_n;
    
    /// Solution at n (past) state
    TPZFMatrix<STATE> m_X;

    /// Solution at n (past) state
    TPZFMatrix<STATE> m_X_tilde;
    
    /// Modified Thomas cceleration factor
    STATE m_alpha = 1.0;
    
    /// Residue error
    STATE m_error;
    
    /// Correction variation
    STATE m_dx_norm;
    
    /// number of Newton iterations
    int m_k_iterations;
    
    /// number of Quase-Newton iterations for the next Jacobian update
    int m_n_update_jac;
    
    /// Post-processor object
    TPZPostProcAnalysis * m_post_processor;
    
    /// Variables being postprocessed
    TPZStack<std::string> m_var_names;
    
    /// Geomechanics states for nonlinear acceleration
    TPZManVector<TPZFMatrix<REAL>,10> m_x_u;
    
public:
    
    /// Default constructor
    TPMRSGeomechanicAnalysis();
    
    /// Destructor
    ~TPMRSGeomechanicAnalysis();
    
    /// Copy constructor
    TPMRSGeomechanicAnalysis(const TPMRSGeomechanicAnalysis & other);
    
    /// Set the pointer of Simulation data object
    void SetSimulationData(TPMRSSimulationData * simulation_data){
        m_simulation_data = simulation_data;
    }
    
    /// Configurate the solver being used to compute the approximation
    void ConfigurateAnalysis(DecomposeType decomposition, TPMRSSimulationData * simulation_data);
    
    /// Execute a single iteration
    void ExecuteInteration(REAL & norm_dx);
    
    /// Execute a single M1 iteration
    void ExecuteM1Interation(REAL & norm_dx);
    
    /// Execute a single Modified Thomas MT iteration
    void ExecuteMTInteration(REAL & norm_dx);
    
    /// Execute a single M3 iteration
    void ExecuteM3Interation(REAL & norm_dx);
    
    /// Execute a single M6 iteration
    void ExecuteM6Interation(REAL & norm_dx);
    
    /// Execute a single ninth order newton iteration
    void ExecuteNinthOrderNewtonInteration(REAL & norm_dx);
    
    /// Execute the evolution for a single pseudo time step
    bool ExecuteOneTimeStep(bool enforced_execution_Q);
    
    /// Execute the undrained respose and update the initial sigma and pressure
    void ExecuteUndrainedResponseStep();
    
    /// Post-processing the variables for a single pseudo time step
    void PostProcessTimeStep(std::string & file);
    
    /// Update the memory with the converged pseudo time step solution
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
    
    /// Set number of Quase-Newton iterations for jacobian update
    void Set_n_update_jac(int n_update_jac)
    {
        m_n_update_jac = n_update_jac;
    }
    
    /// Get number of Quase-Newton iterations for jacobian update
    int Get_n_update_jac()
    {
        return m_n_update_jac;
    }
    
    /// Set Solution at n+1 state
    void SetX_n(TPZFMatrix<STATE> & X_n)
    {
        m_X_n = X_n;
    }
    
    /// Get Solution at n+1 state
    TPZFMatrix<STATE> & X_n()
    {
        return m_X_n;
    }
    
    // Applying the selected nonlinear acceleration
    void ApplyAcceleration();
    
    /// Define an acceleration method for the outer loop k iteration for geomechanics module
    void AccelerationGeo(int k, int n);
    
    /// Apply a tranformation formula based on three states
    TPZFMatrix<REAL> ApplyTransformation(TPZFMatrix<REAL> & An_p_1, TPZFMatrix<REAL> & An, TPZFMatrix<REAL> & An_m_1);
    
    /// Apply a Atiken Delta-2 tranformation formula based on three states
    TPZFMatrix<REAL> FDMTransformation(TPZFMatrix<REAL> & An_p_1, TPZFMatrix<REAL> & An, TPZFMatrix<REAL> & An_m_1);
    
    /// Apply a Anderson tranformation formula based on three states
    TPZFMatrix<REAL> SDMTransformation(TPZFMatrix<REAL> & An_p_1, TPZFMatrix<REAL> & An, TPZFMatrix<REAL> & An_m_1);
    
};

#endif /* TPMRSGeomechanicAnalysis_h */
