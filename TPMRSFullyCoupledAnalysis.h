//
//  TPMRSFullyCoupledAnalysis.h
//  PZ
//
//  Created by Omar and Manouchehr on 9/11/18.
//
//

#ifndef TPMRSFullyCoupledAnalysis_h
#define TPMRSFullyCoupledAnalysis_h

#include <stdio.h>
#include "pzbndcond.h"
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

#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif

class TPMRSFullyCoupledAnalysis : public TPZAnalysis
{
    
private:
    
    /// define the simulation data
    TPMRSSimulationData * m_simulation_data;
    
    /// Vector of compmesh pointers. fmeshvec[0] = flowHdiv, fmeshvec[1] = PressureL2
    TPZManVector<TPZCompMesh * , 2> m_meshvec;
    
    /// Solution at n+1 state
    TPZFMatrix<STATE> m_X_n;
    
    /// Solution at n (past) state
    TPZFMatrix<STATE> m_X;
    
    /// Residue error
    STATE m_error;
    
    /// Correction variation
    STATE m_dx_norm;
    
    /// number of newton corrections
    int m_k_iterations;
    
    /// Post-processor object
    TPZPostProcAnalysis * m_post_processor;
    
    /// Variables being postprocessed
    TPZStack<std::string> m_var_names;
    
    
public:
    
    /// default Constructor
    TPMRSFullyCoupledAnalysis();
    
    /// default Destructor
    ~TPMRSFullyCoupledAnalysis();
    
    /// Copy constructor
    TPMRSFullyCoupledAnalysis(const TPMRSFullyCoupledAnalysis &other);
    
    /// Copy assignemnt operator
    TPMRSFullyCoupledAnalysis &operator=(const TPMRSFullyCoupledAnalysis &other);
    
    
    /// Set Solution at n+1 (current) state
    void SetX_n(TPZFMatrix<STATE> &x)
    {
        m_X_n = x;
    }
    
    /// Get Solution at n+1 (current) state
    TPZFMatrix<STATE> & X_n()
    {
        return m_X_n;
    }
    
    /// Set Solution at n (last) state
    void SetX(TPZFMatrix<STATE> &x)
    {
        m_X = x;
    }
    
    /// Get Solution at n (last) state
    TPZFMatrix<STATE> & X()
    {
        return m_X;
    }
    
    
    /// Set the simulation data
    void SetSimulationData(TPMRSSimulationData * SimulationData)
    {
        m_simulation_data = SimulationData;
        m_meshvec.Resize(2);
    }
    
    /// Get the space generator
    TPMRSSimulationData * SimulationData()
    {
        return m_simulation_data;
    }

    /// Set vector of compmesh pointers. fmeshvec[0] = flux, fmeshvec[1] = Pressure
    void SetMeshvec(TPZManVector<TPZCompMesh * , 2> &Meshvec)
    {
        m_meshvec = Meshvec;
    }
    
    /// Get Vector of compmesh pointers. fmeshvec[0] = flux, fmeshvec[1] = Pressure
    TPZManVector<TPZCompMesh * , 2> & Meshvec()
    {
        return m_meshvec;
    }
    
    /// Resize and fill residue and solution vectors
    void AdjustVectors();
    
    /// Set k iterations
    void Set_k_ietrarions(int k)
    {
        m_k_iterations = k;
    }
    
    /// Get k iterations
    int k_ietrarions()
    {
        return m_k_iterations;
    }
    
    /// Configurate postprocessor
    void ConfiguratePostProcessor();
    
    /// Execute an time step with supstepping control
    void ExecuteOneTimeStep(int i_time_step);
    
    /// Execute an step approximation
    bool ExcecuteOneStepApproximation(bool enforced_execution_Q);
    
    /// Execute a single newton iteration
    void ExecuteNewtonInteration();
    
    /// Execute a single ninth order newton iteration
    void ExecuteNinthOrderNewtonInteration(REAL & norm_dx);
    
    /// Function to decide the postprocess time directive
    bool ShouldPostprocessQ(REAL time);
    
    /// Post-processing the variables for a single time step
    void PostProcessTimeStep(std::string file);
    
    /// Update solution state x = x_n
    void UpdateState();
    
    /// Configurate boundary conditions
    void ConfigureFCBC(REAL t);
    
    /// Update the memory with the converged pseudo time step solution
    void LoadMemorySolution(TPZFMatrix<REAL> & x);
    
    /// Update state x on cmesh solution
    void LoadState(TPZFMatrix<REAL> & x);

    /// execute the evolutionary problem
    void ExecuteTimeEvolution();
    
    /// Auxiliary function for compute power of integers
    int power(int base, int exp)
    {
        int result = 1;
        while (exp != 0)
        {
            if ((exp & 1) == 1)
                result *= base;
            exp >>= 1;
            base *= base;
        }
        
        return result;
    }
    
};


#endif /* TPMRSFullyCoupledAnalysis_h */
