//
//  TPMRSSegregatedAnalysis.h
//  PMRS
//
//  Created by Omar and Manouchehr on 9/13/18.
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
#include "pzfmatrix.h"

#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif

class TPMRSSegregatedAnalysis {
    
private:
    
    /// Pointer to Simulation data object
    TPMRSSimulationData * m_simulation_data;
    
    /// Pointer to geomechanic analysis object
    TPMRSGeomechanicAnalysis * m_geomechanic_analysis;
    
    /// Pointer to reservoir analysis object
    TPMRSMonoPhasicAnalysis * m_reservoir_analysis;
    
    /// Object that store the history of segregated and internal iteraions
    TPZFMatrix<REAL> m_iterations_summary;
    
    /// Object that store the cpu time of segregated and internal process
    TPZFMatrix<REAL> m_cpu_time_summary;
    
    /// Object that store the residuals of segregated and internal process
    TPZFMatrix<REAL> m_residuals_summary;
    
    TPZFMatrix<REAL> m_xp_m;
    
    TPZFMatrix<REAL> m_xp_m_1;
    
    TPZFMatrix<REAL> m_xp_m_2;
    
    TPZFMatrix<REAL> m_xu_m;
    
    TPZFMatrix<REAL> m_xu_m_1;
    
    TPZFMatrix<REAL> m_xu_m_2;
    
    /// Reservoir state at last fixed stress cycle
    TPZFMatrix<REAL> m_p_m;
    
    /// Geomechanic state at last fixed stress cycle
    TPZFMatrix<REAL> m_u_m;
    
public:
    
    /// Default constructor
    TPMRSSegregatedAnalysis();
    
    /// Destructor
    ~TPMRSSegregatedAnalysis();
    
    /// Copy constructor
    TPMRSSegregatedAnalysis(const TPMRSSegregatedAnalysis & other);
    
    /// Set the pointer of Simulation data object
    void SetSimulationData(TPMRSSimulationData * simulation_data);

    /// Attach materials with the same share pointer (For now just volumeric linking)
    void ApplyMemoryLink(TPZCompMesh * cmesh_o, TPZCompMesh * cmesh_d);
    
    /// Adjust integration orders
    void AdjustIntegrationOrder(TPZCompMesh * cmesh_o, TPZCompMesh * cmesh_d);
    
    /// Configurate internal analysis objects and linked them through the memory shared pointer
    void ConfigurateAnalysis(DecomposeType decompose_geo, DecomposeType decompose_res, TPMRSSimulationData * simulation_data, TPZCompMesh * cmesh_geomechanics, TPZCompMesh * cmesh_reservoir, TPZManVector<TPZCompMesh * , 2> & mesh_vec);
    
    /// Fill Memory with alpha, phi_0, k_0 and Se
    void FillMemory(TPZCompMesh * cmesh);
    
    /// Execute the evolution for a single time step
    void ExecuteOneTimeStep(int i_time_step, int k);
    
    /// Post-processing the variables for a single time step
    void PostProcessTimeStep(std::string & geo_file, std::string & res_file);
    
    /// Execute the transient evolution using Fixed Stress Split Iteration
    void ExecuteTimeEvolution();
    
    /// Update solution state x = x_n
    void UpdateState();
    
    /// Update the initial solution for sigma and pressure
    void UpdateInitialSigmaAndPressure();
    
    /// Configurate boundary conditions (IsInitialConditionsQ is false set recurrent BC's)
    void ConfigureGeomechanicsBC(REAL t, bool IsInitialConditionsQ = false);
    
    /// Configurate boundary conditions (IsInitialConditionsQ is false set recurrent BC's)
    void ConfigureReservoirBC(REAL t, bool IsInitialConditionsQ = false);
    
    /// Execute initial problem
    void ExecuteStaticSolution();
    
    /// Get the object that store the history of segregated and internal iteraions
    TPZFMatrix<REAL> & IterationsSummary();
    
    /// Get the object that store the cpu time of segregated and internal process
    TPZFMatrix<REAL> & TimeSummary();

    /// Get the object that store the residuals of segregated and internal process
    TPZFMatrix<REAL> & ResidualsSummary();
    
    /// Resize and storage the positions to track iteraions and cpu time during the entire simulation
    void ConfigurateHistorySummaries();
    
    REAL linear_extrapolation(REAL & f_1, REAL & f_2, REAL & t_1, REAL & t_2, REAL & t);
    
    REAL quadratic_extrapolation(REAL & f_1, REAL & f_2, REAL & f_3, REAL & t_1, REAL & t_2,  REAL & t_3, REAL & t);
    
    /// Define a quase Newton acceleration for the outer loop k iteration for reservoir module
    void QNAccelerationRes(int k);
    
    /// Define a quase Newton acceleration for the outer loop k iteration for geomechanic module
    void QNAccelerationGeo(int k);
    
    /// Define a  Aitken acceleration for the outer loop k iteration for reservoir module
    void AitkenAccelerationRes(int k);
    
    /// Define a  Aitken acceleration for the outer loop k iteration for geomechanic module
    void AitkenAccelerationGeo(int k);
    
    /// Define a  Gauss-Seidel acceleration for the outer loop k iteration for reservoir module
    void GaussSeidelAccelerationRes(int k);
    
    /// Define a  Gauss-Seidel acceleration for the outer loop k iteration for geomechanic module
    void GaussSeidelAccelerationGeo(int k);
    
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

#endif /* TPMRSSegregatedAnalysis_h */
