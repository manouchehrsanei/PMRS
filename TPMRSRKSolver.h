//
//  TPMRSRKSolver.h
//  PMRS
//
//  Created by Manouchehr Sanei and Omar on 3/8/19.
//
//

#ifndef TPMRSRKSolver_h
#define TPMRSRKSolver_h

#include <stdio.h>
#include "pzerror.h"
#include "pzfmatrix.h"
#include "TPZTensor.h"
#include "TPMRSSimulationData.h"
#include "TPZElasticResponse.h"
#include "TPZElasticCriterion.h"


template <class T, class TMEM>
class TPMRSRKSolver {
    
private:
    
    /// Pointer of Simulation data
    TPMRSSimulationData * m_simulation_data;
    
    /// Material dimension
    int m_dimension;
    
    /// Number of state variables
    int m_n_state;
    
    /// Number of spatial steps
    int m_n_steps;
    
    /// Initial state for the vector of states variables
    std::vector<REAL> y_0;
    
    /// Spatial step size
    REAL m_dr;
    
    /// Wellbore region radius
    REAL m_re;
    
    /// Wellbore radius
    REAL m_rw;
    
    /// Array that storage the Runge-Kutta approximated solution
    TPZFMatrix<REAL> m_r_y;
    
    /// Plastic integrator object composed by a yield function and a elastic predictor
    T m_plastic_integrator;
    
    /// Fluid viscosity
    REAL m_eta;
    
public:
    
    /// Default constructor
    TPMRSRKSolver();
    
    /// Destructor
    ~TPMRSRKSolver();
    
    /// Copy constructor
    TPMRSRKSolver(const TPMRSRKSolver & other);
    
    /// Assignment constructor
    TPMRSRKSolver & operator=(const TPMRSRKSolver & other);
    
    /// Class name
    const std::string Name() const;
    
    /// Print class attributes
    void Print(std::ostream &out = std::cout) const;
    
    /// Set the plastic integrator
    void SetPlasticIntegrator(T & plastic_integrator){
        m_plastic_integrator = plastic_integrator;
    }
    
    /// Get the plastic integrator
    T & GetPlasticIntegrator(){
        return m_plastic_integrator;
    }
    
    /// Set the fluid viscosity
    void SetFluidViscosity(REAL eta){
        m_eta = eta;
    }
    
    /// Get the fluid viscosity
    REAL GetFluidViscosity(){
        return m_eta;
    }
    
    /// Set the reservoir dimensions
    void SetReservoirDimensions(REAL rw, REAL re){
        m_rw = rw;
        m_re = re;
    }
    
    /// Set the reservoir dimensions
    void SetNumberOfSteps(int n_steps){
        m_n_steps = ;
    }
    
    /// Right hand side function
    std::vector<REAL> f(std::vector<REAL> & y);
    
    /// Runge-Kutta method with second order accuracy
    std::vector<REAL> RK2Approximation(std::vector<REAL> & y_n);
    
    /// Runge-Kutta method with fourth order accuracy
    std::vector<REAL> RK4Approximation(std::vector<REAL> & y_n);
    
    /// Perform the Runge-Kutta approximation
    void ExecuteRKApproximation();
    
    /// Print the Runge-Kutta approximation
    void PrintRKApproximation();
    
    /// Function that computes total strain state
    TPZTensor<REAL> Epsilon(std::vector<REAL> & y, std::vector<REAL> & y_n);
    
    /// Function that computes stress state
    TPZTensor<REAL> Sigma(TPZTensor<REAL> & epsilon, TPZFMatrix<REAL> * Dep);
    
    /// Bulk modulus function
    REAL K(std::vector<REAL> & y, std::vector<REAL> & y_n);
    
    /// Bulk modulus function
    REAL alpha(REAL & K, REAL & K_s);
    
    /// Porosity function
    REAL Porosity(std::vector<REAL> & y, std::vector<REAL> & y_n);
    
    /// Permeability function
    REAL Permeability(REAL & phi);

};



#endif /* TPMRSRKSolver_h */
