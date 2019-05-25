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
#include "TPMRSMemory.h"

#include "TPMRSKappaParameters.h"


template <class T, class TMEM = TPMRSMemory>
class TPMRSRKSolver {
    
private:
    
    /// Number of spatial steps
    int m_n_steps;
    
    /// Initial state for the vector of states variables
    std::vector<REAL> m_y_0;
    
    /// Defines the permeability model
    TPMRSKappaParameters m_kappa_models;

    /// Wellbore region radius
    REAL m_re;
    
    /// Wellbore radius
    REAL m_rw;
    
    /// Fluid viscosity
    REAL m_eta;
    
    /// Fluid density
    REAL m_rho;
    
    /// Fluid compressibility
    REAL m_cf;
    
    /// Grain bulk modulus
    REAL m_K_s;
    
    /// Spatial step size
    REAL m_dr;
    
    /// Default memory item
    TMEM m_default_memory;
    
    /// Directive that stands for fourth order Runge-Kutta approximation
    bool m_is_RK4_Q;
    
    /// Directive that stands for first order Euler approximation
    bool m_is_Euler_Q;
    
    /// Directive that stands for where y_0 is defined
    bool m_is_Re_Q;
    
    /// Number of state variables
    int m_n_state;
    
    /// Plastic integrator object composed by a yield function and a elastic predictor
    T m_plastic_integrator;
    
    /// Vector of memory items
    std::vector<TMEM> m_memory;
    
    /// Array that storage the Runge-Kutta approximated solution
    TPZFMatrix<REAL> m_r_y;
    
    /// Vector of first Lame parameter
    std::vector<REAL> m_lambda;
    
    /// Vector of second Lame parameter
    std::vector<REAL> m_mu;
    
    /// Effecttive stress tensor associated to the initial data
    TPZTensor<REAL> m_sigma;
    
    /// Elastic strain tensor associated to the initial data
    TPZTensor<REAL> m_eps_e;
    
    /// Plastic strain tensor associated to the initial data
    TPZTensor<REAL> m_eps_p;
    
    /// Directive to load memory vector entry
    bool m_accept_solution_Q;
    
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
    
    /// Set the fluid data
    void SetFluidData(REAL eta, REAL cf, REAL rho){
        m_eta = eta;
        m_cf  = cf;
        m_rho = rho;
    }

    
    /// Set the Grain bulk modulus
    void SetGrainBulkModulus(REAL K_s){
        m_K_s = K_s;
    }
    
    /// Set the discretization parameters
    void SetDiscretization(REAL rw, REAL re, int n_steps){
        m_rw = rw;
        m_re = re;
        m_n_steps = n_steps;
        if(m_is_Re_Q){
            m_dr = (m_rw - m_re)/REAL(n_steps);
        }else{
            m_dr = (m_re - m_rw)/REAL(n_steps);
        }
    }
    
    /// Set the initial data parameters
    void SetInitialData(std::vector<REAL> y_0){
        m_y_0 = y_0;
    }

    /// Set the initial data parameters
    void SetElastoPlasticInitialData(TPZTensor<REAL> & eps_e, TPZTensor<REAL> & eps_p, TPZTensor<REAL> & sigma){
        m_eps_e = eps_e;
        m_eps_p = eps_p;
        m_sigma = sigma;
    }
    
    /// Set the Default memory item
    void SetDefaultMemory(TMEM & default_memory){
        m_default_memory = default_memory;
    }
    
    /// Set the directive that stands for fourth Runge-Kutta approximation
    void SetFourthOrderApproximation(){
        m_is_RK4_Q = true;
    }
    
    /// Set the directive that stands for first order Euler approximation
    void SetFirstOrderApproximation(){
        m_is_Euler_Q = true;
    }
    
    /// Set the directive that stands for first order Euler approximation
    void SetDefineDataAtRe(){
        m_is_Re_Q = true;
    }
    
    /// Set the Permeability Parameters
    void SetKappaParameters(TPMRSKappaParameters kappa_models){
        m_kappa_models = kappa_models;
    }
    
    /// Get the Permeability Parameters
    TPMRSKappaParameters GetKappaParameters(){
        return m_kappa_models;
    }
    
    
    /// Synchronize all the information
    void Synchronize();
    
    /// Perform the Runge-Kutta approximation
    void ExecuteRKApproximation();
    
    /// Print the Runge-Kutta approximation
    void PrintRKApproximation(std::ostream &out = std::cout);
    
    /// Print the secondary variables (s_r,s_t,eps_t_r,eps_t_t,eps_p_r,eps_p_t,phi,kappa)
    void PrintSecondaryVariables(std::ostream &out = std::cout);
    
private:
    
    /// Right hand side function
    std::vector<REAL> f(int i, REAL & r, std::vector<REAL> & y);
    
    /// Reconstruc an approximated eps for coherence with mem and accept the point approximation
    void ReconstructAndAcceptPoint(int i, REAL & r, std::vector<REAL> & y, bool initial_data_Q = false);
    
    /// Euler method with first order accuracy
    std::vector<REAL> EulerApproximation(int i, REAL & r, std::vector<REAL> & y);
    
    /// Runge-Kutta method with second order accuracy
    std::vector<REAL> RK2Approximation(int i, REAL & r, std::vector<REAL> & y);
    
    /// Runge-Kutta method with fourth order accuracy
    std::vector<REAL> RK4Approximation(int i, REAL & r, std::vector<REAL> & y);
    
    /// Function that computes total strain state
    TPZTensor<REAL> Epsilon(int i, REAL & r, std::vector<REAL> & y);
    
    /// Function that computes stress state
    TPZTensor<REAL> Sigma(int i, TPZTensor<REAL> & epsilon, TPZFMatrix<REAL> * Dep);
    
    /// First lame parameter
    REAL lambda(int i);
    
    /// Second lame parameter
    REAL mu(int i);
    
    /// Bulk modulus parameter
    REAL K(REAL & lambda, REAL & mu);
    
    /// Biot coefficient function
    REAL Alpha(REAL & K);
    
    /// Porosity function
    REAL Porosity(int i, REAL & r, std::vector<REAL> & y);
    
    /// Permeability function
    REAL Permeability(int i, REAL & phi);
    
    /// Append solution to m_r_y
    void AppendTo(int i, std::vector<REAL> y);
    
    /// Zaxpy for std::vector
    std::vector<REAL> a_times_v(const REAL a, std::vector<REAL> & v);
    
    /// Sumation of two vector a and b
    std::vector<REAL> a_add_b(std::vector<REAL> & a, std::vector<REAL> & b);

};



#endif /* TPMRSRKSolver_h */
