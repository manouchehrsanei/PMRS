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

#include "TPMRSSimulationData.h"
#include "TPZElasticResponse.h"
#include "TPZElasticCriterion.h"


class TPMRSRKSolver {
    
private:
    
    /// Pointer of Simulation data
    TPMRSSimulationData * m_simulation_data;

    
    /// Material dimension
    int m_dimension;
    
public:
    
    /// Default constructor
    TPMRSRKSolver();
    
    /// Destructor
    ~TPMRSRKSolver();
    
    /// Copy constructor
    TPMRSRKSolver(const TPMRSRKSolver & other);
    
    /// Class name
    const std::string Name() const;
    
    /// Print class attributes
    virtual void Print(std::ostream &out = std::cout) const;
    
    /// Compute the analytic solution of displacement over space
    REAL DisplacementAnalytic(REAL &r, REAL &young, REAL &nu, REAL &biotcof, REAL &u_w, REAL &u_r, REAL &sig_0, REAL &p_0, REAL &p_w, REAL &r_o, REAL &r_w);
    
    
    /// Configurate sigma
    void ConfigureSigma();
    
    /// Configurate displacement
    void ConfigureDisplcement();
    
    /// Configurate pore pressure
    void ConfigurePressure();
    
    /// Configurate flux
    void ConfigureFlux();
    
    /// Runge Kutta solver
    void RungeKuttaSolver();


};



#endif /* TPMRSRKSolver_h */
