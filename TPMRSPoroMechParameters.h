//
//  TPMRSPoroMechParameters.h
//  PMRS
//
//  Created by Omar and Manouchehr on 10/14/18.
//

#ifndef TPMRSPoroMechParameters_h
#define TPMRSPoroMechParameters_h

#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include "pzreal.h"
#include "pzerror.h"

class TPMRSPoroMechParameters {
    
private:
    
    /// model parameters
    std::vector<REAL> m_parameters;
    
public:
    
    /// Default constructor
    TPMRSPoroMechParameters();
    
    /// Copy constructor
    TPMRSPoroMechParameters(const TPMRSPoroMechParameters & other);
    
    /// Assignement constructor
    const TPMRSPoroMechParameters & operator=(const TPMRSPoroMechParameters & other);
    
    /// Destructor
    virtual ~TPMRSPoroMechParameters();
    
    /// Class name
    const std::string Name() const;
    
    /// Print class attributes
    virtual void Print(std::ostream &out = std::cout) const;
    
    /// Set parameters vector {Eyoung, nu, ks, c_f, eta, rho_f, rho_s}
    void SetParameters(std::vector<REAL> & parameters);
    
    /// Get parameters vector {Eyoung, nu, ks, c_f, eta, rho_f, rho_s}
    std::vector<REAL> & GetParameters();
    
};

#endif /* TPMRSPoroMechParameters_h */
