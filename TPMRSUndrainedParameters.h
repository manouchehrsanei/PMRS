//
//  TPMRSUndrainedParameters.h
//  PMRS
//
//  Created by Omar Durán on 10/14/18.
//

#ifndef TPMRSUndrainedParameters_h
#define TPMRSUndrainedParameters_h

#include <stdio.h>
#include <iostream>
#include <string>
#include "pzreal.h"

class TPMRSUndrainedParameters {
  
private:
    
    /// Undrained Young modulus
    REAL m_e_young_u;
    
    /// Undrained Poisson ratio
    REAL m_nu_u;
    
    /// Initial constant porosity
    REAL m_phi_0;
    
    /// Initial constant abolute permeability
    REAL m_kappa_0;
    
public:
    
    /// Default constructor
    TPMRSUndrainedParameters();
    
    /// Copy constructor
    TPMRSUndrainedParameters(const TPMRSUndrainedParameters & other);
    
    /// Assignement constructor
    const TPMRSUndrainedParameters & operator=(const TPMRSUndrainedParameters & other);
    
    /// Desconstructor
    virtual ~TPMRSUndrainedParameters();
    
    /// Class name
    const std::string Name() const;
    
    /// Print class attributes
    virtual void Print(std::ostream &out = std::cout) const;
    
    ///  Set undrained parameters
    void SetUndrainedParameters(REAL e_young_u, REAL nu_u, REAL phi_0, REAL kappa_0);

};

#endif /* TPMRSUndrainedParameters_h */
