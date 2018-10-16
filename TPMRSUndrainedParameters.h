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
#include <vector>
#include "pzreal.h"
#include "pzerror.h"

class TPMRSUndrainedParameters {
  
private:
    
    /// Undrained parameters
    std::vector<REAL> m_parameters;
    
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
    
    /// Set parameters vector {Eyoung_u, nu_u, phi_0, kappa_0}
    void SetParameters(std::vector<REAL> & parameters);
    
    /// Get parameters vector {Eyoung_u, nu_u, phi_0, kappa_0}
    std::vector<REAL> & GetParameters();

};

#endif /* TPMRSUndrainedParameters_h */
