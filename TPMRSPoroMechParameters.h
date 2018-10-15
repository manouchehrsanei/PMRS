//
//  TPMRSPoroMechParameters.h
//  PMRS
//
//  Created by Omar Durán on 10/14/18.
//

#ifndef TPMRSPoroMechParameters_h
#define TPMRSPoroMechParameters_h

#include <stdio.h>
#include <iostream>
#include <string>
#include "pzreal.h"

class TPMRSPoroMechParameters {
    
    // @TODO:: MS, please document this class.
private:
    
    REAL m_e_young;
    
    REAL m_nu;
    
    REAL m_alpha;
    
    REAL m_se;
    
    REAL m_eta;
    
    REAL m_rho_f;
    
    REAL m_rho_s;

    
public:
    
    /// Default constructor
    TPMRSPoroMechParameters();
    
    /// Copy constructor
    TPMRSPoroMechParameters(const TPMRSPoroMechParameters & other);
    
    /// Assignement constructor
    const TPMRSPoroMechParameters & operator=(const TPMRSPoroMechParameters & other);
    
    /// Desconstructor
    virtual ~TPMRSPoroMechParameters();
    
    /// Class name
    const std::string Name() const;
    
    /// Print class attributes
    virtual void Print(std::ostream &out = std::cout) const;
    
    
};

#endif /* TPMRSPoroMechParameters_h */
