//
//  TPMRSPhiParameters.h
//  PMRS
//
//  Created by Omar Durán on 10/15/18.
//

#ifndef TPMRSPhiParameters_h
#define TPMRSPhiParameters_h

#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include "pzreal.h"
#include "pzerror.h"

#endif /* TPMRSPhiParameters_h */


class TPMRSPhiParameters {
    
private:
    
    /// Enumerate defining the porosity model
    enum EModel { e_none = 0, e_linear = 1 };
    
    /// model parameters
    std::vector<REAL> m_parameters;
    
    /// model definition
    EModel m_model;
    
public:
    
    /// Default constructor
    TPMRSPhiParameters();
    
    /// Copy constructor
    TPMRSPhiParameters(const TPMRSPhiParameters & other);
    
    /// Assignement constructor
    const TPMRSPhiParameters & operator=(const TPMRSPhiParameters & other);
    
    /// Desconstructor
    virtual ~TPMRSPhiParameters();
    
    /// Class name
    const std::string Name() const;
    
    /// Print class attributes
    virtual void Print(std::ostream &out = std::cout) const;
    
    /// Configurate the model according to the string and number of parameters
    void ConfigurateModel(std::string model, std::vector<REAL> & parameters);
    
};
