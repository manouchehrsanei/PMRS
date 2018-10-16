//
//  TPMRSPlasticityParameters.h
//  PMRS
//
//  Created by Omar Durán on 10/14/18.
//

#ifndef TPMRSPlasticityParameters_h
#define TPMRSPlasticityParameters_h

#include <stdio.h>
#include <map>
#include <iostream>
#include <string>
#include <vector>
#include "pzreal.h"
#include "pzerror.h"

/// Enumerate defining the porosity model
enum EModel { e_none = 0, e_mc = 1, e_ds = 2, e_cc = 3, e_dp = 4 };

// Map to associate the strings with the enum values
static std::map<std::string, EModel> m_name_to_e_model;

class TPMRSPlasticityParameters {
    
private:
    
    /// model definition
    EModel m_model;
    
    /// model parameters
    std::vector<REAL> m_parameters;
    
    /// Initialization for names_to_e_model map
    static void Initialize();
    
    
public:
    
    /// Default constructor
    TPMRSPlasticityParameters();
    
    /// Copy constructor
    TPMRSPlasticityParameters(const TPMRSPlasticityParameters & other);
    
    /// Assignement constructor
    const TPMRSPlasticityParameters & operator=(const TPMRSPlasticityParameters & other);
    
    /// Desconstructor
    virtual ~TPMRSPlasticityParameters();
    
    /// Class name
    const std::string Name() const;
    
    /// Print class attributes
    virtual void Print(std::ostream &out = std::cout) const;
    
};


#endif /* TPMRSPlasticityParameters_h */
