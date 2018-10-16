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

/// Enumerate defining the plastic model
enum EEPModel { ep_none = 0, ep_mc = 1, ep_ds = 2, ep_cc = 3, ep_dp = 4 };

// Map to associate the strings with the enum values
static std::map<std::string, EEPModel> m_name_to_ep_model;

class TPMRSPlasticityParameters {
    
private:
    
    /// model definition
    EEPModel m_model;
    
    /// model parameters
    std::vector<REAL> m_parameters;
    
    /// Initialization for names_to_ep_model map
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
    
    /// Configurate the model according to the string and number of parameters
    void ConfigurateModel(std::string model, std::vector<REAL> & parameters);
    
    /// Set the parameters
    void SetParameters(std::vector<REAL> parameters);
    
    /// Get the parameters
    std::vector<REAL> & GetParameters();
    
    /// Get the model index
    EEPModel GetModel();
    
};


#endif /* TPMRSPlasticityParameters_h */
