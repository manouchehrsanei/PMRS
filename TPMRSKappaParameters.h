//
//  TPMRSKappaParameters.h
//  PMRS
//
//  Created by Omar Durán on 10/15/18.
//

#ifndef TPMRSKappaParameters_h
#define TPMRSKappaParameters_h

#include <stdio.h>
#include <map>
#include <iostream>
#include <string>
#include <vector>
#include "pzreal.h"
#include "pzerror.h"

/// Enumerate defining the porosity model
enum EModel { e_none = 0, e_constant = 1, e_petunin = 2, e_davies = 3 };

// Map to associate the strings with the enum values
static std::map<std::string, EModel> m_name_to_e_model;

class TPMRSKappaParameters {
    
private:
    
    /// model definition
    EModel m_model;
    
    /// model parameters
    std::vector<REAL> m_parameters;
    
    /// Initialization for names_to_e_model map
    static void Initialize();
    
    
public:
    
    /// Default constructor
    TPMRSKappaParameters();
    
    /// Copy constructor
    TPMRSKappaParameters(const TPMRSKappaParameters & other);
    
    /// Assignement constructor
    const TPMRSKappaParameters & operator=(const TPMRSKappaParameters & other);
    
    /// Desconstructor
    virtual ~TPMRSKappaParameters();
    
    /// Class name
    const std::string Name() const;
    
    /// Print class attributes
    virtual void Print(std::ostream &out = std::cout) const;
    
    /// Configurate the model according to the string and number of parameters
    void ConfigurateModel(std::string model, std::vector<REAL> & parameters);
    
    /// Computes the porosity using the selected model
    void Permeability(REAL &kappa, REAL &dkappa_dphi, REAL &kappa_0, REAL &phi, REAL &phi_0);
    
};

#endif /* TPMRSKappaParameters_h */
