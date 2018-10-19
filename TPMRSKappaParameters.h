//
//  TPMRSKappaParameters.h
//  PMRS
//
//  Created by Omar and Manouchehr on 10/15/18.
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

class TPMRSKappaParameters {
    
public:
    
    /// Enumerate defining the permeability model
    enum EKappaModel { k_none = 0, k_constant = 1, k_petunin = 2, k_davies = 3 };
    
private:
    
    /// model definition
    EKappaModel m_model;
    
    /// model parameters
    std::vector<REAL> m_parameters;
    
    // Map to associate the strings with the enum values
    std::map<std::string, EKappaModel> m_name_to_kappa_model;
    
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
    
    /// Initialization for names_to_k_model map
    void Initialize();
    
    /// Configurate the model according to the string
    void SetModel(std::string model);
    
    /// Set the parameters
    void SetParameters(std::vector<REAL> parameters);
    
    /// Get the parameters
    std::vector<REAL> & GetParameters();
    
    /// Get the model index
    EKappaModel GetModel();
    
    /// Computes the porosity using the selected model
    void Permeability(REAL &kappa, REAL &dkappa_dphi, REAL &kappa_0, REAL &phi, REAL &phi_0);
    
};

#endif /* TPMRSKappaParameters_h */
