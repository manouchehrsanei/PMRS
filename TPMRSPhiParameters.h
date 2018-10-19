//
//  TPMRSPhiParameters.h
//  PMRS
//
//  Created by Omar and Manouchehr on 10/15/18.
//

#ifndef TPMRSPhiParameters_h
#define TPMRSPhiParameters_h

#include <stdio.h>
#include <map>
#include <iostream>
#include <string>
#include <vector>
#include "pzreal.h"
#include "pzerror.h"

class TPMRSPhiParameters {
    
public:
    
    /// Enumerate defining the porosity model
    enum EPhiModel { p_none = 0, p_constant = 1, p_linear = 2 };
    
private:
    
    /// model definition
    EPhiModel m_model;
    
    /// model parameters
    std::vector<REAL> m_parameters;

    // Map to associate the strings with the enum values
    std::map<std::string, EPhiModel> m_name_to_p_model;
    
    /// Rock bulk modulus
    REAL m_K;
    
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
    
    /// Initialization for names_to_p_model map
    void Initialize();
    
    /// Configurate the model according to the string
    void SetModel(std::string model);
    
    /// Set the parameters
    void SetParameters(std::vector<REAL> parameters);
    
    /// Get the parameters
    std::vector<REAL> & GetParameters();
    
    /// Get the model index
    EPhiModel GetModel();
    
    /// Computes the porosity using the selected model
    void Porosity(REAL &phi, REAL &dphi_dp, REAL &phi_0, REAL &p, REAL &p_0, REAL &sigma_v, REAL &sigma_v_0, REAL &alpha, REAL &Se);
    
    /// Set the bulk modulus
    void SetBulkModulus(REAL E, REAL nu){
        REAL K = E/(3.0*(1.0-2.0*nu));
        m_K = K;
    }
    
    /// Get the bulk modulus
    REAL GetBulkModulus(){
        return m_K;
    }
    
};

#endif /* TPMRSPhiParameters_h */
