//
//  TPMRSPlasticityParameters.cpp
//  PMRS
//
//  Created by Omar Durán on 10/14/18.
//

#include "TPMRSPlasticityParameters.h"

TPMRSPlasticityParameters::TPMRSPlasticityParameters(){
    
    m_model         = e_none;
    m_parameters.resize(0);
    Initialize();
    
}

TPMRSPlasticityParameters::TPMRSPlasticityParameters(const TPMRSPlasticityParameters & other){
    m_model         = other.m_model;
    m_parameters    = other.m_parameters;
    Initialize();
    
}

const TPMRSPlasticityParameters & TPMRSPlasticityParameters::operator=(const TPMRSPlasticityParameters & other){
    
    // check for self-assignment
    if(&other == this){
        return *this;
    }
    
    m_model         = other.m_model;
    m_parameters    = other.m_parameters;
    Initialize();
    return *this;
    
}

TPMRSPlasticityParameters::~TPMRSPlasticityParameters(){
    
}

const std::string TPMRSPlasticityParameters::Name() const{
    return "TPMRSPlasticityParameters";
}

void TPMRSPlasticityParameters::Print(std::ostream &out) const {
    out << Name() << std::endl;
    out << "m_model = " << m_model << std::endl;
    for (int i  = 0; i < m_parameters.size(); i++) {
        out << "parameter number " << i << " = " << m_parameters[i] << std::endl;
    }
}

void TPMRSPlasticityParameters::Initialize()
{
    m_name_to_e_model["none"] = e_none;
    m_name_to_e_model["MC"] = e_mc;
    m_name_to_e_model["DS"] = e_ds;
    m_name_to_e_model["CC"] = e_cc;
    m_name_to_e_model["DP"] = e_dp;
    
    std::cout << "TPMRSPlasticityParameters::Initialization." << std::endl;
}
