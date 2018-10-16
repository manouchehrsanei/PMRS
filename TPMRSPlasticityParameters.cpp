//
//  TPMRSPlasticityParameters.cpp
//  PMRS
//
//  Created by Omar Durán on 10/14/18.
//

#include "TPMRSPlasticityParameters.h"

TPMRSPlasticityParameters::TPMRSPlasticityParameters(){
    
    m_model         = ep_none;
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
    m_name_to_ep_model["none"] = ep_none;
    m_name_to_ep_model["MC"] = ep_mc;
    m_name_to_ep_model["DS"] = ep_ds;
    m_name_to_ep_model["CC"] = ep_cc;
    m_name_to_ep_model["DP"] = ep_dp;
    
//    std::cout << "TPMRSPlasticityParameters::Initialization." << std::endl;
}

void TPMRSPlasticityParameters::ConfigurateModel(std::string model, std::vector<REAL> & parameters){
    
    switch (m_name_to_ep_model[model])
    {
        case ep_mc : {
            m_model = ep_mc;
            if (parameters.size()!=2) {
                DebugStop();
            }
        }
            break;
        case ep_ds : {
            m_model = ep_ds;
            if (parameters.size()!=8) {
                DebugStop();
            }
        }
            break;
        default : {
            DebugStop();
        }
    }
    m_parameters = parameters;
}

void TPMRSPlasticityParameters::SetParameters(std::vector<REAL> parameters){
    m_parameters = parameters;
}

std::vector<REAL> & TPMRSPlasticityParameters::GetParameters(){
    return m_parameters;
}

EEPModel TPMRSPlasticityParameters::GetModel(){
    return m_model;
    
}
