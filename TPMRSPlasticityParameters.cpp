//
//  TPMRSPlasticityParameters.cpp
//  PMRS
//
//  Created by Omar and Manouchehr on 10/14/18.
//

#include "TPMRSPlasticityParameters.h"

TPMRSPlasticityParameters::TPMRSPlasticityParameters(){
    
    m_model         = ep_none;
    m_parameters.resize(0);
    Initialize();
    
}

TPMRSPlasticityParameters::TPMRSPlasticityParameters(const TPMRSPlasticityParameters & other){
    m_model            = other.m_model;
    m_parameters       = other.m_parameters;
    m_name_to_ep_model = other.m_name_to_ep_model;
    
}

const TPMRSPlasticityParameters & TPMRSPlasticityParameters::operator=(const TPMRSPlasticityParameters & other){
    
    /// check for self-assignment
    if(&other == this){
        return *this;
    }
    
    m_model            = other.m_model;
    m_parameters       = other.m_parameters;
    m_name_to_ep_model = other.m_name_to_ep_model;
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
    for (int i  = 0; i < m_parameters.size(); i++)
    {
        out << "parameter number " << i << " = " << m_parameters[i] << std::endl;
    }
}

void TPMRSPlasticityParameters::Initialize()
{
    m_name_to_ep_model["none"] = ep_none;
    m_name_to_ep_model["MC"]   = ep_mc;
    m_name_to_ep_model["DS"]   = ep_ds;
    m_name_to_ep_model["CC"]   = ep_cc;
    m_name_to_ep_model["DP"]   = ep_dp;
}

void TPMRSPlasticityParameters::SetModel(std::string model){
    
    switch (m_name_to_ep_model[model])
    {
        case ep_mc : {
            m_model = ep_mc;
        }
            break;
        case ep_ds : {
            m_model = ep_ds;
        }
            break;
        default : {
            DebugStop();
        }
    }
}

void TPMRSPlasticityParameters::SetParameters(std::vector<REAL> parameters){
    m_parameters = parameters;
}

std::vector<REAL> & TPMRSPlasticityParameters::GetParameters(){
    return m_parameters;
}

TPMRSPlasticityParameters::EEPModel TPMRSPlasticityParameters::GetModel(){
    return m_model;
    
}
