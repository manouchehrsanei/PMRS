//
//  TPMRSKappaParameters.cpp
//  PMRS
//
//  Created by Omar Durán on 10/15/18.
//

#include "TPMRSKappaParameters.h"


TPMRSKappaParameters::TPMRSKappaParameters(){
    
    m_model = k_none;
    m_parameters.resize(0);
    Initialize();
}

TPMRSKappaParameters::TPMRSKappaParameters(const TPMRSKappaParameters & other){
    m_model         = other.m_model;
    m_parameters    =   other.m_parameters;
    m_name_to_kappa_model = other.m_name_to_kappa_model;
}

const TPMRSKappaParameters & TPMRSKappaParameters::operator=(const TPMRSKappaParameters & other){
    
    // check for self-assignment
    if(&other == this){
        return *this;
    }
    
    m_model         = other.m_model;
    m_parameters    =   other.m_parameters;
    m_name_to_kappa_model = other.m_name_to_kappa_model;
    return *this;
}


TPMRSKappaParameters::~TPMRSKappaParameters(){
    
}


const std::string TPMRSKappaParameters::Name() const{
    return "TPMRSKappaParameters";
}


void TPMRSKappaParameters::Print(std::ostream &out) const{
    out << Name() << std::endl;
    out << "m_model = " << m_model << std::endl;
    for (int i  = 0; i < m_parameters.size(); i++) {
        out << "parameter number " << i << " = " << m_parameters[i] << std::endl;
    }
}


void TPMRSKappaParameters::SetModel(std::string model){
    
    switch (m_name_to_kappa_model[model])
    {
        case k_constant : {
            m_model = k_constant;
        }
            break;
        case k_petunin : {
            m_model = k_petunin;
        }
            break;
        case k_davies : {
            m_model = k_davies;
        }
            break;
        default : {
            DebugStop();
        }
    }
    
}

void TPMRSKappaParameters::Initialize()
{
    m_name_to_kappa_model["none"] = k_none;
    m_name_to_kappa_model["Petunin"] = k_petunin;
    m_name_to_kappa_model["Davies"] = k_davies;
}

void TPMRSKappaParameters::Permeability(REAL &kappa, REAL &dkappa_dphi, REAL &kappa_0, REAL &phi, REAL &phi_0){
    
    switch (m_model)
    {
        case k_constant : {
            kappa = kappa_0;
            dkappa_dphi = 0;
        }
            break;
        case k_petunin : {
            REAL A = m_parameters[0];
            kappa = kappa_0*pow(phi/phi_0,A);
            dkappa_dphi = (A*kappa_0*pow(phi/phi_0,-1.0 + A))/phi_0;
        }
            break;
        case k_davies : {
            REAL C = m_parameters[0];
            kappa = exp(C*(-1.0 + phi/phi_0))*kappa_0;
            dkappa_dphi = (C*exp(C*(-1.0 + phi/phi_0))*kappa_0)/phi_0;
        }
            break;
        default : {
            DebugStop();
        }
    }
    
}

void TPMRSKappaParameters::SetParameters(std::vector<REAL> parameters){
    m_parameters = parameters;
}

std::vector<REAL> & TPMRSKappaParameters::GetParameters(){
    return m_parameters;
}

TPMRSKappaParameters::EKappaModel TPMRSKappaParameters::GetModel(){
    return m_model;
}
