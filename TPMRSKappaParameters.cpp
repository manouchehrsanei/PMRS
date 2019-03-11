//
//  TPMRSKappaParameters.cpp
//  PMRS
//
//  Created by Omar and Manouchehr on 10/15/18.
//

#include "TPMRSKappaParameters.h"
#include <math.h>


TPMRSKappaParameters::TPMRSKappaParameters(){
    
    m_model = k_none;
    m_parameters.resize(0);
    Initialize();
}

TPMRSKappaParameters::TPMRSKappaParameters(const TPMRSKappaParameters & other){
    m_model               = other.m_model;
    m_parameters          = other.m_parameters;
    m_name_to_kappa_model = other.m_name_to_kappa_model;
}

const TPMRSKappaParameters & TPMRSKappaParameters::operator=(const TPMRSKappaParameters & other){
    
    /// check for self-assignment
    if(&other == this){
        return *this;
    }
    
    m_model               = other.m_model;
    m_parameters          = other.m_parameters;
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
    for (int i  = 0; i < m_parameters.size(); i++)
    {
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
        case k_costa : {
            m_model = k_costa;
        }
            break;
        case k_nelson : {
            m_model = k_nelson;
        }
            break;
        case k_bayles : {
            m_model = k_bayles;
        }
            
            break;
            
        default : {
            DebugStop();
        }
    }
    
}

void TPMRSKappaParameters::Initialize()
{
    m_name_to_kappa_model["none"]     = k_none;
    m_name_to_kappa_model["Constant"] = k_constant;
    m_name_to_kappa_model["Petunin"]  = k_petunin;
    m_name_to_kappa_model["Davies"]   = k_davies;
    m_name_to_kappa_model["Costa"]    = k_costa;
    m_name_to_kappa_model["Nelson"]   = k_nelson;
    m_name_to_kappa_model["Bayles"]   = k_bayles;
    
}

void TPMRSKappaParameters::Permeability(REAL &kappa, REAL &dkappa_dphi, REAL &kappa_0, REAL &phi, REAL &phi_0){
    
    switch (m_model)
    {
        case k_constant : {
            kappa       = kappa_0;
            dkappa_dphi = 0;
        }
            break;
        case k_petunin : {
            REAL A      = m_parameters[0];
            
            kappa       = kappa_0*pow(phi/phi_0,A);
            dkappa_dphi = (A*kappa_0*pow(phi/phi_0,-1.0 + A))/phi_0;
        }
            break;
        case k_davies : {
            REAL C      = m_parameters[0];
            
            kappa       = exp(C*(-1.0 + phi/phi_0))*kappa_0;
            dkappa_dphi = (C*exp(C*(-1.0 + phi/phi_0))*kappa_0)/phi_0;
        }
            break;
        case k_costa : {
            REAL A      = m_parameters[0];
            REAL C      = m_parameters[1];

            kappa       = kappa_0*A*((1-phi_0)/(1-phi))*pow(phi/phi_0,C);
            dkappa_dphi = (A*kappa_0*(1 - phi_0)*pow(phi/phi_0,C))/pow(1 - phi,2) +
            (A*C*kappa_0*(1 - phi_0)*pow(phi/phi_0,-1 + C))/((1 - phi)*phi_0);
        }
            break;
        case k_nelson : {
            REAL A      = m_parameters[0];
            REAL C      = m_parameters[1];
            
            kappa       = kappa_0*pow(10,(C + A*(phi - phi_0)));
            dkappa_dphi = A*kappa_0*log(10.0)*pow(10,(C + A*(phi - phi_0)));
        }
            break;
            
        case k_bayles : {
            REAL A      = m_parameters[0];
            REAL C      = m_parameters[1];
            
            kappa       = (A*kappa_0*pow(1 - phi_0,2)*pow(phi/phi_0,2 + C))/pow(1 - phi,2);
            dkappa_dphi = (2*A*kappa_0*pow(1 - phi_0,2)*pow(phi/phi_0,2 + C))/pow(1 - phi,3) +
            (A*(2 + C)*kappa_0*pow(1 - phi_0,2)*pow(phi/phi_0,1 + C))/(pow(1 - phi,2)*phi_0);
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
