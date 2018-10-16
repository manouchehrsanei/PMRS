//
//  TPMRSKappaParameters.cpp
//  PMRS
//
//  Created by Omar Durán on 10/15/18.
//

#include "TPMRSKappaParameters.h"


TPMRSKappaParameters::TPMRSKappaParameters(){
    
    m_model = e_none;
    m_parameters.resize(0);
    Initialize();
}

TPMRSKappaParameters::TPMRSKappaParameters(const TPMRSKappaParameters & other){
    m_model         = other.m_model;
    m_parameters    =   other.m_parameters;
    Initialize();
}

const TPMRSKappaParameters & TPMRSKappaParameters::operator=(const TPMRSKappaParameters & other){
    
    // check for self-assignment
    if(&other == this){
        return *this;
    }
    
    m_model         = other.m_model;
    m_parameters    =   other.m_parameters;
    Initialize();
    
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


void TPMRSKappaParameters::ConfigurateModel(std::string model, std::vector<REAL> & parameters){
    
    switch (m_name_to_e_model[model])
    {
        case e_constant : {
            m_model = e_constant;
        }
            break;
        case e_petunin : {
            m_model = e_petunin;
            if (parameters.size()!=1) {
                DebugStop();
            }
        }
            break;
        case e_davies : {
            m_model = e_davies;
            if (parameters.size()!=1) {
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

void TPMRSKappaParameters::Initialize()
{
    m_name_to_e_model["none"] = e_none;
    m_name_to_e_model["Petunin"] = e_petunin;
    m_name_to_e_model["Davies"] = e_davies;
    
    std::cout << "TPMRSPhiParameters::Initialization." << std::endl;
}

void TPMRSKappaParameters::Permeability(REAL &kappa, REAL &dkappa_dphi, REAL &kappa_0, REAL &phi, REAL &phi_0){
    
    switch (m_model)
    {
        case e_constant : {
            kappa = kappa_0;
            dkappa_dphi = 0;
        }
            break;
        case e_petunin : {
            REAL A = m_parameters[0];
            kappa = kappa_0*pow(phi/phi_0,A);
            dkappa_dphi = (A*kappa_0*pow(phi/phi_0,-1.0 + A))/phi_0;
        }
            break;
        case e_davies : {
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
