//
//  TPMRSPoroMechParameters.cpp
//  PMRS
//
//  Created by Omar Durán on 10/14/18.
//

#include "TPMRSPoroMechParameters.h"

TPMRSPoroMechParameters::TPMRSPoroMechParameters(){
    
    m_e_young = 0.0;
    m_nu = 0.0;
    m_alpha = 0.0;
    m_se = 0.0;
    m_eta = 0.0;
    m_rho_f = 0.0;
    m_rho_s = 0.0;
}


TPMRSPoroMechParameters::TPMRSPoroMechParameters(const TPMRSPoroMechParameters & other){
    
    m_e_young = other.m_e_young;
    m_nu = other.m_nu;
    m_alpha = other.m_alpha;
    m_se = other.m_se;
    m_eta = other.m_eta;
    m_rho_f = other.m_rho_f;
    m_rho_s = other.m_rho_s;
}

const TPMRSPoroMechParameters & TPMRSPoroMechParameters::operator=(const TPMRSPoroMechParameters & other){
    
    // check for self-assignment
    if(&other == this){
        return *this;
    }
    
    m_e_young = other.m_e_young;
    m_nu = other.m_nu;
    m_alpha = other.m_alpha;
    m_se = other.m_se;
    m_eta = other.m_eta;
    m_rho_f = other.m_rho_f;
    m_rho_s = other.m_rho_s;
    
    return *this;
}


TPMRSPoroMechParameters::~TPMRSPoroMechParameters(){
    
}

const std::string  TPMRSPoroMechParameters::Name() const{
    return "TPMRSPoroMechParameters";
}


void  TPMRSPoroMechParameters::Print(std::ostream &out) const{
    out << Name() << std::endl;
    out << "m_e_young = " << m_e_young << std::endl;
    out << "m_nu = " << m_nu << std::endl;
    out << "m_alpha = " << m_alpha << std::endl;
    out << "m_se =" << m_se << std::endl;
    out << "m_eta =" << m_eta << std::endl;
    out << "m_rho_f =" << m_rho_f << std::endl;
    out << "m_rho_s =" << m_rho_s << std::endl;
}

