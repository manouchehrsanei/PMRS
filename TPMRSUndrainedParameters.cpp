//
//  TPMRSUndrainedParameters.cpp
//  PMRS
//
//  Created by Omar Durán on 10/14/18.
//

#include "TPMRSUndrainedParameters.h"


TPMRSUndrainedParameters::TPMRSUndrainedParameters(){
    
    m_e_young_u = 0.0;
    m_nu_u      = 0.0;
    m_phi_0     = 0.0;
    m_kappa_0   = 0.0;
    
}

TPMRSUndrainedParameters::TPMRSUndrainedParameters(const TPMRSUndrainedParameters & other){
    
    m_e_young_u = other.m_e_young_u;
    m_nu_u      = other.m_nu_u;
    m_phi_0     = other.m_phi_0;
    m_kappa_0   = other.m_kappa_0;
}

const TPMRSUndrainedParameters & TPMRSUndrainedParameters::operator=(const TPMRSUndrainedParameters & other){
    
    // check for self-assignment
    if(&other == this){
        return *this;
    }
    
    m_e_young_u = other.m_e_young_u;
    m_nu_u      = other.m_nu_u;
    m_phi_0     = other.m_phi_0;
    m_kappa_0   = other.m_kappa_0;
    
    return *this;
}

TPMRSUndrainedParameters::~TPMRSUndrainedParameters(){
    
}

const std::string TPMRSUndrainedParameters::Name() const {
    return "TPMRSUndrainedParameters";
}

void TPMRSUndrainedParameters::Print(std::ostream &out) const {
    out << Name() << std::endl;
    out << "m_e_young_u = " << m_e_young_u << std::endl;
    out << "m_nu_u = " << m_nu_u << std::endl;
    out << "m_phi_0 = " << m_phi_0 << std::endl;
    out << "m_kappa_0 =" << m_kappa_0 << std::endl;
}

void TPMRSUndrainedParameters::SetUndrainedParameters(REAL e_young_u, REAL nu_u, REAL phi_0, REAL kappa_0){
    m_e_young_u = e_young_u;
    m_nu_u      = nu_u;
    m_phi_0     = phi_0;
    m_kappa_0   = kappa_0;
}

REAL TPMRSUndrainedParameters::E_young_u(){
    return m_e_young_u;
}

REAL TPMRSUndrainedParameters::nu_u(){
    return m_nu_u;
}

REAL TPMRSUndrainedParameters::phi_0(){
    return m_phi_0;
}

REAL TPMRSUndrainedParameters::kappa_0(){
    return m_kappa_0;
}

