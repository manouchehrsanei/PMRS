//
//  TPMRSElastoPlasticMemory.cpp
//  PMRS
//
//  Created by Omar Durán on 9/11/18.
//

#include "TPMRSElastoPlasticMemory.h"


TPMRSElastoPlasticMemory::TPMRSElastoPlasticMemory(){
    
    m_sigma_n.Zero();
    m_u_n.resize(0);
    m_sigma.Zero();
    m_u.resize(0);
}

TPMRSElastoPlasticMemory::TPMRSElastoPlasticMemory(const TPMRSElastoPlasticMemory & other){
    
    m_sigma_n           = other.m_sigma_n;
    m_plastic_strain_n  = other.m_plastic_strain_n;
    m_u_n               = other.m_u_n;
    m_sigma             = other.m_sigma;
    m_plastic_strain    = other.m_plastic_strain;
    m_u                 = other.m_u;
}

const TPMRSElastoPlasticMemory & TPMRSElastoPlasticMemory::operator=(const TPMRSElastoPlasticMemory & other){
    
    // check for self-assignment
    if(&other == this){
        return *this;
    }
    
    m_sigma_n           = other.m_sigma_n;
    m_plastic_strain_n  = other.m_plastic_strain_n;
    m_u_n               = other.m_u_n;
    m_sigma             = other.m_sigma;
    m_plastic_strain    = other.m_plastic_strain;
    m_u                 = other.m_u;
    
    return *this;
}

TPMRSElastoPlasticMemory::~TPMRSElastoPlasticMemory(){
    
}

const std::string TPMRSElastoPlasticMemory::Name() const{
    return "TPMRSElastoPlasticMemory";
}

void TPMRSElastoPlasticMemory::Write(TPZStream &buf, int withclassid) const {
    m_sigma_n.Write(buf, withclassid);
    m_plastic_strain_n.Write(buf, withclassid);
    buf.Write(m_u_n);
    m_sigma.Write(buf, withclassid);
    m_plastic_strain.Write(buf, withclassid);
    buf.Write(m_u);
}


void TPMRSElastoPlasticMemory::Read(TPZStream &buf, void *context){
    m_sigma_n.Read(buf, context);
    m_plastic_strain_n.Read(buf, context);
    buf.Read(m_u_n);
    m_sigma.Read(buf, context);
    m_plastic_strain.Read(buf, context);
    buf.Read(m_u);
}

void TPMRSElastoPlasticMemory::Print(std::ostream &out) const{
    out << Name();
    out << "\n Current state stress = " << m_sigma_n;
    out << "\n Current Plastic strasin state = " << m_plastic_strain_n;
    out << "\n Current displacement field = " << m_u_n;
    out << "\n Last state stress = " << m_sigma;
    out << "\n Last Plastic strasin state = " << m_plastic_strain;
    out << "\n Last displacement field = " << m_u;
}

int TPMRSElastoPlasticMemory::ClassId() const{
    return Hash("TPMRSElastoPlasticMemory");
}
