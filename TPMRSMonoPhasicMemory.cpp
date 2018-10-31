//
//  TPMRSMonoPhasicMemory.cpp
//  PMRS
//
//  Created by Omar and Manouchehr on 9/11/18.
//

#include "TPMRSMonoPhasicMemory.h"


TPMRSMonoPhasicMemory::TPMRSMonoPhasicMemory(){
    
    m_p_0       = 0.0;
    m_p         = 0.0;
    m_p_n       = 0.0;
    m_kappa_0   = 1.0e-13;
    m_kappa     = 1.0e-13;
    m_phi_0     = 0.1;
    m_phi       = 0.1;
}

TPMRSMonoPhasicMemory::TPMRSMonoPhasicMemory(const TPMRSMonoPhasicMemory & other){
    
    m_p_0       = other.m_p_0;
    m_p         = other.m_p;
    m_p_n       = other.m_p_n;
    m_kappa_0   = other.m_kappa_0;
    m_kappa     = other.m_kappa;
    m_phi_0     = other.m_phi_0;
    m_phi       = other.m_phi;
}

const TPMRSMonoPhasicMemory & TPMRSMonoPhasicMemory::operator=(const TPMRSMonoPhasicMemory & other){
    
    /// check for self-assignment
    if(&other == this){
        return *this;
    }
    
    m_p_0       = other.m_p_0;
    m_p         = other.m_p;
    m_p_n       = other.m_p_n;
    m_kappa_0   = other.m_kappa_0;
    m_kappa     = other.m_kappa;
    m_phi_0     = other.m_phi_0;
    m_phi       = other.m_phi;
    
    return *this;
}

TPMRSMonoPhasicMemory::~TPMRSMonoPhasicMemory(){
    
}

const std::string TPMRSMonoPhasicMemory::Name() const{
    return "TPMRSMonoPhasicMemory";
}

void TPMRSMonoPhasicMemory::Write(TPZStream &buf, int withclassid) const {
    buf.Write(&m_p_0);
    buf.Write(&m_p);
    buf.Write(&m_p_n);
    buf.Write(&m_kappa_0);
    buf.Write(&m_kappa);
    buf.Write(&m_phi_0);
    buf.Write(&m_phi);
}


void TPMRSMonoPhasicMemory::Read(TPZStream &buf, void *context){
    buf.Read(&m_p_0);
    buf.Read(&m_p);
    buf.Read(&m_p_n);
    buf.Read(&m_kappa_0);
    buf.Read(&m_kappa);
    buf.Read(&m_phi_0);
    buf.Read(&m_phi);
}

void TPMRSMonoPhasicMemory::Print(std::ostream &out) const{
    out << Name();
    out << "\n Initial pressure = " << m_p_0;
    out << "\n Pressure at last state = " << m_p;
    out << "\n Pressure at current state = " << m_p_n;
    out << "\n Initial absolute permeability = " << m_kappa_0;
    out << "\n Current absolute permeability = " << m_kappa;
    out << "\n Initial porosity = " << m_phi_0;
    out << "\n Current porosity = " << m_phi;
}

int TPMRSMonoPhasicMemory::ClassId() const{
    return Hash("TPMRSMonoPhasicMemory");
}
