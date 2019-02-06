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
    m_q.Resize(3,0.);
    m_q_n.Resize(3,0.);
    m_kappa_0   = 1.0e-13;
    m_kappa_n   = 1.0e-13;
    m_phi_0     = 0.1;
    m_delta_phi       = 0.1;
    m_phi_n     = 0.1;
}

TPMRSMonoPhasicMemory::TPMRSMonoPhasicMemory(const TPMRSMonoPhasicMemory & other){
    
    m_p_0       = other.m_p_0;
    m_p         = other.m_p;
    m_p_n       = other.m_p_n;
    m_q         = other.m_q;
    m_q_n       = other.m_q_n;
    m_kappa_0   = other.m_kappa_0;
    m_kappa_n   = other.m_kappa_n;
    m_phi_0     = other.m_phi_0;
    m_delta_phi       = other.m_delta_phi;
    m_phi_n     = other.m_phi_n;
}

const TPMRSMonoPhasicMemory & TPMRSMonoPhasicMemory::operator=(const TPMRSMonoPhasicMemory & other){
    
    /// check for self-assignment
    if(&other == this){
        return *this;
    }
    
    m_p_0       = other.m_p_0;
    m_p         = other.m_p;
    m_p_n       = other.m_p_n;
    m_q         = other.m_q;
    m_q_n       = other.m_q_n;
    m_kappa_0   = other.m_kappa_0;
    m_kappa_n   = other.m_kappa_n;
    m_phi_0     = other.m_phi_0;
    m_delta_phi       = other.m_delta_phi;
    m_phi_n     = other.m_phi_n;
    
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
    buf.Write(m_q);
    buf.Write(m_q_n);
    buf.Write(&m_kappa_0);
    buf.Write(&m_kappa_n);
    buf.Write(&m_phi_0);
    buf.Write(&m_delta_phi);
    buf.Write(&m_phi_n);
}


void TPMRSMonoPhasicMemory::Read(TPZStream &buf, void *context){
    buf.Read(&m_p_0);
    buf.Read(&m_p);
    buf.Read(&m_p_n);
    buf.Read(m_q);
    buf.Read(m_q_n);
    buf.Read(&m_kappa_0);
    buf.Read(&m_kappa_n);
    buf.Read(&m_phi_0);
    buf.Read(&m_delta_phi);
    buf.Read(&m_phi_n);
}

void TPMRSMonoPhasicMemory::Print(std::ostream &out) const{
    out << Name();
    out << "\n Initial pressure              = " << m_p_0;
    out << "\n Pressure at last state        = " << m_p;
    out << "\n Pressure at current state     = " << m_p_n;
    out << "\n Last flux vector              = " << m_q;
    out << "\n Current flux vector           = " << m_q_n;
    out << "\n Initial absolute permeability = " << m_kappa_0;
    out << "\n Current absolute permeability = " << m_kappa_n;
    out << "\n Initial porosity              = " << m_phi_0;
    out << "\n Porosity correction           = " << m_delta_phi;
    out << "\n Current porosity              = " << m_phi_n;
}

int TPMRSMonoPhasicMemory::ClassId() const{
    return Hash("TPMRSMonoPhasicMemory");
}
