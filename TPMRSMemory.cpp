//
//  TPMRSMemory.cpp
//  PMRS
//
//  Created by Omar and Manouchehr on 9/11/18.
//

#include "TPMRSMemory.h"


TPMRSMemory::TPMRSMemory() : TPMRSMonoPhasicMemory() , TPMRSElastoPlasticMemory() {
    m_alpha = 0.0;
    m_Se = 0.0;
}

TPMRSMemory::TPMRSMemory(const TPMRSMemory & other): TPMRSMonoPhasicMemory(other), TPMRSElastoPlasticMemory(other) {
    m_alpha = other.m_alpha;
    m_Se    = other.m_Se;
}

const TPMRSMemory & TPMRSMemory::operator=(const TPMRSMemory & other) {
    
    /// check for self-assignment
    if(&other == this){
        return *this;
    }
    
    m_alpha = other.m_alpha;
    m_Se    = other.m_Se;
    
    return *this;
}

TPMRSMemory::~TPMRSMemory(){
    
}

const std::string TPMRSMemory::Name() const {
    return "TPMRSMemory";
}

void TPMRSMemory::Write(TPZStream &buf, int withclassid) const {
    TPMRSMonoPhasicMemory::Write(buf, withclassid);
    TPMRSElastoPlasticMemory::Write(buf, withclassid);
}

void TPMRSMemory::Read(TPZStream &buf, void *context){
    TPMRSMonoPhasicMemory::Read(buf, context);
    TPMRSElastoPlasticMemory::Read(buf, context);
}

void TPMRSMemory::Print(std::ostream &out) const {
    TPMRSMonoPhasicMemory::Print(out);
    TPMRSElastoPlasticMemory::Print(out);
}

int TPMRSMemory::ClassId() const {
    return Hash("TPMRSMemory");
}
