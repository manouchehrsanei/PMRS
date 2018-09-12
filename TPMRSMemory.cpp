//
//  TPMRSMemory.cpp
//  PMRS
//
//  Created by Omar Durán on 9/9/18.
//

#include "TPMRSMemory.h"


TPMRSMemory::TPMRSMemory() : TPMRSMonoPhasicMemory() , TPMRSElastoPlasticMemory() {
    
}

TPMRSMemory::TPMRSMemory(const TPMRSMemory & other): TPMRSMonoPhasicMemory(other), TPMRSElastoPlasticMemory(other) {
    
}

const TPMRSMemory & TPMRSMemory::operator=(const TPMRSMemory & other) {
    
    // check for self-assignment
    if(&other == this){
        return *this;
    }
    
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
