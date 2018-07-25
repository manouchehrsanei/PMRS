//
//  TPZPMRSMemoryEP.cpp
//  PZ
//
//  Created by Manouchehr on Jun 27, 2018.
//
//

#include "TPZCouplElasPlastMem.h"


TPZCouplElasPlastMem::TPZCouplElasPlastMem(): TPZElastoPlasticMem()
{ 

    /** @brief Gradient of deformation at n (last) state */
    m_grad_u_n.Resize(3, 3);
    m_grad_u_n.Zero();
    
    
    /** @brief Elastic strain at n (last) state */
    m_epsilon_e_n.Resize(3, 3);
    m_epsilon_e_n.Zero();
    
    
    /** @brief Plastic strain at n (last) state */
    m_epsilon_p_n.Resize(3, 3);
    m_epsilon_p_n.Zero();

}
    
TPZCouplElasPlastMem::TPZCouplElasPlastMem(const TPZCouplElasPlastMem & source):  TPZElastoPlasticMem(source)
{
}

TPZCouplElasPlastMem::~TPZCouplElasPlastMem(){ }

void TPZCouplElasPlastMem::Write(TPZStream &buf, int withclassid) const
{
    TPZElastoPlasticMem::Write(buf,withclassid);

}

void TPZCouplElasPlastMem::Read(TPZStream &buf, void *context)
{
    TPZElastoPlasticMem::Read(buf,context);
}

void TPZCouplElasPlastMem::Print(std::ostream &out)const
{
    out << Name();
    out << "\n Parent Class Data:";
    TPZElastoPlasticMem::Print(out);
}

const std::string TPZCouplElasPlastMem::Name()const
{
    return "TPZCouplElasPlastMem";   
}

int TPZCouplElasPlastMem::ClassId() const{
    return Hash("TPZCouplElasPlastMem") ^ TPZElastoPlasticMem::ClassId() << 1;
}

const TPZCouplElasPlastMem & TPZCouplElasPlastMem::operator=(const TPZCouplElasPlastMem & source)
{
    TPZElastoPlasticMem::operator=(source);
    return *this;
}
