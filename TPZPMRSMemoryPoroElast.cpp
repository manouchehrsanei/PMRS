//
//  TPZPMRSMemoryPoroElast.cpp
//  PZ
//
//  Created by Omar and Manouchehr on 9/6/16.
//
//

#include "TPZPMRSMemoryPoroElast.h"

/** @brief Default constructor */
TPZPMRSMemoryPoroElast::TPZPMRSMemoryPoroElast():  fSigma() {
    
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

/** @brief Default destructor */
TPZPMRSMemoryPoroElast::~TPZPMRSMemoryPoroElast(){
    
}

int TPZPMRSMemoryPoroElast::ClassId() const{
    return Hash("TPZPMRSMemoryPoroElast");
}

