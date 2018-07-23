//
//  TPZPMRSMemory.cpp
//  PZ
//
//  Created by Omar and Manouchehr on 9/6/16.
//
//

#include "TPZPMRSMemory.h"

/** @brief Default constructor */
TPZPMRSMemory::TPZPMRSMemory():  fSigma() {
    
   
    
    /** @brief Gradient of deformation at n (last) state */
    m_grad_u_n.Resize(3, 3);
    m_grad_u_n.Zero();
    
    
    /** @brief Elastic strain at n (last) state */
    m_epsilon_e_n.Resize(3, 3);
    m_epsilon_e_n.Zero();
    
    
    /** @brief Plastic strain at n (last) state */
    m_epsilon_p_n.Resize(3, 3);
    m_epsilon_p_n.Zero();
    
    /** @brief Gradient of pore pressure at n (last) state */
    m_grad_p_n.Resize(3, 3);
    m_grad_p_n.Zero();
    
    /** @brief Porosity at n (last) state */
    m_phi_n.Resize(3, 3);
    m_phi_n.Zero();
    
    /** @brief Permeability at n (last) state */
    m_k_n.Resize(3, 3);
    m_k_n.Zero();

    
}

/** @brief Default destructor */
TPZPMRSMemory::~TPZPMRSMemory(){
    
}

int TPZPMRSMemory::ClassId() const{
    return Hash("TPZPMRSMemory");
}

