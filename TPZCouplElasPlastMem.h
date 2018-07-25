//
//  TPZCouplElasPlastMem.h
//  PZ
//
//  Created by Manouchehr on Jun 27, 2018.
//
//

#ifndef TPZCouplElasPlastMem_H
#define TPZCouplElasPlastMem_H

#include "TPZMaterial.h"
#include "TPZTensor.h"
#include "TPZPlasticState.h"
#include "TPZElastoPlasticMem.h"


  /**
   * This class defines the material memory that should be stored at each integration point
   * for the purposes of an elastoplastic material.
   */

class TPZCouplElasPlastMem :public TPZElastoPlasticMem
{  
    
    /** @brief Gradient of deformation at at n (last) state */
    TPZFNMatrix<9,REAL> m_grad_u_n;
    
   /** @brief elastic strain at n (last) state */
    TPZFMatrix<REAL> m_epsilon_e_n;

   /** @brief plastic strain at n (last) state */
    TPZFMatrix<REAL> m_epsilon_p_n;
    


public:
    TPZCouplElasPlastMem();
    
    TPZCouplElasPlastMem(const TPZCouplElasPlastMem & source);
    
    const TPZCouplElasPlastMem & operator=(const TPZCouplElasPlastMem & source);
    
    virtual ~TPZCouplElasPlastMem();
    
    const std::string Name()const;
    
    public:
virtual int ClassId() const;

    
    virtual void Write(TPZStream &buf, int withclassid) const;

    void Read(TPZStream &buf, void *context);

    virtual void Print(std::ostream &out = std::cout)const;
    
    /**
     * Operator<<
     */
    friend std::ostream& operator<<( std::ostream& Out, const TPZCouplElasPlastMem & s )
    {
        s.Print(Out);
        return Out;
    }
    



        /** @brief Set Gradient of deformation at n (last) state */
    void Set_grad_u_n(TPZFMatrix<STATE> & grad_u_n)
    {
        m_grad_u_n = grad_u_n;
    }
    
    /** @brief Get Gradient of deformation at n (last) state */
    TPZFMatrix<STATE> & grad_u_n()
    {
        return m_grad_u_n;
    }
    
    
    /** @brief Set elastic strain at n (last) state */
    void Set_epsilon_e_n(TPZFMatrix<REAL> & epsilon_e_n)
    {
        m_epsilon_e_n = epsilon_e_n;
    }
    
    /** @brief Get elastic strain at n (last) state */
    TPZFMatrix<REAL> epsilon_e_n()
    {
        return m_epsilon_e_n;
    }
    
    
    /** @brief Set plastic strain at n (last) state */
    void Set_epsilon_p_n(TPZFMatrix<REAL> & epsilon_p_n)
    {
        m_epsilon_p_n = epsilon_p_n;
    }
    
    /** @brief Get plastic strain at n (last) state */
    TPZFMatrix<REAL> epsilon_p_n(){
        return m_epsilon_p_n;
    }
    
};


#endif /* TPZCouplElasPlastMem_H */
