//
//  TPMRSCoupPoElaMemory.h
//  PZ
//
//  Created by Omar and Manouchehr on 9/6/16.
//
//

#ifndef TPMRSCoupPoElaMemory_h
#define TPMRSCoupPoElaMemory_h

#include <stdio.h>
#include "pzreal.h"
#include "pzfmatrix.h"
#include "TPZElastoPlasticMem.h"


class TPMRSCoupPoElaMemory
{
    
private:
    /// Brief Gradient of deformation at n (last) state
    TPZFNMatrix<9,REAL> m_grad_u_n;
    
   /// Brief elastic strain at n (last) state
    TPZFMatrix<REAL> m_epsilon_e_n;

   /// Brief plastic strain at n (last) state
    TPZFMatrix<REAL> m_epsilon_p_n;
    
    
    
public:
    
    /// Default constructor
    TPMRSCoupPoElaMemory();
    
    /// Copy constructor
    TPMRSCoupPoElaMemory(const TPMRSCoupPoElaMemory & other);
    
    /// Assignement constructor
    const TPMRSCoupPoElaMemory & operator=(const TPMRSCoupPoElaMemory & other);
    
    /// Default destructor
    ~TPMRSCoupPoElaMemory();
    
    /// Class Name
    const std::string Name() const;
    
    /// Write class attributes
    virtual void Write(TPZStream &buf, int withclassid) const;
    
    /// Read class attributes
    virtual void Read(TPZStream &buf, void *context);
    
    /// Print class attributes
    virtual void Print(std::ostream &out = std::cout) const;
    
    /// Print class attributes
    friend std::ostream & operator<<( std::ostream& out, const TPMRSCoupPoElaMemory & memory ){
        memory.Print(out);
        return out;
    }
    
    /// Brief Class unique identifier
    virtual int ClassId() const;
    
    
    /// ********* Set and Get some functions ********************************************

    
    /// Brief Set Gradient of deformation at n (last) state
    void Set_grad_u_n(TPZFMatrix<STATE> & grad_u_n)
    {
        m_grad_u_n = grad_u_n;
    }
    
    /// Brief Get Gradient of deformation at n (last) state
    TPZFMatrix<STATE> & grad_u_n()
    {
        return m_grad_u_n;
    }
    
    
    /// Brief Set elastic strain at n (last) state
    void Set_epsilon_e_n(TPZFMatrix<REAL> & epsilon_e_n)
    {
        m_epsilon_e_n = epsilon_e_n;
    }
    
    /// Brief Get elastic strain at n (last) state
    TPZFMatrix<REAL> epsilon_e_n()
    {
        return m_epsilon_e_n;
    }
    
    /// Brief Set plastic strain at n (last) state
    void Set_epsilon_p_n(TPZFMatrix<REAL> & epsilon_p_n)
    {
        m_epsilon_p_n = epsilon_p_n;
    }
    
    /// Brief Get plastic strain at n (last) state
    TPZFMatrix<REAL> epsilon_p_n(){
        
        return m_epsilon_p_n;
    }
    
};


#endif /* TPMRSCoupPoElaMemory_h */
