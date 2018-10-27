//
//  TPMRSMemoryPoroElast.h
//  PZ
//
//  Created by Omar and Manouchehr on 9/6/16.
//
//

#ifndef TPMRSMemoryPoroElast_h
#define TPMRSMemoryPoroElast_h

#include <stdio.h>
#include "pzreal.h"
#include "pzfmatrix.h"
#include "TPZElastoPlasticMem.h"


class TPMRSMemoryPoroElast
{
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  Memory : It is written for TPMRSMemoryPoroElast
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /// Basis function
    
    /// Brief Gradient of deformation at at n (last) state
    TPZFNMatrix<9,REAL> m_grad_u_n;
    
   /// Brief elastic strain at n (last) state
    TPZFMatrix<REAL> m_epsilon_e_n;

   /// Brief plastic strain at n (last) state
    TPZFMatrix<REAL> m_epsilon_p_n;
    
    
    
public:
    
    /// Brief Default constructor
    TPMRSMemoryPoroElast();
    
    /// Brief Default destructor
    ~TPMRSMemoryPoroElast();
    
    TPMRSMemoryPoroElast(const TPMRSMemoryPoroElast &copy)
    {
        
        m_grad_u_n    = copy.m_grad_u_n;
        m_epsilon_e_n = copy.m_epsilon_e_n;
        m_epsilon_p_n = copy.m_epsilon_p_n;

    }
    
    TPMRSMemoryPoroElast &operator=(const TPMRSMemoryPoroElast &other)
    {
        
        m_grad_u_n    = other.m_grad_u_n;
        m_epsilon_e_n = other.m_epsilon_e_n;
        m_epsilon_p_n = other.m_epsilon_p_n;
        
        return *this;
    }
    
    
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
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  Memory :
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    /// Brief Write method
    void Write(TPZStream &buf, int withclassid) const
    {
//                buf.Write(&m_Pressure_n);
//                buf.Write(&m_Pressure);
        DebugStop();
    }
    
    /// Brief Read method
    void Read(TPZStream &buf, void *context)
    {
        //        buf.Read(&m_Pressure_n);
        //        buf.Read(&m_Pressure);
        DebugStop();
    }
    
    void Print(std::ostream &out) const
    {
        //        out << m_Pressure_n;
        //        out << m_Pressure;
        DebugStop();
    }

    
    /// Brief Class unique identifier
    virtual int ClassId() const;
    
};

inline std::ostream &operator<<(std::ostream &out,const TPMRSMemoryPoroElast &mem)
{
    mem.Print(out);
    return out;
}


#endif /* TPMRSMemoryPoroElast_h */
