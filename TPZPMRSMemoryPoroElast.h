//
//  TPZPMRSMemoryPoroElast.h
//  PZ
//
//  Created by Omar and Manouchehr on 9/6/16.
//
//

#ifndef TPZPMRSMemoryPoroElast_h
#define TPZPMRSMemoryPoroElast_h

#include <stdio.h>
#include "pzreal.h"
#include "pzfmatrix.h"
#include "TPZElastoPlasticMem.h"


class TPZPMRSMemoryPoroElast
{
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  Memory : It is written for TPZPMRSMemoryPoroElast
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // Basis functions
    TPZElastoPlasticMem m_Mem;
    
    
    /** @brief Gradient of deformation at at n (last) state */
    TPZFNMatrix<9,REAL> m_grad_u_n;
    
   /** @brief elastic strain at n (last) state */
    TPZFMatrix<REAL> m_epsilon_e_n;

   /** @brief plastic strain at n (last) state */
    TPZFMatrix<REAL> m_epsilon_p_n;
    
    
    
public:
    
    /** @brief Default constructor */
    TPZPMRSMemoryPoroElast();
    
    /** @brief Default destructor */
    ~TPZPMRSMemoryPoroElast();
    
    TPZPMRSMemoryPoroElast(const TPZPMRSMemoryPoroElast &copy)
    {
        
        m_grad_u_n    = copy.m_grad_u_n;
        m_epsilon_e_n = copy.m_epsilon_e_n;
        m_epsilon_p_n = copy.m_epsilon_p_n;

    }
    
    TPZPMRSMemoryPoroElast &operator=(const TPZPMRSMemoryPoroElast &other)
    {
        
        m_grad_u_n    = other.m_grad_u_n;
        m_epsilon_e_n = other.m_epsilon_e_n;
        m_epsilon_p_n = other.m_epsilon_p_n;
        
        return *this;
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
    
    

    
    /**
     * Total (elastic+plastic) stress
     */
    TPZTensor<REAL> fSigma;
    
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  Memory :
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    /**
     Write method not implemented.
     
     @param buf TPZStream buffer
     @param withclassid obsolete
     */
    void Write(TPZStream &buf, int withclassid) const
    {
//                buf.Write(&m_Pressure_n);
//                buf.Write(&m_Pressure);
        DebugStop();
    }
    
    /**
     Read method not implemented.
     
     @param buf TPZStream buffer
     @param context obsolete
     */
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

    
    /**
     Class unique identifier

     @return integer class id
     */
    virtual int ClassId() const;
    
};

inline std::ostream &operator<<(std::ostream &out,const TPZPMRSMemoryPoroElast &mem)
{
    mem.Print(out);
    return out;
}


#endif /* TPZPMRSMemoryPoroElast_h */
