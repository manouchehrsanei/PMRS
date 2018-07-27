//
//  TPZPMRSCouplPoroPlast.cpp
//  PZ
//
//  Created by Manouchehr on Jun 27, 2018.
//
//


#include "TPZPMRSCouplPoroPlast.h"
#include <iostream>
#include <string>
#include "pzelasmat.h"
#include "pzbndcond.h"
#include "pzaxestools.h"
#include "pzintel.h"
#include "TPZElasticResponse.h"


#include "pzfmatrix.h"
#include "TPZTensor.h"


#include "TPZCouplElasPlastMem.h"


#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.TPZPMRSCouplPoroPlast"));
#endif



/** @brief default costructor */
template<class T,class TMEM>
TPZPMRSCouplPoroPlast<T,TMEM>::TPZPMRSCouplPoroPlast():TPZMatElastoPlastic2D<T,TMEM>()
{
    m_Dim = 3;
    m_b.resize(3);
    
    m_b[0]=0.;
    m_b[1]=0.;
    m_b[2]=0.;
    
    m_k_model = 0;
    
    m_SetRunPlasticity = false;
}

/** @brief costructor based on a material id */
template<class T,class TMEM>
TPZPMRSCouplPoroPlast<T,TMEM>::TPZPMRSCouplPoroPlast(int matid, int dim):TPZMatElastoPlastic2D<T,TMEM>(matid,1)
{
    m_Dim = dim;
    m_b.resize(3);
    
    m_b[0]=0.;
    m_b[1]=0.;
    m_b[2]=0.;

    m_k_model = 0;

    m_SetRunPlasticity = false;
}

/** @brief default destructor */
template<class T,class TMEM>
TPZPMRSCouplPoroPlast<T,TMEM>::~TPZPMRSCouplPoroPlast()
{
}


/** @brief copy constructor $ */
template<class T,class TMEM>
TPZPMRSCouplPoroPlast<T,TMEM>::TPZPMRSCouplPoroPlast(const TPZPMRSCouplPoroPlast& other)
{
    this->m_Dim    = other.m_Dim;
    this->m_SimulationData    = other.m_SimulationData;
}


/** @brief Copy assignemnt operator $ */
template<class T,class TMEM>
TPZPMRSCouplPoroPlast<T,TMEM> &TPZPMRSCouplPoroPlast<T,TMEM>::operator = (const TPZPMRSCouplPoroPlast& other)
{
    
    if (this != & other) // prevent self-assignment
    {
        this->m_Dim    = other.m_Dim;
        this->m_SimulationData    = other.m_SimulationData;
    }
    return *this;
}


/** @brief number of state variables */
template<class T,class TMEM>
int TPZPMRSCouplPoroPlast<T,TMEM>::NStateVariables()
{
    return 2;
}

/** @brief a computation of delta strain */
template <class T, class TMEM>
void TPZPMRSCouplPoroPlast<T,TMEM>::ComputeDeltaStrainVector(TPZMaterialData & data, TPZFMatrix<REAL> &DeltaStrain)
{
    TPZMatElastoPlastic2D<T,TMEM>::ComputeDeltaStrainVector(data, DeltaStrain);
}

/** @brief a computation of stress and tangent */
template <class T, class TMEM>
void TPZPMRSCouplPoroPlast<T,TMEM>::ApplyDeltaStrainComputeDep(TPZMaterialData & data, TPZFMatrix<REAL> & DeltaStrain,TPZFMatrix<REAL> & Stress, TPZFMatrix<REAL> & Dep)
{
#ifdef PZDEBUG
    if (DeltaStrain.Rows() != 6)
    {
        DebugStop();
    }
#endif
    TPZMatElastoPlastic<T,TMEM>::ApplyDeltaStrainComputeDep(data,DeltaStrain,Stress,Dep);
}

/** @brief a computation of stress */
template <class T, class TMEM>
void TPZPMRSCouplPoroPlast<T,TMEM>::ApplyDeltaStrain(TPZMaterialData & data, TPZFMatrix<REAL> & DeltaStrain,TPZFMatrix<REAL> & Stress)
{
#ifdef PZDEBUG
    if (DeltaStrain.Rows() != 6)
    {
        DebugStop();
    }
#endif
    TPZMatElastoPlastic2D<T,TMEM>::ApplyDeltaStrain(data,DeltaStrain,Stress);
}


template<class T,class TMEM>
void TPZPMRSCouplPoroPlast<T,TMEM>::SetUpdateToUseFullU(bool update)
{
    m_UpdateToUseFullDiplacement = update;
}



/** @brief permeability coupling models  */
template<class T,class TMEM>
REAL TPZPMRSCouplPoroPlast<T,TMEM>::k_permeability(REAL &phi, REAL &k)
{
    switch (m_k_model)
    {
        case 0:
        {
            k = m_k_0;
        }
            break;
            
            
        case 1: // Petunin et al. (2011), A = 2.0
        {
            k = m_k_0*pow((phi/m_porosity_0),2.0);
        }
            break;
            
            
        case 2: // Santos et al. (2014): Unloading/Reloading, Virgin Loading: A = 4.60
        {
            k = m_k_0*pow((phi/m_porosity_0),2.44);
        }
            break;
            
            
        case 3: // Santos et al. (2014): Unloading/Reloading, Virgin Loading: A = 7.19
        {
            k = m_k_0*pow((phi/m_porosity_0),4.62);
        }
            break;
        
        default:
        {
            DebugStop();
        }
            break;
    }
    

    return k;
}

/** @brief Poroelastic porosity correction */
template<class T,class TMEM>
REAL TPZPMRSCouplPoroPlast<T,TMEM>::porosity_corrected_2D(TPZVec<TPZMaterialData> &datavec)
{
    
    int u_b = 0;
    int p_b = 1;
    
    // Getting the space functions
    TPZFNMatrix <9,REAL>	&axes_u	=	datavec[u_b].axes;
    
    // Getting the solutions and derivatives
    TPZManVector<REAL,2> u = datavec[u_b].sol[0];
    TPZManVector<REAL,1> p = datavec[p_b].sol[0];
    
    TPZFNMatrix <6,REAL> du = datavec[u_b].dsol[0];
    TPZFNMatrix <6,REAL> dp = datavec[p_b].dsol[0];
    
    
    TPZFNMatrix<6,REAL> Grad_u(2,2,0.0);
    
    // Computing Gradient of the Solution
    Grad_u(0,0) = du(0,0)*axes_u(0,0)+du(1,0)*axes_u(1,0); // dux/dx
    Grad_u(1,0) = du(0,0)*axes_u(0,1)+du(1,0)*axes_u(1,1); // dux/dy
    
    Grad_u(0,1) = du(0,1)*axes_u(0,0)+du(1,1)*axes_u(1,0); // duy/dx
    Grad_u(1,1) = du(0,1)*axes_u(0,1)+du(1,1)*axes_u(1,1); // duy/dy
    
    REAL div_u = Grad_u(0,0) + Grad_u(1,1);
    REAL phi = m_porosity_0 + m_alpha * div_u + m_Se * p[0];
    
    return phi;

}


/** @brief Poroelastic porosity correction */
template<class T,class TMEM>
REAL TPZPMRSCouplPoroPlast<T,TMEM>::porosity_corrected_3D(TPZVec<TPZMaterialData> &datavec)
{
    
    int u_b = 0;
    int p_b = 1;
    
    // Getting the space functions
    TPZFNMatrix <9,REAL>	&axes_u	=	datavec[u_b].axes;
    
    // Getting the solutions and derivatives
    TPZManVector<REAL,3> u = datavec[u_b].sol[0];
    TPZManVector<REAL,1> p = datavec[p_b].sol[0];
    
    TPZFNMatrix <9,REAL> du = datavec[u_b].dsol[0];
    TPZFNMatrix <6,REAL> dp = datavec[p_b].dsol[0];
    
    TPZFNMatrix<9,REAL> Grad_u(3,3,0.0);
    
    // Computing Gradient of the Solution
    Grad_u(0,0) = du(0,0)*axes_u(0,0)+du(1,0)*axes_u(1,0)+du(2,0)*axes_u(2,0); // dux/dx
    Grad_u(1,0) = du(0,0)*axes_u(0,1)+du(1,0)*axes_u(1,1)+du(2,0)*axes_u(2,1); // dux/dy
    Grad_u(2,0) = du(0,0)*axes_u(0,2)+du(1,0)*axes_u(1,2)+du(2,0)*axes_u(2,2); // dux/dz
    
    Grad_u(0,1) = du(0,1)*axes_u(0,0)+du(1,1)*axes_u(1,0)+du(2,1)*axes_u(2,0); // duy/dx
    Grad_u(1,1) = du(0,1)*axes_u(0,1)+du(1,1)*axes_u(1,1)+du(2,1)*axes_u(2,1); // duy/dy
    Grad_u(2,1) = du(0,1)*axes_u(0,2)+du(1,1)*axes_u(1,2)+du(2,1)*axes_u(2,2); // duy/dz
    
    Grad_u(0,2) = du(0,2)*axes_u(0,0)+du(1,2)*axes_u(1,0)+du(2,2)*axes_u(2,0); // duz/dx
    Grad_u(1,2) = du(0,2)*axes_u(0,1)+du(1,2)*axes_u(1,1)+du(2,2)*axes_u(2,1); // duz/dy
    Grad_u(2,2) = du(0,2)*axes_u(0,2)+du(1,2)*axes_u(1,2)+du(2,2)*axes_u(2,2); // duz/dz
    
    REAL div_u = Grad_u(0,0) + Grad_u(1,1) + Grad_u(2,2);
    REAL phi = m_porosity_0 + m_alpha * div_u + m_Se * p[0];
    
    return phi;
    
}

/// Transform a voight notation to a tensor
template<class T,class TMEM>
void TPZPMRSCouplPoroPlast<T,TMEM>::FromVoight(TPZVec<STATE> &Svoight, TPZFMatrix<STATE> &S)
{
    S(0,0) = Svoight[TPZPMRSCouplPoroPlast<T,TMEM>::Exx];
    S(0,1) = Svoight[TPZPMRSCouplPoroPlast<T,TMEM>::Exy];
    S(0,2) = Svoight[TPZPMRSCouplPoroPlast<T,TMEM>::Exz];
    S(1,0) = Svoight[TPZPMRSCouplPoroPlast<T,TMEM>::Exy];
    S(1,1) = Svoight[TPZPMRSCouplPoroPlast<T,TMEM>::Eyy];
    S(1,2) = Svoight[TPZPMRSCouplPoroPlast<T,TMEM>::Eyz];
    S(2,0) = Svoight[TPZPMRSCouplPoroPlast<T,TMEM>::Exz];
    S(2,1) = Svoight[TPZPMRSCouplPoroPlast<T,TMEM>::Eyz];
    S(2,2) = Svoight[TPZPMRSCouplPoroPlast<T,TMEM>::Ezz];
}


/** @brief of contribute of BC */
template<class T,class TMEM>
void TPZPMRSCouplPoroPlast<T,TMEM>::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE>  &ek, TPZFMatrix<STATE> &ef)
{

    if (m_Dim == 3)
    {
        this->Contribute_3D(datavec, weight, ek, ef);
    }
    else
    {
        this->Contribute_2D(datavec, weight, ek, ef);
    }
    
}



/** @brief of contribute in 2 dimensional */
template<class T,class TMEM>
void TPZPMRSCouplPoroPlast<T,TMEM>::Contribute_2D(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE>  &ek, TPZFMatrix<STATE> &ef)
{
    
    int u_b = 0;
    int p_b = 1;
    
    // Getting the space functions
    TPZFMatrix<REAL>    &phiu   =   datavec[u_b].phi;
    TPZFMatrix<REAL>    &phip   =   datavec[p_b].phi;
    
    TPZFMatrix<REAL>    &grad_phi_u   =   datavec[u_b].dphix;
    TPZFMatrix<REAL>    &grad_phi_p   =   datavec[p_b].dphix;
    
    TPZFNMatrix <9,REAL>	&axes_u	=	datavec[u_b].axes;
    TPZFNMatrix <9,REAL>	&axes_p	=	datavec[p_b].axes;
    
    // Getting the solutions and derivatives
    TPZManVector<REAL,3> u = datavec[u_b].sol[0];
    TPZManVector<REAL,1> p = datavec[p_b].sol[0];
    
    TPZFNMatrix <6,REAL> du = datavec[u_b].dsol[0];
    TPZFNMatrix <6,REAL> dp = datavec[p_b].dsol[0];
    
    TPZFNMatrix<6,REAL> Grad_p(2,1,0.0),Grad_phi_i(2,1,0.0),Grad_phi_j(2,1,0.0);
    Grad_p(0,0) = dp(0,0)*axes_p(0,0)+dp(1,0)*axes_p(1,0);
    Grad_p(1,0) = dp(0,0)*axes_p(0,1)+dp(1,0)*axes_p(1,1);
    
    
    int nphi_u = phiu.Rows();
    int nphi_p = phip.Rows();
    
    int first_u = 0;
    int first_p = 2*nphi_u;
    
    // Compute porosity poroelastic correction
    REAL phi_poro = porosity_corrected_2D(datavec);
    
    
    // @brief of checking whether the time of diffusion is in the current state or not
    REAL dt = m_SimulationData->dt();
    if (!m_SimulationData->IsCurrentStateQ())
    {
        

        // Darcy mono-phascis flow
        for (int ip = 0; ip < nphi_p; ip++)
        {
            
            ef(ip + first_p, 0)		+= - weight * (phi_poro/dt) * phip(ip,0);
        }
        
        return;
    }
    
    REAL rho_avg = (1.0-phi_poro)*m_rho_s+phi_poro*m_rho_f;
    m_b[0] = rho_avg*m_SimulationData->Gravity()[0];
    m_b[1] = rho_avg*m_SimulationData->Gravity()[1];

    // Computing Gradient of the Solution
    TPZFNMatrix<6,REAL> S(3,3,0.0);

    TPZFNMatrix<6,REAL> Grad_vx_i(2,1,0.0);
    TPZFNMatrix<6,REAL> Grad_vy_i(2,1,0.0);

    TPZFNMatrix<6,REAL> Grad_v(2,2,0.0);
    TPZFNMatrix<6,REAL> Grad_vx_j(2,1,0.0);
    TPZFNMatrix<6,REAL> Grad_vy_j(2,1,0.0);

    TPZFMatrix<REAL> & Sigma_0 = m_SimulationData->PreStress();
    
    Sigma_0.Zero();
    
    S = Sigma_0; // Applying prestress

    for (int iu = 0; iu < nphi_u; iu++) {
        
        // Computing Gradient of the test function for each component
        Grad_vx_i(0,0) = grad_phi_u(0,iu)*axes_u(0,0)+grad_phi_u(1,iu)*axes_u(1,0); // dvx/dx
        Grad_vx_i(1,0) = grad_phi_u(0,iu)*axes_u(0,1)+grad_phi_u(1,iu)*axes_u(1,1); // dvx/dy
        
        Grad_vy_i(0,0) = grad_phi_u(0,iu)*axes_u(0,0)+grad_phi_u(1,iu)*axes_u(1,0); // dvy/dx
        Grad_vy_i(1,0) = grad_phi_u(0,iu)*axes_u(0,1)+grad_phi_u(1,iu)*axes_u(1,1); // dvy/dy
        
        ef(2*iu   + first_u, 0) += weight * ((S(0,0)-m_alpha * p[0])*Grad_vx_i(0,0)+S(0,1)*Grad_vx_i(1,0)-(m_b[0])*phiu(iu,0));
        ef(2*iu+1 + first_u, 0)	+= weight * (S(1,0)*Grad_vy_i(0,0)+(S(1,1)-m_alpha*p[0])*Grad_vy_i(1,0)-(m_b[1])*phiu(iu,0));
        
        
        for (int ju = 0; ju < nphi_u; ju++) {
            
           
            // Computing Gradient of the test function
            Grad_vx_j(0,0) = grad_phi_u(0,ju)*axes_u(0,0)+grad_phi_u(1,ju)*axes_u(1,0); // dvx/dx
            Grad_vx_j(1,0) = grad_phi_u(0,ju)*axes_u(0,1)+grad_phi_u(1,ju)*axes_u(1,1); // dvx/dy
            
            Grad_vy_j(0,0) = grad_phi_u(0,ju)*axes_u(0,0)+grad_phi_u(1,ju)*axes_u(1,0); // dvy/dx
            Grad_vy_j(1,0) = grad_phi_u(0,ju)*axes_u(0,1)+grad_phi_u(1,ju)*axes_u(1,1); // dvy/dy
            
            
            ek(2*iu   + first_u, 2*ju   + first_u)  += (1.)*weight*(((2.0*m_mu+m_lambda)*Grad_vx_j(0,0))*Grad_vx_i(0,0)+m_mu*Grad_vx_j(1,0)*Grad_vx_i(1,0));
            ek(2*iu   + first_u, 2*ju+1 + first_u)  += (1.)*weight*((m_lambda*Grad_vy_j(1,0))*Grad_vx_i(0,0)+m_mu*Grad_vy_j(0,0)*Grad_vx_i(1,0));
            ek(2*iu+1 + first_u, 2*ju   + first_u)	+= (1.)*weight*(m_mu*Grad_vx_j(1,0)*Grad_vy_i(0,0)+m_lambda*Grad_vx_j(0,0)*Grad_vy_i(1,0));
            ek(2*iu+1 + first_u, 2*ju+1 + first_u)	+= (1.)*weight*((2.0*m_mu+m_lambda)*Grad_vy_j(1,0)*Grad_vy_i(1,0)+m_mu*Grad_vy_j(0,0)*Grad_vy_i(0,0));
            
        }
        
    }
    
    TPZFNMatrix<6,REAL> dv(2,1,0.0);
    
    //	Matrix -Qc
    //	Coupling matrix
    for(int iu = 0; iu < nphi_u; iu++ )
    {
        
        // Computing Gradient of the test function for each component
        Grad_vx_i(0,0) = grad_phi_u(0,iu)*axes_u(0,0)+grad_phi_u(1,iu)*axes_u(1,0); // dvx/dx
        Grad_vx_i(1,0) = grad_phi_u(0,iu)*axes_u(0,1)+grad_phi_u(1,iu)*axes_u(1,1); // dvx/dy
        
        Grad_vy_i(0,0) = grad_phi_u(0,iu)*axes_u(0,0)+grad_phi_u(1,iu)*axes_u(1,0); // dvy/dx
        Grad_vy_i(1,0) = grad_phi_u(0,iu)*axes_u(0,1)+grad_phi_u(1,iu)*axes_u(1,1); // dvy/dy
        
        for(int jp = 0; jp < nphi_p; jp++)
        {
            
            ek(2*iu,first_p+jp) += (-1.)* weight * m_alpha * phip(jp,0) * Grad_vx_i(0,0);
            ek(2*iu+1,first_p+jp) += (-1.)* weight * m_alpha * phip(jp,0) * Grad_vy_i(1,0);
        }
    }
    
    //	Matrix QcˆT
    //	Coupling matrix transpose
    for(int ip = 0; ip < nphi_p; ip++ )
    {
        
        for(int ju = 0; ju < nphi_u; ju++)
        {
            
            dv(0,0) = grad_phi_u(0,ju)*axes_u(0,0)+grad_phi_u(1,ju)*axes_u(1,0);
            dv(1,0) = grad_phi_u(0,ju)*axes_u(0,1)+grad_phi_u(1,ju)*axes_u(1,1);
            
            ek(first_p+ip,2*ju  ) += (1./dt) * weight * m_alpha * dv(0,0) * phip(ip,0);
            ek(first_p+ip,2*ju+1) += (1./dt) * weight * m_alpha * dv(1,0) * phip(ip,0);
            
        }
    }
    
    /** @brief Rudnicki diffusion coefficient */
    /** J. W. Rudnicki. Fluid mass sources and point forces in linear elastic diffusive solids. Journal of Mechanics of Materials, 5:383–393, 1986. */
    REAL k = 0.0;
    m_k_model = 1;
    k_permeability(phi_poro,k);
    m_lambdau = 1.1 * m_lambda;
    
    REAL c = 1; // (k/m_eta); // (k/m_eta)*(m_lambdau-m_lambda)*(m_lambda + 2.0*m_mu)/(m_alpha*m_alpha*(m_lambdau + 2.0*m_mu));

    // Darcy mono-phascis flow
    for (int ip = 0; ip < nphi_p; ip++)
   {
        
        Grad_phi_i(0,0) = grad_phi_p(0,ip)*axes_p(0,0)+grad_phi_p(1,ip)*axes_p(1,0);
        Grad_phi_i(1,0) = grad_phi_p(0,ip)*axes_p(0,1)+grad_phi_p(1,ip)*axes_p(1,1);
        
        REAL dot = 0.0;
        for (int i = 0;  i < m_Dim; i++)
        {
            dot += Grad_p(i,0) * Grad_phi_i(i,0);
        }
        
        ef(ip + first_p, 0)		+=  weight * ( c * dot + (phi_poro/dt) * phip(ip,0));
        
        for (int jp = 0; jp < nphi_p; jp++)
        {
            
            Grad_phi_j(0,0) = grad_phi_p(0,jp)*axes_p(0,0)+grad_phi_p(1,jp)*axes_p(1,0);
            Grad_phi_j(1,0) = grad_phi_p(0,jp)*axes_p(0,1)+grad_phi_p(1,jp)*axes_p(1,1);
            
            REAL dot = 0.0;
            for (int i = 0;  i < m_Dim; i++)
            {
                dot += Grad_phi_j(i,0) * Grad_phi_i(i,0);
            }
            
            ek(ip + first_p, jp + first_p)		+= weight * (c * dot + (m_Se/dt) * phip(jp,0) * phip(ip,0));
        }
        
    }
    
    
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << "<<< TPZPMRSCouplPoroPlast<T,TMEM>::Contribute_2D ***";
        sout << " Resultant rhs vector:\n" << ef;
        sout << " Resultant stiff vector:\n" << ek;
        LOGPZ_DEBUG(logger,sout.str().c_str());
    }
#endif
    
    
    
    if (m_SetRunPlasticity)
    {
        ContributePlastic_2D(datavec[0],weight,ek,ef);
        
#ifdef LOG4CXX
        if(logger->isDebugEnabled())
        {
            std::stringstream sout;
            sout << "<<< TPZPMRSCouplPoroPlast<T,TMEM>::ContributePlastic_2D ***";
            sout << " Resultant rhs vector:\n" << ef;
            sout << " Resultant stiff vector:\n" << ek;
            LOGPZ_DEBUG(logger,sout.str().c_str());
        }
#endif
        
        return;
    }
    
}



/** @brief of contribute of plasticity in 2 dimensional */
template<class T,class TMEM>
void TPZPMRSCouplPoroPlast<T,TMEM>::ContributePlastic_2D(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef)
{

    TPZFMatrix<REAL> &dphi = data.dphix, dphiXY;
    TPZFMatrix<REAL> &phi  = data.phi;
    dphiXY = dphi;
    
    
    int first_u = 0;
    const int n_phi_u = phi.Rows();
    
    TPZFNMatrix<4>  Deriv(2,2);
    TPZFNMatrix<9> Dep(3,3);
    TPZFNMatrix<3>  DeltaStrain(3,1);
    TPZFNMatrix<3>  Stress(3,1);
    int ptindex = data.intGlobPtIndex;
    

        // Loop over the solutions if update memory is true
    TPZSolVec locsol(data.sol);
    TPZGradSolVec locdsol(data.dsol);
    int numsol = locsol.size();
        
    for (int is=0; is<numsol; is++)
    {
        data.sol[0] = locsol[is];
        data.dsol[0] = locdsol[is];
        
        this->ComputeDeltaStrainVector(data, DeltaStrain);
        this->ApplyDeltaStrainComputeDep(data, DeltaStrain, Stress, Dep);
    }


#ifdef MACOS
    feclearexcept(FE_ALL_EXCEPT);
    if(fetestexcept(/*FE_DIVBYZERO*/ FE_ALL_EXCEPT    )) {
        std::cout << "division by zero reported\n";
        DebugStop();
    }
#endif
    
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << ">>> TPZPMRSCouplPoroPlast<T,TMEM>::ContributePlastic_2D ***";
        sout << "\nIntegration Local Point index = " << data.intGlobPtIndex;
        sout << "\nIntegration Global Point index = " << data.intGlobPtIndex;
        sout << "\ndata.axes = " << data.axes;
        sout << "\nDep " <<endl;
        sout << Dep(0,0) << "\t" << Dep(0,1) << "\t" << Dep(0,2) <<"\n";
        sout << Dep(1,0) << "\t" << Dep(1,1) << "\t" << Dep(1,2) <<"\n";
        sout << Dep(2,0) << "\t" << Dep(2,1) << "\t" << Dep(2,2) <<"\n";
        
        sout << "\nStress " <<endl;
        sout << Stress(0,0) << "\t" << Stress(1,0) << "\t" << Stress(2,0) <<"\n";
        
        sout << "\nDELTA STRAIN " <<endl;
        sout << DeltaStrain(0,0) << "\t" << DeltaStrain(1,0) << "\t" << DeltaStrain(2,0) <<"\n";
        sout << "data.phi" << data.phi;
        
        LOGPZ_DEBUG(logger,sout.str().c_str());
    }
#endif
    
    ptindex = 0;
    REAL val;
    
    for(int iu = 0; iu < n_phi_u; iu++)
    {
        
        val -= m_b[0] * phi(0,iu);
        val += Stress(0,0) * dphiXY(0,iu); //dphixdx
        val += Stress(2,0) * dphiXY(1,iu); //dphixdy
        ef(2*iu+0 + first_u,0) += weight * val;
        
        val -= m_b[1] * phi(1,iu);
        val += Stress(2,0) * dphiXY(2,iu); //dphiydx
        val += Stress(1,0) * dphiXY(3,iu); //dphiydy
        ef(2*iu+1 + first_u,0) += weight * val;
        
        
        for (int ju = 0; ju < n_phi_u; ju++)
        {
            for (int ud = 0; ud < 2; ud++)
            {
                for (int vd = 0; vd < 2; vd++)
                {
                    Deriv(vd, ud) = dphiXY(vd, iu) * dphiXY(ud, ju);
                }
            }
        
        
        
        val  = 2. * Dep(0,0) * Deriv(0, 0); //dvdx*dudx
        val +=      Dep(0,2) * Deriv(0, 1); //dvdx*dudy
        val += 2. * Dep(2,0) * Deriv(1, 0); //dvdy*dudx
        val +=      Dep(2,2) * Deriv(1, 1); //dvdy*dudy
        val *= 0.5;
        ek(2*iu+0 + first_u, 2*ju+0 + first_u) += weight * val;
        
        val  =      Dep(0,2) * Deriv(0, 0); //dvdx*dudx
        val += 2. * Dep(0,1) * Deriv(0, 1); //dvdx*dudy
        val +=      Dep(2,2) * Deriv(1, 0); //dvdy*dudx
        val += 2. * Dep(2,1) * Deriv(1, 1); //dvdy*dudy
        val *= 0.5;
        ek(2*iu+0 + first_u, 2*ju+1 + first_u) += weight * val;
        
        val  = 2. * Dep(2,0) * Deriv(0, 0); //dvdx*dudx
        val +=      Dep(2,2) * Deriv(0, 1); //dvdx*dudy
        val += 2. * Dep(1,0) * Deriv(1, 0); //dvdy*dudx
        val +=      Dep(1,2) * Deriv(1, 1); //dvdy*dudy
        val *= 0.5;
        ek(2*iu+1 + first_u, 2*ju+0 + first_u) += weight * val;
        
        val  =      Dep(2,2) * Deriv(0, 0); //dvdx*dudx
        val += 2. * Dep(2,1) * Deriv(0, 1); //dvdx*dudy
        val +=      Dep(1,2) * Deriv(1, 0); //dvdy*dudx
        val += 2. * Dep(1,1) * Deriv(1, 1); //dvdy*dudy
        val *= 0.5;
        ek(2*iu+1 + first_u, 2*ju+1 + first_u) += weight * val;
            
        }
        
    }
    
}


// Contribute Methods being used
template<class T,class TMEM>
void TPZPMRSCouplPoroPlast<T,TMEM>::Contribute_3D(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    
    int u_b = 0;
    int p_b = 1;
    
    // Getting the space functions
    TPZFMatrix<REAL>    &phiu   =   datavec[u_b].phi;
    TPZFMatrix<REAL>    &phip   =   datavec[p_b].phi;
    
    TPZFMatrix<REAL>    &grad_phi_u   =   datavec[u_b].dphix;
    TPZFMatrix<REAL>    &grad_phi_p   =   datavec[p_b].dphix;
    
    TPZFNMatrix <9,REAL>	&axes_u	=	datavec[u_b].axes;
    TPZFNMatrix <9,REAL>	&axes_p	=	datavec[p_b].axes;
    
    // Getting the solutions and derivatives
    TPZManVector<REAL,3> u = datavec[u_b].sol[0];
    TPZManVector<REAL,1> p = datavec[p_b].sol[0];
    
    TPZFNMatrix <9,REAL> du = datavec[u_b].dsol[0];
    TPZFNMatrix <9,REAL> dp = datavec[p_b].dsol[0];
    
    TPZFNMatrix<9,REAL> Grad_p(3,1,0.0),Grad_phi_i(3,1,0.0),Grad_phi_j(3,1,0.0);
    Grad_p(0,0) = dp(0,0)*axes_p(0,0)+dp(1,0)*axes_p(1,0)+dp(2,0)*axes_p(2,0);
    Grad_p(1,0) = dp(0,0)*axes_p(0,1)+dp(1,0)*axes_p(1,1)+dp(2,0)*axes_p(2,1);
    Grad_p(2,0) = dp(0,0)*axes_p(0,2)+dp(1,0)*axes_p(1,2)+dp(2,0)*axes_p(2,2);
    
    
    int nphi_u = phiu.Rows();
    int nphi_p = phip.Rows();
    
    int first_u = 0;
    int first_p = 3*nphi_u;
    
    // Compute porosity poroelastic correction
    REAL phi_poro = porosity_corrected_3D(datavec);
    
    REAL dt = m_SimulationData->dt();
    if (!m_SimulationData->IsCurrentStateQ()) {
        
        
        // Darcy mono-phascis flow
        for (int ip = 0; ip < nphi_p; ip++) {
            
            ef(ip + first_p, 0)		+= - weight * (phi_poro/dt) * phip(ip,0);
        }
        
        return;
    }
    
    
    REAL rho_avg = (1.0-phi_poro)*m_rho_s+phi_poro*m_rho_f;
    m_b[0] = rho_avg*m_SimulationData->Gravity()[0];
    m_b[1] = rho_avg*m_SimulationData->Gravity()[1];
    m_b[2] = rho_avg*m_SimulationData->Gravity()[2];

    
    // Computing Gradient of the Solution
    TPZFNMatrix<9,REAL> Grad_u(3,3,0.0),Grad_u_n,e_e,e_p,S;
//    Grad_u(0,0) = du(0,0)*axes_u(0,0)+du(1,0)*axes_u(1,0)+du(2,0)*axes_u(2,0); // dux/dx
//    Grad_u(0,1) = du(0,0)*axes_u(0,1)+du(1,0)*axes_u(1,1)+du(2,0)*axes_u(2,1); // dux/dy
//    Grad_u(0,2) = du(0,0)*axes_u(0,2)+du(1,0)*axes_u(1,2)+du(2,0)*axes_u(2,2); // dux/dz
//    
//    Grad_u(1,0) = du(0,1)*axes_u(0,0)+du(1,1)*axes_u(1,0)+du(2,1)*axes_u(2,0); // duy/dx
//    Grad_u(1,1) = du(0,1)*axes_u(0,1)+du(1,1)*axes_u(1,1)+du(2,1)*axes_u(2,1); // duy/dy
//    Grad_u(1,2) = du(0,1)*axes_u(0,2)+du(1,1)*axes_u(1,2)+du(2,1)*axes_u(2,2); // duy/dz
//    
//    Grad_u(2,0) = du(0,2)*axes_u(0,0)+du(1,2)*axes_u(1,0)+du(2,2)*axes_u(2,0); // duz/dx
//    Grad_u(2,1) = du(0,2)*axes_u(0,1)+du(1,2)*axes_u(1,1)+du(2,2)*axes_u(2,1); // duz/dy
//    Grad_u(2,2) = du(0,2)*axes_u(0,2)+du(1,2)*axes_u(1,2)+du(2,2)*axes_u(2,2); // duz/dz
    
    // Get the solution at the integrations points
//    long global_point_index = datavec[0].intGlobPtIndex;
//    TMEM &point_memory = TPZMatWithMem<TMEM>::fMemory[global_point_index];
    
//    TMEM &point_memory = TPZMatWithMem<TMEM>::GetMemory()[global_point_index];

//    e_e = point_memory.epsilon_e_n();
//    e_p = point_memory.epsilon_p_n();
//    Grad_u_n = point_memory.grad_u_n();
//    
//    Compute_Sigma_n(Grad_u_n, Grad_u, e_e, e_p, S);
    
    TPZFNMatrix<9,REAL> Grad_vx_i(3,1,0.0);
    TPZFNMatrix<9,REAL> Grad_vy_i(3,1,0.0);
    TPZFNMatrix<9,REAL> Grad_vz_i(3,1,0.0);
    
    TPZFNMatrix<9,REAL> Grad_v(3,3,0.0);
    TPZFNMatrix<9,REAL> Grad_vx_j(3,1,0.0);
    TPZFNMatrix<9,REAL> Grad_vy_j(3,1,0.0);
    TPZFNMatrix<9,REAL> Grad_vz_j(3,1,0.0);
    
    REAL dvxdx, dvxdy, dvxdz;
    REAL dvydx, dvydy, dvydz;
    REAL dvzdx, dvzdy, dvzdz;
    
    REAL duxdx, duxdy, duxdz;
    REAL duydx, duydy, duydz;
    REAL duzdx, duzdy, duzdz;
    
    TPZFMatrix<REAL> & Sigma_0 = m_SimulationData->PreStress();
    
    Sigma_0.Zero();
    S = Sigma_0; // Applying prestress
    
    for (int iu = 0; iu < nphi_u; iu++)
    {
        
        // Computing Gradient of the test function for each component
        for (int d = 0; d < m_Dim; d++)
        {
            Grad_vx_i(d,0) = grad_phi_u(d,iu);
            Grad_vy_i(d,0) = grad_phi_u(d,iu);
            Grad_vz_i(d,0) = grad_phi_u(d,iu);
        }
        
        ef(3*iu   + first_u, 0)    += weight * ((S(0,0) - m_alpha*p[0]) * Grad_vx_i(0,0) + S(0,1) * Grad_vx_i(1,0) + S(0,2) * Grad_vx_i(2,0) - m_b[0] * phiu(iu,0));
        ef(3*iu+1 + first_u, 0)    += weight * (S(1,0) * Grad_vy_i(0,0) + (S(1,1) - m_alpha*p[0]) * Grad_vy_i(1,0) + S(1,2) * Grad_vy_i(2,0) - m_b[1] * phiu(iu,0));
        ef(3*iu+2 + first_u, 0)    += weight * (S(2,0) * Grad_vz_i(0,0) + S(2,1) * Grad_vz_i(1,0) + (S(2,2) - m_alpha*p[0]) * Grad_vz_i(2,0) - m_b[2] * phiu(iu,0));
        
        //x
        dvxdx = Grad_vx_i(0,0);
        dvxdy = Grad_vx_i(1,0);
        dvxdz = Grad_vx_i(2,0);
        
        //y
        dvydx = Grad_vy_i(0,0);
        dvydy = Grad_vy_i(1,0);
        dvydz = Grad_vy_i(2,0);
        
        //z
        dvzdx = Grad_vz_i(0,0);
        dvzdy = Grad_vz_i(1,0);
        dvzdz = Grad_vz_i(2,0);
        
        
        for (int ju = 0; ju < nphi_u; ju++)
        {
            
            // Computing Gradient of the test function for each component
            for (int d = 0; d < m_Dim; d++)
            {
                Grad_vx_j(d,0) = grad_phi_u(d,ju);
                Grad_vy_j(d,0) = grad_phi_u(d,ju);
                Grad_vz_j(d,0) = grad_phi_u(d,ju);
            }
            
            //x
            duxdx = Grad_vx_j(0,0);
            duxdy = Grad_vx_j(1,0);
            duxdz = Grad_vx_j(2,0);
            
            //y
            duydx = Grad_vy_j(0,0);
            duydy = Grad_vy_j(1,0);
            duydz = Grad_vy_j(2,0);
            
            //z
            duzdx = Grad_vz_j(0,0);
            duzdy = Grad_vz_j(1,0);
            duzdz = Grad_vz_j(2,0);
            
            // Gradient 1
            ek(3*iu   + first_u, 3*ju    + first_u) += weight * ((m_lambda + 2.*m_mu)*duxdx*dvxdx + m_mu*duxdy*dvxdy + m_mu*duxdz*dvxdz);
            ek(3*iu   + first_u, 3*ju+1  + first_u) += weight * (m_lambda*duydy*dvxdx + m_mu*duydx*dvxdy);
            ek(3*iu   + first_u, 3*ju+2  + first_u) += weight * (m_lambda*duzdz*dvxdx + m_mu*duzdx*dvxdz);
            
            // Gradient 2
            ek(3*iu+1 + first_u, 3*ju    + first_u) += weight * (m_lambda*duxdx*dvydy + m_mu*duxdy*dvydx);
            ek(3*iu+1 + first_u, 3*ju+1  + first_u) += weight * ((m_lambda + 2.*m_mu)*duydy*dvydy + m_mu*duydx*dvydx + m_mu*duydz*dvydz);
            ek(3*iu+1 + first_u, 3*ju+2  + first_u) += weight * (m_lambda*duzdz*dvydy + m_mu*duzdy*dvydz);
            
            // Gradient 3
            ek(3*iu+2 + first_u, 3*ju    + first_u) += weight * (m_lambda*duxdx*dvzdz + m_mu*duxdz*dvzdx);
            ek(3*iu+2 + first_u, 3*ju+1  + first_u) += weight * (m_lambda*duydy*dvzdz + m_mu*duydz*dvzdy);
            ek(3*iu+2 + first_u, 3*ju+2  + first_u) += weight * ((m_lambda + 2.*m_mu)*duzdz*dvzdz + m_mu*duzdx*dvzdx + m_mu*duzdy*dvzdy);
            
        }
    }
    
    TPZFNMatrix<9,REAL> dv(3,1,0.0);
    
    //    Matrix -Qc
    //    Coupling matrix
    for(int iu = 0; iu < nphi_u; iu++ )
    {
        
        // Computing Gradient of the test function for each component
        Grad_vx_i(0,0) = grad_phi_u(0,iu)*axes_u(0,0)+grad_phi_u(1,iu)*axes_u(1,0)+grad_phi_u(2,iu)*axes_u(2,0); // dvx/dx
        Grad_vx_i(1,0) = grad_phi_u(0,iu)*axes_u(0,1)+grad_phi_u(1,iu)*axes_u(1,1)+grad_phi_u(2,iu)*axes_u(2,1); // dvx/dy
        Grad_vx_i(2,0) = grad_phi_u(0,iu)*axes_u(0,2)+grad_phi_u(1,iu)*axes_u(1,2)+grad_phi_u(2,iu)*axes_u(2,2); // dvx/dz
        
        Grad_vy_i(0,0) = grad_phi_u(0,iu)*axes_u(0,0)+grad_phi_u(1,iu)*axes_u(1,0)+grad_phi_u(2,iu)*axes_u(2,0); // dvy/dx
        Grad_vy_i(1,0) = grad_phi_u(0,iu)*axes_u(0,1)+grad_phi_u(1,iu)*axes_u(1,1)+grad_phi_u(2,iu)*axes_u(2,1); // dvy/dy
        Grad_vy_i(2,0) = grad_phi_u(0,iu)*axes_u(0,2)+grad_phi_u(1,iu)*axes_u(1,2)+grad_phi_u(2,iu)*axes_u(2,2); // dvy/dz
        
        Grad_vz_i(0,0) = grad_phi_u(0,iu)*axes_u(0,0)+grad_phi_u(1,iu)*axes_u(1,0)+grad_phi_u(2,iu)*axes_u(2,0); // dvz/dx
        Grad_vz_i(1,0) = grad_phi_u(0,iu)*axes_u(0,1)+grad_phi_u(1,iu)*axes_u(1,1)+grad_phi_u(2,iu)*axes_u(2,1); // dvz/dy
        Grad_vz_i(2,0) = grad_phi_u(0,iu)*axes_u(0,2)+grad_phi_u(1,iu)*axes_u(1,2)+grad_phi_u(2,iu)*axes_u(2,2); // dvz/dz
        
        for(int jp = 0; jp < nphi_p; jp++)
        {
            
            ek(3*iu+0,first_p+jp) += (-1.)* weight * m_alpha * phip(jp,0) * Grad_vx_i(0,0);
            ek(3*iu+1,first_p+jp) += (-1.)* weight * m_alpha * phip(jp,0) * Grad_vy_i(1,0);
            ek(3*iu+2,first_p+jp) += (-1.)* weight * m_alpha * phip(jp,0) * Grad_vz_i(2,0);
        }
    }
    
    //    Matrix QcˆT
    //    Coupling matrix transpose
    for(int ip = 0; ip < nphi_p; ip++ )
    {
        
        
        for(int ju = 0; ju < nphi_u; ju++)
        {
            
            dv(0,0) = grad_phi_u(0,ju)*axes_u(0,0)+grad_phi_u(1,ju)*axes_u(1,0)+grad_phi_u(2,ju)*axes_u(2,0);
            dv(1,0) = grad_phi_u(0,ju)*axes_u(0,1)+grad_phi_u(1,ju)*axes_u(1,1)+grad_phi_u(2,ju)*axes_u(2,1);
            dv(2,0) = grad_phi_u(0,ju)*axes_u(0,2)+grad_phi_u(1,ju)*axes_u(1,2)+grad_phi_u(2,ju)*axes_u(2,2);
            
            ek(first_p+ip,3*ju+0) += (1./dt) * weight * m_alpha * dv(0,0) * phip(ip,0);
            ek(first_p+ip,3*ju+1) += (1./dt) * weight * m_alpha * dv(1,0) * phip(ip,0);
            ek(first_p+ip,3*ju+2) += (1./dt) * weight * m_alpha * dv(2,0) * phip(ip,0);
            
        }
    }
    
    
    /** @brief Rudnicki diffusion coefficient */
    /** J. W. Rudnicki. Fluid mass sources and point forces in linear elastic diffusive solids. Journal of Mechanics of Materials, 5:383–393, 1986. */
    REAL k = 0.0;
    m_k_model = 1;
    k_permeability(phi_poro,k);
    m_lambdau = 1.1 * m_lambda;
    REAL c = (k/m_eta)*(m_lambdau-m_lambda)*(m_lambda + 2.0*m_mu)/(m_alpha*m_alpha*(m_lambdau + 2.0*m_mu));
    
    // Darcy mono-phascis flow
    for (int ip = 0; ip < nphi_p; ip++)
    {
        
        Grad_phi_i(0,0) = grad_phi_p(0,ip)*axes_p(0,0)+grad_phi_p(1,ip)*axes_p(1,0)+grad_phi_p(2,ip)*axes_p(2,0);
        Grad_phi_i(1,0) = grad_phi_p(0,ip)*axes_p(0,1)+grad_phi_p(1,ip)*axes_p(1,1)+grad_phi_p(2,ip)*axes_p(2,1);
        Grad_phi_i(2,0) = grad_phi_p(0,ip)*axes_p(0,2)+grad_phi_p(1,ip)*axes_p(1,2)+grad_phi_p(2,ip)*axes_p(2,2);
        
        REAL dot = 0.0;
        for (int i = 0;  i < m_Dim; i++)
        {
            dot += Grad_p(i,0) * Grad_phi_i(i,0);
        }
        
        ef(ip + first_p, 0)		+=  weight * ( c * dot + (phi_poro/dt) * phip(ip,0));
        
        for (int jp = 0; jp < nphi_p; jp++)
        {
            
            Grad_phi_j(0,0) = grad_phi_p(0,jp)*axes_p(0,0)+grad_phi_p(1,jp)*axes_p(1,0)+grad_phi_p(2,jp)*axes_p(2,0);
            Grad_phi_j(1,0) = grad_phi_p(0,jp)*axes_p(0,1)+grad_phi_p(1,jp)*axes_p(1,1)+grad_phi_p(2,jp)*axes_p(2,1);
            Grad_phi_j(2,0) = grad_phi_p(0,jp)*axes_p(0,2)+grad_phi_p(1,jp)*axes_p(1,2)+grad_phi_p(2,jp)*axes_p(2,2);
            
            REAL dot = 0.0;
            for (int i = 0;  i < m_Dim; i++)
            {
                dot += Grad_phi_j(i,0) * Grad_phi_i(i,0);
            }
            
            ek(ip + first_p, jp + first_p)		+= weight * (c * dot + (m_Se/dt) * phip(jp,0) * phip(ip,0) );
        }
        
    }
    
}



/** @brief of contribute  */
template<class T,class TMEM>
void TPZPMRSCouplPoroPlast<T,TMEM>::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef){
    
    TPZFMatrix<STATE>  ek_fake(ef.Rows(),ef.Rows(),0.0);
    this->Contribute(datavec, weight, ek_fake, ef);
    
}

/** @brief of contribute of BC */
template<class T,class TMEM>
void TPZPMRSCouplPoroPlast<T,TMEM>::ContributeBC(TPZVec<TPZMaterialData> &datavec,REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{
    
    if (!m_SimulationData->IsCurrentStateQ())
    {
        return;
    }
    
    if (m_Dim == 3)
    {
        this->ContributeBC_3D(datavec, weight, ek, ef, bc);
    }
    else
    {
         this->ContributeBC_2D(datavec, weight, ek, ef, bc);
    }
    
    
}


/** @brief of contribute of BC_2D */
template<class T,class TMEM>
void TPZPMRSCouplPoroPlast<T,TMEM>::ContributeBC_2D(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{

    int u_b = 0;
    int p_b = 1;

    TPZFMatrix<REAL>  &phiu = datavec[u_b].phi;
    TPZFMatrix<REAL>  &phip = datavec[p_b].phi;

    // Getting the solutions and derivatives
    TPZManVector<REAL,2> u = datavec[u_b].sol[0];
    TPZManVector<REAL,1> p = datavec[p_b].sol[0];

    int phru = phiu.Rows();
    int phrp = phip.Rows();
    short in,jn;
    const REAL BIGNUMBER = TPZMaterial::gBigNumber;

    
    // The boundary is time dependent
    REAL v[3];
    REAL time = this->SimulationData()->t();
    
    if (bc.HasTimedependentBCForcingFunction())
    {
        TPZManVector<REAL,3> f(3);
        TPZFMatrix<REAL> gradf;
        bc.TimedependentBCForcingFunction()->Execute(datavec[u_b].x, time, f, gradf);
        v[0] = f[0];    //  Ux displacement or Tnx
        v[1] = f[1];    //  Uy displacement or Tny
        v[2] = f[2];    //  Pressure or Qn
    }
    
    
    // Boundaries

    // Dirichlet in Pressure
    switch (bc.Type())
    {
        case 0 : // Du_Dp
        {
            REAL v[3];
            v[0] = bc.Val2()(0,0);    //    Ux displacement
            v[1] = bc.Val2()(1,0);    //    Uy displacement
            v[2] = bc.Val2()(2,0);    //    Pressure

            for(in = 0 ; in < phru; in++)
            {
                //    Contribution for load Vector
                ef(2*in  ,0)      += BIGNUMBER*(u[0] - v[0])*phiu(in,0)*weight;    // X displacement Value
                ef(2*in+1,0)      += BIGNUMBER*(u[1] - v[1])*phiu(in,0)*weight;    // y displacement Value

                for (jn = 0 ; jn < phru; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(2*in,2*jn    )    += BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;    // X displacement
                    ek(2*in+1,2*jn+1)    += BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;    // Y displacement
                }
            }

            //    Diffusion Equation
            REAL p_s = p[0];
            REAL d_p = (p_s-v[2]);
            for(in = 0 ; in < phrp; in++)
            {
                //    Contribution for load Vector
                ef(in+2*phru,0)        += BIGNUMBER*(d_p)*phip(in,0)*weight;    // P Pressure

                for (jn = 0 ; jn < phrp; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(in+2*phru,jn+2*phru)        += BIGNUMBER*phip(in,0)*phip(jn,0)*weight;    // P Pressure
                }
            }
            break;
        }

        case 1 : // Dux_Dp
        {
            REAL v[2];
            v[0] = bc.Val2()(0,0);    //    Ux displacement
            v[1] = bc.Val2()(1,0);    //    Pressure

            //    Diffusion Equation
            REAL ux_s = u[0];
            REAL d_ux = (ux_s-v[0]);
            for(in = 0 ; in < phru; in++)
            {
                //    Contribution for load Vector
                ef(2*in,0)        += BIGNUMBER*(d_ux)*phiu(in,0)*weight;    // X displacement Value

                for (jn = 0 ; jn < phru; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(2*in,2*jn)        += BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;    // X displacement
                }
            }

            //    Diffusion Equation
            REAL p_s = p[0];
            REAL d_p = (p_s-v[1]);
            for(in = 0 ; in < phrp; in++)
            {
                //    Contribution for load Vector
                ef(in+2*phru,0)        += BIGNUMBER*(d_p)*phip(in,0)*weight;    // P Pressure

                for (jn = 0 ; jn < phrp; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(in+2*phru,jn+2*phru)        += BIGNUMBER*phip(in,0)*phip(jn,0)*weight;    // P Pressure
                }
            }
            break;
        }

        case 2 : //Duy_Dp
        {

            REAL v[2];
            v[0] = bc.Val2()(0,0);    //    Uy displacement
            v[1] = bc.Val2()(1,0);    //    Pressure

            REAL uy_s = u[1];
            REAL d_uy = (uy_s-v[0]);
            for(in = 0 ; in < phru; in++)
            {
                //    Contribution for load Vector
                ef(2*in+1,0)      += BIGNUMBER*(d_uy)*phiu(in,0)*weight;    // y displacement

                for (jn = 0 ; jn < phru; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(2*in+1,2*jn+1)    += BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;    // Y displacement
                }
            }

            //    Diffusion Equation
            REAL p_s = p[0];
            REAL d_p = (p_s-v[1]);
            for(in = 0 ; in < phrp; in++)
            {
                //    Contribution for load Vector
                ef(in+2*phru,0)        += BIGNUMBER*(d_p)*phip(in,0)*weight;    // P Pressure

                for (jn = 0 ; jn < phrp; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(in+2*phru,jn+2*phru)        += BIGNUMBER*phip(in,0)*phip(jn,0)*weight;    // P Pressure
                }
            }
            break;
        }

        case 3 : // Nt_Dp
        {

            REAL v[3];
            v[0] = bc.Val2()(0,0);    //    Tnx
            v[1] = bc.Val2()(1,0);    //    Tny
            v[2] = bc.Val2()(2,0);    //    Pressure

            //    Neumann condition for each state variable
            //    Elasticity Equation
            for(in = 0 ; in <phru; in++)
            {
                //    Normal Tension Components on neumman boundary
                ef(2*in  ,0)    += -1.0 * weight * v[0] * phiu(in,0);        //    Tnx
                ef(2*in+1,0)    += -1.0 * weight * v[1] * phiu(in,0);        //    Tny
            }

            //    Diffusion Equation
            REAL p_s = p[0];
            REAL d_p = (p_s-v[2]);
            for(in = 0 ; in < phrp; in++)
            {
                //    Contribution for load Vector
                ef(in+2*phru,0)        += BIGNUMBER*(d_p)*phip(in,0)*weight;    // P Pressure

                for (jn = 0 ; jn < phrp; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(in+2*phru,jn+2*phru)        += BIGNUMBER*phip(in,0)*phip(jn,0)*weight;    // P Pressure
                }
            }
            break;
        }

        case 4 : // Ntn_Dp
        {

            REAL v[2];
            v[0] = bc.Val2()(0,0);    //    Tn normal traction
            v[1] = bc.Val2()(1,0);    //    Pressure

            REAL tn = v[0];
            TPZManVector<REAL,2> n = datavec[u_b].normal;
            //    Neumann condition for each state variable
            //    Elasticity Equation
            for(in = 0 ; in <phru; in++)
            {
                //    Normal Tension Components on neumman boundary
                ef(2*in  ,0)    += -1.0 * weight * tn * n[0] * phiu(in,0);        //    Tnx
                ef(2*in+1,0)    += -1.0 * weight * tn * n[1] * phiu(in,0);        //    Tny
            }

            //    Diffusion Equation
            REAL p_s = p[0];
            REAL d_p = (p_s-v[1]);
            for(in = 0 ; in < phrp; in++)
            {
                //    Contribution for load Vector
                ef(in+2*phru,0)        += BIGNUMBER*(d_p)*phip(in,0)*weight;    // P Pressure

                for (jn = 0 ; jn < phrp; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(in+2*phru,jn+2*phru)        += BIGNUMBER*phip(in,0)*phip(jn,0)*weight;    // P Pressure
                }
            }
            break;
        }

        case 5 : // Du_Nq
        {

            REAL v[3];
            v[0] = bc.Val2()(0,0);    //    Ux displacement
            v[1] = bc.Val2()(1,0);    //    Uy displacement
            v[2] = bc.Val2()(2,0);    //    Qn

            for(in = 0 ; in < phru; in++)
            {
                //    Contribution for load Vector
                ef(2*in,0)        += BIGNUMBER*(u[0] - v[0])*phiu(in,0)*weight;    // X displacement Value
                ef(2*in+1,0)      += BIGNUMBER*(u[1] - v[1])*phiu(in,0)*weight;    // y displacement Value

                for (jn = 0 ; jn < phru; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(2*in,2*jn)        += BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;    // X displacement
                    ek(2*in+1,2*jn+1)    += BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;    // Y displacement
                }
            }

            //    Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //    Normal Flux on neumman boundary
                ef(in+2*phru,0)    += 1.0 *  weight * v[2] * phip(in,0);    // Qnormal
            }
            break;
        }

        case 6 : // Dux_Nq
        {


            REAL v[2];
            v[0] = bc.Val2()(0,0);    //    Ux displacement
            v[1] = bc.Val2()(1,0);    //    Qn

            REAL ux_s = u[0];
            REAL d_ux = (ux_s-v[0]);
            for(in = 0 ; in < phru; in++)
            {
                //    Contribution for load Vector
                ef(2*in,0)        += BIGNUMBER*(d_ux)*phiu(in,0)*weight;    // X displacement Value

                for (jn = 0 ; jn < phru; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(2*in,2*jn)        += BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;    // X displacement
                }
            }

            //    Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //    Normal Flux on neumman boundary
                ef(in+2*phru,0)    += 1.0 * weight * v[1] * phip(in,0);    // Qnormal
            }
            break;
        }

        case 7 : // Duy_Nq
        {

            REAL v[2];
            v[0] = bc.Val2()(0,0);    //    Uy displacement
            v[1] = bc.Val2()(1,0);    //    Qn

            REAL uy_s = u[1];
            REAL d_uy = (uy_s-v[0]);
            for(in = 0 ; in < phru; in++)
            {
                //    Contribution for load Vector
                ef(2*in+1,0)      += BIGNUMBER*(d_uy)*phiu(in,0)*weight;    // y displacement Value

                for (jn = 0 ; jn < phru; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(2*in+1,2*jn+1)    += BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;    // Y displacement
                }
            }

            //    Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //    Normal Flux on neumman boundary
                ef(in+2*phru,0)    += 1.0 * weight * v[1] * phip(in,0);    // Qnormal
            }
            break;
        }

        case 8 : // Nt_Nq
        {

            REAL v[3];
            v[0] = bc.Val2()(0,0);    //    Tnx
            v[1] = bc.Val2()(1,0);    //    Tny
            v[2] = bc.Val2()(2,0);    //    Qn

            //    Neumann condition for each state variable
            //    Elasticity Equation
            for(in = 0 ; in <phru; in++)
            {
                //    Normal Tension Components on neumman boundary
                ef(2*in  ,0)     += -1.0 * weight * v[0] * phiu(in,0);        //    Tnx
                ef(2*in+1,0)     += -1.0 * weight * v[1] * phiu(in,0);        //    Tny
            }

            //    Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //    Normal Flux on neumman boundary
                ef(in+2*phru,0)    += 1.0 * weight * v[2] * phip(in,0);    // Qnormal
            }
            break;
        }

        case 9 : // Ntn_Nq
        {

            REAL v[2];
            v[0] = bc.Val2()(0,0);    //    Tn normal traction
            v[1] = bc.Val2()(2,0);    //    Qn

            REAL tn = v[0];
            TPZManVector<REAL,2> n = datavec[u_b].normal;
            //    Neumann condition for each state variable
            //    Elasticity Equation
            for(in = 0 ; in <phru; in++)
            {
                //    Normal Tension Components on neumman boundary
                ef(2*in  ,0)      += -1.0 * weight * tn * n[0] * phiu(in,0);        //    Tnx
                ef(2*in+1,0)      += -1.0 * weight * tn * n[1] * phiu(in,0);        //    Tny
            }

            //    Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //    Normal Flux on neumman boundary
                ef(in+2*phru,0)    += 1.0 * weight * v[1] * phip(in,0);    // Qnormal
            }
            break;
        }

        case 10 : //Du_time_Dp
        {
            
            //    Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //    Contribution for load Vector
                ef(2*in  ,0)      += BIGNUMBER*(u[0] - v[0])*phiu(in,0)*weight;    // X displacement Value
                ef(2*in+1,0)      += BIGNUMBER*(u[1] - v[1])*phiu(in,0)*weight;    // y displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(2*in,2*jn    )    += BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;    // X displacement
                    ek(2*in+1,2*jn+1)    += BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;    // Y displacement
                }
            }
            
            //    Diffusion Equation
            REAL p_s = p[0];
            REAL d_p = (p_s-v[2]);
            for(in = 0 ; in < phrp; in++)
            {
                //    Contribution for load Vector
                ef(in+2*phru,0)        += BIGNUMBER*(d_p)*phip(in,0)*weight;    // P Pressure
                
                for (jn = 0 ; jn < phrp; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(in+2*phru,jn+2*phru)        += BIGNUMBER*phip(in,0)*phip(jn,0)*weight;    // P Pressure
                }
            }
            break;
        }
            
        case 11 : // Ntn_time_Dp
        {
            
            //    Neumann condition for each state variable
            //    Elasticity Equation
            for(in = 0 ; in <phru; in++)
            {
                //    Normal Tension Components on neumman boundary
                ef(2*in  ,0)    += -1.0 * weight * v[0] * phiu(in,0);        //    Tnx
                ef(2*in+1,0)    += -1.0 * weight * v[1] * phiu(in,0);        //    Tny
            }
            
            //    Diffusion Equation
            REAL p_s = p[0];
            REAL d_p = (p_s-v[2]);
            for(in = 0 ; in < phrp; in++)
            {
                //    Contribution for load Vector
                ef(in+2*phru,0)        += BIGNUMBER*(d_p)*phip(in,0)*weight;    // P Pressure
                
                for (jn = 0 ; jn < phrp; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(in+2*phru,jn+2*phru)        += BIGNUMBER*phip(in,0)*phip(jn,0)*weight;    // P Pressure
                }
            }
            break;
        }
            

        default:
        {
            DebugStop();
        }
            break;
    }

}

template<class T,class TMEM>
void TPZPMRSCouplPoroPlast<T,TMEM>::ContributeBC_3D(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    int u_b = 0;
    int p_b = 1;
    
    // Getting the test (trial) function
    TPZFMatrix<REAL>  & phiu = datavec[u_b].phi;
    TPZFMatrix<REAL>  & phip = datavec[p_b].phi;
   
    // Getting the solutions (base function) and derivatives
    TPZManVector<REAL,3> u = datavec[u_b].sol[0];
    TPZManVector<REAL,1> p = datavec[p_b].sol[0];
    
    // Getting the size of test function
    int phru = phiu.Rows();
    int phrp = phip.Rows();
    short in,jn;
    const REAL BIGNUMBER = TPZMaterial::gBigNumber;
    
    switch (bc.Type())
    {
            
        case 0 : // Du_Dp
            // Dirichlet in displacement and Pressure
        {
            REAL v[4];
            v[0] = bc.Val2()(0,0);    //    Ux displacement
            v[1] = bc.Val2()(1,0);    //    Uy displacement
            v[2] = bc.Val2()(2,0);	  //	Uz displacement
            v[3] = bc.Val2()(3,0);    //    Pressure
            
            
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                ef(3*in  ,0)	+= BIGNUMBER*(u[0] - v[0])*phiu(in,0)*weight;	// X displacement Value
                ef(3*in+1,0)	+= BIGNUMBER*(u[1] - v[1])*phiu(in,0)*weight;	// y displacement Value
                ef(3*in+2,0)	+= BIGNUMBER*(u[2] - v[2])*phiu(in,0)*weight;	// Z displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in  ,3*jn  )	+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
                    ek(3*in+1,3*jn+1)	+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// Y displacement
                    ek(3*in+2,3*jn+2)	+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// Z displacement
                }
            }
            
            //    Diffusion Equation
            REAL p_s = p[0];
            REAL d_p = (p_s-v[3]);
            for(in = 0 ; in < phrp; in++)
            {
                //    Contribution for load Vector
                ef(in+3*phru,0)        += BIGNUMBER*(d_p)*phip(in,0)*weight;    // P Pressure
                
                for (jn = 0 ; jn < phrp; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(in+3*phru,jn+3*phru)        += BIGNUMBER*phip(in,0)*phip(jn,0)*weight;    // P Pressure
                }
            }
            break;
            
        }
            
            
        case 1 : // Duxy_Dp
            // Dirichlet in x and y direction of displacement and Pressure
        {
            REAL v[3];
            v[0] = bc.Val2()(0,0);    //    Ux displacement
            v[1] = bc.Val2()(1,0);    //    Uy displacement
            v[2] = bc.Val2()(2,0);    //    Pressure
            
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                ef(3*in  ,0)	+= BIGNUMBER*(u[0] - v[0])*phiu(in,0)*weight;	// X displacement Value
                ef(3*in+1,0)	+= BIGNUMBER*(u[1] - v[1])*phiu(in,0)*weight;	// y displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in  ,3*jn  )	+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
                    ek(3*in+1,3*jn+1)	+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// Y displacement
                }
            }
            //    Diffusion Equation
            REAL p_s = p[0];
            REAL d_p = (p_s-v[2]);
            for(in = 0 ; in < phrp; in++)
            {
                //    Contribution for load Vector
                ef(in+3*phru,0)        += BIGNUMBER*(d_p)*phip(in,0)*weight;    // P Pressure
                
                for (jn = 0 ; jn < phrp; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(in+3*phru,jn+3*phru)        += BIGNUMBER*phip(in,0)*phip(jn,0)*weight;    // P Pressure
                }
            }
            break;
            
        }
            
            
        case 2 : // Duxz_Dp
            // Dirichlet in x and z direction of displacement and Pressure
        {
            REAL v[3];
            v[0] = bc.Val2()(0,0);    //    Ux displacement
            v[1] = bc.Val2()(1,0);    //    Uz displacement
            v[2] = bc.Val2()(2,0);    //    Pressure
            
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                ef(3*in  ,0)	+= BIGNUMBER*(u[0] - v[0])*phiu(in,0)*weight;	// X displacement Value
                ef(3*in+2,0)	+= BIGNUMBER*(u[2] - v[1])*phiu(in,0)*weight;	// Z displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in  ,3*jn  )	+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
                    ek(3*in+2,3*jn+2)	+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// Z displacement
                }
            }
            //    Diffusion Equation
            REAL p_s = p[0];
            REAL d_p = (p_s-v[2]);
            for(in = 0 ; in < phrp; in++)
            {
                //    Contribution for load Vector
                ef(in+3*phru,0)        += BIGNUMBER*(d_p)*phip(in,0)*weight;    // P Pressure
                
                for (jn = 0 ; jn < phrp; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(in+3*phru,jn+3*phru)        += BIGNUMBER*phip(in,0)*phip(jn,0)*weight;    // P Pressure
                }
            }
            break;
            
        }
            
            
        case 3 : // Duyz_Dp
            // Dirichlet in y and z direction of displacement and Pressure
        {
            REAL v[3];
            v[0] = bc.Val2()(0,0);    //    Uy displacement
            v[1] = bc.Val2()(1,0);    //    Uz displacement
            v[2] = bc.Val2()(2,0);    //    Pressure
            
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                ef(3*in+1,0)	+= BIGNUMBER*(u[1] - v[0])*phiu(in,0)*weight;	// Y displacement Value
                ef(3*in+2,0)	+= BIGNUMBER*(u[2] - v[1])*phiu(in,0)*weight;	// Z displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in+1,3*jn+1)	+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// Y displacement
                    ek(3*in+2,3*jn+2)	+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// Z displacement
                }
            }
            //    Diffusion Equation
            REAL p_s = p[0];
            REAL d_p = (p_s-v[2]);
            for(in = 0 ; in < phrp; in++)
            {
                //    Contribution for load Vector
                ef(in+3*phru,0)        += BIGNUMBER*(d_p)*phip(in,0)*weight;    // P Pressure
                
                for (jn = 0 ; jn < phrp; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(in+3*phru,jn+3*phru)        += BIGNUMBER*phip(in,0)*phip(jn,0)*weight;    // P Pressure
                }
            }
            break;
            
        }
            
            
        case 4 : // Dux_Dp
            // Dirichlet in x direction of displacement and Pressure
        {
            REAL v[2];
            v[0] = bc.Val2()(0,0);    //    Ux displacement
            v[1] = bc.Val2()(1,0);    //    Pressure
            
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                ef(3*in,0)		+= BIGNUMBER*(u[0] - v[0])*phiu(in,0)*weight;	// X displacement Value
                
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in,3*jn)		+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
                    
                }
            }
            //    Diffusion Equation
            REAL p_s = p[0];
            REAL d_p = (p_s-v[1]);
            for(in = 0 ; in < phrp; in++)
            {
                //    Contribution for load Vector
                ef(in+3*phru,0)        += BIGNUMBER*(d_p)*phip(in,0)*weight;    // P Pressure
                
                for (jn = 0 ; jn < phrp; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(in+3*phru,jn+3*phru)        += BIGNUMBER*phip(in,0)*phip(jn,0)*weight;    // P Pressure
                }
            }
            break;
            
        }
            
            
        case 5 : // Duy_Dp
            // Dirichlet in y direction of displacement and Pressure
        {
            REAL v[2];
            v[0] = bc.Val2()(0,0);    //    Uy displacement
            v[1] = bc.Val2()(1,0);    //    Pressure
            
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                ef(3*in+1,0)		+= BIGNUMBER*(u[1] - v[0])*phiu(in,0)*weight;	// Y displacement Value
                
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in+1,3*jn+1)		+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// Y displacement
                    
                }
            }
            //    Diffusion Equation
            REAL p_s = p[0];
            REAL d_p = (p_s-v[1]);
            for(in = 0 ; in < phrp; in++)
            {
                //    Contribution for load Vector
                ef(in+3*phru,0)        += BIGNUMBER*(d_p)*phip(in,0)*weight;    // P Pressure
                
                for (jn = 0 ; jn < phrp; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(in+3*phru,jn+3*phru)        += BIGNUMBER*phip(in,0)*phip(jn,0)*weight;    // P Pressure
                }
            }
            break;
            
        }
            
            
        case 6 : // Duz_Dp
            // Dirichlet in z direction of displacement and Pressure
        {
            REAL v[2];
            v[0] = bc.Val2()(0,0);    //    Uz displacement
            v[1] = bc.Val2()(1,0);    //    Pressure
            
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                ef(3*in+2,0)		+= BIGNUMBER*(u[2] - v[0])*phiu(in,0)*weight;	// Z displacement Value
                
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in+2,3*jn+2)		+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// Z displacement
                    
                }
            }
            //    Diffusion Equation
            REAL p_s = p[0];
            REAL d_p = (p_s-v[1]);
            for(in = 0 ; in < phrp; in++)
            {
                //    Contribution for load Vector
                ef(in+3*phru,0)        += BIGNUMBER*(d_p)*phip(in,0)*weight;    // P Pressure
                
                for (jn = 0 ; jn < phrp; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(in+3*phru,jn+3*phru)        += BIGNUMBER*phip(in,0)*phip(jn,0)*weight;    // P Pressure
                }
            }
            break;
            
        }
            
            
        case 7 : //Nt_Dp
            // Neumann condition in displacement and Dirichlet in Pressure
        {
            REAL v[4];
            v[0] = bc.Val2()(0,0);    //    Tnx
            v[1] = bc.Val2()(1,0);    //    Tny
            v[2] = bc.Val2()(2,0);	  //	Tnz
            v[3] = bc.Val2()(3,0);    //    Pressure
            
            
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Normal Tension Components on neumman boundary
                ef(3*in  ,0)	+= -1.0 * weight * v[0] * phiu(in,0);		//	Tnx
                ef(3*in+1,0)	+= -1.0 * weight * v[1] * phiu(in,0);		//	Tny
                ef(3*in+2,0)	+= -1.0 * weight * v[2] * phiu(in,0);		//	Tnz
            }
            
            //    Diffusion Equation
            REAL p_s = p[0];
            REAL d_p = (p_s-v[3]);
            for(in = 0 ; in < phrp; in++)
            {
                //    Contribution for load Vector
                ef(in+3*phru,0)        += BIGNUMBER*(d_p)*phip(in,0)*weight;    // P Pressure
                
                for (jn = 0 ; jn < phrp; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(in+3*phru,jn+3*phru)        += BIGNUMBER*phip(in,0)*phip(jn,0)*weight;    // P Pressure
                }
            }
            break;
        }
            
            
        case 8 : //Ntn_Dp
            // Neumann condition in displacement and Dirichlet in Pressure
        {
            REAL v[2];
            v[0] = bc.Val2()(0,0);    //    Tn normal traction
            v[1] = bc.Val2()(1,0);    //    Pressure
            
            REAL tn = v[0];
            TPZManVector<REAL,3> n = datavec[u_b].normal;
            
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Normal Tension Components on neumman boundary
                ef(3*in  ,0)	+= -1.0 * weight * tn * n[0] * phiu(in,0);		//	Tnx
                ef(3*in+1,0)	+= -1.0 * weight * tn * n[1] * phiu(in,0);		//	Tny
                ef(3*in+2,0)	+= -1.0 * weight * tn * n[2] * phiu(in,0);		//	Tnz
            }
            
            //    Diffusion Equation
            REAL p_s = p[0];
            REAL d_p = (p_s-v[1]);
            for(in = 0 ; in < phrp; in++)
            {
                //    Contribution for load Vector
                ef(in+3*phru,0)        += BIGNUMBER*(d_p)*phip(in,0)*weight;    // P Pressure
                
                for (jn = 0 ; jn < phrp; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(in+3*phru,jn+3*phru)        += BIGNUMBER*phip(in,0)*phip(jn,0)*weight;    // P Pressure
                }
            }
            break;
        }
            
            
        case 9 : // Du_Nq
            // Dirichlet in displacement and Neumann in Pressure
        {
            REAL v[4];
            v[0] = bc.Val2()(0,0);    //    Ux displacement
            v[1] = bc.Val2()(1,0);    //    Uy displacement
            v[2] = bc.Val2()(2,0);	  //	Uz displacement
            v[3] = bc.Val2()(3,0);    //    Qn
            
            
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                ef(3*in  ,0)	+= BIGNUMBER*(u[0] - v[0])*phiu(in,0)*weight;	// X displacement Value
                ef(3*in+1,0)	+= BIGNUMBER*(u[1] - v[1])*phiu(in,0)*weight;	// Y displacement Value
                ef(3*in+2,0)	+= BIGNUMBER*(u[2] - v[2])*phiu(in,0)*weight;	// Z displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in  ,3*jn  )	+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
                    ek(3*in+1,3*jn+1)	+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// Y displacement
                    ek(3*in+2,3*jn+2)	+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// Z displacement
                }
            }
            
            //    Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //    Normal Flux on neumman boundary
                ef(in+3*phru,0)    += 1.0 * weight * v[3]*phip(in,0);    // Qnormal
            }
            break;
            
        }
            
        case 10 : // Duxy_Nq
            // Dirichlet in x and y direction of displacement and Neumann in Pressure
        {
            REAL v[3];
            v[0] = bc.Val2()(0,0);    //    Ux displacement
            v[1] = bc.Val2()(1,0);    //    Uy displacement
            v[2] = bc.Val2()(2,0);    //    Qn
            
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                ef(3*in  ,0)	+= BIGNUMBER*(u[0] - v[0])*phiu(in,0)*weight;	// X displacement Value
                ef(3*in+1,0)	+= BIGNUMBER*(u[1] - v[1])*phiu(in,0)*weight;	// Y displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in  ,3*jn  )	+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
                    ek(3*in+1,3*jn+1)	+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// Y displacement
                }
            }
            //    Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //    Normal Flux on neumman boundary
                ef(in+3*phru,0)    += weight * v[2] * phip(in,0);    // Qnormal
            }
            break;
            
        }
            
        case 11 : // Duxz_Nq
            // Dirichlet in x and z direction of displacement and Neumann in Pressure
        {
            REAL v[3];
            v[0] = bc.Val2()(0,0);    //    Ux displacement
            v[1] = bc.Val2()(1,0);    //    Uz displacement
            v[2] = bc.Val2()(2,0);    //    Qn
            
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                ef(3*in  ,0)	+= BIGNUMBER*(u[0] - v[0])*phiu(in,0)*weight;	// X displacement Value
                ef(3*in+2,0)	+= BIGNUMBER*(u[2] - v[1])*phiu(in,0)*weight;	// Z displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in  ,3*jn  )	+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
                    ek(3*in+2,3*jn+2)	+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// Z displacement
                }
            }
            //    Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //    Normal Flux on neumman boundary
                ef(in+3*phru,0)    += 1.0 * weight * v[2] * phip(in,0);    // Qnormal
            }
            break;
            
        }
            
        case 12 : // Duyz_Nq
            // Dirichlet in y and z direction of displacement and Neumann in Pressure
        {
            REAL v[3];
            v[0] = bc.Val2()(0,0);    //    Uy displacement
            v[1] = bc.Val2()(1,0);    //    Uz displacement
            v[2] = bc.Val2()(2,0);    //    Qn
            
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                ef(3*in+1,0)	+= BIGNUMBER*(u[1] - v[0])*phiu(in,0)*weight;	// Y displacement Value
                ef(3*in+2,0)	+= BIGNUMBER*(u[2] - v[1])*phiu(in,0)*weight;	// Z displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in+1,3*jn+1)	+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// Y displacement
                    ek(3*in+2,3*jn+2)	+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// Z displacement
                }
            }
            //    Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //    Normal Flux on neumman boundary
                ef(in+3*phru,0)    += 1.0 * weight * v[2] * phip(in,0);    // Qnormal
            }
            break;
            
        }
            
        case 13 : // Dux_Nq
            // Dirichlet in x direction of displacement and Neumann in Pressure
        {
            REAL v[2];
            v[0] = bc.Val2()(0,0);    //    Ux displacement
            v[1] = bc.Val2()(1,0);    //    Qn
            
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                ef(3*in,0)		+= BIGNUMBER*(u[0] - v[0])*phiu(in,0)*weight;	// X displacement Value
                
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in,3*jn)		+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
                    
                }
            }
            //    Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //    Normal Flux on neumman boundary
                ef(in+3*phru,0)    += 1.0 * weight * v[1] * phip(in,0);    // Qnormal
            }
            break;
            
        }
            
            
        case 14 : // Duy_Nq
            // Dirichlet in y direction of displacement and Neumann in Pressure
        {
            REAL v[2];
            v[0] = bc.Val2()(0,0);    //    Uy displacement
            v[1] = bc.Val2()(1,0);    //    Qn
            
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                ef(3*in+1,0)		+= BIGNUMBER*(u[1] - v[0])*phiu(in,0)*weight;	// Y displacement Value
                
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in+1,3*jn+1)		+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// Y displacement
                    
                }
            }
            //    Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //    Normal Flux on neumman boundary
                ef(in+3*phru,0)    += 1.0 * weight * v[1] * phip(in,0);    // Qnormal
            }
            break;
            
        }
            
            
        case 15 : // Duz_Nq
            // Dirichlet in z direction of displacement and Neumann in Pressure
        {
            REAL v[2];
            v[0] = bc.Val2()(0,0);    //    Uz displacement
            v[1] = bc.Val2()(1,0);    //    Qn
            
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                ef(3*in+2,0)		+= BIGNUMBER*(u[2] - v[0])*phiu(in,0)*weight;	// Z displacement Value
                
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in+2,3*jn+2)		+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// Z displacement
                    
                }
            }
            //    Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //    Normal Flux on neumman boundary
                ef(in+3*phru,0)    += 1.0 * weight * v[1] * phip(in,0);    // Qnormal
            }
            break;
            
        }
            
            
        case 16 : // Nt_Nq
            // Neumann condition in displacement and Pressure
        {
            REAL v[4];
            v[0] = bc.Val2()(0,0);    //    Tnx
            v[1] = bc.Val2()(1,0);    //    Tny
            v[2] = bc.Val2()(2,0);	  //	Tnz
            v[3] = bc.Val2()(3,0);    //    Qn
            
            
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Normal Tension Components on neumman boundary
                ef(3*in  ,0)	+= -1.0 * weight * v[0] * phiu(in,0);		//	Tnx
                ef(3*in+1,0)	+= -1.0 * weight * v[1] * phiu(in,0);		//	Tny
                ef(3*in+2,0)	+= -1.0 * weight * v[2] * phiu(in,0);		//	Tnz
            }
            
            //    Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //    Normal Flux on neumman boundary
                ef(in+3*phru,0)    += 1.0 * weight * v[3] * phip(in,0);    // Qnormal
            }
            break;
        }
            
            
        case 17 : // Ntn_Nq
            // Neumann condition in displacement and Pressure
        {
            REAL v[2];
            v[0] = bc.Val2()(0,0);    //    Tn normal traction
            v[1] = bc.Val2()(3,0);    //    Qn
            
            REAL tn = v[0];
            TPZManVector<REAL,3> n = datavec[u_b].normal;
            
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Normal Tension Components on neumman boundary
                ef(3*in  ,0)	+= -1.0 * weight * tn * n[0] * phiu(in,0);		//	Tnx
                ef(3*in+1,0)	+= -1.0 * weight * tn * n[1] * phiu(in,0);		//	Tny
                ef(3*in+2,0)	+= -1.0 * weight * tn * n[2] * phiu(in,0);		//	Tnz
            }
            
            //    Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //    Normal Flux on neumman boundary
                ef(in+3*phru,0)    += 1.0 * weight * v[1] * phip(in,0);    // Qnormal
            }
            break;
        }
            
            
        default:
        {
            DebugStop();
        }
            break;
    }
}

template<class T,class TMEM>
void TPZPMRSCouplPoroPlast<T,TMEM>::FillDataRequirements(TPZVec<TPZMaterialData > &datavec)

{
    int nref = datavec.size();
    for(int i = 0; i<nref; i++)
    {
        datavec[i].SetAllRequirements(false);
        datavec[i].fNeedsSol = true;
        datavec[i].fNeedsNeighborSol = true;
        datavec[i].fNeedsNeighborCenter = false;
        datavec[i].fNeedsNormal = true;
    }
}

template<class T,class TMEM>
void TPZPMRSCouplPoroPlast<T,TMEM>::FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec)
{
    int nref = datavec.size();
    for(int i = 0; i<nref; i++)
    {
        datavec[i].fNeedsSol = true;
        datavec[i].fNeedsNormal = true;
        datavec[i].fNeedsNeighborSol = true;
    }
}

template<class T,class TMEM>
void TPZPMRSCouplPoroPlast<T,TMEM>::Print(std::ostream &out)
{
    out << "Material Name : " << Name() << "\n";
    out << "Properties for TPZPMRSCouplPoroPlast: \n";
    out << "\t Poisson Ratio   = "											<< m_nu		<< std::endl;
    out << "\t Undarined Poisson Ratio   = "								<< m_nuu		<< std::endl;
    out << "\t First Lamé Parameter   = "									<< m_lambda	<< std::endl;
    out << "\t Second Lamé Parameter   = "									<< m_mu		<< std::endl;
    out << "\t Undrained First Lamé Parameter   = "							<< m_lambdau	<< std::endl;
    out << "\t Biot coefficient   = "										<< m_alpha	<< std::endl;
    out << "\t Body force vector B {X-direction, Y-direction}   = "			<< m_b[0] << ' ' << m_b[1]   << std::endl;
    out << "Properties for Diffusion: \n";
    out << "\t Initial Permeability   = "											<< m_k_0		<< std::endl;
    out << "\t Fluid Viscosity   = "										<< m_eta	<< std::endl;
    out << "\t Constrained specific storage at constant strain Se = "		<< m_Se		<< std::endl;
    out << "Class properties :";
    TPZMaterial::Print(out);
    out << "\n";
    
}

/** Returns the variable index associated with the name */
template<class T,class TMEM>
int TPZPMRSCouplPoroPlast<T,TMEM>::VariableIndex(const std::string &name)
{
    //	Total Strain Variables
    if(!strcmp("et_v",name.c_str()))          return TPZPMRSCouplPoroPlast<T,TMEM>::InxTotStrainVol;
    if(!strcmp("et_x",name.c_str()))          return TPZPMRSCouplPoroPlast<T,TMEM>::InxTotStrainXX;
    if(!strcmp("et_y",name.c_str()))          return TPZPMRSCouplPoroPlast<T,TMEM>::InxTotStrainYY;
    if(!strcmp("et_z",name.c_str()))          return TPZPMRSCouplPoroPlast<T,TMEM>::InxTotStrainZZ;
    if(!strcmp("et_xy",name.c_str()))         return TPZPMRSCouplPoroPlast<T,TMEM>::InxTotStrainXY;
    if(!strcmp("et_xz",name.c_str()))         return TPZPMRSCouplPoroPlast<T,TMEM>::InxTotStrainXZ;
    if(!strcmp("et_yz",name.c_str()))         return TPZPMRSCouplPoroPlast<T,TMEM>::InxTotStrainYZ;
    
    //	Elastic Strain Variables
    if(!strcmp("e_v",name.c_str()))           return TPZPMRSCouplPoroPlast<T,TMEM>::InxElStrainVol;
    if(!strcmp("e_x",name.c_str()))           return TPZPMRSCouplPoroPlast<T,TMEM>::InxElStrainXX;
    if(!strcmp("e_y",name.c_str()))           return TPZPMRSCouplPoroPlast<T,TMEM>::InxElStrainYY;
    if(!strcmp("e_z",name.c_str()))           return TPZPMRSCouplPoroPlast<T,TMEM>::InxElStrainZZ;
    if(!strcmp("e_xy",name.c_str()))          return TPZPMRSCouplPoroPlast<T,TMEM>::InxElStrainXY;
    if(!strcmp("e_xz",name.c_str()))          return TPZPMRSCouplPoroPlast<T,TMEM>::InxElStrainXZ;
    if(!strcmp("e_yz",name.c_str()))          return TPZPMRSCouplPoroPlast<T,TMEM>::InxElStrainYZ;
    
    //	Plastic Strain Variables
    if(!strcmp("ep_v",name.c_str()))          return TPZPMRSCouplPoroPlast<T,TMEM>::InxPlStrainVol;
    if(!strcmp("ep_x",name.c_str()))          return TPZPMRSCouplPoroPlast<T,TMEM>::InxPlStrainXX;
    if(!strcmp("ep_y",name.c_str()))          return TPZPMRSCouplPoroPlast<T,TMEM>::InxPlStrainYY;
    if(!strcmp("ep_z",name.c_str()))          return TPZPMRSCouplPoroPlast<T,TMEM>::InxPlStrainZZ;
    if(!strcmp("ep_xy",name.c_str()))         return TPZPMRSCouplPoroPlast<T,TMEM>::InxPlStrainXY;
    if(!strcmp("ep_xz",name.c_str()))         return TPZPMRSCouplPoroPlast<T,TMEM>::InxPlStrainXZ;
    if(!strcmp("ep_yz",name.c_str()))         return TPZPMRSCouplPoroPlast<T,TMEM>::InxPlStrainYZ;
    
    //	Displacement Variables
    if(!strcmp("u",name.c_str()))		      return TPZPMRSCouplPoroPlast<T,TMEM>::InxDisplacement;
    
    //	Diffusion Variables
    if(!strcmp("p",name.c_str()))		      return TPZPMRSCouplPoroPlast<T,TMEM>::InxPorePressure;
    if(!strcmp("v",name.c_str()))			  return TPZPMRSCouplPoroPlast<T,TMEM>::InxVelocity;
    if(!strcmp("phi",name.c_str()))			  return TPZPMRSCouplPoroPlast<T,TMEM>::InxPorosity;
    if(!strcmp("k_x",name.c_str()))		      return TPZPMRSCouplPoroPlast<T,TMEM>::InxPermeabilityXX;
    if(!strcmp("k_y",name.c_str()))		      return TPZPMRSCouplPoroPlast<T,TMEM>::InxPermeabilityYY;
    if(!strcmp("k_z",name.c_str()))		      return TPZPMRSCouplPoroPlast<T,TMEM>::InxPermeabilityZZ;
    
    //	Total Stress Variables
    if(!strcmp("s_x",name.c_str()))           return TPZPMRSCouplPoroPlast<T,TMEM>::InxTotStressXX;
    if(!strcmp("s_y",name.c_str()))           return TPZPMRSCouplPoroPlast<T,TMEM>::InxTotStressYY;
    if(!strcmp("s_z",name.c_str()))           return TPZPMRSCouplPoroPlast<T,TMEM>::InxTotStressZZ;
    if(!strcmp("t_xy",name.c_str()))          return TPZPMRSCouplPoroPlast<T,TMEM>::InxTotStressXY;
    if(!strcmp("t_xz",name.c_str()))          return TPZPMRSCouplPoroPlast<T,TMEM>::InxTotStressXZ;
    if(!strcmp("t_yz",name.c_str()))          return TPZPMRSCouplPoroPlast<T,TMEM>::InxTotStressYZ;
    
    //	Stress Ratio Variable
    if(!strcmp("K_0",name.c_str()))           return TPZPMRSCouplPoroPlast<T,TMEM>::InxStressRatio;
    
    //	Yield Surface Variable
    if(!strcmp("YS_1",name.c_str()))         return TPZPMRSCouplPoroPlast<T,TMEM>::InxYieldSurface1;
    if(!strcmp("YS_2",name.c_str()))         return TPZPMRSCouplPoroPlast<T,TMEM>::InxYieldSurface2;
    if(!strcmp("YS_3",name.c_str()))         return TPZPMRSCouplPoroPlast<T,TMEM>::InxYieldSurface3;
    
    return TPZMaterial::VariableIndex(name);
}

template<class T,class TMEM>
int TPZPMRSCouplPoroPlast<T,TMEM>::NSolutionVariables(int var)
{
    if(TPZPMRSCouplPoroPlast<T,TMEM>::InxTotStrainVol)          return 1;
    if(TPZPMRSCouplPoroPlast<T,TMEM>::InxTotStrainXX)           return 1;
    if(TPZPMRSCouplPoroPlast<T,TMEM>::InxTotStrainYY)           return 1;
    if(TPZPMRSCouplPoroPlast<T,TMEM>::InxTotStrainZZ)           return 1;
    if(TPZPMRSCouplPoroPlast<T,TMEM>::InxTotStrainXY)           return 1;
    if(TPZPMRSCouplPoroPlast<T,TMEM>::InxTotStrainXZ)           return 1;
    if(TPZPMRSCouplPoroPlast<T,TMEM>::InxTotStrainYZ)           return 1;
    
    //	Elastic Strain Variables
    if(TPZPMRSCouplPoroPlast<T,TMEM>::InxElStrainVol)           return 1;
    if(TPZPMRSCouplPoroPlast<T,TMEM>::InxElStrainXX)            return 1;
    if(TPZPMRSCouplPoroPlast<T,TMEM>::InxElStrainYY)            return 1;
    if(TPZPMRSCouplPoroPlast<T,TMEM>::InxElStrainZZ)            return 1;
    if(TPZPMRSCouplPoroPlast<T,TMEM>::InxElStrainXY)            return 1;
    if(TPZPMRSCouplPoroPlast<T,TMEM>::InxElStrainXZ)            return 1;
    if(TPZPMRSCouplPoroPlast<T,TMEM>::InxElStrainYZ)            return 1;
    
    //	Plastic Strain Variables
    if(TPZPMRSCouplPoroPlast<T,TMEM>::InxPlStrainVol)           return 1;
    if(TPZPMRSCouplPoroPlast<T,TMEM>::InxPlStrainXX)            return 1;
    if(TPZPMRSCouplPoroPlast<T,TMEM>::InxPlStrainYY)            return 1;
    if(TPZPMRSCouplPoroPlast<T,TMEM>::InxPlStrainZZ)            return 1;
    if(TPZPMRSCouplPoroPlast<T,TMEM>::InxPlStrainXY)            return 1;
    if(TPZPMRSCouplPoroPlast<T,TMEM>::InxPlStrainXZ)            return 1;
    if(TPZPMRSCouplPoroPlast<T,TMEM>::InxPlStrainYZ)            return 1;
    
    //	Displacement Variables
    if(TPZPMRSCouplPoroPlast<T,TMEM>::InxDisplacement)		    return m_Dim;
    
    //	Diffusion Variables
    if(TPZPMRSCouplPoroPlast<T,TMEM>::InxPorePressure)		   return 1;
    if(TPZPMRSCouplPoroPlast<T,TMEM>::InxVelocity)			   return m_Dim;
    if(TPZPMRSCouplPoroPlast<T,TMEM>::InxPorosity)			   return 1;
    if(TPZPMRSCouplPoroPlast<T,TMEM>::InxPermeabilityXX)	   return 1;
    if(TPZPMRSCouplPoroPlast<T,TMEM>::InxPermeabilityYY)	   return 1;
    if(TPZPMRSCouplPoroPlast<T,TMEM>::InxPermeabilityZZ)	   return 1;
    
    //	Total Stress Variables
    if(TPZPMRSCouplPoroPlast<T,TMEM>::InxTotStressXX)          return 1;
    if(TPZPMRSCouplPoroPlast<T,TMEM>::InxTotStressYY)          return 1;
    if(TPZPMRSCouplPoroPlast<T,TMEM>::InxTotStressZZ)          return 1;
    if(TPZPMRSCouplPoroPlast<T,TMEM>::InxTotStressXY)          return 1;
    if(TPZPMRSCouplPoroPlast<T,TMEM>::InxTotStressXZ)          return 1;
    if(TPZPMRSCouplPoroPlast<T,TMEM>::InxTotStressYZ)          return 1;
    
    //	Stress Ratio Variable
    if(TPZPMRSCouplPoroPlast<T,TMEM>::InxStressRatio)          return 1;
    
    //	Yield Surface Variable
    if(TPZPMRSCouplPoroPlast<T,TMEM>::InxYieldSurface1)        return 1;
    if(TPZPMRSCouplPoroPlast<T,TMEM>::InxYieldSurface2)        return 1;
    if(TPZPMRSCouplPoroPlast<T,TMEM>::InxYieldSurface3)        return 1;
    
    
    return TPZMaterial::NSolutionVariables(var);
}

//	Calculate Secondary variables based on ux, uy, Pore pressure and their derivatives
template<class T,class TMEM>
void TPZPMRSCouplPoroPlast<T,TMEM>::Solution(TPZVec<TPZMaterialData > &datavec, int var, TPZVec<STATE> &Solout)
{

    
    
    
}



/** @brief computation of effective sigma */
template<class T,class TMEM>
void TPZPMRSCouplPoroPlast<T,TMEM>::Compute_Sigma_n(TPZFMatrix<REAL> Grad_u_n, TPZFMatrix<REAL> Grad_u, TPZFMatrix<REAL> &e_e, TPZFMatrix<REAL> &e_p, TPZFMatrix<REAL> &S)
{
    
#ifdef PZDEBUG

    if (Grad_u.Rows() != 3 && Grad_u.Cols() != 3)
    {
        DebugStop();
    }
    
    if (Grad_u_n.Rows() != 3 && Grad_u_n.Cols() != 3)
    {
        DebugStop();
    }

    if (e_p.Rows() != 3 && e_p.Cols() != 3)
    {
        DebugStop();
    }

    if (e_e.Rows() != 3 && e_e.Cols() != 3)
    {
        DebugStop();
    }
#endif
    
   
    TPZFNMatrix<9,REAL> Grad_du, Grad_du_Transpose = Grad_u, delta_e;
    
    //
    Grad_u_n = Grad_u;
    Grad_du = Grad_u_n; // Linear case
    Grad_du.Transpose(&Grad_du_Transpose);
    delta_e = Grad_du + Grad_du_Transpose;
    delta_e *= 0.5;
    
    
    TPZFNMatrix<9,REAL> e_t, e_tn;
    TPZFNMatrix<9,REAL> S_tn,s_tn, I(delta_e.Rows(),delta_e.Cols(),0.0);
    I.Identity();
    
    /** Trial strain */
    e_t = e_e + e_p;
    e_tn = e_t + delta_e;

    /** Trial stress */
    REAL trace = (e_tn(0,0) + e_tn(1,1) + e_tn(2,2));
    s_tn = 2.0 * m_mu * e_tn + m_lambda * trace * I;
    
    // convert to principal stresses
    Principal_Stress(s_tn, S_tn);
    
    /** Elastic update */
    e_e = e_tn;
    S = s_tn;
    
    return;
    
}

/** @brief Principal Stress */
template<class T,class TMEM>
void TPZPMRSCouplPoroPlast<T,TMEM>::Principal_Stress(TPZFMatrix<REAL> S, TPZFMatrix<REAL> & PrinStres)
{
    
#ifdef PZDEBUG
    
    if (S.Rows() != 3 && S.Cols() != 3)
    {
        DebugStop();
    }
    
#endif
    
    S += 1.0e-18;
    
    REAL a,b,c,d;
    a = 1.0;
    b = - S(0,0) - S(1,1) - S(2,2);
    c = - S(0,1)*S(0,1) - 2.0*S(0,2)*S(0,2) + S(0,0) * S(1,1) + S(0,0) * S(2,2) + S(1,1) * S(2,2);
    d = S(0,0) * S(0,2)*S(0,2) - 2.0* S(0,1) * S(0,2)*S(0,2) + S(0,2)*S(0,2) * S(1,1) + S(0,1)*S(0,1)*S(2,2) - S(0,0) * S(1,1) * S(2,2);
    REAL p,q;
    
    p = (3.0 * a * c - b * b) / (3.0 * a * a);
    q = (2.0 * b * b * b - 9.0 * a * b * c + 27.0 * a * a * d) / (27.0 * a * a * a);
    
    REAL A ,B, C;
    A = 2.0*sqrt(-p/3.0);
    B = -b/(3.0*a);
    C = acos(3.0*(q/(A*p)));
             
    TPZManVector<REAL,3> r(3,0.0);
    r[0] = A*cos((1.0/3.0) * (C+0.0*M_PI))+B;
    r[1] = A*cos((1.0/3.0) * (C+2.0*M_PI))+B;
    r[2] = A*cos((1.0/3.0) * (C+4.0*M_PI))+B;
    
    // sorting
    REAL prins1 = std::max(r[0], std::max(r[1], r[2]));
    REAL prins3 = std::min(r[0], std::min(r[1], r[2]));
    REAL prins2 = 0.0;
    for (int i = 0; i < 3 ; i++)
    {
        if(fabs(r[i]  - prins1) <= 1.0e-10 || fabs(r[i] - prins3) <= 1.0e-10)
        {
            continue;
        }
        prins2 = r[i];
    }
    
    PrinStres.Resize(3, 3);
    PrinStres.Zero();
    PrinStres(0,0) = prins1;
    PrinStres(1,1) = prins2;
    PrinStres(2,2) = prins3;
    
}

// ****************************************************************** plasticity ***********************

template<class T,class TMEM>
void TPZPMRSCouplPoroPlast<T,TMEM>::SetRunPlasticity(bool IsPlasticity)
{
    m_SetRunPlasticity = IsPlasticity;
}


#include "TPZSandlerExtended.h"
#include "TPZPlasticStepPV.h"
#include "TPZYCMohrCoulombPV.h"
#include "TPZSandlerDimaggio.h"

template class TPZPMRSCouplPoroPlast<TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2>, TPZElastoPlasticMem>;
template class TPZPMRSCouplPoroPlast<TPZPlasticStepPV<TPZYCMohrCoulombPV,TPZElasticResponse> , TPZElastoPlasticMem>;