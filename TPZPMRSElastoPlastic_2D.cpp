//
//  TPZPMRSCoupling.cpp
//  PZ
//
//  Created by Manouchehr on Jun 27, 2018.
//
//

#include "TPZPMRSElastoPlastic_2D.h"
#include <iostream>
#include <string>
#include "pzelasmat.h"
#include "pzbndcond.h"
#include "pzaxestools.h"
#include "pzintel.h"
#include "TPZElasticResponse.h"


#include "pzfmatrix.h"
#include "TPZTensor.h"


#ifdef LOG4CXX
#include "pzlog.h"
static LoggerPtr logger(Logger::getLogger("pz.reducedspace.data"));
#endif




template<class T,class TMEM>
TPZPMRSElastoPlastic_2D<T,TMEM>::TPZPMRSElastoPlastic_2D() : TPZMatElastoPlastic2D<T,TMEM>()
{
    m_matId = 0;
    m_Dim = 2;
    m_b.resize(2);
    
    m_b[0]=0.;
    m_b[1]=0.;
    
    m_PlaneStress = 1.;
    
    m_rho_s = 2700.0;
    m_rho_f = 1000.0;
    
    m_k_model = 0;
    
    this->SetCurrentState();

    
    m_SetRunPlasticity = false;

    
}

template<class T,class TMEM>
TPZPMRSElastoPlastic_2D<T,TMEM>::TPZPMRSElastoPlastic_2D(int matid, int dim) : TPZMatElastoPlastic2D<T,TMEM>(matid,1)
{
   
    m_matId = matid;
    m_Dim = dim;
    m_b.resize(2);
    
    m_b[0]=0.;
    m_b[1]=0.;
    
    m_PlaneStress = 1;
    
    m_rho_s = 2700.0;
    m_rho_f = 1000.0;
    
    m_k_model = 0;
 
    
    this->SetCurrentState();
    m_SetRunPlasticity = false;
    
    
}

template<class T,class TMEM>
TPZPMRSElastoPlastic_2D<T,TMEM>::~TPZPMRSElastoPlastic_2D()
{
}


/** @brief copy constructor $ */
template<class T,class TMEM>
TPZPMRSElastoPlastic_2D<T,TMEM>::TPZPMRSElastoPlastic_2D(const TPZPMRSElastoPlastic_2D& other)
{
    this->m_Dim    = other.m_Dim;
    this->m_SimulationData    = other.m_SimulationData;
}


/** @brief Copy assignemnt operator $ */
template<class T,class TMEM>
TPZPMRSElastoPlastic_2D<T,TMEM> &TPZPMRSElastoPlastic_2D<T,TMEM>::operator = (const TPZPMRSElastoPlastic_2D& other)
{
    
    if (this != & other) // prevent self-assignment
    {
        this->m_Dim    = other.m_Dim;
        this->m_SimulationData    = other.m_SimulationData;
    }
    return *this;
}

template<class T,class TMEM>
int TPZPMRSElastoPlastic_2D<T,TMEM>::NStateVariables()
{
    return 2;
}


/** @brief Poroelastic porosity correction */
template<class T,class TMEM>
REAL TPZPMRSElastoPlastic_2D<T,TMEM>::porosity_corrected(TPZVec<TPZMaterialData> &datavec)
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
REAL TPZPMRSElastoPlastic_2D<T,TMEM>::porosity_corrected_3D(TPZVec<TPZMaterialData> &datavec)
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

/** @brief permeability coupling models  */
template<class T,class TMEM>
REAL TPZPMRSElastoPlastic_2D<T,TMEM>::k_permeability(REAL &phi, REAL &k)
{
    
    
    k = 0.0;
    REAL tom2 = 9.869233e-16;
    switch (m_k_model)
    {
        case 0:
        {
            k = m_k_0;
        }
            break;
            
        case 1:
        {
            k = m_k_0*pow((phi/m_porosity_0),4.0);
        }
            break;
            
        case 2:
        {
            k = 0.136*(pow(phi,1.4))*tom2;
        }
            break;
            
        case 3:
        {
            k = (100.0*pow(phi,2.25))*(100.0*pow(phi,2.25))*tom2;
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



template<class T,class TMEM>
void TPZPMRSElastoPlastic_2D<T,TMEM>::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef)
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
    TPZManVector<REAL,2> u = datavec[u_b].sol[0];
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
    REAL phi_poro = porosity_corrected(datavec);
    
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
    TPZFNMatrix<6,REAL> S;

    TPZFNMatrix<6,REAL> Grad_vx_i(2,1,0.0);
    TPZFNMatrix<6,REAL> Grad_vy_i(2,1,0.0);
    
    TPZFNMatrix<6,REAL> Grad_v(2,2,0.0);
    TPZFNMatrix<6,REAL> Grad_vx_j(2,1,0.0);
    TPZFNMatrix<6,REAL> Grad_vy_j(2,1,0.0);
    
    TPZFMatrix<REAL> & S_0 = m_SimulationData->Sigma_0();
    
    S_0.Zero();
    S -= S_0; // Applying prestress
    
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
    k_permeability(phi_poro,k);
    m_lambdau *=1.1;
    REAL c = 1.0;//(k/feta)*(flambdau-flambda)*(flambda + 2.0*fmu)/(falpha*falpha*(flambdau + 2.0*fmu));
    
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
        sout << "<<< TPZMatElastoPlastic2D<T,TMEM>::Contribute ***";
        sout << " Resultant rhs vector:\n" << ef;
        sout << " Resultant stiff vector:\n" << ek;
        LOGPZ_DEBUG(logger,sout.str().c_str());
    }
#endif
    
    
    
    if (m_SetRunPlasticity)
    {
        ContributePlastic(datavec[0],weight,ek,ef);
        
#ifdef LOG4CXX
        if(logger->isDebugEnabled())
        {
            std::stringstream sout;
            sout << "<<< TPZH1PlasticFrac2D<T,TMEM>::Contribute ***";
            sout << " Resultant rhs vector:\n" << ef;
            sout << " Resultant stiff vector:\n" << ek;
            LOGPZ_DEBUG(logger,sout.str().c_str());
        }
#endif
        
        return;
    } 

}

template<class T,class TMEM>
void TPZPMRSElastoPlastic_2D<T,TMEM>::ContributePlastic(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef)
{
    TPZMatElastoPlastic2D<T,TMEM>::Contribute(data,weight,ek,ef);
    return;
    
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
  
  

  if (TPZMatWithMem<TMEM>::fUpdateMem && data.sol.size() > 1)
  {
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
  }
  else
  {
      
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
        sout << ">>> TPZMatElastoPlastic<T,TMEM>::Contribute ***";
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
        
        for( int ju = 0; ju < n_phi_u; ju++)
        {
            
            val  = 2. * Dep(0,0) * dphiXY(0,iu)*dphiXY(0,ju);//dphixdxI*dphixdxJ
            val +=      Dep(0,2) * dphiXY(0,iu)*dphiXY(1,ju);//dphixdxI*dphixdyJ
            val += 2. * Dep(2,0) * dphiXY(1,iu)*dphiXY(0,ju);//dphixdyI*dphixdxJ
            val +=      Dep(2,2) * dphiXY(1,iu)*dphiXY(1,ju);//dphixdyI*dphixdyJ
            val *= 0.5;
            ek(2*iu+0 + first_u, 2*ju+0 + first_u) += weight * val;
            
            val  =      Dep(0,2) * dphiXY(0,iu)*dphiXY(2,ju);//dphixdxI*dphiydxJ
            val += 2. * Dep(0,1) * dphiXY(0,iu)*dphiXY(3,ju);//dphixdxI*dphiydyJ
            val +=      Dep(2,2) * dphiXY(1,iu)*dphiXY(2,ju);//dphixdyI*dphiydxJ
            val += 2. * Dep(2,1) * dphiXY(1,iu)*dphiXY(3,ju);//dphixdyI*dphiydyJ
            val *= 0.5;
            ek(2*iu+0 + first_u, 2*ju+1 + first_u) += weight * val;
      
            val  = 2. * Dep(2,0) * dphiXY(2,iu)*dphiXY(0,ju);//dphiydxI*dphixdxJ
            val +=      Dep(2,2) * dphiXY(2,iu)*dphiXY(1,ju);//dphiydxI*dphixdyJ
            val += 2. * Dep(1,0) * dphiXY(3,iu)*dphiXY(0,ju);//dphiydyI*dphixdxJ
            val +=      Dep(1,2) * dphiXY(3,iu)*dphiXY(1,ju);//dphiydyI*dphixdyJ
            val *= 0.5;
            ek(2*iu+1 + first_u, 2*ju   + first_u) += weight * val;
      
            val  =      Dep(2,2) * dphiXY(2,iu)*dphiXY(2,ju);//dphiydxI*dphiydxJ
            val += 2. * Dep(2,1) * dphiXY(2,iu)*dphiXY(3,ju);//dphiydxI*dphiydyJ
            val +=      Dep(1,2) * dphiXY(3,iu)*dphiXY(2,ju);//dphiydyI*dphiydxJ
            val += 2. * Dep(1,1) * dphiXY(3,iu)*dphiXY(3,ju);//dphiydyI*dphiydyJ
            val *= 0.5;
            ek(2*iu+1 + first_u, 2*ju+1 + first_u) += weight * val;
            
        }
    }
  
    
#ifdef LOG4CXX
  if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << "<<< TPZH1PlasticFrac2D<T,TMEM>::ContributePlastic ***";
        sout << " Resultant rhs vector:\n" << ef;
        sout << " Resultant stiff vector:\n" << ek;
        LOGPZ_DEBUG(logger,sout.str().c_str());
    }
#endif

}

/** @brief of contribute  */
template<class T,class TMEM>
void TPZPMRSElastoPlastic_2D<T,TMEM>::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef){
    
    TPZFMatrix<STATE>  ek_fake(ef.Rows(),ef.Rows(),0.0);
    this->Contribute(datavec, weight, ek_fake, ef);
    
}

template<class T,class TMEM>
void TPZPMRSElastoPlastic_2D<T,TMEM>::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
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


    
    REAL v[3];
    v[0] = bc.Val2()(0,0);	//	Ux displacement
    v[1] = bc.Val2()(1,0);	//	Uy displacement
    v[2] = bc.Val2()(2,0);	//	Pressure
    
    REAL time = this->m_SimulationData->t();
    REAL Value = bc.Val2()(0,0);
    if (bc.HasTimedependentBCForcingFunction())
    {
        TPZManVector<REAL,3> f(3);
        TPZFMatrix<REAL> gradf;
        bc.TimedependentBCForcingFunction()->Execute(datavec[u_b].x, time, f, gradf);
        v[0] = f[0];	//	Ux displacement or Tnx
        v[1] = f[1];	//	Uy displacement or Tny
        v[2] = f[2];	//	Pressure or Qn
    }
    else{
        Value = bc.Val2()(0,0);
    }
    
    
    
    // Dirichlet in Pressure
    switch (bc.Type())
    {
        case 0 : // Du_Dp
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
            for(in = 0 ; in < phrp; in++)
            {
                //    Contribution for load Vector
                ef(in+2*phru,0)        += BIGNUMBER*(p[0]-v[2])*phip(in,0)*weight;    // P Pressure
                
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
            
            //    Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //    Contribution for load Vector
                ef(2*in,0)        += BIGNUMBER*(u[0]-v[0])*phiu(in,0)*weight;    // X displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(2*in,2*jn)        += BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;    // X displacement
                }
            }
            
            //    Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //    Contribution for load Vector
                ef(in+2*phru,0)        += BIGNUMBER*(p[0]-v[2])*phip(in,0)*weight;    // P Pressure
                
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
            
            //    Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //    Contribution for load Vector
                ef(2*in+1,0)      += BIGNUMBER*(u[1]-v[1])*phiu(in,0)*weight;    // y displacement
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(2*in+1,2*jn+1)    += BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;    // Y displacement
                }
            }
            
            //    Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //    Contribution for load Vector
                ef(in+2*phru,0)        += BIGNUMBER*(p[0]-v[2])*phip(in,0)*weight;    // P Pressure
                
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
            
            //    Neumann condition for each state variable
            //    Elasticity Equation
            for(in = 0 ; in <phru; in++)
            {
                //    Normal Tension Components on neumman boundary
                ef(2*in  ,0)    += -1.0 * weight * v[0] * phiu(in,0);        //    Tnx
                ef(2*in+1,0)    += -1.0 * weight * v[1] * phiu(in,0);        //    Tny
            }
            
            //    Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //    Contribution for load Vector
                ef(in+2*phru,0)        += BIGNUMBER*(p[0]-v[2])*phip(in,0)*weight;    // P Pressure
                
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
            
            REAL tn = v[0];  //    Tn normal traction
            
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
            for(in = 0 ; in < phrp; in++)
            {
                //    Contribution for load Vector
                ef(in+2*phru,0)        += BIGNUMBER*(p[0]-v[2])*phip(in,0)*weight;    // P Pressure
                
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
            
            //    Elasticity Equation
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
            
            //    Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //    Contribution for load Vector
                ef(2*in,0)        += BIGNUMBER*(u[0]-v[0])*phiu(in,0)*weight;    // X displacement Value
                
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
                ef(in+2*phru,0)    += 1.0 * weight * v[2] * phip(in,0);    // Qnormal
            }
            break;
        }
            
        case 7 : // Duy_Nq
        {
            
            
            for(in = 0 ; in < phru; in++)
            {
                //    Contribution for load Vector
                ef(2*in+1,0)      += BIGNUMBER*(u[1]-v[1])*phiu(in,0)*weight;    // y displacement Value
                
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
                ef(in+2*phru,0)    += 1.0 * weight * v[2] * phip(in,0);    // Qnormal
            }
            break;
        }
            
        case 8 : // Nt_Nq
        {
            
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
            
            REAL tn = v[0];  //    Tn normal traction
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
                ef(in+2*phru,0)    += 1.0 * weight * v[2] * phip(in,0);    // Qnormal
            }
            break;
        }
            
        case 10 : //Duy_time_Dp
        {
            
            //    Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //    Contribution for load Vector
                ef(2*in+1,0)      += BIGNUMBER*(u[1]-v[1])*phiu(in,0)*weight;    // y displacement
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(2*in+1,2*jn+1)    += BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;    // Y displacement
                }
            }
            
            //    Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //    Contribution for load Vector
                ef(in+2*phru,0)        += BIGNUMBER*(p[0]-v[2])*phip(in,0)*weight;    // P Pressure
                
                for (jn = 0 ; jn < phrp; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(in+2*phru,jn+2*phru)        += BIGNUMBER*phip(in,0)*phip(jn,0)*weight;    // P Pressure
                }
            }
            break;
        }
            
        case 11 : // Ntny_time_Dp
        {
            
            //    Neumann condition for each state variable
            //    Elasticity Equation
            for(in = 0 ; in <phru; in++)
            {
                //    Normal Tension Components on neumman boundary
                ef(2*in+1,0)    += -1.0 * weight * v[1] * phiu(in,0);        //    Tny
            }
            
            //    Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //    Contribution for load Vector
                ef(in+2*phru,0)        += BIGNUMBER*(p[0]-v[2])*phip(in,0)*weight;    // P Pressure
                
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
void TPZPMRSElastoPlastic_2D<T,TMEM>::FillDataRequirements(TPZVec<TPZMaterialData > &datavec)
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
void TPZPMRSElastoPlastic_2D<T,TMEM>::FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec)
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
void TPZPMRSElastoPlastic_2D<T,TMEM>::Print(std::ostream &out)
{
    out << "Material Name : " << Name() << "\n";
    out << "Plane Problem (fPlaneStress = 0, for Plane Strain conditions) " << m_PlaneStress << std::endl;
    out << "Properties for TPZPMRSCoupling: \n";
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
int TPZPMRSElastoPlastic_2D<T,TMEM>::VariableIndex(const std::string &name)
{
    //	Elasticity Variables
    if(!strcmp("u",name.c_str()))				return	1;
    if(!strcmp("s_x",name.c_str()))             return	2;
    if(!strcmp("s_y",name.c_str()))             return	3;
    if(!strcmp("s_z",name.c_str()))             return	4;
    if(!strcmp("t_xy",name.c_str()))            return	5;
    if(!strcmp("t_xz",name.c_str()))            return	6;
    if(!strcmp("t_yz",name.c_str()))            return	7;
    
    //	Diffusion Variables
    if(!strcmp("p",name.c_str()))				return	8;
    if(!strcmp("v",name.c_str()))				return	9;
    if(!strcmp("k_x",name.c_str()))				return	10;
    if(!strcmp("k_y",name.c_str()))				return	11;
    if(!strcmp("k_z",name.c_str()))				return	12;
    if(!strcmp("phi",name.c_str()))				return	13;
    
    //	Defformation Variables
    if(!strcmp("e_x",name.c_str()))             return	14;
    if(!strcmp("e_y",name.c_str()))             return	15;
    if(!strcmp("e_z",name.c_str()))             return	16;
    if(!strcmp("e_xy",name.c_str()))            return	17;
    if(!strcmp("e_xz",name.c_str()))            return	18;
    if(!strcmp("e_yz",name.c_str()))            return	19;
    if(!strcmp("ep_x",name.c_str()))            return	20;
    if(!strcmp("ep_y",name.c_str()))            return	21;
    if(!strcmp("ep_z",name.c_str()))            return	22;
    if(!strcmp("ep_xy",name.c_str()))           return	23;
    if(!strcmp("ep_xz",name.c_str()))           return	24;
    if(!strcmp("ep_yz",name.c_str()))           return	25;
    
    //	Stress Ratio Variable
    if(!strcmp("K_0",name.c_str()))             return	26;
    
    return TPZMaterial::VariableIndex(name);
}


template<class T,class TMEM>
int TPZPMRSElastoPlastic_2D<T,TMEM>::NSolutionVariables(int var)
{
    if(var == 1)	return m_Dim;
    if(var == 2)	return 1;
    if(var == 3)	return 1;
    if(var == 4)	return 1;
    if(var == 5)	return 1;
    if(var == 6)	return 1;
    if(var == 7)	return 1;
    if(var == 8)	return 1;
    if(var == 9)	return m_Dim;
    if(var == 10)	return 1;
    if(var == 11)	return 1;
    if(var == 12)	return 1;
    if(var == 13)	return 1;
    if(var == 14)	return 1;
    if(var == 15)	return 1;
    if(var == 16)	return 1;
    if(var == 17)	return 1;
    if(var == 18)	return 1;
    if(var == 19)	return 1;
    if(var == 20)	return 1;
    if(var == 21)	return 1;
    if(var == 22)	return 1;
    if(var == 23)	return 1;
    if(var == 24)	return 1;
    if(var == 25)	return 1;
    if(var == 26)	return 1;
    
    return TPZMaterial::NSolutionVariables(var);
}


//	Calculate Secondary variables based on ux, uy, Pore pressure and their derivatives
template<class T,class TMEM>
void TPZPMRSElastoPlastic_2D<T,TMEM>::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout)
{
    
    Solout.Resize( this->NSolutionVariables(var));
    
    int u_b = 0;
    int p_b = 1;
    
    // Getting the space functions
    TPZFNMatrix <9,REAL>	&axes_u	=	datavec[u_b].axes;
    TPZFNMatrix <9,REAL>	&axes_p	=	datavec[p_b].axes;
    
    // Getting the solutions and derivatives
    TPZManVector<REAL,3> u = datavec[u_b].sol[0];
    TPZManVector<REAL,1> p = datavec[p_b].sol[0];
    
    TPZFNMatrix <9,REAL> du = datavec[u_b].dsol[0];
    TPZFNMatrix <9,REAL> dp = datavec[p_b].dsol[0];
    
    
    REAL to_Mpa     = 1.0e-6;
    REAL to_Darcy   = 1.013249966e+12;
    
    
    // Computing Gradient of the Solution
    TPZFNMatrix<9,REAL> Grad_p(3,1,0.0),Grad_u(3,3,0.0),Grad_u_n(3,3,0.0),e_e(3,3,0.0),e_p(3,3,0.0),S;
    
    
    // Computing Gradient of deformation in 2D for corrector_DP function
    if (m_Dim != 3)
    {
        Grad_u(0,0) = du(0,0)*axes_u(0,0)+du(1,0)*axes_u(1,0); // dux/dx
        Grad_u(0,1) = du(0,0)*axes_u(0,1)+du(1,0)*axes_u(1,1); // dux/dy
        
        Grad_u(1,0) = du(0,1)*axes_u(0,0)+du(1,1)*axes_u(1,0); // duy/dx
        Grad_u(1,1) = du(0,1)*axes_u(0,1)+du(1,1)*axes_u(1,1); // duy/dy
    }
    else
    {
        Grad_u(0,0) = du(0,0)*axes_u(0,0)+du(1,0)*axes_u(1,0)+du(2,0)*axes_u(2,0); // dux/dx
        Grad_u(0,1) = du(0,0)*axes_u(0,1)+du(1,0)*axes_u(1,1)+du(2,0)*axes_u(2,1); // dux/dy
        Grad_u(0,2) = du(0,0)*axes_u(0,2)+du(1,0)*axes_u(1,2)+du(2,0)*axes_u(2,2); // dux/dz
        
        Grad_u(1,0) = du(0,1)*axes_u(0,0)+du(1,1)*axes_u(1,0)+du(2,1)*axes_u(2,0); // duy/dx
        Grad_u(1,1) = du(0,1)*axes_u(0,1)+du(1,1)*axes_u(1,1)+du(2,1)*axes_u(2,1); // duy/dy
        Grad_u(1,2) = du(0,1)*axes_u(0,2)+du(1,1)*axes_u(1,2)+du(2,1)*axes_u(2,2); // duy/dz
        
        Grad_u(2,0) = du(0,2)*axes_u(0,0)+du(1,2)*axes_u(1,0)+du(2,2)*axes_u(2,0); // duz/dx
        Grad_u(2,1) = du(0,2)*axes_u(0,1)+du(1,2)*axes_u(1,1)+du(2,2)*axes_u(2,1); // duz/dy
        Grad_u(2,2) = du(0,2)*axes_u(0,2)+du(1,2)*axes_u(1,2)+du(2,2)*axes_u(2,2); // duz/dz
    }
    
    
    
    
    //	Displacements
    if(var == 1)
    {
        Solout[0] = u[0];
        Solout[1] = u[1];
        if (m_Dim == 3)
        {
            Solout[2] = u[2];
        }
        return;
    }
    
    //	sigma_x
    if(var == 2)
    {
        Solout[0] = S(0,0)*to_Mpa;
        return;
    }
    
    //	sigma_y
    if(var == 3)
    {
        Solout[0] = S(1,1)*to_Mpa;
        return;
    }
    
    //	sigma_z
    if(var == 4)
    {
        Solout[0] = S(2,2)*to_Mpa;
        return;
    }
    
    //	tau_xy
    if(var == 5)
    {
        Solout[0] = S(0,1)*to_Mpa;
        return;
    }
    
    //	tau_xz
    if(var == 6)
    {
        Solout[0] = S(0,2)*to_Mpa;
        return;
    }
    
    //	tau_yz
    if(var == 7)
    {
        Solout[0] = S(1,2)*to_Mpa;
        return;
    }
    
    //	pore pressure
    if(var == 8)
    {
        Solout[0] = p[0]*to_Mpa;
        return;
    }
    
    //	Darcy's velocity
    if(var == 9)
    {
        if (m_Dim != 3)
        {
            Grad_p(0,0) = dp(0,0)*axes_p(0,0)+dp(1,0)*axes_p(1,0); // dp/dx
            Grad_p(1,0) = dp(0,0)*axes_p(0,1)+dp(1,0)*axes_p(1,1); // dp/dy
            
            REAL phi = porosity_corrected(datavec);
            REAL k;
            k_permeability(phi, k);
            Solout[0] = -(k/m_eta) * Grad_p(0,0);
            Solout[1] = -(k/m_eta) * Grad_p(1,0);
            return;
        }
        else
        {
            Grad_p(0,0) = dp(0,0)*axes_p(0,0)+dp(1,0)*axes_p(1,0)+dp(2,0)*axes_p(2,0); // dp/dx
            Grad_p(1,0) = dp(0,0)*axes_p(0,1)+dp(1,0)*axes_p(1,1)+dp(2,0)*axes_p(2,1); // dp/dy
            Grad_p(2,0) = dp(0,0)*axes_p(0,2)+dp(1,0)*axes_p(1,2)+dp(2,0)*axes_p(2,2); // dp/dz
            
            REAL phi = porosity_corrected_3D(datavec);
            REAL k;
            k_permeability(phi, k);
            Solout[0] = -(k/m_eta) * Grad_p(0,0);
            Solout[1] = -(k/m_eta) * Grad_p(1,0);
            Solout[2] = -(k/m_eta) * Grad_p(2,0);
            return;
        }
    }
    
    //	k_x
    if(var == 10)
    {
        if (m_Dim != 3)
        {
            REAL phi = porosity_corrected(datavec);
            REAL k = 0.0;
            k_permeability(phi, k);
            Solout[0] = k*to_Darcy;
            return;
        }
        
        else
        {
            REAL phi = porosity_corrected_3D(datavec);
            REAL k = 0.0;
            k_permeability(phi, k);
            Solout[0] = k*to_Darcy;
            return;
        }
    }
    
    //	k_y
    if(var == 11)
    {
        if (m_Dim != 3)
        {
            REAL phi = porosity_corrected(datavec);
            REAL k = 0.0;
            k_permeability(phi, k);
            Solout[0] = k*to_Darcy;
            return;
        }
        else
        {
            REAL phi = porosity_corrected_3D(datavec);
            REAL k = 0.0;
            k_permeability(phi, k);
            Solout[0] = k*to_Darcy;
            return;
        }
    }
    
    //	k_z
    if(var == 12)
    {
        if (m_Dim != 3)
        {
            REAL phi = porosity_corrected(datavec);
            REAL k = 0.0;
            k_permeability(phi, k);
            Solout[0] = k*to_Darcy;
            return;
        }
        else
        {
            REAL phi = porosity_corrected_3D(datavec);
            REAL k = 0.0;
            k_permeability(phi, k);
            Solout[0] = k*to_Darcy;
            return;
        }
    }
    
    //	Porosity form poroelastic correction
    if(var == 13)
    {
        if (m_Dim != 3)
        {
            Solout[0] = porosity_corrected(datavec);
            return;
        }
        else
        {
            Solout[0] = porosity_corrected_3D(datavec);
            return;
        }
    }
    
    //	epsilon_xx
    if(var == 14)
    {
        Solout[0] = e_e(0,0);
        return;
    }
    
    //	epsilon_yy
    if(var == 15)
    {
        Solout[0] = e_e(1,1);
        return;
    }
    
    //	epsilon_zz
    if(var == 16)
    {
        Solout[0] = e_e(2,2);
        return;
    }
    
    //	epsilon_xy
    if(var == 17)
    {
        Solout[0] = e_e(0,1);
        return;
    }
    
    //	epsilon_xz
    if(var == 18)
    {
        Solout[0] = e_e(0,2);
        return;
    }
    
    //	epsilon_yz
    if(var == 19)
    {
        Solout[0] = e_e(1,2);
        return;
    }
    
    //	epsilon_p_xx
    if(var == 20)
    {
        Solout[0] = e_p(0,0);
        return;
    }
    
    //	epsilon_p_yy
    if(var == 21)
    {
        Solout[0] = e_p(1,1);
        return;
    }
    
    //	epsilon_p_zz
    if(var == 22)
    {
        Solout[0] = e_p(2,2);
        return;
    }
    
    //	epsilon_p_xy
    if(var == 23)
    {
        Solout[0] = e_p(0,1);
        return;
    }
    
    //	epsilon_p_xz
    if(var == 24)
    {
        Solout[0] = e_p(0,2);
        return;
    }
    
    //	epsilon_p_yz
    if(var == 25)
    {
        Solout[0] = e_p(1,2);
        return;
    }
    
    //	K_0
    if(var == 26)
    {
        Solout[0] = S(0,0)/S(1,1);
        return;
    }
    
}



template<class T,class TMEM>
void TPZPMRSElastoPlastic_2D<T,TMEM>::SetRunPlasticity(bool IsPlasticity)
{
  m_SetRunPlasticity = IsPlasticity;
}


#include "TPZSandlerExtended.h"
#include "TPZPlasticStepPV.h"
#include "TPZYCMohrCoulombPV.h"
#include "TPZSandlerDimaggio.h"

template class TPZPMRSElastoPlastic_2D<TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2>, TPZElastoPlasticMem>;
template class TPZPMRSElastoPlastic_2D<TPZPlasticStepPV<TPZYCMohrCoulombPV,TPZElasticResponse> , TPZElastoPlasticMem>;
//template class TPZH1PlasticFrac2D<TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse> , TPZElastoPlasticMem>;
