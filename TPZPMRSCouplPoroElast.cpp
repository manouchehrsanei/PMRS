//
//  TPZPMRSCouplPoroElast.cpp
//  PZ
//
//  Created by Omar and Manouchehr on 8/28/16.
//
//

#include "TPZPMRSCouplPoroElast.h"
#include <iostream>
#include <string>
#include "pzbndcond.h"
#include "pzaxestools.h"
#include <algorithm>
#include "pzlog.h"

#include "pzfmatrix.h"
#include "TPZTensor.h"



#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.TPZPMRSCouplPoroElast"));
#endif



/** @brief default costructor */
TPZPMRSCouplPoroElast::TPZPMRSCouplPoroElast():TPZMatWithMem<TPZPMRSMemoryPoroElast,TPZDiscontinuousGalerkin>(), m_nu(0.), m_alpha(0.), m_k_0(0.), m_eta(0.)
{
    m_Dim = 3;
    m_b.resize(3);
    
    m_b[0]=0.;
    m_b[1]=0.;
    m_b[2]=0.;

    m_k_model = 0;
}

/** @brief costructor based on a material id */
TPZPMRSCouplPoroElast::TPZPMRSCouplPoroElast(int matid, int dim):TPZMatWithMem<TPZPMRSMemoryPoroElast,TPZDiscontinuousGalerkin>(matid), m_nu(0.), m_alpha(0.), m_k_0(0.), m_eta(0.)
{
    m_Dim = dim;
    m_b.resize(3);
    
    m_b[0]=0.;
    m_b[1]=0.;
    m_b[2]=0.;
    
    m_k_model = 0;
}

/** @brief default destructor */
TPZPMRSCouplPoroElast::~TPZPMRSCouplPoroElast()
{
}


/** @brief copy constructor $ */
TPZPMRSCouplPoroElast::TPZPMRSCouplPoroElast(const TPZPMRSCouplPoroElast& other)
{
    this->m_Dim    = other.m_Dim;
    this->m_SimulationData    = other.m_SimulationData;
}


/** @brief Copy assignemnt operator $ */
TPZPMRSCouplPoroElast& TPZPMRSCouplPoroElast::operator = (const TPZPMRSCouplPoroElast& other)
{
    
    if (this != & other) // prevent self-assignment
    {
        this->m_Dim    = other.m_Dim;
        this->m_SimulationData    = other.m_SimulationData;
    }
    return *this;
}


/** @brief number of state variables */
int TPZPMRSCouplPoroElast::NStateVariables()
{
    return 2;
}

/** @brief permeability coupling models  */
REAL TPZPMRSCouplPoroElast::k_permeability(REAL &phi, REAL &k)
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
            
        case 4: // Davies and Davies (2001): Exponential function: C = 1;
        {
            k = m_k_0*exp((phi/m_porosity_0)-1);
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
REAL TPZPMRSCouplPoroElast::porosity_corrected_2D(TPZVec<TPZMaterialData> &datavec)
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
REAL TPZPMRSCouplPoroElast::porosity_corrected_3D(TPZVec<TPZMaterialData> &datavec)
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


/** @brief of contribute of BC */
void TPZPMRSCouplPoroElast::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE>  &ek, TPZFMatrix<STATE> &ef)
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
void TPZPMRSCouplPoroElast::Contribute_2D(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE>  &ek, TPZFMatrix<STATE> &ef)
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
    REAL phi_poro = porosity_corrected_2D(datavec);
    
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
    TPZFNMatrix<9,REAL> Grad_u(3,3,0.0),Grad_u_n,e_e,e_p,S;
    Grad_u(0,0) = du(0,0)*axes_u(0,0)+du(1,0)*axes_u(1,0); // dux/dx
    Grad_u(0,1) = du(0,0)*axes_u(0,1)+du(1,0)*axes_u(1,1); // dux/dy
    
    Grad_u(1,0) = du(0,1)*axes_u(0,0)+du(1,1)*axes_u(1,0); // duy/dx
    Grad_u(1,1) = du(0,1)*axes_u(0,1)+du(1,1)*axes_u(1,1); // duy/dy
    
    // Get the solution at the integrations points
    long global_point_index = datavec[0].intGlobPtIndex;
    TPZAdmChunkVector<TPZPMRSMemoryPoroElast> memory_vec = *GetMemory();
    TPZPMRSMemoryPoroElast &memory = memory_vec[global_point_index];
    e_e = memory.epsilon_e_n();
    e_p = memory.epsilon_p_n();
    Grad_u_n = memory.grad_u_n();
    
    Compute_Sigma_n(Grad_u_n, Grad_u, e_e, e_p, S);
    
    TPZFNMatrix<6,REAL> Grad_vx_i(2,1,0.0);
    TPZFNMatrix<6,REAL> Grad_vy_i(2,1,0.0);

    TPZFNMatrix<6,REAL> Grad_v(2,2,0.0);
    TPZFNMatrix<6,REAL> Grad_vx_j(2,1,0.0);
    TPZFNMatrix<6,REAL> Grad_vy_j(2,1,0.0);

    TPZFMatrix<REAL> Sigma_0;// = m_SimulationData->PreStress();
    
    
    Sigma_0.Zero();
    S -= Sigma_0; // Applying prestress

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
    
    REAL c = (k/m_eta); // (k/m_eta)*(m_lambdau-m_lambda)*(m_lambda + 2.0*m_mu)/(m_alpha*m_alpha*(m_lambdau + 2.0*m_mu));

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
    
}


// Contribute Methods being used

void TPZPMRSCouplPoroElast::Contribute_3D(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
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
    Grad_u(0,0) = du(0,0)*axes_u(0,0)+du(1,0)*axes_u(1,0)+du(2,0)*axes_u(2,0); // dux/dx
    Grad_u(0,1) = du(0,0)*axes_u(0,1)+du(1,0)*axes_u(1,1)+du(2,0)*axes_u(2,1); // dux/dy
    Grad_u(0,2) = du(0,0)*axes_u(0,2)+du(1,0)*axes_u(1,2)+du(2,0)*axes_u(2,2); // dux/dz
    
    Grad_u(1,0) = du(0,1)*axes_u(0,0)+du(1,1)*axes_u(1,0)+du(2,1)*axes_u(2,0); // duy/dx
    Grad_u(1,1) = du(0,1)*axes_u(0,1)+du(1,1)*axes_u(1,1)+du(2,1)*axes_u(2,1); // duy/dy
    Grad_u(1,2) = du(0,1)*axes_u(0,2)+du(1,1)*axes_u(1,2)+du(2,1)*axes_u(2,2); // duy/dz
    
    Grad_u(2,0) = du(0,2)*axes_u(0,0)+du(1,2)*axes_u(1,0)+du(2,2)*axes_u(2,0); // duz/dx
    Grad_u(2,1) = du(0,2)*axes_u(0,1)+du(1,2)*axes_u(1,1)+du(2,2)*axes_u(2,1); // duz/dy
    Grad_u(2,2) = du(0,2)*axes_u(0,2)+du(1,2)*axes_u(1,2)+du(2,2)*axes_u(2,2); // duz/dz
    
    // Get the solution at the integrations points
    long global_point_index = datavec[0].intGlobPtIndex;
    TPZAdmChunkVector<TPZPMRSMemoryPoroElast> memory_vec = *GetMemory();
    TPZPMRSMemoryPoroElast &memory = memory_vec[global_point_index];
    e_e = memory.epsilon_e_n();
    e_p = memory.epsilon_p_n();
    Grad_u_n = memory.grad_u_n();
    
    Compute_Sigma_n(Grad_u_n, Grad_u, e_e, e_p, S);
    
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
    
    TPZFMatrix<REAL> Sigma_0;// = m_SimulationData->PreStress();
    
    Sigma_0.Zero();
    S -= Sigma_0; // Applying prestress
    
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
void TPZPMRSCouplPoroElast::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef){
    
    TPZFMatrix<STATE>  ek_fake(ef.Rows(),ef.Rows(),0.0);
    this->Contribute(datavec, weight, ek_fake, ef);
    
}

/** @brief of contribute of BC */
void TPZPMRSCouplPoroElast::ContributeBC(TPZVec<TPZMaterialData> &datavec,REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
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
void TPZPMRSCouplPoroElast::ContributeBC_2D(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
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

    
    
    // Boundaries

    // Dirichlet in Pressure
    switch (bc.Type())
    {

        case 0 :  // Dp
            // Dirichlet BC  PD
        {
            REAL v[1];
            v[0] = bc.Val2()(0,0);    //    QN
            
            for (int ip = 0; ip < phrp; ip++)
            {
                ef(ip+2*phru,0) +=  weight * gBigNumber * (p[0] - v[0]) * phip(ip,0);
                
                for (int jp = 0; jp < phrp; jp++)
                {
                    ek(ip+2*phru,jp+2*phru) +=  weight * gBigNumber * phip(jp,0) * phip(ip,0);
                }
            }
            break;
        }
            
            
        case 1 :    // Nq
            // Neumann BC  QN
        {
            REAL v[1];
            v[0] = bc.Val2()(0,0);    //    Pressure
            
            for (int ip = 0; ip < phrp; ip++)
            {
                ef(ip+2*phru,0) +=  weight * v[0] * phip(ip,0);
            }
            
            break;
        }
            
            
        case 2 : // Du
            // Dirichlet of displacement
            
        {
            REAL v[2];
            v[0] = bc.Val2()(0,0);    //    Ux displacement
            v[1] = bc.Val2()(1,0);    //    Uy displacement
            
            for(in = 0 ; in < phru; in++)
            {
                //    Contribution for load Vector
                ef(2*in+0,0)      += gBigNumber*(u[0] - v[0])*phiu(in,0)*weight;    // X displacement Value
                ef(2*in+1,0)      += gBigNumber*(u[1] - v[1])*phiu(in,0)*weight;    // Y displacement Value
                
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(2*in+0,2*jn+0)    += gBigNumber*phiu(in,0)*phiu(jn,0)*weight;    // X displacement
                    ek(2*in+1,2*jn+1)    += gBigNumber*phiu(in,0)*phiu(jn,0)*weight;    // Y displacement
                    
                }
            }
            
            break;
        }
            
            
        case 3 : // Dux
            // Dirichlet in x direction of displacement
            
        {
            REAL v[1];
            v[0] = bc.Val2()(0,0);    //    Ux displacement
            
            for(in = 0 ; in < phru; in++)
            {
                //    Contribution for load Vector
                ef(2*in,0)        += gBigNumber*(u[0] - v[0])*phiu(in,0)*weight;    // X displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(2*in,2*jn)        += gBigNumber*phiu(in,0)*phiu(jn,0)*weight;    // X displacement
                }
            }
            
            break;
        }
            
        case 4 : //Duy
            // Dirichlet in y direction of displacement
            
        {
            REAL v[1];
            v[0] = bc.Val2()(0,0);    //    Uy displacement
            
            for(in = 0 ; in < phru; in++)
            {
                //    Contribution for load Vector
                ef(2*in+1,0)      += gBigNumber*(u[1] - v[0])*phiu(in,0)*weight;    // Y displacement
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(2*in+1,2*jn+1)    += gBigNumber*phiu(in,0)*phiu(jn,0)*weight;    // Y displacement
                }
            }
            
            break;
        }
            
        case 5 : // Nt
            // Neumann of traction
            
        {
            REAL v[2];
            v[0] = bc.Val2()(0,0);    //    Tnx
            v[1] = bc.Val2()(1,0);    //    Tny
            
            //    Neumann condition for each state variable
            //    Elasticity Equation
            for(in = 0 ; in <phru; in++)
            {
                //    Normal Tension Components on neumman boundary
                ef(2*in+0,0)    += -1.0 * weight * v[0] * phiu(in,0);        //    Tnx
                ef(2*in+1,0)    += -1.0 * weight * v[1] * phiu(in,0);        //    Tny
            }
            
            break;
        }
            
        case 6 : // Ntn
            // Neumann of traction
            
        {
            REAL v[1];
            v[0] = bc.Val2()(0,0);    //    Tn normal traction
            
            REAL tn = v[0];
            TPZManVector<REAL,2> n = datavec[u_b].normal;
            //    Neumann condition for each state variable
            //    Elasticity Equation
            for(in = 0 ; in <phru; in++)
            {
                //    Normal Tension Components on neumman boundary
                ef(2*in+0,0)    += -1.0 * weight * tn * n[0] * phiu(in,0);        //    Tnx
                ef(2*in+1,0)    += -1.0 * weight * tn * n[1] * phiu(in,0);        //    Tny
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
void TPZPMRSCouplPoroElast::ContributeBC_3D(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    
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
    
    switch (bc.Type())
    {
        case 0 :  // Dp
            // Dirichlet BC  PD
        {
            REAL v[1];
            v[0] = bc.Val2()(0,0);    //    QN
            
            for (int ip = 0; ip < phrp; ip++)
            {
                ef(ip+2*phru,0) +=  weight * gBigNumber * (p[0] - v[0]) * phip(ip,0);
                
                for (int jp = 0; jp < phrp; jp++)
                {
                    ek(ip+2*phru,jp+2*phru) +=  weight * gBigNumber * phip(jp,0) * phip(ip,0);
                }
            }
            break;
        }
            
            
        case 1 :    // Nq
            // Neumann BC  QN
        {
            REAL v[1];
            v[0] = bc.Val2()(0,0);    //    Pressure
            
            for (int ip = 0; ip < phrp; ip++)
            {
                ef(ip+2*phru,0) +=  weight * v[0] * phip(ip,0);
            }
            
            break;
        }
            
            
        case 2 : // Du
            // Dirichlet of displacement
            
        {
            REAL v[3];
            v[0] = bc.Val2()(0,0);    //    Ux displacement
            v[1] = bc.Val2()(1,0);    //    Uy displacement
            v[2] = bc.Val2()(2,0);    //    Uz displacement
            
            
            for(in = 0 ; in < phru; in++)
            {
                //    Contribution for load Vector
                ef(3*in+0,0)      += gBigNumber*(u[0] - v[0])*phiu(in,0)*weight;    // X displacement Value
                ef(3*in+1,0)      += gBigNumber*(u[1] - v[1])*phiu(in,0)*weight;    // Y displacement Value
                ef(3*in+2,0)      += gBigNumber*(u[2] - v[2])*phiu(in,0)*weight;    // Z displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(3*in+0,3*jn+0)    += gBigNumber*phiu(in,0)*phiu(jn,0)*weight;    // X displacement
                    ek(3*in+1,3*jn+1)    += gBigNumber*phiu(in,0)*phiu(jn,0)*weight;    // Y displacement
                    ek(3*in+2,3*jn+2)    += gBigNumber*phiu(in,0)*phiu(jn,0)*weight;    // Z displacement
                    
                }
            }
            
            break;
        }
            
            
        case 3 : // Dux
            // Dirichlet in x direction of displacement
            
        {
            REAL v[1];
            v[0] = bc.Val2()(0,0);    //    Ux displacement
            
            for(in = 0 ; in < phru; in++)
            {
                //    Contribution for load Vector
                ef(3*in,0)        += gBigNumber*(u[0] - v[0])*phiu(in,0)*weight;    // X displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(3*in,3*jn)        += gBigNumber*phiu(in,0)*phiu(jn,0)*weight;    // X displacement
                }
            }
            
            break;
        }
            
        case 4 : //Duy
            // Dirichlet in y direction of displacement
            
        {
            REAL v[1];
            v[0] = bc.Val2()(0,0);    //    Uy displacement
            
            for(in = 0 ; in < phru; in++)
            {
                //    Contribution for load Vector
                ef(3*in+1,0)      += gBigNumber*(u[1] - v[0])*phiu(in,0)*weight;    // Y displacement
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(3*in+1,3*jn+1)    += gBigNumber*phiu(in,0)*phiu(jn,0)*weight;    // Y displacement
                }
            }
            
            break;
        }
            
        case 5 : // Nt
            // Neumann of traction
            
        {
            REAL v[3];
            v[0] = bc.Val2()(0,0);    //    Tnx
            v[1] = bc.Val2()(1,0);    //    Tny
            v[2] = bc.Val2()(2,0);    //    Tnz
            
            
            //    Neumann condition for each state variable
            //    Elasticity Equation
            for(in = 0 ; in <phru; in++)
            {
                //    Normal Tension Components on neumman boundary
                ef(3*in+0,0)    += -1.0 * weight * v[0] * phiu(in,0);        //    Tnx
                ef(3*in+1,0)    += -1.0 * weight * v[1] * phiu(in,0);        //    Tny
                ef(3*in+2,0)    += -1.0 * weight * v[2] * phiu(in,0);        //    Tnz
                
            }
            
            break;
        }
            
        case 6 : // Ntn
            // Neumann of traction
            
        {
            REAL v[1];
            v[0] = bc.Val2()(0,0);    //    Tn normal traction
            
            REAL tn = v[0];
            TPZManVector<REAL,3> n = datavec[u_b].normal;
            //    Neumann condition for each state variable
            //    Elasticity Equation
            for(in = 0 ; in <phru; in++)
            {
                //    Normal Tension Components on neumman boundary
                ef(3*in+0,0)    += -1.0 * weight * tn * n[0] * phiu(in,0);        //    Tnx
                ef(3*in+1,0)    += -1.0 * weight * tn * n[1] * phiu(in,0);        //    Tny
                ef(3*in+2,0)    += -1.0 * weight * tn * n[2] * phiu(in,0);        //    Tny
            }
            
            break;
        }
            
            
            
        case 7 : // Duz
            // Dirichlet in z direction of displacement
        {
            REAL v[1];
            v[0] = bc.Val2()(0,0);    //    Uz displacement
            
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                ef(3*in+2,0)		+= gBigNumber*(u[2] - v[0])*phiu(in,0)*weight;	// Z displacement Value
                
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in+2,3*jn+2)		+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// Z displacement
                    
                }
            }
            
            break;
            
        }
            
            
        case 8 : // Duxy
            // Dirichlet in x and y direction of displacement
        {
            REAL v[2];
            v[0] = bc.Val2()(0,0);    //    Ux displacement
            v[1] = bc.Val2()(1,0);    //    Uy displacement
            
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                ef(3*in+0,0)	+= gBigNumber*(u[0] - v[0])*phiu(in,0)*weight;	// X displacement Value
                ef(3*in+1,0)	+= gBigNumber*(u[1] - v[1])*phiu(in,0)*weight;	// Y displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in+0,3*jn+0)	+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
                    ek(3*in+1,3*jn+1)	+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// Y displacement
                }
            }
            
            break;
            
        }
            
        case 9 : // Duxz
            // Dirichlet in x and z direction of displacement
        {
            REAL v[2];
            v[0] = bc.Val2()(0,0);    //    Ux displacement
            v[1] = bc.Val2()(1,0);    //    Uz displacement
            
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                ef(3*in+0,0)	+= gBigNumber*(u[0] - v[0])*phiu(in,0)*weight;	// X displacement Value
                ef(3*in+2,0)	+= gBigNumber*(u[2] - v[1])*phiu(in,0)*weight;	// Z displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in+0,3*jn+0)	+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
                    ek(3*in+2,3*jn+2)	+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// Z displacement
                }
            }
            
            break;
            
        }
            
        case 10 : // Duyz
            // Dirichlet in y and z direction of displacement
        {
            REAL v[2];
            v[0] = bc.Val2()(0,0);    //    Uy displacement
            v[1] = bc.Val2()(1,0);    //    Uz displacement
            
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                ef(3*in+1,0)	+= gBigNumber*(u[1] - v[0])*phiu(in,0)*weight;	// Y displacement Value
                ef(3*in+2,0)	+= gBigNumber*(u[2] - v[1])*phiu(in,0)*weight;	// Z displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in+1,3*jn+1)	+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// Y displacement
                    ek(3*in+2,3*jn+2)	+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// Z displacement
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
void TPZPMRSCouplPoroElast::FillDataRequirements(TPZVec<TPZMaterialData > &datavec)

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

void TPZPMRSCouplPoroElast::FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec)
{
    int nref = datavec.size();
    for(int i = 0; i<nref; i++)
    {
        datavec[i].fNeedsSol = true;
        datavec[i].fNeedsNormal = true;
        datavec[i].fNeedsNeighborSol = true;
    }
}

void TPZPMRSCouplPoroElast::Print(std::ostream &out)
{
    out << "Material Name : " << Name() << "\n";
    out << "Properties for TPZPMRSCouplPoroElast: \n";
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
int TPZPMRSCouplPoroElast::VariableIndex(const std::string &name)
{
    //	Total Strain Variables
    if(!strcmp("et_v",name.c_str()))             return	18;
    if(!strcmp("et_x",name.c_str()))             return	19;
    if(!strcmp("et_y",name.c_str()))             return	20;
    if(!strcmp("et_z",name.c_str()))             return	21;
    if(!strcmp("et_xy",name.c_str()))            return	22;
    if(!strcmp("et_xz",name.c_str()))            return	23;
    if(!strcmp("et_yz",name.c_str()))            return	24;
    
    //	Elastic Strain Variables
    if(!strcmp("e_v",name.c_str()))             return	25;
    if(!strcmp("e_x",name.c_str()))             return	26;
    if(!strcmp("e_y",name.c_str()))             return	27;
    if(!strcmp("e_z",name.c_str()))             return	28;
    if(!strcmp("e_xy",name.c_str()))            return	29;
    if(!strcmp("e_xz",name.c_str()))            return	30;
    if(!strcmp("e_yz",name.c_str()))            return	31;
    
    //	Plastic Strain Variables
    if(!strcmp("ep_v",name.c_str()))            return	32;
    if(!strcmp("ep_x",name.c_str()))            return	33;
    if(!strcmp("ep_y",name.c_str()))            return	34;
    if(!strcmp("ep_z",name.c_str()))            return	35;
    if(!strcmp("ep_xy",name.c_str()))           return	36;
    if(!strcmp("ep_xz",name.c_str()))           return	37;
    if(!strcmp("ep_yz",name.c_str()))           return	38;
    
    //	Displacement Variables
    if(!strcmp("u",name.c_str()))				return	39;
    
    //	Diffusion Variables
    if(!strcmp("p",name.c_str()))				return	40;
    if(!strcmp("v",name.c_str()))				return	41;
    if(!strcmp("phi",name.c_str()))				return	42;
    if(!strcmp("k_x",name.c_str()))				return	43;
    if(!strcmp("k_y",name.c_str()))				return	44;
    if(!strcmp("k_z",name.c_str()))				return	45;
    
    //	Total Stress Variables
    if(!strcmp("s_x",name.c_str()))             return	46;
    if(!strcmp("s_y",name.c_str()))             return	47;
    if(!strcmp("s_z",name.c_str()))             return	48;
    if(!strcmp("t_xy",name.c_str()))            return	49;
    if(!strcmp("t_xz",name.c_str()))            return	50;
    if(!strcmp("t_yz",name.c_str()))            return	51;
    
    //	Stress Ratio Variable
    if(!strcmp("K_0",name.c_str()))             return	52;
    
    //	Yield Surface Variable
    if(!strcmp("YS_1",name.c_str()))             return	53;
    if(!strcmp("YS_2",name.c_str()))             return	54;
    if(!strcmp("YS_3",name.c_str()))             return	55;
    
    return TPZMaterial::VariableIndex(name);
}

int TPZPMRSCouplPoroElast::NSolutionVariables(int var)
{
    if(var == 18)	return 1;
    if(var == 19)	return 1;
    if(var == 20)	return 1;
    if(var == 21)	return 1;
    if(var == 22)	return 1;
    if(var == 23)	return 1;
    if(var == 24)	return 1;
    if(var == 25)	return 1;
    if(var == 26)	return 1;
    if(var == 27)	return 1;
    if(var == 28)	return 1;
    if(var == 29)	return 1;
    if(var == 30)	return 1;
    if(var == 31)	return 1;
    if(var == 32)	return 1;
    if(var == 33)	return 1;
    if(var == 34)	return 1;
    if(var == 35)	return 1;
    if(var == 36)	return 1;
    if(var == 37)	return 1;
    if(var == 38)	return 1;
    if(var == 39)	return m_Dim;
    if(var == 40)	return 1;
    if(var == 41)	return m_Dim;
    if(var == 42)	return 1;
    if(var == 43)	return 1;
    if(var == 44)	return 1;
    if(var == 45)	return 1;
    if(var == 46)	return 1;
    if(var == 47)	return 1;
    if(var == 48)	return 1;
    if(var == 49)	return 1;
    if(var == 50)	return 1;
    if(var == 51)	return 1;
    if(var == 52)	return 1;
    if(var == 53)	return 1;
    if(var == 54)	return 1;
    if(var == 55)	return 1;
    
    return TPZMaterial::NSolutionVariables(var);
}

//	Calculate Secondary variables based on ux, uy, Pore pressure and their derivatives
void TPZPMRSCouplPoroElast::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout)
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
    
    
    REAL to_Mpa     = 1; // 1.0e-6;
    REAL to_Darcy   = 1; // 1.01327e+12; // 1.013249966e+12;
    
    
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
    
    
    Compute_Sigma_n(Grad_u_n, Grad_u, e_e, e_p, S);
    
    // ************************************** The value of parameters ************************
    
    // ************************	Total Strain Variables ************************
    //	epsilon_t_v
    if(var == 18)
    {
        Solout[0] = e_e(0,0) + e_e(1,1) + e_e(2,2);
        return;
    }
    
    //	epsilon_t_xx
    if(var == 19)
    {
        Solout[0] = e_e(0,0);
        return;
    }
    
    //	epsilon_t_yy
    if(var == 20)
    {
        Solout[0] = e_e(1,1);
        return;
    }
    
    //	epsilon_t_zz
    if(var == 21)
    {
        Solout[0] = e_e(2,2);
        return;
    }
    
    //	epsilon_t_xy
    if(var == 22)
    {
        Solout[0] = e_e(0,1);
        return;
    }
    
    //	epsilon_t_xz
    if(var == 23)
    {
        Solout[0] = e_e(0,2);
        return;
    }
    
    //	epsilon_t_yz
    if(var == 24)
    {
        Solout[0] = e_e(1,2);
        return;
    }
    
    // ************************	Elastic Strain Variables ************************
    //	epsilon_e_v
    if(var == 25)
    {
        Solout[0] = e_e(0,0) + e_e(1,1) + e_e(2,2);
        return;
    }
    
    //	epsilon_e_xx
    if(var == 26)
    {
        Solout[0] = e_e(0,0);
        return;
    }
    
    //	epsilon_e_yy
    if(var == 27)
    {
        Solout[0] = e_e(1,1);
        return;
    }
    
    //	epsilon_e_zz
    if(var == 28)
    {
        Solout[0] = e_e(2,2);
        return;
    }
    
    //	epsilon_e_xy
    if(var == 29)
    {
        Solout[0] = e_e(0,1);
        return;
    }
    
    //	epsilon_e_xz
    if(var == 30)
    {
        Solout[0] = e_e(0,2);
        return;
    }
    
    //	epsilon_e_yz
    if(var == 31)
    {
        Solout[0] = e_e(1,2);
        return;
    }
    
    // ************************	Plastic Strain Variables ************************
    //	epsilon_p_v
    if(var == 32)
    {
        Solout[0] = e_p(0,0) + e_p(1,1) + e_p(2,2);
        return;
    }
    
    //	epsilon_p_xx
    if(var == 33)
    {
        Solout[0] = e_p(0,0);
        return;
    }
    
    //	epsilon_p_yy
    if(var == 34)
    {
        Solout[0] = e_p(1,1);
        return;
    }
    
    //	epsilon_p_zz
    if(var == 35)
    {
        Solout[0] = e_p(2,2);
        return;
    }
    
    //	epsilon_p_xy
    if(var == 36)
    {
        Solout[0] = e_p(0,1);
        return;
    }
    
    //	epsilon_p_xz
    if(var == 37)
    {
        Solout[0] = e_p(0,2);
        return;
    }
    
    //	epsilon_p_yz
    if(var == 38)
    {
        Solout[0] = e_p(1,2);
        return;
    }
    
    // ************************	Displacement Variables ************************
    //  Displacement Variable
    if(var == 39)
    {
        Solout[0] = u[0];
        Solout[1] = u[1];
        if (m_Dim == 3)
        {
            Solout[2] = u[2];
        }
        return;
    }
    
    
    // ************************	Diffusion Variables ************************
    //	pore pressure
    if(var == 40)
    {
        Solout[0] = p[0]*to_Mpa;
        return;
    }
    
    //	Darcy's velocity
    if(var == 41)
    {
        if (m_Dim != 3)
        {
            Grad_p(0,0) = dp(0,0)*axes_p(0,0)+dp(1,0)*axes_p(1,0); // dp/dx
            Grad_p(1,0) = dp(0,0)*axes_p(0,1)+dp(1,0)*axes_p(1,1); // dp/dy
            
            REAL phi = porosity_corrected_2D(datavec);
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
    
    
    //	Porosity
    if(var == 42)
    {
        if (m_Dim != 3)
        {
            Solout[0] = porosity_corrected_2D(datavec);
            return;
        }
        else
        {
            Solout[0] = porosity_corrected_3D(datavec);
            return;
        }
    }
    
    
    //	k_x
    if(var == 43)
    {
        if (m_Dim != 3)
        {
        REAL phi = porosity_corrected_2D(datavec);
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
    if(var == 44)
    {
        if (m_Dim != 3)
        {
            REAL phi = porosity_corrected_2D(datavec);
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
    if(var == 45)
    {
        if (m_Dim != 3)
        {
            REAL phi = porosity_corrected_2D(datavec);
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
    
    
    // ************************	Total Stress Variables ************************
    //	sigma_x
    if(var == 46)
    {
        Solout[0] = S(0,0)*to_Mpa;
        return;
    }
    
    //	sigma_y
    if(var == 47)
    {
        Solout[0] = S(1,1)*to_Mpa;
        return;
    }
    
    //	sigma_z
    if(var == 48)
    {
        Solout[0] = S(2,2)*to_Mpa;
        return;
    }
    
    //	tau_xy
    if(var == 49)
    {
        Solout[0] = S(0,1)*to_Mpa;
        return;
    }
    
    //	tau_xz
    if(var == 50)
    {
        Solout[0] = S(0,2)*to_Mpa;
        return;
    }
    
    //	tau_yz
    if(var == 51)
    {
        Solout[0] = S(1,2)*to_Mpa;
        return;
    }
    
    // ************************	Stress Ratio Variable ************************
    //	K_0
    if(var == 52)
    {
        Solout[0] = S(0,0)/S(1,1);
        return;
    }
    
    // ************************	Yield Surface Variable ************************
    //  YS_1
    if(var == 53)
    {
       Solout[0] = 0;
        return;
    }
    
    //  YS_2
    if(var == 54)
    {
        Solout[0] = 0;
        return;
    }
    
    //  YS_3
    if(var == 55)
    {
        Solout[0] = 0;
        return;
    }
    
}



/** @brief computation of effective sigma */
void TPZPMRSCouplPoroElast::Compute_Sigma_n(TPZFMatrix<REAL> Grad_u_n, TPZFMatrix<REAL> Grad_u, TPZFMatrix<REAL> &e_e, TPZFMatrix<REAL> &e_p, TPZFMatrix<REAL> &S)
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
    
    Grad_u_n = Grad_u;
    Grad_du = Grad_u_n; // Linear case
    Grad_du.Transpose(&Grad_du_Transpose);
    delta_e = Grad_du + Grad_du_Transpose;
    delta_e *= 0.5;
    
    TPZFNMatrix<9,REAL> e_t, e_tn;
    TPZFNMatrix<9,REAL> S_tn,s_tn, I(delta_e.Rows(),delta_e.Cols(),0.0);
    I.Identity();
    
    /** Total strain */
    e_t = e_e + e_p;
    e_tn = e_t + delta_e;

    /** Total stress */
    REAL trace = (e_tn(0,0) + e_tn(1,1) + e_tn(2,2));
    s_tn = 2.0 * m_mu * e_tn + m_lambda * trace * I;
    
    // Convert to principal stresses
    Principal_Stress(s_tn, S_tn);
    
    /** Update the parameters */
    e_e = e_tn;
    S = s_tn;
    
    return;
    
}

/** @brief Principal Stress */
void TPZPMRSCouplPoroElast::Principal_Stress(TPZFMatrix<REAL> T, TPZFMatrix<REAL> & S)
{
    
#ifdef PZDEBUG
    
    if (T.Rows() != 3 && T.Cols() != 3)
    {
        DebugStop();
    }
    
#endif
    
    T += 1.0e-18;
    
    REAL a,b,c,d;
    a = 1.0;
    b = - T(0,0) - T(1,1) - T(2,2);
    c = - T(0,1)*T(0,1) - 2.0*T(0,2)*T(0,2) + T(0,0) * T(1,1) + T(0,0) * T(2,2) + T(1,1) * T(2,2);
    d = T(0,0) * T(0,2)*T(0,2) - 2.0* T(0,1) * T(0,2)*T(0,2) + T(0,2)*T(0,2) * T(1,1) + T(0,1)*T(0,1)*T(2,2) - T(0,0) * T(1,1) * T(2,2);
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
    
    // Sorting of Stress
    REAL s1 = std::max(r[0], std::max(r[1], r[2]));
    REAL s3 = std::min(r[0], std::min(r[1], r[2]));
    REAL s2 = 0.0;
    for (int i = 0; i < 3 ; i++)
    {
        if(fabs(r[i]  - s1) <= 1.0e-10 || fabs(r[i] - s3) <= 1.0e-10)
        {
            continue;
        }
        s2 = r[i];
    }
    
    S.Resize(3, 3);
    S.Zero();
    S(0,0) = s1;
    S(1,1) = s2;
    S(2,2) = s3;
    
}

