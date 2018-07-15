//
//  TPZPMRSCoupling.cpp
//  PZ
//
//  Created by Omar and Manouchehr on 8/28/16.
//
//

#include "TPZPMRSCoupling.h"
#include <iostream>
#include <string>
#include "pzbndcond.h"
#include "pzaxestools.h"
#include <algorithm>

#include "pzfmatrix.h"
#include "TPZTensor.h"

// Elasticity
#include "TPZElasticCriterion.h"

// Plasticity
#include "TPZPlasticStepPV.h"
#include "TPZSandlerExtended.h"
#include "TPZYCMohrCoulombPV.h"
#include "TPZMatElastoPlastic.h"
#include "TPZElasticCriterion.h"
#include "TPZSandlerDimaggio.h"



/** @brief default costructor */
TPZPMRSCoupling::TPZPMRSCoupling():TPZMatWithMem<TPZPMRSMemory,TPZDiscontinuousGalerkin>(), m_nu(0.), m_alpha(0.), m_k_0(0.), m_eta(0.), m_PlaneStress(0)
{

    m_Dim = 3;
    m_b.resize(3);
    
    m_b[0]=0.;
    m_b[1]=0.;
    m_b[2]=0.;
    
    m_PlaneStress = 1.;
    
    m_rho_s = 2700.0; 
    m_rho_f = 1000.0;
    
    m_k_model = 0;
    m_eta_dp = 0.0;
    m_xi_dp = 0.0;
    
}

/** @brief costructor based on a material id */
TPZPMRSCoupling::TPZPMRSCoupling(int matid, int dim):TPZMatWithMem<TPZPMRSMemory,TPZDiscontinuousGalerkin>(matid), m_nu(0.), m_alpha(0.), m_k_0(0.), m_eta(0.),m_PlaneStress(0)
{

    m_Dim = dim;
    m_b.resize(3);
    
    m_b[0]=0.;
    m_b[1]=0.;
    m_b[2]=0.;
    
    m_PlaneStress = 1;
    
    m_rho_s = 2700.0;
    m_rho_f = 1000.0;
    
    m_k_model = 0;
    m_eta_dp = 0.0;
    m_xi_dp = 0.0;
    
}

/** @brief default destructor */
TPZPMRSCoupling::~TPZPMRSCoupling()
{
}


/** @brief copy constructor $ */
TPZPMRSCoupling::TPZPMRSCoupling(const TPZPMRSCoupling& other)
{
    this->m_Dim    = other.m_Dim;
    this->m_SimulationData    = other.m_SimulationData;
}


/** @brief Copy assignemnt operator $ */
TPZPMRSCoupling& TPZPMRSCoupling::operator = (const TPZPMRSCoupling& other)
{
    
    if (this != & other) // prevent self-assignment
    {
        this->m_Dim    = other.m_Dim;
        this->m_SimulationData    = other.m_SimulationData;
    }
    return *this;
}


/** @brief number of state variables */
int TPZPMRSCoupling::NStateVariables()
{
    return 1;
}

/** @brief permeability coupling models  */
REAL TPZPMRSCoupling::k_permeability(REAL &phi, REAL &k)
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
REAL TPZPMRSCoupling::porosity_corrected(TPZVec<TPZMaterialData> &datavec)
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
REAL TPZPMRSCoupling::porosity_corrected_3D(TPZVec<TPZMaterialData> &datavec)
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



/** @brief computation of effective sigma 2D */
void TPZPMRSCoupling::Compute_Sigma(TPZFMatrix<REAL> & S_eff,TPZFMatrix<REAL> & Grad_u, REAL p_ex)
{
    
    TPZFNMatrix<6,REAL> Grad_ut(2,2,0.0), epsilon(2,2,0.0), I(2,2,0.0);
    Grad_u.Transpose(&Grad_ut);
    
    epsilon = Grad_u + Grad_ut;
    epsilon *= 0.5;
    
    I(0,0) = 1.0;
    I(1,1) = 1.0;
    
    REAL trace = (epsilon(0,0) + epsilon(1,1));
    
    S_eff = 2.0 * m_mu * epsilon + m_lambda * trace * I - 0.0 * m_alpha * p_ex * I;
    
}

/** @brief computation of sigma 2D */
void TPZPMRSCoupling::Compute_Sigma(TPZFMatrix<REAL> & S,TPZFMatrix<REAL> & Grad_v)
{
    
    TPZFNMatrix<6,REAL> Grad_vt(3,3,0.0), epsilon(3,3,0.0), I(3,3,0.0);
    Grad_v.Transpose(&Grad_vt);
    
    epsilon = Grad_v + Grad_vt;
    epsilon *= 0.5;
    
    I.Identity();
    
    REAL trace = (epsilon(0,0) + epsilon(1,1));
    
    S = 2.0 * m_mu * epsilon + m_lambda * trace * I;
    
}


/** @brief Compute effective stress 3D */
void TPZPMRSCoupling::Compute_Sigma(REAL & l, REAL & mu, TPZFMatrix<REAL> & S,TPZFMatrix<REAL> & Grad_u)
{
    
    
    REAL trace;
    for (int i = 0; i < 3; i++)
    {
        trace = 0.0;
        for (int j = 0; j < 3; j++)
        {
            S(i,j) = mu * (Grad_u(i,j) + Grad_u(j,i));
            trace +=  Grad_u(j,j);
        }
        S(i,i) += l * trace;
    }
    
    return;
}



/** @brief of inner product in 3D */
REAL TPZPMRSCoupling::Inner_Product(TPZFMatrix<REAL> & S,TPZFMatrix<REAL> & T)

{

    REAL inner_product = S(0,0) * T(0,0) + S(0,1) * T(0,1) + S(0,2) * T(0,2) + S(1,0) * T(1,0) + S(1,1) * T(1,1) + S(1,2) * T(1,2) + S(2,0) * T(2,0) + S(2,1) * T(2,1) + S(2,2) * T(2,2); //     S11 T11 + S12 T12 + S13 T13 + S21 T21 + S22 T22 + S23 T23 + S31 T31 + S32 T32 + S33 T33

    return inner_product;
}


/** @brief of contribute of BC */
void TPZPMRSCoupling::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE>  &ek, TPZFMatrix<STATE> &ef)
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
void TPZPMRSCoupling::Contribute_2D(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE>  &ek, TPZFMatrix<STATE> &ef)
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
    TPZFNMatrix<9,REAL> Grad_u(3,3,0.0),Grad_u_n,e_e,e_p,S;
    Grad_u(0,0) = du(0,0)*axes_u(0,0)+du(1,0)*axes_u(1,0); // dux/dx
    Grad_u(0,1) = du(0,0)*axes_u(0,1)+du(1,0)*axes_u(1,1); // dux/dy
    
    Grad_u(1,0) = du(0,1)*axes_u(0,0)+du(1,1)*axes_u(1,0); // duy/dx
    Grad_u(1,1) = du(0,1)*axes_u(0,1)+du(1,1)*axes_u(1,1); // duy/dy
    
    // Get the solution at the integrations points
    long global_point_index = datavec[0].intGlobPtIndex;
    TPZPMRSMemory &point_memory = GetMemory()[global_point_index];
    e_e = point_memory.epsilon_e_n();
    e_p = point_memory.epsilon_p_n();
    Grad_u_n = point_memory.grad_u_n();
    
    corrector_DP(Grad_u_n, Grad_u, e_e, e_p, S);
    
    TPZFNMatrix<6,REAL> Grad_vx_i(2,1,0.0),Si_x;
    TPZFNMatrix<6,REAL> Grad_vy_i(2,1,0.0),Si_y;

    TPZFNMatrix<6,REAL> Grad_v(2,2,0.0),T(2,2,0.0);
    TPZFNMatrix<6,REAL> Grad_vx_j(2,1,0.0),Tj_x;
    TPZFNMatrix<6,REAL> Grad_vy_j(2,1,0.0),Tj_y;

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
    m_k_model = 0;
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
    
}


// Contribute Methods being used

void TPZPMRSCoupling::Contribute_3D(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
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
    TPZPMRSMemory &point_memory = GetMemory()[global_point_index];
    e_e = point_memory.epsilon_e_n();
    e_p = point_memory.epsilon_p_n();
    Grad_u_n = point_memory.grad_u_n();
    
    corrector_DP(Grad_u_n, Grad_u, e_e, e_p, S);
    
    TPZFNMatrix<9,REAL> Grad_vx_i(3,1,0.0),Si_x;
    TPZFNMatrix<9,REAL> Grad_vy_i(3,1,0.0),Si_y;
    TPZFNMatrix<9,REAL> Grad_vz_i(3,1,0.0),Si_z;
    
    TPZFNMatrix<9,REAL> Grad_v(3,3,0.0),T(3,3,0.0);
    TPZFNMatrix<9,REAL> Grad_vx_j(3,1,0.0),Tj_x;
    TPZFNMatrix<9,REAL> Grad_vy_j(3,1,0.0),Tj_y;
    TPZFNMatrix<9,REAL> Grad_vz_j(3,1,0.0),Tj_z;
    
    REAL dvxdx, dvxdy, dvxdz;
    REAL dvydx, dvydy, dvydz;
    REAL dvzdx, dvzdy, dvzdz;
    
    REAL duxdx, duxdy, duxdz;
    REAL duydx, duydy, duydz;
    REAL duzdx, duzdy, duzdz;
    
    TPZFMatrix<REAL> & S_0 = m_SimulationData->Sigma_0();
    
    S_0.Zero();
    S -= S_0; // Applying prestress
    
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
void TPZPMRSCoupling::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef){
    
    TPZFMatrix<STATE>  ek_fake(ef.Rows(),ef.Rows(),0.0);
    this->Contribute(datavec, weight, ek_fake, ef);
    
}

/** @brief of contribute of BC */
void TPZPMRSCoupling::ContributeBC(TPZVec<TPZMaterialData> &datavec,REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
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
void TPZPMRSCoupling::ContributeBC_2D(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
    
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
    
    
    REAL v[3];
    v[0] = bc.Val2()(0,0);	//	Ux displacement
    v[1] = bc.Val2()(1,0);	//	Uy displacement
    v[2] = bc.Val2()(2,0);	//	Pressure
    
    REAL time = this->SimulationData()->t();
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
                ef(2*in  ,0)      += gBigNumber*(u[0] - v[0])*phiu(in,0)*weight;    // X displacement Value
                ef(2*in+1,0)      += gBigNumber*(u[1] - v[1])*phiu(in,0)*weight;    // y displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(2*in,2*jn    )    += gBigNumber*phiu(in,0)*phiu(jn,0)*weight;    // X displacement
                    ek(2*in+1,2*jn+1)    += gBigNumber*phiu(in,0)*phiu(jn,0)*weight;    // Y displacement
                }
            }
            
            //    Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //    Contribution for load Vector
                ef(in+2*phru,0)        += gBigNumber*(p[0]-v[2])*phip(in,0)*weight;    // P Pressure
                
                for (jn = 0 ; jn < phrp; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(in+2*phru,jn+2*phru)        += gBigNumber*phip(in,0)*phip(jn,0)*weight;    // P Pressure
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
                ef(2*in,0)        += gBigNumber*(u[0]-v[0])*phiu(in,0)*weight;    // X displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(2*in,2*jn)        += gBigNumber*phiu(in,0)*phiu(jn,0)*weight;    // X displacement
                }
            }
            
            //    Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //    Contribution for load Vector
                ef(in+2*phru,0)        += gBigNumber*(p[0]-v[2])*phip(in,0)*weight;    // P Pressure
                
                for (jn = 0 ; jn < phrp; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(in+2*phru,jn+2*phru)        += gBigNumber*phip(in,0)*phip(jn,0)*weight;    // P Pressure
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
                ef(2*in+1,0)      += gBigNumber*(u[1]-v[1])*phiu(in,0)*weight;    // y displacement
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(2*in+1,2*jn+1)    += gBigNumber*phiu(in,0)*phiu(jn,0)*weight;    // Y displacement
                }
            }
            
            //    Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //    Contribution for load Vector
                ef(in+2*phru,0)        += gBigNumber*(p[0]-v[2])*phip(in,0)*weight;    // P Pressure
                
                for (jn = 0 ; jn < phrp; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(in+2*phru,jn+2*phru)        += gBigNumber*phip(in,0)*phip(jn,0)*weight;    // P Pressure
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
                ef(in+2*phru,0)        += gBigNumber*(p[0]-v[2])*phip(in,0)*weight;    // P Pressure
                
                for (jn = 0 ; jn < phrp; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(in+2*phru,jn+2*phru)        += gBigNumber*phip(in,0)*phip(jn,0)*weight;    // P Pressure
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
                ef(in+2*phru,0)        += gBigNumber*(p[0]-v[2])*phip(in,0)*weight;    // P Pressure
                
                for (jn = 0 ; jn < phrp; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(in+2*phru,jn+2*phru)        += gBigNumber*phip(in,0)*phip(jn,0)*weight;    // P Pressure
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
                ef(2*in,0)        += gBigNumber*(u[0] - v[0])*phiu(in,0)*weight;    // X displacement Value
                ef(2*in+1,0)      += gBigNumber*(u[1] - v[1])*phiu(in,0)*weight;    // y displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(2*in,2*jn)        += gBigNumber*phiu(in,0)*phiu(jn,0)*weight;    // X displacement
                    ek(2*in+1,2*jn+1)    += gBigNumber*phiu(in,0)*phiu(jn,0)*weight;    // Y displacement
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
                ef(2*in,0)        += gBigNumber*(u[0]-v[0])*phiu(in,0)*weight;    // X displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(2*in,2*jn)        += gBigNumber*phiu(in,0)*phiu(jn,0)*weight;    // X displacement
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
                ef(2*in+1,0)      += gBigNumber*(u[1]-v[1])*phiu(in,0)*weight;    // y displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(2*in+1,2*jn+1)    += gBigNumber*phiu(in,0)*phiu(jn,0)*weight;    // Y displacement
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
                ef(2*in+1,0)      += gBigNumber*(u[1]-v[1])*phiu(in,0)*weight;    // y displacement
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(2*in+1,2*jn+1)    += gBigNumber*phiu(in,0)*phiu(jn,0)*weight;    // Y displacement
                }
            }
            
            //    Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //    Contribution for load Vector
                ef(in+2*phru,0)        += gBigNumber*(p[0]-v[2])*phip(in,0)*weight;    // P Pressure
                
                for (jn = 0 ; jn < phrp; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(in+2*phru,jn+2*phru)        += gBigNumber*phip(in,0)*phip(jn,0)*weight;    // P Pressure
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
                ef(in+2*phru,0)        += gBigNumber*(p[0]-v[2])*phip(in,0)*weight;    // P Pressure
                
                for (jn = 0 ; jn < phrp; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(in+2*phru,jn+2*phru)        += gBigNumber*phip(in,0)*phip(jn,0)*weight;    // P Pressure
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


void TPZPMRSCoupling::ContributeBC_3D(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
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
                ef(3*in  ,0)	+= gBigNumber*(u[0] - v[0])*phiu(in,0)*weight;	// X displacement Value
                ef(3*in+1,0)	+= gBigNumber*(u[1] - v[1])*phiu(in,0)*weight;	// y displacement Value
                ef(3*in+2,0)	+= gBigNumber*(u[2] - v[2])*phiu(in,0)*weight;	// Z displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in  ,3*jn  )	+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
                    ek(3*in+1,3*jn+1)	+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// Y displacement
                    ek(3*in+2,3*jn+2)	+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// Z displacement
                }
            }
            
            //    Diffusion Equation
            REAL p_s = p[0];
            REAL d_p = (p_s-v[3]);
            for(in = 0 ; in < phrp; in++)
            {
                //    Contribution for load Vector
                ef(in+3*phru,0)        += gBigNumber*(d_p)*phip(in,0)*weight;    // P Pressure
                
                for (jn = 0 ; jn < phrp; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(in+3*phru,jn+3*phru)        += gBigNumber*phip(in,0)*phip(jn,0)*weight;    // P Pressure
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
                ef(3*in  ,0)	+= gBigNumber*(u[0] - v[0])*phiu(in,0)*weight;	// X displacement Value
                ef(3*in+1,0)	+= gBigNumber*(u[1] - v[1])*phiu(in,0)*weight;	// y displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in  ,3*jn  )	+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
                    ek(3*in+1,3*jn+1)	+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// Y displacement
                }
            }
            //    Diffusion Equation
            REAL p_s = p[0];
            REAL d_p = (p_s-v[2]);
            for(in = 0 ; in < phrp; in++)
            {
                //    Contribution for load Vector
                ef(in+3*phru,0)        += gBigNumber*(d_p)*phip(in,0)*weight;    // P Pressure
                
                for (jn = 0 ; jn < phrp; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(in+3*phru,jn+3*phru)        += gBigNumber*phip(in,0)*phip(jn,0)*weight;    // P Pressure
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
                ef(3*in  ,0)	+= gBigNumber*(u[0] - v[0])*phiu(in,0)*weight;	// X displacement Value
                ef(3*in+2,0)	+= gBigNumber*(u[2] - v[1])*phiu(in,0)*weight;	// Z displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in  ,3*jn  )	+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
                    ek(3*in+2,3*jn+2)	+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// Z displacement
                }
            }
            //    Diffusion Equation
            REAL p_s = p[0];
            REAL d_p = (p_s-v[2]);
            for(in = 0 ; in < phrp; in++)
            {
                //    Contribution for load Vector
                ef(in+3*phru,0)        += gBigNumber*(d_p)*phip(in,0)*weight;    // P Pressure
                
                for (jn = 0 ; jn < phrp; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(in+3*phru,jn+3*phru)        += gBigNumber*phip(in,0)*phip(jn,0)*weight;    // P Pressure
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
                ef(3*in+1,0)	+= gBigNumber*(u[1] - v[0])*phiu(in,0)*weight;	// Y displacement Value
                ef(3*in+2,0)	+= gBigNumber*(u[2] - v[1])*phiu(in,0)*weight;	// Z displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in+1,3*jn+1)	+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// Y displacement
                    ek(3*in+2,3*jn+2)	+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// Z displacement
                }
            }
            //    Diffusion Equation
            REAL p_s = p[0];
            REAL d_p = (p_s-v[2]);
            for(in = 0 ; in < phrp; in++)
            {
                //    Contribution for load Vector
                ef(in+3*phru,0)        += gBigNumber*(d_p)*phip(in,0)*weight;    // P Pressure
                
                for (jn = 0 ; jn < phrp; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(in+3*phru,jn+3*phru)        += gBigNumber*phip(in,0)*phip(jn,0)*weight;    // P Pressure
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
                ef(3*in,0)		+= gBigNumber*(u[0] - v[0])*phiu(in,0)*weight;	// X displacement Value
                
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in,3*jn)		+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
                    
                }
            }
            //    Diffusion Equation
            REAL p_s = p[0];
            REAL d_p = (p_s-v[1]);
            for(in = 0 ; in < phrp; in++)
            {
                //    Contribution for load Vector
                ef(in+3*phru,0)        += gBigNumber*(d_p)*phip(in,0)*weight;    // P Pressure
                
                for (jn = 0 ; jn < phrp; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(in+3*phru,jn+3*phru)        += gBigNumber*phip(in,0)*phip(jn,0)*weight;    // P Pressure
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
                ef(3*in+1,0)		+= gBigNumber*(u[1] - v[0])*phiu(in,0)*weight;	// Y displacement Value
                
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in+1,3*jn+1)		+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// Y displacement
                    
                }
            }
            //    Diffusion Equation
            REAL p_s = p[0];
            REAL d_p = (p_s-v[1]);
            for(in = 0 ; in < phrp; in++)
            {
                //    Contribution for load Vector
                ef(in+3*phru,0)        += gBigNumber*(d_p)*phip(in,0)*weight;    // P Pressure
                
                for (jn = 0 ; jn < phrp; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(in+3*phru,jn+3*phru)        += gBigNumber*phip(in,0)*phip(jn,0)*weight;    // P Pressure
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
                ef(3*in+2,0)		+= gBigNumber*(u[2] - v[0])*phiu(in,0)*weight;	// Z displacement Value
                
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in+2,3*jn+2)		+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// Z displacement
                    
                }
            }
            //    Diffusion Equation
            REAL p_s = p[0];
            REAL d_p = (p_s-v[1]);
            for(in = 0 ; in < phrp; in++)
            {
                //    Contribution for load Vector
                ef(in+3*phru,0)        += gBigNumber*(d_p)*phip(in,0)*weight;    // P Pressure
                
                for (jn = 0 ; jn < phrp; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(in+3*phru,jn+3*phru)        += gBigNumber*phip(in,0)*phip(jn,0)*weight;    // P Pressure
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
                ef(in+3*phru,0)        += gBigNumber*(d_p)*phip(in,0)*weight;    // P Pressure
                
                for (jn = 0 ; jn < phrp; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(in+3*phru,jn+3*phru)        += gBigNumber*phip(in,0)*phip(jn,0)*weight;    // P Pressure
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
                ef(in+3*phru,0)        += gBigNumber*(d_p)*phip(in,0)*weight;    // P Pressure
                
                for (jn = 0 ; jn < phrp; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(in+3*phru,jn+3*phru)        += gBigNumber*phip(in,0)*phip(jn,0)*weight;    // P Pressure
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
                ef(3*in  ,0)	+= gBigNumber*(u[0] - v[0])*phiu(in,0)*weight;	// X displacement Value
                ef(3*in+1,0)	+= gBigNumber*(u[1] - v[1])*phiu(in,0)*weight;	// Y displacement Value
                ef(3*in+2,0)	+= gBigNumber*(u[2] - v[2])*phiu(in,0)*weight;	// Z displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in  ,3*jn  )	+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
                    ek(3*in+1,3*jn+1)	+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// Y displacement
                    ek(3*in+2,3*jn+2)	+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// Z displacement
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
                ef(3*in  ,0)	+= gBigNumber*(u[0] - v[0])*phiu(in,0)*weight;	// X displacement Value
                ef(3*in+1,0)	+= gBigNumber*(u[1] - v[1])*phiu(in,0)*weight;	// Y displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in  ,3*jn  )	+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
                    ek(3*in+1,3*jn+1)	+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// Y displacement
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
                ef(3*in  ,0)	+= gBigNumber*(u[0] - v[0])*phiu(in,0)*weight;	// X displacement Value
                ef(3*in+2,0)	+= gBigNumber*(u[2] - v[1])*phiu(in,0)*weight;	// Z displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in  ,3*jn  )	+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
                    ek(3*in+2,3*jn+2)	+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// Z displacement
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
                ef(3*in+1,0)	+= gBigNumber*(u[1] - v[0])*phiu(in,0)*weight;	// Y displacement Value
                ef(3*in+2,0)	+= gBigNumber*(u[2] - v[1])*phiu(in,0)*weight;	// Z displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in+1,3*jn+1)	+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// Y displacement
                    ek(3*in+2,3*jn+2)	+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// Z displacement
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
                ef(3*in,0)		+= gBigNumber*(u[0] - v[0])*phiu(in,0)*weight;	// X displacement Value
                
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in,3*jn)		+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
                    
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
                ef(3*in+1,0)		+= gBigNumber*(u[1] - v[0])*phiu(in,0)*weight;	// Y displacement Value
                
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in+1,3*jn+1)		+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// Y displacement
                    
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
                ef(3*in+2,0)		+= gBigNumber*(u[2] - v[0])*phiu(in,0)*weight;	// Z displacement Value
                
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in+2,3*jn+2)		+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// Z displacement
                    
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



void TPZPMRSCoupling::FillDataRequirements(TPZVec<TPZMaterialData > &datavec)

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

void TPZPMRSCoupling::FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec)
{
    int nref = datavec.size();
    for(int i = 0; i<nref; i++)
    {
        datavec[i].fNeedsSol = true;
        datavec[i].fNeedsNormal = true;
        datavec[i].fNeedsNeighborSol = true;
    }
}

void TPZPMRSCoupling::Print(std::ostream &out)
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
int TPZPMRSCoupling::VariableIndex(const std::string &name)
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

int TPZPMRSCoupling::NSolutionVariables(int var)
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
void TPZPMRSCoupling::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout)
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
    REAL to_Darcy   = 1; // 1.013249966e+12;
    
    
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
    
    
    corrector_DP(Grad_u_n, Grad_u, e_e, e_p, S);
    
    
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
//        Solout[0] = S(0,0)*to_Mpa-p[0]*to_Mpa;
        Solout[0] = S(0,0)*to_Mpa;
        return;
    }
    
    //	sigma_y
    if(var == 3)
    {
//        Solout[0] = S(1,1)*to_Mpa-p[0]*to_Mpa;
        Solout[0] = S(1,1)*to_Mpa;
        return;
    }
    
    //	sigma_z
    if(var == 4)
    {
//        Solout[0] = S(2,2)*to_Mpa-p[0]*to_Mpa;
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

/** @brief mean stress */
REAL TPZPMRSCoupling::p_m(TPZFMatrix<REAL> T)
{
    
#ifdef PZDEBUG
    if (T.Rows() != 3 && T.Cols() != 3)
    {
        DebugStop();
    }
#endif
    
    REAL mean_stress = 0.0;
    mean_stress = fabs((T(0,0) + T(1,1) + T(2,2)))/3.0;
    return mean_stress;
}

/** @brief deviatoric stress */
TPZFMatrix<REAL> TPZPMRSCoupling::s(TPZFMatrix<REAL> T)
{

#ifdef PZDEBUG
    if (T.Rows() != 3 && T.Cols() != 3)
    {
        DebugStop();
    }
#endif
    
    REAL  mean_stress = p_m(T);
    TPZFMatrix<REAL> H = T;
    TPZFNMatrix<9,REAL> I(3,3,0.0);
    I.Identity();
    H = T - mean_stress * I;
    return H;
}

/** @brief J2 invariant stress */
REAL TPZPMRSCoupling::J2(TPZFMatrix<REAL> T)
{
    
#ifdef PZDEBUG
    if (T.Rows() != 3 && T.Cols() != 3)
    {
        DebugStop();
    }
#endif
    TPZFMatrix<REAL> S = T;
    TPZFMatrix<REAL> S_inner = T;
    S.Transpose(&S);
    S.Multiply(T,S_inner);
    
    REAL j2 = 0.5*(S_inner(0,0) + S_inner(1,1) + S_inner(2,2));
    
    return j2
    ;
}

/** @brief J3 invariant stress */
REAL TPZPMRSCoupling::J3(TPZFMatrix<REAL> T)
{
    
#ifdef PZDEBUG
    if (T.Rows() != 3 && T.Cols() != 3)
    {
        DebugStop();
    }
#endif
    TPZFMatrix<REAL> Tinv = T;
    REAL det = T(0,0)*T(1,1)-T(1,0)*T(0,1);
    Tinv.Resize(3, 3);
    Tinv(2,2) = Tinv(1,1);
    Tinv.DeterminantInverse(det, Tinv);
    
    return det;
    
}

/** @brief Lode angle, theta */
REAL TPZPMRSCoupling::theta(TPZFMatrix<REAL> T)
{
    
#ifdef PZDEBUG
    if (T.Rows() != 3 && T.Cols() != 3)
    {
        DebugStop();
    }
#endif
    
    REAL theta;
    REAL arg = -3.0*sqrt(3.0)*J3(s(T))/(2.0*pow(J2(s(T)), 1.5)) + 1.0e-14;
    theta = (1.0/3.0)*asin(arg);
    
    return theta;
}

/** @brief Phi Mohr-Coulomb */
REAL TPZPMRSCoupling::Phi_MC(TPZFMatrix<REAL> T)
{

#ifdef PZDEBUG
    if (T.Rows() != 3 && T.Cols() != 3) {
        DebugStop();
    }
#endif
    
    REAL theta_v = theta(T);
    
    REAL phi = (cos(theta_v) - (1.0/sqrt(3.0))*sin(theta_v)*sin(m_phi_f)) * sqrt(J2(s(T))) + p_m(T)*sin(m_phi_f) - m_c * cos(m_phi_f);
    return phi;
}

/** @brief Phi Drucker-Prager */
REAL TPZPMRSCoupling::Phi_DP(TPZFMatrix<REAL> T)
{
    
#ifdef PZDEBUG
    if (T.Rows() != 3 && T.Cols() != 3)
    {
        DebugStop();
    }
#endif
    REAL eta = 6.0*(sin(m_phi_f))/(sqrt(3.0)*(3.0-sin(m_phi_f)));
    REAL xi = 6.0*(cos(m_phi_f))/(sqrt(3.0)*(3.0-sin(m_phi_f)));
    REAL phi = sqrt(J2(s(T))) + eta *  p_m(T) - xi * m_c ;
    return phi;

}

/** @brief plasticity multiplier delta_gamma */
REAL TPZPMRSCoupling::Phi_tilde_DP(TPZFMatrix<REAL> T, REAL d_gamma_guest)
{
    
#ifdef PZDEBUG
    if (T.Rows() != 3 && T.Cols() != 3)
    {
        DebugStop();
    }
#endif
    REAL phi = sqrt(J2(s(T))) - m_mu * d_gamma_guest + m_eta_dp *  (p_m(T) - m_K * m_eta_dp * d_gamma_guest) - m_xi_dp * m_c ;
    return phi;
    
}

/** @brief plasticity multiplier delta_gamma */
REAL TPZPMRSCoupling::Phi_tilde_DP_delta_gamma(TPZFMatrix<REAL> T, REAL d_gamma_guest)
{
    
#ifdef PZDEBUG
    if (T.Rows() != 3 && T.Cols() != 3)
    {
        DebugStop();
    }
#endif

    REAL d_phi_d_delta_gamma = - m_mu - m_K * m_eta_dp * m_eta_dp;
    
    if (fabs(d_phi_d_delta_gamma) <= 1.0e-18)
    {
        d_phi_d_delta_gamma = 1.0e-18;
    }
    
    return d_phi_d_delta_gamma;
    
}

/** @brief plasticity multiplier delta_gamma using newton iterations */
REAL TPZPMRSCoupling::delta_gamma_finder(TPZFMatrix<REAL> T, REAL d_gamma_guest)
{
    
    REAL tol = 1.0e-10;
    REAL error = 1.0;
    int n_iter = 20;
    REAL d_gamma_converged = d_gamma_guest;
    
    for (int i = 0; i < n_iter; i++) {
        d_gamma_converged = d_gamma_converged - Phi_tilde_DP(T,d_gamma_converged) / Phi_tilde_DP_delta_gamma(T,d_gamma_converged);
        error = Phi_tilde_DP(T,d_gamma_converged);
        if (error <= tol) {
            break;
        }
        
    }
    
    return d_gamma_converged;
    
}

/** @brief Drucker prager strain update */
TPZFMatrix<REAL> TPZPMRSCoupling::strain_DP(TPZFMatrix<REAL> T)
{
    DebugStop();
	return T;
}

/** @brief Drucker prager stress update */
TPZFMatrix<REAL> TPZPMRSCoupling::stress_DP(TPZFMatrix<REAL> T)
{
    DebugStop();
	return T;
}

/** @brief Drucker prager elastoplastic corrector  */
void TPZPMRSCoupling::corrector_DP(TPZFMatrix<REAL> Grad_u_n, TPZFMatrix<REAL> Grad_u, TPZFMatrix<REAL> &e_e, TPZFMatrix<REAL> &e_p, TPZFMatrix<REAL> &S)
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
    
    
    TPZFNMatrix<9,REAL> e_t, e_trial;
    TPZFNMatrix<9,REAL> S_trial,s_trial, I(delta_e.Rows(),delta_e.Cols(),0.0);
    I.Identity();
    
    /** Trial strain */
    e_t = e_e + e_p;
    e_trial = e_t + delta_e;

    /** Trial stress */
    REAL trace = (e_trial(0,0) + e_trial(1,1) + e_trial(2,2));
    s_trial = 2.0 * m_mu * e_trial + m_lambda * trace * I;
    
    // convert to principal stresses
    Principal_Stress(s_trial, S_trial);
    
    /** Elastic update */
    e_e = e_trial;
    S = s_trial;
    
    return;
    
    if (Phi_DP(s_trial) < 0.0)
    {
        /** Elastic update */
        e_e = e_trial;
        S = s_trial;
    }
    else{
        /** Plastic update */
        REAL delta_gamma = 0.0;
        delta_gamma = delta_gamma_finder(s_trial, delta_gamma);
        e_e = e_trial;
        e_p = delta_gamma * ( ( 1.0/(2.0*sqrt(J2(s(s_trial))) ) ) * s(s_trial) + (m_eta_dp/3.0)* I);
        S = s_trial - delta_gamma * ( ( m_mu/(2.0*sqrt(J2(s(s_trial))) ) ) * s(s_trial) + (m_K * m_eta_dp/3.0)* I);

//        e_e.Print("e_e = ");
//        e_p.Print("e_p = ");
//        s_trial.Print("s_trial = ");
//        S.Print("s = ");
    }
    
}

/** @brief Principal Stress */
void TPZPMRSCoupling::Principal_Stress(TPZFMatrix<REAL> T, TPZFMatrix<REAL> & S)
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
    
    // sorting
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

