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

#include "TPZMatElastoPlastic.h"


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
    
    m_VariableYoung = false;
    m_SetRunPlasticity = false;
    m_UpdateToUseFullDiplacement = false;
    
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

    m_VariableYoung = false;
    m_SetRunPlasticity = false;
    m_UpdateToUseFullDiplacement = false;

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

/** @brief a computation of strain */
template <class T, class TMEM>
void TPZPMRSCouplPoroPlast<T,TMEM>::ComputeStrainVector(TPZMaterialData & data, TPZFMatrix<REAL> &Strain)
{
    TPZMatElastoPlastic<T,TMEM>::ComputeStrainVector(data, Strain);
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
void TPZPMRSCouplPoroPlast<T,TMEM>::ComputeStressVector(TPZMaterialData & data, TPZFMatrix<REAL> &Stress)
{
    TPZMatElastoPlastic<T,TMEM>::ComputeStressVector(data,Stress);
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

/** @brief the 2D Poroelastic porosity correction from strains and pressure */
template<class T,class TMEM>
REAL TPZPMRSCouplPoroPlast<T,TMEM>::porosity_corrected_2D(TPZTensor<STATE> & eps_elastic, TPZTensor<STATE> & eps_plastic, STATE & pressure){
    REAL eps_volumetric = eps_elastic.I1() + eps_plastic.I1();
    REAL phi = m_porosity_0 + m_alpha * eps_volumetric + m_Se * pressure;
    return phi;
}


/** @brief the 3D Poroelastic porosity correction from strains and pressure */
template<class T,class TMEM>
REAL TPZPMRSCouplPoroPlast<T,TMEM>::porosity_corrected_3D(TPZTensor<STATE> & eps_elastic, TPZTensor<STATE> & eps_plastic, STATE & pressure){
    REAL eps_volumetric = eps_elastic.I1() + eps_plastic.I1();
    REAL phi = m_porosity_0 + m_alpha * eps_volumetric + m_Se * pressure;
    return phi;
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
    TPZFMatrix<REAL>    &phiu      =   datavec[u_b].phi;
    TPZFMatrix<REAL>    &phip      =   datavec[p_b].phi;
    
    TPZFMatrix<REAL> &grad_phi_u   =   datavec[u_b].dphix;
    TPZFMatrix<REAL> &grad_phi_p   =   datavec[p_b].dphix;
    
    TPZFNMatrix <9,REAL> &axes_u   =	datavec[u_b].axes;
    TPZFNMatrix <9,REAL> &axes_p   =	datavec[p_b].axes;
    
    // Getting the solutions and derivatives
    TPZManVector<REAL,3> u  = datavec[u_b].sol[0];
    TPZManVector<REAL,1> p  = datavec[p_b].sol[0];
    
    TPZFNMatrix <6,REAL> du = datavec[u_b].dsol[0];
    TPZFNMatrix <6,REAL> dp = datavec[p_b].dsol[0];
    
    TPZFNMatrix<6,REAL> Grad_p(2,1,0.0),Grad_phi_i(2,1,0.0),Grad_phi_j(2,1,0.0);
    Grad_p(0,0) = dp(0,0)*axes_p(0,0)+dp(1,0)*axes_p(1,0);
    Grad_p(1,0) = dp(0,0)*axes_p(0,1)+dp(1,0)*axes_p(1,1);
    
    int nphi_u = phiu.Rows();
    int nphi_p = phip.Rows();
    
    int first_u = 0;
    int first_p = 2*nphi_u;
    
    // Computing Gradient of the Solution
    TPZFNMatrix<6,REAL> Grad_vx_i(2,1,0.0);
    TPZFNMatrix<6,REAL> Grad_vy_i(2,1,0.0);
    
    TPZFNMatrix<6,REAL> Grad_v(2,2,0.0);
    TPZFNMatrix<6,REAL> Grad_vx_j(2,1,0.0);
    TPZFNMatrix<6,REAL> Grad_vy_j(2,1,0.0);
    
    
    // ElastoPlastic Strain Stress Parameters from material memory
    int global_point_index =  datavec[u_b].intGlobPtIndex;
    TPZAdmChunkVector<TMEM> memory_vec = *TPZMatWithMem<TMEM>::fMemory;
    TMEM &point_memory = memory_vec[global_point_index];
    
    // Setting for the elastoplastic integrator
    T elasto_plastic_integrator(this->m_plasticity_model);
    elasto_plastic_integrator.SetState(point_memory.fPlasticState);
    TPZElasticResponse ER;
    ER.SetUp(m_SimulationData->Get_young(), m_SimulationData->Get_nu());
    elasto_plastic_integrator.SetElasticResponse(ER);
    
    
    // obtaining the total strain
    TPZTensor<STATE> Stress ;//= point_memory.fSigma;
    TPZTensor<REAL> EpsT;
    TPZFNMatrix<6,STATE> deltastrain(6,1,0.);
    ComputeDeltaStrainVector(datavec[u_b], deltastrain);
    EpsT.CopyFrom(deltastrain);
    EpsT.Add(elasto_plastic_integrator.GetState().m_eps_t, 1.);// Adding the last point total strain state
    
    
    // Perform the return mapping algorithm
    elasto_plastic_integrator.ApplyStrainComputeSigma(EpsT, Stress);
    
    TPZTensor<REAL> eps_total_last_state = elasto_plastic_integrator.GetState().m_eps_t;
    TPZTensor<REAL> eps_plastic_last_state = elasto_plastic_integrator.GetState().m_eps_p;
    TPZTensor<REAL> eps_elastic_last_state = eps_total_last_state - eps_plastic_last_state;
    
    // Compute porosity poroelastic correction
    REAL phi_poro = porosity_corrected_2D(eps_elastic_last_state,eps_plastic_last_state,p[0]);
    
    // Computing the density and gravity
    REAL rho_avg = (1.0-phi_poro)*m_rho_s+phi_poro*m_rho_f;
    m_b[0] = rho_avg*m_SimulationData->Gravity()[0];
    m_b[1] = rho_avg*m_SimulationData->Gravity()[1];
    
    // @brief of checking whether the time of diffusion is in the current state or not
    REAL dt = m_SimulationData->dt();
    if (!m_SimulationData->IsCurrentStateQ())
    {
        // Darcy mono-phascis flow
        for (int ip = 0; ip < nphi_p; ip++)
        {
            ef(ip + first_p, 0)        += - weight * (phi_poro/dt) * phip(ip,0);
        }
        return;
    }
    
    TPZFMatrix<REAL> & Sigma_0 = m_SimulationData->PreStress();
    Sigma_0.Zero();
    TPZFNMatrix<9,REAL> delta_S(2,2,0.0);
    
    delta_S(0,0) = ((Stress(_XX_,0))-(Sigma_0(0,0)));
    delta_S(0,1) = ((Stress(_XY_,0))-(Sigma_0(0,1)));
    
    delta_S(1,0) = ((Stress(_XY_,0))-(Sigma_0(1,0)));
    delta_S(1,1) = ((Stress(_YY_,0))-(Sigma_0(1,1)));
       

    for (int iu = 0; iu < nphi_u; iu++) {
        
        // Computing Gradient of the test function for each component
        Grad_vx_i(0,0) = grad_phi_u(0,iu)*axes_u(0,0)+grad_phi_u(1,iu)*axes_u(1,0); // dvx/dx
        Grad_vx_i(1,0) = grad_phi_u(0,iu)*axes_u(0,1)+grad_phi_u(1,iu)*axes_u(1,1); // dvx/dy
        
        Grad_vy_i(0,0) = grad_phi_u(0,iu)*axes_u(0,0)+grad_phi_u(1,iu)*axes_u(1,0); // dvy/dx
        Grad_vy_i(1,0) = grad_phi_u(0,iu)*axes_u(0,1)+grad_phi_u(1,iu)*axes_u(1,1); // dvy/dy
        
        ef(2*iu   + first_u, 0) += weight * ((delta_S(0,0)-m_alpha * p[0])*Grad_vx_i(0,0)+delta_S(0,1)*Grad_vx_i(1,0)-(m_b[0])*phiu(iu,0));
        ef(2*iu+1 + first_u, 0)	+= weight * (delta_S(1,0)*Grad_vy_i(0,0)+(delta_S(1,1)-m_alpha*p[0])*Grad_vy_i(1,0)-(m_b[1])*phiu(iu,0));
        
        
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
    
    // When the residuals expression are solved the memory items are accepted and updated
    if (m_SimulationData->Get_must_accept_solution_Q()) {
//        point_memory.fSigma = Stress;
//        point_memory.fDisplacement = datavec[u_b].sol[0];
//        point_memory.fPlasticState = elasto_plastic_integrator.GetState();
//        point_memory.fPorePressure = p[0];
//        point_memory.fv.resize(m_Dim);
        for (int i = 0;  i < m_Dim; i++)
        {
//            point_memory.fv[i] = - (k/m_eta) * Grad_p(i,0);
        }
        memory_vec[global_point_index] = point_memory;
    }
    
    
     // @brief of checking whether the plasticity is necessary
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
    TPZFNMatrix<36> Dep(6,6,0.0);
    TPZFNMatrix<6>  DeltaStrain(6,1);
    TPZFNMatrix<6>  Stress(6,1);
    
    int ptindex = data.intGlobPtIndex;
    
    
    if (m_UpdateToUseFullDiplacement)
    {
        TPZAdmChunkVector<TMEM> memory_vec = *TPZMatWithMem<TMEM>::fMemory;
        TMEM &point_memory = memory_vec[ptindex];
        point_memory.fPlasticState.m_eps_t.Zero();
        int solsize = data.sol[0].size();
        for(int i=0; i<solsize; i++)
        {
            point_memory.fDisplacement[i] = 0.;
        }
        return;
    }
    
    
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
        // Tenho que fazer esse calculo vetorial tambem????? Acho que nao, o datasol={ux,uy} e data.dsol={{duxdx,duxdy},{duydx,duydy}} aqui
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
    
    
    //#define _XX_ 0
    //#define _YY_ 1
    //#define _XY_ 2
    
    
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << ">>> TPZPMRSCouplPoroPlast<T,TMEM>::ContributePlastic_2D ***";
        sout << "\nIntegration Local Point index = " << data.intGlobPtIndex;
        sout << "\nIntegration Global Point index = " << data.intGlobPtIndex;
        sout << "\ndata.axes = " << data.axes;
        sout << "\nDep " <<endl;
        sout << Dep(_XX_, _XX_) << "\t" << Dep(_XX_, _YY_) << "\t" << Dep(_XX_, _XY_) << "\n";
        sout << Dep(_YY_, _XX_) << "\t" << Dep(_YY_, _YY_) << "\t" << Dep(_YY_, _XY_) << "\n";
        sout << Dep(_XY_, _XX_) << "\t" << Dep(_XY_, _YY_) << "\t" << Dep(_XY_, _XY_) << "\n";
        sout << "\nStress " <<endl;
        sout << Stress(_XX_, 0) << "\t" << Stress(_YY_, 0) << "\t" << Stress(_XY_, 0) << "\n";
        sout << "\nDeltaStrain " <<endl;
        sout << DeltaStrain(_XX_, 0) << "\t" << DeltaStrain(_YY_, 0) << "\t" << DeltaStrain(_XY_, 0) << "\n";
        sout << "data.phi" << data.phi;
        LOGPZ_DEBUG(logger,sout.str().c_str());
    }
#endif
    
    REAL val;
    
    for(int iu = 0; iu < n_phi_u; iu++)
    {
        val = - m_b[0] * phi(iu,0);
        val += Stress(_XX_, 0) * dphiXY(0, iu); //dphixdx
        val += Stress(_XY_, 0) * dphiXY(1, iu); //dphixdy
        ef(2*iu+0 + first_u,0) += weight * val;
        
        val = - m_b[1] * phi(iu,0);
        val += Stress(_XY_, 0) * dphiXY(0, iu); //dphiydx
        val += Stress(_YY_, 0) * dphiXY(1, iu); //dphiydy
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
        
        val  = 2. * Dep(_XX_, _XX_) * Deriv(0, 0); //dvdx*dudx
        val +=      Dep(_XX_, _XY_) * Deriv(0, 1); //dvdx*dudy
        val += 2. * Dep(_XY_, _XX_) * Deriv(1, 0); //dvdy*dudx
        val +=      Dep(_XY_, _XY_) * Deriv(1, 1); //dvdy*dudy
        val *= 0.5;
        ek(2*iu+0 + first_u, 2*ju+0 + first_u) += weight * val;
        
        val  =      Dep(_XX_, _XY_) * Deriv(0, 0); //dvdx*dudx
        val += 2. * Dep(_XX_, _YY_) * Deriv(0, 1); //dvdx*dudy
        val +=      Dep(_XY_, _XY_) * Deriv(1, 0); //dvdy*dudx
        val += 2. * Dep(_XY_, _YY_) * Deriv(1, 1); //dvdy*dudy
        val *= 0.5;
        ek(2*iu+0 + first_u, 2*ju+1 + first_u) += weight * val;
        
        val  = 2. * Dep(_XY_, _XX_) * Deriv(0, 0); //dvdx*dudx
        val +=      Dep(_XY_, _XY_) * Deriv(0, 1); //dvdx*dudy
        val += 2. * Dep(_YY_, _XX_) * Deriv(1, 0); //dvdy*dudx
        val +=      Dep(_YY_, _XY_) * Deriv(1, 1); //dvdy*dudy
        val *= 0.5;
        ek(2*iu+1 + first_u, 2*ju+0 + first_u) += weight * val;
        
        val  =      Dep(_XY_, _XY_) * Deriv(0, 0); //dvdx*dudx
        val += 2. * Dep(_XY_, _YY_) * Deriv(0, 1); //dvdx*dudy
        val +=      Dep(_YY_, _XY_) * Deriv(1, 0); //dvdy*dudx
        val += 2. * Dep(_YY_, _YY_) * Deriv(1, 1); //dvdy*dudy
        val *= 0.5;
        ek(2*iu+1 + first_u, 2*ju+1 + first_u) += weight * val;
            
        }
    }
}


// Contribute Methods being used
template<class T,class TMEM>
void TPZPMRSCouplPoroPlast<T,TMEM>::Contribute_3D(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
//    int u_b = 0;
//    int p_b = 1;
//    
//    // Getting the space functions
//    TPZFMatrix<REAL>     &phiu     =   datavec[u_b].phi;
//    TPZFMatrix<REAL>     &phip     =   datavec[p_b].phi;
//    
//    TPZFMatrix<REAL> &grad_phi_u   =   datavec[u_b].dphix;
//    TPZFMatrix<REAL> &grad_phi_p   =   datavec[p_b].dphix;
//    
//    TPZFNMatrix <9,REAL> &axes_u   =    datavec[u_b].axes;
//    TPZFNMatrix <9,REAL> &axes_p   =    datavec[p_b].axes;
//    
//    // Getting the solutions and derivatives
//    TPZManVector<REAL,3> u  = datavec[u_b].sol[0];
//    TPZManVector<REAL,1> p  = datavec[p_b].sol[0];
//    
//    TPZFNMatrix <9,REAL> du = datavec[u_b].dsol[0];
//    TPZFNMatrix <9,REAL> dp = datavec[p_b].dsol[0];
//    
//    TPZFNMatrix<9,REAL> Grad_p(3,1,0.0),Grad_phi_i(3,1,0.0),Grad_phi_j(3,1,0.0);
//    Grad_p(0,0) = dp(0,0)*axes_p(0,0)+dp(1,0)*axes_p(1,0)+dp(2,0)*axes_p(2,0);
//    Grad_p(1,0) = dp(0,0)*axes_p(0,1)+dp(1,0)*axes_p(1,1)+dp(2,0)*axes_p(2,1);
//    Grad_p(2,0) = dp(0,0)*axes_p(0,2)+dp(1,0)*axes_p(1,2)+dp(2,0)*axes_p(2,2);
//    
//    int nphi_u = phiu.Rows();
//    int nphi_p = phip.Rows();
//    
//    int first_u = 0;
//    int first_p = 3*nphi_u;
//    
//    
//    // Computing Gradient of the Solution
//    TPZFNMatrix<9,REAL> Grad_vx_i(3,1,0.0);
//    TPZFNMatrix<9,REAL> Grad_vy_i(3,1,0.0);
//    TPZFNMatrix<9,REAL> Grad_vz_i(3,1,0.0);
//    
//    TPZFNMatrix<9,REAL> Grad_v(3,3,0.0);
//    TPZFNMatrix<9,REAL> Grad_vx_j(3,1,0.0);
//    TPZFNMatrix<9,REAL> Grad_vy_j(3,1,0.0);
//    TPZFNMatrix<9,REAL> Grad_vz_j(3,1,0.0);
//    
//    REAL dvxdx, dvxdy, dvxdz;
//    REAL dvydx, dvydy, dvydz;
//    REAL dvzdx, dvzdy, dvzdz;
//    
//    REAL duxdx, duxdy, duxdz;
//    REAL duydx, duydy, duydz;
//    REAL duzdx, duzdy, duzdz;
//    
//    
//    // ElastoPlastic Strain Stress Parameters from material memory
//    int global_point_index =  datavec[u_b].intGlobPtIndex;
//    TPZAdmChunkVector<TMEM> & memory_vec = TPZMatWithMem<TMEM>::fMemory;
//    TMEM &point_memory = memory_vec[global_point_index];
//    
//    // Setting for the elastoplastic integrator
//    T elasto_plastic_integrator(this->fPlasticity);
//    elasto_plastic_integrator.SetState(point_memory.fPlasticState);
//    TPZElasticResponse ER;
//    ER.SetUp(m_SimulationData->Get_young(), m_SimulationData->Get_nu());
//    elasto_plastic_integrator.SetElasticResponse(ER);
//    
//    
//    // obtaining the total strain
//    TPZTensor<STATE> Stress = point_memory.fSigma;
//    TPZTensor<REAL> EpsT;
//    TPZFNMatrix<9,STATE> deltastrain(9,1,0.);
//    ComputeDeltaStrainVector(datavec[u_b], deltastrain);
//    EpsT.CopyFrom(deltastrain);
//    EpsT.Add(elasto_plastic_integrator.GetState().fEpsT, 1.);// Adding the last point total strain state
//    
//    // Perform the return mapping algorithm
//    elasto_plastic_integrator.ApplyStrainComputeSigma(EpsT, Stress);
//    
//    TPZTensor<REAL> eps_total_last_state = elasto_plastic_integrator.GetState().fEpsT;
//    TPZTensor<REAL> eps_plastic_last_state = elasto_plastic_integrator.GetState().fEpsP;
//    TPZTensor<REAL> eps_elastic_last_state = eps_total_last_state - eps_plastic_last_state;
//    
//    // Compute porosity poroelastic correction
//    REAL phi_poro = porosity_corrected_3D(eps_elastic_last_state,eps_plastic_last_state,p[0]);
//    
//    // Computing the density and gravity
//    REAL rho_avg = (1.0-phi_poro)*m_rho_s+phi_poro*m_rho_f;
//    m_b[0] = rho_avg*m_SimulationData->Gravity()[0];
//    m_b[1] = rho_avg*m_SimulationData->Gravity()[1];
//    m_b[2] = rho_avg*m_SimulationData->Gravity()[2];
//
//    
//    REAL dt = m_SimulationData->dt();
//    if (!m_SimulationData->IsCurrentStateQ()) {
//        
//        
//        // Darcy mono-phascis flow
//        for (int ip = 0; ip < nphi_p; ip++) {
//            
//            ef(ip + first_p, 0)        += - weight * (phi_poro/dt) * phip(ip,0);
//        }
//        
//        return;
//    }
//    
//    TPZFMatrix<REAL> & Sigma_0 = m_SimulationData->PreStress();
//    Sigma_0.Zero();
//    TPZFNMatrix<9,REAL> S(3,3,0.0);
//    
//    S(0,0) = ((Stress(_XX_,0))-(Sigma_0(0,0)));
//    S(0,1) = ((Stress(_XY_,0))-(Sigma_0(0,1)));
//    S(0,2) = ((Stress(_XZ_,0))-(Sigma_0(0,2)));
//    
//    S(1,0) = ((Stress(_XY_,0))-(Sigma_0(1,0)));
//    S(1,1) = ((Stress(_YY_,0))-(Sigma_0(1,1)));
//    S(1,2) = ((Stress(_YZ_,0))-(Sigma_0(1,2)));
//    
//    S(2,0) = ((Stress(_XZ_,0))-(Sigma_0(2,0)));
//    S(2,1) = ((Stress(_YZ_,0))-(Sigma_0(2,1)));
//    S(2,2) = ((Stress(_ZZ_,0))-(Sigma_0(2,2)));
//    
//    
//    for (int iu = 0; iu < nphi_u; iu++)
//    {
//        // Computing Gradient of the test function for each component
//        for (int d = 0; d < m_Dim; d++)
//        {
//            Grad_vx_i(d,0) = grad_phi_u(d,iu);
//            Grad_vy_i(d,0) = grad_phi_u(d,iu);
//            Grad_vz_i(d,0) = grad_phi_u(d,iu);
//        }
//        
//        ef(3*iu   + first_u, 0)    += weight * ((S(0,0) - m_alpha*p[0]) * Grad_vx_i(0,0) + S(0,1) * Grad_vx_i(1,0) + S(0,2) * Grad_vx_i(2,0) - m_b[0] * phiu(iu,0));
//        ef(3*iu+1 + first_u, 0)    += weight * (S(1,0) * Grad_vy_i(0,0) + (S(1,1) - m_alpha*p[0]) * Grad_vy_i(1,0) + S(1,2) * Grad_vy_i(2,0) - m_b[1] * phiu(iu,0));
//        ef(3*iu+2 + first_u, 0)    += weight * (S(2,0) * Grad_vz_i(0,0) + S(2,1) * Grad_vz_i(1,0) + (S(2,2) - m_alpha*p[0]) * Grad_vz_i(2,0) - m_b[2] * phiu(iu,0));
//        
//        //x
//        dvxdx = Grad_vx_i(0,0);
//        dvxdy = Grad_vx_i(1,0);
//        dvxdz = Grad_vx_i(2,0);
//        
//        //y
//        dvydx = Grad_vy_i(0,0);
//        dvydy = Grad_vy_i(1,0);
//        dvydz = Grad_vy_i(2,0);
//        
//        //z
//        dvzdx = Grad_vz_i(0,0);
//        dvzdy = Grad_vz_i(1,0);
//        dvzdz = Grad_vz_i(2,0);
//        
//        
//        for (int ju = 0; ju < nphi_u; ju++)
//        {
//            
//            // Computing Gradient of the test function for each component
//            for (int d = 0; d < m_Dim; d++)
//            {
//                Grad_vx_j(d,0) = grad_phi_u(d,ju);
//                Grad_vy_j(d,0) = grad_phi_u(d,ju);
//                Grad_vz_j(d,0) = grad_phi_u(d,ju);
//            }
//            
//            //x
//            duxdx = Grad_vx_j(0,0);
//            duxdy = Grad_vx_j(1,0);
//            duxdz = Grad_vx_j(2,0);
//            
//            //y
//            duydx = Grad_vy_j(0,0);
//            duydy = Grad_vy_j(1,0);
//            duydz = Grad_vy_j(2,0);
//            
//            //z
//            duzdx = Grad_vz_j(0,0);
//            duzdy = Grad_vz_j(1,0);
//            duzdz = Grad_vz_j(2,0);
//            
//            // Gradient 1
//            ek(3*iu   + first_u, 3*ju    + first_u) += weight * ((m_lambda + 2.*m_mu)*duxdx*dvxdx + m_mu*duxdy*dvxdy + m_mu*duxdz*dvxdz);
//            ek(3*iu   + first_u, 3*ju+1  + first_u) += weight * (m_lambda*duydy*dvxdx + m_mu*duydx*dvxdy);
//            ek(3*iu   + first_u, 3*ju+2  + first_u) += weight * (m_lambda*duzdz*dvxdx + m_mu*duzdx*dvxdz);
//            
//            // Gradient 2
//            ek(3*iu+1 + first_u, 3*ju    + first_u) += weight * (m_lambda*duxdx*dvydy + m_mu*duxdy*dvydx);
//            ek(3*iu+1 + first_u, 3*ju+1  + first_u) += weight * ((m_lambda + 2.*m_mu)*duydy*dvydy + m_mu*duydx*dvydx + m_mu*duydz*dvydz);
//            ek(3*iu+1 + first_u, 3*ju+2  + first_u) += weight * (m_lambda*duzdz*dvydy + m_mu*duzdy*dvydz);
//            
//            // Gradient 3
//            ek(3*iu+2 + first_u, 3*ju    + first_u) += weight * (m_lambda*duxdx*dvzdz + m_mu*duxdz*dvzdx);
//            ek(3*iu+2 + first_u, 3*ju+1  + first_u) += weight * (m_lambda*duydy*dvzdz + m_mu*duydz*dvzdy);
//            ek(3*iu+2 + first_u, 3*ju+2  + first_u) += weight * ((m_lambda + 2.*m_mu)*duzdz*dvzdz + m_mu*duzdx*dvzdx + m_mu*duzdy*dvzdy);
//        }
//    }
//    
//    TPZFNMatrix<9,REAL> dv(3,1,0.0);
//    
//    //    Matrix -Qc
//    //    Coupling matrix
//    for(int iu = 0; iu < nphi_u; iu++ )
//    {
//        
//        // Computing Gradient of the test function for each component
//        Grad_vx_i(0,0) = grad_phi_u(0,iu)*axes_u(0,0)+grad_phi_u(1,iu)*axes_u(1,0)+grad_phi_u(2,iu)*axes_u(2,0); // dvx/dx
//        Grad_vx_i(1,0) = grad_phi_u(0,iu)*axes_u(0,1)+grad_phi_u(1,iu)*axes_u(1,1)+grad_phi_u(2,iu)*axes_u(2,1); // dvx/dy
//        Grad_vx_i(2,0) = grad_phi_u(0,iu)*axes_u(0,2)+grad_phi_u(1,iu)*axes_u(1,2)+grad_phi_u(2,iu)*axes_u(2,2); // dvx/dz
//        
//        Grad_vy_i(0,0) = grad_phi_u(0,iu)*axes_u(0,0)+grad_phi_u(1,iu)*axes_u(1,0)+grad_phi_u(2,iu)*axes_u(2,0); // dvy/dx
//        Grad_vy_i(1,0) = grad_phi_u(0,iu)*axes_u(0,1)+grad_phi_u(1,iu)*axes_u(1,1)+grad_phi_u(2,iu)*axes_u(2,1); // dvy/dy
//        Grad_vy_i(2,0) = grad_phi_u(0,iu)*axes_u(0,2)+grad_phi_u(1,iu)*axes_u(1,2)+grad_phi_u(2,iu)*axes_u(2,2); // dvy/dz
//        
//        Grad_vz_i(0,0) = grad_phi_u(0,iu)*axes_u(0,0)+grad_phi_u(1,iu)*axes_u(1,0)+grad_phi_u(2,iu)*axes_u(2,0); // dvz/dx
//        Grad_vz_i(1,0) = grad_phi_u(0,iu)*axes_u(0,1)+grad_phi_u(1,iu)*axes_u(1,1)+grad_phi_u(2,iu)*axes_u(2,1); // dvz/dy
//        Grad_vz_i(2,0) = grad_phi_u(0,iu)*axes_u(0,2)+grad_phi_u(1,iu)*axes_u(1,2)+grad_phi_u(2,iu)*axes_u(2,2); // dvz/dz
//        
//        for(int jp = 0; jp < nphi_p; jp++)
//        {
//            
//            ek(3*iu+0,first_p+jp) += (-1.)* weight * m_alpha * phip(jp,0) * Grad_vx_i(0,0);
//            ek(3*iu+1,first_p+jp) += (-1.)* weight * m_alpha * phip(jp,0) * Grad_vy_i(1,0);
//            ek(3*iu+2,first_p+jp) += (-1.)* weight * m_alpha * phip(jp,0) * Grad_vz_i(2,0);
//        }
//    }
//    
//    //    Matrix QcˆT
//    //    Coupling matrix transpose
//    for(int ip = 0; ip < nphi_p; ip++ )
//    {
//        
//        for(int ju = 0; ju < nphi_u; ju++)
//        {
//            dv(0,0) = grad_phi_u(0,ju)*axes_u(0,0)+grad_phi_u(1,ju)*axes_u(1,0)+grad_phi_u(2,ju)*axes_u(2,0);
//            dv(1,0) = grad_phi_u(0,ju)*axes_u(0,1)+grad_phi_u(1,ju)*axes_u(1,1)+grad_phi_u(2,ju)*axes_u(2,1);
//            dv(2,0) = grad_phi_u(0,ju)*axes_u(0,2)+grad_phi_u(1,ju)*axes_u(1,2)+grad_phi_u(2,ju)*axes_u(2,2);
//            
//            ek(first_p+ip,3*ju+0) += (1./dt) * weight * m_alpha * dv(0,0) * phip(ip,0);
//            ek(first_p+ip,3*ju+1) += (1./dt) * weight * m_alpha * dv(1,0) * phip(ip,0);
//            ek(first_p+ip,3*ju+2) += (1./dt) * weight * m_alpha * dv(2,0) * phip(ip,0);
//        }
//    }
//    
//    
//    /** @brief Rudnicki diffusion coefficient */
//    /** J. W. Rudnicki. Fluid mass sources and point forces in linear elastic diffusive solids. Journal of Mechanics of Materials, 5:383–393, 1986. */
//    REAL k = 0.0;
//    m_k_model = 1;
//    k_permeability(phi_poro,k);
//    m_lambdau = 1.1 * m_lambda;
//    REAL c = (k/m_eta)*(m_lambdau-m_lambda)*(m_lambda + 2.0*m_mu)/(m_alpha*m_alpha*(m_lambdau + 2.0*m_mu));
//    
//    // Darcy mono-phascis flow
//    for (int ip = 0; ip < nphi_p; ip++)
//    {
//        
//        Grad_phi_i(0,0) = grad_phi_p(0,ip)*axes_p(0,0)+grad_phi_p(1,ip)*axes_p(1,0)+grad_phi_p(2,ip)*axes_p(2,0);
//        Grad_phi_i(1,0) = grad_phi_p(0,ip)*axes_p(0,1)+grad_phi_p(1,ip)*axes_p(1,1)+grad_phi_p(2,ip)*axes_p(2,1);
//        Grad_phi_i(2,0) = grad_phi_p(0,ip)*axes_p(0,2)+grad_phi_p(1,ip)*axes_p(1,2)+grad_phi_p(2,ip)*axes_p(2,2);
//        
//        REAL dot = 0.0;
//        for (int i = 0;  i < m_Dim; i++)
//        {
//            dot += Grad_p(i,0) * Grad_phi_i(i,0);
//        }
//        
//        ef(ip + first_p, 0)        +=  weight * ( c * dot + (phi_poro/dt) * phip(ip,0));
//        
//        for (int jp = 0; jp < nphi_p; jp++)
//        {
//            
//            Grad_phi_j(0,0) = grad_phi_p(0,jp)*axes_p(0,0)+grad_phi_p(1,jp)*axes_p(1,0)+grad_phi_p(2,jp)*axes_p(2,0);
//            Grad_phi_j(1,0) = grad_phi_p(0,jp)*axes_p(0,1)+grad_phi_p(1,jp)*axes_p(1,1)+grad_phi_p(2,jp)*axes_p(2,1);
//            Grad_phi_j(2,0) = grad_phi_p(0,jp)*axes_p(0,2)+grad_phi_p(1,jp)*axes_p(1,2)+grad_phi_p(2,jp)*axes_p(2,2);
//            
//            REAL dot = 0.0;
//            for (int i = 0;  i < m_Dim; i++)
//            {
//                dot += Grad_phi_j(i,0) * Grad_phi_i(i,0);
//            }
//            
//            ek(ip + first_p, jp + first_p)        += weight * (c * dot + (m_Se/dt) * phip(jp,0) * phip(ip,0) );
//        }
//    }
//    
//#ifdef LOG4CXX
//    if(logger->isDebugEnabled())
//    {
//        std::stringstream sout;
//        sout << "<<< TPZPMRSCouplPoroPlast<T,TMEM>::Contribute_2D ***";
//        sout << " Resultant rhs vector:\n" << ef;
//        sout << " Resultant stiff vector:\n" << ek;
//        LOGPZ_DEBUG(logger,sout.str().c_str());
//    }
//#endif
//    
//    // When the residuals expression are solved the memory items are accepted and updated
//    if (m_SimulationData->Get_must_accept_solution_Q()) {
//        point_memory.fSigma = Stress;
//        point_memory.fDisplacement = datavec[u_b].sol[0];
//        point_memory.fPlasticState = elasto_plastic_integrator.GetState();
//        point_memory.fPorePressure = p[0];
//        point_memory.fv.resize(m_Dim);
//        for (int i = 0;  i < m_Dim; i++)
//        {
//            point_memory.fv[i] = - (k/m_eta) * Grad_p(i,0);
//        }
//        memory_vec[global_point_index] = point_memory;
//    }
//    
//    
//    // @brief of checking whether the plasticity is necessary
//    if (m_SetRunPlasticity)
//    {
//        ContributePlastic_3D(datavec[0],weight,ek,ef);
//        
//#ifdef LOG4CXX
//        if(logger->isDebugEnabled())
//        {
//            std::stringstream sout;
//            sout << "<<< TPZPMRSCouplPoroPlast<T,TMEM>::ContributePlastic_2D ***";
//            sout << " Resultant rhs vector:\n" << ef;
//            sout << " Resultant stiff vector:\n" << ek;
//            LOGPZ_DEBUG(logger,sout.str().c_str());
//        }
//#endif
//        return;
//    }
}


/** @brief of contribute of plasticity in 2 dimensional */
template<class T,class TMEM>
void TPZPMRSCouplPoroPlast<T,TMEM>::ContributePlastic_3D(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef)
{
    TPZFMatrix<REAL> &dphi = data.dphix, dphiXY;
    TPZFMatrix<REAL> &phi  = data.phi;
    dphiXY = dphi;
    
    int first_u = 0;
    const int n_phi_u = phi.Rows();
    
    TPZFNMatrix<9>  Deriv(3,3);
    TPZFNMatrix<36> Dep(6,6);
    TPZFNMatrix<6>  DeltaStrain(6,1);
    TPZFNMatrix<6>  Stress(6,1);
    
    int ptindex = data.intGlobPtIndex;
    
    
    if (m_UpdateToUseFullDiplacement)
    {
        TPZAdmChunkVector<TMEM> memory_vec = *TPZMatWithMem<TMEM>::fMemory;
        TMEM &point_memory = memory_vec[ptindex];
//        point_memory.fPlasticState.fEpsT.Zero();
        int solsize = data.sol[0].size();
        for(int i=0; i<solsize; i++)
        {
//           point_memory.fDisplacement[i] = 0.;
        }
        return;
    }
    
    
    
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
        // Tenho que fazer esse calculo vetorial tambem????? Acho que nao, o datasol={ux,uy} e data.dsol={{duxdx,duxdy},{duydx,duydy}} aqui
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
    
    //#define _XX_ 0
    //#define _XY_ 1
    //#define _XZ_ 2
    //#define _YY_ 3
    //#define _YZ_ 4
    //#define _ZZ_ 5
    
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << ">>> TPZPMRSCouplPoroPlast<T,TMEM>::ContributePlastic_3D ***";
        sout << "\nIntegration Local Point index = " << data.intGlobPtIndex;
        sout << "\nIntegration Global Point index = " << data.intGlobPtIndex;
        sout << "\ndata.axes = " << data.axes;
        sout << "\nDep " <<endl;
        sout << Dep(_XX_,_XX_) << "\t" << Dep(_XX_,_XY_) << "\t" << Dep(_XX_,_XZ_)  << "\t" << Dep(_XX_,_YY_) << "\t" << Dep(_XX_,_YZ_) << "\t" << Dep(_XX_,_ZZ_) << "\n";
        sout << Dep(_XY_,_XX_) << "\t" << Dep(_XY_,_XY_) << "\t" << Dep(_XY_,_XZ_)  << "\t" << Dep(_XY_,_YY_) << "\t" << Dep(_XY_,_YZ_) << "\t" << Dep(_XY_,_ZZ_) << "\n";
        sout << Dep(_XZ_,_XX_) << "\t" << Dep(_XZ_,_XY_) << "\t" << Dep(_XZ_,_XZ_)  << "\t" << Dep(_XZ_,_YY_) << "\t" << Dep(_XZ_,_YZ_) << "\t" << Dep(_XZ_,_ZZ_) << "\n";
        sout << Dep(_YY_,_XX_) << "\t" << Dep(_YY_,_XY_) << "\t" << Dep(_YY_,_XZ_)  << "\t" << Dep(_YY_,_YY_) << "\t" << Dep(_YY_,_YZ_) << "\t" << Dep(_YY_,_ZZ_) << "\n";
        sout << Dep(_YZ_,_XX_) << "\t" << Dep(_YZ_,_XY_) << "\t" << Dep(_YZ_,_XZ_)  << "\t" << Dep(_YZ_,_YY_) << "\t" << Dep(_YZ_,_YZ_) << "\t" << Dep(_YZ_,_ZZ_) << "\n";
        sout << Dep(_ZZ_,_XX_) << "\t" << Dep(_ZZ_,_XY_) << "\t" << Dep(_ZZ_,_XZ_)  << "\t" << Dep(_ZZ_,_YY_) << "\t" << Dep(_ZZ_,_YZ_) << "\t" << Dep(_ZZ_,_ZZ_) << "\n";
        sout << "\nStress " <<endl;
        sout << Stress(_XX_,0) << "\t" << Stress(_XY_,0) << "\t" << Stress(_XZ_,0) << "\t" << Stress(_YY_,0) << "\t" << Stress(_YZ_,0) << "\t" << Stress(_ZZ_,0) <<"\n";
        sout << "\nDeltaStrain " <<endl;
        sout << DeltaStrain(_XX_,0)<<"\t"<<DeltaStrain(_XY_,0) <<"\t"<< DeltaStrain(_XZ_,0)<<"\t"<<DeltaStrain(_YY_,0)<<"\t"<<DeltaStrain(_YZ_,0)<<"\t"<<DeltaStrain(_ZZ_,0)<<"\n";
        sout << "data.phi" << data.phi;
        LOGPZ_DEBUG(logger,sout.str().c_str());
    }
#endif
    
    
    REAL val,val2,val3,val4,val5,val6,val7,val8,val9,val10;
    
    for(int iu = 0; iu < n_phi_u; iu++)
    {
        val = - m_b[0] * phi(iu,0);
        val += Stress(_XX_,0) * dphiXY(0,iu); //dphixdx
        val += Stress(_XY_,0) * dphiXY(1,iu); //dphixdy
        val += Stress(_XZ_,0) * dphiXY(2,iu); //dphixdz
        ef(2*iu+0 + first_u,0) += weight * val;
        
        val = - m_b[1] * phi(iu,0);
        val += Stress(_XY_,0) * dphiXY(0,iu); //dphiydx
        val += Stress(_YY_,0) * dphiXY(1,iu); //dphiydy
        val += Stress(_YZ_,0) * dphiXY(2,iu); //dphiydz
        ef(2*iu+1 + first_u,0) += weight * val;
        
        val = - m_b[2] * phi(iu,0);
        val += Stress(_XZ_,0) * dphiXY(0,iu); //dphizdx
        val += Stress(_YZ_,0) * dphiXY(1,iu); //dphizdy
        val += Stress(_ZZ_,0) * dphiXY(2,iu); //dphizdz
        ef(2*iu+2 + first_u,0) += weight * val;
        
        
        for (int ju = 0; ju < n_phi_u; ju++)
        {
            for (int ud = 0; ud < 3; ud++)
            {
                for (int vd = 0; vd < 3; vd++)
                {
                    Deriv(vd, ud) = dphiXY(vd, iu) * dphiXY(ud, ju);
                }//ud
            }//vd
            
            
            // The matrix is:
            //  {{dvdx*dudx, dvdx*dudy, dvdx*dudz},
            //  {dvdy*dudx , dvdy*dudy, dvdy*dudz},
            //  {dvdz*dudx , dvdz*dudy, dvdz*dudz}}
            //
            
            //First equation Dot[Sigma1, gradV1]
            val2  = 2. * Dep(_XX_,_XX_) * Deriv(0,0);//dvdx*dudx
            val2 +=      Dep(_XX_,_XY_) * Deriv(0,1);//dvdx*dudy
            val2 +=	     Dep(_XX_,_XZ_) * Deriv(0,2);//dvdx*dudz
            val2 += 2. * Dep(_XY_,_XX_) * Deriv(1,0);//dvdy*dudx
            val2 +=      Dep(_XY_,_XY_) * Deriv(1,1);//dvdy*dudy
            val2 +=      Dep(_XY_,_XZ_) * Deriv(1,2);//dvdy*dudz
            val2 += 2. * Dep(_XZ_,_XX_) * Deriv(2,0);//dvdz*dudx
            val2 +=      Dep(_XZ_,_XY_) * Deriv(2,1);//dvdz*dudy
            val2 +=      Dep(_XZ_,_XZ_) * Deriv(2,2);//dvdz*dudz
            val2 *= 0.5;
            ek(2*iu+0 + first_u, 2*ju+0 + first_u) += weight * val2;
            
            val3  =      Dep(_XX_,_XY_) * Deriv(0,0);
            val3 += 2. * Dep(_XX_,_YY_) * Deriv(0,1);
            val3 +=      Dep(_XX_,_YZ_) * Deriv(0,2);
            val3 +=      Dep(_XY_,_XY_) * Deriv(1,0);
            val3 += 2. * Dep(_XY_,_YY_) * Deriv(1,1);
            val3 +=      Dep(_XY_,_YZ_) * Deriv(1,2);
            val3 +=      Dep(_XZ_,_XY_) * Deriv(2,0);
            val3 += 2. * Dep(_XZ_,_YY_) * Deriv(2,1);
            val3 +=      Dep(_XZ_,_YZ_) * Deriv(2,2);
            val3 *= 0.5;
            ek(2*iu+0 + first_u, 2*ju+1 + first_u) += weight * val3;
            
            val4  =      Dep(_XX_,_XZ_) * Deriv(0,0);
            val4 +=      Dep(_XX_,_YZ_) * Deriv(0,1);
            val4 += 2. * Dep(_XX_,_ZZ_) * Deriv(0,2);//
            val4 +=      Dep(_XY_,_XZ_) * Deriv(1,0);
            val4 +=      Dep(_XY_,_YZ_) * Deriv(1,1);
            val4 += 2. * Dep(_XY_,_ZZ_) * Deriv(1,2);//
            val4 +=      Dep(_XZ_,_XZ_) * Deriv(2,0);
            val4 +=      Dep(_XZ_,_YZ_) * Deriv(2,1);
            val4 += 2. * Dep(_XZ_,_ZZ_) * Deriv(2,2);
            val4 *= 0.5;
            ek(2*iu+0 + first_u, 2*ju+2 + first_u) += weight * val4;
            
            //Second equation Dot[Sigma2, gradV2]
            val5  = 2. * Dep(_XY_,_XX_) * Deriv(0,0);
            val5 +=      Dep(_XY_,_XY_) * Deriv(0,1);
            val5 +=      Dep(_XY_,_XZ_) * Deriv(0,2);
            val5 += 2. * Dep(_YY_,_XX_) * Deriv(1,0);
            val5 +=      Dep(_YY_,_XY_) * Deriv(1,1);
            val5 +=      Dep(_YY_,_XZ_) * Deriv(1,2);
            val5 += 2. * Dep(_YZ_,_XX_) * Deriv(2,0);
            val5 +=      Dep(_YZ_,_XY_) * Deriv(2,1);
            val5 +=      Dep(_YZ_,_XZ_) * Deriv(2,2);
            val5 *= 0.5;
            ek(2*iu+1 + first_u, 2*ju+0 + first_u) += weight * val5;
            
            val6  =      Dep(_XY_,_XY_) * Deriv(0,0);
            val6 += 2. * Dep(_XY_,_YY_) * Deriv(0,1);
            val6 +=      Dep(_XY_,_YZ_) * Deriv(0,2);
            val6 +=      Dep(_YY_,_XY_) * Deriv(1,0);
            val6 += 2. * Dep(_YY_,_YY_) * Deriv(1,1);
            val6 +=      Dep(_YY_,_YZ_) * Deriv(1,2);
            val6 +=      Dep(_YZ_,_XY_) * Deriv(2,0);
            val6 += 2. * Dep(_YZ_,_YY_) * Deriv(2,1);
            val6 +=      Dep(_YZ_,_YZ_) * Deriv(2,2);
            val6 *= 0.5;
            ek(2*iu+1 + first_u, 2*ju+1 + first_u) += weight * val6;
            
            val7  =      Dep(_XY_,_XZ_) * Deriv(0,0);
            val7 +=      Dep(_XY_,_YZ_) * Deriv(0,1);
            val7 += 2. * Dep(_XY_,_ZZ_) * Deriv(0,2);//
            val7 +=      Dep(_YY_,_XZ_) * Deriv(1,0);
            val7 +=      Dep(_YY_,_YZ_) * Deriv(1,1);
            val7 += 2. * Dep(_YY_,_ZZ_) * Deriv(1,2);//
            val7 +=      Dep(_YZ_,_XZ_) * Deriv(2,0);
            val7 +=      Dep(_YZ_,_YZ_) * Deriv(2,1);
            val7 += 2. * Dep(_YZ_,_ZZ_) * Deriv(2,2);
            val7 *= 0.5;
            ek(2*iu+1 + first_u, 2*ju+2 + first_u) += weight * val7;
            
            
            //Third equation Dot[Sigma3, gradV3]
            val8  = 2. * Dep(_XZ_,_XX_) * Deriv(0,0);
            val8 +=      Dep(_XZ_,_XY_) * Deriv(0,1);
            val8 +=      Dep(_XZ_,_XZ_) * Deriv(0,2);
            val8 += 2. * Dep(_YZ_,_XX_) * Deriv(1,0);
            val8 +=      Dep(_YZ_,_XY_) * Deriv(1,1);
            val8 +=      Dep(_YZ_,_XZ_) * Deriv(1,2);
            val8 += 2. * Dep(_ZZ_,_XX_) * Deriv(2,0);//
            val8 +=      Dep(_ZZ_,_XY_) * Deriv(2,1);//
            val8 +=      Dep(_ZZ_,_XZ_) * Deriv(2,2);
            val8 *= 0.5;
            ek(2*iu+2 + first_u, 2*ju+0 + first_u) += weight * val8;
            
            val9  =      Dep(_XZ_,_XY_) * Deriv(0,0);
            val9 += 2. * Dep(_XZ_,_YY_) * Deriv(0,1);
            val9 +=      Dep(_XZ_,_YZ_) * Deriv(0,2);
            val9 +=      Dep(_YZ_,_XY_) * Deriv(1,0);
            val9 += 2. * Dep(_YZ_,_YY_) * Deriv(1,1);
            val9 +=      Dep(_YZ_,_YZ_) * Deriv(1,2);
            val9 +=      Dep(_ZZ_,_XY_) * Deriv(2,0);//
            val9 += 2. * Dep(_ZZ_,_YY_) * Deriv(2,1);//
            val9 +=      Dep(_ZZ_,_YZ_) * Deriv(2,2);
            val9 *= 0.5;
            ek(2*iu+2 + first_u, 2*ju+1 + first_u) += weight * val9;
            
            val10  =      Dep(_XZ_,_XZ_) * Deriv(0,0);
            val10 +=      Dep(_XZ_,_YZ_) * Deriv(0,1);
            val10 += 2. * Dep(_XZ_,_ZZ_) * Deriv(0,2);
            val10 +=      Dep(_YZ_,_XZ_) * Deriv(1,0);
            val10 +=      Dep(_YZ_,_YZ_) * Deriv(1,1);
            val10 += 2. * Dep(_YZ_,_ZZ_) * Deriv(1,2);
            val10 +=      Dep(_ZZ_,_XZ_) * Deriv(2,0);
            val10 +=      Dep(_ZZ_,_YZ_) * Deriv(2,1);
            val10 += 2. * Dep(_ZZ_,_ZZ_) * Deriv(2,2);//
            val10 *= 0.5;
            ek(2*iu+2 + first_u, 2*ju+2 + first_u) += weight * val10;
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
    
    // *********************** Boundaries *******************************

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
                ef(2*in+0,0)      += BIGNUMBER*(u[0] - v[0])*phiu(in,0)*weight;    // X displacement Value
                ef(2*in+1,0)      += BIGNUMBER*(u[1] - v[1])*phiu(in,0)*weight;    // y displacement Value

                for (jn = 0 ; jn < phru; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(2*in+0,2*jn+0)    += BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;    // X displacement
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
                ef(2*in+0,0)    += -1.0 * weight * v[0] * phiu(in,0);        //    Tnx
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
                ef(2*in+0,0)    += -1.0 * weight * tn * n[0] * phiu(in,0);        //    Tnx
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
                ef(2*in+0,0)        += BIGNUMBER*(u[0] - v[0])*phiu(in,0)*weight;    // X displacement Value
                ef(2*in+1,0)      += BIGNUMBER*(u[1] - v[1])*phiu(in,0)*weight;    // y displacement Value

                for (jn = 0 ; jn < phru; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(2*in+0,2*jn+0)        += BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;    // X displacement
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
                ef(2*in+0,0)     += -1.0 * weight * v[0] * phiu(in,0);        //    Tnx
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
                ef(2*in+0,0)      += -1.0 * weight * tn * n[0] * phiu(in,0);        //    Tnx
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
                ef(2*in+0,0)      += BIGNUMBER*(u[0] - v[0])*phiu(in,0)*weight;    // X displacement Value
                ef(2*in+1,0)      += BIGNUMBER*(u[1] - v[1])*phiu(in,0)*weight;    // y displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(2*in+0,2*jn+0)    += BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;    // X displacement
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
                ef(3*in+0,0)	+= BIGNUMBER*(u[0] - v[0])*phiu(in,0)*weight;	// X displacement Value
                ef(3*in+1,0)	+= BIGNUMBER*(u[1] - v[1])*phiu(in,0)*weight;	// y displacement Value
                ef(3*in+2,0)	+= BIGNUMBER*(u[2] - v[2])*phiu(in,0)*weight;	// Z displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in+0,3*jn+0)	+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
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
            
        // ************************* Boundaries *******************************

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
                ef(3*in+0,0)	+= BIGNUMBER*(u[0] - v[0])*phiu(in,0)*weight;	// X displacement Value
                ef(3*in+1,0)	+= BIGNUMBER*(u[1] - v[1])*phiu(in,0)*weight;	// y displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in+0,3*jn+0)	+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
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
                ef(3*in+0,0)	+= BIGNUMBER*(u[0] - v[0])*phiu(in,0)*weight;	// X displacement Value
                ef(3*in+2,0)	+= BIGNUMBER*(u[2] - v[1])*phiu(in,0)*weight;	// Z displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in+0,3*jn+0)	+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
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
                ef(3*in+0,0)	+= -1.0 * weight * v[0] * phiu(in,0);		//	Tnx
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
                ef(3*in+0,0)	+= -1.0 * weight * tn * n[0] * phiu(in,0);		//	Tnx
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
                ef(3*in+0,0)	+= BIGNUMBER*(u[0] - v[0])*phiu(in,0)*weight;	// X displacement Value
                ef(3*in+1,0)	+= BIGNUMBER*(u[1] - v[1])*phiu(in,0)*weight;	// Y displacement Value
                ef(3*in+2,0)	+= BIGNUMBER*(u[2] - v[2])*phiu(in,0)*weight;	// Z displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in+0,3*jn+0)	+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
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
                ef(3*in+0,0)	+= BIGNUMBER*(u[0] - v[0])*phiu(in,0)*weight;	// X displacement Value
                ef(3*in+1,0)	+= BIGNUMBER*(u[1] - v[1])*phiu(in,0)*weight;	// Y displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in+0,3*jn+0)	+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
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
                ef(3*in+0,0)	+= BIGNUMBER*(u[0] - v[0])*phiu(in,0)*weight;	// X displacement Value
                ef(3*in+2,0)	+= BIGNUMBER*(u[2] - v[1])*phiu(in,0)*weight;	// Z displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in+0,3*jn+0)	+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
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
          
            
        case 18 : //Du_time_Dp
        {
            REAL v[4];
            v[0] = bc.Val2()(0,0);    //    Ux displacement
            v[1] = bc.Val2()(1,0);    //    Uy displacement
            v[2] = bc.Val2()(2,0);	  //	Uz displacement
            v[3] = bc.Val2()(3,0);    //    Pressure
            //    Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //    Contribution for load Vector
                ef(2*in+0,0)      += BIGNUMBER*(u[0] - v[0])*phiu(in,0)*weight;    // X displacement Value
                ef(2*in+1,0)      += BIGNUMBER*(u[1] - v[1])*phiu(in,0)*weight;    // Y displacement Value
                ef(2*in+2,0)      += BIGNUMBER*(u[2] - v[2])*phiu(in,0)*weight;    // Z displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(2*in+0,2*jn+0)    += BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;    // X displacement
                    ek(2*in+1,2*jn+1)    += BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;    // Y displacement
                    ek(2*in+2,2*jn+2)    += BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;    // Z displacement
                }
            }
            
            //    Diffusion Equation
            REAL p_s = p[0];
            REAL d_p = (p_s-v[3]);
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
            
        case 19 : // Ntn_time_Dp
        {
            //    Neumann condition for each state variable
            REAL v[4];
            v[0] = bc.Val2()(0,0);    //    Tnx
            v[1] = bc.Val2()(1,0);    //    Tny
            v[2] = bc.Val2()(2,0);	  //	Tnz
            v[3] = bc.Val2()(3,0);    //    Pressure
            
            //    Elasticity Equation
            for(in = 0 ; in <phru; in++)
            {
                //    Normal Tension Components on neumman boundary
                ef(2*in+0,0)    += -1.0 * weight * v[0] * phiu(in,0);        //    Tnx
                ef(2*in+1,0)    += -1.0 * weight * v[1] * phiu(in,0);        //    Tny
                ef(2*in+2,0)    += -1.0 * weight * v[2] * phiu(in,0);        //    Tnz
            }
            
            //    Diffusion Equation
            REAL p_s = p[0];
            REAL d_p = (p_s-v[3]);
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
    out << "Properties for this material: \n";
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
    //**** The VariableIndex must be selected more than 17 because of the previous VariableIndex associated with TPZMatElastoPlastic ****
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
    
    return TPZMatWithMem<TMEM>::VariableIndex(name);
}

template<class T,class TMEM>
int TPZPMRSCouplPoroPlast<T,TMEM>::NSolutionVariables(int var)
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
    
    return TPZMatWithMem<TMEM>::NSolutionVariables(var);
}


//	Calculate Secondary variables based on ux, uy, Pore pressure and their derivatives
template<class T,class TMEM>
void TPZPMRSCouplPoroPlast<T,TMEM>::Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout)
{
    int global_point = data.intGlobPtIndex;
    TPZAdmChunkVector<TMEM> memory_vec = *TPZMatWithMem<TMEM>::fMemory;
    TMEM &Memory = memory_vec[global_point];
    
    T plasticloc(this->m_plasticity_model);
    plasticloc.SetState(Memory.m_elastoplastic_state);
    
    TPZTensor<STATE> Sigma ;//= Memory.fSigma;
    
    STATE normdsol = Norm(data.dsol[0]);
    
    if (normdsol != 0.)
    {
        TPZTensor<REAL> EpsT;
        TPZFNMatrix<6,STATE> deltastrain(6,1,0.);
        ComputeDeltaStrainVector(data, deltastrain);
        
        EpsT.CopyFrom(deltastrain);
        EpsT.Add(plasticloc.GetState().m_eps_t, 1.);
        plasticloc.ApplyStrainComputeSigma(EpsT, Sigma);
    }
    
    TPZPlasticState<STATE> PState = plasticloc.GetState();
    TPZTensor<REAL> totalStrain = PState.m_eps_t;
    TPZTensor<REAL> plasticStrain = PState.m_eps_p;
    TPZTensor<REAL> elasticStrain = totalStrain - plasticStrain;
    
    // The values of displacement and total stress
    TPZManVector<REAL,3>  u ;//= Memory.fDisplacement;
    TPZTensor<REAL> totalStress = Sigma;
    
    // The values of pore pressure and darcy velocity
    STATE p ;//= Memory.fPorePressure;
    TPZManVector<STATE,3> v ;//= Memory.fv;


    // ************************************** The value of parameters ************************
    
    // ************************	Total Strain Variables ************************
    //	Total Volumetric Strain
    if(var == 18)
    {
        Solout[0] = totalStrain.XX() + totalStrain.YY() + totalStrain.ZZ();
        return;
    }
    
    //	Total Strain in XX Direction
    if(var == 19)
    {
        Solout[0] = totalStrain.XX();
        return;
    }
    
    //	Total Strain in YY Direction
    if(var == 20)
    {
        Solout[0] = totalStrain.YY();
        return;
    }
    
    //	Total Strain in ZZ Direction
    if(var == 21)
    {
        Solout[0] = totalStrain.ZZ();
        return;
    }
    
    //	Total Strain in XY Direction
    if(var == 22)
    {
        Solout[0] = totalStrain.XY();
        return;
    }
    
    //	Total Strain in XZ Direction
    if(var == 23)
    {
        Solout[0] = totalStrain.XZ();
        return;
    }
    
    //	Total Strain in YZ Direction
    if(var == 24)
    {
        Solout[0] = totalStrain.YZ();
        return;
    }
    
    // ************************	Elastic Strain Variables ************************
    //	Elastic Volumetric Strain
    if(var == 25)
    {
        Solout[0] = elasticStrain.XX() + elasticStrain.YY() + elasticStrain.ZZ();;
        return;
    }
    
    //	Elastic Strain in XX Direction
    if(var == 26)
    {
        Solout[0] = elasticStrain.XX();
        return;
    }
    
    //	Elastic Strain in YY Direction
    if(var == 27)
    {
        Solout[0] = elasticStrain.YY();
        return;
    }
    
    //	Elastic Strain in ZZ Direction
    if(var == 28)
    {
        Solout[0] = elasticStrain.ZZ();
        return;
    }
    
    //	Elastic Strain in XY Direction
    if(var == 29)
    {
        Solout[0] = elasticStrain.XY();
        return;
    }
    
    //	Elastic Strain in XZ Direction
    if(var == 30)
    {
        Solout[0] = elasticStrain.XZ();
        return;
    }
    
    //	Elastic Strain in YZ Direction
    if(var == 31)
    {
        Solout[0] = elasticStrain.YZ();
        return;
    }
    
    // ************************	Plastic Strain Variables ************************
    //	Plastic Volumetric Strain
    if(var == 32)
    {
        Solout[0] = plasticStrain.XX() + plasticStrain.YY() + plasticStrain.ZZ();;
        return;
    }
    
    //	Plastic Strain in XX Direction
    if(var == 33)
    {
        Solout[0] = plasticStrain.XX();
        return;
    }
    
    //	Plastic Strain in YY Direction
    if(var == 34)
    {
        Solout[0] = plasticStrain.YY();
        return;
    }
    
    //	Plastic Strain in ZZ Direction
    if(var == 35)
    {
        Solout[0] = plasticStrain.ZZ();
        return;
    }
    
    //	Plastic Strain in XY Direction
    if(var == 36)
    {
        Solout[0] = plasticStrain.XY();
        return;
    }
    
    //	Plastic Strain in XZ Direction
    if(var == 37)
    {
        Solout[0] = plasticStrain.XZ();
        return;
    }
    
    //	Plastic Strain in YZ Direction
    if(var == 38)
    {
        Solout[0] = plasticStrain.YZ();
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
    //	Pore Pressure
    if(var == 40)
    {
        Solout[0] = p;
        return;
    }
    
    //	Darcy's velocity
    if(var == 41)
    {
        Solout[0] = v[0];
        Solout[1] = v[1];
        if (m_Dim == 3)
        {
            Solout[2] = v[2];
        }
        return;
    }
    
    
    //	Porosity
    if(var == 42)
    {
        Solout[0] = m_porosity_0 + (m_alpha * (totalStrain.XX() + totalStrain.YY() + totalStrain.ZZ())) + (m_Se * p);
        return;
    }
    
    
    //	Permeability in XX Direction
    if(var == 43)
    {
        REAL A = 2.0;
        REAL porosity = m_porosity_0 + (m_alpha * (totalStrain.XX() + totalStrain.YY() + totalStrain.ZZ())) + (m_Se * p);
        Solout[0] = m_k_0*pow((porosity/m_porosity_0),A);
        return;
    }
    
    //	Permeability in YY Direction
    if(var == 44)
    {
        REAL A = 2.0;
        REAL porosity = m_porosity_0 + (m_alpha * (totalStrain.XX() + totalStrain.YY() + totalStrain.ZZ())) + (m_Se * p);
        Solout[0] = m_k_0*pow((porosity/m_porosity_0),A);
        return;
    }
    
    //	Permeability in ZZ Direction
    if(var == 45)
    {
        REAL A = 2.0;
        REAL porosity = m_porosity_0 + (m_alpha * (totalStrain.XX() + totalStrain.YY() + totalStrain.ZZ())) + (m_Se * p);
        Solout[0] = m_k_0*pow((porosity/m_porosity_0),A);
        return;
    }
    
    // ************************	Total Stress Variables ************************
    //	Total Stress in XX Direction
    if(var == 46)
    {
        Solout[0] = Sigma.XX();
        return;
    }
    
    //	Total Stress in YY Direction
    if(var == 47)
    {
        Solout[0] = Sigma.YY();
        return;
    }
    
    //	Total Stress in ZZ Direction
    if(var == 48)
    {
        Solout[0] = Sigma.ZZ();
        return;
    }
    
    //	Total Stress in XY Direction
    if(var == 49)
    {
        Solout[0] = Sigma.XY();
        return;
    }
    
    //	Total Stress in XZ Direction
    if(var == 50)
    {
        Solout[0] = Sigma.XZ();
        return;
    }
    
    //	Total Stress in YZ Direction
    if(var == 51)
    {
        Solout[0] = Sigma.YZ();
        return;
    }
    
    // ************************	Stress Ratio Variable ************************
    //	Stress Ratio in XZ Direction
    if(var == 52)
    {
        Solout[0] = Sigma.XX()/Sigma.ZZ();
        return;
    }
    
    // ************************	Yield Surface Variable ************************
    //	Yield Surface 1
    if(var == 53)
    {
        TPZManVector<STATE,3> yieldVal(3,0.);
        plasticloc.Phi(elasticStrain,yieldVal);
        Solout[0] = yieldVal[0];
        return;
    }
    
    //	Yield Surface 2
    if(var == 54)
    {
        TPZManVector<STATE,3> yieldVal(3,0.);
        plasticloc.Phi(elasticStrain,yieldVal);
        Solout[0] = yieldVal[1];
        return;
    }
    
    //	Yield Surface 3
    if(var == 55)
    {
        TPZManVector<STATE,3> yieldVal(3,0.);
        plasticloc.Phi(elasticStrain,yieldVal);
        Solout[0] = yieldVal[2];
        return;
    }
}


// ****************************************************************** plasticity ***********************

template<class T,class TMEM>
void TPZPMRSCouplPoroPlast<T,TMEM>::SetRunPlasticity(bool IsPlasticity)
{
    m_SetRunPlasticity = IsPlasticity;
}


#include "TPZElasticCriterion.h"
//#include "TPZSandlerExtended.h"
//#include "TPZPlasticStepPV.h"
//#include "TPZYCMohrCoulombPV.h"
//#include "TPZSandlerDimaggio.h"

template class TPZPMRSCouplPoroPlast<TPZElasticCriterion , TPZPMRSMemoryPoroPlast>;
//template class TPZPMRSCouplPoroPlast<TPZPlasticStepPV<TPZYCMohrCoulombPV,TPZElasticResponse> , TPZPoroElastoPlasticMem>;
//template class TPZPMRSCouplPoroPlast<TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse> , TPZPoroElastoPlasticMem>;


