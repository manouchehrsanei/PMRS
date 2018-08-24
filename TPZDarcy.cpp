//
//  TPZDarcyMem.cpp
//  PZ
//
//  Created by Manouchehr on Agust 24, 2018.
//
//

#include "TPZDarcy.h"
#include <iostream>
#include <string>
#include "pzbndcond.h"
#include "pzaxestools.h"
#include <algorithm>
#include "pzlog.h"
#include "pzfmatrix.h"


#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.TPZDarcy"));
#endif



/** @brief default costructor */
TPZDarcy::TPZDarcy() : TPZMatWithMem<TPZFMatrix<STATE>, TPZDiscontinuousGalerkin >()
{
    m_Dim = 0;
}

/** @brief costructor based on a material id */
TPZDarcy::TPZDarcy(int matid, int dim): TPZMatWithMem<TPZFMatrix<STATE>, TPZDiscontinuousGalerkin >(matid)
{
    m_Dim = dim;
}

/** @brief copy constructor $ */
TPZDarcy::TPZDarcy(const TPZDarcy& other): TPZMatWithMem<TPZFMatrix<STATE>, TPZDiscontinuousGalerkin >(other)
{
    this->m_Dim    = other.m_Dim;
    this->m_SimulationData    = other.m_SimulationData;
}


/** @brief default destructor */
TPZDarcy::~TPZDarcy()
{
}


/** @brief Copy assignemnt operator $ */
TPZDarcy& TPZDarcy::operator = (const TPZDarcy& other)
{
    if (this != & other) // prevent self-assignment
    {
        this->m_Dim    = other.m_Dim;
        this->m_SimulationData    = other.m_SimulationData;
    }
    return *this;
}


/** @brief of compute permeability (Kappa) */
void TPZDarcy::Compute_Kappa(TPZMaterialData &data, REAL &kappa)
{
    kappa = m_k_0;
}


/** @brief of contribute in 2 dimensional */
void TPZDarcy::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE>  &ek, TPZFMatrix<STATE> &ef)
{
    int p_b = 0;
    
    // Getting the space functions
    TPZFMatrix<REAL>        &phip         =   datavec[p_b].phi;
    TPZFMatrix<REAL>        &grad_phi_p   =   datavec[p_b].dphix;
    TPZFNMatrix <9,REAL>    &axes_p	      =	  datavec[p_b].axes;
    
    // Getting the solutions and derivatives
    TPZManVector<REAL,1> p = datavec[p_b].sol[0];
    TPZFNMatrix <6,REAL> dp = datavec[p_b].dsol[0];
    
    TPZFNMatrix<6,REAL> Grad_p(2,1,0.0),Grad_phi_i(2,1,0.0),Grad_phi_j(2,1,0.0);
    Grad_p(0,0) = dp(0,0)*axes_p(0,0)+dp(1,0)*axes_p(1,0);
    Grad_p(1,0) = dp(0,0)*axes_p(0,1)+dp(1,0)*axes_p(1,1);
    
    int nphi_p = phip.Rows();
    
    // Compute permeability
    REAL k = 0.0;
    Compute_Kappa(datavec[p_b], k);
    
    REAL c = (k/m_eta);

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
        
            ef(ip, 0)		+=  weight *  c * dot;
        
        for (int jp = 0; jp < nphi_p; jp++)
        {
            
            Grad_phi_j(0,0) = grad_phi_p(0,jp)*axes_p(0,0)+grad_phi_p(1,jp)*axes_p(1,0);
            Grad_phi_j(1,0) = grad_phi_p(0,jp)*axes_p(0,1)+grad_phi_p(1,jp)*axes_p(1,1);
            
            REAL dot = 0.0;
            for (int i = 0;  i < m_Dim; i++)
            {
                dot += Grad_phi_j(i,0) * Grad_phi_i(i,0);
            }
            
            ek(ip, jp)		+= weight * c * dot;
        }
    }
}


/** @brief of contribute of BC_2D */
void TPZDarcy::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{
    int p_b = 0;

    TPZFMatrix<REAL>  &phip = datavec[p_b].phi;
    TPZManVector<REAL,1> p  = datavec[p_b].sol[0];

    int phrp = phip.Rows();
    short in,jn;

    const REAL BIGNUMBER = TPZMaterial::gBigNumber;

    
    // Boundaries

    // Dirichlet in Pressure
    switch (bc.Type())
    {
        case 0 : // Dp
        {
            REAL v[1];
            v[0] = bc.Val2()(0,0);    //    Pressure

            //    Darcy Equation
            for(in = 0 ; in < phrp; in++)
            {
                //    Contribution for load Vector
                ef(in,0)        += BIGNUMBER*(p[0]-v[0])*phip(in,0)*weight;    // P Pressure

                for (jn = 0 ; jn < phrp; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(in,jn)        += BIGNUMBER*phip(in,0)*phip(jn,0)*weight;    // P Pressure
                }
            }
            break;
        }
            
        // Neumman in Flux
        case 1 : // Nq
        {
            REAL v[1];
            v[0] = bc.Val2()(0,0);    //    Qn

            //    Darcy Equation
            for(in = 0 ; in < phrp; in++)
            {
                //    Normal Flux on neumman boundary
                ef(in,0)    += 1.0 *  weight * v[0] * phip(in,0);    // Qnormal
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


void TPZDarcy::FillDataRequirements(TPZVec<TPZMaterialData> &datavec)

{
    int nref = datavec.size();
    for(int i = 0; i<nref; i++)
    {
        datavec[i].SetAllRequirements(false);
        datavec[i].fNeedsSol = true;
    }
}

void TPZDarcy::FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData> &datavec)
{
    int nref = datavec.size();
    for(int i = 0; i<nref; i++)
    {
        datavec[i].SetAllRequirements(false);
        datavec[i].fNeedsSol = true;
    }
}

void TPZDarcy::Print(std::ostream &out)
{
    out << "Material Name : "               << Name()  << "\n";
    out << "Properties for TPZDarcy: \n";
    out << "\t Initial Permeability   = "	<< m_k_0   << std::endl;
    out << "\t Fluid Viscosity   = "		<< m_eta   << std::endl;
    TPZMaterial::Print(out);
    out << "\n";
    
}

/** Returns the variable index associated with the name */
int TPZDarcy::VariableIndex(const std::string &name)
{
    //	Diffusion Variables
    if(!strcmp("p",name.c_str()))				return	0;
    if(!strcmp("v",name.c_str()))				return	1;
    if(!strcmp("k",name.c_str()))				return	2;
    
    return TPZMaterial::VariableIndex(name);
}

int TPZDarcy::NSolutionVariables(int var)
{
    if(var == 0)	return 1;
    if(var == 1)	return m_Dim;
    if(var == 2)	return 1;

    return TPZMaterial::NSolutionVariables(var);
}

//	Calculate Secondary variables based on ux, uy, Pore pressure and their derivatives
void TPZDarcy::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout)
{
    Solout.Resize( this->NSolutionVariables(var));
    
    int p_b = 0;
    
    // Getting the space functions
    TPZFNMatrix <9,REAL>	&axes_p	=	datavec[p_b].axes;
    
    // Getting the solutions and derivatives
    TPZManVector<REAL,1> p  = datavec[p_b].sol[0];
    TPZFNMatrix <9,REAL> dp = datavec[p_b].dsol[0];
    
    REAL to_Mpa     = 1; // 1.0e-6;
    REAL to_Darcy   = 1; // 1.01327e+12;
    
    // Computing Gradient of the Solution
    TPZFNMatrix<3,REAL> Grad_p(3,1,0.0);
    
    // ************************************** The value of parameters ************************
    
    // ************************	Darcy Variables ************************
    //	Pore Pressure
    if(var == 0)
    {
        Solout[0] = p[0]*to_Mpa;
        return;
    }
    
    //	Darcy's velocity
    if(var == 1)
    {
        Grad_p(0,0) = dp(0,0)*axes_p(0,0)+dp(1,0)*axes_p(1,0); // dp/dx
        Grad_p(1,0) = dp(0,0)*axes_p(0,1)+dp(1,0)*axes_p(1,1); // dp/dy
        
        REAL k = 0.0;
        Compute_Kappa(datavec[p_b], k);
        
        Solout[0] = -(k/m_eta) * Grad_p(0,0);
        Solout[1] = -(k/m_eta) * Grad_p(1,0);
        return;
    }
    
    
    //	Permeability
    if(var == 2)
    {
        REAL k = 0.0;
        Compute_Kappa(datavec[p_b], k);
        Solout[0] = k*to_Darcy;
        return;
    }
    
    
}


