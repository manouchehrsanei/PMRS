//
//  TPMRSElastoPlastic_impl.hpp
//  PMRS
//
//  Created by Omar Durán on 9/11/18.
//

#include "TPMRSElastoPlastic.h"

template <class T, class TMEM>
TPMRSElastoPlastic<T,TMEM>::TPMRSElastoPlastic() : TPZMatWithMem<TMEM>(){
    m_simulation_data = NULL;
    m_dimension = 0;
}

template <class T, class TMEM>
TPMRSElastoPlastic<T,TMEM>::~TPMRSElastoPlastic(){
    
}

template <class T, class TMEM>
TPMRSElastoPlastic<T,TMEM>::TPMRSElastoPlastic(int mate_id) : TPZMatWithMem<TMEM>(mate_id) {
    m_simulation_data = NULL;
    m_dimension = 0;
}

template <class T, class TMEM>
TPMRSElastoPlastic<T,TMEM>::TPMRSElastoPlastic(const TPMRSElastoPlastic & other): TPZMatWithMem<TMEM>(other){
    m_simulation_data   =   other.m_simulation_data;
    m_dimension         =   other.m_dimension;
}

template <class T, class TMEM>
void TPMRSElastoPlastic<T,TMEM>::Print(std::ostream &out){
    out << Name() << std::endl;
    out << "Material dimension " << m_dimension << std::endl;
    TPZMatWithMem<TMEM>::Print(out);
}

template <class T, class TMEM>
void TPMRSElastoPlastic<T,TMEM>::Print(std::ostream &out, const int memory){
    out << Name() << std::endl;
    TPZMatWithMem<TMEM>::Print(out,memory);
}

template <class T, class TMEM>
int TPMRSElastoPlastic<T,TMEM>::VariableIndex(const std::string &name){
    DebugStop();
}

template <class T, class TMEM>
int TPMRSElastoPlastic<T,TMEM>::NSolutionVariables(int var){
    DebugStop();
}

template <class T, class TMEM>
void TPMRSElastoPlastic<T,TMEM>::Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout){
    DebugStop();
}

template <class T, class TMEM>
void TPMRSElastoPlastic<T,TMEM>::Epsilon(TPZMaterialData &data, TPZTensor<REAL> & epsilon_t){
    
    
    int gp_index = data.intGlobPtIndex;
    this->MemItem(gp_index).GetPlasticState().Get
    
    TPZFNMatrix<9,STATE> eps, grad_u(3,3,0.0),grad_u_t;
    TPZFMatrix<REAL>  & dsol_u    = data.dsol[0];
    TPZAxesTools<REAL>::Axes2XYZ(dsol_u, grad_u, data.axes);
    
    grad_u.Transpose(&grad_u_t);
    eps = 0.5*(grad_u + grad_u_t);
    
    epsilon_t.XX() = eps(0,0);
    epsilon_t.XY() = eps(0,1);
    epsilon_t.XZ() = eps(0,2);
    epsilon_t.YY() = eps(1,1);
    epsilon_t.YZ() = eps(1,2);
    epsilon_t.ZZ() = eps(2,2);
    
}

template <class T, class TMEM>
void TPMRSElastoPlastic<T,TMEM>::Sigma(TPZTensor<REAL> & epsilon_t, TPZTensor<REAL> & sigma, TPZFMatrix<REAL> * Dep){
    
}

template <class T, class TMEM>
void TPMRSElastoPlastic<T,TMEM>::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef){
    
    DebugStop();
    
}

template <class T, class TMEM>
void TPMRSElastoPlastic<T,TMEM>::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef, TPZBndCond &bc){
    
    TPZFMatrix<REAL>  &phiu = data.phi;
    // Getting the solutions and derivatives
    TPZManVector<REAL,2> u = data.sol[0];
    
    int phru = phiu.Rows();
    short in,jn;
    
    REAL BigNumber = TPZDiscontinuousGalerkin::gBigNumber;
    
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
                ef(2*in+0,0)      += BigNumber*(u[0] - v[0])*phiu(in,0)*weight;    // X displacement Value
                ef(2*in+1,0)      += BigNumber*(u[1] - v[1])*phiu(in,0)*weight;    // y displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(2*in+0,2*jn+0)    += BigNumber*phiu(in,0)*phiu(jn,0)*weight;    // X displacement
                    ek(2*in+1,2*jn+1)    += BigNumber*phiu(in,0)*phiu(jn,0)*weight;    // Y displacement
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
                ef(2*in,0)        += BigNumber*(d_ux)*phiu(in,0)*weight;    // X displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(2*in,2*jn)        += BigNumber*phiu(in,0)*phiu(jn,0)*weight;    // X displacement
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
                ef(2*in+1,0)      += BigNumber*(d_uy)*phiu(in,0)*weight;    // y displacement
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(2*in+1,2*jn+1)    += BigNumber*phiu(in,0)*phiu(jn,0)*weight;    // Y displacement
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
        
            break;
        }
            
        case 4 : // Ntn_Dp
        {
            
            REAL v[2];
            v[0] = bc.Val2()(0,0);    //    Tn normal traction
            v[1] = bc.Val2()(1,0);    //    Pressure
            
            REAL tn = v[0];
            TPZManVector<REAL,2> n = data.normal;
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
            
        case 5 : // Du_Nq
        {
            
            REAL v[3];
            v[0] = bc.Val2()(0,0);    //    Ux displacement
            v[1] = bc.Val2()(1,0);    //    Uy displacement
            v[2] = bc.Val2()(2,0);    //    Qn
            
            for(in = 0 ; in < phru; in++)
            {
                //    Contribution for load Vector
                ef(2*in+0,0)        += BigNumber*(u[0] - v[0])*phiu(in,0)*weight;    // X displacement Value
                ef(2*in+1,0)      += BigNumber*(u[1] - v[1])*phiu(in,0)*weight;    // y displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(2*in+0,2*jn+0)        += BigNumber*phiu(in,0)*phiu(jn,0)*weight;    // X displacement
                    ek(2*in+1,2*jn+1)    += BigNumber*phiu(in,0)*phiu(jn,0)*weight;    // Y displacement
                }
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
                ef(2*in,0)        += BigNumber*(d_ux)*phiu(in,0)*weight;    // X displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(2*in,2*jn)        += BigNumber*phiu(in,0)*phiu(jn,0)*weight;    // X displacement
                }
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
                ef(2*in+1,0)      += BigNumber*(d_uy)*phiu(in,0)*weight;    // y displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(2*in+1,2*jn+1)    += BigNumber*phiu(in,0)*phiu(jn,0)*weight;    // Y displacement
                }
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
            
            break;
        }
            
        case 9 : // Ntn_Nq
        {
            
            REAL v[2];
            v[0] = bc.Val2()(0,0);    //    Tn normal traction
            v[1] = bc.Val2()(2,0);    //    Qn
            
            REAL tn = v[0];
            TPZManVector<REAL,2> n = data.normal;
            //    Neumann condition for each state variable
            //    Elasticity Equation
            for(in = 0 ; in <phru; in++)
            {
                //    Normal Tension Components on neumman boundary
                ef(2*in+0,0)      += -1.0 * weight * tn * n[0] * phiu(in,0);        //    Tnx
                ef(2*in+1,0)      += -1.0 * weight * tn * n[1] * phiu(in,0);        //    Tny
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

template <class T, class TMEM>
void TPMRSElastoPlastic<T,TMEM>::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ef){
    TPZFMatrix<REAL> ek_fake;
    this->Contribute(data, weight, ek_fake, ef);
}

template <class T, class TMEM>
void TPMRSElastoPlastic<T,TMEM>::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ef, TPZBndCond &bc){
    TPZFMatrix<REAL> ek_fake;
    this->ContributeBC(data, weight, ek_fake, ef, bc);
}

template <class T, class TMEM>
void TPMRSElastoPlastic<T,TMEM>::Write(TPZStream &buf, int withclassid) const {
    TPZMatWithMem<TMEM>::Write(buf, withclassid);
}

template <class T, class TMEM>
void TPMRSElastoPlastic<T,TMEM>::Read(TPZStream &buf, void *context){
    TPZMatWithMem<TMEM>::Read(buf, context);
}

template <class T, class TMEM>
void TPMRSElastoPlastic<T,TMEM>::FillDataRequirements(TPZMaterialData &data){
    TPZMatWithMem<TMEM>::FillDataRequirements(data);
    data.fNeedsSol = true;
    data.fNeedsNormal = false;
    data.fNeedsHSize = false;
    data.fNeedsNeighborCenter = false;
}

template <class T, class TMEM>
void TPMRSElastoPlastic<T,TMEM>::FillBoundaryConditionDataRequirement(int type,TPZMaterialData &data){
    data.fNeedsSol = true;
    data.fNeedsNormal = true;
}








