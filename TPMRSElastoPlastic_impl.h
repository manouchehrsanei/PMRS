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
    m_simulation_data       = other.m_simulation_data;
    m_dimension             = other.m_dimension;
    m_plastic_integrator    = other.m_plastic_integrator;
}

template <class T, class TMEM>
void TPMRSElastoPlastic<T,TMEM>::Print(std::ostream &out){
    out << Name() << std::endl;
    out << "Material dimension " << m_dimension << std::endl;
    TPZMatWithMem<TMEM>::Print(out);
    m_plastic_integrator.Print(out);
}

template <class T, class TMEM>
void TPMRSElastoPlastic<T,TMEM>::Print(std::ostream &out, const int memory){
    out << Name() << std::endl;
    TPZMatWithMem<TMEM>::Print(out,memory);
}

template <class T, class TMEM>
int TPMRSElastoPlastic<T,TMEM>::VariableIndex(const std::string &name){
    if (!strcmp("ux", name.c_str())) return 0;
    if (!strcmp("uy", name.c_str())) return 1;
    if (!strcmp("uz", name.c_str())) return 2;
    if (!strcmp("sxx", name.c_str())) return 3;
    if (!strcmp("sxy", name.c_str())) return 4;
    if (!strcmp("sxz", name.c_str())) return 5;
    if (!strcmp("syy", name.c_str())) return 6;
    if (!strcmp("syz", name.c_str())) return 7;
    if (!strcmp("szz", name.c_str())) return 8;
    if (!strcmp("exx", name.c_str())) return 9;
    if (!strcmp("exy", name.c_str())) return 10;
    if (!strcmp("exz", name.c_str())) return 11;
    if (!strcmp("eyy", name.c_str())) return 12;
    if (!strcmp("eyz", name.c_str())) return 13;
    if (!strcmp("ezz", name.c_str())) return 14;
    if (!strcmp("epxx", name.c_str())) return 15;
    if (!strcmp("epxy", name.c_str())) return 16;
    if (!strcmp("epxz", name.c_str())) return 17;
    if (!strcmp("epyy", name.c_str())) return 18;
    if (!strcmp("epyz", name.c_str())) return 19;
    if (!strcmp("epzz", name.c_str())) return 20;
    return TPZMatWithMem<TMEM>::VariableIndex(name);
}

template <class T, class TMEM>
void TPMRSElastoPlastic<T,TMEM>::SetSimulationData(TPZSimulationData * simulation_data){
    m_simulation_data = simulation_data;
}

template <class T, class TMEM>
int TPMRSElastoPlastic<T,TMEM>::NSolutionVariables(int var){
    switch(var) {
        case 0:
            return 1; // Scalar
        case 1:
            return 1; // Scalar
        case 2:
            return 1; // Scalar
        case 3:
            return 1; // Scalar
        case 4:
            return 1; // Scalar
        case 5:
            return 1; // Scalar
        case 6:
            return 1; // Scalar
        case 7:
            return 1; // Scalar
        case 8:
            return 1; // Scalar
        case 9:
            return 1; // Scalar
        case 10:
            return 1; // Scalar
        case 11:
            return 1; // Scalar
        case 12:
            return 1; // Scalar
        case 13:
            return 1; // Scalar
        case 14:
            return 1; // Scalar
        case 15:
            return 1; // Scalar
        case 16:
            return 1; // Scalar
        case 17:
            return 1; // Scalar
        case 18:
            return 1; // Scalar
        case 19:
            return 1; // Scalar
        case 20:
            return 1; // Scalar
    }
    return TPZMatWithMem<TMEM>::NSolutionVariables(var);
}

template <class T, class TMEM>
void TPMRSElastoPlastic<T,TMEM>::Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout){
    
    long gp_index = data.intGlobPtIndex;
    TMEM & memory = this->GetMemory().get()->operator[](gp_index);
    Solout.Resize( this->NSolutionVariables(var));
    
    TPZTensor<REAL> epsilon_t = memory.GetPlasticState_n().fEpsT;
    TPZTensor<REAL> epsilon_p = memory.GetPlasticState_n().fEpsP;
    
    switch (var) {
        case 0:
            {
                Solout[0] = memory.Getu_n()[0];
            }
            break;
        case 1:
        {
            Solout[0] = memory.Getu_n()[1];
        }
            break;
        case 2:
        {
            Solout[0] = memory.Getu_n()[2];
        }
            break;
        case 3:
        {
            Solout[0] = memory.GetSigma_n().XX();
        }
            break;
        case 4:
        {
            Solout[0] = memory.GetSigma_n().XY();
        }
            break;
        case 5:
        {
            Solout[0] = memory.GetSigma_n().XZ();
        }
            break;
        case 6:
        {
            Solout[0] = memory.GetSigma_n().YY();
        }
            break;
        case 7:
        {
            Solout[0] = memory.GetSigma_n().YZ();
        }
            break;
        case 8:
        {
            Solout[0] = memory.GetSigma_n().ZZ();
        }
            break;
        case 9:
        {
            Solout[0] = epsilon_t.XX();
        }
            break;
        case 10:
        {
            Solout[0] = epsilon_t.XY();
        }
            break;
        case 11:
        {
            Solout[0] = epsilon_t.XZ();
        }
            break;
        case 12:
        {
            Solout[0] = epsilon_t.YY();
        }
            break;
        case 13:
        {
            Solout[0] = epsilon_t.YZ();
        }
            break;
        case 14:
        {
            Solout[0] = epsilon_t.ZZ();
        }
            break;
        case 15:
        {
            Solout[0] = epsilon_p.XX();
        }
            break;
        case 16:
        {
            Solout[0] = epsilon_p.XY();
        }
            break;
        case 17:
        {
            Solout[0] = epsilon_p.XZ();
        }
            break;
        case 18:
        {
            Solout[0] = epsilon_p.YY();
        }
            break;
        case 19:
        {
            Solout[0] = epsilon_p.YZ();
        }
            break;
        case 20:
        {
            Solout[0] = epsilon_p.ZZ();
        }
            break;
        default:
        {
            std::cout << "TPMRSElastoPlastic<T,TMEM>:: Variable not implemented." << std::endl;
            DebugStop();
        }
            break;
    }
    
}

template <class T, class TMEM>
void TPMRSElastoPlastic<T,TMEM>::Epsilon(TPZMaterialData &data, TPZTensor<REAL> & epsilon_t){
    
    int gp_index = data.intGlobPtIndex;
    TPZTensor<REAL> last_epsilon = this->MemItem(gp_index).GetPlasticState().fEpsT;
    TPZFNMatrix<9,STATE> delta_eps(3,3,0.0), grad_delta_u,grad_delta_u_t;
    TPZFMatrix<REAL>  & dsol_delta_u    = data.dsol[0];
    TPZAxesTools<REAL>::Axes2XYZ(dsol_delta_u, grad_delta_u, data.axes);
    grad_delta_u.Resize(3, 3);
    grad_delta_u.Transpose(&grad_delta_u_t);
    delta_eps = 0.5*(grad_delta_u + grad_delta_u_t);
    
    epsilon_t.XX() = delta_eps(0,0);
    epsilon_t.XY() = delta_eps(0,1);
    epsilon_t.XZ() = delta_eps(0,2);
    epsilon_t.YY() = delta_eps(1,1);
    epsilon_t.YZ() = delta_eps(1,2);
    epsilon_t.ZZ() = delta_eps(2,2);
    
    epsilon_t += last_epsilon;
}

template <class T, class TMEM>
void TPMRSElastoPlastic<T,TMEM>::Sigma(TPZTensor<REAL> & epsilon_t, TPZTensor<REAL> & sigma, TPZFMatrix<REAL> * Dep){
    T plastic_integrator(m_plastic_integrator);
    plastic_integrator.ApplyStrainComputeSigma(epsilon_t,sigma,Dep);
}

template <class T, class TMEM>
void TPMRSElastoPlastic<T,TMEM>::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef){
    
    
    // Getting weight functions
    TPZFMatrix<REAL>  & phi_u     =  data.phi;
    int n_phi_u = phi_u.Rows();
    int first_u  = 0;
    
    TPZFNMatrix<40,REAL> grad_phi_u(3,n_phi_u);
    TPZAxesTools<REAL>::Axes2XYZ(data.dphix, grad_phi_u, data.axes);

    REAL dvdx,dvdy;
    TPZTensor<STATE> epsilon,sigma;
    
    TPZFNMatrix<36,STATE> De(6,6,0.0);
    
    Epsilon(data,epsilon);
    Sigma(epsilon, sigma, &De);
    
    TPZFNMatrix<4,STATE> Deriv(2, 2);
    STATE val;
    for(int iu = 0; iu < n_phi_u; iu++ )
    {
        dvdx = grad_phi_u(0,iu);
        dvdy = grad_phi_u(1,iu);
        
        ef(2*iu + first_u)     +=    weight * (sigma.XX()*dvdx + sigma.XY()*dvdy);    // x direction
        ef(2*iu+1 + first_u)   +=    weight * (sigma.XY()*dvdx + sigma.YY()*dvdy);    // y direction
        
        for(int ju = 0; ju < n_phi_u; ju++)
        {

            for (int ud = 0; ud < 2; ud++) {
                for (int vd = 0; vd < 2; vd++) {
                    Deriv(vd, ud) = grad_phi_u(vd, iu) * grad_phi_u(ud, ju);
                }
            }
            
            val = 2. * De(_XX_, _XX_) * Deriv(0, 0);
            val += De(_XX_, _XY_) * Deriv(0, 1);
            val += 2. * De(_XY_, _XX_) * Deriv(1, 0);
            val += De(_XY_, _XY_) * Deriv(1, 1);
            val *= 0.5;
            ek(2*iu + first_u,2*ju + first_u) += weight * val;
            
            val = De(_XX_, _XY_) * Deriv(0, 0);
            val += 2. * De(_XX_, _YY_) * Deriv(0, 1);
            val += De(_XY_, _XY_) * Deriv(1, 0);
            val += 2. * De(_XY_, _YY_) * Deriv(1, 1);
            val *= 0.5;
            ek(2*iu + first_u,2*ju+1 + first_u) += weight * val;
            
            val = 2. * De(_XY_, _XX_) * Deriv(0, 0);
            val += De(_XY_, _XY_) * Deriv(0, 1);
            val += 2. * De(_YY_, _XX_) * Deriv(1, 0);
            val += De(_YY_, _XY_) * Deriv(1, 1);
            val *= 0.5;
            ek(2*iu+1 + first_u,2*ju + first_u) += weight * val;
            
            val = De(_XY_, _XY_) * Deriv(0, 0);
            val += 2. * De(_XY_, _YY_) * Deriv(0, 1);
            val += De(_YY_, _XY_) * Deriv(1, 0);
            val += 2. * De(_YY_, _YY_) * Deriv(1, 1);
            val *= 0.5;
            ek(2*iu+1 + first_u,2*ju+1 + first_u) += weight * val;

        }
    }
    
    
}

template <class T, class TMEM>
void TPMRSElastoPlastic<T,TMEM>::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef, TPZBndCond &bc){
    
    TPZFMatrix<REAL>  &phiu = data.phi;
    // Getting the solutions and derivatives
    TPZManVector<REAL,2> u = data.sol[0];
    
    int phru = phiu.Rows();
    int in,jn;
    
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
    
    if (m_simulation_data->Get_must_accept_solution_Q()) {
        
        int gp_index = data.intGlobPtIndex;
        
        if (m_simulation_data->GetTransferCurrentToLastQ()) {
            this->MemItem(gp_index).SetPlasticState(this->MemItem(gp_index).GetPlasticState_n());
            this->MemItem(gp_index).SetSigma(this->MemItem(gp_index).GetSigma_n());
            this->MemItem(gp_index).Setu(this->MemItem(gp_index).Getu_n()) ;
            return;
        }
        

        TPZTensor<STATE> epsilon,sigma;
        Epsilon(data,epsilon);
        T plastic_integrator(m_plastic_integrator);
        plastic_integrator.ApplyStrainComputeSigma(epsilon,sigma);
        
        if (m_simulation_data->IsCurrentStateQ()) {
            this->MemItem(gp_index).SetPlasticState_n(plastic_integrator.fN);
            this->MemItem(gp_index).SetSigma_n(sigma);
            
            TPZManVector<STATE,3> delta_u    = data.sol[0];
            TPZManVector<STATE,3> u_n(m_dimension,0.0);
            TPZManVector<STATE,3> u(this->MemItem(gp_index).Getu());
            for (int i = 0; i < m_dimension; i++) {
                u_n[i] = delta_u[i] + u[i];
            }
            this->MemItem(gp_index).Setu_n(u_n);
            
        }else{
            this->MemItem(gp_index).SetPlasticState(plastic_integrator.fN);
            this->MemItem(gp_index).SetSigma(sigma);
            
            TPZManVector<STATE,3> u    = data.sol[0];
            this->MemItem(gp_index).Setu(u);
        }

        
    }else{
        TPZFMatrix<REAL> ek_fake;
        ek_fake.Resize(ef.Rows(),ef.Rows());
        this->Contribute(data, weight, ek_fake, ef);
    }

    
}

template <class T, class TMEM>
void TPMRSElastoPlastic<T,TMEM>::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ef, TPZBndCond &bc){
    TPZFMatrix<REAL> ek_fake;
    ek_fake.Resize(ef.Rows(),ef.Rows());
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








