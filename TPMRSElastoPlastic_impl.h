//
//  TPMRSElastoPlastic_impl.hpp
//  PMRS
//
//  Created by Omar and Manouchehr on 9/11/18.
//

#include "TPMRSElastoPlastic.h"

template <class T, class TMEM>
TPMRSElastoPlastic<T,TMEM>::TPMRSElastoPlastic() : TPZMatWithMem<TMEM>(){
    m_simulation_data = NULL;
    m_dimension       = 0;
}

template <class T, class TMEM>
TPMRSElastoPlastic<T,TMEM>::~TPMRSElastoPlastic(){
    
}

template <class T, class TMEM>
TPMRSElastoPlastic<T,TMEM>::TPMRSElastoPlastic(int mat_id) : TPZMatWithMem<TMEM>(mat_id) {
    m_simulation_data = NULL;
    m_dimension       = 0;
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
    if (!strcmp("ux"  , name.c_str())) return  0;
    if (!strcmp("uy"  , name.c_str())) return  1;
    if (!strcmp("uz"  , name.c_str())) return  2;
    if (!strcmp("sxx" , name.c_str())) return  3;
    if (!strcmp("sxy" , name.c_str())) return  4;
    if (!strcmp("sxz" , name.c_str())) return  5;
    if (!strcmp("syy" , name.c_str())) return  6;
    if (!strcmp("syz" , name.c_str())) return  7;
    if (!strcmp("szz" , name.c_str())) return  8;
    if (!strcmp("exx" , name.c_str())) return  9;
    if (!strcmp("exy" , name.c_str())) return 10;
    if (!strcmp("exz" , name.c_str())) return 11;
    if (!strcmp("eyy" , name.c_str())) return 12;
    if (!strcmp("eyz" , name.c_str())) return 13;
    if (!strcmp("ezz" , name.c_str())) return 14;
    if (!strcmp("epxx", name.c_str())) return 15;
    if (!strcmp("epxy", name.c_str())) return 16;
    if (!strcmp("epxz", name.c_str())) return 17;
    if (!strcmp("epyy", name.c_str())) return 18;
    if (!strcmp("epyz", name.c_str())) return 19;
    if (!strcmp("epzz", name.c_str())) return 20;
    if (!strcmp("u"   , name.c_str())) return 21;
    if (!strcmp("s"   , name.c_str())) return 22;
    if (!strcmp("e"   , name.c_str())) return 23;
    if (!strcmp("ep"  , name.c_str())) return 24;
    if (!strcmp("s_t" , name.c_str())) return 25;
    if (!strcmp("alpha" , name.c_str())) return 26;
    return TPZMatWithMem<TMEM>::VariableIndex(name);
}

template <class T, class TMEM>
void TPMRSElastoPlastic<T,TMEM>::SetSimulationData(TPMRSSimulationData * simulation_data){
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
        case 21:
            return m_dimension; // Vector
        case 22:
            return 9; // Tensor
        case 23:
            return 9; // Tensor
        case 24:
            return 9; // Tensor
        case 25:
            return 9; // Tensor
        case 26:
            return 1; // Tensor
    }
    return TPZMatWithMem<TMEM>::NSolutionVariables(var);
}

template <class T, class TMEM>
void TPMRSElastoPlastic<T,TMEM>::Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout){
    
    long gp_index = data.intGlobPtIndex;
    TMEM & memory = this->GetMemory().get()->operator[](gp_index);
    Solout.Resize(this->NSolutionVariables(var));
    
    TPZTensor<REAL> epsilon_t = memory.GetPlasticState_n().m_eps_t;
    TPZTensor<REAL> epsilon_p = memory.GetPlasticState_n().m_eps_p;
    
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
        case 21:
        {
            for(int i = 0; i < m_dimension; i++){
                Solout[i] = memory.Getu_n()[i];
            }
        }
            break;
        case 22:
        {
            Solout[0] = memory.GetSigma_n().XX();
            Solout[1] = Solout[3] = memory.GetSigma_n().XY();
            Solout[2] = Solout[6] = memory.GetSigma_n().XZ();
            Solout[4] = memory.GetSigma_n().YY();
            Solout[5] = Solout[7] = memory.GetSigma_n().YZ();
            Solout[8] = memory.GetSigma_n().ZZ();
        }
            break;
        case 23:
        {
            Solout[0] = epsilon_t.XX();
            Solout[1] = Solout[3] = epsilon_t.XY();
            Solout[2] = Solout[6] = epsilon_t.XZ();
            Solout[4] = epsilon_t.YY();
            Solout[5] = Solout[7] = epsilon_t.YZ();
            Solout[8] = epsilon_t.ZZ();
        }
            break;
        case 24:
        {
            Solout[0] = epsilon_p.XX();
            Solout[1] = Solout[3] = epsilon_p.XY();
            Solout[2] = Solout[6] = epsilon_p.XZ();
            Solout[4] = epsilon_p.YY();
            Solout[5] = Solout[7] = epsilon_p.YZ();
            Solout[8] = epsilon_p.ZZ();
        }
            break;
        case 25:
        {
            Solout[0] = memory.GetSigma_n().XX() - memory.Alpha()*memory.p_n();
            Solout[1] = Solout[3] = memory.GetSigma_n().XY();
            Solout[2] = Solout[6] = memory.GetSigma_n().XZ();
            Solout[4] = memory.GetSigma_n().YY() - memory.Alpha()*memory.p_n();
            Solout[5] = Solout[7] = memory.GetSigma_n().YZ();
            Solout[8] = memory.GetSigma_n().ZZ() - memory.Alpha()*memory.p_n();
        }
            break;
        case 26:
        {
            Solout[0] = memory.Alpha();
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
    TPZTensor<REAL> last_epsilon;
    if (m_simulation_data->Get_must_use_sub_stepping_Q()) {
        last_epsilon = this->MemItem(gp_index).GetPlasticStateSubStep().m_eps_t;
    }else{
        last_epsilon = this->MemItem(gp_index).GetPlasticState().m_eps_t;
    }
    TPZFNMatrix<9,STATE> delta_eps(3,3,0.0), grad_delta_u, grad_delta_u_t;
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

//#define MEuler_Q
#define Secant_Q

template <class T, class TMEM>
void TPMRSElastoPlastic<T,TMEM>::Sigma(TPZMaterialData &data, TPZTensor<REAL> & epsilon_t, TPZTensor<REAL> & sigma, TPZFMatrix<REAL> * Dep){
    
#ifdef MEuler_Q
    /// https://www.newcastle.edu.au/__data/assets/pdf_file/0017/22463/12_Substepping-schemes-for-the-numerical-integration-of-elastoplastic-stress-strain-relations.pdf
     int gp_index = data.intGlobPtIndex;
    
//    TPZFNMatrix<6,STATE> delta_eps_t_vec(6,1,0.0),sigma_vec(6,1,0.0);
//    TPZFNMatrix<6,STATE> delta_sigma_1_vec, delta_sigma_2_vec;
//
//    /// step one
//    TPZFNMatrix<36,REAL> & Dep_1 = this->MemItem(gp_index).Dep();
//    TPZTensor<REAL> last_sigma = this->MemItem(gp_index).GetSigma();
//
//    /// Last plastic_strain state
//    TPZPlasticState<REAL> & plastic_strain = this->MemItem(gp_index).GetPlasticState();
//    TPZTensor<REAL> last_epsilon_t = plastic_strain.m_eps_t;
//    TPZTensor<REAL> delta_epsilon_t = epsilon_t - last_epsilon_t;
//
//    delta_epsilon_t.CopyTo(delta_eps_t_vec);
//    Dep_1.Multiply(delta_eps_t_vec, delta_sigma_1_vec);
//
//    TPZTensor<REAL> sigma_star, eps_e;
//    sigma_star.CopyFrom(delta_sigma_1_vec);
//    sigma_star += last_sigma;
//
//    /// step two
//    m_plastic_integrator.GetElasticResponse().ComputeStrain(sigma_star, eps_e);
//    T plastic_integrator(m_plastic_integrator);
//    plastic_integrator.fN.m_eps_p = epsilon_t - eps_e;
//    plastic_integrator.ApplyStrainComputeSigma(epsilon_t,sigma,Dep);
//
//    Dep->Multiply(delta_eps_t_vec, delta_sigma_2_vec);
//    last_sigma.CopyTo(sigma_vec);
//    sigma_vec += 0.5*(delta_sigma_1_vec + delta_sigma_2_vec);
//    Dep->operator*=(0.5);
//    Dep->operator+=(0.5*Dep_1);
//    sigma.CopyFrom(sigma_vec);
    
    TPZFNMatrix<6,STATE> delta_eps_t_vec(6,1,0.0),sigma_vec(6,1,0.0);
    TPZFNMatrix<6,STATE> delta_sigma_star_vec, delta_sigma_2_vec;
    
    /// step one
    T plastic_integrator_1(m_plastic_integrator);
    TPZFNMatrix<36,REAL> Dep_star(6,6);
    TPZTensor<REAL> sigma_star, eps_e;
    plastic_integrator_1.ApplyStrainComputeSigma(epsilon_t,sigma_star,&Dep_star);
//    TPZFNMatrix<36,REAL> & Dep_1 = this->MemItem(gp_index).Dep();
    TPZTensor<REAL> last_sigma = this->MemItem(gp_index).GetSigma();
    
    /// Last plastic_strain state
    TPZPlasticState<REAL> & plastic_strain = this->MemItem(gp_index).GetPlasticState();
    TPZTensor<REAL> last_epsilon_t = plastic_strain.m_eps_t;
    TPZTensor<REAL> delta_epsilon_t = epsilon_t - last_epsilon_t;
    
    delta_epsilon_t.CopyTo(delta_eps_t_vec);
    Dep_star.Multiply(delta_eps_t_vec, delta_sigma_star_vec);
    
//
//    sigma_star.CopyFrom(delta_sigma_1_vec);
//    sigma_star += last_sigma;
    
    /// step two
    m_plastic_integrator.GetElasticResponse().ComputeStrain(sigma_star, eps_e);
    T plastic_integrator(m_plastic_integrator);
    plastic_integrator.fN.m_eps_p = epsilon_t - eps_e;
    plastic_integrator.ApplyStrainComputeSigma(epsilon_t,sigma,Dep);
    
    Dep->Multiply(delta_eps_t_vec, delta_sigma_2_vec);
    last_sigma.CopyTo(sigma_vec);
    sigma_vec += 0.5*(delta_sigma_star_vec + delta_sigma_2_vec);
    Dep->operator*=(0.5);
    Dep->operator+=(0.5*Dep_star);
//    sigma.CopyFrom(sigma_vec);
    
#else
    
    bool SecantQ =  m_simulation_data->Get_is_secant_geomechanics_Q();
    int gp_index = data.intGlobPtIndex;
    T plastic_integrator(m_plastic_integrator);
    
    if(SecantQ)
    {
        plastic_integrator.SetState(this->MemItem(gp_index).GetPlasticState());
        plastic_integrator.ApplyStrainComputeSigma(epsilon_t,sigma);
        
        /// Linear elastic tangent.
        {
            const REAL lambda = plastic_integrator.fER.Lambda();
            const REAL mu = plastic_integrator.fER.G();
            
            // Linha 0
            Dep->PutVal(_XX_,_XX_, lambda + 2. * mu);
            Dep->PutVal(_XX_,_YY_, lambda);
            Dep->PutVal(_XX_,_ZZ_, lambda);
            
            // Linha 1
            Dep->PutVal(_XY_,_XY_, 2. * mu);
            
            // Linha 2
            Dep->PutVal(_XZ_,_XZ_, 2. * mu);
            
            // Linha 3
            Dep->PutVal(_YY_,_XX_, lambda);
            Dep->PutVal(_YY_,_YY_, lambda + 2. * mu);
            Dep->PutVal(_YY_,_ZZ_, lambda);
            
            // Linha 4
            Dep->PutVal(_YZ_,_YZ_, 2. * mu);
            
            // Linha 5
            Dep->PutVal(_ZZ_,_XX_, lambda);
            Dep->PutVal(_ZZ_,_YY_, lambda);
            Dep->PutVal(_ZZ_,_ZZ_, lambda + 2. * mu);
        }
    }
    else{
        plastic_integrator.SetState(this->MemItem(gp_index).GetPlasticState());
        plastic_integrator.ApplyStrainComputeSigma(epsilon_t,sigma,Dep);
    }
    
#endif
    
}

template <class T, class TMEM>
void TPMRSElastoPlastic<T,TMEM>::Contribute_Biot_Stress(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ef){
    
    // Getting weight functions
    TPZFMatrix<REAL>  & phi_u     =  data.phi;
    int n_phi_u = phi_u.Rows();
    int first_u  = 0;
    
    TPZFNMatrix<40,REAL> grad_phi_u(3,n_phi_u);
    TPZAxesTools<REAL>::Axes2XYZ(data.dphix, grad_phi_u, data.axes);
    
    REAL dvdx,dvdy;
    
    int gp_index = data.intGlobPtIndex;
    REAL p_0 = this->MemItem(gp_index).p_0();
    REAL p_n = this->MemItem(gp_index).p_n();
    REAL alpha = this->MemItem(gp_index).Alpha();
    
    if (m_dimension == 2) {
        for(int iu = 0; iu < n_phi_u; iu++ )
        {
            dvdx = grad_phi_u(0,iu);
            dvdy = grad_phi_u(1,iu);
            
            ef(m_dimension*iu+0 + first_u)   +=    -1.0 * weight * (alpha*(p_n-p_0)*dvdx);    // x direction
            ef(m_dimension*iu+1 + first_u)   +=    -1.0 * weight * (alpha*(p_n-p_0)*dvdy);    // y direction
        }
    }
    else{
        REAL dvdz;
        for(int iu = 0; iu < n_phi_u; iu++ )
        {
            dvdx = grad_phi_u(0,iu);
            dvdy = grad_phi_u(1,iu);
            dvdz = grad_phi_u(2,iu);

            ef(m_dimension*iu+0 + first_u)   +=    -1.0 * weight * (alpha*(p_n-p_0)*dvdx);    // x direction
            ef(m_dimension*iu+1 + first_u)   +=    -1.0 * weight * (alpha*(p_n-p_0)*dvdy);    // y direction
            ef(m_dimension*iu+2 + first_u)   +=    -1.0 * weight * (alpha*(p_n-p_0)*dvdz);    // z direction
        }
    }
    
}

template <class T, class TMEM>
void TPMRSElastoPlastic<T,TMEM>::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef){
    
    // Getting weight functions
    TPZFMatrix<REAL>  & phi_u     =  data.phi;
    int n_phi_u = phi_u.Rows();
    int first_u  = 0;
    
    TPZFNMatrix<40,REAL> grad_phi_u(3,n_phi_u);
    TPZAxesTools<REAL>::Axes2XYZ(data.dphix, grad_phi_u, data.axes);

    REAL dvdx,dvdy,dvdz;
    TPZTensor<STATE> epsilon,sigma;
    
    TPZFNMatrix<36,STATE> De(6,6,0.0);

    /// Get initial effective stress state
    int gp_index = data.intGlobPtIndex;
    TPZTensor<REAL> & sigma_0 = this->MemItem(gp_index).GetSigma_0();

    /// Get current effective stress state
    Epsilon(data,epsilon);
    Sigma(data, epsilon, sigma, &De);
    
    TPZFNMatrix<9,STATE> Deriv(m_dimension, m_dimension);
    STATE val;
    
    sigma -= sigma_0;// Applying the intial prestress
    if (m_dimension == 2) { /// Plane strain conditions.
        
        for(int iu = 0; iu < n_phi_u; iu++ )
        {
            dvdx = grad_phi_u(0,iu);
            dvdy = grad_phi_u(1,iu);
            
            ef(m_dimension*iu+0 + first_u)   +=    weight * (sigma.XX()*dvdx + sigma.XY()*dvdy);    // x direction
            ef(m_dimension*iu+1 + first_u)   +=    weight * (sigma.XY()*dvdx + sigma.YY()*dvdy);    // y direction
            
            for(int ju = 0; ju < n_phi_u; ju++)
            {
                
                for (int ud = 0; ud < m_dimension; ud++) {
                    for (int vd = 0; vd < m_dimension; vd++) {
                        Deriv(vd, ud) = grad_phi_u(vd, iu) * grad_phi_u(ud, ju);
                    }
                }
                
                val  = 2. * De(_XX_, _XX_) * Deriv(0, 0);//dvdx*dudx
                val +=      De(_XX_, _XY_) * Deriv(0, 1);//dvdx*dudy
                val += 2. * De(_XY_, _XX_) * Deriv(1, 0);//dvdy*dudx
                val +=      De(_XY_, _XY_) * Deriv(1, 1);//dvdy*dudy
                val *= 0.5;
                ek(m_dimension*iu+0 + first_u, m_dimension*ju+0 + first_u) += weight * val;
                
                val  =      De(_XX_, _XY_) * Deriv(0, 0);
                val += 2. * De(_XX_, _YY_) * Deriv(0, 1);
                val +=      De(_XY_, _XY_) * Deriv(1, 0);
                val += 2. * De(_XY_, _YY_) * Deriv(1, 1);
                val *= 0.5;
                ek(m_dimension*iu+0 + first_u, m_dimension*ju+1 + first_u) += weight * val;
                
                val  = 2. * De(_XY_, _XX_) * Deriv(0, 0);
                val +=      De(_XY_, _XY_) * Deriv(0, 1);
                val += 2. * De(_YY_, _XX_) * Deriv(1, 0);
                val +=      De(_YY_, _XY_) * Deriv(1, 1);
                val *= 0.5;
                ek(m_dimension*iu+1 + first_u, m_dimension*ju+0 + first_u) += weight * val;
                
                val  =      De(_XY_, _XY_) * Deriv(0, 0);
                val += 2. * De(_XY_, _YY_) * Deriv(0, 1);
                val +=      De(_YY_, _XY_) * Deriv(1, 0);
                val += 2. * De(_YY_, _YY_) * Deriv(1, 1);
                val *= 0.5;
                ek(m_dimension*iu+1 + first_u, m_dimension*ju+1 + first_u) += weight * val;
                
            }
        }
    }else{

        
        for(int iu = 0; iu < n_phi_u; iu++ )
        {
            dvdx = grad_phi_u(0,iu);
            dvdy = grad_phi_u(1,iu);
            dvdz = grad_phi_u(2,iu);
            
            ef(m_dimension*iu+0 + first_u)   +=    weight * (sigma.XX()*dvdx + sigma.XY()*dvdy + sigma.XZ()*dvdz);    // x direction
            ef(m_dimension*iu+1 + first_u)   +=    weight * (sigma.XY()*dvdx + sigma.YY()*dvdy + sigma.YZ()*dvdz);    // y direction
            ef(m_dimension*iu+2 + first_u)   +=    weight * (sigma.XZ()*dvdx + sigma.YZ()*dvdy + sigma.ZZ()*dvdz);    // z direction
            
            for(int ju = 0; ju < n_phi_u; ju++)
            {
                
                for (int ud = 0; ud < m_dimension; ud++) {
                    for (int vd = 0; vd < m_dimension; vd++) {
                        Deriv(vd, ud) = grad_phi_u(vd, iu) * grad_phi_u(ud, ju);
                    }
                }
                
                val  = 2. * De(_XX_,_XX_) * Deriv(0,0);//dvdx*dudx
                val +=      De(_XX_,_XY_) * Deriv(0,1);//dvdx*dudy
                val +=      De(_XX_,_XZ_) * Deriv(0,2);//dvdx*dudz
                val += 2. * De(_XY_,_XX_) * Deriv(1,0);//dvdy*dudx
                val +=      De(_XY_,_XY_) * Deriv(1,1);//dvdy*dudy
                val +=      De(_XY_,_XZ_) * Deriv(1,2);//dvdy*dudz
                val += 2. * De(_XZ_,_XX_) * Deriv(2,0);//dvdz*dudx
                val +=      De(_XZ_,_XY_) * Deriv(2,1);//dvdz*dudy
                val +=      De(_XZ_,_XZ_) * Deriv(2,2);//dvdz*dudz
                val *= 0.5;
                ek(m_dimension*iu+0 + first_u, m_dimension*ju+0 + first_u) += weight * val;
                
                val  =      De(_XX_,_XY_) * Deriv(0,0);
                val += 2. * De(_XX_,_YY_) * Deriv(0,1);
                val +=      De(_XX_,_YZ_) * Deriv(0,2);
                val +=      De(_XY_,_XY_) * Deriv(1,0);
                val += 2. * De(_XY_,_YY_) * Deriv(1,1);
                val +=      De(_XY_,_YZ_) * Deriv(1,2);
                val +=      De(_XZ_,_XY_) * Deriv(2,0);
                val += 2. * De(_XZ_,_YY_) * Deriv(2,1);
                val +=      De(_XZ_,_YZ_) * Deriv(2,2);
                val *= 0.5;
                ek(m_dimension*iu+0 + first_u, m_dimension*ju+1 + first_u) += weight * val;
                
                val  =      De(_XX_,_XZ_) * Deriv(0,0);
                val +=      De(_XX_,_YZ_) * Deriv(0,1);
                val += 2. * De(_XX_,_ZZ_) * Deriv(0,2);
                val +=      De(_XY_,_XZ_) * Deriv(1,0);
                val +=      De(_XY_,_YZ_) * Deriv(1,1);
                val += 2. * De(_XY_,_ZZ_) * Deriv(1,2);
                val +=      De(_XZ_,_XZ_) * Deriv(2,0);
                val +=      De(_XZ_,_YZ_) * Deriv(2,1);
                val += 2. * De(_XZ_,_ZZ_) * Deriv(2,2);
                val *= 0.5;
                ek(m_dimension*iu+0 + first_u, m_dimension*ju+2 + first_u) += weight * val;
                
                val  = 2. * De(_XY_,_XX_) * Deriv(0,0);
                val +=      De(_XY_,_XY_) * Deriv(0,1);
                val +=      De(_XY_,_XZ_) * Deriv(0,2);
                val += 2. * De(_YY_,_XX_) * Deriv(1,0);
                val +=      De(_YY_,_XY_) * Deriv(1,1);
                val +=      De(_YY_,_XZ_) * Deriv(1,2);
                val += 2. * De(_YZ_,_XX_) * Deriv(2,0);
                val +=      De(_YZ_,_XY_) * Deriv(2,1);
                val +=      De(_YZ_,_XZ_) * Deriv(2,2);
                val *= 0.5;
                ek(m_dimension*iu+1 + first_u, m_dimension*ju+0 + first_u) += weight * val;
                
                val  =      De(_XY_,_XY_) * Deriv(0,0);
                val += 2. * De(_XY_,_YY_) * Deriv(0,1);
                val +=      De(_XY_,_YZ_) * Deriv(0,2);
                val +=      De(_YY_,_XY_) * Deriv(1,0);
                val += 2. * De(_YY_,_YY_) * Deriv(1,1);
                val +=      De(_YY_,_YZ_) * Deriv(1,2);
                val +=      De(_YZ_,_XY_) * Deriv(2,0);
                val += 2. * De(_YZ_,_YY_) * Deriv(2,1);
                val +=      De(_YZ_,_YZ_) * Deriv(2,2);
                val *= 0.5;
                ek(m_dimension*iu+1 + first_u, m_dimension*ju+1 + first_u) += weight * val;
                
                val  =      De(_XY_,_XZ_) * Deriv(0,0);
                val +=      De(_XY_,_YZ_) * Deriv(0,1);
                val += 2. * De(_XY_,_ZZ_) * Deriv(0,2);
                val +=      De(_YY_,_XZ_) * Deriv(1,0);
                val +=      De(_YY_,_YZ_) * Deriv(1,1);
                val += 2. * De(_YY_,_ZZ_) * Deriv(1,2);
                val +=      De(_YZ_,_XZ_) * Deriv(2,0);
                val +=      De(_YZ_,_YZ_) * Deriv(2,1);
                val += 2. * De(_YZ_,_ZZ_) * Deriv(2,2);
                val *= 0.5;
                ek(m_dimension*iu+1 + first_u, m_dimension*ju+2 + first_u) += weight * val;
                
                val  = 2. * De(_XZ_,_XX_) * Deriv(0,0);
                val +=      De(_XZ_,_XY_) * Deriv(0,1);
                val +=      De(_XZ_,_XZ_) * Deriv(0,2);
                val += 2. * De(_YZ_,_XX_) * Deriv(1,0);
                val +=      De(_YZ_,_XY_) * Deriv(1,1);
                val +=      De(_YZ_,_XZ_) * Deriv(1,2);
                val += 2. * De(_ZZ_,_XX_) * Deriv(2,0);
                val +=      De(_ZZ_,_XY_) * Deriv(2,1);
                val +=      De(_ZZ_,_XZ_) * Deriv(2,2);
                val *= 0.5;
                ek(m_dimension*iu+2 + first_u, m_dimension*ju+0 + first_u) += weight * val;
                
                val  =      De(_XZ_,_XY_) * Deriv(0,0);
                val += 2. * De(_XZ_,_YY_) * Deriv(0,1);
                val +=      De(_XZ_,_YZ_) * Deriv(0,2);
                val +=      De(_YZ_,_XY_) * Deriv(1,0);
                val += 2. * De(_YZ_,_YY_) * Deriv(1,1);
                val +=      De(_YZ_,_YZ_) * Deriv(1,2);
                val +=      De(_ZZ_,_XY_) * Deriv(2,0);
                val += 2. * De(_ZZ_,_YY_) * Deriv(2,1);
                val +=      De(_ZZ_,_YZ_) * Deriv(2,2);
                val *= 0.5;
                ek(m_dimension*iu+2 + first_u, m_dimension*ju+1 + first_u) += weight * val;
                
                val  =      De(_XZ_,_XZ_) * Deriv(0,0);
                val +=      De(_XZ_,_YZ_) * Deriv(0,1);
                val += 2. * De(_XZ_,_ZZ_) * Deriv(0,2);
                val +=      De(_YZ_,_XZ_) * Deriv(1,0);
                val +=      De(_YZ_,_YZ_) * Deriv(1,1);
                val += 2. * De(_YZ_,_ZZ_) * Deriv(1,2);
                val +=      De(_ZZ_,_XZ_) * Deriv(2,0);
                val +=      De(_ZZ_,_YZ_) * Deriv(2,1);
                val += 2. * De(_ZZ_,_ZZ_) * Deriv(2,2);
                val *= 0.5;
                ek(m_dimension*iu+2 + first_u, m_dimension*ju+2 + first_u) += weight * val;
                
            }
        }
    }
    
    /// Biot's effect
    this->Contribute_Biot_Stress(data, weight, ef);

}

template <class T, class TMEM>
void TPMRSElastoPlastic<T,TMEM>::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef, TPZBndCond &bc){
    
    
    if(m_dimension == 3){
        this->ContributeBC_3D(data,weight,ek,ef,bc);
        return;
    }
    
    TPZBndCondWithMem<TPMRSElastoPlasticMemory> & bc_with_memory = dynamic_cast<TPZBndCondWithMem<TPMRSElastoPlasticMemory> &>(bc);
    int gp_index = data.intGlobPtIndex;
    if (m_simulation_data->Get_must_accept_solution_Q())
    {
        
        if (m_simulation_data->GetTransferCurrentToLastQ()) {
            if (m_simulation_data->Get_must_use_sub_stepping_Q()) {
                bc_with_memory.MemItem(gp_index).Setu_sub_step(bc_with_memory.MemItem(gp_index).Getu_n());
            }else{
                bc_with_memory.MemItem(gp_index).Setu_sub_step(bc_with_memory.MemItem(gp_index).Getu_n()) ;
                bc_with_memory.MemItem(gp_index).Setu(bc_with_memory.MemItem(gp_index).Getu_n());
                
            }
            return;
        }
        
        
        if (m_simulation_data->IsCurrentStateQ())
        {
            TPZManVector<STATE,3> delta_u    = data.sol[0];
            TPZManVector<STATE,3> u_n(m_dimension,0.0);
            if (m_simulation_data->Get_must_use_sub_stepping_Q()) {
                TPZManVector<STATE,3> u(bc_with_memory.MemItem(gp_index).Getu_sub_step());
                for (int i = 0; i < m_dimension; i++) {
                    u_n[i] = delta_u[i] + u[i];
                }
            }else{
                TPZManVector<STATE,3> u(bc_with_memory.MemItem(gp_index).Getu());
                for (int i = 0; i < m_dimension; i++) {
                    u_n[i] = delta_u[i] + u[i];
                }
            }
            bc_with_memory.MemItem(gp_index).Setu_n(u_n);
            
        }else
        {
            TPZManVector<STATE,3> u    = data.sol[0];
            bc_with_memory.MemItem(gp_index).Setu(u);
        }
    }
    
    TPZFMatrix<REAL>  &phiu = data.phi;
    TPZManVector<STATE,3> delta_u    = data.sol[0];
    TPZManVector<STATE,3> u_n(m_dimension,0.0);
    if (m_simulation_data->Get_must_use_sub_stepping_Q()) {
        TPZManVector<STATE,3> u(bc_with_memory.MemItem(gp_index).Getu_sub_step());
        for (int i = 0; i < m_dimension; i++) {
            u_n[i] = delta_u[i] + u[i];
        }
    }else{
        TPZManVector<STATE,3> u(bc_with_memory.MemItem(gp_index).Getu());
        for (int i = 0; i < m_dimension; i++) {
            u_n[i] = delta_u[i] + u[i];
        }
    }
    
    int phru = phiu.Rows();
    int in,jn;
    
    REAL BigNumber = TPZDiscontinuousGalerkin::gBigNumber;
    
    switch (bc.Type())
    {
        case 2 : // Du
            /// Dirichlet of displacement

        {
            REAL v[2];
            v[0] = bc.Val2()(0,0);    //    Ux displacement
            v[1] = bc.Val2()(1,0);    //    Uy displacement
            
            for(in = 0 ; in < phru; in++)
            {
                ///    Contribution for load Vector
                ef(2*in+0,0)      += BigNumber*(u_n[0] - v[0])*phiu(in,0)*weight;    // X displacement Value
                ef(2*in+1,0)      += BigNumber*(u_n[1] - v[1])*phiu(in,0)*weight;    // Y displacement Value

                
                for (jn = 0 ; jn < phru; jn++)
                {
                    ///    Contribution for Stiffness Matrix
                    ek(2*in+0,2*jn+0)    += BigNumber*phiu(in,0)*phiu(jn,0)*weight;    /// X displacement
                    ek(2*in+1,2*jn+1)    += BigNumber*phiu(in,0)*phiu(jn,0)*weight;    // Y displacement

                }
            }
            
            break;
        }
            
        
        case 3 : /// Dux
            /// Dirichlet in x direction of displacement

        {
            REAL v[1];
            v[0] = bc.Val2()(0,0);    //    Ux displacement
            
            for(in = 0 ; in < phru; in++)
            {
                ///    Contribution for load Vector
                ef(2*in,0)        += BigNumber*(u_n[0] - v[0])*phiu(in,0)*weight;    // X displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    ///    Contribution for Stiffness Matrix
                    ek(2*in,2*jn)        += BigNumber*phiu(in,0)*phiu(jn,0)*weight;    // X displacement
                }
            }
    
            break;
        }
            
        case 4 : /// Duy
            /// Dirichlet in y direction of displacement

        {
            REAL v[1];
            v[0] = bc.Val2()(0,0);    //    Uy displacement
            
            for(in = 0 ; in < phru; in++)
            {
                ///   Contribution for load Vector
                ef(2*in+1,0)      += BigNumber*(u_n[1] - v[0])*phiu(in,0)*weight;    // Y displacement
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    ///    Contribution for Stiffness Matrix
                    ek(2*in+1,2*jn+1)    += BigNumber*phiu(in,0)*phiu(jn,0)*weight;    // Y displacement
                }
            }

            break;
        }
            
        case 5 : /// Nt
            /// Neumann of traction

        {
            REAL v[2];
            v[0] = bc.Val2()(0,0);    //    Tnx
            v[1] = bc.Val2()(1,0);    //    Tny

            ///    Neumann condition for each state variable
            ///    Elasticity Equation
            for(in = 0 ; in <phru; in++)
            {
                ///    Normal Tension Components on neumman boundary
                ef(2*in+0,0)    += -1.0 * weight * v[0] * phiu(in,0);        //    Tnx
                ef(2*in+1,0)    += -1.0 * weight * v[1] * phiu(in,0);        //    Tny
            }
        
            break;
        }
            
        case 6 : /// Ntn
            /// Neumann of traction

        {
            REAL v[1];
            v[0] = bc.Val2()(0,0);    //    Tn normal traction
            
            REAL tn = v[0];
            TPZManVector<REAL,3> n = data.normal;
            ///    Neumann condition for each state variable
            ///    Elasticity Equation
            for(in = 0 ; in <phru; in++)
            {
                ///   Normal Tension Components on neumman boundary
                ef(2*in+0,0)    += -1.0 * weight * tn * n[0] * phiu(in,0);        //    Tnx
                ef(2*in+1,0)    += -1.0 * weight * tn * n[1] * phiu(in,0);        //    Tny
            }
        
            break;
        }
            
        case 7 : /// NSn
            /// Neumann of normal sigma
            
        {
            TPZManVector<REAL,3> n = data.normal;
            TPZFNMatrix<9,REAL> sigma(3,3),normal(3,1),t(3,1);
            t.Zero();
            sigma(0,0)              = bc.Val2()(0,0);    //    sigma_xx
            sigma(0,1) = sigma(1,0) = bc.Val2()(1,0);    //    sigma_xy
            sigma(0,2) = sigma(2,0) = bc.Val2()(2,0);    //    sigma_xz
            sigma(1,1)              = bc.Val2()(3,0);    //    sigma_yy
            sigma(1,2) = sigma(2,1) = bc.Val2()(4,0);    //    sigma_yz
            sigma(2,2)              = bc.Val2()(5,0);    //    sigma_zz
            normal(0,0) = n[0];
            normal(1,0) = n[1];
            normal(2,0) = n[2];
            sigma.Multiply(normal, t);
            
            ///    Neumann condition for each state variable
            ///    Elasticity Equation
            for(in = 0 ; in <phru; in++)
            {
                ///   Normal Tension Components on neumman boundary
                ef(2*in+0,0)    += -1.0 * weight * t(0,0) * phiu(in,0);        //    Tnx
                ef(2*in+1,0)    += -1.0 * weight * t(1,0) * phiu(in,0);        //    Tny
            }
            
            break;
        }
            
        case 8 : // Dun
            /// Dirichlet of normal displacement
        {
            TPZManVector<REAL,3> n = data.normal;
            REAL v[1];
            v[0] = bc.Val2()(0,0);    //    Un displacement
            
            for(in = 0 ; in < phru; in++)
            {
                ///    Contribution for load Vector
                ef(2*in+0,0)      += BigNumber*((u_n[0] - v[0])*n[0])*phiu(in,0)*weight;    // X displacement Value
                ef(2*in+1,0)      += BigNumber*((u_n[1] - v[0])*n[1])*phiu(in,0)*weight;    // Y displacement Value
                
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    ///    Contribution for Stiffness Matrix
                    ek(2*in+0,2*jn+0)    += BigNumber*n[0]*phiu(in,0)*phiu(jn,0)*weight;    // X displacement
                    ek(2*in+1,2*jn+1)    += BigNumber*n[1]*phiu(in,0)*phiu(jn,0)*weight;    // Y displacement
                    
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


template <class T, class TMEM>
void TPMRSElastoPlastic<T,TMEM>::ContributeBC_3D(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef, TPZBndCond &bc){
    
    TPZBndCondWithMem<TPMRSElastoPlasticMemory> & bc_with_memory = dynamic_cast<TPZBndCondWithMem<TPMRSElastoPlasticMemory> &>(bc);
    int gp_index = data.intGlobPtIndex;
    if (m_simulation_data->Get_must_accept_solution_Q()) {
        
        if (m_simulation_data->GetTransferCurrentToLastQ()) {
            if (m_simulation_data->Get_must_use_sub_stepping_Q()) {
                bc_with_memory.MemItem(gp_index).Setu_sub_step(bc_with_memory.MemItem(gp_index).Getu_n());
            }else{
                bc_with_memory.MemItem(gp_index).Setu_sub_step(bc_with_memory.MemItem(gp_index).Getu_n()) ;
                bc_with_memory.MemItem(gp_index).Setu(bc_with_memory.MemItem(gp_index).Getu_n());
                
            }
            return;
        }
        
        
        if (m_simulation_data->IsCurrentStateQ()) {
            
            TPZManVector<STATE,3> delta_u    = data.sol[0];
            TPZManVector<STATE,3> u_n(m_dimension,0.0);
            if (m_simulation_data->Get_must_use_sub_stepping_Q()) {
                TPZManVector<STATE,3> u(bc_with_memory.MemItem(gp_index).Getu_sub_step());
                for (int i = 0; i < m_dimension; i++) {
                    u_n[i] = delta_u[i] + u[i];
                }
            }else{
                TPZManVector<STATE,3> u(bc_with_memory.MemItem(gp_index).Getu());
                for (int i = 0; i < m_dimension; i++) {
                    u_n[i] = delta_u[i] + u[i];
                }
            }
            bc_with_memory.MemItem(gp_index).Setu_n(u_n);
            
        }else{
            TPZManVector<STATE,3> u    = data.sol[0];
            bc_with_memory.MemItem(gp_index).Setu(u);
        }
        
        
    }
    
    TPZFMatrix<REAL>  &phiu = data.phi;
    TPZManVector<STATE,3> delta_u    = data.sol[0];
    TPZManVector<STATE,3> u_n(m_dimension,0.0);
    if (m_simulation_data->Get_must_use_sub_stepping_Q()) {
        TPZManVector<STATE,3> u(bc_with_memory.MemItem(gp_index).Getu_sub_step());
        for (int i = 0; i < m_dimension; i++) {
            u_n[i] = delta_u[i] + u[i];
        }
    }else{
        TPZManVector<STATE,3> u(bc_with_memory.MemItem(gp_index).Getu());
        for (int i = 0; i < m_dimension; i++) {
            u_n[i] = delta_u[i] + u[i];
        }
    }
    
    int phru = phiu.Rows();
    int in,jn;
    
    REAL BigNumber = TPZDiscontinuousGalerkin::gBigNumber;
    BigNumber *= 0.001;
    
    switch (bc.Type())
    {
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
                ef(3*in+0,0)      += BigNumber*(u_n[0] - v[0])*phiu(in,0)*weight;    // X displacement Value
                ef(3*in+1,0)      += BigNumber*(u_n[1] - v[1])*phiu(in,0)*weight;    // Y displacement Value
                ef(3*in+2,0)      += BigNumber*(u_n[2] - v[2])*phiu(in,0)*weight;    // Z displacement Value
                
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(3*in+0,3*jn+0)    += BigNumber*phiu(in,0)*phiu(jn,0)*weight;    // X displacement
                    ek(3*in+1,3*jn+1)    += BigNumber*phiu(in,0)*phiu(jn,0)*weight;    // Y displacement
                    ek(3*in+2,3*jn+2)    += BigNumber*phiu(in,0)*phiu(jn,0)*weight;    // Z displacement
                    
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
                ef(3*in,0)        += BigNumber*(u_n[0] - v[0])*phiu(in,0)*weight;    // X displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(3*in,3*jn)        += BigNumber*phiu(in,0)*phiu(jn,0)*weight;    // X displacement
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
                ef(3*in+1,0)      += BigNumber*(u_n[1] - v[0])*phiu(in,0)*weight;    // Y displacement
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(3*in+1,3*jn+1)    += BigNumber*phiu(in,0)*phiu(jn,0)*weight;    // Y displacement
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
            v[2] = bc.Val2()(2,0);    //    Tny
            
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
            TPZManVector<REAL,3> n = data.normal;
            //    Neumann condition for each state variable
            //    Elasticity Equation
            for(in = 0 ; in <phru; in++)
            {
                //    Normal Tension Components on neumman boundary
                ef(3*in+0,0)    += -1.0 * weight * tn * n[0] * phiu(in,0);        //    Tnx
                ef(3*in+1,0)    += -1.0 * weight * tn * n[1] * phiu(in,0);        //    Tny
                ef(3*in+2,0)    += -1.0 * weight * tn * n[2] * phiu(in,0);        //    Tnz
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
                ef(3*in+2,0)		+= BigNumber*(u_n[2] - v[0])*phiu(in,0)*weight;	// Z displacement Value
                
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in+2,3*jn+2)		+= BigNumber*phiu(in,0)*phiu(jn,0)*weight;	// Z displacement
                    
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
                ef(3*in+0,0)	+= BigNumber*(u_n[0] - v[0])*phiu(in,0)*weight;	// X displacement Value
                ef(3*in+1,0)	+= BigNumber*(u_n[1] - v[1])*phiu(in,0)*weight;	// Y displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in+0,3*jn+0)	+= BigNumber*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
                    ek(3*in+1,3*jn+1)	+= BigNumber*phiu(in,0)*phiu(jn,0)*weight;	// Y displacement
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
                ef(3*in+0,0)	+= BigNumber*(u_n[0] - v[0])*phiu(in,0)*weight;	// X displacement Value
                ef(3*in+2,0)	+= BigNumber*(u_n[2] - v[1])*phiu(in,0)*weight;	// Z displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in+0,3*jn+0)	+= BigNumber*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
                    ek(3*in+2,3*jn+2)	+= BigNumber*phiu(in,0)*phiu(jn,0)*weight;	// Z displacement
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
                ef(3*in+1,0)	+= BigNumber*(u_n[1] - v[0])*phiu(in,0)*weight;	// Y displacement Value
                ef(3*in+2,0)	+= BigNumber*(u_n[2] - v[1])*phiu(in,0)*weight;	// Z displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in+1,3*jn+1)	+= BigNumber*phiu(in,0)*phiu(jn,0)*weight;	// Y displacement
                    ek(3*in+2,3*jn+2)	+= BigNumber*phiu(in,0)*phiu(jn,0)*weight;	// Z displacement
                }
            }
            
            break;
        }
        case 11 : /// NSn
            /// Neumann of normal sigma
            
        {
            TPZManVector<REAL,3> n = data.normal;
            TPZFNMatrix<9,REAL> sigma(3,3),normal(3,1),t(3,1);
            t.Zero();
            sigma(0,0)              = bc.Val2()(0,0);    //    sigma_xx
            sigma(0,1) = sigma(1,0) = bc.Val2()(1,0);    //    sigma_xy
            sigma(0,2) = sigma(2,0) = bc.Val2()(2,0);    //    sigma_xz
            sigma(1,1)              = bc.Val2()(3,0);    //    sigma_yy
            sigma(1,2) = sigma(2,1) = bc.Val2()(4,0);    //    sigma_yz
            sigma(2,2)              = bc.Val2()(5,0);    //    sigma_zz
            normal(0,0) = n[0];
            normal(1,0) = n[1];
            normal(2,0) = n[2];
            sigma.Multiply(normal, t);
            
            ///    Neumann condition for each state variable
            ///    Elasticity Equation
            for(in = 0 ; in <phru; in++)
            {
                ///   Normal Tension Components on neumman boundary
                ef(3*in+0,0)    += -1.0 * weight * t(0,0) * phiu(in,0);        //    Tnx
                ef(3*in+1,0)    += -1.0 * weight * t(1,0) * phiu(in,0);        //    Tny
                ef(3*in+2,0)    += -1.0 * weight * t(2,0) * phiu(in,0);        //    Tnz
            }
            
            break;
        }
            
        case 12 : // Dun
            /// Dirichlet of normal displacement
            
        {
            TPZManVector<REAL,3> n = data.normal;
            REAL v[1];
            v[0] = bc.Val2()(0,0);    //    Un displacement
            
            for(in = 0 ; in < phru; in++)
            {
                ///    Contribution for load Vector
                ef(3*in+0,0)      += BigNumber*((u_n[0] - v[0])*n[0])*phiu(in,0)*weight;    // X displacement Value
                ef(3*in+1,0)      += BigNumber*((u_n[1] - v[0])*n[1])*phiu(in,0)*weight;    // Y displacement Value
                ef(3*in+2,0)      += BigNumber*((u_n[2] - v[0])*n[2])*phiu(in,0)*weight;    // Y displacement Value
                
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    ///    Contribution for Stiffness Matrix
                    ek(3*in+0,3*jn+0)    += BigNumber*n[0]*phiu(in,0)*phiu(jn,0)*weight;    // X displacement
                    ek(3*in+1,3*jn+1)    += BigNumber*n[1]*phiu(in,0)*phiu(jn,0)*weight;    // Y displacement
                    ek(3*in+2,3*jn+2)    += BigNumber*n[2]*phiu(in,0)*phiu(jn,0)*weight;    // Z displacement
                    
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


template <class T, class TMEM>
void TPMRSElastoPlastic<T,TMEM>::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ef){
    
    if (m_simulation_data->Get_must_accept_solution_Q()) {
        
        int gp_index = data.intGlobPtIndex;
        
        if (m_simulation_data->GetTransferCurrentToLastQ()) {
            
            if (m_simulation_data->Get_must_use_sub_stepping_Q()) {
                this->MemItem(gp_index).SetPlasticStateSubStep(this->MemItem(gp_index).GetPlasticState_n());
                this->MemItem(gp_index).Setu_sub_step(this->MemItem(gp_index).Getu_n());
            }else{
                
                TPZFNMatrix<36,REAL> Dep(6,6,0.0);
                TPZTensor<STATE> epsilon,sigma;
                Epsilon(data,epsilon);
                T plastic_integrator(m_plastic_integrator);
                plastic_integrator.SetState(this->MemItem(gp_index).GetPlasticState());
                plastic_integrator.ApplyStrainComputeSigma(epsilon,sigma,&Dep);
                
//                /// Compute Kep with expression 27
//                TPZTensor<STATE> last_epsilon, last_sigma;
//                if (m_simulation_data->Get_must_use_sub_stepping_Q()) {
//                    last_epsilon = this->MemItem(gp_index).GetPlasticStateSubStep().m_eps_t;
//                }else{
//                    last_epsilon = this->MemItem(gp_index).GetPlasticState().m_eps_t;
//                }
//                last_sigma = this->MemItem(gp_index).GetSigma();
//
//                REAL delta_sigma = (sigma.I1()-last_sigma.I1()) / 3.0;
//                REAL delta_epsilon = (epsilon.I1()-last_epsilon.I1());
//                REAL Kep = delta_sigma / delta_epsilon;
                
                REAL K_ep_xx = (Dep(0,0) + Dep(3,0) + Dep(5,0))/3.0;
                REAL K_ep_yy = (Dep(0,3) + Dep(3,3) + Dep(5,3))/3.0;
                REAL K_ep_zz = (Dep(0,5) + Dep(3,5) + Dep(5,5))/3.0;
                REAL Kep = (K_ep_xx + K_ep_yy + K_ep_zz) / 3.0;
                REAL Ks = this->MemItem(gp_index).Ks();
                REAL alpha = 1.0 - (Kep/Ks);

                this->MemItem(gp_index).SetAlpha(alpha);
                this->MemItem(gp_index).SetPlasticStateSubStep(this->MemItem(gp_index).GetPlasticState_n());
                this->MemItem(gp_index).Setu_sub_step(this->MemItem(gp_index).Getu_n()) ;
                this->MemItem(gp_index).SetPlasticState(this->MemItem(gp_index).GetPlasticState_n());
                this->MemItem(gp_index).SetSigma(this->MemItem(gp_index).GetSigma_n());
                this->MemItem(gp_index).Setu(this->MemItem(gp_index).Getu_n());
            
            }
        
            return;
        }
        
//        TPZFNMatrix<36,REAL> Dep(6,6,0.0);
        TPZTensor<STATE> epsilon,sigma;
        Epsilon(data,epsilon);
        T plastic_integrator(m_plastic_integrator);
        plastic_integrator.SetState(this->MemItem(gp_index).GetPlasticState());
        plastic_integrator.ApplyStrainComputeSigma(epsilon,sigma);
        
//        REAL K_ep_xx = (Dep(0,0) + Dep(3,0) + Dep(5,0))/3.0;
//        REAL K_ep_yy = (Dep(0,3) + Dep(3,3) + Dep(5,3))/3.0;
//        REAL K_ep_zz = (Dep(0,5) + Dep(3,5) + Dep(5,5))/3.0;
//        REAL Kep = (K_ep_xx + K_ep_yy + K_ep_zz) / 3.0;
//        REAL Ks = this->MemItem(gp_index).Ks();
//        REAL alpha = 1.0 - (Kep/Ks);
        
        if (m_simulation_data->IsCurrentStateQ()) {
            
            this->MemItem(gp_index).SetPlasticState_n(plastic_integrator.fN);
            this->MemItem(gp_index).SetSigma_n(sigma);
            
            TPZManVector<STATE,3> delta_u    = data.sol[0];
            TPZManVector<STATE,3> u_n(m_dimension,0.0);
            if (m_simulation_data->Get_must_use_sub_stepping_Q()) {
                TPZManVector<STATE,3> u(this->MemItem(gp_index).Getu_sub_step());
                for (int i = 0; i < m_dimension; i++) {
                    u_n[i] = delta_u[i] + u[i];
                }
            }else{
                TPZManVector<STATE,3> u(this->MemItem(gp_index).Getu());
                for (int i = 0; i < m_dimension; i++) {
                    u_n[i] = delta_u[i] + u[i];
                }
            }
            this->MemItem(gp_index).Setu_n(u_n);
            
            /// Geomechanics porosity correction.
            {
                REAL alpha = this->MemItem(gp_index).Alpha();
                REAL Kdr   = this->MemItem(gp_index).Kdr();
                REAL p      = this->MemItem(gp_index).p();
                REAL p_n   = this->MemItem(gp_index).p_n();
                
                REAL epsilon_p_v   = (this->MemItem(gp_index).GetPlasticState().m_eps_p.I1()/3);
                REAL epsilon_p_v_n = (this->MemItem(gp_index).GetPlasticState_n().m_eps_p.I1()/3);
                REAL phi_p   = alpha * epsilon_p_v;
                REAL phi_p_n = alpha * epsilon_p_v_n;
                
                REAL sigma_t_v   = (this->MemItem(gp_index).GetSigma().I1()/3) - alpha * p;
                REAL sigma_t_v_n = (this->MemItem(gp_index).GetSigma_n().I1()/3)  - alpha * p_n;
                REAL geo_delta_phi_n = (alpha/Kdr)*(sigma_t_v_n-sigma_t_v) + (phi_p_n - phi_p); //  Geomechanic update.
                
                this->MemItem(gp_index).Setdelta_phi(geo_delta_phi_n);
            }
            
            { ///  Check for the need of substeps
                REAL norm = (this->MemItem(gp_index).GetPlasticState_n().m_eps_p - this->MemItem(gp_index).GetPlasticStateSubStep().m_eps_p).Norm();
                if (norm >= m_simulation_data->Get_max_plastic_strain()) {
                    m_simulation_data->Set_must_use_sub_stepping_Q(true);
                }
            }
            
        }else{
            this->MemItem(gp_index).SetPlasticState(plastic_integrator.fN);
            this->MemItem(gp_index).SetSigma(sigma);
            
            TPZManVector<STATE,3> u    = data.sol[0];
            this->MemItem(gp_index).Setu(u);
            
        }

    }
    
    TPZFMatrix<REAL> ek_fake;
    ek_fake.Resize(ef.Rows(),ef.Rows());
    this->Contribute(data, weight, ek_fake, ef);

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
