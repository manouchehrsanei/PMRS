//
//  TPMRSPoroElastoPlastic.h
//  PZ
//
//  Created by Omar and Manouchehr on 8/28/16.
//
//

#include "TPMRSPoroElastoPlastic.h"


template <class T, class TMEM>
TPMRSPoroElastoPlastic<T,TMEM>::TPMRSPoroElastoPlastic() : TPZMatWithMem<TMEM>() {
    
    m_simulation_data   = NULL;
    m_dimension         = 0;
    m_c                 = 0;
    m_eta               = 0;
    m_rho_0             = 0;
    m_scale_factor      = 1;
    m_theta_scheme      = 1;
}

template <class T, class TMEM>
TPMRSPoroElastoPlastic<T,TMEM>::TPMRSPoroElastoPlastic(int matid) : TPZMatWithMem<TMEM>(matid) {
    
    m_simulation_data   = NULL;
    m_dimension         = 0;
    m_c                 = 0;
    m_eta               = 0;
    m_rho_0             = 0;
    m_scale_factor      = 1;
    m_theta_scheme      = 1;
}

template <class T, class TMEM>
TPMRSPoroElastoPlastic<T,TMEM>::~TPMRSPoroElastoPlastic(){
    
}

template <class T, class TMEM>
TPMRSPoroElastoPlastic<T,TMEM>::TPMRSPoroElastoPlastic(const TPMRSPoroElastoPlastic& other): TPZMatWithMem<TMEM>(other){
    m_simulation_data    = other.m_simulation_data;
    m_dimension          = other.m_dimension;
    m_plastic_integrator = other.m_plastic_integrator;
    m_c                  = other.m_c;
    m_eta                = other.m_eta;
    m_rho_0              = other.m_rho_0;
    m_scale_factor       = other.m_scale_factor;
    m_phi_model          = other.m_phi_model;
    m_kappa_model        = other.m_kappa_model;
    m_theta_scheme       = other.m_theta_scheme;
}

template <class T, class TMEM>
TPMRSPoroElastoPlastic<T,TMEM> & TPMRSPoroElastoPlastic<T,TMEM>::operator = (const TPMRSPoroElastoPlastic& other){
    
    // check for self-assignment
    if(&other == this){
        return *this;
    }
    
    m_simulation_data    = other.m_simulation_data;
    m_dimension          = other.m_dimension;
    m_plastic_integrator = other.m_plastic_integrator;
    m_c                  = other.m_c;
    m_eta                = other.m_eta;
    m_rho_0              = other.m_rho_0;
    m_scale_factor       = other.m_scale_factor;
    m_phi_model          = other.m_phi_model;
    m_kappa_model        = other.m_kappa_model;
    m_theta_scheme       = other.m_theta_scheme;
    return *this;
    
}

template <class T, class TMEM>
void TPMRSPoroElastoPlastic<T,TMEM>::Print(std::ostream & out){
    
    out << " Material name : " << this->Name() << "\n";
    out << " Pointer to TPMRSSimulationData : " << m_simulation_data << "\n";
    out << " Material dimension : " << m_dimension << "\n";
    m_plastic_integrator.Print(out);
    out << " Fluid compressibility : " << m_c << "\n";
    out << " Fluid viscosity : " << m_eta << "\n";
    out << " Fluid density : " << m_rho_0 << "\n";
    out << " Scale factor  : " << m_scale_factor << "\n";
    out << " Theta scheme  : " << m_theta_scheme << "\n";
    m_phi_model.Print(out);
    m_kappa_model.Print(out);
    out << "\t Base class print:\n";
    TPZMaterial::Print(out);
}

template <class T, class TMEM>
int TPMRSPoroElastoPlastic<T,TMEM>::NStateVariables() const {
    return 1;
}

template <class T, class TMEM>
void TPMRSPoroElastoPlastic<T,TMEM>::FillDataRequirements(TPZVec<TPZMaterialData > &datavec){
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
    }
}

template <class T, class TMEM>
void TPMRSPoroElastoPlastic<T,TMEM>::FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec){
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
        datavec[idata].fNeedsNormal = true;
    }
}

template <class T, class TMEM>
void TPMRSPoroElastoPlastic<T,TMEM>::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
    // Getting weight functions
    TPZFMatrix<REAL>  & phi_u     =  datavec[m_u_b].phi;
    int n_phi_u = phi_u.Rows();
    int first_u = 0;
    
    TPZFNMatrix<40,REAL> grad_phi_u(3,n_phi_u);
    TPZAxesTools<REAL>::Axes2XYZ(datavec[m_u_b].dphix, grad_phi_u, datavec[m_u_b].axes);
    
    REAL dvdx,dvdy,dvdz;
    TPZTensor<STATE> epsilon,sigma;
    
    TPZFNMatrix<36,STATE> De(6,6,0.0);
    
    /// Get initial effective stress state
    int gp_index = datavec[m_u_b].intGlobPtIndex;
    TPZTensor<REAL> & sigma_0 = this->MemItem(gp_index).GetSigma_0();
    
    /// Get current effective stress state
    Epsilon(datavec[m_u_b],epsilon);
    Sigma(datavec[m_u_b], epsilon, sigma, &De);
    
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
    
    
    
    /// Getting weight functions
    TPZFMatrix<REAL>  & phi_p     =  datavec[m_p_b].phi;
    int n_phi_p = phi_p.Rows();
    
    TPZFNMatrix<40,REAL> grad_phi_p(m_dimension,n_phi_p);
    TPZFNMatrix<9,REAL> grad_p(m_dimension,1);
    TPZAxesTools<REAL>::Axes2XYZ(datavec[m_p_b].dphix, grad_phi_p, datavec[m_p_b].axes);
    TPZAxesTools<REAL>::Axes2XYZ(datavec[m_p_b].dsol[0], grad_p, datavec[m_p_b].axes);
    
    
    /// Get the pressure at the integrations points
    TMEM & memory  = this->MemItem(gp_index);
    
    /// Time
    STATE dt       = m_simulation_data->dt();
    STATE p_n      = datavec[m_p_b].sol[0][0];
    
    STATE p_0      = memory.p_0();
    STATE p        = memory.p();
    REAL alpha = this->MemItem(gp_index).Alpha();
    /// Biot term coupling
    {
    
        if (m_dimension == 2) {
            for(int iu = 0; iu < n_phi_u; iu++ )
            {
                dvdx = grad_phi_u(0,iu);
                dvdy = grad_phi_u(1,iu);
                
                ef(m_dimension*iu+0 + first_u)   +=    -1.0 * weight * (alpha*(p_n-p_0)*dvdx);    // x direction
                ef(m_dimension*iu+1 + first_u)   +=    -1.0 * weight * (alpha*(p_n-p_0)*dvdy);    // y direction
                
                for (int jp = 0; jp < n_phi_p; jp++)
                {
                    ek(m_dimension*iu+0 + first_u,jp+m_dimension*n_phi_u)   +=    -1.0 * weight * (alpha*(phi_p(jp,0))*dvdx);    // x direction
                    ek(m_dimension*iu+1 + first_u,jp+m_dimension*n_phi_u)   +=    -1.0 * weight * (alpha*(phi_p(jp,0))*dvdy);    // y direction
                }
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
                
                for (int jp = 0; jp < n_phi_p; jp++)
                {
                    ek(m_dimension*iu+0 + first_u,jp+m_dimension*n_phi_u)   +=    -1.0 * weight * (alpha*(phi_p(jp,0))*dvdx);    // x direction
                    ek(m_dimension*iu+1 + first_u,jp+m_dimension*n_phi_u)   +=    -1.0 * weight * (alpha*(phi_p(jp,0))*dvdy);    // y direction
                    ek(m_dimension*iu+2 + first_u,jp+m_dimension*n_phi_u)   +=    -1.0 * weight * (alpha*(phi_p(jp,0))*dvdz);    // z direction
                }
                
            }
        }
    }
    
    
    
    TPZManVector<STATE,3> & last_Kl_grad_p        = memory.f_vec();
    
    STATE phi_n,dphi_ndp,phi;
    REAL phi_0 = memory.phi_0();
    this->porosity(gp_index,phi_n,dphi_ndp,phi);
    
    REAL kappa_n;
    REAL dkappa_ndphi,dkappa_ndp;
    this->permeability(gp_index, kappa_n, dkappa_ndphi, phi_n, phi_0);
    dkappa_ndp = dkappa_ndphi * dphi_ndp;
    
    TPZFNMatrix<9,REAL> K(3,3),dKdp(3,3),dKdsigma(3,3);
    
    K.Zero();
    K(0,0) = memory.kappa_n();
    K(1,1) = memory.kappa_n();
    K(2,2) = memory.kappa_n();
    
    dKdp.Zero();
    dKdp(0,0) = dkappa_ndp;
    dKdp(1,1) = dkappa_ndp;
    dKdp(2,2) = dkappa_ndp;
    
    
    
    STATE rho        = m_rho_0 * (1 + (m_c/m_scale_factor)*(p-p_0));
    STATE rho_n      = m_rho_0 * (1 + (m_c/m_scale_factor)*(p_n-p_0));
    STATE drho_ndp_n = m_c/m_scale_factor;
    STATE lambda     = rho_n/m_eta;
    
    /// update for Kdr
    REAL dphi_nds_eff_vol = alpha;
    REAL dkappa_nds_eff_vol = dkappa_ndphi * dphi_nds_eff_vol;
    
    dKdsigma.Zero();
    dKdsigma(0,0) = dkappa_nds_eff_vol;
    dKdsigma(1,1) = dkappa_nds_eff_vol;
    dKdsigma(2,2) = dkappa_nds_eff_vol;
    
    /// Defining local variables
    TPZFNMatrix<3,STATE> Kl_grad_p_(3,1),dKdpl_grad_p_(3,1),dKdsigmal_grad_p_(3,1);
    for (int i = 0; i < Dimension(); i++)
    {
        STATE dot     = 0.0;
        STATE dKdpdot = 0.0;
        STATE dKdsigmadot = 0.0;
        for (int j = 0; j < Dimension(); j++)
        {
            dot            += (1.0/m_scale_factor)*K(i,j)*grad_p(j,0);
            dKdpdot        += dKdp(i,j)*grad_p(j,0);
            dKdsigmadot    += dKdsigma(i,j)*grad_p(j,0);
        }
        
        Kl_grad_p_(i,0)         = lambda * dot;
        dKdpl_grad_p_(i,0)      = lambda * dKdpdot;
        dKdsigmal_grad_p_(i,0)  = lambda * dKdsigmadot;
    }
    
    /// Integration point contribution
    TPZFNMatrix<3,STATE> phi_q_i(3,1), phi_q_j(3,1);
    
    if(!m_simulation_data->IsCurrentStateQ()){
        return;
    }
    
    TPZFNMatrix<3,STATE> Kl_grad_phi_j_(3,1);
    for (int ip = 0; ip < n_phi_p; ip++)
    {
        
        STATE Kl_grad_p_dot_grad_phi        = 0.0;
        STATE last_Kl_grad_p_dot_grad_phi   = 0.0;
        STATE dKdpl_grad_p_dot_grad_phi     = 0.0;
        STATE dKds_vol_l_grad_p_dot_grad_phi     = 0.0;
        for (int i = 0; i < Dimension(); i++)
        {
            Kl_grad_p_dot_grad_phi      += Kl_grad_p_(i,0)*grad_phi_p(i,ip);
            last_Kl_grad_p_dot_grad_phi += last_Kl_grad_p[i]*grad_phi_p(i,ip);
            dKdpl_grad_p_dot_grad_phi   += dKdpl_grad_p_(i,0)*grad_phi_p(i,ip);
            dKds_vol_l_grad_p_dot_grad_phi   += dKdsigmal_grad_p_(i,0)*grad_phi_p(i,ip);
            
        }
        
        
        REAL Kl_grad_p_dot_grad_phi_avg = m_theta_scheme*Kl_grad_p_dot_grad_phi + (1.0-m_theta_scheme)*last_Kl_grad_p_dot_grad_phi;
        
        ef(ip+m_dimension*n_phi_u) +=  weight * ( Kl_grad_p_dot_grad_phi_avg + (1.0/dt) * ( phi_n*rho_n - phi*rho ) * phi_p(ip,0) );
        
        if (m_dimension == 2) {
            for(int ju = 0; ju < n_phi_u; ju++ )
            {
                dvdx = grad_phi_u(0,ju);
                dvdy = grad_phi_u(1,ju);
                
                ek(ip+m_dimension*n_phi_u,m_dimension*ju+0 + first_u)   +=    weight * (m_theta_scheme*dKds_vol_l_grad_p_dot_grad_phi*dvdx+ (1.0/dt) * (dphi_nds_eff_vol*rho_n*dvdx) * phi_p(ip,0));    // x direction
                ek(ip+m_dimension*n_phi_u,m_dimension*ju+1 + first_u)   +=    weight * (m_theta_scheme*dKds_vol_l_grad_p_dot_grad_phi*dvdy+ (1.0/dt) * (dphi_nds_eff_vol*rho_n*dvdy) * phi_p(ip,0));    // y direction
                
            }
        }
        else{
            REAL dvdz;
            for(int ju = 0; ju < n_phi_u; ju++ )
            {
                dvdx = grad_phi_u(0,ju);
                dvdy = grad_phi_u(1,ju);
                dvdz = grad_phi_u(2,ju);
                
                ek(ip+m_dimension*n_phi_u,m_dimension*ju+0 + first_u)   +=    weight * (m_theta_scheme*dKds_vol_l_grad_p_dot_grad_phi*dvdx+ (1.0/dt) * (dphi_nds_eff_vol*rho_n*dvdx) * phi_p(ip,0));    // x direction
                ek(ip+m_dimension*n_phi_u,m_dimension*ju+1 + first_u)   +=    weight * (m_theta_scheme*dKds_vol_l_grad_p_dot_grad_phi*dvdy+ (1.0/dt) * (dphi_nds_eff_vol*rho_n*dvdy) * phi_p(ip,0));    // y direction
                ek(ip+m_dimension*n_phi_u,m_dimension*ju+2 + first_u)   +=    weight * (m_theta_scheme*dKds_vol_l_grad_p_dot_grad_phi*dvdz+ (1.0/dt) * (dphi_nds_eff_vol*rho_n*dvdz) * phi_p(ip,0));    // z direction
                
            }
        }
        
        for (int jp = 0; jp < n_phi_p; jp++)
        {
            
            for (int i = 0; i < Dimension(); i++)
            {
                STATE dot = 0.0;
                for (int j =0; j < Dimension(); j++)
                {
                    dot    += (1.0/m_scale_factor)*K(i,j)*grad_phi_p(j,jp);
                }
                
                Kl_grad_phi_j_(i,0)     = lambda * dot;
            }
            
            STATE Kl_grad_phi_j_dot_grad_phi_j = 0.0;
            for (int i = 0; i < Dimension(); i++)
            {
                Kl_grad_phi_j_dot_grad_phi_j    += Kl_grad_phi_j_(i,0)*grad_phi_p(i,ip);
                
            }
            
            ek(ip+m_dimension*n_phi_u,jp+m_dimension*n_phi_u) +=  weight * ( m_theta_scheme*Kl_grad_phi_j_dot_grad_phi_j + m_theta_scheme*dKdpl_grad_p_dot_grad_phi + (1.0/dt) * ( phi_n * drho_ndp_n + dphi_ndp * rho_n ) * phi_p(jp,0)  * phi_p(ip,0) );
        }
        
    }
    
    if (m_simulation_data->GetTransferCurrentToLastQ()) {
        for (int i = 0; i < Dimension(); i++)
        {
            last_Kl_grad_p[i] = Kl_grad_p_(i,0);
        }
        this->MemItem(gp_index).Setf_vec(last_Kl_grad_p); // Update process for f function.
        this->MemItem(gp_index).Setp(this->MemItem(gp_index).p_n());
        this->MemItem(gp_index).Setq(this->MemItem(gp_index).q_n());
    }
    
}

template <class T, class TMEM>
void TPMRSPoroElastoPlastic<T,TMEM>::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef){
    
    if (m_simulation_data->Get_must_accept_solution_Q()) {
        
        int gp_index = datavec[m_u_b].intGlobPtIndex;
        
        TPZTensor<STATE> epsilon,sigma;
        Epsilon(datavec[m_u_b],epsilon);
        TPZFNMatrix<36,STATE> Dep(6,6,0.0);
        T plastic_integrator(m_plastic_integrator);
        plastic_integrator.SetState(this->MemItem(gp_index).GetPlasticState());
        plastic_integrator.ApplyStrainComputeSigma(epsilon,sigma,&Dep);
        
        REAL K_ep_xx = (Dep(0,0) + Dep(3,0) + Dep(5,0))/3.0;
        REAL K_ep_yy = (Dep(0,3) + Dep(3,3) + Dep(5,3))/3.0;
        REAL K_ep_zz = (Dep(0,5) + Dep(3,5) + Dep(5,5))/3.0;
        REAL Kep = (K_ep_xx + K_ep_yy + K_ep_zz) / 3.0;
        REAL Ks = this->MemItem(gp_index).Ks();
        REAL alpha = 1.0 - (Kep/Ks);
        
        if (m_simulation_data->IsCurrentStateQ()) {
            
            this->MemItem(gp_index).SetPlasticState_n(plastic_integrator.fN);
            this->MemItem(gp_index).SetSigma_n(sigma);
            TPZManVector<STATE,3> delta_u    = datavec[m_u_b].sol[0];
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
            
            // Pressure variable
            STATE p = datavec[m_p_b].sol[0][0];
            
            // flux variable
            TPZManVector<REAL,3> q(3);
            {
                REAL p_0      = this->MemItem(gp_index).p_0();
                REAL p_n      = p;
                REAL rho_n    = m_rho_0 * (1 + m_c/m_scale_factor*(p_n-p_0)); //  Provide the compressibility in MPa
                REAL lambda   = rho_n/m_eta;
                REAL k        = this->MemItem(gp_index).kappa_n();
                TPZFNMatrix<9,REAL> grad_p(m_dimension,1);
                TPZAxesTools<REAL>::Axes2XYZ(datavec[m_p_b].dsol[0], grad_p, datavec[m_p_b].axes);
                for (int i = 0; i < m_dimension; i++)
                {
                    q[i] = (- (1.0/m_scale_factor) * k * lambda * grad_p[i]);
                }
                
                
                this->MemItem(gp_index).Setp_n(p);
                this->MemItem(gp_index).Setq_n(q);
            }
            
            { ///  Check for the need of substeps
                REAL norm = (this->MemItem(gp_index).GetPlasticState_n().m_eps_p - this->MemItem(gp_index).GetPlasticStateSubStep().m_eps_p).Norm();
                if (norm >= m_simulation_data->Get_max_plastic_strain()) {
                    m_simulation_data->Set_must_use_sub_stepping_Q(true);
                }
            }
            
            if (m_simulation_data->GetTransferCurrentToLastQ()) {
                
                if (m_simulation_data->Get_must_use_sub_stepping_Q()) {
                    this->MemItem(gp_index).SetPlasticStateSubStep(this->MemItem(gp_index).GetPlasticState_n());
                    this->MemItem(gp_index).Setu_sub_step(this->MemItem(gp_index).Getu_n());
                }else{
                    this->MemItem(gp_index).SetAlpha(alpha);
                    this->MemItem(gp_index).SetPlasticStateSubStep(this->MemItem(gp_index).GetPlasticState_n());
                    this->MemItem(gp_index).Setu_sub_step(this->MemItem(gp_index).Getu_n()) ;
                    this->MemItem(gp_index).SetPlasticState(this->MemItem(gp_index).GetPlasticState_n());
                    this->MemItem(gp_index).SetSigma(this->MemItem(gp_index).GetSigma_n());
                    this->MemItem(gp_index).Setu(this->MemItem(gp_index).Getu_n());
                }
                
            }
            
        }else{
            this->MemItem(gp_index).SetPlasticState(plastic_integrator.fN);
            this->MemItem(gp_index).SetSigma(sigma);
            
            TPZManVector<STATE,3> u    = datavec[m_u_b].sol[0];
            this->MemItem(gp_index).Setu(u);
            
        }
        
    }
    
    TPZFMatrix<STATE>  ek_fake(ef.Rows(),ef.Rows(),0.0);
    this->Contribute(datavec, weight, ek_fake, ef);
    return;
}

template <class T, class TMEM>
void TPMRSPoroElastoPlastic<T,TMEM>::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
    
    if(m_dimension == 3){
        this->ContributeBC_3D(datavec,weight,ek,ef,bc);
        return;
    }
    
    TPZBndCondWithMem<TMEM> & bc_with_memory = dynamic_cast<TPZBndCondWithMem<TMEM> &>(bc);
    int gp_index = datavec[m_u_b].intGlobPtIndex;
    TMEM & memory = bc_with_memory.MemItem(gp_index);
    
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
            TPZManVector<STATE,3> delta_u    = datavec[m_u_b].sol[0];
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
            TPZManVector<STATE,3> u    = datavec[m_u_b].sol[0];
            bc_with_memory.MemItem(gp_index).Setu(u);
        }
    }
    
    TPZFMatrix<REAL>  &phiu = datavec[m_u_b].phi;
    TPZFMatrix<REAL>  &phip = datavec[m_p_b].phi;
    TPZManVector<STATE,3> delta_u    = datavec[m_u_b].sol[0];
    REAL p  = datavec[m_p_b].sol[0][0];
    
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
    int phrp = phip.Rows();
    int in,jn;
    
    REAL BigNumber = TPZDiscontinuousGalerkin::gBigNumber;
    
    
    
    switch (bc.Type())
    {
    
        case 4 : /// DunNq
            /// Dirichlet of normal displacement and normal flux
            
        {
            TPZManVector<REAL,3> n = datavec[m_u_b].normal;
            REAL v[1];
            v[0] = bc.Val2()(0,0);    //    Un displacement
            
            for(in = 0 ; in < phru; in++)
            {
                ///    Contribution for load Vector
                ef(m_dimension*in+0,0)      += BigNumber*((u_n[0] - v[0])*n[0])*phiu(in,0)*weight;    // X displacement Value
                ef(m_dimension*in+1,0)      += BigNumber*((u_n[1] - v[0])*n[1])*phiu(in,0)*weight;    // Y displacement Value
                
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    ///    Contribution for Stiffness Matrix
                    ek(m_dimension*in+0,m_dimension*jn+0)    += BigNumber*n[0]*phiu(in,0)*phiu(jn,0)*weight;    // X displacement
                    ek(m_dimension*in+1,m_dimension*jn+1)    += BigNumber*n[1]*phiu(in,0)*phiu(jn,0)*weight;    // Y displacement
                    
                }
            }
            
            REAL Value = bc.Val2()(1,0);
            STATE last_qn_N = memory.f();
            STATE current_qn = Value;
            STATE qn_N = m_theta_scheme*current_qn+(1.0-m_theta_scheme)*last_qn_N;
            for (int ip = 0; ip < phrp; ip++)
            {
                ef(ip+m_dimension*phru) += -1.0 * weight * qn_N * phip(ip,0);
            }
            if (m_simulation_data->GetTransferCurrentToLastQ()) {
                memory.Setf(current_qn);
            }
            
            break;
        }
           
            
        case 5 : /// DuxDp
            /// Dirichlet in x direction of displacement and Dirichlet of pore pressure
            
        {
            REAL v[1];
            v[0] = bc.Val2()(0,0);    //    Ux displacement
            
            for(in = 0 ; in < phru; in++)
            {
                ///    Contribution for load Vector
                ef(m_dimension*in+0,0)      += BigNumber*(u_n[0] - v[0])*phiu(in,0)*weight;    // x displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    ///    Contribution for Stiffness Matrix
                    ek(m_dimension*in+0,m_dimension*jn+0)    += BigNumber*phiu(in,0)*phiu(jn,0)*weight;    // x displacement
                    
                }
            }
            
            REAL Value = bc.Val2()(1,0);
            REAL p_D = Value;
            for (int ip = 0; ip < phrp; ip++)
            {
                ef(ip+m_dimension*phru) += weight * m_scale_factor * BigNumber * (p - p_D) * phip(ip,0);
                
                for (int jp = 0; jp < phrp; jp++)
                {
                    
                    ek(ip+m_dimension*phru,jp+m_dimension*phru) += weight * m_scale_factor *  BigNumber * phip(jp,0) * phip(ip,0);
                }
            }

            
            break;
        }
            
            
            
        case 6 : /// DuxNq
            /// Dirichlet in x direction of displacement and normal flux
            
        {
            REAL v[1];
            v[0] = bc.Val2()(0,0);    //    Ux displacement
            
            for(in = 0 ; in < phru; in++)
            {
                ///    Contribution for load Vector
                ef(m_dimension*in+0,0)      += BigNumber*(u_n[0] - v[0])*phiu(in,0)*weight;    // x displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    ///    Contribution for Stiffness Matrix
                    ek(m_dimension*in+0,m_dimension*jn+0)    += BigNumber*phiu(in,0)*phiu(jn,0)*weight;    // x displacement
                    
                }
            }
            
            REAL Value = bc.Val2()(1,0);
            STATE last_qn_N = memory.f();
            STATE current_qn = Value;
            STATE qn_N = m_theta_scheme*current_qn+(1.0-m_theta_scheme)*last_qn_N;
            for (int ip = 0; ip < phrp; ip++)
            {
                ef(ip+m_dimension*phru) += -1.0 * weight * qn_N * phip(ip,0);
            }
            if (m_simulation_data->GetTransferCurrentToLastQ()) {
                memory.Setf(current_qn);
            }
            
            break;
        }
            
            
        case 7 : /// DuyDp
            /// Dirichlet in y direction of displacement and Dirichlet of pore pressure
            
        {
            REAL v[1];
            v[0] = bc.Val2()(0,0);    //    Uy displacement
            
            for(in = 0 ; in < phru; in++)
            {
                ///    Contribution for load Vector
                ef(m_dimension*in+1,0)      += BigNumber*(u_n[1] - v[0])*phiu(in,0)*weight;    // y displacement Value
                
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    ///    Contribution for Stiffness Matrix
                    ek(m_dimension*in+1,m_dimension*jn+1)    += BigNumber*phiu(in,0)*phiu(jn,0)*weight;    // y displacement
                    
                }
            }
            
            REAL Value = bc.Val2()(1,0);
            REAL p_D = Value;
            for (int ip = 0; ip < phrp; ip++)
            {
                ef(ip+m_dimension*phru) += weight * m_scale_factor * BigNumber * (p - p_D) * phip(ip,0);
                
                for (int jp = 0; jp < phrp; jp++)
                {
                    
                    ek(ip+m_dimension*phru,jp+m_dimension*phru) += weight * m_scale_factor *  BigNumber * phip(jp,0) * phip(ip,0);
                }
            }

            
            break;
        }


        case 8 : /// DuyNq
            /// Dirichlet in y direction of displacement and normal flux
            
        {
            REAL v[1];
            v[0] = bc.Val2()(0,0);    //    Uy displacement
            
            for(in = 0 ; in < phru; in++)
            {
                ///    Contribution for load Vector
                ef(m_dimension*in+1,0)      += BigNumber*(u_n[1] - v[0])*phiu(in,0)*weight;    // y displacement Value
                
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    ///    Contribution for Stiffness Matrix
                    ek(m_dimension*in+1,m_dimension*jn+1)    += BigNumber*phiu(in,0)*phiu(jn,0)*weight;    // y displacement
                    
                }
            }
            
            REAL Value = bc.Val2()(1,0);
            STATE last_qn_N = memory.f();
            STATE current_qn = Value;
            STATE qn_N = m_theta_scheme*current_qn+(1.0-m_theta_scheme)*last_qn_N;
            for (int ip = 0; ip < phrp; ip++)
            {
                ef(ip+m_dimension*phru) += -1.0 * weight * qn_N * phip(ip,0);
            }
            if (m_simulation_data->GetTransferCurrentToLastQ()) {
                memory.Setf(current_qn);
            }
            
            break;
        }

            
        case 9 : /// NSDp
            /// Neumann of normal sigma and Dirichlet BC of PD
            
        {
                /// Neumann of normal sigma
                TPZManVector<REAL,3> n = datavec[m_u_b].normal;
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
            
            /// Dirichlet BC of PD
            REAL Value = bc.Val2()(6,0);
            REAL p_D = Value;
            for (int ip = 0; ip < phrp; ip++)
            {
                ef(ip+m_dimension*phru) += weight * m_scale_factor * BigNumber * (p - p_D) * phip(ip,0);
                
                for (int jp = 0; jp < phrp; jp++)
                {
                    
                    ek(ip+m_dimension*phru,jp+m_dimension*phru) += weight * m_scale_factor *  BigNumber * phip(jp,0) * phip(ip,0);
                }
            }
            
            
            break;
        }
            
            
            
        case 13 : /// NtnDp
            /// Neumann of traction and Dirichlet BC of PD
            
        {
            REAL v[1];
            v[0] = bc.Val2()(0,0);    //    Tn normal traction
            
            REAL tn = v[0];
            TPZManVector<REAL,3> n = datavec[m_u_b].normal;
            ///    Neumann condition for each state variable
            ///    Elasticity Equation
            for(in = 0 ; in <phru; in++)
            {
                ///   Normal Tension Components on neumman boundary
                ef(m_dimension*in+0,0)    += -1.0 * weight * tn * n[0] * phiu(in,0);        //    Tnx
                ef(m_dimension*in+1,0)    += -1.0 * weight * tn * n[1] * phiu(in,0);        //    Tny
            }
            
            /// Dirichlet BC of PD
            REAL Value = bc.Val2()(1,0);
            REAL p_D = Value;
            for (int ip = 0; ip < phrp; ip++)
            {
                ef(ip+m_dimension*phru) += weight * m_scale_factor * BigNumber * (p - p_D) * phip(ip,0);
                
                for (int jp = 0; jp < phrp; jp++)
                {
                    
                    ek(ip+m_dimension*phru,jp+m_dimension*phru) += weight * m_scale_factor *  BigNumber * phip(jp,0) * phip(ip,0);
                }
            }

            
            break;
        }
            
        case 14 : /// NtnNq
            /// Neumann of traction and normal flux
            
        {
            REAL v[1];
            v[0] = bc.Val2()(0,0);    //    Tn normal traction
            
            REAL tn = v[0];
            TPZManVector<REAL,3> n = datavec[m_u_b].normal;
            ///    Neumann condition for each state variable
            ///    Elasticity Equation
            for(in = 0 ; in <phru; in++)
            {
                ///   Normal Tension Components on neumman boundary
                ef(m_dimension*in+0,0)    += -1.0 * weight * tn * n[0] * phiu(in,0);        //    Tnx
                ef(m_dimension*in+1,0)    += -1.0 * weight * tn * n[1] * phiu(in,0);        //    Tny
            }
            
            
            REAL Value = bc.Val2()(1,0);
            STATE last_qn_N = memory.f();
            STATE current_qn = Value;
            STATE qn_N = m_theta_scheme*current_qn+(1.0-m_theta_scheme)*last_qn_N;
            for (int ip = 0; ip < phrp; ip++)
            {
                ef(ip+m_dimension*phru) += -1.0 * weight * qn_N * phip(ip,0);
            }
            if (m_simulation_data->GetTransferCurrentToLastQ()) {
                memory.Setf(current_qn);
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
void TPMRSPoroElastoPlastic<T,TMEM>::ContributeBC_3D(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
    
    TPZBndCondWithMem<TMEM> & bc_with_memory = dynamic_cast<TPZBndCondWithMem<TMEM> &>(bc);
    int gp_index = datavec[m_u_b].intGlobPtIndex;
    TMEM & memory = bc_with_memory.MemItem(gp_index);
    
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
            TPZManVector<STATE,3> delta_u    = datavec[m_u_b].sol[0];
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
            TPZManVector<STATE,3> u    = datavec[m_u_b].sol[0];
            bc_with_memory.MemItem(gp_index).Setu(u);
        }
    }
    
    TPZFMatrix<REAL>  &phiu = datavec[m_u_b].phi;
    TPZFMatrix<REAL>  &phip = datavec[m_p_b].phi;
    TPZManVector<STATE,3> delta_u    = datavec[m_u_b].sol[0];
    REAL p  = datavec[m_p_b].sol[0][0];
    
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
    int phrp = phip.Rows();
    int in,jn;
    
    REAL BigNumber = TPZDiscontinuousGalerkin::gBigNumber;
    
    
    
    switch (bc.Type())
    {
            
        case 4 : /// DunNq
            /// Dirichlet of normal displacement and normal flux
            
        {
            TPZManVector<REAL,3> n = datavec[m_u_b].normal;
            REAL v[1];
            v[0] = bc.Val2()(0,0);    //    Un displacement
            
            for(in = 0 ; in < phru; in++)
            {
                ///    Contribution for load Vector
                ef(m_dimension*in+0,0)      += BigNumber*((u_n[0] - v[0])*n[0])*phiu(in,0)*weight;    // X displacement Value
                ef(m_dimension*in+1,0)      += BigNumber*((u_n[1] - v[0])*n[1])*phiu(in,0)*weight;    // Y displacement Value
                ef(m_dimension*in+2,0)      += BigNumber*((u_n[2] - v[0])*n[2])*phiu(in,0)*weight;    // Z displacement Value
                
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    ///    Contribution for Stiffness Matrix
                    ek(m_dimension*in+0,m_dimension*jn+0)    += BigNumber*n[0]*phiu(in,0)*phiu(jn,0)*weight;    // X displacement
                    ek(m_dimension*in+1,m_dimension*jn+1)    += BigNumber*n[1]*phiu(in,0)*phiu(jn,0)*weight;    // Y displacement
                    ek(m_dimension*in+2,m_dimension*jn+2)    += BigNumber*n[2]*phiu(in,0)*phiu(jn,0)*weight;    // Z displacement
                    
                }
            }
            
            REAL Value = bc.Val2()(1,0);
            STATE last_qn_N = memory.f();
            STATE current_qn = Value;
            STATE qn_N = m_theta_scheme*current_qn+(1.0-m_theta_scheme)*last_qn_N;
            for (int ip = 0; ip < phrp; ip++)
            {
                ef(ip+m_dimension*phru) += -1.0 * weight * qn_N * phip(ip,0);
            }
            if (m_simulation_data->GetTransferCurrentToLastQ()) {
                memory.Setf(current_qn);
            }
            
            break;
        }
            
        case 21 : /// NtnDp
            /// Neumann of traction and Dirichlet BC  PD
            
        {
            REAL v[1];
            v[0] = bc.Val2()(0,0);    //    Tn normal traction
            
            REAL tn = v[0];
            TPZManVector<REAL,3> n = datavec[m_u_b].normal;
            ///    Neumann condition for each state variable
            ///    Elasticity Equation
            for(in = 0 ; in <phru; in++)
            {
                ///   Normal Tension Components on neumman boundary
                ef(m_dimension*in+0,0)    += -1.0 * weight * tn * n[0] * phiu(in,0);        //    Tnx
                ef(m_dimension*in+1,0)    += -1.0 * weight * tn * n[1] * phiu(in,0);        //    Tny
                ef(m_dimension*in+2,0)    += -1.0 * weight * tn * n[2] * phiu(in,0);        //    Tnz
            }
            
            REAL Value = bc.Val2()(1,0);
            REAL p_D = Value;
            for (int ip = 0; ip < phrp; ip++)
            {
                ef(ip+m_dimension*phru) += weight * m_scale_factor * BigNumber * (p - p_D) * phip(ip,0);
                
                for (int jp = 0; jp < phrp; jp++)
                {
                    
                    ek(ip+m_dimension*phru,jp+m_dimension*phru) += weight * m_scale_factor *  BigNumber * phip(jp,0) * phip(ip,0);
                }
            }
            
            
            break;
        }
            
        case 22 : /// NtnNq
            /// Neumann of traction and normal flux
            
        {
            REAL v[1];
            v[0] = bc.Val2()(0,0);    //    Tn normal traction
            
            REAL tn = v[0];
            TPZManVector<REAL,3> n = datavec[m_u_b].normal;
            ///    Neumann condition for each state variable
            ///    Elasticity Equation
            for(in = 0 ; in <phru; in++)
            {
                ///   Normal Tension Components on neumman boundary
                ef(m_dimension*in+0,0)    += -1.0 * weight * tn * n[0] * phiu(in,0);        //    Tnx
                ef(m_dimension*in+1,0)    += -1.0 * weight * tn * n[1] * phiu(in,0);        //    Tny
                ef(m_dimension*in+2,0)    += -1.0 * weight * tn * n[2] * phiu(in,0);        //    Tnz
            }
            
            
            REAL Value = bc.Val2()(1,0);
            STATE last_qn_N = memory.f();
            STATE current_qn = Value;
            STATE qn_N = m_theta_scheme*current_qn+(1.0-m_theta_scheme)*last_qn_N;
            for (int ip = 0; ip < phrp; ip++)
            {
                ef(ip+m_dimension*phru) += -1.0 * weight * qn_N * phip(ip,0);
            }
            if (m_simulation_data->GetTransferCurrentToLastQ()) {
                memory.Setf(current_qn);
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
int TPMRSPoroElastoPlastic<T,TMEM>::VariableIndex(const std::string &name){
    if (!strcmp("u"     , name.c_str())) return 0;
    if (!strcmp("s"     , name.c_str())) return 1;
    if (!strcmp("e"     , name.c_str())) return 2;
    if (!strcmp("ep"    , name.c_str())) return 3;
    if (!strcmp("s_t"   , name.c_str())) return 4;
    if (!strcmp("p"     , name.c_str())) return 5;
    if (!strcmp("phi"   , name.c_str())) return 6;
    if (!strcmp("kappa" , name.c_str())) return 7;
    if (!strcmp("q"     , name.c_str())) return 8;
    if (!strcmp("div_q" , name.c_str())) return 9;
    return TPZMatWithMem<TMEM>::VariableIndex(name);
}

template <class T, class TMEM>
int TPMRSPoroElastoPlastic<T,TMEM>::NSolutionVariables(int var){
    switch(var) {
        case 0:
            return m_dimension; /// Vector
        case 1:
            return 9; /// Tensor
        case 2:
            return 9; /// Tensor
        case 3:
            return 9; /// Tensor
        case 4:
            return 9; /// Tensor
        case 5:
            return 1; /// Scalar
        case 6:
            return 1; /// Scalar
        case 7:
            return 1; /// Scalar
        case 8:
            return m_dimension; /// Vector
        case 9:
            return 1; /// Scalar
            
    }
    return TPZMatWithMem<TMEM>::NSolutionVariables(var);
}

template <class T, class TMEM>
void TPMRSPoroElastoPlastic<T,TMEM>::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout){
    
    long gp_index = data.intGlobPtIndex;
    TMEM & memory = this->GetMemory().get()->operator[](gp_index);
    Solout.Resize( this->NSolutionVariables(var));
    TPZTensor<REAL> epsilon_t = memory.GetPlasticState_n().m_eps_t;
    TPZTensor<REAL> epsilon_p = memory.GetPlasticState_n().m_eps_p;
    
    switch (var) {
        case 0:
        {
            for(int i = 0; i < m_dimension; i++){
                Solout[i] = memory.Getu_n()[i];
            }
        }
            break;
        case 1:
        {
            Solout[0] = memory.GetSigma_n().XX();
            Solout[1] = Solout[3] = memory.GetSigma_n().XY();
            Solout[2] = Solout[6] = memory.GetSigma_n().XZ();
            Solout[4] = memory.GetSigma_n().YY();
            Solout[5] = Solout[7] = memory.GetSigma_n().YZ();
            Solout[8] = memory.GetSigma_n().ZZ();
        }
            break;
        case 2:
        {
            Solout[0] = epsilon_t.XX();
            Solout[1] = Solout[3] = epsilon_t.XY();
            Solout[2] = Solout[6] = epsilon_t.XZ();
            Solout[4] = epsilon_t.YY();
            Solout[5] = Solout[7] = epsilon_t.YZ();
            Solout[8] = epsilon_t.ZZ();
        }
            break;
        case 3:
        {
            Solout[0] = epsilon_p.XX();
            Solout[1] = Solout[3] = epsilon_p.XY();
            Solout[2] = Solout[6] = epsilon_p.XZ();
            Solout[4] = epsilon_p.YY();
            Solout[5] = Solout[7] = epsilon_p.YZ();
            Solout[8] = epsilon_p.ZZ();
        }
            break;
        case 4:
        {
            Solout[0] = memory.GetSigma_n().XX() - memory.Alpha()*memory.p_n();
            Solout[1] = Solout[3] = memory.GetSigma_n().XY();
            Solout[2] = Solout[6] = memory.GetSigma_n().XZ();
            Solout[4] = memory.GetSigma_n().YY() - memory.Alpha()*memory.p_n();
            Solout[5] = Solout[7] = memory.GetSigma_n().YZ();
            Solout[8] = memory.GetSigma_n().ZZ() - memory.Alpha()*memory.p_n();
        }
            break;
        case 5:
        {
            Solout[0] = memory.p_n();
        }
            break;
        case 6:
        {
            Solout[0] = memory.phi_n();
        }
            break;
        case 7:
        {
            Solout[0] = memory.kappa_n();
        }
            break;
        case 8:
        {
            for (int i = 0; i < m_dimension; i++)
            {
                Solout[i] = memory.q_n()[i];
            }
        }
            break;
        case 9:
        {
                Solout[0] = -1942; // Meaningless variable.
        }
            break;
        default:
        {
            std::cout << "TPMRSPoroElastoPlastic<T,TMEM>:: Variable not implemented." << std::endl;
            DebugStop();
        }
            break;
    }
    
}

template <class T, class TMEM>
void TPMRSPoroElastoPlastic<T,TMEM>::Epsilon(TPZMaterialData &data, TPZTensor<REAL> & epsilon_t){
    
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

template <class T, class TMEM>
void TPMRSPoroElastoPlastic<T,TMEM>::Sigma(TPZMaterialData &data, TPZTensor<REAL> & epsilon_t, TPZTensor<REAL> & sigma, TPZFMatrix<REAL> * Dep){
    
    int gp_index = data.intGlobPtIndex;
    T plastic_integrator(m_plastic_integrator);
    bool SecantQ = true;
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
    
}

template <class T, class TMEM>
void TPMRSPoroElastoPlastic<T,TMEM>::porosity(long gp_index, REAL &phi_n, REAL &dphi_ndp, REAL &phi){
    
    TMEM & memory = this->MemItem(gp_index);
    
    REAL alpha = memory.Alpha();
    REAL Kdr   = memory.Kdr();
    REAL phi_0 = memory.phi_0();
    
    REAL p_0   = memory.p_0();
    REAL p     = memory.p();
    REAL p_n   = memory.p_n();
    
    REAL S = (1.0-alpha)*(alpha-phi_0)/Kdr;
    
    REAL sigma_v_0 = (memory.GetSigma_0().I1()/3);
    REAL sigma_v   = (memory.GetSigma().I1()/3);
    REAL sigma_v_n = (memory.GetSigma_n().I1()/3);
    
    REAL phi_p_0   = alpha * memory.GetPlasticState_0().m_eps_p.I1();
    REAL phi_p     = alpha * memory.GetPlasticState().m_eps_p.I1();
    REAL phi_p_n   = alpha * memory.GetPlasticState_n().m_eps_p.I1();
    
    phi = phi_0 + (alpha/Kdr)* (sigma_v-sigma_v_0) + (phi_p - phi_p_0) + S * (p - p_0);
    phi_n = phi_0 + (alpha/Kdr)* (sigma_v_n-sigma_v_0) + (phi_p_n - phi_p_0) + S * (p_n - p_0);
    dphi_ndp = S;
    
    this->MemItem(gp_index).Setphi_n(phi_n);
}


template <class T, class TMEM>
void TPMRSPoroElastoPlastic<T,TMEM>::permeability(long gp_index, REAL &kappa_n, REAL &dkappa_ndphi,REAL &phi,REAL &phi_0){
    
    TMEM & memory = this->MemItem(gp_index);
    REAL kappa_0  = memory.kappa_0();
    m_kappa_model.Permeability(kappa_n, dkappa_ndphi, kappa_0, phi, phi_0);
    this->MemItem(gp_index).Setkappa_n(kappa_n);
}
