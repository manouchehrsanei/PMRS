//
//  TPMRSRKSolver_impl.hpp
//  PMRS
//
//  Created by Manouchehr Sanei and Omar on 3/9/19.
//

#include "TPMRSRKSolver.h"


template <class T, class TMEM>
TPMRSRKSolver<T,TMEM>::TPMRSRKSolver(){
    
    m_n_steps   = 0;
    m_y_0.resize(0);
    m_re        = 0.0;
    m_rw        = 0.0;
    m_eta       = 0.0;
    m_cf        = 0.0;
    m_K_s       = 0.0;
    m_dr        = 0.0;
    m_is_RK4_Q  = false;
    m_is_Euler_Q = false;
    m_n_state   = 4;
    m_memory.resize(0);
    m_r_y.Resize(0, 0);
    m_lambda.resize(0);
    m_mu.resize(0);
    m_accept_solution_Q = false;
    m_sigma.Zero();
    m_eps_e.Zero();
    m_eps_p.Zero();
}

template <class T, class TMEM>
TPMRSRKSolver<T,TMEM>::~TPMRSRKSolver(){
    
}

template <class T, class TMEM>
TPMRSRKSolver<T,TMEM>::TPMRSRKSolver(const TPMRSRKSolver & other){
    
    m_n_steps            = other.m_n_steps;
    m_y_0                = other.m_y_0;
    m_re                 = other.m_re;
    m_rw                 = other.m_rw;
    m_eta                = other.m_eta;
    m_cf                 = other.m_cf;
    m_K_s                = other.m_K_s;
    m_dr                 = other.m_dr;
    m_default_memory     = other.m_default_memory;
    m_is_RK4_Q           = other.m_is_RK4_Q;
    m_n_state            = other.m_n_state;
    m_plastic_integrator = other.m_plastic_integrator;
    m_memory             = other.m_memory;
    m_r_y                = other.m_r_y;
    m_lambda             = other.m_lambda;
    m_mu                 = other.m_mu;
    m_accept_solution_Q  = other.m_accept_solution_Q;
    
}

template <class T, class TMEM>
TPMRSRKSolver<T,TMEM> & TPMRSRKSolver<T,TMEM>::operator=(const TPMRSRKSolver & other){
    
    /// check for self-assignment
    if(&other == this){
        return *this;
    }
    
    m_n_steps            = other.m_n_steps;
    m_y_0                = other.m_y_0;
    m_re                 = other.m_re;
    m_rw                 = other.m_rw;
    m_eta                = other.m_eta;
    m_cf                 = other.m_cf;
    m_K_s                = other.m_K_s;
    m_dr                 = other.m_dr;
    m_default_memory     = other.m_default_memory;
    m_is_RK4_Q           = other.m_is_RK4_Q;
    m_n_state            = other.m_n_state;
    m_plastic_integrator = other.m_plastic_integrator;
    m_memory             = other.m_memory;
    m_r_y                = other.m_r_y;
    m_lambda             = other.m_lambda;
    m_mu                 = other.m_mu;
    m_accept_solution_Q  = other.m_accept_solution_Q;
    
    return *this;
}

template <class T, class TMEM>
const std::string TPMRSRKSolver<T,TMEM>::Name() const{
    return "TPMRSRKSolver";
}

template <class T, class TMEM>
void TPMRSRKSolver<T,TMEM>::Print(std::ostream &out) const{
    
    out << Name() << std::endl;
//    out << "m_n_steps = " << m_n_steps << std::endl;
//    out << "m_y_0 = " << m_y_0 << std::endl;
//    out << "m_re = " << m_re << std::endl;
//    out << "m_rw = " << m_rw << std::endl;
//    out << "m_eta = " << m_eta << std::endl;
//    out << "m_cf = " << m_cf << std::endl;
//    out << "m_K_s = " << m_K_s << std::endl;
//    out << "m_dr = " << m_dr << std::endl;
//    out << "m_default_memory = " << m_default_memory << std::endl;
//    out << "m_is_RK4_Q = " << m_is_RK4_Q << std::endl;
//    out << "m_n_state = " << m_n_state << std::endl;
//    out << "m_plastic_integrator = " << m_plastic_integrator << std::endl;
//    out << "m_memory = " << m_memory << std::endl;
//    out << "m_r_y = " << m_r_y << std::endl;
//    out << "m_lambda = " << m_lambda << std::endl;
//    out << "m_mu = " << m_mu << std::endl;
    
}

template <class T, class TMEM>
void TPMRSRKSolver<T,TMEM>::Synchronize(){
    int n_points = m_n_steps + 1;
    m_memory.resize(n_points,m_default_memory);
    m_r_y.Resize(n_points, m_n_state+1);
    m_lambda.resize(n_points);
    m_mu.resize(n_points);
    for (int i = 0; i < n_points; i++) {
        m_lambda[i] = m_plastic_integrator.GetElasticResponse().Lambda();
        m_mu[i] = m_plastic_integrator.GetElasticResponse().G();
    }
}

template <class T, class TMEM>
std::vector<REAL> TPMRSRKSolver<T,TMEM>::f(int i, REAL & r, std::vector<REAL> & y){
    
    REAL ur = y[0];
    REAL sr = y[1];
    REAL qr = y[3];
    
    REAL phi   = Porosity(i,r,y);
    REAL kappa = Permeability(i,phi);
    
    std::vector<REAL> f(4);
    
    REAL l = this->lambda(i);
    REAL mu = this->mu(i);
    REAL alpha = m_memory[i].Alpha();
    REAL sigmar_0 = m_memory[i].GetSigma_0().XX();
    
    REAL eps_t_rr = (-r*sigmar_0+r*sr-l*ur)/(r*(l+2.0*mu));
    f[0] = eps_t_rr;
    f[1] = ((2*mu*(ur/r)) - (2.0*mu*(-r*sigmar_0+r*sr-l*ur)/(r*(l+2.0*mu))))/(r)-alpha*(m_eta/(m_rho*kappa))*qr;
    f[2] = -(m_eta/(m_rho*kappa))*qr;
    f[3] = -qr/r;

    return f;
}

template <class T, class TMEM>
std::vector<REAL> TPMRSRKSolver<T,TMEM>::EulerApproximation(int i, REAL & r, std::vector<REAL> & y){
    
    REAL h = m_dr;
    
    /// k1
    std::vector<REAL> k1;
    k1 = f(i,r,y);
    k1 = a_times_v(h,k1);
    
    std::vector<REAL> y_p_1(4,0.0);
    y_p_1 = a_add_b(y_p_1,y);
    y_p_1 = a_add_b(y_p_1,k1);
    
    return y_p_1;
}

template <class T, class TMEM>
std::vector<REAL> TPMRSRKSolver<T,TMEM>::RK2Approximation(int i, REAL & r, std::vector<REAL> & y){
    
    REAL h = m_dr;
    REAL s = 0.5;
    REAL hhalf = s*h;
    
    REAL half_r = hhalf + r;
    
    /// k1
    std::vector<REAL> k1;
    k1 = f(i,r,y);
    k1 = a_times_v(h,k1);
    
    /// k2
    std::vector<REAL> k2,k1_star;
    k1_star = a_times_v(s,k1);
    k1_star = a_add_b(y,k1_star);
    k2 = f(i,half_r,k1_star);
    k2 = a_times_v(h,k2);
    
    /// construct the approximation y_n+1
    std::vector<REAL> y_p_1(4),s_k2;
    y_p_1 = a_add_b(y_p_1,y);
    y_p_1 = a_add_b(y_p_1,k2);
    
    return y_p_1;
}

template <class T, class TMEM>
std::vector<REAL> TPMRSRKSolver<T,TMEM>::RK4Approximation(int i, REAL & r, std::vector<REAL> & y){
    
    REAL h = m_dr;
    REAL s = 0.5;
    REAL hhalf = s*h;
    
    REAL half_r = hhalf + r;
    REAL next_r = h + r;
    
    /// k1
    std::vector<REAL> k1;
    k1 = f(i,r,y);
    k1 = a_times_v(h,k1);
    
    /// k2
    std::vector<REAL> k2,k1_star;
    k1_star = a_times_v(s,k1);
    k1_star = a_add_b(y,k1_star);
    k2 = f(i,half_r,k1_star);
    k2 = a_times_v(h,k2);
    
    /// k3
    std::vector<REAL> k3,k2_star;
    k2_star = a_times_v(s,k2);
    k2_star = a_add_b(y,k2_star);
    k3 = f(i,half_r,k2_star);
    k3 = a_times_v(h,k3);
    
    /// k4
    std::vector<REAL> k4,k3_star;
    k3_star = a_add_b(y,k3);
    k4 = f(i,next_r,k3_star);
    k4 = a_times_v(h,k4);
    
    /// construct the approximation y_n+1
    std::vector<REAL> y_p_1(4),s_k1,s_k2,s_k3,s_k4;
    s_k1 = a_times_v(1.0/6.0,k1);
    s_k2 = a_times_v(1.0/3.0,k2);
    s_k3 = a_times_v(1.0/3.0,k3);
    s_k4 = a_times_v(1.0/6.0,k4);
    y_p_1 = a_add_b(y_p_1,y);
    y_p_1 = a_add_b(y_p_1,s_k1);
    y_p_1 = a_add_b(y_p_1,s_k2);
    y_p_1 = a_add_b(y_p_1,s_k3);
    y_p_1 = a_add_b(y_p_1,s_k4);
    return y_p_1;
}

template <class T, class TMEM>
std::vector<REAL> TPMRSRKSolver<T,TMEM>::a_times_v(const REAL a, std::vector<REAL> & v){
    std::vector<REAL> y(v.size());
    for (int i = 0; i < v.size(); i++) {
        y[i] = a*v[i];
    }
    return y;
}

template <class T, class TMEM>
std::vector<REAL> TPMRSRKSolver<T,TMEM>::a_add_b(std::vector<REAL> & a, std::vector<REAL> & b){
    if (a.size()!=b.size()) {
        DebugStop();
    }
    std::vector<REAL> y(a.size());
    for (int i = 0; i < a.size(); i++) {
        y[i] = a[i]+b[i];
    }
    return y;
}

template <class T, class TMEM>
void TPMRSRKSolver<T,TMEM>::ExecuteRKApproximation(){
    
    int n_points = m_n_steps + 1;
    
    std::vector<REAL> y = m_y_0;
    if(m_is_Re_Q){
        ReconstructAndAcceptPoint(0,m_re,y,true);
    }else{
        ReconstructAndAcceptPoint(0,m_rw,y,true);
    }
    
    AppendTo(0,y);
    
    for (int i = 1; i < n_points; i++) {
        
        REAL r_ref;
        if(m_is_Re_Q){
            r_ref =  m_re;
        }else{
            r_ref =  m_rw;
        }
        
        REAL r = m_dr*(i-1) + r_ref;
        
        if(m_is_RK4_Q){
            y = RK4Approximation(i-1,r,y);
        }else{
            if(m_is_Euler_Q){
                y = EulerApproximation(i-1,r,y);
            }else{
                y = RK2Approximation(i-1,r,y);
            }
        }
        ReconstructAndAcceptPoint(i,r,y);
        AppendTo(i,y);
    }
    
}

template <class T, class TMEM>
void TPMRSRKSolver<T,TMEM>::AppendTo(int i, std::vector<REAL> y){
    REAL r;
    if(m_is_Re_Q){
        r = m_dr*(i) + m_re;
    }else{
        r = m_dr*(i) + m_rw;
    }
    for (int k = 0; k < m_n_state; k++) {
        m_r_y(i,0) = r;
        m_r_y(i,k+1) = y[k];
    }
}

template <class T, class TMEM>
void TPMRSRKSolver<T,TMEM>::PrintRKApproximation(std::ostream &out){
    m_r_y.Print("rkdata = ",out,EMathematicaInput);
}

/// Print the secondary variables (s_r,s_t,s_z,eps_t_r,eps_t_t,eps_p_r,eps_p_t,phi,kappa)
template <class T, class TMEM>
void TPMRSRKSolver<T,TMEM>::PrintSecondaryVariables(std::ostream &out){
    int n_data = m_memory.size();
    int n_var = 10;
    TPZFMatrix<REAL> s_data(n_data,n_var);
    for (int i = 0; i < n_data; i++) {
        REAL r       = m_r_y(i,0);
        REAL s_r     = m_memory[i].GetSigma_n().XX();
        REAL s_t     = m_memory[i].GetSigma_n().YY();
        REAL s_z     = m_memory[i].GetSigma_n().ZZ();
        REAL eps_t_r = m_memory[i].GetPlasticState_n().m_eps_t.XX();
        REAL eps_t_t = m_memory[i].GetPlasticState_n().m_eps_t.YY();
        REAL eps_p_r = m_memory[i].GetPlasticState_n().m_eps_p.XX();
        REAL eps_p_t = m_memory[i].GetPlasticState_n().m_eps_p.YY();
        REAL phi     = m_memory[i].phi_n();
        REAL kappa   = m_memory[i].kappa_n();
        s_data(i,0) = r;
        s_data(i,1) = s_r;
        s_data(i,2) = s_t;
        s_data(i,3) = s_z;
        s_data(i,4) = eps_t_r;
        s_data(i,5) = eps_t_t;
        s_data(i,6) = eps_p_r;
        s_data(i,7) = eps_p_t;
        s_data(i,8) = phi;
        s_data(i,9) = kappa;
    }
    s_data.Print("rksdata = ",out,EMathematicaInput);
}


template <class T, class TMEM>
REAL TPMRSRKSolver<T,TMEM>::lambda(int i){
    return m_lambda[i];
}

template <class T, class TMEM>
REAL TPMRSRKSolver<T,TMEM>::mu(int i){
    return m_mu[i];
}

template <class T, class TMEM>
REAL TPMRSRKSolver<T,TMEM>::K(REAL & lambda, REAL & mu){
    return (lambda + (2.0/3.0)*mu);
}

template <class T, class TMEM>
REAL TPMRSRKSolver<T,TMEM>::Alpha(REAL & K){
    DebugStop();
}

template <class T, class TMEM>
TPZTensor<REAL> TPMRSRKSolver<T,TMEM>::Epsilon(int i, REAL & r, std::vector<REAL> & y){
    
    REAL l = this->lambda(i);
    REAL mu = this->mu(i);
    
    REAL ur = y[0];
    REAL sr = y[1];
    REAL sigmar_0 = m_memory[i].GetSigma_0().XX();
    REAL d_eps_t_rr = (-r*sigmar_0+r*sr-l*ur)/(r*(l+2.0*mu));
    REAL d_eps_t_tt = ur/r;
    
    TPZTensor<REAL> eps_rr = m_memory[i].GetPlasticState_0().m_eps_t;
    TPZTensor<REAL> eps;
    eps.Zero();
    eps.XX() = d_eps_t_rr;
    eps.YY() = d_eps_t_tt;
    eps += eps_rr;
    return eps;
}

template <class T, class TMEM>
TPZTensor<REAL> TPMRSRKSolver<T,TMEM>::Sigma(int i, TPZTensor<REAL> & epsilon, TPZFMatrix<REAL> * Dep){
    TPZTensor<REAL> sigma;
    T plastic_integrator(m_plastic_integrator);
    plastic_integrator.SetState(m_memory[i].GetPlasticState());
    plastic_integrator.ApplyStrainComputeSigma(epsilon,sigma,Dep);
    
    if(m_accept_solution_Q){
        m_memory[i].GetPlasticState_n() = plastic_integrator.GetState();
        m_memory[i].SetSigma_n(sigma);
        if (!IsZero(plastic_integrator.GetState().m_eps_p.Norm())) {
            std::cout << "Plasticity detected in point = " << i <<std::endl;
        }
    }
    return sigma;
}

template <class T, class TMEM>
void TPMRSRKSolver<T,TMEM>::ReconstructAndAcceptPoint(int i, REAL & r, std::vector<REAL> & y, bool initial_data_Q){
    
    REAL last_l,last_mu;

    TPZTensor<REAL> last_eps_p;
    if(initial_data_Q){
        
//        m_memory[i].GetPlasticState().m_eps_t = m_eps_p;
//        m_memory[i].GetPlasticState().m_eps_t += m_eps_e;
//        m_memory[i].GetPlasticState().m_eps_p = m_eps_p;
        last_eps_p = m_eps_p;
        
        last_l = this->lambda(i);
        last_mu = this->mu(i);
    }else{
//        m_memory[i].GetPlasticState().m_eps_p = m_memory[i-1].GetPlasticState().m_eps_p;
        last_eps_p = m_memory[i-1].GetPlasticState_n().m_eps_p;
        
        last_l = this->lambda(i-1);
        last_mu = this->mu(i-1);
    }
    REAL l = last_l;
    REAL mu = last_mu;
    
    
    bool check_Q = false;
    REAL tol = 1.0e-4;
    int n_iterations = 0;
    m_accept_solution_Q = false;
    REAL l_n = last_l;
    REAL mu_n = last_mu;
    
    for (int k = 0; k < n_iterations; k++) {
        
        m_lambda[i] = l;
        m_mu[i] = mu;
        TPZTensor<REAL> epsilon = Epsilon(i,r,y)  + last_eps_p;
        TPZFNMatrix<36,REAL> Dep(6,6,0.0);
        TPZTensor<REAL> sigma   = Sigma(i,epsilon,&Dep);
//        y[1] = sigma.XX();

        /// lamé data correction
        REAL K_ep_xx = (Dep(0,0) + Dep(3,0) + Dep(5,0))/3.0;
        REAL K_ep_yy = (Dep(0,3) + Dep(3,3) + Dep(5,3))/3.0;
        REAL K_ep_zz = (Dep(0,5) + Dep(3,5) + Dep(5,5))/3.0;
        REAL Kep = (K_ep_xx + K_ep_yy + K_ep_zz) / 3.0;
        mu_n = 0.5*(Dep(1,1)+Dep(2,2)+Dep(4,4))/3.0;
        l_n = Kep - (2.0/3.0)*mu;
        
        REAL error_l  = fabs(l - l_n);
        REAL error_mu  = fabs(mu - mu_n);
        check_Q = error_l < tol && error_mu < tol;
        if (check_Q) {
            break;
        }
        
        mu  = mu_n;
        l   = l_n;
        
    }
    
    m_accept_solution_Q = true;
    /// update secondary variables
    m_lambda[i] = l;
    m_mu[i] = mu;
    TPZTensor<REAL> epsilon = Epsilon(i,r,y) + last_eps_p;
    TPZFNMatrix<36,REAL> Dep(6,6,0.0);
    TPZTensor<REAL> sigma   = Sigma(i,epsilon,&Dep);
    
    /// updating
    REAL error = y[1] - sigma.XX();
    std::cout << "Difference in s_rr = " << error << std::endl;
//    y[1] = sigma.XX();
    
    /// lamé data correction
    REAL K_ep_xx = (Dep(0,0) + Dep(3,0) + Dep(5,0))/3.0;
    REAL K_ep_yy = (Dep(0,3) + Dep(3,3) + Dep(5,3))/3.0;
    REAL K_ep_zz = (Dep(0,5) + Dep(3,5) + Dep(5,5))/3.0;
    REAL Kep = (K_ep_xx + K_ep_yy + K_ep_zz) / 3.0;
    mu = 0.5*(Dep(1,1)+Dep(2,2)+Dep(4,4))/3.0;
    l = Kep - (2.0/3.0)*mu;
    
    REAL Kdr_ep = Kep;
    REAL alpha = 1.0 - Kdr_ep/m_K_s;
    m_memory[i].SetAlpha(alpha);

    
    REAL phi = Porosity(i,r,y);
    REAL kappa = Permeability(i,phi);
    m_memory[i].Setphi_n(phi);
    m_memory[i].Setkappa_n(kappa);
    m_lambda[i] = l;
    m_mu[i] = mu;
    m_accept_solution_Q = false;
    
}

template <class T, class TMEM>
REAL TPMRSRKSolver<T,TMEM>::Porosity(int i, REAL & r, std::vector<REAL> & y){
 
    /// Reconstruction of epsilon, sigma and elastoplastic parameters
    
    TPZTensor<REAL> epsilon = Epsilon(i,r,y);
    
    REAL l = this->lambda(i);
    REAL mu = this->mu(i);
    REAL Kdr = l + (2.0/3.0)*mu;
    REAL alpha = 1.0 - Kdr/m_K_s;
    
    REAL pr = y[2];
    
    REAL phi_0 = m_memory[i].phi_0();
    REAL p_0 = m_memory[i].p_0();
    REAL phi = phi_0;
    
    TPZTensor<REAL> epsilon_0 = m_memory[i].GetPlasticState_0().m_eps_t;
    
    /// Apply geomechanic correction
    REAL eps_v_0 = epsilon_0.I1();
    REAL eps_v = epsilon.I1();
    phi += alpha*(eps_v - eps_v_0);
    
    /// Apply pore correction
    REAL S = (1.0-alpha)*(alpha-phi_0)/Kdr;
    phi += S*(pr - p_0);
    return phi;
}

template <class T, class TMEM>
REAL TPMRSRKSolver<T,TMEM>::Permeability(int i, REAL & phi){

    REAL s = 1.0e6;
    REAL phi_0 = m_memory[i].phi_0();
    REAL kappa_0 = m_memory[i].kappa_0()*s;
    REAL kappa, dkappa_dphi;
    
    TPMRSKappaParameters m_model = GetKappaParameters();
    m_model.Permeability(kappa, dkappa_dphi, kappa_0, phi, phi_0);
    return kappa;
    
}

