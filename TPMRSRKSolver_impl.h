//
//  TPMRSRKSolver_impl.hpp
//  PMRS
//
//  Created by Omar Durán on 3/9/19.
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
    m_n_state   = 4;
    m_memory.resize(0);
    m_r_y.Resize(0, 0);
    m_lambda.resize(0);
    m_mu.resize(0);
    m_accept_solution_Q = false;
}

template <class T, class TMEM>
TPMRSRKSolver<T,TMEM>::~TPMRSRKSolver(){
    
}

template <class T, class TMEM>
TPMRSRKSolver<T,TMEM>::TPMRSRKSolver(const TPMRSRKSolver & other){
    
    m_n_steps       = other.m_n_steps;
    m_y_0           = other.m_y_0;
    m_re            = other.m_re;
    m_rw            = other.m_rw;
    m_eta           = other.m_eta;
    m_cf            = other.m_cf;
    m_K_s           = other.m_K_s;
    m_dr            = other.m_dr;
    m_default_memory = other.m_default_memory;
    m_is_RK4_Q      = other.m_is_RK4_Q;
    m_n_state       = other.m_n_state;
    m_plastic_integrator = other.m_plastic_integrator;
    m_memory        = other.m_memory;
    m_r_y           = other.m_r_y;
    m_lambda        = other.m_lambda;
    m_mu            = other.m_mu;
    m_accept_solution_Q = other.m_accept_solution_Q;
    
}

template <class T, class TMEM>
TPMRSRKSolver<T,TMEM> & TPMRSRKSolver<T,TMEM>::operator=(const TPMRSRKSolver & other){
    
    /// check for self-assignment
    if(&other == this){
        return *this;
    }
    
    m_n_steps       = other.m_n_steps;
    m_y_0           = other.m_y_0;
    m_re            = other.m_re;
    m_rw            = other.m_rw;
    m_eta           = other.m_eta;
    m_cf            = other.m_cf;
    m_K_s           = other.m_K_s;
    m_dr            = other.m_dr;
    m_default_memory = other.m_default_memory;
    m_is_RK4_Q      = other.m_is_RK4_Q;
    m_n_state       = other.m_n_state;
    m_plastic_integrator = other.m_plastic_integrator;
    m_memory        = other.m_memory;
    m_r_y           = other.m_r_y;
    m_lambda        = other.m_lambda;
    m_mu            = other.m_mu;
    m_accept_solution_Q = other.m_accept_solution_Q;
    
    return *this;
}

template <class T, class TMEM>
const std::string TPMRSRKSolver<T,TMEM>::Name() const{
    return "TPMRSRKSolver";
}

template <class T, class TMEM>
void TPMRSRKSolver<T,TMEM>::Print(std::ostream &out) const{
    
    out << Name() << std::endl;
    out << "m_n_steps = " << m_n_steps << std::endl;
    out << "m_y_0 = " << m_y_0 << std::endl;
    out << "m_re = " << m_re << std::endl;
    out << "m_rw = " << m_rw << std::endl;
    out << "m_eta = " << m_eta << std::endl;
    out << "m_cf = " << m_cf << std::endl;
    out << "m_K_s = " << m_K_s << std::endl;
    out << "m_dr = " << m_dr << std::endl;
    out << "m_default_memory = " << m_default_memory << std::endl;
    out << "m_is_RK4_Q = " << m_is_RK4_Q << std::endl;
    out << "m_n_state = " << m_n_state << std::endl;
    out << "m_plastic_integrator = " << m_plastic_integrator << std::endl;
    out << "m_memory = " << m_memory << std::endl;
    out << "m_r_y = " << m_r_y << std::endl;
    out << "m_lambda = " << m_lambda << std::endl;
    out << "m_mu = " << m_mu << std::endl;
    
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

#define new_RK_Q

template <class T, class TMEM>
std::vector<REAL> TPMRSRKSolver<T,TMEM>::f(int i, REAL & r, std::vector<REAL> & y){
    


#ifdef new_RK_Q
    
    REAL qr = y[3];
    REAL sr_0 = m_memory[i].GetSigma_0().XX();
    REAL st_0 = m_memory[i].GetSigma_0().YY();
    
    REAL phi = Porosity(i,r,y);
    REAL kappa = Permeability(i,phi);
    
    std::vector<REAL> f(4);
    REAL alpha = m_memory[i].Alpha();
    TPZTensor<REAL> eps_t = m_memory[i].GetPlasticState_n().m_eps_t;
    TPZTensor<REAL> sigma = m_memory[i].GetSigma_n();
    
    f[0] = eps_t.XX();
    f[1] = -alpha*(m_eta/kappa)*qr + (sr_0-st_0)/(r) + (sigma.YY()-sigma.XX())/(r);
    f[2] = -(m_eta/kappa)*qr;
    f[3] = -qr/r;
    
#else
    
    REAL ur = y[0];
    REAL sr = y[1];
    REAL qr = y[3];
    REAL sr_0 = m_memory[i].GetSigma_0().XX();
    REAL st_0 = m_memory[i].GetSigma_0().YY();
    
    REAL phi = Porosity(i,r,y);
    REAL kappa = Permeability(i,phi);
    
    std::vector<REAL> f(4);
    
    REAL l = this->lambda(i);
    REAL mu = this->mu(i);
    REAL alpha = m_memory[i].Alpha();

    REAL eps_t_rr = (r*sr-l*ur)/(r*(l+2.0*mu));
    f[0] = eps_t_rr;
    f[1] = -alpha*(m_eta/kappa)*qr + (sr_0-st_0)/(r) + 2.0*mu*((ur/(r*r))-eps_t_rr/r);
    f[2] = -(m_eta/kappa)*qr;
    f[3] = -qr/r;
    
#endif

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
    ReconstructAndAcceptPoint(0,m_re,y);
    AppendTo(0,y);
    
    for (int i = 1; i < n_points; i++) {
        
        REAL r = m_dr*(i-1) + m_re;
        if(m_is_RK4_Q){
            y = EulerApproximation(i-1,r,y);
//            y = RK4Approximation(i-1,r,y);
        }else{
            y = RK2Approximation(i-1,r,y);
        }
        ReconstructAndAcceptPoint(i,r,y);
        AppendTo(i,y);
    }
    
}

template <class T, class TMEM>
void TPMRSRKSolver<T,TMEM>::AppendTo(int i, std::vector<REAL> y){
    REAL r = m_dr*(i) + m_re;
    for (int k = 0; k < m_n_state; k++) {
        m_r_y(i,0) = r;
        m_r_y(i,k+1) = y[k];
    }
}

template <class T, class TMEM>
void TPMRSRKSolver<T,TMEM>::PrintRKApproximation(){
    m_r_y.Print("rkdata = ",std::cout,EMathematicaInput);
}

/// Print the secondary variables (s_r,s_t,eps_t_r,eps_t_t,eps_p_r,eps_p_t,phi,kappa)
template <class T, class TMEM>
void TPMRSRKSolver<T,TMEM>::PrintSecondaryVariables(){
    int n_data = m_memory.size();
    int n_var = 9;
    TPZFMatrix<REAL> s_data(n_data,n_var);
    for (int i = 0; i < n_data; i++) {
        REAL r = m_r_y(i,0);
        REAL s_r = m_memory[i].GetSigma_n().XX();
        REAL s_t = m_memory[i].GetSigma_n().YY();
        REAL eps_t_r = m_memory[i].GetPlasticState_n().m_eps_t.XX();
        REAL eps_t_t = m_memory[i].GetPlasticState_n().m_eps_t.YY();
        REAL eps_p_r = m_memory[i].GetPlasticState_n().m_eps_t.XX();
        REAL eps_p_t = m_memory[i].GetPlasticState_n().m_eps_t.YY();
        REAL phi = m_memory[i].phi_n();
        REAL kappa = m_memory[i].kappa_n();
        s_data(i,0) = r;
        s_data(i,1) = s_r;
        s_data(i,2) = s_t;
        s_data(i,3) = eps_t_r;
        s_data(i,4) = eps_t_t;
        s_data(i,5) = eps_p_r;
        s_data(i,6) = eps_p_t;
        s_data(i,7) = phi;
        s_data(i,8) = kappa;
    }
    s_data.Print("rksdata = ",std::cout,EMathematicaInput);
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
    
    REAL eps_t_rr = (r*sr-l*ur)/(r*(l+2.0*mu));
    REAL eps_t_tt = ur/r;
    
    TPZTensor<REAL> eps;
    eps.Zero();
    eps.XX() = eps_t_rr;
    eps.YY() = eps_t_tt;
    
    return eps;
}

template <class T, class TMEM>
TPZTensor<REAL> TPMRSRKSolver<T,TMEM>::Sigma(int i, TPZTensor<REAL> & epsilon, TPZFMatrix<REAL> * Dep){
    TPZTensor<REAL> sigma;
    T plastic_integrator(m_plastic_integrator);
    plastic_integrator.SetState(m_memory[i].GetPlasticState_n());
    plastic_integrator.ApplyStrainComputeSigma(epsilon,sigma,Dep);
    
    if(m_accept_solution_Q){
        m_memory[i].GetPlasticState_n() = plastic_integrator.GetState();
        m_memory[i].SetSigma_n(sigma);
        if (!IsZero(plastic_integrator.GetState().m_eps_p.Norm())) {
            std::cout << "Plasticity " <<std::endl;
        }
    }
    return sigma;
}

template <class T, class TMEM>
void TPMRSRKSolver<T,TMEM>::ReconstructAndAcceptPoint(int i, REAL & r, std::vector<REAL> & y){
    
#ifdef new_RK_Q
    
    m_accept_solution_Q = true;
    /// update secondary variables
    REAL l = this->lambda(i);
    REAL mu = this->mu(i);
    TPZTensor<REAL> epsilon = Epsilon(i,r,y); /// Resconstructed eps
    TPZTensor<REAL> sigma;
    
    sigma.Zero();
    REAL u_r = y[0];
    REAL s_r = y[1];
    REAL eps_r = epsilon.XX();
    REAL s_z = s_r  - 2.0*mu*eps_r;
    REAL s_t = s_z  + 2.0*mu*u_r/r;
    sigma.XX() = s_r;
    sigma.YY() = s_t;
    sigma.ZZ() = s_z;
    
    TPZFNMatrix<36,REAL> Dep(6,6,0.0);
//    T plastic_integrator(m_plastic_integrator);
//    plastic_integrator.SetState(m_memory[i].GetPlasticState_n());
//    plastic_integrator.ApplyLoad(sigma,epsilon);
    
    /// now update all the variables
    sigma = Sigma(i,epsilon,&Dep);
    l = Dep(0,3);
    mu = Dep(1,1)/2.0;
    REAL Kdr_ep = l + (2.0/3.0)*mu;
    REAL alpha = 1.0 - Kdr_ep/m_K_s;
    REAL phi = Porosity(i,r,y);
    REAL kappa = Permeability(i,phi);
    m_lambda[i] = l;
    m_mu[i] = mu;
    m_memory[i].Setphi_n(phi);
    m_memory[i].Setkappa_n(kappa);
    m_memory[i].SetAlpha(alpha);
    m_accept_solution_Q = false;
    
#else
    
    REAL last_l = this->lambda(i);
    REAL last_mu = this->mu(i);
    
    bool check_Q = false;
    REAL tol = 1.0e-8;
    int n_iterations = 50;
    REAL l = last_l;
    REAL mu = last_mu;
    m_accept_solution_Q = false;
    for (int k = 0; k < n_iterations; k++) {
        m_lambda[i] = l;
        m_mu[i] = mu;
        TPZTensor<REAL> epsilon = Epsilon(i,r,y);
        TPZFNMatrix<36,REAL> Dep(6,6,0.0);
        TPZTensor<REAL> sigma   = Sigma(i,epsilon,&Dep);
        REAL error_l  = fabs(l - Dep(0,3));
        REAL error_mu  = fabs(mu - Dep(1,1)/2.0);
        check_Q = error_l < tol && error_mu < tol;
        if (check_Q) {
            break;
        }
        l = Dep(0,3);
        mu = Dep(1,1)/2.0;
    }

    if(!check_Q){
        m_lambda[i] = last_l;
        m_mu[i] = last_mu;
        std::cout << "TPMRSRKSolver<T,TMEM>:: Nonlinear Process does not converge!" << std::endl;
    }
    
    m_accept_solution_Q = true;
    /// update secondary variables
    TPZTensor<REAL> epsilon = Epsilon(i,r,y);
    TPZFNMatrix<36,REAL> Dep(6,6,0.0);
    TPZTensor<REAL> sigma   = Sigma(i,epsilon,&Dep);
    REAL Kdr_ep = l + (2.0/3.0)*mu;
    REAL alpha = 1.0 - Kdr_ep/m_K_s;
    REAL phi = Porosity(i,r,y);
    REAL kappa = Permeability(i,phi);
    m_memory[i].Setphi_n(phi);
    m_memory[i].Setkappa_n(kappa);
    m_memory[i].SetAlpha(alpha);
    m_accept_solution_Q = false;
    
#endif
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
    phi += alpha*(epsilon.I1() - epsilon_0.I1());
    
    /// Apply pore correction
    REAL S = (1.0-alpha)*(alpha-phi_0)/Kdr;
    phi += S*(pr - p_0);
    return phi;
}

template <class T, class TMEM>
REAL TPMRSRKSolver<T,TMEM>::Permeability(int i, REAL & phi){
    REAL A = 0.0;
    REAL s = 1.0e6;
    REAL phi_0 = m_memory[i].phi_0();
    REAL kappa_0 = m_memory[i].kappa_0()*s;
    REAL kappa = kappa_0*pow(phi/phi_0,A);
    return kappa;
}

