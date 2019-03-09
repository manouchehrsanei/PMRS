//
//  TPMRSAxisymRKSolver.cpp
//  PMRS
//
//  Created by Manouchehr Sanei on 3/8/19.
//
//

#include "TPMRSAxisymRKSolver.h"

/// Default constructor
TPMRSAxisymRKSolver::TPMRSAxisymRKSolver(){
    m_simulation_data = NULL;
    m_dimension       = 0;
}

/// Destructor
TPMRSAxisymRKSolver::~TPMRSAxisymRKSolver(){
    
}

/// Copy constructor
TPMRSAxisymRKSolver::TPMRSAxisymRKSolver(const TPMRSAxisymRKSolver & other){
    m_simulation_data       = other.m_simulation_data;
    m_dimension             = other.m_dimension;
}


const std::string TPMRSAxisymRKSolver::Name() const{
    return "TPMRSAxisymRKSolver";
}

void TPMRSAxisymRKSolver::Print(std::ostream &out) const{
    out << Name() << std::endl;
    out << "m_dimension = " << m_dimension << std::endl;

}


REAL TPMRSAxisymRKSolver::DisplacementAnalytic(REAL &r, REAL &young, REAL &nu, REAL &biotcof, REAL &u_w, REAL &u_r, REAL &sig_0, REAL &p_0, REAL &p_w, REAL &r_o, REAL &r_w){
    
    REAL dis = (biotcof*p_0*pow(r,2)*pow(r_o,2)*log(r) - biotcof*p_w*pow(r,2)*pow(r_o,2)*log(r) - biotcof*p_0*pow(r,2)*pow(r_w,2)*log(r) + biotcof*p_w*pow(r,2)*pow(r_w,2)*log(r) -
           biotcof*p_0*pow(r,2)*pow(r_o,2)*log(r_o) + biotcof*p_w*pow(r,2)*pow(r_o,2)*log(r_o) + biotcof*p_0*pow(r_o,2)*pow(r_w,2)*log(r_o) -
           biotcof*p_w*pow(r_o,2)*pow(r_w,2)*log(r_o) + (pow(r,2)*r_o*u_r*young*log(r_o))/(1 + nu) + (4*nu*pow(r,2)*r_o*u_r*young*log(r_o))/((1 - 2*nu)*(1 + nu)) -
           (r_o*pow(r_w,2)*u_r*young*log(r_o))/(1 + nu) - (4*nu*r_o*pow(r_w,2)*u_r*young*log(r_o))/((1 - 2*nu)*(1 + nu)) - (pow(r,2)*r_w*u_w*young*log(r_o))/(1 + nu) -
           (4*nu*pow(r,2)*r_w*u_w*young*log(r_o))/((1 - 2*nu)*(1 + nu)) + (pow(r_o,2)*r_w*u_w*young*log(r_o))/(1 + nu) +
           (4*nu*pow(r_o,2)*r_w*u_w*young*log(r_o))/((1 - 2*nu)*(1 + nu)) + biotcof*p_0*pow(r,2)*pow(r_w,2)*log(r_w) - biotcof*p_w*pow(r,2)*pow(r_w,2)*log(r_w) -
           biotcof*p_0*pow(r_o,2)*pow(r_w,2)*log(r_w) + biotcof*p_w*pow(r_o,2)*pow(r_w,2)*log(r_w) - (pow(r,2)*r_o*u_r*young*log(r_w))/(1 + nu) -
           (4*nu*pow(r,2)*r_o*u_r*young*log(r_w))/((1 - 2*nu)*(1 + nu)) + (r_o*pow(r_w,2)*u_r*young*log(r_w))/(1 + nu) +
           (4*nu*r_o*pow(r_w,2)*u_r*young*log(r_w))/((1 - 2*nu)*(1 + nu)) + (pow(r,2)*r_w*u_w*young*log(r_w))/(1 + nu) +
           (4*nu*pow(r,2)*r_w*u_w*young*log(r_w))/((1 - 2*nu)*(1 + nu)) - (pow(r_o,2)*r_w*u_w*young*log(r_w))/(1 + nu) -
           (4*nu*pow(r_o,2)*r_w*u_w*young*log(r_w))/((1 - 2*nu)*(1 + nu)))/
    (2.*r*(pow(r_o,2) - pow(r_w,2))*(young/(2.*(1 + nu)) + (2*nu*young)/((1 - 2*nu)*(1 + nu)))*(log(r_o) - log(r_w)));
    return dis;
    
}


void TPMRSAxisymRKSolver::ConfigureSigma(){
   
    int n_regions = m_simulation_data->NumberOfRegions();
    
    for (int iregion = 0; iregion < n_regions; iregion++)
    {
        std::tuple<TPMRSUndrainedParameters, TPMRSPoroMechParameters, TPMRSPhiParameters,TPMRSKappaParameters,TPMRSPlasticityParameters> chunk =    m_simulation_data->MaterialProps()[iregion];
        
        REAL E,nu;

        TPMRSPoroMechParameters poroperm_parameters(std::get<1>(chunk));
        std::vector<REAL> poroperm_pars = poroperm_parameters.GetParameters();

        
        E = poroperm_pars[0];
        nu = poroperm_pars[1];
        
        REAL alpha     =  0.75;
        REAL u_w       = -0.01;
        REAL u_r       = -0.001;
        REAL sig_0     = -90.0e6;
        REAL p_0       =  20.0e6;
        REAL p_w       =  10.0e6;
        REAL r_o       =  10.0;
        REAL r_w       =  0.1;
        REAL r;
        REAL udsip;
        
        udsip  = DisplacementAnalytic(r, E, nu, alpha, u_w, u_r, sig_0, p_0, p_w, r_o, r_w);
        

//

        
//        TPZFNMatrix<9,REAL> Grad_du, Grad_du_Transpose, delta_e;
//        
//        Grad_du.Transpose(&Grad_du_Transpose);
//        delta_e = Grad_du + Grad_du_Transpose;
//        delta_e *= 0.5;
//        
//
//        TPZFNMatrix<9,REAL> eps_t;
//        
//        eps_t = delta_e;
        
        
        
        TPZElasticResponse ER;
        ER.SetEngineeringData(E, nu);
        
        TPZElasticCriterion Elastic;
        Elastic.SetElasticResponse(ER);
        
        
        TPZTensor<REAL> epsilon_t,sigma;
        sigma.Zero();
        epsilon_t.Zero();
        
        Elastic.ApplyStrainComputeSigma(epsilon_t, sigma);
        
        
    }

}


void ConfigureDisplcement(){
    DebugStop();
}

void ConfigurePressure(){
    DebugStop();
}

void ConfigureFlux(){
    DebugStop();
}

void RungeKuttaSolver(){
    DebugStop();
}
