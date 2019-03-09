//
//  TPMRSRKSolver.cpp
//  PMRS
//
//  Created by Manouchehr Sanei on 3/8/19.
//
//

#include "TPMRSRKSolver.h"

/// Default constructor
TPMRSRKSolver::TPMRSRKSolver(){
    m_simulation_data = NULL;
    m_dimension       = 0;
}

/// Destructor
TPMRSRKSolver::~TPMRSRKSolver(){
    
}

/// Copy constructor
TPMRSRKSolver::TPMRSRKSolver(const TPMRSRKSolver & other){
    DebugStop();
}

/// Assignment constructor
TPMRSRKSolver & TPMRSRKSolver::operator=(const TPMRSRKSolver & other){
    DebugStop();
}

const std::string TPMRSRKSolver::Name() const{
    return "TPMRSRKSolver";
}

void TPMRSRKSolver::Print(std::ostream &out) const{
    out << Name() << std::endl;
    DebugStop();

}

std::vector<REAL> TPMRSRKSolver::f(std::vector<REAL> & y){
    
    
    
}

std::vector<REAL> TPMRSRKSolver::RK2Approximation(std::vector<REAL> & y_n){
    
}


std::vector<REAL> TPMRSRKSolver::RK4Approximation(std::vector<REAL> & y_n){
    
}

void TPMRSRKSolver::ExecuteRKApproximation(){

}

void TPMRSRKSolver::PrintRKApproximation(){
    
}

REAL TPMRSRKSolver::Porosity(std::vector<REAL> & y){
    
}

REAL TPMRSRKSolver::Permeability(REAL & phi){
    
}

REAL TPMRSRKSolver::K(std::vector<REAL> & y){
    
}

REAL TPMRSRKSolver::alpha(std::vector<REAL> & y){
    
}
