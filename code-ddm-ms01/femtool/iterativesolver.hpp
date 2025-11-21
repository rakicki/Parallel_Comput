#ifndef ITERATIVE_SOLVER_HPP
#define ITERATIVE_SOLVER_HPP

#include <functional>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <type_traits>
#include <vector>
#include "preconditioner.hpp"  // for CholeskyPrec
// #include <femtool.hpp>      // ← remove this from the solver header
#include "coomatrix.hpp"   
std::vector<double>
cgsolve(const CooMatrix<double>&   A,
	const std::vector<double>& b) {
    
  assert((NbCol(A)==NbRow(A)) &&
	 (b.size()==NbCol(A)) );
 
  auto    x   = std::vector<double>(b.size(),0.);
  auto    r   = b-A*x;
  auto    p   = r;
  auto   Ap   = A*p;    
  double r2   = std::pow(Norm(r),2);  
  double eps2 = (1e-8)*r2;
  double eps=1e-6;
  eps2       *= std::abs((b|b));      
  double alpha,beta,pAp;

  Mesh2D Omega;

  // Loading a 2D mesh
  Read(Omega,"exo2h0.0025.mesh");

  // Assembly of a finite element space over Omega
  auto Vh   = FeSpace(Omega);

 
  
  // Function x = (x1,x2) -> cos(omega*x1)
  auto F    = [](const R3& x){return std::cos(10.*M_PI*x[0]);};

  // Manufactured "exact solution" obtained
  // by nodal evaluation of f at the degrees of freedom of Vh
  auto ue   = Vh(F);
  auto err  = ue-x;
  double norm_err = sqrt(( (A(err)|err))/((A(ue)|ue)));
  std::size_t niter = 0;    
  // Function x = (x1,x2) -> cos(omega*x1)

    std::ofstream file("courbe_conv_h0.0025.dat");
  while( norm_err>eps && niter++<2000 ){

    Ap    = A*p;
    pAp   = std::real((Ap|p));    
    alpha = r2/pAp;
    x    += alpha*p;
    r    -= alpha*Ap;    
    r2    = std::pow(Norm(r),2);
    beta  = r2/(alpha*pAp);
    p     = beta*p+r;
    err   = ue - x;
    norm_err = sqrt(( (A(err)|err))/((A(ue)|ue)));
    file << niter << "; " << norm_err  << "\n";
    if((niter%50)==0){
      std::cout << std::left << std::setw(7) << niter << "\t";
      std::cout << Norm(r) << std::endl;
     
   }
      
  }
  std::cout << "CG converged in " << niter << " iterations." << std::endl;
  file.close();
  return x;
    
}

std::vector<double> PCGSolver(const CooMatrix<double>&   A,
	const std::vector<double>& b,const CholeskyPrec& P) {
    
  assert((NbCol(A)==NbRow(A)) &&
   (b.size()==NbCol(A)) );
 
  auto    x   = std::vector<double>(b.size(),0.);
  auto    r   = b-A*x;
  auto    z   = P(r);
  auto    p   = z;
  auto   Ap   = A*p;    
  double rz   = (r|z);  
  double eps2 = (1e-8)*rz;
  double eps=1e-6;
  eps2       *= std::abs((b|b));      
  double alpha,beta,pAp;

  Mesh2D Omega;

  // Loading a 2D mesh
  Read(Omega,"tp1-1.mesh");

  // Assembly of a finite element space over Omega
  auto Vh   = FeSpace(Omega);

 
  
  // Function x = (x1,x2) -> cos(omega*x1)
  auto F    = [](const R3& x){return std::cos(10.*M_PI*x[0]);};

  // Manufactured "exact solution" obtained
  // by nodal evaluation of f at the degrees of freedom of Vh
  auto ue   = Vh(F);
  auto err  = ue-x;
  double norm_err = sqrt(( (A(err)|err))/((A(ue)|ue)));
  std::size_t niter = 0;   
   std::ofstream file("TP2_exo1_h0.025.dat"); 
  while( norm_err>eps && niter++<2000 ){

    Ap    = A*p;
    pAp   = std::real((Ap|p));    
    alpha = rz/pAp;
    x    += alpha*p;
    r    -= alpha*Ap;    
    z     = P(r);
    double rz_new  = (r|z);
    beta  = rz_new/rz;
    p     = beta*p+z;
    rz    = rz_new;
     err   = ue - x;
    norm_err = sqrt(( (A(err)|err))/((A(ue)|ue)));
    file << niter << "; " << norm_err  << "\n";
    if((niter%50)==0){
      std::cout << std::left << std::setw(7) << niter << "\t";
      std::cout << Norm(r) << std::endl;
     
   }
      
  }
  std::cout << "PCG converged in " << niter << " iterations." << std::endl;
  file.close();
  return x;
  }

#endif
