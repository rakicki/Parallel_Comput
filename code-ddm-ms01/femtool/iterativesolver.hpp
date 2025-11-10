#ifndef ITERATIVE_SOLVER_HPP
#define ITERATIVE_SOLVER_HPP

#include <functional>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <type_traits>
#include <vector>

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
  eps2       *= std::abs((b|b));      
  double alpha,beta,pAp;
      std::ofstream file("courbe_conv.dat");
  std::size_t niter = 0;    
  // Function x = (x1,x2) -> cos(omega*x1)
  auto ue    = [](const R3& x){return std::cos(10.*M_PI*x[0]);};
  while( r2>eps2 && niter++<1000 ){
      
    Ap    = A*p;
    pAp   = std::real((Ap|p));    
    alpha = r2/pAp;
    x    += alpha*p;
    r    -= alpha*Ap;    
    r2    = std::pow(Norm(r),2);
    beta  = r2/(alpha*pAp);
    p     = beta*p+r;
    auto err=ue-x;
    file << niter << "; " << (Norm(err))/((ue)) << "\n";
    if((niter%50)==0){
      std::cout << std::left << std::setw(7) << niter << "\t";
      std::cout << Norm(r) << std::endl;
    }
      
  }
    file.close();
  return x;
    
}




#endif
