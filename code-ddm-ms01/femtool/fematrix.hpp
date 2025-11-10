#ifndef FEMATRIX_HPP
#define FEMATRIX_HPP

#include "coomatrix.hpp"


template <std::size_t DIM>
auto Identity(const FeSpace<DIM>& Vh){
  return IdentityMatrix(dim(Vh));
}

template <std::size_t DIM>
auto Mass(const FeSpace<DIM>& Vh){

  CooMatrix<double> M(dim(Vh),dim(Vh));
  for(const auto& I:Vh){
    double h = Vol(I.elt())/((DIM+1.)*(DIM+2.));
    for(std::size_t j=0; j< dim(I); ++j){
      for(std::size_t k=0; k< dim(I); ++k){
	if(j==k){ M.push_back(I[j],I[k],2.*h); }
	else    { M.push_back(I[j],I[k],   h); }
      }
    }
  }
  M.sort();
  return M;
  
}

template <std::size_t DIM>
auto Stiffness(const FeSpace<DIM>& Vh){
  
  CooMatrix<double> K(dim(Vh),dim(Vh));  
  std::size_t d = FeSpace<DIM>::local_space_dim;
  for(const auto& I:Vh){

    const auto& e = I.elt();
    auto  n = BdNormal(e);
    auto  h = Vol(e);
    
    for(std::size_t j=0; j<dim(I); ++j){
      for(std::size_t k=0; k<dim(I); ++k){

	double Kjk = 0.;	
	Kjk  = (n[j]|n[k])*h;
	Kjk /= ((e[j]-e[(j+1)%d])|n[j]);
	Kjk /= ((e[k]-e[(k+1)%d])|n[k]);
	K.push_back(I[j],I[k],Kjk);
	
      }
    }
  }
  K.sort();
  return K;

}




#endif
