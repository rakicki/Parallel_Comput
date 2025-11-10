#ifndef VECTOROPS_HPP
#define VECTOROPS_HPP

#include <assert.h>
#include <vector>
#include <iostream>
#include <random>

template <typename T>
std::ostream& operator<<(std::ostream& o,
			 const std::vector<T>& x){
  for(const auto& xj:x){ o << xj << "\t";}
  return o;}

template <typename T,typename S>
auto& operator+=(std::vector<T>& x,
		const std::vector<S> y){
  assert(x.size() == y.size() );
  for(std::size_t j=0; j<x.size(); ++j){x[j]+=y[j];}
  return x;}

template <typename T,typename S>
auto& operator-=(std::vector<T>& x,
		const std::vector<S> y){
  assert(x.size() == y.size() );
  for(std::size_t j=0; j<x.size(); ++j){x[j]-=y[j];}
  return x;}

template <typename T, typename S>
requires std::same_as<S,double> || std::same_as<S,cplx>
auto& operator*=(std::vector<T>& x,
		 const S& alpha){
  for(auto& xj:x){xj*=alpha;}
  return x;}

template <typename T>
auto operator+(const std::vector<T>& x,
	       const std::vector<T>& y){
  std::vector<T> z = x;
  return z+=y; }

template <typename T>
auto operator-(const std::vector<T>& x,
	       const std::vector<T>& y){
  std::vector<T> z = x;
  return z-=y; }

template <typename T, typename S>
requires std::same_as<S,double> || std::same_as<S,cplx>
auto operator*(const S& alpha,
	       const std::vector<T>& x){
  std::vector<T> z = x;
  return z*=alpha;
}

double operator|(const std::vector<double>& x,
		 const std::vector<double>& y){
   assert(x.size()==y.size());
   double p = 0.;
   for(std::size_t j=0; j<x.size(); ++j){
     p += x[j]*y[j];}
   return p;
}

cplx operator|(const std::vector<cplx>& x,
	       const std::vector<cplx>& y){
  assert(x.size()==y.size());
  cplx p = 0.;
  for(std::size_t j=0; j<x.size(); ++j){
    p += x[j]*std::conj(y[j]);}
  return p;
}

template <typename T>
double Norm(const std::vector<T>& x){
  return std::sqrt( std::abs( x|x ) );}

template <typename S, typename T>
bool Close(const std::vector<S>& x,
	   const std::vector<T>& y,
	   const double& tol = 1e-10){
  return Norm(x-y)<tol;}

std::vector<double> RandomVec(const std::size_t& sz){
  std::vector<double> rv(sz);
  std::random_device rd;
  std::default_random_engine e{rd()};
  std::uniform_real_distribution<double> u(-10.0,10.0);
  for(std::size_t j=0; j<sz; ++j){rv[j]=u(e);}
  rv *= 1./Norm(rv);
  return rv;
}




#endif
