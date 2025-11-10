#ifndef SMALLVECTOR_HPP
#define SMALLVECTOR_HPP

#include <complex>
#include <utility>
#include <initializer_list>
#include <array>
#include <assert.h>

typedef std::complex<double>  cplx;

template <typename T, std::size_t D>
class SmallVector{

public:
  
  SmallVector(){
    for(auto& xj:data){xj=T();}
  }
  
  SmallVector(const SmallVector& x){
    for(std::size_t j=0; j<D; ++j){
      data[j]=x.data[j];}    
  }

  SmallVector(const std::initializer_list<T>& x){
    assert(x.size()==D);
    std::copy(x.begin(),x.end(),std::begin(data));}
  
  SmallVector& operator=(const SmallVector& x){
    for(std::size_t j=0; j<D; ++j){
      data[j]=x.data[j];}
    return *this;
  }

  SmallVector& operator=(const std::initializer_list<T>& x){
    assert(x.size()==D);
    std::copy(x.begin(),x.end(),std::begin(data));
    return *this;
  }
  
  T& operator[](const std::size_t& j){
    assert(j<D);
    return data[j];
  }

  const T& operator[](const std::size_t& j) const {
    assert(j<D);
    return data[j];
  }

  SmallVector& operator+=(const SmallVector& x){
    for(std::size_t j=0; j<D; ++j){
      data[j]+=x.data[j];
    }
    return *this;
  }

  SmallVector& operator-=(const SmallVector& x){
    for(std::size_t j=0; j<D; ++j){
      data[j]-=x.data[j];
    }
    return *this;
  }

  template <typename S> 
  SmallVector& operator*=(const S& alpha){
    for(std::size_t j=0; j<D; ++j){
      data[j]*=alpha;
    }
    return *this;
  }
  
  typedef typename std::array<T,D>::const_iterator const_iterator;
  typedef typename std::array<T,D>::iterator             iterator;
  iterator       begin ()       {return data.begin(); }
  iterator       end   ()       {return data.end();   }
  const_iterator begin () const {return data.cbegin();}
  const_iterator end   () const {return data.cend();  }  
  const_iterator cbegin() const {return data.cbegin();}
  const_iterator cend  () const {return data.cend();  }  

  friend std::ostream& operator<<(std::ostream& o, const SmallVector& x){
    for(const auto& xj:x){o << xj << "\t";} return o;}

  friend std::istream& operator>>(std::istream& i, SmallVector& x){
    for(auto& xj:x){i >> xj;} return i;}
  
  std::size_t size() const {return data.size();}
  
private:

  std::array<T,D> data;
  
};



template <typename T, std::size_t D>
auto operator+(const SmallVector<T,D>& x,
	       const SmallVector<T,D>& y){
  return SmallVector<T,D>(x)+=y;}


template <typename T, std::size_t D>
auto operator-(const SmallVector<T,D>& x,
	       const SmallVector<T,D>& y){
  return SmallVector<T,D>(x)-=y;}


template <typename T, std::size_t D, typename S>
auto operator*(const S& alpha,
	       const SmallVector<T,D>& x){
  return SmallVector(x)*=alpha;}


template <std::size_t D>
double operator|(const SmallVector<double,D>& x,
		 const SmallVector<double,D>& y){
  double p = 0.;
  for(std::size_t j=0; j<D; ++j){
    p+= x[j]*y[j];}
  return p;
}


template <std::size_t D>
cplx operator|(const SmallVector<cplx,D>& x,
	       const SmallVector<cplx,D>& y){
  cplx p = 0.;
  for(std::size_t j=0; j<D; ++j){
    p += x[j]*std::conj(y[j]);}
  return p;
}


template <typename T, std::size_t D>
double Norm(const SmallVector<T,D>& x) {
  return std::sqrt( std::abs( x|x ) );}


template <typename T, std::size_t D>
void Normalize(SmallVector<T,D>& x){
  x *= 1./Norm(x);}



auto
VProd(const SmallVector<double,3>& v,
      const SmallVector<double,3>& w){
  SmallVector<double,3> u;
  u[0] = v[1]*w[2] - v[2]*w[1];
  u[1] = v[2]*w[0] - v[0]*w[2];
  u[2] = v[0]*w[1] - v[1]*w[0];
  return u;
}


template <typename T, std::size_t D>
bool Close(const SmallVector<T,D>& x,
	   const SmallVector<T,D>& y,
	   const double& tol = 1e-10){
  return Norm(x-y)<tol;}


typedef SmallVector<double,1> R1;
typedef SmallVector<double,2> R2;
typedef SmallVector<double,3> R3;
typedef SmallVector<double,4> R4;
typedef SmallVector<double,5> R5;
typedef SmallVector<double,6> R6;
typedef SmallVector<double,7> R7;
typedef SmallVector<double,8> R8;
typedef SmallVector<double,9> R9;

typedef SmallVector<cplx,1>   C1;
typedef SmallVector<cplx,2>   C2;
typedef SmallVector<cplx,3>   C3;
typedef SmallVector<cplx,4>   C4;
typedef SmallVector<cplx,5>   C5;
typedef SmallVector<cplx,6>   C6;
typedef SmallVector<cplx,7>   C7;
typedef SmallVector<cplx,8>   C8;
typedef SmallVector<cplx,9>   C9;

#endif
