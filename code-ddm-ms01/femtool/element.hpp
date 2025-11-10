#ifndef ELEMENT_HPP
#define ELEMENT_HPP

#include <array>
#include <iostream>
#include <algorithm>
#include <assert.h>
#include "smallvector.hpp"

template <std::size_t D>
class Element{

public:

  Element(){data.fill(nullptr);}

  Element(const Element&) = default;

  Element(const R3& x0)
    requires(D==0)
  {data[0]=&x0;}

  Element(const R3& x0, const R3& x1)
    requires(D==1)
  {data[0]=&x0; data[1]=&x1; sort();}

  Element(const R3& x0, const R3& x1, const R3& x2)
    requires(D==2)
  {data[0]=&x0; data[1]=&x1; data[2]=&x2; sort();}

  Element(const R3& x0, const R3& x1, const R3& x2, const R3& x3)
    requires(D==3)
  {data[0]=&x0; data[1]=&x1; data[2]=&x2; data[3]=&x3; sort();}
  
  bool initialized() const {
    return !(data[D]==nullptr);}

  Element& operator=(const Element& e){
    assert(e.initialized());
    data = e.data;
    return *this;
  }

  void push_back(const R3& vtx){
    assert(!initialized());
    for(auto& vtx_ptr:data){
      if(vtx_ptr==nullptr){
	vtx_ptr=&vtx;
	if(initialized()){sort();}
	return;
      }
    }
  }
  
  const R3& operator[](const std::size_t& j) const {
    assert(j<D+1); return *data[j]; }
  
  friend std::ostream& operator<<(std::ostream& o, const Element& e){
    for(std::size_t j=0; j<D+1; ++j){o << e[j] << "\n";}
    return o;}
  
  void sort(){
    assert(initialized());
    std::sort(data.begin(),data.end());}
  
  friend bool operator==(const Element& e1, const Element& e2){
    return (e1.data == e2.data); }
  
  friend bool operator!=(const Element& e1, const Element& e2){
    return (e1.data != e2.data); }
  
  friend bool operator<(const Element& e1, const Element& e2){
    return (e1.data < e2.data); }
  
  friend bool operator>(const Element& e1, const Element& e2){
    return (e1.data < e2.data); }
  
  friend auto Boundary(const Element& e){
    assert(e.initialized());
    std::array<Element<D-1>,D+1> faces;
    for(std::size_t j=0; j<D+1; ++j){
      for(std::size_t k=0; k<D+1; ++k){
	if(k!=D-j){faces[j].push_back(e[k]);}
      }
    }
    return faces;
  }

  void clear(){data.fill(nullptr);}

private:

  using VtxPtr = R3 const *;
  std::array<VtxPtr,D+1> data;

};


using Elt1D = Element<1>;
using Elt2D = Element<2>;
using Elt3D = Element<3>;


template <std::size_t DIM>
R3 Ctr(const Element<DIM>& e){

  R3 ctr;
  for(std::size_t j=0; j<DIM+1; ++j){
    ctr += e[j];}
  ctr *= (1./double(DIM));
  return ctr;

}



//############################//
//     Volume d'un element    //
//############################//
double Vol(const Elt1D& e){
  return Norm(e[1]-e[0]);}

double Vol(const Elt2D& e){
  return 0.5*Norm(VProd(e[1]-e[0],e[2]-e[1]));}

double Vol(const Elt3D& e){
  R3 u1 = e[1]-e[0];
  R3 u2 = e[2]-e[0];
  R3 u3 = e[3]-e[0];  
  return std::abs(( u1 | VProd(u2,u3)))/6.;
}

//############################//
//       Normale au bord      //
//############################//
auto BdNormal(const Elt1D& e){
  std::array<R3,2> n;
  for(std::size_t j=0; j<2; ++j){
    n[j] = e[(j+1)%2]-e[j];
    Normalize(n[j]);
  }  
  return n;
}

auto BdNormal(const Elt2D& e){
  std::array<R3,3> u;
  for(std::size_t j=0; j<3; ++j){
    u[j] = e[j]-e[(j+1)%3];}
  
  std::array<R3,3> n;
  R3 nr;
  for(std::size_t j=0; j<3; ++j){
    nr   = VProd(u[j],u[(j+1)%3]);
    nr  *= 1./Norm(nr);
    n[j] = VProd(nr,u[(j+1)%3]);
    Normalize(n[j]);
    if( (n[j]|u[j])>0 ){ n[j]*= -1.; }
  } 
  return n;
}

auto BdNormal(const Elt3D& e){
  std::array<R3,4> u;
  for(std::size_t j=0; j<4; ++j){
    u[j] = e[j]-e[(j+1)%4];}

  std::array<R3,4> n;
  for(std::size_t j=0; j<4; ++j){
    n[j] = VProd( u[(j+1)%4],u[(j+2)%4] );
    if( (n[j]|u[j])>0 ){n[j]*=-1.;}
    Normalize(n[j]);
  }
  return n;  

}

#endif
