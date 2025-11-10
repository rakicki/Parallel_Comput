#ifndef DENSEMATRIX_HPP
#define DENSEMATRIX_HPP

#include <cassert>
#include <iostream>
#include <iomanip>
#include <memory>
#include <vector>
#include <functional>
#include <type_traits>
  
template <typename VALUE_TYPE>
class DenseMatrix{
  
public:
  
  using ValueType      = VALUE_TYPE;
  using ThisType       = DenseMatrix<ValueType>;
  using ItemType       = ValueType;
  using ContainerType  = std::vector<ItemType>;

  //Member functions
  DenseMatrix(const std::size_t& nr0 = 0,
	      const std::size_t& nc0 = 0):  
    nr(nr0), nc(nc0),
    data_ptr(std::make_shared<ContainerType>(nr0*nc0,ValueType())) {};

  template <typename OtherValueType>
  DenseMatrix(const DenseMatrix<OtherValueType>& m):
    nr(NbRow(m)), nc(NbCol(m)),
    data_ptr(std::make_shared<ContainerType>()) {
    const auto& m_data = GetData(m); 
    data_ptr->clear();
    data_ptr->assign(m_data.begin(), m_data.end());
  }
  
  DenseMatrix(const DenseMatrix&)            = default;
  DenseMatrix(DenseMatrix&&)                 = default;
  DenseMatrix& operator=(const DenseMatrix&) = default;
  DenseMatrix& operator=(DenseMatrix&&)      = default;

  template <typename OtherValueType>
  DenseMatrix& operator=(const DenseMatrix<OtherValueType>& m){
    nr = NbRow(m);
    nc = NbCol(m);
    const auto& m_data = GetData(m); 
    data_ptr->clear();
    data_ptr->assign(m_data.begin(), m_data.end());
    return *this;}  
  
  ValueType&
  operator[](const std::size_t& j, const std::size_t& k){
    assert(j<nr && k<nc);
    return (*data_ptr)[j*nc+k];}

  const ValueType&
  operator[](const std::size_t& j, const std::size_t& k) const {
    assert(j<nr && k<nc);    
    return (*data_ptr)[j*nc+k];}

  friend std::size_t
  NbRow(const ThisType& m){return m.nr;}

  friend std::size_t
  NbCol(const ThisType& m){return m.nc;}
  
  friend const ContainerType&
  GetData(const ThisType& m) {return *m.data_ptr;}

  friend ContainerType&
  GetData(ThisType& m) {return *m.data_ptr;}

  std::size_t
  use_count() const {return data_ptr.use_count();}
  
  friend ThisType
  Copy(const ThisType& m){
    ThisType new_m(NbRow(m),NbCol(m));
    GetData(new_m) = GetData(m);
    return new_m;
  }
  
  friend double
  Norm(const ThisType& m){return Norm(GetData(m));}
  
  template <typename OtherValueType>
  friend double
  Close(const DenseMatrix<OtherValueType>& m1,
	const ThisType& m2,
	const double& tol = 1.e-10){
    return Norm(m1-m2)<tol;}
  
  friend std::ostream&
  operator<<(std::ostream& o,
	     const ThisType&  m){
    for(std::size_t j=0; j<NbRow(m); ++j){
      for(std::size_t k=0; k<NbCol(m); ++k){
	o << std::setw(10) << std::showpos << std::left;
	o << m[j,k] << "\t";}
      o << "\n";
    }
    return o;
  }

  template<typename InputValueType>
  auto operator()(const std::vector<InputValueType>& x) const {
    assert(x.size()==nc);
    using OutputValueType = std::common_type_t<ValueType,InputValueType>;
    std::vector<OutputValueType> y(nr,OutputValueType());
    for(std::size_t k=0; k<data_ptr->size(); ++k){
      y[k/nc] += (*data_ptr)[k]*x[k%nc];}
    return y;
  }
  
  template <typename InputValueType, typename OutputValueType>
  auto& operator()(const std::vector<InputValueType>& x,
		   std::vector<OutputValueType>& y) const {
    assert(x.size()==nc);
    assert(y.size()==nr);    
    for(std::size_t k=0; k<data_ptr->size(); ++k){
      y[k/nc] += (*data_ptr)[k]*x[k%nc];}
    return y; 
  }
  
  template <typename InputValueType, typename OutputValueType>
  auto operator()(const InputValueType* u,
		  OutputValueType* v) const {
    const auto& m = *this;
    for(std::size_t j=0; j<NbRow(m); ++j){v[j]=0.;}    
    for(std::size_t j=0; j<m.nr; ++j){
      for(std::size_t k=0; k<m.nc; ++k){
	v[j]+=m[j,k]*u[k];}}
    return v;
  }

  template <typename InputValueType>
  auto operator*(const std::vector<InputValueType>& x) const {
    assert(x.size()==nc);
    return (*this)(x);}
  
  auto T(){
    const auto& m = *this;
    ThisType mt(m.nc,m.nr);
    for(std::size_t j=0; j<mt.nr; ++j){
      for(std::size_t k=0; k<mt.nc; ++k){
	mt[j,k]=m[k,j];}}
    return mt;
  }
  
  auto H(){
    const auto& m = *this;
    ThisType mt(m.nc,m.nr);
    for(int j=0; j<mt.nr; ++j){
      for(int k=0; k<mt.nc; ++k){
	mt[j,k]=Conj(m[k,j]);}}
    return mt;
  }
  
  // Additions and substractions
  template <typename S>
  ThisType& operator+=(const DenseMatrix<S>& m){
    assert(nr==NbRow(m) && nc==NbCol(m));
    GetData(*this)+=GetData(m);
    return *this;
  }
  
  template <typename S>
  ThisType& operator-=(const DenseMatrix<S>& m){
    assert(nr==m.nr && nc==m.nc);
    GetData(*this)-=GetData(m);
    return *this;}

  template <typename S>
  friend auto
  operator+(const ThisType& u,
	    const DenseMatrix<S>& v){
    DenseMatrix<std::common_type_t<ValueType,S>> w(NbRow(u),NbCol(u));
    return (w+=u)+=v;}

  template <typename S>
  friend auto
  operator-(const ThisType& u,
	    const DenseMatrix<S>& v){
    DenseMatrix<std::common_type_t<ValueType,S>> w(NbRow(u),NbCol(u));
    return (w+=u)-=v;}

  //Scalar-matrix product
  template <typename S>
  requires std::same_as<S,double> || std::same_as<S,cplx>
  ThisType& operator*=(const S& a){
    (*data_ptr)=a*(*data_ptr); return *this;}

  template <typename S>
  requires std::same_as<S,double> || std::same_as<S,cplx>
  friend auto operator*(const S& a, const ThisType& m){
    DenseMatrix<std::common_type_t<ValueType,S>> m2(NbRow(m),NbCol(m));
    return (m2+=m)*=a;}
    
  template <typename S1, typename S2> friend void
  MvAdd(const ThisType& m, S2* v, const S1* u){
    for(std::size_t j=0; j<m.nr; ++j){
      for(std::size_t k=0; k<m.nc; ++k){
	v[j]+=m[j,k]*u[k];}}}
  
  template <typename InputValueType>
  auto T(const std::vector<InputValueType>& v) const {
    assert( nr==v.size() );
    const auto& m = *this;
    using OutputValueType = std::common_type_t<ValueType,InputValueType>; 
    std::vector<OutputValueType> w(nc,OutputValueType());
    for(std::size_t j=0; j<nr; ++j){
      for(std::size_t k=0; k<nc; ++k){
	w[k]+=m[j,k]*v[j];}}
    return w;}

  //Matrix-Matrix product
  template <typename S>
  friend auto
  operator*(const DenseMatrix<ValueType>& m1,
	    const DenseMatrix<S>& m2){
    assert(NbCol(m1)==NbRow(m2));
    DenseMatrix<std::common_type_t<ValueType,S>> m3(m1.nr,NbCol(m2));
    for(std::size_t j=0; j<NbRow(m1); ++j){
      for(std::size_t k=0; k<NbCol(m2); ++k){
	for(std::size_t l=0; l<NbCol(m1); ++l){	
	  m3[j,k] += m1[j,l]*m2[l,k];
	}
      }
    }
    return m3;}
  
  template <typename S>
  auto T(const DenseMatrix<S>& m2) const {
    const auto& m1 = *this;
    assert(NbRow(m1)==NbRow(m2));
    DenseMatrix<std::common_type_t<ValueType,S>> m3(m1.nc,NbCol(m2));
    for(std::size_t j=0; j<NbCol(m1); ++j){
      for(std::size_t k=0; k<NbCol(m2); ++k){
	for(std::size_t l=0; l<NbRow(m1); ++l){	
	  m3[j,k] += m1[l,j]*m2[l,k];
	}
      }
    }
    return m3;
  }
  
private:
  
  //Data members
  std::size_t                       nr,nc;
  std::shared_ptr<ContainerType> data_ptr;
  
};

#endif
