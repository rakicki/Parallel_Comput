#ifndef COOMATRIX_HPP
#define COOMATRIX_HPP

#include <cassert>
#include <iostream>
#include <memory>
#include <vector>
#include <utility>
#include <algorithm>

template <typename ValueType>
bool less_col(const std::tuple<std::size_t, std::size_t, ValueType>& a,
	      const std::tuple<std::size_t, std::size_t, ValueType>& b){
  return ( (std::get<1>(a) <  std::get<1>(b)) ||
	   ( (std::get<1>(a) == std::get<1>(b)) &&
	     (std::get<0>(a) <  std::get<0>(b)) ) );}

template <typename ValueType>
bool less_row(const std::tuple<std::size_t, std::size_t, ValueType>& a,
	      const std::tuple<std::size_t, std::size_t, ValueType>& b){
  return ( (std::get<0>(a) <  std::get<0>(b)) ||
	   ( (std::get<0>(a) == std::get<0>(b)) &&
	     (std::get<1>(a) <  std::get<1>(b)) ) );}


template <typename VALUE_TYPE>
class CooMatrix{
  
public:
  
  using ValueType      = VALUE_TYPE;
  using ThisType       = CooMatrix<ValueType>;  
  using ItemType       = std::tuple<std::size_t,std::size_t,ValueType>;
  using ComparisonType = std::function<bool(const ItemType&,const ItemType&)>;
  using ContainerType  = std::vector<ItemType>;
  
  CooMatrix(const std::size_t& nr0 = 0,
	    const std::size_t& nc0 = 0):
    nr(nr0), nc(nc0),
    data_ptr(std::make_shared<ContainerType>()) {};
  
  CooMatrix(const CooMatrix&) = default;
  CooMatrix(CooMatrix&&)      = default;

  template <typename OtherValueType>
  CooMatrix(const CooMatrix<OtherValueType>& m):
    nr(NbRow(m)), nc(NbCol(m)),
    data_ptr(std::make_shared<ContainerType>()) {
    const auto& m_data = GetData(m); 
    data_ptr->clear();
    data_ptr->assign(m_data.begin(), m_data.end());
  }

  CooMatrix& operator=(const CooMatrix&) = default; 
  CooMatrix& operator=(CooMatrix&&)      = default;   

  template <typename OtherValueType>
  CooMatrix& operator=(const CooMatrix<OtherValueType>& m){
    nr = NbRow(m);
    nc = NbCol(m);
    const auto& m_data = GetData(m); 
    data_ptr->clear();
    data_ptr->assign(m_data.begin(), m_data.end());
    return *this;
  }
  
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
    
  typedef typename ContainerType::const_iterator const_iterator;
  typedef typename ContainerType::iterator             iterator;
  iterator       begin ()       {return data_ptr->begin(); }
  iterator       end   ()       {return data_ptr->end();   }
  const_iterator begin () const {return data_ptr->cbegin();}
  const_iterator end   () const {return data_ptr->cend();  }  
  const_iterator cbegin() const {return data_ptr->cbegin();}
  const_iterator cend  () const {return data_ptr->cend();  }  
  
  void reserve(const std::size_t& new_size){
    data_ptr->reserve(new_size);}
  
  void push_back(const ItemType& jkv){
    assert( (std::get<0>(jkv)<nr) && (std::get<1>(jkv)<nc) );
    data_ptr->push_back(jkv);}

  void push_back(const ContainerType& new_data){
    data_ptr->insert(data_ptr->end(),new_data.begin(),new_data.end());}
  
  void push_back(const std::size_t& j,
		 const std::size_t& k,
		 const ValueType& v){
    assert( (j<nr) && (k<nc) );
    data_ptr->emplace_back(j,k,v);}
  
  template <typename OtherValueType> auto
  operator()(const std::vector<OtherValueType>& x) const {
    assert(x.size()==nc);
    using CommonType = std::common_type_t<ValueType,OtherValueType>;
    std::vector<CommonType> y(nr,CommonType());
    for(const auto& [j,k,m_jk]:*data_ptr){
      y[j]+=m_jk*x[k];}
    return y;
  }
  
  template <typename ValueType1, typename ValueType2>
  auto& operator()(const std::vector<ValueType1>& x,
	           std::vector<ValueType2>& y) const {
    assert(x.size()==nc);
    assert(y.size()==nr);
    for(const auto& [j,k,m_jk]:*data_ptr){
      y[j]+=m_jk*x[k];}
    return y;
  }

  template <typename InputValueType, typename OutputValueType>
  auto operator()(const InputValueType* u,
		  OutputValueType* v) const {
    const auto& m = *this;
    for(std::size_t j=0; j<NbRow(m); ++j){v[j]=0.;}
    for(const auto& [j,k,mjk]:m){v[j]+=mjk*u[k];}
    return v;
  }
  
  template <typename OtherValueType>
  auto operator*(const std::vector<OtherValueType>& x) const {
    assert(x.size()==nc);
    return (*this)(x);}
  
  friend DenseMatrix<ValueType>
  MakeDense(const CooMatrix& m){
    DenseMatrix<ValueType> new_m(NbRow(m),NbCol(m));
    for(const auto& [j,k,m_jk]: m){new_m[j,k]+=m_jk;}
    return new_m;
  }
  
  friend std::ostream&
  operator<<(std::ostream& o, const CooMatrix& m){
    return o << MakeDense(m);}

  friend std::size_t Nnz(ThisType m){
    m.sort(); return m.data_ptr->size();}

  void sort(const ComparisonType& comp = less_row<ValueType>) const {
    if(!data_ptr->empty()){
      
      std::sort(data_ptr->begin(), data_ptr->end(), comp);
      ContainerType new_data;
      new_data.reserve(data_ptr->size());
      new_data.push_back((*data_ptr)[0]);
      for(auto it0=data_ptr->begin(),
	    it1=std::next(data_ptr->begin());
	  it1!=data_ptr->end(); ++it0, ++it1){
	if(std::get<0>(*it0)==std::get<0>(*it1) &&
	   std::get<1>(*it0)==std::get<1>(*it1) ){
	  std::get<2>(new_data.back()) += std::get<2>(*it1);}
	else{new_data.push_back(*it1);}
      }    
      std::swap(*data_ptr,new_data);
    }    
  }
  
  auto T() const {
    ThisType mt(nc,nr);
    mt.reserve(data_ptr->size());
    for(const auto& [j,k,v]:*data_ptr){mt.push_back(k,j,v);}
    return mt;
  }
  
  auto H() const {
    ThisType mt(nc,nr);
    mt.reserve(data_ptr->size());
    for(const auto& [j,k,v]:*data_ptr){
      mt.push_back(k,j,Conj(v));}
    return mt;
  }
    
  friend double
  Norm(ThisType m){
    m.sort(); double nrm=0.;
    for(const auto& [j,k,v]:m){
      nrm+=std::abs(v)*std::abs(v);}
    return std::sqrt(nrm);
  }

  template <typename OtherValueType>
  friend bool
  Close(const CooMatrix<OtherValueType>& m1,
	const ThisType& m2,
	const double& tol = 1.e-10){
    return Norm(m1-m2)<tol;}

  template <typename OtherValueType>
  ThisType& operator+=(const CooMatrix<OtherValueType>& m){
    assert((nr==NbRow(m)) && (nc==NbCol(m)));
    data_ptr->insert(data_ptr->end(),m.begin(),m.end());
    sort();
    return *this;
  }
  
  template <typename OtherValueType>
  friend auto
  operator+(const ThisType& m1,
	    const CooMatrix<OtherValueType>& m2){
    using CommonType = std::common_type_t<ValueType,OtherValueType>; 
    CooMatrix<CommonType> m3(NbRow(m1),NbCol(m1));
    return (m3+=m1)+=m2;
  }
  
  template <typename OtherValueType>
  ThisType& operator-=(const CooMatrix<OtherValueType>& m){
    assert((nr==NbRow(m)) && (nc==NbCol(m)));    
    data_ptr->reserve(data_ptr->size()+m.data_ptr->size());
    for(const auto& [j,k,v]:m){push_back(j,k,-v);}
    sort();
    return *this;
  }
  
  template <typename OtherValueType>
  friend auto
  operator-(const ThisType& m1,
	    const CooMatrix<OtherValueType>& m2){
    using CommonType = std::common_type_t<ValueType,OtherValueType>; 
    CooMatrix<CommonType> m3(NbRow(m1),NbCol(m1));    
    return (m3+=m1)-=m2;
  }

  template <typename S>
  requires std::same_as<S,double> || std::same_as<S,std::complex<double>>
  auto& operator*=(const S& a){
    for(auto& [j,k,v]:*data_ptr){v*=a;}
    return *this;
  }
  
  template <typename S>
  requires std::same_as<S,double> || std::same_as<S,std::complex<double>>
  friend auto operator*(const S& a, const ThisType& m){
    using CommonType = std::common_type_t<ValueType,S>;
    CooMatrix<CommonType> m2(NbRow(m),NbCol(m));
    return (m2+=m)*=a;
  }

  template <typename OtherValueType>
  auto T(const std::vector<OtherValueType>& u) const {
    assert( nr==u.size() );
    const auto& m = *this;
    using CommonType = std::common_type_t<ValueType,OtherValueType>;
    std::vector<CommonType> v(nc,CommonType());
    for(const auto& [j,k,mjk]:m){v[k]+=mjk*u[j];}
    return v;
  }
    
  template <typename OtherValueType>
  friend auto
  operator*(const ThisType& lhs,
	    const CooMatrix<OtherValueType>& rhs){
    assert( NbCol(lhs)==NbRow(rhs) );
    std::size_t nn = NbCol(lhs);
    
    //######################################//    
    auto data0 = GetData(lhs);

    std::sort(data0.begin(), data0.end(), less_col<ValueType>);
    std::vector<std::size_t> col(nn+1,0);    
    std::size_t c=0;
    for(std::size_t p=0; p<data0.size(); ++p){
      const auto& [j,k,v] = data0[p];
      while(c<k){col[++c]=p;}
    }
    while(c<nn){col[++c]=data0.size();}
    
    //######################################//    
    auto data1 = GetData(rhs);

    std::sort(data1.begin(), data1.end(), less_row<OtherValueType>);
    std::vector<std::size_t> row(nn+1,0);
    std::size_t r=0;
    for(std::size_t p=0; p<data1.size(); ++p){
      const auto& [j,k,v] = data1[p];
      while(r<j){row[++r]=p;}
    }
    while(r<nn){row[++r]=data1.size();}
    
    //######################################//    
    using CommonType = std::common_type_t<ValueType,OtherValueType>;
    CooMatrix<CommonType> new_mat(NbRow(lhs),NbCol(rhs));
    for(std::size_t p=0; p<nn; ++p){
      for(std::size_t col_p=col[p]; col_p<col[p+1]; ++col_p){
	for(std::size_t row_p=row[p]; row_p<row[p+1]; ++row_p){
	  const auto& [j0,k0,v0] = data0[col_p];
	  const auto& [j1,k1,v1] = data1[row_p];
	  new_mat.push_back(j0,k1,v0*v1);
	}
      }
    }

    //######################################//    
    new_mat.sort();
    return new_mat;    
  }
  
private:
  
  //Data members
  std::size_t                       nr,nc;
  std::shared_ptr<ContainerType> data_ptr;
  
};


template <typename Pattern> CooMatrix<double> 
BooleanMatrix(const std::size_t& nr,
	      const std::size_t& nc,
	      const Pattern& pattern){  
  CooMatrix<double> B(nr,nc);
  for(const auto& [j,k]: pattern){
    B.push_back(j,k,1.);}
  return B;}

CooMatrix<double>
IdentityMatrix(const std::size_t& sz){
  CooMatrix<double> Id(sz,sz);
  for(std::size_t j=0; j<sz; ++j){
    Id.push_back(j,j,1.);}
  return Id;}  


template <typename ValueType>
auto Transpose(const std::vector<CooMatrix<ValueType>>& M){
  std::vector<CooMatrix<ValueType>> MT(M.size());
  for(std::size_t j=0; j<MT.size(); ++j){
    MT[j] = M[j].T();}
  return MT;}

template <typename ValueType>
CooMatrix<ValueType>
Diagonal(const CooMatrix<ValueType>& m){
  CooMatrix<ValueType> D(NbRow(m),NbCol(m));
  for(const auto& [j,k,mjk]:m){
    if(j==k){D.push_back(j,k,mjk);}}
  return D;
}



#endif
