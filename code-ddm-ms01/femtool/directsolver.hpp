#ifndef DIRECT_SOLVER_HPP
#define DIRECT_SOLVER_HPP

#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include "coomatrix.hpp"

template <typename T>
void Copy(const CooMatrix<T>& A,
	  Eigen::SparseMatrix<T>& Ae){

  using NxNxT = Eigen::Triplet<double>;
  std::vector<NxNxT> jkv;
  for(const auto& [j,k,v]:A){
    jkv.push_back(NxNxT(int(j),int(k),v));}
  Ae.setFromTriplets(jkv.begin(),jkv.end());
}

template <typename T>
void Copy(const std::vector<T>& u,
	  Eigen::Matrix<T,Eigen::Dynamic,1>& ue){
  ue.resize(u.size());
  for(std::size_t j=0; j<u.size(); ++j){ue[j] = u[j];}}

template <typename T>
void Copy(const Eigen::Matrix<T,Eigen::Dynamic,1>& ue,
	  std::vector<T>& u){
  u.resize(ue.size());
  for(std::size_t j=0; j<u.size(); ++j){u[j] = ue[j];}}


template <typename VALUE_TYPE>
class InvCooMatrix{

public:
  
  using ValueType       = VALUE_TYPE;
  using ThisType        = InvCooMatrix<ValueType>;  
  using EigenVectorType = Eigen::Matrix<ValueType,Eigen::Dynamic,1>;
  using EigenMatrixType = Eigen::SparseMatrix<ValueType>;
  using EigenOrderType  = Eigen::COLAMDOrdering<int>;
  using ContainerType   = Eigen::SparseLU<EigenMatrixType,
					  EigenOrderType>;

  InvCooMatrix(const CooMatrix<ValueType>& A):
    nr(NbRow(A)), nc(NbCol(A)),
    data_ptr(std::make_shared<ContainerType>()),
    u_ptr(std::make_shared<EigenVectorType>()),
    b_ptr(std::make_shared<EigenVectorType>())
  {
    assert( nr==nc );
    Eigen::SparseMatrix<double> Ae(nr,nc);
    Copy(A,Ae);
    data_ptr->analyzePattern(Ae);
    data_ptr->factorize(Ae);
    u_ptr->resize(nr);
    b_ptr->resize(nc);
  };

  InvCooMatrix()                               = default;
  InvCooMatrix(const InvCooMatrix&)            = default;
  InvCooMatrix(InvCooMatrix&&)                 = default;
  InvCooMatrix& operator=(const InvCooMatrix&) = default; 
  InvCooMatrix& operator=(InvCooMatrix&&)      = default;     
  friend std::size_t
  NbRow(const ThisType& m){return m.nr;}
  
  friend std::size_t
  NbCol(const ThisType& m){return m.nc;}

  auto operator()(const std::vector<ValueType>& b) const {

    auto& ue = *u_ptr;
    auto& be = *b_ptr;
    
    Copy(b,be);
    ue = data_ptr->solve(*b_ptr);
    std::vector<ValueType> u;
    Copy(ue,u);
    
    return u;
  }

  auto operator*(const std::vector<ValueType>& b) const {
    return (*this)(b);}
  
private:
  
  std::size_t                       nr,nc;
  std::shared_ptr<ContainerType> data_ptr;
  std::shared_ptr<EigenVectorType>  u_ptr;
  std::shared_ptr<EigenVectorType>  b_ptr;

};

template <typename VALUE_TYPE>
auto Inv(const CooMatrix<VALUE_TYPE>& A){
  return InvCooMatrix<VALUE_TYPE>(A); }


#endif
