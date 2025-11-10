#ifndef PRECONDITIONER_HPP
#define PRECONDITIONER_HPP

#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include "coomatrix.hpp"
#include "directsolver.hpp"

class CholeskyPrec{

public:
  
  using ValueType       = double;
  using ThisType        = CholeskyPrec;  
  using EigenVectorType = Eigen::Matrix<ValueType,Eigen::Dynamic,1>;
  using EigenMatrixType = Eigen::SparseMatrix<ValueType>;
  using ContainerType   = Eigen::IncompleteCholesky<ValueType>;
  
  CholeskyPrec(const CooMatrix<ValueType>& A):
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

  CholeskyPrec()                               = default;
  CholeskyPrec(const CholeskyPrec&)            = default;
  CholeskyPrec(CholeskyPrec&&)                 = default;
  CholeskyPrec& operator=(const CholeskyPrec&) = default; 
  CholeskyPrec& operator=(CholeskyPrec&&)      = default;     
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
  
  //Data members
  std::size_t                       nr,nc;
  std::shared_ptr<ContainerType> data_ptr;
  std::shared_ptr<EigenVectorType>  u_ptr;
  std::shared_ptr<EigenVectorType>  b_ptr;
  
};





#endif
