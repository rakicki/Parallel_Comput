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



/*
class AdditiveSchwarz
{
public:
    using ItemType = std::pair<InvCooMatrix<double>, CooMatrix<double>>;

    AdditiveSchwarz(const CooMatrix<double>& A,
                    const std::vector<CooMatrix<double>>& Rlist)
    {
        InvAxR.reserve(Rlist.size());

        for (const CooMatrix<double>& Rj : Rlist)
        {
            CooMatrix<double> Aj = (Rj.T()*A) *Rj;

            // important for InvCooMatrix: ensure merged duplicates
            Aj.sort();
            // Factorize Aj
            InvCooMatrix<double> InvAj(Aj);

            // Store (InvAj, Rj)
            InvAxR.emplace_back(std::move(InvAj), Rj);
        }
    }

    // ASM application
    std::vector<double> operator*(const std::vector<double>& u) const
    {
        std::vector<double> y(u.size(), 0.0);

        for (const auto& [InvAj, Rj] : InvAxR)
        {
            // r = Rj^T u
            std::vector<double> r = Rj.T()(u);

            // z = Aj^{-1} r
            std::vector<double> z = InvAj(r);

            // y += Rj z
            auto Rjz = Rj(z);
            for (std::size_t i = 0; i < y.size(); ++i) y[i] += Rjz[i];
        }

        return y;
    }

private:
    std::vector<ItemType> InvAxR;
};
*/
class AdditiveSchwarz
{
public:
    using ItemType = std::pair<InvCooMatrix<double>, CooMatrix<double>>;
    AdditiveSchwarz(const CooMatrix<double>& A,
                    const std::vector<CooMatrix<double>>& Rlist)
    {
        InvAxR.reserve(Rlist.size());
        for (const auto& Rj : Rlist)
        {
            CooMatrix<double> Aj = (Rj.T() * A) * Rj;
            Aj.sort();
            InvCooMatrix<double> InvAj(Aj);
            InvAxR.emplace_back(std::move(InvAj), Rj);
        }
    }

    std::vector<double> operator*(const std::vector<double>& u) const
    {
        std::vector<double> y(u.size(), 0.0);
        for (const auto& [InvAj, Rj] : InvAxR)
        {
            std::vector<double> r = Rj.T(u);
            std::vector<double> z = InvAj(r);
            // y += Rj z
            std::vector<double> add = Rj(z);
            for (std::size_t i = 0; i < y.size(); ++i)
                y[i] += add[i];
        }

        return y;
    }
private:
    std::vector<ItemType> InvAxR;
};

#endif
