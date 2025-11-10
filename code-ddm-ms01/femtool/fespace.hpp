#ifndef FESPACE_HPP
#define FESPACE_HPP

#include <array>
#include <iostream>
#include <algorithm>
#include <assert.h>
#include <filesystem>
#include <type_traits>
#include "element.hpp"

//###########################//
//   Cellule element fini    //
//###########################//

template <std::size_t DIM>
class FeCell{

public:

  static constexpr std::size_t geometry_dim = DIM;
  static constexpr std::size_t space_dim    = DIM+1;
  static constexpr std::size_t value_dim    = 1;
  
  using EltType  = Element<DIM>;
  using DataType = std::array<std::size_t,DIM+1> ;
  
  FeCell(): e_ptr(nullptr) {
    for(std::size_t j=0; j<space_dim; ++j){
      data[j]=j;
    }
  }

  FeCell(const EltType& e): e_ptr(&e) {
    for(std::size_t j=0; j<space_dim; ++j){
      data[j]=j;
    }
  }
  
  FeCell(const FeCell&)            = default;
  FeCell(FeCell&&)                 = default;  
  FeCell& operator=(const FeCell&) = default;
  FeCell& operator=(FeCell&&)      = default;

  std::size_t& operator[](const std::size_t& j){
    return data[j];}

  const std::size_t& operator[](const std::size_t& j) const {
    return data[j];}

  void attach(const EltType& e){e_ptr=&e;}

  bool attached() const {return (e_ptr!=nullptr);}
  
  const EltType& elt() const {return *e_ptr;}  
  
  typedef typename DataType::const_iterator const_iterator;
  typedef typename DataType::iterator             iterator;
  iterator       begin ()       {return data.begin(); }
  iterator       end   ()       {return data.end();   }
  const_iterator begin () const {return data.cbegin();}
  const_iterator end   () const {return data.cend();  }  
  const_iterator cbegin() const {return data.cbegin();}
  const_iterator cend  () const {return data.cend();  }  

  friend std::ostream&
  operator<<(std::ostream& o, const FeCell& c){
    for(const auto& j:c){o << j << "\t";}
    return o;}

  friend constexpr std::size_t
  dim(const FeCell&){return space_dim;}
  
  std::size_t size() const {return data.size();}
  
  friend auto Boundary(const FeCell& c){
    std::array<FeCell<DIM-1>,DIM+1> f;
    std::size_t kk=0;
    for(std::size_t j=0; j<DIM+1; ++j){
      kk=0;
      for(std::size_t k=0; k<DIM+1; ++k){
	if(k!=DIM-j){f[j][kk++] = c[k];}
      }
    }
    return f;
  }
  
  friend bool
  operator==(const FeCell& c1, const FeCell& c2){
    return (c1.e_ptr==c2.e_ptr) && (c1.data==c2.data);}
  
  friend auto
  Points(const FeCell& c){
    assert(c.attached()); 
    const auto& e = c.elt();
    std::array<R3,space_dim> x;
    for(std::size_t j=0; j<c.size(); ++j){
      x[j] = e[j];} 
    return x;
  }

  template <typename FctType> 
  auto operator()(FctType fct) const {    
    std::array<std::invoke_result_t<FctType,R3>,space_dim> v;
    auto x = Points(*this);
    for(std::size_t j=0; j<space_dim; ++j){
      v[j]=fct(x[j]);}
    return v;
  }
  
private:
  
  const EltType*               e_ptr;
  std::array<std::size_t,DIM+1> data;  
  
};

using FeCell1D = FeCell<1>;
using FeCell2D = FeCell<2>;
using FeCell3D = FeCell<3>;


//###########################//
//   Espace element fini     //
//###########################//

template <std::size_t DIM>
class FeSpace{

public:
  
  using CellType      = FeCell <DIM>;
  using EltType       = Element<DIM>;
  using MeshType      = Mesh   <DIM>;
  using DataContainer = std::vector<CellType>;

  static constexpr std::size_t geometry_dim    = DIM;
  static constexpr std::size_t local_space_dim = CellType::space_dim;
  static constexpr std::size_t value_dim       = CellType::value_dim;
  
  FeSpace(const MeshType& mesh0 = MeshType()):
    data_ptr(std::make_shared<DataContainer>(mesh0.size())),
    mesh_(mesh0), space_dim(0)
  {
    attach(mesh_);
    
    std::size_t d = local_space_dim;
    std::vector<std::tuple<const R3*, std::size_t>> tbl;
    tbl.reserve(d*size());

    for(std::size_t j=0; j<size(); ++j){
      for(std::size_t k=0; k<d; ++k){
	tbl.push_back({ &(mesh_[j][k]), k+j*d });
      }
    }
    std::sort(tbl.begin(),tbl.end());

    if(size()>0){
      auto [p0,jk0] = tbl[0];
      for(const auto& [p,jk]:tbl){
	if(p>p0){
	  p0=p;
	  ++space_dim;
	}
	(*data_ptr)[jk/d][jk%d] = space_dim;
      }
      ++space_dim;
    }
  }
  
  FeSpace(const FeSpace&)            = default;
  FeSpace(FeSpace&&)                 = default;
  FeSpace& operator=(const FeSpace&) = default;
  FeSpace& operator=(FeSpace&&)      = default;

  typedef typename DataContainer::const_iterator const_iterator;
  typedef typename DataContainer::iterator             iterator;
  iterator       begin()        {return data_ptr->begin(); }
  iterator       end()          {return data_ptr->end();   }
  const_iterator begin()  const {return data_ptr->cbegin();}
  const_iterator end()    const {return data_ptr->cend();  }  
  const_iterator cbegin() const {return data_ptr->cbegin();}
  const_iterator cend()   const {return data_ptr->cend();  }  
  
  void attach([[maybe_unused]] const MeshType& m){
    auto& data = *data_ptr;
    assert(data.size()==m.size());
    for(std::size_t j=0; j<size(); ++j){
      data[j].attach(mesh_[j]);}
  }
  
  std::size_t size() const {return data_ptr->size();}

  friend std::size_t dim(const FeSpace& Vh){return Vh.space_dim;}

  friend constexpr std::size_t
  local_dim(const FeSpace&){return FeSpace::local_space_dim;}
  
  auto mesh()  const {return mesh_;}
  
  auto nodes() const {return mesh().nodes();}
  
  const auto& operator[](const std::size_t& j) const {
    return (*data_ptr)[j];}

  friend std::ostream&
  operator<<(std::ostream& o, const FeSpace& Vh){
    for(const auto& I:Vh){o << I << "\n";}
    return o;
  }

  friend auto
  Points(const FeSpace& Vh){
    std::vector<R3> x(dim(Vh));
    for(std::size_t j=0; j<Vh.size(); ++j){
      auto xj = Points(Vh[j]);
      for(std::size_t k=0; k<xj.size(); ++k){
	x[Vh[j][k]]=xj[k];
      }
    }
    return x;
  }
    
  template <typename FctType> 
  auto operator()(FctType fct) const {    
    std::vector<std::invoke_result_t<FctType,R3>> v(space_dim);
    auto x = Points(*this);
    for(std::size_t j=0; j<space_dim; ++j){
      v[j]=fct(x[j]);}
    return v;
  }
    
private:

  std::shared_ptr<DataContainer> data_ptr;
  MeshType                          mesh_;
  std::size_t                   space_dim;
  
};


template <std::size_t DIM, typename VALUE_TYPE>
void Plot(const FeSpace<DIM>&             Vh,
	  const std::vector<VALUE_TYPE>&   u,
	  std::filesystem::path     filename){
  
  assert( dim(Vh)==u.size() );
  
  //###############//
  //  Fichier mesh //  
  std::vector<std::string> tag =
    {"Vertices", "Edges", "Triangles", "Tetrahedra"};

  filename.replace_extension(".mesh");
  std::ofstream f;
  f.open(filename.c_str());
  
  // Preambule
  f << "MeshVersionFormatted 3\n\n";
  f << "Dimension\n3\n\n";

  // Section noeuds
  auto x = Points(Vh);
  f << "Vertices\n";
  f << x.size() << "\n";
  for(const auto& xj:x){
    f << xj << "\t1\n";}
  f << "\n";

  // Section elements
  f << tag[DIM]   << "\n";
  f << Vh.size() << "\n";
  for(const auto& I:Vh){
    for(const auto& Ik:I){
      f << 1+Ik << "\t";}
    f << "1\n";}
  f << "\n";

  // Fermeture
  f << "\nEnd";
  f.close();
  
  //##############//
  //  Fichier sol //

  // Ouverture
  filename.replace_extension(".sol");
  f.open(filename.c_str());

  // Preambule
  f << "MeshVersionFormatted 3\n\n";
  f << "Dimension\n3\n\n";

  // Donnees
  f << "SolAtVertices\n";
  f << dim(Vh) << "\n";
  f << "1\t1\n";
  for(std::size_t j=0; j<u.size(); ++j){
    f << u[j] << "\n";}

  // Fermeture
  f << "\nEnd";
  f.close();
  
}


template <std::size_t DIM>
auto Boundary(const FeSpace<DIM>& Vh,
	      const std::tuple<Mesh<DIM-1>,std::vector<std::size_t>>& Gamma_x_tbl){

  const auto& [Gamma,tbl] = Gamma_x_tbl;
  auto Wh = FeSpace(Gamma);
  using NxN = std::pair<std::size_t,std::size_t>;
  std::vector<NxN> Wh_x_Vh(dim(Wh));

  for(std::size_t j=0; j<tbl.size(); ++j){
    std::size_t k = tbl[j]/(DIM+1);
    std::size_t l = tbl[j]%(DIM+1);    

    const auto& I = Wh[j];
    auto II = Boundary(Vh[k])[l];
    for(std::size_t p=0; p<I.size(); ++p){
      Wh_x_Vh[I[p]] = std::make_pair(I[p],II[p]);
    }    
  }  

  CooMatrix<double> B(dim(Wh),dim(Vh));
  for(const auto& [j,k]:Wh_x_Vh){
    B.push_back(j,k,1.);}

  return std::make_pair(Wh,B);  
}

template <std::size_t DIM>
auto Boundary(const FeSpace<DIM>& Vh){
  return Boundary(Vh, Boundary(Vh.mesh()));}

#endif
