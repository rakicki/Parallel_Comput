#ifndef MESH_HPP
#define MESH_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>
#include <vector>
#include <filesystem>

template <std::size_t DIM>
class Mesh{
  
private:

  using EltType = Element<DIM>;
  using DataContainer = std::vector<EltType>;

  std::shared_ptr<DataContainer> data_ptr;
  Nodes nodes_;
  
public:
  
  Mesh(const Nodes& nodes0 = Nodes()):
    data_ptr(std::make_shared<DataContainer>()),
    nodes_(nodes0) {};
  
  Mesh(const Mesh&)            = default;
  Mesh(Mesh&&)                 = default;
  Mesh& operator=(const Mesh&) = default;
  Mesh& operator=(Mesh&&)      = default;

  typedef typename DataContainer::const_iterator const_iterator;
  typedef typename DataContainer::iterator             iterator;
  iterator       begin()        {return data_ptr->begin(); }
  iterator       end()          {return data_ptr->end();   }
  const_iterator begin()  const {return data_ptr->cbegin();}
  const_iterator end()    const {return data_ptr->cend();  }  
  const_iterator cbegin() const {return data_ptr->cbegin();}
  const_iterator cend()   const {return data_ptr->cend();  }  

  const auto& operator[](const std::size_t& j) const {
    return (*data_ptr)[j];}

  void push_back(const EltType& e){
    data_ptr->push_back(e);}

  void push_back(const R3& x0)
    requires(DIM==0)
  {push_back(EltType(x0));}

  void push_back(const R3& x0, const R3& x1)
    requires(DIM==1)
  {push_back(EltType(x0,x1));}

  void push_back(const R3& x0, const R3& x1, const R3& x2)
    requires(DIM==2)
  {push_back(EltType(x0,x1,x2));}

  void push_back(const R3& x0, const R3& x1, const R3& x2, const R3& x3)
    requires(DIM==3)
  {push_back(EltType(x0,x1,x2,x3));}
  
  void reserve(const std::size_t& n){
    data_ptr->reserve(n);}  
  
  std::size_t size() const {
    return data_ptr->size();}

  Nodes nodes() const {
    return nodes_;}
  
  friend bool operator==(const Mesh& m1,const Mesh& m2){
    const auto& data1 = *(m1.data_ptr);
    const auto& data2 = *(m2.data_ptr);    
    return (m1.data_ptr==m2.data_ptr) ||
      ( (m1.nodes()==m1.nodes()) && (data1==data2) );
  }
  
};

using Mesh1D = Mesh<1>;
using Mesh2D = Mesh<2>;
using Mesh3D = Mesh<3>;



template <std::size_t DIM>
auto Boundary(const Mesh<DIM>& m){

  using EltxN = std::tuple<Element<DIM-1>,std::size_t>;  
  std::vector<EltxN> tbl;
  tbl.reserve(2*(DIM+1)*m.size());
  
  for(std::size_t j=0; j<m.size(); ++j){
    auto f = Boundary(m[j]);    

    for(std::size_t k=0; k<f.size(); ++k){
      tbl.push_back({f[k],k+j*(DIM+1)});
    }
  }

  std::sort(tbl.begin(),tbl.end());
  
  Mesh<DIM-1> mb(m.nodes());
  std::vector<std::size_t> cntb;
  assert(tbl.size()>0);
  for(std::size_t j=0; j<tbl.size(); ++j){

    const auto& [f1,n1] = tbl[j];
    if( (j+1)==(tbl.size()) ){
      mb.push_back(f1);
      cntb.push_back(n1);
      break;
    }
    
    const auto& [f2,n2] = tbl[j+1];    
    if(f1==f2){++j;}
    else{
      mb.push_back(f1);
      cntb.push_back(n1);
    }    
  }
  
  return  std::make_tuple(mb,cntb);
}


template <std::size_t DIM>
void Read(Mesh<DIM>& m,
	  std::filesystem::path filename){
  
  std::vector<std::string> tag =
    {"Vertices", "Edges", "Triangles", "Tetrahedra"};
  
  // Ouverture fichier

  filename.replace_extension(".mesh");
  std::ifstream f;
  f.open(filename.c_str());
  std::string line;
  std::istringstream iss;

  // Lecture section noeuds
  auto v = m.nodes();
  if(v.size()==0){
    Read(v,filename);
  }
  
  // Lecture section elements  
  std::size_t ne = 0;  
  while(std::getline(f,line)){
    auto pos = line.find(tag[DIM]);
    if(pos!=std::string::npos){break;}          
  }
  std::getline(f,line);
  iss.str(line);
  iss >> ne;    
  iss.clear();

  Element<DIM>  e;
  std::size_t I = 0;
  for(std::size_t j=0; j<ne; ++j){
    std::getline(f,line);
    iss.str(line);    
    for(std::size_t k=0; k<DIM+1; ++k){
      iss >> I;
      I = I-1;
      e.push_back(v[I]);
    }
    e.sort();
    m.push_back(e);

    e.clear();
    iss.clear();
  }  

  // Fermeture fichier
  f.close();
}




template <std::size_t DIM>
void Write(const Mesh<DIM>&      m,
	   std::filesystem::path filename){

  std::vector<std::string> tag =
    {"Vertices", "Edges", "Triangles", "Tetrahedra"};

  // Ouverture fichier
  filename.replace_extension(".mesh");
  std::ofstream f;
  f.open(filename.c_str());

  // Preambule
  f << "MeshVersionFormatted 3\n\n";
  f << "Dimension\n3\n\n";

  // Section noeuds
  auto v = m.nodes();
  const auto& v0 = v[0];
  f << "Vertices\n";
  f << v.size() << "\n";
  for(const auto& x:v){
    f << x << "\t1\n";}
  f << "\n";

  // Section elements
  f << tag[DIM]   << "\n";
  f << m.size() << "\n";
  for(const auto& e:m){
    for(std::size_t j=0; j<DIM+1; ++j){
      f << 1+int(&e[j]-&v0) << "\t";}
    f << "1\n";}

  // Fermeture
  f << "\nEnd";
  f.close();
  
}



template <std::size_t DIM, typename VALUE_TYPE>
void Plot(const Mesh<DIM>&             mesh,
	  const std::vector<VALUE_TYPE>&  u,
	  std::filesystem::path    filename){

  assert( mesh.size()==u.size() );

  std::vector<std::string> tag =
    {"Vertices", "Edges", "Triangles", "Tetrahedra"};
  
  //###############//
  //  Fichier mesh //  
  Write(mesh,filename);

  //###############//
  //  Fichier sol  //

  // Ouverture
  filename.replace_extension(".sol");
  std::ofstream f;
  f.open(filename.c_str());

  // Preambule
  f << "MeshVersionFormatted 3\n\n";
  f << "Dimension\n3\n\n";
  
  // Donnees
  f << "SolAt";
  f << tag[DIM] << "\n";
  f << mesh.size() << "\n";
  f << "1\t1\n";
  for(std::size_t j=0; j<u.size(); ++j){
    f << u[j] << "\n";}

  // Fermeture
  f << "\nEnd";
  f.close();
  
}



  
#endif
