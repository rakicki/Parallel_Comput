#ifndef NODES_HPP
#define NODES_HPP

#include <vector>
#include <memory>
#include <fstream>
#include <sstream>
#include "smallvector.hpp"

class Nodes{

private:

  using DataContainer = std::vector<R3>;  
  std::shared_ptr<DataContainer> data_ptr;
  
public:
  
  Nodes(): data_ptr(std::make_shared<DataContainer>()) {};

  Nodes(const Nodes&)  = default;

  Nodes(Nodes&&)       = default;
  
  Nodes(const std::initializer_list<R3>& x):
    data_ptr(std::make_shared<DataContainer>()){
    for(const auto& xj:x){data_ptr->push_back(xj);}}
  
  Nodes& operator=(const Nodes&) = default;

  Nodes& operator=(Nodes&&)      = default;
  
  void push_back(const R3& x){
    data_ptr->push_back(x);}

  void reserve(const std::size_t& n){
    data_ptr->reserve(n);}

  std::size_t size() const {
    return data_ptr->size();}

  const auto& data() const {return *data_ptr;}
  
  typedef typename DataContainer::const_iterator const_iterator;
  typedef typename DataContainer::iterator             iterator;
  iterator       begin()        {return data_ptr->begin(); }
  iterator       end()          {return data_ptr->end();   }
  const_iterator begin()  const {return data_ptr->cbegin();}
  const_iterator end()    const {return data_ptr->cend();  }  
  const_iterator cbegin() const {return data_ptr->cbegin();}
  const_iterator cend()   const {return data_ptr->cend();  }  
  
  friend std::ostream&
  operator<<(std::ostream& o, const Nodes& vtx){
    for(const auto& x:vtx){ o << x << "\n"; }
    return o;}
  
  const auto& operator[](const std::size_t& j) const {
    return (*data_ptr)[j];}

  friend bool operator==(const Nodes& v1, const Nodes& v2){
    return v1.data_ptr==v2.data_ptr;}
  
};


void Read(Nodes& v, const std::string& filename){

  // Ouverture fichier
  std::ifstream f;
  f.open(filename);
  std::string line;
  std::istringstream iss;

  // Lecture section noeuds
  std::size_t nv;
  while(std::getline(f,line)){
    auto pos = line.find("Vertices");
    if(pos!=std::string::npos){break;}          
  }

  std::getline(f,line);
  iss.str(line);
  iss >> nv;    
  iss.clear();

  R3 x;       
  for(std::size_t j=0; j<nv; ++j){
    std::getline(f,line);
    iss.str(line);
    iss >> x;
    v.push_back(x);    
    iss.clear();
  }

  // Fermeture fichier
  f.close();
}

bool Close(const Nodes& v1,
	   const Nodes& v2){

  bool test = (v1.size()==v2.size());
  if(test){
    for(std::size_t j=0; j<v1.size(); ++j){
      test = test && Close(v1[j],v2[j]);
    }
  }
  return test;
}



#endif
