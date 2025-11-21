
#include <cmath>
#include <femtool.hpp>

int main(){

  // Instantiation of a 2D domain
  Mesh2D Omega;

  // Loading a 2D mesh
  Read(Omega,"exo2h0.05.mesh");

  // Assembly of a finite element space over Omega
  auto Vh   = FeSpace(Omega);

 
  
  // Function x = (x1,x2) -> cos(omega*x1)
  auto F    = [](const R3& x){return std::cos(10.*M_PI*x[0]);};

  // Manufactured "exact solution" obtained
  // by nodal evaluation of f at the degrees of freedom of Vh
  auto ue   = Vh(F);

  // Assembly of the mass matrix for Vh

  // Assembly of finite element matrix of the operator
  // -Delta + 1 with Neumann BC
  auto A    = Stiffness(Vh)+Mass(Vh);

  // Manufactured "right hand side"
  auto b    = A(ue);
  
  // Nmerical solution obtained by conjugate gradient
  auto uh   = cgsolve(A,b);
 

  // Error between discrete solution and
  // (manufactured) exact solution

  // Plotting exact solution with vizir4
  //Plot(Vh,ue,"ue");

  // Plotting error with vizir4
  //Plot(Vh,err,"err");
  return 0;


}