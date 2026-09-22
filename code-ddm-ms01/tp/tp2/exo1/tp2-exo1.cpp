
#include <cmath>


#include <femtool.hpp>          // meshes/spaces/assembly/types (defines cplx etc.)
#include "preconditioner.hpp"   // preconditioner
#include "iterativesolver.hpp"  // now “pure” linear algebra, no FEM includes

int main(){

  // Instantiation of a 2D domain
  Mesh2D Omega;

  // Loading a 2D mesh
  Read(Omega,"tp1-1.mesh");

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

 
  auto b    = A(ue);
  CholeskyPrec precon(A);
  // Nmerical solution obtained by conjugate gradient
  auto uh   = PCGSolver(A,b,precon);

 

  // Error between discrete solution and
  // (manufactured) exact solution

  // Plotting exact solution with vizir4
  //Plot(Vh,ue,"ue");

  // Plotting error with vizir4
  //Plot(Vh,err,"err");
  return 0;


}
