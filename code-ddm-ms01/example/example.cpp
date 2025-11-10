#include <cmath>
#include <femtool.hpp>

int main(){

  //########################################################//
  //                   EXAMPLE 1                            //
  //########################################################//

  // Instantiation of a 3D domain
  Mesh3D Omega;

  // Loading a 3D mesh
  Read(Omega,"example1.mesh");

  // Assembly of a finite element space over Omega
  auto Vh   = FeSpace(Omega);

  // Frequency parameter
  double k = 20.;
  
  // Function x = (x1,x2) -> cos(omega*x1)
  auto F    = [&k](const R3& x){return std::cos(k*x[0]);};

  // Evaluating nodal values of f at the degrees of freedom of Vh
  auto u    = Vh(F);

  // Plotting F with vizir4
  Plot(Vh,u,"output");
  
  // Assembly of the finite element matrix of the
  // boundary value problem
  // -Delta u + u = rhs on Omega
  // \partial_n u = 0   on the boundary  
  auto A    = Stiffness(Vh)+Mass(Vh);

  // Generation of a random vector of nodal values on Vh
  auto x    = RandomVec(dim(Vh));

  // Computing a corresponding  right-hand side
  auto b    = A(x);

  // Lazy invertion of A (sparse LU factorization) 
  auto InvA = Inv(A);  

  // Computing A^{-1}b
  auto y    = InvA(b);

  // Computing error: x and y should be very close
  std::cout << "Norm(x-y)/Norm(x) = ";
  std::cout <<  Norm(x-y)/Norm(x) << "\n";




  
  //########################################################//
  //                   EXAMPLE 2                            //
  //########################################################//

  
  // // Instantiation of a 2D domain
  // Mesh2D Omega;

  // // Loading a 2D mesh
  // Read(Omega,"example2.mesh");

  // // Assembly of a finite element space over Omega
  // auto Vh   = FeSpace(Omega);

  // // Frequency parameter
  // double k = 20.;
  
  // // Function x = (x1,x2) -> cos(omega*x1)
  // auto F    = [](const R3& x){return std::cos(k*x[0]);};

  // // Manufactured "exact solution" obtained
  // // by nodal evaluation of f at the degrees of freedom of Vh
  // auto ue   = Vh(F);

  // // Assembly of the mass matrix for Vh
  // auto M    = Mass(Vh);

  // // Assembly of finite element matrix of the operator
  // // -Delta + 1 with Neumann BC
  // auto A    = Stiffness(Vh)+Mass(Vh);

  // // Manufactured "right hand side"
  // auto b    = A(ue);

  // // Nmerical solution obtained by conjugate gradient
  // auto uh   = cgsolve(A,b);  

  // // Error between discrete solution and
  // // (manufactured) exact solution
  // auto err  = ue-uh;

  // // Displaying L2 norm of the error
  // std::cout << "|ue-uh|^2/|ue|^2 = ";
  // std::cout <<  (M(err)|err)/(M(ue)|ue) << "\n";  

  // // Plotting exact solution with vizir4 
  // Plot(Vh,ue,"ue");

  // // Plotting error with vizir4
  // Plot(Vh,err,"err");

}
