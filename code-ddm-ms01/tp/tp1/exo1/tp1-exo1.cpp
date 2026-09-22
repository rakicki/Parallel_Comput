#include <cmath>
#include <femtool.hpp>

int main(){
    // Instantiation of a 2D domain
  Mesh2D Omega;

  // Loading a 2D mesh
  Read(Omega,"tp1-1.mesh"); 
    // Assembly of a finite element space over Omega
  auto Vh   = FeSpace(Omega);
    std::vector<double> uex = Vh(
    [](const R3& x){
      return std::cos(10*M_PI*(x[0]+x[1]));
    }
  );
    // Plotting the exact solution with vizir4
  Plot(Vh,uex,"tp1-1-output");

}