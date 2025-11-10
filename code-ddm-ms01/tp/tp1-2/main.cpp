#include <cmath>
#include <femtool.hpp>

int main(){
    // Instantiation of a 3D domain
  Mesh3D Omega;

  // Loading a 2D mesh
  Read(Omega,"tp1-2.mesh"); 
    // Assembly of a finite element space over Omega
  auto Vh   = FeSpace(Omega);
  FeSpace Wh=Boundary(Vh).first;  
  CooMatrix<double> B=Boundary(Vh).second;       

    std::vector<double> uex = Vh(
    [](const R3& x){
      return std::cos(5*M_PI*(x[0]+x[1]+x[2]));
    }
  );
    // Plotting the exact solution with vizir4
  Plot(Wh,B(uex),"tp1-2-output");
return 0;
}