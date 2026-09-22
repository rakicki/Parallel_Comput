#include <cmath>
#include <femtool.hpp>          // meshes/spaces/assembly/types (defines cplx etc.)
#include "partition.hpp"       // mesh partitioning utilities
#include "preconditioner.hpp"   // preconditioner
#include "iterativesolver.hpp" 
int main() {
// Load a 2D mesh
Mesh2D Omega;
Read(Omega, "exo2h0.05.mesh");
FeSpace2D Vh= FeSpace2D(Omega);
auto F    = [](const R3& x){return std::cos(10.*M_PI*x[0]);};
auto ue   = Vh(F);
std::vector<FeSpace2DxCoo> Uh_and_Rtlist = Partition4(Vh,2);

 auto Uh = std::vector<FeSpace2D>(4);
    std::vector<CooMatrix<double>> Rlist(4);

    for (std::size_t p = 0; p < 4; ++p) {
        Uh[p]    = Uh_and_Rtlist[p].first;
        auto Rt  = Uh_and_Rtlist[p].second; // Rt is Rj^T  (local x global)
        Rlist[p] = Rt.T();                  // Rj = (Rj^T)^T (global x locals
    }
auto A    = Stiffness(Vh)+Mass(Vh);

auto b    = A(ue);
std::cout << "Assembled system of size  P" << NbRow(Rlist[0]) << " x " << NbCol(Rlist[0]) << std::endl;
std::cout << "Assembled system of size  A" << NbRow(A) << " x " << NbCol(A) << std::endl;
std::cout << "dim(Vh)=" << dim(Vh)
          << " dim(Uh)=" << dim(Uh[0]) << "\n";
std::cout << "Rj: " << NbRow(Rlist[0]) << " x " << NbCol(Rlist[0]) << "\n";
AdditiveSchwarz precon(A,Rlist);
  // Nmerical solution obtained by conjugate gradien 
  auto uh   = PCG_ASM_Solver(A,b,precon);
// Partition into 4 parts (like your Partition4)

Plot(Vh, uh, "uh_full_domain");
return 0;
}

