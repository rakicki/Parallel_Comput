#include <cmath>
#include <femtool.hpp>         
#include "partition.hpp"     
int main() {
// Load a 2D mesh
Mesh2D Omega;
Read(Omega, "exo2h0.01.mesh");
FeSpace2D Vh= FeSpace2D(Omega);
auto [Sigma, Q] = Partition16(Omega,2);
std::vector<std::size_t> first_col_rows;

for (const auto& [row, col, val] : GetData(Q)) {
    if (col == 0 && val != 0.0) {   // first column only, nonzero entries
        first_col_rows.push_back(row);
    }
}
std::cout<<"First column rows count: "<<Sigma[0].size()<<"\n";
auto [Sigma0,Uh0]= Restrict(Vh,Sigma[0] ,first_col_rows);
auto F    = [](const R3& x){return std::cos(7.*M_PI*(x[0]+x[1]));};
Plot(Sigma0, Sigma0(F), "F_on_submesh_0");
// Write all 16 submeshes into one colored .mesh
Plot(Sigma, "domaine_0.01_colored_16");/**/
return 0;
}  