#include <cmath>
#include <femtool.hpp>          
#include "partition.hpp"   
#include <vector> 
int main() {
// Load a 2D mesh
Mesh2D Omega;
Read(Omega, "exo2h0.01.mesh");
FeSpace2D Vh= FeSpace2D(Omega);
auto [Sigma, Q] = Partition4(Omega,2);
std::vector<std::size_t> first_col_rows;
for (const auto& [row, col, val] : GetData(Q)) {
    if (col == 0 && val != 0.0) {  
        first_col_rows.push_back(row);
    }
}
std::cout<<"First column rows count: "<<Sigma[0].size()<<"\n";
auto [Sigma0,L]= Restrict(Vh,Sigma[0] ,first_col_rows);



auto F    = [](const R3& x){return std::cos(7.*M_PI*(x[0]+x[1]));}; 
auto V= Vh(F); 
auto V_submesh=L*V;
Plot(Sigma0, V_submesh, "F_on_submesh_0");
//Plot(Sigma0, Sigma0(F), "F_on_submesh_0");

//Plot(Sigma, "domaine_0.01_colored_4");/**/
return 0;
}