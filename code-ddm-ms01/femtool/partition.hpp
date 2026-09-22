#include <cmath>
#include <vector>
#include <array>
#include <unordered_set>
#include <femtool.hpp>
#include "coomatrix.hpp"        
#include "mesh.hpp"   
#include "element.hpp"       
using Mesh2DPart =std::vector<Mesh2D>;
static inline std::array<double,2>
tri_centroid(const Element<2>& T) {
  double cx = (T[0][0] + T[1][0] + T[2][0]) / 3.0;
  double cy = (T[0][1] + T[1][1] + T[2][1]) / 3.0;
  return {cx, cy};
}

using FeSpace2D = FeSpace<2>;
using FeSpace2DxCoo = std::pair<FeSpace2D,CooMatrix<double>>;
std::pair<Mesh2DPart, CooMatrix<double>>
Partition4(const Mesh2D& Omega)
{
    const std::size_t ne = Omega.size();

    Mesh2DPart Sigma(4);
    for (int p = 0; p < 4; ++p)
        Sigma[p] = Mesh2D(Omega.nodes());


    std::vector<std::pair<std::size_t,std::size_t>> pattern;
    pattern.reserve(ne);

    for (std::size_t j = 0; j < ne; ++j) {
        const auto& T = Omega[j];
        auto [cx, cy] = tri_centroid(T);

        int jx = (cx < 0.5) ? 0 : 1;
        int ky = (cy < 0.5) ? 0 : 1;
        int p  = 2 * jx + ky;            // p in {0,1,2,3}

        Sigma[p].push_back(T);
        pattern.emplace_back(j, static_cast<std::size_t>(p));
    }

   
    CooMatrix<double> Q = BooleanMatrix(ne, 4, pattern);

    return {Sigma, Q};
}

void Plot(const std::vector<Mesh2D>& Sigma,
           const std::string& filename)
{
    if (Sigma.empty()) return;

    std::vector<std::string> tag =
        {"Vertices", "Edges", "Triangles", "Tetrahedra"};

    // Open output file
    std::ofstream f;
    std::string meshfile = filename;
    if (meshfile.size() < 5 || meshfile.substr(meshfile.size() - 5) != ".mesh")
        meshfile += ".mesh";
    f.open(meshfile.c_str());


    f << "MeshVersionFormatted 3\n\n";
    f << "Dimension\n3\n\n";

 
    const auto& v = Sigma.front().nodes();
    const auto& v0 = v[0];
    f << "Vertices\n";
    f << v.size() << "\n";
    for (const auto& x : v)
        f << x << "\t1\n";
    f << "\n";
    std::size_t ne_total = 0;
    for (const auto& m : Sigma) ne_total += m.size();
    f << tag[2] << "\n";
    f << ne_total << "\n";
    for (std::size_t p = 0; p < Sigma.size(); ++p) {
        const auto& m = Sigma[p];
        const auto& vp = m.nodes();
        const auto& vp0 = vp[0];

        for (const auto& e : m) {
            for (std::size_t j = 0; j < 2 + 1; ++j)
                f << 1 + int(&e[j] - &vp0) << "\t";
            f << int(p + 1) << "\n"; // region tag = color
        }
    }

    f << "\nEnd";
    f.close();
}


std::pair<Mesh2DPart, CooMatrix<double>>
Partition4(const Mesh2D& Omega, const std::size_t& nl)
{
 
    auto [Sigma_part4, Q_part4] = Partition4(Omega);
    const std::size_t ne = Omega.size();

    // in[i] = triangles déjà dans Σ_i
    std::array<std::unordered_set<std::size_t>, 4> in;
    std::vector<std::pair<std::size_t,std::size_t>> pattern;
    pattern.reserve(ne * 2); 

    // init à partir de Q_part4
    for (const auto& [row, col, val] : GetData(Q_part4)) {
        if (val != 0.0) {
            in[col].insert(row);
            pattern.emplace_back(row, col);
        }
    }

    // Grow nl layers
    for (std::size_t l = 0; l < nl; ++l) {
        for (std::size_t i = 0; i < 4; ++i) {

            auto [mbi, cntbi] = Boundary(Sigma_part4[i]);

            for (std::size_t j = 0; j < ne; ++j) {
                if (in[i].count(j)) continue;  // déjà dedans

                const auto& T = Omega[j];
                auto edgesT = Boundary(T);

                bool touches = false;
                for (std::size_t k = 0; !touches && k < edgesT.size(); ++k) {
                    const Element<1>& eT = edgesT[k];
                    for (std::size_t e = 0; e < cntbi.size(); ++e) {
                        if (eT == mbi[e]) { touches = true; break; }
                    }
                }

                if (touches) {
                    Sigma_part4[i].push_back(T);
                    in[i].insert(j);
                    pattern.emplace_back(j, i);   // appartenance booléenne
                }
            }
        }
    }

    // construire R booléenne finale
    CooMatrix<double> R = BooleanMatrix(ne, 4, pattern);

    return {Sigma_part4, R};
}




std::pair<Mesh2DPart, CooMatrix<double>>
Partition16(const Mesh2D& Omega)
{
    const std::size_t ne = Omega.size();

    // 16 submeshes sharing the same node cloud
    Mesh2DPart Sigma(16);
    for (int p = 0; p < 16; ++p)
        Sigma[p] = Mesh2D(Omega.nodes());

    // Collect the (row=j, col=p) pattern for Q
    std::vector<std::pair<std::size_t,std::size_t>> pattern;
    pattern.reserve(ne);

    for (std::size_t j = 0; j < ne; ++j) {
        const auto& T = Omega[j];
        auto [cx, cy] = tri_centroid(T);
        int jx = (cx < 0.25) ? 0 : (cx < 0.50) ? 1 : (cx < 0.75) ? 2 : 3;
        int ky = (cy < 0.25) ? 0 : (cy < 0.50) ? 1 : (cy < 0.75) ? 2 : 3;


        int p  = 4 * jx + ky; 

        Sigma[p].push_back(T);
        pattern.emplace_back(j, static_cast<std::size_t>(p));
    }


    CooMatrix<double> Q = BooleanMatrix(ne, 16, pattern);

    return {Sigma, Q};
}

std::pair<Mesh2DPart, CooMatrix<double>>
Partition16(const Mesh2D& Omega, const std::size_t& nl)
{
         // base partition (sans overlap)
    auto [Sigma_part16, Q_part16] = Partition16(Omega);
    const std::size_t ne = Omega.size();
    // in[i] = triangles déjà dans Σ_i
    std::array<std::unordered_set<std::size_t>, 16 > in;

    std::vector<std::pair<std::size_t,std::size_t>> pattern;
    pattern.reserve(ne * 2); 

    // init à partir de Q_part16
    for (const auto& [row, col, val] : GetData(Q_part16)) {
        if (val != 0.0) {
            in[col].insert(row);
            pattern.emplace_back(row, col);
        }
    }

    // Grow nl layers
    for (std::size_t l = 0; l < nl; ++l) {
        for (std::size_t i = 0; i < 16; ++i) {

            auto [mbi, cntbi] = Boundary(Sigma_part16[i]);

            for (std::size_t j = 0; j < ne; ++j) {
                if (in[i].count(j)) continue;  // déjà dedans

                const auto& T = Omega[j];
                auto edgesT = Boundary(T);

                bool touches = false;
                for (std::size_t k = 0; !touches && k < edgesT.size(); ++k) {
                    const Element<1>& eT = edgesT[k];
                    for (std::size_t e = 0; e < cntbi.size(); ++e) {
                        if (eT == mbi[e]) { touches = true; break; }
                    }
                }

                if (touches) {
                    Sigma_part16[i].push_back(T);
                    in[i].insert(j);
                    pattern.emplace_back(j, i);   // appartenance booléenne
                }
            }
        }
    }

    // construire R booléenne finale
    CooMatrix<double> R = BooleanMatrix(ne, 16, pattern);

    return {Sigma_part16, R};
 }

FeSpace2DxCoo Restrict(const FeSpace2D& Vh,const Mesh2D& Gamma,const std::vector<std::size_t>& tbl){
    FeSpace2D Vh_Gamma(Gamma);
    std::cout<<"Dim Vh: "<<dim(Vh)<<" Dim Vh_Gamma: "<<dim(Vh_Gamma)<<"\n";
    std::cout<<"Restriction table size: "<<tbl.size()<<"\n";
    CooMatrix<double> R(dim(Vh_Gamma),dim(Vh));

    for (std::size_t j=0; j<Gamma.size(); ++j){
        for (size_t k = 0; k < local_dim(Vh_Gamma); ++k)
        {
            R.push_back(Vh_Gamma[j][k], Vh[tbl[j]][k], 1.0);
        }
    }

    return {Vh_Gamma,R};
}

std::vector<FeSpace2DxCoo>
Partition4(const FeSpace2D& Vh, const std::size_t& nl)
{
    auto [Sigma_part4, R_part4] = Partition4(Vh.mesh(), nl);

    std::vector<FeSpace2DxCoo> Vh_part4(4);

    for (std::size_t p = 0; p < 4; ++p) {
        std::vector<std::size_t> tbl;
        tbl.reserve(Sigma_part4[p].size());

        // collect Omega triangle indices belonging to subdomain p
        for (const auto& [j, k, v] : GetData(R_part4)) {
            if (k == p && v != 0.0) tbl.push_back(j);
        }

        Vh_part4[p] = Restrict(Vh, Sigma_part4[p], tbl);
        
    }

    return Vh_part4;
}


std::vector<FeSpace2DxCoo>
Partition16(const FeSpace2D& Vh, const std::size_t& nl){
     auto [Sigma_part16, R_part16] = Partition16(Vh.mesh(), nl);

    std::vector<FeSpace2DxCoo> Vh_part16(16);

    for (std::size_t p = 0; p < 16; ++p) {
        std::vector<std::size_t> tbl;
        tbl.reserve(Sigma_part16[p].size());
        for (const auto& [j, k, v] : GetData(R_part16)) {
            if (k == p && v != 0.0) tbl.push_back(j);
        }

        Vh_part16[p] = Restrict(Vh, Sigma_part16[p], tbl);
        
    }

    return Vh_part16;

}