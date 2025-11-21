#include <cmath>
#include <vector>
#include <array>
#include <unordered_set>
#include <femtool.hpp>
#include "coomatrix.hpp"          // sparse matrix in COO format
#include "mesh.hpp"   
#include "element.hpp"          // 2D mesh data structure
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

    // 4 submeshes sharing the same node cloud
    Mesh2DPart Sigma(4);
    for (int p = 0; p < 4; ++p)
        Sigma[p] = Mesh2D(Omega.nodes());

    // Collect the (row=j, col=p) pattern for Q
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

    // Build Q as a boolean (0/1) sparse matrix from the pattern
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

    // --- Preamble ---
    f << "MeshVersionFormatted 3\n\n";
    f << "Dimension\n3\n\n";

    // --- Vertices (shared across all submeshes) ---
    const auto& v = Sigma.front().nodes();
    const auto& v0 = v[0];
    f << "Vertices\n";
    f << v.size() << "\n";
    for (const auto& x : v)
        f << x << "\t1\n";
    f << "\n";

    // --- Count total number of elements ---
    std::size_t ne_total = 0;
    for (const auto& m : Sigma) ne_total += m.size();

    // --- Triangles (all parts, region = p+1) ---
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

    // --- End ---
    f << "\nEnd";
    f.close();
}





std::pair<Mesh2DPart, CooMatrix<double>>
Partition4(const Mesh2D& Omega, const std::size_t& nl)
{
    // Compute base 4-way once
     auto [Sigma_part4, Q_part4] = Partition4(Omega);
      std::array<std::unordered_set<std::size_t>, 4> in;
 // Fast membership: in[i] holds all triangle indices currently in Σ_i
/* std::vector<std::pair<std::size_t,std::size_t>> pattern;
    pattern.reserve(Omega.size());
  */
    for (const auto& [row, col, val] : GetData(Q_part4)) {
        if (val != 0.0) {in[col].insert(row);
           // pattern.emplace_back(row, col);
        }
    }
        
    // Start R as a copy of Q (same nonzeros to begin with)
    CooMatrix<double> R =Q_part4;//= BooleanMatrix(NbRow(Q_part4), NbCol(Q_part4), pattern);

   

    const std::size_t ne = Omega.size();

    // Grow nl layers
    for (std::size_t l = 0; l < nl; ++l) {
        for (std::size_t i = 0; i < 4; ++i) {
            auto [mbi, cntbi] = Boundary(Sigma_part4[i]); // mbi: Mesh<1> of boundary edges

            // Try to add triangles not yet in Σ_i that touch its boundary
            for (std::size_t j = 0; j < ne; ++j) {
                if (in[i].count(j)) continue; // already in Σ_i
                
                const auto& T = Omega[j];
                auto edgesT   = Boundary(T);  // 3 edges (Element<1>)
                bool touches  = false;
                auto tri = 0;
                // edge of T vs boundary edges of Σ_i
                for (std::size_t k = 0; !touches && k < edgesT.size(); ++k) {
                    const Element<1>& eT = edgesT[k];
                    for (std::size_t e = 0; e < cntbi.size(); ++e) {
                        
                        
                        if (eT == mbi[e]) { touches = true; tri = cntbi[e] / 3; break; } // edge-to-edge compare
                    }
                }
                if (touches) {
                    Sigma_part4[i].push_back(T);
                    // put new mmbership in R
                    R.push_back(j, i, tri);   // append new membership
                    in[i].insert(j);          // mark as in Σ_i
                }
            }
        }
    }
    return {Sigma_part4, R};
}


/*
std::pair< Mesh2DPart,CooMatrix<double> >
Partition4(const Mesh2D& Omega, const std::size_t& nl){
    Mesh2DPart Sigma_part4=Partition4(Omega).first;
    CooMatrix<double> Q_part4=Partition4(Omega).second;
    CooMatrix<double> R;
   R=Q_part4;
    for (size_t l = 0; l < nl; l++)
    {
    for (int i = 0; i < 4; i++)
    {   
        auto [mbi, cntbi] = Boundary(Sigma_part4[i]);
        // loop over all triangles in Omega
        for (std::size_t j = 0; j < Omega.size(); ++j) {
            
            const auto& T = Omega[j];
            auto edgeT=Boundary(T);
            // loop over all edges of triangle T
            for (std::size_t k = 0; k < edgeT.size(); ++k) {
                const auto& eT = edgeT[k];
                // loop over all boundary edges of sigma_i
                for (std::size_t e = 0; e < cntbi.size(); ++e) {
                    //if edge e of T is in boundary of sigma_i
                   // std::size_t tri = cntbi[e] / 3;
                   // std::size_t localEdge = cntbi[e] % 3;
                    if (eT == mbi[e]) {
                        Sigma_part4[i].push_back(T);
                        R.push_back(j, i, 1.);
                
            }
       

        }}

    }
    
}

}
    return {Sigma_part4,R};
}


// Grows the 4 partitions by nl boundary-touching layers.
// Returns (grown Sigma, R) with R a boolean (0/1) ne×4 membership matrix.
std::pair<Mesh2DPart, CooMatrix<double>>
Partition4(const Mesh2D& Omega, const std::size_t& nl)
{
    const std::size_t ne = Omega.size();

    // Base 4-way cut (compute once)
    auto [Sigma_part4, Q_part4] = Partition4(Omega);

    // Start R with the base membership (copy of Q_part4)
    CooMatrix<double> R;
    R.reserve(Q_part4.Nnz()); // if your COO has reserve()
    for (std::size_t j = 0; j < ne; ++j) {
        for (std::size_t p = 0; p < 4; ++p) {
            if (Q_part4(j, p) == 1.0) R.push_back(j, p, 1.0);
        }
    }

    // Grow nl layers
    for (std::size_t l = 0; l < nl; ++l) {
        for (std::size_t i = 0; i < 4; ++i) {
            // boundary of current submesh i (frozen for this layer)
            auto [mbi, cntbi] = Boundary(Sigma_part4[i]);

            // loop over all triangles in Omega
            for (std::size_t j = 0; j < ne; ++j) {
                // already in Sigma_i at current layer state? skip
                if (R(j, i) == 1.0) continue;

                const auto& T   = Omega[j];
                const auto edgeT = Boundary(T);

                bool touches = false;
                // any edge of T equals a boundary edge of Sigma_i?
                for (std::size_t k = 0; !touches && k < edgeT.size(); ++k) {
                    const auto& eT = edgeT[k];
                    for (std::size_t e = 0; e < cntbi.size(); ++e) {
                        std::size_t tri       = cntbi[e] / 3;
                        std::size_t localEdge = cntbi[e] % 3;
                        if (eT == mbi[tri][localEdge]) { touches = true; break; }
                    }
                }

                if (touches) {
                    Sigma_part4[i].push_back(T);
                    R.push_back(j, i, 1.0);
                }
            }
        }
    }

    return {Sigma_part4, R};
}


std::pair<std::vector<Mesh<2>>, CooMatrix<double>>
Partition4(const Mesh2D& Omega, const std::size_t& nl)
{
    constexpr std::size_t P = 4;
    const std::size_t ne = Omega.size();

    // Base 4-way (your existing centroid-based splitter)
    auto [Sigma, Q_part4] = Partition4(Omega);

    // R: same shape as Q_part4
    CooMatrix<double> R(NbRow(Q_part4), NbCol(Q_part4));
    R.reserve(Nnz(Q_part4));

    // Fast membership: in[i] == set of triangle indices in Σ_i
    std::array<std::unordered_set<std::size_t>, P> in;

    // Seed R and membership from Q_part4 by iterating its nonzeros
    for (const auto& jkv : Q_part4) {
        const auto& [j, p, v] = jkv;
        if (v != 0.0) {
            R.push_back(j, p, 1.0);
            in[p].insert(j);
        }
    }

    // Grow nl layers
    for (std::size_t l = 0; l < nl; ++l) {
        for (std::size_t i = 0; i < P; ++i) {

            // Boundary of current Σ_i (frozen this layer)
            auto [mbi, cntbi] = Boundary(Sigma[i]); // mbi: Mesh<1> (edges)

            // Try to add triangles not yet in Σ_i that touch its boundary
            for (std::size_t j = 0; j < ne; ++j) {
                if (in[i].count(j)) continue; // already in Σ_i

                const auto& T = Omega[j];
                auto edgesT = Boundary(T);    // std::array<Element<1>,3>

                bool touches = false;
                // Check each edge of T against EACH boundary edge of Σ_i
                for (std::size_t k = 0; !touches && k < edgesT.size(); ++k) {
                    const Element<1>& eT = edgesT[k];
                    for (std::size_t e = 0; e < cntbi.size(); ++e) {
                        // ✅ mbi[e] is already an Element<1> edge
                        if (eT == mbi[e]) { touches = true; break; }
                    }
                }

                if (touches) {
                    Sigma[i].push_back(T);
                    R.push_back(j, i, 1.0);
                    in[i].insert(j);
                }
            }
        }
    }

    return { Sigma, R };
}
    */
// Partition the mesh into 16 subdomains on a 4x4 grid.
// Sigma_p ~ ]j/4,(j+1)/4[ x ]k/4,(k+1)/4[, with j = p/4, k = p%4.
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

        // Map centroid to grid indices j,k ∈ {0,1,2,3}.
        // Threshold version (robust for open-interval intent).
        int jx = (cx < 0.25) ? 0 : (cx < 0.50) ? 1 : (cx < 0.75) ? 2 : 3;
        int ky = (cy < 0.25) ? 0 : (cy < 0.50) ? 1 : (cy < 0.75) ? 2 : 3;

        // p = 4*j + k  (so that j = p/4, k = p%4)
        int p  = 4 * jx + ky; // p in {0,...,15}

        Sigma[p].push_back(T);
        pattern.emplace_back(j, static_cast<std::size_t>(p));
    }

    // Build Q as a boolean (0/1) sparse matrix from the pattern
    CooMatrix<double> Q = BooleanMatrix(ne, 16, pattern);

    return {Sigma, Q};
}

std::pair<Mesh2DPart, CooMatrix<double>>
Partition16(const Mesh2D& Omega, const std::size_t& nl)
{
      // Compute base 16-way once
     auto [Sigma_part16, Q_part16] = Partition16(Omega);
      std::array<std::unordered_set<std::size_t>, 16> in;

    for (const auto& [row, col, val] : GetData(Q_part16)) {
        if (val != 0.0) {in[col].insert(row);
           // pattern.emplace_back(row, col);
        }
    }
        
    // Start R as a copy of Q (same nonzeros to begin with)
    CooMatrix<double> R =Q_part16;//= BooleanMatrix(NbRow(Q_part4), NbCol(Q_part4), pattern);

   

    const std::size_t ne = Omega.size();

    // Grow nl layers
    for (std::size_t l = 0; l < nl; ++l) {
        for (std::size_t i = 0; i < 16; ++i) {
            auto [mbi, cntbi] = Boundary(Sigma_part16[i]); // mbi: Mesh<1> of boundary edges

            // Try to add triangles not yet in Σ_i that touch its boundary
            for (std::size_t j = 0; j < ne; ++j) {
                if (in[i].count(j)) continue; // already in Σ_i
                
                const auto& T = Omega[j];
                auto edgesT   = Boundary(T);  // 3 edges (Element<1>)
                bool touches  = false;
                auto tri = 0;
                // edge of T vs boundary edges of Σ_i
                for (std::size_t k = 0; !touches && k < edgesT.size(); ++k) {
                    const Element<1>& eT = edgesT[k];
                    for (std::size_t e = 0; e < cntbi.size(); ++e) {
                        
                        
                        if (eT == mbi[e]) { touches = true; tri = cntbi[e] / 3; break; } // edge-to-edge compare
                    }
                }
                if (touches) {
                    Sigma_part16[i].push_back(T);
                    // put new mmbership in R
                    R.push_back(j, i, tri);   // append new membership
                    in[i].insert(j);          // mark as in Σ_i
                }
            }
        }
    }
    return {Sigma_part16, R};
 }

FeSpace2DxCoo Restrict(const FeSpace2D& Vh,const Mesh2D& Gamma,const std::vector<std::size_t>& tbl){
    FeSpace2D Vh_Gamma(Gamma);
    std::cout<<"Dim Vh: "<<dim(Vh)<<" Dim Vh_Gamma: "<<dim(Vh_Gamma)<<"\n";
    std::cout<<"Restriction table size: "<<tbl.size()<<"\n";
    CooMatrix<double> R(dim(Vh_Gamma),dim(Vh));

    for (std::size_t j=0; j<tbl.size(); ++j){
        std::size_t k = tbl[j]/(2+1);
        std::size_t l = tbl[j]%(2+1);    
        const auto& I = Vh_Gamma[j];
        auto II = Vh[k];
        for(std::size_t p=0; p<I.size(); ++p){
            R.push_back(I[p],II[l],1.);
        }        
    }
        
    return {Vh_Gamma,R};
}

std::vector<FeSpace2DxCoo>Partition4(const FeSpace2D& Vh, const std::size_t& nl){
    auto [Sigma_part4, R_part4] = Partition4G(Vh.mesh(),nl);
    std::vector<FeSpace2DxCoo> Vh_part4(4);
    for (std::size_t p=0; p<4; ++p){
        std::vector<std::size_t> tbl;
        for (const auto& [j,k,v] : GetData(R_part4)){
            tbl[v]=j;
        }
        Vh_part4[p] = Restrict(Vh,Sigma_part4[p],tbl);
    }
    return Vh_part4;

}



std::vector<FeSpace2DxCoo>
Partition16(const FeSpace2D& Vh, const std::size_t& nl){
    auto [Sigma_part16, R_part16] = Partition4G(Vh.mesh(),nl);
    std::vector<FeSpace2DxCoo> Vh_part16(16);
    for (std::size_t p=0; p<16; ++p){
        std::vector<std::size_t> tbl;
        for (const auto& [j,k,v] : GetData(R_part16)){
            tbl[v]=j;
        }
        Vh_part16[p] = Restrict(Vh,Sigma_part16[p],tbl);
    }
    return Vh_part16;

}