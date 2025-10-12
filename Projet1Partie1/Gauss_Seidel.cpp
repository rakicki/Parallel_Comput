
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
using namespace std;
// Code de calcul de la méthode de Gauss-Seidel avec une solution manufacturée cos(pi*x)*cos(2pi*y)
// On évalue ici aussi la convergence de la méthode.
double V(double y,double b){
    return 1-cos(2*M_PI*y/b);
}

double f(double x,double y) {
    return (-5*M_PI*M_PI)*cos(M_PI*x)*cos(2*M_PI*y);
}

void saveSolCSV(const vector<vector<double>>& grid,double dx,double dy, const string& filename) {
    ofstream f(filename);
    for (size_t j=0; j<grid.size(); ++j) {
        for (size_t i=0; i<grid[j].size(); ++i) {
            f << i*dx<<","<<j*dy<<","<<grid[j][i]<< "\n";
           
        }
    }
    f.close();
}
int main (int argc, char* argv[])
{
    //DATA
    int Nx=100;
    int Ny=100;
    int Nl=6000;
    double alpha=0.5;
    double a=1.0;
    double b=1.0;
   
    //CL
    double U0=2;
    double dx = a/(Nx+1.);
    double dy = b/(Ny+1.);
     double h=dx;
    vector<double> u_l(Nx+2,0);
    vector<vector<double>> pos(Ny+2,u_l);
    vector<vector<double>> sol_u(Ny+2,u_l);
    vector<vector<double>> sol_ex(Ny+2,u_l);
     vector<double> error_sol_Linf(Nl);
     vector<double> error_sol_L2(Nl);


     // CL et CI solution manufacturée
 for(int j=1;j<Ny+1;j++){
        for(int i=1;i<Nx+1;i++){
            sol_u[j][i]=0;
        }
        
    }
     for(int j=0;j<=Ny+1;j++){
        for(int i=0;i<=Nx+1;i++){
            sol_ex[j][i]=cos(M_PI*i*dx)*cos(2*M_PI*j*dy);
        }
        
    }
     for(int j=0;j<=Ny+1;j++){
        sol_u[j][0]=cos(2*M_PI*j*dy);
        sol_u[j][Nx+1]=cos(M_PI*(Nx+1)*dx)*cos(2*M_PI*j*dy);
    }
    for(int i=0;i<=Nx+1;i++){
        sol_u[0][i]=cos(M_PI*(i)*dx);
        sol_u[Ny+1][i]=cos(M_PI*(i)*dx)*cos(2*M_PI*(Ny+1)*dy);
    }
   /* CL du projet
   
     for(int j=1;j<=Ny;j++){
        for(int i=1;i<=Nx;i++){
            sol_u[j][i]=0;
        }
        
    }
    for(int j=0;j<=Ny+1;j++){
        sol_u[j][0]=U0*(1+alpha*V(j*dy,b));
        sol_u[j][Nx+1]=U0;
    }
    for(int i=0;i<=Nx+1;i++){
        sol_u[0][i]=U0;
        sol_u[Ny+1][i]=U0;
    }
   
   */
    for(int l=0;l<Nl;l++){
        double err_max=abs(sol_u[0][0]-sol_ex[0][0]); 
       double err_L2=0;
        for(int j=1;j<=Ny;j++){
            for(int i=1;i<=Nx;i++){
                sol_u[j][i] = ( dy*dy * (sol_u[j][i+1] + sol_u[j][i-1]) +dx*dx * (sol_u[j+1][i] + sol_u[j-1][i]) -f(dx*i, j*dy) * dx*dx * dy*dy) / (2 * (dx*dx + dy*dy));
                if(abs(sol_u[j][i]-sol_ex[j][i])>err_max){
                    err_max=abs(sol_u[j][i]-sol_ex[j][i]);
                }
                err_L2+=(sol_u[j][i]-sol_ex[j][i])*(sol_u[j][i]-sol_ex[j][i])*dx*dy;
            }
        }
         error_sol_Linf[l]=err_max;
        error_sol_L2[l]=sqrt(err_L2);
    }
       ofstream file;
  file.open("Gauss_Sidel_err.dat");
  for (int n=0; n<Nl; n++){
    file << n << "; " << error_sol_Linf[n] << endl;
  }
  file.close();
   ofstream file2;
    file2.open("Gauss_Sidel_errL2.dat");
  for (int n=0; n<Nl; n++){
    file2 << n << "; " << error_sol_L2[n] << endl;
  }
  file2.close();
    saveSolCSV(sol_u,dx,dy,"sol_Gauss_Sidel.csv");
    return 0;
}