
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
using namespace std;
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
    int Nl=1000;
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
    vector<vector<double>> sol_u_new(Ny+2,u_l);
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
    for(int j=0;j<=Ny+1;j++){
        sol_u_new[j][0]=U0*(1+alpha*V(j*dy,b));
        sol_u_new[j][Nx+1]=U0;
    }
    for(int i=0;i<=Nx+1;i++){
        sol_u_new[0][i]=U0;
        sol_u_new[Ny+1][i]=U0;
    }
    for(int l=0;l<Nl;l++){
        for(int j=1;j<=Ny;j++){
            for(int i=1;i<=Nx;i++){
                sol_u_new[j][i]=(dy*dy/2)*(sol_u[j][i+1]+sol_u[j][i-1])/(dx*dx+dy*dy)+ (dx*dx/2)*(sol_u[j+1][i]+sol_u[j-1][i])/(dx*dx+dy*dy)-dx*dx*dy*dy*f(dx*i,dy*j)/(2*(dx*dx+dy*dy));
            }
        }
        sol_u.swap(sol_u_new);
    }
    saveSolCSV(sol_u,dx,dy,"sol_Jacobi.csv");
    return 0;
}