
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
using namespace std;
//Code de calcul de la méthode de Jacobi avec une solution manufacturée et de evaluation de sa convergence.
double V(double y,double b){
    return 1-cos(2*M_PI*y/b);
}

double f(double x,double y) {
    return (-5*M_PI*M_PI)*cos(M_PI*x)*cos(2*M_PI*y);
}


int main (int argc, char* argv[])
{
    //DATA
    int Nx=100;
    int Ny=100;
    int Nl=5000;
    double alpha=0.5;
    double a=1.0;
    double b=1.0;
     vector<double> pos(Nl);
    //CL
    double U0=2;
    double dx = a/(Nx+1.);
    double dy = b/(Ny+1.);
     double h=dx;
    vector<double> u_l(Nx+2,0);
   
    vector<vector<double>> sol_u(Ny+2,u_l);
    vector<vector<double>> sol_u_new(Ny+2,u_l);
    vector<double> error_sol_Linf(Nl);
     vector<double> error_sol_L2(Nl);
    vector<vector<double>> sol_ex(Ny+2,u_l);
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
     for(int j=0;j<=Ny+1;j++){
        sol_u_new[j][0]=cos(2*M_PI*j*dy);
        sol_u_new[j][Nx+1]=cos(M_PI*(Nx+1)*dx)*cos(2*M_PI*j*dy);
    }
    for(int i=0;i<=Nx+1;i++){
        sol_u_new[0][i]=cos(M_PI*(i)*dx);
        sol_u_new[Ny+1][i]=cos(M_PI*(i)*dx)*cos(2*M_PI*(Ny+1)*dy);
    }

    for(int l=0;l<Nl;l++){
       double err_max=abs(sol_u_new[0][0]-sol_ex[0][0]); 
       double err_L2=0;
        for(int j=1;j<=Ny;j++){
            for(int i=1;i<=Nx;i++){
                sol_u_new[j][i]=(dy*dy/2)*(sol_u[j][i+1]+sol_u[j][i-1])/(dx*dx+dy*dy)+ (dx*dx/2)*(sol_u[j+1][i]+sol_u[j-1][i])/(dx*dx+dy*dy)-dx*dx*dy*dy*f(dx*i,dy*j)/(2*(dx*dx+dy*dy));
                if(abs(sol_u_new[j][i]-sol_ex[j][i])>err_max){
                    err_max=abs(sol_u_new[j][i]-sol_ex[j][i]);
                }
                err_L2+=(sol_u_new[j][i]-sol_ex[j][i])*(sol_u_new[j][i]-sol_ex[j][i])*dx*dy;
            }
        }
        error_sol_Linf[l]=err_max;
        error_sol_L2[l]=sqrt(err_L2);
        sol_u.swap(sol_u_new);
    }
      ofstream file;
  file.open("Jacobi_err.dat");
  for (int n=0; n<Nl; n++){
    file << n << "; " << error_sol_Linf[n] << endl;
  }
  file.close();
   ofstream file2;
    file2.open("Jacobi_errL2.dat");
  for (int n=0; n<Nl; n++){
    file2 << n << "; " << error_sol_L2[n] << endl;
  }
  file2.close();
    return 0;
}