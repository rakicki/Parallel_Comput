#include <iostream>
#include <cmath>
#include <mpi.h>
#include <fstream>
#include <sstream>
#include <vector>
// Code de Calcul Jacobi parallèle avec les conditions aux limites du projet.
// Le résultat de chaque processeur est écrit dans un fichier nommé "sol2rank<myrank>" .
// Pour visuliser le résultat final on concatène dans un fichier globale  les fichiers  "sol2rank<myrank>" que l'on visualise ave le fichier visu.ipynb.
using namespace std;

double V(double y,double b){
    return 1-cos(2*M_PI*y/b);
}

double f(double x,double y) {
    return 0;
    //(-5*M_PI*M_PI)*cos(M_PI*x)*cos(2*M_PI*y);
}
void saveSolCSV(const vector<vector<double>>& grid,double dx,double dy, int nstarty, int Nx, const string& filename) {
    ofstream f(filename);
    if(!f) return;
    int nrows = grid.size();   
    int cols  = grid[0].size(); 
   
    for (int jr = 1; jr <= nrows-2; ++jr) {
        int j_global = nstarty + (jr - 1); 
        double y = j_global * dy;
        for (int i = 0; i < cols; ++i) {
            double x = i * dx;
            f << x << "," << y << "," << grid[jr][i] << "\n";
        }
    }
    f.close();
}

int main(int argc, char* argv[]){
    //data
    int Nx=1000;
    int Ny=1000;
    int Nl=200;
    double alpha=0.5;
    double a=1.0;
    double b=1.0;
    

    //CL
    double U0=2;
    double dx = a/(Nx+1.);
    double dy = b/(Ny+1.);
double h=dx; 
    //MPI
    MPI_Init(&argc, &argv);

    int nbTask;
    int myRank;
    MPI_Comm_size(MPI_COMM_WORLD, &nbTask);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    cout << "I am task " << myRank << " out of " << nbTask  << endl;
    
    



   int NsP=(Ny)/nbTask;
    int nstarty=myRank*NsP +1;
    int nendy=min((myRank+1)*NsP,Ny);
    vector<double> u_l(Nx+2 ,0);
    vector<vector<double>> sol_u( nendy-nstarty+3,u_l);
    vector<vector<double>> sol_u_new(nendy-nstarty+3,u_l);
      
    if(myRank==0){
        
        for(int i=0;i<=Nx+1;i++){
            sol_u_new[0][i]=U0;
            sol_u[0][i]=U0;
        }
        for(int j=1;j<=nendy-nstarty+2;j++){
        for(int i=1;i<=Nx;i++){
            sol_u[j][i]=0;
        }}
    }

    if(myRank==nbTask-1){
        for(int i=0;i<=Nx+1;i++){
            sol_u_new[nendy-nstarty+2][i]=U0;
            sol_u[nendy-nstarty+2][i]=U0;
           
        }
         for(int j=0;j<nendy-nstarty+2;j++){
        for(int i=1;i<=Nx;i++){
            sol_u[j][i]=0;
        }}
    }
  if(myRank>0 && myRank<nbTask-1){
for(int j=0;j<nendy-nstarty+2; j++){
    for(int i=1;i<=Nx;i++){
    sol_u[j][i]=0;
}
}
  }
    for(int j=0;j<=nendy-nstarty+2;j++){
        int j_global = nstarty + (j-1);
    double y = j_global * dy;
        sol_u_new[j][0]=U0*(1+alpha*V(y,b));
        sol_u_new[j][Nx+1]=U0; 
         sol_u[j][0]=U0*(1+alpha*V(y,b));
        sol_u[j][Nx+1]=U0;
    }
   

   for(int l=0;l<Nl;l++){
        MPI_Status status;
        if(myRank>0){
             MPI_Send(&sol_u[1][0],Nx+2 , MPI_DOUBLE, myRank-1, 23, MPI_COMM_WORLD);
            MPI_Recv(&sol_u[0][0], Nx+2, MPI_DOUBLE, myRank-1 , 22, MPI_COMM_WORLD, &status);
        }
 
        if(myRank<nbTask-1){
            MPI_Recv(&sol_u[nendy-nstarty+2][0], Nx+2 , MPI_DOUBLE, myRank+1 , 23, MPI_COMM_WORLD,&status);
            MPI_Send(&sol_u[nendy-nstarty+1][0], Nx+2, MPI_DOUBLE, myRank+1, 22, MPI_COMM_WORLD);
        }
        for(int j=1;j<=nendy-nstarty+1;j++){
            for(int i=1;i<=Nx;i++){
                int j_global = nstarty + (j-1);
                double y = j_global * dy;
               
                sol_u_new[j][i]=(dy*dy/2)*(sol_u[j][i+1]+sol_u[j][i-1])/(dx*dx+dy*dy)+ (dx*dx/2)*(sol_u[j+1][i]+sol_u[j-1][i])/(dx*dx+dy*dy)-dx*dx*dy*dy*f(dx*i,y)/(2*(dx*dx+dy*dy));
            }
        }
        for (int j=1; j<=nendy-nstarty+1;j++) {
        int j_global=nstarty+(j-1);
        double y=j_global*dy;
        sol_u_new[j][0]=U0*(1 +alpha*V(y,b));
        sol_u_new[j][Nx+1]=U0;
    }
        if (myRank==0) {
            for (int i=0;i<=Nx+1;i++) {
                sol_u_new[0][i]=U0;}
        }
        if (myRank==nbTask-1) {
            for (int i=0;i<=Nx+1;i++){ 
                sol_u_new[nendy-nstarty+2][i]=U0;
            }
        }
        sol_u.swap(sol_u_new);
    }

    stringstream ss;
    ss << "sol2_rank" << myRank << ".csv";
    saveSolCSV(sol_u, dx, dy, nstarty, Nx, ss.str());



  return 0;
}