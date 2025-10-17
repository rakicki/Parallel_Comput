#include <iostream>
#include <cmath>
#include <mpi.h>
#include <fstream>
#include <sstream>
#include <vector>
// Code de calcul de la méthode de Gauss Seidel parallèle avec une solution manufacturée cos(pi*x)*cos(2pi*y)
// On évalue ici la convergence de la méthode et sa performance.
using namespace std;

double V(double y,double b){
    return 1-cos(2*M_PI*y/b);
}

double f(double x,double y) {
    return (-5*M_PI*M_PI)*cos(M_PI*x)*cos(2*M_PI*y);
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
void send_even_indices(const vector<double>& sol, int N, int dest, int tag, MPI_Comm comm) {
    int nblocks = (N + 1) / 2; // nombre d'éléments pairs (indices 0,2,4,...)
    MPI_Datatype type_pairs;
    MPI_Type_vector(
        nblocks,     // count: nombre de blocs
        1,           // blocklength: 1 élément par bloc
        2,           // stride: pas de 2 éléments
        MPI_DOUBLE,  // type de base
        &type_pairs
    );
    MPI_Type_commit(&type_pairs);

    // Envoi des éléments d'indice pair
    MPI_Send(&sol[0], 1, type_pairs, dest, tag, comm);

    MPI_Type_free(&type_pairs);
}



void send_odd_indices(const std::vector<double>& sol, int N, int dest, int tag, MPI_Comm comm) {
  
    int nblocks = N / 2;            
    MPI_Datatype odd_type;

    MPI_Type_vector(
        nblocks,       // count
        1,             // blocklength
        2,             // stride in units of MPI_DOUBLE
        MPI_DOUBLE,
        &odd_type
    );
    MPI_Type_commit(&odd_type);

    if (nblocks > 0) {
        MPI_Send(&sol[1], 1, odd_type, dest, tag, comm);
    }

    MPI_Type_free(&odd_type);
}
int main(int argc, char* argv[]){
    //data
    int Nx= 100;
    int Ny=100;
    int Nl=1000;
    double alpha=0.5;
    double a=1.0;
    double b=1.0;
    

    //CL
    double U0=2;
    double dx = a/(Nx+1.);
    double dy = b/(Ny+1.);

    //MPI
    MPI_Init(&argc, &argv);
    int nbTask;
    int myRank;
    MPI_Comm_size(MPI_COMM_WORLD, &nbTask);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
   cout << "I am task " << myRank << " out of " << nbTask  << endl;
    //data per process



   int NsP=(Ny)/nbTask;
    int nstarty=myRank*NsP +1;
    int nendy=min((myRank+1)*NsP,Ny);
    vector<double> u_l(Nx+2 ,0);
    vector<vector<double>> sol_u( nendy-nstarty+3,u_l);
    vector<vector<double>> sol_u_new(nendy-nstarty+3,u_l);
      vector<vector<double>> sol_ex(nendy-nstarty+3,u_l);
    vector<double> error_sol_Linf(Nl);
     vector<double> error_sol_L2(Nl);
     for(int j=0;j<=nendy-nstarty+2;j++){
        for(int i=0;i<=Nx+1;i++){
             int j_global = nstarty + (j-1);
            double y = j_global * dy;
            sol_ex[j][i]=cos(M_PI*i*dx)*cos(2*M_PI*y);
        }
        
    }
    if(myRank==0){
        
        for(int i=0;i<=Nx+1;i++){
            sol_u_new[0][i]=cos(M_PI*(i)*dx);
            sol_u[0][i]=cos(M_PI*(i)*dx);
        }
        for(int j=1;j<=nendy-nstarty+2;j++){
        for(int i=1;i<=Nx;i++){
            sol_u[j][i]=0;
        }}
    }

    if(myRank==nbTask-1){
        for(int i=0;i<=Nx+1;i++){
            sol_u_new[nendy-nstarty+2][i]=cos(M_PI*(i)*dx)*cos(2*M_PI*(Ny+1)*dy);
            sol_u[nendy-nstarty+2][i]=cos(M_PI*(i)*dx)*cos(2*M_PI*(Ny+1)*dy);
           
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
        sol_u_new[j][0]=cos(2*M_PI*y);
        sol_u_new[j][Nx+1]=cos(M_PI*(Nx+1)*dx)*cos(2*M_PI*y); 
         sol_u[j][0]=cos(2*M_PI*y);
        sol_u[j][Nx+1]=cos(M_PI*(Nx+1)*dx)*cos(2*M_PI*y);
    }

   for(int l=0;l<Nl;l++){
    int nblocs = (Nx + 3) / 2;
    MPI_Status status;
        if(myRank>0){

             MPI_Send(&sol_u[1][0],Nx+2 , MPI_DOUBLE, myRank-1, 23, MPI_COMM_WORLD);
            MPI_Recv(&sol_u[0][0], Nx+2  , MPI_DOUBLE, myRank-1 , 22, MPI_COMM_WORLD, &status);
        }
 
        if(myRank<nbTask-1){
            
           MPI_Recv(&sol_u[nendy-nstarty+2][0], Nx+2 , MPI_DOUBLE, myRank+1 , 23, MPI_COMM_WORLD,&status);
      
            MPI_Send(&sol_u[nendy-nstarty+1][0], Nx+2, MPI_DOUBLE, myRank+1, 22, MPI_COMM_WORLD);
        }

        for(int j=1;j<=nendy-nstarty+1;j++){
            for(int i=1;i<=Nx;i++){
                int j_global = nstarty + (j-1);
                if((i+j_global)%2==0){// if red
                    double y = j_global * dy;
                    sol_u_new[j][i]=(dy*dy/2)*(sol_u[j][i+1]+sol_u[j][i-1])/(dx*dx+dy*dy)+ (dx*dx/2)*(sol_u[j+1][i]+sol_u[j-1][i])/(dx*dx+dy*dy)-dx*dx*dy*dy*f(dx*i,y)/(2*(dx*dx+dy*dy));

                }
            }
        }

        if(myRank>0){
           
        
             MPI_Send(&sol_u_new[1][0],Nx+2 , MPI_DOUBLE, myRank-1, 3, MPI_COMM_WORLD);
            MPI_Recv(&sol_u_new[0][0], Nx+2, MPI_DOUBLE, myRank-1 , 2, MPI_COMM_WORLD, &status);
        }
 
        if(myRank<nbTask-1){
            

            
            MPI_Recv(&sol_u_new[nendy-nstarty+2][0], Nx+2 , MPI_DOUBLE, myRank+1 , 3, MPI_COMM_WORLD,&status);

            MPI_Send(&sol_u_new[nendy-nstarty+1][0], Nx+2, MPI_DOUBLE, myRank+1, 2, MPI_COMM_WORLD);
        }

        for(int j=1;j<=nendy-nstarty+1;j++){
            for(int i=1;i<=Nx;i++){
                 int j_global = nstarty + (j-1);
                if((i+j_global)%2!=0){// if black
                    double y = j_global * dy;
                    sol_u_new[j][i]=(dy*dy/2)*(sol_u_new[j][i+1]+sol_u_new[j][i-1])/(dx*dx+dy*dy)+ (dx*dx/2)*(sol_u_new[j+1][i]+sol_u_new[j-1][i])/(dx*dx+dy*dy)-dx*dx*dy*dy*f(dx*i,y)/(2*(dx*dx+dy*dy));
                   
                }

            }
        }
   /* for (int j=1; j<=nendy-nstarty+1; j++) {
        int j_global = nstarty + (j-1);
        double y = j_global * dy;
        sol_u_new[j][0]=cos(2*M_PI*y);
        sol_u_new[j][Nx+1]=cos(M_PI*(Nx+1)*dx)*cos(2*M_PI*y); 
    }
    if (myRank == 0) {
        for (int i=0; i<=Nx+1; i++) {
            sol_u_new[0][i]=cos(M_PI*(i)*dx);}
    }
    if (myRank == nbTask-1) {
        for (int i=0; i<=Nx+1; i++) { 
            sol_u_new[nendy-nstarty+2][i]=cos(M_PI*(i)*dx)*cos(2*M_PI*(Ny+1)*dy);
        }
    }*/
        sol_u=sol_u_new;
    }


 


stringstream ss;
    ss << "sol_rank" << myRank << ".csv";
    saveSolCSV(sol_u, dx, dy, nstarty, Nx, ss.str());

 MPI_Finalize();
  return 0;
}