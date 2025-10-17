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
// send only even or only odd columns of a row
static inline void send_row_parity(const std::vector<double>& row,
                                   bool send_even, int dest, int tag, MPI_Comm comm) {
    const int N = (int)row.size();
    const int count = (N + (send_even ? 1 : 0)) / 2;     // ceil if starting at even
    MPI_Datatype T;
    MPI_Type_vector(count, 1, 2, MPI_DOUBLE, &T);        // stride in elements
    MPI_Type_commit(&T);
    MPI_Send(row.data() + (send_even ? 0 : 1), 1, T, dest, tag, comm);
    MPI_Type_free(&T);
}

// receive into a full row, placing values at even/odd columns directly
static inline void recv_row_parity(std::vector<double>& row,
                                   bool recv_even, int src, int tag, MPI_Comm comm) {
    const int N = (int)row.size();
    const int count = (N + (recv_even ? 1 : 0)) / 2;
    MPI_Datatype T;
    MPI_Type_vector(count, 1, 2, MPI_DOUBLE, &T);
    MPI_Type_commit(&T);
    MPI_Status st;
    MPI_Recv(row.data() + (recv_even ? 0 : 1), 1, T, src, tag, comm, &st);
    MPI_Type_free(&T);
}

int main(int argc, char* argv[]){
    //data
    int Nx= 2000;
    int Ny=32000;
    int Nl=200;
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
  MPI_Barrier(MPI_COMM_WORLD);
    double time1 = MPI_Wtime();
   for(int l=0;l<Nl;l++){
        int nblocs = (Nx + 3) / 2;
        MPI_Status status;
        const int j_top = nstarty;    
        const int j_bot = nendy;      
     
        if (myRank > 0) {
            bool black_even = (1 - (j_top % 2)) == 0;          // which parity is black at j_top
            send_row_parity(sol_u[1], black_even, myRank-1, 23, MPI_COMM_WORLD);

            bool red_even = (j_top % 2) == 0;             // our red parity at j_top
            recv_row_parity(sol_u[0], red_even, myRank-1, 22, MPI_COMM_WORLD);
        }


        if (myRank < nbTask-1) {
            bool red_even = (j_bot % 2) == 0;            
            recv_row_parity(sol_u[nendy-nstarty+2], red_even, myRank+1, 23, MPI_COMM_WORLD);

            bool black_even = (1 - (j_bot % 2)) == 0;       
            send_row_parity(sol_u[nendy-nstarty+1], black_even, myRank+1, 22, MPI_COMM_WORLD);
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


        if (myRank > 0) {
            bool red_even = (j_top % 2) == 0;                 
            send_row_parity(sol_u_new[1], red_even, myRank-1, 3, MPI_COMM_WORLD);

            bool black_even = (1 - (j_top % 2)) == 0;   
            recv_row_parity(sol_u_new[0],black_even, myRank-1, 2, MPI_COMM_WORLD);
        }

        if (myRank < nbTask-1) {
            bool black_even = (1 - (j_bot % 2)) == 0;
            recv_row_parity(sol_u_new[nendy-nstarty+2], black_even, myRank+1, 3, MPI_COMM_WORLD);

            bool red_even = (j_bot % 2) == 0;
            send_row_parity(sol_u_new[nendy-nstarty+1], red_even, myRank+1, 2, MPI_COMM_WORLD);
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


        sol_u=sol_u_new;
    }


MPI_Barrier(MPI_COMM_WORLD);
double time2 = MPI_Wtime();

if(myRank == 0) cout << "Duration: " << time2-time1 << endl;


 MPI_Finalize();
  return 0;
}