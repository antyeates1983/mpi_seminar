/*! MPI TUTORIAL - solve heat equation in 1d --> C++ VERSION
 -
 - A. Yeates, March 2017
 */
#include <mpi.h> // <-- note this
#include <math.h>
#include <iostream>
#include <fstream>
using namespace std;

int main(){
    static const int nx=2048;
    static const double k=1.0, tMax=0.1;
    int i;
    
    // (A) Initialize MPI
    // ==================
    
    MPI_Init(NULL, NULL);
    
    // - Get number of MPI processes and global rank:
    int nProcs, myRank;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    /*if(myRank==0) printf("%d processes\n", nProcs);
    printf("Hello from process %d\n", myRank);*/
    
    // - Set up Cartesian communicator:
    int dims[1], periodic[1];
    dims[0] = nProcs;
    periodic[0] = 0;
    MPI_Comm comm;
    MPI_Cart_create(MPI_COMM_WORLD, 1, dims, periodic, 1, &comm);
    MPI_Comm_rank(comm, &myRank);
    /*printf("Hello from process %d\n", myRank);*/
    
    // - Get ranks of neighbours:
    int prvRank, nxtRank;
    MPI_Cart_shift(comm, 0, 1, &prvRank, &nxtRank);
    /*printf("%d has neighbours %d and %d\n", myRank, prvRank, nxtRank);*/
    
    // (B) Initial condition
    // =====================
    
    // - Check that grid divides evenly between processes:
    if((nx % nProcs) > 0){
        if(myRank==0){
            printf("%d grid points don't divide evenly into %d processes\n", nx, nProcs);
            return 0;
            }
        }
    int nxLocal = nx/nProcs;
    double dx = 2.0/(double(nx) - 1.0);
    /*printf("%g\n",dx);*/
    
    // - One process defines coordinate array (global):
    double* xGlob; // --> need to define pointer so others can call MPI_Scatter.
    if(myRank==0){
        xGlob = new double[nx];
        for(i=0; i<nx; i++){
            xGlob[i] = -1.0 + double(i)*dx;
            }
        }
    
    // - Distribute chunks of coordinate array to each process:
    double x[nxLocal];
    MPI_Scatter(xGlob, nxLocal, MPI_DOUBLE, x, nxLocal, MPI_DOUBLE, 0, comm);
    /*printf("%d %g %g\n",myRank,x[0],x[nxLocal-1]);*/
    
    // - Define initial condition function (locally for each process):
    double g[nxLocal];
    for(i=0; i<nxLocal; i++){
        g[i] = exp(-10.0*x[i]*x[i]);
        }
    
    // - Set timestep (for stability):
    double dt = 0.2*dx*dx/k;
    double mu = k*dt/(dx*dx);
    
    // - Initialize u (locally, including ghost cells):
    double u[nxLocal+2];
    double uNew[nxLocal];
    u[0] = 0.0;
    u[nxLocal+1] = 0.0;
    for(i=1; i<=nxLocal; i++){
        u[i] = g[i-1];
        }

    // (C) Main loop
    // =============

    double start, finish;
    start=MPI_Wtime();
    
    double t = 0.0;
    double uav, uavGlob;
    while(t < tMax){
        // - Enforce global boundary conditions:
        if(prvRank==MPI_PROC_NULL) u[0] = u[2];
        if(nxtRank==MPI_PROC_NULL) u[nxLocal+1] = u[nxLocal-1];
        
        // - Communicate ghost values of u:
        // -- send to left:
        if(prvRank!=MPI_PROC_NULL) MPI_Send(&u[1], 1, MPI_DOUBLE, prvRank, 13, comm);
        if(nxtRank!=MPI_PROC_NULL) MPI_Recv(&u[nxLocal+1], 1, MPI_DOUBLE, nxtRank, 13, comm, MPI_STATUS_IGNORE);
        // -- send to right:
        if(nxtRank!=MPI_PROC_NULL) MPI_Send(&u[nxLocal], 1, MPI_DOUBLE, nxtRank, 13, comm);
        if(prvRank!=MPI_PROC_NULL) MPI_Recv(&u[0], 1, MPI_DOUBLE, prvRank, 13, comm, MPI_STATUS_IGNORE);
        
        // - Update u:
        for(i=1; i<=nxLocal; i++){
            uNew[i-1] = u[i] + mu*(u[i-1] - 2.0*u[i] + u[i+1]);
            }
         for(i=1; i<=nxLocal; i++){
            u[i] = uNew[i-1];
            }
        
        /*
        // - Compute average of u:
        uav = 0.0;
        for(i=1; i<=nxLocal; i++){
            uav = uav + u[i];
            }
        MPI_Allreduce(&uav, &uavGlob, 1, MPI_DOUBLE, MPI_SUM, comm);
        uavGlob = uavGlob/double(nx);
        if(myRank==0) printf("%f %f\n",t, uavGlob);
        */

        // - Update t:
        t += dt;
        }
    
    finish=MPI_Wtime();
    printf("Parallel Elapsed time: %f seconds\n", finish-start);
    
    // (D) Output results
    // ==================
    
    // - Collect results in control process:
    double uInterior[nxLocal];
    for(i=1; i<=nxLocal; i++){
        uInterior[i-1] = u[i];
        }
    double* uGlob;
    if(myRank==0) uGlob = new double[nx];
    MPI_Gather(uInterior, nxLocal, MPI_DOUBLE, uGlob, nxLocal, MPI_DOUBLE, 0, comm);
    
    // - Save to file:
    if(myRank==0){
        ofstream foo;
        foo.open("u.dat", ios::binary | ios::out);
        foo.write((char *)&nx, sizeof(int));
        for(i=0; i<nx; i++){
            foo.write((char *)&xGlob[i], sizeof(xGlob[i]));
            }
        for(i=0; i<nx; i++){
            foo.write((char *)&uGlob[i], sizeof(uGlob[i]));
            }
        foo.close();
        }

    MPI_Finalize();
    }
