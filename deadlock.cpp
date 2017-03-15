/*! MPI TUTORIAL - example where deadlock can occur --> C++ VERSION
 -
 - A. Yeates, March 2017
 */
#include <mpi.h> // <-- note this
#include <iostream>
using namespace std;

int main(){
    MPI_Init(NULL, NULL);

    int nProcs;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);

    // - Set up PERIODIC Cartesian communicator:
    int dims[1], periodic[1];
    dims[0] = nProcs;
    periodic[0] = 1;
    MPI_Comm comm;
    MPI_Cart_create(MPI_COMM_WORLD, 1, dims, periodic, 1, &comm);
    int myRank;
    MPI_Comm_rank(comm, &myRank);
    
    // - Get ranks of neighbours:
    int prvRank, nxtRank;
    MPI_Cart_shift(comm, 0, 1, &prvRank, &nxtRank);
    
    // - Everyone sends their value to the left:
    int myVal, theirVal;
    myVal = myRank;
    
    //MPI_Send(&myVal, 1, MPI_INT, nxtRank, 13, comm);
    MPI_Ssend(&myVal, 1, MPI_INT, nxtRank, 13, comm);
    MPI_Recv(&theirVal, 1, MPI_INT, prvRank, 13, comm, MPI_STATUS_IGNORE);
    
    printf("%i received from %i\n",myRank,nxtRank);
    
    MPI_Finalize();
    }
