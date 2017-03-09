! MPI TUTORIAL - solve heat equation in 1d --> FORTRAN 90 VERSION
!
! - A. Yeates, March 2017
!
!
program heat
    implicit none
    integer, parameter :: d=kind(0.d0)

    include 'mpif.h'  ! <-- note this

    integer, parameter :: nx=2048
    real(d), parameter :: k = 1.0_d
    real(d), parameter :: tMax = 0.1_d
    real(d) :: dx, dt, t, mu, uav, uavGlob
    integer :: ierr, nProcs, myRank, comm, nxtRank, prvRank, nxLocal, i
    real(d), allocatable :: x(:), g(:), u(:), xGlob(:), uGlob(:)
    real(d) :: start, finish

    ! (A) Initialize MPI
    ! ==================

    call mpi_init(ierr)

    ! - Get number of MPI processes and global rank:
    call mpi_comm_size(MPI_COMM_WORLD, nProcs, ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, myRank, ierr)
!    if (myRank==0) print*,nprocs,' processes'
!    print*,'Hello from process ',myRank

    ! - Set up Cartesian communicator:
    call mpi_cart_create(MPI_COMM_WORLD, 1, nProcs, 0, 1, comm, ierr)
    call mpi_comm_rank(comm, myRank, ierr)

    ! Get ranks of neighbours:
    call mpi_cart_shift(comm, 0, 1, prvRank, nxtRank, ierr)
!    print*,myRank,' has neighbours ',prvRank,' and ',nxtrank
!    print*,'MPI_PROC_NULL=',MPI_PROC_NULL

    ! (B) Initial condition
    ! =====================

    ! - Check that grid divides evenly between processes:
    if (mod(nx, nProcs)>0) then
        if (myRank==0) print*,nx," grid points don't divide evenly among",nProcs,"processors"
        stop
    end if
    nxLocal = nx/nProcs
    dx = 2.0_d/(real(nx,d) - 1.0_d)
!    print*,dx

    ! - One process defines coordinate array (global):
    if (myRank==0) then
        allocate(xGlob(1:nx))
        xGlob = (/ (-1.0_d + real(i)*dx, i=0,nx-1) /)
    end if

    ! - distribute chunks of coordinate array to each process:
    allocate(x(1:nxLocal))
    call mpi_scatter(xGlob, nxLocal, MPI_DOUBLE_PRECISION, x, nxLocal, MPI_DOUBLE_PRECISION, 0, comm, ierr)
!    print*,myRank,x

    ! - define initial condition function (locally for each process):
    allocate(g(1:nxLocal))
    g = exp(-10.0_d*x**2)
!    print*,g

    ! - Set timestep (for stability):
    dt = 0.2_d*dx**2/k
    mu = k*dt/(dx**2)

    ! -Initialize u (locally, including ghost cells):
    allocate(u(0:nxLocal+1))
    u = 0.0_d
    u(1:nxLocal) = g

    ! (C) Main loop
    ! =============

    start=MPI_Wtime();

    t = 0.0_d
    do while (t < tMax)
        ! - Enforce global boundary conditions:
        if (prvRank==MPI_PROC_NULL) u(0) = u(2)
        if (nxtRank==MPI_PROC_NULL) u(nxLocal+1) = u(nxLocal-1)

        ! - Communicate ghost values of u:
        ! -- send to left:
        if (prvRank/=MPI_PROC_NULL) call mpi_send(u(1), 1, MPI_DOUBLE_PRECISION, prvRank, 13, comm, ierr)
        if (nxtRank/=MPI_PROC_NULL) call mpi_recv(u(nxLocal+1), 1, MPI_DOUBLE_PRECISION, nxtRank, 13, comm, MPI_STATUS_IGNORE, ierr)
        ! -- send to right:
        if (nxtRank/=MPI_PROC_NULL) call mpi_send(u(nxLocal), 1, MPI_DOUBLE_PRECISION, nxtRank, 13, comm, ierr)
        if (prvRank/=MPI_PROC_NULL) call mpi_recv(u(0), 1, MPI_DOUBLE_PRECISION, prvRank, 13, comm,  MPI_STATUS_IGNORE, ierr)

        ! - Update u:
        u(1:nxLocal) = u(1:nxLocal) + mu*(u(0:nxLocal-1) - 2.0_d*u(1:nxLocal) + u(2:nxLocal+1))

        ! - Compute average of u:
        !uav = sum(u(1:nxLocal))
        !call mpi_allreduce(uav, uavGlob, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
        !uavGlob = uavGlob/real(nx,d)
        !if (myRank==0) print*,t,uavGlob

        ! - Update t:
        t = t + dt
    end do

    finish=MPI_Wtime();
    print*,'Parallel Elapsed time:',finish-start,' seconds'

    ! (D) Output results
    ! ==================

    ! - Collect results in control process:
    if (myRank==0) allocate(uGlob(1:nx))
    call mpi_gather(u(1:nxLocal), nxLocal, MPI_DOUBLE_PRECISION, uGlob, nxLocal, MPI_DOUBLE_PRECISION, 0, comm, ierr)

    ! - Save to file:
    if (myRank==0) then
        open(unit=1, file='u.dat', action='write', status='replace', form='unformatted', access='stream')
        write(1) nx
        write(1) xGlob
        write(1) uGlob
        close(unit=1)
    end if

    call mpi_finalize(ierr)

end program heat
