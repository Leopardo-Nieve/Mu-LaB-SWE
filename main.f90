!----------------------------------------------------------! 
! main.f90 
! The file main.f90 in FORTRAN 90 is adapted from the book 
! to use the module LABSWE.f90 to solve basic test 
! 7.2.1: Steady Flow over a Bump. Any FORTRAN-90 compiler
! may be used to compile, e.g. 
! "f90 LABSWE.f90 main.f90 -o labswe". 
! It simulates a steady flow in straight channel with a 
! defined inflow discharge and outflow depth boundary  
! conditions in the x direction and periodic boundary 
! condition in y direction. 
! Consequently, a steady solution is obtained after 11332 
! time steps, showing a the expected dip in the surface level profile 
! above the bump. 
! S. Fiset, Montreal, 2025
!----------------------------------------------------------! 
! List of Major Variables 
! ho  - Initial water depth 
! epsilon - Convergence criterion
! itera_no - Total iteration number or time steps 
! simTime - Maximum desired amount of simulation time
! time - Amount of time elapsed since start of simulation
! simTime - Maximum desired amount of simulation time
! uo, vo - Initial velocities 
!----------------------------------------------------------!
program main 
    
    ! call the module LABSWE 
    use ieee_arithmetic  ! Module for IEEE functions
    use Mu_LaB_SWE
    implicit none ! had to write in a second line because VSCode was signaling an error
    
    ! declare local working variables 
    integer:: itera_no
    double precision :: uo, vo,simTime, epsilon
    character:: fdate*24, td*24 ! get date for output
    logical:: steadyFlow

    ! define pi
    pi = dacos(-1.0d0)

    steadyFlow = .false. ! if steady define `.true.`, if tidal define `.false.`

    ! Boundary conditions for inflow and outflow MUST BE LOWER CASE
    BCInflow  = "n" ! "i" (inflow) if assigned depth and velocity, otherwise "n" (Neumann) for zero gradient
    BCOutflow = "n" ! "o" (outflow) if assigned depth and velocity, otherwise "n" (Neumann) for zero gradient
    
    ! initialize stopSim and epsilon to let the simulation run
    stopSim = .false.
    if ( steadyFlow ) then
        epsilon = 1.0d-10
    else
        epsilon = 0.0d0
    end if
    consCriter = 1.0d-3
    
    current_iteration = 0
    itera_no = 105.0d3
        
    time = 0
    simTime = 9117.5d0

    ! assign a value for the inlet discharge
    ! q_in = 4.42d0 ! m^2/s

    ! constants for initializing flow field. 
    uo = 0.0d0
    vo = 0.0d0

    ! define total lattice numbers in x and y directions
    domainX = 14.0d3     
    
    ! assign a value of dx and dy
    dx = domainX/800.0d0 ! m
    dy = dx
    domainY = 5.0d0*dx ! dimensions in metres
    Lx = NINT(domainX/dx); Ly = NINT(domainY/dy) ! nodes

    ! allocate dimensions for dynamic arrays
    allocate (f(9,Lx,Ly),feq(9,Lx,Ly),ftemp(9,Lx,Ly),h(Lx,Ly),& 
        & force_x(Lx,Ly),force_y(Lx,Ly),u(Lx,Ly),v(Lx,Ly),H_part(Lx,Ly),zb(Lx,Ly),dzbdx(Lx,Ly), &
        & consInLft(1,Ly),consInRgt(1,Ly),consOutLft(1,Ly),consOutRgt(1,Ly))!, hIn(Ly), uIn(Ly))

    do x = 1, Lx
        H_part(x,:) = 50.5d0 - 40.0d0*x*dx/domainX - 10.0d0*dsin(pi*(4.0d0*x*dx/domainX - 0.5d0))
    end do

    ! initialize the depth 
    h = H_part

    ! define bed geometry
    zb = H_part(1,Ly/2) - H_part


    dzbdx(2:Lx-1,:) = (zb(3:Lx,:) - zb(1:Lx-2,:)) / (2.0d0 * dx)
    dzbdx(1,:) = (-zb(3,:) + 4.0d0 * zb(2,:) - 3.0d0 * zb(1,:)) / (2.0d0 * dx)
    dzbdx(Lx,:) = (3.0d0 * zb(Lx,:) - 4.0d0 * zb(Lx-1,:) + zb(Lx-2,:)) / (2.0d0 * dx)

    ! constants for boundary conditions
    h(1,:) = h_in(time)
    ! hOut = 2.0d0
    u(Lx,:) = 0.0d0
    
    ! assign a value for the molecular viscosity
    ! nu = 1.004d-6 ! m^2/s molecular viscosity of water

    
    ! calculate the minimum possible value of e such that the stationary population is positive
    ! eMin = dsqrt(5.0d0*gacl*ho/6.0d0 + 2.0d0/3.0d0*(q_in/ho)**2)

    ! define the lattice velocity
    e = 200.0d0

    ! define timestep dt
    dt = dx/e !s
    ! print values
    ! print*, "q inlet =", q_in, "m^2/s"
    print*, "e =", e, "m/s"
    print*, "dt =", dt, "s"

    ! calculate the dimensionless relaxation time
    ! tau = 3.0d0*nu*dt/dx**2 + 0.5d0
    tau = 0.6d0 

    ! calculate molecular viscosity 
    nu = (tau-0.5d0)*e*dx/3.0d0

    ! initialize the velocities
    u = uo
    v = vo

    ! Set constant force
    force_x = -h*gacl*dzbdx ! bed slope force m^2/s^2
    force_y = 0.0d0 ! m^2/s^2
    ! prepare the calculations
    call setup
    
    ! main loop for time marching
    timStep: do

        time = time+dt
        current_iteration = current_iteration + 1

        ! Streaming and collision steps
        call collide_stream

        do i=1,Lx
            do j=1,Ly
                do a=1,9
                    if (ieee_is_nan(ftemp(a,i,j))) then
                        print*, "ftemp",a,x,y,"is not a number"
                        stopSim = .true.
                    end if
                end do
            end do
        end do

        ! Apply Inflow and Outflow BC
        call Inflow_Outflow_BC

        ! make sure no population is NaN
        do i = 1, Lx
            do j = 1, Ly
                do a = 1, 9
                    if ( ieee_is_nan(ftemp(a,i,j)) ) then
                        print*, "ftemp",a,x,y,"is not a number"
                        stopSim = .true.
                    end if
                end do
            end do
        end do

        ! Calculate h, u & v
        if (.not. stopSim) call solution

        ! Update the feq
        call compute_feq

        write(6,'(I5,A2,3(ES26.16,A2))') current_iteration,'   ', h(1,Ly/2)

        do i=1,Lx 
            do j = 1, Ly
                
                ! make sure no u is NaN
                if (ieee_is_nan(u(i,j))) then
                    print*, "u",x,y,"is not a number"
                    stopSim = .true.
                end if
                
                ! make sure no v is NaN
                if (ieee_is_nan(v(i,j))) then
                    print*, "v",x,y,"is not a number"
                    stopSim = .true.
                end if

                ! make sure no h is NaN
                if (ieee_is_nan(h(i,j))) then
                    print*, "h",x,y,"is not a number"
                    stopSim = .true.
                end if
            end do
        end do
        if (current_iteration >= itera_no .or. time >= simTime) then
            print*, "Maximum simulation time reached."
            stopSim = .true. ! stop simulation after desired time reached
        end if
        if (stopSim .or. check_convergence(u,h,epsilon)) then
            call end_simulation 
            exit
        end if

    end do timStep

    call ensure_results_directory ! ensures "../results" exists as a directory
    write(6,*) 
    write(6,*)' Writing results in file: result.dat ... ' 
    open(66,file='../results/result.dat',status='unknown') 
    td=fdate() 
    write(66,*) '# Date: ',td 
    write(66,*) '# Fr =' ,u(1,Ly/2)/sqrt(gacl*h(1,Ly/2)) 
    write(66,*) '# tau =',tau,', uO =',uo 
    write(66,*) '# Iteration No.: ',current_iteration 
    write(66,'(1X,A6,I3,A9,I3)') '# Lx = ', Lx, ' Ly = ', Ly 
    write(66,*) '#      Results of the computations' 
    write(66,'(1X,A3,A4,A11,2A12) ') '# x','y','h(i,j)',& 
                                        & 'u(i,j)' , 'v(i,j)' 
    write(66,*) '#------------------------------------------' 

    do x = 1, Lx 
        do y = 1, Ly 
            write(66,'(2i4,3f12.6)')x,y,h(x,y),u(x,y) ,v(x,y) 
        end do 
    end do 
    close(66) 

    ! Add after the existing result.dat write
    write(6,*) ' Writing CSV results in file: result.csv ... '
    call write_csv
    
end program main 