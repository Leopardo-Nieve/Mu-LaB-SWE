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
    double precision :: uo, vo,simTime,position_x!, position_y, x_r, y_r, r !debug
    character:: fdate*24, td*24 ! get date for output
    logical:: steadyFlow

    ! define Manning's coefficient
    nb = 0.012d0

    steadyFlow = .TRUE. ! if steady define `.true.`, if tidal define `.false.`

    ! Boundary conditions for inflow and outflow MUST BE LOWER CASE
    BCInflow  = "i" ! "i" (inflow) if assigned depth and velocity, otherwise "n" (Neumann) for zero gradient
    BCOutflow = "o" ! "o" (outflow) if assigned depth and velocity, otherwise "n" (Neumann) for zero gradient
    
    ! initialize stopSim and epsilon to let the simulation run
    stopSim = .false.
    if ( steadyFlow ) then
        epsilon = 1.0d-3
    else
        epsilon = 0.0d0
    end if
    consCriter = 1.0d-3
    
    current_iteration = 0
    ! itera_no = 1 !debug
    ! itera_no = 1e4 !debug
    itera_no = NINT(4e8)
        
    time = 0
    simTime = 9.0d20 ! s, maximum simulation time, set to a large value for steady flow

    ! assign a value for the inlet discharge
    ! q_in = 4.42d0 ! m^2/s

    ! define total lattice numbers in x and y directions
    domainX = 4.0d0 ! m
    domainY = 2.0d0*0.5d0 ! m, symetric domain
    
    ! assign a value of dx and dy
    dx = 0.00667 ! m, lattice spacing
    dy = dx ! m, lattice spacing
    
    ! define total number of nodes in x and y directions
    Lx = NINT(domainX/dx); Ly = NINT(domainY/dy) ! nodes

    ! allocate dimensions for dynamic arrays
    allocate (f(9,Lx,Ly),feq(9,Lx,Ly),ftemp(9,Lx,Ly),h(Lx,Ly),u(Lx,Ly),v(Lx,Ly),hLast(Lx,Ly),& 
        & hCentered(2*Lx+1,2*Ly+1),uCentered(2*Lx+1,2*Ly+2),vCentered(2*Lx+1,2*Ly+1),&
        ! & C(Lx,Ly),Cz(2*Lx+1,2*Ly+1),Cb(2*Lx+1,2*Ly+1),tau_bx(2*Lx+1,2*Ly+1),&
        & force_x(2*Lx+1,2*Ly+1),force_y(2*Lx+1,2*Ly+1),&
        & H_part(2*Lx+1,2*Ly+1),zb(2*Lx+1,2*Ly+1),dzbdx(2*Lx+1,2*Ly+1), &
        & consInLft(1,Ly),consInRgt(1,Ly),consOutLft(1,Ly),consOutRgt(1,Ly),&
        & hAnal(Lx,Ly),uAnal(Lx,Ly))!, hIn(Ly), uIn(Ly))

    dzbdx = -6.25d-4 ! m/m, slope of the bed
    
    ! define bathymetry and node state array
    ! C = 0.0d0 ! m^2/s, assume all nodes are fluid nodes
    ! x_r = 2.0d0 ! m, position of the cylinder in x direction
    ! y_r = 0.0d0 ! m, position of the cylinder in y
    ! r = 0.11d0 ! m, radius of the cylinder

    do x = 1, 2*Lx+1
        position_x = dx*(DBLE(x-1)*0.5d0)
        zb(x,:) = dzbdx(x,:)*(position_x - domainX) ! m, bed geometry
        ! do y = 1, 2*Ly+1
        !     position_y = dy*(DBLE(y-1)*0.5d0)
        !     if ( dsqrt((position_x - x_r)*(position_x - x_r) + (position_y - y_r)*(position_y - y_r)) <= r) then
        !         C(2*x,2*y) = 1.0d0 ! m^2/s, solid node, different array dimension
        !     end if
        ! end do
    end do

    ! determine boundary nodes
    ! do x = 1, Lx
    !     xf = x + 1
    !     xb = x - 1
    !     do y = 1, Ly
    !         if ( C(x,y) == 1 .OR. C(x,y) == 0.5) then
    !             cycle ! skip solid and boundary nodes
    !         end if

    !         yf = y + 1
    !         yb = y - 1
            
    !         if (C(xf,y) == 1 .OR. &
    !          & C(xf,yf) == 1 .OR. &
    !          & C(x,yf)  == 1 .OR. &
    !          & C(xb,yf) == 1 .OR. &
    !          & C(xb,y)  == 1 .OR. &
    !          & C(xb,yb) == 1 .OR. &
    !          & C(x,yb)  == 1 .OR. &
    !          & C(xf,yb) == 1) then
    !             C(x,y) = 0.5 ! m^2/s, boundary node
    !         end if
    !     end do
    ! end do


    ! constants for initializing flow field. 
    ho = 0.185d0 ! m, initial water depth
    ! uo = 0.0d0
    vo = 0.0d0

    ! initialize the depth 
    do x = 1, Lx
        h(x,:) = ho - zb(2*x,Ly/2) ! different array dimension
    end do

    ! dzbdx(2:2*Lx,:) = (zb(3:2*Lx+1,:) - zb(1:2*Lx-1,:))/(dx)
    ! dzbdx(1,:) = (-zb(3,:) + 4.0d0 * zb(2,:) - 3.0d0 * zb(1,:)) / (dx)
    ! dzbdx(2*Lx+1,:) = (3.0d0 * zb(2*Lx+1,:) - 4.0d0 * zb(2*Lx,:) + zb(2*Lx-1,:)) / (dx)

    ! constants for boundary conditions
    q_in = 0.248*0.5d0 ! m^3/s, inlet discharge, symmetric domain, so divide by 2
    ! h()
    u(1,:) = q_in/(h(1,:)*DBLE(domainY)) ! m/s, inlet velocity
    hOut = 0.185d0 ! m, outflow depth
    h(Lx,:) = hOut ! set outflow depth
    ! u(Lx,:) = 0.0d0
    
    ! assign a value for the molecular viscosity
    ! nu = 1.004d-6 ! m^2/s molecular viscosity of water

    
    ! calculate the minimum possible value of e such that the stationary population is positive
    ! eMin = dsqrt(5.0d0*gacl*ho/6.0d0 + 2.0d0/3.0d0*(q_in/ho)**2)

    ! define timestep dt
    dt = 0.00145d0 !s

    ! define the lattice velocity
    e = dx/dt ! m/s, lattice velocity

    ! print values
    ! print*, "q inlet =", q_in, "m^2/s"
    print*, "e =", e, "m/s"
    print*, "dt =", dt, "s"

    ! calculate the dimensionless relaxation time
    ! tau = 3.0d0*nu*dt/dx**2 + 0.5d0
    tau = 1.982d0 

    ! calculate molecular viscosity 
    nu = (tau-0.5d0)*e*dx/3.0d0

    ! initialize the velocities
    u = q_in/(h*DBLE(domainY)) ! m/s, inlet velocity
    v = vo

    ! prepare the calculations
    call setup
    
    ! main loop for time marching
    timStep: do

        time = time+dt
        current_iteration = current_iteration + 1

        ! Update the body force with the current h
        call update_body_force

        ! Streaming and collision steps
        call collide_stream

        do i=1,Lx
            do j=1,Ly
                do a=1,9
                    if (ieee_is_nan(ftemp(a,i,j))) then
                        print*, "ftemp",a,i,j,"is not a number"
                        stopSim = .true.
                    end if
                end do
            end do
        end do

        ! Apply no slip at solid boundary nodes
        call Noslip_BC
        
        ! Apply Inflow and Outflow BC
        call Inflow_Outflow_BC

        ! make sure no population is NaN
        do i = 1, Lx
            do j = 1, Ly
                do a = 1, 9
                    if ( ieee_is_nan(ftemp(a,i,j)) ) then
                        print*, "ftemp",a,i,j,"is not a number"
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
                    print*, "u",i,j,"is not a number"
                    stopSim = .true.
                end if
                
                ! make sure no v is NaN
                if (ieee_is_nan(v(i,j))) then
                    print*, "v",i,j,"is not a number"
                    stopSim = .true.
                end if

                ! make sure no h is NaN
                if (ieee_is_nan(h(i,j))) then
                    print*, "h",i,j,"is not a number"
                    stopSim = .true.
                end if
            end do
        end do
        if (current_iteration >= itera_no .or. time >= simTime) then
            print*, "Maximum simulation time reached."
            stopSim = .true. ! stop simulation after desired time reached
        end if
        if (stopSim .or. check_convergence(h,hLast,epsilon)) then
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
    write(6,*) ' CSV results written! ... '
    
end program main 