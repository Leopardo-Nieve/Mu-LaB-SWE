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
    use m_config
    implicit none ! had to write in a second line because VSCode was signalling an error

    ! configuration file
    type(CFG_t)           :: my_cfg

    ! declare local working variables
    integer:: itera_no
    double precision :: uo, vo,simTime, x_r, y_r, radius, r
    character:: td*24 ! get date for output
    character(len=3), dimension(2) :: all_forcing_schemes
    logical:: steadyFlow

    ! standard values
    call CFG_add(my_cfg, "physical_constants%nb", 0.012d0, "nb")
    call CFG_add(my_cfg, "simulation_parameters%steadyFlow", .false., "steadyFlow")
    call CFG_add(my_cfg, "simulation_parameters%epsilon", 0.0d0, "epsilon")
    call CFG_add(my_cfg, "simulation_parameters%consCriter", 1.0d-3, "consCriter")
    call CFG_add(my_cfg, "simulation_parameters%all_forcing_schemes", ["BG ", "GZS"], "all_forcing_schemes")
    call CFG_add(my_cfg, "simulation_parameters%forcing_scheme", "BG ", "all_forcing_schemes")
    call CFG_add(my_cfg, "simulation_parameters%simTime", 500.0d0, "simTime")
    call CFG_add(my_cfg, "simulation_parameters%itera_no", NINT(4e8), "itera_no")
    call CFG_add(my_cfg, "simulation_parameters%domainX", 14.0d3, "domainX")
    call CFG_add(my_cfg, "simulation_parameters%domainY", 70.0d0, "domainY")
    call CFG_add(my_cfg, "simulation_parameters%dx", 14.0d0, "dx")
    call CFG_add(my_cfg, "simulation_parameters%dt", 0.07d0, "dt")
    call CFG_add(my_cfg, "simulation_parameters%tau", 1.0d0, "tau")

    ! read parameters.cfg file to update values
    call CFG_read_file(my_cfg, "doc/parameters.cfg")

    r = 1.0d0 ! convergence ratio. start with 1, then 2, 4, 8, 16, 32

    ! define Manning's coefficient
!    nb = 0.012d0
    call CFG_get(my_cfg, "physical_constants%nb", nb)
    call CFG_get(my_cfg, "simulation_parameters%steadyFlow", steadyFlow)
    call CFG_get(my_cfg, "simulation_parameters%consCriter", consCriter)
    call CFG_get(my_cfg, "simulation_parameters%all_forcing_schemes", all_forcing_schemes)
    call CFG_get(my_cfg, "simulation_parameters%forcing_scheme", forcing_scheme)
    forcing_scheme = trim(forcing_scheme)
    call CFG_get(my_cfg, "simulation_parameters%simTime", simTime)
    call CFG_get(my_cfg, "simulation_parameters%itera_no", itera_no)
    call CFG_get(my_cfg, "simulation_parameters%domainX", domainX)
    call CFG_get(my_cfg, "simulation_parameters%domainY", domainY)
    call CFG_get(my_cfg, "simulation_parameters%dx", dx)
    call CFG_get(my_cfg, "simulation_parameters%dt", dt)
    call CFG_get(my_cfg, "simulation_parameters%tau", tau)

    ! initialize stopSim to let the simulation run
    stopSim = .false.

    ! initialize epsilon to let the simulation run
    if ( steadyFlow ) then
        call CFG_get(my_cfg, "simulation_parameters%epsilon", epsilon)
    end if

    ! debugging config fortran
    print*, "Manning constant: ", nb
    print*, "Steady flow: ", steadyFlow
    print*, "Epsilon: ", epsilon
    print*, "Consistency criterion: ", consCriter
    print*, "All forcing schemes: ", all_forcing_schemes
    print*, "Chosen forcing scheme: ", forcing_scheme
    print*, "Simulation time: ", simTime
    print*, "Iteration number: ", itera_no
    print*, "Domain x: ", domainX
    print*, "Domain y: ", domainY
    print*, "dx: ", dx
    print*, "dt: ", dt
    print*, "tau: ", tau

    ! Boundary conditions for inflow and outflow MUST BE LOWER CASE
    BCInflow  = "i" ! "i" (inflow) if assigned depth and velocity, otherwise "n" (Neumann) for zero gradient
    BCOutflow = "o" ! "o" (outflow) if assigned depth and velocity, otherwise "n" (Neumann) for zero gradient

    ! Forcing schemes: "bg": Buick-Greated, "gzs": Guo-Zheng-Shi
!    all_forcing_schemes = [character(len=3) :: "BG", "GZS"]



!    do i=1, size(all_forcing_schemes)
!        if (forcing_scheme == all_forcing_schemes(i)) exit
!        if (i == size(all_forcing_schemes)) then
!!            print*, "Please select a valid forcing scheme from: ", all_forcing_schemes
!!            stopSim = .true.
!        end if
!    end do

    if (forcing_scheme == "BG ") then
        B = 1.0d0-1.0d0/(2.0d0*tau)
        C = 0
    elseif (forcing_scheme == "GZS") then
        B = 1.0d0-1.0d0/(2.0d0*tau)
        C = 1.0d0-1.0d0/(2.0d0*tau)
    else
        print*, "Please select a valid forcing scheme from: ", all_forcing_schemes
        stopSim = .true.
    end if

    current_iteration = 0

    ! Initialise time
    time = 0

    ! assign a value of dy
    dy = dx ! m, lattice spacing

    ! define total number of nodes in x and y directions
    Lx = NINT(domainX/dx); Ly = NINT(domainY/dy) ! nodes

    ! allocate dimensions for dynamic arrays
    allocate (f(9,Lx,Ly),feq(9,Lx,Ly),ftemp(9,Lx,Ly),h(Lx,Ly),u(2,Lx,Ly),hLast(Lx,Ly),uLast(Lx,Ly),vLast(Lx,Ly),&
        & hCentered(2*Lx+1,2*Ly+1),uCentered(2*Lx+1,2*Ly+2),vCentered(2*Lx+1,2*Ly+1),&
        ! & C(Lx,Ly),Cz(2*Lx+1,2*Ly+1),Cb(2*Lx+1,2*Ly+1),tau_bx(2*Lx+1,2*Ly+1),&
        & force_x(2*Lx+1,2*Ly+1),force_y(2*Lx+1,2*Ly+1),&
        & H_part(2*Lx+1,2*Ly+1),zb(Lx,Ly),dzbdx(2*Lx+1,2*Ly+1), &
        & consInLft(1,Ly),consInRgt(1,Ly),consOutLft(1,Ly),consOutRgt(1,Ly),&
        & hAnal(Lx,Ly),uAnal(Lx,Ly),vAnal(Lx,Ly), &
        & force_x_MMS(2*Lx+1,2*Ly+1),force_y_MMS(2*Lx+1,2*Ly+1),S(9,Lx,Ly),force(2,9,Lx,Ly),&
        & u_vecLast(2,Lx,Ly))!, hIn(Ly), uIn(Ly))


    ! define pi
    pi = dacos(-1.0d0)

    ! call MMS_analytic_solution ! calculate analytical solution
    if (stopSim) STOP
    ! define bathymetry and node state array
    ! C = 0.0d0 ! m^2/s, assume all nodes are fluid nodes
    x_r    = 10.0d0 ! m, position of the bump in x direction
    y_r    = 5.0d0 ! m, position of the bump in y
    radius = 4.0d0 ! m, radius of the bump

    ! define partial depth and bathymetry
    zb = 0
    do x = 1, 2*Lx+1 ! to allow for body force scheme to have nodes in between each node
        position_x = dx*(DBLE(x-1)*0.5d0)

        ! partial depth
        H_part(x,:) = 50.5d0 - 40.0d0*position_x/domainX - 10.0d0*DSIN(pi*(4.0d0*position_x/domainX - 0.5d0))

        ! bathymetry
        zb(x,:) = H_part(1,:) - H_part(x,:)

        ! force_x_MMS(x,:) = 0.0d0 ! debug

        ! with bed
        ! force_x_MMS(x,:) = (1.0d0/9.0d0)*pi*gacl/domainX*(0.12d0*dsin(2.0d0*pi*position_x/domainX) &
        ! & + 8.0d0*dcos(2.0d0*pi*position_x/domainX))*&
        ! & (dsin(2.0d0*pi*position_x/domainX) + 3.0d0)

        ! without bed
        force_x_MMS(x,:) = (8.0d0/9.0d0)*pi*gacl/domainX*(dsin(2.0d0*pi*position_x/domainX) + 3.0d0)&
         & * dcos(2.0d0*pi*position_x/domainX)

        force_y_MMS(x,:) = 0.0d0

        ! if ( position_x > 0.8 .and. position_x < 1.2) then
        !     ! zb(x,:) = 0.2d0 - 0.05d0 * (position_x - 10.0d0)**2.0d0 ! bump function
        !     zb(x,:) = 0.2d-1 - 0.05d1 * (position_x - 10.0d-1)**2.0d0 ! bump function resized for 2 m x 2 m domain
        ! end if
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

    ! assign a value for the inlet discharge
    q_in = 4.42d0 ! m^2/s

    ! ho = 2.0d0 ! m, initial water depth

    uo = 0.0d0
    vo = 0.0d0

    dzbdx(2:2*Lx,:) = (zb(3:2*Lx+1,:) - zb(1:2*Lx-1,:))/(dx)
    dzbdx(1,:) = (-zb(3,:) + 4.0d0 * zb(2,:) - 3.0d0 * zb(1,:)) / (dx)
    dzbdx(2*Lx+1,:) = (3.0d0 * zb(2*Lx+1,:) - 4.0d0 * zb(2*Lx,:) + zb(2*Lx-1,:)) / (dx)

    ! constants for boundary conditions
    ! u_in = 0.3125d0
    ! u_out = -0.636d0
    ! q_in = 1d-2 ! debug
    ! q_in = 0.248*0.5d0 ! m^3/s, inlet discharge, symmetric domain, so divide by 2
    ! h()
    ! u_vec(1,1,:) = q_in/(h(1,:)*DBLE(domainY)) ! m/s, inlet velocity
    ! u_vec(1,1,:) = uAnal(1,:)
    ! u_vec(1,Lx,:) = uAnal(Lx,:)
    ! u_vec(1,1,:) = u_in
    ! u_vec(1,Lx,:) = u_out
    ! v(1,:) = vAnal(1,:)
    ! v(Lx,:) = vAnal(Lx,:)
    ! h(1,:) = 2.0d0
    ! hOut = 2.0d0 ! m, outflow depth
    ! h(Lx,:) = hOut ! set outflow depth
    ! u_vec(1,Lx,:) = 0.0d0

    ! assign a value for the molecular viscosity
    ! nu = 1.004d-6 ! m^2/s molecular viscosity of water


    ! calculate the minimum possible value of e such that the stationary population is positive
    ! eMin = dsqrt(5.0d0*gacl*ho/6.0d0 + 2.0d0/3.0d0*(q_in/ho)**2)

    ! define the lattice velocity
    e = dx/dt ! m/s, lattice velocity

    ! print values
    ! print*, "q inlet =", q_in, "m^2/s"
    print*, "e =", e, "m/s"
    print*, "dt =", dt, "s"
    ! print *,  "passed dt print" ! debug

    ! calculate the dimensionless relaxation time
    ! tau = 3.0d0*nu*dt/dx**2 + 0.5d0


    ! calculate molecular viscosity
    nu = (tau-0.5d0)*e*dx/3.0d0

    ! initialize the velocities
    ! u = q_in/(h*DBLE(domainY)) ! m/s, inlet velocity
    ! do i = 2, Lx-2
    !     u_vec(1,i,:) = (u_vec(1,Lx,:) - u_vec(1,1,:))/Lx * i + u_vec(1,1,:)
    !     v(i,:) = (v(Lx,:) - v(1,:))/Lx * i + v(1,:)
    ! end do

    ! define initial water depth
    do x = 1, Lx
        h(x,:) = H_part(2*x,Ly/2) ! m, initial water dept
    end do
    ! print *,  "passed h initializing" ! debug
    u(1,:,:) = uo
    u(2,:,:) = vo
    ! print*, "After initializing" ! debug
    ! print*, "h(1,:) =", h(1,:)" ! debug
    ! print *,  "reached setup" ! debug
    ! prepare the calculations
    call setup
    ! print*, "After setup"" ! debug
    ! print*, "u_vec(1,1,:) =", u_vec(1,1,:)" ! debug
    ! print*, "v(1,:) =", v(1,:)" ! debug

    ! do a=1,9
    !     do y=1,Ly
    !         print*, "f(",a,",1,",y,"): ", f(a,1,y)
    !     end do
    ! end do
    ! print *,  "passed setup" ! debug

    !  initialise errors
    L1_error = 0.0d0
    L2_error = 0.0d0

    ! main loop for time marching
    timStep: do

        time = time+dt
        current_iteration = current_iteration + 1

        call analytical_solution(time, Lx, Ly) ! update the analytical solution for the current timestep
        ! print *,  "passed analytical_solution" ! debug

        ! Update the body force with the current h
        call update_body_force
        ! print *,  "passed update_body_force" ! debug

        ! Apply no slip at solid boundary nodes to use the modified bounceback scheme (precollision)
        ! call Noslip_BC

        ! Streaming and collision steps
        call collide_stream
        ! print *,  "passed collide_stream" ! debug


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



        ! Apply Inflow and Outflow BC
        call Inflow_Outflow_BC
        ! print *,  "passed Inflow_Outflow_BC" ! debug


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

        ! do a=1,9 ! debug
            ! do y= 1, Ly ! debug
                ! print*, "f(a,1,",y,"):", f(a,1,y), "ftemp(a,1,",y,")", ftemp(a,1,y) ! debug
            ! end do ! debug
        ! end do ! debug

        ! Calculate h, u & v
        if (.not. stopSim) call solution
        ! print *,  "passed solution" ! debug

        if (.NOT. steadyFlow) call calculate_errors

        ! Update the feq
        call compute_feq
        ! print *,  "passed compute_feq" ! debug

        write(6,'(I8,A2,F20.14,A2,3(ES26.16,A2))') current_iteration,'   ', time, '   ',&
        & hAnal(1,Ly/2), '   ', h(1,Ly/2) ! debug
        ! & h(1,Ly/2), '   ', u_vec(1,1,Ly/2), '   ', v(1,Ly/2) ! commented for debug

        do i=1,Lx
            do j = 1, Ly

                ! make sure no u is NaN
                if (ieee_is_nan(u(1,i,j))) then
                    print*, "u",i,j,"is not a number"
                    stopSim = .true.
                end if

                ! make sure no v is NaN
                if (ieee_is_nan(u(2,i,j))) then
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
            if (steadyFlow) call calculate_errors ! only calculate final timestep, do not calculate twice if not steady
            call end_simulation
            exit
        end if

    end do timStep

    call ensure_results_directory ! ensures "../results" exists as a directory
    write(6,*)
    ! write(6,*)' Writing results in file: result.dat ... '
    ! open(66,file='../results/result.dat',status='unknown')    ! run from \src
    open(66,file='./results/result.dat',status='unknown')       ! run from \Mu-LaB-SWE
    call fdate(td)
    write(66,*) '# Date: ',td
    write(66,*) '# Fr =' ,u(1,1,Ly/2)/sqrt(gacl*h(1,Ly/2))
    write(66,*) '# tau =',tau,', uO =',uo
    write(66,*) '# Iteration No.: ',current_iteration
    write(66,'(1X,A6,I3,A9,I3)') '# Lx = ', Lx, ' Ly = ', Ly
    write(66,*) '#      Results of the computations'
    write(66,'(1X,A3,A4,A11,2A12) ') '# x','y','h(i,j)',&
                                        & 'u(1,i,j)' , 'u(1,i,j)'
    write(66,*) '#------------------------------------------'

    do x = 1, Lx
        do y = 1, Ly
            write(66,'(2i4,3f12.6)')x,y,h(x,y),u(1,x,y),u(2,x,y)
        end do
    end do
    close(66)

    ! Add after the existing result.dat write
    call write_csv
    write(6,*) ' CSV results written! ... '

end program main
