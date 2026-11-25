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
    implicit none

    ! configuration file
    type(CFG_t)           :: my_cfg

    ! declare local working variables
    integer:: itera_no, per_checkpoint, per_print
    double precision :: uo, vo, simTime, x_r, y_r, radius, r
    double precision :: hAnal_x, dhanal_dx
    double precision :: h_scale
    character:: fdate*24, td*24 ! get date for output
    logical:: steadyFlow

    ! read parameters.cfg file to update values
    call CFG_read_file(my_cfg, "doc/parameters.cfg")

    ! define Manning's coefficient
    !nb = 0.012d0

    call CFG_add_get(my_cfg, "simulation_parameters%project_name", project_name, "Name of this project")
    call CFG_add_get(my_cfg, "simulation_parameters%steadyFlow", steadyFlow, "Whether the flow is steady or tidal")
    call CFG_add_get(my_cfg, "simulation_parameters%epsilon", epsilon, "convergence criterion")
    call CFG_add_get(my_cfg, "simulation_parameters%consCriter", consCriter, "conservation criterion")
    call CFG_add_get(my_cfg, "simulation_parameters%simTime", simTime, "total simulation time allowed [s]")
    call CFG_add_get(my_cfg, "simulation_parameters%itera_no", itera_no, "total number of iterations allowed")
    call CFG_add_get(my_cfg, "simulation_parameters%per_checkpoint_iter", per_checkpoint, "period between each checkpoint")
    call CFG_add_get(my_cfg, "simulation_parameters%per_print", per_print, "period between values printed to screen")
    call CFG_add_get(my_cfg, "simulation_parameters%domainX", domainX, "domain size in x direction [m]")
    call CFG_add_get(my_cfg, "simulation_parameters%domainY", domainY, "domain size in y direction [m]")
    call CFG_add_get(my_cfg, "simulation_parameters%dx", dx, "lattice spacing [m]")
    call CFG_add_get(my_cfg, "simulation_parameters%tau", tau, "dimensionless relaxation time [-]")
    call CFG_add_get(my_cfg, "simulation_parameters%S", S, "lattice velocity safety factor [-]")
    call CFG_add_get(my_cfg, "physical_constants%nb", nb, "Manning's coefficient")




    ! ==========================================================
    ! TEMPORARY DEBUGGING PARAMETERS
    ! ==========================================================

    call CFG_add_get(my_cfg, "simulation_parameters%refinement_ratio", &
        & r, "Grid refinement ratio")

    call CFG_add_get(my_cfg, "simulation_parameters%parameter_mode", &
        & parameter_mode, "Numerical parameter selection mode")

    call CFG_add_get(my_cfg, "simulation_parameters%profile_mode", &
        & profile_mode, "Initial macroscopic profile")

    call CFG_add_get(my_cfg, "simulation_parameters%nu_input", &
        & nu_input, "Prescribed kinematic viscosity [m^2/s]")

    call CFG_add_get(my_cfg, "simulation_parameters%h_test", &
        & h_test, "Uniform initial depth [m]")

    call CFG_add_get(my_cfg, "simulation_parameters%u_test", &
        & u_test, "Uniform initial x velocity [m/s]")

    call CFG_add_get(my_cfg, "simulation_parameters%v_test", &
        & v_test, "Uniform initial y velocity [m/s]")

    call CFG_add_get(my_cfg, "simulation_parameters%wave_amplitude", &
        & wave_amplitude, "Gravity wave depth perturbation [m]")

    call CFG_add_get(my_cfg, "simulation_parameters%alpha_hydro", &
        & alpha_hydro, "Hydrostatic force multiplier")

    call CFG_add_get(my_cfg, "simulation_parameters%alpha_MMS", &
        &  alpha_MMS, "MMS force scaling factor: 0.0=off, 1.0=full [-]")

    call CFG_add_get(my_cfg, "simulation_parameters%h_scale", &
        & h_scale, "Fixed reference depth for mode 2 [m]")

    ! ==========================================================
    ! END TEMPORARY DEBUGGING PARAMETERS
    ! ==========================================================


    ! Boundary conditions for inflow and outflow MUST BE LOWER CASE
    BCInflow  = "i" ! "i" (inflow) if assigned depth and velocity, otherwise "n" (Neumann) for zero gradient
    BCOutflow = "o" ! "o" (outflow) if assigned depth and velocity, otherwise "n" (Neumann) for zero gradient

    if (unix()) then
        filename_separator = "/"
    else
        filename_separator = "\\"
    end if

    ! initialize stopSim and epsilon to let the simulation run
    stopSim = .false.
!    if ( steadyFlow ) then
!        epsilon = max(epsilon_base / DBLE(r)**6.0, 1.0d-13)
!    else
!        epsilon = 0.0d0
!    end if

    current_iteration = 0
    time = 0

    ! assign a value of dx and dy
    dx = dx/r ! m, lattice spacing
    dy = dx ! m, lattice spacing

    if (h_scale <= 0.0d0) then
        print*, "ERROR: h_scale must be positive."
        stop
    end if

    if (tau <= 0.5d0) then
        print*, "ERROR: tau must be greater than 0.5."
        stop
    end if

    if (S <= 1.0d0) then
        print*, "ERROR: Base-grid S must be greater than 1."
        stop
    end if

    ! define total number of nodes in x and y directions
    Lx = NINT(domainX/dx); Ly = NINT(domainY/dy) ! nodes

    ! allocate dimensions for dynamic arrays
    allocate (f(9,Lx,Ly),feq(9,Lx,Ly),ftemp(9,Lx,Ly),h(Lx,Ly),u(Lx,Ly),v(Lx,Ly), &
        & hLast(Lx,Ly),uLast(Lx,Ly),vLast(Lx,Ly), &
        & hLast2(Lx,Ly),uLast2(Lx,Ly),vLast2(Lx,Ly), &
        & hCentered(2*Lx+1,2*Ly+1),uCentered(2*Lx+1,2*Ly+2),vCentered(2*Lx+1,2*Ly+1),&
        & force_x(2*Lx+1,2*Ly+1),force_y(2*Lx+1,2*Ly+1),&
        & H_part(2*Lx+1,2*Ly+1),zb(2*Lx+1,2*Ly+1),dzbdx(2*Lx+1,2*Ly+1), &
        & consInLft(1,Ly),consInRgt(1,Ly),consOutLft(1,Ly),consOutRgt(1,Ly),&
        & hAnal(Lx,Ly),uAnal(Lx,Ly),vAnal(Lx,Ly), &
        & force_x_MMS(2*Lx+1,2*Ly+1),force_y_MMS(2*Lx+1,2*Ly+1))


    ! define pi
    pi = dacos(-1.0d0)

    call MMS_analytic_solution ! calculate analytical solution
    if (stopSim) STOP


    ! ==========================================================
    ! SELECT INITIAL MACROSCOPIC PROFILE
    ! ==========================================================

    select case (profile_mode)

    case (0)
        ! Uniform stationary flow
        h = h_test
        u = 0.0d0
        v = 0.0d0

    case (1)
        ! Uniform moving flow
        h = h_test
        u = u_test
        v = v_test

    case (2)
        ! Sinusoidal depth perturbation: gravity-wave test
        do i = 1, Lx
            position_x = (DBLE(i)-0.5d0)*dx

            h(i,:) = h_test + wave_amplitude * &
                & dsin(2.0d0*pi*position_x/domainX)

            u(i,:) = u_test
            v(i,:) = v_test
        end do

    case (3)
        ! Current analytical MMS solution
        h = hAnal
        u = uAnal
        v = vAnal

    case default
        print*, "ERROR: Invalid profile_mode."
        stop

    end select

    ! Set reference analytical values for diagnostic output
    if (profile_mode /= 3) then
        hAnal = h
        uAnal = u
        vAnal = v
    end if

    ! Check initial conditions
    if (minval(h) <= 0.0d0) then
        print*, "ERROR: Non-positive initial water depth."
        stop
    end if

    print*, "Initial profile mode =", profile_mode
    print*, "Minimum initial depth =", minval(h)
    print*, "Maximum initial depth =", maxval(h)
    print*, "Maximum initial velocity =", &
        & maxval(dsqrt(u*u + v*v))

    ! ==========================================================
    ! END INITIAL PROFILE SELECTION
    ! ==========================================================



    ! define bathymetry and node state array
    ! C = 0.0d0 ! m^2/s, assume all nodes are fluid nodes
    x_r    = 10.0d0 ! m, position of the bump in x direction
    y_r    = 5.0d0 ! m, position of the bump in y
    radius = 4.0d0 ! m, radius of the bump

    ! define bathymetry
    zb = 0.0d0
    dzbdx = 0.0d0

!    do x = 1, 2*Lx+1
!        position_x = dx*(DBLE(x-1)*0.5d0)
!        force_x_MMS(x,:) = +(8.0d0/9.0d0)*pi*gacl/domainX &
!            & *(dsin(2.0d0*pi*position_x/domainX) + 3.0d0) &
!            & *dcos(2.0d0*pi*position_x/domainX)
!        force_y_MMS(x,:) = 0.0d0
!    end do

    ! Define bed geometry on half-node grid
    ! debug commented to identify source of instability
!    do x = 1, 2*Lx+1
!        position_x = dx*(DBLE(x-1)*0.5d0)
!        zb(x,:)    = 0.3d0*dsin(4.0d0*pi*position_x/domainX)
!        dzbdx(x,:) = 0.3d0*(4.0d0*pi/domainX)*dcos(4.0d0*pi*position_x/domainX)
!    end do

    do x = 1, 2*Lx+1
        position_x = dx*(DBLE(x-1)*0.5d0)
        do y = 1,2*Ly+1
            position_y = dx*(DBLE(y-1)*0.5d0)
            force_x_MMS(x,y) = (4.0d0/45.0d0)*pi*(gacl*(dsin(2.0d0*pi*position_x/domainX) &
                & + 3.0d0)**3.0d0*(10.0d0*dcos(2.0d0*pi*position_x/domainX) + 9.0d0*dcos( &
                & 4.0d0*pi*position_x/domainX)) - 135.0d0*dcos(2.0d0*pi*position_x/domainX) &
                & )/(domainX*(dsin(2.0d0*pi*position_x/domainX) + 3.0d0)**2.0d0)
            force_y_MMS(x,y) = 0.0d0
        end do
    end do



    !do x = 1, 2*Lx+1 ! to allow for body force scheme to have nodes in between each node
     !   position_x = dx*(DBLE(x-1)*0.5d0)
        ! commented to debug MMS
        !zb(x,:) =       0.02d0*dsin(3.1415926535897932d0*position_x/domainX)**2.0d0 ! 2 m wide bump function

        ! force_x_MMS(x,:) = 0.0d0 ! debug

	! with bed
        ! force_x_MMS(x,:) = (1.0d0/9.0d0)*pi*gacl/domainX*(0.12d0*dsin(2.0d0*pi*position_x/domainX) &
        ! & + 8.0d0*dcos(2.0d0*pi*position_x/domainX))*&
        ! & (dsin(2.0d0*pi*position_x/domainX) + 3.0d0)

        ! without bed
      !  force_x_MMS(x,:) = -(8.0d0/9.0d0)*pi*gacl/domainX*(dsin(2.0d0*pi*position_x/domainX) + 3.0d0) &
       !     & * dcos(2.0d0*pi*position_x/domainX)

        !force_y_MMS(x,:) = 0.0d0

        ! if ( position_x > 0.8 .and. position_x < 1.2) then
        !     ! zb(x,:) = 0.2d0 - 0.05d0 * (position_x - 10.0d0)**2.0d0 ! bump function
        !     zb(x,:) = 0.2d-1 - 0.05d1 * (position_x - 10.0d-1)**2.0d0 ! bump function resized for 2 m x 2 m domain
        ! end if
!    end do

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


!    do x = 1, Lx
!        do y = 1, Ly
!            h(x,y) = hAnal(x,y) - zb(2*x, 2*y) ! different array dimension
!        end do
!    end do

!    dzbdx(2:2*Lx,:) = (zb(3:2*Lx+1,:) - zb(1:2*Lx-1,:))/(dx)
!    dzbdx(1,:) = (-zb(3,:) + 4.0d0 * zb(2,:) - 3.0d0 * zb(1,:)) / (dx)
!    dzbdx(2*Lx+1,:) = (3.0d0 * zb(2*Lx+1,:) - 4.0d0 * zb(2*Lx,:) + zb(2*Lx-1,:)) / (dx)

    ! constants for boundary conditions
    ! u_in = 0.3125d0
    ! u_out = -0.636d0
    ! q_in = 1d-2 ! debug
    ! q_in = 0.248*0.5d0 ! m^3/s, inlet discharge, symmetric domain, so divide by 2
    ! h()
    ! u(1,:) = q_in/(h(1,:)*DBLE(domainY)) ! m/s, inlet velocity
    ! u(1,:) = uAnal(1,:)
    ! u(Lx,:) = uAnal(Lx,:)
    ! u(1,:) = u_in
    ! u(Lx,:) = u_out
    ! v(1,:) = vAnal(1,:)
    ! v(Lx,:) = vAnal(Lx,:)
    ! h(1,:) = 2.0d0
    ! hOut = 2.0d0 ! m, outflow depth
    ! h(Lx,:) = hOut ! set outflow depth
    ! u(Lx,:) = 0.0d0


    ! calculate the minimum possible value of e such that the stationary population is positive
    ! eMin = dsqrt(5.0d0*gacl*ho/6.0d0 + 2.0d0/3.0d0*(q_in/ho)**2)


    ! ==========================================================
    ! NUMERICAL PARAMETER SELECTION
    ! ==========================================================

    h_max = maxval(h)
    u_max = dsqrt(maxval(u**2 + v**2))

    if (tau <= 0.5d0) then
        print*, "ERROR: tau must be greater than 0.5."
        stop
    end if

    if (parameter_mode == 1) then

        ! MODE 1:
        ! Prescribe physical viscosity.
        ! Derive dt and lattice velocity.

        if (nu_input <= 0.0d0) then
            print*, "ERROR: nu_input must be positive."
            stop
        end if

        nu = nu_input

        dt = dx**2 * (2.0d0*tau - 1.0d0) / &
            & (6.0d0*nu)

        e = dx/dt

    elseif (parameter_mode == 2) then

    ! MODE 2 WITH DIFFUSIVE SCALING:
    ! S is prescribed on the base grid.
    ! The effective S increases with refinement ratio.
    ! Derive e, dt and physical viscosity.

    if (S <= 1.0d0) then
        print*, "ERROR: Base-grid S must be greater than 1."
        stop
    end if

    if (h_scale <= 0.0d0) then
        print*, "ERROR: h_scale must be positive."
        stop
    end if

    ! Scale lattice velocity linearly with refinement.
    e = (S*r)*dsqrt(gacl*h_scale)

    ! Diffusive scaling: dt is proportional to dx^2.
    dt = dx/e

    ! Physical viscosity implied by the LBM parameters.
    nu = dx**2*(2.0d0*tau-1.0d0)/(6.0d0*dt)

    else

        print*, "ERROR: Invalid parameter_mode."
        stop

    end if

    ! ==========================================================
    ! STABILITY DIAGNOSTICS
    ! ==========================================================

    print*, " "
    print*, "========== LBM PARAMETERS =========="

    print*, "dx      =", dx
    print*, "dt      =", dt
    print*, "e       =", e
    print*, "tau     =", tau
    print*, "nu      =", nu

    print*, "h_max   =", h_max
    print*, "u_max   =", u_max

    print*, "u/e     =", u_max/e
    print*, "gh/e^2  =", gacl*h_max/e**2

    print*, "Fr_max  =", u_max/dsqrt(gacl*h_max)

    print*, "Ma_max  =", dsqrt(3.0d0)*u_max/e

    print*, "===================================="

    print*, " "
    print*, "========== REFINEMENT CHECK =========="

    print*, "refinement ratio =", r
    print*, "base-grid S      =", S
    print*, "effective S      =", S*r
    print*, "reference depth  =", h_scale

    print*, "dx*r             =", dx*r
    print*, "dt*r^2           =", dt*r**2
    print*, "e/r              =", e/r
    print*, "nu               =", nu

    print*, "======================================"

    if (u_max/e >= 1.0d0) then
        print*, "WARNING: Velocity stability condition violated."
    end if

    if (gacl*h_max/e**2 >= 1.0d0) then
        print*, "WARNING: Celerity stability condition violated."
    end if

    if (u_max**2/(gacl*h_max) >= 1.0d0) then
        print*, "WARNING: Flow is not subcritical."
    end if

    ! ==========================================================
    ! END NUMERICAL PARAMETER SELECTION
    ! ==========================================================


    if (per_print == 0) then
        per_print = huge(1)
    end if

    if (per_checkpoint == 0) then
        per_checkpoint = huge(per_checkpoint)
!    else
!        per_checkpoint = per_checkpoint/dt ! convert from seconds to number of iterations
    end if

    ! calculate the dimensionless relaxation time
    ! tau = 3.0d0*nu*dt/dx**2 + 0.5d0

    ! initialize the velocities
    ! u = q_in/(h*DBLE(domainY)) ! m/s, inlet velocity
    ! do i = 2, Lx-2
    !     u(i,:) = (u(Lx,:) - u(1,:))/Lx * i + u(1,:)
    !     v(i,:) = (v(Lx,:) - v(1,:))/Lx * i + v(1,:)
    ! end do

    ! prepare the calculations
    call setup

    ! Save the initial populations as the reference solution.
    allocate(f_initial(9,Lx,Ly))
    f_initial = f

    print*, "Initial population reference saved."

    ! main loop for time marching
    timStep: do

        time = time+dt
        current_iteration = current_iteration + 1

        ! Update the body force with the current h
        call update_body_force

        ! Streaming and collision steps
        call collide_stream

        ! Stage 1: verify that stationary equilibrium is preserved.
        if (profile_mode == 0) then

            max_pop_error = maxval(abs(ftemp - f_initial))

            pop_scale = maxval(abs(f_initial))
            pop_tol = 100.0d0 * 2.2204460492503131d-16 &
                      * max(1.0d0, pop_scale)

            if (current_iteration == 1 .or. mod(current_iteration,100) == 0) then
                print*, "Iteration:", current_iteration
                print*, "Maximum population difference:", max_pop_error
                print*, "Allowed tolerance:", pop_tol
            end if

            if (max_pop_error > pop_tol) then
                print*, "ERROR: Stationary populations have changed!"
                stopSim = .true.
            end if

        end if

        ! Apply no slip at solid boundary nodes
        ! call Noslip_BC

        ! Apply Inflow and Outflow BC
        ! call Inflow_Outflow_BC

        ! Check populations for NaN or infinity.
        ! Report the first problematic node and all nine populations.
        check_populations: do j = 1, Ly
            do i = 1, Lx

                ! Check whether any of the nine populations is non-finite.
                do a = 1, 9
                    if (.not. ieee_is_finite(ftemp(a,i,j))) then

                        stopSim = .true.

                        print*, " "
                        print*, "========================================"
                        print*, "ERROR: Non-finite population detected!"
                        print*, "Iteration:", current_iteration
                        print*, "Time:", time, "s"
                        print*, "Node indices: i =", i, ", j =", j
                        print*, "Node coordinates:"
                        print*, "x =", (DBLE(i)-0.5d0)*dx, "m"
                        print*, "y =", (DBLE(j)-0.5d0)*dy, "m"
                        print*, "Problematic population: a =", a
                        print*, "========================================"

                        ! Print all nine populations at this node.
                        print*, "Distribution functions ftemp(a,i,j):"
                        do k = 1, 9
                            print*, "a =", k, ", ftemp =", ftemp(k,i,j)
                        end do

                        print*, "========================================"
                        print*, " "

                        stopSim = .true.

                        ! Stop searching: report only the first problematic node.
                        exit check_populations

                    end if
                end do

            end do
        end do check_populations

        ! Calculate h, u & v
        if (.not. stopSim) call solution

        ! ==========================================================
        ! STAGE 2: UNIFORM MOVING FLOW DIAGNOSTICS
        ! ==========================================================

        if (profile_mode == 1 .and. .not. stopSim) then

            max_h_error = maxval(abs(h - h_test))
            max_u_error = maxval(abs(u - u_test))
            max_v_error = maxval(abs(v - v_test))

            if (current_iteration == 1 .or. &
                & mod(current_iteration,100) == 0) then

                print*, " "
                print*, "Stage 2 - Iteration:", current_iteration
                print*, "Maximum depth error:", max_h_error
                print*, "Maximum u error:", max_u_error
                print*, "Maximum v error:", max_v_error

            end if

        end if

        ! Update the feq
        call compute_feq

        if (mod(current_iteration, per_print) == 0) write(6,'(I0,A2,3(ES26.16,A2))') current_iteration,'   ', h(1,Ly/2)

!        do i=1,Lx
!            do j = 1, Ly
!
!                ! make sure no u is NaN
!                if (ieee_is_nan(u(i,j))) then
!                    print*, "u",i,j,"is not a number"
!                    stopSim = .true.
!                end if
!
!                ! make sure no v is NaN
!                if (ieee_is_nan(v(i,j))) then
!                    print*, "v",i,j,"is not a number"
!                    stopSim = .true.
!                end if
!
!                ! make sure no h is NaN
!                if (ieee_is_nan(h(i,j))) then
!                    print*, "h",i,j,"is not a number"
!                    stopSim = .true.
!                end if
!            end do
!        end do

        if (current_iteration >= itera_no .or. time >= simTime) then
            print*, "Maximum simulation time reached."
            stopSim = .true. ! stop simulation after desired time reached
        elseif (mod(current_iteration, per_checkpoint) == 0) then
            call write_csv
        end if

        if (stopSim .or. &
            & (steadyFlow .and. &
            & check_convergence(h,hLast,hLast2,epsilon) .and. &
            & check_convergence(u,uLast,uLast2,epsilon) .and. &
            & check_convergence(v,vLast,vLast2,epsilon))) then

            call end_simulation
            exit

        end if

    end do timStep

    call ensure_directory("./results", len_trim("./results")) ! ensures "./results" exists as a directory    write(6,*)
    ! write(6,*)' Writing results in file: result.dat ... '
    ! open(66,file='../results/result.dat',status='unknown')    ! run from \src
    open(66,file='./results/result.dat',status='unknown')       ! run from \Mu-LaB-SWE
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
    call write_csv
    write(6,*) ' CSV results written! ... '

end program main
