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
    double precision :: uo, vo,simTime, x_r, y_r, radius, r, hAnal_x, dhanal_dx
    character:: fdate*24, td*24 ! get date for output
    logical :: steadyFlow, diagnosticFailure

    ! read parameters.cfg file to update values
    call CFG_read_file(my_cfg, "doc/parameters.cfg")

    r = 1.0d0 ! convergence ratio. start with 1, then 2, 4, 8

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
    call CFG_add_get(my_cfg, "simulation_parameters%alpha_MMS", alpha_MMS, &
        & "MMS force scaling factor: 0.0=off, 1.0=full [-]")

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

    ! ho = 2.0d0 ! m, initial water depth
    ho = hAnal(1,Ly/2) ! m, initial water depth
    uo = 0.0d0
    vo = 0.0d0

    ! initialize macroscopic variables
!    h = hAnal
!    u = uAnal
!    v = vAnal
    h = 2.0d0; u = 0.0d0; v = 0.0d0

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

    ! --- Smart parameter derivation from stability conditions ---
    ! Scan analytical solution for worst-case flow conditions
    h_max = maxval(hAnal)
    u_max = dsqrt(maxval(uAnal**2 + vAnal**2))

    ! Condition 3 (gh/e^2 < 1) is binding for subcritical flow (Fr < 1).
    ! S > 1 ensures gh/e^2 = 1/S^2, safely inside the stability boundary.
    e   = S * dsqrt(gacl * h_max)
    dt  = dx / e
    nu  = (tau - 0.5d0) * dx**2 / (3.0d0 * dt)

    ! Print derived values and stability diagnostics
    print*, "--- Derived LBM parameters ---"
    print*, "h_max =", h_max,                        "m"
    print*, "u_max =", u_max,                        "m/s"
    print*, "S     =", S,                            "[-]"
    print*, "e     =", e,                            "m/s"
    print*, "dt    =", dt,                           "s"
    print*, "tau   =", tau,                          "[-]"
    print*, "nu    =", nu,                           "m^2/s"
    print*, "u/e   =", u_max/e,                      "[-]  (stable if < 1)"
    print*, "gh/e^2=", gacl*h_max/e**2,              "[-]  (stable if < 1)"
    print*, "Fr    =", u_max/dsqrt(gacl*h_max),      "[-]  (subcritical if < 1)"

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
    ! Open stability diagnostic log
    open(newunit=stabilityUnit, file='stability_log.txt', &
         status='replace', action='write')

    flush(stabilityUnit)

    call setup

    ! main loop for time marching
    timStep: do

        ! Advance simulation counters
        time = time + dt
        current_iteration = current_iteration + 1

        ! --------------------------------------------------
        ! 1. Validate current macroscopic state
        ! --------------------------------------------------

        call audit_stability('PRE-COLLISION', diagnosticFailure)

        if (diagnosticFailure) then

            flush(stabilityUnit)
            stopSim = .true.

        end if

        ! --------------------------------------------------
        ! 2. Collision and streaming
        ! --------------------------------------------------

        if (.not. stopSim) then

            call update_body_force

            call collide_stream

            ! --------------------------------------------------
            ! 3. Validate proposed post-collision state
            ! --------------------------------------------------

            call validate_postcollision(diagnosticFailure)

            if (diagnosticFailure) then

                write(stabilityUnit,'(/,A)') &
                    'POST-COLLISION FAILURE DETECTED'

                ! Audit the previous valid macroscopic state
                call audit_stability('LAST VALID STATE', diagnosticFailure)

                write(stabilityUnit,'(A,I0)') &
                    'SIMULATION STOPPED AT ITERATION: ',current_iteration

                flush(stabilityUnit)

                stopSim = .true.

            end if

        end if

        ! --------------------------------------------------
        ! 4. Accept new state only if valid
        ! --------------------------------------------------

        if (.not. stopSim) then

            call solution

            ! Check reconstructed macroscopic variables
            do i = 1, Lx
                do j = 1, Ly

                    if (.not. ieee_is_finite(h(i,j)) .or. &
                        .not. ieee_is_finite(u(i,j)) .or. &
                        .not. ieee_is_finite(v(i,j))) then

                        write(stabilityUnit,'(A,I0,A,I0,A,I0)') &
                            'ITER=',current_iteration, &
                            ' CELL=(',i,',',j, &
                            '): NON-FINITE MACROSCOPIC VARIABLE'

                        stopSim = .true.

                    else if (h(i,j) <= 0.0d0) then

                        write(stabilityUnit,'(A,I0,A,I0,A,I0,A,ES16.8)') &
                            'ITER=',current_iteration, &
                            ' CELL=(',i,',',j, &
                            '): INVALID RECONSTRUCTED DEPTH = ',h(i,j)

                        stopSim = .true.

                    end if

                end do
            end do

            flush(stabilityUnit)

        end if

        ! --------------------------------------------------
        ! 5. Calculate equilibrium only if state is valid
        ! --------------------------------------------------

        if (.not. stopSim) then

            call compute_feq

            if (.not. all(ieee_is_finite(feq))) then

                write(stabilityUnit,'(A,I0)') &
                    'NON-FINITE EQUILIBRIUM POPULATION AT ITERATION: ', &
                    current_iteration

                flush(stabilityUnit)
                stopSim = .true.

            end if

        end if

        ! --------------------------------------------------
        ! 6. Print/checkpoint only if simulation is valid
        ! --------------------------------------------------

        if (.not. stopSim) then

            if (mod(current_iteration,per_print) == 0) then

                write(6,'(I0,A2,ES26.16)') &
                    current_iteration,'   ',h(1,Ly/2)

            end if

            if (current_iteration >= itera_no .or. time >= simTime) then

                print*, 'Maximum simulation time reached.'
                stopSim = .true.

            else if (mod(current_iteration,per_checkpoint) == 0) then

                call write_csv

            end if

        end if

        ! --------------------------------------------------
        ! 7. Stop or check convergence
        ! --------------------------------------------------

        if (stopSim .or. &
            (check_convergence(h,hLast,hLast2,epsilon) .and. &
             check_convergence(u,uLast,uLast2,epsilon) .and. &
             check_convergence(v,vLast,vLast2,epsilon))) then

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
