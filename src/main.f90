program main 
    use ieee_arithmetic
    use Mu_LaB_SWE
    implicit none
    integer :: itera_no
    double precision :: uo, vo, simTime
    double precision :: dx_ref, dt_ref
    character :: fdate*24, td*24
    logical :: steadyFlow
 
    ! define Manning's coefficient
    nb = 0.012d0
 
    steadyFlow = .TRUE.
 
    ! Boundary conditions for inflow and outflow MUST BE LOWER CASE
    BCInflow  = "i"
    BCOutflow = "o"
    ! initialize stopSim and epsilon to let the simulation run
    stopSim = .false.
    if (steadyFlow) then
        epsilon = 1.0d-3
    else
        epsilon = 0.0d0
    end if
    consCriter = 1.0d-3
    current_iteration = 0
    itera_no = NINT(4e8)
    time = 0.0d0
    simTime = 9.0d20
 
    ! ============================================================
    ! domain
    ! ============================================================
    domainX = 4.0d0
    domainY = 2.0d0*0.5d0
    ! ============================================================
    ! grid spacing
    ! Old very fine debugging mesh:
    ! dx = 0.00667d0
    ! dy = dx
    !
    ! New coarse base mesh for MMS convergence study:
    ! ============================================================
    dx = 0.1d0
    dy = dx
    ! number of nodes
    Lx = NINT(domainX/dx)
    Ly = NINT(domainY/dy)
 
    ! allocate arrays
    allocate (f(9,Lx,Ly), feq(9,Lx,Ly), ftemp(9,Lx,Ly), h(Lx,Ly), u(Lx,Ly), v(Lx,Ly), hLast(Lx,Ly), &
& hCentered(2*Lx+1,2*Ly+1), uCentered(2*Lx+1,2*Ly+2), vCentered(2*Lx+1,2*Ly+1), &
& force_x(2*Lx+1,2*Ly+1), force_y(2*Lx+1,2*Ly+1), &
& H_part(2*Lx+1,2*Ly+1), zb(2*Lx+1,2*Ly+1), dzbdx(2*Lx+1,2*Ly+1), &
& consInLft(1,Ly), consInRgt(1,Ly), consOutLft(1,Ly), consOutRgt(1,Ly), &
& hAnal(Lx,Ly), uAnal(Lx,Ly), vAnal(Lx,Ly))
 
    ! ============================================================
    ! bed slope
    ! ============================================================
    dzbdx = -6.25d-4
 
    ! define bed geometry
    do x = 1, 2*Lx+1
        position_x = dx * (DBLE(x-1) * 0.5d0)
        zb(x,:) = dzbdx(x,:) * (position_x - domainX)
    end do
 
    ! ============================================================
    ! initial flow constants
    ! ============================================================
    ho = 0.185d0
    vo = 0.0d0
 
    ! initialize depth
    do x = 1, Lx
        h(x,:) = ho - zb(2*x,Ly/2)
    end do
 
    ! ============================================================
    ! boundary-condition constants
    ! ============================================================
    q_in = 0.248d0 * 0.5d0
    hOut = 0.185d0
 
    ! initialize boundary values
    u(1,:) = q_in / (h(1,:) * DBLE(domainY))
    h(Lx,:) = hOut
 
    ! ============================================================
    ! time step
    !
    ! Old fixed dt:
    ! dt = 0.00145d0
    !
    ! New choice for MMS convergence:
    ! keep tau constant and scale dt ~ dx^2
    ! so that physical viscosity stays constant across mesh refinement
    ! ============================================================
    dx_ref = 0.1d0
    dt_ref = 0.00145d0
 
    dt = dt_ref * (dx/dx_ref)**2
 
    ! lattice velocity
    e = dx / dt
 
    print*, "dx =", dx, "m"
    print*, "dy =", dy, "m"
    print*, "Lx =", Lx, "nodes"
    print*, "Ly =", Ly, "nodes"
    print*, "e  =", e,  "m/s"
    print*, "dt =", dt, "s"
 
    ! relaxation time
    tau = 1.982d0
 
    ! molecular viscosity
    nu = (tau - 0.5d0) * e * dx / 3.0d0
    print*, "nu =", nu, "m^2/s"
 
    ! initialize full velocity field
    u = q_in / (h * DBLE(domainY))
    v = vo
 
    ! build analytical MMS fields before setup
    call MMS_analytic_solution
 
    ! prepare the calculations
    call setup
    ! ============================================================
    ! main loop for time marching
    ! ============================================================
    timStep: do
 
        time = time + dt
        current_iteration = current_iteration + 1
 
        ! Update the body force with the current h
        call update_body_force
 
        ! Streaming and collision steps
        call collide_stream
 
        do i = 1, Lx
            do j = 1, Ly
                do a = 1, 9
                    if (ieee_is_nan(ftemp(a,i,j))) then
                        print*, "ftemp", a, i, j, "is not a number"
                        stopSim = .true.
                    end if
                end do
            end do
        end do
 
        ! ========================================================
        ! BC in y
        ! Old:
        ! call Noslip_BC
        !
        ! MMS current choice is compatible with slip in y:
        ! ========================================================
        call Slip_BC
        ! Apply Inflow and Outflow BC
        call Inflow_Outflow_BC
 
        do i = 1, Lx
            do j = 1, Ly
                do a = 1, 9
                    if (ieee_is_nan(ftemp(a,i,j))) then
                        print*, "ftemp", a, i, j, "is not a number"
                        stopSim = .true.
                    end if
                end do
            end do
        end do
 
        ! Calculate h, u & v
        if (.not. stopSim) call solution
 
        ! Update feq
        call compute_feq
 
        write(6,'(I5,A2,3(ES26.16,A2))') current_iteration, '   ', h(1,Ly/2)
 
        do i = 1, Lx
            do j = 1, Ly
                if (ieee_is_nan(u(i,j))) then
                    print*, "u", i, j, "is not a number"
                    stopSim = .true.
                end if
                if (ieee_is_nan(v(i,j))) then
                    print*, "v", i, j, "is not a number"
                    stopSim = .true.
                end if
 
                if (ieee_is_nan(h(i,j))) then
                    print*, "h", i, j, "is not a number"
                    stopSim = .true.
                end if
            end do
        end do
 
        if (current_iteration >= itera_no .or. time >= simTime) then
            print*, "Maximum simulation time reached."
            stopSim = .true.
        end if
 
        if (stopSim .or. check_convergence(h, hLast, epsilon)) then
            call end_simulation
            exit
        end if
 
    end do timStep
 
    ! refresh analytical fields before writing results
    call MMS_analytic_solution
 
    call ensure_results_directory
    write(6,*)
 
    open(66, file='./results/result.dat', status='unknown')
    td = fdate()
    write(66,*) '# Date: ', td
    write(66,*) '# Fr =', u(1,Ly/2)/sqrt(gacl*h(1,Ly/2))
    write(66,*) '# tau =', tau, ', uO =', uo
    write(66,*) '# Iteration No.: ', current_iteration
    write(66,'(1X,A6,I3,A9,I3)') '# Lx = ', Lx, ' Ly = ', Ly
    write(66,*) '#      Results of the computations'
    write(66,'(1X,A3,A4,A11,2A12) ') '# x','y','h(i,j)', 'u(i,j)' , 'v(i,j)'
    write(66,*) '#------------------------------------------'
 
    do x = 1, Lx
        do y = 1, Ly
            write(66,'(2i4,3f12.6)') x, y, h(x,y), u(x,y), v(x,y)
        end do
    end do
    close(66)
 
    call write_csv
    write(6,*) ' CSV results written! ... '
end program main