!----------------------------------------------------------! 
! main.f90 
! The file main.f90 in FORTRAN 90 is included in the book 
! to show how to use the module LABSWE.f90. Any FORTRAN-90 
! compiler may be used to compile, e.g. 
! "f90 LABSWE.f90 main.f90 -o labswe". 
! It simulates a flow in straight channel under constant 
! force with the periodic boundary conditions in the x 
! direction and no-slip boundary condition in y direction. 
! Consequently, a steady solution is obtained after 9000th 
! time steps, showing a typical laminar flow where a 
! parabolic distribution in velocity across the channel is 
! well developed. 
! J.G. Zhou, Peterborough, 2003 
!----------------------------------------------------------! 
! List of Major Variables 
! Fr - Froude number
! ho  - Initial water depth 
! itera no - Total iteration number or time steps 
! time - Time step or iteration counter 
! uo, vo - initial velocities 
!----------------------------------------------------------!
program main 
    
    ! call the module LABSWE 
    use ieee_arithmetic  ! Module for IEEE functions
    use LABSWE
    implicit none ! had to write in a second line because VSCode was signaling an error
    
    ! declare local working variables 
    integer:: itera_no
    double precision :: ho, uo, vo, time, epsilon!, simTime
    character:: fdate*24, td*24 ! get date for output

    ! initialize stopSim to let the simulation run
    stopSim = .false.
    epsilon = 1.0d-7

    ! constants for initializing flow field. 
    ho = 2.0d0
    uo = 0.0d0
    vo = 0.0d0

    ! assign a value for Zhou's relaxation time
    tauZhou = 1.5d0 

    ! assign a value for Zhou's lattice speed
    eZhou = 15.0d0 ! m/s

    ! assign a value for the inlet discharge
    qZhou = 4.42d0 ! m^2/s

    

    ! assign a value of dx and dy
    dx = 1.0d-1 ! m
    dy = dx

    ! define total lattice numbers in x and y directions
    domainX = 25.0d0     ; domainY = 5.0d0*dx ! dimensions in metres
    Lx = NINT(domainX/dx); Ly = NINT(domainY/dy) ! nodes

    

        ! calculate molecular viscosity 
    nuZhou = (tauZhou-0.5d0)*eZhou*dx/3.0d0
    print*, "nu Zhou =", nuZhou, "m^2/s"

    ! assign a value for the molecular viscosity
    ! nu = 1.004d-6 ! m^2/s molecular viscosity of water

    ! calculate the Reynolds number in Zhou's case
    ReZhou = qZhou/nuZhou
    print*, "ReZhou =", ReZhou

    ! calculate equivalent inlet discharge
    ! q_in = ReZhou*nu
    q_in=qZhou ! to replicate Zhou's results
    print*, "q inlet =", q_in, "m^2/s"
    
    ! calculate the minimum possible value of e such that the stationary population is positive
    eMin = dsqrt(5.0d0*gacl*ho/6.0d0 + 2.0d0/3.0d0*(q_in/ho)**2)
    print*, "e min =", eMin, "m/s"

    ! calculate the lattice velocity
    ! e = 2.0d0*eMin
    ! e=eZhou ! to replicate Zhou's results
    e=15.0d0 ! test to see if eMax is a valid stability condition
    print*, "e =", e, "m/s"

    ! define timestep dt
    dt = dx/e !s
    print*, "dt =", dt, "s"

    ! calculate the dimensionless relaxation time
    ! tau = 3.0d0*nu*dt/dx**2 + 0.5d0
    tau =tauZhou ! to replicate Zhou's results
    nu = (tau-0.5d0)*e*dx/3.0d0

    ! Total simulation time
    ! itera_no = 1
    ! itera_no = NINT(200/dt) 
    ! simTime = 200 ! assuming steady state after 200 s

    ! allocate dimensions for dynamic arrays
    allocate (f(9,Lx,Ly),feq(9,Lx,Ly),ftemp(9,Lx,Ly),h(Lx,Ly),& 
        & force_x(Lx,Ly),force_y(Lx,Ly),u(Lx,Ly),v(Lx,Ly),zb(Lx,Ly),dzbdx(Lx,Ly)) 
    
    ! initialize the depth and velocities 
    h = ho 
    u = uo 
    v = vo

    !define bed geometry
    zb = 0
    do x = 1, Lx
        if (x*dx > 8 .and. x*dx < 12) zb(x,:) = 0.2d0 - 0.05d0 * (x*dx - 10.0d0)**2.0d0 ! bump function
    end do

    dzbdx(2:Lx-1,:) = (zb(3:Lx,:) - zb(1:Lx-2,:)) / (2.0d0 * dx)
    dzbdx(1,:) = (-zb(3,:) + 4.0d0 * zb(2,:) - 3.0d0 * zb(1,:)) / (2.0d0 * dx)
    dzbdx(Lx,:) = (3.0d0 * zb(Lx,:) - 4.0d0 * zb(Lx-1,:) + zb(Lx-2,:)) / (2.0d0 * dx)


    !apply geometry
    h = h - zb


    !define initial velocity profile
    u(1,:) = q_in/h(1,:) 

    ! do y=1,Ly! debug
    !     print*, "u(1,",y,")=", u(1,y), "u(2,",y,")=", u(2,y)! debug
    ! end do! debug

    ! Set constant force
    force_x = h*gacl*dzbdx ! bed slope force
    force_y = 0.0d0
    ! print*, itera_no
    ! print*, ho
    ! print*, uo
    ! print*, vo
    ! print*, force_x, force_y
    ! prepare the calculations
    ! print*, "starting setup" !debugging
    call setup
    ! print*, "setup done" !debugging
    ! main loop for time marching
    time = 0
    timStep: do

        time = time+dt

        ! Streaming and collision steps
        call collide_stream

        ! debug first node negative populations
        ! do a=4,6
        !     if (ftemp(a,1,3)<0) then
        !         print*, "NEGATIVE population",a
        !         print*, "previous timestep: ",a,f(4,1+1,3-1)
        !         print*, "feq                ",a,feq(4,1+1,3-1)
        !         print*, "force term         ",a,dt_6e2*(ex(4)*force_x(1+1,3-1)+ey(4)*force_y(1+1,3-1))
        !     end if
        !     print*, "ftemp(",a,"1,3)=",ftemp(a,1,3)
        ! end do
        ! print*, "Please press enter to continue"
        ! read(*,*)

        ! Apply Slip BC
        ! call Slip_BC
        ! do a=1,9
        !     if (ftemp(a,1,2)/=ftemp(a,1,4)) then
        !         if (.not. stopSim) print*, "populations diverge after slip bc"
        !         print*, "ftemp(",a,"1,2)=", ftemp(a,1,2),"ftemp(",a,"1,4)=", ftemp(a,1,4)
        !         stopSim = .true.
        !     end if
        !     if (stopSim) stop
        ! end do
        ! print*, "after Slip_BC" !debug
        ! print*, "feq",3,"=", feq(3,1,3) !debug
        ! print*, "feq",7,"=", feq(7,1,3) !debug
        ! print*, "--------After Slip_BC---------" !debug
        ! do a=1,9 !debug
        !     print*,"f",a,"(1,3) =",ftemp(a,1,3) !debug
        ! end do !debug

        ! Apply Inflow and Outflow BC
        call Inflow_Outflow_BC

        ! Apply BCs to corners
        ! call Four_Corners_BC
        ! do a=1,9
        !     if (ftemp(a,1,2)/=ftemp(a,1,4)) then
        !         if (.not. stopSim) print*, "populations diverge after 4 corners"
        !         print*, "ftemp(",a,"1,2)=", ftemp(a,1,2),"ftemp(",a,"1,4)=", ftemp(a,1,4)
        !         stopSim = .true.
        !     end if
        !     if (stopSim) stop
        ! end do
        ! print*, "after Four_Corners_BC" !debug
        ! print*, "feq",3,"=", feq(3,1,3) !debug
        ! print*, "feq",7,"=", feq(7,1,3) !debug
        ! print*, "-----After Four_Corners_BC----" !debug
        ! do a=1,9 !debug
        !     print*,"f",a,"(1,3) =",ftemp(a,1,3) !debug
        ! end do !debug

        ! debug to determine which populations are negative
        do x=1,Lx
            do y=1,Ly
                do a=1,9
                    if (ftemp(a,x,y) < -1.0d-23) then
                        if (.not. stopSim) print*, "Error: negative population"
                        print*, "ftemp", a,x,y, "=", ftemp(a,x,y)
                        stopSim=.true.
                    end if
                end do
            end do
        end do

        ! Calculate h, u & v
        call solution

        ! Update the feq
        call compute_feq

        write(6,'(ES26.16,A2,3(ES26.16,A2))') time,'   ', h(100, Ly/2)

        ! if (time == 488 .or. time == 489) then
        !     do a = 1, 9
        !         print*, a, u(1,24), f(a,1,24), ex(a), ey(a)
        !         ! h(:,:) = h(:,:) + f(a,:,:)
        !         ! u(:,:) = u(:,:) + ex(a)*f(a,:,:)
        !         ! v(:,:) = v(:,:) +  ey(a)*f(a,:,:)
        !     end do
        ! end if

        ! do i=1,Lx 
        !     do j = 1, Ly
                
        !         ! Vérifier si u est un NaN
        !         if (ieee_is_nan(u(i,j))) then
        !             print *, "La variable u est un NaN.", i,j
                    
        !             stop
        !         end if
                
        !         ! Vérifier si v est un NaN
        !         if (ieee_is_nan(v(i,j))) then
        !             print *, "La variable v est un NaN.", i,j
        !             stop
        !         end if

        !         ! Vérifier si h est un NaN
        !         if (ieee_is_nan(h(i,j))) then
        !             print *, "La variable h est un NaN.", i,j
        !             stop
        !         end if
        !     end do
        ! end do
        
        ! if (time >= simTime) exit
        if (stopSim .or. check_convergence(u,epsilon)) then
            call end_simulation
            exit
        end if

    end do timStep
    print*, " h after loop =", h(100,Ly/2) !debugging

    write(6,*) 
    write(6,*)' Writing results in file: result.dat ... ' 
    open(66,file='../results_7.2.1/result.dat',status='unknown') 
    td=fdate() 
    write(66,*) '# Date: ',td 
    write(66,*) '# Fr =' ,u(1,Ly/2)/sqrt(gacl*h(1,Ly/2)) 
    write(66,*) '# tau =',tau,', uO =',uo 
    write(66,*) '# Iteration No.: ',itera_no 
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