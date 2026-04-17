!----------------------------------------------------------! 
!                       LABSWE.f90 
! This module written in FORTRAN 90 implements the lattice 
! Boltzmann method for shallow water equations (LABSWE), 
! based on the 9-speed square lattices. It is adapted from 
! the sample code presented in appendix B of
! Zhou, J. G. (2004). Lattice Boltzmann methods for shallow 
! water flows. Springer. 
! https://doi.org/10.1007/978-3-662-08276-8
! Contrary to the sample code, the time step and lattice
! spacing do not have to be taken as units. Also the 
! module provides periodic boundary conditions, 
! inflow-outflow boundary conditions in the x direction 
! boundaries and no-slip boundary conditions in y direction 
! boundaries. 
! S. Fiset, Montreal, 2025 
! ----------------------------------------------------------!
!                   List of Major Variables
! a, x, y, b, i, j, k - Loop integers
! domainX, domainY - Domain's size in x and y directions [m]
! dt - time step [s]
! dx, dy - lattice spacing in x and y directions [m]
! e - lattice velocity [m/s]
! ex, ey - x and y components of particles' velocities
! f  - Distribution function [m]
! feq  - local equilibrium distribution function [m]
! force_x - x-direction component of force term [m^2/s^2]
! force_y - y-direction component of force term [m^2/s^2]
! ftemp - Temple distribution function [m^2/s^2]
! gacl - Gravitational acceleration [m/s^2]
! h - depth [m]
! hOut - fixed depth at outflow [m]
! H_part - partial depth [m]
! Lx, Ly - Total lattice numbers in x and y directions [-]
! nu - Molecular viscosity  [m^2/s]
! q_in - Inflow discharge [m^2/s]
! tau - Relaxation time [-]
! time - Amount of time elapsed since start of simulation [s]
! u, v - x and y components of flow velocity [m/s]
! uOut - fixed velocity at outflow [m/s]
! zb - bed geometry
module Mu_LaB_SWE

        implicit none 

        integer:: Lx,Ly,x,y,a,current_iteration, b,i,j,k,xf,yf,xb,yb
        integer, dimension(2):: hIndex
        logical:: stopSim, tauOk, velOk, celOk, FrOk
        character:: BCInflow, BCOutflow
        double precision:: ho,q_in,dx,dy,domainX,domainY,time,dt,eMin,e,tau,nu,hOut,&!,uOut & !necessary?
        &dt_6e2,one_8th_e4,one_3rd_e2,one_6th_e2,one_12th_e2, one_24th_e2,five_6th_g_e2,two_3rd_e2,gacl = 9.81,&
        & hMax, uMax2, FrMax, Fr, Ma, consCriter,pi,epsilon, nb, position_x, position_y, nu_MMs
        double precision, dimension(9):: ex,ey, eMax
        ! double precision, allocatable, dimension(:):: hIn,uIn ! not necessary?
        double precision, allocatable, dimension(:,:):: u,v,h,hLast,uLast,vLAst,hCentered,uCentered,vCentered,&
        & force_x,force_y,H_part,zb,dzbdx,consInLft,consInRgt,consOutLft,consOutRgt,hAnal,uAnal,vAnal,&
        & force_x_MMS,force_y_MMS!&
        ! &,C,Cz,Cb,tau_bx,tau_by,& !debug
        double precision, allocatable, dimension(:,:,:):: f,feq,ftemp 
    
contains 

subroutine setup 
    

    ! D2Q9 directions:
    ! 1 = E, 2 = NE, 3 = N, 4 = NW, 5 = W, 6 = SW, 7 = S, 8 = SE, 9 = Still

    ex = (/ 1.0d0,  1.0d0,  0.0d0, -1.0d0, -1.0d0, -1.0d0,  0.0d0,  1.0d0, 0.0d0 /)
    ey = (/ 0.0d0,  1.0d0,  1.0d0,  1.0d0,  0.0d0, -1.0d0, -1.0d0, -1.0d0, 0.0d0 /)
    
    ex(9) = 0.0d0; ey(9) = 0.0d0
    eMax = gacl*h(1,3)/3.0d0
    eMax = eMax + ex*ex*u(1,3)*u(1,3) + 2.0d0*ex*ey*u(1,3)*v(1,3) + ey*ey*v(1,3)*v(1,3)
    eMax = eMax - 1.0d0/3.0d0*(u(1,3)*u(1,3)+v(1,3)*v(1,3))
    eMax = eMax/-(ex*u(1,3) + ey*v(1,3))
    eMax = 6.0d0*eMax
    do a=1,9
        if (mod(a,2) == 0) eMax(a) = 2.5d-1*eMax(a) ! if even number index
    end do
    ex(:) = e*ex(:); ey(:) = e*ey(:) !scale for non unit lattice velocity

    ! constants to limit random error
    one_24th_e2=1.0d0/(24.0d0*e*e)
    one_12th_e2=2.0d0*one_24th_e2
    one_6th_e2 =2.0d0*one_12th_e2
    one_3rd_e2 =2.0d0*one_6th_e2
    one_8th_e4 = 1.0d0/(8.0d0*e*e*e*e)
    five_6th_g_e2 = 5.0d0*gacl*one_6th_e2
    two_3rd_e2 = 2.0d0*one_3rd_e2

    dt_6e2=dt/(6.0d0*e*e)

    ! determine initial inlet depth and velocity
    ! h(1,:) = 1.0d-3 ! m, initial depth at inlet
    ! h(1,:) = 2.0d0 ! m, initial depth at inlet
    ! u(1,:) = q_in/h(1,:) ! m/s, initial velocity at inlet

    ! initialize the depth over the entire domain
    ! do x = 1, Lx
    !     h(x,:) = ho - zb(2*x,Ly/2) ! different array dimension
    ! end do

    ! commented do loop and moved compute_feq out of it for MMS
    
    ! compute the equilibrium distribution function feq 
    call compute_feq
    
    ! initDepth: do
    !     ! compute the equilibrium distribution function feq 
    !     call compute_feq

    !     ftemp = feq ! initialize the temporary distribution function

    !     if ( .NOT. check_consistency("east",h,u,e,consCriter,Ly,.FALSE.) ) then 
    !         h = h + 1.0d-3 ! m, increase the depth at inlet by 1 mm
    !         ! commented because MMS requires new velocity inlet condition
    !         ! u(1,:) = q_in/(h(1,:)*DBLE(domainY)) ! m/s, update the velocity at inlet
    !     else
    !         print*, "Inlet consistency satisfied."
    !         exit initDepth ! exit loop if consistency is satisfied
    !     end if
    !     print*, "Initial depth at inlet:", h(1,Ly/2)
    ! end do initDepth

    stopSim = .false. ! reset stopSim flag before starting the simulation

    ! Set the initial distribution function to feq 
    f = feq
    return
end subroutine setup

subroutine update_body_force   
    ! interpolate values of h centred between each nodes to evaluate centred slope body force
    ! do x = 2, 2*Lx
    !     if (mod(x,2) == 0) then
    !         hCentered(x,:) = h(x/2,Ly/2)
    !     else
    !         hCentered(x,:) = (-h((x-3)/2,Ly/2)+6.0d0*h((x-1)/2,Ly/2) + 3.0d0*h((x+1)/2,Ly/2))/8.0d0
    !     end if

    !     ! if ( x>190 .AND. x<210 ) then !debug
    !     !     print*,"h_centred(",x,") =",hCentered(x,Ly/2) !debug
    !     ! end if!debug
    ! end do

    ! hCentered(1,:)    = (15.0d0*h(1,Ly/2) - 10.0d0*h(2,Ly/2) + 3.0d0*h(3,Ly/2))/8.0d0 
    ! hCentered(2*Lx+1,:) = (15.0d0*h(Lx,Ly/2) - 10.0d0*h(Lx-1,Ly/2) + 3.0d0*h(Lx-2,Ly/2))/8.0d0 

    hCentered = centred_interpolation(h,Lx,Ly)
    uCentered = centred_interpolation(u,Lx,Ly)
    vCentered = centred_interpolation(v,Lx,Ly)

    ! ! bed shear stress    
    ! Cz = hCentered**(1.0d0/6.0d0)/nb ! Chezy coefficient
    ! Cb = gacl/(Cz*Cz) ! bed friction coefficient
    ! tau_bx = Cb*uCentered*dsqrt(uCentered*uCentered + vCentered*vCentered) ! x-direction bed shear stress
    ! tau_by = Cb*vCentered*dsqrt(uCentered*uCentered + vCentered*vCentered) ! y-direction bed shear stress

    ! Set body force
    force_x = -hCentered*gacl*dzbdx + force_x_MMS !- tau_bx !debug ! m^2/s^2, bed slope force and bed shear stress
    force_y = 0.0d0 + force_y_MMS !-tau_by !debug ! m^2/s^2, bed shear stress
end subroutine update_body_force

subroutine collide_stream

    ! This calculates distribution function with the LABSWE 

    do y = 1, Ly 
        yf = y + 1
        yb = y -1
    
        do x = 1, Lx
            xf = x + 1
            xb = x -1
            ! if (C(x,y) == 0 .OR. C(x,y) == 0.5) cycle ! skip solid and boundary nodes

            ! Following 4 lines Implement periodic BCs in x or y directions
            if (xf > Lx) xf = xf - Lx !remove outlet periodic boundary
            if (xb < 1) xb = Lx + xb !remove inlet periodic boundary
            if (yf > Ly) yf = yf - Ly
            if (yb < 1) yb = Ly + yb 

            ! start streaming and collision 
            if (xf<=Lx) then ! periodic in y direction
                ftemp(1,xf,y) = f(1,x,y)-(f(1,x,y)-feq(1,x,y))/tau&
                & + dt_6e2*(ex(1)*force_x(2*x+1,2*y)+ey(1)*force_y(2*x+1,2*y))
            end if
            if (xf<=Lx) then !if (xf<=Lx .and. yf<=Ly) ! periodic in y direction
                ftemp(2,xf,yf) = f(2,x,y)-(f(2,x,y)-feq(2,x,y))/tau& 
                & + dt_6e2*(ex(2)*force_x(2*x+1,2*y+1)+ey(2)*force_y(2*x+1,2*y+1))
            end if
            ! if (yf<=Ly) ! periodic in y direction
            ftemp(3,x,yf) = f(3,x,y)-(f(3,x,y)-feq(3,x,y))/tau& 
                & + dt_6e2*(ex(3)*force_x(2*x,2*y+1)+ey(3)*force_y(2*x,2*y+1))
            if (xb>=1) then !if (xb>=1 .and. yf<=Ly) ! periodic in y direction
                ftemp(4,xb,yf) = f(4,x,y)-(f(4,x,y)-feq(4,x,y))/tau& 
                & + dt_6e2*(ex(4)*force_x(2*x-1,2*y+1)+ey(4)*force_y(2*x-1,2*y+1))
            end if
            if (xb>=1) ftemp(5,xb,y) = f(5,x,y)-(f(5,x,y)-feq(5,x,y))/tau& 
                & + dt_6e2*(ex(5)*force_x(2*x-1,2*y)+ey(5)*force_y(2*x-1,2*y))
            if (xb>=1) then !if (xb>=1 .and. yb>=1) ! periodic in y direction
                ftemp(6,xb,yb) = f(6,x,y)-(f(6,x,y)-feq(6,x,y))/tau& 
                & + dt_6e2*(ex(6)*force_x(2*x-1,2*y-1)+ey(6)*force_y(2*x-1,2*y-1))
            end if
            ! if (yb>=1) ! periodic in y direction
            ftemp(7,x,yb) = f(7,x,y)-(f(7,x,y)-feq(7,x,y))/tau& 
                & + dt_6e2*(ex(7)*force_x(2*x,2*y-1)+ey(7)*force_y(2*x,2*y-1))
            if (xf<=Lx) then !if (xf<=Lx .and. yb>=1) ! periodic in y direction
                ftemp(8,xf,yb) = f(8,x,y)-(f(8,x,y)-feq(8,x,y))/tau& 
                & + dt_6e2*(ex(8)*force_x(2*x+1,2*y-1)+ey(8)*force_y(2*x+1,2*y-1))
            end if
            ftemp(9,x,y) = f(9,x,y) - (f(9,x,y)-feq(9,x,y))/tau 
            
        end do 
    end do

return
end subroutine collide_stream

subroutine solution
    
    ! save last timestep
    hLast = h
    uLast = u
    vLast = v
    
    ! compute physical variables h, u and v

    ! Set the distribution function f
    f = ftemp

    ! compute the velocity and depth
    h = 0.0d0
    u = 0.0d0
    v = 0.0d0
    do a = 1, 9
        h(:,:) = h(:,:) + f(a,:,:)
        u(:,:) = u(:,:) + ex(a)*f(a,:,:)
        v(:,:) = v(:,:) +  ey(a)*f(a,:,:)
    end do
    u = u/h
    v = v/h

    return

end subroutine solution

subroutine compute_feq
    ! this computes the local equilibrium distribution function

    do a = 1, 8 
        ! if (mod(a,2) == 0) then 
        feq(a,:,:) = gacl*h(:,:)*h(:,:)*one_24th_e2 +& 
            & h(:,:)*one_12th_e2*(ex(a)*u(:,:)+& 
            & ey(a)*v(:,:))+h(:,:)*one_8th_e4& 
            & *(ex(a)*u(:,:)*ex(a)*u(:,:)+& 
            & 2.0d0*ex(a)*u(:,:)*ey(a)*v(:,:)+& 
            & ey(a)*v(:,:)*ey(a)*v(:,:))-& 
            & h(:,:)*one_24th_e2*(u(:,:)*u(:,:)+& 
            & v(:,:)*v(:,:)) 
        ! end if

        if (mod(a,2) /= 0) feq(a,:,:) = 4.0d0*feq(a,:,:) ! if odd number index
    end do
    feq(9,:,:) = h(:,:) - five_6th_g_e2*h(:,:)*h(:,:) - &
             & two_3rd_e2*h(:,:)*(u(:,:)*u(:,:) + v(:,:)*v(:,:))
    return
end subroutine compute_feq

subroutine Noslip_BC

    ! this is for noslip boundary with Bounce back scheme

    ! for lower boundary
    do a = 2, 4 
        ftemp(a,:,1) = ftemp(a+4,:,1) 
    end do
    
    ! for upper boundary
    do a = 6, 8 
        ftemp(a,:,Ly) = ftemp(a-4,:,Ly) 
    end do

    ! for cylinder boundary
    ! do x = 1, Lx !debug
    !     do y = 1, Ly !debug
    !         if (C(x,y) == 0.5) then ! boundary node !debug
    !             xf = x + 1
    !             xb = x - 1
    !             yf = y + 1
    !             yb = y - 1

    !             ! apply bounce back scheme
    !             if (C(xf,y)  == 1) ftemp(5,x,y) = ftemp(1,x,y) ! E
    !             if (C(xf,yf) == 1) ftemp(6,x,y) = ftemp(2,x,y) ! NE
    !             if (C(x,yf)  == 1) ftemp(7,x,y) = ftemp(3,x,y) ! N
    !             if (C(xb,yf) == 1) ftemp(8,x,y) = ftemp(4,x,y) ! NW
    !             if (C(xb,y)  == 1) ftemp(1,x,y) = ftemp(5,x,y) ! W
    !             if (C(xb,yb) == 1) ftemp(2,x,y) = ftemp(6,x,y) ! SW
    !             if (C(x,yb)  == 1) ftemp(3,x,y) = ftemp(7,x,y) ! S
    !             if (C(xf,yb) == 1) ftemp(4,x,y) = ftemp(8,x,y) ! SE
    !         end if
    !     end do
    ! end do

    return 
end subroutine Noslip_BC 

subroutine Slip_BC

    ! this is for slip boundary with Bounce back scheme

    ! for lower boundary
    ftemp(2,:,1) = ftemp(8,:,1) 
    ftemp(3,:,1) = ftemp(7,:,1) 
    ftemp(4,:,1) = ftemp(6,:,1) 
    
    ! for upper boundary
    ftemp(8,:,Ly) = ftemp(2,:,Ly) 
    ftemp(7,:,Ly) = ftemp(3,:,Ly)
    ftemp(6,:,Ly) = ftemp(4,:,Ly)

    return 
end subroutine Slip_BC 

subroutine Inflow_Outflow_BC
    ! macroscopic values
    ! h(1,:) = h(2,:)
    ! h(Lx,:) = hOut ! m, fixed depth at outflow
    h(1,:) = hAnal(1,:)
    ! u(1,:) = q_in/(h(1,:)*DBLE(domainY)) ! m/s, inflow velocity
    ! u(1,:) = uAnal(1,:)
    ! u(Lx,:) = uAnal(Lx,:)
    ! u(1,:) = q_in/h(1,:)
    ! u(Lx,:) = q_in/h(Lx,:)
    v(1,:) = vAnal(1,:)
    v(Lx,:) = vAnal(Lx,:)
    u(1,:) = e - e/h(1,:)*(ftemp(3,1,:)+ftemp(7,1,:)+ftemp(9,1,:)+2.0d0*(ftemp(4,1,:)+ftemp(5,1,:)+ftemp(6,1,:)))
    ! uAnal = u_analytical(time, Lx, Ly)
    ! u(1,:) = uAnal(1,:)

    h(Lx,:) = hAnal(Lx,:) ! m, fixed depth at outflow
    ! h(Lx,:) = ftemp(3,Lx,:) + ftemp(7,Lx,:) + ftemp(9,Lx,:) + 2.0d0*(ftemp(1,Lx,:) + ftemp(2,Lx,:) + ftemp(8,Lx,:))/(1+u(Lx,:)/e)
    u(Lx,:) = -e + e/h(Lx,:)*(ftemp(3,Lx,:)+ftemp(7,Lx,:)+ftemp(9,Lx,:)&
        &+2.0d0*(ftemp(1,Lx,:)+ftemp(2,Lx,:)+ftemp(8,Lx,:))) ! consistency check equation

    if ( BCInflow == "i" ) then
        ! consInLft(1,:) = h(1,:)-ftemp(9,1,:) ! left side of the consistence equation
        ! do a = 3, 7
        !     consInLft(1,:) = consInLft(1,:) - ftemp(a,1,:)
        ! end do
        ! consInRgt(1,:) = h(1,:)*u(1,:)/e + ftemp(4,1,:) + ftemp(5,1,:) + ftemp(6,1,:) ! right side of the consistence equation
        ! do j = 1, Ly
        !     if ( abs(consInLft(1,j) - consInRgt(1,j)) > consCriter ) then
        !         print*, "consistency fails at node",1,j
        !         print*,consInLft(1,j),"/=", consInRgt(1,j)
        !         stopSim = .true.
        !     end if
        ! end do

        ! if ( .not. stopSim ) then
        
        ! consistence check
        ! if ( check_consistency("east",h,u,e,consCriter,Ly)) then
        if ( .TRUE. ) then !debug to omit consistency errors (continuity equation optional?)
            ! Following lines implement inflow BC (Zhou, p.59)
            ftemp(1,1,:) = ftemp(5,1,:) + 2.0d0*h(1,:)*u(1,:)/(3.0d0*e)
            ftemp(2,1,:) = h(1,:)*u(1,:)/(6.0d0*e) + ftemp(6,1,:) + 0.5d0*(ftemp(7,1,:) - ftemp(3,1,:))
            ftemp(8,1,:) = h(1,:)*u(1,:)/(6.0d0*e) + ftemp(4,1,:) + 0.5d0*(ftemp(3,1,:) - ftemp(7,1,:))
        end if
    elseif (BCInflow == "n") then
            ! Neumann BC at the inflow (p. 58)
            ftemp(1,1,:) = ftemp(1,2,:) ! neigbouring population
            ftemp(2,1,:) = ftemp(2,2,:) ! neigbouring population
            ftemp(8,1,:) = ftemp(8,2,:) ! neigbouring population
    else
        print*, "`BCInflow` variable incorrectly defined as:", BCInflow
        print*, "***Hint: the condition must be written all in lower case.***"
    end if

    if ( BCOutflow == "o" ) then
        ! consistence check
        consOutLft(1,:) = h(Lx,:)
        do a = 1, 9
            if (a >= 4 .and. a <=6 ) cycle
            consOutLft(1,:) = consOutLft(1,:) - ftemp(a,Lx,:)
        end do
        consOutRgt(1,:) = -h(Lx,:)*u(Lx,:)/e + ftemp(1,Lx,:) + ftemp(2,Lx,:) + ftemp(8,Lx,:)
        do j = 1, Ly
            if ( abs(consOutLft(1,j) - consOutRgt(1,j)) > consCriter ) then
                ! print*, "consistence fails at node",Lx,j ! commented because not important if continuous for MMS (debug)
                ! print*,consOutLft(1,j),"/=", consOutRgt(1,j) ! "" debug
                ! stopSim = .true. "" ! debug
            end if
        end do

        ! if ( .not. stopSim ) then
        if ( .TRUE. ) then ! debug
            ! Following lines implement outflow BC (Zhou, p.60)
            ftemp(5,Lx,:) = ftemp(1,Lx,:) - 2.0d0*h(Lx,:)*u(Lx,:)/(3.0d0*e)
            ftemp(4,Lx,:) = -h(Lx,:)*u(Lx,:)/(6.0d0*e) + ftemp(8,Lx,:) + 0.5d0*(ftemp(7,Lx,:) - ftemp(3,Lx,:))
            ftemp(6,Lx,:) = -h(Lx,:)*u(Lx,:)/(6.0d0*e) + ftemp(2,Lx,:) + 0.5d0*(ftemp(3,Lx,:) - ftemp(7,Lx,:))
        end if
    elseif (BCOutflow == "n") then
        ! Neumann BC at outflow (p. 58)
        ftemp(4,Lx,:) = ftemp(4,Lx-1,:) ! neigbouring population
        ftemp(5,Lx,:) = ftemp(5,Lx-1,:) ! neigbouring population
        ftemp(6,Lx,:) = ftemp(6,Lx-1,:) ! neigbouring population
    else
        print*, "`BCOutflow` variable incorrectly defined as:", BCOutflow
        print*, "***Hint: the condition must be written all in lower case.***"
    end if
end subroutine Inflow_Outflow_BC

subroutine ensure_results_directory
    implicit none
    logical :: exists
    integer :: ierr
    character(len=100) :: cmd

    inquire(file='./results', exist=exists)

    if (.not. exists) then
        cmd = 'mkdir ".\results"'
        ierr = system(cmd)
        if (ierr /= 0) then
            print *, 'Error: Could not create the directory.'
            stop
        else
            print *, 'Directory ".\results" created successfully.'
        end if
    ! else
    !     print *, 'Directory "..\results" already exists.'
    end if
end subroutine ensure_results_directory

subroutine write_csv
    implicit none
    integer :: io, try
    integer, parameter :: max_tries = 2
    character(len=100) :: fpath
    character(len=8)     :: date
    character(len=10)    :: t ! time
    character(len=5)     :: zone
    ! integer, dimension(8) :: d
    character(len=19) :: formatted_time
    ! call date_and_time (values=d)
    
    
    call date_and_time(DATE=date, TIME=t, ZONE=zone)
    formatted_time = date(1:4)//'-'//date(5:6)//'-'//date(7:8) &
               & //'T'//                             &
               & t(1:2)//'-'//t(3:4)//'-'//t(5:6)

    fpath = './results/'//formatted_time//'.csv'
    write(6,*) ' Writing CSV results in file: '//fpath//" ... "

    try = 1
    do
        ! Try to open file explicitly for writing
        open(67, file=fpath, status='unknown', action='write', iostat=io)

        if (io == 0) exit  ! File opened successfully

        if (try >= max_tries) then
            print *, "ERROR: Cannot open file: ", trim(fpath)
            print *, "The file may be open in another program like Excel."
            stop
        end if

        print *, "WARNING: File is locked. Please close it (e.g. in Excel)."
        print *, "Retrying in 15 seconds..."
        call sleep_seconds(15)
        try = try + 1
    end do

    ! Write simulation parameters
    ! write(67, '(A)') 'Date,Iteration No.,tau,ujuj/e^2,gh/e^2,Fr'
    ! write(67, '(A,",",I10,",",F10.5,",",F10.5,",",F10.5,",",F10.5)') &
    !     trim(fdate()), current_iteration, tau, uMax2/(e*e), gacl*hMax/(e*e), FrMax

    ! Write CSV header
    write(67, '(A)') 'x (nodes),y (nodes),x (m),y (m),h + zb (m),zb (m),h (m),u (m/s),v (m/s), q (m^2/s)&
    &, h analytical (m), u analytical (m/s), h analytical + zb (m) ,&
    ! to add MMS source terms
    & Fx MMS, Fy MMS&
    &'

    ! Write data points
    do x = 1, Lx
        do y = 1, Ly
            write(67,'(2(I5,","),2(F17.14,","),3(F17.14,","),7(F17.14,","), &
            ! to add MMS source terms
            & 2(F17.14,",")   &
            & )') &
                x, y, &
                dx*(DBLE(x)-0.5d0), dy*(DBLE(y)-0.5d0), &
                h(x,y) + zb(2*x,2*y), zb(2*x,2*y), h(x,y), &
                u(x,y), v(x,y), h(x,y)*u(x,y), hAnal(x,y), uAnal(x,y),  hAnal(x,y) + zb(2*x,2*y) ,&
                & force_x_MMS(2*x,2*y), force_y_MMS(2*x,2*y)
        end do
    end do
    close(67)
end subroutine write_csv

subroutine sleep_seconds(n)
    integer, intent(in) :: n
    character(len=20) :: command
    write(command, '(A,I0,A)') 'timeout /T ', n, ' >nul'
    call system(command)
end subroutine sleep_seconds

subroutine end_simulation
    tauOk = .false.
    velOk = .false.
    celOk = .false.
    FrOk  = .false.

    hIndex = maxloc(h)
    hMax = h(hIndex(1),hIndex(2))
    uMax2 = 0
    FrMax = 0
    do x=1,Lx
        do y=1,Ly
            if (u(x,y)**2+v(x,y)**2>uMax2) uMax2 = u(x,y)**2+v(x,y)**2 
            Fr = (u(x,y)**2 + v(x,y)**2)/h(x,y)
            if (Fr>FrMax) then
                FrMax=Fr
                i=x;j=y
            end if
        end do
    end do
    FrMax = dsqrt(FrMax/gacl)

    if (tau     > 0.5) tauOk = .true. ! stability condition 1
    if (uMax2   < e*e) velOk = .true. ! stability condition 2
    if (gacl*hMax<e*e) celOk = .true. ! stability condition 3
    if (FrMax   < 1)   FrOk  = .true. ! stability condition 4
    Ma = sqrt(uMax2)/(1.0d0/sqrt(3.0d0)*e)
    print*, "tau =",tau 
    if (tauOk) then
        print*, "tau ok!"
    else
        print*, "tau NOT OKAY!!"
    end if
    print*, "u_ju_j/e^2 =", uMax2/(e*e)
    if (velOk) then
        print*, "velocity ok!"
    else
        print*, "velocity NOT OKAY!!"
    end if
    print*, "gh/e^2 =", gacl*hMax/(e*e)
    if (celOk) then
        print*, "celerity ok!"
    else
        print*, "celerity NOT OKAY!!"
    end if
    print*, "Fr max = ",FrMax, "at node:",i,j
    if (FrOk) then
        print*, "Froude ok!"
    else
        print*, "Froude NOT OKAY!!"
    end if
    print*, "Ma max = ",Ma!, "at node:",uIndex
    if (Ma<0.3) then
        print*, "Mach ok!"
    else 
        print*, "Mach NOT OKAY!!"
    end if
end subroutine end_simulation

! double precision function h_in(currentTime)
!     implicit none
!     double precision, intent(in)    :: currentTime
!     h_in = H_part(1,Ly/2) + 4.0d0 - 4.0d0*dsin(pi*(4.0d0*currentTime/86.4d3 + 0.5d0))
! end function h_in

function centred_interpolation(originalArray, dimX, dimY) result(outputArray)
    implicit none
    integer,          intent(in)    :: dimX, dimY ! dimesions of the original array (Lx, Ly)
    double precision, intent(in)    :: originalArray(:,:)
    double precision                :: outputArray(2*dimX+1, 2*dimY+1)

    do i = 2, 2*dimX
        if (mod(i,2) == 0) then
            outputArray(i,:) = originalArray(i/2,:)
        else
            outputArray(i,:) = (-originalArray((i-3)/2,dimY/2)+6.0d0*originalArray((i-1)/2,dimY/2)&
             + 3.0d0*originalArray((i+1)/2,dimY/2))/8.0d0
        end if
    end do

    outputArray(1,:)        = (15.0d0*originalArray(1,dimY/2) - 10.0d0*originalArray(2,dimY/2)&
        & + 3.0d0*originalArray(3,dimY/2))/8.0d0 
    outputArray(2*dimX+1,:) = (15.0d0*originalArray(dimX,dimY/2)-10.0d0*originalArray(dimX-1,dimY/2)&
         + 3.0d0*originalArray(dimX-2,dimY/2))/8.0d0
end function centred_interpolation

! function h_analytical(currentTime, dimX, dimY) result(h_a)
!     implicit none
!     integer,          intent(in)    :: dimX, dimY
!     double precision, intent(in)    :: currentTime
!     double precision                :: h_a(dimX,dimY)
!     h_a = H_part + 4.0d0 - 4.0d0*dsin(pi*(4.0d0*currentTime/86.4d3+0.5d0))
! end function h_analytical

! function u_analytical(currentTime, dimX, dimY) result(u_a)
!     implicit none
!     integer,          intent(in)    :: dimX, dimY
!     double precision, intent(in)    :: currentTime
!     double precision                :: u_a(dimX,dimY)

!     hAnal = h_analytical(currentTime, Lx, Ly)
!     do i = 1, Lx
!         u_a(i,:) = (i*dx - 14.0d3)*pi/(5.4d3*hAnal(i,:))*dcos(pi*(4.0d0*time/86.4d3 + 0.5d0))
!     end do
! end function u_analytical

logical function check_consistency(direction, hCheck, uCheck, eCheck, criterionCheck, dim, printErrors)
    implicit none
    character(len=*), intent(in) :: direction ! north, south, east, west
    double precision, intent(in) :: hCheck(:,:), uCheck(:,:)
    double precision, intent(in) :: eCheck, criterionCheck
    integer, intent(in)          :: dim
    logical, optional, intent(in):: printErrors
    logical                      :: localPrintErrors
    double precision             :: consLft(dim), consRgt(dim)

    ! If printErrors not specified, default to TRUE
    if (MOD(current_iteration,100) == 0) then
        if (present(printErrors)) then ! print errors every 100 iterations if specified
            localPrintErrors = printErrors
        else
            localPrintErrors = .true.
        end if
    else
        localPrintErrors = .false.
    end if
    check_consistency = .false. ! Initialize to false

    if (direction == "east") then
        ! consistence check
        consLft(:) = hCheck(1,:)-ftemp(9,1,:) ! left side of the consistence equation
        do a = 3, 7
            consLft(:) = consLft(:) - ftemp(a,1,:)
        end do
        consRgt(:) = hCheck(1,:)*uCheck(1,:)/eCheck &
        & + ftemp(4,1,:) + ftemp(5,1,:) + ftemp(6,1,:) ! right side of the consistence equation

        do i = 1, dim
            if ( abs(consLft(i) - consRgt(i)) > criterionCheck ) then
                stopSim = .true.
                ! stopSim = .false. ! debug
                ! if (localPrintErrors) then
                print*, "Inlet consistency check failed at node", i
                print*, consLft(i), "=/=", consRgt(i)
                ! end if
            end if
        end do
    end if

    if (.NOT. stopSim) then
        check_consistency = .true. ! If no inconsistency found, return true
        ! print*, "Left side:", consLft(dim/2), "Right side:", consRgt(dim/2) !debug
    else

    end if
    
end function check_consistency

logical function check_convergence(phiCheck, phiPrev, epsilonCheck)
    implicit none
    real(8), intent(in)  :: phiCheck(:,:), phiPrev(:,:) !uCheck(:,:)
    real(8), intent(in)  :: epsilonCheck
    double precision :: R
    ! real(8), save        :: u_nMinus2 = 0.0d0, u_nMinus1 = 0.0d0, u_n = 0.0d0, h_nMinus2 = 0.0d0, h_nMinus1 = 0.0d0, h_n = 0.0d0
    ! real(8)              :: u_avg, h_avg, u_diff1, u_diff2, h_diff1, h_diff2

    ! u_avg = sum(uCheck) / size(uCheck); h_avg = sum(phiCheck)/size(phiCheck)

    ! ! Shift average history
    ! u_nMinus2 = u_nMinus1; h_nMinus2 = h_nMinus1
    ! u_nMinus1 = u_n;       h_nMinus1 = h_n
    ! u_n     = u_avg;       h_n = h_avg

    ! ! Avoid check on first 2 calls
    ! if (u_nMinus2 == 0.0d0 .and. u_nMinus1 == 0.0d0) then
    !   check_convergence = .false.
    !   return
    ! end if

    ! u_diff1 = u_n       - u_nMinus1;    h_diff1 = h_n       - h_nMinus1
    ! u_diff2 = u_nMinus1 - u_nMinus2;    h_diff2 = h_nMinus1 - h_nMinus2

    ! if (abs(u_diff1 - u_diff2) < epsilonCheck .and. abs(h_diff1 - h_diff2) < epsilonCheck) then
    R = 0
    do x = 1, Lx
        do y = 1, Ly
            R = R + ((phiCheck(x,y) - phiPrev(x,y))/phiCheck(x,y))*((phiCheck(x,y) - phiPrev(x,y))/phiCheck(x,y))
        end do
    end do

    R = dsqrt(R)
    if ( R < epsilonCheck ) then
        print*, "Solution converges."
        check_convergence = .true.
    else
        check_convergence = .false.
    end if
  end function check_convergence

subroutine MMS_analytic_solution
    ! double precision:: phi
    double precision, dimension(3):: phi_0,phi_k,phi_x,phi_y,phi_xy,a_phix,a_phiy,a_phixy,BC_coeff ! MMS constants
    double precision, allocatable, dimension(:,:,:) :: phi,phi_1
    
    allocate (phi(3,Lx,Ly), phi_1(3,Lx,Ly))
    
    ! indices: 1-depth (h); 2-horizontal velocity (u); 3-vertical velocity (v)
    nu_MMS = 1.0d-2
    
    phi_0   = [2.0d0, 0.0d0, 0.0d0]
    phi_k   = [0.0d0, 0.0d0, 0.0d0]
    phi_x   = [0.0d0, 0.0d0, 0.0d0]
    phi_y   = [0.0d0, 0.0d0, 0.0d0]
    phi_xy  = [0.0d0, 0.0d0, 0.0d0]
    a_phix  = [0.75d0, 5.0d0/3.0d0, 1.5d0]
    a_phiy  = [1.0d0, 1.2d0, 1.0d0]
    a_phixy = [1.25d0, 0.2d0, 0.9d0]

    
    do i = 1, Lx
        position_x = dx*i
        do j = 1, Ly
            position_y = dy*j
            
            ! phi_1(:,i,j) = phi_k(:) &
            ! & + phi_x(:)  * DSIN(a_phix(:)  * pi * position_x / domainX) &
            ! & + phi_y(:)  * DSIN(a_phiy(:)  * pi * position_y / domainY) &
            ! & + phi_xy(:) * DSIN(a_phixy(:) * pi * position_x * position_y / (domainX * domainY))
            
            ! BC_coeff(1) = (position_x - 0.0d0)*(position_x - 0.0d0) * (domainX - position_x) ! depth
            ! ! BC_coeff(2) = (position_y - 0.0d0) * (domainY - position_y) ! u-velocity
            ! BC_coeff(2) = 1 ! u-velocity
            ! BC_coeff(3) = 1 ! v-velocity
            
            ! phi(:,i,j) = phi_0(:) + phi_1(:,i,j) * BC_coeff(:)

            ! if (phi(1,i,j) >= 0) then
            !     hAnal(i,j) = phi(1,i,j)
            ! else
            !     print *, "Error: Analytic solution has negative depth in i =", i, " j =", j
            !     return
            ! end if
            ! uAnal(i,j) = phi(2,i,j)
            ! vAnal(i,j) = phi(3,i,j)

            ! directly from Sympy code
            hAnal(i,j) = (2.0d0/3.0d0)*sin(6.2831853071795865d0*position_x/domainX) + 2
            if (hAnal(i,j) < 0) then 
                print *, "Error: Analytic solution has negative depth in i =", i, " j =", j
                stopSim = .TRUE.
                exit 
            end if
            uAnal(i,j) = 0.0d0
            vAnal(i,j) = 0.0d0
        end do
    end do
end subroutine MMS_analytic_solution

end module Mu_LaB_SWE


