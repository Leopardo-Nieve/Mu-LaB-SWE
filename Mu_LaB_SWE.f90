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
! a, x, y, b, i, j - Loop integers
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
! h_out - fixed depth at outflow [m]
! H_part - partial depth [m]
! Lx, Ly - Total lattice numbers in x and y directions [-]
! nu - Molecular viscosity  [m^2/s]
! q_in - Inflow discharge [m^2/s]
! tau - Relaxation time [-]
! time - Amount of time elapsed since start of simulation [s]
! u, v - x and y components of flow velocity [m/s]
! u_out - fixed velocity at outflow [m/s]
! zb - bed geometry
module Mu_LaB_SWE

        implicit none 

        integer:: Lx,Ly,x,y,a,current_iteration, b,i,j
        integer, dimension(2):: hIndex
        logical:: stopSim, tauOk, velOk, celOk, FrOk
        double precision:: q_in,h_out,u_out, dx,dy,domainX,domainY,time,dt,eMin,e,tau,nu,&
        &dt_6e2,one_8th_e4,one_3rd_e2,one_6th_e2,one_12th_e2, one_24th_e2,five_6th_g_e2,two_3rd_e2,gacl = 9.81,&
        & hMax, uMax2, FrMax, Fr, Ma, consCriter,pi 
        double precision, dimension(9):: ex,ey, eMax
        double precision, allocatable, dimension(:,:):: u,v,h,force_x,force_y,H_part,zb,dzbdx,consInLft,consInRgt,consOutLft,consOutRgt
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

    ! compute the equilibrium distribution function feq 
    call compute_feq

    ! Set the initial distribution function to feq 
    f = feq
    return
end subroutine setup

subroutine collide_stream

    ! This calculates distribution function with the LABSWE 

    ! local working integers 
    integer:: xf,yf,xb,yb

    do y = 1, Ly 
        yf = y + 1
        yb = y -1
    
        do x = 1, Lx
            xf = x + 1
            xb = x -1

            ! Following 2 lines Implement periodic BCs in y-direction 
            ! if (xf > Lx) xf = xf - Lx !remove outlet periodic boundary
            ! if (xb < 1) xb = Lx + xb !remove inlet periodic boundary
            if (yf > Ly) yf = yf - Ly
            if (yb < 1) yb = Ly + yb 

            ! start streaming and collision 
            if (xf<=Lx) then ! periodic in y direction
                ftemp(1,xf,y) = f(1,x,y)-(f(1,x,y)-feq(1,x,y))/tau&
                & + dt_6e2*(ex(1)*0.5d0*(force_x(x,y)+force_x(xf,y))+ey(1)*0.5d0*(force_y(x,y)+force_y(xf,y)))
            end if
            if (xf<=Lx) then !if (xf<=Lx .and. yf<=Ly) ! periodic in y direction
            ftemp(2,xf,yf) = f(2,x,y)-(f(2,x,y)-feq(2,x,y))/tau& 
                & + dt_6e2*(ex(2)*0.5d0*(force_x(x,y)+force_x(xf,yf))+ey(2)*0.5d0*(force_y(x,y)+force_y(xf,yf)))
            end if
            ! if (yf<=Ly) ! periodic in y direction
            ftemp(3,x,yf) = f(3,x,y)-(f(3,x,y)-feq(3,x,y))/tau& 
                & + dt_6e2*(ex(3)*0.5d0*(force_x(x,y)+force_x(x,yf))+ey(3)*0.5d0*(force_y(x,y)+force_y(x,yf)))
            if (xb>=1) then !if (xb>=1 .and. yf<=Ly) ! periodic in y direction
                ftemp(4,xb,yf) = f(4,x,y)-(f(4,x,y)-feq(4,x,y))/tau& 
                & + dt_6e2*(ex(4)*0.5d0*(force_x(x,y)+force_x(xb,yf))+ey(4)*0.5d0*(force_y(x,y)+force_y(xb,yf)))
            end if
            if (xb>=1) ftemp(5,xb,y) = f(5,x,y)-(f(5,x,y)-feq(5,x,y))/tau& 
                & + dt_6e2*(ex(5)*0.5d0*(force_x(x,y)+force_x(xb,y))+ey(5)*0.5d0*(force_y(x,y)+force_y(xb,y)))
            if (xb>=1) then !if (xb>=1 .and. yb>=1) ! periodic in y direction
                ftemp(6,xb,yb) = f(6,x,y)-(f(6,x,y)-feq(6,x,y))/tau& 
                & + dt_6e2*(ex(6)*0.5d0*(force_x(x,y)+force_x(xb,yb))+ey(6)*0.5d0*(force_y(x,y)+force_y(xb,yb)))
            end if
            ! if (yb>=1) ! periodic in y direction
            ftemp(7,x,yb) = f(7,x,y)-(f(7,x,y)-feq(7,x,y))/tau& 
                & + dt_6e2*(ex(7)*0.5d0*(force_x(x,y)+force_x(x,yb))+ey(7)*0.5d0*(force_y(x,y)+force_y(x,yb)))
            if (xf<=Lx) then !if (xf<=Lx .and. yb>=1) ! periodic in y direction
                ftemp(8,xf,yb) = f(8,x,y)-(f(8,x,y)-feq(8,x,y))/tau& 
                & + dt_6e2*(ex(8)*0.5d0*(force_x(x,y)+force_x(xf,yb))+ey(8)*0.5d0*(force_y(x,y)+force_y(xf,yb)))
            end if
            ftemp(9,x,y) = f(9,x,y) - (f(9,x,y)-feq(9,x,y))/tau 
            
        end do 
    end do

return
end subroutine collide_stream

subroutine solution
    
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
    ! macroscopic inflow values
    h(1,:) = h_in(time)

    ! consistence check
    consInLft(1,:) = h(1,:)-ftemp(9,1,:) ! left side of the consistence equation
    do a = 3, 7
        consInLft(1,:) = consInLft(1,:) - ftemp(a,1,:)
    end do
    consInRgt(1,:) = h(1,:)*u(1,:)/e + ftemp(4,1,:) + ftemp(5,1,:) + ftemp(6,1,:) ! right side of the consistence equation
    do j = 1, Ly
        if ( abs(consInLft(1,j) - consInRgt(1,j)) > consCriter ) then
            print*, "consistency fails at node",1,j
            print*,consInLft(1,j),"/=", consInRgt(1,j)
            stopSim = .true.
        end if
    end do
    
    if ( .not. stopSim ) then
        ! Following lines implement inflow BC (Zhou, p.59)
        ftemp(1,1,:) = ftemp(5,1,:) + 2.0d0*q_in/(3.0d0*e)
        ftemp(2,1,:) = q_in/(6.0d0*e) + ftemp(6,1,:) + 0.5d0*(ftemp(7,1,:) - ftemp(3,1,:))
        ftemp(8,1,:) = q_in/(6.0d0*e) + ftemp(4,1,:) + 0.5d0*(ftemp(3,1,:) - ftemp(7,1,:))    
    end if

    ! macroscopic outflow values
    u_out = 0.0d0
    do a = 1, 9
        u_out = u_out + ex(a)*ftemp(a,Lx-1,Ly/2)
    end do
    u_out = u_out/h_out
    
    ! consistence check
    consOutLft(1,:) = h_out
    do a = 1, 9
        if (a >= 4 .and. a <=6 ) cycle
        consOutLft(1,:) = consOutLft(1,:) - ftemp(a,Lx,:)
    end do
    consOutRgt(1,:) = -h_out*u_out/e + ftemp(1,Lx,:) + ftemp(2,Lx,:) + ftemp(8,Lx,:)
    do j = 1, Ly
        if ( abs(consOutLft(1,j) - consOutRgt(1,j)) > consCriter ) then
            print*, "consistence fails at node",Lx,j
            print*,consOutLft(1,j),"/=", consOutRgt(1,j)
            stopSim = .true.
        end if
    end do
    if ( .not. stopSim ) then
        ! Following lines implement outflow BC (Zhou, p.60) and Neumann  
        ftemp(5,Lx,:) = ftemp(1,Lx,:) - 2.0d0*h_out*u_out/(3.0d0*e)
        ftemp(4,Lx,:) = -h_out*u_out/(6.0d0*e) + ftemp(8,Lx,:) + 0.5d0*(ftemp(7,Lx,:) - ftemp(3,Lx,:))
        ftemp(6,Lx,:) = -h_out*u_out/(6.0d0*e) + ftemp(2,Lx,:) + 0.5d0*(ftemp(3,Lx,:) - ftemp(7,Lx,:))    
    end if
end subroutine Inflow_Outflow_BC


subroutine write_csv
    ! Write simulation results to a CSV file for ParaView visualization
    open(67, file='../results_7.2.1/result.csv', status='unknown')
    
    ! write simulation parameters
    write(67, '(A)') 'Date,Iteration No.,tau,ujuj/e^2,gh/e^2,Fr'
    write(67, '(A,",",I5,",",F10.5,",",F10.5,",",F10.5,",",F10.5)') &
     trim(fdate()), current_iteration, tau, uMax2/(e*e), gacl*hMax/(e*e), FrMax


    
    ! Write CSV header
    write(67, '(A)') 'x (nodes),y (nodes),x (m),y (m),h + zb (m),zb (m),h (m),u (m/s),v (m/s), q (m^2/s)'

    ! Write data points
    do x = 1, Lx
        do y = 1, Ly
            write(67,'(2(I5,","),2(F12.4,","),3(F12.4,","),3(F12.4,","),F12.4)') &
                x, y, &
                x * dx, y * dy, &
                h(x,y) + zb(x,y), zb(x,y), h(x,y), &
                u(x,y), v(x,y), h(x,y)*u(x,y)
        end do
    end do
    
    close(67)
end subroutine write_csv

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

double precision function h_in(currentTime)
    implicit none
    double precision, intent(in)    :: currentTime
    pi = dacos(-1.0d0)

    h_in = H_part(0,Ly/2) + 4.0d0 - 4.0d0*dsin(pi*(4.0d0*currentTime/86.4d3 + 0.5d0))
end function h_in

logical function check_convergence(uCheck, hCheck, epsilonCheck)
    implicit none
    real(8), intent(in)  :: uCheck(:,:), hCheck(:,:)
    real(8), intent(in)  :: epsilonCheck
    real(8), save        :: u_nMinus2 = 0.0d0, u_nMinus1 = 0.0d0, u_n = 0.0d0, h_nMinus2 = 0.0d0, h_nMinus1 = 0.0d0, h_n = 0.0d0
    real(8)              :: u_avg, h_avg, u_diff1, u_diff2, h_diff1, h_diff2

    u_avg = sum(uCheck) / size(uCheck); h_avg = sum(hCheck)/size(hCheck)

    ! Shift average history
    u_nMinus2 = u_nMinus1; h_nMinus2 = h_nMinus1
    u_nMinus1 = u_n;       h_nMinus1 = h_n
    u_n     = u_avg;       h_n = h_avg

    ! Avoid check on first 2 calls
    if (u_nMinus2 == 0.0d0 .and. u_nMinus1 == 0.0d0) then
      check_convergence = .false.
      return
    end if

    u_diff1 = u_n       - u_nMinus1;    h_diff1 = h_n       - h_nMinus1
    u_diff2 = u_nMinus1 - u_nMinus2;    h_diff2 = h_nMinus1 - h_nMinus2

    if (abs(u_diff1 - u_diff2) < epsilonCheck .and. abs(h_diff1 - h_diff2) < epsilonCheck) then
        print*, "Solution converges."
        check_convergence = .true.
    else
        check_convergence = .false.
    end if
  end function check_convergence

end module Mu_LaB_SWE