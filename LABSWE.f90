!----------------------------------------------------------! 
!                       LABSWE.f90 
! This module written in FORTRAN 90 implements the lattice 
! Boltzmann method for shallow water equations (LABSWE), 
! based on the 9-speed square lattices. It is included as 
! a sample code in the book and hence only highlights the 
! major procedure in the LABSWE under conditions that time 
! step and lattice spacing are taken as units. Also the 
! module provides periodic boundary conditions and no-slip 
! boundary conditions in y direction boundaries. It should 
! be pointed that the module can easily be adapted into a 
! practical code by changing these settings or adding more 
! boundary conditions. A sample file main.f90 is also 
! included to show how to use the module. 
! J.G. Zhou, Peterborough, 2003 
! ----------------------------------------------------------!
!                   List of Major Variables
! a, x, y - Loop integers
! ex, ey - x and y components of particles' velocities
! f  - Distribution function
! feq  - local equilibrium distribution function
! force_x - x-direction component of force term
! force_y - y-direction component of force term
! ftemp - Temple distribution function
! gael - Gravitational acceleration
! h - water depth
! Lx, Ly - Total lattice numbers in x and y directions
! nu - Molecular viscosity  
! tau - Relaxation time
! u, v - x and y components of flow velocity 
module LABSWE 

        implicit none 

        integer:: Lx,Ly,x,y,a,b,i,j !b,i,j for debug
        integer, dimension(2):: hIndex
        logical:: stopSim, tauOk, velOk, celOk, FrOk
        double precision:: q_in,h_out,dx,dy,domainX,domainY,dt,eMin,e,eZhou,tau,tauZhou,nu,nuZhou,qZhou,ReZhou,&
        &dt_6e2,one_8th_e4,one_3rd_e2,one_6th_e2,one_12th_e2, one_24th_e2,five_6th_g_e2,two_3rd_e2,gacl = 9.81,&
        & hMax, uMax2, FrMax, Fr, Ma 
        double precision, dimension(9):: ex,ey, eMax
        double precision, allocatable, dimension(:,:):: u,v,h,force_x,force_y,zb,dzbdx
        double precision, allocatable, dimension(:,:,:):: f,feq,ftemp 
    
contains 

subroutine setup 
    ! declare a local double precision for the quarter of PI 
    ! double precision:: quarter_pi 

    ! ! set constant PI 
    ! quarter_pi = datan(1.0d0) 

    ! ! compute the particle velocities 
    ! do a = 1, 8 
    !     if (mod(a,2) == 0) then 
    !         ex(a) = dsqrt(2.0d0)*dcos(quarter_pi*dble(a-1)) 
    !         ey(a) = dsqrt(2.0d0)*dsin(quarter_pi*dble(a-1)) 
    !     else
    !         ex(a) = dcos(quarter_pi*dble(a-1)) 
    !         ey(a) = dsin(quarter_pi*dble(a-1))

    !     end if 
    ! end do
     ! hard-code lattice velocities    

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
    print*, "eMax = ", eMax


    ! constants to limit random error
    one_24th_e2=1.0d0/(24.0d0*e*e)
    one_12th_e2=2.0d0*one_24th_e2
    one_6th_e2 =2.0d0*one_12th_e2
    one_3rd_e2 =2.0d0*one_6th_e2
    one_8th_e4 = 1.0d0/(8.0d0*e*e*e*e)
    five_6th_g_e2 = 5.0d0*gacl*one_6th_e2
    two_3rd_e2 = 2.0d0*one_3rd_e2

    dt_6e2=dt/(6.0d0*e*e)


    ! print*, "particle velocities defined" !debug
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

            ! print*, x, y !debug
            ! Following 4 lines Implement periodic BCs 
            ! if (xf > Lx) xf = xf - Lx
            ! if (xb < 1) xb = Lx + xb 
            if (yf > Ly) yf = yf - Ly
            if (yb < 1) yb = Ly + yb 

            ! start streaming and collision 
            if (xf<=Lx) then
                ftemp(1,xf,y) = f(1,x,y)-(f(1,x,y)-feq(1,x,y))/tau&
                & + dt_6e2*(ex(1)*0.05d0*(force_x(x,y)+force_x(xf,y))+ey(1)*0.05d0*(force_y(x,y)+force_y(xf,y)))
            end if
            if (xf<=Lx) then !if (xf<=Lx .and. yf<=Ly) 
            ftemp(2,xf,yf) = f(2,x,y)-(f(2,x,y)-feq(2,x,y))/tau& 
                & + dt_6e2*(ex(2)*0.05d0*(force_x(x,y)+force_x(xf,yf))+ey(2)*0.05d0*(force_y(x,y)+force_y(xf,yf)))
            end if
            ! if (yf<=Ly) 
            ftemp(3,x,yf) = f(3,x,y)-(f(3,x,y)-feq(3,x,y))/tau& 
                & + dt_6e2*(ex(3)*0.05d0*(force_x(x,y)+force_x(x,yf))+ey(3)*0.05d0*(force_y(x,y)+force_y(x,yf)))
            if (xb>=1) then !if (xb>=1 .and. yf<=Ly) 
                ftemp(4,xb,yf) = f(4,x,y)-(f(4,x,y)-feq(4,x,y))/tau& 
                & + dt_6e2*(ex(4)*0.05d0*(force_x(x,y)+force_x(xb,yf))+ey(4)*0.05d0*(force_y(x,y)+force_y(xb,yf)))
            end if
            if (xb>=1) ftemp(5,xb,y) = f(5,x,y)-(f(5,x,y)-feq(5,x,y))/tau& 
                & + dt_6e2*(ex(5)*0.05d0*(force_x(x,y)+force_x(xb,y))+ey(5)*0.05d0*(force_y(x,y)+force_y(xb,y)))
            if (xb>=1) then !if (xb>=1 .and. yb>=1) 
                ftemp(6,xb,yb) = f(6,x,y)-(f(6,x,y)-feq(6,x,y))/tau& 
                & + dt_6e2*(ex(6)*0.05d0*(force_x(x,y)+force_x(xb,yb))+ey(6)*0.05d0*(force_y(x,y)+force_y(xb,yb)))
            end if
            ! if (yb>=1) 
            ftemp(7,x,yb) = f(7,x,y)-(f(7,x,y)-feq(7,x,y))/tau& 
                & + dt_6e2*(ex(7)*0.05d0*(force_x(x,y)+force_x(x,yb))+ey(7)*0.05d0*(force_y(x,y)+force_y(x,yb)))
            if (xf<=Lx) then !if (xf<=Lx .and. yb>=1) 
                ftemp(8,xf,yb) = f(8,x,y)-(f(8,x,y)-feq(8,x,y))/tau& 
                & + dt_6e2*(ex(8)*0.05d0*(force_x(x,y)+force_x(xf,yb))+ey(8)*0.05d0*(force_y(x,y)+force_y(xf,yb)))
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

    ! print*, "beginning equilibrium populations" !debug
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
    ! do a=1,9
    !     ! i=1+1
    !     ! if (a==4) then
    !     !     j=3-1
    !     ! elseif(a==5) then
    !     !     j=3
    !     ! else
    !     !     j=3+1
    !     ! end if
    !     i=1
    !     j=3
    !     print*, "Population:",a
    !     print*, "ex(",a,") =",ex(a)
    !     print*, "1st term", gacl*h(i,j)*h(i,j)*one_24th_e2
    !     print*, "2nd term", h(i,j)*one_12th_e2*(ex(a)*u(i,j)+ ey(a)*v(i,j))
    !     print*, "3rd term",h(i,j)*one_8th_e4& 
    !         & *(ex(a)*u(i,j)*ex(a)*u(i,j)+& 
    !         & 2.0d0*ex(a)*u(i,j)*ey(a)*v(i,j)+& 
    !         & ey(a)*v(i,j)*ey(a)*v(i,j))
    !     print*, "4th term",- h(i,j)*one_24th_e2*(u(i,j)*u(i,j)+ v(i,j)*v(i,j)) 
    !     print*,"feq",a,feq(a,i,j)
    !     print*, "Please press enter to continue"
    ! read(*,*)
    ! end do

    ! do a=1,9
    !     print*, "f(",a,",1,3) =",f(a,1,3)
    ! end do
    ! read(*,*)

             ! print*, "feq_0 = ",feq(9,1,1) !debug
    ! if (feq(a,Lx/2,Ly/2) < -1.0d-23) then !debug
    !     print*, "feq", a, "=", feq(a,Lx/2,Ly/2) !debug
    !     if (a == 9) then !debug
    !         print*, "term 1 =", h(Lx/2,Ly/2) !debug
    !         print*, "term 2 =", - 5.0d0*gacl*h(Lx/2,Ly/2)*h(Lx/2,Ly/2)/(6.0d0*e**2) !debug
    !         print*, "term 3 =", - 2.0d0*h(Lx/2,Ly/2)/(3.0d0*e**2)*(u(Lx/2,Ly/2)**2 + v(Lx/2,Ly/2)**2.0d0) !debug
    !     end if !debug
    ! end if !debug
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
    ! Following lines implement inflow BC (Zhou, p.59)
    ftemp(1,1,:) = ftemp(5,1,:) + 2.0d0*q_in/(3.0d0*e)
    ftemp(2,1,:) = q_in/(6.0d0*e) + ftemp(6,1,:) + 0.5d0*(ftemp(7,1,:) - ftemp(3,1,:))
    ftemp(8,1,:) = q_in/(6.0d0*e) + ftemp(4,1,:) + 0.5d0*(ftemp(3,1,:) - ftemp(7,1,:))


    ! do a=1,9 !debug
    !     print*, "f after inflow", a, "=", ftemp(a,1,1)
    ! end do ! debug

    ! Following lines implement outflow BC (Zhou, p.60) and Neumann  
    ! h(Lx,:) = h_out
    ftemp(5,Lx,:) = ftemp(1,Lx,:) - 2.0d0*h(Lx,:)*u(Lx-1,:)/(3.0d0*e)
    ftemp(4,Lx,:) = -h(Lx,:)*u(Lx-1,:)/(6.0d0*e) + ftemp(8,Lx,:) + 0.5d0*(ftemp(7,Lx,:) - ftemp(3,Lx,:))
    ftemp(6,Lx,:) = -h(Lx,:)*u(Lx-1,:)/(6.0d0*e) + ftemp(2,Lx,:) + 0.5d0*(ftemp(3,Lx,:) - ftemp(7,Lx,:))

    ! do a=1,9 !debug
    !     print*, "f final", a, "=", ftemp(a,1,1)
    ! end do ! debug
end subroutine Inflow_Outflow_BC

subroutine Four_Corners_BC
    ! node (1,1)
    ! ftemp(3,1,1) = ftemp(7,1,1) !slip (should be already done)
    ! ftemp(4,1,1) = ftemp(6,1,1) !slip (should be already done)
    ! ftemp(2,1,1) = q_in/(6*e) + ftemp(6,1,1) + 0.5*(ftemp(7,1,1) - ftemp(3,1,1)) !inflow
    ftemp(8,1,1) = ftemp(2,1,1) !slip

    ! node (1,Ly)
    ! ftemp(6,1,Ly) = ftemp(4,1,Ly) !slip (should be already done)
    ! ftemp(7,1,Ly) = ftemp(3,1,Ly) !slip (should be already done)
    ! ftemp(8,1,Ly) = q_in/(6*e) + ftemp(4,1,Ly) + 0.5*(ftemp(3,1,Ly) - ftemp(7,1,Ly)) !inflow
    ftemp(2,1,Ly) = ftemp(8,1,Ly) !slip

    ! node (Lx,1)
    ! ftemp(2,1,1) = ftemp(8,1,1) !slip (should be already done)
    ! ftemp(3,1,1) = ftemp(7,1,1) !slip (should be already done)
    ! ftemp(4,Lx,1) = -h(Lx,1)*u(Lx-1,1)/(6*e) + ftemp(8,Lx,1) + 0.5*(ftemp(7,Lx,1) - ftemp(3,Lx,1))
    ftemp(6,Lx,1) = ftemp(4,Lx,1)

    ! node (Lx,Ly)
    ftemp(4,Lx,Ly) = ftemp(6,Lx,Ly)

end subroutine Four_Corners_BC


subroutine write_csv
    ! Write simulation results to a CSV file for ParaView visualization
    open(67, file='../results_7.2.1/result.csv', status='unknown')
    
    ! Write CSV header
    write(67, '(A)') 'x (nodes),y (nodes),x (m),y (m),h + zb (m),zb (m),h (m),u (m/s),v (m/s)'

    ! Write data points
    do x = 1, Lx
        do y = 1, Ly
            write(67,'(2(I5,","),2(F12.4,","),3(F12.4,","),2(F12.4,","),F12.4)') &
                x, y, &
                x * dx, y * dy, &
                h(x,y) + zb(x,y), zb(x,y), h(x,y), &
                u(x,y), v(x,y)
        end do
    end do
    
    close(67)
end subroutine write_csv

logical function check_convergence(uCheck, epsilonCheck)
    implicit none
    real(8), intent(in)  :: uCheck(:,:)
    real(8), intent(in)  :: epsilonCheck
    real(8), save        :: u_nMinus2 = 0.0d0, u_nMinus1 = 0.0d0, u_n = 0.0d0
    real(8)              :: current_avg, diff1, diff2

    current_avg = sum(uCheck) / size(uCheck)

    ! Shift average history
    u_nMinus2 = u_nMinus1
    u_nMinus1 = u_n
    u_n = current_avg

    ! Avoid check on first 2 calls
    if (u_nMinus2 == 0.0d0 .and. u_nMinus1 == 0.0d0) then
      check_convergence = .false.
      return
    end if

    diff1 = u_n - u_nMinus1
    diff2 = u_nMinus1 - u_nMinus2

    if (abs(diff1 - diff2) < epsilonCheck) then
      check_convergence = .true.
    else
      check_convergence = .false.
    end if
  end function check_convergence

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

end module LABSWE