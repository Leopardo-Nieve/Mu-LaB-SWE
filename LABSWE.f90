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

        integer:: Lx,Ly,x,y,a!,b !b for debugging
        double precision:: q_in,h_out,dx,dy,domainX,domainY,dt,eMin,e,eZhou,tau,tauZhou,nu,nuZhou,qZhou,ReZhou,gacl = 9.81 
        double precision, dimension(9):: ex,ey 
        double precision, allocatable, dimension(:,:):: u,v,h,force_x,force_y
        double precision, allocatable, dimension(:,:,:):: f,feq,ftemp 
    
contains 

subroutine setup 
    ! declare a local double precision for the quarter of PI 
    double precision:: quarter_pi 

    ! set constant PI 
    quarter_pi = datan(1.0d0) 

    ! compute the particle velocities 
    do a = 1, 8 
        if (mod(a,2) == 0) then 
            ex(a) = dsqrt(2.0d0)*dcos(quarter_pi*dble(a-1)) 
            ey(a) = dsqrt(2.0d0)*dsin(quarter_pi*dble(a-1)) 
        else
            ex(a) = dcos(quarter_pi*dble(a-1)) 
            ey(a) = dsin(quarter_pi*dble(a-1))

        end if 
    end do
    ex(9) = 0.0d0; ey(9) = 0.0d0

    ! print*, "particle velocities defined" !debugging
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
            if (xf > Lx) xf = xf - Lx
            if (xb < 1) xb = Lx + xb 
            if (yf > Ly) yf = yf - Ly
            if (yb < 1) yb = Ly + yb 
            
            ! Following lines implement inflow BC (Zhou, p.59)
            ! ftemp(1,0,:) = f(5,0,:) + 2*q_in/(3*e)
            ! ftemp(2,0,:) = q_in/(6*e) + f(6,0,:) + 0.5*(f(7,0,:) - f(3,0,:))
            ! ftemp(8,0,:) = q_in/(6*e) + f(4,0,:) + 0.5*(f(3,0,:) - f(7,0,:))

            ! Following lines implement outflow BC (Zhou, p.60) and Neumann  
            ! h(Lx,:) = h_out
            ! ftemp(5,Lx,:) = f(1,Lx,:) - 2*h(Lx,:)*u(Lx-1,:)/(3*e)
            ! ftemp(4,Lx,:) = -h(Lx,:)*u(Lx-1,:)/(6*e) + f(8,Lx,:) + 0.5*(f(7,Lx,:) - f(3,Lx,:))
            ! ftemp(6,Lx,:) = -h(Lx,:)*u(Lx-1,:)/(6*e) + f(2,Lx,:) + 0.5*(f(3,Lx,:) - f(7,Lx,:))

            ! start streaming and collision 
            ftemp(1,xf,y) = f(1,x,y)-(f(1,x,y)-feq(1,x,y))/tau&
                & + 1.0d0/6.0d0*(ex(1)*force_x(x,y)+ey(1)*force_y(x,y))
            ftemp(2,xf,yf) = f(2,x,y)-(f(2,x,y)-feq(2,x,y))/tau& 
                & + 1.0d0/6.0d0*(ex(2)*force_x(x,y)+ey(2)*force_y(x,y)) 
            ftemp(3,x,yf) = f(3,x,y)-(f(3,x,y)-feq(3,x,y))/tau& 
                & + 1.0d0/6.0d0*(ex(3)*force_x(x,y)+ey(3)*force_y(x,y)) 
            ftemp(4,xb,yf) = f(4,x,y)-(f(4,x,y)-feq(4,x,y))/tau& 
                & + 1.0d0/6.0d0*(ex(4)*force_x(x,y)+ey(4)*force_y(x,y)) 
            ftemp(5,xb,y) = f(5,x,y)-(f(5,x,y)-feq(5,x,y))/tau& 
                & + 1.0d0/6.0d0*(ex(5)*force_x(x,y)+ey(5)*force_y(x,y)) 
            ftemp(6,xb,yb) = f(6,x,y)-(f(6,x,y)-feq(6,x,y))/tau& 
                & + 1.0d0/6.0d0*(ex(6)*force_x(x,y)+ey(6)*force_y(x,y)) 
            ftemp(7,x,yb) = f(7,x,y)-(f(7,x,y)-feq(7,x,y))/tau& 
                & + 1.0d0/6.0d0*(ex(7)*force_x(x,y)+ey(7)*force_y(x,y)) 
            ftemp(8,xf,yb) = f(8,x,y)-(f(8,x,y)-feq(8,x,y))/tau& 
                & + 1.0d0/6.0d0*(ex(8)*force_x(x,y)+ey(8)*force_y(x,y)) 
            ftemp(9,x,y) = f(9,x,y) - (f(9,x,y)-feq(9,x,y))/tau 
            
            ! debug to determine which populations are negative
            ! do a=1,9
            !     if (ftemp(a,Lx/2,Ly/2) < 0) print*, "f", a, "=", ftemp(a,Lx/2,Ly/2)
            ! end do
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

    ! print*, "beginning equilibrium populations" !debugging
    do a = 1, 8 
        ! if (mod(a,2) == 0) then 
        feq(a,:,:) = gacl*h(:,:)*h(:,:)/24.0d0 +& 
            & h(:,:)/12.0d0*(ex(a)*u(: ,:)+& 
            & ey(a)*v(:,:))+h(:,:)/8.0d0& 
            & *(ex(a)*u(:,:)*ex(a)*u(:,:)+& 
            & 2.0d0*ex(a)*u(:,:)*ey(a)*v(:,:)+& 
            & ey(a)*v(:,:)*ey(a)*v(:,:))-& 
            & h(:,:)/24.0d0*(u(:,:)*u(:,:)+& 
            & V ( : , : ) *V ( : , : ) ) 
        ! end if

        if (mod(a,2) /= 0) feq(a,:,:) = 4.0d0*feq(a,:,:) ! if odd number index
        print*, "feq", a, "=", feq(a,Lx/2,Ly/2) !debugging
    end do
    feq(9,:,:) = h(:,:) - 5.0d0*gacl*h(:,:)*h(:,:)/(6.0d0*e**2) - &
             & 2.0d0*h(:,:)/(3.0d0*e**2)*(u(:,:)**2 + v(:,:)**2)
    
    print*, "feq", a, "=", feq(a,Lx/2,Ly/2) !debugging
    ! if (feq(a,Lx/2,Ly/2) < 0) then !debugging
    !     print*, "feq", a, "=", feq(a,Lx/2,Ly/2) !debugging
    !     if (a == 9) then !debugging
    !         print*, "term 1 =", h(Lx/2,Ly/2) !debugging
    !         print*, "term 2 =", - 5.0d0*gacl*h(Lx/2,Ly/2)*h(Lx/2,Ly/2)/(6.0d0*e**2) !debugging
    !         print*, "term 3 =", - 2.0d0*h(Lx/2,Ly/2)/(3.0d0*e**2)*(u(Lx/2,Ly/2)**2 + v(Lx/2,Ly/2)**2.0d0) !debugging
    !     end if !debugging
    ! end if
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
    ftemp(2,:,1) = f(8,:,1) 
    ftemp(3,:,1) = f(7,:,1) 
    ftemp(4,:,1) = f(6,:,1) 
    
    ! for upper boundary
    ftemp(8,:,Ly) = f(2,:,Ly) 
    ftemp(7,:,Ly) = f(3,:,Ly)
    ftemp(6,:,Ly) = f(4,:,Ly)

    return 
end subroutine Slip_BC 

subroutine write_csv
    ! Write simulation results to a CSV file for ParaView visualization
    open(67, file='result.csv', status='unknown')
    
    ! Write CSV header
    write(67, '(A)') 'x,y,h,u,v'
    
    ! Write data points
    do x = 1, Lx
        do y = 1, Ly
            write(67, '(I0,",",I0,3(",",ES15.6))') x, y, h(x,y), u(x,y), v(x,y)
        end do
    end do
    
    close(67)
end subroutine write_csv

end module LABSWE