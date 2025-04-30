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
    integer:: time,itera_no!, i, j commented because only used to debug
    double precision :: ho, uo, vo, simTime
    character:: fdate*24, td*24 ! get date for output
    
    ! Total iteration numbers 
    itera_no = 1
    ! itera_no = 200

    ! constants for initializing flow field. 
    ho = 2.0d0
    uo = 0.0d0
    vo = 0.0d0

    ! define total lattice numbers in x and y directions
    Lx = 250; Ly = 5

    ! assign a value for the relaxation time
    tau = 1.5d0 ! peut-être tendre plus vers 1

    ! allocate dimensions for dynamic arrays
    allocate (f(9,Lx,Ly),feq(9,Lx,Ly),ftemp(9,Lx,Ly),h(Lx,Ly),& 
        & force_x(Lx,Ly),force_y(Lx,Ly),u(Lx,Ly),v(Lx,Ly)) 
    
    ! initialize the depth and velocities 
    h = ho 
    u = uo 
    v = vo

    ! Set constant force
    force_x = 2.4d-8 ! modified from book example
    ! force_x = 2.4d-6 ! original book example
    force_y = 0.0d0
    print*, itera_no
    print*, ho
    print*, uo
    print*, vo
    ! print*, force_x, force_y
    ! prepare the calculations
    call setup

    ! main loop for time marching
    time = 0
    timStep: do

        time = time+1

        ! Streaming and collision steps
        call collide_stream

        ! Apply Noslip BC
        call Noslip_BC

        ! Calculate h, u & v
        call solution

        ! Update the feq
        call compute_feq

        write(6,'(I6,A2,3(ES26.16,A2))') time,'   ', u(Lx/2, Ly/2)

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
        
        if (time == itera_no) exit

    end do timStep

    write(6,*) 
    write(6,*)' Writing results in file: result.dat ... ' 
    open(66,file='result.dat',status='unknown') 
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