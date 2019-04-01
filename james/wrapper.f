!**********************************************************************!
!                          Charlie Anderson                            !
!                  University of Bristol, March 2019                   !
!**********************************************************************!
!    
      program wrapper
!    
      integer :: nx, ny, nl
!    
      integer :: na = 2
      integer :: np = 5
      double precision, dimension(1,2) :: alpha
      double precision, dimension(1,5) :: params
!    
!    
!    
!    
!    
!     mesh information, numbers define x and y direction and number of
!     Halo cells 
!    
      nx = 100
      ny = 40
      nl = 2

!      
!    
!     design variables: currently only throat ratio but will increase at
!     some point  
!    
      alpha(1,1) = 2.0d0
      alpha(1,2) = 5.0d0
!    
!    
!    
!     information to run the CFD code, imput as a vector for simplicity
!     Mach_in   Pres_rat   gam   CFL    tmax
!    
      params(1,1) = 0.3d0
      params(1,2) = 0.8d0
      params(1,3) = 1.4d0
      params(1,4) = 0.2d0
      params(1,5) = DBLE(300)
!    
!    
!    
!    
      call main(nx,ny,na,nl,np,alpha,params)
!    
!    
!    
!    
      end program wrapper
