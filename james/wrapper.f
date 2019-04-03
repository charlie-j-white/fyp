!**********************************************************************!
!                          Charlie Anderson                            !
!                  University of Bristol, March 2019                   !
!**********************************************************************!
!    
      program wrapper
!    
      integer :: na = 40
      integer :: np = 12
      double precision, dimension(1,40) :: alpha
      double precision, dimension(1,12) :: params
!    
      integer :: nx, ny, nl
      double precision :: AR, pi, frac
!    
!    
!    
!    
!    
!     information about the structured mesh
!-----------------------------------------------------------------------
!    
!     x cells      y cells      halo cells (must be 2 for JST)
!    
      nx = 100
      ny = 40
      nl = 2

!    
!    
!      
!    
!    
!      
!    
!    
!     information to run the CFD code, imput as a vector for simplicity
!-----------------------------------------------------------------------
!    
!    
!     Mach_in      Pres_rat       gam       CFL      tmax   
!    
      params(1,1) = 0.3d0
      params(1,2) = 0.8d0
      params(1,3) = 1.4d0
      params(1,4) = 0.3d0
      params(1,5) = DBLE(10000)
!    
!    
!     K2_jst      K4_jst
!    
      params(1,6) = 1.0d0/16.0d0
      params(1,7) = 1.0d0/32.0d0
!    
!    
!     adaptive meshing ON/OFF (1/0)     
!    
      params(1,8) = 1.0d0
!    
!    
!     s_pos       s_hgt         s_wdt      Ly
!    
      params(1,9) = 0.7d0
      params(1,10) = 5.0d0
      params(1,11) = 0.05d0
      params(1,12) = 0.4d0
!    
!    
!      
!    
!    
!      
!      
!    
!    
!      
!    
!     design variables: currently only throat ratio but will increase at
!     some point  
!-----------------------------------------------------------------------
!    
      AR = 2.0d0
      pi = 3. 1415926535897932d0
!
!     set the height of the mesh along the length of the duct
      do i = 1,na
      frac = DBLE(i)/(na+1.0d0)
      alpha(1,i) = 1.0d0/AR + (1.0d0 - 1.0d0/AR)*cos(pi*frac)**2.0d0
!      alpha(1,i)=1.0d0/AR+(1.0d0-1.0d0/AR)*cos(2.0d0*pi*frac)**2.0d0
      end do
!
!     add a perturbation
!      alpha(1,10) = alpha(1,10) + 0.1
!
!
!    
!    
!    
!    
!    
!    
!    
!    
!    
!    
!    
!-----------------------------------------------------------------------
!    
      call main(nx,ny,na,nl,np,alpha,params)
!    
!    
!    
!    
      end program wrapper
!**********************************************************************!
