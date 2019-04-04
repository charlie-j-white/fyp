!**********************************************************************!
!                          Charlie Anderson                            !
!                  University of Bristol, March 2019                   !
!**********************************************************************!
!    
      program wrapper
!    
      integer :: na = 5
      integer :: np = 13
      double precision, dimension(1,5) :: alpha, alpha_fd,
     &   dJda_fd, dJda_ad
      double precision, dimension(1,13) :: params
!    
      integer :: nx, ny, nl, runtype, i, j
      double precision :: AR, pi, frac, h, Jcost, Jcost0
!    
      print*
      print*, "START PROGRAM"
      print*
      print*, '                 Discrete Adjoint Solver'
      print*
      print*
      print*
!      
!    
!    
!    
!     define what type of run the code will have 
!-----------------------------------------------------------------------
!
      runtype = 0
!    
!             = 0: Normal operation    INPUT: alpha;
!                                     OUTPUT: cost; solution &
!                                             convergence files;
!    
!             = 1: Finite difference   INPUT: alpha; 
!                                     OUTPUT: cost; sensitivity;
!      
!             = 2: Adjoint solver      INPUT: alpha;
!                                     OUTPUT: cost; sensitivity;
!    
!             < 0: No operation - displays help message 
!    
!    
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
      nx = 10
      ny = 10
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
      params(1,5) = DBLE(20000)
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
      params(1,9) = 0.71d0
      params(1,10) = 5.0d0
      params(1,11) = 1.0d0/25.0d0
      params(1,12) = 0.2d0
!    
!    
!     residual_tol
!    
      params(1,13) = 0.00000001d0
!      params(1,13) = 0.000001d0
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
      if (runtype .EQ. 0) then
      call main(nx,ny,na,nl,np,alpha,params,runtype,Jcost)
      end if
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
!    
!     finite difference section
!-----------------------------------------------------------------------
!
      if (runtype .EQ. 1) then
!    
!   
!   
!
!     set up finite difference calculations 
      print*, "Start Finite difference calcs . . ."
      print*
!   
!   
!     step size of 0.001 seems to work best   
      h = 0.001d0
      print*, "      h = ", h
!   
!    
!     calculate initial cost function for forward difference
      call main(nx,ny,na,nl,np,alpha,params,runtype,Jcost0)
!    
!
!
!
!     begin looping over design variables   
      do i = 1,na
      print*, "      difference",i,"of",na
!
!
!     re-initialise variable vector
      do j = 1,na
      alpha_fd(1,j) = alpha(1,j)
      end do
!    
!    
!     add perturbation
      alpha_fd(1,i) = alpha(1,i) + h
!     
!    
!     call function with relevant perturbation    
      call main(nx,ny,na,nl,np,alpha_fd,params,runtype,Jcost)
!    
!    
!     calculate actual finite difference for this design variable 
      dJda_fd(1,i) = (Jcost - Jcost0)/h
!    
!    
      end do
      end if
!    
!    
!    
!    
!    
!    
!    
!-----------------------------------------------------------------------
!    
!    
!    
!    
!    
!    
!    
!     display results    
!-----------------------------------------------------------------------
!
      if (runtype .NE. 0) then    
!
      print*    
      print*    
      print*, "Results of sensitivity calcs . . ."
      print*    
      print*, "                 FD                      Adjoint" 
!    
!    
!    
!
      do i = 1,na
!
      print*, "      ", dJda_fd(1,i)
!
      end do    
!    
!    
!    
!    
!    
!    
!    
      end if
!    
!-----------------------------------------------------------------------
!
      print*
      print*
      print*, "END PROGRAM"    
      print*
!    
      end program wrapper
!**********************************************************************!
