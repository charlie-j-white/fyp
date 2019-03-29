!**********************************************************************!
!                          Charlie Anderson                            !
!                  University of Bristol, March 2019                   !
!**********************************************************************!
!                                                          
      subroutine main(nx,ny,na,nl,np,alpha,params)
!
!
      integer :: nx,ny,na,nl,np
      integer :: i,j,tmax
!
      double precision :: dt,rmsmax
!
      double precision, dimension(1,na) :: alpha
      double precision, dimension(1,np) :: params
      double precision, dimension(1,4*(nx+2*nl)*(ny+2*nl)) ::
     & flow,flow0,residual
      double precision, dimension(1-nl:nx+nl,1-nl:ny+nl) ::
     & u1,u2,u3,u5,r1,r2,r3,r5
      double precision, dimension(0-nl:nx+nl,0-nl:ny+nl) ::
     & meshX,meshY
!
!
!
!
!     A cell is defined as the space between mesh nodes. The labelling 
!     for a general cell is as follows:
!
!
!                     mesh(i-1,j)             mesh(i,j)
!                          *-------------------*
!                          |                   |
!                          |                   |
!                          |                   |
!                          |     cell(i,j)     |
!                          |                   |
!                          |                   |
!                          |                   |
!                          *-------------------*
!                     mesh(i-1,j-1)           mesh(i,j-1)
!                   
!
!     The number of halo cells in each direction (i,j) is defined by nl.
!     This is shown for on layer of halo cells:
!
!
!     Ny+1  *------*------*------*------*------*------*------*------*
!           |                                                       |
!       Ny  *      *------*------*------*------*------*------*      *
!           |      | 1,Ny                              Nx,Ny |      |
!     Ny-1  *      *      *      *      *      *      *      *      *
!           |      |                                         |      |
!           *      *      *      *      *      *      *      *      *
!           |      |                                         |      |
!           *      *      *      *      *      *      *      *      *
!           |      |                                         |      |
!        1  *      *      *      *      *      *      *      *      *
!           |      | 1,1                                Nx,1 |      |
!        0  *      *------*------*------*------*------*------*      *
!           | 0,0                                                   |
!       -1  *------*------*------*------*------*------*------*------*
!          -1      0      1                         Nx-1    Nx    Nx+1        
!
!
!     And so begins the main function:
!
!      call split_fwd(nx,ny,nl,flow,u1,u2,u3,u5)
!      call split_rev(nx,ny,nl,flow,u1,u2,u3,u5)
!
      print*, "START PROGRAM main(...)"
      print*, "  "
!
!
!
!
!     start meshing and initialisation processes
!-----------------------------------------------------------------------
!
      call meshing(nx,ny,na,nl,alpha,meshX,meshY)
!
      call initialise(nx,ny,nl,flow,np,params,meshX,meshY)
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
!     start timestepping section 
!-----------------------------------------------------------------------
!
      rmsmax = -10.0d0
      tmax = 1
      do i = 1,tmax
!
      call timestep(nx,ny,nl,flow,flow0,dt,i,rmsmax)
!
      do j = 1,4*(nx+2*nl)*(ny+2*nl)
      flow0(1,j) = flow(1,j)
      end do
!
      call resid(nx,ny,nl,flow,residual,np,params,meshX,meshY)
!
      call update(nx,ny,nl,flow,residual,np,params,meshX,meshY,dt)
!
      call itinfo(nx,ny,nl,flow,flow0,dt,i,rmsmax)
!
      end do
!
!
!
!
!
!     end of main calcs, debugging and post-processing
!-----------------------------------------------------------------------
!
      call split_fwd(nx,ny,nl,flow,u1,u2,u3,u5)
!
!      call debug_cell(nx,ny,nl,u1)
!      call debug_cell(nx,ny,nl,u2)
!      call debug_cell(nx,ny,nl,u3)
!      call debug_cell(nx,ny,nl,u5)
!
      call postprocess(nx,ny,nl,flow,u1,u2,u3,u5,meshX,meshY)
!
!
!
!
!
!
!
!
!
      end subroutine main
!-----------------------------------------------------------------------
