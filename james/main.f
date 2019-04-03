!**********************************************************************!
!                          Charlie Anderson                            !
!                  University of Bristol, March 2019                   !
!**********************************************************************!
!                                                          
      subroutine main(nx,ny,na,nl,np,alpha,params,runtype,Jcost)
!
!
      integer :: nx,ny,na,nl,np
      integer :: i,j,tmax,runtype
!
      double precision :: rmsmax,rms,Jcost
!
      double precision, dimension(1,na) :: alpha
      double precision, dimension(1,np) :: params
      double precision, dimension(1,4*(nx+2*nl)*(ny+2*nl)) ::
     & flow,flow0,residual,residual0
      double precision, dimension(1-nl:nx+nl,1-nl:ny+nl) ::
     & u1,u2,u3,u5,dt,
     & r1,r2,r3,r5
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
!     And so begins the main function:
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
!     start meshing and initialisation processes
!-----------------------------------------------------------------------
!
!
!     open file for convergence information
      if (runtype .EQ. 0) then
      open(456, file='convergence.plt')
      write(456,*) 'VARIABLES = "it" "log(rms)"'
457   format(2f15.8)
      end if
!
!
!     initialise solution and residual vectors
      do j = 1,4*(nx+2*nl)*(ny+2*nl)
      flow(1,j) = 0.0d0
      residual(1,j) = 0.0d0
      end do
!
!
      call meshing(nx,ny,na,nl,alpha,meshX,meshY,np,params)
      call initialise(nx,ny,nl,flow,np,params,meshX,meshY)
!
!
!
!
!      do j = 0-nl,ny+nl
!      print*, meshX(:,j)
!      end do
!      print*
!      do j = 0-nl,ny+nl
!      print*, meshY(:,j)
!      end do
!      print*
!
!
      print*
!
!
      call meshing(nx,ny,na,nl,alpha,meshX,meshY,np,params)
!
!      do j = 0-nl,ny+nl
!      print*, meshX(:,j)
!      end do
!      print*
!      do j = 0-nl,ny+nl
!      print*, meshY(:,j)
!      end do
!      print*
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
!
!
!
!     start timestepping 
!-----------------------------------------------------------------------
!
!
!
      rmsmax = -10.0d0
      tmax = INT(params(1,5))
      do i = 1,tmax
!
!
      call timestep(nx,ny,nl,flow,np,params,meshX,meshY,dt)
!
      do j = 1,4*(nx+2*nl)*(ny+2*nl)
      flow0(1,j) = flow(1,j)
      end do
!
      call resid(nx,ny,nl,flow,residual,np,params,meshX,meshY,dt)
!
      call update(nx,ny,nl,flow,residual,np,params,dt,flow0)
!
      call itinfo(nx,ny,nl,flow,flow0,dt,i,rmsmax,rms,runtype)
!
!
!     write to file if a runtype of 0
      if (runtype .EQ. 0) then
      write(456,457) DBLE(i), DLOG10(rms/rmsmax)
      end if
!
!
!     escape loop if RMS error is low enough
      if (rms/rmsmax .LE. params(1,13)) then
        goto 99
      end if
!
!
      end do
99    continue
!
!
!
!
!
!
!
!     end of main calcs, post-processing
!-----------------------------------------------------------------------
!
!     post-process the results if runtype is 0
      if (runtype .EQ. 0) then
      call postprocess(nx,ny,nl,flow,u1,u2,u3,u5,meshX,meshY,rmsmax,
     & flow0)
      end if
!
!     close convergence file if a runtype of 0
      if (runtype .EQ. 0) then
      close(456)
      end if
!
!
!
!
!
!
!
!     call cost function from here
!-----------------------------------------------------------------------
!
!     bear in mind relationship between flow variables w and design
!     variables alpha for the AD. define:
!
!            inheritance: w inherits effects from alpha
!              w == w(a)
!
!              isolation: w is not affected by alpha
!              w /= w(a)
!
!
!
!     1) Create functions that will give the correct partial derivatives
!
!     resid.f          aresid.f          aMresid.f          aMresid_d.f
!     cost.f           acost.f           aMcost.f           aMcost_d.f
!           ----------->      ----------->       ----------->
!          isolate w and       | merge.sh          Tapenade
!            a by hand         |
!                              `--> edit 
!                             subroutine names
!
!
!     2) Take the partial derivatives from the functions in hold arrays
!        and delete rows with boundary conditions.
!
!
!     3) Solve adjoint equation.
!
!
!     4) Form total derivative.
!
!
!
!
!
!
!
!      call split_fwd(nx,ny,nl,residual,r1,r2,r3,r5)
!      call debug_cell(nx,ny,nl,r1)
!      call debug_cell(nx,ny,nl,r2)
!      call debug_cell(nx,ny,nl,r3)
!      call debug_cell(nx,ny,nl,r5)
!
!
!
!
!
!
!
      Jcost = 0.0d0;
!
      call cost_is(nx,ny,nl,na,np,params,flow,alpha,Jcost)
      call aresid(nx,ny,nl,flow0,residual,np,params,dt,na,alpha)
!
!
!      call split_fwd(nx,ny,nl,residual,r1,r2,r3,r5)
!      call debug_cell(nx,ny,nl,r1)
!      call debug_cell(nx,ny,nl,r2)
!      call debug_cell(nx,ny,nl,r3)
!      call debug_cell(nx,ny,nl,r5)
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
!      call adjoint function from here
!-----------------------------------------------------------------------
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
!      de-bugging functions, if needed 
!-----------------------------------------------------------------------
!
!      call split_fwd(nx,ny,nl,flow,u1,u2,u3,u5)
!      call debug_cell(nx,ny,nl,u1)
!      call debug_cell(nx,ny,nl,u2)
!      call debug_cell(nx,ny,nl,u3)
!      call debug_cell(nx,ny,nl,u5)
!!
!      call split_fwd(nx,ny,nl,residual,r1,r2,r3,r5)
!      call debug_cell(nx,ny,nl,r1)
!      call debug_cell(nx,ny,nl,r2)
!      call debug_cell(nx,ny,nl,r3)
!      call debug_cell(nx,ny,nl,r5)
!
!      call split_fwd(nx,ny,nl,flow,u1,u2,u3,u5)
!      call split_rev(nx,ny,nl,flow,u1,u2,u3,u5)
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
