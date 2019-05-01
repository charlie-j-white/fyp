!**********************************************************************!
!                          Charlie Anderson                            !
!                  University of Bristol, March 2019                   !
!**********************************************************************!
!    
      subroutine adjoint(nx,ny,nl,flow,residual,np,params,dt,na,alpha,
     & Jcost,SENS,timeM)
!
      integer, external :: colour
!
      integer :: nx,ny,nl,np,na
      integer :: nh,nb,nt
      integer :: i,j,k,Rx,Ry,INFO,tmax
      integer :: cerr
!
      double precision :: Jcost,ans,ansT,ansN,ansI
      double precision :: dtS
      double precision :: con, conmin, conmax, contest
      double precision :: numzer, numnon
      double precision :: time0, timeDA, timeDB, timeDC, timeDD,
     &    timeLF, timeLS, timeS, timeM, time00, timeA
      double precision, dimension(1,np) :: params
!      double precision, dimension(0-nl:nx+nl,0-nl:ny+nl) :: 
!     &         meshX,meshY
!      double precision, dimension(1-nl:nx+nl,1-nl:ny+nl) :: 
!     &         u1,u2,u3,u5,
!     &         r1,r2,r3,r5
!
!
!
      double precision, dimension(1,na) :: 
     &         alpha,alpha_s
      double precision, dimension(1,4*(nx+2*nl)*(ny+2*nl)) :: 
     &         flow,flow_s,
     &         residual,
     &         dt
!
      double precision, dimension(4*nx*ny,na) :: dRda
      double precision, dimension(4*nx*ny,4*nx*ny) :: fluxjac
      integer, dimension(4*nx*ny,4*nx*ny) :: IPIV
      double precision, dimension(1,na) :: dJda,SENS
      double precision, dimension(1,4*nx*ny) :: adj, dJdw
!
      double precision, dimension(4*(nx+2*nl)*(ny+2*nl),na) :: 
     &             dRda_h
      double precision, dimension(4*(nx+2*nl)*(ny+2*nl),
     &                            4*(nx+2*nl)*(ny+2*nl)) :: 
     &             fluxjac_h
      double precision, dimension(1,4*(nx+2*nl)*(ny+2*nl)) :: 
     &             dJdw_h
!
!
!
!
!
      nb = nx*ny
      nh = 2*nl*(nx+ny)
      nt = 4*(nx+2*nl)*(ny+2*nl)
      print*, "      adjoint called from main . . ."
!
!
!
!
!
!
!-----------------------------------------------------------------------
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
!      call cost_is(nx,ny,nl,na,np,params,flow,alpha,Jcost)
!!
!      call aresid(nx,ny,nl,flow,residual,np,params,dt,na,alpha)
!!
!      call cost_is_d(nx, ny, nl, na, np, params, flow, flowd,
!     & alpha, alphad, jcost, jcostd)
!!
!      call aresid_d(nx, ny, nl, flow, flowd, residual, residuald, 
!     & np, params, dt, na, alpha, alphad)
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
!     get derivatives from the differentiated cost and resdifual
!     functions. should be using
!
!             dJda
!
!             dJdw_h
!
!     which are both multiple input single output functions, and
!
!             dRda_h
!
!             fluxjac_h
!
!     which are both multiple input multiple output functions.
!
!
!
      print*, "        calculate partial derivatives . . ."
!
!
!     ------- get dJda -------
!
!
!
      call CPU_TIME(time0)
      call CPU_TIME(time00)
      do i = 1,na
!
!
!     initialise all seeds as zero
      do j = 1,na
      alpha_s(1,j) = 0.0d0
      end do
      do j = 1,nt
      flow_s(1,j) = 0.d0
      end do
!
!
!     set relevant seed to one
      alpha_s(1,i) = 1.0d0
!
!
!     call cost function to relevant part of hold array
      call cost_is_d(nx, ny, nl, na, np, params, flow, flow_s,
     & alpha, alpha_s, jcost, dJda(1,i))
!
!
      end do
      call CPU_TIME(timeDA)
      timeDA = timeDA - time0
!
!
!
!
!
!
!
!     ------- get dJdw -------
!
!
!
      call CPU_TIME(time0)
      do i = 1,nt
!
!
!     initialise all seeds as zero
      do j = 1,na
      alpha_s(1,j) = 0.0d0
      end do
      do j = 1,nt
      flow_s(1,j) = 0.d0
      end do
!
!
!     set relevant seed to one
      flow_s(1,i) = 1.0d0
!
!
!     call cost function to relevant part of hold array
      call cost_is_d(nx, ny, nl, na, np, params, flow, flow_s,
     & alpha, alpha_s, jcost, dJdw_h(1,i))
!
!
!
!
      end do
      call CPU_TIME(timeDB)
      timeDB = timeDB - time0
!
!
!
!
!
!
!     ------- get dRda -------
!
!
!
      call CPU_TIME(time0)
      do i = 1,na
!
!
!     initialise all seeds as zero
      do j = 1,na
      alpha_s(1,j) = 0.0d0
      end do
      do j = 1,nt
      flow_s(1,j) = 0.d0
      end do
!
!
!     set relevant seed to one
      alpha_s(1,i) = 1.0d0
!
!
!     call cost function to relevant part of hold array
      call aresid_d(nx, ny, nl, flow, flow_s, residual, dRda_h(:,i), 
     & np, params, dt, na, alpha, alpha_s)
!
!
      end do
      call CPU_TIME(timeDC)
      timeDC = timeDC - time0
!
!
!
!
!
!
!     ------- get dRdw -------
!
!
!
      call CPU_TIME(time0)
      do i = 1,nt
!
!
!     initialise all seeds as zero
      do j = 1,na
      alpha_s(1,j) = 0.0d0
      end do
      do j = 1,nt
      flow_s(1,j) = 0.d0
      end do
!
!
!     set relevant seed to one
      flow_s(1,i) = 1.0d0
!
!
!     call cost function to relevant part of hold array
      call aresid_d(nx, ny, nl, flow, flow_s, residual, fluxjac_h(:,i), 
     & np, params, dt, na, alpha, alpha_s)
!
!
      end do
      call CPU_TIME(timeDD)
      timeDD = timeDD - time0
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
!     now delete all of the rows with halo or corner cells. note that nt
!     is defined as the storage vector including halo and corner cells,
!     while nw is just the body cells = 4*nx*ny
!
!               dJda     --->   dJda     *no action required*
!              ( 1,na)         ( 1,na)
!
!              dJdw_h    --->    adj
!              ( 1,nt)         ( 1,nw)
!
!              dRda_h    --->   dRda
!              (nt,na)         (nw,na)
!
!             fluxjac_h  --->  fluxjac
!              (nt,nt)         (nw,nw)
!
      nw = 4*nb
!
!
      Rx = 4*nh
      do i = 1,nw
      adj(1,i) = dJdw_h(1,i+Rx)
      end do
!
!
      Rx = 4*nh
      do i = 1,nw
      do j = 1,na
      dRda(i,j) = dRda_h(i+Rx,j)
      end do
      end do
!
!
      Rx = 4*nh
      Ry = 4*nh
      do j = 1,nw
      do i = 1,nw
      fluxjac(i,j) = fluxjac_h(i+Rx,j+Ry)
      end do
      end do
!
!
      print*, "        partial derivatives calculated."
!
!
!
!
!
!     also, make a function that counts the number of non-zero elements 
!     on each row or column of flux jacobian for debugging purposes
!
!
!
!
      if (0 == 1) then
!-------------------------------------------------
      print*, "      "
      print*, "      "
      print*, "flux jacobian visualisation"
      print*, "   index            non-zero                 non-zeroT",
     &          "                  NaN                      infty"
      do j=1,nw
!
      ans  = 0.0d0
      ansT = 0.0d0
      ansN = 0.0d0
      ansI = 0.0d0
      do i=1,nw
!        non-zero elements
        if( fluxjac(i,j) /= 0.0d0 ) then
        ans = ans + 1.0d0
        else
        continue
        end if
!
!        non-zero elements transpose
        if( fluxjac(j,i) /= 0.0d0 ) then
        ansT = ansT + 1.0d0
        else
        continue
        end if
!
!        search for NaN
        if( fluxjac(i,j) /= fluxjac(i,j) ) then
        ansN = ansN + 1.0d0
        else
        continue
        end if
!
!        search for infinities
        if( fluxjac(i,j) == fluxjac(i,j) - 1.0d0 ) then
        ansI = ansI + 1.0d0
        else
        continue
        end if
      end do
!
      print*, j, ans, ansT, ansN, ansI
      end do
!-------------------------------------------------
      end if
!
!
!
!
!
!
!
!
      conmin = 10000.0d0
      conmax = 0.0d0
      numzer = 0.0d0
      numnon = 0.0d0
!
      do i = 1,nw
      do j = 1,nw
!
!
      if (DABS(fluxjac(i,j)) .NE. 0.0d0) then
!
          contest = DABS(fluxjac(i,j))
!
          if (contest .GT. conmax) then
              conmax = contest
          end if
          if (contest .LT. conmin) then
              conmin = contest
          end if
!
          numnon = numnon + 1
!
      end if
!
!
      if (DABS(fluxjac(i,j)) .EQ. 0.0d0) then
          numzer = numzer + 1
      end if
!
!
      end do
      end do
!
!
      con = conmax/conmin
      print*
      print*, "conmax = ", conmax
      print*, "conmin = ", conmin
      print*, "   con = ", con
      print*
      print*, "numnon = ", numnon
      print*, "numzer = ", numzer
      print*, "numtot = ", numzer+numnon
      print*, "sparsity=", numzer/(numzer+numnon)
      print*
!
!      cerr = colour(nw,fluxjac)
!
!
!
!
!
!-----------------------------------------------------------------------
!
!     solve the adjoint equation. for simplicity LAPack will be used
!     but solution of the sparse flux jacobian can be done more
!     efficiently with different software.
!
!                   fluxjac * adjoint = dJdw
!
!     for greater compatibility with the solver software dJdw has
!     already been stored in the vector call adj.
!
!
!      print*, "   index         adjoint vector"
!
!
      do i = 1,nw
      dJdw(1,i) = adj(1,i)
      end do
!
!
!
!
      INFO = -5
      call CPU_TIME(time0)
      print*, "        begin flux Jacobian factorisation . . ."
      call DGETRF(nw,nw,fluxjac,nw,IPIV,INFO)
      print*, "        matrix solved with status", 
     &  INFO, "; success = '0';"
      call CPU_TIME(timeLF)
      timeLF = timeLF - time0
!
!
!
      call CPU_TIME(time0)
      print*, "        begin system solution . . ."
      call DGETRS('T',nw,1,fluxjac,nw,IPIV,adj,nw,INFO)
      print*, "        system solved with status", 
     &  INFO, "; success = '0';"
      call CPU_TIME(timeLS)
      timeLS = timeLS - time0
!
!
!      print*, "   index         adjoint vector"
!      do i = 1,nw
!      print*, i, adj(1,i)
!      end do
!
!
!
!
!
!
!
!!     define constants for iterative scheme
!      tmax = 500
!      dtS = 0.0001d0
!      Rx = 4*nh
!!
!!     start pseudo-timestepping
!      do j = 1,tmax
!!
!!     update scheme
!      do i = 1,nw
!!
!!        flux jacobian multiplication
!         ans = 0.0d0
!         do k = 1,nw
!            ans = ans + adj(1,k)*fluxjac(k,i)
!         end do
!!
!!        rest of calcs
!         adj(1,i) = adj(1,i) + dtS*dt(1,i+Rx)*(ans - dJdw(1,i))
!!
!      end do
!!
!      end do
!
!
!
!      stop
!
!
!
!
!
!
!
!
!
!     verify system solve
!      do i = 1,nw
!!
!      ans = 0.0d0
!      do j = 1,nw
!         ans = ans + adj(1,j)*fluxjac(i,j)
!      end do
!!
!      print*, i, dJdw(1,i), ans, dJdw(1,i)-ans
!!
!      end do
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
!     finally, form the total derivative by combining the remaining
!     partial derivatives.
!
!                 SENS = dJda + adjoint * dRda
!
!
!
      call CPU_TIME(time0)
      do i = 1,na
!
      ans = 0.0d0
      do j = 1,nw
      ans = ans + adj(1,j)*dRda(j,i)
      end do
!
!
      SENS(1,i) = dJda(1,i) - ans
!      print*, SENS(1,i)
      call CPU_TIME(timeS)
      call CPU_TIME(timeA)
      timeS = timeS - time0
      timeA = timeA - time00
!
!
      end do
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
      print*, "      adjoint completed in main."
!
!
      print*
      print*, "primal   =", timeM
      print*, "adjoint' =", timeA
      print*
      print*, "time DA  =", timeDA
      print*, "time DB  =", timeDB
      print*, "time DC  =", timeDC
      print*, "time DD  =", timeDD
      print*
      print*, "time LF  =", timeLF
      print*, "time LS  =", timeLS
      print*
      print*, "time S   =", timeS
!
!
      open(444, file='timeadjoint.dat')
!
!
      write(444,*)  na, ";",
     &              nx*ny, ";",
     &              nx, ";",
     &              ny, ";",
     &              timeM, ";",
     &              timeA, ";",
     &              timeDA, ";",
     &              timeDB, ";",
     &              timeDC, ";",
     &              timeDD, ";",
     &              timeLF, ";",
     &              timeLS, ";",
     &              timeS, ";"
!
!
!
      close(444)
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
!     create image of the flux jacobian
!-----------------------------------------------------------------------
!
!      cerr = colour(nw,fluxjac)
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
      end subroutine adjoint
