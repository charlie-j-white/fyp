!**********************************************************************!
!                          Charlie Anderson                            !
!                  University of Bristol, March 2019                   !
!**********************************************************************!
!    
      subroutine adjoint(nx,ny,nl,flow,residual,np,params,dt,na,alpha,
     & Jcost)
!
!
      integer :: nx,ny,nl,np,na
      integer :: nh,nb,nt
      integer :: i,j,Rx,Ry,INFO
!
      double precision :: Jcost,ans,ansT,ansN,ansI
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
      double precision, dimension(1,4*nx*ny) :: adj
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
      nt = 4*nb + 4*nh + 16*nl*nl
      print*
      print*, "      adjoint called from main . . ."
      print*
      print*, "  4*nh = ", 4*nh
      print*, "  4*nb = ", 4*nb
      print*, "16*nl2 = ", 16*nl*nl
      print*, "    nt = ", nt
      print*
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
!
!
!     ------- get dJda -------
!
!
!
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
!      call cost_is_d(nx, ny, nl, na, np, params, flow, flow_s,
!     & alpha, alpha_s, jcost, dJda(1,i))
!
!
      end do
!
!
!
      print*, "dsfsdafdsafdsfadsafa"
!
!
      stop
!
!
!     ------- get dJdw -------
!
!
!
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
!      print*, dJdw_h(1,i)
!
!
      end do
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
!
!
!
!
!
!     also, make a function that counts the number of non-zero elements 
!     on each row or column of flux jacobian for debugging purposes
!
!
      if (1 == 0) then
!-------------------------------------------------
      print*, "      "
      print*, "      "
      print*, "flux jacobian visualisation"
      print*, "   index            non-zero                 non-zeroT",
     &          "                  NaN                      infty"
      do j=1,nt
!
      ans  = 0.0d0
      ansT = 0.0d0
      ansN = 0.0d0
      ansI = 0.0d0
      do i=1,nt
!        non-zero elements
        if( fluxjac_h(i,j) /= 0.0d0 ) then
        ans = ans + 1.0d0
        else
        continue
        end if
!
!        non-zero elements transpose
        if( fluxjac_h(j,i) /= 0.0d0 ) then
        ansT = ansT + 1.0d0
        else
        continue
        end if
!
!        search for NaN
        if( fluxjac_h(i,j) /= fluxjac_h(i,j) ) then
        ansN = ansN + 1.0d0
        else
        continue
        end if
!
!        search for infinities
        if( fluxjac_h(i,j) == fluxjac_h(i,j) - 1.0d0 ) then
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
!
!
!
      INFO = -5
      print*, "      "
      print*, "      begin flux Jacobian factorisation . . ."
      print*,"      matrix solved with status", INFO, "; success = '0';"
      call DGETRF(nw,nw,fluxjac,nw,IPIV,INFO)
      print*, "      begin system solution . . ."
      call DGETRS('T',nw,1,fluxjac,nw,IPIV,adj,nw,INFO)
      print*,"      system solved with status", INFO, "; success = '0';"
      print*, "      "
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
!-----------------------------------------------------------------------
!
!     finally, form the total derivative by combining the remaining
!     partial derivatives.
!
!                 SENS = dJda + adjoint * dRda
!
!
!
!
!
      print*
!
!
      do i = 1,na
!
      ans = 0.0d0
      do j = 1,nw
      ans = ans + adj(1,j)*dRda(j,i)
      end do
!
!
      SENS(1,i) = dJda(1,i) + ans
      print*, SENS(1,i)
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
      print*, "      adjoint completed in main . . ."
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