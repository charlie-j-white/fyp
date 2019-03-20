!***********************************************************************
!  CHRISTIAN ALLEN. UNIVERSITY OF BRISTOL.
!  LAST UPDATE: NOVEMBER 2014.
!***********************************************************************
!
       subroutine readconfig(
     &      pmachin,arearat,prat,CFLscale,pk,nstages,maxit,ni,restol)
!
      IMPLICIT NONE
!
      LOGICAL,DIMENSION(9):: READIN
      CHARACTER*80:: linestring,linestring1
      CHARACTER*50,DIMENSION(9):: label
!
      DOUBLE PRECISION:: restol,pmachin,CFLscale,arearat,prat,pk
!
      INTEGER:: maxit,ni,nstages,i,j,n
!
        do n=1,9
          READIN(n)=.FALSE.
        enddo
!
       label(1)='FLOW.MACH'
       label(2)='FLOW.AREARATIO'
       label(3)='FLOW.PRESSURERATIO'
       label(4)='FLOW.NPOINTS'
       label(5)='SOLVER.CFL_SCALE'
       label(6)='SOLVER.RESID_TOLERANCE'
       label(7)='SOLVER.MAX_CYCLES'
       label(8)='SOLVER.STAGES'
       label(9)='SOLVER.KAPPA'
!
        open(10,file='solver.conf')
!
88      linestring(:)=' '
        linestring1(:)=' '
        read(10,101,ERR=2000,END=2000) linestring
        linestring1=linestring
        call TOUC(linestring) 
        do i=1,79
          if(linestring(i:i).eq.'#') then
            do j=i,80
              linestring(j:j)=' '
            enddo
          endif
        enddo
!
        do i=1,50 
!
          if(linestring(i:i+8).eq.'FLOW.MACH') then 
            do j=i,i+8
              linestring(j:j)=' '
            enddo
            do j=1,65
              if(linestring(j:j).eq.'=') then 
                linestring(j:j)=' '
              endif
            enddo
            read(linestring,*) pmachin 
            READIN(1)=.TRUE.
!
          elseif(linestring(i:i+13).eq.'FLOW.AREARATIO') then 
            do j=i,i+13
              linestring(j:j)=' '
            enddo
            do j=1,65
              if(linestring(j:j).eq.'=') then 
                linestring(j:j)=' '
              endif
            enddo
            read(linestring,*) arearat
            READIN(2)=.TRUE.
!
          elseif(linestring(i:i+17).eq.'FLOW.PRESSURERATIO') then 
            do j=i,i+17
              linestring(j:j)=' '
            enddo
            do j=1,65
              if(linestring(j:j).eq.'=') then 
                linestring(j:j)=' '
              endif
            enddo
            read(linestring,*) prat 
            READIN(3)=.TRUE.
!
          elseif(linestring(i:i+11).eq.'FLOW.NPOINTS') then 
            do j=i,i+11
              linestring(j:j)=' '
            enddo
            do j=1,65
              if(linestring(j:j).eq.'=') then 
                linestring(j:j)=' '
              endif
            enddo
            read(linestring,*) ni 
            READIN(4)=.TRUE.
!
          elseif(linestring(i:i+15).eq.'SOLVER.CFL_SCALE') then 
            do j=i,i+15
              linestring(j:j)=' '
            enddo
            do j=1,65
              if(linestring(j:j).eq.'=') then 
                linestring(j:j)=' '
              endif
            enddo
            read(linestring,*) cflscale
            READIN(5)=.TRUE.
!
          elseif(linestring(i:i+21).eq.'SOLVER.RESID_TOLERANCE') then 
            do j=i,i+21
              linestring(j:j)=' '
            enddo
            do j=1,65
              if(linestring(j:j).eq.'=') then 
                linestring(j:j)=' '
              endif
            enddo
            read(linestring,*) restol 
            READIN(6)=.TRUE.
!
          elseif(linestring(i:i+16).eq.'SOLVER.MAX_CYCLES') then 
            do j=i,i+16 
              linestring(j:j)=' '
            enddo
            do j=1,65
              if(linestring(j:j).eq.'=') then 
                linestring(j:j)=' '
              endif
            enddo
            read(linestring,*) maxit 
            READIN(7)=.TRUE.
!
          elseif(linestring(i:i+12).eq.'SOLVER.STAGES') then 
            do j=i,i+12
              linestring(j:j)=' '
            enddo
            do j=1,65
              if(linestring(j:j).eq.'=') then 
                linestring(j:j)=' '
              endif
            enddo
            read(linestring,*) nstages 
            READIN(8)=.TRUE.
!
          elseif(linestring(i:i+11).eq.'SOLVER.KAPPA') then 
            do j=i,i+11
              linestring(j:j)=' '
            enddo
            do j=1,65
              if(linestring(j:j).eq.'=') then 
                linestring(j:j)=' '
              endif
            enddo
            read(linestring,*) pk 
            READIN(9)=.TRUE.
!
        ENDIF
!
       ENDDO 
!
       go to 88 
!
2000    continue
!
        close(10)
!
!
        do n=1,9
          if(READIN(n).eqv..FALSE.) THEN
            print*,'ERROR: VARIABLE:' 
            print*,label(n)
            print*,'NOT READ IN CORRECTLY'
            stop
          endif
        enddo 
!
101     format((a))
!
        return
        end subroutine readconfig 
!
!**********************************************************************
!***********************************************************************
!  CHRISTIAN ALLEN. UNIVERSITY OF BRISTOL.
!  LAST UPDATE: NOVEMBER 2014.
!***********************************************************************
!
      subroutine touc(AIN)
      IMPLICIT NONE
      CHARACTER*80:: ain
      INTEGER:: l,lenain,ascbot,asctop,asclow,iain
!
      DATA ascbot/97/,asctop/122/,asclow/32/
!
      lenain = len(ain)
      do l=1,lenain
         iain = ichar(ain(l:l))
         if (iain.ge.ascbot.and.iain.le.asctop) then
            ain(l:l) = char(iain-asclow)
         endif
      enddo
      return
      end subroutine touc 
!
!**********************************************************************

