!     this program was adapted with help from 'pmy' on stack overflow
!     URL: https://stackoverflow.com/questions/22028571/track-memory-usage-in-fortran-90
!
!
!
!
      program main 
!
      implicit none
!
      integer :: valueRSS
      character(len=200):: filename=' '
      character(len=80) :: line
      character(len=8)  :: pid_char=' '
      integer :: pid
      logical :: ifxst
!
!
      integer :: i, j, tmax, sle
      double precision :: sLen
!
!
      valueRSS=-1    ! return negative number if not found
!
      pid = 9999999
      tmax = 100
      sle = 1
!
!
!
!
!
!     remove traces of old programs
      call EXECUTE_COMMAND_LINE('rm PROCESS')
!
!
!
!     stuff to add to process being studied
!______________________________________
!
!     get process ID
      pid=getpid()
!
!     write to output file
      open(222, file='PROCESS')
      write(222,*) pid
      close(222)
!
!
!______________________________________
!
!
!
      print*, pid
!
!
!
!     do some long process until cancelled
      do i = 1,tmax
!
      call SLEEP(3)
!
      end do
!
!
!
      end program main 
