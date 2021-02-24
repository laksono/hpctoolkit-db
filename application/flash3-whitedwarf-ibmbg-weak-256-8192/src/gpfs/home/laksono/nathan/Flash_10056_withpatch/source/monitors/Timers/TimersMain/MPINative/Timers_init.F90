!!****if* source/monitors/Timers/TimersMain/MPINative/Timers_init
!!
!! NAME
!!  Timers_init
!!
!! SYNOPSIS
!!
!!  Timers_init(integer(IN) :: myPE,
!!              integer(IN) :: numProcs,
!!              real(OUT) :: initialWCTime)
!!
!! DESCRIPTION 
!!  
!!  Initialize the timer data structures.  This will
!!  essentially delete all information previously gathered by all timers
!!  and make it safe to start timers from scratch.  In the middle of 
!!  a run, for instance, this could be called once per timestep along with
!!  Timers_getSummary to get timer summary information for each timestep. 
!!  
!! ARGUMENTS
!!
!!  initialWCTime -- the initial wall clock time when this was called.
!!  numProcs -- the global number of processors.
!!  myPE -- local processor number.
!!  
!!
!!***

subroutine Timers_init(myPE,numProcs,initialWCTime)

  use Timers_data, ONLY: tmr_initDate, tmr_initTime, tmr_writeStatSummary, &
       tmr_eachProcWritesSummary, tmr_myPE, tmr_numProcs
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none
  integer, intent(in) :: myPE, numProcs
  real, intent(out) :: initialWCTime

  ! Everybody should know this
  tmr_myPE = myPE
  tmr_numProcs = numProcs

  call RuntimeParameters_get("writeStatSummary", tmr_writeStatSummary)
  call RuntimeParameters_get("eachProcWritesSummary", tmr_eachProcWritesSummary)

  call tmr_etime(initialWCTime)
  tmr_initTime = initialWCTime
  call current_date_time(tmr_initDate)

  call tmr_init()

end subroutine Timers_init
