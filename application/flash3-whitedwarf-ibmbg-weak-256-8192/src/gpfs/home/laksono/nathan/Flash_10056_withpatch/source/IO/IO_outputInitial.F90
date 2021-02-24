!!****f* source/IO/IO_outputInitial
!!
!! NAME
!!
!!  IO_outputInitial
!!
!!
!! SYNOPSIS
!!
!!  IO_outputInitial(integer(in) :: myPE, 
!!                   integer(in) :: numProcs,
!!                   integer(in) :: nbegin,
!!                   real(in) :: initialSimTime  
!!                  
!!
!!
!! DESCRIPTION
!!
!!  This routine is called before the main timestep loop.  It outputs the 
!!  initial data to a checkpoint file and plotfile, and particle plotfiles
!!
!!  If particles are not included a stub (empty) routine will be called.
!!
!!
!! ARGUMENTS
!!
!!  myPE - current processor
!!  numProcs - number of processors running the simulation
!!  nbegin - initial step of simulation
!!  initialSimTime - initial simulation time
!!
!!
!!***

subroutine IO_outputInitial(myPE, numProcs, nbegin, initialSimTime)

  use IO_interface, ONLY : IO_writeIntegralQuantities

  implicit none
  integer, intent(in) :: myPE, numProcs, nbegin
  real, intent(in)    :: initialSimTime

  call IO_writeIntegralQuantities(myPE, 1, initialSimTime)

end subroutine IO_outputInitial
