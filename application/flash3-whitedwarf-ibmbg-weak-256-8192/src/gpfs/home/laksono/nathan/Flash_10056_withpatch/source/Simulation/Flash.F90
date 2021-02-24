!!****f* source/Simulation/Flash
!!
!! NAME
!!
!!  Flash
!!
!!
!! SYNOPSIS
!!
!!  N/A
!!
!!
!! DESCRIPTION
!!
!!  The source file Flash.F90 in the Simulation unit contains the Fortran
!!  PROGRAM. As such it can be considered the top-level "driver" of an application.
!!  By default it is set up to drive the simulation of a time-dependent
!!  problem by calling:
!!  - Driver_initFlash  for initializations,
!!  - Driver_evolveFlash  for managing the computation, and
!!  - Driver_finalizeFlash  for cleaning up.
!!
!! SEE ALSO
!!
!!  Driver_initFlash
!!  Driver_evolveFlash
!!  Driver_finalizeFlash
!!
!!***

program Flash

  implicit none


  call Driver_initFlash()

  call Driver_evolveFlash( )

  call Driver_finalizeFlash ( )
  

end program Flash
