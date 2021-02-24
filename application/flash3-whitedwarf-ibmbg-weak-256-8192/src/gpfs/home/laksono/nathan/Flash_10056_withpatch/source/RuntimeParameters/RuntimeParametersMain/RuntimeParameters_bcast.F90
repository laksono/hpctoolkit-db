!!****if* source/RuntimeParameters/RuntimeParametersMain/RuntimeParameters_bcast
!!
!! NAME
!!  RuntimeParameters_bcast
!!
!! SYNOPSIS
!!
!!  RuntimeParameters_bcast(integer(in) :: myPE)
!!
!! DESCRIPTION
!!
!! Broadcasts parameters from the current processor to the other processors.
!! Only the master processor MASTER_PE reads in runtime parameters.
!! 
!! ARGUMENTS
!!
!!        
!! myPE:      current processor number
!!
!!
!!
!!***

subroutine RuntimeParameters_bcast(myPE)

  use RuntimeParameters_data, ONLY : parameter

  implicit none
  
  integer, intent(in)                      :: myPE

  call nameValueLL_bcast(parameter, myPE)

  return

end subroutine RuntimeParameters_bcast


  

