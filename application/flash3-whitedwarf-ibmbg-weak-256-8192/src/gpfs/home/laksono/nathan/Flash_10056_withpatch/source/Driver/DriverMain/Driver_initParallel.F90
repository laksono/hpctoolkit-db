!!****if* source/Driver/DriverMain/Driver_initParallel
!!
!! NAME
!!
!!  Driver_initParallel
!!
!! SYNOPSIS
!!
!!  Driver_initParallel(integer(OUT) :: myPE,
!!                      integer(OUT) :: numProcs)
!!
!! DESCRIPTION
!!
!!  Initialize the parallel message-passing interface,
!!  the number of processors in a run and each processing
!!  element
!!
!!
!!  ARGUMENTS
!!    myPE : current processor
!!    numProcs : number of processors
!!  
!!
!!
!!***

subroutine Driver_initParallel (myPE, numProcs)

  implicit none             

  include "Flash_mpi.h"

  integer, intent(out) :: myPE, numProcs

  integer   error

  call MPI_Init (error)

  call MPI_Comm_Rank (MPI_COMM_WORLD, myPE, error)
  call MPI_Comm_Size (MPI_COMM_WORLD, numProcs, error)

  return
end subroutine Driver_initParallel
