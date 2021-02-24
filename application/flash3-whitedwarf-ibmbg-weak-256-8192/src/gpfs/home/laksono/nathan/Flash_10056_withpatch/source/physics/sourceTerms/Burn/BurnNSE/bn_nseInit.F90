!!****if* source/physics/sourceTerms/Burn/BurnNSE/bn_nseInit
!!
!! NAME
!!
!!  bn_nseInit
!!
!! SYNOPSIS
!!
!!  bn_nseInit(  integer, intent(in)     :: myPE)
!!
!! DESCRIPTION
!!
!!  Initialize NSE Data.
!!  C.A.Meakin, 26 Feb. 2007 (Chicago)
!!
!! ARGUMENTS
!!
!!   myPE:  current processor
!!
!!***


subroutine bn_nseInit(myPE)

  use bn_nseInterface, ONLY : bn_nseInitTables
  use Driver_interface, ONLY : Driver_abortFlash
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, ONLY : Logfile_stampMessage
  use Multispecies_interface, ONLY : Multispecies_getProperty
  use PhysicalConstants_interface, ONLY:  PhysicalConstants_get

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent(IN)     :: myPE
  character (len=50),save :: denstablename,prestablename


  !READ IN NSE TABLE DATA
  call RuntimeParameters_get( "pbPresTableName", prestablename)
  call RuntimeParameters_get( "pbDensTableName", denstablename)

  call bn_nseInitTables(myPE,prestablename,denstablename)



end subroutine bn_nseInit

