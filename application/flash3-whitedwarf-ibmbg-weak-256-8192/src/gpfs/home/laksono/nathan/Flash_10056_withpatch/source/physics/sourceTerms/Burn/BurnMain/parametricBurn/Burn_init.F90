!!****if* source/physics/sourceTerms/Burn/BurnMain/parametricBurn/Burn_init
!!
!! NAME
!!  
!!  Burn_init
!!
!! SYNOPSIS
!! 
!!  call Burn_init(integer(IN) :: myPE)
!!  
!! DESCRIPTION
!!
!!  Initial data and runtime parameteres used by the
!!  multi-stage parameteric burning routines.
!!
!! ARGUMENTS
!!
!!  myPE -- current local processor
!!
!! PARAMETERS
!!
!!  useBurn -- Boolean, True.  Turns on burning module
!!  useShockBurn -- Boolean, FALSE.  Controls whether burning is allowed inside
!!                a regime experiencing shocks
!!
!!
!!
!!***
subroutine Burn_init(myPE)
  
  use Burn_data, ONLY                   : bn_useBurn, bn_myPE, bn_useShockBurn, &
                                          bn_smallx,bn_enucDtFactor
  use Driver_interface, ONLY            : Driver_abortFlash
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  
  implicit none
#include "constants.h"
#include "Flash.h"
  integer, intent(IN) :: myPE

  ! This strange initialization is to keep interface 
  ! consistent with other sourceTerms
  bn_myPE = myPE
  
  call RuntimeParameters_get("useBurn", bn_useBurn )
  if (.not. bn_useBurn) then
     write(6,*)' WARNING:  You have included the Burn unit but have set '
     write(6,*)'           the runtime parameter useBurn to FALSE'
     write(6,*)'           No burning will occur but Burn_init will continue.'
  end if
  call RuntimeParameters_get('useShockBurn', bn_useShockBurn)
  call RuntimeParameters_get('smallx', bn_smallx)
  call RuntimeParameters_get('enucDtFactor', bn_enucDtFactor)
  
  !!  Now initialize the network and nse things
  call bn_nseInit(myPE)
  call bn_paraInit(myPE)
  
end subroutine Burn_init
