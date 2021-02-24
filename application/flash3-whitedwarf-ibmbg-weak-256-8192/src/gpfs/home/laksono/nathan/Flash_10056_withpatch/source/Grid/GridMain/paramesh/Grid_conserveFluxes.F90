!!****if* source/Grid/GridMain/paramesh/Grid_conserveFluxes
!!
!! NAME
!!  Grid_conserveFluxes
!!
!! SYNOPSIS
!!
!!  Grid_conserveFluxes(integer(IN) :: myPE, 
!!                      integer(IN) :: axis,
!!                      integer(IN) :: level)
!!  
!! DESCRIPTION 
!!  
!!  Flux conservation is necessary when 2 blocks of differing
!!  levels (meaning having different grid spacings) border 
!!  one another. 
!!  
!!  This routine can perform flux conservation on the finest
!!  blocks, the most typical usage for the Paramesh Grid or on
!!  blocks of a certain level.
!!  
!!  The routine overwrites the flux arrays maintained by the Grid
!!  
!! ARGUMENTS 
!!
!!  myPE - current processor number
!!
!!  axis - conserve fluxes in just one direction if 
!!         IAXIS, JAXIS, KAXIS, or in all directions if ALLDIR.
!!         These constants are defined in constants.h.
!!
!!  level - refinement level. Ignored.
!!
!!***
subroutine Grid_conserveFluxes(MyPE, axis, level)
  use paramesh_interfaces, ONLY : amr_flux_conserve
#include "Flash.h"
#include "constants.h"
#ifndef FLASH_GRID_PARAMESH2
  use physicaldata, ONLY: no_permanent_guardcells
#endif
implicit none
  integer, intent(in) :: MyPE, axis, level
  integer :: gridDataStruct

  call amr_flux_conserve(MyPE, 0, axis)

#ifndef FLASH_GRID_PARAMESH2
  gridDataStruct = CENTER
#if NFACE_VARS > 0
  gridDataStruct = CENTER_FACES
#endif
  if (no_permanent_guardcells) then
     call gr_commSetup(gridDataStruct)
  end if
#endif
end subroutine Grid_conserveFluxes
