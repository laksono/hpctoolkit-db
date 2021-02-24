!!****if* source/Grid/GridMain/paramesh/paramesh4/Grid_getFluxData
!!
!! NAME
!!  Grid_getFluxData
!!
!! SYNOPSIS
!!
!!
!!  Grid_getFluxData(integer(IN) :: blockID,
!!                   integer(IN) :: axis,
!!                   real(INOUT) :: fluxes(NFLUXES,dataSize(1),dataSize(2),dataSize(3)),
!!                   integer(IN) :: dataSize(3),
!!          OPTIONAL,integer(IN) :: pressureSlot,
!!          OPTIONAL,real(IN)    :: areaLeft(:,:,:))
!!
!! DESCRIPTION 
!!
!!  Get the fluxes in a direction specified by axis for boundary cells
!!  for block blockID. This routine needs to be used when using adaptive mesh
!!  since fluxes calculated by the two blocks that are at a fine/coarse
!!  boundary have different accuracy.
!!  
!!  This should be called after Grid_conserveFluxes, which gets the 
!!  fluxes updated for cells at fine/coarse boundaries to make them consistent.
!!
!! ARGUMENTS
!!
!!  blockID : The local blockid
!!
!!
!!  axis : integer value specifying on which cell faces to get fluxes. 
!!         The options are IAXIS, JAXIS, or KAXIS defined in constants.h
!!
!!
!!  fluxes :  real array with space for fluxes, through one axis, 
!!            for all cells of a block and for all flux variables.
!!            fluxes(VAR, i, j, k) is VAR's flux through 
!!            the left cell face for cell i, j, k.  The 
!!            fluxes of the boundary cells of coarse blocks
!!            will have been appropriately changed, if the flux matching
!!            between fine and coarse boundaries (normally by a call to
!!            Grid_conserveFluxes), has occured since the last call
!!            of Grid_putFluxData.
!!
!!
!!  dataSize : integer array specifying the dimensions for fluxes
!!
!!             dataSize (1) holds the number of cells returned in the i direction
!!
!!             dataSize (2) holds the number of cells returned in the j direction
!!                          if 1 d problem, set datasize(2) = 1
!!
!!             dataSize (3) holds the number of cells returned in the k direction
!!                          if 1 or 2 d problem, set datasize(3) = 1
!!
!!             fluxes should contain space for fluxes of all cells in the block, 
!!             including guardcells.
!!
!!  pressureSlot : If present and greater than zero, this indicates one flux variable
!!                 in the fluxes array that may need special handling because it
!!                 really scales like a flux; normally this would be pressure,
!!                 but it could be another flux variable that the caller keeps in
!!                 flux density form.
!!
!!  areaLeft :     areas of left and right faces, only used if special scaling is
!!                 requested with the pressureSlot argument.
!!
!! NOTES 
!!
!!   Any code calling this subroutine needs to know the explicit interface,
!!   since this interface contains optional dummy arguments and assumed-shape
!!   dummy arrays. Calling FORTRAN units should therefore contain a line like
!!       use Grid_interface, ONLY: Grid_getFluxData
!!
!!   This implementation is specific to Paramesh 4.
!!
!! SEE ALSO
!!
!!   Grid_putFluxData
!!   Grid_conserveFluxes
!!   hy_sweep
!!***

subroutine Grid_getFluxData(blockID, axis, fluxes, dataSize, pressureSlot, areaLeft)

  use physicaldata, ONLY : flux_x, flux_y, flux_z, nfluxes
  use tree, ONLY : surr_blks, nodetype
  use Grid_data, ONLY : gr_xflx, gr_yflx, gr_zflx
  implicit none

#include "Flash.h"
#include "constants.h"

  integer, intent(IN) :: blockID
  integer, intent(IN) :: axis
  integer, intent(IN), dimension(3) :: dataSize
  real, intent(INOUT), dimension(NFLUXES,dataSize(1),dataSize(2),dataSize(3)) :: fluxes
  integer, intent(IN), OPTIONAL :: pressureSlot
  real, intent(IN), OPTIONAL :: areaLeft(:,:,:)

#if NFLUXES > 0
  integer :: presVar
  integer :: sx,ex,sy,ey,sz,ez

  presVar = -1
  if (present(pressureSlot)) presVar = pressureSlot

  sx = NGUARD+1
  sy = NGUARD*K2D+1
  sz = NGUARD*K3D+1
  ex = dataSize(1)-NGUARD
  ey = dataSize(2)-NGUARD*K2D
  ez = dataSize(3)-NGUARD*K3D

  select case(axis)
  case(IAXIS)
     fluxes(:,sx,sy:ey,sz:ez) = flux_x(:nfluxes,1,:,:,blockID) 
     fluxes(:,ex+1,sy:ey,sz:ez) = flux_x(:nfluxes,2,:,:,blockID)
     fluxes(:,sx+1,sy:ey,sz:ez) = gr_xflx(:,1,:,:,blockID) 
     fluxes(:,ex,sy:ey,sz:ez) = gr_xflx(:,2,:,:,blockID)
     if (presVar > 0) then
        if (.NOT.(surr_blks(1,1,1+K2D,1+K3D,blockID) > 0 .AND. &
             surr_blks(3,1,1+K2D,1+K3D,blockID) == nodetype(blockID))) then
           where (areaLeft(sx,sy:ey,sz:ez).NE.0.0) &
             fluxes(presVar,sx,sy:ey,sz:ez) = fluxes(presVar,sx,sy:ey,sz:ez) / areaLeft(sx,sy:ey,sz:ez)
        end if
        if (.NOT.(surr_blks(1,3,1+K2D,1+K3D,blockID) > 0 .AND. &
             surr_blks(3,3,1+K2D,1+K3D,blockID) == nodetype(blockID))) then
           fluxes(presVar,ex+1,sy:ey,sz:ez) = fluxes(presVar,ex+1,sy:ey,sz:ez) / areaLeft(ex+1,sy:ey,sz:ez)
        end if
     end if

  case(JAXIS)
     fluxes(:,sx:ex,sy,sz:ez) = flux_y(:nfluxes,:,1,:,blockID) 
     fluxes(:,sx:ex,ey+1,sz:ez) =flux_y(:nfluxes,:,2,:,blockID) 
     fluxes(:,sx:ex,sy+1,sz:ez) = gr_yflx(:,:,1,:,blockID) 
     fluxes(:,sx:ex,ey,sz:ez) = gr_yflx(:,:,2,:,blockID)
#if NDIM > 1
     if (presVar > 0) then
        if (.NOT.(surr_blks(1,2,1,1+K3D,blockID) > 0 .AND. &
             surr_blks(3,2,1,1+K3D,blockID) == nodetype(blockID))) then
           where (areaLeft(sx:ex,sy,sz:ez).NE.0.0) &
                fluxes(presVar,sx:ex,sy,sz:ez) = fluxes(presVar,sx:ex,sy,sz:ez) / areaLeft(sx:ex,sy,sz:ez)
        end if
        if (.NOT.(surr_blks(1,2,3,1+K3D,blockID) > 0 .AND. &
             surr_blks(3,2,3,1+K3D,blockID) == nodetype(blockID))) then
           where (areaLeft(sx:ex,ey+1,sz:ez).NE.0.0) &
             fluxes(presVar,sx:ex,ey+1,sz:ez) = fluxes(presVar,sx:ex,ey+1,sz:ez) / areaLeft(sx:ex,ey+1,sz:ez)
        end if
     end if
#endif

  case(KAXIS)
     fluxes(:,sx:ex,sy:ey,sz) = flux_z(:nfluxes,:,:,1,blockID) 
     fluxes(:,sx:ex,sy:ey,ez+1) = flux_z(:nfluxes,:,:,2,blockID) 
     fluxes(:,sx:ex,sy:ey,sz+1) = gr_zflx(:,:,:,1,blockID) 
     fluxes(:,sx:ex,sy:ey,ez) = gr_zflx(:,:,:,2,blockID) 
#if NDIM > 2
     if (presVar > 0) then
!!$        where (areaLeft(sx:ex,sy:ey,sz).NE.0.0) &  ! should not happen in any supported geometry
        if (.NOT.(surr_blks(1,2,2,1,blockID) > 0 .AND. &
             surr_blks(3,2,2,1,blockID) == nodetype(blockID))) then
           fluxes(presVar,sx:ex,sy:ey,sz) = fluxes(presVar,sx:ex,sy:ey,sz) / areaLeft(sx:ex,sy:ey,sz)
        end if
        if (.NOT.(surr_blks(1,2,2,3,blockID) > 0 .AND. &
             surr_blks(3,2,2,3,blockID) == nodetype(blockID))) then
           fluxes(presVar,sx:ex,sy:ey,ez+1) = fluxes(presVar,sx:ex,sy:ey,ez+1) / areaLeft(sx:ex,sy:ey,ez+1)
        end if
     end if
#endif
  end select
#endif

  return
end subroutine Grid_getFluxData


