!!****if* source/Grid/GridBoundaryConditions/gr_bcPutRegion
!!
!! NAME
!!  gr_bcPutRegion
!!
!! SYNOPSIS
!!
!!  gr_bcPutRegion(integer(IN)   :: gridDataStruct,
!!                 integer(IN)   :: axis,
!!                 integer(IN)   :: endPoints(LOW:HIGH,MDIM),
!!                 integer(IN)   :: regionSize(REGION_DIM),
!!                 integer(IN)   :: mask(regionSize(STRUCTSIZE)),
!!                 real(IN)      :: region(:,:,:,:),
!!                 integer(IN)   :: blockID,
!!                 integer(IN)   :: idest)
!!  
!! DESCRIPTION 
!!  This routine puts the data with boundary conditions applied back into
!!  the Grid-maintained storage for the specified Grid data structure for 
!!  all supported data structures, described for argument gridDataStruct below.
!!  The data are returned in a four-dimensional array, region, where the fourth
!!  dimension represents the individual variables of the data structure.
!!  The first three dimensions store a set of rows with relevant sections of
!!  the block on which the boundary conditions have been applied. Each row
!!  contains complete set of data points to correctly apply the boundary
!!  conditions along the specified axis. The endPoints argument specifies
!!  the bounding box of the regions being selected. For more details, see
!!  the example in gr_bcGetRegion:
!!  
!!
!! ARGUMENTS 
!!
!!
!!  gridDataStruct : optional integer value specifying data structure. 
!!                   The options are defined in constants.h and they are :
!!                   CENTER cell centered variables (default)
!!                   FACEX  face centered variable on faces along IAXIS
!!                   FACEY  face centered variable on faces along JAXIS
!!                   FACEZ  face centered variable on faces along IAXIS
!!                   WORK   work array specific to paramesh
!!  axis           : The axis on which boundary condition is being applied
!!  endPoints      : the boundaries of the region to be extracted from the 
!!                   Grid block
!!  regionSize     : regionSize(BC_DIR) contains the size of the each row;
!!                   regionSize(SECOND_DIR) contains the number of rows along the second
!!                   direction, and regionSize(THIRD_DIR) has the number of rows
!!                   along the third direction. regionSize(STRUCTSIZE) contains the
!!                   number of variables in the data structure
!!  mask           : mask if selected variables are getting boundary filled.
!!                   currently this has meaning for only PM3. Currently ignored here.
!!  region         : the extracted region
!!  blockID        : the local blockid
!!  idest          : has meaning only for PM3, where it distinguishes 
!!                   between leaf and parent nodes; see NOTES below.
!!
!!
!! NOTES
!!  Beginning with PARAMESH3: The updated solution data in the
!!  region array are not actually copied directly to "permanent" storage
!!  (UNK,WORK,etc.), but to the one-block arrays (UNK1,WORK1,etc.)
!!  that are used by PARAMESH while it is processing a block's data.
!!  Calls to gr_bcGetRegion are therefore only meaningful in certain
!!  contexts.  The idest argument is then taken as an index into these
!!  one-block arrays. It distinguishes between the slots available and
!!  must match the slot that is actually being used by PARAMESH.
!!
!! SEE ALSO
!!   gr_bcGetRegion
!!***
#include "Flash.h"

#ifdef FLASH_GRID_PARAMESH3OR4
!!REORDER(5): unk1, facevar1[xyz]
#endif
#ifdef FLASH_GRID_UG
!!REORDER(5): unk, facevar[xyz]
#endif
#ifdef FLASH_GRID_PARAMESH2
!!REORDER(5) : unk
#endif

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

subroutine gr_bcPutRegion(gridDataStruct,axis,endPoints,regionSize,mask,&
     region,blockID,idest)
  
#include "constants.h"
  
  use Driver_interface, ONLY : Driver_abortFlash
  
#ifdef FLASH_GRID_UG
  use physicaldata, ONLY: unk,facevarx,facevary,facevarz
#endif
#ifdef FLASH_GRID_PARAMESH3OR4
  use physicaldata, ONLY: unk1,facevarx1, facevary1, facevarz1
  use workspace, ONLY : work1
#endif

#ifdef FLASH_GRID_PARAMESH2
  use workspace, ONLY : work
  use physicaldata, ONLY : unk
#endif

  implicit none
  
  integer, intent(in) :: gridDataStruct,axis
  integer,dimension(LOW:HIGH,MDIM),intent(IN) :: endPoints
  integer,intent(IN) :: regionSize(REGION_DIM)
  logical,dimension(regionSize(STRUCTSIZE)),intent(IN) :: mask
  real,dimension(regionSize(BC_DIR),regionSize(SECOND_DIR),&
       regionSize(THIRD_DIR),regionSize(STRUCTSIZE)),intent(IN) :: region
  integer, intent(in) :: blockID
  integer,intent(IN) :: idest

  integer :: var,i,j,k,n,m,strt,fin
  logical :: validGridDataStruct
  integer :: varCount, bcVecEnd

  strt = endPoints(LOW,axis)
  fin  = endPoints(HIGH,axis)
  varCount = regionSize(STRUCTSIZE)
  bcVecEnd = regionSize(BC_DIR)
#ifdef DEBUG_GRID

  validGridDataStruct = .false.
  validGridDataStruct= (gridDataStruct == CENTER).or.validGridDataStruct
  validGridDataStruct= (gridDataStruct == FACEX).or.validGridDataStruct
  validGridDataStruct= (gridDataStruct == FACEY).or.validGridDataStruct
  validGridDataStruct= (gridDataStruct == FACEZ).or.validGridDataStruct
  validGridDataStruct= (gridDataStruct == WORK).or.validGridDataStruct
  
  if(.not.validGridDataStruct) then
     print *, "gr_bcPutRegion: gridDataStruct set to improper value"
     print *, "gridDataStruct must = CENTER,FACEX,FACEY,FACEZ,WORK " // &
          " (defined in constants.h)"
     call Driver_abortFlash("gr_bcPutRegion gridDataStruct must be one of CENTER,FACEX,FACEY,FACEZ,WORK(see constants.h)")
  end if

  if((gridDataStruct==WORK).and.(varCount/=1)) &
       call Driver_abortFlash("gr_bcPutRegion: varCount be 1 for work array")

  if((fin-strt+1)/=bcVecEnd)&
       call Driver_abortFlash("gr_bcPutRegion: mismatch between rowSize and the region size")
       

#endif
  
  if(axis==IAXIS) then
     do k=endPoints(LOW,KAXIS),endPoints(HIGH,KAXIS)
        m=k-endPoints(LOW,KAXIS)+1
        do j=endPoints(LOW,JAXIS),endPoints(HIGH,JAXIS)
           n=j-endPoints(LOW,JAXIS)+1
           do var=1,varCount
              select case(gridDataStruct)
#ifdef FLASH_GRID_PARAMESH3OR4
                 !! since PM3 insists on using unk1 etc
              case(CENTER)
                 unk1(var,strt:fin,j,k,idest)=region(1:bcVecEnd,n,m,var)
              case(FACEX)
                 facevarx1(var,strt:fin,j,k,idest)=region(1:bcVecEnd,n,m,var)
              case(FACEY)
                 facevary1(var,strt:fin,j,k,idest)=region(1:bcVecEnd,n,m,var)
              case(FACEZ)
                 facevarz1(var,strt:fin,j,k,idest)=region(1:bcVecEnd,n,m,var)
              case(WORK)
                 work1(strt:fin,j,k,idest)=region(1:bcVecEnd,n,m,varCount)
#else
                 !! this section in play if the grid is UG or PM2
              case(CENTER)
                 unk(var,strt:fin,j,k,blockID)=region(1:bcVecEnd,n,m,var)
#if NFACE_VARS>0
              case(FACEX)
                 facevarx(var,strt:fin,j,k,blockID)=region(1:bcVecEnd,n,m,var)
              case(FACEY)
                 facevary(var,strt:fin,j,k,blockID)=region(1:bcVecEnd,n,m,var)
              case(FACEZ)
                 facevarz(var,strt:fin,j,k,blockID)=region(1:bcVecEnd,n,m,var)
#endif
#ifdef FLASH_GRID_PARAMESH2
              case(WORK)
                 work(strt:fin,j,k,blockID,1)=region(1:bcVecEnd,n,m,varCount)
#endif
#endif
              end select
           end do
        end do
     end do

  elseif(axis==JAXIS) then
     do k=endPoints(LOW,KAXIS),endPoints(HIGH,KAXIS)
        m=k+1-endPoints(LOW,KAXIS)
        do i=endPoints(LOW,IAXIS),endPoints(HIGH,IAXIS)
           n=i+1-endPoints(LOW,IAXIS)
           do var=1,varCount
              select case(gridDataStruct)
#ifdef FLASH_GRID_PARAMESH3OR4
                 !! since PM3 insists on using unk1 etc
              case(CENTER)
                 unk1(var,i,strt:fin,k,idest)=region(1:bcVecEnd,n,m,var)
              case(FACEX)
                 facevarx1(var,i,strt:fin,k,idest)=region(1:bcVecEnd,n,m,var)
              case(FACEY)
                 facevary1(var,i,strt:fin,k,idest)=region(1:bcVecEnd,n,m,var)
              case(FACEZ)
                 facevarz1(var,i,strt:fin,k,idest)=region(1:bcVecEnd,n,m,var)
              case(WORK)
                 work1(i,strt:fin,k,idest)=region(1:bcVecEnd,n,m,varCount)
#else
                 !! this section in play if the grid is UG or PM2
              case(CENTER)
                 unk(var,i,strt:fin,k,blockID)=region(1:bcVecEnd,n,m,var)
#if NFACE_VARS>0
              case(FACEX)
                 facevarx(var,i,strt:fin,k,blockID)=region(1:bcVecEnd,n,m,var)
              case(FACEY)
                 facevary(var,i,strt:fin,k,blockID)=region(1:bcVecEnd,n,m,var)
              case(FACEZ)
                 facevarz(var,i,strt:fin,k,blockID)=region(1:bcVecEnd,n,m,var)
#endif
#ifdef FLASH_GRID_PARAMESH2
              case(WORK)
                 work(i,strt:fin,k,blockID,1)=region(1:bcVecEnd,n,m,varCount)
#endif
#endif
              end select
           end do
        end do
     end do
  elseif(axis==KAXIS) then
     do j=endPoints(LOW,JAXIS),endPoints(HIGH,JAXIS)
        m=j+1-endPoints(LOW,JAXIS)
        do i=endPoints(LOW,IAXIS),endPoints(HIGH,IAXIS)
           n=i+1-endPoints(LOW,IAXIS)
           do var=1,varCount
              select case(gridDataStruct)
#ifdef FLASH_GRID_PARAMESH3OR4
                 !! since PM3 insists on using unk1 etc
              case(CENTER)
                 unk1(var,i,j,strt:fin,idest)=region(1:bcVecEnd,n,m,var)
              case(FACEX)
                 facevarx1(var,i,j,strt:fin,idest)=region(1:bcVecEnd,n,m,var)
              case(FACEY)
                 facevary1(var,i,j,strt:fin,idest)=region(1:bcVecEnd,n,m,var)
              case(FACEZ)
                 facevarz1(var,i,j,strt:fin,idest)=region(1:bcVecEnd,n,m,var)
              case(WORK)
                 work1(i,j,strt:fin,idest)=region(1:bcVecEnd,n,m,varCount)
#else
                 !! this section in play if the grid is UG or PM2
              case(CENTER)
                 unk(var,i,j,strt:fin,blockID)=region(1:bcVecEnd,n,m,var)
#if NFACE_VARS>0
              case(FACEX)
                 facevarx(var,i,j,strt:fin,blockID)=region(1:bcVecEnd,n,m,var)
              case(FACEY)
                 facevary(var,i,j,strt:fin,blockID)=region(1:bcVecEnd,n,m,var)
              case(FACEZ)
                 facevarz(var,i,j,strt:fin,blockID)=region(1:bcVecEnd,n,m,var)
#endif
#ifdef FLASH_GRID_PARAMESH2
              case(WORK)
                 work(i,j,strt:fin,blockID,1)=region(1:bcVecEnd,n,m,varCount)
#endif
#endif
              end select
           end do
        end do
     end do
  end if
  return
end subroutine gr_bcPutRegion








