!!****if* source/Grid/GridBoundaryConditions/gr_bcGetRegion
!!
!! NAME
!!  gr_bcGetRegion
!!
!! SYNOPSIS
!!
!!  gr_bcGetRegion(integer(IN)   :: gridDataStruct,
!!                 integer(IN)   :: axis,
!!                 integer(IN)   :: endPoints(LOW:HIGH,MDIM),
!!                 integer(IN)   :: regionSize(REGION_DIM),
!!                 integer(OUT)  :: mask(regionSize(4)),
!!                 real(out)     :: region(regionSize(1),regionSize(2),regionSize(3),regionSize(4)),
!!                 integer(IN)   :: blockID,
!!                 integer(IN)   :: idest)
!!  
!! DESCRIPTION 
!!  This routine creates a region for the application of boundary condition
!!  to all supported data structures, described for argument 
!!  gridDataStruct below.
!!  The region is stored in a four-dimensional array, where the fourth
!!  dimension represents the individual variables of the data structure.
!!  The other dimensions store a set of rows containing relevant sections of
!!  the block on which the boundary conditions are being applied. Each row
!!  contains a complete set of data points to correctly apply the boundary
!!  conditions along the specified axis. The endPoints argument specifies
!!  the bounding box of the regions being selected. For more details, see
!!  the example below:
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
!!  regionSize     : regionSize(BC_DIR) contains the size of the each row
!!                   regionSize(SECOND_DIR) contains the number of rows along the
!!                   second direction, and regionSize(THIRD_DIR) has the number of rows
!!                   along the third direction. regionSize(STRUCTSIZE) contains the
!!                   number of variables in the data structure
!!  mask           : Mask to be used if selected variables are getting boundary
!!                   filled. Currently this has meaning for only PM3 and PM4.
!!  region         : the extracted region
!!  blockID        : the local block ID.
!!                   (With Paramesh3 f. this may actually be a blockHandle that refers
!!                   to a remote block, but this implementation does not actually use
!!                   the blockID at all if the grid is Paramesh3 f. - see idest instead).
!!  idest          : has meaning only for PM3 and PM4, where it distinguishes between 
!!                   leaf and parent nodes, should be 1 or 2; see NOTES below.
!!
!!
!! EXAMPLE 
!!   In this example with 2D data on a LEAF block, 
!!   we want to apply boundary conditions on the right face of IAXIS, 
!!   and we wish to fetch columns 7 and 8 of the interior
!!   data and all the columns of the guardcell data. Since this example has
!!   4 guardcells on each side, the actual column numbers are 11,12 for the
!!   interior and 13-16 for the guardcells to be filled.
!!
!!       ---- - - - - - - - - ----
!!       ---- - - - - - - - - ----
!!       ---- - - - - - - - - ----
!!       ---- - - - - - - - - ----
!!     8 ----|-|-|-|-|-|-|*|*|----
!!     7 ----|-|-|-|-|-|-|*|*|----
!!     6 ----|-|-|-|-|-|-|*|*|----
!!     5 ----|-|-|-|-|-|-|*|*|----
!!     4 ----|-|-|-|-|-|-|*|*|----
!!     3 ----|-|-|-|-|-|-|*|*|----
!!     2 ----|-|-|-|-|-|-|*|*|----
!!     1 ----|-|-|-|-|-|-|*|*|----
!!     j ---- - - - - - - - - ----
!!       ---- - - - - - - - - ----
!!       ---- - - - - - - - - ----
!!       ---- - - - - - - - - ----
!!         i  1-2-3-4 5-6-7-8 
!!     
!!     Then the values in the arguement endpoints should be
!!
!!          endPoint(LOW,IAXIS) = 11; endPoint(HIGH,IAXIS)=16
!!          endPoint(LOW,JAXIS) = 1 ; endPoint(HIGH,JAXIS)=16
!!          endPoint(LOW:HIGH,KAXIS) = 1
!!
!!     RegionSize argument should have 
!!         RegionSize(:)=endPoint(HIGH,:)-endPoint(LOW,:)+1
!!
!!     The argument Region will contain the data being fetched, so
!!     should be allocated as 
!!     Region(RegionSize(IAXIS),RegionSize(JAXIS),RegionSize(KAXIS),vars)
!!     where vars is the number of variables in the data structure;
!!     NUNK_VARS for cell centered, NFACE_VARS for face centered along IAXIS etc.
!!
!!     Please Note that if we were interested in rows (in other words
!!     the top regions) then the allocation would have been
!!     Region(RegionSize(JAXIS),RegionSize(IAXIS),RegionSize(KAXIS),vars)
!!
!!     The call will have the following syntax:
!!
!!     call gr_bcGetRegion(CENTER,IAXIS,endPoint,regionSize,mask,Region,blockID,
!!     idest)
!!
!! NOTES
!!  Beginning with PARAMESH3: The solution data used to fill the
!!  region array are not copied directly from "permanent" storage
!!  (UNK,WORK,etc.), but from the one-block arrays (UNK1,WORK1,etc.)
!!  that are filled by PARAMESH while it is processing a block's data.
!!  Calls to gr_bcGetRegion are therefore only valid in certain
!!  contexts.  The idest argument is then taken as an index into these
!!  one-block arrays. It distinguishes between the slots available and
!!  must match the slot that has actually been filled by PARAMESH.
!!  
!!
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

subroutine gr_bcGetRegion(gridDataStruct,axis,endPoints,regionSize,mask,&
     region,blockID,idest)
  
#include "constants.h"
  
  use Driver_interface, ONLY : Driver_abortFlash
  
#ifdef FLASH_GRID_UG
  use physicaldata, ONLY: unk,facevarx,facevary,facevarz
#endif
#ifdef FLASH_GRID_PARAMESH3OR4
  use physicaldata, ONLY: unk1,facevarx1, facevary1, facevarz1, gcell_on_cc, gcell_on_fc
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
  logical,dimension(regionSize(STRUCTSIZE)),intent(OUT) :: mask
  real,dimension(regionSize(BC_DIR),regionSize(SECOND_DIR),&
       regionSize(THIRD_DIR),regionSize(STRUCTSIZE)),intent(OUT) :: region
  integer, intent(in) :: blockID
  integer,intent(IN) :: idest

  integer :: var,i,j,k,n,m,strt,fin, varCount,bcVecEnd
  logical :: validGridDataStruct


  strt = endPoints(LOW,axis)
  fin  = endPoints(HIGH,axis)
  varCount=regionSize(STRUCTSIZE)
  bcVecEnd=regionSize(BC_DIR)

#ifdef DEBUG_GRID

  validGridDataStruct = .false.
  validGridDataStruct= (gridDataStruct == CENTER).or.validGridDataStruct
  validGridDataStruct= (gridDataStruct == FACEX).or.validGridDataStruct
  validGridDataStruct= (gridDataStruct == FACEY).or.validGridDataStruct
  validGridDataStruct= (gridDataStruct == FACEZ).or.validGridDataStruct
  validGridDataStruct= (gridDataStruct == WORK).or.validGridDataStruct
  
  if(.not.validGridDataStruct) then
     print *, "gr_bcGetRegion: gridDataStruct set to improper value"
     print *, "gridDataStruct must = CENTER,FACEX,FACEY,FACEZ,WORK " // &
          " (defined in constants.h)"
     call Driver_abortFlash("gr_bcGetRegion gridDataStruct must be one of CENTER,FACEX,FACEY,FACEZ,WORK(see constants.h)")
  end if

  if((gridDataStruct==WORK).and.(varCount/=1)) &
       call Driver_abortFlash("gr_bcGetRegion: varCount be 1 for work array")

  if((fin-strt+1)/=bcVecEnd)&
       call Driver_abortFlash("gr_bcGetRegion: mismatch between rowSize and the region size")
       

#endif
  
  mask=.true.
#ifdef FLASH_GRID_PARAMESH3OR4
  if (gridDataStruct==CENTER) then
     mask(1:varCount)=gcell_on_cc(1:varCount)
  elseif(gridDataStruct==FACEX) then
     mask(1:varCount)=gcell_on_fc(1,1:varCount)
  elseif(gridDataStruct==FACEY) then
     mask(1:varCount)=gcell_on_fc(2,1:varCount)
  elseif(gridDataStruct==FACEZ) then
     mask(1:varCount)=gcell_on_fc(3,1:varCount)
  end if
#endif
  if(axis==IAXIS) then
     do k=endPoints(LOW,KAXIS),endPoints(HIGH,KAXIS)
        m=k-endPoints(LOW,KAXIS)+1
        do j=endPoints(LOW,JAXIS),endPoints(HIGH,JAXIS)
           n=j-endPoints(LOW,JAXIS)+1
           do var=1,varCount
              select case(gridDataStruct)
#ifdef FLASH_GRID_PARAMESH3OR4
                 !! since PM3 f. insists on using unk1 etc
              case(CENTER)
                 region(1:bcVecEnd,n,m,var)=unk1(var,strt:fin,j,k,idest)
              case(FACEX)
                 region(1:bcVecEnd,n,m,var)=facevarx1(var,strt:fin,j,k,idest)
              case(FACEY)
                 region(1:bcVecEnd,n,m,var)=facevary1(var,strt:fin,j,k,idest)
              case(FACEZ)
                 region(1:bcVecEnd,n,m,var)=facevarz1(var,strt:fin,j,k,idest)
              case(WORK)
                 region(1:bcVecEnd,n,m,varCount)=work1(strt:fin,j,k,idest)
#else
                 !! this section in play if the grid is UG or PM2
              case(CENTER)
                 region(1:bcVecEnd,n,m,var)=unk(var,strt:fin,j,k,blockID)
#if NFACE_VARS>0
              case(FACEX)
                 region(1:bcVecEnd,n,m,var)=facevarx(var,strt:fin,j,k,blockID)
              case(FACEY)
                 region(1:bcVecEnd,n,m,var)=facevary(var,strt:fin,j,k,blockID)
              case(FACEZ)
                 region(1:bcVecEnd,n,m,var)=facevarz(var,strt:fin,j,k,blockID)
#endif
#ifdef FLASH_GRID_PARAMESH2
              case(WORK)
                 region(1:bcVecEnd,n,m,varCount)=work(strt:fin,j,k,blockID,1)
#endif
#endif
              end select
           end do
        end do
     end do

  elseif(axis==JAXIS) then
     do k=endPoints(LOW,KAXIS),endPoints(HIGH,KAXIS)
        m=k-endPoints(LOW,KAXIS)+1
        do i=endPoints(LOW,IAXIS),endPoints(HIGH,IAXIS)
           n=i-endPoints(LOW,IAXIS)+1
           do var=1,varCount
              select case(gridDataStruct)
#ifdef FLASH_GRID_PARAMESH3OR4
                 !! since PM3 f. insists on using unk1 etc
              case(CENTER)
                 region(1:bcVecEnd,n,m,var)=unk1(var,i,strt:fin,k,idest)
              case(FACEX)
                 region(1:bcVecEnd,n,m,var)=facevarx1(var,i,strt:fin,k,idest)
              case(FACEY)
                 region(1:bcVecEnd,n,m,var)=facevary1(var,i,strt:fin,k,idest)
              case(FACEZ)
                 region(1:bcVecEnd,n,m,var)=facevarz1(var,i,strt:fin,k,idest)
              case(WORK)
                 region(1:bcVecEnd,n,m,varCount)=work1(i,strt:fin,k,idest)
#else
                 !! this section in play if the grid is UG or PM2
              case(CENTER)
                 region(1:bcVecEnd,n,m,var)=unk(var,i,strt:fin,k,blockID)
#if NFACE_VARS>0
              case(FACEX)
                 region(1:bcVecEnd,n,m,var)=facevarx(var,i,strt:fin,k,blockID)
              case(FACEY)
                 region(1:bcVecEnd,n,m,var)=facevary(var,i,strt:fin,k,blockID)
              case(FACEZ)
                 region(1:bcVecEnd,n,m,var)=facevarz(var,i,strt:fin,k,blockID)
#endif
#ifdef FLASH_GRID_PARAMESH2
              case(WORK)
                 region(1:bcVecEnd,n,m,varCount)=work(i,strt:fin,k,blockID,1)
#endif
#endif
              end select
           end do
        end do
     end do
  elseif(axis==KAXIS) then
     do j=endPoints(LOW,JAXIS),endPoints(HIGH,JAXIS)
        m=j-endPoints(LOW,JAXIS)+1
        do i=endPoints(LOW,IAXIS),endPoints(HIGH,IAXIS)
           n=i-endPoints(LOW,IAXIS)+1
           do var=1,varCount
              select case(gridDataStruct)
#ifdef FLASH_GRID_PARAMESH3OR4
                 !! since PM3 f. insists on using unk1 etc
              case(CENTER)
                 region(1:bcVecEnd,n,m,var)=unk1(var,i,j,strt:fin,idest)
              case(FACEX)
                 region(1:bcVecEnd,n,m,var)=facevarx1(var,i,j,strt:fin,idest)
              case(FACEY)
                 region(1:bcVecEnd,n,m,var)=facevary1(var,i,j,strt:fin,idest)
              case(FACEZ)
                 region(1:bcVecEnd,n,m,var)=facevarz1(var,i,j,strt:fin,idest)
              case(WORK)
                 region(1:bcVecEnd,n,m,varCount)=work1(i,j,strt:fin,idest)
#else
                 !! this section in play if the grid is UG or PM2
              case(CENTER)
                 region(1:bcVecEnd,n,m,var)=unk(var,i,j,strt:fin,blockID)
#if NFACE_VARS>0
              case(FACEX)
                 region(1:bcVecEnd,n,m,var)=facevarx(var,i,j,strt:fin,blockID)
              case(FACEY)
                 region(1:bcVecEnd,n,m,var)=facevary(var,i,j,strt:fin,blockID)
              case(FACEZ)
                 region(1:bcVecEnd,n,m,var)=facevarz(var,i,j,strt:fin,blockID)
#endif
#ifdef FLASH_GRID_PARAMESH2
              case(WORK)
                 region(1:bcVecEnd,n,m,varCount)=work(i,j,strt:fin,blockID,1)
#endif
#endif
              end select
           end do
        end do
     end do
  end if
  return
end subroutine gr_bcGetRegion








