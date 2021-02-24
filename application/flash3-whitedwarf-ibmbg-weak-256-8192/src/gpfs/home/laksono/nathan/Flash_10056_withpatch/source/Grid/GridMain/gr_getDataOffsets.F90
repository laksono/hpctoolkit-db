!!****if* source/Grid/GridMain/gr_getDataOffsets
!!
!! NAME
!!  gr_getDataOffsets
!!
!! SYNOPSIS
!!
!!  gr_getDataOffsets(integer(IN) :: blockID,
!!                    integer(IN) :: gridDataStruct,
!!                    integer(IN) :: beginCount, 
!!                    integer(IN) :: startingPos(MDIM),
!!                    integer(IN) :: length(MDIM),
!!                    integer(OUT):: begOffset(MDIM)
!!                    logical(OUT):: getBlkPtr
!!  
!! DESCRIPTION 
!!  
!!  This routine determines the offset of the data withing the block from
!!  where to fetch. These offsets are used by data get/put routines to 
!!  correctly determine the starting index in all dimensions
!!
!! ARGUMENTS 
!!
!!  blockID   : my block number
!!
!!  gridDataStruct : integer value specifying the type of data desired.
!!
!!  beginCount : tells the routine where to start index counting.  beginCount can
!!               be set to INTERIOR or EXTERIOR.  If INTERIOR is specified
!!               guardcell indices are not included and index 1 is the first interior cell. 
!!               If EXTERIOR is specified
!!               the first index, 1, is the left most guardcell.  See examples
!!               in get/put data routines for more details
!!
!!  startingPos(MDIM):
!!           specifies the starting position in each dimension of 
!!           the plane of data being fetched.
!!   
!!           startingPos(1) = i
!!           startingPos(2) = j
!!           startingPos(3) = k
!!
!!           If a problem is 2 dimensions startingPos(3) is irrelevant and
!!           ignored.  If a problem is only 1 dimension this routine doesn't
!!           make any sense and an error is returned.
!!
!! length : the length of the data in each dimension
!!
!! begOffset - the calculated offset values
!! getBlkPtr - indicates to the data get/put routines whether blk pointer needs to 
!!          fetched. It is always false in permanent GC mode, and in UG. Depending
!!          upon the starting index, the block may or may not need to be fetched in
!!          the non permanent mode when using PM3
!!
!!***
subroutine gr_getDataOffsets(blockID,gridDataStruct,startingPos,length,beginCount,begOffset,getBlkPtr)
  use Grid_interface, ONLY : Grid_getBlkIndexLimits
#include "constants.h"
#include "Flash.h"
  implicit none
  integer, intent(IN)  :: blockID
  integer, intent(IN)  :: gridDataStruct
  integer, intent(IN)  :: beginCount
  integer,dimension(MDIM),intent(IN) :: startingPos
  integer,dimension(MDIM),intent(IN) :: length
  integer,dimension(MDIM),intent(OUT) :: begOffset
  logical,intent(OUT) :: getBlkPtr

  integer,dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  integer :: i

  begOffset = 0
  !! We fetch data from the main data structure (unk etc in most cases) unless
  !! in non-permanent mode, when some of the guardcell need to be fetched too
  
  getBlkPtr=.false.
  
#ifdef FL_NON_PERMANENT_GUARDCELLS 
  !! If running with PM3 and operating the non permanent guardcell mode, and
  !! the grid data structure is not SCRATCH
  if((gridDataStruct == CENTER).or.(gridDataStruct == FACEX).or.(gridDataStruct ==FACEY).or.&
       (gridDataStruct == FACEZ)) then
     if(beginCount == EXTERIOR) then
        call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC,gridDataStruct)
        do i = 1,NDIM
           getBlkPtr = getBlkPtr.or.(startingPos(i) < blkLimits(LOW,i))
           getBlkPtr = getBlkPtr.or.((startingPos(i)+length(i)-1)>blkLimits(HIGH,i))
        end do
        !! if all the data is from the interior of the block, we don't need to
        !! form the block, just subtracting the number of guardcells from the starting
        !! position gives the right offset.
        if(.not.getBlkPtr) begOffset(1:NDIM) = blkLimitsGC(LOW,1:NDIM)-blkLimits(LOW,1:NDIM)
     end if !! If the beginCount is interior, then nothing needs to be done to the offset
  end if
#else
  !! If operating with permanent guardcells, then the interior implies offset by gcell
  if(beginCount == INTERIOR) then
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC,gridDataStruct)
     begOffset(1:NDIM) = blkLimits(LOW,1:NDIM)-1
  end if
#endif  
  return
end subroutine gr_getDataOffsets
