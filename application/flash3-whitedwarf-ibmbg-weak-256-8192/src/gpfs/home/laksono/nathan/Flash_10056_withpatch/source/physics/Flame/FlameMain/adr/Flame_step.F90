!!****if* source/physics/Flame/FlameMain/adr/Flame_step
!!
!! NAME
!!
!!  Flame_step
!!
!! SYNOPSIS
!!
!!  Flame_step(integer, INTENT(in)  :: num_blocks,
!!                 integer, INTENT(in), DIMENSION(num_blocks)  :: blocklist,
!!                 real, INTENT(in)  :: dt)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   num_blocks : 
!!
!!   blocklist : 
!!
!!   dt : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

#ifdef DEBUG_ALL
!#define DEBUG_INTERFACES
#endif

subroutine Flame_step( num_blocks, blockList, dt  )    

  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Flame_data, ONLY : fl_myPE,&
                         fl_shockBurning, fl_xrenormBurning,fl_meshGeom,&
                         fl_useFlame, fl_useBurn, fl_eps, &
                         fl_tnucmin,fl_tnucmax,&
                         fl_dnucmin, fl_dnucmax, fl_gcMaskSize,fl_gcMask
  use Driver_interface, ONLY : Driver_getSimTime, Driver_abortFlash
  use Eos_interface, ONLY : Eos_wrapped
  use Multispecies_interface, ONLY : Multispecies_getSumInv, &
    Multispecies_getSumFrac
  use Grid_interface, ONLY : Grid_fillGuardCells, &
    Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_getCellCoords, &
    Grid_releaseBlkPtr
  use fl_fsInterface, ONLY : fl_flameSpeed
  use fl_trackerInterface, ONLY : fl_advanceTrackerDR,fl_advanceTrackerDRCyl,&
                                   fl_advanceTrackerDRSph
  use Hydro_interface, ONLY : Hydro_detectShock
  use Flame_interface, ONLY : Flame_Effects

  implicit none

#include "constants.h"
#include "Flash.h"
    
! args these are currently not used (except for dt)
  integer, INTENT(in)                        :: num_blocks
  integer, INTENT(in), DIMENSION(num_blocks) :: blockList
  real,    INTENT(in)                        :: dt

  real :: time
  integer :: blockId
    
  real, pointer, dimension(:,:,:,:) :: solnData

  integer,dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
!.. old burn_block declarations
#ifdef FIXEDBLOCKSIZE

  real, dimension(GRID_ILO_GC:GRID_IHI_GC,&
                  GRID_JLO_GC:GRID_JHI_GC,&
                  GRID_KLO_GC:GRID_KHI_GC) :: fspeed,rpv1,rpv1dot,shock
  real,dimension(GRID_ILO_GC:GRID_IHI_GC) :: iCenter
  real,dimension(GRID_JLO_GC:GRID_JHI_GC) :: jCenter
  real,dimension(GRID_KLO_GC:GRID_KHI_GC) :: kCenter
#else
  integer :: istat
  real,allocatable,dimension(:,:,:):: fspeed,rpv1,rpv1dot,shock
  real, allocatable,dimension(:) :: iCenter, jCenter, kCenter
#endif
  integer :: iSize,jSize, kSize

  real          :: dti
  integer       :: i, j, k, n, blk
  logical       :: gcell = .true.

!==========================================


#if defined(FSPD_SCRATCH_GRID_VAR) || defined(ENUC_SCRATCH_GRID_VAR)
  ! Zero-fill the fspd and/or enuc GRIDVAR starage if it is defined - KW
  do blk = 1, num_blocks
     blockID = blockList(blk)
     call Grid_getBlkPtr(blockID, solnData, SCRATCH)
#ifdef FSPD_SCRATCH_GRID_VAR
     solnData(:,:,:, FSPD_SCRATCH_GRID_VAR) = 0.0
#endif
#ifdef ENUC_SCRATCH_GRID_VAR
     solnData(:,:,:, ENUC_SCRATCH_GRID_VAR) = 0.0
#endif
     call Grid_releaseBlkPtr(blockID, solnData, SCRATCH)
  enddo
#endif

  if( .not. fl_useFlame ) return
    
  call Driver_getSimTime(time)
  ! initialize the timers
  
  call Timers_start ("flame_step")
  
  ! make sure the guardcells are up-to-date

  call Grid_fillGuardCells(fl_myPE, CENTER, ALLDIR,&
       maskSize=fl_gcMaskSize, mask=fl_gcMask,makeMaskConsistent=.true.)
  
  ! the whole array will be loaded with data in case it is needed,
  ! so we do not need to initialize it each time

  shock = 0.e0

  dti = 1.e0/dt
    
  ! loop over all the blocks -- if we are a leaf, update progress variable

  do blk = 1, num_blocks
     blockID = blockList(blk)

     call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)
     iSize = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
     jSize = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
     kSize = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

#ifndef FIXEDBLOCKSIZE
     allocate(fspeed(iSize,jSize,kSize),STAT=istat)
      if (istat /= 0) call Driver_abortFlash("Cannot allocate fspeed in Flame_step")
     allocate(rpv1(iSize,jSize,kSize),STAT=istat)
      if (istat /= 0) call Driver_abortFlash("Cannot allocate rpv1 in Flame_step")
     allocate(rpv1dot(iSize,jSize,kSize),STAT=istat)
      if (istat /= 0) call Driver_abortFlash("Cannot allocate rpv1dot in Flame_step")
     allocate(shock(iSize,jSize,kSize),STAT=istat)
      if (istat /= 0) call Driver_abortFlash("Cannot allocate shock in Flame_step")
 
     allocate(iCenter(iSize),STAT=istat)
      if (istat /= 0) call Driver_abortFlash("Cannot allocate iCenter in Flame_step")
     allocate(jCenter(jSize),STAT=istat)
      if (istat /= 0) call Driver_abortFlash("Cannot allocate jCenter in Flame_step")
     allocate(kCenter(kSize),STAT=istat)
      if (istat /= 0) call Driver_abortFlash("Cannot allocate kCenter in Flame_step")
#endif

     call Grid_getBlkPtr(blockID, solnData, CENTER)
     !.. get reaction progress variables from database

     rpv1(:,:,:) = solnData(FLAM_MSCALAR,:,:,:)

     !.. find flame speed
     
     call fl_flameSpeed(solnData, fspeed, blockID)
     
     !.. turn off nuclear burning in shocks
     !  This isn't really used anymore
     
     if ( fl_useBurn  .and. (.not. fl_shockBurning) ) then
        call Grid_getCellCoords(IAXIS,blockID,CENTER,gcell,iCenter,iSize)
        call Grid_getCellCoords(JAXIS,blockID,CENTER,gcell,jCenter,jSize)
        call Grid_getCellCoords(KAXIS,blockID,CENTER,gcell,kCenter,kSize)
        call Hydro_detectShock(solnData,  shock, blkLimits, blkLimitsGC,&
             iCenter,jCenter,kCenter)
     endif

     !.. find time derivative of progress variable
     !   due to both reaction and diffusion
     
     select case (fl_meshGeom)
     case (CARTESIAN)
        call fl_advanceTrackerDR(fspeed, rpv1, rpv1dot, blockID,&
             blkLimits,blkLimitsGC)
     case (CYLINDRICAL)
        call fl_advanceTrackerDRcyl(fspeed, rpv1, rpv1dot, blockID,&
             blkLimits,blkLimitsGC)
     case (SPHERICAL)
        call fl_advanceTrackerDRsph(fspeed, rpv1, rpv1dot, blockID,&
             blkLimits,blkLimitsGC)
     case default
        call Driver_abortFlash("burn_block: invalid geometry")
     end select
     
     !.. update progress variable
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

           rpv1(i,j,k) = max(0.e0,min(1.e0, solnData(FLAM_MSCALAR,i,j,k) + rpv1dot(i,j,k) * dt ))
           rpv1dot(i,j,k) = ( rpv1(i,j,k) - solnData(FLAM_MSCALAR,i,j,k) ) * dti
           solnData(FLAM_MSCALAR,i,j,k) = rpv1(i,j,k)
           enddo
        enddo
     enddo

     call Flame_Effects( solnData, rpv1dot, blkLimits, blkLimitsGC , time, dt, blockID)

     ! work is done

     call Grid_releaseBlkPtr(blockID, solnData)
     


#ifndef FIXEDBLOCKSIZE
     deallocate(fspeed,STAT=istat)
       if (istat /= 0) call Driver_abortFlash("Cannot deallocate fspeed in Flame_step")
     deallocate(rpv1,STAT=istat)
       if (istat /= 0) call Driver_abortFlash("Cannot deallocate rpv1 in Flame_step")
     deallocate(rpv1dot,STAT=istat)
       if (istat /= 0) call Driver_abortFlash("Cannot deallocate rpv1dot in Flame_step")
     deallocate(shock,STAT=istat)
       if (istat /= 0) call Driver_abortFlash("Cannot deallocate shock in Flame_step")

     deallocate(iCenter,STAT=istat)
       if (istat /= 0) call Driver_abortFlash("Cannot deallocate iCenter in Flame_step")
     deallocate(jCenter,STAT=istat)
       if (istat /= 0) call Driver_abortFlash("Cannot deallocate jCenter in Flame_step")
     deallocate(kCenter,STAT=istat)
       if (istat /= 0) call Driver_abortFlash("Cannot deallocate kCenter in Flame_step")
#endif

  enddo
  
  call Timers_stop ("flame_step")
  
  return
  
end subroutine Flame_step

