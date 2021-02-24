
subroutine fl_flameSpeed(solnData, fspeed, blockID)

#include "constants.h"
#include "Flash.h"
  
  use Flame_data, ONLY :fl_useConstFlameSpeed, fl_useQuenching,&
       fl_approxAtwood, fl_meshGeom,&
       fl_constFspeed, fl_fspeedMult,&
       fl_quenchingDens0,fl_quenchingDens1,fl_quenchingDdensi,&
       fl_adrDebug
  use fl_fsData, ONLY : fl_gcdFlameSuppress, fl_gcdFlameSuppressTime, &
       fl_gcdFlameSuppressAngle

  use Driver_interface, ONLY : Driver_abortFlash, Driver_getSimTime
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, &
    Grid_releaseBlkPtr, Grid_getCellCoords
#ifdef FSPD_SCRATCH_GRID_VAR
  use Grid_interface, ONLY : Grid_getBlkPtr,Grid_releaseBlkPtr
#endif
  use fl_fsInterface, ONLY : fl_laminarFlameSpeedBlock,&
                                fl_subgridFlameSpeedCar,&
                                fl_subgridFlameSpeedCyl,&
                                fl_subgridFlameSpeedSph
  use fl_effInterface, ONLY : fl_heatReleaseBlock

  implicit none
  
  integer, intent(in) :: blockID
  real,dimension(:,:,:,:),pointer :: solnData
  real,dimension(:,:,:),intent(INOUT) ::fspeed
  
  integer :: i,j,k
  real :: time, costheta
  integer :: sizeX, sizeY, sizeZ
  integer, dimension(LOW:HIGH,MDIM)::blkLimits, blkLimitsGC
#ifdef FSPD_SCRATCH_GRID_VAR
  real, dimension(:,:,:,:), pointer :: scratchData
#endif
#ifdef FIXEDBLOCKSIZE 
  real, dimension(NSPECIES, GRID_ILO_GC:GRID_IHI_GC,&
       GRID_JLO_GC:GRID_JHI_GC,&
       GRID_KLO_GC:GRID_KHI_GC) :: mf
  
  real, dimension(GRID_ILO_GC:GRID_IHI_GC,&
                  GRID_JLO_GC:GRID_JHI_GC,&
                  GRID_KLO_GC:GRID_KHI_GC) :: atwood,factor,phi,qburn,dens,temp,gpot
  real, dimension(GRID_ILO_GC:GRID_IHI_GC) :: x, y, z
#else
  real, allocatable, dimension(:,:,:,:) :: mf
  real, allocatable, dimension(:,:,:):: atwood,factor,phi,qburn,dens,temp,gpot
  integer :: istat
  real,allocatable,dimension(:) :: x,y,z
#endif

  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

  sizeX=blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  sizeY=blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
  sizeZ=blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

#ifndef FIXEDBLOCKSIZE 
  allocate(atwood(sizeX,sizeY,sizeZ),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate atwood in fl_flameSpeed")
  allocate(factor(sizeX,sizeY,sizeZ),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate factor in fl_flameSpeed")
  allocate(phi(sizeX,sizeY,sizeZ),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate phi in fl_flameSpeed")
  allocate(qburn(sizeX,sizeY,sizeZ),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate qburn in fl_flameSpeed")
  allocate(dens(sizeX,sizeY,sizeZ),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate dens in fl_flameSpeed")
  allocate(temp(sizeX,sizeY,sizeZ),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate temp in fl_flameSpeed")
  allocate(gpot(sizeX,sizeY,sizeZ),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate gpot in fl_flameSpeed")
  allocate(x(sizeX),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate x in fl_flameSpeed")
  allocate(y(sizeY),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate y in fl_flameSpeed")
  allocate(z(sizeZ),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate z in fl_flameSpeed")
#endif

  
#ifndef FSPD_VAR
#ifndef FSPD_SCRATCH_GRID_VAR
  if (fl_adrDebug)call Driver_abortFlash&
       ('fl_flameSpeed: fl_adrDebug option needs \"fspd\" variable or GRIDVAR')
#endif
#endif
  
  
  !.. constant flame speed
  
  if (fl_useConstFlameSpeed) then
     fspeed = fl_constFspeed
#ifdef FSPD_VAR
     if (fl_adrDebug) solndata(FSPD_VAR, :,:,:) = fspeed
#endif
#ifdef FSPD_SCRATCH_GRID_VAR
     if (fl_adrDebug) then
        call Grid_getBlkPtr(blockID, scratchData, SCRATCH)
        scratchData(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
             blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
             blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS), FSPD_SCRATCH_GRID_VAR) = fspeed
        call Grid_releaseBlkPtr(blockID, scratchData, SCRATCH)
     end if
#endif
     return
  endif
  
  !.. get variables from the Grid
  !! Dev : should be possible to improve performance here
  dens = solnData(DENS_VAR,:,:,:)
  temp = solnData(TEMP_VAR,:,:,:)
#ifdef RPV1_MSCALAR
  phi  = solnData(RPV1_MSCALAR, :,:,:)
#else
  phi  = solnData(FLAM_MSCALAR, :,:,:)
#endif
  mf   = solnData(SPECIES_BEGIN:SPECIES_END,:,:,:)

!.. laminar

  call fl_heatReleaseBlock(qburn, 1)
 
  factor = solnData(GAME_VAR,:,:,:)
  factor = (factor - 1.e0)/factor
  factor = factor*qburn*dens/solnData(PRES_VAR,:,:,:)

  if (fl_approxAtwood) then
     atwood = 0.5e0*factor
     call fl_laminarFlameSpeedBlock(dens, temp, fspeed,blkLimitsGC)
  else
     !.estimate unburned density
     dens   = max(dens, dens/(1.e0-factor*phi))
     !.unburned to burned density ratio
     factor = 1.e0 + max(0.e0, factor/(1.e0 - factor*phi))
     !.atwood number
     atwood = (factor - 1.e0)/(factor + 1.e0)
     call fl_laminarFlameSpeedBlock(dens, temp, fspeed,blkLimitsGC)
     !.get local density again, since we just overwrote it
     dens   = solnData(DENS_VAR,:,:,:)
  endif

!.. if gravity present, overwrite laminar flame speed with turbulent 
#ifdef GPOT_VAR
  
  gpot(:,:,:) = solnData(GPOT_VAR,:,:,:)
  select case (fl_meshGeom)
  case (CARTESIAN)
     call fl_subgridFlameSpeedCar(atwood, gpot, phi, fspeed, blockID)
  case (CYLINDRICAL)
     call fl_subgridFlameSpeedCyl(atwood, gpot, phi, fspeed, blockID)
  case (SPHERICAL)
     call fl_subgridFlameSpeedSph(atwood, gpot, phi, fspeed, blockID)
  case default
     call Driver_abortFlash("burn_block: invalid geometry")
  end select
  
#endif
  
  ! multiply flame speed by a constant
  
  if ( fl_fspeedMult /= 1.e0 ) fspeed = fl_fspeedMult * fspeed
  
  !.. quench flame, if density falls below threshold
  
  if (fl_useQuenching) then
     
     factor = (dens - fl_quenchingDens0) * fl_quenchingDDensi
     factor = 3.e0*factor**2 - 2.e0*factor**3
     
     where ( dens < fl_quenchingDens1)
        fspeed = fspeed * factor
     endwhere
     
     where ( dens < fl_quenchingDens0)
        fspeed = 0.e0
     endwhere
     
  endif

  !  Suppress flame in selected region

  if (fl_gcdFlameSuppress) then

     call Driver_getSimTime(time)
     if (time > fl_gcdFlameSuppressTime) then
        call Grid_getCellCoords(IAXIS,blockID,CENTER,.true.,x,sizeX)
        call Grid_getCellCoords(JAXIS,blockID,CENTER,.true.,y,sizeY)
        call Grid_getCellCoords(KAXIS,blockID,CENTER,.true.,z,sizeZ)

        do k = blklimitsGC(LOW,KAXIS),blklimitsGC(HIGH,KAXIS)
           do j = blklimitsGC(LOW,JAXIS),blklimitsGC(HIGH,JAXIS)
              do i = blklimitsGC(LOW,IAXIS),blklimitsGC(HIGH,IAXIS)
                 select case (fl_meshGeom)
                 case (CARTESIAN)
                    costheta = z(k)/sqrt(x(i)**2+y(j)**2+z(k)**2)
                 case (CYLINDRICAL)
                    costheta = y(j)/sqrt(x(i)**2+y(j)**2)
                 case default
                    costheta = 1.0
                 end select
                 if (costheta <= -cos(fl_gcdFlameSuppressAngle)) then
                    fspeed(i,j,k) = 0.0
                 endif
              enddo
           enddo
        enddo

     endif
  endif
  
  !.. save flame speed for debugging
  
#ifdef FSPD_VAR
  if (fl_adrDebug)solndata(FSPD_VAR, :,:,:) = fspeed
#endif
#ifdef FSPD_SCRATCH_GRID_VAR
  if (fl_adrDebug) then
     call Grid_getBlkPtr(blockID, scratchData, SCRATCH)
        scratchData(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
             blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
             blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS), FSPD_SCRATCH_GRID_VAR) = fspeed
     call Grid_releaseBlkPtr(blockID, scratchData, SCRATCH)
  end if
#endif
  
#ifndef FIXEDBLOCKSIZE
  deallocate(atwood,STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot deallocate atwood in fl_flameSpeed")
  deallocate(factor,STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot deallocate factor in fl_flameSpeed")
  deallocate(phi,STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot deallocate phi in fl_flameSpeed")
  deallocate(qburn,STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot deallocate qburn in fl_flameSpeed")
  deallocate(dens,STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot deallocate dens in fl_flameSpeed")
  deallocate(temp,STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot deallocate temp in fl_flameSpeed")
  deallocate(gpot,STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot deallocate gpot in fl_flameSpeed")
  deallocate(x,STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot deallocate x in fl_flameSpeed")
  deallocate(y,STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot deallocate y in fl_flameSpeed")
  deallocate(z,STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot deallocate z in fl_flameSpeed")
#endif
  
  return
  
end subroutine fl_flameSpeed
