!!****if* source/physics/Flame/FlameSpeed/fl_laminarFlameSpeedBlock
!!
!! NAME
!!
!!  fl_laminarFlameSpeedBlock
!!
!! SYNOPSIS
!!
!!  fl_laminarFlameSpeedBlock(real,dimension(:,:,:),intent(IN) :: dens,
!!                            real,dimension(:,:,:),intent(IN) :: pres,
!!                            real,dimension(:,:,:),intent(OUT)  :: s,
!!                            integer,dimension(:,:) :: blklimitsgc)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   dens : density
!!
!!   pres : pressure
!!
!!   s :  flame speed
!!
!!   blklimitsgc : integer index bounds for the block including guardcells
!!
!!
!!
!!***

subroutine fl_laminarFlameSpeedBlock(dens, pres, s,blkLimitsGC)
  use Flame_data,ONLY :  fl_useConstFlameSpeed, fl_constFspeed
  use Flame_data,ONLY : II_CPOWER,II_CSPEED, II_CDENSI,II_C1,II_C2,II_C3,II_C4

  use Driver_interface, ONLY: Driver_abortFlash

  implicit none
#include "constants.h"

  real,dimension(:,:,:),intent(IN):: dens, pres
  real,dimension(:,:,:),intent(OUT) :: s
  integer,dimension(LOW:HIGH,MDIM),intent(IN):: blkLimitsGC

  ! Parameters for Shimon's fit of Laminar flame speed as function of fuel density
  ! Fit only for carbon fraction Xc = 0.5 

  real, parameter    :: c1 = -43.0e0
  real, parameter    :: c2 =  4.534e0
  real, parameter    :: c3 = -0.08333e0

#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_ILO_GC:GRID_IHI_GC,&
       GRID_JLO_GC:GRID_JHI_GC,&
       GRID_KLO_GC:GRID_KHI_GC) :: ldens
#else
  real, allocatable, dimension(:,:,:):: ldens
  integer :: sizeX, sizeY, sizeZ, istat
  sizeX=blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  sizeY=blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
  sizeZ=blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

  allocate(ldens(sizeX,sizeY,sizeZ),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot deallocate ldens in fl_laminarFlameSpeedBlock")

#endif


  if (fl_useConstFlameSpeed) then
     s = fl_constFspeed
  else
     ldens = log(dens)
     s = exp(c1 + c2*ldens + c3*(ldens**2)) 

  endif

#ifdef FIXEDBLOCKSIZE
#else
  deallocate(ldens,STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot deallocate ldens in fl_laminarFlameSpeedBlock")
#endif


end subroutine fl_laminarFlameSpeedBlock

!------------------------------------------------------------------------
