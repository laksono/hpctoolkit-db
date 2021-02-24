
subroutine fl_subgridFlameSpeedCar(atwood, gpot, phi, flamespeed, blockID)

  use Flame_data, ONLY : fl_minCellSize, fl_s1, fl_s2, &
    fl_subgridSuppress, fl_subgridSuppressTime, fl_subgridSuppressAngle
  use Grid_interface, ONLY : Grid_getDeltas, Grid_getBlkIndexLimits, &
    Grid_getCellCoords
  use Driver_interface, ONLY:  Driver_abortFlash, Driver_getSimTime
  implicit none
#include "constants.h"
#include "Flash.h"
      
  integer, intent(in) :: blockID

  real, dimension(:,:,:),intent(IN) :: atwood,gpot,phi
  real, dimension(:,:,:),intent(inout) :: flamespeed

  integer             :: i, j, k
  real                :: f, df, dg, dfdg, dgdg, modg, ndotg
  real                :: deltas(MDIM), dx, over_twelve_dx, over_sq_twelve_dx 
  real                :: subgrid
  integer,dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  integer             :: sizeX, sizeY, sizeZ
  integer             :: istat
  real                :: time
!==============================================================================

#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_ILO_GC:GRID_IHI_GC) :: x
  real, dimension(GRID_JLO_GC:GRID_JHI_GC) :: y
  real, dimension(GRID_KLO_GC:GRID_KHI_GC) :: z
#else
  real,allocatable,dimension(:) :: x,y,z
#endif

  call Grid_getBlkIndexLimits(blockID, blkLimits,blkLimitsGC)
  sizeX=blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  sizeY=blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
  sizeZ=blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

#ifndef FIXEDBLOCKSIZE
  allocate(x(sizeX),STAT=istat)
  if (istat /= 0) call Driver_abortFlash("Cannot allocate x in fl_subgridFlameSpeedCar")
  allocate(y(sizeY),STAT=istat)
  if (istat /= 0) call Driver_abortFlash("Cannot allocate y in fl_subgridFlameSpeedCar")
  allocate(z(sizeZ),STAT=istat)
  if (istat /= 0) call Driver_abortFlash("Cannot allocate z in fl_subgridFlameSpeedCar")
#endif

  call Driver_getSimTime(time)

  call Grid_getCellCoords(IAXIS,blockID,CENTER,.true.,x,sizeX)
  call Grid_getCellCoords(JAXIS,blockID,CENTER,.true.,y,sizeY)
  call Grid_getCellCoords(KAXIS,blockID,CENTER,.true.,z,sizeZ)

!.. for now, assume dx = dy = dz = const across the block

  call Grid_getDeltas(blockID, deltas)
  dx = deltas(1)
  over_twelve_dx    = 1.e0/(12.e0*dx)
  over_sq_twelve_dx = over_twelve_dx**2


  do k = blklimits(LOW,KAXIS)-2*K3D,blklimits(HIGH,KAXIS)+2*K3D
      do j = blklimits(LOW,JAXIS)-2*K2D,blklimits(HIGH,JAXIS)+2*K2D
        do i = blklimits(LOW,IAXIS)-2,blklimits(HIGH,IAXIS)+2

           df =  phi(i-2,j,k) - 8.e0*phi(i-1,j,k) + &
                8.e0*phi(i+1,j,k) - phi(i+2,j,k) 
           
           dg = gpot(i-2,j,k) - 8.e0*gpot(i-1,j,k) + &
                8.e0*gpot(i+1,j,k) - gpot(i+2,j,k) 
           
           dfdg = df*dg
           dgdg = dg*dg
           
           if (NDIM >= 2) then
              
              df =  phi(i,j-2,k) - 8.e0*phi(i,j-1,k) + &
                    8.e0*phi(i,j+1,k)  - phi(i,j+2,k) 

              dg = gpot(i,j-2,k) - 8.e0*gpot(i,j-1,k) + &
                                   8.e0*gpot(i,j+1,k) - gpot(i,j+2,k) 

              dfdg = dfdg + df*dg
              dgdg = dgdg + dg*dg
              
           endif
           
           if (NDIM == 3) then
              
              df =  phi(i,j,k-2) - 8.e0*phi(i,j,k-1) + &
                   8.e0*phi(i,j,k+1)  - phi(i,j,k+2) 
              
              dg = gpot(i,j,k-2) - 8.e0*gpot(i,j,k-1) + &
                   8.e0*gpot(i,j,k+1) - gpot(i,j,k+2) 
              
              dfdg = dfdg + df*dg
              dgdg = dgdg + dg*dg
              
           endif
           
           modg  = sqrt(dgdg) * over_twelve_dx 
           
           f = phi(i,j,k)
           
           if (f <= 0.e0 .or. f >= 1.e0) then
              ndotg = 0.e0
           else
              ndotg = dfdg * over_sq_twelve_dx / f
           endif
           
           subgrid = sqrt( atwood(i,j,k)*( fl_s1*modg + fl_s2*max(0.e0, ndotg) ) )

           if (fl_subgridSuppress) then
              if (NDIM /= 3) call Driver_abortFlash("subgrid suppression only valid for cartesian in 3D")
              if (.not. (time > fl_subgridSuppressTime .and. &
                    z(k)/sqrt(x(i)**2+y(j)**2+z(k)**2) <= -cos(fl_subgridSuppressAngle)) ) then
                 flamespeed(i,j,k) = max(flamespeed(i,j,k), subgrid)
              endif
           else
              flamespeed(i,j,k) = max(flamespeed(i,j,k), subgrid)
           endif
           
        enddo
     enddo
  enddo

#ifndef FIXEDBLOCKSIZE
  deallocate(x,STAT=istat)
  if (istat /= 0) call Driver_abortFlash("Cannot deallocate x in fl_subgridFlameSpeedCar")
  deallocate(y,STAT=istat)
  if (istat /= 0) call Driver_abortFlash("Cannot deallocate y in fl_subgridFlameSpeedCar")
  deallocate(z,STAT=istat)
  if (istat /= 0) call Driver_abortFlash("Cannot deallocate z in fl_subgridFlameSpeedCar")
#endif

 
  !------------------------------------------------------------------
  
end subroutine fl_subgridFlameSpeedCar


