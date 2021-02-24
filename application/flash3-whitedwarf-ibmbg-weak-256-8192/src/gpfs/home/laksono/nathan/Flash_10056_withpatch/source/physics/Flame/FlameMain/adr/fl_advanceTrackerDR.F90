
subroutine fl_advanceTrackerDR(flamespeed, phi, phidot,blockID,&
                                         blkLimits,blkLimitsGC)

  use Flame_data, ONLY :fl_adrLimitPhi 
  use Driver_interface, ONLY : Driver_getDt, Driver_abortFlash
  use Grid_interface, ONLY : Grid_getDeltas
  use fl_adrInterface, ONLY : fl_adrReactionRateBlock,fl_adrDiffusivityBlock
 
  implicit none
#include "Flash.h"
#include "constants.h"
      
  integer,intent(IN) :: blockID
  real, dimension(:,:,:),intent(INOUT) :: flamespeed,phidot

  real, dimension(:,:,:),intent(INOUT) :: phi
  integer,dimension(LOW:HIGH,MDIM),intent(IN) ::blkLimits,blkLimitsGC

  real, dimension(MDIM) :: delta

#ifdef FIXEDBLOCKSIZE

  real, dimension(GRID_ILO_GC:GRID_IHI_GC, &
       GRID_JLO_GC:GRID_JHI_GC, GRID_KLO_GC:GRID_KHI_GC) :: phidotrr,kappa

#else
  integer :: istat
  real,allocatable,dimension(:,:,:)::phidotrr,kappa
#endif
  
  integer       :: i, j, k,sizeX, sizeY,sizeZ
  real          :: dx, df, ddf, dt
  real          :: over_twelve_sq_dx

!==============================================================================

!------------------------------------------------------------------------------

  sizeX=blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  sizeY=blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
  sizeZ=blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

#ifndef FIXEDBLOCKSIZE
  allocate(phidotrr(sizeX,sizeY,sizeZ),STAT=istat) 
     if (istat /= 0) call Driver_abortFlash("Cannot allocate phidotrr in fl_advanceTrackerDR")
 allocate(kappa(sizeX,sizeY,sizeZ),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate kappa in fl_advanceTrackerDR")
#endif

  call Driver_getDt(dt)
  call Grid_getDeltas(blockID,delta)
!.. for now, assume dx = dy = dz = const across the block

  dx = delta(IAXIS)
  over_twelve_sq_dx = 1.e0/(12.e0*dx**2)

!.. compute reaction term
  call fl_adrReactionRateBlock(flamespeed, phi, phidotrr)

  phi = phi + phidotrr*dt

!.. add diffusion term
  call fl_adrDiffusivityBlock(flamespeed, kappa)

  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

           ddf =    - phi(i-2,j,k) + 16.e0*phi(i-1,j,k) - &
                30.e0*phi(i,j,k)   + 16.e0*phi(i+1,j,k) - phi(i+2,j,k) 

           if (NDIM >= 2) then
              ddf =  ddf - phi(i,j-2,k) + 16.e0*phi(i,j-1,k) - &
                     30.e0*phi(i,j,k)   + 16.e0*phi(i,j+1,k) - phi(i,j+2,k) 
           endif

           if (NDIM == 3) then
              ddf =  ddf - phi(i,j,k-2) + 16.e0*phi(i,j,k-1) - &
                     30.e0*phi(i,j,k)   + 16.e0*phi(i,j,k+1) - phi(i,j,k+2) 
           endif

           phidot(i,j,k) = kappa(i,j,k) * ddf*over_twelve_sq_dx

           phidot(i,j,k) = phidot(i,j,k) + phidotrr(i,j,k)

        enddo
     enddo
  enddo

  phi = phi - phidotrr*dt
  
!------------------------------------------------------------------------------
#ifndef FIXEDBLOCKSIZE
  deallocate(phidotrr,STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot deallocate phidotrr in fl_advanceTrackerDR")
  deallocate(kappa,STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot deallocate kappa in fl_advanceTrackerDR")
#endif

end subroutine fl_advanceTrackerDR
