!!****if* source/Simulation/SimulationMain/WD_def/varminmax
!!
!!  NAME
!!   varminmax
!!
!! SYNOPSIS
!!  varminmax
!!        ( real(in)    :: varData,
!!          real(out)   :: varmin, 
!!          real(out)   :: varmax,
!!          real(in)    :: blkLimitsGC)
!!
!!  DESCRIPTION
!!
!!
!!  ARGUMENTS
!!
!!***
!       

!
subroutine varminmax(varData, varmin, varmax, blkLimitsGC)
  
  implicit none
#include "constants.h"
#include "Flash.h"
  
  !..DECLAR LOCAL VARIABLES
  integer :: i, j, k, ii, jj, kk
  integer, intent(IN),dimension(2,MDIM)::blkLimitsGC
  real, intent(IN), dimension(GRID_ILO_GC:GRID_IHI_GC,GRID_JLO_GC:GRID_JHI_GC,GRID_KLO_GC:GRID_KHI_GC) ::  &
       varData
  real, intent(OUT),dimension(GRID_ILO_GC:GRID_IHI_GC,GRID_JLO_GC:GRID_JHI_GC,GRID_KLO_GC:GRID_KHI_GC) ::  &
       varmin, varmax
  
  do k = blkLimitsGC(LOW,KAXIS)+1*K3D, blkLimitsGC(HIGH,KAXIS)-1*K3D
     do j = blkLimitsGC(LOW,JAXIS)+1*K2D, blkLimitsGC(HIGH,JAXIS)-1*K2D
        do i = blkLimitsGC(LOW,IAXIS)+1, blkLimitsGC(HIGH,IAXIS)-1
           varmin(i,j,k) = varData(i,j,k)
           varmax(i,j,k) = varData(i,j,k)
        end do
     end do
  end do
  
  do kk = -K3D,K3D
     do jj = -K2D,+K2D
        do ii = -1,+1
           do k = blkLimitsGC(LOW,KAXIS)+1*K3D, blkLimitsGC(HIGH,KAXIS)-1*K3D
              do j = blkLimitsGC(LOW,JAXIS)+1*K2D, blkLimitsGC(HIGH,JAXIS)-1*K2D
                 do i = blkLimitsGC(LOW,IAXIS)+1, blkLimitsGC(HIGH,IAXIS)-1
                    varmin(i,j,k) = min(varmin(i,j,k),varData(i+ii,j+jj,k+kk))
                    varmax(i,j,k) = max(varmax(i,j,k),varData(i+ii,j+jj,k+kk))
                 end do
              end do
           end do
        end do
     end do
  end do
  
  return
  
end subroutine varminmax


