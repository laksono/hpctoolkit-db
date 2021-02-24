!!****if* source/physics/sourceTerms/Burn/BurnMain/parametricBurn/Co2nse/Burn
!!
!! NAME
!!
!!  Burn
!!
!!
!! SYNOPSIS
!!
!!   call Burn ( integer, intent(IN)    :: blockCount, 
!!               integer(:), intent(IN) :: blockList, 
!!               real, intent(IN)       ::  dt  )    
!!
!! DESCRIPTION
!!
!!  Apply parameterized, multi-stage progress variable burner 
!!  to all blocks in specified list.
!!
!! ARGUMENTS
!!
!!   blockCount -- dimension of blockList
!!   blockList -- array of blocks which should receive burning
!!   dt  --       passed to the internal bn_burner module  
!!
!! PARAMETERS
!!
!!  useBurn -- Boolean, True.  Turns on burning module
!!
!!  useShockBurn -- Boolean, FALSE.  Controls whether burning is allowed inside
!!                a regime experiencing shocks
!!
!! NOTES
!!
!! This burning unit adds the following new variables to the unk array:
!! 
!! VARIABLE    ENUC
!! MASS_SCALAR RPV1
!! MASS_SCALAR RPV2
!! MASS_SCALAR RPV3
!! MASS_SCALAR YE
!! MASS_SCALAR SUMY
!! MASS_SCALAR QBAR
!!
!!***

#define DEBUG_GRID_GCMASK

subroutine Burn ( blockCount, blockList, dt )
  
  use Grid_interface, ONLY  : Grid_fillGuardCells, &
       Grid_getBlkIndexLimits, Grid_getCellCoords, Grid_getBlkPtr, &
       Grid_releaseBlkPtr, Grid_getMyPE
  use Eos_interface, ONLY   : Eos_wrapped
  use Hydro_interface, ONLY : Hydro_detectShock
  use Burn_data
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use bn_paraInterface, ONLY : bn_paraBurn
  use Logfile_interface, ONLY : Logfile_stampVarMask
  
  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"

  !args
  integer, INTENT(in)                        :: blockCount
  integer, INTENT(in), DIMENSION(blockCount)  :: blockList
  real,    INTENT(in)                        :: dt

  ! locals
  integer                    :: i, j, k
  integer                    :: blockID, thisBlock, myPE
  real                       :: temp, dens, eint, pres
  real                       :: phi1, phi2, phi3, sumyi, ye
  real                       :: flame
  real                       :: qdot, qbar
  logical                    :: burnedZone

  integer, parameter         :: NNVARS = 11
  integer, parameter         :: nin  = NSPECIES+NNVARS 
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  logical okBurnTemp, okBurnDens, okBurnShock, okBurnNSE
  logical :: getGuardCells = .true.
  real, allocatable, dimension(:)         :: xCoord, yCoord, zCoord
  integer                                 :: xSizeCoord, ySizeCoord, zSizeCoord
  real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: shock
  real, pointer, dimension(:,:,:,:)                    :: solnData
#ifndef EINT_VAR
  real :: energyKinetic
#endif

  ! variables that are implemented in this file need values in guardcells
  logical :: gcMask(NUNK_VARS)  




  ! CHECK BURN FLAG
  if (.not. bn_useBurn) return
  

  ! START TIMERS
  call Timers_start("burn")
  

  ! UPDATE GUARDCELLS WHEN USING SHOCK DETECTION
  if (.NOT. bn_useShockBurn) then
     ! Make sure guardcells are up-to-date for Hydro_detectShock
     gcMask = .FALSE.
     gcMask(PRES_VAR) = .TRUE.
#ifdef VELX_VAR
     gcMask(VELX_VAR) = .TRUE.
#endif
#if NDIM > 1
#ifdef VELY_VAR
     gcMask(VELY_VAR) = .TRUE.
#endif
#endif
#if NDIM > 2
#ifdef VELZ_VAR
     gcMask(VELZ_VAR) = .TRUE.
#endif
#endif


     call Grid_getMyPE(myPE)
     call Grid_fillGuardCells(myPE,CENTER, ALLDIR, doEos=.false.,&
       maskSize=NUNK_VARS, mask=gcMask,makeMaskConsistent=.true.)
  endif


  ! BEING LOOP OVER BLOCKS PASSED IN     gcMask(:) = gcNeededMask(:)

  do thisBlock = 1, blockCount
     
     blockID = blockList(thisBlock)
     burnedZone = .FALSE.
     
     !GET DIMENSION AND COORD POSITIONS
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     xSizeCoord = blkLimitsGC(HIGH,IAXIS)
     ySizeCoord = blkLimitsGC(HIGH,JAXIS)
     zSizeCoord = blkLimitsGC(HIGH,KAXIS)
     allocate(xCoord(xSizeCoord))
     allocate(yCoord(ySizeCoord))
     allocate(zCoord(zSizeCoord))
     call Grid_getCellCoords(IAXIS,blockID,CENTER,getGuardCells,xCoord,xSizeCoord)
     call Grid_getCellCoords(JAXIS,blockID,CENTER,getGuardCells,yCoord,ySizeCoord)
     call Grid_getCellCoords(KAXIS,blockID,CENTER,getGuardCells,zCoord,zSizeCoord)
     
     !GET POINTER TO SOLUTION DATA
     call Grid_getBlkPtr(blockID,solnData)
     
     !SHOCK DETECTION IF REQUESTED
     if (.NOT. bn_useShockBurn) then
        call Hydro_detectShock(solnData, shock, blkLimits, blkLimitsGC, &
             xCoord,yCoord,zCoord)
     else
        shock(:,:,:) = 0
     endif
     
     !LOOP OVER CURRENT BLOCK ZONES
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              
              dens    = solnData(DENS_VAR,i,j,k)
              temp    = solnData(TEMP_VAR,i,j,k)

#ifdef EINT_VAR
              eint    = solnData(EINT_VAR,i,j,k)
#else
              energyKinetic = solnData(VELX_VAR,i,j,k)**2
#if NDIM >= 2
              energyKinetic = energyKinetic + solnData(VELY_VAR,i,j,k)**2
#endif
#if NDIM > 2
              energyKinetic = energyKinetic + solnData(VELZ_VAR,i,j,k)**2
#endif
              eint    = solnData(ENER_VAR,i,j,k) - 0.5*energyKinetic
#endif
              pres    = solnData(PRES_VAR,i,j,k)
              
#ifdef FLAM_MSCALAR
              flame   = solnData(FLAM_MSCALAR,i,j,k)
#else
              flame   = 0.e0
#endif


              phi1    = solnData(RPV1_MSCALAR,i,j,k)
              phi2    = solnData(RPV2_MSCALAR,i,j,k)
              phi3    = solnData(RPV3_MSCALAR,i,j,k)
              qdot    = 0.e0
              
              ye    = solnData(YE_MSCALAR,  i,j,k)
              sumyi = solnData(SUMY_MSCALAR,i,j,k)
              qbar  = solnData(QBAR_MSCALAR,i,j,k)
              
              ! CHECK FOR BURN CONDITIONS
              okBurnTemp   = .FALSE.
              okBurnDens   = .FALSE.
              okBurnShock  = .FALSE.
              okBurnNSE    = .FALSE.
              

              ! EXCLUDE BURNING IN SHOCK IF DESIRED
              if(shock(i,j,k) .eq. 0)      okBurnShock = .true.
              ! AVOID NSE CALCULATION AT LOW DENSITY
              !if(dens > pbquenchingDens0) okBurnNSE   = .true.
              okBurnNSE = .true.
              
              ! ADVANCE PROGRESS VARIABLES
              if (okBurnShock .and. okBurnNSE) then
                 burnedZone = .TRUE.
                 call bn_paraBurn(dens, temp, eint, pres, 0.5, &
                      phi1, phi2, phi3, flame, &
                      ye, sumyi, qbar, qdot, dt)
                 !call ii_nseBurn(dens, temp, eint, pres, &
                 !     0.5, phi1, phi2, phi3, phi1dotv,   &
                 !     ye, sumyi, qbar, qdot, time, dt)
              endif
              
              ! UPDATE SOLUTION DATA
              solnData(RPV1_MSCALAR,i,j,k) = phi1
              solnData(RPV2_MSCALAR,i,j,k) = phi2
              solnData(RPV3_MSCALAR,i,j,k) = phi3
              

              solnData(ENUC_VAR,i,j,k)     = qdot 
              solnData(ENER_VAR,i,j,k)     = solnData(ENER_VAR,i,j,k) + qdot*dt
#ifdef EINT_VAR
              solnData(EINT_VAR,i,j,k)     = solnData(EINT_VAR,i,j,k) + qdot*dt
#endif
              solnData(YE_MSCALAR,  i,j,k) = ye
              solnData(SUMY_MSCALAR,i,j,k) = sumyi
              solnData(QBAR_MSCALAR,i,j,k) = qbar
              
           enddo
           
        enddo
     enddo
     
     
     ! MAKE HYDRO CONSISTEN WITH UPDATED INTERNAL ENERGY
     if (burnedZone) then
        call Eos_wrapped(MODE_DENS_EI,blkLimits,blockID)
     end if
     
     ! RELEASE MEMORY/POINTERS
     call Grid_releaseBlkPtr(blockID,solnData)
     deallocate(xCoord)
     deallocate(yCoord)
     deallocate(zCoord)
  end do
  
  call Timers_stop("burn")
  
  return
  
end subroutine Burn
