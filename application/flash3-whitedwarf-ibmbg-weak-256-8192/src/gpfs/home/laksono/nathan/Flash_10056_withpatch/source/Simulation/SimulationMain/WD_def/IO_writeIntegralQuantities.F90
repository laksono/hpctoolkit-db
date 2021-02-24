!!****if* source/Simulation/SimulationMain/WD_def/IO_writeIntegralQuantities
!!
!!
!!  NAME
!!    IO_writeIntegralQuantities
!!
!!  SYNOPSIS
!!    call IO_writeIntegralQuantities(integer(in) :: myPE, 
!!                                    integer(in) :: isFirst,
!!                                    real(in)    :: simTime)
!!
!!  DESCRIPTION
!!
!!   Compute the values of integral quantities (eg. total energy)
!!   and write them to an ASCII file.  If this is the initial step,
!!   create the file and write a header to it before writing the data.
!!
!!   Presently, this supports 1, 2, and 3-d Cartesian geometry and 2-d
!!   cylindrical geometry (r,z).  More geometries can be added by
!!   modifying the volume of each zone (dvol).
!!
!!   Users should modify this routine if they want to store any
!!   quantities other than default values in the flash.dat.  Make sure
!!   to modify the nGlobalSum parameter to match the number of
!!   quantities written.  Also make sure to modify the header to match
!!   the names of quantities with those calculated in the lsum and
!!   gsum arrays.
!!  
!!  ARGUMENTS
!!    
!!   myPE - current processor number
!!   isFirst - if 1 then write header info plus data, otherwise just write data
!!   simTime - simulation time
!!
!!
!!***

!!REORDER(4):solnData

subroutine IO_writeIntegralQuantities (myPE, isFirst, simTime)

  use Simulation_data, ONLY: sim_gridGeom, sim_restart
  use fl_adrData, ONLY : fl_adrEps
  use Hydro_data, ONLY: hy_gravMassXYZ, hy_gravMassZYX
  use Grid_interface, ONLY : Grid_getListOfBlocks, &
    Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_getSingleCellVol, &
    Grid_releaseBlkPtr, Grid_getCellCoords
  use fl_adrInterface, ONLY : fl_adrReactionRate


  implicit none

#include "mpif.h"
#include "constants.h"
#include "Flash.h"
  
  integer, intent(in) :: myPE
  real, intent(in) :: simTime

  integer, intent(in) :: isFirst

  integer :: lb, count
  
  integer :: funit = 99
  integer :: error
  
  character (len=MAX_STRING_LENGTH), save :: statsfilename = 'flash.dat'
  character (len=MAX_STRING_LENGTH), save :: fname 

  integer :: blockList(MAXBLOCKS)

  integer :: blkLimits(HIGH, MDIM), blkLimitsGC(HIGH, MDIM)

  integer, parameter ::  nGlobalSum = 36          ! Number of globally-summed quantities
  real :: gsum(nGlobalSum) !Global summed quantities
  real :: lsum(nGlobalSum) !Global summed quantities
  real :: gmax  ! Added by Cal to hold Maximum quantites from MPI_REDUCE command. Just one quantity at the moment
  integer :: isizeGC, jsizeGC, ksizeGC
  integer :: i, j, k
  real :: dvol, dm, hdm, dgrav, dedis, s, phi, phidot, &
       eps, momx, momy, momz, ekin, phi_dm, nse_threshold

  real, DIMENSION(:,:,:,:), POINTER :: solnData
  
!!  real,    dimension(GRID_KLO_GC:GRID_KHI_GC) :: x,y,z  !  NOTE -- this won't work in 1D
  integer :: istat
  real, allocatable, dimension(:) :: x,y,z
  integer :: point(MDIM)

  ! Sum quantities over all locally held leaf-node blocks.
  gsum  = 0.
  lsum = 0.
  gmax = 0.

  nse_threshold = 1. - fl_adrEps

  call Grid_getListOfBlocks(LEAF, blockList, count)
  
  do lb = 1, count
     !get the index limits of the block
     call Grid_getBlkIndexLimits(blockList(lb), blkLimits, blkLimitsGC)
     isizeGC=blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
     jsizeGC=blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
     ksizeGC=blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

     ! allocate some space for the axes
     allocate(x(isizeGC),STAT=istat)
       if (istat /= 0) print *, 'ERROR cannot allocate x in IO_writeIntegralQuantities'
     allocate(y(jsizeGC),STAT=istat)
       if (istat /= 0) print *, 'ERROR cannot allocate y in IO_writeIntegralQuantities'
     allocate(z(ksizeGC),STAT=istat)
       if (istat /= 0) print *, 'ERROR cannot allocate z in IO_writeIntegralQuantities'
     x = 0.        ! Just initialize them here.....
     y = 0.
     z = 0.

     ! get a pointer to the current block of data
     call Grid_getBlkPtr(blockList(lb), solnData)

     x = 0.

     call Grid_getCellCoords(IAXIS, blockList(lb), CENTER, .true., x, isizeGC)
     if ( NDIM > 1) call Grid_getCellCoords(JAXIS, blockList(lb), CENTER, .true., y, jsizeGC)
     if ( NDIM > 2) call Grid_getCellCoords(KAXIS, blockList(lb), CENTER, .true., z, ksizeGC)

     ! Sum contributions from the indicated blkLimits of cells.
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

              point(IAXIS) = i
              point(JAXIS) = j
              point(KAXIS) = k


              !! Get the cell volume for a single cell
              call Grid_getSingleCellVol(blockList(lb), EXTERIOR, point, dvol)

              ! mass   

#ifdef DENS_VAR
              dm = solnData(DENS_VAR,i,j,k)*dvol
              hdm = 0.5e0*dm
              lsum(1) = lsum(1) + dm
#endif           

#ifdef DENS_VAR
#ifdef VELX_VAR      
              ! momentum

              momx = dm * solnData(VELX_VAR,i,j,k)
              lsum(2) = lsum(2) + momx

#endif
#ifdef VELY_VAR      
              momy = dm * solnData(VELY_VAR,i,j,k)
              lsum(3) = lsum(3) + momy

#endif
#ifdef VELZ_VAR      
              momz = dm * solnData(VELZ_VAR,i,j,k) 
              lsum(4) = lsum(4) + momz

#endif

              ! total energy
#ifdef ENER_VAR
              lsum(5) = lsum(5) + solnData(ENER_VAR,i,j,k) * dm
#endif


#ifdef VELX_VAR      
#ifdef VELY_VAR      
#ifdef VELZ_VAR      
              ! kinetic energy
              ekin = hdm * ( &  
                   solnData(VELX_VAR,i,j,k)**2 + & 
                   solnData(VELY_VAR,i,j,k)**2 + & 
                   solnData(VELZ_VAR,i,j,k)**2 )           

              lsum(6) = lsum(6) + ekin
#endif
#endif
#endif


#ifdef EINT_VAR
              ! internal energy
              lsum(7) = lsum(7) + solnData(EINT_VAR,i,j,k)*dm
#endif
#endif ! ifdef DENS_VAR

              ! potential energy

#ifdef GPOT_VAR

              dgrav = hdm*solnData(GPOT_VAR,i,j,k)

              lsum(8) = lsum(8) + dgrav

              dedis = dm*solnData(ENER_VAR,i,j,k) + dgrav

              if ( dedis > 0.e0 ) then
                 lsum(9)  = lsum(9)  + dedis
                 lsum(10) = lsum(10) + dm
              else
                 lsum(11) = lsum(11) + dedis
                 lsum(12) = lsum(12) + dm
              end if

#endif

              ! burning quantities

#ifdef RPV1_MSCALAR 

              phi = solnData(RPV1_MSCALAR,i,j,k)

              ! burned mass, momenta and kinetic energy

              phi_dm   = phi*dm

              lsum(13) = lsum(13) + phi_dm
              lsum(14) = lsum(14) + phi*momx
              lsum(15) = lsum(15) + phi*momy
              lsum(16) = lsum(16) + phi*momz
              lsum(17) = lsum(17) + phi*ekin

              ! burning rate and area

#ifdef FSPD_VAR

              s   = solnData(FSPD_VAR,i,j,k)
              call fl_adrReactionRate(s, phi, phidot)

              lsum(18) = lsum(18) + phidot*dm                    ! mass
              if (s > 0.e0) lsum(19) = lsum(19) + dvol*phidot/s  ! area

#endif

              lsum(20) = lsum(20) + x(i)*phi_dm
              lsum(21) = lsum(21) + y(j)*phi_dm
              lsum(22) = lsum(22) + z(k)*phi_dm

#endif

              ! mass in NSQE
#ifdef RPV2_MSCALAR
              if ( solnData(RPV2_MSCALAR,i,j,k) < nse_threshold) then
                 lsum(23) = lsum(23) + solnData(RPV2_MSCALAR,i,j,k)*dm 
              end if
#endif
              ! density distribution histogram

              if ( solnData(DENS_VAR,i,j,k) > 1.e9 ) lsum(24) = lsum(24) + dm
              if ( solnData(DENS_VAR,i,j,k) > 3.e8 ) lsum(25) = lsum(25) + dm
              if ( solnData(DENS_VAR,i,j,k) > 1.e8 ) lsum(26) = lsum(26) + dm
              if ( solnData(DENS_VAR,i,j,k) > 3.e7 ) lsum(27) = lsum(27) + dm
              if ( solnData(DENS_VAR,i,j,k) > 1.e7 ) lsum(28) = lsum(28) + dm
              if ( solnData(DENS_VAR,i,j,k) > 3.e6 ) lsum(29) = lsum(29) + dm
              if ( solnData(DENS_VAR,i,j,k) > 1.e6 ) lsum(30) = lsum(30) + dm

              !! Some extra quantities added by Cal

              ! Mass of material with density > 1.5E7 which should burn to
              ! intermediate mass elements during the detonation.
              if ( solnData(DENS_VAR,i,j,k) > 1.5e7 ) lsum(34) = lsum(34) + dm

              ! Maximum density  (should be associated with the central density)
              if (solnData(DENS_VAR,i,j,k) .gt. lsum(35)) lsum(35) = solnData(DENS_VAR,i,j,k)

              ! Bubble Volume calculated by adding integrating phi*dvol
              lsum(36) = lsum(36) + phi*dvol

              !! End of Cals extra quantities
 
           enddo
        enddo
     enddo

     !! Gravity diagnostics requested by Cal/Carlo
     lsum(31) =  hy_gravMassXYZ(IAXIS)  ! acceleration*mass in I direction
     lsum(32) =  hy_gravMassXYZ(JAXIS)
     lsum(33) =  hy_gravMassXYZ(KAXIS)

     !  Cleanup and deallocate arrays

     call Grid_releaseBlkPtr(blockList(lb), solnData)

     deallocate(x,STAT=istat)
       if (istat /= 0) print *, 'ERROR cannot deallocate x in IO_writeIntegralQuantities'
     deallocate(y,STAT=istat)
       if (istat /= 0) print *, 'ERROR cannot deallocate y in IO_writeIntegralQuantities'
     deallocate(z,STAT=istat)
       if (istat /= 0) print *, 'ERROR cannot deallocate z in IO_writeIntegralQuantities'


  enddo
  
  
  ! Now the MASTER_PE sums the local contributions from all of
  ! the processors and writes the total to a file.

!  An extra call added by Cal to get the maximum density should be tested

  call MPI_Reduce (lsum(35), gmax, 1, MPI_Double_Precision, MPI_MAX, &
       &                MASTER_PE, MPI_Comm_World, error)

  
  call MPI_Reduce (lsum, gsum, nGlobalSum, MPI_Double_Precision, MPI_Sum, & 
       &                MASTER_PE, MPI_Comm_World, error)

!Ugly hack by Cal to facilitate output

  gsum(35) = gmax
  

  if (MyPE == MASTER_PE) then
     !!  Decrease gravitational acceleration terms by total mass 
     gsum(31) = gsum(31) / gsum(1) ! net acceleration in I direction
     gsum(32) = gsum(32) / gsum(1) 
     gsum(33) = gsum(33) / gsum(1)   
     
     ! create the file from scratch if it is a not a restart simulation, 
     ! otherwise append to the end of the file
     if (isfirst == 0) then
        open (funit, file=trim(statsFileName), position='APPEND')
     else 
        if (.NOT. sim_restart) then
           open (funit, file=trim(statsFileName)) 
           write (funit, 10)               &
           '#time                    ', &
           'mass                     ', &
           'x-momentum               ', &
           'y-momentum               ', & 
           'z-momentum               ', &
           'E_total                  ', &
           'E_kinetic                ', &
           'E_internal               ', & 
           'E_potential              ', &
           'E_binding                ', &
           'E_bin pos                ', &
           'mass(E_bin pos)          ', &
           'E_bin neg                ', &
           'mass(E_bin neg)          ', &
           'burnt mass               ', &
           'burnt mass x-momentum    ', &
           'burnt mass y-momentum    ', &
           'burnt mass z-momentum    ', &
           'burnt mass E_kinetic     ', &
           'burning rate             ', &
           'burning area             ', &
           'burnt mass x-centroid    ', &
           'burnt mass y-centroid    ', &
           'burnt mass z-centroid    ', &
           'mass (NSQE)              ', &
           'mass (dens > 1e9)        ', &
           'mass (dens > 3e8)        ', &
           'mass (dens > 1e8)        ', &
           'mass (dens > 3e7)        ', &
           'mass (dens > 1e7)        ', &
           'mass (dens > 3e6)        ', &
           'mass (dens > 1e6)        ', &
           'net grav accel in I      ', &
           'net grav accel in J      ', &
           'net grav accel in K      ', &
           'mass (dens > 1.5e7)      ', &
           'max dens                 ', &
           'bubble volume            '  
           
           
10         format (2x,50(a22, :, 1X))
           
        else
           open (funit, file=trim(statsFileName), position='APPEND')
           write (funit, 11) 
11         format('# simulation restarted')
        endif
     endif
     
     write (funit, 12) simtime, (gsum(i),i=1,8), gsum(5)+gsum(8), (gsum(i),i=9,19),   &
                                    (gsum(i)/gsum(13),i=20,22), (gsum(i),i=23,30), (gsum(i),i=31,nglobalSum) !,&
!                                    enuctint, (gsum(i),i=33,size(gsum))     ! Write the global sums to the file.
12   format (1x, 1P, 50(E22.15, :, 1x))
 
     close (funit)          ! Close the file.
     
  endif
  
  call MPI_Barrier (MPI_Comm_World, error)

  

  !=============================================================================
  
  return
end subroutine IO_writeIntegralQuantities



