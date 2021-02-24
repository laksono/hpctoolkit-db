!!****if* source/Simulation/SimulationMain/WD_def/Grid_markRefineDerefine
!!
!! NAME
!!  Grid_markRefineDerefine
!!
!! SYNOPSIS
!!
!!  Grid_markRefineDerefine(integer(IN) :: myPE)
!!  
!! DESCRIPTION 
!!  Mark blocks for refinement or derefinement
!!  This routine has been customized for the WD detonation problem.
!!
!! ARGUMENTS
!!  myPE : my processor id.
!! 
!! NOTES
!!
!!  WARNING this version of the routine refines the propagating flame
!!          higher than anything else.  This is not yet well tested!!!!!
!!  
!!  Several kinds of conditions are enforced:
!!
!!  1. gradient checks at high density
!!     checks on rpv1, dens, |velecity|
!!     these are weighted with a log density scale so that they are
!!     stronger at higher density (enables us to de-emphasise the stars edge)
!!  2. forced refinement/derefinement
!!     enforce an elevated minimum refinement for r < refine_uni_rad or dens > uni_dens
!!     further enforce an elevated minimum refinement of inner_dens_dx for rho > inner_dens_min.
!!       Used to spread out particles over more processors.
!!     enforce a reduced maximum refinement for r > refine_max_rad
!!       for GCD problem refine_max_rad can be made larger near south pole
!!  3. enforce full refinement where flame is actively propagating
!!  4. GCD collision region override
!!     make resolution close to gcd_focus_dx near south pole
!!
!!  The above description 2007-02-22, there are several older checks
!!  related to GCD and derfinement with time that are commented out or otherwise disabled
!!
!!  Dean Townsley 2007/02/22
!!  Spring cleaning: Casey Meakin, 21 June 2007
!!
!!***


subroutine Grid_markRefineDerefine(MyPE)

  use Grid_data, ONLY : gr_refine_cutoff, gr_derefine_cutoff,&
       gr_refine_filter, &
       gr_numRefineVars,gr_refine_var,gr_refineOnParticleCount
  use Driver_interface, ONLY : Driver_getSimTime
  use tree, ONLY : newchild, refine, derefine, stay,lrefine, lrefine_max,lrefine_min, &
                   nodetype, parent, child, nchild
  use Simulation_data, ONLY : sim_refineDDens, sim_refineDtVel, &
       sim_refineXtVel, sim_refineDphi1, sim_refineXphi1, &
       sim_refineSphi1, sim_refineXenuc, sim_refineDensMin, &
       sim_refineDensMax, sim_refineUniDens, sim_refineUniDx, &
       sim_refineUniRadius, sim_refineMaxRadius, &
       sim_refineIgnitionTime, sim_refineIgnitionRadius, &
       sim_lrefineNonFlameDecrement, sim_smallU, &
       sim_xMatch, sim_yMatch, sim_zMatch, sim_gridGeom, &
       sim_lrefineMax, sim_lrefineMin, sim_refine_inner_dens_min, &
       sim_refine_inner_dens_dx, sim_GCDRefineAngle, sim_GCDRefineMaxRadius, &
       sim_GCDFocusTime, sim_GCDFocusAngle, sim_GCDFocusMinRadius, &
       sim_GCDFocusMaxRadius, sim_GCDFocusDx, sim_refine_allphi
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_getListOfBlocks, &
    Grid_releaseBlkPtr, Grid_getBlkIndexLimits, Grid_getCellCoords, &
    Grid_getDeltas, Grid_fillGuardCells

  implicit none

  include 'mpif.h'
#include "constants.h"
#include "Flash.h"

  integer,intent(IN) :: MyPE

  real :: ref_cut,deref_cut,ref_filter
  integer       :: l,iref

  integer, dimension(2,MDIM) :: blkLimits
  integer, dimension(2,MDIM) :: blkLimitsGC
  integer ,dimension(MAXBLOCKS) :: blkList
  integer :: blkCount

  logical :: addLevel
  real    :: d_, dj, dj_, dref, maxphi1, mdens, mpres, errh, errl, dx, nxbh, flameerrh, flameerrl
  real    :: ldwg_min, ldwg_max, ldwg_
  integer :: i, j, k, iblk, ii, jj, kk
  integer :: isizeGC, jsizeGC, ksizeGC
  integer :: isize, jsize, ksize
  real    :: time

  real    :: xll, yll, zll, rll, height, costheta
  real    :: xli, yli, zli, rli

  real    :: refineMaxRadius_, tderef
  real    :: vphi1l, vphi2l, vphi3l, vphi1h, vphi2h, vphi3h
  real    :: v1_min, v1_max, v2_min, v2_max, v3_min, v3_max
  
  integer, dimension(MAXBLOCKS) :: lrefineDel_

  real, dimension(:,:,:,:), pointer :: solnData
#ifdef FSPD_SCRATCH_GRID_VAR
  real, dimension(:,:,:,:), pointer :: scratchData
#endif

  logical :: gcMask(NUNK_VARS) ! indicates for which variables the refinement criteria
  ! that are implemented in this file need values in guardcells

  logical, dimension(MAXBLOCKS) :: refine_par
  integer :: nsend, nrecv, ierr
  integer, dimension(MAXBLOCKS) :: reqr, reqs
  integer, dimension(MPI_STATUS_SIZE,MAXBLOCKS) :: statr, stats

#ifdef FIXEDBLOCKSIZE
  real,    dimension(GRID_ILO_GC:GRID_IHI_GC, GRID_JLO_GC:GRID_JHI_GC, GRID_KLO_GC:GRID_KHI_GC) ::  &
       phi1, phi2, phi3, tvel, dwg, varmin, varmax
  real,    dimension(GRID_ILO_GC:GRID_IHI_GC) :: xLeft
  real,    dimension(GRID_JLO_GC:GRID_JHI_GC) :: yLeft 
  real,    dimension(GRID_KLO_GC:GRID_KHI_GC) :: zLeft
#else
  integer :: istat
  real, allocatable, dimension(:,:,:) :: phi1, tvel, dwg, varmin, varmax, phi2, phi3
  real, allocatable, dimension(:) :: xLeft, yLeft, zLeft
#endif

  newchild(:) = .FALSE.
  refine(:)   = .FALSE.
  derefine(:) = .FALSE.
  stay(:)     = .FALSE.

!  the code below would allow addition of standard refinement conditions vi flash.par.
!  However, the logic at the end here does not respect the values of refine(:),
!    derefine(:) and stay(:) set by them, so they won't work anyway.
!  We could rewrite logic to respect them so this can be enabled.
!  We would also have to call Grid_fillGuardCells before the next lines
!  to ensure that gr_markRefineDerefine sees valid data in guard cells.
!
!  do l = 1,gr_numRefineVars
!     iref = gr_refine_var(l)
!     ref_cut = gr_refine_cutoff(l)
!     deref_cut = gr_derefine_cutoff(l)
!     ref_filter = gr_refine_filter(l)
!     call gr_markRefineDerefine(MyPE,iref,ref_cut,deref_cut,ref_filter)
!  end do

  
  !---------------------------------------------------------------!
  ! START Applying problem-specific refinement criteria.
  !---------------------------------------------------------------!
  


  !---------------------------!
  ! GET SIMULATION TIME FOR
  ! TIME DEPENDENT REFINEMENT
  ! CRITERIA
  !---------------------------!
  call Driver_getSimTime(time)
  
  !---------------------!
  ! GUARD CELL FILLING
  !---------------------!
  gcMask = .FALSE.
#ifdef DENS_VAR
  gcMask(DENS_VAR) = .TRUE.
#endif
#ifdef RPV1_MSCALAR
  gcMask(RPV1_MSCALAR) = .TRUE.
#endif
#ifdef FSPD_VAR
  gcMask(FSPD_VAR) = .TRUE.
#endif
#ifdef ENUC_VAR
  !!  gcMask(ENUC_VAR) = .TRUE.
#endif
#ifdef VELX_VAR
  gcMask(VELX_VAR) = .TRUE.
#endif
#ifdef VELY_VAR
#if NDIM > 1
  gcMask(VELY_VAR) = .TRUE.
#endif
#endif
#ifdef VELZ_VAR
#if NDIM > 2
  gcMask(VELZ_VAR) = .TRUE.
#endif
#if defined(FLASH_MHD_DIM)
  if (NDIM .GE. 2) gcMask(VELZ_VAR) = .TRUE.
#endif
#if defined(FLASH_2P5DIM)
  if (NDIM .GE. 2) gcMask(VELZ_VAR) = .TRUE.
#endif
#endif


  call Grid_fillGuardCells(MyPE,CENTER,ALLDIR,doEos=.true.,&
       maskSize=NUNK_VARS, mask=gcMask, makeMaskConsistent=.true.)

    !--------------------------!
  ! THE APPROPRIATE GUARD
  ! CELLS ARE NOW FILLED
  !--------------------------!
  



  !--------------------------!
  ! GET LIST OF BLOCKS AND
  ! BEGIN LOOP OVER THEM
  !--------------------------!  
  call Grid_getListOfBlocks(ACTIVE_BLKS, blkList,blkCount)

  do iblk = 1, blkCount
     
     if (lrefine(blkList(iblk)) < lrefine_min) then
        refine(blkList(iblk)) = .TRUE.
        derefine(blkList(iblk)) = .FALSE.
        CYCLE                   ! Don't even try to analyze contents of this block - KW
     end if


     call Grid_getBlkPtr(blkList(iblk),solnData,CENTER)
#ifdef FSPD_SCRATCH_GRID_VAR
     call Grid_getBlkPtr(blkList(iblk),scratchData,SCRATCH)
#endif
     call Grid_getBlkIndexLimits(blkList(iblk), blkLimits, blkLimitsGC)     

     isizeGC=blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
     jsizeGC=blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
     ksizeGC=blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1
     
     isize=blkLimits(HIGH,IAXIS)-blkLimits(LOW,IAXIS)+1
     jsize=blkLimits(HIGH,JAXIS)-blkLimits(LOW,JAXIS)+1
     ksize=blkLimits(HIGH,KAXIS)-blkLimits(LOW,KAXIS)+1
     
     nxbh = 0.5e0*isize
     
#ifndef FIXEDBLOCKSIZE
     allocate(phi1(isizeGC,jsizeGC,ksizeGC),STAT=istat)
     if (istat /= 0) print *,' ERROR in allocating phi1 in Grid_markRefineDerefine'
     allocate(phi2(isizeGC,jsizeGC,ksizeGC),STAT=istat)
     if (istat /= 0) print *,' ERROR in allocating phi2 in Grid_markRefineDerefine'
     allocate(phi3(isizeGC,jsizeGC,ksizeGC),STAT=istat)
     if (istat /= 0) print *,' ERROR in allocating phi3 in Grid_markRefineDerefine'
     allocate(tvel(isizeGC,jsizeGC,ksizeGC),STAT=istat)
     if (istat /= 0) print *,' ERROR in allocating tvel in Grid_markRefineDerefine'
     allocate(dwg(isizeGC,jsizeGC,ksizeGC),STAT=istat)
     if (istat /= 0) print *,' ERROR in allocating dwg in Grid_markRefineDerefine'
     allocate(varmin(isizeGC,jsizeGC,ksizeGC)),STAT=istat)
     if (istat /= 0) print *,' ERROR in allocating varmin in Grid_markRefineDerefine'
     allocate(varmax(isizeGC,jsizeGC,ksizeGC),STAT=istat)
     if (istat /= 0) print *,' ERROR in allocating varmax in Grid_markRefineDerefine'
     allocate(xLeft(isizeGC),STAT=istat)
     if (istat /= 0) print *,' ERROR in allocating xLeft in Grid_markRefineDerefine'
     allocate(yLeft(jsizeGC),STAT=istat)
     if (istat /= 0) print *,' ERROR in allocating yLeft in Grid_markRefineDerefine'
     allocate(zLeft(ksizeGC),STAT=istat)
     if (istat /= 0) print *,' ERROR in allocating zLeft in Grid_markRefineDerefine'
#endif
     


     !------------------------------------------!
     ! INITIALIZE VARIABLES TO BE USED TO KEEP 
     ! TRACK OF REFINMENT STATUS
     !------------------------------------------!
     errh = 0.e0
     errl = 0.e0
     

     !------------------------------------!
     ! STORE MAXIMUM DENSITY FOR BLOCK
     ! (not including guardcells)
     !------------------------------------!
     mdens = 0.e0     
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              mdens = max(mdens, solnData(DENS_VAR,i,j,k))
           end do
        end do
     end do


     !--------------------------------------------------!
     ! READ IN FIRST PROGRESS VARIABLES IF NEEDED LATER
     !--------------------------------------------------!
     !     if( sim_refineIgnitionTime > time .or. &
     !          (min(sim_refineDensMin,sim_refineDensMax) > 0.e0 .AND.  &
     !          mdens > sim_refineDensMin)) then
     do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
        do j = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
           do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)
              phi1(i,j,k)  = solnData(RPV1_MSCALAR,i,j,k)
              phi2(i,j,k) = solnData(RPV2_MSCALAR,i,j,k)
              phi3(i,j,k) = solnData(RPV3_MSCALAR,i,j,k)
           end do
        end do
     end do
     ! (OBSOLETE: endif)
     
     
     !-------------------------------------------------------------------------------!
     ! PART 1: REFINEMENT CONDITION BASED ON VARIABLE GRADIENTS (LOG DENSITY WEIGHTED).
     ! (these conditions are only checked for blocks with large enough density)
     !-------------------------------------------------------------------------------!
     if ( min(sim_refineDensMin,sim_refineDensMax) > 0.e0 .AND.  &
          mdens > sim_refineDensMin ) then
        
        ! for weighting indicators with the first volume scalar
        ! and density
        ! density weights for refinement
        
        ldwg_min = log(sim_refineDensMin)
        ldwg_max = log(sim_refineDensMax)
        ldwg_    = 1.e0/(ldwg_max-ldwg_min)
        
        maxphi1 = 0.e0
        ! store weighting factor, phi1 and get max phi1 in block&gc
        do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
           do j = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
              do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)                 
                 if (phi1(i,j,k) > maxphi1) maxphi1 = phi1(i,j,k)
                 dwg(i,j,k) = max(0.e0,(log(solnData(DENS_VAR,i,j,k))-ldwg_min)*ldwg_)
              end do
           end do
        end do
        

        ! DENSITY GRADIENT CONDITION
        iref = DENS_VAR
        call varminmax(solnData(iref,:,:,:), varmin, varmax, blkLimitsGC)
        dj = 0.e0
        do k = blkLimitsGC(LOW,KAXIS)+1*K3D, blkLimitsGC(HIGH,KAXIS)-1*K3D
           do j = blkLimitsGC(LOW,JAXIS)+1*K2D, blkLimitsGC(HIGH,JAXIS)-1*K2D
              do i = blkLimitsGC(LOW,IAXIS)+1, blkLimitsGC(HIGH,IAXIS)-1
                 dj_ = varmax(i,j,k)/varmin(i,j,k) - 1.e0
                 dj  = max(dj, dj_*dwg(i,j,k))
              end do
           end do
        end do
        dref = sim_refineDDens
        errh = errh + max(0.e0, dj-dref)
        errl = errl + max(0.e0, dj-0.33e0*dref)
        

        
        ! TOTAL VELOCITY GRADIENT CONDITION
        do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
           do j = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
              do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)
                 dj_ =       solnData(VELX_VAR,i,j,k)**2
#if NDIM >= 2
                 dj_ = dj_ + solnData(VELY_VAR,i,j,k)**2
#endif
#if NDIM == 3
                 dj_ = dj_ + solnData(VELZ_VAR,i,j,k)**2
#endif
                 tvel(i,j,k) = sqrt(dj_+sim_smallU)
              end do
           end do
        end do
        call varminmax(tvel, varmin, varmax, blkLimitsGC)
        dj = 0.e0
        do k = blkLimitsGC(LOW,KAXIS)+1*K3D, blkLimitsGC(HIGH,KAXIS)-1*K3D
           do j = blkLimitsGC(LOW,JAXIS)+1*K2D, blkLimitsGC(HIGH,JAXIS)-1*K2D
              do i = blkLimitsGC(LOW,IAXIS)+1, blkLimitsGC(HIGH,IAXIS)-1
                 if ( varmin(i,j,k) > sim_refineXtvel ) then
                    dj_ = varmax(i,j,k)/varmin(i,j,k) - 1.e0
                    dj  = max(dj, dj_*dwg(i,j,k)*phi1(i,j,k))
                 end if
              end do
           end do
        end do
        dref = sim_refineDtvel
        errh = errh + max(0.e0, dj-dref)
        errl = errl + max(0.e0, dj-0.33e0*dref)


        ! FIRST PROGRESS VARIABLE (WHERE VALUE EXCEDES THRESHOLD)
        if ( maxphi1 > sim_refineXphi1 ) then           
           call varminmax(phi1, varmin, varmax, blkLimitsGC)
           dj = 0.e0
           do k = blkLimitsGC(LOW,KAXIS)+1*K3D, blkLimitsGC(HIGH,KAXIS)-1*K3D
              do j = blkLimitsGC(LOW,JAXIS)+1*K2D, blkLimitsGC(HIGH,JAXIS)-1*K2D
                 do i = blkLimitsGC(LOW,IAXIS)+1, blkLimitsGC(HIGH,IAXIS)-1
                    dj_ = varmax(i,j,k)/(varmin(i,j,k)+sim_refineSphi1) - 1.e0
                    dj  = max(dj, dj_)
                 end do
              end do
           end do
           dref = sim_refineDphi1
           errh = errh + max(0.e0, dj-dref)
           errl = errl + max(0.e0, dj-0.33e0*dref)
        end if
        
     end if
     !--------------------------------------!
     ! END DENSE REGION CONDITIONAL AND
     ! GRADIENT-BASED REFINMENT CONDITIONS
     !--------------------------------------!






     ! -----------------------------------------------------------------------!
     ! PART 2: BEGIN UNIFORM AND FORCED REFINMENT AND DEREFINMENT CONDITIONS
     !------------------------------------------------------------------------!

     !--------------------------------------------------!
     ! CALCULATE GEOMETRIC DISTANCES AND ANGLES:
     ! rll is dist from center to closest block corner
     ! rli is dist from ignition point to closest cell 
     !   corner in this block
     !--------------------------------------------------!
     call Grid_getCellCoords(IAXIS, blkList(iblk), LEFT_EDGE, .true., xLeft, isizeGC)
     
     dx = xLeft(blkLimits(LOW,IAXIS)) - xLeft(blkLimits(LOW,IAXIS)-1)
     xll = min(abs(xLeft(blkLimits(LOW,IAXIS))),abs(xLeft(blkLimits(HIGH,IAXIS)+1)))
     xli = abs(xLeft(blkLimits(LOW,IAXIS)) - sim_xMatch)
     do i = blkLimits(LOW,IAXIS)+1, blkLimits(HIGH,IAXIS)+1
        xli = min(xli, abs(xLeft(i) - sim_xMatch))
     end do
     
     if ( NDIM > 1) then 
        call Grid_getCellCoords(JAXIS, blkList(iblk), LEFT_EDGE, .true., yLeft, jsizeGC)
        yll = min(abs(yLeft(blkLimits(LOW,JAXIS))),abs(yLeft(blkLimits(HIGH,JAXIS)+1)))
        yli = abs(yLeft(blkLimits(LOW,JAXIS)) - sim_yMatch)
        do i = blkLimits(LOW,JAXIS)+1, blkLimits(HIGH,JAXIS)+1
           yli = min(yli, abs(yLeft(i) - sim_yMatch))
        end do
     else
        yli = 0.0
     endif
     
     if ( NDIM > 2 ) then
        call Grid_getCellCoords(KAXIS, blkList(iblk), LEFT_EDGE, .true., zLeft, ksizeGC)
        zll = min(abs(zLeft(blkLimits(LOW,KAXIS))),abs(zLeft(blkLimits(HIGH,KAXIS)+1)))
        zli = abs(zLeft(blkLimits(LOW,KAXIS)) - sim_zMatch)
        do i = blkLimits(LOW,KAXIS)+1, blkLimits(HIGH,KAXIS)+1
           zli = min(zli, abs(zLeft(i) - sim_zMatch))
        end do
     else
        zli = 0.e0
     endif
     
     select case(sim_gridGeom)
     case (CARTESIAN)
        rll = sqrt( xll**2 + yll**2 + zll**2 )
        rli = sqrt( xli**2 + yli**2 + zli**2 )

        ! assume bubble rising along z-direction
        if ( abs(zLeft(blkLimits(LOW,KAXIS))) < abs(zLeft(blkLimits(HIGH,KAXIS)+1)) ) then
           height = zLeft(blkLimits(LOW,KAXIS))
        else
           height = zLeft(blkLimits(HIGH,KAXIS)+1)
        endif
     case (CYLINDRICAL)
        rll = sqrt( xll**2 + yll**2 )
        rli = sqrt( xli**2 + yli**2 )

        ! bubble rising along y-direction
        if ( abs(yLeft(blkLimits(LOW,JAXIS))) < abs(yLeft(blkLimits(HIGH,JAXIS)+1)) ) then
           height = yLeft(blkLimits(LOW,JAXIS))
        else
           height = yLeft(blkLimits(HIGH,JAXIS)+1)
        endif
     case (SPHERICAL)
        rll = xll
        rli = xli
        ! no bubble
        height = 0.e0
     end select
     
     ! COS (POLAR ANGLE)
     if (rll == 0.0) then
        costheta = 0.0
     else
        costheta = height/rll
     end if



     !----------------------------------------------!
     ! ENFORCE REGION OF ENHANCED UNIFORM REFINMENT
     !----------------------------------------------!
     if ( sim_refineUniDx > 0.e0 ) then
        addLevel = .false.
        ! resolve dense regions
        if ( sim_refineUniDens > 0.e0 .AND. mdens > sim_refineUniDens ) &
             addLevel = .true.
        ! ensure uniform resolution up to certain distance from (0,0,0)
        if ( sim_refineUniRadius > 0.e0 ) then
           if ( rll < sim_refineUniRadius ) addLevel = .true.
        end if
        ! if the additional refinement criteria are met...
        if ( addLevel ) then
           ! effect  0.5*uni_dx < dx < uni_dx
           ! force refinement if the resolution is coarse
           if ( dx > 1.01e0*sim_refineUniDx ) errh = 1.e0
           ! prevent derefinement if we are just at the right resolution
           if ( dx > 0.5e0*sim_refineUniDx ) errl = 1.e0
        end if
     end if

     
     
     !-----------------------------------------!
     ! ADD COMMENT 
     !-----------------------------------------!
     ! further enhanced refinement inside WD 
     ! to better balance particles    
     !-----------------------------------------!
     if (mdens .gt. sim_refine_inner_dens_min) then
        ! effect  0.5*inner_dens_dx < dx < inner_dens_dx
        ! force refinement if the resolution is coarse
        if ( dx > 1.01e0*sim_refine_inner_dens_dx ) errh = 1.e0
        ! prevent derefinement if we are just at the right resolution
        if ( dx > 0.5e0*sim_refine_inner_dens_dx ) errl = 1.e0
     end if
     
     
     !------------------------------------------!
     ! ADD COMMENT 
     !------------------------------------------!
     ! forced derefinement at large radii
     !   i.e. depress max refine at large radii
     !------------------------------------------!
     if ( sim_refineMaxRadius > 0.e0 ) then
        !  for GCD studies
        !  push out derefine boundary near south pole
        !  (note units are taken care of in initialization of GCDRefineAngle)
        if ( sim_GCDRefineAngle > 0.e0 .and. costheta <= -cos(sim_GCDRefineAngle) ) then
           refineMaxRadius_ = sim_GCDRefineMaxRadius
        else
           refineMaxRadius_ = sim_refineMaxRadius
        end if
        
        if ( rll > refineMaxRadius_ ) then
           ! force resolution to be > 2 * refine_uni_dx
           if ( dx < 2.e0 * sim_refineUniDx ) then
              ! force derefine
              errh = 0.e0
              errl = 0.e0
           else if (dx < 4.001e0 * sim_refineUniDx ) then
              ! prevent refine
              errh = 0.e0
           end if           
        end if
     end if


     !---------------------------------------------------------------!
     ! ADD COMMENT 
     !---------------------------------------------------------------!
     ! ensure full refinment in areas where the flame is propagating
     !
     ! this is necessary because the flame is "tuned" to have the 
     ! right width and propagation speed on this grid resolution
     ! this comes after the force derefine above so that it can 
     ! override if necessary
     !---------------------------------------------------------------!
     
     flameerrh=0.e0
     flameerrl=0.e0
     if ( sim_refineIgnitionTime > time ) then
        ! ensure full refinement "inside" flame front
        v1_min = 1.e0
        v1_max = 0.e0

        v2_min = 1.e0
        v2_max = 0.e0

        v3_min = 1.e0
        v3_max = 0.e0

        do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
           do j = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
              do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)
                 
                 if( sim_refine_allphi) then
                    v1_min = min(v1_min, phi1(i,j,k))
                    v1_max = max(v1_max, phi1(i,j,k))
                    
                    v2_min = min(v2_min, phi2(i,j,k))
                    v2_max = max(v2_max, phi2(i,j,k))
                    
                    v3_min = min(v3_min, phi3(i,j,k))
                    v3_max = max(v3_max, phi3(i,j,k))
                    
                 else                 
#ifdef FSPD_SCRATCH_GRID_VAR
                    if (scratchData(FSPD_VAR,i,j,k) > 0.0) then
#else
                    if (solnData(FSPD_VAR,i,j,k) > 0.0) then
#endif
                       v1_min = min(v1_min, phi1(i,j,k) )
                       v1_max = max(v1_max, phi1(i,j,k) )
                    endif
                    
                 endif

                    
              end do
           end do
        end do
        
        ! want to trigger refinement when v_max > 0.1 AND v_min < 0.9
        !  that is            when (started burning)  AND (not done burning)
        ! but we need to have some hysteresis
        
        vphi1h = min( max(0.e0, v1_max - 0.1e0),max(0.e0, (1.0-v1_min) - 0.1e0) )
        vphi1l = min( max(0.e0, v1_max - 0.33e0*0.1e0),max(0.e0, (1.0-v1_min) - 0.33e0*0.1e0) )
        
        if( sim_refine_allphi) then 
           vphi2h = min( max(0.e0, v2_max - 0.1e0),max(0.e0, (1.0-v2_min) - 0.1e0) )
           vphi3h = min( max(0.e0, v3_max - 0.1e0),max(0.e0, (1.0-v3_min) - 0.1e0) )
           vphi2l = min( max(0.e0, v2_max - 0.33e0*0.1e0),max(0.e0, (1.0-v2_min) - 0.33e0*0.1e0) )
           vphi3l = min( max(0.e0, v3_max - 0.33e0*0.1e0),max(0.e0, (1.0-v3_min) - 0.33e0*0.1e0) )
           errh = errh + max(max(vphi1h,vphi2h),vphi3h)
           errl = errl + max(max(vphi1l,vphi2l),vphi3l)
        else
           errh = errh + vphi1h
           errl = errl + vphi1l
        endif
        
        vphi1l = min( max(0.e0, v1_max - 0.66e0*0.1e0),max(0.e0, (1.0-v1_min) - 0.66e0*0.1e0) )
        
        if( sim_refine_allphi) then 
           vphi2l = min( max(0.e0, v2_max - 0.66e0*0.1e0),max(0.e0, (1.0-v2_min) - 0.66e0*0.1e0) )
           vphi3l = min( max(0.e0, v3_max - 0.66e0*0.1e0),max(0.e0, (1.0-v3_min) - 0.66e0*0.1e0) )
           flameerrh = flameerrh + max(max(vphi1h,vphi2h),vphi3h)
           flameerrl = flameerrl + max(max(vphi1l,vphi2l),vphi3l)
        else
           flameerrh = flameerrh + vphi1h
           flameerrl = flameerrl + vphi1l
        endif
        
     end if



     !--------------------------------------------------------!
     ! ADD COMMENT 
     !--------------------------------------------------------!
     ! force specified resolution in a predefined region
     !  near the south pole after a certain time
     !  this is an override and therefore must be at the end
     !--------------------------------------------------------!

     if ( time > sim_GCDFocusTime ) then
        if ( sim_GCDFocusAngle > 0.0 ) then
           if ( (rll > sim_GCDFocusMinRadius) .and. (rll < sim_GCDFocusMaxRadius)  &
                .and. (costheta < -cos(sim_GCDFocusAngle)) ) then
              ! we want the refinement level which is closest to gcd_focus_dx
              ! this is done by checking increments
              if ( abs(0.5*dx-sim_GCDFocusDx) < abs(dx-sim_GCDFocusDx) ) then
                 ! refine
                 errh = 1.e0
              else if ( abs(2*dx)-sim_GCDFocusDx < abs(dx-sim_GCDFocusDx) ) then
                 ! derefine
                 errh = 0.e0
                 errl = 0.e0
              else
                 ! stay
                 errh = 0.e0
                 errl = 1.e0
              endif
           endif
        endif
     endif
     



     !---------------------------!
     ! COLLECT WISDOM: MARK FOR
     ! REFINMENT/DEREFINMENT
     !---------------------------!
     ! check all jumps
     ! refine or derefine as allowed if it is allowed
     ! only flameerr tests can push us over lrefinemax-nonflamedecrement
     !defaults...
     refine(blkList(iblk)) = .FALSE.
     derefine(blkList(iblk)) = .FALSE.
     !now test
     if (     ( lrefine(blkList(iblk)) < sim_lrefineMin )  &
         .or. ( errh > 0.e0      .and. lrefine(blkList(iblk)) < (sim_lrefineMax-sim_lrefineNonflameDecrement) ) &
         .or. ( flameerrh > 0.e0 .and. lrefine(blkList(iblk)) < sim_lrefineMax ) ) then
        refine(blkList(iblk)) = .TRUE.
     endif
     if (     ( lrefine(blkList(iblk)) > sim_lrefineMax ) &
         .or. ( lrefine(blkList(iblk)) > (sim_lrefineMax-sim_lrefineNonflameDecrement) .and. flameerrl <= 0.e0 )  &
         .or. ( errl <= 0.e0 .and. lrefine(blkList(iblk)) > sim_lrefineMin) ) then
        derefine(blkList(iblk)) = .TRUE.
     endif


     !---------------------------------!
     ! CLEANUP:
     ! REALEASE CURRENT BLOCK POINTER 
     ! AND DEALLOCATE TEMP DATA.
     !---------------------------------!
     call Grid_releaseBlkPtr(blkList(iblk),solnData)
#ifdef FSPD_SCRATCH_GRID_VAR
     call Grid_releaseBlkPtr(blkList(iblk),scratchData, SCRATCH)
#endif
     
     
#ifndef FIXEDBLOCKSIZE
     deallocate(phi1,STAT=istat)
     if (istat /= 0) print *,' ERROR deallocating phi1 in Grid_markRefineDerefine'
     deallocate(phi2,STAT=istat)
     if (istat /= 0) print *,' ERROR deallocating phi2 in Grid_markRefineDerefine'
     deallocate(phi3,STAT=istat)
     if (istat /= 0) print *,' ERROR deallocating phi3 in Grid_markRefineDerefine'
     deallocate(tvel,STAT=istat)
     if (istat /= 0) print *,' ERROR deallocating tvel in Grid_markRefineDerefine'
     deallocate(dwg,STAT=istat)
     if (istat /= 0) print *,' ERROR deallocating dwg in Grid_markRefineDerefine'
     deallocate(varmin,STAT=istat)
     if (istat /= 0) print *,' ERROR deallocating varmin in Grid_markRefineDerefine'
     deallocate(varmax,STAT=ista)
     if (istat /= 0) print *,' ERROR deallocating varmax in Grid_markRefineDerefine'
     deallocate(xLeft,STAT=istat)
     if (istat /= 0) print *,' ERROR deallocating xLeft in Grid_markRefineDerefine'
     deallocate(yLeft,STAT=istat)
     if (istat /= 0) print *,' ERROR deallocating yLeft in Grid_markRefineDerefine'
     deallocate(zLeft,STAT=istat)
     if (istat /= 0) print *,' ERROR deallocating zLeft in Grid_markRefineDerefine'
#endif
     
  end do !----END LOOP OVER BLOCKS----!
  
  



  !------------------------------------------------------------------------------


  !---------------------------------------------------------------!
  ! CONSITENCY CHECK: UPDATE REFINMENT CONDITION BASED ON PARENTS:
  ! for children that are marked derefine, check if parent is 
  ! marked refine and if so unmark derefine
  !
  !---------------------------------------------------------------!
  call Grid_getListOfBlocks(ALL_BLKS, blkList,blkCount)
  refine_par(:) = .false.
  ! children collect messages from parents
  !    message id is child block number on local processor
  nrecv = 0
  do iblk = 1, blkCount
     i = blkList(iblk)
     if (parent(1,i) > 0) then
        if (parent(2,i)/=myPE) then
           nrecv = nrecv+1
           call MPI_IRecv(refine_par(i), 1, MPI_LOGICAL, parent(2,i), &
                                    i, MPI_COMM_WORLD, reqr(nrecv), ierr)
        else
           refine_par(i) = refine(parent(1,i))
        endif
     endif
  end do
  ! parents send error
  nsend = 0
  do iblk = 1, blkCount
     i = blkList(iblk)
     do j = 1,nchild
        if (child(1,j,i) > 0) then
           if (child(2,j,i) /= myPE) then
              nsend = nsend + 1
              call MPI_ISend(refine(i), 1, MPI_LOGICAL, child(2,j,i), &
                                  child(1,j,i), MPI_COMM_WORLD, reqs(nsend), ierr)
           endif
        endif
     enddo
  enddo
  if (nsend > 0) then
     call MPI_Waitall(nsend,reqs,stats,ierr)
  endif
  if (nrecv > 0) then
     call MPI_Waitall(nrecv,reqr,statr,ierr)
  endif
  do iblk = 1, blkCount
     i = blkList(iblk)
     if (nodetype(i) == LEAF .and. derefine(i) .and. refine_par(i) ) then
        derefine(i) = .false.
     endif
  enddo

  !---------!
  ! DONE 
  !---------!

  if(gr_refineOnParticleCount)call gr_ptMarkRefineDerefine()


  return

end subroutine Grid_markRefineDerefine
