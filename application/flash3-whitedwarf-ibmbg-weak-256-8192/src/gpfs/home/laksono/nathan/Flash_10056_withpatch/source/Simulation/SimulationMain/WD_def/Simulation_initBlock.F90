!!****if* source/Simulation/SimulationMain/WD_def/Simulation_initBlock
!!
!! NAME
!! 
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(IN) :: blockID,
!!                       integer(IN) :: myPE  )
!!
!! DESCRIPTION
!!
!!  Provide initial conditions for a Type Ia supernova problem
!!  with deflagration initialized in the center of the WD.
!!
!! ARGUMENTS
!!
!!  blockID - my block number
!!  myPE    - local processor number
!!
!!***
subroutine Simulation_initBlock(blockID,myPE)

  use Simulation_data
  use Driver_interface, ONLY : Driver_abortFlash
  use Flame_interface, ONLY : Flame_getProfile, &
    Flame_heatRelease, Flame_laminarSpeed, &
    Flame_rhJump
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_getCellCoords, Grid_getDeltas, Grid_putRowData
  use Eos_interface, ONLY : Eos
  use Flame_data, ONLY : fl_quenchingDens0

  implicit none
#include "constants.h"
#include "Flash.h"
#include "Eos.h"
  integer, intent(IN) :: blockID
  integer, intent(IN) :: myPE

  real, allocatable, dimension(:) :: xCenter, yCenter, zCenter
  real, allocatable, dimension(:) :: dx, dy, dz
  real, allocatable, dimension(:) :: rho, p, e, eint, t, game, gamc, &
       vx, vy, vz, f, trcr, ye_row, sumy, qbar_row, fspeed
  real, allocatable, dimension(:) :: xnAny
  
  logical                 :: contains_match 
  
  !front tracking stuff
  !.. runtime parameters and values calculated from them
  
  real                    :: tsim_dens_b, tpres_b, ttemp_b, teint_b
  
  
  !.. local variables
  real,allocatable, dimension(:) :: xvector, yvector, zvector
  real,allocatable, dimension(:) :: velxvector, velyvector, velzvector
  real,allocatable, dimension(:) :: densvector, presvector, &
       tempvector, enervector
  
  real                    :: dens,   pres,   temp,   ener,   phi
  real                    :: r, a, vr, ksi, r0 !,z <-- this is doubly defined
  real                    :: fronts, deltae, alpha
  
  real                    :: x_conv,  y_conv,  z_conv
  real                    :: vx_conv, vy_conv, vz_conv

  integer                 :: i, j, k, n, err, ign
  integer                 :: istat
  character(len=16)       :: ionname
  integer                 :: vecLen
  real                    :: tdens_b

  !------------------------------------------------------------------

  !..for the interpolations

  integer ::  op, kat
  real    ::  dist, xc12, xo16, dyi
  real    ::  distsq, min_distsq

  !..for calling the eos

  real :: pel, eel, sel, uint
  real :: dpt, dpd, ded, det, c_v, c_p, cs1, gamma, xalfa, xxni, xxne, xxnp
  real :: entropy, dst, dsd, abar, zbar
  real :: ye, sumyi, qbar
#if NSPECIES > 0
  real :: mfrac
#endif

  real, dimension(EOS_NUM) :: eosData
  logical,dimension(EOS_VARS+1:EOS_NUM) :: eosMask

  real, dimension(1) :: xn_dummy

  ! random number generator

  integer, save  :: i_seed
  real           :: x_seed
  integer, dimension(MDIM) :: pos
  real,dimension(MDIM) :: del
  integer,dimension(LOW:HIGH,MDIM)::blkLimits,blkLimitsGC
  integer :: iSize, jSize, kSize
  integer :: iSizeGC, jSizeGC, kSizeGC
  integer :: ilo, ihi

  !..set the number of points in the interpolation
  !..2 points = linear, 3 points = quadratic, 4 points = cubic and so on
  !..kat is a first guess at where in the table we are


  op  = 2
  kat = 1

  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  iSizeGC = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  jSizeGC = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
  kSizeGC = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

  iSize = blkLimits(HIGH,IAXIS)-blkLimits(LOW,IAXIS)+1
  jSize = blkLimits(HIGH,JAXIS)-blkLimits(LOW,JAXIS)+1
  kSize = blkLimits(HIGH,KAXIS)-blkLimits(LOW,KAXIS)+1

  ilo = blkLimits(LOW,IAXIS)
  ihi = blkLimits(HIGH,IAXIS)

  !! allocate all needed space
  allocate(xCenter(iSizeGC),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate xCenter in Simulation_initBlock")
  allocate(yCenter(jSizeGC),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate yCenter in Simulation_initBlock")
  allocate(zCenter(kSizeGC),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate zCenter in Simulation_initBlock")
  allocate(dx(iSizeGC),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate dx in Simulation_initBlock")
  allocate(dy(jSizeGC),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate dy in Simulation_initBlock")
  allocate(dz(kSizeGC),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate dz in Simulation_initBlock")
  allocate(rho(iSizeGC),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate rho in Simulation_initBlock")
  allocate(p(iSizeGC),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate p in Simulation_initBlock")
  allocate(e(iSizeGC),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate e in Simulation_initBlock")
  allocate(eint(iSizeGC),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate eint in Simulation_initBlock")
  allocate(t(iSizeGC),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate t in Simulation_initBlock")
  allocate(game(iSizeGC),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate game in Simulation_initBlock")
  allocate(gamc(iSizeGC),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate gamc in Simulation_initBlock")
  allocate(vx(iSizeGC),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate vx in Simulation_initBlock")
  allocate(vy(iSizeGC),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate vy in Simulation_initBlock")
  allocate(vz(iSizeGC),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate vz in Simulation_initBlock")
  allocate(f(iSizeGC),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate f in Simulation_initBlock")
  allocate(ye_row(iSizeGC),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate ye_row in Simulation_initBlock")
  allocate(sumy(iSizeGC),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate sumy in Simulation_initBlock")
  allocate(qbar_row(iSizeGC),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate qbar_row in Simulation_initBlock")
  allocate(trcr(iSizeGC),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate trcr in Simulation_initBlock")
#if NSPECIES > 0
  allocate( xnAny(iSizeGC),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate xnAny in Simulation_initBlock")
#endif
  allocate(xvector(iSizeGC),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate xvector in Simulation_initBlock")
  allocate(yvector(jSizeGC),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate yvector in Simulation_initBlock")
  allocate(zvector(kSizeGC),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate zvector in Simulation_initBlock")
  allocate(velxvector(iSizeGC),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate velxvector in Simulation_initBlock")
  allocate(velyvector(iSizeGC),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate velyvector in Simulation_initBlock")
  allocate(velzvector(iSizeGC),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate velzvector in Simulation_initBlock")
  allocate(densvector(iSizeGC),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate densvector in Simulation_initBlock")
  allocate(presvector(iSizeGC),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate presvector in Simulation_initBlock")
  allocate(tempvector(iSizeGC),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate tempvector in Simulation_initBlock")
  allocate(enervector(iSizeGC),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate enervector in Simulation_initBlock")
  allocate(fspeed(iSizeGC),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate fspeed in Simulation_initBlock")


  xCenter(:) = 0.e0
  yCenter(:) = 0.e0
  zCenter(:) = 0.e0

  dx(:) = 0.e0
  dy(:) = 0.e0
  dz(:) = 0.e0
  call Grid_getDeltas(blockId, del)

  if(NDIM > 2) then
     call Grid_getCellCoords(KAXIS,blockID, CENTER, .true.,zCenter,kSizeGC)
     dx=del(KAXIS)
  end if

  if(NDIM > 1) then
     call Grid_getCellCoords(JAXIS,blockID, CENTER, .true.,yCenter,jSizeGC)
     dx=del(JAXIS)
  end if

  call Grid_getCellCoords(IAXIS,blockID, CENTER, .true.,xCenter,iSizeGC)
  dx=del(IAXIS)

  t  (:) = 0.e0
  rho(:) = 0.e0
  p  (:) = 0.e0
  e  (:) = 0.e0
  vx (:) = 0.e0
  vy (:) = 0.e0
  vz (:) = 0.e0
  ye_row (:) = 0.e0
  sumy (:) = 0.e0
  qbar_row (:) = 0.e0
#if NSPECIES > 0
!  xnAny (:) = 0.e0
  mfrac = 1.0/NSPECIES
#endif
  fspeed(:) = 0.e0

 

  eosMask=.false.

  ! now fill the master arrays
  pos(IAXIS)=blkLimits(LOW,IAXIS)
  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     pos(KAXIS)=k
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        pos(JAXIS)=j
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           
           ! compute the distance from the center
           
           if ( NDIM == 1 ) then
              
              dist = abs(xCenter(i))
              
           else if ( NDIM == 2 ) then
              
              dist = sqrt(xCenter(i)**2 + yCenter(j)**2)
              
           else if ( NDIM == 3 ) then
              
              dist = sqrt(xCenter(i)**2 + yCenter(j)**2 + zCenter(k)**2)
              
           end if
           
           !..for the interior and atmosphere
           
           if ( dist <=  sim_1dRad(sim_npoints) ) then
              
              !..locate this zone in the model
              
              call ut_hunt(sim_1dRad,sim_npoints,dist,kat)
              
              kat = max(1, min(kat - op/2 + 1, sim_npoints - op + 1))
              
              ! actually we should identify the first zone in the 
              ! damped model which has density
              ! lower thatn sim_densFluff and use this to determine 
              ! numerical stellar radius; this
              ! makes us independent of possible density fluctuations
              ! far in the ambient medium
              ! or once one specifies sim_densFluff much smaller than
              ! used in damping

              !..get the interpolated values at this distance

              call ut_polint(sim_1dRad(kat),sim_1dDens(kat),op,dist,dens,dyi)
              call ut_polint(sim_1dRad(kat),sim_1dTemp(kat),op,dist,temp,dyi)
              call ut_polint(sim_1dRad(kat),sim_1dXC12(kat),op,dist,xc12,dyi)
              call ut_polint(sim_1dRad(kat),sim_1dXO16(kat),op,dist,xo16,dyi)
              
              !..If the interpolated density is less than that of the fluff,
              !..use the fluff density
              
              if ( dens < sim_densFluff ) then
                 dens    = sim_densFluff
                 temp    = sim_tempFluff
                 xc12    = sim_xc12Fluff
                 xo16    = sim_xo16Fluff
                 trcr(i) = 0.e0
              else
                 trcr(i) = 1.e0
              end if
              
           else
              
              !..for the fluff
              
              if ( sim_densFluff <= 0.e0 ) then 
                 dens   = sim_1dDens(sim_npoints)
                 temp   = sim_1dTemp(sim_npoints)
                 xc12   = sim_1dXC12(sim_npoints)
                 xo16   = sim_1dXO16(sim_npoints)
              else
                 dens = sim_densFluff
                 temp = sim_tempFluff
                 xc12 = sim_xc12Fluff
                 xo16 = sim_xo16Fluff
              end if
              
              trcr(i) = 0.e0
              
           end if
           
           ! set initial composition to unburned
           
           phi = 0.e0
           eosData(EOS_TEMP) = temp
           eosData(EOS_DENS) = dens
           eosData(EOS_ABAR) = sim_eosDataU(EOS_ABAR)
           eosData(EOS_ZBAR) = sim_eosDataU(EOS_ZBAR)
           qbar=sim_qbarU

           vecLen=1

           call Eos(MODE_DENS_TEMP, vecLen,eosData,xn_dummy,eosMask)
           abar = eosData(EOS_ABAR)
           zbar = eosData(EOS_ZBAR)
           pres = eosData(EOS_PRES)
           uint = eosData(EOS_EINT)
           gamma = eosData(EOS_GAMC)
!!$           call eos(dens,temp,pres,uint,sim_mfrac,entropy,  &
!!$                abar,zbar,dpt,dpd,det,ded,dst,dsd,  c_v,c_p,gamma,pel, &
!!$                xxne,xalfa,1)

           !-- begin if ignite -------------------------------------------
           !
           
           !if igniting the model set a spherical region to the burned state
           !note that this is really a perturbation on the above. If ignite
           !is set to .FALSE., then one should wind up with just the WD model.
           
           if ( sim_ignite ) then
              
              if ( sim_ignitionFile ) then
                 ! multi point ignition,  points were read from file
                 ! find the minimum distance to any ignition point
                 min_distsq = HUGE(1.0)
                 do ign = 1, sim_num_ign
                    if (NDIM == 1) then
                       distsq = (xCenter(i) - sim_ign_x(ign))**2
                    else if (NDIM == 2) then
                       distsq = (xCenter(i) - sim_ign_x(ign))**2 + &
                                (yCenter(j) - sim_ign_y(ign))**2
                    else if (NDIM == 3) then
                       distsq = (xCenter(i) - sim_ign_x(ign))**2 + &
                                (yCenter(j) - sim_ign_y(ign))**2 + &
                                (zCenter(k) - sim_ign_z(ign))**2
                    endif
                    if (distsq < min_distsq) then
                       min_distsq = distsq
                       sim_rMatch = sim_ign_r(ign)
                    endif
                 enddo
                 dist = sqrt(min_distsq)

              else 
                 ! single point ignition according to location and radius in flash.par
                 if ( sim_rMatch <= 0.e0 ) then
                    
                    contains_match = (xCenter(i)-0.5e0*dx(i) <= sim_xMatch).and.&
                                     (xCenter(i)+0.5e0*dx(i) >= sim_xMatch)
                    
                    if ( NDIM >= 2 ) then
                       contains_match = contains_match .and. &
                                        (yCenter(j)-0.5e0*dy(j) <= sim_yMatch)&
                                        .and.&
                                        (yCenter(j)+0.5e0*dy(j) >= sim_yMatch)
                    end if

                    if ( NDIM == 3 ) then
                       contains_match = contains_match .and. &
                                        (zCenter(k)-0.5e0*dz(k) <= sim_zMatch)&
                                        .and.   &
                                        (zCenter(k)+0.5e0*dz(k) >= sim_zMatch)
                    end if
                 
                    if ( contains_match ) then
                       call Driver_abortFlash('[SIMULATION_INITBLOCK] ERROR: not implemented')
                    endif

                 else

                    dist = (xCenter(i)-sim_xMatch)**2
                 
                    if (NDIM > 1) then
                       dist = dist + (yCenter(j)-sim_yMatch)**2
                    end if

                    if(NDIM > 2) then
                       dist = dist + (zCenter(k)-sim_zMatch)**2
                    end if

                    dist = sqrt(dist)
                 end if
              end if
              
              r0 = dist
              
                 
              ! set the value of phi
              if ( r0 < sim_rMatch - sim_drMin ) then
                 phi = 1.e0
              else if ( r0 > sim_rMatch + sim_drMin ) then
                 phi = 0.e0
              else
                 call Flame_getProfile( r0-sim_rMatch, phi )
              end if
              
              ! If radius is less than that of matchhead + width of the flame
              ! use the flame stuff. Otherwise, leave things alone.
              
              if ( r0 <= sim_rMatch + sim_drMin)then
                 
                 ! set the new abundances
                 
                 
                 ye = (1.e0-phi)* sim_eosDataU(EOS_ZBAR)/sim_eosDataU(EOS_ABAR) +&
                      phi* sim_eosDataNse(EOS_ZBAR)/sim_eosDataNse(EOS_ABAR)
                 sumyi = (1.e0-phi)/sim_eosDataU(EOS_ABAR) + phi/sim_eosDataNse(EOS_ABAR) 
                 qbar = (1.e0-phi)*sim_qbarU + phi*sim_qbarNse
                 
                 deltae = sim_deltaeNse
                 
                 call Flame_laminarSpeed(dens, temp, fronts)
                 
                 ! For now, and this is really a test, I am faking the
                 ! partial burning by treating phi as a factor of the
                 ! heat release.  As I see it now, the fact that the
                 ! constant burned and unburned states are globals,
                 ! sim_uMfrac and sim_bMfrac, keeps me from having
                 ! less than a complete burn in the jump
                 ! calculation. But, scaling the energy release by the
                 ! factor phi should produce the equivalent amount of
                 ! energy for the state. The jump calculation should
                 ! not adjust the abundances, so those will stay set
                 ! as they should be.
                 
                 call Flame_rhJump(dens, pres, temp, uint, zbar/abar, &
                      1.e0/abar, tdens_b, tpres_b, ttemp_b, teint_b, ye, &
                      sumyi, phi*deltae, 0.0, MODE_DENS_TEMP)
                 
                 temp = ttemp_b
                 dens = tdens_b
                 uint = teint_b
                 
                 ! This gives the 'temporary' burned state given the
                 ! input state
                 ! To keep hydrostatic equlibrium, something we have
                 ! been striving so
                 ! dilligently for, use the orignial pressure and
                 ! burned density to get
                 ! the temperature and energy
                 
                 eosData(EOS_TEMP)=temp
                 eosData(EOS_DENS)=dens
                 eosData(EOS_EINT)=teint_b
                 eosData(EOS_ABAR)=1.e0/sumyi
                 eosData(EOS_ZBAR)=ye*eosData(EOS_ABAR)
                 vecLen=1
                 
                 call Eos(MODE_DENS_TEMP, vecLen,eosData,xn_dummy,eosMask)
                 dens = eosData(EOS_DENS)
                 temp = eosData(EOS_TEMP)
                 pres = eosData(EOS_PRES)
                 abar = eosData(EOS_ABAR)
                 zbar = eosData(EOS_ZBAR)
                 uint = eosData(EOS_EINT)
                 gamma = eosData(EOS_GAMC)
                 
!!$                 call eos(dens,temp,pres,uint,sim_mfrac,entropy,  &
!!$                    abar,zbar,dpt,dpd,det,ded,dst,dsd,  c_v,c_p,gamma,pel, &
!!$                    xxne,xalfa,3)

                 if( sim_detonate ) then  
                    ! reset phi1 = 0.0 to start detonation instead of flame
                    phi = 0.0
                 endif
                 
              end if
           end if
           !-- end if ignite -----------------------------------------



           !-- initial velocity field
           
           if ( sim_vPert > 0.e0 ) then
              
              ! random fluctuations
              
              call random_number (harvest=x_seed)
              vx(i) = sim_vPert * 2.e0 * (x_seed - 0.5e0)
              vx(i) = vx(i)*sim_sdimI
              
              if ( NDIM == 2 ) then
                 call random_number (harvest=x_seed)
                 vy(i) = sim_vPert * 2.e0 * (x_seed - 0.5e0)
                 vy(i) = vy(i)*sim_sdimI
              end if
              
              if ( NDIM == 3 ) then
                 call random_number (harvest=x_seed)
                 vz(i) = sim_vPert * 2.e0 * (x_seed - 0.5e0)
                 vz(i) = vz(i)*sim_sdimI
              end if
              
           else
              
              ! assume the star is at rest
              
              vx(i) = 0.e0
              vy(i) = 0.e0
              vz(i) = 0.e0
              
           end if
              
           ! add initial convection, if needed
           
           if( sim_vConv > 0.e0 .AND. sim_rConv > 0.e0  ) then
              if( NDIM == 2 .AND. xCenter(i)**2+yCenter(j)**2 <= sim_rc2)then
                 vx(i) = vx(i) + sim_vConv*xCenter(i)*yCenter(j)/sim_rc2
                 vy(i) = vy(i) + sim_vConv*(1.e0-(2.e0*xCenter(i)**2+&
                         yCenter(j)**2)/sim_rc2)
              else if( NDIM == 3 .AND. xCenter(i)**2+yCenter(j)**2+&
                   zCenter(k)**2 <= sim_rc2 ) then
                    ! Rotate into the frame of reference of the dipole
                    ! phi:   rotation around z axis
                    ! theta: rotation around y axis (left-handed)
                 x_conv =  sim_cosPhi*sim_cosTheta*xCenter(i) + &
                           sim_sinPhi*sim_cosTheta*yCenter(j) + &
                           sim_sinTheta*zCenter(k)
                 y_conv = -  sim_sinPhi*xCenter(i) +&
                             sim_cosPhi*yCenter(j)
                 z_conv = -sim_cosPhi*sim_sinTheta*xCenter(i) - &
                           sim_sinPhi*sim_sinTheta*yCenter(j) + &
                           sim_cosTheta*zCenter(k)

                 ! Compute velocities in this frame of reference
                 vx_conv = sim_vConv*x_conv*z_conv/sim_rc2
                 vy_conv = sim_vConv*y_conv*z_conv/sim_rc2
                 vz_conv = sim_vConv*(1.e0-(2.e0*(x_conv**2+y_conv**2)+&
                           z_conv**2)/sim_rc2)

                 ! Convert velocities to the lab frame
                 vx(i) = vx(i) + sim_cosPhi*sim_cosTheta*vx_conv - &
                         sim_sinPhi*vy_conv - &
                         sim_cosPhi*sim_sinTheta*vz_conv
                 vy(i) = vy(i) + sim_sinPhi*sim_cosTheta*vx_conv + &
                         sim_cosPhi*vy_conv - &
                         sim_sinPhi*sim_sinTheta*vz_conv
                 vz(i) = vz(i) + sim_sinTheta*vx_conv + sim_cosTheta*vz_conv
              end if
           end if
              
           ! fill in the flash arrays
           
           t   (i) = temp
           rho (i) = dens
           p   (i) = pres
           e   (i) = uint + 0.5e0*(vx(i)**2+vy(i)**2+vz(i)**2)
           eint(i) = uint
           gamc(i) = gamma
           game(i) = pres/(dens*uint) + 1.e0
           ye_row(i) = zbar/abar
           sumy(i) = 1.e0/abar
           qbar_row(i) = qbar

           f(i) = phi

#if NSPECIES > 0
           xnAny(i) = mfrac
#endif

           ! this is a hack to get the refinement conditions to act correctly
           ! the real flamespeed cannot be calculated until the
           ! gravitational potential is initialized, which hasn't happened yet
           ! these values will be overwritten before they are used for anything real
           if (rho(i) > fl_quenchingDens0) then
              fspeed(i) = 1.e0
           else
              fspeed(i) = 0.e0
           endif
           
           ! end of 3d loops
        end do
        
        call Grid_putRowData(blockID, CENTER, DENS_VAR, EXTERIOR,IAXIS,&
             pos,rho(ilo:ihi),isize)
        call Grid_putRowData(blockID, CENTER, PRES_VAR, EXTERIOR,IAXIS,&
             pos,p(ilo:ihi),iSize)
        call Grid_putRowData(blockID, CENTER, TEMP_VAR, EXTERIOR,IAXIS,&
             pos,t(ilo:ihi),iSize)
        call Grid_putRowData(blockID, CENTER, VELX_VAR, EXTERIOR,IAXIS,&
             pos,vx(ilo:ihi),iSize)
        call Grid_putRowData(blockID, CENTER, VELY_VAR, EXTERIOR,IAXIS,&
             pos,vy(ilo:ihi),iSize)
        call Grid_putRowData(blockID, CENTER, VELZ_VAR, EXTERIOR,IAXIS,&
             pos,vz(ilo:ihi),iSize)
        call Grid_putRowData(blockID, CENTER, GAMC_VAR, EXTERIOR,IAXIS,&
             pos,gamc(ilo:ihi),iSize)
        call Grid_putRowData(blockID, CENTER, GAME_VAR, EXTERIOR,IAXIS,&
             pos,game(ilo:ihi),iSize)
        call Grid_putRowData(blockID, CENTER, ENER_VAR, EXTERIOR,IAXIS,&
             pos,e(ilo:ihi),iSize)
#ifdef EINT_VAR
        call Grid_putRowData(blockID, CENTER, EINT_VAR, EXTERIOR,IAXIS,&
             pos,eint(ilo:ihi),iSize)
#endif
        call Grid_putRowData(blockID, CENTER, FLAM_MSCALAR, EXTERIOR,IAXIS,&
             pos,f(ilo:ihi),iSize)
        call Grid_putRowData(blockID, CENTER, RPV1_MSCALAR, EXTERIOR,IAXIS,&
             pos,f(ilo:ihi),iSize)
        call Grid_putRowData(blockID, CENTER, RPV2_MSCALAR, EXTERIOR,IAXIS,&
             pos,f(ilo:ihi),iSize)
        call Grid_putRowData(blockID, CENTER, YE_MSCALAR, EXTERIOR,IAXIS,&
             pos,ye_row(ilo:ihi),iSize)
        call Grid_putRowData(blockID, CENTER, SUMY_MSCALAR, EXTERIOR,IAXIS,&
             pos,sumy(ilo:ihi),iSize)
        call Grid_putRowData(blockID, CENTER, QBAR_MSCALAR, EXTERIOR,IAXIS,&
             pos,qbar_row(ilo:ihi),iSize)
        
#ifdef TRCR_VAR
        call Grid_putRowData(blockID, CENTER, TRCR_VAR, EXTERIOR,IAXIS,&
             pos,trcr(ilo:ihi),iSize)
#endif
        do n = SPECIES_BEGIN,SPECIES_END
           call Grid_putRowData(blockID, CENTER,n, EXTERIOR,IAXIS,&
                pos,xnAny(ilo:ihi),iSize)
        end do

        call Grid_putRowData(blockID, CENTER, FSPD_VAR, EXTERIOR,IAXIS,&
             pos,fspeed(ilo:ihi),iSize)
        
     end do
  end do
  
  deallocate(xCenter,STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot deallocate xCenter in Simulation_initBlock")
  deallocate(yCenter,STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot deallocate yCenter in Simulation_initBlock")
  deallocate(zCenter,STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot deallocate zCenter in Simulation_initBlock")
  deallocate(dx,STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot deallocate dx in Simulation_initBlock")
  deallocate(dy,STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot deallocate dy in Simulation_initBlock")
  deallocate(dz,STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot deallocate dz in Simulation_initBlock")
  deallocate(rho,STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot deallocate rho in Simulation_initBlock")
  deallocate(p,STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot deallocate p in Simulation_initBlock")
  deallocate(e,STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot deallocate e in Simulation_initBlock")
  deallocate(eint,STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot deallocate eint in Simulation_initBlock")
  deallocate(t,STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot deallocate t in Simulation_initBlock")
  deallocate(game,STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot deallocate game in Simulation_initBlock")
  deallocate(gamc,STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot deallocate gamc in Simulation_initBlock")
  deallocate(vx,STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot deallocate vx in Simulation_initBlock")
  deallocate(vy,STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot deallocate vy in Simulation_initBlock")
  deallocate(vz,STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot deallocate vz in Simulation_initBlock")
  deallocate(f,STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot deallocate f in Simulation_initBlock")
  deallocate(ye_row,STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot deallocate ye_row in Simulation_initBlock")
  deallocate(sumy,STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot deallocate sumy in Simulation_initBlock")
  deallocate(qbar_row,STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot deallocate qbar_row in Simulation_initBlock")
  deallocate(trcr,STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot deallocate trcr in Simulation_initBlock")
#if NSPECIES > 0
  deallocate( xnAny,STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot deallocate xnAny in Simulation_initBlock")
#endif
  deallocate(xvector,STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot deallocate xvector in Simulation_initBlock")
  deallocate(yvector,STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot deallocate yvector in Simulation_initBlock")
  deallocate(zvector,STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot deallocate zvector in Simulation_initBlock")
  deallocate(velxvector,STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot deallocate velxvector in Simulation_initBlock")
  deallocate(velyvector,STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot deallocate velyvector in Simulation_initBlock")
  deallocate(velzvector,STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot deallocate velzvector in Simulation_initBlock")
  deallocate(densvector,STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot deallocate densvector in Simulation_initBlock")
  deallocate(presvector,STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot deallocate presvector in Simulation_initBlock")
  deallocate(tempvector,STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot deallocate tempvector in Simulation_initBlock")
  deallocate(enervector,STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot deallocate enervector in Simulation_initBlock")
  deallocate(fspeed,STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot deallocate fspeed in Simulation_initBlock")


  return
end subroutine Simulation_initBlock

