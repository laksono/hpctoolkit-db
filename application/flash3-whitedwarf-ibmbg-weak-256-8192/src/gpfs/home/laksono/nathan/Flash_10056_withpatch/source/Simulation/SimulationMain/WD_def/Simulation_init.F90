!!****if* source/Simulation/SimulationMain/WD_def/Simulation_init
!!
!! NAME
!! 
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  Simulation_init(integer(IN) :: myPE)
!!
!! DESCRIPTION
!!
!!  Provide initial conditions for a Type Ia supernova problem
!!  with deflagration initialized in the center of the WD.
!!
!! ARGUMENTS
!!
!!    myPE - my processor number
!!
!! PARAMETERS
!!
!!   dens_fluff = density of the fluff
!!   temp_fluff = temperature of fluff
!!   xc12_fluff = c12 mass fraction of fluff
!!   xo16_fluff = o16 mass fraction of fluff
!!   ignite     = whether to create ignition region
!!   x_match    = x-coordinate of the ignition region
!!   y_match    = y-coordinate of the ignition region
!!   z_match    = z-coordinate of the ignition region
!!   r_match    = radius of the ignition region
!!   r_conv     = radius of initial convection region
!!   v_conv     = velocity in initial convection region
!!
!!***

subroutine Simulation_init(myPE)

  use Simulation_data
  use bn_paraData, ONLY : pbEbC12,pbAC12,pbEbO16,pbAO16

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, ONLY : Logfile_stampMessage
  use Grid_interface, ONLY : Grid_getGeometry
  use Flame_interface, ONLY : Flame_heatRelease, &
                                   Flame_rhJump,&
                                   Flame_laminarSpeed,&
                                   Flame_getWidth, &
                                   Flame_rhJumpReactive
  use Burn_data, ONLY : bn_useBurn 
  use Driver_interface, ONLY : Driver_abortFlash
  use Driver_interface, ONLY : Driver_abortFlash


  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"

  integer, intent(IN) :: myPE

  real :: v0, p0, p1, p2
  real :: dens_u,pres_u, temp_u, eint_u, uint, pres, laminarWidth
  real, save :: pi
  real :: alpha
  integer :: n_pert
  integer    :: jdens, jtemp, jxc12, jxo16
  integer, parameter ::  NVARS_MODEL = 24
  real,dimension( NVARS_MODEL)  :: tv
  character (len=4)       :: unklabels(NVARS_MODEL)
  integer :: i, j, i_seed, istat
  character (len=40) :: infile, string, ignitionFileName
  integer :: inlen
  real :: fronts, deltae
  integer :: op, kat
  real :: dyi, deltaeNse
  real :: dist

  uint = 0.e0
  pres = 0.e0

  if (.not. bn_useBurn) then 
     call Driver_abortFlash("WD_def now requires useBurn=true for normal operation") 
  endif

  call RuntimeParameters_get('smallu',     sim_smallU)

  call RuntimeParameters_get('dens_fluff', sim_densFluff)

  call RuntimeParameters_get('temp_fluff', sim_tempFluff)
  
  call RuntimeParameters_get('xc12_fluff', sim_xc12Fluff)
  
  call RuntimeParameters_get('xo16_fluff', sim_xo16Fluff)
  
  call RuntimeParameters_get('rstar', sim_rstar)
  
  call RuntimeParameters_get('x_match', sim_xMatch)
  
  call RuntimeParameters_get('y_match', sim_yMatch)
  
  call RuntimeParameters_get('z_match', sim_zMatch)
  
  call RuntimeParameters_get('r_match', sim_rMatch)
  
  call RuntimeParameters_get('ignite', sim_ignite)
  call RuntimeParameters_get('detonate', sim_detonate)

  call RuntimeParameters_get('ignition_file', sim_ignitionFile)
  call RuntimeParameters_get('ignition_file_name', ignitionFileName)
  
  call RuntimeParameters_get('v_pert', sim_vPert)
  
  call RuntimeParameters_get('n_pert', n_pert)
  call RuntimeParameters_get('r_conv', sim_rConv)
  call RuntimeParameters_get('v_conv', sim_vConv)

  call RuntimeParameters_get('refine_ddens',           sim_refineDDens)
  call RuntimeParameters_get('refine_dtvel',           sim_refineDtVel)
  call RuntimeParameters_get('refine_xtvel',           sim_refineXtVel)
  call RuntimeParameters_get('refine_dphi1',           sim_refineDphi1)
  call RuntimeParameters_get('refine_xphi1',           sim_refineXphi1)
  call RuntimeParameters_get('refine_sphi1',           sim_refineSphi1)
  call RuntimeParameters_get('refine_xenuc',           sim_refineXenuc)

  call RuntimeParameters_get('refine_allphi',           sim_refine_allphi)

  call RuntimeParameters_get('refine_dens_min',        sim_refineDensMin)
  call RuntimeParameters_get('refine_dens_max',        sim_refineDensMax)
  call RuntimeParameters_get('refine_uni_dens',        sim_refineUniDens)
  call RuntimeParameters_get('refine_uni_dx',          sim_refineUniDx)
  call RuntimeParameters_get('refine_uni_radius',      sim_refineUniRadius)
  call RuntimeParameters_get('refine_max_radius',      sim_refineMaxRadius)

  call RuntimeParameters_get('refine_inner_dens_min',  sim_refine_inner_dens_min)
  call RuntimeParameters_get('refine_inner_dens_dx',   sim_refine_inner_dens_dx)

  call RuntimeParameters_get('gcd_refine_max_radius',  sim_GCDRefineMaxRadius)
  call RuntimeParameters_get('gcd_refine_angle',       sim_GCDRefineAngle)
  sim_GCDRefineAngle = PI/180.0*sim_GCDRefineAngle

  call RuntimeParameters_get('gcd_focus_time',        sim_GCDFocusTime)
  call RuntimeParameters_get('gcd_focus_max_radius',  sim_GCDFocusMaxRadius)
  call RuntimeParameters_get('gcd_focus_min_radius',  sim_GCDFocusMinRadius)
  call RuntimeParameters_get('gcd_focus_angle',       sim_GCDFocusAngle)
  call RuntimeParameters_get('gcd_focus_dx',          sim_GCDFocusDx)
  sim_GCDFocusAngle = PI/180.0*sim_GCDFocusAngle

  call RuntimeParameters_get('refine_ignition_time',   sim_refineIgnitionTime)
  call RuntimeParameters_get('refine_ignition_radius', sim_refineIgnitionRadius)

  call RuntimeParameters_get('lrefine_min',   sim_lrefineMin)
  call RuntimeParameters_get('lrefine_max',   sim_lrefineMax)
  call RuntimeParameters_get('lrefine_nonflamedecrement',   sim_lrefineNonflameDecrement)
  if (sim_lrefineNonflameDecrement > (sim_lrefineMax-sim_lrefineMin) ) then
     call Driver_abortFlash("refinement max-min must be greater than nonflame decrement")
  endif

  call RuntimeParameters_get('restart',   sim_restart)
  
  call Grid_getGeometry(sim_gridGeom)

  if( NDIM == 3 ) then
     call RuntimeParameters_get('phi_conv', sim_phiConv)
     call RuntimeParameters_get('theta_conv', sim_thetaConv)
     pi = PI
     sim_phiConv   = sim_phiConv  *  pi/180.
     sim_thetaConv = sim_thetaConv * pi/180.
     
     sim_cosPhi   = cos(sim_phiConv)
     sim_sinPhi   = sin(sim_phiConv)
     sim_cosTheta = cos(sim_thetaConv)
     sim_sinTheta = sin(sim_thetaConv)
  end if

  sim_rc2 = sim_rConv*sim_rConv
  
  if ( sim_vPert > 0.e0 ) then
     
     if ( n_pert < 0 ) then
        n_pert = 0
        write(*,*)  '[SIMULATION_INITBLOCK] Warning: n_pert reset to 0'
        call Logfile_stampMessage(myPE, '[SIMULATION_INITBLOCK] Warning: n_pert reset to 0')
     end if
     
     do i = 1, n_pert*MyPE + MyPE + 1
        call random_seed (size=i_seed)
     end do
     
     sim_sdimI = 1.e0/sqrt(float(NDIM))
     
  end if

  !-----------------------------
  ! READ LIST OF IGNITION POINTS
  !-----------------------------
  if (sim_ignitionFile) then
     call Logfile_stampMessage(myPE, '[SIMULATION_INITBLOCK] Reading ignition points from file')
     open(unit=2,file=ignitionFileName,status='OLD',iostat=istat)
     if (istat /= 0) call Driver_abortFlash("Unable to open ignition points file")

     ! eat header
     read(2,*)
     read(2,*) sim_num_ign
     do i = 1, sim_num_ign
        read(2,*) sim_ign_x(i), sim_ign_y(i), sim_ign_z(i), sim_ign_r(i)
     enddo
!     if (MyPE == MASTER_PE) then
!        write(*,*)  'ignition points x y z r'
!        do i = 1, sim_num_ign
!           write(*,*) sim_ign_x(i), sim_ign_y(i), sim_ign_z(i), sim_ign_r(i)
!        enddo
!     endif
  endif
  
  !-------------------------------------
  !  READ INITIAL PROFILE FROM DATA FILE
  !-------------------------------------
  !..open the data file
  
  ! for original unrelaxed cold wd model
  ! infile = 'coldwd_mchandra.dat'
  ! nvar_ = 4
  
  ! for damped cold wd model
  infile = 'coldwd_mchandra_damped.dat'
  !nvar_ = 24

  inlen  = index(infile,' ') - 1
  
  open(unit=2,file=infile,status='old')
  
  !..strip the header
  
  read(2,'(a)') string
  read(2,'(a)') string
  
  ! identify variables
  
  do i = 1,NVARS_MODEL
     read(2,'(a)') unklabels(i)
     if ( unklabels(i) == 'dens' ) jdens = i
     if ( unklabels(i) == 'temp' ) jtemp = i
     if ( unklabels(i) == 'C12 ' ) jxc12 = i
     if ( unklabels(i) == 'O16 ' ) jxo16 = i
  end do

  if ( MyPE == MASTER_PE ) then
     write(*,*) 'density     : ',jdens
     write(*,*) 'temperature : ',jtemp
     write(*,*) 'XC12        : ',jxc12
     write(*,*) 'XO16        : ',jxo16
  end if

     !..read the radius, density, temperature, xc12 and xo16

  sim_npoints = 0
  
  do i = 1,SIM_NMAX
     read(2,*,end=11) sim_1dRad(i), (tv(j), j=1,NVARS_MODEL)
  
     sim_1dDens(i) = tv(jdens)
     sim_1dTemp(i) = tv(jtemp)
     sim_1dXC12(i) = tv(jxc12)
     sim_1dXO16(i) = tv(jxo16)
     
!!$     read(2,*,end=11) &
!!$          sim_1dRad(i),tv(1), tv(2), tv(3), tv(4), tv(5), sim_1dDens(i),         &
!!$          tv(7), tv(8), tv(9), tv(10), sim_1dTemp(i), tv(12), tv(13), tv(14), &
!!$          tv(15), sim_1dXC12(i), sim_1dXO16(i)
     
     sim_npoints = i
     
     ! preserving temperature profile resulting from the relaxation
     ! seems producing more stable models
     
     ! sim_1dTemp(i) = temp_fluff
  end do
  
11 close(unit=2)
  
  if (MyPE == MASTER_PE) then
     write(*,'(1x,a,a,a,i4)')' done reading file = ',infile(1:inlen), & 
          ' number of data points =',sim_npoints
  end if
  
  ! front tracking stuff
  !.. find unburned compositions: composition at the center of the cold star
  

  if (NDIM == 1) then
     dist = (sim_xMatch)**2
  else if (NDIM == 2) then
     dist = (sim_xMatch)**2 +                 &
          (sim_yMatch)**2
  else
     dist = (sim_xMatch)**2 +                 &
          (sim_yMatch)**2 +                 &
          (sim_zMatch)**2
  end if
  dist = sqrt(dist)

  ! find index in 1d model
  ! assumes that ignition point is not some insane spot outside WD
  
  !..set the number of points in the interpolation
  !..2 points = linear, 3 points = quadratic, 4 points = cubic and so on
  !..kat is a first guess at where in the table we are
  op  = 2
  kat = 1

  call ut_hunt(sim_1dRad, sim_nPoints, dist, kat)

  kat = max(1, min(kat - op/2 + 1, sim_npoints - op + 1))
  !..get the interpolated values at this distance

  call ut_polint(sim_1dRad(kat),sim_1dDens(kat),op,dist,sim_eosDataU(EOS_DENS),dyi)
  call ut_polint(sim_1dRad(kat),sim_1dTemp(kat),op,dist,sim_eosDataU(EOS_TEMP),dyi)
  ! composition

  call ut_polint(sim_1dRad(kat),sim_1dXC12(kat),op,dist,sim_cFrac,dyi)
  
  sim_eosDataU(EOS_ABAR) = 1.e0 / (sim_cFrac/pbAC12 + (1.e0-sim_cFrac)/pbAO16)
  sim_eosDataU(EOS_ZBAR) = 0.5e0 * sim_eosDataU(EOS_ABAR)
  sim_qbarU = sim_cFrac*pbEbC12/pbAC12 + (1.e0-sim_cFrac)*pbEbO16/pbAO16

  !.. find burned composition: nse at pressure of unburned material

  call Flame_laminarSpeed(sim_eosDataU(EOS_DENS), sim_eosDataU(EOS_TEMP), fronts)

  call Flame_rhJumpReactive(sim_eosDataU, sim_qbarU, sim_eosDataNse, sim_qbarNse, MODE_DENS_TEMP) 
  sim_deltaeNse = (sim_qbarNse - sim_qbarU) * 9.6485e17
  sim_eosDataB = sim_eosDataNse

  !.. some constants to initialize pressure and velocity
  
  alpha = sim_eosDataU(EOS_DENS)/sim_eosDataB(EOS_DENS)
  
  sim_v0 = fronts * (alpha - 1.e0)
  sim_p0 = sim_eosDataU(EOS_PRES)
  sim_p1 = 2.e0 * alpha * (alpha - 1.e0) * sim_eosDataU(EOS_DENS) * fronts*fronts
  sim_p2 = sim_p1 * (alpha - 1.e0) / (4.e0*alpha)
  
  !.. width of the "narrow band"
  call Flame_getWidth(laminarWidth)
  sim_drMin = 1.5e0*laminarWidth
  




  
end subroutine Simulation_init


