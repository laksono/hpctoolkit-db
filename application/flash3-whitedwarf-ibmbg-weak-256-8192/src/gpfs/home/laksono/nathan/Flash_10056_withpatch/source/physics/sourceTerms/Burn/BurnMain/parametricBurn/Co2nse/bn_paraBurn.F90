!!****if* source/physics/sourceTerms/Burn/BurnMain/parametricBurn/Co2nse/bn_paraBurn
!!
!! NAME
!!
!!  bn_paraBurn
!!
!! SYNOPSIS
!!
!!  bn_paraBurn(  real, intent(in)     :: dens,
!!             real, intent(in)     :: temp,
!!             real, intent(in)     :: eint,
!!             real, intent(in)     :: xc12init,
!!             real, intent(inout)  :: phi1,
!!             real, intent(inout)  :: phi2,
!!             real, intent(inout)  :: phi3,
!!             real, intent(inout)  :: flame,
!!             real, intent(inout)  :: ye,
!!             real, intent(inout)  :: sumyi,
!!             real, intent(inout)  :: qbar,
!!             real, intent(out)    :: qdot,
!!             real, intent(in)     :: dt)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   dens : density 
!!
!!   temp : temperature
!!
!!   eint : internal energy
!!
!!   xc12init : fraction of initial carbon abundance
!!
!!   phi1 : progress variable 1
!!
!!   phi2 : progress variable 2
!!
!!   phi3 : progress variable 3
!!
!!   ye : electron fraction
!!
!!   sumyi : 
!!
!!   qbar : 
!!
!!   qdot : 
!!
!!   dt : timestep over which to burn material
!!
!!***

subroutine bn_paraBurn(dens, temp, eint, pres, xc12init, phi1, phi2, phi3, flame, &
     ye, sumyi, qbar, qdot, dt)
  
  use bn_paraData, ONLY : &
       pbEbC12, pbAC12, pbEbO16, pbAO16, pbEbMg24, pbAMg24,&
       pbFqdotMult, pbNseWarnLevel, &
       pbClight, pbNA, pbMp, pbMn, pbElectronMass, pbReact

  use Burn_interface, ONLY:  Burn_nseAtDens, Burn_nseAtPres  

  implicit none
  
#include "Flash.h"
#include "constants.h"
  
  real, intent(in)    :: dens, temp, eint, pres, xc12init
  real, intent(inout) :: phi1, phi2, phi3
  real, intent(in)    :: flame
  real, intent(inout) :: ye, sumyi, qbar
  real, intent(out)   :: qdot
  real, intent(in)    :: dt

  
  !----------!
  !LOCAL VARS!
  !----------!
  real, parameter :: ash_update_threshold = 1.e-6
  integer :: i, j
  real    :: dti
  
  !INITIAL STATE OF PLASMA
  real :: phi1_i, phi2_i, phi3_i
  real :: phi1dot, phi2dot, phi3dot
  real :: eint_i, pres_i
  real :: qbar_i, qbar_ash, qbar_oldash, qbar_oldash_i
  real :: sumyi_ash, sumyi_oldash, sumyi_oldash_i
  real :: ye_i, ye_ash,  ye_oldash, ye_oldash_i, sumyi_i
  real :: xc12_i, xo16_i, xmg24_i, xcarbon

  !RELAXATION TIMESCALE PARAMETERS
  real :: tau_nsqe, tau_nse, maxExp, testExp  
  
  !NSE CALCULATION - TABLE DATA
  real :: s
  real :: qbar_finnse_d, sumyi_finnse_d, temp_finnse_d, edot_d, yedot_d
  real :: qbar_finnse_p, sumyi_finnse_p, temp_finnse_p, edot_p, yedot_p
  real :: qbar_finnse, sumyi_finnse, temp_finnse, edot, yedot

  !CARBON BURNING REACTION RATE PARAMETERS
  real :: nasigmav, t9, t9a
  integer :: icburn

  !-------------------------------------------------------------------------------------
  !inverse timestep
  dti    = 1.e0/dt
  
  !initial plasma properties
  qbar_i  = qbar
  ye_i    = ye  
  sumyi_i = sumyi
  eint_i  = eint
  pres_i  = pres
  phi1_i  = phi1
  phi2_i  = phi2
  phi3_i  = phi3
  
  
  xc12_i  = (1.e0-phi1)*xc12init
  xo16_i  = (1.e0-phi2)*(1.e0-xc12init)
  xmg24_i = (phi1-phi2)*xc12init


  ! Properties of ash from pevious timestep:
  ! Used at end of subroutine to find new state after burn.
  if (phi2 > ash_update_threshold) then !Divide by zero test.
     ye_oldash_i   = 0.5e0 + (ye-0.5e0)/phi2
     
     qbar_oldash_i = (qbar - xc12_i/pbAC12*pbEbC12 &
          - xo16_i/pbAO16*pbEbO16 &
          - xmg24_i/pbAMg24*pbEbMg24)/phi2
     
     sumyi_oldash_i = (sumyi - xc12_i/pbAC12 &
          - xo16_i/pbAO16 &
          - xmg24_i/pbAMg24)/phi2
  else
     ye_oldash_i    = ye
     sumyi_oldash_i = sumyi
     qbar_oldash_i  = qbar
  endif
  
  
  
  ! CONSTRAIN YE>=0.42, WARNING MESSAGE IF TRIGGERED
  if (ye_i <0.42) then
     if (pbNseWarnLevel >=1) then
        write(6,*) 'pbburnNse: ye_i bad on entry, ye_i= ', ye_i, ', resetting to 0.42'
     endif
     ye_i = 0.42
  endif

  

  !--------------------------------------------------------------------------
  ! NSE DATA: qbar, sumyi, temp, edot, yedot for ISOCHORIC and ISOBARIC BURN
  !--------------------------------------------------------------------------
  ! Retrieve qbar, sumyi, temp, edot, and, yedot in NSE state from tables.
  ! Tables available for isochoric and isobaric burn. Choice of isochoric,
  ! isobaric, and transition (linear interp. between) deteremined by second
  ! progress variable (phi2) when we are in a flame front (deflagration).
  !--------------------------------------------------------------------------
  if ( flame >= 0.9999 .or. flame < ash_update_threshold ) then !ISOCHORIC BURN DATA
     call Burn_nseAtDens(qbar_finnse, sumyi_finnse, temp_finnse, edot, yedot, &
       ye_i, dens, eint_i - qbar_i*9.6485e17)
  
  else if (flame >= 0.99) then !LINEAR INTERPOLATION BETWEEN ISOCHORIC AND ISOBARIC DATA
     call Burn_nseAtDens(qbar_finnse_d, sumyi_finnse_d, temp_finnse_d, edot_d, yedot_d, &
            ye_i, dens, eint_i - qbar_i*9.6485e17)
     call Burn_nseAtPres(qbar_finnse_p, sumyi_finnse_p, temp_finnse_p, edot_p, yedot_p, &
          ye_i, pres_i, eint_i + pres_i/dens - qbar_i*9.6485e17)
     s = (phi2 - 0.99)/(0.9999-0.99)
     qbar_finnse = s*qbar_finnse_d + (1.0-s)*qbar_finnse_p
     sumyi_finnse = 1.0/(s/sumyi_finnse_d + (1.0-s)/sumyi_finnse_p)
     temp_finnse = s*temp_finnse_d + (1.0-s)*temp_finnse_p
     edot = s*edot_d + (1.0-s)*edot_p
     yedot = s*yedot_d + (1.0-s)*yedot_p
       
  else !ISOBARIC BURN DATA
     call Burn_nseAtPres(qbar_finnse, sumyi_finnse, temp_finnse, edot, yedot, &
          ye_i, pres_i, eint_i + pres_i/dens - qbar_i*9.6485e17)
  endif
  
  
  
  !-----------------------------------------------------------------------
  !                   MAKE EXPLICIT ESTIMATES FOR BURNING RATES 
  !                     AND TIME ADVANCE PROGRESS VARIABLES
  !-----------------------------------------------------------------------
  
  !--------------------------
  ! STAGE 1: CARBON BURNING
  !--------------------------
  
  phi1dot = 0.e0

  !---- UPDATE PHI1 -----
  if( (.not. pbReact) .or. ( flame .ge. ash_update_threshold ) ) then
     ! follow flame tracking variable
     phi1 = flame
     phi1dot = (phi1 - phi1_i)*dti
  else
     ! burn based on nuclear reaction rate
     icburn = 1
     t9 = temp/1.d9
  
     !USE REACTION RATE: icburn = 1
     if(icburn.eq.1.and.t9.gt.0.8.and.dens.gt.0.5e5) then
     !   Caughlin & Fowler 1975
     !   t9a      = t9/(1.0 + 0.067*t9)
     !   nalambda = 1.26e27*t9a**(5.0/6.0)*t9**(-1.5d0)*exp(-84.165*t9a**(-1.0/3.0))/ &
     !        (exp(-0.01*t9a**4.0)+5.56e-3*exp(1.685*t9a**(2.0/3.0)))
     !   phi1dot  = 12.0*7.0/36.0*dens*nalambda*(1-phi1)**2.0*xc12init**2
        ! From Caughlin & Fowler 1988, ADNDT 40 283
        ! note this is UNSCREENED, should be screened at some point
        t9a = t9/(1.0+0.0396*t9)
        nasigmav = 4.27e26*t9a**(5.0/6)*t9**(-1.5)*exp(-84.165*t9a**(-1.0/3.0)-2.12e-3*t9**3)
        phi1dot = dens*xc12init*(1.0-phi1)**2*0.5/12.0*nasigmav
     endif

     !EXPLOSIVELY BURN CARBON WITH STEP FUNCTION: icburn = 0
     if(icburn.eq.0.and.t9.gt.4.d0) then 
        !9901 format(A20,4e16.8)
        !     write(*,9901) 'phi1,phi1dot,dphi1,dt = ',phi1,phi1dot,phi1dot*dt,t9
        phi1dot = dti
     endif
     phi1    = max(min(phi1_i + phi1dot*dt,1.0),0.e0)
  endif

 

  !--------------------------
  ! STAGE 2: NSQE RELAXATION
  !---------------------------
  phi2dot = 0.e0
  testExp = 182.06e9/temp_finnse - 46.054e0
  maxExp = log(HUGE(1.0)/10.0)
  if (testExp .GE. maxExp) then
     tau_nsqe = HUGE(1.0)
  else
     tau_nsqe = exp(testExp)
  endif
  if (tau_nsqe < dt) then
     phi2dot = (phi1-phi2)*dti
  else
     phi2dot = (phi1-phi2)/tau_nsqe
  endif
  ! slow down approach to NSQE except at the beginning of the simulation
  ! (REMOVED)
  !phi2dot = sign(min(abs(phi2dot),pbDphiMax*dti),phi2dot)
  !phi2dot = sign(min(abs(phi2dot),1.0*dti),phi2dot)

   !---- UPDATE PHI2 -----
  phi2 = max(min(phi2_i + phi2dot*dt,phi1),0.0e0)
  phi2dot = (phi2-phi2_i)*dti
  ! evolve already present ash toward NSE qbar and 1/abar on the NSQE timescale
  if (tau_nsqe < dt) then
     qbar_oldash  = qbar_finnse
     sumyi_oldash = sumyi_finnse
  else
     qbar_oldash  = qbar_oldash_i + (qbar_finnse-qbar_oldash_i)/tau_nsqe*dt
     sumyi_oldash = sumyi_oldash_i + (sumyi_finnse-sumyi_oldash_i)/tau_nsqe*dt
  endif


 

  !EXPLOSIVELY BURN OXYGEN WITH STEP FUNCTION 
  !if(t9 .ge. 1.d0) then 
  !   phi2dot = dti
  !endif


  !--------------------------------------------------------------------------
  !STAGE 3: NSE RELAXATION:
  ! silicon into nickel, start neutronization
  !--------------------------------------------------------------------------
  ! ajk     tau_nse    = (dens**0.2e0)*exp(179.7e9/temp - 40.5e0)
  ! fang1   tau_nse    = (dens**0.2e0)*exp(207.76e9/temp - 47.262e0)
  ! lbr     tau_nse    = exp(196.02e9/temp_finnse - 41.645e0)
  !--------------------------------------------------------------------------
  phi3dot = 0.d0
  testExp = 196.02e9/temp_finnse - 41.645e0
  if (testExp .GE. maxExp) then
     tau_nse = HUGE(1.0)
  else
     tau_nse = exp(testExp)
  endif
  if (tau_nse < dt) then
     phi3dot = (phi2-phi3)*dti
  else
     phi3dot = (phi2-phi3)/tau_nse
  endif
  !phi3dot = sign(min(abs(phi3dot),pbDphiMax*dti),phi3dot)
  !---- UPDATE PHI3 -----
  phi3 = max(min(phi3_i + phi3dot*dt,phi2),0.0e0)
  phi3dot = (phi3- phi3_i)*dti


  !INSTANTLY RELAX TO NSE
  !if(t9 .gt. 1) then
  !   phi3dot = dti
  !endif
  



  !--------------------------------------------------------------------------------
  ! Update abundances by mixing NS(Q)E ash and unburned material.
  ! Scale so that NSE abundances carry entire neutron excess and
  ! other material is all Ye=0.5.
  ! Also include neutronization here...
  ! can't do if phi2 is too small (we're dividing by phi2 in several places here)
  !--------------------------------------------------------------------------------
  if (phi2 > ash_update_threshold ) then
     qbar_ash  = (phi2_i*qbar_oldash + phi2dot*dt*qbar_finnse)/phi2
     sumyi_ash = (phi2_i*sumyi_oldash + phi2dot*dt*sumyi_finnse)/phi2
     ye_oldash = ye_oldash_i + phi3*yedot*dt
     ye_ash    = (phi2_i*ye_oldash +phi2dot*dt*0.5)/phi2
     
     if (ye_ash <0.42) then
        if (pbNseWarnLevel >= 2) then
           write(6,*)'pbnseBurn: low ye_ash in evolution: rho, temp', dens, temp
           write(6,*)'   phi1 phi2 phi3', phi1, phi2, phi3
           write(6,*)'   ye_ash, ye_i, ye_oldash_i,  yedot', &
                ye_ash, ye_i, ye_oldash_i, yedot
           write(6,*)'   changing ye_ash to 0.42'
        endif
        ye_ash = 0.42
     endif
  else
     yedot     = 0.e0
     edot      = 0.e0
     qbar_ash  = qbar_finnse
     sumyi_ash = sumyi_finnse
     ye_ash    = ye_i
  endif

  ! a sanity check because we've had some weird behavior lately (2007/04/02)
  if (sumyi_ash < 0.e0) then
!     write (6,*) 'sumyi_ash < 0   sumyi_ash phi1 phi2 phi3', phi1, phi2, phi3
     sumyi_ash = 1.e-2 ! really shouldn't be much less than 1/56
  endif

  ! Construct final state from carbon burn and NS(Q)E pieces
  qbar = (1.e0-phi1)*pbEbC12/pbAC12*xc12init &
       + (1.e0-phi2)*(1.e0-xc12init)*pbEbO16/pbAO16 &
       + (phi1-phi2)*xc12init*pbEbMg24/pbAMg24 &
       + phi2*qbar_ash

  sumyi = (1.e0-phi1)/pbAC12*xc12init &
       + (1.e0-phi2)*(1.e0-xc12init)/pbAO16 &
       + (phi1-phi2)*xc12init/pbAMg24 &
       + phi2*sumyi_ash
  
  ye = (1.e0-phi2)*0.5 + phi2*ye_ash
  
  qdot = 9.6485e17*(qbar - qbar_i)*dti  &
       - phi3*( yedot*pbNA*pbClight*pbClight*(pbMp + pbElectronMass - pbMn) + edot )
  
  qdot = pbFqdotMult * qdot

  
  !DONE
end subroutine bn_paraBurn
