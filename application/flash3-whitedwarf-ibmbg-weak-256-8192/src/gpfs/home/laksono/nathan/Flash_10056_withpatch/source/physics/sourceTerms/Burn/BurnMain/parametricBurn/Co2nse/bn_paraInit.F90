!!****if* source/physics/sourceTerms/Burn/BurnMain/parametricBurn/Co2nse/bn_paraInit
!!
!! NAME
!!
!!  bn_paraInit
!!
!! SYNOPSIS
!!
!!  bn_paraInit(  integer, intent(in)     :: myPE)
!!
!! DESCRIPTION
!!
!!  Initialize paraBurn Data.
!!  C.A.Meakin, 26 Feb. 2007 (Chicago)
!!
!! ARGUMENTS
!!
!!   myPE:  current processor
!!
!!***
!
!

subroutine bn_paraInit(myPE)

  use bn_paraData 

  use Driver_interface, ONLY : Driver_abortFlash
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, ONLY : Logfile_stampMessage
  use PhysicalConstants_interface, ONLY:  PhysicalConstants_get

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent(IN)     :: myPE
  integer, parameter      :: SPEC_UNIT=43,SPEC_NUM=238,SPEC_TO_READ=5
  character(len=4)        :: isotopeName
  real                    :: abar,zbar,bindEnergy,dummy
  integer                 :: i,isotope,count

  !WHAT IS THIS?
  call RuntimeParameters_get("pbFqdotMult", pbFqdotMult)
  call RuntimeParameters_get("pbUseConstComposition", pbUseConstComposition)

  !SETUP PARAMS FOR EVOLING NSE STATE
  call RuntimeParameters_get("eos_tolerance", pbTol)
  !call RuntimeParameters_get("eos_maxNewton", pbMaxNewton)
  !call RuntimeParameters_get("irenorm", pbRenorm)
  !.. find parameters for limiting
  call RuntimeParameters_get("smallx", pbSmallx)
  !Dev : find out the default value for this parameter
  call RuntimeParameters_get("pbNseWarnLevel", pbNseWarnLevel)
  call RuntimeParameters_get("pbReact", pbReact)
  call RuntimeParameters_get("c_frac", pbC12Init)

  !PHYSICAL CONSTANTS
  call PhysicalConstants_get("electron mass", pbElectronMass)
  call PhysicalConstants_get("Avogadro", pbNA)
  call PhysicalConstants_get("speed of light", pbClight)
  call PhysicalConstants_get("proton mass", pbMp)
  pbMn = 1.67492716e-24 !neutron mass

  !READ IN SOME RELEVANT NUCLIDE DATA
  open(unit=SPEC_UNIT,file="SpeciesList.txt")
  count=0
  i=0
  do while((count<SPEC_TO_READ).and.(i<=SPEC_NUM))
     i=i+1
     read(SPEC_UNIT,*)isotopeName,zbar,abar,dummy,bindEnergy
     if (trim(isotopeName) .eq. 'mg24') then
        pbAMg24 = abar
        pbEbMg24 = bindEnergy
        count = count + 1
     else if (trim(isotopeName) .eq. 'c12') then
        pbAC12 = abar
        pbEbC12 = bindEnergy
        count = count + 1
     else if (trim(isotopeName) .eq. 'o16') then
        pbAO16 = abar
        pbEbO16 = bindEnergy
        count = count + 1
     else if (trim(isotopeName) .eq. 'si28') then
        pbASi28 = abar
        pbEbSi28 = bindEnergy
        count = count + 1
     else if (trim(isotopeName) .eq. 'ni56') then
        pbANi56 = abar
        pbEbNi56 = bindEnergy
        count = count + 1
     end if
  end do
  close(unit=SPEC_UNIT)


  !CALCULATE SOME ADDITIONAL BURNING RELATED QUANTITIES
  !(find energy release during burning, nsqe, and nse stages) 
  pbMev2ergMp = pbMev2erg/pbMp
  pbQc2mg = (pbEbMg24/pbAMg24 - pbEbC12/pbAC12) * pbMev2ergMp
  pbQo2si  = (pbEbSi28/pbASi28 - pbEbO16/pbAO16) * pbMev2ergMp
  pbQmg2si = (pbEbSi28/pbASi28 - pbEbMg24/pbAMg24) * pbMev2ergMp
  pbQsi2ni = (pbEbNi56/pbANi56 - pbEbSi28/pbASi28) * pbMev2ergMp



end subroutine bn_paraInit

