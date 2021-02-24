
subroutine Flame_rhJump(dens_u, pres_u, temp_u, ener_u, ye_u, sumy_u, &
                  dens_b, pres_b, temp_b, ener_b, ye_b, sumy_b, q, s, mode)

  use Flame_data, ONLY : fl_smallt, fl_smlrho
  use Eos_interface, ONLY : Eos
  
  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"  

  real,    intent(inout)  :: dens_u, pres_u, temp_u, ener_u
  real,    intent(out)    :: dens_b, pres_b, temp_b, ener_b
  real,    intent(in)     :: q, s
  integer, intent(in)     :: mode  !! This is the Eos mode
  real, intent(in)        :: ye_u, sumy_u 
  real, intent(in)        :: ye_b, sumy_b
  real, dimension(NSPECIES) :: xn_b, xn_u
  real                    :: dfd, dft, dgd, dgt, f, g
  real                    :: sq_s_dens, d_inv_dens, determinant_inv
  real                    :: dens_b_old, temp_b_old
  real                    ::  error
  integer                 :: niters
  real, dimension(EOS_NUM) :: eosData
  logical, dimension(EOS_VARS+1:EOS_NUM) :: mask 
  integer                  :: vecLen=1
1 format(5(2X, E10.4))

!------------------------------------------------------------------
! given: unburned state (dens_u, pres_u, velx_u),  s, q
!
 
! unburned state:
  eosData(EOS_DENS)=dens_u
  eosData(EOS_PRES)=pres_u
  eosData(EOS_TEMP)=temp_u
  eosData(EOS_EINT)=ener_u
  
     
!!!
!
! GCJIV
!
!!!

  eosData(EOS_ABAR) = 1.0e0/sumy_u
  eosData(EOS_ZBAR) = ye_u/sumy_u
  
  mask = .false.
  mask(EOS_DPT)=.true.
  mask(EOS_DPD)=.true.
  mask(EOS_DET)=.true.
  mask(EOS_DED)=.true.

  call Eos(mode,vecLen,eosData,xn_u,mask)


  pres_u = eosData(EOS_PRES)
  temp_u = eosData(EOS_TEMP)
  ener_u = eosData(EOS_EINT)


#if NSPECIES > 0
!!!
!
! Terrible!! Please fix to be general.
!
! GCJIV

  xn_b(:) = 1.0E-10
  xn_u(:) = 1.0E-10

  xn_b(1) = 1.0e0
  xn_u(1)  = 1.0e0
!
!!!
#endif


!  write (*,*) '------------------------------------------------------'
!  write (*,*) ''
!  write (*,*) 'rhjump: flame speed, heat release:'
!  write (*,1) s, q
!  write (*,*) 'rhjump: unburned temp, dens, ener, pres:'
!  write (*,1) temp_u, dens_u, ener_u, pres_u
!  write (*,*) ''

!.. burned state:

  dens_b = dens_u
!  eosData(EOS_DENS) = dens_b
  eosData(EOS_EINT) = ener_u + q
  eosData(EOS_ZBAR) = ye_b / sumy_b
  eosData(EOS_ABAR) = 1.0e0/sumy_b

  call Eos(MODE_DENS_EI,vecLen,eosData,xn_b,mask)


!!$  call eos (dens_b, temp_b, pres_b, ener_b, xn_b, dmy1, dmy2, dmy3, &
!!$            dpt, dpd, det, ded, & 
!!$            dmy8, dmy9, dmy10, dmy11, dmy12, dmy13, dmy14, dmy15, 2)

  ener_b = eosData(EOS_EINT)
  pres_b = eosData(EOS_PRES)
  temp_b = eosData(EOS_TEMP) 
  error = 1.
  niters = 0

  sq_s_dens = (s*dens_u)**2

  do while (error > 1.e-8 .and. niters < 100)


     d_inv_dens = 1./dens_b - 1./dens_u

     f = pres_b - pres_u + sq_s_dens * d_inv_dens 
     g = ener_b - ener_u - q + 0.5*(pres_b + pres_u) * d_inv_dens

     dfd = eosData(EOS_DPD) - sq_s_dens/dens_b**2
     dft = eosData(EOS_DPT)
     dgd = eosData(EOS_DED) + 0.5*d_inv_dens*eosData(EOS_DPD) - 0.5*(pres_b + pres_u) / dens_b**2
     dgt = eosData(EOS_DET) + 0.5*d_inv_dens*eosData(EOS_DPT)

     determinant_inv = 1./(dfd*dgt - dft*dgd)

     dens_b_old = dens_b
     temp_b_old = temp_b

     dens_b = dens_b - (f*dgt - g*dft) * determinant_inv
     temp_b = temp_b + (f*dgd - g*dfd) * determinant_inv

     if (dens_b .lt. fl_smlrho) then

        dens_b = 0.5*dens_b_old
        temp_b = temp_b_old
    
     elseif (temp_b .lt. temp_u) then

        temp_b = temp_u
        dens_b = dens_b_old

     endif

     eosData(EOS_DENS)=dens_b
     eosData(EOS_TEMP)=temp_b

     call Eos(MODE_DENS_TEMP,vecLen,eosData,xn_b,mask)

     pres_b=eosData(EOS_PRES)
     ener_b=eosData(EOS_EINT)

     error = abs(f/pres_u) + abs(g/ener_u)

!     write(*,1) temp_b, dens_b, ener_b, pres_b, error

     niters = niters + 1

  enddo

!  write (*,*) ''
!  write (*,*) 'rhjump: burned temp, dens, ener, pres:'
!  write (*,1) temp_b, dens_b, ener_b, pres_b
!  write (*,*) 'rhjump:   niters, error:', niters, error
!  write (*,*) ''
!  write (*,*) '-------------------------------------------------'
   if (niters >= 100) write (6,*) 'rhjump did not converge'

!------------------------------------------------------------------

end subroutine Flame_rhJump
