!!****if* source/physics/Flame/FlameSpeed/Flame_laminarSpeed
!!
!! NAME
!!
!!  Flame_laminarSpeed
!!
!! SYNOPSIS
!!
!!  Flame_laminarSpeed(real, intent(in)  :: dens,
!!                               real, intent(in)  :: pres,
!!                               real, intent(out)  :: s)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   dens : density 
!!
!!   pres : pressure
!!
!!   s : 
!!
!! AUTOGENROBODOC
!!
!!
!!***


!------------------------------------------------------------------------

subroutine Flame_laminarSpeed(dens, pres, s)

  use Flame_data,ONLY : fl_useConstFlameSpeed,fl_constFspeed,&
                             II_CSPEED,II_CPOWER,II_CDENSI,II_C1,II_C2,II_C3,II_C4

  implicit none

  real, intent(in)   :: dens
  real, intent(in)   :: pres
  real, intent(out)  :: s

  real, parameter    :: c1 = -43.0e0
  real, parameter    :: c2 =  4.534e0
  real, parameter    :: c3 = -0.08333e0

  real               :: ldens

  if (fl_useConstFlameSpeed) then
     s = fl_constFspeed
  else
     ldens = log(dens)
     s = exp(c1 + c2*ldens + c3*(ldens**2)) 
  endif

end subroutine Flame_laminarSpeed

!------------------------------------------------------------------------

