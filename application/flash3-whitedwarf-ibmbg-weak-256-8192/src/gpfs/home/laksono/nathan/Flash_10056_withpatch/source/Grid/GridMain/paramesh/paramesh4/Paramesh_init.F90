!!****if* source/Grid/GridMain/paramesh/paramesh4/Paramesh_init
!!
!! NAME
!!  Paramesh_init
!! 
!! SYNOPSIS
!!  call Paramesh_init()
!!
!! DESCRIPTION
!!  Simple interface to mesh initialization routine.
!!  Each mesh package will have a different copy of this function.
!!
!!  This Paramesh3 / Paramesh4 implementation:
!!  does some geometry-related checks and (if possible) readjusts
!!  geometry-related paramesh variables to the geometry requested by
!!  Flash runtime parameter.
!!
!!  Paramesh runtime parameters (these are different from Flash runtime
!!  parameters!) will be dumped to a file.
!!
!! USED BY
!!  Grid_init
!!***
subroutine Paramesh_init()
  use paramesh_interfaces , ONLY : amr_initialize
  use Driver_interface, ONLY: Driver_abortFlash
  use physicaldata, ONLY: curvilinear, cartesian_pm, cylindrical_pm, &
       spherical_pm, polar_pm
!  use physicaldata, ONLY: interp_mask_unk
!  use workspace, ONLY: interp_mask_work
  use Grid_data, ONLY: gr_geometry
  implicit none

#include "Flash.h"
! Undefine these symbols if defined by Flash.h, we want the ones from constants.h! - KW
#ifdef CARTESIAN
#undef CARTESIAN
#endif
#ifdef SPHERICAL
#undef SPHERICAL
#endif
#ifdef CYLINDRICAL
#undef CYLINDRICAL
#endif
#ifdef POLAR
#undef POLAR
#endif
#include "constants.h"


  logical :: testvar

  call amr_initialize()

#ifdef LIBRARY
  select case (gr_geometry)
     case (CARTESIAN)
        if (curvilinear) then
           cartesian_pm = .true.
        else
           cartesian_pm = .false.
        end if
        spherical_pm = .false.
        cylindrical_pm = .false.
        polar_pm = .false.
     case (SPHERICAL)
        cartesian_pm = .false.
        spherical_pm = .true.
        cylindrical_pm = .false.
        polar_pm = .false.
     case (CYLINDRICAL)
        cartesian_pm = .false.
        spherical_pm = .false.
        cylindrical_pm = .true.
        polar_pm = .false.
     case (POLAR)
        cartesian_pm = .false.
        spherical_pm = .false.
        cylindrical_pm = .false.
        polar_pm = .true.
     end select
  call gr_amr_dump_runtime_parameters
#else
  call gr_amr_dump_runtime_parameters
  select case (gr_geometry)
     case (CARTESIAN)
        testvar = (cartesian_pm .or. .not.curvilinear)
     case (SPHERICAL)
        testvar = spherical_pm
     case (CYLINDRICAL)
        testvar = cylindrical_pm
     case (POLAR)
        testvar = polar_pm
     case default
        testvar = .false.
     end select
     if (.not. testvar) then
        call Driver_abortFlash('[Paramesh_init] The value of the geometry runtime parameter is incompatible    &
             & with the geometry compiled into the AMR code.')
     end if
#endif


  if (cylindrical_pm .or. spherical_pm .or. polar_pm) then
     if (.not. curvilinear) then
        call Driver_abortFlash('[Paramesh_init] Curvilinear support must be enabled in order to use the requested geometry!')
     end if
  end if


!! NOTE: Maybe this is the place to put stuff like this?
!  interp_mask_unk = 2
!  interp_mask_work = 2
!!  call init_flash_physicaldata()
end subroutine Paramesh_init
