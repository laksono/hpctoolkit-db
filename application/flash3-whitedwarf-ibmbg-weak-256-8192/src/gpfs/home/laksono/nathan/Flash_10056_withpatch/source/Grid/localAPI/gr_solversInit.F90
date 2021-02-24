!!****if* source/Grid/localAPI/gr_solversInit
!!
!! NAME
!!
!!  gr_solversInit
!!
!! 
!! SYNOPSIS
!!
!!  gr_solversInit()
!!
!!
!! DESCRIPTION
!!
!!  This routine initializes all grid solvers; solvers in use will
!!    have something more exciting than a stub
!!
!!***

subroutine gr_solversInit()    

  use Grid_data, ONLY : gr_isolatedBoundaries
  use gr_pfftInterface, ONLY : gr_pfftInit
  use gr_hgInterface, ONLY : gr_hgInit, gr_hgPfftInit

  implicit none 

  ! parallel FFT
  call gr_pfftInit()

  ! Multipole
  call gr_mpoleInit()

  ! Multigrid
  call gr_hgInit()

  ! Isolated multipole
  if (gr_isolatedBoundaries) then
     call gr_isoMpoleInit() 
  end if

  ! Parallel FFT at coarse level in Multigrid solver.
  call gr_hgPfftInit()

  return
end subroutine gr_solversInit     
