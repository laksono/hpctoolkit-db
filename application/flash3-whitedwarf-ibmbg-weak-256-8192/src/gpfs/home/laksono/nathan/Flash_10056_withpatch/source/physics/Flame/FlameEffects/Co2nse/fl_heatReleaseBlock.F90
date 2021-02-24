
!------------------------------------------------------------------------

subroutine fl_heatReleaseBlock(q, flag)

  use bn_paraData, ONLY : pbQc2mg, pbQmg2si, pbQo2si, pbQsi2ni, pbC12Init
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none
#include "Flash.h"

  integer, parameter :: BURN = 1
  integer, parameter :: NSQE = 2
  integer, parameter :: NSE  = 3

  real, dimension(:,:,:),   intent(out) :: q
  integer,                  intent(in)  :: flag

  select case(flag)
  case(BURN)
     q = pbQc2mg * pbC12Init
  case(NSQE)
     q = pbQo2si * (1.e0-pbC12Init) + pbQmg2si * pbC12Init
  case(NSE)
     q = pbQsi2ni
  case default
     call Driver_abortFlash("heatRelease: invalid input parameter")
  end select

end subroutine fl_heatReleaseBlock
