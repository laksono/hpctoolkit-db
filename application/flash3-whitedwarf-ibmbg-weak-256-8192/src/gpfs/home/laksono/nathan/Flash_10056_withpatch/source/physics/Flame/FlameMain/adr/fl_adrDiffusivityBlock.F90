
!! Dev : is there a need for this function ?
subroutine fl_adrDiffusivityBlock(s, D)

  use fl_adrData, ONLY : fl_adrDelta
  implicit none

  real, dimension(:,:,:), intent(in)  :: s
  real, dimension(:,:,:), intent(out) :: D

  D(:,:,:) = fl_adrDelta * s(:,:,:)

  return

end subroutine fl_adrDiffusivityBlock

