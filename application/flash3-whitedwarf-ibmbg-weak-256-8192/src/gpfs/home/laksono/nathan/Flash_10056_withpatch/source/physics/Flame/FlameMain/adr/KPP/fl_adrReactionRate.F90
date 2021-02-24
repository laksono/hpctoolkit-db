
subroutine fl_adrReactionRate(s, phi, phidot)
  use fl_adrData, ONLY : fl_adrEps0, fl_adrEps1, fl_adrTauInv, fl_adrKppFact
  implicit none

  real, intent(in)  :: s
  real, intent(in)  :: phi
  real, intent(out) :: phidot


  phidot = fl_adrKppFact * fl_adrTauInv * s * (1.e0+fl_adrEps1-phi)*(phi-fl_adrEps0)
  
  return
  
end subroutine fl_adrReactionRate
