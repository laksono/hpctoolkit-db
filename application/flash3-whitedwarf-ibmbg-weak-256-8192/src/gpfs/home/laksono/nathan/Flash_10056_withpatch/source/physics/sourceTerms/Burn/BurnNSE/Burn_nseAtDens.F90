!!****ih* source/physics/sourceTerms/Burn/BurnNSE/Burn_nseAtDens
!!
!! NAME
!!
!!    Burn_nseAtDens
!!
!!
!! SYNOPSIS
!!
!!    Burn_nseAtDens( real(OUT)   :: qbar_nse,
!!                          real(OUT)   :: sumyi_nse,
!!                          real(OUT)   :: approxtemp,
!!                          real(OUT)   :: edot,
!!                          real(OUT)   :: Yedot,
!!                          real(IN)    :: Ye,
!!                          real(IN) :: dens,
!!                          real(INOUT) :: emq)
!!
!!
!! DESCRIPTION
!!  
!! This function gives the NSE final state for a constant density burn
!! via table lookup.  qbar is returned in MeV/nucleon.
!! Only the compositional information (qbar, 1/Abar, edot, Yedot) should be
!! used for hydrodynamics purposes.
!! Accurate thermodynamic prperties of the NSE final state (rho, T) should be
!! obtained by solving
!!        eint - qbar_nse = emq .
!! This is important due to the limited accuracy of the table interpolation
!! (i.e., the interpolated values cannot satisfy any constraints on their own.)
!!
!! ARGUMENTS
!!
!!    qbar_nse    --  Average binding energy per nucleon (averaged over nse)
!!    sumyi_nse   --  Sum over the abundances. (Y = X/A = X/(#of nucleons))
!!    approxtemp  --  ?? Ask Dean
!!    edot        --  Neutrino and anti-neutrino loss rates [ergs/g]
!!    Yedot       --  Change in electron fraction wrt time (neutronization rate)
!!    Ye          --  Electron fraction
!!    dens        --  Density in/out
!!    emq         --  E-Q = Internal energy minus binding energy (per volume?)
!!
!! NOTES
!!   Takes place after burning is complete.  NSE expands after burning
!!   and the equilibrium state is a function of the density and the
!!   temperature.  Been burned, flame is done, but then the hot stuff
!!   expands.  This changes the composition and releases energy.
!!
!!   LBR changed usage after discussion with Casey of a possible bug.
!!   There was not symmetry between the nseAtDens and the nseAtPres routines.
!!   Casey said it was irrelevant to results because they are really the same grid.
!!     p_nYe -> d_nYe
!!     p_Ye_grid -> d_Ye_grid
!!     p_dYe -> d_dYe
!!
!!***


subroutine Burn_nseAtDens(qbar_nse,sumyi_nse,approxtemp,edot,Yedot, Ye, dens, emq)
  
  use bn_nseData, ONLY: nemq, emq_grid, ldens_grid, demq, dldens, nldens, &
       & d_dYe, d_Ye_grid, d_nYe, &
       & d_ltemp_tab, d_lqbar_tab, d_ledot_tab, d_lmYedot_tab, d_lAbar_tab
  
  implicit none
  
  real, intent(IN) :: Ye, dens, emq
  real, intent(OUT) :: qbar_nse,sumyi_nse,approxtemp,edot,Yedot
  
  integer :: emq_a, ldens_a, Ye_a
  
  real :: ldens, emq_v
  real :: te111, te211, te121, te112, te221, te212, te122, te222
  real :: qb111, qb211, qb121, qb112, qb221, qb212, qb122, qb222
  real :: ed111, ed211, ed121, ed112, ed221, ed212, ed122, ed222
  real :: yd111, yd211, yd121, yd112, yd221, yd212, yd122, yd222
  real :: ab111, ab211, ab121, ab112, ab221, ab212, ab122, ab222
  
  real :: c1, c2, c3, abar

  
  !! find the location in the table grid
  !! remember that we use the log of density
  
  emq_v = emq
  ldens = log10(dens) 

  
  ! we want to truncate at the "bottom" edges (in density and emq)
  !instead of extrapolate.  Should be pure nickel down here.
  ! emq is stored in drecrasing order

  if ( emq_v < emq_grid(nemq)   )   emq_v = emq_grid(nemq)
  if ( ldens < ldens_grid(1) ) ldens = ldens_grid(1)

  ! find v_a  such that v_grid(v_a) <= v < v_grid(v_a+1)
  emq_a   = floor((emq_v - emq_grid(1))/demq)+1
  ldens_a = floor((ldens - ldens_grid(1))/dldens)+1
  Ye_a   = floor((Ye - d_Ye_grid(1))/d_dYe)+1

  ! at upper edges we will extrapolate
  if (emq_a < 1) emq_a = 1
  if (emq_a >= nemq) emq_a = nemq - 1
  ! lower limit handled above
  if (ldens_a >= nldens) ldens_a = nldens - 1

  ! extrapolate both ways in Ye
  if (Ye_a >= d_nYe) Ye_a = d_nYe - 1
  if (Ye_a < 1) Ye_a = 1

  !! calculate the coefficients
  
  c1 = ( emq_v - emq_grid(emq_a) )/                &
       ( emq_grid(emq_a+1) - emq_grid(emq_a) )
  c2 = ( ldens - ldens_grid(ldens_a)) /              &
       ( ldens_grid(ldens_a+1) - ldens_grid(ldens_a) )
  c3 = ( Ye - d_Ye_grid(Ye_a) )/                  &
       ( d_Ye_grid(Ye_a+1) - d_Ye_grid(Ye_a) )

  !! build the local cubes
  te111 = d_ltemp_tab(emq_a,ldens_a,Ye_a)
  te211 = d_ltemp_tab(emq_a+1,ldens_a,Ye_a)
  te121 = d_ltemp_tab(emq_a,ldens_a+1,Ye_a)
  te112 = d_ltemp_tab(emq_a,ldens_a,Ye_a+1)
  te221 = d_ltemp_tab(emq_a+1,ldens_a+1,Ye_a)
  te212 = d_ltemp_tab(emq_a+1,ldens_a,Ye_a+1)
  te122 = d_ltemp_tab(emq_a,ldens_a+1,Ye_a+1)
  te222 = d_ltemp_tab(emq_a+1,ldens_a+1,Ye_a+1)
  
  qb111 = d_lqbar_tab(emq_a,ldens_a,Ye_a)
  qb211 = d_lqbar_tab(emq_a+1,ldens_a,Ye_a)
  qb121 = d_lqbar_tab(emq_a,ldens_a+1,Ye_a)
  qb112 = d_lqbar_tab(emq_a,ldens_a,Ye_a+1)
  qb221 = d_lqbar_tab(emq_a+1,ldens_a+1,Ye_a)
  qb212 = d_lqbar_tab(emq_a+1,ldens_a,Ye_a+1)
  qb122 = d_lqbar_tab(emq_a,ldens_a+1,Ye_a+1)
  qb222 = d_lqbar_tab(emq_a+1,ldens_a+1,Ye_a+1)
  
  ed111 = d_ledot_tab(emq_a,ldens_a,Ye_a)
  ed211 = d_ledot_tab(emq_a+1,ldens_a,Ye_a)
  ed121 = d_ledot_tab(emq_a,ldens_a+1,Ye_a)
  ed112 = d_ledot_tab(emq_a,ldens_a,Ye_a+1)
  ed221 = d_ledot_tab(emq_a+1,ldens_a+1,Ye_a)
  ed212 = d_ledot_tab(emq_a+1,ldens_a,Ye_a+1)
  ed122 = d_ledot_tab(emq_a,ldens_a+1,Ye_a+1)
  ed222 = d_ledot_tab(emq_a+1,ldens_a+1,Ye_a+1)
  
  yd111 = d_lmYedot_tab(emq_a,ldens_a,Ye_a)
  yd211 = d_lmYedot_tab(emq_a+1,ldens_a,Ye_a)
  yd121 = d_lmYedot_tab(emq_a,ldens_a+1,Ye_a)
  yd112 = d_lmYedot_tab(emq_a,ldens_a,Ye_a+1)
  yd221 = d_lmYedot_tab(emq_a+1,ldens_a+1,Ye_a)
  yd212 = d_lmYedot_tab(emq_a+1,ldens_a,Ye_a+1)
  yd122 = d_lmYedot_tab(emq_a,ldens_a+1,Ye_a+1)
  yd222 = d_lmYedot_tab(emq_a+1,ldens_a+1,Ye_a+1)
  
  ab111 = d_lAbar_tab(emq_a,ldens_a,Ye_a)
  ab211 = d_lAbar_tab(emq_a+1,ldens_a,Ye_a)
  ab121 = d_lAbar_tab(emq_a,ldens_a+1,Ye_a)
  ab112 = d_lAbar_tab(emq_a,ldens_a,Ye_a+1)
  ab221 = d_lAbar_tab(emq_a+1,ldens_a+1,Ye_a)
  ab212 = d_lAbar_tab(emq_a+1,ldens_a,Ye_a+1)
  ab122 = d_lAbar_tab(emq_a,ldens_a+1,Ye_a+1)
  ab222 = d_lAbar_tab(emq_a+1,ldens_a+1,Ye_a+1)

  
  !! now interpolate

   approxtemp  =    c3                          &
       *( (1.0-c1)*(1.0-c2)    &
       *te112    &
       +      c1*(1.0-c2)    &
       *te212    &
       +      c2*(1.0-c1)    &
       *te122    &
       +      c1*c2          &
       *te222    &
       )                      &
       +(1.0-c3)                    &
       *( (1.0-c1)*(1.0-c2)    &
       *te111    &
       +      c1*(1.0-c2)    &
       *te211    &
       +      c2*(1.0-c1)    &
       *te121    &
       +      c1*c2          &
       *te221    &
       )
  
  approxtemp = 10.e0**approxtemp
  
  qbar_nse  =    c3                          &
       *( (1.0-c1)*(1.0-c2)    &
       *qb112    &
       +      c1*(1.0-c2)    &
       *qb212    &
       +      c2*(1.0-c1)    &
       *qb122    &
       +      c1*c2          &
       *qb222    &
       )                      &
       +(1.0-c3)                    &
       *( (1.0-c1)*(1.0-c2)    &
       *qb111    &
       +      c1*(1.0-c2)    &
       *qb211    &
       +      c2*(1.0-c1)    &
       *qb121    &
       +      c1*c2          &
       *qb221    &
       )
  
  qbar_nse = 10.e0**qbar_nse
  
  edot  =    c3                          &
       *( (1.0-c1)*(1.0-c2)    &
       *ed112    &
       +      c1*(1.0-c2)    &
       *ed212    &
       +      c2*(1.0-c1)    &
       *ed122    &
       +      c1*c2          &
       *ed222    &
       )                      &
       +(1.0-c3)                    &
       *( (1.0-c1)*(1.0-c2)    &
       *ed111    &
       +      c1*(1.0-c2)    &
       *ed211    &
       +      c2*(1.0-c1)    &
       *ed121    &
       +      c1*c2          &
       *ed221    &
       )
  
  edot = 10.e0**edot
  
  Yedot  =    c3                         &
       *( (1.0-c1)*(1.0-c2)    &
       *yd112    &
       +      c1*(1.0-c2)    &
       *yd212    &
       +      c2*(1.0-c1)    &
       *yd122    &
       +      c1*c2          &
       *yd222    &
       )                      &
       +(1.0-c3)                    &
       *( (1.0-c1)*(1.0-c2)    &
       *yd111    &
       +      c1*(1.0-c2)    &
       *yd211    &
       +      c2*(1.0-c1)    &
       *yd121    &
       +      c1*c2          &
       *yd221    &
       )
  
  ! note that Yedot is a negative quantity
  Yedot = - 10.e0**Yedot
  
  abar  =    c3                          &
       *( (1.0-c1)*(1.0-c2)    &
       *ab112    &
       +      c1*(1.0-c2)    &
       *ab212    &
       +      c2*(1.0-c1)    &
       *ab122    &
       +      c1*c2          &
       *ab222    &
       )                      &
       +(1.0-c3)                    &
       *( (1.0-c1)*(1.0-c2)    &
       *ab111    &
       +      c1*(1.0-c2)    &
       *ab211    &
       +      c2*(1.0-c1)    &
       *ab121    &
       +      c1*c2          &
       *ab221    &
       )
  
  abar = 10.e0**abar

  sumyi_nse = 1.0/abar

  !write (6,*) 'returning qbar_nse,sumyi_nse,approxtemp,edot,Yedot', &
  !              qbar_nse,sumyi_nse,approxtemp,edot,Yedot
  return
end subroutine Burn_nseAtDens
