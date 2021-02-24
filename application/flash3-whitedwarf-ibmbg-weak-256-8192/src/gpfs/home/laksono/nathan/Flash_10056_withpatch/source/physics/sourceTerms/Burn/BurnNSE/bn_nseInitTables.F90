!!****ih* source/physics/sourceTerms/Burn/BurnNSE/bn_nseInitTables
!!
!! NAME
!!
!!    bn_nseInitTables 
!!
!!
!! SYNOPSIS
!!
!!    call bn_nseInitTables(integer(IN) :: mype,
!!                        character(IN) :: prestablename, 
!!                        character(IN) :: denstablename)
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
!!    mype -- current local processor number
!!    prestablename -- filename containing the interpolated pressure
!!    denstablename -- filename containing the interpolated density
!!
!! NOTES
!!
!!***

subroutine bn_nseInitTables(mype, prestablename, denstablename)
  
  use bn_nseData, ONLY: nemq, emq_grid, ldens_grid, demq, dldens, nldens, &
       & d_dYe, d_Ye_grid, d_nYe, &
       & d_ltemp_tab, d_lqbar_tab, d_ledot_tab, d_lmYedot_tab, d_lAbar_tab, &
       &  nhmq, hmq_grid, lpres_grid, dhmq, dlpres, nlpres, &
       & p_dYe, p_Ye_grid, p_nYe, p_ltemp_tab, p_lqbar_tab, p_ledot_tab, &
       & p_lmYedot_tab, p_lAbar_tab

  use Driver_interface, ONLY : Driver_abortFlash
  use Logfile_interface, ONLY : Logfile_stampMessage
  
#include "constants.h"
  
  implicit none
  
  integer, intent(IN)           :: mype
  character (len=*), intent(IN) :: prestablename, denstablename
  integer :: i, j, k, istat
  real    :: rtemp1, rtemp2, rtemp3, rtemp4, rtemp5
  real    :: rtemp6, rtemp7, rtemp8, rtemp9


  ! -----------------------------!
  !  Tables of NSE final state   !
  !------------------------------!

  !------------------------!
  ! READ IN PRESSURE TABLE !
  !------------------------!
  open (unit=21,file=prestablename,status='OLD',iostat=istat)
  if (istat /= 0) call Driver_abortFlash("Unable to open nse pressure table")

  read(21,*) p_nYe
  read(21,*) nlpres
  read(21,*) nhmq
  read(21,*) 
  
  ! space for coordinate grid
  allocate(p_Ye_grid(p_nYe),STAT=istat)
  if (istat /= 0) call Driver_abortFlash("Cannot allocate p_Ye_grid in InitNseTable")
  allocate(lpres_grid(nlpres),STAT=istat)
  if (istat /= 0) call Driver_abortFlash("Cannot allocate lpres_grid in InitNseTable")
  allocate(hmq_grid(nhmq),STAT=istat)
  if (istat /= 0) call Driver_abortFlash("Cannot allocate hmq_grid in InitNseTable")
  
  ! space for tables
  allocate(p_ltemp_tab(nhmq,nlpres,p_nYe),STAT=istat)
  if (istat /= 0) call Driver_abortFlash("Cannot allocate p_qbartab in InitNseTable")
  allocate(p_lqbar_tab(nhmq,nlpres,p_nYe),STAT=istat)
  if (istat /= 0) call Driver_abortFlash("Cannot allocate p_qbartab in InitNseTable")
  allocate(p_ledot_tab(nhmq,nlpres,p_nYe),STAT=istat)
  if (istat /= 0) call Driver_abortFlash("Cannot allocate p_edottab in InitNseTable")
  allocate(p_lmYedot_tab(nhmq,nlpres,p_nYe),STAT=istat)
  if (istat /= 0) call Driver_abortFlash("Cannot allocate p_Yedottab in InitNseTable")
  allocate(p_lAbar_tab(nhmq,nlpres,p_nYe),STAT=istat)
  if (istat /= 0) call Driver_abortFlash("Cannot allocate p_Abartab in InitNseTable")
  
  !! read the table, taking logs will be done separately
  do k = 1, p_nYe
     do j = 1, nlpres
        do i = 1, nhmq
           read(21,*) rtemp1, rtemp2, hmq_grid(i), rtemp4, p_ltemp_tab(i,j,k), &
                rtemp6, p_lqbar_tab(i,j,k), p_lAbar_tab(i,j,k), rtemp9, p_ledot_tab(i,j,k), &
                p_lmYedot_tab(i,j,k)
           
           ! avoid posttive Yedot
           if (p_lmYedot_tab(i,j,k).gt.0.0) p_lmYedot_tab(i,j,k) = -1.e-20
           
           ! empty portions of table are filled with values at the last good
           ! point up the column in hmq (actually done by propagating the value)
           if (p_lqbar_tab(i,j,k) == 0.0) then
              p_ltemp_tab(i,j,k) = p_ltemp_tab(i-1,j,k)
              p_lqbar_tab(i,j,k) = p_lqbar_tab(i-1,j,k)
              p_lAbar_tab(i,j,k) = p_lAbar_tab(i-1,j,k)
              p_ledot_tab(i,j,k) = p_ledot_tab(i-1,j,k)
              p_lmYedot_tab(i,j,k) = p_lmYedot_tab(i-1,j,k)
           endif
           
        enddo
        lpres_grid(j) = log10(rtemp2)
     enddo
     p_Ye_grid(k) = rtemp1
  enddo
  close(21)
  
  !----------------!
  ! CALCUALTE LOGS !
  !----------------!
  !! work with logs for temp, qbar, edot, -Yedot, abartables
  ! working with logs of Yedot. Change the sign to make it a positive
  ! quantity, then change it back after interpolating
  p_ltemp_tab = log10(p_ltemp_tab)
  p_lqbar_tab = log10(p_lqbar_tab)
  p_ledot_tab = log10(p_ledot_tab)
  p_lmYedot_tab = log10(-p_lmYedot_tab)
  p_lAbar_tab = log10(p_lAbar_tab)

  !-----------------!
  !CALCULATE DELTAS !
  !-----------------!
  !! store the deltas of grid for fast index calculation
  p_dYe = (p_Ye_grid(p_nYe)-p_Ye_grid(1))/(p_nYe-1)
  dlpres = (lpres_grid(nlpres)-lpres_grid(1))/(nlpres-1)
  dhmq = (hmq_grid(nhmq)-hmq_grid(1))/(nhmq-1)
  
  if (mype == MASTER_PE) then
     write(*,*)'Table of NSE final state at pressure initialized.'
  endif
  call Logfile_stampMessage(myPE,'[bn_nseInitTables] Table of NSE final state at pressure initialized.')


  !---------------------------!
  !  READ IN DENSITY TABLE    !
  !---------------------------!
  open (unit=21,file=denstablename,status='OLD',iostat=istat)
  if (istat /= 0) call Driver_abortFlash("Unable to open nse density table")
  
  read(21,*) d_nYe
  read(21,*) nldens
  read(21,*) nemq
  read(21,*) 
  
  ! space for coordinate grid
  allocate(d_Ye_grid(d_nYe),STAT=istat)
  if (istat /= 0) call Driver_abortFlash("Cannot allocate d_Ye_grid in InitNseTable")
  allocate(ldens_grid(nldens),STAT=istat)
  if (istat /= 0) call Driver_abortFlash("Cannot allocate ldens_grid in InitNseTable")
  allocate(emq_grid(nemq),STAT=istat)
  if (istat /= 0) call Driver_abortFlash("Cannot allocate emq_grid in InitNseTable")
  
  ! space for tables
  allocate(d_ltemp_tab(nemq,nldens,d_nYe),STAT=istat)
  if (istat /= 0) call Driver_abortFlash("Cannot allocate d_qbartab in InitNseTable")
  allocate(d_lqbar_tab(nemq,nldens,d_nYe),STAT=istat)
  if (istat /= 0) call Driver_abortFlash("Cannot allocate d_qbartab in InitNseTable")
  allocate(d_ledot_tab(nemq,nldens,d_nYe),STAT=istat)
  if (istat /= 0) call Driver_abortFlash("Cannot allocate d_edottab in InitNseTable")
  allocate(d_lmYedot_tab(nemq,nldens,d_nYe),STAT=istat)
  if (istat /= 0) call Driver_abortFlash("Cannot allocate d_Yedottab in InitNseTable")
  allocate(d_lAbar_tab(nemq,nldens,d_nYe),STAT=istat)
  if (istat /= 0) call Driver_abortFlash("Cannot allocate d_Abartab in InitNseTable")
  
  !! read the table data
  do k = 1, d_nYe
     do j = 1, nldens
        do i = 1, nemq
           read(21,*) rtemp1, rtemp2, emq_grid(i), rtemp4, d_ltemp_tab(i,j,k), &
                rtemp6, d_lqbar_tab(i,j,k), d_lAbar_tab(i,j,k), rtemp9, d_ledot_tab(i,j,k), &
                d_lmYedot_tab(i,j,k)
           
           ! avoid posttive Yedot
           if (d_lmYedot_tab(i,j,k).gt.0.0) d_lmYedot_tab(i,j,k) = -1.e-20
           
           ! empty portions of table are filled with values at the last good
           ! point up the column in hmq (actually done by propagating the value)
           if (d_lqbar_tab(i,j,k) == 0.0) then
              d_ltemp_tab(i,j,k) = d_ltemp_tab(i-1,j,k)
              d_lqbar_tab(i,j,k) = d_lqbar_tab(i-1,j,k)
              d_lAbar_tab(i,j,k) = d_lAbar_tab(i-1,j,k)
              d_ledot_tab(i,j,k) = d_ledot_tab(i-1,j,k)
              d_lmYedot_tab(i,j,k) = d_lmYedot_tab(i-1,j,k)
           endif
           
        enddo
        ldens_grid(j) = log10(rtemp2)
     enddo
     d_Ye_grid(k) = rtemp1
  enddo
  close(21)
  

  !----------------!
  ! CALCUALTE LOGS !
  !----------------!
  ! work with logs for temp, qbar, edot, -Yedot, abartables
  ! working with logs of Yedot. Change the sign to make it a positive
  ! quantity, then change it back after interpolating

  d_ltemp_tab   = log10(d_ltemp_tab)
  d_lqbar_tab   = log10(d_lqbar_tab)
  d_ledot_tab   = log10(d_ledot_tab)
  d_lmYedot_tab = log10(-d_lmYedot_tab)
  d_lAbar_tab   = log10(d_lAbar_tab)
  
  !-----------------!
  !CALCULATE DELTAS !
  !-----------------!
  !! store the deltas of grid for fast index calculation
  d_dYe  = (d_Ye_grid(d_nYe)-d_Ye_grid(1))/(d_nYe-1)
  dldens = (ldens_grid(nldens)-ldens_grid(1))/(nldens-1)
  demq   = (emq_grid(nemq)-emq_grid(1))/(nemq-1)
  
  
  if (mype == MASTER_PE) then
     write(*,*)'[bn_nseInitTables] Table of NSE final state at density initialized.'
  endif
  call Logfile_stampMessage(myPE,'[bn_nseInitTables] Table of NSE final state at density initialized.')

  return
end subroutine bn_nseInitTables
