!#define DEBUG 1
!> @author  Vladimir Mironov
!
!> @brief This module contains subroutines for 1-electron integrals
!>  calculation. Based on `int1.src` code. It supports all calculation
!>  types except `RUNTYP=MOROKUMA`.
!
!> @note all subroutines in this file uses SP-free basis set which can
!>  be prepared simply with call to `SPLIT_SP_BASIS` from `MOD_NOSP_BASIS` module
!>  All low-level 1e-integral subroutines were moved to `MOD_1E_PRIMITIVES` module
!
!  REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!
MODULE omp_int1

    USE prec, ONLY: dp, fp
    USE mx_limits, ONLY: MXSH, MXGTOT, MXATM, MAXSH, MXABC, MAXDEN
    USE mod_nosp_basis, ONLY: basis_set, nosp_basis, split_sp_basis, &
            bas_norm_matrix, &
            bas_denorm_matrix
    USE mod_1e_primitives, ONLY: &
        update_triang_matrix, &
        density_ordered, &
        comp_coulomb_int1_prim, &
        comp_kin_ovl_int1_prim, &
        comp_lz_int1_prim, &
        comp_fmoesp_int1_prim, &
        comp_coulomb_dampch_int1_prim, &
        comp_coulpot_prim

    USE mod_shell_tools, ONLY: shell_t, shpair_t

    USE mod_gauss_hermite, ONLY: prepQuadGaussHermite
    IMPLICIT NONE

    REAL(KIND=fp), PARAMETER :: ZERO = 0.0_fp           !< Floating point zero
    REAL(KIND=fp), PARAMETER :: RLN10 = 2.30258_fp      !< \f$ ln(10) \f$

    INTEGER, PARAMETER :: MXCHRM = 1

    !<   size of shell pair block (square of max.num. basis functions in max.ang.m.)
    INTEGER, PARAMETER :: BLOCKSIZE = 28*28

    REAL(KIND=dp), PARAMETER :: &
        RLUT = transfer('LUT-IOTC', 1.0_dp)

    CHARACTER(LEN=*), PARAMETER :: STR_NONE = 'NONE'

    INTERFACE int1_coul_ext_chg
        MODULE PROCEDURE chg_ints_soa
        MODULE PROCEDURE chg_ints_aos
    END INTERFACE

    INTERFACE int1_coul
       MODULE PROCEDURE int1_coul_xyzc
       MODULE PROCEDURE int1_coul_x_y_z_c
       MODULE PROCEDURE int1_coul_xyz_c
    END INTERFACE

    INTERFACE int1_fmoesp
       MODULE PROCEDURE int1_fmoesp_ffcxyz
       MODULE PROCEDURE int1_fmoesp_fcall
    END INTERFACE

    PRIVATE
    PUBLIC omp_hst
    PUBLIC grid_fmo_ints
    PUBLIC omp_epoten
    PUBLIC omp_fmo_epoten

CONTAINS
!> @brief Driver for conventional h, S, and T integrals
!
!> @details  Compute one electron integrals and core Hamiltonian,
!>  - S is evaluated by Gauss-Hermite quadrature,
!>  - T is an overlap with -2,0,+2 angular momentum shifts,
!>  - V is evaluated by Gauss-Rys quadrature, then \f$ h = T+V \f$
!>  Also, do \f$ L_z \f$ integrals for atoms and linear cases.
!>  Also, do FMO ESP integrals if needed.
!>  This subroutine is capable to do integrals in parallel using
!>  both OpenMP and MPI. It it helpful when running large FMO jobs.
!
!> @note Based on `HSANDT` subroutine from file `INT1.SRC`
!
!> @author Vladimir Mironov
!
!   REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!>
!> @param[in,out]   h       one-electron Hamiltonian matrix in packet format
!> @param[in,out]   s       packed matrix of overlap integrals
!> @param[in,out]   t       packed matrix of kinetic energy integrals
!> @param[in,out]   z       packed matrix of z-angular momentum (Lz) integrals
!> @param[in,out]   esp1e   packed matrix of FMO ESP 1e integrals
!> @param[in]       nlct    an array for Divide-and-Conquer (DC) method screening
!> @param[in]       atchrg  array of atomic charges used in FMO
!> @param[in]       ll2     size of matrices `H`, `S`, `T`,..
!> @param[in]       nat1e   number of atoms for FMO runs
!> @param[in]       lfmoc   pointer to coordinates of FMO atoms in global (`FMCOM`) memory array
!> @param[in]       lzint   flag to compute Lz integrals
!> @param[in]       doesp0  flag to run FMO ESP0 job
!> @param[in]       doesp1  flag to run FMO ESP1 job
!> @param[in]       dccut   flag specifying DC job
!> @param[in]       dbug    flag for debug output
 SUBROUTINE omp_hst(h,s,t,z,esp1e,nlct,atchrg,ll2,nat1e,lfmoc,lzint,doesp0,doesp1,dccut,dbug)

    REAL(KIND=fp), CONTIGUOUS, INTENT(INOUT) :: h(:), s(:), t(:), z(:), esp1e(:) ! All dims. LL2
    INTEGER, INTENT(IN) :: nlct(:) ! All dims. LL2
    REAL(KIND=fp), INTENT(IN) :: atchrg(*)
    INTEGER, INTENT(IN) :: ll2, lfmoc, nat1e
    LOGICAL, INTENT(IN) :: dbug, lzint, doesp0, doesp1, dccut

    COMMON /CHMGMS/ xchm(mxchrm),ychm(mxchrm),zchm(mxchrm), &
                    dxelmm(mxchrm),dyelmm(mxchrm),dzelmm(mxchrm), &
                    qchm(mxchrm),nchmat,kchrmm
      REAL(KIND=fp) :: xchm,ychm,zchm,dxelmm,dyelmm,dzelmm,qchm
      INTEGER :: nchmat,kchrmm

    COMMON /CONV  / dentol,en,etot,ehf,ehf0,diff,iter,icalp,icbet
      REAL(KIND=fp) :: dentol,en,etot,ehf,ehf0,diff
      INTEGER :: iter,icalp,icbet

    COMMON /FMOINF/ nfg,nlayer,natfmo,nbdfg,naotyp,nbody
      INTEGER :: nfg,nlayer,natfmo,nbdfg,naotyp,nbody

    COMMON /INFOA / nat,ich,mul,num,nqmt,ne,na,nb, &
                    zan(mxatm),c(3,mxatm),ian(mxatm)
      INTEGER :: nat,ich,mul,num,nqmt,ne,na,nb,ian
      REAL(KIND=fp) :: zan,c

    COMMON /IOFILE/ ir,iw,ip,is,ipk,idaf,nav,ioda(950)
      INTEGER :: ir,iw,ip,is,ipk,idaf,nav,ioda

    COMMON /OUTPUT/ nprint,itol,icut,normf,normp,nopk
      INTEGER :: nprint,itol,icut,normf,normp,nopk

    COMMON /PAR   / me,master,nproc,ibtyp,iptim,goparr,dskwrk,maswrk
      INTEGER :: me,master,nproc,ibtyp,iptim
      LOGICAL :: goparr,dskwrk,maswrk

    COMMON /RUNOPT/ runtyp,exetyp,nevals,nglevl,nhlevl
      REAL(KIND=fp) :: runtyp,exetyp
      INTEGER :: nevals,nglevl,nhlevl

!    COMMON /SYMIND/ tol,ii,jj,lit,ljt,mini,minj,maxi,maxj,iandj
      REAL(KIND=fp) :: tol

    COMMON /RELWFN/ rmethod,qrqmt,clig,clig2,qrtol,tau, &
                    iqrord,modqr,nesoc,nratom, &
                    numu,nqmtr,nqrdaf,morda,ndarelb
      REAL(KIND=fp) :: rmethod,qrqmt,clig,clig2,qrtol,tau
      INTEGER :: iqrord,modqr,nesoc,nratom,numu,nqmtr,nqrdaf,morda,ndarelb

    COMMON /WFNOPT/ scftyp,vbtyp,dftype,tddftyp,cityp,cctyp,mplevl,mpctyp
      REAL(KIND=fp) :: scftyp,vbtyp,dftype,tddftyp,cityp,cctyp
      INTEGER :: mplevl,mpctyp

    COMMON /FMCOM / x(1)
      REAL(KIND=fp), TARGET :: x

    COMMON /TINOPT/ mparti,mmonly,qmmm
      INTEGER :: mparti
      LOGICAL :: mmonly, qmmm

    COMMON /ISEPS / iseps,useps
      LOGICAL :: iseps, useps

    COMMON /COSDAT/ se2,secorr,qvcosmo,elast,emp2cos,emp2last, &
                    cosvol,cossar,ediel,eoc1,deoc_rs,sumqsc, &
                    sumqscold,zsum,zsum2,zsum3,fepsi,rds,disex2, &
                    epsi,cosrad,disex,outchg,ediel_save, &
                    maxnps,icorr,itrip,nqs,mp2trip,mp2iter, &
                    icfreq,nspa,nsph,npsd,nps,nps2,nden,npspher, &
                    cosbug,coswrt,dcosmo,prfcnd,ioutch
      REAL(KIND=fp) :: se2,secorr,qvcosmo,elast,emp2cos,emp2last, &
          cosvol,cossar,ediel,eoc1,deoc_rs,sumqsc,sumqscold,zsum,zsum2,zsum3,fepsi,rds,disex2, &
          epsi,cosrad,disex,outchg,ediel_save
      INTEGER :: maxnps,icorr,itrip,nqs,mp2trip,mp2iter, &
          icfreq,nspa,nsph,npsd,nps,nps2,nden,npspher
      LOGICAL :: cosbug,coswrt,dcosmo,prfcnd,ioutch

    COMMON /COSVCE/ qscnet(mxabc),corzan(3,mxabc), &
                    qden(maxden),qscnet_save(mxabc),iatsp(mxabc+1)
      REAL(KIND=fp) :: qscnet,corzan,qden,qscnet_save
      INTEGER :: iatsp

    LOGICAL :: uncon, is_lut

    REAL(KIND=fp) :: dummy
    INTEGER :: l1, l2, nqmmatm

    REAL(KIND=fp), POINTER :: fmoc(:,:)

    CALL prepQuadGaussHermite
    ! TODO: find a way to skip unnecessary basis set recreation step in FMO
    !IF (nfg==0) CALL split_sp_basis
    CALL split_sp_basis

    uncon = transfer(rmethod,STR_NONE)/=STR_NONE .AND. mod(modqr,2)==1
    is_lut = (rmethod==RLUT)

    IF (is_lut) uncon = .FALSE.

    tol = RLN10*itol
!   norm = normf /= 1 .OR. normp /= 1

    IF (doesp1.AND.lzint) THEN
!      ZBLK is reused for some other purpose.
       WRITE(iw,*) 'STORAGE NOT ALLOCATED IN HSANDT'
       CALL abrt
    END IF

!   Exclude 1e potential in ESDIM, because density is used,
!   not point charges

!   ----- MOPAC INTEGRALS ARE DONE ELSEWHERE -----
    IF (transfer(mpctyp,STR_NONE)/=STR_NONE) THEN
       CALL mpcint
       RETURN
    END IF

!   For COSMO direct SCF, we need the difference of H between the SCF
!   iterations in RHFCL, so hold is copied into section 87
!   and HNEW is in section 11 as usual, beginning in SCF-cycle 2
    IF (iseps .AND. iter>0) THEN
       CALL daread(idaf,ioda,h,ll2,11,0)
       CALL dawrit(idaf,ioda,h,ll2,87,0)
    END IF

    l1 = num
    IF (uncon) l1 = numu
    l2 = l1*(l1+1)/2

    nqmmatm = 0
    IF (nfg/=0) THEN
        IF (qmmm) CALL getmmchg(nqmmatm,dummy)
    END IF

!    Zero out all arrays
     s(1:ll2) = 0.0
     t(1:ll2) = 0.0
     h(1:ll2) = 0.0
     IF (doesp1) esp1e(1:ll2) = 0.0
     IF (lzint) z(1:ll2) = 0.0

    IF (dccut) CALL daread(idaf,ioda,nlct,ll2,272,1)

    CALL kin_ovl_ints(s, t, nlct, nosp_basis, dccut, is_lut, tol)

    IF (.NOT.doesp1) CALL nuc_ints(h, nlct, nosp_basis, dccut, is_lut, tol)

!   `NCHMAT` is nonzero if there are external charges which
!   perturb the system, such as if CHARMM is in use. Note
!   that there is also a nuclear repulsion term which is not
!   included here, it is in the charmm interface code.

!   COSMO or CHARMM point charges
    IF (iseps) CALL int1_coul_ext_chg(h, nosp_basis, nps, &
                corzan, qscnet, tol, 1.0d-8)

    IF (doesp0) THEN
        fmoc(1:3,1:nat1e) => x(lfmoc:lfmoc+3*nat1e-1)
        CALL int1_coul_ext_chg(h, nosp_basis, nat1e, &
                fmoc, atchrg(1:nat1e), tol, 1.0d-8)
    END IF

    IF (nqmmatm/=0) CALL int1_coul_ext_chg(h, nosp_basis, nqmmatm, &
                xchm, ychm, zchm, qchm, tol, 1.0d-8)

    IF (doesp1) CALL hstesp(h, esp1e, nosp_basis, nqmmatm, tol)

    IF (lzint) CALL lzints(z, nlct, nosp_basis, dccut, tol)

!   Normalize 1-e integrals all at once
    CALL bas_norm_matrix(h, nosp_basis%bfnrm, l1)
    CALL bas_norm_matrix(s, nosp_basis%bfnrm, l1)
    CALL bas_norm_matrix(t, nosp_basis%bfnrm, l1)
    IF (lzint)  CALL bas_norm_matrix(z,     nosp_basis%bfnrm, l1)
    IF (doesp1) CALL bas_norm_matrix(esp1e, nosp_basis%bfnrm, l1)

!   Sum up partial contributions if parallel
    IF (goparr) THEN
       CALL ddi_gsumf(910,h,l2)
       CALL ddi_gsumf(911,s,l2)
       CALL ddi_gsumf(912,t,l2)
       IF (lzint)  CALL ddi_gsumf(913,z,l2)
       IF (doesp1) CALL ddi_gsumf(913,esp1e,l2)
    END IF

!   Form one electron Hamiltonian
!   Hcore = Vne + Te
    h(1:l2) = h(1:l2) + t(1:l2)

!   Save H, S, and T matrices on the DAF
    CALL dawrit(idaf,ioda,h,ll2,11,0)
    CALL dawrit(idaf,ioda,s,ll2,12,0)
    CALL dawrit(idaf,ioda,t,ll2,13,0)

!   Various solvent models unfortunately mess up record 11,
!   so keep a pristine copy.
!   It will lack any corrections for core potentials or
!   for scalar relativity, however!
    CALL dawrit(idaf,ioda,h,ll2,290,0)
    IF (lzint) CALL dawrit(idaf,ioda,z,ll2,379,0)
    IF (doesp1) CALL dawrit(idaf,ioda,esp1e,l2,311,0)
    IF (is_lut) CALL dawrit(idaf,ioda,h,ll2,396,0)

!   Optional debug printout
    IF (dbug) THEN
       WRITE(iw,*) 'OVERLAP MATRIX'
       CALL prtril(s,l1)
       WRITE(iw,*) 'BARE NUCLEUS HAMILTONIAN INTEGRALS (H=T+V)'
       CALL prtril(h,l1)
       WRITE(iw,*) 'KINETIC ENERGY INTEGRALS'
       CALL prtril(t,l1)
       IF (lzint) THEN
          WRITE(iw,*) 'Z-ANGULAR MOMENTUM INTEGRALS'
          CALL prtril(z,l1)
       END IF
       IF (doesp1) THEN
          WRITE(iw,*) 'FMO ESP INTEGRALS'
          CALL prtril(esp1e,l1)
       END IF
    END IF

 END SUBROUTINE

!-------------------------------------------------------------------------------

!> @brief Compute 1e ESP on grid in FMO/FMO-PCM jobs
!> @details MPI/OpenMP version of `PCMPOT` subroutine
!> @note mode=0 compute V and update the Fock matrix H (PCM)
!> @note mode=1 only compute V (PCM)
!> @note mode=2 only compute V (1e ESP)
!> @note `SOME` must be set identically on master and slaves (never set it to true on masters only).
!
!> @author Vladimir Mironov
!
!   REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!>
!> @param[in,out]   h       one-electron Hamiltonian matrix in packet format
!> @param[in,out]   vpcm    packed matrix of one-electron Coulomb integrals
!> @param[in]       qse     array of particle charges
!> @param[in]       xyzcts  array of particle coordinates in AoS or SoA format
!> @param[in]       l2      size of matrices `H` and `VPCM`
!> @param[in]       mode    determines AoS (`==2`) or SoA (`/=2`) layout of `XYZCTS`; also, if `mode==0` flush `VPCM` on disk and add it to `H`
!> @param[in]       some    flag for timing and output
 SUBROUTINE grid_fmo_ints(h,vpcm,npt,qse,xyzcts,l2,mode,some)

    REAL(KIND=fp), CONTIGUOUS, INTENT(INOUT) :: h(:), vpcm(:)
    REAL(KIND=fp), INTENT(IN) :: qse(*)
    REAL(KIND=fp), TARGET, CONTIGUOUS, INTENT(IN) :: xyzcts(:,:) !xyzcts(mxts,3)
    INTEGER, INTENT(IN) :: npt, l2, mode
    LOGICAL, INTENT(IN) :: some

    COMMON /DFTB  / dftbfl,scc,srscc,dftb3,dampxh,lrdftb
        LOGICAL :: dftbfl,scc,srscc,dftb3,dampxh,lrdftb

    COMMON /OUTPUT/ nprint,itol,icut,normf,normp,nopk
      INTEGER :: nprint,itol,icut,normf,normp,nopk

    COMMON /PAR   / me,master,nproc,ibtyp,iptim,goparr,dskwrk,maswrk
      INTEGER :: me,master,nproc,ibtyp,iptim
      LOGICAL :: goparr,dskwrk,maswrk

    COMMON /PCMDIM/ mxsp,mxts,mempcm1,mempcm2,nts
        INTEGER :: mxsp,mxts,mempcm1,mempcm2,nts

    COMMON /IOFILE/ ir,iw,ip,is,ipk,idaf,nav,ioda(950)
      INTEGER :: ir,iw,ip,is,ipk,idaf,nav,ioda

    REAL(KIND=fp) :: tol
    REAL(KIND=fp), POINTER :: xyz(:,:)
    INTEGER :: l1

!    CALL prepQuadGaussHermite

    IF (some) CALL timit(1)

!    WRITE(*,'("OMP PCMPOT(MODE=",I1,")")') mode

    l1 = floor(sqrt(2.0*l2))

    tol = RLN10*(itol-4)

    vpcm(1:l2) = 0.0

    IF (dftbfl) THEN
        CALL dftb_pcmpot(vpcm,qse,xyzcts,l2)
    ELSE IF (mode==2) THEN
        tol = RLN10*itol
        xyz(1:3,1:npt) => xyzcts !AoS layout
        CALL int1_coul_ext_chg(vpcm(1:l2), nosp_basis, npt, &
                               xyz, qse(1:npt), tol, 0.0d0)
    ELSE
        xyz(1:mxts,1:3) => xyzcts ! SoA layout
        CALL int1_coul_ext_chg(vpcm(1:l2), nosp_basis, npt, &
                               xyz(:,1), xyz(:,2), xyz(:,3), qse(1:npt), tol, 0.0d0)
    END IF

    IF (goparr) THEN
       !IF (nxt) CALL ddi_dlbreset
       CALL ddi_gsumf(914,vpcm,l2)
    END IF

    CALL bas_norm_matrix(vpcm, nosp_basis%bfnrm, l1)

    IF (mode==0) THEN
!nb   This write may cause trouble to FMO/F?!
      CALL dawrit(idaf,ioda,vpcm,l2,313,0)
!     Add VPCM to the core (aka 1e) Hamiltonian.
      CALL daread(idaf,ioda,h,l2,11,0)
      h(1:l2) = h(1:l2) + vpcm(1:l2)
      CALL dawrit(idaf,ioda,h,l2,11,0)
    END IF

    IF (some) THEN
      IF (mode/=2.AND.maswrk) WRITE(iw,*) 'Done VPCM.'
      CALL timit(1)
    endif

 END SUBROUTINE

!-------------------------------------------------------------------------------

!> @brief Compute nuclear attraciton 1e ESP in FMO jobs
!> @details This subroutine is a wrapper to two computational
!>  kernels: FASTVESP or general case
!> @param[in,out]   h       one-electron Hamiltonian matrix in packet format
!> @param[in,out]   esp1e   packed matrix of one-electron Coulomb integrals
!> @param[in]       basis   basis w/ SP-shells separated
!> @param[in]       nqmmatm number of QMMM atoms in job
!> @param[in]       tol     1-e exponential prefactor tolerance
 SUBROUTINE hstesp(h,esp1e,basis,nqmmatm,tol)

    REAL(KIND=fp), CONTIGUOUS, INTENT(INOUT) :: h(:), esp1e(:) ! All dims. LL2
    INTEGER, INTENT(IN) :: nqmmatm
    REAL(KIND=fp), INTENT(IN) :: tol
    TYPE(basis_set), INTENT(IN) :: basis ! basis without sp-shells
!    LOGICAL, INTENT(IN) :: fastvesp

    COMMON /EFMO  / efmoetot,efmoeserg,efmopolerg,efmodiserg, &
                    efmorepnrg,efmochtnrg,efmoepen,efmopcmg, &
                    iefmorun,imodefp,iefmort,iefmocfrg,iefmonfrg, &
                    iefmodim(2),imodefe,natefmo,imodefd,imodefer, &
                    imodefct,idimtyp,iefmo_agrad
      REAL(KIND=fp) :: efmoetot,efmoeserg,efmopolerg,efmodiserg, &
          efmorepnrg,efmochtnrg,efmoepen,efmopcmg
      INTEGER :: iefmorun,imodefp,iefmort,iefmocfrg,iefmonfrg, &
          iefmodim,imodefe,natefmo,imodefd,imodefer, &
          imodefct,idimtyp,iefmo_agrad

    COMMON /FMOINF/ nfg,nlayer,natfmo,nbdfg,naotyp,nbody
      INTEGER :: nfg,nlayer,natfmo,nbdfg,naotyp,nbody

    COMMON /FMOOPT/ espsca(9),respap(2),resppc(2),resdim,restri(4), &
                    rcorsd,respct,convfg,cnvdmp,coroff,rflmo(4), &
                    orshft,orshft2,cnvafo,ascreen(4),ixesp,mxitfg, &
                    nguess,nbsse,modorb,modpar,irststp,irstlay,nprfmo, &
                    nfmopal,modprp,maxl1c,ipieda,modgrd,modesp,ivmul, &
                    modlmo,nopden,mofock,modfd,modfmm,ncentm,ndualb
      REAL(KIND=fp) :: espsca,respap,resppc,resdim,restri, &
          rcorsd,respct,convfg,cnvdmp,coroff,rflmo, &
          orshft,orshft2,cnvafo,ascreen
      INTEGER :: ixesp,mxitfg,nguess,nbsse,modorb,modpar,&
          irststp,irstlay,nprfmo,nfmopal,modprp,maxl1c,ipieda,&
          modgrd,modesp,ivmul,modlmo,nopden,mofock,modfd,&
          modfmm,ncentm,ndualb

    COMMON /FMORUN/ espscf,e0scf(2),emp2s,idafmo,icurfg,jcurfg,kcurfg, &
                    icurlay,icurunt,nat1e,ncursh,ngau,icurpop,ifmostp, &
                    moncor,needr,modrst,norbproj,nunesp,iskipesp, &
                    iesdppc,idoprop,mp2run,icurit,idmfmo,iddfmo, &
                    iddcur,nddleft,ivmfmo,nzmtfmo,ifmobas,itmfmo(2)
      REAL(KIND=fp) :: espscf,e0scf,emp2s
      INTEGER :: idafmo,icurfg,jcurfg,kcurfg, &
          icurlay,icurunt,nat1e,ncursh,ngau,icurpop,ifmostp, &
          moncor,needr,modrst,norbproj,nunesp,iskipesp, &
          iesdppc,idoprop,mp2run,icurit,idmfmo,iddfmo, &
          iddcur,nddleft,ivmfmo,nzmtfmo,ifmobas,itmfmo

    COMMON /RELWFN/ rmethod,qrqmt,clig,clig2,qrtol,tau, &
                    iqrord,modqr,nesoc,nratom, &
                    numu,nqmtr,nqrdaf,morda,ndarelb
      REAL(KIND=fp) :: rmethod,qrqmt,clig,clig2,qrtol,tau
      INTEGER :: iqrord,modqr,nesoc,nratom,numu,nqmtr,nqrdaf,morda,ndarelb

    LOGICAL :: dampall, dodamp, fastvesp
    INTEGER :: needfv, lfvesp

!   norm = normf /= 1 .OR. normp /= 1

    dampall = iand(modesp,64)/=0
    dodamp = ascreen(1)/=0.AND.(nbdfg/=0.OR.dampall).AND. &
             ifmostp/=6.AND.ifmostp/=7.AND.iefmorun==0

!   SETATZ allocates memory. One should be careful in returning it.
    CALL setatz(1,nat1e,fastvesp,lfvesp,needfv)

    IF (fastvesp) THEN
        CALL fmo_fastvesp(h,esp1e,basis,&
                nat1e,ascreen(1),ascreen(2),nqmmatm,tol,lfvesp, &
                dodamp, dampall)
    ELSE
        CALL fmo_esp(h,esp1e,basis,&
                nat1e,ascreen(1),ascreen(2),tol, &
                dodamp, dampall, rmethod==RLUT)
    END IF

    !IF (nfg/=0.AND.fastvesp) CALL retfm(needfv)
    IF (fastvesp) CALL retfm(needfv)

 END SUBROUTINE

!-------------------------------------------------------------------------------

!> @brief Compute nuclear attraciton 1e ESP in FMO jobs, FASTVESP case
!> @details The integrals are comptued over all shell pairs of current
!>  fragment and all FMO atoms. The coordinates and scaling coefficients
!>  for 1e potential has to be precomputed and stored in global
!>  memory (`FMCOM`) at address `LFVESP`
!> @param[in,out]   h       one-electron Hamiltonian matrix in packet format
!> @param[in,out]   esp1e   packed matrix of one-electron Coulomb integrals
!> @param[in]       basis   basis w/ SP-shells separated
!> @param[in]       nat1e   number of FMO atoms for which ESP is computed
!> @param[in]       alpha   dumping exponent
!> @param[in]       beta    exponential dumping factor
!> @param[in]       nqmmatm number of QMMM atoms in job
!> @param[in]       tol     1-e exponential prefactor tolerance
!> @param[in]       lfvesp  pointer to coordinates, nuclear charges and scaling coefficients in global (`FMCOM`) memory
!> @param[in]       dodamp  dumping flag, dump only `FRACV/={0,1}` integrals
!> @param[in]       dampall dumping flag, dump all Coulomb integrals
 SUBROUTINE fmo_fastvesp(h, esp1e, basis,&
                 nat1e, alpha, beta, nqmmatm, tol, lfvesp,&
                 dodamp, dampall)

    REAL(KIND=fp), CONTIGUOUS, INTENT(INOUT) :: h(:), esp1e(:) ! All dims. LL2
    INTEGER, INTENT(IN) :: nat1e, nqmmatm, lfvesp
    REAL(KIND=fp), INTENT(IN) :: alpha, beta, tol
    TYPE(basis_set), INTENT(IN) :: basis ! basis without sp-shells
    LOGICAL, INTENT(IN) :: dampall, dodamp

    COMMON /PAR   / me,master,nproc,ibtyp,iptim,goparr,dskwrk,maswrk
      INTEGER :: me,master,nproc,ibtyp,iptim
      LOGICAL :: goparr,dskwrk,maswrk

    COMMON /FMCOM / x(1)
      REAL(KIND=fp) :: x

    REAL(KIND=fp) :: fracesp, fracv
    REAL(KIND=fp), DIMENSION(BLOCKSIZE) :: vblk, zblk
!dir$ attributes align : 64 :: vblk, zblk
    INTEGER :: icind, ii, jj
    TYPE(shell_t) :: shi, shj
    TYPE(shpair_t) :: cntp

!$omp parallel &
!$omp   private( &
!$omp       ii, jj, &
!$omp       icind, fracv, fracesp, &
!$omp       vblk, zblk, &
!$omp       shi, shj, cntp &
!$omp   )

    CALL cntp%alloc

    fracesp = 0.0
    fracv = 0.0

!   I shell
    DO ii = basis%nshell, 1, -1
!       Go parallel (SLB)
        IF (goparr) THEN
            IF (mod(ii,nproc)/=me) CYCLE
        END IF

        CALL shi%fetch_by_id(ii)

!       J shell
!$omp do schedule(dynamic)
        DO jj = 1, ii

            CALL shj%fetch_by_id(jj)

            CALL cntp%shell_pair(shi, shj, tol)
            IF (cntp%numpairs==0) CYCLE

            vblk = 0.0
            zblk = 0.0

            CALL int1_fmoesp_ffcxyz(cntp, x(lfvesp:lfvesp+6*nat1e), &
                             nat1e, alpha, beta, &
                             dodamp, dampall, vblk, zblk)

            IF (nqmmatm>0) THEN
                icind = lfvesp + nat1e*6
                CALL int1_coul(cntp, x(icind:icind+nqmmatm*4), nqmmatm, 0.0_fp, vblk)
            END IF

            CALL update_triang_matrix(shi, shj, vblk, h)
            CALL update_triang_matrix(shi, shj, zblk, esp1e)

        END DO
!$omp end do nowait
    END DO
!$omp end parallel
!   End of shell loops

 END SUBROUTINE

!-------------------------------------------------------------------------------

!> @brief Compute nuclear attraciton 1e ESP in FMO jobs, general case
!> @details The integrals are comptued over all shell pairs of current
!>  fragment and all FMO atoms. The coordinates and scaling coefficients
!>  for 1e potential are computed at each step with `FMATFRG` subroutine
!> @param[in,out]   h       one-electron Hamiltonian matrix in packet format
!> @param[in,out]   esp1e   packed matrix of one-electron Coulomb integrals
!> @param[in]       basis   basis w/ SP-shells separated
!> @param[in]       nat1e   number of FMO atoms for which ESP is computed
!> @param[in]       alpha   dumping exponent
!> @param[in]       beta    exponential dumping factor
!> @param[in]       tol     1-e exponential prefactor tolerance
!> @param[in]       dodamp  dumping flag, dump only `FRACV/={0,1}` integrals
!> @param[in]       dampall dumping flag, dump all Coulomb integrals
!> @param[in]       is_lut  .TRUE. if `RMETHOD=='LUT'`
 SUBROUTINE fmo_esp(h, esp1e, basis,&
                 nat1e, alpha, beta, tol, &
                 dodamp, dampall, is_lut)

    REAL(KIND=fp), CONTIGUOUS, INTENT(INOUT) :: h(:), esp1e(:) ! All dims. LL2
    INTEGER, INTENT(IN) :: nat1e
    REAL(KIND=fp), INTENT(IN) :: alpha, beta, tol
    TYPE(basis_set), INTENT(IN) :: basis ! basis without sp-shells
    LOGICAL, INTENT(IN) :: dampall, dodamp, is_lut

    COMMON /PAR   / me,master,nproc,ibtyp,iptim,goparr,dskwrk,maswrk
      INTEGER :: me,master,nproc,ibtyp,iptim
      LOGICAL :: goparr,dskwrk,maswrk

    REAL(KIND=fp), DIMENSION(BLOCKSIZE) :: vblk, zblk
!dir$ attributes align : 64 :: vblk, zblk
    INTEGER :: ic, ii, jj
    TYPE(shell_t) :: shi, shj
    TYPE(shpair_t) :: cntp

!$omp parallel &
!$omp   private( &
!$omp       ii, jj, ic, &
!$omp       vblk, zblk, &
!$omp       shi, shj, cntp &
!$omp   )

    CALL cntp%alloc

!   I shell
    DO ii = basis%nshell, 1, -1
!       Go parallel (SLB)
        IF (goparr) THEN
            IF (mod(ii,nproc)/=me) CYCLE
        END IF

        CALL shi%fetch_by_id(ii)

!       J shell
!$omp do schedule(dynamic)
        DO jj = 1, ii

            CALL shj%fetch_by_id(jj)

            CALL cntp%shell_pair(shi, shj, tol)
            IF (cntp%numpairs==0) CYCLE

            vblk = 0.0
            zblk = 0.0

            CALL int1_fmoesp_fcall(cntp, shi%atid, shj%atid, nat1e, alpha, beta, &
                             dodamp, dampall, is_lut, vblk, zblk)

            CALL update_triang_matrix(shi, shj, vblk, h)
            CALL update_triang_matrix(shi, shj, zblk, esp1e)

        END DO
!$omp end do nowait
    END DO
!$omp end parallel
!   End of shell loops

 END SUBROUTINE

!-------------------------------------------------------------------------------

!> @brief Compute overlap and kinetic integrals
!
!> @details Overlap and electron kinetic energy integrals
!>  are computed using Gauss-Hermite quadrature formula
!>  Kinetic energy integrals are actually overlap integrals with +2, -2 angular
!>  momentum shifts
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!>
!> @param[in,out]   s       packed matrix of overlap integrals
!> @param[in,out]   t       packed matrix of kinetic energy integrals
!> @param[in]       nlct    an array for Divide-and-Conquer (DC) method screening
!> @param[in]       basis   basis w/ SP-shells separated
!> @param[in]       dccut   flag specifying DC job
!> @param[in]       is_lut  .TRUE. if `RMETHOD=='LUT'`
!> @param[in]       tol     1-e exponential prefactor tolerance
 SUBROUTINE kin_ovl_ints(s, t, nlct, basis, dccut, is_lut, tol)

    REAL(KIND=fp), CONTIGUOUS,  INTENT(INOUT)  :: s(:), t(:)
    TYPE(basis_set), INTENT(IN)     :: basis ! basis without sp-shells
    INTEGER, CONTIGUOUS,        INTENT(IN)     :: nlct(:)
    LOGICAL,         INTENT(IN)     :: dccut, is_lut
    REAL(KIND=fp),   INTENT(IN)     :: tol

    COMMON /PAR   / me,master,nproc,ibtyp,iptim,goparr,dskwrk,maswrk
      INTEGER :: me,master,nproc,ibtyp,iptim
      LOGICAL :: goparr,dskwrk,maswrk

    INTEGER :: &
        ii, jj, ipos

    LOGICAL :: dokinetic

    REAL(KIND=fp), DIMENSION(BLOCKSIZE) :: tblk, sblk
!dir$ attributes align : 64 :: tblk, sblk

    LOGICAL :: dcfskp

    TYPE(shell_t) :: shi, shj
    TYPE(shpair_t) :: cntp

    dcfskp = .FALSE.

!$omp parallel &
!$omp   private( &
!$omp       ii, jj, &
!$omp       sblk, tblk, &
!$omp       shi, shj, cntp, dokinetic, ipos &
!$omp   ) &
!$omp   firstprivate( &
!$omp       dcfskp &
!$omp   )

    CALL cntp%alloc

!   I shell
    DO ii = basis%nshell, 1, -1
!       Go parallel (SLB)
        IF (goparr.AND.(mod(ii,nproc)/=me)) CYCLE

        CALL shi%fetch_by_id(ii)
        ipos = (shi%locao*(shi%locao-1))/2

!       J shell
!$omp do schedule(dynamic)
        DO jj = 1, ii

            CALL shj%fetch_by_id(jj)
            IF (dccut) dcfskp = (nlct(ipos+shj%locao)==0)

!           Skip calculation of elements
!           which are not included in DC-Fock (T and V)
            dokinetic = (.NOT.dcfskp).AND. &
                        (.NOT.(is_lut .AND. shi%atid==shj%atid))

            CALL cntp%shell_pair(shi, shj, tol)
            IF (cntp%numpairs==0) CYCLE

            sblk = 0.0
            tblk = 0.0

            CALL int1_kin_ovl(cntp, dokinetic, sblk, tblk)

            CALL update_triang_matrix(shi, shj, sblk, s)
            IF (.NOT.dcfskp) CALL update_triang_matrix(shi, shj, tblk, t)

        END DO
!$omp end do nowait
    END DO
!$omp end parallel
!   End of shell loops
 END SUBROUTINE

!-------------------------------------------------------------------------------

!> @brief Compute \f$ L_z \f$ integrals
!
!> @details \f$ L_z \f$ are actually overlap integrals with +1, -1 angular
!>  momentum shifts
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!>
!> @param[in,out]   z       packed matrix of Lz integrals
!> @param[in]       nlct    an array for Divide-and-Conquer (DC) method screening
!> @param[in]       basis   basis w/ SP-shells separated
!> @param[in]       dccut   flag specifying DC job
!> @param[in]       tol     1-e exponential prefactor tolerance
 SUBROUTINE lzints(z, nlct, basis, dccut, tol)

    TYPE(basis_set), INTENT(IN)  :: basis ! basis without sp-shells
    REAL(KIND=fp), CONTIGUOUS,  INTENT(OUT) :: z(:)
    LOGICAL,         INTENT(IN)  :: dccut
    REAL(KIND=fp),   INTENT(IN)  :: tol
    INTEGER,  CONTIGUOUS,       INTENT(IN)  :: nlct(:)

    COMMON /PAR   / me,master,nproc,ibtyp,iptim,goparr,dskwrk,maswrk
      INTEGER :: me,master,nproc,ibtyp,iptim
      LOGICAL :: goparr,dskwrk,maswrk

    INTEGER :: &
        ii, jj, ipos

    REAL(KIND=fp), DIMENSION(BLOCKSIZE) :: zblk
!dir$ attributes align : 64 :: zblk

    TYPE(shell_t) :: shi, shj
    TYPE(shpair_t)  :: cntp

    CALL cntp%alloc

!   I shell
    DO ii = 1, basis%nshell
        CALL shi%fetch_by_id(ii)
        ipos = (shi%locao*(shi%locao-1))/2

!       J shell
        DO jj = 1, ii

!           Go parallel (SLB)
            IF (goparr) THEN
                IF (mod(jj,nproc)/=me) CYCLE
            END IF

            CALL shj%fetch_by_id(jj)
            IF (dccut .AND. nlct(ipos+shj%locao)==0) CYCLE

            CALL cntp%shell_pair(shi, shj, tol)
            IF (cntp%numpairs==0) CYCLE

            zblk = 0.0
            CALL int1_lz(cntp, zblk)
            CALL update_triang_matrix(shi, shj, zblk, z)

        END DO
    END DO
!   End of shell loops
 END SUBROUTINE

!-------------------------------------------------------------------------------

!> @brief Compute nuclear attraction integrals
!
!> @details Nuclear attaction integrals are computed using Gauss-Rys quadrature
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!>
!> @param[in,out]   h       core Hamiltonian matrix
!> @param[in]       nlct    an array for Divide-and-Conquer (DC) method screening
!> @param[in]       basis   basis w/ SP-shells separated
!> @param[in]       dccut   flag specifying DC job
!> @param[in]       is_lut  .TRUE. if `RMETHOD=='LUT'`
!> @param[in]       tol     1-e exponential prefactor tolerance
 SUBROUTINE nuc_ints(h, nlct, basis, dccut, is_lut, tol)

    REAL(KIND=fp), CONTIGUOUS,  INTENT(INOUT)  :: h(:)
    TYPE(basis_set), INTENT(IN)     :: basis ! basis without sp-shells
    INTEGER,  CONTIGUOUS,       INTENT(IN)     :: nlct(:)
    LOGICAL,         INTENT(IN)     :: dccut, is_lut
    REAL(KIND=fp),   INTENT(IN)     :: tol

    COMMON /INFOA / nat,ich,mul,num,nqmt,ne,na,nb, &
                    zan(mxatm),c(3,mxatm),ian(mxatm)
      INTEGER :: nat,ich,mul,num,nqmt,ne,na,nb,ian
      REAL(KIND=fp) :: zan,c

    COMMON /PAR   / me,master,nproc,ibtyp,iptim,goparr,dskwrk,maswrk
      INTEGER :: me,master,nproc,ibtyp,iptim
      LOGICAL :: goparr,dskwrk,maswrk

    INTEGER :: &
        ii, jj, ipos, iskp

    REAL(KIND=fp), DIMENSION(BLOCKSIZE) :: vblk
!dir$ attributes align : 64 :: vblk

    TYPE(shell_t) :: shi, shj
    TYPE(shpair_t) :: cntp

!$omp parallel &
!$omp   private( &
!$omp       ii, jj, iskp, &
!$omp       vblk, &
!$omp       shi, shj, cntp, ipos &
!$omp   )

    CALL cntp%alloc
    iskp = -1

!   I shell
    DO ii = basis%nshell, 1, -1
!       Go parallel (SLB)
        IF (goparr.AND.(mod(ii,nproc)/=me)) CYCLE

        CALL shi%fetch_by_id(ii)
        ipos = (shi%locao*(shi%locao-1))/2

!       J shell
!$omp do schedule(dynamic)
        DO jj = 1, ii

            CALL shj%fetch_by_id(jj)
            IF (dccut .AND. nlct(ipos+shj%locao)==0) CYCLE

            CALL cntp%shell_pair(shi, shj, tol)
            IF (cntp%numpairs==0) CYCLE

            IF (is_lut .AND. shi%atid==shj%atid) iskp = shi%atid

            vblk = 0.0

            CALL int1_coul(cntp, c, zan, nat, 0.0_fp, vblk, iskp)

            CALL update_triang_matrix(shi, shj, vblk, h)


        END DO
!$omp end do nowait
    END DO
!$omp end parallel
!   End of shell loops

 END SUBROUTINE

!-------------------------------------------------------------------------------

!> @brief General way to compute integrals of charge interaction, charge
!>  data are stored in array of structures like xyz(3,:), charge(:)
!
!> @details Electron-charge interaction integrals are computed using
!>  Gauss-Rys quadrature
!> @note This case is a variation of nuclear attraction case. They differ in
!>  data representation and in screening logic
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!>
!> @param[in,out]   h       packed matrix of one-electon Coulomb integrals
!> @param[in]       basis   basis w/ SP-shells separated
!> @param[in]       nat     number of atoms
!> @param[in]       xyz     array of particle coordinates in AoS format
!> @param[in]       chg     array of particle charges
!> @param[in]       tol     1-e exponential prefactor tolerance
!> @param[in]       chgtol  tolerance for particle charge
 SUBROUTINE chg_ints_aos(h, basis, nat, xyz, chg, tol, chgtol)

    REAL(KIND=fp), CONTIGUOUS,  INTENT(INOUT)  :: h(:)
    TYPE(basis_set), INTENT(IN)     :: basis ! basis without sp-shells
    INTEGER,         INTENT(IN)     :: nat
    REAL(KIND=fp), CONTIGUOUS,  INTENT(IN)     :: xyz(:,:), chg(:)
    REAL(KIND=fp),   INTENT(IN)     :: tol, chgtol

    COMMON /PAR   / me,master,nproc,ibtyp,iptim,goparr,dskwrk,maswrk
      INTEGER :: me,master,nproc,ibtyp,iptim
      LOGICAL :: goparr,dskwrk,maswrk

    INTEGER :: &
        ii, jj

    REAL(KIND=fp), DIMENSION(BLOCKSIZE) :: vblk
!dir$ attributes align : 64 :: vblk

    TYPE(shell_t) :: shi, shj
    TYPE(shpair_t) :: cntp

!$omp parallel &
!$omp   private( &
!$omp       ii, jj, &
!$omp       vblk, &
!$omp       shi, shj, cntp &
!$omp   )

    CALL cntp%alloc

!   I shell
    DO ii = basis%nshell, 1, -1
!       Go parallel (SLB)
        IF (goparr.AND.(mod(ii,nproc)/=me)) CYCLE

        CALL shi%fetch_by_id(ii)

!       J shell
!$omp do schedule(dynamic)
        DO jj = 1, ii

            CALL shj%fetch_by_id(jj)

            CALL cntp%shell_pair(shi, shj, tol)
            IF (cntp%numpairs==0) CYCLE

            vblk = 0.0

            CALL int1_coul(cntp, xyz, chg, nat, chgtol, vblk, -1)

            CALL update_triang_matrix(shi, shj, vblk, h)

        END DO
!$omp end do nowait
    END DO
!$omp end parallel
!   End of shell loops

 END SUBROUTINE

!-------------------------------------------------------------------------------

!> @brief General way to compute integrals of charge interaction, charge
!>  data are stored in a structure of arrays x(:),y(:),z(:),charge(:)
!
!> @details Electron-charge interaction integrals are computed using
!>  Gauss-Rys quadrature
!> @note This case is a variation of nuclear attraction case. They differ in
!>  data representation and in screening logic
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!>
!> @param[in,out]   h       packed matrix of one-electon Coulomb integrals
!> @param[in]       basis   basis w/ SP-shells separated
!> @param[in]       nat     number of atoms
!> @param[in]       x       array of X particle coordinates
!> @param[in]       y       array of Y particle coordinates
!> @param[in]       z       array of Z particle coordinates
!> @param[in]       chg     array of particle charges
!> @param[in]       tol     1-e exponential prefactor tolerance
!> @param[in]       chgtol  tolerance for particle charge
 SUBROUTINE chg_ints_soa(h, basis, nat, x, y, z, chg, tol, chgtol)

    REAL(KIND=fp), CONTIGUOUS,  INTENT(INOUT)  :: h(:)
    TYPE(basis_set), INTENT(IN)     :: basis ! basis without sp-shells
    INTEGER,         INTENT(IN)     :: nat
    REAL(KIND=fp), CONTIGUOUS,  INTENT(IN)     :: x(:), y(:), z(:), chg(:)
    REAL(KIND=fp),   INTENT(IN)     :: tol, chgtol

    COMMON /PAR   / me,master,nproc,ibtyp,iptim,goparr,dskwrk,maswrk
      INTEGER :: me,master,nproc,ibtyp,iptim
      LOGICAL :: goparr,dskwrk,maswrk

    INTEGER :: &
        ii, jj

    REAL(KIND=fp), DIMENSION(BLOCKSIZE) :: vblk
!dir$ attributes align : 64 :: vblk

    TYPE(shell_t) :: shi, shj
    TYPE(shpair_t) :: cntp

!$omp parallel &
!$omp   private( &
!$omp       ii, jj, &
!$omp       vblk, &
!$omp       shi, shj, cntp &
!$omp   )
    CALL cntp%alloc

!   I shell
    DO ii = basis%nshell, 1, -1
!       Go parallel (SLB)
        IF (goparr.AND.(mod(ii,nproc)/=me)) CYCLE

        CALL shi%fetch_by_id(ii)

!       J shell
!$omp do schedule(dynamic)
        DO jj = 1, ii

            CALL shj%fetch_by_id(jj)

            CALL cntp%shell_pair(shi, shj, tol)
            IF (cntp%numpairs==0) CYCLE

            vblk = 0.0

            CALL int1_coul(cntp, x, y, z, chg, nat, chgtol, vblk)

            CALL update_triang_matrix(shi, shj, vblk, h)

        END DO
!$omp end do nowait
    END DO
!$omp end parallel
!   End of shell loops

 END SUBROUTINE

!-------------------------------------------------------------------------------

!> @brief Compute electrostatic potential of a frament
!>  on atoms of other fragments
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Oct, 2018_ Initial release
!>
!> @param[in]       basis   basis w/ SP-shells separated
!> @param[in]       aconst  scaling value for `VAL`
!> @param[in,out]   val     array of electron-charge repulsion energies
!> @param[in,out]   da      density matrix in packed format, remains unchanged on return
!> @param[in]       np      number of particles
!> @param[in]       l2      dimension of `DA` array
!> @param[in]       norm    flag for integrals normalization
!> @param[in]       tol     1-e exponential prefactor tolerance
 SUBROUTINE omp_fmo_epoten(basis,aconst,val,da,np,l2,norm,tol)

    TYPE(basis_set), INTENT(IN) :: basis ! basis without sp-shells
    INTEGER, INTENT(IN) :: np, l2
    REAL(KIND=fp), INTENT(IN) :: aconst, tol
    LOGICAL, INTENT(IN) :: norm
    REAL(KIND=fp), CONTIGUOUS, INTENT(INOUT) :: da(:)
    REAL(KIND=fp), CONTIGUOUS, INTENT(OUT) :: val(:,:)

    COMMON /FMOINF/ nfg,nlayer,natfmo,nbdfg,naotyp,nbody
      INTEGER :: nfg,nlayer,natfmo,nbdfg,naotyp,nbody

    COMMON /FMOPNT/ lichfg,lmulfg,lidmrec,lfrgnam,llayfrg,lindat,lncbs,lfmozan,lfmoc,lfmomas,lizbas,liaglob,libdgh, &
          liabdfg,ljabdfg,lncao,lidxcao,liaprjo,ljaprjo,lcoreao,locccor,lshiftb,liodfmo,lfmoda,lfmodb, &
          lfmoespa,lfmoespb,llocfmo,lscffrg,lfmoscf,lrij,lpopmul,lpopmat,lialoc,lindbd,liatfrg,lindfrg, &
          lindgfrg,lnatfrg,lnat0frg,lianfrg,lzanfrg,lcfrg,llibish,llibnsh,llibng,lindatg,lfmobuf(3),lfmode, &
          lnumfrg,lloctat,liaoglob,lloadm,lfmoge,ldgrid,liodcfmo,ljob2grp,lfmopg,lemocdr,luntxyz,luntrot, &
          lstonep,lmapsu,lfrgmul,lclmo,lialmo,lindlmo,latclmo,llmobdf,lfgflmo,lnfglmo,llfglmo,lpfglmo, &
          lpopdmat,lidmpnt,liddpnt,livmpnt,liactfg,lcrfrg,lzlmfrgv,lylmfrgv,lndtfrg,lf_mm,lg_mm,lmaxl30,libuffg
      INTEGER :: lichfg,lmulfg,lidmrec,lfrgnam,llayfrg,lindat,lncbs,lfmozan,lfmoc,lfmomas,lizbas,liaglob,libdgh, &
          liabdfg,ljabdfg,lncao,lidxcao,liaprjo,ljaprjo,lcoreao,locccor,lshiftb,liodfmo,lfmoda,lfmodb, &
          lfmoespa,lfmoespb,llocfmo,lscffrg,lfmoscf,lrij,lpopmul,lpopmat,lialoc,lindbd,liatfrg,lindfrg, &
          lindgfrg,lnatfrg,lnat0frg,lianfrg,lzanfrg,lcfrg,llibish,llibnsh,llibng,lindatg,lfmobuf,lfmode, &
          lnumfrg,lloctat,liaoglob,lloadm,lfmoge,ldgrid,liodcfmo,ljob2grp,lfmopg,lemocdr,luntxyz,luntrot, &
          lstonep,lmapsu,lfrgmul,lclmo,lialmo,lindlmo,latclmo,llmobdf,lfgflmo,lnfglmo,llfglmo,lpfglmo, &
          lpopdmat,lidmpnt,liddpnt,livmpnt,liactfg,lcrfrg,lzlmfrgv,lylmfrgv,lndtfrg,lf_mm,lg_mm,lmaxl30,libuffg

    COMMON /PAR   / me,master,nproc,ibtyp,iptim,goparr,dskwrk,maswrk
      INTEGER :: me,master,nproc,ibtyp,iptim
      LOGICAL :: goparr,dskwrk,maswrk

    COMMON /FMCOM / x(1)
      REAL(KIND=fp) :: x

    REAL(KIND=fp) :: fracesp, fracv, znuc, cx, cy, cz, vsum
    INTEGER :: l1
    INTEGER :: ic, ii, jj, kfg, kat1, kat2, kfg1, kfg2
    TYPE(shell_t) :: shi, shj
    TYPE(shpair_t) :: cntp

    REAL(KIND=fp) :: dij(225) ! increase to 784 for I functions

    IF (norm) THEN
        l1 = floor(sqrt(2.0*l2))
        CALL bas_norm_matrix(da, nosp_basis%bfnrm, l1)
    END IF

!$omp parallel &
!$omp   private( &
!$omp       ii, jj, ic, &
!$omp       znuc, cx, cy, cz, fracv, fracesp, &
!$omp       dij, &
!$omp       vsum, &
!$omp       kat1, kat2, &
!$omp       kfg1, kfg2, &
!$omp       shi, shj, cntp &
!$omp   ) &
!$omp   reduction(+:val)

    CALL cntp%alloc

    fracesp = 0.0
    fracv = 0.0
    kat1 = 0
    kat2 = 0
    kfg1 = 0
    kfg2 = 0

!   I shell
    !DO ii = basis%nshell, 1, -1
    DO ii = 1, basis%nshell
!       Go parallel (SLB)
        IF (goparr) THEN
            IF (mod(ii,nproc)/=me) CYCLE
        END IF

        CALL shi%fetch_by_id(ii)

!       J shell
!$omp do schedule(dynamic)
        DO jj = 1, ii

            CALL shj%fetch_by_id(jj)

            CALL cntp%shell_pair(shi, shj, tol)
            IF (cntp%numpairs==0) CYCLE

            CALL density_ordered(shi, shj, dij, da)

            DO ic = 1, np

                CALL fmoatfrg2(ic,x(lindat),x(lindatg),x(liaglob), &
                               x(lialoc),x(liabdfg),x(ljabdfg), &
                               x(lindbd),x(lfmozan),x(lfmoc), &
                               natfmo+nbdfg,x(luntxyz),x(lpopmat), &
                               shi%atid,shj%atid,fracv,fracesp, &
                               kfg,znuc,cx,cy,cz,kfg1,kat1,kfg2,kat2)

                IF (kfg1==0 .AND. kfg2==0) CYCLE

                CALL int1_epoten(cntp, cx, cy, cz, dij, vsum)

                IF (kfg1/=0) THEN
                    !!$omp atomic
                    val(kat1,kfg1) = val(kat1,kfg1) + vsum*aconst
                END IF
                IF (kfg2/=0) THEN
                    !!$omp atomic
                    val(kat2,kfg2) = val(kat2,kfg2) + vsum*aconst
                END IF
            END DO

        END DO
!$omp end do nowait
    END DO
!$omp end parallel
!   End of shell loops
    IF (norm) THEN
        CALL bas_denorm_matrix(da, nosp_basis%bfnrm, l1)
    END IF

 END SUBROUTINE

!-------------------------------------------------------------------------------

!> @brief Compute electrostatic potential on a grid
!>  on atoms of other fragments
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Oct, 2018_ Initial release
!>
!> @param[in]       basis   basis w/ SP-shells separated
!> @param[in]       aconst  scaling value for `VAL`
!> @param[in]       xcoord  `X` coordinates of grid points
!> @param[in]       ycoord  `Y` coordinates of grid points
!> @param[in]       zcoord  `Z` coordinates of grid points
!> @param[in,out]   val     array of electron-charge repulsion energies
!> @param[in,out]   da      density matrix in packed format, remains unchanged on return
!> @param[in]       np      number of particles
!> @param[in]       l2      dimension of `DA` array
!> @param[in]       norm    flag for integrals normalization
!> @param[in]       tol     1-e exponential prefactor tolerance
 SUBROUTINE omp_epoten(basis,aconst,xcoord,ycoord,zcoord,val,da,np,l2,norm,tol)

    TYPE(basis_set), INTENT(IN) :: basis ! basis without sp-shells
    INTEGER, INTENT(IN) :: np, l2
    REAL(KIND=fp), INTENT(IN) :: aconst, tol
    REAL(KIND=fp), CONTIGUOUS, INTENT(IN) :: xcoord(:), ycoord(:), zcoord(:)
    LOGICAL, INTENT(IN) :: norm
    REAL(KIND=fp), CONTIGUOUS, INTENT(INOUT) :: da(:)
    REAL(KIND=fp), CONTIGUOUS, INTENT(OUT) :: val(:)

    COMMON /FMOINF/ nfg,nlayer,natfmo,nbdfg,naotyp,nbody
      INTEGER :: nfg,nlayer,natfmo,nbdfg,naotyp,nbody

    COMMON /PAR   / me,master,nproc,ibtyp,iptim,goparr,dskwrk,maswrk
      INTEGER :: me,master,nproc,ibtyp,iptim
      LOGICAL :: goparr,dskwrk,maswrk

    COMMON /FMCOM / x(1)
      REAL(KIND=fp) :: x

    REAL(KIND=fp) :: fracesp, fracv, znuc, cx, cy, cz, vsum
    INTEGER :: l1
    INTEGER :: ic, ii, jj, kat1, kat2, kfg1, kfg2
    TYPE(shell_t) :: shi, shj
    TYPE(shpair_t) :: cntp
    INTEGER :: icount

    REAL(KIND=fp) :: dij(225) ! increase to 784 for I functions

    IF (norm) THEN
        l1 = floor(sqrt(2.0*l2))
        CALL bas_norm_matrix(da, nosp_basis%bfnrm, l1)
    END IF

!$omp parallel &
!$omp   private( &
!$omp       ii, jj, ic, &
!$omp       znuc, cx, cy, cz, fracv, fracesp, &
!$omp       dij, &
!$omp       vsum, &
!$omp       icount, &
!$omp       kat1, kat2, &
!$omp       kfg1, kfg2, &
!$omp       shi, shj, cntp &
!$omp   ) &
!$omp   reduction(+:val)

    icount = 0

    CALL cntp%alloc

!   I shell
    !DO ii = basis%nshell, 1, -1
    DO ii = 1, basis%nshell
        CALL shi%fetch_by_id(ii)

!       J shell
        DO jj = 1, ii
            icount = icount + 1
!           Go parallel (SLB)
            IF (goparr) THEN
                IF (mod(icount,nproc)/=me) CYCLE
            END IF

            CALL shj%fetch_by_id(jj)

            CALL cntp%shell_pair(shi, shj, tol)
            IF (cntp%numpairs==0) CYCLE

            CALL density_ordered(shi, shj, dij, da)

!$omp do schedule(dynamic,8)
            DO ic = 1, np

                cx = xcoord(ic)
                cy = ycoord(ic)
                cz = zcoord(ic)

                CALL int1_epoten(cntp, cx, cy, cz, dij, vsum)

!                !$omp atomic
                val(ic) = val(ic) + vsum*aconst
            END DO
!$omp end do nowait

        END DO
    END DO
!$omp end parallel
!   End of shell loops
    IF (norm) THEN
        CALL bas_denorm_matrix(da, nosp_basis%bfnrm, l1)
    END IF

 END SUBROUTINE

!--------------------------------------------------------------------------------
!       ONE-ELECTRON INTEGRALS CALCULATION (CONTRACTED SHELLS)
!--------------------------------------------------------------------------------

!> @brief Compute contracted block of kinetic energy and overlap 1e integrals
!> @param[in]       cntp        shell pair data
!> @param[in]       dokinetic   if `.FALSE.` compute only overlap integrals
!> @param[out]      sblk        block of overlap integrals
!> @param[out]      tblk        block of kinetic energy integrals
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!
 SUBROUTINE int1_kin_ovl(cntp, dokinetic, sblk, tblk)
!dir$ attributes inline :: int1_kin_ovl
    TYPE(shpair_t), INTENT(IN) :: cntp
    LOGICAL, INTENT(IN) :: dokinetic
    REAL(KIND=fp), CONTIGUOUS, INTENT(INOUT) :: sblk(:), tblk(:)

    INTEGER :: ig

!dir$ assume_aligned sblk : 64
!dir$ assume_aligned tblk : 64

    DO ig = 1, cntp%numpairs
        CALL comp_kin_ovl_int1_prim(cntp, ig, dokinetic, sblk, tblk)
    END DO

 END SUBROUTINE

!--------------------------------------------------------------------------------

!> @brief Compute contracted block of Coulomb 1e integrals
!> @param[in]       cntp        shell pair data
!> @param[in]       xyz         coordinates of particles
!> @param[in]       c           charges of particles
!> @param[in]       nat         number of particles
!> @param[in]       chgtol      cut-off for charge
!> @param[inout]    blk         block of 1e Coulomb integrals
!> @param[in]       iskp        [opt] particle with this index will be skipped
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!
 SUBROUTINE int1_coul_xyz_c(cntp, xyz, c, nat, chgtol, blk, iskp)
!dir$ attributes inline :: int1_coul_xyz_c
    TYPE(shpair_t), INTENT(IN) :: cntp
    REAL(KIND=fp), CONTIGUOUS, INTENT(IN) :: xyz(:,:), c(:)
    REAL(KIND=fp), INTENT(IN) :: chgtol
    INTEGER, INTENT(IN) :: nat, iskp
    REAL(KIND=fp), CONTIGUOUS, INTENT(INOUT) :: blk(:)

    COMMON /ROOT  / xx,u(13),w(13),nroots
      REAL(KIND=fp) :: xx,u,w
      INTEGER :: nroots
!$omp threadprivate(/ROOT  /)

    INTEGER :: ig, iat

!dir$ assume_aligned blk : 64

    nroots = (cntp%iang+cntp%jang-2)/2 + 1

!   Interaction with point charge
    DO iat = 1, nat
        IF (iat==iskp) CYCLE
        IF (abs(c(iat))<chgtol) CYCLE
        DO ig = 1, cntp%numpairs
            CALL comp_coulomb_int1_prim(cntp, ig, xyz(:,iat), -c(iat), blk)
        END DO
    END DO

 END SUBROUTINE

!--------------------------------------------------------------------------------

!> @brief Compute contracted block of Coulomb 1e integrals
!> @param[in]       cntp        shell pair data
!> @param[in]       xyzc        coordinates and charges of particles
!> @param[in]       nat         number of particles
!> @param[in]       chgtol      cut-off for charge
!> @param[inout]    blk         block of 1e Coulomb integrals
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!
 SUBROUTINE int1_coul_xyzc(cntp, xyzc, nat, chgtol, blk)
!dir$ attributes inline :: int1_coul_xyzc
    TYPE(shpair_t), INTENT(IN) :: cntp
    REAL(KIND=fp), CONTIGUOUS, INTENT(IN) :: xyzc(:)
    REAL(KIND=fp), INTENT(IN) :: chgtol
    INTEGER, INTENT(IN) :: nat
    REAL(KIND=fp), CONTIGUOUS, INTENT(INOUT) :: blk(:)

    COMMON /ROOT  / xx,u(13),w(13),nroots
      REAL(KIND=fp) :: xx,u,w
      INTEGER :: nroots
!$omp threadprivate(/ROOT  /)

    INTEGER :: ig, iat
    REAL(KIND=fp) :: c(3), znuc
!dir$ assume_aligned blk : 64

    nroots = (cntp%iang+cntp%jang-2)/2 + 1

!   Interaction with point charge
    DO iat = 1, nat
        c = xyzc((iat-1)*4+1:iat*4-1)
        znuc = -xyzc(iat*4)
        IF (abs(znuc)<chgtol) CYCLE
        DO ig = 1, cntp%numpairs
            CALL comp_coulomb_int1_prim(cntp, ig, c(:), znuc, blk)
        END DO
    END DO

 END SUBROUTINE

!--------------------------------------------------------------------------------

!> @brief Compute contracted block of Coulomb 1e integrals
!> @param[in]       cntp        shell pair data
!> @param[in]       x           `X` coordinates of charged particles
!> @param[in]       y           `Y` coordinates of charged particles
!> @param[in]       z           `Z` coordinates of charged particles
!> @param[in]       c           charges of particles
!> @param[in]       nat         number of particles
!> @param[in]       chgtol      cut-off for charge
!> @param[inout]    blk         block of 1e Coulomb integrals
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!
 SUBROUTINE int1_coul_x_y_z_c(cntp, x, y, z, c, nat, chgtol, blk)
!dir$ attributes inline :: int1_coul_x_y_z_c
    TYPE(shpair_t), INTENT(IN) :: cntp
    REAL(KIND=fp), CONTIGUOUS, INTENT(IN) :: x(:), y(:), z(:), c(:)
    REAL(KIND=fp), INTENT(IN) :: chgtol
    INTEGER, INTENT(IN) :: nat
    REAL(KIND=fp), CONTIGUOUS, INTENT(INOUT) :: blk(:)

    COMMON /ROOT  / xx,u(13),w(13),nroots
      REAL(KIND=fp) :: xx,u,w
      INTEGER :: nroots
!$omp threadprivate(/ROOT  /)

    INTEGER :: ig, iat
    REAL(KIND=fp) :: crd(3)
!dir$ assume_aligned blk : 64

    nroots = (cntp%iang+cntp%jang-2)/2 + 1

!   Interaction with point charge
    DO iat = 1, nat
        IF (abs(c(iat))<chgtol) CYCLE
        crd(1) = x(iat)
        crd(2) = y(iat)
        crd(3) = z(iat)
        DO ig = 1, cntp%numpairs
            CALL comp_coulomb_int1_prim(cntp, ig, crd, -c(iat), blk)
        END DO
    END DO

 END SUBROUTINE

!--------------------------------------------------------------------------------

!> @brief Compute contracted block of (possibly dumped) Coulomb 1e integrals in FMO method
!> @param[in]       cntp        shell pair data
!> @param[in]       at_info     coordinates charges and scaling factors of particles
!> @param[in]       nat         number of FMO atoms for ESP calculation
!> @param[in]       alpha       dumping exponent
!> @param[in]       beta        dumping factor
!> @param[in]       dodamp      `.true.` to enable damping (based of scaling factors value)
!> @param[in]       dampall     `.true.` to turn dumping of all integrals
!> @param[inout]    vblk        block of 1e Coulomb integrals
!> @param[inout]    zblk        block of 1e Coulomb ESP integrals
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!
 SUBROUTINE int1_fmoesp_ffcxyz(cntp, at_info, nat, &
                 alpha, beta, dodamp, dampall, vblk, zblk)
!dir$ attributes inline :: int1_fmoesp_ffcxyz
    TYPE(shpair_t), INTENT(IN) :: cntp
    REAL(KIND=fp), CONTIGUOUS, INTENT(IN) :: at_info(:)
    INTEGER, INTENT(IN) :: nat
    REAL(KIND=fp), INTENT(IN) :: alpha, beta
    LOGICAL, INTENT(IN) :: dodamp, dampall
    REAL(KIND=fp), CONTIGUOUS, INTENT(INOUT) :: vblk(:), zblk(:)

    COMMON /ROOT  / xx,u(13),w(13),nroots
      REAL(KIND=fp) :: xx,u,w
      INTEGER :: nroots
!$omp threadprivate(/ROOT  /)

    INTEGER :: ig, n
    !fracv - FMO nuclear attraction scaling factor
    !fracesp - FMO ESP scaling factor
    REAL(KIND=fp) :: fracv, fracesp
    REAL(KIND=fp) :: c(3), znuc
    LOGICAL :: damp
!dir$ assume_aligned vblk : 64
!dir$ assume_aligned zblk : 64

    nroots = (cntp%iang+cntp%jang-2)/2 + 1

    DO n = 1, nat
        fracv   = at_info((n-1)*6+1)
        fracesp = at_info((n-1)*6+2)
        znuc    = at_info((n-1)*6+3)
        c(1)    = at_info((n-1)*6+4)
        c(2)    = at_info((n-1)*6+5)
        c(3)    = at_info((n-1)*6+6)

        IF (fracv==ZERO .AND. fracesp==ZERO) CYCLE
!       Damp only BDA atoms unless DAMPALL
        damp = dodamp.AND.fracesp/=ZERO.AND.&
               (fracv/=ZERO .AND. fracv/=1.0 .OR. dampall)

!       Interaction with point charge
        DO ig = 1, cntp%numpairs
            CALL comp_fmoesp_int1_prim(cntp, ig, c, znuc, fracv, fracesp, vblk, zblk)
        END DO

!           The FMO damping term always goes to the ESP (ZBLK)
        IF (damp) THEN
            DO ig = 1, cntp%numpairs
                CALL comp_coulomb_dampch_int1_prim(cntp,ig,alpha,beta,c,znuc*fracesp,zblk)
            END DO
        END IF
    END DO

 END SUBROUTINE

!--------------------------------------------------------------------------------

!> @brief Compute contracted block of (possibly dumped) Coulomb 1e integrals in FMO method
!> @param[in]       cntp        shell pair data
!> @param[in]       ishat       atomic index of the first shell
!> @param[in]       jshat       atomic index of the first shell
!> @param[in]       nat         number of FMO atoms for ESP calculation
!> @param[in]       alpha       dumping exponent
!> @param[in]       beta        dumping factor
!> @param[in]       dodamp      `.true.` to enable damping (based of scaling factors value)
!> @param[in]       dampall     `.true.` to turn dumping of all integrals
!> @param[in]       is_lut      `.true.` if LUT is used
!> @param[inout]    vblk        block of 1e Coulomb integrals
!> @param[inout]    zblk        block of 1e Coulomb ESP integrals
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!
 SUBROUTINE int1_fmoesp_fcall(cntp, ishat, jshat, nat,&
                 alpha, beta, dodamp, dampall, is_lut, vblk, zblk)
!dir$ attributes inline :: int1_fmoesp_fcall
    TYPE(shpair_t), INTENT(IN) :: cntp
    INTEGER, INTENT(IN) :: ishat, jshat, nat
    REAL(KIND=fp), INTENT(IN) :: alpha, beta
    LOGICAL, INTENT(IN) :: dodamp, dampall, is_lut
    REAL(KIND=fp), CONTIGUOUS, INTENT(INOUT) :: vblk(:), zblk(:)

    COMMON /FMOINF/ nfg,nlayer,natfmo,nbdfg,naotyp,nbody
      INTEGER :: nfg,nlayer,natfmo,nbdfg,naotyp,nbody

    COMMON /FMOPNT/ lichfg,lmulfg,lidmrec,lfrgnam,llayfrg,lindat,lncbs,lfmozan,lfmoc,lfmomas,lizbas,liaglob,libdgh, &
          liabdfg,ljabdfg,lncao,lidxcao,liaprjo,ljaprjo,lcoreao,locccor,lshiftb,liodfmo,lfmoda,lfmodb, &
          lfmoespa,lfmoespb,llocfmo,lscffrg,lfmoscf,lrij,lpopmul,lpopmat,lialoc,lindbd,liatfrg,lindfrg, &
          lindgfrg,lnatfrg,lnat0frg,lianfrg,lzanfrg,lcfrg,llibish,llibnsh,llibng,lindatg,lfmobuf(3),lfmode, &
          lnumfrg,lloctat,liaoglob,lloadm,lfmoge,ldgrid,liodcfmo,ljob2grp,lfmopg,lemocdr,luntxyz,luntrot, &
          lstonep,lmapsu,lfrgmul,lclmo,lialmo,lindlmo,latclmo,llmobdf,lfgflmo,lnfglmo,llfglmo,lpfglmo, &
          lpopdmat,lidmpnt,liddpnt,livmpnt,liactfg,lcrfrg,lzlmfrgv,lylmfrgv,lndtfrg,lf_mm,lg_mm,lmaxl30,libuffg
      INTEGER :: lichfg,lmulfg,lidmrec,lfrgnam,llayfrg,lindat,lncbs,lfmozan,lfmoc,lfmomas,lizbas,liaglob,libdgh, &
          liabdfg,ljabdfg,lncao,lidxcao,liaprjo,ljaprjo,lcoreao,locccor,lshiftb,liodfmo,lfmoda,lfmodb, &
          lfmoespa,lfmoespb,llocfmo,lscffrg,lfmoscf,lrij,lpopmul,lpopmat,lialoc,lindbd,liatfrg,lindfrg, &
          lindgfrg,lnatfrg,lnat0frg,lianfrg,lzanfrg,lcfrg,llibish,llibnsh,llibng,lindatg,lfmobuf,lfmode, &
          lnumfrg,lloctat,liaoglob,lloadm,lfmoge,ldgrid,liodcfmo,ljob2grp,lfmopg,lemocdr,luntxyz,luntrot, &
          lstonep,lmapsu,lfrgmul,lclmo,lialmo,lindlmo,latclmo,llmobdf,lfgflmo,lnfglmo,llfglmo,lpfglmo, &
          lpopdmat,lidmpnt,liddpnt,livmpnt,liactfg,lcrfrg,lzlmfrgv,lylmfrgv,lndtfrg,lf_mm,lg_mm,lmaxl30,libuffg

    COMMON /ROOT  / xx,u(13),w(13),nroots
      REAL(KIND=fp) :: xx,u,w
      INTEGER :: nroots
!$omp threadprivate(/ROOT  /)

    COMMON /FMCOM / x(1)
      REAL(KIND=fp) :: x

    REAL(KIND=fp) :: c, fracv, fracesp
    INTEGER :: ig, ic, kfg
    REAL(KIND=fp) :: crd(3)
    LOGICAL :: damp
!dir$ assume_aligned vblk : 64
!dir$ assume_aligned zblk : 64

    fracesp = 0.0
    fracv = 0.0

    nroots = (cntp%iang+cntp%jang-2)/2 + 1

    DO ic = 1, nat

        IF (is_lut .AND. ishat==jshat .AND. ishat==ic) CYCLE

        CALL fmoatfrg(ic,x(lindat),x(lindatg),x(liaglob), &
                      x(lialoc),x(liabdfg),x(ljabdfg), &
                      x(lindbd),x(lfmozan),x(lfmoc), &
                      natfmo+nbdfg,x(luntxyz),x(lpopmat), &
                      ishat,jshat,fracv,fracesp, &
                      kfg,c,crd(1),crd(2),crd(3))

        IF (fracv==ZERO .AND. fracesp==ZERO) CYCLE

!       Damp only BDA atoms unless DAMPALL
        damp = dodamp.AND.fracesp/=ZERO.AND.&
               (fracv/=ZERO .AND. fracv/=1.0 .OR. dampall)

!       Interaction with point charge
        DO ig = 1, cntp%numpairs
            CALL comp_fmoesp_int1_prim(cntp, ig, crd, c, fracv, fracesp, vblk, zblk)
        END DO

!       The FMO damping term always goes to the ESP (ZBLK)
        IF (damp) THEN
            DO ig = 1, cntp%numpairs
                CALL comp_coulomb_dampch_int1_prim(cntp,ig,alpha,beta,crd,c*fracesp,zblk)
            END DO
        END IF

    END DO

 END SUBROUTINE

!--------------------------------------------------------------------------------

!> @brief Compute sum of Coulomb integrals over pair of contracted shells
!> @param[in]       cntp        shell pair data
!> @param[in]       x           `X` coordinate of the charged particle
!> @param[in]       y           `Y` coordinate of the charged particle
!> @param[in]       z           `Z` coordinate of the charged particle
!> @param[in]       den         normalized density matrix block
!> @param[inout]    vsum        sum of Coulomb integrals over shell pair
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Oct, 2018_ Initial release
!
 SUBROUTINE int1_epoten(cntp, x, y, z, den, vsum)
!dir$ attributes inline :: int1_epoten
    TYPE(shpair_t), INTENT(IN) :: cntp
    REAL(KIND=fp), INTENT(IN) :: x, y, z, den(:)
    REAL(KIND=fp), INTENT(OUT) :: vsum

    COMMON /ROOT  / xx,u(13),w(13),nroots
      REAL(KIND=fp) :: xx,u,w
      INTEGER :: nroots
!$omp threadprivate(/ROOT  /)

    INTEGER :: ig
    REAL(KIND=fp) :: crd(3)

    vsum = 0.0
    crd(1) = x
    crd(2) = y
    crd(3) = z

    nroots = (cntp%iang+cntp%jang-2)/2 + 1

    DO ig = 1, cntp%numpairs
!       Interaction with point charge
        CALL comp_coulpot_prim(cntp, ig, crd, den, vsum)
    END DO

 END SUBROUTINE

!--------------------------------------------------------------------------------

!> @brief Compute contracted block of Z-angular momentum integrals
!> @param[in]       cntp        shell pair data
!> @param[inout]    blk         block of 1e Coulomb Lz-integrals
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!
 SUBROUTINE int1_lz(cntp, blk)
!dir$ attributes inline :: int1_lz
    TYPE(shpair_t), INTENT(IN) :: cntp
    REAL(KIND=fp), CONTIGUOUS, INTENT(INOUT) :: blk(:)

    INTEGER :: ig
!dir$ assume_aligned blk : 64

    DO ig = 1, cntp%numpairs
        CALL comp_lz_int1_prim(cntp, ig, blk)
    END DO

 END SUBROUTINE

END MODULE
