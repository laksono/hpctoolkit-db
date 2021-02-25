#define __TIMING 0
!---------------------------------------------------------------------
!=====================================================================
!---------------------------------------------------------------------
!*MODULE OMPMOD
!>  @brief   In this module OpenMP versions of Fock matrixIJPRIM_gpu_2 calculation routines are collected.
!>           They are called explicitly from corresponding routines from int2a.src module,
!>           namely TWOEI and EXCHNG, using conditional OpenMP compiling. They are never
!>           called if compiler doesn't support OpenMP or -openmp/-fopenmp/-qopenmp flag
!>           is not set at compile time for the all files.
!>  @details A lot of the work on thread-safety is done by small script "addomp.sh" which
!>           is located in ../tools directory. It sets <threadprivate> attribute
!>           for sensitive common blocks throughout all files in GAMESS source directory.
!>  @author  Vladimir Mironov
!>  @date    2016-2017


!---------------------------------------------------------------------
!=====================================================================
!---------------------------------------------------------------------




MODULE ompmod

    USE omp_lib
    USE prec, ONLY: fp

    USE mx_limits, ONLY: &
        mxsh, mxgtot, mxatm, mxang, mxang2, mxgsh, mxg2

    IMPLICIT NONE

    PRIVATE

    PUBLIC &
        ompmod_twoei_jk, &
        ompmod_twoei_kl, &
        ompmod_twoei_shf_kl_rhf, &
        ompmod_exchng, &
        ompmod_rysint, &
        ompmod_raxintsp, &
        ompmod_raxintspd, &
        ompmod_shellquart

    CHARACTER(LEN=*), PARAMETER :: &
        dbgfmt1 = '(/2x,&
                    &"Thread | Number of |",19X,"Timing",&
                    &/1x," number | quartets  |  Integrals |   F update  ",&
                    &"|   Schwartz  |    Total    |")', &
        dbgfmt2 = '(i5,4x,"|",i10," |",4(f9.2," s | "))', &

        dbgfmt_exch ='(1X,"SCHWARZ INEQUALITY OVERHEAD:",I10," INTEGRALS, &
                  &T=",F12.2)'

    INTEGER, PARAMETER :: &
        angm(0:6) = (/4,4,6,10,15,21,28/)

    REAL(KIND=fp), PARAMETER :: &
        zero=0.0_fp, one=1.0_fp

    INTERFACE
        REAL(KIND=fp) FUNCTION schwdn(dsh,ish,jsh,ksh,lsh,ia)
            USE prec
            INTEGER, INTENT(IN) :: ish, jsh, ksh, lsh, ia(*)
            REAL(KIND=fp), INTENT(IN) :: dsh(*)
        END FUNCTION
    END INTERFACE

CONTAINS

!*MODULE OMPMOD   *DECK OMPMOD_TWOEI_JK
!
!>    @brief   Threaded version of two-electron integral
!>             calculation routine. Load balance of OpenMP code
!>             is done over J and K shell indices,
!>             MPI load balance - over I index.
!>             Related input file option:
!>             `$INTGRL INTOMP=1 $END`
!>
!>    @details Calculates two-electron contribution to the Fock
!>             matrix. Based on `TWOEI` subroutine from `int2a.src`
!>
!>    @author  Vladimir Mironov
!
!     REVISION HISTORY:
!>    @date _May, 2016_ Initial release
!>    @date _Jun, 2017_ Bug fixes
!
!     PARAMETERS:
!
!>    @param[in]     typscf    hollerith, the kind of SCF (=RHF,UHF, etc.)
!>    @param[in]     schwrz    determines whether to use integral screening
!>    @param[out]    nint      number of calculated ERIs
!>    @param[out]    nschwz    number of ERIs that were screened out
!>    @param[in]     l1        size of index array
!>    @param[in]     l2a       size of Fock and density matrices (alpha)
!>    @param[in]     l2b       size of Fock and density matrices (beta)
!>    @param[in]     xints(:)  array of exchange integrals
!>    @param[in]     nsh2      size of `xints(:)` array
!>    @param[in]     maxg      size of temporary array for Rys code
!>    @param[in]     ia(:)     index array, contains "triangular numbers"
!>    @param[in,out] da(:)     density matrix for alpha electrons
!>    @param[out]    fa(:)     Fock matrix for alpha electrons
!>    @param[in,out] db(:)     density matrix for beta electrons,
!>                                 used in open-shell calculations only
!>    @param[out]    fb(:)     Fock matrix for beta electrons,
!>                                 used in open-shell calculations only
!>    @param[in]     dsh       density matrix packed in shells
!>                                 for screening purposes
!>    @param[in]     nflmat    when >1 selects CPHF calculation
!>    @param[in]     cutoff    cutoff for integral screening
!>    @param[in]     oflag     logical parameter for debug timing output
      SUBROUTINE ompmod_twoei_jk(typscf,schwrz,nint,nschwz, &
                             l1,l2a,l2b,xints, &
                             nsh2,maxg, &
                             ia,da,fa,db,fb,dsh,nflmat, &
                             cutoff, oflag)
      use omp_lib
      REAL(KIND=8),PARAMETER :: ZER=0.0D00
      LOGICAL, INTENT(IN) :: schwrz 
      INTEGER, INTENT(OUT) :: nint, nschwz
      INTEGER, INTENT(IN) :: &
        l1, l2a, l2b, nsh2, nflmat, maxg, &
        ia(l1)

!   typscf is an 8-byte hollerith constant
      REAL(KIND=8), INTENT(IN) :: typscf  
      REAL(KIND=fp), INTENT(IN) :: cutoff, xints(nsh2), dsh(nsh2)  
      REAL(KIND=fp), INTENT(INOUT) :: da(l2a), db(l2b)  
      REAL(KIND=fp), INTENT(OUT) :: fa(l2a*nflmat), fb(l2b*nflmat)
      LOGICAL, INTENT(IN) :: oflag
 
      COMMON /par   / me,master,nproc,ibtyp,iptim,goparr,dskwrk,maswrk
          INTEGER :: me,master,nproc,ibtyp,iptim
          LOGICAL :: dskwrk, maswrk, goparr  
      COMMON /iofile/ ir,iw,ip,is,ipk,idaf,nav,ioda(950)
          INTEGER :: ir,iw,ip,is,ipk,idaf,nav,ioda
      COMMON /fmcom / xx(1)
          REAL(KIND=fp), TARGET :: xx
      COMMON /NSHEL / ex(mxgtot),cs(mxgtot),cp(mxgtot),cd(mxgtot),     &
                      cf(mxgtot),cg(mxgtot),ch(mxgtot),ci(mxgtot),     &
                      kstart(mxsh),katom(mxsh),ktype(mxsh),kng(mxsh),  &
                      kloc(mxsh),kmin(mxsh),kmax(mxsh),nshell
          INTEGER :: kstart,katom,ktype,kng,kloc,kmin,kmax,nshell
          REAL(KIND=fp) :: ex,cs,cp,cd,cf,cg,ch,ci

      !to get rid of nshell
      INTEGER,ALLOCATABLE,DIMENSION(:) :: kstart_tmp,katom_tmp,ktype_tmp,kng_tmp,kloc_tmp,kmin_tmp,kmax_tmp
      INTEGER :: nshell_tmp
      REAL(KIND=fp),ALLOCATABLE,DIMENSION(:) :: ex_tmp,cs_tmp,cp_tmp,cd_tmp,cf_tmp,cg_tmp,ch_tmp,ci_tmp

      COMMON /shlnos/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl, &
                    mini,minj,mink,minl,maxi,maxj,maxk,maxl, &
                    nij,ij,kl,ijkl
        INTEGER :: lit,ljt,lkt,llt,loci,locj,lock,locl, &
                   mini,minj,mink,minl,maxi,maxj,maxk,maxl, &
                   nij,ij,kl,ijkl
        REAL(KIND=fp) :: qq4
      COMMON /maxc  / cmax(mxgtot),cmaxa(mxgsh),cmaxb(mxgsh), &
                      cmaxc(mxgsh),cmaxd(mxgsh),ismlp(mxg2),ismlq
          REAL(KIND=fp) :: cmax,cmaxa,cmaxb,cmaxc,cmaxd
          INTEGER :: ismlp,ismlq
      COMMON /nlrcf / lrint
          LOGICAL :: lrint
      COMMON /dftpar/ dfttyp(20),exena,exenb,exenc, &
                      idft34,nauxfun,nauxshl
          REAL(KIND=fp) :: dfttyp,exena,exenb,exenc
          INTEGER :: idft34,nauxfun,nauxshl

!   Common blocks /MAXC  /, /SHLNOS/ and /GOUT/ contain
!   both global (cmax(:),qq4,norgp) and thread-local
!   data (the rest). Set them threadprivate and copyin
!   global data later in $omp parallel section.

!$omp threadprivate(/shlnos/,/maxc  /,/nlrcf /,/dftpar/)
!!$omp threadprivate(/maxc  /,/nlrcf /,/dftpar/)

      INTEGER :: &
        next, num_threads, ithread, thr_nshq, &
        ii, jj, kk, ll, ijij, klkl, jork,MM
  
      LOGICAL :: schskp, dirscf
      REAL(KIND=fp) :: &
        denmax, tim, tim0, tim1, tim2, tim3, tim4, test, tim2222s, tim2222e
    
      REAL(KIND=fp),DIMENSION(:),ALLOCATABLE :: ddij, ghondo
  
      INTEGER :: shlmxang, shlsumang, shltotcon
      LOGICAL :: lshell

      LOGICAL :: DO_raxintsp,DO_raxintspd,DO_eric_ts,DO_rysint        
      INTEGER :: NCOUNT_raxintspd,NCOUNT_eric_ts,NCOUNT_rysint
      INTEGER,ALLOCATABLE,DIMENSION(:,:) :: &
         INDEX_raxintspd,INDEX_eric_ts,INDEX_rysint

      INTEGER :: ICOUNT_raxintspd,ICOUNT_eric_ts,ICOUNT_rysint

      INTEGER :: NSA,ISA !sum angular momentum
      INTEGER :: NCP,ICP,ACP(100) !product contraction
      INTEGER :: MAXSP !max number of sp quartets
      INTEGER,ALLOCATABLE :: SPTYP(:,:,:) !type of SP quartets
      INTEGER,ALLOCATABLE :: SPSHINDEX(:,:,:,:,:) !SP shell index

      !rys
      INTEGER :: NSA_rys,ISA_rys !sum angular momentum
      INTEGER :: NCP_rys,ICP_rys,ACP_rys(100) !product contraction
      INTEGER :: MAXRYS !max number of rys quartets
      INTEGER,ALLOCATABLE :: RYSTYP(:,:) !type of rys quartets
      INTEGER,ALLOCATABLE :: RYSINDEX(:,:,:,:) !rys shell index
      INTEGER :: MM_rys
      INTEGER :: INEW,JNEW,KNEW,LNEW
      INTEGER :: IS1,JS1,KS1,LS1
! threadprivate common blocks ==========================================
      !shlnos
      INTEGER :: LOCI_rys,LOCJ_rys,LOCK_rys,LOCL_rys, &
                 MINI_rys,MINJ_rys,MINK_rys,MINL_rys,&
                 MAXI_rys,MAXJ_rys,MAXK_rys,MAXL_rys,&
                 lit_rys,ljt_rys,lkt_rys,llt_rys,&
                 nij_rys,ij_rys,kl_rys,ijkl_rys
      !REAL(KIND=fp) :: qq4
      !shlinf
      REAL(KIND=fp) :: GA(MXGSH),CSA1(MXGSH),CPA1(MXGSH),CDA(MXGSH),&
                       CFA(MXGSH),CGA(MXGSH),CHA(MXGSH),CIA(MXGSH),&
                       GB(MXGSH),CSB1(MXGSH),CPB1(MXGSH),CDB(MXGSH),&
                       CFB(MXGSH),CGB(MXGSH),CHB(MXGSH),CIB(MXGSH),&
                       GC(MXGSH),CSC1(MXGSH),CPC1(MXGSH),CDC(MXGSH),&
                       CFC(MXGSH),CGC(MXGSH),CHC(MXGSH),CIC(MXGSH),&
                       GD(MXGSH),CSD1(MXGSH),CPD1(MXGSH),CDD(MXGSH),&
                       CFD(MXGSH),CGD(MXGSH),CHD(MXGSH),CID(MXGSH)
      REAL(KIND=fp) :: AX,AY,AZ,BX,BY,BZ,RAB,CX,CY,CZ,DX,DY,DZ,RCD1
      INTEGER :: NGA1,NGB1,NGC1,NGD1    
! END threadprivate common blocks =======================================
          REAL(KIND=fp),DIMENSION(84) :: &
     IX,IY,IZ,JX,JY,JZ,KX,KY,KZ,LX,LY,LZ

     DATA JX /   0,  49,   0,   0,  98,   0,   0,  49,  49,   0, &
         147,   0,   0,  98,  98,  49,   0,  49,   0,  49, &
         196,   0,   0, 147, 147,  49,   0,  49,   0,  98, &
          98,   0,  98,  49,  49, &
         245,   0,   0, 196, 196,  49,   0,  49,   0, 147, &
         147,  98,   0,  98,   0, 147,  49,  49,  98,  98, &
          49, &
         294,   0,   0, 245, 245,  49,   0,  49,   0, 196, &
         196,  98,   0,  98,   0, 196,  49,  49, 147, 147, &
           0, 147, 147,  98,  49,  98,  49,  98/

      DATA IX /   1, 344,   1,   1, 687,   1,   1, 344, 344,   1, &
              1030,   1,   1, 687, 687, 344,   1, 344,   1, 344, &
              1373,   1,   1,1030,1030, 344,   1, 344,   1, 687, &
               687,   1, 687, 344, 344, &
              1716,   1,   1,1373,1373, 344,   1, 344,   1,1030, &
              1030, 687,   1, 687,   1,1030, 344, 344, 687, 687, &
               344, &
              2059,   1,   1,1716,1716, 344,   1, 344,   1,1373, &
              1373, 687,   1, 687,   1,1373, 344, 344,1030,1030, &
                 1,1030,1030, 687, 344, 687, 344, 687/
      DATA JY /   0,   0,  49,   0,   0,  98,   0,  49,   0,  49, &
                 0, 147,   0,  49,   0,  98,  98,   0,  49,  49, &
                 0, 196,   0,  49,   0, 147, 147,   0,  49,  98, &
                 0,  98,  49,  98,  49, &
                 0, 245,   0,  49,   0, 196, 196,   0,  49,  98, &
                 0, 147, 147,   0,  98,  49, 147,  49,  98,  49, &
                98, &
                 0, 294,   0,  49,   0, 245, 245,   0,  49,  98, &
                 0, 196, 196,   0,  98,  49, 196,  49, 147,   0, &
               147,  98,  49, 147, 147,  49,  98,  98/
      DATA IY /   1,   1, 344,   1,   1, 687,   1, 344,   1, 344, &
                 1,1030,   1, 344,   1, 687, 687,   1, 344, 344, &
                 1,1373,   1, 344,   1,1030,1030,   1, 344, 687, &
                 1, 687, 344, 687, 344, &
                 1,1716,   1, 344,   1,1373,1373,   1, 344, 687, &
                 1,1030,1030,   1, 687, 344,1030, 344, 687, 344, &
               687, &
                 1,2059,   1, 344,   1,1716,1716,   1, 344, 687, &
                 1,1373,1373,   1, 687, 344,1373, 344,1030,   1, &
              1030, 687, 344,1030,1030, 344, 687, 687/
      DATA JZ /   0,   0,   0,  49,   0,   0,  98,   0,  49,  49, &
                 0,   0, 147,   0,  49,   0,  49,  98,  98,  49, &
                 0,   0, 196,   0,  49,   0,  49, 147, 147,   0, &
                98,  98,  49,  49,  98, &
                 0,   0, 245,   0,  49,   0,  49, 196, 196,   0, &
                98,   0,  98, 147, 147,  49,  49, 147,  49,  98, &
                98, &
                 0,   0, 294,   0,  49,   0,  49, 245, 245,   0, &
                98,   0,  98, 196, 196,  49,  49, 196,   0, 147, &
               147,  49,  98,  49,  98, 147, 147,  98/
      DATA IZ /   1,   1,   1, 344,   1,   1, 687,   1, 344, 344, &
                 1,   1,1030,   1, 344,   1, 344, 687, 687, 344, &
                 1,   1,1373,   1, 344,   1, 344,1030,1030,   1, &
               687, 687, 344, 344, 687, &
                 1,   1,1716,   1, 344,   1, 344,1373,1373,   1, &
               687,   1, 687,1030,1030, 344, 344,1030, 344, 687, &
               687, &
                 1,   1,2059,   1, 344,   1, 344,1716,1716,   1, &
               687,   1, 687,1373,1373, 344, 344,1373,   1,1030, &
              1030, 344, 687, 344, 687,1030,1030, 687/
      DATA LX /   0,   1,   0,   0,   2,   0,   0,   1,   1,   0,&
                 3,   0,   0,   2,   2,   1,   0,   1,   0,   1,&
                 4,   0,   0,   3,   3,   1,   0,   1,   0,   2,&
                 2,   0,   2,   1,   1,&
                 5,   0,   0,   4,   4,   1,   0,   1,   0,   3,&
                 3,   2,   0,   2,   0,   3,   1,   1,   2,   2,&
                 1,&
                 6,   0,   0,   5,   5,   1,   0,   1,   0,   4,&
                 4,   2,   0,   2,   0,   4,   1,   1,   3,   3,&
                 0,   3,   3,   2,   1,   2,   1,   2/
      DATA KX /   0,   7,   0,   0,  14,   0,   0,   7,   7,   0,&
                21,   0,   0,  14,  14,   7,   0,   7,   0,   7,&
                28,   0,   0,  21,  21,   7,   0,   7,   0,  14,&
                14,   0,  14,   7,   7,&
                35,   0,   0,  28,  28,   7,   0,   7,   0,  21,&
                21,  14,   0,  14,   0,  21,   7,   7,  14,  14,&
                 7,&
                42,   0,   0,  35,  35,   7,   0,   7,   0,  28,&
                28,  14,   0,  14,   0,  28,   7,   7,  21,  21,&
                 0,  21,  21,  14,   7,  14,   7,  14/
     DATA LY /   0,   0,   1,   0,   0,   2,   0,   1,   0,   1,&
                 0,   3,   0,   1,   0,   2,   2,   0,   1,   1,&
                 0,   4,   0,   1,   0,   3,   3,   0,   1,   2,&
                 0,   2,   1,   2,   1,&
                 0,   5,   0,   1,   0,   4,   4,   0,   1,   2,&
                 0,   3,   3,   0,   2,   1,   3,   1,   2,   1,&
                 2,&
                 0,   6,   0,   1,   0,   5,   5,   0,   1,   2,&
                 0,   4,   4,   0,   2,   1,   4,   1,   3,   0,&
                 3,   2,   1,   3,   3,   1,   2,   2/
      DATA KY /   0,   0,   7,   0,   0,  14,   0,   7,   0,   7,&
                 0,  21,   0,   7,   0,  14,  14,   0,   7,   7,&
                 0,  28,   0,   7,   0,  21,  21,   0,   7,  14,&
                 0,  14,   7,  14,   7,&
                 0,  35,   0,   7,   0,  28,  28,   0,   7,  14,&
                 0,  21,  21,   0,  14,   7,  21,   7,  14,   7,&
                14,&
                 0,  42,   0,   7,   0,  35,  35,   0,   7,  14,&
                 0,  28,  28,   0,  14,   7,  28,   7,  21,   0,&
                21,  14,   7,  21,  21,   7,  14,  14/
      DATA LZ /   0,   0,   0,   1,   0,   0,   2,   0,   1,   1,&
                 0,   0,   3,   0,   1,   0,   1,   2,   2,   1,&
                 0,   0,   4,   0,   1,   0,   1,   3,   3,   0,&
                 2,   2,   1,   1,   2,&
                 0,   0,   5,   0,   1,   0,   1,   4,   4,   0,&
                 2,   0,   2,   3,   3,   1,   1,   3,   1,   2,&
                 2,&
                 0,   0,   6,   0,   1,   0,   1,   5,   5,   0,&
                 2,   0,   2,   4,   4,   1,   1,   4,   0,   3,&
                 3,   1,   2,   1,   2,   3,   3,   2/
      DATA KZ /   0,   0,   0,   7,   0,   0,  14,   0,   7,   7,&
                 0,   0,  21,   0,   7,   0,   7,  14,  14,   7,&
                 0,   0,  28,   0,   7,   0,   7,  21,  21,   0,&
                14,  14,   7,   7,  14,&
                 0,   0,  35,   0,   7,   0,   7,  28,  28,   0,&
                14,   0,  14,  21,  21,   7,   7,  21,   7,  14,&
                14,&
                 0,   0,  42,   0,   7,   0,   7,  35,  35,   0,&
                14,   0,  14,  28,  28,   7,   7,  28,   0,  21,&
                21,   7,  14,   7,  14,  21,  21,  14/

      COMMON /shlexc/ norgsh(3),norgsp(3),iexch,nangm,ngth(4)
        INTEGER :: norgsh,norgsp,iexch,nangm,ngth
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB, &
                    ZAN(MXATM),C(3,MXATM),IAN(MXATM)
        INTEGER :: NAT,ICH,MUL,NUM,NQMT,NE,NA,NB
        REAL(KIND=fp) :: ZAN,IAN,C

      ! ERIOUT
      INTEGER :: lstri,lstrj,lstrk,lstrl

      COMMON /SHLT  / TOL,CUTOFFAO,ICOUNT,OUT
      REAL(KIND=fp) :: TOL,CUTOFFAO,OUT
      INTEGER :: ICOUNT

      REAL(KIND=fp) :: HFSCAL,CSCALT
      INTEGER :: JTYPE
      INTEGER :: max_threads

! SP additional data ===================================================
      COMMON /SHLLFO/ NGA,LA,EXA(MXGSH),CSA(MXGSH),CPA(MXGSH), &
                      NGB,LB,EXB(MXGSH),CSB(MXGSH),CPB(MXGSH), &
                      NGC,LC,EXC(MXGSH),CSC(MXGSH),CPC(MXGSH), &
                      NGD,LD,EXD(MXGSH),CSD(MXGSH),CPD(MXGSH)
      REAL(KIND=fp) :: EXA,CSA,CPA,EXB,CSB,CPB,EXC,CSC,CPC,EXD,CSD,CPD      
      INTEGER :: NGA,LA,NGB,LB,NGC,LC,NGD,LD
!$omp threadprivate(/SHLLFO/)   
      COMMON /KI3 / R00(25),R01(120),R02(156),R03(80),R04(15)
         REAL(KIND=fp) :: R00,R01,R02,R03,R04
      COMMON /eridat/ len1,len2,len3,len4
         INTEGER :: len1,len2,len3,len4
      REAL(KIND=fp) :: gpople(768)
      !INTEGER :: INEW,JNEW,KNEW,LNEW
      INTEGER :: norgp
      INTEGER :: LAT,LBT,LCT,LDT,ITYPE
      REAL(KIND=fp) :: RCD,SING,COSG
      INTEGER :: lpopi,lpopj,lpopk,lpopl

      REAL(KIND=fp) :: P12(3,3),P34(3,3),P(3,3),T(3)
      REAL(KIND=fp) :: R34,ACX,ACZ

!$omp threadprivate(/KI3/)
      INTEGER :: IKL,IG,I,J
      REAL(KIND=fp) :: QX,QZ,TMP
! END SP additional data ===============================================
!NEW RYS STUFF
      INTEGER,PARAMETER :: NUM_QUART=9169597 !number of quartets
      INTEGER :: NUM_QUART_0000
      INTEGER,ALLOCATABLE :: IDX_RYS_2(:,:,:),N_RYS_2(:)

      INTEGER,ALLOCATABLE :: IDX_RYS_2000(:,:,:),N_RYS_2000(:)
      INTEGER,ALLOCATABLE :: IDX_RYS_2100(:,:,:),N_RYS_2100(:)

      INTEGER,ALLOCATABLE :: IDX_RYS_3(:,:,:),N_RYS_3(:)
      INTEGER,ALLOCATABLE :: IDX_RYS_4(:,:,:),N_RYS_4(:)
      INTEGER,ALLOCATABLE :: IDX_RYS_5(:,:,:),N_RYS_5(:)

      DOUBLE PRECISION :: varx(31213,1000),vary(31213,1000),varz(31213,1000)
      DOUBLE PRECISION :: XIN(31213),YIN(31213),ZIN(31213)

!
!   --- Initialization of variables ---
!

CALL OMP_SET_DYNAMIC(.FALSE.)

      nint   = 0
      nschwz = 0
      schskp = .FALSE.
      denmax = zero
      dirscf = .TRUE.

!     no symmetry supported yet
      qq4 = 1
      norgp = 0

      !allocate temporary stuff
      ALLOCATE(ex_tmp(mxgtot))
      ALLOCATE(cs_tmp(mxgtot))
      ALLOCATE(cp_tmp(mxgtot))
      ALLOCATE(cd_tmp(mxgtot))
      ALLOCATE(cf_tmp(mxgtot))
      ALLOCATE(cg_tmp(mxgtot))
      ALLOCATE(ch_tmp(mxgtot))
      ALLOCATE(ci_tmp(mxgtot))
      ALLOCATE(kstart_tmp(mxsh))
      ALLOCATE(katom_tmp(mxsh))
      ALLOCATE(ktype_tmp(mxsh))
      ALLOCATE(kng_tmp(mxsh))
      ALLOCATE(kloc_tmp(mxsh))
      ALLOCATE(kmin_tmp(mxsh))
      ALLOCATE(kmax_tmp(mxsh))
      !new
      nshell_tmp=nshell
      ex_tmp=ex
      cs_tmp=cs
      cp_tmp=cp
      cd_tmp=cd
      cf_tmp=cf
      cg_tmp=cg
      ch_tmp=ch
      ci_tmp=ci
      kstart_tmp=kstart
      katom_tmp=katom
      ktype_tmp=ktype
      kng_tmp=kng
      kloc_tmp=kloc
      kmin_tmp=kmin
      kmax_tmp=kmax

      ! allocate temporary memory for integrals
      ALLOCATE(ddij(mxang2*mxg2),ghondo(maxg))


      ! ================================================================
      ! START ===== rotated axis SP ============================== START
      ! ================================================================


      HFSCAL=DFTTYP(3)
      CSCALT=1.0D00

      ! number of quartets for SP, SPD, ERIC and RYS
      CALL RHF_NCOUNT_raxintsp_gpu &
          (MAXSP,ACP,NCP,dirscf,schwrz,cutoff,nschwz, &
           nsh2,l1,typscf,xints,dsh,nint,ia)


      ! do SP?
      DO_raxintsp=.FALSE.
      IF(MAXSP.GT.0) DO_raxintsp=.TRUE.

      ! allocate SP index array      
      IF(DO_raxintsp) THEN
         ALLOCATE(SPTYP(6,NCP,5))
         ALLOCATE(SPSHINDEX(4,MAXSP,6,NCP,5))
      ENDIF


      ! shell information to index arrays
      IF(DO_raxintsp) THEN      
        CALL RHF_INDEX_raxintsp_gpu &
            (SPSHINDEX,ACP,NCP,SPTYP,MAXSP,     &
             dirscf,schwrz,cutoff,nschwz,nsh2,  &
             l1,typscf,xints,dsh,nint,ia)
      ENDIF


      ! calculate SP quartets (ghondo)
      ! and fock contribution (fa)
      NSA=5
      IF(DO_raxintsp) THEN
         DO ISA=1,NSA
            DO ICP=1,NCP
               DO JTYPE=1,6

!$omp parallel do default(NONE) &
!$omp reduction(+:fa) &
!$omp shared(SPTYP,SPSHINDEX,da,ia) &
!$omp shared(HFSCAL,CSCALT,CUTOFFAO) &
!$omp shared(maxg,norgp,IKL) &
!$omp shared(JTYPE,ICP,ISA) &
!$omp shared(kmax,kmin,ktype) &
!$omp private(ghondo,gpople,MM,ii,jj,kk,ll,nint) &
!$omp private(lstri,lstrj,lstrk,lstrl) &
!$omp private(INEW,JNEW,KNEW,LNEW) &
!$omp private(LAT,LBT,LCT,LDT,ITYPE) &
!$omp private(RCD,SING,COSG,P) &
!$omp private(lpopi,lpopj,lpopk,lpopl) &
!$omp private(QX,QZ,TMP,IG) &
!$omp private(P12,P34,R34,T,ACX,ACZ)

                  DO MM=1,SPTYP(JTYPE,ICP,ISA)

                     ii=SPSHINDEX(1,MM,JTYPE,ICP,ISA)
                     jj=SPSHINDEX(2,MM,JTYPE,ICP,ISA)
                     kk=SPSHINDEX(3,MM,JTYPE,ICP,ISA)
                     ll=SPSHINDEX(4,MM,JTYPE,ICP,ISA)

                     mini = kmin(ii)
                     maxi = kmax(ii)
                     minj = kmin(jj)
                     maxj = kmax(jj)
                     mink = kmin(kk)
                     maxk = kmax(kk)
                     minl = kmin(ll)
                     maxl = kmax(ll)

                     LAT= ktype(ii)-1
                     LBT= ktype(jj)-1
                     LCT= ktype(kk)-1
                     LDT= ktype(ll)-1
                     ITYPE= 1+LDT+2*(LCT+2*(LBT+2*LAT))

                     CALL genr70_switch_gpu &
                         (ii,jj,kk,ll, &
                          ITYPE,LAT,LBT,LCT,LDT, &
                          LA,LB,LC,LD, &
                          INEW,JNEW,KNEW,LNEW, &
                          LPOPI,LPOPJ,LPOPK,LPOPL)

                     IKL= 0
                     IF(JTYPE.EQ. 1) THEN
                        R00(1)= ZER
                     ELSEIF(JTYPE.EQ. 2) THEN
                        CALL INTK2_1_gpu(IKL,R00,R01)
                     ELSEIF(JTYPE.EQ. 3) THEN
                        CALL INTK3_1_gpu(IKL,R00,R01,R02)
                     ELSEIF(JTYPE.EQ. 4) THEN
                        CALL INTK4_1_gpu(IKL,R00,R01,R02)
                     ELSEIF(JTYPE.EQ. 5) THEN
                        CALL INTK5_1_gpu(IKL,R00,R01,R02,R03)
                     ELSEIF(JTYPE.EQ. 6) THEN
                        CALL INTK6_1_gpu(IKL,R00,R01,R02,R03,R04)
                     ENDIF

                     CALL genr70_pre_gpu &
                         (JTYPE,ii,jj,kk,ll, &
                          INEW,JNEW,KNEW,LNEW, &
                          P12,P34,R34,ACX,ACZ, &
                          RCD,SING,COSG,P,T)

                     CALL genr70_P_gpu &
                         (JTYPE,ii,jj,kk,ll, &
                          INEW,JNEW,KNEW,LNEW, &
                          P12,P34,R34,ACX,ACZ, &
                          RCD,SING,COSG,P,T)

                     IF(JTYPE.EQ.1) THEN
                        CALL genr70_Q1_gpu(JTYPE,P12,R34,ACX,ACZ,RCD,SING,COSG)
                     ELSEIF(JTYPE.EQ.2) THEN
                        CALL genr70_Q2_gpu(JTYPE,P12,R34,ACX,ACZ,RCD,SING,COSG)
                     ELSEIF(JTYPE.EQ.3) THEN
                        CALL genr70_Q3_gpu(JTYPE,P12,R34,ACX,ACZ,RCD,SING,COSG)
                     ELSEIF(JTYPE.EQ.4) THEN
                        CALL genr70_Q4_gpu(JTYPE,P12,R34,ACX,ACZ,RCD,SING,COSG)
                     ELSEIF(JTYPE.EQ.5) THEN
                        CALL genr70_Q5_gpu(JTYPE,P12,R34,ACX,ACZ,RCD,SING,COSG)
                     ELSEIF(JTYPE.EQ.6) THEN
                        CALL genr70_Q6_gpu(JTYPE,P12,R34,ACX,ACZ,RCD,SING,COSG)
                     ENDIF

                     IG= 1+NORGP
                     IF(JTYPE.EQ.1) THEN
                        GPOPLE(IG)= R00(1)
                     ENDIF

                     QX= RCD*SING
                     QZ= RCD*COSG
                     IF(JTYPE.EQ. 2)THEN
                        CALL MCDV2_gpu(GPOPLE(IG),QX,QZ,R00,R01)
                     ELSEIF(JTYPE.EQ. 3) THEN
                        CALL MCDV3_gpu(GPOPLE(IG),QX,QZ,R00,R01,R02)
                     ELSEIF(JTYPE.EQ. 4) THEN
                        CALL MCDV4_gpu(GPOPLE(IG),QX,QZ,R00,R01,R02)
                     ELSEIF(JTYPE.EQ. 5) THEN
                        CALL MCDV5_gpu(GPOPLE(IG),QX,QZ,R00,R01,R02,R03)
                     ELSEIF(JTYPE.EQ. 6) THEN
                        CALL MCDV6_gpu(GPOPLE(IG),QX,QZ,R00,R01,R02,R03,R04)
                     ENDIF
!
!JMS  NOW, THE TRANSPOSE OF P TO BE USED FOR COMPUTATIONAL EFFICIENCY
!
                     DO J=1,2
                        DO I=J+1,3
                           TMP= P(I,J)
                           P(I,J)= P(J,I)
                           P(J,I)= TMP
                        ENDDO
                     ENDDO

                     IF(JTYPE.EQ.2) THEN
                        CALL R30S1S_J2_gpu(JTYPE,GPOPLE(IG),P)
                     ELSEIF(JTYPE.EQ.3) THEN
                        CALL R30S1S_J3_gpu(JTYPE,GPOPLE(IG),P)
                     ELSEIF(JTYPE.EQ.4) THEN
                        CALL R30S1S_J4_gpu(JTYPE,GPOPLE(IG),P)
                     ELSEIF(JTYPE.EQ.5) THEN
                        CALL R30S1S_J5_gpu(JTYPE,GPOPLE(IG),P)
                     ELSE
                        CALL R30S1S_J6_gpu(JTYPE,GPOPLE(IG),P)
                     ENDIF

                     CALL raxintsp_cpyint_gpu &
                         (lstri,lstrj,lstrk,lstrl, &
                          ghondo,gpople,maxg, &
                          lpopi,lpopj,lpopk,lpopl, &
                          mini,maxi,minj,maxj,mink,maxk,minl,maxl)

                     CALL dirfck_rhf_gpu &
                         (ia,da,fa,ii,jj,kk,ll, &
                          ghondo,lstri,lstrj,lstrk,lstrl, &
                          HFSCAL,CSCALT,CUTOFFAO,nint)

                  ENDDO !MM
!$omp end parallel do

               ENDDO !JTYPE
            ENDDO !ICP
         ENDDO !ISA
      ENDIF !DO_raxintsp

      ! ===============================================================
      ! END ===== rotated axis SP ================================= END
      ! ===============================================================



     !SPD AND ERIC HERE
      CALL RHF_NCOUNT_else_gpu &
          (NCOUNT_raxintspd,NCOUNT_eric_ts, NCOUNT_rysint,&
           dirscf,schwrz,cutoff,nschwz,nsh2,l1,typscf,xints,dsh,nint,ia)

      !which integral packages to be used
      DO_raxintspd=.FALSE.
      DO_eric_ts=.FALSE.
      DO_rysint=.FALSE.

      IF(NCOUNT_raxintspd.GT.0) DO_raxintspd=.TRUE.
      IF(NCOUNT_eric_ts.GT.0) DO_eric_ts=.TRUE.
      IF(NCOUNT_rysint.GT.0) DO_rysint=.TRUE.
      !allocate shell index arrays
      IF(DO_raxintspd) ALLOCATE(INDEX_raxintspd(4,NCOUNT_raxintspd))
      IF(DO_eric_ts) ALLOCATE(INDEX_eric_ts(4,NCOUNT_eric_ts))

      CALL RHF_INDEX_new_gpu &
          (INDEX_raxintspd,INDEX_eric_ts, &
           NCOUNT_raxintspd,NCOUNT_eric_ts, &
           dirscf,schwrz,cutoff,nschwz,nsh2,l1,typscf,xints,dsh,nint,ia)

      ! SPD inegrals
      IF(DO_raxintspd) THEN
         DO ICOUNT_raxintspd=1,NCOUNT_raxintspd
            ii=INDEX_raxintspd(1,ICOUNT_raxintspd)
            jj=INDEX_raxintspd(2,ICOUNT_raxintspd)
            kk=INDEX_raxintspd(3,ICOUNT_raxintspd)
            ll=INDEX_raxintspd(4,ICOUNT_raxintspd)
            CALL ompmod_raxintspd(ii, jj, kk, ll, ghondo)

            CALL dirfck_gpu(typscf,ia,da,fa,db,fb,ghondo, &
                        l2a,nint,nflmat)
         ENDDO
      ENDIF


      ! ERIC integrals
      IF(DO_eric_ts) THEN
         DO ICOUNT_eric_ts=1,NCOUNT_eric_ts
            ii=INDEX_eric_ts(1,ICOUNT_eric_ts)
            jj=INDEX_eric_ts(2,ICOUNT_eric_ts)
            kk=INDEX_eric_ts(3,ICOUNT_eric_ts)
            ll=INDEX_eric_ts(4,ICOUNT_eric_ts)
            CALL eric_ts(ii, jj, kk, ll, ghondo)

            CALL dirfck_gpu(typscf,ia,da,fa,db,fb,ghondo, &
                        l2a,nint,nflmat)
         ENDDO
      ENDIF

  133 CONTINUE


      ! number of quartets for SP, SPD, ERIC and RYS
      CALL CONTRACTION_PROD_gpu (ACP_rys,NCP_rys)
      !RYS QUADRATURE HERE
DO_rysint=.TRUE.
IF(DO_rysint) THEN
      ALLOCATE(N_RYS_2(NCP_rys))
      ALLOCATE(IDX_RYS_2(4,NUM_QUART,NCP_rys))

      ALLOCATE(N_RYS_2000(NCP_rys))
      ALLOCATE(IDX_RYS_2000(4,NUM_QUART,NCP_rys))
      ALLOCATE(N_RYS_2100(NCP_rys))
      ALLOCATE(IDX_RYS_2100(4,NUM_QUART,NCP_rys))

      ALLOCATE(N_RYS_3(NCP_rys))
      ALLOCATE(IDX_RYS_3(4,NUM_QUART,NCP_rys))
      ALLOCATE(N_RYS_4(NCP_rys))
      ALLOCATE(IDX_RYS_4(4,NUM_QUART,NCP_rys))
      ALLOCATE(N_RYS_5(NCP_rys))
      ALLOCATE(IDX_RYS_5(4,NUM_QUART,NCP_rys))
ENDIF

IF(DO_rysint) THEN
      ! write(*,*) "original"
      CALL RHF_INDEX_rysint_gpu &
          (ACP_rys,NCP_rys, &
          dirscf,schwrz,cutoff,nschwz,nsh2, &
          l1,typscf,xints,dsh,nint,ia,&
          N_RYS_2000,N_RYS_2100,&
          IDX_RYS_2000,IDX_RYS_2100,&
          N_RYS_2,IDX_RYS_2, &
          N_RYS_3,IDX_RYS_3, &
          N_RYS_4,IDX_RYS_4, &
          N_RYS_5,IDX_RYS_5)
      ! write(*,*) "new"
      ! CALL RHF_INDEX_rysint_dddd_gpu &
      !     (ACP_rys,NCP_rys, &
      !     dirscf,schwrz,cutoff,nschwz,nsh2, &
      !     l1,typscf,xints,dsh,nint,ia,&
      !     N_RYS_2000,N_RYS_2100,&
      !     IDX_RYS_2000,IDX_RYS_2100,&
      !     N_RYS_2,IDX_RYS_2, &
      !     N_RYS_3,IDX_RYS_3, &
      !     N_RYS_4,IDX_RYS_4, &
      !     N_RYS_5,IDX_RYS_5)

ENDIF

IF(DO_rysint) THEN

!$omp target enter data &
!$omp map(TO:FA) &
!$omp map(TO:DA,IA) &
!$omp map(TO:ktype_tmp,kstart_tmp,katom_tmp,kng_tmp,kloc_tmp,kmin_tmp,kmax_tmp) &
!$omp map(TO:ex_tmp,cs_tmp,cp_tmp,cd_tmp) &
!$omp map(TO:IX,IY,IZ,JX,JY,JZ) &
!$omp map(TO:KX,KY,KZ,LX,LY,LZ) &
!$omp map(to:IDX_RYS_2,N_RYS_2,NCP_rys)
tim0 = omp_get_wtime()
      CALL rysint_gpu_2(ii,jj,kk,ll,ktype_tmp,ghondo,ddij,ngth,c,&
                     kstart_tmp,katom_tmp,kng_tmp,kloc_tmp,kmin_tmp,kmax_tmp,&
                     ex_tmp,cs_tmp,cp_tmp,cd_tmp,&
                     lstri,lstrj,lstrk,lstrl,&
                     nij_rys,ij_rys,kl_rys,ijkl_rys, &
                     qq4,&
                     AX,AY,AZ,BX,BY,BZ,RAB,CX,CY,CZ,DX,DY,DZ,RCD1,&
                     dfttyp,&
                     ia,da,fa,nint,&
                     l1,l2a,nflmat,&
                     IDX_RYS_2,N_RYS_2,NCP_rys)
tim1 = omp_get_wtime()
write(*,*) "time in rysint_gpu_2 is", tim1-tim0
!$omp target update from(fa)

!$omp target exit data &
!$omp map(delete:fa) &
!$omp map(delete:da,ia) &
!$omp map(delete:ktype_tmp,kstart_tmp,katom_tmp,kng_tmp,kloc_tmp,kmin_tmp,kmax_tmp) &
!$omp map(delete:ex_tmp,cs_tmp,cp_tmp,cd_tmp) &
!$omp map(delete:IX,IY,IZ,JX,JY,JZ) &
!$omp map(delete:KX,KY,KZ,LX,LY,LZ) &
!$omp map(delete:IDX_RYS_2,N_RYS_2,NCP_rys)


! ! !$omp target enter data &
! ! !$omp map(TO:FA) &
! ! !$omp map(TO:DA,IA) &
! ! !$omp map(TO:ktype_tmp,kstart_tmp,katom_tmp,kng_tmp,kloc_tmp,kmin_tmp,kmax_tmp) &
! ! !$omp map(TO:ex_tmp,cs_tmp,cp_tmp,cd_tmp) &
! ! !$omp map(TO:IX,IY,IZ,JX,JY,JZ) &
! ! !$omp map(TO:KX,KY,KZ,LX,LY,LZ) &
! ! !$omp map(to:IDX_RYS_2000,N_RYS_2000,IDX_RYS_2100,N_RYS_2100,NCP_rys)
! ! tim0 = omp_get_wtime()
! !       CALL rysint_gpu_2000(ii,jj,kk,ll,ktype_tmp,ghondo,ddij,ngth,c,&
! !                      kstart_tmp,katom_tmp,kng_tmp,kloc_tmp,kmin_tmp,kmax_tmp,&
! !                      ex_tmp,cs_tmp,cp_tmp,cd_tmp,&
! !                      lstri,lstrj,lstrk,lstrl,&
! !                      nij_rys,ij_rys,kl_rys,ijkl_rys, &
! !                      qq4,&
! !                      AX,AY,AZ,BX,BY,BZ,RAB,CX,CY,CZ,DX,DY,DZ,RCD1,&
! !                      dfttyp,&
! !                      ia,da,fa,nint,&
! !                      l1,l2a,nflmat,&
! !                      IDX_RYS_2000,N_RYS_2000,NCP_rys)
! !       ! CALL rysint_gpu_0200(ii,jj,kk,ll,ktype_tmp,ghondo,ddij,ngth,c,&
! !       !                kstart_tmp,katom_tmp,kng_tmp,kloc_tmp,kmin_tmp,kmax_tmp,&
! !       !                ex_tmp,cs_tmp,cp_tmp,cd_tmp,&
! !       !                lstri,lstrj,lstrk,lstrl,&
! !       !                nij_rys,ij_rys,kl_rys,ijkl_rys, &
! !       !                qq4,&
! !       !                AX,AY,AZ,BX,BY,BZ,RAB,CX,CY,CZ,DX,DY,DZ,RCD1,&
! !       !                dfttyp,&
! !       !                ia,da,fa,nint,&
! !       !                l1,l2a,nflmat,&
! !       !                IDX_RYS_0200,N_RYS_0200,NCP_rys)
! !       ! CALL rysint_gpu_0020(ii,jj,kk,ll,ktype_tmp,ghondo,ddij,ngth,c,&
! !       !                kstart_tmp,katom_tmp,kng_tmp,kloc_tmp,kmin_tmp,kmax_tmp,&
! !       !                ex_tmp,cs_tmp,cp_tmp,cd_tmp,&
! !       !                lstri,lstrj,lstrk,lstrl,&
! !       !                nij_rys,ij_rys,kl_rys,ijkl_rys, &
! !       !                qq4,&
! !       !                AX,AY,AZ,BX,BY,BZ,RAB,CX,CY,CZ,DX,DY,DZ,RCD1,&
! !       !                dfttyp,&
! !       !                ia,da,fa,nint,&
! !       !                l1,l2a,nflmat,&
! !       !                IDX_RYS_0020,N_RYS_0020,NCP_rys)
! !       ! CALL rysint_gpu_0002(ii,jj,kk,ll,ktype_tmp,ghondo,ddij,ngth,c,&
! !       !                kstart_tmp,katom_tmp,kng_tmp,kloc_tmp,kmin_tmp,kmax_tmp,&
! !       !                ex_tmp,cs_tmp,cp_tmp,cd_tmp,&
! !       !                lstri,lstrj,lstrk,lstrl,&
! !       !                nij_rys,ij_rys,kl_rys,ijkl_rys, &
! !       !                qq4,&
! !       !                AX,AY,AZ,BX,BY,BZ,RAB,CX,CY,CZ,DX,DY,DZ,RCD1,&
! !       !                dfttyp,&
! !       !                ia,da,fa,nint,&
! !       !                l1,l2a,nflmat,&
! !       !                IDX_RYS_0002,N_RYS_0002,NCP_rys)
! !       ! CALL rysint_gpu_2100(ii,jj,kk,ll,ktype_tmp,ghondo,ddij,ngth,c,&
! !       !                kstart_tmp,katom_tmp,kng_tmp,kloc_tmp,kmin_tmp,kmax_tmp,&
! !       !                ex_tmp,cs_tmp,cp_tmp,cd_tmp,&
! !       !                lstri,lstrj,lstrk,lstrl,&
! !       !                nij_rys,ij_rys,kl_rys,ijkl_rys, &
! !       !                qq4,&
! !       !                AX,AY,AZ,BX,BY,BZ,RAB,CX,CY,CZ,DX,DY,DZ,RCD1,&
! !       !                dfttyp,&
! !       !                ia,da,fa,nint,&
! !       !                l1,l2a,nflmat,&
! !       !                IDX_RYS_2100,N_RYS_2100,NCP_rys)
! ! tim1 = omp_get_wtime()
! ! write(*,*) "time in rysint_gpu_2 is", tim1-tim0
! ! !$omp target update from(fa)

! ! !$omp target exit data &
! ! !$omp map(delete:fa) &
! ! !$omp map(delete:da,ia) &
! ! !$omp map(delete:ktype_tmp,kstart_tmp,katom_tmp,kng_tmp,kloc_tmp,kmin_tmp,kmax_tmp) &
! ! !$omp map(delete:ex_tmp,cs_tmp,cp_tmp,cd_tmp) &
! ! !$omp map(delete:IX,IY,IZ,JX,JY,JZ) &
! ! !$omp map(delete:KX,KY,KZ,LX,LY,LZ) &
! ! !$omp map(delete:IDX_RYS_2000,N_RYS_2000,IDX_RYS_2100,N_RYS_2100,NCP_rys)
      ! !write(*,*) "calling rysint gpu 3"
      CALL rysint_gpu_3(ii,jj,kk,ll,ktype_tmp,ghondo,ddij,ngth,c,&
                     kstart_tmp,katom_tmp,kng_tmp,kloc_tmp,kmin_tmp,kmax_tmp,&
                     ex_tmp,cs_tmp,cp_tmp,cd_tmp,cf_tmp,cg_tmp,ch_tmp,ci_tmp,&
                     lstri,lstrj,lstrk,lstrl,&
                     nij_rys,ij_rys,kl_rys,ijkl_rys, &
                     qq4,&
                     AX,AY,AZ,BX,BY,BZ,RAB,CX,CY,CZ,DX,DY,DZ,RCD1,&
                     dfttyp,&
                     ia,da,fa,nint,&
                     l1,l2a,nflmat,&
                     IDX_RYS_3,N_RYS_3,NCP_rys, &
                     XIN,YIN,ZIN)

      !write(*,*) "calling rysint gpu 4"
      CALL rysint_gpu_4(ii,jj,kk,ll,ktype_tmp,ghondo,ddij,ngth,c,&
                     kstart_tmp,katom_tmp,kng_tmp,kloc_tmp,kmin_tmp,kmax_tmp,&
                     ex_tmp,cs_tmp,cp_tmp,cd_tmp,cf_tmp,cg_tmp,ch_tmp,ci_tmp,&
                     lstri,lstrj,lstrk,lstrl,&
                     nij_rys,ij_rys,kl_rys,ijkl_rys, &
                     qq4,&
                     AX,AY,AZ,BX,BY,BZ,RAB,CX,CY,CZ,DX,DY,DZ,RCD1,&
                     dfttyp,&
                     ia,da,fa,nint,&
                     l1,l2a,nflmat,&
                     IDX_RYS_4,N_RYS_4,NCP_rys)
!$omp target enter data &
!$omp map(to:IDX_RYS_5,N_RYS_5,NCP_rys)
tim2222s = omp_get_wtime()
write(*,*) "calling rys5"
      !write(*,*) "calling rysint gpu 5"
      CALL rysint_gpu_5(ii,jj,kk,ll,ghondo,ddij,ngth,c,&
                     kstart_tmp,katom_tmp,kloc_tmp,kmin_tmp,kmax_tmp,&
                     ex_tmp,cs_tmp,cp_tmp,cd_tmp,&
                     lstri,lstrj,lstrk,lstrl,&
                     nij_rys,ij_rys,kl_rys,ijkl_rys, &
                     qq4,&
                     AX,AY,AZ,BX,BY,BZ,RAB,CX,CY,CZ,DX,DY,DZ,RCD1,&
                     dfttyp,&
                     ia,da,fa,nint,&
                     l1,l2a,nflmat,&
                     IDX_RYS_5,N_RYS_5,NCP_rys)
tim2222e = omp_get_wtime()
write(*,*) "time in rysint_gpu_5 is", tim2222e-tim2222s
!$omp target update from(fa)

!$omp target exit data &
!$omp map(delete:IDX_RYS_5,N_RYS_5,NCP_rys)

ENDIF

      ! deallocate temporary memory for integrals
      DEALLOCATE(ddij,ghondo)

      ! allocate SP index array      
      IF(MAXSP.GT.0) THEN
         DEALLOCATE(SPTYP)
         DEALLOCATE(SPSHINDEX)
      ENDIF
            ! DEallocate shell index arrays
      IF(DO_raxintspd) DEALLOCATE(INDEX_raxintspd)
      IF(DO_eric_ts) DEALLOCATE(INDEX_eric_ts)


      END SUBROUTINE ompmod_twoei_jk

      SUBROUTINE RHF_INDEX_rysint_gpu &
               (ACP_rys,NCP_rys, &
                dirscf,schwrz,cutoff,nschwz,nsh2, &
                l1,typscf,xints,dsh,nint,ia,&
                N_RYS_2000,N_RYS_2100,&
                IDX_RYS_2000,IDX_RYS_2100,&
                N_RYS_2,IDX_RYS_2, &
                N_RYS_3,IDX_RYS_3, &
                N_RYS_4,IDX_RYS_4, &
                N_RYS_5,IDX_RYS_5)

      use mx_limits, only: mxgtot,mxsh
      USE prec, ONLY: fp

      COMMON /NSHEL / ex(mxgtot),cs(mxgtot),cp(mxgtot),cd(mxgtot),     &
                      cf(mxgtot),cg(mxgtot),ch(mxgtot),ci(mxgtot),     &
                      kstart(mxsh),katom(mxsh),ktype(mxsh),kng(mxsh),  &
                      kloc(mxsh),kmin(mxsh),kmax(mxsh),nshell

      INTEGER :: kstart,katom,ktype,kng,kloc,kmin,kmax,nshell
      REAL(KIND=fp) :: ex,cs,cp,cd,cf,cg,ch,ci

!   typscf is an 8-byte hollerith constant
      REAL(KIND=8), INTENT(IN) :: typscf
      REAL(KIND=fp), INTENT(IN) :: cutoff,xints(nsh2), dsh(nsh2)
      LOGICAL, INTENT(IN) :: schwrz
      LOGICAL :: schskp, dirscf
      LOGICAL :: lshell
      INTEGER :: ii, jj, kk, ll, ijij, klkl, jork
      INTEGER :: shlmxang, shlsumang, shltotcon   
      INTEGER, INTENT(IN) :: l1, nsh2,ia(l1)
      INTEGER, INTENT(OUT) :: nint, nschwz
      REAL(KIND=fp) :: denmax, test


      INTEGER :: NSA_rys,ISA_rys !sum angular momentum
      INTEGER :: NCP_rys,ICP_rys,ACP_rys(100) !product contraction
      INTEGER :: LAI,LAJ,LAK,LAL
      INTEGER :: LAT,LBT,LCT,LDT

      INTEGER,PARAMETER :: NUM_QUART=9169597 !number of quartets
      INTEGER :: N_RYS_2(NCP_rys)
      INTEGER :: IDX_RYS_2(4,NUM_QUART,NCP_rys)

      INTEGER :: N_RYS_2000(NCP_rys)
      INTEGER :: IDX_RYS_2000(4,NUM_QUART,NCP_rys)
      INTEGER :: N_RYS_2100(NCP_rys)
      INTEGER :: IDX_RYS_2100(4,NUM_QUART,NCP_rys)

      INTEGER :: N_RYS_3(NCP_rys)
      INTEGER :: IDX_RYS_3(4,NUM_QUART,NCP_rys)
      INTEGER :: N_RYS_4(NCP_rys)
      INTEGER :: IDX_RYS_4(4,NUM_QUART,NCP_rys)
      INTEGER :: N_RYS_5(NCP_rys)
      INTEGER :: IDX_RYS_5(4,NUM_QUART,NCP_rys)

      INTEGER :: NROOTS


      ! INTEGER :: MAXRYS        
      ! INTEGER,INTENT(OUT) :: RYSINDEX(4,MAXRYS,NCP_rys,9)
      ! INTEGER :: RYSTYP(NCP_rys,9)

      !RYSTYP=0
      N_RYS_2000=0
      N_RYS_2100=0
      N_RYS_2=0
      N_RYS_3=0
      N_RYS_4=0
      N_RYS_5=0

      DO ii = 1,nshell
         DO jj = 1,ii
            DO kk = 1,ii
               IF (ii.NE.kk) THEN
                   jork=kk
               ELSE
                   jork=jj
               ENDIF
               DO ll = 1,jork

                  ! Schwartz screening:
                  IF (schwrz) THEN
                     ijij = (ii*ii-ii)/2 + jj
                     klkl = (kk*kk-kk)/2 + ll
                     test = xints(ijij)*xints(klkl)
                     IF (dirscf) THEN
                        denmax = schwdn(dsh,ii,jj,kk,ll,ia)
                        test = test*denmax
                     ENDIF
                     schskp = test.LT.cutoff
                  ENDIF

                  IF(.NOT.schskp) THEN

                     LAI=ktype(ii)
                     LAJ=ktype(jj)
                     LAK=ktype(kk)
                     LAL=ktype(ll)
                  !i do not want to do s and p functions for rys quadrature
                  !IF(LAI.GT.2 .AND. LAJ.GT.2 .AND. LAK.GT.2 .AND. LAL.GT.2) THEN
                     shlmxang  = max(LAI,LAJ,LAK,LAL) ! max(ktype(ii), ktype(jj), ktype(kk), ktype(ll))
                     !write(*,*) "shlmxang", shlmxang
                     IF (shlmxang.GT.2) THEN
                     !write(*,*) "shlmxang in", shlmxang
                     shlsumang = LAI+LAJ+LAK+LAL-3 ! ktype(ii) + ktype(jj) + ktype(kk) + ktype(ll) - 3                     
                     shltotcon = kng(ii)*kng(jj)*kng(kk)*kng(ll)
                     !write(*,*) "shltotcon in", shltotcon
                     LAT= LAI-1
                     LBT= LAJ-1
                     LCT= LAK-1
                     LDT= LAL-1
                     NROOTS = (LAT+LBT+LCT+LDT)/2 +1
                     !NROOTS = (LAI+LAJ+LAK+LAT-2)/2
                     !write(*,*) "nroots in index is", NROOTS 

                         IF (shltotcon.EQ.1) THEN
                          !new stuff
                          IF (NROOTS.EQ.2) THEN
                            DO ICP_rys=1,NCP_rys
                               IF(ACP_rys(ICP_rys).EQ.shltotcon) THEN
                                  N_RYS_2(ICP_rys)=N_RYS_2(ICP_rys)+1 
                                  IDX_RYS_2(1,N_RYS_2(ICP_rys),ICP_rys) = ii
                                  IDX_RYS_2(2,N_RYS_2(ICP_rys),ICP_rys) = jj
                                  IDX_RYS_2(3,N_RYS_2(ICP_rys),ICP_rys) = kk
                                  IDX_RYS_2(4,N_RYS_2(ICP_rys),ICP_rys) = ll
                                  EXIT
                               ENDIF
                            ENDDO
                          
                            ! IF (lat.EQ.2 .AND. lbt.EQ.0 .AND. lct.EQ.0 .AND. ldt.EQ.0) THEN
                            !   !&
                            !   !lat.EQ.0 .AND. lbt.EQ.2 .AND. lct.EQ.0 .AND. ldt.EQ.0 .OR. &
                            !   !lat.EQ.0 .AND. lbt.EQ.0 .AND. lct.EQ.2 .AND. ldt.EQ.0 .OR. &
                            !   !lat.EQ.0 .AND. lbt.EQ.0 .AND. lct.EQ.0 .AND. ldt.EQ.2) THEN
                            !   DO ICP_rys=1,NCP_rys
                            !      IF(ACP_rys(ICP_rys).EQ.shltotcon) THEN
                            !         N_RYS_2000(ICP_rys)=N_RYS_2000(ICP_rys)+1 
                            !         IDX_RYS_2000(1,N_RYS_2000(ICP_rys),ICP_rys) = ii
                            !         IDX_RYS_2000(2,N_RYS_2000(ICP_rys),ICP_rys) = jj
                            !         IDX_RYS_2000(3,N_RYS_2000(ICP_rys),ICP_rys) = kk
                            !         IDX_RYS_2000(4,N_RYS_2000(ICP_rys),ICP_rys) = ll
                            !         EXIT
                            !      ENDIF
                            !   ENDDO

                            ! ELSE IF (lat + lbt + lct + ldt .EQ. 3 .AND. & 
                            !          lat.EQ.2 .OR. lbt.EQ.2 .OR. lct.EQ.2 .OR. ldt.EQ.2) THEN  
                            !   DO ICP_rys=1,NCP_rys
                            !      IF(ACP_rys(ICP_rys).EQ.shltotcon) THEN
                            !         N_RYS_2100(ICP_rys)=N_RYS_2100(ICP_rys)+1 
                            !         IDX_RYS_2100(1,N_RYS_2100(ICP_rys),ICP_rys) = ii
                            !         IDX_RYS_2100(2,N_RYS_2100(ICP_rys),ICP_rys) = jj
                            !         IDX_RYS_2100(3,N_RYS_2100(ICP_rys),ICP_rys) = kk
                            !         IDX_RYS_2100(4,N_RYS_2100(ICP_rys),ICP_rys) = ll
                            !         EXIT
                            !      ENDIF
                            !   ENDDO
                            ! ENDIF

                          ELSEIF (NROOTS.EQ.3) THEN
                            DO ICP_rys=1,NCP_rys
                               IF(ACP_rys(ICP_rys).EQ.shltotcon) THEN
                                  N_RYS_3(ICP_rys)=N_RYS_3(ICP_rys)+1
                                  !write(*,*) "n rys3 in index", N_RYS_3(ICP_rys) 
                                  IDX_RYS_3(1,N_RYS_3(ICP_rys),ICP_rys) = ii
                                  IDX_RYS_3(2,N_RYS_3(ICP_rys),ICP_rys) = jj
                                  IDX_RYS_3(3,N_RYS_3(ICP_rys),ICP_rys) = kk
                                  IDX_RYS_3(4,N_RYS_3(ICP_rys),ICP_rys) = ll
                                  EXIT
                               ENDIF
                            ENDDO
                          ELSEIF (NROOTS.EQ.4) THEN
                            DO ICP_rys=1,NCP_rys
                               IF(ACP_rys(ICP_rys).EQ.shltotcon) THEN
                                  N_RYS_4(ICP_rys)=N_RYS_4(ICP_rys)+1
                                  !write(*,*) "n rys4 in index", N_RYS_4(ICP_rys) 
                                  IDX_RYS_4(1,N_RYS_4(ICP_rys),ICP_rys) = ii
                                  IDX_RYS_4(2,N_RYS_4(ICP_rys),ICP_rys) = jj
                                  IDX_RYS_4(3,N_RYS_4(ICP_rys),ICP_rys) = kk
                                  IDX_RYS_4(4,N_RYS_4(ICP_rys),ICP_rys) = ll
                                  EXIT
                               ENDIF
                            ENDDO
                          ELSEIF (NROOTS.EQ.5) THEN
                            DO ICP_rys=1,NCP_rys
                               IF(ACP_rys(ICP_rys).EQ.shltotcon) THEN
                                  N_RYS_5(ICP_rys)=N_RYS_5(ICP_rys)+1
                                  !write(*,*) "n rys5 in index", N_RYS_5(ICP_rys) 
                                  IDX_RYS_5(1,N_RYS_5(ICP_rys),ICP_rys) = ii
                                  IDX_RYS_5(2,N_RYS_5(ICP_rys),ICP_rys) = jj
                                  IDX_RYS_5(3,N_RYS_5(ICP_rys),ICP_rys) = kk
                                  IDX_RYS_5(4,N_RYS_5(ICP_rys),ICP_rys) = ll
                                  EXIT
                               ENDIF
                            ENDDO
                          ENDIF
                         ENDIF
                     ENDIF
                  !ENDIF
                  ELSE
                     nschwz = nschwz + 1
                  ENDIF
                  
               ENDDO !ll
            ENDDO !kk
         ENDDO !jj
      ENDDO !ii
      END SUBROUTINE RHF_INDEX_rysint_gpu

      SUBROUTINE RHF_INDEX_rysint_dddd_gpu &
               (ACP_rys,NCP_rys, &
                dirscf,schwrz,cutoff,nschwz,nsh2, &
                l1,typscf,xints,dsh,nint,ia,&
                N_RYS_2000,N_RYS_2100,&
                IDX_RYS_2000,IDX_RYS_2100,&
                N_RYS_2,IDX_RYS_2, &
                N_RYS_3,IDX_RYS_3, &
                N_RYS_4,IDX_RYS_4, &
                N_RYS_5,IDX_RYS_5)

      use mx_limits, only: mxgtot,mxsh
      USE prec, ONLY: fp

      COMMON /NSHEL / ex(mxgtot),cs(mxgtot),cp(mxgtot),cd(mxgtot),     &
                      cf(mxgtot),cg(mxgtot),ch(mxgtot),ci(mxgtot),     &
                      kstart(mxsh),katom(mxsh),ktype(mxsh),kng(mxsh),  &
                      kloc(mxsh),kmin(mxsh),kmax(mxsh),nshell

      INTEGER :: kstart,katom,ktype,kng,kloc,kmin,kmax,nshell
      REAL(KIND=fp) :: ex,cs,cp,cd,cf,cg,ch,ci

!   typscf is an 8-byte hollerith constant
      REAL(KIND=8), INTENT(IN) :: typscf
      REAL(KIND=fp), INTENT(IN) :: cutoff,xints(nsh2), dsh(nsh2)
      LOGICAL, INTENT(IN) :: schwrz
      LOGICAL :: schskp, dirscf
      LOGICAL :: lshell
      INTEGER :: ii, jj, kk, ll, ijij, klkl, jork
      INTEGER :: shlmxang, shlsumang, shltotcon   
      INTEGER, INTENT(IN) :: l1, nsh2,ia(l1)
      INTEGER, INTENT(OUT) :: nint, nschwz
      REAL(KIND=fp) :: denmax, test


      INTEGER :: NSA_rys,ISA_rys !sum angular momentum
      INTEGER :: NCP_rys,ICP_rys,ACP_rys(100) !product contraction
      INTEGER :: LAI,LAJ,LAK,LAL
      INTEGER :: LAT,LBT,LCT,LDT

      INTEGER,PARAMETER :: NUM_QUART=9169597 !number of quartets
      INTEGER :: N_RYS_2(NCP_rys)
      INTEGER :: IDX_RYS_2(4,NUM_QUART,NCP_rys)

      INTEGER :: N_RYS_2000(NCP_rys)
      INTEGER :: IDX_RYS_2000(4,NUM_QUART,NCP_rys)
      INTEGER :: N_RYS_2100(NCP_rys)
      INTEGER :: IDX_RYS_2100(4,NUM_QUART,NCP_rys)

      INTEGER :: N_RYS_3(NCP_rys)
      INTEGER :: IDX_RYS_3(4,NUM_QUART,NCP_rys)
      INTEGER :: N_RYS_4(NCP_rys)
      INTEGER :: IDX_RYS_4(4,NUM_QUART,NCP_rys)
      INTEGER :: N_RYS_5(NCP_rys)
      INTEGER :: IDX_RYS_5(4,NUM_QUART,NCP_rys)

      INTEGER :: NROOTS

      !RYSTYP=0
      N_RYS_2000=0
      N_RYS_2100=0
      N_RYS_2=0
      N_RYS_3=0
      N_RYS_4=0
      N_RYS_5=0

      DO ii = 1,nshell
         LAI=ktype(ii)
         IF(LAI.NE.3) CYCLE
         DO jj = 1,ii
            LAJ=ktype(jj)
            IF(LAJ.NE.3) CYCLE            
            DO kk = 1,ii
               LAK=ktype(kk)
               IF(LAK.NE.3) CYCLE               
               IF (ii.NE.kk) THEN
                   jork=kk
               ELSE
                   jork=jj
               ENDIF
               DO ll = 1,jork
                  LAL=ktype(ll)
                  IF(LAL.NE.3) CYCLE
                  ! Schwartz screening:
                  IF (schwrz) THEN
                     ijij = (ii*ii-ii)/2 + jj
                     klkl = (kk*kk-kk)/2 + ll
                     test = xints(ijij)*xints(klkl)
                     IF (dirscf) THEN
                        denmax = schwdn(dsh,ii,jj,kk,ll,ia)
                        test = test*denmax
                     ENDIF
                     schskp = test.LT.cutoff
                  ENDIF

                  IF(.NOT.schskp) THEN

                     ! LAI=ktype(ii)
                     ! LAJ=ktype(jj)
                     ! LAK=ktype(kk)
                     ! LAL=ktype(ll)
                  !i do not want to do s and p functions for rys quadrature
                  !IF(LAI.GT.2 .AND. LAJ.GT.2 .AND. LAK.GT.2 .AND. LAL.GT.2) THEN
                     shlmxang  = max(LAI,LAJ,LAK,LAL) ! max(ktype(ii), ktype(jj), ktype(kk), ktype(ll))
                     !write(*,*) "shlmxang", shlmxang
                     IF (shlmxang.GT.2) THEN
                     !write(*,*) "shlmxang in", shlmxang
                     shlsumang = LAI+LAJ+LAK+LAL-3 ! ktype(ii) + ktype(jj) + ktype(kk) + ktype(ll) - 3                     
                     shltotcon = kng(ii)*kng(jj)*kng(kk)*kng(ll)
                     !write(*,*) "shltotcon in", shltotcon
                     LAT= LAI-1
                     LBT= LAJ-1
                     LCT= LAK-1
                     LDT= LAL-1
                     !write(*,*) "LAT",LAT
                     NROOTS = (LAT+LBT+LCT+LDT)/2 +1
                     !NROOTS = (LAI+LAJ+LAK+LAT-2)/2
                     !write(*,*) "nroots in index is", NROOTS 

                         IF (shltotcon.EQ.1) THEN
                          !new stuff
                          IF (NROOTS.EQ.2) THEN
                            !write(*,*) "nroot 2"
                            DO ICP_rys=1,NCP_rys
                               IF(ACP_rys(ICP_rys).EQ.shltotcon) THEN
                                  N_RYS_2(ICP_rys)=N_RYS_2(ICP_rys)+1 
                                  IDX_RYS_2(1,N_RYS_2(ICP_rys),ICP_rys) = ii
                                  IDX_RYS_2(2,N_RYS_2(ICP_rys),ICP_rys) = jj
                                  IDX_RYS_2(3,N_RYS_2(ICP_rys),ICP_rys) = kk
                                  IDX_RYS_2(4,N_RYS_2(ICP_rys),ICP_rys) = ll
                                  EXIT
                               ENDIF
                            ENDDO
                          
                            ! IF (lat.EQ.2 .AND. lbt.EQ.0 .AND. lct.EQ.0 .AND. ldt.EQ.0) THEN
                            !   !&
                            !   !lat.EQ.0 .AND. lbt.EQ.2 .AND. lct.EQ.0 .AND. ldt.EQ.0 .OR. &
                            !   !lat.EQ.0 .AND. lbt.EQ.0 .AND. lct.EQ.2 .AND. ldt.EQ.0 .OR. &
                            !   !lat.EQ.0 .AND. lbt.EQ.0 .AND. lct.EQ.0 .AND. ldt.EQ.2) THEN
                            !   DO ICP_rys=1,NCP_rys
                            !      IF(ACP_rys(ICP_rys).EQ.shltotcon) THEN
                            !         N_RYS_2000(ICP_rys)=N_RYS_2000(ICP_rys)+1 
                            !         IDX_RYS_2000(1,N_RYS_2000(ICP_rys),ICP_rys) = ii
                            !         IDX_RYS_2000(2,N_RYS_2000(ICP_rys),ICP_rys) = jj
                            !         IDX_RYS_2000(3,N_RYS_2000(ICP_rys),ICP_rys) = kk
                            !         IDX_RYS_2000(4,N_RYS_2000(ICP_rys),ICP_rys) = ll
                            !         EXIT
                            !      ENDIF
                            !   ENDDO

                            ! ELSE IF (lat + lbt + lct + ldt .EQ. 3 .AND. & 
                            !          lat.EQ.2 .OR. lbt.EQ.2 .OR. lct.EQ.2 .OR. ldt.EQ.2) THEN  
                            !   DO ICP_rys=1,NCP_rys
                            !      IF(ACP_rys(ICP_rys).EQ.shltotcon) THEN
                            !         N_RYS_2100(ICP_rys)=N_RYS_2100(ICP_rys)+1 
                            !         IDX_RYS_2100(1,N_RYS_2100(ICP_rys),ICP_rys) = ii
                            !         IDX_RYS_2100(2,N_RYS_2100(ICP_rys),ICP_rys) = jj
                            !         IDX_RYS_2100(3,N_RYS_2100(ICP_rys),ICP_rys) = kk
                            !         IDX_RYS_2100(4,N_RYS_2100(ICP_rys),ICP_rys) = ll
                            !         EXIT
                            !      ENDIF
                            !   ENDDO
                            ! ENDIF

                          ELSEIF (NROOTS.EQ.3) THEN
                            !write(*,*) "nroot 3"
                            DO ICP_rys=1,NCP_rys
                               IF(ACP_rys(ICP_rys).EQ.shltotcon) THEN
                                  N_RYS_3(ICP_rys)=N_RYS_3(ICP_rys)+1
                                  !write(*,*) "n rys3 in index", N_RYS_3(ICP_rys) 
                                  IDX_RYS_3(1,N_RYS_3(ICP_rys),ICP_rys) = ii
                                  IDX_RYS_3(2,N_RYS_3(ICP_rys),ICP_rys) = jj
                                  IDX_RYS_3(3,N_RYS_3(ICP_rys),ICP_rys) = kk
                                  IDX_RYS_3(4,N_RYS_3(ICP_rys),ICP_rys) = ll
                                  EXIT
                               ENDIF
                            ENDDO
                          ELSEIF (NROOTS.EQ.4) THEN
                            !write(*,*) "nroot 4"
                            DO ICP_rys=1,NCP_rys
                               IF(ACP_rys(ICP_rys).EQ.shltotcon) THEN
                                  N_RYS_4(ICP_rys)=N_RYS_4(ICP_rys)+1
                                  !write(*,*) "n rys4 in index", N_RYS_4(ICP_rys) 
                                  IDX_RYS_4(1,N_RYS_4(ICP_rys),ICP_rys) = ii
                                  IDX_RYS_4(2,N_RYS_4(ICP_rys),ICP_rys) = jj
                                  IDX_RYS_4(3,N_RYS_4(ICP_rys),ICP_rys) = kk
                                  IDX_RYS_4(4,N_RYS_4(ICP_rys),ICP_rys) = ll
                                  EXIT
                               ENDIF
                            ENDDO
                          ELSEIF (NROOTS.EQ.5) THEN
                            !write(*,*) "nroot 5"
                            DO ICP_rys=1,NCP_rys
                               IF(ACP_rys(ICP_rys).EQ.shltotcon) THEN
                                  N_RYS_5(ICP_rys)=N_RYS_5(ICP_rys)+1
                                  !write(*,*) "n rys5 in index", N_RYS_5(ICP_rys) 
                                  IDX_RYS_5(1,N_RYS_5(ICP_rys),ICP_rys) = ii
                                  IDX_RYS_5(2,N_RYS_5(ICP_rys),ICP_rys) = jj
                                  IDX_RYS_5(3,N_RYS_5(ICP_rys),ICP_rys) = kk
                                  IDX_RYS_5(4,N_RYS_5(ICP_rys),ICP_rys) = ll
                                  EXIT
                               ENDIF
                            ENDDO
                          ENDIF
                         ENDIF
                     ENDIF
                  !ENDIF
                  ELSE
                     nschwz = nschwz + 1
                  ENDIF
                  
               ENDDO !ll
            ENDDO !kk
         ENDDO !jj
      ENDDO !ii
      END SUBROUTINE RHF_INDEX_rysint_dddd_gpu

      SUBROUTINE CONTRACTION_PROD_gpu (ACP,NCP)

      use mx_limits, only: mxgtot,mxsh
      USE prec, ONLY: fp

      COMMON /NSHEL / ex(mxgtot),cs(mxgtot),cp(mxgtot),cd(mxgtot),     &
                      cf(mxgtot),cg(mxgtot),ch(mxgtot),ci(mxgtot),     &
                      kstart(mxsh),katom(mxsh),ktype(mxsh),kng(mxsh),  &
                      kloc(mxsh),kmin(mxsh),kmax(mxsh),nshell

      INTEGER :: kstart,katom,ktype,kng,kloc,kmin,kmax,nshell
      REAL(KIND=fp) :: ex,cs,cp,cd,cf,cg,ch,ci


      INTEGER :: ACP(100),NCP,ICP,MM
      INTEGER :: NCONTYP(100)
      INTEGER :: ish, ii, jj, kk, ll, ijij, klkl, jork
      INTEGER :: NCT,NGI

      NCONTYP(1) = KNG(1)
      NCT=1
      DO ISH=2,NSHELL
         NGI=KNG(ISH)
         MM=0
         DO KK=1,NCT
            IF(NCONTYP(KK).EQ.NGI) EXIT
            MM=MM+1
         ENDDO
         IF(MM.EQ.NCT) THEN
            NCT=NCT+1
            NCONTYP(NCT)=NGI
         ENDIF
      ENDDO

      NCP=0
      DO II=1,NCT
         DO JJ=1,II
            DO KK=1,II
               IF (II.NE.KK) THEN
                   JORK=KK
               ELSE
                   JORK=JJ
               ENDIF
               DO LL = 1,JORK
                  NCP=NCP+1
                  ACP(NCP)=NCONTYP(II)*NCONTYP(JJ)*NCONTYP(KK)*NCONTYP(LL)
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      END SUBROUTINE CONTRACTION_PROD_gpu


     SUBROUTINE RHF_NCOUNT_rysint_gpu &
               (MAXRYS,ACP_rys,NCP_rys,dirscf,schwrz,cutoff,nschwz, &
                nsh2,l1,typscf,xints,dsh,nint,ia)


      COMMON /NSHEL / ex(mxgtot),cs(mxgtot),cp(mxgtot),cd(mxgtot),     &
                      cf(mxgtot),cg(mxgtot),ch(mxgtot),ci(mxgtot),     &
                      kstart(mxsh),katom(mxsh),ktype(mxsh),kng(mxsh),  &
                      kloc(mxsh),kmin(mxsh),kmax(mxsh),nshell
      INTEGER :: kstart,katom,ktype,kng,kloc,kmin,kmax,nshell
      REAL(KIND=fp) :: ex,cs,cp,cd,cf,cg,ch,ci

!   typscf is an 8-byte hollerith constant
      REAL(KIND=8), INTENT(IN) :: typscf
      REAL(KIND=fp), INTENT(IN) :: cutoff,xints(nsh2), dsh(nsh2)
      LOGICAL, INTENT(IN) :: schwrz
      LOGICAL :: schskp, dirscf
      LOGICAL :: lshell
      INTEGER :: ii, jj, kk, ll, ijij, klkl, jork
      INTEGER :: shlmxang, shlsumang, shltotcon   
      INTEGER, INTENT(IN) :: l1, nsh2,ia(l1)
      INTEGER, INTENT(OUT) :: nint, nschwz
      REAL(KIND=fp) :: denmax, test

      INTEGER :: NSA_rys,ISA_rys !sum angular momentum
      INTEGER :: NCP_rys,ICP_rys,ACP_rys(100) !product contraction
      INTEGER :: NSP0_rys,NSP1_rys,NSP2_rys,NSP3_rys,NSP4_rys
      iNTEGER :: NSP5_rys,NSP6_rys,NSP7_rys,NSP8_rys,NSP9_rys,NSP10_rys
      !INTEGER :: NCOUNT_rysint   
      INTEGER :: MAXRYS 
      INTEGER :: MM_rys    

      !MAXRYS=0
      DO ii = 1,nshell
         DO jj = 1,ii
            DO kk = 1,ii
               IF (ii.NE.kk) THEN
                   jork=kk
               ELSE
                   jork=jj
               ENDIF
               DO ll = 1,jork

                  ! Schwartz screening:
                  IF (schwrz) THEN
                     ijij = (ii*ii-ii)/2 + jj
                     klkl = (kk*kk-kk)/2 + ll
                     test = xints(ijij)*xints(klkl)
                     IF (dirscf) THEN
                        denmax = schwdn(dsh,ii,jj,kk,ll,ia)
                        test = test*denmax
                     ENDIF
                     schskp = test.LT.cutoff
                  ENDIF

                  IF(.NOT.schskp) THEN

                     shlmxang  = max(ktype(ii), ktype(jj), ktype(kk), ktype(ll))
                     shltotcon = kng(ii)*kng(jj)*kng(kk)*kng(ll)                     
                     IF (shlmxang.GT.2) THEN
                        ACP_rys(1)=shltotcon
                        GOTO 150
                     ENDIF
                  ELSE
                     nschwz = nschwz + 1
                  ENDIF
                  
               ENDDO !ll
            ENDDO !kk
         ENDDO !jj
      ENDDO !ii


  150 CONTINUE

      NSP0_rys=0
      NSP1_rys=0
      NSP2_rys=0
      NSP3_rys=0
      NSP4_rys=0
      NSP5_rys=0
      NSP6_rys=0
      NSP7_rys=0
      NSP8_rys=0
      NSP9_rys=0
      NSP10_rys=0
      NCP_rys=1
      DO ii = 1,nshell
         DO jj = 1,ii
            DO kk = 1,ii
               IF (ii.NE.kk) THEN
                   jork=kk
               ELSE
                   jork=jj
               ENDIF

               DO ll = 1,jork

                  ! Schwartz screening:
                  IF (schwrz) THEN
                     ijij = (ii*ii-ii)/2 + jj
                     klkl = (kk*kk-kk)/2 + ll
                     test = xints(ijij)*xints(klkl)
                     IF (dirscf) THEN
                        denmax = schwdn(dsh,ii,jj,kk,ll,ia)
                        test = test*denmax
                     ENDIF
                     schskp = test.LT.cutoff
                  ENDIF

                  IF(.NOT.schskp) THEN

                     !write(*,*) "ktype(ii)", ktype(ii)
                     !write(*,*) "ktype(jj)", ktype(jj)
                     !write(*,*) "ktype(kk)", ktype(kk)
                     !write(*,*) "ktype(ll)", ktype(ll)
                     shlmxang  = max(ktype(ii), ktype(jj), ktype(kk), ktype(ll)) 
                     shlsumang = ktype(ii) + ktype(jj) + ktype(kk) + ktype(ll) - 4
                     !write(*,*) "shlsumang is", shlsumang
                     shltotcon = kng(ii)*kng(jj)*kng(kk)*kng(ll)
                     !write(*,*) "shltotcon is", shltotcon
                     lshell = (kmax(ii)-kmin(ii)).EQ.3 .OR. &
                              (kmax(jj)-kmin(jj)).EQ.3 .OR. &
                              (kmax(kk)-kmin(kk)).EQ.3 .OR. &
                              (kmax(ll)-kmin(ll)).EQ.3

                     
                     IF (shlmxang.GT.2) THEN

                            IF(shlsumang.EQ.0) THEN
                              NSP0_rys=NSP0_rys+1
                            ELSEIF(shlsumang.EQ.1) THEN
                              NSP1_rys=NSP1_rys+1
                            ELSEIF(shlsumang.EQ.2) THEN
                              NSP2_rys=NSP2_rys+1
                            ELSEIF(shlsumang.EQ.3) THEN
                              NSP3_rys=NSP3_rys+1
                            ELSEIF(shlsumang.EQ.4) THEN
                              NSP4_rys=NSP4_rys+1
                            ELSEIF(shlsumang.EQ.5) THEN
                              NSP5_rys=NSP5_rys+1
                            ELSEIF(shlsumang.EQ.6) THEN
                              NSP6_rys=NSP6_rys+1
                            ELSEIF(shlsumang.EQ.7) THEN
                              NSP7_rys=NSP7_rys+1
                            ELSEIF(shlsumang.EQ.8) THEN
                              NSP8_rys=NSP8_rys+1
                              !write(*,*) "NSP8_rys ", NSP8_rys
                            ELSEIF(shlsumang.EQ.9) THEN
                              NSP9_rys=NSP9_rys+1
                            ENDIF
                            MM_rys=0
                            DO ICP_rys=1,NCP_rys
                              MM_rys=MM_rys+1
                              IF(ACP_rys(ICP_rys).EQ.shltotcon) EXIT
                            ENDDO
                            IF(MM_rys.EQ.NCP_rys) THEN
                              NCP_rys=NCP_rys+1
                              ACP_rys(NCP_rys)=shltotcon
                            ENDIF

                     ENDIF
                  ELSE
                     nschwz = nschwz + 1
                  ENDIF
                  
               ENDDO !ll
            ENDDO !kk
         ENDDO !jj
      ENDDO !ii

      MAXRYS=MAX(NSP0_rys,NSP1_rys,NSP2_rys,NSP3_rys,NSP4_rys, &
                NSP5_rys,NSP6_rys,NSP7_rys,NSP8_rys,NSP9_rys)
      !write(*,*) "maxrys is", MAXRYS

      END SUBROUTINE RHF_NCOUNT_rysint_gpu



      SUBROUTINE RHF_NCOUNT_else_gpu &
               (NCOUNT_raxintspd,NCOUNT_eric_ts, NCOUNT_rysint,&
                dirscf,schwrz,cutoff,nschwz,nsh2,l1,typscf,xints,dsh,nint,ia)


      COMMON /NSHEL / ex(mxgtot),cs(mxgtot),cp(mxgtot),cd(mxgtot),     &
                      cf(mxgtot),cg(mxgtot),ch(mxgtot),ci(mxgtot),     &
                      kstart(mxsh),katom(mxsh),ktype(mxsh),kng(mxsh),  &
                      kloc(mxsh),kmin(mxsh),kmax(mxsh),nshell
      INTEGER :: kstart,katom,ktype,kng,kloc,kmin,kmax,nshell
      REAL(KIND=fp) :: ex,cs,cp,cd,cf,cg,ch,ci

!   typscf is an 8-byte hollerith constant
      REAL(KIND=8), INTENT(IN) :: typscf
      REAL(KIND=fp), INTENT(IN) :: cutoff,xints(nsh2), dsh(nsh2)
      LOGICAL, INTENT(IN) :: schwrz
      LOGICAL :: schskp, dirscf
      LOGICAL :: lshell
      INTEGER :: ii, jj, kk, ll, ijij, klkl, jork
      INTEGER :: shlmxang, shlsumang, shltotcon   
      INTEGER, INTENT(IN) :: l1, nsh2,ia(l1)
      INTEGER, INTENT(OUT) :: nint, nschwz
      REAL(KIND=fp) :: denmax, test

      INTEGER :: NCOUNT_raxintsp,NCOUNT_raxintspd,NCOUNT_eric_ts,NCOUNT_rysint        



      NCOUNT_raxintspd=0
      NCOUNT_eric_ts=0
      NCOUNT_rysint=0

      DO ii = 1,nshell
         DO jj = 1,ii

            DO kk = 1,ii

               IF (ii.NE.kk) THEN
                   jork=kk
               ELSE
                   jork=jj
               ENDIF

               DO ll = 1,jork

                  ! Schwartz screening:
                  IF (schwrz) THEN
                     ijij = (ii*ii-ii)/2 + jj
                     klkl = (kk*kk-kk)/2 + ll
                     test = xints(ijij)*xints(klkl)
                     IF (dirscf) THEN
                        denmax = schwdn(dsh,ii,jj,kk,ll,ia)
                        test = test*denmax
                     ENDIF
                     schskp = test.LT.cutoff
                  ENDIF

                  IF(.NOT.schskp) THEN

                     shlmxang  = max(ktype(ii), ktype(jj), ktype(kk), ktype(ll))
                     shlsumang = ktype(ii) + ktype(jj) + ktype(kk) + ktype(ll) - 4
                     shltotcon = kng(ii)*kng(jj)*kng(kk)*kng(ll)
                     lshell = (kmax(ii)-kmin(ii)).EQ.3 .OR. &
                              (kmax(jj)-kmin(jj)).EQ.3 .OR. &
                              (kmax(kk)-kmin(kk)).EQ.3 .OR. &
                              (kmax(ll)-kmin(ll)).EQ.3

                    
                     IF (shlmxang.LE.2) THEN
                     ELSE
                         IF (shltotcon.NE.1) THEN
                             IF (shlmxang.EQ.3) THEN
                                 NCOUNT_raxintspd=NCOUNT_raxintspd+1
                             ELSE IF (shlsumang.LE.5 .AND. .NOT.lshell.AND.shlmxang.LE.5) THEN
                                 NCOUNT_eric_ts=NCOUNT_eric_ts+1
                             ELSE 
                                 NCOUNT_rysint=NCOUNT_rysint+1
                             ENDIF
                          ELSE
                             NCOUNT_rysint=NCOUNT_rysint+1
                          ENDIF
                     ENDIF
                  ELSE
                     nschwz = nschwz + 1
                  ENDIF
                  
               ENDDO !ll
            ENDDO !kk
         ENDDO !jj
      ENDDO !ii

      END SUBROUTINE RHF_NCOUNT_else_gpu

     SUBROUTINE RHF_INDEX_new_gpu &
               (INDEX_raxintspd,INDEX_eric_ts, &
                NCOUNT_raxintspd,NCOUNT_eric_ts, &
                dirscf,schwrz,cutoff,nschwz,nsh2,l1,typscf,xints,dsh,nint,ia)

      COMMON /NSHEL / ex(mxgtot),cs(mxgtot),cp(mxgtot),cd(mxgtot),     &
                      cf(mxgtot),cg(mxgtot),ch(mxgtot),ci(mxgtot),     &
                      kstart(mxsh),katom(mxsh),ktype(mxsh),kng(mxsh),  &
                      kloc(mxsh),kmin(mxsh),kmax(mxsh),nshell
      INTEGER :: kstart,katom,ktype,kng,kloc,kmin,kmax,nshell
      REAL(KIND=fp) :: ex,cs,cp,cd,cf,cg,ch,ci

!   typscf is an 8-byte hollerith constant
      REAL(KIND=8), INTENT(IN) :: typscf
      REAL(KIND=fp), INTENT(IN) :: cutoff,xints(nsh2), dsh(nsh2)
      LOGICAL, INTENT(IN) :: schwrz
      LOGICAL :: schskp, dirscf
      LOGICAL :: lshell
      INTEGER :: ii, jj, kk, ll, ijij, klkl, jork
      INTEGER :: shlmxang, shlsumang, shltotcon   
      INTEGER, INTENT(IN) :: l1, nsh2,ia(l1)
      INTEGER, INTENT(OUT) :: nint, nschwz
      REAL(KIND=fp) :: denmax, test

      INTEGER :: NCOUNT_raxintsp,NCOUNT_raxintspd,NCOUNT_eric_ts      
      INTEGER :: INDEX_raxintspd(4,NCOUNT_raxintspd), &
                 INDEX_eric_ts(4,NCOUNT_eric_ts)
      INTEGER :: ICOUNT_raxintspd,ICOUNT_eric_ts


      ICOUNT_raxintspd=0
      ICOUNT_eric_ts=0
      !ICOUNT_rysint=0

      DO ii = 1,nshell
         DO jj = 1,ii

            DO kk = 1,ii

               IF (ii.NE.kk) THEN
                   jork=kk
               ELSE
                   jork=jj
               ENDIF

               DO ll = 1,jork

                  ! Schwartz screening:
                  IF (schwrz) THEN
                     ijij = (ii*ii-ii)/2 + jj
                     klkl = (kk*kk-kk)/2 + ll
                     test = xints(ijij)*xints(klkl)
                     IF (dirscf) THEN
                        denmax = schwdn(dsh,ii,jj,kk,ll,ia)
                        test = test*denmax
                     ENDIF
                     schskp = test.LT.cutoff
                  ENDIF

                  IF(.NOT.schskp) THEN

                     shlmxang  = max(ktype(ii), ktype(jj), ktype(kk), ktype(ll))
                     shlsumang = ktype(ii) + ktype(jj) + ktype(kk) + ktype(ll) - 4
                     shltotcon = kng(ii)*kng(jj)*kng(kk)*kng(ll)
                     lshell = (kmax(ii)-kmin(ii)).EQ.3 .OR. &
                              (kmax(jj)-kmin(jj)).EQ.3 .OR. &
                              (kmax(kk)-kmin(kk)).EQ.3 .OR. &
                              (kmax(ll)-kmin(ll)).EQ.3
                    
                     IF (shlmxang.LE.2) THEN
                     ELSE
                         IF (shltotcon.NE.1) THEN
                             IF (shlmxang.EQ.3) THEN
                                 !write(*,*) "ktype(ii)",ktype(ii)
                                 !write(*,*) "ktype(jj)",ktype(jj) 
                                 !write(*,*) "ktype(kk)",ktype(kk) 
                                 !write(*,*) "ktype(ll)",ktype(ll) 

                                 ICOUNT_raxintspd=ICOUNT_raxintspd+1
                                 INDEX_raxintspd(1,ICOUNT_raxintspd)=ii
                                 INDEX_raxintspd(2,ICOUNT_raxintspd)=jj
                                 INDEX_raxintspd(3,ICOUNT_raxintspd)=kk
                                 INDEX_raxintspd(4,ICOUNT_raxintspd)=ll
                              ELSE IF (shlsumang.LE.5 .AND. .NOT.lshell.AND.shlmxang.LE.5) THEN
                                 ICOUNT_eric_ts=ICOUNT_eric_ts+1
                                 INDEX_eric_ts(1,ICOUNT_eric_ts)=ii
                                 INDEX_eric_ts(2,ICOUNT_eric_ts)=jj
                                 INDEX_eric_ts(3,ICOUNT_eric_ts)=kk
                                 INDEX_eric_ts(4,ICOUNT_eric_ts)=ll
                             ENDIF
                         ELSE
                         ENDIF
                     ENDIF
                  ELSE
                     nschwz = nschwz + 1
                  ENDIF
                  
               ENDDO !ll
            ENDDO !kk
         ENDDO !jj
      ENDDO !ii
      END SUBROUTINE RHF_INDEX_new_gpu



      SUBROUTINE RHF_NCOUNT_gpu &
               (NCOUNT_raxintspd,NCOUNT_eric_ts,NCOUNT_rysint, &
                dirscf,schwrz,cutoff,nschwz,nsh2,l1,typscf,xints,dsh,nint,ia)


      COMMON /NSHEL / ex(mxgtot),cs(mxgtot),cp(mxgtot),cd(mxgtot),     &
                      cf(mxgtot),cg(mxgtot),ch(mxgtot),ci(mxgtot),     &
                      kstart(mxsh),katom(mxsh),ktype(mxsh),kng(mxsh),  &
                      kloc(mxsh),kmin(mxsh),kmax(mxsh),nshell
      INTEGER :: kstart,katom,ktype,kng,kloc,kmin,kmax,nshell
      REAL(KIND=fp) :: ex,cs,cp,cd,cf,cg,ch,ci

!   typscf is an 8-byte hollerith constant
      REAL(KIND=8), INTENT(IN) :: typscf
      REAL(KIND=fp), INTENT(IN) :: cutoff,xints(nsh2), dsh(nsh2)
      LOGICAL, INTENT(IN) :: schwrz
      LOGICAL :: schskp, dirscf
      LOGICAL :: lshell
      INTEGER :: ii, jj, kk, ll, ijij, klkl, jork
      INTEGER :: shlmxang, shlsumang, shltotcon   
      INTEGER, INTENT(IN) :: l1, nsh2,ia(l1)
      INTEGER, INTENT(OUT) :: nint, nschwz
      REAL(KIND=fp) :: denmax, test

      INTEGER :: NCOUNT_raxintsp,NCOUNT_raxintspd,NCOUNT_eric_ts,NCOUNT_rysint        



      NCOUNT_raxintspd=0
      NCOUNT_eric_ts=0
      NCOUNT_rysint=0

      DO ii = 1,nshell
         DO jj = 1,ii

            DO kk = 1,ii

               IF (ii.NE.kk) THEN
                   jork=kk
               ELSE
                   jork=jj
               ENDIF

               DO ll = 1,jork

                  ! Schwartz screening:
                  IF (schwrz) THEN
                     ijij = (ii*ii-ii)/2 + jj
                     klkl = (kk*kk-kk)/2 + ll
                     test = xints(ijij)*xints(klkl)
                     IF (dirscf) THEN
                        denmax = schwdn(dsh,ii,jj,kk,ll,ia)
                        test = test*denmax
                     ENDIF
                     schskp = test.LT.cutoff
                  ENDIF

                  IF(.NOT.schskp) THEN

                     shlmxang  = max(ktype(ii), ktype(jj), ktype(kk), ktype(ll))
                     shlsumang = ktype(ii) + ktype(jj) + ktype(kk) + ktype(ll) - 4
                     shltotcon = kng(ii)*kng(jj)*kng(kk)*kng(ll)
                     lshell = (kmax(ii)-kmin(ii)).EQ.3 .OR. &
                              (kmax(jj)-kmin(jj)).EQ.3 .OR. &
                              (kmax(kk)-kmin(kk)).EQ.3 .OR. &
                              (kmax(ll)-kmin(ll)).EQ.3

                     
                     IF (shlmxang.LE.2) THEN
                     ELSE
                         IF (shltotcon.NE.1) THEN
                             IF (shlmxang.EQ.3) THEN
                                 NCOUNT_raxintspd=NCOUNT_raxintspd+1
                             ELSE IF (shlsumang.LE.5 .AND. .NOT.lshell.AND.shlmxang.LE.5) THEN
                                 NCOUNT_eric_ts=NCOUNT_eric_ts+1
                             ELSE
                                 NCOUNT_rysint=NCOUNT_rysint+1
                             ENDIF
                         ELSE
                              NCOUNT_rysint=NCOUNT_rysint+1
                         ENDIF
                     ENDIF
                  ELSE
                     nschwz = nschwz + 1
                  ENDIF
                  
               ENDDO !ll
            ENDDO !kk
         ENDDO !jj
      ENDDO !ii

      END SUBROUTINE RHF_NCOUNT_gpu



      SUBROUTINE RHF_NCOUNT_raxintsp_gpu &
               (MAXSP,ACP,NCP,dirscf,schwrz,cutoff,nschwz, &
                nsh2,l1,typscf,xints,dsh,nint,ia)


      COMMON /NSHEL / ex(mxgtot),cs(mxgtot),cp(mxgtot),cd(mxgtot),     &
                      cf(mxgtot),cg(mxgtot),ch(mxgtot),ci(mxgtot),     &
                      kstart(mxsh),katom(mxsh),ktype(mxsh),kng(mxsh),  &
                      kloc(mxsh),kmin(mxsh),kmax(mxsh),nshell
      INTEGER :: kstart,katom,ktype,kng,kloc,kmin,kmax,nshell
      REAL(KIND=fp) :: ex,cs,cp,cd,cf,cg,ch,ci

!   typscf is an 8-byte hollerith constant
      REAL(KIND=8), INTENT(IN) :: typscf
      REAL(KIND=fp), INTENT(IN) :: cutoff,xints(nsh2), dsh(nsh2)
      LOGICAL, INTENT(IN) :: schwrz
      LOGICAL :: schskp, dirscf
      LOGICAL :: lshell
      INTEGER :: ii, jj, kk, ll, ijij, klkl, jork
      INTEGER :: shlmxang, shlsumang, shltotcon   
      INTEGER, INTENT(IN) :: l1, nsh2,ia(l1)
      INTEGER, INTENT(OUT) :: nint, nschwz
      REAL(KIND=fp) :: denmax, test

      INTEGER :: ACP(100),NCP,ICP,MM
      INTEGER :: MAXSP      
      INTEGER :: NSP0,NSP1,NSP2,NSP3,NSP4



      ! initialization ACP
      DO ii = 1,nshell
         IF(ktype(ii).GT.2) CYCLE
         DO jj = 1,ii
            IF(ktype(jj).GT.2) CYCLE
            DO kk = 1,ii
               IF(ktype(kk).GT.2) CYCLE
               IF (ii.NE.kk) THEN
                   jork=kk
               ELSE
                   jork=jj
               ENDIF
               DO ll = 1,jork
                  IF(ktype(ll).GT.2) CYCLE

                  ! Schwartz screening:
                  IF (schwrz) THEN
                     ijij = (ii*ii-ii)/2 + jj
                     klkl = (kk*kk-kk)/2 + ll
                     test = xints(ijij)*xints(klkl)
                     IF (dirscf) THEN
                        denmax = schwdn(dsh,ii,jj,kk,ll,ia)
                        test = test*denmax
                     ENDIF
                     schskp = test.LT.cutoff
                  ENDIF

                  IF(.NOT.schskp) THEN

                     shlmxang  = max(ktype(ii), ktype(jj), ktype(kk), ktype(ll))
                     shltotcon = kng(ii)*kng(jj)*kng(kk)*kng(ll)                     
                     IF (shlmxang.LE.2) THEN
                        ACP(1)=shltotcon
                        GOTO 100
                     ENDIF
                  ELSE
                     nschwz = nschwz + 1
                  ENDIF
                  
               ENDDO !ll
            ENDDO !kk
         ENDDO !jj
      ENDDO !ii


  100 CONTINUE
      NSP0=0
      NSP1=0
      NSP2=0
      NSP3=0
      NSP4=0

      NCP=1    

      DO ii = 1,nshell
         IF(ktype(ii).GT.2) CYCLE
         DO jj = 1,ii
            IF(ktype(jj).GT.2) CYCLE
            DO kk = 1,ii
               IF(ktype(kk).GT.2) CYCLE
               IF (ii.NE.kk) THEN
                   jork=kk
               ELSE
                   jork=jj
               ENDIF
               DO ll = 1,jork
                  IF(ktype(ll).GT.2) CYCLE

                  ! Schwartz screening:
                  IF (schwrz) THEN
                     ijij = (ii*ii-ii)/2 + jj
                     klkl = (kk*kk-kk)/2 + ll
                     test = xints(ijij)*xints(klkl)
                     IF (dirscf) THEN
                        denmax = schwdn(dsh,ii,jj,kk,ll,ia)
                        test = test*denmax
                     ENDIF
                     schskp = test.LT.cutoff
                  ENDIF

                  IF(.NOT.schskp) THEN

                     shlmxang  = max(ktype(ii), ktype(jj), ktype(kk), ktype(ll))
                     shlsumang = ktype(ii) + ktype(jj) + ktype(kk) + ktype(ll) - 4
                     shltotcon = kng(ii)*kng(jj)*kng(kk)*kng(ll)
                     lshell = (kmax(ii)-kmin(ii)).EQ.3 .OR. &
                              (kmax(jj)-kmin(jj)).EQ.3 .OR. &
                              (kmax(kk)-kmin(kk)).EQ.3 .OR. &
                              (kmax(ll)-kmin(ll)).EQ.3
                    
                     IF (shlmxang.LE.2) THEN
                        IF(shlsumang.EQ.0) THEN
                          NSP0=NSP0+1
                        ELSEIF(shlsumang.EQ.1) THEN
                          NSP1=NSP1+1
                        ELSEIF(shlsumang.EQ.2) THEN
                          NSP2=NSP2+1
                        ELSEIF(shlsumang.EQ.3) THEN
                          NSP3=NSP3+1
                        ELSEIF(shlsumang.EQ.4) THEN
                          NSP4=NSP4+1
                        ENDIF
                        MM=0
                        DO ICP=1,NCP
                           MM=MM+1
                           IF(ACP(ICP).EQ.shltotcon) EXIT
                        ENDDO
                        IF(MM.EQ.NCP) THEN
                           NCP=NCP+1
                           ACP(NCP)=shltotcon
                        ENDIF
                     ENDIF
                  ELSE
                     nschwz = nschwz + 1
                  ENDIF
                  
               ENDDO !ll
            ENDDO !kk
         ENDDO !jj
      ENDDO !ii


      MAXSP=MAX(NSP0,NSP1,NSP2,NSP3,NSP4)


      END SUBROUTINE RHF_NCOUNT_raxintsp_gpu




      SUBROUTINE RHF_INDEX_gpu &
               (INDEX_raxintspd,INDEX_eric_ts,INDEX_rysint, &
                NCOUNT_raxintspd,NCOUNT_eric_ts,NCOUNT_rysint, &
                dirscf,schwrz,cutoff,nschwz,nsh2,l1,typscf,xints,dsh,nint,ia)

      COMMON /NSHEL / ex(mxgtot),cs(mxgtot),cp(mxgtot),cd(mxgtot),     &
                      cf(mxgtot),cg(mxgtot),ch(mxgtot),ci(mxgtot),     &
                      kstart(mxsh),katom(mxsh),ktype(mxsh),kng(mxsh),  &
                      kloc(mxsh),kmin(mxsh),kmax(mxsh),nshell
      INTEGER :: kstart,katom,ktype,kng,kloc,kmin,kmax,nshell
      REAL(KIND=fp) :: ex,cs,cp,cd,cf,cg,ch,ci

!   typscf is an 8-byte hollerith constant
      REAL(KIND=8), INTENT(IN) :: typscf
      REAL(KIND=fp), INTENT(IN) :: cutoff,xints(nsh2), dsh(nsh2)
      LOGICAL, INTENT(IN) :: schwrz
      LOGICAL :: schskp, dirscf
      LOGICAL :: lshell
      INTEGER :: ii, jj, kk, ll, ijij, klkl, jork
      INTEGER :: shlmxang, shlsumang, shltotcon   
      INTEGER, INTENT(IN) :: l1, nsh2,ia(l1)
      INTEGER, INTENT(OUT) :: nint, nschwz
      REAL(KIND=fp) :: denmax, test

      INTEGER :: NCOUNT_raxintsp,NCOUNT_raxintspd,NCOUNT_eric_ts,NCOUNT_rysint        
      INTEGER :: INDEX_raxintspd(4,NCOUNT_raxintspd), &
                 INDEX_eric_ts(4,NCOUNT_eric_ts),INDEX_rysint(4,NCOUNT_rysint)
      INTEGER :: ICOUNT_raxintspd,ICOUNT_eric_ts,ICOUNT_rysint


      ICOUNT_raxintspd=0
      ICOUNT_eric_ts=0
      ICOUNT_rysint=0

      DO ii = 1,nshell
         DO jj = 1,ii

            DO kk = 1,ii

               IF (ii.NE.kk) THEN
                   jork=kk
               ELSE
                   jork=jj
               ENDIF

               DO ll = 1,jork

                  ! Schwartz screening:
                  IF (schwrz) THEN
                     ijij = (ii*ii-ii)/2 + jj
                     klkl = (kk*kk-kk)/2 + ll
                     test = xints(ijij)*xints(klkl)
                     IF (dirscf) THEN
                        denmax = schwdn(dsh,ii,jj,kk,ll,ia)
                        test = test*denmax
                     ENDIF
                     schskp = test.LT.cutoff
                  ENDIF

                  IF(.NOT.schskp) THEN

                     shlmxang  = max(ktype(ii), ktype(jj), ktype(kk), ktype(ll))
                     shlsumang = ktype(ii) + ktype(jj) + ktype(kk) + ktype(ll) - 4
                     shltotcon = kng(ii)*kng(jj)*kng(kk)*kng(ll)
                     lshell = (kmax(ii)-kmin(ii)).EQ.3 .OR. &
                              (kmax(jj)-kmin(jj)).EQ.3 .OR. &
                              (kmax(kk)-kmin(kk)).EQ.3 .OR. &
                              (kmax(ll)-kmin(ll)).EQ.3
                    
                     IF (shlmxang.LE.2) THEN
                     ELSE
                         IF (shltotcon.NE.1) THEN
                             IF (shlmxang.EQ.3) THEN
                                 ICOUNT_raxintspd=ICOUNT_raxintspd+1
                                 INDEX_raxintspd(1,ICOUNT_raxintspd)=ii
                                 INDEX_raxintspd(2,ICOUNT_raxintspd)=jj
                                 INDEX_raxintspd(3,ICOUNT_raxintspd)=kk
                                 INDEX_raxintspd(4,ICOUNT_raxintspd)=ll
                              ELSE IF (shlsumang.LE.5 .AND. .NOT.lshell.AND.shlmxang.LE.5) THEN
                                 ICOUNT_eric_ts=ICOUNT_eric_ts+1
                                 INDEX_eric_ts(1,ICOUNT_eric_ts)=ii
                                 INDEX_eric_ts(2,ICOUNT_eric_ts)=jj
                                 INDEX_eric_ts(3,ICOUNT_eric_ts)=kk
                                 INDEX_eric_ts(4,ICOUNT_eric_ts)=ll
                             ELSE
                                 ICOUNT_rysint=ICOUNT_rysint+1
                                 INDEX_rysint(1,ICOUNT_rysint)=ii
                                 INDEX_rysint(2,ICOUNT_rysint)=jj
                                 INDEX_rysint(3,ICOUNT_rysint)=kk
                                 INDEX_rysint(4,ICOUNT_rysint)=ll
                             ENDIF
                         ELSE
                            ICOUNT_rysint=ICOUNT_rysint+1
                            INDEX_rysint(1,ICOUNT_rysint)=ii
                            INDEX_rysint(2,ICOUNT_rysint)=jj
                            INDEX_rysint(3,ICOUNT_rysint)=kk
                            INDEX_rysint(4,ICOUNT_rysint)=ll
                         ENDIF
                     ENDIF
                  ELSE
                     nschwz = nschwz + 1
                  ENDIF
                  
               ENDDO !ll
            ENDDO !kk
         ENDDO !jj
      ENDDO !ii
      END SUBROUTINE RHF_INDEX_gpu





      SUBROUTINE RHF_INDEX_raxintsp_gpu &
               (SPSHINDEX,ACP,NCP,SPTYP,MAXSP, &
                dirscf,schwrz,cutoff,nschwz,nsh2,l1, &
                typscf,xints,dsh,nint,ia)

      COMMON /NSHEL / ex(mxgtot),cs(mxgtot),cp(mxgtot),cd(mxgtot),     &
                      cf(mxgtot),cg(mxgtot),ch(mxgtot),ci(mxgtot),     &
                      kstart(mxsh),katom(mxsh),ktype(mxsh),kng(mxsh),  &
                      kloc(mxsh),kmin(mxsh),kmax(mxsh),nshell


      INTEGER :: kstart,katom,ktype,kng,kloc,kmin,kmax,nshell
      REAL(KIND=fp) :: ex,cs,cp,cd,cf,cg,ch,ci

      ! typscf is an 8-byte hollerith constant
      REAL(KIND=8), INTENT(IN) :: typscf
      REAL(KIND=fp), INTENT(IN) :: cutoff,xints(nsh2), dsh(nsh2)
      LOGICAL, INTENT(IN) :: schwrz
      LOGICAL :: schskp, dirscf
      LOGICAL :: lshell
      INTEGER :: ii, jj, kk, ll, ijij, klkl, jork
      INTEGER :: shlmxang, shlsumang, shltotcon   
      INTEGER, INTENT(IN) :: l1, nsh2,ia(l1)
      INTEGER, INTENT(OUT) :: nint, nschwz
      REAL(KIND=fp) :: denmax, test


      INTEGER :: MAXSP
      INTEGER :: SPTYP(6,NCP,5)
      INTEGER,INTENT(OUT) :: SPSHINDEX(4,MAXSP,6,NCP,5)
      INTEGER :: ACP(100),NCP,ICP

      INTEGER :: LAI,LAJ,LAK,LAL, LAT,LBT,LCT,LDT

      INTEGER :: ITYPE,JTYPE


      SPTYP=0

      DO ii = 1,nshell
         DO jj = 1,ii
            DO kk = 1,ii
               IF (ii.NE.kk) THEN
                   jork=kk
               ELSE
                   jork=jj
               ENDIF
               DO ll = 1,jork

                  ! Schwartz screening:
                  IF (schwrz) THEN
                     ijij = (ii*ii-ii)/2 + jj
                     klkl = (kk*kk-kk)/2 + ll
                     test = xints(ijij)*xints(klkl)
                     IF (dirscf) THEN
                        denmax = schwdn(dsh,ii,jj,kk,ll,ia)
                        test = test*denmax
                     ENDIF
                     schskp = test.LT.cutoff
                  ENDIF

                  IF(.NOT.schskp) THEN

                     ! check max angular momentum
                     LAI=ktype(ii)
                     LAJ=ktype(jj)
                     LAK=ktype(kk)
                     LAL=ktype(ll)
                     shlmxang  = max(LAI,LAJ,LAK,LAL)

                     ! SP quartets only
                     IF (shlmxang.LE.2) THEN

                        ! six type of SP quartets (JTYPE)
                        ! 0000; 0001; 0011; 0101; 0111; 1111
                        LAT= LAI-1
                        LBT= LAJ-1
                        LCT= LAK-1
                        LDT= LAL-1
                        ITYPE= 1+LDT+2*(LCT+2*(LBT+2*LAT))

                        IF(ITYPE.EQ. 1) THEN
                           JTYPE= 1
                        ELSEIF(ITYPE.EQ. 2 .OR. ITYPE.EQ. 3 .OR. &
                               ITYPE.EQ. 5 .OR. ITYPE.EQ. 9) THEN
                           JTYPE= 2
                        ELSEIF(ITYPE.EQ. 4 .OR. ITYPE.EQ.13) THEN
                           JTYPE= 3
                        ELSEIF(ITYPE.EQ. 6 .OR. ITYPE.EQ. 7 .OR. &
                               ITYPE.EQ.10 .OR. ITYPE.EQ.11) THEN
                           JTYPE= 4
                        ELSEIF(ITYPE.EQ. 8 .OR. ITYPE.EQ.12 .OR. &
                               ITYPE.EQ.14 .OR. ITYPE.EQ.15) THEN
                           JTYPE= 5
                        ELSEIF(ITYPE.EQ.16) THEN
                           JTYPE= 6
                        ENDIF

                        ! the sum of angular momentum of four shells
                        shlsumang = LAI+LAJ+LAK+LAL-3

                        ! contraction product
                        shltotcon = kng(ii)*kng(jj)*kng(kk)*kng(ll)

                        ! number of quartets that has the same:
                        ! * Total angular momentum
                        ! ** Contraction product
                        ! *** Integral type JTYPE
                        DO ICP=1,NCP
                           IF(ACP(ICP).EQ.shltotcon) THEN
                              SPTYP(JTYPE,ICP,shlsumang)=SPTYP(JTYPE,ICP,shlsumang)+1
                              EXIT
                           ENDIF
                        ENDDO

                        ! get shell index of quartets satisfied
                        ! the three conditions (*), (**) & (***)
                        SPSHINDEX(1,SPTYP(JTYPE,ICP,shlsumang),JTYPE,ICP,shlsumang)=ii
                        SPSHINDEX(2,SPTYP(JTYPE,ICP,shlsumang),JTYPE,ICP,shlsumang)=jj
                        SPSHINDEX(3,SPTYP(JTYPE,ICP,shlsumang),JTYPE,ICP,shlsumang)=kk
                        SPSHINDEX(4,SPTYP(JTYPE,ICP,shlsumang),JTYPE,ICP,shlsumang)=ll
                     ENDIF                    
                  ELSE
                     nschwz = nschwz + 1
                  ENDIF
                  
               ENDDO !ll
            ENDDO !kk
         ENDDO !jj
      ENDDO !ii

      END SUBROUTINE RHF_INDEX_raxintsp_gpu


!---------------------------------------------------------------------
!=====================================================================
!---------------------------------------------------------------------

!*MODULE OMPMOD   *DECK OMPMOD_TWOEI_KL
!
!>    @brief   Threaded version of two-electron integral
!>             calculation routine. For closed shell only.
!>             OpenMP parallelization over K and L shell indices,
!>             MPI -  over I and J indices.
!>             Related input file options:
!>             `$INTGRL INTOMP=2 $END`
!>
!>    @details Calculates two-electron contribution to the Fock
!>             matrix. Based on `TWOEI` subroutine from `int2a.src`
!>
!>    @author  Vladimir Mironov
!
!     REVISION HISTORY:
!>    @date _Jun, 2017_ Initial release
!
!     PARAMETERS:
!
!>    @param[in]     typscf    hollerith, the kind of SCF (=RHF,UHF, etc.)
!>    @param[in]     schwrz    determines whether to use integral screening
!>    @param[out]    nint      number of calculated ERIs
!>    @param[out]    nschwz    number of ERIs that were screened out
!>    @param[in]     l1        size of index array
!>    @param[in]     l2a       size of Fock and density matrices (alpha)
!>    @param[in]     l2b       size of Fock and density matrices (beta)
!>    @param[in]     xints(:)  array of exchange integrals
!>    @param[in]     nsh2      size of `xints(:)` array
!>    @param[in]     maxg      size of temporary array for Rys code
!>    @param[in]     ia(:)     index array, contains "triangular numbers"
!>    @param[in,out] da(:)     density matrix for alpha electrons
!>    @param[out]    fa(:)     Fock matrix for alpha electrons
!>    @param[in,out] db(:)     density matrix for beta electrons,
!>                                 used in open-shell calculations only
!>    @param[out]    fb(:)     Fock matrix for beta electrons,
!>                                 used in open-shell calculations only
!>    @param[in]     dsh       density matrix packed in shells
!>                                 for screening purposes
!>    @param[in]     nflmat    when >1 selects CPHF calculation
!>    @param[in]     cutoff    cutoff for integral screening
!>    @param[in]     oflag     logical parameter for debug timing output
  SUBROUTINE ompmod_twoei_kl(typscf,schwrz,nint,nschwz, &
                             l1,l2a,l2b,xints, &
                             nsh2,maxg, &
                             ia,da,fa,db,fb,dsh,nflmat, &
                             cutoff, oflag)

    LOGICAL, INTENT(IN) :: &
      schwrz

    INTEGER, INTENT(OUT) :: &
      nint, nschwz

    INTEGER, INTENT(IN) :: &
      l1, l2a, l2b, nsh2, nflmat, maxg, &
      ia(l1)

!   typscf is an 8-byte hollerith constant
    REAL(KIND=8), INTENT(IN) :: &
        typscf

    REAL(KIND=fp), INTENT(IN) :: &
      cutoff, &
      xints(nsh2), dsh(nsh2)

    REAL(KIND=fp), INTENT(INOUT) :: &
      da(l2a), db(l2b)

    REAL(KIND=fp), INTENT(OUT) :: &
!      fa(l2a), fb(l2b)
      fa(l2a*nflmat), fb(l2b*nflmat)

    LOGICAL, INTENT(IN) :: &
        oflag
!
!
!
    COMMON /NSHEL / ex(mxgtot),cs(mxgtot),cp(mxgtot),cd(mxgtot),     &
                    cf(mxgtot),cg(mxgtot),ch(mxgtot),ci(mxgtot),     &
                    kstart(mxsh),katom(mxsh),ktype(mxsh),kng(mxsh),  &
                    kloc(mxsh),kmin(mxsh),kmax(mxsh),nshell
        INTEGER :: kstart,katom,ktype,kng,kloc,kmin,kmax,nshell
        REAL(KIND=fp) :: ex,cs,cp,cd,cf,cg,ch,ci

    COMMON /par   / me,master,nproc,ibtyp,iptim,goparr,dskwrk,maswrk
        INTEGER :: me,master,nproc,ibtyp,iptim
        LOGICAL :: dskwrk, maswrk, goparr

    COMMON /iofile/ ir,iw,ip,is,ipk,idaf,nav,ioda(950)
        INTEGER :: ir,iw,ip,is,ipk,idaf,nav,ioda

    COMMON /fmcom / xx(1)
        REAL(KIND=fp), TARGET :: xx

    COMMON /shlnos/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl, &
                    mini,minj,mink,minl,maxi,maxj,maxk,maxl, &
                    nij,ij,kl,ijkl
        INTEGER :: lit,ljt,lkt,llt,loci,locj,lock,locl, &
                   mini,minj,mink,minl,maxi,maxj,maxk,maxl, &
                   nij,ij,kl,ijkl
        REAL(KIND=fp) :: qq4

    COMMON /maxc  / cmax(mxgtot),cmaxa(mxgsh),cmaxb(mxgsh), &
                    cmaxc(mxgsh),cmaxd(mxgsh),ismlp(mxg2),ismlq
        REAL(KIND=fp) :: cmax,cmaxa,cmaxb,cmaxc,cmaxd
        INTEGER :: ismlp,ismlq

    COMMON /gout  / gpople(768),norgp
        REAL(KIND=fp) :: gpople
        INTEGER :: norgp

    COMMON /nlrcf / lrint
        LOGICAL :: lrint

    COMMON /dftpar/ dfttyp(20),exena,exenb,exenc, &
                    idft34,nauxfun,nauxshl
        REAL(KIND=fp) :: dfttyp,exena,exenb,exenc
        INTEGER :: idft34,nauxfun,nauxshl

!   Common blocks /MAXC  /, /SHLNOS/ and /GOUT/ contain
!   both global (cmax(:),qq4,norgp) and thread-local
!   data (the rest). Set them threadprivate and copyin
!   global data later in $omp parallel section.

!$omp threadprivate(/shlnos/,/maxc  /,/gout  /,/nlrcf /,/dftpar/)

    INTEGER :: &
      next, num_threads, ithread, thr_nshq, &
      ii, jj, kk, ll, ijij, klkl, ijmax, klmax

    LOGICAL :: &
      schskp, dirscf

    REAL(KIND=fp) :: &
      denmax, tim, tim0, tim1, tim2, tim3, tim4, test

    REAL(KIND=fp),DIMENSION(:),ALLOCATABLE :: &
      ddij_t, ghondo_t

!
!   --- Initialization of variables ---
!

    tim = zero
    CALL tsecnd(tim)

    nint   = 0
    nschwz = 0
    schskp = .FALSE.
    denmax = zero
    dirscf = .TRUE.

!   no symmetry supported yet
    qq4 = 1
    norgp = 0


#if __TIMING==1
!   text for debut output

    IF (oflag) WRITE(iw,dbgfmt1)
#endif
!
!   --- Initiate OpenMP parallel ---
!
!   early thread invocation decrease OpenMP lib overhead
!
!$omp parallel                                  &
!$omp   private(ithread, thr_nshq,              &
!$omp     kk,ll, klkl, test,                    &
!$omp     ijij, ii, jj, ijmax, klmax,           &
!$omp     tim, tim0, tim1, tim2, tim3, tim4,    &
!$omp     ddij_t, ghondo_t                      &
!$omp   )                                       &
!$omp   firstprivate(denmax, schskp)            &
!$omp   shared(ia,da,l1,l2a,l2b,num_threads,    &
!$omp     next, maxg, xx,                       &
!$omp     cutoff, dirscf, xints, dsh, nflmat,   &
!$omp     typscf,schwrz,oflag)                  &
!$omp   reduction(+:fa,fb,nint,nschwz)          &
!$omp   copyin(cmax,qq4,norgp,lrint,/dftpar/)


!
!   --- variable initialization for OpenMP ---
!

    num_threads = omp_get_num_threads()

!   reset all info stats
    ithread  = omp_get_thread_num()
    thr_nshq = 0
#if __TIMING==1
    tim1     = zero
    tim2     = zero
    tim3     = zero
    tim4     = zero
!   main timer initialization
    tim0     = omp_get_wtime()
#endif


!   allocate temporary memory for integrals
    allocate(ddij_t(mxang2*mxg2),ghondo_t(maxg))

    ijmax = nshell*(nshell+1)/2
!   go over shells from largest value (nshell)
!   to 1 for better load balance
    ijij = ijmax + 1
ijc:DO
!   --- MPI load balance ---
        IF (ibtyp.eq.1) THEN
!           Dynamic load balance
!$omp barrier
!$omp master
            CALL ddi_dlbnext(next)
!$omp flush (next)
!$omp end master
!$omp barrier
            ijij = ijmax - next
            IF (ijij.LE.0) EXIT ijc
        ELSE
!           Static load balance
            ijij = ijij - 1
            IF (ijij.LE.0) EXIT ijc
            IF (MOD(ijij,nproc).NE.me) CYCLE ijc
        ENDIF

        ii = ceiling((sqrt(8.0*ijij+1)-1)/2)
        jj = ijij - ii*(ii-1)/2

        klmax = ii*(ii+1)/2

!   --- Begin OpenMP integral calculation ---

!$omp do schedule(dynamic,1)
klc:    DO klkl = 1, klmax

            kk = ceiling((sqrt(8.0*klkl+1)-1)/2)
            ll = klkl - kk*(kk-1)/2

            IF (ii.EQ.kk.AND.ll.GT.jj) CYCLE klc

!           Schwartz screening:
#if __TIMING==1
            tim1 = omp_get_wtime()
#endif
            IF (schwrz) THEN
                test = xints(ijij)*xints(klkl)
                IF (dirscf) THEN
                    denmax = schwdn(dsh,ii,jj,kk,ll,ia)
                    test = test*denmax
                ENDIF
                schskp = test.LT.cutoff
            ENDIF

#if __TIMING==1
            tim2 = tim2 + omp_get_wtime() - tim1
#endif

            IF(.NOT.schskp) THEN

                    CALL ompmod_shellquart(ii,jj,kk,ll,ghondo_t,ddij_t)

                    thr_nshq = thr_nshq + 1

#if __TIMING==1
                    tim3 = tim3 + omp_get_wtime() - tim1
#endif
!          ---  Fock matrix update ---

                    CALL dirfck(typscf,ia,da,fa,db,fb,ghondo_t, &
                                l2a,nint,nflmat)

#if __TIMING==1
                    tim4 = tim4 + omp_get_wtime() - tim1
#endif

            ELSE
                nschwz = nschwz + 1
            ENDIF
        ENDDO klc
!$omp enddo nowait
    ENDDO ijc
!$omp barrier

#if __TIMING==1
    tim0 = omp_get_wtime() - tim0

!   debug timing output
    IF (oflag) THEN
!$omp do ordered
        DO ii = 0, num_threads-1
!$omp ordered
            WRITE(iw,dbgfmt2) ithread,thr_nshq,tim3-tim2,tim4-tim3,tim2,tim0
!$omp end ordered
        ENDDO
!$omp enddo
    ENDIF
#endif

!   deallocate temporary memory for integrals
    deallocate(ddij_t,ghondo_t)

!$omp end parallel

    CALL ddi_dlbreset

  END SUBROUTINE ompmod_twoei_kl

!---------------------------------------------------------------------
!=====================================================================
!---------------------------------------------------------------------

!*MODULE OMPMOD   *DECK OMPMOD_TWOEI_KL_RHF_SHF
!
!>    @brief   Threaded version of two-electron integral
!>             calculation routine. Fock matrix is shared
!>             among threads. For closed shell only.
!>             Load balance of OpenMP code
!>             is done over K and L shell indices,
!>             MPI - over I and J indices.
!>             Related input file options:
!>             `$INTGRL INTOMP=2 SHFOCK=.TRUE. $END`
!>
!>    @details Calculates two-electron contribution to the Fock
!>             matrix. Based on `TWOEI` subroutine from `int2a.src`
!>
!>    @author  Vladimir Mironov
!
!     REVISION HISTORY:
!>    @date _Jun, 2017_ Initial release
!
!     PARAMETERS:
!
!>    @param[in]     typscf    hollerith, the kind of SCF (=RHF,UHF, etc.)
!>    @param[in]     schwrz    determines whether to use integral screening
!>    @param[out]    nint      number of calculated ERIs
!>    @param[out]    nschwz    number of ERIs that were screened out
!>    @param[in]     l1        size of index array
!>    @param[in]     l2a       size of Fock and density matrices (alpha)
!>    @param[in]     l2b       size of Fock and density matrices (beta)
!>    @param[in]     xints(:)  array of exchange integrals
!>    @param[in]     nsh2      size of `xints(:)` array
!>    @param[in]     maxg      size of temporary array for Rys code
!>    @param[in]     ia(:)     index array, contains "triangular numbers"
!>    @param[in,out] da(:)     density matrix for alpha electrons
!>    @param[out]    fa(:)     Fock matrix for alpha electrons
!>    @param[in]     dsh       density matrix packed in shells
!>                                 for screening purposes
!>    @param[in]     nflmat    when >1 selects CPHF calculation
!>    @param[in]     cutoff    cutoff for integral screening
!>    @param[in]     oflag     logical parameter for debug timing output
  SUBROUTINE ompmod_twoei_shf_kl_rhf(typscf,schwrz,nint,nschwz, &
                                     l1,l2a,l2b,xints, &
                                     nsh2,maxg, &
                                     ia,da,fa,dsh,nflmat, &
                                     cutoff, oflag)
    USE blkint

    USE camdft, ONLY: alphac => cam_alpha, betac => cam_beta, &
      cammu => cam_mu, camvwn => cam_vwn5, camlyp => cam_lyp, camflag
    USE lrcdft, ONLY: lcflag, emu, emu2, lrfile
    USE ompmod_tools, ONLY: &
        ompmod_tools_col_reduce, &
        ompmod_dirfck_rhf

    LOGICAL, INTENT(IN) :: &
      schwrz

    INTEGER, INTENT(OUT) :: &
      nint, nschwz

    INTEGER, INTENT(IN) :: &
      l1, l2a, l2b, nsh2, nflmat, maxg, &
      ia(l1)

!   typscf is an 8-byte hollerith constant
    REAL(KIND=8), INTENT(IN) :: &
        typscf

    REAL(KIND=fp), INTENT(INOUT) :: &
      cutoff, &
      xints(nsh2), dsh(nsh2)

    REAL(KIND=fp), INTENT(INOUT) :: &
      da(l2a)

    REAL(KIND=fp), INTENT(OUT) :: &
      fa(l2a)

    LOGICAL, INTENT(IN) :: &
        oflag
!
!
!
    COMMON /par   / me,master,nproc,ibtyp,iptim,goparr,dskwrk,maswrk
        INTEGER :: me,master,nproc,ibtyp,iptim
        LOGICAL :: dskwrk, maswrk, goparr

    COMMON /iofile/ ir,iw,ip,is,ipk,idaf,nav,ioda(950)
        INTEGER :: ir,iw,ip,is,ipk,idaf,nav,ioda

    COMMON /fmcom / xx(1)
        REAL(KIND=fp), TARGET :: xx

    COMMON /NSHEL / ex(mxgtot),cs(mxgtot),cp(mxgtot),cd(mxgtot),     &
                    cf(mxgtot),cg(mxgtot),ch(mxgtot),ci(mxgtot),     &
                    kstart(mxsh),katom(mxsh),ktype(mxsh),kng(mxsh),  &
                    kloc(mxsh),kmin(mxsh),kmax(mxsh),nshell
        INTEGER :: kstart,katom,ktype,kng,kloc,kmin,kmax,nshell
        REAL(KIND=fp) :: ex,cs,cp,cd,cf,cg,ch,ci

    COMMON /shlnos/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl, &
                    mini,minj,mink,minl,maxi,maxj,maxk,maxl, &
                    nij,ij,kl,ijkl
        INTEGER :: lit,ljt,lkt,llt,loci,locj,lock,locl, &
                   mini,minj,mink,minl,maxi,maxj,maxk,maxl, &
                   nij,ij,kl,ijkl
        REAL(KIND=fp) :: qq4

    COMMON /maxc  / cmax(mxgtot),cmaxa(mxgsh),cmaxb(mxgsh), &
                    cmaxc(mxgsh),cmaxd(mxgsh),ismlp(mxg2),ismlq
        REAL(KIND=fp) :: cmax,cmaxa,cmaxb,cmaxc,cmaxd
        INTEGER :: ismlp,ismlq

    COMMON /gout  / gpople(768),norgp
        REAL(KIND=fp) :: gpople
        INTEGER :: norgp

    COMMON /shlexc/ norgsh(3),norgsp(3),iexch,nangm,ngth(4)
        INTEGER :: norgsh,norgsp,iexch,nangm,ngth

    COMMON /nlrcf / lrint
        LOGICAL :: lrint

    COMMON /dftpar/ dfttyp(20),exena,exenb,exenc, &
                    idft34,nauxfun,nauxshl
        REAL(KIND=fp) :: dfttyp,exena,exenb,exenc
        INTEGER :: idft34,nauxfun,nauxshl

!   Common blocks /MAXC  /, /SHLNOS/ and /GOUT/ contain
!   both global (cmax(:),qq4,norgp) and thread-local
!   data (the rest). Set them threadprivate and copyin
!   global data later in $omp parallel section.

!$omp threadprivate(/shlnos/,/maxc  /,/gout  /,/nlrcf /,/dftpar/)

    INTEGER :: &
      next, max_threads, ithread, nthreads, thr_nshq, &
      ii, jj, kk, ll, ijij, klkl, ijmax, klmax, &
      iprev

    LOGICAL :: &
      schskp, dirscf

    REAL(KIND=fp) :: &
      denmax, tim, tim0, tim1, tim2, tim3, tim4, test, &
      hfscal, cscalt

    REAL(KIND=fp), ALLOCATABLE :: &
      ddij_t(:), ghondo_t(:), &
      faomp(:,:), daomp(:,:), &
      Fi(:,:,:), Fj(:,:,:)

!dir$ attributes align : 64 :: faomp, daomp, Fi, Fj

!
!   --- Initialization of variables ---
!

    tim = zero
    CALL tsecnd(tim)

    nint   = 0
    nschwz = 0
    schskp = .FALSE.
    denmax = zero
    dirscf = .TRUE.

!   DFT scaling factors
    hfscal=dfttyp(3)
    cscalt=1.0d+00
    if(lcflag) then
        if(lrint) then
            hfscal=1.0d+00
            cscalt=0.0d+00
        else
            hfscal=0.0d+00
            cscalt=1.0d+00
        endif
    endif
    if(camflag.and.lrint) cscalt=0.0d+00

!   no symmetry supported yet
    qq4 = 1
    norgp = 0

!   use square fock and density matrices
    allocate( faomp(l1,l1), &
              daomp(l1,l1))

    faomp(:,:) = 0

    CALL blkint_shlsrt

    CALL blkint_mtrx_unpack(da,daomp,l1)
    CALL blkint_aomtrx_fwd(daomp)
    CALL blkint_shmtrx_fwd(xints)
    CALL blkint_shmtrx_fwd(dsh)

    max_threads = omp_get_max_threads()

!   allocate temporary storage for I and J columns of
!   the Fock matrix
    allocate(Fi(l1,nangm,0:max_threads-1))
    allocate(Fj(l1,nangm,0:max_threads-1))


#if __TIMING==1
!   text for timing output
    IF (oflag) WRITE(iw,dbgfmt1)
#endif
!
!   --- Initiate OpenMP parallel ---
!
!   early thread invocation decrease OpenMP lib overhead
!
!$omp parallel                                  &
!$omp   private(ithread, thr_nshq,              &
!$omp     kk,ll, klkl, test,                    &
!$omp     ijij, ii, jj, ijmax, klmax,           &
!$omp     tim, tim0, tim1, tim2, tim3, tim4,    &
!$omp     ddij_t, ghondo_t                      &
!$omp   )                                       &
!$omp   firstprivate(denmax, schskp)            &
!$omp   shared(ia,da,l1,l2a,l2b,max_threads,    &
!$omp     next, maxg, xx, nthreads,             &
!$omp     cutoff, dirscf, xints, dsh, nflmat,   &
!$omp     typscf,schwrz,oflag, iprev,           &
!$omp     daomp, Fi, Fj, faomp                  &
!$omp   )                                       &
!$omp   reduction(+:nint,nschwz)                &
!$omp   copyin(cmax,qq4,norgp,lrint,/dftpar/)
!!$omp   reduction(+:faomp)                      &
!!$omp   reduction(+:fa,fb)                      &


!
!   --- Variable initialization for OpenMP ---
!

!   allocate temporary memory for integrals
    allocate(ddij_t(mxang2*mxg2),ghondo_t(maxg))

!  Reset all info stats
    ithread  = omp_get_thread_num()
    nthreads = omp_get_num_threads()
    thr_nshq = 0
#if __TIMING==1
    tim1     = zero
    tim2     = zero
    tim3     = zero
    tim4     = zero
!   main timer initialization
    tim0     = omp_get_wtime()
#endif

    Fi(:,:,ithread) = 0.0_fp
    Fj(:,:,ithread) = 0.0_fp

    ijmax = nshell*(nshell+1)/2
!   go over shells from largest value (nshell)
!   to 1 for better load balance
    ijij = ijmax + 1
    iprev = nshell+1
ijc:DO
!   --- MPI load balance ---
        IF (ibtyp.eq.1) THEN
!           Dynamic load balance
!$omp barrier
!$omp master
            CALL ddi_dlbnext(next)
!$omp flush (next)
!$omp end master
!$omp barrier
            ijij = ijmax - next
            IF (ijij.LE.0) EXIT ijc
        ELSE
!           static load balance
            ijij = ijij - 1
            IF (ijij.LE.0) EXIT ijc
            IF (MOD(ijij,nproc).NE.me) CYCLE ijc
        ENDIF

        ii = ceiling((sqrt(8.0*ijij+1)-1)/2)
        jj = ijij - ii*(ii-1)/2

        klmax = ii*(ii+1)/2

        IF (schwrz) THEN
            test = xints(ijij)
            IF (dirscf) THEN
                denmax = schwdn(dsh,ii,jj,ii,jj,ia)
                test = test*denmax
            ENDIF
            schskp = test.LT.cutoff
        ENDIF
        IF (schskp) cycle ijc

!       if we switch to new I subblock, flush Fi to the Fock matrix
        IF (iprev/=ii.AND.iprev<=nshell) THEN
            CALL ompmod_tools_col_reduce(kmax(iprev)-kmin(iprev)+1, &
                                         kloc(iprev),ithread,nthreads,Fi,faomp)
        ENDIF
!   --- Begin OpenMP integral calculation ---

!$omp do schedule(dynamic,1)
klc:    DO klkl = 1, klmax

            kk = ceiling((sqrt(8.0*klkl+1)-1)/2)
            ll = klkl - kk*(kk-1)/2

            IF (ii.EQ.kk.AND.ll.GT.jj) CYCLE klc

!           Schwartz screening:
#if __TIMING==1
            tim1 = omp_get_wtime()
#endif
            IF (schwrz) THEN
                test = xints(ijij)*xints(klkl)
                IF (dirscf) THEN
                    denmax = schwdn(dsh,ii,jj,kk,ll,ia)
                    test = test*denmax
                ENDIF
                schskp = test.LT.cutoff
            ENDIF

#if __TIMING==1
            tim2 = tim2 + omp_get_wtime() - tim1
#endif
            IF(.NOT.schskp) THEN

                    CALL ompmod_shellquart(ii,jj,kk,ll,ghondo_t,ddij_t)

                    thr_nshq = thr_nshq + 1

#if __TIMING==1
                    tim3 = tim3 + omp_get_wtime() - tim1
#endif
!          ---  Fock matrix update ---

!                    CALL dirfck(typscf,ia,da,fa,db,fb,ghondo_t, &
!                                l2a,nint,nflmat)
!                    write(*,'(4i4)') ii,jj,kk,ll
                    CALL ompmod_dirfck_rhf(ii,jj,kk,ll,ghondo_t, &
                                           hfscal,cscalt,cutoff, &
                                           nint,ia,daomp, &
                                           Fi(:,:,ithread), Fj(:,:,ithread), faomp)

#if __TIMING==1
                    tim4 = tim4 + omp_get_wtime() - tim1
#endif
            ELSE
                nschwz = nschwz + 1
            ENDIF
        ENDDO klc
!!$omp enddo
!$omp enddo nowait

!   always need to add J subblock contribution
    CALL ompmod_tools_col_reduce(kmax(jj)-kmin(jj)+1, &
                                 kloc(jj),ithread,nthreads,Fj,faomp)

        iprev = ii
    ENDDO ijc

!$omp barrier

!   add last I subblock contribution
    IF (iprev<=nshell) THEN
        CALL ompmod_tools_col_reduce(kmax(iprev)-kmin(iprev)+1, &
                                     kloc(iprev),ithread,nthreads,Fi,faomp)
    ENDIF
#if __TIMING==1
    tim0 = omp_get_wtime() - tim0

!   debug timing output
    IF (oflag) THEN
!$omp do ordered
        DO ii = 0, max_threads-1
!$omp ordered
            WRITE(iw,dbgfmt2) ithread,thr_nshq,tim3-tim2,tim4-tim3,tim2,tim0
!$omp end ordered
        ENDDO
!$omp enddo
    ENDIF
#endif


!   deallocate temporary memory for integrals
    deallocate(ddij_t,ghondo_t)

!$omp end parallel

    CALL blkint_shmtrx_rev(dsh)
    CALL blkint_shmtrx_rev(xints)
    CALL blkint_aomtrx_rev(faomp)
    CALL blkint_mtrx_spack(faomp,fa,l1)

!   deallocate temporary Fock matrix subblocks
    deallocate(Fi,Fj)

!   deallocate square Fock and density matrices
    deallocate(faomp, daomp)

    CALL ddi_dlbreset

    CALL blkint_shrevert

  END SUBROUTINE ompmod_twoei_shf_kl_rhf

!---------------------------------------------------------------------
!=====================================================================
!---------------------------------------------------------------------

!*MODULE OMPMOD   *DECK OMPMOD_EXCHNG
!
!>    @brief   Calculates integrals for screening purposes.
!>             Threaded version of EXCHNG routine
!>
!>    @details Calculates exhcange ERIs \f$ (i,j|i,j) \f$.
!>             Largest ERI in shell pair is stored in
!>             array `xints(:)` to be used later for screening.
!>             Based on `EXCHNG` subroutine from `int2a.src`
!>             OpenMP algorithm is turned on by default when `INTOMP`
!>             option of the `$INTGRL` group in the input file
!>             is nonzero.
!>
!>    @author  Vladimir Mironov
!
!     REVISION HISTORY:
!>    @date _May, 2016_ Initial release
!>    @date _Jun, 2017_ Bug fixes
!
!     PARAMETERS:
!
!>    @param[out] xints(:) array to store max exchange
!>                         integrals over shells
!>    @param[in]  nsh2     size of `xints(:)`
!>    @param[in]  maxg     dimension of temporary array for
!>                         Rys integral code
!>    @param[in]  inttyp   selects ERI scheme
  SUBROUTINE ompmod_exchng(xints,nsh2,maxg,inttyp)

    REAL(KIND=fp), INTENT(OUT) :: &
        xints(nsh2)

    INTEGER, INTENT(IN) :: &
        nsh2, maxg, inttyp

    COMMON /intac2/ ei1,ei2,cux
        REAL(KIND=fp) :: ei1,ei2,cux

    COMMON /flips / ib(4,3)
        INTEGER :: ib

    COMMON /fmoinf/ nfg,nlayer,natfmo,nbdfg,naotyp,nbody
        INTEGER :: nfg,nlayer,natfmo,nbdfg,naotyp,nbody

    COMMON /fmorun/ espscf,e0scf(2),emp2s,idafmo,icurfg,jcurfg,kcurfg, &
                    icurlay,icurunt,nat1e,ncursh,ngau,icurpop,ifmostp, &
                    moncor,needr,modrst,norbproj,nunesp,iskipesp,      &
                    iesdppc,idoprop,mp2run,icurit,idmfmo,iddfmo,       &
                    iddcur,nddleft,ivmfmo,nzmtfmo,ifmobas,itmfmo(2)
        REAL(KIND=fp) :: espscf,e0scf,emp2s
        INTEGER :: idafmo,icurfg,jcurfg,kcurfg, &
                   icurlay,icurunt,nat1e,ncursh,ngau,icurpop,ifmostp, &
                   moncor,needr,modrst,norbproj,nunesp,iskipesp,      &
                   iesdppc,idoprop,mp2run,icurit,idmfmo,iddfmo,       &
                   iddcur,nddleft,ivmfmo,nzmtfmo,ifmobas,itmfmo

    COMMON /gout  / gpople(768),norgp
        REAL(KIND=fp) :: gpople
        INTEGER :: norgp

    COMMON /iofile/ ir,iw,ip,is,ipk,idaf,nav,ioda(950)
        INTEGER :: ir,iw,ip,is,ipk,idaf,nav,ioda

    COMMON /elgidx/ lcut
        LOGICAL :: lcut

    COMMON /nshel / ex(mxgtot),cs(mxgtot),cp(mxgtot),cd(mxgtot),       &
                    cf(mxgtot),cg(mxgtot),ch(mxgtot),ci(mxgtot),       &
                    kstart(mxsh),katom(mxsh),ktype(mxsh),kng(mxsh),    &
                    kloc(mxsh),kmin(mxsh),kmax(mxsh),nshell
        REAL(KIND=fp) :: ex,cs,cp,cd,cf,cg,ch,ci
        INTEGER :: kstart,katom,ktype,kng,kloc,kmin,kmax,nshell

    COMMON /output/ nprint,itol,icut,normf,normp,nopk
        INTEGER :: nprint,itol,icut,normf,normp,nopk

    COMMON /par   / me,master,nproc,ibtyp,iptim,goparr,dskwrk,maswrk
        INTEGER :: me,master,nproc,ibtyp,iptim
        LOGICAL :: dskwrk, maswrk, goparr

    COMMON /shlexc/ norgsh(3),norgsp(3),iexch,nangm,ngth(4)
        INTEGER :: norgsh,norgsp,iexch,nangm,ngth

    COMMON /shlg70/ ish,jsh,ksh,lsh,ijklxx(4)
        INTEGER :: ish,jsh,ksh,lsh,ijklxx

    COMMON /shlnos/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,           &
                    mini,minj,mink,minl,maxi,maxj,maxk,maxl,           &
                    nij,ij,kl,ijkl
        INTEGER :: lit,ljt,lkt,llt,loci,locj,lock,locl, &
                   mini,minj,mink,minl,maxi,maxj,maxk,maxl, &
                   nij,ij,kl,ijkl
        REAL(KIND=fp) :: qq4

    COMMON /maxc  / cmax(mxgtot),cmaxa(mxgsh),cmaxb(mxgsh), &
                    cmaxc(mxgsh),cmaxd(mxgsh),ismlp(mxg2),ismlq
        REAL(KIND=fp) :: cmax,cmaxa,cmaxb,cmaxc,cmaxd
        INTEGER :: ismlp,ismlq

    COMMON /nlrcf / lrint
        LOGICAL :: lrint

    COMMON /dftpar/ dfttyp(20),exena,exenb,exenc, &
                    idft34,nauxfun,nauxshl
        REAL(KIND=fp) :: dfttyp,exena,exenb,exenc
        INTEGER :: idft34,nauxfun,nauxshl

    COMMON /shlt  / tol,cutoff,icount,out
        REAL(KIND=fp) :: tol, cutoff
        INTEGER :: icount
        LOGICAL :: out

    COMMON /fmcom / xx(1)
        REAL(KIND=fp), TARGET :: xx

!$omp threadprivate(/flips /,/gout  /,/shlnos/,/shlg70/, &
!$omp               /maxc  /,/nlrcf /,/dftpar/)

    REAL(KIND=fp) :: &
        tim, tim0, texch, tolsv, ei1sv, ei2sv, cuxsv, vmax

    REAL(KIND=fp),DIMENSION(:),ALLOCATABLE :: &
        ddij_t, ghondo_t

    LOGICAL :: &
        pople,some,saveint, doesp

    INTEGER :: &
        ithread, lmax, nint, ijij, ii, jj, ish2 

    INTEGER :: LEVELEFP,LVLEFP

    doesp=nfg.NE.0.AND.ncursh.NE.0
    IF (goparr.OR.doesp) xints(1:nsh2)=0

    some = nprint.NE.-5 .AND. maswrk
    saveint=nfg.EQ.0.OR.ncursh.EQ.0.OR.ifmostp.EQ.6

    IF (some) THEN
        tim = zero
        CALL tsecnd(tim)
        tim0 = tim
    ENDIF

    CALL baschk(lmax)
    nangm = angm(lmax)
    ngth(4) = 1
    ngth(3) = ngth(4) * nangm
    ngth(2) = ngth(3) * nangm
    ngth(1) = ngth(2) * nangm

!   The idea is to do even small integrals, below the usual
!   cutoff threshholds, by resetting tolerances tightly.

    tolsv = tol
    tol = 75.0d+00

    ei1sv = ei1
    ei2sv = ei2
    cuxsv = cux
    ei1 = 1.0d-17
    ei2 = 1.0d-17
    cux = 50.0d+00

    iexch = 1
    norgp = 0
    qq4   = one
    nint  = 0

    norgsh(1:3) = 0
    norgsp(1:3) = 0

    ijij = 0


!$omp parallel                                      &
!$omp   private(ii,jj,ijij,ish2,vmax,               &
!$omp     pople,ithread,ghondo_t,ddij_t)            &
!$omp   shared(doesp,ncursh,ktype)                  &
!$omp   reduction(+:nint)                           &
!$omp   copyin(cmax,qq4,norgp,lrint,/dftpar/)


    ithread = omp_get_thread_num()

    allocate(ddij_t(mxang2*mxg2),ghondo_t(maxg))

ic: DO ii = 1,nshell
        ish = ii
!
!   ----- go parallel! -----
!
        IF (goparr) THEN
            IF (mod(ish,nproc).NE.me) CYCLE ic
        ENDIF

        ish2 = ish*(ish-1)/2
!$omp do schedule(dynamic,1)
jc:     DO jj = 1,ii
            jsh = jj
            ijij = ish2 + jj

!           skip unneeded off-diagonal blocks for FMO ESP screening
            IF (doesp.AND.(ish.GT.ncursh).AND.(jsh.LE.ncursh)) CYCLE jc

!           use Pople code for any pure SP integral blocks,
!           use HONDO Rys polynomial code for other blocks


            pople = inttyp.LT.2.AND.max(ktype(ish),ktype(jsh)).LE.2

            IF (pople) THEN
                ksh=ish
                lsh=jsh
                gpople(1:256)=0.0_fp
                CALL genr70(1,.FALSE.)
            ELSE
                ghondo_t(:) = 0.0
                CALL ompmod_rysint(ish,jsh,ish,jsh,ktype,ghondo_t,ddij_t)
            ENDIF

! -----     pick out largest exchange integral for this block -----

            IF (pople) THEN
                vmax = maxval(abs(gpople(1:256)))
                nint = nint+count(abs(gpople(1:256)).gt.zero)
            ELSE
                vmax = maxval(abs(ghondo_t))
                nint = nint+count(abs(ghondo_t).gt.zero)
            ENDIF

            xints(ijij) = sqrt(vmax)
        ENDDO jc
!$omp enddo nowait
    ENDDO ic

!   deallocate temporary memory for integrals
    deallocate(ddij_t,ghondo_t)

!$omp end parallel


! ----- sum up partial contributions if parallel -----

    IF (goparr) THEN
        CALL ddi_sync(1052)
        CALL ddi_gsumf(1050,xints,nsh2)
        CALL ddi_gsumi(1051,nint ,1)
    ENDIF
!
    IF (out) THEN
        WRITE(iw,*) 'MAX EXCHANGE INTEGRAL IN SHELL'
        CALL prtri(xints,nshell)
    ENDIF
!
    LVLEFP=LEVELEFP()
    IF (some) THEN
        CALL tsecnd(tim)
        texch = tim-tim0
        IF(LVLEFP.NE.2) WRITE(iw,dbgfmt_exch) nint,texch
    ENDIF
!
    tol = tolsv
    ei1 = ei1sv
    ei2 = ei2sv
    cux = cuxsv

!   During FMO ESP runs, exchange integrals have different size so
!   one cannot write them to the same record. the only exception is
!   the separated dimer energies where there is just one set of 2e
!   integrals.
!   Elongation method also must decide on this.

    IF (saveint.AND.(.NOT.lcut)) &
        CALL dawrit(idaf,ioda,xints,nsh2,54,0)

  END SUBROUTINE ompmod_exchng

!---------------------------------------------------------------------
!=====================================================================
!---------------------------------------------------------------------

!*MODULE OMPMOD   *DECK OMPMOD_RYSINT
!>    @brief   Rys quadrature ERI calculation
!>
!>    @details Calculates two-electron integrals over
!>             shell quartet \f$ (ii,jj|kk,ll) \f$ using Rys quadrature
!>             code of GAMESS (`SHELLS, IJPRIM, GENRAL` and `S0000`
!>             routines). Based on `SHELLQUART` routine from `int2a.src`
!>             Supports S, P, D, F, G, H and I shells.
!>
!>    @author  Vladimir Mironov
!
!     REVISION HISTORY:
!>    @date _May, 2016_ Initial release
!>    @date _Jun, 2017_ Bug fixes
!
!     PARAMETERS:
!
!>    @param[in]     ii        first shell index
!>    @param[in]     jj        second shell index
!>    @param[in]     kk        third shell index
!>    @param[in]     ll        fourth shell index
!>    @param[in]     ktype(:)  array of shell types according
!>                                 to their angular momentum
!>    @param[out]    ghondo(:) array of integrals over shell quartet
!>    @param[in]     ddij(:)   temporary array for ingeral normalization
  SUBROUTINE ompmod_rysint(ii,jj,kk,ll,ktype,ghondo,ddij)

    INTEGER, INTENT(IN) :: &
      ii,jj,kk,ll,ktype(mxsh)

    REAL(KIND=fp),INTENT(INOUT),ALLOCATABLE :: &
      ghondo(:)

    REAL(KIND=fp),INTENT(INOUT) :: &
      ddij(:)

    INTEGER :: &
      lqsum

!   fill indices for Rys quadrature
    CALL shells(1,ii,jj,kk,ll,.TRUE.)
    CALL shells(2,ii,jj,kk,ll,.TRUE.)
    CALL ijprim(ddij)

!    ghondo(:) = 0.0_fp
    call zqout(ghondo)

    lqsum = ktype(ii) + ktype(jj) +     &
            ktype(kk) + ktype(ll) - 4

!   calculate integrals
    IF (lqsum.EQ.0) THEN
      CALL s0000(ghondo,ddij)
    ELSE
      CALL genral(ghondo,ddij)
    ENDIF

  END SUBROUTINE ompmod_rysint

  SUBROUTINE rysint_gpu(ii,jj,kk,ll,ktype,ghondo,ddij,ngth,c,&
                        kstart,katom,kng,kloc,kmin,kmax,&
                        ex,cs,cp,cd,cf,cg,ch,ci,&
                        NGTI,NGTJ,NGTK,NGTL,&
                        NIJ,IJ,KL,IJKL,&
                        qq4,&
                        AX,AY,AZ,BX,BY,BZ,RAB,CX,CY,CZ,DX,DY,DZ,RCD1,&
                        dfttyp,&
                        ia,da,fa,nint,&
                        l1,l2a,nflmat,&
                        RYSINDEX,MM_rys,ICP_rys,ISA_rys,RYSTYP,&
                        NSA_rys,NCP_rys)

    USE omp_lib
    USE prec, ONLY: fp
    USE mx_limits, ONLY: &
    mxsh,mxg2,mxatm,mxgtot,mxgsh

    !INTEGER, INTENT(IN) :: &
    !l1, l2a,nflmat, ia(l1)

    !for dirfck
    REAL(KIND=fp),INTENT(OUT) :: fa(l2a*nflmat) !FOCK MATRIX
    REAL(KIND=fp),INTENT(INOUT) :: da(l2a) !DENSITY MATRIX
    INTEGER,INTENT(IN) :: ia(l1) !TRIANGULAR MATRIX 
    INTEGER :: l2a,nflmat,l1,nint 

    REAL(KIND=fp),DIMENSION(20) :: DFTTYP

    ! INTEGER, INTENT(IN) :: &
    !   ii,jj,kk,ll,ktype(mxsh)

    INTEGER, INTENT(IN) :: ktype(mxsh)
    INTEGER :: ii,jj,kk,ll

    REAL(KIND=fp),INTENT(INOUT),ALLOCATABLE :: &
      ghondo(:)

    REAL(KIND=fp),INTENT(INOUT) :: &
      ddij(:)

    REAL(KIND=fp),DIMENSION(784) :: DKL,DIJ
    REAL(KIND=fp),DIMENSION(784) :: &
     IJGT,IJX,IJY,IJZ,IK,KLGT,KLX,KLY,KLZ

    REAL(KIND=fp),DIMENSION(MXG2) :: A,R,X1,Y1,Z1
    REAL(KIND=fp),DIMENSION(784) :: IJD

      !for shlexc
      INTEGER,DIMENSION(4) :: NGTH
      !for infoa
      REAL(KIND=fp) :: C(3,mxatm)
      !for roots
      REAL(KIND=fp),DIMENSION(13) :: U,W
      REAL(KIND=fp) :: XX
      integer :: NROOTS
      !for eriout
      INTEGER :: NGTI,NGTJ,NGTK,NGTL
      !shlnos
      INTEGER ::mini,maxi,minj,maxj,mink,maxk,minl,maxl,&
                LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,&
                NIJ,IJ,KL,IJKL
      REAL(KIND=fp) :: qq4  
      !nshell
      INTEGER,DIMENSION(mxsh) :: kstart,katom,kng,kloc,kmin,kmax,nshell
      REAL(KIND=fp),DIMENSION(mxgtot) :: ex,cs,cp,cd,cf,cg,ch,ci

      !shlinf
      double precision :: ga(mxgsh),csa1(mxgsh),cpa1(mxgsh),cda(mxgsh)
      double precision :: cfa(mxgsh),cga(mxgsh),cha(mxgsh),cia(mxgsh)
      double precision :: gb(mxgsh),csb1(mxgsh),cpb1(mxgsh),cdb(mxgsh)
      double precision :: cfb(mxgsh),cgb(mxgsh),chb(mxgsh),cib(mxgsh)
      double precision :: gc(mxgsh),csc1(mxgsh),cpc1(mxgsh),cdc(mxgsh)
      double precision :: cfc(mxgsh),cgc(mxgsh),chc(mxgsh),cic(mxgsh)
      double precision :: gd(mxgsh),csd1(mxgsh),cpd1(mxgsh),cdd(mxgsh)
      double precision :: cfd(mxgsh),cgd(mxgsh),chd(mxgsh),cid(mxgsh)

      double precision :: AX,AY,AZ,BX,BY,BZ,RAB,CX,CY,CZ,DX,DY,DZ,RCD1
      integer :: NGA1,NGB1,NGC1,NGD1
      LOGICAL IANDJ,KANDL,SAME
      integer :: IS1,JS1,KS1,LS1
      !IS1=kstart(ii)
            !rys
      INTEGER :: NSA_rys,ISA_rys !sum angular momentum
      INTEGER :: NCP_rys,ICP_rys !product contraction
      INTEGER,ALLOCATABLE :: RYSINDEX(:,:,:,:) !rys shell index
      INTEGER,ALLOCATABLE :: RYSTYP(:,:) !rystype
      INTEGER :: MM_rys
      
      NSA_rys=9
      DO ISA_rys=1,NSA_rys
      DO ICP_rys=1,NCP_rys

      DO MM_rys=1,RYSTYP(ICP_rys,ISA_rys)

      ii=RYSINDEX(1,MM_rys,ICP_rys,ISA_rys)
      jj=RYSINDEX(2,MM_rys,ICP_rys,ISA_rys)
      kk=RYSINDEX(3,MM_rys,ICP_rys,ISA_rys)
      ll=RYSINDEX(4,MM_rys,ICP_rys,ISA_rys)

      IS1= kstart(ii)
      JS1= kstart(jj)
      KS1= kstart(kk)
      LS1= kstart(ll)

      MINI = kmin(ii)
      MAXI = kmax(ii)
      MINJ = kmin(jj)
      MAXJ = kmax(jj)
      MINK = kmin(kk)
      MAXK = kmax(kk)
      MINL = kmin(ll)
      MAXL = kmax(ll)

      lit= ktype(ii)-1
      ljt= ktype(jj)-1
      lkt= ktype(kk)-1
      llt= ktype(ll)-1

      NGA1= kng(ii)
      NGB1= kng(jj)
      NGC1= kng(kk)
      NGD1= kng(ll)

      LOCI = kloc(ii)-MINI
      LOCJ = kloc(jj)-MINJ
      LOCK = kloc(kk)-MINK
      LOCL = kloc(ll)-MINL

      ga=ex(is1)
      csa1=cs(is1)
      cpa1=cp(is1)
      cda=cd(is1)
      cfa=cf(is1)
      cga=cg(is1)
      cha=ch(is1)
      cia=ci(is1)

      gb=ex(js1)
      csb1=cs(js1)
      cpb1=cp(js1)
      cdb=cd(js1)
      cfb=cf(js1)
      cgb=cg(js1)
      chb=ch(js1)
      cib=ci(js1)

      gc=ex(ks1)
      csc1=cs(ks1)
      cpc1=cp(ks1)
      cdc=cd(ks1)
      cfc=cf(ks1)
      cgc=cg(ks1)
      chc=ch(ks1)
      cic=ci(ks1)

      gd=ex(ls1)
      csd1=cs(ls1)
      cpd1=cp(ls1)
      cdd=cd(ls1)
      cfd=cf(ls1)
      cgd=cg(ls1)
      chd=ch(ls1)
      cid=ci(ls1)

      IANDJ = II .EQ. JJ
      KANDL = KK .EQ. LL
      SAME = (II .EQ. KK) .AND. (JJ .EQ. LL)

!   fill indices for Rys quadrature
   CALL shells_gpu(1,ii,jj,kk,ll,.TRUE.,&
IJGT,IJX,IJY,IJZ,IK,KLGT,KLX,KLY,KLZ,&
IANDJ,KANDL,SAME,&
ii,jj,kk,ll,NGTI,NGTJ,NGTK,NGTL,&
LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,&
mini,maxi,minj,maxj,mink,maxk,minl,maxl,&
NIJ,IJ,KL,IJKL,&
NROOTS,NGTH,C,&
kstart,katom,ktype,kng,kloc,kmin,kmax,nshell,&
ex,cs,cp,cd,cf,cg,ch,ci,&
ga,csa1,cpa1,cda,cfa,cga,cha,cia,&
gb,csb1,cpb1,cdb,cfb,cgb,chb,cib,&
gc,csc1,cpc1,cdc,cfc,cgc,chc,cic,&
gd,csd1,cpd1,cdd,cfd,cgd,chd,cid,&
AX,AY,AZ,BX,BY,BZ,RAB,CX,CY,CZ,DX,DY,DZ,RCD1,&
NGA1,NGB1,NGC1,NGD1)

   CALL shells_gpu(2,ii,jj,kk,ll,.TRUE.,&
IJGT,IJX,IJY,IJZ,IK,KLGT,KLX,KLY,KLZ,&
IANDJ,KANDL,SAME,&
ii,jj,kk,ll,NGTI,NGTJ,NGTK,NGTL,&
LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,&
mini,maxi,minj,maxj,mink,maxk,minl,maxl,&
NIJ,IJ,KL,IJKL,&
NROOTS,NGTH,C,&
kstart,katom,ktype,kng,kloc,kmin,kmax,nshell,&
ex,cs,cp,cd,cf,cg,ch,ci,&
ga,csa1,cpa1,cda,cfa,cga,cha,cia,&
gb,csb1,cpb1,cdb,cfb,cgb,chb,cib,&
gc,csc1,cpc1,cdc,cfc,cgc,chc,cic,&
gd,csd1,cpd1,cdd,cfd,cgd,chd,cid,&
AX,AY,AZ,BX,BY,BZ,RAB,CX,CY,CZ,DX,DY,DZ,RCD1,&
NGA1,NGB1,NGC1,NGD1)

   CALL ijprim_gpu(ddij,a,r,x1,y1,z1,ijd,&
IANDJ,KANDL,SAME,&
LIT,mini,maxi,minj,maxj,mink,maxk,minl,maxl,NIJ,&
ga,csa1,cpa1,cda,cfa,cga,cha,cia,&
gb,csb1,cpb1,cdb,cfb,cgb,chb,cib,&
gc,csc1,cpc1,cdc,cfc,cgc,chc,cic,&
gd,csd1,cpd1,cdd,cfd,cgd,chd,cid,&
AX,AY,AZ,BX,BY,BZ,RAB,CX,CY,CZ,DX,DY,DZ,RCD1,&
NGA1,NGB1,NGC1,NGD1)

   call zqout_gpu(ghondo,&
IJGT,IJX,IJY,IJZ,IK,KLGT,KLX,KLY,KLZ,&
IANDJ,KANDL,SAME,&
mini,maxi,minj,maxj,mink,maxk,minl,maxl)

! calculate integrals
     CALL genral_gpu(ghondo,ddij,dkl,dij,&
IJGT,IJX,IJY,IJZ,IK,KLGT,KLX,KLY,KLZ, &
A,R,X1,Y1,Z1,IJD,iandj,kandl,same,&
LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,&
qq4,mini,maxi,minj,maxj,mink,maxk,minl,maxl,&
NIJ,IJ,KL,IJKL,&
NROOTS,&
ga,csa1,cpa1,cda,cfa,cga,cha,cia,&
gb,csb1,cpb1,cdb,cfb,cgb,chb,cib,&
gc,csc1,cpc1,cdc,cfc,cgc,chc,cic,&
gd,csd1,cpd1,cdd,cfd,cgd,chd,cid,&
AX,AY,AZ,BX,BY,BZ,RAB,CX,CY,CZ,DX,DY,DZ,RCD1,&
NGA1,NGB1,NGC1,NGD1)

     CALL dirfck_gpu_rys(ia,da,fa,ghondo,nint,&
kstart,katom,kng,kloc,kmin,kmax,&
ii,jj,kk,ll,&
NGTI,NGTJ,NGTK,NGTL,dfttyp)

  ENDDO
  ENDDO
  ENDDO

  END SUBROUTINE rysint_gpu

SUBROUTINE rysint_gpu_5(ii,jj,kk,ll,ghondo,ddij,ngth,c,&
                        kstart,katom,kloc,kmin,kmax,&
                        ex,cs,cp,cd,&
                        NGTI,NGTJ,NGTK,NGTL,&
                        NIJ,IJ,KL,IJKL,&
                        qq4,&
                        AX,AY,AZ,BX,BY,BZ,RAB,CX,CY,CZ,DX,DY,DZ,RCD1,&
                        dfttyp,&
                        ia,da,fa,nint,&
                        l1,l2a,nflmat,&
                        IDX_RYS_5,N_RYS_5,NCP_rys)

    USE omp_lib
    USE prec, ONLY: fp
    USE mx_limits, ONLY: &
    mxsh,mxatm,mxgtot

    !for dirfck
    REAL(KIND=fp),INTENT(OUT) :: fa(l2a*nflmat) !FOCK MATRIX
    REAL(KIND=fp),INTENT(INOUT) :: da(l2a) !DENSITY MATRIX
    INTEGER,INTENT(IN) :: ia(l1) !TRIANGULAR MATRIX 
    INTEGER :: l2a,nflmat,l1,nint 

    REAL(KIND=fp),DIMENSION(20) :: DFTTYP

    INTEGER :: ii,jj,kk,ll

    REAL(KIND=fp),INTENT(INOUT),ALLOCATABLE :: &
      ghondo(:)

    ! REAL(KIND=fp),INTENT(INOUT) :: &
    !   ddij(:)

    ! REAL(KIND=fp),DIMENSION(784) :: DKL,DIJ
    ! REAL(KIND=fp),DIMENSION(784) :: ddij
    REAL(KIND=fp),DIMENSION(36) :: DKL,DIJ
    REAL(KIND=fp),DIMENSION(36) :: ddij
    ! REAL(KIND=fp),DIMENSION(784) :: &
    !  IJGT,IJX,IJY,IJZ,IK,KLGT,KLX,KLY,KLZ

     ! integer, DIMENSION(784) :: &
     ! IJGT,IJX,IJY,IJZ,IK,KLGT,KLX,KLY,KLZ
     integer, DIMENSION(36) :: IJGT,IJX,IJY,IJZ,IK,KLGT,KLX,KLY,KLZ

    !REAL(KIND=fp),DIMENSION(MXG2) :: A,R,X1,Y1,Z1
    ! REAL(KIND=fp),DIMENSION(784) :: IJD
    REAL(KIND=fp),DIMENSION(36) :: IJD

      !for shlexc
      INTEGER,DIMENSION(4) :: NGTH
      !for infoa
      REAL(KIND=fp) :: C(3,mxatm)
      !for eriout
      INTEGER :: NGTI,NGTJ,NGTK,NGTL
      !shlnos
      INTEGER ::mini,maxi,minj,maxj,mink,maxk,minl,maxl,&
                LOCI,LOCJ,LOCK,LOCL,&
                NIJ,IJ,KL,IJKL
      REAL(KIND=fp) :: qq4  
      !nshell
      INTEGER,DIMENSION(mxsh) :: kstart,katom,kloc,kmin,kmax,nshell
      REAL(KIND=fp),DIMENSION(mxgtot) :: ex,cs,cp,cd

      double precision :: AX,AY,AZ,BX,BY,BZ,RAB,CX,CY,CZ,DX,DY,DZ,RCD1
      integer :: NGA1,NGB1,NGC1,NGD1
      LOGICAL IANDJ,KANDL,SAME
      integer :: IS1,JS1,KS1,LS1
      !IS1=kstart(ii)
      !rys
      INTEGER,PARAMETER :: NUM_QUART=9169597 !number of quartets
      INTEGER :: IDX_RYS_5(4,NUM_QUART,NCP_rys)
      INTEGER :: N_RYS_5(NCP_rys)

      INTEGER :: NCP_rys,ICP_rys !product contraction

      iNTEGER :: MM_rys
      integer :: NROOTS

      double precision,parameter :: HALF=0.5D00
      double precision,parameter :: FOUR=4.0D00
      !stuff for digestion
      INTEGER :: II2,JJ2,KK2,IL,JK,JL
      DOUBLE PRECISION :: XVAL4, VAL,VAL1,VAL4,XVAL1
      !REAL(KIND=FP) :: CUTOFFAO
      integer :: i_index,ij_index,ijk_index,itmp,nkl,ijkl_index
      integer :: I1_index,I2_index,J1_index,J2_index,K1_index,K2_index,L1_index,L2_index
      integer :: I_count,J_count,K_count,L_count
      !integer :: IJJ,IKK,ILL,JKK,JLL,KLL
      integer :: IKK
      integer :: i,j,k,l,I1,I2,J1,J2,NX,NY,NZ,K1,K2,L11,L2
      !stuff for zqout
      integer :: ijn,kln,n1,nn
      !stuff for shells
      integer :: inu,jnu,knu,lnu,max
      ! REAL(KIND=fp),DIMENSION(84) :: IX,IY,IZ,JX,JY,JZ,KX,KY,KZ,LX,LY,LZ
      integer,DIMENSION(84) :: IX,IY,IZ,JX,JY,JZ,KX,KY,KZ,LX,LY,LZ
      DATA JX /   0,  49,   0,   0,  98,   0,   0,  49,  49,   0, &
               147,   0,   0,  98,  98,  49,   0,  49,   0,  49, &
               196,   0,   0, 147, 147,  49,   0,  49,   0,  98, &
                98,   0,  98,  49,  49, &
               245,   0,   0, 196, 196,  49,   0,  49,   0, 147, &
               147,  98,   0,  98,   0, 147,  49,  49,  98,  98, &
                49, &
               294,   0,   0, 245, 245,  49,   0,  49,   0, 196, &
               196,  98,   0,  98,   0, 196,  49,  49, 147, 147, &
                 0, 147, 147,  98,  49,  98,  49,  98/
      DATA IX /   1, 344,   1,   1, 687,   1,   1, 344, 344,   1, &
              1030,   1,   1, 687, 687, 344,   1, 344,   1, 344, &
              1373,   1,   1,1030,1030, 344,   1, 344,   1, 687, &
               687,   1, 687, 344, 344, &
              1716,   1,   1,1373,1373, 344,   1, 344,   1,1030, &
              1030, 687,   1, 687,   1,1030, 344, 344, 687, 687, &
               344, &
              2059,   1,   1,1716,1716, 344,   1, 344,   1,1373, &
              1373, 687,   1, 687,   1,1373, 344, 344,1030,1030, &
                 1,1030,1030, 687, 344, 687, 344, 687/
      DATA JY /   0,   0,  49,   0,   0,  98,   0,  49,   0,  49, &
                 0, 147,   0,  49,   0,  98,  98,   0,  49,  49, &
                 0, 196,   0,  49,   0, 147, 147,   0,  49,  98, &
                 0,  98,  49,  98,  49, &
                 0, 245,   0,  49,   0, 196, 196,   0,  49,  98, &
                 0, 147, 147,   0,  98,  49, 147,  49,  98,  49, &
                98, &
                 0, 294,   0,  49,   0, 245, 245,   0,  49,  98, &
                 0, 196, 196,   0,  98,  49, 196,  49, 147,   0, &
               147,  98,  49, 147, 147,  49,  98,  98/
      DATA IY /   1,   1, 344,   1,   1, 687,   1, 344,   1, 344, &
                 1,1030,   1, 344,   1, 687, 687,   1, 344, 344, &
                 1,1373,   1, 344,   1,1030,1030,   1, 344, 687, &
                 1, 687, 344, 687, 344, &
                 1,1716,   1, 344,   1,1373,1373,   1, 344, 687, &
                 1,1030,1030,   1, 687, 344,1030, 344, 687, 344, &
               687, &
                 1,2059,   1, 344,   1,1716,1716,   1, 344, 687, &
                 1,1373,1373,   1, 687, 344,1373, 344,1030,   1, &
              1030, 687, 344,1030,1030, 344, 687, 687/
      DATA JZ /   0,   0,   0,  49,   0,   0,  98,   0,  49,  49, &
                 0,   0, 147,   0,  49,   0,  49,  98,  98,  49, &
                 0,   0, 196,   0,  49,   0,  49, 147, 147,   0, &
                98,  98,  49,  49,  98, &
                 0,   0, 245,   0,  49,   0,  49, 196, 196,   0, &
                98,   0,  98, 147, 147,  49,  49, 147,  49,  98, &
                98, &
                 0,   0, 294,   0,  49,   0,  49, 245, 245,   0, &
                98,   0,  98, 196, 196,  49,  49, 196,   0, 147, &
               147,  49,  98,  49,  98, 147, 147,  98/
      DATA IZ /   1,   1,   1, 344,   1,   1, 687,   1, 344, 344, &
                 1,   1,1030,   1, 344,   1, 344, 687, 687, 344, &
                 1,   1,1373,   1, 344,   1, 344,1030,1030,   1, &
               687, 687, 344, 344, 687, &
                 1,   1,1716,   1, 344,   1, 344,1373,1373,   1, &
               687,   1, 687,1030,1030, 344, 344,1030, 344, 687, &
               687, &
                 1,   1,2059,   1, 344,   1, 344,1716,1716,   1, &
               687,   1, 687,1373,1373, 344, 344,1373,   1,1030, &
              1030, 344, 687, 344, 687,1030,1030, 687/
      DATA LX /   0,   1,   0,   0,   2,   0,   0,   1,   1,   0,&
                 3,   0,   0,   2,   2,   1,   0,   1,   0,   1,&
                 4,   0,   0,   3,   3,   1,   0,   1,   0,   2,&
                 2,   0,   2,   1,   1,&
                 5,   0,   0,   4,   4,   1,   0,   1,   0,   3,&
                 3,   2,   0,   2,   0,   3,   1,   1,   2,   2,&
                 1,&
                 6,   0,   0,   5,   5,   1,   0,   1,   0,   4,&
                 4,   2,   0,   2,   0,   4,   1,   1,   3,   3,&
                 0,   3,   3,   2,   1,   2,   1,   2/
      DATA KX /   0,   7,   0,   0,  14,   0,   0,   7,   7,   0,&
                21,   0,   0,  14,  14,   7,   0,   7,   0,   7,&
                28,   0,   0,  21,  21,   7,   0,   7,   0,  14,&
                14,   0,  14,   7,   7,&
                35,   0,   0,  28,  28,   7,   0,   7,   0,  21,&
                21,  14,   0,  14,   0,  21,   7,   7,  14,  14,&
                 7,&
                42,   0,   0,  35,  35,   7,   0,   7,   0,  28,&
                28,  14,   0,  14,   0,  28,   7,   7,  21,  21,&
                 0,  21,  21,  14,   7,  14,   7,  14/
     DATA LY /   0,   0,   1,   0,   0,   2,   0,   1,   0,   1,&
                 0,   3,   0,   1,   0,   2,   2,   0,   1,   1,&
                 0,   4,   0,   1,   0,   3,   3,   0,   1,   2,&
                 0,   2,   1,   2,   1,&
                 0,   5,   0,   1,   0,   4,   4,   0,   1,   2,&
                 0,   3,   3,   0,   2,   1,   3,   1,   2,   1,&
                 2,&
                 0,   6,   0,   1,   0,   5,   5,   0,   1,   2,&
                 0,   4,   4,   0,   2,   1,   4,   1,   3,   0,&
                 3,   2,   1,   3,   3,   1,   2,   2/
      DATA KY /   0,   0,   7,   0,   0,  14,   0,   7,   0,   7,&
                 0,  21,   0,   7,   0,  14,  14,   0,   7,   7,&
                 0,  28,   0,   7,   0,  21,  21,   0,   7,  14,&
                 0,  14,   7,  14,   7,&
                 0,  35,   0,   7,   0,  28,  28,   0,   7,  14,&
                 0,  21,  21,   0,  14,   7,  21,   7,  14,   7,&
                14,&
                 0,  42,   0,   7,   0,  35,  35,   0,   7,  14,&
                 0,  28,  28,   0,  14,   7,  28,   7,  21,   0,&
                21,  14,   7,  21,  21,   7,  14,  14/
      DATA LZ /   0,   0,   0,   1,   0,   0,   2,   0,   1,   1,&
                 0,   0,   3,   0,   1,   0,   1,   2,   2,   1,&
                 0,   0,   4,   0,   1,   0,   1,   3,   3,   0,&
                 2,   2,   1,   1,   2,&
                 0,   0,   5,   0,   1,   0,   1,   4,   4,   0,&
                 2,   0,   2,   3,   3,   1,   1,   3,   1,   2,&
                 2,&
                 0,   0,   6,   0,   1,   0,   1,   5,   5,   0,&
                 2,   0,   2,   4,   4,   1,   1,   4,   0,   3,&
                 3,   1,   2,   1,   2,   3,   3,   2/
      DATA KZ /   0,   0,   0,   7,   0,   0,  14,   0,   7,   7,&
                 0,   0,  21,   0,   7,   0,   7,  14,  14,   7,&
                 0,   0,  28,   0,   7,   0,   7,  21,  21,   0,&
                14,  14,   7,   7,  14,&
                 0,   0,  35,   0,   7,   0,   7,  28,  28,   0,&
                14,   0,  14,  21,  21,   7,   7,  21,   7,  14,&
                14,&
                 0,   0,  42,   0,   7,   0,   7,  35,  35,   0,&
                14,   0,  14,  28,  28,   7,   7,  28,   0,  21,&
                21,   7,  14,   7,  14,  21,  21,  14/
      !stuff for IJPRIM
      integer :: n,nm
      double precision :: ai,arri,axi,xi,ayi,yi,azi,zi
      double precision :: aj,aainv,axj,xj,ayj,yj,azj,zj
      double precision :: A,R,X1,Y1,Z1
      double precision :: aa,tol
      double precision :: csi,cpi,cdi
      double precision :: csj,cpj,cdj
      double precision :: csk,cpk,cdk
      double precision :: csl,cpl,cdl
      double precision :: ak,brrk,akxk,akyk,akzk,axak,ayak,azak,axai,ayai,azai
      double precision :: al,b,bbrrk,bxbk,bybk,bzbk,bxbi,bybi,bzbi,binv
      double precision :: xa,ya,za,xb,yb,zb
      double precision :: dum,dum1,dum2
      double precision :: dxij,dyij,dzij,dxkl,dykl,dzkl
      double precision :: factor
      !integer :: iaa,jb
      double precision, PARAMETER :: SQRT3=1.73205080756888D+00
      double precision, PARAMETER :: SQRT5=2.23606797749979D+00
      double precision, PARAMETER :: SQRT7=2.64575131106459D+00
      double precision, PARAMETER :: SQRT9=3.0D+00
      double precision, PARAMETER :: SQRT11=3.3166247903553998D+00
      double precision, PARAMETER :: PI252=34.986836655250D+00
      !integer, parameter :: lit=3,ljt=3,lkt=3,llt=3
      !DOUBLE PRECISION :: XIN(31213,1000),YIN(31213,1000),ZIN(31213,1000)
      !DOUBLE PRECISION :: XIN(15000,1000),YIN(15000,1000),ZIN(15000,1000)
      !DOUBLE PRECISION :: XIN(15000),YIN(15000),ZIN(15000)
      !size of the array is d*d*d*d*nroots
      DOUBLE PRECISION :: XIN(5000),YIN(5000),ZIN(5000)
      !DOUBLE PRECISION :: XIN(410),YIN(410),ZIN(410)
      !integer :: ia
      !genral
      double precision :: c1x,c2x,c3x,c4x,c1y,c2y,c3y,c4y,c1z,c2z,c3z,c4z
      double precision :: U(5),W(5),XX
      double precision :: RT1,RT2,RT3,RT4,RT5,WW1,WW2,WW3,WW4,WW5
      integer :: m,mx,my,mz,nnn
      double precision :: d1
      integer :: IN1(5) 
      DATA IN1 /   1,  344,   687,   736,  785 /
      integer :: KN(5) 
      DATA KN /   0,  7,   14,   15,  16 /
      integer :: mm,in(5)
      double precision :: expe,u2,f00
      double precision :: duminv,DM2INV,bp01,b00,b10,XCP00,XC00,ycp00,yc00,zcp00,zc00
      integer :: ID,ITEAM,ITHREAD,NTHREADS
      !REAL(KIND=fp) :: HFSCAL,CSCALT
      !delete this later
      integer :: kg,lg
      logical :: double
      !also delete
      integer :: maxj2,maxl2
      !XYZINT_gpu_2222
      logical :: first1,first2,first3,first4
      integer :: i3,i4,i5,k3,k4,mmm,minn,km,ni,nj,nk,nl,iaa,ibb
      double precision :: cp10,c10,cp01,c01
      !ROOT5TS_GPU
      double precision :: xxx,y,e
      double precision,parameter :: R15=1.17581320211778D-01
      double precision,parameter :: PIE4=7.85398163397448D-01
      double precision,parameter :: R25=1.07456201243690D+00
      double precision,parameter :: W25=2.70967405960535D-01
      double precision,parameter :: R35=3.08593744371754D+00
      double precision,parameter :: W35=3.82231610015404D-02
      double precision,parameter :: R45=6.41472973366203D+00
      double precision,parameter :: W45=1.51614186862443D-03
      double precision,parameter :: R55=1.18071894899717D+01
      double precision,parameter :: W55=8.62130526143657D-06


      NROOTS=5

      DO ICP_rys=1,NCP_rys

!$omp target teams distribute parallel do default(none) &
!$omp map(tofrom:fa) &
!$omp map(to:da,ia) &
!$omp map(to:IDX_RYS_5,N_RYS_5) &
!$omp map(to:ngth,c) &
!$omp map(to:dfttyp) &
!$omp map(to:KMAX,KMIN,KLOC,KSTART,KATOM) &
!$omp map(to:EX,CS,CP,CD) &
!$omp map(to:NROOTS,ICP_rys) &
!$omp map(to:ix,iy,iz,jx,jy,jz,kx,ky,kz,lx,ly,lz) &
!$omp map(to:qq4) &
!$omp private(factor) &
!$omp private(xin,yin,zin)&
!$omp private(ii,jj,kk,ll,is1,js1,ks1,ls1) &
!$omp private(mini,minj,mink,minl,maxi,maxj,maxk,maxl) &
!$omp private(loci,locj,lock,locl) &
!$omp private(iandj,kandl,same,max) &
!$omp private(ngti,ngtj,ngtk,ngtl) &
!$omp private(i,j,k,l) &
!it used to be private
!$omp private(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz) &
!!$omp private(lit,ljt,lkt,llt) &
!$omp private(rab,rcd1)&
!$omp private(nx,ny,nz)&
!$omp private(ij,kl,ijgt,klgt,ijkl)&
!$omp private(ijx,ijy,ijz,klx,kly,klz)&
!$omp private(ik)&
!$omp private(ijn,kln,n1,nn)&
!$omp private(ghondo)&
!this is for ijprim, can i simplify this?
!!$omp private(x1,y1,z1,nij,ijd,a,r,ddij)&
!$omp private(x1,y1,z1,nij,ijd,r)&
!!$omp shared(ddij)&
!$omp private(ddij)&
!$omp private(nm,n)&
!$omp private(csi,cpi,cdi)&
!$omp private(csj,cpj,cdj)&
!$omp private(csk,cpk,cdk)&
!$omp private(csl,cpl,cdl)&
!$omp private(ai,arri,axi,ayi,azi)&
!$omp private(aj,aa,aainv,dum1,dum2)&
!genral
!$omp private(dkl)&
!$omp private(dij)&
!$omp private(dxij,dyij,dzij,dxkl,dykl,dzkl)&
!$omp private(ak,brrk)&
!$omp private(al,binv,bbrrk)&
!next two were shared
!$omp private(bxbk,bybk,bzbk,bxbi,bybi,bzbi)&
!$omp private(akxk,akyk,akzk,axak,ayak,azak,axai,ayai,azai)&
!$omp private(a,b,xa,ya,za,xb,yb,zb)&
!next one was shared
!$omp private(c1x,c2x,c3x,c4x,c1y,c2y,c3y,c4y,c1z,c2z,c3z,c4z)&
!$omp private(RT1,RT2,RT3,RT4,RT5,WW1,WW2,WW3,WW4,WW5)&
!!$omp shared(RT1,RT2,RT3,RT4,RT5,WW1,WW2,WW3,WW4,WW5)&
!$omp private(xx,f00)&
!$omp private(expe)&
!!$omp shared(u,w)&
!$omp private(u,w)&
!$omp private(mx,my,mz,nnn,d1)&
!next two were shared
!$omp private(u2)&
!$omp private(duminv,DM2INV,bp01,b00,b10,XCP00,XC00,ycp00,yc00,zcp00,zc00)&
!$omp private(in,mm)&
!ROOT5TS_GPU
!$omp private(xxx,y,e)&
!!$omp shared(R15,R25,R35,R45,R55,PIE4,W25,W35,W45,W55)&
!these have to be shared
!$omp shared(in1,kn)&
!XYZINT_gpu_2222
!$omp private(cp10,c10,cp01,c01,i3,i4,i5,k3,k4,mmm)&
!$omp private(first1,first2,first3,first4,minn,km)&
!$omp private(iaa,ibb,ni,nj,nk,nl)&
!dirfck
!$omp private(nint,val,val1,val4,xval1,xval4)&
!$omp private(ii2,jj2,kk2)&
!$omp private(i1,i2,j1,j2,k1,k2,l1,l2)&
!figure out itmp
!$omp private(nkl,itmp,ikk,jl,jk,il)&
!$omp private(i_index,ij_index,ijk_index,ijkl_index)&
!!$omp private(ID,ITEAM,NTHREADS,ITHREAD)&
!!$omp private(hfscal,cscalt)&
!delete later
!$omp private(kg,lg,double,maxj2,maxl2)

       DO MM_rys=1,N_RYS_5(ICP_rys)
      !write(*,*) "N_RYS_5(ICP_rys)",N_RYS_5(ICP_rys)
      ! ITEAM=omp_get_team_num()
      ! ITHREAD=omp_get_thread_num()
      ! NTHREADS=omp_get_num_threads()

      ! ID=ITEAM*NTHREADS+ITHREAD+1
      !write(*,*) "ID",ID

      ii=IDX_RYS_5(1,MM_rys,ICP_rys)
      jj=IDX_RYS_5(2,MM_rys,ICP_rys)
      kk=IDX_RYS_5(3,MM_rys,ICP_rys)
      ll=IDX_RYS_5(4,MM_rys,ICP_rys)
      !write(*,*) "in loop"
      MINI = kmin(ii)
      MAXI = kmax(ii)
      MINJ = kmin(jj)
      MAXJ = kmax(jj)
      MINK = kmin(kk)
      MAXK = kmax(kk)
      MINL = kmin(ll)
      MAXL = kmax(ll)

      IS1= kstart(ii)
      JS1= kstart(jj)
      KS1= kstart(kk)
      LS1= kstart(ll)

      LOCI = kloc(ii)-MINI
      LOCJ = kloc(jj)-MINJ
      LOCK = kloc(kk)-MINK
      LOCL = kloc(ll)-MINL

      IANDJ = II .EQ. JJ
      KANDL = KK .EQ. LL
      SAME = (II .EQ. KK) .AND. (JJ .EQ. LL)

!   fill indices for Rys quadrature
      NGTI = NGTH(1)
      NGTJ = NGTH(2)
      NGTK = NGTH(3)
      NGTL = NGTH(4)

      I = KATOM(ii)
      J = KATOM(jj)
      K = KATOM(kk)
      L = KATOM(ll)

      AX = C(1,I)
      AY = C(2,I)
      AZ = C(3,I)
      BX = C(1,J)
      BY = C(2,J)
      BZ = C(3,J)
      CX = C(1,K)
      CY = C(2,K)
      CZ = C(3,K)
      DX = C(1,L)
      DY = C(2,L)
      DZ = C(3,L)

      RAB = ((AX-BX)*(AX-BX) + (AY-BY)*(AY-BY) + (AZ-BZ)*(AZ-BZ))
      !write(*,*) "RAB", RAB
      RCD1 = ((CX-DX)*(CX-DX) + (CY-DY)*(CY-DY) + (CZ-DZ)*(CZ-DZ))
      !write(*,*) "RCD", RCD1

      IJ = 0
      DO I = MINI,MAXI
         NX = IX(I)
         NY = IY(I)
         NZ = IZ(I)
         IF (IANDJ) MAXJ = I
         DO J = MINJ,MAXJ
            IJ = IJ+1
            !write(*,*) "IJ", IJ
            IJX(IJ) = NX+JX(J)
            IJY(IJ) = NY+JY(J)
            IJZ(IJ) = NZ+JZ(J)
            IJGT(IJ) = NGTI*(I-MINI)+NGTJ*(J-MINJ)+1
         ENDDO
      ENDDO
   
      KL = 0
      DO K = MINK,MAXK
         NX = KX(K)
         NY = KY(K)
         NZ = KZ(K)
         IF (KANDL) MAXL = K
         DO L = MINL,MAXL
            KL = KL+1
            KLX(KL) = NX+LX(L)
            KLY(KL) = NY+LY(L)
            KLZ(KL) = NZ+LZ(L)
            KLGT(KL) = NGTK*(K-MINK)+NGTL*(L-MINL)
         ENDDO
      ENDDO
      !time is 0.07 here

      MAX = KL
      DO I = 1,IJ
      IF (SAME) MAX = I
      IK(I) = MAX
      ENDDO
      IJKL = IJ*KL
      IF (SAME) IJKL = IJ*(IJ+1)/2

!====================get density for IJ =================
      MAX = MAXJ
      N = 0
      NN = 0
      NM = -2**20
      DO 180 I = MINI,MAXI
         GO TO (100,100,120,120,100,120,120,100,120,120,&
               100,120,120,100,120,120,120,120,120,100,&
               100,120,120,100,120,120,120,120,120,100,&
               120,120,100,120,120,&
               100,120,120,100,120,120,120,120,120,100,&
               120,120,120,120,120,100,120,120,100,120,&
               120,&
               100,120,120,100,120,120,120,120,120,100,&
               120,120,120,120,120,100,120,120,100,120,&
               120,100,120,120,120,120,120,100),I
  100    NM = NN
  120    NN = NM
         IF (IANDJ) MAX = I
         DO 170 J = MINJ,MAX
            GO TO (140,140,160,160,140,160,160,140,160,160,&
                  140,160,160,140,160,160,160,160,160,140,&
                  140,160,160,140,160,160,160,160,160,140,&
                  160,160,140,160,160,&
                  140,160,160,140,160,160,160,160,160,140,&
                  160,160,160,160,160,140,160,160,140,160,&
                  160,&
                  140,160,160,140,160,160,160,160,160,140,&
                  160,160,160,160,160,140,160,160,140,160,&
                  160,140,160,160,160,160,160,140),J
  140       NN = NN+1
  160       N = N+1
            IJD(N) = NN
  170    CONTINUE
  180 CONTINUE

          AI = ex(is1)
          !write(*,*) "AI"
          ARRI = AI*RAB
          AXI = AI*AX
          AYI = AI*AY
          AZI = AI*AZ
          CSI = cs(is1)
          CPI = cp(is1)
          CDI = cd(is1)
          AJ = ex(js1)
          AA = AI+AJ
          AAINV = ONE/AA
          !DUM = AJ*ARRI*AAINV
          !IF (DUM .GT. TOL) GO TO 520
          CSJ = cs(js1)
          CPJ = cp(js1)
          CDJ = cd(js1)

          NN = 0
          R = AJ*ARRI*AAINV
          A = AA
          X1 = (AXI+AJ*BX)*AAINV
          Y1 = (AYI+AJ*BY)*AAINV
          Z1 = (AZI+AJ*BZ)*AAINV


            DUM1 = ZERO
            DUM2 = ZERO
            DO 420 I = MINI,MAXI
               GO TO (200,220,420,420,240,420,420,260,420,420),I
  200          DUM1 = CSI*AAINV
               GO TO 280
  220          DUM1 = CPI*AAINV
               GO TO 280
  240          DUM1 = CDI*AAINV
               GO TO 280
  260          DUM1 = DUM1*SQRT3
               GO TO 280

  280          IF (IANDJ) MAX = I
                 DO 400 J = MINJ,MAX
                  GO TO (300,320,400,400,340,400,400,360,400,400),J
  300             DUM2 = DUM1*CSJ
                  GO TO 380
  320             DUM2 = DUM1*CPJ
                  GO TO 380
  340             DUM2 = DUM1*CDJ
                  GO TO 380
  360             DUM2 = DUM2*SQRT3
                  GO TO 380

  380             NN = NN+1
                  DDIJ(NN) = DUM2
  400          CONTINUE
420          CONTINUE

!====================zero out ghondo===========
      IJN = 0
      DO I = MINI,MAXI
         IF (IANDJ) MAXJ = I
         DO J = MINJ,MAXJ
            IJN = IJN+1
            N1 = IJGT(IJN)
            KLN = 0
            DO K =  MINK,MAXK
               IF (KANDL) MAXL = K
               DO L = MINL,MAXL
                  KLN = KLN+1
                  NN = N1+KLGT(KLN)
                  GHONDO(NN) = 0
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      !0.27 here

!====================calculate integrals===========
      FACTOR = PI252*QQ4
      DXIJ = AX-BX
      DYIJ = AY-BY
      DZIJ = AZ-BZ
      DXKL = CX-DX
      DYKL = CY-DY
      DZKL = CZ-DZ

      AK = ex(ks1)
      BRRK = AK*RCD1
      AKXK = AK*CX
      AKYK = AK*CY
      AKZK = AK*CZ
      CSK = cs(ks1)*FACTOR
      CPK = cp(ks1)*FACTOR
      CDK = cd(ks1)*FACTOR

      AL = ex(ls1)
      B = AK+AL
      BINV = ONE/B
      BBRRK = AL*BRRK*BINV
!             !IF (BBRRK .GT. TOL) GO TO 460
      CSL = cs(ls1)
      CPL = cp(ls1)
      CDL = cd(ls1)
      XB = (AKXK+AL*DX)*BINV
      YB = (AKYK+AL*DY)*BINV
      ZB = (AKZK+AL*DZ)*BINV
      BXBK = B*(XB-CX)
      BYBK = B*(YB-CY)
      BZBK = B*(ZB-CZ)
      BXBI = B*(XB-AX)
      BYBI = B*(YB-AY)
      BZBI = B*(ZB-AZ)

            KG=1
            LG=1

            DOUBLE=KANDL.AND.KG.NE.LG
            N = 0
            MAX = MAXL
            DUM1 = ZERO
            DUM2 = ZERO
            DO K = MINK,MAXK
               GO TO (540,560,620,620,580,620,620,600,620,620,&
                     601),K
  540          DUM1 = CSK*BINV
               GO TO 620
  560          DUM1 = CPK*BINV
               GO TO 620
  580          DUM1 = CDK*BINV
               GO TO 620
  600          DUM1 = DUM1*SQRT3
               GO TO 620
  601          GO TO 620

  620          IF (KANDL) MAX = K
               DO L = MINL,MAX
                  GO TO (740,780,840,840,800,840,840,820,840,840,&
                        821),L
  740             DUM2 = DUM1*CSL
                  GO TO 840
                  IF ( .NOT. DOUBLE) GO TO 840
                  IF (K .GT. 1) GO TO 760
                  DUM2 = DUM2+DUM2
                  GO TO 840
  760             DUM2 = DUM2+CSK*CPL*BINV
                  GO TO 840
  780             DUM2 = DUM1*CPL
                  IF (DOUBLE) DUM2 = DUM2+DUM2
                  GO TO 840
  800             DUM2 = DUM1*CDL
                  IF (DOUBLE) DUM2 = DUM2+DUM2
                  GO TO 840
  820             DUM2 = DUM2*SQRT3
                  GO TO 840
  821             GO TO 840

  840             N = N+1
                  DKL(N) = DUM2
              ENDDO
          ENDDO
!write(*,*) "here"
            NN = 0
               DO I = 1,IJ
                  DIJ(I) = DDIJ(IJD(I)+NN)
                  !write(*,*) "DIJ", DIJ(I)
               ENDDO
               A = AA
               !AB = A*B
               !AANDB = A+B
               ! EXPE = EXP(-DUM)/SQRT(AANDB)
               !write(*,*) "R is", R
               EXPE = EXP(-(BBRRK+R))/SQRT((A+B))
               XA = X1
               YA = Y1
               ZA = Z1
               !XX = RHO*((XA-XB)*(XA-XB) + (YA-YB)*(YA-YB)+ (ZA-ZB)*(ZA-ZB))
               XX = ((A*B)/(A+B))*((XA-XB)*(XA-XB) + (YA-YB)*(YA-YB)+ (ZA-ZB)*(ZA-ZB))
               !write(*,*) "XX in 2222 is", XX
               !make these shared?
               AXAK = A*(XA-CX)
               AYAK = A*(YA-CY)
               AZAK = A*(ZA-CZ)
               AXAI = A*(XA-AX)
               AYAI = A*(YA-AY)
               AZAI = A*(ZA-AZ)
               C1X = BXBK+AXAK
               C2X = A*BXBK
               C3X = BXBI+AXAI
               C4X = B*AXAI
               C1Y = BYBK+AYAK
               C2Y = A*BYBK
               C3Y = BYBI+AYAI
               C4Y = B*AYAI
               C1Z = BZBK+AZAK
               C2Z = A*BZBK
               C3Z = BZBI+AZAI
               C4Z = B*AZAI
               !write(*,*) "C4Z", C4Z
!===================calculate roots===============
      IF (XX .GT. 15.0D+00) GO TO 1800
      IF (XX .GT. 5.0D+00) GO TO 1400
      IF (XX .GT. 1.0D+00) GO TO 1200
      IF (XX .GT. 3.0D-07) GO TO 1000
!C     XX IS APPROXXIMATELY ZERO.                   NROOTS = 5
      RT1 = 2.26659266316985D-02 -2.15865967920897D-03 *XX
      RT2 = 2.31271692140903D-01 -2.20258754389745D-02 *XX
      RT3 = 8.57346024118836D-01 -8.16520023025515D-02 *XX
      RT4 = 2.97353038120346D+00 -2.83193369647137D-01 *XX
      RT5 = 1.84151859759051D+01 -1.75382723579439D+00 *XX
      WW1 = 2.95524224714752D-01 -1.96867576909777D-02 *XX
      WW2 = 2.69266719309995D-01 -5.61737590184721D-02 *XX
      WW3 = 2.19086362515981D-01 -9.71152726793658D-02 *XX
      WW4 = 1.49451349150580D-01 -1.02979262193565D-01 *XX
      WW5 = 6.66713443086877D-02 -5.73782817488315D-02 *XX
      !RETURN
      GO TO 7777


!C     XX=0.0 TO 1.0                               NROOTS = 5
  1000 RT1 = ((((((-4.46679165328413D-11*XX+1.21879111988031D-09)*XX- &
          2.62975022612104D-08 )*XX+5.15106194905897D-07 )*XX- &
          9.27933625824749D-06 )*XX+1.51794097682482D-04 )*XX- &
          2.15865967920301D-03 )*XX+2.26659266316985D-02
      RT2 = (((((( 1.93117331714174D-10*XX-4.57267589660699D-09)*XX+ &
          2.48339908218932D-08 )*XX+1.50716729438474D-06 )*XX- &
          6.07268757707381D-05 )*XX+1.37506939145643D-03 )*XX- &
          2.20258754419939D-02 )*XX+2.31271692140905D-01
      RT3 = ((((( 4.84989776180094D-09*XX+1.31538893944284D-07)*XX- &
          2.766753852879D-06)*XX-7.651163510626D-05)*XX+ &
          4.033058545972D-03)*XX-8.16520022916145D-02 )*XX+ &
          8.57346024118779D-01
      RT4 = ((((-2.48581772214623D-07*XX-4.34482635782585D-06)*XX- &
          7.46018257987630D-07 )*XX+1.01210776517279D-02 )*XX- &
          2.83193369640005D-01 )*XX+2.97353038120345D+00
      RT5 = (((((-8.92432153868554D-09*XX+1.77288899268988D-08)*XX+ &
          3.040754680666D-06)*XX+1.058229325071D-04)*XX+ &
          4.596379534985D-02)*XX-1.75382723579114D+00 )*XX+ &
          1.84151859759049D+01
      WW1 = ((((((-2.03822632771791D-09*XX+3.89110229133810D-08)*XX- &
          5.84914787904823D-07 )*XX+8.30316168666696D-06 )*XX- &
          1.13218402310546D-04 )*XX+1.49128888586790D-03 )*XX- &
          1.96867576904816D-02 )*XX+2.95524224714749D-01
      WW2 = ((((((( 8.62848118397570D-09*XX-1.38975551148989D-07)*XX+ &
          1.602894068228D-06)*XX-1.646364300836D-05)*XX+ &
          1.538445806778D-04)*XX-1.28848868034502D-03 )*XX+ &
          9.38866933338584D-03 )*XX-5.61737590178812D-02 )*XX+ &
          2.69266719309991D-01
      WW3 = ((((((((-9.41953204205665D-09*XX+1.47452251067755D-07)*XX- &
          1.57456991199322D-06 )*XX+1.45098401798393D-05 )*XX- &
          1.18858834181513D-04 )*XX+8.53697675984210D-04 )*XX- &
          5.22877807397165D-03 )*XX+2.60854524809786D-02 )*XX- &
          9.71152726809059D-02 )*XX+2.19086362515979D-01
      WW4 = ((((((((-3.84961617022042D-08*XX+5.66595396544470D-07)*XX- &
          5.52351805403748D-06 )*XX+4.53160377546073D-05 )*XX- &
          3.22542784865557D-04 )*XX+1.95682017370967D-03 )*XX- &
          9.77232537679229D-03 )*XX+3.79455945268632D-02 )*XX- &
          1.02979262192227D-01 )*XX+1.49451349150573D-01
      WW5 = ((((((((( 4.09594812521430D-09*XX-6.47097874264417D-08)*XX+ &
          6.743541482689D-07)*XX-5.917993920224D-06)*XX+ &
          4.531969237381D-05)*XX-2.99102856679638D-04 )*XX+ &
          1.65695765202643D-03 )*XX-7.40671222520653D-03 )*XX+ &
          2.50889946832192D-02 )*XX-5.73782817487958D-02 )*XX+ &
          6.66713443086877D-02
          GO TO 7777
!C     XXXX=1.0 TO 5.0                               NROOTS = 5
  1200 Y = XX-3.0D+00
      RT1 = ((((((((-2.58163897135138D-14*Y+8.14127461488273D-13)*Y- &
          2.11414838976129D-11 )*Y+5.09822003260014D-10 )*Y- &
          1.16002134438663D-08 )*Y+2.46810694414540D-07 )*Y- &
          4.92556826124502D-06 )*Y+9.02580687971053D-05 )*Y- &
          1.45190025120726D-03 )*Y+1.73416786387475D-02
      RT2 = ((((((((( 1.04525287289788D-14*Y+5.44611782010773D-14)*Y- &
          4.831059411392D-12)*Y+1.136643908832D-10)*Y- &
          1.104373076913D-09)*Y-2.35346740649916D-08 )*Y+ &
          1.43772622028764D-06 )*Y-4.23405023015273D-05 )*Y+ &
          9.12034574793379D-04 )*Y-1.52479441718739D-02 )*Y+ &
          1.76055265928744D-01
      RT3 = (((((((((-6.89693150857911D-14*Y+5.92064260918861D-13)*Y+ &
          1.847170956043D-11)*Y-3.390752744265D-10)*Y- &
          2.995532064116D-09)*Y+1.57456141058535D-07 )*Y- &
          3.95859409711346D-07 )*Y-9.58924580919747D-05 )*Y+ &
          3.23551502557785D-03 )*Y-5.97587007636479D-02 )*Y+ &
          6.46432853383057D-01
      RT4 = ((((((((-3.61293809667763D-12*Y-2.70803518291085D-11)*Y+ &
          8.83758848468769D-10 )*Y+1.59166632851267D-08 )*Y- &
          1.32581997983422D-07 )*Y-7.60223407443995D-06 )*Y- &
          7.41019244900952D-05 )*Y+9.81432631743423D-03 )*Y- &
          2.23055570487771D-01 )*Y+2.21460798080643D+00
      RT5 = ((((((((( 7.12332088345321D-13*Y+3.16578501501894D-12)*Y- &
          8.776668218053D-11)*Y-2.342817613343D-09)*Y- &
          3.496962018025D-08)*Y-3.03172870136802D-07 )*Y+ &
          1.50511293969805D-06 )*Y+1.37704919387696D-04 )*Y+ &
          4.70723869619745D-02 )*Y-1.47486623003693D+00 )*Y+ &
          1.35704792175847D+01
      WW1 = ((((((((( 1.04348658616398D-13*Y-1.94147461891055D-12)*Y+ &
          3.485512360993D-11)*Y-6.277497362235D-10)*Y+ &
          1.100758247388D-08)*Y-1.88329804969573D-07 )*Y+ &
          3.12338120839468D-06 )*Y-5.04404167403568D-05 )*Y+ &
          8.00338056610995D-04 )*Y-1.30892406559521D-02 )*Y+ &
          2.47383140241103D-01
      WW2 = ((((((((((( 3.23496149760478D-14*Y-5.24314473469311D-13)*Y+ &
          7.743219385056D-12)*Y-1.146022750992D-10)*Y+ &
          1.615238462197D-09)*Y-2.15479017572233D-08 )*Y+ &
          2.70933462557631D-07 )*Y-3.18750295288531D-06 )*Y+ &
          3.47425221210099D-05 )*Y-3.45558237388223D-04 )*Y+ &
          3.05779768191621D-03 )*Y-2.29118251223003D-02 )*Y+ &
          1.59834227924213D-01
      WW3 = ((((((((((((-3.42790561802876D-14*Y+5.26475736681542D-13)*Y- &
          7.184330797139D-12)*Y+9.763932908544D-11)*Y- &
          1.244014559219D-09)*Y+1.472744068942D-08)*Y- &
          1.611749975234D-07)*Y+1.616487851917D-06)*Y- &
          1.46852359124154D-05 )*Y+1.18900349101069D-04 )*Y- &
          8.37562373221756D-04 )*Y+4.93752683045845D-03 )*Y- &
          2.25514728915673D-02 )*Y+6.95211812453929D-02
      WW4 = ((((((((((((( 1.04072340345039D-14*Y-1.60808044529211D-13)* &
          Y+2.183534866798D-12)*Y-2.939403008391D-11)*Y+ &
          3.679254029085D-10)*Y-4.23775673047899D-09 )*Y+ &
          4.46559231067006D-08 )*Y-4.26488836563267D-07 )*Y+ &
          3.64721335274973D-06 )*Y-2.74868382777722D-05 )*Y+ &
          1.78586118867488D-04 )*Y-9.68428981886534D-04 )*Y+ &
          4.16002324339929D-03 )*Y-1.28290192663141D-02 )*Y+ &
          2.22353727685016D-02
      WW5 = ((((((((((((((-8.16770412525963D-16*Y+1.31376515047977D-14)* &
          Y-1.856950818865D-13)*Y+2.596836515749D-12)*Y- &
          3.372639523006D-11)*Y+4.025371849467D-10)*Y- &
          4.389453269417D-09)*Y+4.332753856271D-08)*Y- &
          3.82673275931962D-07 )*Y+2.98006900751543D-06 )*Y- &
          2.00718990300052D-05 )*Y+1.13876001386361D-04 )*Y- &
          5.23627942443563D-04 )*Y+1.83524565118203D-03 )*Y- &
          4.37785737450783D-03 )*Y+5.36963805223095D-03
! !       RETURN
          GO TO 7777

  1400 IF (XX .GT. 10.0D+00) GO TO 1600
!C     XXXX=5.0 TO 10.0                              NROOTS = 5
      Y = XX-7.5D+00
      RT1 = ((((((((-1.13825201010775D-14*Y+1.89737681670375D-13)*Y- &
          4.81561201185876D-12 )*Y+1.56666512163407D-10 )*Y- &
          3.73782213255083D-09 )*Y+9.15858355075147D-08 )*Y- &
          2.13775073585629D-06 )*Y+4.56547356365536D-05 )*Y- &
          8.68003909323740D-04 )*Y+1.22703754069176D-02
      RT2 = (((((((((-3.67160504428358D-15*Y+1.27876280158297D-14)*Y- &
          1.296476623788D-12)*Y+1.477175434354D-11)*Y+ &
          5.464102147892D-10)*Y-2.42538340602723D-08 )*Y+ &
          8.20460740637617D-07 )*Y-2.20379304598661D-05 )*Y+ &
          4.90295372978785D-04 )*Y-9.14294111576119D-03 )*Y+ &
          1.22590403403690D-01
      RT3 = ((((((((( 1.39017367502123D-14*Y-6.96391385426890D-13)*Y+ &
          1.176946020731D-12)*Y+1.725627235645D-10)*Y- &
          3.686383856300D-09)*Y+2.87495324207095D-08 )*Y+ &
          1.71307311000282D-06 )*Y-7.94273603184629D-05 )*Y+ &
          2.00938064965897D-03 )*Y-3.63329491677178D-02 )*Y+ &
          4.34393683888443D-01
      RT4 = ((((((((((-1.27815158195209D-14*Y+1.99910415869821D-14)*Y+ &
          3.753542914426D-12)*Y-2.708018219579D-11)*Y- &
          1.190574776587D-09)*Y+1.106696436509D-08)*Y+ &
          3.954955671326D-07)*Y-4.398596059588D-06)*Y- &
          2.01087998907735D-04 )*Y+7.89092425542937D-03 )*Y- &
          1.42056749162695D-01 )*Y+1.39964149420683D+00
      RT5 = ((((((((((-1.19442341030461D-13*Y-2.34074833275956D-12)*Y+ &
          6.861649627426D-12)*Y+6.082671496226D-10)*Y+ &
          5.381160105420D-09)*Y-6.253297138700D-08)*Y- &
          2.135966835050D-06)*Y-2.373394341886D-05)*Y+ &
          2.88711171412814D-06 )*Y+4.85221195290753D-02 )*Y- &
          1.04346091985269D+00 )*Y+7.89901551676692D+00
      WW1 = ((((((((( 7.95526040108997D-15*Y-2.48593096128045D-13)*Y+ &
          4.761246208720D-12)*Y-9.535763686605D-11)*Y+ &
          2.225273630974D-09)*Y-4.49796778054865D-08 )*Y+ &
          9.17812870287386D-07 )*Y-1.86764236490502D-05 )*Y+ &
          3.76807779068053D-04 )*Y-8.10456360143408D-03 )*Y+ &
          2.01097936411496D-01
      WW2 = ((((((((((( 1.25678686624734D-15*Y-2.34266248891173D-14)*Y+ & 
          3.973252415832D-13)*Y-6.830539401049D-12)*Y+ &
          1.140771033372D-10)*Y-1.82546185762009D-09 )*Y+ &
          2.77209637550134D-08 )*Y-4.01726946190383D-07 )*Y+ &
          5.48227244014763D-06 )*Y-6.95676245982121D-05 )*Y+ &
          8.05193921815776D-04 )*Y-8.15528438784469D-03 )*Y+ &
          9.71769901268114D-02
      WW3 = ((((((((((((-8.20929494859896D-16*Y+1.37356038393016D-14)*Y- &
          2.022863065220D-13)*Y+3.058055403795D-12)*Y- &
          4.387890955243D-11)*Y+5.923946274445D-10)*Y- &
          7.503659964159D-09)*Y+8.851599803902D-08)*Y- &
          9.65561998415038D-07 )*Y+9.60884622778092D-06 )*Y- &
          8.56551787594404D-05 )*Y+6.66057194311179D-04 )*Y- &
          4.17753183902198D-03 )*Y+2.25443826852447D-02
      WW4 = ((((((((((((((-1.08764612488790D-17*Y+1.85299909689937D-16)* &
          Y-2.730195628655D-15)*Y+4.127368817265D-14)*Y- &
          5.881379088074D-13)*Y+7.805245193391D-12)*Y- &
          9.632707991704D-11)*Y+1.099047050624D-09)*Y- &
          1.15042731790748D-08 )*Y+1.09415155268932D-07 )*Y- &
          9.33687124875935D-07 )*Y+7.02338477986218D-06 )*Y- &
          4.53759748787756D-05 )*Y+2.41722511389146D-04 )*Y- &
          9.75935943447037D-04 )*Y+2.57520532789644D-03
      WW5 = ((((((((((((((( 7.28996979748849D-19*Y-1.26518146195173D-17) &
          *Y+1.886145834486D-16)*Y-2.876728287383D-15)*Y+ &
          4.114588668138D-14)*Y-5.44436631413933D-13 )*Y+ &
          6.64976446790959D-12 )*Y-7.44560069974940D-11 )*Y+ &
          7.57553198166848D-10 )*Y-6.92956101109829D-09 )*Y+ &
          5.62222859033624D-08 )*Y-3.97500114084351D-07 )*Y+ &
          2.39039126138140D-06 )*Y-1.18023950002105D-05 )*Y+ &
          4.52254031046244D-05 )*Y-1.21113782150370D-04 )*Y+ &
          1.75013126731224D-04
!       RETURN
          GO TO 7777

! !     XXXX=10.0 TO 15.0                             NROOTS = 5
  1600 Y = XX-12.5D+00
      RT1 = ((((((((((-4.16387977337393D-17*Y+7.20872997373860D-16)*Y+ & 
          1.395993802064D-14)*Y+3.660484641252D-14)*Y- &
          4.154857548139D-12)*Y+2.301379846544D-11)*Y- &
          1.033307012866D-09)*Y+3.997777641049D-08)*Y- &
          9.35118186333939D-07 )*Y+2.38589932752937D-05 )*Y- &
          5.35185183652937D-04 )*Y+8.85218988709735D-03
      RT2 = ((((((((((-4.56279214732217D-16*Y+6.24941647247927D-15)*Y+ &
          1.737896339191D-13)*Y+8.964205979517D-14)*Y- &
          3.538906780633D-11)*Y+9.561341254948D-11)*Y- &
          9.772831891310D-09)*Y+4.240340194620D-07)*Y- &
          1.02384302866534D-05 )*Y+2.57987709704822D-04 )*Y- &
          5.54735977651677D-03 )*Y+8.68245143991948D-02
      RT3 = ((((((((((-2.52879337929239D-15*Y+2.13925810087833D-14)*Y+ &
          7.884307667104D-13)*Y-9.023398159510D-13)*Y- &
          5.814101544957D-11)*Y-1.333480437968D-09)*Y- &
          2.217064940373D-08)*Y+1.643290788086D-06)*Y- &
          4.39602147345028D-05 )*Y+1.08648982748911D-03 )*Y- &
          2.13014521653498D-02 )*Y+2.94150684465425D-01
      RT4 = ((((((((((-6.42391438038888D-15*Y+5.37848223438815D-15)*Y+ &
          8.960828117859D-13)*Y+5.214153461337D-11)*Y- &
          1.106601744067D-10)*Y-2.007890743962D-08)*Y+ &
          1.543764346501D-07)*Y+4.520749076914D-06)*Y- &
          1.88893338587047D-04 )*Y+4.73264487389288D-03 )*Y- &
          7.91197893350253D-02 )*Y+8.60057928514554D-01
      RT5 = (((((((((((-2.24366166957225D-14*Y+4.87224967526081D-14)*Y+ &
          5.587369053655D-12)*Y-3.045253104617D-12)*Y- &
          1.223983883080D-09)*Y-2.05603889396319D-09 )*Y+ &
          2.58604071603561D-07 )*Y+1.34240904266268D-06 )*Y- &
          5.72877569731162D-05 )*Y-9.56275105032191D-04 )*Y+ &
          4.23367010370921D-02 )*Y-5.76800927133412D-01 )*Y+ &
          3.87328263873381D+00
      WW1 = ((((((((( 8.98007931950169D-15*Y+7.25673623859497D-14)*Y+ &
          5.851494250405D-14)*Y-4.234204823846D-11)*Y+ &
          3.911507312679D-10)*Y-9.65094802088511D-09 )*Y+ &
          3.42197444235714D-07 )*Y-7.51821178144509D-06 )*Y+ &
          1.94218051498662D-04 )*Y-5.38533819142287D-03 )*Y+ &
          1.68122596736809D-01
      WW2 = ((((((((((-1.05490525395105D-15*Y+1.96855386549388D-14)*Y- &
          5.500330153548D-13)*Y+1.003849567976D-11)*Y- &
          1.720997242621D-10)*Y+3.533277061402D-09)*Y- &
          6.389171736029D-08)*Y+1.046236652393D-06)*Y- &
          1.73148206795827D-05 )*Y+2.57820531617185D-04 )*Y- &
          3.46188265338350D-03 )*Y+7.03302497508176D-02
      WW3 = ((((((((((( 3.60020423754545D-16*Y-6.24245825017148D-15)*Y+ &
          9.945311467434D-14)*Y-1.749051512721D-12)*Y+ &
          2.768503957853D-11)*Y-4.08688551136506D-10 )*Y+ &
          6.04189063303610D-09 )*Y-8.23540111024147D-08 )*Y+ &
          1.01503783870262D-06 )*Y-1.20490761741576D-05 )*Y+ &
          1.26928442448148D-04 )*Y-1.05539461930597D-03 )*Y+ &
          1.15543698537013D-02
      WW4 = ((((((((((((( 2.51163533058925D-18*Y-4.31723745510697D-17)* &
          Y+6.557620865832D-16)*Y-1.016528519495D-14)*Y+ &
          1.491302084832D-13)*Y-2.06638666222265D-12 )*Y+ &
          2.67958697789258D-11 )*Y-3.23322654638336D-10 )*Y+ &
          3.63722952167779D-09 )*Y-3.75484943783021D-08 )*Y+ &
          3.49164261987184D-07 )*Y-2.92658670674908D-06 )*Y+ &
          2.12937256719543D-05 )*Y-1.19434130620929D-04 )*Y+ &
          6.45524336158384D-04
      WW5 = ((((((((((((((-1.29043630202811D-19*Y+2.16234952241296D-18)* &
          Y-3.107631557965D-17)*Y+4.570804313173D-16)*Y- &
          6.301348858104D-15)*Y+8.031304476153D-14)*Y- &
          9.446196472547D-13)*Y+1.018245804339D-11)*Y- &
          9.96995451348129D-11 )*Y+8.77489010276305D-10 )*Y- &
          6.84655877575364D-09 )*Y+4.64460857084983D-08 )*Y- & 
          2.66924538268397D-07 )*Y+1.24621276265907D-06 )*Y- &
          4.30868944351523D-06 )*Y+9.94307982432868D-06
! !       RETURN
          GO TO 7777

  1800 IF (XX .GT. 25.0D+00) GO TO 2200
      IF (XX .GT. 20.0D+00) GO TO 2000
!     XXXX=15.0 TO 20.0                             NROOTS = 5
      Y = XX-17.5D+00
      RT1 = (((((((((( 1.91875764545740D-16*Y+7.8357401095707D-16)*Y- &
          3.260875931644D-14)*Y-1.186752035569D-13)*Y+ & 
          4.275180095653D-12)*Y+3.357056136731D-11)*Y- &
          1.123776903884D-09)*Y+1.231203269887D-08)*Y- &
          3.99851421361031D-07 )*Y+1.45418822817771D-05 )*Y- &
          3.49912254976317D-04 )*Y+6.67768703938812D-03
      RT2 = (((((((((( 2.02778478673555D-15*Y+1.01640716785099D-14)*Y- &
          3.385363492036D-13)*Y-1.615655871159D-12)*Y+ &
          4.527419140333D-11)*Y+3.853670706486D-10)*Y- &
          1.184607130107D-08)*Y+1.347873288827D-07)*Y- &
          4.47788241748377D-06 )*Y+1.54942754358273D-04 )*Y- &
          3.55524254280266D-03 )*Y+6.44912219301603D-02
      RT3 = (((((((((( 7.79850771456444D-15*Y+6.00464406395001D-14)*Y- &
          1.249779730869D-12)*Y-1.020720636353D-11)*Y+ &
          1.814709816693D-10)*Y+1.766397336977D-09)*Y- &
          4.603559449010D-08)*Y+5.863956443581D-07)*Y- &
          2.03797212506691D-05 )*Y+6.31405161185185D-04 )*Y- &
          1.30102750145071D-02 )*Y+2.10244289044705D-01
      RT4 = (((((((((((-2.92397030777912D-15*Y+1.94152129078465D-14)*Y+ &
          4.859447665850D-13)*Y-3.217227223463D-12)*Y- &
          7.484522135512D-11)*Y+7.19101516047753D-10 )*Y+ &
          6.88409355245582D-09 )*Y-1.44374545515769D-07 )*Y+ &
          2.74941013315834D-06 )*Y-1.02790452049013D-04 )*Y+ &
          2.59924221372643D-03 )*Y-4.35712368303551D-02 )*Y+ &
          5.62170709585029D-01
      RT5 = ((((((((((( 1.17976126840060D-14*Y+1.24156229350669D-13)*Y- &
          3.892741622280D-12)*Y-7.755793199043D-12)*Y+ &
          9.492190032313D-10)*Y-4.98680128123353D-09 )*Y- &
          1.81502268782664D-07 )*Y+2.69463269394888D-06 )*Y+ &
          2.50032154421640D-05 )*Y-1.33684303917681D-03 )*Y+ &
          2.29121951862538D-02 )*Y-2.45653725061323D-01 )*Y+ &
          1.89999883453047D+00
      WW1 = (((((((((( 1.74841995087592D-15*Y-6.95671892641256D-16)*Y- &
          3.000659497257D-13)*Y+2.021279817961D-13)*Y+ &
          3.853596935400D-11)*Y+1.461418533652D-10)*Y- &
          1.014517563435D-08)*Y+1.132736008979D-07)*Y- &
          2.86605475073259D-06 )*Y+1.21958354908768D-04 )*Y- &
          3.86293751153466D-03 )*Y+1.45298342081522D-01
      WW2 = ((((((((((-1.11199320525573D-15*Y+1.85007587796671D-15)*Y+ &
          1.220613939709D-13)*Y+1.275068098526D-12)*Y- &
          5.341838883262D-11)*Y+6.161037256669D-10)*Y- &
          1.009147879750D-08)*Y+2.907862965346D-07)*Y- &
          6.12300038720919D-06 )*Y+1.00104454489518D-04 )*Y- &
          1.80677298502757D-03 )*Y+5.78009914536630D-02
      WW3 = ((((((((((-9.49816486853687D-16*Y+6.67922080354234D-15)*Y+ &
          2.606163540537D-15)*Y+1.983799950150D-12)*Y- &
          5.400548574357D-11)*Y+6.638043374114D-10)*Y- &
          8.799518866802D-09)*Y+1.791418482685D-07)*Y- &
          2.96075397351101D-06 )*Y+3.38028206156144D-05 )*Y- &
          3.58426847857878D-04 )*Y+8.39213709428516D-03
      WW4 = ((((((((((( 1.33829971060180D-17*Y-3.44841877844140D-16)*Y+ &
          4.745009557656D-15)*Y-6.033814209875D-14)*Y+ &
          1.049256040808D-12)*Y-1.70859789556117D-11 )*Y+ &
          2.15219425727959D-10 )*Y-2.52746574206884D-09 )*Y+ &
          3.27761714422960D-08 )*Y-3.90387662925193D-07 )*Y+ &
          3.46340204593870D-06 )*Y-2.43236345136782D-05 )*Y+ &
          3.54846978585226D-04
      WW5 = ((((((((((((( 2.69412277020887D-20*Y-4.24837886165685D-19)* &
          Y+6.030500065438D-18)*Y-9.069722758289D-17)*Y+ &
          1.246599177672D-15)*Y-1.56872999797549D-14 )*Y+ & 
          1.87305099552692D-13 )*Y-2.09498886675861D-12 )*Y+ &
          2.11630022068394D-11 )*Y-1.92566242323525D-10 )*Y+ &
          1.62012436344069D-09 )*Y-1.23621614171556D-08 )*Y+ &
          7.72165684563049D-08 )*Y-3.59858901591047D-07 )*Y+ &
          2.43682618601000D-06
!     XXXX=20.0 TO 25.0                             NROOTS = 5
  2000 Y = XX-22.5D+00
      RT1 = (((((((((-1.13927848238726D-15*Y+7.39404133595713D-15)*Y+ &
          1.445982921243D-13)*Y-2.676703245252D-12)*Y+ &
          5.823521627177D-12)*Y+2.17264723874381D-10 )*Y+ &
          3.56242145897468D-09 )*Y-3.03763737404491D-07 )*Y+ &
          9.46859114120901D-06 )*Y-2.30896753853196D-04 )*Y+ &
          5.24663913001114D-03
      RT2 = (((((((((( 2.89872355524581D-16*Y-1.22296292045864D-14)*Y+ &
          6.184065097200D-14)*Y+1.649846591230D-12)*Y- &
          2.729713905266D-11)*Y+3.709913790650D-11)*Y+ &
          2.216486288382D-09)*Y+4.616160236414D-08)*Y- &
          3.32380270861364D-06 )*Y+9.84635072633776D-05 )*Y- &
          2.30092118015697D-03 )*Y+5.00845183695073D-02
      RT3 = (((((((((( 1.97068646590923D-15*Y-4.89419270626800D-14)*Y+ &
          1.136466605916D-13)*Y+7.546203883874D-12)*Y- &
          9.635646767455D-11)*Y-8.295965491209D-11)*Y+ &
          7.534109114453D-09)*Y+2.699970652707D-07)*Y- &
          1.42982334217081D-05 )*Y+3.78290946669264D-04 )*Y- &
          8.03133015084373D-03 )*Y+1.58689469640791D-01
      RT4 = (((((((((( 1.33642069941389D-14*Y-1.55850612605745D-13)*Y- &
          7.522712577474D-13)*Y+3.209520801187D-11)*Y- &
          2.075594313618D-10)*Y-2.070575894402D-09)*Y+ &
          7.323046997451D-09)*Y+1.851491550417D-06)*Y- &
          6.37524802411383D-05 )*Y+1.36795464918785D-03 )*Y- &
          2.42051126993146D-02 )*Y+3.97847167557815D-01
      RT5 = ((((((((((-6.07053986130526D-14*Y+1.04447493138843D-12)*Y- &
          4.286617818951D-13)*Y-2.632066100073D-10)*Y+ &
          4.804518986559D-09)*Y-1.835675889421D-08)*Y- &
          1.068175391334D-06)*Y+3.292234974141D-05)*Y- &
          5.94805357558251D-04 )*Y+8.29382168612791D-03 )*Y- &
          9.93122509049447D-02 )*Y+1.09857804755042D+00
      WW1 = (((((((((-9.10338640266542D-15*Y+1.00438927627833D-13)*Y+ &
          7.817349237071D-13)*Y-2.547619474232D-11)*Y+ &
          1.479321506529D-10)*Y+1.52314028857627D-09 )*Y+ &
          9.20072040917242D-09 )*Y-2.19427111221848D-06 )*Y+ &
          8.65797782880311D-05 )*Y-2.82718629312875D-03 )*Y+ &
          1.28718310443295D-01
      WW2 = ((((((((( 5.52380927618760D-15*Y-6.43424400204124D-14)*Y- &
          2.358734508092D-13)*Y+8.261326648131D-12)*Y+ &
          9.229645304956D-11)*Y-5.68108973828949D-09 )*Y+ &
          1.22477891136278D-07 )*Y-2.11919643127927D-06 )*Y+ &
          4.23605032368922D-05 )*Y-1.14423444576221D-03 )*Y+ &
          5.06607252890186D-02
      WW3 = ((((((((( 3.99457454087556D-15*Y-5.11826702824182D-14)*Y- &
          4.157593182747D-14)*Y+4.214670817758D-12)*Y+ &
          6.705582751532D-11)*Y-3.36086411698418D-09 )*Y+ &
          6.07453633298986D-08 )*Y-7.40736211041247D-07 )*Y+ &
          8.84176371665149D-06 )*Y-1.72559275066834D-04 )*Y+ &
          7.16639814253567D-03
      WW4 = (((((((((((-2.14649508112234D-18*Y-2.45525846412281D-18)*Y+ &
          6.126212599772D-16)*Y-8.526651626939D-15)*Y+ &
          4.826636065733D-14)*Y-3.39554163649740D-13 )*Y+ &
          1.67070784862985D-11 )*Y-4.42671979311163D-10 )*Y+ &
          6.77368055908400D-09 )*Y-7.03520999708859D-08 )*Y+ &
          6.04993294708874D-07 )*Y-7.80555094280483D-06 )*Y+ &
          2.85954806605017D-04
      WW5 = ((((((((((((-5.63938733073804D-21*Y+6.92182516324628D-20)*Y- &
          1.586937691507D-18)*Y+3.357639744582D-17)*Y- &
          4.810285046442D-16)*Y+5.386312669975D-15)*Y- &
          6.117895297439D-14)*Y+8.441808227634D-13)*Y- &
          1.18527596836592D-11 )*Y+1.36296870441445D-10 )*Y- &
          1.17842611094141D-09 )*Y+7.80430641995926D-09 )*Y- &
          5.97767417400540D-08 )*Y+1.65186146094969D-06
! ! !       !RETUR
          GO TO 7777
  2200 WW1 = SQRT(PIE4/XX)
      !write(*,*) "PIE4",PIE4
      !write(*,*) "XX",XX
      !write(*,*) "sqrt stuff",SQRT(PIE4/XX)
      IF (XX .GT. 40.0D+00) GO TO 2400
!     XXXX=25.0 TO 40.0                             NROOTS = 5
      !write(*,*) "R15",R15
      E = EXP(-XX)
      !write(*,*) "E",E
      !write(*,*) "WW1",WW1
      !WW1 is problematic
      RT1 = ((((((((-1.73363958895356D-06*XX+1.19921331441483D-04)*XX - &
          1.59437614121125D-02)*XX+1.13467897349442D+00)*XX - &
          4.47216460864586D+01)*XX+1.06251216612604D+03)*XX - &
          1.52073917378512D+04)*XX+1.20662887111273D+05)*XX - &
          4.07186366852475D+05)*E + R15/(XX-R15)
      RT2 = ((((((((-1.60102542621710D-05*XX+1.10331262112395D-03)*XX - &
          1.50043662589017D-01)*XX+1.05563640866077D+01)*XX - &
          4.10468817024806D+02)*XX+9.62604416506819D+03)*XX - &
          1.35888069838270D+05)*XX+1.06107577038340D+06)*XX - &
          3.51190792816119D+06)*E + R25/(XX-R25)
      RT3 = ((((((((-4.48880032128422D-05*XX+2.69025112122177D-03)*XX - &
          4.01048115525954D-01)*XX+2.78360021977405D+01)*XX - &
          1.04891729356965D+03)*XX+2.36985942687423D+04)*XX - &
          3.19504627257548D+05)*XX+2.34879693563358D+06)*XX - &
          7.16341568174085D+06)*E + R35/(XX-R35)
      RT4 = ((((((((-6.38526371092582D-05*XX-2.29263585792626D-03)*XX - &
          7.65735935499627D-02)*XX+9.12692349152792D+00)*XX - &
          2.32077034386717D+02)*XX+2.81839578728845D+02)*XX + &
          9.59529683876419D+04)*XX-1.77638956809518D+06)*XX + &
          1.02489759645410D+07)*E + R45/(XX-R45)
      RT5 = ((((((((-3.59049364231569D-05*XX-2.25963977930044D-02)*XX + &
          1.12594870794668D+00)*XX-4.56752462103909D+01)*XX + &
          1.05804526830637D+03)*XX-1.16003199605875D+04)*XX - &
          4.07297627297272D+04)*XX+2.22215528319857D+06)*XX - &
          1.61196455032613D+07)*E + R55/(XX-R55)
      WW5 = (((((((((-4.61100906133970D-10*XX+1.43069932644286D-07)*XX - &
          1.63960915431080D-05)*XX+1.15791154612838D-03)*XX - &
          5.30573476742071D-02)*XX+1.61156533367153D+00)*XX - &
          3.23248143316007D+01)*XX+4.12007318109157D+02)*XX - &
          3.02260070158372D+03)*XX+9.71575094154768D+03)*E + W55*WW1
      WW4 = (((((((((-2.40799435809950D-08*XX+8.12621667601546D-06)*XX - &
          9.04491430884113D-04)*XX+6.37686375770059D-02)*XX - &
          2.96135703135647D+00)*XX+9.15142356996330D+01)*XX - &
          1.86971865249111D+03)*XX+2.42945528916947D+04)*XX - &
          1.81852473229081D+05)*XX+5.96854758661427D+05)*E + W45*WW1
      WW3 = (((((((( 1.83574464457207D-05*XX-1.54837969489927D-03)*XX + &
          1.18520453711586D-01)*XX-6.69649981309161D+00)*XX + &
          2.44789386487321D+02)*XX-5.68832664556359D+03)*XX + &
          8.14507604229357D+04)*XX-6.55181056671474D+05)*XX + &
          2.26410896607237D+06)*E + W35*WW1
      WW2 = (((((((( 2.77778345870650D-05*XX-2.22835017655890D-03)*XX + &
          1.61077633475573D-01)*XX-8.96743743396132D+00)*XX + &
          3.28062687293374D+02)*XX-7.65722701219557D+03)*XX + &
          1.10255055017664D+05)*XX-8.92528122219324D+05)*XX + &
          3.10638627744347D+06)*E + W25*WW1
      WW1 = WW1-0.01962D+00*E-WW2-WW3-WW4-WW5
      GO TO 7777

  2400 IF (XX .GT. 59.0D+00) GO TO 2600
!     X=40.0 TO 59.0                             NROOTS = 5
      XXX = XX**3
      E = XXX*EXP(-XX)
      RT1 = (((-2.43758528330205D-02*XX+2.07301567989771D+00)*XX -6.45964225381113D+01)*XX+7.14160088655470D+02)*E + R15/(XX-R15)
      RT2 = (((-2.28861955413636D-01*XX+1.93190784733691D+01)*XX -5.99774730340912D+02)*XX+6.61844165304871D+03)*E + R25/(XX-R25)
      RT3 = (((-6.95053039285586D-01*XX+5.76874090316016D+01)*XX -1.77704143225520D+03)*XX+1.95366082947811D+04)*E + R35/(XX-R35)
      RT4 = (((-1.58072809087018D+00*XX+1.27050801091948D+02)*XX -3.86687350914280D+03)*XX+4.23024828121420D+04)*E + R45/(XX-R45)
      RT5 = (((-3.33963830405396D+00*XX+2.51830424600204D+02)*XX -7.57728527654961D+03)*XX+8.21966816595690D+04)*E + R55/(XX-R55)
      E = XXX*E
      WW5 = (( 1.35482430510942D-08*XX-3.27722199212781D-07)*XX +2.41522703684296D-06)*E + W55*WW1
      WW4 = (( 1.23464092261605D-06*XX-3.55224564275590D-05)*XX +3.03274662192286D-04)*E + W45*WW1
      WW3 = (( 1.34547929260279D-05*XX-4.19389884772726D-04)*XX +3.87706687610809D-03)*E + W35*WW1
      WW2 = (( 2.09539509123135D-05*XX-6.87646614786982D-04)*XX +6.68743788585688D-03)*E + W25*WW1
      WW1 = WW1-WW2-WW3-WW4-WW5
      GO TO 7777

  2600 RT1 = R15/(XX-R15)
      RT2 = R25/(XX-R25)
      RT3 = R35/(XX-R35)
      RT4 = R45/(XX-R45)
      RT5 = R55/(XX-R55)
      WW2 = W25*WW1
      WW3 = W35*WW1
      WW4 = W45*WW1
      WW5 = W55*WW1
      WW1 = WW1-WW2-WW3-WW4-WW5
      !GO TO 7777

7777 CONTINUE

                   U(1)=RT1
                   !write(*,*) "U1",U(1)
                   U(2)=RT2
                   U(3)=RT3
                   U(4)=RT4
                   U(5)=RT5
                   W(1)=WW1
                   !write(*,*) "W1",W(1)
                   W(2)=WW2
                   W(3)=WW3
                   W(4)=WW4
                   W(5)=WW5

                MM = 0
                DO M = 1,5
                   !U2 = U(M)*RHO
                   U2 = U(M)*((A*B)/(A+B)) !maybe bring back rho?
                   F00 = EXPE*W(M)
                    DO I = 1,5
                       IN(I) = IN1(I)+MM
                    ENDDO

                     DUMINV = ONE/((A*B)+U2*(A+B))
                     DM2INV = HALF*DUMINV
                     BP01 = (A+U2)*DM2INV
                     B00 = U2*DM2INV
                     B10 = (B+U2)*DM2INV
                     XCP00 = (U2*C1X+C2X)*DUMINV
                     XC00 = (U2*C3X+C4X)*DUMINV
                     YCP00 = (U2*C1Y+C2Y)*DUMINV
                     YC00 = (U2*C3Y+C4Y)*DUMINV
                     ZCP00 = (U2*C1Z+C2Z)*DUMINV
                     ZC00 = (U2*C3Z+C4Z)*DUMINV

! !======================exchange integrals==============
! !     ----- I(0,0) -----
!       I1 = IN(1)

!       XIN(I1,ID) = ONE
!       YIN(I1,ID) = ONE
!       ZIN(I1,ID) = F00      

!       I2 = IN(2)
!       K2 = KN(2)
!       CP10 = B00

!     ! ----- I(1,0) -----
!       XIN(I2,ID) = XC00
!       YIN(I2,ID) = YC00
!       ZIN(I2,ID) = ZC00*F00

! !     ----- I(0,1) -----
!       I3 = I1+K2

!       XIN(I3,ID) = XCP00
!       YIN(I3,ID) = YCP00
!       ZIN(I3,ID) = ZCP00*F00

! !     ----- I(1,1) -----
!       I3 = I2+K2

!       XIN(I3,ID) = XCP00*XIN(I2,ID)+CP10
!       YIN(I3,ID) = YCP00*YIN(I2,ID)+CP10
!       ZIN(I3,ID) = ZCP00*ZIN(I2,ID)+CP10*F00

!         C10 = ZERO
!         I3 = I1
!         I4 = I2
!         DO N = 2,4
!           C10 = C10+B10

! !     ----- I(N,0) -----

!           I5 = IN(N+1)

!           XIN(I5,ID) = C10*XIN(I3,ID)+XC00*XIN(I4,ID)
!           YIN(I5,ID) = C10*YIN(I3,ID)+YC00*YIN(I4,ID)
!           ZIN(I5,ID) = C10*ZIN(I3,ID)+ZC00*ZIN(I4,ID)

!             CP10 = CP10+B00

! !     ----- I(N,1) -----

!             I3 = I5+K2

!             XIN(I3,ID) = XCP00*XIN(I5,ID)+CP10*XIN(I4,ID)
!             YIN(I3,ID) = YCP00*YIN(I5,ID)+CP10*YIN(I4,ID)
!             ZIN(I3,ID) = ZCP00*ZIN(I5,ID)+CP10*ZIN(I4,ID)

!           I3 = I4
!           I4 = I5
!       ENDDO

!         CP01 = ZERO
!         C01 = B00
!         I3 = I1
!         I4 = I1+K2
!         DO MMM = 2,4
!           CP01 = CP01+BP01

! !     ----- I(0,M) -----
!           I5 = I1+KN(MMM+1)

!           XIN(I5,ID) = CP01*XIN(I3,ID)+XCP00*XIN(I4,ID)
!           YIN(I5,ID) = CP01*YIN(I3,ID)+YCP00*YIN(I4,ID)
!           ZIN(I5,ID) = CP01*ZIN(I3,ID)+ZCP00*ZIN(I4,ID)
! !     ----- I(1,M) -----
!             C01 = C01+B00
!             I3 = I2+KN(MMM+1)

!             XIN(I3,ID) = XC00*XIN(I5,ID)+C01*XIN(I4,ID)
!             YIN(I3,ID) = YC00*YIN(I5,ID)+C01*YIN(I4,ID)
!             ZIN(I3,ID) = ZC00*ZIN(I5,ID)+C01*ZIN(I4,ID)

!           I3 = I4
!           I4 = I5
!         ENDDO

! !     ----- I(N,M) -----
!         C01 = B00
!         K3 = K2
!         DO MMM = 2,4
!           K4 = KN(MMM+1)
!           C01 = C01+B00
!           I3 = I1
!           I4 = I2
!           C10 = B10
!           DO N = 2,4
!             I5 = IN(N+1)

!             XIN(I5+K4,ID) = C10*XIN(I3+K4,ID)+XC00*XIN(I4+K4,ID)+C01*XIN(I4+K3,ID)
!             YIN(I5+K4,ID) = C10*YIN(I3+K4,ID)+YC00*YIN(I4+K4,ID)+C01*YIN(I4+K3,ID)
!             ZIN(I5+K4,ID) = C10*ZIN(I3+K4,ID)+ZC00*ZIN(I4+K4,ID)+C01*ZIN(I4+K3,ID)

!             C10 = C10+B10
!             I3 = I4
!             I4 = I5
!           ENDDO
!           K3 = K4
!         ENDDO

! !     ----- I(NI,NJ,M) -----
!         MMM = 0
!         I5 = IN(5)
!         FIRST1 = .TRUE.
!         DO 4300 WHILE (FIRST1 .OR. MMM .LE. 4)
!           MINN = 2
!           KM = KN(MMM+1)
!           FIRST2 = .TRUE.
!           DO 3600 WHILE (FIRST2 .OR. MINN .LT. 4)
!             N = 4
!             I3 = I5+KM
!             FIRST3 = .TRUE.
!             DO 3400 WHILE (FIRST3 .OR. N .GT. MINN)
!               I4 = IN(N)+KM

!               XIN(I3,ID) = XIN(I3,ID)+DXIJ*XIN(I4,ID)
!               YIN(I3,ID) = YIN(I3,ID)+DYIJ*YIN(I4,ID)
!               ZIN(I3,ID) = ZIN(I3,ID)+DZIJ*ZIN(I4,ID)

!               I3 = I4
!               N = N-1
!               FIRST3 = .FALSE.
!   3400       END DO
!             MINN = MINN+1
!             FIRST2 = .FALSE.
!   3600     END DO

!             I3 = 49+KM+I1
!             DO 4000 NJ = 1,2
!               I4 = I3
!               DO 3800 NI = 1,2

!                 XIN(I4,ID) = XIN(I4+294,ID)+DXIJ*XIN(I4-49,ID)
!                 YIN(I4,ID) = YIN(I4+294,ID)+DYIJ*YIN(I4-49,ID)
!                 ZIN(I4,ID) = ZIN(I4+294,ID)+DZIJ*ZIN(I4-49,ID)
!                 I4 = I4+343
!   3800         CONTINUE
!               I3 = I3+49
!   4000       CONTINUE

!           MMM = MMM+1
!           FIRST1 = .FALSE.
!   4300   END DO

! !     ----- I(NI,NJ,NK,NL) -----

!         I5 = KN(5)
!         IAA = I1
!         NI = 0
!         FIRST4 = .TRUE.
!         DO 5800 WHILE (FIRST4 .OR. NI .LE. 2)
!           NJ = 0
!           IBB = IAA
!           FIRST1 = .TRUE.
!           DO 5700 WHILE (FIRST1 .OR. NJ .LE. 2)
!             MINN = 2
!             FIRST2 = .TRUE.
!             DO 5300 WHILE (FIRST2 .OR. MINN .LT. 4)
!               MMM = 4
!               I3 = IBB+I5
!               FIRST3 = .TRUE.
!               DO 5200 WHILE (FIRST3 .OR. MMM .GT. MINN)
!                 I4 = IBB+KN(MMM)

!                 XIN(I3,ID) = XIN(I3,ID)+DXKL*XIN(I4,ID)
!                 YIN(I3,ID) = YIN(I3,ID)+DYKL*YIN(I4,ID)
!                 ZIN(I3,ID) = ZIN(I3,ID)+DZKL*ZIN(I4,ID)

!                 I3 = I4
!                 MMM = MMM-1
!                 FIRST3 = .FALSE.
!   5200         END DO
!               MINN = MINN+1
!               FIRST2 = .FALSE.
!   5300       END DO
!               I3 = IBB+1
!               DO NL = 1,2
!                 I4 = I3
!                 DO NK = 1,2

!                   XIN(I4,ID) = XIN(I4+6,ID)+DXKL*XIN(I4-1,ID)
!                   YIN(I4,ID) = YIN(I4+6,ID)+DYKL*YIN(I4-1,ID)
!                   ZIN(I4,ID) = ZIN(I4+6,ID)+DZKL*ZIN(I4-1,ID)

!                   I4 = I4+7
!                 ENDDO
!               I3 = I3+1
!               ENDDO
            
!             NJ = NJ+1
!             IBB = IBB+49
!             FIRST1 = .FALSE.
!   5700     END DO
!           NI = NI+1
!           IAA = IAA+343
!           FIRST4 = .FALSE.
!   5800   END DO
!     ----- I(0,0) -----
      I1 = IN(1)

      XIN(I1) = ONE
      YIN(I1) = ONE
      ZIN(I1) = F00      

      I2 = IN(2)
      K2 = KN(2)
      CP10 = B00

    ! ----- I(1,0) -----
      !write(*,*) "size of XIN array", SIZEOF(XIN(I2))
      XIN(I2) = XC00
      YIN(I2) = YC00
      ZIN(I2) = ZC00*F00

!     ----- I(0,1) -----
      I3 = I1+K2

      XIN(I3) = XCP00
      YIN(I3) = YCP00
      ZIN(I3) = ZCP00*F00

!     ----- I(1,1) -----
      I3 = I2+K2

      XIN(I3) = XCP00*XIN(I2)+CP10
      YIN(I3) = YCP00*YIN(I2)+CP10
      ZIN(I3) = ZCP00*ZIN(I2)+CP10*F00

        C10 = ZERO
        I3 = I1
        I4 = I2
        DO N = 2,4
          C10 = C10+B10

!     ----- I(N,0) -----

          I5 = IN(N+1)

          XIN(I5) = C10*XIN(I3)+XC00*XIN(I4)
          YIN(I5) = C10*YIN(I3)+YC00*YIN(I4)
          ZIN(I5) = C10*ZIN(I3)+ZC00*ZIN(I4)

            CP10 = CP10+B00

!     ----- I(N,1) -----

            I3 = I5+K2

            XIN(I3) = XCP00*XIN(I5)+CP10*XIN(I4)
            YIN(I3) = YCP00*YIN(I5)+CP10*YIN(I4)
            ZIN(I3) = ZCP00*ZIN(I5)+CP10*ZIN(I4)

          I3 = I4
          I4 = I5
      ENDDO

        CP01 = ZERO
        C01 = B00
        I3 = I1
        I4 = I1+K2
        DO MMM = 2,4
          CP01 = CP01+BP01

!     ----- I(0,M) -----
          I5 = I1+KN(MMM+1)

          XIN(I5) = CP01*XIN(I3)+XCP00*XIN(I4)
          YIN(I5) = CP01*YIN(I3)+YCP00*YIN(I4)
          ZIN(I5) = CP01*ZIN(I3)+ZCP00*ZIN(I4)
!     ----- I(1,M) -----
            C01 = C01+B00
            I3 = I2+KN(MMM+1)

            XIN(I3) = XC00*XIN(I5)+C01*XIN(I4)
            YIN(I3) = YC00*YIN(I5)+C01*YIN(I4)
            ZIN(I3) = ZC00*ZIN(I5)+C01*ZIN(I4)

          I3 = I4
          I4 = I5
        ENDDO

!     ----- I(N,M) -----
        C01 = B00
        K3 = K2
        DO MMM = 2,4
          K4 = KN(MMM+1)
          C01 = C01+B00
          I3 = I1
          I4 = I2
          C10 = B10
          DO N = 2,4
            I5 = IN(N+1)

            XIN(I5+K4) = C10*XIN(I3+K4)+XC00*XIN(I4+K4)+C01*XIN(I4+K3)
            YIN(I5+K4) = C10*YIN(I3+K4)+YC00*YIN(I4+K4)+C01*YIN(I4+K3)
            ZIN(I5+K4) = C10*ZIN(I3+K4)+ZC00*ZIN(I4+K4)+C01*ZIN(I4+K3)

            C10 = C10+B10
            I3 = I4
            I4 = I5
          ENDDO
          K3 = K4
        ENDDO

!     ----- I(NI,NJ,M) -----
        MMM = 0
        I5 = IN(5)
        FIRST1 = .TRUE.
        DO 4300 WHILE (FIRST1 .OR. MMM .LE. 4)
          MINN = 2
          KM = KN(MMM+1)
          FIRST2 = .TRUE.
          DO 3600 WHILE (FIRST2 .OR. MINN .LT. 4)
            N = 4
            I3 = I5+KM
            FIRST3 = .TRUE.
            DO 3400 WHILE (FIRST3 .OR. N .GT. MINN)
              I4 = IN(N)+KM

              XIN(I3) = XIN(I3)+DXIJ*XIN(I4)
              YIN(I3) = YIN(I3)+DYIJ*YIN(I4)
              ZIN(I3) = ZIN(I3)+DZIJ*ZIN(I4)

              I3 = I4
              N = N-1
              FIRST3 = .FALSE.
  3400       END DO
            MINN = MINN+1
            FIRST2 = .FALSE.
  3600     END DO

            I3 = 49+KM+I1
            DO 4000 NJ = 1,2
              I4 = I3
              DO 3800 NI = 1,2

                XIN(I4) = XIN(I4+294)+DXIJ*XIN(I4-49)
                YIN(I4) = YIN(I4+294)+DYIJ*YIN(I4-49)
                ZIN(I4) = ZIN(I4+294)+DZIJ*ZIN(I4-49)
                I4 = I4+343
  3800         CONTINUE
              I3 = I3+49
  4000       CONTINUE

          MMM = MMM+1
          FIRST1 = .FALSE.
  4300   END DO

!     ----- I(NI,NJ,NK,NL) -----

        I5 = KN(5)
        IAA = I1
        NI = 0
        FIRST4 = .TRUE.
        DO 5800 WHILE (FIRST4 .OR. NI .LE. 2)
          NJ = 0
          IBB = IAA
          FIRST1 = .TRUE.
          DO 5700 WHILE (FIRST1 .OR. NJ .LE. 2)
            MINN = 2
            FIRST2 = .TRUE.
            DO 5300 WHILE (FIRST2 .OR. MINN .LT. 4)
              MMM = 4
              I3 = IBB+I5
              FIRST3 = .TRUE.
              DO 5200 WHILE (FIRST3 .OR. MMM .GT. MINN)
                I4 = IBB+KN(MMM)

                XIN(I3) = XIN(I3)+DXKL*XIN(I4)
                YIN(I3) = YIN(I3)+DYKL*YIN(I4)
                ZIN(I3) = ZIN(I3)+DZKL*ZIN(I4)

                I3 = I4
                MMM = MMM-1
                FIRST3 = .FALSE.
  5200         END DO
              MINN = MINN+1
              FIRST2 = .FALSE.
  5300       END DO
              I3 = IBB+1
              DO NL = 1,2
                I4 = I3
                DO NK = 1,2

                  XIN(I4) = XIN(I4+6)+DXKL*XIN(I4-1)
                  YIN(I4) = YIN(I4+6)+DYKL*YIN(I4-1)
                  ZIN(I4) = ZIN(I4+6)+DZKL*ZIN(I4-1)

                  I4 = I4+7
                ENDDO
              I3 = I3+1
              ENDDO
            
            NJ = NJ+1
            IBB = IBB+49
            FIRST1 = .FALSE.
  5700     END DO
          NI = NI+1
          IAA = IAA+343
          FIRST4 = .FALSE.
  5800   END DO


                   !MM = MM+2401
                   MM = MM+1000
                   !MM = MM+800
                   !MM = MM+410
                   !what is 2401????
                   !write(*,*) "MM is ",MM
                ENDDO

      DO I = 1,IJ
      D1 = DIJ(I)
      NX = IJX(I)
      !write(*,*) "NX",NX
      NY = IJY(I)
      NZ = IJZ(I)
      N1 = IJGT(I)
      MAX = IK(I)
      DO K = 1,MAX
      MX = NX+KLX(K)
      MY = NY+KLY(K)
      MZ = NZ+KLZ(K)
      ! write(*,*) "MX",MX
      NNN = N1+KLGT(K)
      !GHONDO(NNN) = GHONDO(NNN) + DKL(K)*D1*( XIN(MX)*YIN(MY)*ZIN(MZ)+ XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401) &      
      !      + XIN(MX+ 4802)*YIN(MY+ 4802)*ZIN(MZ+ 4802)+ XIN(MX+ 7203)*YIN(MY+ 7203)*ZIN(MZ+ 7203)+XIN(MX+ 9604)*YIN(MY+ 9604)*ZIN(MZ+ 9604))
      GHONDO(NNN) = GHONDO(NNN) + DKL(K)*D1*( XIN(MX)*YIN(MY)*ZIN(MZ)+ XIN(MX+ 1000)*YIN(MY+ 1000)*ZIN(MZ+ 1000) &      
            + XIN(MX+ 2000)*YIN(MY+ 2000)*ZIN(MZ+ 2000)+ XIN(MX+ 3000)*YIN(MY+ 3000)*ZIN(MZ+ 3000)+XIN(MX+ 4000)*YIN(MY+ 4000)*ZIN(MZ+ 4000))
      !GHONDO(NNN) = GHONDO(NNN) + DKL(K)*D1*( XIN(MX)*YIN(MY)*ZIN(MZ)+ XIN(MX+ 800)*YIN(MY+ 800)*ZIN(MZ+ 800) &      
      !      + XIN(MX+ 1600)*YIN(MY+ 1600)*ZIN(MZ+ 1600)+ XIN(MX+ 2400)*YIN(MY+ 2400)*ZIN(MZ+ 2400)+XIN(MX+ 3200)*YIN(MY+ 3200)*ZIN(MZ+ 3200))
      !GHONDO(NNN) = GHONDO(NNN) + DKL(K)
      ENDDO
      ENDDO
!================ end calculation of integrals ================

!==================== digestion of integrals ================
      XVAL1 = DFTTYP(3)
      !write(*,*) "XVAL1", XVAL1
      XVAL4 = FOUR

      NIJ = 0
      MAXJ2 = MAXJ
      I_INDEX = 1
      DO I = MINI,MAXI
         IF (IANDJ) MAXJ2 = I

         I1 = I+LOCI

         IJ_INDEX = I_INDEX
         !write(*,*) "IJ_INDEX", IJ_INDEX
         I_INDEX = I_INDEX + NGTI

         DO J = MINJ,MAXJ2
            NIJ = NIJ+1
            MAXL2 = MAXL

            J1 = J+LOCJ
            I2 = I1
            J2 = J1
            IF (I1.LT.J1) THEN ! SORT <IJ|
               I2 = J1
               J2 = I1
            ENDIF

            IJK_INDEX = IJ_INDEX
            IJ_INDEX = IJ_INDEX + NGTJ

            NKL = NIJ

            DO K =  MINK,MAXK
               IF (KANDL) MAXL2 = K

               K1 = K + LOCK

               IF(SAME) THEN ! ACCOUNT FOR NON-UNIQUE PERMUTATIONS
                  ITMP = MIN(MAXL2-MINL+1,NKL)
                  IF (ITMP.EQ.0) CYCLE
                  MAXL2 = MINL + ITMP - 1
                  NKL = NKL - ITMP
               ENDIF

               IJKL_INDEX = IJK_INDEX
               IJK_INDEX = IJK_INDEX + NGTK

               DO L=MINL,MAXL2

                  VAL = GHONDO( IJKL_INDEX )
                  !IF(VAL .NE. 0.00) write(*,*) "val is", VAL
                  !write(*,*) "val is", VAL
                  IJKL_INDEX = IJKL_INDEX + NGTL
                  !IF(ABS(VAL).LT.CUTOFF) CYCLE ! GOTO 300
                  NINT = NINT + 1
                  !write(*,*) "NINT", NINT

                  L1 = L + LOCL
                  K2 = K1
                  L2 = L1

                  IF (K2.LT.L2) THEN ! SORT |KL>
                     K2 = L1
                     L2 = K1
                  ENDIF

                  II = I2
                  JJ = J2
                  KK = K2
                  LL = L2

                  IF (II.LT.KK) THEN ! SORT <IJ|KL>
                     II = K2
                     JJ = L2
                     KK = I2
                     LL = J2
                  ELSE IF (II.EQ.KK.AND.JJ.LT.LL) THEN ! SORT <IJ|IL>
                     JJ = L2
                     LL = J2
                  ENDIF

                  II2 = IA(II)
                  JJ2 = IA(JJ)
                  KK2 = IA(KK)

                  IJ = II2 + JJ
                  IKK = II2 + KK
                  IL = II2 + LL
                  JK = JJ2 + KK
                  JL = JJ2 + LL
                  KL = KK2 + LL
                  IF (JJ.LT.KK) JK = KK2 + JJ
                  IF (JJ.LT.LL) JL = IA(LL) + JJ
!
!       ACCOUNT FOR IDENTICAL PERMUTATIONS.
!
                  IF(II.EQ.JJ) VAL = VAL*HALF
                  IF(KK.EQ.LL) VAL = VAL*HALF
                  IF(II.EQ.KK.AND.JJ.EQ.LL) VAL = VAL*HALF
                  VAL1 = VAL*XVAL1
                  !write(*,*) "VAL1", VAL1
                  VAL4 = VAL*XVAL4
                  ! write(*,*) "VAL4", VAL4
                  ! write(*,*) "DA(KL)", DA(KL)
!$omp atomic
                  FA(IJ) = FA(IJ) + VAL4*DA(KL)
                    !write(*,*) "IJ", IJ
!$omp atomic
                  FA(KL) = FA(KL) + VAL4*DA(IJ)
!$omp atomic
                  FA(IKK) = FA(IKK) - VAL1*DA(JL)
!$omp atomic
                  FA(JL) = FA(JL) - VAL1*DA(IKk)
!$omp atomic
                  FA(IL) = FA(IL) - VAL1*DA(JK)
!$omp atomic
                  FA(JK) = FA(JK) - VAL1*DA(IL)
               ENDDO
            ENDDO
         ENDDO
      ENDDO

   ENDDO
!$omp end target teams distribute parallel do
  ENDDO

  END SUBROUTINE rysint_gpu_5

SUBROUTINE rysint_gpu_4(ii,jj,kk,ll,ktype,ghondo,ddij,ngth,c,&
                        kstart,katom,kng,kloc,kmin,kmax,&
                        ex,cs,cp,cd,cf,cg,ch,ci,&
                        NGTI,NGTJ,NGTK,NGTL,&
                        NIJ,IJ,KL,IJKL,&
                        qq4,&
                        AX,AY,AZ,BX,BY,BZ,RAB,CX,CY,CZ,DX,DY,DZ,RCD1,&
                        dfttyp,&
                        ia,da,fa,nint,&
                        l1,l2a,nflmat,&
                        IDX_RYS_4,N_RYS_4,NCP_rys)

    USE omp_lib
    USE prec, ONLY: fp
    USE mx_limits, ONLY: &
    mxsh,mxg2,mxatm,mxgtot,mxgsh

    !INTEGER, INTENT(IN) :: &
    !l1, l2a,nflmat, ia(l1)

    !for dirfck
    REAL(KIND=fp),INTENT(OUT) :: fa(l2a*nflmat) !FOCK MATRIX
    REAL(KIND=fp),INTENT(INOUT) :: da(l2a) !DENSITY MATRIX
    INTEGER,INTENT(IN) :: ia(l1) !TRIANGULAR MATRIX 
    INTEGER :: l2a,nflmat,l1,nint 

    REAL(KIND=fp),DIMENSION(20) :: DFTTYP

    INTEGER, INTENT(IN) :: ktype(mxsh)
    INTEGER :: ii,jj,kk,ll

    REAL(KIND=fp),INTENT(INOUT),ALLOCATABLE :: &
      ghondo(:)

    REAL(KIND=fp),INTENT(INOUT) :: &
      ddij(:)

    REAL(KIND=fp),DIMENSION(784) :: DKL,DIJ
    REAL(KIND=fp),DIMENSION(784) :: &
     IJGT,IJX,IJY,IJZ,IK,KLGT,KLX,KLY,KLZ

    REAL(KIND=fp),DIMENSION(MXG2) :: A,R,X1,Y1,Z1
    REAL(KIND=fp),DIMENSION(784) :: IJD

      !for shlexc
      INTEGER,DIMENSION(4) :: NGTH
      !for infoa
      REAL(KIND=fp) :: C(3,mxatm)
      !for roots
      REAL(KIND=fp),DIMENSION(13) :: U,W
      REAL(KIND=fp) :: XX
      !for eriout
      INTEGER :: NGTI,NGTJ,NGTK,NGTL
      !shlnos
      INTEGER ::mini,maxi,minj,maxj,mink,maxk,minl,maxl,&
                LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,&
                NIJ,IJ,KL,IJKL
      REAL(KIND=fp) :: qq4  
      !nshell
      INTEGER,DIMENSION(mxsh) :: kstart,katom,kng,kloc,kmin,kmax,nshell
      REAL(KIND=fp),DIMENSION(mxgtot) :: ex,cs,cp,cd,cf,cg,ch,ci

      !shlinf
      double precision :: ga(mxgsh),csa1(mxgsh),cpa1(mxgsh),cda(mxgsh)
      double precision :: cfa(mxgsh),cga(mxgsh),cha(mxgsh),cia(mxgsh)
      double precision :: gb(mxgsh),csb1(mxgsh),cpb1(mxgsh),cdb(mxgsh)
      double precision :: cfb(mxgsh),cgb(mxgsh),chb(mxgsh),cib(mxgsh)
      double precision :: gc(mxgsh),csc1(mxgsh),cpc1(mxgsh),cdc(mxgsh)
      double precision :: cfc(mxgsh),cgc(mxgsh),chc(mxgsh),cic(mxgsh)
      double precision :: gd(mxgsh),csd1(mxgsh),cpd1(mxgsh),cdd(mxgsh)
      double precision :: cfd(mxgsh),cgd(mxgsh),chd(mxgsh),cid(mxgsh)

      double precision :: AX,AY,AZ,BX,BY,BZ,RAB,CX,CY,CZ,DX,DY,DZ,RCD1
      integer :: NGA1,NGB1,NGC1,NGD1
      LOGICAL IANDJ,KANDL,SAME
      integer :: IS1,JS1,KS1,LS1
      !IS1=kstart(ii)
      !rys
      INTEGER,PARAMETER :: NUM_QUART=9169597 !number of quartets
      INTEGER :: IDX_RYS_4(4,NUM_QUART,NCP_rys)
      INTEGER :: N_RYS_4(NCP_rys)

      INTEGER :: NSA_rys,ISA_rys !sum angular momentum
      INTEGER :: NCP_rys,ICP_rys !product contraction
      !INTEGER,ALLOCATABLE :: RYSINDEX(:,:,:,:) !rys shell index
      !INTEGER,ALLOCATABLE :: RYSTYP(:,:) !rystype
      INTEGER :: MAXRYS
      iNTEGER :: MM_rys
      integer :: NROOTS

      
      NSA_rys=9
      ISA_rys=1
      NROOTS=4
      !DO ISA_rys=1,NSA_rys
      !write(*,*) "NROOTS in 4", NROOTS
      DO ICP_rys=1,NCP_rys

      DO MM_rys=1,N_RYS_4(ICP_rys)

      ii=IDX_RYS_4(1,MM_rys,ICP_rys)
      jj=IDX_RYS_4(2,MM_rys,ICP_rys)
      kk=IDX_RYS_4(3,MM_rys,ICP_rys)
      ll=IDX_RYS_4(4,MM_rys,ICP_rys)

      MINI = kmin(ii)
      MAXI = kmax(ii)
      MINJ = kmin(jj)
      MAXJ = kmax(jj)
      MINK = kmin(kk)
      MAXK = kmax(kk)
      MINL = kmin(ll)
      MAXL = kmax(ll)

      IS1= kstart(ii)
      JS1= kstart(jj)
      KS1= kstart(kk)
      LS1= kstart(ll)

      lit= ktype(ii)-1
      ljt= ktype(jj)-1
      lkt= ktype(kk)-1
      llt= ktype(ll)-1

      NGA1= kng(ii)
      NGB1= kng(jj)
      NGC1= kng(kk)
      NGD1= kng(ll)

      LOCI = kloc(ii)-MINI
      LOCJ = kloc(jj)-MINJ
      LOCK = kloc(kk)-MINK
      LOCL = kloc(ll)-MINL

      ga=ex(is1)
      csa1=cs(is1)
      cpa1=cp(is1)
      cda=cd(is1)
      cfa=cf(is1)
      cga=cg(is1)
      cha=ch(is1)
      cia=ci(is1)

      gb=ex(js1)
      csb1=cs(js1)
      cpb1=cp(js1)
      cdb=cd(js1)
      cfb=cf(js1)
      cgb=cg(js1)
      chb=ch(js1)
      cib=ci(js1)

      gc=ex(ks1)
      csc1=cs(ks1)
      cpc1=cp(ks1)
      cdc=cd(ks1)
      cfc=cf(ks1)
      cgc=cg(ks1)
      chc=ch(ks1)
      cic=ci(ks1)

      gd=ex(ls1)
      csd1=cs(ls1)
      cpd1=cp(ls1)
      cdd=cd(ls1)
      cfd=cf(ls1)
      cgd=cg(ls1)
      chd=ch(ls1)
      cid=ci(ls1)

      IANDJ = II .EQ. JJ
      KANDL = KK .EQ. LL
      SAME = (II .EQ. KK) .AND. (JJ .EQ. LL)

!   fill indices for Rys quadrature
   CALL shells_gpu(1,ii,jj,kk,ll,.TRUE.,&
IJGT,IJX,IJY,IJZ,IK,KLGT,KLX,KLY,KLZ,&
IANDJ,KANDL,SAME,&
ii,jj,kk,ll,NGTI,NGTJ,NGTK,NGTL,&
LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,&
mini,maxi,minj,maxj,mink,maxk,minl,maxl,&
NIJ,IJ,KL,IJKL,&
NROOTS,NGTH,C,&
kstart,katom,ktype,kng,kloc,kmin,kmax,nshell,&
ex,cs,cp,cd,cf,cg,ch,ci,&
ga,csa1,cpa1,cda,cfa,cga,cha,cia,&
gb,csb1,cpb1,cdb,cfb,cgb,chb,cib,&
gc,csc1,cpc1,cdc,cfc,cgc,chc,cic,&
gd,csd1,cpd1,cdd,cfd,cgd,chd,cid,&
AX,AY,AZ,BX,BY,BZ,RAB,CX,CY,CZ,DX,DY,DZ,RCD1,&
NGA1,NGB1,NGC1,NGD1)

   CALL shells_gpu(2,ii,jj,kk,ll,.TRUE.,&
IJGT,IJX,IJY,IJZ,IK,KLGT,KLX,KLY,KLZ,&
IANDJ,KANDL,SAME,&
ii,jj,kk,ll,NGTI,NGTJ,NGTK,NGTL,&
LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,&
mini,maxi,minj,maxj,mink,maxk,minl,maxl,&
NIJ,IJ,KL,IJKL,&
NROOTS,NGTH,C,&
kstart,katom,ktype,kng,kloc,kmin,kmax,nshell,&
ex,cs,cp,cd,cf,cg,ch,ci,&
ga,csa1,cpa1,cda,cfa,cga,cha,cia,&
gb,csb1,cpb1,cdb,cfb,cgb,chb,cib,&
gc,csc1,cpc1,cdc,cfc,cgc,chc,cic,&
gd,csd1,cpd1,cdd,cfd,cgd,chd,cid,&
AX,AY,AZ,BX,BY,BZ,RAB,CX,CY,CZ,DX,DY,DZ,RCD1,&
NGA1,NGB1,NGC1,NGD1)

   CALL ijprim_gpu(ddij,a,r,x1,y1,z1,ijd,&
IANDJ,KANDL,SAME,&
LIT,mini,maxi,minj,maxj,mink,maxk,minl,maxl,NIJ,&
ga,csa1,cpa1,cda,cfa,cga,cha,cia,&
gb,csb1,cpb1,cdb,cfb,cgb,chb,cib,&
gc,csc1,cpc1,cdc,cfc,cgc,chc,cic,&
gd,csd1,cpd1,cdd,cfd,cgd,chd,cid,&
AX,AY,AZ,BX,BY,BZ,RAB,CX,CY,CZ,DX,DY,DZ,RCD1,&
NGA1,NGB1,NGC1,NGD1)

   call zqout_gpu(ghondo,&
IJGT,IJX,IJY,IJZ,IK,KLGT,KLX,KLY,KLZ,&
IANDJ,KANDL,SAME,&
mini,maxi,minj,maxj,mink,maxk,minl,maxl)

! calculate integrals
     CALL genral_gpu_4(ghondo,ddij,dkl,dij,&
IJGT,IJX,IJY,IJZ,IK,KLGT,KLX,KLY,KLZ, &
A,R,X1,Y1,Z1,IJD,iandj,kandl,same,&
LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,&
qq4,mini,maxi,minj,maxj,mink,maxk,minl,maxl,&
NIJ,IJ,KL,IJKL,&
NROOTS,&
ga,csa1,cpa1,cda,cfa,cga,cha,cia,&
gb,csb1,cpb1,cdb,cfb,cgb,chb,cib,&
gc,csc1,cpc1,cdc,cfc,cgc,chc,cic,&
gd,csd1,cpd1,cdd,cfd,cgd,chd,cid,&
AX,AY,AZ,BX,BY,BZ,RAB,CX,CY,CZ,DX,DY,DZ,RCD1,&
NGA1,NGB1,NGC1,NGD1)

     CALL dirfck_gpu_rys(ia,da,fa,ghondo,nint,&
kstart,katom,kng,kloc,kmin,kmax,&
ii,jj,kk,ll,&
NGTI,NGTJ,NGTK,NGTL,dfttyp)

  ENDDO
  ENDDO


  END SUBROUTINE rysint_gpu_4


  SUBROUTINE rysint_gpu_2(ii,jj,kk,ll,ktype,ghondo,ddij,ngth,c,&
                        kstart,katom,kng,kloc,kmin,kmax,&
                        ex,cs,cp,cd,&
                        NGTI,NGTJ,NGTK,NGTL,&
                        NIJ,IJ,KL,IJKL,&
                        qq4,&
                        AX,AY,AZ,BX,BY,BZ,RAB,CX,CY,CZ,DX,DY,DZ,RCD1,&
                        dfttyp,&
                        ia,da,fa,nint,&
                        l1,l2a,nflmat,&
                        IDX_RYS_2,N_RYS_2,NCP_rys)

    USE omp_lib
    USE prec, ONLY: fp
    USE mx_limits, ONLY: &
    mxsh,mxg2,mxatm,mxgtot,mxgsh

    IMPLICIT DOUBLE PRECISION (A-H,O-Z)

    !INTEGER, INTENT(IN) :: &
    !l1, l2a,nflmat, ia(l1)

    !for dirfck
    REAL(KIND=fp),INTENT(OUT) :: fa(l2a*nflmat) !FOCK MATRIX
    REAL(KIND=fp),INTENT(INOUT) :: da(l2a) !DENSITY MATRIX
    INTEGER,INTENT(IN) :: ia(l1) !TRIANGULAR MATRIX 
    INTEGER :: l2a,nflmat,l1,nint 

    REAL(KIND=fp),DIMENSION(20) :: DFTTYP

    INTEGER, INTENT(IN) :: ktype(mxsh)
    INTEGER :: ii,jj,kk,ll,inu,jnu,knu,lnu

    REAL(KIND=fp),INTENT(INOUT),ALLOCATABLE :: &
      ghondo(:)

    REAL(KIND=fp),INTENT(INOUT) :: &
      ddij(:)

    REAL(KIND=fp),DIMENSION(784) :: DKL,DIJ
    REAL(KIND=fp),DIMENSION(784) :: &
     IJGT,IJX,IJY,IJZ,IK,KLGT,KLX,KLY,KLZ
    REAL(KIND=fp),DIMENSION(84) :: &
     IX,IY,IZ,JX,JY,JZ,KX,KY,KZ,LX,LY,LZ
     INTEGER:: IJ_count

    ! REAL(KIND=fp),DIMENSION(MXG2) :: A,R,X1,Y1,Z1
    REAL(KIND=fp),DIMENSION(784) :: IJD

      !for shlexc
      INTEGER,DIMENSION(4) :: NGTH
      !for infoa
      REAL(KIND=fp) :: C(3,mxatm)
      !for roots
      REAL(KIND=fp),DIMENSION(13) :: U,W
      REAL(KIND=fp) :: XX
      !for eriout
      INTEGER :: NGTI,NGTJ,NGTK,NGTL
      !shlnos
      INTEGER ::mini,maxi,minj,maxj,mink,maxk,minl,maxl,&
                LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,&
                NIJ,IJ,KL,IJKL, max,IKK
      REAL(KIND=fp) :: qq4

      !nshell
      INTEGER,DIMENSION(mxsh) :: kstart,katom,kng,kloc,kmin,kmax
      REAL(KIND=fp),DIMENSION(mxgtot) :: ex,cs,cp,cd,cf,cg,ch,ci

      !shlinf
      double precision :: ga(mxgsh),csa1(mxgsh),cpa1(mxgsh),cda(mxgsh)
      double precision :: gb(mxgsh),csb1(mxgsh),cpb1(mxgsh),cdb(mxgsh)
      double precision :: gc(mxgsh),csc1(mxgsh),cpc1(mxgsh),cdc(mxgsh)
      double precision :: gd(mxgsh),csd1(mxgsh),cpd1(mxgsh),cdd(mxgsh)

      double precision :: AX,AY,AZ,BX,BY,BZ,RAB,CX,CY,CZ,DX,DY,DZ,RCD1
      integer :: NGA1,NGB1,NGC1,NGD1
      LOGICAL IANDJ,KANDL,SAME,NORM,NORMF,NORMP
      integer :: IS1,JS1,KS1,LS1,I1,I2,J1,J2,NX,NY,NZ,K1,K2,L11,L2
      !IS1=kstart(ii)
      !rys
      INTEGER,PARAMETER :: NUM_QUART=9169597 !number of quartets
      INTEGER :: IDX_RYS_2(4,NUM_QUART,NCP_rys)
      INTEGER :: N_RYS_2(NCP_rys)

      INTEGER :: NSA_rys,ISA_rys !sum angular momentum
      INTEGER :: NCP_rys,ICP_rys !product contraction

      INTEGER :: MAXRYS
      iNTEGER :: MM_rys
      integer :: NROOTS
      REAL(KIND=fp) :: HFSCAL,CSCALT
      double precision,parameter :: HALF=0.5D00
      double precision,parameter :: FOUR=4.0D00
      !stuff for digestion
      !INTEGER :: II1,JJ1,KK1,LL1, I2,J2,K2,L2
      INTEGER :: II2,JJ2,KK2,IL,JK,JL
      DOUBLE PRECISION :: XVAL4, VAL,VAL1,VAL4
      REAL(KIND=FP) :: CUTOFFAO

      INTEGER::IJN,KLN,NN,I,J,K,L
      double precision :: N1
      ! !new
      ! double precision :: N,NM

      DOUBLE PRECISION :: XIN(31213),YIN(31213),ZIN(31213)
      !added for dirfck
      integer :: i_index,ij_index,ijk_index,itmp,nkl,ijkl_index
      !stuff for IJPRIM
      integer :: n,nm
      double precision :: ai,arri,axi,xi,ayi,yi,azi,zi
      double precision :: aj,aainv,axj,xj,ayj,yj,azj,zj
      double precision :: A,R,X1,Y1,Z1
      double precision :: aa,tol
      double precision :: csi,cpi,cdi
      double precision :: csj,cpj,cdj
      double precision :: csk,cpk,cdk
      double precision :: csl,cpl,cdl
      double precision, PARAMETER :: SQRT3=1.73205080756888D+00

      NROOTS=2

      !write(*,*) "N_RYS_2", N_RYS_2

       DO ICP_rys=1,NCP_rys


!$omp target teams distribute parallel do default(none) &
!$omp shared(fa,da,ia) &
!$omp shared(IDX_RYS_2,N_RYS_2) &
!$omp shared(ngth,c) &
!$omp shared(dfttyp) &
!$omp shared(KMAX,KMIN,KTYPE,KLOC,KNG,KSTART,KATOM) &
!$omp shared(EX,CS,CP,CD) &
!$omp shared(NROOTS,ICP_rys) &
!$omp shared(qq4) &
!$omp shared(ix,iy,iz,jx,jy,jz,kx,ky,kz,lx,ly,lz) &
!!$omp private(XIN,YIN,ZIN) &
!$omp private(XIN,YIN,ZIN) &
!$omp private(bx,by,bz,i1,j1,i2,j2)  &
!$omp private(az,rab1,knu,ay,kandl)  &
!$omp private(same,iandj,cpd1,cdd) &
!$omp private(ngtj,ax,ngti,inu,jnu)  &
!$omp private(kl,nx,rcd1,l11,l2,nz)  &
!$omp private(ny,dz,cx,cy)  &
!$omp private(ngtl,lnu,ngtk,dx,dy,k2)  &
!$omp private(cz,k1,ks1,ls1,js1,maxl)  &
!$omp private(is1,llt,nga1,lkt,lit,ljt) &
!$omp private(minl,ll,mini,kk,ii,jj)   &

!ijprim
!zqout
!$omp private(ijn,kln,n1)&
!$omp private(ddij)&
!$omp private(nm,n,nn)&
!$omp private(csi,cpi,cdi)&
!$omp private(csj,cpj,cdj)&
!$omp private(csk,cpk,cdk)&
!$omp private(csl,cpl,cdl)&
!$omp private(ai,arri,axi,ayi,azi)&
!$omp private(aj,aa,aainv,dum1,dum2)&

!$omp private(mink,maxk,maxj,maxi,minj)  &
!$omp private(cdb,gc,cpb1,gb,csb1,gd)  &
!$omp private(csd1,cdc,csc1,cpc1,cda)  &
!$omp private(loci,locj,ngd1,ngb1,ngc1)  &
!$omp private(csa1,cpa1,ga,lock,locl) &
!$omp private(ijx,ijy,ijz,ijgt,ij) &
!$omp private(klx,kly,klz,klgt,ijkl,ik,max,nij) &
!$omp private(ijd,z1,y1,r,a,x1,ghondo,dij,dkl) &
!!!!added stuff for dirfck
!$omp private(jk,il,jl,kk2,val4,val1,val,xval4,xval1,jj2,ii2,l1,ikk) &
!$omp private(nkl,itmp,ijkl_index,i_index,ij_index,ijk_index) &
!$omp private(nint,cscalt,hfscal)

! ! !FFF
       DO MM_rys=1,N_RYS_2(ICP_rys)
      !write(*,*) "i am in rys2"


      ii=IDX_RYS_2(1,MM_rys,ICP_rys)
      jj=IDX_RYS_2(2,MM_rys,ICP_rys)
      kk=IDX_RYS_2(3,MM_rys,ICP_rys)
      ll=IDX_RYS_2(4,MM_rys,ICP_rys)

      ! MINI = kmin(ii)
      ! MAXI = kmax(ii)
      ! MINJ = kmin(jj)
      ! MAXJ = kmax(jj)
      ! MINK = kmin(kk)
      ! MAXK = kmax(kk)
      ! MINL = kmin(ll)
      ! MAXL = kmax(ll)

      IS1= kstart(ii)
      JS1= kstart(jj)
      KS1= kstart(kk)
      LS1= kstart(ll)

      ! lit= ktype(ii)-1
      ! !write(*,*) "lit is", lit
      ! ljt= ktype(jj)-1
      ! !write(*,*) "ljt is", ljt
      ! lkt= ktype(kk)-1
      ! !write(*,*) "lkt is", lkt
      ! llt= ktype(ll)-1     
      ! !write(*,*) "llt is", llt

      ! IF (lit.EQ.2 .AND. ljt.EQ.0 .AND. lkt.EQ.0 .AND. llt.EQ.0) then
      ! write(*,*) "rys 2 does 2000"
      ! ENDIF

      ! IF (lit.EQ.2 .AND. ljt.EQ.1 .AND. lkt.EQ.0 .AND. llt.EQ.0) then
      ! write(*,*) "rys 2 does 2100"
      ! ENDIF


      ! LOCI = kloc(ii)-MINI
      ! LOCJ = kloc(jj)-MINJ
      ! LOCK = kloc(kk)-MINK
      ! LOCL = kloc(ll)-MINL

      ! ga=ex(is1)
      ! csa1=cs(is1)
      ! cpa1=cp(is1)
      ! cda=cd(is1)

      ! gb=ex(js1)
      ! csb1=cs(js1)
      ! cpb1=cp(js1)
      ! cdb=cd(js1)

      ! gc=ex(ks1)
      ! csc1=cs(ks1)
      ! cpc1=cp(ks1)
      ! cdc=cd(ks1)

      ! gd=ex(ls1)
      ! csd1=cs(ls1)
      ! cpd1=cp(ls1)
      ! cdd=cd(ls1)

      IANDJ = II .EQ. JJ
      KANDL = KK .EQ. LL
      SAME = (II .EQ. KK) .AND. (JJ .EQ. LL)

!=================first subroutine============
      IF (KTYPE(ii) .LT. KTYPE(jj)) THEN
         INU = jj
         JNU = ii
         NGTI = NGTH(2)
         NGTJ = NGTH(1)
      ELSE
         INU = ii
         JNU = jj
         NGTI = NGTH(1)
         NGTJ = NGTH(2)
      END IF
      IF (KTYPE(kk) .LT. KTYPE(ll)) THEN
         KNU = ll
         LNU = kk
         NGTK = NGTH(4)
         NGTL = NGTH(3)
      ELSE
         KNU = kk
         LNU = ll
         NGTK = NGTH(3)
         NGTL = NGTH(4)
      END IF

      MINI = kmin(inu)
      MAXI = kmax(inu)
      MINJ = kmin(jnu)
      MAXJ = kmax(jnu)
      MINK = kmin(knu)
      MAXK = kmax(knu)
      MINL = kmin(lnu)
      MAXL = kmax(lnu)


      LOCI = kloc(inu)-MINI
      LOCJ = kloc(jnu)-MINJ
      LOCK = kloc(knu)-MINK
      LOCL = kloc(lnu)-MINL

      I = KATOM(INU)
      J = KATOM(JNU)
      K = KATOM(KNU)
      L = KATOM(LNU)
      !write(*,*) "I is ", I
      AX = C(1,I)
      AY = C(2,I)
      AZ = C(3,I)
      BX = C(1,J)
      BY = C(2,J)
      BZ = C(3,J)
      CX = C(1,K)
      CY = C(2,K)
      CZ = C(3,K)
      DX = C(1,L)
      DY = C(2,L)
      DZ = C(3,L)
      !I1 = KSTART(INU)
      !I2 = I1+KNG(INU)-1
      LIT = KTYPE(INU)
      LJT = KTYPE(JNU)
      LKT = KTYPE(KNU)
      LLT = KTYPE(LNU)

      RAB1 = ((AX-BX)*(AX-BX) + (AY-BY)*(AY-BY) + (AZ-BZ)*(AZ-BZ))
      RCD1 = ((CX-DX)*(CX-DX) + (CY-DY)*(CY-DY) + (CZ-DZ)*(CZ-DZ))
!this is the same
      IJ = 0
      DO I = MINI,MAXI
         NX = IX(I)
         NY = IY(I)
         NZ = IZ(I)
         IF (IANDJ) MAXJ = I
         DO J = MINJ,MAXJ
            IJ = IJ+1
            IJX(IJ) = NX+JX(J)
            !write(*,*) "IJX(IJ)", IJX(IJ)
            IJY(IJ) = NY+JY(J)
            IJZ(IJ) = NZ+JZ(J)
            IJGT(IJ) = NGTI*(I-MINI)+NGTJ*(J-MINJ)+1
         ENDDO
      ENDDO
    
      KL = 0
      DO K = MINK,MAXK
         NX = KX(K)
         NY = KY(K)
         NZ = KZ(K)
         IF (KANDL) MAXL = K
         DO L = MINL,MAXL
            KL = KL+1
            KLX(KL) = NX+LX(L)
            KLY(KL) = NY+LY(L)
            KLZ(KL) = NZ+LZ(L)
            KLGT(KL) = NGTK*(K-MINK)+NGTL*(L-MINL)
         ENDDO
      ENDDO
      MAX = KL
      DO I = 1,IJ
      IF (SAME) MAX = I
      IK(I) = MAX
      ENDDO
      IJKL = IJ*KL
      IF (SAME) IJKL = IJ*(IJ+1)/2

!====================get density for IJ =================

         ! CALL shells_gpu2_ij(1,ii,jj,.TRUE.,&
         !      IJGT,IJX,IJY,IJZ,&
         !      IANDJ,SAME,&
         !      ii,jj,NGTI,NGTJ,&
         !      LIT,LJT,LOCI,LOCJ,&
         !      mini,maxi,minj,maxj,&
         !      NIJ,IJ,KL,IJKL,&
         !      NGTH,C,&
         !      kstart,katom,ktype,kng,kloc,kmin,kmax,&
         !      ex,cs,cp,cd,&
         !      ga,csa1,cpa1,cda,&
         !      gb,csb1,cpb1,cdb,&
         !      AX,AY,AZ,BX,BY,BZ,RAB1,&
         !      NGA1,NGB1)

         ! CALL shells_gpu2_kl(2,ii,jj,kk,ll,.TRUE.,&
         !      IK,KLGT,KLX,KLY,KLZ,&
         !      IANDJ,KANDL,SAME,&
         !      kk,ll,NGTK,NGTL,&
         !      LKT,LLT,LOCK,LOCL,&
         !      mink,maxk,minl,maxl,&
         !      IJ,KL,IJKL,&
         !      NGTH,C,&
         !      kstart,katom,ktype,kng,kloc,kmin,kmax,&
         !      ex,cs,cp,cd,&
         !      gc,csc1,cpc1,cdc,&
         !      gd,csd1,cpd1,cdd,&
         !      CX,CY,CZ,DX,DY,DZ,RCD1,&
         !      NGC1,NGD1)
!       ii=IDX_RYS_2(1,MM_rys,ICP_rys)
!       jj=IDX_RYS_2(2,MM_rys,ICP_rys)
!       kk=IDX_RYS_2(3,MM_rys,ICP_rys)
!       ll=IDX_RYS_2(4,MM_rys,ICP_rys)
!       !write(*,*) "in loop"
!       MINI = kmin(ii)
!       MAXI = kmax(ii)
!       MINJ = kmin(jj)
!       MAXJ = kmax(jj)
!       MINK = kmin(kk)
!       MAXK = kmax(kk)
!       MINL = kmin(ll)
!       MAXL = kmax(ll)

!       IS1= kstart(ii)
!       JS1= kstart(jj)
!       KS1= kstart(kk)
!       LS1= kstart(ll)

!       NGA1= kng(ii)
!       NGB1= kng(jj)
!       NGC1= kng(kk)
!       NGD1= kng(ll)

!       LOCI = kloc(ii)-MINI
!       LOCJ = kloc(jj)-MINJ
!       LOCK = kloc(kk)-MINK
!       LOCL = kloc(ll)-MINL

!       ga=ex(is1)
!       csa1=cs(is1)
!       cpa1=cp(is1)
!       cda=cd(is1)

!       gb=ex(js1)
!       csb1=cs(js1)
!       cpb1=cp(js1)
!       cdb=cd(js1)

!       gc=ex(ks1)
!       csc1=cs(ks1)
!       cpc1=cp(ks1)
!       cdc=cd(ks1)

!       gd=ex(ls1)
!       csd1=cs(ls1)
!       cpd1=cp(ls1)
!       cdd=cd(ls1)

!       IANDJ = II .EQ. JJ
!       KANDL = KK .EQ. LL
!       SAME = (II .EQ. KK) .AND. (JJ .EQ. LL)

! !   fill indices for Rys quadrature
!       NGTI = NGTH(1)
!       NGTJ = NGTH(2)
!       NGTK = NGTH(3)
!       NGTL = NGTH(4)
!       !write(*,*) "here"

!       I = KATOM(ii)
!       J = KATOM(jj)
!       K = KATOM(kk)
!       L = KATOM(ll)
!       !write(*,*) "here"

      ! AX = C(1,I)
      ! AY = C(2,I)
      ! AZ = C(3,I)
      ! BX = C(1,J)
      ! BY = C(2,J)
      ! BZ = C(3,J)
      ! CX = C(1,K)
      ! CY = C(2,K)
      ! CZ = C(3,K)
      ! DX = C(1,L)
      ! DY = C(2,L)
      ! DZ = C(3,L)

      ! LIT = KTYPE(ii)
      ! LJT = KTYPE(jj)
      ! LKT = KTYPE(kk)
      ! LLT = KTYPE(ll)

      ! RAB1 = ((AX-BX)*(AX-BX) + (AY-BY)*(AY-BY) + (AZ-BZ)*(AZ-BZ))
      ! RCD1 = ((CX-DX)*(CX-DX) + (CY-DY)*(CY-DY) + (CZ-DZ)*(CZ-DZ))
      ! !write(*,*) "here"
      
      ! IJ = 0
      ! DO I = MINI,MAXI
      !    NX = IX(I)
      !    NY = IY(I)
      !    NZ = IZ(I)
      !    IF (IANDJ) MAXJ = I
      !    DO J = MINJ,MAXJ
      !       IJ = IJ+1
      !       !write(*,*) "here"
      !       IJX(IJ) = NX+JX(J)
      !       IJY(IJ) = NY+JY(J)
      !       IJZ(IJ) = NZ+JZ(J)
      !       IJGT(IJ) = NGTI*(I-MINI)+NGTJ*(J-MINJ)+1
      !    ENDDO
      ! ENDDO
   
      ! KL = 0
      ! DO K = MINK,MAXK
      !    NX = KX(K)
      !    NY = KY(K)
      !    NZ = KZ(K)
      !    IF (KANDL) MAXL = K
      !    DO L = MINL,MAXL
      !       KL = KL+1
      !       KLX(KL) = NX+LX(L)
      !       KLY(KL) = NY+LY(L)
      !       KLZ(KL) = NZ+LZ(L)
      !       KLGT(KL) = NGTK*(K-MINK)+NGTL*(L-MINL)
      !    ENDDO
      ! ENDDO
      ! !time is 0.07 here

      ! MAX = KL
      ! DO I = 1,IJ
      ! IF (SAME) MAX = I
      ! IK(I) = MAX
      ! ENDDO
      ! IJKL = IJ*KL
      ! IF (SAME) IJKL = IJ*(IJ+1)/2
      ! !time is 0.11 here
!====================get density for IJ =================
      MAX = MAXJ
      N = 0
      NN = 0
      NM = -2**20
      DO 180 I = MINI,MAXI
         GO TO (100,100,120,120,100,120,120,100,120,120,&
               100,120,120,100,120,120,120,120,120,100,&
               100,120,120,100,120,120,120,120,120,100,&
               120,120,100,120,120,&
               100,120,120,100,120,120,120,120,120,100,&
               120,120,120,120,120,100,120,120,100,120,&
               120,&
               100,120,120,100,120,120,120,120,120,100,&
               120,120,120,120,120,100,120,120,100,120,&
               120,100,120,120,120,120,120,100),I
  100    NM = NN
  120    NN = NM
         IF (IANDJ) MAX = I
         DO 170 J = MINJ,MAX
            GO TO (140,140,160,160,140,160,160,140,160,160,&
                  140,160,160,140,160,160,160,160,160,140,&
                  140,160,160,140,160,160,160,160,160,140,&
                  160,160,140,160,160,&
                  140,160,160,140,160,160,160,160,160,140,&
                  160,160,160,160,160,140,160,160,140,160,&
                  160,&
                  140,160,160,140,160,160,160,160,160,140,&
                  160,160,160,160,160,140,160,160,140,160,&
                  160,140,160,160,160,160,160,140),J
  140       NN = NN+1
  160       N = N+1
            IJD(N) = NN
  170    CONTINUE
  180 CONTINUE

          AI = ex(is1)
          !write(*,*) "AI"
          ARRI = AI*RAB1
          AXI = AI*AX
          AYI = AI*AY
          AZI = AI*AZ
          CSI = cs(is1)
          CPI = cp(is1)
          CDI = cd(is1)
          AJ = ex(js1)
          AA = AI+AJ
          AAINV = ONE/AA
          !DUM = AJ*ARRI*AAINV
          !IF (DUM .GT. TOL) GO TO 520
          CSJ = cs(js1)
          CPJ = cp(js1)
          CDJ = cd(js1)

          NN = 0
          R = AJ*ARRI*AAINV
          A = AA
          X1 = (AXI+AJ*BX)*AAINV
          Y1 = (AYI+AJ*BY)*AAINV
          Z1 = (AZI+AJ*BZ)*AAINV


            DUM1 = ZERO
            DUM2 = ZERO
            DO 420 I = MINI,MAXI
               GO TO (200,220,420,420,240,420,420,260,420,420),I
  200          DUM1 = CSI*AAINV
               GO TO 280
  220          DUM1 = CPI*AAINV
               GO TO 280
  240          DUM1 = CDI*AAINV
               GO TO 280
  260          DUM1 = DUM1*SQRT3
               GO TO 280

  280          IF (IANDJ) MAX = I
                 DO 400 J = MINJ,MAX
                  GO TO (300,320,400,400,340,400,400,360,400,400),J
  300             DUM2 = DUM1*CSJ
                  GO TO 380
  320             DUM2 = DUM1*CPJ
                  GO TO 380
  340             DUM2 = DUM1*CDJ
                  GO TO 380
  360             DUM2 = DUM2*SQRT3
                  GO TO 380

  380             NN = NN+1
                  DDIJ(NN) = DUM2
  400          CONTINUE
420          CONTINUE
         ! CALL ijprim_gpu_2(ddij,a,r,x1,y1,z1,ijd,&
         !      IANDJ,KANDL,SAME,&
         !      LIT,mini,maxi,minj,maxj,NIJ,&
         !      ga,csa1,cpa1,cda,&
         !      gb,csb1,cpb1,cdb,&
         !      AX,AY,AZ,BX,BY,BZ,RAB1,&
         !      NGA1,NGB1)

!====================zero out ghondo===========
      IJN = 0
      DO I = MINI,MAXI
         IF (IANDJ) MAXJ = I
         DO J = MINJ,MAXJ
            IJN = IJN+1
            N1 = IJGT(IJN)
            KLN = 0
            DO K =  MINK,MAXK
               IF (KANDL) MAXL = K
               DO L = MINL,MAXL
                  KLN = KLN+1
                  NN = N1+KLGT(KLN)
                  GHONDO(NN) = 0
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      ! !0.27 here
      !    call zqout_gpu_2(ghondo,IJGT,KLGT,&
      !         IANDJ,KANDL,SAME,&
      !         mini,maxi,minj,maxj,mink,maxk,minl,maxl)

         CALL GENRAL_GPU_2(ghondo,ddij,dkl,dij,&
             IJGT,IJX,IJY,IJZ,IK,KLGT,KLX,KLY,KLZ, &
             A,R,X1,Y1,Z1,IJD,iandj,kandl,same,&
             LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,&
             qq4,mini,maxi,minj,maxj,mink,maxk,minl,maxl,&
             NIJ,IJ,KL,IJKL,&
             ga,csa1,cpa1,cda,&
             gb,csb1,cpb1,cdb,&
             gc,csc1,cpc1,cdc,&
             gd,csd1,cpd1,cdd,&
             AX,AY,AZ,BX,BY,BZ,RAB1,CX,CY,CZ,DX,DY,DZ,RCD1,&
             NGA1,NGB1,NGC1,NGD1,&
             XIN,YIN,ZIN)

       HFSCAL=DFTTYP(3)
       CSCALT=1.0D+00

      XVAL1 = HFSCAL
      XVAL4 = FOUR*CSCALT

      NIJ = 0
      I_INDEX = 1
      DO I = MINI,MAXI
         IF (IANDJ) MAXJ = I

         I1 = I+LOCI

         IJ_INDEX = I_INDEX
         !write(*,*) "IJ_INDEX", IJ_INDEX
         I_INDEX = I_INDEX + NGTI

         DO J = MINJ,MAXJ
            NIJ = NIJ+1

            J1 = J+LOCJ
            I2 = I1
            J2 = J1
            IF (I1.LT.J1) THEN ! SORT <IJ|
               I2 = J1
               J2 = I1
            ENDIF

            IJK_INDEX = IJ_INDEX
            IJ_INDEX = IJ_INDEX + NGTJ

            NKL = NIJ

            DO K =  MINK,MAXK
               IF (KANDL) MAXL = K

               K1 = K + LOCK

               IF(SAME) THEN ! ACCOUNT FOR NON-UNIQUE PERMUTATIONS
                  ITMP = MIN(MAXL-MINL+1,NKL)
                  IF (ITMP.EQ.0) CYCLE
                  MAXL = MINL + ITMP - 1
                  NKL = NKL - ITMP
               ENDIF

               IJKL_INDEX = IJK_INDEX
               IJK_INDEX = IJK_INDEX + NGTK

               DO L=MINL,MAXL

                  VAL = GHONDO( IJKL_INDEX )
                  !IF(VAL .NE. 0.00) write(*,*) "val is", VAL
                  !write(*,*) "val is", VAL
                  IJKL_INDEX = IJKL_INDEX + NGTL
                  !IF(ABS(VAL).LT.CUTOFF) CYCLE ! GOTO 300
                  NINT = NINT + 1

                  L1 = L + LOCL
                  K2 = K1
                  L2 = L1

                  IF (K2.LT.L2) THEN ! SORT |KL>
                     K2 = L1
                     L2 = K1
                  ENDIF

                  II = I2
                  JJ = J2
                  KK = K2
                  LL = L2

                  IF (II.LT.KK) THEN ! SORT <IJ|KL>
                     II = K2
                     JJ = L2
                     KK = I2
                     LL = J2
                  ELSE IF (II.EQ.KK.AND.JJ.LT.LL) THEN ! SORT <IJ|IL>
                     JJ = L2
                     LL = J2
                  ENDIF

                  II2 = IA(II)
                  JJ2 = IA(JJ)
                  KK2 = IA(KK)

                  IJ = II2 + JJ
                  IKK = II2 + KK
                  IL = II2 + LL
                  JK = JJ2 + KK
                  JL = JJ2 + LL
                  KL = KK2 + LL
                  IF (JJ.LT.KK) JK = KK2 + JJ
                  IF (JJ.LT.LL) JL = IA(LL) + JJ
!
!       ACCOUNT FOR IDENTICAL PERMUTATIONS.
!
                  IF(II.EQ.JJ) VAL = VAL*HALF
                  IF(KK.EQ.LL) VAL = VAL*HALF
                  IF(II.EQ.KK.AND.JJ.EQ.LL) VAL = VAL*HALF
                  VAL1 = VAL*XVAL1
                  VAL4 = VAL*XVAL4
!$omp atomic
                  FA(IJ) = FA(IJ) + VAL4*DA(KL)
!$omp atomic
                  FA(KL) = FA(KL) + VAL4*DA(IJ)
!$omp atomic
                  FA(IKK) = FA(IKK) - VAL1*DA(JL)
!$omp atomic
                  FA(JL) = FA(JL) - VAL1*DA(IKK)
!$omp atomic
                  FA(IL) = FA(IL) - VAL1*DA(JK)
!$omp atomic
                  FA(JK) = FA(JK) - VAL1*DA(IL)
               ENDDO
            ENDDO
         ENDDO
      ENDDO


      ENDDO
!$omp end target teams distribute parallel do
  ENDDO


  END SUBROUTINE rysint_gpu_2

  SUBROUTINE rysint_gpu_2000(ii,jj,kk,ll,ktype,ghondo,ddij,ngth,c,&
                        kstart,katom,kng,kloc,kmin,kmax,&
                        ex,cs,cp,cd,&
                        NGTI,NGTJ,NGTK,NGTL,&
                        NIJ,IJ,KL,IJKL,&
                        qq4,&
                        AX,AY,AZ,BX,BY,BZ,RAB,CX,CY,CZ,DX,DY,DZ,RCD1,&
                        dfttyp,&
                        ia,da,fa,nint,&
                        l1,l2a,nflmat,&
                        IDX_RYS_2000,N_RYS_2000,NCP_rys)

    USE omp_lib
    USE prec, ONLY: fp
    USE mx_limits, ONLY: &
    mxsh,mxg2,mxatm,mxgtot,mxgsh

    IMPLICIT DOUBLE PRECISION (A-H,O-Z)

    !INTEGER, INTENT(IN) :: &
    !l1, l2a,nflmat, ia(l1)

    !for dirfck
    REAL(KIND=fp),INTENT(OUT) :: fa(l2a*nflmat) !FOCK MATRIX
    REAL(KIND=fp),INTENT(INOUT) :: da(l2a) !DENSITY MATRIX
    INTEGER,INTENT(IN) :: ia(l1) !TRIANGULAR MATRIX 
    INTEGER :: l2a,nflmat,l1,nint 

    REAL(KIND=fp),DIMENSION(20) :: DFTTYP

    INTEGER, INTENT(IN) :: ktype(mxsh)
    INTEGER :: ii,jj,kk,ll,inu,jnu,knu,lnu

    REAL(KIND=fp),INTENT(INOUT),ALLOCATABLE :: &
      ghondo(:)

    REAL(KIND=fp),INTENT(INOUT) :: &
      ddij(:)

    REAL(KIND=fp),DIMENSION(784) :: DKL,DIJ
    REAL(KIND=fp),DIMENSION(784) :: &
     IJGT,IJX,IJY,IJZ,IK,KLGT,KLX,KLY,KLZ
    REAL(KIND=fp),DIMENSION(84) :: &
     IX,IY,IZ,JX,JY,JZ,KX,KY,KZ,LX,LY,LZ
     INTEGER:: IJ_count

    REAL(KIND=fp),DIMENSION(MXG2) :: A,R,X1,Y1,Z1
    REAL(KIND=fp),DIMENSION(784) :: IJD

      !for shlexc
      INTEGER,DIMENSION(4) :: NGTH
      !for infoa
      REAL(KIND=fp) :: C(3,mxatm)
      !for roots
      REAL(KIND=fp),DIMENSION(13) :: U,W
      REAL(KIND=fp) :: XX
      !for eriout
      INTEGER :: NGTI,NGTJ,NGTK,NGTL
      !shlnos
      INTEGER ::mini,maxi,minj,maxj,mink,maxk,minl,maxl,&
                LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,&
                NIJ,IJ,KL,IJKL, max,IKK
      REAL(KIND=fp) :: qq4

      !nshell
      INTEGER,DIMENSION(mxsh) :: kstart,katom,kng,kloc,kmin,kmax
      REAL(KIND=fp),DIMENSION(mxgtot) :: ex,cs,cp,cd,cf,cg,ch,ci

      !shlinf
      double precision :: ga(mxgsh),csa1(mxgsh),cpa1(mxgsh),cda(mxgsh)
      double precision :: gb(mxgsh),csb1(mxgsh),cpb1(mxgsh),cdb(mxgsh)
      double precision :: gc(mxgsh),csc1(mxgsh),cpc1(mxgsh),cdc(mxgsh)
      double precision :: gd(mxgsh),csd1(mxgsh),cpd1(mxgsh),cdd(mxgsh)

      double precision :: AX,AY,AZ,BX,BY,BZ,RAB,CX,CY,CZ,DX,DY,DZ,RCD1
      integer :: NGA1,NGB1,NGC1,NGD1
      LOGICAL IANDJ,KANDL,SAME,NORM,NORMF,NORMP
      integer :: IS1,JS1,KS1,LS1,I1,I2,J1,J2,NX,NY,NZ,K1,K2,L11,L2
      !IS1=kstart(ii)
      !rys
      INTEGER,PARAMETER :: NUM_QUART=9169597 !number of quartets
      INTEGER :: IDX_RYS_2000(4,NUM_QUART,NCP_rys)
      INTEGER :: N_RYS_2000(NCP_rys)

      !INTEGER :: NSA_rys,ISA_rys !sum angular momentum
      INTEGER :: NCP_rys,ICP_rys !product contraction

      INTEGER :: MAXRYS
      iNTEGER :: MM_rys
      integer :: NROOTS
      REAL(KIND=fp) :: HFSCAL,CSCALT
      double precision,parameter :: HALF=0.5D00
      double precision,parameter :: FOUR=4.0D00
      !stuff for digestion
      !INTEGER :: II1,JJ1,KK1,LL1, I2,J2,K2,L2
      INTEGER :: II2,JJ2,KK2,IL,JK,JL
      DOUBLE PRECISION :: XVAL4, VAL,VAL1,VAL4
      REAL(KIND=FP) :: CUTOFFAO

      INTEGER::IJN,KLN,NN,I,J,K,L
      double precision :: N1
      !new
      double precision :: N,NM

      DOUBLE PRECISION :: XIN(31213),YIN(31213),ZIN(31213)
      !added for dirfck
      integer :: i_index,ij_index,ijk_index,itmp,nkl,ijkl_index

      !NSA_rys=9
      !ISA_rys=1
      !write(*,*) "N_RYS_2000", N_RYS_2000
      NROOTS=2

      !write(*,*) "in N_RYS_2000"

      DO ICP_rys=1,NCP_rys

!$omp target teams distribute parallel do default(none) &
!$omp shared(fa,da,ia) &
!$omp shared(IDX_RYS_2000,N_RYS_2000) &
!$omp shared(ngth,c) &
!$omp shared(dfttyp) &
!$omp shared(KMAX,KMIN,KTYPE,KLOC,KNG,KSTART,KATOM) &
!$omp shared(EX,CS,CP,CD,CF,CG,CH,CI) &
!$omp shared(NROOTS,ICP_rys) &
!$omp shared(qq4) &
!!$omp private(XIN,YIN,ZIN) &
!$omp private(XIN,YIN,ZIN) &
!$omp private(bx,by,bz,i1,j1,i2,j2)  &
!$omp private(az,rab1,knu,ay,kandl)  &
!$omp private(same,iandj,cpd1,cdd) &
!$omp private(ngtj,ax,ngti,inu,jnu)  &
!$omp private(kl,nx,rcd1,l11,l2,nz)  &
!$omp private(kz,ky,kx,ny,dz,cx,cy)  &
!$omp private(ngtl,lnu,ngtk,dx,dy,k2)  &
!$omp private(cz,k1,ks1,ls1,js1,maxl)  &
!$omp private(is1,llt,nga1,lkt,lit,ljt) &
!$omp private(minl,ll,mini,kk,ii,jj)   &
!$omp private(mink,maxk,maxj,maxi,minj)  &
!$omp private(cdb,gc,cpb1,gb,csb1,gd)  &
!$omp private(csd1,cdc,csc1,cpc1,cda)  &
!$omp private(loci,locj,ngd1,ngb1,ngc1)  &
!$omp private(csa1,cpa1,ga,lock,locl) &
!$omp private(jy,ijx,ijy,ijz,ijgt,jz,jx,ix,ij,iy,iz) &
!$omp private(klx,kly,klz,klgt,lz,lx,ly,ijkl,ik,max,nij) &
!$omp private(ijd,z1,y1,a,ddij,x1,r,ghondo,dij,dkl) &
!!!!added stuff for dirfck
!$omp private(MM_rys) &
!$omp private(jk,il,jl,kk2,val4,val1,val,xval4,xval1,jj2,ii2,l1,ikk) &
!$omp private(nkl,itmp,ijkl_index,i_index,ij_index,ijk_index) &
!$omp private(nint,cscalt,hfscal)


      DO MM_rys=1,N_RYS_2000(ICP_rys)
      !write(*,*) "i am in rys2"
      ! write(*,*) "N_RYS_2000", N_RYS_2000



      ii=IDX_RYS_2000(1,MM_rys,ICP_rys)
      jj=IDX_RYS_2000(2,MM_rys,ICP_rys)
      kk=IDX_RYS_2000(3,MM_rys,ICP_rys)
      ll=IDX_RYS_2000(4,MM_rys,ICP_rys)

      MINI = kmin(ii)
      MAXI = kmax(ii)
      MINJ = kmin(jj)
      MAXJ = kmax(jj)
      MINK = kmin(kk)
      MAXK = kmax(kk)
      MINL = kmin(ll)
      MAXL = kmax(ll)

      IS1= kstart(ii)
      JS1= kstart(jj)
      KS1= kstart(kk)
      LS1= kstart(ll)

      lit= ktype(ii)-1
      write(*,*) "lit is", lit
      ljt= ktype(jj)-1
      !write(*,*) "ljt is", ljt
      lkt= ktype(kk)-1
      !write(*,*) "lkt is", lkt
      llt= ktype(ll)-1     
      !write(*,*) "llt is", llt

      ! IF (lit.EQ.2 .AND. ljt.EQ.0 .AND. lkt.EQ.0 .AND. llt.EQ.0) then
      ! write(*,*) "rys 2 does 2000"
      ! ENDIF

      ! IF (lit.EQ.2 .AND. ljt.EQ.1 .AND. lkt.EQ.0 .AND. llt.EQ.0) then
      ! write(*,*) "rys 2 does 2100"
      ! ENDIF

      NGA1= kng(ii)
      NGB1= kng(jj)
      NGC1= kng(kk)
      NGD1= kng(ll)

      LOCI = kloc(ii)-MINI
      LOCJ = kloc(jj)-MINJ
      LOCK = kloc(kk)-MINK
      LOCL = kloc(ll)-MINL

      ga=ex(is1)
      csa1=cs(is1)
      cpa1=cp(is1)
      cda=cd(is1)

      gb=ex(js1)
      csb1=cs(js1)
      cpb1=cp(js1)
      cdb=cd(js1)

      gc=ex(ks1)
      csc1=cs(ks1)
      cpc1=cp(ks1)
      cdc=cd(ks1)

      gd=ex(ls1)
      csd1=cs(ls1)
      cpd1=cp(ls1)
      cdd=cd(ls1)

      IANDJ = II .EQ. JJ
      KANDL = KK .EQ. LL
      SAME = (II .EQ. KK) .AND. (JJ .EQ. LL)


      ! IF (KTYPE(ii) .LT. KTYPE(jj)) THEN
      !    INU = jj
      !    JNU = ii
      !    NGTI = NGTH(2)
      !    NGTJ = NGTH(1)
      ! ELSE
      !    INU = ii
      !    JNU = jj
      !    NGTI = NGTH(1)
      !    NGTJ = NGTH(2)
      ! END IF

      ! I = KATOM(INU)
      ! !write(*,*) "I is ", I
      ! AX = C(1,I)
      ! AY = C(2,I)
      ! AZ = C(3,I)
      ! I1 = KSTART(INU)
      ! I2 = I1+KNG(INU)-1
      ! LIT = KTYPE(INU)
      ! MINI = KMIN(INU)
      ! MAXI = KMAX(INU)
      ! LOCI = KLOC(INU)-MINI
      ! NGA1 = 0
      ! !write(*,*) "I1 is ", I1
      ! !write(*,*) "I2 is ", I2
      ! DO I = I1,I2
      !    NGA1 = NGA1+1
      !    GA(NGA1) = EX(I)
      !    CSA1(NGA1) = CS(I)
      !    CPA1(NGA1) = CP(I)
      !    CDA(NGA1) = CD(I)
      ! ENDDO

      ! J = KATOM(JNU)
      ! BX = C(1,J)
      ! BY = C(2,J)
      ! BZ = C(3,J)
      ! J1 = KSTART(JNU)
      ! J2 = J1+KNG(JNU)-1
      ! LJT = KTYPE(JNU)
      ! MINJ = KMIN(JNU)
      ! MAXJ = KMAX(JNU)
      ! LOCJ = KLOC(JNU)-MINJ
      ! NGB1 = 0
      !  DO J = J1,J2
      !    NGB1 = NGB1+1
      !    GB(NGB1) = EX(J)
      !    CSB1(NGB1) = CS(J)
      !    CPB1(NGB1) = CP(J)
      !    CDB(NGB1) = CD(J)
      !  ENDDO
      ! RAB1 = ((AX-BX)*(AX-BX) + (AY-BY)*(AY-BY) + (AZ-BZ)*(AZ-BZ))

      ! IJ = 0
      ! DO I = MINI,MAXI
      !    NX = IX(I)
      !    NY = IY(I)
      !    NZ = IZ(I)
      !    IF (IANDJ) MAXJ = I
      !    DO J = MINJ,MAXJ
      !       IJ = IJ+1
      !       IJX(IJ) = NX+JX(J)
      !       !write(*,*) "IJX(IJ)", IJX(IJ)
      !       IJY(IJ) = NY+JY(J)
      !       IJZ(IJ) = NZ+JZ(J)
      !       IJGT(IJ) = NGTI*(I-MINI)+NGTJ*(J-MINJ)+1
      !    ENDDO
      ! ENDDO



         CALL shells_gpu2_ij(1,ii,jj,.TRUE.,&
              IJGT,IJX,IJY,IJZ,&
              IANDJ,SAME,&
              ii,jj,NGTI,NGTJ,&
              LIT,LJT,LOCI,LOCJ,&
              mini,maxi,minj,maxj,&
              NIJ,IJ,KL,IJKL,&
              NGTH,C,&
              kstart,katom,ktype,kng,kloc,kmin,kmax,&
              ex,cs,cp,cd,&
              ga,csa1,cpa1,cda,&
              gb,csb1,cpb1,cdb,&
              AX,AY,AZ,BX,BY,BZ,RAB1,&
              NGA1,NGB1)

         CALL shells_gpu2_kl(2,ii,jj,kk,ll,.TRUE.,&
              IK,KLGT,KLX,KLY,KLZ,&
              IANDJ,KANDL,SAME,&
              kk,ll,NGTK,NGTL,&
              LKT,LLT,LOCK,LOCL,&
              mink,maxk,minl,maxl,&
              IJ,KL,IJKL,&
              NGTH,C,&
              kstart,katom,ktype,kng,kloc,kmin,kmax,&
              ex,cs,cp,cd,&
              gc,csc1,cpc1,cdc,&
              gd,csd1,cpd1,cdd,&
              CX,CY,CZ,DX,DY,DZ,RCD1,&
              NGC1,NGD1)

         CALL ijprim_gpu_2(ddij,a,r,x1,y1,z1,ijd,&
              IANDJ,KANDL,SAME,&
              LIT,mini,maxi,minj,maxj,NIJ,&
              ga,csa1,cpa1,cda,&
              gb,csb1,cpb1,cdb,&
              AX,AY,AZ,BX,BY,BZ,RAB1,&
              NGA1,NGB1)

         call zqout_gpu_2(ghondo,IJGT,KLGT,&
              IANDJ,KANDL,SAME,&
              mini,maxi,minj,maxj,mink,maxk,minl,maxl)

         CALL GENRAL_GPU_2(ghondo,ddij,dkl,dij,&
             IJGT,IJX,IJY,IJZ,IK,KLGT,KLX,KLY,KLZ, &
             A,R,X1,Y1,Z1,IJD,iandj,kandl,same,&
             LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,&
             qq4,mini,maxi,minj,maxj,mink,maxk,minl,maxl,&
             NIJ,IJ,KL,IJKL,&
             ga,csa1,cpa1,cda,&
             gb,csb1,cpb1,cdb,&
             gc,csc1,cpc1,cdc,&
             gd,csd1,cpd1,cdd,&
             AX,AY,AZ,BX,BY,BZ,RAB1,CX,CY,CZ,DX,DY,DZ,RCD1,&
             NGA1,NGB1,NGC1,NGD1,&
             XIN,YIN,ZIN)

         HFSCAL=DFTTYP(3)
         CSCALT=1.0D+00

         ! CALL dirfck_rhf_rys_2(IA,DA,FA,ii,jj,kk,ll, &
         !      GHONDO,NGTI,NGTJ,NGTK,NGTL, &
         !      HFSCAL,CSCALT,NINT,&
         !      kstart,katom,kng,kloc,kmin,kmax)
    
      ! SAME  = ii.EQ.kk.AND.jj.EQ.ll
      ! IANDJ = ii.EQ.jj
      ! KANDL = kk.EQ.ll

      ! MINI = KMIN(ii)
      ! MINJ = KMIN(jj)
      ! MINK = KMIN(kk)
      ! MINL = KMIN(ll)
      ! MAXI = KMAX(ii)
      ! MAXJ = KMAX(jj)
      ! MAXK = KMAX(kk)
      ! MAXL = KMAX(ll)
      ! LOCI = KLOC(ii)-MINI
      ! LOCJ = KLOC(jj)-MINJ
      ! LOCK = KLOC(kk)-MINK
      ! LOCL = KLOC(ll)-MINL

      XVAL1 = HFSCAL
      XVAL4 = FOUR*CSCALT

      NIJ = 0
      I_INDEX = 1
      DO I = MINI,MAXI
         IF (IANDJ) MAXJ = I

         I1 = I+LOCI

         IJ_INDEX = I_INDEX
         !write(*,*) "IJ_INDEX", IJ_INDEX
         I_INDEX = I_INDEX + NGTI

         DO J = MINJ,MAXJ
            NIJ = NIJ+1

            J1 = J+LOCJ
            I2 = I1
            J2 = J1
            IF (I1.LT.J1) THEN ! SORT <IJ|
               I2 = J1
               J2 = I1
            ENDIF

            IJK_INDEX = IJ_INDEX
            IJ_INDEX = IJ_INDEX + NGTJ

            NKL = NIJ

            DO K =  MINK,MAXK
               IF (KANDL) MAXL = K

               K1 = K + LOCK

               IF(SAME) THEN ! ACCOUNT FOR NON-UNIQUE PERMUTATIONS
                  ITMP = MIN(MAXL-MINL+1,NKL)
                  IF (ITMP.EQ.0) CYCLE
                  MAXL = MINL + ITMP - 1
                  NKL = NKL - ITMP
               ENDIF

               IJKL_INDEX = IJK_INDEX
               IJK_INDEX = IJK_INDEX + NGTK

               DO L=MINL,MAXL

                  VAL = GHONDO( IJKL_INDEX )
                  !IF(VAL .NE. 0.00) write(*,*) "val is", VAL
                  !write(*,*) "val is", VAL
                  IJKL_INDEX = IJKL_INDEX + NGTL
                  !IF(ABS(VAL).LT.CUTOFF) CYCLE ! GOTO 300
                  NINT = NINT + 1

                  L1 = L + LOCL
                  K2 = K1
                  L2 = L1

                  IF (K2.LT.L2) THEN ! SORT |KL>
                     K2 = L1
                     L2 = K1
                  ENDIF

                  II = I2
                  JJ = J2
                  KK = K2
                  LL = L2

                  IF (II.LT.KK) THEN ! SORT <IJ|KL>
                     II = K2
                     JJ = L2
                     KK = I2
                     LL = J2
                  ELSE IF (II.EQ.KK.AND.JJ.LT.LL) THEN ! SORT <IJ|IL>
                     JJ = L2
                     LL = J2
                  ENDIF

                  II2 = IA(II)
                  JJ2 = IA(JJ)
                  KK2 = IA(KK)

                  IJ = II2 + JJ
                  IKK = II2 + KK
                  IL = II2 + LL
                  JK = JJ2 + KK
                  JL = JJ2 + LL
                  KL = KK2 + LL
                  IF (JJ.LT.KK) JK = KK2 + JJ
                  IF (JJ.LT.LL) JL = IA(LL) + JJ
!
!       ACCOUNT FOR IDENTICAL PERMUTATIONS.
!
                  IF(II.EQ.JJ) VAL = VAL*HALF
                  IF(KK.EQ.LL) VAL = VAL*HALF
                  IF(II.EQ.KK.AND.JJ.EQ.LL) VAL = VAL*HALF
                  VAL1 = VAL*XVAL1
                  VAL4 = VAL*XVAL4
!$omp atomic
                  FA(IJ) = FA(IJ) + VAL4*DA(KL)
!$omp atomic
                  FA(KL) = FA(KL) + VAL4*DA(IJ)
!$omp atomic
                  FA(IKK) = FA(IKK) - VAL1*DA(JL)
!$omp atomic
                  FA(JL) = FA(JL) - VAL1*DA(IKK)
!$omp atomic
                  FA(IL) = FA(IL) - VAL1*DA(JK)
!$omp atomic
                  FA(JK) = FA(JK) - VAL1*DA(IL)
               ENDDO
            ENDDO
         ENDDO
      ENDDO


      ENDDO
!$omp end target teams distribute parallel do
  ENDDO
  !ENDDO

  END SUBROUTINE rysint_gpu_2000

  SUBROUTINE rysint_gpu_2100(ii,jj,kk,ll,ktype,ghondo,ddij,ngth,c,&
                        kstart,katom,kng,kloc,kmin,kmax,&
                        ex,cs,cp,cd,&
                        NGTI,NGTJ,NGTK,NGTL,&
                        NIJ,IJ,KL,IJKL,&
                        qq4,&
                        AX,AY,AZ,BX,BY,BZ,RAB,CX,CY,CZ,DX,DY,DZ,RCD1,&
                        dfttyp,&
                        ia,da,fa,nint,&
                        l1,l2a,nflmat,&
                        IDX_RYS_2100,N_RYS_2100,NCP_rys)

    USE omp_lib
    USE prec, ONLY: fp
    USE mx_limits, ONLY: &
    mxsh,mxg2,mxatm,mxgtot,mxgsh

    IMPLICIT DOUBLE PRECISION (A-H,O-Z)

    !INTEGER, INTENT(IN) :: &
    !l1, l2a,nflmat, ia(l1)

    !for dirfck
    REAL(KIND=fp),INTENT(OUT) :: fa(l2a*nflmat) !FOCK MATRIX
    REAL(KIND=fp),INTENT(INOUT) :: da(l2a) !DENSITY MATRIX
    INTEGER,INTENT(IN) :: ia(l1) !TRIANGULAR MATRIX 
    INTEGER :: l2a,nflmat,l1,nint 

    REAL(KIND=fp),DIMENSION(20) :: DFTTYP

    INTEGER, INTENT(IN) :: ktype(mxsh)
    INTEGER :: ii,jj,kk,ll,inu,jnu,knu,lnu

    REAL(KIND=fp),INTENT(INOUT),ALLOCATABLE :: &
      ghondo(:)

    REAL(KIND=fp),INTENT(INOUT) :: &
      ddij(:)

    REAL(KIND=fp),DIMENSION(784) :: DKL,DIJ
    REAL(KIND=fp),DIMENSION(784) :: &
     IJGT,IJX,IJY,IJZ,IK,KLGT,KLX,KLY,KLZ
    REAL(KIND=fp),DIMENSION(84) :: &
     IX,IY,IZ,JX,JY,JZ,KX,KY,KZ,LX,LY,LZ
     INTEGER:: IJ_count

    REAL(KIND=fp),DIMENSION(MXG2) :: A,R,X1,Y1,Z1
    REAL(KIND=fp),DIMENSION(784) :: IJD

      !for shlexc
      INTEGER,DIMENSION(4) :: NGTH
      !for infoa
      REAL(KIND=fp) :: C(3,mxatm)
      !for roots
      REAL(KIND=fp),DIMENSION(13) :: U,W
      REAL(KIND=fp) :: XX
      !for eriout
      INTEGER :: NGTI,NGTJ,NGTK,NGTL
      !shlnos
      INTEGER ::mini,maxi,minj,maxj,mink,maxk,minl,maxl,&
                LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,&
                NIJ,IJ,KL,IJKL, max,IKK
      REAL(KIND=fp) :: qq4

      !nshell
      INTEGER,DIMENSION(mxsh) :: kstart,katom,kng,kloc,kmin,kmax
      REAL(KIND=fp),DIMENSION(mxgtot) :: ex,cs,cp,cd,cf,cg,ch,ci

      !shlinf
      double precision :: ga(mxgsh),csa1(mxgsh),cpa1(mxgsh),cda(mxgsh)
      double precision :: gb(mxgsh),csb1(mxgsh),cpb1(mxgsh),cdb(mxgsh)
      double precision :: gc(mxgsh),csc1(mxgsh),cpc1(mxgsh),cdc(mxgsh)
      double precision :: gd(mxgsh),csd1(mxgsh),cpd1(mxgsh),cdd(mxgsh)

      double precision :: AX,AY,AZ,BX,BY,BZ,RAB,CX,CY,CZ,DX,DY,DZ,RCD1
      integer :: NGA1,NGB1,NGC1,NGD1
      LOGICAL IANDJ,KANDL,SAME,NORM,NORMF,NORMP
      integer :: IS1,JS1,KS1,LS1,I1,I2,J1,J2,NX,NY,NZ,K1,K2,L11,L2
      !IS1=kstart(ii)
      !rys
      INTEGER,PARAMETER :: NUM_QUART=9169597 !number of quartets
      INTEGER :: IDX_RYS_2100(4,NUM_QUART,NCP_rys)
      INTEGER :: N_RYS_2100(NCP_rys)

      !INTEGER :: NSA_rys,ISA_rys !sum angular momentum
      INTEGER :: NCP_rys,ICP_rys !product contraction

      INTEGER :: MAXRYS
      iNTEGER :: MM_rys
      integer :: NROOTS
      REAL(KIND=fp) :: HFSCAL,CSCALT
      double precision,parameter :: HALF=0.5D00
      double precision,parameter :: FOUR=4.0D00
      !stuff for digestion
      !INTEGER :: II1,JJ1,KK1,LL1, I2,J2,K2,L2
      INTEGER :: II2,JJ2,KK2,IL,JK,JL
      DOUBLE PRECISION :: XVAL4, VAL,VAL1,VAL4
      REAL(KIND=FP) :: CUTOFFAO

      INTEGER::IJN,KLN,NN,I,J,K,L
      double precision :: N1
      !new
      double precision :: N,NM

      DOUBLE PRECISION :: XIN(31213),YIN(31213),ZIN(31213)
      !added for dirfck
      integer :: i_index,ij_index,ijk_index,itmp,nkl,ijkl_index

      !NSA_rys=9
      !ISA_rys=1
      NROOTS=2

      !write(*,*) "N_RYS_2", N_RYS_2

      DO ICP_rys=1,NCP_rys


!$omp target teams distribute parallel do default(none) &
!$omp shared(fa,da,ia) &
!$omp shared(IDX_RYS_2100,N_RYS_2100) &
!$omp shared(ngth,c) &
!$omp shared(dfttyp) &
!$omp shared(KMAX,KMIN,KTYPE,KLOC,KNG,KSTART,KATOM) &
!$omp shared(EX,CS,CP,CD,CF,CG,CH,CI) &
!$omp shared(NROOTS,ICP_rys) &
!$omp shared(qq4) &
!!$omp private(XIN,YIN,ZIN) &
!$omp private(XIN,YIN,ZIN) &
!$omp private(bx,by,bz,i1,j1,i2,j2)  &
!$omp private(az,rab1,knu,ay,kandl)  &
!$omp private(same,iandj,cpd1,cdd) &
!$omp private(ngtj,ax,ngti,inu,jnu)  &
!$omp private(kl,nx,rcd1,l11,l2,nz)  &
!$omp private(kz,ky,kx,ny,dz,cx,cy)  &
!$omp private(ngtl,lnu,ngtk,dx,dy,k2)  &
!$omp private(cz,k1,ks1,ls1,js1,maxl)  &
!$omp private(is1,llt,nga1,lkt,lit,ljt) &
!$omp private(minl,ll,mini,kk,ii,jj)   &
!$omp private(mink,maxk,maxj,maxi,minj)  &
!$omp private(cdb,gc,cpb1,gb,csb1,gd)  &
!$omp private(csd1,cdc,csc1,cpc1,cda)  &
!$omp private(loci,locj,ngd1,ngb1,ngc1)  &
!$omp private(csa1,cpa1,ga,lock,locl) &
!$omp private(jy,ijx,ijy,ijz,ijgt,jz,jx,ix,ij,iy,iz) &
!$omp private(klx,kly,klz,klgt,lz,lx,ly,ijkl,ik,max,nij) &
!$omp private(ijd,z1,y1,a,ddij,x1,r,ghondo,dij,dkl) &
!!!!added stuff for dirfck
!$omp private(jk,il,jl,kk2,val4,val1,val,xval4,xval1,jj2,ii2,l1,ikk) &
!$omp private(nkl,itmp,ijkl_index,i_index,ij_index,ijk_index) &
!$omp private(nint,cscalt,hfscal)


      DO MM_rys=1,N_RYS_2100(ICP_rys)
      !write(*,*) "i am in rys2"


      ii=IDX_RYS_2100(1,MM_rys,ICP_rys)
      jj=IDX_RYS_2100(2,MM_rys,ICP_rys)
      kk=IDX_RYS_2100(3,MM_rys,ICP_rys)
      ll=IDX_RYS_2100(4,MM_rys,ICP_rys)

      MINI = kmin(ii)
      MAXI = kmax(ii)
      MINJ = kmin(jj)
      MAXJ = kmax(jj)
      MINK = kmin(kk)
      MAXK = kmax(kk)
      MINL = kmin(ll)
      MAXL = kmax(ll)

      IS1= kstart(ii)
      JS1= kstart(jj)
      KS1= kstart(kk)
      LS1= kstart(ll)

      lit= ktype(ii)-1
      !write(*,*) "lit is", lit
      ljt= ktype(jj)-1
      !write(*,*) "ljt is", ljt
      lkt= ktype(kk)-1
      !write(*,*) "lkt is", lkt
      llt= ktype(ll)-1     
      !write(*,*) "llt is", llt

      ! IF (lit.EQ.2 .AND. ljt.EQ.0 .AND. lkt.EQ.0 .AND. llt.EQ.0) then
      ! write(*,*) "rys 2 does 2000"
      ! ENDIF

      ! IF (lit.EQ.2 .AND. ljt.EQ.1 .AND. lkt.EQ.0 .AND. llt.EQ.0) then
      ! write(*,*) "rys 2 does 2100"
      ! ENDIF

      NGA1= kng(ii)
      NGB1= kng(jj)
      NGC1= kng(kk)
      NGD1= kng(ll)

      LOCI = kloc(ii)-MINI
      LOCJ = kloc(jj)-MINJ
      LOCK = kloc(kk)-MINK
      LOCL = kloc(ll)-MINL

      ga=ex(is1)
      csa1=cs(is1)
      cpa1=cp(is1)
      cda=cd(is1)

      gb=ex(js1)
      csb1=cs(js1)
      cpb1=cp(js1)
      cdb=cd(js1)

      gc=ex(ks1)
      csc1=cs(ks1)
      cpc1=cp(ks1)
      cdc=cd(ks1)

      gd=ex(ls1)
      csd1=cs(ls1)
      cpd1=cp(ls1)
      cdd=cd(ls1)

      IANDJ = II .EQ. JJ
      KANDL = KK .EQ. LL
      SAME = (II .EQ. KK) .AND. (JJ .EQ. LL)


      ! IF (KTYPE(ii) .LT. KTYPE(jj)) THEN
      !    INU = jj
      !    JNU = ii
      !    NGTI = NGTH(2)
      !    NGTJ = NGTH(1)
      ! ELSE
      !    INU = ii
      !    JNU = jj
      !    NGTI = NGTH(1)
      !    NGTJ = NGTH(2)
      ! END IF

      ! I = KATOM(INU)
      ! !write(*,*) "I is ", I
      ! AX = C(1,I)
      ! AY = C(2,I)
      ! AZ = C(3,I)
      ! I1 = KSTART(INU)
      ! I2 = I1+KNG(INU)-1
      ! LIT = KTYPE(INU)
      ! MINI = KMIN(INU)
      ! MAXI = KMAX(INU)
      ! LOCI = KLOC(INU)-MINI
      ! NGA1 = 0
      ! !write(*,*) "I1 is ", I1
      ! !write(*,*) "I2 is ", I2
      ! DO I = I1,I2
      !    NGA1 = NGA1+1
      !    GA(NGA1) = EX(I)
      !    CSA1(NGA1) = CS(I)
      !    CPA1(NGA1) = CP(I)
      !    CDA(NGA1) = CD(I)
      ! ENDDO

      ! J = KATOM(JNU)
      ! BX = C(1,J)
      ! BY = C(2,J)
      ! BZ = C(3,J)
      ! J1 = KSTART(JNU)
      ! J2 = J1+KNG(JNU)-1
      ! LJT = KTYPE(JNU)
      ! MINJ = KMIN(JNU)
      ! MAXJ = KMAX(JNU)
      ! LOCJ = KLOC(JNU)-MINJ
      ! NGB1 = 0
      !  DO J = J1,J2
      !    NGB1 = NGB1+1
      !    GB(NGB1) = EX(J)
      !    CSB1(NGB1) = CS(J)
      !    CPB1(NGB1) = CP(J)
      !    CDB(NGB1) = CD(J)
      !  ENDDO
      ! RAB1 = ((AX-BX)*(AX-BX) + (AY-BY)*(AY-BY) + (AZ-BZ)*(AZ-BZ))

      ! IJ = 0
      ! DO I = MINI,MAXI
      !    NX = IX(I)
      !    NY = IY(I)
      !    NZ = IZ(I)
      !    IF (IANDJ) MAXJ = I
      !    DO J = MINJ,MAXJ
      !       IJ = IJ+1
      !       IJX(IJ) = NX+JX(J)
      !       !write(*,*) "IJX(IJ)", IJX(IJ)
      !       IJY(IJ) = NY+JY(J)
      !       IJZ(IJ) = NZ+JZ(J)
      !       IJGT(IJ) = NGTI*(I-MINI)+NGTJ*(J-MINJ)+1
      !    ENDDO
      ! ENDDO



         CALL shells_gpu2_ij(1,ii,jj,.TRUE.,&
              IJGT,IJX,IJY,IJZ,&
              IANDJ,SAME,&
              ii,jj,NGTI,NGTJ,&
              LIT,LJT,LOCI,LOCJ,&
              mini,maxi,minj,maxj,&
              NIJ,IJ,KL,IJKL,&
              NGTH,C,&
              kstart,katom,ktype,kng,kloc,kmin,kmax,&
              ex,cs,cp,cd,&
              ga,csa1,cpa1,cda,&
              gb,csb1,cpb1,cdb,&
              AX,AY,AZ,BX,BY,BZ,RAB1,&
              NGA1,NGB1)

         CALL shells_gpu2_kl(2,ii,jj,kk,ll,.TRUE.,&
              IK,KLGT,KLX,KLY,KLZ,&
              IANDJ,KANDL,SAME,&
              kk,ll,NGTK,NGTL,&
              LKT,LLT,LOCK,LOCL,&
              mink,maxk,minl,maxl,&
              IJ,KL,IJKL,&
              NGTH,C,&
              kstart,katom,ktype,kng,kloc,kmin,kmax,&
              ex,cs,cp,cd,&
              gc,csc1,cpc1,cdc,&
              gd,csd1,cpd1,cdd,&
              CX,CY,CZ,DX,DY,DZ,RCD1,&
              NGC1,NGD1)

         CALL ijprim_gpu_2(ddij,a,r,x1,y1,z1,ijd,&
              IANDJ,KANDL,SAME,&
              LIT,mini,maxi,minj,maxj,NIJ,&
              ga,csa1,cpa1,cda,&
              gb,csb1,cpb1,cdb,&
              AX,AY,AZ,BX,BY,BZ,RAB1,&
              NGA1,NGB1)

         call zqout_gpu_2(ghondo,IJGT,KLGT,&
              IANDJ,KANDL,SAME,&
              mini,maxi,minj,maxj,mink,maxk,minl,maxl)

         CALL GENRAL_GPU_2(ghondo,ddij,dkl,dij,&
             IJGT,IJX,IJY,IJZ,IK,KLGT,KLX,KLY,KLZ, &
             A,R,X1,Y1,Z1,IJD,iandj,kandl,same,&
             LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,&
             qq4,mini,maxi,minj,maxj,mink,maxk,minl,maxl,&
             NIJ,IJ,KL,IJKL,&
             ga,csa1,cpa1,cda,&
             gb,csb1,cpb1,cdb,&
             gc,csc1,cpc1,cdc,&
             gd,csd1,cpd1,cdd,&
             AX,AY,AZ,BX,BY,BZ,RAB1,CX,CY,CZ,DX,DY,DZ,RCD1,&
             NGA1,NGB1,NGC1,NGD1,&
             XIN,YIN,ZIN)

         HFSCAL=DFTTYP(3)
         CSCALT=1.0D+00

         ! CALL dirfck_rhf_rys_2(IA,DA,FA,ii,jj,kk,ll, &
         !      GHONDO,NGTI,NGTJ,NGTK,NGTL, &
         !      HFSCAL,CSCALT,NINT,&
         !      kstart,katom,kng,kloc,kmin,kmax)
    
      ! SAME  = ii.EQ.kk.AND.jj.EQ.ll
      ! IANDJ = ii.EQ.jj
      ! KANDL = kk.EQ.ll

      ! MINI = KMIN(ii)
      ! MINJ = KMIN(jj)
      ! MINK = KMIN(kk)
      ! MINL = KMIN(ll)
      ! MAXI = KMAX(ii)
      ! MAXJ = KMAX(jj)
      ! MAXK = KMAX(kk)
      ! MAXL = KMAX(ll)
      ! LOCI = KLOC(ii)-MINI
      ! LOCJ = KLOC(jj)-MINJ
      ! LOCK = KLOC(kk)-MINK
      ! LOCL = KLOC(ll)-MINL

      XVAL1 = HFSCAL
      XVAL4 = FOUR*CSCALT

      NIJ = 0
      I_INDEX = 1
      DO I = MINI,MAXI
         IF (IANDJ) MAXJ = I

         I1 = I+LOCI

         IJ_INDEX = I_INDEX
         !write(*,*) "IJ_INDEX", IJ_INDEX
         I_INDEX = I_INDEX + NGTI

         DO J = MINJ,MAXJ
            NIJ = NIJ+1

            J1 = J+LOCJ
            I2 = I1
            J2 = J1
            IF (I1.LT.J1) THEN ! SORT <IJ|
               I2 = J1
               J2 = I1
            ENDIF

            IJK_INDEX = IJ_INDEX
            IJ_INDEX = IJ_INDEX + NGTJ

            NKL = NIJ

            DO K =  MINK,MAXK
               IF (KANDL) MAXL = K

               K1 = K + LOCK

               IF(SAME) THEN ! ACCOUNT FOR NON-UNIQUE PERMUTATIONS
                  ITMP = MIN(MAXL-MINL+1,NKL)
                  IF (ITMP.EQ.0) CYCLE
                  MAXL = MINL + ITMP - 1
                  NKL = NKL - ITMP
               ENDIF

               IJKL_INDEX = IJK_INDEX
               IJK_INDEX = IJK_INDEX + NGTK

               DO L=MINL,MAXL

                  VAL = GHONDO( IJKL_INDEX )
                  !IF(VAL .NE. 0.00) write(*,*) "val is", VAL
                  !write(*,*) "val is", VAL
                  IJKL_INDEX = IJKL_INDEX + NGTL
                  !IF(ABS(VAL).LT.CUTOFF) CYCLE ! GOTO 300
                  NINT = NINT + 1

                  L1 = L + LOCL
                  K2 = K1
                  L2 = L1

                  IF (K2.LT.L2) THEN ! SORT |KL>
                     K2 = L1
                     L2 = K1
                  ENDIF

                  II = I2
                  JJ = J2
                  KK = K2
                  LL = L2

                  IF (II.LT.KK) THEN ! SORT <IJ|KL>
                     II = K2
                     JJ = L2
                     KK = I2
                     LL = J2
                  ELSE IF (II.EQ.KK.AND.JJ.LT.LL) THEN ! SORT <IJ|IL>
                     JJ = L2
                     LL = J2
                  ENDIF

                  II2 = IA(II)
                  JJ2 = IA(JJ)
                  KK2 = IA(KK)

                  IJ = II2 + JJ
                  IKK = II2 + KK
                  IL = II2 + LL
                  JK = JJ2 + KK
                  JL = JJ2 + LL
                  KL = KK2 + LL
                  IF (JJ.LT.KK) JK = KK2 + JJ
                  IF (JJ.LT.LL) JL = IA(LL) + JJ
!
!       ACCOUNT FOR IDENTICAL PERMUTATIONS.
!
                  IF(II.EQ.JJ) VAL = VAL*HALF
                  IF(KK.EQ.LL) VAL = VAL*HALF
                  IF(II.EQ.KK.AND.JJ.EQ.LL) VAL = VAL*HALF
                  VAL1 = VAL*XVAL1
                  VAL4 = VAL*XVAL4
!$omp atomic
                  FA(IJ) = FA(IJ) + VAL4*DA(KL)
!$omp atomic
                  FA(KL) = FA(KL) + VAL4*DA(IJ)
!$omp atomic
                  FA(IKK) = FA(IKK) - VAL1*DA(JL)
!$omp atomic
                  FA(JL) = FA(JL) - VAL1*DA(IKK)
!$omp atomic
                  FA(IL) = FA(IL) - VAL1*DA(JK)
!$omp atomic
                  FA(JK) = FA(JK) - VAL1*DA(IL)
               ENDDO
            ENDDO
         ENDDO
      ENDDO


      ENDDO
!$omp end target teams distribute parallel do
  ENDDO
  !ENDDO

  END SUBROUTINE rysint_gpu_2100


SUBROUTINE rysint_gpu_3(ii,jj,kk,ll,ktype,ghondo,ddij,ngth,c,&
                        kstart,katom,kng,kloc,kmin,kmax,&
                        ex,cs,cp,cd,cf,cg,ch,ci,&
                        NGTI,NGTJ,NGTK,NGTL,&
                        NIJ,IJ,KL,IJKL,&
                        qq4,&
                        AX,AY,AZ,BX,BY,BZ,RAB,CX,CY,CZ,DX,DY,DZ,RCD1,&
                        dfttyp,&
                        ia,da,fa,nint,&
                        l1,l2a,nflmat,&
                        IDX_RYS_3,N_RYS_3,NCP_rys, &
                        XIN,YIN,ZIN)

    USE omp_lib
    USE prec, ONLY: fp
    USE mx_limits, ONLY: &
    mxsh,mxg2,mxatm,mxgtot,mxgsh

    !INTEGER, INTENT(IN) :: &
    !l1, l2a,nflmat, ia(l1)

    !for dirfck
    REAL(KIND=fp),INTENT(OUT) :: fa(l2a*nflmat) !FOCK MATRIX
    REAL(KIND=fp),INTENT(INOUT) :: da(l2a) !DENSITY MATRIX
    INTEGER,INTENT(IN) :: ia(l1) !TRIANGULAR MATRIX 
    INTEGER :: l2a,nflmat,l1,nint 

    REAL(KIND=fp),DIMENSION(20) :: DFTTYP

    INTEGER, INTENT(IN) :: ktype(mxsh)
    INTEGER :: ii,jj,kk,ll

    REAL(KIND=fp),INTENT(INOUT),ALLOCATABLE :: &
      ghondo(:)

    REAL(KIND=fp),INTENT(INOUT) :: &
      ddij(:)

    REAL(KIND=fp),DIMENSION(784) :: DKL,DIJ
    REAL(KIND=fp),DIMENSION(784) :: &
     IJGT,IJX,IJY,IJZ,IK,KLGT,KLX,KLY,KLZ

    REAL(KIND=fp),DIMENSION(MXG2) :: A,R,X1,Y1,Z1
    REAL(KIND=fp),DIMENSION(784) :: IJD

      !for shlexc
      INTEGER,DIMENSION(4) :: NGTH
      !for infoa
      REAL(KIND=fp) :: C(3,mxatm)
      !for roots
      REAL(KIND=fp),DIMENSION(13) :: U,W
      REAL(KIND=fp) :: XX
      !for eriout
      INTEGER :: NGTI,NGTJ,NGTK,NGTL
      !shlnos
      INTEGER ::mini,maxi,minj,maxj,mink,maxk,minl,maxl,&
                LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,&
                NIJ,IJ,KL,IJKL
      REAL(KIND=fp) :: qq4  
      !nshell
      INTEGER,DIMENSION(mxsh) :: kstart,katom,kng,kloc,kmin,kmax,nshell
      REAL(KIND=fp),DIMENSION(mxgtot) :: ex,cs,cp,cd,cf,cg,ch,ci

      !shlinf
      double precision :: ga(mxgsh),csa1(mxgsh),cpa1(mxgsh),cda(mxgsh)
      double precision :: cfa(mxgsh),cga(mxgsh),cha(mxgsh),cia(mxgsh)
      double precision :: gb(mxgsh),csb1(mxgsh),cpb1(mxgsh),cdb(mxgsh)
      double precision :: cfb(mxgsh),cgb(mxgsh),chb(mxgsh),cib(mxgsh)
      double precision :: gc(mxgsh),csc1(mxgsh),cpc1(mxgsh),cdc(mxgsh)
      double precision :: cfc(mxgsh),cgc(mxgsh),chc(mxgsh),cic(mxgsh)
      double precision :: gd(mxgsh),csd1(mxgsh),cpd1(mxgsh),cdd(mxgsh)
      double precision :: cfd(mxgsh),cgd(mxgsh),chd(mxgsh),cid(mxgsh)

      double precision :: AX,AY,AZ,BX,BY,BZ,RAB,CX,CY,CZ,DX,DY,DZ,RCD1
      integer :: NGA1,NGB1,NGC1,NGD1
      LOGICAL IANDJ,KANDL,SAME
      integer :: IS1,JS1,KS1,LS1
      !IS1=kstart(ii)
      !rys
      INTEGER,PARAMETER :: NUM_QUART=9169597 !number of quartets
      INTEGER :: IDX_RYS_3(4,NUM_QUART,NCP_rys)
      INTEGER :: N_RYS_3(NCP_rys)

      INTEGER :: NSA_rys,ISA_rys !sum angular momentum
      INTEGER :: NCP_rys,ICP_rys !product contraction
      !INTEGER,ALLOCATABLE :: RYSINDEX(:,:,:,:) !rys shell index
      !INTEGER,ALLOCATABLE :: RYSTYP(:,:) !rystype
      INTEGER :: MAXRYS
      iNTEGER :: MM_rys
      integer :: NROOTS

      DOUBLE PRECISION :: XIN(31213),YIN(31213),ZIN(31213)

      
      NSA_rys=9
      ISA_rys=1
      NROOTS=3
      !DO ISA_rys=1,NSA_rys
      !write(*,*) "NROOTS in 3", NROOTS
      DO ICP_rys=1,NCP_rys

! $omp target teams distribute parallel do default(none) &
! !$omp target teams

! !$omp parallel &
! $omp shared(fa,da,ia) &
! $omp shared(IDX_RYS_3,N_RYS_2) &
! $omp shared(ngth,c) &
! $omp shared(dfttyp) &
! $omp shared(KMAX,KMIN,KTYPE,KLOC,KNG,KSTART,KATOM) &
! $omp shared(EX,CS,CP,CD,CF,CG,CH,CI) &
! $omp shared(NROOTS,ICP_rys) &
! $omp shared(qq4,nshell) &
! $omp private(XIN,YIN,ZIN) &
! $omp private(bx,by,bz,i1,j1,i2,j2)  &
! $omp private(az,rab1,knu,ay,kandl)  &
! $omp private(same,iandj,cpd1,cdd) &
! $omp private(ngtj,ax,ngti,inu,jnu)  &
! $omp private(kl,nx,rcd1,l11,l2,nz)  &
! $omp private(kz,ky,kx,ny,dz,cx,cy)  &
! $omp private(ngtl,lnu,ngtk,dx,dy,k2)  &
! $omp private(cz,k1,ks1,ls1,js1,maxl)  &
! $omp private(is1,llt,nga1,lkt,lit,ljt) &
! $omp private(minl,ll,mini,kk,ii,jj)   &
! $omp private(mink,maxk,maxj,maxi,minj)  &
! $omp private(cdb,gc,cpb1,gb,csb1,gd)  &
! $omp private(csd1,cdc,csc1,cpc1,cda)  &
! $omp private(loci,locj,ngd1,ngb1,ngc1)  &
! $omp private(csa1,cpa1,ga,lock,locl) &
! $omp private(jy,ijx,ijy,ijz,ijgt,jz,jx,ix,ij,iy,iz) &
! $omp private(klx,kly,klz,klgt,lz,lx,ly,ijkl,ik,max,nij) &
! $omp private(ijd,z1,y1,a,ddij,x1,r,ghondo,dij,dkl) &
! $omp private(nint,cscalt,hfscal)

      DO MM_rys=1,N_RYS_3(ICP_rys)

      ii=IDX_RYS_3(1,MM_rys,ICP_rys)
      jj=IDX_RYS_3(2,MM_rys,ICP_rys)
      kk=IDX_RYS_3(3,MM_rys,ICP_rys)
      ll=IDX_RYS_3(4,MM_rys,ICP_rys)

      MINI = kmin(ii)
      MAXI = kmax(ii)
      MINJ = kmin(jj)
      MAXJ = kmax(jj)
      MINK = kmin(kk)
      MAXK = kmax(kk)
      MINL = kmin(ll)
      MAXL = kmax(ll)

      IS1= kstart(ii)
      JS1= kstart(jj)
      KS1= kstart(kk)
      LS1= kstart(ll)

      lit= ktype(ii)-1
      ljt= ktype(jj)-1
      lkt= ktype(kk)-1
      llt= ktype(ll)-1

      NGA1= kng(ii)
      NGB1= kng(jj)
      NGC1= kng(kk)
      NGD1= kng(ll)

      LOCI = kloc(ii)-MINI
      LOCJ = kloc(jj)-MINJ
      LOCK = kloc(kk)-MINK
      LOCL = kloc(ll)-MINL

      ga=ex(is1)
      csa1=cs(is1)
      cpa1=cp(is1)
      cda=cd(is1)
      cfa=cf(is1)
      cga=cg(is1)
      cha=ch(is1)
      cia=ci(is1)

      gb=ex(js1)
      csb1=cs(js1)
      cpb1=cp(js1)
      cdb=cd(js1)
      cfb=cf(js1)
      cgb=cg(js1)
      chb=ch(js1)
      cib=ci(js1)

      gc=ex(ks1)
      csc1=cs(ks1)
      cpc1=cp(ks1)
      cdc=cd(ks1)
      cfc=cf(ks1)
      cgc=cg(ks1)
      chc=ch(ks1)
      cic=ci(ks1)

      gd=ex(ls1)
      csd1=cs(ls1)
      cpd1=cp(ls1)
      cdd=cd(ls1)
      cfd=cf(ls1)
      cgd=cg(ls1)
      chd=ch(ls1)
      cid=ci(ls1)

      IANDJ = II .EQ. JJ
      KANDL = KK .EQ. LL
      SAME = (II .EQ. KK) .AND. (JJ .EQ. LL)

!write(*,*) "calling shells gpu in 3"
!   fill indices for Rys quadrature
   CALL shells_gpu(1,ii,jj,kk,ll,.TRUE.,&
IJGT,IJX,IJY,IJZ,IK,KLGT,KLX,KLY,KLZ,&
IANDJ,KANDL,SAME,&
ii,jj,kk,ll,NGTI,NGTJ,NGTK,NGTL,&
LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,&
mini,maxi,minj,maxj,mink,maxk,minl,maxl,&
NIJ,IJ,KL,IJKL,&
NROOTS,NGTH,C,&
kstart,katom,ktype,kng,kloc,kmin,kmax,nshell,&
ex,cs,cp,cd,cf,cg,ch,ci,&
ga,csa1,cpa1,cda,cfa,cga,cha,cia,&
gb,csb1,cpb1,cdb,cfb,cgb,chb,cib,&
gc,csc1,cpc1,cdc,cfc,cgc,chc,cic,&
gd,csd1,cpd1,cdd,cfd,cgd,chd,cid,&
AX,AY,AZ,BX,BY,BZ,RAB,CX,CY,CZ,DX,DY,DZ,RCD1,&
NGA1,NGB1,NGC1,NGD1)

   CALL shells_gpu(2,ii,jj,kk,ll,.TRUE.,&
IJGT,IJX,IJY,IJZ,IK,KLGT,KLX,KLY,KLZ,&
IANDJ,KANDL,SAME,&
ii,jj,kk,ll,NGTI,NGTJ,NGTK,NGTL,&
LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,&
mini,maxi,minj,maxj,mink,maxk,minl,maxl,&
NIJ,IJ,KL,IJKL,&
NROOTS,NGTH,C,&
kstart,katom,ktype,kng,kloc,kmin,kmax,nshell,&
ex,cs,cp,cd,cf,cg,ch,ci,&
ga,csa1,cpa1,cda,cfa,cga,cha,cia,&
gb,csb1,cpb1,cdb,cfb,cgb,chb,cib,&
gc,csc1,cpc1,cdc,cfc,cgc,chc,cic,&
gd,csd1,cpd1,cdd,cfd,cgd,chd,cid,&
AX,AY,AZ,BX,BY,BZ,RAB,CX,CY,CZ,DX,DY,DZ,RCD1,&
NGA1,NGB1,NGC1,NGD1)

   CALL ijprim_gpu(ddij,a,r,x1,y1,z1,ijd,&
IANDJ,KANDL,SAME,&
LIT,mini,maxi,minj,maxj,mink,maxk,minl,maxl,NIJ,&
ga,csa1,cpa1,cda,cfa,cga,cha,cia,&
gb,csb1,cpb1,cdb,cfb,cgb,chb,cib,&
gc,csc1,cpc1,cdc,cfc,cgc,chc,cic,&
gd,csd1,cpd1,cdd,cfd,cgd,chd,cid,&
AX,AY,AZ,BX,BY,BZ,RAB,CX,CY,CZ,DX,DY,DZ,RCD1,&
NGA1,NGB1,NGC1,NGD1)

   call zqout_gpu(ghondo,&
IJGT,IJX,IJY,IJZ,IK,KLGT,KLX,KLY,KLZ,&
IANDJ,KANDL,SAME,&
mini,maxi,minj,maxj,mink,maxk,minl,maxl)

! calculate integrals
     CALL genral_gpu_3(ghondo,ddij,dkl,dij,&
IJGT,IJX,IJY,IJZ,IK,KLGT,KLX,KLY,KLZ, &
A,R,X1,Y1,Z1,IJD,iandj,kandl,same,&
LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,&
qq4,mini,maxi,minj,maxj,mink,maxk,minl,maxl,&
NIJ,IJ,KL,IJKL,&
NROOTS,&
ga,csa1,cpa1,cda,cfa,cga,cha,cia,&
gb,csb1,cpb1,cdb,cfb,cgb,chb,cib,&
gc,csc1,cpc1,cdc,cfc,cgc,chc,cic,&
gd,csd1,cpd1,cdd,cfd,cgd,chd,cid,&
AX,AY,AZ,BX,BY,BZ,RAB,CX,CY,CZ,DX,DY,DZ,RCD1,&
NGA1,NGB1,NGC1,NGD1)

     CALL dirfck_gpu_rys(ia,da,fa,ghondo,nint,&
kstart,katom,kng,kloc,kmin,kmax,&
ii,jj,kk,ll,&
NGTI,NGTJ,NGTK,NGTL,dfttyp)

  ! CALL dirfck_rhf_rys_3(IA,DA,FA,ii,jj,kk,ll, &
  !                 GHONDO,NGTI,NGTJ,NGTK,NGTL, &
  !                 HFSCAL,CSCALT,NINT,&
  !                 kstart,katom,kng,kloc,kmin,kmax)

  ENDDO
!!$omp end target teams distribute parallel do
  ENDDO
  !ENDDO

  END SUBROUTINE rysint_gpu_3



!---------------------------------------------------------------------
!=====================================================================
!---------------------------------------------------------------------

!*MODULE OMPMOD   *DECK RAXINTSP
!>    @brief   Rotated axis method ERI calculation for S,P,L shells
!>
!>    @details Calculates two-electron integrals over
!>             shell quartet \f$ (ii,jj|kk,ll) \f$ using Pople algorithm
!>             implemented in GAMESS (`GENR70` and `GENR03` routines).
!>             Supports S, P and L shells.
!>             Based on `SHELLQUART` routine from `int2a.src`
!>
!>    @author  Vladimir Mironov
!
!     REVISION HISTORY:
!>    @date _May, 2016_ Initial release
!>    @date _Jun, 2017_ Bug fixes
!
!     PARAMETERS:
!
!>    @param[in]     ii        first shell index
!>    @param[in]     jj        second shell index
!>    @param[in]     kk        third shell index
!>    @param[in]     ll        fourth shell index
!>    @param[out]    ghondo(:) array of integrals over shell quartet
  SUBROUTINE ompmod_raxintsp(ish,jsh,ksh,lsh,ghondo)


    INTEGER, INTENT(IN) :: &
      ish,jsh,ksh,lsh

    REAL(KIND=fp),INTENT(INOUT),ALLOCATABLE :: &
      ghondo(:)

    COMMON /gout  / gpople(768),norgp
        REAL(KIND=fp) :: gpople
        INTEGER :: norgp

    COMMON /popout/ lpopi,lpopj,lpopk,lpopl
        INTEGER :: lpopi,lpopj,lpopk,lpopl

    COMMON /shlg70/ ipl,jpl,kpl,lpl,inew,jnew,knew,lnew
        INTEGER :: ipl,jpl,kpl,lpl,inew,jnew,knew,lnew

    COMMON /eridat/ len1,len2,len3,len4
        INTEGER :: len1,len2,len3,len4

    COMMON /eriout/ inw,jnw,knw,lnw,lstri,lstrj,lstrk,lstrl
        INTEGER :: inw,jnw,knw,lnw,lstri,lstrj,lstrk,lstrl

! required for thread-safety:
!$omp threadprivate(/gout  /,/popout/,/shlg70/,/eriout/)

    COMMON /nshel / ex(mxgtot),cs(mxgtot),cp(mxgtot),cd(mxgtot),    &
                    cf(mxgtot),cg(mxgtot),ch(mxgtot),ci(mxgtot),    &
                    kstart(mxsh),katom(mxsh),ktype(mxsh),kng(mxsh), &
                    kloc(mxsh),kmin(mxsh),kmax(mxsh),nshell
        REAL(KIND=fp) :: ex,cs,cp,cd,cf,cg,ch,ci
        INTEGER :: kstart,katom,ktype,kng,kloc,kmin,kmax,nshell

    INTEGER :: &
      mini,maxi,minj,maxj,mink,maxk,minl,maxl, &
      ii,i,ij,ijk,ijkl,ip,ijp,ijkp,ijklp,j,k,l

    ipl = ish
    jpl = jsh
    kpl = ksh
    lpl = lsh
    inw = ish
    jnw = jsh
    knw = ksh
    lnw = lsh

    CALL genr70(1,.FALSE.)

!  save to output array with HONDO indexing

    mini = kmin(inw)
    maxi = kmax(inw)
    minj = kmin(jnw)
    maxj = kmax(jnw)
    mink = kmin(knw)
    maxk = kmax(knw)
    minl = kmin(lnw)
    maxl = kmax(lnw)

    ii = 1
    DO i = mini, maxi
      ip = (i-1)*lpopi + 1
      ij  = ii
      DO j = minj, maxj
        ijp = (j-1)*lpopj + ip
        ijk  = ij
        DO k = mink, maxk
          ijkp = (k-1)*lpopk + ijp
          ijkl  = ijk
          DO l = minl, maxl
            ijklp = (l-1)*lpopl + ijkp
            ghondo(ijkl) = gpople(ijklp)
            ijkl = ijkl  + len1
          ENDDO
          ijk  = ijk  + len2
        ENDDO
        ij  = ij  + len3
      ENDDO
      ii = ii + len4
    ENDDO


    lstri = len4
    lstrj = len3
    lstrk = len2
    lstrl = len1
  END SUBROUTINE ompmod_raxintsp


      SUBROUTINE raxintsp_cpyint_gpu &
                (lstri,lstrj,lstrk,lstrl, &
                 ghondo,gpople,maxg, &
                 lpopi,lpopj,lpopk,lpopl, &
                 mini,maxi,minj,maxj,mink,maxk,minl,maxl)
 
      REAL(KIND=fp),INTENT(INOUT) :: ghondo(maxg)
      INTEGER :: maxg
	   REAL(KIND=fp) :: gpople(768)

      INTEGER :: lpopi,lpopj,lpopk,lpopl
  
      COMMON /eridat/ len1,len2,len3,len4
      INTEGER :: len1,len2,len3,len4
  
      INTEGER :: lstri,lstrj,lstrk,lstrl      

      INTEGER :: &
        mini,maxi,minj,maxj,mink,maxk,minl,maxl, &
        ii,i,ij,ijk,ijkl,ip,ijp,ijkp,ijklp,j,k,l


      ii = 1
      DO i = mini, maxi
        ip = (i-1)*lpopi + 1
        ij  = ii
        DO j = minj, maxj
          ijp = (j-1)*lpopj + ip
          ijk  = ij
          DO k = mink, maxk
            ijkp = (k-1)*lpopk + ijp
            ijkl  = ijk
            DO l = minl, maxl
              ijklp = (l-1)*lpopl + ijkp
              ghondo(ijkl) = gpople(ijklp)
              ijkl = ijkl  + len1
            ENDDO
            ijk  = ijk  + len2
          ENDDO
          ij  = ij  + len3
        ENDDO
        ii = ii + len4
      ENDDO


      lstri = len4
      lstrj = len3
      lstrk = len2
      lstrl = len1
      END SUBROUTINE raxintsp_cpyint_gpu




!---------------------------------------------------------------------
!=====================================================================
!---------------------------------------------------------------------

!*MODULE OMPMOD   *DECK OMPMOD_RAXINTSPD
!>    @brief   Rotated axis method ERI calculation for \f$ L\le 3 \f$ shells
!>
!>    @details Calculates two-electron integrals over
!>             shell quartet \f$ (ii,jj|kk,ll) \f$ using Pople algorithm
!>             implemented in GAMESS (`GENR70` and `GENR03` routines).
!>             Supports S, P, D and L shells.
!>             Based on `SHELLQUART` routine from `int2a.src`
!>
!>    @author  Vladimir Mironov
!
!     REVISION HISTORY:
!>    @date _May, 2016_ Initial release
!>    @date _Jun, 2017_ Bug fixes
!
!     PARAMETERS:
!
!>    @param[in]     ii        first shell index
!>    @param[in]     jj        second shell index
!>    @param[in]     kk        third shell index
!>    @param[in]     ll        fourth shell index
!>    @param[out]    ghondo(:) array of integrals over shell quartet
  SUBROUTINE ompmod_raxintspd(ish,jsh,ksh,lsh,ghondo)

    INTEGER, INTENT(IN) :: &
      ish,jsh,ksh,lsh

    REAL(KIND=fp),INTENT(INOUT),ALLOCATABLE :: &
      ghondo(:)

    COMMON /flips / ib(4,3)
        INTEGER :: ib

    COMMON /shlg70/ ipl,jpl,kpl,lpl,inew,jnew,knew,lnew
        INTEGER :: ipl,jpl,kpl,lpl,inew,jnew,knew,lnew

    COMMON /eridat/ len1,len2,len3,len4
        INTEGER :: len1,len2,len3,len4

    COMMON /eriout/ inw,jnw,knw,lnw,lstri,lstrj,lstrk,lstrl
        INTEGER :: inw,jnw,knw,lnw,lstri,lstrj,lstrk,lstrl

! required for thread-safety:
!$omp threadprivate(/flips /,/eriout/,/shlg70/)

    COMMON /nshel / ex(mxgtot),cs(mxgtot),cp(mxgtot),cd(mxgtot),    &
                    cf(mxgtot),cg(mxgtot),ch(mxgtot),ci(mxgtot),    &
                    kstart(mxsh),katom(mxsh),ktype(mxsh),kng(mxsh), &
                    kloc(mxsh),kmin(mxsh),kmax(mxsh),nshell
        REAL(KIND=fp) :: ex,cs,cp,cd,cf,cg,ch,ci
        INTEGER :: kstart,katom,ktype,kng,kloc,kmin,kmax,nshell

    COMMON /output/ nprint,itol,icut,normf,normp,nopk
        INTEGER :: nprint,itol,icut,normf,normp,nopk

    REAL(KIND=fp) :: &
      grotspd(1296)

    LOGICAL :: &
      iandj,kandl,same

    INTEGER :: &
      mini,maxi,minj,maxj,mink,maxk,minl,maxl, &
      i,j,k,l,ijn,kln, &
      ibb, jbb, kbb, lbb, jmax, lmax, &
      ihondo,irotax,ijklhondo,ijklrotax, &
      ijhondo, ijrotax, ijkhondo, ijkrotax, &
      idpop(4,10)

    DATA idpop /0,0,0,0,216,36,6,1,432,72,12,2,648,108,18,3, &
                0,0,0,0,216,36,6,1,432,72,12,2,648,108,18,3, &
                864,144,24,4,1080,180,30,5/

    ipl = ish
    jpl = jsh
    kpl = ksh
    lpl = lsh

    inw = ish
    jnw = jsh
    knw = ksh
    lnw = lsh

    CALL genr03(grotspd)

!   save to output array with HONDO indexing

    iandj = ish.EQ.jsh
    kandl = ksh.EQ.lsh
    same  = ish.EQ.ksh  .AND.  jsh.EQ.lsh
    IF(nopk.EQ.0) same = .FALSE.

    ibb = ib(1,1)
    jbb = ib(2,1)
    kbb = ib(3,1)
    lbb = ib(4,1)

    mini = kmin(inw)
    maxi = kmax(inw)
    minj = kmin(jnw)
    maxj = kmax(jnw)
    mink = kmin(knw)
    maxk = kmax(knw)
    minl = kmin(lnw)
    maxl = kmax(lnw)

    ijn = 0
    jmax = maxj

ic: DO i = mini, maxi

      ihondo = (i-mini)*len4 + 1
      irotax = idpop(ibb,i)  + 1

      IF(iandj) jmax=i

jc:   DO j = minj, jmax

        ijhondo = (j-minj)*len3 + ihondo
        ijrotax = idpop(jbb,j)  + irotax
        ijn = ijn+1
        lmax=maxl
        kln=0

kc:     DO k = mink, maxk

          ijkhondo = (k-mink)*len2 + ijhondo
          ijkrotax = idpop(kbb,k)  + ijrotax

          IF(kandl) lmax=k

lc:       DO l = minl, lmax

            kln = kln+1

            IF(same .AND. kln.GT.ijn) EXIT kc

            ijklhondo = (l-minl)*len1 + ijkhondo
            ijklrotax = idpop(lbb,l)  + ijkrotax
            ghondo(ijklhondo) = grotspd(ijklrotax)

          ENDDO lc
        ENDDO kc
      ENDDO jc
    ENDDO ic


    lstri = len4
    lstrj = len3
    lstrk = len2
    lstrl = len1

  END SUBROUTINE ompmod_raxintspd

!---------------------------------------------------------------------
!=====================================================================
!---------------------------------------------------------------------

!*MODULE OMPMOD   *DECK OMPMOD_SHELLQUART
!>    @brief   Thread-safe analog of the `SHELLQUART` subroutine for
!>             2e integrals algorithm selection
!>
!>    @details Calculates two-electron integrals over
!>     shell quartet \f$ (ii,jj|kk,ll) \f$
!>     Currently implements `INTTYP=0` ("best") scheme:
!>     - LMAX = 2 (SP(L)D) shells : rotated axis code
!>     - LMAX > 2, uncontracted shells : Rys quadrature
!>     - LMAX > 2, no L-type shells, L(i)+L(j)+L(k)+L(l)<=5 : ERIC code
!>     - otherwise : Rys quadrature
!>
!>    @author  Vladimir Mironov
!
!     REVISION HISTORY:
!>    @date _Jan, 2019_ Initial release
!
!     PARAMETERS:
!
!>    @param[in]     ii        first shell index
!>    @param[in]     jj        second shell index
!>    @param[in]     kk        third shell index
!>    @param[in]     ll        fourth shell index
!>    @param[out]    ghondo(:) array of integrals over shell quartet
!>    @param[out]    ddij(:)   temporary array for 'densities' of integrals
      SUBROUTINE ompmod_shellquart(ii,jj,kk,ll,ghondo,ddij)

      INTEGER, INTENT(IN) :: ii, jj, kk, ll
      REAL(KIND=fp), ALLOCATABLE, INTENT(INOUT) :: ghondo(:), ddij(:)

      COMMON /nshel / ex(mxgtot),cs(mxgtot),cp(mxgtot),cd(mxgtot),    &
                    cf(mxgtot),cg(mxgtot),ch(mxgtot),ci(mxgtot),    &
                    kstart(mxsh),katom(mxsh),ktype(mxsh),kng(mxsh), &
                    kloc(mxsh),kmin(mxsh),kmax(mxsh),nshell
      REAL(KIND=fp) :: ex,cs,cp,cd,cf,cg,ch,ci
      INTEGER :: kstart,katom,ktype,kng,kloc,kmin,kmax,nshell

      INTEGER :: shlmxang, shlsumang, shltotcon
      LOGICAL :: lshell

      shlmxang  = max(ktype(ii), ktype(jj), ktype(kk), ktype(ll))

      shlsumang = ktype(ii) + ktype(jj) + ktype(kk) + ktype(ll) - 4

      shltotcon = kng(ii)*kng(jj)*kng(kk)*kng(ll)

      lshell = (kmax(ii)-kmin(ii)).EQ.3 .OR. &
               (kmax(jj)-kmin(jj)).EQ.3 .OR. &
               (kmax(kk)-kmin(kk)).EQ.3 .OR. &
               (kmax(ll)-kmin(ll)).EQ.3

      IF (shlmxang.LE.2) THEN
          CALL ompmod_raxintsp(ii, jj, kk, ll, ghondo)
      ELSE
          IF (shltotcon.NE.1) THEN
              IF (shlmxang.EQ.3) THEN
                  CALL ompmod_raxintspd(ii, jj, kk, ll, ghondo)
              ELSE IF (shlsumang.LE.5 .AND. .NOT.lshell.AND.shlmxang.LE.5) THEN
                  CALL eric_ts(ii, jj, kk, ll, ghondo)
              ELSE
                  CALL ompmod_rysint(ii, jj, kk, ll, ktype, ghondo, ddij)
              ENDIF
          ELSE
              CALL ompmod_rysint(ii, jj, kk, ll, ktype, ghondo, ddij)
          ENDIF
      ENDIF

      END SUBROUTINE ompmod_shellquart







END MODULE ompmod



!
!*MODULE SCFLIB  *DECK DIRFCK
!>
!>    @brief form fock operator directly from integrals
!>
!>    @author Frank Jensen, MWS
!>
!>    @date January, 2017 - C.Bertoni
!>          - Added changes for EFMO gradient
!>
!>    @param scftyp : [in] scf type 8 byte hollerith constant.
!>    @param ia : [in] triangular index addresses
!>    @param da : [in] rhf or rohf/uhf alpha density matrix.
!>    @param fa : [in,out] rhf or rohf/uhf alpha fock matrix.
!>    @param db : [in] rohf/uhf beta density matrix.
!>    @param fb : [in,out] rohf/uhf beta fock matrix.
!>    @param ghondo : [in] hondo indexed integrals.
!>    @param l2 : [in] size of symmetric packed density or fock matrix.
!>    @param nint : [in,out] number of non-zero integrals.
!>    @param nxyz : [in] nxyz > 1 selects cphf calculation.
!>           nxyz < 0 selects code for efmo dispersion gradient.
!>


      SUBROUTINE dirfck_gpu(SCFTYP,IA,DA,FA,DB,FB,GHONDO,L2,NINT,NXYZ)

      USE camdft, ONLY: CAMFLAG
      USE lrcdft, ONLY: LCFLAG, EMU, EMU2, LRFILE
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION IA(*),DA(*),FA(*),DB(*),FB(*),GHONDO(*)

      LOGICAL RHF,CPHF,UROHF,GVB,CPUHF, dyndisp
      LOGICAL LRINT,OUT

      COMMON /DFTPAR/ DFTTYP(20),EXENA,EXENB,EXENC, &
                      IDFT34,NAUXFUN,NAUXSHL
!$omp threadprivate(/DFTPAR/)                      
      COMMON /ERIOUT/ ISH,JSH,KSH,LSH,LSTRI,LSTRJ,LSTRK,LSTRL
!$omp threadprivate(/ERIOUT/)      
      COMMON /NLRCF / LRINT
!$omp threadprivate(/NLRCF /)      
      COMMON /SHLT  / TOL,CUTOFFAO,ICOUNT,OUT

      DATA DHRHF,DHUHF,DHROHF /8HRHF     ,8HUHF     ,8HROHF     /
      DATA DHGVB /8HGVB     /
!
!     ----- FORM FOCK OPERATOR DIRECTLY FROM INTEGRALS -----
!     THIS ROUTINE WAS PIECED TOGETHER FROM QOUT AND HSTAR
!     BY FRANK JENSEN AT ODENSE UNIVERSITY IN MARCH 1990.
!     SCF FUNCTIONS BESIDES RHF ADDED BY MWS IN AUGUST 1991.
!
!     NOTE THAT OFF-DIAGONAL ELEMENTS WILL NEED TO BE HALVED LATER.
!
!     DIRFCK CALLS ONE OF THE FOLLOWING:
!         DIRFCK_RHF FOR RHF CALCULATIONS
!         DIRFCK_CPHF FOR CPHF CALCULATIONS
!         DIRFCK_UROHF FOR UHF OR ROHF CALCULATIONS
!         DIRFCK_GVB FOR GVB CALCULATIONS
!         DIRFCK_DYNDISP FOR EFMO GRADIENT CALCULATIONS
!     THE DFT CASES ARE HANDLED VIA HFSCAL AND CSCALT FACTORS.
!
!        SCFTYP    [IN] SCF TYPE 8 BYTE HOLLERITH CONSTANT.
!        IA(*)     [IN] TRIANGULAR INDEX ADDRESSES
!        DA(*)     [IN] RHF OR ROHF/UHF ALPHA DENSITY MATRIX.
!        FA(*)     [IN,OUT] RHF OR ROHF/UHF ALPHA FOCK MATRIX.
!        DB(*)     [IN] ROHF/UHF BETA DENSITY MATRIX.
!        FB(*)     [IN,OUT] ROHF/UHF BETA FOCK MATRIX.
!        GHONDO(*) [IN] HONDO INDEXED INTEGRALS.
!        L2        [IN] SIZE OF SYMMETRIC PACKED DENSITY OR FOCK MATRIX.
!        NINT      [IN,OUT] NUMBER OF NON-ZERO INTEGRALS.
!        NXYZ      [IN] NXYZ > 1 SELECTS CPHF CALCULATION.
!
      DYNDISP = NXYZ .lt. 0
      CPHF = NXYZ.GT.1.AND.SCFTYP.NE.DHUHF.AND.SCFTYP.NE.DHROHF
      UROHF = SCFTYP.EQ.DHUHF.OR.SCFTYP.EQ.DHROHF.AND..NOT.CPHF
      CPUHF= NXYZ.GT.1.AND.SCFTYP.EQ.DHUHF
      UROHF = UROHF.AND..NOT.(CPUHF.or. dyndisp)
      GVB = SCFTYP.EQ.DHGVB
      RHF = SCFTYP.EQ.DHRHF.AND..NOT.(CPHF.OR.UROHF.OR.GVB.or.dyndisp)

      HFSCAL=DFTTYP(3)
      CSCALT=1.0D+00
      IF(LCFLAG) THEN
         IF(LRINT) THEN
            HFSCAL=1.0D+00
            CSCALT=0.0D+00
         ELSE
            HFSCAL=0.0D+00
            CSCALT=1.0D+00
         ENDIF
      ENDIF
      IF(CAMFLAG.AND.LRINT) CSCALT=0.0D+00
!
!               FIRST THREE ARE SCF EQUATIONS, WITH THE FIRST TWO
!               POSSIBLY INVOLVING THE INTEGRALS PART OF DFT.
      IF (RHF) THEN
         CALL dirfck_rhf_gpu(IA,DA,FA,ISH,JSH,KSH,LSH, &
                        GHONDO,LSTRI,LSTRJ,LSTRK,LSTRL, &
                        HFSCAL,CSCALT,CUTOFFAO,NINT)
      ELSE IF (UROHF) THEN
         CALL DIRFCK_UROHF(IA,DA,FA,DB,FB,ISH,JSH,KSH,LSH, &
                          GHONDO,LSTRI,LSTRJ,LSTRK,LSTRL, &
                          HFSCAL,CSCALT,CUTOFFAO,NINT,NXYZ,L2)
      ELSE IF (GVB) THEN
         CALL DIRFCK_GVB(IA,DA,FA,L2,ISH,JSH,KSH,LSH, &
                        GHONDO,LSTRI,LSTRJ,LSTRK,LSTRL, &
                        CUTOFFAO,NINT)
!               NEXT CASE IS CLOSED SHELL RESPONSE EQN. USING AO INTS
      ELSE IF (CPHF) THEN
         CALL DIRFCK_CPHF(IA,DA,FA,NXYZ,ISH,JSH,KSH,LSH, &
                         GHONDO,LSTRI,LSTRJ,LSTRK,LSTRL, &
                         CUTOFFAO,NINT,HFSCAL,CSCALT)
      ELSE IF (CPUHF) THEN
         CALL DIRFCK_CPUHF(IA,DA,FA,NXYZ,ISH,JSH,KSH,LSH, &
                         GHONDO,LSTRI,LSTRJ,LSTRK,LSTRL, &
                         CUTOFFAO,NINT,DB,FB,HFSCAL,CSCALT,L2)
      ELSE IF (DYNDISP) THEN
         CALL DIRFCK_DYNDISP(IA,DA,FA,nxyz,ISH,JSH,KSH,LSH, &
                        GHONDO,LSTRI,LSTRJ,LSTRK,LSTRL, &
                        CUTOFFAO,NINT)
      ENDIF

      RETURN
      END



!*MODULE SCFLIB  *DECK DIRFCK_RHF
      SUBROUTINE dirfck_rhf_gpu(IA,DA,FA,ISH,JSH,KSH,LSH, &
                            GHONDO,ISTRIDE,JSTRIDE,KSTRIDE,LSTRIDE, &
                            HFSCAL,CSCALT,CUTOFF,NINT)
      use mx_limits, only: mxsh,mxgtot

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION IA(*),DA(*),FA(*),GHONDO(*)

      LOGICAL :: IANDJ,KANDL,SAME

      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT), &
                      CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT), &
                      KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH), &
                      KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL

      double precision,parameter :: HALF=0.5D00
      double precision,parameter :: FOUR=4.0D00
!
!     COMPUTES FOCK MATRIX ELEMENTS OF RHF OR RHF/DFT WAVEFUNCTION.
!     SHOULD NOT BE CALLED DIRECTLY, USE THE -DIRFCK- INTERFACE.
!        IA(*)                  [IN] TRIANGULAR INDEX ARRAY
!        DA(*)                  [IN] DENSITY MATRIX.
!        FA(*)                  [IN,OUT] FOCK MATRIX.
!        ISH,JSH,KSH,LSH        [IN] SHELL INDICES
!        GHONDO(*)              [IN] HONDO INDEXED INTEGRALS.
!        ISTRIDE,JSTRIDE,KSTRIDE,LSTRIDE -
!             [IN] HONDO INDEXED INTEGRALS SHELL STRIDES.
!        HFSCAL,CSCALT          [IN] DFT SCALING FACTORS.
!        CUTOFF                 [IN] NON-ZERO INTEGRAL CUTOFF
!        NINT                   [IN,OUT] NUMBER OF NONZERO INTEGRALS
!
      SAME  = ISH.EQ.KSH.AND.JSH.EQ.LSH
      IANDJ = ISH.EQ.JSH
      KANDL = KSH.EQ.LSH

      MINI = KMIN(ISH)
      MINJ = KMIN(JSH)
      MINK = KMIN(KSH)
      MINL = KMIN(LSH)
      MAXI = KMAX(ISH)
      MAXJ = KMAX(JSH)
      MAXK = KMAX(KSH)
      MAXL = KMAX(LSH)
      LOCI = KLOC(ISH)-MINI
      LOCJ = KLOC(JSH)-MINJ
      LOCK = KLOC(KSH)-MINK
      LOCL = KLOC(LSH)-MINL

      XVAL1 = HFSCAL
      XVAL4 = FOUR*CSCALT

      NIJ = 0
      MAXJ2 = MAXJ
      I_INDEX = 1
      DO I = MINI,MAXI
         IF (IANDJ) MAXJ2 = I

         I1 = I+LOCI

         IJ_INDEX = I_INDEX
         I_INDEX = I_INDEX + ISTRIDE

         DO J = MINJ,MAXJ2
            NIJ = NIJ+1
            MAXL2 = MAXL

            J1 = J+LOCJ
            I2 = I1
            J2 = J1
            IF (I1.LT.J1) THEN ! SORT <IJ|
               I2 = J1
               J2 = I1
            ENDIF

            IJK_INDEX = IJ_INDEX
            IJ_INDEX = IJ_INDEX + JSTRIDE

            NKL = NIJ

            DO K =  MINK,MAXK
               IF (KANDL) MAXL2 = K

               K1 = K + LOCK

               IF(SAME) THEN ! ACCOUNT FOR NON-UNIQUE PERMUTATIONS
                  ITMP = MIN(MAXL2-MINL+1,NKL)
                  IF (ITMP.EQ.0) CYCLE
                  MAXL2 = MINL + ITMP - 1
                  NKL = NKL - ITMP
               ENDIF

               IJKL_INDEX = IJK_INDEX
               IJK_INDEX = IJK_INDEX + KSTRIDE

               DO L=MINL,MAXL2

                  VAL = GHONDO( IJKL_INDEX )
                  IJKL_INDEX = IJKL_INDEX + LSTRIDE
                  !IF(ABS(VAL).LT.CUTOFF) CYCLE ! GOTO 300
                  NINT = NINT + 1

                  L1 = L + LOCL
                  K2 = K1
                  L2 = L1

                  IF (K2.LT.L2) THEN ! SORT |KL>
                     K2 = L1
                     L2 = K1
                  ENDIF

                  II = I2
                  JJ = J2
                  KK = K2
                  LL = L2

                  IF (II.LT.KK) THEN ! SORT <IJ|KL>
                     II = K2
                     JJ = L2
                     KK = I2
                     LL = J2
                  ELSE IF (II.EQ.KK.AND.JJ.LT.LL) THEN ! SORT <IJ|IL>
                     JJ = L2
                     LL = J2
                  ENDIF

                  II2 = IA(II)
                  JJ2 = IA(JJ)
                  KK2 = IA(KK)

                  IJ = II2 + JJ
                  IK = II2 + KK
                  IL = II2 + LL
                  JK = JJ2 + KK
                  JL = JJ2 + LL
                  KL = KK2 + LL
                  IF (JJ.LT.KK) JK = KK2 + JJ
                  IF (JJ.LT.LL) JL = IA(LL) + JJ
!
!       ACCOUNT FOR IDENTICAL PERMUTATIONS.
!
                  IF(II.EQ.JJ) VAL = VAL*HALF
                  IF(KK.EQ.LL) VAL = VAL*HALF
                  IF(II.EQ.KK.AND.JJ.EQ.LL) VAL = VAL*HALF
                  VAL1 = VAL*XVAL1
                  VAL4 = VAL*XVAL4
                  FA(IJ) = FA(IJ) + VAL4*DA(KL)
                  FA(KL) = FA(KL) + VAL4*DA(IJ)
                  FA(IK) = FA(IK) - VAL1*DA(JL)
                  FA(JL) = FA(JL) - VAL1*DA(IK)
                  FA(IL) = FA(IL) - VAL1*DA(JK)
                  FA(JK) = FA(JK) - VAL1*DA(IL)
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      END
      
      SUBROUTINE dirfck_gpu_rys(IA,DA,FA,GHONDO,NINT,&
                                kstart,katom,kng,kloc,kmin,kmax,&
                                ish,jsh,ksh,lsh,&
                                LSTRI,LSTRJ,LSTRK,LSTRL,dfttyp)

      USE mx_limits, ONLY: mxsh
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION IA(*),DA(*),FA(*),GHONDO(*)                   

      DIMENSION :: DFTTYP(20)
      LOGICAL :: IANDJ,KANDL,SAME

      double precision,parameter :: HALF=0.5D00
      double precision,parameter :: FOUR=4.0D00

      INTEGER,DIMENSION(mxsh) :: kstart,katom,kng,kloc,kmin,kmax
      INTEGER, INTENT(IN) :: ish,jsh,ksh,lsh
      INTEGER, INTENT(IN) :: LSTRI,LSTRJ,LSTRK,LSTRL
!!$omp declare target

      HFSCAL=DFTTYP(3)
      CSCALT=1.0D+00

         CALL dirfck_rhf_rys(IA,DA,FA,ISH,JSH,KSH,LSH, &
                        GHONDO,LSTRI,LSTRJ,LSTRK,LSTRL, &
                        HFSCAL,CSCALT,NINT,&
                        kstart,katom,kng,kloc,kmin,kmax)

      RETURN
      END
      SUBROUTINE dirfck_gpu_rys_2(IA,DA,FA,GHONDO,NINT,&
                                kstart,katom,kng,kloc,kmin,kmax,&
                                ish,jsh,ksh,lsh,&
                                LSTRI,LSTRJ,LSTRK,LSTRL,dfttyp)

      USE mx_limits, ONLY: mxsh
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION IA(*),DA(*),FA(*),GHONDO(*)                   

      DIMENSION :: DFTTYP(20)
      LOGICAL :: IANDJ,KANDL,SAME

      double precision,parameter :: HALF=0.5D00
      double precision,parameter :: FOUR=4.0D00

      INTEGER,DIMENSION(mxsh) :: kstart,katom,kng,kloc,kmin,kmax
      INTEGER, INTENT(IN) :: ish,jsh,ksh,lsh
      INTEGER, INTENT(IN) :: LSTRI,LSTRJ,LSTRK,LSTRL
!!$omp declare target

      HFSCAL=DFTTYP(3)
      CSCALT=1.0D+00

         CALL dirfck_rhf_rys_2(IA,DA,FA,ISH,JSH,KSH,LSH, &
                        GHONDO,LSTRI,LSTRJ,LSTRK,LSTRL, &
                        HFSCAL,CSCALT,NINT,&
                        kstart,katom,kng,kloc,kmin,kmax)

      RETURN
      END

      SUBROUTINE dirfck_rhf_rys(IA,DA,FA,ISH,JSH,KSH,LSH, &
                            GHONDO,ISTRIDE,JSTRIDE,KSTRIDE,LSTRIDE, &
                            HFSCAL,CSCALT,NINT,&
                            kstart,katom,kng,kloc,kmin,kmax)
      use mx_limits, only: mxsh

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION IA(*),DA(*),FA(*),GHONDO(*)

      LOGICAL :: IANDJ,KANDL,SAME

      double precision,parameter :: HALF=0.5D00
      double precision,parameter :: FOUR=4.0D00

      INTEGER,DIMENSION(mxsh) :: kstart,katom,kng,kloc,kmin,kmax
      INTEGER, INTENT(IN) :: ish,jsh,ksh,lsh
      INTEGER, INTENT(IN) :: ISTRIDE,JSTRIDE,KSTRIDE,LSTRIDE
!
!!$omp declare target
      SAME  = ISH.EQ.KSH.AND.JSH.EQ.LSH
      IANDJ = ISH.EQ.JSH
      KANDL = KSH.EQ.LSH

      MINI = KMIN(ISH)
      MINJ = KMIN(JSH)
      MINK = KMIN(KSH)
      MINL = KMIN(LSH)
      MAXI = KMAX(ISH)
      MAXJ = KMAX(JSH)
      MAXK = KMAX(KSH)
      MAXL = KMAX(LSH)
      LOCI = KLOC(ISH)-MINI
      LOCJ = KLOC(JSH)-MINJ
      LOCK = KLOC(KSH)-MINK
      LOCL = KLOC(LSH)-MINL

      XVAL1 = HFSCAL
      XVAL4 = FOUR*CSCALT

      NIJ = 0
      MAXJ2 = MAXJ
      I_INDEX = 1
      DO I = MINI,MAXI
         IF (IANDJ) MAXJ2 = I

         I1 = I+LOCI

         IJ_INDEX = I_INDEX
         I_INDEX = I_INDEX + ISTRIDE

         DO J = MINJ,MAXJ2
            NIJ = NIJ+1
            MAXL2 = MAXL

            J1 = J+LOCJ
            I2 = I1
            J2 = J1
            IF (I1.LT.J1) THEN ! SORT <IJ|
               I2 = J1
               J2 = I1
            ENDIF

            IJK_INDEX = IJ_INDEX
            IJ_INDEX = IJ_INDEX + JSTRIDE

            NKL = NIJ

            DO K =  MINK,MAXK
               IF (KANDL) MAXL2 = K

               K1 = K + LOCK

               IF(SAME) THEN ! ACCOUNT FOR NON-UNIQUE PERMUTATIONS
                  ITMP = MIN(MAXL2-MINL+1,NKL)
                  IF (ITMP.EQ.0) CYCLE
                  MAXL2 = MINL + ITMP - 1
                  NKL = NKL - ITMP
               ENDIF

               IJKL_INDEX = IJK_INDEX
               IJK_INDEX = IJK_INDEX + KSTRIDE

               DO L=MINL,MAXL2

                  VAL = GHONDO( IJKL_INDEX )
                  IJKL_INDEX = IJKL_INDEX + LSTRIDE
                  !IF(ABS(VAL).LT.CUTOFF) CYCLE ! GOTO 300
                  NINT = NINT + 1

                  L1 = L + LOCL
                  K2 = K1
                  L2 = L1

                  IF (K2.LT.L2) THEN ! SORT |KL>
                     K2 = L1
                     L2 = K1
                  ENDIF

                  II = I2
                  JJ = J2
                  KK = K2
                  LL = L2

                  IF (II.LT.KK) THEN ! SORT <IJ|KL>
                     II = K2
                     JJ = L2
                     KK = I2
                     LL = J2
                  ELSE IF (II.EQ.KK.AND.JJ.LT.LL) THEN ! SORT <IJ|IL>
                     JJ = L2
                     LL = J2
                  ENDIF

                  II2 = IA(II)
                  JJ2 = IA(JJ)
                  KK2 = IA(KK)

                  IJ = II2 + JJ
                  IK = II2 + KK
                  IL = II2 + LL
                  JK = JJ2 + KK
                  JL = JJ2 + LL
                  KL = KK2 + LL
                  IF (JJ.LT.KK) JK = KK2 + JJ
                  IF (JJ.LT.LL) JL = IA(LL) + JJ
!
!       ACCOUNT FOR IDENTICAL PERMUTATIONS.
!
                  IF(II.EQ.JJ) VAL = VAL*HALF
                  IF(KK.EQ.LL) VAL = VAL*HALF
                  IF(II.EQ.KK.AND.JJ.EQ.LL) VAL = VAL*HALF
                  VAL1 = VAL*XVAL1
                  VAL4 = VAL*XVAL4
                  FA(IJ) = FA(IJ) + VAL4*DA(KL)
                  FA(KL) = FA(KL) + VAL4*DA(IJ)
                  FA(IK) = FA(IK) - VAL1*DA(JL)
                  FA(JL) = FA(JL) - VAL1*DA(IK)
                  FA(IL) = FA(IL) - VAL1*DA(JK)
                  FA(JK) = FA(JK) - VAL1*DA(IL)
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      END

SUBROUTINE dirfck_rhf_rys_2(IA,DA,FA,ISH,JSH,KSH,LSH, &
                            GHONDO,ISTRIDE,JSTRIDE,KSTRIDE,LSTRIDE, &
                            HFSCAL,CSCALT,NINT,&
                            kstart,katom,kng,kloc,kmin,kmax)
      USE OMP_LIB 
      use mx_limits, only: mxsh

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

!$omp declare target
      DIMENSION IA(*),DA(*),FA(*),GHONDO(*)

      LOGICAL :: IANDJ,KANDL,SAME

      double precision,parameter :: HALF=0.5D00
      double precision,parameter :: FOUR=4.0D00

      INTEGER,DIMENSION(mxsh) :: kstart,katom,kng,kloc,kmin,kmax
      INTEGER, INTENT(IN) :: ish,jsh,ksh,lsh
      INTEGER, INTENT(IN) :: ISTRIDE,JSTRIDE,KSTRIDE,LSTRIDE
!!$omp declare target
      SAME  = ISH.EQ.KSH.AND.JSH.EQ.LSH
      IANDJ = ISH.EQ.JSH
      KANDL = KSH.EQ.LSH

      MINI = KMIN(ISH)
      MINJ = KMIN(JSH)
      MINK = KMIN(KSH)
      MINL = KMIN(LSH)
      MAXI = KMAX(ISH)
      MAXJ = KMAX(JSH)
      MAXK = KMAX(KSH)
      MAXL = KMAX(LSH)
      LOCI = KLOC(ISH)-MINI
      LOCJ = KLOC(JSH)-MINJ
      LOCK = KLOC(KSH)-MINK
      LOCL = KLOC(LSH)-MINL

      XVAL1 = HFSCAL
      XVAL4 = FOUR*CSCALT

      ! write(*,*) "ISTRIDE", ISTRIDE
      ! write(*,*) "JSTRIDE", JSTRIDE
      ! write(*,*) "KSTRIDE", KSTRIDE
      ! write(*,*) "LSTRIDE", LSTRIDE

      NIJ = 0
      MAXJ2 = MAXJ
      I_INDEX = 1
      DO I = MINI,MAXI
         IF (IANDJ) MAXJ2 = I

         I1 = I+LOCI

         IJ_INDEX = I_INDEX
         !write(*,*) "IJ_INDEX", IJ_INDEX
         I_INDEX = I_INDEX + ISTRIDE

         DO J = MINJ,MAXJ2
            NIJ = NIJ+1
            MAXL2 = MAXL

            J1 = J+LOCJ
            I2 = I1
            J2 = J1
            IF (I1.LT.J1) THEN ! SORT <IJ|
               I2 = J1
               J2 = I1
            ENDIF

            IJK_INDEX = IJ_INDEX
            IJ_INDEX = IJ_INDEX + JSTRIDE

            NKL = NIJ

            DO K =  MINK,MAXK
               IF (KANDL) MAXL2 = K

               K1 = K + LOCK

               IF(SAME) THEN ! ACCOUNT FOR NON-UNIQUE PERMUTATIONS
                  ITMP = MIN(MAXL2-MINL+1,NKL)
                  IF (ITMP.EQ.0) CYCLE
                  MAXL2 = MINL + ITMP - 1
                  NKL = NKL - ITMP
               ENDIF

               IJKL_INDEX = IJK_INDEX
               IJK_INDEX = IJK_INDEX + KSTRIDE

               DO L=MINL,MAXL2

                  VAL = GHONDO( IJKL_INDEX )
                  !IF(VAL .NE. 0.00) write(*,*) "val is", VAL
                  !write(*,*) "val is", VAL
                  IJKL_INDEX = IJKL_INDEX + LSTRIDE
                  !IF(ABS(VAL).LT.CUTOFF) CYCLE ! GOTO 300
                  NINT = NINT + 1

                  L1 = L + LOCL
                  K2 = K1
                  L2 = L1

                  IF (K2.LT.L2) THEN ! SORT |KL>
                     K2 = L1
                     L2 = K1
                  ENDIF

                  II = I2
                  JJ = J2
                  KK = K2
                  LL = L2

                  IF (II.LT.KK) THEN ! SORT <IJ|KL>
                     II = K2
                     JJ = L2
                     KK = I2
                     LL = J2
                  ELSE IF (II.EQ.KK.AND.JJ.LT.LL) THEN ! SORT <IJ|IL>
                     JJ = L2
                     LL = J2
                  ENDIF

                  II2 = IA(II)
                  JJ2 = IA(JJ)
                  KK2 = IA(KK)

                  IJ = II2 + JJ
                  IK = II2 + KK
                  IL = II2 + LL
                  JK = JJ2 + KK
                  JL = JJ2 + LL
                  KL = KK2 + LL
                  IF (JJ.LT.KK) JK = KK2 + JJ
                  IF (JJ.LT.LL) JL = IA(LL) + JJ
!
!       ACCOUNT FOR IDENTICAL PERMUTATIONS.
!
                  IF(II.EQ.JJ) VAL = VAL*HALF
                  IF(KK.EQ.LL) VAL = VAL*HALF
                  IF(II.EQ.KK.AND.JJ.EQ.LL) VAL = VAL*HALF
                  VAL1 = VAL*XVAL1
                  VAL4 = VAL*XVAL4
!$omp atomic
                  FA(IJ) = FA(IJ) + VAL4*DA(KL)
!$omp atomic
                  FA(KL) = FA(KL) + VAL4*DA(IJ)
!$omp atomic
                  FA(IK) = FA(IK) - VAL1*DA(JL)
!$omp atomic
                  FA(JL) = FA(JL) - VAL1*DA(IK)
!$omp atomic
                  FA(IL) = FA(IL) - VAL1*DA(JK)
!$omp atomic
                  FA(JK) = FA(JK) - VAL1*DA(IL)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      END

SUBROUTINE dirfck_rhf_rys_3(IA,DA,FA,ISH,JSH,KSH,LSH, &
                            GHONDO,ISTRIDE,JSTRIDE,KSTRIDE,LSTRIDE, &
                            HFSCAL,CSCALT,NINT,&
                            kstart,katom,kng,kloc,kmin,kmax)
      USE OMP_LIB 
      use mx_limits, only: mxsh

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

!$omp declare target
      DIMENSION IA(*),DA(*),FA(*),GHONDO(*)

      LOGICAL :: IANDJ,KANDL,SAME

      double precision,parameter :: HALF=0.5D00
      double precision,parameter :: FOUR=4.0D00

      INTEGER,DIMENSION(mxsh) :: kstart,katom,kng,kloc,kmin,kmax
      INTEGER, INTENT(IN) :: ish,jsh,ksh,lsh
      INTEGER, INTENT(IN) :: ISTRIDE,JSTRIDE,KSTRIDE,LSTRIDE
!!$omp declare target
      SAME  = ISH.EQ.KSH.AND.JSH.EQ.LSH
      IANDJ = ISH.EQ.JSH
      KANDL = KSH.EQ.LSH

      MINI = KMIN(ISH)
      MINJ = KMIN(JSH)
      MINK = KMIN(KSH)
      MINL = KMIN(LSH)
      MAXI = KMAX(ISH)
      MAXJ = KMAX(JSH)
      MAXK = KMAX(KSH)
      MAXL = KMAX(LSH)
      LOCI = KLOC(ISH)-MINI
      LOCJ = KLOC(JSH)-MINJ
      LOCK = KLOC(KSH)-MINK
      LOCL = KLOC(LSH)-MINL

      XVAL1 = HFSCAL
      XVAL4 = FOUR*CSCALT

      ! write(*,*) "ISTRIDE", ISTRIDE
      ! write(*,*) "JSTRIDE", JSTRIDE
      ! write(*,*) "KSTRIDE", KSTRIDE
      ! write(*,*) "LSTRIDE", LSTRIDE

      NIJ = 0
      MAXJ2 = MAXJ
      I_INDEX = 1
      DO I = MINI,MAXI
         IF (IANDJ) MAXJ2 = I

         I1 = I+LOCI

         IJ_INDEX = I_INDEX
         !write(*,*) "IJ_INDEX", IJ_INDEX
         I_INDEX = I_INDEX + ISTRIDE

         DO J = MINJ,MAXJ2
            NIJ = NIJ+1
            MAXL2 = MAXL

            J1 = J+LOCJ
            I2 = I1
            J2 = J1
            IF (I1.LT.J1) THEN ! SORT <IJ|
               I2 = J1
               J2 = I1
            ENDIF

            IJK_INDEX = IJ_INDEX
            IJ_INDEX = IJ_INDEX + JSTRIDE

            NKL = NIJ

            DO K =  MINK,MAXK
               IF (KANDL) MAXL2 = K

               K1 = K + LOCK

               IF(SAME) THEN ! ACCOUNT FOR NON-UNIQUE PERMUTATIONS
                  ITMP = MIN(MAXL2-MINL+1,NKL)
                  IF (ITMP.EQ.0) CYCLE
                  MAXL2 = MINL + ITMP - 1
                  NKL = NKL - ITMP
               ENDIF

               IJKL_INDEX = IJK_INDEX
               IJK_INDEX = IJK_INDEX + KSTRIDE

               DO L=MINL,MAXL2

                  VAL = GHONDO( IJKL_INDEX )
                  !IF(VAL .NE. 0.00) write(*,*) "val is", VAL
                  !write(*,*) "val is", VAL
                  IJKL_INDEX = IJKL_INDEX + LSTRIDE
                  !IF(ABS(VAL).LT.CUTOFF) CYCLE ! GOTO 300
                  NINT = NINT + 1

                  L1 = L + LOCL
                  K2 = K1
                  L2 = L1

                  IF (K2.LT.L2) THEN ! SORT |KL>
                     K2 = L1
                     L2 = K1
                  ENDIF

                  II = I2
                  JJ = J2
                  KK = K2
                  LL = L2

                  IF (II.LT.KK) THEN ! SORT <IJ|KL>
                     II = K2
                     JJ = L2
                     KK = I2
                     LL = J2
                  ELSE IF (II.EQ.KK.AND.JJ.LT.LL) THEN ! SORT <IJ|IL>
                     JJ = L2
                     LL = J2
                  ENDIF

                  II2 = IA(II)
                  JJ2 = IA(JJ)
                  KK2 = IA(KK)

                  IJ = II2 + JJ
                  IK = II2 + KK
                  IL = II2 + LL
                  JK = JJ2 + KK
                  JL = JJ2 + LL
                  KL = KK2 + LL
                  IF (JJ.LT.KK) JK = KK2 + JJ
                  IF (JJ.LT.LL) JL = IA(LL) + JJ
!
!       ACCOUNT FOR IDENTICAL PERMUTATIONS.
!
                  IF(II.EQ.JJ) VAL = VAL*HALF
                  IF(KK.EQ.LL) VAL = VAL*HALF
                  IF(II.EQ.KK.AND.JJ.EQ.LL) VAL = VAL*HALF
                  VAL1 = VAL*XVAL1
                  VAL4 = VAL*XVAL4
!$omp atomic
                  FA(IJ) = FA(IJ) + VAL4*DA(KL)
                  !write(*,*) "FA(IJ) is", FA(IJ)
                  !write(*,*) "VAL4", VAL4
                  !write(*,*) "DA(KL)", DA(KL)
                  !write(*,*) "second term", VAL4*DA(KL)
!$omp atomic
                  FA(KL) = FA(KL) + VAL4*DA(IJ)
!$omp atomic
                  FA(IK) = FA(IK) - VAL1*DA(JL)
!$omp atomic
                  FA(JL) = FA(JL) - VAL1*DA(IK)
!$omp atomic
                  FA(IL) = FA(IL) - VAL1*DA(JK)
!$omp atomic
                  FA(JK) = FA(JK) - VAL1*DA(IL)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      END


      SUBROUTINE SHELLS_gpu2_ij(NELEC,ISH,JSH,FLIP,&
      IJGT,IJX,IJY,IJZ,&
      IANDJ,SAME,&
      INU,JNU,NGTI,NGTJ,&
      LIT,LJT,LOCI,LOCJ,&
      mini,maxi,minj,maxj,&
      NIJ,IJ,KL,IJKL,&
      NGTH,C,&
      kstart,katom,ktype,kng,kloc,kmin,kmax,&
      ex,cs,cp,cd,&
      ga,csa,cpa,cda,&
      gb,csb,cpb,cdb,&
      AX,AY,AZ,BX,BY,BZ,RAB,&
      NGA,NGB)
      use omp_lib
      use mx_limits, only: mxsh,mxgsh,mxgtot,mxatm

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      LOGICAL FLIP
      LOGICAL IANDJ,KANDL,SAME

      INTEGER :: NGTI,NGTJ
      double precision :: C(3,MXATM)

      DIMENSION :: EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT)
      INTEGER,DIMENSION(mxsh) :: kstart,katom,kng,kloc,kmin,kmax,ktype
      INTEGER,DIMENSION(4) :: NGTH
      INTEGER :: NELEC 
      ! DIMENSION :: NGTH(4)

      DIMENSION :: GA(MXGSH),CSA(MXGSH),CPA(MXGSH),CDA(MXGSH)
      DIMENSION :: GB(MXGSH),CSB(MXGSH),CPB(MXGSH),CDB(MXGSH)
      double precision :: AX,AY,AZ,BX,BY,BZ,RAB
      integer :: NGA,NGB

      double precision :: qq4
      INTEGER :: LIT,LJT,LOCI,LOCJ
      INTEGER :: MINI,MINJ,MAXI,MAXJ
      INTEGER :: NIJ,IJ,KL,IJKL
      !ADDED
      INTEGER :: IJN,JMAX,I,J,N1,LMAX,KLN,K,L,NN
      INTEGER :: INU,JNU,ISH,JSH
      INTEGER, DIMENSION(84) :: IX,IY,IZ,JX,JY,JZ,IJGT,IJX,IJY,IJZ
      double precision :: NX,NY,NZ
      INTEGER:: I1,J1,I2,J2  
      INTEGER :: NROOTS
      DATA JX /   0,  49,   0,   0,  98,   0,   0,  49,  49,   0, &
               147,   0,   0,  98,  98,  49,   0,  49,   0,  49, &
               196,   0,   0, 147, 147,  49,   0,  49,   0,  98, &
                98,   0,  98,  49,  49, &
               245,   0,   0, 196, 196,  49,   0,  49,   0, 147, &
               147,  98,   0,  98,   0, 147,  49,  49,  98,  98, &
                49, &
               294,   0,   0, 245, 245,  49,   0,  49,   0, 196, &
               196,  98,   0,  98,   0, 196,  49,  49, 147, 147, &
                 0, 147, 147,  98,  49,  98,  49,  98/
      DATA IX /   1, 344,   1,   1, 687,   1,   1, 344, 344,   1, &
              1030,   1,   1, 687, 687, 344,   1, 344,   1, 344, &
              1373,   1,   1,1030,1030, 344,   1, 344,   1, 687, &
               687,   1, 687, 344, 344, &
              1716,   1,   1,1373,1373, 344,   1, 344,   1,1030, &
              1030, 687,   1, 687,   1,1030, 344, 344, 687, 687, &
               344, &
              2059,   1,   1,1716,1716, 344,   1, 344,   1,1373, &
              1373, 687,   1, 687,   1,1373, 344, 344,1030,1030, &
                 1,1030,1030, 687, 344, 687, 344, 687/
      DATA JY /   0,   0,  49,   0,   0,  98,   0,  49,   0,  49, &
                 0, 147,   0,  49,   0,  98,  98,   0,  49,  49, &
                 0, 196,   0,  49,   0, 147, 147,   0,  49,  98, &
                 0,  98,  49,  98,  49, &
                 0, 245,   0,  49,   0, 196, 196,   0,  49,  98, &
                 0, 147, 147,   0,  98,  49, 147,  49,  98,  49, &
                98, &
                 0, 294,   0,  49,   0, 245, 245,   0,  49,  98, &
                 0, 196, 196,   0,  98,  49, 196,  49, 147,   0, &
               147,  98,  49, 147, 147,  49,  98,  98/
      DATA IY /   1,   1, 344,   1,   1, 687,   1, 344,   1, 344, &
                 1,1030,   1, 344,   1, 687, 687,   1, 344, 344, &
                 1,1373,   1, 344,   1,1030,1030,   1, 344, 687, &
                 1, 687, 344, 687, 344, &
                 1,1716,   1, 344,   1,1373,1373,   1, 344, 687, &
                 1,1030,1030,   1, 687, 344,1030, 344, 687, 344, &
               687, &
                 1,2059,   1, 344,   1,1716,1716,   1, 344, 687, &
                 1,1373,1373,   1, 687, 344,1373, 344,1030,   1, &
              1030, 687, 344,1030,1030, 344, 687, 687/
      DATA JZ /   0,   0,   0,  49,   0,   0,  98,   0,  49,  49, &
                 0,   0, 147,   0,  49,   0,  49,  98,  98,  49, &
                 0,   0, 196,   0,  49,   0,  49, 147, 147,   0, &
                98,  98,  49,  49,  98, &
                 0,   0, 245,   0,  49,   0,  49, 196, 196,   0, &
                98,   0,  98, 147, 147,  49,  49, 147,  49,  98, &
                98, &
                 0,   0, 294,   0,  49,   0,  49, 245, 245,   0, &
                98,   0,  98, 196, 196,  49,  49, 196,   0, 147, &
               147,  49,  98,  49,  98, 147, 147,  98/
      DATA IZ /   1,   1,   1, 344,   1,   1, 687,   1, 344, 344, &
                 1,   1,1030,   1, 344,   1, 344, 687, 687, 344, &
                 1,   1,1373,   1, 344,   1, 344,1030,1030,   1, &
               687, 687, 344, 344, 687, &
                 1,   1,1716,   1, 344,   1, 344,1373,1373,   1, &
               687,   1, 687,1030,1030, 344, 344,1030, 344, 687, &
               687, &
                 1,   1,2059,   1, 344,   1, 344,1716,1716,   1, &
               687,   1, 687,1373,1373, 344, 344,1373,   1,1030, &
              1030, 344, 687, 344, 687,1030,1030, 687/

!$omp declare target
      IANDJ = ISH .EQ. JSH
      IF (KTYPE(ISH) .LT. KTYPE(JSH)) THEN
         INU = JSH
         JNU = ISH
         NGTI = NGTH(2)
         NGTJ = NGTH(1)
      ELSE
         INU = ISH
         JNU = JSH
         NGTI = NGTH(1)
         NGTJ = NGTH(2)
      END IF

      I = KATOM(INU)
      !write(*,*) "I is ", I
      AX = C(1,I)
      AY = C(2,I)
      AZ = C(3,I)
      I1 = KSTART(INU)
      I2 = I1+KNG(INU)-1
      LIT = KTYPE(INU)
      MINI = KMIN(INU)
      MAXI = KMAX(INU)
      LOCI = KLOC(INU)-MINI
      NGA = 0
      !write(*,*) "I1 is ", I1
      !write(*,*) "I2 is ", I2
      DO I = I1,I2
         NGA = NGA+1
         GA(NGA) = EX(I)
         CSA(NGA) = CS(I)
         CPA(NGA) = CP(I)
         CDA(NGA) = CD(I)
      ENDDO

      J = KATOM(JNU)
      BX = C(1,J)
      BY = C(2,J)
      BZ = C(3,J)
      J1 = KSTART(JNU)
      J2 = J1+KNG(JNU)-1
      LJT = KTYPE(JNU)
      MINJ = KMIN(JNU)
      MAXJ = KMAX(JNU)
      LOCJ = KLOC(JNU)-MINJ
      NGB = 0
       DO J = J1,J2
         NGB = NGB+1
         GB(NGB) = EX(J)
         CSB(NGB) = CS(J)
         CPB(NGB) = CP(J)
         CDB(NGB) = CD(J)
       ENDDO
      RAB = ((AX-BX)*(AX-BX) + (AY-BY)*(AY-BY) + (AZ-BZ)*(AZ-BZ))
      !write(*,*) "RAB", RAB

      IJ = 0
      JMAX = MAXJ
      DO I = MINI,MAXI
         NX = IX(I)
         NY = IY(I)
         NZ = IZ(I)
         IF (IANDJ) JMAX = I
         DO J = MINJ,JMAX
            IJ = IJ+1
            IJX(IJ) = NX+JX(J)
            !write(*,*) "IJX(IJ)", IJX(IJ)
            IJY(IJ) = NY+JY(J)
            IJZ(IJ) = NZ+JZ(J)
            IJGT(IJ) = NGTI*(I-MINI)+NGTJ*(J-MINJ)+1
         ENDDO
      ENDDO
      RETURN
      END
      SUBROUTINE IJPRIM_gpu_2(DDIJ,A,R,X1,Y1,Z1,IJD,IANDJ,KANDL,SAME,&
      LIT,mini,maxi,minj,maxj,NIJ,&
      ag,csa,cpa,cda,&
      bg,csb,cpb,cdb,&
      XI,YI,ZI,XJ,YJ,ZJ,RRI,&
      NGA,NGB)

      use mx_limits, only: mxgsh,mxg2

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      LOGICAL IANDJ,KANDL,SAME,OUT,NORM
      DIMENSION DDIJ(49*MXG2)
      DIMENSION :: A(MXG2),R(MXG2),X1(MXG2),Y1(MXG2),Z1(MXG2),IJD(10)
      DIMENSION :: AG(MXGSH),CSA(MXGSH),CPA(MXGSH),CDA(MXGSH)
      DIMENSION :: BG(MXGSH),CSB(MXGSH),CPB(MXGSH),CDB(MXGSH)

      double precision :: XI,YI,ZI,XJ,YJ,ZJ,RRI
      integer :: NGA,NGB
      double precision :: qq4
      INTEGER :: LIT
      INTEGER :: MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL
      INTEGER :: NIJ,IJ,KL,IJKL

      PARAMETER (SQRT3=1.73205080756888D+00)
      PARAMETER (SQRT5=2.23606797749979D+00)
      PARAMETER (SQRT7=2.64575131106459D+00)
      PARAMETER (SQRT9=3.0D+00)
      PARAMETER (SQRT11=3.3166247903553998D+00)
      PARAMETER (ZERO=0.0D+00)
      PARAMETER (ONE=1.0D+00)
!$omp declare target

      NORM = NORMF .NE. 1 .OR. NORMP .NE. 1
      MAX = MAXJ
      N = 0
      NN = 0
      NM = -2**20
      DO 180 I = MINI,MAXI
         GO TO (100,100,120,120,100,120,120,100,120,120,&
               100,120,120,100,120,120,120,120,120,100,&
               100,120,120,100,120,120,120,120,120,100,&
               120,120,100,120,120,&
               100,120,120,100,120,120,120,120,120,100,&
               120,120,120,120,120,100,120,120,100,120,&
               120,&
               100,120,120,100,120,120,120,120,120,100,&
               120,120,120,120,120,100,120,120,100,120,&
               120,100,120,120,120,120,120,100),I
  100    NM = NN
  120    NN = NM
         IF (IANDJ) MAX = I
         DO 170 J = MINJ,MAX
            GO TO (140,140,160,160,140,160,160,140,160,160,&
                  140,160,160,140,160,160,160,160,160,140,&
                  140,160,160,140,160,160,160,160,160,140,&
                  160,160,140,160,160,&
                  140,160,160,140,160,160,160,160,160,140,&
                  160,160,160,160,160,140,160,160,140,160,&
                  160,&
                  140,160,160,140,160,160,160,160,160,140,&
                  160,160,160,160,160,140,160,160,140,160,&
                  160,140,160,160,160,160,160,140),J
  140       NN = NN+1
  160       N = N+1
            IJD(N) = NN
  170    CONTINUE
  180 CONTINUE

      NIJ = 0
      JBMAX = NGB
      DO IA = 1,NGA
         AI = AG(IA)
         ARRI = AI*RRI
         AXI = AI*XI
         AYI = AI*YI
         AZI = AI*ZI
         CSI = CSA(IA)
         CPI = CPA(IA)
         CDI = CDA(IA)

         IF (IANDJ) JBMAX = IA
         DO 520 JB = 1,JBMAX
            AJ = BG(JB)
            AA = AI+AJ
            AAINV = ONE/AA
            DUM = AJ*ARRI*AAINV
            !IF (DUM .GT. TOL) GO TO 520
            CSJ = CSB(JB)
            CPJ = CPB(JB)
            CDJ = CDB(JB)

            NM = 49*NIJ
            NN = NM
            NIJ = NIJ+1
            R(NIJ) = DUM
            A(NIJ) = AA
            X1(NIJ) = (AXI+AJ*XJ)*AAINV
            Y1(NIJ) = (AYI+AJ*YJ)*AAINV
            Z1(NIJ) = (AZI+AJ*ZJ)*AAINV

            DUM1 = ZERO
            DUM2 = ZERO
            DO 420 I = MINI,MAXI

               GO TO (200,220,420,420,240,420,420,260,420,420),I
  200          DUM1 = CSI*AAINV
               GO TO 280
  220          DUM1 = CPI*AAINV
               GO TO 280
  240          DUM1 = CDI*AAINV
               GO TO 280
  260          DUM1 = DUM1*SQRT3
               GO TO 280


  280          IF (IANDJ) MAX = I
                 DO 400 J = MINJ,MAX
                  GO TO (300,320,400,400,340,400,400,360,400,400),J
  300             DUM2 = DUM1*CSJ
                  GO TO 380
  320             DUM2 = DUM1*CPJ
                  GO TO 380
  340             DUM2 = DUM1*CDJ
                  GO TO 380
  360             DUM2 = DUM2*SQRT3
                  GO TO 380


  380             NN = NN+1
                  DDIJ(NN) = DUM2
!                  write(*,*) "DDIJ is", DDIJ(NN)
  400          CONTINUE
   420       CONTINUE
            IF ( .NOT. IANDJ) GO TO 520
            IF (IA .EQ. JB) GO TO 520
            !GO TO (500,440,460,455,450,445,444),LIT
            GO TO (500,440),LIT
  440       IF (MINI .EQ. 2) GO TO 500
            DDIJ(NM+2) = DDIJ(NM+2)+CSI*CPJ*AAINV
            GO TO 480
  ! 444       DDIJ(NM+28) = DDIJ(NM+28)+DDIJ(NM+28)
  !           DDIJ(NM+27) = DDIJ(NM+27)+DDIJ(NM+27)
  !           DDIJ(NM+26) = DDIJ(NM+26)+DDIJ(NM+26)
  !           DDIJ(NM+25) = DDIJ(NM+25)+DDIJ(NM+25)
  !           DDIJ(NM+24) = DDIJ(NM+24)+DDIJ(NM+24)
  !           DDIJ(NM+23) = DDIJ(NM+23)+DDIJ(NM+23)
  !           DDIJ(NM+22) = DDIJ(NM+22)+DDIJ(NM+22)
  !           DDIJ(NM+21) = DDIJ(NM+21)+DDIJ(NM+21)
  !           DDIJ(NM+20) = DDIJ(NM+20)+DDIJ(NM+20)
  !           DDIJ(NM+19) = DDIJ(NM+19)+DDIJ(NM+19)
  !           DDIJ(NM+18) = DDIJ(NM+18)+DDIJ(NM+18)
  !           DDIJ(NM+17) = DDIJ(NM+17)+DDIJ(NM+17)
  !           DDIJ(NM+16) = DDIJ(NM+16)+DDIJ(NM+16)
  ! 445       DDIJ(NM+15) = DDIJ(NM+15)+DDIJ(NM+15)
  !           DDIJ(NM+14) = DDIJ(NM+14)+DDIJ(NM+14)
  !           DDIJ(NM+13) = DDIJ(NM+13)+DDIJ(NM+13)
  !           DDIJ(NM+12) = DDIJ(NM+12)+DDIJ(NM+12)
  !           DDIJ(NM+11) = DDIJ(NM+11)+DDIJ(NM+11)
  ! 450       DDIJ(NM+10) = DDIJ(NM+10)+DDIJ(NM+10)
  !           DDIJ(NM+9) = DDIJ(NM+9)+DDIJ(NM+9)
  !           DDIJ(NM+8) = DDIJ(NM+8)+DDIJ(NM+8)
  !           DDIJ(NM+7) = DDIJ(NM+7)+DDIJ(NM+7)
  ! 455       DDIJ(NM+6) = DDIJ(NM+6)+DDIJ(NM+6)
  !           DDIJ(NM+5) = DDIJ(NM+5)+DDIJ(NM+5)
  !           DDIJ(NM+4) = DDIJ(NM+4)+DDIJ(NM+4)
  ! 460       DDIJ(NM+2) = DDIJ(NM+2)+DDIJ(NM+2)
  480       DDIJ(NM+3) = DDIJ(NM+3)+DDIJ(NM+3)
  500       DDIJ(NM+1) = DDIJ(NM+1)+DDIJ(NM+1)
  520 CONTINUE
      ENDDO
      RETURN
      END
      SUBROUTINE IJPRIM_gpu_2222(DDIJ,A,R,X1,Y1,Z1,IJD,IANDJ,KANDL,SAME,&
      LIT,mini,maxi,minj,maxj,NIJ,&
      ag,csa,cpa,cda,&
      bg,csb,cpb,cdb,&
      XI,YI,ZI,XJ,YJ,ZJ,RRI)

      use mx_limits, only: mxgsh,mxg2

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      LOGICAL IANDJ,KANDL,SAME,OUT,NORM
      DIMENSION DDIJ(49*MXG2)
      DIMENSION :: A(MXG2),R(MXG2),X1(MXG2),Y1(MXG2),Z1(MXG2),IJD(10)
      DIMENSION :: AG(MXGSH),CSA(MXGSH),CPA(MXGSH),CDA(MXGSH)
      DIMENSION :: BG(MXGSH),CSB(MXGSH),CPB(MXGSH),CDB(MXGSH)

      double precision :: XI,YI,ZI,XJ,YJ,ZJ,RRI
      !integer :: NGA,NGB
      double precision :: qq4
      INTEGER :: LIT
      INTEGER :: MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL
      INTEGER :: NIJ,IJ,KL,IJKL

      PARAMETER (SQRT3=1.73205080756888D+00)
      PARAMETER (SQRT5=2.23606797749979D+00)
      PARAMETER (SQRT7=2.64575131106459D+00)
      PARAMETER (SQRT9=3.0D+00)
      PARAMETER (SQRT11=3.3166247903553998D+00)
      PARAMETER (ZERO=0.0D+00)
      PARAMETER (ONE=1.0D+00)
!$omp declare target

      NORM = NORMF .NE. 1 .OR. NORMP .NE. 1
      MAX = MAXJ
      N = 0
      NN = 0
      NM = -2**20
      DO 180 I = MINI,MAXI
         GO TO (100,100,120,120,100,120,120,100,120,120,&
               100,120,120,100,120,120,120,120,120,100,&
               100,120,120,100,120,120,120,120,120,100,&
               120,120,100,120,120,&
               100,120,120,100,120,120,120,120,120,100,&
               120,120,120,120,120,100,120,120,100,120,&
               120,&
               100,120,120,100,120,120,120,120,120,100,&
               120,120,120,120,120,100,120,120,100,120,&
               120,100,120,120,120,120,120,100),I
  100    NM = NN
  120    NN = NM
         IF (IANDJ) MAX = I
         DO 170 J = MINJ,MAX
            GO TO (140,140,160,160,140,160,160,140,160,160,&
                  140,160,160,140,160,160,160,160,160,140,&
                  140,160,160,140,160,160,160,160,160,140,&
                  160,160,140,160,160,&
                  140,160,160,140,160,160,160,160,160,140,&
                  160,160,160,160,160,140,160,160,140,160,&
                  160,&
                  140,160,160,140,160,160,160,160,160,140,&
                  160,160,160,160,160,140,160,160,140,160,&
                  160,140,160,160,160,160,160,140),J
  140       NN = NN+1
  160       N = N+1
            IJD(N) = NN
  170    CONTINUE
  180 CONTINUE

      NIJ = 0

         AI = AG(IA)
         ARRI = AI*RRI
         AXI = AI*XI
         AYI = AI*YI
         AZI = AI*ZI
         CSI = CSA(IA)
         CPI = CPA(IA)
         CDI = CDA(IA)

         IF (IANDJ) JBMAX = IA
         !DO 520 JB = 1,JBMAX
            AJ = BG(JB)
            AA = AI+AJ
            AAINV = ONE/AA
            DUM = AJ*ARRI*AAINV
            !IF (DUM .GT. TOL) GO TO 520
            CSJ = CSB(JB)
            CPJ = CPB(JB)
            CDJ = CDB(JB)

            NM = 49*NIJ
            !write(*,*) "NM",NM
            NN = NM
            !write(*,*) "NN",NN
            NIJ = NIJ+1
            !write(*,*) "NIJ",NIJ
            R(NIJ) = DUM
            A(NIJ) = AA
            X1(NIJ) = (AXI+AJ*XJ)*AAINV
            Y1(NIJ) = (AYI+AJ*YJ)*AAINV
            Z1(NIJ) = (AZI+AJ*ZJ)*AAINV

            DUM1 = ZERO
            DUM2 = ZERO
            DO 420 I = MINI,MAXI

               GO TO (200,220,420,420,240,420,420,260,420,420),I
  200          DUM1 = CSI*AAINV
               GO TO 280
  220          DUM1 = CPI*AAINV
               GO TO 280
  240          DUM1 = CDI*AAINV
               GO TO 280
  260          DUM1 = DUM1*SQRT3
               GO TO 280


  280          IF (IANDJ) MAX = I
                 DO 400 J = MINJ,MAX
                  GO TO (300,320,400,400,340,400,400,360,400,400),J
  300             DUM2 = DUM1*CSJ
                  GO TO 380
  320             DUM2 = DUM1*CPJ
                  GO TO 380
  340             DUM2 = DUM1*CDJ
                  GO TO 380
  360             DUM2 = DUM2*SQRT3
                  GO TO 380


  380             NN = NN+1
                  DDIJ(NN) = DUM2
                 !write(*,*) "DDIJ is", DDIJ(NN)
  400          CONTINUE
420          CONTINUE
      RETURN
      END


      SUBROUTINE SHELLS_gpu2_kl(NELEC,ISH,JSH,KSH,LSH,FLIP,&
      IK,KLGT,KLX,KLY,KLZ,&
      IANDJ,KANDL,SAME,&
      KNU,LNU,NGTK,NGTL,&
      LKT,LLT,LOCK,LOCL,&
      mink,maxk,minl,maxl,&
      IJ,KL,IJKL,&
      NGTH,C,&
      kstart,katom,ktype,kng,kloc,kmin,kmax,&
      ex,cs,cp,cd,&
      gc,csc,cpc,cdc,&
      gd,csd,cpd,cdd,&
      CX,CY,CZ,DX,DY,DZ,RCD,&
      NGC,NGD)
      use mx_limits, only: mxsh,mxgsh,mxgtot,mxatm
      use omp_lib

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      LOGICAL FLIP
      LOGICAL IANDJ,KANDL,SAME

      DIMENSION :: KX(84),KY(84),KZ(84),LX(84),LY(84),LZ(84)

      INTEGER :: INU,JNU,KNU,LNU,NGTI,NGTJ,NGTK,NGTL
      DIMENSION :: IJGT(784),IJX(784),IJY(784),IJZ(784),IK(784)
      DIMENSION :: KLGT(784),KLX(784),KLY(784),KLZ(784)
      double precision :: C(3,MXATM)

      DIMENSION :: EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT)
      DIMENSION :: KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH)
      DIMENSION :: KLOC(MXSH),KMIN(MXSH),KMAX(MXSH)

      DIMENSION :: NGTH(4)

      DIMENSION :: GC(MXGSH),CSC(MXGSH),CPC(MXGSH),CDC(MXGSH)
      DIMENSION :: GD(MXGSH),CSD(MXGSH),CPD(MXGSH),CDD(MXGSH)
      double precision :: CX,CY,CZ,DX,DY,DZ,RCD
      integer :: NGC,NGD

      double precision :: qq4
      INTEGER :: LOCI,LOCJ,LOCK,LOCL
      INTEGER :: MINK,MINL,MAXK,MAXL
      INTEGER :: NIJ,IJ,KL,IJKL
      DATA LX /   0,   1,   0,   0,   2,   0,   0,   1,   1,   0,&
                 3,   0,   0,   2,   2,   1,   0,   1,   0,   1,&
                 4,   0,   0,   3,   3,   1,   0,   1,   0,   2,&
                 2,   0,   2,   1,   1,&
                 5,   0,   0,   4,   4,   1,   0,   1,   0,   3,&
                 3,   2,   0,   2,   0,   3,   1,   1,   2,   2,&
                 1,&
                 6,   0,   0,   5,   5,   1,   0,   1,   0,   4,&
                 4,   2,   0,   2,   0,   4,   1,   1,   3,   3,&
                 0,   3,   3,   2,   1,   2,   1,   2/
      DATA KX /   0,   7,   0,   0,  14,   0,   0,   7,   7,   0,&
                21,   0,   0,  14,  14,   7,   0,   7,   0,   7,&
                28,   0,   0,  21,  21,   7,   0,   7,   0,  14,&
                14,   0,  14,   7,   7,&
                35,   0,   0,  28,  28,   7,   0,   7,   0,  21,&
                21,  14,   0,  14,   0,  21,   7,   7,  14,  14,&
                 7,&
                42,   0,   0,  35,  35,   7,   0,   7,   0,  28,&
                28,  14,   0,  14,   0,  28,   7,   7,  21,  21,&
                 0,  21,  21,  14,   7,  14,   7,  14/
     DATA LY /   0,   0,   1,   0,   0,   2,   0,   1,   0,   1,&
                 0,   3,   0,   1,   0,   2,   2,   0,   1,   1,&
                 0,   4,   0,   1,   0,   3,   3,   0,   1,   2,&
                 0,   2,   1,   2,   1,&
                 0,   5,   0,   1,   0,   4,   4,   0,   1,   2,&
                 0,   3,   3,   0,   2,   1,   3,   1,   2,   1,&
                 2,&
                 0,   6,   0,   1,   0,   5,   5,   0,   1,   2,&
                 0,   4,   4,   0,   2,   1,   4,   1,   3,   0,&
                 3,   2,   1,   3,   3,   1,   2,   2/
      DATA KY /   0,   0,   7,   0,   0,  14,   0,   7,   0,   7,&
                 0,  21,   0,   7,   0,  14,  14,   0,   7,   7,&
                 0,  28,   0,   7,   0,  21,  21,   0,   7,  14,&
                 0,  14,   7,  14,   7,&
                 0,  35,   0,   7,   0,  28,  28,   0,   7,  14,&
                 0,  21,  21,   0,  14,   7,  21,   7,  14,   7,&
                14,&
                 0,  42,   0,   7,   0,  35,  35,   0,   7,  14,&
                 0,  28,  28,   0,  14,   7,  28,   7,  21,   0,&
                21,  14,   7,  21,  21,   7,  14,  14/
      DATA LZ /   0,   0,   0,   1,   0,   0,   2,   0,   1,   1,&
                 0,   0,   3,   0,   1,   0,   1,   2,   2,   1,&
                 0,   0,   4,   0,   1,   0,   1,   3,   3,   0,&
                 2,   2,   1,   1,   2,&
                 0,   0,   5,   0,   1,   0,   1,   4,   4,   0,&
                 2,   0,   2,   3,   3,   1,   1,   3,   1,   2,&
                 2,&
                 0,   0,   6,   0,   1,   0,   1,   5,   5,   0,&
                 2,   0,   2,   4,   4,   1,   1,   4,   0,   3,&
                 3,   1,   2,   1,   2,   3,   3,   2/
      DATA KZ /   0,   0,   0,   7,   0,   0,  14,   0,   7,   7,&
                 0,   0,  21,   0,   7,   0,   7,  14,  14,   7,&
                 0,   0,  28,   0,   7,   0,   7,  21,  21,   0,&
                14,  14,   7,   7,  14,&
                 0,   0,  35,   0,   7,   0,   7,  28,  28,   0,&
                14,   0,  14,  21,  21,   7,   7,  21,   7,  14,&
                14,&
                 0,   0,  42,   0,   7,   0,   7,  35,  35,   0,&
                14,   0,  14,  28,  28,   7,   7,  28,   0,  21,&
                21,   7,  14,   7,  14,  21,  21,  14/

!$omp declare target
      KANDL = KSH .EQ. LSH
      SAME = ISH .EQ. KSH .AND. JSH .EQ. LSH
      !write(*,*) "here"
      IF (KTYPE(KSH) .LT. KTYPE(LSH)  .AND.  FLIP) THEN
         KNU = LSH
         LNU = KSH
         NGTK = NGTH(4)
         NGTL = NGTH(3)
      ELSE
         KNU = KSH
         LNU = LSH
         NGTK = NGTH(3)
         NGTL = NGTH(4)
      END IF

      K = KATOM(KNU)
      CX = C(1,K)
      CY = C(2,K)
      CZ = C(3,K)
      K1 = KSTART(KNU)
      K2 = K1+KNG(KNU)-1
      LKT = KTYPE(KNU)
      MINK = KMIN(KNU)
      MAXK = KMAX(KNU)
      LOCK = KLOC(KNU)-MINK
      NGC = 0
      DO K = K1,K2
         NGC = NGC+1
         GC(NGC) = EX(K)
         CSC(NGC) = CS(K)
         CPC(NGC) = CP(K)
         CDC(NGC) = CD(K)
      ENDDO
      L = KATOM(LNU)
      DX = C(1,L)
      DY = C(2,L)
      DZ = C(3,L)
      L1 = KSTART(LNU)
      L2 = L1+KNG(LNU)-1
      LLT = KTYPE(LNU)
      MINL = KMIN(LNU)
      MAXL = KMAX(LNU)
      LOCL = KLOC(LNU)-MINL
      NGD = 0
      DO L = L1,L2
         NGD = NGD+1
         GD(NGD) = EX(L)
         CSD(NGD) = CS(L)
         CPD(NGD) = CP(L)
         CDD(NGD) = CD(L)
      ENDDO
      RCD = ((CX-DX)*(CX-DX) + (CY-DY)*(CY-DY) + (CZ-DZ)*(CZ-DZ))
      !write(*,*) "RCD", RCD

      KL = 0
      LMAX = MAXL
      DO K = MINK,MAXK
         NX = KX(K)
         NY = KY(K)
         NZ = KZ(K)
         IF (KANDL) LMAX = K
         DO L = MINL,LMAX
            KL = KL+1
            KLX(KL) = NX+LX(L)
            KLY(KL) = NY+LY(L)
            KLZ(KL) = NZ+LZ(L)
            KLGT(KL) = NGTK*(K-MINK)+NGTL*(L-MINL)
         ENDDO
      ENDDO
      MAX = KL
      !write(*,*) "IJ", IJ
!      DO 320 I = 1,IJ
      DO I = 1,IJ
      !write(*,*) "IJ", IJ
      IF (SAME) MAX = I
      IK(I) = MAX
      ENDDO
      !write(*,*) "IK", IK(I)
!320 IK(I) = MAX
      !write(*,*) "IK(I)", IK(I)
      IJKL = IJ*KL
      IF (SAME) IJKL = IJ*(IJ+1)/2
      RETURN
      END



      SUBROUTINE ZQOUT_gpu_2(GHONDO,&
      IJGT,KLGT,&
      IANDJ,KANDL,SAME,&
      mini,maxi,minj,maxj,mink,maxk,minl,maxl)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION GHONDO(*)
      LOGICAL IANDJ,KANDL,SAME
      DIMENSION :: IJGT(784),KLGT(784)
      INTEGER :: MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL
      PARAMETER (ZERO=0.0D+00)

!$omp declare target
      IJN = 0
      JMAX = MAXJ
      DO 260 I = MINI,MAXI
         IF (IANDJ) JMAX = I
         DO 240 J = MINJ,JMAX
            IJN = IJN+1
            N1 = IJGT(IJN)
            LMAX = MAXL
            KLN = 0
            DO 220 K =  MINK,MAXK
               IF (KANDL) LMAX = K
               DO 200 L = MINL,LMAX
                  KLN = KLN+1
                  !IF (SAME .AND. KLN .GT. IJN) GO TO 240
                  NN = N1+KLGT(KLN)
                  GHONDO(NN) = ZERO
  200          CONTINUE
  220       CONTINUE
  240    CONTINUE
  260 CONTINUE
      RETURN
      END

SUBROUTINE GENRAL_GPU_2(GHONDO,DDIJ,DKL,DIJ,&
      IJGT,IJX,IJY,IJZ,IK,KLGT,KLX,KLY,KLZ,&
      AA,R,X1,Y1,Z1,IJD,IANDJ,KANDL,SAME,&
      LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,&
      qq4,mini,maxi,minj,maxj,mink,maxk,minl,maxl,&
      NIJ,IJ,KL,IJKL,&
      ag,csa,cpa,cda,&
      bg,csb,cpb,cdb,&
      cg,csc,cpc,cdc,&
      dg,csd,cpd,cdd,&
      XI,YI,ZI,XJ,YJ,ZJ,RRI,XK,YK,ZK,XL,YL,ZL,RRK,&
      NGA,NGB,NGC,NGD, &
      XIN,YIN,ZIN)

      use omp_lib
      USE lrcdft, ONLY: EMU2
      use mx_limits, only: mxgsh,mxg2

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION GHONDO(*),DDIJ(*)
      LOGICAL IANDJ,KANDL,SAME,OUT,NORM,DOUBLE
      double precision :: XIN(31213),YIN(31213),ZIN(31213)
  
      DIMENSION :: DKL(784),DIJ(784)
      DIMENSION :: IJGT(784),IJX(784),IJY(784),IJZ(784),IK(784)
      DIMENSION :: KLGT(784),KLX(784),KLY(784),KLZ(784)
      DIMENSION :: AA(MXG2),R(MXG2),X1(MXG2),Y1(MXG2),Z1(MXG2),IJD(784)

      double precision :: XX
      double precision :: U(13)
      double precision :: W(13)
      DIMENSION :: IN(13),KN(13)
      INTEGER :: NI,NJ,NK,NL,NMAX,MMAX
      double precision :: BP01,B00,B10,XCP00,XC00,YCP00,YC00,ZCP00,ZC00
      double precision :: F00,DXIJ,DYIJ,DZIJ,DXKL,DYKL,DZKL

      DIMENSION :: AG(MXGSH),CSA(MXGSH),CPA(MXGSH),CDA(MXGSH)
      DIMENSION :: BG(MXGSH),CSB(MXGSH),CPB(MXGSH),CDB(MXGSH)
      DIMENSION :: CG(MXGSH),CSC(MXGSH),CPC(MXGSH),CDC(MXGSH)
      DIMENSION :: DG(MXGSH),CSD(MXGSH),CPD(MXGSH),CDD(MXGSH)

      double precision :: XI,YI,ZI,XJ,YJ,ZJ,RRI,XK,YK,ZK,XL,YL,ZL,RRK
      integer :: NGA,NGB,NGC,NGD

      double precision :: qq4
      INTEGER :: LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL
      INTEGER :: MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL
      INTEGER :: NIJ,IJ,KL,IJKL
      DIMENSION IN1(13)
      PARAMETER (SQRT3=1.73205080756888D+00, SQRT5=2.23606797749979D+00,&
                SQRT7=2.64575131106459D+00, PI252=34.986836655250D+00,&
                SQRT9=3.0D+00,SQRT11=3.3166247903553998D+00,&
                ZERO=0.0D+00, HALF=0.5D+00, ONE=1.0D+00,&
                TOL=46.0515999999999934 )
!$omp declare target

!new
! teamid = omp_get_num_teams()
! threadid = omp_get_num_threads()
! new = omp_get_thread_num()
! !write(*,*) "new", NEW
!ID = omp_get_thread_num()

      !write(*,*) "ID is", ID
      FACTOR = PI252*QQ4
      NI = LIT-1
      NJ = LJT-1
      NK = LKT-1
      NL = LLT-1
      DXIJ = XI-XJ
      DYIJ = YI-YJ
      DZIJ = ZI-ZJ
      DXKL = XK-XL
      DYKL = YK-YL
      DZKL = ZK-ZL
      NMAX = NI+NJ
      MMAX = NK+NL
      MAX = NMAX+1
      DO I = 1,MAX
         N = I-1
         IF (N .LE. NI) IN1(I) = 343*N+1
         IF (N .GT. NI) IN1(I) = 343*NI+49*(N-NI)+1
      ENDDO
      MAX = MMAX+1
      DO K = 1,MAX
         N = K-1
         IF (N .LE. NK) KN(K) = 7*N
         IF (N .GT. NK) KN(K) = 7*NK+N-NK
      ENDDO

      LGMAX = NGD
      DO 480 KG = 1,NGC
         AK = CG(KG)
         BRRK = AK*RRK
         AKXK = AK*XK
         AKYK = AK*YK
         AKZK = AK*ZK
         CSK = CSC(KG)*FACTOR
         CPK = CPC(KG)*FACTOR
         CDK = CDC(KG)*FACTOR
         IF (KANDL) LGMAX = KG
         DO 460 LG = 1,LGMAX
            AL = DG(LG)
            B = AK+AL
            BINV = ONE/B
            BBRRK = AL*BRRK*BINV
            CSL = CSD(LG)
            CPL = CPD(LG)
            CDL = CDD(LG)
            XB = (AKXK+AL*XL)*BINV
            YB = (AKYK+AL*YL)*BINV
            ZB = (AKZK+AL*ZL)*BINV
            BXBK = B*(XB-XK)
            BYBK = B*(YB-YK)
            BZBK = B*(ZB-ZK)
            BXBI = B*(XB-XI)
            BYBI = B*(YB-YI)
            BZBI = B*(ZB-ZI)

            DOUBLE=KANDL.AND.KG.NE.LG
            N = 0
            MAX = MAXL
            DUM1 = ZERO
            DUM2 = ZERO


            DO 370 K = MINK,MAXK
               GO TO (140,160,220,220,180,220,220,200,220,220,201),K

  140          DUM1 = CSK*BINV
               GO TO 220
  160          DUM1 = CPK*BINV
               GO TO 220
  180          DUM1 = CDK*BINV
               GO TO 220
  200          DUM1 = DUM1*SQRT3
               GO TO 220
  201          GO TO 220

  220          IF (KANDL) MAX = K
               DO 360 L = MINL,MAX
                  GO TO (240,280,340,340,300,340,340,320,340),L
  240             DUM2 = DUM1*CSL
                  IF ( .NOT. DOUBLE) GO TO 340
                  IF (K .GT. 1) GO TO 260
                  DUM2 = DUM2+DUM2
                  GO TO 340
  260             DUM2 = DUM2+CSK*CPL*BINV
                  GO TO 340
  280             DUM2 = DUM1*CPL
                  IF (DOUBLE) DUM2 = DUM2+DUM2
                  GO TO 340
  300             DUM2 = DUM1*CDL
                  IF (DOUBLE) DUM2 = DUM2+DUM2
                  GO TO 340
  320             DUM2 = DUM2*SQRT3
                  GO TO 340
  340             N = N+1
                  DKL(N) = DUM2
                  ! write(*,*) "DKL(N)", DKL(N)

  360          CONTINUE
  370       CONTINUE

            NN = 0
            DO 440 N = 1,NIJ
               DUM = BBRRK+R(N)
               DO I = 1,IJ
                  DIJ(I) = DDIJ(IJD(I)+NN)
               ENDDO
               A = AA(N)
               AB = A*B
               AANDB = A+B
               EXPE = EXP(-DUM)/SQRT(AANDB)
               RHO = AB/AANDB
               XA = X1(N)
               YA = Y1(N)
               ZA = Z1(N)
               XX = RHO*((XA-XB)*(XA-XB)+(YA-YB)*(YA-YB)+(ZA-ZB)*(ZA-ZB))
               AXAK = A*(XA-XK)
               AYAK = A*(YA-YK)
               AZAK = A*(ZA-ZK)
               AXAI = A*(XA-XI)
               AYAI = A*(YA-YI)
               AZAI = A*(ZA-ZI)
               C1X = BXBK+AXAK
               C2X = A*BXBK
               C3X = BXBI+AXAI
               C4X = B*AXAI
               C1Y = BYBK+AYAK
               C2Y = A*BYBK
               C3Y = BYBI+AYAI
               C4Y = B*AYAI
               C1Z = BZBK+AZAK
               C2Z = A*BZBK
               C3Z = BZBI+AZAI
               C4Z = B*AZAI

               CALL RT123TS_gpu(XX,2,RT1,RT2,RT3,WW1,WW2,WW3)
               U(1)=RT1
               U(2)=RT2
               U(3)=RT3
               W(1)=WW1
               W(2)=WW2
               W(3)=WW3

               MM = 0
               MAX = NMAX+1
               !THIS IS OK
               !write(*,*) "MAX is", MAX

               DO M = 1,2 !ROOTS VALUE 1 AND 2
                  U2 = U(M)*RHO
                  F00 = EXPE*W(M)
                  DO I = 1,MAX
                     IN(I) = IN1(I)+MM
                  ENDDO

                     DUMINV = ONE/(AB+U2*AANDB)
                     DM2INV = HALF*DUMINV
                     BP01 = (A+U2)*DM2INV
                     B00 = U2*DM2INV
                     B10 = (B+U2)*DM2INV
                     XCP00 = (U2*C1X+C2X)*DUMINV
                     XC00 = (U2*C3X+C4X)*DUMINV
                     YCP00 = (U2*C1Y+C2Y)*DUMINV
                     YC00 = (U2*C3Y+C4Y)*DUMINV
                     ZCP00 = (U2*C1Z+C2Z)*DUMINV
                     ZC00 = (U2*C3Z+C4Z)*DUMINV
                      !newer version 


                  CALL XYZINT_gpu_2(XIN,YIN,ZIN,&
                  IN,KN,NI,NJ,NK,NL,NMAX,MMAX,&
                  BP01,B00,B10,XCP00,XC00,YCP00,YC00,ZCP00,ZC00,F00,&
                  DXIJ,DYIJ,DZIJ,DXKL,DYKL,DZKL)

                  MM = MM+2401
               ENDDO

      DO I = 1,IJ
      D1 = DIJ(I)
      !write(*,*) "D1 is", D1
      NX = IJX(I)
      !write(*,*) "NX is", NX
      NY = IJY(I)
      NZ = IJZ(I)
      N1 = IJGT(I)
      MAX = IK(I)
      !write(*,*) "MAX is", MAX
      DO K = 1,MAX
      !write(*,*) "MAX is", MAX
      MX = NX+KLX(K)
      !write(*,*) "MX is", MX
      !write(*,*) "KLX is", KLX(K)
      !these should be good
      MY = NY+KLY(K)
      MZ = NZ+KLZ(K)
      NNN = N1+KLGT(K)
      !write(*,*) "NNN is", NNN
      !write(*,*) "DKL(K) is", DKL(K)
      !write(*,*) "D1*DKL(K) is", D1*DKL(K)
      !write(*,*) "XIN is", XIN(MX,ID)
      !!$omp atomic
      !GHONDO(NNN) = GHONDO(NNN) + D1*DKL(K)*( varx(MX,ID)*vary(MY,ID)*varz(MZ,ID)+ varx(MX+ 2401,ID)*vary(MY+ 2401,ID)*varz(MZ+ 2401,ID))
      GHONDO(NNN) = GHONDO(NNN) + D1*DKL(K)*( XIN(MX)*YIN(MY)*ZIN(MZ)+ XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401))
      !GHONDO(NNN) = GHONDO(NNN) + D1*DKL(K)*( XIN(MX,ID)*YIN(MY,ID)*ZIN(MZ,ID)+ XIN(MX+ 2401,ID)*YIN(MY+ 2401,ID)*ZIN(MZ+ 2401,ID))
      !GHONDO(NNN) = GHONDO(NNN) + D1*DKL(K)
      !if (GHONDO(NNN) .NE. 0.00) write(*,*) GHONDO(NNN)
      !write(*,*) "GHONDO  is", GHONDO(NNN)
      ENDDO
      ENDDO

  440       NN = NN+49
  460    CONTINUE
  480 CONTINUE
      RETURN
      END


      SUBROUTINE GENRAL_gpu_2222(GHONDO,DDIJ,DKL,DIJ,&
      IJGT,IJX,IJY,IJZ,IK,KLGT,KLX,KLY,KLZ,&
      AA,R,X1,Y1,Z1,IJD,IANDJ,KANDL,SAME,&
      LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,&
      qq4,mini,maxi,minj,maxj,mink,maxk,minl,maxl,&
      NIJ,IJ,KL,IJKL,&
      NROOTS,&
      cg,csc,cpc,cdc,&
      dg,csd,cpd,cdd,&
      XI,YI,ZI,XJ,YJ,ZJ,RRI,XK,YK,ZK,XL,YL,ZL,RRK)


      USE lrcdft, ONLY: EMU2
      use mx_limits, only: mxgsh,mxg2

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION GHONDO(*),DDIJ(*)

      LOGICAL IANDJ,KANDL,SAME,OUT,NORM,DOUBLE
      DIMENSION :: XIN(31213),YIN(31213),ZIN(31213)
      DIMENSION :: DKL(784),DIJ(784)
      DIMENSION :: IJGT(784),IJX(784),IJY(784),IJZ(784),IK(784)
      DIMENSION :: KLGT(784),KLX(784),KLY(784),KLZ(784)
      DIMENSION :: AA(MXG2),R(MXG2),X1(MXG2),Y1(MXG2),Z1(MXG2),IJD(784)

      double precision :: XX
      double precision :: U(13)
      double precision :: W(13)
      ! double precision :: U(5)
      ! double precision :: W(5)
      integer :: NROOTS
      DIMENSION :: IN(13),KN(13)
      INTEGER :: NI,NJ,NK,NL,NMAX,MMAX
      double precision :: BP01,B00,B10,XCP00,XC00,YCP00,YC00,ZCP00,ZC00
      double precision :: F00,DXIJ,DYIJ,DZIJ,DXKL,DYKL,DZKL

      DIMENSION :: CG(MXGSH),CSC(MXGSH),CPC(MXGSH),CDC(MXGSH)
      DIMENSION :: DG(MXGSH),CSD(MXGSH),CPD(MXGSH),CDD(MXGSH)

      double precision :: XI,YI,ZI,XJ,YJ,ZJ,RRI,XK,YK,ZK,XL,YL,ZL,RRK
      !integer :: NGA,NGB,NGC,NGD

      double precision :: qq4
      INTEGER :: LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL
      INTEGER :: MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL
      INTEGER :: NIJ,IJ,KL,IJKL

      DIMENSION IN1(13)

      PARAMETER (SQRT3=1.73205080756888D+00, SQRT5=2.23606797749979D+00,&
                SQRT7=2.64575131106459D+00, PI252=34.986836655250D+00,&
                SQRT9=3.0D+00,SQRT11=3.3166247903553998D+00,&
                ZERO=0.0D+00, HALF=0.5D+00, ONE=1.0D+00,&
                TOL=46.0515999999999934 )

      ! DATA R15,PIE4/1.17581320211778D-01, 7.85398163397448D-01/
      ! DATA R25,W25/ 1.07456201243690D+00, 2.70967405960535D-01/
      ! DATA R35,W35/ 3.08593744371754D+00, 3.82231610015404D-02/
      ! DATA R45,W45/ 6.41472973366203D+00, 1.51614186862443D-03/
      ! DATA R55,W55/ 1.18071894899717D+01, 8.62130526143657D-06/
!$omp declare target
      FACTOR = PI252*QQ4
      NI = LIT-1
      NJ = LJT-1
      NK = LKT-1
      NL = LLT-1
      DXIJ = XI-XJ
      DYIJ = YI-YJ
      DZIJ = ZI-ZJ
      DXKL = XK-XL
      DYKL = YK-YL
      DZKL = ZK-ZL
      NMAX = NI+NJ
      MMAX = NK+NL
      MAX = NMAX+1
      DO 100 I = 1,MAX
         N = I-1
         IF (N .LE. NI) IN1(I) = 343*N+1
         IF (N .GT. NI) IN1(I) = 343*NI+49*(N-NI)+1
  100 CONTINUE
      MAX = MMAX+1
      DO 120 K = 1,MAX
         N = K-1
         IF (N .LE. NK) KN(K) = 7*N
         IF (N .GT. NK) KN(K) = 7*NK+N-NK
  120 CONTINUE

      !LGMAX = NGD
      !DO 480 KG = 1,NGC
         AK = CG(KG)
         BRRK = AK*RRK
         AKXK = AK*XK
         AKYK = AK*YK
         AKZK = AK*ZK
         CSK = CSC(KG)*FACTOR
         CPK = CPC(KG)*FACTOR
         CDK = CDC(KG)*FACTOR

         !IF (KANDL) LGMAX = KG
         !DO 460 LG = 1,LGMAX
            AL = DG(LG)
            B = AK+AL
            BINV = ONE/B
            BBRRK = AL*BRRK*BINV
            !IF (BBRRK .GT. TOL) GO TO 460
            CSL = CSD(LG)
            CPL = CPD(LG)
            CDL = CDD(LG)

            XB = (AKXK+AL*XL)*BINV
            YB = (AKYK+AL*YL)*BINV
            ZB = (AKZK+AL*ZL)*BINV
            BXBK = B*(XB-XK)
            BYBK = B*(YB-YK)
            BZBK = B*(ZB-ZK)
            BXBI = B*(XB-XI)
            BYBI = B*(YB-YI)
            BZBI = B*(ZB-ZI)

            DOUBLE=KANDL.AND.KG.NE.LG
            N = 0
            MAX = MAXL
            DUM1 = ZERO
            DUM2 = ZERO

            DO 370 K = MINK,MAXK
               GO TO (140,160,220,220,180,220,220,200,220,220,&
                     201),K
  140          DUM1 = CSK*BINV
               GO TO 220
  160          DUM1 = CPK*BINV
               GO TO 220
  180          DUM1 = CDK*BINV
               GO TO 220
  200          DUM1 = DUM1*SQRT3
               GO TO 220
  201          GO TO 220

  220          IF (KANDL) MAX = K
               DO 360 L = MINL,MAX
                  GO TO (240,280,340,340,300,340,340,320,340,340,&
                        321),L
  240             DUM2 = DUM1*CSL
                  IF ( .NOT. DOUBLE) GO TO 340
                  IF (K .GT. 1) GO TO 260
                  DUM2 = DUM2+DUM2
                  GO TO 340
  260             DUM2 = DUM2+CSK*CPL*BINV
                  GO TO 340
  280             DUM2 = DUM1*CPL
                  IF (DOUBLE) DUM2 = DUM2+DUM2
                  GO TO 340
  300             DUM2 = DUM1*CDL
                  IF (DOUBLE) DUM2 = DUM2+DUM2
                  GO TO 340
  320             DUM2 = DUM2*SQRT3
                  GO TO 340
  321             GO TO 340

  340             N = N+1
                  DKL(N) = DUM2
  360          CONTINUE
  370       CONTINUE

            NN = 0
            DO 440 N = 1,NIJ
               DUM = BBRRK+R(N)
               IF (DUM .GT. TOL) GO TO 440
               DO 380 I = 1,IJ
                  DIJ(I) = DDIJ(IJD(I)+NN)
  380          CONTINUE
               A = AA(N)
               AB = A*B
               AANDB = A+B
               EXPE = EXP(-DUM)/SQRT(AANDB)
               RHO = AB/AANDB
               XA = X1(N)
               YA = Y1(N)
               ZA = Z1(N)
               XX = RHO*((XA-XB)*(XA-XB) + (YA-YB)*(YA-YB)+ (ZA-ZB)*(ZA-ZB))
               !write(*,*) "XX in 2222 is", XX
               AXAK = A*(XA-XK)
               AYAK = A*(YA-YK)
               AZAK = A*(ZA-ZK)
               AXAI = A*(XA-XI)
               AYAI = A*(YA-YI)
               AZAI = A*(ZA-ZI)
               C1X = BXBK+AXAK
               C2X = A*BXBK
               C3X = BXBI+AXAI
               C4X = B*AXAI
               C1Y = BYBK+AYAK
               C2Y = A*BYBK
               C3Y = BYBI+AYAI
               C4Y = B*AYAI
               C1Z = BZBK+AZAK
               C2Z = A*BZBK
               C3Z = BZBI+AZAI
               C4Z = B*AZAI

                   CALL ROOT5TS_GPU(XX,NROOTS,RT1,RT2,RT3,RT4,RT5,&
      WW1,WW2,WW3,WW4,WW5)
                   U(1)=RT1
                   U(2)=RT2
                   U(3)=RT3
                   U(4)=RT4
                   U(5)=RT5
                   W(1)=WW1
                   W(2)=WW2
                   W(3)=WW3
                   W(4)=WW4
                   W(5)=WW5

               MM = 0
               MAX = NMAX+1

               DO 420 M = 1,NROOTS
! C               DO 420 M = 1,5
                  U2 = U(M)*RHO
                  F00 = EXPE*W(M)
                  DO 400 I = 1,MAX
                     IN(I) = IN1(I)+MM
  400             CONTINUE

                     DUMINV = ONE/(AB+U2*AANDB)
                     DM2INV = HALF*DUMINV
                     BP01 = (A+U2)*DM2INV
                     B00 = U2*DM2INV
                     B10 = (B+U2)*DM2INV
                     XCP00 = (U2*C1X+C2X)*DUMINV
                     XC00 = (U2*C3X+C4X)*DUMINV
                     YCP00 = (U2*C1Y+C2Y)*DUMINV
                     YC00 = (U2*C3Y+C4Y)*DUMINV
                     ZCP00 = (U2*C1Z+C2Z)*DUMINV
                     ZC00 = (U2*C3Z+C4Z)*DUMINV

                  CALL XYZINT_gpu_2222(XIN,YIN,ZIN,&
      IN,KN,NI,NJ,NK,NL,NMAX,MMAX,&
      BP01,B00,B10,XCP00,XC00,YCP00,YC00,ZCP00,ZC00,F00,&
      DXIJ,DYIJ,DZIJ,DXKL,DYKL,DZKL)
                  MM = MM+2401
  420          CONTINUE

      DO I = 1,IJ
      D1 = DIJ(I)
      NX = IJX(I)
      NY = IJY(I)
      NZ = IJZ(I)
      N1 = IJGT(I)
      MAX = IK(I)
      DO K = 1,MAX
      MX = NX+KLX(K)
      MY = NY+KLY(K)
      MZ = NZ+KLZ(K)
      NNN = N1+KLGT(K)
      !GHONDO(NNN) = GHONDO(NNN) + D1*DKL(K)*( XIN(MX)*YIN(MY)*ZIN(MZ)+ XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401) &      
      !      + XIN(MX+ 4802)*YIN(MY+ 4802)*ZIN(MZ+ 4802)+ XIN(MX+ 7203)*YIN(MY+ 7203)*ZIN(MZ+ 7203)+XIN(MX+ 9604)*YIN(MY+ 9604)*ZIN(MZ+ 9604))
      ENDDO
      ENDDO


  440       NN = NN+49
  !460    CONTINUE
  !480 CONTINUE

      RETURN
      END
!       SUBROUTINE GENRAL_gpu_2222(GHONDO,DDIJ,DKL,DIJ,&
!       IJGT,IJX,IJY,IJZ,IK,KLGT,KLX,KLY,KLZ,&
!       AA,R,X1,Y1,Z1,IJD,IANDJ,KANDL,SAME,&
!       LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,&
!       qq4,mink,maxk,minl,maxl,&
!       NIJ,IJ,KL,IJKL,&
!       NROOTS,&
!       cg,csc,cpc,cdc,&
!       dg,csd,cpd,cdd,&
!       XI,YI,ZI,XJ,YJ,ZJ,RRI,XK,YK,ZK,XL,YL,ZL,RRK,&
!       NGA,NGB,NGC,NGD)

!       use mx_limits, only: mxgsh,mxg2

!       IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!       DIMENSION GHONDO(*),DDIJ(*)

!       LOGICAL IANDJ,KANDL,SAME,OUT,NORM,DOUBLE
!       DIMENSION :: XIN(31213),YIN(31213),ZIN(31213)
!       DIMENSION :: DKL(784),DIJ(784)
!       DIMENSION :: IJGT(784),IJX(784),IJY(784),IJZ(784),IK(784)
!       DIMENSION :: KLGT(784),KLX(784),KLY(784),KLZ(784)
!       DIMENSION :: AA(MXG2),R(MXG2),X1(MXG2),Y1(MXG2),Z1(MXG2),IJD(784)

!       double precision :: XX
!       double precision :: U(13)
!       double precision :: W(13)
!       integer :: NROOTS
!       DIMENSION :: IN(13),KN(13)
!       INTEGER :: NI,NJ,NK,NL,NMAX,MMAX
!       double precision :: BP01,B00,B10,XCP00,XC00,YCP00,YC00,ZCP00,ZC00
!       double precision :: F00,DXIJ,DYIJ,DZIJ,DXKL,DYKL,DZKL

!       !DIMENSION :: AG(MXGSH),CSA(MXGSH),CPA(MXGSH),CDA(MXGSH)
!       !DIMENSION :: BG(MXGSH),CSB(MXGSH),CPB(MXGSH),CDB(MXGSH)
!       DIMENSION :: CG(MXGSH),CSC(MXGSH),CPC(MXGSH),CDC(MXGSH)
!       DIMENSION :: DG(MXGSH),CSD(MXGSH),CPD(MXGSH),CDD(MXGSH)

!       double precision :: XI,YI,ZI,XJ,YJ,ZJ,RRI,XK,YK,ZK,XL,YL,ZL,RRK
!       integer :: NGA,NGB,NGC,NGD

!       double precision :: qq4
!       INTEGER :: LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL
!       INTEGER :: MINK,MINL,MAXK,MAXL
!       INTEGER :: NIJ,IJ,KL,IJKL

!       !LOGICAL FIRST1,FIRST2,FIRST3,FIRST4

!       DIMENSION IN1(13)

!       PARAMETER (SQRT3=1.73205080756888D+00, SQRT5=2.23606797749979D+00,&
!                 SQRT7=2.64575131106459D+00, PI252=34.986836655250D+00,&
!                 SQRT9=3.0D+00,SQRT11=3.3166247903553998D+00,&
!                 ZERO=0.0D+00, HALF=0.5D+00, ONE=1.0D+00,&
!                 TOL=46.0515999999999934 )
      
!       ! DATA R15,PIE4/1.17581320211778D-01, 7.85398163397448D-01/
!       ! DATA R25,W25/ 1.07456201243690D+00, 2.70967405960535D-01/
!       ! DATA R35,W35/ 3.08593744371754D+00, 3.82231610015404D-02/
!       ! DATA R45,W45/ 6.41472973366203D+00, 1.51614186862443D-03/
!       ! DATA R55,W55/ 1.18071894899717D+01, 8.62130526143657D-06/
! !$omp declare target
!       FACTOR = PI252*QQ4
!       NI = LIT-1
!       NJ = LJT-1
!       NK = LKT-1
!       NL = LLT-1
!       DXIJ = XI-XJ
!       DYIJ = YI-YJ
!       DZIJ = ZI-ZJ
!       DXKL = XK-XL
!       DYKL = YK-YL
!       DZKL = ZK-ZL
!       NMAX = NI+NJ
!       MMAX = NK+NL
!       MAX = NMAX+1
!       DO 100 I = 1,MAX
!          N = I-1
!          IF (N .LE. NI) IN1(I) = 343*N+1
!          IF (N .GT. NI) IN1(I) = 343*NI+49*(N-NI)+1
!   100 CONTINUE
!       MAX = MMAX+1
!       DO 120 K = 1,MAX
!          N = K-1
!          IF (N .LE. NK) KN(K) = 7*N
!          IF (N .GT. NK) KN(K) = 7*NK+N-NK
!   120 CONTINUE

!       LGMAX = NGD
!       DO 480 KG = 1,NGC
!          AK = CG(KG)
!          BRRK = AK*RRK
!          AKXK = AK*XK
!          AKYK = AK*YK
!          AKZK = AK*ZK
!          CSK = CSC(KG)*FACTOR
!          CPK = CPC(KG)*FACTOR
!          CDK = CDC(KG)*FACTOR


!          IF (KANDL) LGMAX = KG
!          DO 460 LG = 1,LGMAX
!             AL = DG(LG)
!             B = AK+AL
!             BINV = ONE/B
!             BBRRK = AL*BRRK*BINV
!             IF (BBRRK .GT. TOL) GO TO 460
!             CSL = CSD(LG)
!             CPL = CPD(LG)
!             CDL = CDD(LG)

!             XB = (AKXK+AL*XL)*BINV
!             YB = (AKYK+AL*YL)*BINV
!             ZB = (AKZK+AL*ZL)*BINV
!             BXBK = B*(XB-XK)
!             BYBK = B*(YB-YK)
!             BZBK = B*(ZB-ZK)
!             BXBI = B*(XB-XI)
!             BYBI = B*(YB-YI)
!             BZBI = B*(ZB-ZI)

!             DOUBLE=KANDL.AND.KG.NE.LG
!             N = 0
!             MAX = MAXL
!             DUM1 = ZERO
!             DUM2 = ZERO

!             DO 370 K = MINK,MAXK
!                GO TO (140,160,220,220,180,220,220,200,220,220,&
!                      201),K
!   140          DUM1 = CSK*BINV
!                GO TO 220
!   160          DUM1 = CPK*BINV
!                GO TO 220
!   180          DUM1 = CDK*BINV
!                GO TO 220
!   200          DUM1 = DUM1*SQRT3
!                GO TO 220
!   201          GO TO 220

!   220          IF (KANDL) MAX = K
!                DO 360 L = MINL,MAX
!                   GO TO (240,280,340,340,300,340,340,320,340,340,&
!                         321),L
!   240             DUM2 = DUM1*CSL
!                   IF ( .NOT. DOUBLE) GO TO 340
!                   IF (K .GT. 1) GO TO 260
!                   DUM2 = DUM2+DUM2
!                   GO TO 340
!   260             DUM2 = DUM2+CSK*CPL*BINV
!                   GO TO 340
!   280             DUM2 = DUM1*CPL
!                   IF (DOUBLE) DUM2 = DUM2+DUM2
!                   GO TO 340
!   300             DUM2 = DUM1*CDL
!                   IF (DOUBLE) DUM2 = DUM2+DUM2
!                   GO TO 340
!   320             DUM2 = DUM2*SQRT3
!                   GO TO 340
!   321             GO TO 340

!   340             N = N+1
!                   DKL(N) = DUM2
!   360          CONTINUE
!   370       CONTINUE

!             NN = 0
!             DO 440 N = 1,NIJ
!                DUM = BBRRK+R(N)
!                IF (DUM .GT. TOL) GO TO 440
!                DO 380 I = 1,IJ
!                   DIJ(I) = DDIJ(IJD(I)+NN)
!   380          CONTINUE
!                A = AA(N)
!                AB = A*B
!                AANDB = A+B
!                EXPE = EXP(-DUM)/SQRT(AANDB)
!                RHO = AB/AANDB
!                XA = X1(N)
!                YA = Y1(N)
!                ZA = Z1(N)
!                XX = RHO*((XA-XB)*(XA-XB) + (YA-YB)*(YA-YB)+ (ZA-ZB)*(ZA-ZB))
!                !write(*,*) "XX in 2222 is", XX
!                AXAK = A*(XA-XK)
!                AYAK = A*(YA-YK)
!                AZAK = A*(ZA-ZK)
!                AXAI = A*(XA-XI)
!                AYAI = A*(YA-YI)
!                AZAI = A*(ZA-ZI)
!                C1X = BXBK+AXAK
!                C2X = A*BXBK
!                C3X = BXBI+AXAI
!                C4X = B*AXAI
!                C1Y = BYBK+AYAK
!                C2Y = A*BYBK
!                C3Y = BYBI+AYAI
!                C4Y = B*AYAI
!                C1Z = BZBK+AZAK
!                C2Z = A*BZBK
!                C3Z = BZBI+AZAI
!                C4Z = B*AZAI


! ! !!$omp declare target
! !       IF (XX .GT. 15.0D+00) GO TO 1800
! !       IF (XX .GT. 5.0D+00) GO TO 1400
! !       IF (XX .GT. 1.0D+00) GO TO 1200
! !       IF (XX .GT. 3.0D-07) GO TO 1000
! ! !C     X IS APPROXIMATELY ZERO.                   NROOTS = 5
! !       RT1 = 2.26659266316985D-02 -2.15865967920897D-03 *XX
! !       RT2 = 2.31271692140903D-01 -2.20258754389745D-02 *XX
! !       RT3 = 8.57346024118836D-01 -8.16520023025515D-02 *XX
! !       RT4 = 2.97353038120346D+00 -2.83193369647137D-01 *XX
! !       RT5 = 1.84151859759051D+01 -1.75382723579439D+00 *XX
! !       WW1 = 2.95524224714752D-01 -1.96867576909777D-02 *XX
! !       WW2 = 2.69266719309995D-01 -5.61737590184721D-02 *XX
! !       WW3 = 2.19086362515981D-01 -9.71152726793658D-02 *XX
! !       WW4 = 1.49451349150580D-01 -1.02979262193565D-01 *XX
! !       WW5 = 6.66713443086877D-02 -5.73782817488315D-02 *XX
! !       RETURN

! ! !C     X=0.0 TO 1.0                               NROOTS = 5
! !   1000 RT1 = ((((((-4.46679165328413D-11*XX+1.21879111988031D-09)*XX- &
! !           2.62975022612104D-08 )*XX+5.15106194905897D-07 )*XX- &
! !           9.27933625824749D-06 )*XX+1.51794097682482D-04 )*XX- &
! !           2.15865967920301D-03 )*XX+2.26659266316985D-02
! !       RT2 = (((((( 1.93117331714174D-10*XX-4.57267589660699D-09)*XX+ &
! !           2.48339908218932D-08 )*XX+1.50716729438474D-06 )*XX- &
! !           6.07268757707381D-05 )*XX+1.37506939145643D-03 )*XX- &
! !           2.20258754419939D-02 )*XX+2.31271692140905D-01
! !       RT3 = ((((( 4.84989776180094D-09*XX+1.31538893944284D-07)*XX- &
! !           2.766753852879D-06)*XX-7.651163510626D-05)*XX+ &
! !           4.033058545972D-03)*XX-8.16520022916145D-02 )*XX+ &
! !           8.57346024118779D-01
! !       RT4 = ((((-2.48581772214623D-07*XX-4.34482635782585D-06)*XX- &
! !           7.46018257987630D-07 )*XX+1.01210776517279D-02 )*XX- &
! !           2.83193369640005D-01 )*XX+2.97353038120345D+00
! !       RT5 = (((((-8.92432153868554D-09*XX+1.77288899268988D-08)*XX+ &
! !           3.040754680666D-06)*XX+1.058229325071D-04)*XX+ &
! !           4.596379534985D-02)*XX-1.75382723579114D+00 )*XX+ &
! !           1.84151859759049D+01
! !       WW1 = ((((((-2.03822632771791D-09*XX+3.89110229133810D-08)*XX- &
! !           5.84914787904823D-07 )*XX+8.30316168666696D-06 )*XX- &
! !           1.13218402310546D-04 )*XX+1.49128888586790D-03 )*XX- &
! !           1.96867576904816D-02 )*XX+2.95524224714749D-01
! !       WW2 = ((((((( 8.62848118397570D-09*XX-1.38975551148989D-07)*XX+ &
! !           1.602894068228D-06)*XX-1.646364300836D-05)*XX+ &
! !           1.538445806778D-04)*XX-1.28848868034502D-03 )*XX+ &
! !           9.38866933338584D-03 )*XX-5.61737590178812D-02 )*XX+ &
! !           2.69266719309991D-01
! !       WW3 = ((((((((-9.41953204205665D-09*XX+1.47452251067755D-07)*XX- &
! !           1.57456991199322D-06 )*XX+1.45098401798393D-05 )*XX- &
! !           1.18858834181513D-04 )*XX+8.53697675984210D-04 )*XX- &
! !           5.22877807397165D-03 )*XX+2.60854524809786D-02 )*XX- &
! !           9.71152726809059D-02 )*XX+2.19086362515979D-01
! !       WW4 = ((((((((-3.84961617022042D-08*XX+5.66595396544470D-07)*XX- &
! !           5.52351805403748D-06 )*XX+4.53160377546073D-05 )*XX- &
! !           3.22542784865557D-04 )*XX+1.95682017370967D-03 )*XX- &
! !           9.77232537679229D-03 )*XX+3.79455945268632D-02 )*XX- &
! !           1.02979262192227D-01 )*XX+1.49451349150573D-01
! !       WW5 = ((((((((( 4.09594812521430D-09*XX-6.47097874264417D-08)*XX+ &
! !           6.743541482689D-07)*XX-5.917993920224D-06)*XX+ &
! !           4.531969237381D-05)*XX-2.99102856679638D-04 )*XX+ &
! !           1.65695765202643D-03 )*XX-7.40671222520653D-03 )*XX+ &
! !           2.50889946832192D-02 )*XX-5.73782817487958D-02 )*XX+ &
! !           6.66713443086877D-02
! !       RETURN

! ! !C     XX=1.0 TO 5.0                               NROOTS = 5
! !   1200 Y = XX-3.0D+00
! !       RT1 = ((((((((-2.58163897135138D-14*Y+8.14127461488273D-13)*Y- &
! !           2.11414838976129D-11 )*Y+5.09822003260014D-10 )*Y- &
! !           1.16002134438663D-08 )*Y+2.46810694414540D-07 )*Y- &
! !           4.92556826124502D-06 )*Y+9.02580687971053D-05 )*Y- &
! !           1.45190025120726D-03 )*Y+1.73416786387475D-02
! !       RT2 = ((((((((( 1.04525287289788D-14*Y+5.44611782010773D-14)*Y- &
! !           4.831059411392D-12)*Y+1.136643908832D-10)*Y- &
! !           1.104373076913D-09)*Y-2.35346740649916D-08 )*Y+ &
! !           1.43772622028764D-06 )*Y-4.23405023015273D-05 )*Y+ &
! !           9.12034574793379D-04 )*Y-1.52479441718739D-02 )*Y+ &
! !           1.76055265928744D-01
! !       RT3 = (((((((((-6.89693150857911D-14*Y+5.92064260918861D-13)*Y+ &
! !           1.847170956043D-11)*Y-3.390752744265D-10)*Y- &
! !           2.995532064116D-09)*Y+1.57456141058535D-07 )*Y- &
! !           3.95859409711346D-07 )*Y-9.58924580919747D-05 )*Y+ &
! !           3.23551502557785D-03 )*Y-5.97587007636479D-02 )*Y+ &
! !           6.46432853383057D-01
! !       RT4 = ((((((((-3.61293809667763D-12*Y-2.70803518291085D-11)*Y+ &
! !           8.83758848468769D-10 )*Y+1.59166632851267D-08 )*Y- &
! !           1.32581997983422D-07 )*Y-7.60223407443995D-06 )*Y- &
! !           7.41019244900952D-05 )*Y+9.81432631743423D-03 )*Y- &
! !           2.23055570487771D-01 )*Y+2.21460798080643D+00
! !       RT5 = ((((((((( 7.12332088345321D-13*Y+3.16578501501894D-12)*Y- &
! !           8.776668218053D-11)*Y-2.342817613343D-09)*Y- &
! !           3.496962018025D-08)*Y-3.03172870136802D-07 )*Y+ &
! !           1.50511293969805D-06 )*Y+1.37704919387696D-04 )*Y+ &
! !           4.70723869619745D-02 )*Y-1.47486623003693D+00 )*Y+ &
! !           1.35704792175847D+01
! !       WW1 = ((((((((( 1.04348658616398D-13*Y-1.94147461891055D-12)*Y+ &
! !           3.485512360993D-11)*Y-6.277497362235D-10)*Y+ &
! !           1.100758247388D-08)*Y-1.88329804969573D-07 )*Y+ &
! !           3.12338120839468D-06 )*Y-5.04404167403568D-05 )*Y+ &
! !           8.00338056610995D-04 )*Y-1.30892406559521D-02 )*Y+ &
! !           2.47383140241103D-01
! !       WW2 = ((((((((((( 3.23496149760478D-14*Y-5.24314473469311D-13)*Y+ &
! !           7.743219385056D-12)*Y-1.146022750992D-10)*Y+ &
! !           1.615238462197D-09)*Y-2.15479017572233D-08 )*Y+ &
! !           2.70933462557631D-07 )*Y-3.18750295288531D-06 )*Y+ &
! !           3.47425221210099D-05 )*Y-3.45558237388223D-04 )*Y+ &
! !           3.05779768191621D-03 )*Y-2.29118251223003D-02 )*Y+ &
! !           1.59834227924213D-01
! !       WW3 = ((((((((((((-3.42790561802876D-14*Y+5.26475736681542D-13)*Y- &
! !           7.184330797139D-12)*Y+9.763932908544D-11)*Y- &
! !           1.244014559219D-09)*Y+1.472744068942D-08)*Y- &
! !           1.611749975234D-07)*Y+1.616487851917D-06)*Y- &
! !           1.46852359124154D-05 )*Y+1.18900349101069D-04 )*Y- &
! !           8.37562373221756D-04 )*Y+4.93752683045845D-03 )*Y- &
! !           2.25514728915673D-02 )*Y+6.95211812453929D-02
! !       WW4 = ((((((((((((( 1.04072340345039D-14*Y-1.60808044529211D-13)* &
! !           Y+2.183534866798D-12)*Y-2.939403008391D-11)*Y+ &
! !           3.679254029085D-10)*Y-4.23775673047899D-09 )*Y+ &
! !           4.46559231067006D-08 )*Y-4.26488836563267D-07 )*Y+ &
! !           3.64721335274973D-06 )*Y-2.74868382777722D-05 )*Y+ &
! !           1.78586118867488D-04 )*Y-9.68428981886534D-04 )*Y+ &
! !           4.16002324339929D-03 )*Y-1.28290192663141D-02 )*Y+ &
! !           2.22353727685016D-02
! !       WW5 = ((((((((((((((-8.16770412525963D-16*Y+1.31376515047977D-14)* &
! !           Y-1.856950818865D-13)*Y+2.596836515749D-12)*Y- &
! !           3.372639523006D-11)*Y+4.025371849467D-10)*Y- &
! !           4.389453269417D-09)*Y+4.332753856271D-08)*Y- &
! !           3.82673275931962D-07 )*Y+2.98006900751543D-06 )*Y- &
! !           2.00718990300052D-05 )*Y+1.13876001386361D-04 )*Y- &
! !           5.23627942443563D-04 )*Y+1.83524565118203D-03 )*Y- &
! !           4.37785737450783D-03 )*Y+5.36963805223095D-03
! !       RETURN

! !   1400 IF (XX .GT. 10.0D+00) GO TO 1600
! ! !C     XX=5.0 TO 10.0                              NROOTS = 5
! !       Y = XX-7.5D+00
! !       RT1 = ((((((((-1.13825201010775D-14*Y+1.89737681670375D-13)*Y- &
! !           4.81561201185876D-12 )*Y+1.56666512163407D-10 )*Y- &
! !           3.73782213255083D-09 )*Y+9.15858355075147D-08 )*Y- &
! !           2.13775073585629D-06 )*Y+4.56547356365536D-05 )*Y- &
! !           8.68003909323740D-04 )*Y+1.22703754069176D-02
! !       RT2 = (((((((((-3.67160504428358D-15*Y+1.27876280158297D-14)*Y- &
! !           1.296476623788D-12)*Y+1.477175434354D-11)*Y+ &
! !           5.464102147892D-10)*Y-2.42538340602723D-08 )*Y+ &
! !           8.20460740637617D-07 )*Y-2.20379304598661D-05 )*Y+ &
! !           4.90295372978785D-04 )*Y-9.14294111576119D-03 )*Y+ &
! !           1.22590403403690D-01
! !       RT3 = ((((((((( 1.39017367502123D-14*Y-6.96391385426890D-13)*Y+ &
! !           1.176946020731D-12)*Y+1.725627235645D-10)*Y- &
! !           3.686383856300D-09)*Y+2.87495324207095D-08 )*Y+ &
! !           1.71307311000282D-06 )*Y-7.94273603184629D-05 )*Y+ &
! !           2.00938064965897D-03 )*Y-3.63329491677178D-02 )*Y+ &
! !           4.34393683888443D-01
! !       RT4 = ((((((((((-1.27815158195209D-14*Y+1.99910415869821D-14)*Y+ &
! !           3.753542914426D-12)*Y-2.708018219579D-11)*Y- &
! !           1.190574776587D-09)*Y+1.106696436509D-08)*Y+ &
! !           3.954955671326D-07)*Y-4.398596059588D-06)*Y- &
! !           2.01087998907735D-04 )*Y+7.89092425542937D-03 )*Y- &
! !           1.42056749162695D-01 )*Y+1.39964149420683D+00
! !       RT5 = ((((((((((-1.19442341030461D-13*Y-2.34074833275956D-12)*Y+ &
! !           6.861649627426D-12)*Y+6.082671496226D-10)*Y+ &
! !           5.381160105420D-09)*Y-6.253297138700D-08)*Y- &
! !           2.135966835050D-06)*Y-2.373394341886D-05)*Y+ &
! !           2.88711171412814D-06 )*Y+4.85221195290753D-02 )*Y- &
! !           1.04346091985269D+00 )*Y+7.89901551676692D+00
! !       WW1 = ((((((((( 7.95526040108997D-15*Y-2.48593096128045D-13)*Y+ &
! !           4.761246208720D-12)*Y-9.535763686605D-11)*Y+ &
! !           2.225273630974D-09)*Y-4.49796778054865D-08 )*Y+ &
! !           9.17812870287386D-07 )*Y-1.86764236490502D-05 )*Y+ &
! !           3.76807779068053D-04 )*Y-8.10456360143408D-03 )*Y+ &
! !           2.01097936411496D-01
! !       WW2 = ((((((((((( 1.25678686624734D-15*Y-2.34266248891173D-14)*Y+ & 
! !           3.973252415832D-13)*Y-6.830539401049D-12)*Y+ &
! !           1.140771033372D-10)*Y-1.82546185762009D-09 )*Y+ &
! !           2.77209637550134D-08 )*Y-4.01726946190383D-07 )*Y+ &
! !           5.48227244014763D-06 )*Y-6.95676245982121D-05 )*Y+ &
! !           8.05193921815776D-04 )*Y-8.15528438784469D-03 )*Y+ &
! !           9.71769901268114D-02
! !       WW3 = ((((((((((((-8.20929494859896D-16*Y+1.37356038393016D-14)*Y- &
! !           2.022863065220D-13)*Y+3.058055403795D-12)*Y- &
! !           4.387890955243D-11)*Y+5.923946274445D-10)*Y- &
! !           7.503659964159D-09)*Y+8.851599803902D-08)*Y- &
! !           9.65561998415038D-07 )*Y+9.60884622778092D-06 )*Y- &
! !           8.56551787594404D-05 )*Y+6.66057194311179D-04 )*Y- &
! !           4.17753183902198D-03 )*Y+2.25443826852447D-02
! !       WW4 = ((((((((((((((-1.08764612488790D-17*Y+1.85299909689937D-16)* &
! !           Y-2.730195628655D-15)*Y+4.127368817265D-14)*Y- &
! !           5.881379088074D-13)*Y+7.805245193391D-12)*Y- &
! !           9.632707991704D-11)*Y+1.099047050624D-09)*Y- &
! !           1.15042731790748D-08 )*Y+1.09415155268932D-07 )*Y- &
! !           9.33687124875935D-07 )*Y+7.02338477986218D-06 )*Y- &
! !           4.53759748787756D-05 )*Y+2.41722511389146D-04 )*Y- &
! !           9.75935943447037D-04 )*Y+2.57520532789644D-03
! !       WW5 = ((((((((((((((( 7.28996979748849D-19*Y-1.26518146195173D-17) &
! !           *Y+1.886145834486D-16)*Y-2.876728287383D-15)*Y+ &
! !           4.114588668138D-14)*Y-5.44436631413933D-13 )*Y+ &
! !           6.64976446790959D-12 )*Y-7.44560069974940D-11 )*Y+ &
! !           7.57553198166848D-10 )*Y-6.92956101109829D-09 )*Y+ &
! !           5.62222859033624D-08 )*Y-3.97500114084351D-07 )*Y+ &
! !           2.39039126138140D-06 )*Y-1.18023950002105D-05 )*Y+ &
! !           4.52254031046244D-05 )*Y-1.21113782150370D-04 )*Y+ &
! !           1.75013126731224D-04
! !       RETURN

! ! !     XX=10.0 TO 15.0                             NROOTS = 5
! !   1600 Y = XX-12.5D+00
! !       RT1 = ((((((((((-4.16387977337393D-17*Y+7.20872997373860D-16)*Y+ & 
! !           1.395993802064D-14)*Y+3.660484641252D-14)*Y- &
! !           4.154857548139D-12)*Y+2.301379846544D-11)*Y- &
! !           1.033307012866D-09)*Y+3.997777641049D-08)*Y- &
! !           9.35118186333939D-07 )*Y+2.38589932752937D-05 )*Y- &
! !           5.35185183652937D-04 )*Y+8.85218988709735D-03
! !       RT2 = ((((((((((-4.56279214732217D-16*Y+6.24941647247927D-15)*Y+ &
! !           1.737896339191D-13)*Y+8.964205979517D-14)*Y- &
! !           3.538906780633D-11)*Y+9.561341254948D-11)*Y- &
! !           9.772831891310D-09)*Y+4.240340194620D-07)*Y- &
! !           1.02384302866534D-05 )*Y+2.57987709704822D-04 )*Y- &
! !           5.54735977651677D-03 )*Y+8.68245143991948D-02
! !       RT3 = ((((((((((-2.52879337929239D-15*Y+2.13925810087833D-14)*Y+ &
! !           7.884307667104D-13)*Y-9.023398159510D-13)*Y- &
! !           5.814101544957D-11)*Y-1.333480437968D-09)*Y- &
! !           2.217064940373D-08)*Y+1.643290788086D-06)*Y- &
! !           4.39602147345028D-05 )*Y+1.08648982748911D-03 )*Y- &
! !           2.13014521653498D-02 )*Y+2.94150684465425D-01
! !       RT4 = ((((((((((-6.42391438038888D-15*Y+5.37848223438815D-15)*Y+ &
! !           8.960828117859D-13)*Y+5.214153461337D-11)*Y- &
! !           1.106601744067D-10)*Y-2.007890743962D-08)*Y+ &
! !           1.543764346501D-07)*Y+4.520749076914D-06)*Y- &
! !           1.88893338587047D-04 )*Y+4.73264487389288D-03 )*Y- &
! !           7.91197893350253D-02 )*Y+8.60057928514554D-01
! !       RT5 = (((((((((((-2.24366166957225D-14*Y+4.87224967526081D-14)*Y+ &
! !           5.587369053655D-12)*Y-3.045253104617D-12)*Y- &
! !           1.223983883080D-09)*Y-2.05603889396319D-09 )*Y+ &
! !           2.58604071603561D-07 )*Y+1.34240904266268D-06 )*Y- &
! !           5.72877569731162D-05 )*Y-9.56275105032191D-04 )*Y+ &
! !           4.23367010370921D-02 )*Y-5.76800927133412D-01 )*Y+ &
! !           3.87328263873381D+00
! !       WW1 = ((((((((( 8.98007931950169D-15*Y+7.25673623859497D-14)*Y+ &
! !           5.851494250405D-14)*Y-4.234204823846D-11)*Y+ &
! !           3.911507312679D-10)*Y-9.65094802088511D-09 )*Y+ &
! !           3.42197444235714D-07 )*Y-7.51821178144509D-06 )*Y+ &
! !           1.94218051498662D-04 )*Y-5.38533819142287D-03 )*Y+ &
! !           1.68122596736809D-01
! !       WW2 = ((((((((((-1.05490525395105D-15*Y+1.96855386549388D-14)*Y- &
! !           5.500330153548D-13)*Y+1.003849567976D-11)*Y- &
! !           1.720997242621D-10)*Y+3.533277061402D-09)*Y- &
! !           6.389171736029D-08)*Y+1.046236652393D-06)*Y- &
! !           1.73148206795827D-05 )*Y+2.57820531617185D-04 )*Y- &
! !           3.46188265338350D-03 )*Y+7.03302497508176D-02
! !       WW3 = ((((((((((( 3.60020423754545D-16*Y-6.24245825017148D-15)*Y+ &
! !           9.945311467434D-14)*Y-1.749051512721D-12)*Y+ &
! !           2.768503957853D-11)*Y-4.08688551136506D-10 )*Y+ &
! !           6.04189063303610D-09 )*Y-8.23540111024147D-08 )*Y+ &
! !           1.01503783870262D-06 )*Y-1.20490761741576D-05 )*Y+ &
! !           1.26928442448148D-04 )*Y-1.05539461930597D-03 )*Y+ &
! !           1.15543698537013D-02
! !       WW4 = ((((((((((((( 2.51163533058925D-18*Y-4.31723745510697D-17)* &
! !           Y+6.557620865832D-16)*Y-1.016528519495D-14)*Y+ &
! !           1.491302084832D-13)*Y-2.06638666222265D-12 )*Y+ &
! !           2.67958697789258D-11 )*Y-3.23322654638336D-10 )*Y+ &
! !           3.63722952167779D-09 )*Y-3.75484943783021D-08 )*Y+ &
! !           3.49164261987184D-07 )*Y-2.92658670674908D-06 )*Y+ &
! !           2.12937256719543D-05 )*Y-1.19434130620929D-04 )*Y+ &
! !           6.45524336158384D-04
! !       WW5 = ((((((((((((((-1.29043630202811D-19*Y+2.16234952241296D-18)* &
! !           Y-3.107631557965D-17)*Y+4.570804313173D-16)*Y- &
! !           6.301348858104D-15)*Y+8.031304476153D-14)*Y- &
! !           9.446196472547D-13)*Y+1.018245804339D-11)*Y- &
! !           9.96995451348129D-11 )*Y+8.77489010276305D-10 )*Y- &
! !           6.84655877575364D-09 )*Y+4.64460857084983D-08 )*Y- & 
! !           2.66924538268397D-07 )*Y+1.24621276265907D-06 )*Y- &
! !           4.30868944351523D-06 )*Y+9.94307982432868D-06
! !       RETURN

! !   1800 IF (XX .GT. 25.0D+00) GO TO 2200
! !       IF (XX .GT. 20.0D+00) GO TO 2000
! ! !     XX=15.0 TO 20.0                             NROOTS = 5
! !       Y = XX-17.5D+00
! !       RT1 = (((((((((( 1.91875764545740D-16*Y+7.8357401095707D-16)*Y- &
! !           3.260875931644D-14)*Y-1.186752035569D-13)*Y+ & 
! !           4.275180095653D-12)*Y+3.357056136731D-11)*Y- &
! !           1.123776903884D-09)*Y+1.231203269887D-08)*Y- &
! !           3.99851421361031D-07 )*Y+1.45418822817771D-05 )*Y- &
! !           3.49912254976317D-04 )*Y+6.67768703938812D-03
! !       RT2 = (((((((((( 2.02778478673555D-15*Y+1.01640716785099D-14)*Y- &
! !           3.385363492036D-13)*Y-1.615655871159D-12)*Y+ &
! !           4.527419140333D-11)*Y+3.853670706486D-10)*Y- &
! !           1.184607130107D-08)*Y+1.347873288827D-07)*Y- &
! !           4.47788241748377D-06 )*Y+1.54942754358273D-04 )*Y- &
! !           3.55524254280266D-03 )*Y+6.44912219301603D-02
! !       RT3 = (((((((((( 7.79850771456444D-15*Y+6.00464406395001D-14)*Y- &
! !           1.249779730869D-12)*Y-1.020720636353D-11)*Y+ &
! !           1.814709816693D-10)*Y+1.766397336977D-09)*Y- &
! !           4.603559449010D-08)*Y+5.863956443581D-07)*Y- &
! !           2.03797212506691D-05 )*Y+6.31405161185185D-04 )*Y- &
! !           1.30102750145071D-02 )*Y+2.10244289044705D-01
! !       RT4 = (((((((((((-2.92397030777912D-15*Y+1.94152129078465D-14)*Y+ &
! !           4.859447665850D-13)*Y-3.217227223463D-12)*Y- &
! !           7.484522135512D-11)*Y+7.19101516047753D-10 )*Y+ &
! !           6.88409355245582D-09 )*Y-1.44374545515769D-07 )*Y+ &
! !           2.74941013315834D-06 )*Y-1.02790452049013D-04 )*Y+ &
! !           2.59924221372643D-03 )*Y-4.35712368303551D-02 )*Y+ &
! !           5.62170709585029D-01
! !       RT5 = ((((((((((( 1.17976126840060D-14*Y+1.24156229350669D-13)*Y- &
! !           3.892741622280D-12)*Y-7.755793199043D-12)*Y+ &
! !           9.492190032313D-10)*Y-4.98680128123353D-09 )*Y- &
! !           1.81502268782664D-07 )*Y+2.69463269394888D-06 )*Y+ &
! !           2.50032154421640D-05 )*Y-1.33684303917681D-03 )*Y+ &
! !           2.29121951862538D-02 )*Y-2.45653725061323D-01 )*Y+ &
! !           1.89999883453047D+00
! !       WW1 = (((((((((( 1.74841995087592D-15*Y-6.95671892641256D-16)*Y- &
! !           3.000659497257D-13)*Y+2.021279817961D-13)*Y+ &
! !           3.853596935400D-11)*Y+1.461418533652D-10)*Y- &
! !           1.014517563435D-08)*Y+1.132736008979D-07)*Y- &
! !           2.86605475073259D-06 )*Y+1.21958354908768D-04 )*Y- &
! !           3.86293751153466D-03 )*Y+1.45298342081522D-01
! !       WW2 = ((((((((((-1.11199320525573D-15*Y+1.85007587796671D-15)*Y+ &
! !           1.220613939709D-13)*Y+1.275068098526D-12)*Y- &
! !           5.341838883262D-11)*Y+6.161037256669D-10)*Y- &
! !           1.009147879750D-08)*Y+2.907862965346D-07)*Y- &
! !           6.12300038720919D-06 )*Y+1.00104454489518D-04 )*Y- &
! !           1.80677298502757D-03 )*Y+5.78009914536630D-02
! !       WW3 = ((((((((((-9.49816486853687D-16*Y+6.67922080354234D-15)*Y+ &
! !           2.606163540537D-15)*Y+1.983799950150D-12)*Y- &
! !           5.400548574357D-11)*Y+6.638043374114D-10)*Y- &
! !           8.799518866802D-09)*Y+1.791418482685D-07)*Y- &
! !           2.96075397351101D-06 )*Y+3.38028206156144D-05 )*Y- &
! !           3.58426847857878D-04 )*Y+8.39213709428516D-03
! !       WW4 = ((((((((((( 1.33829971060180D-17*Y-3.44841877844140D-16)*Y+ &
! !           4.745009557656D-15)*Y-6.033814209875D-14)*Y+ &
! !           1.049256040808D-12)*Y-1.70859789556117D-11 )*Y+ &
! !           2.15219425727959D-10 )*Y-2.52746574206884D-09 )*Y+ &
! !           3.27761714422960D-08 )*Y-3.90387662925193D-07 )*Y+ &
! !           3.46340204593870D-06 )*Y-2.43236345136782D-05 )*Y+ &
! !           3.54846978585226D-04
! !       WW5 = ((((((((((((( 2.69412277020887D-20*Y-4.24837886165685D-19)* &
! !           Y+6.030500065438D-18)*Y-9.069722758289D-17)*Y+ &
! !           1.246599177672D-15)*Y-1.56872999797549D-14 )*Y+ & 
! !           1.87305099552692D-13 )*Y-2.09498886675861D-12 )*Y+ &
! !           2.11630022068394D-11 )*Y-1.92566242323525D-10 )*Y+ &
! !           1.62012436344069D-09 )*Y-1.23621614171556D-08 )*Y+ &
! !           7.72165684563049D-08 )*Y-3.59858901591047D-07 )*Y+ &
! !           2.43682618601000D-06
! !       RETURN

! ! !     XX=20.0 TO 25.0                             NROOTS = 5
! !   2000 Y = XX-22.5D+00
! !       RT1 = (((((((((-1.13927848238726D-15*Y+7.39404133595713D-15)*Y+ &
! !           1.445982921243D-13)*Y-2.676703245252D-12)*Y+ &
! !           5.823521627177D-12)*Y+2.17264723874381D-10 )*Y+ &
! !           3.56242145897468D-09 )*Y-3.03763737404491D-07 )*Y+ &
! !           9.46859114120901D-06 )*Y-2.30896753853196D-04 )*Y+ &
! !           5.24663913001114D-03
! !       RT2 = (((((((((( 2.89872355524581D-16*Y-1.22296292045864D-14)*Y+ &
! !           6.184065097200D-14)*Y+1.649846591230D-12)*Y- &
! !           2.729713905266D-11)*Y+3.709913790650D-11)*Y+ &
! !           2.216486288382D-09)*Y+4.616160236414D-08)*Y- &
! !           3.32380270861364D-06 )*Y+9.84635072633776D-05 )*Y- &
! !           2.30092118015697D-03 )*Y+5.00845183695073D-02
! !       RT3 = (((((((((( 1.97068646590923D-15*Y-4.89419270626800D-14)*Y+ &
! !           1.136466605916D-13)*Y+7.546203883874D-12)*Y- &
! !           9.635646767455D-11)*Y-8.295965491209D-11)*Y+ &
! !           7.534109114453D-09)*Y+2.699970652707D-07)*Y- &
! !           1.42982334217081D-05 )*Y+3.78290946669264D-04 )*Y- &
! !           8.03133015084373D-03 )*Y+1.58689469640791D-01
! !       RT4 = (((((((((( 1.33642069941389D-14*Y-1.55850612605745D-13)*Y- &
! !           7.522712577474D-13)*Y+3.209520801187D-11)*Y- &
! !           2.075594313618D-10)*Y-2.070575894402D-09)*Y+ &
! !           7.323046997451D-09)*Y+1.851491550417D-06)*Y- &
! !           6.37524802411383D-05 )*Y+1.36795464918785D-03 )*Y- &
! !           2.42051126993146D-02 )*Y+3.97847167557815D-01
! !       RT5 = ((((((((((-6.07053986130526D-14*Y+1.04447493138843D-12)*Y- &
! !           4.286617818951D-13)*Y-2.632066100073D-10)*Y+ &
! !           4.804518986559D-09)*Y-1.835675889421D-08)*Y- &
! !           1.068175391334D-06)*Y+3.292234974141D-05)*Y- &
! !           5.94805357558251D-04 )*Y+8.29382168612791D-03 )*Y- &
! !           9.93122509049447D-02 )*Y+1.09857804755042D+00
! !       WW1 = (((((((((-9.10338640266542D-15*Y+1.00438927627833D-13)*Y+ &
! !           7.817349237071D-13)*Y-2.547619474232D-11)*Y+ &
! !           1.479321506529D-10)*Y+1.52314028857627D-09 )*Y+ &
! !           9.20072040917242D-09 )*Y-2.19427111221848D-06 )*Y+ &
! !           8.65797782880311D-05 )*Y-2.82718629312875D-03 )*Y+ &
! !           1.28718310443295D-01
! !       WW2 = ((((((((( 5.52380927618760D-15*Y-6.43424400204124D-14)*Y- &
! !           2.358734508092D-13)*Y+8.261326648131D-12)*Y+ &
! !           9.229645304956D-11)*Y-5.68108973828949D-09 )*Y+ &
! !           1.22477891136278D-07 )*Y-2.11919643127927D-06 )*Y+ &
! !           4.23605032368922D-05 )*Y-1.14423444576221D-03 )*Y+ &
! !           5.06607252890186D-02
! !       WW3 = ((((((((( 3.99457454087556D-15*Y-5.11826702824182D-14)*Y- &
! !           4.157593182747D-14)*Y+4.214670817758D-12)*Y+ &
! !           6.705582751532D-11)*Y-3.36086411698418D-09 )*Y+ &
! !           6.07453633298986D-08 )*Y-7.40736211041247D-07 )*Y+ &
! !           8.84176371665149D-06 )*Y-1.72559275066834D-04 )*Y+ &
! !           7.16639814253567D-03
! !       WW4 = (((((((((((-2.14649508112234D-18*Y-2.45525846412281D-18)*Y+ &
! !           6.126212599772D-16)*Y-8.526651626939D-15)*Y+ &
! !           4.826636065733D-14)*Y-3.39554163649740D-13 )*Y+ &
! !           1.67070784862985D-11 )*Y-4.42671979311163D-10 )*Y+ &
! !           6.77368055908400D-09 )*Y-7.03520999708859D-08 )*Y+ &
! !           6.04993294708874D-07 )*Y-7.80555094280483D-06 )*Y+ &
! !           2.85954806605017D-04
! !       WW5 = ((((((((((((-5.63938733073804D-21*Y+6.92182516324628D-20)*Y- &
! !           1.586937691507D-18)*Y+3.357639744582D-17)*Y- &
! !           4.810285046442D-16)*Y+5.386312669975D-15)*Y- &
! !           6.117895297439D-14)*Y+8.441808227634D-13)*Y- &
! !           1.18527596836592D-11 )*Y+1.36296870441445D-10 )*Y- &
! !           1.17842611094141D-09 )*Y+7.80430641995926D-09 )*Y- &
! !           5.97767417400540D-08 )*Y+1.65186146094969D-06
! !       RETURN

! !   2200 WW1 = SQRT(PIE4/XX)
! !       IF (XX .GT. 40.0D+00) GO TO 2400
! ! !     XX=25.0 TO 40.0                             NROOTS = 5
! !       E = EXP(-XX)
! !       RT1 = ((((((((-1.73363958895356D-06*XX+1.19921331441483D-04)*XX - &
! !           1.59437614121125D-02)*XX+1.13467897349442D+00)*XX - &
! !           4.47216460864586D+01)*XX+1.06251216612604D+03)*XX - &
! !           1.52073917378512D+04)*XX+1.20662887111273D+05)*XX - &
! !           4.07186366852475D+05)*E + R15/(XX-R15)
! !       RT2 = ((((((((-1.60102542621710D-05*XX+1.10331262112395D-03)*XX - &
! !           1.50043662589017D-01)*XX+1.05563640866077D+01)*XX - &
! !           4.10468817024806D+02)*XX+9.62604416506819D+03)*XX - &
! !           1.35888069838270D+05)*XX+1.06107577038340D+06)*XX - &
! !           3.51190792816119D+06)*E + R25/(XX-R25)
! !       RT3 = ((((((((-4.48880032128422D-05*XX+2.69025112122177D-03)*XX - &
! !           4.01048115525954D-01)*XX+2.78360021977405D+01)*XX - &
! !           1.04891729356965D+03)*XX+2.36985942687423D+04)*XX - &
! !           3.19504627257548D+05)*XX+2.34879693563358D+06)*XX - &
! !           7.16341568174085D+06)*E + R35/(XX-R35)
! !       RT4 = ((((((((-6.38526371092582D-05*XX-2.29263585792626D-03)*XX - &
! !           7.65735935499627D-02)*XX+9.12692349152792D+00)*XX - &
! !           2.32077034386717D+02)*XX+2.81839578728845D+02)*XX + &
! !           9.59529683876419D+04)*XX-1.77638956809518D+06)*XX + &
! !           1.02489759645410D+07)*E + R45/(XX-R45)
! !       RT5 = ((((((((-3.59049364231569D-05*XX-2.25963977930044D-02)*XX + &
! !           1.12594870794668D+00)*XX-4.56752462103909D+01)*XX + &
! !           1.05804526830637D+03)*XX-1.16003199605875D+04)*XX - &
! !           4.07297627297272D+04)*XX+2.22215528319857D+06)*XX - &
! !           1.61196455032613D+07)*E + R55/(XX-R55)
! !       WW5 = (((((((((-4.61100906133970D-10*XX+1.43069932644286D-07)*XX - &
! !           1.63960915431080D-05)*XX+1.15791154612838D-03)*XX - &
! !           5.30573476742071D-02)*XX+1.61156533367153D+00)*XX - &
! !           3.23248143316007D+01)*XX+4.12007318109157D+02)*XX - &
! !           3.02260070158372D+03)*XX+9.71575094154768D+03)*E + W55*WW1
! !       WW4 = (((((((((-2.40799435809950D-08*XX+8.12621667601546D-06)*XX - &
! !           9.04491430884113D-04)*XX+6.37686375770059D-02)*XX - &
! !           2.96135703135647D+00)*XX+9.15142356996330D+01)*XX - &
! !           1.86971865249111D+03)*XX+2.42945528916947D+04)*XX - &
! !           1.81852473229081D+05)*XX+5.96854758661427D+05)*E + W45*WW1
! !       WW3 = (((((((( 1.83574464457207D-05*XX-1.54837969489927D-03)*XX + &
! !           1.18520453711586D-01)*XX-6.69649981309161D+00)*XX + &
! !           2.44789386487321D+02)*XX-5.68832664556359D+03)*XX + &
! !           8.14507604229357D+04)*XX-6.55181056671474D+05)*XX + &
! !           2.26410896607237D+06)*E + W35*WW1
! !       WW2 = (((((((( 2.77778345870650D-05*XX-2.22835017655890D-03)*XX + &
! !           1.61077633475573D-01)*XX-8.96743743396132D+00)*XX + &
! !           3.28062687293374D+02)*XX-7.65722701219557D+03)*XX + &
! !           1.10255055017664D+05)*XX-8.92528122219324D+05)*XX + &
! !           3.10638627744347D+06)*E + W25*WW1
! !       WW1 = WW1-0.01962D+00*E-WW2-WW3-WW4-WW5
! !       RETURN

! !   2400 IF (XX .GT. 59.0D+00) GO TO 2600
! ! !     X=40.0 TO 59.0                             NROOTS = 5
! !       XXX = XX**3
! !       E = XXX*EXP(-XX)
! !       RT1 = (((-2.43758528330205D-02*XX+2.07301567989771D+00)*XX -6.45964225381113D+01)*XX+7.14160088655470D+02)*E + R15/(XX-R15)
! !       RT2 = (((-2.28861955413636D-01*XX+1.93190784733691D+01)*XX -5.99774730340912D+02)*XX+6.61844165304871D+03)*E + R25/(XX-R25)
! !       RT3 = (((-6.95053039285586D-01*XX+5.76874090316016D+01)*XX -1.77704143225520D+03)*XX+1.95366082947811D+04)*E + R35/(XX-R35)
! !       RT4 = (((-1.58072809087018D+00*XX+1.27050801091948D+02)*XX -3.86687350914280D+03)*XX+4.23024828121420D+04)*E + R45/(XX-R45)
! !       RT5 = (((-3.33963830405396D+00*XX+2.51830424600204D+02)*XX -7.57728527654961D+03)*XX+8.21966816595690D+04)*E + R55/(XX-R55)
! !       E = XXX*E
! !       WW5 = (( 1.35482430510942D-08*XX-3.27722199212781D-07)*XX +2.41522703684296D-06)*E + W55*WW1
! !       WW4 = (( 1.23464092261605D-06*XX-3.55224564275590D-05)*XX +3.03274662192286D-04)*E + W45*WW1
! !       WW3 = (( 1.34547929260279D-05*XX-4.19389884772726D-04)*XX +3.87706687610809D-03)*E + W35*WW1
! !       WW2 = (( 2.09539509123135D-05*XX-6.87646614786982D-04)*XX +6.68743788585688D-03)*E + W25*WW1
! !       WW1 = WW1-WW2-WW3-WW4-WW5
! !       RETURN

! ! !     XX=59.0 TO INFINITY                         NROOTS = 5
! !   2600 RT1 = R15/(XX-R15)
! !       RT2 = R25/(XX-R25)
! !       RT3 = R35/(XX-R35)
! !       RT4 = R45/(XX-R45)
! !       RT5 = R55/(XX-R55)
! !       WW2 = W25*WW1
! !       WW3 = W35*WW1
! !       WW4 = W45*WW1
! !       WW5 = W55*WW1
! !       WW1 = WW1-WW2-WW3-WW4-WW5
! !       !RETURN


!                    CALL ROOT5TS_GPU(XX,NROOTS,RT1,RT2,RT3,RT4,RT5,&
!       WW1,WW2,WW3,WW4,WW5)
!                    U(1)=RT1
!                    U(2)=RT2
!                    U(3)=RT3
!                    U(4)=RT4
!                    U(5)=RT5
!                    W(1)=WW1
!                    W(2)=WW2
!                    W(3)=WW3
!                    W(4)=WW4
!                    W(5)=WW5
             
!                MM = 0
!                MAX = NMAX+1

!                DO 420 M = 1,NROOTS
! ! C               DO 420 M = 1,5
!                   U2 = U(M)*RHO
!                   F00 = EXPE*W(M)
!                   DO I = 1,MAX
!                      IN(I) = IN1(I)+MM
!                   ENDDO

!                      DUMINV = ONE/(AB+U2*AANDB)
!                      DM2INV = HALF*DUMINV
!                      BP01 = (A+U2)*DM2INV
!                      B00 = U2*DM2INV
!                      B10 = (B+U2)*DM2INV
!                      XCP00 = (U2*C1X+C2X)*DUMINV
!                      XC00 = (U2*C3X+C4X)*DUMINV
!                      YCP00 = (U2*C1Y+C2Y)*DUMINV
!                      YC00 = (U2*C3Y+C4Y)*DUMINV
!                      ZCP00 = (U2*C1Z+C2Z)*DUMINV
!                      ZC00 = (U2*C3Z+C4Z)*DUMINV

!                   CALL XYZINT_gpu_2222(XIN,YIN,ZIN,&
!       IN,KN,NI,NJ,NK,NL,NMAX,MMAX,&
!       BP01,B00,B10,XCP00,XC00,YCP00,YC00,ZCP00,ZC00,F00,&
!       DXIJ,DYIJ,DZIJ,DXKL,DYKL,DZKL)

! ! ! !     ----- I(0,0) -----
! !       I1 = IN(1)
! !       XIN(I1) = ONE
! !       YIN(I1) = ONE
! !       ZIN(I1) = F00
! !       I2 = IN(2)
! !       K2 = KN(2)
! !       CP10 = B00
! !     ! ----- I(1,0) -----
! !       XIN(I2) = XC00
! !       YIN(I2) = YC00
! !       ZIN(I2) = ZC00*F00
! ! !     ----- I(0,1) -----
! !       I3 = I1+K2
! !       XIN(I3) = XCP00
! !       YIN(I3) = YCP00
! !       ZIN(I3) = ZCP00*F00
! ! !     ----- I(1,1) -----
! !       I3 = I2+K2
! !       XIN(I3) = XCP00*XIN(I2)+CP10
! !       YIN(I3) = YCP00*YIN(I2)+CP10
! !       ZIN(I3) = ZCP00*ZIN(I2)+CP10*F00
! !       C10 = ZERO
! !       I3 = I1
! !       I4 = I2
! !       DO N_NMAX = 2,NMAX
! !         C10 = C10+B10
! ! ! ----- I(N,0) -----
! !         I5 = IN(N_NMAX+1)

! !         XIN(I5) = C10*XIN(I3)+XC00*XIN(I4)
! !         YIN(I5) = C10*YIN(I3)+YC00*YIN(I4)
! !           ZIN(I5) = C10*ZIN(I3)+ZC00*ZIN(I4)

! !             CP10 = CP10+B00

! ! ! ! !     ----- I(N,1) -----

! !             I3 = I5+K2

! !             XIN(I3) = XCP00*XIN(I5)+CP10*XIN(I4)
! !             YIN(I3) = YCP00*YIN(I5)+CP10*YIN(I4)
! !             ZIN(I3) = ZCP00*ZIN(I5)+CP10*ZIN(I4)

! !           I3 = I4
! !           I4 = I5
! !         ENDDO

! !         CP01 = ZERO
! !         C01 = B00
! !         I3 = I1
! !         I4 = I1+K2
! !         DO M_MMAX = 2,MMAX
! !           CP01 = CP01+BP01

! ! ! !     ----- I(0,M) -----

! !           I5 = I1+KN(M_MMAX+1)

! !           XIN(I5) = CP01*XIN(I3)+XCP00*XIN(I4)
! !           YIN(I5) = CP01*YIN(I3)+YCP00*YIN(I4)
! !           ZIN(I5) = CP01*ZIN(I3)+ZCP00*ZIN(I4)
! ! ! !     ----- I(1,M) -----
! !             C01 = C01+B00
! !             I3 = I2+KN(M_MMAX+1)

! !             XIN(I3) = XC00*XIN(I5)+C01*XIN(I4)
! !             YIN(I3) = YC00*YIN(I5)+C01*YIN(I4)
! !             ZIN(I3) = ZC00*ZIN(I5)+C01*ZIN(I4)
! !           I3 = I4
! !           I4 = I5
! !         ENDDO

! ! !     ----- I(N,M) -----
! !         C01 = B00
! !         K3 = K2
! !         DO M_MMAX = 2,MMAX
! !           K4 = KN(M_MMAX+1)
! !           C01 = C01+B00
! !           I3 = I1
! !           I4 = I2
! !           C10 = B10
! !           DO N_NMAX = 2,NMAX
! !             I5 = IN(N_NMAX+1)

! !             XIN(I5+K4) = C10*XIN(I3+K4)+XC00*XIN(I4+K4)+C01*XIN(I4+K3)
! !             YIN(I5+K4) = C10*YIN(I3+K4)+YC00*YIN(I4+K4)+C01*YIN(I4+K3)
! !             ZIN(I5+K4) = C10*ZIN(I3+K4)+ZC00*ZIN(I4+K4)+C01*ZIN(I4+K3)

! !             C10 = C10+B10
! !             I3 = I4
! !             I4 = I5
! !           ENDDO
! !           K3 = K4
! !         ENDDO

! ! !     ----- I(NI,NJ,M) -----
! !         M_MMAX = 0
! !         I5 = IN(NMAX+1)
! !         FIRST1 = .TRUE.
! !         DO WHILE (FIRST1 .OR. M_MMAX .LE. MMAX)
! !           MIN = NI
! !           KM = KN(M+1)
! !           FIRST2 = .TRUE.
! !           DO WHILE (FIRST2 .OR. MIN .LT. NMAX)
! !             N_NMAX = NMAX
! !             I3 = I5+KM
! !             FIRST3 = .TRUE.
! !             DO WHILE (FIRST3 .OR. N_NMAX .GT. MIN)
! !               I4 = IN(N_NMAX)+KM

! !               XIN(I3) = XIN(I3)+DXIJ*XIN(I4)
! !               YIN(I3) = YIN(I3)+DYIJ*YIN(I4)
! !               ZIN(I3) = ZIN(I3)+DZIJ*ZIN(I4)

! !               I3 = I4
! !               N_NMAX = N_NMAX-1
! !               FIRST3 = .FALSE.
! !             ENDDO
! !             MIN = MIN+1
! !             FIRST2 = .FALSE.
! !           ENDDO
! !             I3 = 49+KM+I1
! !             DO NJ = 1,NJ
! !               I4 = I3
! !               DO NI = 1,NI

! !                 XIN(I4) = XIN(I4+294)+DXIJ*XIN(I4-49)
! !                 YIN(I4) = YIN(I4+294)+DYIJ*YIN(I4-49)
! !                 ZIN(I4) = ZIN(I4+294)+DZIJ*ZIN(I4-49)
! !                 I4 = I4+343
! !               ENDDO
! !               I3 = I3+49
! !             ENDDO
! !           M_MMAX = M_MMAX+1
! !           FIRST1 = .FALSE.
! !         ENDDO

! ! !     ----- I(NI,NJ,NK,NL) -----
! !         I5 = KN(MMAX+1)
! !         IA = I1
! !         NI = 0
! !         FIRST4 = .TRUE.
! !         DO 580 WHILE (FIRST4 .OR. NI .LE. NIMAX)
! !           NJ = 0
! !           IB = IA
! !           FIRST1 = .TRUE.
! !           DO 570 WHILE (FIRST1 .OR. NJ .LE. NJ)
! !             MIN = NK
! !             FIRST2 = .TRUE.
! !             DO 530 WHILE (FIRST2 .OR. MIN .LT. MMAX)
! !               M_MMAX = MMAX
! !               I3 = IB+I5
! !               FIRST3 = .TRUE.
! !               DO 520 WHILE (FIRST3 .OR. M_MMAX .GT. MIN)
! !                 I4 = IB+KN(M_MMAX)

! !                 XIN(I3) = XIN(I3)+DXKL*XIN(I4)
! !                 YIN(I3) = YIN(I3)+DYKL*YIN(I4)
! !                 ZIN(I3) = ZIN(I3)+DZKL*ZIN(I4)

! !                 I3 = I4
! !                 M_MMAX = M_MMAX-1
! !                 FIRST3 = .FALSE.
! !   520         END DO
! !               MIN = MIN+1
! !               FIRST2 = .FALSE.
! !   530       END DO

! !               I3 = IB+1
! !               DO NL = 1,NL
! !                 I4 = I3
! !                 DO NK = 1,NK

! !                   XIN(I4) = XIN(I4+6)+DXKL*XIN(I4-1)
! !                   YIN(I4) = YIN(I4+6)+DYKL*YIN(I4-1)
! !                   ZIN(I4) = ZIN(I4+6)+DZKL*ZIN(I4-1)

! !                   I4 = I4+7
! !                 ENDDO
! !               I3 = I3+1
! !               ENDDO
! !             NJ = NJ+1
! !             IB = IB+49
! !             FIRST1 = .FALSE.
! !   570     END DO
! !           NI = NI+1
! !           IA = IA+343
! !           FIRST4 = .FALSE.
! !   580   END DO





!                   MM = MM+2401
!   420          CONTINUE

!       DO I = 1,IJ
!       D1 = DIJ(I)
!       NX = IJX(I)
!       NY = IJY(I)
!       NZ = IJZ(I)
!       N1 = IJGT(I)
!       MAX = IK(I)
!       DO K = 1,MAX
!       MX = NX+KLX(K)
!       MY = NY+KLY(K)
!       MZ = NZ+KLZ(K)
!       NNN = N1+KLGT(K)
!       GHONDO(NNN) = GHONDO(NNN) + D1*DKL(K)*( XIN(MX)*YIN(MY)*ZIN(MZ)+ XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401) &      
!             + XIN(MX+ 4802)*YIN(MY+ 4802)*ZIN(MZ+ 4802)+ XIN(MX+ 7203)*YIN(MY+ 7203)*ZIN(MZ+ 7203)+XIN(MX+ 9604)*YIN(MY+ 9604)*ZIN(MZ+ 9604))
!       ENDDO
!       ENDDO
!   440       NN = NN+49
!   460    CONTINUE
!   480 CONTINUE








!       RETURN
!       END


      SUBROUTINE RT123TS_gpu(X,NROOTS,&
                            RT1,RT2,RT3,&
                            WW1,WW2,WW3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DATA R12,PIE4/2.75255128608411D-01, 7.85398163397448D-01/
      DATA R22,W22/ 2.72474487139158D+00, 9.17517095361369D-02/
      DATA R13/     1.90163509193487D-01/
      DATA R23,W23/ 1.78449274854325D+00, 1.77231492083829D-01/
      DATA R33,W33/ 5.52534374226326D+00, 5.11156880411248D-03/

!$omp declare target
      IF (X .GT. 5.0D+00) GO TO 400
      IF (X .GT. 1.0D+00) GO TO 280
      IF (X .GT. 3.0D-07) GO TO 180
!     X IS APPROXIMATELY ZERO.         NROOTS=1,2, OR 3
  !     IF (NROOTS-2) 120,140,160
  ! 120 RT1 = 0.5D+00 -X/5.0D+00
  !     WW1 = 1.0D+00 -X/3.0D+00
  !     RETURN
  140 RT1 = 1.30693606237085D-01 -2.90430236082028D-02 *X
      RT2 = 2.86930639376291D+00 -6.37623643058102D-01 *X
      WW1 = 6.52145154862545D-01 -1.22713621927067D-01 *X
      WW2 = 3.47854845137453D-01 -2.10619711404725D-01 *X
      RETURN
  ! 160 RT1 = 6.03769246832797D-02 -9.28875764357368D-03 *X
  !     RT2 = 7.76823355931043D-01 -1.19511285527878D-01 *X
  !     RT3 = 6.66279971938567D+00 -1.02504611068957D+00 *X
  !     WW1 = 4.67913934572691D-01 -5.64876917232519D-02 *X
  !     WW2 = 3.60761573048137D-01 -1.49077186455208D-01 *X
  !     WW3 = 1.71324492379169D-01 -1.27768455150979D-01 *X
  !     RETURN
!     X = 0.0 TO 1.0                   NROOTS=1,2, OR 3
  180 F1 = ((((((((-8.36313918003957D-08*X+1.21222603512827D-06 )*X-1.15662609053481D-05 )*X+9.25197374512647D-05 )*X-6.40994113129432D-04 )*X+3.78787044215009D-03 )*X-1.85185172458485D-02 )*X+7.14285713298222D-02 )*X-1.99999999997023D-01 )*X+3.33333333333318D-01
      WW1 = (X+X)*F1+EXP(-X)
      ! IF (NROOTS .EQ. 2) GO TO 200
      ! RT1 = F1/(WW1-F1)
      ! RETURN
  200 RT1 = (((((((-2.35234358048491D-09*X+2.49173650389842D-08)*X-4.558315364581D-08)*X-2.447252174587D-06)*X+4.743292959463D-05)*X-5.33184749432408D-04 )*X+ 4.44654947116579D-03 )*X-2.90430236084697D-02 )*X+1.30693606237085D-01
      RT2 = (((((((-2.47404902329170D-08*X+2.36809910635906D-07)*X+1.835367736310D-06)*X-2.066168802076D-05)*X-1.345693393936D-04)*X-5.88154362858038D-05 )*X+5.32735082098139D-02 )*X-6.37623643056745D-01 )*X+2.86930639376289D+00
      WW2 = ((F1-WW1)*RT1+F1)*(1.0D+00+RT2)/(RT2-RT1)
      WW1 = WW1-WW2
      RETURN
  260 T1 = RT1/(RT1+1.0D+00)
      T2 = RT2/(RT2+1.0D+00)
      T3 = RT3/(RT3+1.0D+00)
      A2 = F2-T1*F1
      A1 = F1-T1*WW1
      WW3 = (A2-T2*A1)/((T3-T2)*(T3-T1))
      WW2 = (T3*A1-A2)/((T3-T2)*(T2-T1))
      WW1 = WW1-WW2-WW3
      RETURN
  280 IF (X .GT. 3.0D+00) GO TO 340
!     X = 1.0 TO 3.0                   NROOTS=1,2, OR 3
      Y = X-2.0D+00
      F1 = ((((((((((-1.61702782425558D-10*Y+1.96215250865776D-09 )*Y-2.14234468198419D-08 )*Y+2.17216556336318D-07 )*Y-1.98850171329371D-06 )*Y+1.62429321438911D-05 )*Y-1.16740298039895D-04 )*Y+7.24888732052332D-04 )*Y-3.79490003707156D-03 )*Y+1.61723488664661D-02 )*Y-5.29428148329736D-02 )*Y+1.15702180856167D-01
      WW1 = (X+X)*F1+EXP(-X)
      ! IF (NROOTS .EQ. 2) GO TO 300
      ! RT1 = F1/(WW1-F1)
      ! RETURN
  300 RT1 = (((((((((-6.36859636616415D-12*Y+8.47417064776270D-11)*Y-5.152207846962D-10)*Y-3.846389873308D-10)*Y+8.472253388380D-08)*Y-1.85306035634293D-06 )*Y+2.47191693238413D-05 )*Y-2.49018321709815D-04 )*Y+2.19173220020161D-03 )*Y-1.63329339286794D-02 )*Y+8.68085688285261D-02
      RT2 = ((((((((( 1.45331350488343D-10*Y+2.07111465297976D-09)*Y-1.878920917404D-08)*Y-1.725838516261D-07)*Y+2.247389642339D-06)*Y+9.76783813082564D-06 )*Y-1.93160765581969D-04 )*Y-1.58064140671893D-03 )*Y+4.85928174507904D-02 )*Y-4.30761584997596D-01 )*Y+1.80400974537950D+00
      WW2 = ((F1-WW1)*RT1+F1)*(1.0D+00+RT2)/(RT2-RT1)
      WW1 = WW1-WW2
      RETURN
!     X = 3.0 TO 5.0                   NROOTS =1,2, OR 3
  340 Y = X-4.0D+00
      F1 = ((((((((((-2.62453564772299D-11*Y+3.24031041623823D-10 )*Y-3.614965656163D-09)*Y+3.760256799971D-08)*Y-3.553558319675D-07)*Y+3.022556449731D-06)*Y-2.290098979647D-05)*Y+1.526537461148D-04)*Y-8.81947375894379D-04 )*Y+4.33207949514611D-03 )*Y-1.75257821619926D-02 )*Y+5.28406320615584D-02
      WW1 = (X+X)*F1+EXP(-X)
  ! GO TO 360
      ! IF (NROOTS .EQ. 2) GO TO 360
      ! RT1 = F1/(WW1-F1)
      ! RETURN
  360 RT1 = ((((((((-4.11560117487296D-12*Y+7.10910223886747D-11)*Y-1.73508862390291D-09 )*Y+5.93066856324744D-08 )*Y- 9.76085576741771D-07 )*Y+1.08484384385679D-05 )*Y-1.12608004981982D-04 )*Y+1.16210907653515D-03 )*Y-9.89572595720351D-03 )*Y+6.12589701086408D-02
      RT2 = (((((((((-1.80555625241001D-10*Y+5.44072475994123D-10)*Y+1.603498045240D-08)*Y-1.497986283037D-07)*Y-7.017002532106D-07)*Y+1.85882653064034D-05 )*Y-2.04685420150802D-05 )*Y-2.49327728643089D-03 )*Y+3.56550690684281D-02 )*Y-2.60417417692375D-01 )*Y+ 1.12155283108289D+00
      WW2 = ((F1-WW1)*RT1+F1)*(1.0D+00+RT2)/(RT2-RT1)
      WW1 = WW1-WW2
      RETURN
  400 CONTINUE
  !400 IF (X .GT. 15.0D+00) GO TO 560
      E = EXP(-X)
      IF (X .GT. 10.0D+00) GO TO 480
!     X = 5.0 TO 10.0                  NROOTS =1,2, OR 3
      WW1 = (((((( 4.6897511375022D-01/X-6.9955602298985D-01)/X +5.3689283271887D-01)/X-3.2883030418398D-01)/X +2.4645596956002D-01)/X-4.9984072848436D-01)/X -3.1501078774085D-06)*E + SQRT(PIE4/X)
      F1 = (WW1-E)/(X+X)
      IF (NROOTS-2) 420,440,460
  420 RT1 = F1/(WW1-F1)
      RETURN
  440 Y = X-7.5D+00
      RT1 = (((((((((((((-1.43632730148572D-16*Y+2.38198922570405D-16)*Y+1.358319618800D-14)*Y-7.064522786879D-14)*Y-7.719300212748D-13)*Y+7.802544789997D-12)*Y+6.628721099436D-11)*Y-1.775564159743D-09)*Y+1.713828823990D-08)*Y-1.497500187053D-07)*Y+2.283485114279D-06)*Y-3.76953869614706D-05 )*Y+4.74791204651451D-04 )*Y-4.60448960876139D-03 )*Y+3.72458587837249D-02
      RT2 = (((((((((((( 2.48791622798900D-14*Y-1.36113510175724D-13)*Y-2.224334349799D-12)*Y+4.190559455515D-11)*Y-2.222722579924D-10)*Y-2.624183464275D-09)*Y+6.128153450169D-08)*Y-4.383376014528D-07)*Y-2.49952200232910D-06 )*Y+1.03236647888320D-04 )*Y-1.44614664924989D-03 )*Y+1.35094294917224D-02 )*Y-9.53478510453887D-02 )*Y+5.44765245686790D-01
      WW2 = ((F1-WW1)*RT1+F1)*(1.0D+00+RT2)/(RT2-RT1)
      WW1 = WW1-WW2
      RETURN
  460 F2 = (F1+F1+F1-E)/(X+X)
      Y = X-7.5D+00
      RT1 = ((((((((((( 5.74429401360115D-16*Y+7.11884203790984D-16)*Y-6.736701449826D-14)*Y-6.264613873998D-13)*Y+1.315418927040D-11)*Y-4.23879635610964D-11 )*Y+1.39032379769474D-09 )*Y-4.65449552856856D-08 )*Y+7.34609900170759D-07 )*Y-1.08656008854077D-05 )*Y+1.77930381549953D-04 )*Y-2.39864911618015D-03 )*Y+2.39112249488821D-02
      RT2 = ((((((((((( 1.13464096209120D-14*Y+6.99375313934242D-15)*Y-8.595618132088D-13)*Y-5.293620408757D-12)*Y-2.492175211635D-11)*Y+2.73681574882729D-09 )*Y-1.06656985608482D-08 )*Y-4.40252529648056D-07 )*Y+9.68100917793911D-06 )*Y-1.68211091755327D-04 )*Y+2.69443611274173D-03 )*Y-3.23845035189063D-02 )*Y+2.75969447451882D-01
      RT3 = (((((((((((( 6.66339416996191D-15*Y+1.84955640200794D-13)*Y-1.985141104444D-12)*Y-2.309293727603D-11)*Y+3.917984522103D-10)*Y+1.663165279876D-09)*Y-6.205591993923D-08)*Y+8.769581622041D-09)*Y+8.97224398620038D-06 )*Y-3.14232666170796D-05 )*Y-1.83917335649633D-03 )*Y+3.51246831672571D-02 )*Y-3.22335051270860D-01 )*Y+1.73582831755430D+00
      GO TO 260
!     X = 10.0 TO 15.0                 NROOTS=1,2, OR 3
  480 WW1 = (((-1.8784686463512D-01/X+2.2991849164985D-01)/X -4.9893752514047D-01)/X-2.1916512131607D-05)*E + SQRT(PIE4/X)
      F1 = (WW1-E)/(X+X)
      IF (NROOTS-2) 500,520,540
  500 RT1 = F1/(WW1-F1)
      RETURN
  520 RT1 = ((((-1.01041157064226D-05*X+1.19483054115173D-03)*X -6.73760231824074D-02)*X+1.25705571069895D+00)*X + (((-8.57609422987199D+03/X+5.91005939591842D+03)/X -1.70807677109425D+03)/X+2.64536689959503D+02)/X -2.38570496490846D+01)*E + R12/(X-R12)
      RT2 = ((( 3.39024225137123D-04*X-9.34976436343509D-02)*X -4.22216483306320D+00)*X + (((-2.08457050986847D+03/X -1.04999071905664D+03)/X+3.39891508992661D+02)/X -1.56184800325063D+02)/X+8.00839033297501D+00)*E + R22/(X-R22)
      WW2 = ((F1-WW1)*RT1+F1)*(1.0D+00+RT2)/(RT2-RT1)
      WW1 = WW1-WW2
      RETURN
  540 F2 = (F1+F1+F1-E)/(X+X)
      Y = X-12.5D+00
      RT1 = ((((((((((( 4.42133001283090D-16*Y-2.77189767070441D-15)*Y-4.084026087887D-14)*Y+5.379885121517D-13)*Y+1.882093066702D-12)*Y-8.67286219861085D-11 )*Y+7.11372337079797D-10 )*Y-3.55578027040563D-09 )*Y+1.29454702851936D-07 )*Y-4.14222202791434D-06 )*Y+8.04427643593792D-05 )*Y-1.18587782909876D-03 )*Y+1.53435577063174D-02
      RT2 = ((((((((((( 6.85146742119357D-15*Y-1.08257654410279D-14)*Y-8.579165965128D-13)*Y+6.642452485783D-12)*Y+4.798806828724D-11)*Y-1.13413908163831D-09 )*Y+7.08558457182751D-09 )*Y-5.59678576054633D-08 )*Y+2.51020389884249D-06 )*Y-6.63678914608681D-05 )*Y+1.11888323089714D-03 )*Y-1.45361636398178D-02 )*Y+1.65077877454402D-01
      RT3 = (((((((((((( 3.20622388697743D-15*Y-2.73458804864628D-14)*Y-3.157134329361D-13)*Y+8.654129268056D-12)*Y-5.625235879301D-11)*Y-7.718080513708D-10)*Y+2.064664199164D-08)*Y-1.567725007761D-07)*Y-1.57938204115055D-06 )*Y+6.27436306915967D-05 )*Y-1.01308723606946D-03 )*Y+1.13901881430697D-02 )*Y-1.01449652899450D-01 )*Y+7.77203937334739D-01
      GO TO 260
      !560 CONTINUE
  ! 560 IF (X .GT. 33.0D+00) GO TO 660
!     X = 15.0 TO 33.0                 NROOTS=1,2, OR 3
      E = EXP(-X)
      WW1 = (( 1.9623264149430D-01/X-4.9695241464490D-01)/X -6.0156581186481D-05)*E + SQRT(PIE4/X)
      F1 = (WW1-E)/(X+X)
      IF (NROOTS-2) 580,600,620
  580 RT1 = F1/(WW1-F1)
      RETURN
  600 RT1 = ((((-1.14906395546354D-06*X+1.76003409708332D-04)*X -1.71984023644904D-02)*X-1.37292644149838D-01)*X + (-4.75742064274859D+01/X+9.21005186542857D+00)/X -2.31080873898939D-02)*E + R12/(X-R12)
      RT2 = ((( 3.64921633404158D-04*X-9.71850973831558D-02)*X -4.02886174850252D+00)*X + (-1.35831002139173D+02/X -8.66891724287962D+01)/X+2.98011277766958D+00)*E + R22/(X-R22)
      WW2 = ((F1-WW1)*RT1+F1)*(1.0D+00+RT2)/(RT2-RT1)
      WW1 = WW1-WW2
      RETURN
  620 F2 = (F1+F1+F1-E)/(X+X)
      !IF (X .GT. 20.0D+00) GO TO 640
      RT1 = ((((((-2.43270989903742D-06*X+3.57901398988359D-04)*X -2.34112415981143D-02)*X+7.81425144913975D-01)*X -1.73209218219175D+01)*X+2.43517435690398D+02)*X + (-1.97611541576986D+04/X+9.82441363463929D+03)/X -2.07970687843258D+03)*E + R13/(X-R13)
      RT2 = (((((-2.62627010965435D-04*X+3.49187925428138D-02)*X -3.09337618731880D+00)*X+1.07037141010778D+02)*X -2.36659637247087D+03)*X + ((-2.91669113681020D+06/X +1.41129505262758D+06)/X-2.91532335433779D+05)/X +3.35202872835409D+04)*E + R23/(X-R23)
      RT3 = ((((( 9.31856404738601D-05*X-2.87029400759565D-02)*X -7.83503697918455D-01)*X-1.84338896480695D+01)*X +4.04996712650414D+02)*X + (-1.89829509315154D+05/X +5.11498390849158D+04)/X-6.88145821789955D+03)*E + R33/(X-R33)
      GO TO 260
  ! 640 RT1 = ((((-4.97561537069643D-04*X-5.00929599665316D-02)*X +1.31099142238996D+00)*X-1.88336409225481D+01)*X -6.60344754467191D+02 /X+1.64931462413877D+02)*E + R13/(X-R13)
  !     RT2 = ((((-4.48218898474906D-03*X-5.17373211334924D-01)*X +1.13691058739678D+01)*X-1.65426392885291D+02)*X -6.30909125686731D+03 /X+1.52231757709236D+03)*E + R23/(X-R23)
  !     RT3 = ((((-1.38368602394293D-02*X-1.77293428863008D+00)*X +1.73639054044562D+01)*X-3.57615122086961D+02)*X -1.45734701095912D+04 /X+2.69831813951849D+03)*E + R33/(X-R33)
  !     GO TO 260
!     X = 33.0 TO INFINITY             NROOTS=1,2, OR 3
  ! 660 WW1 = SQRT(PIE4/X)
  !     IF (NROOTS-2) 680,700,720
  ! 680 RT1 = 0.5D+00/(X-0.5D+00)
  !     RETURN
  ! ! 700 IF (X .GT. 40.0D+00) GO TO 740
  ! 700 E = EXP(-X)
  !     RT1 = (-8.78947307498880D-01*X+1.09243702330261D+01)*E + R12/(X-R12)
  !     RT2 = (-9.28903924275977D+00*X+8.10642367843811D+01)*E + R22/(X-R22)
  !     WW2 = ( 4.46857389308400D+00*X-7.79250653461045D+01)*E + W22*WW1
  !     WW1 = WW1-WW2
  !     RETURN
  ! 720 E = EXP(-X)
  !     RT1 = ((-7.39058467995275D+00*X+3.21318352526305D+02)*X -3.99433696473658D+03)*E + R13/(X-R13)
  !     RT2 = ((-7.38726243906513D+01*X+3.13569966333873D+03)*X -3.86862867311321D+04)*E + R23/(X-R23)
  !     RT3 = ((-2.63750565461336D+02*X+1.04412168692352D+04)*X -1.28094577915394D+05)*E + R33/(X-R33)
  !     WW3 = ((( 1.52258947224714D-01*X-8.30661900042651D+00)*X +1.92977367967984D+02)*X-1.67787926005344D+03)*E + W33*WW1
  !     WW2 = (( 6.15072615497811D+01*X-2.91980647450269D+03)*X +3.80794303087338D+04)*E + W23*WW1
  !     WW1 = WW1-WW2-WW3
  !     RETURN
  ! 740 RT1 = R12/(X-R12)
  !     RT2 = R22/(X-R22)
  !     WW2 = W22*WW1
  !     WW1 = WW1-WW2
  !     RETURN
  ! 760 RT1 = R13/(X-R13)
  !     RT2 = R23/(X-R23)
  !     RT3 = R33/(X-R33)
  !     WW2 = W23*WW1
  !     WW3 = W33*WW1
  !     WW1 = WW1-WW2-WW3
  !     RETURN
      END
      SUBROUTINE ROOT5TS_GPU(X,NROOTS,&
          RT1,RT2,RT3,RT4,RT5,&
          WW1,WW2,WW3,WW4,WW5)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
! C
! C          *****   VERSION  FEBRUARY 27,1975   *****
! C
      DATA R15,PIE4/1.17581320211778D-01, 7.85398163397448D-01/
      DATA R25,W25/ 1.07456201243690D+00, 2.70967405960535D-01/
      DATA R35,W35/ 3.08593744371754D+00, 3.82231610015404D-02/
      DATA R45,W45/ 6.41472973366203D+00, 1.51614186862443D-03/
      DATA R55,W55/ 1.18071894899717D+01, 8.62130526143657D-06/
! C
! C     X LARGE (ASYMPTOTIC SOLUTION)              NROOTS=5
! C
! CA    IF(X.GT.59.0D+00) THEN
! CA       CALL ROOTSA
! CA       RETURN
! CA    ENDIF
! C
!$omp declare target
      IF (X .GT. 15.0D+00) GO TO 180
      IF (X .GT. 5.0D+00) GO TO 140
      IF (X .GT. 1.0D+00) GO TO 120
      IF (X .GT. 3.0D-07) GO TO 100
!C     X IS APPROXIMATELY ZERO.                   NROOTS = 5
      RT1 = 2.26659266316985D-02 -2.15865967920897D-03 *X
      RT2 = 2.31271692140903D-01 -2.20258754389745D-02 *X
      RT3 = 8.57346024118836D-01 -8.16520023025515D-02 *X
      RT4 = 2.97353038120346D+00 -2.83193369647137D-01 *X
      RT5 = 1.84151859759051D+01 -1.75382723579439D+00 *X
      WW1 = 2.95524224714752D-01 -1.96867576909777D-02 *X
      WW2 = 2.69266719309995D-01 -5.61737590184721D-02 *X
      WW3 = 2.19086362515981D-01 -9.71152726793658D-02 *X
      WW4 = 1.49451349150580D-01 -1.02979262193565D-01 *X
      WW5 = 6.66713443086877D-02 -5.73782817488315D-02 *X
      RETURN

!C     X=0.0 TO 1.0                               NROOTS = 5
  100 RT1 = ((((((-4.46679165328413D-11*X+1.21879111988031D-09)*X- &
          2.62975022612104D-08 )*X+5.15106194905897D-07 )*X- &
          9.27933625824749D-06 )*X+1.51794097682482D-04 )*X- &
          2.15865967920301D-03 )*X+2.26659266316985D-02
      RT2 = (((((( 1.93117331714174D-10*X-4.57267589660699D-09)*X+ &
          2.48339908218932D-08 )*X+1.50716729438474D-06 )*X- &
          6.07268757707381D-05 )*X+1.37506939145643D-03 )*X- &
          2.20258754419939D-02 )*X+2.31271692140905D-01
      RT3 = ((((( 4.84989776180094D-09*X+1.31538893944284D-07)*X- &
          2.766753852879D-06)*X-7.651163510626D-05)*X+ &
          4.033058545972D-03)*X-8.16520022916145D-02 )*X+ &
          8.57346024118779D-01
      RT4 = ((((-2.48581772214623D-07*X-4.34482635782585D-06)*X- &
          7.46018257987630D-07 )*X+1.01210776517279D-02 )*X- &
          2.83193369640005D-01 )*X+2.97353038120345D+00
      RT5 = (((((-8.92432153868554D-09*X+1.77288899268988D-08)*X+ &
          3.040754680666D-06)*X+1.058229325071D-04)*X+ &
          4.596379534985D-02)*X-1.75382723579114D+00 )*X+ &
          1.84151859759049D+01
      WW1 = ((((((-2.03822632771791D-09*X+3.89110229133810D-08)*X- &
          5.84914787904823D-07 )*X+8.30316168666696D-06 )*X- &
          1.13218402310546D-04 )*X+1.49128888586790D-03 )*X- &
          1.96867576904816D-02 )*X+2.95524224714749D-01
      WW2 = ((((((( 8.62848118397570D-09*X-1.38975551148989D-07)*X+ &
          1.602894068228D-06)*X-1.646364300836D-05)*X+ &
          1.538445806778D-04)*X-1.28848868034502D-03 )*X+ &
          9.38866933338584D-03 )*X-5.61737590178812D-02 )*X+ &
          2.69266719309991D-01
      WW3 = ((((((((-9.41953204205665D-09*X+1.47452251067755D-07)*X- &
          1.57456991199322D-06 )*X+1.45098401798393D-05 )*X- &
          1.18858834181513D-04 )*X+8.53697675984210D-04 )*X- &
          5.22877807397165D-03 )*X+2.60854524809786D-02 )*X- &
          9.71152726809059D-02 )*X+2.19086362515979D-01
      WW4 = ((((((((-3.84961617022042D-08*X+5.66595396544470D-07)*X- &
          5.52351805403748D-06 )*X+4.53160377546073D-05 )*X- &
          3.22542784865557D-04 )*X+1.95682017370967D-03 )*X- &
          9.77232537679229D-03 )*X+3.79455945268632D-02 )*X- &
          1.02979262192227D-01 )*X+1.49451349150573D-01
      WW5 = ((((((((( 4.09594812521430D-09*X-6.47097874264417D-08)*X+ &
          6.743541482689D-07)*X-5.917993920224D-06)*X+ &
          4.531969237381D-05)*X-2.99102856679638D-04 )*X+ &
          1.65695765202643D-03 )*X-7.40671222520653D-03 )*X+ &
          2.50889946832192D-02 )*X-5.73782817487958D-02 )*X+ &
          6.66713443086877D-02
      RETURN

!C     XX=1.0 TO 5.0                               NROOTS = 5
  120 Y = X-3.0D+00
      RT1 = ((((((((-2.58163897135138D-14*Y+8.14127461488273D-13)*Y- &
          2.11414838976129D-11 )*Y+5.09822003260014D-10 )*Y- &
          1.16002134438663D-08 )*Y+2.46810694414540D-07 )*Y- &
          4.92556826124502D-06 )*Y+9.02580687971053D-05 )*Y- &
          1.45190025120726D-03 )*Y+1.73416786387475D-02
      RT2 = ((((((((( 1.04525287289788D-14*Y+5.44611782010773D-14)*Y- &
          4.831059411392D-12)*Y+1.136643908832D-10)*Y- &
          1.104373076913D-09)*Y-2.35346740649916D-08 )*Y+ &
          1.43772622028764D-06 )*Y-4.23405023015273D-05 )*Y+ &
          9.12034574793379D-04 )*Y-1.52479441718739D-02 )*Y+ &
          1.76055265928744D-01
      RT3 = (((((((((-6.89693150857911D-14*Y+5.92064260918861D-13)*Y+ &
          1.847170956043D-11)*Y-3.390752744265D-10)*Y- &
          2.995532064116D-09)*Y+1.57456141058535D-07 )*Y- &
          3.95859409711346D-07 )*Y-9.58924580919747D-05 )*Y+ &
          3.23551502557785D-03 )*Y-5.97587007636479D-02 )*Y+ &
          6.46432853383057D-01
      RT4 = ((((((((-3.61293809667763D-12*Y-2.70803518291085D-11)*Y+ &
          8.83758848468769D-10 )*Y+1.59166632851267D-08 )*Y- &
          1.32581997983422D-07 )*Y-7.60223407443995D-06 )*Y- &
          7.41019244900952D-05 )*Y+9.81432631743423D-03 )*Y- &
          2.23055570487771D-01 )*Y+2.21460798080643D+00
      RT5 = ((((((((( 7.12332088345321D-13*Y+3.16578501501894D-12)*Y- &
          8.776668218053D-11)*Y-2.342817613343D-09)*Y- &
          3.496962018025D-08)*Y-3.03172870136802D-07 )*Y+ &
          1.50511293969805D-06 )*Y+1.37704919387696D-04 )*Y+ &
          4.70723869619745D-02 )*Y-1.47486623003693D+00 )*Y+ &
          1.35704792175847D+01
      WW1 = ((((((((( 1.04348658616398D-13*Y-1.94147461891055D-12)*Y+ &
          3.485512360993D-11)*Y-6.277497362235D-10)*Y+ &
          1.100758247388D-08)*Y-1.88329804969573D-07 )*Y+ &
          3.12338120839468D-06 )*Y-5.04404167403568D-05 )*Y+ &
          8.00338056610995D-04 )*Y-1.30892406559521D-02 )*Y+ &
          2.47383140241103D-01
      WW2 = ((((((((((( 3.23496149760478D-14*Y-5.24314473469311D-13)*Y+ &
          7.743219385056D-12)*Y-1.146022750992D-10)*Y+ &
          1.615238462197D-09)*Y-2.15479017572233D-08 )*Y+ &
          2.70933462557631D-07 )*Y-3.18750295288531D-06 )*Y+ &
          3.47425221210099D-05 )*Y-3.45558237388223D-04 )*Y+ &
          3.05779768191621D-03 )*Y-2.29118251223003D-02 )*Y+ &
          1.59834227924213D-01
      WW3 = ((((((((((((-3.42790561802876D-14*Y+5.26475736681542D-13)*Y- &
          7.184330797139D-12)*Y+9.763932908544D-11)*Y- &
          1.244014559219D-09)*Y+1.472744068942D-08)*Y- &
          1.611749975234D-07)*Y+1.616487851917D-06)*Y- &
          1.46852359124154D-05 )*Y+1.18900349101069D-04 )*Y- &
          8.37562373221756D-04 )*Y+4.93752683045845D-03 )*Y- &
          2.25514728915673D-02 )*Y+6.95211812453929D-02
      WW4 = ((((((((((((( 1.04072340345039D-14*Y-1.60808044529211D-13)* &
          Y+2.183534866798D-12)*Y-2.939403008391D-11)*Y+ &
          3.679254029085D-10)*Y-4.23775673047899D-09 )*Y+ &
          4.46559231067006D-08 )*Y-4.26488836563267D-07 )*Y+ &
          3.64721335274973D-06 )*Y-2.74868382777722D-05 )*Y+ &
          1.78586118867488D-04 )*Y-9.68428981886534D-04 )*Y+ &
          4.16002324339929D-03 )*Y-1.28290192663141D-02 )*Y+ &
          2.22353727685016D-02
      WW5 = ((((((((((((((-8.16770412525963D-16*Y+1.31376515047977D-14)* &
          Y-1.856950818865D-13)*Y+2.596836515749D-12)*Y- &
          3.372639523006D-11)*Y+4.025371849467D-10)*Y- &
          4.389453269417D-09)*Y+4.332753856271D-08)*Y- &
          3.82673275931962D-07 )*Y+2.98006900751543D-06 )*Y- &
          2.00718990300052D-05 )*Y+1.13876001386361D-04 )*Y- &
          5.23627942443563D-04 )*Y+1.83524565118203D-03 )*Y- &
          4.37785737450783D-03 )*Y+5.36963805223095D-03
      RETURN

  140 IF (X .GT. 10.0D+00) GO TO 160
!C     XX=5.0 TO 10.0                              NROOTS = 5
      Y = X-7.5D+00
      RT1 = ((((((((-1.13825201010775D-14*Y+1.89737681670375D-13)*Y- &
          4.81561201185876D-12 )*Y+1.56666512163407D-10 )*Y- &
          3.73782213255083D-09 )*Y+9.15858355075147D-08 )*Y- &
          2.13775073585629D-06 )*Y+4.56547356365536D-05 )*Y- &
          8.68003909323740D-04 )*Y+1.22703754069176D-02
      RT2 = (((((((((-3.67160504428358D-15*Y+1.27876280158297D-14)*Y- &
          1.296476623788D-12)*Y+1.477175434354D-11)*Y+ &
          5.464102147892D-10)*Y-2.42538340602723D-08 )*Y+ &
          8.20460740637617D-07 )*Y-2.20379304598661D-05 )*Y+ &
          4.90295372978785D-04 )*Y-9.14294111576119D-03 )*Y+ &
          1.22590403403690D-01
      RT3 = ((((((((( 1.39017367502123D-14*Y-6.96391385426890D-13)*Y+ &
          1.176946020731D-12)*Y+1.725627235645D-10)*Y- &
          3.686383856300D-09)*Y+2.87495324207095D-08 )*Y+ &
          1.71307311000282D-06 )*Y-7.94273603184629D-05 )*Y+ &
          2.00938064965897D-03 )*Y-3.63329491677178D-02 )*Y+ &
          4.34393683888443D-01
      RT4 = ((((((((((-1.27815158195209D-14*Y+1.99910415869821D-14)*Y+ &
          3.753542914426D-12)*Y-2.708018219579D-11)*Y- &
          1.190574776587D-09)*Y+1.106696436509D-08)*Y+ &
          3.954955671326D-07)*Y-4.398596059588D-06)*Y- &
          2.01087998907735D-04 )*Y+7.89092425542937D-03 )*Y- &
          1.42056749162695D-01 )*Y+1.39964149420683D+00
      RT5 = ((((((((((-1.19442341030461D-13*Y-2.34074833275956D-12)*Y+ &
          6.861649627426D-12)*Y+6.082671496226D-10)*Y+ &
          5.381160105420D-09)*Y-6.253297138700D-08)*Y- &
          2.135966835050D-06)*Y-2.373394341886D-05)*Y+ &
          2.88711171412814D-06 )*Y+4.85221195290753D-02 )*Y- &
          1.04346091985269D+00 )*Y+7.89901551676692D+00
      WW1 = ((((((((( 7.95526040108997D-15*Y-2.48593096128045D-13)*Y+ &
          4.761246208720D-12)*Y-9.535763686605D-11)*Y+ &
          2.225273630974D-09)*Y-4.49796778054865D-08 )*Y+ &
          9.17812870287386D-07 )*Y-1.86764236490502D-05 )*Y+ &
          3.76807779068053D-04 )*Y-8.10456360143408D-03 )*Y+ &
          2.01097936411496D-01
      WW2 = ((((((((((( 1.25678686624734D-15*Y-2.34266248891173D-14)*Y+ & 
          3.973252415832D-13)*Y-6.830539401049D-12)*Y+ &
          1.140771033372D-10)*Y-1.82546185762009D-09 )*Y+ &
          2.77209637550134D-08 )*Y-4.01726946190383D-07 )*Y+ &
          5.48227244014763D-06 )*Y-6.95676245982121D-05 )*Y+ &
          8.05193921815776D-04 )*Y-8.15528438784469D-03 )*Y+ &
          9.71769901268114D-02
      WW3 = ((((((((((((-8.20929494859896D-16*Y+1.37356038393016D-14)*Y- &
          2.022863065220D-13)*Y+3.058055403795D-12)*Y- &
          4.387890955243D-11)*Y+5.923946274445D-10)*Y- &
          7.503659964159D-09)*Y+8.851599803902D-08)*Y- &
          9.65561998415038D-07 )*Y+9.60884622778092D-06 )*Y- &
          8.56551787594404D-05 )*Y+6.66057194311179D-04 )*Y- &
          4.17753183902198D-03 )*Y+2.25443826852447D-02
      WW4 = ((((((((((((((-1.08764612488790D-17*Y+1.85299909689937D-16)* &
          Y-2.730195628655D-15)*Y+4.127368817265D-14)*Y- &
          5.881379088074D-13)*Y+7.805245193391D-12)*Y- &
          9.632707991704D-11)*Y+1.099047050624D-09)*Y- &
          1.15042731790748D-08 )*Y+1.09415155268932D-07 )*Y- &
          9.33687124875935D-07 )*Y+7.02338477986218D-06 )*Y- &
          4.53759748787756D-05 )*Y+2.41722511389146D-04 )*Y- &
          9.75935943447037D-04 )*Y+2.57520532789644D-03
      WW5 = ((((((((((((((( 7.28996979748849D-19*Y-1.26518146195173D-17) &
          *Y+1.886145834486D-16)*Y-2.876728287383D-15)*Y+ &
          4.114588668138D-14)*Y-5.44436631413933D-13 )*Y+ &
          6.64976446790959D-12 )*Y-7.44560069974940D-11 )*Y+ &
          7.57553198166848D-10 )*Y-6.92956101109829D-09 )*Y+ &
          5.62222859033624D-08 )*Y-3.97500114084351D-07 )*Y+ &
          2.39039126138140D-06 )*Y-1.18023950002105D-05 )*Y+ &
          4.52254031046244D-05 )*Y-1.21113782150370D-04 )*Y+ &
          1.75013126731224D-04
      RETURN

!     XX=10.0 TO 15.0                             NROOTS = 5
  160 Y = X-12.5D+00
      RT1 = ((((((((((-4.16387977337393D-17*Y+7.20872997373860D-16)*Y+ & 
          1.395993802064D-14)*Y+3.660484641252D-14)*Y- &
          4.154857548139D-12)*Y+2.301379846544D-11)*Y- &
          1.033307012866D-09)*Y+3.997777641049D-08)*Y- &
          9.35118186333939D-07 )*Y+2.38589932752937D-05 )*Y- &
          5.35185183652937D-04 )*Y+8.85218988709735D-03
      RT2 = ((((((((((-4.56279214732217D-16*Y+6.24941647247927D-15)*Y+ &
          1.737896339191D-13)*Y+8.964205979517D-14)*Y- &
          3.538906780633D-11)*Y+9.561341254948D-11)*Y- &
          9.772831891310D-09)*Y+4.240340194620D-07)*Y- &
          1.02384302866534D-05 )*Y+2.57987709704822D-04 )*Y- &
          5.54735977651677D-03 )*Y+8.68245143991948D-02
      RT3 = ((((((((((-2.52879337929239D-15*Y+2.13925810087833D-14)*Y+ &
          7.884307667104D-13)*Y-9.023398159510D-13)*Y- &
          5.814101544957D-11)*Y-1.333480437968D-09)*Y- &
          2.217064940373D-08)*Y+1.643290788086D-06)*Y- &
          4.39602147345028D-05 )*Y+1.08648982748911D-03 )*Y- &
          2.13014521653498D-02 )*Y+2.94150684465425D-01
      RT4 = ((((((((((-6.42391438038888D-15*Y+5.37848223438815D-15)*Y+ &
          8.960828117859D-13)*Y+5.214153461337D-11)*Y- &
          1.106601744067D-10)*Y-2.007890743962D-08)*Y+ &
          1.543764346501D-07)*Y+4.520749076914D-06)*Y- &
          1.88893338587047D-04 )*Y+4.73264487389288D-03 )*Y- &
          7.91197893350253D-02 )*Y+8.60057928514554D-01
      RT5 = (((((((((((-2.24366166957225D-14*Y+4.87224967526081D-14)*Y+ &
          5.587369053655D-12)*Y-3.045253104617D-12)*Y- &
          1.223983883080D-09)*Y-2.05603889396319D-09 )*Y+ &
          2.58604071603561D-07 )*Y+1.34240904266268D-06 )*Y- &
          5.72877569731162D-05 )*Y-9.56275105032191D-04 )*Y+ &
          4.23367010370921D-02 )*Y-5.76800927133412D-01 )*Y+ &
          3.87328263873381D+00
      WW1 = ((((((((( 8.98007931950169D-15*Y+7.25673623859497D-14)*Y+ &
          5.851494250405D-14)*Y-4.234204823846D-11)*Y+ &
          3.911507312679D-10)*Y-9.65094802088511D-09 )*Y+ &
          3.42197444235714D-07 )*Y-7.51821178144509D-06 )*Y+ &
          1.94218051498662D-04 )*Y-5.38533819142287D-03 )*Y+ &
          1.68122596736809D-01
      WW2 = ((((((((((-1.05490525395105D-15*Y+1.96855386549388D-14)*Y- &
          5.500330153548D-13)*Y+1.003849567976D-11)*Y- &
          1.720997242621D-10)*Y+3.533277061402D-09)*Y- &
          6.389171736029D-08)*Y+1.046236652393D-06)*Y- &
          1.73148206795827D-05 )*Y+2.57820531617185D-04 )*Y- &
          3.46188265338350D-03 )*Y+7.03302497508176D-02
      WW3 = ((((((((((( 3.60020423754545D-16*Y-6.24245825017148D-15)*Y+ &
          9.945311467434D-14)*Y-1.749051512721D-12)*Y+ &
          2.768503957853D-11)*Y-4.08688551136506D-10 )*Y+ &
          6.04189063303610D-09 )*Y-8.23540111024147D-08 )*Y+ &
          1.01503783870262D-06 )*Y-1.20490761741576D-05 )*Y+ &
          1.26928442448148D-04 )*Y-1.05539461930597D-03 )*Y+ &
          1.15543698537013D-02
      WW4 = ((((((((((((( 2.51163533058925D-18*Y-4.31723745510697D-17)* &
          Y+6.557620865832D-16)*Y-1.016528519495D-14)*Y+ &
          1.491302084832D-13)*Y-2.06638666222265D-12 )*Y+ &
          2.67958697789258D-11 )*Y-3.23322654638336D-10 )*Y+ &
          3.63722952167779D-09 )*Y-3.75484943783021D-08 )*Y+ &
          3.49164261987184D-07 )*Y-2.92658670674908D-06 )*Y+ &
          2.12937256719543D-05 )*Y-1.19434130620929D-04 )*Y+ &
          6.45524336158384D-04
      WW5 = ((((((((((((((-1.29043630202811D-19*Y+2.16234952241296D-18)* &
          Y-3.107631557965D-17)*Y+4.570804313173D-16)*Y- &
          6.301348858104D-15)*Y+8.031304476153D-14)*Y- &
          9.446196472547D-13)*Y+1.018245804339D-11)*Y- &
          9.96995451348129D-11 )*Y+8.77489010276305D-10 )*Y- &
          6.84655877575364D-09 )*Y+4.64460857084983D-08 )*Y- & 
          2.66924538268397D-07 )*Y+1.24621276265907D-06 )*Y- &
          4.30868944351523D-06 )*Y+9.94307982432868D-06
      RETURN

  180 IF (X .GT. 25.0D+00) GO TO 220
      !IF (X .GT. 20.0D+00) GO TO 200
!     XX=15.0 TO 20.0                             NROOTS = 5
      Y = X-17.5D+00
      RT1 = (((((((((( 1.91875764545740D-16*Y+7.8357401095707D-16)*Y- &
          3.260875931644D-14)*Y-1.186752035569D-13)*Y+ & 
          4.275180095653D-12)*Y+3.357056136731D-11)*Y- &
          1.123776903884D-09)*Y+1.231203269887D-08)*Y- &
          3.99851421361031D-07 )*Y+1.45418822817771D-05 )*Y- &
          3.49912254976317D-04 )*Y+6.67768703938812D-03
      RT2 = (((((((((( 2.02778478673555D-15*Y+1.01640716785099D-14)*Y- &
          3.385363492036D-13)*Y-1.615655871159D-12)*Y+ &
          4.527419140333D-11)*Y+3.853670706486D-10)*Y- &
          1.184607130107D-08)*Y+1.347873288827D-07)*Y- &
          4.47788241748377D-06 )*Y+1.54942754358273D-04 )*Y- &
          3.55524254280266D-03 )*Y+6.44912219301603D-02
      RT3 = (((((((((( 7.79850771456444D-15*Y+6.00464406395001D-14)*Y- &
          1.249779730869D-12)*Y-1.020720636353D-11)*Y+ &
          1.814709816693D-10)*Y+1.766397336977D-09)*Y- &
          4.603559449010D-08)*Y+5.863956443581D-07)*Y- &
          2.03797212506691D-05 )*Y+6.31405161185185D-04 )*Y- &
          1.30102750145071D-02 )*Y+2.10244289044705D-01
      RT4 = (((((((((((-2.92397030777912D-15*Y+1.94152129078465D-14)*Y+ &
          4.859447665850D-13)*Y-3.217227223463D-12)*Y- &
          7.484522135512D-11)*Y+7.19101516047753D-10 )*Y+ &
          6.88409355245582D-09 )*Y-1.44374545515769D-07 )*Y+ &
          2.74941013315834D-06 )*Y-1.02790452049013D-04 )*Y+ &
          2.59924221372643D-03 )*Y-4.35712368303551D-02 )*Y+ &
          5.62170709585029D-01
      RT5 = ((((((((((( 1.17976126840060D-14*Y+1.24156229350669D-13)*Y- &
          3.892741622280D-12)*Y-7.755793199043D-12)*Y+ &
          9.492190032313D-10)*Y-4.98680128123353D-09 )*Y- &
          1.81502268782664D-07 )*Y+2.69463269394888D-06 )*Y+ &
          2.50032154421640D-05 )*Y-1.33684303917681D-03 )*Y+ &
          2.29121951862538D-02 )*Y-2.45653725061323D-01 )*Y+ &
          1.89999883453047D+00
      WW1 = (((((((((( 1.74841995087592D-15*Y-6.95671892641256D-16)*Y- &
          3.000659497257D-13)*Y+2.021279817961D-13)*Y+ &
          3.853596935400D-11)*Y+1.461418533652D-10)*Y- &
          1.014517563435D-08)*Y+1.132736008979D-07)*Y- &
          2.86605475073259D-06 )*Y+1.21958354908768D-04 )*Y- &
          3.86293751153466D-03 )*Y+1.45298342081522D-01
      WW2 = ((((((((((-1.11199320525573D-15*Y+1.85007587796671D-15)*Y+ &
          1.220613939709D-13)*Y+1.275068098526D-12)*Y- &
          5.341838883262D-11)*Y+6.161037256669D-10)*Y- &
          1.009147879750D-08)*Y+2.907862965346D-07)*Y- &
          6.12300038720919D-06 )*Y+1.00104454489518D-04 )*Y- &
          1.80677298502757D-03 )*Y+5.78009914536630D-02
      WW3 = ((((((((((-9.49816486853687D-16*Y+6.67922080354234D-15)*Y+ &
          2.606163540537D-15)*Y+1.983799950150D-12)*Y- &
          5.400548574357D-11)*Y+6.638043374114D-10)*Y- &
          8.799518866802D-09)*Y+1.791418482685D-07)*Y- &
          2.96075397351101D-06 )*Y+3.38028206156144D-05 )*Y- &
          3.58426847857878D-04 )*Y+8.39213709428516D-03
      WW4 = ((((((((((( 1.33829971060180D-17*Y-3.44841877844140D-16)*Y+ &
          4.745009557656D-15)*Y-6.033814209875D-14)*Y+ &
          1.049256040808D-12)*Y-1.70859789556117D-11 )*Y+ &
          2.15219425727959D-10 )*Y-2.52746574206884D-09 )*Y+ &
          3.27761714422960D-08 )*Y-3.90387662925193D-07 )*Y+ &
          3.46340204593870D-06 )*Y-2.43236345136782D-05 )*Y+ &
          3.54846978585226D-04
      WW5 = ((((((((((((( 2.69412277020887D-20*Y-4.24837886165685D-19)* &
          Y+6.030500065438D-18)*Y-9.069722758289D-17)*Y+ &
          1.246599177672D-15)*Y-1.56872999797549D-14 )*Y+ & 
          1.87305099552692D-13 )*Y-2.09498886675861D-12 )*Y+ &
          2.11630022068394D-11 )*Y-1.92566242323525D-10 )*Y+ &
          1.62012436344069D-09 )*Y-1.23621614171556D-08 )*Y+ &
          7.72165684563049D-08 )*Y-3.59858901591047D-07 )*Y+ &
          2.43682618601000D-06
!       RETURN

! !     XX=20.0 TO 25.0                             NROOTS = 5
!   200 Y = X-22.5D+00
!       RT1 = (((((((((-1.13927848238726D-15*Y+7.39404133595713D-15)*Y+ &
!           1.445982921243D-13)*Y-2.676703245252D-12)*Y+ &
!           5.823521627177D-12)*Y+2.17264723874381D-10 )*Y+ &
!           3.56242145897468D-09 )*Y-3.03763737404491D-07 )*Y+ &
!           9.46859114120901D-06 )*Y-2.30896753853196D-04 )*Y+ &
!           5.24663913001114D-03
!       RT2 = (((((((((( 2.89872355524581D-16*Y-1.22296292045864D-14)*Y+ &
!           6.184065097200D-14)*Y+1.649846591230D-12)*Y- &
!           2.729713905266D-11)*Y+3.709913790650D-11)*Y+ &
!           2.216486288382D-09)*Y+4.616160236414D-08)*Y- &
!           3.32380270861364D-06 )*Y+9.84635072633776D-05 )*Y- &
!           2.30092118015697D-03 )*Y+5.00845183695073D-02
!       RT3 = (((((((((( 1.97068646590923D-15*Y-4.89419270626800D-14)*Y+ &
!           1.136466605916D-13)*Y+7.546203883874D-12)*Y- &
!           9.635646767455D-11)*Y-8.295965491209D-11)*Y+ &
!           7.534109114453D-09)*Y+2.699970652707D-07)*Y- &
!           1.42982334217081D-05 )*Y+3.78290946669264D-04 )*Y- &
!           8.03133015084373D-03 )*Y+1.58689469640791D-01
!       RT4 = (((((((((( 1.33642069941389D-14*Y-1.55850612605745D-13)*Y- &
!           7.522712577474D-13)*Y+3.209520801187D-11)*Y- &
!           2.075594313618D-10)*Y-2.070575894402D-09)*Y+ &
!           7.323046997451D-09)*Y+1.851491550417D-06)*Y- &
!           6.37524802411383D-05 )*Y+1.36795464918785D-03 )*Y- &
!           2.42051126993146D-02 )*Y+3.97847167557815D-01
!       RT5 = ((((((((((-6.07053986130526D-14*Y+1.04447493138843D-12)*Y- &
!           4.286617818951D-13)*Y-2.632066100073D-10)*Y+ &
!           4.804518986559D-09)*Y-1.835675889421D-08)*Y- &
!           1.068175391334D-06)*Y+3.292234974141D-05)*Y- &
!           5.94805357558251D-04 )*Y+8.29382168612791D-03 )*Y- &
!           9.93122509049447D-02 )*Y+1.09857804755042D+00
!       WW1 = (((((((((-9.10338640266542D-15*Y+1.00438927627833D-13)*Y+ &
!           7.817349237071D-13)*Y-2.547619474232D-11)*Y+ &
!           1.479321506529D-10)*Y+1.52314028857627D-09 )*Y+ &
!           9.20072040917242D-09 )*Y-2.19427111221848D-06 )*Y+ &
!           8.65797782880311D-05 )*Y-2.82718629312875D-03 )*Y+ &
!           1.28718310443295D-01
!       WW2 = ((((((((( 5.52380927618760D-15*Y-6.43424400204124D-14)*Y- &
!           2.358734508092D-13)*Y+8.261326648131D-12)*Y+ &
!           9.229645304956D-11)*Y-5.68108973828949D-09 )*Y+ &
!           1.22477891136278D-07 )*Y-2.11919643127927D-06 )*Y+ &
!           4.23605032368922D-05 )*Y-1.14423444576221D-03 )*Y+ &
!           5.06607252890186D-02
!       WW3 = ((((((((( 3.99457454087556D-15*Y-5.11826702824182D-14)*Y- &
!           4.157593182747D-14)*Y+4.214670817758D-12)*Y+ &
!           6.705582751532D-11)*Y-3.36086411698418D-09 )*Y+ &
!           6.07453633298986D-08 )*Y-7.40736211041247D-07 )*Y+ &
!           8.84176371665149D-06 )*Y-1.72559275066834D-04 )*Y+ &
!           7.16639814253567D-03
!       WW4 = (((((((((((-2.14649508112234D-18*Y-2.45525846412281D-18)*Y+ &
!           6.126212599772D-16)*Y-8.526651626939D-15)*Y+ &
!           4.826636065733D-14)*Y-3.39554163649740D-13 )*Y+ &
!           1.67070784862985D-11 )*Y-4.42671979311163D-10 )*Y+ &
!           6.77368055908400D-09 )*Y-7.03520999708859D-08 )*Y+ &
!           6.04993294708874D-07 )*Y-7.80555094280483D-06 )*Y+ &
!           2.85954806605017D-04
!       WW5 = ((((((((((((-5.63938733073804D-21*Y+6.92182516324628D-20)*Y- &
!           1.586937691507D-18)*Y+3.357639744582D-17)*Y- &
!           4.810285046442D-16)*Y+5.386312669975D-15)*Y- &
!           6.117895297439D-14)*Y+8.441808227634D-13)*Y- &
!           1.18527596836592D-11 )*Y+1.36296870441445D-10 )*Y- &
!           1.17842611094141D-09 )*Y+7.80430641995926D-09 )*Y- &
!           5.97767417400540D-08 )*Y+1.65186146094969D-06
!       RETURN
  220 WW1 = SQRT(PIE4/X)
      write(*,*) "PIE4",PIE4
      write(*,*) "X",X
      !IF (X .GT. 40.0D+00) GO TO 240
!     XX=25.0 TO 40.0                             NROOTS = 5
      E = EXP(-X)
      RT1 = ((((((((-1.73363958895356D-06*X+1.19921331441483D-04)*X - &
          1.59437614121125D-02)*X+1.13467897349442D+00)*X - &
          4.47216460864586D+01)*X+1.06251216612604D+03)*X - &
          1.52073917378512D+04)*X+1.20662887111273D+05)*X - &
          4.07186366852475D+05)*E + R15/(X-R15)
      RT2 = ((((((((-1.60102542621710D-05*X+1.10331262112395D-03)*X - &
          1.50043662589017D-01)*X+1.05563640866077D+01)*X - &
          4.10468817024806D+02)*X+9.62604416506819D+03)*X - &
          1.35888069838270D+05)*X+1.06107577038340D+06)*X - &
          3.51190792816119D+06)*E + R25/(X-R25)
      RT3 = ((((((((-4.48880032128422D-05*X+2.69025112122177D-03)*X - &
          4.01048115525954D-01)*X+2.78360021977405D+01)*X - &
          1.04891729356965D+03)*X+2.36985942687423D+04)*X - &
          3.19504627257548D+05)*X+2.34879693563358D+06)*X - &
          7.16341568174085D+06)*E + R35/(X-R35)
      RT4 = ((((((((-6.38526371092582D-05*X-2.29263585792626D-03)*X - &
          7.65735935499627D-02)*X+9.12692349152792D+00)*X - &
          2.32077034386717D+02)*X+2.81839578728845D+02)*X + &
          9.59529683876419D+04)*X-1.77638956809518D+06)*X + &
          1.02489759645410D+07)*E + R45/(X-R45)
      RT5 = ((((((((-3.59049364231569D-05*X-2.25963977930044D-02)*X + &
          1.12594870794668D+00)*X-4.56752462103909D+01)*X + &
          1.05804526830637D+03)*X-1.16003199605875D+04)*X - &
          4.07297627297272D+04)*X+2.22215528319857D+06)*X - &
          1.61196455032613D+07)*E + R55/(X-R55)
      WW5 = (((((((((-4.61100906133970D-10*X+1.43069932644286D-07)*X - &
          1.63960915431080D-05)*X+1.15791154612838D-03)*X - &
          5.30573476742071D-02)*X+1.61156533367153D+00)*X - &
          3.23248143316007D+01)*X+4.12007318109157D+02)*X - &
          3.02260070158372D+03)*X+9.71575094154768D+03)*E + W55*WW1
      WW4 = (((((((((-2.40799435809950D-08*X+8.12621667601546D-06)*X - &
          9.04491430884113D-04)*X+6.37686375770059D-02)*X - &
          2.96135703135647D+00)*X+9.15142356996330D+01)*X - &
          1.86971865249111D+03)*X+2.42945528916947D+04)*X - &
          1.81852473229081D+05)*X+5.96854758661427D+05)*E + W45*WW1
      WW3 = (((((((( 1.83574464457207D-05*X-1.54837969489927D-03)*X + &
          1.18520453711586D-01)*X-6.69649981309161D+00)*X + &
          2.44789386487321D+02)*X-5.68832664556359D+03)*X + &
          8.14507604229357D+04)*X-6.55181056671474D+05)*X + &
          2.26410896607237D+06)*E + W35*WW1
      WW2 = (((((((( 2.77778345870650D-05*X-2.22835017655890D-03)*X + &
          1.61077633475573D-01)*X-8.96743743396132D+00)*X + &
          3.28062687293374D+02)*X-7.65722701219557D+03)*X + &
          1.10255055017664D+05)*X-8.92528122219324D+05)*X + &
          3.10638627744347D+06)*E + W25*WW1
      WW1 = WW1-0.01962D+00*E-WW2-WW3-WW4-WW5
!       RETURN

!   240 IF (X .GT. 59.0D+00) GO TO 260
! !     X=40.0 TO 59.0                             NROOTS = 5
!       XXX = X**3
!       E = XXX*EXP(-X)
!       RT1 = (((-2.43758528330205D-02*X+2.07301567989771D+00)*X -6.45964225381113D+01)*X+7.14160088655470D+02)*E + R15/(X-R15)
!       RT2 = (((-2.28861955413636D-01*X+1.93190784733691D+01)*X -5.99774730340912D+02)*X+6.61844165304871D+03)*E + R25/(X-R25)
!       RT3 = (((-6.95053039285586D-01*X+5.76874090316016D+01)*X -1.77704143225520D+03)*X+1.95366082947811D+04)*E + R35/(X-R35)
!       RT4 = (((-1.58072809087018D+00*X+1.27050801091948D+02)*X -3.86687350914280D+03)*X+4.23024828121420D+04)*E + R45/(X-R45)
!       RT5 = (((-3.33963830405396D+00*X+2.51830424600204D+02)*X -7.57728527654961D+03)*X+8.21966816595690D+04)*E + R55/(X-R55)
!       E = XXX*E
!       WW5 = (( 1.35482430510942D-08*X-3.27722199212781D-07)*X +2.41522703684296D-06)*E + W55*WW1
!       WW4 = (( 1.23464092261605D-06*X-3.55224564275590D-05)*X +3.03274662192286D-04)*E + W45*WW1
!       WW3 = (( 1.34547929260279D-05*X-4.19389884772726D-04)*X +3.87706687610809D-03)*E + W35*WW1
!       WW2 = (( 2.09539509123135D-05*X-6.87646614786982D-04)*X +6.68743788585688D-03)*E + W25*WW1
!       WW1 = WW1-WW2-WW3-WW4-WW5
!       RETURN

! !     XX=59.0 TO INFINITY                         NROOTS = 5
!   260 RT1 = R15/(X-R15)
!       RT2 = R25/(X-R25)
!       RT3 = R35/(X-R35)
!       RT4 = R45/(X-R45)
!       RT5 = R55/(X-R55)
!       WW2 = W25*WW1
!       WW3 = W35*WW1
!       WW4 = W45*WW1
!       WW5 = W55*WW1
!       WW1 = WW1-WW2-WW3-WW4-WW5
      RETURN
      END

      SUBROUTINE XYZINT_gpu_2(XINT,YINT,ZINT,&
      I,K,NIMAX,NJMAX,NKMAX,NLMAX,NMAX,MMAX,&
      BP01,B00,B10,XCP00,XC00,YCP00,YC00,ZCP00,ZC00,F00,&
      DXIJ,DYIJ,DZIJ,DXKL,DYKL,DZKL)


      !use integral_tmp_mod, only: XINT
      use omp_lib

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      LOGICAL N0,N1,M0,M1,FIRST1,FIRST2,FIRST3,FIRST4
      DIMENSION :: I(13),K(13)
      INTEGER :: NIMAX,NJMAX,NKMAX,NLMAX,NMAX,MMAX
      double precision :: BP01,B00,B10,XCP00,XC00,YCP00,YC00,ZCP00,ZC00
      double precision :: F00,DXIJ,DYIJ,DZIJ,DXKL,DYKL,DZKL
      double precision :: XINT(31213),YINT(31213),ZINT(31213)

      PARAMETER (ZERO=0.0D+00)
      PARAMETER (ONE=1.0D+00)
!$omp declare target

!ID=omp_get_thread_num()

      N0 = NMAX .EQ. 0
      N1 = NMAX .LE. 1
      M0 = MMAX .EQ. 0
      M1 = MMAX .LE. 1
! !     ----- I(0,0) -----
      I1 = I(1)

      XINT(I1) = ONE
      YINT(I1) = ONE
      ZINT(I1) = F00
      IF (N0 .AND. M0) RETURN
      I2 = I(2)
      K2 = K(2)
      CP10 = B00

    ! ----- I(1,0) -----
      IF (.NOT. N0) THEN
        XINT(I2) = XC00
        YINT(I2) = YC00
        ZINT(I2) = ZC00*F00
        IF (M0) GO TO 120
      END IF

! ! !     ----- I(0,1) -----
      I3 = I1+K2

      XINT(I3) = XCP00
      YINT(I3) = YCP00
      ZINT(I3) = ZCP00*F00

! ! ! ! !     ----- I(1,1) -----

      IF (.NOT. N0) THEN
        I3 = I2+K2

        XINT(I3) = XCP00*XINT(I2)+CP10
        YINT(I3) = YCP00*YINT(I2)+CP10
        ZINT(I3) = ZCP00*ZINT(I2)+CP10*F00
        !write(*,*) "ZINT(I3)", ZINT(I3)

        !if (varx(I3,ID) .NE. XINT(I3)) write(*,*) "ugh33"

      END IF

  120 CONTINUE

      IF (.NOT. N1) THEN
        C10 = ZERO
        I3 = I1
        I4 = I2
        DO 160 N = 2,NMAX
          C10 = C10+B10

! ! !     ----- I(N,0) -----

          I5 = I(N+1)

          XINT(I5) = C10*XINT(I3)+XC00*XINT(I4)
          YINT(I5) = C10*YINT(I3)+YC00*YINT(I4)
          ZINT(I5) = C10*ZINT(I3)+ZC00*ZINT(I4)
          IF ( .NOT. M0) THEN
            CP10 = CP10+B00

! ! !     ----- I(N,1) -----

            I3 = I5+K2

            XINT(I3) = XCP00*XINT(I5)+CP10*XINT(I4)
            YINT(I3) = YCP00*YINT(I5)+CP10*YINT(I4)
            ZINT(I3) = ZCP00*ZINT(I5)+CP10*ZINT(I4)

          END IF
          I3 = I4
          I4 = I5
  160     CONTINUE
      END IF

      IF ( .NOT. M1) THEN
        CP01 = ZERO
        C01 = B00
        I3 = I1
        I4 = I1+K2
        DO 220 M = 2,MMAX
          CP01 = CP01+BP01

! !     ----- I(0,M) -----

          I5 = I1+K(M+1)

          XINT(I5) = CP01*XINT(I3)+XCP00*XINT(I4)
          YINT(I5) = CP01*YINT(I3)+YCP00*YINT(I4)
          ZINT(I5) = CP01*ZINT(I3)+ZCP00*ZINT(I4)
! !     ----- I(1,M) -----

          IF (.NOT. N0) THEN
            C01 = C01+B00
            I3 = I2+K(M+1)

            XINT(I3) = XC00*XINT(I5)+C01*XINT(I4)
            YINT(I3) = YC00*YINT(I5)+C01*YINT(I4)
            ZINT(I3) = ZC00*ZINT(I5)+C01*ZINT(I4)
          END IF
          I3 = I4
          I4 = I5
  220   CONTINUE
      END IF

!     ----- I(N,M) -----

      IF (.NOT. N1 .AND. .NOT. M1) THEN
        C01 = B00
        K3 = K2
        DO 280 M = 2,MMAX
          K4 = K(M+1)
          C01 = C01+B00
          I3 = I1
          I4 = I2
          C10 = B10
          DO 260 N = 2,NMAX
            I5 = I(N+1)

            XINT(I5+K4) = C10*XINT(I3+K4)+XC00*XINT(I4+K4)+C01*XINT(I4+K3)
            YINT(I5+K4) = C10*YINT(I3+K4)+YC00*YINT(I4+K4)+C01*YINT(I4+K3)
            ZINT(I5+K4) = C10*ZINT(I3+K4)+ZC00*ZINT(I4+K4)+C01*ZINT(I4+K3)

            C10 = C10+B10
            I3 = I4
            I4 = I5
  260     CONTINUE
          K3 = K4
  280   CONTINUE
      END IF

!     ----- I(NI,NJ,M) -----

      IF (NJMAX .GT. 0) THEN
        M = 0
        I5 = I(NMAX+1)
        FIRST1 = .TRUE.
        DO 430 WHILE (FIRST1 .OR. M .LE. MMAX)
          MIN = NIMAX
          KM = K(M+1)
          FIRST2 = .TRUE.
          DO 360 WHILE (FIRST2 .OR. MIN .LT. NMAX)
            N = NMAX
            I3 = I5+KM
            FIRST3 = .TRUE.
            DO 340 WHILE (FIRST3 .OR. N .GT. MIN)
              I4 = I(N)+KM

              XINT(I3) = XINT(I3)+DXIJ*XINT(I4)
              YINT(I3) = YINT(I3)+DYIJ*YINT(I4)
              ZINT(I3) = ZINT(I3)+DZIJ*ZINT(I4)

              I3 = I4
              N = N-1
              FIRST3 = .FALSE.
  340       END DO
            MIN = MIN+1
            FIRST2 = .FALSE.
  360     END DO
          IF (NIMAX .GT. 0) THEN
            I3 = 49+KM+I1
            DO 400 NJ = 1,NJMAX
              I4 = I3
              DO 380 NI = 1,NIMAX

                XINT(I4) = XINT(I4+294)+DXIJ*XINT(I4-49)
                YINT(I4) = YINT(I4+294)+DYIJ*YINT(I4-49)
                ZINT(I4) = ZINT(I4+294)+DZIJ*ZINT(I4-49)
                I4 = I4+343
  380         CONTINUE
              I3 = I3+49
  400       CONTINUE
          END IF
          M = M+1
          FIRST1 = .FALSE.
  430   END DO
      END IF

! !     ----- I(NI,NJ,NK,NL) -----

      IF (NLMAX .GT. 0) THEN

        I5 = K(MMAX+1)
        IA = I1
        NI = 0
        FIRST4 = .TRUE.
        DO 580 WHILE (FIRST4 .OR. NI .LE. NIMAX)
          NJ = 0
          IB = IA
          FIRST1 = .TRUE.
          DO 570 WHILE (FIRST1 .OR. NJ .LE. NJMAX)
            MIN = NKMAX
            FIRST2 = .TRUE.
            DO 530 WHILE (FIRST2 .OR. MIN .LT. MMAX)
              M = MMAX
              I3 = IB+I5
              FIRST3 = .TRUE.
              DO 520 WHILE (FIRST3 .OR. M .GT. MIN)
                I4 = IB+K(M)

                XINT(I3) = XINT(I3)+DXKL*XINT(I4)
                YINT(I3) = YINT(I3)+DYKL*YINT(I4)
                ZINT(I3) = ZINT(I3)+DZKL*ZINT(I4)

                I3 = I4
                M = M-1
                FIRST3 = .FALSE.
  520         END DO
              MIN = MIN+1
              FIRST2 = .FALSE.
  530       END DO
            IF (NKMAX .GT. 0) THEN
              I3 = IB+1
              DO 560 NL = 1,NLMAX
                I4 = I3
                DO 540 NK = 1,NKMAX

                  XINT(I4) = XINT(I4+6)+DXKL*XINT(I4-1)
                  YINT(I4) = YINT(I4+6)+DYKL*YINT(I4-1)
                  ZINT(I4) = ZINT(I4+6)+DZKL*ZINT(I4-1)

                  I4 = I4+7
  540           END DO
              I3 = I3+1
  560         END DO
            END IF
            NJ = NJ+1
            IB = IB+49
            FIRST1 = .FALSE.
  570     END DO
          NI = NI+1
          IA = IA+343
          FIRST4 = .FALSE.
  580   END DO
      END IF
      RETURN
      END
      SUBROUTINE XYZINT_gpu_2222(XINT,YINT,ZINT,&
      I,K,NIMAX,NJMAX,NKMAX,NLMAX,NMAX,MMAX,&
      BP01,B00,B10,XCP00,XC00,YCP00,YC00,ZCP00,ZC00,F00,&
      DXIJ,DYIJ,DZIJ,DXKL,DYKL,DZKL)


      !use integral_tmp_mod, only: XINT
      use omp_lib

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      !LOGICAL N0,N1,M0,M1,FIRST1,FIRST2,FIRST3,FIRST4
      LOGICAL FIRST1,FIRST2,FIRST3,FIRST4
      DIMENSION :: I(13),K(13)
      INTEGER :: NIMAX,NJMAX,NKMAX,NLMAX,NMAX,MMAX
      double precision :: BP01,B00,B10,XCP00,XC00,YCP00,YC00,ZCP00,ZC00
      double precision :: F00,DXIJ,DYIJ,DZIJ,DXKL,DYKL,DZKL
      !double precision :: XINT(31213,ID),YINT(31213,ID),ZINT(31213,ID)
      double precision :: XINT(31213),YINT(31213),ZINT(31213)
      integer :: ID

      PARAMETER (ZERO=0.0D+00)
      PARAMETER (ONE=1.0D+00)
!$omp declare target

!ID=omp_get_thread_num()
      !write(*,*) "ID inside", ID


      ! N0 = NMAX .EQ. 0
      ! N1 = NMAX .LE. 1
      ! M0 = MMAX .EQ. 0
      ! M1 = MMAX .LE. 1
! !     ----- I(0,0) -----
      I1 = I(1)

      ! XINT(I1,ID) = ONE
      ! YINT(I1,ID) = ONE
      ! ZINT(I1,ID) = F00
      XINT(I1) = ONE
      YINT(I1) = ONE
      ZINT(I1) = F00      
      !IF (N0 .AND. M0) RETURN
      I2 = I(2)
      K2 = K(2)
      CP10 = B00

    ! ----- I(1,0) -----
      !IF (.NOT. N0) THEN
        XINT(I2) = XC00
        YINT(I2) = YC00
        ZINT(I2) = ZC00*F00

        ! XINT(I2,ID) = XC00
        ! YINT(I2,ID) = YC00
        ! ZINT(I2,ID) = ZC00*F00
        !IF (M0) GO TO 120
      !END IF

! ! !     ----- I(0,1) -----
      I3 = I1+K2

      ! XINT(I3,ID) = XCP00
      ! YINT(I3,ID) = YCP00
      ! ZINT(I3,ID) = ZCP00*F00
      XINT(I3) = XCP00
      YINT(I3) = YCP00
      ZINT(I3) = ZCP00*F00

! ! ! ! !     ----- I(1,1) -----

      !IF (.NOT. N0) THEN
        I3 = I2+K2

        ! XINT(I3,ID) = XCP00*XINT(I2,ID)+CP10
        ! YINT(I3,ID) = YCP00*YINT(I2,ID)+CP10
        ! ZINT(I3,ID) = ZCP00*ZINT(I2,ID)+CP10*F00
        XINT(I3) = XCP00*XINT(I2)+CP10
        YINT(I3) = YCP00*YINT(I2)+CP10
        ZINT(I3) = ZCP00*ZINT(I2)+CP10*F00
        !write(*,*) "ZINT(I3) inside", ZINT(I3)
        !write(*,*) "XINT(I3) inside", XINT(I3)

        !if (varx(I3,ID) .NE. XINT(I3)) write(*,*) "ugh33"

      !END IF

  !120 CONTINUE

      !IF (.NOT. N1) THEN
        C10 = ZERO
        I3 = I1
        I4 = I2
        DO N = 2,NMAX
          C10 = C10+B10

! ! !     ----- I(N,0) -----

          I5 = I(N+1)

          XINT(I5) = C10*XINT(I3)+XC00*XINT(I4)
          YINT(I5) = C10*YINT(I3)+YC00*YINT(I4)
          ZINT(I5) = C10*ZINT(I3)+ZC00*ZINT(I4)
          !write(*,*) "XINT(I5) inside", XINT(I5)
          !IF ( .NOT. M0) THEN
            CP10 = CP10+B00

! ! !     ----- I(N,1) -----

            I3 = I5+K2

            XINT(I3) = XCP00*XINT(I5)+CP10*XINT(I4)
            YINT(I3) = YCP00*YINT(I5)+CP10*YINT(I4)
            ZINT(I3) = ZCP00*ZINT(I5)+CP10*ZINT(I4)

          !END IF
          I3 = I4
          I4 = I5
      ENDDO
      !END IF

      !IF ( .NOT. M1) THEN
        CP01 = ZERO
        C01 = B00
        I3 = I1
        I4 = I1+K2
        DO M = 2,MMAX
          CP01 = CP01+BP01

! !     ----- I(0,M) -----

          I5 = I1+K(M+1)

          XINT(I5) = CP01*XINT(I3)+XCP00*XINT(I4)
          YINT(I5) = CP01*YINT(I3)+YCP00*YINT(I4)
          ZINT(I5) = CP01*ZINT(I3)+ZCP00*ZINT(I4)
! !     ----- I(1,M) -----

          !IF (.NOT. N0) THEN
            C01 = C01+B00
            I3 = I2+K(M+1)

            XINT(I3) = XC00*XINT(I5)+C01*XINT(I4)
            YINT(I3) = YC00*YINT(I5)+C01*YINT(I4)
            ZINT(I3) = ZC00*ZINT(I5)+C01*ZINT(I4)
          !END IF
          I3 = I4
          I4 = I5
        ENDDO
      !END IF

!     ----- I(N,M) -----

      !IF (.NOT. N1 .AND. .NOT. M1) THEN
        C01 = B00
        K3 = K2
        DO M = 2,MMAX
          K4 = K(M+1)
          C01 = C01+B00
          I3 = I1
          I4 = I2
          C10 = B10
          DO N = 2,NMAX
            I5 = I(N+1)

            XINT(I5+K4) = C10*XINT(I3+K4)+XC00*XINT(I4+K4)+C01*XINT(I4+K3)
            YINT(I5+K4) = C10*YINT(I3+K4)+YC00*YINT(I4+K4)+C01*YINT(I4+K3)
            ZINT(I5+K4) = C10*ZINT(I3+K4)+ZC00*ZINT(I4+K4)+C01*ZINT(I4+K3)

            C10 = C10+B10
            I3 = I4
            I4 = I5
          ENDDO
          K3 = K4
        ENDDO
      !END IF

!     ----- I(NI,NJ,M) -----

      !IF (NJMAX .GT. 0) THEN
        M = 0
        I5 = I(NMAX+1)
        FIRST1 = .TRUE.
        DO 430 WHILE (FIRST1 .OR. M .LE. MMAX)
          MIN = NIMAX
          KM = K(M+1)
          FIRST2 = .TRUE.
          DO 360 WHILE (FIRST2 .OR. MIN .LT. NMAX)
            N = NMAX
            I3 = I5+KM
            FIRST3 = .TRUE.
            DO 340 WHILE (FIRST3 .OR. N .GT. MIN)
              I4 = I(N)+KM

              XINT(I3) = XINT(I3)+DXIJ*XINT(I4)
              YINT(I3) = YINT(I3)+DYIJ*YINT(I4)
              ZINT(I3) = ZINT(I3)+DZIJ*ZINT(I4)

              I3 = I4
              N = N-1
              FIRST3 = .FALSE.
  340       END DO
            MIN = MIN+1
            FIRST2 = .FALSE.
  360     END DO
          !IF (NIMAX .GT. 0) THEN
            I3 = 49+KM+I1
            DO 400 NJ = 1,NJMAX
              I4 = I3
              DO 380 NI = 1,NIMAX

                XINT(I4) = XINT(I4+294)+DXIJ*XINT(I4-49)
                YINT(I4) = YINT(I4+294)+DYIJ*YINT(I4-49)
                ZINT(I4) = ZINT(I4+294)+DZIJ*ZINT(I4-49)
                I4 = I4+343
  380         CONTINUE
              I3 = I3+49
  400       CONTINUE
          !END IF
          M = M+1
          FIRST1 = .FALSE.
  430   END DO
      !END IF

! !     ----- I(NI,NJ,NK,NL) -----

      !IF (NLMAX .GT. 0) THEN

        I5 = K(MMAX+1)
        IA = I1
        NI = 0
        FIRST4 = .TRUE.
        DO 580 WHILE (FIRST4 .OR. NI .LE. NIMAX)
          NJ = 0
          IB = IA
          FIRST1 = .TRUE.
          DO 570 WHILE (FIRST1 .OR. NJ .LE. NJMAX)
            MIN = NKMAX
            FIRST2 = .TRUE.
            DO 530 WHILE (FIRST2 .OR. MIN .LT. MMAX)
              M = MMAX
              I3 = IB+I5
              FIRST3 = .TRUE.
              DO 520 WHILE (FIRST3 .OR. M .GT. MIN)
                I4 = IB+K(M)

                XINT(I3) = XINT(I3)+DXKL*XINT(I4)
                YINT(I3) = YINT(I3)+DYKL*YINT(I4)
                ZINT(I3) = ZINT(I3)+DZKL*ZINT(I4)

                I3 = I4
                M = M-1
                FIRST3 = .FALSE.
  520         END DO
              MIN = MIN+1
              FIRST2 = .FALSE.
  530       END DO
            !IF (NKMAX .GT. 0) THEN
              I3 = IB+1
              DO 560 NL = 1,NLMAX
                I4 = I3
                DO 540 NK = 1,NKMAX

                  XINT(I4) = XINT(I4+6)+DXKL*XINT(I4-1)
                  YINT(I4) = YINT(I4+6)+DYKL*YINT(I4-1)
                  ZINT(I4) = ZINT(I4+6)+DZKL*ZINT(I4-1)
                  !write(*,*) "XINT(I4) inside", XINT(I4)

                  I4 = I4+7
  540           END DO
              I3 = I3+1
  560         END DO
            !END IF
            NJ = NJ+1
            IB = IB+49
            FIRST1 = .FALSE.
  570     END DO
          NI = NI+1
          IA = IA+343
          FIRST4 = .FALSE.
  580   END DO
      !END IF
      RETURN
      END

