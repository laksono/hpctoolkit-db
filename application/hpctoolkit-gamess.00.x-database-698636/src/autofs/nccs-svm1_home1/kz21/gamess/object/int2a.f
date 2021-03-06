C*MODULE INT2A   *DECK BASCHK                                           
      SUBROUTINE BASCHK(LMAX)                                           
      use mx_limits, only: mxsh,mxgtot                                  
C                                                                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C                                                                       
C                                                                       
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),      
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),      
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),   
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL           
C                                                                       
C     RETURN THE HIGHEST ANGULAR MOMENTUM PRESENT IN THE BASIS.         
C     NOTE THAT KTYPE=1,2,3,4,5 MEANS S, P(L), D, F, G FUNCTION.        
C                                                                       
      KANG = 0                                                          
      DO 100 N=1,NSHELL                                                 
          IF(KTYPE(N).GT.KANG) KANG = KTYPE(N)                          
  100 CONTINUE                                                          
      LMAX = KANG-1                                                     
      RETURN                                                            
      END                                                               
C*MODULE INT2A   *DECK EXCHNG                                           
C>    @brief exchange integral for Schwarz inequality screening         
C>    @ author unknonw                                                  
C>    @date Oct 2018 Peng Xu and Tosaporn Sattasathuchana               
C>    - not print Schwarz inequality overhead for QM-EFP2               
C>                                                                      
      SUBROUTINE EXCHNG(XINTS,GHONDO,DDIJ,NSH2,MAXG,INTTYP)             
      use mx_limits, only: mxsh,mxgtot                                  
C                                                                       
      USE params, ONLY: intomp                                          
C$    USE ompmod                                                        
C                                                                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C                                                                       
      DIMENSION XINTS(NSH2),GHONDO(MAXG),DDIJ(*)                        
      DIMENSION IBPOP(4,4)                                              
C                                                                       
      LOGICAL POPLE,OUT,SOME,IANDJ,GOPARR,DSKWRK,MASWRK,SAVEINT,NXT,    
     *        DOESP,LCUT                                                
C                                                                       
C                                                                       
      COMMON /INTAC2/ EI1,EI2,CUX                                       
      COMMON /FLIPS / IB(4,3)                                           
C$omp threadprivate(/FLIPS/)
      COMMON /FMOINF/ NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY              
      common /fmoopt/ espsca(9),RESPAP(2),rESPPC(2),rESDIM,restri(4),   
     *                rcorsd,respct,convfg,cnvdmp,coroff,rflmo(4),      
     *                orshft,orshft2,cnvafo,ascreen(4),IXESP,mxitfg,    
     *                nguess,NBSSE,modorb,modpar,irststp,irstlay,nprfmo,
     *                nfmopal,modprp,maxl1c,ipieda,modgrd,modesp,ivmul, 
     *                modlmo,nopden,mofock,modfd,modfmm,ncentm,ndualb   
      COMMON /FMORUN/ ESPSCF,E0SCF(2),EMP2S,IDAFMO,ICURFG,JCURFG,KCURFG,
     *                ICURLAY,ICURUNT,NAT1E,NCURSH,NGAU,ICURPOP,IFMOSTP,
     *                MONCOR,NEEDR,MODRST,NORBPROJ,NUNESP,ISKIPESP,     
     *                IESDPPC,IDOPROP,MP2RUN,ICURIT,IDMFMO,IDDFMO,      
     *                IDDCUR,NDDLEFT,IVMFMO,nzmtfmo,ifmobas,itmfmo(2)   
      COMMON /GOUT  / GPOPLE(768),NORGP                                 
C$omp threadprivate(/GOUT/)
      COMMON /INTDEX/ IJGT(784),IJX(784),IJY(784),IJZ(784),IK(784),     
     *                KLGT(784),KLX(784),KLY(784),KLZ(784)              
C$omp threadprivate(/INTDEX/)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)                
      COMMON /ELGIDX/ LCUT                                              
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),      
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),      
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),   
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL           
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK                 
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
      COMMON /SHLEXC/ NORGSH(3),NORGSP(3),IEXCH,NANGM,NGTH(4)           
      COMMON /SHLG70/ ISH,JSH,KSH,LSH,IJKLXX(4)                         
C$omp threadprivate(/SHLG70/)
      COMMON /SHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,          
     *                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,          
     *                NIJ,IJ,KL,IJKL                                    
C$omp threadprivate(/SHLNOS/)
      COMMON /SHLT  / TOL,CUTOFF,ICOUNT,OUT                             
C                                                                       
      PARAMETER (ZERO=0.0D+00)                                          
      PARAMETER (ONE=1.0D+00)                                           
C                                                                       
      DATA IBPOP/0,0,0,0,64,16,4,1,128,32,8,2,192,48,12,3/              
C                                                                       
C     ----- COMPUTE ALL EXCHANGE INTEGRALS [II,JJ/II,JJ] -----          
C     THE LARGEST EXCHANGE INTEGRAL FROM EACH INTEGRAL BLOCK            
C     IS PICKED OUT AND SAVED, SO THAT THE SCHWARZ INEQUALITY           
C     CAN BE USED LATER TO AVOID ENTIRE INTEGRAL BLOCKS.                
C                                                                       
C                                                                       
C     CALL THREADED VERSION IF POSSIBLE                                 
C                                                                       
C       The code below is only active                                   
C       with -fopenmp/-qopenmp compiler option                          
C                                                                       
C$    IF (intomp.ne.0) THEN                                             
C$      CALL ompmod_exchng(xints,nsh2,maxg,inttyp)                      
C$      RETURN                                                          
C$    END IF                                                            
C                                                                       
C                                                                       
C     ----- INTIALIZE PARALLEL -----                                    
C     DYNAMIC LOAD BALANCING PROVED AN ENTIRE FAILURE, IN BOTH 1ST      
C     AND 2ND NESTED SHELL LOOP, PRESUMABLY DUE TO TOO HIGH RATIO       
C     OF EXTRA EXPENSES/COMPUTIONS (ON GIGABIT NETWORK).  THE CODE      
C     IS KEPT FOR FUTURE AND LOAD BALANCING IS HARDWIRED TO STATIC.     
C     NXT = IBTYP.EQ.1                                                  
      NXT=.FALSE.                                                       
      IPCOUNT = ME - 1                                                  
      IF(NXT) IPCOUNT = - 1                                             
      NEXT = -1                                                         
      DOESP=NFG.NE.0.AND.NCURSH.NE.0                                    
      IF (GOPARR.OR.DOESP) CALL VCLR(XINTS,1,NSH2)                      
C                                                                       
      SOME = NPRINT.NE.-5 .AND. MASWRK                                  
      SAVEINT=NFG.EQ.0.OR.NCURSH.EQ.0.OR.                               
     *        (IFMOSTP.EQ.6.and.iand(mofock,128).eq.0)                  
C                                                                       
      IF(SOME) THEN                                                     
         TIM = ZERO                                                     
         CALL TSECND(TIM)                                               
         TIM0 = TIM                                                     
      ENDIF                                                             
C                                                                       
      CALL BASCHK(LMAX)                                                 
                    NANGM =  4                                          
      IF(LMAX.EQ.2) NANGM =  6                                          
      IF(LMAX.EQ.3) NANGM = 10                                          
      IF(LMAX.EQ.4) NANGM = 15                                          
      IF(LMAX.EQ.5) NANGM = 21                                          
      IF(LMAX.EQ.6) NANGM = 28                                          
      NGTH(4) = 1                                                       
      NGTH(3) = NGTH(4) * NANGM                                         
      NGTH(2) = NGTH(3) * NANGM                                         
      NGTH(1) = NGTH(2) * NANGM                                         
      IF(NOPK.EQ.0) THEN                                                
         NORGSH(1) = 0                                                  
         NORGSH(2) = NORGSH(1) + NANGM**4                               
         NORGSH(3) = NORGSH(2) + NANGM**4                               
         NORGSP(1) = 0                                                  
         NORGSP(2) = 256                                                
         NORGSP(3) = 512                                                
      ELSE                                                              
         DO I=1,3                                                       
            NORGSH(I) = 0                                               
            NORGSP(I) = 0                                               
         ENDDO                                                          
      END IF                                                            
C                                                                       
C        THE IDEA IS TO DO EVEN SMALL INTEGRALS, BELOW THE USUAL        
C        CUTOFF THRESHHOLDS, BY RESETTING TOLERANCES TIGHTLY.           
C                                                                       
      TOLSV = TOL                                                       
      TOL = 75.0D+00                                                    
C                                                                       
      EI1SV = EI1                                                       
      EI2SV = EI2                                                       
      CUXSV = CUX                                                       
      EI1 = 1.0D-17                                                     
      EI2 = 1.0D-17                                                     
      CUX = 50.0D+00                                                    
C                                                                       
      IEXCH = 1                                                         
      NORGP = 0                                                         
      QQ4   = ONE                                                       
      NINT  = 0                                                         
C                                                                       
C     ----- LOOP OVER ALL SHELL BLOCKS -----                            
C                                                                       
      IJIJ = 0                                                          
      DO 600 ISH = 1,NSHELL                                             
         IF (GOPARR.AND.NXT) THEN                                       
            IPCOUNT = IPCOUNT + 1                                       
            IF (IPCOUNT.GT.NEXT) CALL DDI_DLBNEXT(NEXT)                 
            IF (NEXT.NE.IPCOUNT) THEN                                   
               IJIJ = IJIJ+ISH                                          
               GO TO 600                                                
            ENDIF                                                       
         END IF                                                         
         DO 500 JSH = 1,ISH                                             
            IJIJ = IJIJ+1                                               
C                                                                       
C           SKIP UNNEEDED OFF-DIAGONAL BLOCKS FOR FMO ESP SCREENING.    
C                                                                       
            IF(DOESP.AND.ISH.GT.NCURSH.AND.JSH.LE.NCURSH) GO TO 500     
C                                                                       
C           ----- GO PARALLEL! -----                                    
C                                                                       
            IF (GOPARR.AND..NOT.NXT) THEN                               
               IPCOUNT = IPCOUNT + 1                                    
               IF (MOD(IPCOUNT,NPROC).NE.0) GO TO 500                   
            END IF                                                      
C                                                                       
C     USE POPLE CODE FOR ANY PURE SP INTEGRAL BLOCKS,                   
C     USE HONDO RYS POLYNOMIAL CODE FOR OTHER BLOCKS                    
C                                                                       
            POPLE=.TRUE.                                                
            IF(INTTYP.GE.2)     POPLE=.FALSE.                           
            IF(KTYPE(ISH).GT.2) POPLE=.FALSE.                           
            IF(KTYPE(JSH).GT.2) POPLE=.FALSE.                           
C                                                                       
            IF(POPLE) THEN                                              
               KSH=ISH                                                  
               LSH=JSH                                                  
               CALL GENR70(IEXCH,.FALSE.)                               
            ELSE                                                        
               CALL SHELLS(1,ISH,JSH,ISH,JSH,.TRUE.)                    
            CALL IJPRIM(DDIJ)                                           
            CALL SHELLS(2,ISH,JSH,ISH,JSH,.TRUE.)                       
            CALL ZQOUT(GHONDO)                                          
            IF(IJKL.EQ.1) CALL S0000(GHONDO,DDIJ)                       
            IF(IJKL.GT.1) CALL GENRAL(GHONDO,DDIJ)                      
         END IF                                                         
C                                                                       
C     ----- PICK OUT LARGEST EXCHANGE INTEGRAL FOR THIS BLOCK -----     
C                                                                       
         VMAX = ZERO                                                    
         MINI = KMIN(ISH)                                               
         MINJ = KMIN(JSH)                                               
         MAXI = KMAX(ISH)                                               
         JMAX = KMAX(JSH)                                               
         IANDJ=ISH.EQ.JSH                                               
         IBB = IB(1,IEXCH)                                              
         JBB = IB(2,IEXCH)                                              
         KBB = IB(3,IEXCH)                                              
         LBB = IB(4,IEXCH)                                              
         IJN = 0                                                        
         DO 300 I=MINI,MAXI                                             
            IF(IANDJ) JMAX = I                                          
            DO 200 J=MINJ,JMAX                                          
               IF(POPLE) THEN                                           
                  NN = IBPOP(IBB,I) + IBPOP(JBB,J)                      
     *               + IBPOP(KBB,I) + IBPOP(LBB,J) + 1                  
                  VAL = GPOPLE(NN)                                      
               ELSE                                                     
                  IJN = IJN+1                                           
                  NN = IJGT(IJN) + KLGT(IJN)                            
                  VAL = GHONDO(NN)                                      
               END IF                                                   
               IF(VAL.GT.ZERO) NINT=NINT+1                              
               IF(VAL.GT.VMAX) VMAX=VAL                                 
  200       CONTINUE                                                    
  300    CONTINUE                                                       
         XINTS(IJIJ)=SQRT(VMAX)                                         
  500    CONTINUE                                                       
  600 CONTINUE                                                          
C                                                                       
C     ----- SUM UP PARTIAL CONTRIBUTIONS IF PARALLEL -----              
C                                                                       
      IF (GOPARR) THEN                                                  
         CALL DDI_GSUMF(1050,XINTS,NSH2)                                
         CALL DDI_GSUMI(1051,NINT ,1)                                   
         IF(NXT) CALL DDI_DLBRESET                                      
      END IF                                                            
C                                                                       
      IF(OUT) THEN                                                      
         WRITE(IW,*) 'MAX EXCHANGE INTEGRAL IN SHELL'                   
         CALL PRTRI(XINTS,NSHELL)                                       
      END IF                                                            
C                                                                       
      LVLEFP=LEVELEFP()                                                 
C THIS WILL NOT BE PRINT OUT FOR QM-EFP2                                
C      IF(SOME) THEN                                                    
      IF(SOME.AND.(LVLEFP.NE.2)) THEN                                   
         CALL TSECND(TIM)                                               
         TEXCH = TIM-TIM0                                               
         WRITE(IW,9000) NINT,TEXCH                                      
      ENDIF                                                             
C                                                                       
      TOL = TOLSV                                                       
      EI1 = EI1SV                                                       
      EI2 = EI2SV                                                       
      CUX = CUXSV                                                       
C                                                                       
C     DURING FMO ESP RUNS, EXCHANGE INTEGRALS HAVE DIFFERENT SIZE SO    
C     ONE CANNOT WRITE THEM TO THE SAME RECORD. THE ONLY EXCEPTION IS   
C     THE SEPARATED DIMER ENERGIES WHERE THERE IS JUST ONE SET OF 2E    
C     INTEGRALS.                                                        
C     ELONGATION METHOD ALSO MUST DECIDE ON THIS.                       
C                                                                       
      IF(SAVEINT.AND.(.NOT.LCUT)) CALL DAWRIT(IDAF,IODA,XINTS,NSH2,54,0)
      RETURN                                                            
 9000 FORMAT(1X,'SCHWARZ INEQUALITY OVERHEAD:',I10,' INTEGRALS, T=',    
     *       F12.2)                                                     
      END                                                               
C*MODULE INT2A   *DECK DEBUT                                            
      SUBROUTINE DEBUT(DIRSCF,BUFP,BUFK,IX,NINTMX,NEED,DIRTRF)          
      USE camdft, ONLY: CAMFLAG                                         
      USE lrcdft, ONLY: LCFLAG, EMU, EMU2, LRFILE                       
      use mx_limits, only: mxsh,mxgtot,mxatm,mxao                       
C                                                                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C                                                                       
      LOGICAL DIRSCF,DIRTRF                                             
      LOGICAL DIR,GOPARR,DSKWRK,MASWRK,PK,PANDK,BLOCK,OUT               
      LOGICAL LTRMST                                                    
      LOGICAL LRINT                                                     
C                                                                       
      DIMENSION BUFP(NINTMX),BUFK(NINTMX),IX(*)                         
C                                                                       
C                                                                       
      COMMON /B     / CO(MXSH,3)                                        
      COMMON /ELGPMT/ NELONG,NATM,NASPIN,NCT,NBNDAB,NTMLB,IPRI,LDOS     
      COMMON /ELGRST/ IRSTRT,I2EA,IGOOD                                 
      COMMON /ELGTRM/ LTRMST,NFLTRM,NRCTRM,NPSTRM,NHTSHL                
      COMMON /IJPAIR/ IA(MXAO)                                          
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,                    
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)                  
      COMMON /INTPR / QINT(2),VALINT(2),JCINT(11)                       
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)                
      COMMON /NLRCF / LRINT                                             
C$omp threadprivate(/NLRCF /)
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),      
     2                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),      
     3                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),   
     4                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL           
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK                 
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
      COMMON /PCKLAB/ LABSIZ                                            
      COMMON /PKFIL / PK,PANDK,BLOCK                                    
      COMMON /RESTAR/ TIMLIM,IREST,NREC,INTLOC,IST,JST,KST,LST          
      COMMON /INT2IC/ NINTIC,ININTIC,NXXIC,LBUFPIC,LIXIC,LABSIX,NINTIX  
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL                
      COMMON /SHLT  / TOL,CUTOFF,ICOUNT,OUT                             
      COMMON /WFNOPT/ SCFTYP,VBTYP,DFTYPE,TDDFTYP,CITYP,CCTYP,          
     *                MPLEVL,MPCTYP                                     
C                                                                       
      PARAMETER (ONE=1.0D+00, TEN=10.0D+00, RLN10=2.30258D+00)          
C                                                                       
      DATA CHECK/8HCHECK   /                                            
      DATA NONE/4HNONE/                                                 
      ISSAVE=IS                                                         
      IF(LRINT)IS=LRFILE                                                
C                                                                       
C     ----- INITIALIZE TWO ELECTRON INTEGRAL CALCULATION -----          
C                                                                       
      TOL= ITOL*RLN10                                                   
      CUTOFF= ONE/(TEN**ICUT)                                           
      DIR= DIRSCF .OR. DIRTRF                                           
C                                                                       
      DO 100 I=1,NUM                                                    
         IA(I)=(I*I-I)/2                                                
  100 CONTINUE                                                          
C                                                                       
C     ----- ERIC (AND OTHER INTEGRAL PACKAGE) INITIALIZATIONS -----     
C                                                                       
C      CALL ERIPRE                                                      
      IF(LCFLAG.OR.CAMFLAG)THEN                                         
         IF(LRINT)CALL ERIPRE                                           
      ELSE                                                              
         CALL ERIPRE                                                    
      END IF                                                            
C                                                                       
      IF((NPRINT.NE.-5.AND.NPRINT.NE.-23) .AND. MASWRK) THEN            
         IF(DIR) THEN                                                   
C           IF(DIRSCF) WRITE(IW,9000)                                   
            IF(.NOT.LRINT.AND.DIRSCF)WRITE(IW,9000)                     
            IF(     LRINT.AND.DIRSCF)WRITE(IW,9001)                     
            IF(DIRSCF) WRITE(IW,9040)                                   
            IF(DIRTRF) WRITE(IW,9045)                                   
         ELSE                                                           
C           WRITE(IW,9000)                                              
            IF(.NOT.LRINT)WRITE(IW,9000)                                
            IF(     LRINT)WRITE(IW,9001)                                
            IF(PK) THEN                                                 
               IF(     PANDK .AND. NPRINT.NE.-5) WRITE(IW,9010)         
               IF(.NOT.PANDK) WRITE(IW,9020)                            
            ELSE                                                        
               WRITE(IW,9030)                                           
            END IF                                                      
C  J (OR P) INTEGRAL, 4 LABELS OF 1 OR 2 BYTES, PLUS MAYBE A K INTEGRAL 
            NBYTES = 8 + 4*LABSIZ                                       
            IF(PANDK) NBYTES = NBYTES+8                                 
            IF(NINTIC.NE.0) WRITE(IW,9055) NINTIC                       
            WRITE(IW,9050) NINTMX,NBYTES                                
            WRITE(IW,9060) NEED                                         
         END IF                                                         
      END IF                                                            
C                                                                       
      OUT = NPRINT.EQ.4 .AND. MASWRK                                    
      JCINT(1) = 0                                                      
C                                                                       
      DO 200 I=1,NSHELL                                                 
         ICC = KATOM(I)                                                 
         CO(I,1)= C(1,ICC)                                              
         CO(I,2)= C(2,ICC)                                              
         CO(I,3)= C(3,ICC)                                              
  200 CONTINUE                                                          
C                                                                       
      IF(DIRSCF .OR. DIRTRF .OR. EXETYP.EQ.CHECK) GO TO 400             
C                                                                       
      CALL SEQREW(IS)                                                   
C                                                                       
      IF(IREST.LT.1 .OR. NREC.LE.1 .OR. INTLOC.LE.1) GO TO 400          
  300 CONTINUE                                                          
C                                                                       
C     ----- POSITION THE INTEGRAL FILE FOR A RESTART JOB -----          
C                                                                       
      ICOUNT = INTLOC                                                   
      N = NREC-1                                                        
      IF(MASWRK) THEN                                                   
         DO 310 I=1,N                                                   
            READ(IS)                                                    
  310    CONTINUE                                                       
      END IF                                                            
      IF(NINTIC.NE.0) CALL ABRT                                         
C     NO RESTART SUPPORT YET.                                           
      IF(.NOT.PANDK) CALL PREAD(IS,BUFP,IX,NXX,NINTMX)                  
      IF(     PANDK) CALL PKREAD(IS,BUFP,BUFK,IX,NXX,NINTMX)            
      CALL SEQREW(IS)                                                   
      IF(MASWRK) THEN                                                   
         DO 320 I=1,N                                                   
            READ(IS)                                                    
  320    CONTINUE                                                       
      END IF                                                            
      GO TO 999                                                         
  400 CONTINUE                                                          
C                                                                       
C     ----- NORMAL START -----                                          
C                                                                       
      IF(IST.LT.1) IST = 1                                              
      IF(JST.LT.1) JST = 1                                              
      IF(KST.LT.1) KST = 1                                              
      IF(LST.LT.1) LST = 1                                              
      IF(IST.GT.NSHELL) GO TO 999                                       
      IF(IST.NE.1 .OR. JST.NE.1 .OR. KST.NE.1 .OR. LST.NE.1) GO TO 300  
      NREC   = 1                                                        
      INTLOC = 1                                                        
      ICOUNT = 1                                                        
  999 CONTINUE                                                          
C                                                                       
C         ELONGATION METHOD INTEGRAL FILES                              
      IF(NELONG.GE.2) THEN                                              
         IF(MPCTYP.NE.NONE) I2EA = 0                                    
         IF(DIRSCF) I2EA = 0                                            
         CALL ELGINT(IS,BUFP,IX,NINTMX)                                 
         IF(LTRMST) THEN                                                
            NHTSHL = NSHELL                                             
            DO I=1,NSHELL                                               
               IF(KATOM(I).EQ.(NAT-NTMLB+1)) GOTO 1500                  
            ENDDO                                                       
 1500       NHTSHL = I                                                  
         ELSE                                                           
            NFLTRM = 1                                                  
            NRCTRM = 0                                                  
            NPSTRM = 0                                                  
            NHTSHL = NSHELL                                             
         ENDIF                                                          
      ENDIF                                                             
C                                                                       
      IF(LRINT)IS=ISSAVE                                                
      RETURN                                                            
C                                                                       
 9000 FORMAT(/10X,20(1H-)/10X,'2 ELECTRON INTEGRALS'/10X,20(1H-)/)      
 9001 FORMAT(/10X,32(1H-)/10X,'2 ELECTRON LC EXCHANGE INTEGRALS'        
     * /10X,32(1H-)/)                                                   
 9010 FORMAT(' THE -PK- OPTION IS ON, CREATING -P- AND -K-',            
     *       ' SUPERMATRICES.')                                         
 9020 FORMAT(' THE -PK- OPTION IS ON, CREATING A -P- SUPERMATRIX.')     
 9030 FORMAT(' THE -PK- OPTION IS OFF, THE INTEGRALS ARE NOT IN',       
     *       ' SUPERMATRIX FORM.')                                      
 9040 FORMAT(' DIRECT SCF METHOD SKIPS INTEGRAL STORAGE ON DISK.')      
 9045 FORMAT(' DIRECT TRANSFORMATION SKIPS AO INTEGRAL STORAGE',        
     *       ' ON DISK.')                                               
 9050 FORMAT(' STORING',I8,' INTEGRALS/RECORD ON DISK, USING',I3,       
     *       ' BYTES/INTEGRAL.')                                        
 9055 FORMAT(' STORING',I9,' INTEGRALS IN MEMORY,')                     
 9060 FORMAT(' TWO ELECTRON INTEGRAL EVALUATION REQUIRES',I8,           
     *       ' WORDS OF MEMORY.')                                       
      END                                                               
C*MODULE INT2A   *DECK FINAL                                            
      SUBROUTINE FINAL(INDEX,II,JJ,KK,LL,PANDK,BUFP,BUFK,IX,NINTMX)     
      USE lrcdft, ONLY: LRFILE                                          
      use mx_limits, only: mxsh,mxgtot                                  
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      LOGICAL PANDK,OUT,GOPARR,DSKWRK,MASWRK,LTRMST                     
      LOGICAL LRINT                                                     
C                                                                       
      DIMENSION BUFP(NINTMX),BUFK(NINTMX),IX(*)                         
C                                                                       
C                                                                       
      COMMON /ELGPMT/ NELONG,NATM,NASPIN,NCT,NBNDAB,NTMLB,IPRI,LDOS     
      COMMON /ELGTRM/ LTRMST,NFLTRM,NRCTRM,NPSTRM,NHTSHL                
      COMMON /FMOINF/ NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY              
      COMMON /INT2IC/ NINTIC,ININTIC,NXXIC,LBUFPIC,LIXIC,LABSIX,NINTIX  
      COMMON /INTPR / Q(2),V(2),JC,N1(2),J1(2),J2(2),J3(2),J4(2)        
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)                
      COMMON /NLRCF / LRINT                                             
C$omp threadprivate(/NLRCF /)
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),      
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),      
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),   
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL           
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK                 
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
      COMMON /RESTAR/ TIMLIM,IREST,NREC,INTLOC,IST,JST,KST,LST          
      COMMON /SHLT  / TOL,CUTOFF,ICOUNT,OUT                             
C                                                                       
      ISSAVE=IS                                                         
      IF(LRINT)IS=LRFILE                                                
C                                                                       
      IF(OUT .AND. JC.GT.0)                                             
     * WRITE(IW,9088) (J1(M),J2(M),J3(M),J4(M),Q(M),V(M),M=1,JC)        
C                                                                       
      IF (INDEX.EQ.1) GO TO 140                                         
C                                                                       
C         RAN OUT OF TIME, FLUSH PARTIAL BUFFER, PRINT RESTART DATA     
C                                                                       
      IREST = 1                                                         
      IST = II                                                          
      JST = JJ                                                          
      KST = KK                                                          
      LST = LL+1                                                        
      IF(LST.LE.KK) GO TO 120                                           
      LST = 1                                                           
      KST = KK+1                                                        
      IF(KST.LE.JJ) GO TO 120                                           
      KST = 1                                                           
      JST = JJ+1                                                        
      IF(JST.LE.II) GO TO 120                                           
      JST = 1                                                           
      IST = II+1                                                        
      IF(IST.GT.NSHELL) GO TO 140                                       
C                                                                       
  120 CONTINUE                                                          
CNB   DO SOMETHING FOR IN CORE                                          
      NXX = ICOUNT-1                                                    
      IF(.NOT.PANDK) CALL PWRIT(IS,BUFP,IX,NXX,NINTMX)                  
      IF(     PANDK) CALL PKWRIT(IS,BUFP,BUFK,IX,NXX,NINTMX)            
      NINT = NINTMX*(NREC-1)+ICOUNT-1                                   
      IF (MASWRK) THEN                                                  
         WRITE(IW,9010) NINT,NREC,IS                                    
         WRITE(IW,9020) NREC,ICOUNT,IST,JST,KST,LST                     
      END IF                                                            
      IF(LRINT)IS=ISSAVE                                                
      RETURN                                                            
C                                                                       
C        DONE WITH INTEGRALS, WRITE LAST BUFFER, PRINT STATISTICS       
C                                                                       
  140 CONTINUE                                                          
      IREST = 0                                                         
      IST = 1                                                           
      JST = 1                                                           
      KST = 1                                                           
      LST = 1                                                           
      NXX = ICOUNT-1                                                    
      IF(NXX.GE.NINTIC) THEN                                            
C                                                                       
C     GE INSTEAD OF GT ENFORCES WRITING A RECORD WITH 0 INTEGRALS       
C     (IF NINTIC IS EQUAL TO NXX) TO COMPLY WITH THE INTEGRAL FILE      
C     STRUCTURE                                                         
C                                                                       
         NXX=NXX-NINTIC                                                 
         NXX = -NXX                                                     
C        WRITE(6,*) 'SAVING',NXX                                        
         IF(.NOT.PANDK) CALL PWRIT(IS,BUFP(NINTIC+1),IX(ININTIC+1),     
     *                             NXX,NINTMX)                          
         IF(     PANDK) CALL PKWRIT(IS,BUFP,BUFK,IX,NXX,NINTMX)         
         NINT = NINTMX*(NREC-1)+ICOUNT-1                                
C        NOTE THAT ICOUNT-1 INCLUDES NINTIC                             
      ELSE                                                              
         NINT = NXX                                                     
      ENDIF                                                             
C     THE # OF INTEGRALS MUST BE STORED ON EACH NODE BEFORE GSUMI.      
      IF(NINTIC.NE.0) NXXIC=MIN(NINT,NINTIC)                            
      NINTMY=NINT                                                       
C                                                                       
C          ELONGATION METHOD INTEGRAL FILES                             
C                                                                       
      IF(NELONG.GT.1.AND.(.NOT.LTRMST)) THEN                            
         NRCTRM = NREC                                                  
         NPSTRM = ICOUNT                                                
      END IF                                                            
C                                                                       
      IF (GOPARR) THEN                                                  
         CALL DDI_GSUMI(1056,NINT,1)                                    
         CALL DDI_GSUMI(1057,NREC,1)                                    
      END IF                                                            
      IF (MASWRK.AND.(NFG.EQ.0.OR.NPRINT.NE.-5)) THEN                   
         IF(NINTIC.EQ.0) THEN                                           
            WRITE(IW,9010) NINT,NREC,IS                                 
         ELSE IF(NINTMY.LE.NINTIC) THEN                                 
            WRITE(IW,9015) NINT,(NINTMY*1.0D+02)/NINTIC                 
         ELSE                                                           
            WRITE(IW,9017) NINT,NINTIC,(NINTIC*1.0D+02)/NINTMY,         
     *                     NREC,IS                                      
         ENDIF                                                          
      ENDIF                                                             
      IF(LRINT)IS=ISSAVE                                                
      RETURN                                                            
C                                                                       
 9010 FORMAT(1X,'TOTAL NUMBER OF NONZERO TWO-ELECTRON INTEGRALS =',I20/ 
     *       1X,I10,' INTEGRAL RECORDS WERE STORED ON DISK FILE',I3,'.')
 9015 FORMAT(1X,'TOTAL NUMBER OF NONZERO TWO-ELECTRON INTEGRALS =',I20/ 
     *       1X,'ALL INTEGRALS FITTED INTO MEMORY (',F5.1,              
     *          '% OF BUFFER USED)!')                                   
 9017 FORMAT(1X,'TOTAL NUMBER OF NONZERO TWO-ELECTRON INTEGRALS =',I20  
     *      /1X,'ON NODE 0',I13,' INTEGRALS STORED IN MEMORY(',F4.1,'%)'
     *     ,/I10,' INTEGRAL RECORDS WERE STORED ON DISK FILE',I3,'.')   
 9020 FORMAT(/1X,'...... WARNING   .......   WARNING   .......'/        
     *        1X,'TIME LIMIT HAS EXPIRED.  THIS JOB MUST BE RESTARTED.'/
     *        1X,'IF YOU SAVED THE INTEGRALS, RESTART WITH IREST=1,'/   
     *        1X,'NREC=',I8,' INTLOC=',I6,' IST,JST,KST,LST=',4I6)      
 9088 FORMAT(2(4I4,F5.1,F17.9,1X))                                      
      END                                                               
C*MODULE INT2A   *DECK FORMS_gpu2                                       
      SUBROUTINE FORMS_gpu2(GHONDO,NROOTS,DKL,DIJ,XIN,YIN,ZIN,          
     * IJGT,IJX,IJY,IJZ,IK,KLGT,KLX,KLY,KLZ,IJ)                         
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      DIMENSION GHONDO(*)                                               
C                                                                       
C      COMMON /DENS  / DKL(784),DIJ(784)                                
                                                                        
      DIMENSION :: DKL(784),DIJ(784)                                    
      DIMENSION :: IJGT(784),IJX(784),IJY(784),IJZ(784),IK(784)         
      DIMENSION :: KLGT(784),KLX(784),KLY(784),KLZ(784)                 
      DIMENSION :: XIN(31213),YIN(31213),ZIN(31213)                     
                                                                        
C      COMMON /INTDEX/ IJGT(784),IJX(784),IJY(784),IJZ(784),IK(784),    
C     *                KLGT(784),KLX(784),KLY(784),KLZ(784)             
C      COMMON /SHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,         
C     *                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,         
C     *                NIJ,IJ,KL,IJKL                                   
C      COMMON /ROOT  / XX,U(13),W(13),NROOTS                            
C      COMMON /XYZ   / XIN(31213),YIN(31213),ZIN(31213)                 
C                                                                       
C     ----- FORM INTEGRALS OVER FUNCTIONS -----                         
C     DIMENSIONING XIN(81,5), AND ROLLING UP THE COMPUTATION            
C     OF GHONDO IN A LOOP OF LENGTH NROOTS ADDS 33 SECONDS TO           
C     A 240 SECOND INTEGRAL COMPUTATION JOB.  LEAVE IT UNROLLED.        
C                                                                       
      GO TO (10,20,30,40,50,60,70,80,90,100,110,120,130),NROOTS         
C                                                                       
C          CODE FOR NROOTS=1                                            
C                                                                       
   10 CONTINUE                                                          
      DO 12 I = 1,IJ                                                    
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO 11 K = 1,MAX                                                   
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*                                
     *            ( XIN(MX      )*YIN(MY      )*ZIN(MZ      ))          
   11 CONTINUE                                                          
   12 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C          CODE FOR NROOTS=2                                            
C                                                                       
   20 CONTINUE                                                          
      DO 22 I = 1,IJ                                                    
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO 21 K = 1,MAX                                                   
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*                                
     *            ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )           
     *          +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401))          
   21 CONTINUE                                                          
   22 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C          CODE FOR NROOTS=3                                            
C                                                                       
   30 CONTINUE                                                          
      DO 32 I = 1,IJ                                                    
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO 31 K = 1,MAX                                                   
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*                                
     *            ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )           
     *          +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401)           
     *          +   XIN(MX+ 4802)*YIN(MY+ 4802)*ZIN(MZ+ 4802))          
   31 CONTINUE                                                          
   32 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C          CODE FOR NROOTS=4                                            
C                                                                       
   40 CONTINUE                                                          
      DO 42 I = 1,IJ                                                    
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO 41 K = 1,MAX                                                   
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*                                
     *            ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )           
     *          +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401)           
     *          +   XIN(MX+ 4802)*YIN(MY+ 4802)*ZIN(MZ+ 4802)           
     *          +   XIN(MX+ 7203)*YIN(MY+ 7203)*ZIN(MZ+ 7203))          
   41 CONTINUE                                                          
   42 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C          CODE FOR NROOTS=5                                            
C                                                                       
   50 CONTINUE                                                          
      DO 52 I = 1,IJ                                                    
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO 51 K = 1,MAX                                                   
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*                                
     *            ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )           
     *          +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401)           
     *          +   XIN(MX+ 4802)*YIN(MY+ 4802)*ZIN(MZ+ 4802)           
     *          +   XIN(MX+ 7203)*YIN(MY+ 7203)*ZIN(MZ+ 7203)           
     *          +   XIN(MX+ 9604)*YIN(MY+ 9604)*ZIN(MZ+ 9604))          
   51 CONTINUE                                                          
   52 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C          CODE FOR NROOTS=6                                            
C                                                                       
   60 CONTINUE                                                          
      DO 62 I = 1,IJ                                                    
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO 61 K = 1,MAX                                                   
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*                                
     *            ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )           
     *          +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401)           
     *          +   XIN(MX+ 4802)*YIN(MY+ 4802)*ZIN(MZ+ 4802)           
     *          +   XIN(MX+ 7203)*YIN(MY+ 7203)*ZIN(MZ+ 7203)           
     *          +   XIN(MX+ 9604)*YIN(MY+ 9604)*ZIN(MZ+ 9604)           
     *          +   XIN(MX+12005)*YIN(MY+12005)*ZIN(MZ+12005))          
   61 CONTINUE                                                          
   62 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C          CODE FOR NROOTS=7                                            
C                                                                       
   70 CONTINUE                                                          
      DO 72 I = 1,IJ                                                    
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO 71 K = 1,MAX                                                   
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*                                
     *            ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )           
     *          +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401)           
     *          +   XIN(MX+ 4802)*YIN(MY+ 4802)*ZIN(MZ+ 4802)           
     *          +   XIN(MX+ 7203)*YIN(MY+ 7203)*ZIN(MZ+ 7203)           
     *          +   XIN(MX+ 9604)*YIN(MY+ 9604)*ZIN(MZ+ 9604)           
     *          +   XIN(MX+12005)*YIN(MY+12005)*ZIN(MZ+12005)           
     *          +   XIN(MX+14406)*YIN(MY+14406)*ZIN(MZ+14406))          
   71 CONTINUE                                                          
   72 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C          CODE FOR NROOTS=8                                            
C                                                                       
   80 CONTINUE                                                          
      DO 82 I = 1,IJ                                                    
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO 81 K = 1,MAX                                                   
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*                                
     *            ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )           
     *          +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401)           
     *          +   XIN(MX+ 4802)*YIN(MY+ 4802)*ZIN(MZ+ 4802)           
     *          +   XIN(MX+ 7203)*YIN(MY+ 7203)*ZIN(MZ+ 7203)           
     *          +   XIN(MX+ 9604)*YIN(MY+ 9604)*ZIN(MZ+ 9604)           
     *          +   XIN(MX+12005)*YIN(MY+12005)*ZIN(MZ+12005)           
     *          +   XIN(MX+14406)*YIN(MY+14406)*ZIN(MZ+14406)           
     *          +   XIN(MX+16807)*YIN(MY+16807)*ZIN(MZ+16807))          
   81 CONTINUE                                                          
   82 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C          CODE FOR NROOTS=9                                            
C                                                                       
   90 CONTINUE                                                          
      DO 92 I = 1,IJ                                                    
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO 91 K = 1,MAX                                                   
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*                                
     *            ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )           
     *          +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401)           
     *          +   XIN(MX+ 4802)*YIN(MY+ 4802)*ZIN(MZ+ 4802)           
     *          +   XIN(MX+ 7203)*YIN(MY+ 7203)*ZIN(MZ+ 7203)           
     *          +   XIN(MX+ 9604)*YIN(MY+ 9604)*ZIN(MZ+ 9604)           
     *          +   XIN(MX+12005)*YIN(MY+12005)*ZIN(MZ+12005)           
     *          +   XIN(MX+14406)*YIN(MY+14406)*ZIN(MZ+14406)           
     *          +   XIN(MX+16807)*YIN(MY+16807)*ZIN(MZ+16807)           
     *          +   XIN(MX+19208)*YIN(MY+19208)*ZIN(MZ+19208))          
   91 CONTINUE                                                          
   92 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C          CODE FOR NROOTS=10                                           
C                                                                       
  100 CONTINUE                                                          
      DO 102 I = 1,IJ                                                   
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO 101 K = 1,MAX                                                  
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*                                
     *            ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )           
     *          +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401)           
     *          +   XIN(MX+ 4802)*YIN(MY+ 4802)*ZIN(MZ+ 4802)           
     *          +   XIN(MX+ 7203)*YIN(MY+ 7203)*ZIN(MZ+ 7203)           
     *          +   XIN(MX+ 9604)*YIN(MY+ 9604)*ZIN(MZ+ 9604)           
     *          +   XIN(MX+12005)*YIN(MY+12005)*ZIN(MZ+12005)           
     *          +   XIN(MX+14406)*YIN(MY+14406)*ZIN(MZ+14406)           
     *          +   XIN(MX+16807)*YIN(MY+16807)*ZIN(MZ+16807)           
     *          +   XIN(MX+19208)*YIN(MY+19208)*ZIN(MZ+19208)           
     *          +   XIN(MX+21609)*YIN(MY+21609)*ZIN(MZ+21609))          
  101 CONTINUE                                                          
  102 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C          CODE FOR NROOTS=11                                           
C                                                                       
  110 CONTINUE                                                          
      DO 112 I = 1,IJ                                                   
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO 111 K = 1,MAX                                                  
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*                                
     *            ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )           
     *          +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401)           
     *          +   XIN(MX+ 4802)*YIN(MY+ 4802)*ZIN(MZ+ 4802)           
     *          +   XIN(MX+ 7203)*YIN(MY+ 7203)*ZIN(MZ+ 7203)           
     *          +   XIN(MX+ 9604)*YIN(MY+ 9604)*ZIN(MZ+ 9604)           
     *          +   XIN(MX+12005)*YIN(MY+12005)*ZIN(MZ+12005)           
     *          +   XIN(MX+14406)*YIN(MY+14406)*ZIN(MZ+14406)           
     *          +   XIN(MX+16807)*YIN(MY+16807)*ZIN(MZ+16807)           
     *          +   XIN(MX+19208)*YIN(MY+19208)*ZIN(MZ+19208)           
     *          +   XIN(MX+21609)*YIN(MY+21609)*ZIN(MZ+21609)           
     *          +   XIN(MX+24010)*YIN(MY+24010)*ZIN(MZ+24010))          
  111 CONTINUE                                                          
  112 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C          CODE FOR NROOTS=12                                           
C                                                                       
  120 CONTINUE                                                          
      DO 122 I = 1,IJ                                                   
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO 121 K = 1,MAX                                                  
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*                                
     *            ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )           
     *          +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401)           
     *          +   XIN(MX+ 4802)*YIN(MY+ 4802)*ZIN(MZ+ 4802)           
     *          +   XIN(MX+ 7203)*YIN(MY+ 7203)*ZIN(MZ+ 7203)           
     *          +   XIN(MX+ 9604)*YIN(MY+ 9604)*ZIN(MZ+ 9604)           
     *          +   XIN(MX+12005)*YIN(MY+12005)*ZIN(MZ+12005)           
     *          +   XIN(MX+14406)*YIN(MY+14406)*ZIN(MZ+14406)           
     *          +   XIN(MX+16807)*YIN(MY+16807)*ZIN(MZ+16807)           
     *          +   XIN(MX+19208)*YIN(MY+19208)*ZIN(MZ+19208)           
     *          +   XIN(MX+21609)*YIN(MY+21609)*ZIN(MZ+21609)           
     *          +   XIN(MX+24010)*YIN(MY+24010)*ZIN(MZ+24010)           
     *          +   XIN(MX+26411)*YIN(MY+26411)*ZIN(MZ+26411))          
  121 CONTINUE                                                          
  122 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C          CODE FOR NROOTS=13                                           
C                                                                       
  130 DO 132 I = 1,IJ                                                   
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO 131 K = 1,MAX                                                  
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*                                
     *            ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )           
     *          +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401)           
     *          +   XIN(MX+ 4802)*YIN(MY+ 4802)*ZIN(MZ+ 4802)           
     *          +   XIN(MX+ 7203)*YIN(MY+ 7203)*ZIN(MZ+ 7203)           
     *          +   XIN(MX+ 9604)*YIN(MY+ 9604)*ZIN(MZ+ 9604)           
     *          +   XIN(MX+12005)*YIN(MY+12005)*ZIN(MZ+12005)           
     *          +   XIN(MX+14406)*YIN(MY+14406)*ZIN(MZ+14406)           
     *          +   XIN(MX+16807)*YIN(MY+16807)*ZIN(MZ+16807)           
     *          +   XIN(MX+19208)*YIN(MY+19208)*ZIN(MZ+19208)           
     *          +   XIN(MX+21609)*YIN(MY+21609)*ZIN(MZ+21609)           
     *          +   XIN(MX+24010)*YIN(MY+24010)*ZIN(MZ+24010)           
     *          +   XIN(MX+26411)*YIN(MY+26411)*ZIN(MZ+26411)           
     *          +   XIN(MX+28812)*YIN(MY+28812)*ZIN(MZ+28812))          
  131 CONTINUE                                                          
  132 CONTINUE                                                          
      RETURN                                                            
      END                                                               
C*MODULE INT2A   *DECK FORMS                                            
      SUBROUTINE FORMS(GHONDO)                                          
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      DIMENSION GHONDO(*)                                               
C                                                                       
      COMMON /DENS  / DKL(784),DIJ(784)                                 
C$omp threadprivate(/DENS/)
      COMMON /INTDEX/ IJGT(784),IJX(784),IJY(784),IJZ(784),IK(784),     
     *                KLGT(784),KLX(784),KLY(784),KLZ(784)              
C$omp threadprivate(/INTDEX/)
      COMMON /SHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,          
     *                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,          
     *                NIJ,IJ,KL,IJKL                                    
C$omp threadprivate(/SHLNOS/)
      COMMON /ROOT  / XX,U(13),W(13),NROOTS                             
C$omp threadprivate(/ROOT/)
      COMMON /XYZ   / XIN(31213),YIN(31213),ZIN(31213)                  
C$omp threadprivate(/XYZ/)
C                                                                       
C     ----- FORM INTEGRALS OVER FUNCTIONS -----                         
C     DIMENSIONING XIN(81,5), AND ROLLING UP THE COMPUTATION            
C     OF GHONDO IN A LOOP OF LENGTH NROOTS ADDS 33 SECONDS TO           
C     A 240 SECOND INTEGRAL COMPUTATION JOB.  LEAVE IT UNROLLED.        
C                                                                       
      GO TO (10,20,30,40,50,60,70,80,90,100,110,120,130),NROOTS         
C                                                                       
C          CODE FOR NROOTS=1                                            
C                                                                       
   10 CONTINUE                                                          
      DO 12 I = 1,IJ                                                    
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO 11 K = 1,MAX                                                   
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*                                
     *            ( XIN(MX      )*YIN(MY      )*ZIN(MZ      ))          
   11 CONTINUE                                                          
   12 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C          CODE FOR NROOTS=2                                            
C                                                                       
   20 CONTINUE                                                          
      DO 22 I = 1,IJ                                                    
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO 21 K = 1,MAX                                                   
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*                                
     *            ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )           
     *          +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401))          
   21 CONTINUE                                                          
   22 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C          CODE FOR NROOTS=3                                            
C                                                                       
   30 CONTINUE                                                          
      DO 32 I = 1,IJ                                                    
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO 31 K = 1,MAX                                                   
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*                                
     *            ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )           
     *          +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401)           
     *          +   XIN(MX+ 4802)*YIN(MY+ 4802)*ZIN(MZ+ 4802))          
   31 CONTINUE                                                          
   32 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C          CODE FOR NROOTS=4                                            
C                                                                       
   40 CONTINUE                                                          
      DO 42 I = 1,IJ                                                    
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO 41 K = 1,MAX                                                   
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*                                
     *            ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )           
     *          +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401)           
     *          +   XIN(MX+ 4802)*YIN(MY+ 4802)*ZIN(MZ+ 4802)           
     *          +   XIN(MX+ 7203)*YIN(MY+ 7203)*ZIN(MZ+ 7203))          
   41 CONTINUE                                                          
   42 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C          CODE FOR NROOTS=5                                            
C                                                                       
   50 CONTINUE                                                          
      DO 52 I = 1,IJ                                                    
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO 51 K = 1,MAX                                                   
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*                                
     *            ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )           
     *          +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401)           
     *          +   XIN(MX+ 4802)*YIN(MY+ 4802)*ZIN(MZ+ 4802)           
     *          +   XIN(MX+ 7203)*YIN(MY+ 7203)*ZIN(MZ+ 7203)           
     *          +   XIN(MX+ 9604)*YIN(MY+ 9604)*ZIN(MZ+ 9604))          
   51 CONTINUE                                                          
   52 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C          CODE FOR NROOTS=6                                            
C                                                                       
   60 CONTINUE                                                          
      DO 62 I = 1,IJ                                                    
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO 61 K = 1,MAX                                                   
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*                                
     *            ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )           
     *          +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401)           
     *          +   XIN(MX+ 4802)*YIN(MY+ 4802)*ZIN(MZ+ 4802)           
     *          +   XIN(MX+ 7203)*YIN(MY+ 7203)*ZIN(MZ+ 7203)           
     *          +   XIN(MX+ 9604)*YIN(MY+ 9604)*ZIN(MZ+ 9604)           
     *          +   XIN(MX+12005)*YIN(MY+12005)*ZIN(MZ+12005))          
   61 CONTINUE                                                          
   62 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C          CODE FOR NROOTS=7                                            
C                                                                       
   70 CONTINUE                                                          
      DO 72 I = 1,IJ                                                    
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO 71 K = 1,MAX                                                   
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*                                
     *            ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )           
     *          +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401)           
     *          +   XIN(MX+ 4802)*YIN(MY+ 4802)*ZIN(MZ+ 4802)           
     *          +   XIN(MX+ 7203)*YIN(MY+ 7203)*ZIN(MZ+ 7203)           
     *          +   XIN(MX+ 9604)*YIN(MY+ 9604)*ZIN(MZ+ 9604)           
     *          +   XIN(MX+12005)*YIN(MY+12005)*ZIN(MZ+12005)           
     *          +   XIN(MX+14406)*YIN(MY+14406)*ZIN(MZ+14406))          
   71 CONTINUE                                                          
   72 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C          CODE FOR NROOTS=8                                            
C                                                                       
   80 CONTINUE                                                          
      DO 82 I = 1,IJ                                                    
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO 81 K = 1,MAX                                                   
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*                                
     *            ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )           
     *          +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401)           
     *          +   XIN(MX+ 4802)*YIN(MY+ 4802)*ZIN(MZ+ 4802)           
     *          +   XIN(MX+ 7203)*YIN(MY+ 7203)*ZIN(MZ+ 7203)           
     *          +   XIN(MX+ 9604)*YIN(MY+ 9604)*ZIN(MZ+ 9604)           
     *          +   XIN(MX+12005)*YIN(MY+12005)*ZIN(MZ+12005)           
     *          +   XIN(MX+14406)*YIN(MY+14406)*ZIN(MZ+14406)           
     *          +   XIN(MX+16807)*YIN(MY+16807)*ZIN(MZ+16807))          
   81 CONTINUE                                                          
   82 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C          CODE FOR NROOTS=9                                            
C                                                                       
   90 CONTINUE                                                          
      DO 92 I = 1,IJ                                                    
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO 91 K = 1,MAX                                                   
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*                                
     *            ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )           
     *          +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401)           
     *          +   XIN(MX+ 4802)*YIN(MY+ 4802)*ZIN(MZ+ 4802)           
     *          +   XIN(MX+ 7203)*YIN(MY+ 7203)*ZIN(MZ+ 7203)           
     *          +   XIN(MX+ 9604)*YIN(MY+ 9604)*ZIN(MZ+ 9604)           
     *          +   XIN(MX+12005)*YIN(MY+12005)*ZIN(MZ+12005)           
     *          +   XIN(MX+14406)*YIN(MY+14406)*ZIN(MZ+14406)           
     *          +   XIN(MX+16807)*YIN(MY+16807)*ZIN(MZ+16807)           
     *          +   XIN(MX+19208)*YIN(MY+19208)*ZIN(MZ+19208))          
   91 CONTINUE                                                          
   92 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C          CODE FOR NROOTS=10                                           
C                                                                       
  100 CONTINUE                                                          
      DO 102 I = 1,IJ                                                   
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO 101 K = 1,MAX                                                  
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*                                
     *            ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )           
     *          +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401)           
     *          +   XIN(MX+ 4802)*YIN(MY+ 4802)*ZIN(MZ+ 4802)           
     *          +   XIN(MX+ 7203)*YIN(MY+ 7203)*ZIN(MZ+ 7203)           
     *          +   XIN(MX+ 9604)*YIN(MY+ 9604)*ZIN(MZ+ 9604)           
     *          +   XIN(MX+12005)*YIN(MY+12005)*ZIN(MZ+12005)           
     *          +   XIN(MX+14406)*YIN(MY+14406)*ZIN(MZ+14406)           
     *          +   XIN(MX+16807)*YIN(MY+16807)*ZIN(MZ+16807)           
     *          +   XIN(MX+19208)*YIN(MY+19208)*ZIN(MZ+19208)           
     *          +   XIN(MX+21609)*YIN(MY+21609)*ZIN(MZ+21609))          
  101 CONTINUE                                                          
  102 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C          CODE FOR NROOTS=11                                           
C                                                                       
  110 CONTINUE                                                          
      DO 112 I = 1,IJ                                                   
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO 111 K = 1,MAX                                                  
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*                                
     *            ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )           
     *          +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401)           
     *          +   XIN(MX+ 4802)*YIN(MY+ 4802)*ZIN(MZ+ 4802)           
     *          +   XIN(MX+ 7203)*YIN(MY+ 7203)*ZIN(MZ+ 7203)           
     *          +   XIN(MX+ 9604)*YIN(MY+ 9604)*ZIN(MZ+ 9604)           
     *          +   XIN(MX+12005)*YIN(MY+12005)*ZIN(MZ+12005)           
     *          +   XIN(MX+14406)*YIN(MY+14406)*ZIN(MZ+14406)           
     *          +   XIN(MX+16807)*YIN(MY+16807)*ZIN(MZ+16807)           
     *          +   XIN(MX+19208)*YIN(MY+19208)*ZIN(MZ+19208)           
     *          +   XIN(MX+21609)*YIN(MY+21609)*ZIN(MZ+21609)           
     *          +   XIN(MX+24010)*YIN(MY+24010)*ZIN(MZ+24010))          
  111 CONTINUE                                                          
  112 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C          CODE FOR NROOTS=12                                           
C                                                                       
  120 CONTINUE                                                          
      DO 122 I = 1,IJ                                                   
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO 121 K = 1,MAX                                                  
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*                                
     *            ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )           
     *          +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401)           
     *          +   XIN(MX+ 4802)*YIN(MY+ 4802)*ZIN(MZ+ 4802)           
     *          +   XIN(MX+ 7203)*YIN(MY+ 7203)*ZIN(MZ+ 7203)           
     *          +   XIN(MX+ 9604)*YIN(MY+ 9604)*ZIN(MZ+ 9604)           
     *          +   XIN(MX+12005)*YIN(MY+12005)*ZIN(MZ+12005)           
     *          +   XIN(MX+14406)*YIN(MY+14406)*ZIN(MZ+14406)           
     *          +   XIN(MX+16807)*YIN(MY+16807)*ZIN(MZ+16807)           
     *          +   XIN(MX+19208)*YIN(MY+19208)*ZIN(MZ+19208)           
     *          +   XIN(MX+21609)*YIN(MY+21609)*ZIN(MZ+21609)           
     *          +   XIN(MX+24010)*YIN(MY+24010)*ZIN(MZ+24010)           
     *          +   XIN(MX+26411)*YIN(MY+26411)*ZIN(MZ+26411))          
  121 CONTINUE                                                          
  122 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C          CODE FOR NROOTS=13                                           
C                                                                       
  130 DO 132 I = 1,IJ                                                   
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO 131 K = 1,MAX                                                  
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*                                
     *            ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )           
     *          +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401)           
     *          +   XIN(MX+ 4802)*YIN(MY+ 4802)*ZIN(MZ+ 4802)           
     *          +   XIN(MX+ 7203)*YIN(MY+ 7203)*ZIN(MZ+ 7203)           
     *          +   XIN(MX+ 9604)*YIN(MY+ 9604)*ZIN(MZ+ 9604)           
     *          +   XIN(MX+12005)*YIN(MY+12005)*ZIN(MZ+12005)           
     *          +   XIN(MX+14406)*YIN(MY+14406)*ZIN(MZ+14406)           
     *          +   XIN(MX+16807)*YIN(MY+16807)*ZIN(MZ+16807)           
     *          +   XIN(MX+19208)*YIN(MY+19208)*ZIN(MZ+19208)           
     *          +   XIN(MX+21609)*YIN(MY+21609)*ZIN(MZ+21609)           
     *          +   XIN(MX+24010)*YIN(MY+24010)*ZIN(MZ+24010)           
     *          +   XIN(MX+26411)*YIN(MY+26411)*ZIN(MZ+26411)           
     *          +   XIN(MX+28812)*YIN(MY+28812)*ZIN(MZ+28812))          
  131 CONTINUE                                                          
  132 CONTINUE                                                          
      RETURN                                                            
      END                                                               
C*MODULE INT2A   *DECK GENRAL                                           
      SUBROUTINE GENRAL(GHONDO,DDIJ)                                    
      USE lrcdft, ONLY: EMU2                                            
      use mx_limits, only: mxgsh,mxg2                                   
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      DIMENSION GHONDO(*),DDIJ(*)                                       
C                                                                       
      LOGICAL IANDJ,KANDL,SAME,OUT,NORM,DOUBLE                          
      LOGICAL LRINT                                                     
C                                                                       
C                                                                       
      COMMON /DENS  / DKL(784),DIJ(784)                                 
C$omp threadprivate(/DENS/)
      COMMON /IJGNRL/ AA(MXG2),R(MXG2),X1(MXG2),Y1(MXG2),Z1(MXG2),      
     *                IJD(784)                                          
C$omp threadprivate(/IJGNRL/)
      COMMON /MISC  / IANDJ,KANDL,SAME                                  
C$omp threadprivate(/MISC/)
      COMMON /NLRCF / LRINT                                             
C$omp threadprivate(/NLRCF /)
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK                 
      COMMON /ROOT  / XX,U(13),W(13),NROOTS                             
C$omp threadprivate(/ROOT/)
      COMMON /SETINT/ IN(13),KN(13),NI,NJ,NK,NL,NMAX,MMAX,              
     +                BP01,B00,B10,XCP00,XC00,YCP00,YC00,ZCP00,ZC00,F00,
     +                DXIJ,DYIJ,DZIJ,DXKL,DYKL,DZKL                     
C$omp threadprivate(/SETINT/)
      COMMON /SHLINF/  AG(MXGSH),CSA(MXGSH),CPA(MXGSH),CDA(MXGSH),      
     *                CFA(MXGSH),CGA(MXGSH),CHA(MXGSH),CIA(MXGSH),      
     *                 BG(MXGSH),CSB(MXGSH),CPB(MXGSH),CDB(MXGSH),      
     *                CFB(MXGSH),CGB(MXGSH),CHB(MXGSH),CIB(MXGSH),      
     *                 CG(MXGSH),CSC(MXGSH),CPC(MXGSH),CDC(MXGSH),      
     *                CFC(MXGSH),CGC(MXGSH),CHC(MXGSH),CIC(MXGSH),      
     *                 DG(MXGSH),CSD(MXGSH),CPD(MXGSH),CDD(MXGSH),      
     *                CFD(MXGSH),CGD(MXGSH),CHD(MXGSH),CID(MXGSH),      
     *                XI,YI,ZI,XJ,YJ,ZJ,RRI,XK,YK,ZK,XL,YL,ZL,RRK,      
     *                NGA,NGB,NGC,NGD                                   
C$omp threadprivate(/SHLINF/)
      COMMON /SHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,          
     +                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,          
     +                NIJ,IJ,KL,IJKL                                    
C$omp threadprivate(/SHLNOS/)
      COMMON /SHLT  / TOL,CUTOFF,ICOUNT,OUT                             
C                                                                       
      DIMENSION IN1(13)                                                 
C                                                                       
      PARAMETER (SQRT3=1.73205080756888D+00, SQRT5=2.23606797749979D+00,
     *           SQRT7=2.64575131106459D+00, PI252=34.986836655250D+00, 
     *           SQRT9=3.0D+00,SQRT11=3.3166247903553998D+00,           
     *           ZERO=0.0D+00, HALF=0.5D+00, ONE=1.0D+00)               
C                                                                       
C     GENERAL INTEGRAL ROUTINE FOR SPDFGHI AND L FUNCTIONS              
C                                                                       
      FACTOR = PI252*QQ4                                                
      NORM = NORMF .NE. 1 .OR. NORMP .NE. 1                             
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
C                                                                       
C     ----- K PRIMITIVE                                                 
C                                                                       
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
         CFK = CFC(KG)*FACTOR                                           
         CGK = CGC(KG)*FACTOR                                           
         CHK = CHC(KG)*FACTOR                                           
         CIK = CIC(KG)*FACTOR                                           
C                                                                       
C        ----- L PRIMITIVE                                              
C                                                                       
         IF (KANDL) LGMAX = KG                                          
         DO 460 LG = 1,LGMAX                                            
            AL = DG(LG)                                                 
            B = AK+AL                                                   
            BINV = ONE/B                                                
            BBRRK = AL*BRRK*BINV                                        
            IF (BBRRK .GT. TOL) GO TO 460                               
            CSL = CSD(LG)                                               
            CPL = CPD(LG)                                               
            CDL = CDD(LG)                                               
            CFL = CFD(LG)                                               
            CGL = CGD(LG)                                               
            CHL = CHD(LG)                                               
            CIL = CID(LG)                                               
            XB = (AKXK+AL*XL)*BINV                                      
            YB = (AKYK+AL*YL)*BINV                                      
            ZB = (AKZK+AL*ZL)*BINV                                      
            BXBK = B*(XB-XK)                                            
            BYBK = B*(YB-YK)                                            
            BZBK = B*(ZB-ZK)                                            
            BXBI = B*(XB-XI)                                            
            BYBI = B*(YB-YI)                                            
            BZBI = B*(ZB-ZI)                                            
C                                                                       
C           ----- DENSITY FACTOR                                        
C                                                                       
            DOUBLE=KANDL.AND.KG.NE.LG                                   
            N = 0                                                       
            MAX = MAXL                                                  
            DUM1 = ZERO                                                 
            DUM2 = ZERO                                                 
            DO 370 K = MINK,MAXK                                        
               GO TO (140,160,220,220,180,220,220,200,220,220,          
     1                201,220,220,202,220,220,220,220,220,203,          
     1                204,220,220,205,220,220,220,220,220,206,          
     1                220,220,207,220,220,                              
     1                208,220,220,209,220,220,220,220,220,210,          
     1                220,220,220,220,220,211,220,220,212,220,          
     1                220,                                              
     1                213,220,220,214,220,220,220,220,220,215,          
     1                220,220,220,220,220,216,220,220,217,220,          
     1                220,218,220,220,220,220,220,219),K                
  140          DUM1 = CSK*BINV                                          
               GO TO 220                                                
  160          DUM1 = CPK*BINV                                          
               GO TO 220                                                
  180          DUM1 = CDK*BINV                                          
               GO TO 220                                                
  200          IF (NORM) DUM1 = DUM1*SQRT3                              
               GO TO 220                                                
  201          DUM1 = CFK*BINV                                          
               GO TO 220                                                
  202          IF (NORM) DUM1 = DUM1*SQRT5                              
               GO TO 220                                                
  203          IF (NORM) DUM1 = DUM1*SQRT3                              
               GO TO 220                                                
  204          DUM1 = CGK*BINV                                          
               GO TO 220                                                
  205          IF (NORM) DUM1 = DUM1*SQRT7                              
               GO TO 220                                                
  206          IF (NORM) DUM1 = DUM1*SQRT5/SQRT3                        
               GO TO 220                                                
  207          IF (NORM) DUM1 = DUM1*SQRT3                              
               GO TO 220                                                
  208          DUM1 = CHK*BINV                                          
               GO TO 220                                                
  209          IF (NORM) DUM1 = DUM1*SQRT9                              
               GO TO 220                                                
  210          IF (NORM) DUM1 = DUM1*SQRT7/SQRT3                        
               GO TO 220                                                
  211          IF (NORM) DUM1 = DUM1*SQRT3                              
               GO TO 220                                                
  212          IF (NORM) DUM1 = DUM1*SQRT5/SQRT3                        
               GO TO 220                                                
  213          DUM1 = CIK*BINV                                          
               GO TO 220                                                
  214          IF (NORM) DUM1 = DUM1*SQRT11                             
               GO TO 220                                                
  215          IF (NORM) DUM1 = DUM1*SQRT3                              
               GO TO 220                                                
  216          IF (NORM) DUM1 = DUM1*SQRT3                              
               GO TO 220                                                
  217          IF (NORM) DUM1 = DUM1*SQRT7/(SQRT5*SQRT3)                
               GO TO 220                                                
  218          IF (NORM) DUM1 = DUM1*SQRT5                              
               GO TO 220                                                
  219          IF (NORM) DUM1 = DUM1*SQRT5/SQRT3                        
C                                                                       
  220          IF (KANDL) MAX = K                                       
               DO 360 L = MINL,MAX                                      
                  GO TO (240,280,340,340,300,340,340,320,340,340,       
     1                   321,340,340,322,340,340,340,340,340,323,       
     1                   324,340,340,325,340,340,340,340,340,326,       
     1                   340,340,327,340,340,                           
     1                   328,340,340,329,340,340,340,340,340,330,       
     1                   340,340,340,340,340,331,340,340,332,340,       
     1                   340,                                           
     1                   333,340,340,334,340,340,340,340,340,335,       
     1                   340,340,340,340,340,336,340,340,337,340,       
     1                   340,338,340,340,340,340,340,339),L             
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
  320             IF (NORM) DUM2 = DUM2*SQRT3                           
                  GO TO 340                                             
  321             DUM2 = DUM1*CFL                                       
                  IF (DOUBLE) DUM2 = DUM2+DUM2                          
                  GO TO 340                                             
  322             IF (NORM) DUM2 = DUM2*SQRT5                           
                  GO TO 340                                             
  323             IF (NORM) DUM2 = DUM2*SQRT3                           
                  GO TO 340                                             
  324             DUM2 = DUM1*CGL                                       
                  IF (DOUBLE) DUM2 = DUM2+DUM2                          
                  GO TO 340                                             
  325             IF (NORM) DUM2 = DUM2*SQRT7                           
                  GO TO 340                                             
  326             IF (NORM) DUM2 = DUM2*SQRT5/SQRT3                     
                  GO TO 340                                             
  327             IF (NORM) DUM2 = DUM2*SQRT3                           
                  GO TO 340                                             
  328             DUM2 = DUM1*CHL                                       
                  IF (DOUBLE) DUM2 = DUM2+DUM2                          
                  GO TO 340                                             
  329             IF (NORM) DUM2 = DUM2*SQRT9                           
                  GO TO 340                                             
  330             IF (NORM) DUM2 = DUM2*SQRT7/SQRT3                     
                  GO TO 340                                             
  331             IF (NORM) DUM2 = DUM2*SQRT3                           
                  GO TO 340                                             
  332             IF (NORM) DUM2 = DUM2*SQRT5/SQRT3                     
                  GO TO 340                                             
  333             DUM2 = DUM1*CIL                                       
                  IF (DOUBLE) DUM2 = DUM2+DUM2                          
                  GO TO 340                                             
  334             IF (NORM) DUM2 = DUM2*SQRT11                          
                  GO TO 340                                             
  335             IF (NORM) DUM2 = DUM2*SQRT3                           
                  GO TO 340                                             
  336             IF (NORM) DUM2 = DUM2*SQRT3                           
                  GO TO 340                                             
  337             IF (NORM) DUM2 = DUM2*SQRT7/(SQRT5*SQRT3)             
                  GO TO 340                                             
  338             IF (NORM) DUM2 = DUM2*SQRT5                           
                  GO TO 340                                             
  339             IF (NORM) DUM2 = DUM2*SQRT5/SQRT3                     
C                                                                       
  340             N = N+1                                               
                  DKL(N) = DUM2                                         
  360          CONTINUE                                                 
  370       CONTINUE                                                    
C                                                                       
C           ----- PAIR OF I,J PRIMITIVES                                
C                                                                       
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
               IF(LRINT) THEN                                           
                 RHO0 = RHO                                             
                 RHO  = RHO0*EMU2/(RHO0+EMU2)                           
               ENDIF                                                    
               XA = X1(N)                                               
               YA = Y1(N)                                               
               ZA = Z1(N)                                               
               XX = RHO*((XA-XB)*(XA-XB) + (YA-YB)*(YA-YB)              
     *                                   + (ZA-ZB)*(ZA-ZB))             
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
C                                                                       
C              ----- ROOTS AND WEIGHTS FOR QUADRATURE                   
C                                                                       
               !IF (NROOTS .LE. 3) CALL RT123                           
               !IF (NROOTS .EQ. 4) CALL ROOT4                           
               !IF (NROOTS .EQ. 5) CALL ROOT5                           
               !IF (NROOTS .GE. 6) CALL ROOT6                           
                                                                        
               MM = 0                                                   
               MAX = NMAX+1                                             
C                                                                       
C              COMPUTE TWO-ELECTRON INTEGRALS FOR EACH ROOT             
C                                                                       
               DO 420 M = 1,NROOTS                                      
                  U2 = U(M)*RHO                                         
                  F00 = EXPE*W(M)                                       
                  IF(LRINT) F00 = F00*SQRT(EMU2/(RHO0+EMU2))            
                  DO 400 I = 1,MAX                                      
                     IN(I) = IN1(I)+MM                                  
  400             CONTINUE                                              
                  IF(.NOT.LRINT)THEN                                    
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
                  ELSE                                                  
                     T2    = U2/(U2+RHO)                                
                     T2AR  = T2*RHO/A                                   
                     T2BR  = T2*RHO/B                                   
                     BP01  = HALF/B*(ONE-T2BR)                          
                     B00   = HALF/AANDB*T2*RHO/RHO0                     
                     B10   = HALF/A*(ONE-T2AR)                          
                     XCP00 = (XB-XK)+T2BR*(XA-XB)                       
                     XC00  = (XA-XI)-T2AR*(XA-XB)                       
                     YCP00 = (YB-YK)+T2BR*(YA-YB)                       
                     YC00  = (YA-YI)-T2AR*(YA-YB)                       
                     ZCP00 = (ZB-ZK)+T2BR*(ZA-ZB)                       
                     ZC00  = (ZA-ZI)-T2AR*(ZA-ZB)                       
                  END IF                                                
                  CALL XYZINT                                           
                  MM = MM+2401                                          
  420          CONTINUE                                                 
C                                                                       
C              ----- FORM (I,J//K,L) INTEGRALS OVER FUNCTIONS           
C                                                                       
               CALL FORMS(GHONDO)                                       
  440       NN = NN+49                                                  
  460    CONTINUE                                                       
  480 CONTINUE                                                          
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2A   *DECK IJPRIM                                           
      SUBROUTINE IJPRIM(DDIJ)                                           
      use mx_limits, only: mxgsh,mxg2                                   
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      LOGICAL IANDJ,KANDL,SAME,OUT,NORM                                 
C                                                                       
C                                                                       
      DIMENSION DDIJ(49*MXG2)                                           
C                                                                       
      COMMON /IJGNRL/ A(MXG2),R(MXG2),X1(MXG2),Y1(MXG2),Z1(MXG2),       
     *                IJD(784)                                          
C$omp threadprivate(/IJGNRL/)
      COMMON /MISC  / IANDJ,KANDL,SAME                                  
C$omp threadprivate(/MISC/)
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK                 
      COMMON /SHLINF/  AG(MXGSH),CSA(MXGSH),CPA(MXGSH),CDA(MXGSH),      
     *                CFA(MXGSH),CGA(MXGSH),CHA(MXGSH),CIA(MXGSH),      
     *                 BG(MXGSH),CSB(MXGSH),CPB(MXGSH),CDB(MXGSH),      
     *                CFB(MXGSH),CGB(MXGSH),CHB(MXGSH),CIB(MXGSH),      
     *                 CG(MXGSH),CSC(MXGSH),CPC(MXGSH),CDC(MXGSH),      
     *                CFC(MXGSH),CGC(MXGSH),CHC(MXGSH),CIC(MXGSH),      
     *                 DG(MXGSH),CSD(MXGSH),CPD(MXGSH),CDD(MXGSH),      
     *                CFD(MXGSH),CGD(MXGSH),CHD(MXGSH),CID(MXGSH),      
     *                XI,YI,ZI,XJ,YJ,ZJ,RRI,XK,YK,ZK,XL,YL,ZL,RRK,      
     *                NGA,NGB,NGC,NGD                                   
C$omp threadprivate(/SHLINF/)
      COMMON /SHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,          
     *                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,          
     *                NIJ,IJ,KL,IJKL                                    
C$omp threadprivate(/SHLNOS/)
      COMMON /SHLT  / TOL,CUTOFF,ICOUNT,OUT                             
C                                                                       
      PARAMETER (SQRT3=1.73205080756888D+00)                            
      PARAMETER (SQRT5=2.23606797749979D+00)                            
      PARAMETER (SQRT7=2.64575131106459D+00)                            
      PARAMETER (SQRT9=3.0D+00)                                         
      PARAMETER (SQRT11=3.3166247903553998D+00)                         
      PARAMETER (ZERO=0.0D+00)                                          
      PARAMETER (ONE=1.0D+00)                                           
C                                                                       
      NORM = NORMF .NE. 1 .OR. NORMP .NE. 1                             
      MAX = MAXJ                                                        
      N = 0                                                             
      NN = 0                                                            
      NM = -2**20                                                       
      DO 180 I = MINI,MAXI                                              
         GO TO (100,100,120,120,100,120,120,100,120,120,                
     1          100,120,120,100,120,120,120,120,120,100,                
     1          100,120,120,100,120,120,120,120,120,100,                
     1          120,120,100,120,120,                                    
     1          100,120,120,100,120,120,120,120,120,100,                
     1          120,120,120,120,120,100,120,120,100,120,                
     1          120,                                                    
     1          100,120,120,100,120,120,120,120,120,100,                
     1          120,120,120,120,120,100,120,120,100,120,                
     1          120,100,120,120,120,120,120,100),I                      
  100    NM = NN                                                        
  120    NN = NM                                                        
         IF (IANDJ) MAX = I                                             
         DO 170 J = MINJ,MAX                                            
            GO TO (140,140,160,160,140,160,160,140,160,160,             
     1             140,160,160,140,160,160,160,160,160,140,             
     1             140,160,160,140,160,160,160,160,160,140,             
     1             160,160,140,160,160,                                 
     1             140,160,160,140,160,160,160,160,160,140,             
     1             160,160,160,160,160,140,160,160,140,160,             
     1             160,                                                 
     1             140,160,160,140,160,160,160,160,160,140,             
     1             160,160,160,160,160,140,160,160,140,160,             
     1             160,140,160,160,160,160,160,140),J                   
  140       NN = NN+1                                                   
  160       N = N+1                                                     
            IJD(N) = NN                                                 
  170    CONTINUE                                                       
  180 CONTINUE                                                          
C                                                                       
C     ----- I PRIMITIVE                                                 
C                                                                       
      NIJ = 0                                                           
      JBMAX = NGB                                                       
      DO 540 IA = 1,NGA                                                 
         AI = AG(IA)                                                    
         ARRI = AI*RRI                                                  
         AXI = AI*XI                                                    
         AYI = AI*YI                                                    
         AZI = AI*ZI                                                    
         CSI = CSA(IA)                                                  
         CPI = CPA(IA)                                                  
         CDI = CDA(IA)                                                  
         CFI = CFA(IA)                                                  
         CGI = CGA(IA)                                                  
         CHI = CHA(IA)                                                  
         CII = CIA(IA)                                                  
C                                                                       
C        ----- J PRIMITIVE                                              
C                                                                       
         IF (IANDJ) JBMAX = IA                                          
         DO 520 JB = 1,JBMAX                                            
            AJ = BG(JB)                                                 
            AA = AI+AJ                                                  
            AAINV = ONE/AA                                              
            DUM = AJ*ARRI*AAINV                                         
            IF (DUM .GT. TOL) GO TO 520                                 
            CSJ = CSB(JB)                                               
            CPJ = CPB(JB)                                               
            CDJ = CDB(JB)                                               
            CFJ = CFB(JB)                                               
            CGJ = CGB(JB)                                               
            CHJ = CHB(JB)                                               
            CIJ = CIB(JB)                                               
            NM = 49*NIJ                                                 
            NN = NM                                                     
            NIJ = NIJ+1                                                 
            R(NIJ) = DUM                                                
            A(NIJ) = AA                                                 
            X1(NIJ) = (AXI+AJ*XJ)*AAINV                                 
            Y1(NIJ) = (AYI+AJ*YJ)*AAINV                                 
            Z1(NIJ) = (AZI+AJ*ZJ)*AAINV                                 
C                                                                       
C           ----- DENSITY FACTOR                                        
C                                                                       
            DUM1 = ZERO                                                 
            DUM2 = ZERO                                                 
            DO 420 I = MINI,MAXI                                        
               GO TO (200,220,420,420,240,420,420,260,420,420,          
     1                261,420,420,262,420,420,420,420,420,263,          
     1                264,420,420,265,420,420,420,420,420,266,          
     1                420,420,267,420,420,                              
     1                268,420,420,269,420,420,420,420,420,270,          
     1                420,420,420,420,420,271,420,420,272,420,          
     1                420,                                              
     1                273,420,420,274,420,420,420,420,420,275,          
     1                420,420,420,420,420,276,420,420,277,420,          
     1                420,278,420,420,420,420,420,279),I                
  200          DUM1 = CSI*AAINV                                         
               GO TO 280                                                
  220          DUM1 = CPI*AAINV                                         
               GO TO 280                                                
  240          DUM1 = CDI*AAINV                                         
               GO TO 280                                                
  260          IF (NORM) DUM1 = DUM1*SQRT3                              
               GO TO 280                                                
  261          DUM1 = CFI*AAINV                                         
               GO TO 280                                                
  262          IF (NORM) DUM1 = DUM1*SQRT5                              
               GO TO 280                                                
  263          IF (NORM) DUM1 = DUM1*SQRT3                              
               GO TO 280                                                
  264          DUM1 = CGI*AAINV                                         
               GO TO 280                                                
  265          IF (NORM) DUM1 = DUM1*SQRT7                              
               GO TO 280                                                
  266          IF (NORM) DUM1 = DUM1*SQRT5/SQRT3                        
               GO TO 280                                                
  267          IF (NORM) DUM1 = DUM1*SQRT3                              
               GO TO 280                                                
  268          DUM1 = CHI*AAINV                                         
               GO TO 280                                                
  269          IF (NORM) DUM1 = DUM1*SQRT9                              
               GO TO 280                                                
  270          IF (NORM) DUM1 = DUM1*SQRT7/SQRT3                        
               GO TO 280                                                
  271          IF (NORM) DUM1 = DUM1*SQRT3                              
               GO TO 280                                                
  272          IF (NORM) DUM1 = DUM1*SQRT5/SQRT3                        
               GO TO 280                                                
  273          DUM1 = CII*AAINV                                         
               GO TO 280                                                
  274          IF (NORM) DUM1 = DUM1*SQRT11                             
               GO TO 280                                                
  275          IF (NORM) DUM1 = DUM1*SQRT3                              
               GO TO 280                                                
  276          IF (NORM) DUM1 = DUM1*SQRT3                              
               GO TO 280                                                
  277          IF (NORM) DUM1 = DUM1*SQRT7/(SQRT5*SQRT3)                
               GO TO 280                                                
  278          IF (NORM) DUM1 = DUM1*SQRT5                              
               GO TO 280                                                
  279          IF (NORM) DUM1 = DUM1*SQRT5/SQRT3                        
C                                                                       
  280          IF (IANDJ) MAX = I                                       
               DO 400 J = MINJ,MAX                                      
                  GO TO (300,320,400,400,340,400,400,360,400,400,       
     1                   361,400,400,362,400,400,400,400,400,363,       
     1                   364,400,400,365,400,400,400,400,400,366,       
     1                   400,400,367,400,400,                           
     1                   368,400,400,369,400,400,400,400,400,370,       
     1                   400,400,400,400,400,371,400,400,372,400,       
     1                   400,                                           
     1                   373,400,400,374,400,400,400,400,400,375,       
     1                   400,400,400,400,400,376,400,400,377,400,       
     1                   400,378,400,400,400,400,400,379),J             
  300             DUM2 = DUM1*CSJ                                       
                  GO TO 380                                             
  320             DUM2 = DUM1*CPJ                                       
                  GO TO 380                                             
  340             DUM2 = DUM1*CDJ                                       
                  GO TO 380                                             
  360             IF (NORM) DUM2 = DUM2*SQRT3                           
                  GO TO 380                                             
  361             DUM2 = DUM1*CFJ                                       
                  GO TO 380                                             
  362             IF (NORM) DUM2 = DUM2*SQRT5                           
                  GO TO 380                                             
  363             IF (NORM) DUM2 = DUM2*SQRT3                           
                  GO TO 380                                             
  364             DUM2 = DUM1*CGJ                                       
                  GO TO 380                                             
  365             IF (NORM) DUM2 = DUM2*SQRT7                           
                  GO TO 380                                             
  366             IF (NORM) DUM2 = DUM2*SQRT5/SQRT3                     
                  GO TO 380                                             
  367             IF (NORM) DUM2 = DUM2*SQRT3                           
                  GO TO 380                                             
  368             DUM2 = DUM1*CHJ                                       
                  GO TO 380                                             
  369             IF (NORM) DUM2 = DUM2*SQRT9                           
                  GO TO 380                                             
  370             IF (NORM) DUM2 = DUM2*SQRT7/SQRT3                     
                  GO TO 380                                             
  371             IF (NORM) DUM2 = DUM2*SQRT3                           
                  GO TO 380                                             
  372             IF (NORM) DUM2 = DUM2*SQRT5/SQRT3                     
                  GO TO 380                                             
  373             DUM2 = DUM1*CIJ                                       
                  GO TO 380                                             
  374             IF (NORM) DUM2 = DUM2*SQRT11                          
                  GO TO 380                                             
  375             IF (NORM) DUM2 = DUM2*SQRT3                           
                  GO TO 380                                             
  376             IF (NORM) DUM2 = DUM2*SQRT3                           
                  GO TO 380                                             
  377             IF (NORM) DUM2 = DUM2*SQRT7/(SQRT5*SQRT3)             
                  GO TO 380                                             
  378             IF (NORM) DUM2 = DUM2*SQRT5                           
                  GO TO 380                                             
  379             IF (NORM) DUM2 = DUM2*SQRT5/SQRT3                     
C                                                                       
  380             NN = NN+1                                             
                  DDIJ(NN) = DUM2                                       
  400          CONTINUE                                                 
  420       CONTINUE                                                    
            IF ( .NOT. IANDJ) GO TO 520                                 
            IF (IA .EQ. JB) GO TO 520                                   
            GO TO (500,440,460,455,450,445,444),LIT                     
  440       IF (MINI .EQ. 2) GO TO 500                                  
            DDIJ(NM+2) = DDIJ(NM+2)+CSI*CPJ*AAINV                       
            GO TO 480                                                   
  444       DDIJ(NM+28) = DDIJ(NM+28)+DDIJ(NM+28)                       
            DDIJ(NM+27) = DDIJ(NM+27)+DDIJ(NM+27)                       
            DDIJ(NM+26) = DDIJ(NM+26)+DDIJ(NM+26)                       
            DDIJ(NM+25) = DDIJ(NM+25)+DDIJ(NM+25)                       
            DDIJ(NM+24) = DDIJ(NM+24)+DDIJ(NM+24)                       
            DDIJ(NM+23) = DDIJ(NM+23)+DDIJ(NM+23)                       
            DDIJ(NM+22) = DDIJ(NM+22)+DDIJ(NM+22)                       
            DDIJ(NM+21) = DDIJ(NM+21)+DDIJ(NM+21)                       
            DDIJ(NM+20) = DDIJ(NM+20)+DDIJ(NM+20)                       
            DDIJ(NM+19) = DDIJ(NM+19)+DDIJ(NM+19)                       
            DDIJ(NM+18) = DDIJ(NM+18)+DDIJ(NM+18)                       
            DDIJ(NM+17) = DDIJ(NM+17)+DDIJ(NM+17)                       
            DDIJ(NM+16) = DDIJ(NM+16)+DDIJ(NM+16)                       
  445       DDIJ(NM+15) = DDIJ(NM+15)+DDIJ(NM+15)                       
            DDIJ(NM+14) = DDIJ(NM+14)+DDIJ(NM+14)                       
            DDIJ(NM+13) = DDIJ(NM+13)+DDIJ(NM+13)                       
            DDIJ(NM+12) = DDIJ(NM+12)+DDIJ(NM+12)                       
            DDIJ(NM+11) = DDIJ(NM+11)+DDIJ(NM+11)                       
  450       DDIJ(NM+10) = DDIJ(NM+10)+DDIJ(NM+10)                       
            DDIJ(NM+9) = DDIJ(NM+9)+DDIJ(NM+9)                          
            DDIJ(NM+8) = DDIJ(NM+8)+DDIJ(NM+8)                          
            DDIJ(NM+7) = DDIJ(NM+7)+DDIJ(NM+7)                          
  455       DDIJ(NM+6) = DDIJ(NM+6)+DDIJ(NM+6)                          
            DDIJ(NM+5) = DDIJ(NM+5)+DDIJ(NM+5)                          
            DDIJ(NM+4) = DDIJ(NM+4)+DDIJ(NM+4)                          
  460       DDIJ(NM+2) = DDIJ(NM+2)+DDIJ(NM+2)                          
  480       DDIJ(NM+3) = DDIJ(NM+3)+DDIJ(NM+3)                          
  500       DDIJ(NM+1) = DDIJ(NM+1)+DDIJ(NM+1)                          
  520    CONTINUE                                                       
  540 CONTINUE                                                          
      RETURN                                                            
      END                                                               
C*MODULE INT2A   *DECK IJPRIM_gpu2                                      
      SUBROUTINE IJPRIM_gpu2(DDIJ)                                      
      use mx_limits, only: mxgsh,mxg2                                   
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      LOGICAL IANDJ,KANDL,SAME,OUT,NORM                                 
C                                                                       
C                                                                       
      DIMENSION DDIJ(49*MXG2)                                           
C                                                                       
      COMMON /IJGNRL/ A(MXG2),R(MXG2),X1(MXG2),Y1(MXG2),Z1(MXG2),       
     *                IJD(784)                                          
C$omp threadprivate(/IJGNRL/)
      COMMON /MISC  / IANDJ,KANDL,SAME                                  
C$omp threadprivate(/MISC/)
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK                 
      COMMON /SHLINF/  AG(MXGSH),CSA(MXGSH),CPA(MXGSH),CDA(MXGSH),      
     *                CFA(MXGSH),CGA(MXGSH),CHA(MXGSH),CIA(MXGSH),      
     *                 BG(MXGSH),CSB(MXGSH),CPB(MXGSH),CDB(MXGSH),      
     *                CFB(MXGSH),CGB(MXGSH),CHB(MXGSH),CIB(MXGSH),      
     *                 CG(MXGSH),CSC(MXGSH),CPC(MXGSH),CDC(MXGSH),      
     *                CFC(MXGSH),CGC(MXGSH),CHC(MXGSH),CIC(MXGSH),      
     *                 DG(MXGSH),CSD(MXGSH),CPD(MXGSH),CDD(MXGSH),      
     *                CFD(MXGSH),CGD(MXGSH),CHD(MXGSH),CID(MXGSH),      
     *                XI,YI,ZI,XJ,YJ,ZJ,RRI,XK,YK,ZK,XL,YL,ZL,RRK,      
     *                NGA,NGB,NGC,NGD                                   
C$omp threadprivate(/SHLINF/)
      COMMON /SHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,          
     *                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,          
     *                NIJ,IJ,KL,IJKL                                    
C$omp threadprivate(/SHLNOS/)
      COMMON /SHLT  / TOL,CUTOFF,ICOUNT,OUT                             
C                                                                       
      PARAMETER (SQRT3=1.73205080756888D+00)                            
      PARAMETER (SQRT5=2.23606797749979D+00)                            
      PARAMETER (SQRT7=2.64575131106459D+00)                            
      PARAMETER (SQRT9=3.0D+00)                                         
      PARAMETER (SQRT11=3.3166247903553998D+00)                         
      PARAMETER (ZERO=0.0D+00)                                          
      PARAMETER (ONE=1.0D+00)                                           
C                                                                       
      NORM = NORMF .NE. 1 .OR. NORMP .NE. 1                             
      MAX = MAXJ                                                        
      N = 0                                                             
      NN = 0                                                            
      NM = -2**20                                                       
      DO 180 I = MINI,MAXI                                              
         GO TO (100,100,120,120,100,120,120,100,120,120,                
     1          100,120,120,100,120,120,120,120,120,100,                
     1          100,120,120,100,120,120,120,120,120,100,                
     1          120,120,100,120,120,                                    
     1          100,120,120,100,120,120,120,120,120,100,                
     1          120,120,120,120,120,100,120,120,100,120,                
     1          120,                                                    
     1          100,120,120,100,120,120,120,120,120,100,                
     1          120,120,120,120,120,100,120,120,100,120,                
     1          120,100,120,120,120,120,120,100),I                      
  100    NM = NN                                                        
  120    NN = NM                                                        
         IF (IANDJ) MAX = I                                             
         DO 170 J = MINJ,MAX                                            
            GO TO (140,140,160,160,140,160,160,140,160,160,             
     1             140,160,160,140,160,160,160,160,160,140,             
     1             140,160,160,140,160,160,160,160,160,140,             
     1             160,160,140,160,160,                                 
     1             140,160,160,140,160,160,160,160,160,140,             
     1             160,160,160,160,160,140,160,160,140,160,             
     1             160,                                                 
     1             140,160,160,140,160,160,160,160,160,140,             
     1             160,160,160,160,160,140,160,160,140,160,             
     1             160,140,160,160,160,160,160,140),J                   
  140       NN = NN+1                                                   
  160       N = N+1                                                     
            IJD(N) = NN                                                 
  170    CONTINUE                                                       
  180 CONTINUE                                                          
C                                                                       
C     ----- I PRIMITIVE                                                 
C                                                                       
      NIJ = 0                                                           
      JBMAX = NGB                                                       
      DO 540 IA = 1,NGA                                                 
         AI = AG(IA)                                                    
         ARRI = AI*RRI                                                  
         AXI = AI*XI                                                    
         AYI = AI*YI                                                    
         AZI = AI*ZI                                                    
         CSI = CSA(IA)                                                  
         CPI = CPA(IA)                                                  
         CDI = CDA(IA)                                                  
         CFI = CFA(IA)                                                  
         CGI = CGA(IA)                                                  
         CHI = CHA(IA)                                                  
         CII = CIA(IA)                                                  
C                                                                       
C        ----- J PRIMITIVE                                              
C                                                                       
         IF (IANDJ) JBMAX = IA                                          
         DO 520 JB = 1,JBMAX                                            
            AJ = BG(JB)                                                 
            AA = AI+AJ                                                  
            AAINV = ONE/AA                                              
            DUM = AJ*ARRI*AAINV                                         
            IF (DUM .GT. TOL) GO TO 520                                 
            CSJ = CSB(JB)                                               
            CPJ = CPB(JB)                                               
            CDJ = CDB(JB)                                               
            CFJ = CFB(JB)                                               
            CGJ = CGB(JB)                                               
            CHJ = CHB(JB)                                               
            CIJ = CIB(JB)                                               
            NM = 49*NIJ                                                 
            NN = NM                                                     
            NIJ = NIJ+1                                                 
            R(NIJ) = DUM                                                
            A(NIJ) = AA                                                 
            X1(NIJ) = (AXI+AJ*XJ)*AAINV                                 
            Y1(NIJ) = (AYI+AJ*YJ)*AAINV                                 
            Z1(NIJ) = (AZI+AJ*ZJ)*AAINV                                 
C                                                                       
C           ----- DENSITY FACTOR                                        
C                                                                       
            DUM1 = ZERO                                                 
            DUM2 = ZERO                                                 
            DO 420 I = MINI,MAXI                                        
               GO TO (200,220,420,420,240,420,420,260,420,420,          
     1                261,420,420,262,420,420,420,420,420,263,          
     1                264,420,420,265,420,420,420,420,420,266,          
     1                420,420,267,420,420,                              
     1                268,420,420,269,420,420,420,420,420,270,          
     1                420,420,420,420,420,271,420,420,272,420,          
     1                420,                                              
     1                273,420,420,274,420,420,420,420,420,275,          
     1                420,420,420,420,420,276,420,420,277,420,          
     1                420,278,420,420,420,420,420,279),I                
  200          DUM1 = CSI*AAINV                                         
               GO TO 280                                                
  220          DUM1 = CPI*AAINV                                         
               GO TO 280                                                
  240          DUM1 = CDI*AAINV                                         
               GO TO 280                                                
  260          IF (NORM) DUM1 = DUM1*SQRT3                              
               GO TO 280                                                
  261          DUM1 = CFI*AAINV                                         
               GO TO 280                                                
  262          IF (NORM) DUM1 = DUM1*SQRT5                              
               GO TO 280                                                
  263          IF (NORM) DUM1 = DUM1*SQRT3                              
               GO TO 280                                                
  264          DUM1 = CGI*AAINV                                         
               GO TO 280                                                
  265          IF (NORM) DUM1 = DUM1*SQRT7                              
               GO TO 280                                                
  266          IF (NORM) DUM1 = DUM1*SQRT5/SQRT3                        
               GO TO 280                                                
  267          IF (NORM) DUM1 = DUM1*SQRT3                              
               GO TO 280                                                
  268          DUM1 = CHI*AAINV                                         
               GO TO 280                                                
  269          IF (NORM) DUM1 = DUM1*SQRT9                              
               GO TO 280                                                
  270          IF (NORM) DUM1 = DUM1*SQRT7/SQRT3                        
               GO TO 280                                                
  271          IF (NORM) DUM1 = DUM1*SQRT3                              
               GO TO 280                                                
  272          IF (NORM) DUM1 = DUM1*SQRT5/SQRT3                        
               GO TO 280                                                
  273          DUM1 = CII*AAINV                                         
               GO TO 280                                                
  274          IF (NORM) DUM1 = DUM1*SQRT11                             
               GO TO 280                                                
  275          IF (NORM) DUM1 = DUM1*SQRT3                              
               GO TO 280                                                
  276          IF (NORM) DUM1 = DUM1*SQRT3                              
               GO TO 280                                                
  277          IF (NORM) DUM1 = DUM1*SQRT7/(SQRT5*SQRT3)                
               GO TO 280                                                
  278          IF (NORM) DUM1 = DUM1*SQRT5                              
               GO TO 280                                                
  279          IF (NORM) DUM1 = DUM1*SQRT5/SQRT3                        
C                                                                       
  280          IF (IANDJ) MAX = I                                       
               DO 400 J = MINJ,MAX                                      
                  GO TO (300,320,400,400,340,400,400,360,400,400,       
     1                   361,400,400,362,400,400,400,400,400,363,       
     1                   364,400,400,365,400,400,400,400,400,366,       
     1                   400,400,367,400,400,                           
     1                   368,400,400,369,400,400,400,400,400,370,       
     1                   400,400,400,400,400,371,400,400,372,400,       
     1                   400,                                           
     1                   373,400,400,374,400,400,400,400,400,375,       
     1                   400,400,400,400,400,376,400,400,377,400,       
     1                   400,378,400,400,400,400,400,379),J             
  300             DUM2 = DUM1*CSJ                                       
                  GO TO 380                                             
  320             DUM2 = DUM1*CPJ                                       
                  GO TO 380                                             
  340             DUM2 = DUM1*CDJ                                       
                  GO TO 380                                             
  360             IF (NORM) DUM2 = DUM2*SQRT3                           
                  GO TO 380                                             
  361             DUM2 = DUM1*CFJ                                       
                  GO TO 380                                             
  362             IF (NORM) DUM2 = DUM2*SQRT5                           
                  GO TO 380                                             
  363             IF (NORM) DUM2 = DUM2*SQRT3                           
                  GO TO 380                                             
  364             DUM2 = DUM1*CGJ                                       
                  GO TO 380                                             
  365             IF (NORM) DUM2 = DUM2*SQRT7                           
                  GO TO 380                                             
  366             IF (NORM) DUM2 = DUM2*SQRT5/SQRT3                     
                  GO TO 380                                             
  367             IF (NORM) DUM2 = DUM2*SQRT3                           
                  GO TO 380                                             
  368             DUM2 = DUM1*CHJ                                       
                  GO TO 380                                             
  369             IF (NORM) DUM2 = DUM2*SQRT9                           
                  GO TO 380                                             
  370             IF (NORM) DUM2 = DUM2*SQRT7/SQRT3                     
                  GO TO 380                                             
  371             IF (NORM) DUM2 = DUM2*SQRT3                           
                  GO TO 380                                             
  372             IF (NORM) DUM2 = DUM2*SQRT5/SQRT3                     
                  GO TO 380                                             
  373             DUM2 = DUM1*CIJ                                       
                  GO TO 380                                             
  374             IF (NORM) DUM2 = DUM2*SQRT11                          
                  GO TO 380                                             
  375             IF (NORM) DUM2 = DUM2*SQRT3                           
                  GO TO 380                                             
  376             IF (NORM) DUM2 = DUM2*SQRT3                           
                  GO TO 380                                             
  377             IF (NORM) DUM2 = DUM2*SQRT7/(SQRT5*SQRT3)             
                  GO TO 380                                             
  378             IF (NORM) DUM2 = DUM2*SQRT5                           
                  GO TO 380                                             
  379             IF (NORM) DUM2 = DUM2*SQRT5/SQRT3                     
C                                                                       
  380             NN = NN+1                                             
                  DDIJ(NN) = DUM2                                       
  400          CONTINUE                                                 
  420       CONTINUE                                                    
            IF ( .NOT. IANDJ) GO TO 520                                 
            IF (IA .EQ. JB) GO TO 520                                   
            GO TO (500,440,460,455,450,445,444),LIT                     
  440       IF (MINI .EQ. 2) GO TO 500                                  
            DDIJ(NM+2) = DDIJ(NM+2)+CSI*CPJ*AAINV                       
            GO TO 480                                                   
  444       DDIJ(NM+28) = DDIJ(NM+28)+DDIJ(NM+28)                       
            DDIJ(NM+27) = DDIJ(NM+27)+DDIJ(NM+27)                       
            DDIJ(NM+26) = DDIJ(NM+26)+DDIJ(NM+26)                       
            DDIJ(NM+25) = DDIJ(NM+25)+DDIJ(NM+25)                       
            DDIJ(NM+24) = DDIJ(NM+24)+DDIJ(NM+24)                       
            DDIJ(NM+23) = DDIJ(NM+23)+DDIJ(NM+23)                       
            DDIJ(NM+22) = DDIJ(NM+22)+DDIJ(NM+22)                       
            DDIJ(NM+21) = DDIJ(NM+21)+DDIJ(NM+21)                       
            DDIJ(NM+20) = DDIJ(NM+20)+DDIJ(NM+20)                       
            DDIJ(NM+19) = DDIJ(NM+19)+DDIJ(NM+19)                       
            DDIJ(NM+18) = DDIJ(NM+18)+DDIJ(NM+18)                       
            DDIJ(NM+17) = DDIJ(NM+17)+DDIJ(NM+17)                       
            DDIJ(NM+16) = DDIJ(NM+16)+DDIJ(NM+16)                       
  445       DDIJ(NM+15) = DDIJ(NM+15)+DDIJ(NM+15)                       
            DDIJ(NM+14) = DDIJ(NM+14)+DDIJ(NM+14)                       
            DDIJ(NM+13) = DDIJ(NM+13)+DDIJ(NM+13)                       
            DDIJ(NM+12) = DDIJ(NM+12)+DDIJ(NM+12)                       
            DDIJ(NM+11) = DDIJ(NM+11)+DDIJ(NM+11)                       
  450       DDIJ(NM+10) = DDIJ(NM+10)+DDIJ(NM+10)                       
            DDIJ(NM+9) = DDIJ(NM+9)+DDIJ(NM+9)                          
            DDIJ(NM+8) = DDIJ(NM+8)+DDIJ(NM+8)                          
            DDIJ(NM+7) = DDIJ(NM+7)+DDIJ(NM+7)                          
  455       DDIJ(NM+6) = DDIJ(NM+6)+DDIJ(NM+6)                          
            DDIJ(NM+5) = DDIJ(NM+5)+DDIJ(NM+5)                          
            DDIJ(NM+4) = DDIJ(NM+4)+DDIJ(NM+4)                          
  460       DDIJ(NM+2) = DDIJ(NM+2)+DDIJ(NM+2)                          
  480       DDIJ(NM+3) = DDIJ(NM+3)+DDIJ(NM+3)                          
  500       DDIJ(NM+1) = DDIJ(NM+1)+DDIJ(NM+1)                          
  520    CONTINUE                                                       
  540 CONTINUE                                                          
      RETURN                                                            
      END                                                               
C*MODULE INT2A   *DECK INTIN                                            
C>    @brief    read in $INTGRL group                                   
C>    @author   ??                                                      
C>    @date     Sep, 2019 Peng Xu and Tosaporn Sattasathuchana          
C>             - adding comments for QM-EFP2 MPI/OMP runs               
C>                                                                      
      SUBROUTINE INTIN                                                  
      USE camdft, ONLY: CAMFLAG                                         
      USE lrcdft, ONLY: LCFLAG, EMU, EMU2, LRFILE                       
      use mx_limits, only: mxatm,mxgrid                                 
C                                                                       
      USE params, ONLY: intomp, shfock                                  
C$    USE omp_lib                                                       
C                                                                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C                                                                       
      DOUBLE PRECISION INTGRL,MCSCF,MOROKM                              
      LOGICAL PK,PANDK,BLOCK,PACK2E,SCHWRZ,DIRSCF,FDIFF,                
     *        GOPARR,DSKWRK,MASWRK,SVDSKW,DIRTRF,SCREEN,                
     *        MOIDON,EDCOMP,DIPDCM,DEPRNT,QADDCM,ZDO,                   
     *        POLDCM,POLANG,POLAPP,KMIDPT,POLDYN,                       
     *        QOPS,QFMM,OK,SG1                                          
C                                                                       
      PARAMETER (NNAM=22)                                               
      PARAMETER (NMO=500)                                               
C                                                                       
      DIMENSION QNAM(NNAM),KQNAM(NNAM)                                  
C                                                                       
      COMMON /DFGRID/ DFTTHR,DFTGTHR,SWOFF,SW0,BSLRD(137),NDFTFG,       
     *                NRAD,NTHE,NPHI,NRAD0,NTHE0,NPHI0,                 
     *                NANGPT(MXGRID),NANGPT0(MXGRID),SG1,JANS           
      COMMON /DFTPAR/ DFTTYP(20),EXENA,EXENB,EXENC,IDFT34,NAUXFUN,      
     *                                                    NAUXSHL       
C$omp threadprivate(/DFTPAR/)
      COMMON /EDCMP / ZIJ(NMO),ZMO(5,NMO),OCCUP(NMO),DPFREQ(50),        
     *                MOIDNO(5,NMO),IJMO(2,NMO),MOIJ(NMO),NMOIJ(NMO),   
     *                NMOAT(NMO),NDPFREQ,IPROT(5),NPROT,                
     *                MOIDON,EDCOMP,DIPDCM,DEPRNT,QADDCM,ZDO,POLDCM,    
     *                POLANG,POLAPP,KMIDPT,POLDYN                       
      COMMON /ELGPMT/ NELONG,NATM,NASPIN,NCT,NBNDAB,NTMLB,IPRI,LDOS     
      COMMON /FMOINF/ NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY              
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,                    
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)                  
      COMMON /INTFIL/ NINTMX,NHEX,NTUPL,PACK2E,INTTYP,IGRDTYP           
      COMMON /INT2IC/ NINTIC,ININTIC,NXXIC,LBUFPIC,LIXIC,LABSIX,NINTIX  
      COMMON /INTOPT/ ISCHWZ,IECP,NECP,IEFLD                            
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)                
      COMMON /ORDOPT/ NORDER,NDAR  ,LDAR  ,NBOXMX,NWORD,NOMEM,NSQUAR    
      COMMON /OPTSCF/ DIRSCF,FDIFF                                      
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK                 
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
      COMMON /PCKLAB/ LABSIZ                                            
      COMMON /PKFIL / PK,PANDK,BLOCK                                    
      COMMON /PRPOPT/ ETOLLZ,ILOCAL,IAHARD                              
      COMMON /QMFM  / SIZE,EPS,DPGD,QFMM,NP,NS,IWS,NPGP,MPMTHD,NUMRD,   
     *                ITERMS,QOPS,ISCUT                                 
      COMMON /RESTAR/ TIMLIM,IREST,NREC,INTLOC,IST,JST,KST,LST          
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL                
      COMMON /SCINP / VLAMB,SCREEN                                      
      COMMON /TRFOPT/ CUTTRF,NWDTRF,MPTRAN,ITRFAO,NOSYMT,IPURTF,DIRTRF  
      COMMON /WFNOPT/ SCFTYP,VBTYP,DFTYPE,TDDFTYP,CITYP,CCTYP,          
     *                MPLEVL,MPCTYP                                     
C                                                                       
      DATA IJKO/24/                                                     
      DATA RHF,UHF,ROHF,GVB,MCSCF                                       
     *    /8HRHF     ,8HUHF     ,8HROHF    ,8HGVB     ,8HMCSCF   /      
      DATA RNONE/8HNONE    /                                            
      DATA MOROKM/8HMOROKUMA/                                           
C                                                                       
      DATA INTGRL/8HINTGRL  /                                           
      DATA QNAM/8HNORDER  ,8HNDAR    ,8HLDAR    ,8HNBOXMX  ,            
     *          8HNWORD   ,8HNINTMX  ,8HNOMEM   ,8HNSQUAR  ,            
     *          8HNOPK    ,                                             
     *          8HIST     ,8HJST     ,8HKST     ,8HLST     ,            
     *          8HNREC    ,8HINTLOC  ,8HSCHWRZ  ,                       
     *          8HSCREEN  ,8HVLAMB   ,8HQFMM    ,8HNINTIC  ,            
     *          8HINTOMP  ,8HSHFOCK  /                                  
      DATA KQNAM/1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,3,0,1   ,1,0/        
C                                                                       
C     ----- INITIALIZE VARIABLES FOR NAMELIST $INTGRL -----             
C                                                                       
      NORDER = 0                                                        
      NDAR   = 2000                                                     
      CALL RASIZE(LDAR)                                                 
      NBOXMX = 200                                                      
      NWORD  = 0                                                        
      NINTMX = 0                                                        
      NOMEM  = 0                                                        
      NSQUAR = 0                                                        
      QFMM  =.FALSE.                                                    
      NINTIC = 0                                                        
      ININTIC= 0                                                        
C                    INITIAL VALUE OF NOPK IS ILLEGAL                   
      NOPK   = -27                                                      
      IST    = 1                                                        
      JST    = 1                                                        
      KST    = 1                                                        
      LST    = 1                                                        
      NREC   = 1                                                        
      INTLOC = 1                                                        
      SCHWRZ =.FALSE.                                                   
      SCREEN=.FALSE.                                                    
      VLAMB=0.0D+00                                                     
      IF(NAT.GT.5.OR.NFG.NE.0) SCHWRZ=.TRUE.                            
C     SETTING SCHWRZ TO .FALSE. FOR FMO FRANKENSTEINIZES A BIZARRE BUG  
C     OF THE FOLLOWING MOST MACABRE NATURE:                             
C     SCHWRZ IS SET HERE FOR THE FAKE $DATA (THAT USUALLY HAS 2-4 ATOMS,
C     WHEREAS REAL FRAGMENTS HAVE OFTEN MORE).                          
C                                                                       
      IF(DIRSCF) SCHWRZ=.TRUE.                                          
                                                                        
      intomp = 0                                                        
C$    intomp = 1                                                        
      shfock = .false.                                                  
C                                                                       
C    ----- READ NAMELIST $INTGRL -----                                  
C                                                                       
      CALL NAMEIO(IR,JRET,INTGRL,NNAM,QNAM,KQNAM,                       
     *            NORDER,NDAR,LDAR,NBOXMX,NWORD,NINTMX,NOMEM,NSQUAR,    
     *            NOPK,IST,JST,KST,LST,NREC,INTLOC,SCHWRZ,SCREEN,VLAMB, 
     *            QFMM,NINTIC,                                          
     *            intomp, shfock,                                       
     *            0,0,                                                  
     *    0,0,0,0,0,  0,0,0,0,0,    0,0,0,0,0,  0,0,0,0,0,              
     *    0,0,0,0,0,  0,0,0,0,0,    0,0,0,0,0,  0,0,0,0,0,              
     *   0,0,0,0,0, 0)                                                  
      IF(JRET.EQ.2) THEN                                                
         IF (MASWRK) WRITE(IW,9000)                                     
         CALL ABRT                                                      
      END IF                                                            
C                                                                       
      IF (QFMM) CALL QFMMIN                                             
C                                                                       
C     LONGRANGE !!                                                      
C     IF(LCFLAG)SCHWRZ =.FALSE.                                         
C                                                                       
C     IF SCREEN THEN SCHWARZ INEQUALITIES ARE SKIPPED                   
C                                                                       
      IF (SCREEN) SCHWRZ=.FALSE.                                        
C                                                                       
      ISCHWZ = 0                                                        
      IF(SCHWRZ) ISCHWZ = 1                                             
C                                                                       
      IF(QFMM) THEN                                                     
          OK = .TRUE.                                                   
          CALL DERCHK(NDER)                                             
          IF(SCFTYP.EQ.UHF    .AND.  NDER.GT.1)     OK=.FALSE.          
          IF(SCFTYP.EQ.ROHF   .AND.  NDER.GT.1)     OK=.FALSE.          
          IF(SCFTYP.EQ.GVB)                         OK=.FALSE.          
          IF(SCFTYP.EQ.MCSCF)                       OK=.FALSE.          
          IF(DFTYPE.NE.RNONE  .AND.  SCFTYP.NE.RHF) OK=.FALSE.          
          IF(.NOT.OK) THEN                                              
             IF(MASWRK) WRITE(IW,9030)                                  
             CALL ABRT                                                  
          END IF                                                        
      END IF                                                            
C                                                                       
      IF((LCFLAG.or.camflag).AND.NINTIC.NE.0) THEN                      
         IF(MASWRK) WRITE(IW,9040)                                      
         CALL ABRT                                                      
      END IF                                                            
C                                                                       
C$    IF (intomp.eq.0) shfock = .FALSE.                                 
C                                                                       
C$    IF (MASWRK)                                                       
C$   *      WRITE(IW,'(/10X,14(1H-)/10X,"OPENMP OPTIONS"/10X,14(1H-)/   
C$   *          1X,"INTOMP  =",I8,/                                     
C$   *          1X,"SHFOCK  =",L8)')                                    
C$   *              intomp, shfock                                      
C$    IF (intomp.GT.2.OR.intomp.LT.0) THEN                              
C$      IF(MASWRK) WRITE(IW,95)                                         
C$      CALL ABRT                                                       
C$    END IF                                                            
C                                                                       
C$    IF (.NOT.DIRSCF) THEN                                             
C$      IF(MASWRK) WRITE(IW,96)                                         
C$      intomp=0                                                        
C$      shfock=.false.                                                  
C$    END IF                                                            
C                                                                       
C$    IF (SCFTYP.EQ.GVB) THEN                                           
C$      IF(MASWRK) WRITE(IW,97)                                         
C$      intomp=0                                                        
C$      shfock=.false.                                                  
C$    END IF                                                            
C                                                                       
C$    IF ((omp_get_max_threads().lt.2).AND.                             
C$   *    (intomp.ge.0).AND.MASWRK) THEN                                
C$      WRITE(IW,'(/," WARNING!",                                       
C$   *              " YOU HAVE SELECTED AN OPENMP CODE,",               
C$   *              " BUT ONLY ONE THREAD IS AVAILABLE" )')             
C$      WRITE(IW,*)                                                     
C$   *   "MAYBE YOU FOGOT TO SET OMP_NUM_THREADS ENVIRONMENT VARIABLE?" 
CC$      intomp=0                                                       
CC$      shfock=.false.                                                 
C$    END IF                                                            
C                                                                       
C$    IF ((intomp.NE.0).AND.MASWRK) THEN                                
C$      WRITE(IW,98) omp_get_max_threads()                              
C$    END IF                                                            
C                                                                       
C$    IF (shfock.AND.(intomp.NE.2)) THEN                                
C$      intomp = 2                                                      
C$      IF (MASWRK) THEN                                                
C$      WRITE(IW,'(/," WARNING!",                                       
C$   *               " USING SHARED FOCK MATRIX CODE ASSUMES INTOMP=2", 
C$   *             /," SETTING INTOMP=2")')                             
C$      END IF                                                          
C$    END IF                                                            
C                                                                       
C$    IF (MASWRK) THEN                                                  
C                                                                       
C$      SELECT CASE (intomp)                                            
C$      CASE (0)                                                        
C$        WRITE(IW,'(/," 2E INTS ALGORITHM: MPI(IJ), NO OPENMP")')      
C$      CASE (1)                                                        
C$        WRITE(IW,'(/," 2E INTS ALGORITHM: MPI(I)  + OPENMP(JK)")')    
C$      CASE (2)                                                        
C$        WRITE(IW,'(/," 2E INTS ALGORITHM: MPI(IJ) + OPENMP(KL)")')    
C$      END SELECT                                                      
C                                                                       
C$      LVLEFP = LEVELEFP() ! To see if EFP2 fragments present          
C$                                                                      
C$      IF (intomp.gt.0) THEN                                           
C$          IF (shfock) THEN                                            
C$              WRITE(IW,*) "FOCK MATRIX WILL BE SHARED AMONG THREADS"  
C$                                                                      
C$              IF(LVLEFP.EQ.2) THEN                                    
C$                WRITE(IW,*)                                           
C$                WRITE(IW,*) "CURRENTLY, QM-EFP2 2E INTS (MODES 1,2,3)"
C$                WRITE(IW,*) "DO NOT SHARE FOCK AMONG THREADS"         
C$                WRITE(IW,*) "SHFOCK OPTION ONLY APPLIES TO QM REGION" 
C$              END IF                                                  
C$          ELSE                                                        
C$              WRITE(IW,*) "PRIVATE FOCK MATRIX CODE WILL BE USED"     
C$          END IF                                                      
C$      END IF                                                          
C                                                                       
C$    END IF                                                            
C                                                                       
C     ----- DETERMINE WHICH 2E- INTEGRAL SCHEME TO USE -----            
C     NOPK=0 MEANS A -P- SUPERMATRIX, AND POSSIBLY -K- SUPERMATRIX.     
C     NOPK=1 MEANS A CONVENTIONAL -J- INTEGRAL LIST WILL BE USED.       
C     SUPERMATRIX IS INCORRECT FOR MOROKUMA ANALYSIS.                   
C     ANYTHING THAT DOES AN INTEGRAL TRANSFORMATION (E.G. MCSCF,        
C     CI, MP2, RUEDENBERG LOCALIZATION, ANALYTIC HESSIAN, LOCALIZED     
C     ORBITAL POLARIZABILITIES) MAY NOT USE A SUPERMATRIX.              
C     IN CORE RUNS ARE NOT YET TAUGHT TO USE NOPK (THEY COULD BE).      
C                                                                       
C     IN CASE OF FMO PK OPTION IS DOOMED IN MANY PLACES. THE HARD       
C     ONES TO FIND ARE ZPKOUT AND PKFILE, WHERE THE FIRST ROUTINE       
C     IS SET FOR THE FIRST FRAGMENT ONLY (SEE SHELLS FOR A SOLUTION).   
C     OTHER PLACES ARE PROBABLY WRONG TOO (NEVER CHECKED).              
C                                                                       
      CALL DERCHK(NDER)                                                 
      IF(NDFTFG.NE.0)              NOPK = 1                             
      IF(DFTTYP(1) .NE. 0.0D+00)   NOPK = 1                             
      IF(RUNTYP.EQ.MOROKM)         NOPK = 1                             
      IF(SCFTYP.EQ.MCSCF)          NOPK = 1                             
      IF(CITYP.NE.RNONE)           NOPK = 1                             
      IF(ILOCAL.EQ.2)              NOPK = 1                             
      IF(POLDCM)                   NOPK = 1                             
      IF(POLDYN)                   NOPK = 1                             
      IF(NDER.EQ.2)                NOPK = 1                             
      IF(NINTIC.NE.0)              NOPK = 1                             
      IF(NFG.NE.0)                 NOPK = 1                             
C                                                                       
C         SINCE J FILE IS MUCH SHORTER THAN P OR PK FILE,               
C         AND I/O IS STILL QUITE BAD EVEN IF STRIPING,                  
C         WE MAKE ALL INTEGRAL FILES NON-SUPERMATRIX,                   
C         UNLESS THE USER HAS SPECIFICALLY ASKED FOR ONE.               
C                                                                       
      IF(NOPK.EQ.-27)              NOPK = 1                             
      IF(NOPK.NE.1)                NOPK = 0                             
C                                                                       
      PK = NOPK.EQ.0                                                    
      PANDK = (PK .AND. SCFTYP.EQ. UHF) .OR.                            
     *        (PK .AND. SCFTYP.EQ.ROHF) .OR.                            
     *        (PK .AND. SCFTYP.EQ. GVB)                                 
C                                                                       
      IF(      PANDK .AND. NINTMX.EQ.0) NINTMX=10000                    
      IF(.NOT. PANDK .AND. NINTMX.EQ.0) NINTMX=15000                    
C                                                                       
C     ----- OPEN INTEGRAL FILES -----                                   
C                                                                       
      IF(DIRSCF  .OR.  DIRTRF) THEN                                     
         PK    = .FALSE.                                                
         PANDK = .FALSE.                                                
         NOPK  = 1                                                      
      ELSE                                                              
         SVDSKW = DSKWRK                                                
         DSKWRK = .TRUE.                                                
         CALL SEQOPN(IS,'AOINTS','UNKNOWN',.FALSE.,'UNFORMATTED')       
         IF(LCFLAG.OR.CAMFLAG) THEN                                     
           CALL SEQOPN(LRFILE,'DFTINTS','UNKNOWN',.FALSE.,'UNFORMATTED')
         END IF                                                         
C                                                                       
C     EXCHANGE -IJKO- AND -IS- SO WE DON'T GET TWO COPIES OF 2E-INTS    
C                                                                       
         IF (RUNTYP.EQ.MOROKM) THEN                                     
            CALL SEQOPN(IJKO,'ORDINT','UNKNOWN',.FALSE.,'UNFORMATTED')  
            ISK = IJKO                                                  
            IJKO = IS                                                   
            IS = ISK                                                    
         END IF                                                         
         DSKWRK = SVDSKW                                                
      END IF                                                            
C                                                                       
C     ARRAY DIMENSIONS THROUGH OUT THE PROGRAM ALLOW -MXAO- AO-S,       
C     EXCEPT IN THE CI CODE WHERE THE LIMIT IS -MXAOCI-.                
C     THESE LIMITS ARE TESTED NEAR WHERE THE MOLECULE IS READ IN,       
C     AND NEED NOT CONCERN US FURTHER.                                  
C     FOR CONVENTIONAL SCF, WITH A P OR PK SUPERMATRIX, THE 16          
C     BIT PACKING LIMITS THE AOS TO 361.  ELSEWHERE WE USE 8            
C     BITS TO HANDLE CASES SMALLER THAN 256, 16 BITS OTHERWISE.         
C                                                                       
      MAXAO=255                                                         
      IF(PK) MAXAO=361                                                  
      LABSIZ = 1                                                        
      IF(NUM.GT.MAXAO.OR.NFG.NE.0) LABSIZ = 2                           
C                                                                       
C     DURING IN-CORE FMO RUNS LABSIZ SHOULD NOT BE SET BASED            
C     ON THE LIBRARY $DATA (THAT DELIBERATELY RESULTS IN LABSIZ=1).     
C     THE SOLUTION IS TO SET LABSIZ TO 2 ALWAYS AND RECOMPUTE FOR       
C     EACH FRAGMENT LATER. POSITIVE INPUT VALUES OF NINTIC MAY          
C     RESULT IN NOT FULLY USING AVAILABLE MEMORY AND SHOULD BE          
C     AVOIDED (USE NEGATIVE VALUES IN $INTGRL NINTIC INSTEAD).          
C                                                                       
      IF(PK  .AND.  NUM.GT.MAXAO) THEN                                  
         IF(MASWRK) WRITE(IW,*) 'TOO MANY AOS TO USE PK SUPERMATRIX.'   
         IF(MASWRK) WRITE(IW,*) 'USE J INTEGRAL LIST, OR DIRECT SCF'    
         CALL ABRT                                                      
      END IF                                                            
C                                                                       
C        ELONGATION METHOD                                              
C                                                                       
      IF(NELONG.GT.0) LABSIZ = 2                                        
C                                                                       
      IF (MASWRK) WRITE(IW,9010) NOPK,NORDER,SCHWRZ                     
C                                                                       
C     CHECK IF IT IS AN SCREEN CALCULATION. IF SO, TELL IT ,            
C     MAKE SCHWZ FALSE AND                                              
C     CALL TO THE INITFCTS SUBROUTINE TO GENERATE THE VALUES            
C     FOR COMMON /FCTS/                                                 
C                                                                       
      IF (SCREEN) THEN                                                  
          WRITE(IW,9011) SCREEN,VLAMB                                   
          CALL INITFCTS                                                 
      END IF                                                            
C                                                                       
      RETURN                                                            
C                                                                       
 9000 FORMAT(1X,'ERROR IN $INTGRL INPUT')                               
 9010 FORMAT(/10X,22(1H-)/10X,'INTEGRAL INPUT OPTIONS'/10X,22(1H-)/     
     *       1X,'NOPK  =',I8,' NORDER=',I8,' SCHWRZ=',L8)               
 9011 FORMAT(/,'SCREENED INTEGRALS',10X,                                
     *       ' SCREEN=',L8,' VLAMB= ',F6.3,/)                           
 9030 FORMAT(1X,'*** ERROR ***'/                                        
     *       1X,'THE QFMM=.TRUE. OPTION IS AVAILABLE ONLY FOR'/         
     *       1X,'CLOSED SHELL RHF OR DFT ENERGIES AND GRADIENTS,'/      
     *       1X,'OR FOR UHF AND ROHF NON-DFT ENERGIES.')                
 9040 FORMAT(/1X,'Set NINTIC to 0 for DFT/LC or CAMB3LYP.')             
C$ 95 FORMAT(/1X,'ERROR, INTOMP MUST BE 0, 1 OR 2'/)                    
C$ 96 FORMAT(/1X,'*** WARNING ***'/,                                    
C$   *        1X,'THREADED 2-ELECTRON INTEGRAL CODE ONLY WORKS',        
C$   *        1X,'WITH DIRECT SCF'/,                                    
C$   *        1X,'INTOMP OPTION WAS IGNORED')                           
C$ 97 FORMAT(/1X,'*** WARNING ***'/,                                    
C$   *        1X,'THREADED 2-ELECTRON INTEGRAL CODE DOES NOT',          
C$   *        1X,'SUPPORT GVB METHOD'/,                                 
C$   *        1X,'INTOMP OPTION WAS IGNORED')                           
C$ 98 FORMAT(/1X,'THREADED 2-ELECTRON INTEGRAL CODE WILL BE USED',      
C$   *       /,' NUMBER OF OPENMP THREADS:', I6)                        
      END                                                               
C*MODULE INT2A   *DECK INTOUT                                           
      SUBROUTINE INTOUT(I1,I2,I3,I4,Q4,NN,VAL)                          
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      LOGICAL GOPARR,DSKWRK,MASWRK                                      
C                                                                       
      COMMON /INTPR / Q(2),V(2),JC,N1(2),J1(2),J2(2),J3(2),J4(2)        
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)                
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
C                                                                       
      JC = JC+1                                                         
      J1(JC) = I1                                                       
      J2(JC) = I2                                                       
      J3(JC) = I3                                                       
      J4(JC) = I4                                                       
      Q(JC) = Q4                                                        
      N1(JC) = NN                                                       
      V(JC) = VAL                                                       
      IF (JC.LT.2) GO TO 100                                            
      JC = 0                                                            
      IF (MASWRK) WRITE (IW,9088)                                       
     *   (J1(M),J2(M),J3(M),J4(M),Q(M),V(M),M = 1,2)                    
  100 CONTINUE                                                          
      RETURN                                                            
C                                                                       
 9088 FORMAT(2(4I4,F5.1,F17.9,1X))                                      
      END                                                               
C*MODULE INT2A   *DECK JANDK                                            
!> @brief Compute 2-electron integrals                                  
!>                                                                      
!> @author Unknown                                                      
!> @date September 2010 - Albert DeFusco                                
!> - make sure exchange integrals are computed for correlated subsystems
      SUBROUTINE JANDK                                                  
      use mx_limits, only: mxsh,mxgtot,mxgsh,mxg2                       
C                                                                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C                                                                       
      CHARACTER*8 INAOFL                                                
C                                                                       
      LOGICAL SCHWRZ,PACK2E,PK,PANDK,BLOCK,DIRTRF,FDIFF                 
      LOGICAL GOPARR,DSKWRK,MASWRK,TDSKWRK,DIRSCF                       
      LOGICAL LCUT                                                      
C                                                                       
C                                                                       
      COMMON /ELGIDX/ LCUT                                              
      COMMON /ELGPMT/ NELONG,NATM,NASPIN,NCT,NBNDAB,NTMLB,IPRI,LDOS     
      COMMON /ELGFIL/ NFILE,INAOFL                                      
      COMMON /FMCOM / XX(1)                                             
      COMMON /FMOINF/ NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY              
      COMMON /FMOOPT/ ESPSCA(9),RESPAP(2),RESPPC(2),RESDIM,RESTRI(4),   
     *                RCORSD,RESPCT,CONVFG,CNVDMP,COROFF,RFLMO(4),      
     *                ORSHFT,ORSHFT2,CNVAFO,ASCREEN(4),IXESP,MXITFG,    
     *                NGUESS,NBSSE,MODORB,MODPAR,IRSTSTP,IRSTLAY,NPRFMO,
     *                NFMOPAL,MODPRP,MAXL1C,IPIEDA,MODGRD,MODESP,IVMUL, 
     *                MODLMO,NOPDEN,MOFOCK,MODFD,modfmm,ncentm,ndualb   
      COMMON /INT2IC/ NINTIC,ININTIC,NXXIC,LBUFPIC,LIXIC,LABSIX,NINTIX  
      COMMON /INTFIL/ NINTMX,NHEX,NTUPL,PACK2E,INTTYP,IGRDTYP           
      COMMON /INTOPT/ ISCHWZ,IECP,NECP,IEFLD                            
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)                
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),      
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),      
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),   
     *                KLOC(MXSH),MIN(MXSH),MAX(MXSH),NSHELL             
      COMMON /OPTSCF/ DIRSCF,FDIFF                                      
      COMMON /ORDOPT/ NORDER,NDAR,LDAR,NBOXMX,NWORD,NOMEM,NSQUAR        
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK                 
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
      COMMON /SIMDAT/ NACC,NREJ,IGOMIN,NRPA,IBWM,NACCT,NREJT,NRPAT,     
     *                NPRTGO,IDPUNC,IGOFLG                              
      COMMON /PKFIL / PK,PANDK,BLOCK                                    
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL                
      COMMON /TRFOPT/ CUTTRF,NWDTRF,MPTRAN,ITRFAO,NOSYMT,IPURTF,DIRTRF  
      COMMON /WFNOPT/ SCFTYP,VBTYP,DFTYPE,TDDFTYP,CITYP,CCTYP,          
     *                MPLEVL,MPCTYP                                     
C                                                                       
      DATA CHECK/8HCHECK   /                                            
      DATA NONE/4HNONE/                                                 
      DATA CIMSUB/8HCIMSUB  /                                           
C                                                                       
C     ----- MAIN DRIVER FOR CALCULATION OF 2E- INTEGRALS -----          
C                                                                       
C     ----- MOPAC 2-ELECTRON INTEGRALS ALREADY CALCULATED -----         
C                                                                       
      IF(MPCTYP.NE.NONE) RETURN                                         
C                                                                       
C            ELONGATION METHOD                                          
      IF(NELONG.GT.1) THEN                                              
         INAOFL = 'AOINTS  '                                            
         NFILE = 1                                                      
      END IF                                                            
C                                                                       
C     HONDO INTEGRAL PACKAGE REQUIRES A BUFFER DIMENSIONED FOR          
C     THE MAXIMUM ANGULAR MOMENTUM OCCURING IN THE BASIS SET.           
C     A PURE SP BASIS USING THE POPLE INTEGRALS DOESN'T NEED THIS.      
C                                                                       
      CALL BASCHK(LMAX)                                                 
                    NANGM =  4                                          
      IF(LMAX.EQ.2) NANGM =  6                                          
      IF(LMAX.EQ.3) NANGM = 10                                          
      IF(LMAX.EQ.4) NANGM = 15                                          
      IF(LMAX.EQ.5) NANGM = 21                                          
      IF(LMAX.EQ.6) NANGM = 28                                          
                 MAXG = NANGM**4                                        
      IF(PK)     MAXG = MAXG*3                                          
      IF(DIRSCF) MAXG = 1                                               
C                                                                       
      NSH2 = (NSHELL*NSHELL+NSHELL)/2                                   
      MINTMX=NINTMX                                                     
      IF(NINTIC.NE.0) MINTMX=0                                          
C     THE MEMORY IS ALLOCATED ELSEWHERE                                 
C                                                                       
      CALL VALFM(LOADFM)                                                
      LBUFP  = LOADFM + 1                                               
      LBUFK  = LBUFP  + MINTMX                                          
      LIX    = LBUFK  + MINTMX                                          
      LXINTS = LIX    + MINTMX                                          
      LGHOND = LXINTS + NSH2                                            
      LDDIJ  = LGHOND + MAXG                                            
      LAST   = LDDIJ  + 49*MXG2                                         
      NEED = LAST - LOADFM - 1                                          
      CALL GETFM(NEED)                                                  
      IF(NINTIC.NE.0) THEN                                              
         LBUFP=LBUFPIC                                                  
         LIX=LIXIC                                                      
      ENDIF                                                             
C                                                                       
      IF(EXETYP.EQ.CHECK) THEN                                          
         IF(ISCHWZ.GT.0) THEN                                           
            CALL VCLR(XX(LXINTS),1,NSH2)                                
            CALL DAWRIT(IDAF,IODA,XX(LXINTS),NSH2,54,0)                 
         END IF                                                         
      END IF                                                            
C                                                                       
C           INITIALIZE THE CALCULATION                                  
C                                                                       
      TDSKWRK = DSKWRK                                                  
      DSKWRK  = .TRUE.                                                  
      CALL DEBUT(DIRSCF,XX(LBUFP),XX(LBUFK),XX(LIX),NINTMX,NEED,DIRTRF) 
      IF(EXETYP.EQ.CHECK)  GO TO 200                                    
      IF(runtyp.eq.cimsub.and.dirtrf)                                   
     *  CALL EXCHNG(XX(LXINTS),XX(LGHOND),XX(LDDIJ),                    
     *                       NSH2,MAXG,INTTYP)                          
      IF(DIRSCF.OR.DIRTRF) GO TO 200                                    
C                                                                       
C     ----- PACKING PARAMETERS: THIS IS AN INACTIVE OPTION.             
C           NHEX = DESIRED HEXADECIMAL ACCURACY.                        
C           NTUPL= # OF BYTES PER INTEGER WORD TO BE PACKED.            
C                  THE NEGATIVE VALUE OF -NTUPL- MEANS THAT THE         
C                  INTEGER LABELS ARE NOT TO BE PACKED.                 
C                                                                       
      PACK2E = ICUT .LT. 0                                              
      NTUPL = -4                                                        
      IF (PACK2E .AND. MASWRK) WRITE (IW,9088) NHEX,NTUPL               
 9088 FORMAT(41H THE INTEGRALS ARE PACKED WITH A.D.MCLEAN,              
     *     26H PACKING UTILITIES (1977).,8H NHEX = ,I5,9H NTUPL = ,I5/) 
C                                                                       
C           GENERATE ALL EXCHANGE INTEGRALS                             
C                                                                       
      SCHWRZ = ISCHWZ.GT.0                                              
      IF(SCHWRZ) CALL EXCHNG(XX(LXINTS),XX(LGHOND),XX(LDDIJ),           
     *                       NSH2,MAXG,INTTYP)                          
C                                                                       
C           CRUNCH OUT THE INTEGRALS                                    
C                                                                       
      IDUMMY=0                                                          
      DUMMY=0.0D+00                                                     
      CALL TWOEI(SCFTYP,.FALSE.,.FALSE.,.FALSE.,.FALSE.,                
     *           INTTYP,SCHWRZ,NINT,NSCHWZ,1,1,                         
     *           XX(LBUFP),XX(LBUFK),XX(LIX),NINTMX,                    
     *           XX(LXINTS),NSH2,XX(LGHOND),MAXG,XX(LDDIJ),             
     *           IDUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,1)    
C                                                                       
C            ELONGATION METHOD INTEGRAL FILES                           
      IF(NELONG.GT.1.AND.NFILE.GT.1) THEN                               
         CALL SEQCLO(IS,'KEEP')                                         
         CALL SEQOPN(IS,INAOFL(1:6),'UNKNOWN',.FALSE.,'UNFORMATTED')    
         CALL SEQREW(IS)                                                
      END IF                                                            
C                                                                       
  200 CONTINUE                                                          
      DSKWRK  = TDSKWRK                                                 
      CALL RETFM(NEED)                                                  
C                                                                       
C            ELONGATION METHOD                                          
      IF(LCUT) THEN                                                     
         CALL DEBTCT(DIRSCF)                                            
         IF(DIRSCF) GOTO 300                                            
C                                                                       
         INAOFL = 'EGINTA  '                                            
         NFILE = 1                                                      
C                                                                       
         NSH2 = (NSHELL*NSHELL+NSHELL)/2                                
         CALL VALFM(LOADFM)                                             
         LBUFP  = LOADFM + 1                                            
         LBUFK  = LBUFP  + NINTMX                                       
         LIX    = LBUFK  + NINTMX                                       
         LXINTS = LIX    + NINTMX                                       
         LGHOND = LXINTS + NSH2                                         
         LDDIJ  = LGHOND + MAXG                                         
         LAST   = LDDIJ  + 49*MXG2                                      
         NEED = LAST - LOADFM - 1                                       
         CALL GETFM(NEED)                                               
C                                                                       
         IDUMMY=0                                                       
         DUMMY=0.0D+00                                                  
         IF(SCHWRZ) CALL EXCHNG(XX(LXINTS),XX(LGHOND),XX(LDDIJ),        
     *                          NSH2,MAXG,INTTYP)                       
C                                                                       
         CALL ADDINT(SCFTYP,.FALSE.,.FALSE.,.FALSE.,.FALSE.,            
     *               SCHWRZ,NINT,NSCHWZ,1,                              
     *               XX(LBUFP),XX(LBUFK),XX(LIX),NINTMX,                
     *               XX(LXINTS),NSH2,XX(LGHOND),MAXG,                   
     *               IDUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,              
     *               .FALSE.,1,1,1)                                     
         CALL RETFM(NEED)                                               
300      CONTINUE                                                       
      ENDIF                                                             
C                                                                       
      IF (MASWRK.AND.NPRTGO.NE.2.AND.(NFG.EQ.0.OR.IAND(NPRFMO,3).EQ.0)) 
     *   WRITE(IW,*) ' ...... END OF TWO-ELECTRON INTEGRALS .....'      
      IF (MASWRK.AND.(NFG.EQ.0.OR.IAND(NPRFMO,3).EQ.0)) CALL TEXIT(1,1) 
C                                                                       
C     ----- INTEGRAL ORDERING -----                                     
C                                                                       
      IF (NORDER.EQ.1) CALL ORDRJK                                      
      RETURN                                                            
      END                                                               
C*MODULE INT2A   *DECK PKFILE                                           
      SUBROUTINE PKFILE(II,JJ,KK,LL,SKIPA,SKIPB,SKIPC,NPSYM,            
     *                  BUFP,BUFK,IX,NINTMX,GHONDO)                     
      USE lrcdft, ONLY: LRFILE                                          
      use mx_limits, only: mxsh,mxgtot,mxao                             
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      DIMENSION BUFP(NINTMX),BUFK(NINTMX),IX(NINTMX),GHONDO(*)          
      DIMENSION IB(15),JB(15),KB(15),LB(15)                             
C                                                                       
      LOGICAL SKIPA,SKIPB,SKIPC,NPSYM                                   
      LOGICAL PKFL,PANDK,BLOCK,IANDJ,KANDL,SAME,OUT,FIRST               
      LOGICAL   LRINT                                                   
C                                                                       
C                                                                       
      COMMON /IJPAIR/ IJADD(MXAO)                                       
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)                
      COMMON /MISC  / IANDJ,KANDL,SAME                                  
C$omp threadprivate(/MISC/)
      COMMON /NLRCF / LRINT                                             
C$omp threadprivate(/NLRCF /)
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),      
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),      
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),   
     *                KLOC(MXSH),MIN(MXSH),MAX(MXSH),NSHELL             
      COMMON /PCKLAB/ LABSIZ                                            
      COMMON /PKFIL / PKFL,PANDK,BLOCK                                  
      COMMON /RESTAR/ TIMLIM,IREST,NREC,INTLOC,IST,JST,KST,LST          
      COMMON /SHLEXC/ NORGSH(3),NORGSP(3),IEXCH,NANGM,NGTH(4)           
      COMMON /SHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,          
     *                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,          
     *                NIJ,IJ,KL,IJKL                                    
C$omp threadprivate(/SHLNOS/)
      COMMON /SHLT  / TOL,CUTOFF,ICOUNT,OUT                             
C                                                                       
      PARAMETER (ZERO=0.0D+00)                                          
      PARAMETER (HALF=0.5D+00)                                          
C                                                                       
      SAVE FIRST,IB,JB,KB,LB                                            
      DATA FIRST/.TRUE./                                                
      ISSAVE=IS                                                         
      IF(LRINT)IS=LRFILE                                                
C                                                                       
C     ----- WRITE THIS SHELL'S P OR PK INTEGRAL FILE -----              
C                                                                       
      IF(FIRST) THEN                                                    
         FIRST=.FALSE.                                                  
         DO 50 I=1,NANGM                                                
            LB(I) = I-1                                                 
            KB(I) = LB(I) * NANGM                                       
            JB(I) = KB(I) * NANGM                                       
            IB(I) = JB(I) * NANGM                                       
  50     CONTINUE                                                       
      END IF                                                            
C                                                                       
      IND = 1                                                           
      IF (SKIPA .AND. NPSYM) IND = 2                                    
      IF (SKIPB .AND. NPSYM) IND = 3                                    
      IF (SKIPC .AND. NPSYM) IND = 3                                    
      IF (SKIPA .AND. SKIPB .AND. NPSYM) IND = 4                        
C                                                                       
      NORG1 = NORGSH(1)+1                                               
      NORG2 = NORGSH(2)+1                                               
      NORG3 = NORGSH(3)+1                                               
C                                                                       
      LIT = KTYPE(II)                                                   
      LKT = KTYPE(KK)                                                   
      MINI = MIN(II)                                                    
      MINJ = MIN(JJ)                                                    
      MINK = MIN(KK)                                                    
      MINL = MIN(LL)                                                    
      MAXI = MAX(II)                                                    
      MAXJ = MAX(JJ)                                                    
      MAXK = MAX(KK)                                                    
      MAXL = MAX(LL)                                                    
      LOCI = KLOC(II)-MINI                                              
      LOCJ = KLOC(JJ)-MINJ                                              
      LOCK = KLOC(KK)-MINK                                              
      LOCL = KLOC(LL)-MINL                                              
C                                                                       
      IANDJ = II .EQ. JJ                                                
      KANDL = KK .EQ. LL                                                
      SAME = (II .EQ. KK) .AND. (JJ .EQ. LL)                            
C                                                                       
C     TYPE = 1 FOR (II II II II)                                        
C            2     (II JJ JJ JJ)                                        
C            3     (II II KK KK) AND  LIT.GE.LKT                        
C            4     (II II KK KK) AND  LIT.LT.LKT                        
C            5     (II II II LL)                                        
C            6     (II JJ KK KK)                                        
C            7     (II JJ JJ LL)                                        
C            8     (II II KK LL)                                        
C            9     (II JJ KK LL)                                        
C                                                                       
      NTYP = 0                                                          
      IF (II.EQ.JJ .AND. JJ.EQ.KK .AND. KK.EQ.LL) NTYP = 1              
      IF (II.GT.JJ .AND. JJ.EQ.KK .AND. KK.EQ.LL) NTYP = 2              
      IF (II.EQ.JJ .AND. JJ.GT.KK .AND. KK.EQ.LL                        
     *                          .AND. LIT.GE.LKT) NTYP = 3              
      IF (II.EQ.JJ .AND. JJ.GT.KK .AND. KK.EQ.LL                        
     *                          .AND. LIT.LT.LKT) NTYP = 4              
      IF (II.EQ.JJ .AND. JJ.EQ.KK .AND. KK.GT.LL) NTYP = 5              
      IF (II.GT.JJ .AND. JJ.GT.KK .AND. KK.EQ.LL) NTYP = 6              
      IF (II.GT.JJ .AND. JJ.EQ.KK .AND. KK.GT.LL) NTYP = 7              
      IF (II.EQ.JJ .AND. JJ.GT.KK .AND. KK.GT.LL) NTYP = 8              
      IF (II.GT.JJ .AND. JJ.GT.KK .AND. KK.GT.LL) NTYP = 9              
      IF (SKIPA .AND. .NOT. NPSYM) NORG2 = 1                            
      IF (SKIPC .AND. .NOT. NPSYM) NORG3 = NORGSH(2)+1                  
      IF (SKIPB .AND. .NOT. NPSYM) NORG3 = 1                            
C                                                                       
C     ----- BEGIN LOOPS OVER PRIMITIVES IN THIS SHELL -----             
C                                                                       
C     INTEGRAL TYPES N1,G1 FOR (I,J//K,L)                               
C                    N2,G2 FOR (I,K//J,L)                               
C                    N3,G3 FOR (I,L//J,K)                               
C                                                                       
C        N1 = IB(IA)+JB(JA)+KB(KA)+LB(LA)+NORG1                         
C        N2 = IB(IA)+JB(KA)+KB(JA)+LB(LA)+NORG2                         
C        N3 = IB(IA)+JB(LA)+KB(JA)+LB(KA)+NORG3                         
C                                                                       
      JMAX = MAXJ                                                       
      KMAX = MAXK                                                       
      LMAX = MAXL                                                       
      DO 860 I = MINI,MAXI                                              
      IAO = LOCI + I                                                    
      IA = I-MINI+1                                                     
      N1I = NORG1 + IB(IA)                                              
      N2I = NORG2 + IB(IA)                                              
C                                                                       
      IF (IANDJ) JMAX = I                                               
      DO 840 J = MINJ,JMAX                                              
      IF (JJ .EQ. KK) KMAX = J                                          
      JAO = LOCJ + J                                                    
      JA = J-MINJ+1                                                     
      N1IJ = N1I + JB(JA)                                               
      N2IJ = N2I + KB(JA)                                               
C                                                                       
      DO 820 K = MINK,KMAX                                              
      KAO = LOCK + K                                                    
      KA = K-MINK+1                                                     
      N1IJK = N1IJ + KB(KA)                                             
      N2IJK = N2IJ + JB(KA)                                             
C                                                                       
      IF (KANDL) LMAX = K                                               
      DO 800 L = MINL,LMAX                                              
         LAO = LOCL + L                                                 
         LA = L-MINL+1                                                  
         N1 = N1IJK + LB(LA)                                            
         N2 = N2IJK + LB(LA)                                            
C                                                                       
         GO TO (200,220,230,250,270,280,290,300,310),NTYP               
  200    IF (IA .EQ. JA) GO TO 210                                      
         N3 = IB(IA)+JB(LA)+KB(JA)+LB(KA)+NORG3                         
         GO TO 400                                                      
  210    N3 = IB(JA)+JB(KA)+KB(IA)+LB(LA)+NORG3                         
         GO TO 400                                                      
  220    N3 = IB(IA)+JB(LA)+KB(JA)+LB(KA)+NORG3                         
         GO TO 400                                                      
  230    IF (IA .EQ. JA) GO TO 240                                      
         N3 = IB(IA)+JB(LA)+KB(JA)+LB(KA)+NORG3                         
         GO TO 400                                                      
  240    N3 = IB(JA)+JB(KA)+KB(IA)+LB(LA)+NORG3                         
         GO TO 400                                                      
  250    IF (KA .EQ. LA) GO TO 260                                      
         N3 = IB(JA)+JB(KA)+KB(IA)+LB(LA)+NORG3                         
         GO TO 400                                                      
  260    N3 = IB(IA)+JB(LA)+KB(JA)+LB(KA)+NORG3                         
         GO TO 400                                                      
  270    N3 = IB(JA)+JB(KA)+KB(IA)+LB(LA)+NORG3                         
         GO TO 400                                                      
  280    N3 = IB(IA)+JB(LA)+KB(JA)+LB(KA)+NORG3                         
         GO TO 400                                                      
  290    N3 = IB(IA)+JB(LA)+KB(JA)+LB(KA)+NORG3                         
         GO TO 400                                                      
  300    N3 = IB(JA)+JB(KA)+KB(IA)+LB(LA)+NORG3                         
         GO TO 400                                                      
  310    N3 = IB(IA)+JB(LA)+KB(JA)+LB(KA)+NORG3                         
C                                                                       
C     ----- FORM FIRST LINEAR COMBINATION -----                         
C                                                                       
  400 CONTINUE                                                          
      G1 = GHONDO(N1)                                                   
      G2 = GHONDO(N2)                                                   
      G3 = GHONDO(N3)                                                   
C                                                                       
      JUMP = 1                                                          
      I1 = IAO                                                          
      I2 = JAO                                                          
      I3 = KAO                                                          
      I4 = LAO                                                          
      IF (I2 .EQ. I3) JUMP = 2                                          
      IF ((I2 .EQ. I4) .OR. (I1 .EQ. I3)) JUMP = 3                      
C                                                                       
      GO TO (410,420,430,440),IND                                       
  410 VALK = G2+G3                                                      
      VALP = (G1+G1)+(G1+G1)-VALK                                       
      GO TO 460                                                         
  420 VALK = G3                                                         
      VALP = (G1+G1)+(G1+G1)-VALK                                       
      GO TO 460                                                         
  430 VALK = G2                                                         
      VALP = (G1+G1)+(G1+G1)-VALK                                       
      GO TO 460                                                         
  440 VALK = ZERO                                                       
      VALP = (G1+G1)+(G1+G1)                                            
  460 CONTINUE                                                          
      NN = N1                                                           
      GO TO 700                                                         
C                                                                       
C     ----- FORM SECOND LINEAR COMBINATION -----                        
C                                                                       
  500 CONTINUE                                                          
      GO TO (510,520,530,540),IND                                       
  510 VALK = G3+G1                                                      
      VALP = (G2+G2)+(G2+G2)-VALK                                       
      GO TO 560                                                         
  520 VALK = G1+G3                                                      
      VALP = -VALK                                                      
      GO TO 560                                                         
  530 VALK = G1                                                         
      VALP = (G2+G2)+(G2+G2)-VALK                                       
      GO TO 560                                                         
  540 VALK = G1                                                         
      VALP = -VALK                                                      
  560 CONTINUE                                                          
      NN = N2                                                           
      JUMP = 2                                                          
      IF ((I1 .EQ. I2) .OR. (I3 .EQ. I4)) JUMP = 3                      
      I2 = KAO                                                          
      I3 = JAO                                                          
      GO TO 700                                                         
C                                                                       
C     ----- FORM THIRD LINEAR COMBINATION -----                         
C                                                                       
  600 CONTINUE                                                          
      GO TO (610,620,630,640),IND                                       
  610 VALK = G1+G2                                                      
      VALP = (G3+G3)+(G3+G3)-VALK                                       
      GO TO 660                                                         
  620 VALK = G1                                                         
      VALP = (G3+G3)+(G3+G3)-VALK                                       
      GO TO 660                                                         
  630 VALK = G1+G2                                                      
      VALP = -VALK                                                      
      GO TO 660                                                         
  640 VALK = G1                                                         
      VALP = -VALK                                                      
  660 CONTINUE                                                          
      NN = N3                                                           
      I2 = LAO                                                          
      I3 = JAO                                                          
      I4 = KAO                                                          
      JUMP = 3                                                          
C                                                                       
C     ----- STORE INTEGRAL AND INDICES -----                            
C                                                                       
  700 CONTINUE                                                          
      IF (PANDK) GO TO 740                                              
C                                                                       
C     ----- -P- SUPERMATRIX ONLY -----                                  
C                                                                       
      IF (ABS(VALP) .LT. CUTOFF) GO TO 780                              
      IF (OUT) CALL INTOUT(I1,I2,I3,I4,QQ4,NN,VALP)                     
      IF (I1.EQ.I3 .AND. I2.EQ.I4) VALP = VALP*HALF                     
      I1I2 = IJADD(I1) + I2                                             
      I3I4 = IJADD(I3) + I4                                             
C                                                                       
                 NPACK = ICOUNT                                         
                 IPACK = I1I2                                           
                 JPACK = I3I4                                           
                 IF (LABSIZ .EQ. 2) THEN                                
*I32               IX( 2*NPACK-1 ) = IPACK                              
*I32               IX( 2*NPACK   ) = JPACK                              
                   LABEL = ISHFT( IPACK, 32 ) + JPACK                   
                   IX(NPACK) = LABEL                                    
                 ELSE IF (LABSIZ .EQ. 1) THEN                           
*I32               LABEL = ISHFT( IPACK, 16 ) + JPACK                   
*I32               IX(NPACK) = LABEL                                    
                   IF ( MOD(NPACK,2) .EQ. 0 ) THEN                      
                     LABEL = ISHFT( IPACK, 16 ) + JPACK                 
                     IX( NPACK/2 ) = IX( NPACK/2 ) + LABEL              
                   ELSE                                                 
                     LABEL = ISHFT( IPACK, 48 ) + ISHFT( JPACK, 32 )    
                     IX( (NPACK/2)+1 ) = LABEL                          
                   END IF                                               
                 END IF                                                 
C                                                                       
      BUFP(ICOUNT) = VALP                                               
      ICOUNT = ICOUNT+1                                                 
      IF (ICOUNT .GT. NINTMX) THEN                                      
         NXX = NINTMX                                                   
         CALL PWRIT(IS,BUFP,IX,NXX,NINTMX)                              
         ICOUNT = 1                                                     
         NREC = NREC+1                                                  
      END IF                                                            
      GO TO 780                                                         
C                                                                       
C     ----- -P- AND -K- SUPERMATRICES -----                             
C                                                                       
  740 CONTINUE                                                          
      IF (ABS(VALP) .LT. CUTOFF .AND. ABS(VALK) .LT. CUTOFF) GO TO 780  
      IF(OUT) THEN                                                      
         CALL INTOUT(I1,I2,I3,I4,QQ4,NN,VALP)                           
         CALL INTOUT(I1,I2,I3,I4,QQ4,NN,VALK)                           
      END IF                                                            
      IF (I1.EQ.I3 .AND. I2.EQ.I4) THEN                                 
         VALP = VALP*HALF                                               
         VALK = VALK*HALF                                               
      END IF                                                            
      I1I2 = IJADD(I1) + I2                                             
      I3I4 = IJADD(I3) + I4                                             
C                                                                       
                 NPACK = ICOUNT                                         
                 IPACK = I1I2                                           
                 JPACK = I3I4                                           
                 IF (LABSIZ .EQ. 2) THEN                                
*I32               IX( 2*NPACK-1 ) = IPACK                              
*I32               IX( 2*NPACK   ) = JPACK                              
                   LABEL = ISHFT( IPACK, 32 ) + JPACK                   
                   IX(NPACK) = LABEL                                    
                 ELSE IF (LABSIZ .EQ. 1) THEN                           
*I32               LABEL = ISHFT( IPACK, 16 ) + JPACK                   
*I32               IX(NPACK) = LABEL                                    
                   IF ( MOD(NPACK,2) .EQ. 0 ) THEN                      
                     LABEL = ISHFT( IPACK, 16 ) + JPACK                 
                     IX( NPACK/2 ) = IX( NPACK/2 ) + LABEL              
                   ELSE                                                 
                     LABEL = ISHFT( IPACK, 48 ) + ISHFT( JPACK, 32 )    
                     IX( (NPACK/2)+1 ) = LABEL                          
                   END IF                                               
                 END IF                                                 
C                                                                       
      BUFP(ICOUNT) = VALP                                               
      BUFK(ICOUNT) = VALK                                               
      ICOUNT = ICOUNT+1                                                 
      IF (ICOUNT .GT. NINTMX) THEN                                      
         NXX = NINTMX                                                   
         CALL PKWRIT(IS,BUFP,BUFK,IX,NXX,NINTMX)                        
         ICOUNT = 1                                                     
         NREC = NREC+1                                                  
      END IF                                                            
  780 CONTINUE                                                          
      GO TO (500,600,800),JUMP                                          
  800 CONTINUE                                                          
  820 CONTINUE                                                          
  840 CONTINUE                                                          
  860 CONTINUE                                                          
      IF(LRINT)IS=ISSAVE                                                
      RETURN                                                            
      END                                                               
C*MODULE INT2A   *DECK QOUT                                             
      SUBROUTINE QOUT(BUFP,IX,NINTMX,GHONDO)                            
      USE lrcdft, ONLY: LRFILE                                          
      use mx_limits, only: mxsh,mxgtot                                  
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      LOGICAL IANDJ,KANDL,SAME,OUT,LRINT                                
C                                                                       
*I32  CHARACTER*8 INAOFL                                                
C                                                                       
      DIMENSION BUFP(NINTMX),IX(*),GHONDO(*)                            
C                                                                       
C                                                                       
*I32  COMMON /ELGFIL/ NFILE,INAOFL                                      
      COMMON /ELGPMT/ NELONG,NATM,NASPIN,NCT,NBNDAB,NTMLB,IPRI,LDOS     
      COMMON /ERIOUT/ ISH,JSH,KSH,LSH,LSTRI,LSTRJ,LSTRK,LSTRL           
C$omp threadprivate(/ERIOUT/)
      COMMON /INT2IC/ NINTIC,ININTIC,NXXIC,LBUFPIC,LIXIC,LABSIX,NINTIX  
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)                
      COMMON /MISC  / IANDJ,KANDL,SAME                                  
C$omp threadprivate(/MISC/)
      COMMON /NLRCF / LRINT                                             
C$omp threadprivate(/NLRCF /)
      COMMON /NSHEL / EXX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),     
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),      
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),   
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL           
      COMMON /PCKLAB/ LABSIZ                                            
      COMMON /RESTAR/ TIMLIM,IREST,NREC,INTLOC,IST,JST,KST,LST          
      COMMON /SHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,          
     *                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,          
     *                NIJ,IJ,KL,IJKL                                    
C$omp threadprivate(/SHLNOS/)
      COMMON /SHLT  / TOL,CUTOFF,ICOUNT,OUT                             
C                                                                       
C-NEXT STATEMENT IS FOR VARIOUS IBM XLF 3.X AND 5.X COMPILERS-          
C                                                                       
      SAVE IJN,KLN                                                      
C                                                                       
      DATA HALF /0.5D+00/                                               
      ISSAVE=IS                                                         
      IF(LRINT)IS=LRFILE                                                
C                                                                       
C     ----- PACK THE 4 INDICES OF INTEGRAL INTO ONE WORD                
C     ----- WRITE LABEL + INTEGRAL ON TAPE (IS)                         
C                                                                       
      SAME  = ISH .EQ. KSH .AND. JSH .EQ. LSH                           
      IANDJ = ISH .EQ. JSH                                              
      KANDL = KSH .EQ. LSH                                              
C                                                                       
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
C                                                                       
      IJN = 0                                                           
      JMAX = MAXJ                                                       
      DO 260 I = MINI,MAXI                                              
         I_INDEX = (I-MINI)*LSTRI + 1                                   
         IF (IANDJ) JMAX = I                                            
         DO 240 J = MINJ,JMAX                                           
            IJ_INDEX = (J-MINJ)*LSTRJ + I_INDEX                         
            IJN = IJN+1                                                 
            LMAX = MAXL                                                 
            KLN = 0                                                     
            DO 220 K =  MINK,MAXK                                       
               IJK_INDEX = (K-MINK)*LSTRK + IJ_INDEX                    
               IF (KANDL) LMAX = K                                      
               DO 200 L = MINL,LMAX                                     
                  KLN = KLN+1                                           
                  IF(SAME  .AND.  KLN.GT.IJN) GO TO 240                 
                  IJKL_INDEX = (L-MINL)*LSTRL + IJK_INDEX               
C                                                                       
                  VAL = GHONDO( IJKL_INDEX )                            
                  IF(ABS(VAL).LT.CUTOFF) GO TO 200                      
C                                                                       
                  I1 = LOCI+I                                           
                  I2 = LOCJ+J                                           
                  I3 = LOCK+K                                           
                  I4 = LOCL+L                                           
                  IF (I1 .GE. I2) GO TO 100                             
                  N = I1                                                
                  I1 = I2                                               
                  I2 = N                                                
  100             IF (I3 .GE. I4) GO TO 120                             
                  N = I3                                                
                  I3 = I4                                               
                  I4 = N                                                
  120             IF (I1-I3) 140,160,180                                
  140             N = I1                                                
                  I1 = I3                                               
                  I3 = N                                                
                  N = I2                                                
                  I2 = I4                                               
                  I4 = N                                                
                  GO TO 180                                             
  160             IF (I2 .LT. I4) GO TO 140                             
  180             CONTINUE                                              
C                                                                       
                  IF (OUT) CALL INTOUT(I1,I2,I3,I4,QQ4,IJKL_INDEX,VAL)  
                  IF (I1 .EQ. I2) VAL = VAL*HALF                        
                  IF (I3 .EQ. I4) VAL = VAL*HALF                        
                  IF (I1 .EQ. I3 .AND. I2 .EQ. I4) VAL = VAL*HALF       
C                                                                       
                  NPACK = ICOUNT                                        
                  IPACK = I1                                            
                  JPACK = I2                                            
                  KPACK = I3                                            
                  LPACK = I4                                            
                  IF(LABSIZ .EQ. 2) THEN                                
*I32                 LABEL1 = ISHFT( IPACK, 16 ) + JPACK                
*I32                 LABEL2 = ISHFT( KPACK, 16 ) + LPACK                
*I32                 IX( 2*NPACK-1 ) = LABEL1                           
*I32                 IX( 2*NPACK   ) = LABEL2                           
                     LABEL = ISHFT( IPACK, 48 ) + ISHFT( JPACK, 32 ) +  
     *                       ISHFT( KPACK, 16 ) + LPACK                 
                     IX(NPACK) = LABEL                                  
                  ELSE IF (LABSIZ .EQ. 1) THEN                          
*I32                 LABEL = ISHFT( IPACK, 24 ) + ISHFT( JPACK, 16 ) +  
*I32 *                       ISHFT( KPACK,  8 ) + LPACK                 
*I32                 IX(NPACK) = LABEL                                  
                     IF ( MOD(NPACK,2) .EQ. 0 ) THEN                    
                       LABEL = ISHFT( IPACK, 24 ) + ISHFT( JPACK, 16 ) +
     *                         ISHFT( KPACK,  8 ) + LPACK               
                       IX( NPACK/2 ) = IX( NPACK/2 ) + LABEL            
                     ELSE                                               
                       LABEL = ISHFT( IPACK, 56 ) + ISHFT( JPACK, 48 ) +
     *                         ISHFT( KPACK, 40 ) + ISHFT( LPACK, 32 )  
                       IX( (NPACK/2)+1 ) = LABEL                        
                     END IF                                             
                  END IF                                                
C                                                                       
                  BUFP(ICOUNT) = VAL                                    
                  ICOUNT = ICOUNT+1                                     
                  IF(ICOUNT .GT. NINTIC) THEN                           
                     JCOUNT=ICOUNT-NINTIC                               
                     IF(JCOUNT .GT. NINTMX) THEN                        
                        NXX = NINTMX                                    
                        CALL PWRIT(IS,BUFP(NINTIC+1),IX(ININTIC+1),     
     *                             NXX,NINTMX)                          
                        ICOUNT = NINTIC+1                               
                        NREC = NREC+1                                   
C                                                                       
C             ELONGATION METHOD                                         
                        IF(NELONG.GT.1) THEN                            
*I32                    IF(NREC.GT.8650) THEN                           
*I32                       WRITE(IW,*)                                  
*I32 *                       'NEXT FILE WITH INTEGRALS IS OPENED'       
*I32                       NREC = 1                                     
*I32                       NFILE = NFILE + 1                            
*I32                       WRITE(IW,*)' NFILE = ', NFILE                
*I32                       CALL SEQCLO(IS,'KEEP')                       
*I32                       IF(NFILE.GT.9) THEN                          
*I32                          WRITE(INAOFL(7:8),'(I2)') NFILE           
*I32                       ELSE                                         
*I32                          WRITE(INAOFL(7:8),'(I1,'' '')') NFILE     
*I32                       ENDIF                                        
*I32                       CALL SEQOPN(IS,INAOFL,'UNKNOWN',.FALSE.,     
*I32 *                        'UNFORMATTED')                            
*I32                    ENDIF                                           
                       ENDIF                                            
C                                                                       
                     END IF ! JCOUNT                                    
                  END IF ! ICOUNT                                       
  200          CONTINUE                                                 
  220       CONTINUE                                                    
  240    CONTINUE                                                       
  260 CONTINUE                                                          
      IF(LRINT)IS=ISSAVE                                                
      RETURN                                                            
      END                                                               
C*MODULE INT2A   *DECK S0000                                            
      SUBROUTINE S0000(GHONDO,DDIJ)                                     
      USE lrcdft, ONLY: EMU2                                            
      use mx_limits, only: mxgsh,mxg2                                   
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      LOGICAL IANDJ,KANDL,SAME,OUT,LRINT                                
C                                                                       
C                                                                       
      DIMENSION GHONDO(*),DDIJ(49*MXG2)                                 
C                                                                       
      COMMON /IJGNRL/ A(MXG2),R(MXG2),X1(MXG2),Y1(MXG2),Z1(MXG2),       
     *                IJD(784)                                          
C$omp threadprivate(/IJGNRL/)
      COMMON /MISC  / IANDJ,KANDL,SAME                                  
C$omp threadprivate(/MISC/)
      COMMON /NLRCF / LRINT                                             
C$omp threadprivate(/NLRCF /)
      COMMON /SHLINF/  AG(MXGSH),CSA(MXGSH),CPA(MXGSH),CDA(MXGSH),      
     *                CFA(MXGSH),CGA(MXGSH),CHA(MXGSH),CIA(MXGSH),      
     *                 BG(MXGSH),CSB(MXGSH),CPB(MXGSH),CDB(MXGSH),      
     *                CFB(MXGSH),CGB(MXGSH),CHB(MXGSH),CIB(MXGSH),      
     *                 CG(MXGSH),CSC(MXGSH),CPC(MXGSH),CDC(MXGSH),      
     *                CFC(MXGSH),CGC(MXGSH),CHC(MXGSH),CIC(MXGSH),      
     *                 DG(MXGSH),CSD(MXGSH),CPD(MXGSH),CDD(MXGSH),      
     *                CFD(MXGSH),CGD(MXGSH),CHD(MXGSH),CID(MXGSH),      
     *                XI,YI,ZI,XJ,YJ,ZJ,RRI,XK,YK,ZK,XL,YL,ZL,RRK,      
     *                NGA,NGB,NGC,NGD                                   
C$omp threadprivate(/SHLINF/)
      COMMON /SHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,          
     *                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,          
     *                NIJ,IJ,KL,IJKL                                    
C$omp threadprivate(/SHLNOS/)
      COMMON /SHLT  / TOL,CUTOFF,ICOUNT,OUT                             
C                                                                       
      PARAMETER (PI252=34.986836655250D+00)                             
      PARAMETER (PIE4=7.85398163397448D-01)                             
      PARAMETER (ZERO=0.0D+00)                                          
      PARAMETER (ONE=1.0D+00)                                           
C                                                                       
C     SPECIAL SSSS INTEGRAL ROUTINE WHEN USING HONDO INTEGRALS          
C                                                                       
      GGOUT = ZERO                                                      
      LGMAX = NGD                                                       
      DO 300 KG = 1,NGC                                                 
      BK = CG(KG)                                                       
      BRRK = BK*RRK                                                     
      BXK = BK*XK                                                       
      BYK = BK*YK                                                       
      BZK = BK*ZK                                                       
      CSK = CSC(KG)                                                     
      IF (KANDL) LGMAX = KG                                             
      DO 280 LG = 1,LGMAX                                               
      BL = DG(LG)                                                       
      BB = BK+BL                                                        
      BBINV = ONE/BB                                                    
      DUM = BL*BRRK*BBINV                                               
      IF (DUM .GT. TOL) GO TO 280                                       
      BBRRK = DUM                                                       
      D2 = CSD(LG)*CSK*BBINV                                            
      IF (KANDL .AND. LG .NE. KG) D2 = D2+D2                            
      BBX = (BXK+BL*XL)*BBINV                                           
      BBY = (BYK+BL*YL)*BBINV                                           
      BBZ = (BZK+BL*ZL)*BBINV                                           
      SUM = ZERO                                                        
      NN = 1                                                            
      DO 260 N = 1,NIJ                                                  
      DUM = BBRRK+R(N)                                                  
      IF (DUM .GT. TOL) GO TO 260                                       
      EXPE = EXP(-DUM)                                                  
      AA = A(N)                                                         
      AB = AA+BB                                                        
      DUM = X1(N)-BBX                                                   
      XX = DUM*DUM                                                      
      DUM = Y1(N)-BBY                                                   
      XX = DUM*DUM+XX                                                   
      DUM = Z1(N)-BBZ                                                   
      XX = DUM*DUM+XX                                                   
      X = XX*AA*BB/AB                                                   
      IF(LRINT) THEN                                                    
         RHO = AA*BB*EMU2/(AB*EMU2+AA*BB)                               
         X   = XX*RHO                                                   
      ENDIF                                                             
C                                                                       
      IF (X .GT. 5.0D+00) GO TO 160                                     
      IF (X .GT. 1.0D+00) GO TO 120                                     
      IF (X .GT. 3.0D-07) GO TO 100                                     
      WW1 = 1.0D+00-X/3.0D+00                                           
      GO TO 240                                                         
C                                                                       
  100 CONTINUE                                                          
      F1 = ((((((((-8.36313918003957D-08*X+1.21222603512827D-06 )*X-    
     +     1.15662609053481D-05 )*X+9.25197374512647D-05 )*X-           
     +     6.40994113129432D-04 )*X+3.78787044215009D-03 )*X-           
     +     1.85185172458485D-02 )*X+7.14285713298222D-02 )*X-           
     +     1.99999999997023D-01 )*X+3.33333333333318D-01                
      WW1 = (X+X)*F1+EXP(-X)                                            
      GO TO 240                                                         
C                                                                       
  120 CONTINUE                                                          
      IF (X .GT. 3.0D+00) GO TO 140                                     
      Y = X-2.0D+00                                                     
      F1 = ((((((((((-1.61702782425558D-10*Y+1.96215250865776D-09 )*Y-  
     +     2.14234468198419D-08 )*Y+2.17216556336318D-07 )*Y-           
     +     1.98850171329371D-06 )*Y+1.62429321438911D-05 )*Y-           
     +     1.16740298039895D-04 )*Y+7.24888732052332D-04 )*Y-           
     +     3.79490003707156D-03 )*Y+1.61723488664661D-02 )*Y-           
     +     5.29428148329736D-02 )*Y+1.15702180856167D-01                
      WW1 = (X+X)*F1+EXP(-X)                                            
      GO TO 240                                                         
C                                                                       
  140 CONTINUE                                                          
      Y = X-4.0D+00                                                     
      F1 = ((((((((((-2.62453564772299D-11*Y+3.24031041623823D-10 )*Y-  
     +     3.614965656163D-09)*Y+3.760256799971D-08)*Y-                 
     +     3.553558319675D-07)*Y+3.022556449731D-06)*Y-                 
     +     2.290098979647D-05)*Y+1.526537461148D-04)*Y-                 
     +     8.81947375894379D-04 )*Y+4.33207949514611D-03 )*Y-           
     +     1.75257821619926D-02 )*Y+5.28406320615584D-02                
      WW1 = (X+X)*F1+EXP(-X)                                            
      GO TO 240                                                         
C                                                                       
  160 CONTINUE                                                          
      IF (X .GT. 15.0D+00) GO TO 200                                    
      E = EXP(-X)                                                       
      IF (X .GT. 10.0D+00) GO TO 180                                    
      XINV = ONE/X                                                      
      WW1 = (((((( 4.6897511375022D-01*XINV-6.9955602298985D-01)*XINV + 
     +     5.3689283271887D-01)*XINV-3.2883030418398D-01)*XINV +        
     +     2.4645596956002D-01)*XINV-4.9984072848436D-01)*XINV -        
     +     3.1501078774085D-06)*E + SQRT(PIE4*XINV)                     
      GO TO 240                                                         
C                                                                       
  180 CONTINUE                                                          
      XINV = ONE/X                                                      
      WW1 = (((-1.8784686463512D-01*XINV+2.2991849164985D-01)*XINV      
     +         -4.9893752514047D-01)*XINV-2.1916512131607D-05)*E        
     +         + SQRT(PIE4*XINV)                                        
      GO TO 240                                                         
C                                                                       
  200 CONTINUE                                                          
      IF (X .GT. 33.0D+00) GO TO 220                                    
      XINV = ONE/X                                                      
      E = EXP(-X)                                                       
      WW1 = (( 1.9623264149430D-01*XINV-4.9695241464490D-01)*XINV -     
     +     6.0156581186481D-05)*E + SQRT(PIE4*XINV)                     
      GO TO 240                                                         
C                                                                       
  220 WW1 = SQRT(PIE4/X)                                                
C                                                                       
  240 CONTINUE                                                          
      IF(.NOT.LRINT)SUM = SUM+DDIJ(NN)*WW1*EXPE/SQRT(AB)                
      IF(     LRINT)SUM = SUM+DDIJ(NN)*WW1*EXPE/SQRT(AA*BB/RHO)         
  260 NN = NN+49                                                        
      GGOUT = GGOUT+D2*SUM                                              
  280 CONTINUE                                                          
  300 CONTINUE                                                          
      GHONDO(1) = GGOUT*PI252*QQ4                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2A   *DECK SCHWDN                                           
      DOUBLE PRECISION FUNCTION SCHWDN(DSH,ISH,JSH,KSH,LSH,IA)          
C                                                                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C                                                                       
      DIMENSION DSH(*),IA(*)                                            
C                                                                       
      PARAMETER (FOUR=4.0D+00)                                          
C                                                                       
C     ----- FIND MAXIMUM DENSITY CONTRIBUTION TO THIS SHELL SET -----   
C     -DSH- IS THE DENSITY MATRIX ALREADY COMPRESSED TO SHELLS          
C                                                                       
      IJ = IA(ISH)+JSH                                                  
      IK = IA(ISH)+KSH                                                  
      IL = IA(ISH)+LSH                                                  
      KL = IA(KSH)+LSH                                                  
      JK = IA(JSH)+KSH                                                  
      JL = IA(JSH)+LSH                                                  
      IF(JSH.LT.KSH) JK=IA(KSH)+JSH                                     
      IF(JSH.LT.LSH) JL=IA(LSH)+JSH                                     
      SCHWDN=MAX(FOUR*DSH(IJ),FOUR*DSH(KL),                             
     *           DSH(JL),DSH(JK),DSH(IL),DSH(IK))                       
      RETURN                                                            
      END                                                               
C*MODULE INT2A   *DECK SHELLS                                           
      SUBROUTINE SHELLS(NELEC,ISH,JSH,KSH,LSH,FLIP)                     
      use mx_limits, only: mxsh,mxgsh,mxgtot,mxatm                      
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      LOGICAL FLIP                                                      
      LOGICAL IANDJ,KANDL,SAME                                          
C                                                                       
      DIMENSION IX(84),IY(84),IZ(84),                                   
     *          JX(84),JY(84),JZ(84),                                   
     *          KX(84),KY(84),KZ(84),                                   
     *          LX(84),LY(84),LZ(84)                                    
C                                                                       
C                                                                       
      COMMON /ERIOUT/ INU,JNU,KNU,LNU,NGTI,NGTJ,NGTK,NGTL               
C$omp threadprivate(/ERIOUT/)
      COMMON /INTDEX/ IJGT(784),IJX(784),IJY(784),IJZ(784),IK(784),     
     *                KLGT(784),KLX(784),KLY(784),KLZ(784)              
C$omp threadprivate(/INTDEX/)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,                    
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)                  
      COMMON /MISC  / IANDJ,KANDL,SAME                                  
C$omp threadprivate(/MISC/)
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),      
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),      
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),   
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL           
      COMMON /ROOT  / XX,U(13),W(13),NROOTS                             
C$omp threadprivate(/ROOT/)
      COMMON /SHLEXC/ NORGSH(3),NORGSP(3),IEXCH,NANGM,NGTH(4)           
      COMMON /SHLINF/  GA(MXGSH),CSA(MXGSH),CPA(MXGSH),CDA(MXGSH),      
     *                CFA(MXGSH),CGA(MXGSH),CHA(MXGSH),CIA(MXGSH),      
     *                 GB(MXGSH),CSB(MXGSH),CPB(MXGSH),CDB(MXGSH),      
     *                CFB(MXGSH),CGB(MXGSH),CHB(MXGSH),CIB(MXGSH),      
     *                 GC(MXGSH),CSC(MXGSH),CPC(MXGSH),CDC(MXGSH),      
     *                CFC(MXGSH),CGC(MXGSH),CHC(MXGSH),CIC(MXGSH),      
     *                 GD(MXGSH),CSD(MXGSH),CPD(MXGSH),CDD(MXGSH),      
     *                CFD(MXGSH),CGD(MXGSH),CHD(MXGSH),CID(MXGSH),      
     *                AX,AY,AZ,BX,BY,BZ,RAB,CX,CY,CZ,DX,DY,DZ,RCD,      
     *                NGA,NGB,NGC,NGD                                   
C$omp threadprivate(/SHLINF/)
      COMMON /SHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,          
     +                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,          
     +                NIJ,IJ,KL,IJKL                                    
C$omp threadprivate(/SHLNOS/)
C                                                                       
      DATA LX /   0,   1,   0,   0,   2,   0,   0,   1,   1,   0,       
     *            3,   0,   0,   2,   2,   1,   0,   1,   0,   1,       
     *            4,   0,   0,   3,   3,   1,   0,   1,   0,   2,       
     *            2,   0,   2,   1,   1,                                
     *            5,   0,   0,   4,   4,   1,   0,   1,   0,   3,       
     *            3,   2,   0,   2,   0,   3,   1,   1,   2,   2,       
     *            1,                                                    
     *            6,   0,   0,   5,   5,   1,   0,   1,   0,   4,       
     *            4,   2,   0,   2,   0,   4,   1,   1,   3,   3,       
     *            0,   3,   3,   2,   1,   2,   1,   2/                 
      DATA KX /   0,   7,   0,   0,  14,   0,   0,   7,   7,   0,       
     *           21,   0,   0,  14,  14,   7,   0,   7,   0,   7,       
     *           28,   0,   0,  21,  21,   7,   0,   7,   0,  14,       
     *           14,   0,  14,   7,   7,                                
     *           35,   0,   0,  28,  28,   7,   0,   7,   0,  21,       
     *           21,  14,   0,  14,   0,  21,   7,   7,  14,  14,       
     *            7,                                                    
     *           42,   0,   0,  35,  35,   7,   0,   7,   0,  28,       
     *           28,  14,   0,  14,   0,  28,   7,   7,  21,  21,       
     *            0,  21,  21,  14,   7,  14,   7,  14/                 
      DATA JX /   0,  49,   0,   0,  98,   0,   0,  49,  49,   0,       
     *          147,   0,   0,  98,  98,  49,   0,  49,   0,  49,       
     *          196,   0,   0, 147, 147,  49,   0,  49,   0,  98,       
     *           98,   0,  98,  49,  49,                                
     *          245,   0,   0, 196, 196,  49,   0,  49,   0, 147,       
     *          147,  98,   0,  98,   0, 147,  49,  49,  98,  98,       
     *           49,                                                    
     *          294,   0,   0, 245, 245,  49,   0,  49,   0, 196,       
     *          196,  98,   0,  98,   0, 196,  49,  49, 147, 147,       
     *            0, 147, 147,  98,  49,  98,  49,  98/                 
      DATA IX /   1, 344,   1,   1, 687,   1,   1, 344, 344,   1,       
     *         1030,   1,   1, 687, 687, 344,   1, 344,   1, 344,       
     *         1373,   1,   1,1030,1030, 344,   1, 344,   1, 687,       
     *          687,   1, 687, 344, 344,                                
     *         1716,   1,   1,1373,1373, 344,   1, 344,   1,1030,       
     *         1030, 687,   1, 687,   1,1030, 344, 344, 687, 687,       
     *          344,                                                    
     *         2059,   1,   1,1716,1716, 344,   1, 344,   1,1373,       
     *         1373, 687,   1, 687,   1,1373, 344, 344,1030,1030,       
     *            1,1030,1030, 687, 344, 687, 344, 687/                 
      DATA LY /   0,   0,   1,   0,   0,   2,   0,   1,   0,   1,       
     *            0,   3,   0,   1,   0,   2,   2,   0,   1,   1,       
     *            0,   4,   0,   1,   0,   3,   3,   0,   1,   2,       
     *            0,   2,   1,   2,   1,                                
     *            0,   5,   0,   1,   0,   4,   4,   0,   1,   2,       
     *            0,   3,   3,   0,   2,   1,   3,   1,   2,   1,       
     *            2,                                                    
     *            0,   6,   0,   1,   0,   5,   5,   0,   1,   2,       
     *            0,   4,   4,   0,   2,   1,   4,   1,   3,   0,       
     *            3,   2,   1,   3,   3,   1,   2,   2/                 
      DATA KY /   0,   0,   7,   0,   0,  14,   0,   7,   0,   7,       
     *            0,  21,   0,   7,   0,  14,  14,   0,   7,   7,       
     *            0,  28,   0,   7,   0,  21,  21,   0,   7,  14,       
     *            0,  14,   7,  14,   7,                                
     *            0,  35,   0,   7,   0,  28,  28,   0,   7,  14,       
     *            0,  21,  21,   0,  14,   7,  21,   7,  14,   7,       
     *           14,                                                    
     *            0,  42,   0,   7,   0,  35,  35,   0,   7,  14,       
     *            0,  28,  28,   0,  14,   7,  28,   7,  21,   0,       
     *           21,  14,   7,  21,  21,   7,  14,  14/                 
      DATA JY /   0,   0,  49,   0,   0,  98,   0,  49,   0,  49,       
     *            0, 147,   0,  49,   0,  98,  98,   0,  49,  49,       
     *            0, 196,   0,  49,   0, 147, 147,   0,  49,  98,       
     *            0,  98,  49,  98,  49,                                
     *            0, 245,   0,  49,   0, 196, 196,   0,  49,  98,       
     *            0, 147, 147,   0,  98,  49, 147,  49,  98,  49,       
     *           98,                                                    
     *            0, 294,   0,  49,   0, 245, 245,   0,  49,  98,       
     *            0, 196, 196,   0,  98,  49, 196,  49, 147,   0,       
     *          147,  98,  49, 147, 147,  49,  98,  98/                 
      DATA IY /   1,   1, 344,   1,   1, 687,   1, 344,   1, 344,       
     *            1,1030,   1, 344,   1, 687, 687,   1, 344, 344,       
     *            1,1373,   1, 344,   1,1030,1030,   1, 344, 687,       
     *            1, 687, 344, 687, 344,                                
     *            1,1716,   1, 344,   1,1373,1373,   1, 344, 687,       
     *            1,1030,1030,   1, 687, 344,1030, 344, 687, 344,       
     *          687,                                                    
     *            1,2059,   1, 344,   1,1716,1716,   1, 344, 687,       
     *            1,1373,1373,   1, 687, 344,1373, 344,1030,   1,       
     *         1030, 687, 344,1030,1030, 344, 687, 687/                 
      DATA LZ /   0,   0,   0,   1,   0,   0,   2,   0,   1,   1,       
     *            0,   0,   3,   0,   1,   0,   1,   2,   2,   1,       
     *            0,   0,   4,   0,   1,   0,   1,   3,   3,   0,       
     *            2,   2,   1,   1,   2,                                
     *            0,   0,   5,   0,   1,   0,   1,   4,   4,   0,       
     *            2,   0,   2,   3,   3,   1,   1,   3,   1,   2,       
     *            2,                                                    
     *            0,   0,   6,   0,   1,   0,   1,   5,   5,   0,       
     *            2,   0,   2,   4,   4,   1,   1,   4,   0,   3,       
     *            3,   1,   2,   1,   2,   3,   3,   2/                 
      DATA KZ /   0,   0,   0,   7,   0,   0,  14,   0,   7,   7,       
     *            0,   0,  21,   0,   7,   0,   7,  14,  14,   7,       
     *            0,   0,  28,   0,   7,   0,   7,  21,  21,   0,       
     *           14,  14,   7,   7,  14,                                
     *            0,   0,  35,   0,   7,   0,   7,  28,  28,   0,       
     *           14,   0,  14,  21,  21,   7,   7,  21,   7,  14,       
     *           14,                                                    
     *            0,   0,  42,   0,   7,   0,   7,  35,  35,   0,       
     *           14,   0,  14,  28,  28,   7,   7,  28,   0,  21,       
     *           21,   7,  14,   7,  14,  21,  21,  14/                 
      DATA JZ /   0,   0,   0,  49,   0,   0,  98,   0,  49,  49,       
     *            0,   0, 147,   0,  49,   0,  49,  98,  98,  49,       
     *            0,   0, 196,   0,  49,   0,  49, 147, 147,   0,       
     *           98,  98,  49,  49,  98,                                
     *            0,   0, 245,   0,  49,   0,  49, 196, 196,   0,       
     *           98,   0,  98, 147, 147,  49,  49, 147,  49,  98,       
     *           98,                                                    
     *            0,   0, 294,   0,  49,   0,  49, 245, 245,   0,       
     *           98,   0,  98, 196, 196,  49,  49, 196,   0, 147,       
     *          147,  49,  98,  49,  98, 147, 147,  98/                 
      DATA IZ /   1,   1,   1, 344,   1,   1, 687,   1, 344, 344,       
     *            1,   1,1030,   1, 344,   1, 344, 687, 687, 344,       
     *            1,   1,1373,   1, 344,   1, 344,1030,1030,   1,       
     *          687, 687, 344, 344, 687,                                
     *            1,   1,1716,   1, 344,   1, 344,1373,1373,   1,       
     *          687,   1, 687,1030,1030, 344, 344,1030, 344, 687,       
     *          687,                                                    
     *            1,   1,2059,   1, 344,   1, 344,1716,1716,   1,       
     *          687,   1, 687,1373,1373, 344, 344,1373,   1,1030,       
     *         1030, 344, 687, 344, 687,1030,1030, 687/                 
C                                                                       
C     PREPARE SHELL INFORMATION/FOR HONDO INTEGRATION                   
C                                                                       
      IF(NELEC.EQ.2) GO TO 200                                          
C                                                                       
C     ----- PERMUTE ISH AND JSH SHELLS, FOR THEIR TYPE                  
C     THIS IS DONE FOR SPEED REASONS.  THE CODE GETS THE RIGHT ANSWER   
C     WITHOUT THE ANGULAR MOMENTUM FLIPPING, AND THEREFORE A CALLING    
C     ARGUMENT ALLOWS ONE DO EXACTLY THE INTEGRAL BLOCK AS SPECIFIED,   
C     SHOULD THAT BE DESIRED.                                           
C                                                                       
      IANDJ = ISH .EQ. JSH                                              
      IF (KTYPE(ISH) .LT. KTYPE(JSH)  .AND.  FLIP) THEN                 
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
C                                                                       
C     ----- ISHELL                                                      
C                                                                       
      I = KATOM(INU)                                                    
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
      DO 140 I = I1,I2                                                  
         NGA = NGA+1                                                    
         GA(NGA) = EX(I)                                                
         CSA(NGA) = CS(I)                                               
         CPA(NGA) = CP(I)                                               
         CDA(NGA) = CD(I)                                               
         CFA(NGA) = CF(I)                                               
         CGA(NGA) = CG(I)                                               
         CHA(NGA) = CH(I)                                               
         CIA(NGA) = CI(I)                                               
  140 CONTINUE                                                          
C                                                                       
C     ----- JSHELL                                                      
C                                                                       
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
      DO 160 J = J1,J2                                                  
         NGB = NGB+1                                                    
         GB(NGB) = EX(J)                                                
         CSB(NGB) = CS(J)                                               
         CPB(NGB) = CP(J)                                               
         CDB(NGB) = CD(J)                                               
         CFB(NGB) = CF(J)                                               
         CGB(NGB) = CG(J)                                               
         CHB(NGB) = CH(J)                                               
         CIB(NGB) = CI(J)                                               
  160 CONTINUE                                                          
      RAB = ((AX-BX)*(AX-BX) + (AY-BY)*(AY-BY) + (AZ-BZ)*(AZ-BZ))       
C                                                                       
C     ----- PREPARE INDICES FOR PAIRS OF (I,J) FUNCTIONS                
C                                                                       
      IJ = 0                                                            
      JMAX = MAXJ                                                       
      DO 190 I = MINI,MAXI                                              
         NX = IX(I)                                                     
         NY = IY(I)                                                     
         NZ = IZ(I)                                                     
         IF (IANDJ) JMAX = I                                            
         DO 180 J = MINJ,JMAX                                           
            IJ = IJ+1                                                   
            IJX(IJ) = NX+JX(J)                                          
            IJY(IJ) = NY+JY(J)                                          
            IJZ(IJ) = NZ+JZ(J)                                          
            IJGT(IJ) = NGTI*(I-MINI)+NGTJ*(J-MINJ)+1                    
  180    CONTINUE                                                       
  190 CONTINUE                                                          
      RETURN                                                            
C     ******                                                            
C                                                                       
C        K AND L SHELL                                                  
C                                                                       
  200 CONTINUE                                                          
      KANDL = KSH .EQ. LSH                                              
      SAME = ISH .EQ. KSH .AND. JSH .EQ. LSH                            
C                                                                       
C     ----- PERMUTE KSH AND LSH SHELLS, FOR THEIR TYPE                  
C                                                                       
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
C                                                                       
C     ----- K SHELL                                                     
C                                                                       
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
      DO 260 K = K1,K2                                                  
         NGC = NGC+1                                                    
         GC(NGC) = EX(K)                                                
         CSC(NGC) = CS(K)                                               
         CPC(NGC) = CP(K)                                               
         CDC(NGC) = CD(K)                                               
         CFC(NGC) = CF(K)                                               
         CGC(NGC) = CG(K)                                               
         CHC(NGC) = CH(K)                                               
         CIC(NGC) = CI(K)                                               
  260 CONTINUE                                                          
C                                                                       
C     ----- LSHELL                                                      
C                                                                       
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
      DO 280 L = L1,L2                                                  
         NGD = NGD+1                                                    
         GD(NGD) = EX(L)                                                
         CSD(NGD) = CS(L)                                               
         CPD(NGD) = CP(L)                                               
         CDD(NGD) = CD(L)                                               
         CFD(NGD) = CF(L)                                               
         CGD(NGD) = CG(L)                                               
         CHD(NGD) = CH(L)                                               
         CID(NGD) = CI(L)                                               
  280 CONTINUE                                                          
      NROOTS = (LIT+LJT+LKT+LLT-2)/2                                    
      RCD = ((CX-DX)*(CX-DX) + (CY-DY)*(CY-DY) + (CZ-DZ)*(CZ-DZ))       
C                                                                       
C     ----- PREPARE INDICES FOR PAIRS OF (K,L) FUNCTIONS                
C                                                                       
      KL = 0                                                            
      LMAX = MAXL                                                       
      DO 310 K = MINK,MAXK                                              
         NX = KX(K)                                                     
         NY = KY(K)                                                     
         NZ = KZ(K)                                                     
         IF (KANDL) LMAX = K                                            
         DO 300 L = MINL,LMAX                                           
            KL = KL+1                                                   
            KLX(KL) = NX+LX(L)                                          
            KLY(KL) = NY+LY(L)                                          
            KLZ(KL) = NZ+LZ(L)                                          
            KLGT(KL) = NGTK*(K-MINK)+NGTL*(L-MINL)                      
  300    CONTINUE                                                       
  310 CONTINUE                                                          
      MAX = KL                                                          
      DO 320 I = 1,IJ                                                   
      IF (SAME) MAX = I                                                 
  320 IK(I) = MAX                                                       
      IJKL = IJ*KL                                                      
      IF (SAME) IJKL = IJ*(IJ+1)/2                                      
      RETURN                                                            
      END                                                               
                                                                        
C*MODULE INT2A   *DECK SHELLS                                           
      SUBROUTINE SHELLS_gpu2(NELEC,ISH,JSH,KSH,LSH,FLIP)                
      use mx_limits, only: mxsh,mxgsh,mxgtot,mxatm                      
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      LOGICAL FLIP                                                      
      LOGICAL IANDJ,KANDL,SAME                                          
C                                                                       
      DIMENSION IX(84),IY(84),IZ(84),                                   
     *          JX(84),JY(84),JZ(84),                                   
     *          KX(84),KY(84),KZ(84),                                   
     *          LX(84),LY(84),LZ(84)                                    
C                                                                       
C                                                                       
      COMMON /ERIOUT/ INU,JNU,KNU,LNU,NGTI,NGTJ,NGTK,NGTL               
C$omp threadprivate(/ERIOUT/)
      COMMON /INTDEX/ IJGT(784),IJX(784),IJY(784),IJZ(784),IK(784),     
     *                KLGT(784),KLX(784),KLY(784),KLZ(784)              
C$omp threadprivate(/INTDEX/)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,                    
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)                  
      COMMON /MISC  / IANDJ,KANDL,SAME                                  
C$omp threadprivate(/MISC/)
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),      
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),      
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),   
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL           
      COMMON /ROOT  / XX,U(13),W(13),NROOTS                             
C$omp threadprivate(/ROOT/)
      COMMON /SHLEXC/ NORGSH(3),NORGSP(3),IEXCH,NANGM,NGTH(4)           
      COMMON /SHLINF/  GA(MXGSH),CSA(MXGSH),CPA(MXGSH),CDA(MXGSH),      
     *                CFA(MXGSH),CGA(MXGSH),CHA(MXGSH),CIA(MXGSH),      
     *                 GB(MXGSH),CSB(MXGSH),CPB(MXGSH),CDB(MXGSH),      
     *                CFB(MXGSH),CGB(MXGSH),CHB(MXGSH),CIB(MXGSH),      
     *                 GC(MXGSH),CSC(MXGSH),CPC(MXGSH),CDC(MXGSH),      
     *                CFC(MXGSH),CGC(MXGSH),CHC(MXGSH),CIC(MXGSH),      
     *                 GD(MXGSH),CSD(MXGSH),CPD(MXGSH),CDD(MXGSH),      
     *                CFD(MXGSH),CGD(MXGSH),CHD(MXGSH),CID(MXGSH),      
     *                AX,AY,AZ,BX,BY,BZ,RAB,CX,CY,CZ,DX,DY,DZ,RCD,      
     *                NGA,NGB,NGC,NGD                                   
C$omp threadprivate(/SHLINF/)
      COMMON /SHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,          
     +                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,          
     +                NIJ,IJ,KL,IJKL                                    
C$omp threadprivate(/SHLNOS/)
C                                                                       
      DATA LX /   0,   1,   0,   0,   2,   0,   0,   1,   1,   0,       
     *            3,   0,   0,   2,   2,   1,   0,   1,   0,   1,       
     *            4,   0,   0,   3,   3,   1,   0,   1,   0,   2,       
     *            2,   0,   2,   1,   1,                                
     *            5,   0,   0,   4,   4,   1,   0,   1,   0,   3,       
     *            3,   2,   0,   2,   0,   3,   1,   1,   2,   2,       
     *            1,                                                    
     *            6,   0,   0,   5,   5,   1,   0,   1,   0,   4,       
     *            4,   2,   0,   2,   0,   4,   1,   1,   3,   3,       
     *            0,   3,   3,   2,   1,   2,   1,   2/                 
      DATA KX /   0,   7,   0,   0,  14,   0,   0,   7,   7,   0,       
     *           21,   0,   0,  14,  14,   7,   0,   7,   0,   7,       
     *           28,   0,   0,  21,  21,   7,   0,   7,   0,  14,       
     *           14,   0,  14,   7,   7,                                
     *           35,   0,   0,  28,  28,   7,   0,   7,   0,  21,       
     *           21,  14,   0,  14,   0,  21,   7,   7,  14,  14,       
     *            7,                                                    
     *           42,   0,   0,  35,  35,   7,   0,   7,   0,  28,       
     *           28,  14,   0,  14,   0,  28,   7,   7,  21,  21,       
     *            0,  21,  21,  14,   7,  14,   7,  14/                 
      DATA JX /   0,  49,   0,   0,  98,   0,   0,  49,  49,   0,       
     *          147,   0,   0,  98,  98,  49,   0,  49,   0,  49,       
     *          196,   0,   0, 147, 147,  49,   0,  49,   0,  98,       
     *           98,   0,  98,  49,  49,                                
     *          245,   0,   0, 196, 196,  49,   0,  49,   0, 147,       
     *          147,  98,   0,  98,   0, 147,  49,  49,  98,  98,       
     *           49,                                                    
     *          294,   0,   0, 245, 245,  49,   0,  49,   0, 196,       
     *          196,  98,   0,  98,   0, 196,  49,  49, 147, 147,       
     *            0, 147, 147,  98,  49,  98,  49,  98/                 
      DATA IX /   1, 344,   1,   1, 687,   1,   1, 344, 344,   1,       
     *         1030,   1,   1, 687, 687, 344,   1, 344,   1, 344,       
     *         1373,   1,   1,1030,1030, 344,   1, 344,   1, 687,       
     *          687,   1, 687, 344, 344,                                
     *         1716,   1,   1,1373,1373, 344,   1, 344,   1,1030,       
     *         1030, 687,   1, 687,   1,1030, 344, 344, 687, 687,       
     *          344,                                                    
     *         2059,   1,   1,1716,1716, 344,   1, 344,   1,1373,       
     *         1373, 687,   1, 687,   1,1373, 344, 344,1030,1030,       
     *            1,1030,1030, 687, 344, 687, 344, 687/                 
      DATA LY /   0,   0,   1,   0,   0,   2,   0,   1,   0,   1,       
     *            0,   3,   0,   1,   0,   2,   2,   0,   1,   1,       
     *            0,   4,   0,   1,   0,   3,   3,   0,   1,   2,       
     *            0,   2,   1,   2,   1,                                
     *            0,   5,   0,   1,   0,   4,   4,   0,   1,   2,       
     *            0,   3,   3,   0,   2,   1,   3,   1,   2,   1,       
     *            2,                                                    
     *            0,   6,   0,   1,   0,   5,   5,   0,   1,   2,       
     *            0,   4,   4,   0,   2,   1,   4,   1,   3,   0,       
     *            3,   2,   1,   3,   3,   1,   2,   2/                 
      DATA KY /   0,   0,   7,   0,   0,  14,   0,   7,   0,   7,       
     *            0,  21,   0,   7,   0,  14,  14,   0,   7,   7,       
     *            0,  28,   0,   7,   0,  21,  21,   0,   7,  14,       
     *            0,  14,   7,  14,   7,                                
     *            0,  35,   0,   7,   0,  28,  28,   0,   7,  14,       
     *            0,  21,  21,   0,  14,   7,  21,   7,  14,   7,       
     *           14,                                                    
     *            0,  42,   0,   7,   0,  35,  35,   0,   7,  14,       
     *            0,  28,  28,   0,  14,   7,  28,   7,  21,   0,       
     *           21,  14,   7,  21,  21,   7,  14,  14/                 
      DATA JY /   0,   0,  49,   0,   0,  98,   0,  49,   0,  49,       
     *            0, 147,   0,  49,   0,  98,  98,   0,  49,  49,       
     *            0, 196,   0,  49,   0, 147, 147,   0,  49,  98,       
     *            0,  98,  49,  98,  49,                                
     *            0, 245,   0,  49,   0, 196, 196,   0,  49,  98,       
     *            0, 147, 147,   0,  98,  49, 147,  49,  98,  49,       
     *           98,                                                    
     *            0, 294,   0,  49,   0, 245, 245,   0,  49,  98,       
     *            0, 196, 196,   0,  98,  49, 196,  49, 147,   0,       
     *          147,  98,  49, 147, 147,  49,  98,  98/                 
      DATA IY /   1,   1, 344,   1,   1, 687,   1, 344,   1, 344,       
     *            1,1030,   1, 344,   1, 687, 687,   1, 344, 344,       
     *            1,1373,   1, 344,   1,1030,1030,   1, 344, 687,       
     *            1, 687, 344, 687, 344,                                
     *            1,1716,   1, 344,   1,1373,1373,   1, 344, 687,       
     *            1,1030,1030,   1, 687, 344,1030, 344, 687, 344,       
     *          687,                                                    
     *            1,2059,   1, 344,   1,1716,1716,   1, 344, 687,       
     *            1,1373,1373,   1, 687, 344,1373, 344,1030,   1,       
     *         1030, 687, 344,1030,1030, 344, 687, 687/                 
      DATA LZ /   0,   0,   0,   1,   0,   0,   2,   0,   1,   1,       
     *            0,   0,   3,   0,   1,   0,   1,   2,   2,   1,       
     *            0,   0,   4,   0,   1,   0,   1,   3,   3,   0,       
     *            2,   2,   1,   1,   2,                                
     *            0,   0,   5,   0,   1,   0,   1,   4,   4,   0,       
     *            2,   0,   2,   3,   3,   1,   1,   3,   1,   2,       
     *            2,                                                    
     *            0,   0,   6,   0,   1,   0,   1,   5,   5,   0,       
     *            2,   0,   2,   4,   4,   1,   1,   4,   0,   3,       
     *            3,   1,   2,   1,   2,   3,   3,   2/                 
      DATA KZ /   0,   0,   0,   7,   0,   0,  14,   0,   7,   7,       
     *            0,   0,  21,   0,   7,   0,   7,  14,  14,   7,       
     *            0,   0,  28,   0,   7,   0,   7,  21,  21,   0,       
     *           14,  14,   7,   7,  14,                                
     *            0,   0,  35,   0,   7,   0,   7,  28,  28,   0,       
     *           14,   0,  14,  21,  21,   7,   7,  21,   7,  14,       
     *           14,                                                    
     *            0,   0,  42,   0,   7,   0,   7,  35,  35,   0,       
     *           14,   0,  14,  28,  28,   7,   7,  28,   0,  21,       
     *           21,   7,  14,   7,  14,  21,  21,  14/                 
      DATA JZ /   0,   0,   0,  49,   0,   0,  98,   0,  49,  49,       
     *            0,   0, 147,   0,  49,   0,  49,  98,  98,  49,       
     *            0,   0, 196,   0,  49,   0,  49, 147, 147,   0,       
     *           98,  98,  49,  49,  98,                                
     *            0,   0, 245,   0,  49,   0,  49, 196, 196,   0,       
     *           98,   0,  98, 147, 147,  49,  49, 147,  49,  98,       
     *           98,                                                    
     *            0,   0, 294,   0,  49,   0,  49, 245, 245,   0,       
     *           98,   0,  98, 196, 196,  49,  49, 196,   0, 147,       
     *          147,  49,  98,  49,  98, 147, 147,  98/                 
      DATA IZ /   1,   1,   1, 344,   1,   1, 687,   1, 344, 344,       
     *            1,   1,1030,   1, 344,   1, 344, 687, 687, 344,       
     *            1,   1,1373,   1, 344,   1, 344,1030,1030,   1,       
     *          687, 687, 344, 344, 687,                                
     *            1,   1,1716,   1, 344,   1, 344,1373,1373,   1,       
     *          687,   1, 687,1030,1030, 344, 344,1030, 344, 687,       
     *          687,                                                    
     *            1,   1,2059,   1, 344,   1, 344,1716,1716,   1,       
     *          687,   1, 687,1373,1373, 344, 344,1373,   1,1030,       
     *         1030, 344, 687, 344, 687,1030,1030, 687/                 
C                                                                       
C     PREPARE SHELL INFORMATION/FOR HONDO INTEGRATION                   
C                                                                       
      IF(NELEC.EQ.2) GO TO 200                                          
C                                                                       
C     ----- PERMUTE ISH AND JSH SHELLS, FOR THEIR TYPE                  
C     THIS IS DONE FOR SPEED REASONS.  THE CODE GETS THE RIGHT ANSWER   
C     WITHOUT THE ANGULAR MOMENTUM FLIPPING, AND THEREFORE A CALLING    
C     ARGUMENT ALLOWS ONE DO EXACTLY THE INTEGRAL BLOCK AS SPECIFIED,   
C     SHOULD THAT BE DESIRED.                                           
C                                                                       
      IANDJ = ISH .EQ. JSH                                              
      IF (KTYPE(ISH) .LT. KTYPE(JSH)  .AND.  FLIP) THEN                 
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
C                                                                       
C     ----- ISHELL                                                      
C                                                                       
      I = KATOM(INU)                                                    
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
      DO 140 I = I1,I2                                                  
         NGA = NGA+1                                                    
         GA(NGA) = EX(I)                                                
         CSA(NGA) = CS(I)                                               
         CPA(NGA) = CP(I)                                               
         CDA(NGA) = CD(I)                                               
         CFA(NGA) = CF(I)                                               
         CGA(NGA) = CG(I)                                               
         CHA(NGA) = CH(I)                                               
         CIA(NGA) = CI(I)                                               
  140 CONTINUE                                                          
C                                                                       
C     ----- JSHELL                                                      
C                                                                       
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
      DO 160 J = J1,J2                                                  
         NGB = NGB+1                                                    
         GB(NGB) = EX(J)                                                
         CSB(NGB) = CS(J)                                               
         CPB(NGB) = CP(J)                                               
         CDB(NGB) = CD(J)                                               
         CFB(NGB) = CF(J)                                               
         CGB(NGB) = CG(J)                                               
         CHB(NGB) = CH(J)                                               
         CIB(NGB) = CI(J)                                               
  160 CONTINUE                                                          
      RAB = ((AX-BX)*(AX-BX) + (AY-BY)*(AY-BY) + (AZ-BZ)*(AZ-BZ))       
C                                                                       
C     ----- PREPARE INDICES FOR PAIRS OF (I,J) FUNCTIONS                
C                                                                       
      IJ = 0                                                            
      JMAX = MAXJ                                                       
      DO 190 I = MINI,MAXI                                              
         NX = IX(I)                                                     
         NY = IY(I)                                                     
         NZ = IZ(I)                                                     
         IF (IANDJ) JMAX = I                                            
         DO 180 J = MINJ,JMAX                                           
            IJ = IJ+1                                                   
            IJX(IJ) = NX+JX(J)                                          
            IJY(IJ) = NY+JY(J)                                          
            IJZ(IJ) = NZ+JZ(J)                                          
            IJGT(IJ) = NGTI*(I-MINI)+NGTJ*(J-MINJ)+1                    
  180    CONTINUE                                                       
  190 CONTINUE                                                          
      RETURN                                                            
C     ******                                                            
C                                                                       
C        K AND L SHELL                                                  
C                                                                       
  200 CONTINUE                                                          
      KANDL = KSH .EQ. LSH                                              
      SAME = ISH .EQ. KSH .AND. JSH .EQ. LSH                            
C                                                                       
C     ----- PERMUTE KSH AND LSH SHELLS, FOR THEIR TYPE                  
C                                                                       
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
C                                                                       
C     ----- K SHELL                                                     
C                                                                       
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
      DO 260 K = K1,K2                                                  
         NGC = NGC+1                                                    
         GC(NGC) = EX(K)                                                
         CSC(NGC) = CS(K)                                               
         CPC(NGC) = CP(K)                                               
         CDC(NGC) = CD(K)                                               
         CFC(NGC) = CF(K)                                               
         CGC(NGC) = CG(K)                                               
         CHC(NGC) = CH(K)                                               
         CIC(NGC) = CI(K)                                               
  260 CONTINUE                                                          
C                                                                       
C     ----- LSHELL                                                      
C                                                                       
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
      DO 280 L = L1,L2                                                  
         NGD = NGD+1                                                    
         GD(NGD) = EX(L)                                                
         CSD(NGD) = CS(L)                                               
         CPD(NGD) = CP(L)                                               
         CDD(NGD) = CD(L)                                               
         CFD(NGD) = CF(L)                                               
         CGD(NGD) = CG(L)                                               
         CHD(NGD) = CH(L)                                               
         CID(NGD) = CI(L)                                               
  280 CONTINUE                                                          
      NROOTS = (LIT+LJT+LKT+LLT-2)/2                                    
      RCD = ((CX-DX)*(CX-DX) + (CY-DY)*(CY-DY) + (CZ-DZ)*(CZ-DZ))       
C                                                                       
C     ----- PREPARE INDICES FOR PAIRS OF (K,L) FUNCTIONS                
C                                                                       
      KL = 0                                                            
      LMAX = MAXL                                                       
      DO 310 K = MINK,MAXK                                              
         NX = KX(K)                                                     
         NY = KY(K)                                                     
         NZ = KZ(K)                                                     
         IF (KANDL) LMAX = K                                            
         DO 300 L = MINL,LMAX                                           
            KL = KL+1                                                   
            KLX(KL) = NX+LX(L)                                          
            KLY(KL) = NY+LY(L)                                          
            KLZ(KL) = NZ+LZ(L)                                          
            KLGT(KL) = NGTK*(K-MINK)+NGTL*(L-MINL)                      
  300    CONTINUE                                                       
  310 CONTINUE                                                          
      MAX = KL                                                          
      DO 320 I = 1,IJ                                                   
      IF (SAME) MAX = I                                                 
  320 IK(I) = MAX                                                       
      IJKL = IJ*KL                                                      
      IF (SAME) IJKL = IJ*(IJ+1)/2                                      
      RETURN                                                            
      END                                                               
C*MODULE INT2C   *DECK SHELLQUART                                       
      SUBROUTINE SHELLQUART(ISH,JSH,KSH,LSH,GHONDO)                     
      use mx_limits, only: mxgsh,mxg2,mxsh,mxgtot                       
C                                                                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C                                                                       
      DIMENSION GHONDO(*)                                               
      DIMENSION GROTSPD(1296)                                           
C                                                                       
      LOGICAL PACK2E,IANDJ,KANDL,SAME                                   
C                                                                       
C                                                                       
      COMMON /FLIPS / IB(4,3)                                           
C$omp threadprivate(/FLIPS/)
      COMMON /INTFIL/ NINTMX,NHEX,NTUPL,PACK2E,INTTYP,IGRDTYP           
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),      
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),      
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),   
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL           
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK                 
      COMMON /SHLEXC/ NORGSH(3),NORGSP(3),IEXCH,NANGM,NGTH(4)           
C                                                                       
C  FOR RYS QUADRATURE CODE                                              
C                                                                       
      LOGICAL LSHEL                                                     
      COMMON /FMCOM /XX(1)                                              
C                                                                       
C  FOR ERIC CODE, BUT ALL CODES MUST SET ERIOUT FOR OUTPUT ROUTINES     
C                                                                       
      LOGICAL ERICQT                                                    
      COMMON /ERIDAT/ LEN1,LEN2,LEN3,LEN4                               
      COMMON /ERIOUT/ INW,JNW,KNW,LNW,LSTRI,LSTRJ,LSTRK,LSTRL           
C$omp threadprivate(/ERIOUT/)
C                                                                       
C  FOR THE TWO ROTATED AXIS CODES                                       
C                                                                       
      INTEGER IDPOP(4,10)                                               
      LOGICAL SP,SPD,ROTSP,ROTSPD,SPDFG                                 
      COMMON /GOUT  / GPOPLE(768),NORGP                                 
C$omp threadprivate(/GOUT/)
      COMMON /POPOUT/ LPOPI,LPOPJ,LPOPK,LPOPL                           
C$omp threadprivate(/POPOUT/)
      COMMON /SHLG70/ IPL,JPL,KPL,LPL,INEW,JNEW,KNEW,LNEW               
C$omp threadprivate(/SHLG70/)
      DATA IDPOP/0,0,0,0,216,36,6,1,432,72,12,2,648,108,18,3,           
     *           0,0,0,0,216,36,6,1,432,72,12,2,648,108,18,3,           
     *           864,144,24,4,1080,180,30,5/                            
C                                                                       
C     ----- SELECT THE INTEGRAL CODE FOR THIS SHELL QUARTET -----       
C     THE USER INPUT SELECTION -INTTYP- HAS THE FOLLOWING MEANING:      
C      INTTYP=0    BEST TIMING:                                         
C                  USE ROTATED AXIS CODES FOR ANY S,P,D,L SHELLS,       
C                  OTHERWISE PICK ERIC CODE WHENEVER POSSIBLE, BUT      
C                  USE RYS CODE FOR L SHELL OR IF TOO MUCH ANG.MOM.     
C      INTTYP=1    USE THE S,P,L OR S,P,D,L ROTATED AXIS CODE WHENEVER  
C                  POSSIBLE, OTHERWISE RYS QUADRATURE (NO ERIC).        
C      INTTYP=2    USE ERIC CODE AS MUCH AS POSSIBLE, OTHERWISE         
C                  USE THE RYS QUADRATURE (NO ROTATED AXIS).            
C      INTTYP=3    USE RYS POLYNOMIAL QUADRATURE FOR EVERYTHING.        
C                                                                       
      SP    = KTYPE(ISH).LE.2.AND.                                      
     *        KTYPE(JSH).LE.2.AND.                                      
     *        KTYPE(KSH).LE.2.AND.                                      
     *        KTYPE(LSH).LE.2                                           
      SPD   = KTYPE(ISH).LE.3.AND.                                      
     *        KTYPE(JSH).LE.3.AND.                                      
     *        KTYPE(KSH).LE.3.AND.                                      
     *        KTYPE(LSH).LE.3                                           
      IF(SP) SPD=.FALSE.                                                
      LSHEL = (KMAX(ISH)-KMIN(ISH)+1).EQ.4.OR.                          
     *        (KMAX(JSH)-KMIN(JSH)+1).EQ.4.OR.                          
     *        (KMAX(KSH)-KMIN(KSH)+1).EQ.4.OR.                          
     *        (KMAX(LSH)-KMIN(LSH)+1).EQ.4                              
      LQSUM = KTYPE(ISH) + KTYPE(JSH) + KTYPE(KSH) + KTYPE(LSH) - 4     
C                                                                       
      ROTSP  = SP                                                       
      ROTSPD = SPD                                                      
      SPDFG = KTYPE(ISH).LE.5.AND.                                      
     *        KTYPE(JSH).LE.5.AND.                                      
     *        KTYPE(KSH).LE.5.AND.                                      
     *        KTYPE(LSH).LE.5                                           
      ERICQT = .NOT.LSHEL .AND. LQSUM.LE.5 .AND. SPDFG                  
C                                                                       
C        RYS QUADRATURE IS A BIT FASTER AT UNCONTRACTED QUARTETS        
C                                                                       
      KQCON = KNG(ISH) * KNG(JSH) * KNG(KSH) * KNG(LSH) + INTTYP        
      IF(KQCON.EQ.1  .AND.  .NOT.SP) THEN                               
         ROTSP  = .FALSE.                                               
         ROTSPD = .FALSE.                                               
         ERICQT = .FALSE.                                               
      END IF                                                            
C                                                                       
C        INPUT OVERRIDES                                                
C                                                                       
      IF(INTTYP.EQ.1) ERICQT = .FALSE.                                  
      IF(INTTYP.GE.2) ROTSP  = .FALSE.                                  
      IF(INTTYP.GE.2) ROTSPD = .FALSE.                                  
      IF(INTTYP.EQ.3) ERICQT = .FALSE.                                  
C                                                                       
C        THE VARIOUS PACKAGES ARE CALLED BELOW IN THE ORDER OF          
C        FIRST ROTATED AXIS, THEN ERIC, FINALLY RYS QUADRATURE.         
C        SINCE EACH CODE RETURNS, A QUARTET IS NEVER DONE TWICE.        
C                                                                       
C  ROTATED AXIS CODE FOR PURE SP SHELL QUARTET                          
C                                                                       
      IF (ROTSP) THEN                                                   
        IPL = ISH                                                       
        JPL = JSH                                                       
        KPL = KSH                                                       
        LPL = LSH                                                       
        INW = ISH                                                       
        JNW = JSH                                                       
        KNW = KSH                                                       
        LNW = LSH                                                       
        NORGP = NORGSP(IEXCH)                                           
        NORGH = NORGSH(IEXCH)                                           
C                                                                       
        CALL GENR70(1,.FALSE.)                                          
C                                                                       
C  SAVE TO OUTPUT ARRAY WITH HONDO INDEXING                             
C                                                                       
        MINI = KMIN(INW)                                                
        MAXI = KMAX(INW)                                                
        MINJ = KMIN(JNW)                                                
        MAXJ = KMAX(JNW)                                                
        MINK = KMIN(KNW)                                                
        MAXK = KMAX(KNW)                                                
        MINL = KMIN(LNW)                                                
        MAXL = KMAX(LNW)                                                
C                                                                       
        II = 1                                                          
        DO I = MINI, MAXI                                               
          IP = (I-1)*LPOPI + 1                                          
          IJ  = II                                                      
          DO J = MINJ, MAXJ                                             
            IJP = (J-1)*LPOPJ + IP                                      
            IJK  = IJ                                                   
            DO K = MINK, MAXK                                           
              IJKP = (K-1)*LPOPK + IJP                                  
              IJKL  = IJK                                               
              DO L = MINL, MAXL                                         
                IJKLP = (L-1)*LPOPL + IJKP                              
                GHONDO(IJKL+NORGH) = GPOPLE(IJKLP+NORGP)                
                IJKL = IJKL  + LEN1                                     
              END DO                                                    
              IJK  = IJK  + LEN2                                        
            END DO                                                      
            IJ  = IJ  + LEN3                                            
          END DO                                                        
          II = II + LEN4                                                
        END DO                                                          
        LSTRI = LEN4                                                    
        LSTRJ = LEN3                                                    
        LSTRK = LEN2                                                    
        LSTRL = LEN1                                                    
        RETURN                                                          
C       ******                                                          
C                                                                       
C  ROTATED AXIS CODE FOR QUARTET CONTAINING AT LEAST ONE D FUNCTION     
C                                                                       
      ELSE IF (ROTSPD) THEN                                             
        INW = ISH                                                       
        JNW = JSH                                                       
        KNW = KSH                                                       
        LNW = LSH                                                       
C                                                                       
        CALL GENR03(GROTSPD)                                            
C                                                                       
C  SAVE TO OUTPUT ARRAY WITH HONDO INDEXING                             
C                                                                       
        NORGH = NORGSH(IEXCH)                                           
        IANDJ = ISH.EQ.JSH                                              
        KANDL = KSH.EQ.LSH                                              
        SAME  = ISH.EQ.KSH  .AND.  JSH.EQ.LSH                           
        IF(NOPK.EQ.0) SAME=.FALSE.                                      
C                                                                       
        IEX=1                                                           
        IBB = IB(1,IEX)                                                 
        JBB = IB(2,IEX)                                                 
        KBB = IB(3,IEX)                                                 
        LBB = IB(4,IEX)                                                 
C                                                                       
        MINI = KMIN(INW)                                                
        MAXI = KMAX(INW)                                                
        MINJ = KMIN(JNW)                                                
        MAXJ = KMAX(JNW)                                                
        MINK = KMIN(KNW)                                                
        MAXK = KMAX(KNW)                                                
        MINL = KMIN(LNW)                                                
        MAXL = KMAX(LNW)                                                
C                                                                       
        IJN = 0                                                         
        JMAX = MAXJ                                                     
        DO I = MINI, MAXI                                               
          IHONDO = (I-MINI)*LEN4 + 1                                    
          IROTAX = IDPOP(IBB,I)  + 1                                    
          IF(IANDJ) JMAX=I                                              
          DO 340 J = MINJ, JMAX                                         
            IJHONDO = (J-MINJ)*LEN3 + IHONDO                            
            IJROTAX = IDPOP(JBB,J)  + IROTAX                            
            IJN = IJN+1                                                 
            LMAX=MAXL                                                   
            KLN=0                                                       
            DO K = MINK, MAXK                                           
              IJKHONDO = (K-MINK)*LEN2 + IJHONDO                        
              IJKROTAX = IDPOP(KBB,K)  + IJROTAX                        
              IF(KANDL) LMAX=K                                          
              DO L = MINL, LMAX                                         
                KLN = KLN+1                                             
                IF(SAME .AND. KLN.GT.IJN) GO TO 340                     
                IJKLHONDO = (L-MINL)*LEN1 + IJKHONDO                    
                IJKLROTAX = IDPOP(LBB,L)  + IJKROTAX                    
                GHONDO(IJKLHONDO+NORGH) = GROTSPD(IJKLROTAX)            
              END DO                                                    
            END DO                                                      
  340     CONTINUE                                                      
        END DO                                                          
        LSTRI = LEN4                                                    
        LSTRJ = LEN3                                                    
        LSTRK = LEN2                                                    
        LSTRL = LEN1                                                    
        RETURN                                                          
C       ******                                                          
C                                                                       
C  USE ERIC FAST CODES, REQUIRES THAT LQSYM.LE.5 AND NO L-SHELLS        
C  NOTE THAT THERE IS CODE COPYING ERIC BUFFERS INTO THE HONDO          
C  FORMAT BUFFER AT THE END OF THE -ERIC- ROUTINE.                      
C                                                                       
      ELSE IF (ERICQT) THEN                                             
        NORGH = NORGSH(IEXCH)                                           
        CALL ERIC(ISH,JSH,KSH,LSH,GHONDO(1+NORGH))                      
        RETURN                                                          
C       ******                                                          
C                                                                       
C  GENERAL CASE = HONDO/RYS QUADRATURE: ANY S,P,D,F,G,H,I OR L SHELLS   
C                                                                       
      ELSE                                                              
        CALL VALFM(LOADFM)                                              
        IDDIJ = LOADFM + 1                                              
        NEED  = 49*MXG2                                                 
        CALL GETFM(NEED)                                                
        CALL SHELLS(1,ISH,JSH,KSH,LSH,.TRUE.)                           
        CALL SHELLS(2,ISH,JSH,KSH,LSH,.TRUE.)                           
        CALL IJPRIM(XX(IDDIJ))                                          
        NORGH = NORGSH(IEXCH)                                           
        IF(NOPK.NE.0) CALL ZQOUT(GHONDO)                                
        IF(LQSUM.EQ.0) THEN                                             
           CALL S0000(GHONDO(1+NORGH),XX(IDDIJ))                        
        ELSE                                                            
           CALL GENRAL(GHONDO(1+NORGH),XX(IDDIJ))                       
        END IF                                                          
        CALL RETFM(NEED)                                                
        RETURN                                                          
C       ******                                                          
C                                                                       
      END IF                                                            
      END                                                               
C*MODULE INT2A   *DECK TWOEI                                            
      SUBROUTINE TWOEI(TYPSCF,DIRSCF,DIRNLO,DIRTRF,DIRCIS,              
     *                 INTTYP,SCHWRZ,NINT,NSCHWZ,L1,L2,                 
     *                 BUFP,BUFK,IX,NINTMX,                             
     *                 XINTS,NSH2,GHONDO,MAXG,DDIJ,                     
     *                 IA,DA,FA,DB,FB,DSH,DNLO,FNLO,NFLMAT)             
      use mx_limits, only: mxsh,mxgtot,mxatm                            
C                                                                       
C       MANY ARGUMENTS ARE OPTIONAL, YOU MUST ALLOCATE STORAGE FOR      
C              ALL CALLS: GHONDO, DDIJ, XINTS                           
C           CONVENTIONAL: BUFP, IX, AND POSSIBLY BUFK                   
C             DIRECT SCF: IA, DSH, DA, FA, AND POSSIBLY DB, FB          
C             DIRECT CIS: IA, DSH, DA, FA, AND POSSIBLY DB, FB          
C             DIRECT NLO: DNLO, FNLO, NFLMAT                            
C  DIRECT TRANSFORMATION: BUFP, IX                                      
C     RESPONSE EQUATIONS: MUST DEFINE NFLMAT.NE.1                       
C  NOTE THAT THE TYPE OF FOCK MATRIX BUILT (TYPSCF) DOES NOT            
C  NECESSARILY HAVE TO MATCH THIS RUNS SCFTYP VALUE.                    
C                                                                       
C$    USE params, ONLY: intomp, shfock                                  
C$    USE ompmod                                                        
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      LOGICAL DIRSCF,DIRNLO,DIRTRF,DIRCIS,CMBDIR                        
      LOGICAL OUT,SCHWRZ,SCHSKP,GOPARR,DSKWRK,MASWRK,DLB,SLB,C1GRP      
      LOGICAL SKIPA,SKIPB,SKIPC,NPSYM                                   
      LOGICAL PK,PANDK,NOTPK,BLOCK,GPSAVE,SCREEN                        
      LOGICAL LTRMST                                                    
C                                                                       
      CHARACTER*8 INAOFL                                                
C                                                                       
      DIMENSION BUFP(NINTMX),BUFK(NINTMX),IX(*),XINTS(NSH2),            
     *          GHONDO(MAXG),IA(L1),DA(L2),FA(L2),DB(L2),FB(L2),        
     *          DSH(NSH2),DDIJ(*),DNLO(L1,L1),FNLO(L1,L1)               
      DIMENSION MI(48),MJ(48),MK(48),M0(48)                             
C                                                                       
C                                                                       
      COMMON /ELGFIL/ NFILE,INAOFL                                      
      COMMON /ELGTRM/ LTRMST,NFLTRM,NRCTRM,NPSTRM,NHTSHL                
      COMMON /INT2IC/ NINTIC,ININTIC,NXXIC,LBUFPIC,LIXIC,LABSIX,NINTIX  
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)                
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),      
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),      
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),   
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL           
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK                 
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
      COMMON /PKFIL / PK,PANDK,BLOCK                                    
      COMMON /RESTAR/ TIMLIM,IREST,NREC,INTLOC,IST,JST,KST,LST          
      COMMON /SCINP / VLAMB,SCREEN                                      
      COMMON /SHLEXC/ NORGSH(3),NORGSP(3),IEXCH,NANGM,NGTH(4)           
      COMMON /SHLG70/ ISH,JSH,KSH,LSH,IJKLXX(4)                         
C$omp threadprivate(/SHLG70/)
      COMMON /SHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,          
     *                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,          
     *                NIJ,IJ,KL,IJKL                                    
C$omp threadprivate(/SHLNOS/)
      COMMON /SHLT  / TOL,CUTOFF,ICOUNT,OUT                             
      COMMON /SYMTRY/ MAPSHL(MXSH,48),MAPCTR(MXATM,48),                 
     *                T(432),INVT(48),NT                                
C                                                                       
      PARAMETER (ZERO=0.0D+00)                                          
C                                                                       
C$    DOUBLE PRECISION, ALLOCATABLE :: fb_cphf(:)                       
C$    DATA DHRHF,DHUHF,DHROHF /8HRHF     ,8HUHF     ,8HROHF     /       
C$    DATA DHGVB /8HGVB     /                                           
C                                                                       
C          ----- TWO-ELECTRON INTEGRALS -----                           
C     THIS VERSION CAN HANDLE S,P,D,F,G,H,I AND L SHELLS                
C                                                                       
      TIM = ZERO                                                        
      CALL TSECND(TIM)                                                  
C                                                                       
C           PERHAPS CALL XABI LOPEZ/JOSE UGALDE'S SCREENED INTEGRAL CODE
C                                                                       
      IF (SCREEN) THEN                                                  
         WRITE(IW,36)                                                   
  36     FORMAT (/,'#########',/,'SCREEN TWO-ELECTRON INTEGRAL',/)      
         CALL STWOEI(TYPSCF,DIRSCF,DIRNLO,DIRTRF,                       
     *               INTTYP,SCHWRZ,NINT,NSCHWZ,L1,L2,                   
     *               BUFP,BUFK,IX,NINTMX,                               
     *               XINTS,NSH2,GHONDO,MAXG,DDIJ,                       
     *               IA,DA,FA,DB,FB,DSH,DNLO,FNLO)                      
         RETURN                                                         
      END IF                                                            
C                                                                       
      CMBDIR= DIRSCF .OR. DIRNLO .OR. DIRTRF .OR. DIRCIS                
C                                                                       
C        THE OLD FASHIONED PARALLEL INTEGRAL TRANSFORMATIONS DO NOT     
C        ALLOW THE AO INTEGRAL WORK TO BE RUN IN PARALLEL.  THE MODERN  
C        DISTRIBUTED MEMORY TRANSFORMATIONS DO NOT PASS THROUGH HERE.   
C                                                                       
      ICONT=  0                                                         
      GPSAVE = GOPARR                                                   
      IF(DIRTRF) GOPARR=.FALSE.                                         
C                                                                       
C     ----- INITIALIZATION FOR PARALLEL WORK -----                      
C     BOTH STATIC AND DYNAMIC LOAD BALANCING ARE IMPLEMENTED BELOW      
C                                                                       
      SLB = GOPARR  .AND.  IBTYP.EQ.0                                   
      DLB = GOPARR  .AND.  IBTYP.EQ.1                                   
      NEXT = -1                                                         
      MINE = -1                                                         
      IPCOUNT = ME - 1                                                  
C                                                                       
      C1GRP = NT.EQ.1                                                   
C                                                                       
      CALL BASCHK(LMAX)                                                 
                    NANGM =  4                                          
      IF(LMAX.EQ.2) NANGM =  6                                          
      IF(LMAX.EQ.3) NANGM = 10                                          
      IF(LMAX.EQ.4) NANGM = 15                                          
      IF(LMAX.EQ.5) NANGM = 21                                          
      IF(LMAX.EQ.6) NANGM = 28                                          
      NGTH(4) = 1                                                       
      NGTH(3) = NGTH(4) * NANGM                                         
      NGTH(2) = NGTH(3) * NANGM                                         
      NGTH(1) = NGTH(2) * NANGM                                         
      IF(NOPK.EQ.0) THEN                                                
         NORGSH(1) = 0                                                  
         NORGSH(2) = NORGSH(1) + NANGM**4                               
         NORGSH(3) = NORGSH(2) + NANGM**4                               
         NORGSP(1) = 0                                                  
         NORGSP(2) = 256                                                
         NORGSP(3) = 512                                                
      ELSE                                                              
         DO I=1,3                                                       
            NORGSH(I) = 0                                               
            NORGSP(I) = 0                                               
         ENDDO                                                          
      END IF                                                            
C                                                                       
C     CALL THREADED VERSION IF POSSIBLE                                 
C                                                                       
C       The code below is only active                                   
C       with openmp enabled by compiler                                 
C                                                                       
!     Temporary workaround for CPHF case                                
C$    IF (intomp.GT.0.AND.typscf.EQ.dhrhf) ALLOCATE(fb_cphf(nflmat))    
C$                                                                      
C$!    IF (.not.(DIRNLO.OR.DIRTRF).and.DIRSCF.AND.(NFLMAT.LE.1)) THEN   
C$    IF (.not.(DIRNLO.OR.DIRTRF).and.DIRSCF.and.(intomp.ne.0)) THEN    
C$      IF (intomp.EQ.1) THEN                                           
C$          IF (typscf.EQ.dhrhf) THEN                                   
C$              CALL ompmod_twoei_jk(typscf,schwrz,nint,nschwz,         
C$   +                               l1,l2,1,xints,                     
C$   +                               nsh2,maxg,                         
C$   +                               ia,da,fa,db,fb_cphf,dsh,nflmat,    
C$   +                               cutoff,out)                        
C$          ELSE                                                        
C$              CALL ompmod_twoei_jk(typscf,schwrz,nint,nschwz,         
C$   +                               l1,l2,l2,xints,                    
C$   +                               nsh2,maxg,                         
C$   +                               ia,da,fa,db,fb,dsh,nflmat,         
C$   +                               cutoff,out)                        
C$          END IF                                                      
C$                                                                      
C$      ELSEIF (intomp.EQ.2) THEN                                       
C$                                                                      
C$          IF (shfock.AND.(typscf.EQ.dhrhf).and.(nflmat.le.1)) THEN    
C$              CALL ompmod_twoei_shf_kl_rhf(typscf,schwrz,nint,nschwz, 
C$   +                                       l1,l2,1,xints,             
C$   +                                       nsh2,maxg,                 
C$   +                                       ia,da,fa,dsh,nflmat,       
C$   +                                       cutoff,out)                
C$          ELSE                                                        
C$              IF (typscf.EQ.dhrhf) THEN                               
C$                  CALL ompmod_twoei_kl(typscf,schwrz,nint,nschwz,     
C$   +                                   l1,l2,1,xints,                 
C$   +                                   nsh2,maxg,                     
C$   +                                   ia,da,fa,db,fb_cphf,dsh,nflmat,
C$   +                                   cutoff,out)                    
C$              ELSE                                                    
C$                  CALL ompmod_twoei_kl(typscf,schwrz,nint,nschwz,     
C$   +                                   l1,l2,l2,xints,                
C$   +                                   nsh2,maxg,                     
C$   +                                   ia,da,fa,db,fb,dsh,nflmat,     
C$   +                                   cutoff,out)                    
C$              END IF                                                  
C$          END IF                                                      
C$      ELSE                                                            
C$!       Should never get here, just in case                           
C$        IF(MASWRK) WRITE(IW,'(" INCONSISTENT OPENMP PARAMETERS")')    
C$        IF(MASWRK) WRITE(IW,'(" ABORTING CALCULATION")')              
C$        CALL ABRT                                                     
C$      END IF                                                          
C$      RETURN                                                          
C$    END IF                                                            
C                                                                       
C     ----- END OPENMP CALCULATION -----                                
C                                                                       
      NOTPK = .NOT.PK                                                   
      NINT  = 0                                                         
      NSCHWZ= 0                                                         
      SCHSKP=.FALSE.                                                    
      DENMAX = ZERO                                                     
C                                                                       
C        NOW WE ARE READY TO LOOP OVER ALL NSHELL**4 SHELL QUARTETS     
C                                                                       
C     ----- I SHELL -----                                               
C                                                                       
      DO 920 II = IST,NSHELL                                            
C                                                                       
C     ----- CHECK CPU TIME -----                                        
C                                                                       
      CALL TSECND(TIM)                                                  
      IF(TIM.GE.TIMLIM) THEN                                            
C        NOTHING CAN BE DONE FOR IN-CORE INTEGRALS: JUST FORGET THEM    
         IF(.NOT.CMBDIR)                                                
     *      CALL FINAL(0,II,1,1,1,PANDK,BUFP,BUFK,IX,NINTMX)            
         IF(MASWRK) WRITE(IW,9030)                                      
         CALL ABRT                                                      
      END IF                                                            
C                                                                       
C         ELONGATION METHOD                                             
      IF(LTRMST) THEN                                                   
         IF(II.EQ.NHTSHL) THEN                                          
            NFLTRM = NFILE                                              
            NRCTRM = NREC                                               
            NPSTRM = ICOUNT                                             
         ENDIF                                                          
      ENDIF                                                             
C                                                                       
C     ----- PRINT INTERMEDIATE RESTART DATA -----                       
C                                                                       
      IF(NPRINT.NE.-5  .AND.  .NOT.CMBDIR .AND. MASWRK) THEN            
         IF(ICOUNT.LE.NINTIC) THEN                                      
            WRITE(IW,9015) II,JST,KST,LST,ICOUNT                        
         ELSE                                                           
            WRITE(IW,9010) II,JST,KST,LST,NREC,ICOUNT-NINTIC            
         ENDIF                                                          
      ENDIF                                                             
C                                                                       
C     ----- SKIP I SHELL IF NOT SYMMETRY UNIQUE -----                   
C     THIS, AND THE SIMILAR BRANCHINGS FOR THE J, K, AND L LOOPS IS     
C     WHAT GENERATES THE "PETITE" RATHER THAN "GRANDE" INTEGRAL LIST.   
C                                                                       
      IF(C1GRP) THEN                                                    
         MI(1)=II                                                       
      ELSE                                                              
         DO 120 IT = 1,NT                                               
            ID = MAPSHL(II,IT)                                          
            IF (ID .GT. II) GO TO 920                                   
            MI(IT) = ID                                                 
  120    CONTINUE                                                       
      END IF                                                            
C                                                                       
C     ----- J SHELL -----                                               
C                                                                       
      J0 = JST                                                          
      DO 900 JJ = J0,II                                                 
      JST = 1                                                           
C                                                                       
      IF(C1GRP) THEN                                                    
         MJ(1)=JJ                                                       
      ELSE                                                              
         DO 200 IT = 1,NT                                               
            ID = MI(IT)                                                 
            JD = MAPSHL(JJ,IT)                                          
            MJ(IT) = JD                                                 
            IF (ID .GE. JD) GO TO 160                                   
            ND = ID                                                     
            ID = JD                                                     
            JD = ND                                                     
  160       IF (ID-II) 200,180,900                                      
  180       IF (JD-JJ) 200,200,900                                      
  200    CONTINUE                                                       
      END IF                                                            
C                                                                       
C     ----- GO PARALLEL! -----                                          
C                                                                       
      IF (DLB) THEN                                                     
         MINE = MINE + 1                                                
         IF (MINE.GT.NEXT) CALL DDI_DLBNEXT(NEXT)                       
         IF (NEXT.NE.MINE) GO TO 900                                    
      END IF                                                            
C                                                                       
C     ----- K SHELL -----                                               
C                                                                       
      K0 = KST                                                          
      DO 880 KK = K0,JJ                                                 
      KST = 1                                                           
C                                                                       
      IF(C1GRP) THEN                                                    
         MK(1)=KK                                                       
      ELSE                                                              
         DO 340 IT = 1,NT                                               
            ID = MI(IT)                                                 
            JD = MJ(IT)                                                 
            KD = MAPSHL(KK,IT)                                          
            MK(IT) = KD                                                 
  240       IF (ID .GE. JD) GO TO 260                                   
            ND = ID                                                     
            ID = JD                                                     
            JD = ND                                                     
  260       IF (JD .GE. KD) GO TO 280                                   
            ND = JD                                                     
            JD = KD                                                     
            KD = ND                                                     
            GO TO 240                                                   
  280       IF (ID-II) 340,300,880                                      
  300       IF (JD-JJ) 340,320,880                                      
  320       IF (KD-KK) 340,340,880                                      
  340    CONTINUE                                                       
      END IF                                                            
C                                                                       
C     ----- GO PARALLEL! -----                                          
C                                                                       
      IF(SLB) THEN                                                      
         IPCOUNT = IPCOUNT + 1                                          
         IF (MOD(IPCOUNT,NPROC).NE.0) GO TO 880                         
      END IF                                                            
C                                                                       
C     ----- L SHELL ----                                                
C                                                                       
      L0 = LST                                                          
      DO 860 LL = L0,KK                                                 
      LST = 1                                                           
C                                                                       
      IF(C1GRP) THEN                                                    
         M0(1)=1                                                        
         N4=1                                                           
      ELSE                                                              
         N4 = 0                                                         
         DO 540 IT = 1,NT                                               
            ID = MI(IT)                                                 
            JD = MJ(IT)                                                 
            KD = MK(IT)                                                 
            LD = MAPSHL(LL,IT)                                          
  380       IF (ID .GE. JD) GO TO 400                                   
            ND = ID                                                     
            ID = JD                                                     
            JD = ND                                                     
  400       IF (JD .GE. KD) GO TO 420                                   
            ND = JD                                                     
            JD = KD                                                     
            KD = ND                                                     
            GO TO 380                                                   
  420       IF (KD .GE. LD) GO TO 440                                   
            ND = KD                                                     
            KD = LD                                                     
            LD = ND                                                     
            GO TO 400                                                   
  440       IF (ID-II) 540,460,860                                      
  460       IF (JD-JJ) 540,480,860                                      
  480       IF (KD-KK) 540,500,860                                      
  500       IF (LD-LL) 540,520,860                                      
  520       N4 = N4+1                                                   
            M0(N4) = IT                                                 
  540    CONTINUE                                                       
      END IF                                                            
C                                                                       
C         THE LOOP STRUCTURE IN THIS ROUTINE IS DESIGNED TO FACILITATE  
C         SUPERMATRIX CONSTRUCTION BY HAVING UP TO THREE "EXCHANGED"    
C         QUARTETS AVAILABLE AT ONCE.  THE LOOP STRUCTURE TO GENERATE   
C         A MORE NORMAL CANONICAL ORDERING OF THE QUARTETS HITS THE     
C         SAME QUARTETS IN A SLIGHTLY DIFFERENT ORDER, BUT BOTH LOOPS   
C         WILL DO EXACTLY THE SAME QUARTETS.                            
C                                                                       
C             CANONICAL                      SUPERMATRIX                
C         DO ISH=1,NSHELL                 DO II=1,NSHELL                
C           DO JSH=1,ISH                    DO JJ=1,II                  
C             IJSH = IA(ISH)+JSH                                        
C             DO KSH=1,ISH                    DO KK=1,JJ                
C               DO LSH=1,KSH                    DO LL=1,KK              
C                 KLSH=IA(KSH)+LSH                                      
C                 IF(IJSH.LT.KLSH),               [II JJ|KK LL],        
C                    CYCLE KSH LOOP               [II KK|JJ LL],        
C                 [ISH JSH|KSH LSH]               [II LL|JJ KK]         
C               ENDDO                           ENDDO                   
C             ENDDO                           ENDDO                     
C           ENDDO                           ENDDO                       
C         ENDDO                           ENDDO                         
C                                                                       
C     ----- CHECK FOR REDUNDANIES BETWEEN THE 3 COMBINATIONS -----      
C            (II,JJ//KK,LL), (II,KK//JJ,LL), (II,LL//JJ,KK)             
C                                                                       
      SKIPA =  JJ.EQ.KK                                                 
      SKIPB = (II.EQ.KK) .OR. (JJ.EQ.LL)                                
      SKIPC = (II.EQ.JJ) .OR. (KK.EQ.LL)                                
      NPSYM = .FALSE.                                                   
      IF (SKIPA .OR. SKIPB .OR. SKIPC) GO TO 720                        
      NPSYM = .TRUE.                                                    
      DO 640 M = 1,N4                                                   
         IT = M0(M)                                                     
         IH = MI(IT)                                                    
         JH = MJ(IT)                                                    
         IF(JH.LE.IH) THEN                                              
            ID = IH                                                     
            JD = JH                                                     
         ELSE                                                           
            ID = JH                                                     
            JD = IH                                                     
         END IF                                                         
         IF(.NOT.SKIPA) SKIPA = (ID.EQ.II .AND. JD.EQ.KK) .OR.          
     *                          (ID.EQ.JJ .AND. JD.EQ.LL)               
         IF(.NOT.SKIPB) SKIPB = (ID.EQ.II .AND. JD.EQ.LL) .OR.          
     *                          (ID.EQ.JJ .AND. JD.EQ.KK)               
         IF (SKIPA .AND. SKIPB) GO TO 660                               
         KH = MK(IT)                                                    
         IF(KH.LE.IH) THEN                                              
            ID = IH                                                     
            KD = KH                                                     
         ELSE                                                           
            ID = KH                                                     
            KD = IH                                                     
         END IF                                                         
         IF(.NOT.SKIPC) SKIPC = (ID.EQ.II .AND. KD.EQ.LL) .OR.          
     *                          (ID.EQ.JJ .AND. KD.EQ.KK)               
         IF(SKIPA .AND. SKIPC) GO TO 680                                
         IF(SKIPB .AND. SKIPC) GO TO 700                                
  640 CONTINUE                                                          
      GO TO 720                                                         
C                                                                       
  660 SKIPC = .TRUE.                                                    
      GO TO 720                                                         
  680 SKIPB = .TRUE.                                                    
      GO TO 720                                                         
  700 SKIPA = .TRUE.                                                    
C                                                                       
C        GENERATE SYMMETRY FACTOR -Q4- FOR THIS QUARTET IN PETITE LIST  
C                                                                       
  720 CONTINUE                                                          
      Q4 = NT                                                           
      Q4 = Q4 / N4                                                      
C                                                                       
C     ----- (II,JJ//KK,LL) -----                                        
C                                                                       
      IEXCH = 1                                                         
      ISH = II                                                          
      JSH = JJ                                                          
      KSH = KK                                                          
      LSH = LL                                                          
      QQ4 = Q4                                                          
      IF(SKIPA .AND. NPSYM) QQ4 = QQ4+Q4                                
      IF(SKIPB .AND. NPSYM) QQ4 = QQ4+Q4                                
      GO TO 780                                                         
C                                                                       
C     ----- (II,KK//JJ,LL) -----                                        
C                                                                       
  740 IF (SKIPA) GO TO 760                                              
      IEXCH = 2                                                         
      ISH = II                                                          
      JSH = KK                                                          
      KSH = JJ                                                          
      LSH = LL                                                          
      QQ4 = Q4                                                          
      IF (SKIPC .AND. NPSYM) QQ4 = QQ4+Q4                               
      GO TO 780                                                         
C                                                                       
C     ----- (II,LL//JJ,KK) -----                                        
C                                                                       
  760 IF (SKIPB .OR. SKIPC) GO TO 840                                   
      IEXCH = 3                                                         
      ISH = II                                                          
      JSH = LL                                                          
      KSH = JJ                                                          
      LSH = KK                                                          
      QQ4 = Q4                                                          
C                                                                       
C        ----- COMPUTE TWO-ELECTRON INTEGRALS ----                      
C                                                                       
  780 CONTINUE                                                          
      IF(PK .AND. IEXCH.EQ.1) CALL ZPKOUT(ISH,JSH,KSH,LSH,GHONDO,       
     *                                    SKIPA,SKIPB,SKIPC,NPSYM)      
C                                                                       
C     APPLY THE SCHWARZ INEQUALITY TO SCREEN OUT SMALL INTEGRALS,       
C      (II,JJ//KK,LL) .LE.  SQRT( (II,JJ//II,JJ)*(KK,LL//KK,LL) )       
C     SEE, FOR EXAMPLE, J.L.WHITTEN, J.CHEM.PHYS. 58,4496-4501(1973)    
C                                                                       
      IF(SCHWRZ) THEN                                                   
         IJIJ = (ISH*ISH-ISH)/2 + JSH                                   
         KLKL = (KSH*KSH-KSH)/2 + LSH                                   
         TEST = QQ4*XINTS(IJIJ)*XINTS(KLKL)                             
         IF(DIRSCF) THEN                                                
            DENMAX = SCHWDN(DSH,ISH,JSH,KSH,LSH,IA)                     
            TEST = TEST*DENMAX                                          
         END IF                                                         
         SCHSKP = TEST.LT.CUTOFF                                        
         IF(SCHSKP) NSCHWZ = NSCHWZ + 1                                 
      END IF                                                            
      IF(SCHSKP) GO TO 820                                              
C                                                                       
C        ----- ELECTRON REPULSION INTEGRAL CALCULATION -----            
C     THIS MAY USE ROTATED AXIS, ERIC, OR RYS QUADRATURE METHODS        
C                                                                       
      CALL SHELLQUART(ISH,JSH,KSH,LSH,GHONDO)                           
C                                                                       
C        USE THE INTEGRALS JUST FORMED.  AT MOST, 1 OF THESE IS CALLED  
C                                                                       
           IF(DIRSCF) THEN                                              
         CALL DIRFCK(TYPSCF,IA,DA,FA,DB,FB,GHONDO,L2,NINT,NFLMAT)       
      ELSE IF(DIRCIS) THEN                                              
         CALL DRFCIS(DA,FA,DB,FB,GHONDO,L1,NINT,NFLMAT)                 
      ELSE IF(DIRNLO) THEN                                              
         CALL DFCKNS(DNLO,FNLO,GHONDO,L1,NINT,NFLMAT)                   
      ELSE IF(DIRTRF) THEN                                              
         CALL DIRTRN(BUFP,IX,NINTMX,GHONDO,ICONT,NINT)                  
      ELSE                                                              
         IF(NOTPK) CALL QOUT(BUFP,IX,NINTMX,GHONDO)                     
      END IF                                                            
C                                                                       
  820 CONTINUE                                                          
      GO TO (740,760,840),IEXCH                                         
C                                                                       
C     ----- WRITE THE P (OR PK) SUPERMATRIX TO DISK FILE -----          
C                                                                       
  840 CONTINUE                                                          
      IF(PK) CALL PKFILE(II,JJ,KK,LL,SKIPA,SKIPB,SKIPC,NPSYM,           
     *                   BUFP,BUFK,IX,NINTMX,GHONDO)                    
C                                                                       
C     ----- END OF SHELL LOOPS -----                                    
C                                                                       
  860 CONTINUE                                                          
  880 CONTINUE                                                          
  900 CONTINUE                                                          
  920 CONTINUE                                                          
C                                                                       
      IF(DLB) CALL DDI_DLBRESET                                         
C                                                                       
C     ----- OUTPUT THE LAST BITS OF INTEGRALS -----                     
C                                                                       
      IF(DIRTRF) CALL ONEIDX(BUFP,IX,ICONT)                             
      IF(.NOT.CMBDIR) THEN                                              
         IF(SCHWRZ) THEN                                                
            IF(GOPARR) CALL DDI_GSUMI(1055,NSCHWZ,1)                    
            IF(NPRINT.NE.-5 .AND. MASWRK) WRITE(IW,9020) NSCHWZ         
         END IF                                                         
         CALL FINAL(1,II,II,II,II,PANDK,BUFP,BUFK,IX,NINTMX)            
      END IF                                                            
C                                                                       
      GOPARR = GPSAVE                                                   
      RETURN                                                            
C                                                                       
 9010 FORMAT(1X,'II,JST,KST,LST =',4I3,' NREC =',I10,' INTLOC =',I5)    
 9015 FORMAT(1X,'II,JST,KST,LST =',4I3,' IN CORE, INTLOC =',I12)        
 9020 FORMAT(1X,'SCHWARZ INEQUALITY TEST SKIPPED',I12,                  
     *        ' INTEGRAL BLOCKS.')                                      
 9030 FORMAT(//1X,'*** THIS JOB HAS EXHAUSTED ITS CPU TIME ***'/        
     *         1X,'     (WHILE COMPUTING 2E- INTEGRALS)'///)            
      END                                                               
C*MODULE INT2A   *DECK XYZINT                                           
      SUBROUTINE XYZINT                                                 
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      LOGICAL N0,N1,M0,M1,FIRST1,FIRST2,FIRST3,FIRST4                   
C                                                                       
      COMMON /SETINT/ I(13),K(13),NIMAX,NJMAX,NKMAX,NLMAX,NMAX,MMAX     
     +               ,BP01,B00,B10,XCP00,XC00,YCP00,YC00,ZCP00,ZC00,F00 
     +               ,DXIJ,DYIJ,DZIJ,DXKL,DYKL,DZKL                     
C$omp threadprivate(/SETINT/)
      COMMON /XYZ   / XINT(31213),YINT(31213),ZINT(31213)               
C$omp threadprivate(/XYZ/)
C                                                                       
      PARAMETER (ZERO=0.0D+00)                                          
      PARAMETER (ONE=1.0D+00)                                           
C                                                                       
      N0 = NMAX .EQ. 0                                                  
      N1 = NMAX .LE. 1                                                  
      M0 = MMAX .EQ. 0                                                  
      M1 = MMAX .LE. 1                                                  
C                                                                       
C     ----- I(0,0) -----                                                
C                                                                       
      I1 = I(1)                                                         
      XINT(I1) = ONE                                                    
      YINT(I1) = ONE                                                    
      ZINT(I1) = F00                                                    
      IF (N0 .AND. M0) RETURN                                           
      I2 = I(2)                                                         
      K2 = K(2)                                                         
      CP10 = B00                                                        
C                                                                       
C     ----- I(1,0) -----                                                
C                                                                       
      IF (.NOT. N0) THEN                                                
        XINT(I2) = XC00                                                 
        YINT(I2) = YC00                                                 
        ZINT(I2) = ZC00*F00                                             
        IF (M0) GO TO 120                                               
      END IF                                                            
C                                                                       
C     ----- I(0,1) -----                                                
C                                                                       
      I3 = I1+K2                                                        
      XINT(I3) = XCP00                                                  
      YINT(I3) = YCP00                                                  
      ZINT(I3) = ZCP00*F00                                              
C                                                                       
C     ----- I(1,1) -----                                                
C                                                                       
      IF (.NOT. N0) THEN                                                
        I3 = I2+K2                                                      
        XINT(I3) = XCP00*XINT(I2)+CP10                                  
        YINT(I3) = YCP00*YINT(I2)+CP10                                  
        ZINT(I3) = ZCP00*ZINT(I2)+CP10*F00                              
      END IF                                                            
C                                                                       
  120 CONTINUE                                                          
      IF (.NOT. N1) THEN                                                
        C10 = ZERO                                                      
        I3 = I1                                                         
        I4 = I2                                                         
        DO 160 N = 2,NMAX                                               
          C10 = C10+B10                                                 
C                                                                       
C     ----- I(N,0) -----                                                
C                                                                       
          I5 = I(N+1)                                                   
          XINT(I5) = C10*XINT(I3)+XC00*XINT(I4)                         
          YINT(I5) = C10*YINT(I3)+YC00*YINT(I4)                         
          ZINT(I5) = C10*ZINT(I3)+ZC00*ZINT(I4)                         
          IF ( .NOT. M0) THEN                                           
            CP10 = CP10+B00                                             
C                                                                       
C     ----- I(N,1) -----                                                
C                                                                       
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
C                                                                       
C     ----- I(0,M) -----                                                
C                                                                       
          I5 = I1+K(M+1)                                                
          XINT(I5) = CP01*XINT(I3)+XCP00*XINT(I4)                       
          YINT(I5) = CP01*YINT(I3)+YCP00*YINT(I4)                       
          ZINT(I5) = CP01*ZINT(I3)+ZCP00*ZINT(I4)                       
C                                                                       
C     ----- I(1,M) -----                                                
C                                                                       
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
C                                                                       
C     ----- I(N,M) -----                                                
C                                                                       
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
            XINT(I5+K4) = C10*XINT(I3+K4)+XC00*XINT(I4+K4)              
     *                    +C01*XINT(I4+K3)                              
            YINT(I5+K4) = C10*YINT(I3+K4)+YC00*YINT(I4+K4)              
     *                    +C01*YINT(I4+K3)                              
            ZINT(I5+K4) = C10*ZINT(I3+K4)+ZC00*ZINT(I4+K4)              
     *                    +C01*ZINT(I4+K3)                              
            C10 = C10+B10                                               
            I3 = I4                                                     
            I4 = I5                                                     
  260     CONTINUE                                                      
          K3 = K4                                                       
  280   CONTINUE                                                        
      END IF                                                            
C                                                                       
C     ----- I(NI,NJ,M) -----                                            
C                                                                       
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
C                                                                       
C     ----- I(NI,NJ,NK,NL) -----                                        
C                                                                       
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
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2A   *DECK XYZINT_gpu2                                      
      SUBROUTINE XYZINT_gpu2(XINT,YINT,ZINT,                            
     * I,K,NIMAX,NJMAX,NKMAX,NLMAX,NMAX,MMAX,                           
     * BP01,B00,B10,XCP00,XC00,YCP00,YC00,ZCP00,ZC00,F00,               
     * DXIJ,DYIJ,DZIJ,DXKL,DYKL,DZKL)                                   
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
C      COMMON /SETINT/ I(13),K(13),NIMAX,NJMAX,NKMAX,NLMAX,NMAX,MMAX    
C     +               ,BP01,B00,B10,XCP00,XC00,YCP00,YC00,ZCP00,ZC00,F00
C     +               ,DXIJ,DYIJ,DZIJ,DXKL,DYKL,DZKL                    
C      COMMON /XYZ   / XINT(31213),YINT(31213),ZINT(31213)              
                                                                        
      LOGICAL N0,N1,M0,M1,FIRST1,FIRST2,FIRST3,FIRST4                   
      DIMENSION :: I(13),K(13)                                          
      INTEGER :: NIMAX,NJMAX,NKMAX,NLMAX,NMAX,MMAX                      
      double precision :: BP01,B00,B10,XCP00,XC00,YCP00,YC00,ZCP00,ZC00 
      double precision :: F00,DXIJ,DYIJ,DZIJ,DXKL,DYKL,DZKL             
      DIMENSION :: XINT(31213),YINT(31213),ZINT(31213)                  
C                                                                       
      PARAMETER (ZERO=0.0D+00)                                          
      PARAMETER (ONE=1.0D+00)                                           
C                                                                       
      N0 = NMAX .EQ. 0                                                  
      N1 = NMAX .LE. 1                                                  
      M0 = MMAX .EQ. 0                                                  
      M1 = MMAX .LE. 1                                                  
C                                                                       
C     ----- I(0,0) -----                                                
C                                                                       
      I1 = I(1)                                                         
      XINT(I1) = ONE                                                    
      YINT(I1) = ONE                                                    
      ZINT(I1) = F00                                                    
      IF (N0 .AND. M0) RETURN                                           
      I2 = I(2)                                                         
      K2 = K(2)                                                         
      CP10 = B00                                                        
C                                                                       
C     ----- I(1,0) -----                                                
C                                                                       
      IF (.NOT. N0) THEN                                                
        XINT(I2) = XC00                                                 
        YINT(I2) = YC00                                                 
        ZINT(I2) = ZC00*F00                                             
        IF (M0) GO TO 120                                               
      END IF                                                            
C                                                                       
C     ----- I(0,1) -----                                                
C                                                                       
      I3 = I1+K2                                                        
      XINT(I3) = XCP00                                                  
      YINT(I3) = YCP00                                                  
      ZINT(I3) = ZCP00*F00                                              
C                                                                       
C     ----- I(1,1) -----                                                
C                                                                       
      IF (.NOT. N0) THEN                                                
        I3 = I2+K2                                                      
        XINT(I3) = XCP00*XINT(I2)+CP10                                  
        YINT(I3) = YCP00*YINT(I2)+CP10                                  
        ZINT(I3) = ZCP00*ZINT(I2)+CP10*F00                              
      END IF                                                            
C                                                                       
  120 CONTINUE                                                          
      IF (.NOT. N1) THEN                                                
        C10 = ZERO                                                      
        I3 = I1                                                         
        I4 = I2                                                         
        DO 160 N = 2,NMAX                                               
          C10 = C10+B10                                                 
C                                                                       
C     ----- I(N,0) -----                                                
C                                                                       
          I5 = I(N+1)                                                   
          XINT(I5) = C10*XINT(I3)+XC00*XINT(I4)                         
          YINT(I5) = C10*YINT(I3)+YC00*YINT(I4)                         
          ZINT(I5) = C10*ZINT(I3)+ZC00*ZINT(I4)                         
          IF ( .NOT. M0) THEN                                           
            CP10 = CP10+B00                                             
C                                                                       
C     ----- I(N,1) -----                                                
C                                                                       
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
C                                                                       
C     ----- I(0,M) -----                                                
C                                                                       
          I5 = I1+K(M+1)                                                
          XINT(I5) = CP01*XINT(I3)+XCP00*XINT(I4)                       
          YINT(I5) = CP01*YINT(I3)+YCP00*YINT(I4)                       
          ZINT(I5) = CP01*ZINT(I3)+ZCP00*ZINT(I4)                       
C                                                                       
C     ----- I(1,M) -----                                                
C                                                                       
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
C                                                                       
C     ----- I(N,M) -----                                                
C                                                                       
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
            XINT(I5+K4) = C10*XINT(I3+K4)+XC00*XINT(I4+K4)              
     *                    +C01*XINT(I4+K3)                              
            YINT(I5+K4) = C10*YINT(I3+K4)+YC00*YINT(I4+K4)              
     *                    +C01*YINT(I4+K3)                              
            ZINT(I5+K4) = C10*ZINT(I3+K4)+ZC00*ZINT(I4+K4)              
     *                    +C01*ZINT(I4+K3)                              
            C10 = C10+B10                                               
            I3 = I4                                                     
            I4 = I5                                                     
  260     CONTINUE                                                      
          K3 = K4                                                       
  280   CONTINUE                                                        
      END IF                                                            
C                                                                       
C     ----- I(NI,NJ,M) -----                                            
C                                                                       
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
C                                                                       
C     ----- I(NI,NJ,NK,NL) -----                                        
C                                                                       
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
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2A   *DECK ZPKOUT                                           
      SUBROUTINE ZPKOUT(II,JJ,KK,LL,GHONDO,SKIPA,SKIPB,SKIPC,NPSYM)     
      use mx_limits, only: mxsh,mxgtot                                  
C                                                                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C                                                                       
      DIMENSION GHONDO(*)                                               
      DIMENSION IB(28),JB(28),KB(28),LB(28)                             
C                                                                       
      LOGICAL SKIPA,SKIPB,SKIPC,NPSYM,IANDJ,KANDL,FIRST                 
C                                                                       
C                                                                       
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),      
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),      
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),   
     *                KLOC(MXSH),MIN(MXSH),MAX(MXSH),NSHELL             
      COMMON /SHLEXC/ NORGSH(3),NORGSP(3),IEXCH,NANGM,NGTH(4)           
C                                                                       
      PARAMETER (ZERO=0.0D+00)                                          
C                                                                       
      SAVE FIRST,IB,JB,KB,LB                                            
      DATA FIRST/.TRUE./                                                
C                                                                       
C     ----- ZERO HONDO PK INTEGRAL FORMATION REGION -----               
C                                                                       
      IF(FIRST) THEN                                                    
         FIRST=.FALSE.                                                  
         DO 100 I=1,NANGM                                               
            LB(I) = I-1                                                 
            KB(I) = LB(I) * NANGM                                       
            JB(I) = KB(I) * NANGM                                       
            IB(I) = JB(I) * NANGM                                       
 100     CONTINUE                                                       
      END IF                                                            
C                                                                       
      NORG1 = NORGSH(1)+1                                               
      NORG2 = NORGSH(2)+1                                               
      NORG3 = NORGSH(3)+1                                               
C                                                                       
      LIT = KTYPE(II)                                                   
      MINI = MIN(II)                                                    
      MAXI = MAX(II)                                                    
      MINJ = MIN(JJ)                                                    
      MAXJ = MAX(JJ)                                                    
      LKT = KTYPE(KK)                                                   
      MINK = MIN(KK)                                                    
      MAXK = MAX(KK)                                                    
      MINL = MIN(LL)                                                    
      MAXL = MAX(LL)                                                    
      IANDJ = II .EQ. JJ                                                
      KANDL = KK .EQ. LL                                                
C                                                                       
C     TYPE = 1 FOR (II II II II)                                        
C            2     (II JJ JJ JJ)                                        
C            3     (II II KK KK) AND  LIT.GE.LKT                        
C            4     (II II KK KK) AND  LIT.LT.LKT                        
C            5     (II II II LL)                                        
C            6     (II JJ KK KK)                                        
C            7     (II JJ JJ LL)                                        
C            8     (II II KK LL)                                        
C            9     (II JJ KK LL)                                        
C                                                                       
      NTYP = 0                                                          
      IF(II.EQ.JJ .AND. JJ.EQ.KK .AND. KK.EQ.LL) NTYP = 1               
      IF(II.GT.JJ .AND. JJ.EQ.KK .AND. KK.EQ.LL) NTYP = 2               
      IF(II.EQ.JJ .AND. JJ.GT.KK .AND. KK.EQ.LL                         
     *                         .AND. LIT.GE.LKT) NTYP = 3               
      IF(II.EQ.JJ .AND. JJ.GT.KK .AND. KK.EQ.LL                         
     *                         .AND. LIT.LT.LKT) NTYP = 4               
      IF(II.EQ.JJ .AND. JJ.EQ.KK .AND. KK.GT.LL) NTYP = 5               
      IF(II.GT.JJ .AND. JJ.GT.KK .AND. KK.EQ.LL) NTYP = 6               
      IF(II.GT.JJ .AND. JJ.EQ.KK .AND. KK.GT.LL) NTYP = 7               
      IF(II.EQ.JJ .AND. JJ.GT.KK .AND. KK.GT.LL) NTYP = 8               
      IF(II.GT.JJ .AND. JJ.GT.KK .AND. KK.GT.LL) NTYP = 9               
C                                                                       
      IF(SKIPA .AND. .NOT.NPSYM) NORG2 = 1                              
      IF(SKIPC .AND. .NOT.NPSYM) NORG3 = NORGSH(2)+1                    
      IF(SKIPB .AND. .NOT.NPSYM) NORG3 = 1                              
C                                                                       
C     ----- BEGIN LOOPS OVER PRIMITIVES IN THIS SHELL -----             
C                                                                       
C     INTEGRAL TYPES N1,G1 FOR (I,J//K,L)                               
C                    N2,G2 FOR (I,K//J,L)                               
C                    N3,G3 FOR (I,L//J,K)                               
C                                                                       
C                      HONDO INTEGRALS                                  
C        N1 = IB(IA)+JB(JA)+KB(KA)+LB(LA)+NORG1                         
C        N2 = IB(IA)+JB(KA)+KB(JA)+LB(LA)+NORG2                         
C                                                                       
      JMAX = MAXJ                                                       
      KMAX = MAXK                                                       
      LMAX = MAXL                                                       
      DO 860 I = MINI,MAXI                                              
      IA = I-MINI+1                                                     
      N1I = NORG1 + IB(IA)                                              
      N2I = NORG2 + IB(IA)                                              
C                                                                       
      IF (IANDJ) JMAX = I                                               
      DO 840 J = MINJ,JMAX                                              
      IF (JJ .EQ. KK) KMAX = J                                          
      JA = J-MINJ+1                                                     
      N1IJ = N1I + JB(JA)                                               
      N2IJ = N2I + KB(JA)                                               
C                                                                       
      DO 820 K = MINK,KMAX                                              
      KA = K-MINK+1                                                     
      N1IJK = N1IJ + KB(KA)                                             
      N2IJK = N2IJ + JB(KA)                                             
C                                                                       
      IF (KANDL) LMAX = K                                               
      DO 800 L = MINL,LMAX                                              
      LA = L-MINL+1                                                     
      N1 = N1IJK + LB(LA)                                               
      N2 = N2IJK + LB(LA)                                               
C                                                                       
      GO TO (200,220,230,250,270,280,290,300,310),NTYP                  
  200 IF (IA .EQ. JA) GO TO 210                                         
      N3 = IB(IA)+JB(LA)+KB(JA)+LB(KA)+NORG3                            
      GO TO 400                                                         
  210 N3 = IB(JA)+JB(KA)+KB(IA)+LB(LA)+NORG3                            
      GO TO 400                                                         
  220 N3 = IB(IA)+JB(LA)+KB(JA)+LB(KA)+NORG3                            
      GO TO 400                                                         
  230 IF (IA .EQ. JA) GO TO 240                                         
      N3 = IB(IA)+JB(LA)+KB(JA)+LB(KA)+NORG3                            
      GO TO 400                                                         
  240 N3 = IB(JA)+JB(KA)+KB(IA)+LB(LA)+NORG3                            
      GO TO 400                                                         
  250 IF (KA .EQ. LA) GO TO 260                                         
      N3 = IB(JA)+JB(KA)+KB(IA)+LB(LA)+NORG3                            
      GO TO 400                                                         
  260 N3 = IB(IA)+JB(LA)+KB(JA)+LB(KA)+NORG3                            
      GO TO 400                                                         
  270 N3 = IB(JA)+JB(KA)+KB(IA)+LB(LA)+NORG3                            
      GO TO 400                                                         
  280 N3 = IB(IA)+JB(LA)+KB(JA)+LB(KA)+NORG3                            
      GO TO 400                                                         
  290 N3 = IB(IA)+JB(LA)+KB(JA)+LB(KA)+NORG3                            
      GO TO 400                                                         
  300 N3 = IB(JA)+JB(KA)+KB(IA)+LB(LA)+NORG3                            
      GO TO 400                                                         
  310 N3 = IB(IA)+JB(LA)+KB(JA)+LB(KA)+NORG3                            
C                                                                       
  400 CONTINUE                                                          
      GHONDO(N1) = ZERO                                                 
      GHONDO(N2) = ZERO                                                 
      GHONDO(N3) = ZERO                                                 
  800 CONTINUE                                                          
  820 CONTINUE                                                          
  840 CONTINUE                                                          
  860 CONTINUE                                                          
      RETURN                                                            
      END                                                               
C*MODULE INT2A   *DECK ZQOUT                                            
      SUBROUTINE ZQOUT(GHONDO)                                          
C                                                                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C                                                                       
      DIMENSION GHONDO(*)                                               
C                                                                       
      LOGICAL IANDJ,KANDL,SAME                                          
C                                                                       
      COMMON /INTDEX/ IJGT(784),IJX(784),IJY(784),IJZ(784),IK(784),     
     *                KLGT(784),KLX(784),KLY(784),KLZ(784)              
C$omp threadprivate(/INTDEX/)
      COMMON /MISC  / IANDJ,KANDL,SAME                                  
C$omp threadprivate(/MISC/)
      COMMON /SHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,          
     *                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,          
     *                NIJX,IJ,KL,IJKL                                   
C$omp threadprivate(/SHLNOS/)
C                                                                       
      PARAMETER (ZERO=0.0D+00)                                          
C                                                                       
C     ----- ZERO HONDO CONVENTIONAL INTEGRAL OUTPUT REGION -----        
C                                                                       
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
                  IF (SAME .AND. KLN .GT. IJN) GO TO 240                
                  NN = N1+KLGT(KLN)                                     
                  GHONDO(NN) = ZERO                                     
  200          CONTINUE                                                 
  220       CONTINUE                                                    
  240    CONTINUE                                                       
  260 CONTINUE                                                          
      RETURN                                                            
      END                                                               
C*MODULE INT2A   *DECK ZQOUT_gpu2                                       
      SUBROUTINE ZQOUT_gpu2(GHONDO)                                     
C                                                                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C                                                                       
      DIMENSION GHONDO(*)                                               
C                                                                       
      LOGICAL IANDJ,KANDL,SAME                                          
C                                                                       
      COMMON /INTDEX/ IJGT(784),IJX(784),IJY(784),IJZ(784),IK(784),     
     *                KLGT(784),KLX(784),KLY(784),KLZ(784)              
C$omp threadprivate(/INTDEX/)
      COMMON /MISC  / IANDJ,KANDL,SAME                                  
C$omp threadprivate(/MISC/)
      COMMON /SHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,          
     *                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,          
     *                NIJX,IJ,KL,IJKL                                   
C$omp threadprivate(/SHLNOS/)
C                                                                       
      PARAMETER (ZERO=0.0D+00)                                          
C                                                                       
C     ----- ZERO HONDO CONVENTIONAL INTEGRAL OUTPUT REGION -----        
C                                                                       
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
                  IF (SAME .AND. KLN .GT. IJN) GO TO 240                
                  NN = N1+KLGT(KLN)                                     
                  GHONDO(NN) = ZERO                                     
  200          CONTINUE                                                 
  220       CONTINUE                                                    
  240    CONTINUE                                                       
  260 CONTINUE                                                          
      RETURN                                                            
      END                                                               
C*MODULE INT2A   *DECK INT2EIC                                          
      SUBROUTINE INT2EIC(IMODE)                                         
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      LOGICAL PACK2E                                                    
      LOGICAL GOPARR,DSKWRK,MASWRK                                      
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)                
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
      COMMON /INTFIL/ NINTMX,NHEX,NTUPL,PACK2E,INTTYP,IGRDTYP           
      COMMON /INT2IC/ NINTIC,ININTIC,NXXIC,LBUFPIC,LIXIC,LABSIX,NINTIX  
      COMMON /MACHIN/ NWDVAR,MAXFM,MAXSM,LIMFM,LIMSM                    
      COMMON /PCKLAB/ LABSIZ                                            
      SAVE NEED                                                         
C                                                                       
C     INITIALISE(IMODE=0)/FINALISE(IMODE=1) IN CORE 2E INTEGRALS        
C     IF NINTIC IS NEGATIVE,                                            
C         IT IS ASSUMED TO BE THE MEMORY AMOUNT TO BE USED              
C     IF NINTIC IS POSITIVE,                                            
C         IT IS THE NUMBER OF INTEGRALS TO BE STORED IN CORE            
C     IF NINTIC IS 0,                                                   
C         NONE OF THE INTEGRALS WILL BE STORED IN CORE.                 
C                                                                       
C     GAMESS SEEMS TO PAY NO ATTENTION ALLOCATING INTEGRAL INDEX        
C     BUFFER IGNORING LABSIZ. BELOW WE DO STRICTLY CORRECT USAGE        
C     ALLOCATING MINIMUM POSSIBLE SIZE. IT MAY NOT WORK IF WHEREVER     
C     IN-CORE INTEGRALS ARE USED THE BUFFERS ARE REUSED FOR SOMETHING   
C     ELSE WITH INCORRECT ASSUMPTION ABOUT THE SIZE SINCE THE BUFFERS   
C     ALLOCATED BELOW WILL BE USED INSTEAD OF WHATEVER IS ALLOCATED     
C     IN THE MAIN BODY.                                                 
C     SEE ALSO OPTFMOX IF MEMORY ALLOCATION IN THIS SUBROUTINE CHANGES. 
C                                                                       
      IF(NINTIC.EQ.0) RETURN                                            
      LABSIZ2=2/LABSIZ                                                  
      IF(LABSIZ.NE.1.AND.LABSIZ.NE.2) CALL ABRT                         
C     THIS WILL WORK FOR LABSIZ=1 AND 2 (ONLY).                         
C                                                                       
      IF(NINTIC.LT.0) NINTIC=(-NINTIC*LABSIZ2)/(LABSIZ2+1)              
C                                                                       
C     ININTIC IS USED TO ADDRESS INTEGER ARRAYS TO PROPERLY SHIFT       
C     INDICES, HAVING THEM POINTING WHERE DISK DRIVEN INTEGRALS START   
C                                                                       
      ININTIC=NINTIC                                                    
      IF(LABSIZ.EQ.2.AND.NWDVAR.EQ.2) ININTIC=NINTIC*2                  
      IF(LABSIZ.EQ.1.AND.NWDVAR.EQ.1) THEN                              
        IF(MOD(NINTIC,2).EQ.1) THEN                                     
C         IT IS NOT THAT DIFFICULT TO IMPLEMENT THIS                    
C         (ADD SHIFTING 4 BYTES IN ADDRESSING 8-BYTE INTEGERS INDICES), 
C         BUT WHY BOTHER?                                               
          NINTIC=NINTIC-1                                               
C         IF(MASWRK) WRITE(IW,9020)                                     
C         CALL ABRT                                                     
C       ELSE                                                            
        ENDIF                                                           
        ININTIC=NINTIC/2                                                
      ENDIF                                                             
      IF(IMODE.EQ.0) THEN                                               
        CALL VALFM(LOADFM)                                              
        LBUFPIC= LOADFM + 1                                             
        LIXIC=   LBUFPIC+ NINTMX+NINTIC                                 
        LAST=    LIXIC  + (NINTMX+NINTIC-1)/LABSIZ2+1                   
        NEED=    LAST - LOADFM - 1                                      
        CALL GETFM(NEED)                                                
        IF(MASWRK) WRITE(IW,9000) NEED                                  
      ELSE                                                              
        CALL RETFM(NEED)                                                
        IF(MASWRK) WRITE(IW,9010)                                       
      ENDIF                                                             
 9000 FORMAT(/1X,'ALLOCATED ',I12,' WORDS FOR IN CORE 2E INTEGRALS.',/) 
 9010 FORMAT(1X,'RETURNED IN CORE 2E INTEGRAL BUFFER.',/)               
C9020 FORMAT(/1X,'PLEASE SELECT AN EVEN NINTIC IF ',                    
C    *          'LABSIZ.EQ.1.AND.NWDVAR.EQ.1',/)                        
      RETURN                                                            
      END                                                               
C*MODULE INT2A   *DECK GENRAL_gpu                                       
      SUBROUTINE GENRAL_gpu(GHONDO,DDIJ,DKL,DIJ,                        
     * IJGT,IJX,IJY,IJZ,IK,KLGT,KLX,KLY,KLZ,                            
     * AA,R,X1,Y1,Z1,IJD,IANDJ,KANDL,SAME,                              
     * LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,                             
     * qq4,mini,maxi,minj,maxj,mink,maxk,minl,maxl,                     
     * NIJ,IJ,KL,IJKL,                                                  
     * NROOTS,                                                          
     * ag,csa,cpa,cda,cfa,cga,cha,cia,                                  
     * bg,csb,cpb,cdb,cfb,cgb,chb,cib,                                  
     * cg,csc,cpc,cdc,cfc,cgc,chc,cic,                                  
     * dg,csd,cpd,cdd,cfd,cgd,chd,cid,                                  
     * XI,YI,ZI,XJ,YJ,ZJ,RRI,XK,YK,ZK,XL,YL,ZL,RRK,                     
     * NGA,NGB,NGC,NGD)                                                 
                                                                        
                                                                        
      USE lrcdft, ONLY: EMU2                                            
      use mx_limits, only: mxgsh,mxg2                                   
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      DIMENSION GHONDO(*),DDIJ(*)                                       
C                                                                       
      LOGICAL IANDJ,KANDL,SAME,OUT,NORM,DOUBLE                          
      DIMENSION :: XIN(31213),YIN(31213),ZIN(31213)                     
      DIMENSION :: DKL(784),DIJ(784)                                    
      DIMENSION :: IJGT(784),IJX(784),IJY(784),IJZ(784),IK(784)         
      DIMENSION :: KLGT(784),KLX(784),KLY(784),KLZ(784)                 
      DIMENSION :: AA(MXG2),R(MXG2),X1(MXG2),Y1(MXG2),Z1(MXG2),IJD(784) 
                                                                        
                                                                        
C      COMMON /MISC  / IANDJ,KANDL,SAME                                 
C      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK                
C      COMMON /ROOT  / XX,U(13),W(13),NROOTS                            
      double precision :: XX                                            
      double precision :: U(13)                                         
      double precision :: W(13)                                         
      integer :: NROOTS                                                 
C      COMMON /SETINT/ IN(13),KN(13),NI,NJ,NK,NL,NMAX,MMAX,             
C     +                BP01,B00,B10,XCP00,XC00,YCP00,YC00,ZCP00,ZC00,F00
C     +                DXIJ,DYIJ,DZIJ,DXKL,DYKL,DZKL                    
      DIMENSION :: IN(13),KN(13)                                        
      INTEGER :: NI,NJ,NK,NL,NMAX,MMAX                                  
      double precision :: BP01,B00,B10,XCP00,XC00,YCP00,YC00,ZCP00,ZC00 
      double precision :: F00,DXIJ,DYIJ,DZIJ,DXKL,DYKL,DZKL             
C      COMMON /SHLINF/  AG(MXGSH),CSA(MXGSH),CPA(MXGSH),CDA(MXGSH),     
C     *                CFA(MXGSH),CGA(MXGSH),CHA(MXGSH),CIA(MXGSH),     
C     *                 BG(MXGSH),CSB(MXGSH),CPB(MXGSH),CDB(MXGSH),     
C     *                CFB(MXGSH),CGB(MXGSH),CHB(MXGSH),CIB(MXGSH),     
C     *                 CG(MXGSH),CSC(MXGSH),CPC(MXGSH),CDC(MXGSH),     
C     *                CFC(MXGSH),CGC(MXGSH),CHC(MXGSH),CIC(MXGSH),     
C     *                 DG(MXGSH),CSD(MXGSH),CPD(MXGSH),CDD(MXGSH),     
C     *                CFD(MXGSH),CGD(MXGSH),CHD(MXGSH),CID(MXGSH),     
C     *                XI,YI,ZI,XJ,YJ,ZJ,RRI,XK,YK,ZK,XL,YL,ZL,RRK,     
C     *                NGA,NGB,NGC,NGD                                  
                                                                        
      DIMENSION :: AG(MXGSH),CSA(MXGSH),CPA(MXGSH),CDA(MXGSH)           
      DIMENSION :: CFA(MXGSH),CGA(MXGSH),CHA(MXGSH),CIA(MXGSH)          
      DIMENSION :: BG(MXGSH),CSB(MXGSH),CPB(MXGSH),CDB(MXGSH)           
      DIMENSION :: CFB(MXGSH),CGB(MXGSH),CHB(MXGSH),CIB(MXGSH)          
      DIMENSION :: CG(MXGSH),CSC(MXGSH),CPC(MXGSH),CDC(MXGSH)           
      DIMENSION :: CFC(MXGSH),CGC(MXGSH),CHC(MXGSH),CIC(MXGSH)          
      DIMENSION :: DG(MXGSH),CSD(MXGSH),CPD(MXGSH),CDD(MXGSH)           
      DIMENSION :: CFD(MXGSH),CGD(MXGSH),CHD(MXGSH),CID(MXGSH)          
      double precision :: XI,YI,ZI,XJ,YJ,ZJ,RRI,XK,YK,ZK,XL,YL,ZL,RRK   
      integer :: NGA,NGB,NGC,NGD                                        
                                                                        
                                                                        
C      COMMON /SHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,         
C     +                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,         
C     +                NIJ,IJ,KL,IJKL                                   
      double precision :: qq4                                           
      INTEGER :: LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL                    
      INTEGER :: MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL                
      INTEGER :: NIJ,IJ,KL,IJKL                                         
C                                                                       
      DIMENSION IN1(13)                                                 
C                                                                       
      PARAMETER (SQRT3=1.73205080756888D+00, SQRT5=2.23606797749979D+00,
     *           SQRT7=2.64575131106459D+00, PI252=34.986836655250D+00, 
     *           SQRT9=3.0D+00,SQRT11=3.3166247903553998D+00,           
     *           ZERO=0.0D+00, HALF=0.5D+00, ONE=1.0D+00,               
     *           TOL=46.0515999999999934 )                              
C                                                                       
C     GENERAL INTEGRAL ROUTINE FOR SPDFGHI AND L FUNCTIONS              
C                                                                       
!!$omp declare target                                                   
      FACTOR = PI252*QQ4                                                
C      NORM = NORMF .NE. 1 .OR. NORMP .NE. 1                            
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
C                                                                       
C     ----- K PRIMITIVE                                                 
C                                                                       
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
         CFK = CFC(KG)*FACTOR                                           
         CGK = CGC(KG)*FACTOR                                           
         CHK = CHC(KG)*FACTOR                                           
         CIK = CIC(KG)*FACTOR                                           
C                                                                       
C        ----- L PRIMITIVE                                              
C                                                                       
         IF (KANDL) LGMAX = KG                                          
         DO 460 LG = 1,LGMAX                                            
            AL = DG(LG)                                                 
            B = AK+AL                                                   
            BINV = ONE/B                                                
            BBRRK = AL*BRRK*BINV                                        
            IF (BBRRK .GT. TOL) GO TO 460                               
            CSL = CSD(LG)                                               
            CPL = CPD(LG)                                               
            CDL = CDD(LG)                                               
            CFL = CFD(LG)                                               
            CGL = CGD(LG)                                               
            CHL = CHD(LG)                                               
            CIL = CID(LG)                                               
            XB = (AKXK+AL*XL)*BINV                                      
            YB = (AKYK+AL*YL)*BINV                                      
            ZB = (AKZK+AL*ZL)*BINV                                      
            BXBK = B*(XB-XK)                                            
            BYBK = B*(YB-YK)                                            
            BZBK = B*(ZB-ZK)                                            
            BXBI = B*(XB-XI)                                            
            BYBI = B*(YB-YI)                                            
            BZBI = B*(ZB-ZI)                                            
C                                                                       
C           ----- DENSITY FACTOR                                        
C                                                                       
            DOUBLE=KANDL.AND.KG.NE.LG                                   
            N = 0                                                       
            MAX = MAXL                                                  
            DUM1 = ZERO                                                 
            DUM2 = ZERO                                                 
            DO 370 K = MINK,MAXK                                        
               GO TO (140,160,220,220,180,220,220,200,220,220,          
     1                201,220,220,202,220,220,220,220,220,203,          
     1                204,220,220,205,220,220,220,220,220,206,          
     1                220,220,207,220,220,                              
     1                208,220,220,209,220,220,220,220,220,210,          
     1                220,220,220,220,220,211,220,220,212,220,          
     1                220,                                              
     1                213,220,220,214,220,220,220,220,220,215,          
     1                220,220,220,220,220,216,220,220,217,220,          
     1                220,218,220,220,220,220,220,219),K                
  140          DUM1 = CSK*BINV                                          
               GO TO 220                                                
  160          DUM1 = CPK*BINV                                          
               GO TO 220                                                
  180          DUM1 = CDK*BINV                                          
               GO TO 220                                                
  200          DUM1 = DUM1*SQRT3                                        
               GO TO 220                                                
  201          DUM1 = CFK*BINV                                          
               GO TO 220                                                
  202          DUM1 = DUM1*SQRT5                                        
               GO TO 220                                                
  203          DUM1 = DUM1*SQRT3                                        
               GO TO 220                                                
  204          DUM1 = CGK*BINV                                          
               GO TO 220                                                
  205          DUM1 = DUM1*SQRT7                                        
               GO TO 220                                                
  206          DUM1 = DUM1*SQRT5/SQRT3                                  
               GO TO 220                                                
  207          DUM1 = DUM1*SQRT3                                        
               GO TO 220                                                
  208          DUM1 = CHK*BINV                                          
               GO TO 220                                                
  209          DUM1 = DUM1*SQRT9                                        
               GO TO 220                                                
  210          DUM1 = DUM1*SQRT7/SQRT3                                  
               GO TO 220                                                
  211          DUM1 = DUM1*SQRT3                                        
               GO TO 220                                                
  212          DUM1 = DUM1*SQRT5/SQRT3                                  
               GO TO 220                                                
  213          DUM1 = CIK*BINV                                          
               GO TO 220                                                
  214          DUM1 = DUM1*SQRT11                                       
               GO TO 220                                                
  215          DUM1 = DUM1*SQRT3                                        
               GO TO 220                                                
  216          DUM1 = DUM1*SQRT3                                        
               GO TO 220                                                
  217          DUM1 = DUM1*SQRT7/(SQRT5*SQRT3)                          
               GO TO 220                                                
  218          DUM1 = DUM1*SQRT5                                        
               GO TO 220                                                
  219          DUM1 = DUM1*SQRT5/SQRT3                                  
C                                                                       
  220          IF (KANDL) MAX = K                                       
               DO 360 L = MINL,MAX                                      
                  GO TO (240,280,340,340,300,340,340,320,340,340,       
     1                   321,340,340,322,340,340,340,340,340,323,       
     1                   324,340,340,325,340,340,340,340,340,326,       
     1                   340,340,327,340,340,                           
     1                   328,340,340,329,340,340,340,340,340,330,       
     1                   340,340,340,340,340,331,340,340,332,340,       
     1                   340,                                           
     1                   333,340,340,334,340,340,340,340,340,335,       
     1                   340,340,340,340,340,336,340,340,337,340,       
     1                   340,338,340,340,340,340,340,339),L             
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
  321             DUM2 = DUM1*CFL                                       
                  IF (DOUBLE) DUM2 = DUM2+DUM2                          
                  GO TO 340                                             
  322             DUM2 = DUM2*SQRT5                                     
                  GO TO 340                                             
  323             DUM2 = DUM2*SQRT3                                     
                  GO TO 340                                             
  324             DUM2 = DUM1*CGL                                       
                  IF (DOUBLE) DUM2 = DUM2+DUM2                          
                  GO TO 340                                             
  325             DUM2 = DUM2*SQRT7                                     
                  GO TO 340                                             
  326             DUM2 = DUM2*SQRT5/SQRT3                               
                  GO TO 340                                             
  327             DUM2 = DUM2*SQRT3                                     
                  GO TO 340                                             
  328             DUM2 = DUM1*CHL                                       
                  IF (DOUBLE) DUM2 = DUM2+DUM2                          
                  GO TO 340                                             
  329             DUM2 = DUM2*SQRT9                                     
                  GO TO 340                                             
  330             DUM2 = DUM2*SQRT7/SQRT3                               
                  GO TO 340                                             
  331             DUM2 = DUM2*SQRT3                                     
                  GO TO 340                                             
  332             DUM2 = DUM2*SQRT5/SQRT3                               
                  GO TO 340                                             
  333             DUM2 = DUM1*CIL                                       
                  IF (DOUBLE) DUM2 = DUM2+DUM2                          
                  GO TO 340                                             
  334             DUM2 = DUM2*SQRT11                                    
                  GO TO 340                                             
  335             DUM2 = DUM2*SQRT3                                     
                  GO TO 340                                             
  336             DUM2 = DUM2*SQRT3                                     
                  GO TO 340                                             
  337             DUM2 = DUM2*SQRT7/(SQRT5*SQRT3)                       
                  GO TO 340                                             
  338             DUM2 = DUM2*SQRT5                                     
                  GO TO 340                                             
  339             DUM2 = DUM2*SQRT5/SQRT3                               
C                                                                       
  340             N = N+1                                               
                  DKL(N) = DUM2                                         
  360          CONTINUE                                                 
  370       CONTINUE                                                    
C                                                                       
C           ----- PAIR OF I,J PRIMITIVES                                
C                                                                       
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
C               IF(LRINT) THEN                                          
C                 RHO0 = RHO                                            
C                 RHO  = RHO0*EMU2/(RHO0+EMU2)                          
C               ENDIF                                                   
               XA = X1(N)                                               
               YA = Y1(N)                                               
               ZA = Z1(N)                                               
               XX = RHO*((XA-XB)*(XA-XB) + (YA-YB)*(YA-YB)              
     *                                   + (ZA-ZB)*(ZA-ZB))             
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
C                                                                       
C              ----- ROOTS AND WEIGHTS FOR QUADRATURE                   
               IF (NROOTS .LE. 3) THEN                                  
                   CALL RT123TS(XX,NROOTS,RT1,RT2,RT3,WW1,WW2,WW3)      
                   U(1)=RT1                                             
                   U(2)=RT2                                             
                   U(3)=RT3                                             
                   W(1)=WW1                                             
                   W(2)=WW2                                             
                   W(3)=WW3                                             
               ENDIF                                                    
               IF (NROOTS .EQ. 4) THEN                                  
                   CALL ROOT4TS(XX,NROOTS,RT1,RT2,RT3,RT4,              
     * WW1,WW2,WW3,WW4)                                                 
                   U(1)=RT1                                             
                   U(2)=RT2                                             
                   U(3)=RT3                                             
                   U(4)=RT4                                             
                   W(1)=WW1                                             
                   W(2)=WW2                                             
                   W(3)=WW3                                             
                   W(4)=WW4                                             
               ENDIF                                                    
               IF (NROOTS .EQ. 5) THEN                                  
                   CALL ROOT5TS(XX,NROOTS,RT1,RT2,RT3,RT4,RT5,          
     * WW1,WW2,WW3,WW4,WW5)                                             
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
               ENDIF                                                    
               MM = 0                                                   
               MAX = NMAX+1                                             
C                                                                       
C              COMPUTE TWO-ELECTRON INTEGRALS FOR EACH ROOT             
C                                                                       
               DO 420 M = 1,NROOTS                                      
                  U2 = U(M)*RHO                                         
                  F00 = EXPE*W(M)                                       
C                  IF(LRINT) F00 = F00*SQRT(EMU2/(RHO0+EMU2))           
                  DO 400 I = 1,MAX                                      
                     IN(I) = IN1(I)+MM                                  
  400             CONTINUE                                              
C                  IF(.NOT.LRINT)THEN                                   
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
C                  ELSE                                                 
C                     T2    = U2/(U2+RHO)                               
C                     T2AR  = T2*RHO/A                                  
C                     T2BR  = T2*RHO/B                                  
C                     BP01  = HALF/B*(ONE-T2BR)                         
C                     B00   = HALF/AANDB*T2*RHO/RHO0                    
C                     B10   = HALF/A*(ONE-T2AR)                         
C                     XCP00 = (XB-XK)+T2BR*(XA-XB)                      
C                     XC00  = (XA-XI)-T2AR*(XA-XB)                      
C                     YCP00 = (YB-YK)+T2BR*(YA-YB)                      
C                     YC00  = (YA-YI)-T2AR*(YA-YB)                      
C                     ZCP00 = (ZB-ZK)+T2BR*(ZA-ZB)                      
C                     ZC00  = (ZA-ZI)-T2AR*(ZA-ZB)                      
C                  END IF                                               
                  CALL XYZINT_gpu(XIN,YIN,ZIN,                          
     * IN,KN,NI,NJ,NK,NL,NMAX,MMAX,                                     
     * BP01,B00,B10,XCP00,XC00,YCP00,YC00,ZCP00,ZC00,F00,               
     * DXIJ,DYIJ,DZIJ,DXKL,DYKL,DZKL)                                   
                  MM = MM+2401                                          
  420          CONTINUE                                                 
C                                                                       
C              ----- FORM (I,J//K,L) INTEGRALS OVER FUNCTIONS           
C                                                                       
               CALL FORMS_gpu(GHONDO,NROOTS,DKL,DIJ,XIN,YIN,ZIN,        
     * IJGT,IJX,IJY,IJZ,IK,KLGT,KLX,KLY,KLZ,IJ)                         
  440       NN = NN+49                                                  
  460    CONTINUE                                                       
  480 CONTINUE                                                          
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2A   *DECK GENRAL_gpu_5                                     
      SUBROUTINE GENRAL_gpu_5(GHONDO,DDIJ,DKL,DIJ,                      
     * IJGT,IJX,IJY,IJZ,IK,KLGT,KLX,KLY,KLZ,                            
     * AA,R,X1,Y1,Z1,IJD,IANDJ,KANDL,SAME,                              
     * LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,                             
     * qq4,mini,maxi,minj,maxj,mink,maxk,minl,maxl,                     
     * NIJ,IJ,KL,IJKL,                                                  
     * NROOTS,                                                          
     * ag,csa,cpa,cda,cfa,cga,cha,cia,                                  
     * bg,csb,cpb,cdb,cfb,cgb,chb,cib,                                  
     * cg,csc,cpc,cdc,cfc,cgc,chc,cic,                                  
     * dg,csd,cpd,cdd,cfd,cgd,chd,cid,                                  
     * XI,YI,ZI,XJ,YJ,ZJ,RRI,XK,YK,ZK,XL,YL,ZL,RRK,                     
     * NGA,NGB,NGC,NGD)                                                 
                                                                        
                                                                        
      USE lrcdft, ONLY: EMU2                                            
      use mx_limits, only: mxgsh,mxg2                                   
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      DIMENSION GHONDO(*),DDIJ(*)                                       
C                                                                       
      LOGICAL IANDJ,KANDL,SAME,OUT,NORM,DOUBLE                          
      DIMENSION :: XIN(31213),YIN(31213),ZIN(31213)                     
      DIMENSION :: DKL(784),DIJ(784)                                    
      DIMENSION :: IJGT(784),IJX(784),IJY(784),IJZ(784),IK(784)         
      DIMENSION :: KLGT(784),KLX(784),KLY(784),KLZ(784)                 
      DIMENSION :: AA(MXG2),R(MXG2),X1(MXG2),Y1(MXG2),Z1(MXG2),IJD(784) 
                                                                        
      double precision :: XX                                            
      double precision :: U(13)                                         
      double precision :: W(13)                                         
      integer :: NROOTS                                                 
      DIMENSION :: IN(13),KN(13)                                        
      INTEGER :: NI,NJ,NK,NL,NMAX,MMAX                                  
      double precision :: BP01,B00,B10,XCP00,XC00,YCP00,YC00,ZCP00,ZC00 
      double precision :: F00,DXIJ,DYIJ,DZIJ,DXKL,DYKL,DZKL             
                                                                        
      DIMENSION :: AG(MXGSH),CSA(MXGSH),CPA(MXGSH),CDA(MXGSH)           
      DIMENSION :: CFA(MXGSH),CGA(MXGSH),CHA(MXGSH),CIA(MXGSH)          
      DIMENSION :: BG(MXGSH),CSB(MXGSH),CPB(MXGSH),CDB(MXGSH)           
      DIMENSION :: CFB(MXGSH),CGB(MXGSH),CHB(MXGSH),CIB(MXGSH)          
      DIMENSION :: CG(MXGSH),CSC(MXGSH),CPC(MXGSH),CDC(MXGSH)           
      DIMENSION :: CFC(MXGSH),CGC(MXGSH),CHC(MXGSH),CIC(MXGSH)          
      DIMENSION :: DG(MXGSH),CSD(MXGSH),CPD(MXGSH),CDD(MXGSH)           
      DIMENSION :: CFD(MXGSH),CGD(MXGSH),CHD(MXGSH),CID(MXGSH)          
      double precision :: XI,YI,ZI,XJ,YJ,ZJ,RRI,XK,YK,ZK,XL,YL,ZL,RRK   
      integer :: NGA,NGB,NGC,NGD                                        
                                                                        
      double precision :: qq4                                           
      INTEGER :: LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL                    
      INTEGER :: MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL                
      INTEGER :: NIJ,IJ,KL,IJKL                                         
C                                                                       
      DIMENSION IN1(13)                                                 
C                                                                       
      PARAMETER (SQRT3=1.73205080756888D+00, SQRT5=2.23606797749979D+00,
     *           SQRT7=2.64575131106459D+00, PI252=34.986836655250D+00, 
     *           SQRT9=3.0D+00,SQRT11=3.3166247903553998D+00,           
     *           ZERO=0.0D+00, HALF=0.5D+00, ONE=1.0D+00,               
     *           TOL=46.0515999999999934 )                              
                                                                        
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
C                                                                       
C     ----- K PRIMITIVE                                                 
C                                                                       
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
         CFK = CFC(KG)*FACTOR                                           
         CGK = CGC(KG)*FACTOR                                           
         CHK = CHC(KG)*FACTOR                                           
         CIK = CIC(KG)*FACTOR                                           
C                                                                       
C        ----- L PRIMITIVE                                              
C                                                                       
         IF (KANDL) LGMAX = KG                                          
         DO 460 LG = 1,LGMAX                                            
            AL = DG(LG)                                                 
            B = AK+AL                                                   
            BINV = ONE/B                                                
            BBRRK = AL*BRRK*BINV                                        
            IF (BBRRK .GT. TOL) GO TO 460                               
            CSL = CSD(LG)                                               
            CPL = CPD(LG)                                               
            CDL = CDD(LG)                                               
            CFL = CFD(LG)                                               
            CGL = CGD(LG)                                               
            CHL = CHD(LG)                                               
            CIL = CID(LG)                                               
            XB = (AKXK+AL*XL)*BINV                                      
            YB = (AKYK+AL*YL)*BINV                                      
            ZB = (AKZK+AL*ZL)*BINV                                      
            BXBK = B*(XB-XK)                                            
            BYBK = B*(YB-YK)                                            
            BZBK = B*(ZB-ZK)                                            
            BXBI = B*(XB-XI)                                            
            BYBI = B*(YB-YI)                                            
            BZBI = B*(ZB-ZI)                                            
C                                                                       
C           ----- DENSITY FACTOR                                        
C                                                                       
            DOUBLE=KANDL.AND.KG.NE.LG                                   
            N = 0                                                       
            MAX = MAXL                                                  
            DUM1 = ZERO                                                 
            DUM2 = ZERO                                                 
            DO 370 K = MINK,MAXK                                        
               GO TO (140,160,220,220,180,220,220,200,220,220,          
     1                201,220,220,202,220,220,220,220,220,203,          
     1                204,220,220,205,220,220,220,220,220,206,          
     1                220,220,207,220,220,                              
     1                208,220,220,209,220,220,220,220,220,210,          
     1                220,220,220,220,220,211,220,220,212,220,          
     1                220,                                              
     1                213,220,220,214,220,220,220,220,220,215,          
     1                220,220,220,220,220,216,220,220,217,220,          
     1                220,218,220,220,220,220,220,219),K                
  140          DUM1 = CSK*BINV                                          
               GO TO 220                                                
  160          DUM1 = CPK*BINV                                          
               GO TO 220                                                
  180          DUM1 = CDK*BINV                                          
               GO TO 220                                                
  200          DUM1 = DUM1*SQRT3                                        
               GO TO 220                                                
  201          DUM1 = CFK*BINV                                          
               GO TO 220                                                
  202          DUM1 = DUM1*SQRT5                                        
               GO TO 220                                                
  203          DUM1 = DUM1*SQRT3                                        
               GO TO 220                                                
  204          DUM1 = CGK*BINV                                          
               GO TO 220                                                
  205          DUM1 = DUM1*SQRT7                                        
               GO TO 220                                                
  206          DUM1 = DUM1*SQRT5/SQRT3                                  
               GO TO 220                                                
  207          DUM1 = DUM1*SQRT3                                        
               GO TO 220                                                
  208          DUM1 = CHK*BINV                                          
               GO TO 220                                                
  209          DUM1 = DUM1*SQRT9                                        
               GO TO 220                                                
  210          DUM1 = DUM1*SQRT7/SQRT3                                  
               GO TO 220                                                
  211          DUM1 = DUM1*SQRT3                                        
               GO TO 220                                                
  212          DUM1 = DUM1*SQRT5/SQRT3                                  
               GO TO 220                                                
  213          DUM1 = CIK*BINV                                          
               GO TO 220                                                
  214          DUM1 = DUM1*SQRT11                                       
               GO TO 220                                                
  215          DUM1 = DUM1*SQRT3                                        
               GO TO 220                                                
  216          DUM1 = DUM1*SQRT3                                        
               GO TO 220                                                
  217          DUM1 = DUM1*SQRT7/(SQRT5*SQRT3)                          
               GO TO 220                                                
  218          DUM1 = DUM1*SQRT5                                        
               GO TO 220                                                
  219          DUM1 = DUM1*SQRT5/SQRT3                                  
C                                                                       
  220          IF (KANDL) MAX = K                                       
               DO 360 L = MINL,MAX                                      
                  GO TO (240,280,340,340,300,340,340,320,340,340,       
     1                   321,340,340,322,340,340,340,340,340,323,       
     1                   324,340,340,325,340,340,340,340,340,326,       
     1                   340,340,327,340,340,                           
     1                   328,340,340,329,340,340,340,340,340,330,       
     1                   340,340,340,340,340,331,340,340,332,340,       
     1                   340,                                           
     1                   333,340,340,334,340,340,340,340,340,335,       
     1                   340,340,340,340,340,336,340,340,337,340,       
     1                   340,338,340,340,340,340,340,339),L             
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
  321             DUM2 = DUM1*CFL                                       
                  IF (DOUBLE) DUM2 = DUM2+DUM2                          
                  GO TO 340                                             
  322             DUM2 = DUM2*SQRT5                                     
                  GO TO 340                                             
  323             DUM2 = DUM2*SQRT3                                     
                  GO TO 340                                             
  324             DUM2 = DUM1*CGL                                       
                  IF (DOUBLE) DUM2 = DUM2+DUM2                          
                  GO TO 340                                             
  325             DUM2 = DUM2*SQRT7                                     
                  GO TO 340                                             
  326             DUM2 = DUM2*SQRT5/SQRT3                               
                  GO TO 340                                             
  327             DUM2 = DUM2*SQRT3                                     
                  GO TO 340                                             
  328             DUM2 = DUM1*CHL                                       
                  IF (DOUBLE) DUM2 = DUM2+DUM2                          
                  GO TO 340                                             
  329             DUM2 = DUM2*SQRT9                                     
                  GO TO 340                                             
  330             DUM2 = DUM2*SQRT7/SQRT3                               
                  GO TO 340                                             
  331             DUM2 = DUM2*SQRT3                                     
                  GO TO 340                                             
  332             DUM2 = DUM2*SQRT5/SQRT3                               
                  GO TO 340                                             
  333             DUM2 = DUM1*CIL                                       
                  IF (DOUBLE) DUM2 = DUM2+DUM2                          
                  GO TO 340                                             
  334             DUM2 = DUM2*SQRT11                                    
                  GO TO 340                                             
  335             DUM2 = DUM2*SQRT3                                     
                  GO TO 340                                             
  336             DUM2 = DUM2*SQRT3                                     
                  GO TO 340                                             
  337             DUM2 = DUM2*SQRT7/(SQRT5*SQRT3)                       
                  GO TO 340                                             
  338             DUM2 = DUM2*SQRT5                                     
                  GO TO 340                                             
  339             DUM2 = DUM2*SQRT5/SQRT3                               
C                                                                       
  340             N = N+1                                               
                  DKL(N) = DUM2                                         
  360          CONTINUE                                                 
  370       CONTINUE                                                    
C                                                                       
C           ----- PAIR OF I,J PRIMITIVES                                
C                                                                       
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
               XX = RHO*((XA-XB)*(XA-XB) + (YA-YB)*(YA-YB)              
     *                                   + (ZA-ZB)*(ZA-ZB))             
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
C                                                                       
C              ----- ROOTS AND WEIGHTS FOR QUADRATURE                   
C                IF (NROOTS .LE. 3) THEN                                
C                    CALL RT123TS(XX,NROOTS,RT1,RT2,RT3,WW1,WW2,WW3)    
C                    U(1)=RT1                                           
C                    U(2)=RT2                                           
C                    U(3)=RT3                                           
C                    W(1)=WW1                                           
C                    W(2)=WW2                                           
C                    W(3)=WW3                                           
C                ENDIF                                                  
C                IF (NROOTS .EQ. 4) THEN                                
C                    CALL ROOT4TS(XX,NROOTS,RT1,RT2,RT3,RT4,            
C      * WW1,WW2,WW3,WW4)                                               
C                    U(1)=RT1                                           
C                    U(2)=RT2                                           
C                    U(3)=RT3                                           
C                    U(4)=RT4                                           
C                    W(1)=WW1                                           
C                    W(2)=WW2                                           
C                    W(3)=WW3                                           
C                    W(4)=WW4                                           
C                ENDIF                                                  
               !IF (NROOTS .EQ. 5) THEN                                 
                   CALL ROOT5TS(XX,NROOTS,RT1,RT2,RT3,RT4,RT5,          
     * WW1,WW2,WW3,WW4,WW5)                                             
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
               !ENDIF                                                   
               MM = 0                                                   
               MAX = NMAX+1                                             
C                                                                       
C              COMPUTE TWO-ELECTRON INTEGRALS FOR EACH ROOT             
C                                                                       
               DO 420 M = 1,NROOTS                                      
C               DO 420 M = 1,5                                          
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
                                                                        
                  CALL XYZINT_gpu(XIN,YIN,ZIN,                          
     * IN,KN,NI,NJ,NK,NL,NMAX,MMAX,                                     
     * BP01,B00,B10,XCP00,XC00,YCP00,YC00,ZCP00,ZC00,F00,               
     * DXIJ,DYIJ,DZIJ,DXKL,DYKL,DZKL)                                   
                  MM = MM+2401                                          
  420          CONTINUE                                                 
C                                                                       
C              ----- FORM (I,J//K,L) INTEGRALS OVER FUNCTIONS           
C                                                                       
               CALL FORMS_gpu_5(GHONDO,NROOTS,DKL,DIJ,XIN,YIN,ZIN,      
     * IJGT,IJX,IJY,IJZ,IK,KLGT,KLX,KLY,KLZ,IJ)                         
  440       NN = NN+49                                                  
  460    CONTINUE                                                       
  480 CONTINUE                                                          
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2A   *DECK GENRAL_gpu_4                                     
      SUBROUTINE GENRAL_gpu_4(GHONDO,DDIJ,DKL,DIJ,                      
     * IJGT,IJX,IJY,IJZ,IK,KLGT,KLX,KLY,KLZ,                            
     * AA,R,X1,Y1,Z1,IJD,IANDJ,KANDL,SAME,                              
     * LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,                             
     * qq4,mini,maxi,minj,maxj,mink,maxk,minl,maxl,                     
     * NIJ,IJ,KL,IJKL,                                                  
     * NROOTS,                                                          
     * ag,csa,cpa,cda,cfa,cga,cha,cia,                                  
     * bg,csb,cpb,cdb,cfb,cgb,chb,cib,                                  
     * cg,csc,cpc,cdc,cfc,cgc,chc,cic,                                  
     * dg,csd,cpd,cdd,cfd,cgd,chd,cid,                                  
     * XI,YI,ZI,XJ,YJ,ZJ,RRI,XK,YK,ZK,XL,YL,ZL,RRK,                     
     * NGA,NGB,NGC,NGD)                                                 
                                                                        
                                                                        
      USE lrcdft, ONLY: EMU2                                            
      use mx_limits, only: mxgsh,mxg2                                   
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      DIMENSION GHONDO(*),DDIJ(*)                                       
C                                                                       
      LOGICAL IANDJ,KANDL,SAME,OUT,NORM,DOUBLE                          
      DIMENSION :: XIN(31213),YIN(31213),ZIN(31213)                     
      DIMENSION :: DKL(784),DIJ(784)                                    
      DIMENSION :: IJGT(784),IJX(784),IJY(784),IJZ(784),IK(784)         
      DIMENSION :: KLGT(784),KLX(784),KLY(784),KLZ(784)                 
      DIMENSION :: AA(MXG2),R(MXG2),X1(MXG2),Y1(MXG2),Z1(MXG2),IJD(784) 
                                                                        
      double precision :: XX                                            
      double precision :: U(13)                                         
      double precision :: W(13)                                         
      integer :: NROOTS                                                 
      DIMENSION :: IN(13),KN(13)                                        
      INTEGER :: NI,NJ,NK,NL,NMAX,MMAX                                  
      double precision :: BP01,B00,B10,XCP00,XC00,YCP00,YC00,ZCP00,ZC00 
      double precision :: F00,DXIJ,DYIJ,DZIJ,DXKL,DYKL,DZKL             
                                                                        
      DIMENSION :: AG(MXGSH),CSA(MXGSH),CPA(MXGSH),CDA(MXGSH)           
      DIMENSION :: CFA(MXGSH),CGA(MXGSH),CHA(MXGSH),CIA(MXGSH)          
      DIMENSION :: BG(MXGSH),CSB(MXGSH),CPB(MXGSH),CDB(MXGSH)           
      DIMENSION :: CFB(MXGSH),CGB(MXGSH),CHB(MXGSH),CIB(MXGSH)          
      DIMENSION :: CG(MXGSH),CSC(MXGSH),CPC(MXGSH),CDC(MXGSH)           
      DIMENSION :: CFC(MXGSH),CGC(MXGSH),CHC(MXGSH),CIC(MXGSH)          
      DIMENSION :: DG(MXGSH),CSD(MXGSH),CPD(MXGSH),CDD(MXGSH)           
      DIMENSION :: CFD(MXGSH),CGD(MXGSH),CHD(MXGSH),CID(MXGSH)          
      double precision :: XI,YI,ZI,XJ,YJ,ZJ,RRI,XK,YK,ZK,XL,YL,ZL,RRK   
      integer :: NGA,NGB,NGC,NGD                                        
                                                                        
      double precision :: qq4                                           
      INTEGER :: LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL                    
      INTEGER :: MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL                
      INTEGER :: NIJ,IJ,KL,IJKL                                         
C                                                                       
      DIMENSION IN1(13)                                                 
C                                                                       
      PARAMETER (SQRT3=1.73205080756888D+00, SQRT5=2.23606797749979D+00,
     *           SQRT7=2.64575131106459D+00, PI252=34.986836655250D+00, 
     *           SQRT9=3.0D+00,SQRT11=3.3166247903553998D+00,           
     *           ZERO=0.0D+00, HALF=0.5D+00, ONE=1.0D+00,               
     *           TOL=46.0515999999999934 )                              
                                                                        
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
C                                                                       
C     ----- K PRIMITIVE                                                 
C                                                                       
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
         CFK = CFC(KG)*FACTOR                                           
         CGK = CGC(KG)*FACTOR                                           
         CHK = CHC(KG)*FACTOR                                           
         CIK = CIC(KG)*FACTOR                                           
C                                                                       
C        ----- L PRIMITIVE                                              
C                                                                       
         IF (KANDL) LGMAX = KG                                          
         DO 460 LG = 1,LGMAX                                            
            AL = DG(LG)                                                 
            B = AK+AL                                                   
            BINV = ONE/B                                                
            BBRRK = AL*BRRK*BINV                                        
            IF (BBRRK .GT. TOL) GO TO 460                               
            CSL = CSD(LG)                                               
            CPL = CPD(LG)                                               
            CDL = CDD(LG)                                               
            CFL = CFD(LG)                                               
            CGL = CGD(LG)                                               
            CHL = CHD(LG)                                               
            CIL = CID(LG)                                               
            XB = (AKXK+AL*XL)*BINV                                      
            YB = (AKYK+AL*YL)*BINV                                      
            ZB = (AKZK+AL*ZL)*BINV                                      
            BXBK = B*(XB-XK)                                            
            BYBK = B*(YB-YK)                                            
            BZBK = B*(ZB-ZK)                                            
            BXBI = B*(XB-XI)                                            
            BYBI = B*(YB-YI)                                            
            BZBI = B*(ZB-ZI)                                            
C                                                                       
C           ----- DENSITY FACTOR                                        
C                                                                       
            DOUBLE=KANDL.AND.KG.NE.LG                                   
            N = 0                                                       
            MAX = MAXL                                                  
            DUM1 = ZERO                                                 
            DUM2 = ZERO                                                 
            DO 370 K = MINK,MAXK                                        
               GO TO (140,160,220,220,180,220,220,200,220,220,          
     1                201,220,220,202,220,220,220,220,220,203,          
     1                204,220,220,205,220,220,220,220,220,206,          
     1                220,220,207,220,220,                              
     1                208,220,220,209,220,220,220,220,220,210,          
     1                220,220,220,220,220,211,220,220,212,220,          
     1                220,                                              
     1                213,220,220,214,220,220,220,220,220,215,          
     1                220,220,220,220,220,216,220,220,217,220,          
     1                220,218,220,220,220,220,220,219),K                
  140          DUM1 = CSK*BINV                                          
               GO TO 220                                                
  160          DUM1 = CPK*BINV                                          
               GO TO 220                                                
  180          DUM1 = CDK*BINV                                          
               GO TO 220                                                
  200          DUM1 = DUM1*SQRT3                                        
               GO TO 220                                                
  201          DUM1 = CFK*BINV                                          
               GO TO 220                                                
  202          DUM1 = DUM1*SQRT5                                        
               GO TO 220                                                
  203          DUM1 = DUM1*SQRT3                                        
               GO TO 220                                                
  204          DUM1 = CGK*BINV                                          
               GO TO 220                                                
  205          DUM1 = DUM1*SQRT7                                        
               GO TO 220                                                
  206          DUM1 = DUM1*SQRT5/SQRT3                                  
               GO TO 220                                                
  207          DUM1 = DUM1*SQRT3                                        
               GO TO 220                                                
  208          DUM1 = CHK*BINV                                          
               GO TO 220                                                
  209          DUM1 = DUM1*SQRT9                                        
               GO TO 220                                                
  210          DUM1 = DUM1*SQRT7/SQRT3                                  
               GO TO 220                                                
  211          DUM1 = DUM1*SQRT3                                        
               GO TO 220                                                
  212          DUM1 = DUM1*SQRT5/SQRT3                                  
               GO TO 220                                                
  213          DUM1 = CIK*BINV                                          
               GO TO 220                                                
  214          DUM1 = DUM1*SQRT11                                       
               GO TO 220                                                
  215          DUM1 = DUM1*SQRT3                                        
               GO TO 220                                                
  216          DUM1 = DUM1*SQRT3                                        
               GO TO 220                                                
  217          DUM1 = DUM1*SQRT7/(SQRT5*SQRT3)                          
               GO TO 220                                                
  218          DUM1 = DUM1*SQRT5                                        
               GO TO 220                                                
  219          DUM1 = DUM1*SQRT5/SQRT3                                  
C                                                                       
  220          IF (KANDL) MAX = K                                       
               DO 360 L = MINL,MAX                                      
                  GO TO (240,280,340,340,300,340,340,320,340,340,       
     1                   321,340,340,322,340,340,340,340,340,323,       
     1                   324,340,340,325,340,340,340,340,340,326,       
     1                   340,340,327,340,340,                           
     1                   328,340,340,329,340,340,340,340,340,330,       
     1                   340,340,340,340,340,331,340,340,332,340,       
     1                   340,                                           
     1                   333,340,340,334,340,340,340,340,340,335,       
     1                   340,340,340,340,340,336,340,340,337,340,       
     1                   340,338,340,340,340,340,340,339),L             
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
  321             DUM2 = DUM1*CFL                                       
                  IF (DOUBLE) DUM2 = DUM2+DUM2                          
                  GO TO 340                                             
  322             DUM2 = DUM2*SQRT5                                     
                  GO TO 340                                             
  323             DUM2 = DUM2*SQRT3                                     
                  GO TO 340                                             
  324             DUM2 = DUM1*CGL                                       
                  IF (DOUBLE) DUM2 = DUM2+DUM2                          
                  GO TO 340                                             
  325             DUM2 = DUM2*SQRT7                                     
                  GO TO 340                                             
  326             DUM2 = DUM2*SQRT5/SQRT3                               
                  GO TO 340                                             
  327             DUM2 = DUM2*SQRT3                                     
                  GO TO 340                                             
  328             DUM2 = DUM1*CHL                                       
                  IF (DOUBLE) DUM2 = DUM2+DUM2                          
                  GO TO 340                                             
  329             DUM2 = DUM2*SQRT9                                     
                  GO TO 340                                             
  330             DUM2 = DUM2*SQRT7/SQRT3                               
                  GO TO 340                                             
  331             DUM2 = DUM2*SQRT3                                     
                  GO TO 340                                             
  332             DUM2 = DUM2*SQRT5/SQRT3                               
                  GO TO 340                                             
  333             DUM2 = DUM1*CIL                                       
                  IF (DOUBLE) DUM2 = DUM2+DUM2                          
                  GO TO 340                                             
  334             DUM2 = DUM2*SQRT11                                    
                  GO TO 340                                             
  335             DUM2 = DUM2*SQRT3                                     
                  GO TO 340                                             
  336             DUM2 = DUM2*SQRT3                                     
                  GO TO 340                                             
  337             DUM2 = DUM2*SQRT7/(SQRT5*SQRT3)                       
                  GO TO 340                                             
  338             DUM2 = DUM2*SQRT5                                     
                  GO TO 340                                             
  339             DUM2 = DUM2*SQRT5/SQRT3                               
C                                                                       
  340             N = N+1                                               
                  DKL(N) = DUM2                                         
  360          CONTINUE                                                 
  370       CONTINUE                                                    
C                                                                       
C           ----- PAIR OF I,J PRIMITIVES                                
C                                                                       
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
               XX = RHO*((XA-XB)*(XA-XB) + (YA-YB)*(YA-YB)              
     *                                   + (ZA-ZB)*(ZA-ZB))             
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
C                                                                       
C              ----- ROOTS AND WEIGHTS FOR QUADRATURE                   
C                IF (NROOTS .LE. 3) THEN                                
C                    CALL RT123TS(XX,NROOTS,RT1,RT2,RT3,WW1,WW2,WW3)    
C                    U(1)=RT1                                           
C                    U(2)=RT2                                           
C                    U(3)=RT3                                           
C                    W(1)=WW1                                           
C                    W(2)=WW2                                           
C                    W(3)=WW3                                           
C                ENDIF                                                  
C               IF (NROOTS .EQ. 4) THEN                                 
                   CALL ROOT4TS(XX,NROOTS,RT1,RT2,RT3,RT4,              
     * WW1,WW2,WW3,WW4)                                                 
                   U(1)=RT1                                             
                   U(2)=RT2                                             
                   U(3)=RT3                                             
                   U(4)=RT4                                             
                   W(1)=WW1                                             
                   W(2)=WW2                                             
                   W(3)=WW3                                             
                   W(4)=WW4                                             
C               ENDIF                                                   
C                IF (NROOTS .EQ. 5) THEN                                
C                    CALL ROOT5TS(XX,NROOTS,RT1,RT2,RT3,RT4,RT5,        
C      * WW1,WW2,WW3,WW4,WW5)                                           
C                    U(1)=RT1                                           
C                    U(2)=RT2                                           
C                    U(3)=RT3                                           
C                    U(4)=RT4                                           
C                    U(5)=RT5                                           
C                    W(1)=WW1                                           
C                    W(2)=WW2                                           
C                    W(3)=WW3                                           
C                    W(4)=WW4                                           
C                    W(5)=WW5                                           
C                ENDIF                                                  
               MM = 0                                                   
               MAX = NMAX+1                                             
C                                                                       
C              COMPUTE TWO-ELECTRON INTEGRALS FOR EACH ROOT             
C                                                                       
               DO 420 M = 1,NROOTS                                      
C               DO 420 M = 1,5                                          
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
                                                                        
                  CALL XYZINT_gpu(XIN,YIN,ZIN,                          
     * IN,KN,NI,NJ,NK,NL,NMAX,MMAX,                                     
     * BP01,B00,B10,XCP00,XC00,YCP00,YC00,ZCP00,ZC00,F00,               
     * DXIJ,DYIJ,DZIJ,DXKL,DYKL,DZKL)                                   
                  MM = MM+2401                                          
  420          CONTINUE                                                 
C                                                                       
C              ----- FORM (I,J//K,L) INTEGRALS OVER FUNCTIONS           
C                                                                       
C               write(*,*) "nroots before forms4", NROOTS               
               CALL FORMS_gpu_4(GHONDO,NROOTS,DKL,DIJ,XIN,YIN,ZIN,      
     * IJGT,IJX,IJY,IJZ,IK,KLGT,KLX,KLY,KLZ,IJ)                         
  440       NN = NN+49                                                  
  460    CONTINUE                                                       
  480 CONTINUE                                                          
C                                                                       
      RETURN                                                            
      END                                                               
C C*MODULE INT2A   *DECK GENRAL_gpu_2                                   
C       SUBROUTINE GENRAL_gpu_2(GHONDO,DDIJ,DKL,DIJ,                    
C      * IJGT,IJX,IJY,IJZ,IK,KLGT,KLX,KLY,KLZ,                          
C      * AA,R,X1,Y1,Z1,IJD,IANDJ,KANDL,SAME,                            
C      * LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,                           
C      * qq4,mini,maxi,minj,maxj,mink,maxk,minl,maxl,                   
C      * NIJ,IJ,KL,IJKL,                                                
C      * NROOTS,                                                        
C      * ag,csa,cpa,cda,cfa,cga,cha,cia,                                
C      * bg,csb,cpb,cdb,cfb,cgb,chb,cib,                                
C      * cg,csc,cpc,cdc,cfc,cgc,chc,cic,                                
C      * dg,csd,cpd,cdd,cfd,cgd,chd,cid,                                
C      * XI,YI,ZI,XJ,YJ,ZJ,RRI,XK,YK,ZK,XL,YL,ZL,RRK,                   
C      * NGA,NGB,NGC,NGD)                                               
                                                                        
                                                                        
C       USE lrcdft, ONLY: EMU2                                          
C       use mx_limits, only: mxgsh,mxg2                                 
C C                                                                     
C       IMPLICIT DOUBLE PRECISION (A-H,O-Z)                             
C C                                                                     
C       DIMENSION GHONDO(*),DDIJ(*)                                     
C C                                                                     
C       LOGICAL IANDJ,KANDL,SAME,OUT,NORM,DOUBLE                        
C       DIMENSION :: XIN(31213),YIN(31213),ZIN(31213)                   
C       DIMENSION :: DKL(784),DIJ(784)                                  
C       DIMENSION :: IJGT(784),IJX(784),IJY(784),IJZ(784),IK(784)       
C       DIMENSION :: KLGT(784),KLX(784),KLY(784),KLZ(784)               
C       DIMENSION :: AA(MXG2),R(MXG2),X1(MXG2),Y1(MXG2),Z1(MXG2),IJD(784
                                                                        
C       double precision :: XX                                          
C       double precision :: U(13)                                       
C       double precision :: W(13)                                       
C       integer :: NROOTS                                               
C       DIMENSION :: IN(13),KN(13)                                      
C       INTEGER :: NI,NJ,NK,NL,NMAX,MMAX                                
C       double precision :: BP01,B00,B10,XCP00,XC00,YCP00,YC00,ZCP00,ZC0
C       double precision :: F00,DXIJ,DYIJ,DZIJ,DXKL,DYKL,DZKL           
                                                                        
C       DIMENSION :: AG(MXGSH),CSA(MXGSH),CPA(MXGSH),CDA(MXGSH)         
C       DIMENSION :: CFA(MXGSH),CGA(MXGSH),CHA(MXGSH),CIA(MXGSH)        
C       DIMENSION :: BG(MXGSH),CSB(MXGSH),CPB(MXGSH),CDB(MXGSH)         
C       DIMENSION :: CFB(MXGSH),CGB(MXGSH),CHB(MXGSH),CIB(MXGSH)        
C       DIMENSION :: CG(MXGSH),CSC(MXGSH),CPC(MXGSH),CDC(MXGSH)         
C       DIMENSION :: CFC(MXGSH),CGC(MXGSH),CHC(MXGSH),CIC(MXGSH)        
C       DIMENSION :: DG(MXGSH),CSD(MXGSH),CPD(MXGSH),CDD(MXGSH)         
C       DIMENSION :: CFD(MXGSH),CGD(MXGSH),CHD(MXGSH),CID(MXGSH)        
C       double precision :: XI,YI,ZI,XJ,YJ,ZJ,RRI,XK,YK,ZK,XL,YL,ZL,RRK 
C       integer :: NGA,NGB,NGC,NGD                                      
                                                                        
C       double precision :: qq4                                         
C       INTEGER :: LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL                  
C       INTEGER :: MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL              
C       INTEGER :: NIJ,IJ,KL,IJKL                                       
C C                                                                     
C       DIMENSION IN1(13)                                               
C C                                                                     
C       PARAMETER (SQRT3=1.73205080756888D+00, SQRT5=2.23606797749979D+0
C      *           SQRT7=2.64575131106459D+00, PI252=34.986836655250D+00
C      *           SQRT9=3.0D+00,SQRT11=3.3166247903553998D+00,         
C      *           ZERO=0.0D+00, HALF=0.5D+00, ONE=1.0D+00,             
C      *           TOL=46.0515999999999934 )                            
                                                                        
C       FACTOR = PI252*QQ4                                              
C       NI = LIT-1                                                      
C       NJ = LJT-1                                                      
C       NK = LKT-1                                                      
C       NL = LLT-1                                                      
C       DXIJ = XI-XJ                                                    
C       DYIJ = YI-YJ                                                    
C       DZIJ = ZI-ZJ                                                    
C       DXKL = XK-XL                                                    
C       DYKL = YK-YL                                                    
C       DZKL = ZK-ZL                                                    
C       NMAX = NI+NJ                                                    
C       MMAX = NK+NL                                                    
C       MAX = NMAX+1                                                    
C       DO I = 1,MAX                                                    
C          N = I-1                                                      
C          IF (N .LE. NI) IN1(I) = 343*N+1                              
C          IF (N .GT. NI) IN1(I) = 343*NI+49*(N-NI)+1                   
C       ENDDO                                                           
C       MAX = MMAX+1                                                    
C       DO K = 1,MAX                                                    
C          N = K-1                                                      
C          IF (N .LE. NK) KN(K) = 7*N                                   
C          IF (N .GT. NK) KN(K) = 7*NK+N-NK                             
C       ENDDO                                                           
C C                                                                     
C C     ----- K PRIMITIVE                                               
C C                                                                     
C       LGMAX = NGD                                                     
C       DO 480 KG = 1,NGC                                               
C          AK = CG(KG)                                                  
C          BRRK = AK*RRK                                                
C          AKXK = AK*XK                                                 
C          AKYK = AK*YK                                                 
C          AKZK = AK*ZK                                                 
C          CSK = CSC(KG)*FACTOR                                         
C C          write(*,*) "CSK in genral", CSK                            
C          CPK = CPC(KG)*FACTOR                                         
C C          write(*,*) "CPK in genral", CPK                            
C          CDK = CDC(KG)*FACTOR                                         
C C         write(*,*) "CDK in genral", CDK                             
C C          CFK = CFC(KG)*FACTOR                                       
C C          CGK = CGC(KG)*FACTOR                                       
C C          CHK = CHC(KG)*FACTOR                                       
C C          CIK = CIC(KG)*FACTOR                                       
C C                                                                     
C C        ----- L PRIMITIVE                                            
C C                                                                     
C          IF (KANDL) LGMAX = KG                                        
C          DO 460 LG = 1,LGMAX                                          
C             AL = DG(LG)                                               
C             B = AK+AL                                                 
C             BINV = ONE/B                                              
C             BBRRK = AL*BRRK*BINV                                      
C             CSL = CSD(LG)                                             
C             CPL = CPD(LG)                                             
C             CDL = CDD(LG)                                             
C             CFL = CFD(LG)                                             
C             CGL = CGD(LG)                                             
C             CHL = CHD(LG)                                             
C             CIL = CID(LG)                                             
C             XB = (AKXK+AL*XL)*BINV                                    
C             YB = (AKYK+AL*YL)*BINV                                    
C             ZB = (AKZK+AL*ZL)*BINV                                    
C             BXBK = B*(XB-XK)                                          
C             BYBK = B*(YB-YK)                                          
C             BZBK = B*(ZB-ZK)                                          
C             BXBI = B*(XB-XI)                                          
C             BYBI = B*(YB-YI)                                          
C             BZBI = B*(ZB-ZI)                                          
C C                                                                     
C C           ----- DENSITY FACTOR                                      
C C                                                                     
C             DOUBLE=KANDL.AND.KG.NE.LG                                 
C             N = 0                                                     
C             MAX = MAXL                                                
C             DUM1 = ZERO                                               
C             DUM2 = ZERO                                               
C             DO 370 K = MINK,MAXK                                      
C                GO TO (140,160,220,220,180,220,220,200,220,220,        
C      1                201),K                                          
                                                                        
C   140          DUM1 = CSK*BINV                                        
C                GO TO 220                                              
C   160          DUM1 = CPK*BINV                                        
C                GO TO 220                                              
C   180          DUM1 = CDK*BINV                                        
C                GO TO 220                                              
C   200          DUM1 = DUM1*SQRT3                                      
C                GO TO 220                                              
C   201          GO TO 220                                              
C C   202          DUM1 = DUM1*SQRT5                                    
C C                GO TO 220                                            
C C   203          DUM1 = DUM1*SQRT3                                    
C C                GO TO 220                                            
C C   204          GO TO 220                                            
C C   205          DUM1 = DUM1*SQRT7                                    
C C                GO TO 220                                            
C C   206          DUM1 = DUM1*SQRT5/SQRT3                              
C C                GO TO 220                                            
C C   207          DUM1 = DUM1*SQRT3                                    
C C                GO TO 220                                            
C C   208          GO TO 220                                            
C C   209          DUM1 = DUM1*SQRT9                                    
C C                GO TO 220                                            
C C   210          DUM1 = DUM1*SQRT7/SQRT3                              
C C                GO TO 220                                            
C C   211          DUM1 = DUM1*SQRT3                                    
C C                GO TO 220                                            
C C   212          DUM1 = DUM1*SQRT5/SQRT3                              
C C                GO TO 220                                            
C C   213          GO TO 220                                            
C C   214          DUM1 = DUM1*SQRT11                                   
C C                GO TO 220                                            
C C   215          DUM1 = DUM1*SQRT3                                    
C C                GO TO 220                                            
C C   216          DUM1 = DUM1*SQRT3                                    
C C                GO TO 220                                            
C C   217          DUM1 = DUM1*SQRT7/(SQRT5*SQRT3)                      
C C                GO TO 220                                            
C C   218          DUM1 = DUM1*SQRT5                                    
C C                GO TO 220                                            
C C   219          DUM1 = DUM1*SQRT5/SQRT3                              
C C                                                                     
C   220          IF (KANDL) MAX = K                                     
C                DO 360 L = MINL,MAX                                    
C                   GO TO (240,280,340,340,300,340,340,320,340,340,     
C      1                   321,340,340,322,340,340,340,340,340,323,     
C      1                   324,340,340,325,340,340,340,340,340,326,     
C      1                   340,340,327,340,340,                         
C      1                   328,340,340,329,340,340,340,340,340,330,     
C      1                   340,340,340,340,340,331,340,340,332,340,     
C      1                   340,                                         
C      1                   333,340,340,334,340,340,340,340,340,335,     
C      1                   340,340,340,340,340,336,340,340,337,340,     
C      1                   340,338,340,340,340,340,340,339),L           
C   240             DUM2 = DUM1*CSL                                     
C                   IF ( .NOT. DOUBLE) GO TO 340                        
C                   IF (K .GT. 1) GO TO 260                             
C                   DUM2 = DUM2+DUM2                                    
C                   GO TO 340                                           
C   260             DUM2 = DUM2+CSK*CPL*BINV                            
C                   GO TO 340                                           
C   280             DUM2 = DUM1*CPL                                     
C                   IF (DOUBLE) DUM2 = DUM2+DUM2                        
C                   GO TO 340                                           
C   300             DUM2 = DUM1*CDL                                     
C                   IF (DOUBLE) DUM2 = DUM2+DUM2                        
C                   GO TO 340                                           
C   320             DUM2 = DUM2*SQRT3                                   
C                   GO TO 340                                           
C   321             DUM2 = DUM1*CFL                                     
C                   IF (DOUBLE) DUM2 = DUM2+DUM2                        
C                   GO TO 340                                           
C   322             DUM2 = DUM2*SQRT5                                   
C                   GO TO 340                                           
C   323             DUM2 = DUM2*SQRT3                                   
C                   GO TO 340                                           
C   324             DUM2 = DUM1*CGL                                     
C                   IF (DOUBLE) DUM2 = DUM2+DUM2                        
C                   GO TO 340                                           
C   325             DUM2 = DUM2*SQRT7                                   
C                   GO TO 340                                           
C   326             DUM2 = DUM2*SQRT5/SQRT3                             
C                   GO TO 340                                           
C   327             DUM2 = DUM2*SQRT3                                   
C                   GO TO 340                                           
C   328             DUM2 = DUM1*CHL                                     
C                   IF (DOUBLE) DUM2 = DUM2+DUM2                        
C                   GO TO 340                                           
C   329             DUM2 = DUM2*SQRT9                                   
C                   GO TO 340                                           
C   330             DUM2 = DUM2*SQRT7/SQRT3                             
C                   GO TO 340                                           
C   331             DUM2 = DUM2*SQRT3                                   
C                   GO TO 340                                           
C   332             DUM2 = DUM2*SQRT5/SQRT3                             
C                   GO TO 340                                           
C   333             DUM2 = DUM1*CIL                                     
C                   IF (DOUBLE) DUM2 = DUM2+DUM2                        
C                   GO TO 340                                           
C   334             DUM2 = DUM2*SQRT11                                  
C                   GO TO 340                                           
C   335             DUM2 = DUM2*SQRT3                                   
C                   GO TO 340                                           
C   336             DUM2 = DUM2*SQRT3                                   
C                   GO TO 340                                           
C   337             DUM2 = DUM2*SQRT7/(SQRT5*SQRT3)                     
C                   GO TO 340                                           
C   338             DUM2 = DUM2*SQRT5                                   
C                   GO TO 340                                           
C   339             DUM2 = DUM2*SQRT5/SQRT3                             
C C                                                                     
C   340             N = N+1                                             
C                   DKL(N) = DUM2                                       
C   360          CONTINUE                                               
C   370       CONTINUE                                                  
C C                                                                     
C C           ----- PAIR OF I,J PRIMITIVES                              
C C                                                                     
C             NN = 0                                                    
C             DO 440 N = 1,NIJ                                          
C                DUM = BBRRK+R(N)                                       
C                DO I = 1,IJ                                            
C                   DIJ(I) = DDIJ(IJD(I)+NN)                            
C                ENDDO                                                  
C                A = AA(N)                                              
C                AB = A*B                                               
C                AANDB = A+B                                            
C                EXPE = EXP(-DUM)/SQRT(AANDB)                           
C                RHO = AB/AANDB                                         
C                XA = X1(N)                                             
C                YA = Y1(N)                                             
C                ZA = Z1(N)                                             
C                XX = RHO*((XA-XB)*(XA-XB) + (YA-YB)*(YA-YB)            
C      *                                   + (ZA-ZB)*(ZA-ZB))           
C                AXAK = A*(XA-XK)                                       
C                AYAK = A*(YA-YK)                                       
C                AZAK = A*(ZA-ZK)                                       
C                AXAI = A*(XA-XI)                                       
C                AYAI = A*(YA-YI)                                       
C                AZAI = A*(ZA-ZI)                                       
C                C1X = BXBK+AXAK                                        
C                C2X = A*BXBK                                           
C                C3X = BXBI+AXAI                                        
C                C4X = B*AXAI                                           
C                C1Y = BYBK+AYAK                                        
C                C2Y = A*BYBK                                           
C                C3Y = BYBI+AYAI                                        
C                C4Y = B*AYAI                                           
C                C1Z = BZBK+AZAK                                        
C                C2Z = A*BZBK                                           
C                C3Z = BZBI+AZAI                                        
C                C4Z = B*AZAI                                           
C C                                                                     
C C              ----- ROOTS AND WEIGHTS FOR QUADRATURE                 
C                CALL RT123TS(XX,NROOTS,RT1,RT2,RT3,WW1,WW2,WW3)        
C                U(1)=RT1                                               
C                U(2)=RT2                                               
C                U(3)=RT3                                               
C                W(1)=WW1                                               
C                W(2)=WW2                                               
C                W(3)=WW3                                               
                                                                        
C                MM = 0                                                 
C                MAX = NMAX+1                                           
C C                                                                     
C C              COMPUTE TWO-ELECTRON INTEGRALS FOR EACH ROOT           
C C                                                                     
C                DO M = 1,2                                             
C                   U2 = U(M)*RHO                                       
C                   F00 = EXPE*W(M)                                     
C                   DO I = 1,MAX                                        
C                      IN(I) = IN1(I)+MM                                
C                   ENDDO                                               
                                                                        
C                      DUMINV = ONE/(AB+U2*AANDB)                       
C                      DM2INV = HALF*DUMINV                             
C                      BP01 = (A+U2)*DM2INV                             
C                      B00 = U2*DM2INV                                  
C                      B10 = (B+U2)*DM2INV                              
C                      XCP00 = (U2*C1X+C2X)*DUMINV                      
C                      XC00 = (U2*C3X+C4X)*DUMINV                       
C                      YCP00 = (U2*C1Y+C2Y)*DUMINV                      
C                      YC00 = (U2*C3Y+C4Y)*DUMINV                       
C                      ZCP00 = (U2*C1Z+C2Z)*DUMINV                      
C                      ZC00 = (U2*C3Z+C4Z)*DUMINV                       
                                                                        
C                   CALL XYZINT_gpu_2(XIN,YIN,ZIN,                      
C      * IN,KN,NI,NJ,NK,NL,NMAX,MMAX,                                   
C      * BP01,B00,B10,XCP00,XC00,YCP00,YC00,ZCP00,ZC00,F00,             
C      * DXIJ,DYIJ,DZIJ,DXKL,DYKL,DZKL)                                 
C                   MM = MM+2401                                        
C                ENDDO                                                  
C C                                                                     
C C              ----- FORM (I,J//K,L) INTEGRALS OVER FUNCTIONS         
C C                                                                     
C       DO I = 1,IJ                                                     
C       D1 = DIJ(I)                                                     
C       NX = IJX(I)                                                     
C       NY = IJY(I)                                                     
C       NZ = IJZ(I)                                                     
C       N1 = IJGT(I)                                                    
C       MAX = IK(I)                                                     
C       DO K = 1,MAX                                                    
C       MX = NX+KLX(K)                                                  
C       MY = NY+KLY(K)                                                  
C       MZ = NZ+KLZ(K)                                                  
C       NNN = N1+KLGT(K)                                                
C       GHONDO(NNN) = GHONDO(NNN) + D1*DKL(K)*                          
C      *            ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )         
C      *          +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401))        
C       ENDDO                                                           
C       ENDDO                                                           
C C      write(*,*) "GHONDO is ", GHONDO(NNN)                           
                                                                        
C   440       NN = NN+49                                                
C   460    CONTINUE                                                     
C   480 CONTINUE                                                        
C C                                                                     
C       RETURN                                                          
C       END                                                             
C*MODULE INT2A   *DECK GENRAL_gpu_3                                     
      SUBROUTINE GENRAL_gpu_3(GHONDO,DDIJ,DKL,DIJ,                      
     * IJGT,IJX,IJY,IJZ,IK,KLGT,KLX,KLY,KLZ,                            
     * AA,R,X1,Y1,Z1,IJD,IANDJ,KANDL,SAME,                              
     * LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,                             
     * qq4,mini,maxi,minj,maxj,mink,maxk,minl,maxl,                     
     * NIJ,IJ,KL,IJKL,                                                  
     * NROOTS,                                                          
     * ag,csa,cpa,cda,cfa,cga,cha,cia,                                  
     * bg,csb,cpb,cdb,cfb,cgb,chb,cib,                                  
     * cg,csc,cpc,cdc,cfc,cgc,chc,cic,                                  
     * dg,csd,cpd,cdd,cfd,cgd,chd,cid,                                  
     * XI,YI,ZI,XJ,YJ,ZJ,RRI,XK,YK,ZK,XL,YL,ZL,RRK,                     
     * NGA,NGB,NGC,NGD)                                                 
                                                                        
                                                                        
      USE lrcdft, ONLY: EMU2                                            
      use mx_limits, only: mxgsh,mxg2                                   
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      DIMENSION GHONDO(*),DDIJ(*)                                       
C                                                                       
      LOGICAL IANDJ,KANDL,SAME,OUT,NORM,DOUBLE                          
      DIMENSION :: XIN(31213),YIN(31213),ZIN(31213)                     
      DIMENSION :: DKL(784),DIJ(784)                                    
      DIMENSION :: IJGT(784),IJX(784),IJY(784),IJZ(784),IK(784)         
      DIMENSION :: KLGT(784),KLX(784),KLY(784),KLZ(784)                 
      DIMENSION :: AA(MXG2),R(MXG2),X1(MXG2),Y1(MXG2),Z1(MXG2),IJD(784) 
                                                                        
      double precision :: XX                                            
      double precision :: U(13)                                         
      double precision :: W(13)                                         
      integer :: NROOTS                                                 
      DIMENSION :: IN(13),KN(13)                                        
      INTEGER :: NI,NJ,NK,NL,NMAX,MMAX                                  
      double precision :: BP01,B00,B10,XCP00,XC00,YCP00,YC00,ZCP00,ZC00 
      double precision :: F00,DXIJ,DYIJ,DZIJ,DXKL,DYKL,DZKL             
                                                                        
      DIMENSION :: AG(MXGSH),CSA(MXGSH),CPA(MXGSH),CDA(MXGSH)           
      DIMENSION :: CFA(MXGSH),CGA(MXGSH),CHA(MXGSH),CIA(MXGSH)          
      DIMENSION :: BG(MXGSH),CSB(MXGSH),CPB(MXGSH),CDB(MXGSH)           
      DIMENSION :: CFB(MXGSH),CGB(MXGSH),CHB(MXGSH),CIB(MXGSH)          
      DIMENSION :: CG(MXGSH),CSC(MXGSH),CPC(MXGSH),CDC(MXGSH)           
      DIMENSION :: CFC(MXGSH),CGC(MXGSH),CHC(MXGSH),CIC(MXGSH)          
      DIMENSION :: DG(MXGSH),CSD(MXGSH),CPD(MXGSH),CDD(MXGSH)           
      DIMENSION :: CFD(MXGSH),CGD(MXGSH),CHD(MXGSH),CID(MXGSH)          
      double precision :: XI,YI,ZI,XJ,YJ,ZJ,RRI,XK,YK,ZK,XL,YL,ZL,RRK   
      integer :: NGA,NGB,NGC,NGD                                        
                                                                        
      double precision :: qq4                                           
      INTEGER :: LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL                    
      INTEGER :: MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL                
      INTEGER :: NIJ,IJ,KL,IJKL                                         
C                                                                       
      DIMENSION IN1(13)                                                 
C                                                                       
      PARAMETER (SQRT3=1.73205080756888D+00, SQRT5=2.23606797749979D+00,
     *           SQRT7=2.64575131106459D+00, PI252=34.986836655250D+00, 
     *           SQRT9=3.0D+00,SQRT11=3.3166247903553998D+00,           
     *           ZERO=0.0D+00, HALF=0.5D+00, ONE=1.0D+00,               
     *           TOL=46.0515999999999934 )                              
                                                                        
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
C                                                                       
C     ----- K PRIMITIVE                                                 
C                                                                       
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
         CFK = CFC(KG)*FACTOR                                           
         CGK = CGC(KG)*FACTOR                                           
         CHK = CHC(KG)*FACTOR                                           
         CIK = CIC(KG)*FACTOR                                           
C                                                                       
C        ----- L PRIMITIVE                                              
C                                                                       
         IF (KANDL) LGMAX = KG                                          
         DO 460 LG = 1,LGMAX                                            
            AL = DG(LG)                                                 
            B = AK+AL                                                   
            BINV = ONE/B                                                
            BBRRK = AL*BRRK*BINV                                        
C            IF (BBRRK .GT. TOL) GO TO 460                              
            CSL = CSD(LG)                                               
            CPL = CPD(LG)                                               
            CDL = CDD(LG)                                               
            CFL = CFD(LG)                                               
            CGL = CGD(LG)                                               
            CHL = CHD(LG)                                               
            CIL = CID(LG)                                               
            XB = (AKXK+AL*XL)*BINV                                      
            YB = (AKYK+AL*YL)*BINV                                      
            ZB = (AKZK+AL*ZL)*BINV                                      
            BXBK = B*(XB-XK)                                            
            BYBK = B*(YB-YK)                                            
            BZBK = B*(ZB-ZK)                                            
            BXBI = B*(XB-XI)                                            
            BYBI = B*(YB-YI)                                            
            BZBI = B*(ZB-ZI)                                            
C                                                                       
C           ----- DENSITY FACTOR                                        
C                                                                       
            DOUBLE=KANDL.AND.KG.NE.LG                                   
            N = 0                                                       
            MAX = MAXL                                                  
            DUM1 = ZERO                                                 
            DUM2 = ZERO                                                 
            DO 370 K = MINK,MAXK                                        
               GO TO (140,160,220,220,180,220,220,200,220,220,          
     1                201,220,220,202,220,220,220,220,220,203,          
     1                204,220,220,205,220,220,220,220,220,206,          
     1                220,220,207,220,220,                              
     1                208,220,220,209,220,220,220,220,220,210,          
     1                220,220,220,220,220,211,220,220,212,220,          
     1                220,                                              
     1                213,220,220,214,220,220,220,220,220,215,          
     1                220,220,220,220,220,216,220,220,217,220,          
     1                220,218,220,220,220,220,220,219),K                
  140          DUM1 = CSK*BINV                                          
               GO TO 220                                                
  160          DUM1 = CPK*BINV                                          
               GO TO 220                                                
  180          DUM1 = CDK*BINV                                          
               GO TO 220                                                
  200          DUM1 = DUM1*SQRT3                                        
               GO TO 220                                                
  201          DUM1 = CFK*BINV                                          
               GO TO 220                                                
  202          DUM1 = DUM1*SQRT5                                        
               GO TO 220                                                
  203          DUM1 = DUM1*SQRT3                                        
               GO TO 220                                                
  204          DUM1 = CGK*BINV                                          
               GO TO 220                                                
  205          DUM1 = DUM1*SQRT7                                        
               GO TO 220                                                
  206          DUM1 = DUM1*SQRT5/SQRT3                                  
               GO TO 220                                                
  207          DUM1 = DUM1*SQRT3                                        
               GO TO 220                                                
  208          DUM1 = CHK*BINV                                          
               GO TO 220                                                
  209          DUM1 = DUM1*SQRT9                                        
               GO TO 220                                                
  210          DUM1 = DUM1*SQRT7/SQRT3                                  
               GO TO 220                                                
  211          DUM1 = DUM1*SQRT3                                        
               GO TO 220                                                
  212          DUM1 = DUM1*SQRT5/SQRT3                                  
               GO TO 220                                                
  213          DUM1 = CIK*BINV                                          
               GO TO 220                                                
  214          DUM1 = DUM1*SQRT11                                       
               GO TO 220                                                
  215          DUM1 = DUM1*SQRT3                                        
               GO TO 220                                                
  216          DUM1 = DUM1*SQRT3                                        
               GO TO 220                                                
  217          DUM1 = DUM1*SQRT7/(SQRT5*SQRT3)                          
               GO TO 220                                                
  218          DUM1 = DUM1*SQRT5                                        
               GO TO 220                                                
  219          DUM1 = DUM1*SQRT5/SQRT3                                  
C                                                                       
  220          IF (KANDL) MAX = K                                       
               DO 360 L = MINL,MAX                                      
                  GO TO (240,280,340,340,300,340,340,320,340,340,       
     1                   321,340,340,322,340,340,340,340,340,323,       
     1                   324,340,340,325,340,340,340,340,340,326,       
     1                   340,340,327,340,340,                           
     1                   328,340,340,329,340,340,340,340,340,330,       
     1                   340,340,340,340,340,331,340,340,332,340,       
     1                   340,                                           
     1                   333,340,340,334,340,340,340,340,340,335,       
     1                   340,340,340,340,340,336,340,340,337,340,       
     1                   340,338,340,340,340,340,340,339),L             
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
  321             DUM2 = DUM1*CFL                                       
                  IF (DOUBLE) DUM2 = DUM2+DUM2                          
                  GO TO 340                                             
  322             DUM2 = DUM2*SQRT5                                     
                  GO TO 340                                             
  323             DUM2 = DUM2*SQRT3                                     
                  GO TO 340                                             
  324             DUM2 = DUM1*CGL                                       
                  IF (DOUBLE) DUM2 = DUM2+DUM2                          
                  GO TO 340                                             
  325             DUM2 = DUM2*SQRT7                                     
                  GO TO 340                                             
  326             DUM2 = DUM2*SQRT5/SQRT3                               
                  GO TO 340                                             
  327             DUM2 = DUM2*SQRT3                                     
                  GO TO 340                                             
  328             DUM2 = DUM1*CHL                                       
                  IF (DOUBLE) DUM2 = DUM2+DUM2                          
                  GO TO 340                                             
  329             DUM2 = DUM2*SQRT9                                     
                  GO TO 340                                             
  330             DUM2 = DUM2*SQRT7/SQRT3                               
                  GO TO 340                                             
  331             DUM2 = DUM2*SQRT3                                     
                  GO TO 340                                             
  332             DUM2 = DUM2*SQRT5/SQRT3                               
                  GO TO 340                                             
  333             DUM2 = DUM1*CIL                                       
                  IF (DOUBLE) DUM2 = DUM2+DUM2                          
                  GO TO 340                                             
  334             DUM2 = DUM2*SQRT11                                    
                  GO TO 340                                             
  335             DUM2 = DUM2*SQRT3                                     
                  GO TO 340                                             
  336             DUM2 = DUM2*SQRT3                                     
                  GO TO 340                                             
  337             DUM2 = DUM2*SQRT7/(SQRT5*SQRT3)                       
                  GO TO 340                                             
  338             DUM2 = DUM2*SQRT5                                     
                  GO TO 340                                             
  339             DUM2 = DUM2*SQRT5/SQRT3                               
C                                                                       
  340             N = N+1                                               
                  DKL(N) = DUM2                                         
  360          CONTINUE                                                 
  370       CONTINUE                                                    
C                                                                       
C           ----- PAIR OF I,J PRIMITIVES                                
C                                                                       
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
               XX = RHO*((XA-XB)*(XA-XB) + (YA-YB)*(YA-YB)              
     *                                   + (ZA-ZB)*(ZA-ZB))             
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
C                                                                       
               CALL RT123TS(XX,NROOTS,RT1,RT2,RT3,WW1,WW2,WW3)          
               U(1)=RT1                                                 
               U(2)=RT2                                                 
               U(3)=RT3                                                 
               W(1)=WW1                                                 
               W(2)=WW2                                                 
               W(3)=WW3                                                 
                                                                        
                                                                        
               MM = 0                                                   
               MAX = NMAX+1                                             
C                                                                       
C              COMPUTE TWO-ELECTRON INTEGRALS FOR EACH ROOT             
C                                                                       
               DO 420 M = 1,3                                           
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
                                                                        
                  CALL XYZINT_gpu(XIN,YIN,ZIN,                          
     * IN,KN,NI,NJ,NK,NL,NMAX,MMAX,                                     
     * BP01,B00,B10,XCP00,XC00,YCP00,YC00,ZCP00,ZC00,F00,               
     * DXIJ,DYIJ,DZIJ,DXKL,DYKL,DZKL)                                   
                  MM = MM+2401                                          
  420          CONTINUE                                                 
C                                                                       
C              ----- FORM (I,J//K,L) INTEGRALS OVER FUNCTIONS           
C                                                                       
                                                                        
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
      GHONDO(NNN) = GHONDO(NNN) + D1*DKL(K)*                            
     *            ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )           
     *          +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401)           
     *          +   XIN(MX+ 4802)*YIN(MY+ 4802)*ZIN(MZ+ 4802))          
      ENDDO                                                             
      ENDDO                                                             
                                                                        
C                CALL FORMS_gpu_3(GHONDO,NROOTS,DKL,DIJ,XIN,YIN,ZIN,    
C      * IJGT,IJX,IJY,IJZ,IK,KLGT,KLX,KLY,KLZ,IJ)                       
  440       NN = NN+49                                                  
  460    CONTINUE                                                       
  480 CONTINUE                                                          
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2A   *DECK SHELLS_gpu                                       
      SUBROUTINE SHELLS_gpu(NELEC,ISH,JSH,KSH,LSH,FLIP,                 
     * IJGT,IJX,IJY,IJZ,IK,KLGT,KLX,KLY,KLZ,                            
     * IANDJ,KANDL,SAME,                                                
     * INU,JNU,KNU,LNU,NGTI,NGTJ,NGTK,NGTL,                             
     * LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,                             
     * mini,maxi,minj,maxj,mink,maxk,minl,maxl,                         
     * NIJ,IJ,KL,IJKL,                                                  
     * NROOTS,NGTH,C,                                                   
     * kstart,katom,ktype,kng,kloc,kmin,kmax,nshell,                    
     * ex,cs,cp,cd,cf,cg,ch,ci,                                         
     * ga,csa,cpa,cda,cfa,cga,cha,cia,                                  
     * gb,csb,cpb,cdb,cfb,cgb,chb,cib,                                  
     * gc,csc,cpc,cdc,cfc,cgc,chc,cic,                                  
     * gd,csd,cpd,cdd,cfd,cgd,chd,cid,                                  
     * AX,AY,AZ,BX,BY,BZ,RAB,CX,CY,CZ,DX,DY,DZ,RCD,                     
     * NGA,NGB,NGC,NGD)                                                 
      use mx_limits, only: mxsh,mxgsh,mxgtot,mxatm                      
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      LOGICAL FLIP                                                      
      LOGICAL IANDJ,KANDL,SAME                                          
C                                                                       
      DIMENSION IX(84),IY(84),IZ(84),                                   
     *          JX(84),JY(84),JZ(84),                                   
     *          KX(84),KY(84),KZ(84),                                   
     *          LX(84),LY(84),LZ(84)                                    
                                                                        
      INTEGER :: INU,JNU,KNU,LNU,NGTI,NGTJ,NGTK,NGTL                    
      DIMENSION :: IJGT(784),IJX(784),IJY(784),IJZ(784),IK(784)         
      DIMENSION :: KLGT(784),KLX(784),KLY(784),KLZ(784)                 
      double precision :: C(3,MXATM)                                    
                                                                        
      DIMENSION :: EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT)          
      DIMENSION :: CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT)          
      DIMENSION :: KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH)       
      DIMENSION :: KLOC(MXSH),KMIN(MXSH),KMAX(MXSH)                     
      INTEGER :: NSHELL                                                 
                                                                        
      DIMENSION :: NGTH(4)                                              
                                                                        
      DIMENSION :: GA(MXGSH),CSA(MXGSH),CPA(MXGSH),CDA(MXGSH)           
      DIMENSION :: CFA(MXGSH),CGA(MXGSH),CHA(MXGSH),CIA(MXGSH)          
      DIMENSION :: GB(MXGSH),CSB(MXGSH),CPB(MXGSH),CDB(MXGSH)           
      DIMENSION :: CFB(MXGSH),CGB(MXGSH),CHB(MXGSH),CIB(MXGSH)          
      DIMENSION :: GC(MXGSH),CSC(MXGSH),CPC(MXGSH),CDC(MXGSH)           
      DIMENSION :: CFC(MXGSH),CGC(MXGSH),CHC(MXGSH),CIC(MXGSH)          
      DIMENSION :: GD(MXGSH),CSD(MXGSH),CPD(MXGSH),CDD(MXGSH)           
      DIMENSION :: CFD(MXGSH),CGD(MXGSH),CHD(MXGSH),CID(MXGSH)          
      double precision :: AX,AY,AZ,BX,BY,BZ,RAB,CX,CY,CZ,DX,DY,DZ,RCD   
      integer :: NGA,NGB,NGC,NGD                                        
                                                                        
      double precision :: qq4                                           
      INTEGER :: LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL                    
      INTEGER :: MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL                
      INTEGER :: NIJ,IJ,KL,IJKL                                         
                                                                        
      DATA LX /   0,   1,   0,   0,   2,   0,   0,   1,   1,   0,       
     *            3,   0,   0,   2,   2,   1,   0,   1,   0,   1,       
     *            4,   0,   0,   3,   3,   1,   0,   1,   0,   2,       
     *            2,   0,   2,   1,   1,                                
     *            5,   0,   0,   4,   4,   1,   0,   1,   0,   3,       
     *            3,   2,   0,   2,   0,   3,   1,   1,   2,   2,       
     *            1,                                                    
     *            6,   0,   0,   5,   5,   1,   0,   1,   0,   4,       
     *            4,   2,   0,   2,   0,   4,   1,   1,   3,   3,       
     *            0,   3,   3,   2,   1,   2,   1,   2/                 
      DATA KX /   0,   7,   0,   0,  14,   0,   0,   7,   7,   0,       
     *           21,   0,   0,  14,  14,   7,   0,   7,   0,   7,       
     *           28,   0,   0,  21,  21,   7,   0,   7,   0,  14,       
     *           14,   0,  14,   7,   7,                                
     *           35,   0,   0,  28,  28,   7,   0,   7,   0,  21,       
     *           21,  14,   0,  14,   0,  21,   7,   7,  14,  14,       
     *            7,                                                    
     *           42,   0,   0,  35,  35,   7,   0,   7,   0,  28,       
     *           28,  14,   0,  14,   0,  28,   7,   7,  21,  21,       
     *            0,  21,  21,  14,   7,  14,   7,  14/                 
      DATA JX /   0,  49,   0,   0,  98,   0,   0,  49,  49,   0,       
     *          147,   0,   0,  98,  98,  49,   0,  49,   0,  49,       
     *          196,   0,   0, 147, 147,  49,   0,  49,   0,  98,       
     *           98,   0,  98,  49,  49,                                
     *          245,   0,   0, 196, 196,  49,   0,  49,   0, 147,       
     *          147,  98,   0,  98,   0, 147,  49,  49,  98,  98,       
     *           49,                                                    
     *          294,   0,   0, 245, 245,  49,   0,  49,   0, 196,       
     *          196,  98,   0,  98,   0, 196,  49,  49, 147, 147,       
     *            0, 147, 147,  98,  49,  98,  49,  98/                 
      DATA IX /   1, 344,   1,   1, 687,   1,   1, 344, 344,   1,       
     *         1030,   1,   1, 687, 687, 344,   1, 344,   1, 344,       
     *         1373,   1,   1,1030,1030, 344,   1, 344,   1, 687,       
     *          687,   1, 687, 344, 344,                                
     *         1716,   1,   1,1373,1373, 344,   1, 344,   1,1030,       
     *         1030, 687,   1, 687,   1,1030, 344, 344, 687, 687,       
     *          344,                                                    
     *         2059,   1,   1,1716,1716, 344,   1, 344,   1,1373,       
     *         1373, 687,   1, 687,   1,1373, 344, 344,1030,1030,       
     *            1,1030,1030, 687, 344, 687, 344, 687/                 
      DATA LY /   0,   0,   1,   0,   0,   2,   0,   1,   0,   1,       
     *            0,   3,   0,   1,   0,   2,   2,   0,   1,   1,       
     *            0,   4,   0,   1,   0,   3,   3,   0,   1,   2,       
     *            0,   2,   1,   2,   1,                                
     *            0,   5,   0,   1,   0,   4,   4,   0,   1,   2,       
     *            0,   3,   3,   0,   2,   1,   3,   1,   2,   1,       
     *            2,                                                    
     *            0,   6,   0,   1,   0,   5,   5,   0,   1,   2,       
     *            0,   4,   4,   0,   2,   1,   4,   1,   3,   0,       
     *            3,   2,   1,   3,   3,   1,   2,   2/                 
      DATA KY /   0,   0,   7,   0,   0,  14,   0,   7,   0,   7,       
     *            0,  21,   0,   7,   0,  14,  14,   0,   7,   7,       
     *            0,  28,   0,   7,   0,  21,  21,   0,   7,  14,       
     *            0,  14,   7,  14,   7,                                
     *            0,  35,   0,   7,   0,  28,  28,   0,   7,  14,       
     *            0,  21,  21,   0,  14,   7,  21,   7,  14,   7,       
     *           14,                                                    
     *            0,  42,   0,   7,   0,  35,  35,   0,   7,  14,       
     *            0,  28,  28,   0,  14,   7,  28,   7,  21,   0,       
     *           21,  14,   7,  21,  21,   7,  14,  14/                 
      DATA JY /   0,   0,  49,   0,   0,  98,   0,  49,   0,  49,       
     *            0, 147,   0,  49,   0,  98,  98,   0,  49,  49,       
     *            0, 196,   0,  49,   0, 147, 147,   0,  49,  98,       
     *            0,  98,  49,  98,  49,                                
     *            0, 245,   0,  49,   0, 196, 196,   0,  49,  98,       
     *            0, 147, 147,   0,  98,  49, 147,  49,  98,  49,       
     *           98,                                                    
     *            0, 294,   0,  49,   0, 245, 245,   0,  49,  98,       
     *            0, 196, 196,   0,  98,  49, 196,  49, 147,   0,       
     *          147,  98,  49, 147, 147,  49,  98,  98/                 
      DATA IY /   1,   1, 344,   1,   1, 687,   1, 344,   1, 344,       
     *            1,1030,   1, 344,   1, 687, 687,   1, 344, 344,       
     *            1,1373,   1, 344,   1,1030,1030,   1, 344, 687,       
     *            1, 687, 344, 687, 344,                                
     *            1,1716,   1, 344,   1,1373,1373,   1, 344, 687,       
     *            1,1030,1030,   1, 687, 344,1030, 344, 687, 344,       
     *          687,                                                    
     *            1,2059,   1, 344,   1,1716,1716,   1, 344, 687,       
     *            1,1373,1373,   1, 687, 344,1373, 344,1030,   1,       
     *         1030, 687, 344,1030,1030, 344, 687, 687/                 
      DATA LZ /   0,   0,   0,   1,   0,   0,   2,   0,   1,   1,       
     *            0,   0,   3,   0,   1,   0,   1,   2,   2,   1,       
     *            0,   0,   4,   0,   1,   0,   1,   3,   3,   0,       
     *            2,   2,   1,   1,   2,                                
     *            0,   0,   5,   0,   1,   0,   1,   4,   4,   0,       
     *            2,   0,   2,   3,   3,   1,   1,   3,   1,   2,       
     *            2,                                                    
     *            0,   0,   6,   0,   1,   0,   1,   5,   5,   0,       
     *            2,   0,   2,   4,   4,   1,   1,   4,   0,   3,       
     *            3,   1,   2,   1,   2,   3,   3,   2/                 
      DATA KZ /   0,   0,   0,   7,   0,   0,  14,   0,   7,   7,       
     *            0,   0,  21,   0,   7,   0,   7,  14,  14,   7,       
     *            0,   0,  28,   0,   7,   0,   7,  21,  21,   0,       
     *           14,  14,   7,   7,  14,                                
     *            0,   0,  35,   0,   7,   0,   7,  28,  28,   0,       
     *           14,   0,  14,  21,  21,   7,   7,  21,   7,  14,       
     *           14,                                                    
     *            0,   0,  42,   0,   7,   0,   7,  35,  35,   0,       
     *           14,   0,  14,  28,  28,   7,   7,  28,   0,  21,       
     *           21,   7,  14,   7,  14,  21,  21,  14/                 
      DATA JZ /   0,   0,   0,  49,   0,   0,  98,   0,  49,  49,       
     *            0,   0, 147,   0,  49,   0,  49,  98,  98,  49,       
     *            0,   0, 196,   0,  49,   0,  49, 147, 147,   0,       
     *           98,  98,  49,  49,  98,                                
     *            0,   0, 245,   0,  49,   0,  49, 196, 196,   0,       
     *           98,   0,  98, 147, 147,  49,  49, 147,  49,  98,       
     *           98,                                                    
     *            0,   0, 294,   0,  49,   0,  49, 245, 245,   0,       
     *           98,   0,  98, 196, 196,  49,  49, 196,   0, 147,       
     *          147,  49,  98,  49,  98, 147, 147,  98/                 
      DATA IZ /   1,   1,   1, 344,   1,   1, 687,   1, 344, 344,       
     *            1,   1,1030,   1, 344,   1, 344, 687, 687, 344,       
     *            1,   1,1373,   1, 344,   1, 344,1030,1030,   1,       
     *          687, 687, 344, 344, 687,                                
     *            1,   1,1716,   1, 344,   1, 344,1373,1373,   1,       
     *          687,   1, 687,1030,1030, 344, 344,1030, 344, 687,       
     *          687,                                                    
     *            1,   1,2059,   1, 344,   1, 344,1716,1716,   1,       
     *          687,   1, 687,1373,1373, 344, 344,1373,   1,1030,       
     *         1030, 344, 687, 344, 687,1030,1030, 687/                 
C                                                                       
C     PREPARE SHELL INFORMATION/FOR HONDO INTEGRATION                   
C                                                                       
!!$omp declare target                                                   
      IF(NELEC.EQ.2) GO TO 200                                          
C                                                                       
C     ----- PERMUTE ISH AND JSH SHELLS, FOR THEIR TYPE                  
C     THIS IS DONE FOR SPEED REASONS.  THE CODE GETS THE RIGHT ANSWER   
C     WITHOUT THE ANGULAR MOMENTUM FLIPPING, AND THEREFORE A CALLING    
C     ARGUMENT ALLOWS ONE DO EXACTLY THE INTEGRAL BLOCK AS SPECIFIED,   
C     SHOULD THAT BE DESIRED.                                           
C                                                                       
      IANDJ = ISH .EQ. JSH                                              
C      IF (KTYPE(ISH) .LT. KTYPE(JSH)  .AND.  FLIP) THEN                
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
C                                                                       
C     ----- ISHELL                                                      
C                                                                       
      I = KATOM(INU)                                                    
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
      DO 140 I = I1,I2                                                  
         NGA = NGA+1                                                    
         GA(NGA) = EX(I)                                                
C         write(*,*) "GA", GA                                           
         CSA(NGA) = CS(I)                                               
         CPA(NGA) = CP(I)                                               
         CDA(NGA) = CD(I)                                               
         CFA(NGA) = CF(I)                                               
         CGA(NGA) = CG(I)                                               
         CHA(NGA) = CH(I)                                               
C         write(*,*) "CHA", CHA                                         
         CIA(NGA) = CI(I)                                               
  140 CONTINUE                                                          
C                                                                       
C     ----- JSHELL                                                      
C                                                                       
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
      DO 160 J = J1,J2                                                  
         NGB = NGB+1                                                    
         GB(NGB) = EX(J)                                                
         CSB(NGB) = CS(J)                                               
         CPB(NGB) = CP(J)                                               
         CDB(NGB) = CD(J)                                               
         CFB(NGB) = CF(J)                                               
         CGB(NGB) = CG(J)                                               
         CHB(NGB) = CH(J)                                               
         CIB(NGB) = CI(J)                                               
  160 CONTINUE                                                          
      RAB = ((AX-BX)*(AX-BX) + (AY-BY)*(AY-BY) + (AZ-BZ)*(AZ-BZ))       
C                                                                       
C     ----- PREPARE INDICES FOR PAIRS OF (I,J) FUNCTIONS                
C                                                                       
      IJ = 0                                                            
      JMAX = MAXJ                                                       
      DO 190 I = MINI,MAXI                                              
         NX = IX(I)                                                     
         NY = IY(I)                                                     
         NZ = IZ(I)                                                     
         IF (IANDJ) JMAX = I                                            
         DO 180 J = MINJ,JMAX                                           
            IJ = IJ+1                                                   
            IJX(IJ) = NX+JX(J)                                          
            IJY(IJ) = NY+JY(J)                                          
            IJZ(IJ) = NZ+JZ(J)                                          
            IJGT(IJ) = NGTI*(I-MINI)+NGTJ*(J-MINJ)+1                    
  180    CONTINUE                                                       
  190 CONTINUE                                                          
      RETURN                                                            
C     ******                                                            
C                                                                       
C        K AND L SHELL                                                  
C                                                                       
  200 CONTINUE                                                          
      KANDL = KSH .EQ. LSH                                              
      SAME = ISH .EQ. KSH .AND. JSH .EQ. LSH                            
C                                                                       
C     ----- PERMUTE KSH AND LSH SHELLS, FOR THEIR TYPE                  
C                                                                       
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
C                                                                       
C     ----- K SHELL                                                     
C                                                                       
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
      DO 260 K = K1,K2                                                  
         NGC = NGC+1                                                    
         GC(NGC) = EX(K)                                                
         CSC(NGC) = CS(K)                                               
         CPC(NGC) = CP(K)                                               
         CDC(NGC) = CD(K)                                               
         CFC(NGC) = CF(K)                                               
         CGC(NGC) = CG(K)                                               
         CHC(NGC) = CH(K)                                               
         CIC(NGC) = CI(K)                                               
  260 CONTINUE                                                          
C                                                                       
C     ----- LSHELL                                                      
C                                                                       
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
      DO 280 L = L1,L2                                                  
         NGD = NGD+1                                                    
         GD(NGD) = EX(L)                                                
         CSD(NGD) = CS(L)                                               
         CPD(NGD) = CP(L)                                               
         CDD(NGD) = CD(L)                                               
         CFD(NGD) = CF(L)                                               
         CGD(NGD) = CG(L)                                               
         CHD(NGD) = CH(L)                                               
         CID(NGD) = CI(L)                                               
  280 CONTINUE                                                          
C      NROOTS = (LIT+LJT+LKT+LLT-2)/2                                   
      RCD = ((CX-DX)*(CX-DX) + (CY-DY)*(CY-DY) + (CZ-DZ)*(CZ-DZ))       
C                                                                       
C     ----- PREPARE INDICES FOR PAIRS OF (K,L) FUNCTIONS                
C                                                                       
      KL = 0                                                            
      LMAX = MAXL                                                       
      DO 310 K = MINK,MAXK                                              
         NX = KX(K)                                                     
         NY = KY(K)                                                     
         NZ = KZ(K)                                                     
         IF (KANDL) LMAX = K                                            
         DO 300 L = MINL,LMAX                                           
            KL = KL+1                                                   
            KLX(KL) = NX+LX(L)                                          
            KLY(KL) = NY+LY(L)                                          
            KLZ(KL) = NZ+LZ(L)                                          
            KLGT(KL) = NGTK*(K-MINK)+NGTL*(L-MINL)                      
  300    CONTINUE                                                       
  310 CONTINUE                                                          
      MAX = KL                                                          
      DO 320 I = 1,IJ                                                   
      IF (SAME) MAX = I                                                 
  320 IK(I) = MAX                                                       
      IJKL = IJ*KL                                                      
      IF (SAME) IJKL = IJ*(IJ+1)/2                                      
      RETURN                                                            
      END                                                               
C*MODULE INT2A   *DECK SHELLS_gpu2_ij                                   
C       SUBROUTINE SHELLS_gpu2_ij(NELEC,ISH,JSH,FLIP,                   
C      * IJGT,IJX,IJY,IJZ,IK,                                           
C      * IANDJ,SAME,                                                    
C      * INU,JNU,KNU,LNU,NGTI,NGTJ,                                     
C      * LIT,LJT,LOCI,LOCJ,                                             
C      * mini,maxi,minj,maxj,                                           
C      * NIJ,IJ,KL,IJKL,                                                
C      * NROOTS,NGTH,C,                                                 
C      * kstart,katom,ktype,kng,kloc,kmin,kmax,nshell,                  
C      * ex,cs,cp,cd,                                                   
C      * ga,csa,cpa,cda,                                                
C      * gb,csb,cpb,cdb,                                                
C      * AX,AY,AZ,BX,BY,BZ,RAB,                                         
C      * NGA,NGB)                                                       
C       use mx_limits, only: mxsh,mxgsh,mxgtot,mxatm                    
C C                                                                     
C       IMPLICIT DOUBLE PRECISION (A-H,O-Z)                             
C C                                                                     
C       LOGICAL FLIP                                                    
C       LOGICAL IANDJ,KANDL,SAME                                        
C C                                                                     
C       DIMENSION IX(84),IY(84),IZ(84),                                 
C      *          JX(84),JY(84),JZ(84),                                 
C      *          KX(84),KY(84),KZ(84),                                 
C      *          LX(84),LY(84),LZ(84)                                  
                                                                        
C       INTEGER :: INU,JNU,KNU,LNU,NGTI,NGTJ,NGTK,NGTL                  
C       DIMENSION :: IJGT(784),IJX(784),IJY(784),IJZ(784),IK(784)       
C       DIMENSION :: KLGT(784),KLX(784),KLY(784),KLZ(784)               
C       double precision :: C(3,MXATM)                                  
                                                                        
C       DIMENSION :: EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT)        
C       DIMENSION :: CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT)        
C       DIMENSION :: KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH)     
C       DIMENSION :: KLOC(MXSH),KMIN(MXSH),KMAX(MXSH)                   
C       INTEGER :: NSHELL                                               
                                                                        
C       DIMENSION :: NGTH(4)                                            
                                                                        
C       DIMENSION :: GA(MXGSH),CSA(MXGSH),CPA(MXGSH),CDA(MXGSH)         
C       DIMENSION :: CFA(MXGSH),CGA(MXGSH),CHA(MXGSH),CIA(MXGSH)        
C       DIMENSION :: GB(MXGSH),CSB(MXGSH),CPB(MXGSH),CDB(MXGSH)         
C       DIMENSION :: CFB(MXGSH),CGB(MXGSH),CHB(MXGSH),CIB(MXGSH)        
C       DIMENSION :: GC(MXGSH),CSC(MXGSH),CPC(MXGSH),CDC(MXGSH)         
C       DIMENSION :: CFC(MXGSH),CGC(MXGSH),CHC(MXGSH),CIC(MXGSH)        
C       DIMENSION :: GD(MXGSH),CSD(MXGSH),CPD(MXGSH),CDD(MXGSH)         
C       DIMENSION :: CFD(MXGSH),CGD(MXGSH),CHD(MXGSH),CID(MXGSH)        
C       double precision :: AX,AY,AZ,BX,BY,BZ,RAB,CX,CY,CZ,DX,DY,DZ,RCD 
C       integer :: NGA,NGB,NGC,NGD                                      
                                                                        
C       double precision :: qq4                                         
C       INTEGER :: LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL                  
C       INTEGER :: MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL              
C       INTEGER :: NIJ,IJ,KL,IJKL                                       
C       DATA JX /   0,  49,   0,   0,  98,   0,   0,  49,  49,   0,     
C      *          147,   0,   0,  98,  98,  49,   0,  49,   0,  49,     
C      *          196,   0,   0, 147, 147,  49,   0,  49,   0,  98,     
C      *           98,   0,  98,  49,  49,                              
C      *          245,   0,   0, 196, 196,  49,   0,  49,   0, 147,     
C      *          147,  98,   0,  98,   0, 147,  49,  49,  98,  98,     
C      *           49,                                                  
C      *          294,   0,   0, 245, 245,  49,   0,  49,   0, 196,     
C      *          196,  98,   0,  98,   0, 196,  49,  49, 147, 147,     
C      *            0, 147, 147,  98,  49,  98,  49,  98/               
C       DATA IX /   1, 344,   1,   1, 687,   1,   1, 344, 344,   1,     
C      *         1030,   1,   1, 687, 687, 344,   1, 344,   1, 344,     
C      *         1373,   1,   1,1030,1030, 344,   1, 344,   1, 687,     
C      *          687,   1, 687, 344, 344,                              
C      *         1716,   1,   1,1373,1373, 344,   1, 344,   1,1030,     
C      *         1030, 687,   1, 687,   1,1030, 344, 344, 687, 687,     
C      *          344,                                                  
C      *         2059,   1,   1,1716,1716, 344,   1, 344,   1,1373,     
C      *         1373, 687,   1, 687,   1,1373, 344, 344,1030,1030,     
C      *            1,1030,1030, 687, 344, 687, 344, 687/               
C       DATA JY /   0,   0,  49,   0,   0,  98,   0,  49,   0,  49,     
C      *            0, 147,   0,  49,   0,  98,  98,   0,  49,  49,     
C      *            0, 196,   0,  49,   0, 147, 147,   0,  49,  98,     
C      *            0,  98,  49,  98,  49,                              
C      *            0, 245,   0,  49,   0, 196, 196,   0,  49,  98,     
C      *            0, 147, 147,   0,  98,  49, 147,  49,  98,  49,     
C      *           98,                                                  
C      *            0, 294,   0,  49,   0, 245, 245,   0,  49,  98,     
C      *            0, 196, 196,   0,  98,  49, 196,  49, 147,   0,     
C      *          147,  98,  49, 147, 147,  49,  98,  98/               
C       DATA IY /   1,   1, 344,   1,   1, 687,   1, 344,   1, 344,     
C      *            1,1030,   1, 344,   1, 687, 687,   1, 344, 344,     
C      *            1,1373,   1, 344,   1,1030,1030,   1, 344, 687,     
C      *            1, 687, 344, 687, 344,                              
C      *            1,1716,   1, 344,   1,1373,1373,   1, 344, 687,     
C      *            1,1030,1030,   1, 687, 344,1030, 344, 687, 344,     
C      *          687,                                                  
C      *            1,2059,   1, 344,   1,1716,1716,   1, 344, 687,     
C      *            1,1373,1373,   1, 687, 344,1373, 344,1030,   1,     
C      *         1030, 687, 344,1030,1030, 344, 687, 687/               
C       DATA JZ /   0,   0,   0,  49,   0,   0,  98,   0,  49,  49,     
C      *            0,   0, 147,   0,  49,   0,  49,  98,  98,  49,     
C      *            0,   0, 196,   0,  49,   0,  49, 147, 147,   0,     
C      *           98,  98,  49,  49,  98,                              
C      *            0,   0, 245,   0,  49,   0,  49, 196, 196,   0,     
C      *           98,   0,  98, 147, 147,  49,  49, 147,  49,  98,     
C      *           98,                                                  
C      *            0,   0, 294,   0,  49,   0,  49, 245, 245,   0,     
C      *           98,   0,  98, 196, 196,  49,  49, 196,   0, 147,     
C      *          147,  49,  98,  49,  98, 147, 147,  98/               
C       DATA IZ /   1,   1,   1, 344,   1,   1, 687,   1, 344, 344,     
C      *            1,   1,1030,   1, 344,   1, 344, 687, 687, 344,     
C      *            1,   1,1373,   1, 344,   1, 344,1030,1030,   1,     
C      *          687, 687, 344, 344, 687,                              
C      *            1,   1,1716,   1, 344,   1, 344,1373,1373,   1,     
C      *          687,   1, 687,1030,1030, 344, 344,1030, 344, 687,     
C      *          687,                                                  
C      *            1,   1,2059,   1, 344,   1, 344,1716,1716,   1,     
C      *          687,   1, 687,1373,1373, 344, 344,1373,   1,1030,     
C      *         1030, 344, 687, 344, 687,1030,1030, 687/               
C !$omp declare target                                                  
C       IANDJ = ISH .EQ. JSH                                            
C       IF (KTYPE(ISH) .LT. KTYPE(JSH)) THEN                            
C          INU = JSH                                                    
C          JNU = ISH                                                    
C          NGTI = NGTH(2)                                               
C          NGTJ = NGTH(1)                                               
C       ELSE                                                            
C          INU = ISH                                                    
C          JNU = JSH                                                    
C          NGTI = NGTH(1)                                               
C          NGTJ = NGTH(2)                                               
C       END IF                                                          
C C                                                                     
C C     ----- ISHELL                                                    
C       I = KATOM(INU)                                                  
C       AX = C(1,I)                                                     
C       AY = C(2,I)                                                     
C       AZ = C(3,I)                                                     
C       I1 = KSTART(INU)                                                
C       I2 = I1+KNG(INU)-1                                              
C       LIT = KTYPE(INU)                                                
C       MINI = KMIN(INU)                                                
C       MAXI = KMAX(INU)                                                
C       LOCI = KLOC(INU)-MINI                                           
C       NGA = 0                                                         
C       DO I = I1,I2                                                    
C          NGA = NGA+1                                                  
C          GA(NGA) = EX(I)                                              
C          CSA(NGA) = CS(I)                                             
C          CPA(NGA) = CP(I)                                             
C          CDA(NGA) = CD(I)                                             
C C          FA(NGA) = CF(I)                                            
C C          GA(NGA) = CG(I)                                            
C C          HA(NGA) = CH(I)                                            
C C          IA(NGA) = CI(I)                                            
C       ENDDO                                                           
C C                                                                     
C C     ----- JSHELL                                                    
C       J = KATOM(JNU)                                                  
C       BX = C(1,J)                                                     
C       BY = C(2,J)                                                     
C       BZ = C(3,J)                                                     
C       J1 = KSTART(JNU)                                                
C       J2 = J1+KNG(JNU)-1                                              
C       LJT = KTYPE(JNU)                                                
C       MINJ = KMIN(JNU)                                                
C       MAXJ = KMAX(JNU)                                                
C       LOCJ = KLOC(JNU)-MINJ                                           
C       NGB = 0                                                         
C       DO J = J1,J2                                                    
C          NGB = NGB+1                                                  
C          GB(NGB) = EX(J)                                              
C          CSB(NGB) = CS(J)                                             
C          CPB(NGB) = CP(J)                                             
C          CDB(NGB) = CD(J)                                             
C C          FB(NGB) = CF(J)                                            
C C          GB(NGB) = CG(J)                                            
C C          HB(NGB) = CH(J)                                            
C C          IB(NGB) = CI(J)                                            
C       ENDDO                                                           
C       RAB = ((AX-BX)*(AX-BX) + (AY-BY)*(AY-BY) + (AZ-BZ)*(AZ-BZ))     
C C                                                                     
C C     ----- PREPARE INDICES FOR PAIRS OF (I,J) FUNCTIONS              
C C                                                                     
C       IJ = 0                                                          
C       JMAX = MAXJ                                                     
C       DO I = MINI,MAXI                                                
C          NX = IX(I)                                                   
C C         write(*,*) "NX is", NX                                      
C          NY = IY(I)                                                   
C          NZ = IZ(I)                                                   
C          IF (IANDJ) JMAX = I                                          
C          DO J = MINJ,JMAX                                             
C C             IJ = IJ+1                                               
C C             IJX(IJ) = NX+JX(J)                                      
C C             IJY(IJ) = NY+JY(J)                                      
C C             IJZ(IJ) = NZ+JZ(J)                                      
C C             IJGT(IJ) = NGTI*(I-MINI)+NGTJ*(J-MINJ)+1                
C          ENDDO                                                        
C       ENDDO                                                           
C       RETURN                                                          
C       END                                                             
C C*MODULE INT2A   *DECK SHELLS_gpu2_kl                                 
C       SUBROUTINE SHELLS_gpu2_kl(NELEC,ISH,JSH,KSH,LSH,FLIP,           
C      * IJGT,IJX,IJY,IJZ,IK,KLGT,KLX,KLY,KLZ,                          
C      * IANDJ,KANDL,SAME,                                              
C      * INU,JNU,KNU,LNU,NGTI,NGTJ,NGTK,NGTL,                           
C      * LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,                           
C      * mini,maxi,minj,maxj,mink,maxk,minl,maxl,                       
C      * NIJ,IJ,KL,IJKL,                                                
C      * NROOTS,NGTH,C,                                                 
C      * kstart,katom,ktype,kng,kloc,kmin,kmax,nshell,                  
C      * ex,cs,cp,cd,cf,cg,ch,ci,                                       
C      * ga,csa,cpa,cda,cfa,cga,cha,cia,                                
C      * gb,csb,cpb,cdb,cfb,cgb,chb,cib,                                
C      * gc,csc,cpc,cdc,cfc,cgc,chc,cic,                                
C      * gd,csd,cpd,cdd,cfd,cgd,chd,cid,                                
C      * AX,AY,AZ,BX,BY,BZ,RAB,CX,CY,CZ,DX,DY,DZ,RCD,                   
C      * NGA,NGB,NGC,NGD)                                               
C       use mx_limits, only: mxsh,mxgsh,mxgtot,mxatm                    
C C                                                                     
C       IMPLICIT DOUBLE PRECISION (A-H,O-Z)                             
C C                                                                     
C       LOGICAL FLIP                                                    
C       LOGICAL IANDJ,KANDL,SAME                                        
C C                                                                     
C       DIMENSION IX(84),IY(84),IZ(84),                                 
C      *          JX(84),JY(84),JZ(84),                                 
C      *          KX(84),KY(84),KZ(84),                                 
C      *          LX(84),LY(84),LZ(84)                                  
                                                                        
C       INTEGER :: INU,JNU,KNU,LNU,NGTI,NGTJ,NGTK,NGTL                  
C       DIMENSION :: IJGT(784),IJX(784),IJY(784),IJZ(784),IK(784)       
C       DIMENSION :: KLGT(784),KLX(784),KLY(784),KLZ(784)               
C       double precision :: C(3,MXATM)                                  
                                                                        
C       DIMENSION :: EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT)        
C       DIMENSION :: CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT)        
C       DIMENSION :: KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH)     
C       DIMENSION :: KLOC(MXSH),KMIN(MXSH),KMAX(MXSH)                   
C       INTEGER :: NSHELL                                               
                                                                        
C       DIMENSION :: NGTH(4)                                            
                                                                        
C       DIMENSION :: GA(MXGSH),CSA(MXGSH),CPA(MXGSH),CDA(MXGSH)         
C       DIMENSION :: CFA(MXGSH),CGA(MXGSH),CHA(MXGSH),CIA(MXGSH)        
C       DIMENSION :: GB(MXGSH),CSB(MXGSH),CPB(MXGSH),CDB(MXGSH)         
C       DIMENSION :: CFB(MXGSH),CGB(MXGSH),CHB(MXGSH),CIB(MXGSH)        
C       DIMENSION :: GC(MXGSH),CSC(MXGSH),CPC(MXGSH),CDC(MXGSH)         
C       DIMENSION :: CFC(MXGSH),CGC(MXGSH),CHC(MXGSH),CIC(MXGSH)        
C       DIMENSION :: GD(MXGSH),CSD(MXGSH),CPD(MXGSH),CDD(MXGSH)         
C       DIMENSION :: CFD(MXGSH),CGD(MXGSH),CHD(MXGSH),CID(MXGSH)        
C       double precision :: AX,AY,AZ,BX,BY,BZ,RAB,CX,CY,CZ,DX,DY,DZ,RCD 
C       integer :: NGA,NGB,NGC,NGD                                      
                                                                        
C       double precision :: qq4                                         
C       INTEGER :: LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL                  
C       INTEGER :: MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL              
C       INTEGER :: NIJ,IJ,KL,IJKL                                       
                                                                        
C       DATA LX /   0,   1,   0,   0,   2,   0,   0,   1,   1,   0,     
C      *            3,   0,   0,   2,   2,   1,   0,   1,   0,   1,     
C      *            4,   0,   0,   3,   3,   1,   0,   1,   0,   2,     
C      *            2,   0,   2,   1,   1,                              
C      *            5,   0,   0,   4,   4,   1,   0,   1,   0,   3,     
C      *            3,   2,   0,   2,   0,   3,   1,   1,   2,   2,     
C      *            1,                                                  
C      *            6,   0,   0,   5,   5,   1,   0,   1,   0,   4,     
C      *            4,   2,   0,   2,   0,   4,   1,   1,   3,   3,     
C      *            0,   3,   3,   2,   1,   2,   1,   2/               
C       DATA KX /   0,   7,   0,   0,  14,   0,   0,   7,   7,   0,     
C      *           21,   0,   0,  14,  14,   7,   0,   7,   0,   7,     
C      *           28,   0,   0,  21,  21,   7,   0,   7,   0,  14,     
C      *           14,   0,  14,   7,   7,                              
C      *           35,   0,   0,  28,  28,   7,   0,   7,   0,  21,     
C      *           21,  14,   0,  14,   0,  21,   7,   7,  14,  14,     
C      *            7,                                                  
C      *           42,   0,   0,  35,  35,   7,   0,   7,   0,  28,     
C      *           28,  14,   0,  14,   0,  28,   7,   7,  21,  21,     
C      *            0,  21,  21,  14,   7,  14,   7,  14/               
C       DATA JX /   0,  49,   0,   0,  98,   0,   0,  49,  49,   0,     
C      *          147,   0,   0,  98,  98,  49,   0,  49,   0,  49,     
C      *          196,   0,   0, 147, 147,  49,   0,  49,   0,  98,     
C      *           98,   0,  98,  49,  49,                              
C      *          245,   0,   0, 196, 196,  49,   0,  49,   0, 147,     
C      *          147,  98,   0,  98,   0, 147,  49,  49,  98,  98,     
C      *           49,                                                  
C      *          294,   0,   0, 245, 245,  49,   0,  49,   0, 196,     
C      *          196,  98,   0,  98,   0, 196,  49,  49, 147, 147,     
C      *            0, 147, 147,  98,  49,  98,  49,  98/               
C       DATA IX /   1, 344,   1,   1, 687,   1,   1, 344, 344,   1,     
C      *         1030,   1,   1, 687, 687, 344,   1, 344,   1, 344,     
C      *         1373,   1,   1,1030,1030, 344,   1, 344,   1, 687,     
C      *          687,   1, 687, 344, 344,                              
C      *         1716,   1,   1,1373,1373, 344,   1, 344,   1,1030,     
C      *         1030, 687,   1, 687,   1,1030, 344, 344, 687, 687,     
C      *          344,                                                  
C      *         2059,   1,   1,1716,1716, 344,   1, 344,   1,1373,     
C      *         1373, 687,   1, 687,   1,1373, 344, 344,1030,1030,     
C      *            1,1030,1030, 687, 344, 687, 344, 687/               
C       DATA LY /   0,   0,   1,   0,   0,   2,   0,   1,   0,   1,     
C      *            0,   3,   0,   1,   0,   2,   2,   0,   1,   1,     
C      *            0,   4,   0,   1,   0,   3,   3,   0,   1,   2,     
C      *            0,   2,   1,   2,   1,                              
C      *            0,   5,   0,   1,   0,   4,   4,   0,   1,   2,     
C      *            0,   3,   3,   0,   2,   1,   3,   1,   2,   1,     
C      *            2,                                                  
C      *            0,   6,   0,   1,   0,   5,   5,   0,   1,   2,     
C      *            0,   4,   4,   0,   2,   1,   4,   1,   3,   0,     
C      *            3,   2,   1,   3,   3,   1,   2,   2/               
C       DATA KY /   0,   0,   7,   0,   0,  14,   0,   7,   0,   7,     
C      *            0,  21,   0,   7,   0,  14,  14,   0,   7,   7,     
C      *            0,  28,   0,   7,   0,  21,  21,   0,   7,  14,     
C      *            0,  14,   7,  14,   7,                              
C      *            0,  35,   0,   7,   0,  28,  28,   0,   7,  14,     
C      *            0,  21,  21,   0,  14,   7,  21,   7,  14,   7,     
C      *           14,                                                  
C      *            0,  42,   0,   7,   0,  35,  35,   0,   7,  14,     
C      *            0,  28,  28,   0,  14,   7,  28,   7,  21,   0,     
C      *           21,  14,   7,  21,  21,   7,  14,  14/               
C       DATA JY /   0,   0,  49,   0,   0,  98,   0,  49,   0,  49,     
C      *            0, 147,   0,  49,   0,  98,  98,   0,  49,  49,     
C      *            0, 196,   0,  49,   0, 147, 147,   0,  49,  98,     
C      *            0,  98,  49,  98,  49,                              
C      *            0, 245,   0,  49,   0, 196, 196,   0,  49,  98,     
C      *            0, 147, 147,   0,  98,  49, 147,  49,  98,  49,     
C      *           98,                                                  
C      *            0, 294,   0,  49,   0, 245, 245,   0,  49,  98,     
C      *            0, 196, 196,   0,  98,  49, 196,  49, 147,   0,     
C      *          147,  98,  49, 147, 147,  49,  98,  98/               
C       DATA IY /   1,   1, 344,   1,   1, 687,   1, 344,   1, 344,     
C      *            1,1030,   1, 344,   1, 687, 687,   1, 344, 344,     
C      *            1,1373,   1, 344,   1,1030,1030,   1, 344, 687,     
C      *            1, 687, 344, 687, 344,                              
C      *            1,1716,   1, 344,   1,1373,1373,   1, 344, 687,     
C      *            1,1030,1030,   1, 687, 344,1030, 344, 687, 344,     
C      *          687,                                                  
C      *            1,2059,   1, 344,   1,1716,1716,   1, 344, 687,     
C      *            1,1373,1373,   1, 687, 344,1373, 344,1030,   1,     
C      *         1030, 687, 344,1030,1030, 344, 687, 687/               
C       DATA LZ /   0,   0,   0,   1,   0,   0,   2,   0,   1,   1,     
C      *            0,   0,   3,   0,   1,   0,   1,   2,   2,   1,     
C      *            0,   0,   4,   0,   1,   0,   1,   3,   3,   0,     
C      *            2,   2,   1,   1,   2,                              
C      *            0,   0,   5,   0,   1,   0,   1,   4,   4,   0,     
C      *            2,   0,   2,   3,   3,   1,   1,   3,   1,   2,     
C      *            2,                                                  
C      *            0,   0,   6,   0,   1,   0,   1,   5,   5,   0,     
C      *            2,   0,   2,   4,   4,   1,   1,   4,   0,   3,     
C      *            3,   1,   2,   1,   2,   3,   3,   2/               
C       DATA KZ /   0,   0,   0,   7,   0,   0,  14,   0,   7,   7,     
C      *            0,   0,  21,   0,   7,   0,   7,  14,  14,   7,     
C      *            0,   0,  28,   0,   7,   0,   7,  21,  21,   0,     
C      *           14,  14,   7,   7,  14,                              
C      *            0,   0,  35,   0,   7,   0,   7,  28,  28,   0,     
C      *           14,   0,  14,  21,  21,   7,   7,  21,   7,  14,     
C      *           14,                                                  
C      *            0,   0,  42,   0,   7,   0,   7,  35,  35,   0,     
C      *           14,   0,  14,  28,  28,   7,   7,  28,   0,  21,     
C      *           21,   7,  14,   7,  14,  21,  21,  14/               
C       DATA JZ /   0,   0,   0,  49,   0,   0,  98,   0,  49,  49,     
C      *            0,   0, 147,   0,  49,   0,  49,  98,  98,  49,     
C      *            0,   0, 196,   0,  49,   0,  49, 147, 147,   0,     
C      *           98,  98,  49,  49,  98,                              
C      *            0,   0, 245,   0,  49,   0,  49, 196, 196,   0,     
C      *           98,   0,  98, 147, 147,  49,  49, 147,  49,  98,     
C      *           98,                                                  
C      *            0,   0, 294,   0,  49,   0,  49, 245, 245,   0,     
C      *           98,   0,  98, 196, 196,  49,  49, 196,   0, 147,     
C      *          147,  49,  98,  49,  98, 147, 147,  98/               
C       DATA IZ /   1,   1,   1, 344,   1,   1, 687,   1, 344, 344,     
C      *            1,   1,1030,   1, 344,   1, 344, 687, 687, 344,     
C      *            1,   1,1373,   1, 344,   1, 344,1030,1030,   1,     
C      *          687, 687, 344, 344, 687,                              
C      *            1,   1,1716,   1, 344,   1, 344,1373,1373,   1,     
C      *          687,   1, 687,1030,1030, 344, 344,1030, 344, 687,     
C      *          687,                                                  
C      *            1,   1,2059,   1, 344,   1, 344,1716,1716,   1,     
C      *          687,   1, 687,1373,1373, 344, 344,1373,   1,1030,     
C      *         1030, 344, 687, 344, 687,1030,1030, 687/               
C C                                                                     
C C     PREPARE SHELL INFORMATION/FOR HONDO INTEGRATION                 
                                                                        
C       KANDL = KSH .EQ. LSH                                            
C       SAME = ISH .EQ. KSH .AND. JSH .EQ. LSH                          
C C                                                                     
C C     ----- PERMUTE KSH AND LSH SHELLS, FOR THEIR TYPE                
C C                                                                     
C       IF (KTYPE(KSH) .LT. KTYPE(LSH)  .AND.  FLIP) THEN               
C          KNU = LSH                                                    
C          LNU = KSH                                                    
C          NGTK = NGTH(4)                                               
C          NGTL = NGTH(3)                                               
C       ELSE                                                            
C          KNU = KSH                                                    
C          LNU = LSH                                                    
C          NGTK = NGTH(3)                                               
C          NGTL = NGTH(4)                                               
C       END IF                                                          
C C                                                                     
C C     ----- K SHELL                                                   
C C                                                                     
C       K = KATOM(KNU)                                                  
C       CX = C(1,K)                                                     
C       CY = C(2,K)                                                     
C       CZ = C(3,K)                                                     
C       K1 = KSTART(KNU)                                                
C       K2 = K1+KNG(KNU)-1                                              
C       LKT = KTYPE(KNU)                                                
C       MINK = KMIN(KNU)                                                
C       MAXK = KMAX(KNU)                                                
C       LOCK = KLOC(KNU)-MINK                                           
C       NGC = 0                                                         
C       DO K = K1,K2                                                    
C          NGC = NGC+1                                                  
C          GC(NGC) = EX(K)                                              
C          CSC(NGC) = CS(K)                                             
C          CPC(NGC) = CP(K)                                             
C          CDC(NGC) = CD(K)                                             
C          CFC(NGC) = CF(K)                                             
C          CGC(NGC) = CG(K)                                             
C          CHC(NGC) = CH(K)                                             
C          CIC(NGC) = CI(K)                                             
C       ENDDO                                                           
C C                                                                     
C C     ----- LSHELL                                                    
C C                                                                     
C       L = KATOM(LNU)                                                  
C       DX = C(1,L)                                                     
C       DY = C(2,L)                                                     
C       DZ = C(3,L)                                                     
C       L1 = KSTART(LNU)                                                
C       L2 = L1+KNG(LNU)-1                                              
C       LLT = KTYPE(LNU)                                                
C       MINL = KMIN(LNU)                                                
C       MAXL = KMAX(LNU)                                                
C       LOCL = KLOC(LNU)-MINL                                           
C       NGD = 0                                                         
C       DO L = L1,L2                                                    
C          NGD = NGD+1                                                  
C          GD(NGD) = EX(L)                                              
C          CSD(NGD) = CS(L)                                             
C          CPD(NGD) = CP(L)                                             
C          CDD(NGD) = CD(L)                                             
C          CFD(NGD) = CF(L)                                             
C          CGD(NGD) = CG(L)                                             
C          CHD(NGD) = CH(L)                                             
C          CID(NGD) = CI(L)                                             
C       ENDDO                                                           
C C      NROOTS = (LIT+LJT+LKT+LLT-2)/2                                 
C       RCD = ((CX-DX)*(CX-DX) + (CY-DY)*(CY-DY) + (CZ-DZ)*(CZ-DZ))     
C C                                                                     
C C     ----- PREPARE INDICES FOR PAIRS OF (K,L) FUNCTIONS              
C       KL = 0                                                          
C       LMAX = MAXL                                                     
C       DO K = MINK,MAXK                                                
C          NX = KX(K)                                                   
C          NY = KY(K)                                                   
C          NZ = KZ(K)                                                   
C          IF (KANDL) LMAX = K                                          
C          DO L = MINL,LMAX                                             
C             KL = KL+1                                                 
C             KLX(KL) = NX+LX(L)                                        
C             KLY(KL) = NY+LY(L)                                        
C             KLZ(KL) = NZ+LZ(L)                                        
C             KLGT(KL) = NGTK*(K-MINK)+NGTL*(L-MINL)                    
C          ENDDO                                                        
C       ENDDO                                                           
C       MAX = KL                                                        
C       DO 320 I = 1,IJ                                                 
C       IF (SAME) MAX = I                                               
C   320 IK(I) = MAX                                                     
C       IJKL = IJ*KL                                                    
C       IF (SAME) IJKL = IJ*(IJ+1)/2                                    
C       RETURN                                                          
C       END                                                             
C*MODULE INT2A   *DECK FORMS_gpu                                        
      SUBROUTINE FORMS_gpu(GHONDO,NROOTS,DKL,DIJ,XIN,YIN,ZIN,           
     * IJGT,IJX,IJY,IJZ,IK,KLGT,KLX,KLY,KLZ,IJ)                         
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      DIMENSION GHONDO(*)                                               
C                                                                       
      DIMENSION :: DKL(784),DIJ(784)                                    
      DIMENSION :: IJGT(784),IJX(784),IJY(784),IJZ(784),IK(784)         
      DIMENSION :: KLGT(784),KLX(784),KLY(784),KLZ(784)                 
      DIMENSION :: XIN(31213),YIN(31213),ZIN(31213)                     
C                                                                       
C     ----- FORM INTEGRALS OVER FUNCTIONS -----                         
C     DIMENSIONING XIN(81,5), AND ROLLING UP THE COMPUTATION            
C     OF GHONDO IN A LOOP OF LENGTH NROOTS ADDS 33 SECONDS TO           
C     A 240 SECOND INTEGRAL COMPUTATION JOB.  LEAVE IT UNROLLED.        
C                                                                       
!!$omp declare target                                                   
      GO TO (10,20,30,40,50,60,70,80,90,100,110,120,130),NROOTS         
C                                                                       
C          CODE FOR NROOTS=1                                            
C                                                                       
   10 CONTINUE                                                          
      DO 12 I = 1,IJ                                                    
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO 11 K = 1,MAX                                                   
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*                                
     *            ( XIN(MX      )*YIN(MY      )*ZIN(MZ      ))          
   11 CONTINUE                                                          
   12 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C          CODE FOR NROOTS=2                                            
C                                                                       
   20 CONTINUE                                                          
      DO 22 I = 1,IJ                                                    
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO 21 K = 1,MAX                                                   
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*                                
     *            ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )           
     *          +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401))          
   21 CONTINUE                                                          
   22 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C          CODE FOR NROOTS=3                                            
C                                                                       
   30 CONTINUE                                                          
      DO 32 I = 1,IJ                                                    
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO 31 K = 1,MAX                                                   
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*                                
     *            ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )           
     *          +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401)           
     *          +   XIN(MX+ 4802)*YIN(MY+ 4802)*ZIN(MZ+ 4802))          
   31 CONTINUE                                                          
   32 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C          CODE FOR NROOTS=4                                            
C                                                                       
   40 CONTINUE                                                          
      DO 42 I = 1,IJ                                                    
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO 41 K = 1,MAX                                                   
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*                                
     *            ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )           
     *          +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401)           
     *          +   XIN(MX+ 4802)*YIN(MY+ 4802)*ZIN(MZ+ 4802)           
     *          +   XIN(MX+ 7203)*YIN(MY+ 7203)*ZIN(MZ+ 7203))          
   41 CONTINUE                                                          
   42 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C          CODE FOR NROOTS=5                                            
C                                                                       
   50 CONTINUE                                                          
      DO 52 I = 1,IJ                                                    
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO 51 K = 1,MAX                                                   
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*                                
     *            ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )           
     *          +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401)           
     *          +   XIN(MX+ 4802)*YIN(MY+ 4802)*ZIN(MZ+ 4802)           
     *          +   XIN(MX+ 7203)*YIN(MY+ 7203)*ZIN(MZ+ 7203)           
     *          +   XIN(MX+ 9604)*YIN(MY+ 9604)*ZIN(MZ+ 9604))          
   51 CONTINUE                                                          
   52 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C          CODE FOR NROOTS=6                                            
C                                                                       
   60 CONTINUE                                                          
      DO 62 I = 1,IJ                                                    
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO 61 K = 1,MAX                                                   
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*                                
     *            ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )           
     *          +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401)           
     *          +   XIN(MX+ 4802)*YIN(MY+ 4802)*ZIN(MZ+ 4802)           
     *          +   XIN(MX+ 7203)*YIN(MY+ 7203)*ZIN(MZ+ 7203)           
     *          +   XIN(MX+ 9604)*YIN(MY+ 9604)*ZIN(MZ+ 9604)           
     *          +   XIN(MX+12005)*YIN(MY+12005)*ZIN(MZ+12005))          
   61 CONTINUE                                                          
   62 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C          CODE FOR NROOTS=7                                            
C                                                                       
   70 CONTINUE                                                          
      DO 72 I = 1,IJ                                                    
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO 71 K = 1,MAX                                                   
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*                                
     *            ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )           
     *          +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401)           
     *          +   XIN(MX+ 4802)*YIN(MY+ 4802)*ZIN(MZ+ 4802)           
     *          +   XIN(MX+ 7203)*YIN(MY+ 7203)*ZIN(MZ+ 7203)           
     *          +   XIN(MX+ 9604)*YIN(MY+ 9604)*ZIN(MZ+ 9604)           
     *          +   XIN(MX+12005)*YIN(MY+12005)*ZIN(MZ+12005)           
     *          +   XIN(MX+14406)*YIN(MY+14406)*ZIN(MZ+14406))          
   71 CONTINUE                                                          
   72 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C          CODE FOR NROOTS=8                                            
C                                                                       
   80 CONTINUE                                                          
      DO 82 I = 1,IJ                                                    
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO 81 K = 1,MAX                                                   
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*                                
     *            ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )           
     *          +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401)           
     *          +   XIN(MX+ 4802)*YIN(MY+ 4802)*ZIN(MZ+ 4802)           
     *          +   XIN(MX+ 7203)*YIN(MY+ 7203)*ZIN(MZ+ 7203)           
     *          +   XIN(MX+ 9604)*YIN(MY+ 9604)*ZIN(MZ+ 9604)           
     *          +   XIN(MX+12005)*YIN(MY+12005)*ZIN(MZ+12005)           
     *          +   XIN(MX+14406)*YIN(MY+14406)*ZIN(MZ+14406)           
     *          +   XIN(MX+16807)*YIN(MY+16807)*ZIN(MZ+16807))          
   81 CONTINUE                                                          
   82 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C          CODE FOR NROOTS=9                                            
C                                                                       
   90 CONTINUE                                                          
      DO 92 I = 1,IJ                                                    
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO 91 K = 1,MAX                                                   
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*                                
     *            ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )           
     *          +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401)           
     *          +   XIN(MX+ 4802)*YIN(MY+ 4802)*ZIN(MZ+ 4802)           
     *          +   XIN(MX+ 7203)*YIN(MY+ 7203)*ZIN(MZ+ 7203)           
     *          +   XIN(MX+ 9604)*YIN(MY+ 9604)*ZIN(MZ+ 9604)           
     *          +   XIN(MX+12005)*YIN(MY+12005)*ZIN(MZ+12005)           
     *          +   XIN(MX+14406)*YIN(MY+14406)*ZIN(MZ+14406)           
     *          +   XIN(MX+16807)*YIN(MY+16807)*ZIN(MZ+16807)           
     *          +   XIN(MX+19208)*YIN(MY+19208)*ZIN(MZ+19208))          
   91 CONTINUE                                                          
   92 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C          CODE FOR NROOTS=10                                           
C                                                                       
  100 CONTINUE                                                          
      DO 102 I = 1,IJ                                                   
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO 101 K = 1,MAX                                                  
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*                                
     *            ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )           
     *          +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401)           
     *          +   XIN(MX+ 4802)*YIN(MY+ 4802)*ZIN(MZ+ 4802)           
     *          +   XIN(MX+ 7203)*YIN(MY+ 7203)*ZIN(MZ+ 7203)           
     *          +   XIN(MX+ 9604)*YIN(MY+ 9604)*ZIN(MZ+ 9604)           
     *          +   XIN(MX+12005)*YIN(MY+12005)*ZIN(MZ+12005)           
     *          +   XIN(MX+14406)*YIN(MY+14406)*ZIN(MZ+14406)           
     *          +   XIN(MX+16807)*YIN(MY+16807)*ZIN(MZ+16807)           
     *          +   XIN(MX+19208)*YIN(MY+19208)*ZIN(MZ+19208)           
     *          +   XIN(MX+21609)*YIN(MY+21609)*ZIN(MZ+21609))          
  101 CONTINUE                                                          
  102 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C          CODE FOR NROOTS=11                                           
C                                                                       
  110 CONTINUE                                                          
      DO 112 I = 1,IJ                                                   
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO 111 K = 1,MAX                                                  
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*                                
     *            ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )           
     *          +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401)           
     *          +   XIN(MX+ 4802)*YIN(MY+ 4802)*ZIN(MZ+ 4802)           
     *          +   XIN(MX+ 7203)*YIN(MY+ 7203)*ZIN(MZ+ 7203)           
     *          +   XIN(MX+ 9604)*YIN(MY+ 9604)*ZIN(MZ+ 9604)           
     *          +   XIN(MX+12005)*YIN(MY+12005)*ZIN(MZ+12005)           
     *          +   XIN(MX+14406)*YIN(MY+14406)*ZIN(MZ+14406)           
     *          +   XIN(MX+16807)*YIN(MY+16807)*ZIN(MZ+16807)           
     *          +   XIN(MX+19208)*YIN(MY+19208)*ZIN(MZ+19208)           
     *          +   XIN(MX+21609)*YIN(MY+21609)*ZIN(MZ+21609)           
     *          +   XIN(MX+24010)*YIN(MY+24010)*ZIN(MZ+24010))          
  111 CONTINUE                                                          
  112 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C          CODE FOR NROOTS=12                                           
C                                                                       
  120 CONTINUE                                                          
      DO 122 I = 1,IJ                                                   
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO 121 K = 1,MAX                                                  
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*                                
     *            ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )           
     *          +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401)           
     *          +   XIN(MX+ 4802)*YIN(MY+ 4802)*ZIN(MZ+ 4802)           
     *          +   XIN(MX+ 7203)*YIN(MY+ 7203)*ZIN(MZ+ 7203)           
     *          +   XIN(MX+ 9604)*YIN(MY+ 9604)*ZIN(MZ+ 9604)           
     *          +   XIN(MX+12005)*YIN(MY+12005)*ZIN(MZ+12005)           
     *          +   XIN(MX+14406)*YIN(MY+14406)*ZIN(MZ+14406)           
     *          +   XIN(MX+16807)*YIN(MY+16807)*ZIN(MZ+16807)           
     *          +   XIN(MX+19208)*YIN(MY+19208)*ZIN(MZ+19208)           
     *          +   XIN(MX+21609)*YIN(MY+21609)*ZIN(MZ+21609)           
     *          +   XIN(MX+24010)*YIN(MY+24010)*ZIN(MZ+24010)           
     *          +   XIN(MX+26411)*YIN(MY+26411)*ZIN(MZ+26411))          
  121 CONTINUE                                                          
  122 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C          CODE FOR NROOTS=13                                           
C                                                                       
  130 DO 132 I = 1,IJ                                                   
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO 131 K = 1,MAX                                                  
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*                                
     *            ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )           
     *          +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401)           
     *          +   XIN(MX+ 4802)*YIN(MY+ 4802)*ZIN(MZ+ 4802)           
     *          +   XIN(MX+ 7203)*YIN(MY+ 7203)*ZIN(MZ+ 7203)           
     *          +   XIN(MX+ 9604)*YIN(MY+ 9604)*ZIN(MZ+ 9604)           
     *          +   XIN(MX+12005)*YIN(MY+12005)*ZIN(MZ+12005)           
     *          +   XIN(MX+14406)*YIN(MY+14406)*ZIN(MZ+14406)           
     *          +   XIN(MX+16807)*YIN(MY+16807)*ZIN(MZ+16807)           
     *          +   XIN(MX+19208)*YIN(MY+19208)*ZIN(MZ+19208)           
     *          +   XIN(MX+21609)*YIN(MY+21609)*ZIN(MZ+21609)           
     *          +   XIN(MX+24010)*YIN(MY+24010)*ZIN(MZ+24010)           
     *          +   XIN(MX+26411)*YIN(MY+26411)*ZIN(MZ+26411)           
     *          +   XIN(MX+28812)*YIN(MY+28812)*ZIN(MZ+28812))          
  131 CONTINUE                                                          
  132 CONTINUE                                                          
      RETURN                                                            
      END                                                               
C*MODULE INT2A   *DECK FORMS_gpu_4                                      
      SUBROUTINE FORMS_gpu_4(GHONDO,NROOTS,DKL,DIJ,XIN,YIN,ZIN,         
     * IJGT,IJX,IJY,IJZ,IK,KLGT,KLX,KLY,KLZ,IJ)                         
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      DIMENSION GHONDO(*)                                               
C                                                                       
      DIMENSION :: DKL(784),DIJ(784)                                    
      DIMENSION :: IJGT(784),IJX(784),IJY(784),IJZ(784),IK(784)         
      DIMENSION :: KLGT(784),KLX(784),KLY(784),KLZ(784)                 
      DIMENSION :: XIN(31213),YIN(31213),ZIN(31213)                     
C                                                                       
   40 CONTINUE                                                          
      DO 42 I = 1,IJ                                                    
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO 41 K = 1,MAX                                                   
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*                                
     *            ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )           
     *          +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401)           
     *          +   XIN(MX+ 4802)*YIN(MY+ 4802)*ZIN(MZ+ 4802)           
     *          +   XIN(MX+ 7203)*YIN(MY+ 7203)*ZIN(MZ+ 7203))          
   41 CONTINUE                                                          
   42 CONTINUE                                                          
C      write(*,*) "GHONDO in 4 is", GHONDO(N)                           
      RETURN                                                            
      END                                                               
C*MODULE INT2A   *DECK FORMS_gpu_5                                      
      SUBROUTINE FORMS_gpu_5(GHONDO,NROOTS,DKL,DIJ,XIN,YIN,ZIN,         
     * IJGT,IJX,IJY,IJZ,IK,KLGT,KLX,KLY,KLZ,IJ)                         
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      DIMENSION GHONDO(*)                                               
C                                                                       
      DIMENSION :: DKL(784),DIJ(784)                                    
      DIMENSION :: IJGT(784),IJX(784),IJY(784),IJZ(784),IK(784)         
      DIMENSION :: KLGT(784),KLX(784),KLY(784),KLZ(784)                 
      DIMENSION :: XIN(31213),YIN(31213),ZIN(31213)                     
C                                                                       
   50 CONTINUE                                                          
      DO 52 I = 1,IJ                                                    
      D1 = DIJ(I)                                                       
      NX = IJX(I)                                                       
      NY = IJY(I)                                                       
      NZ = IJZ(I)                                                       
      N1 = IJGT(I)                                                      
      MAX = IK(I)                                                       
      DO 51 K = 1,MAX                                                   
      MX = NX+KLX(K)                                                    
      MY = NY+KLY(K)                                                    
      MZ = NZ+KLZ(K)                                                    
      N = N1+KLGT(K)                                                    
      GHONDO(N) = GHONDO(N) + D1*DKL(K)*                                
     *            ( XIN(MX      )*YIN(MY      )*ZIN(MZ      )           
     *          +   XIN(MX+ 2401)*YIN(MY+ 2401)*ZIN(MZ+ 2401)           
     *          +   XIN(MX+ 4802)*YIN(MY+ 4802)*ZIN(MZ+ 4802)           
     *          +   XIN(MX+ 7203)*YIN(MY+ 7203)*ZIN(MZ+ 7203)           
     *          +   XIN(MX+ 9604)*YIN(MY+ 9604)*ZIN(MZ+ 9604))          
   51 CONTINUE                                                          
   52 CONTINUE                                                          
C      write(*,*) "GHONDO in 5 is", GHONDO(N)                           
      RETURN                                                            
      END                                                               
C*MODULE INT2A   *DECK XYZINT_gpu                                       
      SUBROUTINE XYZINT_gpu(XINT,YINT,ZINT,                             
     * I,K,NIMAX,NJMAX,NKMAX,NLMAX,NMAX,MMAX,                           
     * BP01,B00,B10,XCP00,XC00,YCP00,YC00,ZCP00,ZC00,F00,               
     * DXIJ,DYIJ,DZIJ,DXKL,DYKL,DZKL)                                   
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      LOGICAL N0,N1,M0,M1,FIRST1,FIRST2,FIRST3,FIRST4                   
      DIMENSION :: I(13),K(13)                                          
      INTEGER :: NIMAX,NJMAX,NKMAX,NLMAX,NMAX,MMAX                      
      double precision :: BP01,B00,B10,XCP00,XC00,YCP00,YC00,ZCP00,ZC00 
      double precision :: F00,DXIJ,DYIJ,DZIJ,DXKL,DYKL,DZKL             
      DIMENSION :: XINT(31213),YINT(31213),ZINT(31213)                  
C                                                                       
      PARAMETER (ZERO=0.0D+00)                                          
      PARAMETER (ONE=1.0D+00)                                           
C                                                                       
      N0 = NMAX .EQ. 0                                                  
      N1 = NMAX .LE. 1                                                  
      M0 = MMAX .EQ. 0                                                  
      M1 = MMAX .LE. 1                                                  
C                                                                       
C     ----- I(0,0) -----                                                
      I1 = I(1)                                                         
      XINT(I1) = ONE                                                    
      YINT(I1) = ONE                                                    
      ZINT(I1) = F00                                                    
C      IF (N0 .AND. M0) RETURN                                          
      I2 = I(2)                                                         
      K2 = K(2)                                                         
      CP10 = B00                                                        
C                                                                       
C     ----- I(1,0) -----                                                
C                                                                       
C      IF (.NOT. N0) THEN                                               
        XINT(I2) = XC00                                                 
        YINT(I2) = YC00                                                 
        ZINT(I2) = ZC00*F00                                             
        IF (M0) GO TO 120                                               
C      END IF                                                           
C                                                                       
C     ----- I(0,1) -----                                                
C                                                                       
      I3 = I1+K2                                                        
      XINT(I3) = XCP00                                                  
      YINT(I3) = YCP00                                                  
      ZINT(I3) = ZCP00*F00                                              
C                                                                       
C     ----- I(1,1) -----                                                
C                                                                       
      IF (.NOT. N0) THEN                                                
        I3 = I2+K2                                                      
        XINT(I3) = XCP00*XINT(I2)+CP10                                  
        YINT(I3) = YCP00*YINT(I2)+CP10                                  
        ZINT(I3) = ZCP00*ZINT(I2)+CP10*F00                              
      END IF                                                            
C                                                                       
  120 CONTINUE                                                          
      IF (.NOT. N1) THEN                                                
        C10 = ZERO                                                      
        I3 = I1                                                         
        I4 = I2                                                         
        DO 160 N = 2,NMAX                                               
          C10 = C10+B10                                                 
C                                                                       
C     ----- I(N,0) -----                                                
C                                                                       
          I5 = I(N+1)                                                   
          XINT(I5) = C10*XINT(I3)+XC00*XINT(I4)                         
          YINT(I5) = C10*YINT(I3)+YC00*YINT(I4)                         
          ZINT(I5) = C10*ZINT(I3)+ZC00*ZINT(I4)                         
          IF ( .NOT. M0) THEN                                           
            CP10 = CP10+B00                                             
C                                                                       
C     ----- I(N,1) -----                                                
C                                                                       
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
C                                                                       
C     ----- I(0,M) -----                                                
C                                                                       
          I5 = I1+K(M+1)                                                
          XINT(I5) = CP01*XINT(I3)+XCP00*XINT(I4)                       
          YINT(I5) = CP01*YINT(I3)+YCP00*YINT(I4)                       
          ZINT(I5) = CP01*ZINT(I3)+ZCP00*ZINT(I4)                       
C                                                                       
C     ----- I(1,M) -----                                                
C                                                                       
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
C                                                                       
C     ----- I(N,M) -----                                                
C                                                                       
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
            XINT(I5+K4) = C10*XINT(I3+K4)+XC00*XINT(I4+K4)              
     *                    +C01*XINT(I4+K3)                              
            YINT(I5+K4) = C10*YINT(I3+K4)+YC00*YINT(I4+K4)              
     *                    +C01*YINT(I4+K3)                              
            ZINT(I5+K4) = C10*ZINT(I3+K4)+ZC00*ZINT(I4+K4)              
     *                    +C01*ZINT(I4+K3)                              
            C10 = C10+B10                                               
            I3 = I4                                                     
            I4 = I5                                                     
  260     CONTINUE                                                      
          K3 = K4                                                       
  280   CONTINUE                                                        
      END IF                                                            
C                                                                       
C     ----- I(NI,NJ,M) -----                                            
C                                                                       
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
C                                                                       
C     ----- I(NI,NJ,NK,NL) -----                                        
C                                                                       
C      write(*,*) "NLMAX is", NLMAX                                     
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
C                                                                       
      RETURN                                                            
      END                                                               
C C*MODULE INT2A   *DECK XYZINT_gpu_2                                   
C       SUBROUTINE XYZINT_gpu_2(XINT,YINT,ZINT,                         
C      * I,K,NIMAX,NJMAX,NKMAX,NLMAX,NMAX,MMAX,                         
C      * BP01,B00,B10,XCP00,XC00,YCP00,YC00,ZCP00,ZC00,F00,             
C      * DXIJ,DYIJ,DZIJ,DXKL,DYKL,DZKL)                                 
C C                                                                     
C       IMPLICIT DOUBLE PRECISION (A-H,O-Z)                             
C C                                                                     
C       LOGICAL N0,N1,M0,M1,FIRST1,FIRST2,FIRST3,FIRST4                 
C       DIMENSION :: I(13),K(13)                                        
C       INTEGER :: NIMAX,NJMAX,NKMAX,NLMAX,NMAX,MMAX                    
C       double precision :: BP01,B00,B10,XCP00,XC00,YCP00,YC00,ZCP00,ZC0
C       double precision :: F00,DXIJ,DYIJ,DZIJ,DXKL,DYKL,DZKL           
C       DIMENSION :: XINT(31213),YINT(31213),ZINT(31213)                
C C                                                                     
C       PARAMETER (ZERO=0.0D+00)                                        
C       PARAMETER (ONE=1.0D+00)                                         
C C                                                                     
C       N0 = NMAX .EQ. 0                                                
C       N1 = NMAX .LE. 1                                                
C       M0 = MMAX .EQ. 0                                                
C       M1 = MMAX .LE. 1                                                
C C                                                                     
C C     ----- I(0,0) -----                                              
C C                                                                     
C       I1 = I(1)                                                       
C       XINT(I1) = ONE                                                  
C       YINT(I1) = ONE                                                  
C       ZINT(I1) = F00                                                  
C       IF (N0 .AND. M0) RETURN                                         
C       I2 = I(2)                                                       
C       K2 = K(2)                                                       
C       CP10 = B00                                                      
C C                                                                     
C C     ----- I(1,0) -----                                              
C C                                                                     
C       IF (.NOT. N0) THEN                                              
C C        write(*,*) "here1"                                           
C         XINT(I2) = XC00                                               
C         YINT(I2) = YC00                                               
C         ZINT(I2) = ZC00*F00                                           
C         IF (M0) GO TO 120                                             
C       END IF                                                          
C C                                                                     
C C     ----- I(0,1) -----                                              
C C                                                                     
C       I3 = I1+K2                                                      
C       XINT(I3) = XCP00                                                
C       YINT(I3) = YCP00                                                
C       ZINT(I3) = ZCP00*F00                                            
C C                                                                     
C C     ----- I(1,1) -----                                              
C C                                                                     
C       IF (.NOT. N0) THEN                                              
C C        write(*,*) "here2"                                           
C         I3 = I2+K2                                                    
C         XINT(I3) = XCP00*XINT(I2)+CP10                                
C         YINT(I3) = YCP00*YINT(I2)+CP10                                
C         ZINT(I3) = ZCP00*ZINT(I2)+CP10*F00                            
C       END IF                                                          
C C                                                                     
C   120 CONTINUE                                                        
C       IF (.NOT. N1) THEN                                              
C C        write(*,*) "here3"                                           
C         C10 = ZERO                                                    
C         I3 = I1                                                       
C         I4 = I2                                                       
C         DO 160 N = 2,NMAX                                             
C           C10 = C10+B10                                               
C C                                                                     
C C     ----- I(N,0) -----                                              
C C                                                                     
C           I5 = I(N+1)                                                 
C           XINT(I5) = C10*XINT(I3)+XC00*XINT(I4)                       
C           YINT(I5) = C10*YINT(I3)+YC00*YINT(I4)                       
C           ZINT(I5) = C10*ZINT(I3)+ZC00*ZINT(I4)                       
C           IF ( .NOT. M0) THEN                                         
C             CP10 = CP10+B00                                           
C C                                                                     
C C     ----- I(N,1) -----                                              
C C                                                                     
C             I3 = I5+K2                                                
C             XINT(I3) = XCP00*XINT(I5)+CP10*XINT(I4)                   
C             YINT(I3) = YCP00*YINT(I5)+CP10*YINT(I4)                   
C             ZINT(I3) = ZCP00*ZINT(I5)+CP10*ZINT(I4)                   
C           END IF                                                      
C           I3 = I4                                                     
C           I4 = I5                                                     
C   160     CONTINUE                                                    
C       END IF                                                          
C       IF ( .NOT. M1) THEN                                             
C C        write(*,*) "here4"                                           
C         CP01 = ZERO                                                   
C         C01 = B00                                                     
C         I3 = I1                                                       
C         I4 = I1+K2                                                    
C         DO 220 M = 2,MMAX                                             
C           CP01 = CP01+BP01                                            
C C                                                                     
C C     ----- I(0,M) -----                                              
C C                                                                     
C           I5 = I1+K(M+1)                                              
C           XINT(I5) = CP01*XINT(I3)+XCP00*XINT(I4)                     
C           YINT(I5) = CP01*YINT(I3)+YCP00*YINT(I4)                     
C           ZINT(I5) = CP01*ZINT(I3)+ZCP00*ZINT(I4)                     
C C                                                                     
C C     ----- I(1,M) -----                                              
C C                                                                     
C           IF (.NOT. N0) THEN                                          
C             C01 = C01+B00                                             
C             I3 = I2+K(M+1)                                            
C             XINT(I3) = XC00*XINT(I5)+C01*XINT(I4)                     
C             YINT(I3) = YC00*YINT(I5)+C01*YINT(I4)                     
C             ZINT(I3) = ZC00*ZINT(I5)+C01*ZINT(I4)                     
C           END IF                                                      
C           I3 = I4                                                     
C           I4 = I5                                                     
C   220   CONTINUE                                                      
C       END IF                                                          
C C                                                                     
C C     ----- I(N,M) -----                                              
C C                                                                     
C       IF (.NOT. N1 .AND. .NOT. M1) THEN                               
C C        write(*,*) "here5"                                           
C         C01 = B00                                                     
C         K3 = K2                                                       
C         DO 280 M = 2,MMAX                                             
C           K4 = K(M+1)                                                 
C           C01 = C01+B00                                               
C           I3 = I1                                                     
C           I4 = I2                                                     
C           C10 = B10                                                   
C           DO 260 N = 2,NMAX                                           
C             I5 = I(N+1)                                               
C             XINT(I5+K4) = C10*XINT(I3+K4)+XC00*XINT(I4+K4)            
C      *                    +C01*XINT(I4+K3)                            
C             YINT(I5+K4) = C10*YINT(I3+K4)+YC00*YINT(I4+K4)            
C      *                    +C01*YINT(I4+K3)                            
C             ZINT(I5+K4) = C10*ZINT(I3+K4)+ZC00*ZINT(I4+K4)            
C      *                    +C01*ZINT(I4+K3)                            
C             C10 = C10+B10                                             
C             I3 = I4                                                   
C             I4 = I5                                                   
C   260     CONTINUE                                                    
C           K3 = K4                                                     
C   280   CONTINUE                                                      
C       END IF                                                          
C C                                                                     
C C     ----- I(NI,NJ,M) -----                                          
C C                                                                     
C       IF (NJMAX .GT. 0) THEN                                          
C C        write(*,*) "here6"                                           
C         M = 0                                                         
C         I5 = I(NMAX+1)                                                
C         FIRST1 = .TRUE.                                               
C         DO 430 WHILE (FIRST1 .OR. M .LE. MMAX)                        
C           MIN = NIMAX                                                 
C           KM = K(M+1)                                                 
C           FIRST2 = .TRUE.                                             
C           DO 360 WHILE (FIRST2 .OR. MIN .LT. NMAX)                    
C             N = NMAX                                                  
C             I3 = I5+KM                                                
C             FIRST3 = .TRUE.                                           
C             DO 340 WHILE (FIRST3 .OR. N .GT. MIN)                     
C               I4 = I(N)+KM                                            
C               XINT(I3) = XINT(I3)+DXIJ*XINT(I4)                       
C               YINT(I3) = YINT(I3)+DYIJ*YINT(I4)                       
C               ZINT(I3) = ZINT(I3)+DZIJ*ZINT(I4)                       
C               I3 = I4                                                 
C               N = N-1                                                 
C               FIRST3 = .FALSE.                                        
C   340       END DO                                                    
C             MIN = MIN+1                                               
C             FIRST2 = .FALSE.                                          
C   360     END DO                                                      
C           IF (NIMAX .GT. 0) THEN                                      
C C            write(*,*) "here7"                                       
C             I3 = 49+KM+I1                                             
C             DO 400 NJ = 1,NJMAX                                       
C               I4 = I3                                                 
C               DO 380 NI = 1,NIMAX                                     
C                 XINT(I4) = XINT(I4+294)+DXIJ*XINT(I4-49)              
C                 YINT(I4) = YINT(I4+294)+DYIJ*YINT(I4-49)              
C                 ZINT(I4) = ZINT(I4+294)+DZIJ*ZINT(I4-49)              
C                 I4 = I4+343                                           
C   380         CONTINUE                                                
C               I3 = I3+49                                              
C   400       CONTINUE                                                  
C           END IF                                                      
C           M = M+1                                                     
C           FIRST1 = .FALSE.                                            
C   430   END DO                                                        
C       END IF                                                          
C C                                                                     
C C     ----- I(NI,NJ,NK,NL) -----                                      
C C                                                                     
C       IF (NLMAX .GT. 0) THEN                                          
C C        write(*,*) "here8"                                           
C         I5 = K(MMAX+1)                                                
C         IA = I1                                                       
C         NI = 0                                                        
C         FIRST4 = .TRUE.                                               
C         DO 580 WHILE (FIRST4 .OR. NI .LE. NIMAX)                      
C           NJ = 0                                                      
C           IB = IA                                                     
C           FIRST1 = .TRUE.                                             
C           DO 570 WHILE (FIRST1 .OR. NJ .LE. NJMAX)                    
C             MIN = NKMAX                                               
C             FIRST2 = .TRUE.                                           
C             DO 530 WHILE (FIRST2 .OR. MIN .LT. MMAX)                  
C               M = MMAX                                                
C               I3 = IB+I5                                              
C               FIRST3 = .TRUE.                                         
C               DO 520 WHILE (FIRST3 .OR. M .GT. MIN)                   
C                 I4 = IB+K(M)                                          
C                 XINT(I3) = XINT(I3)+DXKL*XINT(I4)                     
C                 YINT(I3) = YINT(I3)+DYKL*YINT(I4)                     
C                 ZINT(I3) = ZINT(I3)+DZKL*ZINT(I4)                     
C                 I3 = I4                                               
C                 M = M-1                                               
C                 FIRST3 = .FALSE.                                      
C   520         END DO                                                  
C               MIN = MIN+1                                             
C               FIRST2 = .FALSE.                                        
C   530       END DO                                                    
C             IF (NKMAX .GT. 0) THEN                                    
C               I3 = IB+1                                               
C               DO 560 NL = 1,NLMAX                                     
C                 I4 = I3                                               
C                 DO 540 NK = 1,NKMAX                                   
C                   XINT(I4) = XINT(I4+6)+DXKL*XINT(I4-1)               
C                   YINT(I4) = YINT(I4+6)+DYKL*YINT(I4-1)               
C                   ZINT(I4) = ZINT(I4+6)+DZKL*ZINT(I4-1)               
C                   I4 = I4+7                                           
C   540           END DO                                                
C               I3 = I3+1                                               
C   560         END DO                                                  
C             END IF                                                    
C             NJ = NJ+1                                                 
C             IB = IB+49                                                
C             FIRST1 = .FALSE.                                          
C   570     END DO                                                      
C           NI = NI+1                                                   
C           IA = IA+343                                                 
C           FIRST4 = .FALSE.                                            
C   580   END DO                                                        
C       END IF                                                          
C C                                                                     
C       RETURN                                                          
C       END                                                             
C*MODULE INT2A   *DECK ZQOUT_gpu                                        
      SUBROUTINE ZQOUT_gpu(GHONDO,                                      
     * IJGT,IJX,IJY,IJZ,IK,KLGT,KLX,KLY,KLZ,                            
     * IANDJ,KANDL,SAME,                                                
     * mini,maxi,minj,maxj,mink,maxk,minl,maxl)                         
C                                                                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C                                                                       
      DIMENSION GHONDO(*)                                               
C                                                                       
      LOGICAL IANDJ,KANDL,SAME                                          
      DIMENSION :: IJGT(784),IJX(784),IJY(784),IJZ(784),IK(784)         
      DIMENSION :: KLGT(784),KLX(784),KLY(784),KLZ(784)                 
      INTEGER :: MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL                
C                                                                       
      PARAMETER (ZERO=0.0D+00)                                          
C                                                                       
C     ----- ZERO HONDO CONVENTIONAL INTEGRAL OUTPUT REGION -----        
C                                                                       
!!$omp declare target                                                   
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
                  IF (SAME .AND. KLN .GT. IJN) GO TO 240                
                  NN = N1+KLGT(KLN)                                     
                  GHONDO(NN) = ZERO                                     
  200          CONTINUE                                                 
  220       CONTINUE                                                    
  240    CONTINUE                                                       
  260 CONTINUE                                                          
      RETURN                                                            
      END                                                               
C C*MODULE INT2A   *DECK ZQOUT_gpu_2                                    
C       SUBROUTINE ZQOUT_gpu_2(GHONDO,                                  
C      * IJGT,KLGT,                                                     
C      * IANDJ,KANDL,SAME,                                              
C      * mini,maxi,minj,maxj,mink,maxk,minl,maxl)                       
C C                                                                     
C       IMPLICIT DOUBLE PRECISION(A-H,O-Z)                              
C C                                                                     
C       DIMENSION GHONDO(*)                                             
C C                                                                     
C       LOGICAL IANDJ,KANDL,SAME                                        
C       DIMENSION :: IJGT(784),KLGT(784)                                
C       INTEGER :: MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL              
C C                                                                     
C       PARAMETER (ZERO=0.0D+00)                                        
C C                                                                     
C C     ----- ZERO HONDO CONVENTIONAL INTEGRAL OUTPUT REGION -----      
C C                                                                     
C !!$omp declare target                                                 
C       IJN = 0                                                         
C       JMAX = MAXJ                                                     
C       DO 260 I = MINI,MAXI                                            
C          IF (IANDJ) JMAX = I                                          
C          DO 240 J = MINJ,JMAX                                         
C             IJN = IJN+1                                               
C             N1 = IJGT(IJN)                                            
C             LMAX = MAXL                                               
C             KLN = 0                                                   
C             DO 220 K =  MINK,MAXK                                     
C                IF (KANDL) LMAX = K                                    
C                DO 200 L = MINL,LMAX                                   
C                   KLN = KLN+1                                         
C                   !IF (SAME .AND. KLN .GT. IJN) GO TO 240             
C                   NN = N1+KLGT(KLN)                                   
C                   GHONDO(NN) = ZERO                                   
C   200          CONTINUE                                               
C   220       CONTINUE                                                  
C   240    CONTINUE                                                     
C   260 CONTINUE                                                        
C       RETURN                                                          
C       END                                                             
C*MODULE INT2A   *DECK IJPRIM_gpu                                       
      SUBROUTINE IJPRIM_gpu(DDIJ,A,R,X1,Y1,Z1,IJD,IANDJ,KANDL,SAME,     
     * LIT,mini,maxi,minj,maxj,mink,maxk,minl,maxl,NIJ,                 
     * ag,csa,cpa,cda,cfa,cga,cha,cia,                                  
     * bg,csb,cpb,cdb,cfb,cgb,chb,cib,                                  
     * cg,csc,cpc,cdc,cfc,cgc,chc,cic,                                  
     * dg,csd,cpd,cdd,cfd,cgd,chd,cid,                                  
     * XI,YI,ZI,XJ,YJ,ZJ,RRI,XK,YK,ZK,XL,YL,ZL,RRK,                     
     * NGA,NGB,NGC,NGD)                                                 
                                                                        
      use mx_limits, only: mxgsh,mxg2                                   
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      LOGICAL IANDJ,KANDL,SAME,OUT,NORM                                 
      DIMENSION DDIJ(49*MXG2)                                           
      DIMENSION :: A(MXG2),R(MXG2),X1(MXG2),Y1(MXG2),Z1(MXG2),IJD(784)  
      DIMENSION :: AG(MXGSH),CSA(MXGSH),CPA(MXGSH),CDA(MXGSH)           
      DIMENSION :: CFA(MXGSH),CGA(MXGSH),CHA(MXGSH),CIA(MXGSH)          
      DIMENSION :: BG(MXGSH),CSB(MXGSH),CPB(MXGSH),CDB(MXGSH)           
      DIMENSION :: CFB(MXGSH),CGB(MXGSH),CHB(MXGSH),CIB(MXGSH)          
      DIMENSION :: CG(MXGSH),CSC(MXGSH),CPC(MXGSH),CDC(MXGSH)           
      DIMENSION :: CFC(MXGSH),CGC(MXGSH),CHC(MXGSH),CIC(MXGSH)          
      DIMENSION :: DG(MXGSH),CSD(MXGSH),CPD(MXGSH),CDD(MXGSH)           
      DIMENSION :: CFD(MXGSH),CGD(MXGSH),CHD(MXGSH),CID(MXGSH)          
      double precision :: XI,YI,ZI,XJ,YJ,ZJ,RRI,XK,YK,ZK,XL,YL,ZL,RRK   
      integer :: NGA,NGB,NGC,NGD                                        
      double precision :: qq4                                           
C      INTEGER :: NROOTS                                                
      INTEGER :: LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL                    
      INTEGER :: MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL                
      INTEGER :: NIJ,IJ,KL,IJKL                                         
      COMMON /SHLT  / TOL,CUTOFF,ICOUNT,OUT                             
C                                                                       
      PARAMETER (SQRT3=1.73205080756888D+00)                            
      PARAMETER (SQRT5=2.23606797749979D+00)                            
      PARAMETER (SQRT7=2.64575131106459D+00)                            
      PARAMETER (SQRT9=3.0D+00)                                         
      PARAMETER (SQRT11=3.3166247903553998D+00)                         
      PARAMETER (ZERO=0.0D+00)                                          
      PARAMETER (ONE=1.0D+00)                                           
!!$omp declare target                                                   
C                                                                       
      NORM = NORMF .NE. 1 .OR. NORMP .NE. 1                             
      MAX = MAXJ                                                        
      N = 0                                                             
      NN = 0                                                            
      NM = -2**20                                                       
      DO 180 I = MINI,MAXI                                              
         GO TO (100,100,120,120,100,120,120,100,120,120,                
     1          100,120,120,100,120,120,120,120,120,100,                
     1          100,120,120,100,120,120,120,120,120,100,                
     1          120,120,100,120,120,                                    
     1          100,120,120,100,120,120,120,120,120,100,                
     1          120,120,120,120,120,100,120,120,100,120,                
     1          120,                                                    
     1          100,120,120,100,120,120,120,120,120,100,                
     1          120,120,120,120,120,100,120,120,100,120,                
     1          120,100,120,120,120,120,120,100),I                      
  100    NM = NN                                                        
  120    NN = NM                                                        
         IF (IANDJ) MAX = I                                             
         DO 170 J = MINJ,MAX                                            
            GO TO (140,140,160,160,140,160,160,140,160,160,             
     1             140,160,160,140,160,160,160,160,160,140,             
     1             140,160,160,140,160,160,160,160,160,140,             
     1             160,160,140,160,160,                                 
     1             140,160,160,140,160,160,160,160,160,140,             
     1             160,160,160,160,160,140,160,160,140,160,             
     1             160,                                                 
     1             140,160,160,140,160,160,160,160,160,140,             
     1             160,160,160,160,160,140,160,160,140,160,             
     1             160,140,160,160,160,160,160,140),J                   
  140       NN = NN+1                                                   
  160       N = N+1                                                     
            IJD(N) = NN                                                 
  170    CONTINUE                                                       
  180 CONTINUE                                                          
C                                                                       
C     ----- I PRIMITIVE                                                 
C                                                                       
      NIJ = 0                                                           
      JBMAX = NGB                                                       
      DO 540 IA = 1,NGA                                                 
         AI = AG(IA)                                                    
         ARRI = AI*RRI                                                  
         AXI = AI*XI                                                    
         AYI = AI*YI                                                    
         AZI = AI*ZI                                                    
         CSI = CSA(IA)                                                  
         CPI = CPA(IA)                                                  
         CDI = CDA(IA)                                                  
         CFI = CFA(IA)                                                  
         CGI = CGA(IA)                                                  
         CHI = CHA(IA)                                                  
         CII = CIA(IA)                                                  
C                                                                       
C        ----- J PRIMITIVE                                              
C                                                                       
         IF (IANDJ) JBMAX = IA                                          
         DO 520 JB = 1,JBMAX                                            
            AJ = BG(JB)                                                 
            AA = AI+AJ                                                  
            AAINV = ONE/AA                                              
            DUM = AJ*ARRI*AAINV                                         
            IF (DUM .GT. TOL) GO TO 520                                 
            CSJ = CSB(JB)                                               
            CPJ = CPB(JB)                                               
            CDJ = CDB(JB)                                               
            CFJ = CFB(JB)                                               
            CGJ = CGB(JB)                                               
            CHJ = CHB(JB)                                               
            CIJ = CIB(JB)                                               
            NM = 49*NIJ                                                 
            NN = NM                                                     
            NIJ = NIJ+1                                                 
            R(NIJ) = DUM                                                
            A(NIJ) = AA                                                 
            X1(NIJ) = (AXI+AJ*XJ)*AAINV                                 
            Y1(NIJ) = (AYI+AJ*YJ)*AAINV                                 
            Z1(NIJ) = (AZI+AJ*ZJ)*AAINV                                 
C                                                                       
C           ----- DENSITY FACTOR                                        
C                                                                       
            DUM1 = ZERO                                                 
            DUM2 = ZERO                                                 
            DO 420 I = MINI,MAXI                                        
               GO TO (200,220,420,420,240,420,420,260,420,420,          
     1                261,420,420,262,420,420,420,420,420,263,          
     1                264,420,420,265,420,420,420,420,420,266,          
     1                420,420,267,420,420,                              
     1                268,420,420,269,420,420,420,420,420,270,          
     1                420,420,420,420,420,271,420,420,272,420,          
     1                420,                                              
     1                273,420,420,274,420,420,420,420,420,275,          
     1                420,420,420,420,420,276,420,420,277,420,          
     1                420,278,420,420,420,420,420,279),I                
  200          DUM1 = CSI*AAINV                                         
               GO TO 280                                                
  220          DUM1 = CPI*AAINV                                         
               GO TO 280                                                
  240          DUM1 = CDI*AAINV                                         
               GO TO 280                                                
  260          DUM1 = DUM1*SQRT3                                        
               GO TO 280                                                
  261          DUM1 = CFI*AAINV                                         
               GO TO 280                                                
  262          DUM1 = DUM1*SQRT5                                        
               GO TO 280                                                
  263          DUM1 = DUM1*SQRT3                                        
               GO TO 280                                                
  264          DUM1 = CGI*AAINV                                         
               GO TO 280                                                
  265          DUM1 = DUM1*SQRT7                                        
               GO TO 280                                                
  266          DUM1 = DUM1*SQRT5/SQRT3                                  
               GO TO 280                                                
  267          DUM1 = DUM1*SQRT3                                        
               GO TO 280                                                
  268          DUM1 = CHI*AAINV                                         
               GO TO 280                                                
  269          DUM1 = DUM1*SQRT9                                        
               GO TO 280                                                
  270          DUM1 = DUM1*SQRT7/SQRT3                                  
               GO TO 280                                                
  271          DUM1 = DUM1*SQRT3                                        
               GO TO 280                                                
  272          DUM1 = DUM1*SQRT5/SQRT3                                  
               GO TO 280                                                
  273          DUM1 = CII*AAINV                                         
               GO TO 280                                                
  274          DUM1 = DUM1*SQRT11                                       
               GO TO 280                                                
  275          DUM1 = DUM1*SQRT3                                        
               GO TO 280                                                
  276          DUM1 = DUM1*SQRT3                                        
               GO TO 280                                                
  277          DUM1 = DUM1*SQRT7/(SQRT5*SQRT3)                          
               GO TO 280                                                
  278          DUM1 = DUM1*SQRT5                                        
               GO TO 280                                                
  279          DUM1 = DUM1*SQRT5/SQRT3                                  
C                                                                       
  280          IF (IANDJ) MAX = I                                       
               DO 400 J = MINJ,MAX                                      
                  GO TO (300,320,400,400,340,400,400,360,400,400,       
     1                   361,400,400,362,400,400,400,400,400,363,       
     1                   364,400,400,365,400,400,400,400,400,366,       
     1                   400,400,367,400,400,                           
     1                   368,400,400,369,400,400,400,400,400,370,       
     1                   400,400,400,400,400,371,400,400,372,400,       
     1                   400,                                           
     1                   373,400,400,374,400,400,400,400,400,375,       
     1                   400,400,400,400,400,376,400,400,377,400,       
     1                   400,378,400,400,400,400,400,379),J             
  300             DUM2 = DUM1*CSJ                                       
                  GO TO 380                                             
  320             DUM2 = DUM1*CPJ                                       
                  GO TO 380                                             
  340             DUM2 = DUM1*CDJ                                       
                  GO TO 380                                             
  360             DUM2 = DUM2*SQRT3                                     
                  GO TO 380                                             
  361             DUM2 = DUM1*CFJ                                       
                  GO TO 380                                             
  362             DUM2 = DUM2*SQRT5                                     
                  GO TO 380                                             
  363             DUM2 = DUM2*SQRT3                                     
                  GO TO 380                                             
  364             DUM2 = DUM1*CGJ                                       
                  GO TO 380                                             
  365             DUM2 = DUM2*SQRT7                                     
                  GO TO 380                                             
  366             DUM2 = DUM2*SQRT5/SQRT3                               
                  GO TO 380                                             
  367             DUM2 = DUM2*SQRT3                                     
                  GO TO 380                                             
  368             DUM2 = DUM1*CHJ                                       
                  GO TO 380                                             
  369             DUM2 = DUM2*SQRT9                                     
                  GO TO 380                                             
  370             DUM2 = DUM2*SQRT7/SQRT3                               
                  GO TO 380                                             
  371             DUM2 = DUM2*SQRT3                                     
                  GO TO 380                                             
  372             DUM2 = DUM2*SQRT5/SQRT3                               
                  GO TO 380                                             
  373             DUM2 = DUM1*CIJ                                       
                  GO TO 380                                             
  374             DUM2 = DUM2*SQRT11                                    
                  GO TO 380                                             
  375             DUM2 = DUM2*SQRT3                                     
                  GO TO 380                                             
  376             DUM2 = DUM2*SQRT3                                     
                  GO TO 380                                             
  377             DUM2 = DUM2*SQRT7/(SQRT5*SQRT3)                       
                  GO TO 380                                             
  378             DUM2 = DUM2*SQRT5                                     
                  GO TO 380                                             
  379             DUM2 = DUM2*SQRT5/SQRT3                               
C                                                                       
  380             NN = NN+1                                             
                  DDIJ(NN) = DUM2                                       
C                  write(*,*) "DDIJ is", DDIJ(NN)                       
  400          CONTINUE                                                 
  420       CONTINUE                                                    
            IF ( .NOT. IANDJ) GO TO 520                                 
            IF (IA .EQ. JB) GO TO 520                                   
            GO TO (500,440,460,455,450,445,444),LIT                     
  440       IF (MINI .EQ. 2) GO TO 500                                  
            DDIJ(NM+2) = DDIJ(NM+2)+CSI*CPJ*AAINV                       
            GO TO 480                                                   
  444       DDIJ(NM+28) = DDIJ(NM+28)+DDIJ(NM+28)                       
            DDIJ(NM+27) = DDIJ(NM+27)+DDIJ(NM+27)                       
            DDIJ(NM+26) = DDIJ(NM+26)+DDIJ(NM+26)                       
            DDIJ(NM+25) = DDIJ(NM+25)+DDIJ(NM+25)                       
            DDIJ(NM+24) = DDIJ(NM+24)+DDIJ(NM+24)                       
            DDIJ(NM+23) = DDIJ(NM+23)+DDIJ(NM+23)                       
            DDIJ(NM+22) = DDIJ(NM+22)+DDIJ(NM+22)                       
            DDIJ(NM+21) = DDIJ(NM+21)+DDIJ(NM+21)                       
            DDIJ(NM+20) = DDIJ(NM+20)+DDIJ(NM+20)                       
            DDIJ(NM+19) = DDIJ(NM+19)+DDIJ(NM+19)                       
            DDIJ(NM+18) = DDIJ(NM+18)+DDIJ(NM+18)                       
            DDIJ(NM+17) = DDIJ(NM+17)+DDIJ(NM+17)                       
            DDIJ(NM+16) = DDIJ(NM+16)+DDIJ(NM+16)                       
  445       DDIJ(NM+15) = DDIJ(NM+15)+DDIJ(NM+15)                       
            DDIJ(NM+14) = DDIJ(NM+14)+DDIJ(NM+14)                       
            DDIJ(NM+13) = DDIJ(NM+13)+DDIJ(NM+13)                       
            DDIJ(NM+12) = DDIJ(NM+12)+DDIJ(NM+12)                       
            DDIJ(NM+11) = DDIJ(NM+11)+DDIJ(NM+11)                       
  450       DDIJ(NM+10) = DDIJ(NM+10)+DDIJ(NM+10)                       
            DDIJ(NM+9) = DDIJ(NM+9)+DDIJ(NM+9)                          
            DDIJ(NM+8) = DDIJ(NM+8)+DDIJ(NM+8)                          
            DDIJ(NM+7) = DDIJ(NM+7)+DDIJ(NM+7)                          
  455       DDIJ(NM+6) = DDIJ(NM+6)+DDIJ(NM+6)                          
            DDIJ(NM+5) = DDIJ(NM+5)+DDIJ(NM+5)                          
            DDIJ(NM+4) = DDIJ(NM+4)+DDIJ(NM+4)                          
  460       DDIJ(NM+2) = DDIJ(NM+2)+DDIJ(NM+2)                          
  480       DDIJ(NM+3) = DDIJ(NM+3)+DDIJ(NM+3)                          
  500       DDIJ(NM+1) = DDIJ(NM+1)+DDIJ(NM+1)                          
  520    CONTINUE                                                       
  540 CONTINUE                                                          
      RETURN                                                            
      END                                                               
C C*MODULE INT2A   *DECK IJPRIM_gpu_2                                   
C       SUBROUTINE IJPRIM_gpu_2(DDIJ,A,R,X1,Y1,Z1,IJD,IANDJ,KANDL,SAME, 
C      * LIT,mini,maxi,minj,maxj,mink,maxk,minl,maxl,NIJ,               
C      * ag,csa,cpa,cda,cfa,cga,cha,cia,                                
C      * bg,csb,cpb,cdb,cfb,cgb,chb,cib,                                
C      * cg,csc,cpc,cdc,cfc,cgc,chc,cic,                                
C      * dg,csd,cpd,cdd,cfd,cgd,chd,cid,                                
C      * XI,YI,ZI,XJ,YJ,ZJ,RRI,XK,YK,ZK,XL,YL,ZL,RRK,                   
C      * NGA,NGB,NGC,NGD)                                               
                                                                        
C       use mx_limits, only: mxgsh,mxg2                                 
C C                                                                     
C       IMPLICIT DOUBLE PRECISION (A-H,O-Z)                             
C C                                                                     
C       LOGICAL IANDJ,KANDL,SAME,OUT,NORM                               
C       DIMENSION DDIJ(49*MXG2)                                         
C       DIMENSION :: A(MXG2),R(MXG2),X1(MXG2),Y1(MXG2),Z1(MXG2),IJD(784)
C       DIMENSION :: AG(MXGSH),CSA(MXGSH),CPA(MXGSH),CDA(MXGSH)         
C       DIMENSION :: CFA(MXGSH),CGA(MXGSH),CHA(MXGSH),CIA(MXGSH)        
C       DIMENSION :: BG(MXGSH),CSB(MXGSH),CPB(MXGSH),CDB(MXGSH)         
C       DIMENSION :: CFB(MXGSH),CGB(MXGSH),CHB(MXGSH),CIB(MXGSH)        
C       DIMENSION :: CG(MXGSH),CSC(MXGSH),CPC(MXGSH),CDC(MXGSH)         
C       DIMENSION :: CFC(MXGSH),CGC(MXGSH),CHC(MXGSH),CIC(MXGSH)        
C       DIMENSION :: DG(MXGSH),CSD(MXGSH),CPD(MXGSH),CDD(MXGSH)         
C       DIMENSION :: CFD(MXGSH),CGD(MXGSH),CHD(MXGSH),CID(MXGSH)        
C       double precision :: XI,YI,ZI,XJ,YJ,ZJ,RRI,XK,YK,ZK,XL,YL,ZL,RRK 
C       integer :: NGA,NGB,NGC,NGD                                      
C       double precision :: qq4                                         
C C      INTEGER :: NROOTS                                              
C       INTEGER :: LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL                  
C       INTEGER :: MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL              
C       INTEGER :: NIJ,IJ,KL,IJKL                                       
C       COMMON /SHLT  / TOL,CUTOFF,ICOUNT,OUT                           
C C                                                                     
C       PARAMETER (SQRT3=1.73205080756888D+00)                          
C       PARAMETER (SQRT5=2.23606797749979D+00)                          
C       PARAMETER (SQRT7=2.64575131106459D+00)                          
C       PARAMETER (SQRT9=3.0D+00)                                       
C       PARAMETER (SQRT11=3.3166247903553998D+00)                       
C       PARAMETER (ZERO=0.0D+00)                                        
C       PARAMETER (ONE=1.0D+00)                                         
C !!$omp declare target                                                 
C C                                                                     
C       NORM = NORMF .NE. 1 .OR. NORMP .NE. 1                           
C       MAX = MAXJ                                                      
C       N = 0                                                           
C       NN = 0                                                          
C       NM = -2**20                                                     
C       DO 180 I = MINI,MAXI                                            
C          GO TO (100,100,120,120,100,120,120,100,120,120,              
C      1          100,120,120,100,120,120,120,120,120,100,              
C      1          100,120,120,100,120,120,120,120,120,100,              
C      1          120,120,100,120,120,                                  
C      1          100,120,120,100,120,120,120,120,120,100,              
C      1          120,120,120,120,120,100,120,120,100,120,              
C      1          120,                                                  
C      1          100,120,120,100,120,120,120,120,120,100,              
C      1          120,120,120,120,120,100,120,120,100,120,              
C      1          120,100,120,120,120,120,120,100),I                    
C   100    NM = NN                                                      
C   120    NN = NM                                                      
C          IF (IANDJ) MAX = I                                           
C          DO 170 J = MINJ,MAX                                          
C             GO TO (140,140,160,160,140,160,160,140,160,160,           
C      1             140,160,160,140,160,160,160,160,160,140,           
C      1             140,160,160,140,160,160,160,160,160,140,           
C      1             160,160,140,160,160,                               
C      1             140,160,160,140,160,160,160,160,160,140,           
C      1             160,160,160,160,160,140,160,160,140,160,           
C      1             160,                                               
C      1             140,160,160,140,160,160,160,160,160,140,           
C      1             160,160,160,160,160,140,160,160,140,160,           
C      1             160,140,160,160,160,160,160,140),J                 
C   140       NN = NN+1                                                 
C   160       N = N+1                                                   
C             IJD(N) = NN                                               
C   170    CONTINUE                                                     
C   180 CONTINUE                                                        
C C                                                                     
C C     ----- I PRIMITIVE                                               
C C                                                                     
C       NIJ = 0                                                         
C       JBMAX = NGB                                                     
C       DO 540 IA = 1,NGA                                               
C          AI = AG(IA)                                                  
C          ARRI = AI*RRI                                                
C          AXI = AI*XI                                                  
C          AYI = AI*YI                                                  
C          AZI = AI*ZI                                                  
C          CSI = CSA(IA)                                                
C          CPI = CPA(IA)                                                
C          CDI = CDA(IA)                                                
C          CFI = CFA(IA)                                                
C          CGI = CGA(IA)                                                
C          CHI = CHA(IA)                                                
C          CII = CIA(IA)                                                
C C                                                                     
C C        ----- J PRIMITIVE                                            
C C                                                                     
C          IF (IANDJ) JBMAX = IA                                        
C          DO 520 JB = 1,JBMAX                                          
C             AJ = BG(JB)                                               
C             AA = AI+AJ                                                
C             AAINV = ONE/AA                                            
C             DUM = AJ*ARRI*AAINV                                       
C             IF (DUM .GT. TOL) GO TO 520                               
C             CSJ = CSB(JB)                                             
C             CPJ = CPB(JB)                                             
C             CDJ = CDB(JB)                                             
C             CFJ = CFB(JB)                                             
C             CGJ = CGB(JB)                                             
C             CHJ = CHB(JB)                                             
C             CIJ = CIB(JB)                                             
C             NM = 49*NIJ                                               
C             NN = NM                                                   
C             NIJ = NIJ+1                                               
C             R(NIJ) = DUM                                              
C             A(NIJ) = AA                                               
C             X1(NIJ) = (AXI+AJ*XJ)*AAINV                               
C             Y1(NIJ) = (AYI+AJ*YJ)*AAINV                               
C             Z1(NIJ) = (AZI+AJ*ZJ)*AAINV                               
C C                                                                     
C C           ----- DENSITY FACTOR                                      
C C                                                                     
C             DUM1 = ZERO                                               
C             DUM2 = ZERO                                               
C             DO 420 I = MINI,MAXI                                      
C                GO TO (200,220,420,420,240,420,420,260,420,420,        
C      1                261,420,420,262,420,420,420,420,420,263,        
C      1                264,420,420,265,420,420,420,420,420,266,        
C      1                420,420,267,420,420,                            
C      1                268,420,420,269,420,420,420,420,420,270,        
C      1                420,420,420,420,420,271,420,420,272,420,        
C      1                420,                                            
C      1                273,420,420,274,420,420,420,420,420,275,        
C      1                420,420,420,420,420,276,420,420,277,420,        
C      1                420,278,420,420,420,420,420,279),I              
C   200          DUM1 = CSI*AAINV                                       
C                GO TO 280                                              
C   220          DUM1 = CPI*AAINV                                       
C                GO TO 280                                              
C   240          DUM1 = CDI*AAINV                                       
C                GO TO 280                                              
C   260          DUM1 = DUM1*SQRT3                                      
C                GO TO 280                                              
C   261          DUM1 = CFI*AAINV                                       
C                GO TO 280                                              
C   262          DUM1 = DUM1*SQRT5                                      
C                GO TO 280                                              
C   263          DUM1 = DUM1*SQRT3                                      
C                GO TO 280                                              
C   264          DUM1 = CGI*AAINV                                       
C                GO TO 280                                              
C   265          DUM1 = DUM1*SQRT7                                      
C                GO TO 280                                              
C   266          DUM1 = DUM1*SQRT5/SQRT3                                
C                GO TO 280                                              
C   267          DUM1 = DUM1*SQRT3                                      
C                GO TO 280                                              
C   268          DUM1 = CHI*AAINV                                       
C                GO TO 280                                              
C   269          DUM1 = DUM1*SQRT9                                      
C                GO TO 280                                              
C   270          DUM1 = DUM1*SQRT7/SQRT3                                
C                GO TO 280                                              
C   271          DUM1 = DUM1*SQRT3                                      
C                GO TO 280                                              
C   272          DUM1 = DUM1*SQRT5/SQRT3                                
C                GO TO 280                                              
C   273          DUM1 = CII*AAINV                                       
C                GO TO 280                                              
C   274          DUM1 = DUM1*SQRT11                                     
C                GO TO 280                                              
C   275          DUM1 = DUM1*SQRT3                                      
C                GO TO 280                                              
C   276          DUM1 = DUM1*SQRT3                                      
C                GO TO 280                                              
C   277          DUM1 = DUM1*SQRT7/(SQRT5*SQRT3)                        
C                GO TO 280                                              
C   278          DUM1 = DUM1*SQRT5                                      
C                GO TO 280                                              
C   279          DUM1 = DUM1*SQRT5/SQRT3                                
C C                                                                     
C   280          IF (IANDJ) MAX = I                                     
C                DO 400 J = MINJ,MAX                                    
C                   GO TO (300,320,400,400,340,400,400,360,400,400,     
C      1                   361,400,400,362,400,400,400,400,400,363,     
C      1                   364,400,400,365,400,400,400,400,400,366,     
C      1                   400,400,367,400,400,                         
C      1                   368,400,400,369,400,400,400,400,400,370,     
C      1                   400,400,400,400,400,371,400,400,372,400,     
C      1                   400,                                         
C      1                   373,400,400,374,400,400,400,400,400,375,     
C      1                   400,400,400,400,400,376,400,400,377,400,     
C      1                   400,378,400,400,400,400,400,379),J           
C   300             DUM2 = DUM1*CSJ                                     
C                   GO TO 380                                           
C   320             DUM2 = DUM1*CPJ                                     
C                   GO TO 380                                           
C   340             DUM2 = DUM1*CDJ                                     
C                   GO TO 380                                           
C   360             DUM2 = DUM2*SQRT3                                   
C                   GO TO 380                                           
C   361             DUM2 = DUM1*CFJ                                     
C                   GO TO 380                                           
C   362             DUM2 = DUM2*SQRT5                                   
C                   GO TO 380                                           
C   363             DUM2 = DUM2*SQRT3                                   
C                   GO TO 380                                           
C   364             DUM2 = DUM1*CGJ                                     
C                   GO TO 380                                           
C   365             DUM2 = DUM2*SQRT7                                   
C                   GO TO 380                                           
C   366             DUM2 = DUM2*SQRT5/SQRT3                             
C                   GO TO 380                                           
C   367             DUM2 = DUM2*SQRT3                                   
C                   GO TO 380                                           
C   368             DUM2 = DUM1*CHJ                                     
C                   GO TO 380                                           
C   369             DUM2 = DUM2*SQRT9                                   
C                   GO TO 380                                           
C   370             DUM2 = DUM2*SQRT7/SQRT3                             
C                   GO TO 380                                           
C   371             DUM2 = DUM2*SQRT3                                   
C                   GO TO 380                                           
C   372             DUM2 = DUM2*SQRT5/SQRT3                             
C                   GO TO 380                                           
C   373             DUM2 = DUM1*CIJ                                     
C                   GO TO 380                                           
C   374             DUM2 = DUM2*SQRT11                                  
C                   GO TO 380                                           
C   375             DUM2 = DUM2*SQRT3                                   
C                   GO TO 380                                           
C   376             DUM2 = DUM2*SQRT3                                   
C                   GO TO 380                                           
C   377             DUM2 = DUM2*SQRT7/(SQRT5*SQRT3)                     
C                   GO TO 380                                           
C   378             DUM2 = DUM2*SQRT5                                   
C                   GO TO 380                                           
C   379             DUM2 = DUM2*SQRT5/SQRT3                             
C C                                                                     
C   380             NN = NN+1                                           
C                   DDIJ(NN) = DUM2                                     
C C                  write(*,*) "DDIJ is", DDIJ(NN)                     
C   400          CONTINUE                                               
C   420       CONTINUE                                                  
C             IF ( .NOT. IANDJ) GO TO 520                               
C             IF (IA .EQ. JB) GO TO 520                                 
C             GO TO (500,440,460,455,450,445,444),LIT                   
C   440       IF (MINI .EQ. 2) GO TO 500                                
C             DDIJ(NM+2) = DDIJ(NM+2)+CSI*CPJ*AAINV                     
C             GO TO 480                                                 
C   444       DDIJ(NM+28) = DDIJ(NM+28)+DDIJ(NM+28)                     
C             DDIJ(NM+27) = DDIJ(NM+27)+DDIJ(NM+27)                     
C             DDIJ(NM+26) = DDIJ(NM+26)+DDIJ(NM+26)                     
C             DDIJ(NM+25) = DDIJ(NM+25)+DDIJ(NM+25)                     
C             DDIJ(NM+24) = DDIJ(NM+24)+DDIJ(NM+24)                     
C             DDIJ(NM+23) = DDIJ(NM+23)+DDIJ(NM+23)                     
C             DDIJ(NM+22) = DDIJ(NM+22)+DDIJ(NM+22)                     
C             DDIJ(NM+21) = DDIJ(NM+21)+DDIJ(NM+21)                     
C             DDIJ(NM+20) = DDIJ(NM+20)+DDIJ(NM+20)                     
C             DDIJ(NM+19) = DDIJ(NM+19)+DDIJ(NM+19)                     
C             DDIJ(NM+18) = DDIJ(NM+18)+DDIJ(NM+18)                     
C             DDIJ(NM+17) = DDIJ(NM+17)+DDIJ(NM+17)                     
C             DDIJ(NM+16) = DDIJ(NM+16)+DDIJ(NM+16)                     
C   445       DDIJ(NM+15) = DDIJ(NM+15)+DDIJ(NM+15)                     
C             DDIJ(NM+14) = DDIJ(NM+14)+DDIJ(NM+14)                     
C             DDIJ(NM+13) = DDIJ(NM+13)+DDIJ(NM+13)                     
C             DDIJ(NM+12) = DDIJ(NM+12)+DDIJ(NM+12)                     
C             DDIJ(NM+11) = DDIJ(NM+11)+DDIJ(NM+11)                     
C   450       DDIJ(NM+10) = DDIJ(NM+10)+DDIJ(NM+10)                     
C             DDIJ(NM+9) = DDIJ(NM+9)+DDIJ(NM+9)                        
C             DDIJ(NM+8) = DDIJ(NM+8)+DDIJ(NM+8)                        
C             DDIJ(NM+7) = DDIJ(NM+7)+DDIJ(NM+7)                        
C   455       DDIJ(NM+6) = DDIJ(NM+6)+DDIJ(NM+6)                        
C             DDIJ(NM+5) = DDIJ(NM+5)+DDIJ(NM+5)                        
C             DDIJ(NM+4) = DDIJ(NM+4)+DDIJ(NM+4)                        
C   460       DDIJ(NM+2) = DDIJ(NM+2)+DDIJ(NM+2)                        
C   480       DDIJ(NM+3) = DDIJ(NM+3)+DDIJ(NM+3)                        
C   500       DDIJ(NM+1) = DDIJ(NM+1)+DDIJ(NM+1)                        
C   520    CONTINUE                                                     
C   540 CONTINUE                                                        
C       RETURN                                                          
C       END                                                             
C*MODULE INT2A   *DECK S0000_gpu                                        
      SUBROUTINE S0000_gpu(GHONDO,DDIJ,IANDJ,KANDL,SAME,                
     * LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,                             
     * qq4,mini,maxi,minj,maxj,mink,maxk,minl,maxl,                     
     * NIJ,IJ,KL,IJKL)                                                  
      USE lrcdft, ONLY: EMU2                                            
      use mx_limits, only: mxgsh,mxg2                                   
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      LOGICAL IANDJ,KANDL,SAME,OUT,LRINT                                
C                                                                       
C                                                                       
      DIMENSION GHONDO(*),DDIJ(49*MXG2)                                 
C                                                                       
      COMMON /IJGNRL/ A(MXG2),R(MXG2),X1(MXG2),Y1(MXG2),Z1(MXG2),       
     *                IJD(784)                                          
C$omp threadprivate(/IJGNRL/)
C      COMMON /MISC  / IANDJ,KANDL,SAME                                 
      COMMON /NLRCF / LRINT                                             
C$omp threadprivate(/NLRCF /)
      COMMON /SHLINF/  AG(MXGSH),CSA(MXGSH),CPA(MXGSH),CDA(MXGSH),      
     *                CFA(MXGSH),CGA(MXGSH),CHA(MXGSH),CIA(MXGSH),      
     *                 BG(MXGSH),CSB(MXGSH),CPB(MXGSH),CDB(MXGSH),      
     *                CFB(MXGSH),CGB(MXGSH),CHB(MXGSH),CIB(MXGSH),      
     *                 CG(MXGSH),CSC(MXGSH),CPC(MXGSH),CDC(MXGSH),      
     *                CFC(MXGSH),CGC(MXGSH),CHC(MXGSH),CIC(MXGSH),      
     *                 DG(MXGSH),CSD(MXGSH),CPD(MXGSH),CDD(MXGSH),      
     *                CFD(MXGSH),CGD(MXGSH),CHD(MXGSH),CID(MXGSH),      
     *                XI,YI,ZI,XJ,YJ,ZJ,RRI,XK,YK,ZK,XL,YL,ZL,RRK,      
     *                NGA,NGB,NGC,NGD                                   
C$omp threadprivate(/SHLINF/)
C      COMMON /SHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,         
C     *                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,         
C     *                NIJ,IJ,KL,IJKL                                   
      double precision :: qq4                                           
      INTEGER :: LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL                    
      INTEGER :: MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL                
      INTEGER :: NIJ,IJ,KL,IJKL                                         
      COMMON /SHLT  / TOL,CUTOFF,ICOUNT,OUT                             
C                                                                       
      PARAMETER (PI252=34.986836655250D+00)                             
      PARAMETER (PIE4=7.85398163397448D-01)                             
      PARAMETER (ZERO=0.0D+00)                                          
      PARAMETER (ONE=1.0D+00)                                           
C                                                                       
C     SPECIAL SSSS INTEGRAL ROUTINE WHEN USING HONDO INTEGRALS          
C                                                                       
      write(*,*) "I am in S0000"                                        
      GGOUT = ZERO                                                      
      LGMAX = NGD                                                       
      DO 300 KG = 1,NGC                                                 
      BK = CG(KG)                                                       
      BRRK = BK*RRK                                                     
      BXK = BK*XK                                                       
      BYK = BK*YK                                                       
      BZK = BK*ZK                                                       
      CSK = CSC(KG)                                                     
      IF (KANDL) LGMAX = KG                                             
      DO 280 LG = 1,LGMAX                                               
      BL = DG(LG)                                                       
      BB = BK+BL                                                        
      BBINV = ONE/BB                                                    
      DUM = BL*BRRK*BBINV                                               
      IF (DUM .GT. TOL) GO TO 280                                       
      BBRRK = DUM                                                       
      D2 = CSD(LG)*CSK*BBINV                                            
      IF (KANDL .AND. LG .NE. KG) D2 = D2+D2                            
      BBX = (BXK+BL*XL)*BBINV                                           
      BBY = (BYK+BL*YL)*BBINV                                           
      BBZ = (BZK+BL*ZL)*BBINV                                           
      SUM = ZERO                                                        
      NN = 1                                                            
      DO 260 N = 1,NIJ                                                  
      DUM = BBRRK+R(N)                                                  
      IF (DUM .GT. TOL) GO TO 260                                       
      EXPE = EXP(-DUM)                                                  
      AA = A(N)                                                         
      AB = AA+BB                                                        
      DUM = X1(N)-BBX                                                   
      XX = DUM*DUM                                                      
      DUM = Y1(N)-BBY                                                   
      XX = DUM*DUM+XX                                                   
      DUM = Z1(N)-BBZ                                                   
      XX = DUM*DUM+XX                                                   
      X = XX*AA*BB/AB                                                   
      IF(LRINT) THEN                                                    
         RHO = AA*BB*EMU2/(AB*EMU2+AA*BB)                               
         X   = XX*RHO                                                   
      ENDIF                                                             
C                                                                       
      IF (X .GT. 5.0D+00) GO TO 160                                     
      IF (X .GT. 1.0D+00) GO TO 120                                     
      IF (X .GT. 3.0D-07) GO TO 100                                     
      WW1 = 1.0D+00-X/3.0D+00                                           
      GO TO 240                                                         
C                                                                       
  100 CONTINUE                                                          
      F1 = ((((((((-8.36313918003957D-08*X+1.21222603512827D-06 )*X-    
     +     1.15662609053481D-05 )*X+9.25197374512647D-05 )*X-           
     +     6.40994113129432D-04 )*X+3.78787044215009D-03 )*X-           
     +     1.85185172458485D-02 )*X+7.14285713298222D-02 )*X-           
     +     1.99999999997023D-01 )*X+3.33333333333318D-01                
      WW1 = (X+X)*F1+EXP(-X)                                            
      GO TO 240                                                         
C                                                                       
  120 CONTINUE                                                          
      IF (X .GT. 3.0D+00) GO TO 140                                     
      Y = X-2.0D+00                                                     
      F1 = ((((((((((-1.61702782425558D-10*Y+1.96215250865776D-09 )*Y-  
     +     2.14234468198419D-08 )*Y+2.17216556336318D-07 )*Y-           
     +     1.98850171329371D-06 )*Y+1.62429321438911D-05 )*Y-           
     +     1.16740298039895D-04 )*Y+7.24888732052332D-04 )*Y-           
     +     3.79490003707156D-03 )*Y+1.61723488664661D-02 )*Y-           
     +     5.29428148329736D-02 )*Y+1.15702180856167D-01                
      WW1 = (X+X)*F1+EXP(-X)                                            
      GO TO 240                                                         
C                                                                       
  140 CONTINUE                                                          
      Y = X-4.0D+00                                                     
      F1 = ((((((((((-2.62453564772299D-11*Y+3.24031041623823D-10 )*Y-  
     +     3.614965656163D-09)*Y+3.760256799971D-08)*Y-                 
     +     3.553558319675D-07)*Y+3.022556449731D-06)*Y-                 
     +     2.290098979647D-05)*Y+1.526537461148D-04)*Y-                 
     +     8.81947375894379D-04 )*Y+4.33207949514611D-03 )*Y-           
     +     1.75257821619926D-02 )*Y+5.28406320615584D-02                
      WW1 = (X+X)*F1+EXP(-X)                                            
      GO TO 240                                                         
C                                                                       
  160 CONTINUE                                                          
      IF (X .GT. 15.0D+00) GO TO 200                                    
      E = EXP(-X)                                                       
      IF (X .GT. 10.0D+00) GO TO 180                                    
      XINV = ONE/X                                                      
      WW1 = (((((( 4.6897511375022D-01*XINV-6.9955602298985D-01)*XINV + 
     +     5.3689283271887D-01)*XINV-3.2883030418398D-01)*XINV +        
     +     2.4645596956002D-01)*XINV-4.9984072848436D-01)*XINV -        
     +     3.1501078774085D-06)*E + SQRT(PIE4*XINV)                     
      GO TO 240                                                         
C                                                                       
  180 CONTINUE                                                          
      XINV = ONE/X                                                      
      WW1 = (((-1.8784686463512D-01*XINV+2.2991849164985D-01)*XINV      
     +         -4.9893752514047D-01)*XINV-2.1916512131607D-05)*E        
     +         + SQRT(PIE4*XINV)                                        
      GO TO 240                                                         
C                                                                       
  200 CONTINUE                                                          
      IF (X .GT. 33.0D+00) GO TO 220                                    
      XINV = ONE/X                                                      
      E = EXP(-X)                                                       
      WW1 = (( 1.9623264149430D-01*XINV-4.9695241464490D-01)*XINV -     
     +     6.0156581186481D-05)*E + SQRT(PIE4*XINV)                     
      GO TO 240                                                         
C                                                                       
  220 WW1 = SQRT(PIE4/X)                                                
C                                                                       
  240 CONTINUE                                                          
      IF(.NOT.LRINT)SUM = SUM+DDIJ(NN)*WW1*EXPE/SQRT(AB)                
      IF(     LRINT)SUM = SUM+DDIJ(NN)*WW1*EXPE/SQRT(AA*BB/RHO)         
  260 NN = NN+49                                                        
      GGOUT = GGOUT+D2*SUM                                              
  280 CONTINUE                                                          
  300 CONTINUE                                                          
      GHONDO(1) = GGOUT*PI252*QQ4                                       
      RETURN                                                            
      END                                                               
