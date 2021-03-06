c 22 Oct 14 - DGF - changes for FMO 5.1                                 
C 19 Mar 14 - MWS - TSECND: spell wall timer as XQTIME for Linux/Mac    
C 27 Apr 12 - MWS - provide MEMGET/MEMREL for HP-UX                     
C 28 DEC 11 - DGF - REDUCE I/O BY NOT FLUSHING                          
C 15 APR 11 - MWS - ADD SIMULATION OF RANDOM_NUMBER INTRINSIC           
C 13 MAY 10 - SS  - ADD *W32/*W64 WINDOWS SUPPORT                       
C 10 MAY 10 - AA  - 3 ROUTINES TO CUSTOMIZE SERIAL EXEC. MATRIX SIZES   
C 18 JUL 08 - MWS - VNAN: NEW ROUTINE TO CREATE BAD VALUES              
C 14 FEB 07 - AA  - ADD "FGE" ROUTINE TO GET ENV. VALUES FROM A DISKFILE
C 22 DEC 06 - MWS,SK - LINUX PORTS *F2C,*AMD,*ITA RENAMED *L32,*L64,*INT
C 22 DEC 06 - DGF - ADD TIMING PRINT-OUT ON SLAVES                      
C 10 JUL 06 - SK,MWS - TSECND: CHANGE *AMD TIMING TO USUAL LINUX WAY    
C 17 JAN 06 - RAK - ADD FLUSHONWRITE ROUTINE FOR IBM/AIX                
C 19 NOV 05 - TJP - ADD *XT3 CODE FOR CRAY XT3                          
C  5 JUL 05 - MWS - SELECT NEW ATOM,BASIS,EFP,PCM,DAF DIMENSIONS        
C  1 JUN 05 - MWS - FLSHBF: IBM/AIX SHOULD PASS 4 BYTE ARG TO FLUSH_    
C  5 FEB 05 - MWS - FLSHBF: ADD ITANIUM CODING                          
C  7 SEP 04 - MWS - PAD COMMON BLOCK INTFIL                             
C 23 JUL 03 - SK,HU,TJP - INCLUDE 64 BIT AMD CHIP VERSION (*AMD)        
C  4 NOV 03 - TJP - INCLUDE ROUTINE FOR EFFECTIVE VECTOR LENGTH         
C  3 SEP 03 - MWS - TSECND: *ITA TIMING CALLS FIXED                     
C  3 JUL 03 - JMM - SUPPRESS PRINTING FOR MONTE CARLO JOBS              
C 26 MAR 03 - TJP - ADD *CX1 CODE FOR CRAY-X1                           
C 28 JAN 03 - MWS - ADAPTATIONS FOR A 64 BIT SUN VERSION                
C 12 DEC 02 - YA  - GETFM: LET ANY NODE PRINT MEMORY ALLOCATION ERROR   
C  7 AUG 02 - ZK  - USE FORTRAN LIBRARY ROUTINES FOR THE HP VERSION     
C 20 JUN 02 - FPA - ADD *ITA ITANIUM LINUX VERSION                      
C 17 APR 02 - KRG - REMOVE A STOP STATEMENT                             
C 16 FEB 02 - MWS - TSECND,TMDATE: CHANGE IBM TIMING CALLS              
C 24 JAN 02 - MWS - ABRT: STALL 2 SECONDS TO ALLOW FOR OUTPUT FLUSH     
C  7 DEC 01 - STE - ADD *WIN LINES FOR WIN32 VERSION                    
C  8 OCT 01 - MWS - PARSET: OMIT BALTYP OVERRIDES, UPDATE ITRFAO VALS   
C  6 SEP 01 - MWS - FLSHBF: DIGITAL CHANGED                             
C 10 JUL 01 - MWS - TMDATE: ADD 1900 TO HP'S YEAR VALUE                 
C  1 MAY 00 - MWS - CHANGES TO MATCH MODIFIED ZUNIX.C                   
C 25 MAR 00 - MWS - REMOVE BRIAN'S GO AT A SEPARATE ALPHA LINUX VERSION 
C 25 SEP 99 - MWS - BEGING: REMOVE ZMIPS CALL                           
C 29 AUG 99 - BSD - ADD STATIC MEMORY ALPHA LINUX CODE                  
C  6 JUN 99 - MWS - ADD CALLING ARGS TO DDI_PBEG AND DDI_PEND           
C 12 NOV 98 - MWS - T3D TO T3E, DEL CONVEX (CVX), LOGSHF, DEC=64 BIT    
C 22 APR 98 - MWS - ENDING: AVOID FINAL FLOATING PT. WARNINGS ON SUN    
C 13 APR 98 - AJL - FLUSH OUTPUT ON FUJITSU                             
C 27 FEB 98 - AR,AR - CORRECT FUJITSU LOGICAL FUNCTIONS, CHANGE TIMING  
C  6 JAN 98 - MY  - INCLUDE HITACHI SR VERSION                          
C 28 SEP 97 - MWS - DEPRECATE MANY OLD MACHINES. NAMELY, WE DROP SUPPORT
C                   FOR *ALL, *AMD, *APO, *CEL, *INT, *KSR, *STR, *TMC  
C  8 SEP 97 - MWS - CONVERT PARALLEL CALLS FROM TCGMSG TO DDI           
C  2 SEP 97 - MWS - ABRT: INSERT PROPER F2C CODE                        
C 24 JUL 97 - KRG - TIMIT: CHANGE FORMATTING TO 2 LINES                 
C 22 JUL 97 - DGF - ADD BIT OPERATIONS LOGNOT, LOGOR, LOGSHF, LOGXOR    
C 19 FEB 97 - FPR - TMDATE: DECLARATION FOR HP                          
C 20 DEC 96 - HPP - BEGING,TSECND: T3E IN BANNER, ADD CPU TIMING        
C 12 MAY 96 - MWS - PARSET,MSGSET: CHANGES FOR MPI INTERFACE            
C 24 MAY 95 - MWS - PARSET: ASSIGN PARALLEL AO TRANSFORMATION STORAGE   
C 29 MAR 95 - MK  - AIX'S XLF VERSION 3 CAN FLUSH_                      
C  6 MAR 95 - KPG - TSECND: USE AN ACTUAL TIMING CALL                   
C 17 NOV 94 - MWS - ADD TEXIT,TIMIT TO CONSOLIDATE ALL TIMING CODE      
C 11 NOV 94 - MWS - REMOVE A FEW FTNCHEK WARNINGS                       
C 28 OCT 94 - MWS - BEGING: USE XDR UNTIL AFTER $SYSTEM IS READ         
C 23 AUG 94 - ML  - INCLUDE CRAY-T3D VERSION                            
C 10 AUG 94 - MWS - INCREASE NUMBER OF DAF RECORDS                      
C 19 JUL 94 - MWS - ADD ROUTINE ASKVEC                                  
C 19 JUN 94 - DM  - TSECND: CHANGE HP TIMING CALL                       
C  4 APR 94 - MWS - MSGSET: XDR USAGE IS OPTIONAL                       
C 18 MAR 94 - MM  - TSECND: CONVEX CPUTIME IS SINGLE PRECISION          
C 29 DEC 93 - PVZ - ADDED F2C OPTION FOR F2C/GCC COMPILATION            
C 31 SEP 93 - GJA - ADD FLUSH FOR AIX                                   
C 24 SEP 93 - JF  - ADD *NEC LINES FOR NEC                              
C 20 SEP 93 - WS,JAB - ADD *KSR LINES                                   
C  7 SEP 93 - SS  - BEGING: CM-5 NEEDS MODE OF STDERR RESET             
C 19 AUG 93 - RHN - TMDATE: CHANGE HP TIME STAMP CALLS                  
C 24 MAY 93 - NAN - MSGSET: DON'T USE MESSAGE MASKS FOR CRAY            
C 10 FEB 93 - TLW - THINKING MACHINES VERSION                           
C 31 OCT 92 - MWS - TMDATE: FIX COUNTING LEAP YEARS FOR AIX             
C 20 JUN 92 - DP,FS,CJ,TH - HP9000-7X0 TIMING CODE                      
C 17 APR 92 - TLW,MWS - ADD PARSET, OTHER IPSC/860 MODIFICATIONS        
C  7 APR 92 - MS  - SETFM: USE MALLOPT(%VAL(1),%VAL(0)) FOR AIX 3.1.5   
C 26 MAR 92 - TLW - MESSAGE MASKS DEFINED IN ROUTINE MSGSET             
C  2 MAR 92 - TLW - FLSHBF: ONLY MASTER WORK                            
C 30 JAN 92 - TLW - TMDATE: FIX PARALLEL CALLS                          
C 29 JAN 92 - TLW - BEGING: PUT IN PARALLEL INITIALIZATION;             
C                   ENDING: PUT IN PARALLEL CLEANUP                     
C 11 JAN 92 - TLW - ABRT: RUN IN PARALLEL                               
C  8 JAN 92 - TLW - MAKE WRITES PARALLEL; ADD COMMON PAR                
C  1 JAN 92 - MWS - FLSHBF: UNICOS VERSION 6 HAS FLUSH CALL             
C 18 DEC 91 - MWS - PUT IN FRED SENESE'S HP CODE                        
C 18 NOV 91 - NM  - ABRT: NO ARGUMENT IN AIX CALL TO ABORT              
C  6 SEP 91 - JAB - ADD INT FOR INTEL IPSC/860 RUNNING NX/2             
C 28 AUG 91 - RHN - ADD FUJITSU UXP/M SUPPORT                           
C 20 JUN 91 - MWS - CHANGE IN TMDATE'S DECSTATION CODE FOR F77 3.0      
C 15 FEB 91 - TLW - PRAGMA ADDED TO ABRT FOR AIX/370                    
C  4 FEB 91 - TLW - ABRT: ADD IF ICORFL LINES                           
C 29 JAN 91 - TLW - ABRT: ADD COMMON MACHSW AND FLUSHES                 
C  8 DEC 90 - UK  - TMDATE: FIX LEAP YEAR COMPUTATION FOR AIX           
C 15 SEP 90 - MWS - MOVE UPCASE BACK INTO FRFMT SECTION                 
C  1 AUG 90 - MWS - AIX CHANGED TO USE DYNAMIC MEMORY                   
C 20 JUN 90 - FRJ - TSECND: CONVEX CPUTIME IS DOUBLE PRECISION          
C  1 MAY 90 - STE - BEGING: ENABLE FP EXCEPT. ON MIPS BASED SYSTEMS     
C  5 APR 90 - JHJ - ADD UPCASE ROUTINE                                  
C 20 MAR 90 - MG  - ADD ETIME SIMULATOR FOR 68000 BASED APOLLOS.        
C 23 FEB 90 - MWS - DECSTN MEMORY ALLOCATION MUST HAVE EXTERNAL LOC,    
C                   FLUSH PUNCH BUFFER IN BIGFM.                        
C  6 FEB 90 - MWS - SETFM: LIMFM SHOULD BE MEMLIM, NOT MEMAVL           
C 22 JAN 90 - MWS - UNICOS VERSION, WITH DYNAMIC MEMORY                 
C 16 JAN 90 - MWS - SETFM: EXPLANATORY COMMENTS, DYNAMIC MEMORY ALA     
C                   STEVE ELBERT FOR ALL,ARD,CEL,CVX,SUN VERSIONS.      
C  3 JAN 90 - STE - DYNAMIC MEMORY SUPPORT FOR DEC, SGI VERSIONS        
C 29 DEC 89 - JAM - CONVEX CAN FLUSH BUFFERS TOO.                       
C 22 SEP 89 - MWS - ADD STEVE'S APO, SGI, STR, SUN VERSIONS             
C 21 SEP 89 - STE - MOVE LOGAND FROM MTHLIB TO HERE                     
C 21 AUG 89 - CFJ - ABRT CALLS TRACEBACK FOR CONVEX VERSION             
C 18 AUG 89 - MWS - ADD IBM AIX SUPPORT                                 
C 11 JUL 89 - STE - ADD DECSTATION SUPPORT (MIPS CHIP)                  
C  5 JUL 89 - FJ  - ADD AMDAHL SUPPORT                                  
C  7 JUN 89 - MWS - ADD ARDENT TITAN SUPPORT                            
C 10 APR 89 - STE - ABORT IF REQUESTED MEMORY EXCEEDS AMOUNT AVAILABLE  
C 22 FEB 89 - STE - PUT ALL THE NASTY MACHINE DEPENDENCIES HERE         
C*MODULE UNPORT                                                         
C-----------------------------------------------------------------------
C     THIS IS THE NASTIEST, TOTALLY UNPORTABLE CODE IN GAMESS.          
C     THE VERSIONS WHICH ARE SUPPORTED HERE ARE                         
C          *AIX - IBM POWER, POWERPC, SP, ETC RUNNING AIX (SEE IRON.DOC)
C          *CRY - CRAY PARALLEL VECTOR PROCESSOR (OLDER MODELS),        
C                 SEE ALSO OTHER CRAY PRODUCTS *T3E, *CX1, *XT3         
C          *CX1 - CRAY X1 OR X1E (THE MODERN VECTOR SYSTEMS) (UNICOS/MP)
C          *DEC - DECSTATION/DECSYSTEM (UNIX)                           
C          *F77 - A GENERIC, AND MOSTLY DO NOTHING VERSION              
C          *FUJ - FUJITSU UXP/M                                         
C          *HIT - HITACHI SR2201(RUNNING HI-UX/MPP)                     
C          *HP  - HEWLETT-PACKARD HP/9000 SERIES                        
C          *IBM - IBM MAINFRAME (VM OR MVS, VS FORTRAN) - SEE ALSO *AIX 
C          *INT - 64 BIT LINUX ON INTEL ITAN2/EMT64, WITH IFORT COMPILER
C          *L32 - 32 BIT LINUX, USING G77, GFORTRAN, PGF77, IFORT ...   
C          *L64 - 64 BIT LINUX, USING      GFORTRAN, PGF77, IFORT ...   
C                 THESE TWO MUST CONFORM TO GNU/GCCLIB CALLS, ONLY.     
C          *NEC - NEC SX (SUPER-UX)                                     
C          *SGI - SILICON GRAPHICS INC.                                 
C          *SUN - SUN WORKSTATIONS                                      
C          *T3E - CRAY T3E (T3D MAY NOT WORK)                           
C          *VMS - VMS SYSTEMS RUNNING ON EITHER VAX OR AXP PROCESSORS   
C          *XT3 - CRAY XT3 (CATAMOUNT)                                  
C          *W64 - WIN64 SYSTEMS USING THE PGI COMPILERS   : PGFORTRAN   
C          *W32 - WIN32 SYSTEMS USING THE PGI COMPILERS   : PGFORTRAN   
C          *I64 - WIN64 SYSTEMS USING THE INTEL COMPILERS : IFORT       
C          *WAZ - AZURE TIMER                             : w/ SERIAL   
C          *64W - WIN64 TIMER                             : PGI w/ MPI  
C          *64I - WIN64 TIMER                             : INTEL w/ MPI
C-----------------------------------------------------------------------
C                                                                       
C*MODULE UNPORT  *DECK ABRT                                             
      SUBROUTINE ABRT                                                   
C                                                                       
C     ----- GENERATE A CALLING SEQUENCE TRACE AND STOP EXECUTION -----  
C                                                                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      DIMENSION TIMSTR(3)                                               
C                                                                       
      LOGICAL GOPARR,DSKWRK,MASWRK                                      
C                                                                       
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)              
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
      COMMON /MACHIN/ NWDVAR,MAXFM,MAXSM,LIMFM,LIMSM                    
      COMMON /MACHSW/ KDIAG,ICORFL,IXDR,modio,mem10,lpnt10,mem10m       
C                                                                       
C     ----- PRINT ACCOUNTING INFO AND GOODBYE MESSAGE -----             
C                                                                       
      CALL TMDATE(TIMSTR)                                               
      IF (MASWRK) WRITE(IW,900) TIMSTR                                  
  900 FORMAT(1X,'EXECUTION OF GAMESS TERMINATED -ABNORMALLY- AT ',3A8)  
      CALL BIGFM(MAXFM)                                                 
      CALL TIMIT(1)                                                     
C                                                                       
C     STALL FOR A SHORT WHILE TO ALLOW A CHANCE OF PRINTING             
C     ANY ERROR MESSAGES, BEFORE DDIKICK.X BLOWS US AWAY.               
C     MAINLY THIS IS IN CASE A NON-MASTER PROCESS ENTERS FIRST,         
C     MOMENTS BEFORE THE MASTER, WHICH SHOULD BE GIVEN A CHANCE         
C     TO PRINT THE ERROR MESSAGE BEFORE THE NON-MASTER EXITS.           
C                                                                       
      CALL FLSHBF(IW)                                                   
      CALL FLSHBF(IP)                                                   
      IDELAY=1                                                          
      IF(ME.GT.0) IDELAY=10                                             
      CALL NAPTIME(IDELAY)                                              
C                                                                       
C     ----- EXIT PARALLEL RUNS GRACELESSLY -----                        
C                                                                       
      CALL DDI_PEND(1)                                                  
C                                                                       
C     ----- DO A TRACEBACK AND/OR CORE DUMP -----                       
C                                                                       
      IF (MASWRK) THEN                                                  
C-----------------------------------------------------------------------
C                          IBM RUNNING AIX                              
      IF (ICORFL.EQ.1) THEN                                             
         WRITE(IW,1) 'DBX','WHERE','QUIT'                               
         CALL ABORT                                                     
      ELSE                                                              
         WRITE(IW,2)                                                    
      END IF                                                            
C-----------------------------------------------------------------------
C                          CRAY UNICOS                                  
*CRY  IF (ICORFL.EQ.1) THEN                                             
*CRY     CALL TRBK(6)                                                   
*CRY     WRITE(IW,1) 'CDBX','WHERE','QUIT'                              
*CRY     CALL ABORT('CORE DEBUG FILE CREATED')                          
*CRY  ELSE                                                              
*CRY     WRITE(IW,2)                                                    
*CRY     CALL TRBK(6)                                                   
*CRY  END IF                                                            
C-----------------------------------------------------------------------
C                          CRAY X1 - UNICOS/MP                          
C          (IDENTICAL TO *CRY, EXCEPT DEBUGGER HAS CHANGED)             
*CX1  IF (ICORFL.EQ.1) THEN                                             
*CX1     WRITE(IW,1) 'GDB','WHERE','QUIT'                               
*CX1     CALL ABORT('CORE DEBUG FILE CREATED')                          
*CX1  ELSE                                                              
*CX1     WRITE(IW,2)                                                    
*CX1  END IF                                                            
C-----------------------------------------------------------------------
C                          DIGITAL EQUIPMENT CORPORATION                
*DEC  IF (ICORFL.EQ.1) THEN                                             
*DEC     WRITE(IW,1) 'DBX','WHERE','QUIT'                               
*DEC     CALL ABORT('CORE DEBUG FILE CREATED')                          
*DEC  ELSE                                                              
*DEC     WRITE(IW,2)                                                    
*DEC  END IF                                                            
C-----------------------------------------------------------------------
C                          GENERIC FORTRAN 77                           
*F77  BYEBYE = SQRT(REAL(-IW))                                          
C-----------------------------------------------------------------------
C                          FUJITSU                                      
*FUJ  CALL ERRTRA                                                       
C-----------------------------------------------------------------------
C                          HITACHI HI-UX/MPP                            
*HIT  IF (ICORFL.EQ.1) THEN                                             
*HIT     CALL ERRTRA                                                    
*HIT     WRITE(IW,1) 'NDB','WHERE','QUIT'                               
*HIT     CALL HF_CFUNC('ABORT',0,ITRASH,0)                              
*HIT  ELSE                                                              
*HIT     WRITE(IW,2)                                                    
*HIT     CALL ERRTRA                                                    
*HIT  END IF                                                            
C-----------------------------------------------------------------------
C                          HEWLETT-PACKARD                              
*HP   IF (ICORFL.EQ.1) THEN                                             
*HP      WRITE(IW,1) 'ADB','C','Q'                                      
*HP      CALL EXIT(1)                                                   
*HP   ELSE                                                              
*HP      WRITE(IW,2)                                                    
*HP   END IF                                                            
C-----------------------------------------------------------------------
C                          INTERNATIONAL BUSINESS MACHINES              
*IBM  CALL ERRTRA                                                       
C-----------------------------------------------------------------------
C                          LINUX ON ITANIUM/EMT64 USING IFORT COMPILER  
*INT  IF (ICORFL.EQ.1) THEN                                             
*INT     WRITE(IW,1) 'GDB','BACKTRACE','QUIT'                           
*INT     CALL ABORT                                                     
*INT  ELSE                                                              
*INT     WRITE(IW,2)                                                    
*INT  END IF                                                            
C-----------------------------------------------------------------------
C                          GENERIC 32 BIT LINUX (X86, AMD, G4, G5...)   
*L32  IF (ICORFL.EQ.1) THEN                                             
*L32     WRITE(IW,1) 'GDB','BACKTRACE','QUIT'                           
*L32     CALL ABORT                                                     
*L32  ELSE                                                              
*L32     WRITE(IW,2)                                                    
*L32  END IF                                                            
C-----------------------------------------------------------------------
C                          GENERIC 64 BIT LINUX (OPTERON, EMT64, ...)   
*L64  IF (ICORFL.EQ.1) THEN                                             
*L64     WRITE(IW,1) 'GDB','BACKTRACE','QUIT'                           
*L64     CALL ABORT                                                     
*L64  ELSE                                                              
*L64     WRITE(IW,2)                                                    
*L64  END IF                                                            
C-----------------------------------------------------------------------
C                          NEC SX SUPER-UX                              
*NEC  IF (ICORFL.EQ.1) THEN                                             
*NEC     WRITE(IW,1) 'DBX','$WHERE','$Q'                                
*NEC     CALL MESPUT('CORE DEBUG FILE CREATED',23,1)                    
*NEC     CALL ABORT()                                                   
*NEC  ELSE                                                              
*NEC     WRITE(IW,2)                                                    
*NEC     CALL MESPUT('ABORT FROM GAMESS',17,1)                          
*NEC  END IF                                                            
C-----------------------------------------------------------------------
C                          SILICON GRAPHICS                             
*SGI  IF (ICORFL.EQ.1) THEN                                             
*SGI     WRITE(IW,1) 'DBX','$WHERE','$Q'                                
*SGI     CALL ABORT('CORE DEBUG FILE CREATED')                          
*SGI  ELSE                                                              
*SGI     WRITE(IW,2)                                                    
*SGI  END IF                                                            
C-----------------------------------------------------------------------
*SUN  IF (ICORFL.EQ.1) THEN                                             
*SUN     WRITE(IW,1) 'DBX','$WHERE','$Q'                                
*SUN     CALL ABORT('CORE DEBUG FILE CREATED')                          
*SUN  ELSE                                                              
*SUN     WRITE(IW,2)                                                    
*SUN  END IF                                                            
C-----------------------------------------------------------------------
C                          CRAY T3E (IDENTICAL TO *CRY)                 
*T3E  IF (ICORFL.EQ.1) THEN                                             
*T3E     CALL TRBK(6)                                                   
*T3E     WRITE(IW,1) 'CDBX','WHERE','QUIT'                              
*T3E     CALL ABORT('CORE DEBUG FILE CREATED')                          
*T3E  ELSE                                                              
*T3E     WRITE(IW,2)                                                    
*T3E     CALL TRBK(6)                                                   
*T3E  END IF                                                            
C-----------------------------------------------------------------------
C                          DIGITAL EQUIPMENT CORPORATION                
*VMS  CALL LIB$STOP(%VAL(0))                                            
C-----------------------------------------------------------------------
C                          CRAY XT3                                     
*XT3  IF (ICORFL.EQ.1) THEN                                             
*XT3     WRITE(IW,1) 'GDB','BACKTRACE','QUIT'                           
*XT3     CALL ABORT('CORE DEBUG FILE CREATED')                          
*XT3  ELSE                                                              
*XT3     WRITE(IW,2)                                                    
*XT3  END IF                                                            
C-----------------------------------------------------------------------
C                          WINDOWS 32 BIT                               
*W32  IF (ICORFL.EQ.1) THEN                                             
*W32     WRITE(IW,1) 'GDB','BACKTRACE','QUIT'                           
*W32     CALL ABORT                                                     
*W32  ELSE                                                              
*W32     WRITE(IW,2)                                                    
*W32  END IF                                                            
C-----------------------------------------------------------------------
C                          WINDOWS 64 BIT                               
*W64  IF (ICORFL.EQ.1) THEN                                             
*W64     WRITE(IW,1) 'GDB','BACKTRACE','QUIT'                           
*W64     CALL ABORT                                                     
*W64  ELSE                                                              
*W64     WRITE(IW,2)                                                    
*W64  END IF                                                            
C-----------------------------------------------------------------------
C              MS WINDOWS & 64 BIT INTEL CHIPS, USING IFORT             
C                             SPECIFIC PORTION                          
*I64  IF (ICORFL.EQ.1) THEN                                             
*I64     WRITE(IW,1) 'GDB','BACKTRACE','QUIT'                           
*I64     CALL ABORT                                                     
*I64  ELSE                                                              
*I64     WRITE(IW,2)                                                    
*I64  END IF                                                            
C-----------------------------------------------------------------------
C                                                                       
      END IF                                                            
C                                                                       
C        ----- GENERIC STOP, IN CASE YOU GET THIS FAR -----             
C                                                                       
      STOP 'IN ABRT'                                                    
C                                                                       
    1 FORMAT(1X,'TRACEBACK CAN BE OBTAINED BY LOCATING',                
     *          ' THE ''CORE'' MEMORY DUMP FILE'/                       
     *   1X,'AND TYPING THE FOLLOWING (IN LOWER CASE)'/                 
     *   1X,A,' /U.../GAMESS/GAMESS.01.X CORE (INVOKES DEBUGGER)'/      
     *   1X,'   ',A,'                         (PRINTS TRACEBACK)'/      
     *   1X,'   ',A,'                         (QUITS DEBUGGER)')        
    2 FORMAT(1X,'IF YOU WANT A CORE FILE, SET COREFL=.TRUE.',           
     *          ' IN $SYSTEM.')                                         
      END                                                               
C*MODULE UNPORT  *DECK BEGING                                           
      SUBROUTINE BEGING(VERSN)                                          
C                                                                       
C------------------------------------------------                       
C     ----- TAKE CARE OF RUN INITIALIZATION -----                       
C------------------------------------------------                       
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      LOGICAL GOPARR,DSKWRK,MASWRK                                      
C                                                                       
      CHARACTER*40 VERSN                                                
      CHARACTER*1  DIRSEP                                               
C                                                                       
C                                                                       
      COMMON /MACHIN/ NWDVAR,MAXFM,MAXSM,LIMFM,LIMSM                    
      COMMON /MACHSW/ KDIAG,ICORFL,IXDR,modio,mem10,lpnt10,mem10m       
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
      COMMON /FILESY/ DIRSEP                                            
C                                                                       
C-----------------------------------------------------------------------
C     ----- SET -NWDVAR- IN COMMON /MACHIN/ -----                       
C     THIS IS THE NUMBER OF INTEGERS IN A WORKING PRECISION             
C     FLOATING POINT NUMBER, WHICH IS MOST COMMONLY 64 BITS.            
C     THUS SYSTEMS WITH 64 BIT INTEGERS USE A VALUE OF 1, WHILE         
C          SYSTEMS WITH 32 BIT INTEGERS USE A VALUE OF 2.               
C                                                                       
      NWDVAR = 1                                                        
*CRY  NWDVAR = 1                                                        
*CX1  NWDVAR = 1                                                        
*DEC  NWDVAR = 1                                                        
*INT  NWDVAR = 1                                                        
*L64  NWDVAR = 1                                                        
*NEC  NWDVAR = 1                                                        
*T3E  NWDVAR = 1                                                        
*XT3  NWDVAR = 1                                                        
*W64  NWDVAR = 1                                                        
*I64  NWDVAR = 1                                                        
c     directory separator symbol in UNIX                                
      DIRSEP='/'                                                        
c     that other operating system                                       
*W32  DIRSEP='\'                                                        
*W64  DIRSEP='\'                                                        
C                                                                       
      MAXFM=0                                                           
      MAXSM=0                                                           
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
C INITIALIZE PARALLEL                                                   
C                                                                       
      CALL DDI_PBEG(NWDVAR)                                             
      CALL DDI_NPROC(NPROC,ME)                                          
      MASTER = 0                                                        
      GOPARR = NPROC.GT.1                                               
      MASWRK = ME.EQ.MASTER                                             
      DSKWRK = .FALSE.                                                  
      IXDR   = 1                                                        
C-----------------------------------------------------------------------
C                                                                       
C     ----- SPECIAL ERROR/CONDITION HANDLING -----                      
C                                                                       
C        KILL UNDERFLOW ERROR LOGGING ON IBM MAINFRAME                  
*IBM  CALL XUFLOW(0)                                                    
C                                                                       
C        KILL OUTPUT CONVERSION ERRORS FOR VAX/VMS SYSTEMS              
C        THIS LINE SHOULD NOT BE USED ON ALPHA VMS SYSTEMS              
CVAX  CALL ERRSET(63,,,,.FALSE.,1000)                                   
C-----------------------------------------------------------------------
C                                                                       
C        ----- DEFINE THE MACHINE VERSION -----                         
C                                                                       
      VERSN = '*********** IBM (AIX) VERSION **********'                
*CRY  VERSN = '*********** CRAY PVP VERSION ***********'                
*CX1  VERSN = '************ CRAY X1 VERSION ***********'                
*DEC  VERSN = '********* COMPAQ (DEC) VERSION *********'                
*F77  VERSN = '****** GENERIC FORTRAN 77 VERSION ******'                
*FUJ  VERSN = '******** FUJITSU UXP/M VERSION *********'                
*HIT  VERSN = '****** HITACHI HI-UX/MPP VERSION *******'                
*HP   VERSN = '*** HEWLETT-PACKARD (HP-UX) VERSION ****'                
*IBM  VERSN = '********* IBM (MVS/VM) VERSION *********'                
*INT  VERSN = '********* 64 BIT INTEL VERSION *********'                
*L32  VERSN = '********* 32 BIT LINUX VERSION *********'                
*L64  VERSN = '********* 64 BIT LINUX VERSION *********'                
*NEC  VERSN = '******* NEC SX (SUPER-UX) VERSION ******'                
*SGI  VERSN = '******** SILICON GRAPHICS VERSION ******'                
*SUN  VERSN = '***** SUN MICROSYSTEMS INC. VERSION ****'                
*T3E  VERSN = '*********** CRAY T3E VERSION ***********'                
*VMS  VERSN = '********* VAX/AXP (VMS) VERSION ********'                
*XT3  VERSN = '*********** CRAY XT3 VERSION ***********'                
*W32  VERSN = '******** 32 BIT WINDOWS VERSION ********'                
*W64  VERSN = '******** 64 BIT WINDOWS VERSION ********'                
*I64  VERSN = '***** 64 BIT WINDOWS INTEL VERSION *****'                
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE UNPORT  *DECK ENDING                                           
      SUBROUTINE ENDING                                                 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)              
C                                                                       
C     ----- TERMINATE EXECUTION SMOOTHLY -----                          
C     THE TIME DELAY IS TO TRY TO ALLOW ALL BUFFERS TO FLUSH PROPERLY   
C                                                                       
      CALL FLSHBF(IW)                                                   
      CALL FLSHBF(IP)                                                   
      IDELAY=2                                                          
      CALL NAPTIME(IDELAY)                                              
C                                                                       
C        CLEAN UP PARALLEL EXECUTION                                    
C                                                                       
      CALL DDI_PEND(0)                                                  
C                                                                       
C        THIS AVOIDS THE "FORTRAN STOP" MESSAGE                         
*VMS  CALL EXIT                                                         
*W32  CALL EXIT                                                         
*W64  CALL EXIT                                                         
      RETURN                                                            
      END                                                               
C*MODULE UNPORT  *DECK FLSHBF                                           
      SUBROUTINE FLSHBF(LUNIT)                                          
C                                                                       
      LOGICAL GOPARR,DSKWRK,MASWRK                                      
      INTEGER*4 LTEMP                                                   
*CX1  INTEGER*4 LTEMP                                                   
*DEC  INTEGER*4 LTEMP                                                   
*INT  INTEGER*4 LTEMP                                                   
*INT  LOGICAL*4 STAT,COMMITQQ                                           
*SGI  INTEGER*4 LTEMP                                                   
*SUN  INTEGER*4 LTEMP,FLUSH,ISTAT                                       
*I64  INTEGER*4 LTEMP                                                   
*I64  LOGICAL*4 STAT,COMMITQQ                                           
C                                                                       
      COMMON /MACHSW/ KDIAG,ICORFL,IXDR,modio,mem10,lpnt10,mem10m       
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
C                                                                       
C        ----- FLUSH THE BUFFER FOR LOGICAL UNIT LUNIT -----            
C        THIS ROUTINE IS MOST IMPORTANT FOR UNIX SYSTEMS,               
C        WHERE OUTPUT OTHERWISE STAYS IN BUFFERS FOREVER.               
C                                                                       
C     NOTE: ON AIX SYSTEMS WITH XLF VERSION 2, THE FLUSH CALL DOES      
C           NOT WORK WELL, AND SHOULD BE COMMENTED OUT.  THE MOST       
C           PROMINENT SIGN OF A BAD "FLUSH_" IS GARBAGE CHARACTERS      
C           IN THE PUNCH FILE.  E.G., CHANGE TO "C---*AIX ..."          
C                                                                       
      if(iand(modio,1).ne.0) return                                     
      IF (MASWRK) THEN                                                  
      LTEMP = LUNIT                                                     
      CALL FLUSH_(LTEMP)                                                
*CRY  CALL FLUSH(LUNIT)                                                 
*CX1  LTEMP = LUNIT                                                     
*CX1  CALL FLUSH(LTEMP)                                                 
*DEC  LTEMP = LUNIT                                                     
*DEC  CALL FLUSH(LTEMP)                                                 
*FUJ  CALL FLUSH(LUNIT)                                                 
*HIT  CALL HF_FLUSH(LUNIT)                                              
*HP   CALL FLUSH(LUNIT)                                                 
*INT  LTEMP = LUNIT                                                     
*INT  STAT = COMMITQQ(LTEMP)                                            
*L32  CALL FLUSH(LUNIT)                                                 
*L64  CALL FLUSH(LUNIT)                                                 
*NEC  CALL FLUSH(LUNIT)                                                 
*SGI  LTEMP = LUNIT                                                     
*SGI  CALL FLUSH(LTEMP)                                                 
*SUN  LTEMP = LUNIT                                                     
*SUN  ISTAT = FLUSH(LTEMP)                                              
*T3E  CALL FLUSH(LUNIT)                                                 
*XT3  CALL FLUSH(LUNIT)                                                 
*W32  CALL FLUSH(LUNIT)                                                 
*W64  CALL FLUSH(LUNIT)                                                 
*I64  LTEMP = LUNIT                                                     
*I64  STAT = COMMITQQ(LTEMP)                                            
      END IF                                                            
      RETURN                                                            
      END                                                               
C*MODULE UNPORT  *DECK FLUSHONWRITE                                     
      SUBROUTINE FLUSHONWRITE(LUNIT)                                    
C                                                                       
C        THE AIX SYSTEM'S FLUSH_ CALL WORKS ONLY ON UNIT 6, WHICH IS    
C        CONSIDERED CONNECTED TO A TERMINAL, BUT NOT OTHER UNITS.       
C                                                                       
      INCLUDE "fiosetup_.h"                                             
      INTEGER*4 LTEMP, IGOOD, IRESULT, FIOSETUP_                        
      COMMON /MACHSW/ KDIAG,ICORFL,IXDR,modio,mem10,lpnt10,mem10m       
      if(iand(modio,1).ne.0) return                                     
      LTEMP = LUNIT                                                     
      IGOOD = 0                                                         
      IRESULT = FIOSETUP_(LTEMP,                                        
     *                    IO_CMD_FLUSH_AFTER_WRITE,                     
     *                    IO_ARG_FLUSH_YES)                             
      IF (IRESULT.EQ.IGOOD) RETURN                                      
      WRITE(6,9) IRESULT                                                
    9 FORMAT(1X,'FLUSHONWRITE: ERROR, FIOSETUP_ RETURNED IRESULT=',I5)  
      CALL ABRT                                                         
      RETURN                                                            
      END                                                               
C*MODULE UNPORT  *DECK IGETGRDVECLEN                                    
      INTEGER FUNCTION IGETGRDVECLEN(MAXVEC)                            
C                                                                       
C        SCALAR MACHINES SHOULD AVOID SPECIAL VECTOR ROUTINES           
C        IN THE GRADIENT PROGRAM BY THE FOLLOWING COMPUTATION.          
C                                                                       
      IGETGRDVECLEN = MAXVEC*3+1                                        
C                                                                       
C        VECTOR MACHINES MAY WISH TO EXPERIMENT WITH THE VECTOR         
C        LENGTH AT WHICH TO CALL THE SPECIAL GRADIENT ROUTINES.         
C        NOTE THAT THIS IS NOT INTENDED TO BE A GENERAL ROUTINE         
C        FOR DECIDING ON PIPELINE LENGTHS, BUT RATHER AS THE NAME       
C        OF THE ROUTINE IMPLIES, IS SPECIFIC TO THE GRADIENT CODES.     
C                                                                       
C        THE CRAY-X1 VALUE WAS SELECTED BY TED PACKWOOD IN 2003         
C                                                                       
*CRY  IGETGRDVECLEN = 24                                                
*CX1  IGETGRDVECLEN = 4                                                 
*FUJ  IGETGRDVECLEN = 24                                                
*IBM  IGETGRDVECLEN = 24                                                
*NEC  IGETGRDVECLEN = 24                                                
      RETURN                                                            
      END                                                               
C*MODULE UNPORT  *DECK LOGAND                                           
      INTEGER FUNCTION LOGAND(IVAL,JVAL)                                
      INTEGER IVAL,JVAL                                                 
C                                                                       
C     RETURN LOGICAL AND OF ALL BITS IN THE INTEGERS IVAL, JVAL         
C                                                                       
*CRY  INTEGER AND                                                       
*CX1  INTEGER AND                                                       
*DEC  INTEGER AND                                                       
*L32  INTEGER AND                                                       
*L64  INTEGER AND                                                       
*NEC  INTEGER AND                                                       
*SGI  INTEGER AND                                                       
*SUN  INTEGER AND                                                       
*T3E  INTEGER AND                                                       
*W32  INTEGER AND                                                       
*W64  INTEGER AND                                                       
C                                                                       
      LOGAND = IAND(IVAL,JVAL)                                          
*CRY  LOGAND =  AND(IVAL,JVAL)                                          
*CX1  LOGAND =  AND(IVAL,JVAL)                                          
*DEC  LOGAND =  AND(IVAL,JVAL)                                          
*F77  LOGAND = 0                                                        
*F77  IF (IVAL .GE. JVAL) LOGAND = 1                                    
*FUJ  LOGAND = IAND(IVAL,JVAL)                                          
*HIT  LOGAND = IAND(IVAL,JVAL)                                          
*HP   LOGAND = IAND(IVAL,JVAL)                                          
*IBM  LOGAND = IAND(IVAL,JVAL)                                          
*INT  LOGAND = IAND(IVAL,JVAL)                                          
*L32  LOGAND =  AND(IVAL,JVAL)                                          
*L64  LOGAND =  AND(IVAL,JVAL)                                          
*NEC  LOGAND =  AND(IVAL,JVAL)                                          
*SGI  LOGAND =  AND(IVAL,JVAL)                                          
*SUN  LOGAND =  AND(IVAL,JVAL)                                          
*T3E  LOGAND =  AND(IVAL,JVAL)                                          
*VMS  LOGAND = IAND(IVAL,JVAL)                                          
*XT3  LOGAND =  AND(IVAL,JVAL)                                          
*W32  LOGAND =  AND(IVAL,JVAL)                                          
*W64  LOGAND =  AND(IVAL,JVAL)                                          
*I64  LOGAND = IAND(IVAL,JVAL)                                          
      RETURN                                                            
      END                                                               
C                                                                       
C*MODULE UNPORT  *DECK LOGOR                                            
      INTEGER FUNCTION LOGOR(IVAL,JVAL)                                 
      INTEGER IVAL,JVAL                                                 
C                                                                       
C     RETURN LOGICAL OR OF ALL BITS IN THE INTEGERS IVAL, JVAL          
C                                                                       
      INTEGER OR                                                        
*CRY  INTEGER OR                                                        
*CX1  INTEGER OR                                                        
*DEC  INTEGER OR                                                        
*L32  INTEGER OR                                                        
*L64  INTEGER OR                                                        
*NEC  INTEGER OR                                                        
*SGI  INTEGER OR                                                        
*SUN  INTEGER OR                                                        
*T3E  INTEGER OR                                                        
*W32  INTEGER OR                                                        
*W64  INTEGER OR                                                        
C                                                                       
      LOGOR =  OR(IVAL,JVAL)                                            
*CRY  LOGOR =  OR(IVAL,JVAL)                                            
*CX1  LOGOR =  OR(IVAL,JVAL)                                            
*DEC  LOGOR =  OR(IVAL,JVAL)                                            
*F77  LOGOR = 0                                                         
*F77  IF (IVAL.NE.0.AND.JVAL.NE.0) LOGOR = 1                            
*FUJ  LOGOR = IOR(IVAL,JVAL)                                            
*HIT  LOGOR = IOR(IVAL,JVAL)                                            
*HP   LOGOR = IOR(IVAL,JVAL)                                            
*IBM  LOGOR = IOR(IVAL,JVAL)                                            
*INT  LOGOR = IOR(IVAL,JVAL)                                            
*L32  LOGOR =  OR(IVAL,JVAL)                                            
*L64  LOGOR =  OR(IVAL,JVAL)                                            
*NEC  LOGOR =  OR(IVAL,JVAL)                                            
*SGI  LOGOR =  OR(IVAL,JVAL)                                            
*SUN  LOGOR =  OR(IVAL,JVAL)                                            
*T3E  LOGOR =  OR(IVAL,JVAL)                                            
*VMS  LOGOR = IOR(IVAL,JVAL)                                            
*XT3  LOGOR =  OR(IVAL,JVAL)                                            
*W32  LOGOR =  OR(IVAL,JVAL)                                            
*W64  LOGOR =  OR(IVAL,JVAL)                                            
*I64  LOGOR = IOR(IVAL,JVAL)                                            
      RETURN                                                            
      END                                                               
C                                                                       
C*MODULE UNPORT  *DECK LOGSHF                                           
      INTEGER FUNCTION LOGSHF(IVAL,JVAL)                                
      LOGSHF = ISHFT(IVAL,JVAL)                                         
      RETURN                                                            
      END                                                               
C                                                                       
C*MODULE UNPORT  *DECK LOGXOR                                           
      INTEGER FUNCTION LOGXOR(IVAL,JVAL)                                
      INTEGER IVAL,JVAL                                                 
C                                                                       
C     RETURN LOGICAL XOR OF ALL BITS IN THE INTEGERS IVAL, JVAL         
C                                                                       
      INTEGER XOR                                                       
*CRY  INTEGER XOR                                                       
*CX1  INTEGER XOR                                                       
*DEC  INTEGER XOR                                                       
*L32  INTEGER XOR                                                       
*L64  INTEGER XOR                                                       
*NEC  INTEGER XOR                                                       
*SGI  INTEGER XOR                                                       
*SUN  INTEGER XOR                                                       
*T3E  INTEGER XOR                                                       
*W32  INTEGER XOR                                                       
*W64  INTEGER XOR                                                       
C                                                                       
      LOGXOR =  XOR(IVAL,JVAL)                                          
*CRY  LOGXOR =  XOR(IVAL,JVAL)                                          
*CX1  LOGXOR =  XOR(IVAL,JVAL)                                          
*DEC  LOGXOR =  XOR(IVAL,JVAL)                                          
*F77  LOGXOR = 0                                                        
*F77  IF (IVAL.NE.JVAL) LOGXOR = 1                                      
*FUJ  LOGXOR =  XOR(IVAL,JVAL)                                          
*HIT  LOGXOR = IEOR(IVAL,JVAL)                                          
*HP   LOGXOR = IXOR(IVAL,JVAL)                                          
*IBM  LOGXOR = IEOR(IVAL,JVAL)                                          
*INT  LOGXOR = IEOR(IVAL,JVAL)                                          
*L32  LOGXOR =  XOR(IVAL,JVAL)                                          
*L64  LOGXOR =  XOR(IVAL,JVAL)                                          
*NEC  LOGXOR =  XOR(IVAL,JVAL)                                          
*SGI  LOGXOR =  XOR(IVAL,JVAL)                                          
*SUN  LOGXOR =  XOR(IVAL,JVAL)                                          
*T3E  LOGXOR =  XOR(IVAL,JVAL)                                          
*VMS  LOGXOR = IEOR(IVAL,JVAL)                                          
*XT3  LOGXOR =  XOR(IVAL,JVAL)                                          
*W32  LOGXOR =  XOR(IVAL,JVAL)                                          
*W64  LOGXOR =  XOR(IVAL,JVAL)                                          
*I64  LOGXOR = IEOR(IVAL,JVAL)                                          
      RETURN                                                            
      END                                                               
C                                                                       
C*MODULE UNPORT  *DECK PARSET                                           
      SUBROUTINE PARSET                                                 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      LOGICAL DIRTRF                                                    
      COMMON /TRFOPT/ CUTTRF,NWDTRF,MPTRAN,ITRFAO,NOSYMT,IPURTF,DIRTRF  
C                                                                       
C     ----- SET CONTROL FOR PARALLEL COMPUTATION -----                  
C                                                                       
C     ITRFAO CONTROLS AO INTEGRAL DISK STORAGE DURING CONVENTIONAL      
C     INTEGRAL TRANSFORMATIONS, 1=DUPLICATE AO LIST IF COMMUNICATION    
C     SPEED IS VERY POOR, E.G. ETHERNET, 2=DISTRIBUTE AO LIST BECAUSE   
C     BROADCAST OF INTEGRALS MAY BE MORE EFFICIENT THAN DISK I/O, E.G.  
C     THE SWITCH IN A SP2 MACHINE.  THIS VARIABLE CAN BE OVERRIDDEN     
C     LATER BY INPUT IN $MP2 OR $TRANS GROUPS.                          
C                                                                       
C     IBM AND COMPAQ/DIGITAL WORKSTATIONS WITH SLOWER NETWORKS WILL     
C     PREFER TO USE DUPLICATED AO INTEGRAL LISTS IN THE INTEGRAL        
C     TRANSFORMATIONS.  HOWEVER, THE IBM SP OR THE COMPAQ SUPERCLUSTER  
C     WITH BETTER COMMUNICATIONS MIGHT LIKE TO HAVE THE DEFAULT BE      
C     DUPLICATED AO INTEGRAL LISTS.  ANY SUCH MACHINE CAN USE A         
C     "SED HACK" TO SELECT THE *TRF LINE BELOW.                         
C                                                                       
      ITRFAO=1                                                          
C                                                                       
*CX1  ITRFAO= 2                                                         
*FUJ  ITRFAO= 2                                                         
*HIT  ITRFAO= 2                                                         
*NEC  ITRFAO= 2                                                         
*T3E  ITRFAO= 2                                                         
*TRF  ITRFAO= 2                                                         
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE UNPORT  *DECK SETFM                                            
      SUBROUTINE SETFM(IPAR)                                            
C-----------------------------------------------------------------------
C     ----- FAST MEMORY (FM) MANAGEMENT ROUTINES -----                  
C     NOTE THAT THERE ARE SIX ROUTINES, NAMELY                          
C        SETFM AND BIGFM ARE CALLED ONLY ONCE EACH,                     
C             TO OBTAIN AND RELEASE THE TOTAL MEMORY POOL.  IN UNIX,    
C             THESE NORMALLY WIND UP CALLING 'MALLOC' AND 'FREE'.       
C        VALFM, GETFM, RETFM ARE USED TOGETHER, AS TRIPLETS OF CALLS:   
C           VALFM TELLS THE LOCATION OF THE AVAILABLE MEMORY            
C           GETFM ASKS TO USE SOME OF THAT MEMORY                       
C           RETFM GIVES THE MEMORY BACK WHEN FINISHED                   
C              THESE TRIPLETS CAN BE NESTED, SO LONG AS THE CALLS       
C              TO RETFM ARE IN REVERSE ORDER OF THE CALLS TO GETFM,     
C              THAT IS TO SAY, THE IMPLEMENTATION IS A STACK.           
C        GOTFM MAY BE USED TO LEARN THE TOTAL AMOUNT OF AVAILABLE       
C              MEMORY, E.G. NOT YET COMMITTED BY GETFM'S.               
C                                                                       
C     ALL MEMORY CALCULATIONS HERE ARE TO BE BASED ON -WORDS-,          
C     WHERE A WORD IS DEFINED AS BEING A 64 BIT QUANTITY.               
C                                                                       
C     SETFM(IPAR) - ON ENTRY, IPAR IS MAXIMUM MEMORY DESIRED BY USER.   
C                   ON EXIT, IPAR IS MAXIMUM MEMORY ACTUALLY AVAILABLE. 
C                   ALLOCATES THE MEMORY POOL FROM THE SYSTEM.          
C                                                                       
C     BIGFM(IPAR) - ON EXIT, IPAR IS MAXIMUM MEMORY EVER USED.          
C                   FREES THE ENTIRE MEMORY POOL TO THE SYSTEM.         
C                                                                       
C     VALFM(IPAR) - ON EXIT, IPAR IS AN OFFSET (LTOP) TO THE HIGHEST    
C                   POSITION IN MEMORY CURRENTLY USED.                  
C                                                                       
C     GETFM(IPAR) - ON ENTRY, IPAR IS THE MEMORY PIECE REQUESTED.       
C                   IT IS ALLOCATED AT LTOP+1, AFTER WHICH THE          
C                   LTOP POINTER IS AUTOMATICALLY ADJUSTED UPWARDS.     
C                                                                       
C     RETFM(IPAR) - ON ENTRY, IPAR IS THE MEMORY PIECE TO BE FREED.     
C                   THE LTOP POINTER IS AUTOMATICALLY ADJUSTED DOWN.    
C                                                                       
C     GOTFM(IPAR) - ON EXIT, IPAR IS THE CURRENTLY UNUSED MEMORY.       
C                                                                       
C--------------------------------------------------------------------   
C                                                                       
C     A TYPICAL MEMORY REQUEST SEQUENCE IS VALFM,GETFM,...,RETFM.       
C     IT IS CRUCIAL THAT RETFM'S BE IN INVERSE ORDER OF GETFM CALLS.    
C     THE DYNAMIC POOL DURING THE MIDDLE ELLIPSIS OF THE SEQUENCE       
C                                                                       
C         CALL VALFM,GETFM,...VALFM,GETFM,...RETFM,...RETFM             
C     WITH ARG=      NEED1          NEED2    NEED2    NEED1             
C                                                                       
C     LOOKS LIKE THIS                                                   
C                                                                       
C            .....RESERVED........      AVAILABLE                       
C            <-------><----------><---------------------->              
C            X  NEED1    NEED2   Y                       Z              
C                                                                       
C     WHERE X=LOFFS, Y=LTOP, Z=MEMLIM                                   
C                                                                       
C     THERE ARE TWO IMPLEMENTATION STRATEGIES BELOW.                    
C                                                                       
C     ONE APPROACH IS TO DECLARE A LARGE FIXED LENGTH ARRAY             
C     X(MEMSIZ) JUST BELOW.  THE USER REQUESTED MEMORY -MEMLIM-         
C     CANNOT EXCEED THE "STATIC" SIZE -MEMSIZ- UNLESS -MEMSIZ-          
C     IS INCREASED HERE IN THE SOURCE, AND THE CODE RECOMPILED.         
C                                                                       
C     THE OTHER ALLOCATES ONLY A SINGLE WORD X(1), AND ACTUALLY         
C     ALLOCATES THE POOL ELSEWHERE.  -LOFFS- MUST BE THE DISTANCE       
C     IN WORDS FROM X(1) TO WHEREVER THIS POOL IS ALLOCATED.            
C     THIS STRATEGY IS TRULY "DYNAMIC", THE CALL TO SETFM CAN REQUEST   
C     ANY AMOUNT UP TO WHATEVER LIMITS THE OPERATING SYSTEM IMPOSES.    
C     IN THIS CASE, -MEMSIZ- IS JUST A DEFAULT, NOT AN UPPER BOUND.     
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C                                                                       
      LOGICAL GOPARR,DSKWRK,MASWRK                                      
C                                                                       
      COMMON /FMPARM/ LTOP,LOFFS,LENHI,LOCMEM,MEMLIM,MEMOK,nalign       
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)              
      COMMON /MACHIN/ NWDVAR,MAXFM,MAXSM,LIMFM,LIMSM                    
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
C                                                                       
      PARAMETER (MEMSIZ= 1 000 000)                                     
      COMMON /FMCOM / X(1)                                              
C                                                                       
*CRY  PARAMETER (MEMSIZ= 1 000 000)                                     
*CRY  COMMON /FMCOM / X(1)                                              
C                                                                       
*CX1  PARAMETER (MEMSIZ= 1 000 000)                                     
*CX1  COMMON /FMCOM / X(1)                                              
C                                                                       
*DEC  PARAMETER (MEMSIZ= 1 000 000)                                     
*DEC  COMMON /FMCOM / X(1)                                              
*DEC  EXTERNAL LOC                                                      
C                                                                       
*F77  PARAMETER (MEMSIZ= 1 000 000)                                     
*F77  COMMON /FMCOM / X(MEMSIZ)                                         
C                                                                       
*FUJ  PARAMETER (MEMSIZ= 1 000 000)                                     
*FUJ  COMMON /FMCOM / X(1)                                              
C                                                                       
*HIT  PARAMETER (MEMSIZ= 1 000 000)                                     
*HIT  COMMON /FMCOM / X(1)                                              
C                                                                       
*HP   PARAMETER (MEMSIZ= 1 000 000)                                     
*HP   COMMON /FMCOM / X(1)                                              
C                                                                       
*IBM  PARAMETER (MEMSIZ= 1 000 000)                                     
*IBM  COMMON /FMCOM / X(MEMSIZ)                                         
C                                                                       
*INT  PARAMETER (MEMSIZ= 1 000 000)                                     
*INT  COMMON /FMCOM / X(1)                                              
C                                                                       
*L32  PARAMETER (MEMSIZ= 1 000 000)                                     
*L32  COMMON /FMCOM / X(1)                                              
C                                                                       
*L64  PARAMETER (MEMSIZ= 1 000 000)                                     
*L64  COMMON /FMCOM / X(1)                                              
C                                                                       
*NEC  PARAMETER (MEMSIZ= 1 000 000)                                     
*NEC  COMMON /FMCOM / X(1)                                              
C                                                                       
*SGI  PARAMETER (MEMSIZ= 1 000 000)                                     
*SGI  COMMON /FMCOM / X(1)                                              
C                                                                       
*SUN  PARAMETER (MEMSIZ= 1 000 000)                                     
*SUN  COMMON /FMCOM / X(1)                                              
C                                                                       
*T3E  PARAMETER (MEMSIZ= 1 000 000)                                     
*T3E  COMMON /FMCOM / X(1)                                              
C                                                                       
*VMS  PARAMETER (MEMSIZ= 1 000 000)                                     
*VMS  COMMON /FMCOM / X(1)                                              
C                                                                       
*XT3  PARAMETER (MEMSIZ= 1 000 000)                                     
*XT3  COMMON /FMCOM / X(1)                                              
C                                                                       
*W32  PARAMETER (MEMSIZ= 1 000 000)                                     
*W32  COMMON /FMCOM / X(1)                                              
C                                                                       
*W64  PARAMETER (MEMSIZ= 1 000 000)                                     
*W64  COMMON /FMCOM / X(1)                                              
C                                                                       
*I64  PARAMETER (MEMSIZ= 1 000 000)                                     
*I64  COMMON /FMCOM / X(1)                                              
C                                                                       
C     ----- INITIALIZE FAST MEMORY (FM) MANAGEMENT -----                
C     -LOCX- IS THE ADDRESS OF THE FIRST ELEMENT OF X.                  
C     FOR STATIC IMPLEMENTATIONS, THIS SHOULD BE JUST 1.                
C                                                                       
      LOCX = 1                                                          
      LOCX = LADDRS(X)                                                  
*CRY  LOCX = LOC(X)                                                     
*CX1  LOCX = LOC(X)                                                     
*DEC  LOCX = LOC(X)                                                     
*FUJ  LOCX = LADDRS(X)                                                  
*HIT  LOCX = LADDRS(X)                                                  
*HP   LOCX = LOC(X)                                                     
*INT  LOCX = LOC(X)                                                     
*L32  LOCX = LADDRS(X)                                                  
*L64  LOCX = LOC(X)                                                     
*NEC  LOCX = LOC(X)                                                     
*NEC  CALL LADDRS(X,LOCX)                                               
*SGI  LOCX = LOC(X)                                                     
*SUN  LOCX = LOC(X)                                                     
*T3E  LOCX = LOC(X)                                                     
*VMS  LOCX = %LOC(X)                                                    
*XT3  LOCX = LOC(X)                                                     
*W32  LOCX = LADDRS(X)                                                  
*W64  LOCX = LOC(X)                                                     
*I64  LOCX = LOC(X)                                                     
C                                                                       
C     -MEMLIM- IS THE DESIRED MAXIMUM SIZE OF THE MEMORY POOL.          
C                                                                       
      MEMLIM = IPAR                                                     
      IF (MEMLIM .LE. 0) MEMLIM = MEMSIZ                                
C                                                                       
C     ----- CHECK FOR AVAILABILITY OF MEMORY -----                      
C     -MEMSIZ- IS AN ABSOLUTE BOUND FOR STATIC IMPLEMENTATIONS          
C                                                                       
      MEMAVL = MEMLIM                                                   
*F77  MEMAVL = MEMSIZ                                                   
*IBM  MEMAVL = MEMSIZ                                                   
C                                                                       
      IF (MEMLIM .GT. MEMAVL) THEN                                      
         IF (MASWRK)                                                    
     *   WRITE(IW,FMT='('' REQUESTED AMOUNT OF MEMORY ('',I10,          
     *      '' WORDS) EXCEEDS THE AMOUNT AVAILABLE ('',I10,'')'')')     
     *               MEMLIM,MEMAVL                                      
         CALL ABRT                                                      
      END IF                                                            
C                                                                       
C     ----- ALLOCATE THE DYNAMIC MEMORY POOL -----                      
C     -LOCMEM- IS THE STARTING ADDRESS OF THE DYNAMIC POOL.             
C     FOR STATIC IMPLEMENTATIONS, THIS SHOULD JUST BE 1.                
C                                                                       
      LOCMEM = 1                                                        
      LOCMEM = MEMGET(MEMLIM)                                           
*CRY  LOCMEM = MEMGET(MEMLIM)                                           
*CX1  LOCMEM = MEMGET(MEMLIM)                                           
*DEC  LOCMEM = MEMGET(MEMLIM)                                           
*FUJ  LOCMEM = MEMGET(MEMLIM)                                           
*HIT  LOCMEM = MEMGET(MEMLIM)                                           
*HP   LOCMEM = MALLOC((MEMLIM+2)*8)                                     
*INT  LOCMEM = MEMGET(MEMLIM)                                           
*L32  LOCMEM = MEMGET(MEMLIM)                                           
*L64  LOCMEM = MEMGET(MEMLIM)                                           
*NEC  CALL MEMGET(MEMLIM,LOCMEM)                                        
*SGI  LOCMEM = MEMGET(MEMLIM)                                           
*SUN  LOCMEM = MEMGET(MEMLIM)                                           
*T3E  LOCMEM = MEMGET(MEMLIM)                                           
*VMS  ISTAT = LIB$GET_VM(8*MEMLIM+MOD(LOCX,8),LOCMEM)                   
*VMS  IF (.NOT.ISTAT) CALL LIB$STOP(%VAL(ISTAT))                        
*XT3  LOCMEM = MEMGET(MEMLIM)                                           
*W32  LOCMEM = MEMGET(MEMLIM)                                           
*W64  LOCMEM = MEMGET(MEMLIM)                                           
*I64  LOCMEM = MEMGET(MEMLIM)                                           
c     if(iand(LOCMEM,4).ne.0) write(6,*) 'mempnt=',LOCMEM               
c     check alignment                                                   
C                                                                       
      IF (LOCMEM .EQ. 0) THEN                                           
         IF (MASWRK) WRITE(IW,*) MEMLIM,' WORDS OF MEMORY UNAVAILABLE'  
         CALL ABRT                                                      
      END IF                                                            
C                                                                       
C         COMPUTE THE OFFSET -LOFFS- FROM X(1) TO THE BEGINNING         
C         OF THE DYNAMIC POOL (E.G. BYTE TO WORD CONVERSION).           
C         -LOFFS- WILL BE ZERO FOR STATIC IMPLEMENTATIONS.              
C                                                                       
      LOFFS = LOCMEM - LOCX                                             
      LOFFS = (LOFFS+7)/8 + 1                                           
*CRY  LOFFS =  LOFFS                                                    
*CX1  LOFFS = (LOFFS+7)/8 + 1                                           
*DEC  LOFFS = (LOFFS+7)/8 + 1                                           
*FUJ  LOFFS = (LOFFS+7)/8 + 1                                           
*HIT  LOFFS = (LOFFS+7)/8 + 1                                           
*HP   LOFFS = (LOFFS+7)/8 + 1                                           
*INT  LOFFS = (LOFFS+7)/8 + 1                                           
*L32  LOFFS = (LOFFS+7)/8 + 1                                           
*L64  LOFFS = (LOFFS+7)/8 + 1                                           
*NEC  LOFFS = (LOFFS+7)/8 + 1                                           
*SGI  LOFFS = (LOFFS+7)/8 + 1                                           
*SUN  LOFFS = (LOFFS+7)/8 + 1                                           
*T3E  LOFFS = (LOFFS+7)/8 + 1                                           
*VMS  LOFFS = (LOFFS+7)/8 + 1                                           
*XT3  LOFFS = (LOFFS+7)/8 + 1                                           
*W32  LOFFS = (LOFFS+7)/8 + 1                                           
*W64  LOFFS = (LOFFS+7)/8 + 1                                           
*I64  LOFFS = (LOFFS+7)/8 + 1                                           
C                                                                       
      LTOP  = 0                                                         
      LENHI = 0                                                         
      MEMOK = 99                                                        
C                                                                       
C     ----- REPORT MEMORY AVAILABILITY -----                            
C                                                                       
      MEMEXT = MEMAVL - MEMLIM                                          
      IF (MASWRK) THEN                                                  
      IF (MEMEXT .NE. 0) THEN                                           
         WRITE(IW,*) MEMLIM,' WORDS REQUESTED ',                        
     *               MEMAVL,' WORDS AVAILABLE ',                        
     *               MEMEXT,' WORDS REMAIN FREE'                        
      ELSE                                                              
         WRITE(IW,FMT='(1X,I10,'' WORDS OF MEMORY AVAILABLE'')') MEMLIM 
      END IF                                                            
      END IF                                                            
C                                                                       
C     SET /MACHIN/ VALUES.                                              
C        MAXFM/SM  = TOTAL FAST/SLOW MEMORY IN DYNAMIC POOL.            
C        LIMFM/SM  = TOP ADDRESS OF THE FAST/SLOW DYNAMIC POOL.         
C     GAMESS DOES NOT CURRENTLY USE SLOW MEMORY.                        
C                                                                       
      MAXFM = MEMLIM                                                    
      LIMFM = MEMLIM + LOFFS                                            
      MAXSM = 0                                                         
      LIMSM = 0                                                         
C                                                                       
      IPAR = MEMAVL                                                     
      RETURN                                                            
      END                                                               
C*MODULE UNPORT  *DECK BIGFM                                            
      SUBROUTINE BIGFM(IPAR)                                            
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      LOGICAL GOPARR,DSKWRK,MASWRK                                      
      COMMON /FMPARM/ LTOP,LOFFS,LENHI,LOCMEM,MEMLIM,MEMOK,nalign       
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)              
C                                                                       
C     ----- RETURN MAXIMUM MEMORY USED -----                            
C                                                                       
      IF (MEMOK .NE. 99) RETURN                                         
      MEMOK = 0                                                         
      IPAR = LENHI                                                      
      IF (MASWRK) WRITE(IW,*) IPAR,' WORDS OF DYNAMIC MEMORY USED'      
      CALL FLSHBF(IW)                                                   
      CALL FLSHBF(IP)                                                   
C                                                                       
C        ----- RELEASE MEMORY TO SYSTEM -----                           
C     STATIC IMPLEMENTATIONS SHOULDN'T DO ANYTHING HERE.                
C                                                                       
      CALL MEMREL(LOCMEM)                                               
*CRY  CALL MEMREL(LOCMEM)                                               
*CX1  CALL MEMREL(LOCMEM)                                               
CDEC  CALL MEMREL(LOCMEM)    NOTE: ALPHA CORE DUMPS IF CALLED           
*FUJ  CALL MEMREL(LOCMEM)                                               
*HIT  CALL MEMREL(LOCMEM)                                               
*HP   CALL FREE(LOCMEM)                                                 
*INT  CALL MEMREL(LOCMEM)                                               
*L32  CALL MEMREL(LOCMEM)                                               
*L64  CALL MEMREL(LOCMEM)                                               
*NEC  CALL MEMREL(LOCMEM)                                               
*SGI  CALL MEMREL(LOCMEM)                                               
*SUN  CALL MEMREL(LOCMEM)                                               
*T3E  CALL MEMREL(LOCMEM)                                               
*VMS  ISTAT = LIB$FREE_VM(8*MEMLIM,LOCMEM,)                             
*VMS  IF (.NOT.ISTAT) CALL LIB$STOP(%VAL(ISTAT))                        
*XT3  CALL MEMREL(LOCMEM)                                               
*W32  CALL MEMREL(LOCMEM)                                               
*W64  CALL MEMREL(LOCMEM)                                               
*I64  CALL MEMREL(LOCMEM)                                               
      RETURN                                                            
      END                                                               
C*MODULE UNPORT  *DECK VALFM                                            
      SUBROUTINE VALFM(IPAR)                                            
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      COMMON /FMPARM/ LTOP,LOFFS,LENHI,LOCMEM,MEMLIM,MEMOK,nalign       
C                                                                       
C     ----- RETURN THE CURRENT TOP OF FM ARRAY -----                    
C                                                                       
      IPAR = LTOP + LOFFS                                               
      RETURN                                                            
      END                                                               
C*MODULE UNPORT  *DECK GETFM                                            
      SUBROUTINE GETFM(IPAR0)                                           
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      LOGICAL GOPARR,DSKWRK,MASWRK                                      
      COMMON /FMPARM/ LTOP,LOFFS,LENHI,LOCMEM,MEMLIM,MEMOK,nalign       
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)              
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
C                                                                       
C         UNCOMMENT NEXT LINE AND ALSO IN -RETFM- FOR PLUMBING JOBS.    
CLEAK IF(MASWRK) WRITE(IW,*) 'GETFM: ALLOCATING',IPAR                   
C                                                                       
C     ----- RESERVE IPAR WORDS OF FM -----                              
C                                                                       
      IPAR=((IPAR0-1)/nalign+1)*nalign                                  
c     write(6,*) 'wwwa',IPAR0,IPAR,nalign                               
      LTOP = LTOP + IPAR                                                
      IF (LTOP .LE. MEMLIM ) THEN                                       
         LENHI = MAX(LENHI,LTOP)                                        
      ELSE                                                              
         WRITE(IW,9000) ME,LTOP,MEMLIM                                  
         CALL ABRT                                                      
      END IF                                                            
      RETURN                                                            
C                                                                       
 9000 FORMAT(1X,'***** ERROR: MEMORY REQUEST EXCEEDS AVAILABLE MEMORY'/ 
     *    1X,'PROCESS NO.',I5,' WORDS REQUIRED=',I10,' AVAILABLE=',I10) 
      END                                                               
C*MODULE UNPORT  *DECK RETFM                                            
      SUBROUTINE RETFM(IPAR0)                                           
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      COMMON /FMPARM/ LTOP,LOFFS,LENHI,LOCMEM,MEMLIM,MEMOK,nalign       
C                                                                       
C         UNCOMMENT NEXT LINES AND ALSO IN -GETFM- FOR PLUMBING JOBS.   
CLEAK LOGICAL GOPARR,DSKWRK,MASWRK                                      
CLEAK COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)              
CLEAK COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
CLEAK IF(MASWRK) WRITE(IW,*) 'RETFM:  RETURNING',IPAR                   
C                                                                       
C     ----- RETURN IPAR WORDS OF FM -----                               
C                                                                       
      IPAR=((IPAR0-1)/nalign+1)*nalign                                  
c     write(6,*) 'wwwb',IPAR0,IPAR,nalign                               
      LTOP = LTOP - IPAR                                                
      RETURN                                                            
      END                                                               
C*MODULE UNPORT  *DECK GOTFM                                            
      SUBROUTINE GOTFM(IPAR)                                            
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      COMMON /FMPARM/ LTOP,LOFFS,LENHI,LOCMEM,MEMLIM,MEMOK,nalign       
C                                                                       
C     ----- RETURN NUMBER OF FREE WORDS IN FM -----                     
C                                                                       
      IPAR0 = MEMLIM - LTOP                                             
      IPAR=((IPAR0-1)/nalign+1)*nalign-nalign                           
      RETURN                                                            
      END                                                               
C*MODULE UNPORT  *DECK TEXIT                                            
      SUBROUTINE TEXIT(NCALL,NREST)                                     
C                                                                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C                                                                       
      LOGICAL PACK,GOPARR,DSKWRK,MASWRK                                 
C                                                                       
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)              
      COMMON /RESTAR/ TIMLIM,IREST,NREC,INTLOC,IST,JST,KST,LST          
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK                 
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
      COMMON /INTFIL/ NINTMX,NHEX,NTUPL,PACK,INTTYP,IGRDTYP             
      COMMON /TMVALS/ TI,TX,TIM                                         
C$omp threadprivate(/TMVALS/)
C                                                                       
      CALL TIMIT(NCALL)                                                 
      IF (TIM .LT. TIMLIM) RETURN                                       
      IF (MASWRK) THEN                                                  
         WRITE(IW,9018)                                                 
         WRITE (IW,9008) TIMLIM,NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK,      
     *       NREST,IST,JST,KST,LST,NREC,INTLOC,NINTMX                   
      END IF                                                            
      CALL ABRT                                                         
      RETURN                                                            
C                                                                       
 9008 FORMAT(F10.0,11I3,I10,2I5)                                        
 9018 FORMAT(1X,'**** JOB HAS EXHAUSTED ITS CPU ALLOTMENT ****')        
      END                                                               
C*MODULE UNPORT  *DECK TIMIT                                            
      SUBROUTINE TIMIT(INDEX)                                           
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C                                                                       
      LOGICAL GOPARR,DSKWRK,MASWRK                                      
C                                                                       
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)              
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
      COMMON /TMVALS/ TI,TX,TIM                                         
C$omp threadprivate(/TMVALS/)
      COMMON /TIMING/ CPU,WALL                                          
      COMMON /SIMDAT/ NACC,NREJ,IGOMIN,NRPA,IBWM,NACCT,NREJT,NRPAT,     
     *                NPRTGO,IDPUNC,IGOFLG                              
C                                                                       
C        COMPUTE AND PRINT INTERVAL CPU TIME                            
C        THIS IS CALLED WITH INDEX=0 ONLY FROM THE MAIN PROGRAM         
C        ALL OTHER CALLS SHOULD PRINT THE INTERVAL TIME.                
C                                                                       
      CALL TSECND(TIM)                                                  
      TX = TIM-TI                                                       
      TI = TIM                                                          
      IF (INDEX .EQ. 0) RETURN                                          
C                                                                       
      RATIO=100.0D+00                                                   
      IF(WALL.GT.0.0D+00) RATIO = 100.0D+00*TIM/WALL                    
      TMINUT = TIM/60.0D+00                                             
      IF((IPTIM.GT.0.OR.MASWRK).AND.NPRTGO.NE.2) THEN                   
         IF(GOPARR) THEN                                                
            WRITE(IW,9010) ME,TX,TIM,TMINUT,WALL,RATIO                  
         ELSE                                                           
            WRITE(IW,9000) TX,TIM,TMINUT,WALL,RATIO                     
         END IF                                                         
      END IF                                                            
      CALL FLSHBF(IW)                                                   
      CALL FLSHBF(IP)                                                   
      RETURN                                                            
C                                                                       
 9000 FORMAT(1X,'STEP CPU TIME =',F9.2,                                 
     *          ' TOTAL CPU TIME =',F13.1,' (',F9.1,' MIN)'/            
     *       1X,'TOTAL WALL CLOCK TIME=',F13.1,                         
     *          ' SECONDS, CPU UTILIZATION IS',F9.2,'%')                
 9010 FORMAT(1X,'CPU',I6,': STEP CPU TIME=',F9.2,                       
     *          ' TOTAL CPU TIME=',F13.1,' (',F9.1,' MIN)'/             
     *       1X,'TOTAL WALL CLOCK TIME=',F13.1,                         
     *          ' SECONDS, CPU UTILIZATION IS',F9.2,'%')                
      END                                                               
C*MODULE UNPORT  *DECK TMDATE                                           
      SUBROUTINE TMDATE(TIMSTR)                                         
C-----------------------------------------------------------------------
C     ----- RETURN REAL ARRAY CONTAINING 24 CHARACTER TIME STAMP -----  
C           ANY FORMAT RESEMBLING "HH:MM:SS DD-MMM-YYYY" WILL DO.       
C-----------------------------------------------------------------------
C                                                                       
C            IBM RS/6000 RUNNING AIX                                    
      DOUBLE PRECISION TIMSTR(3)                                        
      CHARACTER*26 STAMP                                                
      CALL FDATE_(STAMP)                                                
      READ (UNIT=STAMP,FMT='(3A8)') TIMSTR                              
C                                                                       
C        THE ROUTINE 'FDATE' FIRST APPEARED IN IBM'S XL FORTRAN AT      
C        VERSION 2.3, RELEASED SEPTEMBER 1992.  IF YOU HAVE RECEIVED    
C        AN UNRESOLVED EXTERNAL REFERENCE FOR 'FDATE', YOUR COMPILER    
C        MUST BE AN ANTIQUE!  IN THIS CASE, YOU CAN USE OUR ORIGINAL    
C        CODE WRITTEN FOR XL FORTRAN 1.0 RELEASED APRIL 1990, BELOW:    
CAIX  DOUBLE PRECISION WALL,TIMSTR(3)                                   
CAIX  CHARACTER*24 STAMP                                                
CAIX  CHARACTER*3 MONTHS(12)                                            
CAIX  DIMENSION IDAYS(12)                                               
CAIX  INTEGER*4 IWALL                                                   
CAIX  DATA MONTHS/'JAN','FEB','MAR','APR','MAY','JUN',                  
CAIX *            'JUL','AUG','SEP','OCT','NOV','DEC'/                  
CAIX  DATA IDAYS/31,28,31, 30,31,30, 31,31,30, 31,30,31/                
C                                                                       
C     UGH!!!  AIX RETURNS THE NUMBER OF SECONDS SINCE 1-JAN-1970 GMT.   
C     (GMT IS 6 HOURS AHEAD OF AMES'S CST).                             
C                                                                       
CAIX  ISTAT = TIME(IWALL)                                               
CAIX  ITZONE = 6                                                        
CAIX  IWALL = IWALL - ITZONE*3600                                       
CAIX  WALL = IWALL                                                      
C         NLEAP IS LEAP YEARS AFTER 1970 --PRIOR-- TO THIS YEAR         
CAIX  NDAY = WALL/86400.0D+00                                           
CAIX  WALL = WALL - 86400.0D+00*NDAY                                    
CAIX  NYEAR = NDAY/365.25                                               
CAIX  NLEAP = (NYEAR+1)/4                                               
CAIX  WALL = IWALL - 86400.0D+00*(NYEAR*365 + NLEAP)                    
CAIX  NYEAR = NYEAR+1970                                                
CAIX  IF(MOD(NYEAR,4).EQ.0) IDAYS(2) = 29                               
C                                                                       
CAIX  NDAY = WALL/86400.0D+00                                           
CAIX  WALL = WALL - 86400.0D+00*NDAY                                    
CAIX  NDAY = NDAY+1                                                     
CAIX  DO 1 M=1,12                                                       
CAIX     MONTH = M                                                      
CAIX     NDAY = NDAY - IDAYS(MONTH)                                     
CAIX     IF(NDAY.LE.0) GO TO 2                                          
CAIX1 CONTINUE                                                          
CAIX2 CONTINUE                                                          
CAIX  NDAY = NDAY + IDAYS(MONTH)                                        
C                                                                       
CAIX  NHOUR = WALL/3600.0D+00                                           
CAIX  WALL = WALL - 3600.0D+00*NHOUR                                    
CAIX  NMIN = WALL/60.0D+00                                              
CAIX  WALL = WALL - 60.0D+00*NMIN                                       
CAIX  NSEC = WALL                                                       
CAIX  WRITE(UNIT=STAMP,FMT=5) NHOUR,NMIN,NSEC,NDAY,MONTHS(MONTH),NYEAR  
CAIX5 FORMAT(I2,':',I2.2,':',I2.2,' CST ',I2,'-',A3,'-',I4)             
CAIX  READ(UNIT=STAMP,FMT='(3A8)') TIMSTR                               
C-----------------------------------------------------------------------
C                                                                       
*CRY  DIMENSION TIMSTR(3)                                               
*CRY  DATA BLANK/8H        /                                            
*CRY  CALL CLOCK(TIMSTR(1))                                             
*CRY  TIMSTR(2) = BLANK                                                 
*CRY  CALL DATE(TIMSTR(3))                                              
C-----------------------------------------------------------------------
C                                                                       
*CX1  DIMENSION TIMSTR(3)                                               
*CX1  REAL*8 TTIMSTR(3)                                                 
*CX1  CHARACTER*8 DATESTR,TIMESTR                                       
*CX1  DATA BLANK/8H        /                                            
*CX1  EQUIVALENCE (TIMESTR,TTIMSTR(1))                                  
*CX1  EQUIVALENCE (BLANK,TTIMSTR(2))                                    
*CX1  EQUIVALENCE (DATESTR,TTIMSTR(3))                                  
*CX1  CALL CLOCK(TIMESTR)                                               
*CX1  CALL DATE(DATESTR)                                                
*CX1  TIMSTR(1) = TTIMSTR(1)                                            
*CX1  TIMSTR(2) = TTIMSTR(2)                                            
*CX1  TIMSTR(3) = TTIMSTR(3)                                            
C-----------------------------------------------------------------------
C                                                                       
*DEC  EXTERNAL FDATE                                                    
*DEC  DOUBLE PRECISION TIMSTR(3)                                        
*DEC  CHARACTER*24 FDATE,STAMP                                          
*DEC  STAMP = FDATE()                                                   
*DEC  READ (UNIT=STAMP,FMT='(3A8)') TIMSTR                              
C-----------------------------------------------------------------------
C                                                                       
*F77  DOUBLE PRECISION TIMSTR(3)                                        
*F77  CHARACTER*24 STAMP                                                
*F77  STAMP='HH:MM:SS DD-MMM-19YY    '                                  
*F77  READ(UNIT=STAMP,FMT='(3A8)') TIMSTR                               
C-----------------------------------------------------------------------
C                                                                       
*FUJ  EXTERNAL FDATE                                                    
*FUJ  DOUBLE PRECISION TIMSTR(3)                                        
*FUJ  CHARACTER*24 FDATE,STAMP                                          
*FUJ  STAMP = FDATE()                                                   
*FUJ  READ (UNIT=STAMP,FMT='(3A8)') TIMSTR                              
C-----------------------------------------------------------------------
C                                                                       
*HIT  DOUBLE PRECISION BLANK,TIMSTR(3)                                  
*HIT  CHARACTER*8 CLKSTR                                                
*HIT  DATA BLANK/8H        /                                            
C                                                                       
*HIT  CALL CLOCK(IS,2)                                                  
*HIT  IM=IS/60                                                          
*HIT  IH=IM/60                                                          
*HIT  IM=IM-60*IH                                                       
*HIT  IS=IS-60*IM-3600*IH                                               
*HIT  WRITE (UNIT=CLKSTR,FMT=1) IH,IM,IS                                
*HIT1 FORMAT(I2,':',I2,':',I2)                                          
*HIT  READ  (UNIT=CLKSTR,FMT='(1A8)') TIMSTR(1)                         
C                                                                       
*HIT  TIMSTR(2) = BLANK                                                 
*HIT  CALL DATE(TIMSTR(3))                                              
C-----------------------------------------------------------------------
C                                                                       
*HP   DOUBLE PRECISION TIMSTR(3)                                        
*HP   CHARACTER*24 STAMP                                                
*HP   CALL FDATE(STAMP)                                                 
*HP   READ (UNIT=STAMP,FMT='(3A8)') TIMSTR                              
C-----------------------------------------------------------------------
C                                                                       
*IBM  DOUBLE PRECISION TIMSTR(3)                                        
*IBM  CALL ZDATE(TIMSTR)                                                
C-----------------------------------------------------------------------
C                                                                       
*INT  DOUBLE PRECISION TIMSTR(3)                                        
*INT  CHARACTER*24 STAMP                                                
*INT  CALL FDATE(STAMP)                                                 
*INT  READ(UNIT=STAMP,FMT='(3A8)') TIMSTR                               
C-----------------------------------------------------------------------
C        FDATE IS IMPLEMENTED IN ZUNIX.C, NOT A SYSTEM LIBRARY CALL     
*L32  DOUBLE PRECISION TIMSTR(3)                                        
*L32  CHARACTER*24 STAMP                                                
*L32  CALL FDATE(STAMP)                                                 
*L32  READ(UNIT=STAMP,FMT='(3A8)') TIMSTR                               
C-----------------------------------------------------------------------
C        FDATE COMES FROM THE GFORTRAN LIBRARY                          
*L64  DOUBLE PRECISION TIMSTR(3)                                        
*L64  CHARACTER*24 STAMP                                                
*L64  CALL FDATE(STAMP)                                                 
*L64  READ(UNIT=STAMP,FMT='(3A8)') TIMSTR                               
C-----------------------------------------------------------------------
C                                                                       
*NEC  DOUBLE PRECISION TIMSTR(3)                                        
*NEC  CHARACTER*24 DATSTR                                               
*NEC  CALL FDATE(DATSTR)                                                
*NEC  READ (UNIT=DATSTR,FMT='(3A8)') TIMSTR                             
C-----------------------------------------------------------------------
C                                                                       
*SGI  DOUBLE PRECISION TIMSTR(3)                                        
*SGI  CHARACTER*24 FDATE                                                
*SGI  READ (UNIT=FDATE(),FMT='(3A8)') TIMSTR                            
C-----------------------------------------------------------------------
C                                                                       
*SUN  DOUBLE PRECISION TIMSTR(3)                                        
*SUN  CHARACTER*24 DATSTR                                               
*SUN  CALL FDATE(DATSTR)                                                
*SUN  READ (UNIT=DATSTR,FMT='(3A8)') TIMSTR                             
C-----------------------------------------------------------------------
C                                                                       
*T3E  DIMENSION TIMSTR(3)                                               
*T3E  DATA BLANK/8H        /                                            
*T3E  CALL CLOCK(TIMSTR(1))                                             
*T3E  TIMSTR(2) = BLANK                                                 
*T3E  CALL DATE(TIMSTR(3))                                              
C-----------------------------------------------------------------------
C                                                                       
*VMS  DIMENSION TIMSTR(6)                                               
*VMS  DATA BLANK/4H    /                                                
*VMS  CALL TIME(TIMSTR)                                                 
*VMS  TIMSTR(3)=BLANK                                                   
*VMS  TIMSTR(6)=BLANK                                                   
*VMS  CALL DATE(TIMSTR(4))                                              
C-----------------------------------------------------------------------
C                                                                       
*XT3  DOUBLE PRECISION TIMSTR(3)                                        
*XT3  CHARACTER*24 STAMP                                                
*XT3  CALL FDATE(STAMP)                                                 
*XT3  READ(UNIT=STAMP,FMT='(3A8)') TIMSTR                               
C-----------------------------------------------------------------------
C                                                                       
*W32  DOUBLE PRECISION TIMSTR(3)                                        
*W32  CHARACTER*24 STAMP                                                
*W32  CALL FDATE(STAMP)                                                 
*W32  READ (UNIT=STAMP,FMT='(3A8)') TIMSTR                              
C-----------------------------------------------------------------------
C                                                                       
*W64  DOUBLE PRECISION TIMSTR(3)                                        
*W64  CHARACTER*24 STAMP                                                
*W64  CALL FDATE(STAMP)                                                 
*W64  READ (UNIT=STAMP,FMT='(3A8)') TIMSTR                              
C-----------------------------------------------------------------------
C                                                                       
*I64  DOUBLE PRECISION TIMSTR(3)                                        
*I64  CHARACTER*24 STAMP                                                
*I64  CALL FDATE(STAMP)                                                 
*I64  READ(UNIT=STAMP,FMT='(3A8)') TIMSTR                               
C-----------------------------------------------------------------------
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE UNPORT  *DECK TSECND                                           
C>                                                                      
C>    @brief   do CPU and wall clock timings                            
C>                                                                      
C>    @details do CPU and wall clock timings through a name TSECND      
C>             that is not used on any operating system, in order       
C>             to call whatever name each kind of computer uses.        
C>                                                                      
C>    @param   TIM   return value is CPU time since job start           
C>    @param   CPU   CPU time,   since job start (in /TIMING/)          
C>    @param   WALL  wall clock, since job start (in /TIMING/)          
C>                                                                      
C>    @author  Mike Schmidt, many moons ago                             
C>                                                                      
      SUBROUTINE TSECND(TIM)                                            
C-----------------------------------------------------------------------
C       ----- THIS ROUTINE PERFORMS CPU AND WALL CLOCK TIMING -----     
C                                                                       
C       THIS ROUTINE SHOULD SET 'CPU' AND 'WALL' VARIABLES              
C       TO THE TOTAL ELAPSED CPU AND WALL CLOCK TIMES,                  
C       MEASURED IN SECONDS.  IN ADDITION, THE CALLING                  
C       ARGUMENT 'TIM' SHOULD BE SET EQUAL TO 'CPU'.                    
C                                                                       
C       ON THE FIRST ENTRY, 'CPU0' AND 'WALL0' SHOULD BE SET            
C       TO THE APPROPRIATE BASE VALUE ON JOB START.                     
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
C                       >>>>> UNIX NOTE <<<<<                           
C                                                                       
C     DEPENDING ON JUST WHAT UNIX YOU ARE USING, ETIME IS EITHER        
C     SINGLE OR DOUBLE PRECISION (DOUBLE ON THE CELERITY AND ALLIANT).  
C     SOME SYSTEMS ALSO DON'T RETURN THE SUM OF USER AND SYSTEM         
C     TIME AS THE FUNCTION VALUE.  THEREFORE, THE MOST PORTABLE         
C     WAY OF DOING THIS TIMING IS TO ADD THE USER AND SYSTEM TIMES      
C     RETURNED IN TARRAY, WHICH SEEMS TO BE SINGLE PRECISION ON         
C     ALL UNIX SYSTEMS THUS FAR TESTED.                                 
C                                                                       
C-----------------------------------------------------------------------
C      * * * *  IBM RUNNING AIX PORTION  * * * *                        
C                                                                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      REAL*4 ETIME_                                                     
      TYPE TB_TYPE                                                      
        SEQUENCE                                                        
        REAL*4 USRTIME                                                  
        REAL*4 SYSTIME                                                  
      END TYPE                                                          
      TYPE (TB_TYPE) ETIME_STRUCT                                       
      INTEGER*4 TIME_                                                   
      LOGICAL FIRST                                                     
      COMMON /TIMING/ CPU,WALL                                          
      SAVE FIRST,CPU0,WALL0                                             
      DATA FIRST/.TRUE./                                                
C                                                                       
C       ----- INITIALIZE CLOCKS -----                                   
C                                                                       
      IF(FIRST) THEN                                                    
         FIRST = .FALSE.                                                
         CPU0  = ETIME_(ETIME_STRUCT)                                   
         WALL0 = TIME_()                                                
      END IF                                                            
C                                                                       
C       ----- OBTAIN ELAPSED TIMES SINCE JOB STARTED -----              
C                                                                       
      CPU  = ETIME_(ETIME_STRUCT)                                       
      CPU  = CPU - CPU0                                                 
      WALL = TIME_()                                                    
      WALL = WALL - WALL0                                               
      TIM  = CPU                                                        
      RETURN                                                            
C                                                                       
C        THE ROUTINE 'ETIME' FIRST APPEARED IN IBM'S XL FORTRAN AT      
C        VERSION 2.3, RELEASED SEPTEMBER 1992.  IF YOU HAVE RECEIVED    
C        AN UNRESOLVED EXTERNAL REFERENCE FOR 'ETIME', YOUR COMPILER    
C        MUST BE AN ANTIQUE!  IN THIS CASE, YOU CAN USE OUR ORIGINAL    
C        CODE WRITTEN FOR XL FORTRAN 1.0 RELEASED APRIL 1990, BELOW:    
CAIX  DOUBLE PRECISION CPU0,WALL0,CPU,WALL,TIM                          
CAIX  LOGICAL FIRST                                                     
CAIX  INTEGER*4 BUF(4)                                                  
CAIX  INTEGER*4 IWALL                                                   
CAIX  PARAMETER (TCFCTR=100.0D+00)                                      
CAIX  COMMON /TIMING/ CPU,WALL                                          
CAIX  SAVE FIRST,CPU0,WALL0                                             
CAIX  DATA FIRST/.TRUE./                                                
C                                                                       
CAIX  IF(FIRST) THEN                                                    
CAIX     FIRST=.FALSE.                                                  
CAIX     ISTAT = TIME(IWALL)                                            
CAIX     WALL0 = IWALL                                                  
CAIX     CPU0 = 0.0D+00                                                 
CAIX  END IF                                                            
C                                                                       
CAIX  ISTAT = TIME(IWALL)                                               
CAIX  WALL = IWALL                                                      
CAIX  WALL = WALL-WALL0                                                 
CAIX  ISTAT = TIMES(BUF)                                                
CAIX  USER   = BUF(1)/TCFCTR                                            
CAIX  SYSTEM = BUF(2)/TCFCTR                                            
CAIX  CPU = USER+SYSTEM                                                 
CAIX  CPU = CPU-CPU0                                                    
CAIX  TIM = CPU                                                         
CAIX  RETURN                                                            
C-----------------------------------------------------------------------
C      * * * *  CRAY (UNICOS) SPECIFIC PORTION  * * * *                 
C                                                                       
*CRY  REAL CPU0,WALL0,CPU,WALL,TIM                                      
*CRY  LOGICAL FIRST                                                     
*CRY  COMMON /TIMING/ CPU,WALL                                          
*CRY  SAVE FIRST,CPU0,WALL0                                             
*CRY  DATA FIRST/.TRUE./                                                
C                                                                       
C       ----- INITIALIZE CLOCKS -----                                   
C                                                                       
*CRY  IF(FIRST) THEN                                                    
*CRY     FIRST=.FALSE.                                                  
*CRY     CPU0 = 0.0E+00                                                 
*CRY     WALL0 = TIMEF()                                                
*CRY     WALL0 = 0.001E+00*WALL0                                        
*CRY  END IF                                                            
C                                                                       
C       ----- OBTAIN ELAPSED TIMES SINCE JOB STARTED -----              
C                                                                       
*CRY  CALL SECOND(CPU)                                                  
*CRY  WALL = TIMEF()                                                    
*CRY  WALL = 0.001E+00 * WALL - WALL0                                   
*CRY  TIM = CPU                                                         
*CRY  RETURN                                                            
C                                                                       
C-----------------------------------------------------------------------
C      * * * *  CRAY X1 SPECIFIC PORTION  * * * *                       
C                                                                       
*CX1  REAL CPU0,CPU,TIM,WALL0,WALL                                      
*CX1  REAL*4 SECOND                                                     
*CX1  LOGICAL FIRST                                                     
*CX1  COMMON /TIMING/ CPU,WALL                                          
*CX1  SAVE FIRST,CPU0,WALL0                                             
*CX1  DATA FIRST/.TRUE./                                                
C                                                                       
C       ----- INITIALIZE CLOCKS -----                                   
C                                                                       
*CX1  IF(FIRST) THEN                                                    
*CX1     FIRST=.FALSE.                                                  
*CX1     CPU0 = 0.0E+00                                                 
*CX1     WALL0 = TIMEF()                                                
*CX1     WALL0 = 0.001E+00*WALL0                                        
*CX1  END IF                                                            
C                                                                       
C       ----- OBTAIN ELAPSED TIMES SINCE JOB STARTED -----              
C                                                                       
*CX1  CPU = SECOND()                                                    
*CX1  WALL = TIMEF()                                                    
*CX1  WALL = 0.001E+00 * WALL - WALL0                                   
*CX1  TIM = CPU                                                         
*CX1  RETURN                                                            
C                                                                       
C-----------------------------------------------------------------------
C      * * * *  DECSTATION SPECIFIC PORTION  * * * *                    
C                                                                       
*DEC  IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C        --EXTERNAL STATEMENT NEEDED HERE                               
*DEC  EXTERNAL ETIME,TIME                                               
*DEC  DOUBLE PRECISION ETIME                                            
*DEC  REAL TARRAY(2)                                                    
*DEC  INTEGER TIME                                                      
*DEC  LOGICAL FIRST                                                     
*DEC  COMMON /TIMING/ CPU,WALL                                          
*DEC  SAVE FIRST,CPU0,WALL0                                             
*DEC  DATA FIRST/.TRUE./                                                
C                                                                       
C       ----- INITIALIZE CLOCKS -----                                   
C                                                                       
*DEC  IF(FIRST) THEN                                                    
*DEC     FIRST=.FALSE.                                                  
*DEC     DUMMY = ETIME(TARRAY)                                          
*DEC     CPU0 = TARRAY(1)+TARRAY(2)                                     
*DEC     WALL0 = TIME()                                                 
*DEC  END IF                                                            
C                                                                       
C       ----- OBTAIN ELAPSED TIMES SINCE JOB STARTED -----              
C                                                                       
*DEC  DUMMY = ETIME(TARRAY)                                             
*DEC  CPU = TARRAY(1)+TARRAY(2)                                         
*DEC  CPU = CPU - CPU0                                                  
*DEC  WALL = TIME()                                                     
*DEC  WALL = WALL - WALL0                                               
*DEC  TIM = CPU                                                         
*DEC  RETURN                                                            
C                                                                       
C-----------------------------------------------------------------------
C      * * * *  GENERIC FORTRAN 77 SPECIFIC PORTION  * * * *            
C                                                                       
*F77  DOUBLE PRECISION CPU,WALL                                         
*F77  COMMON /TIMING/ CPU,WALL                                          
C                                                                       
*F77  CPU = 0.0D+00                                                     
*F77  WALL =0.0D+00                                                     
*F77  RETURN                                                            
C                                                                       
C-----------------------------------------------------------------------
C      * * * *  FUJITSU SPECIFIC PORTION  * * * *                       
C                                                                       
*FUJ  IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C        --EXTERNAL STATEMENT NEEDED HERE                               
*FUJ  EXTERNAL ETIME,TIME                                               
*FUJ  DOUBLE PRECISION ETIME                                            
*FUJ  REAL TARRAY(2)                                                    
*FUJ  INTEGER TIME                                                      
*FUJ  LOGICAL FIRST                                                     
*FUJ  COMMON /TIMING/ CPU,WALL                                          
*FUJ  SAVE FIRST,CPU0,WALL0                                             
*FUJ  DATA FIRST/.TRUE./                                                
C                                                                       
C       ----- INITIALIZE CLOCKS -----                                   
C                                                                       
*FUJ  IF(FIRST) THEN                                                    
*FUJ     FIRST=.FALSE.                                                  
*FUJ     DUMMY = ETIME(TARRAY)                                          
*FUJ     CPU0 = TARRAY(1)+TARRAY(2)                                     
*FUJ     WALL0 = TIME()                                                 
*FUJ  END IF                                                            
C                                                                       
C       ----- OBTAIN ELAPSED TIMES SINCE JOB STARTED -----              
C                                                                       
*FUJ  DUMMY = ETIME(TARRAY)                                             
*FUJ  CPU = TARRAY(1)+TARRAY(2)                                         
*FUJ  CPU = CPU - CPU0                                                  
*FUJ  WALL = TIME()                                                     
*FUJ  WALL = WALL - WALL0                                               
*FUJ  TIM = CPU                                                         
*FUJ  RETURN                                                            
C                                                                       
C-----------------------------------------------------------------------
C      * * * *  HITACHI RUNNING HI-UX/MPP PORTION  * * * *              
C                                                                       
*HIT  DOUBLE PRECISION CPU,WALL,TIM                                     
*HIT  LOGICAL FIRST                                                     
*HIT  COMMON /TIMING/ CPU,WALL                                          
*HIT  DATA FIRST/.TRUE./                                                
C                                                                       
*HIT  IF(FIRST) THEN                                                    
*HIT     FIRST=.FALSE.                                                  
*HIT     CALL XCLOCK(DUMMY,3)                                           
*HIT     CALL XCLOCK(DUMMY,7)                                           
*HIT  END IF                                                            
C                                                                       
*HIT  CALL XCLOCK(CPU,5)                                                
*HIT  CALL XCLOCK(WALL,8)                                               
*HIT  TIM = CPU                                                         
*HIT  RETURN                                                            
C                                                                       
C-----------------------------------------------------------------------
C      * * * *  HP 9000 SPECIFIC PORTION  * * * *                       
C                                                                       
*HP   DOUBLE PRECISION CPU0,WALL0,CPU,WALL,TIM                          
*HP   EXTERNAL ETIME,TIME                                               
*HP   REAL ETIME                                                        
*HP   REAL TARRAY(2)                                                    
*HP   INTEGER TIME                                                      
*HP   LOGICAL FIRST                                                     
*HP   COMMON /TIMING/ CPU,WALL                                          
*HP   SAVE FIRST,CPU0,WALL0                                             
*HP   DATA FIRST/.TRUE./                                                
*HP   IF(FIRST) THEN                                                    
*HP      FIRST=.FALSE.                                                  
*HP      CPU0 = ETIME(TARRAY)                                           
*HP      WALL0 = TIME()                                                 
*HP   END IF                                                            
*HP   CPU = ETIME(TARRAY)                                               
*HP   CPU = CPU - CPU0                                                  
*HP   WALL = TIME()                                                     
*HP   WALL = WALL - WALL0                                               
*HP   TIM = CPU                                                         
*HP   RETURN                                                            
C                                                                       
C-----------------------------------------------------------------------
C      * * * *  IBM MAINFRAME SPECIFIC PORTION  * * * *                 
C                                                                       
*IBM  DOUBLE PRECISION CPU0,WALL0,CPU,WALL,TIM                          
*IBM  LOGICAL FIRST                                                     
*IBM  COMMON /TIMING/ CPU,WALL                                          
*IBM  SAVE FIRST,CPU0,WALL0                                             
*IBM  DATA FIRST/.TRUE./                                                
C                                                                       
C       ----- INITIALIZE CLOCKS -----                                   
C                                                                       
*IBM  IF(FIRST) THEN                                                    
*IBM     FIRST=.FALSE.                                                  
*IBM     WALL0 =0.0D+00                                                 
*IBM     CPU0 =0.0D+00                                                  
*IBM     CALL ZTIME1(WALL0,CPU0)                                        
*IBM  END IF                                                            
C                                                                       
C       ----- OBTAIN ELAPSED TIMES SINCE JOB STARTED -----              
C                                                                       
*IBM  CPU = CPU0                                                        
*IBM  WALL = WALL0                                                      
*IBM  CALL ZTIME2(WALL,CPU)                                             
*IBM  TIM = CPU                                                         
*IBM  RETURN                                                            
C-----------------------------------------------------------------------
C           * * * * 64 BIT INTEL CHIPS, USING IFORT * * * *             
C         THIS IS AN EXACT COPY OF THE 32 BIT LINUX CODE                
C         SEE THOSE COMMENTS ABOUT THE XQTIME USED ON LINUX...OURS!     
C                                                                       
*INT  DOUBLE PRECISION CPU0,WALL0,CPU,WALL,TIM                          
*INT  DOUBLE PRECISION XQTIME                                           
*INT  REAL TARRAY(2)                                                    
*INT  LOGICAL FIRST                                                     
*INT  COMMON /TIMING/ CPU,WALL                                          
*INT  SAVE FIRST,CPU0,WALL0                                             
*INT  DATA FIRST/.TRUE./                                                
C                                                                       
*INT  IF(FIRST) THEN                                                    
*INT     TARRAY(1)=0.0E+00                                              
*INT     FIRST=.FALSE.                                                  
*INT     WALL0 = XQTIME(TARRAY)                                         
*INT     CPU0 = TARRAY(1)+TARRAY(2)                                     
*INT  END IF                                                            
C                                                                       
*INT  WALL = XQTIME(TARRAY)                                             
*INT  CPU = TARRAY(1)+TARRAY(2)                                         
*INT  CPU = CPU - CPU0                                                  
*INT  WALL = WALL - WALL0                                               
*INT  TIM = CPU                                                         
*INT  RETURN                                                            
C                                                                       
C-----------------------------------------------------------------------
C           * * * *  32 BIT LINUX VERSION * * * *                       
C     THIS USES A CUSTOM VERSION OF -XQTIME-, IMPLEMENTED IN ZUNIX.C.   
C     IN PARTICULAR, THIS CUSTOM XQTIME WILL RETURN THE ELAPSED TIME,   
C     MEANING WALL CLOCK TIME, WHICH DIFFERS FROM THE STANDARD RETURN   
C     VALUE OF THIS ROUTINE WHEN IT EXISTS, WHICH IS A CPU TIME.        
C     THE RETURN VALUE IS DOUBLE PRECISION, BUT THE ARGS ARE NOT.       
C                                                                       
*L32  DOUBLE PRECISION CPU0,WALL0,CPU,WALL,TIM                          
*L32  DOUBLE PRECISION XQTIME                                           
*L32  REAL*4 TARRAY(2)                                                  
*L32  LOGICAL FIRST                                                     
*L32  COMMON /TIMING/ CPU,WALL                                          
*L32  SAVE FIRST,CPU0,WALL0                                             
*L32  DATA FIRST/.TRUE./                                                
C                                                                       
*L32  IF(FIRST) THEN                                                    
*L32     TARRAY(1)=0.0E+00                                              
*L32     FIRST=.FALSE.                                                  
*L32     WALL0 = XQTIME(TARRAY)                                         
*L32     CPU0 = TARRAY(1)+TARRAY(2)                                     
*L32  END IF                                                            
C                                                                       
*L32  WALL = XQTIME(TARRAY)                                             
*L32  CPU = TARRAY(1)+TARRAY(2)                                         
*L32  CPU = CPU - CPU0                                                  
*L32  WALL = WALL - WALL0                                               
*L32  TIM = CPU                                                         
*L32  RETURN                                                            
C-----------------------------------------------------------------------
C           * * * *  64 BIT LINUX VERSION * * * *                       
C         THIS IS AN EXACT COPY OF THE 32 BIT LINUX CODE                
C         SEE THOSE COMMENTS ABOUT THE XQTIME USED ON LINUX...OURS!     
C                                                                       
*L64  DOUBLE PRECISION CPU0,WALL0,CPU,WALL,TIM                          
*L64  DOUBLE PRECISION XQTIME                                           
*L64  REAL*4 TARRAY(2)                                                  
*L64  LOGICAL FIRST                                                     
*L64  COMMON /TIMING/ CPU,WALL                                          
*L64  SAVE FIRST,CPU0,WALL0                                             
*L64  DATA FIRST/.TRUE./                                                
C                                                                       
*L64  IF(FIRST) THEN                                                    
*L64     TARRAY(1)=0.0E+00                                              
*L64     FIRST=.FALSE.                                                  
*L64     WALL0 = XQTIME(TARRAY)                                         
*L64     CPU0 = TARRAY(1)+TARRAY(2)                                     
*L64  END IF                                                            
C                                                                       
*L64  WALL = XQTIME(TARRAY)                                             
*L64  CPU = TARRAY(1)+TARRAY(2)                                         
*L64  CPU = CPU - CPU0                                                  
*L64  WALL = WALL - WALL0                                               
*L64  TIM = CPU                                                         
*L64  RETURN                                                            
C-----------------------------------------------------------------------
C                  * * * NEC SX (SUPER-UX) * * *                        
*NEC  LOGICAL FIRST                                                     
*NEC  COMMON /TIMING/ CPU,WALL                                          
*NEC  SAVE FIRST,CPU0,WALL0                                             
*NEC  DATA FIRST/.TRUE./                                                
*NEC  CHARACTER*8 SXDAT                                                 
C                                                                       
*NEC  IF(FIRST) THEN                                                    
*NEC     FIRST=.FALSE.                                                  
*NEC     ITIC = MCLOCK()                                                
*NEC     CPU0 = ITIC * 0.005                                            
*NEC     CALL DATIM(SXDAT,SXTIM,3)                                      
*NEC     WALL0 = SXTIM*3600.E0                                          
*NEC  END IF                                                            
C                                                                       
*NEC  ITIC = MCLOCK()                                                   
*NEC  CPU  = ITIC * 0.005                                               
*NEC  CPU = CPU - CPU0                                                  
*NEC  CALL DATIM(SXDAT,SXTIM,3)                                         
*NEC  WALL  = SXTIM*3600.E0                                             
*NEC  WALL = WALL - WALL0                                               
*NEC  TIM = CPU                                                         
*NEC  RETURN                                                            
C-----------------------------------------------------------------------
C                                                                       
C                  * * * SILICON GRAPHICS * * *                         
*SGI  IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
*SGI  EXTERNAL TIME                                                     
*SGI  DOUBLE PRECISION ETIME                                            
*SGI  REAL TARRAY(2)                                                    
*SGI  INTEGER TIME                                                      
*SGI  LOGICAL FIRST                                                     
*SGI  COMMON /TIMING/ CPU,WALL                                          
*SGI  SAVE FIRST,CPU0,WALL0                                             
*SGI  DATA FIRST/.TRUE./                                                
C                                                                       
*SGI  IF(FIRST) THEN                                                    
*SGI     FIRST=.FALSE.                                                  
*SGI     DUMMY = ETIME(TARRAY)                                          
*SGI     CPU0 = TARRAY(1)+TARRAY(2)                                     
*SGI     WALL0 = TIME()                                                 
*SGI  END IF                                                            
C                                                                       
*SGI  DUMMY = ETIME(TARRAY)                                             
*SGI  CPU = TARRAY(1)+TARRAY(2)                                         
*SGI  CPU = CPU - CPU0                                                  
*SGI  WALL = TIME()                                                     
*SGI  WALL = WALL - WALL0                                               
*SGI  TIM = CPU                                                         
*SGI  RETURN                                                            
C                                                                       
C-----------------------------------------------------------------------
C                  * * * SUN WORKSTATIONS * * *                         
*SUN  IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
*SUN  REAL ETIME,TARRAY(2)                                              
*SUN  INTEGER TIME                                                      
*SUN  LOGICAL FIRST                                                     
*SUN  COMMON /TIMING/ CPU,WALL                                          
*SUN  SAVE FIRST,CPU0,WALL0                                             
*SUN  DATA FIRST/.TRUE./                                                
C                                                                       
*SUN  IF(FIRST) THEN                                                    
*SUN     FIRST=.FALSE.                                                  
*SUN     DUMMY = ETIME(TARRAY)                                          
*SUN     CPU0 = TARRAY(1)+TARRAY(2)                                     
*SUN     WALL0 = TIME()                                                 
*SUN  END IF                                                            
C                                                                       
*SUN  DUMMY = ETIME(TARRAY)                                             
*SUN  CPU = TARRAY(1)+TARRAY(2)                                         
*SUN  CPU = CPU - CPU0                                                  
*SUN  WALL = TIME()                                                     
*SUN  WALL = WALL - WALL0                                               
*SUN  TIM = CPU                                                         
*SUN  RETURN                                                            
C                                                                       
C-----------------------------------------------------------------------
C      * * * *  CRAY T3E SPECIFIC PORTION  (IDENTICAL TO *CRY) * * * *  
C                                                                       
*T3E  REAL CPU0,WALL0,CPU,WALL,TIM                                      
*T3E  LOGICAL FIRST                                                     
*T3E  COMMON /TIMING/ CPU,WALL                                          
*T3E  SAVE FIRST,CPU0,WALL0                                             
*T3E  DATA FIRST/.TRUE./                                                
C                                                                       
C       ----- INITIALIZE CLOCKS -----                                   
C                                                                       
*T3E  IF(FIRST) THEN                                                    
*T3E     FIRST=.FALSE.                                                  
*T3E     CPU0 = 0.0E+00                                                 
*T3E     WALL0 = TIMEF()                                                
*T3E     WALL0 = 0.001E+00*WALL0                                        
*T3E  END IF                                                            
C                                                                       
C       ----- OBTAIN ELAPSED TIMES SINCE JOB STARTED -----              
C                                                                       
*T3E  CALL SECOND(CPU)                                                  
*T3E  WALL = TIMEF()                                                    
*T3E  WALL = 0.001E+00 * WALL - WALL0                                   
*T3E  TIM = CPU                                                         
*T3E  RETURN                                                            
C                                                                       
C-----------------------------------------------------------------------
C      * * * *  VMS SPECIFIC PORTION  * * * *                           
C                                                                       
*VMS  DOUBLE PRECISION CPU0,WALL0,CPU,WALL,TIM                          
*VMS  CHARACTER*8 BUF                                                   
*VMS  LOGICAL FIRST,GOPARR,DSKWRK,MASWRK                                
*VMS  INTEGER CPUTIM, ITMLST(4), SYS$GETJPI                             
*VMS  PARAMETER JPI$_CPUTIM = '04070004'X   ! ACCUMULATED CPU, 10 MS    
*VMS  COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
*VMS  COMMON /TIMING/ CPU,WALL                                          
*VMS  SAVE FIRST,CPU0,WALL0                                             
*VMS  DATA FIRST/.TRUE./                                                
C                                                                       
C       ----- INITIALIZE CLOCKS -----                                   
C                                                                       
*VMS  IF(FIRST) THEN                                                    
*VMS     FIRST=.FALSE.                                                  
C              INITIALIZE ITEM LIST FOR VMS GETJPI                      
*VMS     ITMLST(1) = JPI$_CPUTIM  ! SET ITEM CODE, INCLUDING LENGTH     
*VMS     ITMLST(2) = %LOC(CPUTIM) ! LOCATION OF CPUTIM VARIABLE         
*VMS     ITMLST(3) = 0            ! LENGTH NOT NEEDED                   
*VMS     ITMLST(4) = 0            ! TERMINATE ITMLST                    
C              GET CPUTIM INFORMATION FROM GETJPI                       
*VMS     IF( .NOT. SYS$GETJPI(,,,ITMLST,,,) .AND. MASWRK)               
*VMS *      WRITE(6,*) 'GETJPI ERROR'                                   
*VMS     CPU0=CPUTIM*0.01D+00                                           
C              NOW GET STARTING WALL TIME                               
*VMS     CALL IDATE(MONTH0,IDAY,IYEAR)                                  
*VMS     CALL TIME(BUF)                                                 
*VMS     READ(UNIT=BUF,FMT='(I2,1X,I2,1X,I2)') IHH,IMM,ISS              
*VMS     WALL0 = ((IDAY*24.0+IHH)*60.0+IMM)*60.0+ISS                    
*VMS  END IF                                                            
C                                                                       
C       ----- OBTAIN ELAPSED TIMES SINCE JOB STARTED -----              
C   (INCORRECT OVER MONTH END, IF YEAR END OR IF PREV. MONTH.NE.31 DAYS)
C                                                                       
*VMS  IF( .NOT. SYS$GETJPI(,,,ITMLST,,,) .AND. MASWRK)                  
*VMS *   WRITE(6,*) 'GETJPI ERROR'                                      
*VMS  CPU=CPUTIM*0.01D+00 - CPU0                                        
*VMS  CALL IDATE(MONTH,IDAY,IYEAR)                                      
*VMS  IF(MONTH.GT.MONTH0) IDAY=IDAY+31   ! WELL, SOME MONTHS ANYWAY     
*VMS  CALL TIME(BUF)                                                    
*VMS  READ(UNIT=BUF,FMT='(I2,1X,I2,1X,I2)') IHH,IMM,ISS                 
*VMS  WALL = ((IDAY*24.0+IHH)*60.0+IMM)*60.0+ISS - WALL0                
*VMS  TIM=CPU                                                           
*VMS  RETURN                                                            
C                                                                       
C-----------------------------------------------------------------------
C           * * * *  MS WINDOWS AZURE SPECIFIC PORTION  * * * *         
C                                                                       
*WAZ  DOUBLE PRECISION CPU,WALL,WALLIN, TIM, CPUIN, CPUOUT              
*WAZ  INTEGER ICPUIN, ICPUOUT, IWALLIN, IWALLOUT, ITICKMAX, IFREQ       
*WAZ  LOGICAL FIRST                                                     
*WAZ  COMMON /TIMING/ CPU, WALL                                         
*WAZ  SAVE FIRST,ICPUIN,IWALLIN,CPUIN                                   
*WAZ  DATA FIRST/.TRUE./                                                
C                                                                       
*WAZ  IF(FIRST) THEN                                                    
*WAZ     FIRST=.FALSE.                                                  
*WAZ     CALL SYSTEM_CLOCK(IWALLIN,IFREQ,ITICKMAX)                      
*WAZ     !CALL SYSTEM_CLOCK(ICPUIN,IFREQ,ITICKMAX)                      
*WAZ     CALL CPU_TIME(CPUIN)                                           
*WAZ  END IF                                                            
C                                                                       
*WAZ  CALL SYSTEM_CLOCK(IWALLOUT,IFREQ,ITICKMAX)                        
*WAZ  !CALL SYSTEM_CLOCK(ICPUOUT,IFREQ,ITICKMAX)                        
*WAZ  CALL CPU_TIME(CPUOUT)                                             
*WAZ  !CPU=DBLE(ICPUOUT-ICPUIN)/DBLE(IFREQ)                             
*WAZ  CPU=CPUOUT-CPUIN                                                  
*WAZ  WALL=DBLE(IWALLOUT-IWALLIN)/DBLE(IFREQ)                           
*WAZ  TIM=CPU                                                           
*WAZ  RETURN                                                            
C-----------------------------------------------------------------------
C           * * * *  MS WINDOWS (WIN32) SPECIFIC PORTION  * * * *       
C                                                                       
*W32  IMPLICIT NONE                                                     
*W32  EXTERNAL WINDOWSTIMER                                             
*W32  DOUBLE PRECISION FUNCTION WINDOWSTIMER                            
*W32  DOUBLE PRECISION CPU0,WALL0,CPU,WALL,TIM,WINDOWSTIMER             
*W32  LOGICAL FIRST                                                     
*W32  COMMON /TIMING/ CPU,WALL                                          
*W32  SAVE FIRST,CPU0,WALL0                                             
*W32  DATA FIRST/.TRUE./                                                
*W32  IF(FIRST) THEN                                                    
*W32     FIRST = .FALSE.                                                
*W32     CPU0  = WINDOWSTIMER(WALL0)                                    
*W32  END IF                                                            
*W32  CPU  = WINDOWSTIMER(WALL)                                         
*W32  CPU  = CPU - CPU0                                                 
*W32  WALL = WALL - WALL0                                               
*W32  TIM  = CPU                                                        
*W32  RETURN                                                            
C                                                                       
C-----------------------------------------------------------------------
C           * * * *  MS WINDOWS (WIN64) SPECIFIC PORTION  * * * *       
C                                                                       
*64W  IMPLICIT NONE                                                     
*64W  EXTERNAL WINDOWSTIMER                                             
*64W  DOUBLE PRECISION FUNCTION WINDOWSTIMER                            
*64W  DOUBLE PRECISION CPU0,WALL0,CPU,WALL,TIM,WINDOWSTIMER             
*64W  LOGICAL FIRST                                                     
*64W  COMMON /TIMING/ CPU,WALL                                          
*64W  SAVE FIRST,CPU0,WALL0                                             
*64W  DATA FIRST/.TRUE./                                                
*64W  IF(FIRST) THEN                                                    
*64W     FIRST = .FALSE.                                                
*64W     CPU0  = WINDOWSTIMER(WALL0)                                    
*64W  END IF                                                            
*64W  CPU  = WINDOWSTIMER(WALL)                                         
*64W  CPU  = CPU - CPU0                                                 
*64W  WALL = WALL - WALL0                                               
*64W  TIM  = CPU                                                        
*64W  RETURN                                                            
C                                                                       
C-----------------------------------------------------------------------
C     * * * *  MS WINDOWS & 64 BIT INTEL CHIPS, USING IFORT  * * * *    
C     * * * *                 SPECIFIC PORTION               * * * *    
C                                                                       
*64I  EXTERNAL WALLTIMER                                                
*64I  DOUBLE PRECISION FUNCTION WALLTIMER                               
*64I  DOUBLE PRECISION CPU0,WALL0,CPU,WALL,TIM,WALLTIMER                
*64I  EXTERNAL ETIME                                                    
*64I  DOUBLE PRECISION ETIME                                            
*64I  REAL TARRAY(2)                                                    
*64I  LOGICAL FIRST                                                     
*64I  COMMON /TIMING/ CPU,WALL                                          
*64I  SAVE FIRST,CPU0,WALL0                                             
*64I  DATA FIRST/.TRUE./                                                
C                                                                       
*64I  IF(FIRST) THEN                                                    
*64I     TARRAY(1)=0.0E+00                                              
*64I     FIRST=.FALSE.                                                  
*64I     WALL0 = ETIME(TARRAY)                                          
*64I     WALL0 = WALLTIMER()                                            
*64I     CPU0 = TARRAY(1)+TARRAY(2)                                     
*64I  END IF                                                            
C                                                                       
*64I  WALL = ETIME(TARRAY)                                              
*64I  WALL = WALLTIMER()                                                
*64I  CPU = TARRAY(1)+TARRAY(2)                                         
*64I  CPU = CPU - CPU0                                                  
*64I  WALL = WALL - WALL0                                               
*64I  TIM = CPU                                                         
*64I  RETURN                                                            
C                                                                       
C-----------------------------------------------------------------------
C           * * * *  CRAY XT3 SPECIFIC PORTION  * * * *                 
C                                                                       
*XT3  DOUBLE PRECISION CPU0,WALL0,CPU,WALL,TIM,MPI_WTIME                
*XT3  LOGICAL FIRST                                                     
*XT3  INTEGER TIME                                                      
*XT3  COMMON /TIMING/ CPU,WALL                                          
*XT3  SAVE FIRST,CPU0,WALL0                                             
*XT3  DATA FIRST/.TRUE./                                                
C                                                                       
*XT3  IF(FIRST) THEN                                                    
*XT3     FIRST=.FALSE.                                                  
*XT3     CALL CPU_TIME(CPU0)                                            
*XT3     WALL0 = MPI_WTIME()                                            
*XT3  END IF                                                            
C                                                                       
*XT3  CALL CPU_TIME(CPU)                                                
*XT3  CPU = CPU - CPU0                                                  
*XT3  WALL = MPI_WTIME()                                                
*XT3  WALL = WALL - WALL0                                               
*XT3  TIM = CPU                                                         
*XT3  RETURN                                                            
      END                                                               
C                                                                       
C*MODULE UNPORT  *DECK VNAN                                             
      SUBROUTINE VNAN(ARRAY,ISTRIDE,LENGTH)                             
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      DIMENSION ARRAY(*)                                                
C                                                                       
      DATA QNAN/Z'7FF8000000000001'/                                    
*CRY  DATA QNAN/9.87654321D+30/                                         
*CX1  DATA QNAN/9.87654321D+30/                                         
*DEC  DATA QNAN/Z'7FF8000000000001'/                                    
*F77  DATA QNAN/9.87654321D+30/                                         
*FUJ  DATA QNAN/9.87654321D+30/                                         
*HIT  DATA QNAN/9.87654321D+30/                                         
*HP   DATA QNAN/Z'7FF8000000000001'/                                    
*IBM  DATA QNAN/9.87654321D+30/                                         
*INT  DATA QNAN/Z'7FF8000000000001'/                                    
*L32  INTEGER*4 NAN(2)                                                  
*L32  EQUIVALENCE (QNAN,NAN(1))                                         
*L32  DATA  NAN/Z'7FF80001',Z'7FF80001'/  ! FOR BIG/LITTLE ENDIAN       
*L64  INTEGER NAN                                                       
*L64  EQUIVALENCE (QNAN,NAN)                                            
*L64  DATA  NAN/Z'7FF8000000000001'/                                    
*NEC  DATA QNAN/9.87654321D+30/                                         
*SGI  DATA QNAN/Z'7FF8000000000001'/                                    
*SUN  DATA QNAN/Z'7FF8000000000001'/                                    
*T3E  DATA QNAN/9.87654321D+30/                                         
*VMS  DATA QNAN/9.87654321D+30/                                         
*XT3  DATA QNAN/9.87654321D+30/                                         
*W32  INTEGER*4 NAN(2)                                                  
*W32  EQUIVALENCE (QNAN,NAN(1))                                         
*W32  DATA  NAN/Z'7FF80001',Z'7FF80001'/  ! FOR BIG/LITTLE ENDIAN       
*W64  INTEGER NAN                                                       
*W64  EQUIVALENCE (QNAN,NAN)                                            
*W64  DATA  NAN/Z'7FF8000000000001'/                                    
*I64  DATA QNAN/Z'7FF8000000000001'/                                    
C                                                                       
C     ---- INITIALIZES -ARRAY- TO QUIET FORM OF "NOT A NUMBER" ----     
C     SET -LENGTH- TOTAL ELEMENTS, SEPARATED BY INCREMENT -ISTRIDE-.    
C     THIS IS MEANT FOR DEBUGGING PURPOSES, SINCE IT IS THE ANTI-VCLR.  
C                                                                       
C     WE CHOOSE THE SO-CALLED QUIET FORM OF NOT-A-NUMBER, WHOSE         
C     IEEE REPRESENTATION REQUIRES SETTING ALL EXPONENT BITS, WITH      
C     A NON-ZERO VALUE IN THE MANTISSA.  A 'QUIET NAN', AS OPPOSED      
C     TO THE SIGNALLING-NOT-A-NUMBER (Z'7FF0000000000001'), HAS THE     
C     MANTISSA'S SIGN BIT SET.  NOTE: INFINITY=Z'7FF0000000000000'.     
C                                                                       
C     USAGE IN THE MIDDLE OF A DODGY BIT OF CODE:                       
C               CALL VALFM(LOADFM)                                      
C               LA   = LOADFM + 1                                       
C               LB   = LA     + ...                                     
C                   ...                                                 
C               LQ   = LP     + ...                                     
C               LAST = LQ     + ...                                     
C               NEED = LAST - LOADFM - 1                                
C               CALL GETFM(NEED)                                        
C         C------------SET BAD VALUES--------                           
C               CALL VNAN(XX(LOADFM+1),1,NEED)                          
C         C------------SET BAD VALUES--------                           
C               ...DO SOMETHING WITH THE MEMORY HERE...                 
C               CALL RETFM(NEED)                                        
C                                                                       
      IF (ISTRIDE.EQ.1) THEN                                            
         DO L=1,LENGTH                                                  
            ARRAY(L) = QNAN                                             
         ENDDO                                                          
      ELSE                                                              
         LA=1-ISTRIDE                                                   
         DO L=1,LENGTH                                                  
            LA=LA+ISTRIDE                                               
            ARRAY(LA) = QNAN                                            
         ENDDO                                                          
      END IF                                                            
      RETURN                                                            
C                                                                       
C-----------------------------------------------------------------------
C     THE END STATEMENT FOR THE ROUTINE ABOVE IS AT THE VERY END OF     
C     THIS FILE, SO THAT WE CAN OBEY THE FORTRAN RULE THAT THERE        
C     SHALL BE NO COMMENTS AFTER THE LAST END STATEMENT.                
C-----------------------------------------------------------------------
C                                                                       
C*MODULE UNPORT  *DECK ERSATZ                                           
C                                                                       
C     ANY LITTLE DUMMY ROUTINES THAT ONLY ONE SYSTEM NEEDS GO HERE,     
C     EACH MUST BEGIN WITH AN END STATEMENT FOR THE ABOVE ROUTINE.      
C                                                                       
C        SIMULATE FORTRAN90 INTRINSIC, BUT FOR ONLY 1 RANDOM NO.        
C                                                                       
*DEC  END                                                               
*DEC  SUBROUTINE RANDOM_NUMBER(XX)                                      
*DEC  DOUBLE PRECISION XX                                               
*DEC  LOGICAL FIRST                                                     
*DEC  DATA FIRST/.TRUE./                                                
*DEC  SAVE FIRST                                                        
*DEC  IRTP=3                                                            
*DEC  IF(FIRST) THEN                                                    
*DEC     FIRST=.FALSE.                                                  
*DEC     CALL RNGEN(XX,1,IRTP,0)                                        
*DEC  ELSE                                                              
*DEC     CALL RNGEN(XX,1,IRTP,1)                                        
*DEC  ENDIF                                                             
*DEC  RETURN                                                            
C                                                                       
C        GETENV USING A POSIX CALL, USEFUL ON THE CRAY T3E              
C                                                                       
*T3E  END                                                               
*T3E  SUBROUTINE GETENV(NAME,VALUE)                                     
*T3E  CHARACTER*(*) NAME,VALUE                                          
*T3E  CALL PXFGETENV(NAME,LENNAM,VALUE,LENVAL,IERROR)                   
*T3E  RETURN                                                            
C                                                                       
C        GETENV USING A POSIX CALL, USEFUL ON THE CRAY X1               
C                                                                       
*CX1  END                                                               
*CX1  SUBROUTINE GETENV(NAME,VALUE)                                     
*CX1  CHARACTER*(*) NAME,VALUE                                          
*CX1  CALL PXFGETENV(NAME,LENNAM,VALUE,LENVAL,IERROR)                   
*CX1  RETURN                                                            
C                                                                       
C        VB2000 always tries to use MEMGET/MEMREL                       
*HP   END                                                               
*HP   INTEGER FUNCTION MEMGET(MEMLIM)                                   
*HP   INTEGER MEMLIM,NBYTES                                             
*HP   NBYTES = (MEMLIM+2)*8                                             
*HP   MEMGET = MALLOC(NBYTES)                                           
*HP   RETURN                                                            
*HP   END                                                               
*HP   SUBROUTINE MEMREL(LOCMEM)                                         
*HP   INTEGER LOCMEM                                                    
*HP   CALL FREE(LOCMEM)                                                 
*HP   RETURN                                                            
C                                                                       
C        THE COMMENT LETTERS FGE STAND FOR FILE/GET/ENVIRONMENT         
C                                                                       
C        GETENV IS REPLACED BY READING VARIABLE=VALUE PAIRS FROM DISK,  
C        WHICH MAY BE USEFUL IF THE COMPUTER IS BRAIN DEAD ABOUT GETENV.
C                                                                       
*FGE  END                                                               
*FGE  SUBROUTINE FILE_GETENV(NAME, VALUE)                               
*FGE  CHARACTER*(*) NAME, VALUE                                         
*FGE  CHARACTER*6   ENVFIL                                              
*FGE  CHARACTER*256 FILENM                                              
C        264=256 PLUS 7 LETTER ENV.VAR.NAME PLUS A SEPARATOR            
*FGE  CHARACTER*264 BUFFER                                              
*FGE  LOGICAL GOPARR,DSKWRK,MASWRK                                      
*FGE  COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
C                                                                       
C   THIS IS INTENDED FOR IBM BLUE GENE, WHERE THE ENVIRONMENTAL         
C   VARIABLES HAVE TO TAKE LESS THAN 1KB OF MEMORY, OR ON OTHER         
C   SYSTEMS THAT FAIL TO BOTHER TO PASS THE ENVIRONMENT PROPERLY.       
C                                                                       
C   SEE THE "FGE" HACK IN THE COMP SCRIPT TO USE THIS ROUTINE.          
C                                                                       
C   THIS ROUTINE TRIES TO GET ENVIRONMENT VALUES BY                     
C   1. ASK FOR THE VALUE FOR -NAME- FROM THE ENVIRONMENT,               
C   2. ASK FOR THE VALUE FOR -NAME- FROM DISK, FROM A FILENAME          
C          A) GIVEN BY ENVIRONMENT VARIABLE -ENVFIL-                    
C          B) ASSUMED TO BE GAMESS.ENVIRON                              
C   IN THE LATTER CASE, EACH USER WILL BE LIMITED TO ONE JOB AT         
C   A TIME ON THE MASTER NODE, AND THE RUNTIME SCRIPT SHOULD            
C   BE SURE THE WORKING DIRECTORY IS THE SCRATCH DIRECTORY, WHERE       
C   THE 'GAMESS.ENVIRON' FILE HAS BEEN CREATED.  SHIKATA GA NAI.        
C                                                                       
C   IF THE VALUES COME FROM A DISK FILE, ONLY THE MASTER NODE           
C   SHOULD PROCESS THE FILE.                                            
C   IT IS THE RESPONSIBILITY OF THE CALLER TO BE SURE THAT ONLY         
C   THE MASTER NODE CALLS THIS ROUTINE, FOR THIS ROUTINE NO LONGER      
C   CONTAINS A RETURN (ANDREY DOESN'T LIKE THIS FOR SOME REASON).       
C                                                                       
C   IT IS THE RESPONSIBILITY OF THE CALLER TO BROADCAST THE RESULT      
C   TO OTHER NODES, IF THAT IS REQUIRED.  OF COURSE, THAT IS ALWAYS     
C   DONE BY THE ROUTINES THAT CALL GMS_GETENV.                          
C                                                                       
C   THE VALUES IN THE FILE SHOULD BE ONE PER LINE,                      
C   VARIABLE=VALUE, WITH THE SEPARATOR ALLOWED TO BE AN EQUALS,         
C   A SPACE, OR A NULL BYTE.  NOTE THAT THE C-SHELL COMMAND             
C        "ENV >& /SCR/MIKE/GAMESS.ENVIRON"                              
C   PROVIDES SUCH AN EQUALS-SEPARATED FILE VERY EASILY.                 
C                                                                       
C-----*FGE  IF(.NOT.MASWRK) RETURN                                      
C                                                                       
*FGE  IUNIT=2                                                           
*FGE  ENVFIL='ENVFIL'                                                   
C                                                                       
C     TRY INTRINSIC GETENV FIRST, DIRECTLY ON THE VARIABLE WE SEEK.     
C                                                                       
*FGE  CALL GETENV(NAME, VALUE)                                          
*FGE  IF ((VALUE(1:1).NE.' ').AND.(VALUE(1:1).NE.CHAR(0))) RETURN       
C                                                                       
C     OTHERWISE, WE NEED TO READ THE VALUE FROM A DISK FILE.            
C                                                                       
C     IF THE FILE NAME IS IN THE ENVIRONMENT, USE THAT FILE NAME,       
C     BUT IF NOT, THEN WE MUST HARD-CODE A CONSTANT FILE NAME.          
C                                                                       
*FGE  CALL GETENV(ENVFIL, FILENM)                                       
*FGE  IF ((FILENM(1:1).EQ.' ').OR.(FILENM(1:1).EQ.CHAR(0))) THEN        
*FGE     FILENM(1:14)='GAMESS.ENVIRON'                                  
*FGE  END IF                                                            
C                                                                       
*FGE  OPEN(UNIT=IUNIT, FILE=FILENM, STATUS='OLD',                       
*FGE *     ACCESS='SEQUENTIAL', FORM='FORMATTED', ERR=7)                
*FGE  REWIND(UNIT=IUNIT)                                                
C                                                                       
*FGE1 CONTINUE                                                          
*FGE  READ(UNIT=IUNIT, FMT='(A)', ERR=9, END=9) BUFFER                  
C                                                                       
C        SKIP OVER LEADING SPACES IN THIS LINE                          
*FGE  DO I=1,LEN(BUFFER),1                                              
*FGE     IF (BUFFER(I:I).NE.' ') GO TO 2                                
*FGE  ENDDO                                                             
*FGE  GO TO 1                                                           
C                                                                       
C        FIND THE VALID SEPARATORS BETWEEN VARIABLE AND VALUE           
*FGE2 CONTINUE                                                          
*FGE  DO J=I,LEN(BUFFER),1                                              
*FGE     IF ((BUFFER(J:J).EQ.' ').OR.                                   
*FGE *       (BUFFER(J:J).EQ.'=').OR.                                   
*FGE *       (BUFFER(J:J).EQ.CHAR(0))) GO TO 3                          
*FGE  ENDDO                                                             
*FGE  GO TO 1                                                           
C                                                                       
C        IS THE VARIABLE NAME ON THIS LINE THE ONE WE ARE LOOKING FOR?  
C        IF SO, RETURN THE VALUE, OTHERWISE KEEP LOOKING.               
*FGE3 CONTINUE                                                          
*FGE  IF (NAME.EQ.BUFFER(I:J-1)) THEN                                   
*FGE     DO K=J,LEN(BUFFER),1                                           
*FGE        IF ((BUFFER(K:K).NE.' ').AND.                               
*FGE *          (BUFFER(K:K).NE.'=').AND.                               
*FGE *          (BUFFER(K:K).NE.CHAR(0))) GO TO 4                       
*FGE     ENDDO                                                          
*FGE4    CONTINUE                                                       
*FGE     VALUE = BUFFER(K:)                                             
*FGE     CLOSE(UNIT=IUNIT, STATUS='KEEP')                               
*FGE     RETURN                                                         
*FGE  ELSE                                                              
*FGE     GO TO 1                                                        
*FGE  END IF                                                            
C                                                                       
C     HANDLE MISSING VARIABLE=VALUE FILE BY OFFICIALLY CRASHING         
C                                                                       
*FGE7 CONTINUE                                                          
*FGE  WRITE(6,8) FILENM                                                 
*FGE  CALL ABRT                                                         
*FGE  STOP                                                              
*FGE8 FORMAT(1X,'FAILED TO OPEN DISK FILE ',A,','/                      
*FGE *    1X,'EXPECTED TO CONTAIN THE ENVIRONMENTAL VARIABLES.'/        
*FGE *    1X,'PLEASE SET -ENVFIL- TO YOUR FILE OF ENVIRONMENT STRINGS'/ 
*FGE *    1X,'OR USE THE DISK FILE NAME -GAMESS.ENVIRON-')              
C                                                                       
C     HANDLE MISSING VARIABLES BY RETURNING NOTHING                     
C                                                                       
*FGE9 CONTINUE                                                          
*FGE  VALUE=' '                                                         
*FGE  CLOSE(UNIT=IUNIT, STATUS='KEEP')                                  
*FGE  RETURN                                                            
C                                                                       
      END                                                               
C                                                                       
C*MODULE UNPORT  *DECK SETMXSEQMTX                                      
      SUBROUTINE SETMXSEQMTX(MXSEQ2IN,MXSEQ3IN)                         
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      COMMON /MTXSIZ/ MXSEQ2,MXSEQ3                                     
C                                                                       
C        VERY SMALL MATRIX SIZES SHOULD NOT ATTEMPT TO RUN THE          
C        LINEAR ALGEBRA STEPS IN PARALLEL.  THESE THREE ROUTINES        
C        RETURN SOME NOTION OF WHAT "VERY SMALL" MIGHT MEAN.            
C        IN PRINCIPLE, ONE COULD USE DIFFERENT VALUES FOR               
C        DIFFERENT TYPES OF MACHINES USING *XXX LINES BELOW:            
C                                                                       
      IF(MXSEQ2IN.GT.0) THEN                                            
         MXSEQ2 = MXSEQ2IN                                              
      ELSE                                                              
         MXSEQ2 = 300                                                   
      END IF                                                            
C                                                                       
      IF(MXSEQ3IN.GT.0) THEN                                            
         MXSEQ3 = MXSEQ3IN                                              
      ELSE                                                              
         MXSEQ3 = 150                                                   
      END IF                                                            
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE UNPORT  *DECK MXSQN2                                           
      INTEGER FUNCTION MXSQN2()                                         
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      COMMON /MTXSIZ/ MXSEQ2,MXSEQ3                                     
C       RETURNS MATRIX THRESHOLD SIZE TO FORCE SEQUENTIAL N**2 LOOPS.   
      MXSQN2 = MXSEQ2                                                   
      RETURN                                                            
      END                                                               
C*MODULE UNPORT  *DECK MXSQN3                                           
      INTEGER FUNCTION MXSQN3()                                         
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      COMMON /MTXSIZ/ MXSEQ2,MXSEQ3                                     
C       RETURNS MATRIX THRESHOLD SIZE TO FORCE SEQUENTIAL N**3 LOOPS.   
      MXSQN3 = MXSEQ3                                                   
      RETURN                                                            
      END                                                               
