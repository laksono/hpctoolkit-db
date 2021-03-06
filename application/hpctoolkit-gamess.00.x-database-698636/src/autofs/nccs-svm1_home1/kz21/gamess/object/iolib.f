c 12 Sep 18 - JLGV - Reimplement Jerry's workaroud for CC file sizes    
C  6 Jun 18 - DGF - tweaks for FMO 5.3                                  
c 18 Apr 16 - DGF - changes for FMO 5.2                                 
c 22 Sep 15 - HL  - ALLOW QUANPOL CHANGE FILE SIZES                     
C 23 Feb 15 - SRP - ADDNANODE: NPUniv added so that calculations using M
C                   can use more than 999 ranks.                        
c 22 Oct 14 - DGF - fileless F1O and F15                                
c 12 Aug 13 - DGF - finish FMO 5.0                                      
C 21 May 13 - DGF - tweaks for SMP file distribution over local disks   
C 19 Oct 12 - MWS - synchronize FRGINF common                           
C 21 Aug 12 - MWS - DAREAD,DAWRIT: I/O adjustments for NBO interface    
C 27 Jul 12 - JAB - CCRED2: use different READ when ioerr=0             
C 24 JUL 12 - DGF - PAD COMMON FOR FMO 4.3                              
C 21 JUN 12 - DGF - EXTCAB PATCH TO OPEN IT IN GDDI                     
C 23 MAR 12 - DGF - handle files for GDDI in disguise, pad commons      
C 28 DEC 11 - DGF,HU,SRP - GMS_GETENV change, avoid closing DA in GDDI, 
C 29 JUN 11 - AG  - use MeUniv rather than MeGlob when creating files (M
C 15 APR 11 - MWS - SAVE ENVIRONMENT VALUE FOR BASPATH                  
C 11 AUG 10 - DGF - ADD MEMORY-BASED I/O FOR FMO                        
C 11 AUG 10 - MWS - ADD DIRECT ACCESS OPEN AND CLOSE ROUTINES           
C 23 JUN 10 - MWS - PRE-STORE CSFSAVE NAME, IMPROVE ADDNANODE ERR MSG   
C 14 OCT 09 - DGF - SYNCHRONISE FMOOPT AND FMORUN                       
C  1 MAY 09 - PFS,HL - ALLOW SIZE-CHANGE IN DAWRIT FOR LMOEDA           
C 12 JAN 09 - DGF - I/O REDUCTION FOR FMO                               
C 15 DEC 08 - DGF,MWS - ALLOW IOSMP FOR NON-FMO, RESTART/TRAJECT FILES  
C 20 NOV 08 - MK  - ADD RAOPDC,RAREDC,RAWRDC FOR DIVIDE-AND-CONQUER     
C 11 APR 08 - JMM - GMS_GETENV: MAKEFP OUTPUT GOES TO -MAKEFP- ENV.VAR. 
C 20 AUG 07 - DGF - PAD FILES IN STORENV                                
C 20 AUG 07 - KI  - ADD MPOPEN AND MPCLOS FOR MP2 DIRECT ACCESS FILES   
C 20 AUG 07 - MWS - CCOPEN: PRINT WARNING IF RECORD LENGTH EXCEEDS 2GB  
C 25 JUN 07 - MWS - PKWRT: FIX ERROR MSG'S PRINTING OF THE I/O DIRECTION
C 24 MAR 07 - MWS - PAD FRGINF COMMON BLOCK                             
C 15 FEB 07 - MWS - ADD GMS_GETENV WRAPPER AROUND CALLS TO UNIX GETENV  
C 22 DEC 06 - DGF - ADD SMP I/O OPTION, PAD MORE FILE NAMES TO STORENV  
C 22 DEC 06 - BN  - DAWRIT: INITGDDI MAKES ALL GROUP MASTERS WRITE DAF  
C  6 NOV 06 - MWS - ADJUST GDDI COMMON BLOCK                            
C 13 MAR 06 - DGF - SEQOPN: MOVE NODEXT INITIALIZATION HERE             
C 22 FEB 06 - TN  - STORENV: STORE ONE MORE FILE NAME FOR EFP+FMO       
C 17 JAN 06 - RAK - SEQOPN: TRY TO SET FLUSH-ON-WRITE FOR "PUNCH" FILES 
C 18 NOV 05 - MWS - ADD SQRDCH/SQWRCH CHUNKED SEQUENTIAL I/O ROUTINES   
C 14 NOV 05 - DGF - STORENV: MORE GDDI FILE NAMES BROADCAST             
C 19 SEP 05 - IA  - SYNCHRONIZE FRGINF COMMON                           
C  6 JUL 05 - DGF - STORENV: ADD NEW CC FILE NAMES                      
C  5 JUL 05 - MWS - SELECT NEW ATOM,BASIS,EFP,PCM,DAF DIMENSIONS        
C 30 APR 05 - DGF - CHANGE GDDI FILE EXTENSIONS TO GLOBAL RANK          
C  7 MAR 05 - IA  - FIX COMMON BLOCK FRGINF                             
C 13 FEB 05 - MWS - PAD COMMON BLOCK FRGINF                             
C  5 FEB 05 - DGF - STORENV: SAVE ERIC FILE, PARENV: ALLOW P>999        
C 23 JUL 04 - MWS - PREAD: MODIFY NINTIC WARNING                        
C 19 MAY 04 - DGF - CHANGES TO ADD THE FMO METHOD                       
C  9 DEC 03 - OQ  - RAISE MXUNIT FROM 100 TO 299 (THIS IS BAD FOR F2C)  
C  4 NOV 03 - DGF - DAREAD,DAWRIT: ALLOW EACH NODE TO HAVE ITS OWN DAF  
C  3 SEP 03 - MWS - RAOPEN,RACLOS,CCOPEN,SEQOPN,SEQCLO: MAX. FILE NUMBER
C  3 JUL 03 - JMM - SUPPRESS PRINTING FOR MONTE CARLO JOBS              
C 12 DEC 02 - MWS - DAREAD: IMPROVED ERROR HANDLING/MESSAGES            
C 17 APR 02 - MWS - SYNCH UP FRGINF COMMON                              
C 26 MAR 02 - MWS - ADD ROUTINE TO OPEN CC DIRECT ACCESS FILES          
C  6 SEP 01 - HU  - PREADP: CORRECT LABEL READING LENGTHS               
C 20 FEB 01 - PND - PUVIB: ADITIONAL PRINT OUT CONTROL                  
C 21 DEC 99 - DGF - RAOPEN,RAOPEN2: PHYSICAL RECORD SIZE CHANGES        
C 12 NOV 98 - GDF - PREAD,PKREAD,PWRIT,PKWRIT: CHOOSE INT. LENGTHS      
C 27 SEP 98 - MWS - RAOPEN: CHANGE FORMAT FIELD WIDTH                   
C 27 FEB 98 - AL  - PARENV: STRIP OF PATH NAME MADE MORE OPTIONAL       
C  6 JAN 98 - MWS - RAOPEN: PASS CORRECT FILE NAMES TO OTHER NODES      
C 14 OCT 97 - DGF - ADD RAWRITE ROUTINE                                 
C 28 SEP 97 - MWS - CONVERT PARALLEL CALLS FROM TCGMSG TO DDI           
C  2 SEP 97 - MWS - PKWRIT,PWRIT1,PWRIT2,PXWRIT,SQWRIT,SEQREW,          
C                   PKREAD,PREAD1,PREAD2,PXREAD: ADD ERROR BRANCHES     
C 16 JUL 97 - GNM - PUVIB: CHANGE FOR FRAGONLY OPTION                   
C 10 JUL 97 - MWS - SQREAD: ERROR BRANCH LEADS TO JOB TERMINATION       
C  7 FEB 97 - MWS - DAREAD,DAWRIT: IMPROVE ERROR MESSAGE CONTENT        
C  6 NOV 96 - MWS - ADD MQOPDA,MQDARE,MQDAWR                            
C 30 OCT 96 - HN  - SEQOPN,SEQCLO: ALLOW UNIT NUMBERS TO 99             
C 24 MAY 96 - GMC - RAOPEN: CHANGE FILE NAME TO MCDIIS FOR UNIT 30      
C  3 JAN 96 - MWS - BRING XABI'S THREE DUMMY ROUTINES HERE              
C 24 MAY 95 - MWS - INCLUDE NEW PREADP ROUTINE                          
C 21 APR 95 - MWS - SEQCLO,SEQOPN: AVOID REPEAT OPEN REQUESTS           
C 17 NOV 94 - WC  - PUVIB: FRAGMENT RELATED CHANGES                     
C 10 NOV 94 - MWS - REMOVE SOME FTNCHEK WARNINGS                        
C 11 AUG 94 - MWS - INCREASE NUMBER OF DAF RECORDS                      
C  4 AUG 94 - MWS - PREAD,PWRIT: ALLOW DOUBLE LABELS                    
C 31 MAY 94 - PRD - PUVIB: ADD FRAGMENT OUTPUT                          
C  4 JAN 94 - MWS - FIX RA ROUTINES TO WORK ON ALL NODES                
C 16 JUL 93 - MWS - ADD PXREAD ROUTINE                                  
C  5 JUN 92 - TLW - ADD ROUTINE PARENV; CHANGE ROUTINE ADVFIL TO SEQFIL;
C                   MAKE SEQENTIAL DISK OPS WORK ON DIFFERENT NODES     
C  3 APR 92 - TLW - DAREAD/WRIT: ACCOUNT FOR PURELY INTEGER RECORDS     
C  2 MAR 92 - TLW - SQREAD: BROADCAST INFORMATION                       
C 21 FEB 92 - TLW - RARD: BROADCAST INFORMATION                         
C  3 FEB 92 - TLW - DARD: BROADCAST DAF INFORMATION                     
C 30 JAN 92 - TLW - OPENDA: BROADCAST OLD DAF INFORMATION IF NEEDED     
C 11 JAN 92 - TLW - MAKE READS PARALLEL                                 
C 10 JAN 92 - TLW - MAKE NEW ROUTINE SEQREW                             
C 10 JAN 92 - TLW - MAKE OPENS AND CLOSES PARALLEL                      
C 10 JAN 92 - TLW,MWS - DELETE OPEN ROUTINES OPENCI,OPENCF,OPENDF,      
C                   OPENFM,OPENIP,OPENIR,OPENIS,OPENIW,OPENJK,OPENPK,   
C                   OPNIRC TO USE GENERIC SEQOPN; SEQCLO REPLACED CLOSJK
C  7 JAN 92 - TLW - MAKE WRITES PARALLEL; ADD COMMON PAR                
C 18 NOV 91 - MWS,JHJ - SEQOPN:INTRODUCED FROM UTILS.CODE.              
C 19 JAN 90 - MWS - INITIALIZE FILENM IN OPENIW IN UNIX SECTION         
C 21 DEC 89 - STE - ADD ERROR MESSAGES TO PKREAD,PREAD                  
C 26 SEP 89 - MWS - ADD NFT13,NFT14 TO /CIFILS/                         
C 12 APR 89 - STE - ADD FILE CIINTS (NFT14) TO OPENCI                   
C 28 MAR 89 - MWS - RAOPEN RETURNS W/O OPENING IF EXETYP.EQ.CHECK       
C 25 FEB 89 - MWS - USE FLSHBF IN PUVIB, LOWER CASE CTSS SUFFIXES       
C 27 JAN 89 - MWS - OPENCI OPENS A FILE ONLY IF ITS LRECL IS NONZERO    
C 18 JAN 89 - MWS - ADD ROUTINE OPENDF, DELETE RDGRD AND WRTGRD         
C 14 NOV 88 - MWS - PKREAD,PKWRIT,PREAD,PWRIT: NX,IX BEFORE XP,XK VALUES
C  7 OCT 88 - MWS - ADD NEW ROUTINE OPENFM                              
C  1 APR 88 - MWS - INCREASE NON-IBM RECORD SIZES FOR DASORT            
C                   TO VAX MAXIMUM OF 2047 W.P. WORDS.                  
C 15 MAR 88 - MWS - ALTER IBM CLOSE IN RAOPEN, FIX VAX LDAR IN RASIZE   
C 29 FEB 88 - STE - ADD FPS INFO TO OPENCF                              
C 18 FEB 88 - MWS - ADD NEW ROUTINE OPENCF                              
C 11 JAN 88 - MWS - ADD DIPOLE TO $VIB OUTPUT GROUP                     
C  2 NOV 87 - STE - CHANGE VAX RECL TO WORK WITH REMOTE DISKS           
C 24 AUG 87 - STE - OPENJK: FIX FPS TO WORK LIKE OTHER FILES            
C  6 AUG 87 - MWS - INCLUDE ETA VERSION                                 
C 26 JUN 87 - MWS - DON'T PRINT FORCE RESTART DATA IN PUVIB             
C  5 MAY 87 - STE - PUVIB: DIMENSION BLANK,VIBWRD                       
C 12 APR 87 - MWS - ROUTINE PUVIB TRANSFERRED FROM FORCE MODULE         
C 15 OCT 86 - MWS - ADD CRAY/CTSS SPECIFIC INSTRUCTIONS IN ORDER        
C                   TO BUILD UNIQUE FILENAMES UNDER CTSS                
C 13 OCT 86 - MWS - USE IDMY IN DAWRIT COMMON TO AVOID ARG COLLISIONS   
C 14 AUG 86 - MWS - REMOVE MULTIPLE FILE OPTION STUFF FOR IBM           
C                   IN PREAD,PWRIT,PKREAD,PKWRIT,   ADD ROUTINE OPNIRC  
C  1 AUG 86 - MWS - OBTAIN FILE NAME FROM THE ENVIRONMENT FOR CELERITY  
C  9 JUL 86 - MWS - OPEN/CLOSE STATEMENTS FOR CELERITY AND CRAY,        
C                   SANITIZE FLOATING POINT CONSTANTS, MOVE ROUTINES    
C                   ABRT,MCHPAR,SECOND,TMDATE TO MODULE GAMESS          
C 11 JUN 86 - MWS - WRITE CORRECT ERROR MESSAGE IN DAWRIT               
C 25 APR 86 - LAM - OPEN IP WITH CC='LIST' IN OPENIP ON THE VAX         
C 16 APR 86 - LAM - CHANGE FORMAT STATEMENT IN RAOPEN                   
C 26 FEB 86 - LAM - FIX PKREAD AND PWRIT TO READ INTEGRALS PROPERLY     
C                   FROM MORE THAN ONE FILE                             
C 13 FEB 86 - LAM - PKREAD, PKWRIT, PREAD, PWRIT CHANGED TO HANDLE      
C                   MULTIPLE (TAPE) FILES FOR IBM                       
C 31 JAN 86 - LAM - DEFINE LOCAL VARIABLE FIRLCL IN SECOND              
C 13 NOV 85 - STE - FPS COMMAND LINE CONTROL OF FILES                   
C 26 OCT 85 - STE - TMDATE: FPS CAN GET DATE NOW                        
C 23 OCT 85 - STE - OPENCI: CHANGE CALL                                 
C  9 SEP 85 - STE - FIX FPS OPEN IN OPENIS; CHANGE 700 IN RACLOS TO 7   
C                   ADD /IOFILE/ FOR FPS IN OPENCI                      
C  6 SEP 85 - LAM - IRECLN=512 IN OPENDA, IRECLN=4094 IN RAOPEN         
C                   AND LDAR=4090 IN RASIZE (IBM)                       
C 12 JUL 85 - MWS - 7 CHAR FILENAMES, JKFILE UNKOWN, KNOWN              
C                   DON'T OPEN IBM FT06F001                             
C 30 MAY 85 - MWS - ADD NPRINT ARGUMENT TO RAOPEN                       
C 15 APR 85 - MWS - PUT ERR= BRANCH IN RACLOS                           
C 10 APR 85 - MWS - REMOVE VAX RECL= FROM OPENCI                        
C 15 MAR 85 - MWS - ADD ROUTINE RASIZE TO SET LDAR                      
C 10 MAR 85 - MWS - ADD ROUTINE OPENCI,OPENJK,SEQCLO,                   
C                   CONVERT IBM VERSION TO VS FORTRAN OPEN/CLOSES,      
C 22 AUG 84 - STE - REMOVE UNUSED PARAMETER N FROM MSREAD/WRITE         
C 21 JUL 84 - STE - ADD IS,IPK TO UNIT10 REC=1 IN DAWRIT,OPENDA         
C 15 MAY 84 - STE - FPS TRACEBACK IN ABRT                               
C 28 APR 84 - STE - FIX MSREAD,MSWRIT TO ALLOW RESTARTS                 
C 26 MAR 84 - STE - ALLOW RUN-TIME ACCESS TO OPEN STATEMENTS            
C 23 FEB 84 - STE - INCLUDE PXWRIT                                      
C  2 FEB 84 - STE - REMOVE REWRITE LENGTH CHECK FROM RAWRIT             
C 21 JAN 84 - STE - OPEN VAX INPUT READONLY, SHARED                     
C 17 JAN 84 - STE - EVALUATE ELAPSED REAL TIME IN SECOND                
C 10 JAN 84 - STE - CHANGE STATUS=NEW TO STATUS=UNKNOWN, EXCEPT IW      
C 28 DEC 83 - STE - REMOVE INTEGER*2 IN P/PK-READ/WRIT                  
C 27 DEC 83 - STE - ADD ROUTINES SQREAD AND SQWRIT                      
C 19 DEC 83 - STE - SAVE TIME FROM INITIAL CALL TO SECOND FOR VAX       
C 15 DEC 83 - STE - ADD LRECL TO OPENIS, OPENPK                         
C  1 DEC 83 - STE - FIX NWDVAR INDEXING IN MSREAD,MSWRIT,RDGRD,WRTGRD   
C 29 NOV 83 - STE - FPS CONVERSION; CHANGE OPENDA,RAOPEN                
C 21 NOV 83 - STE - FPS CONVERSION; ADD OPENIP AND OPENPK               
C  5 NOV 83 - STE - INCLUDE ABRT,MCHPAR,SECOND,TMDATE FROM IBMLIB       
C                   FIX VAX I/O, INITIALIZE IBM TIME CORRECTLY          
C                   INSERT CORRECT TMDATE FOR VAX                       
C                   DELETE BUFINT,DARTRN,REWBF,UNPACK                   
C 13 MAY 83 - STE - INITIALIZE IBM TIMES TO ZERO.                       
C  8 MAY 83 - MWS - IOLIB FORMED FROM PIECES OF                         
C                   GAMESS, IBMLIB, AND CILIB.                          
C  2 NOV 82 - MWS - DEFINE FILE STATEMENT MADE FOR WORST CASE           
C 29 SEP 82 - MWS - CONVERT FOR IBM                                     
C                                                                       
C   IN THIS FILE,                                                       
C       *IBM REFERS TO IBM MAINFRAMES RUNNING MVS OR VM (NOT RS/6000)   
C       *VMS REFERS TO VAX OR AXP RUNNING VMS OPERATING SYSTEM          
C       *UNX REFERS TO EVERYTHING ELSE, RUNNING A UNIX O/S              
C                                                                       
C*MODULE IOLIB   *DECK DABIGIO                                          
!> @brief  This routine controls which direct access files use the new C
!> @author Jerry Boatz                                                  
!>  - June 2015 - subroutine written                                    
!>	- September 2018 - implementation into 2018 gamess                   
!>	@author Jorge Luis Galvez Vallejo                                    
!> @details  This routine was developed as part of the workaround to    
!>           avoid file record sizes larger than 2GB.  It specifies whic
!>           direct access files (according to their assigned unit numbe
!>           been set up to use the new CC I/O routines.                
!>                                                                      
!> @param IRAF Direct access file logical unit number                   
!>                                                                      
      LOGICAL FUNCTION DABIGIO(IRAF)                                    
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      PARAMETER (mxunit=299)                                            
C                                                                       
      COMMON /CCIO  / f70ok,iosize,lrecfl(mxunit)                       
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL                
c                                                                       
      LOGICAL GOPARR,MASWRK,DSKWRK,OUT                                  
      logical f70ok                                                     
C                                                                       
C     USED BY CC DIRECT ACCESS FILE ROUTINES                            
C                                                                       
C     THE RECORD SIZE OF MOST of THE FOLLOWING COUPLED CLUSTER FILES    
C     SCALE AS N**4 (N=# OF AOS), AND THEREFORE CAN POTENTIALLY EXCEED  
C     2GB, WHICH MANY COMPILERS CANNOT HANDLE; e.g., gfortran,          
C     ifort (starting with version 15), Cray compiler, Portland Group,  
C                                                                       
C     UNIT NUMBER     UNIT VARIABLE(S)     LOGICAL FILENAME(S)     SOURC
C                                                                       
C     43              NVT                  BBBB43                  ccqua
C     70              NFT825,NRESF         CCREST,AMPROCC          sever
C     71              NFT858,NFRLE         CCDIIS,ITOPNCC          sever
C     72              NFT833,INTG          CCINTS                  ccaux
C     73              NFT859,NT1           CCT1AMP,LAMB23          sever
C     74              NFT820,NT2           CCT2AMP,VHHAA           sever
C     75              NFT821,NT3           CCT3AMP                 ccaux
C     76              NVM,NFT822           CCVM,VHHAB              ccsdt
C     80              NSTAR                EOMSTAR,VMBA            eomcc
C     82              NC2,NFT830           EOMVEC2,VHPRBB          eomcc
C     84              NHC2,NFT831          EOMHC2,VHPLAA           eomcc
C     85              NHHHH,               EOMHHHH,VHPLBB          ccsdt
C                     NHH,NFT832                                        
C     87              NAX,NFT835           EOMRAMP,VHPLBA          eomcc
C     88              NAXX,NFT836          EOMRTMP,VEAA            eomcc
C     89              NDIAG,NFT837         EOMDG12,VEBB            eomcc
C     92              JCISD                MMCIVEC,VPPPP           mm23 
C     93              JCISDNX              MMCIVC1                 mm23 
C                     NFT891               INTERM1                 rohfc
C     94              NFT892,KCIF          INTERM2,MMCIITR         sever
C     96              NL2                  EOMVL2                  eomcc
C                     n/a                  ITSPACE                 roeom
C                     NEXE                 MMNEXE                  mm23 
C     97              NNAXL                EOMLVEC                 eomcc
C                     n/a                  INSTART                 roeom
C                     NREXM                MMNREXM                 mm23 
C     98              n/a                  ITSPC3                  eaipc
C                     NREXE,NHL1           MMNREXE                 mm23 
C     99              n/a                  INSTRT3                 eaipc
C                     NHL2,NF1             EOMHL2                  eomcc
C                                                                       
C                                                                       
      DATA dbugccio /8HDBUGCCIO/                                        
C                                                                       
      OUT = (exetyp.eq.dbugccio .AND. MASWRK)                           
c                                                                       
      dabigio = .false.                                                 
      if (iraf .eq. 43) dabigio = .true.                                
      if (iraf .eq. 70 .and. f70ok) dabigio = .true.                    
      if (iraf .eq. 71) dabigio = .true.                                
      if (iraf .eq. 72) dabigio = .true.                                
      if (iraf .eq. 73) dabigio = .true.                                
      if (iraf .eq. 74) dabigio = .true.                                
      if (iraf .eq. 75) dabigio = .true.                                
      if (iraf .eq. 76) dabigio = .true.                                
      if (iraf .eq. 77) dabigio = .true.                                
      if (iraf .eq. 80) dabigio = .true.                                
      if (iraf .eq. 82) dabigio = .true.                                
      if (iraf .eq. 84) dabigio = .true.                                
      if (iraf .eq. 85) dabigio = .true.                                
      if (iraf .eq. 87) dabigio = .true.                                
      if (iraf .eq. 88) dabigio = .true.                                
      if (iraf .eq. 89) dabigio = .true.                                
      if (iraf .eq. 92) dabigio = .true.                                
      if (iraf .eq. 93) dabigio = .true.                                
      if (iraf .eq. 94) dabigio = .true.                                
      if (iraf .eq. 96) dabigio = .true.                                
      if (iraf .eq. 97) dabigio = .true.                                
      if (iraf .eq. 98) dabigio = .true.                                
      if (iraf .eq. 99) dabigio = .true.                                
c                                                                       
c     --- disable all file chunking if iosize=0                         
c                                                                       
      if (iosize .eq. 0) dabigio = .false.                              
c                                                                       
      if (out) write(6,*)'DBUGCCIO, in function DABIGIO: iraf,dabigio=',
     *                                                   iraf,dabigio   
c                                                                       
c     --- Uncomment these 4 lines to completely disable new CC I/O routi
c                                                                       
c-2gb-off    dabigio=.false.                                            
c-2gb-off    if (maswrk) then                                           
c-2gb-off       write(6,*) '*** NOTE:  New CC I/O routines are disabled 
c-2gb-off    end if                                                     
      RETURN                                                            
      END                                                               
C*MODULE IOLIB   *DECK CCDAREAD                                         
!> @brief  This routine controls the "chunking" of a single large record
!> @author Jerry Boatz                                                  
!>  - June 2015 - subroutine written                                    
!>                                                                      
!> @details  This routine was developed as part of the workaround to    
!>           avoid file record sizes larger than 2GB.  It calls CCDARED2
!>           times, once for each smaller physical record comprising the
!>           large record passed in as input.                           
!> @param IRAF Direct access file logical unit number                   
!> @param V  Input, double precision array to be read in                
!> @param LENLOG  Logical length of array V to be read in               
!> @param NREC    Logical record number of record to be read in         
!> @param Ioerr  Flag controlling how I/O errors are to be handled      
!> @param from   Character string passed in from calling routine        
!>                                                                      
      SUBROUTINE CCDAREAD(IRAF,V,LENLOG,NREC,ioerr,from)                
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
c                                                                       
      PARAMETER (mxunit=299)                                            
      character*(*) from                                                
      DIMENSION V(lenlog)                                               
C                                                                       
      LOGICAL GOPARR,DSKWRK,MASWRK,out,dabigio,f70ok                    
C                                                                       
C                                                                       
      COMMON /CCIO  / f70ok,iosize,lrecfl(mxunit)                       
c-remove-fm                                                             
c     COMMON /FMCOM / XX(1)                                             
c-remove-fm                                                             
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)                
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL                
C                                                                       
      DATA dbugccio /8HDBUGCCIO/                                        
C                                                                       
      OUT = (exetyp.eq.dbugccio .AND. MASWRK)                           
c                                                                       
c     --- arugment error checking                                       
c                                                                       
      ierr = 0                                                          
      if (maswrk) then                                                  
         IF (IRAF .LT. 0) then                                          
            WRITE(IW,*) 'BOGUS CCDAREAD, called from ',from,            
     *                  ': unit number iraf= ',iraf                     
            ierr = ierr + 1                                             
         end if                                                         
         if (lenlog .lt. 1) then                                        
            WRITE(IW,*) 'BOGUS CCDAREAD, called from ',from,            
     *                  ': record length lenlog= ',lenlog               
            ierr = ierr + 1                                             
         end if                                                         
         if (nrec .lt. 1) then                                          
            WRITE(IW,*) 'BOGUS CCDAREAD, called from ',from,            
     *                  ': logical record number nrec= ',nrec           
            ierr = ierr + 1                                             
         end if                                                         
      end if                                                            
C                      *********                                        
      if (ierr .gt. 0) call abrt                                        
C                      *********                                        
c                                                                       
c     --- retrieve the number of physical records per logical record fro
c                                                                       
      nprlog = lrecfl(iraf)                                             
c                                                                       
c     --- if nprlog=1, this file is not being chunked and therefore can 
c     --- to ccdared2 without having to set up an I/O buffer.           
c                                                                       
      if (nprlog .eq. 1) then                                           
         if (out) write(iw,9000) iraf,from,nrec,lenlog                  
         call ccdared2(v,lenlog,iraf,nrec,ioerr,'CCDAREAD:1')           
C        *********                                                      
         GO TO 999                                                      
C        *********                                                      
      end if                                                            
c                                                                       
c     --- allocate memory for temporary I/O buffer.                     
c                                                                       
c-remove-fm                                                             
c-new                                                                   
c      call valfm(loadfm)                                               
c      liobuf = loadfm + 1                                              
c      needio = iosize                                                  
c      call getfm(needio)                                               
c-new                                                                   
c-remove-fm                                                             
c                                                                       
c     --- set the physical record length to the user-specified value of 
c                                                                       
      lenphz=iosize                                                     
      if (lenphz .gt. lenlog) lenphz=lenlog                             
c                                                                       
c     --- two more error checks                                         
c                                                                       
      if (maswrk) then                                                  
         if (lenphz .lt. 1) then                                        
            WRITE(IW,*) 'BOGUS CCDAREAD, called from ',from,            
     *                  ': physical record length lenphz= ',lenphz      
            ierr = ierr + 1                                             
         end if                                                         
         if (nprlog .lt. 1) then                                        
            WRITE(IW,*) 'BOGUS CCDAREAD, called from ',from,            
     *                  ': number of physical records per logical ',    
     *                  'record nprlog= ',nprlog                        
            ierr = ierr + 1                                             
         end if                                                         
      end if                                                            
C                                                                       
c     --- determine number of physical records per logical record       
C     -intrec- is the number of complete physical records per logical re
C     -lftovr- is the length of the final, partial physical record (if a
c                                                                       
c     --- Note that this allows for the possibility that the logical len
c     --- of this specific record (lenlog) is smaller than the initial l
c     --- length specified when this file was opened.                   
c                                                                       
      intrec = lenlog/lenphz                                            
      lftovr = lenlog - (intrec * lenphz)                               
c                                                                       
c     --- if total number of physical records exceeds original value det
c     --- at file creation, something has gone wrong!                   
c                                                                       
      nphys = intrec                                                    
      if (lftovr .gt. 0) nphys=nphys+1                                  
      if (nphys .gt. nprlog) then                                       
         if (maswrk) then                                               
            WRITE(IW,*) '*** ERROR in CCDAREAD *** '                    
            WRITE(IW,*) '  Number of physical records (',nphys,') for ',
     *                  'logical record ',nrec,' is larger than the ',  
     *                  'number set at file creation (',nprlog,'.)'     
            ierr = ierr + 1                                             
         end if                                                         
      end if                                                            
C                      *********                                        
      if (ierr .gt. 0) call abrt                                        
C                      *********                                        
c                                                                       
c     --- map the logical record number -nrec- to the first physical rec
c                                                                       
      iphz1 = 1 + (nrec-1)*nprlog                                       
c                                                                       
C     ------ READ the -intrec- complete physical records                
c                                                                       
      ibeg = 1                                                          
      iphzrc = iphz1                                                    
      if (out) write(iw,9100) iraf,from,intrec,lenphz,nrec              
      DO Ired=1,intrec                                                  
c-remove-fm                                                             
c        CALL CCDARED2(XX(liobuf),lenphz,IRAF,iphzrc,ioerr,'CCDAREAD:2')
c        call XCOPY(lenphz,xx(liobuf),1,v(ibeg),1)                      
         CALL CCDARED2(V(ibeg),lenphz,IRAF,iphzrc,ioerr,'CCDAREAD:2')   
c-remove-fm                                                             
         ibeg = ibeg + lenphz                                           
         iphzrc = iphzrc+1                                              
      end do                                                            
c                                                                       
c     ----- read the remaining partial physical record (if any), using X
c     ----- instead of reading into V directly.  This avoids the potenti
c     ----- of changing the contents of V immediately following V(ibeg+l
c                                                                       
      if (lftovr .gt. 0) then                                           
         if (out) write(iw,9120) iraf,from,iphzrc,lftovr,nrec           
c-remove-fm                                                             
c        call CCDARED2(XX(liobuf),lftovr,iraf,iphzrc,ioerr,'CCDAREAD:3')
c        call XCOPY(lftovr,xx(liobuf),1,v(ibeg),1)                      
         call CCDARED2(V(ibeg),lftovr,iraf,iphzrc,ioerr,'CCDAREAD:3')   
c-remove-fm                                                             
      end if                                                            
C                                                                       
c-remove-fm                                                             
c-new                                                                   
C                                                                       
C     ---- return the memory for the I/O buffer                         
c                                                                       
c     call retfm(needio)                                                
c-new                                                                   
c-remove-fm                                                             
  999 CONTINUE                                                          
      RETURN                                                            
 9000 FORMAT(1X,'DBUGCCIO:unit=',I4,':  In routine ccdaread:1, called ',
     *          'from ',A,', reading logical record # ',I8,' with ',    
     *          'length ',I10)                                          
 9100 FORMAT(1X,'DBUGCCIO:unit=',I4,':  In routine ccdaread:2, called ',
     *          'from ',A,', preparing to read ',I8,' physical ',       
     *          ' records with length ',I8,', for logical record ',I8)  
 9120 FORMAT(1X,'DBUGCCIO:unit=',I4,':  In routine ccdaread:3, called ',
     *          'from ',A,', reading physical record # ',I8,' with ',   
     *          'length ',I8,' for logical record ',I8)                 
      end                                                               
C*MODULE IOLIB   *DECK CCDARED2                                         
!> @brief  This routine is called by CCDAREAD to read a single physical 
!> @author Jerry Boatz                                                  
!>  - June 2015 - subroutine written                                    
!>                                                                      
!> @details  This routine was developed as part of the workaround to    
!>           avoid file record sizes larger than 2GB.  It is a low-level
!>           counterpart to CCDAREAD, handling reading of just a single 
!>           All "non-chunked" I/O is also handed off to this routine.  
!> @param V  Input, double precision array to be read in                
!> @param LEN  Length of array V to be read in                          
!> @param IRAF Direct access file logical unit number                   
!> @param NS   Record number to be read                                 
!> @param Ioerr  Flag controlling how I/O errors are to be handled      
!> @param from   Character string passed in from calling routine        
!>                                                                      
      SUBROUTINE CCDARED2(V,LEN,IRAF,NS,ioerr,from)                     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      LOGICAL GOPARR,DSKWRK,MASWRK,chk2gb                               
      DIMENSION V(LEN)                                                  
      character*(*) from                                                
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)                
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL                
C                                                                       
      INTEGER IOCODE                                                    
C                                                                       
C     -IOERR-  is a flag to determine how I/O errors are handled.       
C              IOERR = 0 means all errors result in a call to ABRT.     
C              If IOERR.ne.0, then the IOSTAT error code is returned in 
C                                                                       
      logical out                                                       
C                                                                       
      DATA dbugccio /8HDBUGCCIO/                                        
C                                                                       
      OUT = exetyp.eq.dbugccio .AND. MASWRK                             
C                                                                       
C     ----- READ A PHYSICAL RECORD FROM IRAF -----                      
C                                                                       
c     if (out) write(iw,9000) iraf,FROM,NS,LEN                          
c                                                                       
c     ----- Since "len" is in units of words, must check 8*len for 2GB l
c                                                                       
      len8 = len*8                                                      
      if (chk2gb(len8)) then                                            
         if (maswrk) then                                               
            write(iw,*) '*** Warning in ccdared2:  Attempting to read ',
     *                  'direct access record of length ',len8,' bytes,'
            write(iw,*) 'which is .ge. 2GB.  This may cause an abend.'  
            write(iw,*) 'Unit number = ',iraf                           
         end if                                                         
      end if                                                            
      iocode = 0                                                        
      if (ioerr .eq. 0) then                                            
         READ(UNIT=IRAF,REC=NS) V                                       
      else                                                              
         READ(UNIT=IRAF,REC=NS,iostat=iocode) V                         
      end if                                                            
      if (iocode .ne. 0) then                                           
         if (ioerr .eq. 0) then                                         
            if (maswrk)                                                 
     *      write(iw,*) '*** Error reading unit=',iraf,', called from ',
     *                   from,':  iostat error code = ',iocode          
c           *********                                                   
            call abrt                                                   
c           *********                                                   
         end if                                                         
      end if                                                            
      ioerr = iocode                                                    
      RETURN                                                            
c                                                                       
 9000 format(1x,'DBUGCCIO:unit=',I4,': In routine CCDARED2, called from'
     *       ,A,' READING PHYSICAL RECORD # ',I8,' OF LENGTH ',I8)      
      END                                                               
C*MODULE IOLIB   *DECK CCDAWRIT                                         
!> @brief  This routine controls the "chunking" of a single large record
!> @author Jerry Boatz                                                  
!>  - June 2015 - subroutine written                                    
!>                                                                      
!> @details  This routine was developed as part of the workaround to    
!>           avoid file record sizes larger than 2GB.  It calls CCDAWRT2
!>           times, once for each smaller physical record comprising the
!>           large record passed in as input.                           
!> @param IRAF    Direct access file logical unit number                
!> @param V       Output, double precision array to be written to unit I
!> @param LENLOG  Logical length of array V to be written               
!> @param NREC    Logical record number of record to be written         
!> @param Ioerr   Flag controlling how I/O errors are to be handled     
!> @param from    Character string passed in from calling routine       
!>                                                                      
      SUBROUTINE CCDAWRIT(IRAF,V,lenlog,NREC,ioerr,from)                
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
c                                                                       
      PARAMETER (mxunit=299)                                            
      character*(*) from                                                
      DIMENSION V(lenlog)                                               
      LOGICAL GOPARR,DSKWRK,MASWRK,out,dabigio,f70ok                    
C                                                                       
      COMMON /CCIO  / f70ok,iosize,lrecfl(mxunit)                       
c-remove-fm                                                             
c     COMMON /FMCOM / XX(1)                                             
c-remove-fm                                                             
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)                
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL                
C                                                                       
      DATA dbugccio /8HDBUGCCIO/                                        
C                                                                       
      OUT = (exetyp.eq.dbugccio .AND. MASWRK)                           
c                                                                       
c     --- error checking                                                
c                                                                       
      ierr = 0                                                          
      if (maswrk) then                                                  
         IF (IRAF .LT. 0) then                                          
            WRITE(IW,*) '*** BOGUS call to CCDAWRIT, from ',from,       
     *                  ': unit number iraf= ',iraf                     
            ierr = ierr + 1                                             
         end if                                                         
         if (lenlog .lt. 1) then                                        
            WRITE(IW,*) '*** BOGUS call to CCDAWRIT, from ',from,       
     *                  ': unit #, record length lenlog= ',iraf,lenlog  
            ierr = ierr + 1                                             
         end if                                                         
         if (nrec .lt. 1) then                                          
            WRITE(IW,*) '*** BOGUS call to CCDAWRIT, from ',from,       
     *                  ': unit #, record number nrec= ',iraf,nrec      
            ierr = ierr + 1                                             
         end if                                                         
      end if                                                            
C                      *********                                        
      if (ierr .gt. 0) call abrt                                        
C                      *********                                        
c                                                                       
c     --- retrieve the max number of physical records per logical record
c                                                                       
      nprlog = lrecfl(iraf)                                             
c                                                                       
c     --- if nprlog = 1, then this file is not being chunked            
c                                                                       
      if (nprlog .eq. 1) then                                           
         if (out) write(iw,9100) iraf,from,nrec,lenlog                  
         CALL CCDAWRT2(V,lenlog,IRAF,nrec,ioerr,'CCDAWRIT:1')           
C        *********                                                      
         go to 999                                                      
C        *********                                                      
      end if                                                            
c                                                                       
c     --- set the physical record length to the user-specified value of 
c                                                                       
      lenphz=iosize                                                     
      if (lenphz .gt. lenlog) lenphz=lenlog                             
c                                                                       
c     --- map the logical record number -nrec- to the first physical rec
c                                                                       
      iphz1 = 1 + (nrec-1)*nprlog                                       
c                                                                       
c     --- two more error checks                                         
c                                                                       
      ierr = 0                                                          
      if (maswrk) then                                                  
         if (lenphz .lt. 1) then                                        
            WRITE(IW,*) '*** BOGUS call to CCDAWRIT, from ',from,       
     *                  ': unit #, phys. record length lenphz= ',       
     *          iraf,lenphz                                             
            ierr = ierr + 1                                             
         end if                                                         
         if (nprlog .lt. 1) then                                        
            WRITE(IW,*) '*** BOGUS call to CCDAWRIT, from ',from,       
     *                  ': iraf, number of physical records per ',      
     *                  'logical record nprlog= ',iraf,nprlog           
            ierr = ierr + 1                                             
         end if                                                         
      end if                                                            
C                                                                       
c     --- determine number of physical records per logical record; intre
C     -intrec- is the number of complete physical records per logical re
C     -lftovr- is the length of the final, partial physical record (if a
c                                                                       
c     --- note that this allows for the possibility that the logical len
c     --- this specific record (lenlog) is less than the logical record 
c     --- when the file was opened.  Of course, lenlog can never exceed 
c     --- logical record length!                                        
c                                                                       
      intrec = lenlog/lenphz                                            
      lftovr = lenlog - (intrec * lenphz)                               
c                                                                       
c     --- if total number of physical records exceeds original value det
c     --- at file creation, something has gone wrong!                   
c                                                                       
      nphys = intrec                                                    
      if (lftovr .gt. 0) nphys=nphys+1                                  
      if (nphys .gt. nprlog) then                                       
         if (maswrk) then                                               
            WRITE(IW,*) '*** ERROR in CCDAWRIT *** '                    
            WRITE(IW,*) '  Unit # = ',iraf,'; number of physical ',     
     *                  'records (',nphys,') for logical record ',nrec, 
     *                  ' is larger than the number set at file ',      
     *                  'creation (',nprlog,'.)'                        
            ierr = ierr + 1                                             
         end if                                                         
      end if                                                            
C                      *********                                        
      if (ierr .gt. 0) call abrt                                        
C                      *********                                        
c                                                                       
C     ------ WRITE the -intrec- complete physical records               
c                                                                       
      ibeg = 1                                                          
      iphzrc = iphz1                                                    
      if (out) write(iw,9110) iraf,from,intrec,lenphz,nrec              
      DO IWRT=1,intrec                                                  
         CALL CCDAWRT2(V(ibeg),lenphz,IRAF,iphzrc,ioerr,'CCDAWRIT:2')   
         ibeg = ibeg + lenphz                                           
         iphzrc = iphzrc+1                                              
      end do                                                            
c                                                                       
      if (lftovr .gt. 0)  then                                          
c                                                                       
c        --- The last, partial record is written using fast memory for  
c        --- a temporary buffer.  This is done to prevent a possible    
c        --- segmentation fault caused by accessing the vector V beyond 
c        --- the last element V(ibeg+lftovr)...                         
c                                                                       
c-remove-fm                                                             
c        call valfm(loadfm)                                             
c        liobuf = loadfm + 1                                            
c        needio = iosize                                                
c        call getfm(needio)                                             
c                                                                       
c        --- copy the leftover section of V to the buffer and then hand 
c                                                                       
c        call XCOPY(lftovr,v(ibeg),1,xx(liobuf),1)                      
c        CALL CCDAWRT2(xx(liobuf),lftovr,IRAF,iphzrc,ioerr,'CCDAWRIT:3')
c                                                                       
c        call retfm(needio)                                             
c-remove-fm                                                             
         if (out) write(iw,9120) iraf,from,iphzrc,lftovr,nrec           
         CALL CCDAWRT2(V(ibeg),lftovr,IRAF,iphzrc,ioerr,'CCDAWRIT:3')   
      end if                                                            
C                                                                       
  999 CONTINUE                                                          
      RETURN                                                            
 9000 FORMAT(1X,'*** ERROR in iolib:ccdawrit ****',/,                   
     *          '    Direct access file with unit number = ',I5,        
     *          ' not been enabled for chunking.',/,                    
     *          '    Its record length (',I10,') exceeds the maximum',  
     *          'size specified by the IOSIZE (',I10,') input ',        
     *          'variable in $CCINP; calling ABRT...')                  
 9100 FORMAT(1X,'DBUGCCIO:unit=',I4,':  In routine ccdawrit:1, called ',
     *          'from ',A,', writing logical record # ',I8,' with ',    
     *          'length ',I10)                                          
 9110 FORMAT(1X,'DBUGCCIO:unit=',I4,':  In routine ccdawrit:2, called ',
     *          'from ',A,', preparing to write ',I8,' physical ',      
     *          'records with length ',I8,', for logical record ',I8)   
 9120 FORMAT(1X,'DBUGCCIO:unit=',I4,':  In routine ccdawrit:3, called ',
     *          'from ',A,', writing physical record # ',I8,' with ',   
     *          'length ',I8,', for logical record ',I8)                
      end                                                               
C*MODULE IOLIB   *DECK CCDAWRT2                                         
!> @brief  This routine is called by CCDAWRIT to write a single physical
!> @author Jerry Boatz                                                  
!>  - June 2015 - subroutine written                                    
!>                                                                      
!> @details  This routine was developed as part of the workaround to    
!>           avoid file record sizes larger than 2GB.  It is a low-level
!>           counterpart to CCDAWRIT, handling reading of just a single 
!>           All "non-chunked" I/O is also handed off to this routine.  
!> @param V      Output, double precision array to be written to disk   
!> @param LEN    Length of array V to be written                        
!> @param IRAF   Direct access file logical unit number                 
!> @param NS     Record number of the direct access record to be written
!> @param Ioerr  Flag controlling how I/O errors are to be handled      
!> @param from   Character string passed in from calling routine        
!>                                                                      
      SUBROUTINE CCDAWRT2(V,LEN,IRAF,NS,IOERR,from)                     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION V(LEN)                                                  
c                                                                       
      LOGICAL GOPARR,DSKWRK,MASWRK,chk2gb                               
      character*(*) from                                                
C                                                                       
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)                
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
C                                                                       
C     -IOERR-  is a flag to determine how I/O errors are handled.       
C     IOERR = 0 means all errors result in a call to ABRT.              
C     If IOERR.ne.0, then the IOSTAT error code is returned in IOERR.   
c                                                                       
      integer iocode                                                    
C                                                                       
C     ----- WRITE A PHYSICAL RECORD ON IRAF -----                       
C                                                                       
      iocode = 0                                                        
c                                                                       
c     ---- since "len" is in units of words, must check 8*len for 2GB li
c                                                                       
      len8 = len*8                                                      
      if (chk2gb(len8)) then                                            
         if (maswrk) then                                               
            write(iw,*) '*** Warning in ccdawrt2: Attempting to write ',
     *                  'direct access record of length ',len8,' bytes,'
            write(iw,*) 'which is .ge. 2GB.  This may cause an abend.'  
            write(iw,*) 'Unit number = ',iraf                           
         end if                                                         
      end if                                                            
      if (ioerr .eq. 0) then                                            
         WRITE (UNIT=IRAF, REC=NS ) V                                   
      else                                                              
         WRITE (UNIT=IRAF, REC=NS, IOSTAT=IOCODE) V                     
      end if                                                            
      if (iocode .ne. 0) then                                           
         if (ioerr .eq. 0) then                                         
            if (maswrk)                                                 
     *      write(iw,*) '*** Error writing unit=',iraf,', called from ',
     *                   from,':  iostat error code = ',iocode          
c           *********                                                   
            call abrt                                                   
c           *********                                                   
         end if                                                         
      end if                                                            
      ioerr = iocode                                                    
      RETURN                                                            
      END                                                               
C*MODULE IOLIB   *DECK CHGMIU                                           
      SUBROUTINE CHGMIU(IR,IW)                                          
C     SATISFY FTNCHEK'S STALWARTH ZEAL                                  
      IF(IR.LT.0) WRITE(6,*) IR,IW                                      
      RETURN                                                            
      END                                                               
C*MODULE IOLIB   *DECK CHMDAT                                           
      SUBROUTINE CHMDAT(AATOM,AZNUC,CORD,NAT)                           
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      CHARACTER*10 AATOM(*)                                             
      DIMENSION AZNUC(*),CORD(*)                                        
C     SATISFY FTNCHEK'S STALWARTH ZEAL                                  
      IF(NAT.LT.0) WRITE(6,*) AATOM(1),AZNUC(1),CORD(1),NAT             
      RETURN                                                            
      END                                                               
C*MODULE IOLIB   *DECK FDNAI                                            
      DOUBLE PRECISION FUNCTION FDNAI(PLAMBDA,AI,AJ,                    
     *                                L1A,M1A,N1A,L2B,M2B,N2B,          
     *                                XI,YI,ZI,XJ,YJ,ZJ,CX,CY,CZ)       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C                                                                       
C        THESE THREE ROUTINES IMPLEMENT THE SCREENED INTEGRAL           
C        METHOD OF JOSE UGALDE AND XABI LOPEZ, WHICH IS NOT             
C        AT PRESENT IN THE STANDARD GAMESS CODE.                        
C                                                                       
      FDNAI = 0.0D+00                                                   
C     SATISFY FTNCHEK'S STALWARTH ZEAL                                  
      IF(FDNAI.GT.1) WRITE(6,*) PLAMBDA,AI,AJ,L1A,M1A,N1A,L2B,M2B,N2B,  
     *                          XI,YI,ZI,XJ,YJ,ZJ,CX,CY,CZ              
      RETURN                                                            
      END                                                               
C*MODULE IOLIB   *DECK STWOEI                                           
      SUBROUTINE STWOEI(SCFTYP,DIRSCF,DIRNLO,DIRTRF,                    
     *               INTG76,SCHWRZ,NINT,NSCHWZ,L1,L2,                   
     *               BUFP,BUFK,IBUF,NINTMX,                             
     *               XINTS,NSH2,GHONDO,MAXG,DDIJ,                       
     *               IA,DA,FA,DB,FB,DSH,DNLO,FNLO)                      
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C     SATISFY FTNCHEK'S STALWARTH ZEAL                                  
      IF(SCFTYP.EQ.1) WRITE(6,*) SCFTYP,DIRSCF,DIRNLO,DIRTRF,           
     *                           INTG76,SCHWRZ,NINT,NSCHWZ,L1,L2,       
     *                           BUFP,BUFK,IBUF,NINTMX,                 
     *                           XINTS,NSH2,GHONDO,MAXG,DDIJ,           
     *                           IA,DA,FA,DB,FB,DSH,DNLO,FNLO           
      RETURN                                                            
      END                                                               
C*MODULE IOLIB   *DECK INITFCTS                                         
      SUBROUTINE INITFCTS                                               
      RETURN                                                            
      END                                                               
C*MODULE IOLIB   *DECK DARD                                             
      SUBROUTINE DARD(V,LEN,IDAF,NS,IDTYP)                              
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      LOGICAL GOPARR,DSKWRK,MASWRK                                      
      DIMENSION V(LEN)                                                  
      COMMON /MACHIN/ NWDVAR,MAXFM,MAXSM,LIMFM,LIMSM                    
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
C                                                                       
C       READ A PHYSICAL RECORD FROM THE DAF                             
C                                                                       
      IF (MASWRK) READ (UNIT=IDAF, REC=NS) V                            
      ITYP = 100 + IDAF                                                 
      IF (GOPARR) THEN                                                  
         IF (IDTYP.NE.1) THEN                                           
            CALL DDI_BCAST(ITYP,'F',V,LEN,MASTER)                       
         ELSE                                                           
            CALL DDI_BCAST(ITYP,'I',V,NWDVAR*LEN,MASTER)                
         END IF                                                         
      END IF                                                            
      RETURN                                                            
      END                                                               
C*MODULE IOLIB   *DECK DAREAD                                           
!> @brief Read a logical record from the dictionary file                
!> @details The direct access (DA) file is designed to store small data 
!>   quantities, such as one electron integrals, in a manner that       
!>   allows this data to be recovered at whim, in any order.            
!>     DA READ/WRIT arguments are (IDAF,IODA,V,LEN,NREC,IDTYP)          
!>     COMMON /DAIOLN/ IRECLN,IRECST,NRECUS,IFILEN(950),MERF10(950)     
!>     COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)             
!>                                                                      
!>   V is the data item to be written or read, and LEN is V's length.   
!>   Typical data is of the size of the basis set, or a symmetric       
!>   or square matrix of the size of the basis set.  The data may       
!>   be floating point or integer, according to IDTYP being 0 or 1.     
!>   Each unique kind of data is identified with a unique integer,      
!>   namely NREC, running from 1 to 950.  A programmer need worry       
!>   only about the final four arguments.                               
!>                                                                      
!>   Person's using this routine should pass in the correct unit        
!>   number IDAF and a bookkeeping array IODA from /IOFILE/, but        
!>   should never change these values.  IODA should never have been     
!>   a calling argument, but is, like IDAF, so just pass them in.       
!>                                                                      
!>   Since the data length is not fixed, we have "logical" records      
!>   which are written to the disk as one or more "physical" records.   
!>                                                                      
!>   NREC identifies the logical record number, such as 12 for the      
!>   overlap integrals.  The value must be unique, and there is         
!>   a list of all data records in the Programmer's Guide chapter       
!>   of the GAMESS manual.  Data may be written in any order, so        
!>   the only requirement on NREC is it be unique, and less than        
!>   or equal to 950.  It is vital to record in the programmer's        
!>   manual any new logical record, since they must be unique.          
!>                                                                      
!>   Someone using this routine can forget about physical records,      
!>   and not be concerned with anything below about how it works.       
!>                                                                      
!>   Physical records are of fixed size, a language requirement,        
!>   as well as a practical one in making the data location within      
!>   the file an easily computed quantity.                              
!>   A new logical record is placed at the end of all existing          
!>   physical records, using as many physical records to store the      
!>   data as is necessary.  Only part of the final physical             
!>   record will be used, any remainder wasted disk space.              
!>   Since additional logical records will be appended later, the       
!>   data length LEN requested at first write is stored in IFILEN,      
!>   and it is strictly enforced that later writes use the same LEN.    
!>   It is clear that one wants to be sure not to write into the        
!>   trailing data, as the direct access file grows.  Since writing     
!>   less data later is far more likely to be a calling argument bug    
!>   than  anything else, writing must always use the initial LEN.      
!>   It is permitted to read less than was written, such as asking      
!>   only for the occupied orbitals of a full orbital deck, but         
!>   a request to read more than the original LEN is a fatal error.     
!>                                                                      
!>   IODA, which should never be touched by the user, stores the        
!>   first physical record number of the given NREC logical record.     
!>   This, together with IFILEN, allows reading the consecutive         
!>   physical records to string together a complete logical record.     
!>   Anyone using the direct access write/read routines does not        
!>   need to be concerned about the underlying physical records.        
!>   Their number, and their size, has been increased perhaps           
!>   three times over the years, at about 15 year intervals.            
!>                                                                      
!>   In regards parallel computing, the MASTER rank maintains           
!>   the disk file, so it is the only one that writes data out.         
!>   For a read, the MASTER will read the disk file and broadcast       
!>   to all other ranks.  This is the purpose of the data type,         
!>   to generate floating point or integer broadcasts.  One             
!>   consequence of having only the master rank write, is that          
!>   all ranks are supplied with bitwise identical data if the          
!>   data is read back in.  This can help to control roundoff           
!>   issues, by annointing the MASTER rank's data as "best".            
!>                                                                      
!>   There are parallel computer systems sold that are so pathetic      
!>   at I/O that even having only one rank perform I/O is slow.         
!>   An option by Dmitri Federov stores the data into parts             
!>   of the system memory.  Mike Schmidt doesn't really know how        
!>   that works, so read the code.  He guess that MERF10 is             
!>   related to storing the data in memory, and that that is            
!>   also the purpose of various new data in /MACHSW/.                  
!>                                                                      
!> @param IDAF the unit number of the dictionary file                   
!> @param IODA                                                          
!> @param V output double precision array                               
!> @param LEN lenght of V                                               
!> @param NREC                                                          
!> @param IDTYP                                                         
!> @author Unkown                                                       
!> @date September 2010                                                 
!> - Support for the CIMDCT file used in correlated subsystem calculatio
      SUBROUTINE DAREAD(IDAF,IODA,V,LEN,NREC,IDTYP)                     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DOUBLE PRECISION NEDA                                             
C                                                                       
      LOGICAL GOPARR,DSKWRK,MASWRK                                      
C                                                                       
      DIMENSION V(LEN),IODA(950)                                        
C                                                                       
      COMMON /DAIOLN/ IRECLN,IRECST,NRECUS,IFILEN(950),MERF10(950)      
      COMMON /FMCOM / X(1)                                              
      COMMON /IOFILE/ IR,IW,IP,IIS,IPK,IDAFX,NAV,IODAX(950)             
      COMMON /MACHIN/ NWDVAR,MAXFM,MAXSM,LIMFM,LIMSM                    
      COMMON /MACHSW/ KDIAG,ICORFL,IXDR,modio,mem10,lpnt10,mem10m       
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL                
C                                                                       
      DATA NEDA   /8HNEDA    /                                          
C                                                                       
C         READ A LOGICAL RECORD FROM THE DAF DICTIONARY FILE            
C         A LOGICAL RECORD MAY SPAN SEVERAL PHYSICAL RECORDS.           
C                                                                       
C         CALLING ARGUMENT -IDTYP- IS 0 OR 1 FOR FLOATING POINT         
C         OR INTEGER RECORDS.  NO OTHER DATA TYPE IS ALLOWED!           
C         RECORDS MUST BE PURELY FLOATING POINT OR PURELY INTEGER.      
C         NO MATTER WHAT -IDTYP- IS, THE -LEN- OF THE RECORD            
C         MUST BE GIVEN IN TERMS OF FLOATING POINT WORDS.  IN           
C         TURN, THIS MEANS THAT INTEGER RECORDS ON 32 BIT MACHINES      
C         MUST CONTAIN AN EVEN NUMBER OF INTEGERS.                      
C                                                                       
c     write(6,*) 'wwwdaread',NREC,LEN                                   
      IF(IDAF.NE.10) THEN                                               
         if(idaf.ne.194) then                                           
            CALL RAREAD(IDAF,IODA,V,LEN,NREC,NAV)                       
            RETURN                                                      
         endif                                                          
      ENDIF                                                             
C                                                                       
      N = IODA(NREC)                                                    
      IF(N.EQ.-1) GO TO 800                                             
C                                                                       
C         TRAP ANY READS THAT ARE TOO LARGE, BUT PERMIT READS WITH      
C         SIZES BETWEEN 1 AND THE AVAILABLE RECORD TO PROCEED.          
C         NOTE: NBO RUNS MAY USE ZERO RECORD LENGTHS DURING NBO         
C         STEPS, WHEN -RUNTYP- WILL HAVE BEEN MANIPULATED TO -NEDA-.    
C                                                                       
      IF(RUNTYP.EQ.NEDA  .AND.  LEN.EQ.0) GO TO 500                     
C                                                                       
      IF(LEN.LE.0           ) GO TO 810                                 
      IF(LEN.GT.IFILEN(NREC)) GO TO 820                                 
C                                                                       
      if(mem10.ne.0) then                                               
c        MEN=mod(nrec,nproc)                                            
         men=merf10(nrec)                                               
c        write(6,*) 'wwwread',nrec,men,me                               
         if(men.lt.0.or.men.ge.nproc) then                              
           write(6,*) 'MEM10 internal error',men,nproc,nrec             
           call abrt                                                    
         endif                                                          
c        sanity check                                                   
         if(men.eq.me) then                                             
            call XCOPY(len,x(lpnt10+N-1),1,v,1)                         
         endif                                                          
         ITYP = 100 + IDAF                                              
         IF (GOPARR) THEN                                               
            IF (IDTYP.NE.1) THEN                                        
               CALL DDI_BCAST(ITYP,'F',V,LEN,MEN)                       
            ELSE                                                        
               CALL DDI_BCAST(ITYP,'I',V,NWDVAR*LEN,MEN)                
            END IF                                                      
         END IF                                                         
c        write(6,*) 'Read F10 record',me,men,nrec,n                     
         return                                                         
      endif                                                             
C                                                                       
      IS = -IRECLN + 1                                                  
      NS = N                                                            
      LENT = LEN                                                        
  100 CONTINUE                                                          
         IS = IS + IRECLN                                               
         IF = IS + LENT - 1                                             
         IF ((IF-IS+1) .GT. IRECLN) IF = IS + IRECLN - 1                
         NSP = NS                                                       
         LENW = IF - IS + 1                                             
         CALL DARD(V(IS),LENW,IDAF,NSP,IDTYP)                           
         LENT = LENT - IRECLN                                           
         NS = NS + 1                                                    
         N = NS                                                         
      IF (LENT .GE. 1) GO TO 100                                        
C                                                                       
  500 CONTINUE                                                          
c     if(nrec.eq.14) write(6,*) 'FFread',v(1),v(2)                      
      RETURN                                                            
C                                                                       
  800 CONTINUE                                                          
      IF(MASWRK) WRITE(IW,9000) NREC,LEN                                
      GO TO 890                                                         
  810 CONTINUE                                                          
      IF(MASWRK) WRITE(IW,9010) NREC,LEN                                
      GO TO 890                                                         
  820 CONTINUE                                                          
      IF(MASWRK) WRITE(IW,9020) LEN,NREC,IFILEN(NREC)                   
      GO TO 890                                                         
C                                                                       
  890 CONTINUE                                                          
      IF(MASWRK) WRITE(IW,9090)                                         
      CALL ABRT                                                         
      RETURN                                                            
C                                                                       
 9000 FORMAT(/1X,'ERROR *** ATTEMPTING A BOGUS READ OF A DAF RECORD.'/  
     *       1X,'RECORD NUMBER',I5,' OF LENGTH',I10,                    
     *          ' WAS NEVER PREVIOUSLY WRITTEN.')                       
 9010 FORMAT(/1X,'ERROR *** ATTEMPTING A BOGUS READ OF A DAF RECORD.'/  
     *       1X,'RECORD NUMBER',I5,' OF LENGTH',I10,' HAS NO LENGTH.')  
 9020 FORMAT(/1X,'ERROR *** ATTEMPTING A BOGUS READ OF A DAF RECORD.'/  
     *       1X,'ATTEMPTING TO READ',I10,' WORDS FROM RECORD NUMBER',I5/
     *       1X,'BUT THIS RECORD WAS PREVIOUSLY WRITTEN WITH ONLY',     
     *          I10,' WORDS.')                                          
 9090 FORMAT(/1X,'THIS ERROR IS LIKELY TO BE A BUG IN THE PROGRAM, BUT'/
     *       1X,'PLEASE CONSIDER INPUT ERRORS AS A POSSIBLE',           
     *          ' CAUSE, TOO.'/                                         
     *       1X,'MANUAL CHAPTER -PROG.DOC- HAS A LIST OF ALL',          
     *          ' DIRECT ACCESS FILE RECORD NUMBERS,'/                  
     *       1X,'WHICH WILL HELP YOU UNDERSTAND WHAT THE INCORRECT',    
     *          ' DATA IS SUPPOSED TO BE.'/)                            
      END                                                               
C*MODULE IOLIB   *DECK DAWRIT                                           
!> @brief Write a logical record from the dictionary file               
!> @details See DAREAD for details                                      
!>                                                                      
!> @param IDAF the unit number of the dictionary file                   
!> @param IODA                                                          
!> @param V the double precision array to be written                    
!> @param LEN size of the V array                                       
!> @param NREC                                                          
!> @param IDTYP                                                         
!> @author Unknown                                                      
!> @date September 2010                                                 
!> - Support for the CIMDCT file used in correlated subsystem calculatio
      SUBROUTINE DAWRIT(IDAF,IODA,V,LEN,NREC,IDTYP)                     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      LOGICAL NEWREC,GOPARR,DSKWRK,MASWRK,ISGDDI,PAROUT,INITGDDI,       
     *        WASGDDI,EDATYP,MLGDDI                                     
C                                                                       
      DIMENSION V(LEN),IODA(950)                                        
C                                                                       
      INTEGER DDI_WORLD,DDI_GROUP                                       
      DOUBLE PRECISION LMOEDA,NEDA                                      
      PARAMETER (DDI_WORLD=0)                                           
      PARAMETER (DDI_GROUP=1)                                           
C                                                                       
      COMMON /DAIOLN/ IRECLN,IRECST,NRECUS,IFILEN(950),MERF10(950)      
      COMMON /FFPARM/ NFFAT,NBOND,NANGL,NDIHR,NDIHB,NCMAP,NWAGG,        
     *                N1213J,N14J,NLKQMM,IDOCHG,IDOPOL,IDOLJ,IDOCMAP    
      COMMON /FMCOM / X(1)                                              
      COMMON /GDDI  / ISCOPE,NGROUPS,MYGROUP,MEGLOB,NPGLOB,NNGLOB,JBTYP,
     *                ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,NSUBGR,     
     *                MeUniv,NPUniv,numdlb,myworld,nworlds              
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAFX,NAV,IODAX(950)              
      COMMON /MACHSW/ KDIAG,ICORFL,IXDR,modio,mem10,lpnt10,mem10m       
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL                
C                                                                       
      DATA LMOEDA /8HLMOEDA  /                                          
      DATA NEDA   /8HNEDA    /                                          
C                                                                       
C         WRITE A LOGICAL RECORD ON THE DAF DICTIONARY FILE.            
C         A LOGICAL RECORD MAY SPAN SEVERAL PHYSICAL RECORDS.           
C                                                                       
      IF(IDAF.LT.0) WRITE(6,*) 'BOGUS DAWRIT, IDTYP=',IDTYP             
C                                                                       
C     INITILIZE ON ALL MASTERS EVERYTHING THAT HAS BEEN INITIALIZED     
C     ON THE GLOBAL MASTER (IF -INITGDDI-).  IN OTHER WORDS, SWITCH     
C     TO THE GROUP SCOPE, AND AFTER THE WRITE, SWITCH BACK TO WORLD.    
C                                                                       
      IF(ISGDDI.AND.INITGDDI) CALL GDDI_SCOPE(DDI_GROUP)                
C                                                                       
c     if(nrec.eq.14) write(6,*) 'FFsav',v                               
c     write(6,*) 'write',NREC                                           
      IF(IDAF.NE.10) THEN                                               
         if(idaf.ne.194) then                                           
            CALL RAWRIT(IDAF,IODA,V,LEN,NREC,NAV)                       
            IF(ISGDDI.AND.INITGDDI) CALL GDDI_SCOPE(DDI_WORLD)          
            RETURN                                                      
         endif                                                          
      ENDIF                                                             
C                                                                       
C        NEARLY ALWAYS, ENFORCE CONDITION THAT WRITING SHOULD           
C        USE THE SAME RECORD SIZE AS THE FIRST TIME WRITTEN.            
C        ENERGY DECOMPOSITION ANALYSIS MAY WORK WITH SMALLER            
C        MONOMERS, AND MUST FIRST WRITE FULL SIZED RECORDS,             
C        AFTER WHICH SMALLER WRITES ARE PERMISSIBLE.                    
C        NOTE: -NBO- RUNS MAY TEMPORARILY SWITCH THE RUN TYPE           
C        TO -NEDA- WHICH IS OTHERWISE NOT SEEN (E.G. IN INPUT).         
C                                                                       
      N = IODA(NREC)                                                    
C                                                                       
      EDATYP = RUNTYP.EQ.LMOEDA  .OR.  RUNTYP.EQ.NEDA  .OR.             
     *         NFFAT .GT.0                                              
C                                                                       
      IF(.NOT.EDATYP) THEN                                              
         IF (N .GT. 0 .AND. LEN .NE. IFILEN(NREC)) GO TO 800            
      ELSE                                                              
         IF (N .GT. 0 .AND. LEN .GT. IFILEN(NREC)) GO TO 800            
      END IF                                                            
c                                                                       
      if(mem10.ne.0) then                                               
c       merf10(nrec)=mod(nrec,nproc)                                    
        if(merf10(nrec).lt.0) then                                      
          nrecus=nrecus+1                                               
          merf10(nrec)=mod(nrecus,nproc)                                
          if(n.gt.0) call abrt                                          
c         sanity check: IODA and MERF10 should be bot set or unset.     
        endif                                                           
        men=merf10(nrec)                                                
        if(men.eq.me) then                                              
c         Divide storage by record number. A better way is to pick the  
c         core which used less memory!                                  
c         Each core maintains its very own record/memory bookkeeping!   
          if(n.le.0) then                                               
            IFILEN(NREC)=LEN                                            
            IODA(NREC)=IRECST                                           
            IRECST=IRECST+LEN                                           
            if(IRECST.gt.mem10) then                                    
c             Now, slaves toil too, so a slave can run out of memory.   
              write(iw,*) 'Out of memory for F10:',ME,IRECST,mem10      
              call abrt                                                 
            else                                                        
              mem10m=max(mem10m,IRECST-1)                               
            endif                                                       
c           write(6,*) 'New F10 record',me,IRECST,IODA(NREC),LEN        
          endif                                                         
          call XCOPY(len,v,1,x(lpnt10+IODA(NREC)-1),1)                  
c         write(6,*) 'Write F10 record',me,NREC,IODA(NREC)              
        else                                                            
c         chicken out and save the record size and a fake pointer on all
c         cores, just to do error processing in DAREAD on all cores.    
c         One could avoid the data saving below but then DAREAD needs a 
          if(n.le.0) then                                               
            IFILEN(NREC)=LEN                                            
            IODA(NREC)=IRECST                                           
c           DO NOT change IRECST here! IRECST and IODA are core specific
c           pointing to what each core actually saved. IFILEN is duplica
          endif                                                         
        endif                                                           
c       We do not update record 1.                                      
        IF(ISGDDI.AND.INITGDDI) CALL GDDI_SCOPE(DDI_WORLD)              
c       Do we really need this scope change?                            
        return                                                          
      endif                                                             
C                                                                       
      NEWREC = .FALSE.                                                  
      IF (N .GT. 0) GO TO 100                                           
      IODA(NREC) = IRECST                                               
      IFILEN(NREC) = LEN                                                
      NEWREC = .TRUE.                                                   
      IRECST = IRECST + (LEN-1)/IRECLN + 1                              
      N = IODA(NREC)                                                    
  100 CONTINUE                                                          
      IST = -IRECLN + 1                                                 
      NS = N                                                            
      LENT = LEN                                                        
  120 CONTINUE                                                          
         IST = IST + IRECLN                                             
         IF = IST + LENT - 1                                            
         IF ((IF-IST+1) .GT. IRECLN) IF = IST+IRECLN-1                  
         NSP = NS                                                       
         LENW = IF - IST + 1                                            
         CALL DAWRT(V(IST),LENW,IDAF,NSP)                               
         LENT = LENT - IRECLN                                           
         NS = NS + 1                                                    
         N = NS                                                         
      IF (LENT .GE. 1) GO TO 120                                        
      IF (NEWREC .AND. MASWRK)                                          
     *      WRITE(UNIT=IDAF,REC=1) IRECST,IODA,IFILEN,IS,IPK            
      IF(ISGDDI.AND.INITGDDI) CALL GDDI_SCOPE(DDI_WORLD)                
      RETURN                                                            
C                                                                       
  800 CONTINUE                                                          
      WRITE (IW,9008) ME,NREC,LEN,IFILEN(NREC)                          
      CALL ABRT                                                         
      RETURN                                                            
C                                                                       
 9008 FORMAT(1X,'*** ERROR *** IN -DAWRIT- ROUTINE ON NODE',I4/         
     *       1X,'DAWRIT HAS REQUESTED A RECORD WITH LENGTH',            
     *       1X,'DIFFERENT THAN BEFORE - ABORT FORCED.'/                
     *       1X,'DAF RECORD ',I5,' NEW LENGTH =',I10,                   
     *                           ' OLD LENGTH =',I10)                   
      END                                                               
C*MODULE IOLIB   *DECK DAWRT                                            
      SUBROUTINE DAWRT(V,LEN,IDAF,NS)                                   
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      LOGICAL GOPARR,DSKWRK,MASWRK                                      
      DIMENSION V(LEN)                                                  
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
C                                                                       
C     ----- WRITE A PHYSICAL RECORD ON THE DAF -----                    
C                                                                       
      IF (MASWRK) WRITE (UNIT=IDAF, REC=NS) V                           
      RETURN                                                            
      END                                                               
C*MODULE IOLIB   *DECK OPENDA                                           
      SUBROUTINE OPENDA(IREST)                                          
C                                                                       
C     - - - - OPEN MASTER DICTIONARY FILE 10 - - - -                    
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      LOGICAL GOPARR,DSKWRK,MASWRK,ISGDDI,wasgddi,PAROUT,INITGDDI,MLGDDI
C                                                                       
      PARAMETER (MXUNIT=299)                                            
      CHARACTER*256 FILENM,ENVBUF                                       
      COMMON /ENVIR / ENVBUF(-5:MXUNIT)                                 
C                                                                       
      COMMON /DAIOLN/ IRECLN,IRECST,NRECUS,IFILEN(950),MERF10(950)      
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)                
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
      COMMON /GDDI  / ISCOPE,NGROUPS,MYGROUP,MEGLOB,NPGLOB,NNGLOB,JBTYP,
     *                ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,NSUBGR,     
     *                MeUniv,NPUniv,numdlb,myworld,nworlds              
      COMMON /MACHSW/ KDIAG,ICORFL,IXDR,modio,mem10,lpnt10,mem10m       
      COMMON /OPNNFT/ NFTOPN(MXUNIT),NODEXT(MXUNIT),IOSMP(2)            
C                                                                       
      IDAF = 10                                                         
      IRECLN=NRASIZ(IDAF)                                               
C                                                                       
      NRECUS=0                                                          
c     NRECUS is the number of records actually saved to DAF.            
c                                                                       
      IF ((MASWRK.or.mem10.ne.0).and.NFTOPN(idaf).EQ.0) THEN            
      NFTOPN(IDAF) =1                                                   
C                                                                       
      if(mem10.eq.0) then                                               
*IBM  OPEN (UNIT=IDAF, FILE='DICTNRY', STATUS='UNKNOWN',                
*IBM *      ACCESS='DIRECT', FORM='UNFORMATTED', RECL=8*IRECLN)         
C                                                                       
      IF(ISGDDI.or.wasgddi) THEN                                        
         FILENM=ENVBUF(IDAF)                                            
         IF(NODEXT(IDAF).EQ.0)                                          
     *       CALL ADDNANODE(FILENM,MeUniv,IDAF,iout)                    
      ELSE                                                              
         CALL GMS_GETENV('DICTNRY',FILENM)                              
      ENDIF                                                             
      OPEN (UNIT=IDAF, FILE=FILENM, STATUS='UNKNOWN',                   
     *      ACCESS='DIRECT', FORM='UNFORMATTED',                        
     *      RECL=8*IRECLN)                                              
C                                                                       
*VMS  OPEN (UNIT=IDAF, FILE='DICTNRY', STATUS='UNKNOWN',                
*VMS *      ACCESS='DIRECT', FORM='UNFORMATTED', RECL=2*IRECLN)         
      else if(lpnt10.eq.0) then                                         
         CALL VALFM(LOADFM)                                             
         lpnt10= LOADFM + 1                                             
         LAST  = lpnt10 + mem10                                         
         NEED  = LAST- LOADFM -1                                        
         CALL GETFM(NEED)                                               
c        We do not return memory here. It is important to return it     
c        properly later (it is done in PROGRAM GAMESS).                 
         if(maswrk) write(6,*) 'Allocated',mem10,' for F10'             
      endif                                                             
C                                                                       
      END IF                                                            
C                                                                       
C     ----- IS THIS A NEW OR OLD DAF FILE? -----                        
C     EITHER MARK THE NEW DAF RECORDS AS EMPTY                          
C     OR ELSE LOAD THE OLD DAF DIRECTORY                                
C                                                                       
      IF(IREST.EQ.0) THEN                                               
         IRECST = 1                                                     
         DO 100 I = 1,950                                               
            IODA(I) = -1                                                
            MERF10(I) = -1                                              
  100    CONTINUE                                                       
         if(mem10.eq.0) then                                            
           IRECST = IRECST + 1                                          
           if(maswrk) WRITE(UNIT=IDAF, REC=1) IRECST,IODA,IFILEN,IS,IPK 
         else                                                           
c          Do nothing: no restarts!                                     
         endif                                                          
      ELSE                                                              
         if(mem10.ne.0) call abrt                                       
c        Restarts are not supported with in-memory F10!                 
         IF(MASWRK) READ (UNIT=IDAF, REC=1) IRECST,IODA,IFILEN,IS,IPK   
         IF(GOPARR) THEN                                                
            CALL DDI_BCAST(200,'I',IRECST,1,MASTER)                     
            CALL DDI_BCAST(201,'I',IODA,950,MASTER)                     
            CALL DDI_BCAST(202,'I',IFILEN,950,MASTER)                   
            CALL DDI_BCAST(203,'I',IS,1,MASTER)                         
            CALL DDI_BCAST(204,'I',IPK,1,MASTER)                        
         END IF                                                         
      END IF                                                            
      RETURN                                                            
      END                                                               
C*MODULE IOLIB   *DECK MQDARE                                           
      SUBROUTINE MQDARE(IDAF,IODA,V,LEN,NREC,IDTYP)                     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      LOGICAL GOPARR,DSKWRK,MASWRK                                      
C                                                                       
      DIMENSION V(LEN),IODA(400)                                        
C                                                                       
      COMMON /IOFILE/ IR,IW,IP,IIS,IPK,IDAFX,NAV,IODAX(950)             
      COMMON /MQDAIO/ IRECLN,IRECST,IFILEN(400)                         
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
C                                                                       
C         READ A LOGICAL RECORD FROM THE DAF DICTIONARY FILE            
C         A LOGICAL RECORD MAY SPAN SEVERAL PHYSICAL RECORDS.           
C                                                                       
C         CALLING ARGUMENT -IDTYP- IS 0 OR 1 FOR FLOATING POINT         
C         OR INTEGER RECORDS.  NO OTHER DATA TYPE IS ALLOWED!           
C         RECORDS MUST BE PURELY FLOATING POINT OR PURELY INTEGER.      
C         NO MATTER WHAT -IDTYP- IS, THE -LEN- OF THE RECORD            
C         MUST BE GIVEN IN TERMS OF FLOATING POINT WORDS.  IN           
C         TURN, THIS MEANS THAT INTEGER RECORDS ON 32 BIT MACHINES      
C         MUST CONTAIN AN EVEN NUMBER OF INTEGERS.                      
C                                                                       
      N = IODA(NREC)                                                    
      IF(N.EQ.-1) GO TO 800                                             
      IS = -IRECLN + 1                                                  
      NS = N                                                            
      LENT = LEN                                                        
  100 CONTINUE                                                          
         IS = IS + IRECLN                                               
         IF = IS + LENT - 1                                             
         IF ((IF-IS+1) .GT. IRECLN) IF = IS + IRECLN - 1                
         NSP = NS                                                       
         LENW = IF - IS + 1                                             
         CALL DARD(V(IS),LENW,IDAF,NSP,IDTYP)                           
         LENT = LENT - IRECLN                                           
         NS = NS + 1                                                    
         N = NS                                                         
      IF (LENT .GE. 1) GO TO 100                                        
      RETURN                                                            
C                                                                       
  800 CONTINUE                                                          
      IF (MASWRK) WRITE(IW,9000) NREC,LEN                               
      CALL ABRT                                                         
      RETURN                                                            
C                                                                       
 9000 FORMAT(1X,'*** ERROR ***, ATTEMPT TO READ A -MCQD50- RECORD',     
     *         ' THAT WAS NEVER WRITTEN.'/1X,'NREC,LEN=',I5,I10)        
      END                                                               
C*MODULE IOLIB   *DECK MQDAWR                                           
      SUBROUTINE MQDAWR(IDAF,IODA,V,LEN,NREC,IDTYP)                     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      LOGICAL NEWREC,GOPARR,DSKWRK,MASWRK                               
C                                                                       
      DIMENSION V(LEN),IODA(400)                                        
C                                                                       
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAFX,NAV,IODAX(950)              
      COMMON /MQDAIO/ IRECLN,IRECST,IFILEN(400)                         
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
C                                                                       
C         WRITE A LOGICAL RECORD ON THE DAF DICTIONARY FILE             
C         A LOGICAL RECORD MAY SPAN SEVERAL PHYSICAL RECORDS            
C                                                                       
      IF(IDAF.LT.0) WRITE(6,*) 'BOGUS DAWRIT, IDTYP=',IDTYP             
C                                                                       
      N = IODA(NREC)                                                    
      IF (N .GT. 0 .AND. LEN .NE. IFILEN(NREC)) GO TO 800               
      NEWREC = .FALSE.                                                  
      IF (N .GT. 0) GO TO 100                                           
      IODA(NREC) = IRECST                                               
      IFILEN(NREC) = LEN                                                
      NEWREC = .TRUE.                                                   
      IRECST = IRECST + (LEN-1)/IRECLN + 1                              
      N = IODA(NREC)                                                    
  100 CONTINUE                                                          
      IST = -IRECLN + 1                                                 
      NS = N                                                            
      LENT = LEN                                                        
  120 CONTINUE                                                          
         IST = IST + IRECLN                                             
         IF = IST + LENT - 1                                            
         IF ((IF-IST+1) .GT. IRECLN) IF = IST+IRECLN-1                  
         NSP = NS                                                       
         LENW = IF - IST + 1                                            
         CALL DAWRT(V(IST),LENW,IDAF,NSP)                               
         LENT = LENT - IRECLN                                           
         NS = NS + 1                                                    
         N = NS                                                         
      IF (LENT .GE. 1) GO TO 120                                        
      IF (NEWREC .AND. MASWRK)                                          
     *      WRITE(UNIT=IDAF,REC=1) IRECST,IODA,IFILEN,IS,IPK            
      RETURN                                                            
C                                                                       
  800 CONTINUE                                                          
      IF (MASWRK) WRITE (IW,9008) NREC,LEN,IFILEN(NREC)                 
      CALL ABRT                                                         
      RETURN                                                            
C                                                                       
 9008 FORMAT(1X,'MQDAWR HAS REQUESTED A RECORD WITH LENGTH',            
     *       1X,'DIFFERENT THAN BEFORE - ABORT FORCED.'/                
     *       1X,'MCQD50 RECORD ',I5,' NEW LENGTH =',I5,                 
     *          ' OLD LENGTH =',I5)                                     
      END                                                               
C*MODULE IOLIB   *DECK MQOPDA                                           
      SUBROUTINE MQOPDA(IREST)                                          
C=====================================================                  
C     A VARIANT OF OPENDA CREATED BY H. NAKANO 8/14/96                  
C=====================================================                  
C                                                                       
C     - - - - OPEN RANDOM ACCESS FILE - - - -                           
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      LOGICAL GOPARR,DSKWRK,MASWRK                                      
C                                                                       
      CHARACTER*256 FILENM                                              
C                                                                       
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)                
      COMMON /MQDAIO/ IRECLN,IRECST,IFILEN(400)                         
      COMMON /MQIOFI/ IDAF50,NAV50,IODA50(400)                          
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
C                                                                       
*IBM  IRECLN = 512                                                      
      IRECLN = 1023                                                     
*VMS  IRECLN = 512                                                      
C                                                                       
      IF (MASWRK) THEN                                                  
C                                                                       
*IBM  OPEN (UNIT=IDAF50, FILE='MCQD50', STATUS='UNKNOWN',               
*IBM *      ACCESS='DIRECT', FORM='UNFORMATTED', RECL=8*IRECLN)         
C                                                                       
      CALL GMS_GETENV('MCQD50',FILENM)                                  
      OPEN (UNIT=IDAF50, FILE=FILENM, STATUS='UNKNOWN',                 
     *      ACCESS='DIRECT', FORM='UNFORMATTED',                        
     *      RECL=8*IRECLN)                                              
C                                                                       
*VMS  OPEN (UNIT=IDAF50, FILE='MCQD50', STATUS='UNKNOWN',               
*VMS *      ACCESS='DIRECT', FORM='UNFORMATTED', RECL=2*IRECLN)         
C                                                                       
      END IF                                                            
C                                                                       
C     ----- IS THIS A NEW OR OLD DAF FILE? -----                        
C     EITHER MARK THE NEW DAF RECORDS AS EMPTY                          
C     OR ELSE LOAD THE OLD DAF DIRECTORY                                
C                                                                       
      IF(IREST.EQ.0) THEN                                               
         IRECST = 1                                                     
         DO 100 I = 1,400                                               
            IODA50(I) = -1                                              
  100    CONTINUE                                                       
         IRECST = IRECST + 1                                            
         IF(MASWRK) WRITE(UNIT=IDAF50,REC=1) IRECST,IODA50,IFILEN,IS,IPK
      ELSE                                                              
         IF(MASWRK) READ (UNIT=IDAF50,REC=1) IRECST,IODA50,IFILEN,IS,IPK
         IF(GOPARR) THEN                                                
            CALL DDI_BCAST(200,'I',IRECST,1,MASTER)                     
            CALL DDI_BCAST(201,'I',IODA50,400,MASTER)                   
            CALL DDI_BCAST(202,'I',IFILEN,400,MASTER)                   
            CALL DDI_BCAST(203,'I',IS,1,MASTER)                         
            CALL DDI_BCAST(204,'I',IPK,1,MASTER)                        
         END IF                                                         
      END IF                                                            
      RETURN                                                            
      END                                                               
C*MODULE IOLIB   *DECK PARENV                                           
C>                                                                      
C>     @brief environmental variables                                   
C>                                                                      
C>     @details Process environmental variables for parallel runs.      
C>                                                                      
C>     @author unknown                                                  
C>                                                                      
      SUBROUTINE PARENV(FNAME,FILENM,IOUT)                              
C                                                                       
      IMPLICIT INTEGER (A-Z)                                            
C                                                                       
      PARAMETER (MAXLEN=256)                                            
      PARAMETER (MXUNIT=299)                                            
C                                                                       
      LOGICAL GOPARR,DSKWRK,MASWRK,ISGDDI,PAROUT,INITGDDI,wasgddi,      
     *        MLGDDI                                                    
C                                                                       
      INTEGER CFILENM(MAXLEN)                                           
C                                                                       
      CHARACTER*(*) FNAME                                               
      CHARACTER*256 FILENM                                              
      CHARACTER*256 DUMMY                                               
C                                                                       
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
      COMMON /GDDI/   ISCOPE,NGROUPS,MYGROUP,MEGLOB,NPGLOB,NNGLOB,JBTYP,
     *                ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,NSUBGR,     
     *                MeUniv,NPUniv,numdlb,myworld,nworlds              
      COMMON /OPNNFT/ NFTOPN(MXUNIT),NODEXT(MXUNIT),IOSMP(2)            
C                                                                       
C     THIS ROUTINE DECIDES WHETHER OR NOT SLAVES SHOULD OPEN A FILE.    
C     IF A FILE IS TO BE OPENED, THE SLAVES GET THE NAME OF THE FILE    
C     TO OPEN.                                                          
C                                                                       
      IOUT = 0                                                          
C                                                                       
      IF ((FNAME.EQ.'INPUT')   .OR. (FNAME.EQ.'OUTPUT')  .OR.           
     *    (FNAME.EQ.'PUNCH')   .OR. (FNAME.EQ.'RESTART') .OR.           
     *    (FNAME.EQ.'TRAJECT') .OR. (FNAME.EQ.'MAKEFP')  .OR.           
     *       (.NOT.DSKWRK)     .OR.     (NPROC.EQ.1)) THEN              
         IF (.NOT.MASWRK) IOUT = 1                                      
C                                                                       
C        THIS IS USEFUL FOR DEBUGGING. SLAVE NODES ARE NOT SUPPOSED     
C        TO NORMALLY WRITE SOMETHING OUT TO IW. AS USUAL, WRITING TO    
C        THE SAME FILE BY SEVERAL NODES CAN LEAD TO FILE CONTENT        
C        CORRUPTION.                                                    
C                                                                       
         IF(ISGDDI.AND.FNAME.EQ.'OUTPUT') IOUT = 0                      
         RETURN                                                         
      END IF                                                            
C                                                                       
C     IF FILENM IS EMPTY, ABRT ALL SLAVES; ELSE STRIP OFF THE PATH      
C     AND SEND THE NAME TO ALL SLAVES.                                  
C                                                                       
      IF (MASWRK) THEN                                                  
C                                                                       
         IEND = 0                                                       
         IF (FILENM(1:1).EQ.' ') THEN                                   
            IEND = 1                                                    
            CALL DDI_BCAST(205,'I',IEND,1,MASTER)                       
            RETURN                                                      
         END IF                                                         
C                                                                       
         DUMMY = FILENM                                                 
         INDEX = 0                                                      
C                                                                       
C             NEXT CODE STRIPS PATH NAME IF ISTRIP IS 1                 
C             TO PARALLELISE I/O OVER SEVERAL HDDS, ISTRIP MUST BE 0.   
C             A 'SED' HACK CAN BE USED TO PREVENT STRIPPING.            
C                                                                       
         ISTRIP=1                                                       
         IF(IOSMP(2).NE.0) ISTRIP=0                                     
         IF(ISTRIP.EQ.1) THEN                                           
            DO 10 I = MAXLEN,1,-1                                       
               IF (FILENM(I:I).EQ.'/') THEN                             
                  INDEX = I                                             
                  GO TO 20                                              
               END IF                                                   
   10       CONTINUE                                                    
   20       CONTINUE                                                    
C                                                                       
            IF (INDEX.NE.0) THEN                                        
               LDIFF = MAXLEN - INDEX                                   
               DO 30 I = 1,LDIFF                                        
                  INDEX = INDEX + 1                                     
                  FILENM(I:I) = FILENM(INDEX:INDEX)                     
   30          CONTINUE                                                 
            END IF                                                      
         END IF                                                         
C                                                                       
         LWORD = 0                                                      
         DO 40 I = 1,MAXLEN                                             
            IF ((FILENM(I:I).EQ.CHAR(0)).OR.(FILENM(I:I).EQ.' '))       
     *        GO TO 50                                                  
            LWORD = LWORD + 1                                           
   40    CONTINUE                                                       
   50    CONTINUE                                                       
      END IF                                                            
C                                                                       
      CALL DDI_BCAST(205,'I',IEND,1,MASTER)                             
      IF ((.NOT.MASWRK).AND.(IEND.EQ.1)) CALL ABRT                      
C                                                                       
      IF (MASWRK) THEN                                                  
         DO 100 I=1,MAXLEN                                              
            CFILENM(I) = ICHAR(FILENM(I:I))                             
  100    CONTINUE                                                       
      END IF                                                            
C                                                                       
      CALL DDI_BCAST(206,'I',CFILENM,MAXLEN,MASTER)                     
      CALL DDI_BCAST(207,'I',LWORD,1,MASTER)                            
C                                                                       
C     MASTER CAN RETURN NOW.                                            
C                                                                       
      IF (MASWRK) THEN                                                  
         FILENM = DUMMY                                                 
         RETURN                                                         
      END IF                                                            
C                                                                       
      DO 150 I=1,MAXLEN                                                 
         FILENM(I:I) = CHAR(CFILENM(I))                                 
  150 CONTINUE                                                          
C                                                                       
C     APPEND NODE NUMBER, IN ORDER TO OPEN UNIQUE FILE NAMES.           
C                                                                       
      DO 175 I = 1,LWORD                                                
         DUMMY(I:I) = FILENM(I:I)                                       
  175 CONTINUE                                                          
C                                                                       
                          NDIGIT=3                                      
      IF(NPROC.GE.1000)   NDIGIT=4                                      
      IF(NPROC.GE.10000)  NDIGIT=5                                      
      IF(NPROC.GE.100000) NDIGIT=6                                      
C                                                                       
      IF(NDIGIT.EQ.3) THEN                                              
         WRITE(UNIT=DUMMY(LWORD+1:LWORD+4),FMT='(1H.,I3.3)') ME         
         DUMMY(LWORD+5:LWORD+5) = CHAR(0)                               
      ELSE IF(NDIGIT.EQ.4) THEN                                         
         WRITE(UNIT=DUMMY(LWORD+1:LWORD+5),FMT='(1H.,I4.4)') ME         
         DUMMY(LWORD+6:LWORD+6) = CHAR(0)                               
      ELSE IF(NDIGIT.EQ.5) THEN                                         
         WRITE(UNIT=DUMMY(LWORD+1:LWORD+6),FMT='(1H.,I5.5)') ME         
         DUMMY(LWORD+7:LWORD+7) = CHAR(0)                               
      ELSE                                                              
         IF(MASWRK) WRITE(6,*) 'PARENV: TOO MANY CPUS REQUESTED'        
         CALL ABRT                                                      
      END IF                                                            
C                                                                       
      FILENM = DUMMY                                                    
C                                                                       
C     PARALLELISE I/O ON SMPS!                                          
C     SOME (INPUT) FILES ARE EXCEPTED FROM THIS: 2, 3 AND IR.           
C     ISTRIP=1 -> ISTRIP=0 SED HACK IS REQUIRED!!                       
C     MASTER IS NOT AFFECTED.                                           
C                                                                       
      J=IOSMP(2)                                                        
      IF(J.NE.0.AND.J.LE.MAXLEN.AND.ISTRIP.EQ.0) THEN                   
C       WRITE(6,*) ME,FILENM(1:20)                                      
        if(iosmp(1).ge.0) then                                          
          M=MOD(MEGLOB,IOSMP(1))                                        
        else                                                            
          M=MOD(-MEGLOB/IOSMP(1),2)                                     
        endif                                                           
        READ(UNIT=FILENM(J:J),FMT='(I1)') ID                            
        WRITE(UNIT=FILENM(J:J),FMT='(I1)') ID+M                         
C       WRITE(6,*) ME,FILENM(1:20)                                      
      ENDIF                                                             
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE IOLIB   *DECK PKREAD                                           
      SUBROUTINE PKREAD(IS,XP,XK,IX,NX,MAX)                             
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      LOGICAL GOPARR,DSKWRK,MASWRK                                      
      DIMENSION XP(*),XK(*),IX(*)                                       
      COMMON /MACHIN/ NWDVAR,MAXFM,MAXSM,LIMFM,LIMSM                    
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
      COMMON /PCKLAB/ LABSIZ                                            
C                                                                       
      IF (LABSIZ.EQ.1) THEN                                             
         IF(NWDVAR.EQ.2) MAX2 = MAX                                     
         IF(NWDVAR.EQ.1) MAX2 = (MAX+1)/2                               
      ELSE IF (LABSIZ.EQ.2) THEN                                        
         IF(NWDVAR.EQ.2) MAX2 = MAX*2                                   
         IF(NWDVAR.EQ.1) MAX2 = MAX                                     
      ELSE                                                              
         IF(MASWRK) WRITE(6,*) 'PKREAD: CONFUSION WITH -LABSIZ-'        
         CALL ABRT                                                      
      END IF                                                            
      CALL PKRD(IS,XP,XK,IX,NX, MAX, MAX2 )                             
      RETURN                                                            
      END                                                               
C*MODULE IOLIB   *DECK PKRD                                             
      SUBROUTINE PKRD(IS,XP,XK,IX,NX,MAX,MAX2)                          
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      LOGICAL GOPARR,DSKWRK,MASWRK                                      
      DIMENSION XP(MAX),XK(MAX),IX(MAX2)                                
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
C                                                                       
      IF (DSKWRK.OR.MASWRK) READ(IS,END=200,ERR=300) NX,IX,XP,XK        
      RETURN                                                            
C                                                                       
  200 CONTINUE                                                          
      WRITE(6,9010) ME,IS                                               
 9010 FORMAT(1X,'PKRD: NODE',I4,                                        
     *          ' ENCOUNTERED UNEXPECTED END OF FILE READING UNIT',I4)  
      CALL ABRT                                                         
C                                                                       
  300 CONTINUE                                                          
      WRITE(6,9020) ME,IS                                               
 9020 FORMAT(1X,'PKRD: NODE',I4,                                        
     *          ' ENCOUNTERED I/O ERROR READING UNIT',I4)               
      CALL ABRT                                                         
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE IOLIB   *DECK PREAD                                            
      SUBROUTINE PREAD(IS,XX,IX,NX,MAX)                                 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      LOGICAL GOPARR,DSKWRK,MASWRK                                      
      DIMENSION XX(MAX),IX(*)                                           
      COMMON /INT2IC/ NINTIC,ININTIC,NXXIC,LBUFPIC,LIXIC,LABSIX,NINTIX  
      COMMON /MACHIN/ NWDVAR,MAXFM,MAXSM,LIMFM,LIMSM                    
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
      COMMON /PCKLAB/ LABSIZ                                            
C                                                                       
      IF (LABSIZ.EQ.1) THEN                                             
         IF(NWDVAR.EQ.2) MAX2 = MAX                                     
         IF(NWDVAR.EQ.1) MAX2 = (MAX+1)/2                               
      ELSE IF (LABSIZ.EQ.2) THEN                                        
         IF(NWDVAR.EQ.2) MAX2 = MAX*2                                   
         IF(NWDVAR.EQ.1) MAX2 = MAX                                     
      ELSE                                                              
         IF(MASWRK) WRITE(6,*) 'PREAD: CONFUSION WITH -LABSIZ-'         
         CALL ABRT                                                      
      END IF                                                            
C                                                                       
      CALL PRD(IS,XX,IX,NX, MAX, MAX2 )                                 
C                                                                       
C       MOST LIKELY, THE CORRESPONDING CODE HAS NOT BEEN CHANGED.       
C       IT IS ALSO POSSIBLE THAT THE NUMBER OF INTEGRALS EXACTLY MATCHES
C       THE NINTIC CAPACITY, THUS ONE EMPTY PWRITE IS DONE TO INDICATE  
C       NO MORE INTEGRALS.                                              
C                                                                       
      IF(NX.EQ.0 .AND. IS.EQ.8 .AND. NINTIC.NE.0 .AND. MASWRK) THEN     
         WRITE(6,*) 'WARNING: NO INTEGRALS TO READ: NINTIC BUG?'        
      END IF                                                            
      RETURN                                                            
      END                                                               
C*MODULE IOLIB   *DECK PRD                                              
      SUBROUTINE PRD(IS,XX,IX,NX,MAX,MAX2)                              
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      LOGICAL GOPARR,DSKWRK,MASWRK                                      
      DIMENSION XX(MAX),IX(MAX2)                                        
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
C                                                                       
      IF (DSKWRK.OR.MASWRK) READ(IS,END=200,ERR=300) NX,IX,XX           
      RETURN                                                            
C                                                                       
  200 CONTINUE                                                          
      WRITE(6,9010) ME,IS                                               
 9010 FORMAT(1X,'PRD: NODE',I4,                                         
     *          ' ENCOUNTERED UNEXPECTED END OF FILE READING UNIT',I4)  
      CALL ABRT                                                         
C                                                                       
  300 CONTINUE                                                          
      WRITE(6,9020) ME,IS                                               
 9020 FORMAT(1X,'PRD: NODE',I4,                                         
     *          ' ENCOUNTERED I/O ERROR READING UNIT',I4)               
      CALL ABRT                                                         
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE IOLIB   *DECK PREADP                                           
      SUBROUTINE PREADP(NFT,XX,IX,NXX,NINTMX,KAP,IFLG,JFLG,NPROC)       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      LOGICAL GOPARR,DSKWRK,MASWRK                                      
      DIMENSION XX(NINTMX),IX(NINTMX),JFLG(0:NPROC-1)                   
      COMMON /MACHIN/ NWDVAR,MAXFM,MAXSM,LIMFM,LIMSM                    
      COMMON /PAR   / ME,MASTER,NPROCX,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK 
      COMMON /PCKLAB/ LABSIZ                                            
C                                                                       
C     ----- PARALLEL FORM OF -PREAD- ROUTINE -----                      
C     THIS ROUTINE READS A SINGLE RECORD, FROM WHICHEVER NODE WHOSE     
C     TURN IS NEXT, AND BROADCASTS IT TO ALL NODES.  THIS MAKES A       
C     FILE WHICH HAS BEEN DISTRIBUTED ACROSS NODES AVAILABLE IN ITS     
C     ENTIRETY TO ALL NODES.                                            
C                                                                       
C     THIS ROUTINE READS IN A ROUND-ROBIN FASHION, WITH NODE -KAP-      
C     READING THE 1ST RECORD FROM NODES 0,1,...,NPROC-1, THEN THE 2ND   
C     RECORDS FROM EACH, ETC, ETC.  -JFLG- KEEPS TRACK OF END-OF-FILE   
C     ON EACH NODE, WHEREAS -IFLG- KEEPS TRACK OF END-OF-FILE ON ALL    
C     NODES, IFLG.EQ.0 ON EXIT MEANS ALL NODES HAVE PREVIOUSLY HIT      
C     END OF FILE, AND THUS NO DATA RECORD IS RETURNED.                 
C                                                                       
C                                                                       
  100 CONTINUE                                                          
      KAP = KAP+1                                                       
      IF(KAP.EQ.NPROC) THEN                                             
         IF(IFLG.EQ.0) RETURN                                           
         KAP = 0                                                        
         IFLG= 0                                                        
      END IF                                                            
C                                                                       
      IF(JFLG(KAP).GT.0) THEN                                           
         IFLG=1                                                         
         IF(KAP.EQ.ME) CALL PREAD(NFT,XX,IX,NXX,NINTMX)                 
         CALL DDI_BCAST(1205,'I',NXX,1,KAP)                             
         NINT = ABS(NXX)                                                
         IF (LABSIZ.EQ.1) THEN                                          
            IF(NWDVAR.EQ.2) NINT2 = NINT                                
            IF(NWDVAR.EQ.1) NINT2 = (NINT+1)/2                          
         ELSE IF (LABSIZ.EQ.2) THEN                                     
            IF(NWDVAR.EQ.2) NINT2 = NINT*2                              
            IF(NWDVAR.EQ.1) NINT2 = NINT                                
         ELSE                                                           
            IF(MASWRK) WRITE(6,*) 'PREADP: CONFUSION WITH -LABSIZ-'     
            CALL ABRT                                                   
         END IF                                                         
         CALL DDI_BCAST(1210,'F',XX,NINT,KAP)                           
         CALL DDI_BCAST(1215,'I',IX,NINT2,KAP)                          
         JFLG(KAP) = NXX                                                
         RETURN                                                         
      END IF                                                            
C                                                                       
      GO TO 100                                                         
      END                                                               
C*MODULE IOLIB   *DECK PXREAD                                           
      SUBROUTINE PXREAD(NFT,X,XX,MX,NX)                                 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      LOGICAL GOPARR,DSKWRK,MASWRK                                      
      DIMENSION X(NX),XX(NX)                                            
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
C                                                                       
      IF (DSKWRK.OR.MASWRK) READ(NFT,END=200,ERR=300) X,XX,MX           
      RETURN                                                            
C                                                                       
  200 CONTINUE                                                          
      WRITE(6,9010) ME,NFT                                              
 9010 FORMAT(1X,'PXREAD: NODE',I4,                                      
     *          ' ENCOUNTERED UNEXPECTED END OF FILE READING UNIT',I4)  
      CALL ABRT                                                         
C                                                                       
  300 CONTINUE                                                          
      WRITE(6,9020) ME,NFT                                              
 9020 FORMAT(1X,'PXREAD: NODE',I4,                                      
     *          ' ENCOUNTERED I/O ERROR READING UNIT',I4)               
      CALL ABRT                                                         
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE IOLIB   *DECK PUVIB                                            
      SUBROUTINE PUVIB(IFCM,IW,RSTART,NCOORD,IVIB,IATOM,ICOORD,         
     *                 E,EG,DIP)                                        
C                                                                       
      USE MX_LIMITS,ONLY:MXFRG,MXFGPT,MXDFG,MXDPPT,mxatm                
      USE comm_FGRAD                                                    
      USE comm_FRGINF                                                   
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C                                                                       
C                                                                       
      DIMENSION EG(NCOORD), DIP(3), BLANK(1), VIBWRD(1)                 
C                                                                       
      LOGICAL RSTART,GOPARR,DSKWRK,MASWRK,FGONLY                        
C                                                                       
      COMMON /FMOINF/ NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY              
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,                    
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)                  
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
      COMMON /SIMDAT/ NACC,NREJ,IGOMIN,NRPA,IBWM,NACCT,NREJT,NRPAT,     
     *                NPRTGO,IDPUNC,IGOFLG                              
C                                                                       
      DATA BLANK/8H        /, VIBWRD/8H $VIB   /                        
C                                                                       
C     ----- SAVE DATA FOR CURRENT VIBRATION FOR HESSIAN RESTARTS -----  
C     THIS IS CRUCIAL DATA, SO IF POSSIBLE FLUSH THE BUFFER TO DISK.    
C     AND FOR GOOD MEASURE, PUT A COPY IN THE PRINTOUT.                 
C                                                                       
      FGONLY = NUM.EQ.0  .AND.  NFRG.GT.0                               
C                                                                       
      IF(.NOT.MASWRK  .OR.  NPRTGO.LE.0) RETURN                         
C                                                                       
      IF(IVIB.EQ.0) WRITE(IFCM,9000) VIBWRD                             
      WRITE(IFCM,9020) BLANK,IVIB,IATOM,ICOORD,E                        
C                                                                       
      IF (.NOT. FGONLY) WRITE(IFCM,9040) EG                             
      IF(NFRG.GT.0) THEN                                                
         DO IFR=1,NFRG                                                  
            WRITE(IFCM,9040) (DEFT(IFRC,IFR),IFRC=1,3)                  
            WRITE(IFCM,9040) (TORQ(IFRC,IFR),IFRC=1,3)                  
         ENDDO                                                          
      END IF                                                            
      WRITE(IFCM,9040) DIP                                              
      CALL FLSHBF(IFCM)                                                 
C                                                                       
      IF(RSTART) RETURN                                                 
      if(nfg.eq.0) then                                                 
      WRITE(IW,9100) VIBWRD                                             
      WRITE(IW,9120) BLANK,IVIB,IATOM,ICOORD,E                          
      IF (.NOT. FGONLY) WRITE(IW,9140) EG                               
C                                                                       
      IF(NFRG.GT.0) THEN                                                
         DO IFR=1,NFRG                                                  
            WRITE(IW,9040) (DEFT(IFRC,IFR),IFRC=1,3)                    
            WRITE(IW,9040) (TORQ(IFRC,IFR),IFRC=1,3)                    
         ENDDO                                                          
      END IF                                                            
      WRITE(IW,9140) DIP                                                
      endif                                                             
      RETURN                                                            
C                                                                       
 9000 FORMAT(A8)                                                        
 9020 FORMAT(A8,' IVIB=',I4,' IATOM=',I4,' ICOORD=',I4,' E=',F20.10)    
 9040 FORMAT(1P,5E16.9)                                                 
 9100 FORMAT(1X,A8)                                                     
 9120 FORMAT(1X,A8,' IVIB=',I4,' IATOM=',I4,' ICOORD=',I4,' E=',F20.10) 
 9140 FORMAT(1X,1P,5E16.9)                                              
      END                                                               
C*MODULE IOLIB   *DECK PKWRIT                                           
      SUBROUTINE PKWRIT(IS,XP,XK,IX,NX,MAX)                             
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      LOGICAL GOPARR,DSKWRK,MASWRK                                      
      DIMENSION XP(*),XK(*),IX(*)                                       
      COMMON /MACHIN/ NWDVAR,MAXFM,MAXSM,LIMFM,LIMSM                    
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
      COMMON /PCKLAB/ LABSIZ                                            
C                                                                       
      IF (LABSIZ.EQ.1) THEN                                             
         IF(NWDVAR.EQ.2) MAX2 = MAX                                     
         IF(NWDVAR.EQ.1) MAX2 = (MAX+1)/2                               
      ELSE IF (LABSIZ.EQ.2) THEN                                        
         IF(NWDVAR.EQ.2) MAX2 = MAX*2                                   
         IF(NWDVAR.EQ.1) MAX2 = MAX                                     
      ELSE                                                              
         IF(MASWRK) WRITE(6,*) 'PKWRIT: CONFUSION WITH -LABSIZ-'        
         CALL ABRT                                                      
      END IF                                                            
      CALL PKWRT(IS,XP,XK,IX,NX, MAX, MAX2 )                            
      RETURN                                                            
      END                                                               
C*MODULE IOLIB   *DECK PKWRT                                            
      SUBROUTINE PKWRT(IS,XP,XK,IX,NX,MAX,MAX2)                         
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      LOGICAL GOPARR,DSKWRK,MASWRK                                      
      DIMENSION XP(MAX),XK(MAX),IX(MAX2)                                
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
C                                                                       
      IF (DSKWRK.OR.MASWRK) WRITE(IS,ERR=300) NX,IX,XP,XK               
      RETURN                                                            
C                                                                       
  300 CONTINUE                                                          
      WRITE(6,9020) ME,IS                                               
 9020 FORMAT(1X,'PKWRT: NODE',I4,                                       
     *          ' ENCOUNTERED I/O ERROR WRITING UNIT',I4)               
      CALL ABRT                                                         
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE IOLIB   *DECK PWRIT                                            
      SUBROUTINE PWRIT(IS,XX,IX,NX,MAX)                                 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      LOGICAL GOPARR,DSKWRK,MASWRK                                      
      DIMENSION XX(MAX),IX(*)                                           
      COMMON /MACHIN/ NWDVAR,MAXFM,MAXSM,LIMFM,LIMSM                    
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
      COMMON /PCKLAB/ LABSIZ                                            
C                                                                       
      IF (LABSIZ.EQ.1) THEN                                             
         IF(NWDVAR.EQ.2) MAX2 = MAX                                     
         IF(NWDVAR.EQ.1) MAX2 = (MAX+1)/2                               
      ELSE IF (LABSIZ.EQ.2) THEN                                        
         IF(NWDVAR.EQ.2) MAX2 = MAX*2                                   
         IF(NWDVAR.EQ.1) MAX2 = MAX                                     
      ELSE                                                              
         IF(MASWRK) WRITE(6,*) 'PWRIT: CONFUSION WITH -LABSIZ-'         
         CALL ABRT                                                      
      END IF                                                            
C                                                                       
      CALL PWRT(IS,XX,IX,NX, MAX, MAX2 )                                
      RETURN                                                            
      END                                                               
C*MODULE IOLIB   *DECK PWRT                                             
      SUBROUTINE PWRT(IS,XX,IX,NX,MAX,MAX2)                             
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      LOGICAL GOPARR,DSKWRK,MASWRK                                      
      DIMENSION XX(MAX),IX(MAX2)                                        
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
C                                                                       
      IF (DSKWRK.OR.MASWRK) WRITE(IS,ERR=300) NX,IX,XX                  
      RETURN                                                            
C                                                                       
  300 CONTINUE                                                          
      WRITE(6,9020) ME,IS                                               
 9020 FORMAT(1X,'PWRT: NODE',I4,                                        
     *          ' ENCOUNTERED I/O ERROR WRITING UNIT',I4)               
      CALL ABRT                                                         
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE IOLIB   *DECK PXWRIT                                           
      SUBROUTINE PXWRIT(NFT,X,XX,MX,NX)                                 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      LOGICAL GOPARR,DSKWRK,MASWRK                                      
      DIMENSION X(NX),XX(NX)                                            
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
C                                                                       
      IF (DSKWRK.OR.MASWRK) WRITE (NFT,ERR=300) X,XX,MX                 
      RETURN                                                            
C                                                                       
  300 CONTINUE                                                          
      WRITE(6,9020) ME,NFT                                              
 9020 FORMAT(1X,'PXWRIT: NODE',I4,                                      
     *          ' ENCOUNTERED I/O ERROR WRITING UNIT',I4)               
      CALL ABRT                                                         
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE IOLIB   *DECK RACLOS                                           
!> @date November 2015, Jerry Boatz                                     
!> - Add the code and debug output option for the 2GB record size workar
!>                                                                      
      SUBROUTINE RACLOS(IRAF,FSTAT)                                     
C                                                                       
      CHARACTER*(*) FSTAT                                               
C                                                                       
      LOGICAL GOPARR,DSKWRK,MASWRK,out,dabigio,f70ok                    
C                                                                       
      PARAMETER (MXUNIT=299)                                            
C                                                                       
      COMMON /CCIO  / f70ok,iosize,lrecfl(mxunit)                       
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)                
      COMMON /OPNNFT/ NFTOPN(MXUNIT),NODEXT(MXUNIT),IOSMP(2)            
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
C                                                                       
C     ----- CLOSE FILE IRAF -----                                       
C                                                                       
      IF(IRAF.GT.MXUNIT) THEN                                           
         WRITE(IW,900) IRAF,MXUNIT                                      
         CALL ABRT                                                      
      END IF                                                            
      if (fstat .eq. 'DELETE') lrecfl(iraf) = 0                         
      IF(NFTOPN(IRAF).EQ.0) RETURN                                      
  900 FORMAT(1X,'RACLOS: ATTEMPT TO CLOSE FILE',I5,                     
     *          ' GREATER THAN MAXIMUM',I5)                             
C                                                                       
      IF(MASWRK  .OR.  DSKWRK) CLOSE(UNIT=IRAF, STATUS=FSTAT)           
      NFTOPN(IRAF) = 0                                                  
      RETURN                                                            
      END                                                               
C*MODULE IOLIB   *DECK RAOPEN                                           
      SUBROUTINE RAOPEN(IRAF,IORA,LPHYS,NUMREC,LENREC,NPRINT)           
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      LOGICAL GOPARR,DSKWRK,MASWRK,ISGDDI,wasgddi,PAROUT,INITGDDI,MLGDDI
C                                                                       
      PARAMETER (MXUNIT=299)                                            
C                                                                       
      CHARACTER*6 FILENM                                                
      CHARACTER*1 NULL                                                  
      CHARACTER*256 PATHNM,ENVBUF                                       
      COMMON /ENVIR / ENVBUF(-5:MXUNIT)                                 
C                                                                       
      DIMENSION IORA(NUMREC)                                            
C                                                                       
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)                
      COMMON /OPNNFT/ NFTOPN(MXUNIT),NODEXT(MXUNIT),IOSMP(2)            
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
      COMMON /RAIOLN/ JRECLN(10),JRECST(10)                             
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL                
      COMMON /GDDI/   ISCOPE,NGROUPS,MYGROUP,MEGLOB,NPGLOB,NNGLOB,JBTYP,
     *                ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,NSUBGR,     
     *                MeUniv,NPUniv,numdlb,myworld,nworlds              
C                                                                       
      DATA CHECK/8HCHECK    /                                           
C                                                                       
C     THIS ROUTINE, AND OTHERS IN THIS GROUP CONTROL THE RANDOM         
C     ACCESS SORTFILE IRAF, WHICH IS USED FOR SORTING MANY THINGS,      
C     SUCH AS ATOMIC INTEGRALS, TRANSFORMED INTEGRALS, MCSCF            
C     HESSIAN MATRIX, 2-PARTICLE DENSITY, ETC.                          
C     EACH TIME A SORT IS DONE, THE DISK REQUIREMENTS CHANGE.           
C     HENCE RAOPEN AND RACLOS BENEFIT MUCH FROM THE FORTRAN-77          
C     OPEN AND CLOSE STATEMENTS. IN THE ABSENCE OF FORTRAN-77,          
C     ONE MAKES FILE 20 LARGE ENOUGH FOR THE BIGGEST SORT.              
C     IF ONE CANNOT TELL HOW BIG THIS IS, ONE OPENS A FIXED             
C     AND VERY LARGE SORTFILE.                                          
C                                                                       
C     -NUMREC- LOGICAL RECORDS OF LENGTH -LENREC- ARE SPREAD OVER AS    
C     MANY PHYSICAL RECORDS OF LENGTH -IRECLN- AS NEEDED, SO THAT       
C     -MXREC- PHYSICAL RECORDS ARE REQUIRED.                            
C     IF -LPHYS- IS NONZERO, AND THE RECORD SIZE -LENREC- IS SMALLER    
C     THAN THE DEFAULT PHYSICAL LENGTH, THE PHYSICAL LENGTH IS          
C     TRIMMED TO JUST THE SIZE NEEDED.                                  
C                                                                       
      IF(IRAF.GT.MXUNIT) THEN                                           
         WRITE(IW,900) IRAF,MXUNIT                                      
         CALL ABRT                                                      
      END IF                                                            
      IF(NFTOPN(IRAF).EQ.1) RETURN                                      
      NFTOPN(IRAF) = 1                                                  
  900 FORMAT(1X,'RAOPEN: ATTEMPT TO OPEN FILE',I5,                      
     *          ' GREATER THAN MAXIMUM',I5)                             
C                                                                       
C     ----- SET THE PHYSICAL RECORD LENGTH -----                        
C                                                                       
      IRECLN=NRASIZ(IRAF)                                               
      IF(LPHYS.NE.0  .AND.  LENREC.LT.IRECLN) IRECLN=LENREC             
C                                                                       
C     ----- COUNT HOW MANY PHYSICAL RECORDS ARE NEEDED -----            
C                                                                       
      MXREC = (((LENREC-1)/IRECLN)+1)*NUMREC                            
      JRECST(IRAF/10)=1                                                 
      JRECLN(IRAF/10)=IRECLN                                            
      IF(EXETYP.EQ.CHECK) RETURN                                        
C                                                                       
C     ----- OPEN THE RANDOM ACCESS FILE -----                           
C                                                                       
C        BY CONVENTION, FILES ENDING IN 0'S ARE DIRECT ACCESS,          
C        AT PRESENT ONLY FILE 20 AND 30 ARE USED.                       
C                                                                       
      FILENM='      '                                                   
      IF(IRAF.EQ.20) FILENM = 'DASORT'                                  
      IF(IRAF.EQ.30) FILENM = 'DAFL30'                                  
      IF(IRAF.EQ.40) FILENM = 'FMODAT'                                  
      IF(FILENM.EQ.'      ') THEN                                       
         IF(MASWRK) WRITE(IW,*) 'ERROR IN RAOPEN WITH UNIT NUMBER',IRAF 
         CALL ABRT                                                      
      END IF                                                            
C                                                                       
      IF(NPRINT.NE.-5  .AND.  MASWRK)                                   
     *     WRITE(IW,9000) FILENM,NUMREC,LENREC,MXREC,IRECLN             
C                                                                       
      IF (MASWRK  .OR.  DSKWRK) THEN                                    
C                                                                       
*IBM  OPEN (UNIT=IRAF, FILE=FILENM, STATUS='UNKNOWN',                   
*IBM *      ACCESS='DIRECT', FORM='UNFORMATTED', RECL=8*IRECLN)         
C                                                                       
      IF(MASWRK) CALL GMS_GETENV(FILENM,PATHNM)                         
      IF(ISGDDI.or.wasgddi) THEN                                        
         PATHNM=ENVBUF(IRAF)                                            
         IF(NODEXT(IRAF).EQ.0)                                          
     *       CALL ADDNANODE(PATHNM,MeUniv,IRAF,iout)                    
      ELSE                                                              
         CALL PARENV(FILENM,PATHNM,IOUT)                                
         IF (IOUT.EQ.1) RETURN                                          
      ENDIF                                                             
      NULL = CHAR(0)                                                    
      DO 3 KOL=1,256                                                    
         IF(PATHNM(KOL:KOL).EQ.' '  .OR.                                
     *      PATHNM(KOL:KOL).EQ.NULL) GO TO 4                            
    3 CONTINUE                                                          
      KOL=257                                                           
    4 CONTINUE                                                          
      IF(KOL.EQ.1) THEN                                                 
         WRITE(IW,1) FILENM                                             
         CALL ABRT                                                      
      END IF                                                            
      KOL=KOL-1                                                         
      OPEN (UNIT=IRAF, FILE=PATHNM(1:KOL), STATUS='UNKNOWN',            
     *      ACCESS='DIRECT', FORM='UNFORMATTED',                        
     *      RECL=8*IRECLN)                                              
    1 FORMAT(1X,'YOU MUST ASSIGN GENERIC NAME ',A,' WITH A SETENV.')    
C                                                                       
*VMS  OPEN (UNIT=IRAF, FILE=FILENM, STATUS='UNKNOWN',                   
*VMS *      ACCESS='DIRECT', FORM='UNFORMATTED', RECL=2*IRECLN)         
C                                                                       
      END IF                                                            
C                                                                       
C      ----- INITIALIZE FOR EMPTY SORTFILE ----                         
C                                                                       
      DO 100 I = 1,NUMREC                                               
         IORA(I) = -1                                                   
  100 CONTINUE                                                          
      RETURN                                                            
C                                                                       
 9000 FORMAT(1X,'OPENING FILE ',A6,' WITH',I8,' LOGICAL RECORDS OF',I8, 
     *          ' WORDS'/1X,'WITH A MAXIMUM OF',I12,                    
     *          ' PHYSICAL RECORDS OF',I8,' WORDS')                     
      END                                                               
C*MODULE IOLIB   *DECK RAOPEN2                                          
      SUBROUTINE RAOPEN2(IRAF,IORA,LPHYS,NUMREC,LENREC,LDAR,NPRINT)     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      LOGICAL GOPARR,DSKWRK,MASWRK,ISGDDI,wasgddi,PAROUT,INITGDDI,MLGDDI
      PARAMETER (MXUNIT=299)                                            
C                                                                       
      CHARACTER*7 FILENM                                                
      CHARACTER*1 NULL                                                  
      CHARACTER*256 PATHNM,ENVBUF                                       
      COMMON /ENVIR / ENVBUF(-5:MXUNIT)                                 
C                                                                       
      DIMENSION IORA(NUMREC)                                            
C                                                                       
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)                
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
      COMMON /RAIOLN/ JRECLN(10),JRECST(10)                             
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL                
      COMMON /GDDI/   ISCOPE,NGROUPS,MYGROUP,MEGLOB,NPGLOB,NNGLOB,JBTYP,
     *                ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,NSUBGR,     
     *                MeUniv,NPUniv,numdlb,myworld,nworlds              
      COMMON /OPNNFT/ NFTOPN(MXUNIT),NODEXT(MXUNIT),IOSMP(2)            
C                                                                       
      DATA CHECK/8HCHECK    /                                           
C                                                                       
C     THIS ROUTINE, AND OTHERS IN THIS GROUP CONTROL THE RANDOM         
C     ACCESS SORTFILE IRAF, WHICH IS USED FOR SORTING MANY THINGS,      
C     SUCH AS ATOMIC INTEGRALS, TRANSFORMED INTEGRALS, MCSCF            
C     HESSIAN MATRIX, 2-PARTICLE DENSITY, ETC.                          
C     EACH TIME A SORT IS DONE, THE DISK REQUIREMENTS CHANGE.           
C     HENCE RAOPEN AND RACLOS BENEFIT MUCH FROM THE FORTRAN-77          
C     OPEN AND CLOSE STATEMENTS. IN THE ABSENCE OF FORTRAN-77,          
C     ONE MAKES FILE 20 LARGE ENOUGH FOR THE BIGGEST SORT.              
C     IF ONE CANNOT TELL HOW BIG THIS IS, ONE OPENS A FIXED             
C     AND VERY LARGE SORTFILE.                                          
C                                                                       
C     DIFFERS FROM -RAOPEN- IN THAT THERE IS A FILE 40, AND             
C     THE SIZE OF THE RECORDS IS A CALLING ARGUMENT -LDAR-.             
C                                                                       
C     -NUMREC- LOGICAL RECORDS OF LENGTH -LENREC- ARE SPREAD OVER AS    
C     MANY PHYSICAL RECORDS OF LENGTH -IRECLN- AS NEEDED, SO THAT       
C     -MXREC- PHYSICAL RECORDS ARE REQUIRED.                            
C     IF -LPHYS- IS NONZERO, AND THE RECORD SIZE -LENREC- IS SMALLER    
C     THAN THE DEFAULT PHYSICAL LENGTH, THE PHYSICAL LENGTH IS          
C     TRIMMED TO JUST THE SIZE NEEDED.                                  
C                                                                       
C     ----- SET THE PHYSICAL RECORD LENGTH -----                        
C                                                                       
      IF(LDAR.EQ.0) THEN                                                
         IRECLN=NRASIZ(IRAF)                                            
      ELSE                                                              
         IRECLN=LDAR                                                    
      END IF                                                            
      IF(LPHYS.NE.0  .AND.  LENREC.LT.IRECLN) IRECLN=LENREC             
C                                                                       
C     ----- COUNT HOW MANY PHYSICAL RECORDS ARE NEEDED -----            
C                                                                       
      MXREC = (((LENREC-1)/IRECLN)+1)*NUMREC                            
      JRECST(IRAF/10)=1                                                 
      JRECLN(IRAF/10)=IRECLN                                            
      IF(EXETYP.EQ.CHECK) RETURN                                        
C                                                                       
C        BY CONVENTION, FILES ENDING IN 0'S ARE DIRECT ACCESS,          
C        AT PRESENT ONLY FILE 20, 30, AND 40 ARE USED.                  
C                                                                       
      FILENM='      '                                                   
      NFILNM=6                                                          
      IF(IRAF.EQ.20) FILENM = 'DASORT'                                  
      IF(IRAF.EQ.30) FILENM = 'DAFL30'                                  
      IF(IRAF.EQ.40) FILENM = 'SOCCDAT'                                 
      IF(IRAF.EQ.40) NFILNM = 7                                         
      IF(FILENM(1:1).EQ.' ') THEN                                       
         IF(MASWRK) WRITE(IW,*) 'ERROR IN RAOPEN2 WITH UNIT NUMBER',IRAF
         CALL ABRT                                                      
      END IF                                                            
C                                                                       
      IF(NPRINT.NE.-5  .AND.  MASWRK)                                   
     *     WRITE(IW,9000) FILENM(1:NFILNM),NUMREC,LENREC,MXREC,IRECLN   
C                                                                       
C     ----- OPEN THE RANDOM ACCESS FILE -----                           
C                                                                       
      IF (MASWRK  .OR.  DSKWRK) THEN                                    
C                                                                       
*IBM  OPEN (UNIT=IRAF, FILE=FILENM(1:NFILNM), STATUS='UNKNOWN',         
*IBM *      ACCESS='DIRECT', FORM='UNFORMATTED', RECL=8*IRECLN)         
C                                                                       
      IF(MASWRK) CALL GMS_GETENV(FILENM(1:NFILNM),PATHNM)               
      IF(ISGDDI.or.wasgddi) THEN                                        
         PATHNM=ENVBUF(IRAF)                                            
         IF(NODEXT(IRAF).EQ.0)                                          
     *       CALL ADDNANODE(PATHNM,MeUniv,IRAF,iout)                    
      ELSE                                                              
         CALL PARENV(FILENM(1:NFILNM),PATHNM,IOUT)                      
         IF (IOUT.EQ.1) RETURN                                          
      ENDIF                                                             
      NULL = CHAR(0)                                                    
      DO 3 KOL=1,256                                                    
         IF(PATHNM(KOL:KOL).EQ.' '  .OR.                                
     *      PATHNM(KOL:KOL).EQ.NULL) GO TO 4                            
    3 CONTINUE                                                          
      KOL=257                                                           
    4 CONTINUE                                                          
      IF(KOL.EQ.1) THEN                                                 
         WRITE(IW,1) FILENM(1:NFILNM)                                   
         CALL ABRT                                                      
      END IF                                                            
      KOL=KOL-1                                                         
      OPEN (UNIT=IRAF, FILE=PATHNM(1:KOL), STATUS='UNKNOWN',            
     *      ACCESS='DIRECT', FORM='UNFORMATTED',                        
     *      RECL=8*IRECLN)                                              
    1 FORMAT(1X,'YOU MUST ASSIGN GENERIC NAME ',A,' WITH A SETENV.')    
C                                                                       
*VMS  OPEN (UNIT=IRAF, FILE=FILENM(1:NFILNM), STATUS='UNKNOWN',         
*VMS *      ACCESS='DIRECT', FORM='UNFORMATTED', RECL=2*IRECLN)         
C                                                                       
      END IF                                                            
C                                                                       
C      ----- INITIALIZE FOR EMPTY SORTFILE ----                         
C                                                                       
      DO 100 I = 1,NUMREC                                               
         IORA(I) = -1                                                   
  100 CONTINUE                                                          
      RETURN                                                            
C                                                                       
 9000 FORMAT(1X,'OPENING FILE ',A,' WITH',I8,' LOGICAL RECORDS OF',I8,  
     *          ' WORDS'/1X,'WITH A MAXIMUM OF',I12,                    
     *          ' PHYSICAL RECORDS OF',I8,' WORDS')                     
      END                                                               
C*MODULE IOLIB   *DECK RARD                                             
      SUBROUTINE RARD(V,LEN,IRAF,NS)                                    
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      LOGICAL GOPARR,DSKWRK,MASWRK                                      
      DIMENSION V(LEN)                                                  
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
C                                                                       
C     ----- READ A PHYSICAL RECORD FROM IRAF -----                      
C                                                                       
      IF(MASWRK  .OR.  DSKWRK) READ(UNIT=IRAF,REC=NS) V                 
C                                                                       
      IF(GOPARR  .AND.  .NOT.DSKWRK) THEN                               
         CALL DDI_BCAST(225,'F',V,LEN,MASTER)                           
      END IF                                                            
      RETURN                                                            
      END                                                               
C*MODULE IOLIB   *DECK RAREAD                                           
      SUBROUTINE RAREAD(IRAF,IORA,V,LEN,NREC,NAV)                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION V(*),IORA(*)                                            
      COMMON /RAIOLN/ JRECLN(10),JRECST(10)                             
C                                                                       
C      ----- READ A LOGICAL RECORD FROM IRAF -----                      
C                                                                       
      IF(IRAF.LT.0) WRITE(6,*) 'BOGUS RAREAD, NAV=',NAV                 
C                                                                       
      IRECLN=JRECLN(IRAF/10)                                            
      N = IORA(NREC)                                                    
      IS = -IRECLN + 1                                                  
      NS = N                                                            
      LENT = LEN                                                        
C                                                                       
  100 CONTINUE                                                          
         IS = IS + IRECLN                                               
         IF = IS + LENT - 1                                             
         IF ((IF-IS+1) .GT. IRECLN) IF = IS + IRECLN - 1                
         NSP = NS                                                       
         LENW = IF - IS + 1                                             
         CALL RARD(V(IS),LENW,IRAF,NSP)                                 
         LENT = LENT - IRECLN                                           
         NS = NS + 1                                                    
         N = NS                                                         
      IF (LENT .GE. 1) GO TO 100                                        
      RETURN                                                            
      END                                                               
C*MODULE IOLIB   *DECK RASIZE                                           
      SUBROUTINE RASIZE(LDAR)                                           
C                                                                       
C     THIS ROUTINE RETURNS THE LENGTH OF PHYSICAL RECORDS (IN UNITS OF  
C     W.P. FLOATING POINT NOS) USED BY THE DIRECT ACCESS SORT FILE IRAF.
C     SELECT A ROUND NUMBER FOR LEN, A BIT SMALLER THAN THE PHYSICAL    
C     RECORD LENGTH USED IN ROUTINE RAOPEN.  SEE ROUTINE RAOPEN.        
C                                                                       
*VMS  LDAR=2045                                                         
*IBM  LDAR=2930                                                         
      LDAR=2045                                                         
      RETURN                                                            
      END                                                               
C*MODULE IOLIB   *DECK NRASIZ                                           
!> @date June 2015, Jerry Boatz                                         
!> -  Add section to implement "chunking" of CC files                   
!>                                                                      
      INTEGER FUNCTION NRASIZ(IRAF)                                     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
c-remove??                                                              
c-2gb                                                                   
c     logical dabigio                                                   
c-2gb                                                                   
c-remove??                                                              
C                                                                       
C     USED BY DIRECT ACCESS FILE OPENING ROUTINES                       
C                                                                       
      IF(IRAF.EQ.10) THEN                                               
*IBM      NRASIZ = 2048                                                 
          NRASIZ = 4090                                                 
*VAX      NRASIZ = 2048                                                 
      ELSE                                                              
*IBM      NRASIZ = 2934                                                 
          NRASIZ = 2048                                                 
*VAX      NRASIZ = 2046                                                 
C         PLEASE USE AN EVEN NUMBER                                     
      ENDIF                                                             
c-remove??                                                              
c     if (dabigio(iraf)) then                                           
c*IBM      NRASIZ = 293400                                              
c*UNX      NRASIZ = 204800                                              
c*VAX      NRASIZ = 204600                                              
c-2gb-test                                                              
c                                                                       
c     --- try a very small record size, so that all files are           
c     --- chopped into smaller pieces, even for small test jobs         
c                                                                       
c         NRASIZ = 200                                                  
c                                                                       
c--- efficiency testing...                                              
c         NRASIZ = 20480                                                
c         NRASIZ = 2048000                                              
c         NRASIZ = 20480000                                             
c--- efficiency testing...                                              
c         write(6,*)                                                    
c         write(6,*) '*** NOTE ***'                                     
c         write(6,*) 'In NRASIZ:  default record size set to ',nrasiz   
c         write(6,*)                                                    
c-2gb-test                                                              
c     end if                                                            
c-remove??                                                              
      RETURN                                                            
      END                                                               
C*MODULE IOLIB   *DECK RAWRIT                                           
      SUBROUTINE RAWRIT(IRAF,IORA,V,LEN,NREC,NAVM)                      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION V(*),IORA(*)                                            
      LOGICAL GOPARR,DSKWRK,MASWRK                                      
C                                                                       
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
      COMMON /RAIOLN/ JRECLN(10),JRECST(10)                             
C                                                                       
C     ------ WRITE A LOGICAL RECORD ON IRAF -----                       
C                                                                       
      IF(IRAF.LT.0) WRITE(6,*) 'BOGUS RAWRIT, NAVM=',NAVM               
C                                                                       
      IRECST=JRECST(IRAF/10)                                            
      IRECLN=JRECLN(IRAF/10)                                            
      N = IORA(NREC)                                                    
      IF(N.LE.0) THEN                                                   
         IORA(NREC) = IRECST                                            
         IRECST = IRECST + (LEN-1)/IRECLN + 1                           
         N = IORA(NREC)                                                 
      END IF                                                            
C                                                                       
      JRECST(IRAF/10)=IRECST                                            
      IF(GOPARR  .AND.  .NOT.(MASWRK.OR.DSKWRK)) RETURN                 
C                                                                       
      IST = -IRECLN + 1                                                 
      NS = N                                                            
      LENT = LEN                                                        
  120 CONTINUE                                                          
         IST = IST + IRECLN                                             
         IF = IST + LENT - 1                                            
         IF ((IF-IST+1) .GT. IRECLN) IF = IST+IRECLN-1                  
         NSP = NS                                                       
         LENW = IF - IST + 1                                            
         CALL RAWRT(V(IST),LENW,IRAF,NSP)                               
         LENT = LENT - IRECLN                                           
         NS = NS + 1                                                    
         N = NS                                                         
      IF (LENT .GE. 1) GO TO 120                                        
      RETURN                                                            
C                                                                       
      END                                                               
C*MODULE IOLIB   *DECK RAWRITE                                          
      SUBROUTINE RAWRITE(IRAF,IORA,V,LEN1,LEN,NREC,NPHYSOFF)            
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION V(*),IORA(*)                                            
      LOGICAL GOPARR,DSKWRK,MASWRK                                      
C                                                                       
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
      COMMON /RAIOLN/ JRECLN(10),JRECST(10)                             
C                                                                       
C     ------ WRITE A LOGICAL RECORD ON IRAF -----                       
C     THIS IS THE SAME AS RAWRIT EXCEPT IT ENFORCES THAT LEN1 WORDS     
C     BE ALLOCATED ON IRAF FOR THIS RECORD EVEN IF ONE WRITES ONLY      
C     LEN WORDS THIS TIME                                               
C                                                                       
C     NPHYSOFF SETS AN OFFSET IN PHYSICAL RECORDS                       
C     (IN ORDER TO WRITE A PARTIAL LOGICAL ERCORD)                      
C                                                                       
      IRECST=JRECST(IRAF/10)                                            
      IRECLN=JRECLN(IRAF/10)                                            
      N = IORA(NREC)                                                    
      IF(N.LE.0) THEN                                                   
         IORA(NREC) = IRECST                                            
         IRECST = IRECST + (LEN1-1)/IRECLN + 1                          
         N = IORA(NREC)                                                 
      END IF                                                            
C                                                                       
      JRECST(IRAF/10)=IRECST                                            
      IF(GOPARR  .AND.  .NOT.(MASWRK.OR.DSKWRK)) RETURN                 
C                                                                       
      IST = -IRECLN + 1                                                 
      NS = N                                                            
      LENT = LEN                                                        
  120 CONTINUE                                                          
         IST = IST + IRECLN                                             
         IF = IST + LENT - 1                                            
         IF ((IF-IST+1) .GT. IRECLN) IF = IST+IRECLN-1                  
         NSP = NS                                                       
         LENW = IF - IST + 1                                            
         CALL RAWRT(V(IST),LENW,IRAF,NSP+NPHYSOFF)                      
         LENT = LENT - IRECLN                                           
         NS = NS + 1                                                    
         N = NS                                                         
      IF (LENT .GE. 1) GO TO 120                                        
      RETURN                                                            
      END                                                               
C*MODULE IOLIB   *DECK RAWRT                                            
      SUBROUTINE RAWRT(V,LEN,IRAF,NS)                                   
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION V(LEN)                                                  
C                                                                       
C     ----- WRITE A PHYSICAL RECORD ON IRAF -----                       
C                                                                       
      WRITE (UNIT=IRAF, REC=NS) V                                       
      RETURN                                                            
      END                                                               
C*MODULE IOLIB   *DECK SEQADV                                           
      SUBROUTINE SEQADV(LUNIT)                                          
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      LOGICAL GOPARR,DSKWRK,MASWRK                                      
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
C                                                                       
C     ----- ADVANCE ONE RECORD ON FILE LUNIT -----                      
C                                                                       
      IF (DSKWRK.OR.MASWRK) READ (LUNIT)                                
      RETURN                                                            
      END                                                               
C*MODULE IOLIB   *DECK SEQCLO                                           
      SUBROUTINE SEQCLO(IFILE,FSTAT)                                    
C                                                                       
      CHARACTER*(*) FSTAT                                               
C                                                                       
      LOGICAL GOPARR,DSKWRK,MASWRK                                      
      PARAMETER (MXUNIT=299)                                            
C                                                                       
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)                
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
      COMMON /OPNNFT/ NFTOPN(MXUNIT),NODEXT(MXUNIT),IOSMP(2)            
C                                                                       
C     ----- CLOSE THE SEQUENTIAL FILE IFILE -----                       
C     ----- FSTAT MAY BE 'KEEP' OR 'DELETE' -----                       
C                                                                       
      IF(IFILE.GT.MXUNIT) THEN                                          
         WRITE(IW,900) IFILE,MXUNIT                                     
         CALL ABRT                                                      
      END IF                                                            
      IF(NFTOPN(IFILE).EQ.0) RETURN                                     
C                                                                       
      IF(DSKWRK.OR.MASWRK) CLOSE (UNIT=IFILE, STATUS=FSTAT)             
      NFTOPN(IFILE) = 0                                                 
      RETURN                                                            
  900 FORMAT(1X,'SEQCLO: ATTEMPT TO CLOSE FILE',I5,                     
     *          ' GREATER THAN MAXIMUM',I5)                             
      END                                                               
C*MODULE IOLIB   *DECK SEQOPN                                           
!> @brief Open a sequential file                                        
!>                                                                      
!> @author Michael W. Schmidt                                           
!> @date September 2010                                                 
!> - Add the CIMFILE to flush-on-write                                  
      SUBROUTINE SEQOPN(IUNIT,FNAME,FSTAT,RDONLY,FMT)                   
C                                                                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C                                                                       
      CHARACTER*(*) FNAME,FSTAT,FMT                                     
      LOGICAL RDONLY,GOPARR,DSKWRK,MASWRK,ISGDDI,wasgddi,PAROUT,INITGDDI
     *       ,MLGDDI                                                    
      PARAMETER (MXUNIT=299)                                            
C                                                                       
      CHARACTER*256 FILENM,ENVBUF                                       
      CHARACTER*1 NULL                                                  
      COMMON /ENVIR / ENVBUF(-5:MXUNIT)                                 
C                                                                       
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)                
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
      COMMON /OPNNFT/ NFTOPN(MXUNIT),NODEXT(MXUNIT),IOSMP(2)            
      COMMON /GDDI/   ISCOPE,NGROUPS,MYGROUP,MEGLOB,NPGLOB,NNGLOB,JBTYP,
     *                ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,NSUBGR,     
     *                MeUniv,NPUniv,numdlb,myworld,nworlds              
C                                                                       
C      OPEN SEQUENTIAL FILE -IUNIT- WITH GENERIC NAME -FNAME-           
C      WITH STATUS -FSTAT-, AS EITHER READONLY/READWRITE, AND           
C      FOR EITHER FORMATTED OR UNFORMATTED I/O. EXAMPLE:                
C           CALL SEQOPN(IS,'AOINTS','UNKNOWN',.FALSE.,'UNFORMATTED')    
C                                                                       
C                                                                       
C      SOME COMPILERS (NOTABLY F2C OR G77) MAY LIMIT THE UNIT NUMBER    
C      TO A SMALL VALUE (99 IN THE CASE OF THOSE TWO).  A 'SED HACK'    
C      CAN BE USED TO ADJUST -MXUNIT- BEFORE COMPILING, IF SO. ANY      
C      ATTEMPT TO OPEN TOO LARGE A UNIT NUMBER SHOULD BE TERMINATED.    
C      IF THE LIMIT IS 99, ASSUME THIS IS G77, AND PRINT A MESSAGE      
C      ABOUT THE ONLY RUNTYP THAT PRESENTLY GOES OVER 99.               
C                                                                       
      MAXFILES = MXUNIT  ! TEMP VALUE SUPPRESSES A DIAGNOSTIC ON AXP64  
      IF(IUNIT.GT.MAXFILES) THEN                                        
         WRITE(IW,900) IUNIT,MAXFILES                                   
         IF(MAXFILES.EQ.99) WRITE(IW,901)                               
         CALL ABRT                                                      
      END IF                                                            
  900 FORMAT(1X,'SEQOPN: ATTEMPT TO OPEN FILE',I5,                      
     *          ' GREATER THAN MAXIMUM',I5)                             
  901 FORMAT(1X,'G77 HAS AN INTERNAL LIMIT OF 99 FILES,',               
     *          ' AND SO CANNOT EXECUTE RUNTYP=TDHFX')                  
C                                                                       
C         -NFTOPN- IS FOR THE AIX XL FORTRAN VERSION 3, WHICH DOESN'T   
C         SILENTLY IGNORE REQUESTS TO OPEN A FILE WHICH IS ALREADY      
C         OPEN.  IF ANY FILE HAS ALREADY BEEN OPENED, WE JUST EXIT,     
C         WHICH IS THE USUAL FORTRAN BEHAVIOR.  0=NOT OPEN, 1=OPEN.     
C                                                                       
C         -NODEXT- CONTROLS WHETHER TO ADD PROCESS (RANK) EXTENSION     
C         OF .001, .002, ETC TO EACH FILE NAME: 0=ADD, 1=DO NOT ADD.    
C         NOTE THAT -NODEXT- IS POSSIBLE INPUT (LATER) IN $SYSTEM.      
C                                                                       
C         THE INPUT FILE -5- IS THE FIRST FILE OPENED, AND THUS A       
C         GOOD TRIGGER FOR DOING THESE INITIALIZATIONS.  (2=ERICFMT)    
C                                                                       
      IF(IUNIT.EQ.5.and..not.ISGDDI) THEN                               
         DO I=1,MXUNIT                                                  
            NFTOPN(I)=0                                                 
            NODEXT(I)=0                                                 
         ENDDO                                                          
         NODEXT(2)=1                                                    
         NODEXT(5)=1                                                    
         IF(.NOT.RDONLY) WRITE(IW,*) 'CAUTION, INPUT FILE NOT READ ONLY'
      END IF                                                            
C                                                                       
      IF(NFTOPN(IUNIT).EQ.1) RETURN                                     
      NFTOPN(IUNIT) = 1                                                 
C                                                                       
*IBM  IF (MASWRK) THEN                                                  
*IBM  IF(FNAME.NE.'OUTPUT') THEN                                        
*IBM     OPEN(UNIT=IUNIT, FILE=FNAME, STATUS=FSTAT,                     
*IBM *        ACCESS='SEQUENTIAL', FORM=FMT)                            
*IBM  END IF                                                            
*IBM  END IF                                                            
C                                                                       
C     FOR NON-GDDI RUNS,                                                
C         FIRST, TRANSLATE THE ENVIRONMENT VARIABLE TO A FULLY          
C                QUALIFIED FILE NAME, FNAME --> FILENM.                 
C          NEXT, APPEND AN EXTENSION TO THE FILE NAME, GIVING ITS       
C                RANK, SUCH AS .001, ON ALL PROCESSES WHICH ARE NOT     
C                THE MASTER.  THUS SMP SYSTEMS HAVE UNIQUE FILE NAMES.  
C         THIRD, STICK IN A NULL BYTE SO THAT UNIX COMPILERS, WHICH     
C                DON'T UNDERSTAND THE FORTRAN CHARACTER TYPE, HAVE      
C                A C-STYLE CHAR VARIABLE TO WORK WITH (UGH).            
C   AND FINALLY, OPEN THE FILE.                                         
C                                                                       
C     FOR GDDI RUNS,                                                    
C       THE ENVIRONMENT VARIABLES ARE ALREADY TRANSLATED, SEE -STORENV- 
C       THE RANK EXTENSIONS MAY NOT BE WANTED (SEE -NODEXT- IN $SYSTEM) 
C                                                                       
C     MOST PEOPLE DO NOT DEFINE 'OUTPUT' IN -RUNGMS-, WHICH LEADS TO    
C     NO OPEN STATEMENT BEING EXECUTED: PRINT OUT TO "STANDARD OUTPUT". 
C     IN CONTRAST, UNIT 5 IS ALWAYS A DISK FILE, AND IS ALWAYS OPENED.  
C                                                                       
      FILENM=' '                                                        
      IF (MASWRK) CALL GMS_GETENV(FNAME,FILENM)                         
      IF(ISGDDI.or.wasgddi) THEN                                        
        FILENM=ENVBUF(IUNIT)                                            
        IF(NODEXT(IUNIT).EQ.0) then                                     
           CALL ADDNANODE(FILENM,MeUniv,IUNIT,iout)                     
           IF (IOUT.EQ.1) RETURN                                        
        ENDIF                                                           
      ELSE                                                              
        CALL PARENV(FNAME,FILENM,IOUT)                                  
        IF (IOUT.EQ.1) RETURN                                           
      END IF                                                            
      IF(FNAME.EQ.'OUTPUT') THEN                                        
         IF(FILENM(1:1).NE.' '  .AND.  FILENM(1:6).NE.'OUTPUT') THEN    
            OPEN(UNIT=IUNIT, FILE=FILENM, STATUS=FSTAT,                 
     *           ACCESS='SEQUENTIAL', FORM=FMT)                         
         END IF                                                         
      ELSE                                                              
         NULL = CHAR(0)                                                 
         DO 1 KOL=1,256                                                 
            IF(FILENM(KOL:KOL).EQ.' '  .OR.                             
     *         FILENM(KOL:KOL).EQ.NULL) GO TO 2                         
    1    CONTINUE                                                       
         KOL=257                                                        
    2    CONTINUE                                                       
         IF(KOL.EQ.1) THEN                                              
            WRITE(IW,3) FNAME                                           
            CALL ABRT                                                   
         END IF                                                         
         KOL=KOL-1                                                      
c        if(IUNIT.eq.5) write(6,*) 'wwwopen',me,IUNIT,FILENM(1:KOL)     
         OPEN(UNIT=IUNIT, FILE=FILENM(1:KOL), STATUS=FSTAT,             
     *        ACCESS='SEQUENTIAL', FORM=FMT, ERR=4)                     
      END IF                                                            
      RETURN                                                            
C                                                                       
C         IF YOU HAVE A F90 COMPILER, AND ARE TRYING TO KEEP SOME OF    
C         THE FILES LIKE -ERICFMT- AND -MCPDATA- READ-ONLY, YOU COULD   
C         CHANGE THE OPEN STATEMENT ABOVE TO:    (USE *'S IN COLUMN 1)  
CUNX     IF (RDONLY) THEN                                               
CUNX        OPEN(UNIT=IUNIT, FILE=FILENM(1:KOL), STATUS=FSTAT,          
CUNX *           ACCESS='SEQUENTIAL', FORM=FMT, ERR=4, ACTION='READ')   
CUNX     ELSE                                                           
CUNX        OPEN(UNIT=IUNIT, FILE=FILENM(1:KOL), STATUS=FSTAT,          
CUNX *           ACCESS='SEQUENTIAL', FORM=FMT, ERR=4)                  
CUNX     END IF                                                         
C                                                                       
C         ERROR HANDLING (E.G. FOR ERICFMT, ...)                        
C                                                                       
    3 FORMAT(1X,'YOU MUST ASSIGN GENERIC NAME ',A,' WITH A SETENV.')    
C                                                                       
    4 CONTINUE                                                          
      IF(FSTAT.EQ.'OLD') THEN                                           
         IF(MASWRK) WRITE(IW,5) 'PRE-EXISTING',FNAME,FILENM(1:KOL)      
      ELSE                                                              
         IF(MASWRK) WRITE(IW,5) 'NEW',FNAME,FILENM(1:KOL)               
      ENDIF                                                             
      CALL ABRT                                                         
    5 FORMAT(//1X,'ERROR OPENING ',A,' FILE ',A,','/                    
     *         1X,'ASSIGNED TO EXPLICIT FILE NAME ',A,','/              
     *         1X,'PLEASE CHECK THE -SETENV- FILE ASSIGNMENTS',         
     *            ' IN YOUR -RUNGMS- SCRIPT.')                          
C                                                                       
*VMS  IF(RDONLY) THEN                                                   
*VMS     OPEN(UNIT=IUNIT, FILE=FNAME, STATUS=FSTAT, READONLY,           
*VMS *        ACCESS='SEQUENTIAL', FORM=FMT, SHARED)                    
*VMS  ELSE                                                              
*VMS     IF(FNAME.EQ.'OUTPUT') THEN                                     
*VMS        OPEN(UNIT=IUNIT, FILE=FNAME, STATUS=FSTAT,                  
*VMS *           ACCESS='SEQUENTIAL', FORM=FMT, SHARED)                 
*VMS     ELSE IF(FNAME.EQ.'PUNCH'   .OR. FNAME.EQ.'TRAJECT'             
*VMS             FNAME.EQ.'RESTART' .OR. FNAME.EQ.'MAKEFP') THEN        
*VMS        OPEN(UNIT=IUNIT, FILE=FNAME, STATUS=FSTAT,                  
*VMS *           ACCESS='SEQUENTIAL', FORM=FMT,                         
*VMS *           CARRIAGECONTROL='LIST', SHARED)                        
*VMS     ELSE                                                           
*VMS        OPEN(UNIT=IUNIT, FILE=FNAME, STATUS=FSTAT,                  
*VMS *           ACCESS='SEQUENTIAL', FORM=FMT)                         
*VMS     END IF                                                         
*VMS  END IF                                                            
C                                                                       
C       IT MAY BE USEFUL TO SET FLUSH-ON-WRITE FOR ASCII OUTPUT FILES,  
C       NOTE THAT WE WE SKIP HERE UNIT 6 WHICH SHOULD BE OK ALREADY.    
C       THIS IS PRIMARILY AN AIX THING!                                 
C                                                                       
      IF(FNAME.EQ.'PUNCH')   CALL FLUSHONWRITE(IUNIT)                   
      IF(FNAME.EQ.'RESTART') CALL FLUSHONWRITE(IUNIT)                   
      IF(FNAME.EQ.'TRAJECT') CALL FLUSHONWRITE(IUNIT)                   
      IF(FNAME.EQ.'MAKEFP')  CALL FLUSHONWRITE(IUNIT)                   
      if(fname.eq.'CIMFILE') call flushonwrite(iunit)                   
      RETURN                                                            
      END                                                               
C*MODULE IOLIB   *DECK SEQREW                                           
      SUBROUTINE SEQREW(IFILE)                                          
C                                                                       
      LOGICAL GOPARR,DSKWRK,MASWRK                                      
C                                                                       
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
C                                                                       
C   ----- REWIND FILE IFILE ----                                        
C                                                                       
      IF (DSKWRK.OR.MASWRK) REWIND (UNIT=IFILE, ERR=300)                
  300 CONTINUE                                                          
      RETURN                                                            
      END                                                               
C*MODULE IOLIB   *DECK SQRDCH                                           
      SUBROUTINE SQRDCH(LFILE,REGION,LENGTH)                            
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      DIMENSION REGION(LENGTH)                                          
C                                                                       
C         READ A SEQUENTIAL RECORD FROM FILE -LFILE- OF -LENGTH- WORDS. 
C         SEE ALSO THE ROUTINE SQWRCH, WHICH IS SQRDCH'S PARTNER.       
C                                                                       
C         IN CASE OF END OF FILE DETECTION, -LENGTH- IS RETURNED AS 0!  
C                                                                       
C         IN CASE THE -LENGTH- IS VERY LARGE, THE I/O IS DIVIDED INTO   
C         SMALLER SUBRECORDS, OF -MAXSEQ- WORDS.  AS AN EXAMPLE, IF     
C         THE DISK STRIPE SIZE IS 32K, A RECORD OF 250,000 WORDS WILL   
C         CONSIST OF ABOUT 61 STRIPES, AND IS THUS A RATHER BIG         
C         RECORD SIZE, AT LEAST IN THE YEAR 2005.                       
C                                                                       
C         NOTE THAT CHUNKING IS INCOMPATIBLE WITH THE ABILITY TO USE    
C         THE -SEQADV- ROUTINE.                                         
C                                                                       
      MAXSEQ = 250000                                                   
C                                                                       
      IADD = 1                                                          
  100 CONTINUE                                                          
         LEN = MIN(MAXSEQ,LENGTH-IADD+1)                                
         CALL SQREAD(LFILE,REGION(IADD),LEN)                            
C               THIS RETURNS THE END OF FILE HAS BEEN DETECTED          
         IF(LEN.EQ.0) LENGTH=0                                          
         IADD = IADD+MAXSEQ                                             
         IF(IADD.GT.LENGTH) RETURN                                      
      GO TO 100                                                         
      END                                                               
C*MODULE IOLIB   *DECK SQREAD                                           
      SUBROUTINE SQREAD(LFILE,REGION,LENGTH)                            
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      DIMENSION REGION(LENGTH)                                          
      LOGICAL GOPARR,DSKWRK,MASWRK                                      
      COMMON /IOFILE/ IR,IW,IP,IIS,IPK,IDAFX,NAV,IODAX(950)             
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
C                                                                       
C     READ -LENGTH- FLOATING POINT WORDS INTO -REGION- FROM -LFILE-     
C     A POSSIBLE END-OF-FILE RESULTS IN RETURNING -LENGTH- AS ZERO!     
C     NOTE THAT DUE TO ITS POSSIBLE RESET TO 0, -LENGTH- SHOULD         
C     BE A INTEGER VARIABLE, RATHER THAN AN INTEGER CONSTANT.           
C     SEE ALSO -SQWRIT- AND -SEQADV- WHICH ARE PARTNERS TO SQREAD.      
C                                                                       
C         THE FLAG -DSKWRK- DISTINGUISHES THE PARALLEL I/O STRATEGY:    
C           IF DSKWRK=.TRUE.,  ALL PROCESSES READ THEIR OWN FILES.      
C           IF DSKWRK=.FALSE., ONLY THE MASTER WILL READ THE FILE,      
C           BUT THE DATA IS THE BROADCAST TO ALL OTHER PROCESSES.       
C                                                                       
      IF (DSKWRK.OR.MASWRK) READ(LFILE, END=200, ERR=300) REGION        
C                                                                       
C         IF RUNNING IN PARALLEL, AND THE FILE EXISTS ONLY              
C         ON THE MASTER NODE (DSKWRK=.FALSE.), THEN THE DATA            
C         SHOULD BE BROADCAST FROM THE MASTER TO ALL OTHER NODES.       
C                                                                       
      IF (GOPARR.AND.(.NOT.DSKWRK)) THEN                                
         CALL DDI_BCAST(230,'F',REGION,LENGTH,MASTER)                   
      END IF                                                            
      RETURN                                                            
C                                                                       
C                  END OF FILE                                          
C        THIS IS HANDLED BY RETURNING ZERO LENGTH READ, SO THE CALLER   
C        CAN DETERMINE IF THIS IS REALLY AN ERROR, OR WAS EXPECTED.     
C                                                                       
  200 CONTINUE                                                          
      LENGTH=0                                                          
      RETURN                                                            
C                                                                       
C                  ERROR READING FILE, PULL THE PLUG ON THE JOB         
C                                                                       
  300 CONTINUE                                                          
      WRITE(IW,9000) LFILE,ME,LENGTH                                    
      CALL ABRT                                                         
C                                                                       
      RETURN                                                            
 9000 FORMAT(1X,'ERROR READING FILE',I4,' ON NODE',I5,' LENGTH=',I10)   
      END                                                               
C*MODULE IOLIB   *DECK SQWRCH                                           
      SUBROUTINE SQWRCH(LFILE,REGION,LENGTH)                            
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      DIMENSION REGION(LENGTH)                                          
C                                                                       
C         WRITE A SEQUENTIAL RECORD TO FILE -LFILE- OF -LENGTH- WORDS,  
C         CHUNKED INTO -MAXSEQ- SUBRECORDS.                             
C         SEE ALSO THE ROUTINE SQRDCH, WHICH IS SQWRCH'S PARTNER.       
C                                                                       
      MAXSEQ = 250000                                                   
C                                                                       
      IADD = 1                                                          
  100 CONTINUE                                                          
         LEN = MIN(MAXSEQ,LENGTH-IADD+1)                                
         CALL SQWRIT(LFILE,REGION(IADD),LEN)                            
         IADD = IADD+MAXSEQ                                             
         IF(IADD.GT.LENGTH) RETURN                                      
      GO TO 100                                                         
      END                                                               
C*MODULE IOLIB   *DECK SQWRIT                                           
      SUBROUTINE SQWRIT(LFILE,REGION,LENGTH)                            
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      LOGICAL GOPARR,DSKWRK,MASWRK                                      
      DIMENSION REGION(LENGTH)                                          
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
C                                                                       
C     WRITE AN ARRAY -REGION- OF LENGTH -LENGTH- TO UNIT -LFILE-        
C     SEE ALSO -SQREAD- AND -SEQADV- WHICH ARE PARTNERS TO SQWRIT.      
C                                                                       
C         THE FLAG -DSKWRK- DISTINGUISHES THE PARALLEL I/O STRATEGY:    
C           IF DSKWRK=.TRUE.,  ALL PROCESSES WRITE THEIR OWN FILES.     
C           IF DSKWRK=.FALSE., ONLY THE MASTER WILL WRITE THE FILE.     
C                                                                       
      IF (DSKWRK.OR.MASWRK) WRITE(LFILE,ERR=300) REGION                 
      RETURN                                                            
C                                                                       
  300 CONTINUE                                                          
      WRITE(6,9020) ME,LFILE                                            
      CALL ABRT                                                         
      STOP                                                              
C                                                                       
 9020 FORMAT(1X,'SQWRIT: NODE',I4,                                      
     *          ' ENCOUNTERED I/O ERROR WRITING UNIT',I4)               
      END                                                               
C                                                                       
C*MODULE IOLIB   *DECK CCCLOS                                           
!> @date June 2015, Jerry Boatz                                         
!> - Add the code and debug output option for the 2GB record size workar
!>                                                                      
      SUBROUTINE CCCLOS(IFILE,FSTAT)                                    
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      CHARACTER*(*) FSTAT                                               
C                                                                       
      LOGICAL GOPARR,DSKWRK,MASWRK,out,dabigio,f70ok                    
      PARAMETER (mxunit=299)                                            
C                                                                       
      COMMON /CCIO  / f70ok,iosize,lrecfl(mxunit)                       
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)                
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
      COMMON /OPNNFT/ NFTOPN(MXUNIT),NODEXT(MXUNIT),IOSMP(2)            
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL                
c                                                                       
      DATA dbugccio /8HDBUGCCIO/                                        
C                                                                       
      OUT = (exetyp.eq.dbugccio .AND. MASWRK)                           
C                                                                       
C     ----- CLOSE THE SEQUENTIAL (direct access, actually?) FILE IFILE -
C     ----- FSTAT MAY BE 'KEEP' OR 'DELETE' -----                       
C                                                                       
      IF(IFILE.GT.MXUNIT) THEN                                          
         WRITE(IW,900) IFILE,MXUNIT                                     
         CALL ABRT                                                      
      END IF                                                            
      if (fstat .eq. 'DELETE') lrecfl(ifile) = 0                        
      IF(NFTOPN(IFILE).EQ.0) RETURN                                     
C                                                                       
      if (out) write(iw,*) 'DBUGCCIO:unit=',ifile,':  In CCCLOS:1, ',   
     *         'closing file  with status = ',fstat                     
      IF(DSKWRK.OR.MASWRK) CLOSE (UNIT=IFILE, STATUS=FSTAT)             
      NFTOPN(IFILE) = 0                                                 
      RETURN                                                            
  900 FORMAT(1X,'CCCLOS: ATTEMPT TO CLOSE FILE',I5,                     
     *          ' GREATER THAN MAXIMUM',I5)                             
      END                                                               
C                                                                       
C*MODULE IOLIB   *DECK CCOPEN                                           
!> @date June 2015, Jerry Boatz                                         
!> - Add the code and debug output option for the 2GB record size workar
!>                                                                      
      SUBROUTINE CCOPEN(IUNIT,IRECLN,FNAME)                             
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      LOGICAL GOPARR,DSKWRK,MASWRK,ISGDDI,wasgddi,PAROUT,INITGDDI,CHK2GB
     *        ,MLGDDI                                                   
      logical out,dabigio,f70ok,first                                   
C                                                                       
      PARAMETER (mxunit=299)                                            
      PARAMETER (itwo28=268435455)   ! = (2**28)-1                      
C                                                                       
      CHARACTER*256 FILENM,ENVBUF                                       
      COMMON /ENVIR / ENVBUF(-5:MXUNIT)                                 
      CHARACTER*(*) FNAME                                               
C                                                                       
      COMMON /CCIO  / f70ok,iosize,lrecfl(mxunit)                       
      COMMON /GDDI/   ISCOPE,NGROUPS,MYGROUP,MEGLOB,NPGLOB,NNGLOB,JBTYP,
     *                ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,NSUBGR,     
     *                MeUniv,NPUniv,numdlb,myworld,nworlds              
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)                
      COMMON /OPNNFT/ NFTOPN(MXUNIT),NODEXT(MXUNIT),IOSMP(2)            
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL                
c                                                                       
      DATA check,dbugccio /8HCHECK   ,8HDBUGCCIO/                       
      data first /.true./                                               
c                                                                       
      save first                                                        
C                                                                       
      OUT = (exetyp.eq.dbugccio .AND. MASWRK)                           
      if (out) write(iw,9000) iunit,fname,irecln                        
C                                                                       
C        ----- OPEN DIRECT ACCESS FILES USED BY CC PROGRAM -----        
C                                                                       
      IF(IUNIT.GT.MXUNIT) THEN                                          
         WRITE(IW,900) IUNIT,MXUNIT                                     
         CALL ABRT                                                      
      END IF                                                            
c-orig                                                                  
c-moved     IF(NFTOPN(IUNIT).EQ.1) RETURN                               
c-moved     NFTOPN(IUNIT) = 1                                           
c-orig                                                                  
  900 FORMAT(1X,'CCOPEN: ATTEMPT TO OPEN FILE',I5,                      
     *          ' GREATER THAN MAXIMUM',I5)                             
c                                                                       
c     --- ccdaopn2 handles files which are to have large logical records
c     --- broken up into multiple smaller physical records              
c                                                                       
c     --- Unit 70 is a special case.                                    
c     --- When used for logical file CCREST, do not use the new CC I/O r
c     --- The reason for this exception is that CCREST records may conta
c     --- real and integer data, which complicates the I/O chunking.    
c                                                                       
      f70ok = .true.                                                    
      if (iunit.eq.70 .and. fname(1:6).eq.'CCREST') f70ok=.false.       
      if (dabigio(iunit)) then                                          
         if (exetyp .eq. check) then                                    
c                                                                       
c        --- allocate memory for I/O buffer which would be used         
c        --- in a production (EXETYP=RUN) job                           
c                                                                       
c-remove-fm                                                             
c           if (first .and. maswrk) then                                
c              first = .false.                                          
c              write(iw,9110) iosize                                    
c              call getfm(iosize)                                       
c              call retfm(iosize)                                       
c           end if                                                      
c-remove-fm                                                             
         end if                                                         
         call ccdaopn2(iunit,irecln,fname)                              
      else                                                              
         lrecfl(iunit) = 1                                              
C                               ******                                  
         IF(NFTOPN(IUNIT).EQ.1) RETURN                                  
C                               ******                                  
         NFTOPN(IUNIT) = 1                                              
C                                                                       
C    THE OPEN SHELL CC CODES CAN HAVE VERY LARGE RECORD SIZES,          
C    WHICH ARE NOT HANDLED WELL BY GFORTRAN OR PATHF90 (AT LEAST).      
C    WE PRINT A WARNING THAT THE NEXT OPEN MAY FAIL.                    
C                                                                       
         IF (CHK2GB(8*IRECLN)) THEN                                     
            nchunx = int(irecln/itwo28) + 1                             
            iomax  = int(irecln/nchunx) + 1                             
            IF(MASWRK) WRITE(IW,910) IUNIT,FNAME,IRECLN                 
            IF(MASWRK) WRITE(IW,920) iomax                              
            CALL FLSHBF(IW)                                             
         END IF                                                         
  910 FORMAT(/1X,'WARNING: ATTEMPTING TO OPEN FILE',I4,' NAMED ',A/     
     *    1X,'WITH A RECORD SIZE OVER 2 GBYTES, NAMELY,',I20,' WORDS'/  
     *    1X,'YOUR FORTRAN I/O LIBRARY MAY CRASH TRYING THIS!')         
  920 FORMAT(/1X,'THIS CAN BE AVOIDED BY USING THE "IOSIZE" OPTION',/   
     *    1X,'IN $CCINP, WITH A RECOMMENDED VALUE OF (',I10,'/N + 1);', 
     *    1X,'N=1,2,3,...',/)                                           
 9000 format(1x,'DBUGCCIO:unit=',I4,':  In CCOPEN, opening direct ',    
     *          'access file with logical file name ',A,/,              
     *       1x,'and record length ',I10)                               
c-remove-fm                                                             
c9110 FORMAT(/1X,'ALLOCATING IOSIZE=',I10,' WORDS OF MEMORY FOR ',      
c    *      'CC I/O BUFFER.',/)                                         
c-remove-fm                                                             
C                                                                       
         IF (MASWRK) THEN                                               
C                                                                       
*IBM     OPEN (UNIT=IUNIT, FILE=FNAME, STATUS='UNKNOWN',                
*IBM *         ACCESS='DIRECT', FORM='UNFORMATTED', RECL=8*IRECLN)      
C                                                                       
         CALL GMS_GETENV(FNAME,FILENM)                                  
         IF(ISGDDI.or.wasgddi) THEN                                     
           FILENM=ENVBUF(IUNIT)                                         
           IF(NODEXT(IUNIT).EQ.0)                                       
     *         CALL ADDNANODE(FILENM,MeUniv,IUNIT,iout)                 
CUNX     ELSE                                                           
CUNX       CALL PARENV(FNAME,FILENM,IOUT)                               
CUNX       IF (IOUT.EQ.1) RETURN                                        
         END IF                                                         
         OPEN (UNIT=IUNIT, FILE=FILENM, STATUS='UNKNOWN',               
     *         ACCESS='DIRECT', FORM='UNFORMATTED',                     
     *         RECL=8*IRECLN)                                           
C                                                                       
*VMS     OPEN (UNIT=IUNIT, FILE=FNAME, STATUS='UNKNOWN',                
*VMS *         ACCESS='DIRECT', FORM='UNFORMATTED', RECL=2*IRECLN)      
C                                                                       
         END IF                                                         
      end if                                                            
      RETURN                                                            
      END                                                               
C*MODULE IOLIB   *DECK ccdaopn2                                         
!> @brief  This routine is called by CCOPEN to open a 'chunked' direct a
!> @author Jerry Boatz                                                  
!>  - June 2015 - subroutine written                                    
!>                                                                      
!> @details  This routine was developed as part of the workaround to    
!>           avoid file record sizes larger than 2GB.  It is a low-level
!>           counterpart to CCOPEN.  NOTE that the array lrecfl(x) is us
!>           to store the number of physical records per logical record,
!>           logical unit "x".                                          
!> @param IRAF   Direct access file logical unit number                 
!> @param LENLOG Logical record length of each direct access record     
!> @param fname  Character string containing the logical filename assign
!>                                                                      
      SUBROUTINE ccdaopn2(IRAF,lenlog,fname)                            
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      LOGICAL GOPARR,DSKWRK,MASWRK,ISGDDI,wasgddi,PAROUT,INITGDDI       
      logical dabigio,f70ok,out,mlgddi                                  
      PARAMETER (mxunit=299)                                            
      PARAMETER (hundrd=1.0e+02)                                        
C                                                                       
      character*(*) fname                                               
      CHARACTER*7 FILENM                                                
      CHARACTER*1 NULL                                                  
      CHARACTER*256 PATHNM,ENVBUF                                       
      COMMON /CCIO  / f70ok,iosize,lrecfl(mxunit)                       
      COMMON /ENVIR / ENVBUF(-5:MXUNIT)                                 
      COMMON /GDDI  / ISCOPE,NGROUPS,MYGROUP,MEGLOB,NPGLOB,NNGLOB,JBTYP,
     *                ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,NSUBGR,     
     *                MeUniv,NPUniv,numdlb,myworld,nworlds              
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)                
      COMMON /OPNNFT/ NFTOPN(MXUNIT),NODEXT(MXUNIT),IOSMP(2)            
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL                
c                                                                       
      DATA check,dbugccio /8HCHECK   ,8HDBUGCCIO/                       
C                                                                       
      OUT = (exetyp.eq.dbugccio .AND. MASWRK)                           
                                                                        
C     THIS ROUTINE, AND OTHERS IN THIS GROUP CONTROL RANDOM ACCESS      
C     files for the coupled cluster codes, mainly to avoid the          
c     problems encountered under some compilers when using record sizes 
C                                                                       
C     LOGICAL RECORDS OF LENGTH -lenlog- ARE SPREAD OVER AS             
C     MANY PHYSICAL RECORDS OF LENGTH -lenphz- AS NEEDED.               
C                                                                       
C     ----- SET THE PHYSICAL RECORD LENGTH -----                        
C                                                                       
      lenphz=iosize                                                     
      if (lenphz .gt. lenlog) lenphz=lenlog                             
c                                                                       
      intrec = lenlog/lenphz                                            
      lftovr = lenlog - (intrec * lenphz)                               
      nprlog = intrec                                                   
      nextra = 0                                                        
      if (lftovr .gt. 0) then                                           
         nprlog = nprlog + 1                                            
         nextra = iosize - lftovr                                       
      end if                                                            
      percnt = (dble(lenlog)/dble(nprlog*lenphz))*hundrd                
c                                                                       
c     --- Use lrefcl to store the number of physical records per logical
c                                                                       
      lrecfl(iraf) = nprlog                                             
      if (out) write(iw,9000) iraf,fname,lenlog,lenphz,nprlog           
c                                                                       
      FILENM='      '                                                   
      NFILNM=0                                                          
c                                                                       
      if (iraf .eq. 43) then                                            
         FILENM = 'BBBB43'                                              
         NFILNM = 6                                                     
      end if                                                            
c                                                                       
c     --- note:  unit 70 can have the logical name 'CCREST' or 'AMPROCC'
c                chunking is enabled only for the latter case..         
c                                                                       
      if (iraf .eq. 70) then                                            
         FILENM = 'AMPROCC'                                             
         NFILNM = 7                                                     
      end if                                                            
c                                                                       
c     --- note:  unit 71 can have the logical name 'CCDIIS' or 'ITOPNCC'
c                                                                       
      if (iraf .eq. 71) then                                            
         if (fname(1:6) .eq. 'CCDIIS') then                             
            FILENM = fname(1:6)                                         
            NFILNM = 6                                                  
         elseif (fname(1:7) .eq. 'ITOPNCC') then                        
            FILENM = fname(1:7)                                         
            NFILNM = 7                                                  
         end if                                                         
      end if                                                            
c                                                                       
      if (iraf .eq. 72) then                                            
         FILENM = 'CCINTS'                                              
         NFILNM = 6                                                     
      end if                                                            
c                                                                       
c     --- note:  unit 73 can have the logical name 'CCT1AMP' or 'LAMB23'
c                                                                       
      if (iraf .eq. 73) then                                            
         if (fname(1:7) .eq. 'CCT1AMP') then                            
            FILENM = fname(1:7)                                         
            NFILNM = 7                                                  
         elseif (fname(1:6) .eq. 'LAMB23') then                         
            FILENM = fname(1:6)                                         
            NFILNM = 6                                                  
         end if                                                         
      end if                                                            
c                                                                       
c     --- note:  unit 74 can have the logical name 'CCT2AMP' or 'VHHAA' 
c                                                                       
      if (iraf .eq. 74) then                                            
         if (fname(1:7) .eq. 'CCT2AMP') then                            
            FILENM = fname(1:7)                                         
            NFILNM = 7                                                  
         elseif (fname(1:5) .eq. 'VHHAA') then                          
            FILENM = fname(1:5)                                         
            NFILNM = 5                                                  
         end if                                                         
      end if                                                            
c                                                                       
      if (iraf .eq. 75) then                                            
         FILENM = 'CCT3AMP'                                             
         NFILNM = 7                                                     
      end if                                                            
c                                                                       
      if (iraf .eq. 76) then                                            
         FILENM = 'CCVM'                                                
         NFILNM = 4                                                     
      end if                                                            
c                                                                       
      if (iraf .eq. 77) then                                            
         FILENM = 'CCVE'                                                
         NFILNM = 4                                                     
      end if                                                            
c                                                                       
c     --- note:  unit 80 can have the logical name 'EOMSTAR or 'VMBA' ..
c                                                                       
      if (iraf .eq. 80) then                                            
         if (fname(1:7) .eq. 'EOMSTAR') then                            
            FILENM = fname(1:7)                                         
            NFILNM = 7                                                  
         elseif (fname(1:4) .eq. 'VMBA') then                           
            FILENM = fname(1:4)                                         
            NFILNM = 4                                                  
         end if                                                         
      end if                                                            
c                                                                       
c     --- note:  unit 82 can have the logical name 'EOMVEC2' or 'VHPRBB 
c                                                                       
      if (iraf .eq. 82) then                                            
         if (fname(1:7) .eq. 'EOMVEC2') then                            
            FILENM = fname(1:7)                                         
            NFILNM = 7                                                  
         elseif (fname(1:6) .eq. 'VHPRBB') then                         
            FILENM = fname(1:6)                                         
            NFILNM = 6                                                  
         end if                                                         
      end if                                                            
c                                                                       
c     --- note:  unit 84 can have the logical name 'EOMHC2' or 'VHPLAA' 
c                                                                       
      if (iraf .eq. 84) then                                            
         if (fname(1:6).eq.'EOMHC2' .or. fname(1:6).eq.'VHPLAA') then   
            FILENM = fname(1:6)                                         
            NFILNM = 6                                                  
         end if                                                         
      end if                                                            
c                                                                       
c     --- note:  unit 85 can have the logical name 'EOMHHHH' or 'VHPLBB'
c                                                                       
      if (iraf .eq. 85) then                                            
         if (fname(1:7) .eq. 'EOMHHHH') then                            
            FILENM = fname(1:7)                                         
            NFILNM = 7                                                  
         elseif (fname(1:6) .eq. 'VHPLBB') then                         
            FILENM = fname(1:6)                                         
            NFILNM = 6                                                  
         end if                                                         
      end if                                                            
c                                                                       
c     --- note:  unit 87 can have the logical name 'EOMRAMP' or 'VHPLBA 
c                                                                       
      if (iraf .eq. 87) then                                            
         if (fname(1:7) .eq. 'EOMRAMP') then                            
            FILENM = fname(1:7)                                         
            NFILNM = 7                                                  
         elseif (fname(1:6) .eq. 'VHPLBA') then                         
            FILENM = fname(1:6)                                         
            NFILNM = 6                                                  
         end if                                                         
      end if                                                            
c                                                                       
c     --- note:  unit 88 can have the logical name 'EOMRTMP' or 'VEAA ' 
c                                                                       
      if (iraf .eq. 88) then                                            
         if (fname(1:7) .eq. 'EOMRTMP') then                            
            FILENM = fname(1:7)                                         
            NFILNM = 7                                                  
         elseif (fname(1:4) .eq. 'VEAA') then                           
            FILENM = fname(1:4)                                         
            NFILNM = 4                                                  
         end if                                                         
      end if                                                            
c                                                                       
c     --- note:  unit 89 can have the logical name 'EOMDG12' or 'VEBB' .
c                                                                       
      if (iraf .eq. 89) then                                            
         if (fname(1:7) .eq. 'EOMDG12') then                            
            FILENM = fname(1:7)                                         
            NFILNM = 7                                                  
         elseif (fname(1:4) .eq. 'VEBB') then                           
            FILENM = fname(1:4)                                         
            NFILNM = 4                                                  
         end if                                                         
      end if                                                            
c                                                                       
c     --- note:  unit 92 can have the logical name 'MMCIVEC' or 'VPPPP  
c                                                                       
      if (iraf .eq. 92) then                                            
         if (fname(1:7) .eq. 'MMCIVEC') then                            
            FILENM = fname(1:7)                                         
            NFILNM = 7                                                  
         elseif (fname(1:5) .eq. 'VPPPP') then                          
            FILENM = fname(1:5)                                         
            NFILNM = 5                                                  
         end if                                                         
      end if                                                            
c                                                                       
c     --- note:  unit 93 can have the logical name 'MMCIVC1' or 'INTERM1
c                                                                       
      if (iraf .eq. 93) then                                            
         if (fname(1:7).eq.'MMCIVC1' .or. fname(1:7).eq.'INTERM1') then 
            FILENM = fname(1:7)                                         
            NFILNM = 7                                                  
         end if                                                         
      end if                                                            
c                                                                       
c     --- note:  unit 94 can have the logical name 'INTERM2' or 'MMCIITR
c                                                                       
      if (iraf .eq. 94) then                                            
         if (fname(1:7).eq.'INTERM2' .or. fname(1:7).eq.'MMCIITR') then 
            FILENM = fname(1:7)                                         
            NFILNM = 7                                                  
         end if                                                         
      end if                                                            
c                                                                       
c     --- note:  unit 96 can have the logical name 'EOMVL2' or 'ITSPACE'
c                                                                       
      if (iraf .eq. 96) then                                            
         if (fname(1:6) .eq. 'EOMVL2') then                             
            FILENM = fname(1:6)                                         
            NFILNM = 6                                                  
         elseif (fname(1:7) .eq. 'ITSPACE') then                        
            FILENM = fname(1:7)                                         
            NFILNM = 7                                                  
         end if                                                         
      end if                                                            
c                                                                       
c     --- note:  unit 97 can have the logical name 'EOMLVEC' or 'INSTART
c                                                                       
      if (iraf .eq. 97) then                                            
         if (fname(1:7).eq.'EOMLVEC' .or. fname(1:7).eq.'INSTART'       
     *  .or. fname(1:7).eq.'MMNREXM') then                              
            FILENM = fname(1:7)                                         
            NFILNM = 7                                                  
         end if                                                         
      end if                                                            
c                                                                       
c     --- note:  unit 98 can have the logical name 'EOMHL1' or 'ITSPC3'.
c                                                                       
      if (iraf .eq. 98) then                                            
         if (fname(1:6).eq.'EOMHL1' .or. fname(1:6).eq.'ITSPC3') then   
            FILENM = fname(1:6)                                         
            NFILNM = 6                                                  
         end if                                                         
      end if                                                            
c                                                                       
c     --- note:  unit 99 can have the logical name 'EOMHL2' or 'INSTRT3'
c                                                                       
      if (iraf .eq. 99) then                                            
         if (fname(1:6) .eq. 'EOMHL2') then                             
            FILENM = fname(1:6)                                         
            NFILNM = 6                                                  
         elseif (fname(1:7) .eq. 'INSTRT3') then                        
            FILENM = fname(1:7)                                         
            NFILNM = 7                                                  
         end if                                                         
      end if                                                            
c                                                                       
      IF(FILENM(1:1).EQ.' ') THEN                                       
         IF(MASWRK) then                                                
            WRITE(IW,*) 'ERROR IN ccdaopn2 WITH UNIT NUMBER',IRAF       
            WRITE(IW,*) 'Could not find match for filename ',fname      
         end if                                                         
         CALL ABRT                                                      
      END IF                                                            
c                                                                       
c     --- if this file is currently open already, return here before    
c     --- the actual open statement.                                    
c                                                                       
C                           ******                                      
      IF(NFTOPN(IRAF).EQ.1) RETURN                                      
C                           ******                                      
      IF (maswrk) then                                                  
         WRITE(IW,9100) FILENM(1:NFILNM),iraf,lenlog,nprlog,lenphz      
         if (nextra .gt. 0) WRITE(IW,9110) nextra,percnt                
         write(iw,*)                                                    
      end if                                                            
C                                                                       
C     ----- OPEN THE RANDOM ACCESS FILE -----                           
C                                                                       
      IF (MASWRK  .OR.  DSKWRK) THEN                                    
C                                                                       
*IBM  OPEN (UNIT=IRAF, FILE=FILENM(1:NFILNM), STATUS='UNKNOWN',         
*IBM *      ACCESS='DIRECT', FORM='UNFORMATTED', RECL=8*lenphz)         
C                                                                       
      IF(MASWRK) CALL GMS_GETENV(FILENM(1:NFILNM),PATHNM)               
      IF(ISGDDI.or.wasgddi) THEN                                        
         PATHNM=ENVBUF(IRAF)                                            
         IF(NODEXT(IRAF).EQ.0)                                          
     *       CALL ADDNANODE(PATHNM,MEGLOB,IRAF,iout)                    
      ELSE                                                              
         CALL PARENV(FILENM(1:NFILNM),PATHNM,IOUT)                      
         IF (IOUT.EQ.1) RETURN                                          
      ENDIF                                                             
      NULL = CHAR(0)                                                    
      DO 3 KOL=1,256                                                    
         IF(PATHNM(KOL:KOL).EQ.' '  .OR.                                
     *      PATHNM(KOL:KOL).EQ.NULL) GO TO 4                            
    3 CONTINUE                                                          
      KOL=257                                                           
    4 CONTINUE                                                          
      IF(KOL.EQ.1) THEN                                                 
         WRITE(IW,1) FILENM(1:NFILNM)                                   
         CALL ABRT                                                      
      END IF                                                            
      KOL=KOL-1                                                         
      OPEN (UNIT=IRAF, FILE=PATHNM(1:KOL), STATUS='UNKNOWN',            
     *      ACCESS='DIRECT', FORM='UNFORMATTED',                        
     *      RECL=8*lenphz)                                              
    1 FORMAT(1X,'YOU MUST ASSIGN GENERIC NAME ',A,' WITH A SETENV.')    
C                                                                       
*VMS  OPEN (UNIT=IRAF, FILE=FILENM(1:NFILNM), STATUS='UNKNOWN',         
*VMS *      ACCESS='DIRECT', FORM='UNFORMATTED', RECL=2*lenphz)         
C                                                                       
      END IF                                                            
      RETURN                                                            
C                                                                       
 9000 FORMAT(1X,'DBUGCCIO:unit=',I4,':  In routine ccdaopn2, ',         
     *          'opening logical file ',A,' with ',                     
     *          'logical/physical record lengths of ',I10,'/',I8,       
     *          ' and ',I8,' physical records per logical record.')     
 9100 FORMAT(1X,'OPENING LOGICAL FILE ',A8,' (UNIT',I4,'), WITH',/,     
     *       5X,'EACH LOGICAL RECORD OF ',I10,' WORDS spread across',/, 
     *       1X,I6,' PHYSICAL RECORDS OF ',I10,' WORDS EACH.')          
 9110 FORMAT(5X,'NOTE THAT ',13X,I10,' EXCESS WORDS PER LOGICAL',/,     
     *       5X,' RECORD WILL BE WRITTEN TO DISK, FOR AN',              
     *          ' EFFICIENCY OF ',F5.1,'%.')                            
      END                                                               
C*MODULE IOLIB   *DECK RAWRITES                                         
      SUBROUTINE RAWRITES(IRAF,IORA,V,LEN1,LEN,NREC,NPHYSOFF)           
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION V(*),IORA(*)                                            
      LOGICAL GOPARR,DSKWRK,MASWRK,DSKSAV,ISGDDI,wasgddi,               
     *        PAROUT,INITGDDI,MASIO,MLGDDI                              
      INTEGER DDI_WORLD                                                 
      PARAMETER (DDI_WORLD=0)                                           
      COMMON /FMCOM / X(1)                                              
      COMMON /FMOINF/ NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY              
      COMMON /FMOOPT/ ESPSCA(9),RESPAP(2),RESPPC(2),RESDIM,RESTRI(4),   
     *                RCORSD,RESPCT,CONVFG,CNVDMP,COROFF,RFLMO(4),      
     *                ORSHFT,ORSHFT2,CNVAFO,ASCREEN(4),IXESP,MXITFG,    
     *                NGUESS,NBSSE,MODORB,MODPAR,IRSTSTP,IRSTLAY,NPRFMO,
     *                NFMOPAL,MODPRP,MAXL1C,IPIEDA,MODGRD,MODESP,IVMUL, 
     *                MODLMO,NOPDEN,MOFOCK,MODFD,modfmm,ncentm,ndualb   
      COMMON /FMOPNT/ LICHFG,LMULFG,LIDMREC,LFRGNAM,LLAYFRG,LINDAT,     
     *                LNCBS,LFMOZAN,LFMOC,LFMOMAS,LIZBAS,LIAGLOB,LIBDGH,
     *                LIABDFG,LJABDFG,LNCAO,LIDXCAO,LIAPRJO,LJAPRJO,    
     *                LCOREAO,LOCCCOR,LSHIFTB,LIODFMO,LFMODA,LFMODB,    
     *                LFMOESPA,LFMOESPB,LLOCFMO,LSCFFRG,LFMOSCF,LRIJ,   
     *                LPOPMUL,LPOPMAT,LIALOC,LINDBD,LIATFRG,LINDFRG,    
     *                LINDGFRG,LNATFRG,LNAT0FRG,LIANFRG,LZANFRG,LCFRG,  
     *                LLIBISH,LLIBNSH,LLIBNG,LINDATG,LFMOBUF(3),LFMODE, 
     *                LNUMFRG,LLOCTAT,LIAOGLOB,LLOADM,LFMOGE,LDGRID,    
     *                LIODCFMO,LJOB2GRP,LFMOPG,LEMOCDR,LUNTXYZ,LUNTROT, 
     *                LSTONEP,LMAPSU,LFRGMUL,LCLMO,LIALMO,LINDLMO,      
     *                LATCLMO,LLMOBDF,LFGFLMO,LNFGLMO,LLFGLMO,LPFGLMO,  
     *                LPOPDMAT,LIDMPNT,LIDDPNT,LIVMPNT,LIACTFG,lcrfrg,  
     *                lzlmfrgv,lYlmfrgv,lndtfrg,lf_mm,lg_mm,lmaxl30,    
     *                libuffg                                           
      COMMON /FMORUN/ ESPSCF,E0SCF(2),EMP2S,IDAFMO,ICURFG,JCURFG,KCURFG,
     *                ICURLAY,ICURUNT,NAT1E,NCURSH,NGAU,ICURPOP,IFMOSTP,
     *                MONCOR,NEEDR,MODRST,NORBPROJ,NUNESP,ISKIPESP,     
     *                IESDPPC,IDOPROP,MP2RUN,ICURIT,IDMFMO,IDDFMO,      
     *                IDDCUR,NDDLEFT,IVMFMO,nzmtfmo,ifmobas,itmfmo(2)   
      COMMON /GDDI/   ISCOPE,NGROUPS,MYGROUP,MEGLOB,NPGLOB,NNGLOB,JBTYP,
     *                ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,NSUBGR,     
     *                MeUniv,NPUniv,numdlb,myworld,nworlds              
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL                
      DATA CHECK/8HCHECK    /                                           
C                                                                       
      IF(EXETYP.EQ.CHECK) RETURN                                        
c     write(6,*) 'wwwrawrites',iraf,nrec                                
C     IF(NREC.GT.NFG+1.AND.NREC.LE.NFG*2+1) WRITE(6,*) 'WRITEV',NREC    
C     IF(NREC.GT.NFG*3+1) WRITE(6,*) 'WRITEV',NREC                      
      IOPTDM=ISHFT(IAND(MODPAR,512+1024),-9)                            
      IF(IDMFMO.EQ.-1) THEN                                             
        DSKSAV=DSKWRK                                                   
        DSKWRK=IAND(MODPAR,256).EQ.0.OR.IDOPROP.EQ.1                    
        CALL RAWRITE(IRAF,IORA,V,LEN1,LEN,NREC,NPHYSOFF)                
        DSKWRK=DSKSAV                                                   
      ELSE IF(NREC.NE.1.OR.MEGLOB.EQ.0) THEN                            
        MASIO=IAND(MODPAR,2048).EQ.0                                    
C       ONLY GRAND MASTER WRITES RECORD 1 (RESTARTS)                    
C       THERE IS SOME REMOTE POSSIBILITY OF SLAVES TRYING TO READ       
C       WHAT HAS NOT BEEN SAVED BY GLOBAL MASTER YET.                   
c                                                                       
c       Right now no check if writing too long records!!                
c                                                                       
        IST=IXFTCH(X(LIDMPNT),NREC)                                     
        IF(IST.LE.0) CALL ABRT                                          
        IEND=IST+LEN-1                                                  
        IF(.NOT.MASIO.OR.MASWRK) THEN                                   
          ITMP=ISCOPE                                                   
          IF(ISGDDI) CALL GDDI_ASCOPE(DDI_WORLD)                        
C             RESDIM GRADIENT RECORDS;                                  
C             SHOULD BE CONSISTENT WITH ESDIM AND ESDGRD!!              
          IF(NREC.GT.NFG+1.AND.NREC.LE.NFG*2+1) THEN                    
C             WRITE(6,*) 'WWWAHA',IST,IEND,NREC,LEN                     
            IF(IOPTDM.EQ.1) THEN                                        
              CALL DDI_ACC(IDMFMO,1,1,IST,IEND,V)                       
c             WRITE(6,6666) 'WWWacc',IST,IEND,NREC,LEN                  
            ELSE                                                        
              CALL DDI_ACC(IDMFMO,1,LEN,IST,IST,V)                      
            ENDIF                                                       
          ELSE                                                          
            IF(IOPTDM.EQ.1) THEN                                        
              CALL DDI_PUT(IDMFMO,1,1,IST,IEND,V)                       
c             WRITE(6,6666) 'WWWput',IST,IEND,NREC,LEN                  
c6666         format(A8,4I8)                                            
            ELSE                                                        
              CALL DDI_PUT(IDMFMO,1,LEN,IST,IST,V)                      
            ENDIF                                                       
          ENDIF                                                         
          IF(ISGDDI.AND.ITMP.NE.DDI_WORLD) CALL GDDI_ASCOPE(ITMP)       
        ENDIF                                                           
c       IF(NREC.NE.1.and.MASIO) CALL DDI_SYNC(1147)                     
      ENDIF                                                             
      RETURN                                                            
      END                                                               
C*MODULE IOLIB   *DECK RAREADS                                          
      SUBROUTINE RAREADS(IRAF,IORA,V,LEN,NREC,NAV)                      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION V(*),IORA(*)                                            
      LOGICAL GOPARR,DSKWRK,MASWRK,DSKSAV,ISGDDI,PAROUT,INITGDDI,MASIO, 
     *        WASGDDI,MLGDDI                                            
      INTEGER DDI_WORLD                                                 
      PARAMETER (DDI_WORLD=0)                                           
      COMMON /FMCOM / X(1)                                              
C     COMMON /FMOINF/ NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY              
      COMMON /FMOOPT/ ESPSCA(9),RESPAP(2),RESPPC(2),RESDIM,RESTRI(4),   
     *                RCORSD,RESPCT,CONVFG,CNVDMP,COROFF,RFLMO(4),      
     *                ORSHFT,ORSHFT2,CNVAFO,ASCREEN(4),IXESP,MXITFG,    
     *                NGUESS,NBSSE,MODORB,MODPAR,IRSTSTP,IRSTLAY,NPRFMO,
     *                NFMOPAL,MODPRP,MAXL1C,IPIEDA,MODGRD,MODESP,IVMUL, 
     *                MODLMO,NOPDEN,MOFOCK,MODFD,modfmm,ncentm,ndualb   
      COMMON /FMOPNT/ LICHFG,LMULFG,LIDMREC,LFRGNAM,LLAYFRG,LINDAT,     
     *                LNCBS,LFMOZAN,LFMOC,LFMOMAS,LIZBAS,LIAGLOB,LIBDGH,
     *                LIABDFG,LJABDFG,LNCAO,LIDXCAO,LIAPRJO,LJAPRJO,    
     *                LCOREAO,LOCCCOR,LSHIFTB,LIODFMO,LFMODA,LFMODB,    
     *                LFMOESPA,LFMOESPB,LLOCFMO,LSCFFRG,LFMOSCF,LRIJ,   
     *                LPOPMUL,LPOPMAT,LIALOC,LINDBD,LIATFRG,LINDFRG,    
     *                LINDGFRG,LNATFRG,LNAT0FRG,LIANFRG,LZANFRG,LCFRG,  
     *                LLIBISH,LLIBNSH,LLIBNG,LINDATG,LFMOBUF(3),LFMODE, 
     *                LNUMFRG,LLOCTAT,LIAOGLOB,LLOADM,LFMOGE,LDGRID,    
     *                LIODCFMO,LJOB2GRP,LFMOPG,LEMOCDR,LUNTXYZ,LUNTROT, 
     *                LSTONEP,LMAPSU,LFRGMUL,LCLMO,LIALMO,LINDLMO,      
     *                LATCLMO,LLMOBDF,LFGFLMO,LNFGLMO,LLFGLMO,LPFGLMO,  
     *                LPOPDMAT,LIDMPNT,LIDDPNT,LIVMPNT,LIACTFG,lcrfrg,  
     *                lzlmfrgv,lYlmfrgv,lndtfrg,lf_mm,lg_mm,lmaxl30,    
     *                libuffg                                           
      COMMON /FMORUN/ ESPSCF,E0SCF(2),EMP2S,IDAFMO,ICURFG,JCURFG,KCURFG,
     *                ICURLAY,ICURUNT,NAT1E,NCURSH,NGAU,ICURPOP,IFMOSTP,
     *                MONCOR,NEEDR,MODRST,NORBPROJ,NUNESP,ISKIPESP,     
     *                IESDPPC,IDOPROP,MP2RUN,ICURIT,IDMFMO,IDDFMO,      
     *                IDDCUR,NDDLEFT,IVMFMO,nzmtfmo,ifmobas,itmfmo(2)   
      COMMON /GDDI/   ISCOPE,NGROUPS,MYGROUP,MEGLOB,NPGLOB,NNGLOB,JBTYP,
     *                ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,NSUBGR,     
     *                MeUniv,NPUniv,numdlb,myworld,nworlds              
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL                
      DATA CHECK/8HCHECK    /                                           
C                                                                       
      IF(EXETYP.EQ.CHECK) RETURN                                        
C     IF(NREC.GT.NFG+1.AND.NREC.LE.NFG*2+1) WRITE(6,*) 'READV',NREC     
C     IF(NREC.GT.NFG*3+1) WRITE(6,*) 'READV',NREC                       
      IOPTDM=ISHFT(IAND(MODPAR,512+1024),-9)                            
      IF(IDMFMO.EQ.-1) THEN                                             
        DSKSAV=DSKWRK                                                   
        DSKWRK=IAND(MODPAR,256).EQ.0.OR.IFMOSTP.EQ.6.OR.IFMOSTP.EQ.7    
        CALL RAREAD(IRAF,IORA,V,LEN,NREC,NAV)                           
        DSKWRK=DSKSAV                                                   
C     ELSE IF(NREC.NE.1) THEN                                           
      ELSE                                                              
        MASIO=IAND(MODPAR,2048).EQ.0                                    
C       ALL NODES NOW READ RECORD 1                                     
        IST=IXFTCH(X(LIDMPNT),NREC)                                     
        IF(IST.LE.0) CALL ABRT                                          
        IEND=IST+LEN-1                                                  
        IF(.NOT.MASIO.OR.MASWRK) THEN                                   
          ITMP=ISCOPE                                                   
          IF(ISGDDI) CALL GDDI_ASCOPE(DDI_WORLD)                        
          IF(IOPTDM.EQ.1) THEN                                          
            CALL DDI_GET(IDMFMO,1,1,IST,IEND,V)                         
c             WRITE(6,6666) 'WWWget',IST,IEND,NREC,LEN                  
c6666         format(A8,4I8)                                            
          ELSE                                                          
            CALL DDI_GET(IDMFMO,1,LEN,IST,IST,V)                        
          ENDIF                                                         
          IF(ISGDDI.AND.ITMP.NE.DDI_WORLD) CALL GDDI_ASCOPE(ITMP)       
        ENDIF                                                           
        IF(MASIO) CALL DDI_BCAST(2422,'F',V,LEN,MASTER)                 
      ENDIF                                                             
      RETURN                                                            
      END                                                               
C*MODULE IOLIB   *DECK CLOSDA                                           
      SUBROUTINE CLOSDA(FSTAT)                                          
C                                                                       
      CHARACTER*(*) FSTAT                                               
C                                                                       
      PARAMETER (MXUNIT=299)                                            
C                                                                       
      COMMON /DAIOLN/ IRECLN,IRECST,NRECUS,IFILEN(950),MERF10(950)      
      COMMON /OPNNFT/ NFTOPN(MXUNIT),NODEXT(MXUNIT),IOSMP(2)            
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)                
      COMMON /MACHSW/ KDIAG,ICORFL,IXDR,modio,mem10,lpnt10,mem10m       
C                                                                       
C     ----- CLOSE THE SEQUENTIAL FILE IFILE -----                       
C     ----- FSTAT MAY BE 'KEEP' OR 'DELETE' -----                       
C                                                                       
C     For PAROUT=.f. ON SMP MACHINES THERE IS A CONFLICT OF CLOSING DAF,
c     BECAUSE EACH NODE OPENS/CLOSES DAF THAT HAS NO                    
C     NODE-SPECIFIC EXTENSION (THUS FOR HALF NODES THE FILE IS          
C     SNEAKINGLY CLOSED BY ANOTHER HALF). THIS BOMBS OUT ON AIX.        
C     THIS IS SOLVED BY STARTING OFF WITH NGROUPS EQUAL TO THE NUMBER   
C     OF SMP BOXES, NOT THE NUMBER OF CPUS. ALTERNATIVELY, PAROUT=.t.   
C     ADDs node extensions .001 ETC.                                    
C                                                                       
      IF (NFTOPN(IDAF).NE.0.and.iand(modio,2).eq.0) THEN                
        if(mem10.eq.0) CLOSE (UNIT=IDAF, STATUS=FSTAT)                  
c       if(mem10.ne.0) write(6,*) 'close F10'                           
c       one may like to update record 1 for in-memory F10               
        NFTOPN(IDAF) = 0                                                
      ENDIF                                                             
      NRECUS=0                                                          
      RETURN                                                            
      END                                                               
C*MODULE IOLIB   *DECK ADDNANODE                                        
C>                                                                      
C>     @brief node number addition                                      
C>                                                                      
C>     @details Add node numbers to file names.                         
C>                                                                      
C>     @author Dmitri Fedorov                                           
C>                                                                      
      SUBROUTINE ADDNANODE(FILENM,ME,LUFILE,iout)                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      PARAMETER (MXUNIT=299)                                            
      CHARACTER*256 FILENM                                              
      CHARACTER*30 NODFMT                                               
      LOGICAL ISGDDI,PAROUT,INITGDDI,wasgddi,GOPARR,DSKWRK,MASWRK,MLGDDI
      COMMON /GDDI/   ISCOPE,NGROUPS,MYGROUP,MEGLOB,NPGLOB,NNGLOB,JBTYP,
     *                ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,NSUBGR,     
     *                MeUniv,NPUniv,numdlb,myworld,nworlds              
      COMMON /MACHSW/ KDIAG,ICORFL,IXDR,modio,mem10,lpnt10,mem10m       
      COMMON /OPNNFT/ NFTOPN(MXUNIT),NODEXT(MXUNIT),IOSMP(2)            
      COMMON /PAR   / ME0,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK 
C                                                                       
C     APPEND PROCESS RANK TO A FILE NAME, THIS ONLY WORKS FOR UNIX(?)   
C                                                                       
      iout=0                                                            
      if(.not.maswrk.and.(iand(modio,32).ne.0.and.LUFILE.eq.15 .or.     
     *             LUFILE.eq.4.or.LUFILE.eq.35.or.LUFILE.eq.104)) iout=1
c     4 is trajectory and 35 is restart - do not open on slaves.        
      if(LUFILE.eq.104) return                                          
c     For XYZ file do not add node extension.                           
c                                                                       
      IF(FILENM(1:1).EQ.CHAR(0)) THEN                                   
         WRITE(6,190) LUFILE                                            
         CALL ABRT                                                      
      ENDIF                                                             
  190 FORMAT(1X,'ADDNANODE: ERROR, PLEASE PRESTORE YOUR FILE NAME.'/    
     *       1X,'PROBLEM IS WITH UNIT NUMBER',I4/                       
     *       1X,'PROBABLY SOLUTION IS TO UPDATE -STORENV- ROUTINE.')    
C                                                                       
      LENSTR=256                                                        
      DO 175 I = 1,256                                                  
         IF(FILENM(I:I).EQ.' '.OR.FILENM(I:I).EQ.CHAR(0)) THEN          
            LENSTR=I-1                                                  
            GOTO 200                                                    
         ENDIF                                                          
  175 CONTINUE                                                          
  200 CONTINUE                                                          
      IF(LENSTR.GT.256-4) CALL ABRT                                     
C                                                                       
C                                                                       
C     ADDING 1.0D+00 SERVES TWO PURPOSES: TYPE CONVERSION AND ROUND-OFF 
      NFINGERS=3                                                        
      IF(MLGDDI) THEN                                                   
         IF(NPUniv.GT.999) NFINGERS=INT(LOG10(NPUniv+1.0D+00))+1        
      ELSE                                                              
         IF(NPGLOB.GT.999) NFINGERS=INT(LOG10(NPGLOB+1.0D+00))+1        
      ENDIF                                                             
C     USING I1 BELOW MEANS WE CAN ONLY USE AT MOST 10**10-1 NODES.      
C     POOR US! WHAT SHALL WE DO IF WE HAVE GOT MORE?!                   
      WRITE(UNIT=NODFMT,FMT='(6H(1H.,I,I1,1H.,I1,1H))')NFINGERS,NFINGERS
      WRITE(UNIT=FILENM(LENSTR+1:LENSTR+NFINGERS+1),FMT=NODFMT) ME      
      FILENM(LENSTR+NFINGERS+2:LENSTR+NFINGERS+2) = CHAR(0)             
C                                                                       
C     PARALLELISE I/O ON SMPS!                                          
C     SOME (INPUT) FILES ARE EXCEPTED FROM THIS: 2, 3 AND IR.           
C                                                                       
      J=IOSMP(2)                                                        
      IF(J.NE.0.AND.J.LE.LENSTR) THEN                                   
        if(iosmp(1).ge.0) then                                          
          M=MOD(MeUniv,IOSMP(1))                                        
        else                                                            
          M=MOD(-MeUniv/IOSMP(1),2)                                     
        endif                                                           
        READ(UNIT=FILENM(J:J),FMT='(I1)') ID                            
        WRITE(UNIT=FILENM(J:J),FMT='(I1)') ID+M                         
      ENDIF                                                             
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE IOLIB   *DECK STORENV                                          
C>                                                                      
C>     @brief SAVE ENVIRONMENT VARIABLE FOR RUNNING IN SUBGROUPS        
C>                                                                      
C>     @author Dmitri Fedorov                                           
C>                                                                      
      SUBROUTINE STORENV                                                
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      PARAMETER (MXUNIT=299)                                            
      CHARACTER*256 ENVBUF                                              
      LOGICAL GOPARR,DSKWRK,MASWRK                                      
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)              
      COMMON /ENVIR / ENVBUF(-5:MXUNIT)                                 
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
      COMMON /MACHIN/ NWDVAR,MAXFM,MAXSM,LIMFM,LIMSM                    
      COMMON /CIFILS/ NFT11,NFT12,NFT13,NFT14,NFT15,NFT16,IDAF20,NEMEMX 
C                                                                       
C        SAVE ENVIRONMENT VARIABLE FOR RUNNING IN SUBGROUPS.            
C                                                                       
C                THERE ARE TWO KINDS OF BASIS SET FILES:                
C        A SINGLE FILE 3 SPECIFIED BY THE END USER, THEIR OWN DATA,     
C        IS TO BE DEFINED BY -EXTBAS- (DEFAULT IS /DEV/NULL ON UNIX).   
C        TWO DIFFERENT DIRECTORIES, -MCPPATH- AND -BASPATH- CONTAIN     
C        MODEL CORE BASES AND POTENTIALS, AND ALL ELECTRON BASES.       
C        THESE ARE ALSO READ ON UNIT 3, BUT REQUIRE ADDING A SPECIFIC   
C        FILE NAME TO THE PATH NAME, AS APPROPRIATE.  HERE WE STORE     
C        THESE POSSIBLE PATH NAMES AS NEGATIVE UNIT NUMBERS, BUT OF     
C        COURSE A POSITIVE VALUE OF UNIT 3 IS USED IN THE END.          
C                                                                       
      DO I=1,MXUNIT                                                     
         ENVBUF(I)(1:1)=CHAR(0)                                         
      ENDDO                                                             
      IF (MASWRK) THEN                                                  
c        For now, replace MCQDPT as it is not used with GDDI?           
         CALL GMS_GETENV('RIVMAT', ENVBUF(51))                          
         CALL GMS_GETENV('RIT2A',  ENVBUF(52))                          
         CALL GMS_GETENV('RIT3A',  ENVBUF(53))                          
         CALL GMS_GETENV('RIT2B',  ENVBUF(54))                          
         CALL GMS_GETENV('RIT3B',  ENVBUF(55))                          
c                                                                       
         CALL GMS_GETENV('MCPPATH',ENVBUF(-5))                          
         CALL GMS_GETENV('BASPATH',ENVBUF(-4))                          
         CALL GMS_GETENV('EXTCAB', ENVBUF(-3))                          
         CALL GMS_GETENV('DFTBPAR',ENVBUF(-2))                          
C              -2 TO 0 COULD BE USED SIMILARLY TO THE ABOVE TRICKERY.   
         CALL GMS_GETENV('MAKEFP', ENVBUF(1))                           
         CALL GMS_GETENV('ERICFMT',ENVBUF(2))                           
         CALL GMS_GETENV('EXTBAS', ENVBUF(3))                           
         CALL GMS_GETENV('TRAJECT',ENVBUF(4))                           
         CALL GMS_GETENV('INPUT',  ENVBUF(IR))                          
         CALL GMS_GETENV('OUTPUT', ENVBUF(IW))                          
         CALL GMS_GETENV('PUNCH',  ENVBUF(IP))                          
         CALL GMS_GETENV('AOINTS', ENVBUF(IJK))                         
         CALL GMS_GETENV('MOINTS', ENVBUF(IJKT))                        
         CALL GMS_GETENV('DICTNRY',ENVBUF(10))                          
         CALL GMS_GETENV('DRTFILE',ENVBUF(11))                          
         CALL GMS_GETENV('CIVECTR',ENVBUF(12))                          
         CALL GMS_GETENV('CASINTS',ENVBUF(13))                          
         CALL GMS_GETENV('CIINTS', ENVBUF(14))                          
         CALL GMS_GETENV('WORK15', ENVBUF(15))                          
         CALL GMS_GETENV('WORK16', ENVBUF(16))                          
         CALL GMS_GETENV('CSFSAVE',ENVBUF(17))                          
         CALL GMS_GETENV('FOCKDER',ENVBUF(18))                          
         CALL GMS_GETENV('WORK19', ENVBUF(19))                          
         CALL GMS_GETENV('DASORT', ENVBUF(IDAF20))                      
         CALL GMS_GETENV('DFTINTS',ENVBUF(21))                          
         CALL GMS_GETENV('DFTGRID',ENVBUF(22))                          
         CALL GMS_GETENV('JKFILE', ENVBUF(23))                          
         CALL GMS_GETENV('ORDINT', ENVBUF(24))                          
         CALL GMS_GETENV('EFPIND', ENVBUF(25))                          
         CALL GMS_GETENV('PCMDATA',ENVBUF(26))                          
         CALL GMS_GETENV('PCMINTS',ENVBUF(27))                          
         CALL GMS_GETENV('MLTPL',  ENVBUF(28))                          
         CALL GMS_GETENV('MLTPLT', ENVBUF(29))                          
         CALL GMS_GETENV('DAFL30', ENVBUF(30))                          
         CALL GMS_GETENV('RESTART',ENVBUF(35))                          
         CALL GMS_GETENV('HESSIAN',ENVBUF(38))                          
         CALL GMS_GETENV('SOCCDAT',ENVBUF(40))                          
         CALL GMS_GETENV('AABB41', ENVBUF(41))                          
         CALL GMS_GETENV('BBAA42', ENVBUF(42))                          
         CALL GMS_GETENV('BBBB43', ENVBUF(43))                          
         CALL GMS_GETENV('REMD',   ENVBUF(44))                          
         CALL GMS_GETENV('MCQD50', ENVBUF(50))                          
c        CALL GMS_GETENV('MCQD51', ENVBUF(51))                          
c        CALL GMS_GETENV('MCQD52', ENVBUF(52))                          
c        CALL GMS_GETENV('MCQD53', ENVBUF(53))                          
c        CALL GMS_GETENV('MCQD54', ENVBUF(54))                          
c        CALL GMS_GETENV('MCQD55', ENVBUF(55))                          
         CALL GMS_GETENV('MCQD56', ENVBUF(56))                          
         CALL GMS_GETENV('MCQD57', ENVBUF(57))                          
         CALL GMS_GETENV('MCQD58', ENVBUF(58))                          
         CALL GMS_GETENV('MCQD59', ENVBUF(59))                          
         CALL GMS_GETENV('MCQD60', ENVBUF(60))                          
         CALL GMS_GETENV('MCQD61', ENVBUF(61))                          
         CALL GMS_GETENV('MCQD62', ENVBUF(62))                          
         CALL GMS_GETENV('MCQD63', ENVBUF(63))                          
         CALL GMS_GETENV('MCQD64', ENVBUF(64))                          
         CALL GMS_GETENV('DCPHFH2',ENVBUF(67))                          
         CALL GMS_GETENV('NMRINT1',ENVBUF(61))                          
         CALL GMS_GETENV('CCREST', ENVBUF(70))                          
         CALL GMS_GETENV('CCDIIS', ENVBUF(71))                          
         CALL GMS_GETENV('CCINTS', ENVBUF(72))                          
         CALL GMS_GETENV('CCT1AMP',ENVBUF(73))                          
         CALL GMS_GETENV('CCT2AMP',ENVBUF(74))                          
         CALL GMS_GETENV('CCT3AMP',ENVBUF(75))                          
         CALL GMS_GETENV('CCVM',   ENVBUF(76))                          
         CALL GMS_GETENV('CCVE',   ENVBUF(77))                          
         CALL GMS_GETENV('CCQUADS',ENVBUF(78))                          
         CALL GMS_GETENV('QUADSVO',ENVBUF(79))                          
         CALL GMS_GETENV('EOMSTAR',ENVBUF(80))                          
         CALL GMS_GETENV('EOMVEC1',ENVBUF(81))                          
         CALL GMS_GETENV('EOMVEC2',ENVBUF(82))                          
         CALL GMS_GETENV('EOMHC1', ENVBUF(83))                          
         CALL GMS_GETENV('EOMHC2', ENVBUF(84))                          
         CALL GMS_GETENV('EOMHHHH',ENVBUF(85))                          
         CALL GMS_GETENV('EOMPPPP',ENVBUF(86))                          
         CALL GMS_GETENV('EOMRAMP',ENVBUF(87))                          
         CALL GMS_GETENV('EOMRTMP',ENVBUF(88))                          
         CALL GMS_GETENV('EOMDG12',ENVBUF(89))                          
         CALL GMS_GETENV('MMPP',   ENVBUF(90))                          
         CALL GMS_GETENV('MMHPP',  ENVBUF(91))                          
         CALL GMS_GETENV('MMCIVEC',ENVBUF(92))                          
         CALL GMS_GETENV('MMCIVC1',ENVBUF(93))                          
         CALL GMS_GETENV('MMCIITR',ENVBUF(94))                          
         CALL GMS_GETENV('EOMVL1', ENVBUF(95))                          
         CALL GMS_GETENV('EOMVL2', ENVBUF(96))                          
         CALL GMS_GETENV('EOMLVEC',ENVBUF(97))                          
         CALL GMS_GETENV('EOMHL1', ENVBUF(98))                          
         CALL GMS_GETENV('EOMHL2', ENVBUF(99))                          
         CALL GMS_GETENV('EFMOI', ENVBUF(102))                          
         CALL GMS_GETENV('EFMOF', ENVBUF(103))                          
         CALL GMS_GETENV('XYZ', ENVBUF(104))                            
c        The four files below replace those above.                      
c        It is a mistake, each file should have a unique number.        
c        CALL GMS_GETENV('DEN2P1',ENVBUF(70))                           
c        CALL GMS_GETENV('DEN2P2',ENVBUF(71))                           
c        CALL GMS_GETENV('DEN2P3',ENVBUF(72))                           
c        CALL GMS_GETENV('DEN2P4',ENVBUF(73))                           
      END IF                                                            
C                                                                       
C     A LITTLE NASTY TRICK TO BROADCAST CHARACTERS: PRETEND THEY ARE    
C     INTEGERS.  IT IS NOT SAFE TO PRETEND THEY ARE REAL DUE TO SOME    
C     "MAGIC" COMBINATIONS THAT CAN CAUSE SOME CAUTIOUS COMPILERS TO    
C     REPORT FPU ERRORS.                                                
C     DIMENSIONS ASSUME THAT A REAL VARIABLE HOLDS 8 BYTES WHICH SO FAR 
C     HAS ALWAYS BEEN TRUE. THAT IS WHY THERE IS AN "8" BELOW.          
C     IT DOES HELP THAT 256 IS A MULTIPLE OF 8!                         
C     A WORD OF CAUTION: THIS SWINDLING MAY PROVE TO BE UNPORTABLE!     
C     THE REASON IS THAT SOME COMPILERS MAY SECRETLY STORE THE STRING   
C     SIZE.                                                             
C                                                                       
      CALL DDI_BCAST(206,'I',ENVBUF,((MXUNIT*256)/8)*NWDVAR,MASTER)     
      RETURN                                                            
      END                                                               
C                                                                       
C*MODULE IOLIB   *DECK GMS_GETENV                                       
      SUBROUTINE GMS_GETENV(NAME,VALUE)                                 
      CHARACTER*(*) NAME,VALUE                                          
C                                                                       
C        THIS IS JUST A SILLY WRAPPER AROUND THE SYSTEM'S -GETENV-      
C                                                                       
C        HOWEVER, A SED HACK WOULD ALLOW US TO CALL SOME                
C        OTHER NAME, INSTEAD OF GETENV, IF NEEDS MUST.                  
C                                                                       
C        AS AN EXAMPLE, SEE THE FILE_GETENV CODE IN UNPORT.SRC.         
C                                                                       
      LENNAM = LEN(NAME)                                                
      IF(LENNAM.LE.0) THEN                                              
         WRITE(6,*) '--ERROR CALLING GMS_GETENV--'                      
         CALL ABRT                                                      
      END IF                                                            
C                                                                       
      CALL GETENV(NAME,VALUE)                                           
      RETURN                                                            
      END                                                               
C                                                                       
C*MODULE IOLIB   *DECK MPOPEN                                           
      SUBROUTINE MPOPEN(IUNIT,IRECLN,FNAME)                             
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      LOGICAL GOPARR,DSKWRK,MASWRK,ISGDDI,wasgddi,PAROUT,INITGDDI,MLGDDI
C                                                                       
      PARAMETER (MXUNIT=299)                                            
C                                                                       
      CHARACTER*256 FILENM,ENVBUF                                       
      COMMON /ENVIR / ENVBUF(-5:MXUNIT)                                 
      CHARACTER*1 NULL                                                  
      CHARACTER*(*) FNAME                                               
C                                                                       
      COMMON /GDDI/   ISCOPE,NGROUPS,MYGROUP,MEGLOB,NPGLOB,NNGLOB,JBTYP,
     *                ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,NSUBGR,     
     *                MeUniv,NPUniv,numdlb,myworld,nworlds              
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)                
      COMMON /OPNNFT/ NFTOPN(MXUNIT),NODEXT(MXUNIT),IOSMP(2)            
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
C                                                                       
C        ----- OPEN DIRECT ACCESS FILES USED BY IMS MP2 PROGRAM -----   
C                                                                       
      IF(IUNIT.GT.MXUNIT) THEN                                          
         WRITE(IW,900) IUNIT,MXUNIT                                     
         CALL ABRT                                                      
      END IF                                                            
      IF(NFTOPN(IUNIT).EQ.1) RETURN                                     
      NFTOPN(IUNIT) = 1                                                 
  900 FORMAT(1X,'MPOPEN: ATTEMPT TO OPEN FILE',I5,                      
     *          ' GREATER THAN MAXIMUM',I5)                             
C                                                                       
*IBM  IF (MASWRK)                                                       
*IBM *   OPEN (UNIT=IUNIT, FILE=FNAME, STATUS='UNKNOWN',                
*IBM *         ACCESS='DIRECT', FORM='UNFORMATTED', RECL=8*IRECLN)      
C                                                                       
      FILENM=' '                                                        
      IF(MASWRK) CALL GETENV(FNAME,FILENM)                              
      IF(ISGDDI.or.wasgddi) THEN                                        
        FILENM=ENVBUF(IUNIT)                                            
        IF(NODEXT(IUNIT).EQ.0)                                          
     *        CALL ADDNANODE(FILENM,MeUniv,IUNIT,iout)                  
      ELSE                                                              
        CALL PARENV(FNAME,FILENM,IOUT)                                  
        IF (IOUT.EQ.1) RETURN                                           
      END IF                                                            
      NULL = CHAR(0)                                                    
      DO 3 KOL=1,256                                                    
         IF(FILENM(KOL:KOL).EQ.' '  .OR.                                
     *      FILENM(KOL:KOL).EQ.NULL) GO TO 4                            
    3 CONTINUE                                                          
      KOL=257                                                           
    4 CONTINUE                                                          
      IF(KOL.EQ.1) THEN                                                 
         WRITE(IW,1) FNAME                                              
         CALL ABRT                                                      
      END IF                                                            
    1 FORMAT(1X,'YOU MUST ASSIGN GENERIC NAME ',A,' WITH A SETENV.')    
      KOL=KOL-1                                                         
      OPEN (UNIT=IUNIT, FILE=FILENM(1:KOL), STATUS='UNKNOWN',           
     *      ACCESS='DIRECT', FORM='UNFORMATTED',                        
     *      RECL=8*IRECLN)                                              
C                                                                       
*VMS  IF (MASWRK)                                                       
*VMS *    OPEN (UNIT=IUNIT, FILE=FNAME, STATUS='UNKNOWN',               
*VMS *       ACCESS='DIRECT', FORM='UNFORMATTED', RECL=2*IRECLN)        
C                                                                       
      RETURN                                                            
      END                                                               
C                                                                       
C                                                                       
C*MODULE IOLIB   *DECK MPCLOS                                           
      SUBROUTINE MPCLOS(IFILE,FSTAT)                                    
C                                                                       
      CHARACTER*(*) FSTAT                                               
C                                                                       
      LOGICAL GOPARR,DSKWRK,MASWRK                                      
      PARAMETER (MXUNIT=299)                                            
C                                                                       
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)                
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
      COMMON /OPNNFT/ NFTOPN(MXUNIT),NODEXT(MXUNIT),IOSMP(2)            
C                                                                       
C     ----- CLOSE THE DIRECT ACCESS FILE IFILE -----                    
C     ----- FSTAT MAY BE 'KEEP' OR 'DELETE' -----                       
C                                                                       
      IF(IFILE.GT.MXUNIT) THEN                                          
         WRITE(IW,900) IFILE,MXUNIT                                     
         CALL ABRT                                                      
      END IF                                                            
      IF(NFTOPN(IFILE).EQ.0) RETURN                                     
C                                                                       
      IF(DSKWRK.OR.MASWRK) CLOSE (UNIT=IFILE, STATUS=FSTAT)             
      NFTOPN(IFILE) = 0                                                 
      RETURN                                                            
  900 FORMAT(1X,'MPCLOS: ATTEMPT TO CLOSE FILE',I5,                     
     *          ' GREATER THAN MAXIMUM',I5)                             
      END                                                               
C*MODULE IOLIB   *DECK RAOPDC                                           
      SUBROUTINE RAOPDC(IRAF,IORA,NUMREC,NPRINT)                        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      LOGICAL GOPARR,DSKWRK,MASWRK,ISGDDI,wasgddi,PAROUT,INITGDDI,MLGDDI
C                                                                       
      PARAMETER (MXUNIT=299)                                            
C                                                                       
      CHARACTER*6 FILENM                                                
      CHARACTER*1 NULL                                                  
      CHARACTER*256 PATHNM,ENVBUF                                       
      COMMON /ENVIR / ENVBUF(-5:MXUNIT)                                 
C                                                                       
      COMMON /RAIODC/ IRECLN,IRECST                                     
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAFX,NAVX,IODAX(950)             
      COMMON /OPNNFT/ NFTOPN(MXUNIT),NODEXT(MXUNIT),IOSMP(2)            
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL                
      COMMON /GDDI/   ISCOPE,NGROUPS,MYGROUP,MEGLOB,NPGLOB,NNGLOB,JBTYP,
     *                ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,NSUBGR,     
     *                MeUniv,NPUniv,numdlb,myworld,nworlds              
C                                                                       
      DATA CHECK/8HCHECK    /                                           
C                                                                       
      DIMENSION IORA(NUMREC)                                            
C                                                                       
      IF(IRAF.GT.MXUNIT) THEN                                           
         WRITE(IW,900) IRAF,MXUNIT                                      
         CALL ABRT                                                      
      END IF                                                            
      IF(NFTOPN(IRAF).EQ.1) RETURN                                      
      NFTOPN(IRAF) = 1                                                  
  900 FORMAT(1X,'RAOPDC: ATTEMPT TO OPEN FILE',I5,                      
     *          ' GREATER THAN MAXIMUM',I5)                             
C                                                                       
C     ----- SET THE PHYSICAL RECORD LENGTH -----                        
C                                                                       
      IRECLN=NRASIZ(IRAF)                                               
      IF(EXETYP.EQ.CHECK) RETURN                                        
C                                                                       
C     ----- OPEN THE RANDOM ACCESS FILE -----                           
C                                                                       
C        BY CONVENTION, FILES ENDING IN 0'S ARE DIRECT ACCESS,          
C        AT PRESENT ONLY FILE 20 AND 30 ARE USED.                       
C                                                                       
      FILENM='      '                                                   
      IF(IRAF.EQ.20) FILENM = 'DASORT'                                  
      IF(IRAF.EQ.30) FILENM = 'DAFL30'                                  
      IF(IRAF.EQ.40) FILENM = 'FMODAT'                                  
      IF(FILENM.EQ.'      ') THEN                                       
         IF(MASWRK) WRITE(IW,*) 'ERROR IN RAOPEN WITH UNIT NUMBER',IRAF 
         CALL ABRT                                                      
      END IF                                                            
C                                                                       
      IF(NPRINT.NE.-5  .AND.  MASWRK)                                   
     *     WRITE(IW,9000) FILENM,NUMREC,IRECLN                          
C                                                                       
      IF (MASWRK  .OR.  DSKWRK) THEN                                    
C                                                                       
*IBM  OPEN (UNIT=IRAF, FILE=FILENM, STATUS='UNKNOWN',                   
*IBM *      ACCESS='DIRECT', FORM='UNFORMATTED', RECL=8*IRECLN)         
C                                                                       
      IF(MASWRK) CALL GMS_GETENV(FILENM,PATHNM)                         
      IF(ISGDDI.or.wasgddi) THEN                                        
         PATHNM=ENVBUF(IRAF)                                            
         IF(NODEXT(IRAF).EQ.0)                                          
     *       CALL ADDNANODE(PATHNM,MeUniv,IRAF,iout)                    
      ELSE                                                              
         CALL PARENV(FILENM,PATHNM,IOUT)                                
         IF (IOUT.EQ.1) RETURN                                          
      ENDIF                                                             
      NULL = CHAR(0)                                                    
      DO 3 KOL=1,256                                                    
         IF(PATHNM(KOL:KOL).EQ.' '  .OR.                                
     *      PATHNM(KOL:KOL).EQ.NULL) GO TO 4                            
    3 CONTINUE                                                          
      KOL=257                                                           
    4 CONTINUE                                                          
      IF(KOL.EQ.1) THEN                                                 
         WRITE(IW,1) FILENM                                             
         CALL ABRT                                                      
      END IF                                                            
      KOL=KOL-1                                                         
      OPEN (UNIT=IRAF, FILE=PATHNM(1:KOL), STATUS='UNKNOWN',            
     *      ACCESS='DIRECT', FORM='UNFORMATTED',                        
     *      RECL=8*IRECLN)                                              
    1 FORMAT(1X,'YOU MUST ASSIGN GENERIC NAME ',A,' WITH A SETENV.')    
C                                                                       
*VMS  OPEN (UNIT=IRAF, FILE=FILENM, STATUS='UNKNOWN',                   
*VMS *      ACCESS='DIRECT', FORM='UNFORMATTED', RECL=2*IRECLN)         
C                                                                       
      END IF                                                            
C                                                                       
C      ----- INITIALIZE FOR EMPTY SORTFILE ----                         
C                                                                       
      IRECST = 1                                                        
      DO 100 I = 1,NUMREC                                               
         IORA(I) = -1                                                   
  100 CONTINUE                                                          
      RETURN                                                            
C                                                                       
 9000 FORMAT(1X,'OPENING FILE ',A6,' WITH',I8,' LOGICAL RECORDS.',      
     *          /1X,'IRECLN=',I8,' WORDS')                              
      END                                                               
C*MODULE IOLIB   *DECK RAREDC                                           
      SUBROUTINE RAREDC(IRAF,IORA,V,LEN,NREC,NAVM)                      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      LOGICAL GOPARR,DSKWRK,MASWRK                                      
C                                                                       
      DIMENSION V(*),IORA(*)                                            
C                                                                       
      COMMON /RAIODC/ IRECLN,IRECST                                     
      COMMON /IOFILE/ IR,IW,IP,IIS,IPK,IDAFX,NAVX,IODAX(950)            
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
C                                                                       
C         READ A LOGICAL RECORD FROM THE DAF DICTIONARY FILE            
C         A LOGICAL RECORD MAY SPAN SEVERAL PHYSICAL RECORDS.           
C                                                                       
C         CALLING ARGUMENT -IDTYP- IS 0 OR 1 FOR FLOATING POINT         
C         OR INTEGER RECORDS.  NO OTHER DATA TYPE IS ALLOWED!           
C         RECORDS MUST BE PURELY FLOATING POINT OR PURELY INTEGER.      
C         NO MATTER WHAT -IDTYP- IS, THE -LEN- OF THE RECORD            
C         MUST BE GIVEN IN TERMS OF FLOATING POINT WORDS.  IN           
C         TURN, THIS MEANS THAT INTEGER RECORDS ON 32 BIT MACHINES      
C         MUST CONTAIN AN EVEN NUMBER OF INTEGERS.                      
C                                                                       
      N = IORA(NREC)                                                    
      IF(N.EQ.-1) GO TO 800                                             
      IS = -IRECLN + 1                                                  
      NS = N                                                            
      LENT = LEN                                                        
  100 CONTINUE                                                          
         IS = IS + IRECLN                                               
         IF = IS + LENT - 1                                             
         IF ((IF-IS+1) .GT. IRECLN) IF = IS + IRECLN - 1                
         NSP = NS                                                       
         LENW = IF - IS + 1                                             
         CALL DARD(V(IS),LENW,IRAF,NSP,NAVM)                            
         LENT = LENT - IRECLN                                           
         NS = NS + 1                                                    
         N = NS                                                         
      IF (LENT .GE. 1) GO TO 100                                        
      RETURN                                                            
C                                                                       
  800 CONTINUE                                                          
      IF (MASWRK) WRITE(IW,9000) NREC,LEN                               
      CALL ABRT                                                         
      RETURN                                                            
C                                                                       
 9000 FORMAT(1X,'*** ERROR ***, ATTEMPT TO READ A -DASORT- RECORD',     
     *         ' THAT WAS NEVER WRITTEN.'/1X,'NREC,LEN=',I5,I10)        
      END                                                               
C*MODULE IOLIB   *DECK RAWRDC                                           
      SUBROUTINE RAWRDC(IRAF,IORA,IRALEN,V,LEN,NREC,NAVM)               
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      LOGICAL GOPARR,DSKWRK,MASWRK                                      
C                                                                       
      DIMENSION V(*),IORA(*),IRALEN(*)                                  
C                                                                       
      COMMON /RAIODC/ IRECLN,IRECST                                     
      COMMON /IOFILE/ IR,IW,IP,IIS,IPK,IDAFX,NAVX,IODAX(950)            
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
C                                                                       
C         WRITE A LOGICAL RECORD ON THE DAF DICTIONARY FILE             
C         A LOGICAL RECORD MAY SPAN SEVERAL PHYSICAL RECORDS            
C                                                                       
      IF(IRAF.LT.0) WRITE(6,*) 'BOGUS RAWRIT, NAVM=',NAVM               
C                                                                       
      N = IORA(NREC)                                                    
      IF (N .GT. 0 .AND. LEN .NE. IRALEN(NREC)) GO TO 800               
      IF (N .GT. 0) GO TO 100                                           
      IORA(NREC) = IRECST                                               
      IRALEN(NREC) = LEN                                                
      IRECST = IRECST + (LEN-1)/IRECLN + 1                              
      N = IORA(NREC)                                                    
  100 CONTINUE                                                          
      IST = -IRECLN + 1                                                 
      NS = N                                                            
      LENT = LEN                                                        
  120 CONTINUE                                                          
         IST = IST + IRECLN                                             
         IF = IST + LENT - 1                                            
         IF ((IF-IST+1) .GT. IRECLN) IF = IST+IRECLN-1                  
         NSP = NS                                                       
         LENW = IF - IST + 1                                            
         CALL DAWRT(V(IST),LENW,IRAF,NSP)                               
         LENT = LENT - IRECLN                                           
         NS = NS + 1                                                    
         N = NS                                                         
      IF (LENT .GE. 1) GO TO 120                                        
      RETURN                                                            
C                                                                       
  800 CONTINUE                                                          
      IF (MASWRK) WRITE (IW,9008) NREC,LEN,IRALEN(NREC)                 
      CALL ABRT                                                         
      RETURN                                                            
C                                                                       
 9008 FORMAT(1X,'RAWRDC HAS REQUESTED A RECORD WITH LENGTH',            
     *       1X,'DIFFERENT THAN BEFORE - ABORT FORCED.'/                
     *       1X,'RECORD ',I5,' NEW LENGTH =',I5,                        
     *          ' OLD LENGTH =',I5)                                     
      END                                                               
C*MODULE IOLIB   *DECK DIROPN                                           
      SUBROUTINE DIROPN(IUNIT,FNAME,FSTAT,IRECLN)                       
C                                                                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C                                                                       
      CHARACTER*(*) FNAME,FSTAT                                         
      LOGICAL GOPARR,DSKWRK,MASWRK,ISGDDI,PAROUT,INITGDDI,wasgddi,      
     *        MLGDDI                                                    
      PARAMETER (MXUNIT=299)                                            
C                                                                       
      CHARACTER*256 FILENM,ENVBUF                                       
      COMMON /ENVIR / ENVBUF(-5:MXUNIT)                                 
      CHARACTER*1 NULL                                                  
C                                                                       
      COMMON /GDDI/   ISCOPE,NGROUPS,MYGROUP,MEGLOB,NPGLOB,NNGLOB,JBTYP,
     *                ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,NSUBGR,     
     *                MeUniv,NPUniv,numdlb,myworld,nworlds              
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)                
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
      COMMON /OPNNFT/ NFTOPN(MXUNIT),NODEXT(MXUNIT),IOSMP(2)            
C                                                                       
C      OPEN DIRECT ACCESS FILE -IUNIT- WITH GENERIC NAME -FNAME-,       
C      WITH STATUS -FSTAT-, AND WITH RECORD LENGTH -IRECLN-.            
C      THE LATTER QUANTITY SHOULD BE MEASURED IN 64 BIT WORDS.          
C                                                                       
      IF (.NOT.MASWRK  .AND.  .NOT.DSKWRK) RETURN                       
C                                                                       
      IF(NFTOPN(IUNIT).EQ.1) RETURN                                     
      NFTOPN(IUNIT) = 1                                                 
C                                                                       
      MAXFILES = MXUNIT  ! TEMP VALUE SUPPRESSES A DIAGNOSTIC ON AXP64  
      IF(IUNIT.GT.MAXFILES) THEN                                        
         WRITE(IW,900) IUNIT,MAXFILES                                   
         IF(MAXFILES.EQ.99) WRITE(IW,901)                               
         CALL ABRT                                                      
      END IF                                                            
  900 FORMAT(1X,'DIROPN: ATTEMPT TO OPEN FILE',I5,                      
     *          ' GREATER THAN MAXIMUM',I5)                             
  901 FORMAT(1X,'G77 HAS AN INTERNAL LIMIT OF 99 FILES,',               
     *          ' AND SO CANNOT EXECUTE RUNTYP=TDHFX')                  
C                                                                       
*IBM  OMITTED.                                                          
C                                                                       
C         FIRST, TRANSLATE THE ENVIRONMENT VARIABLE TO A FULLY          
C                QUALIFIED FILE NAME, FNAME --> FILENM.                 
C          NEXT, APPEND AN EXTENSION TO THE FILE NAME, GIVING ITS       
C                RANK, SUCH AS .001, ON ALL PROCESSES WHICH ARE NOT     
C                THE MASTER.  THUS SMP SYSTEMS HAVE UNIQUE FILE NAMES.  
C         THIRD, STICK IN A NULL BYTE SO THAT UNIX COMPILERS, WHICH     
C                DON'T UNDERSTAND THE FORTRAN CHARACTER TYPE, HAVE      
C                A C-STYLE CHAR VARIABLE TO WORK WITH (UGH).            
C   AND FINALLY, OPEN THE FILE.                                         
C                                                                       
      FILENM=' '                                                        
      IF (MASWRK) CALL GMS_GETENV(FNAME,FILENM)                         
      IF(ISGDDI.or.wasgddi) THEN                                        
        FILENM=ENVBUF(IUNIT)                                            
        IF(NODEXT(IUNIT).EQ.0)                                          
     *     CALL ADDNANODE(FILENM,MeUniv,IUNIT,iout)                     
      ELSE                                                              
        CALL PARENV(FNAME,FILENM,IOUT)                                  
        IF (IOUT.EQ.1) RETURN                                           
      END IF                                                            
      NULL = CHAR(0)                                                    
      DO 1 KOL=1,256                                                    
         IF(FILENM(KOL:KOL).EQ.' '  .OR.                                
     *      FILENM(KOL:KOL).EQ.NULL) GO TO 2                            
    1 CONTINUE                                                          
      KOL=257                                                           
    2 CONTINUE                                                          
      IF(KOL.EQ.1) THEN                                                 
         WRITE(IW,3) FNAME                                              
         CALL ABRT                                                      
      END IF                                                            
      KOL=KOL-1                                                         
      OPEN(UNIT=IUNIT, FILE=FILENM(1:KOL), STATUS=FSTAT,                
     *      ACCESS='DIRECT', FORM='UNFORMATTED', ERR=4,                 
     *      RECL=8*IRECLN)                                              
      RETURN                                                            
C                                                                       
C         ERROR HANDLING (E.G. FOR ERICFMT, ...)                        
C                                                                       
    3 FORMAT(1X,'YOU MUST ASSIGN GENERIC NAME ',A,' WITH A SETENV.')    
C                                                                       
    4 CONTINUE                                                          
      IF(FSTAT.EQ.'OLD') THEN                                           
         IF(MASWRK) WRITE(IW,5) 'PRE-EXISTING',FNAME,FILENM(1:KOL)      
      ELSE                                                              
         IF(MASWRK) WRITE(IW,5) 'NEW',FNAME,FILENM(1:KOL)               
      ENDIF                                                             
      CALL ABRT                                                         
    5 FORMAT(//1X,'ERROR OPENING ',A,' FILE ',A,','/                    
     *         1X,'ASSIGNED TO EXPLICIT FILE NAME ',A,','/              
     *         1X,'PLEASE CHECK THE -SETENV- FILE ASSIGNMENTS',         
     *            ' IN YOUR -RUNGMS- SCRIPT.')                          
C                                                                       
*VMS  OMITTED.                                                          
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE IOLIB   *DECK DIRCLO                                           
      SUBROUTINE DIRCLO(IFILE,FSTAT)                                    
C                                                                       
      CHARACTER*(*) FSTAT                                               
C                                                                       
      LOGICAL GOPARR,DSKWRK,MASWRK                                      
      PARAMETER (MXUNIT=299)                                            
C                                                                       
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)                
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
      COMMON /OPNNFT/ NFTOPN(MXUNIT),NODEXT(MXUNIT),IOSMP(2)            
C                                                                       
C     ----- CLOSE THE DIRECT ACCESS FILE IFILE -----                    
C     ----- FSTAT MAY BE 'KEEP' OR 'DELETE' -----                       
C                                                                       
      IF(IFILE.GT.MXUNIT) THEN                                          
         WRITE(IW,900) IFILE,MXUNIT                                     
         CALL ABRT                                                      
      END IF                                                            
      IF(NFTOPN(IFILE).EQ.0) RETURN                                     
C                                                                       
      IF(DSKWRK.OR.MASWRK) CLOSE (UNIT=IFILE, STATUS=FSTAT)             
      NFTOPN(IFILE) = 0                                                 
      RETURN                                                            
  900 FORMAT(1X,'DIRCLO: ATTEMPT TO CLOSE FILE',I5,                     
     *          ' GREATER THAN MAXIMUM',I5)                             
      END                                                               
C*MODULE IOLIB   *DECK OPENDAC                                          
!> @brief Clone of OPENDA to control the CIMDCT file                    
!>                                                                      
!> @param IREST the restart value                                       
!> @author Albert DeFusco                                               
!> @date September 2010                                                 
      SUBROUTINE OPENDAC(IREST)                                         
C                                                                       
C     - - - - OPEN CIMSUB DICTIONARY FILE 194 - - - -                   
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      LOGICAL GOPARR,DSKWRK,MASWRK,ISGDDI,wasgddi,PAROUT,INITGDDI,MLGDDI
C                                                                       
      PARAMETER (MXUNIT=299)                                            
      CHARACTER*256 FILENM,ENVBUF                                       
      COMMON /ENVIR / ENVBUF(-5:MXUNIT)                                 
C                                                                       
      COMMON /DAIOLN/ IRECLN,IRECST,NRECUS,IFILEN(950),MERF10(950)      
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)                
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
      COMMON /GDDI/   ISCOPE,NGROUPS,MYGROUP,MEGLOB,NPGLOB,NNGLOB,JBTYP,
     *                ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,NSUBGR,     
     *                MeUniv,NPUniv,numdlb,myworld,nworlds              
      COMMON /OPNNFT/ NFTOPN(MXUNIT),NODEXT(MXUNIT),IOSMP(2)            
C                                                                       
      IDAF = 194                                                        
      IRECLN=NRASIZ(IDAF)                                               
C                                                                       
      IF (MASWRK) THEN                                                  
      NFTOPN(IDAF) =1                                                   
C                                                                       
*IBM  OPEN (UNIT=IDAF, FILE='CIMDCT', STATUS='UNKNOWN',                 
*IBM *      ACCESS='DIRECT', FORM='UNFORMATTED', RECL=8*IRECLN)         
C                                                                       
      IF(ISGDDI.or.wasgddi) THEN                                        
         FILENM=ENVBUF(IDAF)                                            
         IF(PAROUT.AND.NODEXT(IDAF).EQ.0)                               
     *       CALL ADDNANODE(FILENM,MeUniv,IDAF,iout)                    
      ELSE                                                              
         CALL GMS_GETENV('CIMDCT',FILENM)                               
      ENDIF                                                             
      OPEN (UNIT=IDAF, FILE=FILENM, STATUS='UNKNOWN',                   
     *      ACCESS='DIRECT', FORM='UNFORMATTED',                        
     *      RECL=8*IRECLN)                                              
C                                                                       
*VMS  OPEN (UNIT=IDAF, FILE='CIMDCT', STATUS='UNKNOWN',                 
*VMS *      ACCESS='DIRECT', FORM='UNFORMATTED', RECL=2*IRECLN)         
C                                                                       
      END IF                                                            
C                                                                       
C     ----- IS THIS A NEW OR OLD DAF FILE? -----                        
C     EITHER MARK THE NEW DAF RECORDS AS EMPTY                          
C     OR ELSE LOAD THE OLD DAF DIRECTORY                                
C                                                                       
      IF(IREST.EQ.0) THEN                                               
         IRECST = 1                                                     
         DO 100 I = 1,950                                               
            IODA(I) = -1                                                
  100    CONTINUE                                                       
         IRECST = IRECST + 1                                            
         IF(MASWRK) WRITE(UNIT=IDAF, REC=1) IRECST,IODA,IFILEN,IS,IPK   
      ELSE                                                              
         IF(MASWRK) READ (UNIT=IDAF, REC=1) IRECST,IODA,IFILEN,IS,IPK   
         IF(GOPARR) THEN                                                
            CALL DDI_BCAST(200,'I',IRECST,1,MASTER)                     
            CALL DDI_BCAST(201,'I',IODA,950,MASTER)                     
            CALL DDI_BCAST(202,'I',IFILEN,950,MASTER)                   
            CALL DDI_BCAST(203,'I',IS,1,MASTER)                         
            CALL DDI_BCAST(204,'I',IPK,1,MASTER)                        
         END IF                                                         
      END IF                                                            
      RETURN                                                            
      END                                                               
