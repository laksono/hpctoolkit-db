c 17 Oct 19 - DGF - patch LAPACK KDIAG bug                              
C  6 Jun 18 - DGF - tweaks for FMO 5.3                                  
C 22 Oct 14 - DGF - PAD MACHSW COMMON                                   
C  7 Sep 12 - MWS - RENAME EPSLON TO DEPSLON                            
C  7 MAR 12 - MWS - PAD MACHSW COMMON                                   
C 28 DEC 11 - DGF - PAD COMMONS FOR FMO 4.2                             
C 11 AUG 11 - AA  - GLDIAG: ADD, AS COMMENTS, AN LAPACK CALL            
C 10 MAY 10 - AA  - PARALLEL EVVRSP (ETRED3,ELAU,FREDA,EINVIT,ETRBK3)   
C 22 DEC 06 - MWS - USE IMPLICIT STATEMENTS TO COVER BLAS NAME HACK     
C  5 JUL 05 - MWS - SELECT NEW ATOM,BASIS,EFP,PCM,DAF DIMENSIONS        
C  1 JUN 05 - DGF - GLDIAG: OPTION FORCING ALL DIAGS TO MASTER (FOR FMO)
C  3 FEB 97 - MWS - MAKE EINPI A ROUTINE, AS IT MODIFIES ARG ANORM.     
C 16 JAN 97 - STE - EINVIT: AVOID ACCIDENTAL ROW EXCHANGES              
C 10 AUG 94 - MWS - INCREASE NUMBER OF DAF RECORDS                      
C 31 MAR 94 - MWS - ADD A VARIABLE TO END OF MACHSW COMMON              
C 26 JUN 93 - MWS - ETRED3: ADD RETURN FOR SPECIAL CASE N=1             
C  4 JAN 92 - TLW - MAKE WRITES PARALLEL;ADD COMMON PAR                 
C 30 AUG 91 - MWS - JACDIA: LIMIT ITERATIONS, USE EPSLON IN TEST.       
C 14 JUL 91 - MWS - JACOBI DIAGONALIZATION ALLOWS FOR LDVEC.NE.N        
C 29 JAN 91 - TLW - GLDIAG: CHANGED COMMON DIAGSW TO MACHSW             
C 29 OCT 90 - STE - FIX JACDIA UNDEFINED VARIABLE BUG                   
C 14 SEP 90 - MK  - NEW JACOBI DIAGONALIZATION (KDIAG=3)                
C 27 MAR 88 - MWS - ALLOW FOR VECTOR ROUTINE IN GLDIAG                  
C 11 AUG 87 - MWS - SANITIZE CONSTANTS IN EQLRAT                        
C 15 FEB 87 - STE - FIX EINVIT SUB-MATRIX LOOP LIMIT                    
C                   SCRATCH ARRAYS ARE N*8 REAL AND N INTEGER           
C  8 DEC 86 - STE - USE PERF INDEX FROM EINPI1 TO JUDGE EINVIT FAILURE  
C 30 NOV 86 - STE - DELETE LIGENB, MAKE EVVRSP DEFAULT                  
C                   (GIVEIS FAILS ON CRAY FOR BENCHMC AND BENCHCI)      
C  7 JUL 86 - JAB - SANITIZE FLOATING POINT CONSTANTS                   
C 11 OCT 85 - STE - LIGENB,TQL2: USE DROT,DSWAP; TINVTB: SCALE VECTOR   
C                   BEFORE NORMALIZING; GENERIC FUNCTIONS               
C 24 FEB 84 - STE - INITIALIZE INDEX ARRAY FOR LIGENB IN GLDIAG         
C  1 DEC 83 - STE - CHANGE MACHEP FROM 2**-54 TO 2**-50                 
C 28 SEP 82 - MWS - CONVERT TO IBM                                      
C                                                                       
C*MODULE EIGEN   *DECK GLDIAG                                           
      SUBROUTINE GLDIAG(LDVECT,NVECT,N,H,WRK,EIG,VECTOR,IERR,IWRK)      
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      LOGICAL GOPARR,DSKWRK,MASWRK                                      
C                                                                       
      DIMENSION H(*),WRK(N,8),EIG(N),VECTOR(LDVECT,NVECT),IWRK(N)       
c     ALLOCATABLE :: IWRK2(:)                                           
C                                                                       
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)              
      COMMON /MACHSW/ KDIAGG,ICORFL,IXDR,MODIO,mem10,lpnt10,mem10m      
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
c     COMMON /MTXSIZ/ MXSEQ2,MXSEQ3                                     
C                                                                       
C     ----- GENERAL ROUTINE TO DIAGONALIZE A SYMMETRIC MATRIX -----     
C     IF KDIAG = 0, USE A ROUTINE FROM THE VECTOR LIBRARY,              
C                   IF AVAILABLE (SEE THE ROUTINE 'GLDIAG'              
C                   IN VECTOR.SRC), OR EVVRSP OTHERWISE                 
C              = 1, USE EVVRSP                                          
C              = 2, USE GIVEIS                                          
C              = 3, USE JACOBI                                          
C     100 ADDED TO KDIAG FORCES DIAGONALISING ONLY ON MASTER AND        
C     BROADCASTING THE RESULTS ON SLAVES (USEFUL FOR MIXED NODE CLUSTERS
C     AND/OR FOR CLUSTERS WITH CHARGES ACRUED BY CPU TIME).             
C                                                                       
C           N      = DIMENSION (ORDER) OF MATRIX TO BE SOLVED           
C           LDVECT = LEADING DIMENSION OF VECTOR                        
C           NVECT  = NUMBER OF VECTORS DESIRED                          
C           H      = MATRIX TO BE DIAGONALIZED                          
C           WRK    = N*8 W.P. REAL WORDS OF SCRATCH SPACE               
C           EIG    = EIGENVALUES  (OUTPUT)                              
C           VECTOR = EIGENVECTORS (OUTPUT)                              
C           IERR   = ERROR FLAG (OUTPUT)                                
C           IWRK   = N INTEGER WORDS OF SCRATCH SPACE                   
C                                                                       
C        TRADITIONAL RUNS WITH -KDIAG- WILL HAVE LDIAG=0 AND KDIAG SAME.
C        ADDING 100 TO KDIAG BEFORE CALLING THIS WILL RESULT IN THE     
C        DIAGONALIZATION OCCURING ENTIRELY IN THE MASTER PROCESS.       
C                                                                       
      KDIAG=MOD(KDIAGG,100)                                             
      IF(KDIAG.GE.4.AND.(NVECT.NE.N.OR.LDVECT.LT.N)) KDIAG=0            
      KDIAGS=KDIAG                                                      
      IF(KDIAG.GE.4) THEN                                               
        KDIAG=0                                                         
C*LAP   KDIAG=KDIAGS                                                    
c       Use KDIAG=0 if asked for LAPACK but it is not activated.        
      ENDIF                                                             
      LDIAG=KDIAGG/100                                                  
      IF (KDIAG.LE.1) LDIAG = 0 !EVVRSP IS PARALLEL.                    
C                                                                       
      IF(LDIAG.EQ.0.OR.MASWRK) THEN                                     
         IERR = 0                                                       
C                                                                       
C         ----- USE STEVE ELBERT'S ROUTINE -----                        
C                                                                       
         IF(KDIAG.LE.1) THEN                                            
            LENH = (N*N+N)/2                                            
            KORDER =0                                                   
c           write(6,*) 'wwwiii',MXSEQ2,MXSEQ3                           
            if(iand(modio,128).ne.0) call matpar(3,0,max2sav,max3sav)   
c           write(6,*) 'www000',MXSEQ2,MXSEQ3                           
            CALL EVVRSP(IW,N,NVECT,LENH,LDVECT,H,WRK,IWRK,EIG,          
     *                  VECTOR,KORDER,IERR)                             
            if(iand(modio,128).ne.0) call matpar(3,1,max2sav,max3sav)   
c           write(6,*) 'wwwfff',MXSEQ2,MXSEQ3                           
         END IF                                                         
C                                                                       
C         ----- USE MODIFIED EISPAK ROUTINE -----                       
C                                                                       
         IF(KDIAG.EQ.2)                                                 
     *      CALL GIVEIS(N,NVECT,LDVECT,H,WRK,IWRK,EIG,VECTOR,IERR)      
C                                                                       
C         ----- USE JACOBI ROTATION ROUTINE -----                       
C                                                                       
         IF(KDIAG.EQ.3) THEN                                            
            IF(NVECT.EQ.N) THEN                                         
               CALL JACDG(H,VECTOR,EIG,IWRK,WRK,LDVECT,N)               
            ELSE                                                        
               IF (MASWRK) WRITE(IW,9000) N,NVECT,LDVECT                
               CALL ABRT                                                
            END IF                                                      
         END IF                                                         
C                                                                       
C           THE INACTIVE CODE BELOW WOULD LET A THREADED LAPACK         
C           LIBRARY ROUTINE HANDLE THE DIAGONALIZATION.                 
C           THE "KDIAG" CONDITIONS HAVE BEEN SET TO REPLACE EVVRSP,     
C           SO YOU'D PROBABLY WANT TO COMMENT THAT OUT, ABOVE.          
C           EVVRSP IS PREFERRED IN A PARALLEL PROCESS-BASED PROGRAM,    
C           AS IT CONTAINS MESSAGE-PASSING PARALLELIZATION.             
C                                                                       
*L64     IF (KDIAG.EQ.4) THEN                                           
*L64         LDIAG = 1                                                  
*L64         IF (NVECT.EQ.N .AND. LDVECT.GE.N) THEN                     
*L64             IF (MASWRK) THEN                                       
*L64                 CALL DSPEV('V','U',N,H,EIG,VECTOR,LDVECT,WRK,INFO) 
*L64                 IF (INFO.NE.0) THEN                                
*L64                     WRITE(IW,*) 'DSPEV FAILED!'                    
*L64                     CALL ABRT()                                    
*L64                 END IF                                             
*L64             END IF                                                 
c             ELSE IF (NVECT.LT.N) THEN                                 
c                IF (MASWRK) THEN                                       
c                    ALLOCATE(IWRK2(5*N))                               
c                    CALL DSPEVX('V','I','U',N,H,DUM,DUM,               
c    *                   1,NVECT,1.0D-12,NEIGFOUND,EIG,VECTOR,LDVECT,   
c    *                   WRK,IWRK2,IWRK,INFO)                           
c                    IF (INFO.NE.0.OR.NEIGFOUND.LT.NVECT) THEN          
c                        WRITE(IW,*) 'DSPEV FAILED!'                    
c                        CALL ABRT()                                    
c                    END IF                                             
c                    DEALLOCATE(IWRK2)                                  
c                END IF                                                 
*L64            ELSE                                                    
*L64                IF (MASWRK) WRITE(IW,*) 'ERROR, TRY DIFFERENT KDIAG'
*L64                CALL ABRT                                           
*L64            END IF                                                  
*L64        END IF                                                      
C---     IF(KDIAG.LE.1 .OR. KDIAG.GT.3) THEN                            
C---        IF(NVECT.EQ.N .AND. LDVECT.GE.N) THEN                       
C---           LDIAG = 1                                                
C---           IF (MASWRK) THEN                                         
C---              IJ = 1                                                
C---              DO J = 1,N                                            
C---                 DO I = 1,J                                         
C---                    VECTOR(I,J) = H(IJ)                             
C---                    IJ = IJ+1                                       
C---                 ENDDO                                              
C---              ENDDO                                                 
C---              INFO = 0                                              
C---              IF (N.LE.4) THEN                                      
C---                 CALL DSYEV('V', 'U', N, VECTOR, LDVECT,            
C--- *                    EIG, WRK, N*8, INFO)                          
C---              ELSE                                                  
C---                 CALL DSYEV('V', 'U', N, VECTOR, LDVECT,            
C--- *                    EIG, H, (N*N+N)/2, INFO)                      
C---              ENDIF                                                 
C---              IF (INFO.NE.0) CALL ABRT()                            
C---           ENDIF                                                    
C---        ENDIF                                                       
C---     ENDIF                                                          
C                                                                       
         IF(KDIAG.EQ.5) THEN                                            
c           Because DSYEV cannot treat this, use KDIAG=0 instead...     
c           IF(NVECT.EQ.N .AND. LDVECT.GE.N) THEN                       
               IJ = 1                                                   
               DO J = 1,N                                               
                  DO I = 1,J                                            
                     VECTOR(I,J) = H(IJ)                                
                     IJ = IJ+1                                          
                  ENDDO                                                 
               ENDDO                                                    
C*LAP          CALL DSYEV('V','U',N,VECTOR,LDVECT,EIG,WRK,N*8,INFO)     
               IF (INFO.NE.0) CALL ABRT                                 
c           ELSE                                                        
c              IF (MASWRK) WRITE(IW,*) 'SWITCH to KDIAG<4'              
c              CALL ABRT                                                
c           ENDIF                                                       
         ENDIF                                                          
      ENDIF                                                             
C                                                                       
C        BROADCAST THE RESULTS                                          
C                                                                       
      IF(LDIAG.NE.0.AND.GOPARR) THEN                                    
         CALL DDI_BCAST(21028,'I',IERR,1,MASTER)                        
         CALL DDI_BCAST(21029,'F',EIG,N,MASTER)                         
         CALL DDI_BCAST(21030,'F',VECTOR,LDVECT*NVECT,MASTER)           
      ENDIF                                                             
C                                                                       
      RETURN                                                            
C                                                                       
 9000 FORMAT(1X,'IN -GLDIAG-, N,NVECT,LDVECT=',3I8/                     
     *       1X,'THE JACOBI CODE CANNOT COPE WITH N.NE.NVECT!'/         
     *       1X,'SO THIS RUN DOES NOT PERMIT KDIAG=3.')                 
      END                                                               
C                                                                       
C         START OF -EVVRSP- PACKAGE, THE ONLY PARALLEL DIAGONALIZER HERE
C                                                                       
C*MODULE EIGEN   *DECK EVVRSP                                           
      SUBROUTINE EVVRSP(MSGFL,N,NVECT,LENA,NV,A,B,IND,ROOT,             
     *                  VECT,IORDER,IERR)                               
C*                                                                      
C*    AUTHOR:  S. T. ELBERT, AMES LABORATORY-USDOE, JUNE 1985           
C*                                                                      
C*    PURPOSE -                                                         
C*       FINDS   (ALL) EIGENVALUES    AND    (SOME OR ALL) EIGENVECTORS 
C*                     *    *                                   *       
C*       OF A REAL SYMMETRIC PACKED MATRIX.                             
C*            *    *         *                                          
C*                                                                      
C*    METHOD -                                                          
C*       THE METHOD AS PRESENTED IN THIS ROUTINE CONSISTS OF FOUR STEPS:
C*       FIRST, THE INPUT MATRIX IS REDUCED TO TRIDIAGONAL FORM BY THE  
C*       HOUSEHOLDER TECHNIQUE (ORTHOGONAL SIMILARITY TRANSFORMATIONS). 
C*       SECOND, THE ROOTS ARE LOCATED USING THE RATIONAL QL METHOD.    
C*       THIRD, THE VECTORS OF THE TRIDIAGONAL FORM ARE EVALUATED BY THE
C*       INVERSE ITERATION TECHNIQUE.  VECTORS FOR DEGENERATE OR NEAR-  
C*       DEGENERATE ROOTS ARE FORCED TO BE ORTHOGONAL.                  
C*       FOURTH, THE TRIDIAGONAL VECTORS ARE ROTATED TO VECTORS OF THE  
C*       ORIGINAL ARRAY.                                                
C*                                                                      
C*       THESE ROUTINES ARE MODIFICATIONS OF THE EISPACK 3              
C*       ROUTINES TRED3, TQLRAT, TINVIT AND TRBAK3                      
C*                                                                      
C*       FOR FURTHER DETAILS, SEE EISPACK USERS GUIDE, B. T. SMITH      
C*       ET AL, SPRINGER-VERLAG, LECTURE NOTES IN COMPUTER SCIENCE,     
C*       VOL. 6, 2-ND EDITION, 1976.  ANOTHER GOOD REFERENCE IS         
C*       THE SYMMETRIC EIGENVALUE PROBLEM BY B. N. PARLETT              
C*       PUBLISHED BY PRENTICE-HALL, INC., ENGLEWOOD CLIFFS, N.J. (1980)
C*                                                                      
C*    ON ENTRY -                                                        
C*       MSGFL  - INTEGER (LOGICAL UNIT NO.)                            
C*                FILE WHERE ERROR MESSAGES WILL BE PRINTED.            
C*                IF MSGFL IS 0, ERROR MESSAGES WILL BE PRINTED ON LU 6.
C*                IF MSGFL IS NEGATIVE, NO ERROR MESSAGES PRINTED.      
C*       N      - INTEGER                                               
C*                ORDER OF MATRIX A.                                    
C*       NVECT  - INTEGER                                               
C*                NUMBER OF VECTORS DESIRED.  0 .LE. NVECT .LE. N.      
C*       LENA   - INTEGER                                               
C*                DIMENSION OF  A  IN CALLING ROUTINE.  MUST NOT BE LESS
C*                THAN (N*N+N)/2.                                       
C*       NV     - INTEGER                                               
C*                ROW DIMENSION OF VECT IN CALLING ROUTINE.   N .LE. NV.
C*       A      - WORKING PRECISION REAL (LENA)                         
C*                INPUT MATRIX, ROWS OF THE LOWER TRIANGLE PACKED INTO  
C*                LINEAR ARRAY OF DIMENSION N*(N+1)/2.  THE PACKED ORDER
C*                IS A(1,1), A(2,1), A(2,2), A(3,1), A(3,2), ...        
C*       B      - WORKING PRECISION REAL (N,8)                          
C*                SCRATCH ARRAY, 8*N ELEMENTS                           
C*       IND    - INTEGER (N)                                           
C*                SCRATCH ARRAY OF LENGTH N.                            
C*       IORDER - INTEGER                                               
C*                ROOT ORDERING FLAG.                                   
C*                = 0, ROOTS WILL BE PUT IN ASCENDING ORDER.            
C*                = 2, ROOTS WILL BE PUT IN DESCENDING ORDER.           
C*                                                                      
C*    ON EXIT -                                                         
C*       A      - DESTORYED.  NOW HOLDS REFLECTION OPERATORS.           
C*       ROOT   - WORKING PRECISION REAL (N)                            
C*                ALL EIGENVALUES IN ASCENDING OR DESCENDING ORDER.     
C*                  IF IORDER = 0, ROOT(1) .LE. ... .LE. ROOT(N)        
C*                  IF IORDER = 2, ROOT(1) .GE. ... .GE. ROOT(N)        
C*       VECT   - WORKING PRECISION REAL (NV,NVECT)                     
C*                EIGENVECTORS FOR ROOT(1), ..., ROOT(NVECT).           
C*       IERR   - INTEGER                                               
C*                = 0 IF NO ERROR DETECTED,                             
C*                = K IF ITERATION FOR K-TH EIGENVALUE FAILED,          
C*                = -K IF ITERATION FOR K-TH EIGENVECTOR FAILED.        
C*                (FAILURES SHOULD BE VERY RARE.  CONTACT C. MOLER.)    
C*                                                                      
C                                                                       
C                                                                       
      DOUBLE PRECISION A(LENA)                                          
      DOUBLE PRECISION B(N,8)                                           
      DOUBLE PRECISION ROOT(N)                                          
      DOUBLE PRECISION T                                                
      DOUBLE PRECISION VECT(NV,*)                                       
C                                                                       
      INTEGER IND(N)                                                    
C                                                                       
      DOUBLE PRECISION XX                                               
      LOGICAL GOPARR,DSKWRK,MASWRK                                      
      COMMON /FMCOM / XX(1)                                             
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
C                                                                       
  900 FORMAT(26H0*** EVVRSP PARAMETERS ***/                             
     +       14H ***      N = ,I8,4H ***/                               
     +       14H ***  NVECT = ,I8,4H ***/                               
     +       14H ***   LENA = ,I8,4H ***/                               
     +       14H ***     NV = ,I8,4H ***/                               
     +       14H *** IORDER = ,I8,4H ***/                               
     +       14H ***   IERR = ,I8,4H ***)                               
  901 FORMAT(37H VALUE OF LENA IS LESS THAN (N*N+N)/2)                  
  902 FORMAT(39H EQLRAT HAS FAILED TO CONVERGE FOR ROOT,I5)             
  903 FORMAT(18H NV IS LESS THAN N)                                     
  904 FORMAT(41H EINVIT HAS FAILED TO CONVERGE FOR VECTOR,I5)           
  905 FORMAT(51H VALUE OF IORDER MUST BE 0 (SMALLEST ROOT FIRST) OR     
     *      ,23H 2 (LARGEST ROOT FIRST))                                
  906 FORMAT(' VALUE OF N IS LESS THAN OR EQUAL ZERO')                  
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
      LMSGFL=MSGFL                                                      
      IF (MSGFL .EQ. 0) LMSGFL=6                                        
      IERR = N - 1                                                      
      IF (N .LE. 0) GO TO 800                                           
      IERR = N + 1                                                      
      IF ( (N*N+N)/2 .GT. LENA) GO TO 810                               
C                                                                       
C        REDUCE REAL SYMMETRIC MATRIX A TO TRIDIAGONAL FORM             
C                                                                       
      CALL ETRED3(N,LENA,A,B(1,1),B(1,2),B(1,3))                        
C                                                                       
C        FIND ALL EIGENVALUES OF TRIDIAGONAL MATRIX                     
C                                                                       
      CALL EQLRAT(N,B(1,1),B(1,2),B(1,3),ROOT,IND,IERR,B(1,4))          
      IF (IERR .NE. 0) GO TO 820                                        
C                                                                       
C         CHECK THE DESIRED ORDER OF THE EIGENVALUES                    
C                                                                       
      B(1,3) = IORDER                                                   
      IF (IORDER .EQ. 0) GO TO 300                                      
         IF (IORDER .NE. 2) GO TO 850                                   
C                                                                       
C         ORDER ROOTS IN DESCENDING ORDER (LARGEST FIRST)...            
C        TURN ROOT AND IND ARRAYS END FOR END                           
C                                                                       
         DO 210 I = 1, N/2                                              
            J = N+1-I                                                   
            T = ROOT(I)                                                 
            ROOT(I) = ROOT(J)                                           
            ROOT(J) = T                                                 
            L = IND(I)                                                  
            IND(I) = IND(J)                                             
            IND(J) = L                                                  
  210    CONTINUE                                                       
C                                                                       
C           FIND I AND J MARKING THE START AND END OF A SEQUENCE        
C           OF DEGENERATE ROOTS                                         
C                                                                       
         I=0                                                            
  220    CONTINUE                                                       
            I = I+1                                                     
            IF (I .GT. N) GO TO 300                                     
            DO 230 J=I,N                                                
               IF (ROOT(J) .NE. ROOT(I)) GO TO 240                      
  230       CONTINUE                                                    
            J = N+1                                                     
  240       CONTINUE                                                    
            J = J-1                                                     
            IF (J .EQ. I) GO TO 220                                     
C                                                                       
C                    TURN AROUND IND BETWEEN I AND J                    
C                                                                       
            JSV = J                                                     
            KLIM = (J-I+1)/2                                            
            DO 250 K=1,KLIM                                             
               L = IND(J)                                               
               IND(J) = IND(I)                                          
               IND(I) = L                                               
               I = I+1                                                  
               J = J-1                                                  
  250       CONTINUE                                                    
            I = JSV                                                     
         GO TO 220                                                      
C                                                                       
  300 CONTINUE                                                          
C                                                                       
      IF (NVECT .LE. 0) RETURN                                          
      IF (NV .LT. N) GO TO 830                                          
C                                                                       
C        FIND EIGENVECTORS OF TRI-DIAGONAL MATRIX VIA INVERSE ITERATION 
C        PARALLEL REQUIRES EXTRA MEMORY, BEYOND EARLIER VERSION NEEDS.  
C        IT WOULD BE BETTER TO REQUIRE 12 WORK VECTORS FOR GLDIAG!      
C        THIS ISN'T LOCATED INSIDE -EVVRSP- AS THERE ARE EVEN DIRECT    
C        CALLS TO THIS SPECIFIC PACKAGE INSIDE GAMESS.                  
C                                                                       
      CALL VALFM(LOADFM)                                                
      LRV7 = LOADFM + 1                                                 
      LRV8 = LRV7   + 2*N                                               
      LAST = LRV8   + N                                                 
      NEED = LAST - LOADFM - 1                                          
      CALL GETFM(NEED)                                                  
C                                                                       
      IERR = LMSGFL                                                     
      CALL EINVIT(NV,N,B(1,1),B(1,2),B(1,3),NVECT,ROOT,IND,VECT,        
     *            IERR,B(1,4),B(1,5),B(1,6),B(1,7),B(1,8),              
     *            XX(LRV7),XX(LRV8))                                    
      CALL RETFM(NEED)                                                  
      IF (IERR .NE. 0) GO TO 840                                        
C                                                                       
C        FIND EIGENVECTORS OF SYMMETRIC MATRIX VIA BACK TRANSFORMATION  
C                                                                       
  400 CONTINUE                                                          
      CALL ETRBK3(NV,N,LENA,A,NVECT,VECT)                               
      RETURN                                                            
C                                                                       
C        ERROR MESSAGE SECTION                                          
C                                                                       
  800 IF (LMSGFL .LT. 0) RETURN                                         
      IF (MASWRK) WRITE(LMSGFL,906)                                     
      GO TO 890                                                         
C                                                                       
  810 IF (LMSGFL .LT. 0) RETURN                                         
      IF (MASWRK) WRITE(LMSGFL,901)                                     
      GO TO 890                                                         
C                                                                       
  820 IF (LMSGFL .LT. 0) RETURN                                         
      IF (MASWRK) WRITE(LMSGFL,902) IERR                                
      GO TO 890                                                         
C                                                                       
  830 IF (LMSGFL .LT. 0) RETURN                                         
      IF (MASWRK) WRITE(LMSGFL,903)                                     
      GO TO 890                                                         
C                                                                       
  840 CONTINUE                                                          
      IF ((LMSGFL .GT. 0).AND.MASWRK) WRITE(LMSGFL,904) -IERR           
      GO TO 400                                                         
C                                                                       
  850 IERR=-1                                                           
      IF (LMSGFL .LT. 0) RETURN                                         
      IF (MASWRK) WRITE(LMSGFL,905)                                     
      GO TO 890                                                         
C                                                                       
  890 CONTINUE                                                          
      IF (MASWRK) WRITE(LMSGFL,900) N,NVECT,LENA,NV,IORDER,IERR         
      RETURN                                                            
      END                                                               
C                                                                       
C*MODULE EIGEN   *DECK DEPSLON                                          
      DOUBLE PRECISION FUNCTION DEPSLON(X)                              
C*                                                                      
C*    AUTHORS -                                                         
C*       THIS ROUTINE WAS TAKEN FROM EISPACK EDITION 3 DATED 4/6/83     
C*       THIS VERSION IS BY S. T. ELBERT, AMES LABORATORY-USDOE NOV 1986
C*                                                                      
C*    PURPOSE -                                                         
C*       ESTIMATE UNIT ROUNDOFF IN QUANTITIES OF SIZE X.                
C*                                                                      
C*    ON ENTRY -                                                        
C*       X      - WORKING PRECISION REAL                                
C*                VALUES TO FIND DEPSLON FOR                            
C*                                                                      
C*    ON EXIT -                                                         
C*       DEPSLON - WORKING PRECISION REAL                               
C*                SMALLEST POSITIVE VALUE SUCH THAT X+DEPSLON.NE.ZERO   
C*                                                                      
C*    QUALIFICATIONS -                                                  
C*       THIS ROUTINE SHOULD PERFORM PROPERLY ON ALL SYSTEMS            
C*       SATISFYING THE FOLLOWING TWO ASSUMPTIONS,                      
C*          1.  THE BASE USED IN REPRESENTING FLOATING POINT            
C*              NUMBERS IS NOT A POWER OF THREE.                        
C*          2.  THE QUANTITY  A  IN STATEMENT 10 IS REPRESENTED TO      
C*              THE ACCURACY USED IN FLOATING POINT VARIABLES           
C*              THAT ARE STORED IN MEMORY.                              
C*       THE STATEMENT NUMBER 10 AND THE GO TO 10 ARE INTENDED TO       
C*       FORCE OPTIMIZING COMPILERS TO GENERATE CODE SATISFYING         
C*       ASSUMPTION 2.                                                  
C*       UNDER THESE ASSUMPTIONS, IT SHOULD BE TRUE THAT,               
C*              A  IS NOT EXACTLY EQUAL TO FOUR-THIRDS,                 
C*              B  HAS A ZERO FOR ITS LAST BIT OR DIGIT,                
C*              C  IS NOT EXACTLY EQUAL TO ONE,                         
C*              EPS  MEASURES THE SEPARATION OF 1.0 FROM                
C*                   THE NEXT LARGER FLOATING POINT NUMBER.             
C*       THE DEVELOPERS OF EISPACK WOULD APPRECIATE BEING INFORMED      
C*       ABOUT ANY SYSTEMS WHERE THESE ASSUMPTIONS DO NOT HOLD.         
C*                                                                      
C*    DIFFERENCES FROM EISPACK 3 -                                      
C*       USE IS MADE OF PARAMETER STATEMENTS AND INTRINSIC FUNCTIONS    
C*       --NO EXECUTEABLE CODE CHANGES--                                
C*                                                                      
C*    NOTE -                                                            
C*       QUESTIONS AND COMMENTS CONCERNING EISPACK SHOULD BE DIRECTED TO
C*       B. S. GARBOW, APPLIED MATH. DIVISION, ARGONNE NATIONAL LAB.    
C                                                                       
      DOUBLE PRECISION A,B,C,EPS,X                                      
      DOUBLE PRECISION ZERO, ONE, THREE, FOUR                           
C                                                                       
      PARAMETER (ZERO=0.0D+00)                                          
      PARAMETER (ONE=1.0D+00)                                           
      PARAMETER (THREE=3.0D+00)                                         
      PARAMETER (FOUR=4.0D+00)                                          
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
      A = FOUR/THREE                                                    
   10 B = A - ONE                                                       
      C = B + B + B                                                     
      EPS = ABS(C - ONE)                                                
      IF (EPS .EQ. ZERO) GO TO 10                                       
      DEPSLON = EPS*ABS(X)                                              
      RETURN                                                            
      END                                                               
C                                                                       
C*MODULE EIGEN   *DECK ETRED3                                           
      SUBROUTINE ETRED3(N,NV,A,D,E,E2)                                  
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C*                                                                      
C*    AUTHORS -                                                         
C*       THIS IS A MODIFICATION OF ROUTINE TRED3 FROM EISPACK EDITION 3 
C*       DATED AUGUST 1983.                                             
C*       EISPACK TRED3 IS A TRANSLATION OF THE ALGOL PROCEDURE TRED3,   
C*       NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
C*       HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C*       THIS VERSION IS BY S. T. ELBERT, AMES LABORATORY-USDOE JUN 1986
C*                                                                      
C*    PURPOSE -                                                         
C*       THIS ROUTINE REDUCES A REAL SYMMETRIC (PACKED) MATRIX, STORED  
C*       AS A ONE-DIMENSIONAL ARRAY, TO A SYMMETRIC TRIDIAGONAL MATRIX  
C*       USING ORTHOGONAL SIMILARITY TRANSFORMATIONS, PRESERVING THE    
C*       INFORMATION ABOUT THE TRANSFORMATIONS IN  A.                   
C*                                                                      
C*    METHOD -                                                          
C*       THE TRIDIAGONAL REDUCTION IS PERFORMED IN THE FOLLOWING WAY.   
C*       STARTING WITH J=N, THE ELEMENTS IN THE J-TH ROW TO THE         
C*       LEFT OF THE DIAGONAL ARE FIRST SCALED, TO AVOID POSSIBLE       
C*       UNDERFLOW IN THE TRANSFORMATION THAT MIGHT RESULT IN SEVERE    
C*       DEPARTURE FROM ORTHOGONALITY.  THE SUM OF SQUARES  SIGMA  OF   
C*       THESE SCALED ELEMENTS IS NEXT FORMED.  THEN, A VECTOR  U  AND  
C*       A SCALAR                                                       
C*                      H = U(TRANSPOSE) * U / 2                        
C*       DEFINE A REFLECTION OPERATOR                                   
C*                      P = I - U * U(TRANSPOSE) / H                    
C*       WHICH IS ORTHOGONAL AND SYMMETRIC AND FOR WHICH THE            
C*       SIMILIARITY TRANSFORMATION  PAP  ELIMINATES THE ELEMENTS IN    
C*       THE J-TH ROW OF  A  TO THE LEFT OF THE SUBDIAGONAL AND THE     
C*       SYMMETRICAL ELEMENTS IN THE J-TH COLUMN.                       
C*                                                                      
C*       THE NON-ZERO COMPONENTS OF  U  ARE THE ELEMENTS OF THE J-TH    
C*       ROW TO THE LEFT OF THE DIAGONAL WITH THE LAST OF THEM          
C*       AUGMENTED BY THE SQUARE ROOT OF  SIGMA  PREFIXED BY THE SIGN   
C*       OF THE SUBDIAGONAL ELEMENT.  BY STORING THE TRANSFORMED SUB-   
C*       DIAGONAL ELEMENT IN  E(J)  AND NOT OVERWRITING THE ROW         
C*       ELEMENTS ELIMINATED IN THE TRANSFORMATION, FULL INFORMATION    
C*       ABOUT  P  IS SAVE FOR LATER USE IN  ETRBK3.                    
C*                                                                      
C*       THE TRANSFORMATION SETS  E2(J)  EQUAL TO  SIGMA  AND  E(J)     
C*       EQUAL TO THE SQUARE ROOT OF  SIGMA  PREFIXED BY THE SIGN       
C*       OF THE REPLACED SUBDIAGONAL ELEMENT.                           
C*                                                                      
C*       THE ABOVE STEPS ARE REPEATED ON FURTHER ROWS OF THE            
C*       TRANSFORMED  A  IN REVERSE ORDER UNTIL  A  IS REDUCED TO TRI-  
C*       DIAGONAL FORM, THAT IS, REPEATED FOR  J = N-1,N-2,...,3.       
C*                                                                      
C*    COMPLEXITY -                                                      
C*       2/3 N**3                                                       
C*                                                                      
C*    ON ENTRY-                                                         
C*       N      - INTEGER                                               
C*                THE ORDER OF THE MATRIX.                              
C*       NV     - INTEGER                                               
C*                MUST BE SET TO THE DIMENSION OF THE ARRAY PARAMETER A 
C*                AS DECLARED IN THE CALLING ROUTINE DIMENSION STATEMENT
C*       A      - W.P. REAL (NV)                                        
C*                CONTAINS THE LOWER TRIANGLE OF THE REAL SYMMETRIC     
C*                INPUT MATRIX, STORED ROW-WISE AS A ONE-DIMENSIONAL    
C*                ARRAY, IN ITS FIRST N*(N+1)/2 POSITIONS.              
C*                                                                      
C*    ON EXIT-                                                          
C*       A      - W.P. REAL (NV)                                        
C*                CONTAINS INFORMATION ABOUT THE ORTHOGONAL             
C*                TRANSFORMATIONS USED IN THE REDUCTION.                
C*       D      - W.P. REAL (N)                                         
C*                CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL     
C*                MATRIX.                                               
C*       E      - W.P. REAL (N)                                         
C*                CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL  
C*                MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO
C*       E2     - W.P. REAL (N)                                         
C*                CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF 
C*                E. MAY COINCIDE WITH E IF THE SQUARES ARE NOT NEEDED. 
C*                                                                      
C*    DIFFERENCES FROM EISPACK 3 -                                      
C*       OUTER LOOP CHANGED FROM II=1,N TO I=N,3,-1                     
C*       PARAMETER STATEMENT AND GENERIC INTRINSIC FUNCTIONS USED       
C*       SCALE.NE.0 TEST NOW SPOTS TRI-DIAGONAL FORM                    
C*       VALUES LESS THAN DEPSLON CLEARED TO ZERO                       
C*       USE BLAS(1)                                                    
C*       U NOT COPIED TO D, LEFT IN A                                   
C*       E2 COMPUTED FROM E                                             
C*       INNER LOOPS SPLIT INTO ROUTINES ELAU AND A CALL TO DSPR2       
C*       INVERSE OF H STORED INSTEAD OF H                               
C*                                                                      
C*    NOTE -                                                            
C*       QUESTIONS AND COMMENTS CONCERNING EISPACK SHOULD BE DIRECTED TO
C*       B. S. GARBOW, APPLIED MATH. DIVISION, ARGONNE NATIONAL LAB.    
C                                                                       
      INTEGER I,IIA,IZ0,L,N,NV                                          
C                                                                       
      DOUBLE PRECISION A(NV),D(N),E(N),E2(N)                            
      DOUBLE PRECISION AIIMAX,F,G,H,HROOT,SCALE,SCALEI                  
      DOUBLE PRECISION ONE, ZERO                                        
C                                                                       
      PARAMETER (ZERO = 0.0D+00)                                        
      PARAMETER (ONE = 1.0D+00)                                         
C                                                                       
      LOGICAL GOPARR,DSKWRK,MASWRK,PARALL2                              
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
      IF (N .LE. 2) GO TO 310                                           
C                                                                       
C     THE COMMUNICATION COST OF PARALLEL HOUSEHOLDER IS BCAST + SUM,    
C     PER ITERATION.  ITERATIONS ABOVE SEQUENTIAL THESHOLD ARE DONE     
C     IN PARALLEL, ITERATIONS BELOW THE THRESHOLD ARE DONE IN SERIAL.   
C     THE THRESHOLD SHOULD BE AT LEAST 3 TO AVOID COMMUNICATION FOR     
C     FREQUENT 3X3 MATRICES.                                            
C                                                                       
      PARALL2 = GOPARR  .AND.  N-1.GT.MXSQN2()                          
      ICPU = ME                                                         
      NCPU = NPROC                                                      
C                                                                       
      IF (PARALL2) THEN                                                 
         IJ = 1                                                         
         DO I=1,N-1 ! N-1 AVOIDS REDUNDANT BROADCAST                    
            IF (ICPU.NE.MOD(I-1,NCPU)) CALL XCOPY(I,ZERO,0,A(IJ),1)     
            IJ = IJ + I                                                 
         ENDDO                                                          
      ENDIF                                                             
C                                                                       
      IZ0 = (N*N+N)/2                                                   
      AIIMAX = ABS(A(IZ0))                                              
      DO 300 I = N, 3, -1                                               
         L = I - 1                                                      
         IIA = IZ0                                                      
         IZ0 = IZ0 - I                                                  
C                                                                       
         IF (PARALL2) THEN                                              
            PARALL2 = L.GT.MXSQN2()                                     
            IF (PARALL2) THEN                                           
               IF (N.NE.I) THEN ! AVOID REDUNDANT BROADCAST             
                  CALL DDI_BCAST(21026,'F',A(IZ0+1),I,MOD(I-1,NCPU))    
               ENDIF                                                    
            ELSE                                                        
               CALL DDI_GSUMF(21027,A,(I*I+I)/2) ! SUM COMPUTED A       
            ENDIF                                                       
         ENDIF                                                          
C                                                                       
         AIIMAX = MAX(AIIMAX, ABS(A(IIA)))                              
         SCALE = XASUM(L, A(IZ0+1), 1)                                  
C                                                                       
         IF(SCALE .EQ. ABS(A(IIA-1)) .OR. AIIMAX+SCALE .EQ. AIIMAX) THEN
C                                                                       
C           THIS ROW IS ALREADY IN TRI-DIAGONAL FORM                    
C                                                                       
            D(I) = A(IIA)                                               
            IF (AIIMAX+D(I) .EQ. AIIMAX) D(I) = ZERO                    
            E(I) = A(IIA-1)                                             
            IF (AIIMAX+E(I) .EQ. AIIMAX) E(I) = ZERO                    
            E2(I) = E(I)*E(I)                                           
            A(IIA) = ZERO                                               
            GO TO 300                                                   
C                                                                       
         END IF                                                         
C                                                                       
         SCALEI = ONE / SCALE                                           
         CALL XSCAL(L,SCALEI,A(IZ0+1),1)                                
         HROOT = XNRM2(L,A(IZ0+1),1)                                    
C                                                                       
         F = A(IIA-1)                                                   
         G = -SIGN(HROOT,F)                                             
         E(I) = SCALE * G                                               
         E2(I) = E(I)*E(I)                                              
         H = HROOT*HROOT - F * G                                        
         A(IIA-1) = F - G                                               
         D(I) = A(IIA)                                                  
         A(IIA) = ONE / SQRT(H)                                         
C           .......... FORM P THEN Q IN E(1:L) ..........               
         CALL ELAU(ONE/H,L,A(IZ0+1),A,E,PARALL2,ICPU,NCPU)              
C           .......... FORM REDUCED A ..........                        
         CALL FREDA(L,A(IZ0+1),A,E,PARALL2,ICPU,NCPU)                   
C                                                                       
  300 CONTINUE                                                          
C                                                                       
  310 CONTINUE                                                          
      E(1) = ZERO                                                       
      E2(1)= ZERO                                                       
      D(1) = A(1)                                                       
      IF(N.EQ.1) RETURN                                                 
C                                                                       
      E(2) = A(2)                                                       
      E2(2)= A(2)*A(2)                                                  
      D(2) = A(3)                                                       
      RETURN                                                            
      END                                                               
C                                                                       
C*MODULE EIGEN   *DECK ELAU                                             
      SUBROUTINE ELAU(HINV,L,D,A,E,PARALL,IPROC,NPROC)                  
C                                                                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C                                                                       
      DOUBLE PRECISION HINV                                             
      DOUBLE PRECISION A(*)                                             
      DOUBLE PRECISION D(L)                                             
      DOUBLE PRECISION E(L)                                             
      LOGICAL PARALL                                                    
      INTEGER IPROC                                                     
      INTEGER NPROC                                                     
C                                                                       
      CALL XCOPY(L,0.0D+00,0,E,1)                                       
      IF (.NOT.PARALL.OR.IPROC.EQ.0) E(1) = A(1) * D(1)                 
C                                                                       
      JK = 2                                                            
      DO J=2,L                                                          
         IF (.NOT.PARALL.OR.IPROC.EQ.MOD(J-1,NPROC)) THEN               
            CALL XAXPY(J-1,D(J),A(JK),1,E,1)                            
            E(J) = XDOT(J-1,A(JK),1,D,1) + A(JK+J-1)*D(J)               
         ENDIF                                                          
         JK = JK + J                                                    
      ENDDO                                                             
C                                                                       
      IF (PARALL) CALL DDI_GSUMF(21024,E,L)                             
C                                                                       
C        .......... FORM P ..........                                   
C                                                                       
      CALL XSCAL(L,HINV,E,1)                                            
C                                                                       
C     .......... FORM Q ..........                                      
C                                                                       
      HH = 0.5D+00 * HINV * XDOT(L,E,1,D,1)                             
      CALL XAXPY(L,-HH,D(1),1,E(1),1)                                   
C                                                                       
      RETURN                                                            
      END                                                               
C                                                                       
C*MODULE EIGEN   *DECK FREDA                                            
      SUBROUTINE FREDA(L,D,A,E,PARALL,IPROC,NPROC)                      
C                                                                       
      INTEGER L                                                         
      DOUBLE PRECISION A(*)                                             
      DOUBLE PRECISION D(L)                                             
      DOUBLE PRECISION E(L)                                             
      LOGICAL PARALL                                                    
      INTEGER IPROC                                                     
      INTEGER NPROC                                                     
C                                                                       
      IF (PARALL) THEN ! UPDATE LOCAL SEGMENTS OF A                     
         JK = 1                                                         
         DO J=1,L                                                       
            IF (.NOT.PARALL.OR.IPROC.EQ.MOD(J-1,NPROC)) THEN            
               CALL XAXPY(J,-D(J),E,1,A(JK),1) ! A(J) = A(J) - D(J)E    
               CALL XAXPY(J,-E(J),D,1,A(JK),1) ! A(J) = A(J) - E(J)D    
            ENDIF                                                       
            JK = JK + J                                                 
         ENDDO                                                          
      ELSE ! USE BLAS2 OPERATION                                        
         CALL XSPR2('U',L,-1.0D+00,D,1,E,1,A)                           
      ENDIF                                                             
C                                                                       
      RETURN                                                            
      END                                                               
C                                                                       
C*MODULE EIGEN   *DECK EQLRAT                                           
      SUBROUTINE EQLRAT(N,DIAG,E,E2IN,D,IND,IERR,E2)                    
C*                                                                      
C*    AUTHORS -                                                         
C*       THIS IS A MODIFICATION OF ROUTINE EQLRAT FROM EISPACK EDITION 3
C*       DATED AUGUST 1983.                                             
C*       TQLRAT IS A TRANSLATION OF THE ALGOL PROCEDURE TQLRAT,         
C*       ALGORITHM 464, COMM. ACM 16, 689(1973) BY REINSCH.             
C*       THIS VERSION IS BY S. T. ELBERT (AMES LABORATORY-USDOE)        
C*                                                                      
C*    PURPOSE -                                                         
C*       THIS ROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC              
C*       TRIDIAGONAL MATRIX                                             
C*                                                                      
C*    METHOD -                                                          
C*       RATIONAL QL                                                    
C*                                                                      
C*    ON ENTRY -                                                        
C*       N      - INTEGER                                               
C*                THE ORDER OF THE MATRIX.                              
C*       D      - W.P. REAL (N)                                         
C*                CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.   
C*       E2     - W.P. REAL (N)                                         
C*                CONTAINS THE SQUARES OF THE SUBDIAGONAL ELEMENTS OF   
C*                THE INPUT MATRIX IN ITS LAST N-1 POSITIONS.           
C*                E2(1) IS ARBITRARY.                                   
C*                                                                      
C*     ON EXIT -                                                        
C*       D      - W.P. REAL (N)                                         
C*                CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN   
C*                ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND   
C*                ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE     
C*                THE SMALLEST EIGENVALUES.                             
C*       E2     - W.P. REAL (N)                                         
C*                DESTROYED.                                            
C*       IERR   - INTEGER                                               
C*                SET TO                                                
C*                ZERO       FOR NORMAL RETURN,                         
C*                J          IF THE J-TH EIGENVALUE HAS NOT BEEN        
C*                           DETERMINED AFTER 30 ITERATIONS.            
C*                                                                      
C*    DIFFERENCES FROM EISPACK 3 -                                      
C*       G=G+B INSTEAD OF IF(G.EQ.0) G=B ; B=B/4                        
C*       F77 BACKWARD LOOPS INSTEAD OF F66 CONSTRUCT                    
C*       GENERIC INTRINSIC FUNCTIONS                                    
C*       ARRARY  IND  ADDED FOR USE BY EINVIT                           
C*                                                                      
C*    NOTE -                                                            
C*       QUESTIONS AND COMMENTS CONCERNING EISPACK SHOULD BE DIRECTED TO
C*       B. S. GARBOW, APPLIED MATH. DIVISION, ARGONNE NATIONAL LAB.    
C                                                                       
      INTEGER I,J,L,M,N,II,L1,IERR                                      
      INTEGER IND(N)                                                    
C                                                                       
      DOUBLE PRECISION D(N),E(N),E2(N),DIAG(N),E2IN(N)                  
      DOUBLE PRECISION B,C,F,G,H,P,R,S,T,DEPSLON                        
      DOUBLE PRECISION SCALE,ZERO,ONE                                   
C                                                                       
      PARAMETER (ZERO = 0.0D+00)                                        
      PARAMETER (SCALE= 1.0D+00/64.0D+00)                               
      PARAMETER (ONE = 1.0D+00)                                         
C                                                                       
C-----------------------------------------------------------------------
      IERR = 0                                                          
      D(1)=DIAG(1)                                                      
      IND(1) = 1                                                        
      K = 0                                                             
      ITAG = 0                                                          
      IF (N .EQ. 1) GO TO 1001                                          
C                                                                       
      DO 100 I = 2, N                                                   
         D(I)=DIAG(I)                                                   
  100 E2(I-1) = E2IN(I)                                                 
C                                                                       
      F = ZERO                                                          
      T = ZERO                                                          
      B = DEPSLON(ONE)                                                  
      C = B *B                                                          
      B = B * SCALE                                                     
      E2(N) = ZERO                                                      
C                                                                       
      DO 290 L = 1, N                                                   
         H = ABS(D(L)) + ABS(E(L))                                      
         IF (T .GE. H) GO TO 105                                        
            T = H                                                       
            B = DEPSLON(T)                                              
            C = B * B                                                   
            B = B * SCALE                                               
  105    CONTINUE                                                       
C     .......... LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT .......... 
         M = L - 1                                                      
  110    M = M + 1                                                      
         IF (E2(M) .GT. C) GO TO 110                                    
C     .......... E2(N) IS ALWAYS ZERO, SO THERE IS AN EXIT              
C                FROM THE LOOP ..........                               
C                                                                       
         IF (M .LE. K) GO TO 125                                        
            IF (M .NE. N) E2IN(M+1) = ZERO                              
            K = M                                                       
            ITAG = ITAG + 1                                             
  125    CONTINUE                                                       
         IF (M .EQ. L) GO TO 210                                        
C                                                                       
C           ITERATE                                                     
C                                                                       
         DO 205 J = 1, 30                                               
C              .......... FORM SHIFT ..........                         
            L1 = L + 1                                                  
            S = SQRT(E2(L))                                             
            G = D(L)                                                    
            P = (D(L1) - G) / (2.0D+00 * S)                             
            R = SQRT(P*P+1.0D+00)                                       
            D(L) = S / (P + SIGN(R,P))                                  
            H = G - D(L)                                                
C                                                                       
            DO 140 I = L1, N                                            
  140       D(I) = D(I) - H                                             
C                                                                       
            F = F + H                                                   
C              .......... RATIONAL QL TRANSFORMATION ..........         
            G = D(M) + B                                                
            H = G                                                       
            S = ZERO                                                    
            DO 200 I = M-1,L,-1                                         
               P = G * H                                                
               R = P + E2(I)                                            
               E2(I+1) = S * R                                          
               S = E2(I) / R                                            
               D(I+1) = H + S * (H + D(I))                              
               G = D(I) - E2(I) / G   + B                               
               H = G * P / R                                            
  200       CONTINUE                                                    
C                                                                       
            E2(L) = S * G                                               
            D(L) = H                                                    
C              .......... GUARD AGAINST UNDERFLOW IN CONVERGENCE TEST   
            IF (H .EQ. ZERO) GO TO 210                                  
            IF (ABS(E2(L)) .LE. ABS(C/H)) GO TO 210                     
            E2(L) = H * E2(L)                                           
            IF (E2(L) .EQ. ZERO) GO TO 210                              
  205    CONTINUE                                                       
C     .......... SET ERROR -- NO CONVERGENCE TO AN                      
C                EIGENVALUE AFTER 30 ITERATIONS ..........              
      IERR = L                                                          
      GO TO 1001                                                        
C                                                                       
C           CONVERGED                                                   
C                                                                       
  210    P = D(L) + F                                                   
C           .......... ORDER EIGENVALUES ..........                     
         I = 1                                                          
         IF (L .EQ. 1) GO TO 250                                        
            IF (P .LT. D(1)) GO TO 230                                  
               I = L                                                    
C           .......... LOOP TO FIND ORDERED POSITION                    
  220          I = I - 1                                                
               IF (P .LT. D(I)) GO TO 220                               
C                                                                       
               I = I + 1                                                
               IF (I .EQ. L) GO TO 250                                  
  230       CONTINUE                                                    
            DO 240 II = L, I+1, -1                                      
               D(II) = D(II-1)                                          
               IND(II) = IND(II-1)                                      
  240       CONTINUE                                                    
C                                                                       
  250    CONTINUE                                                       
         D(I) = P                                                       
         IND(I) = ITAG                                                  
  290 CONTINUE                                                          
C                                                                       
 1001 RETURN                                                            
      END                                                               
C                                                                       
C*MODULE EIGEN   *DECK EINVIT                                           
      SUBROUTINE EINVIT(NM,N,D,E,E2,M,W,IND,Z,IERR,                     
     *                  RV1,RV2,RV3,RV4,RV6,RV7,RV8)                    
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C*                                                                      
C*    AUTHORS-                                                          
C*       THIS IS A MODIFICATION OF TINVIT FROM EISPACK EDITION 3        
C*       DATED AUGUST 1983.                                             
C*       TINVIT IS A TRANSLATION OF THE INVERSE ITERATION TECHNIQUE     
C*       IN THE ALGOL PROCEDURE TRISTURM BY PETERS AND WILKINSON.       
C*       HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 418-439(1971).
C*       THIS VERSION IS BY S. T. ELBERT (AMES LABORATORY-USDOE)        
C*                                                                      
C*    PURPOSE -                                                         
C*       THIS ROUTINE FINDS THOSE EIGENVECTORS OF A TRIDIAGONAL         
C*       SYMMETRIC MATRIX CORRESPONDING TO SPECIFIED EIGENVALUES.       
C*                                                                      
C*    METHOD -                                                          
C*       INVERSE ITERATION.                                             
C*                                                                      
C*    ON ENTRY -                                                        
C*       NM     - INTEGER                                               
C*                MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL   
C*                ARRAY PARAMETERS AS DECLARED IN THE CALLING ROUTINE   
C*                DIMENSION STATEMENT.                                  
C*       N      - INTEGER                                               
C*       D      - W.P. REAL (N)                                         
C*                CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.   
C*       E      - W.P. REAL (N)                                         
C*                CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX 
C*                IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.        
C*       E2     - W.P. REAL (N)                                         
C*                CONTAINS THE SQUARES OF CORRESPONDING ELEMENTS OF E,  
C*                WITH ZEROS CORRESPONDING TO NEGLIGIBLE ELEMENTS OF E. 
C*                E(I) IS CONSIDERED NEGLIGIBLE IF IT IS NOT LARGER THAN
C*                THE PRODUCT OF THE RELATIVE MACHINE PRECISION AND THE 
C*                SUM OF THE MAGNITUDES OF D(I) AND D(I-1).  E2(1) MUST 
C*                CONTAIN 0.0 IF THE EIGENVALUES ARE IN ASCENDING ORDER,
C*                OR 2.0 IF THE EIGENVALUES ARE IN DESCENDING ORDER.    
C*                IF TQLRAT, BISECT, TRIDIB, OR IMTQLV                  
C*                HAS BEEN USED TO FIND THE EIGENVALUES, THEIR          
C*                OUTPUT E2 ARRAY IS EXACTLY WHAT IS EXPECTED HERE.     
C*       M      - INTEGER                                               
C*                THE NUMBER OF SPECIFIED EIGENVECTORS.                 
C*       W      - W.P. REAL (M)                                         
C*                CONTAINS THE M EIGENVALUES IN ASCENDING               
C*                OR DESCENDING ORDER.                                  
C*       IND    - INTEGER (M)                                           
C*                CONTAINS IN FIRST M POSITIONS THE SUBMATRIX INDICES   
C*                ASSOCIATED WITH THE CORRESPONDING EIGENVALUES IN W -- 
C*                1 FOR EIGENVALUES BELONGING TO THE FIRST SUBMATRIX    
C*                FROM THE TOP, 2 FOR THOSE BELONGING TO THE SECOND     
C*                SUBMATRIX, ETC.                                       
C*       IERR   - INTEGER (LOGICAL UNIT NUMBER)                         
C*                LOGICAL UNIT FOR ERROR MESSAGES                       
C*                                                                      
C*    ON EXIT -                                                         
C*       ALL INPUT ARRAYS ARE UNALTERED.                                
C*       Z      - W.P. REAL (NM,M)                                      
C*                CONTAINS THE ASSOCIATED SET OF ORTHONORMAL            
C*                EIGENVECTORS. ANY VECTOR WHICH WHICH FAILS TO CONVERGE
C*                IS LEFT AS IS (BUT NORMALIZED) WHEN ITERATING STOPPED.
C*       IERR   - INTEGER                                               
C*                SET TO                                                
C*                ZERO    FOR NORMAL RETURN,                            
C*                -R      IF THE EIGENVECTOR CORRESPONDING TO THE R-TH  
C*                        EIGENVALUE FAILS TO CONVERGE IN 5 ITERATIONS. 
C*                        (ONLY LAST FAILURE TO CONVERGE IS REPORTED)   
C*                                                                      
C*       RV1, RV2, RV3, RV4, AND RV6 ARE TEMPORARY STORAGE ARRAYS.      
C*                                                                      
C*       RV1    - W.P. REAL (N)                                         
C*                DIAGONAL ELEMENTS OF U FROM LU DECOMPOSITION          
C*       RV2    - W.P. REAL (N)                                         
C*                SUPER(1)-DIAGONAL ELEMENTS OF U FROM LU DECOMPOSITION 
C*       RV3    - W.P. REAL (N)                                         
C*                SUPER(2)-DIAGONAL ELEMENTS OF U FROM LU DECOMPOSITION 
C*       RV4    - W.P. REAL (N)                                         
C*                ELEMENTS DEFINING L IN LU DECOMPOSITION               
C*       RV6    - W.P. REAL (N)                                         
C*                APPROXIMATE EIGENVECTOR                               
C*       RV7    - W.P. REAL (M,2)                                       
C*                NON-CONVERGED EIGENVECTOR RHO AND NORM                
C*       RV8    - W.P. REAL (M)                                         
C*                                                                      
C*    DIFFERENCES FROM EISPACK 3 -                                      
C*       EPS3 IS SCALED BY  EPSCAL  (ENHANCES CONVERGENCE, BUT          
C*          LOWERS ACCURACY)!                                           
C*       ONE MORE ITERATION (MINIMUM 2) IS PERFORMED AFTER CONVERGENCE  
C*          (ENHANCES ACCURACY)!                                        
C*       REPLACE LOOP WITH PYTHAG WITH SINGLE CALL TO DNRM2!            
C*       IF NOT CONVERGED, USE PERFORMANCE INDEX TO DECIDE ON ERROR     
C*          VALUE SETTING, BUT DO NOT STOP!                             
C*       L.U. FOR ERROR MESSAGES PASSED THROUGH IERR                    
C*       USE PARAMETER STATEMENTS AND GENERIC INTRINSIC FUNCTIONS       
C*       USE LEVEL 1 BLAS                                               
C*       USE IF-THEN-ELSE TO CLARIFY LOGIC                              
C*       LOOP OVER SUBSPACES MADE INTO DO LOOP.                         
C*       LOOP OVER INVERSE ITERATIONS MADE INTO DO LOOP                 
C*       ZERO ONLY REQUIRED PORTIONS OF OUTPUT VECTOR                   
C*                                                                      
C*    NOTE -                                                            
C*       QUESTIONS AND COMMENTS CONCERNING EISPACK SHOULD BE DIRECTED TO
C*       B. S. GARBOW, APPLIED MATH. DIVISION, ARGONNE NATIONAL LAB.    
C*                                                                      
C                                                                       
      LOGICAL CONVGD                                                    
C                                                                       
      INTEGER GROUP,I,IERR,ITS,M,N,NM,P,Q,R,S,SUBMAT,TAG                
      INTEGER RR,RMAX,MM                                                
      INTEGER IND(M)                                                    
C                                                                       
      DOUBLE PRECISION D(N),E(N),E2(N),W(M),Z(NM,M)                     
      DOUBLE PRECISION RV1(N),RV2(N),RV3(N),RV4(N),RV6(N)               
      DOUBLE PRECISION RV7(M,2) ! RHO AND NORM ERROR VECTOR.            
      DOUBLE PRECISION RV8(M)   ! EINICGS TEMPORARY VECTOR              
      DOUBLE PRECISION ANORM,EPS2,EPS3,EPS4,NORM,ORDER,RHO,U,UK,V       
      DOUBLE PRECISION X0,X1,XU                                         
      DOUBLE PRECISION EPSCAL,GRPTOL,HUNDRD,ONE,TEN,ZERO                
      DOUBLE PRECISION DEPSLON                                          
C                                                                       
      PARAMETER (ZERO = 0.0D+00)                                        
      PARAMETER (ONE = 1.0D+00)                                         
      PARAMETER (GRPTOL = 0.001D+00)                                    
      PARAMETER (EPSCAL = 0.5D+00)                                      
      PARAMETER (HUNDRD = 100.0D+00)                                    
      PARAMETER (TEN = 10.0D+00)                                        
C                                                                       
      LOGICAL GOPARR,DSKWRK,MASWRK,PARALL2,PARALL3                      
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
      MXSQGS = MXSQN2()         ! SEQUENTIAL GRAM-SCHMIDT THRESHOLD     
      PARALL3 = GOPARR  .AND.  M.GT.MXSQN3() ! PARALLEL N**3 LOOP       
      PARALL2 = .FALSE.         ! PARALLEL GRAM-SCHMIDT                 
C                                                                       
      CALL XCOPY(M*NM,ZERO,0,Z,1) ! ZERO VECTOR ARRAY                   
      CALL XCOPY(M*2,ZERO,0,RV7,1) ! ZERO ERROR ARRAY                   
C                                                                       
      LUEMSG = IERR                                                     
      IERR = 0                                                          
      X0 = ZERO                                                         
      UK = ZERO                                                         
      NORM = ZERO                                                       
      EPS2 = ZERO                                                       
      EPS3 = ZERO                                                       
      EPS4 = ZERO                                                       
      GROUP = 0                                                         
      TAG = 0                                                           
      ORDER = ONE - E2(1)                                               
      Q = 0                                                             
      RMAX = 0                                                          
      ICPU = ME                                                         
      NCPU = NPROC                                                      
C                                                                       
      DO 930 SUBMAT = 1, N                                              
         P = Q + 1                                                      
C                                                                       
C        .......... ESTABLISH AND PROCESS NEXT SUBMATRIX ..........     
C                                                                       
         DO 120 Q = P, N-1                                              
            IF (E2(Q+1) .EQ. ZERO) GO TO 140                            
  120    CONTINUE                                                       
         Q = N                                                          
C                                                                       
C        .......... FIND VECTORS BY INVERSE ITERATION ..........        
C                                                                       
  140    CONTINUE                                                       
         TAG = TAG + 1                                                  
         ANORM = ZERO                                                   
         S = 0                                                          
         RR = 0                                                         
C                                                                       
         DO 920 R = 1,M                                                 
            IF (IND(R) .NE. TAG) GO TO 920                              
C                                                                       
            IF (R.GT.RMAX) RMAX = R ! TRACK LARGEST R                   
            RR = RR+1 ! TRACK NON-SKIPPED R'S                           
C                                                                       
            ITS = 1                                                     
            X1 = W(R)                                                   
            IF (S .NE. 0) GO TO 510                                     
C                                                                       
C        .......... CHECK FOR ISOLATED ROOT ..........                  
C                                                                       
            XU = ONE                                                    
            IF (P .EQ. Q) THEN                                          
C                 IF NOT ASSIGNED PARALLEL TASK, SKIP OUT               
               IF (PARALL3  .AND.  ICPU.NE.MOD(RR-1,NCPU)) GOTO 910     
               RV6(P) = ONE                                             
               CONVGD = .TRUE.                                          
               GO TO 860                                                
            END IF                                                      
C                                                                       
            NORM = ABS(D(P))                                            
            DO 500 I = P+1, Q                                           
               NORM = MAX( NORM, ABS(D(I)) + ABS(E(I)) )                
  500       CONTINUE                                                    
C                                                                       
C        .......... EPS2 IS THE CRITERION FOR GROUPING,                 
C                   EPS3 REPLACES ZERO PIVOTS AND EQUAL                 
C                   ROOTS ARE MODIFIED BY EPS3,                         
C                   EPS4 IS TAKEN VERY SMALL TO AVOID OVERFLOW .........
C                                                                       
            EPS2 = GRPTOL * NORM                                        
            EPS3 = EPSCAL * DEPSLON(NORM)                               
            UK = Q - P + 1                                              
            EPS4 = UK * EPS3                                            
            UK = EPS4 / SQRT(UK)                                        
            S = P                                                       
            GROUP = 0                                                   
            GO TO 520                                                   
C                                                                       
C        .......... LOOK FOR CLOSE OR COINCIDENT ROOTS ..........       
C                                                                       
  510       IF (ABS(X1-X0) .GE. EPS2) THEN                              
C                                                                       
C                 ROOTS ARE SEPERATE                                    
C                                                                       
               GROUP = 0                                                
            ELSE                                                        
C                                                                       
C                 ROOTS ARE CLOSE                                       
C                                                                       
               GROUP = GROUP + 1                                        
               IF (ORDER * (X1 - X0) .LE. EPS3) X1 = X0 + ORDER * EPS3  
            END IF                                                      
C                                                                       
C        .......... ELIMINATION WITH INTERCHANGES AND                   
C                   INITIALIZATION OF VECTOR ..........                 
C                                                                       
  520       CONTINUE                                                    
C                                                                       
C           WHEN NUMBER OF CLUSTERED EIGENVALUES CROSSES THE THRESHOLD, 
C           TURN OFF N**3 PARALLELIZATION,                              
C           AND USE PARALLEL ITERATED GRAM-SCHMIDT.                     
C           ELSE LOAD BALANCE N**3 LOOP,                                
C           WITH SEQUENTIAL MODIFIED GRAM-SCHMIDT.                      
C                                                                       
C           THERE IS NO SWITCHING BACK FROM N**2 TO N**3 PARALLEL-      
C           IZATION BECAUSE IT IS AWKWARD TO IMPLEMENT DUE TO GENERALLY 
C           NON-CONTIGUOUS TAGS.  MOREOVER, IN MOST QM CASES, THERE IS  
C           ONLY ONE LARGE CLUSTER SPANNING ALMOST THE ENTIRE SPACE.    
C                                                                       
            IF (PARALL3) THEN                                           
               IF (GROUP.GT.MXSQGS) THEN ! SWITCH TO PARALLEL N**2 LOOP 
                  MM = RMAX                                             
C                      CURRENT EIGENVECTOR, HAS NOT BEEN COMPUTED YET   
                  IF (MM.EQ.R) MM = MM - 1                              
C                      SUM UP COMPUTED EIGENVECTORS -Z- AND ERRORS -RV7-
                  CALL DDI_GSUMF(21020,Z,NM*MM)                         
                  CALL DDI_GSUMF(21021,RV7,2*MM)                        
                  PARALL3 = .FALSE.                                     
                  PARALL2 = .TRUE.                                      
               ELSE                                                     
                  IF (ICPU.NE.MOD(RR-GROUP-1,NCPU)) GOTO 910            
               ENDIF                                                    
            ENDIF                                                       
C                                                                       
            U = D(P) - X1                                               
            V = E(P+1)                                                  
            RV6(P) = UK                                                 
            DO 550 I = P+1, Q                                           
               RV6(I) = UK                                              
               IF (ABS(E(I)) .GT. ABS(U)) THEN                          
C                                                                       
C                 EXCHANGE ROWS BEFORE ELIMINATION                      
C                                                                       
C                  *** WARNING -- A DIVIDE CHECK MAY OCCUR HERE IF      
C                      E2 ARRAY HAS NOT BEEN SPECIFIED CORRECTLY .......
C                                                                       
                  XU = U / E(I)                                         
                  RV4(I) = XU                                           
                  RV1(I-1) = E(I)                                       
                  RV2(I-1) = D(I) - X1                                  
                  RV3(I-1) = E(I+1)                                     
                  U = V - XU * RV2(I-1)                                 
                  V = -XU * RV3(I-1)                                    
C                                                                       
               ELSE                                                     
C                                                                       
C                    STRAIGHT ELIMINATION                               
C                                                                       
                  XU = E(I) / U                                         
                  RV4(I) = XU                                           
                  RV1(I-1) = U                                          
                  RV2(I-1) = V                                          
                  RV3(I-1) = ZERO                                       
                  U = D(I) - X1 - XU * V                                
                  V = E(I+1)                                            
               END IF                                                   
  550       CONTINUE                                                    
C                                                                       
            IF (ABS(U) .LE. EPS3) U = EPS3                              
            RV1(Q) = U                                                  
            RV2(Q) = ZERO                                               
            RV3(Q) = ZERO                                               
C                                                                       
C                                                                       
C                                                                       
C              DO INVERSE ITERATIONS                                    
C                                                                       
            CONVGD = .FALSE.                                            
            DO 800 ITS = 1, 5                                           
               IF (ITS .EQ. 1) GO TO 600                                
C                                                                       
C                    .......... FORWARD SUBSTITUTION ..........         
C                                                                       
                  IF (NORM .EQ. ZERO) THEN                              
                     RV6(S) = EPS4                                      
                     S = S + 1                                          
                     IF (S .GT. Q) S = P                                
                  ELSE                                                  
                     XU = EPS4 / NORM                                   
                     CALL XSCAL(Q-P+1, XU, RV6(P), 1)                   
                  END IF                                                
C                                                                       
C                     ... ELIMINATION OPERATIONS ON NEXT VECTOR         
C                                                                       
                  DO 590 I = P+1, Q                                     
                     U = RV6(I)                                         
C                                                                       
C                         IF RV1(I-1) .EQ. E(I), A ROW INTERCHANGE      
C                         WAS PERFORMED EARLIER IN THE                  
C                         TRIANGULARIZATION PROCESS ..........          
C                                                                       
C-STE-               IF (RV1(I-1) .EQ. E(I)) THEN                       
                     IF (RV1(I-1) .EQ. E(I) .AND. RV4(I) .NE. ONE) THEN 
                        U = RV6(I-1)                                    
                        RV6(I-1) = RV6(I)                               
                     ELSE                                               
                        U = RV6(I)                                      
                     END IF                                             
                     RV6(I) = U - RV4(I) * RV6(I-1)                     
  590             CONTINUE                                              
  600          CONTINUE                                                 
C                                                                       
C           .......... BACK SUBSTITUTION                                
C                                                                       
               RV6(Q) = RV6(Q) / RV1(Q)                                 
               V = U                                                    
               U = RV6(Q)                                               
               NORM = ABS(U)                                            
               DO 620 I = Q-1, P, -1                                    
                  RV6(I) = (RV6(I) - U * RV2(I) - V * RV3(I)) / RV1(I)  
                  V = U                                                 
                  U = RV6(I)                                            
                  NORM = NORM + ABS(U)                                  
  620          CONTINUE                                                 
C                                                                       
               IF (GROUP.GT.0) THEN ! CLUSTERED EIGENVALUES             
C                                                                       
C   ..... ORTHOGONALIZE WITH RESPECT TO PREVIOUS MEMBERS OF GROUP ..... 
C              USE PARALLEL ITERATED CLASSICAL GRAM-SCHMIDT             
C                 OR SEQUENTIAL MODIFIED GRAM-SCHMIDT                   
C                                                                       
                  IF (PARALL2.AND.GROUP.GT.MXSQGS) THEN                 
                     CALL EINICGS(Q-P+1,GROUP,Z(P,1),NM,RV6(P),R,TAG,   
     *                            IND,RV8,2,.TRUE.,ICPU,NCPU)           
                  ELSE                                                  
                     CALL EINMGS(Q-P+1,GROUP,Z(P,1),NM,RV6(P),R,TAG,IND)
                  ENDIF                                                 
                  NORM = XASUM(Q-P+1, RV6(P), 1)                        
               ENDIF                                                    
C                                                                       
               IF (CONVGD) GO TO 840                                    
               IF (NORM .GE. ONE) CONVGD = .TRUE.                       
  800       CONTINUE                                                    
C                                                                       
C        .......... NORMALIZE SO THAT SUM OF SQUARES IS                 
C                   1 AND EXPAND TO FULL ORDER ..........               
C                                                                       
  840       CONTINUE                                                    
C                                                                       
            XU = ONE / XNRM2(Q-P+1,RV6(P),1)                            
C                                                                       
  860       CONTINUE                                                    
C                                                                       
            CALL XAXPY(Q-P+1,XU,RV6(P),1,Z(P,R),1)                      
C                                                                       
            IF (.NOT.CONVGD) THEN ! ERRORS ARE DETERMINED IN THE END.   
               CALL EINPI(RHO,Q-P+1,X1,D(P),E(P),Z(P,R),ANORM)          
               RV7(R,1) = RHO                                           
               RV7(R,2) = NORM                                          
            ENDIF                                                       
C                                                                       
 910        CONTINUE                                                    
            X0 = X1                                                     
  920    CONTINUE                                                       
C                                                                       
         IF (Q .EQ. N) GO TO 940                                        
  930 CONTINUE                                                          
  940 CONTINUE                                                          
C                                                                       
      IF (PARALL3) THEN                                                 
         CALL DDI_GSUMF(21022,Z,NM*M)                                   
         CALL DDI_GSUMF(21023,RV7,M*2)                                  
      ENDIF                                                             
C                                                                       
C     CHECK ACCUMULATED ERRORS.                                         
C                                                                       
      DO R=1,M                                                          
         RHO = RV7(R,1)                                                 
         NORM = RV7(R,2)                                                
         IF (RHO.GE.TEN.AND.LUEMSG.GT.0.AND.MASWRK)                     
     *        WRITE(LUEMSG,001) R,NORM,RHO                              
C                                                                       
C        *** SET ERROR -- NON-CONVERGED EIGENVECTOR ..........          
C                                                                       
         IF (RHO.GT.HUNDRD) IERR = -R                                   
      ENDDO                                                             
C                                                                       
      RETURN                                                            
C                                                                       
    1 FORMAT(1X,'WARNING: EIGENVECTOR ROUTINE -EINVIT- DID NOT',        
     *          ' CONVERGE FOR VECTOR',I6/                              
     *       1X,'NORM =',1P,E10.2,' PERFORMANCE INDEX =',E10.2/         
     *       1X,'A HALT WILL OCCUR IF THE PERF.IND. EXCEEDS 100.0')     
      END                                                               
C                                                                       
C*MODULE EIGEN   *DECK EINMGS                                           
      SUBROUTINE EINMGS(M,N,Q,LDQ,R,IR,ITAG,IND)                        
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C                                                                       
C     ARGUMENTS                                                         
C                                                                       
      INTEGER M                 ![IN] VECTOR LENGTH.                    
      INTEGER N                 ![IN] NUMBER OF ORTHOGONAL VECTORS.     
      DOUBLE PRECISION Q(LDQ,*) ![IN] ORTHOGONALIZED VECTORS.           
      DOUBLE PRECISION R(M)     ![IN,OUT] VECTOR TO ORTHOGONALIZE.      
      INTEGER IR                ![IN] VECTOR INDEX.                     
      INTEGER ITAG              ![IN] CLUSTER TAG.                      
      INTEGER IND(*)            ![IN] CLUSTER TAGS BY VECTOR INDICES.   
C                                                                       
C     MODIFIED GRAM-SCHMIDT.                                            
C     SEQUENTIAL GRAM-SCHMIDT ORTHOGONALIZATION ROUTINE USED BY EINVIT. 
C                                                                       
      J = IR                                                            
      DO JJ=1,N                                                         
 100     J = J - 1                                                      
         IF (IND(J).NE.ITAG) GOTO 100                                   
         QR = XDOT(M,Q(1,J),1,R,1)                                      
         CALL XAXPY(M,-QR,Q(1,J),1,R,1)                                 
      ENDDO                                                             
      RETURN                                                            
      END                                                               
C                                                                       
C*MODULE EIGEN   *DECK EINICGS                                          
      SUBROUTINE EINICGS(M,N,Q,LDQ,R,IR,ITAG,IND,                       
     *                   P,NITER,PARALL,IPROC,NPROC)                    
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
C     ARGUMENTS                                                         
C                                                                       
      INTEGER M                 ![IN] VECTOR LENGTH.                    
      INTEGER N                 ![IN] NUMBER OF ORTHOGONAL VECTORS.     
      DOUBLE PRECISION Q(LDQ,*) ![IN] ORTHOGONALIZED VECTORS.           
      INTEGER LDQ               ![IN] Q LEADING DIMENSION.              
      DOUBLE PRECISION R(M)     ![IN,OUT] VECTOR TO ORTHOGONALIZE.      
      INTEGER IR                ![IN] VECTOR INDEX.                     
      INTEGER ITAG              ![IN] CLUSTER TAG.                      
      INTEGER IND(*)            ![IN] CLUSTER TAGS BY VECTOR INDICES.   
      DOUBLE PRECISION P(N)     ![OUT] TEMPORARY VECTOR.                
      INTEGER NITER             ![IN] NUMBER OF ITERATIONS              
      LOGICAL PARALL            ![IN] PARALLEL FLAG.                    
      INTEGER IPROC             ![IN] PROCESS.                          
      INTEGER NPROC             ![IN] NUMBER OF PROCESSES.              
C                                                                       
      PARAMETER (ZERO=0.0D+00)                                          
C                                                                       
C     ITERATED CLASSICAL GRAM-SCHMIDT.                                  
C     PARALLELIZABLE GRAM-SCHMIDT ORTHOGONALIZATION ROUTINE             
C     USED BY EINVIT.  TWO ITERATIONS SHOULD BE ENOUGH TO ACHIEVE       
C     RESULTS ON PAR OR BETTER THAN MODIFIED GRAM-SCHMIDT RESULTS.      
C                                                                       
C     SEE F.J. LINGEN,                                                  
C     EFFICIENT GRAM-SCHMIDT ORTHONORMALISATION ON PARALLEL COMPUTERS.  
C     COMMUNICATIONS IN NUMERICAL METHODS IN ENGINEERING 16, 57 (2000). 
C                                                                       
      DO ITER=1,NITER                                                   
         J = IR                                                         
         JP = 1                                                         
         DO 100 JJ=1,N ! P = QR                                         
 110        J = J - 1                                                   
            IF (IND(J).NE.ITAG) GOTO 110                                
            IF (PARALL.AND.IPROC.NE.MOD(JJ-1,NPROC)) GOTO 100           
            P(JP) = XDOT(M,Q(1,J),1,R,1)                                
            JP = JP + 1                                                 
 100     ENDDO                                                          
C                                                                       
         IF (PARALL.AND.IPROC.NE.0) CALL XCOPY(M,ZERO,0,R,1)            
C                                                                       
         J = IR                                                         
         JP = 1                                                         
         DO 200 JJ=1,N ! R = R - PQ                                     
 210        J = J - 1                                                   
            IF (IND(J).NE.ITAG) GOTO 210                                
            IF (PARALL.AND.IPROC.NE.MOD(JJ-1,NPROC)) GOTO 200           
            CALL XAXPY(M,-P(JP),Q(1,J),1,R,1)                           
            JP = JP + 1                                                 
 200     ENDDO                                                          
C                                                                       
         IF (PARALL) CALL DDI_GSUMF(21031,R,M)                          
C                                                                       
      ENDDO                                                             
      RETURN                                                            
      END                                                               
C                                                                       
C*MODULE EIGEN   *DECK EINPI                                            
      SUBROUTINE EINPI(ESTPI1,N,EVAL,D,E,X,ANORM)                       
C*                                                                      
C*    AUTHOR -                                                          
C*       STEPHEN T. ELBERT (AMES LABORATORY-USDOE) DATE: 5 DEC 1986     
C*                                                                      
C*    PURPOSE -                                                         
C*       EVALUATE SYMMETRIC TRIDIAGONAL MATRIX PERFORMANCE INDEX        
C*       *        *         *                  *           *            
C*       FOR 1 EIGENVECTOR                                              
C*           *                                                          
C*                                                                      
C*    METHOD -                                                          
C*       THIS ROUTINE FORMS THE 1-NORM OF THE RESIDUAL MATRIX A*X-X*EVAL
C*       WHERE  A  IS A SYMMETRIC TRIDIAGONAL MATRIX STORED             
C*       IN THE DIAGONAL (D) AND SUB-DIAGONAL (E) VECTORS, EVAL IS THE  
C*       EIGENVALUE OF AN EIGENVECTOR OF  A,  NAMELY  X.                
C*       THIS NORM IS SCALED BY MACHINE ACCURACY FOR THE PROBLEM SIZE.  
C*       ALL NORMS APPEARING IN THE COMMENTS BELOW ARE 1-NORMS.         
C*                                                                      
C*    ON ENTRY -                                                        
C*       N      - INTEGER                                               
C*                THE ORDER OF THE MATRIX  A.                           
C*       EVAL   - W.P. REAL                                             
C*                THE EIGENVALUE CORRESPONDING TO VECTOR  X.            
C*       D      - W.P. REAL (N)                                         
C*                THE DIAGONAL VECTOR OF  A.                            
C*       E      - W.P. REAL (N)                                         
C*                THE SUB-DIAGONAL VECTOR OF  A.                        
C*       X      - W.P. REAL (N)                                         
C*                AN EIGENVECTOR OF  A.                                 
C*       ANORM  - W.P. REAL                                             
C*                THE NORM OF  A  IF IT HAS BEEN PREVIOUSLY COMPUTED.   
C*                                                                      
C*    ON EXIT -                                                         
C*       ANORM  - W.P. REAL                                             
C*                THE NORM OF  A, COMPUTED IF INITIALLY ZERO.           
C*       ESTPI1 - W.P. REAL                                             
C*          !!A*X-X*EVAL!! / (DEPSLON(10*N)*!!A!!*!!X!!);               
C*          WHERE DEPSLON(X) IS THE SMALLEST NUMBER SUCH THAT           
C*             X + DEPSLON(X) .NE. X                                    
C*                                                                      
C*          ESTPI1 .LT. 1 == SATISFACTORY PERFORMANCE                   
C*                 .GE. 1 AND .LE. 100 == MARGINAL PERFORMANCE          
C*                 .GT. 100 == POOR PERFORMANCE                         
C*          (SEE LECT. NOTES IN COMP. SCI. VOL.6 PP 124-125)            
C                                                                       
      DOUBLE PRECISION ESTPI1,ANORM,EVAL,RNORM,SIZE,XNORM               
      DOUBLE PRECISION D(N), E(N), X(N)                                 
      DOUBLE PRECISION DEPSLON, ONE, ZERO                               
C                                                                       
      PARAMETER (ZERO = 0.0D+00)                                        
      PARAMETER (ONE = 1.0D+00)                                         
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
      ESTPI1 = ZERO                                                     
      IF( N .LE. 1 ) RETURN                                             
      SIZE = 10 * N                                                     
      IF (ANORM .EQ. ZERO) THEN                                         
C                                                                       
C              COMPUTE NORM OF  A                                       
C                                                                       
         ANORM = MAX( ABS(D(1)) + ABS(E(2))                             
     *               ,ABS(D(N)) + ABS(E(N)))                            
         DO 110 I = 2, N-1                                              
            ANORM = MAX( ANORM, ABS(E(I))+ABS(D(I))+ABS(E(I+1)))        
  110    CONTINUE                                                       
         IF(ANORM .EQ. ZERO) ANORM = ONE                                
      END IF                                                            
C                                                                       
C           COMPUTE NORMS OF RESIDUAL AND EIGENVECTOR                   
C                                                                       
      XNORM = ABS(X(1)) + ABS(X(N))                                     
      RNORM = ABS( (D(1)-EVAL)*X(1) + E(2)*X(2))                        
     *       +ABS( (D(N)-EVAL)*X(N) + E(N)*X(N-1))                      
      DO 120 I = 2, N-1                                                 
         XNORM = XNORM + ABS(X(I))                                      
         RNORM = RNORM + ABS(E(I)*X(I-1) + (D(I)-EVAL)*X(I)             
     *                       + E(I+1)*X(I+1))                           
  120 CONTINUE                                                          
C                                                                       
      ESTPI1 = RNORM / (DEPSLON(SIZE)*ANORM*XNORM)                      
      RETURN                                                            
      END                                                               
C                                                                       
C*MODULE EIGEN   *DECK ETRBK3                                           
      SUBROUTINE ETRBK3(NM,N,NV,A,M,Z)                                  
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C*                                                                      
C*    AUTHORS-                                                          
C*       THIS IS A MODIFICATION OF ROUTINE TRBAK3 FROM EISPACK EDITION 3
C*       DATED AUGUST 1983.                                             
C*       EISPACK TRBAK3 IS A TRANSLATION OF THE ALGOL PROCEDURE TRBAK3, 
C*       NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
C*       HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C*       THIS VERSION IS BY S. T. ELBERT (AMES LABORATORY-USDOE)        
C*                                                                      
C*    PURPOSE -                                                         
C*       THIS ROUTINE FORMS THE EIGENVECTORS OF A REAL SYMMETRIC        
C*       MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING         
C*       SYMMETRIC TRIDIAGONAL MATRIX DETERMINED BY  ETRED3.            
C*                                                                      
C*    METHOD -                                                          
C*       THE CALCULATION IS CARRIED OUT BY FORMING THE MATRIX PRODUCT   
C*          Q*Z                                                         
C*       WHERE  Q  IS A PRODUCT OF THE ORTHOGONAL SYMMETRIC MATRICES    
C*                Q = PROD(I)[1 - U(I)*.TRANSPOSE.U(I)*H(I)]            
C*       U  IS THE AUGMENTED SUB-DIAGONAL ROWS OF  A  AND               
C*       Z  IS THE SET OF EIGENVECTORS OF THE TRIDIAGONAL               
C*       MATRIX  F  WHICH WAS FORMED FROM THE ORIGINAL SYMMETRIC        
C*       MATRIX  C  BY THE SIMILARITY TRANSFORMATION                    
C*                F = Q(TRANSPOSE) C Q                                  
C*       NOTE THAT ETRBK3 PRESERVES VECTOR EUCLIDEAN NORMS.             
C*                                                                      
C*                                                                      
C*    COMPLEXITY -                                                      
C*       M*N**2                                                         
C*                                                                      
C*    ON ENTRY-                                                         
C*       NM     - INTEGER                                               
C*                MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL   
C*                ARRAY PARAMETERS AS DECLARED IN THE CALLING ROUTINE   
C*                DIMENSION STATEMENT.                                  
C*       N      - INTEGER                                               
C*                THE ORDER OF THE MATRIX  A.                           
C*       NV     - INTEGER                                               
C*                MUST BE SET TO THE DIMENSION OF THE ARRAY  A  AS      
C*                DECLARED IN THE CALLING ROUTINE DIMENSION STATEMENT.  
C*       A      - W.P. REAL (NV)                                        
C*                CONTAINS INFORMATION ABOUT THE ORTHOGONAL             
C*                TRANSFORMATIONS USED IN THE REDUCTION BY  ETRED3  IN  
C*                ITS FIRST  NV = N*(N+1)/2 POSITIONS.                  
C*       M      - INTEGER                                               
C*                THE NUMBER OF EIGENVECTORS TO BE BACK TRANSFORMED.    
C*       Z      - W.P REAL (NM,M)                                       
C*                CONTAINS THE EIGENVECTORS TO BE BACK TRANSFORMED      
C*                IN ITS FIRST M COLUMNS.                               
C*                                                                      
C*    ON EXIT-                                                          
C*       Z      - W.P. REAL (NM,M)                                      
C*                CONTAINS THE TRANSFORMED EIGENVECTORS                 
C*                IN ITS FIRST M COLUMNS.                               
C*                                                                      
C*    DIFFERENCES WITH EISPACK 3 -                                      
C*       THE TWO INNER LOOPS ARE REPLACED BY DDOT AND DAXPY.            
C*       MULTIPLICATION USED INSTEAD OF DIVISION TO FIND S.             
C*       OUTER LOOP RANGE CHANGED FROM 2,N TO 3,N.                      
C*       ADDRESS POINTERS FOR  A  SIMPLIFIED.                           
C*                                                                      
C*    NOTE -                                                            
C*       QUESTIONS AND COMMENTS CONCERNING EISPACK SHOULD BE DIRECTED TO
C*       B. S. GARBOW, APPLIED MATH. DIVISION, ARGONNE NATIONAL LAB.    
C                                                                       
      INTEGER I,II,IM1,IZ,J,J1,J2,M,N,NM,NV                             
C                                                                       
      DOUBLE PRECISION A(NV),Z(NM,M)                                    
      DOUBLE PRECISION H,S,ZERO                                         
C                                                                       
      PARAMETER (ZERO = 0.0D+00)                                        
C                                                                       
      LOGICAL GOPARR,DSKWRK,MASWRK,PARALL3                              
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
      IF (M .EQ. 0) RETURN                                              
      IF (N .LE. 2) RETURN                                              
C                                                                       
      PARALL3 = GOPARR  .AND.  M.GT.MXSQN3()                            
      ICPU = ME                                                         
      NCPU = NPROC                                                      
C                                                                       
C     DECOMPOSE Z(NM,M)                                                 
      CALL BLK_DECOMP(M,PARALL3,ICPU,NCPU,J1,J2)                        
C                                                                       
C     ZERO SKIPPED ARRAY SEGMENTS                                       
      CALL XCOPY((J1-1)*NM,0.0D+00,0,Z(1,1),1)                          
      CALL XCOPY((M-J2)*NM,0.0D+00,0,Z(1,J2+1),1)                       
C                                                                       
      DO 130 J=J1,J2                                                    
         II=3                                                           
         DO 140 I=3,N                                                   
            IZ=II+1                                                     
            II=II+I                                                     
            H = A(II)                                                   
            IF (H .EQ. ZERO) GO TO 140                                  
            IM1 = I - 1                                                 
            S = -( XDOT(IM1,A(IZ),1,Z(1,J),1) * H) * H                  
            CALL XAXPY(IM1,S,A(IZ),1,Z(1,J),1)                          
 140     CONTINUE                                                       
 130  CONTINUE                                                          
C                                                                       
      IF(PARALL3) CALL DDI_GSUMF(21025,Z,NM*M)                          
C                                                                       
      RETURN                                                            
      END                                                               
C                                                                       
C         END OF -EVVRSP- PACKAGE, AND START OF -GIVEIS- PACKAGE        
C                                                                       
C*MODULE EIGEN   *DECK GIVEIS                                           
      SUBROUTINE GIVEIS(N,NVECT,NV,A,B,INDB,ROOT,VECT,IERR)             
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      DIMENSION A(*),B(N,8),INDB(N),ROOT(N),VECT(NV,NVECT)              
C                                                                       
C     EISPACK-BASED SUBSTITUTE FOR QCPE ROUTINE GIVENS.                 
C     FINDS ALL EIGENVALUES AND SOME EIGENVECTORS OF A REAL SYMMETRIC   
C     MATRIX.   AUTHOR.. C. MOLER AND D. SPANGLER, N.R.C.C., 4/1/79.    
C                                                                       
C     INPUT..                                                           
C     N     = ORDER OF MATRIX .                                         
C     NVECT = NUMBER OF VECTORS DESIRED.  0 .LE. NVECT .LE. N .         
C     NV    = LEADING DIMENSION OF VECT .                               
C     A     = INPUT MATRIX, COLUMNS OF THE UPPER TRIANGLE PACKED INTO   
C             LINEAR ARRAY OF DIMENSION N*(N+1)/2 .                     
C     B     = SCRATCH ARRAY, 8*N ELEMENTS (NOTE THIS IS MORE THAN       
C             PREVIOUS VERSIONS OF GIVENS.)                             
C    IND    = INDEX ARRAY OF N ELEMENTS                                 
C                                                                       
C     OUTPUT..                                                          
C     A       DESTROYED .                                               
C     ROOT  = ALL EIGENVALUES, ROOT(1) .LE. ... .LE. ROOT(N) .          
C             (FOR OTHER ORDERINGS, SEE BELOW.)                         
C     VECT  = EIGENVECTORS FOR ROOT(1),..., ROOT(NVECT) .               
C     IERR  = 0 IF NO ERROR DETECTED,                                   
C           = K IF ITERATION FOR K-TH EIGENVALUE FAILED,                
C           = -K IF ITERATION FOR K-TH EIGENVECTOR FAILED.              
C             (FAILURES SHOULD BE VERY RARE.  CONTACT MOLER.)           
C                                                                       
C     CALLS MODIFIED EISPACK ROUTINES TRED3B, IMTQLV, TINVTB, AND       
C     TRBK3B.  THE ROUTINES TRED3B, TINVTB, AND TRBK3B.                 
C     THE ORIGINAL EISPACK ROUTINES TRED3, TINVIT, AND TRBAK3           
C     WERE MODIFIED BY THE INTRODUCTION OF TWO ROUTINES FROM THE        
C     BLAS LIBRARY - DDOT AND DAXPY.                                    
C                                                                       
C         IF TINVIT FAILS TO CONVERGE, TQL2 IS CALLED                   
C                                                                       
C         SEE EISPACK USERS GUIDE, B. T. SMITH ET AL, SPRINGER-VERLAG   
C     LECTURE NOTES IN COMPUTER SCIENCE, VOL. 6, 2-ND EDITION, 1976 .   
C     NOTE THAT IMTQLV AND TINVTB HAVE INTERNAL MACHINE                 
C     DEPENDENT CONSTANTS.                                              
C                                                                       
      DATA ONE, ZERO /1.0D+00, 0.0D+00/                                 
      CALL TRED3B(N,(N*N+N)/2,A,B(1,1),B(1,2),B(1,3))                   
      CALL IMTQLV(N,B(1,1),B(1,2),B(1,3),ROOT,INDB,IERR,B(1,4))         
      IF (IERR .NE. 0) RETURN                                           
C                                                                       
C     TO REORDER ROOTS...                                               
C     K = N/2                                                           
C     B(1,3) = 2.0D+00                                                  
C     DO 50 I = 1, K                                                    
C        J = N+1-I                                                      
C        T = ROOT(I)                                                    
C        ROOT(I) = ROOT(J)                                              
C        ROOT(J) = T                                                    
C 50  CONTINUE                                                          
C                                                                       
      IF (NVECT .LE. 0) RETURN                                          
      CALL TINVTB(NV,N,B(1,1),B(1,2),B(1,3),NVECT,ROOT,INDB,VECT,IERR,  
     +     B(1,4),B(1,5),B(1,6),B(1,7),B(1,8))                          
      IF (IERR .EQ. 0) GO TO 160                                        
C                                                                       
C      IF INVERSE ITERATION GIVES AN ERROR IN DETERMINING THE           
C      EIGENVECTORS, TRY THE QL ALGORITHM IF ALL THE EIGENVECTORS       
C      ARE DESIRED.                                                     
C                                                                       
      IF (NVECT .NE. N) RETURN                                          
      DO 120 I = 1, NVECT                                               
      DO 100 J = 1, N                                                   
      VECT(I,J) = ZERO                                                  
  100 CONTINUE                                                          
      VECT(I,I) = ONE                                                   
  120 CONTINUE                                                          
      CALL TQL2 (NV,N,B(1,1),B(1,2),VECT,IERR)                          
      DO 140 I = 1, NVECT                                               
      ROOT(I) = B(I,1)                                                  
  140 CONTINUE                                                          
      IF (IERR .NE. 0) RETURN                                           
  160 CALL TRBK3B(NV,N,(N*N+N)/2,A,NVECT,VECT)                          
      RETURN                                                            
      END                                                               
C                                                                       
C*MODULE EIGEN   *DECK TRED3B                                           
      SUBROUTINE TRED3B(N,NV,A,D,E,E2)                                  
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      DIMENSION A(NV),D(N),E(N),E2(N)                                   
C                                                                       
C     THIS ROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED3,       
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.   
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).   
C                                                                       
C     THIS ROUTINE REDUCES A REAL SYMMETRIC MATRIX, STORED AS           
C     A ONE-DIMENSIONAL ARRAY, TO A SYMMETRIC TRIDIAGONAL MATRIX        
C     USING ORTHOGONAL SIMILARITY TRANSFORMATIONS.                      
C                                                                       
C     ON INPUT-                                                         
C                                                                       
C        N IS THE ORDER OF THE MATRIX,                                  
C                                                                       
C        NV MUST BE SET TO THE DIMENSION OF THE ARRAY PARAMETER A       
C          AS DECLARED IN THE CALLING ROUTINE DIMENSION STATEMENT,      
C                                                                       
C        A CONTAINS THE LOWER TRIANGLE OF THE REAL SYMMETRIC            
C          INPUT MATRIX, STORED ROW-WISE AS A ONE-DIMENSIONAL           
C          ARRAY, IN ITS FIRST N*(N+1)/2 POSITIONS.                     
C                                                                       
C     ON OUTPUT-                                                        
C                                                                       
C        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL                    
C          TRANSFORMATIONS USED IN THE REDUCTION,                       
C                                                                       
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX,    
C                                                                       
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL         
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO,      
C                                                                       
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.    
C          E2 MAY COINCIDE WITH E IF THE SQUARES ARE NOT NEEDED.        
C                                                                       
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,        
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY         
C                                                                       
C     ------------------------------------------------------------------
C                                                                       
C     ********** FOR I=N STEP -1 UNTIL 1 DO -- **********               
      DO 300 II = 1, N                                                  
      I = N + 1 - II                                                    
      L = I - 1                                                         
      IZ = (I * L) / 2                                                  
      H = 0.0D+00                                                       
      SCALE = 0.0D+00                                                   
      IF (L .LT. 1) GO TO 120                                           
C     ********** SCALE ROW (ALGOL TOL THEN NOT NEEDED) **********       
      DO 100 K = 1, L                                                   
      IZ = IZ + 1                                                       
      D(K) = A(IZ)                                                      
      SCALE = SCALE + ABS(D(K))                                         
  100 CONTINUE                                                          
C                                                                       
      IF (SCALE .NE. 0.0D+00) GO TO 140                                 
  120 E(I) = 0.0D+00                                                    
      E2(I) = 0.0D+00                                                   
      GO TO 280                                                         
C                                                                       
  140 DO 160 K = 1, L                                                   
      D(K) = D(K) / SCALE                                               
      H = H + D(K) * D(K)                                               
  160 CONTINUE                                                          
C                                                                       
      E2(I) = SCALE * SCALE * H                                         
      F = D(L)                                                          
      G = -SIGN(SQRT(H),F)                                              
      E(I) = SCALE * G                                                  
      H = H - F * G                                                     
      D(L) = F - G                                                      
      A(IZ) = SCALE * D(L)                                              
      IF (L .EQ. 1) GO TO 280                                           
      F = 0.0D+00                                                       
C                                                                       
      JK = 1                                                            
      DO 220 J = 1, L                                                   
      JM1 = J - 1                                                       
      DT = D(J)                                                         
      G = 0.0D+00                                                       
C     ********** FORM ELEMENT OF A*U **********                         
      IF (JM1 .EQ. 0) GO TO 200                                         
      DO 180 K = 1, JM1                                                 
      E(K) = E(K) + DT * A(JK)                                          
      G = G + D(K) * A(JK)                                              
      JK = JK + 1                                                       
  180 CONTINUE                                                          
  200 E(J) = G + A(JK) * DT                                             
      JK = JK + 1                                                       
C     ********** FORM ELEMENT OF P **********                           
  220 CONTINUE                                                          
      F = 0.0D+00                                                       
      DO 240 J = 1, L                                                   
      E(J) = E(J) / H                                                   
      F = F + E(J) * D(J)                                               
  240 CONTINUE                                                          
C                                                                       
      HH = F / (H + H)                                                  
      JK = 0                                                            
C     ********** FORM REDUCED A **********                              
      DO 260 J = 1, L                                                   
      F = D(J)                                                          
      G = E(J) - HH * F                                                 
      E(J) = G                                                          
C                                                                       
      DO 260 K = 1, J                                                   
      JK = JK + 1                                                       
      A(JK) = A(JK) - F * E(K) - G * D(K)                               
  260 CONTINUE                                                          
C                                                                       
  280 D(I) = A(IZ+1)                                                    
      A(IZ+1) = SCALE * SQRT(H)                                         
  300 CONTINUE                                                          
C                                                                       
      RETURN                                                            
      END                                                               
C                                                                       
C*MODULE EIGEN   *DECK IMTQLV                                           
      SUBROUTINE IMTQLV(N,D,E,E2,W,IND,IERR,RV1)                        
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      INTEGER TAG                                                       
      DOUBLE PRECISION MACHEP                                           
      DIMENSION D(N),E(N),E2(N),W(N),RV1(N),IND(N)                      
C                                                                       
C     THIS ROUTINE IS A VARIANT OF  IMTQL1  WHICH IS A TRANSLATION OF   
C     ALGOL PROCEDURE IMTQL1, NUM. MATH. 12, 377-383(1968) BY MARTIN AND
C     WILKINSON, AS MODIFIED IN NUM. MATH. 15, 450(1970) BY DUBRULLE.   
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 241-248(1971).   
C                                                                       
C     THIS ROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC TRIDIAGONAL     
C     MATRIX BY THE IMPLICIT QL METHOD AND ASSOCIATES WITH THEM         
C     THEIR CORRESPONDING SUBMATRIX INDICES.                            
C                                                                       
C     ON INPUT-                                                         
C                                                                       
C        N IS THE ORDER OF THE MATRIX,                                  
C                                                                       
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX,          
C                                                                       
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX        
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY,               
C                                                                       
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.    
C          E2(1) IS ARBITRARY.                                          
C                                                                       
C     ON OUTPUT-                                                        
C                                                                       
C        D AND E ARE UNALTERED,                                         
C                                                                       
C        ELEMENTS OF E2, CORRESPONDING TO ELEMENTS OF E REGARDED        
C          AS NEGLIGIBLE, HAVE BEEN REPLACED BY ZERO CAUSING THE        
C          MATRIX TO SPLIT INTO A DIRECT SUM OF SUBMATRICES.            
C          E2(1) IS ALSO SET TO ZERO,                                   
C                                                                       
C        W CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN          
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND          
C          ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE            
C          THE SMALLEST EIGENVALUES,                                    
C                                                                       
C        IND CONTAINS THE SUBMATRIX INDICES ASSOCIATED WITH THE         
C          CORRESPONDING EIGENVALUES IN W -- 1 FOR EIGENVALUES          
C          BELONGING TO THE FIRST SUBMATRIX FROM THE TOP,               
C          2 FOR THOSE BELONGING TO THE SECOND SUBMATRIX, ETC.,         
C                                                                       
C        IERR IS SET TO                                                 
C          ZERO       FOR NORMAL RETURN,                                
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN               
C                     DETERMINED AFTER 30 ITERATIONS,                   
C                                                                       
C        RV1 IS A TEMPORARY STORAGE ARRAY.                              
C                                                                       
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,        
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY         
C                                                                       
C     ------------------------------------------------------------------
C                                                                       
C     ********** MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING     
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.   
C                                                                       
C                **********                                             
      MACHEP = 2.0D+00**(-50)                                           
C                                                                       
      IERR = 0                                                          
      K = 0                                                             
      TAG = 0                                                           
C                                                                       
      DO 100 I = 1, N                                                   
      W(I) = D(I)                                                       
      IF (I .NE. 1) RV1(I-1) = E(I)                                     
  100 CONTINUE                                                          
C                                                                       
      E2(1) = 0.0D+00                                                   
      RV1(N) = 0.0D+00                                                  
C                                                                       
      DO 360 L = 1, N                                                   
      J = 0                                                             
C     ********** LOOK FOR SMALL SUB-DIAGONAL ELEMENT **********         
  120 DO 140 M = L, N                                                   
      IF (M .EQ. N) GO TO 160                                           
      IF (ABS(RV1(M)) .LE. MACHEP * (ABS(W(M)) + ABS(W(M+1)))) GO TO    
     +     160                                                          
C     ********** GUARD AGAINST UNDERFLOWED ELEMENT OF E2 **********     
      IF (E2(M+1) .EQ. 0.0D+00) GO TO 180                               
  140 CONTINUE                                                          
C                                                                       
  160 IF (M .LE. K) GO TO 200                                           
      IF (M .NE. N) E2(M+1) = 0.0D+00                                   
  180 K = M                                                             
      TAG = TAG + 1                                                     
  200 P = W(L)                                                          
      IF (M .EQ. L) GO TO 280                                           
      IF (J .EQ. 30) GO TO 380                                          
      J = J + 1                                                         
C     ********** FORM SHIFT **********                                  
      G = (W(L+1) - P) / (2.0D+00 * RV1(L))                             
      R = SQRT(G*G+1.0D+00)                                             
      G = W(M) - P + RV1(L) / (G + SIGN(R,G))                           
      S = 1.0D+00                                                       
      C = 1.0D+00                                                       
      P = 0.0D+00                                                       
      MML = M - L                                                       
C     ********** FOR I=M-1 STEP -1 UNTIL L DO -- **********             
      DO 260 II = 1, MML                                                
      I = M - II                                                        
      F = S * RV1(I)                                                    
      B = C * RV1(I)                                                    
      IF (ABS(F) .LT. ABS(G)) GO TO 220                                 
      C = G / F                                                         
      R = SQRT(C*C+1.0D+00)                                             
      RV1(I+1) = F * R                                                  
      S = 1.0D+00 / R                                                   
      C = C * S                                                         
      GO TO 240                                                         
  220 S = F / G                                                         
      R = SQRT(S*S+1.0D+00)                                             
      RV1(I+1) = G * R                                                  
      C = 1.0D+00 / R                                                   
      S = S * C                                                         
  240 G = W(I+1) - P                                                    
      R = (W(I) - G) * S + 2.0D+00 * C * B                              
      P = S * R                                                         
      W(I+1) = G + P                                                    
      G = C * R - B                                                     
  260 CONTINUE                                                          
C                                                                       
      W(L) = W(L) - P                                                   
      RV1(L) = G                                                        
      RV1(M) = 0.0D+00                                                  
      GO TO 120                                                         
C     ********** ORDER EIGENVALUES **********                           
  280 IF (L .EQ. 1) GO TO 320                                           
C     ********** FOR I=L STEP -1 UNTIL 2 DO -- **********               
      DO 300 II = 2, L                                                  
      I = L + 2 - II                                                    
      IF (P .GE. W(I-1)) GO TO 340                                      
      W(I) = W(I-1)                                                     
      IND(I) = IND(I-1)                                                 
  300 CONTINUE                                                          
C                                                                       
  320 I = 1                                                             
  340 W(I) = P                                                          
      IND(I) = TAG                                                      
  360 CONTINUE                                                          
C                                                                       
      GO TO 400                                                         
C     ********** SET ERROR -- NO CONVERGENCE TO AN                      
C                EIGENVALUE AFTER 30 ITERATIONS **********              
  380 IERR = L                                                          
  400 RETURN                                                            
C     ********** LAST CARD OF IMTQLV **********                         
      END                                                               
C                                                                       
C*MODULE EIGEN   *DECK TINVTB                                           
      SUBROUTINE TINVTB(NM,N,D,E,E2,M,W,IND,Z,                          
     *                  IERR,RV1,RV2,RV3,RV4,RV6)                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      DIMENSION D(N),E(N),E2(N),W(M),Z(NM,M),                           
     *          RV1(N),RV2(N),RV3(N),RV4(N),RV6(N),IND(M)               
      DOUBLE PRECISION MACHEP,NORM                                      
      INTEGER P,Q,R,S,TAG,GROUP                                         
C     ------------------------------------------------------------------
C                                                                       
C     THIS ROUTINE IS A TRANSLATION OF THE INVERSE ITERATION TECH-      
C     NIQUE IN THE ALGOL PROCEDURE TRISTURM BY PETERS AND WILKINSON.    
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 418-439(1971).   
C                                                                       
C     THIS ROUTINE FINDS THOSE EIGENVECTORS OF A TRIDIAGONAL            
C     SYMMETRIC MATRIX CORRESPONDING TO SPECIFIED EIGENVALUES,          
C     USING INVERSE ITERATION.                                          
C                                                                       
C     ON INPUT-                                                         
C                                                                       
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL         
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING ROUTINE          
C          DIMENSION STATEMENT,                                         
C                                                                       
C        N IS THE ORDER OF THE MATRIX,                                  
C                                                                       
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX,          
C                                                                       
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX        
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY,               
C                                                                       
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E,    
C          WITH ZEROS CORRESPONDING TO NEGLIGIBLE ELEMENTS OF E.        
C          E(I) IS CONSIDERED NEGLIGIBLE IF IT IS NOT LARGER THAN       
C          THE PRODUCT OF THE RELATIVE MACHINE PRECISION AND THE SUM    
C          OF THE MAGNITUDES OF D(I) AND D(I-1).  E2(1) MUST CONTAIN    
C          0.0 IF THE EIGENVALUES ARE IN ASCENDING ORDER, OR 2.0        
C          IF THE EIGENVALUES ARE IN DESCENDING ORDER.  IF  BISECT,     
C          TRIDIB, OR  IMTQLV  HAS BEEN USED TO FIND THE EIGENVALUES,   
C          THEIR OUTPUT E2 ARRAY IS EXACTLY WHAT IS EXPECTED HERE,      
C                                                                       
C        M IS THE NUMBER OF SPECIFIED EIGENVALUES,                      
C                                                                       
C        W CONTAINS THE M EIGENVALUES IN ASCENDING OR DESCENDING ORDER, 
C                                                                       
C        IND CONTAINS IN ITS FIRST M POSITIONS THE SUBMATRIX INDICES    
C          ASSOCIATED WITH THE CORRESPONDING EIGENVALUES IN W --        
C          1 FOR EIGENVALUES BELONGING TO THE FIRST SUBMATRIX FROM      
C          THE TOP, 2 FOR THOSE BELONGING TO THE SECOND SUBMATRIX, ETC. 
C                                                                       
C     ON OUTPUT-                                                        
C                                                                       
C        ALL INPUT ARRAYS ARE UNALTERED,                                
C                                                                       
C        Z CONTAINS THE ASSOCIATED SET OF ORTHONORMAL EIGENVECTORS.     
C          ANY VECTOR WHICH FAILS TO CONVERGE IS SET TO ZERO,           
C                                                                       
C        IERR IS SET TO                                                 
C          ZERO       FOR NORMAL RETURN,                                
C          -R         IF THE EIGENVECTOR CORRESPONDING TO THE R-TH      
C                     EIGENVALUE FAILS TO CONVERGE IN 5 ITERATIONS,     
C                                                                       
C        RV1, RV2, RV3, RV4, AND RV6 ARE TEMPORARY STORAGE ARRAYS.      
C                                                                       
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,        
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY         
C                                                                       
C     ------------------------------------------------------------------
C                                                                       
C     ********** MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING     
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.   
C                                                                       
C                **********                                             
      MACHEP = 2.0D+00**(-50)                                           
C                                                                       
      IERR = 0                                                          
      IF (M .EQ. 0) GO TO 680                                           
      TAG = 0                                                           
      ORDER = 1.0D+00 - E2(1)                                           
      XU = 0.0D+00                                                      
      UK = 0.0D+00                                                      
      X0 = 0.0D+00                                                      
      U  = 0.0D+00                                                      
      EPS2 = 0.0D+00                                                    
      EPS3 = 0.0D+00                                                    
      EPS4 = 0.0D+00                                                    
      GROUP = 0                                                         
      Q = 0                                                             
C     ********** ESTABLISH AND PROCESS NEXT SUBMATRIX **********        
  100 P = Q + 1                                                         
      IP = P + 1                                                        
C                                                                       
      DO 120 Q = P, N                                                   
      IF (Q .EQ. N) GO TO 140                                           
      IF (E2(Q+1) .EQ. 0.0D+00) GO TO 140                               
  120 CONTINUE                                                          
C     ********** FIND VECTORS BY INVERSE ITERATION **********           
  140 TAG = TAG + 1                                                     
      IQMP = Q - P + 1                                                  
      S = 0                                                             
C                                                                       
      DO 660 R = 1, M                                                   
      IF (IND(R) .NE. TAG) GO TO 660                                    
      ITS = 1                                                           
      X1 = W(R)                                                         
      IF (S .NE. 0) GO TO 220                                           
C     ********** CHECK FOR ISOLATED ROOT **********                     
      XU = 1.0D+00                                                      
      IF (P .NE. Q) GO TO 160                                           
      RV6(P) = 1.0D+00                                                  
      GO TO 600                                                         
  160 NORM = ABS(D(P))                                                  
C                                                                       
      DO 180 I = IP, Q                                                  
  180 NORM = NORM + ABS(D(I)) + ABS(E(I))                               
C     ********** EPS2 IS THE CRITERION FOR GROUPING,                    
C                EPS3 REPLACES ZERO PIVOTS AND EQUAL                    
C                ROOTS ARE MODIFIED BY EPS3,                            
C                EPS4 IS TAKEN VERY SMALL TO AVOID OVERFLOW **********  
      EPS2 = 1.0D-03 * NORM                                             
      EPS3 = MACHEP * NORM                                              
      UK = IQMP                                                         
      EPS4 = UK * EPS3                                                  
      UK = EPS4 / SQRT(UK)                                              
      S = P                                                             
  200 GROUP = 0                                                         
      GO TO 240                                                         
C     ********** LOOK FOR CLOSE OR COINCIDENT ROOTS **********          
  220 IF (ABS(X1-X0) .GE. EPS2) GO TO 200                               
      GROUP = GROUP + 1                                                 
      IF (ORDER * (X1 - X0) .LE. 0.0D+00) X1 = X0 + ORDER * EPS3        
C     ********** ELIMINATION WITH INTERCHANGES AND                      
C                INITIALIZATION OF VECTOR **********                    
  240 V = 0.0D+00                                                       
C                                                                       
      DO 300 I = P, Q                                                   
      RV6(I) = UK                                                       
      IF (I .EQ. P) GO TO 280                                           
      IF (ABS(E(I)) .LT. ABS(U)) GO TO 260                              
C     ********** WARNING -- A DIVIDE CHECK MAY OCCUR HERE IF            
C                E2 ARRAY HAS NOT BEEN SPECIFIED CORRECTLY **********   
      XU = U / E(I)                                                     
      RV4(I) = XU                                                       
      RV1(I-1) = E(I)                                                   
      RV2(I-1) = D(I) - X1                                              
      RV3(I-1) = 0.0D+00                                                
      IF (I .NE. Q) RV3(I-1) = E(I+1)                                   
      U = V - XU * RV2(I-1)                                             
      V = -XU * RV3(I-1)                                                
      GO TO 300                                                         
  260 XU = E(I) / U                                                     
      RV4(I) = XU                                                       
      RV1(I-1) = U                                                      
      RV2(I-1) = V                                                      
      RV3(I-1) = 0.0D+00                                                
  280 U = D(I) - X1 - XU * V                                            
      IF (I .NE. Q) V = E(I+1)                                          
  300 CONTINUE                                                          
C                                                                       
      IF (U .EQ. 0.0D+00) U = EPS3                                      
      RV1(Q) = U                                                        
      RV2(Q) = 0.0D+00                                                  
      RV3(Q) = 0.0D+00                                                  
C     ********** BACK SUBSTITUTION                                      
C                FOR I=Q STEP -1 UNTIL P DO -- **********               
  320 DO 340 II = P, Q                                                  
      I = P + Q - II                                                    
      RV6(I) = (RV6(I) - U * RV2(I) - V * RV3(I)) / RV1(I)              
      V = U                                                             
      U = RV6(I)                                                        
  340 CONTINUE                                                          
C     ********** ORTHOGONALIZE WITH RESPECT TO PREVIOUS                 
C                MEMBERS OF GROUP **********                            
      IF (GROUP .EQ. 0) GO TO 400                                       
      J = R                                                             
C                                                                       
      DO 380 JJ = 1, GROUP                                              
  360 J = J - 1                                                         
      IF (IND(J) .NE. TAG) GO TO 360                                    
      XU = XDOT(IQMP,RV6(P),1,Z(P,J),1)                                 
C                                                                       
      CALL XAXPY(IQMP,-XU,Z(P,J),1,RV6(P),1)                            
C                                                                       
  380 CONTINUE                                                          
C                                                                       
  400 NORM = 0.0D+00                                                    
C                                                                       
      DO 420 I = P, Q                                                   
  420 NORM = NORM + ABS(RV6(I))                                         
C                                                                       
      IF (NORM .GE. 1.0D+00) GO TO 560                                  
C     ********** FORWARD SUBSTITUTION **********                        
      IF (ITS .EQ. 5) GO TO 540                                         
      IF (NORM .NE. 0.0D+00) GO TO 440                                  
      RV6(S) = EPS4                                                     
      S = S + 1                                                         
      IF (S .GT. Q) S = P                                               
      GO TO 480                                                         
  440 XU = EPS4 / NORM                                                  
C                                                                       
      DO 460 I = P, Q                                                   
  460 RV6(I) = RV6(I) * XU                                              
C     ********** ELIMINATION OPERATIONS ON NEXT VECTOR                  
C                ITERATE **********                                     
  480 DO 520 I = IP, Q                                                  
      U = RV6(I)                                                        
C     ********** IF RV1(I-1) .EQ. E(I), A ROW INTERCHANGE               
C                WAS PERFORMED EARLIER IN THE                           
C                TRIANGULARIZATION PROCESS **********                   
      IF (RV1(I-1) .NE. E(I)) GO TO 500                                 
      U = RV6(I-1)                                                      
      RV6(I-1) = RV6(I)                                                 
  500 RV6(I) = U - RV4(I) * RV6(I-1)                                    
  520 CONTINUE                                                          
C                                                                       
      ITS = ITS + 1                                                     
      GO TO 320                                                         
C     ********** SET ERROR -- NON-CONVERGED EIGENVECTOR **********      
  540 IERR = -R                                                         
      XU = 0.0D+00                                                      
      GO TO 600                                                         
C     ********** NORMALIZE SO THAT SUM OF SQUARES IS                    
C                1 AND EXPAND TO FULL ORDER **********                  
  560 U = 0.0D+00                                                       
C                                                                       
      DO 580 I = P, Q                                                   
      RV6(I) = RV6(I) / NORM                                            
  580 U = U + RV6(I)**2                                                 
C                                                                       
      XU = 1.0D+00 / SQRT(U)                                            
C                                                                       
  600 DO 620 I = 1, N                                                   
  620 Z(I,R) = 0.0D+00                                                  
C                                                                       
      DO 640 I = P, Q                                                   
  640 Z(I,R) = RV6(I) * XU                                              
C                                                                       
      X0 = X1                                                           
  660 CONTINUE                                                          
C                                                                       
      IF (Q .LT. N) GO TO 100                                           
  680 RETURN                                                            
C     ********** LAST CARD OF TINVIT **********                         
      END                                                               
C                                                                       
C*MODULE EIGEN   *DECK TQL2                                             
      SUBROUTINE TQL2(NM,N,D,E,Z,IERR)                                  
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      DOUBLE PRECISION MACHEP                                           
      DIMENSION D(N),E(N),Z(NM,N)                                       
C                                                                       
C     THIS ROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQL2,        
C     NUM. MATH. 11, 293-306(1968) BY BOWDLER, MARTIN, REINSCH, AND     
C     WILKINSON.                                                        
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).   
C                                                                       
C     THIS ROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS               
C     OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE QL METHOD.               
C     THE EIGENVECTORS OF A FULL SYMMETRIC MATRIX CAN ALSO              
C     BE FOUND IF  TRED2  HAS BEEN USED TO REDUCE THIS                  
C     FULL MATRIX TO TRIDIAGONAL FORM.                                  
C                                                                       
C     ON INPUT-                                                         
C                                                                       
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL         
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING ROUTINE          
C          DIMENSION STATEMENT,                                         
C                                                                       
C        N IS THE ORDER OF THE MATRIX,                                  
C                                                                       
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX,          
C                                                                       
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX        
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY,               
C                                                                       
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE           
C          REDUCTION BY  TRED2, IF PERFORMED.  IF THE EIGENVECTORS      
C          OF THE TRIDIAGONAL MATRIX ARE DESIRED, Z MUST CONTAIN        
C          THE IDENTITY MATRIX.                                         
C                                                                       
C      ON OUTPUT-                                                       
C                                                                       
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN          
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT          
C          UNORDERED FOR INDICES 1,2,...,IERR-1,                        
C                                                                       
C        E HAS BEEN DESTROYED,                                          
C                                                                       
C        Z CONTAINS ORTHONORMAL EIGENVECTORS OF THE SYMMETRIC           
C          TRIDIAGONAL (OR FULL) MATRIX.  IF AN ERROR EXIT IS MADE,     
C          Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED       
C          EIGENVALUES,                                                 
C                                                                       
C        IERR IS SET TO                                                 
C          ZERO       FOR NORMAL RETURN,                                
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN               
C                     DETERMINED AFTER 30 ITERATIONS.                   
C                                                                       
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,        
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY         
C                                                                       
C     ------------------------------------------------------------------
C                                                                       
C     ********** MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING     
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.   
C                                                                       
C                **********                                             
      MACHEP = 2.0D+00**(-50)                                           
C                                                                       
      IERR = 0                                                          
      IF (N .EQ. 1) GO TO 400                                           
C                                                                       
      DO 100 I = 2, N                                                   
  100 E(I-1) = E(I)                                                     
C                                                                       
      F = 0.0D+00                                                       
      B = 0.0D+00                                                       
      E(N) = 0.0D+00                                                    
C                                                                       
      DO 300 L = 1, N                                                   
      J = 0                                                             
      H = MACHEP * (ABS(D(L)) + ABS(E(L)))                              
      IF (B .LT. H) B = H                                               
C     ********** LOOK FOR SMALL SUB-DIAGONAL ELEMENT **********         
      DO 120 M = L, N                                                   
      IF (ABS(E(M)) .LE. B) GO TO 140                                   
C     ********** E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT               
C                THROUGH THE BOTTOM OF THE LOOP **********              
  120 CONTINUE                                                          
C                                                                       
  140 IF (M .EQ. L) GO TO 280                                           
  160 IF (J .EQ. 30) GO TO 380                                          
      J = J + 1                                                         
C     ********** FORM SHIFT **********                                  
      L1 = L + 1                                                        
      G = D(L)                                                          
      P = (D(L1) - G) / (2.0D+00 * E(L))                                
      R = SQRT(P*P+1.0D+00)                                             
      D(L) = E(L) / (P + SIGN(R,P))                                     
      H = G - D(L)                                                      
C                                                                       
      DO 180 I = L1, N                                                  
  180 D(I) = D(I) - H                                                   
C                                                                       
      F = F + H                                                         
C     ********** QL TRANSFORMATION **********                           
      P = D(M)                                                          
      C = 1.0D+00                                                       
      S = 0.0D+00                                                       
      MML = M - L                                                       
C     ********** FOR I=M-1 STEP -1 UNTIL L DO -- **********             
      DO 260 II = 1, MML                                                
      I = M - II                                                        
      G = C * E(I)                                                      
      H = C * P                                                         
      IF (ABS(P) .LT. ABS(E(I))) GO TO 200                              
      C = E(I) / P                                                      
      R = SQRT(C*C+1.0D+00)                                             
      E(I+1) = S * P * R                                                
      S = C / R                                                         
      C = 1.0D+00 / R                                                   
      GO TO 220                                                         
  200 C = P / E(I)                                                      
      R = SQRT(C*C+1.0D+00)                                             
      E(I+1) = S * E(I) * R                                             
      S = 1.0D+00 / R                                                   
      C = C * S                                                         
  220 P = C * D(I) - S * G                                              
      D(I+1) = H + S * (C * G + S * D(I))                               
C     ********** FORM VECTOR **********                                 
      CALL XROT(N,Z(1,I+1),1,Z(1,I),1,C,S)                              
C                                                                       
  260 CONTINUE                                                          
C                                                                       
      E(L) = S * P                                                      
      D(L) = C * P                                                      
      IF (ABS(E(L)) .GT. B) GO TO 160                                   
  280 D(L) = D(L) + F                                                   
  300 CONTINUE                                                          
C     ********** ORDER EIGENVALUES AND EIGENVECTORS **********          
      DO 360 II = 2, N                                                  
      I = II - 1                                                        
      K = I                                                             
      P = D(I)                                                          
C                                                                       
      DO 320 J = II, N                                                  
      IF (D(J) .GE. P) GO TO 320                                        
      K = J                                                             
      P = D(J)                                                          
  320 CONTINUE                                                          
C                                                                       
      IF (K .EQ. I) GO TO 360                                           
      D(K) = D(I)                                                       
      D(I) = P                                                          
C                                                                       
      CALL XSWAP(N,Z(1,I),1,Z(1,K),1)                                   
C                                                                       
  360 CONTINUE                                                          
C                                                                       
      GO TO 400                                                         
C     ********** SET ERROR -- NO CONVERGENCE TO AN                      
C                EIGENVALUE AFTER 30 ITERATIONS **********              
  380 IERR = L                                                          
  400 RETURN                                                            
C     ********** LAST CARD OF TQL2 **********                           
      END                                                               
C                                                                       
C*MODULE EIGEN   *DECK TRBK3B                                           
      SUBROUTINE TRBK3B(NM,N,NV,A,M,Z)                                  
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      DIMENSION A(NV),Z(NM,M)                                           
C                                                                       
C     THIS ROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRBAK3,      
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.   
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).   
C                                                                       
C     THIS ROUTINE FORMS THE EIGENVECTORS OF A REAL SYMMETRIC           
C     MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING            
C     SYMMETRIC TRIDIAGONAL MATRIX DETERMINED BY  TRED3B.               
C                                                                       
C     ON INPUT-                                                         
C                                                                       
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL         
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING ROUTINE          
C          DIMENSION STATEMENT,                                         
C                                                                       
C        N IS THE ORDER OF THE MATRIX,                                  
C                                                                       
C        NV MUST BE SET TO THE DIMENSION OF THE ARRAY PARAMETER A       
C          AS DECLARED IN THE CALLING ROUTINE DIMENSION STATEMENT,      
C                                                                       
C        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL TRANSFORMATIONS    
C          USED IN THE REDUCTION BY  TRED3B IN ITS FIRST                
C          N*(N+1)/2 POSITIONS,                                         
C                                                                       
C        M IS THE NUMBER OF EIGENVECTORS TO BE BACK TRANSFORMED,        
C                                                                       
C        Z CONTAINS THE EIGENVECTORS TO BE BACK TRANSFORMED             
C          IN ITS FIRST M COLUMNS.                                      
C                                                                       
C     ON OUTPUT-                                                        
C                                                                       
C        Z CONTAINS THE TRANSFORMED EIGENVECTORS                        
C          IN ITS FIRST M COLUMNS.                                      
C                                                                       
C     NOTE THAT TRBAK3 PRESERVES VECTOR EUCLIDEAN NORMS.                
C                                                                       
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,        
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY         
C                                                                       
C     ------------------------------------------------------------------
C                                                                       
      IF (M .EQ. 0) GO TO 140                                           
      IF (N .EQ. 1) GO TO 140                                           
C                                                                       
      DO 120 I = 2, N                                                   
      L = I - 1                                                         
      IZ = (I * L) / 2                                                  
      IK = IZ + I                                                       
      H = A(IK)                                                         
      IF (H .EQ. 0.0D+00) GO TO 120                                     
C                                                                       
      DO 100 J = 1, M                                                   
      S = -XDOT(L,A(IZ+1),1,Z(1,J),1)                                   
C                                                                       
C     ********** DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW **********   
      S = (S / H) / H                                                   
C                                                                       
      CALL XAXPY(L,S,A(IZ+1),1,Z(1,J),1)                                
C                                                                       
  100 CONTINUE                                                          
C                                                                       
  120 CONTINUE                                                          
C                                                                       
  140 RETURN                                                            
C     ********** LAST CARD OF TRBAK3 **********                         
      END                                                               
C                                                                       
C         END OF -GIVEIS- PACKAGE, AND START OF -JACOBI- PACKAGE        
C                                                                       
C*MODULE EIGEN   *DECK JACDG                                            
      SUBROUTINE JACDG(A,VEC,EIG,JBIG,BIG,LDVEC,N)                      
C                                                                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C                                                                       
      DIMENSION A(*),VEC(LDVEC,N),EIG(N),JBIG(N),BIG(N)                 
C                                                                       
      PARAMETER (ONE=1.0D+00)                                           
C                                                                       
C     ----- JACOBI DIAGONALIZATION OF SYMMETRIC MATRIX -----            
C     SYMMETRIC MATRIX -A- OF DIMENSION -N- IS DESTROYED ON EXIT.       
C     ALL EIGENVECTORS ARE FOUND, SO -VEC- MUST BE SQUARE,              
C     UNLESS SOMEONE TAKES THE TROUBLE TO LOOK AT -NMAX- BELOW.         
C     -BIG- AND -JBIG- ARE SCRATCH WORK ARRAYS.                         
C                                                                       
      CALL VCLR(VEC,1,LDVEC*N)                                          
      DO 20 I = 1,N                                                     
        VEC(I,I) = ONE                                                  
   20 CONTINUE                                                          
C                                                                       
      NB1 = N                                                           
      NB2 = (NB1*NB1+NB1)/2                                             
      NMIN = 1                                                          
      NMAX = NB1                                                        
C                                                                       
      CALL JACDIA(A,VEC,NB1,NB2,LDVEC,NMIN,NMAX,BIG,JBIG)               
C                                                                       
      DO 30 I=1,N                                                       
        EIG(I) = A((I*I+I)/2)                                           
   30 CONTINUE                                                          
C                                                                       
      CALL JACORD(VEC,EIG,NB1,LDVEC)                                    
      RETURN                                                            
      END                                                               
C                                                                       
C*MODULE EIGEN   *DECK JACDIA                                           
      SUBROUTINE JACDIA(F,VEC,NB1,NB2,LDVEC,NMIN,NMAX,BIG,JBIG)         
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      LOGICAL GOPARR,DSKWRK,MASWRK                                      
      DIMENSION F(NB2),VEC(LDVEC,NB1),BIG(NB1),JBIG(NB1)                
C                                                                       
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
C                                                                       
      PARAMETER (ROOT2=0.707106781186548D+00 )                          
      PARAMETER (ZERO=0.0D+00)                                          
      PARAMETER (ONE=1.0D+00)                                           
      PARAMETER (D1050=1.05D+00)                                        
      PARAMETER (D1500=1.5D+00)                                         
      PARAMETER (D3875=3.875D+00)                                       
      PARAMETER (D0500=0.5D+00)                                         
      PARAMETER (D1375=1.375D+00)                                       
      PARAMETER (D0250=0.25D+00 )                                       
      PARAMETER (C2=1.0D-12)                                            
      PARAMETER (C3=4.0D-16)                                            
      PARAMETER (C4=2.0D-16)                                            
      PARAMETER (C5=8.0D-09)                                            
      PARAMETER (C6=3.0D-06 )                                           
C                                                                       
C      F IS THE MATRIX TO BE DIAGONALIZED, F IS STORED TRIANGULAR       
C      VEC IS THE ARRAY OF EIGENVECTORS, DIMENSION NB1*NB1              
C      BIG AND JBIG ARE TEMPORARY SCRATCH AREAS OF DIMENSION NB1        
C      THE ROTATIONS AMONG THE FIRST NMIN BASIS FUNCTIONS ARE NOT       
C      ACCOUNTED FOR.                                                   
C      THE ROTATIONS AMONG THE LAST NB1-NMAX BASIS FUNCTIONS ARE NOT    
C      ACCOUNTED FOR.                                                   
C                                                                       
      IEAA=0                                                            
      IEAB=0                                                            
      TT=ZERO                                                           
      EPS = 64.0D+00*DEPSLON(ONE)                                       
C                                                                       
C      LOOP OVER COLUMNS (K) OF TRIANGULAR MATRIX TO DETERMINE          
C      LARGEST OFF-DIAGONAL ELEMENTS IN ROW(I).                         
C                                                                       
      DO 20 I=1,NB1                                                     
         BIG(I)=ZERO                                                    
         JBIG(I)=0                                                      
         IF(I.LT.NMIN  .OR.  I.EQ.1) GO TO 20                           
         II = (I*I-I)/2                                                 
         J=MIN(I-1,NMAX)                                                
         DO 10 K=1,J                                                    
            IF(ABS(BIG(I)).GE.ABS(F(II+K))) GO TO 10                    
            BIG(I)=F(II+K)                                              
            JBIG(I)=K                                                   
   10    CONTINUE                                                       
   20 CONTINUE                                                          
C                                                                       
C     ----- 2X2 JACOBI ITERATIONS BEGIN HERE -----                      
C                                                                       
      MAXIT=MAX(NB2*20,500)                                             
      ITER=0                                                            
   30 CONTINUE                                                          
      ITER=ITER+1                                                       
C                                                                       
C      FIND SMALLEST DIAGONAL ELEMENT                                   
C                                                                       
      SD=D1050                                                          
      JJ=0                                                              
      DO 40 J=1,NB1                                                     
         JJ=JJ+J                                                        
         SD= MIN(SD,ABS(F(JJ)))                                         
   40 CONTINUE                                                          
      TEST = MAX(EPS, C2*MAX(SD,C6))                                    
C                                                                       
C      FIND LARGEST OFF-DIAGONAL ELEMENT                                
C                                                                       
      T=ZERO                                                            
      I1=MAX(2,NMIN)                                                    
      IB = I1                                                           
      DO 50 I=I1,NB1                                                    
         IF(T.GE.ABS(BIG(I))) GO TO 50                                  
         T= ABS(BIG(I))                                                 
         IB=I                                                           
   50 CONTINUE                                                          
C                                                                       
C      TEST FOR CONVERGENCE, THEN DETERMINE ROTATION.                   
C                                                                       
      IF(T.LT.TEST) RETURN                                              
C                   ******                                              
C                                                                       
      IF(ITER.GT.MAXIT) THEN                                            
         IF (MASWRK) THEN                                               
            WRITE(6,*) 'JACOBI DIAGONALIZATION FAILS, DIMENSION=',NB1   
            WRITE(6,9020) ITER,T,TEST,SD                                
         ENDIF                                                          
         CALL ABRT                                                      
         STOP                                                           
      END IF                                                            
C                                                                       
      IA=JBIG(IB)                                                       
      IAA=IA*(IA-1)/2                                                   
      IBB=IB*(IB-1)/2                                                   
      DIF=F(IAA+IA)-F(IBB+IB)                                           
      IF(ABS(DIF).GT.C3*T) GO TO 70                                     
      SX=ROOT2                                                          
      CX=ROOT2                                                          
      GO TO 110                                                         
   70 T2X2=BIG(IB)/DIF                                                  
      T2X25=T2X2*T2X2                                                   
      IF(T2X25 . GT . C4) GO TO 80                                      
      CX=ONE                                                            
      SX=T2X2                                                           
      GO TO 110                                                         
   80 IF(T2X25 . GT . C5) GO TO 90                                      
      SX=T2X2*(ONE-D1500*T2X25)                                         
      CX=ONE-D0500*T2X25                                                
      GO TO 110                                                         
   90 IF(T2X25 . GT . C6) GO TO 100                                     
      CX=ONE+T2X25*(T2X25*D1375 - D0500)                                
      SX= T2X2*(ONE + T2X25*(T2X25*D3875 - D1500))                      
      GO TO 110                                                         
  100 T=D0250  / SQRT(D0250   + T2X25)                                  
      CX= SQRT(D0500   + T)                                             
      SX= SIGN( SQRT(D0500   - T),T2X2)                                 
  110 IEAR=IAA+1                                                        
      IEBR=IBB+1                                                        
C                                                                       
      DO 230 IR=1,NB1                                                   
         T=F(IEAR)*SX                                                   
         F(IEAR)=F(IEAR)*CX+F(IEBR)*SX                                  
         F(IEBR)=T-F(IEBR)*CX                                           
         IF(IR-IA) 220,120,130                                          
  120    TT=F(IEBR)                                                     
         IEAA=IEAR                                                      
         IEAB=IEBR                                                      
         F(IEBR)=BIG(IB)                                                
         IEAR=IEAR+IR-1                                                 
         IF(JBIG(IR)) 200,220,200                                       
  130    T=F(IEAR)                                                      
         IT=IA                                                          
         IEAR=IEAR+IR-1                                                 
         IF(IR-IB) 180,150,160                                          
  150    F(IEAA)=F(IEAA)*CX+F(IEAB)*SX                                  
         F(IEAB)=TT*CX+F(IEBR)*SX                                       
         F(IEBR)=TT*SX-F(IEBR)*CX                                       
         IEBR=IEBR+IR-1                                                 
         GO TO 200                                                      
  160    IF(  ABS(T) . GE .  ABS(F(IEBR))) GO TO 170                    
         IF(IB.GT.NMAX) GO TO 170                                       
         T=F(IEBR)                                                      
         IT=IB                                                          
  170    IEBR=IEBR+IR-1                                                 
  180    IF(  ABS(T) . LT .  ABS(BIG(IR))) GO TO 190                    
         BIG(IR) = T                                                    
         JBIG(IR) = IT                                                  
         GO TO 220                                                      
  190    IF(IA . NE . JBIG(IR) . AND . IB . NE . JBIG(IR))  GO TO 220   
  200    KQ=IEAR-IR-IA+1                                                
         BIG(IR)=ZERO                                                   
         IR1=MIN(IR-1,NMAX)                                             
         DO 210 I=1,IR1                                                 
            K=KQ+I                                                      
            IF(ABS(BIG(IR)) . GE . ABS(F(K)))  GO TO 210                
            BIG(IR) = F(K)                                              
            JBIG(IR)=I                                                  
  210    CONTINUE                                                       
  220    IEAR=IEAR+1                                                    
  230    IEBR=IEBR+1                                                    
C                                                                       
      DO 240 I=1,NB1                                                    
         T1=VEC(I,IA)*CX + VEC(I,IB)*SX                                 
         T2=VEC(I,IA)*SX - VEC(I,IB)*CX                                 
         VEC(I,IA)=T1                                                   
         VEC(I,IB)=T2                                                   
  240 CONTINUE                                                          
      GO TO 30                                                          
C                                                                       
 9020 FORMAT(1X,'ITER=',I6,' T,TEST,SD=',1P,3E20.10)                    
      END                                                               
C                                                                       
C*MODULE EIGEN   *DECK JACORD                                           
      SUBROUTINE JACORD(VEC,EIG,N,LDVEC)                                
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      DIMENSION VEC(LDVEC,N),EIG(N)                                     
C                                                                       
C     ---- SORT EIGENDATA INTO ASCENDING ORDER -----                    
C                                                                       
      DO 290 I = 1, N                                                   
         JJ = I                                                         
         DO 270 J = I, N                                                
            IF (EIG(J) .LT. EIG(JJ)) JJ = J                             
  270    CONTINUE                                                       
         IF (JJ .EQ. I) GO TO 290                                       
         T = EIG(JJ)                                                    
         EIG(JJ) = EIG(I)                                               
         EIG(I) = T                                                     
         DO 280 J = 1, N                                                
            T = VEC(J,JJ)                                               
            VEC(J,JJ) = VEC(J,I)                                        
            VEC(J,I) = T                                                
  280    CONTINUE                                                       
  290 CONTINUE                                                          
      RETURN                                                            
      END                                                               
c                                                                       
C*MODULE EIGEN   *DECK matpar                                           
      SUBROUTINE matpar(n,mode,n2,n3)                                   
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      COMMON /MTXSIZ/ MXSEQ2,MXSEQ3                                     
c                                                                       
      if(mode.eq.0) then                                                
        n2=MXSEQ2                                                       
        MXSEQ2=2**30                                                    
        if(n.gt.2) n3=MXSEQ3                                            
        if(n.gt.2) MXSEQ3=2**30                                         
      else                                                              
        MXSEQ2=n2                                                       
        if(n.gt.2) MXSEQ3=n3                                            
      endif                                                             
      RETURN                                                            
      END                                                               
