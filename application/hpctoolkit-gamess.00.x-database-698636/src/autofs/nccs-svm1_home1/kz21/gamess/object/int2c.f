C  6 Jun 18 - DGF - tweaks for FMO 5.3                                  
C 18 Apr 16 - DGF - pad common blocks                                   
C 22 Nov 13 - JMS - simplify machine generated formulae                 
C 23 MAR 12 - DGF - prevent slaves from opening and closing ERIC in GDDI
C 30 JAN 12 - SAN - RUNTYP=COMP RESET OF ERIC, ADDED /CXTHRM/           
C 25 MAR 10 - GDF - DIMENSION WITH STANDARD PARAMETER -MXGSH-           
C 20 AUG 07 - MWS - ALLOW H,I INDEXING RANGE                            
C 22 DEC 06 - ST,NK,MC - ADD LC EXCHANGE ARGUMENTS                      
C 19 SEP 05 - GDF - ERIPRE: WEE CHANGE FOR H+I IN RYS CODES             
C  5 JUL 05 - MWS - SELECT NEW ATOM,BASIS,EFP,PCM,DAF DIMENSIONS        
C 13 FEB 05 - MWS - PAD COMMON BLOCK NSHEL                              
C  5 FEB 05 - MWS - CHANGE COMMON NAME GAMMAF TO BE FMTTBL INSTEAD      
C 10 NOV 04 - MWS - ENSURE FM(T) BROADCAST ALWAYS DONE, ADJUST CC TOL   
C  7 SEP 04 - GDF - NEW MODULE FOR ERIC, AND SPECIAL S,P ROUTINES       
C                                                                       
C*MODULE INT2C   *DECK ERIC                                             
C>                                                                      
C>    @brief   ERIC: Electron Repulsion Integral Calculator             
C>                                                                      
C>    @details ERIC uses recursion relations for "precursor Hermite     
C>             functions" to evaluate ERI.  The method is described in  
C>             "Recursion Formula for Electron Repulsion Integrals      
C>             Over Hermite Polynomials"                                
C>             G.D.Fletcher  Int.J.Quantum Chem. 106, 355-360(2006)     
C>                                                                      
C>    @author  Graham Fletcher, 2004, using machine generated formulae. 
C>             Jose Sierra, 2013, uses human-created simplifications.   
C>                                                                      
      SUBROUTINE ERIC (ISH,JSH,KSH,LSH,ERI)                             
      use mx_limits, only: mxatm,mxgtot,mxsh,mxgsh,mxg2                 
C ----------------------------------------------------------------------
C                                                                       
C              ELECTRON REPULSION INTEGRAL CALCULATOR                   
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
C  METHOD:                                                              
C     A) FORM SCALED, CONTRACTED 1-CENTER PRECURSOR INTEGRALS.          
C        CONVERT THESE TO 4-CENTER INTEGRALS OVER CARTESIAN             
C        GAUSSIANS USING,                                               
C     B) PRECURSOR-HERMITE TRANSFER EQUATION (PTE)                      
C     C) CONTRACTED TRANSFER  EQUATION (CTE)                            
C     D) HORIZONTAL RECURSION RELATION (HRR)                            
C                                                                       
C                                                                       
C  REFERENCES                                                           
C                                                                       
C     PTE:                                                              
C     "Recursion Formula for Electron Repulsion Integrals Over          
C     Hermite Polynomials"                                              
C         G.D.Fletcher  Int.J.Quantum Chem. 106, 355-360(2006)          
C                                                                       
C     HRR:                                                              
C     M. HEAD-GORDON & J. A. POPLE,                                     
C     J. CHEM. PHYS., 89, 5777-5786 (1988).                             
C                                                                       
C     CTE:                                                              
C     P. M. W. GILL, M. HEAD-GORDON, & J. A. POPLE,                     
C     INT. J. Q. CHEM., SYMP. 23, 269-280 (1989).                       
C                                                                       
C     INTERPOLATION METHOD:                                             
C     P. M. W. GILL, B. G. JOHNSON, & J. A. POPLE,                      
C     INT. J. Q. CHEM., 40, 745-752 (1991).                             
C                                                                       
C                                                                       
C  THIS IS THE DRIVER ROUTINE THAT INTERFACES TO GAMESS COMMON NSHEL    
C                                                                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C                                                                       
      DOUBLE PRECISION ERI(*)                                           
C                                                                       
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,                    
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)                  
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),      
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),      
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),   
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL           
C                                                                       
C  /ERIPRM/ SYMBOLS:                                                    
C     EXI,J,K,L   = PRIMITIVE EXPONENTS OF I,J,K,L SHELLS               
C     CCI,J,K,L   = CONTRACTION COEFFICIENTS OF ISH,JSH,KSH,LSH         
C     XAB,YAB,ZAB = EXPONENT-WEIGHTED DISTANCE FOR BRA                  
C     XCD,YCD,ZCD = EXPONENT-WEIGHTED DISTANCE FOR KET                  
C     CCBRA,CCKET = PRIMITIVE CHARGE-CLOUD FACTOR FOR BRA,KET           
C     SLBRA,SLKET = FACTOR CONVERTING SMALL-T EXPRESSION TO LARGE-T     
C                                                                       
      COMMON /ERIPRM/ EXI(MXGSH),EXJ(MXGSH),EXK(MXGSH),EXL(MXGSH),      
     *                CCI(MXGSH),CCJ(MXGSH),CCK(MXGSH),CCL(MXGSH),      
     *                XAB(MXG2),YAB(MXG2),ZAB(MXG2),                    
     *                XCD(MXG2),YCD(MXG2),ZCD(MXG2),                    
     *                CCBRA(MXG2),CCKET(MXG2),RXB(MXG2),                
     *                SLBRA(MXG2),SLKET(MXG2),RXK(MXG2)                 
C$omp threadprivate(/ERIPRM/)
C                                                                       
      COMMON /ERIDAT/ LEN1,LEN2,LEN3,LEN4                               
      COMMON /ERIOUT/ INW,JNW,KNW,LNW,LSTRI,LSTRJ,LSTRK,LSTRL           
C$omp threadprivate(/ERIOUT/)
C                                                                       
      COMMON /FMCOM /XX(1)                                              
C                                                                       
      LOGICAL IEQJ,KEQL                                                 
C                                                                       
      PARAMETER (PT5=0.5D+00)                                           
      PARAMETER (ONE=1.0D+00)                                           
      PARAMETER (TWO=2.0D+00)                                           
C                                                                       
C     SR3 <==> SQRT( 3)   <==> 1.7320508075688773D+00                   
C     SR5 <==> SQRT( 5)   <==> 2.2360679774997897D+00                   
C     SR7 <==> SQRT( 7)   <==> 2.6457513110645906D+00                   
C     S15 <==> SQRT(15)   <==> 3.8729833462074169D+00                   
C     S35 <==> SQRT(35)   <==> 5.9160797830996160D+00                   
C     S53 <==> SQRT(35/3) <==> 3.4156502553198661D+00                   
C                                                                       
CC    PARAMETER (SR3=1.7320508075688772D+00,                            
CC   *           SR5=2.2360679774997902D+00,                            
CC   *           SR7=2.6457513110645903D+00,                            
CC   *           S15=3.8729833462074170D+00,                            
CC   *           S35=5.9160797830996170D+00,                            
CC   *           S53=3.4156502553198673D+00)                            
      PARAMETER (SR3=1.7320508075688773D+00)                            
      PARAMETER (SR5=2.2360679774997897D+00)                            
      PARAMETER (SR7=2.6457513110645906D+00)                            
      PARAMETER (S15=SR3*SR5)                                           
      PARAMETER (S35=SR5*SR7)                                           
      PARAMETER (S53=S35/SR3)                                           
C                                                                       
CJMS  PARAMETER (PI =3.14159265358979323846D+00)                        
CC    PARAMETER (PI =3.1415926535897932D+00)                            
C                                                                       
CC               PI254 = PI * SQRT( 2 * SQRT(PI ))                      
CC               PI214 = PI254 /(PI +PI )                               
CC    PARAMETER (PI254 =5.914967172795613D+00)                          
CC    PARAMETER (PI214 =0.9413962637767148D+00)                         
      PARAMETER (PI254 =5.9149671727956129D+00)                         
      PARAMETER (PI214 =0.94139626377671481D+00)                        
C                                                                       
      DIMENSION IORD(35),ANGL(35)                                       
      DATA IORD/                                                        
     1       1,                                                         
     2       2,  3,  4,                                                 
     3       5,  7, 10,  6,  8,  9,                                     
     4      11, 14, 20, 12, 15, 13, 17, 18, 19, 16,                     
     5      21, 25, 35, 22, 26, 24, 29, 33, 34, 23, 30, 32, 27, 28, 31/ 
      DATA ANGL/                                                        
     1     ONE,                                                         
     2     ONE,ONE,ONE,                                                 
     3     ONE,SR3,ONE,SR3,SR3,ONE,                                     
     4     ONE,SR5,SR5,ONE,SR5,S15,SR5,SR5,SR5,ONE,                     
     5     ONE,SR7,S53,SR7,ONE,SR7,S35,S35,SR7,S53,S35,S53,SR7,SR7,ONE/ 
C                                                                       
CC    DATA ANGL/1.0D+00,  1.0D+00, 1.0D+00, 1.0D+00,                    
CC   * 1.0D+00, 1.7320508075688772D+00, 1.0D+00, 1.7320508075688772D+00,
CC   * 1.7320508075688772D+00, 1.0D+00,                                 
CC   * 1.0D+00, 2.2360679774997902D+00, 2.2360679774997902D+00, 1.0D+00,
CC   * 2.2360679774997902D+00, 3.8729833462074170D+00,                  
CC   * 2.2360679774997902D+00,                                          
CC   * 2.2360679774997902D+00, 2.2360679774997902D+00, 1.0D+00,         
CC   * 1.0D+00, 2.6457513110645903D+00, 3.4156502553198673D+00,         
CC   * 2.6457513110645903D+00, 1.0D+00, 2.6457513110645903D+00,         
CC   * 5.9160797830996170D+00, 5.9160797830996170D+00,                  
CC   * 2.6457513110645903D+00,                                          
CC   * 3.4156502553198673D+00, 5.9160797830996170D+00,                  
CC   * 3.4156502553198673D+00,                                          
CC   * 2.6457513110645903D+00, 2.6457513110645903D+00, 1.0D+00/         
C                                                                       
C  FAST CODES                                                           
C                                                                       
C  RE-ORDER CHARGE-CLOUD SHELLS                                         
C  (HRR DATA ONLY STORED FOR LI>=LJ)                                    
C                                                                       
      INW = ISH                                                         
      JNW = JSH                                                         
      KNW = KSH                                                         
      LNW = LSH                                                         
      LSTRI = LEN4                                                      
      LSTRJ = LEN3                                                      
      LSTRK = LEN2                                                      
      LSTRL = LEN1                                                      
      LANGI = KTYPE(INW) - 1                                            
      LANGJ = KTYPE(JNW) - 1                                            
      LANGK = KTYPE(KNW) - 1                                            
      LANGL = KTYPE(LNW) - 1                                            
      IF(LANGI.LT.LANGJ) THEN                                           
         INW   = JSH                                                    
         JNW   = ISH                                                    
         ITMP  = LSTRI                                                  
         LSTRI = LSTRJ                                                  
         LSTRJ = ITMP                                                   
      END IF                                                            
      IF(LANGK.LT.LANGL) THEN                                           
         KNW   = LSH                                                    
         LNW   = KSH                                                    
         ITMP  = LSTRK                                                  
         LSTRK = LSTRL                                                  
         LSTRL = ITMP                                                   
      END IF                                                            
C                                                                       
C  SWAP CHARGE CLOUDS FOR EFFICIENCY...                                 
C                                                                       
      LBRA  = LANGI + LANGJ                                             
      LKET  = LANGK + LANGL                                             
      IF(LBRA.GT.LKET) THEN                                             
         ITMP  = INW                                                    
         INW   = KNW                                                    
         KNW   = ITMP                                                   
         ITMP  = LSTRI                                                  
         LSTRI = LSTRK                                                  
         LSTRK = ITMP                                                   
         ITMP  = JNW                                                    
         JNW   = LNW                                                    
         LNW   = ITMP                                                   
         ITMP  = LSTRJ                                                  
         LSTRJ = LSTRL                                                  
         LSTRL = ITMP                                                   
C                                                                       
C  ...OR TO SELECT THE FAST CODE                                        
C                                                                       
      ELSE IF(LBRA.EQ.LKET) THEN                                        
         LANGI = KTYPE(INW)                                             
         LANGJ = KTYPE(JNW)                                             
         LANGK = KTYPE(KNW)                                             
         LANGL = KTYPE(LNW)                                             
         LNGIJ = (LANGI*LANGI-LANGI)/2 + LANGJ                          
         LNGKL = (LANGK*LANGK-LANGK)/2 + LANGL                          
         IF(LNGKL.GT.LNGIJ) THEN                                        
            ITMP  = INW                                                 
            INW   = KNW                                                 
            KNW   = ITMP                                                
            ITMP  = LSTRI                                               
            LSTRI = LSTRK                                               
            LSTRK = ITMP                                                
            ITMP  = JNW                                                 
            JNW   = LNW                                                 
            LNW   = ITMP                                                
            ITMP  = LSTRJ                                               
            LSTRJ = LSTRL                                               
            LSTRL = ITMP                                                
         END IF                                                         
      END IF                                                            
C                                                                       
C  MAIN PARAMETERS                                                      
C                                                                       
      LANGI = KTYPE(INW) - 1                                            
      LANGJ = KTYPE(JNW) - 1                                            
      LANGK = KTYPE(KNW) - 1                                            
      LANGL = KTYPE(LNW) - 1                                            
      LBRA  = LANGI + LANGJ                                             
      LKET  = LANGK + LANGL                                             
      IPRIM = KNG(INW)                                                  
      JPRIM = KNG(JNW)                                                  
      KPRIM = KNG(KNW)                                                  
      LPRIM = KNG(LNW)                                                  
C                                                                       
      IEQJ  = INW.EQ.JNW                                                
      KEQL  = KNW.EQ.LNW                                                
C                                                                       
C  SET UP 1-INDEX PARAMETERS                                            
C                                                                       
      K1 = KSTART(INW)                                                  
      K2 = K1 + IPRIM - 1                                               
      II = 0                                                            
      DO I = K1, K2                                                     
         II = II + 1                                                    
         EXI(II) = EX(I)                                                
      END DO                                                            
      ITYP = KTYPE(INW)                                                 
      IF(ITYP.EQ.1) THEN                                                
         II = 0                                                         
         DO I = K1, K2                                                  
            II = II + 1                                                 
            CCI(II) = CS(I)                                             
         END DO                                                         
      ELSE IF(ITYP.EQ.2) THEN                                           
         II = 0                                                         
         DO I = K1, K2                                                  
            II = II + 1                                                 
            CCI(II) = CP(I)                                             
         END DO                                                         
      ELSE IF(ITYP.EQ.3) THEN                                           
         II = 0                                                         
         DO I = K1, K2                                                  
            II = II + 1                                                 
            CCI(II) = CD(I)                                             
         END DO                                                         
      ELSE IF(ITYP.EQ.4) THEN                                           
         II = 0                                                         
         DO I = K1, K2                                                  
            II = II + 1                                                 
            CCI(II) = CF(I)                                             
         END DO                                                         
      ELSE IF(ITYP.EQ.5) THEN                                           
         II = 0                                                         
         DO I = K1, K2                                                  
            II = II + 1                                                 
            CCI(II) = CG(I)                                             
         END DO                                                         
      END IF                                                            
C                                                                       
      K1 = KSTART(JNW)                                                  
      K2 = K1 + JPRIM - 1                                               
      JJ = 0                                                            
      DO J = K1, K2                                                     
         JJ = JJ + 1                                                    
         EXJ(JJ) = EX(J)                                                
      END DO                                                            
      JTYP = KTYPE(JNW)                                                 
      IF(JTYP.EQ.1) THEN                                                
         JJ = 0                                                         
         DO J = K1, K2                                                  
            JJ = JJ + 1                                                 
            CCJ(JJ) = CS(J)                                             
         END DO                                                         
      ELSE IF(JTYP.EQ.2) THEN                                           
         JJ = 0                                                         
         DO J = K1, K2                                                  
            JJ = JJ + 1                                                 
            CCJ(JJ) = CP(J)                                             
         END DO                                                         
      ELSE IF(JTYP.EQ.3) THEN                                           
         JJ = 0                                                         
         DO J = K1, K2                                                  
            JJ = JJ + 1                                                 
            CCJ(JJ) = CD(J)                                             
         END DO                                                         
      ELSE IF(JTYP.EQ.4) THEN                                           
         JJ = 0                                                         
         DO J = K1, K2                                                  
            JJ = JJ + 1                                                 
            CCJ(JJ) = CF(J)                                             
         END DO                                                         
      ELSE IF(JTYP.EQ.5) THEN                                           
         JJ = 0                                                         
         DO J = K1, K2                                                  
            JJ = JJ + 1                                                 
            CCJ(JJ) = CG(J)                                             
         END DO                                                         
      END IF                                                            
C                                                                       
      K1 = KSTART(KNW)                                                  
      K2 = K1 + KPRIM - 1                                               
      KK = 0                                                            
      DO K = K1, K2                                                     
         KK = KK + 1                                                    
         EXK(KK) = EX(K)                                                
      END DO                                                            
      KTYP = KTYPE(KNW)                                                 
      IF(KTYP.EQ.1) THEN                                                
         KK = 0                                                         
         DO K = K1, K2                                                  
            KK = KK + 1                                                 
            CCK(KK) = CS(K)                                             
         END DO                                                         
      ELSE IF(KTYP.EQ.2) THEN                                           
         KK = 0                                                         
         DO K = K1, K2                                                  
            KK = KK + 1                                                 
            CCK(KK) = CP(K)                                             
         END DO                                                         
      ELSE IF(KTYP.EQ.3) THEN                                           
         KK = 0                                                         
         DO K = K1, K2                                                  
            KK = KK + 1                                                 
            CCK(KK) = CD(K)                                             
         END DO                                                         
      ELSE IF(KTYP.EQ.4) THEN                                           
         KK = 0                                                         
         DO K = K1, K2                                                  
            KK = KK + 1                                                 
            CCK(KK) = CF(K)                                             
         END DO                                                         
      ELSE IF(KTYP.EQ.5) THEN                                           
         KK = 0                                                         
         DO K = K1, K2                                                  
            KK = KK + 1                                                 
            CCK(KK) = CG(K)                                             
         END DO                                                         
      END IF                                                            
C                                                                       
      K1 = KSTART(LNW)                                                  
      K2 = K1 + LPRIM - 1                                               
      LL = 0                                                            
      DO L = K1, K2                                                     
         LL = LL + 1                                                    
         EXL(LL) = EX(L)                                                
      END DO                                                            
      LTYP = KTYPE(LNW)                                                 
      IF(LTYP.EQ.1) THEN                                                
         LL = 0                                                         
         DO L = K1, K2                                                  
            LL = LL + 1                                                 
            CCL(LL) = CS(L)                                             
         END DO                                                         
      ELSE IF(LTYP.EQ.2) THEN                                           
         LL = 0                                                         
         DO L = K1, K2                                                  
            LL = LL + 1                                                 
            CCL(LL) = CP(L)                                             
         END DO                                                         
      ELSE IF(LTYP.EQ.3) THEN                                           
         LL = 0                                                         
         DO L = K1, K2                                                  
            LL = LL + 1                                                 
            CCL(LL) = CD(L)                                             
         END DO                                                         
      ELSE IF(LTYP.EQ.4) THEN                                           
         LL = 0                                                         
         DO L = K1, K2                                                  
            LL = LL + 1                                                 
            CCL(LL) = CF(L)                                             
         END DO                                                         
      ELSE IF(LTYP.EQ.5) THEN                                           
         LL = 0                                                         
         DO L = K1, K2                                                  
            LL = LL + 1                                                 
            CCL(LL) = CG(L)                                             
         END DO                                                         
      END IF                                                            
C                                                                       
C  2-INDEX PARAMETERS                                                   
C                                                                       
      IATM = KATOM(INW)                                                 
      XA   = C(1,IATM)                                                  
      YA   = C(2,IATM)                                                  
      ZA   = C(3,IATM)                                                  
      JATM = KATOM(JNW)                                                 
      XB   = C(1,JATM)                                                  
      YB   = C(2,JATM)                                                  
      ZB   = C(3,JATM)                                                  
      RAB = (XA-XB)**2 + (YA-YB)**2 + (ZA-ZB)**2                        
      IJ = 0                                                            
      DO JJ = 1, JPRIM                                                  
         EJ = EXJ(JJ)                                                   
         CJ = CCJ(JJ)                                                   
         ITOP = IPRIM                                                   
         IF(IEQJ) ITOP = JJ                                             
         DO II = 1, ITOP                                                
            EI = EXI(II)                                                
            CIX = CCI(II)                                               
            EIJ = ONE/(EI+EJ)                                           
            IJ = IJ + 1                                                 
            CCFAC = PI254*CIX*CJ*EIJ*EXP( -EI*EJ*RAB*EIJ )              
            IF(IEQJ .AND. II.NE.JJ) CCFAC = CCFAC*TWO                   
            CCBRA(IJ) = CCFAC                                           
            SLBRA(IJ) = PI214*SQRT(EIJ)                                 
            XAB(IJ) = (EI*XA + EJ*XB)*EIJ                               
            YAB(IJ) = (EI*YA + EJ*YB)*EIJ                               
            ZAB(IJ) = (EI*ZA + EJ*ZB)*EIJ                               
            RXB(IJ) = EIJ*PT5                                           
         END DO                                                         
      END DO                                                            
C                                                                       
      KATM = KATOM(KNW)                                                 
      XC   = C(1,KATM)                                                  
      YC   = C(2,KATM)                                                  
      ZC   = C(3,KATM)                                                  
      LATM = KATOM(LNW)                                                 
      XD   = C(1,LATM)                                                  
      YD   = C(2,LATM)                                                  
      ZD   = C(3,LATM)                                                  
      RCD = (XC-XD)**2 + (YC-YD)**2 + (ZC-ZD)**2                        
      KL = 0                                                            
      DO LL = 1, LPRIM                                                  
         EL = EXL(LL)                                                   
         CL = CCL(LL)                                                   
         KTOP = KPRIM                                                   
         IF(KEQL) KTOP = LL                                             
         DO KK = 1, KTOP                                                
            EK = EXK(KK)                                                
            CK = CCK(KK)                                                
            EKL = ONE/(EK+EL)                                           
            KL = KL + 1                                                 
            CCFAC = PI254*CK*CL*EKL*EXP( -EK*EL*RCD*EKL )               
            IF(KEQL .AND. KK.NE.LL) CCFAC = CCFAC*TWO                   
            CCKET(KL) = CCFAC                                           
            SLKET(KL) = PI214*SQRT(EKL)                                 
            XCD(KL) = (EK*XC + EL*XD)*EKL                               
            YCD(KL) = (EK*YC + EL*YD)*EKL                               
            ZCD(KL) = (EK*ZC + EL*ZD)*EKL                               
            RXK(KL) = EKL*PT5                                           
         END DO                                                         
      END DO                                                            
C                                                                       
      CALL VALFM(LOADFM)                                                
      IPHI = LOADFM + 1                                                 
C                                                                       
C  ANGULAR MOMENTUM 4-INDEX                                             
C                                                                       
      LBGT  = MAX(LANGI,LANGJ)                                          
      LBLT  = MIN(LANGI,LANGJ)                                          
      LNGIJ = (LBGT*LBGT+LBGT)/2 + LBLT                                 
      LKGT  = MAX(LANGK,LANGL)                                          
      LKLT  = MIN(LANGK,LANGL)                                          
      LNGKL = (LKGT*LKGT+LKGT)/2 + LKLT                                 
      LQGT  = MAX(LNGKL,LNGIJ)                                          
      LQLT  = MIN(LNGKL,LNGIJ)                                          
      LIJKL = (LQGT*LQGT+LQGT)/2 + LQLT                                 
C                                                                       
      IF(LIJKL.LE.5) THEN                                               
C                                                                       
C  SP CASES                                                             
C                                                                       
         IF(LIJKL.EQ.0) THEN                                            
C  SSSS                                                                 
            IWK1 = IPHI +    1     ! LPHI                               
            IWK2 = IWK1 +    1     ! LWK1                               
            LAST = IWK2 +    1     ! LWK1                               
            NEED = LAST - LOADFM                                        
            CALL GETFM(NEED)                                            
            IDIM =           1                                          
            CALL SSSS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL,               
     *                 XX(IPHI),XX(IWK2),IDIM)                          
            IOFF = IWK2 +    0                                          
         ELSE IF(LIJKL.EQ.1) THEN                                       
C  PSSS                                                                 
            IWK1 = IPHI +   21                                          
            IWK2 = IWK1 +    1                                          
            LAST = IWK2 +    4                                          
            NEED = LAST - LOADFM                                        
            CALL GETFM(NEED)                                            
            IDIM =           4                                          
            CALL PSSS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL,               
     *                 XC,YC,ZC,XD,YD,ZD,                               
     *                 XX(IPHI),XX(IWK2),IDIM)                          
            IOFF = IWK2 +    1                                          
         ELSE IF(LIJKL.EQ.2) THEN                                       
C  PSPS                                                                 
            IWK1 = IPHI +   65                                          
            IWK2 = IWK1 +    4                                          
            LAST = IWK2 +   12                                          
            NEED = LAST - LOADFM                                        
            CALL GETFM(NEED)                                            
            IDIM =           4                                          
            CALL PSPS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL                
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD              
     *,                XX(IPHI),XX(IWK1),XX(IWK2),IDIM)                 
            IOFF = IWK2 +    1                                          
         ELSE IF(LIJKL.EQ.3) THEN                                       
C  PPSS                                                                 
            IWK1 = IPHI +   65                                          
            IWK2 = IWK1 +    1                                          
            LAST = IWK2 +   24                                          
            NEED = LAST - LOADFM                                        
            CALL GETFM(NEED)                                            
            IDIM =          24                                          
            CALL PPSS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL,               
     *                 XC,YC,ZC,XD,YD,ZD,                               
     *                 XX(IPHI),XX(IWK2),IDIM)                          
            IOFF = IWK2 +   15                                          
         ELSE IF(LIJKL.EQ.4) THEN                                       
C  PPPS                                                                 
            IWK1 = IPHI +  183                                          
            IWK2 = IWK1 +    4                                          
            LAST = IWK2 +   72                                          
            NEED = LAST - LOADFM                                        
            CALL GETFM(NEED)                                            
            IDIM =          24                                          
            CALL PPPS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL                
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD              
     *,                XX(IPHI),XX(IWK1),XX(IWK2),IDIM)                 
            IOFF = IWK2 +   15                                          
         ELSE IF(LIJKL.EQ.5) THEN                                       
C  PPPP                                                                 
            IWK1 = IPHI +  574                                          
            IWK2 = IWK1 +   24                                          
            LAST = IWK2 +  216                                          
            NEED = LAST - LOADFM                                        
            CALL GETFM(NEED)                                            
            IDIM =          24                                          
            CALL PPPP (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL                
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD              
     *,                XX(IPHI),XX(IWK1),XX(IWK2),IDIM)                 
            IOFF = IWK2 +   15                                          
         END IF                                                         
C                                                                       
C  END OF SP CASES                                                      
C                                                                       
      ELSE IF(LIJKL.GE.6 .AND. LIJKL.LE.20) THEN                        
C                                                                       
C  D CASES                                                              
C                                                                       
         IF(LIJKL.EQ.6) THEN                                            
C  DSSS                                                                 
            IWK1 = IPHI +   58                                          
            IWK2 = IWK1 +    1                                          
            LAST = IWK2 +   11                                          
            NEED = LAST - LOADFM                                        
            CALL GETFM(NEED)                                            
            IDIM =          11                                          
            CALL DSSS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL,               
     *                 XC,YC,ZC,XD,YD,ZD,                               
     *                 XX(IPHI),XX(IWK2),IDIM)                          
            IOFF = IWK2 +    5                                          
         ELSE IF(LIJKL.EQ.7) THEN                                       
C  DSPS                                                                 
            IWK1 = IPHI +  159                                          
            IWK2 = IWK1 +    4                                          
            LAST = IWK2 +   33                                          
            NEED = LAST - LOADFM                                        
            CALL GETFM(NEED)                                            
            IDIM =          11                                          
            CALL DSPS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL                
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD              
     *,                XX(IPHI),XX(IWK1),XX(IWK2),IDIM)                 
            IOFF = IWK2 +    5                                          
         ELSE IF(LIJKL.EQ.8) THEN                                       
C  DSPP                                                                 
            IWK1 = IPHI +  453                                          
            IWK2 = IWK1 +   11                                          
            LAST = IWK2 +  144                                          
            NEED = LAST - LOADFM                                        
            CALL GETFM(NEED)                                            
            IDIM =          24                                          
            CALL DSPP (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL                
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD              
     *,                XX(IPHI),XX(IWK1),XX(IWK2),IDIM)                 
            IOFF = IWK2 +   15                                          
         ELSE IF(LIJKL.EQ.9) THEN                                       
C  DSDS                                                                 
            IWK1 = IPHI +  389                                          
            IWK2 = IWK1 +   11                                          
            LAST = IWK2 +   66                                          
            NEED = LAST - LOADFM                                        
            CALL GETFM(NEED)                                            
            IDIM =          11                                          
            CALL DSDS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL                
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD              
     *,                XX(IPHI),XX(IWK1),XX(IWK2),IDIM)                 
            IOFF = IWK2 +    5                                          
         ELSE IF(LIJKL.EQ.10) THEN                                      
C  DPSS                                                                 
            IWK1 = IPHI +  149                                          
            IWK2 = IWK1 +    1                                          
            LAST = IWK2 +   53                                          
            NEED = LAST - LOADFM                                        
            CALL GETFM(NEED)                                            
            IDIM =          53                                          
            CALL DPSS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL,               
     *                 XC,YC,ZC,XD,YD,ZD,                               
     *                 XX(IPHI),XX(IWK2),IDIM)                          
            IOFF = IWK2 +   35                                          
         ELSE IF(LIJKL.EQ.11) THEN                                      
C  DPPS                                                                 
            IWK1 = IPHI +  390                                          
            IWK2 = IWK1 +    4                                          
            LAST = IWK2 +  159                                          
            NEED = LAST - LOADFM                                        
            CALL GETFM(NEED)                                            
            IDIM =          53                                          
            CALL DPPS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL                
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD              
     *,                XX(IPHI),XX(IWK1),XX(IWK2),IDIM)                 
            IOFF = IWK2 +   35                                          
         ELSE IF(LIJKL.EQ.12) THEN                                      
C  DPPP                                                                 
            IWK1 = IPHI + 1181                                          
            IWK2 = IWK1 +   24                                          
            LAST = IWK2 +  477                                          
            NEED = LAST - LOADFM                                        
            CALL GETFM(NEED)                                            
            IDIM =          53                                          
            CALL DPPP (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL                
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD              
     *,                XX(IPHI),XX(IWK1),XX(IWK2),IDIM)                 
            IOFF = IWK2 +   35                                          
         ELSE IF(LIJKL.EQ.13) THEN                                      
C  DPDS                                                                 
            IWK1 = IPHI +  928                                          
            IWK2 = IWK1 +   11                                          
            LAST = IWK2 +  318                                          
            NEED = LAST - LOADFM                                        
            CALL GETFM(NEED)                                            
            IDIM =          53                                          
            CALL DPDS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL                
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD              
     *,                XX(IPHI),XX(IWK1),XX(IWK2),IDIM)                 
            IOFF = IWK2 +   35                                          
         ELSE IF(LIJKL.EQ.15) THEN                                      
C  DDSS                                                                 
            IWK1 = IPHI +  319                                          
            IWK2 = IWK1 +    1                                          
            LAST = IWK2 +  165                                          
            NEED = LAST - LOADFM                                        
            CALL GETFM(NEED)                                            
            IDIM =         165                                          
            CALL DDSS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL,               
     *                 XC,YC,ZC,XD,YD,ZD,                               
     *                 XX(IPHI),XX(IWK2),IDIM)                          
            IOFF = IWK2 +  129                                          
         ELSE IF(LIJKL.EQ.16) THEN                                      
C  DDPS                                                                 
            IWK1 = IPHI +  802                                          
            IWK2 = IWK1 +    4                                          
            LAST = IWK2 +  495                                          
            NEED = LAST - LOADFM                                        
            CALL GETFM(NEED)                                            
            IDIM =         165                                          
            CALL DDPS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL                
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD              
     *,                XX(IPHI),XX(IWK1),XX(IWK2),IDIM)                 
            IOFF = IWK2 +  129                                          
         END IF                                                         
C                                                                       
C  END OF D CASES                                                       
C                                                                       
      ELSE IF(LIJKL.GE.21 .AND. LIJKL.LE.54) THEN                       
C                                                                       
C  F CASES                                                              
C                                                                       
         IF(LIJKL.EQ.21) THEN                                           
C  FSSS                                                                 
            IWK1 = IPHI +  129                                          
            IWK2 = IWK1 +    1                                          
            LAST = IWK2 +   24                                          
            NEED = LAST - LOADFM                                        
            CALL GETFM(NEED)                                            
            IDIM =          24                                          
            CALL FSSS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL,               
     *                 XC,YC,ZC,XD,YD,ZD,                               
     *                 XX(IPHI),XX(IWK2),IDIM)                          
            IOFF = IWK2 +   14                                          
         ELSE IF(LIJKL.EQ.22) THEN                                      
C  FSPS                                                                 
            IWK1 = IPHI +  328                                          
            IWK2 = IWK1 +    4                                          
            LAST = IWK2 +   72                                          
            NEED = LAST - LOADFM                                        
            CALL GETFM(NEED)                                            
            IDIM =          24                                          
            CALL FSPS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL                
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD              
     *,                XX(IPHI),XX(IWK1),XX(IWK2),IDIM)                 
            IOFF = IWK2 +   14                                          
         ELSE IF(LIJKL.EQ.23) THEN                                      
C  FSPP                                                                 
            IWK1 = IPHI +  976                                          
            IWK2 = IWK1 +   24                                          
            LAST = IWK2 +  216                                          
            NEED = LAST - LOADFM                                        
            CALL GETFM(NEED)                                            
            IDIM =          24                                          
            CALL FSPP (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL                
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD              
     *,                XX(IPHI),XX(IWK1),XX(IWK2),IDIM)                 
            IOFF = IWK2 +   14                                          
         ELSE IF(LIJKL.EQ.24) THEN                                      
C  FSDS                                                                 
            IWK1 = IPHI +  768                                          
            IWK2 = IWK1 +   11                                          
            LAST = IWK2 +  144                                          
            NEED = LAST - LOADFM                                        
            CALL GETFM(NEED)                                            
            IDIM =          24                                          
            CALL FSDS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL                
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD              
     *,                XX(IPHI),XX(IWK1),XX(IWK2),IDIM)                 
            IOFF = IWK2 +   14                                          
         ELSE IF(LIJKL.EQ.28) THEN                                      
C  FPSS                                                                 
            IWK1 = IPHI +  299                                          
            IWK2 = IWK1 +    1                                          
            LAST = IWK2 +  100                                          
            NEED = LAST - LOADFM                                        
            CALL GETFM(NEED)                                            
            IDIM =         100                                          
            CALL FPSS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL,               
     *                 XC,YC,ZC,XD,YD,ZD,                               
     *                 XX(IPHI),XX(IWK2),IDIM)                          
            IOFF = IWK2 +   70                                          
         ELSE IF(LIJKL.EQ.29) THEN                                      
C  FPPS                                                                 
            IWK1 = IPHI +  740                                          
            IWK2 = IWK1 +    4                                          
            LAST = IWK2 +  300                                          
            NEED = LAST - LOADFM                                        
            CALL GETFM(NEED)                                            
            IDIM =         100                                          
            CALL FPPS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL                
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD              
     *,                XX(IPHI),XX(IWK1),XX(IWK2),IDIM)                 
            IOFF = IWK2 +   70                                          
         ELSE IF(LIJKL.EQ.36) THEN                                      
C  FDSS                                                                 
            IWK1 = IPHI +  589                                          
            IWK2 = IWK1 +    1                                          
            LAST = IWK2 +  285                                          
            NEED = LAST - LOADFM                                        
            CALL GETFM(NEED)                                            
            IDIM =         285                                          
            CALL FDSS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL,               
     *                 XC,YC,ZC,XD,YD,ZD,                               
     *                 XX(IPHI),XX(IWK2),IDIM)                          
            IOFF = IWK2 +  225                                          
         END IF                                                         
C                                                                       
C  END OF F CASES                                                       
C                                                                       
      ELSE IF(LIJKL.GE.55 .AND. LIJKL.LE.119) THEN                      
C                                                                       
C  G CASES                                                              
C                                                                       
         IF(LIJKL.EQ.55) THEN                                           
C  GSSS                                                                 
            IWK1 = IPHI +  255                                          
            IWK2 = IWK1 +    1                                          
            LAST = IWK2 +   46                                          
            NEED = LAST - LOADFM                                        
            CALL GETFM(NEED)                                            
            IDIM =          46                                          
            CALL GSSS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL,               
     *                 XC,YC,ZC,XD,YD,ZD,                               
     *                 XX(IPHI),XX(IWK2),IDIM)                          
            IOFF = IWK2 +   31                                          
         ELSE IF(LIJKL.EQ.56) THEN                                      
C  GSPS                                                                 
            IWK1 = IPHI +  612                                          
            IWK2 = IWK1 +    4                                          
            LAST = IWK2 +  138                                          
            NEED = LAST - LOADFM                                        
            CALL GETFM(NEED)                                            
            IDIM =          46                                          
            CALL GSPS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL                
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD              
     *,                XX(IPHI),XX(IWK1),XX(IWK2),IDIM)                 
            IOFF = IWK2 +   31                                          
         ELSE IF(LIJKL.EQ.66) THEN                                      
C  GPSS                                                                 
            IWK1 = IPHI +  545                                          
            IWK2 = IWK1 +    1                                          
            LAST = IWK2 +  171                                          
            NEED = LAST - LOADFM                                        
            CALL GETFM(NEED)                                            
            IDIM =         171                                          
            CALL GPSS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL,               
     *                 XC,YC,ZC,XD,YD,ZD,                               
     *                 XX(IPHI),XX(IWK2),IDIM)                          
            IOFF = IWK2 +  126                                          
         END IF                                                         
C                                                                       
C  END OF G CASES                                                       
C                                                                       
      END IF    ! LIJKL                                                 
C                                                                       
C  ANGULAR NORMALIZATION AND SAVE TO OUTPUT ARRAY                       
C  WITH HONDO INDEXING AND REORDERING FOR GAMESS                        
C                                                                       
      MINI = KMIN(INW)                                                  
      MAXI = KMAX(INW)                                                  
      MINJ = KMIN(JNW)                                                  
      MAXJ = KMAX(JNW)                                                  
      MINK = KMIN(KNW)                                                  
      MAXK = KMAX(KNW)                                                  
      MINL = KMIN(LNW)                                                  
      MAXL = KMAX(LNW)                                                  
      LENI = MAXI - MINI + 1                                            
      LENK = MAXK - MINK + 1                                            
      II = 1                                                            
      DO I = MINI, MAXI                                                 
         IO = IORD(I)                                                   
         AI = ANGL(IO)                                                  
         IJ = II                                                        
         DO J = MINJ, MAXJ                                              
            JO = IORD(J)                                                
            AIJ = ANGL(JO)*AI                                           
            JC = ((JO-MINJ)*LENI + IO-MINI)*IDIM + IOFF                 
            IJK = IJ                                                    
            DO K = MINK, MAXK                                           
               KO = IORD(K)                                             
               AIJK = ANGL(KO)*AIJ                                      
               IJKL = IJK                                               
               DO L = MINL, MAXL                                        
                  LO = IORD(L)                                          
                  AIJKL = ANGL(LO)*AIJK                                 
                  IR = (LO-MINL)*LENK + KO-MINK                         
                  ERI(IJKL) = XX(JC+IR)*AIJKL                           
                  IJKL = IJKL + LSTRL                                   
               END DO                                                   
               IJK = IJK + LSTRK                                        
            END DO                                                      
            IJ = IJ + LSTRJ                                             
         END DO                                                         
         II = II + LSTRI                                                
      END DO                                                            
      CALL RETFM(NEED)                                                  
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2C   *DECK ERIC_TS                                          
C>                                                                      
C>    @brief   ERIC: Electron Repulsion Integral Calculator             
C>                                                                      
C>    @details ERIC_TS is a thread-safe verion of ERIC subroutine       
C>             G.D.Fletcher  Int.J.Quantum Chem. 106, 355-360(2006)     
C>                                                                      
C>    @author  Graham Fletcher, 2004, using machine generated formulae. 
C>             Jose Sierra, 2013, uses human-created simplifications.   
C>                                                                      
      SUBROUTINE ERIC_TS (ISH,JSH,KSH,LSH,ERI)                          
      use mx_limits, only: mxatm,mxgtot,mxsh,mxgsh,mxg2                 
C ----------------------------------------------------------------------
C                                                                       
C              ELECTRON REPULSION INTEGRAL CALCULATOR                   
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
C  METHOD:                                                              
C     A) FORM SCALED, CONTRACTED 1-CENTER PRECURSOR INTEGRALS.          
C        CONVERT THESE TO 4-CENTER INTEGRALS OVER CARTESIAN             
C        GAUSSIANS USING,                                               
C     B) PRECURSOR-HERMITE TRANSFER EQUATION (PTE)                      
C     C) CONTRACTED TRANSFER  EQUATION (CTE)                            
C     D) HORIZONTAL RECURSION RELATION (HRR)                            
C                                                                       
C                                                                       
C  REFERENCES                                                           
C                                                                       
C     PTE:                                                              
C     "Recursion Formula for Electron Repulsion Integrals Over          
C     Hermite Polynomials"                                              
C         G.D.Fletcher  Int.J.Quantum Chem. 106, 355-360(2006)          
C                                                                       
C     HRR:                                                              
C     M. HEAD-GORDON & J. A. POPLE,                                     
C     J. CHEM. PHYS., 89, 5777-5786 (1988).                             
C                                                                       
C     CTE:                                                              
C     P. M. W. GILL, M. HEAD-GORDON, & J. A. POPLE,                     
C     INT. J. Q. CHEM., SYMP. 23, 269-280 (1989).                       
C                                                                       
C     INTERPOLATION METHOD:                                             
C     P. M. W. GILL, B. G. JOHNSON, & J. A. POPLE,                      
C     INT. J. Q. CHEM., 40, 745-752 (1991).                             
C                                                                       
C                                                                       
C  THIS IS THE DRIVER ROUTINE THAT INTERFACES TO GAMESS COMMON NSHEL    
C                                                                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C                                                                       
      DOUBLE PRECISION ERI(*)                                           
C                                                                       
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,                    
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)                  
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),      
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),      
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),   
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL           
C                                                                       
C  /ERIPRM/ SYMBOLS:                                                    
C     EXI,J,K,L   = PRIMITIVE EXPONENTS OF I,J,K,L SHELLS               
C     CCI,J,K,L   = CONTRACTION COEFFICIENTS OF ISH,JSH,KSH,LSH         
C     XAB,YAB,ZAB = EXPONENT-WEIGHTED DISTANCE FOR BRA                  
C     XCD,YCD,ZCD = EXPONENT-WEIGHTED DISTANCE FOR KET                  
C     CCBRA,CCKET = PRIMITIVE CHARGE-CLOUD FACTOR FOR BRA,KET           
C     SLBRA,SLKET = FACTOR CONVERTING SMALL-T EXPRESSION TO LARGE-T     
C                                                                       
      COMMON /ERIPRM/ EXI(MXGSH),EXJ(MXGSH),EXK(MXGSH),EXL(MXGSH),      
     *                CCI(MXGSH),CCJ(MXGSH),CCK(MXGSH),CCL(MXGSH),      
     *                XAB(MXG2),YAB(MXG2),ZAB(MXG2),                    
     *                XCD(MXG2),YCD(MXG2),ZCD(MXG2),                    
     *                CCBRA(MXG2),CCKET(MXG2),RXB(MXG2),                
     *                SLBRA(MXG2),SLKET(MXG2),RXK(MXG2)                 
C$omp threadprivate(/ERIPRM/)
C                                                                       
      COMMON /ERIDAT/ LEN1,LEN2,LEN3,LEN4                               
      COMMON /ERIOUT/ INW,JNW,KNW,LNW,LSTRI,LSTRJ,LSTRK,LSTRL           
C$omp threadprivate(/ERIOUT/)
C                                                                       
      COMMON /FMCOM /XX(1)                                              
C                                                                       
      LOGICAL IEQJ,KEQL                                                 
C                                                                       
      PARAMETER (PT5=0.5D+00)                                           
      PARAMETER (ONE=1.0D+00)                                           
      PARAMETER (TWO=2.0D+00)                                           
C                                                                       
C     SR3 <==> SQRT( 3)   <==> 1.7320508075688773D+00                   
C     SR5 <==> SQRT( 5)   <==> 2.2360679774997897D+00                   
C     SR7 <==> SQRT( 7)   <==> 2.6457513110645906D+00                   
C     S15 <==> SQRT(15)   <==> 3.8729833462074169D+00                   
C     S35 <==> SQRT(35)   <==> 5.9160797830996160D+00                   
C     S53 <==> SQRT(35/3) <==> 3.4156502553198661D+00                   
C                                                                       
CC    PARAMETER (SR3=1.7320508075688772D+00,                            
CC   *           SR5=2.2360679774997902D+00,                            
CC   *           SR7=2.6457513110645903D+00,                            
CC   *           S15=3.8729833462074170D+00,                            
CC   *           S35=5.9160797830996170D+00,                            
CC   *           S53=3.4156502553198673D+00)                            
      PARAMETER (SR3=1.7320508075688773D+00)                            
      PARAMETER (SR5=2.2360679774997897D+00)                            
      PARAMETER (SR7=2.6457513110645906D+00)                            
      PARAMETER (S15=SR3*SR5)                                           
      PARAMETER (S35=SR5*SR7)                                           
      PARAMETER (S53=S35/SR3)                                           
C                                                                       
CJMS  PARAMETER (PI =3.14159265358979323846D+00)                        
CC    PARAMETER (PI =3.1415926535897932D+00)                            
C                                                                       
CC               PI254 = PI * SQRT( 2 * SQRT(PI ))                      
CC               PI214 = PI254 /(PI +PI )                               
CC    PARAMETER (PI254 =5.914967172795613D+00)                          
CC    PARAMETER (PI214 =0.9413962637767148D+00)                         
      PARAMETER (PI254 =5.9149671727956129D+00)                         
      PARAMETER (PI214 =0.94139626377671481D+00)                        
C                                                                       
      DIMENSION xxiphi(1182), xxiwk1(25), xxiwk2(496)                   
      SAVE xxiphi, xxiwk1, xxiwk2                                       
!$omp threadprivate(xxiphi, xxiwk1, xxiwk2)                             
      DIMENSION IORD(35),ANGL(35)                                       
      DATA IORD/                                                        
     1       1,                                                         
     2       2,  3,  4,                                                 
     3       5,  7, 10,  6,  8,  9,                                     
     4      11, 14, 20, 12, 15, 13, 17, 18, 19, 16,                     
     5      21, 25, 35, 22, 26, 24, 29, 33, 34, 23, 30, 32, 27, 28, 31/ 
      DATA ANGL/                                                        
     1     ONE,                                                         
     2     ONE,ONE,ONE,                                                 
     3     ONE,SR3,ONE,SR3,SR3,ONE,                                     
     4     ONE,SR5,SR5,ONE,SR5,S15,SR5,SR5,SR5,ONE,                     
     5     ONE,SR7,S53,SR7,ONE,SR7,S35,S35,SR7,S53,S35,S53,SR7,SR7,ONE/ 
C                                                                       
CC    DATA ANGL/1.0D+00,  1.0D+00, 1.0D+00, 1.0D+00,                    
CC   * 1.0D+00, 1.7320508075688772D+00, 1.0D+00, 1.7320508075688772D+00,
CC   * 1.7320508075688772D+00, 1.0D+00,                                 
CC   * 1.0D+00, 2.2360679774997902D+00, 2.2360679774997902D+00, 1.0D+00,
CC   * 2.2360679774997902D+00, 3.8729833462074170D+00,                  
CC   * 2.2360679774997902D+00,                                          
CC   * 2.2360679774997902D+00, 2.2360679774997902D+00, 1.0D+00,         
CC   * 1.0D+00, 2.6457513110645903D+00, 3.4156502553198673D+00,         
CC   * 2.6457513110645903D+00, 1.0D+00, 2.6457513110645903D+00,         
CC   * 5.9160797830996170D+00, 5.9160797830996170D+00,                  
CC   * 2.6457513110645903D+00,                                          
CC   * 3.4156502553198673D+00, 5.9160797830996170D+00,                  
CC   * 3.4156502553198673D+00,                                          
CC   * 2.6457513110645903D+00, 2.6457513110645903D+00, 1.0D+00/         
C                                                                       
C  FAST CODES                                                           
C                                                                       
C  RE-ORDER CHARGE-CLOUD SHELLS                                         
C  (HRR DATA ONLY STORED FOR LI>=LJ)                                    
C                                                                       
      INW = ISH                                                         
      JNW = JSH                                                         
      KNW = KSH                                                         
      LNW = LSH                                                         
      LSTRI = LEN4                                                      
      LSTRJ = LEN3                                                      
      LSTRK = LEN2                                                      
      LSTRL = LEN1                                                      
      LANGI = KTYPE(INW) - 1                                            
      LANGJ = KTYPE(JNW) - 1                                            
      LANGK = KTYPE(KNW) - 1                                            
      LANGL = KTYPE(LNW) - 1                                            
      IF(LANGI.LT.LANGJ) THEN                                           
         INW   = JSH                                                    
         JNW   = ISH                                                    
         ITMP  = LSTRI                                                  
         LSTRI = LSTRJ                                                  
         LSTRJ = ITMP                                                   
      END IF                                                            
      IF(LANGK.LT.LANGL) THEN                                           
         KNW   = LSH                                                    
         LNW   = KSH                                                    
         ITMP  = LSTRK                                                  
         LSTRK = LSTRL                                                  
         LSTRL = ITMP                                                   
      END IF                                                            
C                                                                       
C  SWAP CHARGE CLOUDS FOR EFFICIENCY...                                 
C                                                                       
      LBRA  = LANGI + LANGJ                                             
      LKET  = LANGK + LANGL                                             
      IF(LBRA.GT.LKET) THEN                                             
         ITMP  = INW                                                    
         INW   = KNW                                                    
         KNW   = ITMP                                                   
         ITMP  = LSTRI                                                  
         LSTRI = LSTRK                                                  
         LSTRK = ITMP                                                   
         ITMP  = JNW                                                    
         JNW   = LNW                                                    
         LNW   = ITMP                                                   
         ITMP  = LSTRJ                                                  
         LSTRJ = LSTRL                                                  
         LSTRL = ITMP                                                   
C                                                                       
C  ...OR TO SELECT THE FAST CODE                                        
C                                                                       
      ELSE IF(LBRA.EQ.LKET) THEN                                        
         LANGI = KTYPE(INW)                                             
         LANGJ = KTYPE(JNW)                                             
         LANGK = KTYPE(KNW)                                             
         LANGL = KTYPE(LNW)                                             
         LNGIJ = (LANGI*LANGI-LANGI)/2 + LANGJ                          
         LNGKL = (LANGK*LANGK-LANGK)/2 + LANGL                          
         IF(LNGKL.GT.LNGIJ) THEN                                        
            ITMP  = INW                                                 
            INW   = KNW                                                 
            KNW   = ITMP                                                
            ITMP  = LSTRI                                               
            LSTRI = LSTRK                                               
            LSTRK = ITMP                                                
            ITMP  = JNW                                                 
            JNW   = LNW                                                 
            LNW   = ITMP                                                
            ITMP  = LSTRJ                                               
            LSTRJ = LSTRL                                               
            LSTRL = ITMP                                                
         END IF                                                         
      END IF                                                            
C                                                                       
C  MAIN PARAMETERS                                                      
C                                                                       
      LANGI = KTYPE(INW) - 1                                            
      LANGJ = KTYPE(JNW) - 1                                            
      LANGK = KTYPE(KNW) - 1                                            
      LANGL = KTYPE(LNW) - 1                                            
      LBRA  = LANGI + LANGJ                                             
      LKET  = LANGK + LANGL                                             
      IPRIM = KNG(INW)                                                  
      JPRIM = KNG(JNW)                                                  
      KPRIM = KNG(KNW)                                                  
      LPRIM = KNG(LNW)                                                  
C                                                                       
      IEQJ  = INW.EQ.JNW                                                
      KEQL  = KNW.EQ.LNW                                                
C                                                                       
C  SET UP 1-INDEX PARAMETERS                                            
C                                                                       
      K1 = KSTART(INW)                                                  
      K2 = K1 + IPRIM - 1                                               
      II = 0                                                            
      DO I = K1, K2                                                     
         II = II + 1                                                    
         EXI(II) = EX(I)                                                
      END DO                                                            
      ITYP = KTYPE(INW)                                                 
      IF(ITYP.EQ.1) THEN                                                
         II = 0                                                         
         DO I = K1, K2                                                  
            II = II + 1                                                 
            CCI(II) = CS(I)                                             
         END DO                                                         
      ELSE IF(ITYP.EQ.2) THEN                                           
         II = 0                                                         
         DO I = K1, K2                                                  
            II = II + 1                                                 
            CCI(II) = CP(I)                                             
         END DO                                                         
      ELSE IF(ITYP.EQ.3) THEN                                           
         II = 0                                                         
         DO I = K1, K2                                                  
            II = II + 1                                                 
            CCI(II) = CD(I)                                             
         END DO                                                         
      ELSE IF(ITYP.EQ.4) THEN                                           
         II = 0                                                         
         DO I = K1, K2                                                  
            II = II + 1                                                 
            CCI(II) = CF(I)                                             
         END DO                                                         
      ELSE IF(ITYP.EQ.5) THEN                                           
         II = 0                                                         
         DO I = K1, K2                                                  
            II = II + 1                                                 
            CCI(II) = CG(I)                                             
         END DO                                                         
      END IF                                                            
C                                                                       
      K1 = KSTART(JNW)                                                  
      K2 = K1 + JPRIM - 1                                               
      JJ = 0                                                            
      DO J = K1, K2                                                     
         JJ = JJ + 1                                                    
         EXJ(JJ) = EX(J)                                                
      END DO                                                            
      JTYP = KTYPE(JNW)                                                 
      IF(JTYP.EQ.1) THEN                                                
         JJ = 0                                                         
         DO J = K1, K2                                                  
            JJ = JJ + 1                                                 
            CCJ(JJ) = CS(J)                                             
         END DO                                                         
      ELSE IF(JTYP.EQ.2) THEN                                           
         JJ = 0                                                         
         DO J = K1, K2                                                  
            JJ = JJ + 1                                                 
            CCJ(JJ) = CP(J)                                             
         END DO                                                         
      ELSE IF(JTYP.EQ.3) THEN                                           
         JJ = 0                                                         
         DO J = K1, K2                                                  
            JJ = JJ + 1                                                 
            CCJ(JJ) = CD(J)                                             
         END DO                                                         
      ELSE IF(JTYP.EQ.4) THEN                                           
         JJ = 0                                                         
         DO J = K1, K2                                                  
            JJ = JJ + 1                                                 
            CCJ(JJ) = CF(J)                                             
         END DO                                                         
      ELSE IF(JTYP.EQ.5) THEN                                           
         JJ = 0                                                         
         DO J = K1, K2                                                  
            JJ = JJ + 1                                                 
            CCJ(JJ) = CG(J)                                             
         END DO                                                         
      END IF                                                            
C                                                                       
      K1 = KSTART(KNW)                                                  
      K2 = K1 + KPRIM - 1                                               
      KK = 0                                                            
      DO K = K1, K2                                                     
         KK = KK + 1                                                    
         EXK(KK) = EX(K)                                                
      END DO                                                            
      KTYP = KTYPE(KNW)                                                 
      IF(KTYP.EQ.1) THEN                                                
         KK = 0                                                         
         DO K = K1, K2                                                  
            KK = KK + 1                                                 
            CCK(KK) = CS(K)                                             
         END DO                                                         
      ELSE IF(KTYP.EQ.2) THEN                                           
         KK = 0                                                         
         DO K = K1, K2                                                  
            KK = KK + 1                                                 
            CCK(KK) = CP(K)                                             
         END DO                                                         
      ELSE IF(KTYP.EQ.3) THEN                                           
         KK = 0                                                         
         DO K = K1, K2                                                  
            KK = KK + 1                                                 
            CCK(KK) = CD(K)                                             
         END DO                                                         
      ELSE IF(KTYP.EQ.4) THEN                                           
         KK = 0                                                         
         DO K = K1, K2                                                  
            KK = KK + 1                                                 
            CCK(KK) = CF(K)                                             
         END DO                                                         
      ELSE IF(KTYP.EQ.5) THEN                                           
         KK = 0                                                         
         DO K = K1, K2                                                  
            KK = KK + 1                                                 
            CCK(KK) = CG(K)                                             
         END DO                                                         
      END IF                                                            
C                                                                       
      K1 = KSTART(LNW)                                                  
      K2 = K1 + LPRIM - 1                                               
      LL = 0                                                            
      DO L = K1, K2                                                     
         LL = LL + 1                                                    
         EXL(LL) = EX(L)                                                
      END DO                                                            
      LTYP = KTYPE(LNW)                                                 
      IF(LTYP.EQ.1) THEN                                                
         LL = 0                                                         
         DO L = K1, K2                                                  
            LL = LL + 1                                                 
            CCL(LL) = CS(L)                                             
         END DO                                                         
      ELSE IF(LTYP.EQ.2) THEN                                           
         LL = 0                                                         
         DO L = K1, K2                                                  
            LL = LL + 1                                                 
            CCL(LL) = CP(L)                                             
         END DO                                                         
      ELSE IF(LTYP.EQ.3) THEN                                           
         LL = 0                                                         
         DO L = K1, K2                                                  
            LL = LL + 1                                                 
            CCL(LL) = CD(L)                                             
         END DO                                                         
      ELSE IF(LTYP.EQ.4) THEN                                           
         LL = 0                                                         
         DO L = K1, K2                                                  
            LL = LL + 1                                                 
            CCL(LL) = CF(L)                                             
         END DO                                                         
      ELSE IF(LTYP.EQ.5) THEN                                           
         LL = 0                                                         
         DO L = K1, K2                                                  
            LL = LL + 1                                                 
            CCL(LL) = CG(L)                                             
         END DO                                                         
      END IF                                                            
C                                                                       
C  2-INDEX PARAMETERS                                                   
C                                                                       
      IATM = KATOM(INW)                                                 
      XA   = C(1,IATM)                                                  
      YA   = C(2,IATM)                                                  
      ZA   = C(3,IATM)                                                  
      JATM = KATOM(JNW)                                                 
      XB   = C(1,JATM)                                                  
      YB   = C(2,JATM)                                                  
      ZB   = C(3,JATM)                                                  
      RAB = (XA-XB)**2 + (YA-YB)**2 + (ZA-ZB)**2                        
      IJ = 0                                                            
      DO JJ = 1, JPRIM                                                  
         EJ = EXJ(JJ)                                                   
         CJ = CCJ(JJ)                                                   
         ITOP = IPRIM                                                   
         IF(IEQJ) ITOP = JJ                                             
         DO II = 1, ITOP                                                
            EI = EXI(II)                                                
            CIX = CCI(II)                                               
            EIJ = ONE/(EI+EJ)                                           
            IJ = IJ + 1                                                 
            CCFAC = PI254*CIX*CJ*EIJ*EXP( -EI*EJ*RAB*EIJ )              
            IF(IEQJ .AND. II.NE.JJ) CCFAC = CCFAC*TWO                   
            CCBRA(IJ) = CCFAC                                           
            SLBRA(IJ) = PI214*SQRT(EIJ)                                 
            XAB(IJ) = (EI*XA + EJ*XB)*EIJ                               
            YAB(IJ) = (EI*YA + EJ*YB)*EIJ                               
            ZAB(IJ) = (EI*ZA + EJ*ZB)*EIJ                               
            RXB(IJ) = EIJ*PT5                                           
         END DO                                                         
      END DO                                                            
C                                                                       
      KATM = KATOM(KNW)                                                 
      XC   = C(1,KATM)                                                  
      YC   = C(2,KATM)                                                  
      ZC   = C(3,KATM)                                                  
      LATM = KATOM(LNW)                                                 
      XD   = C(1,LATM)                                                  
      YD   = C(2,LATM)                                                  
      ZD   = C(3,LATM)                                                  
      RCD = (XC-XD)**2 + (YC-YD)**2 + (ZC-ZD)**2                        
      KL = 0                                                            
      DO LL = 1, LPRIM                                                  
         EL = EXL(LL)                                                   
         CL = CCL(LL)                                                   
         KTOP = KPRIM                                                   
         IF(KEQL) KTOP = LL                                             
         DO KK = 1, KTOP                                                
            EK = EXK(KK)                                                
            CK = CCK(KK)                                                
            EKL = ONE/(EK+EL)                                           
            KL = KL + 1                                                 
            CCFAC = PI254*CK*CL*EKL*EXP( -EK*EL*RCD*EKL )               
            IF(KEQL .AND. KK.NE.LL) CCFAC = CCFAC*TWO                   
            CCKET(KL) = CCFAC                                           
            SLKET(KL) = PI214*SQRT(EKL)                                 
            XCD(KL) = (EK*XC + EL*XD)*EKL                               
            YCD(KL) = (EK*YC + EL*YD)*EKL                               
            ZCD(KL) = (EK*ZC + EL*ZD)*EKL                               
            RXK(KL) = EKL*PT5                                           
         END DO                                                         
      END DO                                                            
C                                                                       
      CALL VALFM(LOADFM)                                                
      IPHI = LOADFM + 1                                                 
C                                                                       
C  ANGULAR MOMENTUM 4-INDEX                                             
C                                                                       
      LBGT  = MAX(LANGI,LANGJ)                                          
      LBLT  = MIN(LANGI,LANGJ)                                          
      LNGIJ = (LBGT*LBGT+LBGT)/2 + LBLT                                 
      LKGT  = MAX(LANGK,LANGL)                                          
      LKLT  = MIN(LANGK,LANGL)                                          
      LNGKL = (LKGT*LKGT+LKGT)/2 + LKLT                                 
      LQGT  = MAX(LNGKL,LNGIJ)                                          
      LQLT  = MIN(LNGKL,LNGIJ)                                          
      LIJKL = (LQGT*LQGT+LQGT)/2 + LQLT                                 
C                                                                       
      IF(LIJKL.LE.5) THEN                                               
C                                                                       
C  SP CASES                                                             
C                                                                       
         IF(LIJKL.EQ.0) THEN                                            
C  SSSS                                                                 
            IDIM =           1                                          
            IOFF = 0                                                    
            CALL SSSS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL,               
     *                 XXIPHI,XXIWK2,IDIM)                              
         ELSE IF(LIJKL.EQ.1) THEN                                       
C  PSSS                                                                 
            IDIM =           4                                          
            IOFF = 1                                                    
            CALL PSSS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL,               
     *                 XC,YC,ZC,XD,YD,ZD,                               
     *                 XXIPHI,XXIWK2,IDIM)                              
         ELSE IF(LIJKL.EQ.2) THEN                                       
C  PSPS                                                                 
            IDIM =           4                                          
            IOFF = 1                                                    
            CALL PSPS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL                
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD              
     *,                XXIPHI,XXIWK1,XXIWK2,IDIM)                       
         ELSE IF(LIJKL.EQ.3) THEN                                       
C  PPSS                                                                 
            IDIM =          24                                          
            IOFF = 15                                                   
            CALL PPSS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL,               
     *                 XC,YC,ZC,XD,YD,ZD,                               
     *                 XXIPHI,XXIWK2,IDIM)                              
         ELSE IF(LIJKL.EQ.4) THEN                                       
C  PPPS                                                                 
            IDIM =          24                                          
            IOFF = 15                                                   
            CALL PPPS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL                
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD              
     *,                XXIPHI,XXIWK1,XXIWK2,IDIM)                       
         ELSE IF(LIJKL.EQ.5) THEN                                       
C  PPPP                                                                 
            IDIM =          24                                          
            IOFF = 15                                                   
            CALL PPPP (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL                
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD              
     *,                XXIPHI,XXIWK1,XXIWK2,IDIM)                       
         END IF                                                         
C                                                                       
C  END OF SP CASES                                                      
C                                                                       
      ELSE IF(LIJKL.GE.6 .AND. LIJKL.LE.20) THEN                        
C                                                                       
C  D CASES                                                              
C                                                                       
         IF(LIJKL.EQ.6) THEN                                            
C  DSSS                                                                 
            IDIM =          11                                          
            IOFF = 5                                                    
            CALL DSSS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL,               
     *                 XC,YC,ZC,XD,YD,ZD,                               
     *                 XXIPHI,XXIWK2,IDIM)                              
         ELSE IF(LIJKL.EQ.7) THEN                                       
C  DSPS                                                                 
            IDIM =          11                                          
            IOFF = 5                                                    
            CALL DSPS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL                
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD              
     *,                XXIPHI,XXIWK1,XXIWK2,IDIM)                       
         ELSE IF(LIJKL.EQ.8) THEN                                       
C  DSPP                                                                 
            IDIM =          24                                          
            IOFF = 15                                                   
            CALL DSPP (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL                
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD              
     *,                XXIPHI,XXIWK1,XXIWK2,IDIM)                       
         ELSE IF(LIJKL.EQ.9) THEN                                       
C  DSDS                                                                 
            IDIM =          11                                          
            IOFF = 5                                                    
            CALL DSDS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL                
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD              
     *,                XXIPHI,XXIWK1,XXIWK2,IDIM)                       
         ELSE IF(LIJKL.EQ.10) THEN                                      
C  DPSS                                                                 
            IDIM =          53                                          
            IOFF = 35                                                   
            CALL DPSS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL,               
     *                 XC,YC,ZC,XD,YD,ZD,                               
     *                 XXIPHI,XXIWK2,IDIM)                              
         ELSE IF(LIJKL.EQ.11) THEN                                      
C  DPPS                                                                 
            IDIM =          53                                          
            IOFF = 35                                                   
            CALL DPPS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL                
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD              
     *,                XXIPHI,XXIWK1,XXIWK2,IDIM)                       
         ELSE IF(LIJKL.EQ.12) THEN                                      
C  DPPP                                                                 
            IDIM =          53                                          
            IOFF = 35                                                   
            CALL DPPP (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL                
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD              
     *,                XXIPHI,XXIWK1,XXIWK2,IDIM)                       
         ELSE IF(LIJKL.EQ.13) THEN                                      
C  DPDS                                                                 
            IDIM =          53                                          
            IOFF = 35                                                   
            CALL DPDS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL                
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD              
     *,                XXIPHI,XXIWK1,XXIWK2,IDIM)                       
         ELSE IF(LIJKL.EQ.15) THEN                                      
C  DDSS                                                                 
            IDIM =         165                                          
            IOFF = 129                                                  
            CALL DDSS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL,               
     *                 XC,YC,ZC,XD,YD,ZD,                               
     *                 XXIPHI,XXIWK2,IDIM)                              
         ELSE IF(LIJKL.EQ.16) THEN                                      
C  DDPS                                                                 
            IDIM =         165                                          
            IOFF = 129                                                  
            CALL DDPS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL                
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD              
     *,                XXIPHI,XXIWK1,XXIWK2,IDIM)                       
         END IF                                                         
C                                                                       
C  END OF D CASES                                                       
C                                                                       
      ELSE IF(LIJKL.GE.21 .AND. LIJKL.LE.54) THEN                       
C                                                                       
C  F CASES                                                              
C                                                                       
         IF(LIJKL.EQ.21) THEN                                           
C  FSSS                                                                 
            IDIM =          24                                          
            IOFF = 14                                                   
            CALL FSSS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL,               
     *                 XC,YC,ZC,XD,YD,ZD,                               
     *                 XXIPHI,XXIWK2,IDIM)                              
         ELSE IF(LIJKL.EQ.22) THEN                                      
C  FSPS                                                                 
            IDIM =          24                                          
            IOFF = 14                                                   
            CALL FSPS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL                
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD              
     *,                XXIPHI,XXIWK1,XXIWK2,IDIM)                       
         ELSE IF(LIJKL.EQ.23) THEN                                      
C  FSPP                                                                 
            IDIM =          24                                          
            IOFF = 14                                                   
            CALL FSPP (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL                
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD              
     *,                XXIPHI,XXIWK1,XXIWK2,IDIM)                       
         ELSE IF(LIJKL.EQ.24) THEN                                      
C  FSDS                                                                 
            IDIM =          24                                          
            IOFF = 14                                                   
            CALL FSDS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL                
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD              
     *,                XXIPHI,XXIWK1,XXIWK2,IDIM)                       
         ELSE IF(LIJKL.EQ.28) THEN                                      
C  FPSS                                                                 
            IDIM =         100                                          
            IOFF = 70                                                   
            CALL FPSS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL,               
     *                 XC,YC,ZC,XD,YD,ZD,                               
     *                 XXIPHI,XXIWK2,IDIM)                              
         ELSE IF(LIJKL.EQ.29) THEN                                      
C  FPPS                                                                 
            IDIM =         100                                          
            IOFF = 70                                                   
            CALL FPPS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL                
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD              
     *,                XXIPHI,XXIWK1,XXIWK2,IDIM)                       
         ELSE IF(LIJKL.EQ.36) THEN                                      
C  FDSS                                                                 
            IDIM =         285                                          
            IOFF = 225                                                  
            CALL FDSS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL,               
     *                 XC,YC,ZC,XD,YD,ZD,                               
     *                 XXIPHI,XXIWK2,IDIM)                              
         END IF                                                         
C                                                                       
C  END OF F CASES                                                       
C                                                                       
      ELSE IF(LIJKL.GE.55 .AND. LIJKL.LE.119) THEN                      
C                                                                       
C  G CASES                                                              
C                                                                       
         IF(LIJKL.EQ.55) THEN                                           
C  GSSS                                                                 
            IDIM =          46                                          
            IOFF = 31                                                   
            CALL GSSS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL,               
     *                 XC,YC,ZC,XD,YD,ZD,                               
     *                 XXIPHI,XXIWK2,IDIM)                              
         ELSE IF(LIJKL.EQ.56) THEN                                      
C  GSPS                                                                 
            IDIM =          46                                          
            IOFF = 31                                                   
            CALL GSPS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL                
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD              
     *,                XXIPHI,XXIWK1,XXIWK2,IDIM)                       
         ELSE IF(LIJKL.EQ.66) THEN                                      
C  GPSS                                                                 
            IDIM =         171                                          
            IOFF = 126                                                  
            CALL GPSS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL,               
     *                 XC,YC,ZC,XD,YD,ZD,                               
     *                 XXIPHI,XXIWK2,IDIM)                              
         END IF                                                         
C                                                                       
C  END OF G CASES                                                       
C                                                                       
      END IF    ! LIJKL                                                 
C                                                                       
C  ANGULAR NORMALIZATION AND SAVE TO OUTPUT ARRAY                       
C  WITH HONDO INDEXING AND REORDERING FOR GAMESS                        
C                                                                       
      MINI = KMIN(INW)                                                  
      MAXI = KMAX(INW)                                                  
      MINJ = KMIN(JNW)                                                  
      MAXJ = KMAX(JNW)                                                  
      MINK = KMIN(KNW)                                                  
      MAXK = KMAX(KNW)                                                  
      MINL = KMIN(LNW)                                                  
      MAXL = KMAX(LNW)                                                  
      LENI = MAXI - MINI + 1                                            
      LENK = MAXK - MINK + 1                                            
      II = 1                                                            
      DO I = MINI, MAXI                                                 
         IO = IORD(I)                                                   
         AI = ANGL(IO)                                                  
         IJ = II                                                        
         DO J = MINJ, MAXJ                                              
            JO = IORD(J)                                                
            AIJ = ANGL(JO)*AI                                           
            JC = ((JO-MINJ)*LENI + IO-MINI)*IDIM + IOFF + 1             
            IJK = IJ                                                    
            DO K = MINK, MAXK                                           
               KO = IORD(K)                                             
               AIJK = ANGL(KO)*AIJ                                      
               IJKL = IJK                                               
               DO L = MINL, MAXL                                        
                  LO = IORD(L)                                          
                  AIJKL = ANGL(LO)*AIJK                                 
                  IR = (LO-MINL)*LENK + KO-MINK                         
                  ERI(IJKL) = XXIWK2(JC+IR)*AIJKL                       
                  IJKL = IJKL + LSTRL                                   
               END DO                                                   
               IJK = IJK + LSTRK                                        
            END DO                                                      
            IJ = IJ + LSTRJ                                             
         END DO                                                         
         II = II + LSTRI                                                
      END DO                                                            
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2C   *DECK ERIPRE                                           
C>                                                                      
C>    @brief   ERIC pre-initialization                                  
C>                                                                      
C>    @details SET UP INTEGRAL COMPUTATION DATA FOR LARGEST CASE.       
C>             THIS ROUTINE CAN ONLY BE CALLED WHEN THE HIGHEST         
C>             ANGULAR MOMENTUM (L) OF THE BASIS HAS BEEN ESTABLISHED.  
C>                                                                      
C>    @author  Graham Fletcher, 2004                                    
C>                                                                      
C>    @date Jan 2020  Peng Xu and Tosaporn Sattasathuchana              
C>                  - find the maximum coefficient of s,p,d shells      
C>                    for EFP (CMAXEF), used for screening in           
C>                    rotated axis integral scheme for QMEFP2           
                                                                        
      SUBROUTINE ERIPRE                                                 
      use mx_limits, only: mxgtot,mxsh,mxgsh,mxg2                       
      USE comm_EFPBAS                                                   
      USE comm_FRGINF, only: NFRG                                       
      USE comm_FRGTYP, only: ISET                                       
C                                                                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C                                                                       
      LOGICAL                                     GOPARR,DSKWRK,MASWRK  
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
      COMMON /IOFILE/ IR,IW,IP,IJKO,IJKT,IDAF,NAV,IODA(950)             
C                                                                       
C  /FMTTBL/ SYMBOLS:                                                    
C  LENGTHS                                                              
C     NTX   = MAX ORDER OF INTERPOLATION                                
C     NPF   = NUMBER OF GRID POINTS FOR GAMMA                           
C     NGRD  = NUMBER OF GRIDS FOR GAMMA                                 
C     NPX   = NUMBER OF GRID POINTS FOR EXP(-X)                         
C     MXQT  = MAX L OF QUARTET                                          
C  GRID DATA                                                            
C     FGRID = INTERPOLATION TABLES FOR GAMMA                            
C     RFINC = RECIPROCAL INTERVAL WIDTHS FOR GAMMA                      
C     XGRID = INTERPOLATION TABLE FOR EXP(-X)                           
C     RXINC = RECIPROCAL INTERVAL WIDTH FOR EXP(-X)                     
C     TMAX  = MAXIMUM VALUE TO INTERPOLATE FOR                          
C     NORD  = ORDER OF INTERPOLATION                                    
C  FACTORS                                                              
C     RMR   = 1/(2M-1) FOR FM(T) RECURSION                              
C     TLGM  = (2M-1)!! FACTORS FOR LARGE-T FORMULA                      
C                                                                       
      PARAMETER (NTX=4)                                                 
      PARAMETER (NPF=450)                                               
      PARAMETER (NGRD=7)                                                
      PARAMETER (NPX=1000)                                              
      PARAMETER (MXQT=16)                                               
      COMMON /FMTTBL/ FGRID(0:NTX,0:NPF,0:NGRD),XGRID(0:NTX,0:NPX)      
     *,               TMAX,RFINC(0:NGRD),RXINC                          
     *,               RMR(MXQT),TLGM(0:MXQT),NORD                       
C                                                                       
      LOGICAL                                                SECONDD    
      COMMON /CXTHRM/ CXTHERM(11),CXZPE,METHCX,ICXBAS,ICXPCM,SECONDD    
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),      
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),      
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),   
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL           
C                                                                       
C  VARIOUS TOLERENCES FOR ROTATED AXIS CODE                             
C                                                                       
      COMMON /MAXC  / CMAX(MXGTOT),CMAXA(MXGSH),CMAXB(MXGSH),           
     *                CMAXC(MXGSH),CMAXD(MXGSH),ISMLP(MXG2),ISMLQ       
C$omp threadprivate(/MAXC/)
      LOGICAL                           OUT                             
      COMMON /SHLT  / TOL,CUTOFF,ICOUNT,OUT                             
      COMMON /INTAC2/ EI1,EI2,CUX                                       
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL                
C                                                                       
      COMMON /ERIDAT/ LEN1,LEN2,LEN3,LEN4                               
C                                                                       
      PARAMETER (EC1= 1.0D-02)                                          
      PARAMETER (EC2= 1.0D-04)                                          
      PARAMETER (RLN10=2.30258D+00)                                     
      PARAMETER (CX1=25.0D+00)                                          
C                                                                       
      PARAMETER (ONE=1.0D+00)                                           
C                                                                       
      LOGICAL FIRST,RDGRID                                              
      DATA FIRST/.TRUE./                                                
      DATA COMP/8HCOMP    /                                             
      SAVE FIRST                                                        
C                                                                       
      IF(RUNTYP.EQ.COMP) THEN                                           
         FIRST=.TRUE.                                                   
         IF(SECONDD) FIRST=.FALSE.                                      
      END IF                                                            
      IF(FIRST .AND. MASWRK) WRITE(IW,9900)                             
C                                                                       
C  THE ERIC PACKAGE ONLY GOES UP TO G SHELLS (LMAX=4), AT MOST.         
C                                                                       
      CALL BASCHK(LMAX)                                                 
      LMAX = MIN(LMAX,4)                                                
      MMAX = MAX(4,LMAX*4)                                              
C                                                                       
C  LENGTH OF HIGHEST-L SHELL                                            
C  INDEX DIMENSIONS FOR HONDO-STYLE ADDRESSING                          
C  OTHER INTEGRAL CODES (I.E. HONDO) WILL PERHAPS GO TO HIGHER LMAX.    
C                                                                       
      CALL BASCHK(LMAX)                                                 
      MXLEN = MAX(4, (LMAX*LMAX+3*LMAX+2)/2 )                           
      LEN1 = 1                                                          
      LEN2 = MXLEN                                                      
      LEN3 = MXLEN**2                                                   
      LEN4 = MXLEN**3                                                   
C                                                                       
C  GET INTERPOLATION GRIDS                                              
C  NOTE THAT ERICMP IS HARDWIRED FOR 4TH ORDER                          
C                                                                       
      IF(FIRST) THEN                                                    
         FIRST  = .FALSE.                                               
         SECONDD= .TRUE.                                                
         RDGRID = .TRUE.                                                
         IF(RDGRID) THEN                                                
C                                                                       
C  USE PREVIOUSLY GENERATED GRIDS: 4TH ORDER, 10**-12 ACCURACY, TMAX=25 
C                                                                       
            CALL GRIDIN                                                 
            NORD = 4                                                    
         ELSE                                                           
C                                                                       
C  GENERATE INTERPOLATION GRIDS                                         
C  THIS ADDS A COUPLE OF SECONDS TO THE TOTAL EXECUTION TIME,           
C  AND REQUIRES THE PRESENCE OF QUADRUPLE PRECISION.                    
C                                                                       
            NTMS = 4                 ! 4TH ORDER POLYNOMIAL             
            NTOL = 12                ! 10**-12 ACCURACY                 
            XMAX = 25.0D+00                                             
            CALL CHEBY (MMAX,NTMS,NTOL,XMAX)                            
         END IF                                                         
      END IF                                                            
C                                                                       
C  RECIPROCAL FACTORS USED IN THE DOWNWARDS RECURSION FOR FM(T)         
C                                                                       
      DO M = 1, MXQT                                                    
         RMR(M) = ONE/(2*M-1)                                           
      END DO                                                            
C                                                                       
C  (2M-1)!! FACTORS FOR LARGE-T (IN FM(T))                              
C                                                                       
      DO I = 0, MMAX                                                    
         FII = ONE                                                      
         DO J = 1, 2*I-1, 2                                             
            FII = J*FII                                                 
         END DO                                                         
         TLGM(I) = FII                                                  
      END DO                                                            
C                                                                       
C  SETUP FOR ROTATED AXIS CODE, ORIGINALLY IN GAMGEN                    
C  1) SCREENING PARAMETERS                                              
C                                                                       
      EI1 = EC1*CUTOFF                                                  
      EI2 = EC2*CUTOFF                                                  
      CUX = CX1*RLN10                                                   
C                                                                       
C  2) FIND MAXIMUM VALUE OF S, P, AND/OR D COEFFICIENTS                 
C     FOR SCREENINGS IN THE ROTATED AXIS INTEGRAL CODES.                
C                                                                       
      DO I = 1, NSHELL                                                  
         L = KSTART(I)                                                  
         N = L+KNG(I)-1                                                 
         DO J = L, N                                                    
            CMAX(J)= MAX( ABS(CS(J)), ABS(CP(J)), ABS(CD(J)) )          
         END DO                                                         
      END DO                                                            
C                                                                       
C  3) FOR QMEFP2:FIND MAXIMUM VALUE OF S, P, AND/OR D COEFFICIENTS      
C     FOR SCREENINGS IN THE ROTATED AXIS INTEGRAL CODES.                
C                                                                       
      LVLEFP = LEVELEFP()                                               
      IF(LVLEFP.EQ.2) THEN                                              
      DO MJ = 1, NFRG                                                   
      JM = ISET(MJ)                                                     
      DO I = 1, NSHLEF(JM)                                              
         L = KSTREF(I,JM)                                               
         N = L+KNGEF(I,JM)-1                                            
         DO J = L, N                                                    
            CMAXEF(J,JM)=                                               
     &      MAX( ABS(CSEF(J,JM)), ABS(CPEF(J,JM)),ABS(CDEF(J,JM)) )     
         END DO                                                         
      END DO                                                            
      ENDDO                                                             
      ENDIF                                                             
C                                                                       
      RETURN                                                            
 9900 FORMAT(/20X,22(1H-)/20X,'AO INTEGRAL TECHNOLOGY'/20X,22(1H-)/     
     *  5X,'S,P,L SHELL ROTATED AXIS INTEGRALS, REPROGRAMMED BY'/       
     *  8X,'KAZUYA ISHIMURA (IMS) AND JOSE SIERRA (SYNSTAR).'/          
     *  5X,'S,P,D,L SHELL ROTATED AXIS INTEGRALS PROGRAMMED BY'/        
     *  8X,'KAZUYA ISHIMURA (INSTITUTE FOR MOLECULAR SCIENCE).'/        
     *  5X,'S,P,D,F,G SHELL TO TOTAL QUARTET ANGULAR MOMENTUM SUM 5,'/  
     *  8X,'ERIC PROGRAM BY GRAHAM FLETCHER (ELORET AND NASA ADVANCED'/ 
     *  8X,'SUPERCOMPUTING DIVISION, AMES RESEARCH CENTER).'/           
     *  5X,'S,P,D,F,G,L SHELL GENERAL RYS QUADRATURE PROGRAMMED BY'/    
     *  8X,'MICHEL DUPUIS (PACIFIC NORTHWEST NATIONAL LABORATORY).')    
      END                                                               
C*MODULE INT2C   *DECK CHEBY                                            
C>                                                                      
C>    @brief   incomplete gamma function interpolation setup            
C>                                                                      
C>    @details THIS ROUTINE GENERATES INTERPOLATION GRIDS FOR THE       
C>             INCOMPLETE GAMMA FUNCTION AND EXP(-X) BASED UPON A       
C>             CHEBYSHEV POLYNOMIAL FIT.  IN THIS APPROACH THE GRID     
C>             SPACING IS A FUNCTION OF THE REQUESTED ACCURACY, SO      
C>             HIGHER ORDER (SEE NTMS, BELOW) CORRESPONDS TO FEWER      
C>             GRID POINTS BUT MORE INTERPOLATION FLOPS IN THE INTEGRAL 
C>             CODE, AND VICE VERSA FOR LOWER ORDER. GAMMA FUNCTIONS    
C>             ARE COMPUTED USING THE TAYLOR SERIES FORMULA FOR SMALL   
C>             ARGUMENT. QUAD PRECISION SEEMS ESSENTIAL FOR OBTAINING   
C>             ACCURATE GRID DATA NEAR TO THE LIMIT (25), ESPECIALLY    
C>             FOR HIGH ORDER (SEE MMAX, BELOW).  A TABLE OF RECIPROCALS
C>             IS PRE-COMPUTED FOR SPEED. FOR MMAX HIGHER THAN 16,      
C>             KFX (BELOW) MAY NEED TO BE INCREASED.                    
C>                                                                      
C>    @author  Graham Fletcher, 2004                                    
C>                                                                      
      SUBROUTINE CHEBY (MMAX,NTMS,NTOL,XMAX)                            
C                                                                       
C  INPUT:                                                               
C     MMAX  = HIGHEST ORDER OF GAMMA FUNCTION                           
C     NTMS  = ORDER OF POLYNOMIAL (4TH IS USUALLY SUFFICIENT)           
C     NTOL  = ACCURACY OF INTERPOLATION                                 
C     XMAX  = END OF INTERPOLATION RANGE (25.0)                         
C                                                                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C                                                                       
      PARAMETER (NTX=4)                                                 
      PARAMETER (NPF=450)                                               
      PARAMETER (NGRD=7)                                                
      PARAMETER (NPX=1000)                                              
      PARAMETER (MXQT=16)                                               
      COMMON /FMTTBL/ FGRID(0:NTX,0:NPF,0:NGRD),XGRID(0:NTX,0:NPX)      
     *,               TMAX,RFINC(0:NGRD),RXINC                          
     *,               RMR(MXQT),TLGM(0:MXQT),NORD                       
      LOGICAL                                     GOPARR,DSKWRK,MASWRK  
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
C                                                                       
      LOGICAL WTGRID                                                    
      DOUBLE PRECISION XPW(0:NTX)                                       
C                                                                       
      PARAMETER (KFX=300)                                               
      REAL*16 RMK(KFX),QI                                               
C                                                                       
      PARAMETER (ZER=0.0D+00)                                           
      PARAMETER (ONE=1.0D+00)                                           
      PARAMETER (TWO=2.0D+00)                                           
      PARAMETER (TEN=10.0D+00)                                          
C                                                                       
C  CHEBYSHEV COEFFICIENTS UP TO 12TH ORDER (!)                          
C  GDF 12/09/02: DO NOT KNOW IF 12TH ORDER WORKS, TRIED 8TH ORDER.      
C                                                                       
      INTEGER ICHEB(49)                                                 
      DATA ICHEB/1,1,-1,2,-3,4,1,-8,8,5,-20,16,-1,18,-48,32,-7,56       
     *,         -112,64,1,-32,160,-256,128,9,-120,432,-576,256          
     *,         -1,50,-400,1120,-1280,512,-11,220,-1232,2816,-2816      
     *,          1024,1,-72,840,-3584,6912,-6144,2048/                  
C                                                                       
C  BINOMIAL COEFFICIENTS - ARRANGED TO MATCH CHEBYSHEV TERMS            
C                                                                       
      INTEGER IBINO(252)                                                
      DATA IBINO/1,1,1,1,1,2,1,1,1,1,3,3,1,1,1,2,1,1,4,6,4,1            
     *,          1,1,1,3,3,1,1,5,10,10,5,1,1,1,2,1,1,4,6,4,1            
     *,          1,6,15,20,15,6,1,1,1,1,3,3,1,1,5,10,10,5,1,1           
     *,          7,21,35,35,21,7,1,1,1,2,1,1,4,6,4,1,1,6,15,20          
     *,          15,6,1,1,8,28,56,70,56,28,8,1,1,1,1,3,3,1,1,5          
     *,          10,10,5,1,1,7,21,35,35,21,7,1,1,9,36,84,126,126        
     *,          84,36,9,1,1,1,2,1,1,4,6,4,1,1,6,15,20,15,6,1,1,8       
     *,          28,56,70,56,28,8,1,1,10,45,120,210,252,210,120,45      
     *,          10,1,1,1,1,3,3,1,1,5,10,10,5,1,1,7,21,35,35,21,7,1     
     *,          1,9,36,84,126,126,84,36,9,1,1,11,55,165,330,462,462    
     *,          330,165,55,11,1,1,1,2,1,1,4,6,4,1,1,6,15,20,15,6,1,1   
     *,          8,28,56,70,56,28,8,1,1,10,45,120,210,252,210,120,45    
     *,          10,1,1,12,66,220,495,792,924,792,495,220,66,12,1/      
C                                                                       
      NORD = NTMS                                                       
      TMAX = XMAX                                                       
C                                                                       
      TOL  = TEN**(-NTOL)                                               
      NT2  = 2**NTMS                                                    
      NTF  = 1                                                          
      DO I = 1, NTMS + 1                                                
         NTF = NTF*I                                                    
      END DO                                                            
      EPN  = ONE/( NTMS + 1 )                                           
C                                                                       
C  PRE-COMPUTE RECIPROCALS IN QUAD                                      
C                                                                       
      DO I = 1, KFX                                                     
         QI = I                                                         
         RMK(I) = 1.0Q+00/QI                                            
      END DO                                                            
C                                                                       
C  LOOP OVER SELECTED GAMMA FUNCTION ORDERS, M                          
C  GRIDS ARE COMPUTED FOR M = 0,1,2,3,4, 8, 12, 16, ...                 
C  TO SAVE STORAGE, THE GRIDS FOR 8,12,16 ARE STORED AT FGRID(,,5),     
C  FGRID(,,6),FGRID(,,7) RESPECTIVELY, AND YOU MUST RECUR DOWNWARD      
C  FROM 8,12,16 WASTING 8,7,6 VALUES, IF FOR EXAMPLE YOU WANT M=5.      
C                                                                       
      NGRIDS = 4                                                        
      IF(MMAX.GT.4) NGRIDS = NGRIDS + (MMAX-1)/4                        
      DO IGRD = 0, NGRIDS                                               
         MFT = IGRD                                                     
         IF(IGRD.GT.4) MFT = (IGRD-3)*4                                 
C                                                                       
C  FORMULA FOR OPTIMUM INTERVAL WIDTH                                   
C                                                                       
         DEL  = NT2*NTF*(2*(MFT+NTMS+1)+1)                              
         DEL  = DEL**EPN                                                
         DEL  = DEL*(TOL**EPN)                                          
         FINC = DEL*TWO                                                 
         RFINC(IGRD) = ONE/FINC                                         
         NPTS = NINT(XMAX/FINC)                                         
C                                                                       
C  LOOP OVER GRID POINTS                                                
C                                                                       
         XPT  = ZER                                                     
         DO IPT = 0, NPTS                                               
C                                                                       
C  GENERATE POWERS OF NORMALIZED GRID POINT                             
C                                                                       
            XD = -XPT/DEL                                               
            XPK = ONE                                                   
            DO J = 0, NTMS                                              
               XPW(J) = XPK                                             
               XPK = XPK*XD                                             
               FGRID(J,IPT,IGRD) = ZER                                  
            END DO                                                      
C                                                                       
C  GENERATE TERM COEFFICIENTS                                           
C                                                                       
            IC = 0                                                      
            IB = 0                                                      
            DO KTM = 0, NTMS                                            
               CALL CHEBYG (DEL,XPT,KTM,MFT,RMK,AWT)                    
               NCHEB = (KTM+2)/2    ! TRUNCATE                          
               IN = MOD(KTM,2)                                          
               DO I = 1, NCHEB                                          
                  IC = IC + 1                                           
                  CHEB = ICHEB(IC)*AWT                                  
                  DO J = 0, IN                                          
                     IB = IB + 1                                        
                     II = IN - J                                        
      FGRID(II,IPT,IGRD) = FGRID(II,IPT,IGRD) +                         
     *                    IBINO(IB)*CHEB*XPW(J)*(2**II)                 
                  END DO                                                
                  IN = IN + 2                                           
               END DO                                                   
            END DO                                                      
            XPT = XPT + FINC                                            
         END DO                                                         
      END DO                                                            
C                                                                       
C  INTERPOLATION GRID FOR EXP(-X)                                       
C  FORMULA FOR OPTIMUM INTERVAL WIDTH                                   
C                                                                       
      DEL  = NT2*NTF                                                    
      DEL  = DEL**EPN                                                   
      DEL  = DEL*(TOL**EPN)                                             
      XINC = DEL*TWO                                                    
      RXINC= ONE/XINC                                                   
      NPTS = NINT(XMAX/XINC)                                            
C                                                                       
C  LOOP OVER GRID POINTS                                                
C                                                                       
      XPT  = ZER                                                        
      DO IPT = 0, NPTS                                                  
C                                                                       
C  GENERATE POWERS OF NORMALIZED GRID POINT                             
C                                                                       
         XD = -XPT/DEL                                                  
         XPK = ONE                                                      
         DO J = 0, NTMS                                                 
            XPW(J) = XPK                                                
            XPK = XPK*XD                                                
            XGRID(J,IPT) = ZER                                          
         END DO                                                         
C                                                                       
C  GENERATE TERM COEFFICIENTS                                           
C                                                                       
         IC = 0                                                         
         IB = 0                                                         
         DO KTM = 0, NTMS                                               
            CALL CHEBYX (DEL,XPT,KTM,AWT)                               
            NCHEB = (KTM+2)/2    ! TRUNCATE                             
            IN = MOD(KTM,2)                                             
            DO I = 1, NCHEB                                             
               IC = IC + 1                                              
               CHEB = ICHEB(IC)*AWT                                     
               DO J = 0, IN                                             
                  IB = IB + 1                                           
                  II = IN - J                                           
      XGRID(II,IPT) = XGRID(II,IPT) +                                   
     *             IBINO(IB)*CHEB*XPW(J)*(2**II)                        
               END DO                                                   
               IN = IN + 2                                              
            END DO                                                      
         END DO                                                         
         XPT = XPT + XINC                                               
      END DO                                                            
C                                                                       
C        OPTION TO SAVE THE FM(T) INTERPOLATION DATA JUST GENERATED,    
C        FLIP VARIABLES RDGRID/WTGRID IN ORDER TO REGENERATE THE DATA.  
C                                                                       
C        MACHINES LACKING QUADRUPLE PRECISION MUST READ THIS DATA FROM  
C        DISC, AS THEY CANNOT GENERATE IT SUFFICIENTLY ACCURATELY.      
C                                                                       
      WTGRID = .FALSE.                                                  
      IF(.NOT.WTGRID) RETURN                                            
      IF(.NOT.MASWRK) RETURN                                            
C                                                                       
      OPEN(UNIT=2, FILE='/U1/MIKE/GAMESS/ERICFMT.DAT', STATUS='UNKNOWN',
     *     FORM='FORMATTED', ACCESS='SEQUENTIAL')                       
      WRITE(2,9000)                                                     
      WRITE(2,9010)                                                     
      WRITE(2,9100) (RFINC(III),III=0,NGRD)                             
      WRITE(2,9020)                                                     
      DO KKK=0,NGRD                                                     
         DO JJJ=0,NPF                                                   
            WRITE(2,9100) (FGRID(III,JJJ,KKK),III=0,NTX)                
         ENDDO                                                          
      ENDDO                                                             
      WRITE(2,9030)                                                     
      WRITE(2,9100) RXINC,TMAX                                          
      WRITE(2,9040)                                                     
      DO JJJ=0,NPX                                                      
         WRITE(2,9100) (XGRID(III,JJJ),III=0,NTX)                       
      ENDDO                                                             
      CLOSE(UNIT=2, STATUS='KEEP')                                      
C                                                                       
      RETURN                                                            
 9000 FORMAT('DATA FOR ERIC AND ROTATED AXIS INTEGRAL CODE''S GAMMA',   
     *       ' FUNCTION INTERPOLATION')                                 
 9010 FORMAT('FGRID INCREMENT VALUES')                                  
 9020 FORMAT('FGRID TABLE')                                             
 9030 FORMAT('XGRID INCREMENT, TMAX VALUE')                             
 9040 FORMAT('XGRID TABLE')                                             
 9100 FORMAT(1P,3E25.16)                                                
      END                                                               
C*MODULE INT2C   *DECK CHEBYG                                           
C>                                                                      
C>    @brief   CHEBYSHEV WEIGHT FOR GAMMA FUNCTION                      
C>                                                                      
C>    @details CHEBYSHEV WEIGHT FOR GAMMA FUNCTION                      
C>                                                                      
C>    @author  Graham Fletcher, 2004                                    
C>                                                                      
      SUBROUTINE CHEBYG (DEL,XPT,KTM,MFT,RMK,AWT)                       
C                                                                       
C  CHEBYSHEV WEIGHT FOR GAMMA FUNCTION                                  
C                                                                       
C  DEL = HALF INTERVAL WIDTH     [INPUT]                                
C  XPT = GRID POINT              [INPUT]                                
C  KTM = POLYNOMIAL TERM         [INPUT]                                
C  MFT = ORDER OF GAMMA FUNCTION [INPUT]                                
C  RMK = TABLE OF RECIPROCALS    [INPUT]                                
C        QUAD NOT ESSENTIAL HERE BUT HEY                                
C  AWT = TERM WEIGHT             [OUTPUT]                               
C                                                                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C                                                                       
      PARAMETER (KFX=300)                                               
      REAL*16 RMK(KFX)                                                  
C                                                                       
      LOGICAL CONVERGED                                                 
C                                                                       
      PARAMETER (ZER=0.0D+00)                                           
      PARAMETER (PT5=0.5D+00)                                           
      PARAMETER (ONE=1.0D+00)                                           
      PARAMETER (TWO=2.0D+00)                                           
      PARAMETER (TOL=1.0D-12)                                           
C                                                                       
      DD2 = DEL*PT5                                                     
      DSQ = DD2*DD2                                                     
      AWT = ZER                                                         
      AM  = ONE                                                         
      FK  = ONE                                                         
      DO M = 1, KTM                                                     
         FK = FK*RMK(M)                                                 
      END DO                                                            
C                                                                       
      CONVERGED = .FALSE.                                               
      M   = 0                                                           
      FM  = ONE                                                         
      DO WHILE (.NOT.CONVERGED)                                         
         AP  = AWT                                                      
         MMK = KTM + MFT + 2*M                                          
         CALL FTMVAL (XPT,MMK,RMK,FTM)                                  
         AWT = AWT + AM*FTM*FM*FK                                       
         CONVERGED = ABS(AP-AWT).LT.TOL                                 
         M = M + 1                                                      
         FM = FM*RMK(M)                                                 
         FK = FK*RMK(KTM+M)                                             
         AM = AM*DSQ                                                    
      END DO                                                            
      AWT = AWT*(DD2**KTM)                                              
      IF(KTM.GT.0) AWT = AWT*TWO                                        
      IF(MOD(KTM,2).NE.0) AWT = -AWT                                    
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2C   *DECK FTMVAL                                           
C>                                                                      
C>    @brief   evaluate Fm(t) table                                     
C>                                                                      
C>    @details FORMULA FOR COMPUTING THE GAMMA FUNCTION USING TAYLOR    
C>             EXPANSION FOR TT<25. CARE TAKEN TO DO THIS IN A          
C>             NUMERICALLY STABLE WAY.  NOTE THAT QUADRUPLE PRECISION   
C>             IS USED, ALTHOUGH DOUBLE PRECISION ARGUMENTS ARE PASSED. 
C>             QUADRUPLE PRECISION IS CRITICAL!                         
C>             Because Q.P. is frequently unavailable, this routine     
C>             is not called, instead the data it generates is read     
C>             from a disk file prepared by a machine that has Q.P.     
C>                                                                      
C>    @author  Graham Fletcher, 2004                                    
C>                                                                      
      SUBROUTINE FTMVAL (TT8,MFT,RMK,FTM8)                              
C                                                                       
C  FORMULA FOR COMPUTING THE GAMMA FUNCTION                             
C  USING TAYLOR EXPANSION FOR TT<25. CARE TAKEN                         
C  TO DO THIS IN A NUMERICALLY STABLE WAY.                              
C  NOTE THAT QUADRUPLE PRECISION IS USED, ALTHOUGH DOUBLE PRECISION     
C  ARGUMENTS ARE PASSED. QUADRUPLE PRECISION IS CRITICAL!               
C                                                                       
C  TT8  = GAMMA FUNCTION ARGUMENT  [INPUT]                              
C  MFT  = ORDER OF GAMMA FUNCTION  [INPUT]                              
C  RMK  = TABLE OF RECIPROCALS     [INPUT]                              
C  FTM8 = VALUE OF GAMMA FUNCTION [OUTPUT]                              
C                                                                       
      IMPLICIT REAL*16 (A-H,O-Z)                                        
C                                                                       
      PARAMETER (KFX=300)                                               
      REAL*16 RMK(KFX)                                                  
C                                                                       
      DOUBLE PRECISION TT8,FTM8                                         
C                                                                       
      LOGICAL CONVERGED                                                 
C                                                                       
      PARAMETER (ONE=1.0Q+00)                                           
      PARAMETER (TOL=1.0Q-17)                                           
C                                                                       
      TT = TT8                                                          
C                                                                       
      CONVERGED = .FALSE.                                               
      L   = 1                                                           
      M   = 2*MFT + 1                                                   
      FTM = RMK(M)                                                      
      XK  = ONE                                                         
      DO WHILE (.NOT.CONVERGED)                                         
         M = M + 2                                                      
         A = -TT*RMK(M)                                                 
         FP = FTM                                                       
         FTM = FTM + A*XK                                               
         CONVERGED = ABS(FP-FTM).LT.TOL                                 
         L = L + 1                                                      
         XK = -XK*TT*RMK(L)                                             
      END DO                                                            
C                                                                       
      FTM8 = FTM                                                        
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2C   *DECK CHEBYX                                           
C>                                                                      
C>    @brief   CHEBYSHEV WEIGHT FOR EXP(-X)                             
C>                                                                      
C>    @details CHEBYSHEV WEIGHT FOR EXP(-X)                             
C>                                                                      
C>    @author  Graham Fletcher, 2004                                    
C>                                                                      
      SUBROUTINE CHEBYX (DEL,XPT,KTM,AWT)                               
C                                                                       
C  CHEBYSHEV WEIGHT FOR EXP(-X)                                         
C                                                                       
C  DEL = HALF INTERVAL WIDTH  [INPUT]                                   
C  XPT = GRID POINT           [INPUT]                                   
C  KTM = POLYNOMIAL TERM      [INPUT]                                   
C  AWT = TERM WEIGHT         [OUTPUT]                                   
C                                                                       
C  NOTE THAT APART FROM THE LAST LINE THE WEIGHTS                       
C  ARE INDEPENDENT OF THE VALUE OF XPT                                  
C                                                                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C                                                                       
      LOGICAL CONVERGED                                                 
C                                                                       
      PARAMETER (ZER=0.0D+00)                                           
      PARAMETER (PT5=0.5D+00)                                           
      PARAMETER (ONE=1.0D+00)                                           
      PARAMETER (TWO=2.0D+00)                                           
      PARAMETER (TOL=1.0D-12)                                           
C                                                                       
      DD2 = DEL*PT5                                                     
      DSQ = DD2*DD2                                                     
      AWT = ZER                                                         
      AM  = ONE                                                         
      FK  = ONE                                                         
      DO M = 1, KTM                                                     
         FK = FK*M                                                      
      END DO                                                            
C                                                                       
      CONVERGED = .FALSE.                                               
      M   = 0                                                           
      FM  = ONE                                                         
      DO WHILE (.NOT.CONVERGED)                                         
         AP  = AWT                                                      
         AWT = AWT + AM/(FM*FK)                                         
         CONVERGED = ABS(AP-AWT).LT.TOL                                 
         M = M + 1                                                      
         FM = FM*M                                                      
         FK = FK*(KTM+M)                                                
         AM = AM*DSQ                                                    
      END DO                                                            
      AWT = AWT*(DD2**KTM)                                              
      IF(KTM.GT.0) AWT = AWT*TWO                                        
      IF(MOD(KTM,2).NE.0) AWT = -AWT                                    
      AWT = AWT*EXP(-XPT)                                               
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2C   *DECK GRIDIN                                           
C>                                                                      
C>    @brief   read Fm(t) table                                         
C>                                                                      
C>    @details read Fm(t) table used by ERIC and rotated axis codes     
C>                                                                      
C>    @author  Graham Fletcher, 2004                                    
C>                                                                      
      SUBROUTINE GRIDIN                                                 
C                                                                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C                                                                       
      PARAMETER (NTX=4)                                                 
      PARAMETER (NPF=450)                                               
      PARAMETER (NGRD=7)                                                
      PARAMETER (NPX=1000)                                              
      PARAMETER (MXQT=16)                                               
      COMMON /FMTTBL/ FGRID(0:NTX,0:NPF,0:NGRD),XGRID(0:NTX,0:NPX)      
     *,               TMAX,RFINC(0:NGRD),RXINC                          
     *,               RMR(MXQT),TLGM(0:MXQT),NORD                       
      LOGICAL         ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI             
      COMMON /GDDI/   ISCOPE,NGROUPS,MYGROUP,MEGLOB,npglob,nnglob,JBTYP,
     *                ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,NSUBGR,     
     *                MeUniv,NPUniv,numdlb,myworld,nworlds              
      LOGICAL                                     GOPARR,DSKWRK,MASWRK  
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK  
C                                                                       
      LOGICAL SVDSKW                                                    
C                                                                       
C         FORMATTED READ OF VARIOUS GAMMA FUNCTION INTERPOLATION TABLES.
C         DATA MUST BE GENERATED ON A MACHINE WITH QUADRUPLE PRECISION. 
C                                                                       
      SVDSKW = DSKWRK                                                   
      DSKWRK = .FALSE.                                                  
      LU=2                                                              
      if(.not.isgddi .or. maswrk)                                       
     *   CALL SEQOPN(LU, 'ERICFMT', 'OLD', .TRUE., 'FORMATTED')         
C                                                                       
      IF(MASWRK) THEN                                                   
         READ(LU,9000)                                                  
         READ(LU,9000)                                                  
         READ(LU,9100) (RFINC(III),III=0,NGRD)                          
         READ(LU,9000)                                                  
         DO KKK=0,NGRD                                                  
            DO JJJ=0,NPF                                                
               READ(LU,9100) (FGRID(III,JJJ,KKK),III=0,NTX)             
            ENDDO                                                       
         ENDDO                                                          
         READ(LU,9000)                                                  
         READ(LU,9100) RXINC,TMAX                                       
         READ(LU,9000)                                                  
         DO JJJ=0,NPX                                                   
            READ(LU,9100) (XGRID(III,JJJ),III=0,NTX)                    
         ENDDO                                                          
      END IF                                                            
C                                                                       
      if(.not.isgddi .or. maswrk)                                       
     *   CALL SEQCLO(LU,'KEEP')                                         
      DSKWRK = SVDSKW                                                   
C                                                                       
C         PARALLEL MCSCF RUNS USING DUPLICATED AO INTEGRALS MAY CALL    
C         THE INTEGRAL PACKAGE WITH GOPARR TEMPORARILY SET .FALSE.,     
C         BUT WE MUST BROADCAST THE TABLE INFORMATION ANYWAY.           
C                                                                       
      IF(GOPARR .OR. NPROC.GT.1) THEN                                   
         CALL DDI_BCAST(241,'F',RFINC, NGRD+1                 ,MASTER)  
         CALL DDI_BCAST(242,'F',FGRID,(NGRD+1)*(NPF+1)*(NTX+1),MASTER)  
         CALL DDI_BCAST(243,'F',RXINC,1                       ,MASTER)  
         CALL DDI_BCAST(244,'F',TMAX ,1                       ,MASTER)  
         CALL DDI_BCAST(245,'F',XGRID,         (NPX+1)*(NTX+1),MASTER)  
      END IF                                                            
C                                                                       
      RETURN                                                            
 9000 FORMAT(1X)                                                        
 9100 FORMAT(3E25.16)                                                   
      END                                                               
C                                                                       
C  START OF FAST ROUTINES                                               
C  SP CASES                                                             
C                                                                       
C  LPHI=   1                                                            
C  LWK1=   1                                                            
C  LWK2=   1                                                            
C  LENW=   1                                                            
C*MODULE INT2C   *DECK SSSS                                             
C>                                                                      
C>    @brief   ERIC ssss case                                           
C>                                                                      
C>    @details ERIC [ss|ss] integral quartet                            
C>                                                                      
C>    @author  Graham Fletcher, 2004, modified Jose Sierra, 2013.       
C>                                                                      
      SUBROUTINE SSSS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL,PHI,WK2,LENW)  
      USE lrcdft, ONLY: EMU2                                            
      use mx_limits, only: mxgsh,mxg2                                   
C                                                                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C                                                                       
      LOGICAL    IEQJ,KEQL                                              
      DIMENSION  PHI(*),WK2(LENW,*)                                     
C                                                                       
C                                                                       
      PARAMETER (NTX=4)                                                 
      PARAMETER (NPF=450)                                               
      PARAMETER (NGRD=7)                                                
      PARAMETER (NPX=1000)                                              
      PARAMETER (MXQT=16)                                               
      COMMON /FMTTBL/ FGRID(0:NTX,0:NPF,0:NGRD),XGRID(0:NTX,0:NPX)      
     *,               TMAX,RFINC(0:NGRD),RXINC                          
     *,               RMR(MXQT),TLGM(0:MXQT),NORD                       
      COMMON /ERIPRM/ EXI(MXGSH),EXJ(MXGSH),EXK(MXGSH),EXL(MXGSH),      
     *                CCI(MXGSH),CCJ(MXGSH),CCK(MXGSH),CCL(MXGSH),      
     *                XAB(MXG2),YAB(MXG2),ZAB(MXG2),                    
     *                XCD(MXG2),YCD(MXG2),ZCD(MXG2),                    
     *                CCBRA(MXG2),CCKET(MXG2),RXB(MXG2),                
     *                SLBRA(MXG2),SLKET(MXG2),RXK(MXG2)                 
C$omp threadprivate(/ERIPRM/)
      LOGICAL         LRINT                                             
      COMMON /NLRCF / LRINT                                             
C$omp threadprivate(/NLRCF /)
      COMMON /SHLNOS/ QQ4,IDUMMY(20)                                    
C$omp threadprivate(/SHLNOS/)
C                                                                       
      DIMENSION  FT(0:16)                                               
C                                                                       
      PARAMETER (ZER=0.0D+00)                                           
      PARAMETER (ONE=1.0D+00)                                           
      PARAMETER (CCTOL=1.0D-13)                                         
C                                                                       
      PHI(  1)= ZER                                                     
      KL=0                                                              
      DO LL=1,LPRIM                                                     
         X04= EXL(LL)                                                   
         KTOP=KPRIM                                                     
         IF(KEQL) KTOP=LL                                               
         DO KK=1,KTOP                                                   
            KL=KL+1                                                     
            CFK= CCKET(KL)                                              
            IF(ABS(CFK).GT.CCTOL) THEN                                  
CC             CFK= CFK*QQ4                                             
               X03= EXK(KK)                                             
               X34= X03+X04                                             
               XKL= XCD(KL)                                             
               YKL= YCD(KL)                                             
               ZKL= ZCD(KL)                                             
               IJ=0                                                     
               DO JJ=1,JPRIM                                            
                  X02= EXJ(JJ)                                          
                  ITOP=IPRIM                                            
                  IF(IEQJ) ITOP=JJ                                      
                  DO II=1,ITOP                                          
                     IJ=IJ+1                                            
                     CFB= CCBRA(IJ)                                     
                     IF(ABS(CFB).GT.CCTOL) THEN                         
                        X01= EXI(II)                                    
                        X12= X01+X02                                    
                        X41= ONE/(X12+X34)                              
                        FTZ= CFB*CFK                                    
                        RX = XKL-XAB(IJ)                                
                        RY = YKL-YAB(IJ)                                
                        RZ = ZKL-ZAB(IJ)                                
                        RSQ= RX*RX+RY*RY+RZ*RZ                          
                        RHO= X12*X34*X41                                
                        IF(LRINT) THEN                                  
                           EFR= EMU2/(EMU2+RHO)                         
                           RHO= RHO*EFR                                 
                        ENDIF                                           
                        TT = RSQ*RHO                                    
                        N=0                                             
                        IF(TT.LE.TMAX) THEN                             
C                                                                       
C     FM(T) EVALUATION                                                  
C                                                                       
                           TV= TT*RFINC(N)                              
                           IP= NINT(TV)                                 
                           FX=    FGRID(4,IP,N) *TV                     
                           FX=(FX+FGRID(3,IP,N))*TV                     
                           FX=(FX+FGRID(2,IP,N))*TV                     
                           FX=(FX+FGRID(1,IP,N))*TV                     
                           FX= FX+FGRID(0,IP,N)                         
C                                                                       
                           FT(N)= FX                                    
CC                         T2= TT+TT                                    
CC                         DO M=N,1,-1                                  
CC                            FT(M-1)=(T2*FT(M)+ET)*RMR(M)              
CC                         END DO                                       
C                                                                       
CC                         RHO= RHO+RHO                                 
                           IF(LRINT) FTZ= FTZ*SQRT(EFR)                 
                           FTF= FTZ*SQRT(X41)                           
                              FT(0)= FT(0)*FTF                          
CC                         DO 210 M=0,N                                 
CC                            FT(M)= FT(M)*FTF                          
CC210                      FTF= FTF*RHO                                 
                        ELSE                                            
                           XIN= ONE/RSQ                                 
                           FTF= FTZ*SQRT(XIN)*SLBRA(IJ)*SLKET(KL)       
                           FT(0)= FTF                                   
CC                         DO 220 M=1,N                                 
CC                            FTF= FTF*XIN                              
CC220                      FT(M)= TLGM(M)*FTF                           
                        END IF                                          
                        PHI(  1)= PHI(  1)+FT(0)                        
                     END IF                                             
                  END DO                                                
               END DO                                                   
            END IF                                                      
         END DO                                                         
      END DO                                                            
C  POST-CONTRACTION PHASE                                               
      WK2(  1,1)= PHI(  1)*QQ4                                          
C                                                                       
      RETURN                                                            
      END                                                               
C  LPHI=  21                                                            
C  LWK1=   1                                                            
C  LWK2=   4                                                            
C  LENW=   4                                                            
C*MODULE INT2C   *DECK PSSS                                             
C>                                                                      
C>    @brief   ERIC psss case                                           
C>                                                                      
C>    @details ERIC [ps|ss] integral quartet                            
C>                                                                      
C>    @author  Graham Fletcher, 2004, modified Jose Sierra, 2013.       
C>                                                                      
      SUBROUTINE PSSS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL,               
     *                 XC,YC,ZC,XD,YD,ZD,PHI,WK2,LENW)                  
      USE lrcdft, ONLY: EMU2                                            
      use mx_limits, only: mxgsh,mxg2                                   
C                                                                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C                                                                       
      LOGICAL    IEQJ,KEQL                                              
      DIMENSION  PHI(*),WK2(LENW,*)                                     
C                                                                       
C                                                                       
      PARAMETER (NTX=4)                                                 
      PARAMETER (NPF=450)                                               
      PARAMETER (NGRD=7)                                                
      PARAMETER (NPX=1000)                                              
      PARAMETER (MXQT=16)                                               
      COMMON /FMTTBL/ FGRID(0:NTX,0:NPF,0:NGRD),XGRID(0:NTX,0:NPX)      
     *,               TMAX,RFINC(0:NGRD),RXINC                          
     *,               RMR(MXQT),TLGM(0:MXQT),NORD                       
      COMMON /ERIPRM/ EXI(MXGSH),EXJ(MXGSH),EXK(MXGSH),EXL(MXGSH),      
     *                CCI(MXGSH),CCJ(MXGSH),CCK(MXGSH),CCL(MXGSH),      
     *                XAB(MXG2),YAB(MXG2),ZAB(MXG2),                    
     *                XCD(MXG2),YCD(MXG2),ZCD(MXG2),                    
     *                CCBRA(MXG2),CCKET(MXG2),RXB(MXG2),                
     *                SLBRA(MXG2),SLKET(MXG2),RXK(MXG2)                 
C$omp threadprivate(/ERIPRM/)
      LOGICAL         LRINT                                             
      COMMON /NLRCF / LRINT                                             
C$omp threadprivate(/NLRCF /)
      COMMON /SHLNOS/ QQ4,IDUMMY(20)                                    
C$omp threadprivate(/SHLNOS/)
C                                                                       
      DIMENSION  FT(0:16),SFAC( 9),CNF(11)                              
C                                                                       
      PARAMETER (ZER=0.0D+00)                                           
      PARAMETER (ONE=1.0D+00)                                           
      PARAMETER (TWO=2.0D+00)                                           
      PARAMETER (CCTOL=1.0D-13)                                         
      DATA CNF/0.0D+00,0.0D+00,0.0D+00,1.0D+00,2.0D+00,3.0D+00,         
     &         4.0D+00,5.0D+00,6.0D+00,7.0D+00,8.0D+00/                 
      SAVE CNF                                                          
!$omp threadprivate(cnf)                                                
C                                                                       
      SFAC(1)= ONE                                                      
      PHI(  8)= ZER                                                     
      PHI( 18)= ZER                                                     
      PHI( 19)= ZER                                                     
      PHI( 20)= ZER                                                     
      KL=0                                                              
      DO LL=1,LPRIM                                                     
         X04= EXL(LL)                                                   
         PHI(  7)= ZER                                                  
         PHI( 15)= ZER                                                  
         PHI( 16)= ZER                                                  
         PHI( 17)= ZER                                                  
         KTOP=KPRIM                                                     
         IF(KEQL) KTOP=LL                                               
         DO KK=1,KTOP                                                   
            KL=KL+1                                                     
            CFK= CCKET(KL)                                              
            IF(ABS(CFK).GT.CCTOL) THEN                                  
               CFK= CFK*QQ4                                             
               X03= EXK(KK)                                             
               X34= X03+X04                                             
               XKL= XCD(KL)                                             
               YKL= YCD(KL)                                             
               ZKL= ZCD(KL)                                             
               PHI(  6)= ZER                                            
               PHI( 12)= ZER                                            
               PHI( 13)= ZER                                            
               PHI( 14)= ZER                                            
               IJ=0                                                     
               DO JJ=1,JPRIM                                            
                  X02= EXJ(JJ)                                          
                  PHI(  5)= ZER                                         
                  PHI(  9)= ZER                                         
                  PHI( 10)= ZER                                         
                  PHI( 11)= ZER                                         
                  ITOP=IPRIM                                            
                  IF(IEQJ) ITOP=JJ                                      
                  DO II=1,ITOP                                          
                     IJ=IJ+1                                            
                     CFB= CCBRA(IJ)                                     
                     IF(ABS(CFB).GT.CCTOL) THEN                         
                        X01= EXI(II)                                    
                        X12= X01+X02                                    
                        X41= ONE/(X12+X34)                              
                        FTZ= CFB*CFK                                    
                        RX = XKL-XAB(IJ)                                
                        RY = YKL-YAB(IJ)                                
                        RZ = ZKL-ZAB(IJ)                                
                        RSQ= RX*RX+RY*RY+RZ*RZ                          
                        RHO= X12*X34*X41                                
                        IF(LRINT) THEN                                  
                           EFR= EMU2/(EMU2+RHO)                         
                           RHO= RHO*EFR                                 
                        ENDIF                                           
                        TT = RSQ*RHO                                    
                        N=1                                             
                        IF(TT.LE.TMAX) THEN                             
C                                                                       
C     FM(T) EVALUATION                                                  
C                                                                       
                           TV= TT*RFINC(N)                              
                           IP= NINT(TV)                                 
                           FX=    FGRID(4,IP,N) *TV                     
                           FX=(FX+FGRID(3,IP,N))*TV                     
                           FX=(FX+FGRID(2,IP,N))*TV                     
                           FX=(FX+FGRID(1,IP,N))*TV                     
                           FX= FX+FGRID(0,IP,N)                         
                           TV= TT*RXINC                                 
                           IP= NINT(TV)                                 
                           ET=    XGRID(4,IP) *TV                       
                           ET=(ET+XGRID(3,IP))*TV                       
                           ET=(ET+XGRID(2,IP))*TV                       
                           ET=(ET+XGRID(1,IP))*TV                       
                           ET= ET+XGRID(0,IP)                           
C                                                                       
                           FT(N)= FX                                    
                           T2= TT+TT                                    
                              FT(1-1)=(T2*FT(1)+ET)*RMR(1)              
CC                         DO M=N,1,-1                                  
CC                            FT(M-1)=(T2*FT(M)+ET)*RMR(M)              
CC                         END DO                                       
C                                                                       
                           RHO= RHO+RHO                                 
                           IF(LRINT) FTZ= FTZ*SQRT(EFR)                 
                           FTF= FTZ*SQRT(X41)                           
                              FT(0)= FT(0)*FTF                          
                           FTF= FTF*RHO                                 
                              FT(1)= FT(1)*FTF                          
CC                         DO 210 M=0,N                                 
CC                            FT(M)= FT(M)*FTF                          
CC210                      FTF= FTF*RHO                                 
                        ELSE                                            
                           XIN= ONE/RSQ                                 
                           FTF= FTZ*SQRT(XIN)*SLBRA(IJ)*SLKET(KL)       
                           FT(0)= FTF                                   
                              FTF= FTF*XIN                              
                           FT(1)= TLGM(1)*FTF                           
CC                         DO 220 M=1,N                                 
CC                            FTF= FTF*XIN                              
CC220                      FT(M)= TLGM(M)*FTF                           
                        END IF                                          
C                                                                       
                        CALL PHIFTS(2,RX,RY,RZ,PHI,FT)                  
C                                                                       
C                       FAC= RXB(IJ)                                    
                        PHI(  5)= PHI(  5)+PHI(  1)                     
C                                                                       
                        PHI(  9)= PHI(  9)+PHI(  2)                     
                        PHI( 10)= PHI( 10)+PHI(  3)                     
                        PHI( 11)= PHI( 11)+PHI(  4)                     
                     END IF                                             
                  END DO                                                
C                 FAC= EXJ(JJ)*TWO                                      
                  PHI(  6)= PHI(  6)+PHI(  5)                           
C                                                                       
                  PHI( 12)= PHI( 12)+PHI(  9)                           
                  PHI( 13)= PHI( 13)+PHI( 10)                           
                  PHI( 14)= PHI( 14)+PHI( 11)                           
               END DO                                                   
               FAC= RXK(KL)                                             
               SFAC(2)= FAC                                             
               PHI(  7)= PHI(  7)+PHI(  6)*SFAC(2)                      
C                                                                       
               PHI( 15)= PHI( 15)+PHI( 12)*SFAC(2)                      
               PHI( 16)= PHI( 16)+PHI( 13)*SFAC(2)                      
               PHI( 17)= PHI( 17)+PHI( 14)*SFAC(2)                      
            END IF                                                      
         END DO                                                         
         FAC= EXL(LL)*TWO                                               
         SFAC(2)= FAC                                                   
         PHI(  8)= PHI(  8)+PHI(  7)*SFAC(2)                            
C                                                                       
         PHI( 18)= PHI( 18)+PHI( 15)                                    
         PHI( 19)= PHI( 19)+PHI( 16)                                    
         PHI( 20)= PHI( 20)+PHI( 17)                                    
      END DO                                                            
C  POST-CONTRACTION PHASE                                               
      WK2(  1,1)= PHI(  8)                                              
      WK2(  2,1)=-PHI( 18)                                              
      WK2(  3,1)=-PHI( 19)                                              
      WK2(  4,1)=-PHI( 20)                                              
      CNF(1)= XD-XC                                                     
      CNF(2)= YD-YC                                                     
      CNF(3)= ZD-ZC                                                     
      DO I=1,1                                                          
         WK2(  2,I)= WK2(  2,I)+WK2(  1,I)*CNF( 1)                      
         WK2(  3,I)= WK2(  3,I)+WK2(  1,I)*CNF( 2)                      
         WK2(  4,I)= WK2(  4,I)+WK2(  1,I)*CNF( 3)                      
      END DO                                                            
C                                                                       
      RETURN                                                            
      END                                                               
C  LPHI=  65                                                            
C  LWK1=   4                                                            
C  LWK2=  12                                                            
C  LENW=   4                                                            
C*MODULE INT2C   *DECK PSPS                                             
C>                                                                      
C>    @brief   ERIC psps case                                           
C>                                                                      
C>    @details ERIC [ps|ps] integral quartet                            
C>                                                                      
C>    @author  Graham Fletcher, 2004, modified Jose Sierra, 2013.       
C>                                                                      
      SUBROUTINE PSPS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL,               
     *                 XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD,             
     *                 PHI,WK1,WK2,LENW)                                
      USE lrcdft, ONLY: EMU2                                            
      use mx_limits, only: mxgsh,mxg2                                   
C                                                                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C                                                                       
      LOGICAL    IEQJ,KEQL                                              
      DIMENSION  PHI(*),WK1(*),WK2(LENW,*)                              
C                                                                       
C                                                                       
      PARAMETER (NTX=4)                                                 
      PARAMETER (NPF=450)                                               
      PARAMETER (NGRD=7)                                                
      PARAMETER (NPX=1000)                                              
      PARAMETER (MXQT=16)                                               
      COMMON /FMTTBL/ FGRID(0:NTX,0:NPF,0:NGRD),XGRID(0:NTX,0:NPX)      
     *,               TMAX,RFINC(0:NGRD),RXINC                          
     *,               RMR(MXQT),TLGM(0:MXQT),NORD                       
      COMMON /ERIPRM/ EXI(MXGSH),EXJ(MXGSH),EXK(MXGSH),EXL(MXGSH),      
     *                CCI(MXGSH),CCJ(MXGSH),CCK(MXGSH),CCL(MXGSH),      
     *                XAB(MXG2),YAB(MXG2),ZAB(MXG2),                    
     *                XCD(MXG2),YCD(MXG2),ZCD(MXG2),                    
     *                CCBRA(MXG2),CCKET(MXG2),RXB(MXG2),                
     *                SLBRA(MXG2),SLKET(MXG2),RXK(MXG2)                 
C$omp threadprivate(/ERIPRM/)
      LOGICAL         LRINT                                             
      COMMON /NLRCF / LRINT                                             
C$omp threadprivate(/NLRCF /)
      COMMON /SHLNOS/ QQ4,IDUMMY(20)                                    
C$omp threadprivate(/SHLNOS/)
C                                                                       
      DIMENSION  FT(0:16),SFAC( 9),CNF(11)                              
C                                                                       
      PARAMETER (ZER=0.0D+00)                                           
      PARAMETER (ONE=1.0D+00)                                           
      PARAMETER (TWO=2.0D+00)                                           
      PARAMETER (CCTOL=1.0D-13)                                         
      DATA CNF/0.0D+00,0.0D+00,0.0D+00,1.0D+00,2.0D+00,3.0D+00,         
     &         4.0D+00,5.0D+00,6.0D+00,7.0D+00,8.0D+00/                 
      SAVE CNF                                                          
!$omp threadprivate(cnf)                                                
C                                                                       
      SFAC(1)= ONE                                                      
      PHI( 15)= ZER                                                     
      PHI( 25)= ZER                                                     
      PHI( 26)= ZER                                                     
      PHI( 27)= ZER                                                     
      PHI( 34)= ZER                                                     
      PHI( 35)= ZER                                                     
      PHI( 36)= ZER                                                     
      DO I= 58, 64                                                      
         PHI(I)= ZER                                                    
      END DO                                                            
      KL=0                                                              
      DO LL=1,LPRIM                                                     
         X04= EXL(LL)                                                   
         PHI( 14)= ZER                                                  
         PHI( 22)= ZER                                                  
         PHI( 23)= ZER                                                  
         PHI( 24)= ZER                                                  
         PHI( 31)= ZER                                                  
         PHI( 32)= ZER                                                  
         PHI( 33)= ZER                                                  
         DO I= 51, 57                                                   
            PHI(I)= ZER                                                 
         END DO                                                         
         KTOP=KPRIM                                                     
         IF(KEQL) KTOP=LL                                               
         DO KK=1,KTOP                                                   
            KL=KL+1                                                     
            CFK= CCKET(KL)                                              
            IF(ABS(CFK).GT.CCTOL) THEN                                  
               CFK= CFK*QQ4                                             
               X03= EXK(KK)                                             
               X34= X03+X04                                             
               XKL= XCD(KL)                                             
               YKL= YCD(KL)                                             
               ZKL= ZCD(KL)                                             
               PHI( 13)= ZER                                            
               PHI( 19)= ZER                                            
               PHI( 20)= ZER                                            
               PHI( 21)= ZER                                            
               PHI( 28)= ZER                                            
               PHI( 29)= ZER                                            
               PHI( 30)= ZER                                            
               DO I= 44, 50                                             
                  PHI(I)= ZER                                           
               END DO                                                   
               IJ=0                                                     
               DO JJ=1,JPRIM                                            
                  X02= EXJ(JJ)                                          
                  PHI( 12)= ZER                                         
                  PHI( 16)= ZER                                         
                  PHI( 17)= ZER                                         
                  PHI( 18)= ZER                                         
                  DO I= 37, 43                                          
                     PHI(I)= ZER                                        
                  END DO                                                
                  ITOP=IPRIM                                            
                  IF(IEQJ) ITOP=JJ                                      
                  DO II=1,ITOP                                          
                     IJ=IJ+1                                            
                     CFB= CCBRA(IJ)                                     
                     IF(ABS(CFB).GT.CCTOL) THEN                         
                        X01= EXI(II)                                    
                        X12= X01+X02                                    
                        X41= ONE/(X12+X34)                              
                        FTZ= CFB*CFK                                    
                        RX = XKL-XAB(IJ)                                
                        RY = YKL-YAB(IJ)                                
                        RZ = ZKL-ZAB(IJ)                                
                        RSQ= RX*RX+RY*RY+RZ*RZ                          
                        RHO= X12*X34*X41                                
                        IF(LRINT) THEN                                  
                           EFR= EMU2/(EMU2+RHO)                         
                           RHO= RHO*EFR                                 
                        ENDIF                                           
                        TT = RSQ*RHO                                    
                        N=2                                             
                        IF(TT.LE.TMAX) THEN                             
C                                                                       
C     FM(T) EVALUATION                                                  
C                                                                       
                           TV= TT*RFINC(N)                              
                           IP= NINT(TV)                                 
                           FX=    FGRID(4,IP,N) *TV                     
                           FX=(FX+FGRID(3,IP,N))*TV                     
                           FX=(FX+FGRID(2,IP,N))*TV                     
                           FX=(FX+FGRID(1,IP,N))*TV                     
                           FX= FX+FGRID(0,IP,N)                         
                           TV= TT*RXINC                                 
                           IP= NINT(TV)                                 
                           ET=    XGRID(4,IP) *TV                       
                           ET=(ET+XGRID(3,IP))*TV                       
                           ET=(ET+XGRID(2,IP))*TV                       
                           ET=(ET+XGRID(1,IP))*TV                       
                           ET= ET+XGRID(0,IP)                           
C                                                                       
                           FT(N)= FX                                    
                           T2= TT+TT                                    
                              FT(2-1)=(T2*FT(2)+ET)*RMR(2)              
                              FT(1-1)=(T2*FT(1)+ET)*RMR(1)              
CC                         DO M=N,1,-1                                  
CC                            FT(M-1)=(T2*FT(M)+ET)*RMR(M)              
CC                         END DO                                       
C                                                                       
                           RHO= RHO+RHO                                 
                           IF(LRINT) FTZ= FTZ*SQRT(EFR)                 
                           FTF= FTZ*SQRT(X41)                           
                              FT(0)= FT(0)*FTF                          
                           FTF= FTF*RHO                                 
                              FT(1)= FT(1)*FTF                          
                           FTF= FTF*RHO                                 
                              FT(2)= FT(2)*FTF                          
CC                         DO 210 M=0,N                                 
CC                            FT(M)= FT(M)*FTF                          
CC210                      FTF= FTF*RHO                                 
                        ELSE                                            
                           XIN= ONE/RSQ                                 
                           FTF= FTZ*SQRT(XIN)*SLBRA(IJ)*SLKET(KL)       
                           FT(0)= FTF                                   
                              FTF= FTF*XIN                              
                           FT(1)= TLGM(1)*FTF                           
                              FTF= FTF*XIN                              
                           FT(2)= TLGM(2)*FTF                           
CC                         DO 220 M=1,N                                 
CC                            FTF= FTF*XIN                              
CC220                      FT(M)= TLGM(M)*FTF                           
                        END IF                                          
C                                                                       
                        CALL PHIFTS(3,RX,RY,RZ,PHI,FT)                  
C                                                                       
                        FAC= RXB(IJ)                                    
                        SFAC(2)= FAC                                    
                        PHI( 12)= PHI( 12)+PHI(  1)*SFAC(2)             
C                                                                       
                        PHI( 16)= PHI( 16)+PHI(  2)*SFAC(2)             
                        PHI( 17)= PHI( 17)+PHI(  3)*SFAC(2)             
                        PHI( 18)= PHI( 18)+PHI(  4)*SFAC(2)             
                        J=  5                                           
                        DO I= 37, 43                                    
                           PHI(I)= PHI(I)+PHI(J)*SFAC(2)                
                           J=J+1                                        
                        END DO                                          
                     END IF                                             
                  END DO                                                
                  FAC= EXJ(JJ)*TWO                                      
                  SFAC(2)= FAC                                          
                  PHI( 13)= PHI( 13)+PHI( 12)*SFAC(2)                   
C                                                                       
                  PHI( 19)= PHI( 19)+PHI( 16)                           
                  PHI( 20)= PHI( 20)+PHI( 17)                           
                  PHI( 21)= PHI( 21)+PHI( 18)                           
                  PHI( 28)= PHI( 28)+PHI( 16)*SFAC(2)                   
                  PHI( 29)= PHI( 29)+PHI( 17)*SFAC(2)                   
                  PHI( 30)= PHI( 30)+PHI( 18)*SFAC(2)                   
                  J= 37                                                 
                  DO I= 44, 50                                          
                     PHI(I)= PHI(I)+PHI(J)                              
                     J=J+1                                              
                  END DO                                                
               END DO                                                   
               FAC= RXK(KL)                                             
               SFAC(2)= FAC                                             
               PHI( 14)= PHI( 14)+PHI( 13)*SFAC(2)                      
C                                                                       
               PHI( 22)= PHI( 22)+PHI( 19)*SFAC(2)                      
               PHI( 23)= PHI( 23)+PHI( 20)*SFAC(2)                      
               PHI( 24)= PHI( 24)+PHI( 21)*SFAC(2)                      
               PHI( 31)= PHI( 31)+PHI( 28)*SFAC(2)                      
               PHI( 32)= PHI( 32)+PHI( 29)*SFAC(2)                      
               PHI( 33)= PHI( 33)+PHI( 30)*SFAC(2)                      
               J= 44                                                    
               DO I= 51, 57                                             
                  PHI(I)= PHI(I)+PHI(J)*SFAC(2)                         
                  J=J+1                                                 
               END DO                                                   
            END IF                                                      
         END DO                                                         
         FAC= EXL(LL)*TWO                                               
         SFAC(2)= FAC                                                   
         PHI( 15)= PHI( 15)+PHI( 14)*SFAC(2)                            
C                                                                       
         PHI( 25)= PHI( 25)+PHI( 22)*SFAC(2)                            
         PHI( 26)= PHI( 26)+PHI( 23)*SFAC(2)                            
         PHI( 27)= PHI( 27)+PHI( 24)*SFAC(2)                            
         PHI( 34)= PHI( 34)+PHI( 31)                                    
         PHI( 35)= PHI( 35)+PHI( 32)                                    
         PHI( 36)= PHI( 36)+PHI( 33)                                    
         J= 51                                                          
         DO I= 58, 64                                                   
            PHI(I)= PHI(I)+PHI(J)                                       
            J=J+1                                                       
         END DO                                                         
      END DO                                                            
C  POST-CONTRACTION PHASE                                               
      PHI( 59)= PHI( 59)-PHI( 58)                                       
      PHI( 61)= PHI( 61)-PHI( 58)                                       
      PHI( 64)= PHI( 64)-PHI( 58)                                       
      CNF(1)= XB-XA                                                     
      CNF(2)= YB-YA                                                     
      CNF(3)= ZB-ZA                                                     
      WK1(  2)= PHI( 25)+PHI( 15)*CNF( 1)                               
      WK1(  3)= PHI( 26)+PHI( 15)*CNF( 2)                               
      WK1(  4)= PHI( 27)+PHI( 15)*CNF( 3)                               
      CALL PSPS_BCTE(WK1,WK2,LENW, 1)                                   
      I= 59                                                             
      J= 34                                                             
      L=3                                                               
      DO K=  2,  4                                                      
         IF(K.EQ. 4) I=I+1                                              
         IF(K.EQ. 4) L=L-1                                              
         WK1(  2)=-PHI(I  )-PHI(J)*CNF( 1)                              
         WK1(  3)=-PHI(I+1)-PHI(J)*CNF( 2)                              
         WK1(  4)=-PHI(I+L)-PHI(J)*CNF( 3)                              
         CALL PSPS_BCTE(WK1,WK2,LENW, K)                                
         I=I+1                                                          
         J=J+1                                                          
      END DO                                                            
      CNF(1)= XD-XC                                                     
      CNF(2)= YD-YC                                                     
      CNF(3)= ZD-ZC                                                     
      DO I=1,3                                                          
         WK2(  2,I)= WK2(  2,I)+WK2(  1,I)*CNF( 1)                      
         WK2(  3,I)= WK2(  3,I)+WK2(  1,I)*CNF( 2)                      
         WK2(  4,I)= WK2(  4,I)+WK2(  1,I)*CNF( 3)                      
      END DO                                                            
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2C   *DECK PSPS_BCTE                                        
C>                                                                      
C>    @brief   ERIC psps backtransfer                                   
C>                                                                      
C>    @details ERIC [ps|ps] backtransfer step                           
C>                                                                      
C>    @author  Graham Fletcher, 2004, modified Jose Sierra, 2013.       
C>                                                                      
      SUBROUTINE PSPS_BCTE(WK1,WK2,LENW,JR)                             
C                                                                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C                                                                       
      INTEGER    LENW,I,J,JR                                            
      DIMENSION  WK1(*),WK2(LENW,*)                                     
C                                                                       
      J=  1                                                             
      DO I=1,3                                                          
         J=J+1                                                          
         WK2(JR,I)= WK1(J)                                              
      END DO                                                            
C                                                                       
      RETURN                                                            
      END                                                               
C  LPHI=  65                                                            
C  LWK1=   1                                                            
C  LWK2=  24                                                            
C  LENW=  24                                                            
C*MODULE INT2C   *DECK PPSS                                             
C>                                                                      
C>    @brief   ERIC ppss case                                           
C>                                                                      
C>    @details ERIC [pp|ss] integral quartet                            
C>                                                                      
C>    @author  Graham Fletcher, 2004, modified Jose Sierra, 2013.       
C>                                                                      
      SUBROUTINE PPSS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL,               
     *                 XC,YC,ZC,XD,YD,ZD,PHI,WK2,LENW)                  
      USE lrcdft, ONLY: EMU2                                            
      use mx_limits, only: mxgsh,mxg2                                   
C                                                                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C                                                                       
      LOGICAL    IEQJ,KEQL                                              
      DIMENSION  PHI(*),WK2(LENW,*)                                     
C                                                                       
C                                                                       
      PARAMETER (NTX=4)                                                 
      PARAMETER (NPF=450)                                               
      PARAMETER (NGRD=7)                                                
      PARAMETER (NPX=1000)                                              
      PARAMETER (MXQT=16)                                               
      COMMON /FMTTBL/ FGRID(0:NTX,0:NPF,0:NGRD),XGRID(0:NTX,0:NPX)      
     *,               TMAX,RFINC(0:NGRD),RXINC                          
     *,               RMR(MXQT),TLGM(0:MXQT),NORD                       
      COMMON /ERIPRM/ EXI(MXGSH),EXJ(MXGSH),EXK(MXGSH),EXL(MXGSH),      
     *                CCI(MXGSH),CCJ(MXGSH),CCK(MXGSH),CCL(MXGSH),      
     *                XAB(MXG2),YAB(MXG2),ZAB(MXG2),                    
     *                XCD(MXG2),YCD(MXG2),ZCD(MXG2),                    
     *                CCBRA(MXG2),CCKET(MXG2),RXB(MXG2),                
     *                SLBRA(MXG2),SLKET(MXG2),RXK(MXG2)                 
C$omp threadprivate(/ERIPRM/)
      LOGICAL         LRINT                                             
      COMMON /NLRCF / LRINT                                             
C$omp threadprivate(/NLRCF /)
      COMMON /SHLNOS/ QQ4,IDUMMY(20)                                    
C$omp threadprivate(/SHLNOS/)
C                                                                       
      DIMENSION  FT(0:16),SFAC( 9),CNF(11)                              
C                                                                       
      PARAMETER (ZER=0.0D+00)                                           
      PARAMETER (ONE=1.0D+00)                                           
      PARAMETER (TWO=2.0D+00)                                           
      PARAMETER (CCTOL=1.0D-13)                                         
      DATA CNF/0.0D+00,0.0D+00,0.0D+00,1.0D+00,2.0D+00,3.0D+00,         
     &         4.0D+00,5.0D+00,6.0D+00,7.0D+00,8.0D+00/                 
      SAVE CNF                                                          
!$omp threadprivate(cnf)                                                
C                                                                       
      SFAC(1)= ONE                                                      
      PHI( 15)= ZER                                                     
      PHI( 16)= ZER                                                     
      PHI( 18)= ZER                                                     
      PHI( 28)= ZER                                                     
      PHI( 29)= ZER                                                     
      PHI( 30)= ZER                                                     
      PHI( 34)= ZER                                                     
      PHI( 35)= ZER                                                     
      PHI( 36)= ZER                                                     
      DO I= 58, 64                                                      
         PHI(I)= ZER                                                    
      END DO                                                            
      KL=0                                                              
      DO LL=1,LPRIM                                                     
         X04= EXL(LL)                                                   
         PHI( 14)= ZER                                                  
         PHI( 17)= ZER                                                  
         PHI( 25)= ZER                                                  
         PHI( 26)= ZER                                                  
         PHI( 27)= ZER                                                  
         PHI( 31)= ZER                                                  
         PHI( 32)= ZER                                                  
         PHI( 33)= ZER                                                  
         DO I= 51, 57                                                   
            PHI(I)= ZER                                                 
         END DO                                                         
         KTOP=KPRIM                                                     
         IF(KEQL) KTOP=LL                                               
         DO KK=1,KTOP                                                   
            KL=KL+1                                                     
            CFK= CCKET(KL)                                              
            IF(ABS(CFK).GT.CCTOL) THEN                                  
               CFK= CFK*QQ4                                             
               X03= EXK(KK)                                             
               X34= X03+X04                                             
               XKL= XCD(KL)                                             
               YKL= YCD(KL)                                             
               ZKL= ZCD(KL)                                             
               PHI( 13)= ZER                                            
               PHI( 22)= ZER                                            
               PHI( 23)= ZER                                            
               PHI( 24)= ZER                                            
               DO I= 44, 50                                             
                  PHI(I)= ZER                                           
               END DO                                                   
               IJ=0                                                     
               DO JJ=1,JPRIM                                            
                  X02= EXJ(JJ)                                          
                  PHI( 12)= ZER                                         
                  PHI( 19)= ZER                                         
                  PHI( 20)= ZER                                         
                  PHI( 21)= ZER                                         
                  DO I= 37, 43                                          
                     PHI(I)= ZER                                        
                  END DO                                                
                  ITOP=IPRIM                                            
                  IF(IEQJ) ITOP=JJ                                      
                  DO II=1,ITOP                                          
                     IJ=IJ+1                                            
                     CFB= CCBRA(IJ)                                     
                     IF(ABS(CFB).GT.CCTOL) THEN                         
                        X01= EXI(II)                                    
                        X12= X01+X02                                    
                        X41= ONE/(X12+X34)                              
                        FTZ= CFB*CFK                                    
                        RX = XKL-XAB(IJ)                                
                        RY = YKL-YAB(IJ)                                
                        RZ = ZKL-ZAB(IJ)                                
                        RSQ= RX*RX+RY*RY+RZ*RZ                          
                        RHO= X12*X34*X41                                
                        IF(LRINT) THEN                                  
                           EFR= EMU2/(EMU2+RHO)                         
                           RHO= RHO*EFR                                 
                        ENDIF                                           
                        TT = RSQ*RHO                                    
                        N=2                                             
                        IF(TT.LE.TMAX) THEN                             
C                                                                       
C     FM(T) EVALUATION                                                  
C                                                                       
                           TV= TT*RFINC(N)                              
                           IP= NINT(TV)                                 
                           FX=    FGRID(4,IP,N) *TV                     
                           FX=(FX+FGRID(3,IP,N))*TV                     
                           FX=(FX+FGRID(2,IP,N))*TV                     
                           FX=(FX+FGRID(1,IP,N))*TV                     
                           FX= FX+FGRID(0,IP,N)                         
                           TV= TT*RXINC                                 
                           IP= NINT(TV)                                 
                           ET=    XGRID(4,IP) *TV                       
                           ET=(ET+XGRID(3,IP))*TV                       
                           ET=(ET+XGRID(2,IP))*TV                       
                           ET=(ET+XGRID(1,IP))*TV                       
                           ET= ET+XGRID(0,IP)                           
C                                                                       
                           FT(N)= FX                                    
                           T2= TT+TT                                    
                              FT(2-1)=(T2*FT(2)+ET)*RMR(2)              
                              FT(1-1)=(T2*FT(1)+ET)*RMR(1)              
CC                         DO M=N,1,-1                                  
CC                            FT(M-1)=(T2*FT(M)+ET)*RMR(M)              
CC                         END DO                                       
C                                                                       
                           RHO= RHO+RHO                                 
                           IF(LRINT) FTZ= FTZ*SQRT(EFR)                 
                           FTF= FTZ*SQRT(X41)                           
                              FT(0)= FT(0)*FTF                          
                           FTF= FTF*RHO                                 
                              FT(1)= FT(1)*FTF                          
                           FTF= FTF*RHO                                 
                              FT(2)= FT(2)*FTF                          
CC                         DO 210 M=0,N                                 
CC                            FT(M)= FT(M)*FTF                          
CC210                      FTF= FTF*RHO                                 
                        ELSE                                            
                           XIN= ONE/RSQ                                 
                           FTF= FTZ*SQRT(XIN)*SLBRA(IJ)*SLKET(KL)       
                           FT(0)= FTF                                   
                              FTF= FTF*XIN                              
                           FT(1)= TLGM(1)*FTF                           
                              FTF= FTF*XIN                              
                           FT(2)= TLGM(2)*FTF                           
CC                         DO 220 M=1,N                                 
CC                            FTF= FTF*XIN                              
CC220                      FT(M)= TLGM(M)*FTF                           
                        END IF                                          
C                                                                       
                        CALL PHIFTS(3,RX,RY,RZ,PHI,FT)                  
C                                                                       
C                       FAC= RXB(IJ)                                    
                        PHI( 12)= PHI( 12)+PHI(  1)                     
C                                                                       
                        PHI( 19)= PHI( 19)+PHI(  2)                     
                        PHI( 20)= PHI( 20)+PHI(  3)                     
                        PHI( 21)= PHI( 21)+PHI(  4)                     
                        J=  5                                           
                        DO I= 37, 43                                    
                           PHI(I)= PHI(I)+PHI(J)                        
                           J=J+1                                        
                        END DO                                          
                     END IF                                             
                  END DO                                                
C                 FAC= EXJ(JJ)*TWO                                      
                  PHI( 13)= PHI( 13)+PHI( 12)                           
C                                                                       
                  PHI( 22)= PHI( 22)+PHI( 19)                           
                  PHI( 23)= PHI( 23)+PHI( 20)                           
                  PHI( 24)= PHI( 24)+PHI( 21)                           
                  J= 37                                                 
                  DO I= 44, 50                                          
                     PHI(I)= PHI(I)+PHI(J)                              
                     J=J+1                                              
                  END DO                                                
               END DO                                                   
               FAC= RXK(KL)                                             
               SFAC(2)= FAC                                             
               SFAC(3)= SFAC(2)*FAC                                     
               PHI( 14)= PHI( 14)+PHI( 13)*SFAC(2)                      
               PHI( 17)= PHI( 17)+PHI( 13)*SFAC(3)                      
C                                                                       
               PHI( 25)= PHI( 25)+PHI( 22)*SFAC(2)                      
               PHI( 26)= PHI( 26)+PHI( 23)*SFAC(2)                      
               PHI( 27)= PHI( 27)+PHI( 24)*SFAC(2)                      
               PHI( 31)= PHI( 31)+PHI( 22)*SFAC(3)                      
               PHI( 32)= PHI( 32)+PHI( 23)*SFAC(3)                      
               PHI( 33)= PHI( 33)+PHI( 24)*SFAC(3)                      
               J= 44                                                    
               DO I= 51, 57                                             
                  PHI(I)= PHI(I)+PHI(J)*SFAC(3)                         
                  J=J+1                                                 
               END DO                                                   
            END IF                                                      
         END DO                                                         
         FAC= EXL(LL)*TWO                                               
         SFAC(2)= FAC                                                   
         SFAC(3)= SFAC(2)*FAC                                           
         PHI( 15)= PHI( 15)+PHI( 14)                                    
         PHI( 16)= PHI( 16)+PHI( 14)*SFAC(2)                            
         PHI( 18)= PHI( 18)+PHI( 17)*SFAC(3)                            
C                                                                       
         PHI( 28)= PHI( 28)+PHI( 25)                                    
         PHI( 29)= PHI( 29)+PHI( 26)                                    
         PHI( 30)= PHI( 30)+PHI( 27)                                    
         PHI( 34)= PHI( 34)+PHI( 31)*SFAC(2)                            
         PHI( 35)= PHI( 35)+PHI( 32)*SFAC(2)                            
         PHI( 36)= PHI( 36)+PHI( 33)*SFAC(2)                            
         J= 51                                                          
         DO I= 58, 64                                                   
            PHI(I)= PHI(I)+PHI(J)                                       
            J=J+1                                                       
         END DO                                                         
      END DO                                                            
C  POST-CONTRACTION PHASE                                               
      PHI( 59)= PHI( 59)-PHI( 58)                                       
      PHI( 61)= PHI( 61)-PHI( 58)                                       
      PHI( 64)= PHI( 64)-PHI( 58)                                       
      WK2(  1,1)= PHI( 16)                                              
      WK2(  2,1)=-PHI( 28)                                              
      WK2(  3,1)=-PHI( 29)                                              
      WK2(  4,1)=-PHI( 30)                                              
      WK2(  5,1)= PHI( 18)                                              
      WK2(  6,1)= PHI( 15)                                              
      WK2(  7,1)=-PHI( 34)                                              
      WK2(  8,1)=-PHI( 35)                                              
      WK2(  9,1)=-PHI( 36)                                              
      WK2( 10,1)= PHI( 59)                                              
      WK2( 11,1)= PHI( 60)                                              
      WK2( 12,1)= PHI( 61)                                              
      WK2( 13,1)= PHI( 62)                                              
      WK2( 14,1)= PHI( 63)                                              
      WK2( 15,1)= PHI( 64)                                              
      CNF(1)= XD-XC                                                     
      CNF(2)= YD-YC                                                     
      CNF(3)= ZD-ZC                                                     
      DO I=1,1                                                          
         CALL PPWRK1(I,WK2,LENW,CNF)                                    
      END DO                                                            
C                                                                       
      RETURN                                                            
      END                                                               
C  LPHI= 183                                                            
C  LWK1=   4                                                            
C  LWK2=  72                                                            
C  LENW=  24                                                            
C*MODULE INT2C   *DECK PPPS                                             
C>                                                                      
C>    @brief   ERIC ppps case                                           
C>                                                                      
C>    @details ERIC [pp|ps] integral quartet                            
C>                                                                      
C>    @author  Graham Fletcher, 2004, modified Jose Sierra, 2013.       
C>                                                                      
      SUBROUTINE PPPS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL,               
     *                 XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD,             
     *                 PHI,WK1,WK2,LENW)                                
      USE lrcdft, ONLY: EMU2                                            
      use mx_limits, only: mxgsh,mxg2                                   
C                                                                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C                                                                       
      LOGICAL    IEQJ,KEQL                                              
      DIMENSION  PHI(*),WK1(*),WK2(LENW,*)                              
C                                                                       
C                                                                       
      PARAMETER (NTX=4)                                                 
      PARAMETER (NPF=450)                                               
      PARAMETER (NGRD=7)                                                
      PARAMETER (NPX=1000)                                              
      PARAMETER (MXQT=16)                                               
      COMMON /FMTTBL/ FGRID(0:NTX,0:NPF,0:NGRD),XGRID(0:NTX,0:NPX)      
     *,               TMAX,RFINC(0:NGRD),RXINC                          
     *,               RMR(MXQT),TLGM(0:MXQT),NORD                       
      COMMON /ERIPRM/ EXI(MXGSH),EXJ(MXGSH),EXK(MXGSH),EXL(MXGSH),      
     *                CCI(MXGSH),CCJ(MXGSH),CCK(MXGSH),CCL(MXGSH),      
     *                XAB(MXG2),YAB(MXG2),ZAB(MXG2),                    
     *                XCD(MXG2),YCD(MXG2),ZCD(MXG2),                    
     *                CCBRA(MXG2),CCKET(MXG2),RXB(MXG2),                
     *                SLBRA(MXG2),SLKET(MXG2),RXK(MXG2)                 
C$omp threadprivate(/ERIPRM/)
      LOGICAL         LRINT                                             
      COMMON /NLRCF / LRINT                                             
C$omp threadprivate(/NLRCF /)
      COMMON /SHLNOS/ QQ4,IDUMMY(20)                                    
C$omp threadprivate(/SHLNOS/)
C                                                                       
      DIMENSION  FT(0:16),SFAC( 9),CNF(11)                              
C                                                                       
      PARAMETER (NL=  4)                                                
      PARAMETER (NK=  4)                                                
      PARAMETER (NJ=  3)                                                
      PARAMETER (NI=  2)                                                
      DIMENSION  LSF(NL),KSF(NK),JSF(NJ),ISF(NI)                        
      DIMENSION  LJI(NL),KJI(NK),JJI(NJ),IJI(NI)                        
      DIMENSION  L0I(NL),K0I(NK),J0I(NJ),I0I(NI)                        
      DIMENSION  L0F(NL),K0F(NK),J0F(NJ),I0F(NI)                        
C                                                                       
      PARAMETER (ZER=0.0D+00)                                           
      PARAMETER (ONE=1.0D+00)                                           
      PARAMETER (TWO=2.0D+00)                                           
      PARAMETER (CCTOL=1.0D-13)                                         
      DATA CNF/0.0D+00,0.0D+00,0.0D+00,1.0D+00,2.0D+00,3.0D+00,         
     &         4.0D+00,5.0D+00,6.0D+00,7.0D+00,8.0D+00/                 
      SAVE CNF                                                          
!$omp threadprivate(cnf)                                                
C                                                                       
      DATA LSF/  1,  2,  1,  1/                                         
      DATA LJI/ 82, 96,117,157/                                         
      DATA L0I/ 89,103,124,170/                                         
      DATA L0F/ 95,109,130,182/                                         
C                                                                       
      DATA KSF/  2,  3,  3,  3/                                         
      DATA KJI/ 75, 75,110,144/                                         
      DATA K0I/ 82, 96,117,157/                                         
      DATA K0F/ 88,102,123,169/                                         
C                                                                       
      DATA JSF/  1,  2,  1/                                             
      DATA JJI/ 68, 68,131/                                             
      DATA J0I/ 75,110,144/                                             
      DATA J0F/ 81,116,156/                                             
C                                                                       
      DATA ISF/  2,  2/                                                 
      DATA IJI/  5, 12/                                                 
      DATA I0I/ 68,131/                                                 
      DATA I0F/ 74,143/                                                 
C                                                                       
      SFAC(1)= ONE                                                      
      PHI( 28)= ZER                                                     
      PHI( 29)= ZER                                                     
      PHI( 31)= ZER                                                     
      PHI( 41)= ZER                                                     
      PHI( 42)= ZER                                                     
      PHI( 43)= ZER                                                     
      PHI( 44)= ZER                                                     
      PHI( 45)= ZER                                                     
      PHI( 46)= ZER                                                     
      PHI( 50)= ZER                                                     
      PHI( 51)= ZER                                                     
      PHI( 52)= ZER                                                     
      PHI( 59)= ZER                                                     
      PHI( 60)= ZER                                                     
      PHI( 61)= ZER                                                     
      PHI( 65)= ZER                                                     
      PHI( 66)= ZER                                                     
      PHI( 67)= ZER                                                     
      DO N= 1,NL                                                        
         DO I=L0I(N),L0F(N)                                             
            PHI(I)= ZER                                                 
         END DO                                                         
      END DO                                                            
      KL=0                                                              
      DO LL=1,LPRIM                                                     
         X04= EXL(LL)                                                   
         PHI( 27)= ZER                                                  
         PHI( 30)= ZER                                                  
         PHI( 38)= ZER                                                  
         PHI( 39)= ZER                                                  
         PHI( 40)= ZER                                                  
         PHI( 47)= ZER                                                  
         PHI( 48)= ZER                                                  
         PHI( 49)= ZER                                                  
         PHI( 56)= ZER                                                  
         PHI( 57)= ZER                                                  
         PHI( 58)= ZER                                                  
         PHI( 62)= ZER                                                  
         PHI( 63)= ZER                                                  
         PHI( 64)= ZER                                                  
         DO N= 1,NK                                                     
            DO I=K0I(N),K0F(N)                                          
               PHI(I)= ZER                                              
            END DO                                                      
         END DO                                                         
         KTOP=KPRIM                                                     
         IF(KEQL) KTOP=LL                                               
         DO KK=1,KTOP                                                   
            KL=KL+1                                                     
            CFK= CCKET(KL)                                              
            IF(ABS(CFK).GT.CCTOL) THEN                                  
               CFK= CFK*QQ4                                             
               X03= EXK(KK)                                             
               X34= X03+X04                                             
               XKL= XCD(KL)                                             
               YKL= YCD(KL)                                             
               ZKL= ZCD(KL)                                             
               PHI( 26)= ZER                                            
               PHI( 35)= ZER                                            
               PHI( 36)= ZER                                            
               PHI( 37)= ZER                                            
               PHI( 53)= ZER                                            
               PHI( 54)= ZER                                            
               PHI( 55)= ZER                                            
               DO N= 1,NJ                                               
                  DO I=J0I(N),J0F(N)                                    
                     PHI(I)= ZER                                        
                  END DO                                                
               END DO                                                   
               IJ=0                                                     
               DO JJ=1,JPRIM                                            
                  X02= EXJ(JJ)                                          
                  PHI( 25)= ZER                                         
                  PHI( 32)= ZER                                         
                  PHI( 33)= ZER                                         
                  PHI( 34)= ZER                                         
                  DO N= 1,NI                                            
                     DO I=I0I(N),I0F(N)                                 
                        PHI(I)= ZER                                     
                     END DO                                             
                  END DO                                                
                  ITOP=IPRIM                                            
                  IF(IEQJ) ITOP=JJ                                      
                  DO II=1,ITOP                                          
                     IJ=IJ+1                                            
                     CFB= CCBRA(IJ)                                     
                     IF(ABS(CFB).GT.CCTOL) THEN                         
                        X01= EXI(II)                                    
                        X12= X01+X02                                    
                        X41= ONE/(X12+X34)                              
                        FTZ= CFB*CFK                                    
                        RX = XKL-XAB(IJ)                                
                        RY = YKL-YAB(IJ)                                
                        RZ = ZKL-ZAB(IJ)                                
                        RSQ= RX*RX+RY*RY+RZ*RZ                          
                        RHO= X12*X34*X41                                
                        IF(LRINT) THEN                                  
                           EFR= EMU2/(EMU2+RHO)                         
                           RHO= RHO*EFR                                 
                        ENDIF                                           
                        TT = RSQ*RHO                                    
                        N=3                                             
                        IF(TT.LE.TMAX) THEN                             
C                                                                       
C     FM(T) EVALUATION                                                  
C                                                                       
                           TV= TT*RFINC(N)                              
                           IP= NINT(TV)                                 
                           FX=    FGRID(4,IP,N) *TV                     
                           FX=(FX+FGRID(3,IP,N))*TV                     
                           FX=(FX+FGRID(2,IP,N))*TV                     
                           FX=(FX+FGRID(1,IP,N))*TV                     
                           FX= FX+FGRID(0,IP,N)                         
                           TV= TT*RXINC                                 
                           IP= NINT(TV)                                 
                           ET=    XGRID(4,IP) *TV                       
                           ET=(ET+XGRID(3,IP))*TV                       
                           ET=(ET+XGRID(2,IP))*TV                       
                           ET=(ET+XGRID(1,IP))*TV                       
                           ET= ET+XGRID(0,IP)                           
C                                                                       
                           FT(N)= FX                                    
                           T2= TT+TT                                    
                              FT(3-1)=(T2*FT(3)+ET)*RMR(3)              
                              FT(2-1)=(T2*FT(2)+ET)*RMR(2)              
                              FT(1-1)=(T2*FT(1)+ET)*RMR(1)              
CC                         DO M=N,1,-1                                  
CC                            FT(M-1)=(T2*FT(M)+ET)*RMR(M)              
CC                         END DO                                       
C                                                                       
                           RHO= RHO+RHO                                 
                           IF(LRINT) FTZ= FTZ*SQRT(EFR)                 
                           FTF= FTZ*SQRT(X41)                           
                              FT(0)= FT(0)*FTF                          
                           FTF= FTF*RHO                                 
                              FT(1)= FT(1)*FTF                          
                           FTF= FTF*RHO                                 
                              FT(2)= FT(2)*FTF                          
                           FTF= FTF*RHO                                 
                              FT(3)= FT(3)*FTF                          
CC                         DO 210 M=0,N                                 
CC                            FT(M)= FT(M)*FTF                          
CC210                      FTF= FTF*RHO                                 
                        ELSE                                            
                           XIN= ONE/RSQ                                 
                           FTF= FTZ*SQRT(XIN)*SLBRA(IJ)*SLKET(KL)       
                           FT(0)= FTF                                   
                              FTF= FTF*XIN                              
                           FT(1)= TLGM(1)*FTF                           
                              FTF= FTF*XIN                              
                           FT(2)= TLGM(2)*FTF                           
                              FTF= FTF*XIN                              
                           FT(3)= TLGM(3)*FTF                           
CC                         DO 220 M=1,N                                 
CC                            FTF= FTF*XIN                              
CC220                      FT(M)= TLGM(M)*FTF                           
                        END IF                                          
C                                                                       
                        CALL PHIFTS(4,RX,RY,RZ,PHI,FT)                  
C                                                                       
                        FAC= RXB(IJ)                                    
                        SFAC(2)= FAC                                    
                        PHI( 25)= PHI( 25)+PHI(  1)*SFAC(2)             
C                                                                       
                        PHI( 32)= PHI( 32)+PHI(  2)*SFAC(2)             
                        PHI( 33)= PHI( 33)+PHI(  3)*SFAC(2)             
                        PHI( 34)= PHI( 34)+PHI(  4)*SFAC(2)             
                        DO N= 1,NI                                      
                           FAC= SFAC(ISF(N))                            
                           J=IJI(N)                                     
                           DO I=I0I(N),I0F(N)                           
                              PHI(I)= PHI(I)+PHI(J)*FAC                 
                              J=J+1                                     
                           END DO                                       
                        END DO                                          
                     END IF                                             
                  END DO                                                
                  FAC= EXJ(JJ)*TWO                                      
                  SFAC(2)= FAC                                          
                  PHI( 26)= PHI( 26)+PHI( 25)*SFAC(2)                   
C                                                                       
                  PHI( 35)= PHI( 35)+PHI( 32)                           
                  PHI( 36)= PHI( 36)+PHI( 33)                           
                  PHI( 37)= PHI( 37)+PHI( 34)                           
                  PHI( 53)= PHI( 53)+PHI( 32)*SFAC(2)                   
                  PHI( 54)= PHI( 54)+PHI( 33)*SFAC(2)                   
                  PHI( 55)= PHI( 55)+PHI( 34)*SFAC(2)                   
                  DO N= 1,NJ                                            
                     FAC= SFAC(JSF(N))                                  
                     J=JJI(N)                                           
                     IF(FAC.EQ.ONE) THEN                                
                        DO I=J0I(N),J0F(N)                              
                           PHI(I)= PHI(I)+PHI(J)                        
                           J=J+1                                        
                        END DO                                          
                     ELSE                                               
                        DO I=J0I(N),J0F(N)                              
                           PHI(I)= PHI(I)+PHI(J)*FAC                    
                           J=J+1                                        
                        END DO                                          
                     END IF                                             
                  END DO                                                
               END DO                                                   
               FAC= RXK(KL)                                             
               SFAC(2)= FAC                                             
               SFAC(3)= SFAC(2)*FAC                                     
               PHI( 27)= PHI( 27)+PHI( 26)*SFAC(2)                      
               PHI( 30)= PHI( 30)+PHI( 26)*SFAC(3)                      
C                                                                       
               PHI( 38)= PHI( 38)+PHI( 35)*SFAC(2)                      
               PHI( 39)= PHI( 39)+PHI( 36)*SFAC(2)                      
               PHI( 40)= PHI( 40)+PHI( 37)*SFAC(2)                      
               PHI( 47)= PHI( 47)+PHI( 35)*SFAC(3)                      
               PHI( 48)= PHI( 48)+PHI( 36)*SFAC(3)                      
               PHI( 49)= PHI( 49)+PHI( 37)*SFAC(3)                      
               PHI( 56)= PHI( 56)+PHI( 53)*SFAC(2)                      
               PHI( 57)= PHI( 57)+PHI( 54)*SFAC(2)                      
               PHI( 58)= PHI( 58)+PHI( 55)*SFAC(2)                      
               PHI( 62)= PHI( 62)+PHI( 53)*SFAC(3)                      
               PHI( 63)= PHI( 63)+PHI( 54)*SFAC(3)                      
               PHI( 64)= PHI( 64)+PHI( 55)*SFAC(3)                      
               DO N= 1,NK                                               
                  FAC= SFAC(KSF(N))                                     
                  J=KJI(N)                                              
                  DO I=K0I(N),K0F(N)                                    
                     PHI(I)= PHI(I)+PHI(J)*FAC                          
                     J=J+1                                              
                  END DO                                                
               END DO                                                   
            END IF                                                      
         END DO                                                         
         FAC= EXL(LL)*TWO                                               
         SFAC(2)= FAC                                                   
         SFAC(3)= SFAC(2)*FAC                                           
         PHI( 28)= PHI( 28)+PHI( 27)                                    
         PHI( 29)= PHI( 29)+PHI( 27)*SFAC(2)                            
         PHI( 31)= PHI( 31)+PHI( 30)*SFAC(3)                            
C                                                                       
         PHI( 41)= PHI( 41)+PHI( 38)                                    
         PHI( 42)= PHI( 42)+PHI( 39)                                    
         PHI( 43)= PHI( 43)+PHI( 40)                                    
         PHI( 44)= PHI( 44)+PHI( 38)*SFAC(2)                            
         PHI( 45)= PHI( 45)+PHI( 39)*SFAC(2)                            
         PHI( 46)= PHI( 46)+PHI( 40)*SFAC(2)                            
         PHI( 50)= PHI( 50)+PHI( 47)*SFAC(3)                            
         PHI( 51)= PHI( 51)+PHI( 48)*SFAC(3)                            
         PHI( 52)= PHI( 52)+PHI( 49)*SFAC(3)                            
         PHI( 59)= PHI( 59)+PHI( 56)                                    
         PHI( 60)= PHI( 60)+PHI( 57)                                    
         PHI( 61)= PHI( 61)+PHI( 58)                                    
         PHI( 65)= PHI( 65)+PHI( 62)*SFAC(2)                            
         PHI( 66)= PHI( 66)+PHI( 63)*SFAC(2)                            
         PHI( 67)= PHI( 67)+PHI( 64)*SFAC(2)                            
         DO N= 1,NL                                                     
            FAC= SFAC(LSF(N))                                           
            J=LJI(N)                                                    
            IF(FAC.EQ.ONE) THEN                                         
               DO I=L0I(N),L0F(N)                                       
                  PHI(I)= PHI(I)+PHI(J)                                 
                  J=J+1                                                 
               END DO                                                   
            ELSE                                                        
               DO I=L0I(N),L0F(N)                                       
                  PHI(I)= PHI(I)+PHI(J)*FAC                             
                  J=J+1                                                 
               END DO                                                   
            END IF                                                      
         END DO                                                         
      END DO                                                            
C  POST-CONTRACTION PHASE                                               
      PHI( 90)= PHI( 90)-PHI( 89)                                       
      PHI( 92)= PHI( 92)-PHI( 89)                                       
      PHI( 95)= PHI( 95)-PHI( 89)                                       
      PHI(104)= PHI(104)-PHI(103)                                       
      PHI(106)= PHI(106)-PHI(103)                                       
      PHI(109)= PHI(109)-PHI(103)                                       
      PHI(125)= PHI(125)-PHI(124)                                       
      PHI(127)= PHI(127)-PHI(124)                                       
      PHI(130)= PHI(130)-PHI(124)                                       
      III=130+ 30                                                       
      III=III+ 13                                                       
      CALL PHIIJ3(III,PHI)                                              
C                                                                       
      CNF(1)= XB-XA                                                     
      CNF(2)= YB-YA                                                     
      CNF(3)= ZB-ZA                                                     
      WK1(  2)= PHI( 44)+PHI( 29)*CNF( 1)                               
      WK1(  3)= PHI( 45)+PHI( 29)*CNF( 2)                               
      WK1(  4)= PHI( 46)+PHI( 29)*CNF( 3)                               
      CALL PPPS_BCTE(WK1,WK2,LENW, 1)                                   
      I= 90                                                             
      J= 59                                                             
      L=3                                                               
      DO K=  2,  4                                                      
         IF(K.EQ. 4) I=I+1                                              
         IF(K.EQ. 4) L=L-1                                              
         WK1(  2)=-PHI(I  )-PHI(J)*CNF( 1)                              
         WK1(  3)=-PHI(I+1)-PHI(J)*CNF( 2)                              
         WK1(  4)=-PHI(I+L)-PHI(J)*CNF( 3)                              
         CALL PPPS_BCTE(WK1,WK2,LENW, K)                                
         I=I+1                                                          
         J=J+1                                                          
      END DO                                                            
      WK1(  2)= PHI( 50)+PHI( 31)*CNF( 1)                               
      WK1(  3)= PHI( 51)+PHI( 31)*CNF( 2)                               
      WK1(  4)= PHI( 52)+PHI( 31)*CNF( 3)                               
      CALL PPPS_BCTE(WK1,WK2,LENW, 5)                                   
      WK1(  2)= PHI( 41)+PHI( 28)*CNF( 1)                               
      WK1(  3)= PHI( 42)+PHI( 28)*CNF( 2)                               
      WK1(  4)= PHI( 43)+PHI( 28)*CNF( 3)                               
      CALL PPPS_BCTE(WK1,WK2,LENW, 6)                                   
      I=104                                                             
      J= 65                                                             
      L=3                                                               
      DO K=  7,  9                                                      
         IF(K.EQ. 9) I=I+1                                              
         IF(K.EQ. 9) L=L-1                                              
         WK1(  2)=-PHI(I  )-PHI(J)*CNF( 1)                              
         WK1(  3)=-PHI(I+1)-PHI(J)*CNF( 2)                              
         WK1(  4)=-PHI(I+L)-PHI(J)*CNF( 3)                              
         CALL PPPS_BCTE(WK1,WK2,LENW, K)                                
         I=I+1                                                          
         J=J+1                                                          
      END DO                                                            
      I=173                                                             
      J=125                                                             
      L=4                                                               
      DO K= 10, 15                                                      
         IF(K.EQ.13 .OR. K.EQ.15) I=I+1                                 
         IF(K.EQ.13 .OR. K.EQ.15) L=L-1                                 
         WK1(  2)= PHI(I  )+PHI(J)*CNF( 1)                              
         WK1(  3)= PHI(I+1)+PHI(J)*CNF( 2)                              
         WK1(  4)= PHI(I+L)+PHI(J)*CNF( 3)                              
         CALL PPPS_BCTE(WK1,WK2,LENW, K)                                
         I=I+1                                                          
         J=J+1                                                          
      END DO                                                            
      CNF(1)= XD-XC                                                     
      CNF(2)= YD-YC                                                     
      CNF(3)= ZD-ZC                                                     
      DO I=1,3                                                          
         CALL PPWRK1(I,WK2,LENW,CNF)                                    
      END DO                                                            
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2C   *DECK PPPS_BCTE                                        
C>                                                                      
C>    @brief   ERIC ppps backtransfer                                   
C>                                                                      
C>    @details ERIC [pp|ps] backtransfer step                           
C>                                                                      
C>    @author  Graham Fletcher, 2004, modified Jose Sierra, 2013.       
C>                                                                      
      SUBROUTINE PPPS_BCTE(WK1,WK2,LENW,JR)                             
C                                                                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C                                                                       
      INTEGER    LENW,I,J,JR                                            
      DIMENSION  WK1(*),WK2(LENW,*)                                     
C                                                                       
      J=  1                                                             
      DO I=1,3                                                          
         J=J+1                                                          
         WK2(JR,I)= WK1(J)                                              
      END DO                                                            
C                                                                       
      RETURN                                                            
      END                                                               
C  LPHI= 574                                                            
C  LWK1=  24                                                            
C  LWK2= 216                                                            
C  LENW=  24                                                            
C*MODULE INT2C   *DECK PPPP                                             
C>                                                                      
C>    @brief   ERIC pppp case                                           
C>                                                                      
C>    @details ERIC [pp|pp] integral quartet                            
C>                                                                      
C>    @author  Graham Fletcher, 2004, modified Jose Sierra, 2013.       
C>                                                                      
      SUBROUTINE PPPP (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL,               
     *                 XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD,             
     *                 PHI,WK1,WK2,LENW)                                
      USE lrcdft, ONLY: EMU2                                            
      use mx_limits, only: mxgsh,mxg2                                   
C                                                                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C                                                                       
      LOGICAL    IEQJ,KEQL                                              
      DIMENSION  PHI(*),WK1(*),WK2(LENW,*)                              
C                                                                       
C                                                                       
      PARAMETER (NTX=4)                                                 
      PARAMETER (NPF=450)                                               
      PARAMETER (NGRD=7)                                                
      PARAMETER (NPX=1000)                                              
      PARAMETER (MXQT=16)                                               
      COMMON /FMTTBL/ FGRID(0:NTX,0:NPF,0:NGRD),XGRID(0:NTX,0:NPX)      
     *,               TMAX,RFINC(0:NGRD),RXINC                          
     *,               RMR(MXQT),TLGM(0:MXQT),NORD                       
      COMMON /ERIPRM/ EXI(MXGSH),EXJ(MXGSH),EXK(MXGSH),EXL(MXGSH),      
     *                CCI(MXGSH),CCJ(MXGSH),CCK(MXGSH),CCL(MXGSH),      
     *                XAB(MXG2),YAB(MXG2),ZAB(MXG2),                    
     *                XCD(MXG2),YCD(MXG2),ZCD(MXG2),                    
     *                CCBRA(MXG2),CCKET(MXG2),RXB(MXG2),                
     *                SLBRA(MXG2),SLKET(MXG2),RXK(MXG2)                 
C$omp threadprivate(/ERIPRM/)
      LOGICAL         LRINT                                             
      COMMON /NLRCF / LRINT                                             
C$omp threadprivate(/NLRCF /)
      COMMON /SHLNOS/ QQ4,IDUMMY(20)                                    
C$omp threadprivate(/SHLNOS/)
C                                                                       
      DIMENSION  FT(0:16),SFAC( 9),CNF(11)                              
C                                                                       
      PARAMETER (NL= 15)                                                
      PARAMETER (NK= 13)                                                
      PARAMETER (NJ=  9)                                                
      PARAMETER (NI=  5)                                                
      DIMENSION  LSF(NL),KSF(NK),JSF(NJ),ISF(NI)                        
      DIMENSION  LJI(NL),KJI(NK),JJI(NJ),IJI(NI)                        
      DIMENSION  L0I(NL),K0I(NK),J0I(NJ),I0I(NI)                        
      DIMENSION  L0F(NL),K0F(NK),J0F(NJ),I0F(NI)                        
C                                                                       
      PARAMETER (ZER=0.0D+00)                                           
      PARAMETER (ONE=1.0D+00)                                           
      PARAMETER (TWO=2.0D+00)                                           
      PARAMETER (CCTOL=1.0D-13)                                         
      DATA CNF/0.0D+00,0.0D+00,0.0D+00,1.0D+00,2.0D+00,3.0D+00,         
     &         4.0D+00,5.0D+00,6.0D+00,7.0D+00,8.0D+00/                 
      SAVE CNF                                                          
!$omp threadprivate(cnf)                                                
C                                                                       
      DATA LSF/  1,  1,  2,  1,  1,  2,  3,  1,  2,  1,                 
     *           1,  1,  2,  1,  1/                                     
      DATA LJI/156,170,170,198,226,226,247,268,282,303,                 
     *         343,395,421,460,530/                                     
      DATA L0I/163,177,184,205,233,240,254,275,289,310,                 
     *         356,408,434,473,552/                                     
      DATA L0F/169,183,190,211,239,246,260,281,295,316,                 
     *         368,420,446,485,573/                                     
C                                                                       
      DATA KSF/  2,  3,  3,  2,  3,  2,  3,  3,  3,  2,                 
     *           3,  3,  3/                                             
      DATA KJI/149,149,191,219,219,261,261,296,330,382,                 
     *         382,447,508/                                             
      DATA K0I/156,170,198,226,247,268,282,303,343,395,                 
     *         421,460,530/                                             
      DATA K0F/162,176,204,232,253,274,288,309,355,407,                 
     *         433,472,551/                                             
C                                                                       
      DATA JSF/  1,  2,  1,  2,  3,  1,  1,  2,  1/                     
      DATA JJI/142,142,212,212,212,317,369,369,486/                     
      DATA J0I/149,191,219,261,296,330,382,447,508/                     
      DATA J0F/155,197,225,267,302,342,394,459,529/                     
C                                                                       
      DATA ISF/  2,  3,  2,  3,  3/                                     
      DATA IJI/  5,  5, 12, 12, 25/                                     
      DATA I0I/142,212,317,369,486/                                     
      DATA I0F/148,218,329,381,507/                                     
C                                                                       
      SFAC(1)= ONE                                                      
      PHI( 50)= ZER                                                     
      PHI( 51)= ZER                                                     
      PHI( 53)= ZER                                                     
      PHI( 56)= ZER                                                     
      PHI( 57)= ZER                                                     
      PHI( 59)= ZER                                                     
      PHI( 63)= ZER                                                     
      PHI( 64)= ZER                                                     
      PHI( 66)= ZER                                                     
      PHI( 76)= ZER                                                     
      PHI( 77)= ZER                                                     
      PHI( 78)= ZER                                                     
      PHI( 79)= ZER                                                     
      PHI( 80)= ZER                                                     
      PHI( 81)= ZER                                                     
      PHI( 85)= ZER                                                     
      PHI( 86)= ZER                                                     
      PHI( 87)= ZER                                                     
      PHI( 88)= ZER                                                     
      PHI( 89)= ZER                                                     
      PHI( 90)= ZER                                                     
      PHI( 97)= ZER                                                     
      PHI( 98)= ZER                                                     
      PHI( 99)= ZER                                                     
      PHI(103)= ZER                                                     
      PHI(104)= ZER                                                     
      PHI(105)= ZER                                                     
      PHI(115)= ZER                                                     
      PHI(116)= ZER                                                     
      PHI(117)= ZER                                                     
      PHI(118)= ZER                                                     
      PHI(119)= ZER                                                     
      PHI(120)= ZER                                                     
      PHI(124)= ZER                                                     
      PHI(125)= ZER                                                     
      PHI(126)= ZER                                                     
      PHI(133)= ZER                                                     
      PHI(134)= ZER                                                     
      PHI(135)= ZER                                                     
      PHI(139)= ZER                                                     
      PHI(140)= ZER                                                     
      PHI(141)= ZER                                                     
      DO N= 1,NL                                                        
         DO I=L0I(N),L0F(N)                                             
            PHI(I)= ZER                                                 
         END DO                                                         
      END DO                                                            
      KL=0                                                              
      DO LL=1,LPRIM                                                     
         X04= EXL(LL)                                                   
         PHI( 49)= ZER                                                  
         PHI( 52)= ZER                                                  
         PHI( 55)= ZER                                                  
         PHI( 58)= ZER                                                  
         PHI( 62)= ZER                                                  
         PHI( 65)= ZER                                                  
         PHI( 73)= ZER                                                  
         PHI( 74)= ZER                                                  
         PHI( 75)= ZER                                                  
         PHI( 82)= ZER                                                  
         PHI( 83)= ZER                                                  
         PHI( 84)= ZER                                                  
         PHI( 94)= ZER                                                  
         PHI( 95)= ZER                                                  
         PHI( 96)= ZER                                                  
         PHI(100)= ZER                                                  
         PHI(101)= ZER                                                  
         PHI(102)= ZER                                                  
         PHI(112)= ZER                                                  
         PHI(113)= ZER                                                  
         PHI(114)= ZER                                                  
         PHI(121)= ZER                                                  
         PHI(122)= ZER                                                  
         PHI(123)= ZER                                                  
         PHI(130)= ZER                                                  
         PHI(131)= ZER                                                  
         PHI(132)= ZER                                                  
         PHI(136)= ZER                                                  
         PHI(137)= ZER                                                  
         PHI(138)= ZER                                                  
         DO N= 1,NK                                                     
            DO I=K0I(N),K0F(N)                                          
               PHI(I)= ZER                                              
            END DO                                                      
         END DO                                                         
         KTOP=KPRIM                                                     
         IF(KEQL) KTOP=LL                                               
         DO KK=1,KTOP                                                   
            KL=KL+1                                                     
            CFK= CCKET(KL)                                              
            IF(ABS(CFK).GT.CCTOL) THEN                                  
               CFK= CFK*QQ4                                             
               X03= EXK(KK)                                             
               X34= X03+X04                                             
               XKL= XCD(KL)                                             
               YKL= YCD(KL)                                             
               ZKL= ZCD(KL)                                             
               PHI( 48)= ZER                                            
               PHI( 54)= ZER                                            
               PHI( 61)= ZER                                            
               PHI( 70)= ZER                                            
               PHI( 71)= ZER                                            
               PHI( 72)= ZER                                            
               PHI( 91)= ZER                                            
               PHI( 92)= ZER                                            
               PHI( 93)= ZER                                            
               PHI(109)= ZER                                            
               PHI(110)= ZER                                            
               PHI(111)= ZER                                            
               PHI(127)= ZER                                            
               PHI(128)= ZER                                            
               PHI(129)= ZER                                            
               DO N= 1,NJ                                               
                  DO I=J0I(N),J0F(N)                                    
                     PHI(I)= ZER                                        
                  END DO                                                
               END DO                                                   
               IJ=0                                                     
               DO JJ=1,JPRIM                                            
                  X02= EXJ(JJ)                                          
                  PHI( 47)= ZER                                         
                  PHI( 60)= ZER                                         
                  PHI( 67)= ZER                                         
                  PHI( 68)= ZER                                         
                  PHI( 69)= ZER                                         
                  PHI(106)= ZER                                         
                  PHI(107)= ZER                                         
                  PHI(108)= ZER                                         
                  DO N= 1,NI                                            
                     DO I=I0I(N),I0F(N)                                 
                        PHI(I)= ZER                                     
                     END DO                                             
                  END DO                                                
                  ITOP=IPRIM                                            
                  IF(IEQJ) ITOP=JJ                                      
                  DO II=1,ITOP                                          
                     IJ=IJ+1                                            
                     CFB= CCBRA(IJ)                                     
                     IF(ABS(CFB).GT.CCTOL) THEN                         
                        X01= EXI(II)                                    
                        X12= X01+X02                                    
                        X41= ONE/(X12+X34)                              
                        FTZ= CFB*CFK                                    
                        RX = XKL-XAB(IJ)                                
                        RY = YKL-YAB(IJ)                                
                        RZ = ZKL-ZAB(IJ)                                
                        RSQ= RX*RX+RY*RY+RZ*RZ                          
                        RHO= X12*X34*X41                                
                        IF(LRINT) THEN                                  
                           EFR= EMU2/(EMU2+RHO)                         
                           RHO= RHO*EFR                                 
                        ENDIF                                           
                        TT = RSQ*RHO                                    
                        N=4                                             
                        IF(TT.LE.TMAX) THEN                             
C                                                                       
C     FM(T) EVALUATION                                                  
C                                                                       
                           TV= TT*RFINC(N)                              
                           IP= NINT(TV)                                 
                           FX=    FGRID(4,IP,N) *TV                     
                           FX=(FX+FGRID(3,IP,N))*TV                     
                           FX=(FX+FGRID(2,IP,N))*TV                     
                           FX=(FX+FGRID(1,IP,N))*TV                     
                           FX= FX+FGRID(0,IP,N)                         
                           TV= TT*RXINC                                 
                           IP= NINT(TV)                                 
                           ET=    XGRID(4,IP) *TV                       
                           ET=(ET+XGRID(3,IP))*TV                       
                           ET=(ET+XGRID(2,IP))*TV                       
                           ET=(ET+XGRID(1,IP))*TV                       
                           ET= ET+XGRID(0,IP)                           
C                                                                       
                           FT(N)= FX                                    
                           T2= TT+TT                                    
                           DO M=N,1,-1                                  
                              FT(M-1)=(T2*FT(M)+ET)*RMR(M)              
                           END DO                                       
C                                                                       
                           RHO= RHO+RHO                                 
                           IF(LRINT) FTZ= FTZ*SQRT(EFR)                 
                           FTF= FTZ*SQRT(X41)                           
                           DO 210 M=0,N                                 
                              FT(M)= FT(M)*FTF                          
  210                      FTF= FTF*RHO                                 
                        ELSE                                            
                           XIN= ONE/RSQ                                 
                           FTF= FTZ*SQRT(XIN)*SLBRA(IJ)*SLKET(KL)       
                           FT(0)= FTF                                   
                           DO 220 M=1,N                                 
                              FTF= FTF*XIN                              
  220                      FT(M)= TLGM(M)*FTF                           
                        END IF                                          
C                                                                       
                        CALL PHIFTS(5,RX,RY,RZ,PHI,FT)                  
C                                                                       
                        FAC= RXB(IJ)                                    
                        SFAC(2)= FAC                                    
                        SFAC(3)= SFAC(2)*FAC                            
                        PHI( 47)= PHI( 47)+PHI(  1)*SFAC(2)             
                        PHI( 60)= PHI( 60)+PHI(  1)*SFAC(3)             
C                                                                       
                        PHI( 67)= PHI( 67)+PHI(  2)*SFAC(2)             
                        PHI( 68)= PHI( 68)+PHI(  3)*SFAC(2)             
                        PHI( 69)= PHI( 69)+PHI(  4)*SFAC(2)             
                        PHI(106)= PHI(106)+PHI(  2)*SFAC(3)             
                        PHI(107)= PHI(107)+PHI(  3)*SFAC(3)             
                        PHI(108)= PHI(108)+PHI(  4)*SFAC(3)             
                        DO N= 1,NI                                      
                           FAC= SFAC(ISF(N))                            
                           J=IJI(N)                                     
                           DO I=I0I(N),I0F(N)                           
                              PHI(I)= PHI(I)+PHI(J)*FAC                 
                              J=J+1                                     
                           END DO                                       
                        END DO                                          
                     END IF                                             
                  END DO                                                
                  FAC= EXJ(JJ)*TWO                                      
                  SFAC(2)= FAC                                          
                  SFAC(3)= SFAC(2)*FAC                                  
                  PHI( 48)= PHI( 48)+PHI( 47)                           
                  PHI( 54)= PHI( 54)+PHI( 47)*SFAC(2)                   
                  PHI( 61)= PHI( 61)+PHI( 60)*SFAC(3)                   
C                                                                       
                  PHI( 70)= PHI( 70)+PHI( 67)                           
                  PHI( 71)= PHI( 71)+PHI( 68)                           
                  PHI( 72)= PHI( 72)+PHI( 69)                           
                  PHI( 91)= PHI( 91)+PHI( 67)*SFAC(2)                   
                  PHI( 92)= PHI( 92)+PHI( 68)*SFAC(2)                   
                  PHI( 93)= PHI( 93)+PHI( 69)*SFAC(2)                   
                  PHI(109)= PHI(109)+PHI(106)*SFAC(2)                   
                  PHI(110)= PHI(110)+PHI(107)*SFAC(2)                   
                  PHI(111)= PHI(111)+PHI(108)*SFAC(2)                   
                  PHI(127)= PHI(127)+PHI(106)*SFAC(3)                   
                  PHI(128)= PHI(128)+PHI(107)*SFAC(3)                   
                  PHI(129)= PHI(129)+PHI(108)*SFAC(3)                   
                  DO N= 1,NJ                                            
                     FAC= SFAC(JSF(N))                                  
                     J=JJI(N)                                           
                     IF(FAC.EQ.ONE) THEN                                
                        DO I=J0I(N),J0F(N)                              
                           PHI(I)= PHI(I)+PHI(J)                        
                           J=J+1                                        
                        END DO                                          
                     ELSE                                               
                        DO I=J0I(N),J0F(N)                              
                           PHI(I)= PHI(I)+PHI(J)*FAC                    
                           J=J+1                                        
                        END DO                                          
                     END IF                                             
                  END DO                                                
               END DO                                                   
               FAC= RXK(KL)                                             
               SFAC(2)= FAC                                             
               SFAC(3)= SFAC(2)*FAC                                     
               PHI( 49)= PHI( 49)+PHI( 48)*SFAC(2)                      
               PHI( 52)= PHI( 52)+PHI( 48)*SFAC(3)                      
               PHI( 55)= PHI( 55)+PHI( 54)*SFAC(2)                      
               PHI( 58)= PHI( 58)+PHI( 54)*SFAC(3)                      
               PHI( 62)= PHI( 62)+PHI( 61)*SFAC(2)                      
               PHI( 65)= PHI( 65)+PHI( 61)*SFAC(3)                      
C                                                                       
               PHI( 73)= PHI( 73)+PHI( 70)*SFAC(2)                      
               PHI( 74)= PHI( 74)+PHI( 71)*SFAC(2)                      
               PHI( 75)= PHI( 75)+PHI( 72)*SFAC(2)                      
               PHI( 82)= PHI( 82)+PHI( 70)*SFAC(3)                      
               PHI( 83)= PHI( 83)+PHI( 71)*SFAC(3)                      
               PHI( 84)= PHI( 84)+PHI( 72)*SFAC(3)                      
               PHI( 94)= PHI( 94)+PHI( 91)*SFAC(2)                      
               PHI( 95)= PHI( 95)+PHI( 92)*SFAC(2)                      
               PHI( 96)= PHI( 96)+PHI( 93)*SFAC(2)                      
               PHI(100)= PHI(100)+PHI( 91)*SFAC(3)                      
               PHI(101)= PHI(101)+PHI( 92)*SFAC(3)                      
               PHI(102)= PHI(102)+PHI( 93)*SFAC(3)                      
               PHI(112)= PHI(112)+PHI(109)*SFAC(2)                      
               PHI(113)= PHI(113)+PHI(110)*SFAC(2)                      
               PHI(114)= PHI(114)+PHI(111)*SFAC(2)                      
               PHI(121)= PHI(121)+PHI(109)*SFAC(3)                      
               PHI(122)= PHI(122)+PHI(110)*SFAC(3)                      
               PHI(123)= PHI(123)+PHI(111)*SFAC(3)                      
               PHI(130)= PHI(130)+PHI(127)*SFAC(2)                      
               PHI(131)= PHI(131)+PHI(128)*SFAC(2)                      
               PHI(132)= PHI(132)+PHI(129)*SFAC(2)                      
               PHI(136)= PHI(136)+PHI(127)*SFAC(3)                      
               PHI(137)= PHI(137)+PHI(128)*SFAC(3)                      
               PHI(138)= PHI(138)+PHI(129)*SFAC(3)                      
               DO N= 1,NK                                               
                  FAC= SFAC(KSF(N))                                     
                  J=KJI(N)                                              
                  DO I=K0I(N),K0F(N)                                    
                     PHI(I)= PHI(I)+PHI(J)*FAC                          
                     J=J+1                                              
                  END DO                                                
               END DO                                                   
            END IF                                                      
         END DO                                                         
         FAC= EXL(LL)*TWO                                               
         SFAC(2)= FAC                                                   
         SFAC(3)= SFAC(2)*FAC                                           
         PHI( 50)= PHI( 50)+PHI( 49)                                    
         PHI( 51)= PHI( 51)+PHI( 49)*SFAC(2)                            
         PHI( 53)= PHI( 53)+PHI( 52)*SFAC(3)                            
         PHI( 56)= PHI( 56)+PHI( 55)                                    
         PHI( 57)= PHI( 57)+PHI( 55)*SFAC(2)                            
         PHI( 59)= PHI( 59)+PHI( 58)*SFAC(3)                            
         PHI( 63)= PHI( 63)+PHI( 62)                                    
         PHI( 64)= PHI( 64)+PHI( 62)*SFAC(2)                            
         PHI( 66)= PHI( 66)+PHI( 65)*SFAC(3)                            
C                                                                       
         PHI( 76)= PHI( 76)+PHI( 73)                                    
         PHI( 77)= PHI( 77)+PHI( 74)                                    
         PHI( 78)= PHI( 78)+PHI( 75)                                    
         PHI( 79)= PHI( 79)+PHI( 73)*SFAC(2)                            
         PHI( 80)= PHI( 80)+PHI( 74)*SFAC(2)                            
         PHI( 81)= PHI( 81)+PHI( 75)*SFAC(2)                            
         PHI( 85)= PHI( 85)+PHI( 82)*SFAC(2)                            
         PHI( 86)= PHI( 86)+PHI( 83)*SFAC(2)                            
         PHI( 87)= PHI( 87)+PHI( 84)*SFAC(2)                            
         PHI( 88)= PHI( 88)+PHI( 82)*SFAC(3)                            
         PHI( 89)= PHI( 89)+PHI( 83)*SFAC(3)                            
         PHI( 90)= PHI( 90)+PHI( 84)*SFAC(3)                            
         PHI( 97)= PHI( 97)+PHI( 94)                                    
         PHI( 98)= PHI( 98)+PHI( 95)                                    
         PHI( 99)= PHI( 99)+PHI( 96)                                    
         PHI(103)= PHI(103)+PHI(100)*SFAC(2)                            
         PHI(104)= PHI(104)+PHI(101)*SFAC(2)                            
         PHI(105)= PHI(105)+PHI(102)*SFAC(2)                            
         PHI(115)= PHI(115)+PHI(112)                                    
         PHI(116)= PHI(116)+PHI(113)                                    
         PHI(117)= PHI(117)+PHI(114)                                    
         PHI(118)= PHI(118)+PHI(112)*SFAC(2)                            
         PHI(119)= PHI(119)+PHI(113)*SFAC(2)                            
         PHI(120)= PHI(120)+PHI(114)*SFAC(2)                            
         PHI(124)= PHI(124)+PHI(121)*SFAC(3)                            
         PHI(125)= PHI(125)+PHI(122)*SFAC(3)                            
         PHI(126)= PHI(126)+PHI(123)*SFAC(3)                            
         PHI(133)= PHI(133)+PHI(130)                                    
         PHI(134)= PHI(134)+PHI(131)                                    
         PHI(135)= PHI(135)+PHI(132)                                    
         PHI(139)= PHI(139)+PHI(136)*SFAC(2)                            
         PHI(140)= PHI(140)+PHI(137)*SFAC(2)                            
         PHI(141)= PHI(141)+PHI(138)*SFAC(2)                            
         DO N= 1,NL                                                     
            FAC= SFAC(LSF(N))                                           
            J=LJI(N)                                                    
            IF(FAC.EQ.ONE) THEN                                         
               DO I=L0I(N),L0F(N)                                       
                  PHI(I)= PHI(I)+PHI(J)                                 
                  J=J+1                                                 
               END DO                                                   
            ELSE                                                        
               DO I=L0I(N),L0F(N)                                       
                  PHI(I)= PHI(I)+PHI(J)*FAC                             
                  J=J+1                                                 
               END DO                                                   
            END IF                                                      
         END DO                                                         
      END DO                                                            
C  POST-CONTRACTION PHASE                                               
      PHI(164)= PHI(164)-PHI(163)                                       
      PHI(166)= PHI(166)-PHI(163)                                       
      PHI(169)= PHI(169)-PHI(163)                                       
      PHI(178)= PHI(178)-PHI(177)                                       
      PHI(180)= PHI(180)-PHI(177)                                       
      PHI(183)= PHI(183)-PHI(177)                                       
      PHI(185)= PHI(185)-PHI(184)                                       
      PHI(187)= PHI(187)-PHI(184)                                       
      PHI(190)= PHI(190)-PHI(184)                                       
      PHI(206)= PHI(206)-PHI(205)                                       
      PHI(208)= PHI(208)-PHI(205)                                       
      PHI(211)= PHI(211)-PHI(205)                                       
      PHI(234)= PHI(234)-PHI(233)                                       
      PHI(236)= PHI(236)-PHI(233)                                       
      PHI(239)= PHI(239)-PHI(233)                                       
      PHI(241)= PHI(241)-PHI(240)                                       
      PHI(243)= PHI(243)-PHI(240)                                       
      PHI(246)= PHI(246)-PHI(240)                                       
      PHI(255)= PHI(255)-PHI(254)                                       
      PHI(257)= PHI(257)-PHI(254)                                       
      PHI(260)= PHI(260)-PHI(254)                                       
      PHI(276)= PHI(276)-PHI(275)                                       
      PHI(278)= PHI(278)-PHI(275)                                       
      PHI(281)= PHI(281)-PHI(275)                                       
      PHI(290)= PHI(290)-PHI(289)                                       
      PHI(292)= PHI(292)-PHI(289)                                       
      PHI(295)= PHI(295)-PHI(289)                                       
      PHI(311)= PHI(311)-PHI(310)                                       
      PHI(313)= PHI(313)-PHI(310)                                       
      PHI(316)= PHI(316)-PHI(310)                                       
      III=316+ 30                                                       
      III=III+ 13                                                       
      CALL PHIIJ3(III,PHI)                                              
      III=III+ 13*4                                                     
      CALL PHIIJ3(III,PHI)                                              
      III=III+ 13*2                                                     
      CALL PHIIJ3(III,PHI)                                              
      III=III+ 13*3                                                     
      CALL PHIIJ3(III,PHI)                                              
C                                                                       
      III=III+ 61                                                       
      III=III+ 22                                                       
      CALL PHIIJ4(III,PHI)                                              
C                                                                       
      CNF(1)= XB-XA                                                     
      CNF(2)= YB-YA                                                     
      CNF(3)= ZB-ZA                                                     
      WK1(  2)=          PHI( 79)        +PHI( 57)*CNF( 1)              
      WK1(  3)=          PHI( 80)        +PHI( 57)*CNF( 2)              
      WK1(  4)=          PHI( 81)        +PHI( 57)*CNF( 3)              
      WK1(  7)=          PHI(118)        +PHI( 64)*CNF( 1)              
      WK1(  8)=          PHI(119)        +PHI( 64)*CNF( 2)              
      WK1(  9)=          PHI(120)        +PHI( 64)*CNF( 3)              
      WK1( 10)= PHI(241)+PHI(118)*CNF( 1)+PHI( 51)*CNF( 4)              
      WK1( 11)= PHI(242)+PHI(119)*CNF( 1)                               
      WK1( 12)= PHI(243)+PHI(119)*CNF( 2)+PHI( 51)*CNF( 4)              
      WK1( 13)= PHI(244)+PHI(120)*CNF( 1)                               
      WK1( 14)= PHI(245)+PHI(120)*CNF( 2)                               
      WK1( 15)= PHI(246)+PHI(120)*CNF( 3)+PHI( 51)*CNF( 4)              
      CALL PPPP_BCTE(WK1,WK2,LENW,CNF, 1)                               
      I=164                                                             
      J= 97                                                             
      L=3                                                               
      DO K=  2,  4                                                      
         IF(K.EQ. 4) I=I+1                                              
         IF(K.EQ. 4) L=L-1                                              
         LLL=L+L+1                                                      
         M=I+112                                                        
         N=M+135+3-L                                                    
         WK1(  2)=           -PHI(I   )        -PHI(J    )*CNF( 1)      
         WK1(  3)=           -PHI(I+ 1)        -PHI(J    )*CNF( 2)      
         WK1(  4)=           -PHI(I+ L)        -PHI(J    )*CNF( 3)      
         WK1(  7)=           -PHI(M   )        -PHI(J+ 36)*CNF( 1)      
         WK1(  8)=           -PHI(M+ 1)        -PHI(J+ 36)*CNF( 2)      
         WK1(  9)=           -PHI(M+ L)        -PHI(J+ 36)*CNF( 3)      
         WK1( 10)=-PHI(N    )-PHI(M   )*CNF( 1)-PHI(J- 21)*CNF( 4)      
         WK1( 11)=-PHI(N+  1)-PHI(M+ 1)*CNF( 1)                         
         WK1( 12)=-PHI(N+  2)-PHI(M+ 1)*CNF( 2)-PHI(J- 21)*CNF( 4)      
         WK1( 13)=-PHI(N+1+L)-PHI(M+ L)*CNF( 1)                         
         WK1( 14)=-PHI(N+2+L)-PHI(M+ L)*CNF( 2)                         
         WK1( 15)=-PHI(N+LLL)-PHI(M+ L)*CNF( 3)-PHI(J- 21)*CNF( 4)      
         CALL PPPP_BCTE(WK1,WK2,LENW,CNF, K)                            
         I=I+1                                                          
         J=J+1                                                          
      END DO                                                            
      WK1(  2)=          PHI( 88)        +PHI( 59)*CNF( 1)              
      WK1(  3)=          PHI( 89)        +PHI( 59)*CNF( 2)              
      WK1(  4)=          PHI( 90)        +PHI( 59)*CNF( 3)              
      WK1(  7)=          PHI(124)        +PHI( 66)*CNF( 1)              
      WK1(  8)=          PHI(125)        +PHI( 66)*CNF( 2)              
      WK1(  9)=          PHI(126)        +PHI( 66)*CNF( 3)              
      WK1( 10)= PHI(255)+PHI(124)*CNF( 1)+PHI( 53)*CNF( 4)              
      WK1( 11)= PHI(256)+PHI(125)*CNF( 1)                               
      WK1( 12)= PHI(257)+PHI(125)*CNF( 2)+PHI( 53)*CNF( 4)              
      WK1( 13)= PHI(258)+PHI(126)*CNF( 1)                               
      WK1( 14)= PHI(259)+PHI(126)*CNF( 2)                               
      WK1( 15)= PHI(260)+PHI(126)*CNF( 3)+PHI( 53)*CNF( 4)              
      CALL PPPP_BCTE(WK1,WK2,LENW,CNF, 5)                               
      WK1(  2)=          PHI( 76)        +PHI( 56)*CNF( 1)              
      WK1(  3)=          PHI( 77)        +PHI( 56)*CNF( 2)              
      WK1(  4)=          PHI( 78)        +PHI( 56)*CNF( 3)              
      WK1(  7)=          PHI(115)        +PHI( 63)*CNF( 1)              
      WK1(  8)=          PHI(116)        +PHI( 63)*CNF( 2)              
      WK1(  9)=          PHI(117)        +PHI( 63)*CNF( 3)              
      WK1( 10)= PHI(234)+PHI(115)*CNF( 1)+PHI( 50)*CNF( 4)              
      WK1( 11)= PHI(235)+PHI(116)*CNF( 1)                               
      WK1( 12)= PHI(236)+PHI(116)*CNF( 2)+PHI( 50)*CNF( 4)              
      WK1( 13)= PHI(237)+PHI(117)*CNF( 1)                               
      WK1( 14)= PHI(238)+PHI(117)*CNF( 2)                               
      WK1( 15)= PHI(239)+PHI(117)*CNF( 3)+PHI( 50)*CNF( 4)              
      CALL PPPP_BCTE(WK1,WK2,LENW,CNF, 6)                               
      I=185                                                             
      J=103                                                             
      L=3                                                               
      DO K=  7,  9                                                      
         IF(K.EQ. 9) I=I+1                                              
         IF(K.EQ. 9) L=L-1                                              
         LLL=L+L+1                                                      
         M=I+105                                                        
         N=M+147+3-L                                                    
         WK1(  2)=           -PHI(I   )        -PHI(J    )*CNF( 1)      
         WK1(  3)=           -PHI(I+ 1)        -PHI(J    )*CNF( 2)      
         WK1(  4)=           -PHI(I+ L)        -PHI(J    )*CNF( 3)      
         WK1(  7)=           -PHI(M   )        -PHI(J+ 36)*CNF( 1)      
         WK1(  8)=           -PHI(M+ 1)        -PHI(J+ 36)*CNF( 2)      
         WK1(  9)=           -PHI(M+ L)        -PHI(J+ 36)*CNF( 3)      
         WK1( 10)=-PHI(N    )-PHI(M   )*CNF( 1)-PHI(J- 18)*CNF( 4)      
         WK1( 11)=-PHI(N+  1)-PHI(M+ 1)*CNF( 1)                         
         WK1( 12)=-PHI(N+  2)-PHI(M+ 1)*CNF( 2)-PHI(J- 18)*CNF( 4)      
         WK1( 13)=-PHI(N+1+L)-PHI(M+ L)*CNF( 1)                         
         WK1( 14)=-PHI(N+2+L)-PHI(M+ L)*CNF( 2)                         
         WK1( 15)=-PHI(N+LLL)-PHI(M+ L)*CNF( 3)-PHI(J- 18)*CNF( 4)      
         CALL PPPP_BCTE(WK1,WK2,LENW,CNF, K)                            
         I=I+1                                                          
         J=J+1                                                          
      END DO                                                            
      I=359                                                             
      J=206                                                             
      L=4                                                               
      DO K= 10, 15                                                      
         IF(K.EQ.13 .OR. K.EQ.15) I=I+1                                 
         IF(K.EQ.13 .OR. K.EQ.15) L=L-1                                 
         LLL=L+L+1                                                      
         M=I+117                                                        
         N=M+ 83+4-L                                                    
         WK1(  2)=            PHI(I   )        +PHI(J    )*CNF( 1)      
         WK1(  3)=            PHI(I+ 1)        +PHI(J    )*CNF( 2)      
         WK1(  4)=            PHI(I+ L)        +PHI(J    )*CNF( 3)      
         WK1(  7)=            PHI(M   )        +PHI(J+105)*CNF( 1)      
         WK1(  8)=            PHI(M+ 1)        +PHI(J+105)*CNF( 2)      
         WK1(  9)=            PHI(M+ L)        +PHI(J+105)*CNF( 3)      
         WK1( 10)= PHI(N    )+PHI(M   )*CNF( 1)+PHI(J- 28)*CNF( 4)      
         WK1( 11)= PHI(N+  1)+PHI(M+ 1)*CNF( 1)                         
         WK1( 12)= PHI(N+  2)+PHI(M+ 1)*CNF( 2)+PHI(J- 28)*CNF( 4)      
         WK1( 13)= PHI(N+1+L)+PHI(M+ L)*CNF( 1)                         
         WK1( 14)= PHI(N+2+L)+PHI(M+ L)*CNF( 2)                         
         WK1( 15)= PHI(N+LLL)+PHI(M+ L)*CNF( 3)+PHI(J- 28)*CNF( 4)      
         CALL PPPP_BCTE(WK1,WK2,LENW,CNF, K)                            
         I=I+1                                                          
         J=J+1                                                          
      END DO                                                            
      CNF(1)= XD-XC                                                     
      CNF(2)= YD-YC                                                     
      CNF(3)= ZD-ZC                                                     
      DO I=1,9                                                          
         CALL PPWRK1(I,WK2,LENW,CNF)                                    
      END DO                                                            
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2C   *DECK PPPP_BCTE                                        
C>                                                                      
C>    @brief   ERIC pppp backtransfer                                   
C>                                                                      
C>    @details ERIC [pp|pp] backtransfer step                           
C>                                                                      
C>    @author  Graham Fletcher, 2004, modified Jose Sierra, 2013.       
C>                                                                      
      SUBROUTINE PPPP_BCTE(WK1,WK2,LENW,CNF,JR)                         
C                                                                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C                                                                       
      INTEGER    LENW,I,J,JR                                            
      DIMENSION  WK1(*),WK2(LENW,*),CNF(*)                              
C                                                                       
      WK1( 10)= WK1( 10)+WK1(  7)*CNF( 1)                               
      WK1( 11)= WK1( 11)+WK1(  7)*CNF( 2)                               
      WK1( 12)= WK1( 12)+WK1(  8)*CNF( 2)                               
      WK1( 13)= WK1( 13)+WK1(  7)*CNF( 3)                               
      WK1( 14)= WK1( 14)+WK1(  8)*CNF( 3)                               
      WK1( 15)= WK1( 15)+WK1(  9)*CNF( 3)                               
C                                                                       
      WK1( 16)= WK1( 10)-WK1(  2)*CNF( 1)                               
      WK1( 17)= WK1( 11)-WK1(  3)*CNF( 1)                               
      WK1( 18)= WK1( 13)-WK1(  4)*CNF( 1)                               
      WK1( 19)= WK1( 11)-WK1(  2)*CNF( 2)                               
      WK1( 20)= WK1( 12)-WK1(  3)*CNF( 2)                               
      WK1( 21)= WK1( 14)-WK1(  4)*CNF( 2)                               
      WK1( 22)= WK1( 13)-WK1(  2)*CNF( 3)                               
      WK1( 23)= WK1( 14)-WK1(  3)*CNF( 3)                               
      WK1( 24)= WK1( 15)-WK1(  4)*CNF( 3)                               
      J= 15                                                             
      DO I=1,9                                                          
         J=J+1                                                          
         WK2(JR,I)= WK1(J)                                              
      END DO                                                            
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2C   *DECK PPWRK1                                           
C>                                                                      
C>    @brief   ERIC pp utility                                          
C>                                                                      
C>    @details common expressions for pp products in ERIC               
C>                                                                      
C>    @author  Jose Sierra, 2013                                        
C>                                                                      
      SUBROUTINE PPWRK1(I,WK2,LENW,CNF)                                 
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      INTEGER    I,LENW                                                 
      DIMENSION  WK2(LENW,*),CNF(*)                                     
C                                                                       
         WK2(  2,I)= WK2(  2,I)+WK2(  1,I)*CNF( 1)                      
         WK2(  3,I)= WK2(  3,I)+WK2(  1,I)*CNF( 2)                      
         WK2(  4,I)= WK2(  4,I)+WK2(  1,I)*CNF( 3)                      
C                                                                       
         WK2( 10,I)= WK2( 10,I)+WK2(  7,I)*CNF( 1)+WK2(  6,I)*CNF( 4)   
         WK2( 11,I)= WK2( 11,I)+WK2(  8,I)*CNF( 1)                      
         WK2( 12,I)= WK2( 12,I)+WK2(  8,I)*CNF( 2)+WK2(  6,I)*CNF( 4)   
         WK2( 13,I)= WK2( 13,I)+WK2(  9,I)*CNF( 1)                      
         WK2( 14,I)= WK2( 14,I)+WK2(  9,I)*CNF( 2)                      
         WK2( 15,I)= WK2( 15,I)+WK2(  9,I)*CNF( 3)+WK2(  6,I)*CNF( 4)   
C                                                                       
         WK2(  7,I)= WK2(  7,I)+WK2(  5,I)*CNF( 1)                      
         WK2(  8,I)= WK2(  8,I)+WK2(  5,I)*CNF( 2)                      
         WK2(  9,I)= WK2(  9,I)+WK2(  5,I)*CNF( 3)                      
C                                                                       
         WK2( 10,I)= WK2( 10,I)+WK2(  7,I)*CNF( 1)                      
         WK2( 11,I)= WK2( 11,I)+WK2(  7,I)*CNF( 2)                      
         WK2( 12,I)= WK2( 12,I)+WK2(  8,I)*CNF( 2)                      
         WK2( 13,I)= WK2( 13,I)+WK2(  7,I)*CNF( 3)                      
         WK2( 14,I)= WK2( 14,I)+WK2(  8,I)*CNF( 3)                      
         WK2( 15,I)= WK2( 15,I)+WK2(  9,I)*CNF( 3)                      
C                                                                       
         WK2( 16,I)= WK2( 10,I)-WK2(  2,I)*CNF( 1)                      
         WK2( 17,I)= WK2( 11,I)-WK2(  3,I)*CNF( 1)                      
         WK2( 18,I)= WK2( 13,I)-WK2(  4,I)*CNF( 1)                      
         WK2( 19,I)= WK2( 11,I)-WK2(  2,I)*CNF( 2)                      
         WK2( 20,I)= WK2( 12,I)-WK2(  3,I)*CNF( 2)                      
         WK2( 21,I)= WK2( 14,I)-WK2(  4,I)*CNF( 2)                      
         WK2( 22,I)= WK2( 13,I)-WK2(  2,I)*CNF( 3)                      
         WK2( 23,I)= WK2( 14,I)-WK2(  3,I)*CNF( 3)                      
         WK2( 24,I)= WK2( 15,I)-WK2(  4,I)*CNF( 3)                      
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2C   *DECK PHIFTS                                           
C>                                                                      
C>    @brief   phase array first few initial value set utility          
C>                                                                      
C>    @details common expressions in ERIC                               
C>                                                                      
C>    @author  Jose Sierra, 2013                                        
C>                                                                      
      SUBROUTINE PHIFTS(N,RX,RY,RZ,PHI,FT)                              
C                                                                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C                                                                       
      DIMENSION  PHI(*),FT(0:*)                                         
C                                                                       
      DIMENSION  RXYZ(3),CXYZ(56)                                       
C                                                                       
      PARAMETER (ONE=1.0D+00)                                           
C                                                                       
      RXYZ(1)= RX                                                       
      RXYZ(2)= RY                                                       
      RXYZ(3)= RZ                                                       
C                                                                       
         CXYZ( 1)= ONE                                                  
         CXYZ( 2)= RXYZ(1)                                              
         CXYZ( 3)= RXYZ(2)                                              
         CXYZ( 4)= RXYZ(3)                                              
      IF(N.LE.2) GO TO 100                                              
         CXYZ( 5)= CXYZ( 2)*RXYZ(1)                                     
         CXYZ( 6)= CXYZ( 3)*RXYZ(1)                                     
         CXYZ( 7)= CXYZ( 3)*RXYZ(2)                                     
         CXYZ( 8)= CXYZ( 4)*RXYZ(1)                                     
         CXYZ( 9)= CXYZ( 4)*RXYZ(2)                                     
         CXYZ(10)= CXYZ( 4)*RXYZ(3)                                     
      IF(N.LE.3) GO TO 100                                              
         CXYZ(11)= CXYZ( 5)*RXYZ(1)                                     
         CXYZ(12)= CXYZ( 6)*RXYZ(1)                                     
         CXYZ(13)= CXYZ( 7)*RXYZ(1)                                     
         CXYZ(14)= CXYZ( 7)*RXYZ(2)                                     
         CXYZ(15)= CXYZ( 8)*RXYZ(1)                                     
         CXYZ(16)= CXYZ( 9)*RXYZ(1)                                     
         CXYZ(17)= CXYZ( 9)*RXYZ(2)                                     
         CXYZ(18)= CXYZ(10)*RXYZ(1)                                     
         CXYZ(19)= CXYZ(10)*RXYZ(2)                                     
         CXYZ(20)= CXYZ(10)*RXYZ(3)                                     
      IF(N.LE.4) GO TO 100                                              
         CXYZ(21)= CXYZ(11)*RXYZ(1)                                     
         CXYZ(22)= CXYZ(12)*RXYZ(1)                                     
         CXYZ(23)= CXYZ(13)*RXYZ(1)                                     
         CXYZ(24)= CXYZ(14)*RXYZ(1)                                     
         CXYZ(25)= CXYZ(14)*RXYZ(2)                                     
         CXYZ(26)= CXYZ(15)*RXYZ(1)                                     
         CXYZ(27)= CXYZ(16)*RXYZ(1)                                     
         CXYZ(28)= CXYZ(17)*RXYZ(1)                                     
         CXYZ(29)= CXYZ(17)*RXYZ(2)                                     
         CXYZ(30)= CXYZ(18)*RXYZ(1)                                     
         CXYZ(31)= CXYZ(19)*RXYZ(1)                                     
         CXYZ(32)= CXYZ(19)*RXYZ(2)                                     
         CXYZ(33)= CXYZ(20)*RXYZ(1)                                     
         CXYZ(34)= CXYZ(20)*RXYZ(2)                                     
         CXYZ(35)= CXYZ(20)*RXYZ(3)                                     
      IF(N.LE.5) GO TO 100                                              
         CXYZ(36)= CXYZ(21)*RXYZ(1)                                     
         CXYZ(37)= CXYZ(22)*RXYZ(1)                                     
         CXYZ(38)= CXYZ(23)*RXYZ(1)                                     
         CXYZ(39)= CXYZ(24)*RXYZ(1)                                     
         CXYZ(40)= CXYZ(25)*RXYZ(1)                                     
         CXYZ(41)= CXYZ(25)*RXYZ(2)                                     
         CXYZ(42)= CXYZ(26)*RXYZ(1)                                     
         CXYZ(43)= CXYZ(27)*RXYZ(1)                                     
         CXYZ(44)= CXYZ(28)*RXYZ(1)                                     
         CXYZ(45)= CXYZ(29)*RXYZ(1)                                     
         CXYZ(46)= CXYZ(29)*RXYZ(2)                                     
         CXYZ(47)= CXYZ(30)*RXYZ(1)                                     
         CXYZ(48)= CXYZ(31)*RXYZ(1)                                     
         CXYZ(49)= CXYZ(32)*RXYZ(1)                                     
         CXYZ(50)= CXYZ(32)*RXYZ(2)                                     
         CXYZ(51)= CXYZ(33)*RXYZ(1)                                     
         CXYZ(52)= CXYZ(34)*RXYZ(1)                                     
         CXYZ(53)= CXYZ(34)*RXYZ(2)                                     
         CXYZ(54)= CXYZ(35)*RXYZ(1)                                     
         CXYZ(55)= CXYZ(35)*RXYZ(2)                                     
         CXYZ(56)= CXYZ(35)*RXYZ(3)                                     
  100 CONTINUE                                                          
         PHI(  1)=          FT(0)                                       
         PHI(  2)= CXYZ( 2)*FT(1)                                       
         PHI(  3)= CXYZ( 3)*FT(1)                                       
         PHI(  4)= CXYZ( 4)*FT(1)                                       
      IF(N.LE.2) GO TO 200                                              
         PHI(  5)= CXYZ( 1)*FT(1)                                       
         DO I=  6, 11                                                   
            PHI(I)= CXYZ(I- 1)*FT(2)                                    
         END DO                                                         
      IF(N.LE.3) GO TO 200                                              
         PHI( 12)= CXYZ( 2)*FT(2)                                       
         PHI( 13)= CXYZ( 3)*FT(2)                                       
         PHI( 14)= CXYZ( 4)*FT(2)                                       
         DO I= 15, 24                                                   
            PHI(I)= CXYZ(I- 4)*FT(3)                                    
         END DO                                                         
      IF(N.LE.4) GO TO 200                                              
         PHI( 25)= CXYZ( 1)*FT(2)                                       
         DO I= 26, 31                                                   
            PHI(I)= CXYZ(I-21)*FT(3)                                    
         END DO                                                         
         DO I= 32, 46                                                   
            PHI(I)= CXYZ(I-11)*FT(4)                                    
         END DO                                                         
      IF(N.LE.5) GO TO 200                                              
         PHI( 47)= CXYZ( 2)*FT(3)                                       
         PHI( 48)= CXYZ( 3)*FT(3)                                       
         PHI( 49)= CXYZ( 4)*FT(3)                                       
         DO I= 50, 59                                                   
            PHI(I)= CXYZ(I-39)*FT(4)                                    
         END DO                                                         
         DO I= 60, 80                                                   
            PHI(I)= CXYZ(I-24)*FT(5)                                    
         END DO                                                         
  200 CONTINUE                                                          
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2C   *DECK PHIIJ3                                           
C>                                                                      
C>    @brief   ERIC post-contraction phase utility (3)                  
C>                                                                      
C>    @details common post-contraction expressions in ERIC (3)          
C>                                                                      
C>    @author  Jose Sierra, 2013                                        
C>                                                                      
      SUBROUTINE PHIIJ3(I,PHI)                                          
C                                                                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C                                                                       
      DIMENSION  PHI(*)                                                 
C                                                                       
      PARAMETER (F03=3.0D+00)                                           
C                                                                       
      J=I-3                                                             
      PHI(I  )= PHI(I  )-PHI(J  )*F03                                   
      PHI(I+1)= PHI(I+1)-PHI(J+1)                                       
      PHI(I+2)= PHI(I+2)-PHI(J  )                                       
      PHI(I+3)= PHI(I+3)-PHI(J+1)*F03                                   
      PHI(I+4)= PHI(I+4)-PHI(J+2)                                       
      PHI(I+6)= PHI(I+6)-PHI(J+2)                                       
      PHI(I+7)= PHI(I+7)-PHI(J  )                                       
      PHI(I+8)= PHI(I+8)-PHI(J+1)                                       
      PHI(I+9)= PHI(I+9)-PHI(J+2)*F03                                   
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2C   *DECK PHIIJ4                                           
C>                                                                      
C>    @brief   ERIC post-contraction phase utility (4)                  
C>                                                                      
C>    @details common post-contraction expressions in ERIC (4)          
C>                                                                      
C>    @author  Jose Sierra, 2013                                        
C>                                                                      
      SUBROUTINE PHIIJ4(I,PHI)                                          
C                                                                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C                                                                       
      DIMENSION  PHI(*)                                                 
C                                                                       
      PARAMETER (F03=3.0D+00)                                           
C                                                                       
      J=I-6                                                             
      PHI(I   )= PHI(I   )-PHI(J  )*F03                                 
      PHI(I+ 1)= PHI(I+ 1)-PHI(J+1)*F03                                 
      PHI(I+ 2)= PHI(I+ 2)-PHI(J+2)                                     
      PHI(I+ 3)= PHI(I+ 3)-PHI(J+1)*F03                                 
      PHI(I+ 4)= PHI(I+ 4)-PHI(J+2)*F03                                 
      PHI(I+ 5)= PHI(I+ 5)-PHI(J+3)*F03                                 
      PHI(I+ 6)= PHI(I+ 6)-PHI(J+4)                                     
      PHI(I+ 7)= PHI(I+ 7)-PHI(J+3)                                     
      PHI(I+ 8)= PHI(I+ 8)-PHI(J+4)*F03                                 
      PHI(I+ 9)= PHI(I+ 9)-PHI(J+5)                                     
      PHI(I+10)= PHI(I+10)-PHI(J+1)                                     
      PHI(I+11)= PHI(I+11)-PHI(J+5)                                     
      PHI(I+12)= PHI(I+12)-PHI(J+3)*F03                                 
      PHI(I+13)= PHI(I+13)-PHI(J+4)*F03                                 
      PHI(I+14)= PHI(I+14)-PHI(J+5)*F03                                 
C                                                                       
      PHI(J   )= PHI(J   )-PHI(J-1)                                     
      PHI(J+ 2)= PHI(J +2)-PHI(J-1)                                     
      PHI(J+ 5)= PHI(J +5)-PHI(J-1)                                     
C                                                                       
      PHI(I   )= PHI(I   )-PHI(J  )*F03                                 
      PHI(I+ 2)= PHI(I+ 2)-PHI(J  )                                     
      PHI(I+ 4)= PHI(I+ 4)-PHI(J+2)*F03                                 
      PHI(I+ 9)= PHI(I+ 9)-PHI(J  )                                     
      PHI(I+11)= PHI(I+11)-PHI(J+2)                                     
      PHI(I+14)= PHI(I+14)-PHI(J+5)*F03                                 
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2C   *DECK PHIIJ5                                           
C>                                                                      
C>    @brief   ERIC post-contraction phase utility (5)                  
C>                                                                      
C>    @details common post-contraction expressions in ERIC (5)          
C>                                                                      
C>    @author  Jose Sierra, 2013                                        
C>                                                                      
      SUBROUTINE PHIIJ5(I,PHI)                                          
C                                                                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C                                                                       
      DIMENSION  PHI(*)                                                 
C                                                                       
      PARAMETER (F02=2.0D+00)                                           
      PARAMETER (F03=3.0D+00)                                           
      PARAMETER (F04=4.0D+00)                                           
C                                                                       
      J=I-10                                                            
      PHI(I   )= PHI(I   )-PHI(J  )*F03                                 
      PHI(I+ 1)= PHI(I+ 1)-PHI(J+1)*F03                                 
      PHI(I+ 2)= PHI(I+ 2)-PHI(J+2)*F03                                 
      PHI(I+ 3)= PHI(I+ 3)-PHI(J+3)                                     
      PHI(I+ 4)= PHI(I+ 4)-PHI(J+2)*F03                                 
      PHI(I+ 5)= PHI(I+ 5)-PHI(J+3)*F03                                 
      PHI(I+ 6)= PHI(I+ 6)-PHI(J+4)*F03                                 
      PHI(I+ 7)= PHI(I+ 7)-PHI(J+5)*F03                                 
      PHI(I+ 8)= PHI(I+ 8)-PHI(J+6)                                     
      PHI(I+ 9)= PHI(I+ 9)-PHI(J+5)*F03                                 
      PHI(I+10)= PHI(I+10)-PHI(J+6)*F03                                 
      PHI(I+11)= PHI(I+11)-PHI(J+7)*F03                                 
      PHI(I+12)= PHI(I+12)-PHI(J+8)                                     
      PHI(I+13)= PHI(I+13)-PHI(J+7)                                     
      PHI(I+14)= PHI(I+14)-PHI(J+8)*F03                                 
      PHI(I+15)= PHI(I+15)-PHI(J+9)                                     
      PHI(I+16)= PHI(I+16)-PHI(J+5)*F03                                 
      PHI(I+17)= PHI(I+17)-PHI(J+9)                                     
      PHI(I+18)= PHI(I+18)-PHI(J+7)*F03                                 
      PHI(I+19)= PHI(I+19)-PHI(J+8)*F03                                 
      PHI(I+20)= PHI(I+20)-PHI(J+9)*F03                                 
C                                                                       
      PHI(J   )= PHI(J   )-PHI(J-3)                                     
      PHI(J+ 1)= PHI(J+ 1)-PHI(J-2)                                     
      PHI(J+ 3)= PHI(J+ 3)-PHI(J-2)                                     
      PHI(J+ 4)= PHI(J+ 4)-PHI(J-1)                                     
      PHI(J+ 6)= PHI(J+ 6)-PHI(J-1)                                     
      PHI(J+ 9)= PHI(J+ 9)-PHI(J-1)                                     
C                                                                       
      PHI(I   )= PHI(I   )-PHI(J  )*F03                                 
      PHI(I+ 1)= PHI(I+ 1)-PHI(J+1)*F03                                 
      PHI(I+ 3)= PHI(I+ 3)-PHI(J+1)*F03                                 
      PHI(I+ 5)= PHI(I+ 5)-PHI(J+3)*F03                                 
      PHI(I+ 6)= PHI(I+ 6)-PHI(J+4)*F03                                 
      PHI(I+ 8)= PHI(I+ 8)-PHI(J+4)                                     
      PHI(I+10)= PHI(I+10)-PHI(J+6)*F03                                 
      PHI(I+12)= PHI(I+12)-PHI(J+1)                                     
      PHI(I+15)= PHI(I+15)-PHI(J+4)*F03                                 
      PHI(I+17)= PHI(I+17)-PHI(J+6)*F03                                 
      PHI(I+20)= PHI(I+20)-PHI(J+9)*F03                                 
C                                                                       
      PHI(J   )= PHI(J   )-PHI(J-3)*F02                                 
      PHI(J+ 2)= PHI(J+ 2)-PHI(J-3)                                     
      PHI(J+ 3)= PHI(J+ 3)-PHI(J-2)*F02                                 
      PHI(J+ 7)= PHI(J+ 7)-PHI(J-3)                                     
      PHI(J+ 8)= PHI(J+ 8)-PHI(J-2)                                     
      PHI(J+ 9)= PHI(J+ 9)-PHI(J-1)*F02                                 
C                                                                       
      PHI(I   )= PHI(I   )-PHI(J  )*F04                                 
      PHI(I+ 2)= PHI(I+ 2)-PHI(J  )                                     
      PHI(I+ 4)= PHI(I+ 4)-PHI(J+2)*F03                                 
      PHI(I+ 5)= PHI(I+ 5)-PHI(J+3)*F04                                 
      PHI(I+11)= PHI(I+11)-PHI(J  )                                     
      PHI(I+13)= PHI(I+13)-PHI(J+2)                                     
      PHI(I+14)= PHI(I+14)-PHI(J+3)                                     
      PHI(I+18)= PHI(I+18)-PHI(J+7)*F03                                 
      PHI(I+19)= PHI(I+19)-PHI(J+8)*F03                                 
      PHI(I+20)= PHI(I+20)-PHI(J+9)*F04                                 
C                                                                       
      RETURN                                                            
      END                                                               
