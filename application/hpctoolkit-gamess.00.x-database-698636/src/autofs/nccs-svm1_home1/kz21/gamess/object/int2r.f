C  6 Jun 18 - DGF - tweaks for FMO 5.3                                  
C  8 Aug 13 - JMS - extensive reorganization                            
C 22 DEC 06 - ST,NK,MC - ADD LC EXCHANGE ARGUMENTS                      
C  5 JUL 05 - MWS - SELECT NEW ATOM,BASIS,EFP,PCM,DAF DIMENSIONS        
C 13 FEB 05 - MWS - PAD COMMON BLOCK NSHEL                              
C  5 FEB 05 - MWS - CHANGE COMMON NAME GAMMAF TO BE FMTTBL INSTEAD      
C 10 NOV 04 - KI  - IMPLEMENT S,P,D,L ROTATED AXIS INTEGRAL CODE        
C                                                                       
C*MODULE INT2R   *DECK GENR03                                           
C>                                                                      
C>    @brief   rotated axis integration involving s,p,L,d shells        
C>                                                                      
C>    @details rotated axis integration involving s,p,L,d shells,       
C>             by a mix of rotated axis and McMurchie/Davidson.         
C>                       K.Ishimura, S.Nagase                           
C>                Theoret.Chem.Acc. 120, 185-189(2008)                  
C>                                                                      
C>    @author  Kazuya Ishimura, at the Institute for Molecular Science, 
C>             sponsored by Naregi Nano Science project, in 2004.       
C>             Extensively revised by Jose Sierra in 2013.              
C>                                                                      
      SUBROUTINE GENR03(GHONDO)                                         
      use mx_limits, only: mxgtot,mxsh,mxgsh,mxg2                       
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      DIMENSION GHONDO(*)                                               
C                                                                       
C                                                                       
      COMMON /B     / CO(MXSH,3)                                        
      COMMON /FLIPS / IB(4,3)                                           
C$omp threadprivate(/FLIPS/)
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS,                      
     2                TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2),      
     3                D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2),      
     4                NGANGB                                            
C$omp threadprivate(/GEOMPQ/)
      COMMON /INTAC2/ EI1,EI2,CUX                                       
      COMMON /MAXC  / CMAX(MXGTOT),CMAXA(MXGSH),CMAXB(MXGSH),           
     2                CMAXC(MXGSH),CMAXD(MXGSH),ISMLP(MXG2),ISMLQ       
C$omp threadprivate(/MAXC/)
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),      
     2                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),      
     3                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),   
     4                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL           
      COMMON /SHLG70/ ISHELL,JSHELL,KSHELL,LSHELL,INEW,JNEW,KNEW,LNEW   
C$omp threadprivate(/SHLG70/)
      COMMON /SHLLFO/ NGA,LA,EXA(MXGSH),CSA(MXGSH),CPA(MXGSH),          
     2                NGB,LB,EXB(MXGSH),CSB(MXGSH),CPB(MXGSH),          
     3                NGC,LC,EXC(MXGSH),CSC(MXGSH),CPC(MXGSH),          
     4                NGD,LD,EXD(MXGSH),CSD(MXGSH),CPD(MXGSH)           
C$omp threadprivate(/SHLLFO/)
      COMMON /SHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,          
     2                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,          
     3                NIJ,IJ,KL,IJKL                                    
C$omp threadprivate(/SHLNOS/)
      COMMON /SHLSPD/ CDA(MXGSH),CDB(MXGSH),CDC(MXGSH),CDD(MXGSH),      
     2                D02D(MXG2),D12D(MXG2),D22D(MXG2)                  
C$omp threadprivate(/SHLSPD/)
C                                                                       
C     LABELLED COMMON JMSGYH DEFINED FOR COMPUTATIONAL EFFICIENCY.      
C     IT IS ONLY USED IN THIS MODULE INT2R AND IN MODULE INT2B.         
C                                                                       
      COMMON /JMSGYH/ SQ(4)                                             
C$omp threadprivate(/JMSGYH/)
      COMMON /KI2 / ACY,ACY2,AQX,AQX2,AQXY,Y03,Y04                      
C$omp threadprivate(/KI2/)
C                                                                       
      DIMENSION P12(3,3),P34(3,3),P(3,3),T(3)                           
C                                                                       
      PARAMETER (ZER=0.0D+00)                                           
      PARAMETER (PT5=0.5D+00)                                           
      PARAMETER (ONE=1.0D+00)                                           
      PARAMETER (PT7=0.7D+00)                                           
      PARAMETER (PT9=0.9D+00)                                           
      PARAMETER (ACYCUT=1.0D-10)                                        
      PARAMETER (TENM12=1.0D-12)                                        
C                                                                       
C                PITO52=(PI +PI )*PI * SQRT(PI )                        
      PARAMETER (PITO52=34.98683665524973D+00)                          
C                                                                       
      PIF= PITO52                                                       
C                                                                       
      IEXCH= 1                                                          
      LAT= KTYPE(ISHELL)-1                                              
      LBT= KTYPE(JSHELL)-1                                              
      LCT= KTYPE(KSHELL)-1                                              
      LDT= KTYPE(LSHELL)-1                                              
      ITYPE= 1+LDT+3*(LCT+3*(LBT+3*LAT))                                
      LTYPE= LAT+LBT+LCT+LDT                                            
C                                                                       
      IF(LTYPE.EQ.2) THEN                                               
         JTYPE= 7                                                       
         IF(ITYPE.EQ.3) THEN                                            
            GO TO 1010                                                  
         ELSEIF(ITYPE.EQ. 7) THEN                                       
            GO TO 1030                                                  
         ELSEIF(ITYPE.EQ.19) THEN                                       
            GO TO 1040                                                  
         ELSEIF(ITYPE.EQ.55) THEN                                       
            GO TO 1060                                                  
         ENDIF                                                          
C                                                                       
      ELSEIF(LTYPE.EQ.3) THEN                                           
         IF(ITYPE.EQ. 6) THEN                                           
            JTYPE= 8                                                    
            GO TO 1010                                                  
         ELSEIF(ITYPE.EQ. 8) THEN                                       
            JTYPE= 8                                                    
            GO TO 1030                                                  
         ELSEIF(ITYPE.EQ.46) THEN                                       
            JTYPE= 8                                                    
            GO TO 1040                                                  
         ELSEIF(ITYPE.EQ.64) THEN                                       
            JTYPE= 8                                                    
            GO TO 1060                                                  
         ELSEIF(ITYPE.EQ.12) THEN                                       
            JTYPE= 9                                                    
            GO TO 1010                                                  
         ELSEIF(ITYPE.EQ.30) THEN                                       
            JTYPE= 9                                                    
            GO TO 1020                                                  
         ELSEIF(ITYPE.EQ.16) THEN                                       
            JTYPE= 9                                                    
            GO TO 1030                                                  
         ELSEIF(ITYPE.EQ.20) THEN                                       
            JTYPE= 9                                                    
            GO TO 1040                                                  
         ELSEIF(ITYPE.EQ.34) THEN                                       
            JTYPE= 9                                                    
            GO TO 1050                                                  
         ELSEIF(ITYPE.EQ.56) THEN                                       
            JTYPE= 9                                                    
            GO TO 1060                                                  
         ELSEIF(ITYPE.EQ.22) THEN                                       
            JTYPE= 9                                                    
            GO TO 1070                                                  
         ELSEIF(ITYPE.EQ.58) THEN                                       
            JTYPE= 9                                                    
            GO TO 1080                                                  
         ENDIF                                                          
C                                                                       
      ELSEIF(LTYPE.EQ.4) THEN                                           
         IF(ITYPE.EQ. 9) THEN                                           
            JTYPE= 10                                                   
            GO TO 1010                                                  
         ELSEIF(ITYPE.EQ.73) THEN                                       
            JTYPE= 10                                                   
            GO TO 1040                                                  
         ELSEIF(ITYPE.EQ.15) THEN                                       
            JTYPE= 11                                                   
            GO TO 1010                                                  
         ELSEIF(ITYPE.EQ.33) THEN                                       
            JTYPE= 11                                                   
            GO TO 1020                                                  
         ELSEIF(ITYPE.EQ.17) THEN                                       
            JTYPE= 11                                                   
            GO TO 1030                                                  
         ELSEIF(ITYPE.EQ.47) THEN                                       
            JTYPE= 11                                                   
            GO TO 1040                                                  
         ELSEIF(ITYPE.EQ.35) THEN                                       
            JTYPE= 11                                                   
            GO TO 1050                                                  
         ELSEIF(ITYPE.EQ.65) THEN                                       
            JTYPE= 11                                                   
            GO TO 1060                                                  
         ELSEIF(ITYPE.EQ.49) THEN                                       
            JTYPE= 11                                                   
            GO TO 1070                                                  
         ELSEIF(ITYPE.EQ.67) THEN                                       
            JTYPE= 11                                                   
            GO TO 1080                                                  
         ELSEIF(ITYPE.EQ.21) THEN                                       
            JTYPE= 12                                                   
            GO TO 1010                                                  
         ELSEIF(ITYPE.EQ.57) THEN                                       
            JTYPE= 12                                                   
            GO TO 1020                                                  
         ELSEIF(ITYPE.EQ.25) THEN                                       
            JTYPE= 12                                                   
            GO TO 1030                                                  
         ELSEIF(ITYPE.EQ.61) THEN                                       
            JTYPE= 12                                                   
            GO TO 1050                                                  
         ELSEIF(ITYPE.EQ.39) THEN                                       
            JTYPE= 13                                                   
            GO TO 1010                                                  
         ELSEIF(ITYPE.EQ.43) THEN                                       
            JTYPE= 13                                                   
            GO TO 1030                                                  
         ELSEIF(ITYPE.EQ.23) THEN                                       
            JTYPE= 13                                                   
            GO TO 1040                                                  
         ELSEIF(ITYPE.EQ.59) THEN                                       
            JTYPE= 13                                                   
            GO TO 1060                                                  
         ENDIF                                                          
C                                                                       
      ELSEIF(LTYPE.EQ.5) THEN                                           
         IF(ITYPE.EQ.18) THEN                                           
            JTYPE= 14                                                   
            GO TO 1010                                                  
         ELSEIF(ITYPE.EQ.36) THEN                                       
            JTYPE= 14                                                   
            GO TO 1020                                                  
         ELSEIF(ITYPE.EQ.74) THEN                                       
            JTYPE= 14                                                   
            GO TO 1040                                                  
         ELSEIF(ITYPE.EQ.76) THEN                                       
            JTYPE= 14                                                   
            GO TO 1070                                                  
         ELSEIF(ITYPE.EQ.24) THEN                                       
            JTYPE= 15                                                   
            GO TO 1010                                                  
         ELSEIF(ITYPE.EQ.60) THEN                                       
            JTYPE= 15                                                   
            GO TO 1020                                                  
         ELSEIF(ITYPE.EQ.26) THEN                                       
            JTYPE= 15                                                   
            GO TO 1030                                                  
         ELSEIF(ITYPE.EQ.48) THEN                                       
            JTYPE= 15                                                   
            GO TO 1040                                                  
         ELSEIF(ITYPE.EQ.62) THEN                                       
            JTYPE= 15                                                   
            GO TO 1050                                                  
         ELSEIF(ITYPE.EQ.66) THEN                                       
            JTYPE= 15                                                   
            GO TO 1060                                                  
         ELSEIF(ITYPE.EQ.52) THEN                                       
            JTYPE= 15                                                   
            GO TO 1070                                                  
         ELSEIF(ITYPE.EQ.70) THEN                                       
            JTYPE= 15                                                   
            GO TO 1080                                                  
         ELSEIF(ITYPE.EQ.42) THEN                                       
            JTYPE= 16                                                   
            GO TO 1010                                                  
         ELSEIF(ITYPE.EQ.44) THEN                                       
            JTYPE= 16                                                   
            GO TO 1030                                                  
         ELSEIF(ITYPE.EQ.50) THEN                                       
            JTYPE= 16                                                   
            GO TO 1040                                                  
         ELSEIF(ITYPE.EQ.68) THEN                                       
            JTYPE= 16                                                   
            GO TO 1060                                                  
         ENDIF                                                          
C                                                                       
      ELSEIF(LTYPE.EQ.6) THEN                                           
         IF(ITYPE.EQ.27) THEN                                           
            JTYPE= 17                                                   
            GO TO 1010                                                  
         ELSEIF(ITYPE.EQ.63) THEN                                       
            JTYPE= 17                                                   
            GO TO 1020                                                  
         ELSEIF(ITYPE.EQ.75) THEN                                       
            JTYPE= 17                                                   
            GO TO 1040                                                  
         ELSEIF(ITYPE.EQ.79) THEN                                       
            JTYPE= 17                                                   
            GO TO 1070                                                  
         ELSEIF(ITYPE.EQ.45) THEN                                       
            JTYPE= 18                                                   
            GO TO 1010                                                  
         ELSEIF(ITYPE.EQ.77) THEN                                       
            JTYPE= 18                                                   
            GO TO 1040                                                  
         ELSEIF(ITYPE.EQ.51) THEN                                       
            JTYPE= 19                                                   
            GO TO 1010                                                  
         ELSEIF(ITYPE.EQ.69) THEN                                       
            JTYPE= 19                                                   
            GO TO 1020                                                  
         ELSEIF(ITYPE.EQ.53) THEN                                       
            JTYPE= 19                                                   
            GO TO 1030                                                  
         ELSEIF(ITYPE.EQ.71) THEN                                       
            JTYPE= 19                                                   
            GO TO 1050                                                  
         ENDIF                                                          
C                                                                       
      ELSEIF(LTYPE.EQ.7) THEN                                           
         JTYPE= 20                                                      
         IF(ITYPE.EQ.54) THEN                                           
            GO TO 1010                                                  
         ELSEIF(ITYPE.EQ.72) THEN                                       
            GO TO 1020                                                  
         ELSEIF(ITYPE.EQ.78) THEN                                       
            GO TO 1040                                                  
         ELSEIF(ITYPE.EQ.80) THEN                                       
            GO TO 1070                                                  
         ENDIF                                                          
C                                                                       
      ELSEIF(LTYPE.EQ.8) THEN                                           
         JTYPE= 21                                                      
         GO TO 1010                                                     
      ENDIF                                                             
C                                                                       
 1010 INEW= ISHELL                                                      
      JNEW= JSHELL                                                      
      KNEW= KSHELL                                                      
      LNEW= LSHELL                                                      
      LA= LAT                                                           
      LB= LBT                                                           
      LC= LCT                                                           
      LD= LDT                                                           
      IB(1,IEXCH)= 1                                                    
      IB(2,IEXCH)= 2                                                    
      IB(3,IEXCH)= 3                                                    
      IB(4,IEXCH)= 4                                                    
      GO TO 1090                                                        
C                                                                       
 1020 INEW= JSHELL                                                      
      JNEW= ISHELL                                                      
      KNEW= KSHELL                                                      
      LNEW= LSHELL                                                      
      LA= LBT                                                           
      LB= LAT                                                           
      LC= LCT                                                           
      LD= LDT                                                           
      IB(1,IEXCH)= 2                                                    
      IB(2,IEXCH)= 1                                                    
      IB(3,IEXCH)= 3                                                    
      IB(4,IEXCH)= 4                                                    
      GO TO 1090                                                        
C                                                                       
 1030 INEW= ISHELL                                                      
      JNEW= JSHELL                                                      
      KNEW= LSHELL                                                      
      LNEW= KSHELL                                                      
      LA= LAT                                                           
      LB= LBT                                                           
      LC= LDT                                                           
      LD= LCT                                                           
      IB(1,IEXCH)= 1                                                    
      IB(2,IEXCH)= 2                                                    
      IB(3,IEXCH)= 4                                                    
      IB(4,IEXCH)= 3                                                    
      GO TO 1090                                                        
C                                                                       
 1040 INEW= KSHELL                                                      
      JNEW= LSHELL                                                      
      KNEW= ISHELL                                                      
      LNEW= JSHELL                                                      
      LA= LCT                                                           
      LB= LDT                                                           
      LC= LAT                                                           
      LD= LBT                                                           
      IB(1,IEXCH)= 3                                                    
      IB(2,IEXCH)= 4                                                    
      IB(3,IEXCH)= 1                                                    
      IB(4,IEXCH)= 2                                                    
      GO TO 1090                                                        
C                                                                       
 1050 INEW= JSHELL                                                      
      JNEW= ISHELL                                                      
      KNEW= LSHELL                                                      
      LNEW= KSHELL                                                      
      LA= LBT                                                           
      LB= LAT                                                           
      LC= LDT                                                           
      LD= LCT                                                           
      IB(1,IEXCH)= 2                                                    
      IB(2,IEXCH)= 1                                                    
      IB(3,IEXCH)= 4                                                    
      IB(4,IEXCH)= 3                                                    
      GO TO 1090                                                        
C                                                                       
 1060 INEW= KSHELL                                                      
      JNEW= LSHELL                                                      
      KNEW= JSHELL                                                      
      LNEW= ISHELL                                                      
      LA= LCT                                                           
      LB= LDT                                                           
      LC= LBT                                                           
      LD= LAT                                                           
      IB(1,IEXCH)= 4                                                    
      IB(2,IEXCH)= 3                                                    
      IB(3,IEXCH)= 1                                                    
      IB(4,IEXCH)= 2                                                    
      GO TO 1090                                                        
C                                                                       
 1070 INEW= LSHELL                                                      
      JNEW= KSHELL                                                      
      KNEW= ISHELL                                                      
      LNEW= JSHELL                                                      
      LA= LDT                                                           
      LB= LCT                                                           
      LC= LAT                                                           
      LD= LBT                                                           
      IB(1,IEXCH)= 3                                                    
      IB(2,IEXCH)= 4                                                    
      IB(3,IEXCH)= 2                                                    
      IB(4,IEXCH)= 1                                                    
      GO TO 1090                                                        
C                                                                       
 1080 INEW= LSHELL                                                      
      JNEW= KSHELL                                                      
      KNEW= JSHELL                                                      
      LNEW= ISHELL                                                      
      LA= LDT                                                           
      LB= LCT                                                           
      LC= LBT                                                           
      LD= LAT                                                           
      IB(1,IEXCH)= 4                                                    
      IB(2,IEXCH)= 3                                                    
      IB(3,IEXCH)= 2                                                    
      IB(4,IEXCH)= 1                                                    
 1090 CONTINUE                                                          
C                                                                       
C EMPTY INTEGRAL SUMMATION STORAGE                                      
C                                                                       
      IKL= 0                                                            
      IF(JTYPE.EQ. 7) THEN                                              
         CALL INTK07(IKL)                                               
      ELSEIF(JTYPE.EQ. 8) THEN                                          
         CALL INTK08(IKL)                                               
      ELSEIF(JTYPE.EQ. 9) THEN                                          
         CALL INTK09(IKL)                                               
      ELSEIF(JTYPE.EQ.10) THEN                                          
         CALL INTK10(IKL)                                               
      ELSEIF(JTYPE.EQ.11) THEN                                          
         CALL INTK11(IKL)                                               
      ELSEIF(JTYPE.EQ.12) THEN                                          
         CALL INTK12(IKL)                                               
      ELSEIF(JTYPE.EQ.13) THEN                                          
         CALL INTK13(IKL)                                               
      ELSEIF(JTYPE.EQ.14) THEN                                          
         CALL INTK14(IKL)                                               
      ELSEIF(JTYPE.EQ.15) THEN                                          
         CALL INTK15(IKL)                                               
      ELSEIF(JTYPE.EQ.16) THEN                                          
         CALL INTK16(IKL)                                               
      ELSEIF(JTYPE.EQ.17) THEN                                          
         CALL INTK17(IKL)                                               
      ELSEIF(JTYPE.EQ.18) THEN                                          
         CALL INTK18(IKL)                                               
      ELSEIF(JTYPE.EQ.19) THEN                                          
         CALL INTK19(IKL)                                               
      ELSEIF(JTYPE.EQ.20) THEN                                          
         CALL INTK20(IKL)                                               
      ELSEIF(JTYPE.EQ.21) THEN                                          
         CALL INTK21(IKL)                                               
      ENDIF                                                             
C                                                                       
C OBTAIN INFORMATION ABOUT SHELLS: INEW, KNEW, JNEW, LNEW               
C NUMBER OF GAUSSIANS GO INTO NGA,... IN COMMON SHLLFO                  
C SHELL ANGULAR QUANTUM NUMBERS LA,... GO INTO COMMON SHLLFO            
C GAUSSIAN EXPONENTS GO INTO ARRAYS EXA,EXB,EXC,EXD IN COMMON SHLLFO    
C GAUSSIAN COEFFICIENTS GO INTO ARRAYS CSA,CPA,... IN COMMON SHLLFO     
C                                                                       
C NUMBERS OF GAUSSIAN FUNCTIONS IN SHELLS INEW JNEW KNEW AND LNEW       
C                                                                       
      NGA= KNG(INEW)                                                    
      NGB= KNG(JNEW)                                                    
      NGC= KNG(KNEW)                                                    
      NGD= KNG(LNEW)                                                    
C                                                                       
C STARTING LOCATIONS OF SHELLS INEW JNEW KNEW AND LNEW IN LIST          
C OF GAUSSIAN FUNCTIONS                                                 
C                                                                       
      I= KSTART(INEW)-1                                                 
      J= KSTART(JNEW)-1                                                 
      K= KSTART(KNEW)-1                                                 
      L= KSTART(LNEW)-1                                                 
C                                                                       
      MINI = KMIN(INEW)                                                 
      MAXI = KMAX(INEW)                                                 
      LOCI = KLOC(INEW)-MINI                                            
      MINJ = KMIN(JNEW)                                                 
      MAXJ = KMAX(JNEW)                                                 
      LOCJ = KLOC(JNEW)-MINJ                                            
      MINK = KMIN(KNEW)                                                 
      MAXK = KMAX(KNEW)                                                 
      LOCK = KLOC(KNEW)-MINK                                            
      MINL = KMIN(LNEW)                                                 
      MAXL = KMAX(LNEW)                                                 
      LOCL = KLOC(LNEW)-MINL                                            
C                                                                       
C LOOP OVER GAUSSIANS IN EACH SHELL                                     
C FIRST SHELL INEW                                                      
C                                                                       
      DO NI=1,NGA                                                       
         N=I+NI                                                         
         CMAXA(NI)= CMAX(N)                                             
         EXA(NI)= EX(N)                                                 
         CSA(NI)= CS(N)                                                 
         CPA(NI)= CP(N)                                                 
         CDA(NI)= CD(N)                                                 
      ENDDO                                                             
C                                                                       
      DO NJ=1,NGB                                                       
         N=J+NJ                                                         
         CMAXB(NJ)= CMAX(N)                                             
         EXB(NJ)= EX(N)                                                 
         CSB(NJ)= CS(N)                                                 
         CPB(NJ)= CP(N)                                                 
         CDB(NJ)= CD(N)                                                 
      ENDDO                                                             
C                                                                       
      DO NK=1,NGC                                                       
         N=K+NK                                                         
         CMAXC(NK)= CMAX(N)*QQ4                                         
         EXC(NK)= EX(N)                                                 
         CSC(NK)= CS(N)*QQ4                                             
         CPC(NK)= CP(N)*QQ4                                             
         CDC(NK)= CD(N)*QQ4                                             
      ENDDO                                                             
C                                                                       
      DO NL=1,NGD                                                       
         N=L+NL                                                         
         CMAXD(NL)= CMAX(N)                                             
         EXD(NL)= EX(N)                                                 
         CDD(NL)= CD(N)                                                 
      ENDDO                                                             
C                                                                       
      NGANGB=NGA*NGB                                                    
C                                                                       
C COORDINATES OF ATOMS ASSOCIATED WITH SHELLS INEW JNEW KNEW AND LNEW   
C                                                                       
      R12= ZER                                                          
      R34= ZER                                                          
      DO 150 N=1,3                                                      
         P12(N,1)= CO(INEW,N)                                           
         P12(N,2)= CO(JNEW,N)                                           
         P12(N,3)= P12(N,2)-P12(N,1)                                    
      R12= R12+P12(N,3)*P12(N,3)                                        
         P34(N,1)= CO(KNEW,N)                                           
         P34(N,2)= CO(LNEW,N)                                           
         P34(N,3)= P34(N,2)-P34(N,1)                                    
  150 R34= R34+P34(N,3)*P34(N,3)                                        
C                                                                       
C FIND DIRECTION COSINES OF PENULTIMATE AXES FROM COORDINATES OF AB     
C P(1,1),P(1,2),... ARE DIRECTION COSINES OF AXES AT P.  Z-AXIS ALONG AB
C T(1),T(2),T(3)... ARE DIRECTION COSINES OF AXES AT Q.  Z-AXIS ALONG CD
C                                                                       
C FIND DIRECTION COSINES OF AB AND CD. THESE ARE LOCAL Z-AXES.          
C IF INDETERMINATE TAKE ALONG SPACE Z-AXIS                              
C                                                                       
      P(1,3)= ZER                                                       
      P(2,3)= ZER                                                       
      P(3,3)= ONE                                                       
      RAB= ZER                                                          
      IF(R12.NE.ZER) THEN                                               
         RAB= SQRT(R12)                                                 
         TMP= ONE/RAB                                                   
         P(1,3)= P12(1,3)*TMP                                           
         P(2,3)= P12(2,3)*TMP                                           
         P(3,3)= P12(3,3)*TMP                                           
      ENDIF                                                             
C                                                                       
      T(1)= ZER                                                         
      T(2)= ZER                                                         
      T(3)= ONE                                                         
      RCD= ZER                                                          
      IF(R34.NE.ZER) THEN                                               
         RCD= SQRT(R34)                                                 
         TMP= ONE/RCD                                                   
         T(1)= P34(1,3)*TMP                                             
         T(2)= P34(2,3)*TMP                                             
         T(3)= P34(3,3)*TMP                                             
      ENDIF                                                             
C                                                                       
C FIND LOCAL Y-AXIS AS COMMON PERPENDICULAR TO AB AND CD                
C IF INDETERMINATE TAKE PERPENDICULAR TO AB AND SPACE Z-AXIS            
C IF STILL INDETERMINATE TAKE PERPENDICULAR TO AB AND SPACE X-AXIS      
C                                                                       
      COSG= T(1)*P(1,3)+T(2)*P(2,3)+T(3)*P(3,3)                         
      COSG= MIN( ONE,COSG)                                              
      COSG= MAX(-ONE,COSG)                                              
C     SING= SQRT(ONE-COSG*COSG)                                         
C                                                                       
C MODIFIED ROTATION TESTING.                                            
C THIS FIX CURES THE SMALL ANGLE PROBLEM.                               
C                                                                       
      P(1,2)= T(3)*P(2,3)-T(2)*P(3,3)                                   
      P(2,2)= T(1)*P(3,3)-T(3)*P(1,3)                                   
      P(3,2)= T(2)*P(1,3)-T(1)*P(2,3)                                   
      IF( ABS(COSG).GT.PT9) THEN                                        
         SING= SQRT(P(1,2)*P(1,2)+P(2,2)*P(2,2)+P(3,2)*P(3,2))          
      ELSE                                                              
         SING= SQRT(ONE-COSG*COSG)                                      
      ENDIF                                                             
      IF( ABS(COSG).LE.PT9 .OR. SING.GE.TENM12) THEN                    
         TMP= ONE/SING                                                  
         P(1,2)= P(1,2)*TMP                                             
         P(2,2)= P(2,2)*TMP                                             
         P(3,2)= P(3,2)*TMP                                             
      ELSE                                                              
         I=3                                                            
         IF( ABS(P(1,3)).LE.PT7) I=1                                    
         TMP = P(I,3)*P(I,3)                                            
         TMP = MIN( ONE,TMP)                                            
         TMP = SQRT(ONE-TMP)                                            
         IF(TMP.NE.ZER) TMP= ONE/TMP                                    
         IF( ABS(P(1,3)).LE.PT7) THEN                                   
            P(1,2)= ZER                                                 
            P(2,2)= P(3,3)*TMP                                          
            P(3,2)=-P(2,3)*TMP                                          
         ELSE                                                           
            P(1,2)= P(2,3)*TMP                                          
            P(2,2)=-P(1,3)*TMP                                          
            P(3,2)= ZER                                                 
         ENDIF                                                          
      ENDIF                                                             
C                                                                       
C FIND DIRECTION COSINES OF LOCAL X-AXES                                
C                                                                       
      P(1,1)= P(2,2)*P(3,3)-P(3,2)*P(2,3)                               
      P(2,1)= P(3,2)*P(1,3)-P(1,2)*P(3,3)                               
      P(3,1)= P(1,2)*P(2,3)-P(2,2)*P(1,3)                               
C                                                                       
C FIND COORDINATES OF C RELATIVE TO LOCAL AXES AT A                     
C                                                                       
      T(1)= P34(1,1)-P12(1,1)                                           
      T(2)= P34(2,1)-P12(2,1)                                           
      T(3)= P34(3,1)-P12(3,1)                                           
      ACX = T(1)*P(1,1)+T(2)*P(2,1)+T(3)*P(3,1)                         
      ACY = T(1)*P(1,2)+T(2)*P(2,2)+T(3)*P(3,2)                         
      ACZ = T(1)*P(1,3)+T(2)*P(2,3)+T(3)*P(3,3)                         
C                                                                       
C SET ACY= 0  IF CLOSE                                                  
C                                                                       
      IF( ABS(ACY).LE.ACYCUT) THEN                                      
         ACY = ZER                                                      
         ACY2= ZER                                                      
      ELSE                                                              
         ACY2= ACY*ACY                                                  
      ENDIF                                                             
C                                                                       
C DIRECTION COSINES OF CD LOCAL AXES WITH RESPECT TO AB LOCAL AXES      
C (COSG,0,-SING)  (0,1,0)  (SING,0,COSG)                                
C                                                                       
C PRELIMINARY P LOOP                                                    
C                                                                       
C FILL GEOMPQ WITH INFORMATION ABOUT P IN PRELIMINARY P-LOOP            
C                                                                       
      JI= 1                                                             
      DO 170 I=1,NGA                                                    
         X01= EXA(I)                                                    
         DO 170 J=1,NGB                                                 
            X02= EXB(J)                                                 
            X12= X01+X02                                                
            X21= ONE/X12                                                
            Y01= X01*X21                                                
            Y02= ONE-Y01                                                
            Y12= Y01*X02                                                
            TX12(JI)= X12                                               
            TX21(JI)= X21*PT5                                           
            TY02(JI)= Y02*RAB                                           
            TY01(JI)= TY02(JI)-RAB                                      
            R12Y12= R12*Y12                                             
            IF(R12Y12.GT.CUX) THEN                                      
               ISMLP(JI)=2                                              
               GO TO 160                                                
            ENDIF                                                       
            E12= X21* EXP(-R12Y12)                                      
            TST= E12*CMAXA(I)*CMAXB(J)                                  
            ISMLP(JI)=0                                                 
            IF(TST.LE.EI1) ISMLP(JI)=1                                  
            IF(TST.LE.EI2) ISMLP(JI)=2                                  
            E12= PIF*E12                                                
C                                                                       
C FOR TYPES 0000,0001,0011 ONLY D00P NEEDED                             
C                                                                       
            IF(LB.EQ.2) THEN                                            
               D02D(JI)= E12*CSA(I)*CDB(J)                              
               IF(LA.EQ.1) THEN                                         
                  D12D(JI)= E12*CPA(I)*CDB(J)                           
                  IF(D12D(JI).NE.ZER) D02D(JI)=D02D(JI)/D12D(JI)        
               ELSEIF(LA.EQ.2) THEN                                     
                  D22D(JI)= E12*CDA(I)*CDB(J)                           
               ENDIF                                                    
            ELSEIF(LB.EQ.1) THEN                                        
               D00P(JI)= E12*CSA(I)*CSB(J)                              
               D01P(JI)= E12*CSA(I)*CPB(J)                              
               IF(LA.EQ.0) THEN                                         
                  IF(D01P(JI).NE.ZER) D00P(JI)= D00P(JI)/D01P(JI)       
               ELSE                                                     
                  D10P(JI)= E12*CPA(I)*CSB(J)                           
                  D11P(JI)= E12*CPA(I)*CPB(J)                           
                  IF(D11P(JI).NE.ZER) THEN                              
                     TMP = ONE/D11P(JI)                                 
                     D00P(JI)= D00P(JI)*TMP                             
                     D01P(JI)= D01P(JI)*TMP                             
                     D10P(JI)= D10P(JI)*TMP                             
                  ENDIF                                                 
               ENDIF                                                    
            ELSE                                                        
               D00P(JI)= E12*CSA(I)*CSB(J)                              
            ENDIF                                                       
  160       CONTINUE                                                    
  170 JI=JI+1                                                           
C                                                                       
C BEGIN Q LOOP                                                          
C                                                                       
      IKL= 0                                                            
      DO 190 K=1,NGC                                                    
         X03= EXC(K)                                                    
         DO 180 L=1,NGD                                                 
            X04= EXD(L)                                                 
            X34= X03+X04                                                
            X43= ONE/X34                                                
            Y03= X03*X43                                                
            Y04= ONE-Y03                                                
            Y34= Y03*X04                                                
            R34Y34= R34*Y34                                             
            IF(R34Y34.GT.CUX) GO TO 180                                 
            E34= X43* EXP(-R34Y34)                                      
            TST= E34*CMAXC(K)*CMAXD(L)                                  
            IF(TST.LE.EI2) GO TO 180                                    
            ISMLQ= 0                                                    
            IF(TST.LE.EI1) ISMLQ= 1                                     
C                                                                       
C CQX = COMPONENT OF CQ ALONG PENULTIMATE X-AXIS                        
C CQZ = COMPONENT OF CQ ALONG PENULTIMATE Z-AXIS                        
C                                                                       
            CQ = RCD*Y04                                                
            CQX= CQ*SING                                                
            CQZ= CQ*COSG                                                
C                                                                       
C FIND COORDINATES OF Q RELATIVE TO AXES AT A                           
C QPR IS PERPENDICULAR FROM Q TO AB                                     
C                                                                       
            AQX= ACX+CQX                                                
            AQX2=AQX*AQX                                                
            AQXY=AQX*ACY                                                
            AQZ= ACZ+CQZ                                                
            QPS= AQX2+ACY2                                              
C                                                                       
            SQ(1)= E34*CSC(K)*CDD(L)                                    
            SQ(2)= E34*CPC(K)*CDD(L)                                    
            SQ(3)= E34*CDC(K)*CDD(L)                                    
C                                                                       
C USE SPECIAL FAST ROUTINE FOR INNER LOOPS FOR 0000 ... 1111            
C                                                                       
CJMS  ZEROING OF THE FQDx (x=0..8) ARRAYS OF LABELLED COMMON /FQ08/     
CJMS  TAKES PLACE IN THE INTJxx (xx=07..21) ROUTINES                    
C                                                                       
CJMS  ZEROING OF THE  RDx (x=0..8) ARRAYS OF LABELLED COMMON /KI4 /     
CJMS  TOOK PLACE IN THE CALL (FOR IKL=0) OF THE INTKxx (xx=07..20)      
CJMS  ROUTINES                                                          
C                                                                       
            IKL= IKL+1                                                  
            CALL SPDGEN(JTYPE,IKL)                                      
C                                                                       
  180    CONTINUE                                                       
  190 CONTINUE                                                          
C                                                                       
      QX= RCD*SING                                                      
      QZ= RCD*COSG                                                      
      IF(JTYPE.EQ. 7) THEN                                              
         CALL MCDV07(GHONDO,QX,QZ)                                      
      ELSEIF(JTYPE.EQ. 8) THEN                                          
         CALL MCDV08(GHONDO,QX,QZ)                                      
      ELSEIF(JTYPE.EQ. 9) THEN                                          
         CALL MCDV09(GHONDO,QX,QZ)                                      
      ELSEIF(JTYPE.EQ.10) THEN                                          
         CALL MCDV10(GHONDO,QX,QZ)                                      
      ELSEIF(JTYPE.EQ.11) THEN                                          
         CALL MCDV11(GHONDO,QX,QZ)                                      
      ELSEIF(JTYPE.EQ.12) THEN                                          
         CALL MCDV12(GHONDO,QX,QZ)                                      
      ELSEIF(JTYPE.EQ.13) THEN                                          
         CALL MCDV13(GHONDO,QX,QZ)                                      
      ELSEIF(JTYPE.EQ.14) THEN                                          
         CALL MCDV14(GHONDO,QX,QZ)                                      
      ELSEIF(JTYPE.EQ.15) THEN                                          
         CALL MCDV15(GHONDO,QX,QZ)                                      
      ELSEIF(JTYPE.EQ.16) THEN                                          
         CALL MCDV16(GHONDO,QX,QZ)                                      
      ELSEIF(JTYPE.EQ.17) THEN                                          
         CALL MCDV17(GHONDO,QX,QZ)                                      
      ELSEIF(JTYPE.EQ.18) THEN                                          
         CALL MCDV18(GHONDO,QX,QZ)                                      
      ELSEIF(JTYPE.EQ.19) THEN                                          
         CALL MCDV19(GHONDO,QX,QZ)                                      
      ELSEIF(JTYPE.EQ.20) THEN                                          
         CALL MCDV20(GHONDO,QX,QZ)                                      
      ELSEIF(JTYPE.EQ.21) THEN                                          
         CALL MCDV21(GHONDO,QX,QZ)                                      
      ENDIF                                                             
C                                                                       
CJMS  NOW, THE TRANSPOSE OF P TO BE USED FOR COMPUTATIONAL EFFICIENCY   
C                                                                       
      DO 195 J=1,2                                                      
         DO 195 I=J+1,3                                                 
            TMP= P(I,J)                                                 
            P(I,J)= P(J,I)                                              
            P(J,I)= TMP                                                 
  195 CONTINUE                                                          
C                                                                       
      CALL R30S1D(JTYPE,GHONDO,P)                                       
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2R   *DECK SPDGEN                                           
C>                                                                      
C>    @brief   s,p,d rotated axis type selection                        
C>                                                                      
C>    @details s,p,d rotated axis type selection                        
C>                                                                      
      SUBROUTINE SPDGEN(JTYPE,IKL)                                      
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      IF(JTYPE.EQ. 7) THEN                                              
         CALL INTJ07                                                    
         CALL INTK07(IKL)                                               
      ELSEIF(JTYPE.EQ. 8) THEN                                          
         CALL INTJ08                                                    
         CALL INTK08(IKL)                                               
      ELSEIF(JTYPE.EQ. 9) THEN                                          
         CALL INTJ09                                                    
         CALL INTK09(IKL)                                               
      ELSEIF(JTYPE.EQ.10) THEN                                          
         CALL INTJ10                                                    
         CALL INTK10(IKL)                                               
      ELSEIF(JTYPE.EQ.11) THEN                                          
         CALL INTJ11                                                    
         CALL INTK11(IKL)                                               
      ELSEIF(JTYPE.EQ.12) THEN                                          
         CALL INTJ12                                                    
         CALL INTK12(IKL)                                               
      ELSEIF(JTYPE.EQ.13) THEN                                          
         CALL INTJ13                                                    
         CALL INTK13(IKL)                                               
      ELSEIF(JTYPE.EQ.14) THEN                                          
         CALL INTJ14                                                    
         CALL INTK14(IKL)                                               
      ELSEIF(JTYPE.EQ.15) THEN                                          
         CALL INTJ15                                                    
         CALL INTK15(IKL)                                               
      ELSEIF(JTYPE.EQ.16) THEN                                          
         CALL INTJ16                                                    
         CALL INTK16(IKL)                                               
      ELSEIF(JTYPE.EQ.17) THEN                                          
         CALL INTJ17                                                    
         CALL INTK17(IKL)                                               
      ELSEIF(JTYPE.EQ.18) THEN                                          
         CALL INTJ18                                                    
         CALL INTK18(IKL)                                               
      ELSEIF(JTYPE.EQ.19) THEN                                          
         CALL INTJ19                                                    
         CALL INTK19(IKL)                                               
      ELSEIF(JTYPE.EQ.20) THEN                                          
         CALL INTJ20                                                    
         CALL INTK20(IKL)                                               
      ELSEIF(JTYPE.EQ.21) THEN                                          
         CALL INTJ21                                                    
         CALL INTK21(IKL)                                               
      ENDIF                                                             
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2R   *DECK R30S1D                                           
C>                                                                      
C>    @brief   rotate up to 1296 s,p,d integrals to space fixed axes    
C>                                                                      
C>    @details rotate up to 1296 s,p,d integrals to space fixed axes    
C>             INCOMING AND OUTGOING INTEGRALS IN F, while P(1,1),...   
C>             ARE DIRECTION COSINES OF SPACE FIXED AXES WRT AXES AT P  
C>                                                                      
      SUBROUTINE R30S1D(JTYPE,F,P)                                      
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      DIMENSION F(6,6,6,6),P(3,3)                                       
C                                                                       
      DIMENSION T(6),Q(6,6)                                             
C                                                                       
      PARAMETER (TWO=2.0D+00)                                           
      PARAMETER (SQRT3=1.732050807568877D+00)                           
C                                                                       
      DO 110 I=1,3                                                      
         Q(1,I)= P(1,I)*P(1,I)                                          
         Q(2,I)= P(2,I)*P(2,I)                                          
         Q(3,I)= P(3,I)*P(3,I)                                          
         Q(4,I)= P(1,I)*P(2,I)*TWO                                      
         Q(5,I)= P(1,I)*P(3,I)*TWO                                      
         Q(6,I)= P(2,I)*P(3,I)*TWO                                      
  110 CONTINUE                                                          
      DO 120 I=4,5                                                      
         J=I-2                                                          
         Q(1,I)= P(1,1)*P(1,J)                                          
         Q(2,I)= P(2,1)*P(2,J)                                          
         Q(3,I)= P(3,1)*P(3,J)                                          
         Q(4,I)= P(1,1)*P(2,J)+P(2,1)*P(1,J)                            
         Q(5,I)= P(1,1)*P(3,J)+P(3,1)*P(1,J)                            
         Q(6,I)= P(2,1)*P(3,J)+P(3,1)*P(2,J)                            
  120 CONTINUE                                                          
         Q(1,6)= P(1,2)*P(1,3)                                          
         Q(2,6)= P(2,2)*P(2,3)                                          
         Q(3,6)= P(3,2)*P(3,3)                                          
         Q(4,6)= P(1,2)*P(2,3)+P(2,2)*P(1,3)                            
         Q(5,6)= P(1,2)*P(3,3)+P(3,2)*P(1,3)                            
         Q(6,6)= P(2,2)*P(3,3)+P(3,2)*P(2,3)                            
      DO 130 I=4,6                                                      
         DO 130 J=1,6                                                   
            Q(J,I)= Q(J,I)*SQRT3                                        
  130 CONTINUE                                                          
C                                                                       
      IF(JTYPE.EQ. 7) THEN                                              
         DO 710 I=1,6                                                   
            T(I)= F(I,1,1,1)                                            
  710    CONTINUE                                                       
         DO 720 I=1,6                                                   
            F(I,1,1,1)= T(1)*Q(1,I)+T(2)*Q(2,I)+T(3)*Q(3,I)             
     *                 +T(4)*Q(4,I)+T(5)*Q(5,I)+T(6)*Q(6,I)             
  720    CONTINUE                                                       
C                                                                       
      ELSEIF(JTYPE.EQ. 8) THEN                                          
         DO 810 I=1,6                                                   
            T(1)= F(I,2,1,1)                                            
            T(2)= F(I,3,1,1)                                            
            T(3)= F(I,4,1,1)                                            
            F(I,2,1,1)= T(1)*P(1,1)+T(2)*P(2,1)+T(3)*P(3,1)             
            F(I,3,1,1)= T(1)*P(1,2)+T(2)*P(2,2)+T(3)*P(3,2)             
            F(I,4,1,1)= T(1)*P(1,3)+T(2)*P(2,3)+T(3)*P(3,3)             
  810    CONTINUE                                                       
         DO 830 J=1,4                                                   
            DO 820 I=1,6                                                
               T(I)= F(I,J,1,1)                                         
  820       CONTINUE                                                    
            DO 830 I=1,6                                                
               F(I,J,1,1)= T(1)*Q(1,I)+T(2)*Q(2,I)+T(3)*Q(3,I)          
     *                    +T(4)*Q(4,I)+T(5)*Q(5,I)+T(6)*Q(6,I)          
  830    CONTINUE                                                       
C                                                                       
      ELSEIF(JTYPE.EQ. 9) THEN                                          
         DO 910 I=1,6                                                   
            T(1)= F(I,1,2,1)                                            
            T(2)= F(I,1,3,1)                                            
            T(3)= F(I,1,4,1)                                            
            F(I,1,2,1)= T(1)*P(1,1)+T(2)*P(2,1)+T(3)*P(3,1)             
            F(I,1,3,1)= T(1)*P(1,2)+T(2)*P(2,2)+T(3)*P(3,2)             
            F(I,1,4,1)= T(1)*P(1,3)+T(2)*P(2,3)+T(3)*P(3,3)             
  910    CONTINUE                                                       
         DO 930 J=1,4                                                   
            DO 920 I=1,6                                                
               T(I)= F(I,1,J,1)                                         
  920       CONTINUE                                                    
            DO 930 I=1,6                                                
               F(I,1,J,1)= T(1)*Q(1,I)+T(2)*Q(2,I)+T(3)*Q(3,I)          
     *                    +T(4)*Q(4,I)+T(5)*Q(5,I)+T(6)*Q(6,I)          
  930    CONTINUE                                                       
C                                                                       
      ELSEIF(JTYPE.EQ.10) THEN                                          
         DO 1020 I=1,6                                                  
            DO 1010 J=1,6                                               
               T(J)= F(I,J,1,1)                                         
 1010       CONTINUE                                                    
            DO 1020 J=1,6                                               
               F(I,J,1,1)= T(1)*Q(1,J)+T(2)*Q(2,J)+T(3)*Q(3,J)          
     *                    +T(4)*Q(4,J)+T(5)*Q(5,J)+T(6)*Q(6,J)          
 1020    CONTINUE                                                       
         DO 1040 J=1,6                                                  
            DO 1030 I=1,6                                               
               T(I)= F(I,J,1,1)                                         
 1030       CONTINUE                                                    
            DO 1040 I=1,6                                               
               F(I,J,1,1)= T(1)*Q(1,I)+T(2)*Q(2,I)+T(3)*Q(3,I)          
     *                    +T(4)*Q(4,I)+T(5)*Q(5,I)+T(6)*Q(6,I)          
 1040    CONTINUE                                                       
C                                                                       
      ELSEIF(JTYPE.EQ.11) THEN                                          
         DO 1110 J=1,4                                                  
            DO 1110 I=1,6                                               
               T(1)= F(I,J,2,1)                                         
               T(2)= F(I,J,3,1)                                         
               T(3)= F(I,J,4,1)                                         
               F(I,J,2,1)= T(1)*P(1,1)+T(2)*P(2,1)+T(3)*P(3,1)          
               F(I,J,3,1)= T(1)*P(1,2)+T(2)*P(2,2)+T(3)*P(3,2)          
               F(I,J,4,1)= T(1)*P(1,3)+T(2)*P(2,3)+T(3)*P(3,3)          
 1110    CONTINUE                                                       
         DO 1120 K=1,4                                                  
            DO 1120 I=1,6                                               
               T(1)= F(I,2,K,1)                                         
               T(2)= F(I,3,K,1)                                         
               T(3)= F(I,4,K,1)                                         
               F(I,2,K,1)= T(1)*P(1,1)+T(2)*P(2,1)+T(3)*P(3,1)          
               F(I,3,K,1)= T(1)*P(1,2)+T(2)*P(2,2)+T(3)*P(3,2)          
               F(I,4,K,1)= T(1)*P(1,3)+T(2)*P(2,3)+T(3)*P(3,3)          
 1120    CONTINUE                                                       
         DO 1140 K=1,4                                                  
            DO 1140 J=1,4                                               
               DO 1130 I=1,6                                            
                  T(I)= F(I,J,K,1)                                      
 1130          CONTINUE                                                 
               DO 1140 I=1,6                                            
                  F(I,J,K,1)= T(1)*Q(1,I)+T(2)*Q(2,I)+T(3)*Q(3,I)       
     *                       +T(4)*Q(4,I)+T(5)*Q(5,I)+T(6)*Q(6,I)       
 1140    CONTINUE                                                       
C                                                                       
      ELSEIF(JTYPE.EQ.12) THEN                                          
         DO 1220 I=1,6                                                  
            DO 1210 K=1,6                                               
               T(K)= F(I,1,K,1)                                         
 1210       CONTINUE                                                    
            DO 1220 K=1,6                                               
               F(I,1,K,1)= T(1)*Q(1,K)+T(2)*Q(2,K)+T(3)*Q(3,K)          
     *                    +T(4)*Q(4,K)+T(5)*Q(5,K)+T(6)*Q(6,K)          
 1220    CONTINUE                                                       
         DO 1240 K=1,6                                                  
            DO 1230 I=1,6                                               
               T(I)= F(I,1,K,1)                                         
 1230       CONTINUE                                                    
            DO 1240 I=1,6                                               
               F(I,1,K,1)= T(1)*Q(1,I)+T(2)*Q(2,I)+T(3)*Q(3,I)          
     *                    +T(4)*Q(4,I)+T(5)*Q(5,I)+T(6)*Q(6,I)          
 1240    CONTINUE                                                       
C                                                                       
      ELSEIF(JTYPE.EQ.13) THEN                                          
         DO 1310 K=1,4                                                  
            DO 1310 I=1,6                                               
               T(1)= F(I,1,K,2)                                         
               T(2)= F(I,1,K,3)                                         
               T(3)= F(I,1,K,4)                                         
               F(I,1,K,2)= T(1)*P(1,1)+T(2)*P(2,1)+T(3)*P(3,1)          
               F(I,1,K,3)= T(1)*P(1,2)+T(2)*P(2,2)+T(3)*P(3,2)          
               F(I,1,K,4)= T(1)*P(1,3)+T(2)*P(2,3)+T(3)*P(3,3)          
 1310    CONTINUE                                                       
         DO 1320 L=1,4                                                  
            DO 1320 I=1,6                                               
               T(1)= F(I,1,2,L)                                         
               T(2)= F(I,1,3,L)                                         
               T(3)= F(I,1,4,L)                                         
               F(I,1,2,L)= T(1)*P(1,1)+T(2)*P(2,1)+T(3)*P(3,1)          
               F(I,1,3,L)= T(1)*P(1,2)+T(2)*P(2,2)+T(3)*P(3,2)          
               F(I,1,4,L)= T(1)*P(1,3)+T(2)*P(2,3)+T(3)*P(3,3)          
 1320    CONTINUE                                                       
         DO 1340 L=1,4                                                  
            DO 1340 K=1,4                                               
               DO 1330 I=1,6                                            
                  T(I)= F(I,1,K,L)                                      
 1330          CONTINUE                                                 
               DO 1340 I=1,6                                            
                  F(I,1,K,L)= T(1)*Q(1,I)+T(2)*Q(2,I)+T(3)*Q(3,I)       
     *                       +T(4)*Q(4,I)+T(5)*Q(5,I)+T(6)*Q(6,I)       
 1340    CONTINUE                                                       
C                                                                       
      ELSEIF(JTYPE.EQ.14) THEN                                          
         DO 1410 J=1,6                                                  
            DO 1410 I=1,6                                               
               T(1)= F(I,J,2,1)                                         
               T(2)= F(I,J,3,1)                                         
               T(3)= F(I,J,4,1)                                         
               F(I,J,2,1)= T(1)*P(1,1)+T(2)*P(2,1)+T(3)*P(3,1)          
               F(I,J,3,1)= T(1)*P(1,2)+T(2)*P(2,2)+T(3)*P(3,2)          
               F(I,J,4,1)= T(1)*P(1,3)+T(2)*P(2,3)+T(3)*P(3,3)          
 1410    CONTINUE                                                       
         DO 1430 K=1,4                                                  
            DO 1430 I=1,6                                               
               DO 1420 J=1,6                                            
                  T(J)= F(I,J,K,1)                                      
 1420          CONTINUE                                                 
               DO 1430 J=1,6                                            
                  F(I,J,K,1)= T(1)*Q(1,J)+T(2)*Q(2,J)+T(3)*Q(3,J)       
     *                       +T(4)*Q(4,J)+T(5)*Q(5,J)+T(6)*Q(6,J)       
 1430    CONTINUE                                                       
         DO 1450 K=1,4                                                  
            DO 1450 J=1,6                                               
               DO 1440 I=1,6                                            
                  T(I)= F(I,J,K,1)                                      
 1440          CONTINUE                                                 
               DO 1450 I=1,6                                            
                  F(I,J,K,1)= T(1)*Q(1,I)+T(2)*Q(2,I)+T(3)*Q(3,I)       
     *                       +T(4)*Q(4,I)+T(5)*Q(5,I)+T(6)*Q(6,I)       
 1450    CONTINUE                                                       
C                                                                       
      ELSEIF(JTYPE.EQ.15) THEN                                          
         DO 1520 J=1,4                                                  
            DO 1520 I=1,6                                               
               DO 1510 K=1,6                                            
                  T(K)= F(I,J,K,1)                                      
 1510          CONTINUE                                                 
               DO 1520 K=1,6                                            
                  F(I,J,K,1)= T(1)*Q(1,K)+T(2)*Q(2,K)+T(3)*Q(3,K)       
     *                       +T(4)*Q(4,K)+T(5)*Q(5,K)+T(6)*Q(6,K)       
 1520    CONTINUE                                                       
         DO 1530 K=1,6                                                  
            DO 1530 I=1,6                                               
               T(1)= F(I,2,K,1)                                         
               T(2)= F(I,3,K,1)                                         
               T(3)= F(I,4,K,1)                                         
               F(I,2,K,1)= T(1)*P(1,1)+T(2)*P(2,1)+T(3)*P(3,1)          
               F(I,3,K,1)= T(1)*P(1,2)+T(2)*P(2,2)+T(3)*P(3,2)          
               F(I,4,K,1)= T(1)*P(1,3)+T(2)*P(2,3)+T(3)*P(3,3)          
 1530    CONTINUE                                                       
         DO 1550 K=1,6                                                  
            DO 1550 J=1,4                                               
               DO 1540 I=1,6                                            
                  T(I)= F(I,J,K,1)                                      
 1540          CONTINUE                                                 
               DO 1550 I=1,6                                            
                  F(I,J,K,1)= T(1)*Q(1,I)+T(2)*Q(2,I)+T(3)*Q(3,I)       
     *                       +T(4)*Q(4,I)+T(5)*Q(5,I)+T(6)*Q(6,I)       
 1550    CONTINUE                                                       
C                                                                       
      ELSEIF(JTYPE.EQ.16) THEN                                          
         DO 1610 K=1,4                                                  
            DO 1610 J=1,4                                               
               DO 1610 I=1,6                                            
                  T(1)= F(I,J,K,2)                                      
                  T(2)= F(I,J,K,3)                                      
                  T(3)= F(I,J,K,4)                                      
                  F(I,J,K,2)= T(1)*P(1,1)+T(2)*P(2,1)+T(3)*P(3,1)       
                  F(I,J,K,3)= T(1)*P(1,2)+T(2)*P(2,2)+T(3)*P(3,2)       
                  F(I,J,K,4)= T(1)*P(1,3)+T(2)*P(2,3)+T(3)*P(3,3)       
 1610    CONTINUE                                                       
         DO 1620 L=1,4                                                  
            DO 1620 J=1,4                                               
               DO 1620 I=1,6                                            
                  T(1)= F(I,J,2,L)                                      
                  T(2)= F(I,J,3,L)                                      
                  T(3)= F(I,J,4,L)                                      
                  F(I,J,2,L)= T(1)*P(1,1)+T(2)*P(2,1)+T(3)*P(3,1)       
                  F(I,J,3,L)= T(1)*P(1,2)+T(2)*P(2,2)+T(3)*P(3,2)       
                  F(I,J,4,L)= T(1)*P(1,3)+T(2)*P(2,3)+T(3)*P(3,3)       
 1620    CONTINUE                                                       
         DO 1630 L=1,4                                                  
            DO 1630 K=1,4                                               
               DO 1630 I=1,6                                            
                  T(1)= F(I,2,K,L)                                      
                  T(2)= F(I,3,K,L)                                      
                  T(3)= F(I,4,K,L)                                      
                  F(I,2,K,L)= T(1)*P(1,1)+T(2)*P(2,1)+T(3)*P(3,1)       
                  F(I,3,K,L)= T(1)*P(1,2)+T(2)*P(2,2)+T(3)*P(3,2)       
                  F(I,4,K,L)= T(1)*P(1,3)+T(2)*P(2,3)+T(3)*P(3,3)       
 1630    CONTINUE                                                       
         DO 1650 L=1,4                                                  
            DO 1650 K=1,4                                               
               DO 1650 J=1,4                                            
                  DO 1640 I=1,6                                         
                     T(I)= F(I,J,K,L)                                   
 1640             CONTINUE                                              
                  DO 1650 I=1,6                                         
                     F(I,J,K,L)= T(1)*Q(1,I)+T(2)*Q(2,I)+T(3)*Q(3,I)    
     *                          +T(4)*Q(4,I)+T(5)*Q(5,I)+T(6)*Q(6,I)    
 1650    CONTINUE                                                       
C                                                                       
      ELSEIF(JTYPE.EQ.17) THEN                                          
         DO 1720 J=1,6                                                  
            DO 1720 I=1,6                                               
               DO 1710 K=1,6                                            
                  T(K)= F(I,J,K,1)                                      
 1710          CONTINUE                                                 
               DO 1720 K=1,6                                            
                  F(I,J,K,1)= T(1)*Q(1,K)+T(2)*Q(2,K)+T(3)*Q(3,K)       
     *                       +T(4)*Q(4,K)+T(5)*Q(5,K)+T(6)*Q(6,K)       
 1720    CONTINUE                                                       
         DO 1740 K=1,6                                                  
            DO 1740 I=1,6                                               
               DO 1730 J=1,6                                            
                  T(J)= F(I,J,K,1)                                      
 1730          CONTINUE                                                 
               DO 1740 J=1,6                                            
                  F(I,J,K,1)= T(1)*Q(1,J)+T(2)*Q(2,J)+T(3)*Q(3,J)       
     *                       +T(4)*Q(4,J)+T(5)*Q(5,J)+T(6)*Q(6,J)       
 1740    CONTINUE                                                       
         DO 1760 K=1,6                                                  
            DO 1760 J=1,6                                               
               DO 1750 I=1,6                                            
                  T(I)= F(I,J,K,1)                                      
 1750          CONTINUE                                                 
               DO 1760 I=1,6                                            
                  F(I,J,K,1)= T(1)*Q(1,I)+T(2)*Q(2,I)+T(3)*Q(3,I)       
     *                       +T(4)*Q(4,I)+T(5)*Q(5,I)+T(6)*Q(6,I)       
 1760    CONTINUE                                                       
C                                                                       
      ELSEIF(JTYPE.EQ.18) THEN                                          
         DO 1810 K=1,4                                                  
            DO 1810 J=1,6                                               
               DO 1810 I=1,6                                            
                  T(1)= F(I,J,K,2)                                      
                  T(2)= F(I,J,K,3)                                      
                  T(3)= F(I,J,K,4)                                      
                  F(I,J,K,2)= T(1)*P(1,1)+T(2)*P(2,1)+T(3)*P(3,1)       
                  F(I,J,K,3)= T(1)*P(1,2)+T(2)*P(2,2)+T(3)*P(3,2)       
                  F(I,J,K,4)= T(1)*P(1,3)+T(2)*P(2,3)+T(3)*P(3,3)       
 1810    CONTINUE                                                       
         DO 1820 L=1,4                                                  
            DO 1820 J=1,6                                               
               DO 1820 I=1,6                                            
                  T(1)= F(I,J,2,L)                                      
                  T(2)= F(I,J,3,L)                                      
                  T(3)= F(I,J,4,L)                                      
                  F(I,J,2,L)= T(1)*P(1,1)+T(2)*P(2,1)+T(3)*P(3,1)       
                  F(I,J,3,L)= T(1)*P(1,2)+T(2)*P(2,2)+T(3)*P(3,2)       
                  F(I,J,4,L)= T(1)*P(1,3)+T(2)*P(2,3)+T(3)*P(3,3)       
 1820    CONTINUE                                                       
         DO 1840 L=1,4                                                  
            DO 1840 K=1,4                                               
               DO 1840 I=1,6                                            
                  DO 1830 J=1,6                                         
                     T(J)= F(I,J,K,L)                                   
 1830             CONTINUE                                              
                  DO 1840 J=1,6                                         
                     F(I,J,K,L)= T(1)*Q(1,J)+T(2)*Q(2,J)+T(3)*Q(3,J)    
     *                          +T(4)*Q(4,J)+T(5)*Q(5,J)+T(6)*Q(6,J)    
 1840    CONTINUE                                                       
         DO 1860 L=1,4                                                  
            DO 1860 K=1,4                                               
               DO 1860 J=1,6                                            
                  DO 1850 I=1,6                                         
                     T(I)= F(I,J,K,L)                                   
 1850             CONTINUE                                              
                  DO 1860 I=1,6                                         
                     F(I,J,K,L)= T(1)*Q(1,I)+T(2)*Q(2,I)+T(3)*Q(3,I)    
     *                          +T(4)*Q(4,I)+T(5)*Q(5,I)+T(6)*Q(6,I)    
 1860    CONTINUE                                                       
      ELSEIF(JTYPE.EQ.19) THEN                                          
         DO 1910 K=1,6                                                  
            DO 1910 J=1,4                                               
               DO 1910 I=1,6                                            
                  T(1)= F(I,J,K,2)                                      
                  T(2)= F(I,J,K,3)                                      
                  T(3)= F(I,J,K,4)                                      
                  F(I,J,K,2)= T(1)*P(1,1)+T(2)*P(2,1)+T(3)*P(3,1)       
                  F(I,J,K,3)= T(1)*P(1,2)+T(2)*P(2,2)+T(3)*P(3,2)       
                  F(I,J,K,4)= T(1)*P(1,3)+T(2)*P(2,3)+T(3)*P(3,3)       
 1910    CONTINUE                                                       
         DO 1930 L=1,4                                                  
            DO 1930 J=1,4                                               
               DO 1930 I=1,6                                            
                  DO 1920 K=1,6                                         
                     T(K)= F(I,J,K,L)                                   
 1920             CONTINUE                                              
                  DO 1930 K=1,6                                         
                     F(I,J,K,L)= T(1)*Q(1,K)+T(2)*Q(2,K)+T(3)*Q(3,K)    
     *                          +T(4)*Q(4,K)+T(5)*Q(5,K)+T(6)*Q(6,K)    
 1930    CONTINUE                                                       
         DO 1940 L=1,4                                                  
            DO 1940 K=1,6                                               
               DO 1940 I=1,6                                            
                  T(1)= F(I,2,K,L)                                      
                  T(2)= F(I,3,K,L)                                      
                  T(3)= F(I,4,K,L)                                      
                  F(I,2,K,L)= T(1)*P(1,1)+T(2)*P(2,1)+T(3)*P(3,1)       
                  F(I,3,K,L)= T(1)*P(1,2)+T(2)*P(2,2)+T(3)*P(3,2)       
                  F(I,4,K,L)= T(1)*P(1,3)+T(2)*P(2,3)+T(3)*P(3,3)       
 1940    CONTINUE                                                       
         DO 1960 L=1,4                                                  
            DO 1960 K=1,6                                               
               DO 1960 J=1,4                                            
                  DO 1950 I=1,6                                         
                     T(I)= F(I,J,K,L)                                   
 1950             CONTINUE                                              
                  DO 1960 I=1,6                                         
                     F(I,J,K,L)= T(1)*Q(1,I)+T(2)*Q(2,I)+T(3)*Q(3,I)    
     *                          +T(4)*Q(4,I)+T(5)*Q(5,I)+T(6)*Q(6,I)    
 1960    CONTINUE                                                       
C                                                                       
      ELSEIF(JTYPE.EQ.20) THEN                                          
         DO 2010 K=1,6                                                  
            DO 2010 J=1,6                                               
               DO 2010 I=1,6                                            
                  T(1)= F(I,J,K,2)                                      
                  T(2)= F(I,J,K,3)                                      
                  T(3)= F(I,J,K,4)                                      
                  F(I,J,K,2)= T(1)*P(1,1)+T(2)*P(2,1)+T(3)*P(3,1)       
                  F(I,J,K,3)= T(1)*P(1,2)+T(2)*P(2,2)+T(3)*P(3,2)       
                  F(I,J,K,4)= T(1)*P(1,3)+T(2)*P(2,3)+T(3)*P(3,3)       
 2010    CONTINUE                                                       
         DO 2030 L=1,4                                                  
            DO 2030 J=1,6                                               
               DO 2030 I=1,6                                            
                  DO 2020 K=1,6                                         
                     T(K)= F(I,J,K,L)                                   
 2020             CONTINUE                                              
                  DO 2030 K=1,6                                         
                     F(I,J,K,L)= T(1)*Q(1,K)+T(2)*Q(2,K)+T(3)*Q(3,K)    
     *                          +T(4)*Q(4,K)+T(5)*Q(5,K)+T(6)*Q(6,K)    
 2030    CONTINUE                                                       
         DO 2050 L=1,4                                                  
            DO 2050 K=1,6                                               
               DO 2050 I=1,6                                            
                  DO 2040 J=1,6                                         
                     T(J)= F(I,J,K,L)                                   
 2040             CONTINUE                                              
                  DO 2050 J=1,6                                         
                     F(I,J,K,L)= T(1)*Q(1,J)+T(2)*Q(2,J)+T(3)*Q(3,J)    
     *                          +T(4)*Q(4,J)+T(5)*Q(5,J)+T(6)*Q(6,J)    
 2050    CONTINUE                                                       
         DO 2070 L=1,4                                                  
            DO 2070 K=1,6                                               
               DO 2070 J=1,6                                            
                  DO 2060 I=1,6                                         
                     T(I)= F(I,J,K,L)                                   
 2060             CONTINUE                                              
                  DO 2070 I=1,6                                         
                     F(I,J,K,L)= T(1)*Q(1,I)+T(2)*Q(2,I)+T(3)*Q(3,I)    
     *                          +T(4)*Q(4,I)+T(5)*Q(5,I)+T(6)*Q(6,I)    
 2070    CONTINUE                                                       
C                                                                       
      ELSEIF(JTYPE.EQ.21) THEN                                          
         DO 2120 K=1,6                                                  
            DO 2120 J=1,6                                               
               DO 2120 I=1,6                                            
                  DO 2110 L=1,6                                         
                     T(L)= F(I,J,K,L)                                   
 2110             CONTINUE                                              
                  DO 2120 L=1,6                                         
                     F(I,J,K,L)= T(1)*Q(1,L)+T(2)*Q(2,L)+T(3)*Q(3,L)    
     *                          +T(4)*Q(4,L)+T(5)*Q(5,L)+T(6)*Q(6,L)    
 2120    CONTINUE                                                       
         DO 2140 L=1,6                                                  
            DO 2140 J=1,6                                               
               DO 2140 I=1,6                                            
                  DO 2130 K=1,6                                         
                     T(K)= F(I,J,K,L)                                   
 2130             CONTINUE                                              
                  DO 2140 K=1,6                                         
                     F(I,J,K,L)= T(1)*Q(1,K)+T(2)*Q(2,K)+T(3)*Q(3,K)    
     *                          +T(4)*Q(4,K)+T(5)*Q(5,K)+T(6)*Q(6,K)    
 2140    CONTINUE                                                       
         DO 2160 L=1,6                                                  
            DO 2160 K=1,6                                               
               DO 2160 I=1,6                                            
                  DO 2150 J=1,6                                         
                     T(J)= F(I,J,K,L)                                   
 2150             CONTINUE                                              
                  DO 2160 J=1,6                                         
                     F(I,J,K,L)= T(1)*Q(1,J)+T(2)*Q(2,J)+T(3)*Q(3,J)    
     *                          +T(4)*Q(4,J)+T(5)*Q(5,J)+T(6)*Q(6,J)    
 2160    CONTINUE                                                       
         DO 2180 L=1,6                                                  
            DO 2180 K=1,6                                               
               DO 2180 J=1,6                                            
                  DO 2170 I=1,6                                         
                     T(I)= F(I,J,K,L)                                   
 2170             CONTINUE                                              
                  DO 2180 I=1,6                                         
                     F(I,J,K,L)= T(1)*Q(1,I)+T(2)*Q(2,I)+T(3)*Q(3,I)    
     *                          +T(4)*Q(4,I)+T(5)*Q(5,I)+T(6)*Q(6,I)    
 2180    CONTINUE                                                       
      ENDIF                                                             
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2R   *DECK INTJ07                                           
C>                                                                      
C>    @brief   DSSS case                                                
C>                                                                      
C>    @details integration of the DSSS case                             
C>                                                                      
      SUBROUTINE INTJ07                                                 
      USE lrcdft, ONLY: LCFLAG, EMU, EMU2, LRFILE                       
      use mx_limits, only: mxgtot,mxgsh,mxg2                            
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
C GENERATE JTYPE= 7 INTEGRALS                                           
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
      COMMON /FQ08  / FQD(0:8),FQD0(5),FQD1(26),FQD2(48),FQD3(64),      
     2                FQD4(80),FQD5(66),FQD6(49),FQD7(24),FQD8(9)       
C$omp threadprivate(/FQ08/)
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS,                      
     2                TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2),      
     3                D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2),      
     4                NGANGB                                            
C$omp threadprivate(/GEOMPQ/)
      COMMON /MAXC  / CMAX(MXGTOT),CMAXA(MXGSH),CMAXB(MXGSH),           
     2                CMAXC(MXGSH),CMAXD(MXGSH),ISMLP(MXG2),ISMLQ       
C$omp threadprivate(/MAXC/)
      LOGICAL         LRINT                                             
      COMMON /NLRCF / LRINT                                             
C$omp threadprivate(/NLRCF /)
C                                                                       
      PARAMETER (ZER=0.0D+00)                                           
      PARAMETER (PT5=0.5D+00)                                           
      PARAMETER (ONE=1.0D+00)                                           
      PARAMETER (PI4=0.7853981633974483D+00)                            
C                                                                       
         FQD0(1)= ZER                                                   
         FQD1(1)= ZER                                                   
         FQD1(2)= ZER                                                   
         FQD2(1)= ZER                                                   
         FQD2(2)= ZER                                                   
         FQD2(3)= ZER                                                   
C                                                                       
      DO 300 I=1,NGANGB                                                 
         ISML= ISMLQ+ISMLP(I)                                           
         IF(ISML.GE.2) GO TO 300                                        
         X12= TX12(I)                                                   
         Y02= TY02(I)                                                   
         FQZ= D00P(I)                                                   
         X41= ONE/(X12+X34)                                             
         PQR= Y02-AQZ                                                   
         PQS= PQR*PQR                                                   
         RHO= X12*X34*X41                                               
         IF(LRINT) THEN                                                 
            EFR= EMU2/(EMU2+RHO)                                        
            RHO= RHO*EFR                                                
            FQZ= FQZ*SQRT(EFR)                                          
         ENDIF                                                          
         XVA=(PQS+QPS)*RHO                                              
         RHO= RHO+RHO                                                   
         N=2                                                            
         IF(XVA.LE.TMAX) THEN                                           
C                                                                       
C     FM(T) EVALUATION                                                  
C     SEE NOTES IN ROUTINE INTJ14 ABOUT THE FM(T) INTERPOLATIONS        
C                                                                       
            TV= XVA*RFINC(N)                                            
            IP= NINT(TV)                                                
            FX=    FGRID(4,IP,N) *TV                                    
            FX=(FX+FGRID(3,IP,N))*TV                                    
            FX=(FX+FGRID(2,IP,N))*TV                                    
            FX=(FX+FGRID(1,IP,N))*TV                                    
            FX= FX+FGRID(0,IP,N)                                        
            TV= XVA*RXINC                                               
            IP= NINT(TV)                                                
            ET=    XGRID(4,IP) *TV                                      
            ET=(ET+XGRID(3,IP))*TV                                      
            ET=(ET+XGRID(2,IP))*TV                                      
            ET=(ET+XGRID(1,IP))*TV                                      
            ET= ET+XGRID(0,IP)                                          
C                                                                       
            FQD(N)= FX                                                  
            T2= XVA+XVA                                                 
               FQD(2-1)=(T2*FQD(2)+ET)*RMR(2)                           
               FQD(1-1)=(T2*FQD(1)+ET)*RMR(1)                           
CC          DO M=N,1,-1                                                 
CC             FQD(M-1)=(T2*FQD(M)+ET)*RMR(M)                           
CC          END DO                                                      
C                                                                       
            FQF= FQZ*SQRT(X41)                                          
               FQD(0)= FQD(0)*FQF                                       
            FQF= FQF*RHO                                                
               FQD(1)= FQD(1)*FQF                                       
            FQF= FQF*RHO                                                
               FQD(2)= FQD(2)*FQF                                       
CC          DO 210 M=0,N                                                
CC             FQD(M)= FQD(M)*FQF                                       
CC210       FQF= FQF*RHO                                                
         ELSE                                                           
            XIN= ONE/XVA                                                
            FQD(0)= FQZ*SQRT(PI4*XIN*X41)                               
            ROX= RHO*XIN                                                
            FQF= PT5*ROX                                                
               FQD(1)= FQD(0)*FQF                                       
            FQF= FQF+ROX                                                
               FQD(2)= FQD(1)*FQF                                       
CC          DO 220 M=1,N                                                
CC             FQD(M)= FQD(M-1)*FQF                                     
CC220       FQF= FQF+ROX                                                
         ENDIF                                                          
C                                                                       
         FQD0(1)= FQD0(1)+FQD(0)                                        
         FQD1(1)= FQD1(1)+FQD(1)                                        
         FQD1(2)= FQD1(2)+FQD(1)*PQR                                    
         FQD2(1)= FQD2(1)+FQD(2)                                        
         FQD2(2)= FQD2(2)+FQD(2)*PQR                                    
         FQD2(3)= FQD2(3)+FQD(2)*PQS                                    
  300 CONTINUE                                                          
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2R   *DECK INTJ08                                           
C>                                                                      
C>    @brief   DPSS case                                                
C>                                                                      
C>    @details integration of the DPSS case                             
C>                                                                      
      SUBROUTINE INTJ08                                                 
      USE lrcdft, ONLY: LCFLAG, EMU, EMU2, LRFILE                       
      use mx_limits, only: mxgtot,mxgsh,mxg2                            
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
C GENERATE JTYPE= 8 INTEGRALS                                           
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
      COMMON /FQ08  / FQD(0:8),FQD0(5),FQD1(26),FQD2(48),FQD3(64),      
     2                FQD4(80),FQD5(66),FQD6(49),FQD7(24),FQD8(9)       
C$omp threadprivate(/FQ08/)
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS,                      
     2                TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2),      
     3                D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2),      
     4                NGANGB                                            
C$omp threadprivate(/GEOMPQ/)
      COMMON /MAXC  / CMAX(MXGTOT),CMAXA(MXGSH),CMAXB(MXGSH),           
     2                CMAXC(MXGSH),CMAXD(MXGSH),ISMLP(MXG2),ISMLQ       
C$omp threadprivate(/MAXC/)
      LOGICAL         LRINT                                             
      COMMON /NLRCF / LRINT                                             
C$omp threadprivate(/NLRCF /)
C                                                                       
      PARAMETER (ZER=0.0D+00)                                           
      PARAMETER (PT5=0.5D+00)                                           
      PARAMETER (ONE=1.0D+00)                                           
      PARAMETER (PI4=0.7853981633974483D+00)                            
C                                                                       
         FQD0(1)= ZER                                                   
         FQD1(1)= ZER                                                   
         FQD1(2)= ZER                                                   
         FQD2(1)= ZER                                                   
         FQD2(2)= ZER                                                   
         FQD2(3)= ZER                                                   
         FQD3(1)= ZER                                                   
         FQD3(2)= ZER                                                   
         FQD3(3)= ZER                                                   
         FQD3(4)= ZER                                                   
C                                                                       
      DO 300 I=1,NGANGB                                                 
         ISML= ISMLQ+ISMLP(I)                                           
         IF(ISML.GE.2) GO TO 300                                        
         X12= TX12(I)                                                   
         Y02= TY02(I)                                                   
         FQZ= D00P(I)                                                   
         X41= ONE/(X12+X34)                                             
         PQR= Y02-AQZ                                                   
         PQS= PQR*PQR                                                   
         RHO= X12*X34*X41                                               
         IF(LRINT) THEN                                                 
            EFR= EMU2/(EMU2+RHO)                                        
            RHO= RHO*EFR                                                
            FQZ= FQZ*SQRT(EFR)                                          
         ENDIF                                                          
         XVA=(PQS+QPS)*RHO                                              
         RHO= RHO+RHO                                                   
         N=3                                                            
         IF(XVA.LE.TMAX) THEN                                           
C                                                                       
C     FM(T) EVALUATION...DOWNWARD RECURSION FOR M=3                     
C                                                                       
            TV= XVA*RFINC(N)                                            
            IP= NINT(TV)                                                
            FX=    FGRID(4,IP,N) *TV                                    
            FX=(FX+FGRID(3,IP,N))*TV                                    
            FX=(FX+FGRID(2,IP,N))*TV                                    
            FX=(FX+FGRID(1,IP,N))*TV                                    
            FX= FX+FGRID(0,IP,N)                                        
            TV= XVA*RXINC                                               
            IP= NINT(TV)                                                
            ET=    XGRID(4,IP) *TV                                      
            ET=(ET+XGRID(3,IP))*TV                                      
            ET=(ET+XGRID(2,IP))*TV                                      
            ET=(ET+XGRID(1,IP))*TV                                      
            ET= ET+XGRID(0,IP)                                          
C                                                                       
            FQD(N)= FX                                                  
            T2= XVA+XVA                                                 
               FQD(3-1)=(T2*FQD(3)+ET)*RMR(3)                           
               FQD(2-1)=(T2*FQD(2)+ET)*RMR(2)                           
               FQD(1-1)=(T2*FQD(1)+ET)*RMR(1)                           
CC          DO M=N,1,-1                                                 
CC             FQD(M-1)=(T2*FQD(M)+ET)*RMR(M)                           
CC          END DO                                                      
C                                                                       
            FQF= FQZ*SQRT(X41)                                          
               FQD(0)= FQD(0)*FQF                                       
            FQF= FQF*RHO                                                
               FQD(1)= FQD(1)*FQF                                       
            FQF= FQF*RHO                                                
               FQD(2)= FQD(2)*FQF                                       
            FQF= FQF*RHO                                                
               FQD(3)= FQD(3)*FQF                                       
CC          DO 210 M=0,N                                                
CC             FQD(M)= FQD(M)*FQF                                       
CC210       FQF= FQF*RHO                                                
         ELSE                                                           
            XIN= ONE/XVA                                                
            FQD(0)= FQZ*SQRT(PI4*XIN*X41)                               
            ROX= RHO*XIN                                                
            FQF= PT5*ROX                                                
               FQD(1)= FQD(0)*FQF                                       
            FQF= FQF+ROX                                                
               FQD(2)= FQD(1)*FQF                                       
            FQF= FQF+ROX                                                
               FQD(3)= FQD(2)*FQF                                       
CC          DO 220 M=1,N                                                
CC             FQD(M)= FQD(M-1)*FQF                                     
CC220       FQF= FQF+ROX                                                
         ENDIF                                                          
C                                                                       
         FQD0(1)= FQD0(1)+FQD(0)                                        
         FQD1(1)= FQD1(1)+FQD(1)                                        
         FQD1(2)= FQD1(2)+FQD(1)*PQR                                    
         FQD2(1)= FQD2(1)+FQD(2)                                        
         FQD2(2)= FQD2(2)+FQD(2)*PQR                                    
         FQD2(3)= FQD2(3)+FQD(2)*PQS                                    
         FQD3(1)= FQD3(1)+FQD(3)                                        
         FQD3(2)= FQD3(2)+FQD(3)*PQR                                    
         FQD3(3)= FQD3(3)+FQD(3)*PQS                                    
         FQD3(4)= FQD3(4)+FQD(3)*PQS*PQR                                
  300 CONTINUE                                                          
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2R   *DECK INTJ09                                           
C>                                                                      
C>    @brief   DSPS case                                                
C>                                                                      
C>    @details integration of the DSPS case                             
C>                                                                      
      SUBROUTINE INTJ09                                                 
      USE lrcdft, ONLY: LCFLAG, EMU, EMU2, LRFILE                       
      use mx_limits, only: mxgtot,mxgsh,mxg2                            
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
C GENERATE JTYPE= 9 INTEGRALS                                           
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
      COMMON /FQ08  / FQD(0:8),FQD0(5),FQD1(2,13),FQD2(3,16),FQD3(64),  
     2                FQD4(80),FQD5(66),FQD6(49),FQD7(24),FQD8(9)       
C$omp threadprivate(/FQ08/)
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS,                      
     2                TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2),      
     3                D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2),      
     4                NGANGB                                            
C$omp threadprivate(/GEOMPQ/)
      COMMON /MAXC  / CMAX(MXGTOT),CMAXA(MXGSH),CMAXB(MXGSH),           
     2                CMAXC(MXGSH),CMAXD(MXGSH),ISMLP(MXG2),ISMLQ       
C$omp threadprivate(/MAXC/)
      LOGICAL         LRINT                                             
      COMMON /NLRCF / LRINT                                             
C$omp threadprivate(/NLRCF /)
C                                                                       
      DIMENSION  WORK(3,3)                                              
C                                                                       
      PARAMETER (ZER=0.0D+00)                                           
      PARAMETER (PT5=0.5D+00)                                           
      PARAMETER (ONE=1.0D+00)                                           
      PARAMETER (PI4=0.7853981633974483D+00)                            
C                                                                       
         FQD0(1)= ZER                                                   
         FQD0(2)= ZER                                                   
         DO J= 1, 3                                                     
            FQD1(1,J)= ZER                                              
            FQD1(2,J)= ZER                                              
C                                                                       
            FQD2(1,J)= ZER                                              
            FQD2(2,J)= ZER                                              
            FQD2(3,J)= ZER                                              
         ENDDO                                                          
         FQD3(1)= ZER                                                   
         FQD3(2)= ZER                                                   
         FQD3(3)= ZER                                                   
         FQD3(4)= ZER                                                   
C                                                                       
      DO 300 I=1,NGANGB                                                 
         ISML= ISMLQ+ISMLP(I)                                           
         IF(ISML.GE.2) GO TO 300                                        
         X12= TX12(I)                                                   
         Y02= TY02(I)                                                   
         FQZ= D01P(I)                                                   
         X41= ONE/(X12+X34)                                             
         PQR= Y02-AQZ                                                   
         PQS= PQR*PQR                                                   
         RHO= X12*X34*X41                                               
         IF(LRINT) THEN                                                 
            EFR= EMU2/(EMU2+RHO)                                        
            RHO= RHO*EFR                                                
            FQZ= FQZ*SQRT(EFR)                                          
         ENDIF                                                          
         XVA=(PQS+QPS)*RHO                                              
         RHO= RHO+RHO                                                   
         N=3                                                            
         IF(XVA.LE.TMAX) THEN                                           
C                                                                       
C     FM(T) EVALUATION...DOWNWARD RECURSION FOR M=3                     
C                                                                       
            TV= XVA*RFINC(N)                                            
            IP= NINT(TV)                                                
            FX=    FGRID(4,IP,N) *TV                                    
            FX=(FX+FGRID(3,IP,N))*TV                                    
            FX=(FX+FGRID(2,IP,N))*TV                                    
            FX=(FX+FGRID(1,IP,N))*TV                                    
            FX= FX+FGRID(0,IP,N)                                        
            TV= XVA*RXINC                                               
            IP= NINT(TV)                                                
            ET=    XGRID(4,IP) *TV                                      
            ET=(ET+XGRID(3,IP))*TV                                      
            ET=(ET+XGRID(2,IP))*TV                                      
            ET=(ET+XGRID(1,IP))*TV                                      
            ET= ET+XGRID(0,IP)                                          
C                                                                       
            FQD(N)= FX                                                  
            T2= XVA+XVA                                                 
               FQD(3-1)=(T2*FQD(3)+ET)*RMR(3)                           
               FQD(2-1)=(T2*FQD(2)+ET)*RMR(2)                           
               FQD(1-1)=(T2*FQD(1)+ET)*RMR(1)                           
CC          DO M=N,1,-1                                                 
CC             FQD(M-1)=(T2*FQD(M)+ET)*RMR(M)                           
CC          END DO                                                      
C                                                                       
            FQF= FQZ*SQRT(X41)                                          
               FQD(0)= FQD(0)*FQF                                       
            FQF= FQF*RHO                                                
               FQD(1)= FQD(1)*FQF                                       
            FQF= FQF*RHO                                                
               FQD(2)= FQD(2)*FQF                                       
            FQF= FQF*RHO                                                
               FQD(3)= FQD(3)*FQF                                       
CC          DO 210 M=0,N                                                
CC             FQD(M)= FQD(M)*FQF                                       
CC210       FQF= FQF*RHO                                                
         ELSE                                                           
            XIN= ONE/XVA                                                
            FQD(0)= FQZ*SQRT(PI4*XIN*X41)                               
            ROX= RHO*XIN                                                
            FQF= PT5*ROX                                                
               FQD(1)= FQD(0)*FQF                                       
            FQF= FQF+ROX                                                
               FQD(2)= FQD(1)*FQF                                       
            FQF= FQF+ROX                                                
               FQD(3)= FQD(2)*FQF                                       
CC          DO 220 M=1,N                                                
CC             FQD(M)= FQD(M-1)*FQF                                     
CC220       FQF= FQF+ROX                                                
         ENDIF                                                          
C                                                                       
         WORK(1,1)= D00P(I)                                             
         WORK(1,2)= TY01(I)                                             
         WORK(1,3)= TX21(I)                                             
         DO J= 1, 3                                                     
            WORK(2,J)= WORK(1,J)*PQR                                    
            WORK(3,J)= WORK(1,J)*PQS                                    
         ENDDO                                                          
C                                                                       
         FQD0(1)= FQD0(1)+FQD(0)*WORK(1,1)                              
         FQD0(2)= FQD0(2)+FQD(0)*WORK(1,2)                              
         DO J= 1, 3                                                     
            FQD1(1,J)= FQD1(1,J)+FQD(1)*WORK(1,J)                       
            FQD1(2,J)= FQD1(2,J)+FQD(1)*WORK(2,J)                       
C                                                                       
            FQD2(1,J)= FQD2(1,J)+FQD(2)*WORK(1,J)                       
            FQD2(2,J)= FQD2(2,J)+FQD(2)*WORK(2,J)                       
            FQD2(3,J)= FQD2(3,J)+FQD(2)*WORK(3,J)                       
         ENDDO                                                          
         FQD3(1)= FQD3(1)+FQD(3)*WORK(1,3)                              
         FQD3(2)= FQD3(2)+FQD(3)*WORK(2,3)                              
         FQD3(3)= FQD3(3)+FQD(3)*WORK(3,3)                              
         FQD3(4)= FQD3(4)+FQD(3)*WORK(3,3)*PQR                          
  300 CONTINUE                                                          
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2R   *DECK INTJ10                                           
C>                                                                      
C>    @brief   DDSS case                                                
C>                                                                      
C>    @details integration of the DDSS case                             
C>                                                                      
      SUBROUTINE INTJ10                                                 
      USE lrcdft, ONLY: LCFLAG, EMU, EMU2, LRFILE                       
      use mx_limits, only: mxgtot,mxgsh,mxg2                            
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
C GENERATE JTYPE=10 INTEGRALS                                           
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
      COMMON /FQ08  / FQD(0:8),FQD0(5),FQD1(26),FQD2(48),FQD3(64),      
     2                FQD4(80),FQD5(66),FQD6(49),FQD7(24),FQD8(9)       
C$omp threadprivate(/FQ08/)
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS,                      
     2                TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2),      
     3                D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2),      
     4                NGANGB                                            
C$omp threadprivate(/GEOMPQ/)
      COMMON /MAXC  / CMAX(MXGTOT),CMAXA(MXGSH),CMAXB(MXGSH),           
     2                CMAXC(MXGSH),CMAXD(MXGSH),ISMLP(MXG2),ISMLQ       
C$omp threadprivate(/MAXC/)
      LOGICAL         LRINT                                             
      COMMON /NLRCF / LRINT                                             
C$omp threadprivate(/NLRCF /)
C                                                                       
      PARAMETER (ZER=0.0D+00)                                           
      PARAMETER (PT5=0.5D+00)                                           
      PARAMETER (ONE=1.0D+00)                                           
      PARAMETER (PI4=0.7853981633974483D+00)                            
C                                                                       
         FQD0(1)= ZER                                                   
         FQD1(1)= ZER                                                   
         FQD1(2)= ZER                                                   
         FQD2(1)= ZER                                                   
         FQD2(2)= ZER                                                   
         FQD2(3)= ZER                                                   
         FQD3(1)= ZER                                                   
         FQD3(2)= ZER                                                   
         FQD3(3)= ZER                                                   
         FQD3(4)= ZER                                                   
         FQD4(1)= ZER                                                   
         FQD4(2)= ZER                                                   
         FQD4(3)= ZER                                                   
         FQD4(4)= ZER                                                   
         FQD4(5)= ZER                                                   
C                                                                       
      DO 300 I=1,NGANGB                                                 
         ISML= ISMLQ+ISMLP(I)                                           
         IF(ISML.GE.2) GO TO 300                                        
         X12= TX12(I)                                                   
         Y02= TY02(I)                                                   
         FQZ= D00P(I)                                                   
         X41= ONE/(X12+X34)                                             
         PQR= Y02-AQZ                                                   
         PQS= PQR*PQR                                                   
         RHO= X12*X34*X41                                               
         IF(LRINT) THEN                                                 
            EFR= EMU2/(EMU2+RHO)                                        
            RHO= RHO*EFR                                                
            FQZ= FQZ*SQRT(EFR)                                          
         ENDIF                                                          
         XVA=(PQS+QPS)*RHO                                              
         RHO= RHO+RHO                                                   
         N=4                                                            
         IF(XVA.LE.TMAX) THEN                                           
C                                                                       
C     FM(T) EVALUATION...DOWNWARD RECURSION FOR M=4                     
C                                                                       
            TV= XVA*RFINC(N)                                            
            IP= NINT(TV)                                                
            FX=    FGRID(4,IP,N) *TV                                    
            FX=(FX+FGRID(3,IP,N))*TV                                    
            FX=(FX+FGRID(2,IP,N))*TV                                    
            FX=(FX+FGRID(1,IP,N))*TV                                    
            FX= FX+FGRID(0,IP,N)                                        
            TV= XVA*RXINC                                               
            IP= NINT(TV)                                                
            ET=    XGRID(4,IP) *TV                                      
            ET=(ET+XGRID(3,IP))*TV                                      
            ET=(ET+XGRID(2,IP))*TV                                      
            ET=(ET+XGRID(1,IP))*TV                                      
            ET= ET+XGRID(0,IP)                                          
C                                                                       
            FQD(N)= FX                                                  
            T2= XVA+XVA                                                 
            DO M=N,1,-1                                                 
               FQD(M-1)=(T2*FQD(M)+ET)*RMR(M)                           
            END DO                                                      
C                                                                       
            FQF= FQZ*SQRT(X41)                                          
            DO 210 M=0,N                                                
               FQD(M)= FQD(M)*FQF                                       
  210       FQF= FQF*RHO                                                
         ELSE                                                           
            XIN= ONE/XVA                                                
            FQD(0)= FQZ*SQRT(PI4*XIN*X41)                               
            ROX= RHO*XIN                                                
            FQF= PT5*ROX                                                
            DO 220 M=1,N                                                
               FQD(M)= FQD(M-1)*FQF                                     
  220       FQF= FQF+ROX                                                
         ENDIF                                                          
C                                                                       
         PQT= PQS*PQR                                                   
         PQQ= PQS*PQS                                                   
C                                                                       
         FQD0(1)= FQD0(1)+FQD(0)                                        
         FQD1(1)= FQD1(1)+FQD(1)                                        
         FQD1(2)= FQD1(2)+FQD(1)*PQR                                    
         FQD2(1)= FQD2(1)+FQD(2)                                        
         FQD2(2)= FQD2(2)+FQD(2)*PQR                                    
         FQD2(3)= FQD2(3)+FQD(2)*PQS                                    
         FQD3(1)= FQD3(1)+FQD(3)                                        
         FQD3(2)= FQD3(2)+FQD(3)*PQR                                    
         FQD3(3)= FQD3(3)+FQD(3)*PQS                                    
         FQD3(4)= FQD3(4)+FQD(3)*PQT                                    
         FQD4(1)= FQD4(1)+FQD(4)                                        
         FQD4(2)= FQD4(2)+FQD(4)*PQR                                    
         FQD4(3)= FQD4(3)+FQD(4)*PQS                                    
         FQD4(4)= FQD4(4)+FQD(4)*PQT                                    
         FQD4(5)= FQD4(5)+FQD(4)*PQQ                                    
  300 CONTINUE                                                          
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2R   *DECK INTJ11                                           
C>                                                                      
C>    @brief   DPPS case                                                
C>                                                                      
C>    @details integration of the DPPS case                             
C>                                                                      
      SUBROUTINE INTJ11                                                 
      USE lrcdft, ONLY: LCFLAG, EMU, EMU2, LRFILE                       
      use mx_limits, only: mxgtot,mxgsh,mxg2                            
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
C GENERATE JTYPE=11 INTEGRALS                                           
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
      COMMON /FQ08  / FQD(0:8),FQD0(5),FQD1(2,13),FQD2(3,16),FQD3(4,16),
     2                FQD4(80),FQD5(66),FQD6(49),FQD7(24),FQD8(9)       
C$omp threadprivate(/FQ08/)
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS,                      
     2                TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2),      
     3                D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2),      
     4                NGANGB                                            
C$omp threadprivate(/GEOMPQ/)
      COMMON /MAXC  / CMAX(MXGTOT),CMAXA(MXGSH),CMAXB(MXGSH),           
     2                CMAXC(MXGSH),CMAXD(MXGSH),ISMLP(MXG2),ISMLQ       
C$omp threadprivate(/MAXC/)
      LOGICAL         LRINT                                             
      COMMON /NLRCF / LRINT                                             
C$omp threadprivate(/NLRCF /)
C                                                                       
      DIMENSION  WORK(4,3)                                              
C                                                                       
      PARAMETER (ZER=0.0D+00)                                           
      PARAMETER (PT5=0.5D+00)                                           
      PARAMETER (ONE=1.0D+00)                                           
      PARAMETER (PI4=0.7853981633974483D+00)                            
C                                                                       
         FQD0(1)= ZER                                                   
         FQD0(2)= ZER                                                   
         DO J= 1, 3                                                     
            FQD1(1,J)= ZER                                              
            FQD1(2,J)= ZER                                              
C                                                                       
            FQD2(1,J)= ZER                                              
            FQD2(2,J)= ZER                                              
            FQD2(3,J)= ZER                                              
C                                                                       
            FQD3(1,J)= ZER                                              
            FQD3(2,J)= ZER                                              
            FQD3(3,J)= ZER                                              
            FQD3(4,J)= ZER                                              
         ENDDO                                                          
         FQD4(1)= ZER                                                   
         FQD4(2)= ZER                                                   
         FQD4(3)= ZER                                                   
         FQD4(4)= ZER                                                   
         FQD4(5)= ZER                                                   
C                                                                       
      DO 300 I=1,NGANGB                                                 
         ISML= ISMLQ+ISMLP(I)                                           
         IF(ISML.GE.2) GO TO 300                                        
         X12= TX12(I)                                                   
         Y02= TY02(I)                                                   
         FQZ= D01P(I)                                                   
         X41= ONE/(X12+X34)                                             
         PQR= Y02-AQZ                                                   
         PQS= PQR*PQR                                                   
         RHO= X12*X34*X41                                               
         IF(LRINT) THEN                                                 
            EFR= EMU2/(EMU2+RHO)                                        
            RHO= RHO*EFR                                                
            FQZ= FQZ*SQRT(EFR)                                          
         ENDIF                                                          
         XVA=(PQS+QPS)*RHO                                              
         RHO= RHO+RHO                                                   
         N=4                                                            
         IF(XVA.LE.TMAX) THEN                                           
C                                                                       
C     FM(T) EVALUATION...DOWNWARD RECURSION FOR M=4                     
C                                                                       
            TV= XVA*RFINC(N)                                            
            IP= NINT(TV)                                                
            FX=    FGRID(4,IP,N) *TV                                    
            FX=(FX+FGRID(3,IP,N))*TV                                    
            FX=(FX+FGRID(2,IP,N))*TV                                    
            FX=(FX+FGRID(1,IP,N))*TV                                    
            FX= FX+FGRID(0,IP,N)                                        
            TV= XVA*RXINC                                               
            IP= NINT(TV)                                                
            ET=    XGRID(4,IP) *TV                                      
            ET=(ET+XGRID(3,IP))*TV                                      
            ET=(ET+XGRID(2,IP))*TV                                      
            ET=(ET+XGRID(1,IP))*TV                                      
            ET= ET+XGRID(0,IP)                                          
C                                                                       
            FQD(N)= FX                                                  
            T2= XVA+XVA                                                 
            DO M=N,1,-1                                                 
               FQD(M-1)=(T2*FQD(M)+ET)*RMR(M)                           
            END DO                                                      
C                                                                       
            FQF= FQZ*SQRT(X41)                                          
            DO 210 M=0,N                                                
               FQD(M)= FQD(M)*FQF                                       
  210       FQF= FQF*RHO                                                
         ELSE                                                           
            XIN= ONE/XVA                                                
            FQD(0)= FQZ*SQRT(PI4*XIN*X41)                               
            ROX= RHO*XIN                                                
            FQF= PT5*ROX                                                
            DO 220 M=1,N                                                
               FQD(M)= FQD(M-1)*FQF                                     
  220       FQF= FQF+ROX                                                
         ENDIF                                                          
C                                                                       
         PQT = PQR*PQS                                                  
         WORK(1,1)= D00P(I)                                             
         WORK(1,2)= TY01(I)                                             
         WORK(1,3)= TX21(I)                                             
         DO J= 1, 3                                                     
            WORK(2,J)= WORK(1,J)*PQR                                    
            WORK(3,J)= WORK(1,J)*PQS                                    
            WORK(4,J)= WORK(1,J)*PQT                                    
         ENDDO                                                          
C                                                                       
         FQD0(1)= FQD0(1)+FQD(0)*WORK(1,1)                              
         FQD0(2)= FQD0(2)+FQD(0)*WORK(1,2)                              
         DO J= 1, 3                                                     
            FQD1(1,J)= FQD1(1,J)+FQD(1)*WORK(1,J)                       
            FQD1(2,J)= FQD1(2,J)+FQD(1)*WORK(2,J)                       
C                                                                       
            FQD2(1,J)= FQD2(1,J)+FQD(2)*WORK(1,J)                       
            FQD2(2,J)= FQD2(2,J)+FQD(2)*WORK(2,J)                       
            FQD2(3,J)= FQD2(3,J)+FQD(2)*WORK(3,J)                       
C                                                                       
            FQD3(1,J)= FQD3(1,J)+FQD(3)*WORK(1,J)                       
            FQD3(2,J)= FQD3(2,J)+FQD(3)*WORK(2,J)                       
            FQD3(3,J)= FQD3(3,J)+FQD(3)*WORK(3,J)                       
            FQD3(4,J)= FQD3(4,J)+FQD(3)*WORK(4,J)                       
         ENDDO                                                          
         FQD4(1)= FQD4(1)+FQD(4)*WORK(1,3)                              
         FQD4(2)= FQD4(2)+FQD(4)*WORK(2,3)                              
         FQD4(3)= FQD4(3)+FQD(4)*WORK(3,3)                              
         FQD4(4)= FQD4(4)+FQD(4)*WORK(4,3)                              
         FQD4(5)= FQD4(5)+FQD(4)*WORK(4,3)*PQR                          
  300 CONTINUE                                                          
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2R   *DECK INTJ12                                           
C>                                                                      
C>    @brief   DSDS case                                                
C>                                                                      
C>    @details integration of the DSDS case                             
C>                                                                      
      SUBROUTINE INTJ12                                                 
      USE lrcdft, ONLY: LCFLAG, EMU, EMU2, LRFILE                       
      use mx_limits, only: mxgtot,mxgsh,mxg2                            
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
C GENERATE JTYPE=12 INTEGRALS                                           
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
      COMMON /FQ08  / FQD(0:8),FQD0(5),FQD1(2,13),FQD2(3,16),FQD3(4,16),
     2                FQD4(80),FQD5(66),FQD6(49),FQD7(24),FQD8(9)       
C$omp threadprivate(/FQ08/)
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS,                      
     2                TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2),      
     3                D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2),      
     4                NGANGB                                            
C$omp threadprivate(/GEOMPQ/)
      COMMON /MAXC  / CMAX(MXGTOT),CMAXA(MXGSH),CMAXB(MXGSH),           
     2                CMAXC(MXGSH),CMAXD(MXGSH),ISMLP(MXG2),ISMLQ       
C$omp threadprivate(/MAXC/)
      LOGICAL         LRINT                                             
      COMMON /NLRCF / LRINT                                             
C$omp threadprivate(/NLRCF /)
      COMMON /SHLSPD/ CDA(MXGSH),CDB(MXGSH),CDC(MXGSH),CDD(MXGSH),      
     2                D02D(MXG2),D12D(MXG2),D22D(MXG2)                  
C$omp threadprivate(/SHLSPD/)
C                                                                       
      DIMENSION  WORK(4,4)                                              
C                                                                       
      PARAMETER (ZER=0.0D+00)                                           
      PARAMETER (PT5=0.5D+00)                                           
      PARAMETER (ONE=1.0D+00)                                           
      PARAMETER (PI4=0.7853981633974483D+00)                            
C                                                                       
         FQD0(1)= ZER                                                   
         FQD0(2)= ZER                                                   
         DO J= 1, 4                                                     
            FQD1(1,J)= ZER                                              
            FQD1(2,J)= ZER                                              
C                                                                       
            FQD2(1,J)= ZER                                              
            FQD2(2,J)= ZER                                              
            FQD2(3,J)= ZER                                              
         ENDDO                                                          
         DO J= 1, 2                                                     
            FQD3(1,J)= ZER                                              
            FQD3(2,J)= ZER                                              
            FQD3(3,J)= ZER                                              
            FQD3(4,J)= ZER                                              
         ENDDO                                                          
         FQD4(1)= ZER                                                   
         FQD4(2)= ZER                                                   
         FQD4(3)= ZER                                                   
         FQD4(4)= ZER                                                   
         FQD4(5)= ZER                                                   
C                                                                       
      DO 300 I=1,NGANGB                                                 
         ISML= ISMLQ+ISMLP(I)                                           
         IF(ISML.GE.2) GO TO 300                                        
         X12= TX12(I)                                                   
         Y02= TY02(I)                                                   
         FQZ= D02D(I)                                                   
         X41= ONE/(X12+X34)                                             
         PQR= Y02-AQZ                                                   
         PQS= PQR*PQR                                                   
         RHO= X12*X34*X41                                               
         IF(LRINT) THEN                                                 
            EFR= EMU2/(EMU2+RHO)                                        
            RHO= RHO*EFR                                                
            FQZ= FQZ*SQRT(EFR)                                          
         ENDIF                                                          
         XVA=(PQS+QPS)*RHO                                              
         RHO= RHO+RHO                                                   
         N=4                                                            
         IF(XVA.LE.TMAX) THEN                                           
C                                                                       
C     FM(T) EVALUATION...DOWNWARD RECURSION FOR M=4                     
C                                                                       
            TV= XVA*RFINC(N)                                            
            IP= NINT(TV)                                                
            FX=    FGRID(4,IP,N) *TV                                    
            FX=(FX+FGRID(3,IP,N))*TV                                    
            FX=(FX+FGRID(2,IP,N))*TV                                    
            FX=(FX+FGRID(1,IP,N))*TV                                    
            FX= FX+FGRID(0,IP,N)                                        
            TV= XVA*RXINC                                               
            IP= NINT(TV)                                                
            ET=    XGRID(4,IP) *TV                                      
            ET=(ET+XGRID(3,IP))*TV                                      
            ET=(ET+XGRID(2,IP))*TV                                      
            ET=(ET+XGRID(1,IP))*TV                                      
            ET= ET+XGRID(0,IP)                                          
C                                                                       
            FQD(N)= FX                                                  
            T2= XVA+XVA                                                 
            DO M=N,1,-1                                                 
               FQD(M-1)=(T2*FQD(M)+ET)*RMR(M)                           
            END DO                                                      
C                                                                       
            FQF= FQZ*SQRT(X41)                                          
            DO 210 M=0,N                                                
               FQD(M)= FQD(M)*FQF                                       
  210       FQF= FQF*RHO                                                
         ELSE                                                           
            XIN= ONE/XVA                                                
            FQD(0)= FQZ*SQRT(PI4*XIN*X41)                               
            ROX= RHO*XIN                                                
            FQF= PT5*ROX                                                
            DO 220 M=1,N                                                
               FQD(M)= FQD(M-1)*FQF                                     
  220       FQF= FQF+ROX                                                
         ENDIF                                                          
C                                                                       
         PQT = PQR*PQS                                                  
         XMD1= TX21(I)                                                  
         Y01 = TY01(I)                                                  
         WORK(1,1)= XMD1                                                
         WORK(1,2)= Y01 *Y01                                            
         WORK(1,3)= XMD1*Y01                                            
         WORK(1,4)= XMD1*XMD1                                           
         DO J= 1, 4                                                     
            WORK(2,J)= WORK(1,J)*PQR                                    
            WORK(3,J)= WORK(1,J)*PQS                                    
         ENDDO                                                          
         WORK(4,3)= WORK(1,3)*PQT                                       
         WORK(4,4)= WORK(1,4)*PQT                                       
C                                                                       
         FQD0(1)= FQD0(1)+FQD(0)*WORK(1,1)                              
         FQD0(2)= FQD0(2)+FQD(0)*WORK(1,2)                              
         DO J= 1, 3                                                     
            FQD1(1,J)= FQD1(1,J)+FQD(1)*WORK(1,J)                       
            FQD1(2,J)= FQD1(2,J)+FQD(1)*WORK(2,J)                       
         ENDDO                                                          
            FQD1(1,4)= FQD1(1,4)+FQD(1)*WORK(1,4)                       
         DO J= 1, 4                                                     
            FQD2(1,J)= FQD2(1,J)+FQD(2)*WORK(1,J)                       
            FQD2(2,J)= FQD2(2,J)+FQD(2)*WORK(2,J)                       
            FQD2(3,J)= FQD2(3,J)+FQD(2)*WORK(3,J)                       
         ENDDO                                                          
         DO J= 1, 2                                                     
            FQD3(1,J)= FQD3(1,J)+FQD(3)*WORK(1,J+ 2)                    
            FQD3(2,J)= FQD3(2,J)+FQD(3)*WORK(2,J+ 2)                    
            FQD3(3,J)= FQD3(3,J)+FQD(3)*WORK(3,J+ 2)                    
            FQD3(4,J)= FQD3(4,J)+FQD(3)*WORK(4,J+ 2)                    
         ENDDO                                                          
         FQD4(1)= FQD4(1)+FQD(4)*WORK(1,4)                              
         FQD4(2)= FQD4(2)+FQD(4)*WORK(2,4)                              
         FQD4(3)= FQD4(3)+FQD(4)*WORK(3,4)                              
         FQD4(4)= FQD4(4)+FQD(4)*WORK(4,4)                              
         FQD4(5)= FQD4(5)+FQD(4)*WORK(4,4)*PQR                          
  300 CONTINUE                                                          
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2R   *DECK INTJ13                                           
C>                                                                      
C>    @brief   DSPP case                                                
C>                                                                      
C>    @details integration of the DSPP case                             
C>                                                                      
      SUBROUTINE INTJ13                                                 
      USE lrcdft, ONLY: LCFLAG, EMU, EMU2, LRFILE                       
      use mx_limits, only: mxgtot,mxgsh,mxg2                            
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
C GENERATE JTYPE=13 INTEGRALS                                           
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
      COMMON /FQ08  / FQD(0:8),FQD0(5),FQD1(2,13),FQD2(3,16),FQD3(4,16),
     2                FQD4(80),FQD5(66),FQD6(49),FQD7(24),FQD8(9)       
C$omp threadprivate(/FQ08/)
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS,                      
     2                TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2),      
     3                D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2),      
     4                NGANGB                                            
C$omp threadprivate(/GEOMPQ/)
      COMMON /MAXC  / CMAX(MXGTOT),CMAXA(MXGSH),CMAXB(MXGSH),           
     2                CMAXC(MXGSH),CMAXD(MXGSH),ISMLP(MXG2),ISMLQ       
C$omp threadprivate(/MAXC/)
      LOGICAL         LRINT                                             
      COMMON /NLRCF / LRINT                                             
C$omp threadprivate(/NLRCF /)
C                                                                       
      DIMENSION  WORK(4,9)                                              
C                                                                       
      PARAMETER (ZER=0.0D+00)                                           
      PARAMETER (PT5=0.5D+00)                                           
      PARAMETER (ONE=1.0D+00)                                           
      PARAMETER (PI4=0.7853981633974483D+00)                            
C                                                                       
         DO J= 1, 5                                                     
            FQD0(J)= ZER                                                
         ENDDO                                                          
         DO J= 1, 9                                                     
            FQD1(1,J)= ZER                                              
            FQD1(2,J)= ZER                                              
C                                                                       
            FQD2(1,J)= ZER                                              
            FQD2(2,J)= ZER                                              
            FQD2(3,J)= ZER                                              
         ENDDO                                                          
         DO J= 1, 5                                                     
            FQD3(1,J)= ZER                                              
            FQD3(2,J)= ZER                                              
            FQD3(3,J)= ZER                                              
            FQD3(4,J)= ZER                                              
C                                                                       
            FQD4(J)= ZER                                                
         ENDDO                                                          
C                                                                       
      DO 300 I=1,NGANGB                                                 
         ISML= ISMLQ+ISMLP(I)                                           
         IF(ISML.GE.2) GO TO 300                                        
         X12= TX12(I)                                                   
         Y02= TY02(I)                                                   
         FQZ= D11P(I)                                                   
         X41= ONE/(X12+X34)                                             
         PQR= Y02-AQZ                                                   
         PQS= PQR*PQR                                                   
         RHO= X12*X34*X41                                               
         IF(LRINT) THEN                                                 
            EFR= EMU2/(EMU2+RHO)                                        
            RHO= RHO*EFR                                                
            FQZ= FQZ*SQRT(EFR)                                          
         ENDIF                                                          
         XVA=(PQS+QPS)*RHO                                              
         RHO= RHO+RHO                                                   
         N=4                                                            
         IF(XVA.LE.TMAX) THEN                                           
C                                                                       
C     FM(T) EVALUATION...DOWNWARD RECURSION FOR M=4                     
C                                                                       
            TV= XVA*RFINC(N)                                            
            IP= NINT(TV)                                                
            FX=    FGRID(4,IP,N) *TV                                    
            FX=(FX+FGRID(3,IP,N))*TV                                    
            FX=(FX+FGRID(2,IP,N))*TV                                    
            FX=(FX+FGRID(1,IP,N))*TV                                    
            FX= FX+FGRID(0,IP,N)                                        
            TV= XVA*RXINC                                               
            IP= NINT(TV)                                                
            ET=    XGRID(4,IP) *TV                                      
            ET=(ET+XGRID(3,IP))*TV                                      
            ET=(ET+XGRID(2,IP))*TV                                      
            ET=(ET+XGRID(1,IP))*TV                                      
            ET= ET+XGRID(0,IP)                                          
C                                                                       
            FQD(N)= FX                                                  
            T2= XVA+XVA                                                 
            DO M=N,1,-1                                                 
               FQD(M-1)=(T2*FQD(M)+ET)*RMR(M)                           
            END DO                                                      
C                                                                       
            FQF= FQZ*SQRT(X41)                                          
            DO 210 M=0,N                                                
               FQD(M)= FQD(M)*FQF                                       
  210       FQF= FQF*RHO                                                
         ELSE                                                           
            XIN= ONE/XVA                                                
            FQD(0)= FQZ*SQRT(PI4*XIN*X41)                               
            ROX= RHO*XIN                                                
            FQF= PT5*ROX                                                
            DO 220 M=1,N                                                
               FQD(M)= FQD(M-1)*FQF                                     
  220       FQF= FQF+ROX                                                
         ENDIF                                                          
C                                                                       
         PQT = PQR*PQS                                                  
         XMD1= TX21(I)                                                  
         Y01 = TY01(I)                                                  
         DP00= D00P(I)                                                  
         DP01= D01P(I)                                                  
         DP10= D10P(I)                                                  
         WORK(1,1)= DP00                                                
         WORK(1,2)= Y01 *DP01                                           
         WORK(1,3)= Y02 *DP10                                           
         WORK(1,4)= Y01 *Y02                                            
         WORK(1,5)= XMD1                                                
         WORK(1,6)= XMD1*DP01                                           
         WORK(1,7)= XMD1*DP10                                           
         WORK(1,8)= XMD1*Y01                                            
         WORK(1,9)= XMD1*XMD1                                           
         DO J= 1, 9                                                     
            WORK(2,J)= WORK(1,J)*PQR                                    
            WORK(3,J)= WORK(1,J)*PQS                                    
         ENDDO                                                          
         DO J= 5, 9                                                     
            WORK(4,J)= WORK(1,J)*PQT                                    
         ENDDO                                                          
C                                                                       
         DO J= 1, 5                                                     
            FQD0(J)= FQD0(J)+FQD(0)*WORK(1,J)                           
         ENDDO                                                          
         DO J= 1, 8                                                     
            FQD1(1,J)= FQD1(1,J)+FQD(1)*WORK(1,J)                       
            FQD1(2,J)= FQD1(2,J)+FQD(1)*WORK(2,J)                       
         ENDDO                                                          
            FQD1(1,9)= FQD1(1,9)+FQD(1)*WORK(1,9)                       
         DO J= 1, 9                                                     
            FQD2(1,J)= FQD2(1,J)+FQD(2)*WORK(1,J)                       
            FQD2(2,J)= FQD2(2,J)+FQD(2)*WORK(2,J)                       
            FQD2(3,J)= FQD2(3,J)+FQD(2)*WORK(3,J)                       
         ENDDO                                                          
         DO J= 1, 5                                                     
            FQD3(1,J)= FQD3(1,J)+FQD(3)*WORK(1,J+ 4)                    
            FQD3(2,J)= FQD3(2,J)+FQD(3)*WORK(2,J+ 4)                    
            FQD3(3,J)= FQD3(3,J)+FQD(3)*WORK(3,J+ 4)                    
            FQD3(4,J)= FQD3(4,J)+FQD(3)*WORK(4,J+ 4)                    
         ENDDO                                                          
         FQD4(1)= FQD4(1)+FQD(4)*WORK(1,9)                              
         FQD4(2)= FQD4(2)+FQD(4)*WORK(2,9)                              
         FQD4(3)= FQD4(3)+FQD(4)*WORK(3,9)                              
         FQD4(4)= FQD4(4)+FQD(4)*WORK(4,9)                              
         FQD4(5)= FQD4(5)+FQD(4)*WORK(4,9)*PQR                          
  300 CONTINUE                                                          
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2R   *DECK INTJ14                                           
C>                                                                      
C>    @brief   DDPS case                                                
C>                                                                      
C>    @details integration of the DDPS case                             
C>                                                                      
      SUBROUTINE INTJ14                                                 
      USE lrcdft, ONLY: LCFLAG, EMU, EMU2, LRFILE                       
      use mx_limits, only: mxgtot,mxgsh,mxg2                            
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
C GENERATE JTYPE=14 INTEGRALS                                           
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
      COMMON /FQ08  / FQD(0:8),FQD0(5),FQD1(2,13),FQD2(3,16),FQD3(4,16),
     2                FQD4(5,16),FQD5(66),FQD6(49),FQD7(24),FQD8(9)     
C$omp threadprivate(/FQ08/)
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS,                      
     2                TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2),      
     3                D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2),      
     4                NGANGB                                            
C$omp threadprivate(/GEOMPQ/)
      COMMON /MAXC  / CMAX(MXGTOT),CMAXA(MXGSH),CMAXB(MXGSH),           
     2                CMAXC(MXGSH),CMAXD(MXGSH),ISMLP(MXG2),ISMLQ       
C$omp threadprivate(/MAXC/)
      LOGICAL         LRINT                                             
      COMMON /NLRCF / LRINT                                             
C$omp threadprivate(/NLRCF /)
C                                                                       
      DIMENSION  WORK(5,3)                                              
C                                                                       
      PARAMETER (ZER=0.0D+00)                                           
      PARAMETER (PT5=0.5D+00)                                           
      PARAMETER (ONE=1.0D+00)                                           
      PARAMETER (PI4=0.7853981633974483D+00)                            
C                                                                       
         FQD0(1)= ZER                                                   
         FQD0(2)= ZER                                                   
         DO J= 1, 3                                                     
            FQD1(1,J)= ZER                                              
            FQD1(2,J)= ZER                                              
C                                                                       
            FQD2(1,J)= ZER                                              
            FQD2(2,J)= ZER                                              
            FQD2(3,J)= ZER                                              
C                                                                       
            FQD3(1,J)= ZER                                              
            FQD3(2,J)= ZER                                              
            FQD3(3,J)= ZER                                              
            FQD3(4,J)= ZER                                              
C                                                                       
            FQD4(1,J)= ZER                                              
            FQD4(2,J)= ZER                                              
            FQD4(3,J)= ZER                                              
            FQD4(4,J)= ZER                                              
            FQD4(5,J)= ZER                                              
         ENDDO                                                          
         DO K= 1, 6                                                     
            FQD5(K)= ZER                                                
         ENDDO                                                          
C                                                                       
      DO 300 I=1,NGANGB                                                 
         ISML= ISMLQ+ISMLP(I)                                           
         IF(ISML.GE.2) GO TO 300                                        
         X12= TX12(I)                                                   
         Y02= TY02(I)                                                   
         FQZ= D01P(I)                                                   
         X41= ONE/(X12+X34)                                             
         PQR= Y02-AQZ                                                   
         PQS= PQR*PQR                                                   
         RHO= X12*X34*X41                                               
         IF(LRINT) THEN                                                 
            EFR= EMU2/(EMU2+RHO)                                        
            RHO= RHO*EFR                                                
            FQZ= FQZ*SQRT(EFR)                                          
         ENDIF                                                          
         XVA=(PQS+QPS)*RHO                                              
         RHO= RHO+RHO                                                   
         N=5                                                            
         IF(XVA.LE.TMAX) THEN                                           
C                                                                       
C     FM(T) M=5 INTERPOLATION, GENERATING WASTED M=8,7,6 DATA           
C        FGRID(,,X) FOR  X=0,1,2,3,4,5, 6, 7 HOLDS NECESSARY DATA TO    
C        INTERPOLATE FOR M=0,1,2,3,4,8,12,16.                           
C        HERE M=5, SO WE MUST GENERATE M=8,7,6 VALUES WE DON'T USE.     
C        DOWNWARD RECURSION IS USED FOR GREATER NUMERICAL STABILITY.    
C        NOTE THAT WE ALSO USE AN INTERPOLATION FOR EXP(-T) HERE.       
C                                                                       
            TV= XVA*RFINC(N)                                            
            IP= NINT(TV)                                                
            FX=    FGRID(4,IP,N) *TV                                    
            FX=(FX+FGRID(3,IP,N))*TV                                    
            FX=(FX+FGRID(2,IP,N))*TV                                    
            FX=(FX+FGRID(1,IP,N))*TV                                    
            FX= FX+FGRID(0,IP,N)                                        
            TV= XVA*RXINC                                               
            IP= NINT(TV)                                                
            ET=    XGRID(4,IP) *TV                                      
            ET=(ET+XGRID(3,IP))*TV                                      
            ET=(ET+XGRID(2,IP))*TV                                      
            ET=(ET+XGRID(1,IP))*TV                                      
            ET= ET+XGRID(0,IP)                                          
C                                                                       
            FQD(8)= FX                                                  
            T2= XVA+XVA                                                 
            DO M=8,1,-1                                                 
               FQD(M-1)=(T2*FQD(M)+ET)*RMR(M)                           
            END DO                                                      
C                                                                       
C     THIS IS OTHER PARTS OF THE INTEGRAL, NOT FM(T)                    
C                                                                       
            FQF= FQZ*SQRT(X41)                                          
            DO 210 M=0,N                                                
               FQD(M)= FQD(M)*FQF                                       
  210       FQF= FQF*RHO                                                
         ELSE                                                           
            XIN= ONE/XVA                                                
            FQD(0)= FQZ*SQRT(PI4*XIN*X41)                               
            ROX= RHO*XIN                                                
            FQF= PT5*ROX                                                
            DO 220 M=1,N                                                
               FQD(M)= FQD(M-1)*FQF                                     
  220       FQF= FQF+ROX                                                
         ENDIF                                                          
C                                                                       
         PQT = PQR*PQS                                                  
         PQQ = PQS*PQS                                                  
         WORK(1,1)= D00P(I)                                             
         WORK(1,2)= TY01(I)                                             
         WORK(1,3)= TX21(I)                                             
         DO J= 1, 3                                                     
            WORK(2,J)= WORK(1,J)*PQR                                    
            WORK(3,J)= WORK(1,J)*PQS                                    
            WORK(4,J)= WORK(1,J)*PQT                                    
            WORK(5,J)= WORK(1,J)*PQQ                                    
         ENDDO                                                          
C                                                                       
         FQD0(1)= FQD0(1)+FQD(0)*WORK(1,1)                              
         FQD0(2)= FQD0(2)+FQD(0)*WORK(1,2)                              
         DO J= 1, 3                                                     
            FQD1(1,J)= FQD1(1,J)+FQD(1)*WORK(1,J)                       
            FQD1(2,J)= FQD1(2,J)+FQD(1)*WORK(2,J)                       
C                                                                       
            FQD2(1,J)= FQD2(1,J)+FQD(2)*WORK(1,J)                       
            FQD2(2,J)= FQD2(2,J)+FQD(2)*WORK(2,J)                       
            FQD2(3,J)= FQD2(3,J)+FQD(2)*WORK(3,J)                       
C                                                                       
            FQD3(1,J)= FQD3(1,J)+FQD(3)*WORK(1,J)                       
            FQD3(2,J)= FQD3(2,J)+FQD(3)*WORK(2,J)                       
            FQD3(3,J)= FQD3(3,J)+FQD(3)*WORK(3,J)                       
            FQD3(4,J)= FQD3(4,J)+FQD(3)*WORK(4,J)                       
C                                                                       
            FQD4(1,J)= FQD4(1,J)+FQD(4)*WORK(1,J)                       
            FQD4(2,J)= FQD4(2,J)+FQD(4)*WORK(2,J)                       
            FQD4(3,J)= FQD4(3,J)+FQD(4)*WORK(3,J)                       
            FQD4(4,J)= FQD4(4,J)+FQD(4)*WORK(4,J)                       
            FQD4(5,J)= FQD4(5,J)+FQD(4)*WORK(5,J)                       
         ENDDO                                                          
         DO K= 1, 5                                                     
            FQD5(K)= FQD5(K)+FQD(5)*WORK(K,3)                           
         ENDDO                                                          
            FQD5(6)= FQD5(6)+FQD(5)*WORK(5,3)*PQR                       
  300 CONTINUE                                                          
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2R   *DECK INTJ15                                           
C>                                                                      
C>    @brief   DPDS case                                                
C>                                                                      
C>    @details integration of the DPDS case                             
C>                                                                      
      SUBROUTINE INTJ15                                                 
      USE lrcdft, ONLY: LCFLAG, EMU, EMU2, LRFILE                       
      use mx_limits, only: mxgtot,mxgsh,mxg2                            
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
C GENERATE JTYPE=15 INTEGRALS                                           
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
      COMMON /FQ08  / FQD(0:8),FQD0(5),FQD1(2,13),FQD2(3,16),FQD3(4,16),
     2                FQD4(5,16),FQD5(66),FQD6(49),FQD7(24),FQD8(9)     
C$omp threadprivate(/FQ08/)
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS,                      
     2                TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2),      
     3                D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2),      
     4                NGANGB                                            
C$omp threadprivate(/GEOMPQ/)
      COMMON /MAXC  / CMAX(MXGTOT),CMAXA(MXGSH),CMAXB(MXGSH),           
     2                CMAXC(MXGSH),CMAXD(MXGSH),ISMLP(MXG2),ISMLQ       
C$omp threadprivate(/MAXC/)
      LOGICAL         LRINT                                             
      COMMON /NLRCF / LRINT                                             
C$omp threadprivate(/NLRCF /)
      COMMON /SHLSPD/ CDA(MXGSH),CDB(MXGSH),CDC(MXGSH),CDD(MXGSH),      
     2                D02D(MXG2),D12D(MXG2),D22D(MXG2)                  
C$omp threadprivate(/SHLSPD/)
C                                                                       
      DIMENSION  WORK(5,4)                                              
C                                                                       
      PARAMETER (ZER=0.0D+00)                                           
      PARAMETER (PT5=0.5D+00)                                           
      PARAMETER (ONE=1.0D+00)                                           
      PARAMETER (PI4=0.7853981633974483D+00)                            
C                                                                       
         FQD0(1)= ZER                                                   
         FQD0(2)= ZER                                                   
         DO J= 1, 4                                                     
            FQD1(1,J)= ZER                                              
            FQD1(2,J)= ZER                                              
C                                                                       
            FQD2(1,J)= ZER                                              
            FQD2(2,J)= ZER                                              
            FQD2(3,J)= ZER                                              
C                                                                       
            FQD3(1,J)= ZER                                              
            FQD3(2,J)= ZER                                              
            FQD3(3,J)= ZER                                              
            FQD3(4,J)= ZER                                              
         ENDDO                                                          
         DO J= 1, 2                                                     
            FQD4(1,J)= ZER                                              
            FQD4(2,J)= ZER                                              
            FQD4(3,J)= ZER                                              
            FQD4(4,J)= ZER                                              
            FQD4(5,J)= ZER                                              
         ENDDO                                                          
         DO K= 1, 6                                                     
            FQD5(K)= ZER                                                
         ENDDO                                                          
C                                                                       
      DO 300 I=1,NGANGB                                                 
         ISML= ISMLQ+ISMLP(I)                                           
         IF(ISML.GE.2) GO TO 300                                        
         X12= TX12(I)                                                   
         Y02= TY02(I)                                                   
         FQZ= D02D(I)                                                   
         X41= ONE/(X12+X34)                                             
         PQR= Y02-AQZ                                                   
         PQS= PQR*PQR                                                   
         RHO= X12*X34*X41                                               
         IF(LRINT) THEN                                                 
            EFR= EMU2/(EMU2+RHO)                                        
            RHO= RHO*EFR                                                
            FQZ= FQZ*SQRT(EFR)                                          
         ENDIF                                                          
         XVA=(PQS+QPS)*RHO                                              
         RHO= RHO+RHO                                                   
         N=5                                                            
         IF(XVA.LE.TMAX) THEN                                           
C                                                                       
C     FM(T) M=5 INTERPOLATION, GENERATING WASTED M=8,7,6 DATA           
C                                                                       
            TV= XVA*RFINC(N)                                            
            IP= NINT(TV)                                                
            FX=    FGRID(4,IP,N) *TV                                    
            FX=(FX+FGRID(3,IP,N))*TV                                    
            FX=(FX+FGRID(2,IP,N))*TV                                    
            FX=(FX+FGRID(1,IP,N))*TV                                    
            FX= FX+FGRID(0,IP,N)                                        
            TV= XVA*RXINC                                               
            IP= NINT(TV)                                                
            ET=    XGRID(4,IP) *TV                                      
            ET=(ET+XGRID(3,IP))*TV                                      
            ET=(ET+XGRID(2,IP))*TV                                      
            ET=(ET+XGRID(1,IP))*TV                                      
            ET= ET+XGRID(0,IP)                                          
C                                                                       
            FQD(8)= FX                                                  
            T2= XVA+XVA                                                 
            DO M=8,1,-1                                                 
               FQD(M-1)=(T2*FQD(M)+ET)*RMR(M)                           
            END DO                                                      
C                                                                       
            FQF= FQZ*SQRT(X41)                                          
            DO 210 M=0,N                                                
               FQD(M)= FQD(M)*FQF                                       
  210       FQF= FQF*RHO                                                
         ELSE                                                           
            XIN= ONE/XVA                                                
            FQD(0)= FQZ*SQRT(PI4*XIN*X41)                               
            ROX= RHO*XIN                                                
            FQF= PT5*ROX                                                
            DO 220 M=1,N                                                
               FQD(M)= FQD(M-1)*FQF                                     
  220       FQF= FQF+ROX                                                
         ENDIF                                                          
C                                                                       
         PQT = PQR*PQS                                                  
         XMD1= TX21(I)                                                  
         Y01 = TY01(I)                                                  
         WORK(1,1)= XMD1                                                
         WORK(1,2)= Y01 *Y01                                            
         WORK(1,3)= XMD1*Y01                                            
         WORK(1,4)= XMD1*XMD1                                           
         DO J= 1, 4                                                     
            WORK(2,J)= WORK(1,J)*PQR                                    
            WORK(3,J)= WORK(1,J)*PQS                                    
            WORK(4,J)= WORK(1,J)*PQT                                    
         ENDDO                                                          
         WORK(5,3)= WORK(4,3)*PQR                                       
         WORK(5,4)= WORK(4,4)*PQR                                       
C                                                                       
         FQD0(1)= FQD0(1)+FQD(0)*WORK(1,1)                              
         FQD0(2)= FQD0(2)+FQD(0)*WORK(1,2)                              
         DO J= 1, 3                                                     
            FQD1(1,J)= FQD1(1,J)+FQD(1)*WORK(1,J)                       
            FQD1(2,J)= FQD1(2,J)+FQD(1)*WORK(2,J)                       
         ENDDO                                                          
            FQD1(1,4)= FQD1(1,4)+FQD(1)*WORK(1,4)                       
         DO J= 1, 4                                                     
            FQD2(1,J)= FQD2(1,J)+FQD(2)*WORK(1,J)                       
            FQD2(2,J)= FQD2(2,J)+FQD(2)*WORK(2,J)                       
            FQD2(3,J)= FQD2(3,J)+FQD(2)*WORK(3,J)                       
C                                                                       
            FQD3(1,J)= FQD3(1,J)+FQD(3)*WORK(1,J)                       
            FQD3(2,J)= FQD3(2,J)+FQD(3)*WORK(2,J)                       
            FQD3(3,J)= FQD3(3,J)+FQD(3)*WORK(3,J)                       
            FQD3(4,J)= FQD3(4,J)+FQD(3)*WORK(4,J)                       
         ENDDO                                                          
         DO J= 1, 2                                                     
            FQD4(1,J)= FQD4(1,J)+FQD(4)*WORK(1,J+ 2)                    
            FQD4(2,J)= FQD4(2,J)+FQD(4)*WORK(2,J+ 2)                    
            FQD4(3,J)= FQD4(3,J)+FQD(4)*WORK(3,J+ 2)                    
            FQD4(4,J)= FQD4(4,J)+FQD(4)*WORK(4,J+ 2)                    
            FQD4(5,J)= FQD4(5,J)+FQD(4)*WORK(5,J+ 2)                    
         ENDDO                                                          
         DO K= 1, 5                                                     
            FQD5(K)= FQD5(K)+FQD(5)*WORK(K,4)                           
         ENDDO                                                          
            FQD5(6)= FQD5(6)+FQD(5)*WORK(5,4)*PQR                       
  300 CONTINUE                                                          
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2R   *DECK INTJ16                                           
C>                                                                      
C>    @brief   DPPP case                                                
C>                                                                      
C>    @details integration of the DPPP case                             
C>                                                                      
      SUBROUTINE INTJ16                                                 
      USE lrcdft, ONLY: LCFLAG, EMU, EMU2, LRFILE                       
      use mx_limits, only: mxgtot,mxgsh,mxg2                            
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
C GENERATE JTYPE=16 INTEGRALS                                           
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
      COMMON /FQ08  / FQD(0:8),FQD0(5),FQD1(2,13),FQD2(3,16),FQD3(4,16),
     2                FQD4(5,16),FQD5(66),FQD6(49),FQD7(24),FQD8(9)     
C$omp threadprivate(/FQ08/)
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS,                      
     2                TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2),      
     3                D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2),      
     4                NGANGB                                            
C$omp threadprivate(/GEOMPQ/)
      COMMON /MAXC  / CMAX(MXGTOT),CMAXA(MXGSH),CMAXB(MXGSH),           
     2                CMAXC(MXGSH),CMAXD(MXGSH),ISMLP(MXG2),ISMLQ       
C$omp threadprivate(/MAXC/)
      LOGICAL         LRINT                                             
      COMMON /NLRCF / LRINT                                             
C$omp threadprivate(/NLRCF /)
C                                                                       
      DIMENSION  WORK(5,9)                                              
C                                                                       
      PARAMETER (ZER=0.0D+00)                                           
      PARAMETER (PT5=0.5D+00)                                           
      PARAMETER (ONE=1.0D+00)                                           
      PARAMETER (PI4=0.7853981633974483D+00)                            
C                                                                       
         DO J= 1, 5                                                     
            FQD0(J)= ZER                                                
         ENDDO                                                          
         DO J= 1, 9                                                     
            FQD1(1,J)= ZER                                              
            FQD1(2,J)= ZER                                              
C                                                                       
            FQD2(1,J)= ZER                                              
            FQD2(2,J)= ZER                                              
            FQD2(3,J)= ZER                                              
C                                                                       
            FQD3(1,J)= ZER                                              
            FQD3(2,J)= ZER                                              
            FQD3(3,J)= ZER                                              
            FQD3(4,J)= ZER                                              
         ENDDO                                                          
         DO J= 1, 5                                                     
            FQD4(1,J)= ZER                                              
            FQD4(2,J)= ZER                                              
            FQD4(3,J)= ZER                                              
            FQD4(4,J)= ZER                                              
            FQD4(5,J)= ZER                                              
         ENDDO                                                          
         DO K= 1, 6                                                     
            FQD5(K)= ZER                                                
         ENDDO                                                          
C                                                                       
      DO 300 I=1,NGANGB                                                 
         ISML= ISMLQ+ISMLP(I)                                           
         IF(ISML.GE.2) GO TO 300                                        
         X12= TX12(I)                                                   
         Y02= TY02(I)                                                   
         FQZ= D11P(I)                                                   
         X41= ONE/(X12+X34)                                             
         PQR= Y02-AQZ                                                   
         PQS= PQR*PQR                                                   
         RHO= X12*X34*X41                                               
         IF(LRINT) THEN                                                 
            EFR= EMU2/(EMU2+RHO)                                        
            RHO= RHO*EFR                                                
            FQZ= FQZ*SQRT(EFR)                                          
         ENDIF                                                          
         XVA=(PQS+QPS)*RHO                                              
         RHO= RHO+RHO                                                   
         N=5                                                            
         IF(XVA.LE.TMAX) THEN                                           
C                                                                       
C     FM(T) M=5 INTERPOLATION, GENERATING WASTED M=8,7,6 DATA           
C                                                                       
            TV= XVA*RFINC(N)                                            
            IP= NINT(TV)                                                
            FX=    FGRID(4,IP,N) *TV                                    
            FX=(FX+FGRID(3,IP,N))*TV                                    
            FX=(FX+FGRID(2,IP,N))*TV                                    
            FX=(FX+FGRID(1,IP,N))*TV                                    
            FX= FX+FGRID(0,IP,N)                                        
            TV= XVA*RXINC                                               
            IP= NINT(TV)                                                
            ET=    XGRID(4,IP) *TV                                      
            ET=(ET+XGRID(3,IP))*TV                                      
            ET=(ET+XGRID(2,IP))*TV                                      
            ET=(ET+XGRID(1,IP))*TV                                      
            ET= ET+XGRID(0,IP)                                          
C                                                                       
            FQD(8)= FX                                                  
            T2= XVA+XVA                                                 
            DO M=8,1,-1                                                 
               FQD(M-1)=(T2*FQD(M)+ET)*RMR(M)                           
            END DO                                                      
C                                                                       
            FQF= FQZ*SQRT(X41)                                          
            DO 210 M=0,N                                                
               FQD(M)= FQD(M)*FQF                                       
  210       FQF= FQF*RHO                                                
         ELSE                                                           
            XIN= ONE/XVA                                                
            FQD(0)= FQZ*SQRT(PI4*XIN*X41)                               
            ROX= RHO*XIN                                                
            FQF= PT5*ROX                                                
            DO 220 M=1,N                                                
               FQD(M)= FQD(M-1)*FQF                                     
  220       FQF= FQF+ROX                                                
         ENDIF                                                          
C                                                                       
         PQT = PQR*PQS                                                  
         PQQ = PQS*PQS                                                  
         XMD1= TX21(I)                                                  
         Y01 = TY01(I)                                                  
         DP00= D00P(I)                                                  
         DP01= D01P(I)                                                  
         DP10= D10P(I)                                                  
         WORK(1,1)= DP00                                                
         WORK(1,2)= Y01 *DP01                                           
         WORK(1,3)= Y02 *DP10                                           
         WORK(1,4)= Y01 *Y02                                            
         WORK(1,5)= XMD1                                                
         WORK(1,6)= XMD1*DP01                                           
         WORK(1,7)= XMD1*DP10                                           
         WORK(1,8)= XMD1*Y01                                            
         WORK(1,9)= XMD1*XMD1                                           
         DO J= 1, 9                                                     
            WORK(2,J)= WORK(1,J)*PQR                                    
            WORK(3,J)= WORK(1,J)*PQS                                    
            WORK(4,J)= WORK(1,J)*PQT                                    
         ENDDO                                                          
         DO J= 5, 9                                                     
            WORK(5,J)= WORK(1,J)*PQQ                                    
         ENDDO                                                          
C                                                                       
         DO J= 1, 5                                                     
            FQD0(J)= FQD0(J)+FQD(0)*WORK(1,J)                           
         ENDDO                                                          
         DO J= 1, 8                                                     
            FQD1(1,J)= FQD1(1,J)+FQD(1)*WORK(1,J)                       
            FQD1(2,J)= FQD1(2,J)+FQD(1)*WORK(2,J)                       
         ENDDO                                                          
            FQD1(1,9)= FQD1(1,9)+FQD(1)*WORK(1,9)                       
         DO J= 1, 9                                                     
            FQD2(1,J)= FQD2(1,J)+FQD(2)*WORK(1,J)                       
            FQD2(2,J)= FQD2(2,J)+FQD(2)*WORK(2,J)                       
            FQD2(3,J)= FQD2(3,J)+FQD(2)*WORK(3,J)                       
C                                                                       
            FQD3(1,J)= FQD3(1,J)+FQD(3)*WORK(1,J)                       
            FQD3(2,J)= FQD3(2,J)+FQD(3)*WORK(2,J)                       
            FQD3(3,J)= FQD3(3,J)+FQD(3)*WORK(3,J)                       
            FQD3(4,J)= FQD3(4,J)+FQD(3)*WORK(4,J)                       
         ENDDO                                                          
         DO J= 1, 5                                                     
            FQD4(1,J)= FQD4(1,J)+FQD(4)*WORK(1,J+ 4)                    
            FQD4(2,J)= FQD4(2,J)+FQD(4)*WORK(2,J+ 4)                    
            FQD4(3,J)= FQD4(3,J)+FQD(4)*WORK(3,J+ 4)                    
            FQD4(4,J)= FQD4(4,J)+FQD(4)*WORK(4,J+ 4)                    
            FQD4(5,J)= FQD4(5,J)+FQD(4)*WORK(5,J+ 4)                    
         ENDDO                                                          
         DO K= 1, 5                                                     
            FQD5(K)= FQD5(K)+FQD(5)*WORK(K,9)                           
         ENDDO                                                          
            FQD5(6)= FQD5(6)+FQD(5)*WORK(5,9)*PQR                       
  300 CONTINUE                                                          
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2R   *DECK INTJ17                                           
C>                                                                      
C>    @brief   DDDS case                                                
C>                                                                      
C>    @details integration of the DDDS case                             
C>                                                                      
      SUBROUTINE INTJ17                                                 
      USE lrcdft, ONLY: LCFLAG, EMU, EMU2, LRFILE                       
      use mx_limits, only: mxgtot,mxgsh,mxg2                            
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
C GENERATE JTYPE=17 INTEGRALS                                           
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
      COMMON /FQ08  / FQD(0:8),FQD0(5),FQD1(2,13),FQD2(3,16),FQD3(4,16),
     2                FQD4(5,16),FQD5(6,11),FQD6(49),FQD7(24),FQD8(9)   
C$omp threadprivate(/FQ08/)
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS,                      
     2                TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2),      
     3                D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2),      
     4                NGANGB                                            
C$omp threadprivate(/GEOMPQ/)
      COMMON /MAXC  / CMAX(MXGTOT),CMAXA(MXGSH),CMAXB(MXGSH),           
     2                CMAXC(MXGSH),CMAXD(MXGSH),ISMLP(MXG2),ISMLQ       
C$omp threadprivate(/MAXC/)
      LOGICAL         LRINT                                             
      COMMON /NLRCF / LRINT                                             
C$omp threadprivate(/NLRCF /)
      COMMON /SHLSPD/ CDA(MXGSH),CDB(MXGSH),CDC(MXGSH),CDD(MXGSH),      
     2                D02D(MXG2),D12D(MXG2),D22D(MXG2)                  
C$omp threadprivate(/SHLSPD/)
C                                                                       
      DIMENSION  WORK(6,4)                                              
C                                                                       
      PARAMETER (ZER=0.0D+00)                                           
      PARAMETER (PT5=0.5D+00)                                           
      PARAMETER (ONE=1.0D+00)                                           
      PARAMETER (PI4=0.7853981633974483D+00)                            
C                                                                       
         FQD0(1)= ZER                                                   
         FQD0(2)= ZER                                                   
         DO J= 1, 4                                                     
            FQD1(1,J)= ZER                                              
            FQD1(2,J)= ZER                                              
C                                                                       
            FQD2(1,J)= ZER                                              
            FQD2(2,J)= ZER                                              
            FQD2(3,J)= ZER                                              
C                                                                       
            FQD3(1,J)= ZER                                              
            FQD3(2,J)= ZER                                              
            FQD3(3,J)= ZER                                              
            FQD3(4,J)= ZER                                              
C                                                                       
            FQD4(1,J)= ZER                                              
            FQD4(2,J)= ZER                                              
            FQD4(3,J)= ZER                                              
            FQD4(4,J)= ZER                                              
            FQD4(5,J)= ZER                                              
         ENDDO                                                          
         DO J= 1, 2                                                     
            DO K= 1, 6                                                  
               FQD5(K,J)= ZER                                           
            ENDDO                                                       
         ENDDO                                                          
         DO K= 1, 7                                                     
            FQD6(K)= ZER                                                
         ENDDO                                                          
C                                                                       
      DO 300 I=1,NGANGB                                                 
         ISML= ISMLQ+ISMLP(I)                                           
         IF(ISML.GE.2) GO TO 300                                        
         X12= TX12(I)                                                   
         Y02= TY02(I)                                                   
         FQZ= D02D(I)                                                   
         X41= ONE/(X12+X34)                                             
         PQR= Y02-AQZ                                                   
         PQS= PQR*PQR                                                   
         RHO= X12*X34*X41                                               
         IF(LRINT) THEN                                                 
            EFR= EMU2/(EMU2+RHO)                                        
            RHO= RHO*EFR                                                
            FQZ= FQZ*SQRT(EFR)                                          
         ENDIF                                                          
         XVA=(PQS+QPS)*RHO                                              
         RHO= RHO+RHO                                                   
         N=6                                                            
         IF(XVA.LE.TMAX) THEN                                           
C                                                                       
C     FM(T) M=6 INTERPOLATION, GENERATING WASTED M=8,7 DATA             
C                                                                       
            M=5                                                         
            TV= XVA*RFINC(M)                                            
            IP= NINT(TV)                                                
            FX=    FGRID(4,IP,M) *TV                                    
            FX=(FX+FGRID(3,IP,M))*TV                                    
            FX=(FX+FGRID(2,IP,M))*TV                                    
            FX=(FX+FGRID(1,IP,M))*TV                                    
            FX= FX+FGRID(0,IP,M)                                        
            TV= XVA*RXINC                                               
            IP= NINT(TV)                                                
            ET=    XGRID(4,IP) *TV                                      
            ET=(ET+XGRID(3,IP))*TV                                      
            ET=(ET+XGRID(2,IP))*TV                                      
            ET=(ET+XGRID(1,IP))*TV                                      
            ET= ET+XGRID(0,IP)                                          
C                                                                       
            FQD(8)= FX                                                  
            T2= XVA+XVA                                                 
            DO M=8,1,-1                                                 
               FQD(M-1)=(T2*FQD(M)+ET)*RMR(M)                           
            END DO                                                      
C                                                                       
            FQF= FQZ*SQRT(X41)                                          
            DO 210 M=0,N                                                
               FQD(M)= FQD(M)*FQF                                       
  210       FQF= FQF*RHO                                                
         ELSE                                                           
            XIN= ONE/XVA                                                
            FQD(0)= FQZ*SQRT(PI4*XIN*X41)                               
            ROX= RHO*XIN                                                
            FQF= PT5*ROX                                                
            DO 220 M=1,N                                                
               FQD(M)= FQD(M-1)*FQF                                     
  220       FQF= FQF+ROX                                                
         ENDIF                                                          
C                                                                       
         PQT = PQR*PQS                                                  
         PQQ = PQS*PQS                                                  
         XMD1= TX21(I)                                                  
         Y01 = TY01(I)                                                  
         WORK(1,1)= XMD1                                                
         WORK(1,2)= Y01 *Y01                                            
         WORK(1,3)= XMD1*Y01                                            
         WORK(1,4)= XMD1*XMD1                                           
         DO J= 1, 4                                                     
            WORK(2,J)= WORK(1,J)*PQR                                    
            WORK(3,J)= WORK(1,J)*PQS                                    
            WORK(4,J)= WORK(1,J)*PQT                                    
            WORK(5,J)= WORK(1,J)*PQQ                                    
         ENDDO                                                          
         WORK(6,3)= WORK(5,3)*PQR                                       
         WORK(6,4)= WORK(5,4)*PQR                                       
C                                                                       
         FQD0(1)= FQD0(1)+FQD(0)*WORK(1,1)                              
         FQD0(2)= FQD0(2)+FQD(0)*WORK(1,2)                              
         DO J= 1, 3                                                     
            FQD1(1,J)= FQD1(1,J)+FQD(1)*WORK(1,J)                       
            FQD1(2,J)= FQD1(2,J)+FQD(1)*WORK(2,J)                       
         ENDDO                                                          
            FQD1(1,4)= FQD1(1,4)+FQD(1)*WORK(1,4)                       
         DO J= 1, 4                                                     
            FQD2(1,J)= FQD2(1,J)+FQD(2)*WORK(1,J)                       
            FQD2(2,J)= FQD2(2,J)+FQD(2)*WORK(2,J)                       
            FQD2(3,J)= FQD2(3,J)+FQD(2)*WORK(3,J)                       
C                                                                       
            FQD3(1,J)= FQD3(1,J)+FQD(3)*WORK(1,J)                       
            FQD3(2,J)= FQD3(2,J)+FQD(3)*WORK(2,J)                       
            FQD3(3,J)= FQD3(3,J)+FQD(3)*WORK(3,J)                       
            FQD3(4,J)= FQD3(4,J)+FQD(3)*WORK(4,J)                       
C                                                                       
            FQD4(1,J)= FQD4(1,J)+FQD(4)*WORK(1,J)                       
            FQD4(2,J)= FQD4(2,J)+FQD(4)*WORK(2,J)                       
            FQD4(3,J)= FQD4(3,J)+FQD(4)*WORK(3,J)                       
            FQD4(4,J)= FQD4(4,J)+FQD(4)*WORK(4,J)                       
            FQD4(5,J)= FQD4(5,J)+FQD(4)*WORK(5,J)                       
         ENDDO                                                          
         DO J= 1, 2                                                     
            DO K= 1, 6                                                  
               FQD5(K,J)= FQD5(K,J)+FQD(5)*WORK(K,J+ 2)                 
            ENDDO                                                       
         ENDDO                                                          
         DO K= 1, 6                                                     
            FQD6(K)= FQD6(K)+FQD(6)*WORK(K,4)                           
         ENDDO                                                          
            FQD6(7)= FQD6(7)+FQD(6)*WORK(6,4)*PQR                       
  300 CONTINUE                                                          
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2R   *DECK INTJ18                                           
C>                                                                      
C>    @brief   DDPP case                                                
C>                                                                      
C>    @details integration of the DDPP case                             
C>                                                                      
      SUBROUTINE INTJ18                                                 
      USE lrcdft, ONLY: LCFLAG, EMU, EMU2, LRFILE                       
      use mx_limits, only: mxgtot,mxgsh,mxg2                            
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
C GENERATE JTYPE=18 INTEGRALS                                           
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
      COMMON /FQ08  / FQD(0:8),FQD0(5),FQD1(2,13),FQD2(3,16),FQD3(4,16),
     2                FQD4(5,16),FQD5(6,11),FQD6(49),FQD7(24),FQD8(9)   
C$omp threadprivate(/FQ08/)
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS,                      
     2                TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2),      
     3                D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2),      
     4                NGANGB                                            
C$omp threadprivate(/GEOMPQ/)
      COMMON /MAXC  / CMAX(MXGTOT),CMAXA(MXGSH),CMAXB(MXGSH),           
     2                CMAXC(MXGSH),CMAXD(MXGSH),ISMLP(MXG2),ISMLQ       
C$omp threadprivate(/MAXC/)
      LOGICAL         LRINT                                             
      COMMON /NLRCF / LRINT                                             
C$omp threadprivate(/NLRCF /)
C                                                                       
      DIMENSION  WORK(6,9)                                              
C                                                                       
      PARAMETER (ZER=0.0D+00)                                           
      PARAMETER (PT5=0.5D+00)                                           
      PARAMETER (ONE=1.0D+00)                                           
      PARAMETER (PI4=0.7853981633974483D+00)                            
C                                                                       
         DO J= 1, 5                                                     
            FQD0(J)= ZER                                                
         ENDDO                                                          
         DO J= 1, 9                                                     
            FQD1(1,J)= ZER                                              
            FQD1(2,J)= ZER                                              
C                                                                       
            FQD2(1,J)= ZER                                              
            FQD2(2,J)= ZER                                              
            FQD2(3,J)= ZER                                              
C                                                                       
            FQD3(1,J)= ZER                                              
            FQD3(2,J)= ZER                                              
            FQD3(3,J)= ZER                                              
            FQD3(4,J)= ZER                                              
C                                                                       
            FQD4(1,J)= ZER                                              
            FQD4(2,J)= ZER                                              
            FQD4(3,J)= ZER                                              
            FQD4(4,J)= ZER                                              
            FQD4(5,J)= ZER                                              
         ENDDO                                                          
         DO J= 1, 5                                                     
            DO K= 1, 6                                                  
               FQD5(K,J)= ZER                                           
            ENDDO                                                       
         ENDDO                                                          
         DO K= 1, 7                                                     
            FQD6(K)= ZER                                                
         ENDDO                                                          
C                                                                       
      DO 300 I=1,NGANGB                                                 
         ISML= ISMLQ+ISMLP(I)                                           
         IF(ISML.GE.2) GO TO 300                                        
         X12= TX12(I)                                                   
         Y02= TY02(I)                                                   
         FQZ= D11P(I)                                                   
         X41= ONE/(X12+X34)                                             
         PQR= Y02-AQZ                                                   
         PQS= PQR*PQR                                                   
         RHO= X12*X34*X41                                               
         IF(LRINT) THEN                                                 
            EFR= EMU2/(EMU2+RHO)                                        
            RHO= RHO*EFR                                                
            FQZ= FQZ*SQRT(EFR)                                          
         ENDIF                                                          
         XVA=(PQS+QPS)*RHO                                              
         RHO= RHO+RHO                                                   
         N=6                                                            
         IF(XVA.LE.TMAX) THEN                                           
C                                                                       
C     FM(T) M=6 INTERPOLATION, GENERATING WASTED M=8,7 DATA             
C                                                                       
            M=5                                                         
            TV= XVA*RFINC(M)                                            
            IP= NINT(TV)                                                
            FX=    FGRID(4,IP,M) *TV                                    
            FX=(FX+FGRID(3,IP,M))*TV                                    
            FX=(FX+FGRID(2,IP,M))*TV                                    
            FX=(FX+FGRID(1,IP,M))*TV                                    
            FX= FX+FGRID(0,IP,M)                                        
            TV= XVA*RXINC                                               
            IP= NINT(TV)                                                
            ET=    XGRID(4,IP) *TV                                      
            ET=(ET+XGRID(3,IP))*TV                                      
            ET=(ET+XGRID(2,IP))*TV                                      
            ET=(ET+XGRID(1,IP))*TV                                      
            ET= ET+XGRID(0,IP)                                          
C                                                                       
            FQD(8)= FX                                                  
            T2= XVA+XVA                                                 
            DO M=8,1,-1                                                 
               FQD(M-1)=(T2*FQD(M)+ET)*RMR(M)                           
            END DO                                                      
C                                                                       
            FQF= FQZ*SQRT(X41)                                          
            DO 210 M=0,N                                                
               FQD(M)= FQD(M)*FQF                                       
  210       FQF= FQF*RHO                                                
         ELSE                                                           
            XIN= ONE/XVA                                                
            FQD(0)= FQZ*SQRT(PI4*XIN*X41)                               
            ROX= RHO*XIN                                                
            FQF= PT5*ROX                                                
            DO 220 M=1,N                                                
               FQD(M)= FQD(M-1)*FQF                                     
  220       FQF= FQF+ROX                                                
         ENDIF                                                          
C                                                                       
         PQT = PQR*PQS                                                  
         PQQ = PQS*PQS                                                  
         PQ5 = PQT*PQS                                                  
         XMD1= TX21(I)                                                  
         Y01 = TY01(I)                                                  
         DP00= D00P(I)                                                  
         DP01= D01P(I)                                                  
         DP10= D10P(I)                                                  
         WORK(1,1)= DP00                                                
         WORK(1,2)= Y01 *DP01                                           
         WORK(1,3)= Y02 *DP10                                           
         WORK(1,4)= Y01 *Y02                                            
         WORK(1,5)= XMD1                                                
         WORK(1,6)= XMD1*DP01                                           
         WORK(1,7)= XMD1*DP10                                           
         WORK(1,8)= XMD1*Y01                                            
         WORK(1,9)= XMD1*XMD1                                           
         DO J= 1, 9                                                     
            WORK(2,J)= WORK(1,J)*PQR                                    
            WORK(3,J)= WORK(1,J)*PQS                                    
            WORK(4,J)= WORK(1,J)*PQT                                    
            WORK(5,J)= WORK(1,J)*PQQ                                    
         ENDDO                                                          
         DO J= 5, 9                                                     
            WORK(6,J)= WORK(1,J)*PQ5                                    
         ENDDO                                                          
C                                                                       
         DO J= 1, 5                                                     
            FQD0(J)= FQD0(J)+FQD(0)*WORK(1,J)                           
         ENDDO                                                          
         DO J= 1, 8                                                     
            FQD1(1,J)= FQD1(1,J)+FQD(1)*WORK(1,J)                       
            FQD1(2,J)= FQD1(2,J)+FQD(1)*WORK(2,J)                       
         ENDDO                                                          
            FQD1(1,9)= FQD1(1,9)+FQD(1)*WORK(1,9)                       
         DO J= 1, 9                                                     
            FQD2(1,J)= FQD2(1,J)+FQD(2)*WORK(1,J)                       
            FQD2(2,J)= FQD2(2,J)+FQD(2)*WORK(2,J)                       
            FQD2(3,J)= FQD2(3,J)+FQD(2)*WORK(3,J)                       
C                                                                       
            FQD3(1,J)= FQD3(1,J)+FQD(3)*WORK(1,J)                       
            FQD3(2,J)= FQD3(2,J)+FQD(3)*WORK(2,J)                       
            FQD3(3,J)= FQD3(3,J)+FQD(3)*WORK(3,J)                       
            FQD3(4,J)= FQD3(4,J)+FQD(3)*WORK(4,J)                       
C                                                                       
            FQD4(1,J)= FQD4(1,J)+FQD(4)*WORK(1,J)                       
            FQD4(2,J)= FQD4(2,J)+FQD(4)*WORK(2,J)                       
            FQD4(3,J)= FQD4(3,J)+FQD(4)*WORK(3,J)                       
            FQD4(4,J)= FQD4(4,J)+FQD(4)*WORK(4,J)                       
            FQD4(5,J)= FQD4(5,J)+FQD(4)*WORK(5,J)                       
         ENDDO                                                          
         DO J= 1, 5                                                     
            DO K= 1, 6                                                  
               FQD5(K,J)= FQD5(K,J)+FQD(5)*WORK(K,J+ 4)                 
            ENDDO                                                       
         ENDDO                                                          
         DO K= 1, 6                                                     
            FQD6(K)= FQD6(K)+FQD(6)*WORK(K,9)                           
         ENDDO                                                          
            FQD6(7)= FQD6(7)+FQD(6)*WORK(6,9)*PQR                       
  300 CONTINUE                                                          
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2R   *DECK INTJ19                                           
C>                                                                      
C>    @brief   DPDP case                                                
C>                                                                      
C>    @details integration of the DPDP case                             
C>                                                                      
      SUBROUTINE INTJ19                                                 
      USE lrcdft, ONLY: LCFLAG, EMU, EMU2, LRFILE                       
      use mx_limits, only: mxgtot,mxgsh,mxg2                            
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
C GENERATE JTYPE=19 INTEGRALS                                           
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
      COMMON /FQ08  / FQD(0:8),FQD0(5),FQD1(2,13),FQD2(3,16),FQD3(4,16),
     2                FQD4(5,16),FQD5(6,11),FQD6(49),FQD7(24),FQD8(9)   
C$omp threadprivate(/FQ08/)
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS,                      
     2                TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2),      
     3                D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2),      
     4                NGANGB                                            
C$omp threadprivate(/GEOMPQ/)
      COMMON /MAXC  / CMAX(MXGTOT),CMAXA(MXGSH),CMAXB(MXGSH),           
     2                CMAXC(MXGSH),CMAXD(MXGSH),ISMLP(MXG2),ISMLQ       
C$omp threadprivate(/MAXC/)
      LOGICAL         LRINT                                             
      COMMON /NLRCF / LRINT                                             
C$omp threadprivate(/NLRCF /)
      COMMON /SHLSPD/ CDA(MXGSH),CDB(MXGSH),CDC(MXGSH),CDD(MXGSH),      
     2                D02D(MXG2),D12D(MXG2),D22D(MXG2)                  
C$omp threadprivate(/SHLSPD/)
C                                                                       
      DIMENSION  WORK(6,11)                                             
C                                                                       
      PARAMETER (ZER=0.0D+00)                                           
      PARAMETER (PT5=0.5D+00)                                           
      PARAMETER (ONE=1.0D+00)                                           
      PARAMETER (PI4=0.7853981633974483D+00)                            
C                                                                       
         DO J= 1, 5                                                     
            FQD0(J)= ZER                                                
         ENDDO                                                          
         DO J= 1,10                                                     
            FQD1(1,J)= ZER                                              
            FQD1(2,J)= ZER                                              
         ENDDO                                                          
         DO J= 1,11                                                     
            FQD2(1,J)= ZER                                              
            FQD2(2,J)= ZER                                              
            FQD2(3,J)= ZER                                              
C                                                                       
            FQD3(1,J)= ZER                                              
            FQD3(2,J)= ZER                                              
            FQD3(3,J)= ZER                                              
            FQD3(4,J)= ZER                                              
         ENDDO                                                          
         DO J= 1, 7                                                     
            FQD4(1,J)= ZER                                              
            FQD4(2,J)= ZER                                              
            FQD4(3,J)= ZER                                              
            FQD4(4,J)= ZER                                              
            FQD4(5,J)= ZER                                              
         ENDDO                                                          
         DO J= 1, 4                                                     
            DO K= 1, 6                                                  
               FQD5(K,J)= ZER                                           
            ENDDO                                                       
         ENDDO                                                          
         DO K= 1, 7                                                     
            FQD6(K)= ZER                                                
         ENDDO                                                          
C                                                                       
      DO 300 I=1,NGANGB                                                 
         ISML= ISMLQ+ISMLP(I)                                           
         IF(ISML.GE.2) GO TO 300                                        
         X12= TX12(I)                                                   
         Y02= TY02(I)                                                   
         FQZ= D12D(I)                                                   
         X41= ONE/(X12+X34)                                             
         PQR= Y02-AQZ                                                   
         PQS= PQR*PQR                                                   
         RHO= X12*X34*X41                                               
         IF(LRINT) THEN                                                 
            EFR= EMU2/(EMU2+RHO)                                        
            RHO= RHO*EFR                                                
            FQZ= FQZ*SQRT(EFR)                                          
         ENDIF                                                          
         XVA=(PQS+QPS)*RHO                                              
         RHO= RHO+RHO                                                   
         N=6                                                            
         IF(XVA.LE.TMAX) THEN                                           
C                                                                       
C     FM(T) M=6 INTERPOLATION, GENERATING WASTED M=8,7 DATA             
C                                                                       
            M=5                                                         
            TV= XVA*RFINC(M)                                            
            IP= NINT(TV)                                                
            FX=    FGRID(4,IP,M) *TV                                    
            FX=(FX+FGRID(3,IP,M))*TV                                    
            FX=(FX+FGRID(2,IP,M))*TV                                    
            FX=(FX+FGRID(1,IP,M))*TV                                    
            FX= FX+FGRID(0,IP,M)                                        
            TV= XVA*RXINC                                               
            IP= NINT(TV)                                                
            ET=    XGRID(4,IP) *TV                                      
            ET=(ET+XGRID(3,IP))*TV                                      
            ET=(ET+XGRID(2,IP))*TV                                      
            ET=(ET+XGRID(1,IP))*TV                                      
            ET= ET+XGRID(0,IP)                                          
C                                                                       
            FQD(8)= FX                                                  
            T2= XVA+XVA                                                 
            DO M=8,1,-1                                                 
               FQD(M-1)=(T2*FQD(M)+ET)*RMR(M)                           
            END DO                                                      
C                                                                       
            FQF= FQZ*SQRT(X41)                                          
            DO 210 M=0,N                                                
               FQD(M)= FQD(M)*FQF                                       
  210       FQF= FQF*RHO                                                
         ELSE                                                           
            XIN= ONE/XVA                                                
            FQD(0)= FQZ*SQRT(PI4*XIN*X41)                               
            ROX= RHO*XIN                                                
            FQF= PT5*ROX                                                
            DO 220 M=1,N                                                
               FQD(M)= FQD(M-1)*FQF                                     
  220       FQF= FQF+ROX                                                
         ENDIF                                                          
C                                                                       
         PQT = PQR*PQS                                                  
         PQQ = PQS*PQS                                                  
         PQ5 = PQT*PQS                                                  
         XMD1= TX21(I)                                                  
         XMD2= XMD1*XMD1                                                
         Y01 = TY01(I)                                                  
         Y11 = Y01 *Y01                                                 
         DD02= D02D(I)                                                  
         WORK(1, 1)= XMD1*DD02                                          
         WORK(1, 2)= Y11 *DD02                                          
         WORK(1, 3)= Y11 *Y02                                           
         WORK(1, 4)= XMD1*Y02                                           
         WORK(1, 5)= XMD1*Y01                                           
         WORK(1, 6)= XMD1*Y01*DD02                                      
         WORK(1, 7)= XMD1*Y11                                           
         WORK(1, 8)= XMD2                                               
         WORK(1, 9)= XMD2*DD02                                          
         WORK(1,10)= XMD2*Y01                                           
         WORK(1,11)= XMD2*XMD1                                          
         DO J= 1,11                                                     
            WORK(2,J)= WORK(1,J)*PQR                                    
            WORK(3,J)= WORK(1,J)*PQS                                    
            WORK(4,J)= WORK(1,J)*PQT                                    
         ENDDO                                                          
         DO J= 5,11                                                     
            WORK(5,J)= WORK(1,J)*PQQ                                    
         ENDDO                                                          
         DO J= 8,11                                                     
            WORK(6,J)= WORK(1,J)*PQ5                                    
         ENDDO                                                          
C                                                                       
         DO J= 1, 5                                                     
            FQD0(J)= FQD0(J)+FQD(0)*WORK(1,J)                           
         ENDDO                                                          
         DO J= 1, 8                                                     
            FQD1(1,J)= FQD1(1,J)+FQD(1)*WORK(1,J)                       
            FQD1(2,J)= FQD1(2,J)+FQD(1)*WORK(2,J)                       
         ENDDO                                                          
            FQD1(1, 9)= FQD1(1, 9)+FQD(1)*WORK(1, 9)                    
            FQD1(1,10)= FQD1(1,10)+FQD(1)*WORK(1,10)                    
         DO J= 1,11                                                     
            FQD2(1,J)= FQD2(1,J)+FQD(2)*WORK(1,J)                       
            FQD2(2,J)= FQD2(2,J)+FQD(2)*WORK(2,J)                       
            FQD2(3,J)= FQD2(3,J)+FQD(2)*WORK(3,J)                       
C                                                                       
            FQD3(1,J)= FQD3(1,J)+FQD(3)*WORK(1,J)                       
            FQD3(2,J)= FQD3(2,J)+FQD(3)*WORK(2,J)                       
            FQD3(3,J)= FQD3(3,J)+FQD(3)*WORK(3,J)                       
            FQD3(4,J)= FQD3(4,J)+FQD(3)*WORK(4,J)                       
         ENDDO                                                          
         DO J= 1, 7                                                     
            FQD4(1,J)= FQD4(1,J)+FQD(4)*WORK(1,J+4)                     
            FQD4(2,J)= FQD4(2,J)+FQD(4)*WORK(2,J+4)                     
            FQD4(3,J)= FQD4(3,J)+FQD(4)*WORK(3,J+4)                     
            FQD4(4,J)= FQD4(4,J)+FQD(4)*WORK(4,J+4)                     
            FQD4(5,J)= FQD4(5,J)+FQD(4)*WORK(5,J+4)                     
         ENDDO                                                          
         DO J= 1, 4                                                     
            DO K= 1, 6                                                  
               FQD5(K,J)= FQD5(K,J)+FQD(5)*WORK(K,J+ 7)                 
            ENDDO                                                       
         ENDDO                                                          
         DO K= 1, 6                                                     
            FQD6(K)= FQD6(K)+FQD(6)*WORK(K,11)                          
         ENDDO                                                          
            FQD6(7)= FQD6(7)+FQD(6)*WORK(6,11)*PQR                      
  300 CONTINUE                                                          
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2R   *DECK INTJ20                                           
C>                                                                      
C>    @brief   DDDP case                                                
C>                                                                      
C>    @details integration of the DDDP case                             
C>                                                                      
      SUBROUTINE INTJ20                                                 
      USE lrcdft, ONLY: LCFLAG, EMU, EMU2, LRFILE                       
      use mx_limits, only: mxgtot,mxgsh,mxg2                            
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
C GENERATE JTYPE=20 INTEGRALS                                           
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
      COMMON /FQ08  / FQD(0:8),FQD0(5),FQD1(2,13),FQD2(3,16),FQD3(4,16),
     2                FQD4(5,16),FQD5(6,11),FQD6(7,7),FQD7(24),FQD8(9)  
C$omp threadprivate(/FQ08/)
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS,                      
     2                TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2),      
     3                D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2),      
     4                NGANGB                                            
C$omp threadprivate(/GEOMPQ/)
      COMMON /MAXC  / CMAX(MXGTOT),CMAXA(MXGSH),CMAXB(MXGSH),           
     2                CMAXC(MXGSH),CMAXD(MXGSH),ISMLP(MXG2),ISMLQ       
C$omp threadprivate(/MAXC/)
      LOGICAL         LRINT                                             
      COMMON /NLRCF / LRINT                                             
C$omp threadprivate(/NLRCF /)
      COMMON /SHLSPD/ CDA(MXGSH),CDB(MXGSH),CDC(MXGSH),CDD(MXGSH),      
     2                D02D(MXG2),D12D(MXG2),D22D(MXG2)                  
C$omp threadprivate(/SHLSPD/)
C                                                                       
      DIMENSION  WORK(7,11)                                             
C                                                                       
      PARAMETER (ZER=0.0D+00)                                           
      PARAMETER (PT5=0.5D+00)                                           
      PARAMETER (ONE=1.0D+00)                                           
      PARAMETER (PI4=0.7853981633974483D+00)                            
C                                                                       
         DO J= 1, 5                                                     
            FQD0(J)= ZER                                                
         ENDDO                                                          
         DO J= 1,10                                                     
            FQD1(1,J)= ZER                                              
            FQD1(2,J)= ZER                                              
         ENDDO                                                          
         DO J= 1,11                                                     
            FQD2(1,J)= ZER                                              
            FQD2(2,J)= ZER                                              
            FQD2(3,J)= ZER                                              
C                                                                       
            FQD3(1,J)= ZER                                              
            FQD3(2,J)= ZER                                              
            FQD3(3,J)= ZER                                              
            FQD3(4,J)= ZER                                              
C                                                                       
            FQD4(1,J)= ZER                                              
            FQD4(2,J)= ZER                                              
            FQD4(3,J)= ZER                                              
            FQD4(4,J)= ZER                                              
            FQD4(5,J)= ZER                                              
         ENDDO                                                          
         DO J= 1, 7                                                     
            DO K= 1, 6                                                  
               FQD5(K,J)= ZER                                           
            ENDDO                                                       
         ENDDO                                                          
         DO J= 1, 4                                                     
            DO K= 1, 7                                                  
               FQD6(K,J)= ZER                                           
            ENDDO                                                       
         ENDDO                                                          
         DO K= 1, 8                                                     
            FQD7(K)= ZER                                                
         ENDDO                                                          
C                                                                       
      DO 300 I=1,NGANGB                                                 
         ISML= ISMLQ+ISMLP(I)                                           
         IF(ISML.GE.2) GO TO 300                                        
         X12= TX12(I)                                                   
         Y02= TY02(I)                                                   
         FQZ= D12D(I)                                                   
         X41= ONE/(X12+X34)                                             
         PQR= Y02-AQZ                                                   
         PQS= PQR*PQR                                                   
         RHO= X12*X34*X41                                               
         IF(LRINT) THEN                                                 
            EFR= EMU2/(EMU2+RHO)                                        
            RHO= RHO*EFR                                                
            FQZ= FQZ*SQRT(EFR)                                          
         ENDIF                                                          
         XVA=(PQS+QPS)*RHO                                              
         RHO= RHO+RHO                                                   
         N=7                                                            
         IF(XVA.LE.TMAX) THEN                                           
C                                                                       
C     FM(T) M=7 INTERPOLATION, GENERATING WASTED M=8 DATA               
C                                                                       
            M=5                                                         
            TV= XVA*RFINC(M)                                            
            IP= NINT(TV)                                                
            FX=    FGRID(4,IP,M) *TV                                    
            FX=(FX+FGRID(3,IP,M))*TV                                    
            FX=(FX+FGRID(2,IP,M))*TV                                    
            FX=(FX+FGRID(1,IP,M))*TV                                    
            FX= FX+FGRID(0,IP,M)                                        
            TV= XVA*RXINC                                               
            IP= NINT(TV)                                                
            ET=    XGRID(4,IP) *TV                                      
            ET=(ET+XGRID(3,IP))*TV                                      
            ET=(ET+XGRID(2,IP))*TV                                      
            ET=(ET+XGRID(1,IP))*TV                                      
            ET= ET+XGRID(0,IP)                                          
C                                                                       
            FQD(8)= FX                                                  
            T2= XVA+XVA                                                 
            DO M=8,1,-1                                                 
               FQD(M-1)=(T2*FQD(M)+ET)*RMR(M)                           
            END DO                                                      
C                                                                       
            FQF= FQZ*SQRT(X41)                                          
            DO 210 M=0,N                                                
               FQD(M)= FQD(M)*FQF                                       
  210       FQF= FQF*RHO                                                
         ELSE                                                           
            XIN= ONE/XVA                                                
            FQD(0)= FQZ*SQRT(PI4*XIN*X41)                               
            ROX= RHO*XIN                                                
            FQF= PT5*ROX                                                
            DO 220 M=1,N                                                
               FQD(M)= FQD(M-1)*FQF                                     
  220       FQF= FQF+ROX                                                
         ENDIF                                                          
C                                                                       
         PQT = PQR*PQS                                                  
         PQQ = PQS*PQS                                                  
         PQ5 = PQT*PQS                                                  
         PQ6 = PQQ*PQS                                                  
         XMD1= TX21(I)                                                  
         XMD2= XMD1*XMD1                                                
         Y01 = TY01(I)                                                  
         Y11 = Y01 *Y01                                                 
         DD02= D02D(I)                                                  
         WORK(1, 1)= XMD1*DD02                                          
         WORK(1, 2)= Y11 *DD02                                          
         WORK(1, 3)= Y11 *Y02                                           
         WORK(1, 4)= XMD1*Y02                                           
         WORK(1, 5)= XMD1*Y01                                           
         WORK(1, 6)= XMD1*Y01*DD02                                      
         WORK(1, 7)= XMD1*Y11                                           
         WORK(1, 8)= XMD2                                               
         WORK(1, 9)= XMD2*DD02                                          
         WORK(1,10)= XMD2*Y01                                           
         WORK(1,11)= XMD2*XMD1                                          
         DO J= 1,11                                                     
            WORK(2,J)= WORK(1,J)*PQR                                    
            WORK(3,J)= WORK(1,J)*PQS                                    
            WORK(4,J)= WORK(1,J)*PQT                                    
            WORK(5,J)= WORK(1,J)*PQQ                                    
         ENDDO                                                          
         DO J= 5,11                                                     
            WORK(6,J)= WORK(1,J)*PQ5                                    
         ENDDO                                                          
         DO J= 8,11                                                     
            WORK(7,J)= WORK(1,J)*PQ6                                    
         ENDDO                                                          
C                                                                       
         DO J= 1, 5                                                     
            FQD0(J)= FQD0(J)+FQD(0)*WORK(1,J)                           
         ENDDO                                                          
         DO J= 1, 8                                                     
            FQD1(1,J)= FQD1(1,J)+FQD(1)*WORK(1,J)                       
            FQD1(2,J)= FQD1(2,J)+FQD(1)*WORK(2,J)                       
         ENDDO                                                          
            FQD1(1, 9)= FQD1(1, 9)+FQD(1)*WORK(1, 9)                    
            FQD1(1,10)= FQD1(1,10)+FQD(1)*WORK(1,10)                    
         DO J= 1,11                                                     
            FQD2(1,J)= FQD2(1,J)+FQD(2)*WORK(1,J)                       
            FQD2(2,J)= FQD2(2,J)+FQD(2)*WORK(2,J)                       
            FQD2(3,J)= FQD2(3,J)+FQD(2)*WORK(3,J)                       
C                                                                       
            FQD3(1,J)= FQD3(1,J)+FQD(3)*WORK(1,J)                       
            FQD3(2,J)= FQD3(2,J)+FQD(3)*WORK(2,J)                       
            FQD3(3,J)= FQD3(3,J)+FQD(3)*WORK(3,J)                       
            FQD3(4,J)= FQD3(4,J)+FQD(3)*WORK(4,J)                       
C                                                                       
            FQD4(1,J)= FQD4(1,J)+FQD(4)*WORK(1,J)                       
            FQD4(2,J)= FQD4(2,J)+FQD(4)*WORK(2,J)                       
            FQD4(3,J)= FQD4(3,J)+FQD(4)*WORK(3,J)                       
            FQD4(4,J)= FQD4(4,J)+FQD(4)*WORK(4,J)                       
            FQD4(5,J)= FQD4(5,J)+FQD(4)*WORK(5,J)                       
         ENDDO                                                          
         DO J= 1, 7                                                     
            DO K= 1, 6                                                  
               FQD5(K,J)= FQD5(K,J)+FQD(5)*WORK(K,J+ 4)                 
            ENDDO                                                       
         ENDDO                                                          
         DO J= 1, 4                                                     
            DO K= 1, 7                                                  
               FQD6(K,J)= FQD6(K,J)+FQD(6)*WORK(K,J+ 7)                 
            ENDDO                                                       
         ENDDO                                                          
         DO K= 1, 7                                                     
            FQD7(K)= FQD7(K)+FQD(7)*WORK(K,11)                          
         ENDDO                                                          
            FQD7(8)= FQD7(8)+FQD(7)*WORK(7,11)*PQR                      
  300 CONTINUE                                                          
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2R   *DECK INTJ21                                           
C>                                                                      
C>    @brief   DDDD case                                                
C>                                                                      
C>    @details integration of the DDDD case                             
C>                                                                      
      SUBROUTINE INTJ21                                                 
      USE lrcdft, ONLY: LCFLAG, EMU, EMU2, LRFILE                       
      use mx_limits, only: mxgtot,mxgsh,mxg2                            
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
C GENERATE JTYPE=21 INTEGRALS                                           
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
      COMMON /FQ08  / FQD(0:8),FQD0(5),FQD1(2,13),FQD2(3,16),FQD3(4,16),
     2                FQD4(5,16),FQD5(6,11),FQD6(7,7),FQD7(8,3),FQD8(9) 
C$omp threadprivate(/FQ08/)
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS,                      
     2                TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2),      
     3                D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2),      
     4                NGANGB                                            
C$omp threadprivate(/GEOMPQ/)
      COMMON /MAXC  / CMAX(MXGTOT),CMAXA(MXGSH),CMAXB(MXGSH),           
     2                CMAXC(MXGSH),CMAXD(MXGSH),ISMLP(MXG2),ISMLQ       
C$omp threadprivate(/MAXC/)
      LOGICAL         LRINT                                             
      COMMON /NLRCF / LRINT                                             
C$omp threadprivate(/NLRCF /)
      COMMON /SHLSPD/ CDA(MXGSH),CDB(MXGSH),CDC(MXGSH),CDD(MXGSH),      
     2                D02D(MXG2),D12D(MXG2),D22D(MXG2)                  
C$omp threadprivate(/SHLSPD/)
C                                                                       
      DIMENSION PQWK(2:8),WORK(8,16)                                    
C                                                                       
      PARAMETER (ZER=0.0D+00)                                           
      PARAMETER (PT5=0.5D+00)                                           
      PARAMETER (ONE=1.0D+00)                                           
      PARAMETER (PI4=0.7853981633974483D+00)                            
C                                                                       
         DO J= 1, 5                                                     
            FQD0(J)= ZER                                                
         ENDDO                                                          
         DO J= 1,13                                                     
            FQD1(1,J)= ZER                                              
            FQD1(2,J)= ZER                                              
         ENDDO                                                          
         DO J= 1,16                                                     
            FQD2(1,J)= ZER                                              
            FQD2(2,J)= ZER                                              
            FQD2(3,J)= ZER                                              
C                                                                       
            FQD3(1,J)= ZER                                              
            FQD3(2,J)= ZER                                              
            FQD3(3,J)= ZER                                              
            FQD3(4,J)= ZER                                              
C                                                                       
            FQD4(1,J)= ZER                                              
            FQD4(2,J)= ZER                                              
            FQD4(3,J)= ZER                                              
            FQD4(4,J)= ZER                                              
            FQD4(5,J)= ZER                                              
         ENDDO                                                          
         DO J= 1,11                                                     
            DO K= 1, 6                                                  
               FQD5(K,J)= ZER                                           
            ENDDO                                                       
         ENDDO                                                          
         DO J= 1, 7                                                     
            DO K= 1, 7                                                  
               FQD6(K,J)= ZER                                           
            ENDDO                                                       
         ENDDO                                                          
         DO J= 1, 3                                                     
            DO K= 1, 8                                                  
               FQD7(K,J)= ZER                                           
            ENDDO                                                       
         ENDDO                                                          
         DO K= 1, 9                                                     
            FQD8(K)= ZER                                                
         ENDDO                                                          
C                                                                       
      DO 300 I=1,NGANGB                                                 
         ISML= ISMLQ+ISMLP(I)                                           
         IF(ISML.GE.2) GO TO 300                                        
         X12= TX12(I)                                                   
         Y02= TY02(I)                                                   
         FQZ= D22D(I)                                                   
         X41= ONE/(X12+X34)                                             
         PQR= Y02-AQZ                                                   
         PQS= PQR*PQR                                                   
         RHO= X12*X34*X41                                               
         IF(LRINT) THEN                                                 
            EFR= EMU2/(EMU2+RHO)                                        
            RHO= RHO*EFR                                                
            FQZ= FQZ*SQRT(EFR)                                          
         ENDIF                                                          
         XVA=(PQS+QPS)*RHO                                              
         RHO= RHO+RHO                                                   
         N=8                                                            
         IF(XVA.LE.TMAX) THEN                                           
C                                                                       
C     FM(T) M=8 INTERPOLATION                                           
C                                                                       
            M=5                                                         
            TV= XVA*RFINC(M)                                            
            IP= NINT(TV)                                                
            FX=    FGRID(4,IP,M) *TV                                    
            FX=(FX+FGRID(3,IP,M))*TV                                    
            FX=(FX+FGRID(2,IP,M))*TV                                    
            FX=(FX+FGRID(1,IP,M))*TV                                    
            FX= FX+FGRID(0,IP,M)                                        
            TV= XVA*RXINC                                               
            IP= NINT(TV)                                                
            ET=    XGRID(4,IP) *TV                                      
            ET=(ET+XGRID(3,IP))*TV                                      
            ET=(ET+XGRID(2,IP))*TV                                      
            ET=(ET+XGRID(1,IP))*TV                                      
            ET= ET+XGRID(0,IP)                                          
C                                                                       
            FQD(8)= FX                                                  
            T2= XVA+XVA                                                 
            DO M=8,1,-1                                                 
               FQD(M-1)=(T2*FQD(M)+ET)*RMR(M)                           
            END DO                                                      
C                                                                       
            FQF= FQZ*SQRT(X41)                                          
            DO 210 M=0,N                                                
               FQD(M)= FQD(M)*FQF                                       
  210       FQF= FQF*RHO                                                
         ELSE                                                           
            XIN= ONE/XVA                                                
            FQD(0)= FQZ*SQRT(PI4*XIN*X41)                               
            ROX= RHO*XIN                                                
            FQF= PT5*ROX                                                
            DO 220 M=1,N                                                
               FQD(M)= FQD(M-1)*FQF                                     
  220       FQF= FQF+ROX                                                
         ENDIF                                                          
C                                                                       
         PQWK(2)= PQR                                                   
         PQWK(3)= PQS                                                   
         PQWK(4)= PQS*PQR                                               
         PQWK(5)= PQS*PQS                                               
         PQWK(6)= PQWK(4)*PQS                                           
         PQWK(7)= PQWK(5)*PQS                                           
         PQWK(8)= PQWK(6)*PQS                                           
         XMD1= TX21(I)                                                  
         XMD2= XMD1*XMD1                                                
         XMD3= XMD2*XMD1                                                
         Y01 = TY01(I)                                                  
         Y11 = Y01 *Y01                                                 
         Y12 = Y01 *Y02                                                 
         Y22 = Y02 *Y02                                                 
         WORK(1, 1)= XMD2                                               
         WORK(1, 2)= XMD1*Y11                                           
         WORK(1, 3)= XMD1*Y12                                           
         WORK(1, 4)= XMD1*Y22                                           
         WORK(1, 5)= Y11 *Y22                                           
         WORK(1, 6)= XMD2*Y01                                           
         WORK(1, 7)= XMD2*Y02                                           
         WORK(1, 8)= XMD1*Y11*Y02                                       
         WORK(1, 9)= XMD1*Y12*Y02                                       
         WORK(1,10)= XMD3                                               
         WORK(1,11)= XMD2*Y11                                           
         WORK(1,12)= XMD2*Y12                                           
         WORK(1,13)= XMD2*Y22                                           
         WORK(1,14)= XMD3*Y01                                           
         WORK(1,15)= XMD3*Y02                                           
         WORK(1,16)= XMD3*XMD1                                          
         DO J= 1, 5                                                     
            DO K= 2, 5                                                  
               WORK(K,J)= WORK(1,J)*PQWK(K)                             
            ENDDO                                                       
         ENDDO                                                          
         DO J= 6, 9                                                     
            DO K= 2, 6                                                  
               WORK(K,J)= WORK(1,J)*PQWK(K)                             
            ENDDO                                                       
         ENDDO                                                          
         DO J=10,13                                                     
            DO K= 2, 7                                                  
               WORK(K,J)= WORK(1,J)*PQWK(K)                             
            ENDDO                                                       
         ENDDO                                                          
         DO J=14,16                                                     
            DO K= 2, 8                                                  
               WORK(K,J)= WORK(1,J)*PQWK(K)                             
            ENDDO                                                       
         ENDDO                                                          
C                                                                       
         DO J= 1, 5                                                     
            FQD0(J)= FQD0(J)+FQD(0)*WORK(1,J)                           
         ENDDO                                                          
         DO J= 1,13                                                     
            FQD1(1,J)= FQD1(1,J)+FQD(1)*WORK(1,J)                       
            FQD1(2,J)= FQD1(2,J)+FQD(1)*WORK(2,J)                       
         ENDDO                                                          
         DO J= 1,16                                                     
            FQD2(1,J)= FQD2(1,J)+FQD(2)*WORK(1,J)                       
            FQD2(2,J)= FQD2(2,J)+FQD(2)*WORK(2,J)                       
            FQD2(3,J)= FQD2(3,J)+FQD(2)*WORK(3,J)                       
C                                                                       
            FQD3(1,J)= FQD3(1,J)+FQD(3)*WORK(1,J)                       
            FQD3(2,J)= FQD3(2,J)+FQD(3)*WORK(2,J)                       
            FQD3(3,J)= FQD3(3,J)+FQD(3)*WORK(3,J)                       
            FQD3(4,J)= FQD3(4,J)+FQD(3)*WORK(4,J)                       
C                                                                       
            FQD4(1,J)= FQD4(1,J)+FQD(4)*WORK(1,J)                       
            FQD4(2,J)= FQD4(2,J)+FQD(4)*WORK(2,J)                       
            FQD4(3,J)= FQD4(3,J)+FQD(4)*WORK(3,J)                       
            FQD4(4,J)= FQD4(4,J)+FQD(4)*WORK(4,J)                       
            FQD4(5,J)= FQD4(5,J)+FQD(4)*WORK(5,J)                       
         ENDDO                                                          
         DO J= 1,11                                                     
            DO K= 1, 6                                                  
               FQD5(K,J)= FQD5(K,J)+FQD(5)*WORK(K,J+ 5)                 
            ENDDO                                                       
         ENDDO                                                          
         DO J= 1, 7                                                     
            DO K= 1, 7                                                  
               FQD6(K,J)= FQD6(K,J)+FQD(6)*WORK(K,J+ 9)                 
            ENDDO                                                       
         ENDDO                                                          
         DO J= 1, 3                                                     
            DO K= 1, 8                                                  
               FQD7(K,J)= FQD7(K,J)+FQD(7)*WORK(K,J+13)                 
            ENDDO                                                       
         ENDDO                                                          
         DO K= 1, 8                                                     
            FQD8(K)= FQD8(K)+FQD(8)*WORK(K,16)                          
         ENDDO                                                          
            FQD8(9)= FQD8(9)+FQD(8)*WORK(8,16)*PQR                      
  300 CONTINUE                                                          
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2R   *DECK INTK07                                           
C>                                                                      
C>    @brief   DSSS case                                                
C>                                                                      
C>    @details integration of the DSSS case                             
C>                                                                      
      SUBROUTINE INTK07(IKL)                                            
      use mx_limits, only: mxgsh,mxg2                                   
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
C GENERATE JTYPE= 7 INTEGRALS                                           
C                                                                       
C                                                                       
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS,                      
     2                TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2),      
     3                D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2),      
     4                NGANGB                                            
C$omp threadprivate(/GEOMPQ/)
      COMMON /JMSGYH/ SQ(4)                                             
C$omp threadprivate(/JMSGYH/)
      COMMON /FQ08  / FQD(9),FQD0(5),FQD1(26),FQD2(48),FQD3(64),        
     2                FQD4(80),FQD5(66),FQD6(49),FQD7(24),FQD8( 9)      
C$omp threadprivate(/FQ08/)
      COMMON /KI2 / ACY,ACY2,AQX,AQX2,AQXY,Y03,Y04                      
C$omp threadprivate(/KI2/)
      COMMON /KI4 / RD0(25),RD1(120),RD2(336),RD3(520),RD4(630),        
     *              RD5(504),RD6(336),RD7(144),RD8(45)                  
C$omp threadprivate(/KI4/)
C                                                                       
      PARAMETER (ZER=0.0D+00)                                           
C                                                                       
      IF(IKL.EQ.0) THEN                                                 
         RD0(1)= ZER                                                    
         RD0(2)= ZER                                                    
         RD1(1)= ZER                                                    
         RD1(2)= ZER                                                    
         RD1(3)= ZER                                                    
         RD2(1)= ZER                                                    
         RD2(2)= ZER                                                    
         RD2(3)= ZER                                                    
         RD2(4)= ZER                                                    
         RD2(5)= ZER                                                    
         RD2(6)= ZER                                                    
C                                                                       
         RETURN                                                         
      ENDIF                                                             
C                                                                       
      XMD2= X43 *0.5D+00                                                
      XMD1= XMD2*SQ(1)                                                  
      XMDT= XMD2*XMD1                                                   
      XMD2=-XMD1*Y03                                                    
C                                                                       
      RD0(1)= RD0(1)+FQD0(1)*XMD1                                       
      RD0(2)= RD0(2)+FQD0(1)*Y03 *Y03*SQ(1)                             
C                                                                       
      RD1(1)= RD1(1)-FQD1(1)*AQX *XMD2                                  
      RD1(2)= RD1(2)-FQD1(1)*ACY *XMD2                                  
      RD1(3)= RD1(3)+FQD1(2)     *XMD2                                  
C                                                                       
      RD2(1)= RD2(1)+(FQD2(1)*AQX2-FQD1(1))*XMDT                        
      RD2(2)= RD2(2)+(FQD2(1)*ACY2-FQD1(1))*XMDT                        
      RD2(3)= RD2(3)+(FQD2(3)     -FQD1(1))*XMDT                        
      RD2(4)= RD2(4)+ FQD2(1)*AQXY         *XMDT                        
      RD2(5)= RD2(5)- FQD2(2)*AQX          *XMDT                        
      RD2(6)= RD2(6)- FQD2(2)*ACY          *XMDT                        
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2R   *DECK INTK08                                           
C>                                                                      
C>    @brief   DPSS case                                                
C>                                                                      
C>    @details integration of the DPSS case                             
C>                                                                      
      SUBROUTINE INTK08(IKL)                                            
      use mx_limits, only: mxgsh,mxg2                                   
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
C GENERATE JTYPE= 8 INTEGRALS                                           
C                                                                       
C                                                                       
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS,                      
     2                TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2),      
     3                D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2),      
     4                NGANGB                                            
C$omp threadprivate(/GEOMPQ/)
      COMMON /JMSGYH/ SQ(4)                                             
C$omp threadprivate(/JMSGYH/)
      COMMON /FQ08  / FQD(9),FQD0(5),FQD1(26),FQD2(48),FQD3(64),        
     2                FQD4(80),FQD5(66),FQD6(49),FQD7(24),FQD8( 9)      
C$omp threadprivate(/FQ08/)
      COMMON /KI2 / ACY,ACY2,AQX,AQX2,AQXY,Y03,Y04                      
C$omp threadprivate(/KI2/)
      COMMON /KI4 / RD0(25),RD1(3,40),RD2(6,56),RD3(520),RD4(630),      
     *              RD5(504),RD6(336),RD7(144),RD8(45)                  
C$omp threadprivate(/KI4/)
C                                                                       
      DIMENSION  WORK(4),FWK(6)                                         
C                                                                       
      PARAMETER (ZER=0.0D+00)                                           
      PARAMETER (F03=3.0D+00)                                           
C                                                                       
      IF(IKL.EQ.0) THEN                                                 
         RD0(1)= ZER                                                    
         RD0(2)= ZER                                                    
         RD0(3)= ZER                                                    
         RD0(4)= ZER                                                    
         RD0(5)= ZER                                                    
         DO I= 1, 4                                                     
            RD1(1,I)= ZER                                               
            RD1(2,I)= ZER                                               
            RD1(3,I)= ZER                                               
         ENDDO                                                          
         DO I= 1, 3                                                     
            RD2(1,I)= ZER                                               
            RD2(2,I)= ZER                                               
            RD2(3,I)= ZER                                               
            RD2(4,I)= ZER                                               
            RD2(5,I)= ZER                                               
            RD2(6,I)= ZER                                               
         ENDDO                                                          
         DO J= 1,10                                                     
            RD3(J)= ZER                                                 
         ENDDO                                                          
C                                                                       
         RETURN                                                         
      ENDIF                                                             
C                                                                       
      XMD2= X43 *0.5D+00                                                
      XMD1= XMD2*SQ(1)                                                  
      XMD3= XMD2*SQ(2)                                                  
      XMD4= XMD2*XMD2                                                   
      XMD5= XMD4*SQ(2)                                                  
      XMDT= XMD4*XMD3                                                   
C                                                                       
      XMDTY=-XMDT*ACY                                                   
      XMDTX=-XMDT*AQX                                                   
      XMDTXY=XMDT*AQXY                                                  
C                                                                       
      Y33 = Y03 *Y03                                                    
      Y34 =-Y03 *Y04                                                    
      Y334= Y33 *Y04                                                    
      WORK(1)=-XMD1*Y03                                                 
      WORK(2)= XMD5                                                     
      WORK(3)= XMD3*Y33                                                 
      WORK(4)= XMD3*Y34                                                 
C                                                                       
      RD0(1)= RD0(1)+FQD0(1)*XMD1                                       
      RD0(2)= RD0(2)+FQD0(1)*Y33 *SQ(1)                                 
      RD0(3)= RD0(3)-FQD0(1)*XMD3*Y03                                   
      RD0(4)= RD0(4)+FQD0(1)*XMD3*Y04                                   
      RD0(5)= RD0(5)+FQD0(1)*Y334*SQ(2)                                 
C                                                                       
         FWK(1)=-FQD1(1)*AQX                                            
         FWK(2)=-FQD1(1)*ACY                                            
         FWK(3)= FQD1(2)                                                
      DO I= 1, 4                                                        
         RD1(1,I)= RD1(1,I)+FWK(1)*WORK(I)                              
         RD1(2,I)= RD1(2,I)+FWK(2)*WORK(I)                              
         RD1(3,I)= RD1(3,I)+FWK(3)*WORK(I)                              
      ENDDO                                                             
C                                                                       
      WORK(1)= XMD4*SQ(1)                                               
      WORK(2)=-XMD5*Y03                                                 
      WORK(3)= XMD5*Y04                                                 
         FWK(1)= FQD2(1)*AQX2-FQD1(1)                                   
         FWK(2)= FQD2(1)*ACY2-FQD1(1)                                   
         FWK(3)= FQD2(3)     -FQD1(1)                                   
         FWK(4)= FQD2(1)*AQXY                                           
         FWK(5)=-FQD2(2)*AQX                                            
         FWK(6)=-FQD2(2)*ACY                                            
      DO 210 I=1,3                                                      
         DO 210 J=1,6                                                   
            RD2(J,I)= RD2(J,I)+FWK(J)*WORK(I)                           
  210 CONTINUE                                                          
C                                                                       
      RD3( 1)= RD3( 1)+(FQD3(1)*AQX2-FQD2(1)*F03)*XMDTX                 
      RD3( 2)= RD3( 2)+(FQD3(1)*AQX2-FQD2(1)    )*XMDTY                 
      RD3( 3)= RD3( 3)+(FQD3(2)*AQX2-FQD2(2)    )*XMDT                  
      RD3( 4)= RD3( 4)+(FQD3(1)*ACY2-FQD2(1)    )*XMDTX                 
      RD3( 5)= RD3( 5)+ FQD3(2)                  *XMDTXY                
      RD3( 6)= RD3( 6)+(FQD3(3)     -FQD2(1)    )*XMDTX                 
      RD3( 7)= RD3( 7)+(FQD3(1)*ACY2-FQD2(1)*F03)*XMDTY                 
      RD3( 8)= RD3( 8)+(FQD3(2)*ACY2-FQD2(2)    )*XMDT                  
      RD3( 9)= RD3( 9)+(FQD3(3)     -FQD2(1)    )*XMDTY                 
      RD3(10)= RD3(10)+(FQD3(4)     -FQD2(2)*F03)*XMDT                  
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2R   *DECK INTK09                                           
C>                                                                      
C>    @brief   DSPS case                                                
C>                                                                      
C>    @details integration of the DSPS case                             
C>                                                                      
      SUBROUTINE INTK09(IKL)                                            
      use mx_limits, only: mxgsh,mxg2                                   
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
C GENERATE JTYPE= 9 INTEGRALS                                           
C                                                                       
C                                                                       
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS,                      
     2                TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2),      
     3                D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2),      
     4                NGANGB                                            
C$omp threadprivate(/GEOMPQ/)
      COMMON /JMSGYH/ SQ(4)                                             
C$omp threadprivate(/JMSGYH/)
      COMMON /FQ08  / FQD(9),FQD0(5),FQD1(2,13),FQD2(3,16),FQD3(64),    
     2                FQD4(80),FQD5(66),FQD6(49),FQD7(24),FQD8( 9)      
C$omp threadprivate(/FQ08/)
      COMMON /KI2 / ACY,ACY2,AQX,AQX2,AQXY,Y03,Y04                      
C$omp threadprivate(/KI2/)
      COMMON /KI4 / RD0(25),RD1(3,40),RD2(6,56),RD3(520),RD4(630),      
     *              RD5(504),RD6(336),RD7(144),RD8(45)                  
C$omp threadprivate(/KI4/)
C                                                                       
      DIMENSION  WORK(3)                                                
C                                                                       
      PARAMETER (ZER=0.0D+00)                                           
      PARAMETER (F03=3.0D+00)                                           
C                                                                       
      IF(IKL.EQ.0) THEN                                                 
         RD0(1)= ZER                                                    
         RD0(2)= ZER                                                    
         RD0(3)= ZER                                                    
         RD0(4)= ZER                                                    
         DO I= 1, 4                                                     
            RD1(1,I)= ZER                                               
            RD1(2,I)= ZER                                               
            RD1(3,I)= ZER                                               
         ENDDO                                                          
         DO I= 1, 3                                                     
            RD2(1,I)= ZER                                               
            RD2(2,I)= ZER                                               
            RD2(3,I)= ZER                                               
            RD2(4,I)= ZER                                               
            RD2(5,I)= ZER                                               
            RD2(6,I)= ZER                                               
         ENDDO                                                          
         DO J= 1,10                                                     
            RD3(J)= ZER                                                 
         ENDDO                                                          
C                                                                       
         RETURN                                                         
      ENDIF                                                             
C                                                                       
      XMD2= X43 *0.5D+00                                                
      XMD3= XMD2*SQ(1)                                                  
      XMDT= XMD3*XMD2                                                   
C                                                                       
      XMDTY=-XMDT*ACY                                                   
      XMDTX=-XMDT*AQX                                                   
      XMDTXY=XMDT*AQXY                                                  
C                                                                       
      WORK(1)= XMD3                                                     
      WORK(2)= Y03 *Y03*SQ(1)                                           
      WORK(3)=-XMD3*Y03                                                 
C                                                                       
      RD0(1)= RD0(1)+FQD0(1)*WORK(1)                                    
      RD0(2)= RD0(2)+FQD0(1)*WORK(2)                                    
      RD0(3)= RD0(3)+FQD0(2)*WORK(1)                                    
      RD0(4)= RD0(4)+FQD0(2)*WORK(2)                                    
C                                                                       
      DO I= 1, 2                                                        
         RD1(1,I)= RD1(1,I)-FQD1(1,I)*AQX *WORK(3)                      
         RD1(2,I)= RD1(2,I)-FQD1(1,I)*ACY *WORK(3)                      
         RD1(3,I)= RD1(3,I)+FQD1(2,I)     *WORK(3)                      
      ENDDO                                                             
      DO I= 3, 4                                                        
         RD1(1,I)= RD1(1,I)-FQD1(1,3)*AQX *WORK(I- 2)                   
         RD1(2,I)= RD1(2,I)-FQD1(1,3)*ACY *WORK(I- 2)                   
         RD1(3,I)= RD1(3,I)+FQD1(2,3)     *WORK(I- 2)                   
      ENDDO                                                             
C                                                                       
      WORK(1)= XMDT                                                     
      WORK(2)= XMDT                                                     
      DO I= 1, 3                                                        
         RD2(1,I)= RD2(1,I)+(FQD2(1,I)*AQX2-FQD1(1,I))*WORK(I)          
         RD2(2,I)= RD2(2,I)+(FQD2(1,I)*ACY2-FQD1(1,I))*WORK(I)          
         RD2(3,I)= RD2(3,I)+(FQD2(3,I)     -FQD1(1,I))*WORK(I)          
         RD2(4,I)= RD2(4,I)+ FQD2(1,I)*AQXY           *WORK(I)          
         RD2(5,I)= RD2(5,I)- FQD2(2,I)*AQX            *WORK(I)          
         RD2(6,I)= RD2(6,I)- FQD2(2,I)*ACY            *WORK(I)          
      ENDDO                                                             
C                                                                       
      RD3( 1)= RD3( 1)+(FQD3(1)*AQX2-FQD2(1,3)*F03)*XMDTX               
      RD3( 2)= RD3( 2)+(FQD3(1)*AQX2-FQD2(1,3)    )*XMDTY               
      RD3( 3)= RD3( 3)+(FQD3(2)*AQX2-FQD2(2,3)    )*XMDT                
      RD3( 4)= RD3( 4)+(FQD3(1)*ACY2-FQD2(1,3)    )*XMDTX               
      RD3( 5)= RD3( 5)+ FQD3(2)                    *XMDTXY              
      RD3( 6)= RD3( 6)+(FQD3(3)     -FQD2(1,3)    )*XMDTX               
      RD3( 7)= RD3( 7)+(FQD3(1)*ACY2-FQD2(1,3)*F03)*XMDTY               
      RD3( 8)= RD3( 8)+(FQD3(2)*ACY2-FQD2(2,3)    )*XMDT                
      RD3( 9)= RD3( 9)+(FQD3(3)     -FQD2(1,3)    )*XMDTY               
      RD3(10)= RD3(10)+(FQD3(4)     -FQD2(2,3)*F03)*XMDT                
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2R   *DECK INTK10                                           
C>                                                                      
C>    @brief   DDSS case                                                
C>                                                                      
C>    @details integration of the DDSS case                             
C>                                                                      
      SUBROUTINE INTK10(IKL)                                            
      use mx_limits, only: mxgsh,mxg2                                   
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
C GENERATE JTYPE=10 INTEGRALS                                           
C                                                                       
C                                                                       
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS,                      
     2                TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2),      
     3                D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2),      
     4                NGANGB                                            
C$omp threadprivate(/GEOMPQ/)
      COMMON /JMSGYH/ SQ(4)                                             
C$omp threadprivate(/JMSGYH/)
      COMMON /FQ08  / FQD(9),FQD0(5),FQD1(26),FQD2(48),FQD3(64),        
     2                FQD4(80),FQD5(66),FQD6(49),FQD7(24),FQD8( 9)      
C$omp threadprivate(/FQ08/)
      COMMON /KI2 / ACY,ACY2,AQX,AQX2,AQXY,Y03,Y04                      
C$omp threadprivate(/KI2/)
      COMMON /KI4 / RD0(25),RD1(3,40),RD2(6,56),RD3(10,52),RD4(630),    
     *              RD5(504),RD6(336),RD7(144),RD8(45)                  
C$omp threadprivate(/KI4/)
C                                                                       
      DIMENSION  WORK(4),FWK(10)                                        
C                                                                       
      PARAMETER (ZER=0.0D+00)                                           
      PARAMETER (F03=3.0D+00)                                           
      PARAMETER (F06=6.0D+00)                                           
C                                                                       
      IF(IKL.EQ.0) THEN                                                 
         RD0(1)= ZER                                                    
         RD0(2)= ZER                                                    
         RD0(3)= ZER                                                    
         RD0(4)= ZER                                                    
         RD0(5)= ZER                                                    
         DO I= 1, 4                                                     
            RD1(1,I)= ZER                                               
            RD1(2,I)= ZER                                               
            RD1(3,I)= ZER                                               
C                                                                       
            RD2(1,I)= ZER                                               
            RD2(2,I)= ZER                                               
            RD2(3,I)= ZER                                               
            RD2(4,I)= ZER                                               
            RD2(5,I)= ZER                                               
            RD2(6,I)= ZER                                               
         ENDDO                                                          
         DO I= 1, 2                                                     
            DO J= 1,10                                                  
               RD3(J,I)= ZER                                            
            ENDDO                                                       
         ENDDO                                                          
         DO J= 1,15                                                     
            RD4(J)= ZER                                                 
         ENDDO                                                          
C                                                                       
         RETURN                                                         
      ENDIF                                                             
C                                                                       
      XMD2= X43 *0.5D+00                                                
      XMD3= XMD2*SQ(3)                                                  
      XMD4= XMD3*XMD2                                                   
      XMD6= XMD4*XMD2                                                   
      XMDT= XMD6*XMD2                                                   
      XMD2= XMD3                                                        
C                                                                       
      XMDTY=-XMDT*ACY                                                   
      XMDTX=-XMDT*AQX                                                   
      XMDTXY=XMDT*AQXY                                                  
C                                                                       
      Y33 = Y03 *Y03                                                    
      Y34 =-Y03 *Y04                                                    
      Y44 = Y04 *Y04                                                    
      Y334= Y33 *Y04                                                    
      Y344= Y34 *Y04                                                    
      WORK(1)=-XMD4*Y03                                                 
      WORK(2)= XMD4*Y04                                                 
      WORK(3)= XMD2*Y334                                                
      WORK(4)= XMD2*Y344                                                
C                                                                       
      RD0(1)= RD0(1)+FQD0(1)*XMD4                                       
      RD0(2)= RD0(2)+FQD0(1)*XMD2*Y33                                   
      RD0(3)= RD0(3)+FQD0(1)*XMD2*Y34                                   
      RD0(4)= RD0(4)+FQD0(1)*XMD2*Y44                                   
      RD0(5)= RD0(5)+FQD0(1)*Y33 *Y44*SQ(3)                             
C                                                                       
         FWK(1)=-FQD1(1)*AQX                                            
         FWK(2)=-FQD1(1)*ACY                                            
         FWK(3)= FQD1(2)                                                
      DO I= 1, 4                                                        
         RD1(1,I)= RD1(1,I)+FWK(1)*WORK(I)                              
         RD1(2,I)= RD1(2,I)+FWK(2)*WORK(I)                              
         RD1(3,I)= RD1(3,I)+FWK(3)*WORK(I)                              
      ENDDO                                                             
C                                                                       
      WORK(1)= XMD6                                                     
      WORK(2)= XMD4*Y33                                                 
      WORK(3)= XMD4*Y34                                                 
      WORK(4)= XMD4*Y44                                                 
         FWK(1)= FQD2(1)*AQX2-FQD1(1)                                   
         FWK(2)= FQD2(1)*ACY2-FQD1(1)                                   
         FWK(3)= FQD2(3)     -FQD1(1)                                   
         FWK(4)= FQD2(1)*AQXY                                           
         FWK(5)=-FQD2(2)*AQX                                            
         FWK(6)=-FQD2(2)*ACY                                            
      DO 210 I= 1, 4                                                    
         DO 210 J= 1, 6                                                 
            RD2(J,I)= RD2(J,I)+FWK(J)*WORK(I)                           
  210 CONTINUE                                                          
C                                                                       
      WORK(1)=-XMD6*Y03                                                 
      WORK(2)= XMD6*Y04                                                 
         FWK( 1)=-(FQD3(1)*AQX2-FQD2(1)*F03)*AQX                        
         FWK( 2)=-(FQD3(1)*AQX2-FQD2(1)    )*ACY                        
         FWK( 3)=  FQD3(2)*AQX2-FQD2(2)                                 
         FWK( 4)=-(FQD3(1)*ACY2-FQD2(1)    )*AQX                        
         FWK( 5)=  FQD3(2)                  *AQXY                       
         FWK( 6)=-(FQD3(3)     -FQD2(1)    )*AQX                        
         FWK( 7)=-(FQD3(1)*ACY2-FQD2(1)*F03)*ACY                        
         FWK( 8)=  FQD3(2)*ACY2-FQD2(2)                                 
         FWK( 9)=-(FQD3(3)     -FQD2(1)    )*ACY                        
         FWK(10)=  FQD3(4)     -FQD2(2)*F03                             
      DO 310 I= 1, 2                                                    
         DO 310 J= 1,10                                                 
            RD3(J,I)= RD3(J,I)+FWK(J)*WORK(I)                           
  310 CONTINUE                                                          
C                                                                       
      AQX4= AQX2*AQX2                                                   
      ACY4= ACY2*ACY2                                                   
      X2Y2= AQX2*ACY2                                                   
      Q2C2= AQX2+ACY2                                                   
      RD4( 1)= RD4( 1)+(FQD4(1)*AQX4-FQD3(1)*F06*AQX2                   
     *                              +FQD2(1)*F03         )*XMDT         
      RD4( 2)= RD4( 2)+(FQD4(1)*AQX2-FQD3(1)*F03         )*XMDTXY       
      RD4( 3)= RD4( 3)+(FQD4(2)*AQX2-FQD3(2)*F03         )*XMDTX        
      RD4( 4)= RD4( 4)+(FQD4(1)*X2Y2-FQD3(1)*Q2C2+FQD2(1))*XMDT         
      RD4( 5)= RD4( 5)+(FQD4(2)*AQX2-FQD3(2)             )*XMDTY        
      RD4( 6)= RD4( 6)+(FQD4(3)*AQX2-FQD3(1)*AQX2-FQD3(3)               
     *                              +FQD2(1)             )*XMDT         
      RD4( 7)= RD4( 7)+(FQD4(1)*ACY2-FQD3(1)*F03         )*XMDTXY       
      RD4( 8)= RD4( 8)+(FQD4(2)*ACY2-FQD3(2)             )*XMDTX        
      RD4( 9)= RD4( 9)+(FQD4(3)     -FQD3(1)             )*XMDTXY       
      RD4(10)= RD4(10)+(FQD4(4)     -FQD3(2)*F03         )*XMDTX        
      RD4(11)= RD4(11)+(FQD4(1)*ACY4-FQD3(1)*F06*ACY2                   
     *                              +FQD2(1)*F03         )*XMDT         
      RD4(12)= RD4(12)+(FQD4(2)*ACY2-FQD3(2)*F03         )*XMDTY        
      RD4(13)= RD4(13)+(FQD4(3)*ACY2-FQD3(1)*ACY2-FQD3(3)               
     *                              +FQD2(1)             )*XMDT         
      RD4(14)= RD4(14)+(FQD4(4)     -FQD3(2)*F03         )*XMDTY        
      RD4(15)= RD4(15)+(FQD4(5)  -FQD3(3)*F06+FQD2(1)*F03)*XMDT         
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2R   *DECK INTK11                                           
C>                                                                      
C>    @brief   DPPS case                                                
C>                                                                      
C>    @details integration of the DPPS case                             
C>                                                                      
      SUBROUTINE INTK11(IKL)                                            
      use mx_limits, only: mxgsh,mxg2                                   
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
C GENERATE JTYPE=11 INTEGRALS                                           
C                                                                       
C                                                                       
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS,                      
     2                TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2),      
     3                D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2),      
     4                NGANGB                                            
C$omp threadprivate(/GEOMPQ/)
      COMMON /JMSGYH/ SQ(4)                                             
C$omp threadprivate(/JMSGYH/)
      COMMON /FQ08  / FQD(9),FQD0(5),FQD1(2,13),FQD2(3,16),FQD3(4,16),  
     2                FQD4(80),FQD5(66),FQD6(49),FQD7(24),FQD8( 9)      
C$omp threadprivate(/FQ08/)
      COMMON /KI2 / ACY,ACY2,AQX,AQX2,AQXY,Y03,Y04                      
C$omp threadprivate(/KI2/)
      COMMON /KI4 / RD0(5,5),RD1(3,40),RD2(6,56),RD3(10,52),RD4(630),   
     *              RD5(504),RD6(336),RD7(144),RD8(45)                  
C$omp threadprivate(/KI4/)
C                                                                       
      DIMENSION  WORK(12),FWK(10,3)                                     
C                                                                       
      PARAMETER (ZER=0.0D+00)                                           
      PARAMETER (F03=3.0D+00)                                           
      PARAMETER (F06=6.0D+00)                                           
C                                                                       
      IF(IKL.EQ.0) THEN                                                 
         DO I= 1, 2                                                     
            RD0(1,I)= ZER                                               
            RD0(2,I)= ZER                                               
            RD0(3,I)= ZER                                               
            RD0(4,I)= ZER                                               
            RD0(5,I)= ZER                                               
         ENDDO                                                          
         DO I= 1,13                                                     
            RD1(1,I)= ZER                                               
            RD1(2,I)= ZER                                               
            RD1(3,I)= ZER                                               
         ENDDO                                                          
         DO I= 1,10                                                     
            RD2(1,I)= ZER                                               
            RD2(2,I)= ZER                                               
            RD2(3,I)= ZER                                               
            RD2(4,I)= ZER                                               
            RD2(5,I)= ZER                                               
            RD2(6,I)= ZER                                               
         ENDDO                                                          
         DO I= 1, 5                                                     
            DO J=1,10                                                   
               RD3(J,I)= ZER                                            
            ENDDO                                                       
         ENDDO                                                          
         DO J= 1,15                                                     
            RD4(J)= ZER                                                 
         ENDDO                                                          
C                                                                       
         RETURN                                                         
      ENDIF                                                             
C                                                                       
      XMD2= X43 *0.5D+00                                                
      XMD1= XMD2*SQ(1)                                                  
      XMD3= XMD2*SQ(2)                                                  
      XMD4= XMD2*XMD2                                                   
      XMD5= XMD4*SQ(2)                                                  
      XMDT= XMD4*XMD3                                                   
C                                                                       
      XMDTY=-XMDT*ACY                                                   
      XMDTX=-XMDT*AQX                                                   
      XMDTXY=XMDT*AQXY                                                  
C                                                                       
      Y33 = Y03 *Y03                                                    
      Y34 =-Y03 *Y04*SQ(2)                                              
      Y334= Y33 *Y04*SQ(2)                                              
C                                                                       
      WORK( 1)= XMD1                                                    
      WORK( 2)= Y33 *SQ(1)                                              
      WORK( 3)=-XMD3*Y03                                                
      WORK( 4)= XMD3*Y04                                                
      WORK( 5)= Y334                                                    
      WORK( 6)=-XMD1*Y03                                                
      WORK( 7)= XMD5                                                    
      WORK( 8)= XMD3*Y33                                                
      WORK( 9)= XMD2*Y34                                                
      WORK(10)= XMD4*SQ(1)                                              
      WORK(11)=-XMD5*Y03                                                
      WORK(12)= XMD5*Y04                                                
C                                                                       
      DO 010 I= 1, 2                                                    
         DO 010 J= 1, 5                                                 
            RD0(J,I)= RD0(J,I)+FQD0(I)*WORK(J)                          
  010 CONTINUE                                                          
C                                                                       
      DO I= 1, 3                                                        
         FWK(1,I)=-FQD1(1,I)*AQX                                        
         FWK(2,I)=-FQD1(1,I)*ACY                                        
         FWK(3,I)= FQD1(2,I)                                            
      ENDDO                                                             
      DO I= 1, 4                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,1)*WORK(I+ 5)                         
         RD1(2,I)= RD1(2,I)+FWK(2,1)*WORK(I+ 5)                         
         RD1(3,I)= RD1(3,I)+FWK(3,1)*WORK(I+ 5)                         
      ENDDO                                                             
      DO I= 5, 8                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,2)*WORK(I+ 1)                         
         RD1(2,I)= RD1(2,I)+FWK(2,2)*WORK(I+ 1)                         
         RD1(3,I)= RD1(3,I)+FWK(3,2)*WORK(I+ 1)                         
      ENDDO                                                             
      DO I= 9,13                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,3)*WORK(I- 8)                         
         RD1(2,I)= RD1(2,I)+FWK(2,3)*WORK(I- 8)                         
         RD1(3,I)= RD1(3,I)+FWK(3,3)*WORK(I- 8)                         
      ENDDO                                                             
C                                                                       
      DO I= 1, 3                                                        
         FWK(1,I)= FQD2(1,I)*AQX2-FQD1(1,I)                             
         FWK(2,I)= FQD2(1,I)*ACY2-FQD1(1,I)                             
         FWK(3,I)= FQD2(3,I)     -FQD1(1,I)                             
         FWK(4,I)= FQD2(1,I)*AQXY                                       
         FWK(5,I)=-FQD2(2,I)*AQX                                        
         FWK(6,I)=-FQD2(2,I)*ACY                                        
      ENDDO                                                             
      DO I= 1, 3                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,1)*WORK(I+ 9)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I= 4, 6                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,2)*WORK(I+ 6)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I= 7,10                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,3)*WORK(I- 1)                      
         ENDDO                                                          
      ENDDO                                                             
C                                                                       
      DO I= 1, 2                                                        
         RD3( 1,I)= RD3( 1,I)+(FQD3(1,I)*AQX2-FQD2(1,I)*F03)*XMDTX      
         RD3( 2,I)= RD3( 2,I)+(FQD3(1,I)*AQX2-FQD2(1,I)    )*XMDTY      
         RD3( 3,I)= RD3( 3,I)+(FQD3(2,I)*AQX2-FQD2(2,I)    )*XMDT       
         RD3( 4,I)= RD3( 4,I)+(FQD3(1,I)*ACY2-FQD2(1,I)    )*XMDTX      
         RD3( 5,I)= RD3( 5,I)+ FQD3(2,I)                    *XMDTXY     
         RD3( 6,I)= RD3( 6,I)+(FQD3(3,I)     -FQD2(1,I)    )*XMDTX      
         RD3( 7,I)= RD3( 7,I)+(FQD3(1,I)*ACY2-FQD2(1,I)*F03)*XMDTY      
         RD3( 8,I)= RD3( 8,I)+(FQD3(2,I)*ACY2-FQD2(2,I)    )*XMDT       
         RD3( 9,I)= RD3( 9,I)+(FQD3(3,I)     -FQD2(1,I)    )*XMDTY      
         RD3(10,I)= RD3(10,I)+(FQD3(4,I)     -FQD2(2,I)*F03)*XMDT       
      ENDDO                                                             
         FWK( 1,3)=-(FQD3(1,3)*AQX2-FQD2(1,3)*F03)*AQX                  
         FWK( 2,3)=-(FQD3(1,3)*AQX2-FQD2(1,3)    )*ACY                  
         FWK( 3,3)=  FQD3(2,3)*AQX2-FQD2(2,3)                           
         FWK( 4,3)=-(FQD3(1,3)*ACY2-FQD2(1,3)    )*AQX                  
         FWK( 5,3)=  FQD3(2,3)                    *AQXY                 
         FWK( 6,3)=-(FQD3(3,3)     -FQD2(1,3)    )*AQX                  
         FWK( 7,3)=-(FQD3(1,3)*ACY2-FQD2(1,3)*F03)*ACY                  
         FWK( 8,3)=  FQD3(2,3)*ACY2-FQD2(2,3)                           
         FWK( 9,3)=-(FQD3(3,3)     -FQD2(1,3)    )*ACY                  
         FWK(10,3)=  FQD3(4,3)     -FQD2(2,3)*F03                       
      DO 320 I= 3, 5                                                    
         DO 320 J= 1,10                                                 
            RD3(J,I)= RD3(J,I)+FWK(J,3)*WORK(I+ 7)                      
  320 CONTINUE                                                          
C                                                                       
      AQX4= AQX2*AQX2                                                   
      ACY4= ACY2*ACY2                                                   
      X2Y2= AQX2*ACY2                                                   
      Q2C2= AQX2+ACY2                                                   
      RD4( 1)= RD4( 1)+(FQD4(1)*AQX4-FQD3(1,3)*F06*AQX2                 
     *                              +FQD2(1,3)*F03           )*XMDT     
      RD4( 2)= RD4( 2)+(FQD4(1)*AQX2-FQD3(1,3)*F03           )*XMDTXY   
      RD4( 3)= RD4( 3)+(FQD4(2)*AQX2-FQD3(2,3)*F03           )*XMDTX    
      RD4( 4)= RD4( 4)+(FQD4(1)*X2Y2-FQD3(1,3)*Q2C2+FQD2(1,3))*XMDT     
      RD4( 5)= RD4( 5)+(FQD4(2)*AQX2-FQD3(2,3)               )*XMDTY    
      RD4( 6)= RD4( 6)+(FQD4(3)*AQX2-FQD3(1,3)*AQX2-FQD3(3,3)           
     *                              +FQD2(1,3)               )*XMDT     
      RD4( 7)= RD4( 7)+(FQD4(1)*ACY2-FQD3(1,3)*F03           )*XMDTXY   
      RD4( 8)= RD4( 8)+(FQD4(2)*ACY2-FQD3(2,3)               )*XMDTX    
      RD4( 9)= RD4( 9)+(FQD4(3)     -FQD3(1,3)               )*XMDTXY   
      RD4(10)= RD4(10)+(FQD4(4)     -FQD3(2,3)*F03           )*XMDTX    
      RD4(11)= RD4(11)+(FQD4(1)*ACY4-FQD3(1,3)*F06*ACY2                 
     *                              +FQD2(1,3)*F03           )*XMDT     
      RD4(12)= RD4(12)+(FQD4(2)*ACY2-FQD3(2,3)*F03           )*XMDTY    
      RD4(13)= RD4(13)+(FQD4(3)*ACY2-FQD3(1,3)*ACY2-FQD3(3,3)           
     *                              +FQD2(1,3)               )*XMDT     
      RD4(14)= RD4(14)+(FQD4(4)     -FQD3(2,3)*F03           )*XMDTY    
      RD4(15)= RD4(15)+(FQD4(5)  -FQD3(3,3)*F06+FQD2(1,3)*F03)*XMDT     
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2R   *DECK INTK12                                           
C>                                                                      
C>    @brief   DSDS case                                                
C>                                                                      
C>    @details integration of the DSDS case                             
C>                                                                      
      SUBROUTINE INTK12(IKL)                                            
      use mx_limits, only: mxgsh,mxg2                                   
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
C GENERATE JTYPE=12 INTEGRALS                                           
C                                                                       
C                                                                       
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS,                      
     2                TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2),      
     3                D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2),      
     4                NGANGB                                            
C$omp threadprivate(/GEOMPQ/)
      COMMON /JMSGYH/ SQ(4)                                             
C$omp threadprivate(/JMSGYH/)
      COMMON /FQ08  / FQD(9),FQD0(5),FQD1(2,13),FQD2(3,16),FQD3(4,16),  
     2                FQD4(80),FQD5(66),FQD6(49),FQD7(24),FQD8( 9)      
C$omp threadprivate(/FQ08/)
      COMMON /KI2 / ACY,ACY2,AQX,AQX2,AQXY,Y03,Y04                      
C$omp threadprivate(/KI2/)
      COMMON /KI4 / RD0(25),RD1(3,40),RD2(6,56),RD3(10,52),RD4(630),    
     *              RD5(504),RD6(336),RD7(144),RD8(45)                  
C$omp threadprivate(/KI4/)
C                                                                       
      DIMENSION  WORK(5),FWK(10)                                        
C                                                                       
      PARAMETER (ZER=0.0D+00)                                           
      PARAMETER (F03=3.0D+00)                                           
      PARAMETER (F06=6.0D+00)                                           
C                                                                       
      IF(IKL.EQ.0) THEN                                                 
         RD0(1)= ZER                                                    
         RD0(2)= ZER                                                    
         RD0(3)= ZER                                                    
         RD0(4)= ZER                                                    
         DO I= 1, 4                                                     
            RD1(1,I)= ZER                                               
            RD1(2,I)= ZER                                               
            RD1(3,I)= ZER                                               
         ENDDO                                                          
         DO I= 1, 5                                                     
            RD2(1,I)= ZER                                               
            RD2(2,I)= ZER                                               
            RD2(3,I)= ZER                                               
            RD2(4,I)= ZER                                               
            RD2(5,I)= ZER                                               
            RD2(6,I)= ZER                                               
         ENDDO                                                          
         DO I= 1, 2                                                     
            DO J= 1,10                                                  
               RD3(J,I)= ZER                                            
            ENDDO                                                       
         ENDDO                                                          
         DO J= 1,15                                                     
            RD4(J)= ZER                                                 
         ENDDO                                                          
C                                                                       
         RETURN                                                         
      ENDIF                                                             
C                                                                       
      XMD2= X43 *0.5D+00                                                
      XMD1= XMD2*SQ(1)                                                  
      XMD3=-XMD1*Y03                                                    
      XMDT= XMD1*XMD2                                                   
C                                                                       
      XMD3Y=-XMD3*ACY                                                   
      XMD3X=-XMD3*AQX                                                   
      XMD3XY=XMD3*AQXY                                                  
      XMDTY=-XMDT*ACY                                                   
      XMDTX=-XMDT*AQX                                                   
      XMDTXY=XMDT*AQXY                                                  
C                                                                       
      WORK(1)= XMD1                                                     
      WORK(2)= Y03 *Y03*SQ(1)                                           
      WORK(3)= XMDT                                                     
      WORK(4)= XMDT                                                     
      WORK(5)= XMD3                                                     
C                                                                       
      RD0(1)= RD0(1)+FQD0(1)*WORK(1)                                    
      RD0(2)= RD0(2)+FQD0(1)*WORK(2)                                    
      RD0(3)= RD0(3)+FQD0(2)*WORK(1)                                    
      RD0(4)= RD0(4)+FQD0(2)*WORK(2)                                    
C                                                                       
      DO I= 1, 2                                                        
         RD1(1,I)= RD1(1,I)+FQD1(1,I)*XMD3X                             
         RD1(2,I)= RD1(2,I)+FQD1(1,I)*XMD3Y                             
         RD1(3,I)= RD1(3,I)+FQD1(2,I)*XMD3                              
      ENDDO                                                             
         FWK(1)=-FQD1(1,3)*AQX                                          
         FWK(2)=-FQD1(1,3)*ACY                                          
         FWK(3)= FQD1(2,3)                                              
      DO I= 3, 4                                                        
         RD1(1,I)= RD1(1,I)+FWK(1)*WORK(I- 2)                           
         RD1(2,I)= RD1(2,I)+FWK(2)*WORK(I- 2)                           
         RD1(3,I)= RD1(3,I)+FWK(3)*WORK(I- 2)                           
      ENDDO                                                             
C                                                                       
      DO I= 1, 3                                                        
         RD2(1,I)= RD2(1,I)+(FQD2(1,I)*AQX2-FQD1(1,I))*WORK(I+ 2)       
         RD2(2,I)= RD2(2,I)+(FQD2(1,I)*ACY2-FQD1(1,I))*WORK(I+ 2)       
         RD2(3,I)= RD2(3,I)+(FQD2(3,I)     -FQD1(1,I))*WORK(I+ 2)       
         RD2(4,I)= RD2(4,I)+ FQD2(1,I)*AQXY           *WORK(I+ 2)       
         RD2(5,I)= RD2(5,I)- FQD2(2,I)*AQX            *WORK(I+ 2)       
         RD2(6,I)= RD2(6,I)- FQD2(2,I)*ACY            *WORK(I+ 2)       
      ENDDO                                                             
         FWK(1)= FQD2(1,4)*AQX2-FQD1(1,4)                               
         FWK(2)= FQD2(1,4)*ACY2-FQD1(1,4)                               
         FWK(3)= FQD2(3,4)     -FQD1(1,4)                               
         FWK(4)= FQD2(1,4)*AQXY                                         
         FWK(5)=-FQD2(2,4)*AQX                                          
         FWK(6)=-FQD2(2,4)*ACY                                          
      DO 220 I= 4, 5                                                    
         DO 220 J= 1, 6                                                 
            RD2(J,I)= RD2(J,I)+FWK(J)*WORK(I-3)                         
  220 CONTINUE                                                          
C                                                                       
      RD3( 1,1)= RD3( 1,1)+(FQD3(1,1)*AQX2-FQD2(1,3)*F03)*XMDTX         
      RD3( 2,1)= RD3( 2,1)+(FQD3(1,1)*AQX2-FQD2(1,3)    )*XMDTY         
      RD3( 3,1)= RD3( 3,1)+(FQD3(2,1)*AQX2-FQD2(2,3)    )*XMDT          
      RD3( 4,1)= RD3( 4,1)+(FQD3(1,1)*ACY2-FQD2(1,3)    )*XMDTX         
      RD3( 5,1)= RD3( 5,1)+ FQD3(2,1)                    *XMDTXY        
      RD3( 6,1)= RD3( 6,1)+(FQD3(3,1)     -FQD2(1,3)    )*XMDTX         
      RD3( 7,1)= RD3( 7,1)+(FQD3(1,1)*ACY2-FQD2(1,3)*F03)*XMDTY         
      RD3( 8,1)= RD3( 8,1)+(FQD3(2,1)*ACY2-FQD2(2,3)    )*XMDT          
      RD3( 9,1)= RD3( 9,1)+(FQD3(3,1)     -FQD2(1,3)    )*XMDTY         
      RD3(10,1)= RD3(10,1)+(FQD3(4,1)     -FQD2(2,3)*F03)*XMDT          
C                                                                       
      RD3( 1,2)= RD3( 1,2)+(FQD3(1,2)*AQX2-FQD2(1,4)*F03)*XMD3X         
      RD3( 2,2)= RD3( 2,2)+(FQD3(1,2)*AQX2-FQD2(1,4)    )*XMD3Y         
      RD3( 3,2)= RD3( 3,2)+(FQD3(2,2)*AQX2-FQD2(2,4)    )*XMD3          
      RD3( 4,2)= RD3( 4,2)+(FQD3(1,2)*ACY2-FQD2(1,4)    )*XMD3X         
      RD3( 5,2)= RD3( 5,2)+ FQD3(2,2)                    *XMD3XY        
      RD3( 6,2)= RD3( 6,2)+(FQD3(3,2)     -FQD2(1,4)    )*XMD3X         
      RD3( 7,2)= RD3( 7,2)+(FQD3(1,2)*ACY2-FQD2(1,4)*F03)*XMD3Y         
      RD3( 8,2)= RD3( 8,2)+(FQD3(2,2)*ACY2-FQD2(2,4)    )*XMD3          
      RD3( 9,2)= RD3( 9,2)+(FQD3(3,2)     -FQD2(1,4)    )*XMD3Y         
      RD3(10,2)= RD3(10,2)+(FQD3(4,2)     -FQD2(2,4)*F03)*XMD3          
C                                                                       
      AQX4= AQX2*AQX2                                                   
      ACY4= ACY2*ACY2                                                   
      X2Y2= AQX2*ACY2                                                   
      Q2C2= AQX2+ACY2                                                   
      RD4( 1)= RD4( 1)+(FQD4(1)*AQX4-FQD3(1,2)*F06*AQX2                 
     *                              +FQD2(1,4)*F03           )*XMDT     
      RD4( 2)= RD4( 2)+(FQD4(1)*AQX2-FQD3(1,2)*F03           )*XMDTXY   
      RD4( 3)= RD4( 3)+(FQD4(2)*AQX2-FQD3(2,2)*F03           )*XMDTX    
      RD4( 4)= RD4( 4)+(FQD4(1)*X2Y2-FQD3(1,2)*Q2C2+FQD2(1,4))*XMDT     
      RD4( 5)= RD4( 5)+(FQD4(2)*AQX2-FQD3(2,2)               )*XMDTY    
      RD4( 6)= RD4( 6)+(FQD4(3)*AQX2-FQD3(1,2)*AQX2-FQD3(3,2)           
     *                              +FQD2(1,4)               )*XMDT     
      RD4( 7)= RD4( 7)+(FQD4(1)*ACY2-FQD3(1,2)*F03           )*XMDTXY   
      RD4( 8)= RD4( 8)+(FQD4(2)*ACY2-FQD3(2,2)               )*XMDTX    
      RD4( 9)= RD4( 9)+(FQD4(3)     -FQD3(1,2)               )*XMDTXY   
      RD4(10)= RD4(10)+(FQD4(4)     -FQD3(2,2)*F03           )*XMDTX    
      RD4(11)= RD4(11)+(FQD4(1)*ACY4-FQD3(1,2)*F06*ACY2                 
     *                              +FQD2(1,4)*F03           )*XMDT     
      RD4(12)= RD4(12)+(FQD4(2)*ACY2-FQD3(2,2)*F03           )*XMDTY    
      RD4(13)= RD4(13)+(FQD4(3)*ACY2-FQD3(1,2)*ACY2-FQD3(3,2)           
     *                              +FQD2(1,4)               )*XMDT     
      RD4(14)= RD4(14)+(FQD4(4)     -FQD3(2,2)*F03           )*XMDTY    
      RD4(15)= RD4(15)+(FQD4(5)  -FQD3(3,2)*F06+FQD2(1,4)*F03)*XMDT     
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2R   *DECK INTK13                                           
C>                                                                      
C>    @brief   DSPP case                                                
C>                                                                      
C>    @details integration of the DSPP case                             
C>                                                                      
      SUBROUTINE INTK13(IKL)                                            
      use mx_limits, only: mxgsh,mxg2                                   
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
C GENERATE JTYPE=13 INTEGRALS                                           
C                                                                       
C                                                                       
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS,                      
     2                TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2),      
     3                D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2),      
     4                NGANGB                                            
C$omp threadprivate(/GEOMPQ/)
      COMMON /JMSGYH/ SQ(4)                                             
C$omp threadprivate(/JMSGYH/)
      COMMON /FQ08  / FQD(9),FQD0(5),FQD1(2,13),FQD2(3,16),FQD3(4,16),  
     2                FQD4(80),FQD5(66),FQD6(49),FQD7(24),FQD8( 9)      
C$omp threadprivate(/FQ08/)
      COMMON /KI2 / ACY,ACY2,AQX,AQX2,AQXY,Y03,Y04                      
C$omp threadprivate(/KI2/)
      COMMON /KI4 / RD0(5,5),RD1(3,40),RD2(6,56),RD3(10,52),RD4(630),   
     *              RD5(504),RD6(336),RD7(144),RD8(45)                  
C$omp threadprivate(/KI4/)
C                                                                       
      DIMENSION  WORK(2),FWK(6,4)                                       
C                                                                       
      PARAMETER (ZER=0.0D+00)                                           
      PARAMETER (F03=3.0D+00)                                           
      PARAMETER (F06=6.0D+00)                                           
C                                                                       
      IF(IKL.EQ.0) THEN                                                 
         DO I= 1, 2                                                     
            RD0(1,I)= ZER                                               
            RD0(2,I)= ZER                                               
            RD0(3,I)= ZER                                               
            RD0(4,I)= ZER                                               
            RD0(5,I)= ZER                                               
         ENDDO                                                          
         DO I= 1,13                                                     
            RD1(1,I)= ZER                                               
            RD1(2,I)= ZER                                               
            RD1(3,I)= ZER                                               
         ENDDO                                                          
         DO I= 1,11                                                     
            RD2(1,I)= ZER                                               
            RD2(2,I)= ZER                                               
            RD2(3,I)= ZER                                               
            RD2(4,I)= ZER                                               
            RD2(5,I)= ZER                                               
            RD2(6,I)= ZER                                               
         ENDDO                                                          
         DO I= 1, 5                                                     
            DO J= 1,10                                                  
               RD3(J,I)= ZER                                            
            ENDDO                                                       
         ENDDO                                                          
         DO J= 1,15                                                     
            RD4(J)= ZER                                                 
         ENDDO                                                          
C                                                                       
         RETURN                                                         
      ENDIF                                                             
C                                                                       
      XMD2= X43 *0.5D+00                                                
      XMD1= XMD2*SQ(1)                                                  
      XMD3=-XMD1*Y03                                                    
      XMDT= XMD1*XMD2                                                   
C                                                                       
      XMD3Y=-XMD3*ACY                                                   
      XMD3X=-XMD3*AQX                                                   
      XMD3XY=XMD3*AQXY                                                  
      XMDTY=-XMDT*ACY                                                   
      XMDTX=-XMDT*AQX                                                   
      XMDTXY=XMDT*AQXY                                                  
C                                                                       
      WORK(1)= XMD1                                                     
      WORK(2)= Y03 *Y03*SQ(1)                                           
      DO 010 I= 1, 2                                                    
         DO 010 J= 1, 5                                                 
            RD0(J,I)= RD0(J,I)+FQD0(J)*WORK(I)                          
  010 CONTINUE                                                          
C                                                                       
      DO I= 1, 5                                                        
         RD1(1,I)= RD1(1,I)+FQD1(1,I)*XMD3X                             
         RD1(2,I)= RD1(2,I)+FQD1(1,I)*XMD3Y                             
         RD1(3,I)= RD1(3,I)+FQD1(2,I)*XMD3                              
      ENDDO                                                             
      DO I= 1, 3                                                        
         FWK(1,I)=-FQD1(1,I+5)*AQX                                      
         FWK(2,I)=-FQD1(1,I+5)*ACY                                      
         FWK(3,I)= FQD1(2,I+5)                                          
      ENDDO                                                             
      FQD11 = FQD1(1,5)*RAB +FQD1(1,8)                                  
      FQD12 = FQD1(2,5)*RAB +FQD1(2,8)                                  
         FWK(1,4)=-FQD11*AQX                                            
         FWK(2,4)=-FQD11*ACY                                            
         FWK(3,4)= FQD12                                                
      DO I= 6, 7                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,1)*WORK(I- 5)                         
         RD1(2,I)= RD1(2,I)+FWK(2,1)*WORK(I- 5)                         
         RD1(3,I)= RD1(3,I)+FWK(3,1)*WORK(I- 5)                         
      ENDDO                                                             
      DO I= 8, 9                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,2)*WORK(I- 7)                         
         RD1(2,I)= RD1(2,I)+FWK(2,2)*WORK(I- 7)                         
         RD1(3,I)= RD1(3,I)+FWK(3,2)*WORK(I- 7)                         
      ENDDO                                                             
      DO I=10,11                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,3)*WORK(I- 9)                         
         RD1(2,I)= RD1(2,I)+FWK(2,3)*WORK(I- 9)                         
         RD1(3,I)= RD1(3,I)+FWK(3,3)*WORK(I- 9)                         
      ENDDO                                                             
      DO I=12,13                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,4)*WORK(I-11)                         
         RD1(2,I)= RD1(2,I)+FWK(2,4)*WORK(I-11)                         
         RD1(3,I)= RD1(3,I)+FWK(3,4)*WORK(I-11)                         
      ENDDO                                                             
C                                                                       
      DO I= 1, 5                                                        
         RD2(1,I)= RD2(1,I)+(FQD2(1,I)*AQX2-FQD1(1,I))*XMDT             
         RD2(2,I)= RD2(2,I)+(FQD2(1,I)*ACY2-FQD1(1,I))*XMDT             
         RD2(3,I)= RD2(3,I)+(FQD2(3,I)     -FQD1(1,I))*XMDT             
         RD2(4,I)= RD2(4,I)+ FQD2(1,I)                *XMDTXY           
         RD2(5,I)= RD2(5,I)+ FQD2(2,I)                *XMDTX            
         RD2(6,I)= RD2(6,I)+ FQD2(2,I)                *XMDTY            
      ENDDO                                                             
      DO I= 6, 8                                                        
         RD2(1,I)= RD2(1,I)+(FQD2(1,I)*AQX2-FQD1(1,I))*XMD3             
         RD2(2,I)= RD2(2,I)+(FQD2(1,I)*ACY2-FQD1(1,I))*XMD3             
         RD2(3,I)= RD2(3,I)+(FQD2(3,I)     -FQD1(1,I))*XMD3             
         RD2(4,I)= RD2(4,I)+ FQD2(1,I)                *XMD3XY           
         RD2(5,I)= RD2(5,I)+ FQD2(2,I)                *XMD3X            
         RD2(6,I)= RD2(6,I)+ FQD2(2,I)                *XMD3Y            
      ENDDO                                                             
      FQD2(1,5)= FQD2(1,5)*RAB +FQD2(1,8)                               
      FQD2(2,5)= FQD2(2,5)*RAB +FQD2(2,8)                               
      FQD2(3,5)= FQD2(3,5)*RAB +FQD2(3,8)                               
         RD2(1,9)= RD2(1,9)+(FQD2(1,5)*AQX2-FQD11)*XMD3                 
         RD2(2,9)= RD2(2,9)+(FQD2(1,5)*ACY2-FQD11)*XMD3                 
         RD2(3,9)= RD2(3,9)+(FQD2(3,5)     -FQD11)*XMD3                 
         RD2(4,9)= RD2(4,9)+ FQD2(1,5)            *XMD3XY               
         RD2(5,9)= RD2(5,9)+ FQD2(2,5)            *XMD3X                
         RD2(6,9)= RD2(6,9)+ FQD2(2,5)            *XMD3Y                
         FWK(1,1)= FQD2(1,9)*AQX2-FQD1(1,9)                             
         FWK(2,1)= FQD2(1,9)*ACY2-FQD1(1,9)                             
         FWK(3,1)= FQD2(3,9)     -FQD1(1,9)                             
         FWK(4,1)= FQD2(1,9)*AQXY                                       
         FWK(5,1)=-FQD2(2,9)*AQX                                        
         FWK(6,1)=-FQD2(2,9)*ACY                                        
      DO 240 I=10,11                                                    
         DO 240 J= 1, 6                                                 
            RD2(J,I)= RD2(J,I)+FWK(J,1)*WORK(I- 9)                      
  240 CONTINUE                                                          
C                                                                       
      FQD3(1,1)= FQD3(1,1)*RAB +FQD3(1,4)                               
      FQD3(2,1)= FQD3(2,1)*RAB +FQD3(2,4)                               
      FQD3(3,1)= FQD3(3,1)*RAB +FQD3(3,4)                               
      FQD3(4,1)= FQD3(4,1)*RAB +FQD3(4,4)                               
      DO I= 1, 4                                                        
         RD3( 1,I)= RD3( 1,I)+(FQD3(1,I)*AQX2-FQD2(1,I+4)*F03)*XMDTX    
         RD3( 2,I)= RD3( 2,I)+(FQD3(1,I)*AQX2-FQD2(1,I+4)    )*XMDTY    
         RD3( 3,I)= RD3( 3,I)+(FQD3(2,I)*AQX2-FQD2(2,I+4)    )*XMDT     
         RD3( 4,I)= RD3( 4,I)+(FQD3(1,I)*ACY2-FQD2(1,I+4)    )*XMDTX    
         RD3( 5,I)= RD3( 5,I)+ FQD3(2,I)                      *XMDTXY   
         RD3( 6,I)= RD3( 6,I)+(FQD3(3,I)     -FQD2(1,I+4)    )*XMDTX    
         RD3( 7,I)= RD3( 7,I)+(FQD3(1,I)*ACY2-FQD2(1,I+4)*F03)*XMDTY    
         RD3( 8,I)= RD3( 8,I)+(FQD3(2,I)*ACY2-FQD2(2,I+4)    )*XMDT     
         RD3( 9,I)= RD3( 9,I)+(FQD3(3,I)     -FQD2(1,I+4)    )*XMDTY    
         RD3(10,I)= RD3(10,I)+(FQD3(4,I)     -FQD2(2,I+4)*F03)*XMDT     
      ENDDO                                                             
         RD3( 1,5)= RD3( 1,5)+(FQD3(1,5)*AQX2-FQD2(1,9)*F03)*XMD3X      
         RD3( 2,5)= RD3( 2,5)+(FQD3(1,5)*AQX2-FQD2(1,9)    )*XMD3Y      
         RD3( 3,5)= RD3( 3,5)+(FQD3(2,5)*AQX2-FQD2(2,9)    )*XMD3       
         RD3( 4,5)= RD3( 4,5)+(FQD3(1,5)*ACY2-FQD2(1,9)    )*XMD3X      
         RD3( 5,5)= RD3( 5,5)+ FQD3(2,5)                    *XMD3XY     
         RD3( 6,5)= RD3( 6,5)+(FQD3(3,5)     -FQD2(1,9)    )*XMD3X      
         RD3( 7,5)= RD3( 7,5)+(FQD3(1,5)*ACY2-FQD2(1,9)*F03)*XMD3Y      
         RD3( 8,5)= RD3( 8,5)+(FQD3(2,5)*ACY2-FQD2(2,9)    )*XMD3       
         RD3( 9,5)= RD3( 9,5)+(FQD3(3,5)     -FQD2(1,9)    )*XMD3Y      
         RD3(10,5)= RD3(10,5)+(FQD3(4,5)     -FQD2(2,9)*F03)*XMD3       
C                                                                       
      AQX4= AQX2*AQX2                                                   
      ACY4= ACY2*ACY2                                                   
      X2Y2= AQX2*ACY2                                                   
      Q2C2= AQX2+ACY2                                                   
      RD4( 1)= RD4( 1)+(FQD4(1)*AQX4-FQD3(1,5)*F06*AQX2                 
     *                              +FQD2(1,9)*F03           )*XMDT     
      RD4( 2)= RD4( 2)+(FQD4(1)*AQX2-FQD3(1,5)*F03           )*XMDTXY   
      RD4( 3)= RD4( 3)+(FQD4(2)*AQX2-FQD3(2,5)*F03           )*XMDTX    
      RD4( 4)= RD4( 4)+(FQD4(1)*X2Y2-FQD3(1,5)*Q2C2+FQD2(1,9))*XMDT     
      RD4( 5)= RD4( 5)+(FQD4(2)*AQX2-FQD3(2,5)               )*XMDTY    
      RD4( 6)= RD4( 6)+(FQD4(3)*AQX2-FQD3(1,5)*AQX2                     
     *                              -FQD3(3,5)+FQD2(1,9)     )*XMDT     
      RD4( 7)= RD4( 7)+(FQD4(1)*ACY2-FQD3(1,5)*F03           )*XMDTXY   
      RD4( 8)= RD4( 8)+(FQD4(2)*ACY2-FQD3(2,5)               )*XMDTX    
      RD4( 9)= RD4( 9)+(FQD4(3)     -FQD3(1,5)               )*XMDTXY   
      RD4(10)= RD4(10)+(FQD4(4)     -FQD3(2,5)*F03           )*XMDTX    
      RD4(11)= RD4(11)+(FQD4(1)*ACY4-FQD3(1,5)*F06*ACY2                 
     *                              +FQD2(1,9)*F03           )*XMDT     
      RD4(12)= RD4(12)+(FQD4(2)*ACY2-FQD3(2,5)*F03           )*XMDTY    
      RD4(13)= RD4(13)+(FQD4(3)*ACY2-FQD3(1,5)*ACY2-FQD3(3,5)           
     *                              +FQD2(1,9)               )*XMDT     
      RD4(14)= RD4(14)+(FQD4(4)     -FQD3(2,5)*F03           )*XMDTY    
      RD4(15)= RD4(15)+(FQD4(5)  -FQD3(3,5)*F06+FQD2(1,9)*F03)*XMDT     
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2R   *DECK INTK14                                           
C>                                                                      
C>    @brief   DDPS case                                                
C>                                                                      
C>    @details integration of the DDPS case                             
C>                                                                      
      SUBROUTINE INTK14(IKL)                                            
      use mx_limits, only: mxgsh,mxg2                                   
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
C GENERATE JTYPE=14 INTEGRALS                                           
C                                                                       
C                                                                       
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS,                      
     2                TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2),      
     3                D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2),      
     4                NGANGB                                            
C$omp threadprivate(/GEOMPQ/)
      COMMON /JMSGYH/ SQ(4)                                             
C$omp threadprivate(/JMSGYH/)
      COMMON /FQ08  / FQD(9),FQD0(5),FQD1(2,13),FQD2(3,16),FQD3(4,16),  
     2                FQD4(5,16),FQD5(66),FQD6(49),FQD7(24),FQD8( 9)    
C$omp threadprivate(/FQ08/)
      COMMON /KI2 / ACY,ACY2,AQX,AQX2,AQXY,Y03,Y04                      
C$omp threadprivate(/KI2/)
      COMMON /KI4 / RD0(5,5),RD1(3,40),RD2(6,56),RD3(10,52),RD4(15,42), 
     *              RD5(504),RD6(336),RD7(144),RD8(45)                  
C$omp threadprivate(/KI4/)
C                                                                       
      DIMENSION  WORK(15),FWK(15,3)                                     
C                                                                       
      PARAMETER (ZER=0.0D+00)                                           
      PARAMETER (F03=3.0D+00)                                           
      PARAMETER (F06=6.0D+00)                                           
      PARAMETER (F10=1.0D+01)                                           
      PARAMETER (F15=1.5D+01)                                           
C                                                                       
      IF(IKL.EQ.0) THEN                                                 
         DO I= 1, 2                                                     
            RD0(1,I)= ZER                                               
            RD0(2,I)= ZER                                               
            RD0(3,I)= ZER                                               
            RD0(4,I)= ZER                                               
            RD0(5,I)= ZER                                               
         ENDDO                                                          
         DO I= 1,13                                                     
            RD1(1,I)= ZER                                               
            RD1(2,I)= ZER                                               
            RD1(3,I)= ZER                                               
         ENDDO                                                          
         DO I= 1,12                                                     
            RD2(1,I)= ZER                                               
            RD2(2,I)= ZER                                               
            RD2(3,I)= ZER                                               
            RD2(4,I)= ZER                                               
            RD2(5,I)= ZER                                               
            RD2(6,I)= ZER                                               
         ENDDO                                                          
         DO I= 1, 8                                                     
            DO J= 1,10                                                  
               RD3(J,I)= ZER                                            
            ENDDO                                                       
         ENDDO                                                          
         DO I= 1, 4                                                     
            DO J= 1,15                                                  
               RD4(J,I)= ZER                                            
            ENDDO                                                       
         ENDDO                                                          
         DO J= 1,21                                                     
            RD5(J)= ZER                                                 
         ENDDO                                                          
C                                                                       
         RETURN                                                         
      ENDIF                                                             
C                                                                       
      XMD2= X43 *0.5D+00                                                
      XMD3= XMD2*SQ(3)                                                  
      XMD4= XMD3*XMD2                                                   
      XMD6= XMD4*XMD2                                                   
      XMDT= XMD6*XMD2                                                   
      XMD2= XMD3                                                        
C                                                                       
      XMDTY=-XMDT*ACY                                                   
      XMDTX=-XMDT*AQX                                                   
      XMDTXY=XMDT*AQXY                                                  
C                                                                       
      Y33 = Y03 *Y03                                                    
      Y34 =-Y03 *Y04                                                    
      Y44 = Y04 *Y04                                                    
      Y334= Y33 *Y04                                                    
      Y344= Y34 *Y04                                                    
      WORK( 1)= XMD4                                                    
      WORK( 2)= XMD2*Y33                                                
      WORK( 3)= XMD2*Y34                                                
      WORK( 4)= XMD2*Y44                                                
      WORK( 5)= Y33 *Y44*SQ(3)                                          
C                                                                       
      WORK( 6)=-XMD4*Y03                                                
      WORK( 7)= XMD4*Y04                                                
      WORK( 8)= XMD2*Y334                                               
      WORK( 9)= XMD2*Y344                                               
C                                                                       
      WORK(10)= XMD6                                                    
      WORK(11)= XMD4*Y33                                                
      WORK(12)= XMD4*Y34                                                
      WORK(13)= XMD4*Y44                                                
C                                                                       
      WORK(14)=-XMD6*Y03                                                
      WORK(15)= XMD6*Y04                                                
C                                                                       
      DO 010 I= 1, 2                                                    
         DO 010 J= 1, 5                                                 
            RD0(J,I)= RD0(J,I)+FQD0(I)*WORK(J)                          
  010 CONTINUE                                                          
C                                                                       
      DO I= 1, 3                                                        
         FWK(1,I)=-FQD1(1,I)*AQX                                        
         FWK(2,I)=-FQD1(1,I)*ACY                                        
         FWK(3,I)= FQD1(2,I)                                            
      ENDDO                                                             
      DO I= 1, 4                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,1)*WORK(I+ 5)                         
         RD1(2,I)= RD1(2,I)+FWK(2,1)*WORK(I+ 5)                         
         RD1(3,I)= RD1(3,I)+FWK(3,1)*WORK(I+ 5)                         
      ENDDO                                                             
      DO I= 5, 8                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,2)*WORK(I+ 1)                         
         RD1(2,I)= RD1(2,I)+FWK(2,2)*WORK(I+ 1)                         
         RD1(3,I)= RD1(3,I)+FWK(3,2)*WORK(I+ 1)                         
      ENDDO                                                             
      DO I= 9,13                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,3)*WORK(I- 8)                         
         RD1(2,I)= RD1(2,I)+FWK(2,3)*WORK(I- 8)                         
         RD1(3,I)= RD1(3,I)+FWK(3,3)*WORK(I- 8)                         
      ENDDO                                                             
C                                                                       
      DO I= 1, 3                                                        
         FWK(1,I)= FQD2(1,I)*AQX2-FQD1(1,I)                             
         FWK(2,I)= FQD2(1,I)*ACY2-FQD1(1,I)                             
         FWK(3,I)= FQD2(3,I)     -FQD1(1,I)                             
         FWK(4,I)= FQD2(1,I)*AQXY                                       
         FWK(5,I)=-FQD2(2,I)*AQX                                        
         FWK(6,I)=-FQD2(2,I)*ACY                                        
      ENDDO                                                             
      DO I= 1, 4                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,1)*WORK(I+ 9)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I= 5, 8                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,2)*WORK(I+ 5)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I= 9,12                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,3)*WORK(I- 3)                      
         ENDDO                                                          
      ENDDO                                                             
C                                                                       
      DO I= 1, 3                                                        
         FWK( 1,I)=-(FQD3(1,I)*AQX2-FQD2(1,I)*F03)*AQX                  
         FWK( 2,I)=-(FQD3(1,I)*AQX2-FQD2(1,I)    )*ACY                  
         FWK( 3,I)=  FQD3(2,I)*AQX2-FQD2(2,I)                           
         FWK( 4,I)=-(FQD3(1,I)*ACY2-FQD2(1,I)    )*AQX                  
         FWK( 5,I)=  FQD3(2,I)                    *AQXY                 
         FWK( 6,I)=-(FQD3(3,I)     -FQD2(1,I)    )*AQX                  
         FWK( 7,I)=-(FQD3(1,I)*ACY2-FQD2(1,I)*F03)*ACY                  
         FWK( 8,I)=  FQD3(2,I)*ACY2-FQD2(2,I)                           
         FWK( 9,I)=-(FQD3(3,I)     -FQD2(1,I)    )*ACY                  
         FWK(10,I)=  FQD3(4,I)     -FQD2(2,I)*F03                       
      ENDDO                                                             
      DO I= 1, 2                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,1)*WORK(I+13)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I= 3, 4                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,2)*WORK(I+11)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I= 5, 8                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,3)*WORK(I+ 5)                      
         ENDDO                                                          
      ENDDO                                                             
C                                                                       
      AQX4= AQX2*AQX2                                                   
      ACY4= ACY2*ACY2                                                   
      X2Y2= AQX2*ACY2                                                   
      Q2C2= AQX2+ACY2                                                   
      DO I= 1, 2                                                        
         RD4( 1,I)= RD4( 1,I)+(FQD4(1,I)*AQX4-FQD3(1,I)*F06*AQX2        
     *                        +FQD2(1,I)*F03                 )*XMDT     
         RD4( 2,I)= RD4( 2,I)+(FQD4(1,I)*AQX2-FQD3(1,I)*F03  )*XMDTXY   
         RD4( 3,I)= RD4( 3,I)+(FQD4(2,I)*AQX2-FQD3(2,I)*F03  )*XMDTX    
         RD4( 4,I)= RD4( 4,I)+(FQD4(1,I)*X2Y2-FQD3(1,I)*Q2C2            
     *                        +FQD2(1,I)                     )*XMDT     
         RD4( 5,I)= RD4( 5,I)+(FQD4(2,I)*AQX2-FQD3(2,I)      )*XMDTY    
         RD4( 6,I)= RD4( 6,I)+(FQD4(3,I)*AQX2-FQD3(1,I)*AQX2            
     *                        -FQD3(3,I)+FQD2(1,I)           )*XMDT     
         RD4( 7,I)= RD4( 7,I)+(FQD4(1,I)*ACY2-FQD3(1,I)*F03  )*XMDTXY   
         RD4( 8,I)= RD4( 8,I)+(FQD4(2,I)*ACY2-FQD3(2,I)      )*XMDTX    
         RD4( 9,I)= RD4( 9,I)+(FQD4(3,I)     -FQD3(1,I)      )*XMDTXY   
         RD4(10,I)= RD4(10,I)+(FQD4(4,I)     -FQD3(2,I)*F03  )*XMDTX    
         RD4(11,I)= RD4(11,I)+(FQD4(1,I)*ACY4-FQD3(1,I)*F06*ACY2        
     *                        +FQD2(1,I)*F03                 )*XMDT     
         RD4(12,I)= RD4(12,I)+(FQD4(2,I)*ACY2-FQD3(2,I)*F03  )*XMDTY    
         RD4(13,I)= RD4(13,I)+(FQD4(3,I)*ACY2-FQD3(1,I)*ACY2            
     *                        -FQD3(3,I)+FQD2(1,I)           )*XMDT     
         RD4(14,I)= RD4(14,I)+(FQD4(4,I)     -FQD3(2,I)*F03  )*XMDTY    
         RD4(15,I)= RD4(15,I)+(FQD4(5,I)-FQD3(3,I)*F06                  
     *                        +FQD2(1,I)*F03                 )*XMDT     
      ENDDO                                                             
      FWK( 1,3)=  FQD4(1,3)*AQX4-FQD3(1,3)*F06*AQX2+FQD2(1,3)*F03       
      FWK( 2,3)= (FQD4(1,3)*AQX2-FQD3(1,3)*F03               )*AQXY     
      FWK( 3,3)=-(FQD4(2,3)*AQX2-FQD3(2,3)*F03               )*AQX      
      FWK( 4,3)=  FQD4(1,3)*X2Y2-FQD3(1,3)*Q2C2+FQD2(1,3)               
      FWK( 5,3)=-(FQD4(2,3)*AQX2-FQD3(2,3)                   )*ACY      
      FWK( 6,3)=  FQD4(3,3)*AQX2-FQD3(1,3)*AQX2-FQD3(3,3)+FQD2(1,3)     
      FWK( 7,3)= (FQD4(1,3)*ACY2-FQD3(1,3)*F03               )*AQXY     
      FWK( 8,3)=-(FQD4(2,3)*ACY2-FQD3(2,3)                   )*AQX      
      FWK( 9,3)= (FQD4(3,3)     -FQD3(1,3)                   )*AQXY     
      FWK(10,3)=-(FQD4(4,3)     -FQD3(2,3)*F03               )*AQX      
      FWK(11,3)=  FQD4(1,3)*ACY4-FQD3(1,3)*F06*ACY2+FQD2(1,3)*F03       
      FWK(12,3)=-(FQD4(2,3)*ACY2-FQD3(2,3)*F03               )*ACY      
      FWK(13,3)=  FQD4(3,3)*ACY2-FQD3(1,3)*ACY2-FQD3(3,3)+FQD2(1,3)     
      FWK(14,3)=-(FQD4(4,3)     -FQD3(2,3)*F03               )*ACY      
      FWK(15,3)=  FQD4(5,3)     -FQD3(3,3)*F06+FQD2(1,3)*F03            
      DO 420 I= 3, 4                                                    
         DO 420 J= 1,15                                                 
            RD4(J,I)= RD4(J,I)+FWK(J,3)*WORK(I+11)                      
  420 CONTINUE                                                          
C                                                                       
      RD5( 1)= RD5( 1)+(FQD5(1)*AQX4-FQD4(1,3)*F10*AQX2                 
     *                 +FQD3(1,3)*F15                        )*XMDTX    
      RD5( 2)= RD5( 2)+(FQD5(1)*AQX4-FQD4(1,3)*F06*AQX2                 
     *                 +FQD3(1,3)*F03                        )*XMDTY    
      RD5( 3)= RD5( 3)+(FQD5(2)*AQX4-FQD4(2,3)*F06*AQX2                 
     *                 +FQD3(2,3)*F03                        )*XMDT     
      RD5( 4)= RD5( 4)+(FQD5(1)*X2Y2-FQD4(1,3)*AQX2-FQD4(1,3)*F03*ACY2  
     *                 +FQD3(1,3)*F03                        )*XMDTX    
      RD5( 5)= RD5( 5)+(FQD5(2)*AQX2-FQD4(2,3)*F03           )*XMDTXY   
      RD5( 6)= RD5( 6)+(FQD5(3)*AQX2-FQD4(1,3)*AQX2-FQD4(3,3)*F03       
     *                 +FQD3(1,3)*F03                        )*XMDTX    
      RD5( 7)= RD5( 7)+(FQD5(1)*X2Y2-FQD4(1,3)*F03*AQX2-FQD4(1,3)*ACY2  
     *                 +FQD3(1,3)*F03                        )*XMDTY    
      RD5( 8)= RD5( 8)+(FQD5(2)*X2Y2-FQD4(2,3)*Q2C2+FQD3(2,3))*XMDT     
      RD5( 9)= RD5( 9)+(FQD5(3)*AQX2-FQD4(1,3)*AQX2-FQD4(3,3)           
     *                 +FQD3(1,3)                            )*XMDTY    
      RD5(10)= RD5(10)+(FQD5(4)*AQX2-FQD4(2,3)*F03*AQX2-FQD4(4,3)       
     *                 +FQD3(2,3)*F03                        )*XMDT     
      RD5(11)= RD5(11)+(FQD5(1)*ACY4-FQD4(1,3)*F06*ACY2                 
     *                 +FQD3(1,3)*F03                        )*XMDTX    
      RD5(12)= RD5(12)+(FQD5(2)*ACY2-FQD4(2,3)*F03           )*XMDTXY   
      RD5(13)= RD5(13)+(FQD5(3)*ACY2-FQD4(1,3)*ACY2-FQD4(3,3)           
     *                 +FQD3(1,3)                            )*XMDTX    
      RD5(14)= RD5(14)+(FQD5(4)-FQD4(2,3)*F03                )*XMDTXY   
      RD5(15)= RD5(15)+(FQD5(5)-FQD4(3,3)*F06+FQD3(1,3)*F03  )*XMDTX    
      RD5(16)= RD5(16)+(FQD5(1)*ACY4-FQD4(1,3)*F10*ACY2                 
     *                 +FQD3(1,3)*F15                        )*XMDTY    
      RD5(17)= RD5(17)+(FQD5(2)*ACY4-FQD4(2,3)*F06*ACY2                 
     *                 +FQD3(2,3)*F03                        )*XMDT     
      RD5(18)= RD5(18)+(FQD5(3)*ACY2-FQD4(1,3)*ACY2-FQD4(3,3)*F03       
     *                 +FQD3(1,3)*F03                        )*XMDTY    
      RD5(19)= RD5(19)+(FQD5(4)*ACY2-FQD4(2,3)*F03*ACY2-FQD4(4,3)       
     *                 +FQD3(2,3)*F03                        )*XMDT     
      RD5(20)= RD5(20)+(FQD5(5)-FQD4(3,3)*F06+FQD3(1,3)*F03  )*XMDTY    
      RD5(21)= RD5(21)+(FQD5(6)-FQD4(4,3)*F10+FQD3(2,3)*F15  )*XMDT     
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2R   *DECK INTK15                                           
C>                                                                      
C>    @brief   DPDS case                                                
C>                                                                      
C>    @details integration of the DPDS case                             
C>                                                                      
      SUBROUTINE INTK15(IKL)                                            
      use mx_limits, only: mxgsh,mxg2                                   
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
C GENERATE JTYPE=15 INTEGRALS                                           
C                                                                       
C                                                                       
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS,                      
     2                TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2),      
     3                D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2),      
     4                NGANGB                                            
C$omp threadprivate(/GEOMPQ/)
      COMMON /JMSGYH/ SQ(4)                                             
C$omp threadprivate(/JMSGYH/)
      COMMON /FQ08  / FQD(9),FQD0(5),FQD1(2,13),FQD2(3,16),FQD3(4,16),  
     2                FQD4(5,16),FQD5(66),FQD6(49),FQD7(24),FQD8( 9)    
C$omp threadprivate(/FQ08/)
      COMMON /KI2 / ACY,ACY2,AQX,AQX2,AQXY,Y03,Y04                      
C$omp threadprivate(/KI2/)
      COMMON /KI4 / RD0(5,5),RD1(3,40),RD2(6,56),RD3(10,52),RD4(15,42), 
     *              RD5(504),RD6(336),RD7(144),RD8(45)                  
C$omp threadprivate(/KI4/)
C                                                                       
      DIMENSION  WORK(12),FWK(15,4)                                     
C                                                                       
      PARAMETER (ZER=0.0D+00)                                           
      PARAMETER (F03=3.0D+00)                                           
      PARAMETER (F06=6.0D+00)                                           
      PARAMETER (F10=1.0D+01)                                           
      PARAMETER (F15=1.5D+01)                                           
C                                                                       
      IF(IKL.EQ.0) THEN                                                 
         DO I= 1, 2                                                     
            RD0(1,I)= ZER                                               
            RD0(2,I)= ZER                                               
            RD0(3,I)= ZER                                               
            RD0(4,I)= ZER                                               
            RD0(5,I)= ZER                                               
         ENDDO                                                          
         DO I= 1,13                                                     
            RD1(1,I)= ZER                                               
            RD1(2,I)= ZER                                               
            RD1(3,I)= ZER                                               
         ENDDO                                                          
         DO I= 1,15                                                     
            RD2(1,I)= ZER                                               
            RD2(2,I)= ZER                                               
            RD2(3,I)= ZER                                               
            RD2(4,I)= ZER                                               
            RD2(5,I)= ZER                                               
            RD2(6,I)= ZER                                               
         ENDDO                                                          
         DO I= 1, 9                                                     
            DO J= 1,10                                                  
               RD3(J,I)= ZER                                            
            ENDDO                                                       
         ENDDO                                                          
         DO I= 1, 4                                                     
            DO J= 1,15                                                  
               RD4(J,I)= ZER                                            
            ENDDO                                                       
         ENDDO                                                          
         DO J= 1,21                                                     
            RD5(J)= ZER                                                 
         ENDDO                                                          
C                                                                       
         RETURN                                                         
      ENDIF                                                             
C                                                                       
      XMD2= X43 *0.5D+00                                                
      XMD1= XMD2*SQ(1)                                                  
      XMD3= XMD2*SQ(2)                                                  
      XMD4= XMD2*XMD2                                                   
      XMD5= XMD4*SQ(2)                                                  
      XMDT= XMD4*XMD3                                                   
C                                                                       
      XMDTY=-XMDT*ACY                                                   
      XMDTX=-XMDT*AQX                                                   
      XMDTXY=XMDT*AQXY                                                  
C                                                                       
      Y33 = Y03 *Y03                                                    
      Y34 =-Y03 *Y04                                                    
      WORK( 1)= XMD1                                                    
      WORK( 2)= Y33 *SQ(1)                                              
      WORK( 3)=-XMD3*Y03                                                
      WORK( 4)= XMD3*Y04                                                
      WORK( 5)= Y33 *Y04*SQ(2)                                          
C                                                                       
      WORK( 6)=-XMD1*Y03                                                
      WORK( 7)= XMD5                                                    
      WORK( 8)= XMD3*Y33                                                
      WORK( 9)= XMD3*Y34                                                
C                                                                       
      WORK(10)= XMD4*SQ(1)                                              
      WORK(11)=-XMD5*Y03                                                
      WORK(12)= XMD5*Y04                                                
C                                                                       
      DO 010 I= 1, 2                                                    
         DO 010 J= 1, 5                                                 
            RD0(J,I)= RD0(J,I)+FQD0(I)*WORK(J)                          
  010 CONTINUE                                                          
C                                                                       
      DO I= 1, 3                                                        
         FWK(1,I)=-FQD1(1,I)*AQX                                        
         FWK(2,I)=-FQD1(1,I)*ACY                                        
         FWK(3,I)= FQD1(2,I)                                            
      ENDDO                                                             
      DO I= 1, 4                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,1)*WORK(I+ 5)                         
         RD1(2,I)= RD1(2,I)+FWK(2,1)*WORK(I+ 5)                         
         RD1(3,I)= RD1(3,I)+FWK(3,1)*WORK(I+ 5)                         
      ENDDO                                                             
      DO I= 5, 8                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,2)*WORK(I+ 1)                         
         RD1(2,I)= RD1(2,I)+FWK(2,2)*WORK(I+ 1)                         
         RD1(3,I)= RD1(3,I)+FWK(3,2)*WORK(I+ 1)                         
      ENDDO                                                             
      DO I= 9,13                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,3)*WORK(I- 8)                         
         RD1(2,I)= RD1(2,I)+FWK(2,3)*WORK(I- 8)                         
         RD1(3,I)= RD1(3,I)+FWK(3,3)*WORK(I- 8)                         
      ENDDO                                                             
C                                                                       
      DO I= 1, 4                                                        
         FWK(1,I)= FQD2(1,I)*AQX2-FQD1(1,I)                             
         FWK(2,I)= FQD2(1,I)*ACY2-FQD1(1,I)                             
         FWK(3,I)= FQD2(3,I)     -FQD1(1,I)                             
         FWK(4,I)= FQD2(1,I)*AQXY                                       
         FWK(5,I)=-FQD2(2,I)*AQX                                        
         FWK(6,I)=-FQD2(2,I)*ACY                                        
      ENDDO                                                             
      DO I= 1, 3                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,1)*WORK(I+ 9)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I= 4, 6                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,2)*WORK(I+ 6)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I= 7,10                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,3)*WORK(I- 1)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=11,15                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,4)*WORK(I-10)                      
         ENDDO                                                          
      ENDDO                                                             
C                                                                       
      DO I= 1, 2                                                        
         RD3( 1,I)= RD3( 1,I)+(FQD3(1,I)*AQX2-FQD2(1,I)*F03)*XMDTX      
         RD3( 2,I)= RD3( 2,I)+(FQD3(1,I)*AQX2-FQD2(1,I)    )*XMDTY      
         RD3( 3,I)= RD3( 3,I)+(FQD3(2,I)*AQX2-FQD2(2,I)    )*XMDT       
         RD3( 4,I)= RD3( 4,I)+(FQD3(1,I)*ACY2-FQD2(1,I)    )*XMDTX      
         RD3( 5,I)= RD3( 5,I)+ FQD3(2,I)                    *XMDTXY     
         RD3( 6,I)= RD3( 6,I)+(FQD3(3,I)     -FQD2(1,I)    )*XMDTX      
         RD3( 7,I)= RD3( 7,I)+(FQD3(1,I)*ACY2-FQD2(1,I)*F03)*XMDTY      
         RD3( 8,I)= RD3( 8,I)+(FQD3(2,I)*ACY2-FQD2(2,I)    )*XMDT       
         RD3( 9,I)= RD3( 9,I)+(FQD3(3,I)     -FQD2(1,I)    )*XMDTY      
         RD3(10,I)= RD3(10,I)+(FQD3(4,I)     -FQD2(2,I)*F03)*XMDT       
      ENDDO                                                             
      DO I= 3, 4                                                        
         FWK( 1,I)=-(FQD3(1,I)*AQX2-FQD2(1,I)*F03)*AQX                  
         FWK( 2,I)=-(FQD3(1,I)*AQX2-FQD2(1,I)    )*ACY                  
         FWK( 3,I)=  FQD3(2,I)*AQX2-FQD2(2,I)                           
         FWK( 4,I)=-(FQD3(1,I)*ACY2-FQD2(1,I)    )*AQX                  
         FWK( 5,I)=  FQD3(2,I)                    *AQXY                 
         FWK( 6,I)=-(FQD3(3,I)     -FQD2(1,I)    )*AQX                  
         FWK( 7,I)=-(FQD3(1,I)*ACY2-FQD2(1,I)*F03)*ACY                  
         FWK( 8,I)=  FQD3(2,I)*ACY2-FQD2(2,I)                           
         FWK( 9,I)=-(FQD3(3,I)     -FQD2(1,I)    )*ACY                  
         FWK(10,I)=  FQD3(4,I)     -FQD2(2,I)*F03                       
      ENDDO                                                             
      DO I= 3, 5                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,3)*WORK(I+ 7)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I= 6, 9                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,4)*WORK(I)                         
         ENDDO                                                          
      ENDDO                                                             
C                                                                       
      AQX4= AQX2*AQX2                                                   
      ACY4= ACY2*ACY2                                                   
      X2Y2= AQX2*ACY2                                                   
      Q2C2= AQX2+ACY2                                                   
      RD4( 1,1)= RD4( 1,1)+(FQD4(1,1)*AQX4-FQD3(1,3)*F06*AQX2           
     *                                    +FQD2(1,3)*F03     )*XMDT     
      RD4( 2,1)= RD4( 2,1)+(FQD4(1,1)*AQX2-FQD3(1,3)*F03     )*XMDTXY   
      RD4( 3,1)= RD4( 3,1)+(FQD4(2,1)*AQX2-FQD3(2,3)*F03     )*XMDTX    
      RD4( 4,1)= RD4( 4,1)+(FQD4(1,1)*X2Y2-FQD3(1,3)*Q2C2               
     *                                    +FQD2(1,3)         )*XMDT     
      RD4( 5,1)= RD4( 5,1)+(FQD4(2,1)*AQX2-FQD3(2,3)         )*XMDTY    
      RD4( 6,1)= RD4( 6,1)+(FQD4(3,1)*AQX2-FQD3(1,3)*AQX2               
     *                     -FQD3(3,3)+FQD2(1,3)              )*XMDT     
      RD4( 7,1)= RD4( 7,1)+(FQD4(1,1)*ACY2-FQD3(1,3)*F03     )*XMDTXY   
      RD4( 8,1)= RD4( 8,1)+(FQD4(2,1)*ACY2-FQD3(2,3)         )*XMDTX    
      RD4( 9,1)= RD4( 9,1)+(FQD4(3,1)     -FQD3(1,3)         )*XMDTXY   
      RD4(10,1)= RD4(10,1)+(FQD4(4,1)     -FQD3(2,3)*F03     )*XMDTX    
      RD4(11,1)= RD4(11,1)+(FQD4(1,1)*ACY4-FQD3(1,3)*F06*ACY2           
     *                     +FQD2(1,3)*F03                    )*XMDT     
      RD4(12,1)= RD4(12,1)+(FQD4(2,1)*ACY2-FQD3(2,3)*F03     )*XMDTY    
      RD4(13,1)= RD4(13,1)+(FQD4(3,1)*ACY2-FQD3(1,3)*ACY2               
     *                     -FQD3(3,3)+FQD2(1,3)              )*XMDT     
      RD4(14,1)= RD4(14,1)+(FQD4(4,1)     -FQD3(2,3)*F03     )*XMDTY    
      RD4(15,1)= RD4(15,1)+(FQD4(5,1)-FQD3(3,3)*F06                     
     *                     +FQD2(1,3)*F03                    )*XMDT     
      FWK( 1,2)=  FQD4(1,2)*AQX4-FQD3(1,4)*F06*AQX2+FQD2(1,4)*F03       
      FWK( 2,2)= (FQD4(1,2)*AQX2-FQD3(1,4)*F03               )*AQXY     
      FWK( 3,2)=-(FQD4(2,2)*AQX2-FQD3(2,4)*F03               )*AQX      
      FWK( 4,2)=  FQD4(1,2)*X2Y2-FQD3(1,4)*Q2C2    +FQD2(1,4)           
      FWK( 5,2)=-(FQD4(2,2)*AQX2-FQD3(2,4)                   )*ACY      
      FWK( 6,2)=  FQD4(3,2)*AQX2-FQD3(1,4)*AQX2-FQD3(3,4)+FQD2(1,4)     
      FWK( 7,2)= (FQD4(1,2)*ACY2-FQD3(1,4)*F03               )*AQXY     
      FWK( 8,2)=-(FQD4(2,2)*ACY2-FQD3(2,4)                   )*AQX      
      FWK( 9,2)= (FQD4(3,2)     -FQD3(1,4)                   )*AQXY     
      FWK(10,2)=-(FQD4(4,2)     -FQD3(2,4)*F03               )*AQX      
      FWK(11,2)=  FQD4(1,2)*ACY4-FQD3(1,4)*F06*ACY2+FQD2(1,4)*F03       
      FWK(12,2)=-(FQD4(2,2)*ACY2-FQD3(2,4)*F03               )*ACY      
      FWK(13,2)=  FQD4(3,2)*ACY2-FQD3(1,4)*ACY2-FQD3(3,4)+FQD2(1,4)     
      FWK(14,2)=-(FQD4(4,2)     -FQD3(2,4)*F03               )*ACY      
      FWK(15,2)=  FQD4(5,2)     -FQD3(3,4)*F06+FQD2(1,4)*F03            
      DO 410 I= 2, 4                                                    
         DO 410 J= 1,15                                                 
            RD4(J,I)= RD4(J,I)+FWK(J,2)*WORK(I+ 8)                      
  410 CONTINUE                                                          
C                                                                       
      RD5( 1)= RD5( 1)+(FQD5(1)*AQX4-FQD4(1,2)*F10*AQX2                 
     *                 +FQD3(1,4)*F15                        )*XMDTX    
      RD5( 2)= RD5( 2)+(FQD5(1)*AQX4-FQD4(1,2)*F06*AQX2                 
     *                 +FQD3(1,4)*F03                        )*XMDTY    
      RD5( 3)= RD5( 3)+(FQD5(2)*AQX4-FQD4(2,2)*F06*AQX2                 
     *                 +FQD3(2,4)*F03                        )*XMDT     
      RD5( 4)= RD5( 4)+(FQD5(1)*X2Y2-FQD4(1,2)*AQX2-FQD4(1,2)*F03*ACY2  
     *                 +FQD3(1,4)*F03                        )*XMDTX    
      RD5( 5)= RD5( 5)+(FQD5(2)*AQX2-FQD4(2,2)*F03           )*XMDTXY   
      RD5( 6)= RD5( 6)+(FQD5(3)*AQX2-FQD4(1,2)*AQX2-FQD4(3,2)*F03       
     *                 +FQD3(1,4)*F03                        )*XMDTX    
      RD5( 7)= RD5( 7)+(FQD5(1)*X2Y2-FQD4(1,2)*F03*AQX2-FQD4(1,2)*ACY2  
     *                 +FQD3(1,4)*F03                        )*XMDTY    
      RD5( 8)= RD5( 8)+(FQD5(2)*X2Y2-FQD4(2,2)*Q2C2+FQD3(2,4))*XMDT     
      RD5( 9)= RD5( 9)+(FQD5(3)*AQX2-FQD4(1,2)*AQX2-FQD4(3,2)           
     *                 +FQD3(1,4)                            )*XMDTY    
      RD5(10)= RD5(10)+(FQD5(4)*AQX2-FQD4(2,2)*F03*AQX2-FQD4(4,2)       
     *                 +FQD3(2,4)*F03                        )*XMDT     
      RD5(11)= RD5(11)+(FQD5(1)*ACY4-FQD4(1,2)*F06*ACY2                 
     *                 +FQD3(1,4)*F03                        )*XMDTX    
      RD5(12)= RD5(12)+(FQD5(2)*ACY2-FQD4(2,2)*F03           )*XMDTXY   
      RD5(13)= RD5(13)+(FQD5(3)*ACY2-FQD4(1,2)*ACY2-FQD4(3,2)           
     *                 +FQD3(1,4)                            )*XMDTX    
      RD5(14)= RD5(14)+(FQD5(4)-FQD4(2,2)*F03                )*XMDTXY   
      RD5(15)= RD5(15)+(FQD5(5)-FQD4(3,2)*F06+FQD3(1,4)*F03  )*XMDTX    
      RD5(16)= RD5(16)+(FQD5(1)*ACY4-FQD4(1,2)*F10*ACY2                 
     *                 +FQD3(1,4)*F15                        )*XMDTY    
      RD5(17)= RD5(17)+(FQD5(2)*ACY4-FQD4(2,2)*F06*ACY2                 
     *                 +FQD3(2,4)*F03                        )*XMDT     
      RD5(18)= RD5(18)+(FQD5(3)*ACY2-FQD4(1,2)*ACY2-FQD4(3,2)*F03       
     *                 +FQD3(1,4)*F03                        )*XMDTY    
      RD5(19)= RD5(19)+(FQD5(4)*ACY2-FQD4(2,2)*F03*ACY2-FQD4(4,2)       
     *                 +FQD3(2,4)*F03                        )*XMDT     
      RD5(20)= RD5(20)+(FQD5(5)-FQD4(3,2)*F06+FQD3(1,4)*F03  )*XMDTY    
      RD5(21)= RD5(21)+(FQD5(6)-FQD4(4,2)*F10+FQD3(2,4)*F15  )*XMDT     
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2R   *DECK INTK16                                           
C>                                                                      
C>    @brief   DPPP case                                                
C>                                                                      
C>    @details integration of the DPPP case                             
C>                                                                      
      SUBROUTINE INTK16(IKL)                                            
      use mx_limits, only: mxgsh,mxg2                                   
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
C GENERATE JTYPE=16 INTEGRALS                                           
C                                                                       
C                                                                       
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS,                      
     2                TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2),      
     3                D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2),      
     4                NGANGB                                            
C$omp threadprivate(/GEOMPQ/)
      COMMON /JMSGYH/ SQ(4)                                             
C$omp threadprivate(/JMSGYH/)
      COMMON /FQ08  / FQD(9),FQD0(5),FQD1(2,13),FQD2(3,16),FQD3(4,16),  
     2                FQD4(5,16),FQD5(66),FQD6(49),FQD7(24),FQD8( 9)    
C$omp threadprivate(/FQ08/)
      COMMON /KI2 / ACY,ACY2,AQX,AQX2,AQXY,Y03,Y04                      
C$omp threadprivate(/KI2/)
      COMMON /KI4 / RD0(5,5),RD1(3,40),RD2(6,56),RD3(10,52),RD4(15,42), 
     *              RD5(504),RD6(336),RD7(144),RD8(45)                  
C$omp threadprivate(/KI4/)
C                                                                       
      DIMENSION  WORK(12),FWK(15,10)                                    
C                                                                       
      PARAMETER (ZER=0.0D+00)                                           
      PARAMETER (F03=3.0D+00)                                           
      PARAMETER (F06=6.0D+00)                                           
      PARAMETER (F10=1.0D+01)                                           
      PARAMETER (F15=1.5D+01)                                           
C                                                                       
      IF(IKL.EQ.0) THEN                                                 
         DO I= 1, 5                                                     
            RD0(1,I)= ZER                                               
            RD0(2,I)= ZER                                               
            RD0(3,I)= ZER                                               
            RD0(4,I)= ZER                                               
            RD0(5,I)= ZER                                               
         ENDDO                                                          
         DO I= 1,40                                                     
            RD1(1,I)= ZER                                               
            RD1(2,I)= ZER                                               
            RD1(3,I)= ZER                                               
         ENDDO                                                          
         DO I= 1,36                                                     
            RD2(1,I)= ZER                                               
            RD2(2,I)= ZER                                               
            RD2(3,I)= ZER                                               
            RD2(4,I)= ZER                                               
            RD2(5,I)= ZER                                               
            RD2(6,I)= ZER                                               
         ENDDO                                                          
         DO I= 1,21                                                     
            DO J= 1,10                                                  
               RD3(J,I)= ZER                                            
            ENDDO                                                       
         ENDDO                                                          
         DO I= 1, 7                                                     
            DO J= 1,15                                                  
               RD4(J,I)= ZER                                            
            ENDDO                                                       
         ENDDO                                                          
         DO J= 1,21                                                     
            RD5(J)= ZER                                                 
         ENDDO                                                          
C                                                                       
         RETURN                                                         
      ENDIF                                                             
C                                                                       
      XMD2= X43 *0.5D+00                                                
      XMD1= XMD2*SQ(1)                                                  
      XMD3= XMD2*SQ(2)                                                  
      XMD4= XMD2*XMD2                                                   
      XMD5= XMD4*SQ(2)                                                  
      XMDT= XMD4*XMD3                                                   
C                                                                       
      XMDTY=-XMDT*ACY                                                   
      XMDTX=-XMDT*AQX                                                   
      XMDTXY=XMDT*AQXY                                                  
C                                                                       
      Y33 = Y03 *Y03                                                    
      Y34 =-Y03 *Y04                                                    
      WORK( 1)= XMD1                                                    
      WORK( 2)= Y33 *SQ(1)                                              
      WORK( 3)=-XMD3*Y03                                                
      WORK( 4)= XMD3*Y04                                                
      WORK( 5)= Y33 *Y04*SQ(2)                                          
C                                                                       
      WORK( 6)=-XMD1*Y03                                                
      WORK( 7)= XMD5                                                    
      WORK( 8)= XMD3*Y33                                                
      WORK( 9)= XMD3*Y34                                                
C                                                                       
      WORK(10)= XMD4*SQ(1)                                              
      WORK(11)=-XMD5*Y03                                                
      WORK(12)= XMD5*Y04                                                
C                                                                       
      DO 010 I= 1, 5                                                    
         DO 010 J= 1, 5                                                 
            RD0(J,I)= RD0(J,I)+FQD0(I)*WORK(J)                          
  010 CONTINUE                                                          
C                                                                       
      DO I= 1, 8                                                        
         FWK(1,I)=-FQD1(1,I)*AQX                                        
         FWK(2,I)=-FQD1(1,I)*ACY                                        
         FWK(3,I)= FQD1(2,I)                                            
      ENDDO                                                             
      FQD11 = FQD1(1,5)*RAB +FQD1(1,8)                                  
      FQD12 = FQD1(2,5)*RAB +FQD1(2,8)                                  
         FWK(1,9)=-FQD11*AQX                                            
         FWK(2,9)=-FQD11*ACY                                            
         FWK(3,9)= FQD12                                                
      DO I= 1, 4                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,1)*WORK(I+ 5)                         
         RD1(2,I)= RD1(2,I)+FWK(2,1)*WORK(I+ 5)                         
         RD1(3,I)= RD1(3,I)+FWK(3,1)*WORK(I+ 5)                         
      ENDDO                                                             
      DO I= 5, 8                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,2)*WORK(I+ 1)                         
         RD1(2,I)= RD1(2,I)+FWK(2,2)*WORK(I+ 1)                         
         RD1(3,I)= RD1(3,I)+FWK(3,2)*WORK(I+ 1)                         
      ENDDO                                                             
      DO I= 9,12                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,3)*WORK(I- 3)                         
         RD1(2,I)= RD1(2,I)+FWK(2,3)*WORK(I- 3)                         
         RD1(3,I)= RD1(3,I)+FWK(3,3)*WORK(I- 3)                         
      ENDDO                                                             
      DO I=13,16                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,4)*WORK(I- 7)                         
         RD1(2,I)= RD1(2,I)+FWK(2,4)*WORK(I- 7)                         
         RD1(3,I)= RD1(3,I)+FWK(3,4)*WORK(I- 7)                         
      ENDDO                                                             
      DO I=17,20                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,5)*WORK(I-11)                         
         RD1(2,I)= RD1(2,I)+FWK(2,5)*WORK(I-11)                         
         RD1(3,I)= RD1(3,I)+FWK(3,5)*WORK(I-11)                         
      ENDDO                                                             
      DO I=21,25                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,6)*WORK(I-20)                         
         RD1(2,I)= RD1(2,I)+FWK(2,6)*WORK(I-20)                         
         RD1(3,I)= RD1(3,I)+FWK(3,6)*WORK(I-20)                         
      ENDDO                                                             
      DO I=26,30                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,7)*WORK(I-25)                         
         RD1(2,I)= RD1(2,I)+FWK(2,7)*WORK(I-25)                         
         RD1(3,I)= RD1(3,I)+FWK(3,7)*WORK(I-25)                         
      ENDDO                                                             
      DO I=31,35                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,8)*WORK(I-30)                         
         RD1(2,I)= RD1(2,I)+FWK(2,8)*WORK(I-30)                         
         RD1(3,I)= RD1(3,I)+FWK(3,8)*WORK(I-30)                         
      ENDDO                                                             
      DO I=36,40                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,9)*WORK(I-35)                         
         RD1(2,I)= RD1(2,I)+FWK(2,9)*WORK(I-35)                         
         RD1(3,I)= RD1(3,I)+FWK(3,9)*WORK(I-35)                         
      ENDDO                                                             
C                                                                       
      DO I= 1, 9                                                        
         FWK(1,I)= FQD2(1,I)*AQX2-FQD1(1,I)                             
         FWK(2,I)= FQD2(1,I)*ACY2-FQD1(1,I)                             
         FWK(3,I)= FQD2(3,I)     -FQD1(1,I)                             
         FWK(4,I)= FQD2(1,I)*AQXY                                       
         FWK(5,I)=-FQD2(2,I)*AQX                                        
         FWK(6,I)=-FQD2(2,I)*ACY                                        
      ENDDO                                                             
      FQD21 = FQD2(1,5)*RAB +FQD2(1,8)                                  
      FQD22 = FQD2(2,5)*RAB +FQD2(2,8)                                  
      FQD23 = FQD2(3,5)*RAB +FQD2(3,8)                                  
         FWK(1,10)= FQD21*AQX2-FQD11                                    
         FWK(2,10)= FQD21*ACY2-FQD11                                    
         FWK(3,10)= FQD23     -FQD11                                    
         FWK(4,10)= FQD21*AQXY                                          
         FWK(5,10)=-FQD22*AQX                                           
         FWK(6,10)=-FQD22*ACY                                           
      DO I= 1, 3                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,1)*WORK(I+ 9)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I= 4, 6                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,2)*WORK(I+ 6)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I= 7, 9                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,3)*WORK(I+ 3)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=10,12                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,4)*WORK(I)                         
         ENDDO                                                          
      ENDDO                                                             
      DO I=13,15                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,5)*WORK(I- 3)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=16,19                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,6)*WORK(I-10)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=20,23                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,7)*WORK(I-14)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=24,27                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,8)*WORK(I-18)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=28,31                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,10)*WORK(I-22)                     
         ENDDO                                                          
      ENDDO                                                             
      DO I=32,36                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,9)*WORK(I-31)                      
         ENDDO                                                          
      ENDDO                                                             
C                                                                       
      DO I= 1, 5                                                        
         RD3( 1,I)= RD3( 1,I)+(FQD3(1,I)*AQX2-FQD2(1,I)*F03)*XMDTX      
         RD3( 2,I)= RD3( 2,I)+(FQD3(1,I)*AQX2-FQD2(1,I)    )*XMDTY      
         RD3( 3,I)= RD3( 3,I)+(FQD3(2,I)*AQX2-FQD2(2,I)    )*XMDT       
         RD3( 4,I)= RD3( 4,I)+(FQD3(1,I)*ACY2-FQD2(1,I)    )*XMDTX      
         RD3( 5,I)= RD3( 5,I)+ FQD3(2,I)                    *XMDTXY     
         RD3( 6,I)= RD3( 6,I)+(FQD3(3,I)     -FQD2(1,I)    )*XMDTX      
         RD3( 7,I)= RD3( 7,I)+(FQD3(1,I)*ACY2-FQD2(1,I)*F03)*XMDTY      
         RD3( 8,I)= RD3( 8,I)+(FQD3(2,I)*ACY2-FQD2(2,I)    )*XMDT       
         RD3( 9,I)= RD3( 9,I)+(FQD3(3,I)     -FQD2(1,I)    )*XMDTY      
         RD3(10,I)= RD3(10,I)+(FQD3(4,I)     -FQD2(2,I)*F03)*XMDT       
      ENDDO                                                             
      FQD2(1,5)= FQD21                                                  
      FQD2(2,5)= FQD22                                                  
      FQD3(1,5)= FQD3(1,5)*RAB+FQD3(1,8)                                
      FQD3(2,5)= FQD3(2,5)*RAB+FQD3(2,8)                                
      FQD3(3,5)= FQD3(3,5)*RAB+FQD3(3,8)                                
      FQD3(4,5)= FQD3(4,5)*RAB+FQD3(4,8)                                
      DO I= 5, 9                                                        
         FWK( 1,I)=-(FQD3(1,I)*AQX2-FQD2(1,I)*F03)*AQX                  
         FWK( 2,I)=-(FQD3(1,I)*AQX2-FQD2(1,I)    )*ACY                  
         FWK( 3,I)=  FQD3(2,I)*AQX2-FQD2(2,I)                           
         FWK( 4,I)=-(FQD3(1,I)*ACY2-FQD2(1,I)    )*AQX                  
         FWK( 5,I)=  FQD3(2,I)                    *AQXY                 
         FWK( 6,I)=-(FQD3(3,I)     -FQD2(1,I)    )*AQX                  
         FWK( 7,I)=-(FQD3(1,I)*ACY2-FQD2(1,I)*F03)*ACY                  
         FWK( 8,I)=  FQD3(2,I)*ACY2-FQD2(2,I)                           
         FWK( 9,I)=-(FQD3(3,I)     -FQD2(1,I)    )*ACY                  
         FWK(10,I)=  FQD3(4,I)     -FQD2(2,I)*F03                       
      ENDDO                                                             
      DO I= 6, 8                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,5)*WORK(I+ 4)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I= 9,11                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,6)*WORK(I+ 1)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=12,14                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,7)*WORK(I- 2)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=15,17                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,8)*WORK(I- 5)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=18,21                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,9)*WORK(I-12)                      
         ENDDO                                                          
      ENDDO                                                             
C                                                                       
      DO J= 1, 5                                                        
         FQD4(J,1)= FQD4(J,1)*RAB +FQD4(J,4)                            
      ENDDO                                                             
      AQX4= AQX2*AQX2                                                   
      ACY4= ACY2*ACY2                                                   
      X2Y2= AQX2*ACY2                                                   
      Q2C2= AQX2+ACY2                                                   
      DO I= 1, 4                                                        
         RD4( 1,I)= RD4( 1,I)+(FQD4(1,I)*AQX4-FQD3(1,I+4)*F06*AQX2      
     *                        +FQD2(1,I+4)*F03                )*XMDT    
         RD4( 2,I)= RD4( 2,I)+(FQD4(1,I)*AQX2-FQD3(1,I+4)*F03 )*XMDTXY  
         RD4( 3,I)= RD4( 3,I)+(FQD4(2,I)*AQX2-FQD3(2,I+4)*F03 )*XMDTX   
         RD4( 4,I)= RD4( 4,I)+(FQD4(1,I)*X2Y2-FQD3(1,I+4)*Q2C2          
     *                                       +FQD2(1,I+4)     )*XMDT    
         RD4( 5,I)= RD4( 5,I)+(FQD4(2,I)*AQX2-FQD3(2,I+4)     )*XMDTY   
         RD4( 6,I)= RD4( 6,I)+(FQD4(3,I)*AQX2-FQD3(1,I+4)*AQX2          
     *                        -FQD3(3,I+4)+FQD2(1,I+4)        )*XMDT    
         RD4( 7,I)= RD4( 7,I)+(FQD4(1,I)*ACY2-FQD3(1,I+4)*F03 )*XMDTXY  
         RD4( 8,I)= RD4( 8,I)+(FQD4(2,I)*ACY2-FQD3(2,I+4)     )*XMDTX   
         RD4( 9,I)= RD4( 9,I)+(FQD4(3,I)     -FQD3(1,I+4)     )*XMDTXY  
         RD4(10,I)= RD4(10,I)+(FQD4(4,I)     -FQD3(2,I+4)*F03 )*XMDTX   
         RD4(11,I)= RD4(11,I)+(FQD4(1,I)*ACY4-FQD3(1,I+4)*F06*ACY2      
     *                        +FQD2(1,I+4)*F03                )*XMDT    
         RD4(12,I)= RD4(12,I)+(FQD4(2,I)*ACY2-FQD3(2,I+4)*F03 )*XMDTY   
         RD4(13,I)= RD4(13,I)+(FQD4(3,I)*ACY2-FQD3(1,I+4)*ACY2          
     *                        -FQD3(3,I+4)+FQD2(1,I+4)        )*XMDT    
         RD4(14,I)= RD4(14,I)+(FQD4(4,I)     -FQD3(2,I+4)*F03 )*XMDTY   
         RD4(15,I)= RD4(15,I)+(FQD4(5,I)-FQD3(3,I+4)*F06                
     *                        +FQD2(1,I+4)*F03                )*XMDT    
      ENDDO                                                             
      FWK( 1,5)=  FQD4(1,5)*AQX4-FQD3(1,9)*F06*AQX2+FQD2(1,9)*F03       
      FWK( 2,5)= (FQD4(1,5)*AQX2-FQD3(1,9)*F03                )*AQXY    
      FWK( 3,5)=-(FQD4(2,5)*AQX2-FQD3(2,9)*F03                )*AQX     
      FWK( 4,5)=  FQD4(1,5)*X2Y2-FQD3(1,9)*Q2C2    +FQD2(1,9)           
      FWK( 5,5)=-(FQD4(2,5)*AQX2-FQD3(2,9)                    )*ACY     
      FWK( 6,5)=  FQD4(3,5)*AQX2-FQD3(1,9)*AQX2-FQD3(3,9)+FQD2(1,9)     
      FWK( 7,5)= (FQD4(1,5)*ACY2-FQD3(1,9)*F03                )*AQXY    
      FWK( 8,5)=-(FQD4(2,5)*ACY2-FQD3(2,9)                    )*AQX     
      FWK( 9,5)= (FQD4(3,5)     -FQD3(1,9)                    )*AQXY    
      FWK(10,5)=-(FQD4(4,5)     -FQD3(2,9)*F03                )*AQX     
      FWK(11,5)=  FQD4(1,5)*ACY4-FQD3(1,9)*F06*ACY2+FQD2(1,9)*F03       
      FWK(12,5)=-(FQD4(2,5)*ACY2-FQD3(2,9)*F03                )*ACY     
      FWK(13,5)=  FQD4(3,5)*ACY2-FQD3(1,9)*ACY2-FQD3(3,9)+FQD2(1,9)     
      FWK(14,5)=-(FQD4(4,5)     -FQD3(2,9)*F03                )*ACY     
      FWK(15,5)=  FQD4(5,5)     -FQD3(3,9)*F06+FQD2(1,9)*F03            
      DO 420 I= 5, 7                                                    
         DO 420 J= 1,15                                                 
            RD4(J,I)= RD4(J,I)+FWK(J,5)*WORK(I+ 5)                      
  420 CONTINUE                                                          
C                                                                       
      RD5( 1)= RD5( 1)+(FQD5(1)*AQX4-FQD4(1,5)*F10*AQX2                 
     *                 +FQD3(1,9)*F15                         )*XMDTX   
      RD5( 2)= RD5( 2)+(FQD5(1)*AQX4-FQD4(1,5)*F06*AQX2                 
     *                 +FQD3(1,9)*F03                         )*XMDTY   
      RD5( 3)= RD5( 3)+(FQD5(2)*AQX4-FQD4(2,5)*F06*AQX2                 
     *                 +FQD3(2,9)*F03                         )*XMDT    
      RD5( 4)= RD5( 4)+(FQD5(1)*X2Y2-FQD4(1,5)*AQX2-FQD4(1,5)*F03*ACY2  
     *                 +FQD3(1,9)*F03                         )*XMDTX   
      RD5( 5)= RD5( 5)+(FQD5(2)*AQX2-FQD4(2,5)*F03            )*XMDTXY  
      RD5( 6)= RD5( 6)+(FQD5(3)*AQX2-FQD4(1,5)*AQX2-FQD4(3,5)*F03       
     *                 +FQD3(1,9)*F03                         )*XMDTX   
      RD5( 7)= RD5( 7)+(FQD5(1)*X2Y2-FQD4(1,5)*F03*AQX2-FQD4(1,5)*ACY2  
     *                 +FQD3(1,9)*F03                         )*XMDTY   
      RD5( 8)= RD5( 8)+(FQD5(2)*X2Y2-FQD4(2,5)*Q2C2+FQD3(2,9) )*XMDT    
      RD5( 9)= RD5( 9)+(FQD5(3)*AQX2-FQD4(1,5)*AQX2-FQD4(3,5)           
     *                 +FQD3(1,9)                             )*XMDTY   
      RD5(10)= RD5(10)+(FQD5(4)*AQX2-FQD4(2,5)*F03*AQX2-FQD4(4,5)       
     *                 +FQD3(2,9)*F03                         )*XMDT    
      RD5(11)= RD5(11)+(FQD5(1)*ACY4-FQD4(1,5)*F06*ACY2                 
     *                 +FQD3(1,9)*F03                         )*XMDTX   
      RD5(12)= RD5(12)+(FQD5(2)*ACY2-FQD4(2,5)*F03            )*XMDTXY  
      RD5(13)= RD5(13)+(FQD5(3)*ACY2-FQD4(1,5)*ACY2-FQD4(3,5)           
     *                 +FQD3(1,9)                             )*XMDTX   
      RD5(14)= RD5(14)+(FQD5(4)-FQD4(2,5)*F03                 )*XMDTXY  
      RD5(15)= RD5(15)+(FQD5(5)-FQD4(3,5)*F06+FQD3(1,9)*F03   )*XMDTX   
      RD5(16)= RD5(16)+(FQD5(1)*ACY4-FQD4(1,5)*F10*ACY2                 
     *                 +FQD3(1,9)*F15                         )*XMDTY   
      RD5(17)= RD5(17)+(FQD5(2)*ACY4-FQD4(2,5)*F06*ACY2                 
     *                 +FQD3(2,9)*F03                         )*XMDT    
      RD5(18)= RD5(18)+(FQD5(3)*ACY2-FQD4(1,5)*ACY2-FQD4(3,5)*F03       
     *                 +FQD3(1,9)*F03                         )*XMDTY   
      RD5(19)= RD5(19)+(FQD5(4)*ACY2-FQD4(2,5)*F03*ACY2-FQD4(4,5)       
     *                 +FQD3(2,9)*F03                         )*XMDT    
      RD5(20)= RD5(20)+(FQD5(5)-FQD4(3,5)*F06+FQD3(1,9)*F03   )*XMDTY   
      RD5(21)= RD5(21)+(FQD5(6)-FQD4(4,5)*F10+FQD3(2,9)*F15   )*XMDT    
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2R   *DECK INTK17                                           
C>                                                                      
C>    @brief   DDDS case                                                
C>                                                                      
C>    @details integration of the DDDS case                             
C>                                                                      
      SUBROUTINE INTK17(IKL)                                            
      use mx_limits, only: mxgsh,mxg2                                   
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
C GENERATE JTYPE=17 INTEGRALS                                           
C                                                                       
C                                                                       
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS,                      
     2                TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2),      
     3                D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2),      
     4                NGANGB                                            
C$omp threadprivate(/GEOMPQ/)
      COMMON /JMSGYH/ SQ(4)                                             
C$omp threadprivate(/JMSGYH/)
      COMMON /FQ08  / FQD(9),FQD0(5),FQD1(2,13),FQD2(3,16),FQD3(4,16),  
     2                FQD4(5,16),FQD5(6,11),FQD6(49),FQD7(24),FQD8( 9)  
C$omp threadprivate(/FQ08/)
      COMMON /KI2 / ACY,ACY2,AQX,AQX2,AQXY,Y03,Y04                      
C$omp threadprivate(/KI2/)
      COMMON /KI4 / RD0(5,5),RD1(3,40),RD2(6,56),RD3(10,52),RD4(15,42), 
     *              RD5(21,24),RD6(336),RD7(144),RD8(45)                
C$omp threadprivate(/KI4/)
C                                                                       
      DIMENSION  WORK(15),FWK(21,4),FW6(28)                             
      DIMENSION  FCU(45,8),FCC(45,8)                                    
C                                                                       
      PARAMETER (ZER=0.0D+00)                                           
      PARAMETER (F03=3.0D+00)                                           
      PARAMETER (F06=6.0D+00)                                           
      PARAMETER (F10=1.0D+01)                                           
      PARAMETER (F15=1.5D+01)                                           
C                                                                       
      IF(IKL.EQ.0) THEN                                                 
         DO I= 1, 2                                                     
            RD0(1,I)= ZER                                               
            RD0(2,I)= ZER                                               
            RD0(3,I)= ZER                                               
            RD0(4,I)= ZER                                               
            RD0(5,I)= ZER                                               
         ENDDO                                                          
         DO I= 1,13                                                     
            RD1(1,I)= ZER                                               
            RD1(2,I)= ZER                                               
            RD1(3,I)= ZER                                               
         ENDDO                                                          
         DO I= 1,17                                                     
            RD2(1,I)= ZER                                               
            RD2(2,I)= ZER                                               
            RD2(3,I)= ZER                                               
            RD2(4,I)= ZER                                               
            RD2(5,I)= ZER                                               
            RD2(6,I)= ZER                                               
         ENDDO                                                          
         DO I= 1,12                                                     
            DO J= 1,10                                                  
               RD3(J,I)= ZER                                            
            ENDDO                                                       
         ENDDO                                                          
         DO I= 1, 8                                                     
            DO J= 1,15                                                  
               RD4(J,I)= ZER                                            
            ENDDO                                                       
         ENDDO                                                          
         DO I= 1, 3                                                     
            DO J= 1,21                                                  
               RD5(J,I)= ZER                                            
            ENDDO                                                       
         ENDDO                                                          
         DO J= 1,28                                                     
            RD6(J)= ZER                                                 
         ENDDO                                                          
C                                                                       
         RETURN                                                         
      ENDIF                                                             
C                                                                       
      XMD2= X43 *0.5D+00                                                
      XMD3= XMD2*SQ(3)                                                  
      XMD4= XMD3*XMD2                                                   
      XMD6= XMD4*XMD2                                                   
      XMDT= XMD6*XMD2                                                   
      XMD2= XMD3                                                        
C                                                                       
      XMDTY=-XMDT*ACY                                                   
      XMDTX=-XMDT*AQX                                                   
      XMDTXY=XMDT*AQXY                                                  
C                                                                       
      Y33 = Y03 *Y03                                                    
      Y34 =-Y03 *Y04                                                    
      Y44 = Y04 *Y04                                                    
      WORK( 1)= XMD4                                                    
      WORK( 2)= XMD2*Y33                                                
      WORK( 3)= XMD2*Y34                                                
      WORK( 4)= XMD2*Y44                                                
      WORK( 5)= Y33 *Y44*SQ(3)                                          
C                                                                       
      WORK( 6)=-XMD4*Y03                                                
      WORK( 7)= XMD4*Y04                                                
      WORK( 8)= XMD2*Y33*Y04                                            
      WORK( 9)= XMD2*Y34*Y04                                            
C                                                                       
      WORK(10)= XMD6                                                    
      WORK(11)= XMD4*Y33                                                
      WORK(12)= XMD4*Y34                                                
      WORK(13)= XMD4*Y44                                                
C                                                                       
      WORK(14)=-XMD6*Y03                                                
      WORK(15)= XMD6*Y04                                                
C                                                                       
      CALL FCUFCC(6,XMDT,FCU,FCC)                                       
C                                                                       
      DO 010 I= 1, 2                                                    
         DO 010 J= 1, 5                                                 
            RD0(J,I)= RD0(J,I)+FQD0(I)*WORK(J)                          
  010 CONTINUE                                                          
C                                                                       
      DO I= 1, 3                                                        
         FWK(1,I)=-FQD1(1,I)*AQX                                        
         FWK(2,I)=-FQD1(1,I)*ACY                                        
         FWK(3,I)= FQD1(2,I)                                            
      ENDDO                                                             
      DO I= 1, 4                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,1)*WORK(I+ 5)                         
         RD1(2,I)= RD1(2,I)+FWK(2,1)*WORK(I+ 5)                         
         RD1(3,I)= RD1(3,I)+FWK(3,1)*WORK(I+ 5)                         
      ENDDO                                                             
      DO I= 5, 8                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,2)*WORK(I+ 1)                         
         RD1(2,I)= RD1(2,I)+FWK(2,2)*WORK(I+ 1)                         
         RD1(3,I)= RD1(3,I)+FWK(3,2)*WORK(I+ 1)                         
      ENDDO                                                             
      DO I= 9,13                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,3)*WORK(I- 8)                         
         RD1(2,I)= RD1(2,I)+FWK(2,3)*WORK(I- 8)                         
         RD1(3,I)= RD1(3,I)+FWK(3,3)*WORK(I- 8)                         
      ENDDO                                                             
C                                                                       
      DO I= 1, 4                                                        
         FWK(1,I)= FQD2(1,I)*AQX2-FQD1(1,I)                             
         FWK(2,I)= FQD2(1,I)*ACY2-FQD1(1,I)                             
         FWK(3,I)= FQD2(3,I)     -FQD1(1,I)                             
         FWK(4,I)= FQD2(1,I)*AQXY                                       
         FWK(5,I)=-FQD2(2,I)*AQX                                        
         FWK(6,I)=-FQD2(2,I)*ACY                                        
      ENDDO                                                             
      DO I= 1, 4                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,1)*WORK(I+ 9)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I= 5, 8                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,2)*WORK(I+ 5)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I= 9,12                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,3)*WORK(I- 3)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=13,17                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,4)*WORK(I-12)                      
         ENDDO                                                          
      ENDDO                                                             
C                                                                       
      DO I= 1, 4                                                        
         FWK( 1,I)=-(FQD3(1,I)*AQX2-FQD2(1,I)*F03)*AQX                  
         FWK( 2,I)=-(FQD3(1,I)*AQX2-FQD2(1,I)    )*ACY                  
         FWK( 3,I)=  FQD3(2,I)*AQX2-FQD2(2,I)                           
         FWK( 4,I)=-(FQD3(1,I)*ACY2-FQD2(1,I)    )*AQX                  
         FWK( 5,I)=  FQD3(2,I)                    *AQXY                 
         FWK( 6,I)=-(FQD3(3,I)     -FQD2(1,I)    )*AQX                  
         FWK( 7,I)=-(FQD3(1,I)*ACY2-FQD2(1,I)*F03)*ACY                  
         FWK( 8,I)=  FQD3(2,I)*ACY2-FQD2(2,I)                           
         FWK( 9,I)=-(FQD3(3,I)     -FQD2(1,I)    )*ACY                  
         FWK(10,I)=  FQD3(4,I)     -FQD2(2,I)*F03                       
      ENDDO                                                             
      DO I= 1, 2                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,1)*WORK(I+13)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I= 3, 4                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,2)*WORK(I+11)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I= 5, 8                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,3)*WORK(I+ 5)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I= 9,12                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,4)*WORK(I- 3)                      
         ENDDO                                                          
      ENDDO                                                             
C                                                                       
      AQX4= AQX2*AQX2                                                   
      ACY4= ACY2*ACY2                                                   
      X2Y2= AQX2*ACY2                                                   
      Q2C2= AQX2+ACY2                                                   
      DO I= 1, 2                                                        
         RD4( 1,I)= RD4( 1,I)+(FQD4(1,I)*AQX4-FQD3(1,I)*F06*AQX2        
     *                        +FQD2(1,I)*F03                )*XMDT      
         RD4( 2,I)= RD4( 2,I)+(FQD4(1,I)*AQX2-FQD3(1,I)*F03 )*XMDTXY    
         RD4( 3,I)= RD4( 3,I)+(FQD4(2,I)*AQX2-FQD3(2,I)*F03 )*XMDTX     
         RD4( 4,I)= RD4( 4,I)+(FQD4(1,I)*X2Y2-FQD3(1,I)*Q2C2            
     *                                       +FQD2(1,I)     )*XMDT      
         RD4( 5,I)= RD4( 5,I)+(FQD4(2,I)*AQX2-FQD3(2,I)     )*XMDTY     
         RD4( 6,I)= RD4( 6,I)+(FQD4(3,I)*AQX2-FQD3(1,I)*AQX2            
     *                        -FQD3(3,I)+FQD2(1,I)          )*XMDT      
         RD4( 7,I)= RD4( 7,I)+(FQD4(1,I)*ACY2-FQD3(1,I)*F03 )*XMDTXY    
         RD4( 8,I)= RD4( 8,I)+(FQD4(2,I)*ACY2-FQD3(2,I)     )*XMDTX     
         RD4( 9,I)= RD4( 9,I)+(FQD4(3,I)     -FQD3(1,I)     )*XMDTXY    
         RD4(10,I)= RD4(10,I)+(FQD4(4,I)     -FQD3(2,I)*F03 )*XMDTX     
         RD4(11,I)= RD4(11,I)+(FQD4(1,I)*ACY4-FQD3(1,I)*F06*ACY2        
     *                        +FQD2(1,I)*F03                )*XMDT      
         RD4(12,I)= RD4(12,I)+(FQD4(2,I)*ACY2-FQD3(2,I)*F03 )*XMDTY     
         RD4(13,I)= RD4(13,I)+(FQD4(3,I)*ACY2-FQD3(1,I)*ACY2            
     *                        -FQD3(3,I)+FQD2(1,I)          )*XMDT      
         RD4(14,I)= RD4(14,I)+(FQD4(4,I)     -FQD3(2,I)*F03 )*XMDTY     
         RD4(15,I)= RD4(15,I)+(FQD4(5,I)-FQD3(3,I)*F06                  
     *                        +FQD2(1,I)*F03                )*XMDT      
      ENDDO                                                             
      DO I= 3, 4                                                        
         FWK( 1,I)=  FQD4(1,I)*AQX4-FQD3(1,I)*F06*AQX2+FQD2(1,I)*F03    
         FWK( 2,I)= (FQD4(1,I)*AQX2-FQD3(1,I)*F03           )*AQXY      
         FWK( 3,I)=-(FQD4(2,I)*AQX2-FQD3(2,I)*F03           )*AQX       
         FWK( 4,I)=  FQD4(1,I)*X2Y2-FQD3(1,I)*Q2C2    +FQD2(1,I)        
         FWK( 5,I)=-(FQD4(2,I)*AQX2-FQD3(2,I)               )*ACY       
         FWK( 6,I)=  FQD4(3,I)*AQX2-FQD3(1,I)*AQX2-FQD3(3,I)+FQD2(1,I)  
         FWK( 7,I)= (FQD4(1,I)*ACY2-FQD3(1,I)*F03           )*AQXY      
         FWK( 8,I)=-(FQD4(2,I)*ACY2-FQD3(2,I)               )*AQX       
         FWK( 9,I)= (FQD4(3,I)     -FQD3(1,I)               )*AQXY      
         FWK(10,I)=-(FQD4(4,I)     -FQD3(2,I)*F03           )*AQX       
         FWK(11,I)=  FQD4(1,I)*ACY4-FQD3(1,I)*F06*ACY2+FQD2(1,I)*F03    
         FWK(12,I)=-(FQD4(2,I)*ACY2-FQD3(2,I)*F03           )*ACY       
         FWK(13,I)=  FQD4(3,I)*ACY2-FQD3(1,I)*ACY2-FQD3(3,I)+FQD2(1,I)  
         FWK(14,I)=-(FQD4(4,I)     -FQD3(2,I)*F03           )*ACY       
         FWK(15,I)=  FQD4(5,I)     -FQD3(3,I)*F06+FQD2(1,I)*F03         
      ENDDO                                                             
      DO I= 3, 4                                                        
         DO J= 1,15                                                     
            RD4(J,I)= RD4(J,I)+FWK(J,3)*WORK(I+11)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I= 5, 8                                                        
         DO J= 1,15                                                     
            RD4(J,I)= RD4(J,I)+FWK(J,4)*WORK(I+ 5)                      
         ENDDO                                                          
      ENDDO                                                             
C                                                                       
      RD5( 1,1)= RD5( 1,1)+(FQD5(1,1)*AQX4-FQD4(1,3)*F10*AQX2           
     *                     +FQD3(1,3)*F15                    )*XMDTX    
      RD5( 2,1)= RD5( 2,1)+(FQD5(1,1)*AQX4-FQD4(1,3)*F06*AQX2           
     *                     +FQD3(1,3)*F03                    )*XMDTY    
      RD5( 3,1)= RD5( 3,1)+(FQD5(2,1)*AQX4-FQD4(2,3)*F06*AQX2           
     *                     +FQD3(2,3)*F03                    )*XMDT     
      RD5( 4,1)= RD5( 4,1)+(FQD5(1,1)*X2Y2-FQD4(1,3)*AQX2               
     *                     -FQD4(1,3)*F03*ACY2+FQD3(1,3)*F03 )*XMDTX    
      RD5( 5,1)= RD5( 5,1)+(FQD5(2,1)*AQX2-FQD4(2,3)*F03     )*XMDTXY   
      RD5( 6,1)= RD5( 6,1)+(FQD5(3,1)*AQX2-FQD4(1,3)*AQX2               
     *                     -FQD4(3,3)*F03+FQD3(1,3)*F03      )*XMDTX    
      RD5( 7,1)= RD5( 7,1)+(FQD5(1,1)*X2Y2-FQD4(1,3)*F03*AQX2           
     *                     -FQD4(1,3)*ACY2+FQD3(1,3)*F03     )*XMDTY    
      RD5( 8,1)= RD5( 8,1)+(FQD5(2,1)*X2Y2-FQD4(2,3)*Q2C2               
     *                     +FQD3(2,3)                        )*XMDT     
      RD5( 9,1)= RD5( 9,1)+(FQD5(3,1)*AQX2-FQD4(1,3)*AQX2-FQD4(3,3)     
     *                     +FQD3(1,3)                        )*XMDTY    
      RD5(10,1)= RD5(10,1)+(FQD5(4,1)*AQX2-FQD4(2,3)*F03*AQX2           
     *                     -FQD4(4,3)+FQD3(2,3)*F03          )*XMDT     
      RD5(11,1)= RD5(11,1)+(FQD5(1,1)*ACY4-FQD4(1,3)*F06*ACY2           
     *                     +FQD3(1,3)*F03                    )*XMDTX    
      RD5(12,1)= RD5(12,1)+(FQD5(2,1)*ACY2-FQD4(2,3)*F03     )*XMDTXY   
      RD5(13,1)= RD5(13,1)+(FQD5(3,1)*ACY2-FQD4(1,3)*ACY2-FQD4(3,3)     
     *                     +FQD3(1,3)                        )*XMDTX    
      RD5(14,1)= RD5(14,1)+(FQD5(4,1)-FQD4(2,3)*F03          )*XMDTXY   
      RD5(15,1)= RD5(15,1)+(FQD5(5,1)-FQD4(3,3)*F06                     
     *                     +FQD3(1,3)*F03                    )*XMDTX    
      RD5(16,1)= RD5(16,1)+(FQD5(1,1)*ACY4-FQD4(1,3)*F10*ACY2           
     *                     +FQD3(1,3)*F15                    )*XMDTY    
      RD5(17,1)= RD5(17,1)+(FQD5(2,1)*ACY4-FQD4(2,3)*F06*ACY2           
     *                     +FQD3(2,3)*F03                    )*XMDT     
      RD5(18,1)= RD5(18,1)+(FQD5(3,1)*ACY2-FQD4(1,3)*ACY2               
     *                     -FQD4(3,3)*F03+FQD3(1,3)*F03      )*XMDTY    
      RD5(19,1)= RD5(19,1)+(FQD5(4,1)*ACY2-FQD4(2,3)*F03*ACY2           
     *                     -FQD4(4,3)+FQD3(2,3)*F03          )*XMDT     
      RD5(20,1)= RD5(20,1)+(FQD5(5,1)-FQD4(3,3)*F06                     
     *                     +FQD3(1,3)*F03                    )*XMDTY    
      RD5(21,1)= RD5(21,1)+(FQD5(6,1)-FQD4(4,3)*F10                     
     *                     +FQD3(2,3)*F15                    )*XMDT     
C                                                                       
      FWK( 1,2)=-(FQD5(1,2)*AQX4-FQD4(1,4)*F10*AQX2+FQD3(1,4)*F15)*AQX  
      FWK( 2,2)=-(FQD5(1,2)*AQX4-FQD4(1,4)*F06*AQX2+FQD3(1,4)*F03)*ACY  
      FWK( 3,2)=  FQD5(2,2)*AQX4-FQD4(2,4)*F06*AQX2+FQD3(2,4)*F03       
      FWK( 4,2)=-(FQD5(1,2)*X2Y2-FQD4(1,4)*AQX2-FQD4(1,4)*F03*ACY2      
     *                          +FQD3(1,4)*F03                   )*AQX  
      FWK( 5,2)= (FQD5(2,2)*AQX2-FQD4(2,4)*F03                  )*AQXY  
      FWK( 6,2)=-(FQD5(3,2)*AQX2-FQD4(1,4)*AQX2-FQD4(3,4)*F03           
     *                          +FQD3(1,4)*F03                   )*AQX  
      FWK( 7,2)=-(FQD5(1,2)*X2Y2-FQD4(1,4)*F03*AQX2-FQD4(1,4)*ACY2      
     *                          +FQD3(1,4)*F03                   )*ACY  
      FWK( 8,2)=  FQD5(2,2)*X2Y2-FQD4(2,4)*Q2C2+FQD3(2,4)               
      FWK( 9,2)=-(FQD5(3,2)*AQX2-FQD4(1,4)*AQX2-FQD4(3,4)               
     *                          +FQD3(1,4)                       )*ACY  
      FWK(10,2)=  FQD5(4,2)*AQX2-FQD4(2,4)*F03*AQX2-FQD4(4,4)           
     *                          +FQD3(2,4)*F03                          
      FWK(11,2)=-(FQD5(1,2)*ACY4-FQD4(1,4)*F06*ACY2                     
     *                          +FQD3(1,4)*F03                   )*AQX  
      FWK(12,2)= (FQD5(2,2)*ACY2-FQD4(2,4)*F03                  )*AQXY  
      FWK(13,2)=-(FQD5(3,2)*ACY2-FQD4(1,4)*ACY2-FQD4(3,4)               
     *                          +FQD3(1,4)                       )*AQX  
      FWK(14,2)= (FQD5(4,2)     -FQD4(2,4)*F03                  )*AQXY  
      FWK(15,2)=-(FQD5(5,2)     -FQD4(3,4)*F06+FQD3(1,4)*F03     )*AQX  
      FWK(16,2)=-(FQD5(1,2)*ACY4-FQD4(1,4)*F10*ACY2+FQD3(1,4)*F15)*ACY  
      FWK(17,2)=  FQD5(2,2)*ACY4-FQD4(2,4)*F06*ACY2+FQD3(2,4)*F03       
      FWK(18,2)=-(FQD5(3,2)*ACY2-FQD4(1,4)*ACY2-FQD4(3,4)*F03           
     *                          +FQD3(1,4)*F03                   )*ACY  
      FWK(19,2)=  FQD5(4,2)*ACY2-FQD4(2,4)*F03*ACY2-FQD4(4,4)           
     *                          +FQD3(2,4)*F03                          
      FWK(20,2)=-(FQD5(5,2)     -FQD4(3,4)*F06+FQD3(1,4)*F03     )*ACY  
      FWK(21,2)=  FQD5(6,2)     -FQD4(4,4)*F10+FQD3(2,4)*F15            
      DO 520 I= 2, 3                                                    
         DO 520 J= 1,21                                                 
            RD5(J,I)= RD5(J,I)+FWK(J,2)*WORK(I+12)                      
  520 CONTINUE                                                          
C                                                                       
      CALL FRIKR6( 1, 1,FW6,FQD6, 1,FQD5, 3,FQD4, 3,FQD3)               
C                                                                       
      DO J= 1,28                                                        
         RD6(J)= RD6(J)+FW6(J)*FCC(J,6)                                 
      ENDDO                                                             
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2R   *DECK INTK18                                           
C>                                                                      
C>    @brief   DDPP case                                                
C>                                                                      
C>    @details integration of the DDPP case                             
C>                                                                      
      SUBROUTINE INTK18(IKL)                                            
      use mx_limits, only: mxgsh,mxg2                                   
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
C GENERATE JTYPE=18 INTEGRALS                                           
C                                                                       
C                                                                       
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS,                      
     2                TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2),      
     3                D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2),      
     4                NGANGB                                            
C$omp threadprivate(/GEOMPQ/)
      COMMON /JMSGYH/ SQ(4)                                             
C$omp threadprivate(/JMSGYH/)
      COMMON /FQ08  / FQD(9),FQD0(5),FQD1(2,13),FQD2(3,16),FQD3(4,16),  
     2                FQD4(5,16),FQD5(6,11),FQD6(49),FQD7(24),FQD8( 9)  
C$omp threadprivate(/FQ08/)
      COMMON /KI2 / ACY,ACY2,AQX,AQX2,AQXY,Y03,Y04                      
C$omp threadprivate(/KI2/)
      COMMON /KI4 / RD0(5,5),RD1(3,40),RD2(6,56),RD3(10,52),RD4(15,42), 
     *              RD5(21,24),RD6(336),RD7(144),RD8(45)                
C$omp threadprivate(/KI4/)
C                                                                       
      DIMENSION  WORK(15),FWK(21,10),FW6(28)                            
      DIMENSION  FCU(45,8),FCC(45,8)                                    
C                                                                       
      PARAMETER (ZER=0.0D+00)                                           
      PARAMETER (F03=3.0D+00)                                           
      PARAMETER (F06=6.0D+00)                                           
      PARAMETER (F10=1.0D+01)                                           
      PARAMETER (F15=1.5D+01)                                           
C                                                                       
      IF(IKL.EQ.0) THEN                                                 
         DO I= 1, 5                                                     
            RD0(1,I)= ZER                                               
            RD0(2,I)= ZER                                               
            RD0(3,I)= ZER                                               
            RD0(4,I)= ZER                                               
            RD0(5,I)= ZER                                               
         ENDDO                                                          
         DO I= 1,40                                                     
            RD1(1,I)= ZER                                               
            RD1(2,I)= ZER                                               
            RD1(3,I)= ZER                                               
         ENDDO                                                          
         DO I= 1,41                                                     
            RD2(1,I)= ZER                                               
            RD2(2,I)= ZER                                               
            RD2(3,I)= ZER                                               
            RD2(4,I)= ZER                                               
            RD2(5,I)= ZER                                               
            RD2(6,I)= ZER                                               
         ENDDO                                                          
         DO I= 1,30                                                     
            DO J= 1,10                                                  
               RD3(J,I)= ZER                                            
            ENDDO                                                       
         ENDDO                                                          
         DO I= 1,17                                                     
            DO J= 1,15                                                  
               RD4(J,I)= ZER                                            
            ENDDO                                                       
         ENDDO                                                          
         DO I= 1, 6                                                     
            DO J= 1,21                                                  
               RD5(J,I)= ZER                                            
            ENDDO                                                       
         ENDDO                                                          
         DO J= 1,28                                                     
            RD6(J)= ZER                                                 
         ENDDO                                                          
C                                                                       
         RETURN                                                         
      ENDIF                                                             
C                                                                       
      XMD2= X43 *0.5D+00                                                
      XMD3= XMD2*SQ(3)                                                  
      XMD4= XMD3*XMD2                                                   
      XMD6= XMD4*XMD2                                                   
      XMDT= XMD6*XMD2                                                   
      XMD2= XMD3                                                        
C                                                                       
      XMDTY=-XMDT*ACY                                                   
      XMDTX=-XMDT*AQX                                                   
      XMDTXY=XMDT*AQXY                                                  
C                                                                       
      Y33 = Y03 *Y03                                                    
      Y34 =-Y03 *Y04                                                    
      Y44 = Y04 *Y04                                                    
      WORK( 1)= XMD4                                                    
      WORK( 2)= XMD2*Y33                                                
      WORK( 3)= XMD2*Y34                                                
      WORK( 4)= XMD2*Y44                                                
      WORK( 5)= Y33 *Y44*SQ(3)                                          
C                                                                       
      WORK( 6)=-XMD4*Y03                                                
      WORK( 7)= XMD4*Y04                                                
      WORK( 8)= XMD2*Y33*Y04                                            
      WORK( 9)= XMD2*Y34*Y04                                            
C                                                                       
      WORK(10)= XMD6                                                    
      WORK(11)= XMD4*Y33                                                
      WORK(12)= XMD4*Y34                                                
      WORK(13)= XMD4*Y44                                                
C                                                                       
      WORK(14)=-XMD6*Y03                                                
      WORK(15)= XMD6*Y04                                                
C                                                                       
      CALL FCUFCC(6,XMDT,FCU,FCC)                                       
C                                                                       
      DO 010 I= 1, 5                                                    
         DO 010 J= 1, 5                                                 
            RD0(J,I)= RD0(J,I)+FQD0(I)*WORK(J)                          
  010 CONTINUE                                                          
C                                                                       
      DO I= 1, 8                                                        
         FWK(1,I)=-FQD1(1,I)*AQX                                        
         FWK(2,I)=-FQD1(1,I)*ACY                                        
         FWK(3,I)= FQD1(2,I)                                            
      ENDDO                                                             
      FQD11 = FQD1(1,5)*RAB +FQD1(1,8)                                  
      FQD12 = FQD1(2,5)*RAB +FQD1(2,8)                                  
      FWK(1,9)=-FQD11*AQX                                               
      FWK(2,9)=-FQD11*ACY                                               
      FWK(3,9)= FQD12                                                   
      DO I= 1, 4                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,1)*WORK(I+ 5)                         
         RD1(2,I)= RD1(2,I)+FWK(2,1)*WORK(I+ 5)                         
         RD1(3,I)= RD1(3,I)+FWK(3,1)*WORK(I+ 5)                         
      ENDDO                                                             
      DO I= 5, 8                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,2)*WORK(I+ 1)                         
         RD1(2,I)= RD1(2,I)+FWK(2,2)*WORK(I+ 1)                         
         RD1(3,I)= RD1(3,I)+FWK(3,2)*WORK(I+ 1)                         
      ENDDO                                                             
      DO I= 9,12                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,3)*WORK(I- 3)                         
         RD1(2,I)= RD1(2,I)+FWK(2,3)*WORK(I- 3)                         
         RD1(3,I)= RD1(3,I)+FWK(3,3)*WORK(I- 3)                         
      ENDDO                                                             
      DO I=13,16                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,4)*WORK(I- 7)                         
         RD1(2,I)= RD1(2,I)+FWK(2,4)*WORK(I- 7)                         
         RD1(3,I)= RD1(3,I)+FWK(3,4)*WORK(I- 7)                         
      ENDDO                                                             
      DO I=17,20                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,5)*WORK(I-11)                         
         RD1(2,I)= RD1(2,I)+FWK(2,5)*WORK(I-11)                         
         RD1(3,I)= RD1(3,I)+FWK(3,5)*WORK(I-11)                         
      ENDDO                                                             
      DO I=21,25                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,6)*WORK(I-20)                         
         RD1(2,I)= RD1(2,I)+FWK(2,6)*WORK(I-20)                         
         RD1(3,I)= RD1(3,I)+FWK(3,6)*WORK(I-20)                         
      ENDDO                                                             
      DO I=26,30                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,7)*WORK(I-25)                         
         RD1(2,I)= RD1(2,I)+FWK(2,7)*WORK(I-25)                         
         RD1(3,I)= RD1(3,I)+FWK(3,7)*WORK(I-25)                         
      ENDDO                                                             
      DO I=31,35                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,8)*WORK(I-30)                         
         RD1(2,I)= RD1(2,I)+FWK(2,8)*WORK(I-30)                         
         RD1(3,I)= RD1(3,I)+FWK(3,8)*WORK(I-30)                         
      ENDDO                                                             
      DO I=36,40                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,9)*WORK(I-35)                         
         RD1(2,I)= RD1(2,I)+FWK(2,9)*WORK(I-35)                         
         RD1(3,I)= RD1(3,I)+FWK(3,9)*WORK(I-35)                         
      ENDDO                                                             
C                                                                       
      DO I= 1, 9                                                        
         FWK(1,I)= FQD2(1,I)*AQX2-FQD1(1,I)                             
         FWK(2,I)= FQD2(1,I)*ACY2-FQD1(1,I)                             
         FWK(3,I)= FQD2(3,I)     -FQD1(1,I)                             
         FWK(4,I)= FQD2(1,I)*AQXY                                       
         FWK(5,I)=-FQD2(2,I)*AQX                                        
         FWK(6,I)=-FQD2(2,I)*ACY                                        
      ENDDO                                                             
      FQD21 = FQD2(1,5)*RAB +FQD2(1,8)                                  
      FQD22 = FQD2(2,5)*RAB +FQD2(2,8)                                  
      FQD23 = FQD2(3,5)*RAB +FQD2(3,8)                                  
         FWK(1,10)= FQD21*AQX2-FQD11                                    
         FWK(2,10)= FQD21*ACY2-FQD11                                    
         FWK(3,10)= FQD23     -FQD11                                    
         FWK(4,10)= FQD21*AQXY                                          
         FWK(5,10)=-FQD22*AQX                                           
         FWK(6,10)=-FQD22*ACY                                           
      DO I= 1, 4                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,1)*WORK(I+ 9)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I= 5, 8                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,2)*WORK(I+ 5)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I= 9,12                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,3)*WORK(I+ 1)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=13,16                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,4)*WORK(I- 3)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=17,20                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,5)*WORK(I- 7)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=21,24                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,6)*WORK(I-15)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=25,28                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,7)*WORK(I-19)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=29,32                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,8)*WORK(I-23)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=33,36                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,10)*WORK(I-27)                     
         ENDDO                                                          
      ENDDO                                                             
      DO I=37,41                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,9)*WORK(I-36)                      
         ENDDO                                                          
      ENDDO                                                             
C                                                                       
      DO I= 1, 9                                                        
         FWK( 1,I)=-(FQD3(1,I)*AQX2-FQD2(1,I)*F03)*AQX                  
         FWK( 2,I)=-(FQD3(1,I)*AQX2-FQD2(1,I)    )*ACY                  
         FWK( 3,I)=  FQD3(2,I)*AQX2-FQD2(2,I)                           
         FWK( 4,I)=-(FQD3(1,I)*ACY2-FQD2(1,I)    )*AQX                  
         FWK( 5,I)=  FQD3(2,I)                    *AQXY                 
         FWK( 6,I)=-(FQD3(3,I)     -FQD2(1,I)    )*AQX                  
         FWK( 7,I)=-(FQD3(1,I)*ACY2-FQD2(1,I)*F03)*ACY                  
         FWK( 8,I)=  FQD3(2,I)*ACY2-FQD2(2,I)                           
         FWK( 9,I)=-(FQD3(3,I)     -FQD2(1,I)    )*ACY                  
         FWK(10,I)=  FQD3(4,I)     -FQD2(2,I)*F03                       
      ENDDO                                                             
      FQD31 = FQD3(1,5)*RAB +FQD3(1,8)                                  
      FQD32 = FQD3(2,5)*RAB +FQD3(2,8)                                  
      FQD33 = FQD3(3,5)*RAB +FQD3(3,8)                                  
      FQD34 = FQD3(4,5)*RAB +FQD3(4,8)                                  
      FWK( 1,10)=-(FQD31*AQX2-FQD21*F03)*AQX                            
      FWK( 2,10)=-(FQD31*AQX2-FQD21    )*ACY                            
      FWK( 3,10)=  FQD32*AQX2-FQD22                                     
      FWK( 4,10)=-(FQD31*ACY2-FQD21    )*AQX                            
      FWK( 5,10)=  FQD32                *AQXY                           
      FWK( 6,10)=-(FQD33     -FQD21    )*AQX                            
      FWK( 7,10)=-(FQD31*ACY2-FQD21*F03)*ACY                            
      FWK( 8,10)=  FQD32*ACY2-FQD22                                     
      FWK( 9,10)=-(FQD33     -FQD21    )*ACY                            
      FWK(10,10)=  FQD34     -FQD22*F03                                 
      DO I= 1, 2                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,1)*WORK(I+13)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I= 3, 4                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,2)*WORK(I+11)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I= 5, 6                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,3)*WORK(I+ 9)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I= 7, 8                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,4)*WORK(I+ 7)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I= 9,10                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,5)*WORK(I+ 5)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=11,14                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,6)*WORK(I- 1)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=15,18                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,7)*WORK(I- 5)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=19,22                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,8)*WORK(I- 9)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=23,26                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,10)*WORK(I-13)                     
         ENDDO                                                          
      ENDDO                                                             
      DO I=27,30                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,9)*WORK(I-21)                      
         ENDDO                                                          
      ENDDO                                                             
C                                                                       
      AQX4= AQX2*AQX2                                                   
      ACY4= ACY2*ACY2                                                   
      X2Y2= AQX2*ACY2                                                   
      Q2C2= AQX2+ACY2                                                   
      DO I= 1, 5                                                        
         RD4( 1,I)= RD4( 1,I)+(FQD4(1,I)*AQX4-FQD3(1,I)*F06*AQX2        
     *                                       +FQD2(1,I)*F03  )*XMDT     
         RD4( 2,I)= RD4( 2,I)+(FQD4(1,I)*AQX2-FQD3(1,I)*F03  )*XMDTXY   
         RD4( 3,I)= RD4( 3,I)+(FQD4(2,I)*AQX2-FQD3(2,I)*F03  )*XMDTX    
         RD4( 4,I)= RD4( 4,I)+(FQD4(1,I)*X2Y2-FQD3(1,I)*Q2C2            
     *                                       +FQD2(1,I)      )*XMDT     
         RD4( 5,I)= RD4( 5,I)+(FQD4(2,I)*AQX2-FQD3(2,I)      )*XMDTY    
         RD4( 6,I)= RD4( 6,I)+(FQD4(3,I)*AQX2-FQD3(1,I)*AQX2            
     *                             -FQD3(3,I)+FQD2(1,I)      )*XMDT     
         RD4( 7,I)= RD4( 7,I)+(FQD4(1,I)*ACY2-FQD3(1,I)*F03  )*XMDTXY   
         RD4( 8,I)= RD4( 8,I)+(FQD4(2,I)*ACY2-FQD3(2,I)      )*XMDTX    
         RD4( 9,I)= RD4( 9,I)+(FQD4(3,I)     -FQD3(1,I)      )*XMDTXY   
         RD4(10,I)= RD4(10,I)+(FQD4(4,I)     -FQD3(2,I)*F03  )*XMDTX    
         RD4(11,I)= RD4(11,I)+(FQD4(1,I)*ACY4-FQD3(1,I)*F06*ACY2        
     *                                       +FQD2(1,I)*F03  )*XMDT     
         RD4(12,I)= RD4(12,I)+(FQD4(2,I)*ACY2-FQD3(2,I)*F03  )*XMDTY    
         RD4(13,I)= RD4(13,I)+(FQD4(3,I)*ACY2-FQD3(1,I)*ACY2            
     *                             -FQD3(3,I)+FQD2(1,I)      )*XMDT     
         RD4(14,I)= RD4(14,I)+(FQD4(4,I)     -FQD3(2,I)*F03  )*XMDTY    
         RD4(15,I)= RD4(15,I)+(FQD4(5,I)-FQD3(3,I)*F06                  
     *                                       +FQD2(1,I)*F03  )*XMDT     
      ENDDO                                                             
      FQD2(1,5)= FQD21                                                  
      FQD3(1,5)= FQD31                                                  
      FQD3(2,5)= FQD32                                                  
      FQD3(3,5)= FQD33                                                  
      DO J= 1, 5                                                        
         FQD4(J,5)= FQD4(J,5)*RAB +FQD4(J,8)                            
      ENDDO                                                             
C                                                                       
      DO I= 5, 9                                                        
         FWK( 1,I)=  FQD4(1,I)*AQX4-FQD3(1,I)*F06*AQX2+FQD2(1,I)*F03    
         FWK( 2,I)= (FQD4(1,I)*AQX2-FQD3(1,I)*F03               )*AQXY  
         FWK( 3,I)=-(FQD4(2,I)*AQX2-FQD3(2,I)*F03               )*AQX   
         FWK( 4,I)=  FQD4(1,I)*X2Y2-FQD3(1,I)*Q2C2    +FQD2(1,I)        
         FWK( 5,I)=-(FQD4(2,I)*AQX2-FQD3(2,I)                   )*ACY   
         FWK( 6,I)=  FQD4(3,I)*AQX2-FQD3(1,I)*AQX2-FQD3(3,I)+FQD2(1,I)  
         FWK( 7,I)= (FQD4(1,I)*ACY2-FQD3(1,I)*F03               )*AQXY  
         FWK( 8,I)=-(FQD4(2,I)*ACY2-FQD3(2,I)                   )*AQX   
         FWK( 9,I)= (FQD4(3,I)     -FQD3(1,I)                   )*AQXY  
         FWK(10,I)=-(FQD4(4,I)     -FQD3(2,I)*F03               )*AQX   
         FWK(11,I)=  FQD4(1,I)*ACY4-FQD3(1,I)*F06*ACY2+FQD2(1,I)*F03    
         FWK(12,I)=-(FQD4(2,I)*ACY2-FQD3(2,I)*F03               )*ACY   
         FWK(13,I)=  FQD4(3,I)*ACY2-FQD3(1,I)*ACY2-FQD3(3,I)+FQD2(1,I)  
         FWK(14,I)=-(FQD4(4,I)     -FQD3(2,I)*F03               )*ACY   
         FWK(15,I)=  FQD4(5,I)     -FQD3(3,I)*F06+FQD2(1,I)*F03         
      ENDDO                                                             
      DO I= 6, 7                                                        
         DO J= 1,15                                                     
            RD4(J,I)= RD4(J,I)+FWK(J,5)*WORK(I+ 8)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I= 8, 9                                                        
         DO J= 1,15                                                     
            RD4(J,I)= RD4(J,I)+FWK(J,6)*WORK(I+ 6)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=10,11                                                        
         DO J= 1,15                                                     
            RD4(J,I)= RD4(J,I)+FWK(J,7)*WORK(I+ 4)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=12,13                                                        
         DO J= 1,15                                                     
            RD4(J,I)= RD4(J,I)+FWK(J,8)*WORK(I+ 2)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=14,17                                                        
         DO J= 1,15                                                     
            RD4(J,I)= RD4(J,I)+FWK(J,9)*WORK(I- 4)                      
         ENDDO                                                          
      ENDDO                                                             
C                                                                       
      DO J= 1, 6                                                        
         FQD5(J,1)= FQD5(J,1)*RAB +FQD5(J,4)                            
      ENDDO                                                             
      DO I= 1, 4                                                        
        RD5( 1,I)= RD5( 1,I)+(FQD5(1,I)*AQX4-FQD4(1,I+4)*F10*AQX2       
     *                       +FQD3(1,I+4)*F15                  )*XMDTX  
        RD5( 2,I)= RD5( 2,I)+(FQD5(1,I)*AQX4-FQD4(1,I+4)*F06*AQX2       
     *                       +FQD3(1,I+4)*F03                  )*XMDTY  
        RD5( 3,I)= RD5( 3,I)+(FQD5(2,I)*AQX4-FQD4(2,I+4)*F06*AQX2       
     *                       +FQD3(2,I+4)*F03                  )*XMDT   
        RD5( 4,I)= RD5( 4,I)+(FQD5(1,I)*X2Y2-FQD4(1,I+4)*F03*ACY2       
     *                       -FQD4(1,I+4)*AQX2+FQD3(1,I+4)*F03 )*XMDTX  
        RD5( 5,I)= RD5( 5,I)+(FQD5(2,I)*AQX2-FQD4(2,I+4)*F03   )*XMDTXY 
        RD5( 6,I)= RD5( 6,I)+(FQD5(3,I)*AQX2-FQD4(1,I+4)*AQX2           
     *                       -FQD4(3,I+4)*F03+FQD3(1,I+4)*F03  )*XMDTX  
        RD5( 7,I)= RD5( 7,I)+(FQD5(1,I)*X2Y2-FQD4(1,I+4)*F03*AQX2       
     *                       -FQD4(1,I+4)*ACY2+FQD3(1,I+4)*F03 )*XMDTY  
        RD5( 8,I)= RD5( 8,I)+(FQD5(2,I)*X2Y2-FQD4(2,I+4)*Q2C2           
     *                       +FQD3(2,I+4)                      )*XMDT   
        RD5( 9,I)= RD5( 9,I)+(FQD5(3,I)*AQX2-FQD4(1,I+4)*AQX2           
     *                       -FQD4(3,I+4)+FQD3(1,I+4)          )*XMDTY  
        RD5(10,I)= RD5(10,I)+(FQD5(4,I)*AQX2-FQD4(2,I+4)*F03*AQX2       
     *                       -FQD4(4,I+4)+FQD3(2,I+4)*F03      )*XMDT   
        RD5(11,I)= RD5(11,I)+(FQD5(1,I)*ACY4-FQD4(1,I+4)*F06*ACY2       
     *                       +FQD3(1,I+4)*F03                  )*XMDTX  
        RD5(12,I)= RD5(12,I)+(FQD5(2,I)*ACY2-FQD4(2,I+4)*F03   )*XMDTXY 
        RD5(13,I)= RD5(13,I)+(FQD5(3,I)*ACY2-FQD4(1,I+4)*ACY2           
     *                       -FQD4(3,I+4)+FQD3(1,I+4)          )*XMDTX  
        RD5(14,I)= RD5(14,I)+(FQD5(4,I)-FQD4(2,I+4)*F03        )*XMDTXY 
        RD5(15,I)= RD5(15,I)+(FQD5(5,I)-FQD4(3,I+4)*F06                 
     *                       +FQD3(1,I+4)*F03                  )*XMDTX  
        RD5(16,I)= RD5(16,I)+(FQD5(1,I)*ACY4-FQD4(1,I+4)*F10*ACY2       
     *                       +FQD3(1,I+4)*F15                  )*XMDTY  
        RD5(17,I)= RD5(17,I)+(FQD5(2,I)*ACY4-FQD4(2,I+4)*F06*ACY2       
     *                       +FQD3(2,I+4)*F03                  )*XMDT   
        RD5(18,I)= RD5(18,I)+(FQD5(3,I)*ACY2-FQD4(1,I+4)*ACY2           
     *                       -FQD4(3,I+4)*F03+FQD3(1,I+4)*F03  )*XMDTY  
        RD5(19,I)= RD5(19,I)+(FQD5(4,I)*ACY2-FQD4(2,I+4)*F03*ACY2       
     *                       -FQD4(4,I+4)+FQD3(2,I+4)*F03      )*XMDT   
        RD5(20,I)= RD5(20,I)+(FQD5(5,I)-FQD4(3,I+4)*F06                 
     *                       +FQD3(1,I+4)*F03                  )*XMDTY  
        RD5(21,I)= RD5(21,I)+(FQD5(6,I)-FQD4(4,I+4)*F10                 
     *                       +FQD3(2,I+4)*F15                  )*XMDT   
      ENDDO                                                             
      FWK( 1,5)=-(FQD5(1,5)*AQX4-FQD4(1,9)*F10*AQX2+FQD3(1,9)*F15)*AQX  
      FWK( 2,5)=-(FQD5(1,5)*AQX4-FQD4(1,9)*F06*AQX2+FQD3(1,9)*F03)*ACY  
      FWK( 3,5)=  FQD5(2,5)*AQX4-FQD4(2,9)*F06*AQX2+FQD3(2,9)*F03       
      FWK( 4,5)=-(FQD5(1,5)*X2Y2-FQD4(1,9)*AQX2-FQD4(1,9)*F03*ACY2      
     *                          +FQD3(1,9)*F03                   )*AQX  
      FWK( 5,5)= (FQD5(2,5)*AQX2-FQD4(2,9)*F03                   )*AQXY 
      FWK( 6,5)=-(FQD5(3,5)*AQX2-FQD4(1,9)*AQX2-FQD4(3,9)*F03           
     *                          +FQD3(1,9)*F03                   )*AQX  
      FWK( 7,5)=-(FQD5(1,5)*X2Y2-FQD4(1,9)*F03*AQX2-FQD4(1,9)*ACY2      
     *                          +FQD3(1,9)*F03                   )*ACY  
      FWK( 8,5)=  FQD5(2,5)*X2Y2-FQD4(2,9)*Q2C2+FQD3(2,9)               
      FWK( 9,5)=-(FQD5(3,5)*AQX2-FQD4(1,9)*AQX2-FQD4(3,9)               
     *                          +FQD3(1,9)                       )*ACY  
      FWK(10,5)=  FQD5(4,5)*AQX2-FQD4(2,9)*F03*AQX2-FQD4(4,9)           
     *                          +FQD3(2,9)*F03                          
      FWK(11,5)=-(FQD5(1,5)*ACY4-FQD4(1,9)*F06*ACY2+FQD3(1,9)*F03)*AQX  
      FWK(12,5)= (FQD5(2,5)*ACY2-FQD4(2,9)*F03                   )*AQXY 
      FWK(13,5)=-(FQD5(3,5)*ACY2-FQD4(1,9)*ACY2-FQD4(3,9)               
     *                          +FQD3(1,9)                       )*AQX  
      FWK(14,5)= (FQD5(4,5)     -FQD4(2,9)*F03                   )*AQXY 
      FWK(15,5)=-(FQD5(5,5)     -FQD4(3,9)*F06+FQD3(1,9)*F03     )*AQX  
      FWK(16,5)=-(FQD5(1,5)*ACY4-FQD4(1,9)*F10*ACY2+FQD3(1,9)*F15)*ACY  
      FWK(17,5)=  FQD5(2,5)*ACY4-FQD4(2,9)*F06*ACY2+FQD3(2,9)*F03       
      FWK(18,5)=-(FQD5(3,5)*ACY2-FQD4(1,9)*ACY2-FQD4(3,9)*F03           
     *                          +FQD3(1,9)*F03                   )*ACY  
      FWK(19,5)=  FQD5(4,5)*ACY2-FQD4(2,9)*F03*ACY2-FQD4(4,9)           
     *                          +FQD3(2,9)*F03                          
      FWK(20,5)=-(FQD5(5,5)     -FQD4(3,9)*F06+FQD3(1,9)*F03     )*ACY  
      FWK(21,5)=  FQD5(6,5)     -FQD4(4,9)*F10+FQD3(2,9)*F15            
      DO 520 I= 5, 6                                                    
         DO 520 J= 1,21                                                 
            RD5(J,I)= RD5(J,I)+FWK(J,5)*WORK(I+ 9)                      
  520 CONTINUE                                                          
C                                                                       
      CALL FRIKR6( 1, 1,FW6,FQD6, 4,FQD5, 8,FQD4, 8,FQD3)               
C                                                                       
      DO J= 1,28                                                        
         RD6(J)= RD6(J)+FW6(J)*FCC(J,6)                                 
      ENDDO                                                             
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2R   *DECK INTK19                                           
C>                                                                      
C>    @brief   DPDP case                                                
C>                                                                      
C>    @details integration of the DPDP case                             
C>                                                                      
      SUBROUTINE INTK19(IKL)                                            
      use mx_limits, only: mxgsh,mxg2                                   
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
C GENERATE JTYPE=19 INTEGRALS                                           
C                                                                       
C                                                                       
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS,                      
     2                TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2),      
     3                D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2),      
     4                NGANGB                                            
C$omp threadprivate(/GEOMPQ/)
      COMMON /JMSGYH/ SQ(4)                                             
C$omp threadprivate(/JMSGYH/)
      COMMON /FQ08  / FQD(9),FQD0(5),FQD1(2,13),FQD2(3,16),FQD3(4,16),  
     2                FQD4(5,16),FQD5(6,11),FQD6(49),FQD7(24),FQD8( 9)  
C$omp threadprivate(/FQ08/)
      COMMON /KI2 / ACY,ACY2,AQX,AQX2,AQXY,Y03,Y04                      
C$omp threadprivate(/KI2/)
      COMMON /KI4 / RD0(5,5),RD1(3,40),RD2(6,56),RD3(10,52),RD4(15,42), 
     *              RD5(21,24),RD6(336),RD7(144),RD8(45)                
C$omp threadprivate(/KI4/)
C                                                                       
      DIMENSION  WORK(12),FWK(21,12),FW6(28)                            
      DIMENSION  FCU(45,8),FCC(45,8)                                    
C                                                                       
      PARAMETER (ZER=0.0D+00)                                           
      PARAMETER (F03=3.0D+00)                                           
      PARAMETER (F06=6.0D+00)                                           
      PARAMETER (F10=1.0D+01)                                           
      PARAMETER (F15=1.5D+01)                                           
C                                                                       
      IF(IKL.EQ.0) THEN                                                 
         DO I= 1, 5                                                     
            RD0(1,I)= ZER                                               
            RD0(2,I)= ZER                                               
            RD0(3,I)= ZER                                               
            RD0(4,I)= ZER                                               
            RD0(5,I)= ZER                                               
         ENDDO                                                          
         DO I= 1,40                                                     
            RD1(1,I)= ZER                                               
            RD1(2,I)= ZER                                               
            RD1(3,I)= ZER                                               
         ENDDO                                                          
         DO I= 1,46                                                     
            RD2(1,I)= ZER                                               
            RD2(2,I)= ZER                                               
            RD2(3,I)= ZER                                               
            RD2(4,I)= ZER                                               
            RD2(5,I)= ZER                                               
            RD2(6,I)= ZER                                               
         ENDDO                                                          
         DO I= 1,34                                                     
            DO J= 1,10                                                  
               RD3(J,I)= ZER                                            
            ENDDO                                                       
         ENDDO                                                          
         DO I= 1,17                                                     
            DO J= 1,15                                                  
               RD4(J,I)= ZER                                            
            ENDDO                                                       
         ENDDO                                                          
         DO I= 1, 6                                                     
            DO J= 1,21                                                  
               RD5(J,I)= ZER                                            
            ENDDO                                                       
         ENDDO                                                          
         DO J= 1,28                                                     
            RD6(J)= ZER                                                 
         ENDDO                                                          
C                                                                       
         RETURN                                                         
      ENDIF                                                             
C                                                                       
      XMD2= X43 *0.5D+00                                                
      XMD1= XMD2*SQ(1)                                                  
      XMD3= XMD2*SQ(2)                                                  
      XMD4= XMD2*XMD2                                                   
      XMD5= XMD4*SQ(2)                                                  
      XMDT= XMD4*XMD3                                                   
C                                                                       
      XMDTY=-XMDT*ACY                                                   
      XMDTX=-XMDT*AQX                                                   
      XMDTXY=XMDT*AQXY                                                  
C                                                                       
      Y33 = Y03 *Y03                                                    
      Y34 =-Y03 *Y04                                                    
      WORK( 1)= XMD1                                                    
      WORK( 2)= Y33 *SQ(1)                                              
      WORK( 3)=-XMD3*Y03                                                
      WORK( 4)= XMD3*Y04                                                
      WORK( 5)= Y33 *SQ(2)*Y04                                          
C                                                                       
      WORK( 6)=-XMD1*Y03                                                
      WORK( 7)= XMD5                                                    
      WORK( 8)= XMD3*Y33                                                
      WORK( 9)= XMD3*Y34                                                
C                                                                       
      WORK(10)= XMD4*SQ(1)                                              
      WORK(11)=-XMD5*Y03                                                
      WORK(12)= XMD5*Y04                                                
C                                                                       
      CALL FCUFCC(6,XMDT,FCU,FCC)                                       
C                                                                       
      DO 010 I= 1, 5                                                    
         DO 010 J= 1, 5                                                 
            RD0(J,I)= RD0(J,I)+FQD0(I)*WORK(J)                          
  010 CONTINUE                                                          
C                                                                       
      DO I= 1, 8                                                        
         FWK(1,I)=-FQD1(1,I)*AQX                                        
         FWK(2,I)=-FQD1(1,I)*ACY                                        
         FWK(3,I)= FQD1(2,I)                                            
      ENDDO                                                             
      FQD11 = FQD1(1,5)*RAB +FQD1(1, 7)                                 
      FQD12 = FQD1(2,5)*RAB +FQD1(2, 7)                                 
         FWK(1,9)=-FQD11*AQX                                            
         FWK(2,9)=-FQD11*ACY                                            
         FWK(3,9)= FQD12                                                
      DO I= 1, 4                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,1)*WORK(I+ 5)                         
         RD1(2,I)= RD1(2,I)+FWK(2,1)*WORK(I+ 5)                         
         RD1(3,I)= RD1(3,I)+FWK(3,1)*WORK(I+ 5)                         
      ENDDO                                                             
      DO I= 5, 8                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,2)*WORK(I+ 1)                         
         RD1(2,I)= RD1(2,I)+FWK(2,2)*WORK(I+ 1)                         
         RD1(3,I)= RD1(3,I)+FWK(3,2)*WORK(I+ 1)                         
      ENDDO                                                             
      DO I= 9,12                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,3)*WORK(I- 3)                         
         RD1(2,I)= RD1(2,I)+FWK(2,3)*WORK(I- 3)                         
         RD1(3,I)= RD1(3,I)+FWK(3,3)*WORK(I- 3)                         
      ENDDO                                                             
      DO I=13,16                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,4)*WORK(I- 7)                         
         RD1(2,I)= RD1(2,I)+FWK(2,4)*WORK(I- 7)                         
         RD1(3,I)= RD1(3,I)+FWK(3,4)*WORK(I- 7)                         
      ENDDO                                                             
      DO I=17,20                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,5)*WORK(I-11)                         
         RD1(2,I)= RD1(2,I)+FWK(2,5)*WORK(I-11)                         
         RD1(3,I)= RD1(3,I)+FWK(3,5)*WORK(I-11)                         
      ENDDO                                                             
      DO I=21,25                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,6)*WORK(I-20)                         
         RD1(2,I)= RD1(2,I)+FWK(2,6)*WORK(I-20)                         
         RD1(3,I)= RD1(3,I)+FWK(3,6)*WORK(I-20)                         
      ENDDO                                                             
      DO I=26,30                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,7)*WORK(I-25)                         
         RD1(2,I)= RD1(2,I)+FWK(2,7)*WORK(I-25)                         
         RD1(3,I)= RD1(3,I)+FWK(3,7)*WORK(I-25)                         
      ENDDO                                                             
      DO I=31,35                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,8)*WORK(I-30)                         
         RD1(2,I)= RD1(2,I)+FWK(2,8)*WORK(I-30)                         
         RD1(3,I)= RD1(3,I)+FWK(3,8)*WORK(I-30)                         
      ENDDO                                                             
      DO I=36,40                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,9)*WORK(I-35)                         
         RD1(2,I)= RD1(2,I)+FWK(2,9)*WORK(I-35)                         
         RD1(3,I)= RD1(3,I)+FWK(3,9)*WORK(I-35)                         
      ENDDO                                                             
C                                                                       
      DO I= 1,10                                                        
         FWK(1,I)= FQD2(1,I)*AQX2-FQD1(1,I)                             
         FWK(2,I)= FQD2(1,I)*ACY2-FQD1(1,I)                             
         FWK(3,I)= FQD2(3,I)     -FQD1(1,I)                             
         FWK(4,I)= FQD2(1,I)*AQXY                                       
         FWK(5,I)=-FQD2(2,I)*AQX                                        
         FWK(6,I)=-FQD2(2,I)*ACY                                        
      ENDDO                                                             
      FQD21 = FQD2(1,5)*RAB +FQD2(1, 7)                                 
      FQD22 = FQD2(2,5)*RAB +FQD2(2, 7)                                 
      FQD23 = FQD2(3,5)*RAB +FQD2(3, 7)                                 
         FWK(1,11)= FQD21*AQX2-FQD11                                    
         FWK(2,11)= FQD21*ACY2-FQD11                                    
         FWK(3,11)= FQD23     -FQD11                                    
         FWK(4,11)= FQD21*AQXY                                          
         FWK(5,11)=-FQD22*AQX                                           
         FWK(6,11)=-FQD22*ACY                                           
      FQD13 = FQD1(1,8)*RAB +FQD1(1,10)                                 
      FQD24 = FQD2(1,8)*RAB +FQD2(1,10)                                 
      FQD25 = FQD2(2,8)*RAB +FQD2(2,10)                                 
      FQD26 = FQD2(3,8)*RAB +FQD2(3,10)                                 
         FWK(1,12)= FQD24*AQX2-FQD13                                    
         FWK(2,12)= FQD24*ACY2-FQD13                                    
         FWK(3,12)= FQD26     -FQD13                                    
         FWK(4,12)= FQD24*AQXY                                          
         FWK(5,12)=-FQD25*AQX                                           
         FWK(6,12)=-FQD25*ACY                                           
      DO I= 1, 3                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,1)*WORK(I+ 9)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I= 4, 6                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,2)*WORK(I+ 6)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I= 7, 9                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,3)*WORK(I+ 3)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=10,12                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,4)*WORK(I)                         
         ENDDO                                                          
      ENDDO                                                             
      DO I=13,15                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,5)*WORK(I- 3)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=16,19                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,6)*WORK(I-10)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=20,23                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,7)*WORK(I-14)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=24,27                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,8)*WORK(I-18)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=28,31                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,11)*WORK(I-22)                     
         ENDDO                                                          
      ENDDO                                                             
      DO I=32,36                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,9)*WORK(I-31)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=37,41                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,10)*WORK(I-36)                     
         ENDDO                                                          
      ENDDO                                                             
      DO I=42,46                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,12)*WORK(I-41)                     
         ENDDO                                                          
      ENDDO                                                             
C                                                                       
      DO I= 1, 5                                                        
         RD3( 1,I)= RD3( 1,I)+(FQD3(1,I)*AQX2-FQD2(1,I)*F03)*XMDTX      
         RD3( 2,I)= RD3( 2,I)+(FQD3(1,I)*AQX2-FQD2(1,I)    )*XMDTY      
         RD3( 3,I)= RD3( 3,I)+(FQD3(2,I)*AQX2-FQD2(2,I)    )*XMDT       
         RD3( 4,I)= RD3( 4,I)+(FQD3(1,I)*ACY2-FQD2(1,I)    )*XMDTX      
         RD3( 5,I)= RD3( 5,I)+ FQD3(2,I)                    *XMDTXY     
         RD3( 6,I)= RD3( 6,I)+(FQD3(3,I)     -FQD2(1,I)    )*XMDTX      
         RD3( 7,I)= RD3( 7,I)+(FQD3(1,I)*ACY2-FQD2(1,I)*F03)*XMDTY      
         RD3( 8,I)= RD3( 8,I)+(FQD3(2,I)*ACY2-FQD2(2,I)    )*XMDT       
         RD3( 9,I)= RD3( 9,I)+(FQD3(3,I)     -FQD2(1,I)    )*XMDTY      
         RD3(10,I)= RD3(10,I)+(FQD3(4,I)     -FQD2(2,I)*F03)*XMDT       
      ENDDO                                                             
      FQD2(1,5)= FQD21                                                  
      FQD2(2,5)= FQD22                                                  
      FQD3(1,5)= FQD3(1,5)*RAB +FQD3(1, 7)                              
      FQD3(2,5)= FQD3(2,5)*RAB +FQD3(2, 7)                              
      FQD3(3,5)= FQD3(3,5)*RAB +FQD3(3, 7)                              
      FQD3(4,5)= FQD3(4,5)*RAB +FQD3(4, 7)                              
      DO I= 5,11                                                        
         FWK( 1,I)=-(FQD3(1,I)*AQX2-FQD2(1,I)*F03)*AQX                  
         FWK( 2,I)=-(FQD3(1,I)*AQX2-FQD2(1,I)    )*ACY                  
         FWK( 3,I)=  FQD3(2,I)*AQX2-FQD2(2,I)                           
         FWK( 4,I)=-(FQD3(1,I)*ACY2-FQD2(1,I)    )*AQX                  
         FWK( 5,I)=  FQD3(2,I)                    *AQXY                 
         FWK( 6,I)=-(FQD3(3,I)     -FQD2(1,I)    )*AQX                  
         FWK( 7,I)=-(FQD3(1,I)*ACY2-FQD2(1,I)*F03)*ACY                  
         FWK( 8,I)=  FQD3(2,I)*ACY2-FQD2(2,I)                           
         FWK( 9,I)=-(FQD3(3,I)     -FQD2(1,I)    )*ACY                  
         FWK(10,I)=  FQD3(4,I)     -FQD2(2,I)*F03                       
      ENDDO                                                             
      FQD31 = FQD3(1,8)*RAB +FQD3(1,10)                                 
      FQD32 = FQD3(2,8)*RAB +FQD3(2,10)                                 
      FQD33 = FQD3(3,8)*RAB +FQD3(3,10)                                 
      FQD34 = FQD3(4,8)*RAB +FQD3(4,10)                                 
         FWK( 1,12)=-(FQD31*AQX2-FQD24*F03)*AQX                         
         FWK( 2,12)=-(FQD31*AQX2-FQD24    )*ACY                         
         FWK( 3,12)=  FQD32*AQX2-FQD25                                  
         FWK( 4,12)=-(FQD31*ACY2-FQD24    )*AQX                         
         FWK( 5,12)=  FQD32                *AQXY                        
         FWK( 6,12)=-(FQD33     -FQD24    )*AQX                         
         FWK( 7,12)=-(FQD31*ACY2-FQD24*F03)*ACY                         
         FWK( 8,12)=  FQD32*ACY2-FQD25                                  
         FWK( 9,12)=-(FQD33     -FQD24    )*ACY                         
         FWK(10,12)=  FQD34     -FQD25*F03                              
      DO I= 6, 8                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,5)*WORK(I+ 4)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I= 9,11                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,6)*WORK(I+ 1)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=12,14                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,7)*WORK(I- 2)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=15,17                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,8)*WORK(I- 5)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=18,21                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,9)*WORK(I-12)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=22,25                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,10)*WORK(I-16)                     
         ENDDO                                                          
      ENDDO                                                             
      DO I=26,29                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,12)*WORK(I-20)                     
         ENDDO                                                          
      ENDDO                                                             
      DO I=30,34                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,11)*WORK(I-29)                     
         ENDDO                                                          
      ENDDO                                                             
C                                                                       
      DO J= 1, 5                                                        
         FQD4(J,1)= FQD4(J,1)*RAB +FQD4(J,3)                            
      ENDDO                                                             
      AQX4= AQX2*AQX2                                                   
      ACY4= ACY2*ACY2                                                   
      X2Y2= AQX2*ACY2                                                   
      Q2C2= AQX2+ACY2                                                   
      DO I= 1, 4                                                        
         RD4( 1,I)= RD4( 1,I)+(FQD4(1,I)*AQX4-FQD3(1,I+4)*F06*AQX2      
     *                                       +FQD2(1,I+4)*F03  )*XMDT   
         RD4( 2,I)= RD4( 2,I)+(FQD4(1,I)*AQX2-FQD3(1,I+4)*F03  )*XMDTXY 
         RD4( 3,I)= RD4( 3,I)+(FQD4(2,I)*AQX2-FQD3(2,I+4)*F03  )*XMDTX  
         RD4( 4,I)= RD4( 4,I)+(FQD4(1,I)*X2Y2-FQD3(1,I+4)*Q2C2          
     *                                       +FQD2(1,I+4)      )*XMDT   
         RD4( 5,I)= RD4( 5,I)+(FQD4(2,I)*AQX2-FQD3(2,I+4)      )*XMDTY  
         RD4( 6,I)= RD4( 6,I)+(FQD4(3,I)*AQX2-FQD3(1,I+4)*AQX2          
     *                           -FQD3(3,I+4)+FQD2(1,I+4)      )*XMDT   
         RD4( 7,I)= RD4( 7,I)+(FQD4(1,I)*ACY2-FQD3(1,I+4)*F03  )*XMDTXY 
         RD4( 8,I)= RD4( 8,I)+(FQD4(2,I)*ACY2-FQD3(2,I+4)      )*XMDTX  
         RD4( 9,I)= RD4( 9,I)+(FQD4(3,I)     -FQD3(1,I+4)      )*XMDTXY 
         RD4(10,I)= RD4(10,I)+(FQD4(4,I)     -FQD3(2,I+4)*F03  )*XMDTX  
         RD4(11,I)= RD4(11,I)+(FQD4(1,I)*ACY4-FQD3(1,I+4)*F06*ACY2      
     *                                       +FQD2(1,I+4)*F03  )*XMDT   
         RD4(12,I)= RD4(12,I)+(FQD4(2,I)*ACY2-FQD3(2,I+4)*F03  )*XMDTY  
         RD4(13,I)= RD4(13,I)+(FQD4(3,I)*ACY2-FQD3(1,I+4)*ACY2          
     *                           -FQD3(3,I+4)+FQD2(1,I+4)      )*XMDT   
         RD4(14,I)= RD4(14,I)+(FQD4(4,I)     -FQD3(2,I+4)*F03  )*XMDTY  
         RD4(15,I)= RD4(15,I)+(FQD4(5,I)-FQD3(3,I+4)*F06                
     *                                       +FQD2(1,I+4)*F03  )*XMDT   
      ENDDO                                                             
      FQD2(1,8)= FQD24                                                  
      FQD3(1,8)= FQD31                                                  
      FQD3(2,8)= FQD32                                                  
      FQD3(3,8)= FQD33                                                  
      DO J= 1, 5                                                        
         FQD4(J,4)= FQD4(J,4)*RAB +FQD4(J,6)                            
      ENDDO                                                             
      DO I= 4, 7                                                        
         FWK( 1,I)=  FQD4(1,I)*AQX4-FQD3(1,I+4)*F06*AQX2                
     *                                              +FQD2(1,I+4)*F03    
         FWK( 2,I)= (FQD4(1,I)*AQX2-FQD3(1,I+4)*F03             )*AQXY  
         FWK( 3,I)=-(FQD4(2,I)*AQX2-FQD3(2,I+4)*F03             )*AQX   
         FWK( 4,I)=  FQD4(1,I)*X2Y2-FQD3(1,I+4)*Q2C2+FQD2(1,I+4)        
         FWK( 5,I)=-(FQD4(2,I)*AQX2-FQD3(2,I+4)                 )*ACY   
         FWK( 6,I)=  FQD4(3,I)*AQX2-FQD3(1,I+4)*AQX2-FQD3(3,I+4)        
     *                                              +FQD2(1,I+4)        
         FWK( 7,I)= (FQD4(1,I)*ACY2-FQD3(1,I+4)*F03             )*AQXY  
         FWK( 8,I)=-(FQD4(2,I)*ACY2-FQD3(2,I+4)                 )*AQX   
         FWK( 9,I)= (FQD4(3,I)     -FQD3(1,I+4)                 )*AQXY  
         FWK(10,I)=-(FQD4(4,I)     -FQD3(2,I+4)*F03             )*AQX   
         FWK(11,I)=  FQD4(1,I)*ACY4-FQD3(1,I+4)*F06*ACY2                
     *                                              +FQD2(1,I+4)*F03    
         FWK(12,I)=-(FQD4(2,I)*ACY2-FQD3(2,I+4)*F03             )*ACY   
         FWK(13,I)=  FQD4(3,I)*ACY2-FQD3(1,I+4)*ACY2-FQD3(3,I+4)        
     *                                              +FQD2(1,I+4)        
         FWK(14,I)=-(FQD4(4,I)     -FQD3(2,I+4)*F03             )*ACY   
         FWK(15,I)=  FQD4(5,I)     -FQD3(3,I+4)*F06 +FQD2(1,I+4)*F03    
      ENDDO                                                             
      DO I= 5, 7                                                        
         DO J= 1,15                                                     
            RD4(J,I)= RD4(J,I)+FWK(J,5)*WORK(I+ 5)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I= 8,10                                                        
         DO J= 1,15                                                     
            RD4(J,I)= RD4(J,I)+FWK(J,6)*WORK(I+ 2)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=11,13                                                        
         DO J= 1,15                                                     
            RD4(J,I)= RD4(J,I)+FWK(J,4)*WORK(I- 1)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=14,17                                                        
         DO J= 1,15                                                     
            RD4(J,I)= RD4(J,I)+FWK(J,7)*WORK(I- 8)                      
         ENDDO                                                          
      ENDDO                                                             
C                                                                       
      DO J= 1, 6                                                        
         FQD5(J,1)= FQD5(J,1)*RAB +FQD5(J,3)                            
      ENDDO                                                             
      DO I= 1, 3                                                        
         RD5( 1,I)= RD5( 1,I)+(FQD5(1,I)*AQX4-FQD4(1,I+3)*F10*AQX2      
     *                        +FQD3(1,I+7)*F15                  )*XMDTX 
         RD5( 2,I)= RD5( 2,I)+(FQD5(1,I)*AQX4-FQD4(1,I+3)*F06*AQX2      
     *                        +FQD3(1,I+7)*F03                  )*XMDTY 
         RD5( 3,I)= RD5( 3,I)+(FQD5(2,I)*AQX4-FQD4(2,I+3)*F06*AQX2      
     *                        +FQD3(2,I+7)*F03                  )*XMDT  
         RD5( 4,I)= RD5( 4,I)+(FQD5(1,I)*X2Y2-FQD4(1,I+3)*F03*ACY2      
     *                        -FQD4(1,I+3)*AQX2+FQD3(1,I+7)*F03 )*XMDTX 
         RD5( 5,I)= RD5( 5,I)+(FQD5(2,I)*AQX2-FQD4(2,I+3)*F03   )*XMDTXY
         RD5( 6,I)= RD5( 6,I)+(FQD5(3,I)*AQX2-FQD4(1,I+3)*AQX2          
     *                        -FQD4(3,I+3)*F03+FQD3(1,I+7)*F03  )*XMDTX 
         RD5( 7,I)= RD5( 7,I)+(FQD5(1,I)*X2Y2-FQD4(1,I+3)*F03*AQX2      
     *                        -FQD4(1,I+3)*ACY2+FQD3(1,I+7)*F03 )*XMDTY 
         RD5( 8,I)= RD5( 8,I)+(FQD5(2,I)*X2Y2-FQD4(2,I+3)*Q2C2          
     *                        +FQD3(2,I+7)                      )*XMDT  
         RD5( 9,I)= RD5( 9,I)+(FQD5(3,I)*AQX2-FQD4(1,I+3)*AQX2          
     *                        -FQD4(3,I+3)+FQD3(1,I+7)          )*XMDTY 
         RD5(10,I)= RD5(10,I)+(FQD5(4,I)*AQX2-FQD4(2,I+3)*F03*AQX2      
     *                        -FQD4(4,I+3)+FQD3(2,I+7)*F03      )*XMDT  
         RD5(11,I)= RD5(11,I)+(FQD5(1,I)*ACY4-FQD4(1,I+3)*F06*ACY2      
     *                        +FQD3(1,I+7)*F03                  )*XMDTX 
         RD5(12,I)= RD5(12,I)+(FQD5(2,I)*ACY2-FQD4(2,I+3)*F03   )*XMDTXY
         RD5(13,I)= RD5(13,I)+(FQD5(3,I)*ACY2-FQD4(1,I+3)*ACY2          
     *                        -FQD4(3,I+3)+FQD3(1,I+7)          )*XMDTX 
         RD5(14,I)= RD5(14,I)+(FQD5(4,I)-FQD4(2,I+3)*F03        )*XMDTXY
         RD5(15,I)= RD5(15,I)+(FQD5(5,I)-FQD4(3,I+3)*F06                
     *                        +FQD3(1,I+7)*F03                  )*XMDTX 
         RD5(16,I)= RD5(16,I)+(FQD5(1,I)*ACY4-FQD4(1,I+3)*F10*ACY2      
     *                        +FQD3(1,I+7)*F15                  )*XMDTY 
         RD5(17,I)= RD5(17,I)+(FQD5(2,I)*ACY4-FQD4(2,I+3)*F06*ACY2      
     *                        +FQD3(2,I+7)*F03                  )*XMDT  
         RD5(18,I)= RD5(18,I)+(FQD5(3,I)*ACY2-FQD4(1,I+3)*ACY2          
     *                        -FQD4(3,I+3)*F03+FQD3(1,I+7)*F03  )*XMDTY 
         RD5(19,I)= RD5(19,I)+(FQD5(4,I)*ACY2-FQD4(2,I+3)*F03*ACY2      
     *                        -FQD4(4,I+3)+FQD3(2,I+7)*F03      )*XMDT  
         RD5(20,I)= RD5(20,I)+(FQD5(5,I)-FQD4(3,I+3)*F06                
     *                        +FQD3(1,I+7)*F03                  )*XMDTY 
         RD5(21,I)= RD5(21,I)+(FQD5(6,I)-FQD4(4,I+3)*F10                
     *                        +FQD3(2,I+7)*F15                  )*XMDT  
      ENDDO                                                             
      FWK( 1,4)=-(FQD5(1,4)*AQX4-FQD4(1,7)*F10*AQX2+FQD3(1,11)*F15)*AQX 
      FWK( 2,4)=-(FQD5(1,4)*AQX4-FQD4(1,7)*F06*AQX2+FQD3(1,11)*F03)*ACY 
      FWK( 3,4)=  FQD5(2,4)*AQX4-FQD4(2,7)*F06*AQX2+FQD3(2,11)*F03      
      FWK( 4,4)=-(FQD5(1,4)*X2Y2-FQD4(1,7)*AQX2-FQD4(1,7)*F03*ACY2      
     *                                         +FQD3(1,11)*F03    )*AQX 
      FWK( 5,4)= (FQD5(2,4)*AQX2-FQD4(2,7)*F03                    )*AQXY
      FWK( 6,4)=-(FQD5(3,4)*AQX2-FQD4(1,7)*AQX2-FQD4(3,7)*F03           
     *            +FQD3(1,11)*F03                                 )*AQX 
      FWK( 7,4)=-(FQD5(1,4)*X2Y2-FQD4(1,7)*F03*AQX2-FQD4(1,7)*ACY2      
     *            +FQD3(1,11)*F03                                 )*ACY 
      FWK( 8,4)=  FQD5(2,4)*X2Y2-FQD4(2,7)*Q2C2+FQD3(2,11)              
      FWK( 9,4)=-(FQD5(3,4)*AQX2-FQD4(1,7)*AQX2-FQD4(3,7)               
     *                                             +FQD3(1,11)    )*ACY 
      FWK(10,4)=  FQD5(4,4)*AQX2-FQD4(2,7)*F03*AQX2-FQD4(4,7)           
     *                                             +FQD3(2,11)*F03      
      FWK(11,4)=-(FQD5(1,4)*ACY4-FQD4(1,7)*F06*ACY2                     
     *                                             +FQD3(1,11)*F03)*AQX 
      FWK(12,4)= (FQD5(2,4)*ACY2-FQD4(2,7)*F03                    )*AQXY
      FWK(13,4)=-(FQD5(3,4)*ACY2-FQD4(1,7)*ACY2-FQD4(3,7)               
     *                                             +FQD3(1,11)    )*AQX 
      FWK(14,4)= (FQD5(4,4)     -FQD4(2,7)*F03                    )*AQXY
      FWK(15,4)=-(FQD5(5,4)     -FQD4(3,7)*F06+FQD3(1,11)*F03     )*AQX 
      FWK(16,4)=-(FQD5(1,4)*ACY4-FQD4(1,7)*F10*ACY2                     
     *                                             +FQD3(1,11)*F15)*ACY 
      FWK(17,4)=  FQD5(2,4)*ACY4-FQD4(2,7)*F06*ACY2+FQD3(2,11)*F03      
      FWK(18,4)=-(FQD5(3,4)*ACY2-FQD4(1,7)*ACY2-FQD4(3,7)*F03           
     *                                             +FQD3(1,11)*F03)*ACY 
      FWK(19,4)=  FQD5(4,4)*ACY2-FQD4(2,7)*F03*ACY2-FQD4(4,7)           
     *                                             +FQD3(2,11)*F03      
      FWK(20,4)=-(FQD5(5,4)     -FQD4(3,7)*F06+FQD3(1,11)*F03     )*ACY 
      FWK(21,4)=  FQD5(6,4)     -FQD4(4,7)*F10+FQD3(2,11)*F15           
      DO 520 I= 4, 6                                                    
         DO 520 J= 1,21                                                 
            RD5(J,I)= RD5(J,I)+FWK(J,4)*WORK(I+ 6)                      
  520 CONTINUE                                                          
C                                                                       
      CALL FRIKR6( 1, 1,FW6,FQD6, 3,FQD5, 6,FQD4,10,FQD3)               
C                                                                       
      DO J= 1,28                                                        
         RD6(J)= RD6(J)+FW6(J)*FCC(J,6)                                 
      ENDDO                                                             
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2R   *DECK INTK20                                           
C>                                                                      
C>    @brief   DDDP case                                                
C>                                                                      
C>    @details integration of the DDDP case                             
C>                                                                      
      SUBROUTINE INTK20(IKL)                                            
      use mx_limits, only: mxgsh,mxg2                                   
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
C GENERATE JTYPE=20 INTEGRALS                                           
C                                                                       
C                                                                       
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS,                      
     2                TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2),      
     3                D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2),      
     4                NGANGB                                            
C$omp threadprivate(/GEOMPQ/)
      COMMON /JMSGYH/ SQ(4)                                             
C$omp threadprivate(/JMSGYH/)
      COMMON /FQ08  / FQD(9),FQD0(5),FQD1(2,13),FQD2(3,16),FQD3(4,16),  
     2                FQD4(5,16),FQD5(6,11),FQD6(7,7),FQD7(24),FQD8( 9) 
C$omp threadprivate(/FQ08/)
      COMMON /KI2 / ACY,ACY2,AQX,AQX2,AQXY,Y03,Y04                      
C$omp threadprivate(/KI2/)
      COMMON /KI4 / RD0(5,5),RD1(3,40),RD2(6,56),RD3(10,52),RD4(15,42), 
     *              RD5(21,24),RD6(28,12),RD7(144),RD8(45)              
C$omp threadprivate(/KI4/)
C                                                                       
      DIMENSION  WORK(15),FWK(28,13),FW6(28, 4),FW7(36)                 
      DIMENSION  FCU(45,8),FCC(45,8)                                    
C                                                                       
      PARAMETER (ZER=0.0D+00)                                           
      PARAMETER (F03=3.0D+00)                                           
      PARAMETER (F06=6.0D+00)                                           
      PARAMETER (F10=1.0D+01)                                           
      PARAMETER (F15=1.5D+01)                                           
C                                                                       
      IF(IKL.EQ.0) THEN                                                 
         DO I= 1, 5                                                     
            RD0(1,I)= ZER                                               
            RD0(2,I)= ZER                                               
            RD0(3,I)= ZER                                               
            RD0(4,I)= ZER                                               
            RD0(5,I)= ZER                                               
         ENDDO                                                          
         DO I= 1,40                                                     
            RD1(1,I)= ZER                                               
            RD1(2,I)= ZER                                               
            RD1(3,I)= ZER                                               
         ENDDO                                                          
         DO I= 1,51                                                     
            RD2(1,I)= ZER                                               
            RD2(2,I)= ZER                                               
            RD2(3,I)= ZER                                               
            RD2(4,I)= ZER                                               
            RD2(5,I)= ZER                                               
            RD2(6,I)= ZER                                               
         ENDDO                                                          
         DO I= 1,43                                                     
            DO J= 1,10                                                  
               RD3(J,I)= ZER                                            
            ENDDO                                                       
         ENDDO                                                          
         DO I= 1,29                                                     
            DO J= 1,15                                                  
               RD4(J,I)= ZER                                            
            ENDDO                                                       
         ENDDO                                                          
         DO I= 1,14                                                     
            DO J= 1,21                                                  
               RD5(J,I)= ZER                                            
            ENDDO                                                       
         ENDDO                                                          
         DO I= 1, 5                                                     
            DO J= 1,28                                                  
               RD6(J,I)= ZER                                            
            ENDDO                                                       
         ENDDO                                                          
         DO J= 1,36                                                     
            RD7(J)= ZER                                                 
         ENDDO                                                          
C                                                                       
         RETURN                                                         
      ENDIF                                                             
C                                                                       
      XMD2= X43 *0.5D+00                                                
      XMD3= XMD2*SQ(3)                                                  
      XMD4= XMD3*XMD2                                                   
      XMD6= XMD4*XMD2                                                   
      XMDT= XMD6*XMD2                                                   
      XMD2= XMD3                                                        
C                                                                       
      XMDTY=-XMDT*ACY                                                   
      XMDTX=-XMDT*AQX                                                   
      XMDTXY=XMDT*AQXY                                                  
C                                                                       
      Y33 = Y03 *Y03                                                    
      Y34 =-Y03 *Y04                                                    
      Y44 = Y04 *Y04                                                    
      WORK( 1)= XMD4                                                    
      WORK( 2)= XMD2*Y33                                                
      WORK( 3)= XMD2*Y34                                                
      WORK( 4)= XMD2*Y44                                                
      WORK( 5)= Y33 *Y44*SQ(3)                                          
C                                                                       
      WORK( 6)=-XMD4*Y03                                                
      WORK( 7)= XMD4*Y04                                                
      WORK( 8)= XMD2*Y33*Y04                                            
      WORK( 9)= XMD2*Y34*Y04                                            
C                                                                       
      WORK(10)= XMD6                                                    
      WORK(11)= XMD4*Y33                                                
      WORK(12)= XMD4*Y34                                                
      WORK(13)= XMD4*Y44                                                
C                                                                       
      WORK(14)=-XMD6*Y03                                                
      WORK(15)= XMD6*Y04                                                
C                                                                       
      CALL FCUFCC(7,XMDT,FCU,FCC)                                       
C                                                                       
      DO 010 I= 1, 5                                                    
         DO 010 J= 1, 5                                                 
            RD0(J,I)= RD0(J,I)+FQD0(I)*WORK(J)                          
  010 CONTINUE                                                          
C                                                                       
      DO I= 1, 8                                                        
         FWK(1,I)=-FQD1(1,I)*AQX                                        
         FWK(2,I)=-FQD1(1,I)*ACY                                        
         FWK(3,I)= FQD1(2,I)                                            
      ENDDO                                                             
      FQD1(1,12)= FQD1(1, 8)*RAB +FQD1(1,10)                            
      FQD1(1,13)= FQD1(1, 5)*RAB +FQD1(1, 7)                            
      FQD1(2,13)= FQD1(2, 5)*RAB +FQD1(2, 7)                            
         FWK(1,9)=-FQD1(1,13)*AQX                                       
         FWK(2,9)=-FQD1(1,13)*ACY                                       
         FWK(3,9)= FQD1(2,13)                                           
      DO I= 1, 4                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,1)*WORK(I+ 5)                         
         RD1(2,I)= RD1(2,I)+FWK(2,1)*WORK(I+ 5)                         
         RD1(3,I)= RD1(3,I)+FWK(3,1)*WORK(I+ 5)                         
      ENDDO                                                             
      DO I= 5, 8                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,2)*WORK(I+ 1)                         
         RD1(2,I)= RD1(2,I)+FWK(2,2)*WORK(I+ 1)                         
         RD1(3,I)= RD1(3,I)+FWK(3,2)*WORK(I+ 1)                         
      ENDDO                                                             
      DO I= 9,12                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,3)*WORK(I- 3)                         
         RD1(2,I)= RD1(2,I)+FWK(2,3)*WORK(I- 3)                         
         RD1(3,I)= RD1(3,I)+FWK(3,3)*WORK(I- 3)                         
      ENDDO                                                             
      DO I=13,16                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,4)*WORK(I- 7)                         
         RD1(2,I)= RD1(2,I)+FWK(2,4)*WORK(I- 7)                         
         RD1(3,I)= RD1(3,I)+FWK(3,4)*WORK(I- 7)                         
      ENDDO                                                             
      DO I=17,20                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,5)*WORK(I-11)                         
         RD1(2,I)= RD1(2,I)+FWK(2,5)*WORK(I-11)                         
         RD1(3,I)= RD1(3,I)+FWK(3,5)*WORK(I-11)                         
      ENDDO                                                             
      DO I=21,25                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,6)*WORK(I-20)                         
         RD1(2,I)= RD1(2,I)+FWK(2,6)*WORK(I-20)                         
         RD1(3,I)= RD1(3,I)+FWK(3,6)*WORK(I-20)                         
      ENDDO                                                             
      DO I=26,30                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,7)*WORK(I-25)                         
         RD1(2,I)= RD1(2,I)+FWK(2,7)*WORK(I-25)                         
         RD1(3,I)= RD1(3,I)+FWK(3,7)*WORK(I-25)                         
      ENDDO                                                             
      DO I=31,35                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,8)*WORK(I-30)                         
         RD1(2,I)= RD1(2,I)+FWK(2,8)*WORK(I-30)                         
         RD1(3,I)= RD1(3,I)+FWK(3,8)*WORK(I-30)                         
      ENDDO                                                             
      DO I=36,40                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,9)*WORK(I-35)                         
         RD1(2,I)= RD1(2,I)+FWK(2,9)*WORK(I-35)                         
         RD1(3,I)= RD1(3,I)+FWK(3,9)*WORK(I-35)                         
      ENDDO                                                             
C                                                                       
      DO I= 1,10                                                        
         FWK(1,I)= FQD2(1,I)*AQX2-FQD1(1,I)                             
         FWK(2,I)= FQD2(1,I)*ACY2-FQD1(1,I)                             
         FWK(3,I)= FQD2(3,I)     -FQD1(1,I)                             
         FWK(4,I)= FQD2(1,I)*AQXY                                       
         FWK(5,I)=-FQD2(2,I)*AQX                                        
         FWK(6,I)=-FQD2(2,I)*ACY                                        
      ENDDO                                                             
      DO I= 1, 3                                                        
         FQD2(I,12)= FQD2(I,8)*RAB +FQD2(I,10)                          
         FQD2(I,13)= FQD2(I,5)*RAB +FQD2(I, 7)                          
      ENDDO                                                             
      DO I=12,13                                                        
         FWK(1,I)= FQD2(1,I)*AQX2-FQD1(1,I)                             
         FWK(2,I)= FQD2(1,I)*ACY2-FQD1(1,I)                             
         FWK(3,I)= FQD2(3,I)     -FQD1(1,I)                             
         FWK(4,I)= FQD2(1,I)*AQXY                                       
         FWK(5,I)=-FQD2(2,I)*AQX                                        
         FWK(6,I)=-FQD2(2,I)*ACY                                        
      ENDDO                                                             
      DO I= 1, 4                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,1)*WORK(I+ 9)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I= 5, 8                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,2)*WORK(I+ 5)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I= 9,12                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,3)*WORK(I+ 1)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=13,16                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,4)*WORK(I- 3)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=17,20                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,5)*WORK(I- 7)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=21,24                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,6)*WORK(I-15)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=25,28                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,7)*WORK(I-19)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=29,32                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,8)*WORK(I-23)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=33,36                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,13)*WORK(I-27)                     
         ENDDO                                                          
      ENDDO                                                             
      DO I=37,41                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,9)*WORK(I-36)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=42,46                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,10)*WORK(I-41)                     
         ENDDO                                                          
      ENDDO                                                             
      DO I=47,51                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,12)*WORK(I-46)                     
         ENDDO                                                          
      ENDDO                                                             
C                                                                       
      DO J= 1, 4                                                        
         FQD3(J,12)= FQD3(J,8)*RAB +FQD3(J,10)                          
         FQD3(J,13)= FQD3(J,5)*RAB +FQD3(J, 7)                          
      ENDDO                                                             
      DO I= 1,13                                                        
         FWK( 1,I)=-(FQD3(1,I)*AQX2-FQD2(1,I)*F03)*AQX                  
         FWK( 2,I)=-(FQD3(1,I)*AQX2-FQD2(1,I)    )*ACY                  
         FWK( 3,I)=  FQD3(2,I)*AQX2-FQD2(2,I)                           
         FWK( 4,I)=-(FQD3(1,I)*ACY2-FQD2(1,I)    )*AQX                  
         FWK( 5,I)=  FQD3(2,I)                    *AQXY                 
         FWK( 6,I)=-(FQD3(3,I)     -FQD2(1,I)    )*AQX                  
         FWK( 7,I)=-(FQD3(1,I)*ACY2-FQD2(1,I)*F03)*ACY                  
         FWK( 8,I)=  FQD3(2,I)*ACY2-FQD2(2,I)                           
         FWK( 9,I)=-(FQD3(3,I)     -FQD2(1,I)    )*ACY                  
         FWK(10,I)=  FQD3(4,I)     -FQD2(2,I)*F03                       
      ENDDO                                                             
      DO I= 1, 2                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,1)*WORK(I+13)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I= 3, 4                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,2)*WORK(I+11)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I= 5, 6                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,3)*WORK(I+ 9)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I= 7, 8                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,4)*WORK(I+ 7)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I= 9,10                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,5)*WORK(I+ 5)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=11,14                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,6)*WORK(I- 1)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=15,18                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,7)*WORK(I- 5)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=19,22                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,8)*WORK(I- 9)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=23,26                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,13)*WORK(I-13)                     
         ENDDO                                                          
      ENDDO                                                             
      DO I=27,30                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,9)*WORK(I-21)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=31,34                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,10)*WORK(I-25)                     
         ENDDO                                                          
      ENDDO                                                             
      DO I=35,38                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,12)*WORK(I-29)                     
         ENDDO                                                          
      ENDDO                                                             
      DO I=39,43                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,11)*WORK(I-38)                     
         ENDDO                                                          
      ENDDO                                                             
C                                                                       
      AQX4= AQX2*AQX2                                                   
      ACY4= ACY2*ACY2                                                   
      X2Y2= AQX2*ACY2                                                   
      Q2C2= AQX2+ACY2                                                   
      DO I= 1, 5                                                        
         RD4( 1,I)= RD4( 1,I)+(FQD4(1,I)*AQX4-FQD3(1,I)*F06*AQX2        
     *                                       +FQD2(1,I)*F03   )*XMDT    
         RD4( 2,I)= RD4( 2,I)+(FQD4(1,I)*AQX2-FQD3(1,I)*F03   )*XMDTXY  
         RD4( 3,I)= RD4( 3,I)+(FQD4(2,I)*AQX2-FQD3(2,I)*F03   )*XMDTX   
         RD4( 4,I)= RD4( 4,I)+(FQD4(1,I)*X2Y2-FQD3(1,I)*Q2C2            
     *                                       +FQD2(1,I)       )*XMDT    
         RD4( 5,I)= RD4( 5,I)+(FQD4(2,I)*AQX2-FQD3(2,I)       )*XMDTY   
         RD4( 6,I)= RD4( 6,I)+(FQD4(3,I)*AQX2-FQD3(1,I)*AQX2            
     *                                    -FQD3(3,I)+FQD2(1,I))*XMDT    
         RD4( 7,I)= RD4( 7,I)+(FQD4(1,I)*ACY2-FQD3(1,I)*F03   )*XMDTXY  
         RD4( 8,I)= RD4( 8,I)+(FQD4(2,I)*ACY2-FQD3(2,I)       )*XMDTX   
         RD4( 9,I)= RD4( 9,I)+(FQD4(3,I)     -FQD3(1,I)       )*XMDTXY  
         RD4(10,I)= RD4(10,I)+(FQD4(4,I)     -FQD3(2,I)*F03   )*XMDTX   
         RD4(11,I)= RD4(11,I)+(FQD4(1,I)*ACY4-FQD3(1,I)*F06*ACY2        
     *                                       +FQD2(1,I)*F03   )*XMDT    
         RD4(12,I)= RD4(12,I)+(FQD4(2,I)*ACY2-FQD3(2,I)*F03   )*XMDTY   
         RD4(13,I)= RD4(13,I)+(FQD4(3,I)*ACY2-FQD3(1,I)*ACY2            
     *                                    -FQD3(3,I)+FQD2(1,I))*XMDT    
         RD4(14,I)= RD4(14,I)+(FQD4(4,I)     -FQD3(2,I)*F03   )*XMDTY   
         RD4(15,I)= RD4(15,I)+(FQD4(5,I)     -FQD3(3,I)*F06             
     *                                       +FQD2(1,I)*F03   )*XMDT    
      ENDDO                                                             
      FQD2(1,5)= FQD2(1,13)                                             
      FQD3(1,5)= FQD3(1,13)                                             
      FQD3(2,5)= FQD3(2,13)                                             
      FQD3(3,5)= FQD3(3,13)                                             
      DO J= 1, 5                                                        
         FQD4(J, 5)= FQD4(J,5)*RAB +FQD4(J, 7)                          
         FQD4(J,12)= FQD4(J,8)*RAB +FQD4(J,10)                          
      ENDDO                                                             
      DO I= 5,12                                                        
         FWK( 1,I)=  FQD4(1,I)*AQX4-FQD3(1,I)*F06*AQX2+FQD2(1,I)*F03    
         FWK( 2,I)= (FQD4(1,I)*AQX2-FQD3(1,I)*F03                )*AQXY 
         FWK( 3,I)=-(FQD4(2,I)*AQX2-FQD3(2,I)*F03                )*AQX  
         FWK( 4,I)=  FQD4(1,I)*X2Y2-FQD3(1,I)*Q2C2+FQD2(1,I)            
         FWK( 5,I)=-(FQD4(2,I)*AQX2-FQD3(2,I)                    )*ACY  
         FWK( 6,I)=  FQD4(3,I)*AQX2-FQD3(1,I)*AQX2-FQD3(3,I)+FQD2(1,I)  
         FWK( 7,I)= (FQD4(1,I)*ACY2-FQD3(1,I)*F03                )*AQXY 
         FWK( 8,I)=-(FQD4(2,I)*ACY2-FQD3(2,I)                    )*AQX  
         FWK( 9,I)= (FQD4(3,I)     -FQD3(1,I)                    )*AQXY 
         FWK(10,I)=-(FQD4(4,I)     -FQD3(2,I)*F03                )*AQX  
         FWK(11,I)=  FQD4(1,I)*ACY4-FQD3(1,I)*F06*ACY2+FQD2(1,I)*F03    
         FWK(12,I)=-(FQD4(2,I)*ACY2-FQD3(2,I)*F03                )*ACY  
         FWK(13,I)=  FQD4(3,I)*ACY2-FQD3(1,I)*ACY2-FQD3(3,I)+FQD2(1,I)  
         FWK(14,I)=-(FQD4(4,I)     -FQD3(2,I)*F03                )*ACY  
         FWK(15,I)=  FQD4(5,I)     -FQD3(3,I)*F06 +FQD2(1,I)*F03        
      ENDDO                                                             
      DO I= 6, 7                                                        
         DO J= 1,15                                                     
            RD4(J,I)= RD4(J,I)+FWK(J,5)*WORK(I+ 8)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I= 8, 9                                                        
         DO J= 1,15                                                     
            RD4(J,I)= RD4(J,I)+FWK(J,6)*WORK(I+ 6)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=10,11                                                        
         DO J= 1,15                                                     
            RD4(J,I)= RD4(J,I)+FWK(J,7)*WORK(I+ 4)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=12,13                                                        
         DO J= 1,15                                                     
            RD4(J,I)= RD4(J,I)+FWK(J,8)*WORK(I+ 2)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=14,17                                                        
         DO J= 1,15                                                     
            RD4(J,I)= RD4(J,I)+FWK(J,9)*WORK(I- 4)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=18,21                                                        
         DO J= 1,15                                                     
            RD4(J,I)= RD4(J,I)+FWK(J,10)*WORK(I- 8)                     
         ENDDO                                                          
      ENDDO                                                             
      DO I=22,25                                                        
         DO J= 1,15                                                     
            RD4(J,I)= RD4(J,I)+FWK(J,12)*WORK(I-12)                     
         ENDDO                                                          
      ENDDO                                                             
      DO I=26,29                                                        
         DO J= 1,15                                                     
            RD4(J,I)= RD4(J,I)+FWK(J,11)*WORK(I-20)                     
         ENDDO                                                          
      ENDDO                                                             
C                                                                       
      DO J= 1, 6                                                        
         FQD5(J,1)= FQD5(J,1)*RAB +FQD5(J,3)                            
      ENDDO                                                             
      DO I= 1, 4                                                        
         RD5( 1,I)= RD5( 1,I)+(FQD5(1,I)*AQX4-FQD4(1,I+4)*F10*AQX2      
     *                                       +FQD3(1,I+4)*F15  )*XMDTX  
         RD5( 2,I)= RD5( 2,I)+(FQD5(1,I)*AQX4-FQD4(1,I+4)*F06*AQX2      
     *                                       +FQD3(1,I+4)*F03  )*XMDTY  
         RD5( 3,I)= RD5( 3,I)+(FQD5(2,I)*AQX4-FQD4(2,I+4)*F06*AQX2      
     *                                       +FQD3(2,I+4)*F03  )*XMDT   
         RD5( 4,I)= RD5( 4,I)+(FQD5(1,I)*X2Y2-FQD4(1,I+4)*F03*ACY2      
     *                        -FQD4(1,I+4)*AQX2+FQD3(1,I+4)*F03)*XMDTX  
         RD5( 5,I)= RD5( 5,I)+(FQD5(2,I)*AQX2-FQD4(2,I+4)*F03  )*XMDTXY 
         RD5( 6,I)= RD5( 6,I)+(FQD5(3,I)*AQX2-FQD4(1,I+4)*AQX2          
     *                        -FQD4(3,I+4)*F03+FQD3(1,I+4)*F03 )*XMDTX  
         RD5( 7,I)= RD5( 7,I)+(FQD5(1,I)*X2Y2-FQD4(1,I+4)*F03*AQX2      
     *                        -FQD4(1,I+4)*ACY2+FQD3(1,I+4)*F03)*XMDTY  
         RD5( 8,I)= RD5( 8,I)+(FQD5(2,I)*X2Y2-FQD4(2,I+4)*Q2C2          
     *                                       +FQD3(2,I+4)      )*XMDT   
         RD5( 9,I)= RD5( 9,I)+(FQD5(3,I)*AQX2-FQD4(1,I+4)*AQX2          
     *                        -FQD4(3,I+4)+FQD3(1,I+4)         )*XMDTY  
         RD5(10,I)= RD5(10,I)+(FQD5(4,I)*AQX2-FQD4(2,I+4)*F03*AQX2      
     *                        -FQD4(4,I+4)+FQD3(2,I+4)*F03     )*XMDT   
         RD5(11,I)= RD5(11,I)+(FQD5(1,I)*ACY4-FQD4(1,I+4)*F06*ACY2      
     *                                       +FQD3(1,I+4)*F03  )*XMDTX  
         RD5(12,I)= RD5(12,I)+(FQD5(2,I)*ACY2-FQD4(2,I+4)*F03  )*XMDTXY 
         RD5(13,I)= RD5(13,I)+(FQD5(3,I)*ACY2-FQD4(1,I+4)*ACY2          
     *                        -FQD4(3,I+4)+FQD3(1,I+4)         )*XMDTX  
         RD5(14,I)= RD5(14,I)+(FQD5(4,I)-FQD4(2,I+4)*F03       )*XMDTXY 
         RD5(15,I)= RD5(15,I)+(FQD5(5,I)-FQD4(3,I+4)*F06                
     *                                       +FQD3(1,I+4)*F03  )*XMDTX  
         RD5(16,I)= RD5(16,I)+(FQD5(1,I)*ACY4-FQD4(1,I+4)*F10*ACY2      
     *                                       +FQD3(1,I+4)*F15  )*XMDTY  
         RD5(17,I)= RD5(17,I)+(FQD5(2,I)*ACY4-FQD4(2,I+4)*F06*ACY2      
     *                                       +FQD3(2,I+4)*F03  )*XMDT   
         RD5(18,I)= RD5(18,I)+(FQD5(3,I)*ACY2-FQD4(1,I+4)*ACY2          
     *                        -FQD4(3,I+4)*F03+FQD3(1,I+4)*F03 )*XMDTY  
         RD5(19,I)= RD5(19,I)+(FQD5(4,I)*ACY2-FQD4(2,I+4)*F03*ACY2      
     *                        -FQD4(4,I+4)+FQD3(2,I+4)*F03     )*XMDT   
         RD5(20,I)= RD5(20,I)+(FQD5(5,I)-FQD4(3,I+4)*F06                
     *                                  +FQD3(1,I+4)*F03       )*XMDTY  
         RD5(21,I)= RD5(21,I)+(FQD5(6,I)-FQD4(4,I+4)*F10                
     *                                  +FQD3(2,I+4)*F15       )*XMDT   
      ENDDO                                                             
      FQD3(1,8)= FQD3(1,12)                                             
      FQD3(2,8)= FQD3(2,12)                                             
      FQD4(1,8)= FQD4(1,12)                                             
      FQD4(2,8)= FQD4(2,12)                                             
      FQD4(3,8)= FQD4(3,12)                                             
      FQD4(4,8)= FQD4(4,12)                                             
      DO J= 1, 6                                                        
         FQD5(J,4)= FQD5(J,4)*RAB +FQD5(J,6)                            
      ENDDO                                                             
      DO I= 4, 7                                                        
         FWK( 1,I)=-(FQD5(1,I)*AQX4-FQD4(1,I+4)*F10*AQX2                
     *                             +FQD3(1,I+4)*F15              )*AQX  
         FWK( 2,I)=-(FQD5(1,I)*AQX4-FQD4(1,I+4)*F06*AQX2                
     *                             +FQD3(1,I+4)*F03              )*ACY  
         FWK( 3,I)=  FQD5(2,I)*AQX4-FQD4(2,I+4)*F06*AQX2                
     *                             +FQD3(2,I+4)*F03                     
         FWK( 4,I)=-(FQD5(1,I)*X2Y2-FQD4(1,I+4)*AQX2                    
     *              -FQD4(1,I+4)*F03*ACY2+FQD3(1,I+4)*F03        )*AQX  
         FWK( 5,I)= (FQD5(2,I)*AQX2-FQD4(2,I+4)*F03              )*AQXY 
         FWK( 6,I)=-(FQD5(3,I)*AQX2-FQD4(1,I+4)*AQX2-FQD4(3,I+4)*F03    
     *                             +FQD3(1,I+4)*F03              )*AQX  
         FWK( 7,I)=-(FQD5(1,I)*X2Y2-FQD4(1,I+4)*F03*AQX2                
     *              -FQD4(1,I+4)*ACY2+FQD3(1,I+4)*F03            )*ACY  
         FWK( 8,I)=  FQD5(2,I)*X2Y2-FQD4(2,I+4)*Q2C2+FQD3(2,I+4)        
         FWK( 9,I)=-(FQD5(3,I)*AQX2-FQD4(1,I+4)*AQX2-FQD4(3,I+4)        
     *                             +FQD3(1,I+4)                  )*ACY  
         FWK(10,I)=  FQD5(4,I)*AQX2-FQD4(2,I+4)*F03*AQX2-FQD4(4,I+4)    
     *                             +FQD3(2,I+4)*F03                     
         FWK(11,I)=-(FQD5(1,I)*ACY4-FQD4(1,I+4)*F06*ACY2                
     *                             +FQD3(1,I+4)*F03              )*AQX  
         FWK(12,I)= (FQD5(2,I)*ACY2-FQD4(2,I+4)*F03              )*AQXY 
         FWK(13,I)=-(FQD5(3,I)*ACY2-FQD4(1,I+4)*ACY2-FQD4(3,I+4)        
     *                             +FQD3(1,I+4)                  )*AQX  
         FWK(14,I)= (FQD5(4,I)     -FQD4(2,I+4)*F03              )*AQXY 
         FWK(15,I)=-(FQD5(5,I)     -FQD4(3,I+4)*F06                     
     *                             +FQD3(1,I+4)*F03              )*AQX  
         FWK(16,I)=-(FQD5(1,I)*ACY4-FQD4(1,I+4)*F10*ACY2                
     *                             +FQD3(1,I+4)*F15              )*ACY  
         FWK(17,I)=  FQD5(2,I)*ACY4-FQD4(2,I+4)*F06*ACY2                
     *                             +FQD3(2,I+4)*F03                     
         FWK(18,I)=-(FQD5(3,I)*ACY2-FQD4(1,I+4)*ACY2-FQD4(3,I+4)*F03    
     *                             +FQD3(1,I+4)*F03              )*ACY  
         FWK(19,I)=  FQD5(4,I)*ACY2-FQD4(2,I+4)*F03*ACY2-FQD4(4,I+4)    
     *                             +FQD3(2,I+4)*F03                     
         FWK(20,I)=-(FQD5(5,I)-FQD4(3,I+4)*F06+FQD3(1,I+4)*F03   )*ACY  
         FWK(21,I)=  FQD5(6,I)-FQD4(4,I+4)*F10+FQD3(2,I+4)*F15          
      ENDDO                                                             
      DO I= 5, 6                                                        
         DO J= 1,21                                                     
            RD5(J,I)= RD5(J,I)+FWK(J,5)*WORK(I+ 9)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I= 7, 8                                                        
         DO J= 1,21                                                     
            RD5(J,I)= RD5(J,I)+FWK(J,6)*WORK(I+ 7)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I= 9,10                                                        
         DO J= 1,21                                                     
            RD5(J,I)= RD5(J,I)+FWK(J,4)*WORK(I+ 5)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=11,14                                                        
         DO J= 1,21                                                     
            RD5(J,I)= RD5(J,I)+FWK(J,7)*WORK(I- 1)                      
         ENDDO                                                          
      ENDDO                                                             
C                                                                       
      DO J= 1, 7                                                        
         FQD6(J,1)=FQD6(J,1)*RAB +FQD6(J,3)                             
      ENDDO                                                             
C                                                                       
      CALL FRIKR6( 1, 4,FW6,FQD6, 3,FQD5, 7,FQD4, 7,FQD3)               
C                                                                       
      DO 620 I= 1, 3                                                    
         DO 620 J= 1,28                                                 
            RD6(J,I)= RD6(J,I)+FW6(J,I)*FCC(J,6)                        
  620 CONTINUE                                                          
      I= 4                                                              
C        FW6( 1,I)= FW6( 1,I)*FCU( 1,6)                                 
         FW6( 2,I)= FW6( 2,I)*FCU( 2,6)                                 
         FW6( 3,I)= FW6( 3,I)*FCU( 3,6)                                 
C        FW6( 4,I)= FW6( 4,I)*FCU( 4,6)                                 
         FW6( 5,I)= FW6( 5,I)*FCU( 5,6)                                 
C        FW6( 6,I)= FW6( 6,I)*FCU( 6,6)                                 
         FW6( 7,I)= FW6( 7,I)*FCU( 7,6)                                 
         FW6( 8,I)= FW6( 8,I)*FCU( 8,6)                                 
         FW6( 9,I)= FW6( 9,I)*FCU( 9,6)                                 
         FW6(10,I)= FW6(10,I)*FCU(10,6)                                 
C        FW6(11,I)= FW6(11,I)*FCU(11,6)                                 
         FW6(12,I)= FW6(12,I)*FCU(12,6)                                 
C        FW6(13,I)= FW6(13,I)*FCU(13,6)                                 
         FW6(14,I)= FW6(14,I)*FCU(14,6)                                 
C        FW6(15,I)= FW6(15,I)*FCU(15,6)                                 
         FW6(16,I)= FW6(16,I)*FCU(16,6)                                 
         FW6(17,I)= FW6(17,I)*FCU(17,6)                                 
         FW6(18,I)= FW6(18,I)*FCU(18,6)                                 
         FW6(19,I)= FW6(19,I)*FCU(19,6)                                 
         FW6(20,I)= FW6(20,I)*FCU(20,6)                                 
         FW6(21,I)= FW6(21,I)*FCU(21,6)                                 
C        FW6(22,I)= FW6(22,I)*FCU(22,6)                                 
         FW6(23,I)= FW6(23,I)*FCU(23,6)                                 
C        FW6(24,I)= FW6(24,I)*FCU(24,6)                                 
         FW6(25,I)= FW6(25,I)*FCU(25,6)                                 
C        FW6(26,I)= FW6(26,I)*FCU(26,6)                                 
         FW6(27,I)= FW6(27,I)*FCU(27,6)                                 
C        FW6(28,I)= FW6(28,I)*FCU(28,6)                                 
      DO 630 I= 4, 5                                                    
         DO 630 J= 1,28                                                 
            RD6(J,I)= RD6(J,I)+FW6(J,4)*WORK(I+10)                      
  630 CONTINUE                                                          
C                                                                       
      CALL FRIKR7( 1, 1,FW7,FQD7, 3,FQD6, 6,FQD5,10,FQD4)               
C                                                                       
      DO J= 1,36                                                        
         RD7(J)= RD7(J)+FW7(J)*FCC(J,7)                                 
      ENDDO                                                             
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2R   *DECK INTK21                                           
C>                                                                      
C>    @brief   DDDD case                                                
C>                                                                      
C>    @details integration of the DDDD case                             
C>                                                                      
      SUBROUTINE INTK21(IKL)                                            
      use mx_limits, only: mxgsh,mxg2                                   
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
C GENERATE JTYPE=21 INTEGRALS                                           
C                                                                       
C                                                                       
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS,                      
     2                TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2),      
     3                D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2),      
     4                NGANGB                                            
C$omp threadprivate(/GEOMPQ/)
      COMMON /JMSGYH/ SQ(4)                                             
C$omp threadprivate(/JMSGYH/)
      COMMON /FQ08  / FQD(9),FQD0(5),FQD1(2,13),FQD2(3,16),FQD3(4,16),  
     2                FQD4(5,16),FQD5(6,11),FQD6(7,7),FQD7(8,3),FQD8( 9)
C$omp threadprivate(/FQ08/)
      COMMON /KI2 / ACY,ACY2,AQX,AQX2,AQXY,Y03,Y04                      
C$omp threadprivate(/KI2/)
      COMMON /KI4 / RD0(5,5),RD1(3,40),RD2(6,56),RD3(10,52),RD4(15,42), 
     *              RD5(21,24),RD6(28,12),RD7(36,4),RD8(45)             
C$omp threadprivate(/KI4/)
C                                                                       
      DIMENSION  WORK(15),FWK(36,16),FW6(28, 7),FW7(36, 3),FW8(45)      
      DIMENSION  FCU(45,8),FCC(45,8)                                    
C                                                                       
      PARAMETER (ZER=0.0D+00)                                           
      PARAMETER (F03=3.0D+00)                                           
      PARAMETER (F06=6.0D+00)                                           
      PARAMETER (F10=1.0D+01)                                           
      PARAMETER (F15=1.5D+01)                                           
C                                                                       
      IF(IKL.EQ.0) THEN                                                 
         DO I= 1, 5                                                     
            RD0(1,I)= ZER                                               
            RD0(2,I)= ZER                                               
            RD0(3,I)= ZER                                               
            RD0(4,I)= ZER                                               
            RD0(5,I)= ZER                                               
         ENDDO                                                          
         DO I= 1,40                                                     
            RD1(1,I)= ZER                                               
            RD1(2,I)= ZER                                               
            RD1(3,I)= ZER                                               
         ENDDO                                                          
         DO I= 1,56                                                     
            RD2(1,I)= ZER                                               
            RD2(2,I)= ZER                                               
            RD2(3,I)= ZER                                               
            RD2(4,I)= ZER                                               
            RD2(5,I)= ZER                                               
            RD2(6,I)= ZER                                               
         ENDDO                                                          
         DO I= 1,52                                                     
            DO J= 1,10                                                  
               RD3(J,I)= ZER                                            
            ENDDO                                                       
         ENDDO                                                          
         DO I= 1,42                                                     
            DO J= 1,15                                                  
               RD4(J,I)= ZER                                            
            ENDDO                                                       
         ENDDO                                                          
         DO I= 1,24                                                     
            DO J= 1,21                                                  
               RD5(J,I)= ZER                                            
            ENDDO                                                       
         ENDDO                                                          
         DO I= 1,12                                                     
            DO J= 1,28                                                  
               RD6(J,I)= ZER                                            
            ENDDO                                                       
         ENDDO                                                          
         DO I= 1, 4                                                     
            DO J= 1,36                                                  
               RD7(J,I)= ZER                                            
            ENDDO                                                       
         ENDDO                                                          
         DO J= 1,45                                                     
            RD8(J)= ZER                                                 
         ENDDO                                                          
C                                                                       
         RETURN                                                         
      ENDIF                                                             
C                                                                       
      XMD2= X43 *0.5D+00                                                
      XMD3= XMD2*SQ(3)                                                  
      XMD4= XMD3*XMD2                                                   
      XMD6= XMD4*XMD2                                                   
      XMDT= XMD6*XMD2                                                   
      XMD2= XMD3                                                        
C                                                                       
      XMDTY=-XMDT*ACY                                                   
      XMDTX=-XMDT*AQX                                                   
      XMDTXY=XMDT*AQXY                                                  
C                                                                       
      Y33 = Y03 *Y03                                                    
      Y34 =-Y03 *Y04                                                    
      Y44 = Y04 *Y04                                                    
      WORK( 1)= XMD4                                                    
      WORK( 2)= XMD2*Y33                                                
      WORK( 3)= XMD2*Y34                                                
      WORK( 4)= XMD2*Y44                                                
      WORK( 5)= Y33 *Y44*SQ(3)                                          
C                                                                       
      WORK( 6)=-XMD4*Y03                                                
      WORK( 7)= XMD4*Y04                                                
      WORK( 8)= XMD2*Y33*Y04                                            
      WORK( 9)= XMD2*Y34*Y04                                            
C                                                                       
      WORK(10)= XMD6                                                    
      WORK(11)= XMD4*Y33                                                
      WORK(12)= XMD4*Y34                                                
      WORK(13)= XMD4*Y44                                                
C                                                                       
      WORK(14)=-XMD6*Y03                                                
      WORK(15)= XMD6*Y04                                                
C                                                                       
      CALL FCUFCC(8,XMDT,FCU,FCC)                                       
C                                                                       
      DO 010 I= 1, 5                                                    
         DO 010 J= 1, 5                                                 
            RD0(J,I)= RD0(J,I)+FQD0(I)*WORK(J)                          
  010 CONTINUE                                                          
C                                                                       
      DO I= 1, 9                                                        
         FWK(1,I)=-FQD1(1,I)*AQX                                        
         FWK(2,I)=-FQD1(1,I)*ACY                                        
         FWK(3,I)= FQD1(2,I)                                            
      ENDDO                                                             
      DO I= 1, 4                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,1)*WORK(I+ 5)                         
         RD1(2,I)= RD1(2,I)+FWK(2,1)*WORK(I+ 5)                         
         RD1(3,I)= RD1(3,I)+FWK(3,1)*WORK(I+ 5)                         
      ENDDO                                                             
      DO I= 5, 8                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,2)*WORK(I+ 1)                         
         RD1(2,I)= RD1(2,I)+FWK(2,2)*WORK(I+ 1)                         
         RD1(3,I)= RD1(3,I)+FWK(3,2)*WORK(I+ 1)                         
      ENDDO                                                             
      DO I= 9,12                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,3)*WORK(I- 3)                         
         RD1(2,I)= RD1(2,I)+FWK(2,3)*WORK(I- 3)                         
         RD1(3,I)= RD1(3,I)+FWK(3,3)*WORK(I- 3)                         
      ENDDO                                                             
      DO I=13,16                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,4)*WORK(I- 7)                         
         RD1(2,I)= RD1(2,I)+FWK(2,4)*WORK(I- 7)                         
         RD1(3,I)= RD1(3,I)+FWK(3,4)*WORK(I- 7)                         
      ENDDO                                                             
      DO I=17,20                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,5)*WORK(I-11)                         
         RD1(2,I)= RD1(2,I)+FWK(2,5)*WORK(I-11)                         
         RD1(3,I)= RD1(3,I)+FWK(3,5)*WORK(I-11)                         
      ENDDO                                                             
      DO I=21,25                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,6)*WORK(I-20)                         
         RD1(2,I)= RD1(2,I)+FWK(2,6)*WORK(I-20)                         
         RD1(3,I)= RD1(3,I)+FWK(3,6)*WORK(I-20)                         
      ENDDO                                                             
      DO I=26,30                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,7)*WORK(I-25)                         
         RD1(2,I)= RD1(2,I)+FWK(2,7)*WORK(I-25)                         
         RD1(3,I)= RD1(3,I)+FWK(3,7)*WORK(I-25)                         
      ENDDO                                                             
      DO I=31,35                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,8)*WORK(I-30)                         
         RD1(2,I)= RD1(2,I)+FWK(2,8)*WORK(I-30)                         
         RD1(3,I)= RD1(3,I)+FWK(3,8)*WORK(I-30)                         
      ENDDO                                                             
      DO I=36,40                                                        
         RD1(1,I)= RD1(1,I)+FWK(1,9)*WORK(I-35)                         
         RD1(2,I)= RD1(2,I)+FWK(2,9)*WORK(I-35)                         
         RD1(3,I)= RD1(3,I)+FWK(3,9)*WORK(I-35)                         
      ENDDO                                                             
C                                                                       
      DO I= 1,13                                                        
         FWK(1,I)= FQD2(1,I)*AQX2-FQD1(1,I)                             
         FWK(2,I)= FQD2(1,I)*ACY2-FQD1(1,I)                             
         FWK(3,I)= FQD2(3,I)     -FQD1(1,I)                             
         FWK(4,I)= FQD2(1,I)*AQXY                                       
         FWK(5,I)=-FQD2(2,I)*AQX                                        
         FWK(6,I)=-FQD2(2,I)*ACY                                        
      ENDDO                                                             
      DO I= 1, 4                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,1)*WORK(I+ 9)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I= 5, 8                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,2)*WORK(I+ 5)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I= 9,12                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,3)*WORK(I+ 1)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=13,16                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,4)*WORK(I- 3)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=17,20                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,5)*WORK(I- 7)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=21,24                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,6)*WORK(I-15)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=25,28                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,7)*WORK(I-19)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=29,32                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,8)*WORK(I-23)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=33,36                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,9)*WORK(I-27)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=37,41                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,10)*WORK(I-36)                     
         ENDDO                                                          
      ENDDO                                                             
      DO I=42,46                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,11)*WORK(I-41)                     
         ENDDO                                                          
      ENDDO                                                             
      DO I=47,51                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,12)*WORK(I-46)                     
         ENDDO                                                          
      ENDDO                                                             
      DO I=52,56                                                        
         DO J= 1, 6                                                     
            RD2(J,I)= RD2(J,I)+FWK(J,13)*WORK(I-51)                     
         ENDDO                                                          
      ENDDO                                                             
C                                                                       
      DO I= 1,15                                                        
         FWK( 1,I)=-(FQD3(1,I)*AQX2-FQD2(1,I)*F03)*AQX                  
         FWK( 2,I)=-(FQD3(1,I)*AQX2-FQD2(1,I)    )*ACY                  
         FWK( 3,I)=  FQD3(2,I)*AQX2-FQD2(2,I)                           
         FWK( 4,I)=-(FQD3(1,I)*ACY2-FQD2(1,I)    )*AQX                  
         FWK( 5,I)=  FQD3(2,I)                    *AQXY                 
         FWK( 6,I)=-(FQD3(3,I)     -FQD2(1,I)    )*AQX                  
         FWK( 7,I)=-(FQD3(1,I)*ACY2-FQD2(1,I)*F03)*ACY                  
         FWK( 8,I)=  FQD3(2,I)*ACY2-FQD2(2,I)                           
         FWK( 9,I)=-(FQD3(3,I)     -FQD2(1,I)    )*ACY                  
         FWK(10,I)=  FQD3(4,I)     -FQD2(2,I)*F03                       
      ENDDO                                                             
      DO I= 1, 2                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,1)*WORK(I+13)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I= 3, 4                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,2)*WORK(I+11)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I= 5, 6                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,3)*WORK(I+ 9)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I= 7, 8                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,4)*WORK(I+ 7)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I= 9,10                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,5)*WORK(I+ 5)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=11,14                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,6)*WORK(I- 1)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=15,18                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,7)*WORK(I- 5)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=19,22                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,8)*WORK(I- 9)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=23,26                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,9)*WORK(I-13)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=27,30                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,10)*WORK(I-21)                     
         ENDDO                                                          
      ENDDO                                                             
      DO I=31,34                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,11)*WORK(I-25)                     
         ENDDO                                                          
      ENDDO                                                             
      DO I=35,38                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,12)*WORK(I-29)                     
         ENDDO                                                          
      ENDDO                                                             
      DO I=39,42                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,13)*WORK(I-33)                     
         ENDDO                                                          
      ENDDO                                                             
      DO I=43,47                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,14)*WORK(I-42)                     
         ENDDO                                                          
      ENDDO                                                             
      DO I=48,52                                                        
         DO J= 1,10                                                     
            RD3(J,I)= RD3(J,I)+FWK(J,15)*WORK(I-47)                     
         ENDDO                                                          
      ENDDO                                                             
C                                                                       
      AQX4= AQX2*AQX2                                                   
      ACY4= ACY2*ACY2                                                   
      X2Y2= AQX2*ACY2                                                   
      Q2C2= AQX2+ACY2                                                   
      DO I= 1, 5                                                        
         RD4( 1,I)= RD4( 1,I)+(FQD4(1,I)*AQX4-FQD3(1,I)*F06*AQX2        
     *                                       +FQD2(1,I)*F03   )*XMDT    
         RD4( 2,I)= RD4( 2,I)+(FQD4(1,I)*AQX2-FQD3(1,I)*F03   )*XMDTXY  
         RD4( 3,I)= RD4( 3,I)+(FQD4(2,I)*AQX2-FQD3(2,I)*F03   )*XMDTX   
         RD4( 4,I)= RD4( 4,I)+(FQD4(1,I)*X2Y2-FQD3(1,I)*Q2C2            
     *                                       +FQD2(1,I)       )*XMDT    
         RD4( 5,I)= RD4( 5,I)+(FQD4(2,I)*AQX2-FQD3(2,I)       )*XMDTY   
         RD4( 6,I)= RD4( 6,I)+(FQD4(3,I)*AQX2-FQD3(1,I)*AQX2            
     *                             -FQD3(3,I)+FQD2(1,I)       )*XMDT    
         RD4( 7,I)= RD4( 7,I)+(FQD4(1,I)*ACY2-FQD3(1,I)*F03   )*XMDTXY  
         RD4( 8,I)= RD4( 8,I)+(FQD4(2,I)*ACY2-FQD3(2,I)       )*XMDTX   
         RD4( 9,I)= RD4( 9,I)+(FQD4(3,I)     -FQD3(1,I)       )*XMDTXY  
         RD4(10,I)= RD4(10,I)+(FQD4(4,I)     -FQD3(2,I)*F03   )*XMDTX   
         RD4(11,I)= RD4(11,I)+(FQD4(1,I)*ACY4-FQD3(1,I)*F06*ACY2        
     *                                       +FQD2(1,I)*F03   )*XMDT    
         RD4(12,I)= RD4(12,I)+(FQD4(2,I)*ACY2-FQD3(2,I)*F03   )*XMDTY   
         RD4(13,I)= RD4(13,I)+(FQD4(3,I)*ACY2-FQD3(1,I)*ACY2            
     *                             -FQD3(3,I)+FQD2(1,I)       )*XMDT    
         RD4(14,I)= RD4(14,I)+(FQD4(4,I)     -FQD3(2,I)*F03   )*XMDTY   
         RD4(15,I)= RD4(15,I)+(FQD4(5,I)     -FQD3(3,I)*F06             
     *                                       +FQD2(1,I)*F03   )*XMDT    
      ENDDO                                                             
      DO I= 6,16                                                        
         FWK( 1,I)=  FQD4(1,I)*AQX4-FQD3(1,I)*F06*AQX2+FQD2(1,I)*F03    
         FWK( 2,I)= (FQD4(1,I)*AQX2-FQD3(1,I)*F03               )*AQXY  
         FWK( 3,I)=-(FQD4(2,I)*AQX2-FQD3(2,I)*F03               )*AQX   
         FWK( 4,I)=  FQD4(1,I)*X2Y2-FQD3(1,I)*Q2C2+FQD2(1,I)            
         FWK( 5,I)=-(FQD4(2,I)*AQX2-FQD3(2,I)                   )*ACY   
         FWK( 6,I)=  FQD4(3,I)*AQX2-FQD3(1,I)*AQX2-FQD3(3,I)+FQD2(1,I)  
         FWK( 7,I)= (FQD4(1,I)*ACY2-FQD3(1,I)*F03               )*AQXY  
         FWK( 8,I)=-(FQD4(2,I)*ACY2-FQD3(2,I)                   )*AQX   
         FWK( 9,I)= (FQD4(3,I)     -FQD3(1,I)                   )*AQXY  
         FWK(10,I)=-(FQD4(4,I)     -FQD3(2,I)*F03               )*AQX   
         FWK(11,I)=  FQD4(1,I)*ACY4-FQD3(1,I)*F06*ACY2+FQD2(1,I)*F03    
         FWK(12,I)=-(FQD4(2,I)*ACY2-FQD3(2,I)*F03               )*ACY   
         FWK(13,I)=  FQD4(3,I)*ACY2-FQD3(1,I)*ACY2-FQD3(3,I)+FQD2(1,I)  
         FWK(14,I)=-(FQD4(4,I)     -FQD3(2,I)*F03               )*ACY   
         FWK(15,I)=  FQD4(5,I)     -FQD3(3,I)*F06+FQD2(1,I)*F03         
      ENDDO                                                             
      DO I= 6, 7                                                        
         DO J= 1,15                                                     
            RD4(J,I)= RD4(J,I)+FWK(J,6)*WORK(I+ 8)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I= 8, 9                                                        
         DO J= 1,15                                                     
            RD4(J,I)= RD4(J,I)+FWK(J,7)*WORK(I+ 6)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=10,11                                                        
         DO J= 1,15                                                     
            RD4(J,I)= RD4(J,I)+FWK(J,8)*WORK(I+ 4)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=12,13                                                        
         DO J= 1,15                                                     
            RD4(J,I)= RD4(J,I)+FWK(J,9)*WORK(I+ 2)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=14,17                                                        
         DO J= 1,15                                                     
            RD4(J,I)= RD4(J,I)+FWK(J,10)*WORK(I- 4)                     
         ENDDO                                                          
      ENDDO                                                             
      DO I=18,21                                                        
         DO J= 1,15                                                     
            RD4(J,I)= RD4(J,I)+FWK(J,11)*WORK(I- 8)                     
         ENDDO                                                          
      ENDDO                                                             
      DO I=22,25                                                        
         DO J= 1,15                                                     
            RD4(J,I)= RD4(J,I)+FWK(J,12)*WORK(I-12)                     
         ENDDO                                                          
      ENDDO                                                             
      DO I=26,29                                                        
         DO J= 1,15                                                     
            RD4(J,I)= RD4(J,I)+FWK(J,13)*WORK(I-16)                     
         ENDDO                                                          
      ENDDO                                                             
      DO I=30,33                                                        
         DO J= 1,15                                                     
            RD4(J,I)= RD4(J,I)+FWK(J,14)*WORK(I-24)                     
         ENDDO                                                          
      ENDDO                                                             
      DO I=34,37                                                        
         DO J= 1,15                                                     
            RD4(J,I)= RD4(J,I)+FWK(J,15)*WORK(I-28)                     
         ENDDO                                                          
      ENDDO                                                             
      DO I=38,42                                                        
         DO J= 1,15                                                     
            RD4(J,I)= RD4(J,I)+FWK(J,16)*WORK(I-37)                     
         ENDDO                                                          
      ENDDO                                                             
C                                                                       
      DO I= 1, 4                                                        
         RD5( 1,I)= RD5( 1,I)+(FQD5(1,I)*AQX4-FQD4(1,I+5)*F10*AQX2      
     *                                       +FQD3(1,I+5)*F15  )*XMDTX  
         RD5( 2,I)= RD5( 2,I)+(FQD5(1,I)*AQX4-FQD4(1,I+5)*F06*AQX2      
     *                                       +FQD3(1,I+5)*F03  )*XMDTY  
         RD5( 3,I)= RD5( 3,I)+(FQD5(2,I)*AQX4-FQD4(2,I+5)*F06*AQX2      
     *                                       +FQD3(2,I+5)*F03  )*XMDT   
         RD5( 4,I)= RD5( 4,I)+(FQD5(1,I)*X2Y2-FQD4(1,I+5)*F03*ACY2      
     *                        -FQD4(1,I+5)*AQX2+FQD3(1,I+5)*F03)*XMDTX  
         RD5( 5,I)= RD5( 5,I)+(FQD5(2,I)*AQX2-FQD4(2,I+5)*F03  )*XMDTXY 
         RD5( 6,I)= RD5( 6,I)+(FQD5(3,I)*AQX2-FQD4(1,I+5)*AQX2          
     *                        -FQD4(3,I+5)*F03+FQD3(1,I+5)*F03 )*XMDTX  
         RD5( 7,I)= RD5( 7,I)+(FQD5(1,I)*X2Y2-FQD4(1,I+5)*F03*AQX2      
     *                        -FQD4(1,I+5)*ACY2+FQD3(1,I+5)*F03)*XMDTY  
         RD5( 8,I)= RD5( 8,I)+(FQD5(2,I)*X2Y2-FQD4(2,I+5)*Q2C2          
     *                                       +FQD3(2,I+5)      )*XMDT   
         RD5( 9,I)= RD5( 9,I)+(FQD5(3,I)*AQX2-FQD4(1,I+5)*AQX2          
     *                           -FQD4(3,I+5)+FQD3(1,I+5)      )*XMDTY  
         RD5(10,I)= RD5(10,I)+(FQD5(4,I)*AQX2-FQD4(2,I+5)*F03*AQX2      
     *                           -FQD4(4,I+5)+FQD3(2,I+5)*F03  )*XMDT   
         RD5(11,I)= RD5(11,I)+(FQD5(1,I)*ACY4-FQD4(1,I+5)*F06*ACY2      
     *                                  +FQD3(1,I+5)*F03       )*XMDTX  
         RD5(12,I)= RD5(12,I)+(FQD5(2,I)*ACY2-FQD4(2,I+5)*F03  )*XMDTXY 
         RD5(13,I)= RD5(13,I)+(FQD5(3,I)*ACY2-FQD4(1,I+5)*ACY2          
     *                           -FQD4(3,I+5)+FQD3(1,I+5)      )*XMDTX  
         RD5(14,I)= RD5(14,I)+(FQD5(4,I)-FQD4(2,I+5)*F03       )*XMDTXY 
         RD5(15,I)= RD5(15,I)+(FQD5(5,I)-FQD4(3,I+5)*F06                
     *                                  +FQD3(1,I+5)*F03       )*XMDTX  
         RD5(16,I)= RD5(16,I)+(FQD5(1,I)*ACY4-FQD4(1,I+5)*F10*ACY2      
     *                                  +FQD3(1,I+5)*F15       )*XMDTY  
         RD5(17,I)= RD5(17,I)+(FQD5(2,I)*ACY4-FQD4(2,I+5)*F06*ACY2      
     *                                  +FQD3(2,I+5)*F03       )*XMDT   
         RD5(18,I)= RD5(18,I)+(FQD5(3,I)*ACY2-FQD4(1,I+5)*ACY2          
     *                        -FQD4(3,I+5)*F03+FQD3(1,I+5)*F03 )*XMDTY  
         RD5(19,I)= RD5(19,I)+(FQD5(4,I)*ACY2-FQD4(2,I+5)*F03*ACY2      
     *                           -FQD4(4,I+5)+FQD3(2,I+5)*F03  )*XMDT   
         RD5(20,I)= RD5(20,I)+(FQD5(5,I)-FQD4(3,I+5)*F06                
     *                                  +FQD3(1,I+5)*F03       )*XMDTY  
         RD5(21,I)= RD5(21,I)+(FQD5(6,I)-FQD4(4,I+5)*F10                
     *                                  +FQD3(2,I+5)*F15       )*XMDT   
      ENDDO                                                             
      DO I= 5,11                                                        
         FWK( 1,I)=-(FQD5(1,I)*AQX4-FQD4(1,I+5)*F10*AQX2                
     *                             +FQD3(1,I+5)*F15              )*AQX  
         FWK( 2,I)=-(FQD5(1,I)*AQX4-FQD4(1,I+5)*F06*AQX2                
     *                             +FQD3(1,I+5)*F03              )*ACY  
         FWK( 3,I)=  FQD5(2,I)*AQX4-FQD4(2,I+5)*F06*AQX2                
     *                             +FQD3(2,I+5)*F03                     
         FWK( 4,I)=-(FQD5(1,I)*X2Y2-FQD4(1,I+5)*AQX2                    
     *                      -FQD4(1,I+5)*F03*ACY2+FQD3(1,I+5)*F03)*AQX  
         FWK( 5,I)= (FQD5(2,I)*AQX2-FQD4(2,I+5)*F03              )*AQXY 
         FWK( 6,I)=-(FQD5(3,I)*AQX2-FQD4(1,I+5)*AQX2-FQD4(3,I+5)*F03    
     *                             +FQD3(1,I+5)*F03              )*AQX  
         FWK( 7,I)=-(FQD5(1,I)*X2Y2-FQD4(1,I+5)*F03*AQX2                
     *                    -FQD4(1,I+5)*ACY2+FQD3(1,I+5)*F03      )*ACY  
         FWK( 8,I)=  FQD5(2,I)*X2Y2-FQD4(2,I+5)*Q2C2+FQD3(2,I+5)        
         FWK( 9,I)=-(FQD5(3,I)*AQX2-FQD4(1,I+5)*AQX2-FQD4(3,I+5)        
     *                             +FQD3(1,I+5)                  )*ACY  
         FWK(10,I)=  FQD5(4,I)*AQX2-FQD4(2,I+5)*F03*AQX2-FQD4(4,I+5)    
     *                             +FQD3(2,I+5)*F03                     
         FWK(11,I)=-(FQD5(1,I)*ACY4-FQD4(1,I+5)*F06*ACY2                
     *                             +FQD3(1,I+5)*F03              )*AQX  
         FWK(12,I)= (FQD5(2,I)*ACY2-FQD4(2,I+5)*F03              )*AQXY 
         FWK(13,I)=-(FQD5(3,I)*ACY2-FQD4(1,I+5)*ACY2-FQD4(3,I+5)        
     *                             +FQD3(1,I+5)                  )*AQX  
         FWK(14,I)= (FQD5(4,I)     -FQD4(2,I+5)*F03              )*AQXY 
         FWK(15,I)=-(FQD5(5,I)     -FQD4(3,I+5)*F06                     
     *                             +FQD3(1,I+5)*F03              )*AQX  
         FWK(16,I)=-(FQD5(1,I)*ACY4-FQD4(1,I+5)*F10*ACY2                
     *                             +FQD3(1,I+5)*F15              )*ACY  
         FWK(17,I)=  FQD5(2,I)*ACY4-FQD4(2,I+5)*F06*ACY2                
     *                             +FQD3(2,I+5)*F03                     
         FWK(18,I)=-(FQD5(3,I)*ACY2-FQD4(1,I+5)*ACY2-FQD4(3,I+5)*F03    
     *                             +FQD3(1,I+5)*F03              )*ACY  
         FWK(19,I)=  FQD5(4,I)*ACY2-FQD4(2,I+5)*F03*ACY2-FQD4(4,I+5)    
     *                             +FQD3(2,I+5)*F03                     
         FWK(20,I)=-(FQD5(5,I)-FQD4(3,I+5)*F06+FQD3(1,I+5)*F03   )*ACY  
         FWK(21,I)=  FQD5(6,I)-FQD4(4,I+5)*F10+FQD3(2,I+5)*F15          
      ENDDO                                                             
      DO I= 5, 6                                                        
         DO J= 1,21                                                     
            RD5(J,I)= RD5(J,I)+FWK(J,5)*WORK(I+ 9)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I= 7, 8                                                        
         DO J= 1,21                                                     
            RD5(J,I)= RD5(J,I)+FWK(J,6)*WORK(I+ 7)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I= 9,10                                                        
         DO J= 1,21                                                     
            RD5(J,I)= RD5(J,I)+FWK(J,7)*WORK(I+ 5)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=11,12                                                        
         DO J= 1,21                                                     
            RD5(J,I)= RD5(J,I)+FWK(J,8)*WORK(I+ 3)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=13,16                                                        
         DO J= 1,21                                                     
            RD5(J,I)= RD5(J,I)+FWK(J,9)*WORK(I- 3)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I=17,20                                                        
         DO J= 1,21                                                     
            RD5(J,I)= RD5(J,I)+FWK(J,10)*WORK(I- 7)                     
         ENDDO                                                          
      ENDDO                                                             
      DO I=21,24                                                        
         DO J= 1,21                                                     
            RD5(J,I)= RD5(J,I)+FWK(J,11)*WORK(I-15)                     
         ENDDO                                                          
      ENDDO                                                             
C                                                                       
      CALL FRIKR6( 1, 7,FW6,FQD6, 4,FQD5, 9,FQD4, 9,FQD3)               
C                                                                       
      DO 620 I= 1, 4                                                    
         DO 620 J= 1,28                                                 
            RD6(J,I)= RD6(J,I)+FW6(J,I)*FCC(J,6)                        
  620 CONTINUE                                                          
      DO I= 5, 7                                                        
C        FW6( 1,I)= FW6( 1,I)*FCU( 1,6)                                 
         FW6( 2,I)= FW6( 2,I)*FCU( 2,6)                                 
         FW6( 3,I)= FW6( 3,I)*FCU( 3,6)                                 
C        FW6( 4,I)= FW6( 4,I)*FCU( 4,6)                                 
         FW6( 5,I)= FW6( 5,I)*FCU( 5,6)                                 
C        FW6( 6,I)= FW6( 6,I)*FCU( 6,6)                                 
         FW6( 7,I)= FW6( 7,I)*FCU( 7,6)                                 
         FW6( 8,I)= FW6( 8,I)*FCU( 8,6)                                 
         FW6( 9,I)= FW6( 9,I)*FCU( 9,6)                                 
         FW6(10,I)= FW6(10,I)*FCU(10,6)                                 
C        FW6(11,I)= FW6(11,I)*FCU(11,6)                                 
         FW6(12,I)= FW6(12,I)*FCU(12,6)                                 
C        FW6(13,I)= FW6(13,I)*FCU(13,6)                                 
         FW6(14,I)= FW6(14,I)*FCU(14,6)                                 
C        FW6(15,I)= FW6(15,I)*FCU(15,6)                                 
         FW6(16,I)= FW6(16,I)*FCU(16,6)                                 
         FW6(17,I)= FW6(17,I)*FCU(17,6)                                 
         FW6(18,I)= FW6(18,I)*FCU(18,6)                                 
         FW6(19,I)= FW6(19,I)*FCU(19,6)                                 
         FW6(20,I)= FW6(20,I)*FCU(20,6)                                 
         FW6(21,I)= FW6(21,I)*FCU(21,6)                                 
C        FW6(22,I)= FW6(22,I)*FCU(22,6)                                 
         FW6(23,I)= FW6(23,I)*FCU(23,6)                                 
C        FW6(24,I)= FW6(24,I)*FCU(24,6)                                 
         FW6(25,I)= FW6(25,I)*FCU(25,6)                                 
C        FW6(26,I)= FW6(26,I)*FCU(26,6)                                 
         FW6(27,I)= FW6(27,I)*FCU(27,6)                                 
C        FW6(28,I)= FW6(28,I)*FCU(28,6)                                 
      ENDDO                                                             
      DO I= 5, 6                                                        
         DO J= 1,28                                                     
            RD6(J,I)= RD6(J,I)+FW6(J,5)*WORK(I+ 9)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I= 7, 8                                                        
         DO J= 1,28                                                     
            RD6(J,I)= RD6(J,I)+FW6(J,6)*WORK(I+ 7)                      
         ENDDO                                                          
      ENDDO                                                             
      DO I= 9,12                                                        
         DO J= 1,28                                                     
            RD6(J,I)= RD6(J,I)+FW6(J,7)*WORK(I+ 1)                      
         ENDDO                                                          
      ENDDO                                                             
C                                                                       
      CALL FRIKR7( 1, 3,FW7,FQD7, 4,FQD6, 8,FQD5,13,FQD4)               
C                                                                       
      DO 720 I= 1, 2                                                    
         DO 720 J= 1,36                                                 
            RD7(J,I)= RD7(J,I)+FW7(J,I)*FCC(J,7)                        
  720 CONTINUE                                                          
      I= 3                                                              
         FW7( 1,I)= FW7( 1,I)*FCU( 1,7)                                 
         FW7( 2,I)= FW7( 2,I)*FCU( 2,7)                                 
C        FW7( 3,I)= FW7( 3,I)*FCU( 3,7)                                 
         FW7( 4,I)= FW7( 4,I)*FCU( 4,7)                                 
         FW7( 5,I)= FW7( 5,I)*FCU( 5,7)                                 
         FW7( 6,I)= FW7( 6,I)*FCU( 6,7)                                 
         FW7( 7,I)= FW7( 7,I)*FCU( 7,7)                                 
C        FW7( 8,I)= FW7( 8,I)*FCU( 8,7)                                 
         FW7( 9,I)= FW7( 9,I)*FCU( 9,7)                                 
C        FW7(10,I)= FW7(10,I)*FCU(10,7)                                 
         FW7(11,I)= FW7(11,I)*FCU(11,7)                                 
         FW7(12,I)= FW7(12,I)*FCU(12,7)                                 
         FW7(13,I)= FW7(13,I)*FCU(13,7)                                 
         FW7(14,I)= FW7(14,I)*FCU(14,7)                                 
         FW7(15,I)= FW7(15,I)*FCU(15,7)                                 
         FW7(16,I)= FW7(16,I)*FCU(16,7)                                 
C        FW7(17,I)= FW7(17,I)*FCU(17,7)                                 
         FW7(18,I)= FW7(18,I)*FCU(18,7)                                 
C        FW7(19,I)= FW7(19,I)*FCU(19,7)                                 
         FW7(20,I)= FW7(20,I)*FCU(20,7)                                 
C        FW7(21,I)= FW7(21,I)*FCU(21,7)                                 
         FW7(22,I)= FW7(22,I)*FCU(22,7)                                 
         FW7(23,I)= FW7(23,I)*FCU(23,7)                                 
         FW7(24,I)= FW7(24,I)*FCU(24,7)                                 
         FW7(25,I)= FW7(25,I)*FCU(25,7)                                 
         FW7(26,I)= FW7(26,I)*FCU(26,7)                                 
         FW7(27,I)= FW7(27,I)*FCU(27,7)                                 
         FW7(28,I)= FW7(28,I)*FCU(28,7)                                 
         FW7(29,I)= FW7(29,I)*FCU(29,7)                                 
C        FW7(30,I)= FW7(30,I)*FCU(30,7)                                 
         FW7(31,I)= FW7(31,I)*FCU(31,7)                                 
C        FW7(32,I)= FW7(32,I)*FCU(32,7)                                 
         FW7(33,I)= FW7(33,I)*FCU(33,7)                                 
C        FW7(34,I)= FW7(34,I)*FCU(34,7)                                 
         FW7(35,I)= FW7(35,I)*FCU(35,7)                                 
C        FW7(36,I)= FW7(36,I)*FCU(36,7)                                 
      DO 730 I= 3, 4                                                    
         DO 730 J= 1,36                                                 
            RD7(J,I)= RD7(J,I)+FW7(J,3)*WORK(I+11)                      
  730 CONTINUE                                                          
C                                                                       
      CALL FRIKR8( 1, 1,FW8,FQD8, 2,FQD7, 6,FQD6,10,FQD5,15,FQD4)       
C                                                                       
      DO J= 1,45                                                        
         RD8(J)= RD8(J)+FW8(J)*FCC(J,8)                                 
      ENDDO                                                             
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2R   *DECK FCUFCC                                           
C>                                                                      
C>    @brief   auxiliary routine internal to rot.axis integrations      
C>                                                                      
C>    @details auxiliary routine internal to rot.axis integrations      
C>                                                                      
      SUBROUTINE FCUFCC(N,XMDT,FCU,FCC)                                 
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      DIMENSION  FCU(45,8),FCC(45,8)                                    
C                                                                       
      COMMON /KI2 / ACY,ACY2,AQX,AQX2,AQXY,Y03,Y04                      
C$omp threadprivate(/KI2/)
C                                                                       
      PARAMETER (ONE=1.0D+00)                                           
C                                                                       
      K= 1                                                              
      DO I= 1,N                                                         
         K= K+I+1                                                       
         DO J= 1, K                                                     
            FCU(J,I)= ONE                                               
            FCC(J,I)= XMDT                                              
         ENDDO                                                          
      ENDDO                                                             
C                                                                       
      XMDTY=-XMDT*ACY                                                   
      XMDTX=-XMDT*AQX                                                   
      XMDTXY=XMDT*AQXY                                                  
C                                                                       
         FCU( 1,1)=-AQX                                                 
         FCU( 2,1)=-ACY                                                 
C                                                                       
         FCC( 1,1)= XMDTX                                               
         FCC( 2,1)= XMDTY                                               
C     IF(N.LE.1) GO TO 200                                              
         FCU( 4,2)= AQXY                                                
         FCU( 5,2)=-AQX                                                 
         FCU( 6,2)=-ACY                                                 
C                                                                       
         FCC( 4,2)= XMDTXY                                              
         FCC( 5,2)= XMDTX                                               
         FCC( 6,2)= XMDTY                                               
C     IF(N.LE.2) GO TO 200                                              
         FCU( 1,3)=-AQX                                                 
         FCU( 2,3)=-ACY                                                 
         FCU( 4,3)=-AQX                                                 
         FCU( 5,3)= AQXY                                                
         FCU( 6,3)=-AQX                                                 
         FCU( 7,3)=-ACY                                                 
         FCU( 9,3)=-ACY                                                 
C                                                                       
         FCC( 1,3)= XMDTX                                               
         FCC( 2,3)= XMDTY                                               
         FCC( 4,3)= XMDTX                                               
         FCC( 5,3)= XMDTXY                                              
         FCC( 6,3)= XMDTX                                               
         FCC( 7,3)= XMDTY                                               
         FCC( 9,3)= XMDTY                                               
C     IF(N.LE.3) GO TO 200                                              
         FCU( 2,4)= AQXY                                                
         FCU( 3,4)=-AQX                                                 
         FCU( 5,4)=-ACY                                                 
         FCU( 7,4)= AQXY                                                
         FCU( 8,4)=-AQX                                                 
         FCU( 9,4)= AQXY                                                
         FCU(10,4)=-AQX                                                 
         FCU(12,4)=-ACY                                                 
         FCU(14,4)=-ACY                                                 
C                                                                       
         FCC( 2,4)= XMDTXY                                              
         FCC( 3,4)= XMDTX                                               
         FCC( 5,4)= XMDTY                                               
         FCC( 7,4)= XMDTXY                                              
         FCC( 8,4)= XMDTX                                               
         FCC( 9,4)= XMDTXY                                              
         FCC(10,4)= XMDTX                                               
         FCC(12,4)= XMDTY                                               
         FCC(14,4)= XMDTY                                               
C     IF(N.LE.4) GO TO 200                                              
         DO J= 1,10                                                     
            FCU(J,5)= FCU(J,3)                                          
            FCC(J,5)= FCC(J,3)                                          
         ENDDO                                                          
         FCU(11,5)=-AQX                                                 
         FCU(12,5)= AQXY                                                
         FCU(13,5)=-AQX                                                 
         FCU(14,5)= AQXY                                                
         FCU(15,5)=-AQX                                                 
         FCU(16,5)=-ACY                                                 
         FCU(18,5)=-ACY                                                 
         FCU(20,5)=-ACY                                                 
C                                                                       
         FCC(11,5)= XMDTX                                               
         FCC(12,5)= XMDTXY                                              
         FCC(13,5)= XMDTX                                               
         FCC(14,5)= XMDTXY                                              
         FCC(15,5)= XMDTX                                               
         FCC(16,5)= XMDTY                                               
         FCC(18,5)= XMDTY                                               
         FCC(20,5)= XMDTY                                               
C     IF(N.LE.5) GO TO 200                                              
         DO J= 1,15                                                     
            FCU(J,6)= FCU(J,4)                                          
            FCC(J,6)= FCC(J,4)                                          
         ENDDO                                                          
         FCU(16,6)= AQXY                                                
         FCU(17,6)=-AQX                                                 
         FCU(18,6)= AQXY                                                
         FCU(19,6)=-AQX                                                 
         FCU(20,6)= AQXY                                                
         FCU(21,6)=-AQX                                                 
         FCU(23,6)=-ACY                                                 
         FCU(25,6)=-ACY                                                 
         FCU(27,6)=-ACY                                                 
C                                                                       
         FCC(16,6)= XMDTXY                                              
         FCC(17,6)= XMDTX                                               
         FCC(18,6)= XMDTXY                                              
         FCC(19,6)= XMDTX                                               
         FCC(20,6)= XMDTXY                                              
         FCC(21,6)= XMDTX                                               
         FCC(23,6)= XMDTY                                               
         FCC(25,6)= XMDTY                                               
         FCC(27,6)= XMDTY                                               
      IF(N.LE.6) GO TO 200                                              
         DO J= 1,21                                                     
            FCU(J,7)= FCU(J,5)                                          
            FCC(J,7)= FCC(J,5)                                          
         ENDDO                                                          
         FCU(22,7)=-AQX                                                 
         FCU(23,7)= AQXY                                                
         FCU(24,7)=-AQX                                                 
         FCU(25,7)= AQXY                                                
         FCU(26,7)=-AQX                                                 
         FCU(27,7)= AQXY                                                
         FCU(28,7)=-AQX                                                 
         FCU(29,7)=-ACY                                                 
         FCU(31,7)=-ACY                                                 
         FCU(33,7)=-ACY                                                 
         FCU(35,7)=-ACY                                                 
C                                                                       
         FCC(22,7)= XMDTX                                               
         FCC(23,7)= XMDTXY                                              
         FCC(24,7)= XMDTX                                               
         FCC(25,7)= XMDTXY                                              
         FCC(26,7)= XMDTX                                               
         FCC(27,7)= XMDTXY                                              
         FCC(28,7)= XMDTX                                               
         FCC(29,7)= XMDTY                                               
         FCC(31,7)= XMDTY                                               
         FCC(33,7)= XMDTY                                               
         FCC(35,7)= XMDTY                                               
      IF(N.LE.7) GO TO 200                                              
         DO J= 1,28                                                     
            FCU(J,8)= FCU(J,6)                                          
            FCC(J,8)= FCC(J,6)                                          
         ENDDO                                                          
         FCU(29,8)= AQXY                                                
         FCU(30,8)=-AQX                                                 
         FCU(31,8)= AQXY                                                
         FCU(32,8)=-AQX                                                 
         FCU(33,8)= AQXY                                                
         FCU(34,8)=-AQX                                                 
         FCU(35,8)= AQXY                                                
         FCU(36,8)=-AQX                                                 
         FCU(38,8)=-ACY                                                 
         FCU(40,8)=-ACY                                                 
         FCU(42,8)=-ACY                                                 
         FCU(44,8)=-ACY                                                 
C                                                                       
         FCC(29,8)= XMDTXY                                              
         FCC(30,8)= XMDTX                                               
         FCC(31,8)= XMDTXY                                              
         FCC(32,8)= XMDTX                                               
         FCC(33,8)= XMDTXY                                              
         FCC(34,8)= XMDTX                                               
         FCC(35,8)= XMDTXY                                              
         FCC(36,8)= XMDTX                                               
         FCC(38,8)= XMDTY                                               
         FCC(40,8)= XMDTY                                               
         FCC(42,8)= XMDTY                                               
         FCC(44,8)= XMDTY                                               
C                                                                       
  200 CONTINUE                                                          
      RETURN                                                            
      END                                                               
C*MODULE INT2B   *DECK FRIKR6                                           
C>                                                                      
C>    @brief   auxiliary routine of order 6 for rot.axis integrations   
C>                                                                      
C>    @details auxiliary routine of order 6 for rot.axis integrations   
C>                                                                      
      SUBROUTINE FRIKR6(I1,I2,WRK,QD6,J0,QD5,K0,QD4,L0,QD3)             
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      DIMENSION WRK(28,*)                                               
      DIMENSION QD6( 7,*),QD5( 6,*),QD4( 5,*),QD3( 4,*)                 
C                                                                       
      COMMON /KI2 / ACY,ACY2,AQX,AQX2,AQXY,Y03,Y04                      
C$omp threadprivate(/KI2/)
C                                                                       
      PARAMETER (F03=3.0D+00 )                                          
      PARAMETER (F06=6.0D+00)                                           
      PARAMETER (F10=1.0D+01 )                                          
      PARAMETER (F15=1.5D+01)                                           
C                                                                       
      DO I=I1,I2                                                        
         J=I+J0                                                         
         K=I+K0                                                         
         L=I+L0                                                         
         F31L03= QD3(1,L)*F03                                           
         F31L15= QD3(1,L)*F15                                           
C                                                                       
         F41K03= QD4(1,K)*F03                                           
         F41K06= F41K03+F41K03                                          
         F41K15= QD4(1,K)*F15                                           
         F41K45= F41K15  *F03                                           
         F42K03= QD4(2,K)*F03                                           
         F42K15= QD4(2,K)*F15                                           
         F43K03= QD4(3,K)*F03                                           
         F43K06= F43K03+F43K03                                          
         F43K45= F43K03  *F15                                           
C                                                                       
         F51J03= QD5(1,J)*F03                                           
         F51J06= F51J03+F51J03                                          
         F51J10= QD5(1,J)*F10                                           
         F51J15= QD5(1,J)*F15                                           
         F52J03= QD5(2,J)*F03                                           
         F52J06= F52J03+F52J03                                          
         F52J10= QD5(2,J)*F10                                           
         F53J03= QD5(3,J)*F03                                           
         F53J06= F53J03+F53J03                                          
         F54J03= QD5(4,J)*F03                                           
         F54J10= QD5(4,J)*F10                                           
         F55J15= QD5(5,J)*F15                                           
C                                                                       
         A11   =     QD4(1,K)*AQX2-QD3(1,L)                             
         R11   =     F41K03  *ACY2-F31L03                               
C                                                                       
         B11   =     QD5(1,J)*AQX2-QD4(1,K)                             
         B13   =     QD5(1,J)*AQX2-F41K03                               
         B22   =     F52J03  *AQX2-F42K03                               
         B23   =     QD5(2,J)*AQX2-F42K03                               
         B33   =     QD5(3,J)*AQX2-QD4(3,K)                             
         S13   =     QD5(1,J)*ACY2-F41K03                               
         S22   =     F52J03  *ACY2-F42K03                               
         S23   =     F52J03  *ACY2-F42K03  *F03                         
         S33   =     F53J06  *ACY2-F43K06                               
C                                                                       
         D11   =    (QD5(1,J)*AQX2-F41K06  )*AQX2+F31L03                
         U11   =    (QD5(1,J)*ACY2-F41K06  )*ACY2+F31L03                
C                                                                       
         WRK( 1,I)=((QD6(1,I)*AQX2-F51J15  )*AQX2+F41K45)*AQX2-F31L15   
         WRK( 2,I)= (QD6(1,I)*AQX2-F51J10  )*AQX2+F41K15                
         WRK( 3,I)= (QD6(2,I)*AQX2-F52J10  )*AQX2+F42K15                
         WRK( 4,I)=((QD6(1,I)*AQX2-F51J06  )*AQX2+F41K03)*ACY2-D11      
         WRK( 5,I)= (QD6(2,I)*AQX2-F52J06  )*AQX2+F42K03                
         WRK( 6,I)= (QD6(3,I)*AQX2-F53J06  )*AQX2+F43K03      -D11      
         WRK( 7,I)= (QD6(1,I)*AQX2-F51J03  )*ACY2-B13*F03               
         WRK( 8,I)= (QD6(2,I)*AQX2-F52J03  )*ACY2-B23                   
         WRK( 9,I)=  QD6(3,I)*AQX2-F53J03        -B13                   
         WRK(10,I)=  QD6(4,I)*AQX2-F54J03        -B23*F03               
         WRK(11,I)=((QD6(1,I)*ACY2-F51J06  )*ACY2+F41K03)*AQX2-U11      
         WRK(12,I)= (QD6(2,I)*AQX2-QD5(2,J))*ACY2-B22                   
         WRK(13,I)= (QD6(3,I)*AQX2-QD5(3,J)      -B11   )*ACY2-B33+A11  
         WRK(14,I)=  QD6(4,I)*AQX2-QD5(4,J)      -B22                   
         WRK(15,I)=  QD6(5,I)*AQX2-QD5(5,J)      -B33*F06 +A11*F03      
         WRK(16,I)= (QD6(1,I)*ACY2-F51J10  )*ACY2+F41K15                
         WRK(17,I)= (QD6(2,I)*ACY2-F52J06  )*ACY2+F42K03                
         WRK(18,I)=  QD6(3,I)*ACY2-F53J03        -S13                   
         WRK(19,I)=  QD6(4,I)*ACY2-QD5(4,J)      -S22                   
         WRK(20,I)=  QD6(5,I)     -F53J06        +F41K03                
         WRK(21,I)=  QD6(6,I)     -F54J10        +F42K15                
         WRK(22,I)=((QD6(1,I)*ACY2-F51J15  )*ACY2+F41K45)*ACY2-F31L15   
         WRK(23,I)= (QD6(2,I)*ACY2-F52J10  )*ACY2+F42K15                
         WRK(24,I)= (QD6(3,I)*ACY2-F53J06  )*ACY2+F43K03      -U11      
         WRK(25,I)=  QD6(4,I)*ACY2-F54J03        -S23                   
         WRK(26,I)=  QD6(5,I)*ACY2-QD5(5,J)      -S33         +R11      
         WRK(27,I)=  QD6(6,I)     -F54J10        +F42K15                
         WRK(28,I)=  QD6(7,I)     -F55J15        +F43K45      -F31L15   
      ENDDO                                                             
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2R   *DECK FRIKR7                                           
C>                                                                      
C>    @brief   auxiliary routine of order 7 for rot.axis integrations   
C>                                                                      
C>    @details auxiliary routine of order 7 for rot.axis integrations   
C>                                                                      
      SUBROUTINE FRIKR7(I1,I2,WRK,QD7,J0,QD6,K0,QD5,L0,QD4)             
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      DIMENSION WRK(36,*)                                               
      DIMENSION QD7( 8,*),QD6( 7,*),QD5( 6,*),QD4( 5,*)                 
C                                                                       
      COMMON /KI2 / ACY,ACY2,AQX,AQX2,AQXY,Y03,Y04                      
C$omp threadprivate(/KI2/)
C                                                                       
      PARAMETER (F03=3.0D+00)                                           
      PARAMETER (F05=5.0D+00)                                           
      PARAMETER (F06=6.0D+00)                                           
      PARAMETER (F07=7.0D+00)                                           
      PARAMETER (F10=1.0D+01)                                           
      PARAMETER (F15=1.5D+01)                                           
      PARAMETER (F21=2.1D+01)                                           
C                                                                       
      DO I=I1,I2                                                        
         J=I+J0                                                         
         K=I+K0                                                         
         L=I+L0                                                         
         F41L03= QD4(1,L)*F03                                           
         F41L15= QD4(1,L)*F15                                           
         F41L1H= F41L15  *F07                                           
         F42L03= QD4(2,L)*F03                                           
         F42L15= QD4(2,L)*F15                                           
         F42L1H= F42L15  *F07                                           
C                                                                       
         F51K03= QD5(1,K)*F03                                           
         F51K06= F51K03+F51K03                                          
         F51K10= QD5(1,K)*F10                                           
         F51K15= QD5(1,K)*F15                                           
         F51K45= F51K15  *F03                                           
         F51K1H= F51K15  *F07                                           
         F52K03= QD5(2,K)*F03                                           
         F52K06= F52K03+F52K03                                          
         F52K09= F52K03+F52K06                                          
         F52K15= QD5(2,K)*F15                                           
         F52K45= F52K15  *F03                                           
         F53K03= QD5(3,K)*F03                                           
         F53K15= QD5(3,K)*F15                                           
         F53K45= F53K15  *F03                                           
         F54K03= QD5(4,K)*F03                                           
         F54K1H= F54K03  *F05*F07                                       
C                                                                       
         F61J06= QD6(1,J)*F06                                           
         F61J10= QD6(1,J)*F10                                           
         F61J15= QD6(1,J)*F15                                           
         F61J21= QD6(1,J)*F21                                           
         F62J03= QD6(2,J)*F03                                           
         F62J06= F62J03+F62J03                                          
         F62J10= QD6(2,J)*F10                                           
         F62J15= QD6(2,J)*F15                                           
         F63J03= QD6(3,J)*F03                                           
         F63J06= F63J03+F63J03                                          
         F63J10= QD6(3,J)*F10                                           
         F64J03= QD6(4,J)*F03                                           
         F64J06= F64J03+F64J03                                          
         F64J10= QD6(4,J)*F10                                           
         F65J03= QD6(5,J)*F03                                           
         F65J15= QD6(5,J)*F15                                           
         F66J21= QD6(6,J)*F21                                           
C                                                                       
         A11   =     QD5(1,K)*AQX2-QD4(1,L)                             
         A13   =     QD5(1,K)*AQX2-F41L03                               
         A22   =     QD5(2,K)*AQX2-QD4(2,L)                             
         R11   =     F51K03  *ACY2-F41L03                               
         R13   =     QD5(1,K)*ACY2-F41L03                               
         R22   =     F52K03  *ACY2-F42L03                               
C                                                                       
         B23   =     F62J03  *AQX2-F52K09                               
         B33   =     QD6(3,J)*AQX2-QD5(3,K)                             
         B3T   =     QD6(3,J)*AQX2-F53K03                               
         B44   =     QD6(4,J)*AQX2-QD5(4,K)                             
         S11   =     QD6(1,J)*ACY2-QD5(1,K)                             
         S13   =     QD6(1,J)*ACY2-F51K03                               
         S22   =     F62J03  *ACY2-F52K03                               
         S23   =     F62J03  *ACY2-F52K09                               
         S33   =     F63J03  *ACY2-F53K03                               
         S3T   =     QD6(3,J)*ACY2-F53K03                               
         S44   =     QD6(4,J)*ACY2-QD5(4,K)                             
C                                                                       
         D1S   =    (QD6(1,J)*AQX2-F51K06  )*AQX2+F41L03                
         D1T   =    (QD6(1,J)*AQX2-F51K10  )*AQX2+F41L15                
         D2S   =    (QD6(2,J)*AQX2-F52K06  )*AQX2+F42L03                
         U1S   =    (QD6(1,J)*ACY2-F51K06  )*ACY2+F41L03                
         U1T   =    (QD6(1,J)*ACY2-F51K10  )*ACY2+F41L15                
         U2S   =    (QD6(2,J)*ACY2-F52K06  )*ACY2+F42L03                
C                                                                       
         WRK( 1,I)=((QD7(1,I)*AQX2-F61J21  )*AQX2+F51K1H)*AQX2-F41L1H   
         WRK( 2,I)=((QD7(1,I)*AQX2-F61J15  )*AQX2+F51K45)*AQX2-F41L15   
         WRK( 3,I)=((QD7(2,I)*AQX2-F62J15  )*AQX2+F52K45)*AQX2-F42L15   
         WRK( 4,I)=((QD7(1,I)*AQX2-F61J10  )*AQX2+F51K15)*ACY2-D1T      
         WRK( 5,I)= (QD7(2,I)*AQX2-F62J10  )*AQX2+F52K15                
         WRK( 6,I)= (QD7(3,I)*AQX2-F63J10  )*AQX2+F53K15      -D1T      
         WRK( 7,I)=((QD7(1,I)*AQX2-F61J06  )*AQX2+F51K03)*ACY2-D1S*F03  
         WRK( 8,I)=((QD7(2,I)*AQX2-F62J06  )*AQX2+F52K03)*ACY2-D2S      
         WRK( 9,I)= (QD7(3,I)*AQX2-F63J06  )*AQX2+F53K03      -D1S      
         WRK(10,I)= (QD7(4,I)*AQX2-F64J06  )*AQX2+F54K03      -D2S*F03  
         WRK(11,I)=((QD7(1,I)*ACY2-F61J06  )*ACY2+F51K03)*AQX2-U1S*F03  
         WRK(12,I)= (QD7(2,I)*ACY2-F62J03  )*AQX2-S23                   
         WRK(13,I)= (QD7(3,I)*ACY2-QD6(3,J)      -S11   )*AQX2-S33+R11  
         WRK(14,I)=  QD7(4,I)*AQX2-F64J03        -B23                   
         WRK(15,I)=  QD7(5,I)*AQX2-F65J03        -B3T*F06     +A13*F03  
         WRK(16,I)=((QD7(1,I)*ACY2-F61J10  )*ACY2+F51K15)*AQX2-U1T      
         WRK(17,I)=((QD7(2,I)*ACY2-F62J06  )*ACY2+F52K03)*AQX2-U2S      
         WRK(18,I)= (QD7(3,I)*ACY2-F63J03        -S13   )*AQX2-S3T+R13  
         WRK(19,I)= (QD7(4,I)*ACY2-QD6(4,J)      -S22   )*AQX2-S44+R22  
         WRK(20,I)=  QD7(5,I)*AQX2-QD6(5,J)      -B33*F06     +A11*F03  
         WRK(21,I)=  QD7(6,I)*AQX2-QD6(6,J)      -B44*F10     +A22*F15  
         WRK(22,I)=((QD7(1,I)*ACY2-F61J15  )*ACY2+F51K45)*ACY2-F41L15   
         WRK(23,I)= (QD7(2,I)*ACY2-F62J10  )*ACY2+F52K15                
         WRK(24,I)= (QD7(3,I)*ACY2-F63J06  )*ACY2+F53K03      -U1S      
         WRK(25,I)=  QD7(4,I)*ACY2-F64J03        -S23                   
         WRK(26,I)=  QD7(5,I)*ACY2-QD6(5,J)      -S33-S33     +R11      
         WRK(27,I)=  QD7(6,I)     -F64J10        +F52K15                
         WRK(28,I)=  QD7(7,I)     -F65J15        +F53K45      -F41L15   
         WRK(29,I)=((QD7(1,I)*ACY2-F61J21  )*ACY2+F51K1H)*ACY2-F41L1H   
         WRK(30,I)=((QD7(2,I)*ACY2-F62J15  )*ACY2+F52K45)*ACY2-F42L15   
         WRK(31,I)= (QD7(3,I)*ACY2-F63J10  )*ACY2+F53K15      -U1T      
         WRK(32,I)= (QD7(4,I)*ACY2-F64J06  )*ACY2+F54K03      -U2S*F03  
         WRK(33,I)=  QD7(5,I)*ACY2-F65J03        -S3T*F06     +R13*F03  
         WRK(34,I)=  QD7(6,I)*ACY2-QD6(6,J)      -S44*F10     +R22*F05  
         WRK(35,I)=  QD7(7,I)     -F65J15        +F53K45      -F41L15   
         WRK(36,I)=  QD7(8,I)     -F66J21        +F54K1H      -F42L1H   
      ENDDO                                                             
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2R   *DECK FRIKR8                                           
C>                                                                      
C>    @brief   auxiliary routine of order 8 for rot.axis integrations   
C>                                                                      
C>    @details auxiliary routine of order 8 for rot.axis integrations   
C>                                                                      
      SUBROUTINE FRIKR8(I1,I2,WRK,QD8,J0,QD7,K0,QD6,L0,QD5,M0,QD4)      
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      DIMENSION WRK(45,*)                                               
      DIMENSION QD8( 9,*),QD7( 8,*),QD6( 7,*),QD5( 6,*),QD4( 5,*)       
C                                                                       
      COMMON /KI2 / ACY,ACY2,AQX,AQX2,AQXY,Y03,Y04                      
C$omp threadprivate(/KI2/)
C                                                                       
      PARAMETER (F03=3.0D+00)                                           
      PARAMETER (F05=5.0D+00 )                                          
      PARAMETER (F06=6.0D+00)                                           
      PARAMETER (F07=7.0D+00)                                           
      PARAMETER (F10=1.0D+01 )                                          
      PARAMETER (F15=1.5D+01 )                                          
      PARAMETER (F21=2.1D+01)                                           
      PARAMETER (F28=2.8D+01)                                           
      PARAMETER (F210=2.1D+02)                                          
      PARAMETER (F420=4.2D+02)                                          
C                                                                       
      CY= ACY2                                                          
      QX= AQX2                                                          
      DO I=I1,I2                                                        
         J=I+J0                                                         
         K=I+K0                                                         
         L=I+L0                                                         
         M=I+M0                                                         
         F41M03= QD4(1,M)*F03                                           
         F41M15= QD4(1,M)*F15                                           
         F41M1H= F41M15  *F07                                           
C                                                                       
         F51L03= QD5(1,L)*F03                                           
         F51L06= F51L03+F51L03                                          
         F51L09= F51L03+F51L06                                          
         F51L15= QD5(1,L)*F15                                           
         F51L45= F51L15  *F03                                           
         F51L1H= F51L15  *F07                                           
         F51L4H= QD5(1,L)*F420                                          
         F52L03= QD5(2,L)*F03                                           
         F52L15= QD5(2,L)*F15                                           
         F52L1H= F52L15  *F07                                           
         F53L03= QD5(3,L)*F03                                           
         F53L15= QD5(3,L)*F15                                           
         F53L4H= QD5(3,L)*F420                                          
C                                                                       
         F61K03= QD6(1,K)*F03                                           
         F61K06= F61K03+F61K03                                          
         F61K10= QD6(1,K)*F10                                           
         F61K15= QD6(1,K)*F15                                           
         F61K45= F61K15  *F03                                           
         F61K1H= F61K15  *F07                                           
         F61K2H= QD6(1,K)*F210                                          
         F62K03= QD6(2,K)*F03                                           
         F62K06= F62K03+F62K03                                          
         F62K09= F62K03+F62K06                                          
         F62K10= QD6(2,K)*F10                                           
         F62K15= QD6(2,K)*F15                                           
         F62K45= F62K15  *F03                                           
         F62K1H= F62K15  *F07                                           
         F63K03= QD6(3,K)*F03                                           
         F63K06= F63K03+F63K03                                          
         F63K15= QD6(3,K)*F15                                           
         F63K45= F63K15  *F03                                           
         F64K03= QD6(4,K)*F03                                           
         F64K15= QD6(4,K)*F15                                           
         F64K1H= F64K15  *F07                                           
         F65K03= QD6(5,K)*F03                                           
         F65K2H= QD6(5,K)*F210                                          
C                                                                       
         F71J06= QD7(1,J)*F06                                           
         F71J10= QD7(1,J)*F10                                           
         F71J15= QD7(1,J)*F15                                           
         F71J21= QD7(1,J)*F21                                           
         F71J28= QD7(1,J)*F28                                           
         F72J03= QD7(2,J)*F03                                           
         F72J06= F72J03+F72J03                                          
         F72J10= QD7(2,J)*F10                                           
         F72J15= QD7(2,J)*F15                                           
         F72J21= QD7(2,J)*F21                                           
         F73J03= QD7(3,J)*F03                                           
         F73J06= F73J03+F73J03                                          
         F73J10= QD7(3,J)*F10                                           
         F73J15= QD7(3,J)*F15                                           
         F74J03= QD7(4,J)*F03                                           
         F74J06= F74J03+F74J03                                          
         F74J10= QD7(4,J)*F10                                           
         F75J03= QD7(5,J)*F03                                           
         F75J06= F75J03+F75J03                                          
         F75J15= QD7(5,J)*F15                                           
         F76J03= QD7(6,J)*F03                                           
         F76J21= QD7(6,J)*F21                                           
         F77J28= QD7(7,J)*F28                                           
C                                                                       
         A11   =     QD5(1,L)*QX-QD4(1,M)                               
         R11   =     QD5(1,L)*CY-QD4(1,M)                               
C                                                                       
         B13   =     F61K03  *QX-F51L09                                 
         B22   =     F62K03  *QX-F52L03                                 
         B23   =     F62K15  *QX-F52L15  *F03                           
         B33   =     F63K03  *QX-F53L03                                 
         S11   =     F61K03  *CY-F51L03                                 
         S13   =     F61K03  *CY-F51L09                                 
         S22   =     F62K03  *CY-F52L03                                 
         S23   =     F62K03  *CY-F52L03  *F03                           
         S33   =     F63K03  *CY-F53L03                                 
C                                                                       
         C33   =     F73J06  *QX-F63K03  *F06                           
         C44   =     QD7(4,J)*QX-QD6(4,K)                               
         C4T   =     F74J10  *QX-F64K03  *F10                           
         C55   =     QD7(5,J)*QX-QD6(5,K)                               
         T13   =     QD7(1,J)*CY-F61K03                                 
         T22   =     F72J03  *CY-F62K03                                 
         T23   =     F72J03  *CY-F62K09                                 
         T33   =     QD7(3,J)*CY-QD6(3,K)                               
         T3T   =     QD7(3,J)*CY-F63K03                                 
         T44   =     QD7(4,J)*CY-QD6(4,K)                               
         T43   =     QD7(4,J)*CY-F64K03                                 
         T4D   =     F74J10  *CY-F64K03  *F10                           
         T55   =     QD7(5,J)*CY-QD6(5,K)                               
C                                                                       
         D11   =    (QD6(1,K)*QX-F51L06  )*QX+F41M03                    
         U11   =    (QD6(1,K)*CY-F51L06  )*CY+F41M03                    
C                                                                       
         E16   =    (QD7(1,J)*QX-F61K06  )*QX+F51L03                    
         E1T   =    (QD7(1,J)*QX-F61K10  )*QX+F51L15                    
         E26   =   ((QD7(2,J)*QX-F62K06  )*QX+F52L03)*F03               
         E2T   =    (QD7(2,J)*QX-F62K10  )*QX+F52L15                    
         E36   =    (QD7(3,J)*QX-F63K06  )*QX+F53L03                    
         V16   =    (QD7(1,J)*CY-F61K06  )*CY+F51L03                    
         V1T   =    (QD7(1,J)*CY-F61K10  )*CY+F51L15                    
         V26   =   ((QD7(2,J)*CY-F62K06  )*CY+F52L03)*F03               
         V2T   =    (QD7(2,J)*CY-F62K10  )*CY+F52L15                    
         V36   =    (QD7(3,J)*CY-F63K06  )*CY+F53L03                    
C                                                                       
         G11   =   ((QD7(1,J)*QX-F61K15  )*QX+F51L45)*QX-F41M15         
         W11   =   ((QD7(1,J)*CY-F61K15  )*CY+F51L45)*CY-F41M15         
C                                                                       
         WRK( 1,I)=((QD8(1,I)*QX-F71J28  )*QX+F61K2H)*QX-F51L4H         
         WRK( 1,I)=              WRK( 1,I)*QX+F41M1H                    
         WRK( 2,I)=((QD8(1,I)*QX-F71J21  )*QX+F61K1H)*QX-F51L1H         
         WRK( 3,I)=((QD8(2,I)*QX-F72J21  )*QX+F62K1H)*QX-F52L1H         
         WRK( 4,I)=((QD8(1,I)*QX-F71J15  )*QX+F61K45)*QX-F51L15         
         WRK( 4,I)=              WRK( 4,I)*CY-G11                       
         WRK( 5,I)=((QD8(2,I)*QX-F72J15  )*QX+F62K45)*QX-F52L15         
         WRK( 6,I)=((QD8(3,I)*QX-F73J15  )*QX+F63K45)*QX-F53L15 -G11    
         WRK( 7,I)=((QD8(1,I)*QX-F71J10  )*QX+F61K15)*CY -E1T*F03       
         WRK( 8,I)=((QD8(2,I)*QX-F72J10  )*QX+F62K15)*CY -E2T           
         WRK( 9,I)= (QD8(3,I)*QX-F73J10  )*QX+F63K15     -E1T           
         WRK(10,I)= (QD8(4,I)*QX-F74J10  )*QX+F64K15     -E2T*F03       
         WRK(11,I)= (QD8(1,I)*CY-F71J06  )*CY+F61K03                    
         WRK(11,I)=             (WRK(11,I)*QX-V16*F06)*QX+U11*F03       
         WRK(12,I)=((QD8(2,I)*QX-F72J06  )*QX+F62K03)*CY -E26           
         WRK(13,I)=((QD8(3,I)*QX-F73J06  )*QX+F63K03-E16)*CY-E36+D11    
         WRK(14,I)= (QD8(4,I)*QX-F74J06  )*QX+F64K03-E26                
         WRK(15,I)= (QD8(5,I)*QX-F75J06  )*QX+F65K03-E36*F06+D11*F03    
         WRK(16,I)=((QD8(1,I)*CY-F71J10  )*CY+F61K15)*QX-V1T*F03        
         WRK(17,I)=((QD8(2,I)*CY-F72J06  )*CY+F62K03)*QX-V26            
         WRK(18,I)= (QD8(3,I)*CY-F73J03  -T13)*QX-T3T*F03+S13           
         WRK(19,I)= (QD8(4,I)*CY-QD7(4,J)-T22)*QX-T44*F03+S22*F03       
         WRK(20,I)=  QD8(5,I)*QX-F75J03          -C33+B13               
         WRK(21,I)=  QD8(6,I)*QX-F76J03          -C4T+B23               
         WRK(22,I)=((QD8(1,I)*CY-F71J15  )*CY+F61K45)*CY-F51L15         
         WRK(22,I)=              WRK(22,I)*QX-W11                       
         WRK(23,I)=((QD8(2,I)*CY-F72J10  )*CY+F62K15)*QX -V2T           
         WRK(24,I)=((QD8(3,I)*CY-F73J06  )*CY+F63K03-V16)*QX-V36+U11    
         WRK(25,I)= (QD8(4,I)*CY-F74J03  -T23)*QX-T43+S23               
         WRK(26,I)=  QD8(5,I)*CY-QD7(5,J)-T33*F06+S11                   
         WRK(26,I)=              WRK(26,I)*QX-(T55-S33-S33+R11*F03)     
         WRK(27,I)=  QD8(6,I)*QX-QD7(6,J)    -(C44*F10-B22*F05)         
         WRK(28,I)=  QD8(7,I)*QX-QD7(7,J)    -(C55-B33+A11)*F15         
         WRK(29,I)=((QD8(1,I)*CY-F71J21  )*CY+F61K1H)*CY-F51L1H         
         WRK(30,I)=((QD8(2,I)*CY-F72J15  )*CY+F62K45)*CY-F52L15         
         WRK(31,I)= (QD8(3,I)*CY-F73J10  )*CY+F63K15-V1T                
         WRK(32,I)= (QD8(4,I)*CY-F74J06  )*CY+F64K03-V26                
         WRK(33,I)=  QD8(5,I)*CY-F75J03      -T3T*F06+S13               
         WRK(34,I)=  QD8(6,I)*CY-QD7(6,J)    -T44*F10+S22*F05           
         WRK(35,I)=  QD8(7,I)   -F75J15  +F63K45-F51L15                 
         WRK(36,I)=  QD8(8,I)   -F76J21  +F64K1H-F52L1H                 
         WRK(37,I)=((QD8(1,I)*CY-F71J28  )*CY+F61K2H)*CY-F51L4H         
         WRK(37,I)=              WRK(37,I)*CY+F41M1H                    
         WRK(38,I)=((QD8(2,I)*CY-F72J21  )*CY+F62K1H)*CY-F52L1H         
         WRK(39,I)=((QD8(3,I)*CY-F73J15  )*CY+F63K45)*CY-F53L15-W11     
         WRK(40,I)= (QD8(4,I)*CY-F74J10  )*CY+F64K15   -V2T*F03         
         WRK(41,I)= (QD8(5,I)*CY-F75J06  )*CY+F65K03   -V36*F06+U11*F03 
         WRK(42,I)=  QD8(6,I)*CY-F76J03                -T4D+S23*F05     
         WRK(43,I)=  QD8(7,I)*CY-QD7(7,J)         -(T55-S33+R11)*F15    
         WRK(44,I)=  QD8(8,I)   -F76J21      +F64K1H   -F52L1H          
         WRK(45,I)=  QD8(9,I)   -F77J28      +F65K2H   -F53L4H+F41M1H   
      ENDDO                                                             
C                                                                       
      RETURN                                                            
      END                                                               
