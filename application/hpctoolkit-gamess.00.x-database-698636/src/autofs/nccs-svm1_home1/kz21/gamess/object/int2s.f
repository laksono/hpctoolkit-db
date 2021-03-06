C  8 Aug 13 - JMS - coalesce int2s-2x explicit expressions in one file  
C 10 NOV 04 - KI  - QUADRATURE ROUTINES FOR S,P,D,L ROTATED AXIS CODE   
C     The original program by Kazuya Ishimura in 2004 contained         
C     "explicit expressions" for each kind of integral, comprising      
C     about 28,000 lines of code, divided in files named "s" to "x".    
C     Jose Sierra reorganized this into a looped version, in just       
C     this one file "s", using about 7,000 lines in August 2013.        
C     Execution times are essentially unchanged, but compile time       
C     is greatly improved.                                              
C                                                                       
C*MODULE INT2S   *DECK MCDV07                                           
C>                                                                      
C>    @brief   DSSS case                                                
C>                                                                      
C>    @details integration of a DSSS case                               
C>             Simplified calculation of F(I,J,K,L) for cases where     
C>                 I = 1..6,  J = 1..1,  K = 1..KX,  and  L = 1..LX     
C>             using auxiliary arrays C1, C2 and C3.                    
C>                                                                      
      SUBROUTINE MCDV07(F,QX,QZ)                                        
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      DIMENSION F(6,6,6,6)                                              
C                                                                       
      COMMON /KI4   / R0(   25),R1( 3,40),R2( 6,56),R3(10,52),          
     *                R4(15,42),R5(21,24),R6(28,12),R7(36, 4),R8(   45) 
C$omp threadprivate(/KI4/)
C                                                                       
      PARAMETER (KX=1)                                                  
      PARAMETER (LX=1)                                                  
      DIMENSION       C1(  2,KX,LX),C2(  3,KX,LX),C3(  6,KX,LX)         
C                                                                       
      DIMENSION IN6(6)                                                  
      DATA IN6/ 1, 4, 5, 2, 6, 3/                                       
C                                                                       
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..1,  K = 1  and  L = 1                   
C                                                                       
      DO J=1,6                                                          
         M= IN6(J)                                                      
         C3(J,1,1)=+R2(M,1)                                             
      ENDDO                                                             
C                                                                       
      C2(1,1,1)=+R1(1,1)                                                
C                                                                       
      C2(2,1,1)=+R1(2,1)                                                
      C2(3,1,1)=+R1(3,1)                                                
C                                                                       
      C1(1,1,1)=+R0(  1)                                                
      C1(2,1,1)=+R0(  2)                                                
C                                                                       
CC    DO L=1,LX                                                         
CC       DO K=1,KX                                                      
C                                                                       
         L=1                                                            
            K=1                                                         
C                                                                       
            F(1,1,K,L)=+C3(  1,K,L)+C1(  1,K,L)                         
     *               +(+C2(  1,K,L)+C2(  1,K,L)+C1(  2,K,L)*QX)*QX      
            F(2,1,K,L)=+C3(  4,K,L)+C1(  1,K,L)                         
            F(3,1,K,L)=+C3(  6,K,L)+C1(  1,K,L)                         
     *               +(+C2(  3,K,L)+C2(  3,K,L)+C1(  2,K,L)*QZ)*QZ      
            F(4,1,K,L)=+C3(  2,K,L)+C2(  2,K,L)*QX                      
            F(5,1,K,L)=+C3(  3,K,L)+C2(  3,K,L)*QX                      
     *               +(+C2(  1,K,L)+C1(  2,K,L)*QX)*QZ                  
            F(6,1,K,L)=+C3(  5,K,L)+C2(  2,K,L)*QZ                      
C                                                                       
CC       ENDDO                                                          
CC    ENDDO                                                             
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2S   *DECK MCDV08                                           
C>                                                                      
C>    @brief   DPSS case                                                
C>                                                                      
C>    @details integration of a DPSS case                               
C>                                                                      
      SUBROUTINE MCDV08(F,QX,QZ)                                        
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      DIMENSION F(6,6,6,6)                                              
C                                                                       
      COMMON /KI4   / R0(   25),R1( 3,40),R2( 6,56),R3(10,52),          
     *                R4(15,42),R5(21,24),R6(28,12),R7(36, 4),R8(   45) 
C$omp threadprivate(/KI4/)
C                                                                       
      DIMENSION                 S1( 4, 3),S2( 3, 6)                     
Cjms                                                                    
C     Simplified calculation of F(I,J,K,L) for cases                    
C     where:  I = 1..6,  J = 1..4,  K = 1..KX  and  L = 1..LX           
C     using auxiliary arrays D1, D2, D3 and D4                          
Cjms                                                                    
      PARAMETER (KX=1)                                                  
      PARAMETER (LX=1)                                                  
      DIMENSION       D1(  5,KX,LX),D2(4,3,KX,LX),D3(3,6,KX,LX)         
      DIMENSION       D4( 10,KX,LX)                                     
C                                                                       
      DIMENSION IN6(6)                                                  
      DATA IN6/ 1, 4, 5, 2, 6, 3/                                       
C                                                                       
      DO 101 J=1, 4                                                     
         DO 101 I=1, 3                                                  
  101 S1(J,I)= R1(I,J)                                                  
      DO 102 J=1, 3                                                     
         DO 102 I=1, 6                                                  
  102 S2(J,I)= R2(I,J)                                                  
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..4,  K = 1  and  L = 1                   
Cjms                                                                    
      DO J=1,10                                                         
         D4(  J,1,1)=+R3(J,1)                                           
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         M= IN6(J)                                                      
         D3(1,J,1,1)=+S2(1,M)                                           
         D3(2,J,1,1)=+S2(2,M)                                           
         D3(3,J,1,1)=+S2(3,M)                                           
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         D2(1,J,1,1)=+S1(1,J)                                           
         D2(2,J,1,1)=+S1(2,J)                                           
         D2(3,J,1,1)=+S1(3,J)                                           
         D2(4,J,1,1)=+S1(4,J)                                           
      ENDDO                                                             
C                                                                       
      DO I=1,5                                                          
         D1(I  ,1,1)=+R0(I)                                             
      ENDDO                                                             
C                                                                       
      CALL FIJKL4(KX,LX,F,QX,QZ,D1,D2,D3,D4)                            
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2S   *DECK MCDV09                                           
C>                                                                      
C>    @brief   DSPS case                                                
C>                                                                      
C>    @details integration of a DSPS case                               
C>                                                                      
      SUBROUTINE MCDV09(F,QX,QZ)                                        
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      DIMENSION F(6,6,6,6)                                              
C                                                                       
      COMMON /KI4   / R0(   25),R1( 3,40),R2( 6,56),R3(10,52),          
     *                R4(15,42),R5(21,24),R6(28,12),R7(36, 4),R8(   45) 
C$omp threadprivate(/KI4/)
C                                                                       
      DIMENSION                 S1( 4, 3)                               
Cjms                                                                    
C     Simplified calculation of F(I,J,K,L) for cases                    
C     where:  I = 1..6,  J = 1..1,  K = 1..KX  and  L = 1..LX           
C     using auxiliary arrays C1, C2 and C3                              
Cjms                                                                    
      PARAMETER (KX=4)                                                  
      PARAMETER (LX=1)                                                  
      DIMENSION       C1(  2,KX,LX),C2(  3,KX,LX),C3(  6,KX,LX)         
C                                                                       
      DIMENSION INS(KX,LX),INI(KX,LX),IND( 6,KX,LX)                     
      DIMENSION IN6(6)                                                  
      DATA INS/ 0, 0, 0, 1/                                             
      DATA INI/ 0, 0, 1, 1/                                             
      DATA IN6/ 1, 4, 5, 2, 6, 3/                                       
C                                                                       
      DO L=1,LX                                                         
         DO K=1,KX                                                      
            IJ= 0                                                       
            IN= INS(K,L)                                                
            DO I=1,3                                                    
               IN= IN+INI(K,L)                                          
               DO J=1,I                                                 
                  IJ= IJ+1                                              
                  IN= IN+1                                              
                  IND(IJ,K,L)= IN                                       
               ENDDO                                                    
            ENDDO                                                       
         ENDDO                                                          
      ENDDO                                                             
      DO 101 J=1, 4                                                     
         DO 101 I=1, 3                                                  
  101 S1(J,I)= R1(I,J)                                                  
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..1,  K = 1  and  L = 1                   
Cjms                                                                    
      DO J=1,6                                                          
         M= IN6(J)                                                      
         C3(J,1,1)=+R2(M, 1)                                            
      ENDDO                                                             
C                                                                       
      C2(  1,1,1)=+R1( 1, 1)                                            
C                                                                       
      C2(  2,1,1)=+R1( 2, 1)                                            
      C2(  3,1,1)=+R1( 3, 1)                                            
C                                                                       
      C1(  1,1,1)=+R0(    1)                                            
      C1(  2,1,1)=+R0(    2)                                            
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..1,  K = 2  and  L = 1                   
Cjms                                                                    
      DO J=1,6                                                          
         C3(J,2,1)=-R3(J, 1)                                            
      ENDDO                                                             
C                                                                       
      C2(  1,2,1)=-R2( 1, 3)                                            
C                                                                       
      C2(  2,2,1)=-R2( 4, 3)                                            
      C2(  3,2,1)=-R2( 5, 3)                                            
C                                                                       
      C1(  1,2,1)=-S1( 3, 1)                                            
      C1(  2,2,1)=-S1( 4, 1)                                            
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..1,  K = 3  and  L = 1                   
Cjms                                                                    
      DO J=1,6                                                          
         K= IND(J,3,1)                                                  
         C3(J,3,1)=-R3(K, 1)                                            
      ENDDO                                                             
C                                                                       
      C2(  1,3,1)=-R2( 4, 3)                                            
C                                                                       
      C2(  2,3,1)=-R2( 2, 3)                                            
      C2(  3,3,1)=-R2( 6, 3)                                            
C                                                                       
      C1(  1,3,1)=-S1( 3, 2)                                            
      C1(  2,3,1)=-S1( 4, 2)                                            
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..1,  K = 4  and  L = 1                   
Cjms                                                                    
      DO J=1,6                                                          
         K= IND(J,4,1)                                                  
         M= IN6(J)                                                      
         C3(J,4,1)=-R3(K, 1)+R2( M, 2)                                  
      ENDDO                                                             
C                                                                       
      C2(  1,4,1)=-R2( 5, 3)+R1( 1, 2)                                  
C                                                                       
      C2(  2,4,1)=-R2( 6, 3)+R1( 2, 2)                                  
      C2(  3,4,1)=-R2( 3, 3)+R1( 3, 2)                                  
C                                                                       
      C1(  1,4,1)=-S1( 3, 3)+R0(    3)                                  
      C1(  2,4,1)=-S1( 4, 3)+R0(    4)                                  
C                                                                       
      DO L=1,LX                                                         
         DO K=1,KX                                                      
C                                                                       
            F(1,1,K,L)=+C3(  1,K,L)+C1(  1,K,L)                         
     *               +(+C2(  1,K,L)+C2(  1,K,L)+C1(  2,K,L)*QX)*QX      
            F(2,1,K,L)=+C3(  4,K,L)+C1(  1,K,L)                         
            F(3,1,K,L)=+C3(  6,K,L)+C1(  1,K,L)                         
     *               +(+C2(  3,K,L)+C2(  3,K,L)+C1(  2,K,L)*QZ)*QZ      
            F(4,1,K,L)=+C3(  2,K,L)+C2(  2,K,L)*QX                      
            F(5,1,K,L)=+C3(  3,K,L)+C2(  3,K,L)*QX                      
     *               +(+C2(  1,K,L)+C1(  2,K,L)*QX)*QZ                  
            F(6,1,K,L)=+C3(  5,K,L)+C2(  2,K,L)*QZ                      
C                                                                       
         ENDDO                                                          
      ENDDO                                                             
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2S   *DECK MCDV10                                           
C>                                                                      
C>    @brief   DDSS case                                                
C>                                                                      
C>    @details integration of a DDSS case                               
C>                                                                      
      SUBROUTINE MCDV10(F,QX,QZ)                                        
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      DIMENSION F(6,6,6,6)                                              
C                                                                       
      COMMON /KI4   / R0(   25),R1( 3,40),R2( 6,56),R3(10,52),          
     *                R4(15,42),R5(21,24),R6(28,12),R7(36, 4),R8(   45) 
C$omp threadprivate(/KI4/)
C                                                                       
      DIMENSION                 S1( 4, 3),S2( 4, 6),S3( 2,10)           
Cjms                                                                    
C     Simplified calculation of F(I,J,K,L) for cases                    
C     where:  I = 1..6,  J = 1..6,  K = 1..KX  and  L = 1..LX           
C     using auxiliary arrays E1, E2, E3, E4 and E5                      
Cjms                                                                    
      PARAMETER (KX=1)                                                  
      PARAMETER (LX=1)                                                  
      DIMENSION       E1(  5,KX,LX),E2(4,3,KX,LX),E3(4,6,KX,LX)         
      DIMENSION       E4(2,10,KX,LX),E5( 15,KX,LX)                      
C                                                                       
      DIMENSION IN6(6)                                                  
      DATA IN6/ 1, 4, 5, 2, 6, 3/                                       
C                                                                       
      DO 101 J=1, 4                                                     
         DO 101 I=1, 3                                                  
  101 S1(J,I)= R1(I,J)                                                  
      DO 102 J=1, 4                                                     
         DO 102 I=1, 6                                                  
  102 S2(J,I)= R2(I,J)                                                  
      DO 103 J=1, 2                                                     
         DO 103 I=1,10                                                  
  103 S3(J,I)= R3(I,J)                                                  
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 1  and  L = 1                   
Cjms                                                                    
      DO J=1,15                                                         
         E5(  J,1,1)=+R4(J,1)                                           
      ENDDO                                                             
C                                                                       
      DO J=1,10                                                         
         E4(1,J,1,1)=+S3(1,J)                                           
         E4(2,J,1,1)=+S3(2,J)                                           
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         M= IN6(J)                                                      
         E3(1,J,1,1)=+S2(1,M)                                           
         E3(2,J,1,1)=+S2(2,M)                                           
         E3(3,J,1,1)=+S2(3,M)                                           
         E3(4,J,1,1)=+S2(4,M)                                           
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         E2(1,J,1,1)=+S1(1,J)                                           
         E2(2,J,1,1)=+S1(2,J)                                           
         E2(3,J,1,1)=+S1(3,J)                                           
         E2(4,J,1,1)=+S1(4,J)                                           
      ENDDO                                                             
C                                                                       
      DO I=1,5                                                          
         E1(I  ,1,1)=+R0(I)                                             
      ENDDO                                                             
C                                                                       
      CALL FIJKL5(KX,LX,F,QX,QZ,E1,E2,E3,E4,E5)                         
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2S   *DECK MCDV11                                           
C>                                                                      
C>    @brief   DPPS case                                                
C>                                                                      
C>    @details integration of a DPPS case                               
C>                                                                      
      SUBROUTINE MCDV11(F,QX,QZ)                                        
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      DIMENSION F(6,6,6,6)                                              
C                                                                       
      COMMON /KI4   / R0(   25),R1( 3,40),R2( 6,56),R3(10,52),          
     *                R4(15,42),R5(21,24),R6(28,12),R7(36, 4),R8(   45) 
C$omp threadprivate(/KI4/)
C                                                                       
      DIMENSION                 S1(13, 3),S2(10, 6),S3( 5,10)           
Cjms                                                                    
C     Simplified calculation of F(I,J,K,L) for cases                    
C     where:  I = 1..6,  J = 1..4,  K = 1..KX  and  L = 1..LX           
C     using auxiliary arrays D1, D2, D3 and D4                          
Cjms                                                                    
      PARAMETER (KX=4)                                                  
      PARAMETER (LX=1)                                                  
      DIMENSION       D1(  5,KX,LX),D2(4,3,KX,LX),D3(3,6,KX,LX)         
      DIMENSION       D4( 10,KX,LX)                                     
C                                                                       
      DIMENSION INS(KX,LX),INI(KX,LX),IND(10,KX,LX)                     
      DIMENSION IN6(6)                                                  
      DATA INS/ 0, 0, 0, 1/                                             
      DATA INI/ 0, 0, 1, 1/                                             
      DATA IN6/ 1, 4, 5, 2, 6, 3/                                       
C                                                                       
      DO L=1,LX                                                         
         DO K=1,KX                                                      
            IJ= 0                                                       
            IN= INS(K,L)                                                
            DO I=1,4                                                    
               IN= IN+INI(K,L)                                          
               DO J=1,I                                                 
                  IJ= IJ+1                                              
                  IN= IN+1                                              
                  IND(IJ,K,L)= IN                                       
               ENDDO                                                    
            ENDDO                                                       
         ENDDO                                                          
      ENDDO                                                             
      DO 101 J=1,13                                                     
         DO 101 I=1, 3                                                  
  101 S1(J,I)= R1(I,J)                                                  
      DO 102 J=1,10                                                     
         DO 102 I=1, 6                                                  
  102 S2(J,I)= R2(I,J)                                                  
      DO 103 J=1, 5                                                     
         DO 103 I=1,10                                                  
  103 S3(J,I)= R3(I,J)                                                  
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..4,  K = 1  and  L = 1                   
Cjms                                                                    
      DO J=1,10                                                         
         D4(  J,1,1)=+R3( J,1)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         M= IN6(J)                                                      
         D3(1,J,1,1)=+S2( 1,M)                                          
         D3(2,J,1,1)=+S2( 2,M)                                          
         D3(3,J,1,1)=+S2( 3,M)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         D2(1,J,1,1)=+S1( 1,J)                                          
         D2(2,J,1,1)=+S1( 2,J)                                          
         D2(3,J,1,1)=+S1( 3,J)                                          
         D2(4,J,1,1)=+S1( 4,J)                                          
      ENDDO                                                             
C                                                                       
      DO I=1,5                                                          
         D1(I  ,1,1)=+R0(I)                                             
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..4,  K = 2  and  L = 1                   
Cjms                                                                    
      DO J=1,10                                                         
         D4(  J,2,1)=-R4( J,1)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         D3(1,J,2,1)=-S3( 3,J)                                          
         D3(2,J,2,1)=-S3( 4,J)                                          
         D3(3,J,2,1)=-S3( 5,J)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         M= IN6(J)                                                      
         D2(1,J,2,1)=-S2( 7,M)                                          
         D2(2,J,2,1)=-S2( 8,M)                                          
         D2(3,J,2,1)=-S2( 9,M)                                          
         D2(4,J,2,1)=-S2(10,M)                                          
      ENDDO                                                             
C                                                                       
         J=1                                                            
      DO I=1,5                                                          
         D1(I  ,2,1)=-S1(I+ 8,J)                                        
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..4,  K = 3  and  L = 1                   
Cjms                                                                    
      DO J=1,10                                                         
         K= IND(J,3,1)                                                  
         D4(  J,3,1)=-R4( K,1)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         K= IND(J,3,1)                                                  
         D3(1,J,3,1)=-S3( 3,K)                                          
         D3(2,J,3,1)=-S3( 4,K)                                          
         D3(3,J,3,1)=-S3( 5,K)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         K= IND(J,3,1)                                                  
         M= IN6(K)                                                      
         D2(1,J,3,1)=-S2( 7,M)                                          
         D2(2,J,3,1)=-S2( 8,M)                                          
         D2(3,J,3,1)=-S2( 9,M)                                          
         D2(4,J,3,1)=-S2(10,M)                                          
      ENDDO                                                             
C                                                                       
         J=1                                                            
         K= IND(J,3,1)                                                  
      DO I=1,5                                                          
         D1(I  ,3,1)=-S1(I+ 8,K)                                        
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..4,  K = 4  and  L = 1                   
Cjms                                                                    
      DO J=1,10                                                         
         K= IND(J,4,1)                                                  
         D4(  J,4,1)=-R4( K,1)+R3( J,2)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         K= IND(J,4,1)                                                  
         M= IN6(J)                                                      
         D3(1,J,4,1)=-S3( 3,K)+S2( 4,M)                                 
         D3(2,J,4,1)=-S3( 4,K)+S2( 5,M)                                 
         D3(3,J,4,1)=-S3( 5,K)+S2( 6,M)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         K= IND(J,4,1)                                                  
         M= IN6(K)                                                      
         D2(1,J,4,1)=-S2( 7,M)+S1( 5,J)                                 
         D2(2,J,4,1)=-S2( 8,M)+S1( 6,J)                                 
         D2(3,J,4,1)=-S2( 9,M)+S1( 7,J)                                 
         D2(4,J,4,1)=-S2(10,M)+S1( 8,J)                                 
      ENDDO                                                             
C                                                                       
         J=1                                                            
         K= IND(J,4,1)                                                  
      DO I=1,5                                                          
         D1(I  ,4,1)=-S1(I+ 8,K)+R0(I+ 5)                               
      ENDDO                                                             
C                                                                       
      CALL FIJKL4(KX,LX,F,QX,QZ,D1,D2,D3,D4)                            
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2S   *DECK MCDV12                                           
C>                                                                      
C>    @brief   DSDS case                                                
C>                                                                      
C>    @details integration of a DSDS case                               
C>                                                                      
      SUBROUTINE MCDV12(F,QX,QZ)                                        
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      DIMENSION F(6,6,6,6)                                              
C                                                                       
      COMMON /KI4   / R0(   25),R1( 3,40),R2( 6,56),R3(10,52),          
     *                R4(15,42),R5(21,24),R6(28,12),R7(36, 4),R8(   45) 
C$omp threadprivate(/KI4/)
C                                                                       
      DIMENSION                 S1( 4, 3),S2( 5, 6)                     
C                                                                       
      PARAMETER (F02=2.0D+00)                                           
Cjms                                                                    
C     Simplified calculation of F(I,J,K,L) for cases                    
C     where:  I = 1..6,  J = 1..1,  K = 1..KX  and  L = 1..LX           
C     using auxiliary arrays C1, C2 and C3                              
Cjms                                                                    
      PARAMETER (KX=6)                                                  
      PARAMETER (LX=1)                                                  
      DIMENSION       C1(  2,KX,LX),C2(  3,KX,LX),C3(  6,KX,LX)         
C                                                                       
      DIMENSION INS(KX,LX),INI(KX,LX),IND( 6,KX,LX)                     
      DIMENSION IN6(6)                                                  
      DATA INS/ 0, 1, 3, 0, 1, 2/                                       
      DATA INI/ 0, 2, 2, 1, 1, 2/                                       
      DATA IN6/ 1, 4, 5, 2, 6, 3/                                       
C                                                                       
      DO L=1,LX                                                         
         DO K=1,KX                                                      
            IJ= 0                                                       
            IN= INS(K,L)                                                
            DO I=1,3                                                    
               IN= IN+INI(K,L)                                          
               DO J=1,I                                                 
                  IJ= IJ+1                                              
                  IN= IN+1                                              
                  IND(IJ,K,L)= IN                                       
               ENDDO                                                    
            ENDDO                                                       
         ENDDO                                                          
      ENDDO                                                             
      DO 101 J=1, 4                                                     
         DO 101 I=1, 3                                                  
  101 S1(J,I)= R1(I,J)                                                  
      DO 102 J=1, 5                                                     
         DO 102 I=1, 6                                                  
  102 S2(J,I)= R2(I,J)                                                  
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..1,  K = 1  and  L = 1                   
Cjms                                                                    
      DO J=1,6                                                          
         M= IN6(J)                                                      
         C3(J,1,1)=+R4( J, 1)+R2( M, 1)                                 
      ENDDO                                                             
C                                                                       
      C2(  1,1,1)=+R3( 1, 2)+R1( 1, 1)                                  
C                                                                       
      C2(  2,1,1)=+R3( 2, 2)+R1( 2, 1)                                  
      C2(  3,1,1)=+R3( 3, 2)+R1( 3, 1)                                  
C                                                                       
      C1(  1,1,1)=+S2( 4, 1)+R0(    1)                                  
      C1(  2,1,1)=+S2( 5, 1)+R0(    2)                                  
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..1,  K = 2  and  L = 1                   
Cjms                                                                    
      DO J=1,6                                                          
         K= IND(J,2,1)                                                  
         M= IN6(J)                                                      
         C3(J,2,1)=+R4( K, 1)+R2( M, 1)                                 
      ENDDO                                                             
C                                                                       
      C2(  1,2,1)=+R3( 4, 2)+R1( 1, 1)                                  
C                                                                       
      C2(  2,2,1)=+R3( 7, 2)+R1( 2, 1)                                  
      C2(  3,2,1)=+R3( 8, 2)+R1( 3, 1)                                  
C                                                                       
      C1(  1,2,1)=+S2( 4, 2)+R0(    1)                                  
      C1(  2,2,1)=+S2( 5, 2)+R0(    2)                                  
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..1,  K = 3  and  L = 1                   
Cjms                                                                    
      C3(  1,3,1)=+R4( 6, 1)-R3( 3, 1)*F02+R2( 1, 1)+R2( 1, 2)          
C                                                                       
      C3(  2,3,1)=+R4( 9, 1)-R3( 5, 1)*F02+R2( 4, 1)+R2( 4, 2)          
      C3(  3,3,1)=+R4(10, 1)-R3( 6, 1)*F02+R2( 5, 1)+R2( 5, 2)          
C                                                                       
      C3(  4,3,1)=+R4(13, 1)-R3( 8, 1)*F02+R2( 2, 1)+R2( 2, 2)          
      C3(  5,3,1)=+R4(14, 1)-R3( 9, 1)*F02+R2( 6, 1)+R2( 6, 2)          
      C3(  6,3,1)=+R4(15, 1)-R3(10, 1)*F02+R2( 3, 1)+R2( 3, 2)          
C                                                                       
      C2(  1,3,1)=+R3( 6, 2)-R2( 5, 3)*F02+R1( 1, 1)+R1( 1, 2)          
C                                                                       
      C2(  2,3,1)=+R3( 9, 2)-R2( 6, 3)*F02+R1( 2, 1)+R1( 2, 2)          
      C2(  3,3,1)=+R3(10, 2)-R2( 3, 3)*F02+R1( 3, 1)+R1( 3, 2)          
C                                                                       
      C1(  1,3,1)=+S2( 4, 3)-S1( 3, 3)*F02+R0(    1)+R0(    3)          
      C1(  2,3,1)=+S2( 5, 3)-S1( 4, 3)*F02+R0(    2)+R0(    4)          
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..1,  K = 4  and  L = 1                   
Cjms                                                                    
      DO J=1,6                                                          
         K= IND(J,4,1)                                                  
         C3(J,4,1)=+R4( K, 1)                                           
      ENDDO                                                             
C                                                                       
      C2(  1,4,1)=+R3( 2, 2)                                            
C                                                                       
      C2(  2,4,1)=+R3( 4, 2)                                            
      C2(  3,4,1)=+R3( 5, 2)                                            
C                                                                       
      C1(  1,4,1)=+S2( 4, 4)                                            
      C1(  2,4,1)=+S2( 5, 4)                                            
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..1,  K = 5  and  L = 1                   
Cjms                                                                    
      DO J=1,6                                                          
         K= IND(J,5,1)                                                  
         C3(J,5,1)=+R4( K, 1)-R3( J, 1)                                 
      ENDDO                                                             
C                                                                       
      C2(  1,5,1)=+R3( 3, 2)-R2( 1, 3)                                  
C                                                                       
      C2(  2,5,1)=+R3( 5, 2)-R2( 4, 3)                                  
      C2(  3,5,1)=+R3( 6, 2)-R2( 5, 3)                                  
C                                                                       
      C1(  1,5,1)=+S2( 4, 5)-S1( 3, 1)                                  
      C1(  2,5,1)=+S2( 5, 5)-S1( 4, 1)                                  
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..1,  K = 6  and  L = 1                   
Cjms                                                                    
      C3(  1,6,1)=+R4( 5, 1)-R3( 2, 1)                                  
C                                                                       
      C3(  2,6,1)=+R4( 8, 1)-R3( 4, 1)                                  
      C3(  3,6,1)=+R4( 9, 1)-R3( 5, 1)                                  
C                                                                       
      C3(  4,6,1)=+R4(12, 1)-R3( 7, 1)                                  
      C3(  5,6,1)=+R4(13, 1)-R3( 8, 1)                                  
      C3(  6,6,1)=+R4(14, 1)-R3( 9, 1)                                  
C                                                                       
      C2(  1,6,1)=+R3( 5, 2)-R2( 4, 3)                                  
C                                                                       
      C2(  2,6,1)=+R3( 8, 2)-R2( 2, 3)                                  
      C2(  3,6,1)=+R3( 9, 2)-R2( 6, 3)                                  
C                                                                       
      C1(  1,6,1)=+S2( 4, 6)-S1( 3, 2)                                  
      C1(  2,6,1)=+S2( 5, 6)-S1( 4, 2)                                  
C                                                                       
      DO L=1,LX                                                         
         DO K=1,KX                                                      
C                                                                       
            F(1,1,K,L)=+C3(  1,K,L)+C1(  1,K,L)                         
     *               +(+C2(  1,K,L)+C2(  1,K,L)+C1(  2,K,L)*QX)*QX      
            F(2,1,K,L)=+C3(  4,K,L)+C1(  1,K,L)                         
            F(3,1,K,L)=+C3(  6,K,L)+C1(  1,K,L)                         
     *               +(+C2(  3,K,L)+C2(  3,K,L)+C1(  2,K,L)*QZ)*QZ      
            F(4,1,K,L)=+C3(  2,K,L)+C2(  2,K,L)*QX                      
            F(5,1,K,L)=+C3(  3,K,L)+C2(  3,K,L)*QX                      
     *               +(+C2(  1,K,L)+C1(  2,K,L)*QX)*QZ                  
            F(6,1,K,L)=+C3(  5,K,L)+C2(  2,K,L)*QZ                      
C                                                                       
         ENDDO                                                          
      ENDDO                                                             
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2S   *DECK MCDV13                                           
C>                                                                      
C>    @brief   DSPP case                                                
C>                                                                      
C>    @details integration of a DSPP case                               
C>                                                                      
      SUBROUTINE MCDV13(F,QX,QZ)                                        
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      DIMENSION F(6,6,6,6)                                              
C                                                                       
      COMMON /KI4   / R0(   25),R1( 3,40),R2( 6,56),R3(10,52),          
     *                R4(15,42),R5(21,24),R6(28,12),R7(36, 4),R8(   45) 
C$omp threadprivate(/KI4/)
C                                                                       
      DIMENSION                 S1(13, 3),S2(11, 6)                     
C                                                                       
      LOGICAL   LSYM13                                                  
Cjms                                                                    
C     Simplified calculation of F(I,J,K,L) for cases                    
C     where:  I = 1..6,  J = 1..1,  K = 1..KX  and  L = 1..LX           
C     using auxiliary arrays C1, C2 and C3                              
Cjms                                                                    
      PARAMETER (KX=4)                                                  
      PARAMETER (LX=4)                                                  
      DIMENSION       C1(  2,KX,LX),C2(  3,KX,LX),C3(  6,KX,LX)         
C                                                                       
      DIMENSION INS(KX,LX),INI(KX,LX),IND( 6,KX,LX)                     
      DIMENSION IN6(6)                                                  
      DATA INS/ 0, 0, 0, 1,  0, 0, 0, 1,  0, 0, 1, 2,  1, 1, 2, 3/      
      DATA INI/ 0, 0, 1, 1,  0, 0, 1, 1,  1, 1, 2, 2,  1, 1, 2, 2/      
      DATA IN6/ 1, 4, 5, 2, 6, 3/                                       
C                                                                       
      DO L=1,LX                                                         
         DO K=1,KX                                                      
            IJ= 0                                                       
            IN= INS(K,L)                                                
            DO I=1,3                                                    
               IN= IN+INI(K,L)                                          
               DO J=1,I                                                 
                  IJ= IJ+1                                              
                  IN= IN+1                                              
                  IND(IJ,K,L)= IN                                       
               ENDDO                                                    
            ENDDO                                                       
         ENDDO                                                          
      ENDDO                                                             
      DO 101 J=1,13                                                     
         DO 101 I=1, 3                                                  
  101 S1(J,I)= R1(I,J)                                                  
      DO 102 J=1,11                                                     
         DO 102 I=1, 6                                                  
  102 S2(J,I)= R2(I,J)                                                  
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..1,  K = 1  and  L = 1                   
Cjms                                                                    
      DO J=1,6                                                          
         M= IN6(J)                                                      
         C3(J,1,1)=+R2( M, 1)                                           
      ENDDO                                                             
C                                                                       
      C2(  1,1,1)=+R1( 1, 1)                                            
C                                                                       
      C2(  2,1,1)=+R1( 2, 1)                                            
      C2(  3,1,1)=+R1( 3, 1)                                            
C                                                                       
      C1(  1,1,1)=+R0(    1)                                            
      C1(  2,1,1)=+R0(    6)                                            
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..1,  K = 2  and  L = 1                   
Cjms                                                                    
      DO J=1,6                                                          
         C3(J,2,1)=-R3( J, 2)                                           
      ENDDO                                                             
C                                                                       
      C2(  1,2,1)=-R2( 1, 6)                                            
C                                                                       
      C2(  2,2,1)=-R2( 4, 6)                                            
      C2(  3,2,1)=-R2( 5, 6)                                            
C                                                                       
      C1(  1,2,1)=-S1( 6, 1)                                            
      C1(  2,2,1)=-S1( 7, 1)                                            
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..1,  K = 3  and  L = 1                   
Cjms                                                                    
      DO J=1,6                                                          
         K= IND(J,3,1)                                                  
         C3(J,3,1)=-R3( K, 2)                                           
      ENDDO                                                             
C                                                                       
      C2(  1,3,1)=-R2( 4, 6)                                            
C                                                                       
      C2(  2,3,1)=-R2( 2, 6)                                            
      C2(  3,3,1)=-R2( 6, 6)                                            
C                                                                       
      C1(  1,3,1)=-S1( 6, 2)                                            
      C1(  2,3,1)=-S1( 7, 2)                                            
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..1,  K = 4  and  L = 1                   
Cjms                                                                    
      DO J=1,6                                                          
         K= IND(J,4,1)                                                  
         M= IN6(J)                                                      
         C3(J,4,1)=-R3( K, 2)+R2( M, 2)                                 
      ENDDO                                                             
C                                                                       
      C2(  1,4,1)=-R2( 5, 6)+R1( 1, 2)                                  
C                                                                       
      C2(  2,4,1)=-R2( 6, 6)+R1( 2, 2)                                  
      C2(  3,4,1)=-R2( 3, 6)+R1( 3, 2)                                  
C                                                                       
      C1(  1,4,1)=-S1( 6, 3)+R0(    2)                                  
      C1(  2,4,1)=-S1( 7, 3)+R0(    7)                                  
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..1,  K = 1  and  L = 2                   
Cjms                                                                    
      DO J=1,6                                                          
         C3(J,1,2)=-R3( J, 3)                                           
      ENDDO                                                             
C                                                                       
      C2(  1,1,2)=-R2( 1, 7)                                            
C                                                                       
      C2(  2,1,2)=-R2( 4, 7)                                            
      C2(  3,1,2)=-R2( 5, 7)                                            
C                                                                       
      C1(  1,1,2)=-S1( 8, 1)                                            
      C1(  2,1,2)=-S1( 9, 1)                                            
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..1,  K = 2  and  L = 2                   
Cjms                                                                    
      DO J=1,6                                                          
         M= IN6(J)                                                      
         C3(J,2,2)=+R4( J, 1)+R2( M, 5)                                 
      ENDDO                                                             
C                                                                       
      C2(  1,2,2)=+R3( 1, 5)+R1( 1, 5)                                  
C                                                                       
      C2(  2,2,2)=+R3( 2, 5)+R1( 2, 5)                                  
      C2(  3,2,2)=+R3( 3, 5)+R1( 3, 5)                                  
C                                                                       
      C1(  1,2,2)=+S2(10, 1)+R0(    5)                                  
      C1(  2,2,2)=+S2(11, 1)+R0(   10)                                  
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..1,  K = 3  and  L = 2                   
Cjms                                                                    
      DO J=1,6                                                          
         K= IND(J,3,2)                                                  
         C3(J,3,2)=+R4( K, 1)                                           
      ENDDO                                                             
C                                                                       
      C2(  1,3,2)=+R3( 2, 5)                                            
C                                                                       
      C2(  2,3,2)=+R3( 4, 5)                                            
      C2(  3,3,2)=+R3( 5, 5)                                            
C                                                                       
      C1(  1,3,2)=+S2(10, 4)                                            
      C1(  2,3,2)=+S2(11, 4)                                            
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..1,  K = 4  and  L = 2                   
Cjms                                                                    
      DO J=1,6                                                          
         K= IND(J,4,2)                                                  
         M= IN6(J)                                                      
         C3(J,4,2)=+R4( K, 1)-R3( J, 4)                                 
      ENDDO                                                             
C                                                                       
      C2(  1,4,2)=+R3( 3, 5)-R2( 1, 8)                                  
C                                                                       
      C2(  2,4,2)=+R3( 5, 5)-R2( 4, 8)                                  
      C2(  3,4,2)=+R3( 6, 5)-R2( 5, 8)                                  
C                                                                       
      C1(  1,4,2)=+S2(10, 5)-S1(10, 1)                                  
      C1(  2,4,2)=+S2(11, 5)-S1(11, 1)                                  
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..1,  K = 1  and  L = 3                   
Cjms                                                                    
      DO J=1,6                                                          
         K= IND(J,1,3)                                                  
         C3(J,1,3)=-R3( K, 3)                                           
      ENDDO                                                             
C                                                                       
      C2(  1,1,3)=-R2( 4, 7)                                            
C                                                                       
      C2(  2,1,3)=-R2( 2, 7)                                            
      C2(  3,1,3)=-R2( 6, 7)                                            
C                                                                       
      C1(  1,1,3)=-S1( 8, 2)                                            
      C1(  2,1,3)=-S1( 9, 2)                                            
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..1,  K = 2  and  L = 3                   
C                                                                       
C     next statements commented because F(I,J,2,3)                      
C     will be obtained at the end from  F(I,J,3,2)                      
Cjms                                                                    
C     DO 323 J=1,6                                                      
C 323 C3(  J,2,3)= C3(  J,3,2)                                          
C        DO 322 I=1,3                                                   
C 322 C2(I,1,2,3)= C2(I,1,3,2)                                          
C        DO 321 I=1,2                                                   
C 321 C1(I  ,2,3)= C1(I  ,3,2)                                          
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..1,  K = 3  and  L = 3                   
Cjms                                                                    
      DO J=1,6                                                          
         K= IND(J,3,3)                                                  
         M= IN6(J)                                                      
         C3(J,3,3)=+R4( K, 1)+R2( M, 5)                                 
      ENDDO                                                             
C                                                                       
      C2(  1,3,3)=+R3( 4, 5)+R1( 1, 5)                                  
C                                                                       
      C2(  2,3,3)=+R3( 7, 5)+R1( 2, 5)                                  
      C2(  3,3,3)=+R3( 8, 5)+R1( 3, 5)                                  
C                                                                       
      C1(  1,3,3)=+S2(10, 2)+R0(    5)                                  
      C1(  2,3,3)=+S2(11, 2)+R0(   10)                                  
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..1,  K = 4  and  L = 3                   
Cjms                                                                    
      C3(  1,4,3)=+R4( 5, 1)-R3( 2, 4)                                  
C                                                                       
      C3(  2,4,3)=+R4( 8, 1)-R3( 4, 4)                                  
      C3(  3,4,3)=+R4( 9, 1)-R3( 5, 4)                                  
C                                                                       
      C3(  4,4,3)=+R4(12, 1)-R3( 7, 4)                                  
      C3(  5,4,3)=+R4(13, 1)-R3( 8, 4)                                  
      C3(  6,4,3)=+R4(14, 1)-R3( 9, 4)                                  
C                                                                       
      C2(  1,4,3)=+R3( 5, 5)-R2( 4, 8)                                  
C                                                                       
      C2(  2,4,3)=+R3( 8, 5)-R2( 2, 8)                                  
      C2(  3,4,3)=+R3( 9, 5)-R2( 6, 8)                                  
C                                                                       
      C1(  1,4,3)=+S2(10, 6)-S1(10, 2)                                  
      C1(  2,4,3)=+S2(11, 6)-S1(11, 2)                                  
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..1,  K = 1  and  L = 4                   
Cjms                                                                    
      DO J=1,6                                                          
         K= IND(J,1,4)                                                  
         M= IN6(J)                                                      
         C3(J,1,4)=-R3( K, 3)+R2( M, 3)                                 
      ENDDO                                                             
C                                                                       
      C2(  1,1,4)=-R2( 5, 7)+R1( 1, 3)                                  
C                                                                       
      C2(  2,1,4)=-R2( 6, 7)+R1( 2, 3)                                  
      C2(  3,1,4)=-R2( 3, 7)+R1( 3, 3)                                  
C                                                                       
      C1(  1,1,4)=-S1( 8, 3)+R0(    3)                                  
      C1(  2,1,4)=-S1( 9, 3)+R0(    8)                                  
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..1,  K = 2  and  L = 4                   
Cjms                                                                    
      DO J=1,6                                                          
         K= IND(J,2,4)                                                  
         C3(J,2,4)=+R4( K, 1)-R3( J, 1)                                 
      ENDDO                                                             
C                                                                       
      C2(  1,2,4)=+R3( 3, 5)-R2( 1, 9)                                  
C                                                                       
      C2(  2,2,4)=+R3( 5, 5)-R2( 4, 9)                                  
      C2(  3,2,4)=+R3( 6, 5)-R2( 5, 9)                                  
C                                                                       
      C1(  1,2,4)=+S2(10, 5)-S1(12, 1)                                  
      C1(  2,2,4)=+S2(11, 5)-S1(13, 1)                                  
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..1,  K = 3  and  L = 4                   
Cjms                                                                    
      C3(  1,3,4)=+R4( 5, 1)-R3( 2, 1)                                  
C                                                                       
      C3(  2,3,4)=+R4( 8, 1)-R3( 4, 1)                                  
      C3(  3,3,4)=+R4( 9, 1)-R3( 5, 1)                                  
C                                                                       
      C3(  4,3,4)=+R4(12, 1)-R3( 7, 1)                                  
      C3(  5,3,4)=+R4(13, 1)-R3( 8, 1)                                  
      C3(  6,3,4)=+R4(14, 1)-R3( 9, 1)                                  
C                                                                       
      C2(  1,3,4)=+R3( 5, 5)-R2( 4, 9)                                  
C                                                                       
      C2(  2,3,4)=+R3( 8, 5)-R2( 2, 9)                                  
      C2(  3,3,4)=+R3( 9, 5)-R2( 6, 9)                                  
C                                                                       
      C1(  1,3,4)=+S2(10, 6)-S1(12, 2)                                  
      C1(  2,3,4)=+S2(11, 6)-S1(13, 2)                                  
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..1,  K = 4  and  L = 4                   
Cjms                                                                    
      C3(  1,4,4)=+R4( 6, 1)-R3( 3, 1)-R3( 3, 4)+R2( 1, 4)+R2( 1, 5)    
C                                                                       
      C3(  2,4,4)=+R4( 9, 1)-R3( 5, 1)-R3( 5, 4)+R2( 4, 4)+R2( 4, 5)    
      C3(  3,4,4)=+R4(10, 1)-R3( 6, 1)-R3( 6, 4)+R2( 5, 4)+R2( 5, 5)    
C                                                                       
      C3(  4,4,4)=+R4(13, 1)-R3( 8, 1)-R3( 8, 4)+R2( 2, 4)+R2( 2, 5)    
      C3(  5,4,4)=+R4(14, 1)-R3( 9, 1)-R3( 9, 4)+R2( 6, 4)+R2( 6, 5)    
      C3(  6,4,4)=+R4(15, 1)-R3(10, 1)-R3(10, 4)+R2( 3, 4)+R2( 3, 5)    
C                                                                       
      C2(  1,4,4)=+R3( 6, 5)-R2( 5, 8)-R2( 5, 9)+R1( 1, 4)+R1( 1, 5)    
C                                                                       
      C2(  2,4,4)=+R3( 9, 5)-R2( 6, 8)-R2( 6, 9)+R1( 2, 4)+R1( 2, 5)    
      C2(  3,4,4)=+R3(10, 5)-R2( 3, 8)-R2( 3, 9)+R1( 3, 4)+R1( 3, 5)    
C                                                                       
      C1(  1,4,4)=+S2(10, 3)-S1(10, 3)-S1(12, 3)+R0(    4)+R0(    5)    
      C1(  2,4,4)=+S2(11, 3)-S1(11, 3)-S1(13, 3)+R0(    9)+R0(   10)    
C                                                                       
      LSYM13=(KX.EQ.4 .AND. LX.EQ.4)                                    
      DO L=1,LX                                                         
         DO K=1,KX                                                      
C                                                                       
            IF(LSYM13) THEN                                             
               IF(K.EQ.2 .AND. L.EQ.3) GO TO 001                        
            ENDIF                                                       
C                                                                       
            F(1,1,K,L)=+C3(  1,K,L)+C1(  1,K,L)                         
     *               +(+C2(  1,K,L)+C2(  1,K,L)+C1(  2,K,L)*QX)*QX      
            F(2,1,K,L)=+C3(  4,K,L)+C1(  1,K,L)                         
            F(3,1,K,L)=+C3(  6,K,L)+C1(  1,K,L)                         
     *               +(+C2(  3,K,L)+C2(  3,K,L)+C1(  2,K,L)*QZ)*QZ      
            F(4,1,K,L)=+C3(  2,K,L)+C2(  2,K,L)*QX                      
            F(5,1,K,L)=+C3(  3,K,L)+C2(  3,K,L)*QX                      
     *               +(+C2(  1,K,L)+C1(  2,K,L)*QX)*QZ                  
            F(6,1,K,L)=+C3(  5,K,L)+C2(  2,K,L)*QZ                      
C                                                                       
  001       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
C                                                                       
      IF(.NOT.LSYM13) GO TO 130                                         
CC    DO J=1,1                                                          
         DO I=1,6                                                       
            F(I,1,2,3)= F(I,1,3,2)                                      
         ENDDO                                                          
CC    ENDDO                                                             
  130 CONTINUE                                                          
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2S   *DECK MCDV14                                           
C>                                                                      
C>    @brief   DDPS case                                                
C>                                                                      
C>    @details integration of a DDPS case                               
C>                                                                      
      SUBROUTINE MCDV14(F,QX,QZ)                                        
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      DIMENSION F(6,6,6,6)                                              
C                                                                       
      COMMON /KI4   / R0(   25),R1( 3,40),R2( 6,56),R3(10,52),          
     *                R4(15,42),R5(21,24),R6(28,12),R7(36, 4),R8(   45) 
C$omp threadprivate(/KI4/)
C                                                                       
      DIMENSION                 S1(13, 3),S2(12, 6),S3( 8,10),          
     *                S4( 4,15)                                         
Cjms                                                                    
C     Simplified calculation of F(I,J,K,L) for cases                    
C     where:  I = 1..6,  J = 1..6,  K = 1..KX  and  L = 1..LX           
C     using auxiliary arrays E1, E2, E3, E4 and E5                      
Cjms                                                                    
      PARAMETER (KX=4)                                                  
      PARAMETER (LX=1)                                                  
      DIMENSION       E1(  5,KX,LX),E2(4,3,KX,LX),E3(4,6,KX,LX)         
      DIMENSION       E4(2,10,KX,LX),E5( 15,KX,LX)                      
C                                                                       
      DIMENSION INS(KX,LX),INI(KX,LX),IND(15,KX,LX)                     
      DIMENSION IN6(6)                                                  
      DATA INS/ 0, 0, 0, 1/                                             
      DATA INI/ 0, 0, 1, 1/                                             
      DATA IN6/ 1, 4, 5, 2, 6, 3/                                       
C                                                                       
      DO L=1,LX                                                         
         DO K=1,KX                                                      
            IJ= 0                                                       
            IN= INS(K,L)                                                
            DO I=1,5                                                    
               IN= IN+INI(K,L)                                          
               DO J=1,I                                                 
                  IJ= IJ+1                                              
                  IN= IN+1                                              
                  IND(IJ,K,L)= IN                                       
               ENDDO                                                    
            ENDDO                                                       
         ENDDO                                                          
      ENDDO                                                             
      DO 101 J=1,13                                                     
         DO 101 I=1, 3                                                  
  101 S1(J,I)= R1(I,J)                                                  
      DO 102 J=1,12                                                     
         DO 102 I=1, 6                                                  
  102 S2(J,I)= R2(I,J)                                                  
      DO 103 J=1, 8                                                     
         DO 103 I=1,10                                                  
  103 S3(J,I)= R3(I,J)                                                  
      DO 104 J=1, 4                                                     
         DO 104 I=1,15                                                  
  104 S4(J,I)= R4(I,J)                                                  
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 1  and  L = 1                   
Cjms                                                                    
      DO J=1,15                                                         
         E5(  J,1,1)=+R4(J,1)                                           
      ENDDO                                                             
C                                                                       
      DO J=1,10                                                         
         E4(1,J,1,1)=+S3(1,J)                                           
         E4(2,J,1,1)=+S3(2,J)                                           
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         M= IN6(J)                                                      
         E3(1,J,1,1)=+S2(1,M)                                           
         E3(2,J,1,1)=+S2(2,M)                                           
         E3(3,J,1,1)=+S2(3,M)                                           
         E3(4,J,1,1)=+S2(4,M)                                           
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         E2(1,J,1,1)=+S1(1,J)                                           
         E2(2,J,1,1)=+S1(2,J)                                           
         E2(3,J,1,1)=+S1(3,J)                                           
         E2(4,J,1,1)=+S1(4,J)                                           
      ENDDO                                                             
C                                                                       
      DO I=1,5                                                          
         E1(I  ,1,1)=+R0(I)                                             
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 2  and  L = 1                   
Cjms                                                                    
      DO J=1,15                                                         
         E5(  J,2,1)=-R5( J,1)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,10                                                         
         E4(1,J,2,1)=-S4( 3,J)                                          
         E4(2,J,2,1)=-S4( 4,J)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         E3(1,J,2,1)=-S3( 5,J)                                          
         E3(2,J,2,1)=-S3( 6,J)                                          
         E3(3,J,2,1)=-S3( 7,J)                                          
         E3(4,J,2,1)=-S3( 8,J)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         M= IN6(J)                                                      
         E2(1,J,2,1)=-S2( 9,M)                                          
         E2(2,J,2,1)=-S2(10,M)                                          
         E2(3,J,2,1)=-S2(11,M)                                          
         E2(4,J,2,1)=-S2(12,M)                                          
      ENDDO                                                             
C                                                                       
         J=1                                                            
      DO I=1,5                                                          
         E1(I  ,2,1)=-S1(I+ 8,J)                                        
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 3  and  L = 1                   
Cjms                                                                    
      DO J=1,15                                                         
         K= IND(J,3,1)                                                  
         E5(  J,3,1)=-R5( K,1)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,10                                                         
         K= IND(J,3,1)                                                  
         E4(1,J,3,1)=-S4( 3,K)                                          
         E4(2,J,3,1)=-S4( 4,K)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         K= IND(J,3,1)                                                  
         E3(1,J,3,1)=-S3( 5,K)                                          
         E3(2,J,3,1)=-S3( 6,K)                                          
         E3(3,J,3,1)=-S3( 7,K)                                          
         E3(4,J,3,1)=-S3( 8,K)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         K= IND(J,3,1)                                                  
         M= IN6(K)                                                      
         E2(1,J,3,1)=-S2( 9,M)                                          
         E2(2,J,3,1)=-S2(10,M)                                          
         E2(3,J,3,1)=-S2(11,M)                                          
         E2(4,J,3,1)=-S2(12,M)                                          
      ENDDO                                                             
C                                                                       
         J=1                                                            
         K= IND(J,3,1)                                                  
      DO I=1,5                                                          
         E1(I  ,3,1)=-S1(I+ 8,K)                                        
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 4  and  L = 1                   
Cjms                                                                    
      DO J=1,15                                                         
         K= IND(J,4,1)                                                  
         E5(  J,4,1)=-R5( K,1)+R4(J,2)                                  
      ENDDO                                                             
C                                                                       
      DO J=1,10                                                         
         K= IND(J,4,1)                                                  
         E4(1,J,4,1)=-S4( 3,K)+S3(3,J)                                  
         E4(2,J,4,1)=-S4( 4,K)+S3(4,J)                                  
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         K= IND(J,4,1)                                                  
         M= IN6(J)                                                      
         E3(1,J,4,1)=-S3( 5,K)+S2( 5,M)                                 
         E3(2,J,4,1)=-S3( 6,K)+S2( 6,M)                                 
         E3(3,J,4,1)=-S3( 7,K)+S2( 7,M)                                 
         E3(4,J,4,1)=-S3( 8,K)+S2( 8,M)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         K= IND(J,4,1)                                                  
         M= IN6(K)                                                      
         E2(1,J,4,1)=-S2( 9,M)+S1( 5,J)                                 
         E2(2,J,4,1)=-S2(10,M)+S1( 6,J)                                 
         E2(3,J,4,1)=-S2(11,M)+S1( 7,J)                                 
         E2(4,J,4,1)=-S2(12,M)+S1( 8,J)                                 
      ENDDO                                                             
C                                                                       
         J=1                                                            
         K= IND(J,4,1)                                                  
      DO I=1,5                                                          
         E1(I  ,4,1)=-S1(I+ 8,K)+R0(I+ 5)                               
      ENDDO                                                             
C                                                                       
      CALL FIJKL5(KX,LX,F,QX,QZ,E1,E2,E3,E4,E5)                         
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2S   *DECK MCDV15                                           
C>                                                                      
C>    @brief   DPDS case                                                
C>                                                                      
C>    @details integration of a DPDS case                               
C>                                                                      
      SUBROUTINE MCDV15(F,QX,QZ)                                        
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      DIMENSION F(6,6,6,6)                                              
C                                                                       
      COMMON /KI4   / R0(   25),R1( 3,40),R2( 6,56),R3(10,52),          
     *                R4(15,42),R5(21,24),R6(28,12),R7(36, 4),R8(   45) 
C$omp threadprivate(/KI4/)
C                                                                       
      DIMENSION                 S1(13, 3),S2(15, 6),S3( 9,10),          
     *                S4( 4,15)                                         
C                                                                       
      PARAMETER (F02=2.0D+00)                                           
Cjms                                                                    
C     Simplified calculation of F(I,J,K,L) for cases                    
C     where:  I = 1..6,  J = 1..4,  K = 1..KX  and  L = 1..LX           
C     using auxiliary arrays D1, D2, D3 and D4                          
Cjms                                                                    
      PARAMETER (KX=6)                                                  
      PARAMETER (LX=1)                                                  
      DIMENSION       D1(  5,KX,LX),D2(4,3,KX,LX),D3(3,6,KX,LX)         
      DIMENSION       D4( 10,KX,LX)                                     
C                                                                       
      DIMENSION INS(KX,LX),INI(KX,LX),IND(10,KX,LX)                     
      DIMENSION IN6(6)                                                  
      DATA INS/ 0, 1, 3, 0, 1, 2/                                       
      DATA INI/ 0, 2, 2, 1, 1, 2/                                       
      DATA IN6/ 1, 4, 5, 2, 6, 3/                                       
C                                                                       
      DO L=1,LX                                                         
         DO K=1,KX                                                      
            IJ= 0                                                       
            IN= INS(K,L)                                                
            DO I=1,4                                                    
               IN= IN+INI(K,L)                                          
               DO J=1,I                                                 
                  IJ= IJ+1                                              
                  IN= IN+1                                              
                  IND(IJ,K,L)= IN                                       
               ENDDO                                                    
            ENDDO                                                       
         ENDDO                                                          
      ENDDO                                                             
      DO 101 J=1,13                                                     
         DO 101 I=1, 3                                                  
  101 S1(J,I)= R1(I,J)                                                  
      DO 102 J=1,15                                                     
         DO 102 I=1, 6                                                  
  102 S2(J,I)= R2(I,J)                                                  
      DO 103 J=1, 9                                                     
         DO 103 I=1,10                                                  
  103 S3(J,I)= R3(I,J)                                                  
      DO 104 J=1, 4                                                     
         DO 104 I=1,15                                                  
  104 S4(J,I)= R4(I,J)                                                  
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..4,  K = 1  and  L = 1                   
Cjms                                                                    
      DO J=1,10                                                         
         D4(  J,1,1)=+R5( J,1)+R3( J,1)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         M= IN6(J)                                                      
         D3(1,J,1,1)=+S4( 2,J)+S2( 1,M)                                 
         D3(2,J,1,1)=+S4( 3,J)+S2( 2,M)                                 
         D3(3,J,1,1)=+S4( 4,J)+S2( 3,M)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         D2(1,J,1,1)=+S3( 6,J)+S1( 1,J)                                 
         D2(2,J,1,1)=+S3( 7,J)+S1( 2,J)                                 
         D2(3,J,1,1)=+S3( 8,J)+S1( 3,J)                                 
         D2(4,J,1,1)=+S3( 9,J)+S1( 4,J)                                 
      ENDDO                                                             
C                                                                       
         J=1                                                            
      DO I=1,5                                                          
         D1(I  ,1,1)=+S2(I+10,J)+R0(I)                                  
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..4,  K = 2  and  L = 1                   
Cjms                                                                    
      DO J=1,10                                                         
         K= IND(J,2,1)                                                  
         D4(  J,2,1)=+R5( K,1)+R3( J,1)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         K= IND(J,2,1)                                                  
         M= IN6(J)                                                      
         D3(1,J,2,1)=+S4( 2,K)+S2( 1,M)                                 
         D3(2,J,2,1)=+S4( 3,K)+S2( 2,M)                                 
         D3(3,J,2,1)=+S4( 4,K)+S2( 3,M)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         K= IND(J,2,1)                                                  
         D2(1,J,2,1)=+S3( 6,K)+S1( 1,J)                                 
         D2(2,J,2,1)=+S3( 7,K)+S1( 2,J)                                 
         D2(3,J,2,1)=+S3( 8,K)+S1( 3,J)                                 
         D2(4,J,2,1)=+S3( 9,K)+S1( 4,J)                                 
      ENDDO                                                             
C                                                                       
         J=1                                                            
         K= IND(J,2,1)                                                  
         M= IN6(K)                                                      
      DO I=1,5                                                          
         D1(I  ,2,1)=+S2(I+10,M)+R0(I)                                  
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..4,  K = 3  and  L = 1                   
Cjms                                                                    
      D4(  1,3,1)=+R5( 6, 1)-R4( 3, 1)*F02+R3( 1, 1)+R3( 1, 2)          
C                                                                       
      D4(  2,3,1)=+R5( 9, 1)-R4( 5, 1)*F02+R3( 2, 1)+R3( 2, 2)          
      D4(  3,3,1)=+R5(10, 1)-R4( 6, 1)*F02+R3( 3, 1)+R3( 3, 2)          
C                                                                       
      D4(  4,3,1)=+R5(13, 1)-R4( 8, 1)*F02+R3( 4, 1)+R3( 4, 2)          
      D4(  5,3,1)=+R5(14, 1)-R4( 9, 1)*F02+R3( 5, 1)+R3( 5, 2)          
      D4(  6,3,1)=+R5(15, 1)-R4(10, 1)*F02+R3( 6, 1)+R3( 6, 2)          
C                                                                       
      D4(  7,3,1)=+R5(18, 1)-R4(12, 1)*F02+R3( 7, 1)+R3( 7, 2)          
      D4(  8,3,1)=+R5(19, 1)-R4(13, 1)*F02+R3( 8, 1)+R3( 8, 2)          
      D4(  9,3,1)=+R5(20, 1)-R4(14, 1)*F02+R3( 9, 1)+R3( 9, 2)          
      D4( 10,3,1)=+R5(21, 1)-R4(15, 1)*F02+R3(10, 1)+R3(10, 2)          
C                                                                       
      D3(1,1,3,1)=+S4( 2, 6)-S3( 3, 3)*F02+S2( 1, 1)+S2( 4, 1)          
      D3(2,1,3,1)=+S4( 3, 6)-S3( 4, 3)*F02+S2( 2, 1)+S2( 5, 1)          
      D3(3,1,3,1)=+S4( 4, 6)-S3( 5, 3)*F02+S2( 3, 1)+S2( 6, 1)          
C                                                                       
      D3(1,2,3,1)=+S4( 2, 9)-S3( 3, 5)*F02+S2( 1, 4)+S2( 4, 4)          
      D3(2,2,3,1)=+S4( 3, 9)-S3( 4, 5)*F02+S2( 2, 4)+S2( 5, 4)          
      D3(3,2,3,1)=+S4( 4, 9)-S3( 5, 5)*F02+S2( 3, 4)+S2( 6, 4)          
      D3(1,3,3,1)=+S4( 2,10)-S3( 3, 6)*F02+S2( 1, 5)+S2( 4, 5)          
      D3(2,3,3,1)=+S4( 3,10)-S3( 4, 6)*F02+S2( 2, 5)+S2( 5, 5)          
      D3(3,3,3,1)=+S4( 4,10)-S3( 5, 6)*F02+S2( 3, 5)+S2( 6, 5)          
C                                                                       
      D3(1,4,3,1)=+S4( 2,13)-S3( 3, 8)*F02+S2( 1, 2)+S2( 4, 2)          
      D3(2,4,3,1)=+S4( 3,13)-S3( 4, 8)*F02+S2( 2, 2)+S2( 5, 2)          
      D3(3,4,3,1)=+S4( 4,13)-S3( 5, 8)*F02+S2( 3, 2)+S2( 6, 2)          
      D3(1,5,3,1)=+S4( 2,14)-S3( 3, 9)*F02+S2( 1, 6)+S2( 4, 6)          
      D3(2,5,3,1)=+S4( 3,14)-S3( 4, 9)*F02+S2( 2, 6)+S2( 5, 6)          
      D3(3,5,3,1)=+S4( 4,14)-S3( 5, 9)*F02+S2( 3, 6)+S2( 6, 6)          
      D3(1,6,3,1)=+S4( 2,15)-S3( 3,10)*F02+S2( 1, 3)+S2( 4, 3)          
      D3(2,6,3,1)=+S4( 3,15)-S3( 4,10)*F02+S2( 2, 3)+S2( 5, 3)          
      D3(3,6,3,1)=+S4( 4,15)-S3( 5,10)*F02+S2( 3, 3)+S2( 6, 3)          
C                                                                       
      D2(1,1,3,1)=+S3( 6, 6)-S2( 7, 5)*F02+S1( 1, 1)+S1( 5, 1)          
      D2(2,1,3,1)=+S3( 7, 6)-S2( 8, 5)*F02+S1( 2, 1)+S1( 6, 1)          
      D2(3,1,3,1)=+S3( 8, 6)-S2( 9, 5)*F02+S1( 3, 1)+S1( 7, 1)          
      D2(4,1,3,1)=+S3( 9, 6)-S2(10, 5)*F02+S1( 4, 1)+S1( 8, 1)          
C                                                                       
      D2(1,2,3,1)=+S3( 6, 9)-S2( 7, 6)*F02+S1( 1, 2)+S1( 5, 2)          
      D2(2,2,3,1)=+S3( 7, 9)-S2( 8, 6)*F02+S1( 2, 2)+S1( 6, 2)          
      D2(3,2,3,1)=+S3( 8, 9)-S2( 9, 6)*F02+S1( 3, 2)+S1( 7, 2)          
      D2(4,2,3,1)=+S3( 9, 9)-S2(10, 6)*F02+S1( 4, 2)+S1( 8, 2)          
      D2(1,3,3,1)=+S3( 6,10)-S2( 7, 3)*F02+S1( 1, 3)+S1( 5, 3)          
      D2(2,3,3,1)=+S3( 7,10)-S2( 8, 3)*F02+S1( 2, 3)+S1( 6, 3)          
      D2(3,3,3,1)=+S3( 8,10)-S2( 9, 3)*F02+S1( 3, 3)+S1( 7, 3)          
      D2(4,3,3,1)=+S3( 9,10)-S2(10, 3)*F02+S1( 4, 3)+S1( 8, 3)          
C                                                                       
      DO I=1,5                                                          
         D1(I  ,3,1)=+S2(I+10,3)-S1(I+ 8,3)*F02+R0(I   )+R0(I+ 5)       
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..4,  K = 4  and  L = 1                   
Cjms                                                                    
      DO J=1,10                                                         
         K= IND(J,4,1)                                                  
         D4(  J,4,1)=+R5( K,1)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         K= IND(J,4,1)                                                  
         D3(1,J,4,1)=+S4( 2,K)                                          
         D3(2,J,4,1)=+S4( 3,K)                                          
         D3(3,J,4,1)=+S4( 4,K)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         K= IND(J,4,1)                                                  
         D2(1,J,4,1)=+S3( 6,K)                                          
         D2(2,J,4,1)=+S3( 7,K)                                          
         D2(3,J,4,1)=+S3( 8,K)                                          
         D2(4,J,4,1)=+S3( 9,K)                                          
      ENDDO                                                             
C                                                                       
         J=1                                                            
         K= IND(J,4,1)                                                  
         M= IN6(K)                                                      
      DO I=1,5                                                          
         D1(I  ,4,1)=+S2(I+10,M)                                        
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..4,  K = 5  and  L = 1                   
Cjms                                                                    
      DO J=1,10                                                         
         K= IND(J,5,1)                                                  
         D4(  J,5,1)=+R5( K,1)-R4( J,1)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         K= IND(J,5,1)                                                  
         D3(1,J,5,1)=+S4( 2,K)-S3( 3,J)                                 
         D3(2,J,5,1)=+S4( 3,K)-S3( 4,J)                                 
         D3(3,J,5,1)=+S4( 4,K)-S3( 5,J)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         K= IND(J,5,1)                                                  
         M= IN6(J)                                                      
         D2(1,J,5,1)=+S3( 6,K)-S2( 7,M)                                 
         D2(2,J,5,1)=+S3( 7,K)-S2( 8,M)                                 
         D2(3,J,5,1)=+S3( 8,K)-S2( 9,M)                                 
         D2(4,J,5,1)=+S3( 9,K)-S2(10,M)                                 
      ENDDO                                                             
C                                                                       
         J=1                                                            
         K= IND(J,5,1)                                                  
         M= IN6(K)                                                      
      DO I=1,5                                                          
         D1(I  ,5,1)=+S2(I+10,M)-S1(I+ 8,J)                             
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..4,  K = 6  and  L = 1                   
Cjms                                                                    
      D4(  1,6,1)=+R5( 5, 1)-R4( 2, 1)                                  
C                                                                       
      D4(  2,6,1)=+R5( 8, 1)-R4( 4, 1)                                  
      D4(  3,6,1)=+R5( 9, 1)-R4( 5, 1)                                  
C                                                                       
      D4(  4,6,1)=+R5(12, 1)-R4( 7, 1)                                  
      D4(  5,6,1)=+R5(13, 1)-R4( 8, 1)                                  
      D4(  6,6,1)=+R5(14, 1)-R4( 9, 1)                                  
C                                                                       
      D4(  7,6,1)=+R5(17, 1)-R4(11, 1)                                  
      D4(  8,6,1)=+R5(18, 1)-R4(12, 1)                                  
      D4(  9,6,1)=+R5(19, 1)-R4(13, 1)                                  
      D4( 10,6,1)=+R5(20, 1)-R4(14, 1)                                  
C                                                                       
      D3(1,1,6,1)=+S4( 2, 5)-S3( 3, 2)                                  
      D3(2,1,6,1)=+S4( 3, 5)-S3( 4, 2)                                  
      D3(3,1,6,1)=+S4( 4, 5)-S3( 5, 2)                                  
C                                                                       
      D3(1,2,6,1)=+S4( 2, 8)-S3( 3, 4)                                  
      D3(2,2,6,1)=+S4( 3, 8)-S3( 4, 4)                                  
      D3(3,2,6,1)=+S4( 4, 8)-S3( 5, 4)                                  
      D3(1,3,6,1)=+S4( 2, 9)-S3( 3, 5)                                  
      D3(2,3,6,1)=+S4( 3, 9)-S3( 4, 5)                                  
      D3(3,3,6,1)=+S4( 4, 9)-S3( 5, 5)                                  
C                                                                       
      D3(1,4,6,1)=+S4( 2,12)-S3( 3, 7)                                  
      D3(2,4,6,1)=+S4( 3,12)-S3( 4, 7)                                  
      D3(3,4,6,1)=+S4( 4,12)-S3( 5, 7)                                  
      D3(1,5,6,1)=+S4( 2,13)-S3( 3, 8)                                  
      D3(2,5,6,1)=+S4( 3,13)-S3( 4, 8)                                  
      D3(3,5,6,1)=+S4( 4,13)-S3( 5, 8)                                  
      D3(1,6,6,1)=+S4( 2,14)-S3( 3, 9)                                  
      D3(2,6,6,1)=+S4( 3,14)-S3( 4, 9)                                  
      D3(3,6,6,1)=+S4( 4,14)-S3( 5, 9)                                  
C                                                                       
      D2(1,1,6,1)=+S3( 6, 5)-S2( 7, 4)                                  
      D2(2,1,6,1)=+S3( 7, 5)-S2( 8, 4)                                  
      D2(3,1,6,1)=+S3( 8, 5)-S2( 9, 4)                                  
      D2(4,1,6,1)=+S3( 9, 5)-S2(10, 4)                                  
C                                                                       
      D2(1,2,6,1)=+S3( 6, 8)-S2( 7, 2)                                  
      D2(2,2,6,1)=+S3( 7, 8)-S2( 8, 2)                                  
      D2(3,2,6,1)=+S3( 8, 8)-S2( 9, 2)                                  
      D2(4,2,6,1)=+S3( 9, 8)-S2(10, 2)                                  
      D2(1,3,6,1)=+S3( 6, 9)-S2( 7, 6)                                  
      D2(2,3,6,1)=+S3( 7, 9)-S2( 8, 6)                                  
      D2(3,3,6,1)=+S3( 8, 9)-S2( 9, 6)                                  
      D2(4,3,6,1)=+S3( 9, 9)-S2(10, 6)                                  
C                                                                       
      DO I=1,5                                                          
         D1(I  ,6,1)=+S2(I+10,6)-S1(I+ 8,2)                             
      ENDDO                                                             
C                                                                       
      CALL FIJKL4(KX,LX,F,QX,QZ,D1,D2,D3,D4)                            
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2S   *DECK MCDV16                                           
C>                                                                      
C>    @brief   DPPP case                                                
C>                                                                      
C>    @details integration of a DPPP case                               
C>                                                                      
      SUBROUTINE MCDV16(F,QX,QZ)                                        
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      DIMENSION F(6,6,6,6)                                              
C                                                                       
      COMMON /KI4   / R0(   25),R1( 3,40),R2( 6,56),R3(10,52),          
     *                R4(15,42),R5(21,24),R6(28,12),R7(36, 4),R8(   45) 
C$omp threadprivate(/KI4/)
C                                                                       
      DIMENSION                 S1(40, 3),S2(36, 6),S3(21,10),          
     *                S4(14,15)                                         
Cjms                                                                    
C     Simplified calculation of F(I,J,K,L) for cases                    
C     where:  I = 1..6,  J = 1..4,  K = 1..KX  and  L = 1..LX           
C     using auxiliary arrays D1, D2, D3 and D4                          
Cjms                                                                    
      PARAMETER (KX=4)                                                  
      PARAMETER (LX=4)                                                  
      DIMENSION       D1(  5,KX,LX),D2(4,3,KX,LX),D3(3,6,KX,LX)         
      DIMENSION       D4( 10,KX,LX)                                     
C                                                                       
      DIMENSION INS(KX,LX),INI(KX,LX),IND(10,KX,LX)                     
      DIMENSION IN6(6)                                                  
      DATA INS/ 0, 0, 0, 1,  0, 0, 0, 1,  0, 0, 1, 2,  1, 1, 2, 3/      
      DATA INI/ 0, 0, 1, 1,  0, 0, 1, 1,  1, 1, 2, 2,  1, 1, 2, 2/      
      DATA IN6/ 1, 4, 5, 2, 6, 3/                                       
C                                                                       
      DO L=1,LX                                                         
         DO K=1,KX                                                      
            IJ= 0                                                       
            IN= INS(K,L)                                                
            DO I=1,4                                                    
               IN= IN+INI(K,L)                                          
               DO J=1,I                                                 
                  IJ= IJ+1                                              
                  IN= IN+1                                              
                  IND(IJ,K,L)= IN                                       
               ENDDO                                                    
            ENDDO                                                       
         ENDDO                                                          
      ENDDO                                                             
      DO 101 J=1,40                                                     
         DO 101 I=1, 3                                                  
  101 S1(J,I)= R1(I,J)                                                  
      DO 102 J=1,36                                                     
         DO 102 I=1, 6                                                  
  102 S2(J,I)= R2(I,J)                                                  
      DO 103 J=1,21                                                     
         DO 103 I=1,10                                                  
  103 S3(J,I)= R3(I,J)                                                  
      DO 104 J=1,14                                                     
         DO 104 I=1,15                                                  
  104 S4(J,I)= R4(I,J)                                                  
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..4,  K = 1  and  L = 1                   
Cjms                                                                    
      DO J=1,10                                                         
         D4(  J,1,1)=+R3( J,1)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         M= IN6(J)                                                      
         D3(1,J,1,1)=+S2( 1,M)                                          
         D3(2,J,1,1)=+S2( 2,M)                                          
         D3(3,J,1,1)=+S2( 3,M)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         D2(1,J,1,1)=+S1( 1,J)                                          
         D2(2,J,1,1)=+S1( 2,J)                                          
         D2(3,J,1,1)=+S1( 3,J)                                          
         D2(4,J,1,1)=+S1( 4,J)                                          
      ENDDO                                                             
C                                                                       
      DO I=1,5                                                          
         D1(I  ,1,1)=+R0(I)                                             
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..4,  K = 2  and  L = 1                   
Cjms                                                                    
      DO J=1,10                                                         
         D4(  J,2,1)=-R4( J,2)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         D3(1,J,2,1)=-S3( 9,J)                                          
         D3(2,J,2,1)=-S3(10,J)                                          
         D3(3,J,2,1)=-S3(11,J)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         M= IN6(J)                                                      
         D2(1,J,2,1)=-S2(16,M)                                          
         D2(2,J,2,1)=-S2(17,M)                                          
         D2(3,J,2,1)=-S2(18,M)                                          
         D2(4,J,2,1)=-S2(19,M)                                          
      ENDDO                                                             
C                                                                       
         J=1                                                            
      DO I=1,5                                                          
         D1(I  ,2,1)=-S1(I+20,J)                                        
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..4,  K = 3  and  L = 1                   
Cjms                                                                    
      DO J=1,10                                                         
         K= IND(J,3,1)                                                  
         D4(  J,3,1)=-R4( K,2)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         K= IND(J,3,1)                                                  
         D3(1,J,3,1)=-S3( 9,K)                                          
         D3(2,J,3,1)=-S3(10,K)                                          
         D3(3,J,3,1)=-S3(11,K)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         K= IND(J,3,1)                                                  
         M= IN6(K)                                                      
         D2(1,J,3,1)=-S2(16,M)                                          
         D2(2,J,3,1)=-S2(17,M)                                          
         D2(3,J,3,1)=-S2(18,M)                                          
         D2(4,J,3,1)=-S2(19,M)                                          
      ENDDO                                                             
C                                                                       
         J=1                                                            
         K= IND(J,3,1)                                                  
      DO I=1,5                                                          
         D1(I  ,3,1)=-S1(I+20,K)                                        
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..4,  K = 4  and  L = 1                   
Cjms                                                                    
      DO J=1,10                                                         
         K= IND(J,4,1)                                                  
         D4(  J,4,1)=-R4( K,2)+R3( J,2)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         K= IND(J,4,1)                                                  
         M= IN6(J)                                                      
         D3(1,J,4,1)=-S3( 9,K)+S2( 4,M)                                 
         D3(2,J,4,1)=-S3(10,K)+S2( 5,M)                                 
         D3(3,J,4,1)=-S3(11,K)+S2( 6,M)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         K= IND(J,4,1)                                                  
         M= IN6(K)                                                      
         D2(1,J,4,1)=-S2(16,M)+S1( 5,J)                                 
         D2(2,J,4,1)=-S2(17,M)+S1( 6,J)                                 
         D2(3,J,4,1)=-S2(18,M)+S1( 7,J)                                 
         D2(4,J,4,1)=-S2(19,M)+S1( 8,J)                                 
      ENDDO                                                             
C                                                                       
         J=1                                                            
         K= IND(J,4,1)                                                  
      DO I=1,5                                                          
         D1(I  ,4,1)=-S1(I+20,K)+R0(I+ 5)                               
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..4,  K = 1  and  L = 2                   
Cjms                                                                    
      DO J=1,10                                                         
         D4(  J,1,2)=-R4( J,3)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         D3(1,J,1,2)=-S3(12,J)                                          
         D3(2,J,1,2)=-S3(13,J)                                          
         D3(3,J,1,2)=-S3(14,J)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         M= IN6(J)                                                      
         D2(1,J,1,2)=-S2(20,M)                                          
         D2(2,J,1,2)=-S2(21,M)                                          
         D2(3,J,1,2)=-S2(22,M)                                          
         D2(4,J,1,2)=-S2(23,M)                                          
      ENDDO                                                             
C                                                                       
         J=1                                                            
      DO I=1,5                                                          
         D1(I  ,1,2)=-S1(I+25,J)                                        
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..4,  K = 2  and  L = 2                   
Cjms                                                                    
      DO J=1,10                                                         
         D4(  J,2,2)=+R5( J,1)+R3( J,5)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         M= IN6(J)                                                      
         D3(1,J,2,2)=+S4( 5,J)+S2(13,M)                                 
         D3(2,J,2,2)=+S4( 6,J)+S2(14,M)                                 
         D3(3,J,2,2)=+S4( 7,J)+S2(15,M)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         D2(1,J,2,2)=+S3(18,J)+S1(17,J)                                 
         D2(2,J,2,2)=+S3(19,J)+S1(18,J)                                 
         D2(3,J,2,2)=+S3(20,J)+S1(19,J)                                 
         D2(4,J,2,2)=+S3(21,J)+S1(20,J)                                 
      ENDDO                                                             
C                                                                       
         J=1                                                            
      DO I=1,5                                                          
         D1(I  ,2,2)=+S2(I+31,J)+R0(I+20)                               
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..4,  K = 3  and  L = 2                   
Cjms                                                                    
      DO J=1,10                                                         
         K= IND(J,3,2)                                                  
         D4(  J,3,2)=+R5( K,1)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         K= IND(J,3,2)                                                  
         D3(1,J,3,2)=+S4( 5,K)                                          
         D3(2,J,3,2)=+S4( 6,K)                                          
         D3(3,J,3,2)=+S4( 7,K)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         K= IND(J,3,2)                                                  
         D2(1,J,3,2)=+S3(18,K)                                          
         D2(2,J,3,2)=+S3(19,K)                                          
         D2(3,J,3,2)=+S3(20,K)                                          
         D2(4,J,3,2)=+S3(21,K)                                          
      ENDDO                                                             
C                                                                       
         J=1                                                            
         K= IND(J,3,2)                                                  
         M= IN6(K)                                                      
      DO I=1,5                                                          
         D1(I  ,3,2)=+S2(I+31,M)                                        
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..4,  K = 4  and  L = 2                   
Cjms                                                                    
      DO J=1,10                                                         
         K= IND(J,4,2)                                                  
         D4(  J,4,2)=+R5( K,1)-R4( J,4)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         K= IND(J,4,2)                                                  
         D3(1,J,4,2)=+S4( 5,K)-S3(15,J)                                 
         D3(2,J,4,2)=+S4( 6,K)-S3(16,J)                                 
         D3(3,J,4,2)=+S4( 7,K)-S3(17,J)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         K= IND(J,4,2)                                                  
         M= IN6(J)                                                      
         D2(1,J,4,2)=+S3(18,K)-S2(24,M)                                 
         D2(2,J,4,2)=+S3(19,K)-S2(25,M)                                 
         D2(3,J,4,2)=+S3(20,K)-S2(26,M)                                 
         D2(4,J,4,2)=+S3(21,K)-S2(27,M)                                 
      ENDDO                                                             
C                                                                       
         J=1                                                            
         K= IND(J,4,2)                                                  
         M= IN6(K)                                                      
      DO I=1,5                                                          
         D1(I  ,4,2)=+S2(I+31,M)-S1(I+30,J)                             
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..4,  K = 1  and  L = 3                   
Cjms                                                                    
      DO J=1,10                                                         
         K= IND(J,1,3)                                                  
         D4(  J,1,3)=-R4( K,3)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         K= IND(J,1,3)                                                  
         D3(1,J,1,3)=-S3(12,K)                                          
         D3(2,J,1,3)=-S3(13,K)                                          
         D3(3,J,1,3)=-S3(14,K)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         K= IND(J,1,3)                                                  
         M= IN6(K)                                                      
         D2(1,J,1,3)=-S2(20,M)                                          
         D2(2,J,1,3)=-S2(21,M)                                          
         D2(3,J,1,3)=-S2(22,M)                                          
         D2(4,J,1,3)=-S2(23,M)                                          
      ENDDO                                                             
C                                                                       
         J=1                                                            
         K= IND(J,1,3)                                                  
      DO I=1,5                                                          
         D1(I  ,1,3)=-S1(I+25,K)                                        
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..4,  K = 2  and  L = 3                   
C                                                                       
C     next statements commented because F(I,J,2,3)                      
C     will be obtained in FIJKL4 from   F(I,J,3,2)                      
Cjms                                                                    
C     DO 324 J=1,10                                                     
C 324 D4(  J,2,3)= D4(  J,3,2)                                          
C     DO 323 J=1,6                                                      
C        DO 323 I=1,3                                                   
C 323 D3(I,J,2,3)= D3(I,J,3,2)                                          
C     DO 322 J=1,3                                                      
C        DO 322 I=1,4                                                   
C 322 D2(I,J,2,3)= D2(I,J,3,2)                                          
C        DO 321 I=1,5                                                   
C 321 D1(I  ,2,3)= D1(I  ,3,2)                                          
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..4,  K = 3  and  L = 3                   
Cjms                                                                    
      DO J=1,10                                                         
         K= IND(J,3,3)                                                  
         D4(  J,3,3)=+R5( K,1)+R3( J,5)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         K= IND(J,3,3)                                                  
         M= IN6(J)                                                      
         D3(1,J,3,3)=+S4( 5,K)+S2(13,M)                                 
         D3(2,J,3,3)=+S4( 6,K)+S2(14,M)                                 
         D3(3,J,3,3)=+S4( 7,K)+S2(15,M)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         K= IND(J,3,3)                                                  
         D2(1,J,3,3)=+S3(18,K)+S1(17,J)                                 
         D2(2,J,3,3)=+S3(19,K)+S1(18,J)                                 
         D2(3,J,3,3)=+S3(20,K)+S1(19,J)                                 
         D2(4,J,3,3)=+S3(21,K)+S1(20,J)                                 
      ENDDO                                                             
C                                                                       
         J=1                                                            
         K= IND(J,3,3)                                                  
         M= IN6(K)                                                      
      DO I=1,5                                                          
         D1(I  ,3,3)=+S2(I+31,M)+R0(I+20)                               
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..4,  K = 4  and  L = 3                   
Cjms                                                                    
      D4(  1,4,3)=+R5( 5, 1)-R4( 2, 4)                                  
C                                                                       
      D4(  2,4,3)=+R5( 8, 1)-R4( 4, 4)                                  
      D4(  3,4,3)=+R5( 9, 1)-R4( 5, 4)                                  
C                                                                       
      D4(  4,4,3)=+R5(12, 1)-R4( 7, 4)                                  
      D4(  5,4,3)=+R5(13, 1)-R4( 8, 4)                                  
      D4(  6,4,3)=+R5(14, 1)-R4( 9, 4)                                  
C                                                                       
      D4(  7,4,3)=+R5(17, 1)-R4(11, 4)                                  
      D4(  8,4,3)=+R5(18, 1)-R4(12, 4)                                  
      D4(  9,4,3)=+R5(19, 1)-R4(13, 4)                                  
      D4( 10,4,3)=+R5(20, 1)-R4(14, 4)                                  
C                                                                       
      D3(1,1,4,3)=+S4( 5, 5)-S3(15, 2)                                  
      D3(2,1,4,3)=+S4( 6, 5)-S3(16, 2)                                  
      D3(3,1,4,3)=+S4( 7, 5)-S3(17, 2)                                  
C                                                                       
      D3(1,2,4,3)=+S4( 5, 8)-S3(15, 4)                                  
      D3(2,2,4,3)=+S4( 6, 8)-S3(16, 4)                                  
      D3(3,2,4,3)=+S4( 7, 8)-S3(17, 4)                                  
      D3(1,3,4,3)=+S4( 5, 9)-S3(15, 5)                                  
      D3(2,3,4,3)=+S4( 6, 9)-S3(16, 5)                                  
      D3(3,3,4,3)=+S4( 7, 9)-S3(17, 5)                                  
C                                                                       
      D3(1,4,4,3)=+S4( 5,12)-S3(15, 7)                                  
      D3(2,4,4,3)=+S4( 6,12)-S3(16, 7)                                  
      D3(3,4,4,3)=+S4( 7,12)-S3(17, 7)                                  
      D3(1,5,4,3)=+S4( 5,13)-S3(15, 8)                                  
      D3(2,5,4,3)=+S4( 6,13)-S3(16, 8)                                  
      D3(3,5,4,3)=+S4( 7,13)-S3(17, 8)                                  
      D3(1,6,4,3)=+S4( 5,14)-S3(15, 9)                                  
      D3(2,6,4,3)=+S4( 6,14)-S3(16, 9)                                  
      D3(3,6,4,3)=+S4( 7,14)-S3(17, 9)                                  
C                                                                       
      D2(1,1,4,3)=+S3(18, 5)-S2(24, 4)                                  
      D2(2,1,4,3)=+S3(19, 5)-S2(25, 4)                                  
      D2(3,1,4,3)=+S3(20, 5)-S2(26, 4)                                  
      D2(4,1,4,3)=+S3(21, 5)-S2(27, 4)                                  
C                                                                       
      D2(1,2,4,3)=+S3(18, 8)-S2(24, 2)                                  
      D2(2,2,4,3)=+S3(19, 8)-S2(25, 2)                                  
      D2(3,2,4,3)=+S3(20, 8)-S2(26, 2)                                  
      D2(4,2,4,3)=+S3(21, 8)-S2(27, 2)                                  
      D2(1,3,4,3)=+S3(18, 9)-S2(24, 6)                                  
      D2(2,3,4,3)=+S3(19, 9)-S2(25, 6)                                  
      D2(3,3,4,3)=+S3(20, 9)-S2(26, 6)                                  
      D2(4,3,4,3)=+S3(21, 9)-S2(27, 6)                                  
C                                                                       
      DO I=1,5                                                          
         D1(I  ,4,3)=+S2(I+31,6)-S1(I+30,2)                             
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..4,  K = 1  and  L = 4                   
Cjms                                                                    
      DO J=1,10                                                         
         K= IND(J,1,4)                                                  
         D4(  J,1,4)=-R4( K,3)+R3( J,3)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         K= IND(J,1,4)                                                  
         M= IN6(J)                                                      
         D3(1,J,1,4)=-S3(12,K)+S2( 7,M)                                 
         D3(2,J,1,4)=-S3(13,K)+S2( 8,M)                                 
         D3(3,J,1,4)=-S3(14,K)+S2( 9,M)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         K= IND(J,1,4)                                                  
         M= IN6(K)                                                      
         D2(1,J,1,4)=-S2(20,M)+S1( 9,J)                                 
         D2(2,J,1,4)=-S2(21,M)+S1(10,J)                                 
         D2(3,J,1,4)=-S2(22,M)+S1(11,J)                                 
         D2(4,J,1,4)=-S2(23,M)+S1(12,J)                                 
      ENDDO                                                             
C                                                                       
         J=1                                                            
         K= IND(J,1,4)                                                  
      DO I=1,5                                                          
         D1(I  ,1,4)=-S1(I+25,K)+R0(I+10)                               
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..4,  K = 2  and  L = 4                   
Cjms                                                                    
      DO J=1,10                                                         
         K= IND(J,2,4)                                                  
         D4(  J,2,4)=+R5( K,1)-R4( J,1)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         K= IND(J,2,4)                                                  
         D3(1,J,2,4)=+S4( 5,K)-S3( 6,J)                                 
         D3(2,J,2,4)=+S4( 6,K)-S3( 7,J)                                 
         D3(3,J,2,4)=+S4( 7,K)-S3( 8,J)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         K= IND(J,2,4)                                                  
         M= IN6(J)                                                      
         D2(1,J,2,4)=+S3(18,K)-S2(28,M)                                 
         D2(2,J,2,4)=+S3(19,K)-S2(29,M)                                 
         D2(3,J,2,4)=+S3(20,K)-S2(30,M)                                 
         D2(4,J,2,4)=+S3(21,K)-S2(31,M)                                 
      ENDDO                                                             
C                                                                       
         J=1                                                            
         K= IND(J,2,4)                                                  
         M= IN6(K)                                                      
      DO I=1,5                                                          
         D1(I  ,2,4)=+S2(I+31,M)-S1(I+35,J)                             
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..4,  K = 3  and  L = 4                   
Cjms                                                                    
      D4(  1,3,4)=+R5( 5, 1)-R4( 2, 1)                                  
C                                                                       
      D4(  2,3,4)=+R5( 8, 1)-R4( 4, 1)                                  
      D4(  3,3,4)=+R5( 9, 1)-R4( 5, 1)                                  
C                                                                       
      D4(  4,3,4)=+R5(12, 1)-R4( 7, 1)                                  
      D4(  5,3,4)=+R5(13, 1)-R4( 8, 1)                                  
      D4(  6,3,4)=+R5(14, 1)-R4( 9, 1)                                  
C                                                                       
      D4(  7,3,4)=+R5(17, 1)-R4(11, 1)                                  
      D4(  8,3,4)=+R5(18, 1)-R4(12, 1)                                  
      D4(  9,3,4)=+R5(19, 1)-R4(13, 1)                                  
      D4( 10,3,4)=+R5(20, 1)-R4(14, 1)                                  
C                                                                       
      D3(1,1,3,4)=+S4( 5, 5)-S3( 6, 2)                                  
      D3(2,1,3,4)=+S4( 6, 5)-S3( 7, 2)                                  
      D3(3,1,3,4)=+S4( 7, 5)-S3( 8, 2)                                  
C                                                                       
      D3(1,2,3,4)=+S4( 5, 8)-S3( 6, 4)                                  
      D3(2,2,3,4)=+S4( 6, 8)-S3( 7, 4)                                  
      D3(3,2,3,4)=+S4( 7, 8)-S3( 8, 4)                                  
      D3(1,3,3,4)=+S4( 5, 9)-S3( 6, 5)                                  
      D3(2,3,3,4)=+S4( 6, 9)-S3( 7, 5)                                  
      D3(3,3,3,4)=+S4( 7, 9)-S3( 8, 5)                                  
C                                                                       
      D3(1,4,3,4)=+S4( 5,12)-S3( 6, 7)                                  
      D3(2,4,3,4)=+S4( 6,12)-S3( 7, 7)                                  
      D3(3,4,3,4)=+S4( 7,12)-S3( 8, 7)                                  
      D3(1,5,3,4)=+S4( 5,13)-S3( 6, 8)                                  
      D3(2,5,3,4)=+S4( 6,13)-S3( 7, 8)                                  
      D3(3,5,3,4)=+S4( 7,13)-S3( 8, 8)                                  
      D3(1,6,3,4)=+S4( 5,14)-S3( 6, 9)                                  
      D3(2,6,3,4)=+S4( 6,14)-S3( 7, 9)                                  
      D3(3,6,3,4)=+S4( 7,14)-S3( 8, 9)                                  
C                                                                       
      D2(1,1,3,4)=+S3(18, 5)-S2(28, 4)                                  
      D2(2,1,3,4)=+S3(19, 5)-S2(29, 4)                                  
      D2(3,1,3,4)=+S3(20, 5)-S2(30, 4)                                  
      D2(4,1,3,4)=+S3(21, 5)-S2(31, 4)                                  
C                                                                       
      D2(1,2,3,4)=+S3(18, 8)-S2(28, 2)                                  
      D2(2,2,3,4)=+S3(19, 8)-S2(29, 2)                                  
      D2(3,2,3,4)=+S3(20, 8)-S2(30, 2)                                  
      D2(4,2,3,4)=+S3(21, 8)-S2(31, 2)                                  
      D2(1,3,3,4)=+S3(18, 9)-S2(28, 6)                                  
      D2(2,3,3,4)=+S3(19, 9)-S2(29, 6)                                  
      D2(3,3,3,4)=+S3(20, 9)-S2(30, 6)                                  
      D2(4,3,3,4)=+S3(21, 9)-S2(31, 6)                                  
C                                                                       
      DO I=1,5                                                          
         D1(I  ,3,4)=+S2(I+31,6)-S1(I+35,2)                             
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..4,  K = 4  and  L = 4                   
Cjms                                                                    
      D4(  1,4,4)=+R5( 6, 1)-R4( 3, 1)-R4( 3, 4)+R3( 1, 4)+R3( 1, 5)    
C                                                                       
      D4(  2,4,4)=+R5( 9, 1)-R4( 5, 1)-R4( 5, 4)+R3( 2, 4)+R3( 2, 5)    
      D4(  3,4,4)=+R5(10, 1)-R4( 6, 1)-R4( 6, 4)+R3( 3, 4)+R3( 3, 5)    
C                                                                       
      D4(  4,4,4)=+R5(13, 1)-R4( 8, 1)-R4( 8, 4)+R3( 4, 4)+R3( 4, 5)    
      D4(  5,4,4)=+R5(14, 1)-R4( 9, 1)-R4( 9, 4)+R3( 5, 4)+R3( 5, 5)    
      D4(  6,4,4)=+R5(15, 1)-R4(10, 1)-R4(10, 4)+R3( 6, 4)+R3( 6, 5)    
C                                                                       
      D4(  7,4,4)=+R5(18, 1)-R4(12, 1)-R4(12, 4)+R3( 7, 4)+R3( 7, 5)    
      D4(  8,4,4)=+R5(19, 1)-R4(13, 1)-R4(13, 4)+R3( 8, 4)+R3( 8, 5)    
      D4(  9,4,4)=+R5(20, 1)-R4(14, 1)-R4(14, 4)+R3( 9, 4)+R3( 9, 5)    
      D4( 10,4,4)=+R5(21, 1)-R4(15, 1)-R4(15, 4)+R3(10, 4)+R3(10, 5)    
C                                                                       
      D3(1,1,4,4)=+S4( 5, 6)-S3( 6, 3)-S3(15, 3)+S2(10, 1)+S2(13, 1)    
      D3(2,1,4,4)=+S4( 6, 6)-S3( 7, 3)-S3(16, 3)+S2(11, 1)+S2(14, 1)    
      D3(3,1,4,4)=+S4( 7, 6)-S3( 8, 3)-S3(17, 3)+S2(12, 1)+S2(15, 1)    
C                                                                       
      D3(1,2,4,4)=+S4( 5, 9)-S3( 6, 5)-S3(15, 5)+S2(10, 4)+S2(13, 4)    
      D3(2,2,4,4)=+S4( 6, 9)-S3( 7, 5)-S3(16, 5)+S2(11, 4)+S2(14, 4)    
      D3(3,2,4,4)=+S4( 7, 9)-S3( 8, 5)-S3(17, 5)+S2(12, 4)+S2(15, 4)    
      D3(1,3,4,4)=+S4( 5,10)-S3( 6, 6)-S3(15, 6)+S2(10, 5)+S2(13, 5)    
      D3(2,3,4,4)=+S4( 6,10)-S3( 7, 6)-S3(16, 6)+S2(11, 5)+S2(14, 5)    
      D3(3,3,4,4)=+S4( 7,10)-S3( 8, 6)-S3(17, 6)+S2(12, 5)+S2(15, 5)    
C                                                                       
      D3(1,4,4,4)=+S4( 5,13)-S3( 6, 8)-S3(15, 8)+S2(10, 2)+S2(13, 2)    
      D3(2,4,4,4)=+S4( 6,13)-S3( 7, 8)-S3(16, 8)+S2(11, 2)+S2(14, 2)    
      D3(3,4,4,4)=+S4( 7,13)-S3( 8, 8)-S3(17, 8)+S2(12, 2)+S2(15, 2)    
      D3(1,5,4,4)=+S4( 5,14)-S3( 6, 9)-S3(15, 9)+S2(10, 6)+S2(13, 6)    
      D3(2,5,4,4)=+S4( 6,14)-S3( 7, 9)-S3(16, 9)+S2(11, 6)+S2(14, 6)    
      D3(3,5,4,4)=+S4( 7,14)-S3( 8, 9)-S3(17, 9)+S2(12, 6)+S2(15, 6)    
      D3(1,6,4,4)=+S4( 5,15)-S3( 6,10)-S3(15,10)+S2(10, 3)+S2(13, 3)    
      D3(2,6,4,4)=+S4( 6,15)-S3( 7,10)-S3(16,10)+S2(11, 3)+S2(14, 3)    
      D3(3,6,4,4)=+S4( 7,15)-S3( 8,10)-S3(17,10)+S2(12, 3)+S2(15, 3)    
C                                                                       
      D2(1,1,4,4)=+S3(18, 6)-S2(24, 5)-S2(28, 5)+S1(13, 1)+S1(17, 1)    
      D2(2,1,4,4)=+S3(19, 6)-S2(25, 5)-S2(29, 5)+S1(14, 1)+S1(18, 1)    
      D2(3,1,4,4)=+S3(20, 6)-S2(26, 5)-S2(30, 5)+S1(15, 1)+S1(19, 1)    
      D2(4,1,4,4)=+S3(21, 6)-S2(27, 5)-S2(31, 5)+S1(16, 1)+S1(20, 1)    
C                                                                       
      D2(1,2,4,4)=+S3(18, 9)-S2(24, 6)-S2(28, 6)+S1(13, 2)+S1(17, 2)    
      D2(2,2,4,4)=+S3(19, 9)-S2(25, 6)-S2(29, 6)+S1(14, 2)+S1(18, 2)    
      D2(3,2,4,4)=+S3(20, 9)-S2(26, 6)-S2(30, 6)+S1(15, 2)+S1(19, 2)    
      D2(4,2,4,4)=+S3(21, 9)-S2(27, 6)-S2(31, 6)+S1(16, 2)+S1(20, 2)    
      D2(1,3,4,4)=+S3(18,10)-S2(24, 3)-S2(28, 3)+S1(13, 3)+S1(17, 3)    
      D2(2,3,4,4)=+S3(19,10)-S2(25, 3)-S2(29, 3)+S1(14, 3)+S1(18, 3)    
      D2(3,3,4,4)=+S3(20,10)-S2(26, 3)-S2(30, 3)+S1(15, 3)+S1(19, 3)    
      D2(4,3,4,4)=+S3(21,10)-S2(27, 3)-S2(31, 3)+S1(16, 3)+S1(20, 3)    
C                                                                       
      DO I=1,5                                                          
         D1(I  ,4,4)=+S2(I+31,3)-S1(I+30,3)-S1(I+35,3)+R0(I+15)+R0(I+20)
      ENDDO                                                             
C                                                                       
      CALL FIJKL4(KX,LX,F,QX,QZ,D1,D2,D3,D4)                            
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2S   *DECK MCDV17                                           
C>                                                                      
C>    @brief   DDDS case                                                
C>                                                                      
C>    @details integration of a DDDS case                               
C>                                                                      
      SUBROUTINE MCDV17(F,QX,QZ)                                        
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      DIMENSION F(6,6,6,6)                                              
C                                                                       
      COMMON /KI4   / R0(   25),R1( 3,40),R2( 6,56),R3(10,52),          
     *                R4(15,42),R5(21,24),R6(28,12),R7(36, 4),R8(   45) 
C$omp threadprivate(/KI4/)
C                                                                       
      DIMENSION                 S1(13, 3),S2(17, 6),S3(12,10),          
     *                S4( 8,15),S5( 3,21)                               
C                                                                       
      PARAMETER (F02=2.0D+00)                                           
Cjms                                                                    
C     Simplified calculation of F(I,J,K,L) for cases                    
C     where:  I = 1..6,  J = 1..6,  K = 1..KX  and  L = 1..LX           
C     using auxiliary arrays E1, E2, E3, E4 and E5                      
Cjms                                                                    
      PARAMETER (KX=6)                                                  
      PARAMETER (LX=1)                                                  
      DIMENSION       E1(  5,KX,LX),E2(4,3,KX,LX),E3(4,6,KX,LX)         
      DIMENSION       E4(2,10,KX,LX),E5( 15,KX,LX)                      
C                                                                       
      DIMENSION INS(KX,LX),INI(KX,LX),IND(15,KX,LX)                     
      DIMENSION IN6(6)                                                  
      DATA INS/ 0, 1, 3, 0, 1, 2/                                       
      DATA INI/ 0, 2, 2, 1, 1, 2/                                       
      DATA IN6/ 1, 4, 5, 2, 6, 3/                                       
C                                                                       
      DO L=1,LX                                                         
         DO K=1,KX                                                      
            IJ= 0                                                       
            IN= INS(K,L)                                                
            DO I=1,5                                                    
               IN= IN+INI(K,L)                                          
               DO J=1,I                                                 
                  IJ= IJ+1                                              
                  IN= IN+1                                              
                  IND(IJ,K,L)= IN                                       
               ENDDO                                                    
            ENDDO                                                       
         ENDDO                                                          
      ENDDO                                                             
      DO 101 J=1,13                                                     
         DO 101 I=1, 3                                                  
  101 S1(J,I)= R1(I,J)                                                  
      DO 102 J=1,17                                                     
         DO 102 I=1, 6                                                  
  102 S2(J,I)= R2(I,J)                                                  
      DO 103 J=1,12                                                     
         DO 103 I=1,10                                                  
  103 S3(J,I)= R3(I,J)                                                  
      DO 104 J=1, 8                                                     
         DO 104 I=1,15                                                  
  104 S4(J,I)= R4(I,J)                                                  
      DO 105 J=1, 3                                                     
         DO 105 I=1,21                                                  
  105 S5(J,I)= R5(I,J)                                                  
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 1  and  L = 1                   
Cjms                                                                    
      DO J=1,15                                                         
         E5(  J,1,1)=+R6( J,1)+R4( J,1)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,10                                                         
         E4(1,J,1,1)=+S5( 2,J)+S3( 1,J)                                 
         E4(2,J,1,1)=+S5( 3,J)+S3( 2,J)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         M= IN6(J)                                                      
         E3(1,J,1,1)=+S4( 5,J)+S2( 1,M)                                 
         E3(2,J,1,1)=+S4( 6,J)+S2( 2,M)                                 
         E3(3,J,1,1)=+S4( 7,J)+S2( 3,M)                                 
         E3(4,J,1,1)=+S4( 8,J)+S2( 4,M)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         E2(1,J,1,1)=+S3( 9,J)+S1( 1,J)                                 
         E2(2,J,1,1)=+S3(10,J)+S1( 2,J)                                 
         E2(3,J,1,1)=+S3(11,J)+S1( 3,J)                                 
         E2(4,J,1,1)=+S3(12,J)+S1( 4,J)                                 
      ENDDO                                                             
C                                                                       
         J=1                                                            
      DO I=1,5                                                          
         E1(I  ,1,1)=+S2(I+12,J)+R0(I)                                  
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 2  and  L = 1                   
Cjms                                                                    
      DO J=1,15                                                         
         K= IND(J,2,1)                                                  
         E5(  J,2,1)=+R6( K,1)+R4( J,1)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,10                                                         
         K= IND(J,2,1)                                                  
         E4(1,J,2,1)=+S5( 2,K)+S3( 1,J)                                 
         E4(2,J,2,1)=+S5( 3,K)+S3( 2,J)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         K= IND(J,2,1)                                                  
         M= IN6(J)                                                      
         E3(1,J,2,1)=+S4( 5,K)+S2( 1,M)                                 
         E3(2,J,2,1)=+S4( 6,K)+S2( 2,M)                                 
         E3(3,J,2,1)=+S4( 7,K)+S2( 3,M)                                 
         E3(4,J,2,1)=+S4( 8,K)+S2( 4,M)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         K= IND(J,2,1)                                                  
         E2(1,J,2,1)=+S3( 9,K)+S1( 1,J)                                 
         E2(2,J,2,1)=+S3(10,K)+S1( 2,J)                                 
         E2(3,J,2,1)=+S3(11,K)+S1( 3,J)                                 
         E2(4,J,2,1)=+S3(12,K)+S1( 4,J)                                 
      ENDDO                                                             
C                                                                       
         J=1                                                            
         K= IND(J,2,1)                                                  
         M= IN6(K)                                                      
      DO I=1,5                                                          
         E1(I  ,2,1)=+S2(I+12,M)+R0(I)                                  
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 3  and  L = 1                   
Cjms                                                                    
      J= 0                                                              
      DO JJ=1,5                                                         
         DO II=1,JJ                                                     
            J= J+1                                                      
            K= IND(J,3,1)                                               
            L= K-2-JJ                                                   
C                                                                       
      E5(  J,3,1)=+R6(K,1)-R5(L,1)*F02+S4(1,J)+S4(2,J)                  
C                                                                       
            IF(JJ.GT. 4) GO TO 131                                      
      E4(1,J,3,1)=+S5(2,K)-S4(3,L)*F02+S3(1,J)+S3(3,J)                  
      E4(2,J,3,1)=+S5(3,K)-S4(4,L)*F02+S3(2,J)+S3(4,J)                  
C                                                                       
            IF(JJ.GT. 3) GO TO 131                                      
            M= IN6(J)                                                   
            DO I=1,4                                                    
      E3(I,J,3,1)=+S4(I+ 4,K)-S3(I+ 4,L)*F02+S2(I   ,M)+S2(I+ 4,M)      
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 2) GO TO 131                                      
            M= IN6(L)                                                   
            DO I=1,4                                                    
      E2(I,J,3,1)=+S3(I+ 8,K)-S2(I+ 8,M)*F02+S1(I   ,J)+S1(I+ 4,J)      
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 1) GO TO 131                                      
            M= IN6(K)                                                   
            DO I=1,5                                                    
      E1(I  ,3,1)=+S2(I+12,M)-S1(I+ 8,L)*F02+R0(I   )+R0(I+ 5)          
            ENDDO                                                       
  131       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 4  and  L = 1                   
Cjms                                                                    
      DO J=1,15                                                         
         K= IND(J,4,1)                                                  
         E5(  J,4,1)=+R6( K,1)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,10                                                         
         K= IND(J,4,1)                                                  
         E4(1,J,4,1)=+S5( 2,K)                                          
         E4(2,J,4,1)=+S5( 3,K)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         K= IND(J,4,1)                                                  
         E3(1,J,4,1)=+S4( 5,K)                                          
         E3(2,J,4,1)=+S4( 6,K)                                          
         E3(3,J,4,1)=+S4( 7,K)                                          
         E3(4,J,4,1)=+S4( 8,K)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         K= IND(J,4,1)                                                  
         E2(1,J,4,1)=+S3( 9,K)                                          
         E2(2,J,4,1)=+S3(10,K)                                          
         E2(3,J,4,1)=+S3(11,K)                                          
         E2(4,J,4,1)=+S3(12,K)                                          
      ENDDO                                                             
C                                                                       
         J=1                                                            
         K= IND(J,4,1)                                                  
         M= IN6(K)                                                      
      DO I=1,5                                                          
         E1(I  ,4,1)=+S2(I+12,M)                                        
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 5  and  L = 1                   
Cjms                                                                    
      DO J=1,15                                                         
         K= IND(J,5,1)                                                  
         E5(  J,5,1)=+R6( K,1)-R5( J,1)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,10                                                         
         K= IND(J,5,1)                                                  
         E4(1,J,5,1)=+S5( 2,K)-S4( 3,J)                                 
         E4(2,J,5,1)=+S5( 3,K)-S4( 4,J)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         K= IND(J,5,1)                                                  
         E3(1,J,5,1)=+S4( 5,K)-S3( 5,J)                                 
         E3(2,J,5,1)=+S4( 6,K)-S3( 6,J)                                 
         E3(3,J,5,1)=+S4( 7,K)-S3( 7,J)                                 
         E3(4,J,5,1)=+S4( 8,K)-S3( 8,J)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         K= IND(J,5,1)                                                  
         M= IN6(J)                                                      
         E2(1,J,5,1)=+S3( 9,K)-S2( 9,M)                                 
         E2(2,J,5,1)=+S3(10,K)-S2(10,M)                                 
         E2(3,J,5,1)=+S3(11,K)-S2(11,M)                                 
         E2(4,J,5,1)=+S3(12,K)-S2(12,M)                                 
      ENDDO                                                             
C                                                                       
         J=1                                                            
         K= IND(J,5,1)                                                  
         M= IN6(K)                                                      
      DO I=1,5                                                          
         E1(I  ,5,1)=+S2(I+12,M)-S1(I+ 8,J)                             
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 6  and  L = 1                   
Cjms                                                                    
      J= 0                                                              
      DO JJ=1,5                                                         
         DO II=1,JJ                                                     
            J= J+1                                                      
            K= IND(J,6,1)                                               
            L= K-2-JJ                                                   
C                                                                       
      E5(  J,6,1)=+R6(K,1)-R5(L,1)                                      
C                                                                       
            IF(JJ.GT. 4) GO TO 161                                      
      E4(1,J,6,1)=+S5(2,K)-S4(3,L)                                      
      E4(2,J,6,1)=+S5(3,K)-S4(4,L)                                      
C                                                                       
            IF(JJ.GT. 3) GO TO 161                                      
            DO I=1,4                                                    
      E3(I,J,6,1)=+S4(I+ 4,K)-S3(I+ 4,L)                                
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 2) GO TO 161                                      
            M= IN6(L)                                                   
            DO I=1,4                                                    
      E2(I,J,6,1)=+S3(I+ 8,K)-S2(I+ 8,M)                                
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 1) GO TO 161                                      
            M= IN6(K)                                                   
            DO I=1,5                                                    
      E1(I  ,6,1)=+S2(I+12,M)-S1(I+ 8,L)                                
            ENDDO                                                       
  161       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
C                                                                       
      CALL FIJKL5(KX,LX,F,QX,QZ,E1,E2,E3,E4,E5)                         
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2S   *DECK MCDV18                                           
C>                                                                      
C>    @brief   DDPP case                                                
C>                                                                      
C>    @details integration of a DDPP case                               
C>                                                                      
      SUBROUTINE MCDV18(F,QX,QZ)                                        
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      DIMENSION F(6,6,6,6)                                              
C                                                                       
      COMMON /KI4   / R0(   25),R1( 3,40),R2( 6,56),R3(10,52),          
     *                R4(15,42),R5(21,24),R6(28,12),R7(36, 4),R8(   45) 
C$omp threadprivate(/KI4/)
C                                                                       
      DIMENSION                 S1(40, 3),S2(41, 6),S3(30,10),          
     *                S4(17,15),S5( 6,21)                               
Cjms                                                                    
C     Simplified calculation of F(I,J,K,L) for cases                    
C     where:  I = 1..6,  J = 1..6,  K = 1..KX  and  L = 1..LX           
C     using auxiliary arrays E1, E2, E3, E4 and E5                      
Cjms                                                                    
      PARAMETER (KX=4)                                                  
      PARAMETER (LX=4)                                                  
      DIMENSION       E1(  5,KX,LX),E2(4,3,KX,LX),E3(4,6,KX,LX)         
      DIMENSION       E4(2,10,KX,LX),E5( 15,KX,LX)                      
C                                                                       
      DIMENSION INS(KX,LX),INI(KX,LX),IND(15,KX,LX)                     
      DIMENSION IN6(6)                                                  
      DATA INS/ 0, 0, 0, 1,  0, 0, 0, 1,  0, 0, 1, 2,  1, 1, 2, 3/      
      DATA INI/ 0, 0, 1, 1,  0, 0, 1, 1,  1, 1, 2, 2,  1, 1, 2, 2/      
      DATA IN6/ 1, 4, 5, 2, 6, 3/                                       
C                                                                       
      DO L=1,LX                                                         
         DO K=1,KX                                                      
            IJ= 0                                                       
            IN= INS(K,L)                                                
            DO I=1,5                                                    
               IN= IN+INI(K,L)                                          
               DO J=1,I                                                 
                  IJ= IJ+1                                              
                  IN= IN+1                                              
                  IND(IJ,K,L)= IN                                       
               ENDDO                                                    
            ENDDO                                                       
         ENDDO                                                          
      ENDDO                                                             
      DO 101 J=1,40                                                     
         DO 101 I=1, 3                                                  
  101 S1(J,I)= R1(I,J)                                                  
      DO 102 J=1,41                                                     
         DO 102 I=1, 6                                                  
  102 S2(J,I)= R2(I,J)                                                  
      DO 103 J=1,30                                                     
         DO 103 I=1,10                                                  
  103 S3(J,I)= R3(I,J)                                                  
      DO 104 J=1,17                                                     
         DO 104 I=1,15                                                  
  104 S4(J,I)= R4(I,J)                                                  
      DO 105 J=1, 6                                                     
         DO 105 I=1,21                                                  
  105 S5(J,I)= R5(I,J)                                                  
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 1  and  L = 1                   
Cjms                                                                    
      DO J=1,15                                                         
         E5(  J,1,1)=+R4( J,1)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,10                                                         
         E4(1,J,1,1)=+S3( 1,J)                                          
         E4(2,J,1,1)=+S3( 2,J)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         M= IN6(J)                                                      
         E3(1,J,1,1)=+S2( 1,M)                                          
         E3(2,J,1,1)=+S2( 2,M)                                          
         E3(3,J,1,1)=+S2( 3,M)                                          
         E3(4,J,1,1)=+S2( 4,M)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         E2(1,J,1,1)=+S1( 1,J)                                          
         E2(2,J,1,1)=+S1( 2,J)                                          
         E2(3,J,1,1)=+S1( 3,J)                                          
         E2(4,J,1,1)=+S1( 4,J)                                          
      ENDDO                                                             
C                                                                       
      DO I=1,5                                                          
         E1(I  ,1,1)=+R0(I)                                             
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 2  and  L = 1                   
Cjms                                                                    
      DO J=1,15                                                         
         E5(  J,2,1)=-R5( J,2)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,10                                                         
         E4(1,J,2,1)=-S4( 8,J)                                          
         E4(2,J,2,1)=-S4( 9,J)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         E3(1,J,2,1)=-S3(11,J)                                          
         E3(2,J,2,1)=-S3(12,J)                                          
         E3(3,J,2,1)=-S3(13,J)                                          
         E3(4,J,2,1)=-S3(14,J)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         M= IN6(J)                                                      
         E2(1,J,2,1)=-S2(21,M)                                          
         E2(2,J,2,1)=-S2(22,M)                                          
         E2(3,J,2,1)=-S2(23,M)                                          
         E2(4,J,2,1)=-S2(24,M)                                          
      ENDDO                                                             
C                                                                       
         J=1                                                            
      DO I=1,5                                                          
         E1(I  ,2,1)=-S1(I+20,J)                                        
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 3  and  L = 1                   
Cjms                                                                    
      DO J=1,15                                                         
         K= IND(J,3,1)                                                  
         E5(  J,3,1)=-R5( K,2)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,10                                                         
         K= IND(J,3,1)                                                  
         E4(1,J,3,1)=-S4( 8,K)                                          
         E4(2,J,3,1)=-S4( 9,K)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         K= IND(J,3,1)                                                  
         E3(1,J,3,1)=-S3(11,K)                                          
         E3(2,J,3,1)=-S3(12,K)                                          
         E3(3,J,3,1)=-S3(13,K)                                          
         E3(4,J,3,1)=-S3(14,K)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         K= IND(J,3,1)                                                  
         M= IN6(K)                                                      
         E2(1,J,3,1)=-S2(21,M)                                          
         E2(2,J,3,1)=-S2(22,M)                                          
         E2(3,J,3,1)=-S2(23,M)                                          
         E2(4,J,3,1)=-S2(24,M)                                          
      ENDDO                                                             
C                                                                       
         J=1                                                            
         K= IND(J,3,1)                                                  
      DO I=1,5                                                          
         E1(I  ,3,1)=-S1(I+20,K)                                        
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 4  and  L = 1                   
Cjms                                                                    
      DO J=1,15                                                         
         K= IND(J,4,1)                                                  
         E5(  J,4,1)=-R5( K,2)+R4( J,2)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,10                                                         
         K= IND(J,4,1)                                                  
         E4(1,J,4,1)=-S4( 8,K)+S3( 3,J)                                 
         E4(2,J,4,1)=-S4( 9,K)+S3( 4,J)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         K= IND(J,4,1)                                                  
         M= IN6(J)                                                      
         E3(1,J,4,1)=-S3(11,K)+S2( 5,M)                                 
         E3(2,J,4,1)=-S3(12,K)+S2( 6,M)                                 
         E3(3,J,4,1)=-S3(13,K)+S2( 7,M)                                 
         E3(4,J,4,1)=-S3(14,K)+S2( 8,M)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         K= IND(J,4,1)                                                  
         M= IN6(K)                                                      
         E2(1,J,4,1)=-S2(21,M)+S1( 5,J)                                 
         E2(2,J,4,1)=-S2(22,M)+S1( 6,J)                                 
         E2(3,J,4,1)=-S2(23,M)+S1( 7,J)                                 
         E2(4,J,4,1)=-S2(24,M)+S1( 8,J)                                 
      ENDDO                                                             
C                                                                       
         J=1                                                            
         K= IND(J,4,1)                                                  
      DO I=1,5                                                          
         E1(I  ,4,1)=-S1(I+20,K)+R0(I+ 5)                               
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 1  and  L = 2                   
Cjms                                                                    
      DO J=1,15                                                         
         E5(  J,1,2)=-R5( J,3)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,10                                                         
         E4(1,J,1,2)=-S4(10,J)                                          
         E4(2,J,1,2)=-S4(11,J)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         E3(1,J,1,2)=-S3(15,J)                                          
         E3(2,J,1,2)=-S3(16,J)                                          
         E3(3,J,1,2)=-S3(17,J)                                          
         E3(4,J,1,2)=-S3(18,J)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         M= IN6(J)                                                      
         E2(1,J,1,2)=-S2(25,M)                                          
         E2(2,J,1,2)=-S2(26,M)                                          
         E2(3,J,1,2)=-S2(27,M)                                          
         E2(4,J,1,2)=-S2(28,M)                                          
      ENDDO                                                             
C                                                                       
         J=1                                                            
      DO I=1,5                                                          
         E1(I  ,1,2)=-S1(I+25,J)                                        
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 2  and  L = 2                   
Cjms                                                                    
      DO J=1,15                                                         
         E5(  J,2,2)=+R6( J,1)+R4( J,5)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,10                                                         
         E4(1,J,2,2)=+S5( 5,J)+S3( 9,J)                                 
         E4(2,J,2,2)=+S5( 6,J)+S3(10,J)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         M= IN6(J)                                                      
         E3(1,J,2,2)=+S4(14,J)+S2(17,M)                                 
         E3(2,J,2,2)=+S4(15,J)+S2(18,M)                                 
         E3(3,J,2,2)=+S4(16,J)+S2(19,M)                                 
         E3(4,J,2,2)=+S4(17,J)+S2(20,M)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         E2(1,J,2,2)=+S3(27,J)+S1(17,J)                                 
         E2(2,J,2,2)=+S3(28,J)+S1(18,J)                                 
         E2(3,J,2,2)=+S3(29,J)+S1(19,J)                                 
         E2(4,J,2,2)=+S3(30,J)+S1(20,J)                                 
      ENDDO                                                             
C                                                                       
         J=1                                                            
      DO I=1,5                                                          
         E1(I  ,2,2)=+S2(I+36,J)+R0(I+20)                               
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 3  and  L = 2                   
Cjms                                                                    
      DO J=1,15                                                         
         K= IND(J,3,2)                                                  
         E5(  J,3,2)=+R6( K,1)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,10                                                         
         K= IND(J,3,2)                                                  
         E4(1,J,3,2)=+S5( 5,K)                                          
         E4(2,J,3,2)=+S5( 6,K)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         K= IND(J,3,2)                                                  
         E3(1,J,3,2)=+S4(14,K)                                          
         E3(2,J,3,2)=+S4(15,K)                                          
         E3(3,J,3,2)=+S4(16,K)                                          
         E3(4,J,3,2)=+S4(17,K)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         K= IND(J,3,2)                                                  
         E2(1,J,3,2)=+S3(27,K)                                          
         E2(2,J,3,2)=+S3(28,K)                                          
         E2(3,J,3,2)=+S3(29,K)                                          
         E2(4,J,3,2)=+S3(30,K)                                          
      ENDDO                                                             
C                                                                       
         J=1                                                            
         K= IND(J,3,2)                                                  
         M= IN6(K)                                                      
      DO I=1,5                                                          
         E1(I  ,3,2)=+S2(I+36,M)                                        
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 4  and  L = 2                   
Cjms                                                                    
      DO J=1,15                                                         
         K= IND(J,4,2)                                                  
         E5(  J,4,2)=+R6( K,1)-R5( J,4)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,10                                                         
         K= IND(J,4,2)                                                  
         E4(1,J,4,2)=+S5( 5,K)-S4(12,J)                                 
         E4(2,J,4,2)=+S5( 6,K)-S4(13,J)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         K= IND(J,4,2)                                                  
         E3(1,J,4,2)=+S4(14,K)-S3(19,J)                                 
         E3(2,J,4,2)=+S4(15,K)-S3(20,J)                                 
         E3(3,J,4,2)=+S4(16,K)-S3(21,J)                                 
         E3(4,J,4,2)=+S4(17,K)-S3(22,J)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         K= IND(J,4,2)                                                  
         M= IN6(J)                                                      
         E2(1,J,4,2)=+S3(27,K)-S2(29,M)                                 
         E2(2,J,4,2)=+S3(28,K)-S2(30,M)                                 
         E2(3,J,4,2)=+S3(29,K)-S2(31,M)                                 
         E2(4,J,4,2)=+S3(30,K)-S2(32,M)                                 
      ENDDO                                                             
C                                                                       
         J=1                                                            
         K= IND(J,4,2)                                                  
         M= IN6(K)                                                      
      DO I=1,5                                                          
         E1(I  ,4,2)=+S2(I+36,M)-S1(I+30,J)                             
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 1  and  L = 3                   
Cjms                                                                    
      DO J=1,15                                                         
         K= IND(J,1,3)                                                  
         E5(  J,1,3)=-R5( K,3)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,10                                                         
         K= IND(J,1,3)                                                  
         E4(1,J,1,3)=-S4(10,K)                                          
         E4(2,J,1,3)=-S4(11,K)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         K= IND(J,1,3)                                                  
         E3(1,J,1,3)=-S3(15,K)                                          
         E3(2,J,1,3)=-S3(16,K)                                          
         E3(3,J,1,3)=-S3(17,K)                                          
         E3(4,J,1,3)=-S3(18,K)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         K= IND(J,1,3)                                                  
         M= IN6(K)                                                      
         E2(1,J,1,3)=-S2(25,M)                                          
         E2(2,J,1,3)=-S2(26,M)                                          
         E2(3,J,1,3)=-S2(27,M)                                          
         E2(4,J,1,3)=-S2(28,M)                                          
      ENDDO                                                             
C                                                                       
         J=1                                                            
         K= IND(J,1,3)                                                  
      DO I=1,5                                                          
         E1(I  ,1,3)=-S1(I+25,K)                                        
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 2  and  L = 3                   
C                                                                       
C     next statements commented because F(I,J,2,3)                      
C     will be obtained in FIJKL5 from   F(I,J,3,2)                      
Cjms                                                                    
C     DO 325 J=1,15                                                     
C 325 E5(  J,2,3)= E5(  J,3,2)                                          
C     DO 324 J=1,10                                                     
C        DO 324 I=1,2                                                   
C 324 E4(I,J,2,3)= E4(I,J,3,2)                                          
C     DO 323 J=1,6                                                      
C        DO 323 I=1,4                                                   
C 323 E3(I,J,2,3)= E3(I,J,3,2)                                          
C     DO 322 J=1,3                                                      
C        DO 322 I=1,4                                                   
C 322 E2(I,J,2,3)= E2(I,J,3,2)                                          
C        DO 321 I=1,5                                                   
C 321 E1(I  ,2,3)= E1(I  ,3,2)                                          
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 3  and  L = 3                   
Cjms                                                                    
      DO J=1,15                                                         
         K= IND(J,3,3)                                                  
         E5(  J,3,3)=+R6( K,1)+R4( J,5)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,10                                                         
         K= IND(J,3,3)                                                  
         E4(1,J,3,3)=+S5( 5,K)+S3( 9,J)                                 
         E4(2,J,3,3)=+S5( 6,K)+S3(10,J)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         K= IND(J,3,3)                                                  
         M= IN6(J)                                                      
         E3(1,J,3,3)=+S4(14,K)+S2(17,M)                                 
         E3(2,J,3,3)=+S4(15,K)+S2(18,M)                                 
         E3(3,J,3,3)=+S4(16,K)+S2(19,M)                                 
         E3(4,J,3,3)=+S4(17,K)+S2(20,M)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         K= IND(J,3,3)                                                  
         E2(1,J,3,3)=+S3(27,K)+S1(17,J)                                 
         E2(2,J,3,3)=+S3(28,K)+S1(18,J)                                 
         E2(3,J,3,3)=+S3(29,K)+S1(19,J)                                 
         E2(4,J,3,3)=+S3(30,K)+S1(20,J)                                 
      ENDDO                                                             
C                                                                       
         J=1                                                            
         K= IND(J,3,3)                                                  
         M= IN6(K)                                                      
      DO I=1,5                                                          
         E1(I  ,3,3)=+S2(I+36,M)+R0(I+20)                               
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 4  and  L = 3                   
Cjms                                                                    
      J= 0                                                              
      DO JJ=1,5                                                         
         DO II=1,JJ                                                     
            J= J+1                                                      
            K= IND(J,4,3)                                               
            L= K-2-JJ                                                   
C                                                                       
      E5(  J,4,3)=+R6(K,1)-R5(L,4)                                      
C                                                                       
            IF(JJ.GT. 4) GO TO 343                                      
      E4(1,J,4,3)=+S5( 5,K)-S4(12,L)                                    
      E4(2,J,4,3)=+S5( 6,K)-S4(13,L)                                    
C                                                                       
            IF(JJ.GT. 3) GO TO 343                                      
            DO I=1,4                                                    
      E3(I,J,4,3)=+S4(I+13,K)-S3(I+18,L)                                
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 2) GO TO 343                                      
            M= IN6(L)                                                   
            DO I=1,4                                                    
      E2(I,J,4,3)=+S3(I+26,K)-S2(I+28,M)                                
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 1) GO TO 343                                      
            M= IN6(K)                                                   
            DO I=1,5                                                    
      E1(I  ,4,3)=+S2(I+36,M)-S1(I+30,L)                                
            ENDDO                                                       
  343       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 1  and  L = 4                   
Cjms                                                                    
      DO J=1,15                                                         
         K= IND(J,1,4)                                                  
         E5(  J,1,4)=-R5( K,3)+R4( J,3)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,10                                                         
         K= IND(J,1,4)                                                  
         E4(1,J,1,4)=-S4(10,K)+S3( 5,J)                                 
         E4(2,J,1,4)=-S4(11,K)+S3( 6,J)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         K= IND(J,1,4)                                                  
         M= IN6(J)                                                      
         E3(1,J,1,4)=-S3(15,K)+S2( 9,M)                                 
         E3(2,J,1,4)=-S3(16,K)+S2(10,M)                                 
         E3(3,J,1,4)=-S3(17,K)+S2(11,M)                                 
         E3(4,J,1,4)=-S3(18,K)+S2(12,M)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         K= IND(J,1,4)                                                  
         M= IN6(K)                                                      
         E2(1,J,1,4)=-S2(25,M)+S1( 9,J)                                 
         E2(2,J,1,4)=-S2(26,M)+S1(10,J)                                 
         E2(3,J,1,4)=-S2(27,M)+S1(11,J)                                 
         E2(4,J,1,4)=-S2(28,M)+S1(12,J)                                 
      ENDDO                                                             
C                                                                       
         J=1                                                            
         K= IND(J,1,4)                                                  
      DO I=1,5                                                          
         E1(I  ,1,4)=-S1(I+25,K)+R0(I+10)                               
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 2  and  L = 4                   
Cjms                                                                    
      DO J=1,15                                                         
         K= IND(J,2,4)                                                  
         E5(  J,2,4)=+R6( K,1)-R5( J,1)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,10                                                         
         K= IND(J,2,4)                                                  
         E4(1,J,2,4)=+S5( 5,K)-S4( 6,J)                                 
         E4(2,J,2,4)=+S5( 6,K)-S4( 7,J)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         K= IND(J,2,4)                                                  
         E3(1,J,2,4)=+S4(14,K)-S3(23,J)                                 
         E3(2,J,2,4)=+S4(15,K)-S3(24,J)                                 
         E3(3,J,2,4)=+S4(16,K)-S3(25,J)                                 
         E3(4,J,2,4)=+S4(17,K)-S3(26,J)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         K= IND(J,2,4)                                                  
         M= IN6(J)                                                      
         E2(1,J,2,4)=+S3(27,K)-S2(33,M)                                 
         E2(2,J,2,4)=+S3(28,K)-S2(34,M)                                 
         E2(3,J,2,4)=+S3(29,K)-S2(35,M)                                 
         E2(4,J,2,4)=+S3(30,K)-S2(36,M)                                 
      ENDDO                                                             
C                                                                       
         J=1                                                            
         K= IND(J,2,4)                                                  
         M= IN6(K)                                                      
      DO I=1,5                                                          
         E1(I  ,2,4)=+S2(I+36,M)-S1(I+35,J)                             
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 3  and  L = 4                   
Cjms                                                                    
      J= 0                                                              
      DO JJ=1,5                                                         
         DO II=1,JJ                                                     
            J= J+1                                                      
            K= IND(J,3,4)                                               
            L= K-2-JJ                                                   
C                                                                       
      E5(  J,3,4)=+R6(K,1)-R5(L,1)                                      
C                                                                       
            IF(JJ.GT. 4) GO TO 434                                      
      E4(1,J,3,4)=+S5( 5,K)-S4( 6,L)                                    
      E4(2,J,3,4)=+S5( 6,K)-S4( 7,L)                                    
C                                                                       
            IF(JJ.GT. 3) GO TO 434                                      
            DO I=1,4                                                    
      E3(I,J,3,4)=+S4(I+13,K)-S3(I+22,L)                                
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 2) GO TO 434                                      
            M= IN6(L)                                                   
            DO I=1,4                                                    
      E2(I,J,3,4)=+S3(I+26,K)-S2(I+32,M)                                
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 1) GO TO 434                                      
            M= IN6(K)                                                   
            DO I=1,5                                                    
      E1(I  ,3,4)=+S2(I+36,M)-S1(I+35,L)                                
            ENDDO                                                       
  434       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 4  and  L = 4                   
Cjms                                                                    
      J= 0                                                              
      DO JJ=1,5                                                         
         DO II=1,JJ                                                     
            J= J+1                                                      
            K= IND(J,4,4)                                               
            L= K-2-JJ                                                   
C                                                                       
      E5(  J,4,4)=+R6(K,1)-S5(1,L)-S5(4,L)+S4(4,J)+S4(5,J)              
C                                                                       
            IF(JJ.GT. 4) GO TO 444                                      
      E4(1,J,4,4)=+S5( 5,K)-S4( 6,L)-S4(12,L)+S3( 7,J)+S3( 9,J)         
      E4(2,J,4,4)=+S5( 6,K)-S4( 7,L)-S4(13,L)+S3( 8,J)+S3(10,J)         
C                                                                       
            IF(JJ.GT. 3) GO TO 444                                      
            M= IN6(J)                                                   
            DO I=1,4                                                    
      E3(I,J,4,4)=+S4(I+13,K)-S3(I+18,L)-S3(I+22,L)                     
     *                       +S2(I+12,M)+S2(I+16,M)                     
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 2) GO TO 444                                      
            M= IN6(L)                                                   
            DO I=1,4                                                    
      E2(I,J,4,4)=+S3(I+26,K)-S2(I+28,M)-S2(I+32,M)                     
     *                       +S1(I+12,J)+S1(I+16,J)                     
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 1) GO TO 444                                      
            M= IN6(K)                                                   
            DO I=1,5                                                    
      E1(I  ,4,4)=+S2(I+36,M)-S1(I+30,L)-S1(I+35,L)                     
     *                       +R0(I+15  )+R0(I+20  )                     
            ENDDO                                                       
  444       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
C                                                                       
      CALL FIJKL5(KX,LX,F,QX,QZ,E1,E2,E3,E4,E5)                         
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2S   *DECK MCDV19                                           
C>                                                                      
C>    @brief   DPDP case                                                
C>                                                                      
C>    @details integration of a DPDP case                               
C>                                                                      
      SUBROUTINE MCDV19(F,QX,QZ)                                        
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      DIMENSION F(6,6,6,6)                                              
C                                                                       
      COMMON /KI4   / R0(   25),R1( 3,40),R2( 6,56),R3(10,52),          
     *                R4(15,42),R5(21,24),R6(28,12),R7(36, 4),R8(   45) 
C$omp threadprivate(/KI4/)
C                                                                       
      DIMENSION                 S1(40, 3),S2(46, 6),S3(34,10),          
     *                S4(17,15),S5( 6,21)                               
C                                                                       
      PARAMETER (F02=2.0D+00)                                           
      PARAMETER (F03=3.0D+00)                                           
Cjms                                                                    
C     Simplified calculation of F(I,J,K,L) for cases                    
C     where:  I = 1..6,  J = 1..4,  K = 1..KX  and  L = 1..LX           
C     using auxiliary arrays D1, D2, D3 and D4                          
Cjms                                                                    
      PARAMETER (KX=6)                                                  
      PARAMETER (LX=4)                                                  
      DIMENSION       D1(  5,KX,LX),D2(4,3,KX,LX),D3(3,6,KX,LX)         
      DIMENSION       D4( 10,KX,LX)                                     
C                                                                       
      DIMENSION INS(KX,LX),INI(KX,LX),IND(10,KX,LX)                     
      DIMENSION IN6(6)                                                  
      DATA INS/ 0, 1, 3, 0, 1, 2,  0, 1, 3, 0, 1, 2,  0, 3, 5, 1, 2, 4, 
     *          1, 4, 6, 2, 3, 5/                                       
      DATA INI/ 0, 2, 2, 1, 1, 2,  0, 2, 2, 1, 1, 2,  1, 3, 3, 2, 2, 3, 
     *          1, 3, 3, 2, 2, 3/                                       
      DATA IN6/ 1, 4, 5, 2, 6, 3/                                       
C                                                                       
      DO L=1,LX                                                         
         DO K=1,KX                                                      
            IJ= 0                                                       
            IN= INS(K,L)                                                
            DO I=1,4                                                    
               IN= IN+INI(K,L)                                          
               DO J=1,I                                                 
                  IJ= IJ+1                                              
                  IN= IN+1                                              
                  IND(IJ,K,L)= IN                                       
               ENDDO                                                    
            ENDDO                                                       
         ENDDO                                                          
      ENDDO                                                             
      DO 101 J=1,40                                                     
         DO 101 I=1, 3                                                  
  101 S1(J,I)= R1(I,J)                                                  
      DO 102 J=1,46                                                     
         DO 102 I=1, 6                                                  
  102 S2(J,I)= R2(I,J)                                                  
      DO 103 J=1,34                                                     
         DO 103 I=1,10                                                  
  103 S3(J,I)= R3(I,J)                                                  
      DO 104 J=1,17                                                     
         DO 104 I=1,15                                                  
  104 S4(J,I)= R4(I,J)                                                  
      DO 105 J=1, 6                                                     
         DO 105 I=1,21                                                  
  105 S5(J,I)= R5(I,J)                                                  
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..4,  K = 1  and  L = 1                   
Cjms                                                                    
      DO J=1,10                                                         
         D4(  J,1,1)=+R5( J,2)+R3( J,1)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         M= IN6(J)                                                      
         D3(1,J,1,1)=+S4( 5,J)+S2( 1,M)                                 
         D3(2,J,1,1)=+S4( 6,J)+S2( 2,M)                                 
         D3(3,J,1,1)=+S4( 7,J)+S2( 3,M)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         D2(1,J,1,1)=+S3(18,J)+S1( 1,J)                                 
         D2(2,J,1,1)=+S3(19,J)+S1( 2,J)                                 
         D2(3,J,1,1)=+S3(20,J)+S1( 3,J)                                 
         D2(4,J,1,1)=+S3(21,J)+S1( 4,J)                                 
      ENDDO                                                             
C                                                                       
         J=1                                                            
      DO I=1,5                                                          
         D1(I  ,1,1)=+S2(I+31,J)+R0(I)                                  
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..4,  K = 2  and  L = 1                   
Cjms                                                                    
      DO J=1,10                                                         
         K= IND(J,2,1)                                                  
         D4(  J,2,1)=+R5( K,2)+R3( J,1)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         K= IND(J,2,1)                                                  
         M= IN6(J)                                                      
         D3(1,J,2,1)=+S4( 5,K)+S2( 1,M)                                 
         D3(2,J,2,1)=+S4( 6,K)+S2( 2,M)                                 
         D3(3,J,2,1)=+S4( 7,K)+S2( 3,M)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         K= IND(J,2,1)                                                  
         D2(1,J,2,1)=+S3(18,K)+S1( 1,J)                                 
         D2(2,J,2,1)=+S3(19,K)+S1( 2,J)                                 
         D2(3,J,2,1)=+S3(20,K)+S1( 3,J)                                 
         D2(4,J,2,1)=+S3(21,K)+S1( 4,J)                                 
      ENDDO                                                             
C                                                                       
         J=1                                                            
         K= IND(J,2,1)                                                  
         M= IN6(K)                                                      
      DO I=1,5                                                          
         D1(I  ,2,1)=+S2(I+31,M)+R0(I)                                  
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..4,  K = 3  and  L = 1                   
Cjms                                                                    
      J= 0                                                              
      DO JJ=1,4                                                         
         DO II=1,JJ                                                     
            J= J+1                                                      
            K= IND(J,3,1)                                               
            L= K-2-JJ                                                   
C                                                                       
      D4(  J,3,1)=+R5(K,2)-R4(L,2)*F02+S3(1,J)+S3(2,J)                  
C                                                                       
            IF(JJ.GT. 3) GO TO 131                                      
            M= IN6(J)                                                   
            DO I=1,3                                                    
      D3(I,J,3,1)=+S4(I+ 4,K)-S3(I+ 8,L)*F02+S2(I   ,M)+S2(I+ 3,M)      
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 2) GO TO 131                                      
            M= IN6(L)                                                   
            DO I=1,4                                                    
      D2(I,J,3,1)=+S3(I+17,K)-S2(I+15,M)*F02+S1(I   ,J)+S1(I+ 4,J)      
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 1) GO TO 131                                      
            M= IN6(K)                                                   
            DO I=1,5                                                    
      D1(I  ,3,1)=+S2(I+31,M)-S1(I+20,L)*F02+R0(I   )+R0(I+ 5)          
            ENDDO                                                       
  131       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..4,  K = 4  and  L = 1                   
Cjms                                                                    
      DO J=1,10                                                         
         K= IND(J,4,1)                                                  
         D4(  J,4,1)=+R5( K,2)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         K= IND(J,4,1)                                                  
         D3(1,J,4,1)=+S4( 5,K)                                          
         D3(2,J,4,1)=+S4( 6,K)                                          
         D3(3,J,4,1)=+S4( 7,K)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         K= IND(J,4,1)                                                  
         D2(1,J,4,1)=+S3(18,K)                                          
         D2(2,J,4,1)=+S3(19,K)                                          
         D2(3,J,4,1)=+S3(20,K)                                          
         D2(4,J,4,1)=+S3(21,K)                                          
      ENDDO                                                             
C                                                                       
         J=1                                                            
         K= IND(J,4,1)                                                  
         M= IN6(K)                                                      
      DO I=1,5                                                          
         D1(I  ,4,1)=+S2(I+31,M)                                        
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..4,  K = 5  and  L = 1                   
Cjms                                                                    
      DO J=1,10                                                         
         K= IND(J,5,1)                                                  
         D4(  J,5,1)=+R5( K,2)-R4( J,2)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         K= IND(J,5,1)                                                  
         D3(1,J,5,1)=+S4( 5,K)-S3( 9,J)                                 
         D3(2,J,5,1)=+S4( 6,K)-S3(10,J)                                 
         D3(3,J,5,1)=+S4( 7,K)-S3(11,J)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         K= IND(J,5,1)                                                  
         M= IN6(J)                                                      
         D2(1,J,5,1)=+S3(18,K)-S2(16,M)                                 
         D2(2,J,5,1)=+S3(19,K)-S2(17,M)                                 
         D2(3,J,5,1)=+S3(20,K)-S2(18,M)                                 
         D2(4,J,5,1)=+S3(21,K)-S2(19,M)                                 
      ENDDO                                                             
C                                                                       
         J=1                                                            
         K= IND(J,5,1)                                                  
         M= IN6(K)                                                      
      DO I=1,5                                                          
         D1(I  ,5,1)=+S2(I+31,M)-S1(I+20,J)                             
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..4,  K = 6  and  L = 1                   
Cjms                                                                    
      J= 0                                                              
      DO JJ=1,4                                                         
         DO II=1,JJ                                                     
            J= J+1                                                      
            K= IND(J,6,1)                                               
            L= K-2-JJ                                                   
C                                                                       
      D4(  J,6,1)=+R5(K,2)-R4(L,2)                                      
C                                                                       
            IF(JJ.GT. 3) GO TO 161                                      
            DO I=1,3                                                    
      D3(I,J,6,1)=+S4(I+ 4,K)-S3(I+ 8,L)                                
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 2) GO TO 161                                      
            M= IN6(L)                                                   
            DO I=1,4                                                    
      D2(I,J,6,1)=+S3(I+17,K)-S2(I+15,M)                                
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 1) GO TO 161                                      
            M= IN6(K)                                                   
            DO I=1,5                                                    
      D1(I  ,6,1)=+S2(I+31,M)-S1(I+20,L)                                
            ENDDO                                                       
  161       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..4,  K = 1  and  L = 2                   
Cjms                                                                    
      DO J=1,10                                                         
         D4(  J,1,2)=-R6( J,1)-R4( J,4)*F03                             
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         D3(1,J,1,2)=-S5( 4,J)-S3(15,J)*F03                             
         D3(2,J,1,2)=-S5( 5,J)-S3(16,J)*F03                             
         D3(3,J,1,2)=-S5( 6,J)-S3(17,J)*F03                             
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         M= IN6(J)                                                      
         D2(1,J,1,2)=-S4(14,J)-S2(24,M)*F03                             
         D2(2,J,1,2)=-S4(15,J)-S2(25,M)*F03                             
         D2(3,J,1,2)=-S4(16,J)-S2(26,M)*F03                             
         D2(4,J,1,2)=-S4(17,J)-S2(27,M)*F03                             
      ENDDO                                                             
C                                                                       
         J=1                                                            
      DO I=1,5                                                          
         D1(I  ,1,2)=-S3(I+29,J)-S1(I+30,J)*F03                         
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..4,  K = 2  and  L = 2                   
Cjms                                                                    
      DO J=1,10                                                         
         K= IND(J,2,2)                                                  
         D4(  J,2,2)=-R6( K,1)-R4( J,4)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         K= IND(J,2,2)                                                  
         D3(1,J,2,2)=-S5( 4,K)-S3(15,J)                                 
         D3(2,J,2,2)=-S5( 5,K)-S3(16,J)                                 
         D3(3,J,2,2)=-S5( 6,K)-S3(17,J)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         K= IND(J,2,2)                                                  
         M= IN6(J)                                                      
         D2(1,J,2,2)=-S4(14,K)-S2(24,M)                                 
         D2(2,J,2,2)=-S4(15,K)-S2(25,M)                                 
         D2(3,J,2,2)=-S4(16,K)-S2(26,M)                                 
         D2(4,J,2,2)=-S4(17,K)-S2(27,M)                                 
      ENDDO                                                             
C                                                                       
         J=1                                                            
         K= IND(J,2,2)                                                  
      DO I=1,5                                                          
         D1(I  ,2,2)=-S3(I+29,K)-S1(I+30,J)                             
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..4,  K = 3  and  L = 2                   
Cjms                                                                    
      J= 0                                                              
      DO JJ=1,4                                                         
         DO II=1,JJ                                                     
            J= J+1                                                      
            K= IND(J,3,2)                                               
            L= K-2-JJ                                                   
C                                                                       
      D4(  J,3,2)=-R6(K,1)+R5(L,3)*F02-S4(3,J)-S4(4,J)                  
C                                                                       
            IF(JJ.GT. 3) GO TO 232                                      
            DO I=1,3                                                    
      D3(I,J,3,2)=-S5(I+ 3,K)+S4(I+ 7,L)*F02-S3(I+11,J)-S3(I+14,J)      
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 2) GO TO 232                                      
            M= IN6(J)                                                   
            DO I=1,4                                                    
      D2(I,J,3,2)=-S4(I+13,K)+S3(I+21,L)*F02-S2(I+19,M)-S2(I+23,M)      
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 1) GO TO 232                                      
            M= IN6(L)                                                   
            DO I=1,5                                                    
      D1(I  ,3,2)=-S3(I+29,K)+S2(I+36,M)*F02-S1(I+25,J)-S1(I+30,J)      
            ENDDO                                                       
  232       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..4,  K = 4  and  L = 2                   
Cjms                                                                    
      DO J=1,10                                                         
         K= IND(J,4,2)                                                  
         D4(  J,4,2)=-R6( K,1)-R4( K,4)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         K= IND(J,4,2)                                                  
         D3(1,J,4,2)=-S5( 4,K)-S3(15,K)                                 
         D3(2,J,4,2)=-S5( 5,K)-S3(16,K)                                 
         D3(3,J,4,2)=-S5( 6,K)-S3(17,K)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         K= IND(J,4,2)                                                  
         M= IN6(K)                                                      
         D2(1,J,4,2)=-S4(14,K)-S2(24,M)                                 
         D2(2,J,4,2)=-S4(15,K)-S2(25,M)                                 
         D2(3,J,4,2)=-S4(16,K)-S2(26,M)                                 
         D2(4,J,4,2)=-S4(17,K)-S2(27,M)                                 
      ENDDO                                                             
C                                                                       
         J=1                                                            
         K= IND(J,4,2)                                                  
      DO I=1,5                                                          
         D1(I  ,4,2)=-S3(I+29,K)-S1(I+30,K)                             
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..4,  K = 5  and  L = 2                   
Cjms                                                                    
      DO J=1,10                                                         
         K= IND(J,5,2)                                                  
         D4(  J,5,2)=-R6( K,1)-R4( K,4)+R5( J,3)+R3( J,5)               
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         K= IND(J,5,2)                                                  
         M= IN6(J)                                                      
         D3(1,J,5,2)=-S5( 4,K)-S3(15,K)+S4( 8,J)+S2(13,M)               
         D3(2,J,5,2)=-S5( 5,K)-S3(16,K)+S4( 9,J)+S2(14,M)               
         D3(3,J,5,2)=-S5( 6,K)-S3(17,K)+S4(10,J)+S2(15,M)               
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         K= IND(J,5,2)                                                  
         M= IN6(K)                                                      
         D2(1,J,5,2)=-S4(14,K)-S2(24,M)+S3(22,J)+S1(17,J)               
         D2(2,J,5,2)=-S4(15,K)-S2(25,M)+S3(23,J)+S1(18,J)               
         D2(3,J,5,2)=-S4(16,K)-S2(26,M)+S3(24,J)+S1(19,J)               
         D2(4,J,5,2)=-S4(17,K)-S2(27,M)+S3(25,J)+S1(20,J)               
      ENDDO                                                             
C                                                                       
         J=1                                                            
         K= IND(J,5,2)                                                  
      DO I=1,5                                                          
         D1(I  ,5,2)=-S3(I+29,K)-S1(I+30,K)+S2(I+36,J)+R0(I+20)         
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..4,  K = 6  and  L = 2                   
Cjms                                                                    
      J= 0                                                              
      DO JJ=1,4                                                         
         DO II=1,JJ                                                     
            J= J+1                                                      
            K= IND(J,6,2)                                               
            L= K-2-JJ                                                   
C                                                                       
      D4(  J,6,2)=-R6(K,1)+R5(L,3)                                      
C                                                                       
            IF(JJ.GT. 3) GO TO 262                                      
            DO I=1,3                                                    
      D3(I,J,6,2)=-S5(I+ 3,K)+S4(I+ 7,L)                                
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 2) GO TO 262                                      
            DO I=1,4                                                    
      D2(I,J,6,2)=-S4(I+13,K)+S3(I+21,L)                                
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 1) GO TO 262                                      
            M= IN6(L)                                                   
            DO I=1,5                                                    
      D1(I  ,6,2)=-S3(I+29,K)+S2(I+36,M)                                
            ENDDO                                                       
  262       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..4,  K = 1  and  L = 3                   
C                                                                       
C     next statements commented because F(I,J,1,3)                      
C     will be obtained in FIJKL4 from   F(I,J,4,2)                      
Cjms                                                                    
C     DO 314 J=1,10                                                     
C 314 D4(  J,1,3)= D4(  J,4,2)                                          
C     DO 313 J=1,6                                                      
C        DO 313 I=1,3                                                   
C 313 D3(I,J,1,3)= D3(I,J,4,2)                                          
C     DO 312 J=1,3                                                      
C        DO 312 I=1,4                                                   
C 312 D2(I,J,1,3)= D2(I,J,4,2)                                          
C        DO 311 I=1,5                                                   
C 311 D1(I  ,1,3)= D1(I  ,4,2)                                          
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..4,  K = 2  and  L = 3                   
Cjms                                                                    
      J= 0                                                              
      DO JJ=1,4                                                         
         DO II=1,JJ                                                     
            J= J+1                                                      
            K= IND(J,2,3)                                               
            L= K-3-JJ-JJ                                                
C                                                                       
      D4(  J,2,3)=-R6(K,1)-R4(L,4)*F03                                  
C                                                                       
            IF(JJ.GT. 3) GO TO 323                                      
            DO I=1,3                                                    
      D3(I,J,2,3)=-S5(I+ 3,K)-S3(I+14,L)*F03                            
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 2) GO TO 323                                      
            M= IN6(L)                                                   
            DO I=1,4                                                    
      D2(I,J,2,3)=-S4(I+13,K)-S2(I+23,M)*F03                            
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 1) GO TO 323                                      
            DO I=1,5                                                    
      D1(I  ,2,3)=-S3(I+29,K)-S1(I+30,L)*F03                            
            ENDDO                                                       
  323       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..4,  K = 3  and  L = 3                   
Cjms                                                                    
      J= 0                                                              
      DO JJ=1,4                                                         
         DO II=1,JJ                                                     
            J= J+1                                                      
            K= IND(J,3,3)                                               
            L= K-3-JJ                                                   
            M= L-2-JJ                                                   
C                                                                       
      D4(  J,3,3)=-R6(K,1)+R5(L,3)*F02-S4(3,M)-S4(4,M)                  
C                                                                       
            IF(JJ.GT. 3) GO TO 333                                      
            DO I=1,3                                                    
      D3(I,J,3,3)=-S5(I+ 3,K)+S4(I+ 7,L)*F02-S3(I+11,M)-S3(I+14,M)      
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 2) GO TO 333                                      
            N= IN6(M)                                                   
            DO I=1,4                                                    
      D2(I,J,3,3)=-S4(I+13,K)+S3(I+21,L)*F02-S2(I+19,N)-S2(I+23,N)      
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 1) GO TO 333                                      
            N= IN6(L)                                                   
            DO I=1,5                                                    
      D1(I  ,3,3)=-S3(I+29,K)+S2(I+36,N)*F02-S1(I+25,M)-S1(I+30,M)      
            ENDDO                                                       
  333       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..4,  K = 4  and  L = 3                   
C                                                                       
C     next statements commented because F(I,J,4,3)                      
C     will be obtained in FIJKL4 from   F(I,J,2,2)                      
Cjms                                                                    
C     DO 344 J=1,10                                                     
C 344 D4(  J,4,3)= D4(  J,2,2)                                          
C     DO 343 J=1,6                                                      
C        DO 343 I=1,3                                                   
C 343 D3(I,J,4,3)= D3(I,J,2,2)                                          
C     DO 342 J=1,3                                                      
C        DO 342 I=1,4                                                   
C 342 D2(I,J,4,3)= D2(I,J,2,2)                                          
C        DO 341 I=1,5                                                   
C 341 D1(I  ,4,3)= D1(I  ,2,2)                                          
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..4,  K = 5  and  L = 3                   
C                                                                       
C     next statements commented because F(I,J,5,3)                      
C     will be obtained in FIJKL4 from   F(I,J,6,2)                      
Cjms                                                                    
C     DO 354 J=1,10                                                     
C 354 D4(  J,5,3)= D4(  J,6,2)                                          
C     DO 353 J=1,6                                                      
C        DO 353 I=1,3                                                   
C 353 D3(I,J,5,3)= D3(I,J,6,2)                                          
C     DO 352 J=1,3                                                      
C        DO 352 I=1,4                                                   
C 352 D2(I,J,5,3)= D2(I,J,6,2)                                          
C        DO 351 I=1,5                                                   
C 351 D1(I  ,5,3)= D1(I  ,6,2)                                          
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..4,  K = 6  and  L = 3                   
Cjms                                                                    
      J= 0                                                              
      DO JJ=1,4                                                         
         DO II=1,JJ                                                     
            J= J+1                                                      
            K= IND(J,6,3)                                               
            L= K-3-JJ                                                   
            M= L-JJ                                                     
C                                                                       
      D4(  J,6,3)=-R6(K,1)-R4(M,4)+R5(L,3)+R3(J,5)                      
C                                                                       
            IF(JJ.GT. 3) GO TO 363                                      
            N= IN6(J)                                                   
            DO I=1,3                                                    
      D3(I,J,6,3)=-S5(I+ 3,K)-S3(I+14,M)+S4(I+ 7,L)+S2(I+12,N)          
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 2) GO TO 363                                      
            N= IN6(M)                                                   
            DO I=1,4                                                    
      D2(I,J,6,3)=-S4(I+13,K)-S2(I+23,N)+S3(I+21,L)+S1(I+16,J)          
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 1) GO TO 363                                      
            N= IN6(L)                                                   
            DO I=1,5                                                    
      D1(I  ,6,3)=-S3(I+29,K)-S1(I+30,M)+S2(I+36,N)+R0(I+20)            
            ENDDO                                                       
  363       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..4,  K = 1  and  L = 4                   
Cjms                                                                    
      DO J=1,10                                                         
         K= IND(J,1,4)                                                  
         D4(  J,1,4)=-R6( K,1)-R4( K,4)+R5( J,1)+R3( J,4)               
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         K= IND(J,1,4)                                                  
         M= IN6(J)                                                      
         D3(1,J,1,4)=-S5( 4,K)-S3(15,K)+S4(11,J)+S2(10,M)               
         D3(2,J,1,4)=-S5( 5,K)-S3(16,K)+S4(12,J)+S2(11,M)               
         D3(3,J,1,4)=-S5( 6,K)-S3(17,K)+S4(13,J)+S2(12,M)               
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         K= IND(J,1,4)                                                  
         M= IN6(K)                                                      
         D2(1,J,1,4)=-S4(14,K)-S2(24,M)+S3(26,J)+S1(13,J)               
         D2(2,J,1,4)=-S4(15,K)-S2(25,M)+S3(27,J)+S1(14,J)               
         D2(3,J,1,4)=-S4(16,K)-S2(26,M)+S3(28,J)+S1(15,J)               
         D2(4,J,1,4)=-S4(17,K)-S2(27,M)+S3(29,J)+S1(16,J)               
      ENDDO                                                             
C                                                                       
         J=1                                                            
         K= IND(J,1,4)                                                  
      DO I=1,5                                                          
         D1(I  ,1,4)=-S3(I+29,K)-S1(I+30,K)+S2(I+41,J)+R0(I+15)         
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..4,  K = 2  and  L = 4                   
Cjms                                                                    
      J= 0                                                              
      DO JJ=1,4                                                         
         DO II=1,JJ                                                     
            J= J+1                                                      
            K= IND(J,2,4)                                               
            L= K-3-JJ                                                   
            M= L-JJ                                                     
C                                                                       
      D4(  J,2,4)=-R6(K,1)+R5(L,1)-R4(M,4)+R3(J,4)                      
C                                                                       
            IF(JJ.GT. 3) GO TO 424                                      
            N= IN6(J)                                                   
            DO I=1,3                                                    
      D3(I,J,2,4)=-S5(I+ 3,K)+S4(I+10,L)-S3(I+14,M)+S2(I+ 9,N)          
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 2) GO TO 424                                      
            N= IN6(M)                                                   
            DO I=1,4                                                    
      D2(I,J,2,4)=-S4(I+13,K)+S3(I+25,L)-S2(I+23,N)+S1(I+12,J)          
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 1) GO TO 424                                      
            N= IN6(L)                                                   
            DO I=1,5                                                    
      D1(I  ,2,4)=-S3(I+29,K)+S2(I+41,N)-S1(I+30,M)+R0(I+15)            
            ENDDO                                                       
  424       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..4,  K = 3  and  L = 4                   
Cjms                                                                    
      J= 0                                                              
      DO JJ=1,4                                                         
         DO II=1,JJ                                                     
            J= J+1                                                      
            K= IND(J,3,4)                                               
            L= K-3-JJ                                                   
            M= L-2-JJ                                                   
C                                                                       
      D4(  J,3,4)=-R6(K,1)    +S5(1,L)+S5(3,L)*F02                      
     *            -S4(1,M)*F02-S4(3,M)-S4(4,M)*F03                      
     *            +S3(3,J)    +S3(4,J)+S3(5,J)*F02                      
C                                                                       
            IF(JJ.GT. 3) GO TO 434                                      
            N= IN6(J)                                                   
            DO I=1,3                                                    
      D3(I,J,3,4)=-S5(I+ 3,K)    +S4(I+ 7,L)*F02+S4(I+10,L)             
     *            -S3(I+ 5,M)*F02-S3(I+11,M)    -S3(I+14,M)*F03         
     *            +S2(I+ 6,N)    +S2(I+ 9,N)    +S2(I+12,N)*F02         
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 2) GO TO 434                                      
            N= IN6(M)                                                   
            DO I=1,4                                                    
      D2(I,J,3,4)=-S4(I+13,K)+S3(I+21,L)*F02+S3(I+25,L)                 
     *            -S2(I+19,N)-S2(I+23,N)*F03-S2(I+27,N)*F02             
     *            +S1(I+ 8,J)+S1(I+12,J)    +S1(I+16,J)*F02             
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 1) GO TO 434                                      
            N= IN6(L)                                                   
            DO I=1,5                                                    
      D1(I  ,3,4)=-S3(I+29,K)+S2(I+36,N)*F02+S2(I+41,N)                 
     *            -S1(I+25,M)-S1(I+30,M)*F03-S1(I+35,M)*F02             
     *            +R0(I+10  )+R0(I+15  )    +R0(I+20  )*F02             
            ENDDO                                                       
  434       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..4,  K = 4  and  L = 4                   
Cjms                                                                    
      J= 0                                                              
      DO JJ=1,4                                                         
         DO II=1,JJ                                                     
            J= J+1                                                      
            K= IND(J,4,4)                                               
            L= K-2-JJ                                                   
C                                                                       
      D4(  J,4,4)=-R6(K,1)+R5(L,1)                                      
C                                                                       
            IF(JJ.GT. 3) GO TO 444                                      
            DO I=1,3                                                    
      D3(I,J,4,4)=-S5(I+ 3,K)+S4(I+10,L)                                
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 2) GO TO 444                                      
            DO I=1,4                                                    
      D2(I,J,4,4)=-S4(I+13,K)+S3(I+25,L)                                
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 1) GO TO 444                                      
            M= IN6(L)                                                   
            DO I=1,5                                                    
      D1(I  ,4,4)=-S3(I+29,K)+S2(I+41,M)                                
            ENDDO                                                       
  444       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..4,  K = 5  and  L = 4                   
Cjms                                                                    
      J= 0                                                              
      DO JJ=1,4                                                         
         DO II=1,JJ                                                     
            J= J+1                                                      
            K= IND(J,5,4)                                               
            L= K-2-JJ                                                   
C                                                                       
      D4(  J,5,4)=-R6(K,1)+S5(1,L)+S5(3,L)                              
     *                    -S4(1,J)-S4(4,J)                              
C                                                                       
            IF(JJ.GT. 3) GO TO 454                                      
            DO I=1,3                                                    
      D3(I,J,5,4)=-S5(I+ 3,K)+S4(I+ 7,L)+S4(I+10,L)                     
     *                       -S3(I+ 5,J)-S3(I+14,J)                     
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 2) GO TO 454                                      
            M= IN6(J)                                                   
            DO I=1,4                                                    
      D2(I,J,5,4)=-S4(I+13,K)+S3(I+21,L)+S3(I+25,L)                     
     *                       -S2(I+23,M)-S2(I+27,M)                     
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 1) GO TO 454                                      
            M= IN6(L)                                                   
            DO I=1,5                                                    
      D1(I  ,5,4)=-S3(I+29,K)+S2(I+36,M)+S2(I+41,M)                     
     *                       -S1(I+30,J)-S1(I+35,J)                     
            ENDDO                                                       
  454       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..4,  K = 6  and  L = 4                   
Cjms                                                                    
      J= 0                                                              
      DO JJ=1,4                                                         
         DO II=1,JJ                                                     
            J= J+1                                                      
            K= IND(J,6,4)                                               
            L= K-3-JJ                                                   
            M= L-2-JJ                                                   
C                                                                       
      D4(  J,6,4)=-R6(K,1)+S5(1,L)+S5(3,L)                              
     *                    -S4(1,M)-S4(4,M)                              
C                                                                       
            IF(JJ.GT. 3) GO TO 464                                      
            DO I=1,3                                                    
      D3(I,J,6,4)=-S5(I+ 3,K)+S4(I+ 7,L)+S4(I+10,L)                     
     *                       -S3(I+ 5,M)-S3(I+14,M)                     
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 2) GO TO 464                                      
            N= IN6(M)                                                   
            DO I=1,4                                                    
      D2(I,J,6,4)=-S4(I+13,K)+S3(I+21,L)+S3(I+25,L)                     
     *                       -S2(I+23,N)-S2(I+27,N)                     
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 1) GO TO 464                                      
            N= IN6(L)                                                   
            DO I=1,5                                                    
      D1(I  ,6,4)=-S3(I+29,K)+S2(I+36,N)+S2(I+41,N)                     
     *                       -S1(I+30,M)-S1(I+35,M)                     
            ENDDO                                                       
  464       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
C                                                                       
      CALL FIJKL4(KX,LX,F,QX,QZ,D1,D2,D3,D4)                            
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2S   *DECK MCDV20                                           
C>                                                                      
C>    @brief   DDDP case                                                
C>                                                                      
C>    @details integration of a DDDP case                               
C>                                                                      
      SUBROUTINE MCDV20(F,QX,QZ)                                        
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      DIMENSION F(6,6,6,6)                                              
C                                                                       
      COMMON /KI4   / R0(   25),R1( 3,40),R2( 6,56),R3(10,52),          
     *                R4(15,42),R5(21,24),R6(28,12),R7(36, 4),R8(   45) 
C$omp threadprivate(/KI4/)
C                                                                       
      DIMENSION                 S1(40, 3),S2(51, 6),S3(43,10),          
     *                S4(29,15),S5(14,21),S6( 5,28)                     
C                                                                       
      PARAMETER (F02=2.0D+00)                                           
      PARAMETER (F03=3.0D+00)                                           
Cjms                                                                    
C     Simplified calculation of F(I,J,K,L) for cases                    
C     where:  I = 1..6,  J = 1..6,  K = 1..KX  and  L = 1..LX           
C     using auxiliary arrays E1, E2, E3, E4 and E5                      
Cjms                                                                    
      PARAMETER (KX=6)                                                  
      PARAMETER (LX=4)                                                  
      DIMENSION       E1(  5,KX,LX),E2(4,3,KX,LX),E3(4,6,KX,LX)         
      DIMENSION       E4(2,10,KX,LX),E5( 15,KX,LX)                      
C                                                                       
      DIMENSION INS(KX,LX),INI(KX,LX),IND(15,KX,LX)                     
      DIMENSION IN6(6)                                                  
      DATA INS/ 0, 1, 3, 0, 1, 2,  0, 1, 3, 0, 1, 2,  0, 3, 5, 1, 2, 4, 
     *          1, 4, 6, 2, 3, 5/                                       
      DATA INI/ 0, 2, 2, 1, 1, 2,  0, 2, 2, 1, 1, 2,  1, 3, 3, 2, 2, 3, 
     *          1, 3, 3, 2, 2, 3/                                       
      DATA IN6/ 1, 4, 5, 2, 6, 3/                                       
C                                                                       
      DO L=1,LX                                                         
         DO K=1,KX                                                      
            IJ= 0                                                       
            IN= INS(K,L)                                                
            DO I=1,5                                                    
               IN= IN+INI(K,L)                                          
               DO J=1,I                                                 
                  IJ= IJ+1                                              
                  IN= IN+1                                              
                  IND(IJ,K,L)= IN                                       
               ENDDO                                                    
            ENDDO                                                       
         ENDDO                                                          
      ENDDO                                                             
      DO 101 J=1,40                                                     
         DO 101 I=1, 3                                                  
  101 S1(J,I)= R1(I,J)                                                  
      DO 102 J=1,51                                                     
         DO 102 I=1, 6                                                  
  102 S2(J,I)= R2(I,J)                                                  
      DO 103 J=1,43                                                     
         DO 103 I=1,10                                                  
  103 S3(J,I)= R3(I,J)                                                  
      DO 104 J=1,29                                                     
         DO 104 I=1,15                                                  
  104 S4(J,I)= R4(I,J)                                                  
      DO 105 J=1,14                                                     
         DO 105 I=1,21                                                  
  105 S5(J,I)= R5(I,J)                                                  
      DO 106 J=1, 5                                                     
         DO 106 I=1,28                                                  
  106 S6(J,I)= R6(I,J)                                                  
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 1  and  L = 1                   
Cjms                                                                    
      DO J=1,15                                                         
         E5(  J,1,1)=+R6( J,2)+R4( J,1)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,10                                                         
         E4(1,J,1,1)=+S5( 5,J)+S3( 1,J)                                 
         E4(2,J,1,1)=+S5( 6,J)+S3( 2,J)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         M= IN6(J)                                                      
         E3(1,J,1,1)=+S4(14,J)+S2( 1,M)                                 
         E3(2,J,1,1)=+S4(15,J)+S2( 2,M)                                 
         E3(3,J,1,1)=+S4(16,J)+S2( 3,M)                                 
         E3(4,J,1,1)=+S4(17,J)+S2( 4,M)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         E2(1,J,1,1)=+S3(27,J)+S1( 1,J)                                 
         E2(2,J,1,1)=+S3(28,J)+S1( 2,J)                                 
         E2(3,J,1,1)=+S3(29,J)+S1( 3,J)                                 
         E2(4,J,1,1)=+S3(30,J)+S1( 4,J)                                 
      ENDDO                                                             
C                                                                       
         J=1                                                            
      DO I=1,5                                                          
         E1(I  ,1,1)=+S2(I+36,J)+R0(I)                                  
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 2  and  L = 1                   
Cjms                                                                    
      DO J=1,15                                                         
         K= IND(J,2,1)                                                  
         E5(  J,2,1)=+R6( K,2)+R4( J,1)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,10                                                         
         K= IND(J,2,1)                                                  
         E4(1,J,2,1)=+S5( 5,K)+S3( 1,J)                                 
         E4(2,J,2,1)=+S5( 6,K)+S3( 2,J)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         K= IND(J,2,1)                                                  
         M= IN6(J)                                                      
         E3(1,J,2,1)=+S4(14,K)+S2( 1,M)                                 
         E3(2,J,2,1)=+S4(15,K)+S2( 2,M)                                 
         E3(3,J,2,1)=+S4(16,K)+S2( 3,M)                                 
         E3(4,J,2,1)=+S4(17,K)+S2( 4,M)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         K= IND(J,2,1)                                                  
         E2(1,J,2,1)=+S3(27,K)+S1( 1,J)                                 
         E2(2,J,2,1)=+S3(28,K)+S1( 2,J)                                 
         E2(3,J,2,1)=+S3(29,K)+S1( 3,J)                                 
         E2(4,J,2,1)=+S3(30,K)+S1( 4,J)                                 
      ENDDO                                                             
C                                                                       
         J=1                                                            
         K= IND(J,2,1)                                                  
         M= IN6(K)                                                      
      DO I=1,5                                                          
         E1(I  ,2,1)=+S2(I+36,M)+R0(I)                                  
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 3  and  L = 1                   
Cjms                                                                    
      J= 0                                                              
      DO JJ=1,5                                                         
         DO II=1,JJ                                                     
            J= J+1                                                      
            K= IND(J,3,1)                                               
            L= K-2-JJ                                                   
C                                                                       
      E5(  J,3,1)=+R6(K,2)-R5(L,2)*F02+S4(1,J)+S4(2,J)                  
C                                                                       
            IF(JJ.GT. 4) GO TO 131                                      
      E4(1,J,3,1)=+S5(5,K)-S4(8,L)*F02+S3(1,J)+S3(3,J)                  
      E4(2,J,3,1)=+S5(6,K)-S4(9,L)*F02+S3(2,J)+S3(4,J)                  
C                                                                       
            IF(JJ.GT. 3) GO TO 131                                      
            M= IN6(J)                                                   
            DO I=1,4                                                    
      E3(I,J,3,1)=+S4(I+13,K)-S3(I+10,L)*F02+S2(I   ,M)+S2(I+ 4,M)      
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 2) GO TO 131                                      
            M= IN6(L)                                                   
            DO I=1,4                                                    
      E2(I,J,3,1)=+S3(I+26,K)-S2(I+20,M)*F02+S1(I   ,J)+S1(I+ 4,J)      
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 1) GO TO 131                                      
            M= IN6(K)                                                   
            DO I=1,5                                                    
      E1(I  ,3,1)=+S2(I+36,M)-S1(I+20,L)*F02+R0(I   )+R0(I+ 5)          
            ENDDO                                                       
  131       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 4  and  L = 1                   
Cjms                                                                    
      DO J=1,15                                                         
         K= IND(J,4,1)                                                  
         E5(  J,4,1)=+R6( K,2)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,10                                                         
         K= IND(J,4,1)                                                  
         E4(1,J,4,1)=+S5( 5,K)                                          
         E4(2,J,4,1)=+S5( 6,K)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         K= IND(J,4,1)                                                  
         E3(1,J,4,1)=+S4(14,K)                                          
         E3(2,J,4,1)=+S4(15,K)                                          
         E3(3,J,4,1)=+S4(16,K)                                          
         E3(4,J,4,1)=+S4(17,K)                                          
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         K= IND(J,4,1)                                                  
         E2(1,J,4,1)=+S3(27,K)                                          
         E2(2,J,4,1)=+S3(28,K)                                          
         E2(3,J,4,1)=+S3(29,K)                                          
         E2(4,J,4,1)=+S3(30,K)                                          
      ENDDO                                                             
C                                                                       
         J=1                                                            
         K= IND(J,4,1)                                                  
         M= IN6(K)                                                      
      DO I=1,5                                                          
         E1(I  ,4,1)=+S2(I+36,M)                                        
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 5  and  L = 1                   
Cjms                                                                    
      DO J=1,15                                                         
         K= IND(J,5,1)                                                  
         E5(  J,5,1)=+R6( K,2)-R5( J,2)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,10                                                         
         K= IND(J,5,1)                                                  
         E4(1,J,5,1)=+S5( 5,K)-S4( 8,J)                                 
         E4(2,J,5,1)=+S5( 6,K)-S4( 9,J)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         K= IND(J,5,1)                                                  
         E3(1,J,5,1)=+S4(14,K)-S3(11,J)                                 
         E3(2,J,5,1)=+S4(15,K)-S3(12,J)                                 
         E3(3,J,5,1)=+S4(16,K)-S3(13,J)                                 
         E3(4,J,5,1)=+S4(17,K)-S3(14,J)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         K= IND(J,5,1)                                                  
         M= IN6(J)                                                      
         E2(1,J,5,1)=+S3(27,K)-S2(21,M)                                 
         E2(2,J,5,1)=+S3(28,K)-S2(22,M)                                 
         E2(3,J,5,1)=+S3(29,K)-S2(23,M)                                 
         E2(4,J,5,1)=+S3(30,K)-S2(24,M)                                 
      ENDDO                                                             
C                                                                       
         J=1                                                            
         K= IND(J,5,1)                                                  
         M= IN6(K)                                                      
      DO I=1,5                                                          
         E1(I  ,5,1)=+S2(I+36,M)-S1(I+20,J)                             
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 6  and  L = 1                   
Cjms                                                                    
      J= 0                                                              
      DO JJ=1,5                                                         
         DO II=1,JJ                                                     
            J= J+1                                                      
            K= IND(J,6,1)                                               
            L= K-2-JJ                                                   
C                                                                       
      E5(  J,6,1)=+R6(K,2)-R5(L,2)                                      
C                                                                       
            IF(JJ.GT. 4) GO TO 161                                      
      E4(1,J,6,1)=+S5(5,K)-S4(8,L)                                      
      E4(2,J,6,1)=+S5(6,K)-S4(9,L)                                      
C                                                                       
            IF(JJ.GT. 3) GO TO 161                                      
            DO I=1,4                                                    
      E3(I,J,6,1)=+S4(I+13,K)-S3(I+10,L)                                
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 2) GO TO 161                                      
            M= IN6(L)                                                   
            DO I=1,4                                                    
      E2(I,J,6,1)=+S3(I+26,K)-S2(I+20,M)                                
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 1) GO TO 161                                      
            M= IN6(K)                                                   
            DO I=1,5                                                    
      E1(I  ,6,1)=+S2(I+36,M)-S1(I+20,L)                                
            ENDDO                                                       
  161       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 1  and  L = 2                   
Cjms                                                                    
      DO J=1,15                                                         
         E5(  J,1,2)=-R7( J,1)-R5( J,4)*F03                             
      ENDDO                                                             
C                                                                       
      DO J=1,10                                                         
         E4(1,J,1,2)=-S6( 4,J)-S4(12,J)*F03                             
         E4(2,J,1,2)=-S6( 5,J)-S4(13,J)*F03                             
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         E3(1,J,1,2)=-S5(11,J)-S3(19,J)*F03                             
         E3(2,J,1,2)=-S5(12,J)-S3(20,J)*F03                             
         E3(3,J,1,2)=-S5(13,J)-S3(21,J)*F03                             
         E3(4,J,1,2)=-S5(14,J)-S3(22,J)*F03                             
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         M= IN6(J)                                                      
         E2(1,J,1,2)=-S4(26,J)-S2(29,M)*F03                             
         E2(2,J,1,2)=-S4(27,J)-S2(30,M)*F03                             
         E2(3,J,1,2)=-S4(28,J)-S2(31,M)*F03                             
         E2(4,J,1,2)=-S4(29,J)-S2(32,M)*F03                             
      ENDDO                                                             
C                                                                       
         J=1                                                            
      DO I=1,5                                                          
         E1(I  ,1,2)=-S3(I+38,J)-S1(I+30,J)*F03                         
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 2  and  L = 2                   
Cjms                                                                    
      DO J=1,15                                                         
         K= IND(J,2,2)                                                  
         E5(  J,2,2)=-R7( K,1)-R5( J,4)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,10                                                         
         K= IND(J,2,2)                                                  
         E4(1,J,2,2)=-S6( 4,K)-S4(12,J)                                 
         E4(2,J,2,2)=-S6( 5,K)-S4(13,J)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         K= IND(J,2,2)                                                  
         E3(1,J,2,2)=-S5(11,K)-S3(19,J)                                 
         E3(2,J,2,2)=-S5(12,K)-S3(20,J)                                 
         E3(3,J,2,2)=-S5(13,K)-S3(21,J)                                 
         E3(4,J,2,2)=-S5(14,K)-S3(22,J)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         K= IND(J,2,2)                                                  
         M= IN6(J)                                                      
         E2(1,J,2,2)=-S4(26,K)-S2(29,M)                                 
         E2(2,J,2,2)=-S4(27,K)-S2(30,M)                                 
         E2(3,J,2,2)=-S4(28,K)-S2(31,M)                                 
         E2(4,J,2,2)=-S4(29,K)-S2(32,M)                                 
      ENDDO                                                             
C                                                                       
         J=1                                                            
         K= IND(J,2,2)                                                  
      DO I=1,5                                                          
         E1(I  ,2,2)=-S3(I+38,K)-S1(I+30,J)                             
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 3  and  L = 2                   
Cjms                                                                    
      J= 0                                                              
      DO JJ=1,5                                                         
         DO II=1,JJ                                                     
            J= J+1                                                      
            K= IND(J,3,2)                                               
            L= K-2-JJ                                                   
C                                                                       
      E5(  J,3,2)=-R7(K,1)+R6(L,3)*F02-S5(3,J)-S5(4,J)                  
C                                                                       
            IF(JJ.GT. 4) GO TO 232                                      
      E4(1,J,3,2)=-S6(4,K)+S5(7,L)*F02-S4(10,J)-S4(12,J)                
      E4(2,J,3,2)=-S6(5,K)+S5(8,L)*F02-S4(11,J)-S4(13,J)                
C                                                                       
            IF(JJ.GT. 3) GO TO 232                                      
            DO I=1,4                                                    
      E3(I,J,3,2)=-S5(I+10,K)+S4(I+17,L)*F02-S3(I+14,J)-S3(I+18,J)      
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 2) GO TO 232                                      
            M= IN6(J)                                                   
            DO I=1,4                                                    
      E2(I,J,3,2)=-S4(I+25,K)+S3(I+30,L)*F02-S2(I+24,M)-S2(I+28,M)      
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 1) GO TO 232                                      
            M= IN6(L)                                                   
            DO I=1,5                                                    
      E1(I  ,3,2)=-S3(I+38,K)+S2(I+41,M)*F02-S1(I+25,J)-S1(I+30,J)      
            ENDDO                                                       
  232       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 4  and  L = 2                   
Cjms                                                                    
      DO J=1,15                                                         
         K= IND(J,4,2)                                                  
         E5(  J,4,2)=-R7( K,1)-R5( K,4)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,10                                                         
         K= IND(J,4,2)                                                  
         E4(1,J,4,2)=-S6( 4,K)-S4(12,K)                                 
         E4(2,J,4,2)=-S6( 5,K)-S4(13,K)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         K= IND(J,4,2)                                                  
         E3(1,J,4,2)=-S5(11,K)-S3(19,K)                                 
         E3(2,J,4,2)=-S5(12,K)-S3(20,K)                                 
         E3(3,J,4,2)=-S5(13,K)-S3(21,K)                                 
         E3(4,J,4,2)=-S5(14,K)-S3(22,K)                                 
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         K= IND(J,4,2)                                                  
         M= IN6(K)                                                      
         E2(1,J,4,2)=-S4(26,K)-S2(29,M)                                 
         E2(2,J,4,2)=-S4(27,K)-S2(30,M)                                 
         E2(3,J,4,2)=-S4(28,K)-S2(31,M)                                 
         E2(4,J,4,2)=-S4(29,K)-S2(32,M)                                 
      ENDDO                                                             
C                                                                       
         J=1                                                            
         K= IND(J,4,2)                                                  
      DO I=1,5                                                          
         E1(I  ,4,2)=-S3(I+38,K)-S1(I+30,K)                             
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 5  and  L = 2                   
Cjms                                                                    
      DO J=1,15                                                         
         K= IND(J,5,2)                                                  
         E5(  J,5,2)=-R7( K,1)+R6( J,3)-R5( K,4)+R4( J,5)               
      ENDDO                                                             
C                                                                       
      DO J=1,10                                                         
         K= IND(J,5,2)                                                  
         E4(1,J,5,2)=-S6( 4,K)+S5( 7,J)-S4(12,K)+S3( 9,J)               
         E4(2,J,5,2)=-S6( 5,K)+S5( 8,J)-S4(13,K)+S3(10,J)               
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         K= IND(J,5,2)                                                  
         M= IN6(J)                                                      
         E3(1,J,5,2)=-S5(11,K)+S4(18,J)-S3(19,K)+S2(17,M)               
         E3(2,J,5,2)=-S5(12,K)+S4(19,J)-S3(20,K)+S2(18,M)               
         E3(3,J,5,2)=-S5(13,K)+S4(20,J)-S3(21,K)+S2(19,M)               
         E3(4,J,5,2)=-S5(14,K)+S4(21,J)-S3(22,K)+S2(20,M)               
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         K= IND(J,5,2)                                                  
         M= IN6(K)                                                      
         E2(1,J,5,2)=-S4(26,K)+S3(31,J)-S2(29,M)+S1(17,J)               
         E2(2,J,5,2)=-S4(27,K)+S3(32,J)-S2(30,M)+S1(18,J)               
         E2(3,J,5,2)=-S4(28,K)+S3(33,J)-S2(31,M)+S1(19,J)               
         E2(4,J,5,2)=-S4(29,K)+S3(34,J)-S2(32,M)+S1(20,J)               
      ENDDO                                                             
C                                                                       
         J=1                                                            
         K= IND(J,5,2)                                                  
      DO I=1,5                                                          
         E1(I  ,5,2)=-S3(I+38,K)+S2(I+41,J)-S1(I+30,K)+R0(I+20)         
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 6  and  L = 2                   
Cjms                                                                    
      J= 0                                                              
      DO JJ=1,5                                                         
         DO II=1,JJ                                                     
            J= J+1                                                      
            K= IND(J,6,2)                                               
            L= K-2-JJ                                                   
C                                                                       
      E5(  J,6,2)=-R7(K,1)+R6(L,3)                                      
C                                                                       
            IF(JJ.GT. 4) GO TO 262                                      
      E4(1,J,6,2)=-S6(4,K)+S5(7,L)                                      
      E4(2,J,6,2)=-S6(5,K)+S5(8,L)                                      
C                                                                       
            IF(JJ.GT. 3) GO TO 262                                      
            DO I=1,4                                                    
      E3(I,J,6,2)=-S5(I+10,K)+S4(I+17,L)                                
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 2) GO TO 262                                      
            DO I=1,4                                                    
      E2(I,J,6,2)=-S4(I+25,K)+S3(I+30,L)                                
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 1) GO TO 262                                      
            M= IN6(L)                                                   
            DO I=1,5                                                    
      E1(I  ,6,2)=-S3(I+38,K)+S2(I+41,M)                                
            ENDDO                                                       
  262       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 1  and  L = 3                   
C                                                                       
C     next statements commented because F(I,J,1,3)                      
C     will be obtained in FIJKL5 from   F(I,J,4,2)                      
Cjms                                                                    
C     DO 315 J=1,15                                                     
C 315 E5(  J,1,3)= E5(  J,4,2)                                          
C     DO 314 J=1,10                                                     
C        DO 314 I=1,2                                                   
C 314 E4(I,J,1,3)= E4(I,J,4,2)                                          
C     DO 313 J=1,6                                                      
C        DO 313 I=1,4                                                   
C 313 E3(I,J,1,3)= E3(I,J,4,2)                                          
C     DO 312 J=1,3                                                      
C        DO 312 I=1,4                                                   
C 312 E2(I,J,1,3)= E2(I,J,4,2)                                          
C        DO 311 I=1,5                                                   
C 311 E1(I  ,1,3)= E1(I  ,4,2)                                          
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 2  and  L = 3                   
Cjms                                                                    
      J= 0                                                              
      DO JJ=1,5                                                         
         DO II=1,JJ                                                     
            J= J+1                                                      
            K= IND(J,2,3)                                               
            L= K-3-JJ-JJ                                                
C                                                                       
      E5(  J,2,3)=-R7(K,1)-R5(L,4)*F03                                  
C                                                                       
            IF(JJ.GT. 4) GO TO 323                                      
      E4(1,J,2,3)=-S6(4,K)-S4(12,L)*F03                                 
      E4(2,J,2,3)=-S6(5,K)-S4(13,L)*F03                                 
C                                                                       
            IF(JJ.GT. 3) GO TO 323                                      
            DO I=1,4                                                    
      E3(I,J,2,3)=-S5(I+10,K)-S3(I+18,L)*F03                            
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 2) GO TO 323                                      
            M= IN6(L)                                                   
            DO I=1,4                                                    
      E2(I,J,2,3)=-S4(I+25,K)-S2(I+28,M)*F03                            
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 1) GO TO 323                                      
            DO I=1,5                                                    
      E1(I  ,2,3)=-S3(I+38,K)-S1(I+30,L)*F03                            
            ENDDO                                                       
  323       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 3  and  L = 3                   
Cjms                                                                    
      J= 0                                                              
      DO JJ=1,5                                                         
         DO II=1,JJ                                                     
            J= J+1                                                      
            K= IND(J,3,3)                                               
            L= K-3-JJ                                                   
            M= L-2-JJ                                                   
C                                                                       
      E5(  J,3,3)=-R7(K,1)+R6(L,3)*F02-S5(3,M)-S5(4,M)                  
C                                                                       
            IF(JJ.GT. 4) GO TO 333                                      
      E4(1,J,3,3)=-S6(4,K)+S5(7,L)*F02-S4(10,M)-S4(12,M)                
      E4(2,J,3,3)=-S6(5,K)+S5(8,L)*F02-S4(11,M)-S4(13,M)                
C                                                                       
            IF(JJ.GT. 3) GO TO 333                                      
            DO I=1,4                                                    
      E3(I,J,3,3)=-S5(I+10,K)+S4(I+17,L)*F02-S3(I+14,M)-S3(I+18,M)      
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 2) GO TO 333                                      
            N= IN6(M)                                                   
            DO I=1,4                                                    
      E2(I,J,3,3)=-S4(I+25,K)+S3(I+30,L)*F02-S2(I+24,N)-S2(I+28,N)      
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 1) GO TO 333                                      
            N= IN6(L)                                                   
            DO I=1,5                                                    
      E1(I  ,3,3)=-S3(I+38,K)+S2(I+41,N)*F02-S1(I+25,M)-S1(I+30,M)      
            ENDDO                                                       
  333       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 4  and  L = 3                   
C                                                                       
C     next statements commented because F(I,J,4,3)                      
C     will be obtained in FIJKL5 from   F(I,J,2,2)                      
Cjms                                                                    
C     DO 345 J=1,15                                                     
C 345 E5(  J,4,3)= E5(  J,2,2)                                          
C     DO 344 J=1,10                                                     
C        DO 344 I=1,2                                                   
C 344 E4(I,J,4,3)= E4(I,J,2,2)                                          
C     DO 343 J=1,6                                                      
C        DO 343 I=1,4                                                   
C 343 E3(I,J,4,3)= E3(I,J,2,2)                                          
C     DO 342 J=1,3                                                      
C        DO 342 I=1,4                                                   
C 342 E2(I,J,4,3)= E2(I,J,2,2)                                          
C        DO 341 I=1,5                                                   
C 341 E1(I  ,4,3)= E1(I  ,2,2)                                          
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 5  and  L = 3                   
C                                                                       
C     next statements commented because F(I,J,5,3)                      
C     will be obtained in FIJKL5 from   F(I,J,6,2)                      
Cjms                                                                    
C     DO 355 J=1,15                                                     
C 355 E5(  J,5,3)= E5(  J,6,2)                                          
C     DO 354 J=1,10                                                     
C        DO 354 I=1,2                                                   
C 354 E4(I,J,5,3)= E4(I,J,6,2)                                          
C     DO 353 J=1,6                                                      
C        DO 353 I=1,4                                                   
C 353 E3(I,J,5,3)= E3(I,J,6,2)                                          
C     DO 352 J=1,3                                                      
C        DO 352 I=1,4                                                   
C 352 E2(I,J,5,3)= E2(I,J,6,2)                                          
C        DO 351 I=1,5                                                   
C 351 E1(I  ,5,3)= E1(I  ,6,2)                                          
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 6  and  L = 3                   
Cjms                                                                    
      J= 0                                                              
      DO JJ=1,5                                                         
         DO II=1,JJ                                                     
            J= J+1                                                      
            K= IND(J,6,3)                                               
            L= K-3-JJ                                                   
            M= L-JJ                                                     
C                                                                       
      E5(  J,6,3)=-R7(K,1)+R6(L,3)-R5(M,4)+R4(J,5)                      
C                                                                       
            IF(JJ.GT. 4) GO TO 363                                      
      E4(1,J,6,3)=-S6(4,K)+S5(7,L)-S4(12,M)+S3( 9,J)                    
      E4(2,J,6,3)=-S6(5,K)+S5(8,L)-S4(13,M)+S3(10,J)                    
C                                                                       
            IF(JJ.GT. 3) GO TO 363                                      
            N= IN6(J)                                                   
            DO I=1,4                                                    
      E3(I,J,6,3)=-S5(I+10,K)+S4(I+17,L)-S3(I+18,M)+S2(I+16,N)          
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 2) GO TO 363                                      
            N= IN6(M)                                                   
            DO I=1,4                                                    
      E2(I,J,6,3)=-S4(I+25,K)+S3(I+30,L)-S2(I+28,N)+S1(I+16,J)          
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 1) GO TO 363                                      
            N= IN6(L)                                                   
            DO I=1,5                                                    
      E1(I  ,6,3)=-S3(I+38,K)+S2(I+41,N)-S1(I+30,M)+R0(I+20)            
            ENDDO                                                       
  363       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 1  and  L = 4                   
Cjms                                                                    
      DO J=1,15                                                         
         K= IND(J,1,4)                                                  
         E5(  J,1,4)=-R7( K,1)+R6( J,1)-R5( K,4)+R4( J,4)               
      ENDDO                                                             
C                                                                       
      DO J=1,10                                                         
         K= IND(J,1,4)                                                  
         E4(1,J,1,4)=-S6( 4,K)+S5( 9,J)-S4(12,K)+S3( 7,J)               
         E4(2,J,1,4)=-S6( 5,K)+S5(10,J)-S4(13,K)+S3( 8,J)               
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         K= IND(J,1,4)                                                  
         M= IN6(J)                                                      
         E3(1,J,1,4)=-S5(11,K)+S4(22,J)-S3(19,K)+S2(13,M)               
         E3(2,J,1,4)=-S5(12,K)+S4(23,J)-S3(20,K)+S2(14,M)               
         E3(3,J,1,4)=-S5(13,K)+S4(24,J)-S3(21,K)+S2(15,M)               
         E3(4,J,1,4)=-S5(14,K)+S4(25,J)-S3(22,K)+S2(16,M)               
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         K= IND(J,1,4)                                                  
         M= IN6(K)                                                      
         E2(1,J,1,4)=-S4(26,K)+S3(35,J)-S2(29,M)+S1(13,J)               
         E2(2,J,1,4)=-S4(27,K)+S3(36,J)-S2(30,M)+S1(14,J)               
         E2(3,J,1,4)=-S4(28,K)+S3(37,J)-S2(31,M)+S1(15,J)               
         E2(4,J,1,4)=-S4(29,K)+S3(38,J)-S2(32,M)+S1(16,J)               
      ENDDO                                                             
C                                                                       
         J=1                                                            
         K= IND(J,1,4)                                                  
      DO I=1,5                                                          
         E1(I  ,1,4)=-S3(I+38,K)+S2(I+46,J)-S1(I+30,K)+R0(I+15)         
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 2  and  L = 4                   
Cjms                                                                    
      J= 0                                                              
      DO JJ=1,5                                                         
         DO II=1,JJ                                                     
            J= J+1                                                      
            K= IND(J,2,4)                                               
            L= K-3-JJ                                                   
            M= L-JJ                                                     
C                                                                       
      E5(  J,2,4)=-R7(K,1)+R6(L,1)-R5(M,4)+R4(J,4)                      
C                                                                       
            IF(JJ.GT. 4) GO TO 424                                      
      E4(1,J,2,4)=-S6(4,K)+S5( 9,L)-S4(12,M)+S3(7,J)                    
      E4(2,J,2,4)=-S6(5,K)+S5(10,L)-S4(13,M)+S3(8,J)                    
C                                                                       
            IF(JJ.GT. 3) GO TO 424                                      
            N= IN6(J)                                                   
            DO I=1,4                                                    
      E3(I,J,2,4)=-S5(I+10,K)+S4(I+21,L)-S3(I+18,M)+S2(I+12,N)          
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 2) GO TO 424                                      
            N= IN6(M)                                                   
            DO I=1,4                                                    
      E2(I,J,2,4)=-S4(I+25,K)+S3(I+34,L)-S2(I+28,N)+S1(I+12,J)          
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 1) GO TO 424                                      
            N= IN6(L)                                                   
            DO I=1,5                                                    
      E1(I  ,2,4)=-S3(I+38,K)+S2(I+46,N)-S1(I+30,M)+R0(I+15)            
            ENDDO                                                       
  424       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 3  and  L = 4                   
Cjms                                                                    
      J= 0                                                              
      DO JJ=1,5                                                         
         DO II=1,JJ                                                     
            J= J+1                                                      
            K= IND(J,3,4)                                               
            L= K-3-JJ                                                   
            M= L-2-JJ                                                   
C                                                                       
      E5(  J,3,4)=-R7(K,1)    +S6(1,L)+S6(3,L)*F02                      
     *            -S5(1,M)*F02-S5(3,M)-S5(4,M)*F03                      
     *            +S4(3,J)    +S4(4,J)+S4(5,J)*F02                      
C                                                                       
            IF(JJ.GT. 4) GO TO 434                                      
            DO I=1,2                                                    
      E4(I,J,3,4)=-S6(I+ 3,K)    +S5(I+ 6,L)*F02+S5(I+ 8,L)             
     *            -S4(I+ 5,M)*F02-S4(I+ 9,M)    -S4(I+11,M)*F03         
     *            +S3(I+ 4,J)    +S3(I+ 6,J)    +S3(I+ 8,J)*F02         
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 3) GO TO 434                                      
            N= IN6(J)                                                   
            DO I=1,4                                                    
      E3(I,J,3,4)=-S5(I+10,K)+S4(I+17,L)*F02+S4(I+21,L)                 
     *            -S3(I+14,M)-S3(I+18,M)*F03-S3(I+22,M)*F02             
     *            +S2(I+ 8,N)+S2(I+12,N)    +S2(I+16,N)*F02             
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 2) GO TO 434                                      
            N= IN6(M)                                                   
            DO I=1,4                                                    
      E2(I,J,3,4)=-S4(I+25,K)+S3(I+30,L)*F02+S3(I+34,L)                 
     *            -S2(I+24,N)-S2(I+28,N)*F03-S2(I+32,N)*F02             
     *            +S1(I+ 8,J)+S1(I+12,J)    +S1(I+16,J)*F02             
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 1) GO TO 434                                      
            N= IN6(L)                                                   
            DO I=1,5                                                    
      E1(I  ,3,4)=-S3(I+38,K)+S2(I+41,N)*F02+S2(I+46,N)                 
     *            -S1(I+25,M)-S1(I+30,M)*F03-S1(I+35,M)*F02             
     *            +R0(I+10  )+R0(I+15  )    +R0(I+20  )*F02             
            ENDDO                                                       
  434       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 4  and  L = 4                   
Cjms                                                                    
      J= 0                                                              
      DO JJ=1,5                                                         
         DO II=1,JJ                                                     
            J= J+1                                                      
            K= IND(J,4,4)                                               
            L= K-2-JJ                                                   
C                                                                       
      E5(  J,4,4)=-R7(K,1)+R6(L,1)                                      
C                                                                       
            IF(JJ.GT. 4) GO TO 444                                      
      E4(1,J,4,4)=-S6(4,K)+S5( 9,L)                                     
      E4(2,J,4,4)=-S6(5,K)+S5(10,L)                                     
C                                                                       
            IF(JJ.GT. 3) GO TO 444                                      
            DO I=1,4                                                    
      E3(I,J,4,4)=-S5(I+10,K)+S4(I+21,L)                                
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 2) GO TO 444                                      
            DO I=1,4                                                    
      E2(I,J,4,4)=-S4(I+25,K)+S3(I+34,L)                                
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 1) GO TO 444                                      
            M= IN6(L)                                                   
            DO I=1,5                                                    
      E1(I  ,4,4)=-S3(I+38,K)+S2(I+46,M)                                
            ENDDO                                                       
  444       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 5  and  L = 4                   
Cjms                                                                    
      J= 0                                                              
      DO JJ=1,5                                                         
         DO II=1,JJ                                                     
            J= J+1                                                      
            K= IND(J,5,4)                                               
            L= K-2-JJ                                                   
C                                                                       
      E5(  J,5,4)=-R7(K,1)+S6(1,L)+S6(3,L)                              
     *                    -S5(1,J)-S5(4,J)                              
C                                                                       
            IF(JJ.GT. 4) GO TO 454                                      
      E4(1,J,5,4)=-S6(4,K)+S5(7,L)+S5( 9,L)                             
     *                    -S4(6,J)-S4(12,J)                             
      E4(2,J,5,4)=-S6(5,K)+S5(8,L)+S5(10,L)                             
     *                    -S4(7,J)-S4(13,J)                             
C                                                                       
            IF(JJ.GT. 3) GO TO 454                                      
            DO I=1,4                                                    
      E3(I,J,5,4)=-S5(I+10,K)+S4(I+17,L)+S4(I+21,L)                     
     *                       -S3(I+18,J)-S3(I+22,J)                     
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 2) GO TO 454                                      
            M= IN6(J)                                                   
            DO I=1,4                                                    
      E2(I,J,5,4)=-S4(I+25,K)+S3(I+30,L)+S3(I+34,L)                     
     *                       -S2(I+28,M)-S2(I+32,M)                     
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 1) GO TO 454                                      
            M= IN6(L)                                                   
            DO I=1,5                                                    
      E1(I  ,5,4)=-S3(I+38,K)+S2(I+41,M)+S2(I+46,M)                     
     *                       -S1(I+30,J)-S1(I+35,J)                     
            ENDDO                                                       
  454       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 6  and  L = 4                   
Cjms                                                                    
      J= 0                                                              
      DO JJ=1,5                                                         
         DO II=1,JJ                                                     
            J= J+1                                                      
            K= IND(J,6,4)                                               
            L= K-3-JJ                                                   
            M= L-2-JJ                                                   
C                                                                       
      E5(  J,6,4)=-R7(K,1)+S6(1,L)+S6(3,L)                              
     *                    -S5(1,M)-S5(4,M)                              
C                                                                       
            IF(JJ.GT. 4) GO TO 464                                      
      E4(1,J,6,4)=-S6(4,K)+S5(7,L)+S5( 9,L)                             
     *                    -S4(6,M)-S4(12,M)                             
      E4(2,J,6,4)=-S6(5,K)+S5(8,L)+S5(10,L)                             
     *                    -S4(7,M)-S4(13,M)                             
C                                                                       
            IF(JJ.GT. 3) GO TO 464                                      
            DO I=1,4                                                    
      E3(I,J,6,4)=-S5(I+10,K)+S4(I+17,L)+S4(I+21,L)                     
     *                       -S3(I+18,M)-S3(I+22,M)                     
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 2) GO TO 464                                      
            N= IN6(M)                                                   
            DO I=1,4                                                    
      E2(I,J,6,4)=-S4(I+25,K)+S3(I+30,L)+S3(I+34,L)                     
     *                       -S2(I+28,N)-S2(I+32,N)                     
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 1) GO TO 464                                      
            N= IN6(L)                                                   
            DO I=1,5                                                    
      E1(I  ,6,4)=-S3(I+38,K)+S2(I+41,N)+S2(I+46,N)                     
     *                       -S1(I+30,M)-S1(I+35,M)                     
            ENDDO                                                       
  464       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
C                                                                       
      CALL FIJKL5(KX,LX,F,QX,QZ,E1,E2,E3,E4,E5)                         
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2S   *DECK MCDV21                                           
C>                                                                      
C>    @brief   DDDD case                                                
C>                                                                      
C>    @details integration of a DDDD case                               
C>                                                                      
      SUBROUTINE MCDV21(F,QX,QZ)                                        
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      DIMENSION F(6,6,6,6)                                              
C                                                                       
      COMMON /KI4   / R0(   25),R1( 3,40),R2( 6,56),R3(10,52),          
     *                R4(15,42),R5(21,24),R6(28,12),R7(36, 4),R8(   45) 
C$omp threadprivate(/KI4/)
C                                                                       
      DIMENSION                 S1(40, 3),S2(56, 6),S3(52,10),          
     *                S4(42,15),S5(24,21),S6(12,28),S7( 4,36)           
C                                                                       
      PARAMETER (F02=2.0D+00)                                           
      PARAMETER (F03=3.0D+00)                                           
      PARAMETER (F04=4.0D+00)                                           
      PARAMETER (F06=6.0D+00)                                           
C                                                                       
      DIMENSION       T531E1(  5),T531E2(4,3),T531E3(4,6),T531E4(2,10), 
     *                T531E5( 15),                                      
     *                V531E1(  5),V531E2(4,3),V531E3(4,6),V531E4(2,10), 
     *                V531E5( 15),                                      
     *                U5E1(  5,4),U5E2(4,3,4),U5E3(4,6,4),U5E4(2,10,4), 
     *                U5E5( 15,4),                                      
     *                T632E1(  5),T632E2(4,3),T632E3(4,6),T632E4(2,10), 
     *                T632E5( 15),                                      
     *                V632E1(  5),V632E2(4,3),V632E3(4,6),V632E4(2,10), 
     *                V632E5( 15),                                      
     *                U6E1(  5,4),U6E2(4,3,4),U6E3(4,6,4),U6E4(2,10,4), 
     *                U6E5( 15,4)                                       
Cjms                                                                    
C     Simplified calculation of F(I,J,K,L) for cases                    
C     where:  I = 1..6,  J = 1..6,  K = 1..KX  and  L = 1..LX           
C     using auxiliary arrays E1, E2, E3, E4 and E5                      
Cjms                                                                    
      PARAMETER (KX=6)                                                  
      PARAMETER (LX=6)                                                  
      DIMENSION       E1(  5,KX,LX),E2(4,3,KX,LX),E3(4,6,KX,LX)         
      DIMENSION       E4(2,10,KX,LX),E5( 15,KX,LX)                      
C                                                                       
      DIMENSION INS(KX,LX),INI(KX,LX),IND(15,KX,LX)                     
      DIMENSION IN6(6),JND(15)                                          
      DATA INS/ 0, 1, 3, 0, 1, 2,  1, 6, 8, 3, 4, 7,  3, 8,10, 5, 6, 9, 
     *          0, 3, 5, 1, 2, 4,  1, 4, 6, 2, 3, 5,  2, 7, 9, 4, 5, 8/ 
      DATA INI/ 0, 2, 2, 1, 1, 2,  2, 4, 4, 3, 3, 4,  2, 4, 4, 3, 3, 4, 
     *          1, 3, 3, 2, 2, 3,  1, 3, 3, 2, 2, 3,  2, 4, 4, 3, 3, 4/ 
      DATA IN6/ 1, 4, 5, 2, 6, 3/                                       
C                                                                       
      DO L=1,LX                                                         
         DO K=1,KX                                                      
            IJ= 0                                                       
            IN= INS(K,L)                                                
            DO I=1,5                                                    
               IN= IN+INI(K,L)                                          
               DO J=1,I                                                 
                  IJ= IJ+1                                              
                  IN= IN+1                                              
                  IND(IJ,K,L)= IN                                       
               ENDDO                                                    
            ENDDO                                                       
         ENDDO                                                          
      ENDDO                                                             
      DO 101 J=1,40                                                     
         DO 101 I=1, 3                                                  
  101 S1(J,I)= R1(I,J)                                                  
      DO 102 J=1,56                                                     
         DO 102 I=1, 6                                                  
  102 S2(J,I)= R2(I,J)                                                  
      DO 103 J=1,52                                                     
         DO 103 I=1,10                                                  
  103 S3(J,I)= R3(I,J)                                                  
      DO 104 J=1,42                                                     
         DO 104 I=1,15                                                  
  104 S4(J,I)= R4(I,J)                                                  
      DO 105 J=1,24                                                     
         DO 105 I=1,21                                                  
  105 S5(J,I)= R5(I,J)                                                  
      DO 106 J=1,12                                                     
         DO 106 I=1,28                                                  
  106 S6(J,I)= R6(I,J)                                                  
      DO 107 J=1, 4                                                     
         DO 107 I=1,36                                                  
  107 S7(J,I)= R7(I,J)                                                  
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 1  and  L = 1                   
Cjms                                                                    
      DO J=1,15                                                         
         E5(  J,1,1)=+R8(   J)+R6( J,1)*F06+R4( J,1)*F03                
      ENDDO                                                             
C                                                                       
      DO J=1,10                                                         
         E4(1,J,1,1)=+S7( 3,J)+S5( 5,J)*F06+S3( 1,J)*F03                
         E4(2,J,1,1)=+S7( 4,J)+S5( 6,J)*F06+S3( 2,J)*F03                
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         M= IN6(J)                                                      
         E3(1,J,1,1)=+S6( 9,J)+S4(14,J)*F06+S2( 1,M)*F03                
         E3(2,J,1,1)=+S6(10,J)+S4(15,J)*F06+S2( 2,M)*F03                
         E3(3,J,1,1)=+S6(11,J)+S4(16,J)*F06+S2( 3,M)*F03                
         E3(4,J,1,1)=+S6(12,J)+S4(17,J)*F06+S2( 4,M)*F03                
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         E2(1,J,1,1)=+S5(21,J)+S3(27,J)*F06+S1( 1,J)*F03                
         E2(2,J,1,1)=+S5(22,J)+S3(28,J)*F06+S1( 2,J)*F03                
         E2(3,J,1,1)=+S5(23,J)+S3(29,J)*F06+S1( 3,J)*F03                
         E2(4,J,1,1)=+S5(24,J)+S3(30,J)*F06+S1( 4,J)*F03                
      ENDDO                                                             
C                                                                       
         J=1                                                            
      DO I=1,5                                                          
         E1(I  ,1,1)=+S4(I+37,J)+S2(I+36,J)*F06+R0(I)*F03               
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 2  and  L = 1                   
Cjms                                                                    
      DO J=1,15                                                         
         K= IND(J,2,1)                                                  
         E5(  J,2,1)=+R8(   K)+R6( K,1)+R6( J,1)+R4( J,1)               
      ENDDO                                                             
C                                                                       
      DO J=1,10                                                         
         K= IND(J,2,1)                                                  
         E4(1,J,2,1)=+S7( 3,K)+S5( 5,K)+S5( 5,J)+S3( 1,J)               
         E4(2,J,2,1)=+S7( 4,K)+S5( 6,K)+S5( 6,J)+S3( 2,J)               
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         K= IND(J,2,1)                                                  
         M= IN6(J)                                                      
         E3(1,J,2,1)=+S6( 9,K)+S4(14,K)+S4(14,J)+S2( 1,M)               
         E3(2,J,2,1)=+S6(10,K)+S4(15,K)+S4(15,J)+S2( 2,M)               
         E3(3,J,2,1)=+S6(11,K)+S4(16,K)+S4(16,J)+S2( 3,M)               
         E3(4,J,2,1)=+S6(12,K)+S4(17,K)+S4(17,J)+S2( 4,M)               
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         K= IND(J,2,1)                                                  
         E2(1,J,2,1)=+S5(21,K)+S3(27,K)+S3(27,J)+S1( 1,J)               
         E2(2,J,2,1)=+S5(22,K)+S3(28,K)+S3(28,J)+S1( 2,J)               
         E2(3,J,2,1)=+S5(23,K)+S3(29,K)+S3(29,J)+S1( 3,J)               
         E2(4,J,2,1)=+S5(24,K)+S3(30,K)+S3(30,J)+S1( 4,J)               
      ENDDO                                                             
C                                                                       
         J=1                                                            
         K= IND(J,2,1)                                                  
         M= IN6(K)                                                      
      DO I=1,5                                                          
         E1(I  ,2,1)=+S4(I+37,K)+S2(I+36,M)+S2(I+36,J)+R0(I)            
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 3  and  L = 1                   
C                                                                       
C     Auxiliary arrays T5.., U5.. and V5.. defined here are used for    
C     E1..E5(I,J,3,1), E1..E5(I,J,1,3) and E1..E5(I,J,5,5)              
Cjms                                                                    
      J= 0                                                              
      DO JJ=1,5                                                         
         DO II=1,JJ                                                     
            J= J+1                                                      
            K= IND(J,3,1)                                               
            L= K-2-JJ                                                   
C                                                                       
      T531E5(  J)=+R8(K  )-S7(1,L)-S7(2,L)+R6(J,1)                      
     *            +R6(K,1)-S5(1,L)-S5(2,L)+R4(J,1)                      
C                                                                       
            IF(JJ.GT. 4) GO TO 131                                      
      T531E4(1,J)=+S7(3,K)-S6(5,L)-S6(7,L)+S5(5,J)                      
     *            +S5(5,K)-S4(6,L)-S4(8,L)+S3(1,J)                      
      T531E4(2,J)=+S7(4,K)-S6(6,L)-S6(8,L)+S5(6,J)                      
     *            +S5(6,K)-S4(7,L)-S4(9,L)+S3(2,J)                      
C                                                                       
            IF(JJ.GT. 3) GO TO 131                                      
            M= IN6(J)                                                   
            DO I=1,4                                                    
      T531E3(I,J)=+S6(I+ 8,K)-S5(I+12,L)-S5(I+16,L)+S4(I+13,J)          
     *            +S4(I+13,K)-S3(I+10,L)-S3(I+14,L)+S2(I,M)             
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 2) GO TO 131                                      
            M= IN6(L)                                                   
            DO I=1,4                                                    
      T531E2(I,J)=+S5(I+20,K)-S4(I+29,L)-S4(I+33,L)+S3(I+26,J)          
     *            +S3(I+26,K)-S2(I+20,M)-S2(I+24,M)+S1(I,J)             
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 1) GO TO 131                                      
            DO I=1,5                                                    
      T531E1(I  )=+S4(I+37,K)-S3(I+42,L)-S3(I+47,L)+S2(I+36,J)          
     *            +S2(I+36,L)-S1(I+20,L)-S1(I+25,L)+R0(I  )             
            ENDDO                                                       
  131       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
C                                                                       
      IJ= 0                                                             
      IN= 1                                                             
      DO I=1,5                                                          
         IN= IN+1                                                       
         DO J=1,I                                                       
            IJ= IJ+1                                                    
            IN= IN+1                                                    
            JND(IJ)= IN                                                 
         ENDDO                                                          
      ENDDO                                                             
C                                                                       
      DO J1=2,4                                                         
         DO J=1,15                                                      
               U5E5(  J,J1)= R6(  J,J1)+R4(  J,J1)                      
         ENDDO                                                          
         M1=J1+J1+ 2                                                    
         DO J=1,10                                                      
               U5E4(1,J,J1)= S5(1+M1,J)+S3(1+M1- 4,J)                   
               U5E4(2,J,J1)= S5(2+M1,J)+S3(2+M1- 4,J)                   
         ENDDO                                                          
         M1=M1+M1+ 5                                                    
         DO J=1,6                                                       
            M= IN6(J)                                                   
            DO I=1,4                                                    
               U5E3(I,J,J1)=+S4(I+M1,J)+S2(I+M1-13,M)                   
            ENDDO                                                       
         ENDDO                                                          
         M1=M1+13                                                       
         DO J=1,3                                                       
            DO I=1,4                                                    
               U5E2(I,J,J1)= S3(I+M1,J)+S1(I+M1-26,J)                   
            ENDDO                                                       
         ENDDO                                                          
         M1=M1+J1+ 9                                                    
            DO I=1,5                                                    
               U5E1(I  ,J1)= S2(I+M1,1)+R0(I+M1-36)                     
            ENDDO                                                       
      ENDDO                                                             
C                                                                       
      J= 0                                                              
      DO JJ=1,5                                                         
         DO II=1,JJ                                                     
            J= J+1                                                      
            K= JND(J)                                                   
C                                                                       
      V531E5(  J)=+S7(1,K)-S7(2,K)+S5(1,K)-S5(2,K)                      
C                                                                       
            IF(JJ.GT. 4) GO TO 132                                      
      V531E4(1,J)=+S6(5,K)-S6(7,K)+S4(6,K)-S4(8,K)                      
      V531E4(2,J)=+S6(6,K)-S6(8,K)+S4(7,K)-S4(9,K)                      
C                                                                       
            IF(JJ.GT. 3) GO TO 132                                      
            DO I=1,4                                                    
      V531E3(I,J)=+S5(I+12,K)-S5(I+16,K)+S3(I+10,K)-S3(I+14,K)          
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 2) GO TO 132                                      
            M= IN6(K)                                                   
            DO I=1,4                                                    
      V531E2(I,J)=+S4(I+29,K)-S4(I+33,K)+S2(I+20,M)-S2(I+24,M)          
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 1) GO TO 132                                      
            DO I=1,5                                                    
      V531E1(I  )=+S3(I+42,K)-S3(I+47,K)+S1(I+20,K)-S1(I+25,K)          
            ENDDO                                                       
  132       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
C                                                                       
      J1=2                                                              
      DO J=1,15                                                         
         E5(  J,3,1)= T531E5(  J)+U5E5(  J,J1)-V531E5(  J)              
      ENDDO                                                             
      DO J=1,10                                                         
         E4(1,J,3,1)= T531E4(1,J)+U5E4(1,J,J1)-V531E4(1,J)              
         E4(2,J,3,1)= T531E4(2,J)+U5E4(2,J,J1)-V531E4(2,J)              
      ENDDO                                                             
      DO J=1,6                                                          
         E3(1,J,3,1)= T531E3(1,J)+U5E3(1,J,J1)-V531E3(1,J)              
         E3(2,J,3,1)= T531E3(2,J)+U5E3(2,J,J1)-V531E3(2,J)              
         E3(3,J,3,1)= T531E3(3,J)+U5E3(3,J,J1)-V531E3(3,J)              
         E3(4,J,3,1)= T531E3(4,J)+U5E3(4,J,J1)-V531E3(4,J)              
      ENDDO                                                             
      DO J=1,3                                                          
         E2(1,J,3,1)= T531E2(1,J)+U5E2(1,J,J1)-V531E2(1,J)              
         E2(2,J,3,1)= T531E2(2,J)+U5E2(2,J,J1)-V531E2(2,J)              
         E2(3,J,3,1)= T531E2(3,J)+U5E2(3,J,J1)-V531E2(3,J)              
         E2(4,J,3,1)= T531E2(4,J)+U5E2(4,J,J1)-V531E2(4,J)              
      ENDDO                                                             
      DO I=1,5                                                          
         E1(I  ,3,1)= T531E1(I  )+U5E1(I  ,J1)-V531E1(I  )              
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 4  and  L = 1                   
Cjms                                                                    
      DO J=1,15                                                         
         K= IND(J,4,1)                                                  
         E5(  J,4,1)=+R8(   K)+R6( K,1)*F03                             
      ENDDO                                                             
C                                                                       
      DO J=1,10                                                         
         K= IND(J,4,1)                                                  
         E4(1,J,4,1)=+S7( 3,K)+S5( 5,K)*F03                             
         E4(2,J,4,1)=+S7( 4,K)+S5( 6,K)*F03                             
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         K= IND(J,4,1)                                                  
         E3(1,J,4,1)=+S6( 9,K)+S4(14,K)*F03                             
         E3(2,J,4,1)=+S6(10,K)+S4(15,K)*F03                             
         E3(3,J,4,1)=+S6(11,K)+S4(16,K)*F03                             
         E3(4,J,4,1)=+S6(12,K)+S4(17,K)*F03                             
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         K= IND(J,4,1)                                                  
         E2(1,J,4,1)=+S5(21,K)+S3(27,K)*F03                             
         E2(2,J,4,1)=+S5(22,K)+S3(28,K)*F03                             
         E2(3,J,4,1)=+S5(23,K)+S3(29,K)*F03                             
         E2(4,J,4,1)=+S5(24,K)+S3(30,K)*F03                             
      ENDDO                                                             
C                                                                       
         J=1                                                            
         K= IND(J,4,1)                                                  
         M= IN6(K)                                                      
      DO I=1,5                                                          
         E1(I  ,4,1)=+S4(I+37,K)+S2(I+36,M)*F03                         
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 5  and  L = 1                   
Cjms                                                                    
      DO J=1,15                                                         
         K= IND(J,5,1)                                                  
         E5(  J,5,1)=+R8(   K)-R7( J,1)+(R6( K,1)-R5( J,1))*F03         
      ENDDO                                                             
C                                                                       
      DO J=1,10                                                         
         K= IND(J,5,1)                                                  
         E4(1,J,5,1)=+S7( 3,K)-S6( 5,J)+(S5( 5,K)-S4( 6,J))*F03         
         E4(2,J,5,1)=+S7( 4,K)-S6( 6,J)+(S5( 6,K)-S4( 7,J))*F03         
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         K= IND(J,5,1)                                                  
         E3(1,J,5,1)=+S6( 9,K)-S5(13,J)+(S4(14,K)-S3(11,J))*F03         
         E3(2,J,5,1)=+S6(10,K)-S5(14,J)+(S4(15,K)-S3(12,J))*F03         
         E3(3,J,5,1)=+S6(11,K)-S5(15,J)+(S4(16,K)-S3(13,J))*F03         
         E3(4,J,5,1)=+S6(12,K)-S5(16,J)+(S4(17,K)-S3(14,J))*F03         
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         K= IND(J,5,1)                                                  
         M= IN6(J)                                                      
         E2(1,J,5,1)=+S5(21,K)-S4(30,J)+(S3(27,K)-S2(21,M))*F03         
         E2(2,J,5,1)=+S5(22,K)-S4(31,J)+(S3(28,K)-S2(22,M))*F03         
         E2(3,J,5,1)=+S5(23,K)-S4(32,J)+(S3(29,K)-S2(23,M))*F03         
         E2(4,J,5,1)=+S5(24,K)-S4(33,J)+(S3(30,K)-S2(24,M))*F03         
      ENDDO                                                             
C                                                                       
         J=1                                                            
         K= IND(J,5,1)                                                  
         M= IN6(K)                                                      
      DO I=1,5                                                          
         E1(I  ,5,1)=+S4(I+37,K)-S3(I+42,J)+(S2(I+36,M)-S1(I+20,J))*F03 
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 6  and  L = 1                   
Cjms                                                                    
      J= 0                                                              
      DO JJ=1,5                                                         
         DO II=1,JJ                                                     
            J= J+1                                                      
            K= IND(J,6,1)                                               
            L= K-2-JJ                                                   
C                                                                       
      E5(  J,6,1)=+R8(K  )-R7(L,1)+R6(K,1)-R5(L,1)                      
C                                                                       
            IF(JJ.GT. 4) GO TO 161                                      
      E4(1,J,6,1)=+S7(3,K)-S6(5,L)+S5(5,K)-S4(6,L)                      
      E4(2,J,6,1)=+S7(4,K)-S6(6,L)+S5(6,K)-S4(7,L)                      
C                                                                       
            IF(JJ.GT. 3) GO TO 161                                      
            DO I=1,4                                                    
      E3(I,J,6,1)=+S6(I+ 8,K)-S5(I+12,L)+S4(I+13,K)-S3(I+10,L)          
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 2) GO TO 161                                      
            M= IN6(L)                                                   
            DO I=1,4                                                    
      E2(I,J,6,1)=+S5(I+20,K)-S4(I+29,L)+S3(I+26,K)-S2(I+20,M)          
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 1) GO TO 161                                      
            M= IN6(K)                                                   
            DO I=1,5                                                    
      E1(I  ,6,1)=+S4(I+37,K)-S3(I+42,L)+S2(I+36,M)-S1(I+20,L)          
            ENDDO                                                       
  161       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 1  and  L = 2                   
C                                                                       
C     next statements commented because F(I,J,1,2)                      
C     will be obtained in FIJKL5 from   F(I,J,2,1)                      
Cjms                                                                    
C     DO 215 J=1,15                                                     
C 215 E5(  J,1,2)= E5(  J,2,1)                                          
C     DO 214 J=1,10                                                     
C        DO 214 I=1,2                                                   
C 214 E4(I,J,1,2)= E4(I,J,2,1)                                          
C     DO 213 J=1,6                                                      
C        DO 213 I=1,4                                                   
C 213 E3(I,J,1,2)= E3(I,J,2,1)                                          
C     DO 212 J=1,3                                                      
C        DO 212 I=1,4                                                   
C 212 E2(I,J,1,2)= E2(I,J,2,1)                                          
C        DO 211 I=1,5                                                   
C 211 E1(I  ,1,2)= E1(I  ,2,1)                                          
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 2  and  L = 2                   
Cjms                                                                    
      J= 0                                                              
      DO JJ=1,5                                                         
         DO II=1,JJ                                                     
            J= J+1                                                      
            K= IND(J,2,2)                                               
            L= K-5-JJ-JJ                                                
C                                                                       
      E5(  J,2,2)=+R8(K  )+R6(L,1)*F06+R4(J,1)*F03                      
C                                                                       
            IF(JJ.GT. 4) GO TO 222                                      
      E4(1,J,2,2)=+S7(3,K)+S5(5,L)*F06+S3(1,J)*F03                      
      E4(2,J,2,2)=+S7(4,K)+S5(6,L)*F06+S3(2,J)*F03                      
C                                                                       
            IF(JJ.GT. 3) GO TO 222                                      
            M= IN6(J)                                                   
            DO I=1,4                                                    
      E3(I,J,2,2)=+S6(I+ 8,K)+S4(I+13,L)*F06+S2(I,M)*F03                
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 2) GO TO 222                                      
            DO I=1,4                                                    
      E2(I,J,2,2)=+S5(I+20,K)+S3(I+26,L)*F06+S1(I,J)*F03                
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 1) GO TO 222                                      
            M= IN6(L)                                                   
            DO I=1,5                                                    
      E1(I  ,2,2)=+S4(I+37,K)+S2(I+36,M)*F06+R0(I  )*F03                
            ENDDO                                                       
  222       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 3  and  L = 2                   
C                                                                       
C     Auxiliary arrays T6.., U6.. and V6.. defined here are used for    
C     E1..E5(I,J,3,2), E1..E5(I,J,2,3) and E1..E5(I,J,6,6)              
Cjms                                                                    
      J= 0                                                              
      DO JJ=1,5                                                         
         DO II=1,JJ                                                     
            J= J+1                                                      
            K= IND(J,3,2)                                               
            L= K-4-JJ                                                   
            M= L-3-JJ                                                   
            N= M-JJ                                                     
C                                                                       
      T632E5(  J)=+R8(K  )-S7(1,L)-S7(2,L)+R6(M,1)                      
     *          +R6(M+2,1)-S5(1,N)-S5(2,N)+R4(J,1)                      
C                                                                       
            IF(JJ.GT. 4) GO TO 232                                      
      T632E4(1,J)=+S7(3,K)-S6(5,L)-S6(7,L)+S5(5,M)                      
     *          +S5(5,M+2)-S4(6,N)-S4(8,N)+S3(1,J)                      
      T632E4(2,J)=+S7(4,K)-S6(6,L)-S6(8,L)+S5(6,M)                      
     *          +S5(6,M+2)-S4(7,N)-S4(9,N)+S3(2,J)                      
C                                                                       
            IF(JJ.GT. 3) GO TO 232                                      
            JN= IN6(J)                                                  
            DO I=1,4                                                    
      T632E3(I,J)=+S6(I+ 8,K)-S5(I+12,L)-S5(I+16,L)+S4(I+13,M)          
     *          +S4(I+13,M+2)-S3(I+10,N)-S3(I+14,N)+S2(I,JN)            
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 2) GO TO 232                                      
            N= IN6(N)                                                   
            DO I=1,4                                                    
      T632E2(I,J)=+S5(I+20,K)-S4(I+29,L)-S4(I+33,L)+S3(I+26,M)          
     *          +S3(I+26,M+2)-S2(I+20,N)-S2(I+24,N)+S1(I,J)             
            ENDDO                                                       
            N= M-JJ                                                     
C                                                                       
            IF(JJ.GT. 1) GO TO 232                                      
            M= IN6(M)                                                   
            DO I=1,5                                                    
      T632E1(I  )=+S4(I+37,K)-S3(I+42,L)-S3(I+47,L)+S2(I+36,M)          
     *          +S2(I+36,M+1)-S1(I+20,N)-S1(I+25,N)+R0(I  )             
            ENDDO                                                       
  232       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
C                                                                       
      IJ= 0                                                             
      IN= 1                                                             
      DO I=1,5                                                          
         IN= IN+2                                                       
         DO J=1,I                                                       
            IJ= IJ+1                                                    
            IN= IN+1                                                    
            JND(IJ)= IN                                                 
         ENDDO                                                          
      ENDDO                                                             
C                                                                       
      DO J1=2,4                                                         
         DO J=1,15                                                      
            K= JND(J)                                                   
               U6E5(  J,J1)=+R6(K  ,J1)+R4(  J,J1)                      
         ENDDO                                                          
         M1=J1+J1+ 2                                                    
         DO J=1,10                                                      
            K= JND(J)                                                   
               U6E4(1,J,J1)=+S5(1+M1,K)+S3(1+M1- 4,J)                   
               U6E4(2,J,J1)=+S5(2+M1,K)+S3(2+M1- 4,J)                   
         ENDDO                                                          
         M1=M1+M1+ 5                                                    
         DO J=1,6                                                       
            K= JND(J)                                                   
            M= IN6(J)                                                   
            DO I=1,4                                                    
               U6E3(I,J,J1)=+S4(I+M1,K)+S2(I+M1-13,M)                   
            ENDDO                                                       
         ENDDO                                                          
         M1=M1+13                                                       
         DO J=1,3                                                       
            K= JND(J)                                                   
            DO I=1,4                                                    
               U6E2(I,J,J1)=+S3(I+M1,K)+S1(I+M1-26,J)                   
            ENDDO                                                       
         ENDDO                                                          
         M1=M1+J1+ 9                                                    
            DO I=1,5                                                    
               U6E1(I  ,J1)=+S2(I+M1,2)+R0(I+M1-36)                     
            ENDDO                                                       
      ENDDO                                                             
C                                                                       
      IJ= 0                                                             
      IN= 4                                                             
      DO I=1,5                                                          
         IN= IN+3                                                       
         DO J=1,I                                                       
            IJ= IJ+1                                                    
            IN= IN+1                                                    
            JND(IJ)= IN                                                 
         ENDDO                                                          
      ENDDO                                                             
C                                                                       
      J= 0                                                              
      DO JJ=1,5                                                         
         DO II=1,JJ                                                     
            J= J+1                                                      
            K= JND(J)                                                   
            L= K-3-JJ-JJ                                                
C                                                                       
      V632E5(  J)=+S7(1,K)-S7(2,K)+S5(1,L)-S5(2,L)                      
C                                                                       
            IF(JJ.GT. 4) GO TO 233                                      
      V632E4(1,J)=+S6(5,K)-S6(7,K)+S4(6,L)-S4(8,L)                      
      V632E4(2,J)=+S6(6,K)-S6(8,K)+S4(7,L)-S4(9,L)                      
C                                                                       
            IF(JJ.GT. 3) GO TO 233                                      
            DO I=1,4                                                    
      V632E3(I,J)=+S5(I+12,K)-S5(I+16,K)+S3(I+10,L)-S3(I+14,L)          
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 2) GO TO 233                                      
            M= IN6(L)                                                   
            DO I=1,4                                                    
      V632E2(I,J)=+S4(I+29,K)-S4(I+33,K)+S2(I+20,M)-S2(I+24,M)          
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 1) GO TO 233                                      
            DO I=1,5                                                    
      V632E1(I  )=+S3(I+42,K)-S3(I+47,K)+S1(I+20,L)-S1(I+25,L)          
            ENDDO                                                       
  233       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
C                                                                       
      J1=2                                                              
      DO J=1,15                                                         
         E5(  J,3,2)= T632E5(  J)+U6E5(  J,J1)-V632E5(  J)              
      ENDDO                                                             
      DO J=1,10                                                         
         E4(1,J,3,2)= T632E4(1,J)+U6E4(1,J,J1)-V632E4(1,J)              
         E4(2,J,3,2)= T632E4(2,J)+U6E4(2,J,J1)-V632E4(2,J)              
      ENDDO                                                             
      DO J=1,6                                                          
         E3(1,J,3,2)= T632E3(1,J)+U6E3(1,J,J1)-V632E3(1,J)              
         E3(2,J,3,2)= T632E3(2,J)+U6E3(2,J,J1)-V632E3(2,J)              
         E3(3,J,3,2)= T632E3(3,J)+U6E3(3,J,J1)-V632E3(3,J)              
         E3(4,J,3,2)= T632E3(4,J)+U6E3(4,J,J1)-V632E3(4,J)              
      ENDDO                                                             
      DO J=1,3                                                          
         E2(1,J,3,2)= T632E2(1,J)+U6E2(1,J,J1)-V632E2(1,J)              
         E2(2,J,3,2)= T632E2(2,J)+U6E2(2,J,J1)-V632E2(2,J)              
         E2(3,J,3,2)= T632E2(3,J)+U6E2(3,J,J1)-V632E2(3,J)              
         E2(4,J,3,2)= T632E2(4,J)+U6E2(4,J,J1)-V632E2(4,J)              
      ENDDO                                                             
      DO I=1,5                                                          
         E1(I  ,3,2)= T632E1(I  )+U6E1(I  ,J1)-V632E1(I  )              
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 4  and  L = 2                   
Cjms                                                                    
      J= 0                                                              
      DO JJ=1,5                                                         
         DO II=1,JJ                                                     
            J= J+1                                                      
            K= IND(J,4,2)                                               
            L= K-3-JJ-JJ                                                
C                                                                       
      E5(  J,4,2)=+R8(K  )+R6(L,1)*F03                                  
C                                                                       
            IF(JJ.GT. 4) GO TO 242                                      
      E4(1,J,4,2)=+S7(3,K)+S5(5,L)*F03                                  
      E4(2,J,4,2)=+S7(4,K)+S5(6,L)*F03                                  
C                                                                       
            IF(JJ.GT. 3) GO TO 242                                      
            DO I=1,4                                                    
      E3(I,J,4,2)=+S6(I+ 8,K)+S4(I+13,L)*F03                            
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 2) GO TO 242                                      
            DO I=1,4                                                    
      E2(I,J,4,2)=+S5(I+20,K)+S3(I+26,L)*F03                            
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 1) GO TO 242                                      
            M= IN6(L)                                                   
            DO I=1,5                                                    
      E1(I  ,4,2)=+S4(I+37,K)+S2(I+36,M)*F03                            
            ENDDO                                                       
  242       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 5  and  L = 2                   
Cjms                                                                    
      J= 0                                                              
      DO JJ=1,5                                                         
         DO II=1,JJ                                                     
            J= J+1                                                      
            K= IND(J,5,2)                                               
            L= K-3-JJ                                                   
            M= L-JJ                                                     
C                                                                       
      E5(  J,5,2)=+R8(K  )-R7(L,1)+R6(M,1)-R5(J,1)                      
C                                                                       
            IF(JJ.GT. 4) GO TO 252                                      
      E4(1,J,5,2)=+S7(3,K)-S6(5,L)+S5(5,M)-S4(6,J)                      
      E4(2,J,5,2)=+S7(4,K)-S6(6,L)+S5(6,M)-S4(7,J)                      
C                                                                       
            IF(JJ.GT. 3) GO TO 252                                      
            DO I=1,4                                                    
      E3(I,J,5,2)=+S6(I+ 8,K)-S5(I+12,L)+S4(I+13,M)-S3(I+10,J)          
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 2) GO TO 252                                      
            N= IN6(J)                                                   
            DO I=1,4                                                    
      E2(I,J,5,2)=+S5(I+20,K)-S4(I+29,L)+S3(I+26,M)-S2(I+20,N)          
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 1) GO TO 252                                      
            M= IN6(M)                                                   
            DO I=1,5                                                    
      E1(I  ,5,2)=+S4(I+37,K)-S3(I+42,L)+S2(I+36,M)-S1(I+20,J)          
            ENDDO                                                       
  252       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 6  and  L = 2                   
Cjms                                                                    
      J= 0                                                              
      DO JJ=1,5                                                         
         DO II=1,JJ                                                     
            J= J+1                                                      
            K= IND(J,6,2)                                               
            L= K-4-JJ                                                   
            M= L-1-JJ                                                   
            N= M-2-JJ                                                   
C                                                                       
      E5(  J,6,2)=+R8(K  )-R7(L,1)+R6(M,1)*F03-R5(N,1)*F03              
C                                                                       
            IF(JJ.GT. 4) GO TO 262                                      
      E4(1,J,6,2)=+S7(3,K)-S6(5,L)+S5(5,M)*F03-S4(6,N)*F03              
      E4(2,J,6,2)=+S7(4,K)-S6(6,L)+S5(6,M)*F03-S4(7,N)*F03              
C                                                                       
            IF(JJ.GT. 3) GO TO 262                                      
            DO I=1,4                                                    
            I4= I+10                                                    
      E3(I,J,6,2)=+S6(I+ 8,K)-S5(I+12,L)+S4(I+13,M)*F03-S3(I4,N)*F03    
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 2) GO TO 262                                      
            N= IN6(N)                                                   
            DO I=1,4                                                    
            I4= I+20                                                    
      E2(I,J,6,2)=+S5(I+20,K)-S4(I+29,L)+S3(I+26,M)*F03-S2(I4,N)*F03    
            ENDDO                                                       
            N= M-2-JJ                                                   
C                                                                       
            IF(JJ.GT. 1) GO TO 262                                      
            M= IN6(M)                                                   
            DO I=1,5                                                    
            I4= I+20                                                    
      E1(I  ,6,2)=+S4(I+37,K)-S3(I+42,L)+S2(I+36,M)*F03-S1(I4,N)*F03    
            ENDDO                                                       
  262       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 1  and  L = 3                   
Cjms                                                                    
      J1=4                                                              
      DO J=1,15                                                         
         E5(  J,1,3)= T531E5(  J)+U5E5(  J,J1)+V531E5(  J)              
      ENDDO                                                             
      DO J=1,10                                                         
         E4(1,J,1,3)= T531E4(1,J)+U5E4(1,J,J1)+V531E4(1,J)              
         E4(2,J,1,3)= T531E4(2,J)+U5E4(2,J,J1)+V531E4(2,J)              
      ENDDO                                                             
      DO J=1,6                                                          
         E3(1,J,1,3)= T531E3(1,J)+U5E3(1,J,J1)+V531E3(1,J)              
         E3(2,J,1,3)= T531E3(2,J)+U5E3(2,J,J1)+V531E3(2,J)              
         E3(3,J,1,3)= T531E3(3,J)+U5E3(3,J,J1)+V531E3(3,J)              
         E3(4,J,1,3)= T531E3(4,J)+U5E3(4,J,J1)+V531E3(4,J)              
      ENDDO                                                             
      DO J=1,3                                                          
         E2(1,J,1,3)= T531E2(1,J)+U5E2(1,J,J1)+V531E2(1,J)              
         E2(2,J,1,3)= T531E2(2,J)+U5E2(2,J,J1)+V531E2(2,J)              
         E2(3,J,1,3)= T531E2(3,J)+U5E2(3,J,J1)+V531E2(3,J)              
         E2(4,J,1,3)= T531E2(4,J)+U5E2(4,J,J1)+V531E2(4,J)              
      ENDDO                                                             
      DO I=1,5                                                          
         E1(I  ,1,3)= T531E1(I  )+U5E1(I  ,J1)+V531E1(I  )              
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 2  and  L = 3                   
Cjms                                                                    
      J1=4                                                              
      DO J=1,15                                                         
         E5(  J,2,3)= T632E5(  J)+U6E5(  J,J1)+V632E5(  J)              
      ENDDO                                                             
      DO J=1,10                                                         
         E4(1,J,2,3)= T632E4(1,J)+U6E4(1,J,J1)+V632E4(1,J)              
         E4(2,J,2,3)= T632E4(2,J)+U6E4(2,J,J1)+V632E4(2,J)              
      ENDDO                                                             
      DO J=1,6                                                          
         E3(1,J,2,3)= T632E3(1,J)+U6E3(1,J,J1)+V632E3(1,J)              
         E3(2,J,2,3)= T632E3(2,J)+U6E3(2,J,J1)+V632E3(2,J)              
         E3(3,J,2,3)= T632E3(3,J)+U6E3(3,J,J1)+V632E3(3,J)              
         E3(4,J,2,3)= T632E3(4,J)+U6E3(4,J,J1)+V632E3(4,J)              
      ENDDO                                                             
      DO J=1,3                                                          
         E2(1,J,2,3)= T632E2(1,J)+U6E2(1,J,J1)+V632E2(1,J)              
         E2(2,J,2,3)= T632E2(2,J)+U6E2(2,J,J1)+V632E2(2,J)              
         E2(3,J,2,3)= T632E2(3,J)+U6E2(3,J,J1)+V632E2(3,J)              
         E2(4,J,2,3)= T632E2(4,J)+U6E2(4,J,J1)+V632E2(4,J)              
      ENDDO                                                             
      DO I=1,5                                                          
         E1(I  ,2,3)= T632E1(I  )+U6E1(I  ,J1)+V632E1(I  )              
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 3  and  L = 3                   
Cjms                                                                    
      J= 0                                                              
      DO JJ=1,5                                                         
         DO II=1,JJ                                                     
            J= J+1                                                      
            K= IND(J,3,3)                                               
            L= K-4-JJ                                                   
            M= L-3-JJ                                                   
            N= M-2-JJ                                                   
C                                                                       
      E5(  J,3,3)=    +R8(K  )    -S7(1,L)*F02-S7(2,L)*F02              
     *    +S6(1,M)*F06+S6(2,M)    +S6(3,M)*F04+S6(4,M)+                 
     *   (-S5(1,N)*F03-S5(2,N)*F03-S5(3,N)    -S5(4,N))*F02             
     *    +S4(1,J)*F03+S4(2,J)    +S4(3,J)*F04+S4(4,J)+S4(5,J)          
C                                                                       
            IF(JJ.GT. 4) GO TO 333                                      
            DO I=1,2                                                    
               I4= I+ 4                                                 
               I5= I+ 5                                                 
      E4(I,J,3,3)=    +S7(I+ 2,K)    -S6(I+ 4,L)*F02-S6(I+ 6,L)*F02     
     *   +S5(I4,M)*F06+S5(I4+2,M)    +S5(I4+4,M)*F04+S5(I4+6,M)+        
     *  (-S4(I5,N)*F03-S4(I5+2,N)*F03-S4(I5+4,N)    -S4(I5+6,N))*F02    
     *   +S3(I   ,J)*F03+S3(I+ 2,J)                                     
     *   +S3(I+ 4,J)*F04+S3(I+ 6,J)+S3(I+ 8,J)                          
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 3) GO TO 333                                      
            JN= IN6(J)                                                  
            DO I=1,4                                                    
               I4= I+13                                                 
               I5= I+10                                                 
      E3(I,J,3,3)=    +S6(I+ 8,K)    -S5(I+12,L)*F02-S5(I+16,L)*F02     
     *   +S4(I4,M)*F06+S4(I4+4,M)    +S4(I4+8,M)*F04+S4(I4+12,M)+       
     *  (-S3(I5,N)*F03-S3(I5+4,N)*F03-S3(I5+8,N)    -S3(I5+12,N))*F02   
     *   +S2(I   ,JN)*F03+S2(I+ 4,JN)                                   
     *   +S2(I+ 8,JN)*F04+S2(I+12,JN)+S2(I+16,JN)                       
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 2) GO TO 333                                      
            N= IN6(N)                                                   
            DO I=1,4                                                    
               I4= I+26                                                 
               I5= I+20                                                 
      E2(I,J,3,3)=    +S5(I+20,K)    -S4(I+29,L)*F02-S4(I+33,L)*F02     
     *   +S3(I4,M)*F06+S3(I4+4,M)    +S3(I4+8,M)*F04+S3(I4+12,M)+       
     *  (-S2(I5,N)*F03-S2(I5+4,N)*F03-S2(I5+8,N)    -S2(I5+12,N))*F02   
     *   +S1(I   ,J)*F03+S1(I+ 4,J)                                     
     *   +S1(I+ 8,J)*F04+S1(I+12,J)+S1(I+16,J)                          
            ENDDO                                                       
            N= M-2-JJ                                                   
C                                                                       
            IF(JJ.GT. 1) GO TO 333                                      
            M= IN6(M)                                                   
            DO I=1,5                                                    
               I4= I+36                                                 
               I5= I+20                                                 
      E1(I  ,3,3)=    +S4(I+37,K)    -S3(I+42,L)*F02-S3(I+47,L)*F02     
     *   +S2(I4,M)*F06+S2(I4+5,M)    +S2(I4+10,M)*F04+S2(I4+15,M)+      
     *  (-S1(I5,N)*F03-S1(I5+5,N)*F03-S1(I5+10,N)    -S1(I5+15,N))*F02  
     *   +R0(I   )*F03+R0(I+ 5)      +R0(I+10)*F04   +R0(I+15)+R0(I+20) 
            ENDDO                                                       
  333       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 4  and  L = 3                   
Cjms                                                                    
      J= 0                                                              
      DO JJ=1,5                                                         
         DO II=1,JJ                                                     
            J= J+1                                                      
            K= IND(J,4,3)                                               
            L= K-3-JJ                                                   
            M= L-2-JJ                                                   
C                                                                       
      E5(  J,4,3)=+R8(K  )-R7(L,2)*F02+S6(1,M)+S6(4,M)                  
C                                                                       
            IF(JJ.GT. 4) GO TO 343                                      
      E4(1,J,4,3)=+S7(3,K)-S6(7,L)*F02+S5(5,M)+S5(11,M)                 
      E4(2,J,4,3)=+S7(4,K)-S6(8,L)*F02+S5(6,M)+S5(12,M)                 
C                                                                       
            IF(JJ.GT. 3) GO TO 343                                      
            DO I=1,4                                                    
      E3(I,J,4,3)=+S6(I+ 8,K)-S5(I+16,L)*F02+S4(I+13,M)+S4(I+25,M)      
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 2) GO TO 343                                      
            DO I=1,4                                                    
      E2(I,J,4,3)=+S5(I+20,K)-S4(I+33,L)*F02+S3(I+26,M)+S3(I+38,M)      
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 1) GO TO 343                                      
            M= IN6(M)                                                   
            DO I=1,5                                                    
      E1(I  ,4,3)=+S4(I+37,K)-S3(I+47,L)*F02+S2(I+36,M)+S2(I+51,M)      
            ENDDO                                                       
  343       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 5  and  L = 3                   
Cjms                                                                    
      J= 0                                                              
      DO JJ=1,5                                                         
         DO II=1,JJ                                                     
            J= J+1                                                      
            K= IND(J,5,3)                                               
            L= K-3-JJ                                                   
            M= L-2-JJ                                                   
C                                                                       
      E5(  J,5,3)=+R8(K  )-S7(1,L)    -S7(2,L)*F02                      
     *                    +S6(1,M)*F03+S6(3,M)*F02+S6(4,M)              
     *                    -S5(1,J)    -S5(2,J)*F02-S5(4,J)              
C                                                                       
            IF(JJ.GT. 4) GO TO 353                                      
            DO I=1,2                                                    
      E4(I,J,5,3)=+S7(I+ 2,K)-S6(I+ 4,L)    -S6(I+ 6,L)*F02             
     *                       +S5(I+ 4,M)*F03+S5(I+ 8,M)*F02+S5(I+10,M)  
     *                       -S4(I+ 5,J)    -S4(I+ 7,J)*F02-S4(I+11,J)  
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 3) GO TO 353                                      
            DO I=1,4                                                    
      E3(I,J,5,3)=+S6(I+ 8,K)-S5(I+12,L)    -S5(I+16,L)*F02             
     *                       +S4(I+13,M)*F03+S4(I+21,M)*F02+S4(I+25,M)  
     *                       -S3(I+10,J)    -S3(I+14,J)*F02-S3(I+22,J)  
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 2) GO TO 353                                      
            N= IN6(J)                                                   
            DO I=1,4                                                    
      E2(I,J,5,3)=+S5(I+20,K)-S4(I+29,L)    -S4(I+33,L)*F02             
     *                       +S3(I+26,M)*F03+S3(I+34,M)*F02+S3(I+38,M)  
     *                       -S2(I+20,N)    -S2(I+24,N)*F02-S2(I+32,N)  
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 1) GO TO 353                                      
            M= IN6(M)                                                   
            DO I=1,5                                                    
      E1(I  ,5,3)=+S4(I+37,K)-S3(I+42,L)    -S3(I+47,L)*F02             
     *                       +S2(I+36,M)*F03+S2(I+46,M)*F02+S2(I+51,M)  
     *                       -S1(I+20,J)    -S1(I+25,J)*F02-S1(I+35,J)  
            ENDDO                                                       
  353       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 6  and  L = 3                   
Cjms                                                                    
      J= 0                                                              
      DO JJ=1,5                                                         
         DO II=1,JJ                                                     
            J= J+1                                                      
            K= IND(J,6,3)                                               
            L= K-4-JJ                                                   
            M= L-3-JJ                                                   
            N= M-2-JJ                                                   
C                                                                       
            IF(JJ.LT. 5) THEN                                           
      E5(  J,6,3)= E5(J+JJ,5,3)                                         
            ELSE                                                        
      E5(  J,6,3)=+R8(K  )-S7(1,L)    -S7(2,L)*F02                      
     *                    +S6(1,M)*F03+S6(3,M)*F02+S6(4,M)              
     *                    -S5(1,N)    -S5(2,N)*F02-S5(4,N)              
            ENDIF                                                       
C                                                                       
            IF(JJ.GT. 4) GO TO 363                                      
            IF(JJ.LT. 4) THEN                                           
      E4(1,J,6,3)= E4(1,J+JJ,5,3)                                       
      E4(2,J,6,3)= E4(2,J+JJ,5,3)                                       
            ELSE                                                        
               DO I=1,2                                                 
      E4(I,J,6,3)=+S7(I+ 2,K)-S6(I+ 4,L)    -S6(I+ 6,L)*F02             
     *                       +S5(I+ 4,M)*F03+S5(I+ 8,M)*F02+S5(I+10,M)  
     *                       -S4(I+ 5,N)    -S4(I+ 7,N)*F02-S4(I+11,N)  
               ENDDO                                                    
            ENDIF                                                       
C                                                                       
            IF(JJ.GT. 3) GO TO 363                                      
            IF(JJ.LT. 3) THEN                                           
               DO I=1,4                                                 
      E3(I,J,6,3)= E3(I,J+JJ,5,3)                                       
               ENDDO                                                    
            ELSE                                                        
               DO I=1,4                                                 
      E3(I,J,6,3)=+S6(I+ 8,K)-S5(I+12,L)    -S5(I+16,L)*F02             
     *                       +S4(I+13,M)*F03+S4(I+21,M)*F02+S4(I+25,M)  
     *                       -S3(I+10,N)    -S3(I+14,N)*F02-S3(I+22,N)  
               ENDDO                                                    
            ENDIF                                                       
C                                                                       
            IF(JJ.GT. 2) GO TO 363                                      
            IF(JJ.LT. 2) THEN                                           
               DO I=1,4                                                 
      E2(I,J,6,3)= E2(I,J+JJ,5,3)                                       
               ENDDO                                                    
            ELSE                                                        
               N= IN6(N)                                                
               DO I=1,4                                                 
      E2(I,J,6,3)=+S5(I+20,K)-S4(I+29,L)    -S4(I+33,L)*F02             
     *                       +S3(I+26,M)*F03+S3(I+34,M)*F02+S3(I+38,M)  
     *                       -S2(I+20,N)    -S2(I+24,N)*F02-S2(I+32,N)  
               ENDDO                                                    
               N= M-2-JJ                                                
            ENDIF                                                       
C                                                                       
            IF(JJ.GT. 1) GO TO 363                                      
            M= IN6(M)                                                   
            DO I=1,5                                                    
      E1(I  ,6,3)=+S4(I+37,K)-S3(I+42,L)    -S3(I+47,L)*F02             
     *                       +S2(I+36,M)*F03+S2(I+46,M)*F02+S2(I+51,M)  
     *                       -S1(I+20,N)    -S1(I+25,N)*F02-S1(I+35,N)  
            ENDDO                                                       
  363       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 1  and  L = 4                   
C                                                                       
C     next statements commented because F(I,J,1,4)                      
C     will be obtained in FIJKL5 from   F(I,J,4,1)                      
Cjms                                                                    
C     DO 415 J=1,15                                                     
C 415 E5(  J,1,4)= E5(  J,4,1)                                          
C     DO 414 J=1,10                                                     
C        DO 414 I=1,2                                                   
C 414 E4(I,J,1,4)= E4(I,J,4,1)                                          
C     DO 413 J=1,6                                                      
C        DO 413 I=1,4                                                   
C 413 E3(I,J,1,4)= E3(I,J,4,1)                                          
C     DO 412 J=1,3                                                      
C        DO 412 I=1,4                                                   
C 412 E2(I,J,1,4)= E2(I,J,4,1)                                          
C        DO 411 I=1,5                                                   
C 411 E1(I  ,1,4)= E1(I  ,4,1)                                          
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 2  and  L = 4                   
C                                                                       
C     next statements commented because F(I,J,2,4)                      
C     will be obtained in FIJKL5 from   F(I,J,4,2)                      
Cjms                                                                    
C     DO 425 J=1,15                                                     
C 425 E5(  J,2,4)= E5(  J,4,2)                                          
C     DO 424 J=1,10                                                     
C        DO 424 I=1,2                                                   
C 424 E4(I,J,2,4)= E4(I,J,4,2)                                          
C     DO 423 J=1,6                                                      
C        DO 423 I=1,4                                                   
C 423 E3(I,J,2,4)= E3(I,J,4,2)                                          
C     DO 422 J=1,3                                                      
C        DO 422 I=1,4                                                   
C 422 E2(I,J,2,4)= E2(I,J,4,2)                                          
C        DO 421 I=1,5                                                   
C 421 E1(I  ,2,4)= E1(I  ,4,2)                                          
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 3  and  L = 4                   
Cjms                                                                    
      J= 0                                                              
      DO JJ=1,5                                                         
         DO II=1,JJ                                                     
            J= J+1                                                      
            K= IND(J,3,4)                                               
            L= K-3-JJ                                                   
            M= L-2-JJ                                                   
C                                                                       
      E5(  J,3,4)=+R8(K  )-R7(L,1)*F02+S6(1,M)+S6(2,M)                  
C                                                                       
            IF(JJ.GT. 4) GO TO 434                                      
      E4(1,J,3,4)=+S7(3,K)-S6(5,L)*F02+S5(5,M)+S5(7,M)                  
      E4(2,J,3,4)=+S7(4,K)-S6(6,L)*F02+S5(6,M)+S5(8,M)                  
C                                                                       
            IF(JJ.GT. 3) GO TO 434                                      
            DO I=1,4                                                    
      E3(I,J,3,4)=+S6(I+ 8,K)-S5(I+12,L)*F02+S4(I+13,M)+S4(I+17,M)      
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 2) GO TO 434                                      
            DO I=1,4                                                    
      E2(I,J,3,4)=+S5(I+20,K)-S4(I+29,L)*F02+S3(I+26,M)+S3(I+30,M)      
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 1) GO TO 434                                      
            M= IN6(M)                                                   
            DO I=1,5                                                    
      E1(I  ,3,4)=+S4(I+37,K)-S3(I+42,L)*F02+S2(I+36,M)+S2(I+41,M)      
            ENDDO                                                       
  434       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 4  and  L = 4                   
C                                                                       
C     next statements commented because F(I,J,4,4)                      
C     will be obtained in FIJKL5 from   F(I,J,2,1)                      
Cjms                                                                    
C     DO 445 J=1,15                                                     
C 445 E5(  J,4,4)= E5(  J,2,1)                                          
C     DO 444 J=1,10                                                     
C        DO 444 I=1,2                                                   
C 444 E4(I,J,4,4)= E4(I,J,2,1)                                          
C     DO 443 J=1,6                                                      
C        DO 443 I=1,4                                                   
C 443 E3(I,J,4,4)= E3(I,J,2,1)                                          
C     DO 442 J=1,3                                                      
C        DO 442 I=1,4                                                   
C 442 E2(I,J,4,4)= E2(I,J,2,1)                                          
C        DO 441 I=1,5                                                   
C 441 E1(I  ,4,4)= E1(I  ,2,1)                                          
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 5  and  L = 4                   
C                                                                       
C     next statements commented because F(I,J,5,4)                      
C     will be obtained in FIJKL5 from   F(I,J,6,1)                      
Cjms                                                                    
C     DO 455 J=1,15                                                     
C 455 E5(  J,5,4)= E5(  J,6,1)                                          
C     DO 454 J=1,10                                                     
C        DO 454 I=1,2                                                   
C 454 E4(I,J,5,4)= E4(I,J,6,1)                                          
C     DO 453 J=1,6                                                      
C        DO 453 I=1,4                                                   
C 453 E3(I,J,5,4)= E3(I,J,6,1)                                          
C     DO 452 J=1,3                                                      
C        DO 452 I=1,4                                                   
C 452 E2(I,J,5,4)= E2(I,J,6,1)                                          
C        DO 451 I=1,5                                                   
C 451 E1(I  ,5,4)= E1(I  ,6,1)                                          
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 6  and  L = 4                   
C                                                                       
C     next statements commented because F(I,J,6,4)                      
C     will be obtained in FIJKL5 from   F(I,J,5,2)                      
Cjms                                                                    
C     DO 465 J=1,15                                                     
C 465 E5(  J,6,4)= E5(  J,5,2)                                          
C     DO 464 J=1,10                                                     
C        DO 464 I=1,2                                                   
C 464 E4(I,J,6,4)= E4(I,J,5,2)                                          
C     DO 463 J=1,6                                                      
C        DO 463 I=1,4                                                   
C 463 E3(I,J,6,4)= E3(I,J,5,2)                                          
C     DO 462 J=1,3                                                      
C        DO 462 I=1,4                                                   
C 462 E2(I,J,6,4)= E2(I,J,5,2)                                          
C        DO 461 I=1,5                                                   
C 461 E1(I  ,6,4)= E1(I  ,5,2)                                          
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 1  and  L = 5                   
Cjms                                                                    
      DO J=1,15                                                         
         K= IND(J,1,5)                                                  
         E5(  J,1,5)=+R8(   K)-R7( J,2)+(R6( K,1)-R5( J,2))*F03         
      ENDDO                                                             
C                                                                       
      DO J=1,10                                                         
         K= IND(J,1,5)                                                  
         E4(1,J,1,5)=+S7( 3,K)-S6( 7,J)+(S5( 5,K)-S4( 8,J))*F03         
         E4(2,J,1,5)=+S7( 4,K)-S6( 8,J)+(S5( 6,K)-S4( 9,J))*F03         
      ENDDO                                                             
C                                                                       
      DO J=1,6                                                          
         K= IND(J,1,5)                                                  
         E3(1,J,1,5)=+S6( 9,K)-S5(17,J)+(S4(14,K)-S3(15,J))*F03         
         E3(2,J,1,5)=+S6(10,K)-S5(18,J)+(S4(15,K)-S3(16,J))*F03         
         E3(3,J,1,5)=+S6(11,K)-S5(19,J)+(S4(16,K)-S3(17,J))*F03         
         E3(4,J,1,5)=+S6(12,K)-S5(20,J)+(S4(17,K)-S3(18,J))*F03         
      ENDDO                                                             
C                                                                       
      DO J=1,3                                                          
         K= IND(J,1,5)                                                  
         M= IN6(J)                                                      
         E2(1,J,1,5)=+S5(21,K)-S4(34,J)+(S3(27,K)-S2(25,M))*F03         
         E2(2,J,1,5)=+S5(22,K)-S4(35,J)+(S3(28,K)-S2(26,M))*F03         
         E2(3,J,1,5)=+S5(23,K)-S4(36,J)+(S3(29,K)-S2(27,M))*F03         
         E2(4,J,1,5)=+S5(24,K)-S4(37,J)+(S3(30,K)-S2(28,M))*F03         
      ENDDO                                                             
C                                                                       
         J=1                                                            
         K= IND(J,1,5)                                                  
         M= IN6(K)                                                      
      DO I=1,5                                                          
         E1(I  ,1,5)=+S4(I+37,K)-S3(I+47,J)+(S2(I+36,M)-S1(I+25,J))*F03 
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 2  and  L = 5                   
Cjms                                                                    
      J= 0                                                              
      DO JJ=1,5                                                         
         DO II=1,JJ                                                     
            J= J+1                                                      
            K= IND(J,2,5)                                               
            L= K-3-JJ                                                   
            M= L-JJ                                                     
C                                                                       
      E5(  J,2,5)=+R8(K  )-R7(L,2)+R6(M,1)-R5(J,2)                      
C                                                                       
            IF(JJ.GT. 4) GO TO 525                                      
      E4(1,J,2,5)=+S7(3,K)-S6(7,L)+S5(5,M)-S4(8,J)                      
      E4(2,J,2,5)=+S7(4,K)-S6(8,L)+S5(6,M)-S4(9,J)                      
C                                                                       
            IF(JJ.GT. 3) GO TO 525                                      
            DO I=1,4                                                    
      E3(I,J,2,5)=+S6(I+ 8,K)-S5(I+16,L)+S4(I+13,M)-S3(I+14,J)          
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 2) GO TO 525                                      
            N= IN6(J)                                                   
            DO I=1,4                                                    
      E2(I,J,2,5)=+S5(I+20,K)-S4(I+33,L)+S3(I+26,M)-S2(I+24,N)          
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 1) GO TO 525                                      
            M= IN6(M)                                                   
            DO I=1,5                                                    
      E1(I  ,2,5)=+S4(I+37,K)-S3(I+47,L)+S2(I+36,M)-S1(I+25,J)          
            ENDDO                                                       
  525       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 3  and  L = 5                   
Cjms                                                                    
      J= 0                                                              
      DO JJ=1,5                                                         
         DO II=1,JJ                                                     
            J= J+1                                                      
            K= IND(J,3,5)                                               
            L= K-3-JJ                                                   
            M= L-2-JJ                                                   
C                                                                       
      E5(  J,3,5)=+R8(K  )-S7(1,L)*F02-S7(2,L)                          
     *                    +S6(1,M)*F03+S6(2,M)+S6(3,M)*F02              
     *                    -S5(1,J)*F02-S5(2,J)-S5(3,J)                  
C                                                                       
            IF(JJ.GT. 4) GO TO 535                                      
            DO I=1,2                                                    
      E4(I,J,3,5)=+S7(I+ 2,K)-S6(I+ 4,L)*F02-S6(I+ 6,L)                 
     *                       +S5(I+ 4,M)*F03+S5(I+ 6,M)+S5(I+ 8,M)*F02  
     *                       -S4(I+ 5,J)*F02-S4(I+ 7,J)-S4(I+ 9,J)      
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 3) GO TO 535                                      
            DO I=1,4                                                    
      E3(I,J,3,5)=+S6(I+ 8,K)-S5(I+12,L)*F02-S5(I+16,L)                 
     *                       +S4(I+13,M)*F03+S4(I+17,M)+S4(I+21,M)*F02  
     *                       -S3(I+10,J)*F02-S3(I+14,J)-S3(I+18,J)      
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 2) GO TO 535                                      
            N= IN6(J)                                                   
            DO I=1,4                                                    
      E2(I,J,3,5)=+S5(I+20,K)-S4(I+29,L)*F02-S4(I+33,L)                 
     *                       +S3(I+26,M)*F03+S3(I+30,M)+S3(I+34,M)*F02  
     *                       -S2(I+20,N)*F02-S2(I+24,N)-S2(I+28,N)      
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 1) GO TO 535                                      
            M= IN6(M)                                                   
            DO I=1,5                                                    
      E1(I  ,3,5)=+S4(I+37,K)-S3(I+42,L)*F02-S3(I+47,L)                 
     *                       +S2(I+36,M)*F03+S2(I+41,M)+S2(I+46,M)*F02  
     *                       -S1(I+20,J)*F02-S1(I+25,J)-S1(I+30,J)      
            ENDDO                                                       
  535       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 4  and  L = 5                   
Cjms                                                                    
      J= 0                                                              
      DO JJ=1,5                                                         
         DO II=1,JJ                                                     
            J= J+1                                                      
            K= IND(J,4,5)                                               
            L= K-2-JJ                                                   
C                                                                       
      E5(  J,4,5)=+R8(K  )-R7(L,2)+R6(K,1)-R5(L,2)                      
C                                                                       
            IF(JJ.GT. 4) GO TO 545                                      
      E4(1,J,4,5)=+S7(3,K)-S6(7,L)+S5(5,K)-S4(8,L)                      
      E4(2,J,4,5)=+S7(4,K)-S6(8,L)+S5(6,K)-S4(9,L)                      
C                                                                       
            IF(JJ.GT. 3) GO TO 545                                      
            DO I=1,4                                                    
      E3(I,J,4,5)=+S6(I+ 8,K)-S5(I+16,L)+S4(I+13,K)-S3(I+14,L)          
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 2) GO TO 545                                      
            M= IN6(L)                                                   
            DO I=1,4                                                    
      E2(I,J,4,5)=+S5(I+20,K)-S4(I+33,L)+S3(I+26,K)-S2(I+24,M)          
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 1) GO TO 545                                      
            M= IN6(K)                                                   
            DO I=1,5                                                    
      E1(I  ,4,5)=+S4(I+37,K)-S3(I+47,L)+S2(I+36,M)-S1(I+25,L)          
            ENDDO                                                       
  545       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 5  and  L = 5                   
Cjms                                                                    
C                                                                       
      J1=3                                                              
      DO J=1,15                                                         
         E5(  J,5,5)= T531E5(  J)+U5E5(  J,J1)                          
      ENDDO                                                             
      DO J=1,10                                                         
         E4(1,J,5,5)= T531E4(1,J)+U5E4(1,J,J1)                          
         E4(2,J,5,5)= T531E4(2,J)+U5E4(2,J,J1)                          
      ENDDO                                                             
      DO J=1,6                                                          
         E3(1,J,5,5)= T531E3(1,J)+U5E3(1,J,J1)                          
         E3(2,J,5,5)= T531E3(2,J)+U5E3(2,J,J1)                          
         E3(3,J,5,5)= T531E3(3,J)+U5E3(3,J,J1)                          
         E3(4,J,5,5)= T531E3(4,J)+U5E3(4,J,J1)                          
      ENDDO                                                             
      DO J=1,3                                                          
         E2(1,J,5,5)= T531E2(1,J)+U5E2(1,J,J1)                          
         E2(2,J,5,5)= T531E2(2,J)+U5E2(2,J,J1)                          
         E2(3,J,5,5)= T531E2(3,J)+U5E2(3,J,J1)                          
         E2(4,J,5,5)= T531E2(4,J)+U5E2(4,J,J1)                          
      ENDDO                                                             
      DO I=1,5                                                          
         E1(I  ,5,5)= T531E1(I  )+U5E1(I  ,J1)                          
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 6  and  L = 5                   
Cjms                                                                    
      J= 0                                                              
      DO JJ=1,5                                                         
         DO II=1,JJ                                                     
            J= J+1                                                      
            K= IND(J,6,5)                                               
            L= K-3-JJ                                                   
            M= L-2-JJ                                                   
C                                                                       
      E5(  J,6,5)=+R8(K  )-S7(1,L)-S7(2,L)+S6(1,M)+S6(3,M)              
C                                                                       
            IF(JJ.GT. 4) GO TO 565                                      
      E4(1,J,6,5)=+S7(3,K)-S6(5,L)-S6(7,L)+S5(5,M)+S5( 9,M)             
      E4(2,J,6,5)=+S7(4,K)-S6(6,L)-S6(8,L)+S5(6,M)+S5(10,M)             
C                                                                       
            IF(JJ.GT. 3) GO TO 565                                      
            DO I=1,4                                                    
      E3(I,J,6,5)=+S6(I+ 8,K)-S5(I+12,L)-S5(I+16,L)                     
     *                       +S4(I+13,M)+S4(I+21,M)                     
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 2) GO TO 565                                      
            DO I=1,4                                                    
      E2(I,J,6,5)=+S5(I+20,K)-S4(I+29,L)-S4(I+33,L)                     
     *                       +S3(I+26,M)+S3(I+34,M)                     
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 1) GO TO 565                                      
            M= IN6(M)                                                   
            DO I=1,5                                                    
      E1(I  ,6,5)=+S4(I+37,K)-S3(I+42,L)-S3(I+47,L)                     
     *                       +S2(I+36,M)+S2(I+46,M)                     
            ENDDO                                                       
  565       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 1  and  L = 6                   
C                                                                       
C     next statements commented because F(I,J,1,6)                      
C     will be obtained in FIJKL5 from   F(I,J,4,5)                      
Cjms                                                                    
C     DO 615 J=1,15                                                     
C 615 E5(  J,1,6)= E5(  J,4,5)                                          
C     DO 614 J=1,10                                                     
C        DO 614 I=1,2                                                   
C 614 E4(I,J,1,6)= E4(I,J,4,5)                                          
C     DO 613 J=1,6                                                      
C        DO 613 I=1,4                                                   
C 613 E3(I,J,1,6)= E3(I,J,4,5)                                          
C     DO 612 J=1,3                                                      
C        DO 612 I=1,4                                                   
C 612 E2(I,J,1,6)= E2(I,J,4,5)                                          
C        DO 611 I=1,5                                                   
C 611 E1(I  ,1,6)= E1(I  ,4,5)                                          
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 2  and  L = 6                   
Cjms                                                                    
      J= 0                                                              
      DO JJ=1,5                                                         
         DO II=1,JJ                                                     
            J= J+1                                                      
            K= IND(J,2,6)                                               
            L= K-4-JJ                                                   
            M= L-1-JJ                                                   
            N= M-2-JJ                                                   
C                                                                       
      E5(  J,2,6)=+R8(K  )-R7(L,2)+R6(M,1)*F03-R5(N,2)*F03              
C                                                                       
            IF(JJ.GT. 4) GO TO 626                                      
      E4(1,J,2,6)=+S7(3,K)-S6(7,L)+S5(5,M)*F03-S4(8,N)*F03              
      E4(2,J,2,6)=+S7(4,K)-S6(8,L)+S5(6,M)*F03-S4(9,N)*F03              
C                                                                       
            IF(JJ.GT. 3) GO TO 626                                      
            DO I=1,4                                                    
      E3(I,J,2,6)=+S6(I+8,K)-S5(I+16,L)+S4(I+13,M)*F03-S3(I+14,N)*F03   
            ENDDO                                                       
C                                                                       
            IF(JJ.GT. 2) GO TO 626                                      
            N= IN6(N)                                                   
            DO I=1,4                                                    
            I4= I+20                                                    
      E2(I,J,2,6)=+S5(I4 ,K)-S4(I+33,L)+S3(I+26,M)*F03-S2(I+24,N)*F03   
            ENDDO                                                       
            N= M-2-JJ                                                   
C                                                                       
            IF(JJ.GT. 1) GO TO 626                                      
            M= IN6(M)                                                   
            DO I=1,5                                                    
            I4= I+37                                                    
      E1(I  ,2,6)=+S4(I4 ,K)-S3(I+47,L)+S2(I+36,M)*F03-S1(I+25,N)*F03   
            ENDDO                                                       
  626       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 3  and  L = 6                   
Cjms                                                                    
      J= 0                                                              
      DO JJ=1,5                                                         
         DO II=1,JJ                                                     
            J= J+1                                                      
            K= IND(J,3,6)                                               
            L= K-4-JJ                                                   
            M= L-3-JJ                                                   
            N= M-2-JJ                                                   
C                                                                       
            IF(JJ.LT. 5) THEN                                           
      E5(  J,3,6)= E5(J+JJ,3,5)                                         
            ELSE                                                        
      E5(  J,3,6)=+R8(K  )    -S7(1,L)*F02-S7(2,L)                      
     *            +S6(1,M)*F03+S6(2,M)    +S6(3,M)*F02                  
     *            -S5(1,N)*F02-S5(2,N)    -S5(3,N)                      
            ENDIF                                                       
C                                                                       
            IF(JJ.GT. 4) GO TO 636                                      
            IF(JJ.LT. 4) THEN                                           
      E4(1,J,3,6)= E4(1,J+JJ,3,5)                                       
      E4(2,J,3,6)= E4(2,J+JJ,3,5)                                       
            ELSE                                                        
               DO I=1,2                                                 
      E4(I,J,3,6)=+S7(I+ 2,K)    -S6(I+ 4,L)*F02-S6(I+ 6,L)             
     *            +S5(I+ 4,M)*F03+S5(I+ 6,M)    +S5(I+ 8,M)*F02         
     *            -S4(I+ 5,N)*F02-S4(I+ 7,N)    -S4(I+ 9,N)             
               ENDDO                                                    
            ENDIF                                                       
C                                                                       
            IF(JJ.GT. 3) GO TO 636                                      
            IF(JJ.LT. 3) THEN                                           
               DO I=1,4                                                 
      E3(I,J,3,6)= E3(I,J+JJ,3,5)                                       
               ENDDO                                                    
            ELSE                                                        
               DO I=1,4                                                 
               I4= I+10                                                 
      E3(I,J,3,6)=+S6(I+ 8,K)    -S5(I+12,L)*F02-S5(I+16,L)             
     *            +S4(I+13,M)*F03+S4(I+17,M)    +S4(I+21,M)*F02         
     *            -S3(I+10,N)*F02-S3(I+14,N)    -S3(I+18,N)             
               ENDDO                                                    
            ENDIF                                                       
C                                                                       
            IF(JJ.GT. 2) GO TO 636                                      
            IF(JJ.LT. 2) THEN                                           
               DO I=1,4                                                 
      E2(I,J,3,6)= E2(I,J+JJ,3,5)                                       
               ENDDO                                                    
            ELSE                                                        
               N= IN6(N)                                                
               DO I=1,4                                                 
      E2(I,J,3,6)=+S5(I+20,K)    -S4(I+29,L)*F02-S4(I+33,L)             
     *            +S3(I+26,M)*F03+S3(I+30,M)    +S3(I+34,M)*F02         
     *            -S2(I+20,N)*F02-S2(I+24,N)    -S2(I+28,N)             
               ENDDO                                                    
               N= M-2-JJ                                                
            ENDIF                                                       
C                                                                       
            IF(JJ.GT. 1) GO TO 636                                      
            M= IN6(M)                                                   
            DO I=1,5                                                    
      E1(I  ,3,6)=+S4(I+37,K)    -S3(I+42,L)*F02-S3(I+47,L)             
     *            +S2(I+36,M)*F03+S2(I+41,M)    +S2(I+46,M)*F02         
     *            -S1(I+20,N)*F02-S1(I+25,N)    -S1(I+30,N)             
            ENDDO                                                       
  636       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 4  and  L = 6                   
C                                                                       
C     next statements commented because F(I,J,4,6)                      
C     will be obtained in FIJKL5 from   F(I,J,2,5)                      
Cjms                                                                    
C     DO 645 J=1,15                                                     
C 645 E5(  J,4,6)= E5(  J,2,5)                                          
C     DO 644 J=1,10                                                     
C        DO 644 I=1,2                                                   
C 644 E4(I,J,4,6)= E4(I,J,2,5)                                          
C     DO 643 J=1,6                                                      
C        DO 643 I=1,4                                                   
C 643 E3(I,J,4,6)= E3(I,J,2,5)                                          
C     DO 642 J=1,3                                                      
C        DO 642 I=1,4                                                   
C 642 E2(I,J,4,6)= E2(I,J,2,5)                                          
C        DO 641 I=1,5                                                   
C 641 E1(I  ,4,6)= E1(I  ,2,5)                                          
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 5  and  L = 6                   
C                                                                       
C     next statements commented because F(I,J,5,6)                      
C     will be obtained in FIJKL5 from   F(I,J,6,5)                      
Cjms                                                                    
C     DO 655 J=1,15                                                     
C 655 E5(  J,5,6)= E5(  J,6,5)                                          
C     DO 654 J=1,10                                                     
C        DO 654 I=1,2                                                   
C 654 E4(I,J,5,6)= E4(I,J,6,5)                                          
C     DO 653 J=1,6                                                      
C        DO 653 I=1,4                                                   
C 653 E3(I,J,5,6)= E3(I,J,6,5)                                          
C     DO 652 J=1,3                                                      
C        DO 652 I=1,4                                                   
C 652 E2(I,J,5,6)= E2(I,J,6,5)                                          
C        DO 651 I=1,5                                                   
C 651 E1(I  ,5,6)= E1(I  ,6,5)                                          
Cjms                                                                    
C     Auxiliary arrays to simplify the formulation of F(I,J,K,L)        
C     where:  I = 1..6,  J = 1..6,  K = 6  and  L = 6                   
Cjms                                                                    
      J1=3                                                              
      DO J=1,15                                                         
         E5(  J,6,6)= T632E5(  J)+U6E5(  J,J1)                          
      ENDDO                                                             
      DO J=1,10                                                         
         E4(1,J,6,6)= T632E4(1,J)+U6E4(1,J,J1)                          
         E4(2,J,6,6)= T632E4(2,J)+U6E4(2,J,J1)                          
      ENDDO                                                             
      DO J=1,6                                                          
         E3(1,J,6,6)= T632E3(1,J)+U6E3(1,J,J1)                          
         E3(2,J,6,6)= T632E3(2,J)+U6E3(2,J,J1)                          
         E3(3,J,6,6)= T632E3(3,J)+U6E3(3,J,J1)                          
         E3(4,J,6,6)= T632E3(4,J)+U6E3(4,J,J1)                          
      ENDDO                                                             
      DO J=1,3                                                          
         E2(1,J,6,6)= T632E2(1,J)+U6E2(1,J,J1)                          
         E2(2,J,6,6)= T632E2(2,J)+U6E2(2,J,J1)                          
         E2(3,J,6,6)= T632E2(3,J)+U6E2(3,J,J1)                          
         E2(4,J,6,6)= T632E2(4,J)+U6E2(4,J,J1)                          
      ENDDO                                                             
      DO I=1,5                                                          
         E1(I  ,6,6)= T632E1(I  )+U6E1(I  ,J1)                          
      ENDDO                                                             
C                                                                       
      CALL FIJKL5(KX,LX,F,QX,QZ,E1,E2,E3,E4,E5)                         
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2S   *DECK FIJKL4                                           
C>                                                                      
C>    @brief   auxiliary routine of order 4 for rot.axis integrations   
C>                                                                      
C>    @details auxiliary routine of order 4 for rot.axis integrations   
C>                                                                      
      SUBROUTINE FIJKL4(KX,LX,F,QX,QZ,D1,D2,D3,D4)                      
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      DIMENSION F(6,6,6,6)                                              
C                                                                       
      DIMENSION       D1(  5,KX,LX),D2(4,3,KX,LX),D3(3,6,KX,LX)         
      DIMENSION       D4( 10,KX,LX)                                     
C                                                                       
      LOGICAL   LSYM16,LSYM19                                           
C                                                                       
      PARAMETER (F02=2.0D+00)                                           
      PARAMETER (F03=3.0D+00)                                           
Cjms                                                                    
C     Simplified calculation of F(I,J,K,L) for cases                    
C     where:  I = 1..6,  J = 1..4,  K = 1..KX  and  L = 1..LX           
C     using auxiliary arrays D1, D2, D3 and D4                          
Cjms                                                                    
      XX= QX*QX                                                         
      ZZ= QZ*QZ                                                         
      XZ= QX*QZ                                                         
      XXX= XX*QX                                                        
      XXZ= XX*QZ                                                        
      XZZ= ZZ*QX                                                        
      ZZZ= ZZ*QZ                                                        
C                                                                       
      QXD= QX+QX                                                        
      QZD= QZ+QZ                                                        
C                                                                       
      XZD= XZ+XZ                                                        
C                                                                       
      LSYM16=(KX.EQ.4 .AND. LX.EQ.4)                                    
      LSYM19=(KX.EQ.6 .AND. LX.EQ.4)                                    
      DO L=1,LX                                                         
         DO K=1,KX                                                      
C                                                                       
            IF(LSYM16) THEN                                             
               IF(K.EQ.2 .AND. L.EQ.3) GO TO 001                        
            ENDIF                                                       
C                                                                       
            IF(LSYM19) THEN                                             
               IF(K.EQ.1 .AND. L.EQ.3) GO TO 001                        
               IF(K.EQ.4 .AND. L.EQ.3) GO TO 001                        
               IF(K.EQ.5 .AND. L.EQ.3) GO TO 001                        
            ENDIF                                                       
C                                                                       
            F(1,1,K,L)=+D3(1,1,K,L)+D1(  1,K,L)                         
     *                 +D2(1,1,K,L)*QXD                                 
     *                 +D1(  2,K,L)*XX                                  
            F(2,1,K,L)=+D3(1,4,K,L)+D1(  1,K,L)                         
            F(3,1,K,L)=+D3(1,6,K,L)+D1(  1,K,L)                         
     *                 +D2(1,3,K,L)*QZD                                 
     *                 +D1(  2,K,L)*ZZ                                  
            F(4,1,K,L)=+D3(1,2,K,L)                                     
     *                 +D2(1,2,K,L)*QX                                  
            F(5,1,K,L)=+D3(1,3,K,L)                                     
     *                 +D2(1,3,K,L)*QX                                  
     *                 +D2(1,1,K,L)*QZ                                  
     *                 +D1(  2,K,L)*XZ                                  
            F(6,1,K,L)=+D3(1,5,K,L)                                     
     *                 +D2(1,2,K,L)*QZ                                  
C                                                                       
            F(1,2,K,L)=+D4(  1,K,L)+D2(2,1,K,L)*F03                     
     *               +(+D3(2,1,K,L)*F02+D3(3,1,K,L)                     
     *                 +D1(  3,K,L)*F02+D1(  4,K,L))*QX                 
     *               +(+D2(3,1,K,L)+D2(4,1,K,L)*F02)*XX                 
     *                 +D1(  5,K,L)*XXX                                 
            F(2,2,K,L)=+D4(  4,K,L)+D2(2,1,K,L)                         
     *               +(+D3(3,4,K,L)+D1(  4,K,L))*QX                     
            F(3,2,K,L)=+D4(  6,K,L)+D2(2,1,K,L)                         
     *               +(+D3(3,6,K,L)+D1(  4,K,L))*QX                     
     *                 +D3(2,3,K,L)*QZD                                 
     *                 +D2(4,3,K,L)*XZD                                 
     *                 +D2(3,1,K,L)*ZZ                                  
     *                 +D1(  5,K,L)*XZZ                                 
            F(4,2,K,L)=+D4(  2,K,L)+D2(2,2,K,L)                         
     *               +(+D3(2,2,K,L)+D3(3,2,K,L))*QX                     
     *                 +D2(4,2,K,L)*XX                                  
            F(5,2,K,L)=+D4(  3,K,L)+D2(2,3,K,L)                         
     *               +(+D3(2,3,K,L)+D3(3,3,K,L))*QX                     
     *               +(+D3(2,1,K,L)+D1(  3,K,L))*QZ                     
     *                 +D2(4,3,K,L)*XX                                  
     *               +(+D2(3,1,K,L)+D2(4,1,K,L))*XZ                     
     *                 +D1(  5,K,L)*XXZ                                 
            F(6,2,K,L)=+D4(  5,K,L)                                     
     *                 +D3(3,5,K,L)*QX                                  
     *                 +D3(2,2,K,L)*QZ                                  
     *                 +D2(4,2,K,L)*XZ                                  
C                                                                       
            F(1,3,K,L)=+D4(  2,K,L)+D2(2,2,K,L)                         
     *                 +D3(2,2,K,L)*QXD                                 
     *                 +D2(3,2,K,L)*XX                                  
            F(2,3,K,L)=+D4(  7,K,L)+D2(2,2,K,L)*F03                     
            F(3,3,K,L)=+D4(  9,K,L)+D2(2,2,K,L)                         
     *                 +D3(2,5,K,L)*QZD                                 
     *                 +D2(3,2,K,L)*ZZ                                  
            F(4,3,K,L)=+D4(  4,K,L)+D2(2,1,K,L)                         
     *               +(+D3(2,4,K,L)+D1(  3,K,L))*QX                     
            F(5,3,K,L)=+D4(  5,K,L)                                     
     *                 +D3(2,5,K,L)*QX                                  
     *                 +D3(2,2,K,L)*QZ                                  
     *                 +D2(3,2,K,L)*XZ                                  
            F(6,3,K,L)=+D4(  8,K,L)+D2(2,3,K,L)                         
     *               +(+D3(2,4,K,L)+D1(  3,K,L))*QZ                     
C                                                                       
            F(1,4,K,L)=+D4(  3,K,L)+D2(2,3,K,L)                         
     *                 +D3(2,3,K,L)*QXD                                 
     *               +(+D3(3,1,K,L)+D1(  4,K,L))*QZ                     
     *                 +D2(3,3,K,L)*XX                                  
     *                 +D2(4,1,K,L)*XZD                                 
     *                 +D1(  5,K,L)*XXZ                                 
            F(2,4,K,L)=+D4(  8,K,L)+D2(2,3,K,L)                         
     *               +(+D3(3,4,K,L)+D1(  4,K,L))*QZ                     
            F(3,4,K,L)=+D4( 10,K,L)+D2(2,3,K,L)*F03                     
     *               +(+D3(2,6,K,L)*F02+D3(3,6,K,L)                     
     *                 +D1(  3,K,L)*F02+D1(  4,K,L))*QZ                 
     *               +(+D2(3,3,K,L)+D2(4,3,K,L)*F02)*ZZ                 
     *                 +D1(  5,K,L)*ZZZ                                 
            F(4,4,K,L)=+D4(  5,K,L)                                     
     *                 +D3(2,5,K,L)*QX                                  
     *                 +D3(3,2,K,L)*QZ                                  
     *                 +D2(4,2,K,L)*XZ                                  
            F(5,4,K,L)=+D4(  6,K,L)+D2(2,1,K,L)                         
     *               +(+D3(2,6,K,L)+D1(  3,K,L))*QX                     
     *               +(+D3(2,3,K,L)+D3(3,3,K,L))*QZ                     
     *               +(+D2(3,3,K,L)+D2(4,3,K,L))*XZ                     
     *                 +D2(4,1,K,L)*ZZ                                  
     *                 +D1(  5,K,L)*XZZ                                 
            F(6,4,K,L)=+D4(  9,K,L)+D2(2,2,K,L)                         
     *               +(+D3(2,5,K,L)+D3(3,5,K,L))*QZ                     
     *                 +D2(4,2,K,L)*ZZ                                  
C                                                                       
  001       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
C                                                                       
      IF(.NOT.LSYM16) GO TO 160                                         
      DO J=1,4                                                          
         DO I=1,6                                                       
            F(I,J,2,3)= F(I,J,3,2)                                      
         ENDDO                                                          
      ENDDO                                                             
  160 CONTINUE                                                          
C                                                                       
      IF(.NOT.LSYM19) GO TO 190                                         
      DO J=1,4                                                          
         DO I=1,6                                                       
            F(I,J,1,3)= F(I,J,4,2)                                      
            F(I,J,4,3)= F(I,J,2,2)                                      
            F(I,J,5,3)= F(I,J,6,2)                                      
         ENDDO                                                          
      ENDDO                                                             
  190 CONTINUE                                                          
C                                                                       
      RETURN                                                            
      END                                                               
C*MODULE INT2S   *DECK FIJKL5                                           
C>                                                                      
C>    @brief   auxiliary routine of order 5 for rot.axis integrations   
C>                                                                      
C>    @details auxiliary routine of order 5 for rot.axis integrations   
C>                                                                      
      SUBROUTINE FIJKL5(KX,LX,F,QX,QZ,E1,E2,E3,E4,E5)                   
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      DIMENSION F(6,6,6,6)                                              
C                                                                       
      DIMENSION       E1(  5,KX,LX),E2(4,3,KX,LX),E3(4,6,KX,LX)         
      PARAMETER (N=10)                                                  
      DIMENSION       E4(2,N,KX,LX),E5( 15,KX,LX)                       
C                                                                       
      LOGICAL   LSYM18,LSYM20,LSYM21                                    
C                                                                       
      PARAMETER (F02=2.0D+00)                                           
      PARAMETER (F03=3.0D+00)                                           
      PARAMETER (F04=4.0D+00)                                           
      PARAMETER (F06=6.0D+00)                                           
Cjms                                                                    
C     Simplified calculation of F(I,J,K,L) for cases                    
C     where:  I = 1..6,  J = 1..6,  K = 1..KX  and  L = 1..LX           
C     using auxiliary arrays E1, E2, E3, E4 and E5                      
Cjms                                                                    
      XX= QX*QX                                                         
      ZZ= QZ*QZ                                                         
      XZ= QX*QZ                                                         
      XXX= XX*QX                                                        
      XXZ= XX*QZ                                                        
      XZZ= ZZ*QX                                                        
      ZZZ= ZZ*QZ                                                        
      XXXX= XX*XX                                                       
      XXXZ= XX*XZ                                                       
      XXZZ= XX*ZZ                                                       
      XZZZ= XZ*ZZ                                                       
      ZZZZ= ZZ*ZZ                                                       
C                                                                       
      QXD= QX+QX                                                        
      QZD= QZ+QZ                                                        
C---  XXD= XX+XX                                                        
      XZD= XZ+XZ                                                        
C---  ZZD= ZZ+ZZ                                                        
      XXXD= XXX+XXX                                                     
      XXZD= XXZ+XXZ                                                     
      XZZD= XZZ+XZZ                                                     
      ZZZD= ZZZ+ZZZ                                                     
C                                                                       
      XZQ= XZ+XZ+XZ+XZ                                                  
C                                                                       
      LSYM18=(KX.EQ.4 .AND. LX.EQ.4)                                    
      LSYM20=(KX.EQ.6 .AND. LX.EQ.4)                                    
      LSYM21=(KX.EQ.6 .AND. LX.EQ.6)                                    
      DO L=1,LX                                                         
         DO K=1,KX                                                      
C                                                                       
            IF(LSYM18) THEN                                             
               IF(K.EQ.2 .AND. L.EQ.3) GO TO 001                        
            ENDIF                                                       
C                                                                       
            IF(LSYM20) THEN                                             
               IF(K.EQ.1 .AND. L.EQ.3) GO TO 001                        
               IF(K.EQ.4 .AND. L.EQ.3) GO TO 001                        
               IF(K.EQ.5 .AND. L.EQ.3) GO TO 001                        
            ENDIF                                                       
C                                                                       
            IF(LSYM21) THEN                                             
               IF(K.EQ.1 .AND. L.EQ.2) GO TO 001                        
C                                                                       
               IF(K.NE.3 .AND. L.EQ.4) GO TO 001                        
C                                                                       
               IF(K.EQ.1 .AND. L.EQ.6) GO TO 001                        
               IF(K.EQ.4 .AND. L.EQ.6) GO TO 001                        
               IF(K.EQ.5 .AND. L.EQ.6) GO TO 001                        
            ENDIF                                                       
C                                                                       
            F(1,1,K,L)=+E5(  1,K,L)+E3(1,1,K,L)*F06+E1(  1,K,L)*F03     
     *               +(+E4(1,1,K,L)+E4(2,1,K,L)+                        
     *                (+E2(1,1,K,L)+E2(2,1,K,L))*F03)*QXD               
     *               +(+E3(2,1,K,L)+E3(3,1,K,L)*F04+E3(4,1,K,L)         
     *                 +E1(  2,K,L)+E1(  3,K,L)*F04+E1(  4,K,L))*XX     
     *               +(+E2(3,1,K,L)+E2(4,1,K,L))*XXXD                   
     *                 +E1(  5,K,L)*XXXX                                
            F(2,1,K,L)=+E5(  4,K,L)+E3(1,4,K,L)+E3(1,1,K,L)+E1(  1,K,L) 
     *               +(+E4(2,4,K,L)+E2(2,1,K,L))*QXD                    
     *               +(+E3(4,4,K,L)+E1(  4,K,L))*XX                     
            F(3,1,K,L)=+E5(  6,K,L)+E3(1,6,K,L)+E3(1,1,K,L)+E1(  1,K,L) 
     *               +(+E4(2,6,K,L)+E2(2,1,K,L))*QXD                    
     *               +(+E4(1,3,K,L)+E2(1,3,K,L))*QZD                    
     *               +(+E3(4,6,K,L)+E1(  4,K,L))*XX                     
     *                 +E3(3,3,K,L)*XZQ                                 
     *               +(+E3(2,1,K,L)+E1(  2,K,L))*ZZ                     
     *                 +E2(4,3,K,L)*XXZD                                
     *                 +E2(3,1,K,L)*XZZD                                
     *                 +E1(  5,K,L)*XXZZ                                
            F(4,1,K,L)=+E5(  2,K,L)+E3(1,2,K,L)*F03                     
     *               +(+E4(1,2,K,L)+E4(2,2,K,L)*F02                     
     *                 +E2(1,2,K,L)+E2(2,2,K,L)*F02)*QX                 
     *               +(+E3(3,2,K,L)*F02+E3(4,2,K,L))*XX                 
     *                 +E2(4,2,K,L)*XXX                                 
            F(5,1,K,L)=+E5(  3,K,L)+E3(1,3,K,L)*F03                     
     *               +(+E4(1,3,K,L)+E4(2,3,K,L)*F02                     
     *                 +E2(1,3,K,L)+E2(2,3,K,L)*F02)*QX                 
     *               +(+E4(1,1,K,L)+E2(1,1,K,L)*F03)*QZ                 
     *               +(+E3(3,3,K,L)*F02+E3(4,3,K,L))*XX                 
     *               +(+E3(2,1,K,L)+E3(3,1,K,L)*F02                     
     *                 +E1(  2,K,L)+E1(  3,K,L)*F02)*XZ                 
     *                 +E2(4,3,K,L)*XXX                                 
     *               +(+E2(3,1,K,L)*F02+E2(4,1,K,L))*XXZ                
     *                 +E1(  5,K,L)*XXXZ                                
            F(6,1,K,L)=+E5(  5,K,L)+E3(1,5,K,L)                         
     *                 +E4(2,5,K,L)*QXD                                 
     *               +(+E4(1,2,K,L)+E2(1,2,K,L))*QZ                     
     *                 +E3(4,5,K,L)*XX                                  
     *                 +E3(3,2,K,L)*XZD                                 
     *                 +E2(4,2,K,L)*XXZ                                 
C                                                                       
            F(1,2,K,L)=+E5(  4,K,L)+E3(1,4,K,L)+E3(1,1,K,L)+E1(  1,K,L) 
     *               +(+E4(1,4,K,L)+E2(1,1,K,L))*QXD                    
     *               +(+E3(2,4,K,L)+E1(  2,K,L))*XX                     
            F(2,2,K,L)=+E5( 11,K,L)+E3(1,4,K,L)*F06+E1(  1,K,L)*F03     
            F(3,2,K,L)=+E5( 13,K,L)+E3(1,6,K,L)+E3(1,4,K,L)+E1(  1,K,L) 
     *               +(+E4(1,8,K,L)+E2(1,3,K,L))*QZD                    
     *               +(+E3(2,4,K,L)+E1(  2,K,L))*ZZ                     
            F(4,2,K,L)=+E5(  7,K,L)+E3(1,2,K,L)*F03                     
     *               +(+E4(1,7,K,L)+E2(1,2,K,L)*F03)*QX                 
            F(5,2,K,L)=+E5(  8,K,L)+E3(1,3,K,L)                         
     *               +(+E4(1,8,K,L)+E2(1,3,K,L))*QX                     
     *               +(+E4(1,4,K,L)+E2(1,1,K,L))*QZ                     
     *               +(+E3(2,4,K,L)+E1(  2,K,L))*XZ                     
            F(6,2,K,L)=+E5( 12,K,L)+E3(1,5,K,L)*F03                     
     *               +(+E4(1,7,K,L)+E2(1,2,K,L)*F03)*QZ                 
C                                                                       
            F(1,3,K,L)=+E5(  6,K,L)+E3(1,6,K,L)+E3(1,1,K,L)+E1(  1,K,L) 
     *               +(+E4(1,6,K,L)+E2(1,1,K,L))*QXD                    
     *               +(+E4(2,3,K,L)+E2(2,3,K,L))*QZD                    
     *               +(+E3(2,6,K,L)+E1(  2,K,L))*XX                     
     *                 +E3(3,3,K,L)*XZQ                                 
     *               +(+E3(4,1,K,L)+E1(  4,K,L))*ZZ                     
     *                 +E2(3,3,K,L)*XXZD                                
     *                 +E2(4,1,K,L)*XZZD                                
     *                 +E1(  5,K,L)*XXZZ                                
            F(2,3,K,L)=+E5( 13,K,L)+E3(1,6,K,L)+E3(1,4,K,L)+E1(  1,K,L) 
     *               +(+E4(2,8,K,L)+E2(2,3,K,L))*QZD                    
     *               +(+E3(4,4,K,L)+E1(  4,K,L))*ZZ                     
            F(3,3,K,L)=+E5( 15,K,L)+E3(1,6,K,L)*F06+E1(  1,K,L)*F03     
     *               +(+E4(1,N,K,L)+E4(2,N,K,L)+                        
     *                (+E2(1,3,K,L)+E2(2,3,K,L))*F03)*QZD               
     *               +(+E3(2,6,K,L)+E3(3,6,K,L)*F04+E3(4,6,K,L)         
     *                 +E1(  2,K,L)+E1(  3,K,L)*F04+E1(  4,K,L))*ZZ     
     *               +(+E2(3,3,K,L)+E2(4,3,K,L))*ZZZD                   
     *                 +E1(  5,K,L)*ZZZZ                                
            F(4,3,K,L)=+E5(  9,K,L)+E3(1,2,K,L)                         
     *               +(+E4(1,9,K,L)+E2(1,2,K,L))*QX                     
     *                 +E4(2,5,K,L)*QZD                                 
     *                 +E3(3,5,K,L)*XZD                                 
     *                 +E3(4,2,K,L)*ZZ                                  
     *                 +E2(4,2,K,L)*XZZ                                 
            F(5,3,K,L)=+E5( 10,K,L)+E3(1,3,K,L)*F03                     
     *               +(+E4(1,N,K,L)+E2(1,3,K,L)*F03)*QX                 
     *               +(+E4(1,6,K,L)+E4(2,6,K,L)*F02                     
     *                 +E2(1,1,K,L)+E2(2,1,K,L)*F02)*QZ                 
     *               +(+E3(2,6,K,L)+E3(3,6,K,L)*F02                     
     *                 +E1(  2,K,L)+E1(  3,K,L)*F02)*XZ                 
     *               +(+E3(3,3,K,L)*F02+E3(4,3,K,L))*ZZ                 
     *               +(+E2(3,3,K,L)*F02+E2(4,3,K,L))*XZZ                
     *                 +E2(4,1,K,L)*ZZZ                                 
     *                 +E1(  5,K,L)*XZZZ                                
            F(6,3,K,L)=+E5( 14,K,L)+E3(1,5,K,L)*F03                     
     *               +(+E4(1,9,K,L)+E4(2,9,K,L)*F02                     
     *                 +E2(1,2,K,L)+E2(2,2,K,L)*F02)*QZ                 
     *               +(+E3(3,5,K,L)*F02+E3(4,5,K,L))*ZZ                 
     *                 +E2(4,2,K,L)*ZZZ                                 
C                                                                       
            F(1,4,K,L)=+E5(  2,K,L)+E3(1,2,K,L)*F03                     
     *               +(+E4(1,2,K,L)*F02+E4(2,2,K,L)                     
     *                 +E2(1,2,K,L)*F02+E2(2,2,K,L))*QX                 
     *               +(+E3(2,2,K,L)+E3(3,2,K,L)*F02)*XX                 
     *                 +E2(3,2,K,L)*XXX                                 
            F(2,4,K,L)=+E5(  7,K,L)+E3(1,2,K,L)*F03                     
     *               +(+E4(2,7,K,L)+E2(2,2,K,L)*F03)*QX                 
            F(3,4,K,L)=+E5(  9,K,L)+E3(1,2,K,L)                         
     *               +(+E4(2,9,K,L)+E2(2,2,K,L))*QX                     
     *                 +E4(1,5,K,L)*QZD                                 
     *                 +E3(3,5,K,L)*XZD                                 
     *                 +E3(2,2,K,L)*ZZ                                  
     *                 +E2(3,2,K,L)*XZZ                                 
            F(4,4,K,L)=+E5(  4,K,L)+E3(1,4,K,L)+E3(1,1,K,L)+E1(  1,K,L) 
     *               +(+E4(1,4,K,L)+E4(2,4,K,L)                         
     *                 +E2(1,1,K,L)+E2(2,1,K,L))*QX                     
     *               +(+E3(3,4,K,L)+E1(  3,K,L))*XX                     
            F(5,4,K,L)=+E5(  5,K,L)+E3(1,5,K,L)                         
     *               +(+E4(1,5,K,L)+E4(2,5,K,L))*QX                     
     *               +(+E4(1,2,K,L)+E2(1,2,K,L))*QZ                     
     *                 +E3(3,5,K,L)*XX                                  
     *               +(+E3(2,2,K,L)+E3(3,2,K,L))*XZ                     
     *                 +E2(3,2,K,L)*XXZ                                 
            F(6,4,K,L)=+E5(  8,K,L)+E3(1,3,K,L)                         
     *               +(+E4(2,8,K,L)+E2(2,3,K,L))*QX                     
     *               +(+E4(1,4,K,L)+E2(1,1,K,L))*QZ                     
     *               +(+E3(3,4,K,L)+E1(  3,K,L))*XZ                     
C                                                                       
            F(1,5,K,L)=+E5(  3,K,L)+E3(1,3,K,L)*F03                     
     *               +(+E4(1,3,K,L)*F02+E4(2,3,K,L)                     
     *                 +E2(1,3,K,L)*F02+E2(2,3,K,L))*QX                 
     *               +(+E4(2,1,K,L)+E2(2,1,K,L)*F03)*QZ                 
     *               +(+E3(2,3,K,L)+E3(3,3,K,L)*F02)*XX                 
     *               +(+E3(3,1,K,L)*F02+E3(4,1,K,L)                     
     *                 +E1(  3,K,L)*F02+E1(  4,K,L))*XZ                 
     *                 +E2(3,3,K,L)*XXX                                 
     *               +(+E2(3,1,K,L)+E2(4,1,K,L)*F02)*XXZ                
     *                 +E1(  5,K,L)*XXXZ                                
            F(2,5,K,L)=+E5(  8,K,L)+E3(1,3,K,L)                         
     *               +(+E4(2,8,K,L)+E2(2,3,K,L))*QX                     
     *               +(+E4(2,4,K,L)+E2(2,1,K,L))*QZ                     
     *               +(+E3(4,4,K,L)+E1(  4,K,L))*XZ                     
            F(3,5,K,L)=+E5( 10,K,L)+E3(1,3,K,L)*F03                     
     *               +(+E4(2,N,K,L)+E2(2,3,K,L)*F03)*QX                 
     *               +(+E4(1,6,K,L)*F02+E4(2,6,K,L)                     
     *                 +E2(1,1,K,L)*F02+E2(2,1,K,L))*QZ                 
     *               +(+E3(3,6,K,L)*F02+E3(4,6,K,L)                     
     *                 +E1(  3,K,L)*F02+E1(  4,K,L))*XZ                 
     *               +(+E3(2,3,K,L)+E3(3,3,K,L)*F02)*ZZ                 
     *               +(+E2(3,3,K,L)+E2(4,3,K,L)*F02)*XZZ                
     *                 +E2(3,1,K,L)*ZZZ                                 
     *                 +E1(  5,K,L)*XZZZ                                
            F(4,5,K,L)=+E5(  5,K,L)+E3(1,5,K,L)                         
     *               +(+E4(1,5,K,L)+E4(2,5,K,L))*QX                     
     *               +(+E4(2,2,K,L)+E2(2,2,K,L))*QZ                     
     *                 +E3(3,5,K,L)*XX                                  
     *               +(+E3(3,2,K,L)+E3(4,2,K,L))*XZ                     
     *                 +E2(4,2,K,L)*XXZ                                 
            F(5,5,K,L)=+E5(  6,K,L)+E3(1,6,K,L)+E3(1,1,K,L)+E1(  1,K,L) 
     *               +(+E4(1,6,K,L)+E4(2,6,K,L)                         
     *                 +E2(1,1,K,L)+E2(2,1,K,L))*QX                     
     *               +(+E4(1,3,K,L)+E4(2,3,K,L)                         
     *                 +E2(1,3,K,L)+E2(2,3,K,L))*QZ                     
     *               +(+E3(3,6,K,L)+E1(  3,K,L))*XX                     
     *               +(+E3(2,3,K,L)+E3(3,3,K,L)*F02+E3(4,3,K,L))*XZ     
     *               +(+E3(3,1,K,L)+E1(  3,K,L))*ZZ                     
     *               +(+E2(3,3,K,L)+E2(4,3,K,L))*XXZ                    
     *               +(+E2(3,1,K,L)+E2(4,1,K,L))*XZZ                    
     *                 +E1(  5,K,L)*XXZZ                                
            F(6,5,K,L)=+E5(  9,K,L)+E3(1,2,K,L)                         
     *               +(+E4(2,9,K,L)+E2(2,2,K,L))*QX                     
     *               +(+E4(1,5,K,L)+E4(2,5,K,L))*QZ                     
     *               +(+E3(3,5,K,L)+E3(4,5,K,L))*XZ                     
     *                 +E3(3,2,K,L)*ZZ                                  
     *                 +E2(4,2,K,L)*XZZ                                 
C                                                                       
            F(1,6,K,L)=+E5(  5,K,L)+E3(1,5,K,L)                         
     *                 +E4(1,5,K,L)*QXD                                 
     *               +(+E4(2,2,K,L)+E2(2,2,K,L))*QZ                     
     *                 +E3(2,5,K,L)*XX                                  
     *                 +E3(3,2,K,L)*XZD                                 
     *                 +E2(3,2,K,L)*XXZ                                 
            F(2,6,K,L)=+E5( 12,K,L)+E3(1,5,K,L)*F03                     
     *               +(+E4(2,7,K,L)+E2(2,2,K,L)*F03)*QZ                 
            F(3,6,K,L)=+E5( 14,K,L)+E3(1,5,K,L)*F03                     
     *               +(+E4(1,9,K,L)*F02+E4(2,9,K,L)                     
     *                 +E2(1,2,K,L)*F02+E2(2,2,K,L))*QZ                 
     *               +(+E3(2,5,K,L)+E3(3,5,K,L)*F02)*ZZ                 
     *                 +E2(3,2,K,L)*ZZZ                                 
            F(4,6,K,L)=+E5(  8,K,L)+E3(1,3,K,L)                         
     *               +(+E4(1,8,K,L)+E2(1,3,K,L))*QX                     
     *               +(+E4(2,4,K,L)+E2(2,1,K,L))*QZ                     
     *               +(+E3(3,4,K,L)+E1(  3,K,L))*XZ                     
            F(5,6,K,L)=+E5(  9,K,L)+E3(1,2,K,L)                         
     *               +(+E4(1,9,K,L)+E2(1,2,K,L))*QX                     
     *               +(+E4(1,5,K,L)+E4(2,5,K,L))*QZ                     
     *               +(+E3(2,5,K,L)+E3(3,5,K,L))*XZ                     
     *                 +E3(3,2,K,L)*ZZ                                  
     *                 +E2(3,2,K,L)*XZZ                                 
            F(6,6,K,L)=+E5( 13,K,L)+E3(1,6,K,L)+E3(1,4,K,L)+E1(  1,K,L) 
     *               +(+E4(1,8,K,L)+E4(2,8,K,L)                         
     *                 +E2(1,3,K,L)+E2(2,3,K,L))*QZ                     
     *               +(+E3(3,4,K,L)+E1(  3,K,L))*ZZ                     
C                                                                       
  001       CONTINUE                                                    
         ENDDO                                                          
      ENDDO                                                             
C                                                                       
      IF(.NOT.LSYM18) GO TO 180                                         
      DO J=1,6                                                          
         DO I=1,6                                                       
            F(I,J,2,3)= F(I,J,3,2)                                      
         ENDDO                                                          
      ENDDO                                                             
  180 CONTINUE                                                          
C                                                                       
      IF(.NOT.LSYM20) GO TO 200                                         
      DO J=1,6                                                          
         DO I=1,6                                                       
            F(I,J,1,3)= F(I,J,4,2)                                      
            F(I,J,4,3)= F(I,J,2,2)                                      
            F(I,J,5,3)= F(I,J,6,2)                                      
         ENDDO                                                          
      ENDDO                                                             
  200 CONTINUE                                                          
C                                                                       
      IF(.NOT.LSYM21) GO TO 210                                         
      DO J=1,6                                                          
         DO I=1,6                                                       
            F(I,J,1,2)= F(I,J,2,1)                                      
C                                                                       
            F(I,J,1,4)= F(I,J,4,1)                                      
            F(I,J,2,4)= F(I,J,4,2)                                      
            F(I,J,4,4)= F(I,J,2,1)                                      
            F(I,J,5,4)= F(I,J,6,1)                                      
            F(I,J,6,4)= F(I,J,5,2)                                      
C                                                                       
            F(I,J,1,6)= F(I,J,4,5)                                      
            F(I,J,4,6)= F(I,J,2,5)                                      
            F(I,J,5,6)= F(I,J,6,5)                                      
         ENDDO                                                          
      ENDDO                                                             
  210 CONTINUE                                                          
C                                                                       
      RETURN                                                            
      END                                                               
