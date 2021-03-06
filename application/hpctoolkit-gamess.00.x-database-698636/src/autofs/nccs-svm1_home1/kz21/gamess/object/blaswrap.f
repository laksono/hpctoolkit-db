C  5 APR 10 - MWS - WRAPPERS AROUND SELECTED LEVEL 1,2,3 BLAS CALLS     
C                                                                       
C  THIS IS MEANT TO BE USED ON SYSTEMS WHERE 64 BIT INTEGERS ARE        
C  ALLOWED BY THE FORTRAN COMPILER FLAGS, BUT WHERE THE VENDOR'S        
C  MATHEMATICAL LIBRARIES INSIST UPON HAVING 32 BIT INTEGER ARGS.       
C                                                                       
C  TO USE THIS, ONE MUST HACK ALL CALLS TO -D- BLAS ROUTINES INTO       
C  THE -X- NAMES USED BELOW.  SEE THE -COMP- SCRIPT FOR DETAILS.        
C                                                                       
      DOUBLE PRECISION FUNCTION XASUM(N,DX,INCX)                        
      IMPLICIT NONE                                                     
      DOUBLE PRECISION DX(*),DASUM                                      
      INTEGER   N  ,INCX                                                
      INTEGER*4 N32,INCX32                                              
      N32    = N                                                        
      INCX32 = INCX                                                     
      XASUM = DASUM(N32,DX,INCX32)                                      
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE XAXPY(N,DA,DX,INCX,DY,INCY)                            
      IMPLICIT NONE                                                     
      DOUBLE PRECISION DA,DX(*),DY(*)                                   
      INTEGER   N  ,INCX  ,INCY                                         
      INTEGER*4 N32,INCX32,INCY32                                       
      N32    = N                                                        
      INCX32 = INCX                                                     
      INCY32 = INCY                                                     
      CALL DAXPY(N32,DA,DX,INCX32,DY,INCY32)                            
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE XCOPY(N,DX,INCX,DY,INCY)                               
      IMPLICIT NONE                                                     
      DOUBLE PRECISION DX(*),DY(*)                                      
      INTEGER   N  ,INCX  ,INCY                                         
      INTEGER*4 N32,INCX32,INCY32                                       
      N32    = N                                                        
      INCX32 = INCX                                                     
      INCY32 = INCY                                                     
      CALL DCOPY(N32,DX,INCX32,DY,INCY32)                               
      RETURN                                                            
      END                                                               
C                                                                       
      DOUBLE PRECISION FUNCTION XDOT(N,DX,INCX,DY,INCY)                 
      IMPLICIT NONE                                                     
      DOUBLE PRECISION DX(*),DY(*),DDOT                                 
      INTEGER   N  ,INCX  ,INCY                                         
      INTEGER*4 N32,INCX32,INCY32                                       
      N32    = N                                                        
      INCX32 = INCX                                                     
      INCY32 = INCY                                                     
      XDOT = DDOT(N32,DX,INCX32,DY,INCY32)                              
      RETURN                                                            
      END                                                               
C                                                                       
      DOUBLE PRECISION FUNCTION XNRM2(N,DX,INCX)                        
      IMPLICIT NONE                                                     
      DOUBLE PRECISION DX(*),DNRM2                                      
      INTEGER   N  ,INCX                                                
      INTEGER*4 N32,INCX32                                              
      N32    = N                                                        
      INCX32 = INCX                                                     
      XNRM2 = DNRM2(N32,DX,INCX32)                                      
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE XROT(N,DX,INCX,DY,INCY,C,S)                            
      IMPLICIT NONE                                                     
      DOUBLE PRECISION DX(*),DY(*),C,S                                  
      INTEGER   N  ,INCX  ,INCY                                         
      INTEGER*4 N32,INCX32,INCY32                                       
      N32    = N                                                        
      INCX32 = INCX                                                     
      INCY32 = INCY                                                     
      CALL DROT(N32,DX,INCX32,DY,INCY32,C,S)                            
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE XROTG(DA,DB,C,S)                                       
      IMPLICIT NONE                                                     
      DOUBLE PRECISION DA,DB,C,S                                        
      CALL DROTG(DA,DB,C,S)                                             
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE XSCAL(N,DA,DX,INCX)                                    
      IMPLICIT NONE                                                     
      DOUBLE PRECISION DX(*),DA                                         
      INTEGER   N  ,INCX                                                
      INTEGER*4 N32,INCX32                                              
      N32    = N                                                        
      INCX32 = INCX                                                     
      CALL DSCAL(N32,DA,DX,INCX32)                                      
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE XSWAP(N,DX,INCX,DY,INCY)                               
      IMPLICIT NONE                                                     
      DOUBLE PRECISION DX(*),DY(*)                                      
      INTEGER   N  ,INCX  ,INCY                                         
      INTEGER*4 N32,INCX32,INCY32                                       
      N32    = N                                                        
      INCX32 = INCX                                                     
      INCY32 = INCY                                                     
      CALL DSWAP(N32,DX,INCX32,DY,INCY32)                               
      RETURN                                                            
      END                                                               
C                                                                       
      INTEGER FUNCTION IXAMAX(N,DX,INCX)                                
      IMPLICIT NONE                                                     
      INTEGER IDAMAX                                                    
      DOUBLE PRECISION DX(*)                                            
      INTEGER   N  ,INCX                                                
      INTEGER*4 N32,INCX32                                              
      N32    = N                                                        
      INCX32 = INCX                                                     
      IXAMAX = IDAMAX(N32,DX,INCX32)                                    
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE XGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)                    
      IMPLICIT NONE                                                     
      INTEGER   M  ,N  ,INCX  ,INCY  ,LDA                               
      INTEGER*4 M32,N32,INCX32,INCY32,LDA32                             
      DOUBLE PRECISION ALPHA,X(*),Y(*),A(LDA,*)                         
      M32    = M                                                        
      N32    = N                                                        
      INCX32 = INCX                                                     
      INCY32 = INCY                                                     
      LDA32  = LDA                                                      
      CALL DGER(M32,N32,ALPHA,X,INCX32,Y,INCY32,A,LDA32)                
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE XTRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)                  
      IMPLICIT NONE                                                     
      CHARACTER*1 UPLO,TRANS,DIAG                                       
      INTEGER   N  ,LDA  ,INCX                                          
      INTEGER*4 N32,LDA32,INCX32                                        
      DOUBLE PRECISION A(LDA,*),X(*)                                    
      N32    = N                                                        
      LDA32  = LDA                                                      
      INCX32 = INCX                                                     
      CALL DTRMV(UPLO,TRANS,DIAG,N32,A,LDA32,X,INCX32)                  
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE XGEMV(FORMA,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)        
      IMPLICIT NONE                                                     
      CHARACTER*1 FORMA                                                 
      DOUBLE PRECISION A(LDA,*),X(*),Y(*),ALPHA,BETA                    
      INTEGER   M,  N,  INCX,  INCY,  LDA                               
      INTEGER*4 M32,N32,INCX32,INCY32,LDA32                             
      M32    = M                                                        
      N32    = N                                                        
      INCX32 = INCX                                                     
      INCY32 = INCY                                                     
      LDA32  = LDA                                                      
      CALL DGEMV(FORMA,M32,N32,ALPHA,A,LDA32,X,INCX32,BETA,Y,INCY32)    
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE XSPMV(FORMU,N,ALPHA,A,X,INCX,BETA,Y,INCY)              
      IMPLICIT NONE                                                     
      CHARACTER*1 FORMU                                                 
      DOUBLE PRECISION A(*),X(*),Y(*),ALPHA,BETA                        
      INTEGER   N,  INCX,  INCY                                         
      INTEGER*4 N32,INCX32,INCY32                                       
      N32     = N                                                       
      INCX32  = INCX                                                    
      INCY32  = INCY                                                    
      CALL DSPMV(FORMU,N32,ALPHA,A,X,INCX32,BETA,Y,INCY32)              
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE XSPR2(UPLO,N,ALPHA,X,INCX,Y,INCY,AP)                   
      IMPLICIT NONE                                                     
      CHARACTER*1 UPLO                                                  
      DOUBLE PRECISION ALPHA,X(*),Y(*),AP(*)                            
      INTEGER   N,  INCX,  INCY                                         
      INTEGER*4 N32,INCX32,INCY32                                       
      N32     = N                                                       
      INCX32  = INCX                                                    
      INCY32  = INCY                                                    
      CALL DSPR2(UPLO,N32,ALPHA,X,INCX32,Y,INCY32,AP)                   
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE XGEMM(FORMA,FORMB,L,N,M,ALPHA,A,LDA,B,LDB,BETA,C,LDC)  
      IMPLICIT NONE                                                     
      CHARACTER*1 FORMA,FORMB                                           
      DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*),ALPHA,BETA            
      INTEGER   L,  M,  N,  LDA,  LDB,  LDC                             
      INTEGER*4 L32,M32,N32,LDA32,LDB32,LDC32                           
      L32    = L                                                        
      M32    = M                                                        
      N32    = N                                                        
      LDA32  = LDA                                                      
      LDB32  = LDB                                                      
      LDC32  = LDC                                                      
      CALL DGEMM(FORMA,FORMB,L32,N32,M32,ALPHA,A,LDA32,B,LDB32,BETA,    
     *           C,LDC32)                                               
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE XTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)     
      IMPLICIT NONE                                                     
      CHARACTER*1 SIDE,UPLO,TRANSA,DIAG                                 
      INTEGER   M  ,N  ,LDA  ,LDB                                       
      INTEGER*4 M32,N32,LDA32,LDB32                                     
      DOUBLE PRECISION ALPHA,A(LDA,*),B(LDB,*)                          
      M32   = M                                                         
      N32   = N                                                         
      LDA32 = LDA                                                       
      LDB32 = LDB                                                       
      CALL DTRMM(SIDE,UPLO,TRANSA,DIAG,M32,N32,ALPHA,A,LDA32,B,LDB32)   
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE XTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)     
      CHARACTER*1 SIDE,UPLO,TRANSA,DIAG                                 
      INTEGER   M  ,N  ,LDA  ,LDB                                       
      INTEGER*4 M32,N32,LDA32,LDB32                                     
      DOUBLE PRECISION ALPHA,A(LDA,*),B(LDB,*)                          
      M32   = M                                                         
      N32   = N                                                         
      LDA32 = LDA                                                       
      LDB32 = LDB                                                       
      CALL DTRSM(SIDE,UPLO,TRANSA,DIAG,M32,N32,ALPHA,A,LDA32,B,LDB32)   
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE XSYRK(UPLO, TRANS, N, K, ALPHA, A, LDA, BETA, C, LDC)  
      CHARACTER*1 UPLO,TRANS                                            
      INTEGER   N  ,K  ,LDA  ,LDC                                       
      INTEGER*4 N32,K32,LDA32,LDC32                                     
      DOUBLE PRECISION ALPHA,A(LDA,*),C(LDC,*)                          
      N32   = N                                                         
      K32   = K                                                         
      LDA32 = LDA                                                       
      LDC32 = LDC                                                       
      CALL DSYRK(UPLO, TRANS, N32, K32, ALPHA, A, LDA32, BETA, C, LDC32)
      RETURN                                                            
      END                                                               
