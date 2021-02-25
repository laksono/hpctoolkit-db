      SUBROUTINE genr70_switch_gpu &
                (ISHELL,JSHELL,KSHELL,LSHELL, &
                 ITYPE,LAT,LBT,LCT,LDT, & 
                 LA,LB,LC,LD, &
                 INEW,JNEW,KNEW,LNEW, &
                 LPOPI,LPOPJ,LPOPK,LPOPL)
      use mx_limits, only: mxgtot,mxsh,mxgsh,mxg2

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INTEGER :: LA,LB,LC,LD
      INTEGER :: INEW,JNEW,KNEW,LNEW



      IF(ITYPE.EQ. 1 .OR. ITYPE.EQ. 2 .OR. ITYPE.EQ. 4 .OR. &
         ITYPE.EQ. 6 .OR. ITYPE.EQ. 8 .OR. ITYPE.EQ.16) THEN
!
! TYPES 0000,0001,0011,0101,0111,1111 ARE UNALTERED
!
         INEW= ISHELL
         JNEW= JSHELL
         KNEW= KSHELL
         LNEW= LSHELL
         LA= LAT
         LB= LBT
         LC= LCT
         LD= LDT

         LPOPI = 64
         LPOPJ = 16
         LPOPK =  4
         LPOPL =  1
      ELSEIF(ITYPE.EQ.10 .OR. ITYPE.EQ.12) THEN
!
! TYPES 1001,1011 HAVE IJ SWITCHED
!
         INEW= JSHELL
         JNEW= ISHELL
         KNEW= KSHELL
         LNEW= LSHELL
         LA= LBT
         LB= LAT
         LC= LCT
         LD= LDT

         LPOPI = 16
         LPOPJ = 64
         LPOPK =  4
         LPOPL =  1
      ELSEIF(ITYPE.EQ. 3 .OR. ITYPE.EQ. 7) THEN
!
! TYPES 0010,0110 HAVE KL SWITCHED
!
         INEW= ISHELL
         JNEW= JSHELL
         KNEW= LSHELL
         LNEW= KSHELL
         LA= LAT
         LB= LBT
         LC= LDT
         LD= LCT

         LPOPI = 64
         LPOPJ = 16
         LPOPK =  1
         LPOPL =  4
      ELSEIF(ITYPE.EQ. 5 .OR. ITYPE.EQ.13 .OR. ITYPE.EQ.14) THEN
!
! TYPES 0100,1100,1101 HAVE PAIRS IJ AND KL SWITCHED
!
         INEW= KSHELL
         JNEW= LSHELL
         KNEW= ISHELL
         LNEW= JSHELL
         LA= LCT
         LB= LDT
         LC= LAT
         LD= LBT

         LPOPI =  4
         LPOPJ =  1
         LPOPK = 64
         LPOPL = 16
      ELSEIF(ITYPE.EQ.11) THEN
!
! TYPE 1010 HAS IJ SWITCHED AND KL SWITCHED
!
         INEW= JSHELL
         JNEW= ISHELL
         KNEW= LSHELL
         LNEW= KSHELL
         LA= LBT
         LB= LAT
         LC= LDT
         LD= LCT

         LPOPI = 16
         LPOPJ = 64
         LPOPK =  1
         LPOPL =  4
      ELSEIF(ITYPE.EQ. 9) THEN
!
! TYPE 1000  HAS PAIRS IJ AND KL SWITCHED FOLLOWED BY KL SWITCH
!
         INEW= KSHELL
         JNEW= LSHELL
         KNEW= JSHELL
         LNEW= ISHELL
         LA= LCT
         LB= LDT
         LC= LBT
         LD= LAT

         LPOPI =  1
         LPOPJ =  4
         LPOPK = 64
         LPOPL = 16
      ELSEIF(ITYPE.EQ.15) THEN
!
! TYPE 1110 HAS PAIRS IJ AND KL SWITCHED FOLLOWED BY IJ SWITCH
!
         INEW= LSHELL
         JNEW= KSHELL
         KNEW= ISHELL
         LNEW= JSHELL
         LA= LDT
         LB= LCT
         LC= LAT
         LD= LBT

         LPOPI =  4
         LPOPJ =  1
         LPOPK = 16
         LPOPL = 64
      ENDIF
      
      END




!*MODULE INT2B   *DECK GENR70
!>
!>    @brief   rotated axis integration involving s,p,L shells
!>
!>    @details rotated axis integration involving s,p,L shells,
!>             by a mix of rotated axis manipulations,
!>                        J.A.Pople, W.J.Hehre
!>                    J.Comput.Phys. 27, 161-168(1978)
!>             and McMurchie/Davidson quadrature,
!>                       L.E.McMurchie, E.R.Davidson
!>                    J.Comput.Phys. 26, 218-231(1978)
!>             according to
!>                       K.Ishimura, S.Nagase
!>                Theoret.Chem.Acc. 120, 185-189(2008)
!>             The accuracy parameters were adjusted in 2004 to
!>             make this program as good as Rys quadrature.
!>             Outline:
!>    o  GENR70 EVALUATES INTEGRALS INVOLVING ONLY S AND P FUNCTIONS.
!>       ONLY POSSIBILITIES ALLOWED FOR ANGULAR QUANTUM NUMBERS ARE
!>                0000  0001  0011  0101  0111  1111
!>       DETERMINES TYPE OF INTEGRAL SET BASED ON THE ABOVE NUMBERS
!>       THIS DRIVER CALLS THE FOLLOWING ROUTINES IN THE ORDER GIVEN
!>    o  SINFO, TO PRESET INTEGRAL ACCURACY LIMITS
!>    o  SGEOM, OBTAINS GEOMETRICAL INFORMATION ABOUT THE FOUR CENTERS
!>              AND FINDS TWO SETS OF LOCAL AXES FOR CENTERS:
!>              (A AND B) P SET AND (C AND D) Q SET
!>    o  PINF, OBTAINS INFORMATION ABOUT GAUSSIAN FUNCTIONS CONNECTED
!>             WITH THE P SET OF AXES AT THIS POINT GENR70 OBTAINS
!>             INFORMATION ABOUT THE GAUSSIAN FUNCTIONS CONNECTED WITH
!>             THE Q SET OF AXES
!>    o  SP0000 TO SP1111, OBTAIN UP TO 88 INTEGRALS REFERRED
!>                         TO AXES A B AND Q
!>    o  ROT2, ROTATES THESE INTEGRALS TO UP TO 160 INTEGRALS ON A B AND Q
!>    o  TQ0011 TO TQ1111, TRANSLATE THESE INTEGRALS ON A B AND Q
!>              TO UP TO 256 INTEGRALS ON A B C AND D
!>    o  R30001 TO R31111, ROTATE UP TO 256 INTEGRALS ON A B C AND D
!>              TO THE SAME NUMBER REFERRED TO THE FIXED SPACE AXES
!>
!>    @author  Warren Hehre circa Gaussian70,
!>             Jose Sierra (Synstar, SA, Madrid) revised loop structure
!>             and increased the Fm(t) grid density accuracy in 2001.
!>             Kazuya Ishimura (IMS) replaced inner axis rotation by
!>             direct McMurchie-Davidson quadrature in 2004.
!>             fully accurate Fm(t) and cutoff changes in 2004.
!>




      SUBROUTINE genr70_pre_gpu &
                (JTYPE, &
                 ISHELL,JSHELL,KSHELL,LSHELL, &
                 INEW,JNEW,KNEW,LNEW, &
                 P12,P34,R34,ACX,ACZ, &                 
                 RCD,SING,COSG,P,T)
      use mx_limits, only: mxgtot,mxsh,mxgsh,mxg2

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      double precision,parameter :: ZER=0.0D+00
      double precision,parameter :: PT5=0.5D+00
      double precision,parameter :: ONE=1.0D+00
      double precision,parameter :: PT7=0.7D+00
      double precision,parameter :: PT9=0.9D+00

!        ACYCUT CHANGED FROM A VALUE OF 1E-4 IN 2004,
!        FOR ACCURACY WHEN USING DIFFUSE EXPONENTS.
      double precision,parameter :: ACYCUT=1.0D-10
      double precision,parameter :: TENM12=1.0D-12

      ! PITO52=(PI+PI)*PI*SQRT(PI )
      double precision,parameter :: PITO52=34.986836655249726D+00

      COMMON /B     / CO(MXSH,3)
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS, &
                      TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2), &
                      D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2), &
                      NGANGB
!$omp threadprivate(/GEOMPQ/)

      COMMON /INTAC2/ EI1,EI2,CUX
      COMMON /MAXC  / CMAX(MXGTOT),CMAXA(MXGSH),CMAXB(MXGSH), &
                      CMAXC(MXGSH),CMAXD(MXGSH),ISMLP(MXG2),ISMLQ
!$omp threadprivate(/MAXC/)                      
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT), &
                      CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT), &
                      KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH), &
                      KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /POPOUT/ LPOPI,LPOPJ,LPOPK,LPOPL
!$omp threadprivate(/POPOUT/)      

      INTEGER :: INEW,JNEW,KNEW,LNEW

      COMMON /SHLLFO/ NGA,LA,EXA(MXGSH),CSA(MXGSH),CPA(MXGSH), &
                      NGB,LB,EXB(MXGSH),CSB(MXGSH),CPB(MXGSH), &
                      NGC,LC,EXC(MXGSH),CSC(MXGSH),CPC(MXGSH), &
                      NGD,LD,EXD(MXGSH),CSD(MXGSH),CPD(MXGSH)
!$omp threadprivate(/SHLLFO/)                      
      COMMON /SHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL, &
                      MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL, &
                      NIJ,IJ,KL,IJKL
!$omp threadprivate(/SHLNOS/)                      
!
!     ==================================================================
!
!JMS  LABELLED COMMON JMSGYH DEFINED FOR COMPUTATIONAL EFFICIENCY.
!JMS  IT IS ONLY USED IN THIS MODULE INT2B AND IN MODULE INT2R.
!
      COMMON /JMSGYH/ SQ(0:1,0:1)
!$omp threadprivate(/JMSGYH/)
      COMMON /KI2 / ACY,ACY2,AQX,AQX2,AQXY,Y03,Y04
!$omp threadprivate(/KI2/)      

      double precision :: P12(3,3),P34(3,3),P(3,3),T(3)



      PIF= PITO52


!
! NUMBERS OF GAUSSIAN FUNCTIONS IN SHELLS INEW JNEW KNEW AND LNEW
!
      NGA= KNG(INEW)
      NGB= KNG(JNEW)
      NGC= KNG(KNEW)
      NGD= KNG(LNEW)
!
! STARTING LOCATIONS OF SHELLS INEW JNEW KNEW AND LNEW IN LIST
! OF GAUSSIAN FUNCTIONS
!
      I= KSTART(INEW)-1
      J= KSTART(JNEW)-1
      K= KSTART(KNEW)-1
      L= KSTART(LNEW)-1
!
! LOOP OVER GAUSSIANS IN EACH SHELL
! FIRST SHELL INEW
!
      DO NI=1,NGA
         N=I+NI
         ! THE MAXIMUM COEFFICIENT ASSOCIATED WITH SHELL IS
         ! USED TO DETERMINE IF ANY OF THE INTEGRALS ASSOCIATED WITH A SET
         ! OF SHELLS IS LARGE ENOUGH TO WARRANT EVALUATION OF THE ENTIRE SET
         CMAXA(NI)= CMAX(N)
         ! GAUSSIAN EXPONENTS
         EXA(NI)= EX(N)
         ! S COEFFICIENTS
         CSA(NI)= CS(N)
         ! P COEFFICIENTS
         CPA(NI)= CP(N)
      ENDDO
!
! REPEAT PROCEDURE FOR SHELLS JNEW KNEW AND LNEW
!
      DO NJ=1,NGB
         N=J+NJ
         CMAXB(NJ)= CMAX(N)
         EXB(NJ)= EX(N)
         CSB(NJ)= CS(N)
         CPB(NJ)= CP(N)
      ENDDO

      DO NK=1,NGC
         N=K+NK
         CMAXC(NK)= CMAX(N)*QQ4
         EXC(NK)= EX(N)
         CSC(NK)= CS(N)*QQ4
         CPC(NK)= CP(N)*QQ4
      ENDDO

      DO NL=1,NGD
         N=L+NL
         CMAXD(NL)= CMAX(N)
         EXD(NL)= EX(N)
         CSD(NL)= CS(N)
         CPD(NL)= CP(N)
      ENDDO 

      NGANGB=NGA*NGB
!
! COORDINATES OF ATOMS ASSOCIATED WITH SHELLS INEW JNEW KNEW AND LNEW
!
      R12= ZER
      R34= ZER
      DO N=1,3
         P12(N,1)= CO(INEW,N)
         P12(N,2)= CO(JNEW,N)
         P12(N,3)= P12(N,2)-P12(N,1)
         
         R12= R12+P12(N,3)*P12(N,3)

         P34(N,1)= CO(KNEW,N)
         P34(N,2)= CO(LNEW,N)
         P34(N,3)= P34(N,2)-P34(N,1)
      
         R34= R34+P34(N,3)*P34(N,3)

      ENDDO

!
! FIND DIRECTION COSINES OF PENULTIMATE AXES FROM COORDINATES OF AB
! P(1,1),P(1,2),... ARE DIRECTION COSINES OF AXES AT P.  Z-AXIS ALONG AB
! T(1),T(2),T(3)... ARE DIRECTION COSINES OF AXES AT Q.  Z-AXIS ALONG CD
!
! FIND DIRECTION COSINES OF AB AND CD. THESE ARE LOCAL Z-AXES.
! IF INDETERMINATE TAKE ALONG SPACE Z-AXIS
!
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
!
! FIND LOCAL Y-AXIS AS COMMON PERPENDICULAR TO AB AND CD
! IF INDETERMINATE TAKE PERPENDICULAR TO AB AND SPACE Z-AXIS
! IF STILL INDETERMINATE TAKE PERPENDICULAR TO AB AND SPACE X-AXIS
!
      COSG= T(1)*P(1,3)+T(2)*P(2,3)+T(3)*P(3,3)
      COSG= MIN( ONE,COSG)
      COSG= MAX(-ONE,COSG)
!
! MODIFIED ROTATION TESTING.
! THIS FIX CURES THE SMALL ANGLE PROBLEM.
!
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
!
! FIND DIRECTION COSINES OF LOCAL X-AXES
!
      P(1,1)= P(2,2)*P(3,3)-P(3,2)*P(2,3)
      P(2,1)= P(3,2)*P(1,3)-P(1,2)*P(3,3)
      P(3,1)= P(1,2)*P(2,3)-P(2,2)*P(1,3)
!
! FIND COORDINATES OF C RELATIVE TO LOCAL AXES AT A
!
      T(1)= P34(1,1)-P12(1,1)
      T(2)= P34(2,1)-P12(2,1)
      T(3)= P34(3,1)-P12(3,1)
      ACX = T(1)*P(1,1)+T(2)*P(2,1)+T(3)*P(3,1)
      ACY = T(1)*P(1,2)+T(2)*P(2,2)+T(3)*P(3,2)
      ACZ = T(1)*P(1,3)+T(2)*P(2,3)+T(3)*P(3,3)
!
! SET ACY= 0  IF CLOSE
!
      IF( ABS(ACY).LE.ACYCUT) THEN
         ACY = ZER
         ACY2= ZER
      ELSE
         ACY2= ACY*ACY
      ENDIF

      END



      SUBROUTINE genr70_P_gpu &
                (JTYPE, &
                 ISHELL,JSHELL,KSHELL,LSHELL, &
                 INEW,JNEW,KNEW,LNEW, &
                 P12,P34,R34,ACX,ACZ, &                 
                 RCD,SING,COSG,P,T)
      use mx_limits, only: mxgtot,mxsh,mxgsh,mxg2

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      double precision,parameter :: ZER=0.0D+00
      double precision,parameter :: PT5=0.5D+00
      double precision,parameter :: ONE=1.0D+00
      double precision,parameter :: PT7=0.7D+00
      double precision,parameter :: PT9=0.9D+00

!        ACYCUT CHANGED FROM A VALUE OF 1E-4 IN 2004,
!        FOR ACCURACY WHEN USING DIFFUSE EXPONENTS.
      double precision,parameter :: ACYCUT=1.0D-10
      double precision,parameter :: TENM12=1.0D-12

      ! PITO52=(PI+PI)*PI*SQRT(PI )
      double precision,parameter :: PITO52=34.986836655249726D+00

      COMMON /B     / CO(MXSH,3)
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS, &
                      TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2), &
                      D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2), &
                      NGANGB
!$omp threadprivate(/GEOMPQ/)

      COMMON /INTAC2/ EI1,EI2,CUX
      COMMON /MAXC  / CMAX(MXGTOT),CMAXA(MXGSH),CMAXB(MXGSH), &
                      CMAXC(MXGSH),CMAXD(MXGSH),ISMLP(MXG2),ISMLQ
!$omp threadprivate(/MAXC/)                      
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT), &
                      CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT), &
                      KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH), &
                      KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /POPOUT/ LPOPI,LPOPJ,LPOPK,LPOPL
!$omp threadprivate(/POPOUT/)      

      INTEGER :: INEW,JNEW,KNEW,LNEW

      COMMON /SHLLFO/ NGA,LA,EXA(MXGSH),CSA(MXGSH),CPA(MXGSH), &
                      NGB,LB,EXB(MXGSH),CSB(MXGSH),CPB(MXGSH), &
                      NGC,LC,EXC(MXGSH),CSC(MXGSH),CPC(MXGSH), &
                      NGD,LD,EXD(MXGSH),CSD(MXGSH),CPD(MXGSH)
!$omp threadprivate(/SHLLFO/)                      
      COMMON /SHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL, &
                      MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL, &
                      NIJ,IJ,KL,IJKL
!$omp threadprivate(/SHLNOS/)                      
!
!     ==================================================================
!
!JMS  LABELLED COMMON JMSGYH DEFINED FOR COMPUTATIONAL EFFICIENCY.
!JMS  IT IS ONLY USED IN THIS MODULE INT2B AND IN MODULE INT2R.
!
      COMMON /JMSGYH/ SQ(0:1,0:1)
!$omp threadprivate(/JMSGYH/)
      COMMON /KI2 / ACY,ACY2,AQX,AQX2,AQXY,Y03,Y04
!$omp threadprivate(/KI2/)      

      double precision :: P12(3,3),P34(3,3),P(3,3),T(3)



      PIF= PITO52


!
! DIRECTION COSINES OF CD LOCAL AXES WITH RESPECT TO AB LOCAL AXES
! (COSG,0,-SING)  (0,1,0)  (SING,0,COSG)
!
! PRELIMINARY P LOOP
!
! FILL GEOMPQ WITH INFORMATION ABOUT P IN PRELIMINARY P-LOOP
!
      JI= 1
      DO I=1,NGA
         X01= EXA(I)
         DO J=1,NGB
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
               JI=JI+1
               CYCLE
            ENDIF
            E12= X21* EXP(-R12Y12)
            TST= E12*CMAXA(I)*CMAXB(J)
            ISMLP(JI)=0
            IF(TST.LE.EI1) ISMLP(JI)=1
            IF(TST.LE.EI2) ISMLP(JI)=2
            E12= PIF*E12
!
! FOR TYPES 0000,0001,0011 ONLY D00P NEEDED
!
            D00P(JI)= E12*CSA(I)*CSB(J)
            IF(JTYPE.GT.3) THEN
               D01P(JI)= E12*CSA(I)*CPB(J)
               IF(JTYPE.LE.5) THEN
                  IF(D01P(JI).NE.ZER) THEN
                     D00P(JI)= D00P(JI)/D01P(JI)
                  ENDIF
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
            ENDIF
            JI=JI+1
         ENDDO
      ENDDO

      END




      SUBROUTINE genr70_Q1_gpu(JTYPE,P12,R34,ACX,ACZ,RCD,SING,COSG)
      use mx_limits, only: MXGTOT,MXSH,MXGSH,MXG2

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      double precision,parameter :: ONE=1.0D00

      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS, &
                      TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2), &
                      D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2), &
                      NGANGB
!$omp threadprivate(/GEOMPQ/)
      COMMON /INTAC2/ EI1,EI2,CUX
      COMMON /MAXC  / CMAX(MXGTOT),CMAXA(MXGSH),CMAXB(MXGSH), &
                      CMAXC(MXGSH),CMAXD(MXGSH),ISMLP(MXG2),ISMLQ
!$omp threadprivate(/MAXC/)                      
      COMMON /POPOUT/ LPOPI,LPOPJ,LPOPK,LPOPL
!$omp threadprivate(/POPOUT/)      
      COMMON /SHLLFO/ NGA,LA,EXA(MXGSH),CSA(MXGSH),CPA(MXGSH), &
                      NGB,LB,EXB(MXGSH),CSB(MXGSH),CPB(MXGSH), &
                      NGC,LC,EXC(MXGSH),CSC(MXGSH),CPC(MXGSH), &
                      NGD,LD,EXD(MXGSH),CSD(MXGSH),CPD(MXGSH)
!$omp threadprivate(/SHLLFO/)                      
      COMMON /SHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL, &
                      MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL, &
                      NIJ,IJ,KL,IJKL
!$omp threadprivate(/SHLNOS/)                      
!
!     ==================================================================
!
!JMS  LABELLED COMMON JMSGYH DEFINED FOR COMPUTATIONAL EFFICIENCY.
!JMS  IT IS ONLY USED IN THIS MODULE INT2B AND IN MODULE INT2R.
!
      COMMON /JMSGYH/ SQ(0:1,0:1)
!$omp threadprivate(/JMSGYH/)
      COMMON /KI2 / ACY,ACY2,AQX,AQX2,AQXY,Y03,Y04
!$omp threadprivate(/KI2/)      
      COMMON /KI3 / R00(25),R01(120),R02(156),R03(80),R04(15)
!$omp threadprivate(/KI3/)
      COMMON /FQ04  / FQ(0:4),FQ0(5),FQ1(2,9),FQ2(27),FQ3(20),FQ4(5)
!$omp threadprivate(/FQ04  /)      
      LOGICAL :: LRINT
      COMMON /NLRCF / LRINT
!$omp threadprivate(/NLRCF /)

!
! BEGIN Q LOOP
!
      IKL= 0
      DO K=1,NGC
         X03= EXC(K)
         DO L=1,NGD
            X04= EXD(L)
            X34= X03+X04
            X43= ONE/X34
            Y03= X03*X43
            Y04= ONE-Y03
            Y34= Y03*X04
            R34Y34= R34*Y34
            IF(R34Y34.GT.CUX) CYCLE
            E34= X43* EXP(-R34Y34)
            TST= E34*CMAXC(K)*CMAXD(L)
            IF(TST.LE.EI2) CYCLE
            ISMLQ= 0
            IF(TST.LE.EI1) ISMLQ= 1
!
! CQX = COMPONENT OF CQ ALONG PENULTIMATE X-AXIS
! CQZ = COMPONENT OF CQ ALONG PENULTIMATE Z-AXIS
!
            CQ = RCD*Y04
            CQX= CQ*SING
            CQZ= CQ*COSG
!
! FIND COORDINATES OF Q RELATIVE TO AXES AT A
! QPR IS PERPENDICULAR FROM Q TO AB
!
            AQX= ACX+CQX
            AQX2=AQX*AQX
            AQXY=AQX*ACY
            AQZ= ACZ+CQZ
            QPS= AQX2+ACY2

            SQ(0,0)= E34*CSC(K)*CSD(L)

!
! USE SPECIAL FAST ROUTINE FOR INNER LOOPS FOR 0000 ... 1111
!
!JMS  ZEROING OF THE FQx (x=0..4) ARRAYS OF LABELLED COMMON /FQ04/
!JMS  TAKES PLACE IN THE INTJx (x=1..6) SUBROUTINES
!
!JMS  ZEROING OF THE R0x (x=0..4) ARRAYS OF LABELLED COMMON /KI3 /
!JMS  TOOK PLACE IN THE CALL (FOR IKL=0) OF THE INTKx (x=2..6)
!JMS  SUBROUTINES
!

!
! FROM HERE, MCMURCHIE-DAVIDSON ALGORITHM IS USED.
!
! GENERATE (R] INTEGRALS
!
            IKL= IKL+1
            CALL INTJ1_gpu(LRINT,ISMLP,ISMLQ)
            R00(1)= R00(1)+FQ0(1)*SQ(0,0)

         ENDDO
      ENDDO

      END
      SUBROUTINE genr70_Q2_gpu(JTYPE,P12,R34,ACX,ACZ,RCD,SING,COSG)
      use mx_limits, only: MXGTOT,MXSH,MXGSH,MXG2

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      double precision,parameter :: ONE=1.0D00

      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS, &
                      TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2), &
                      D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2), &
                      NGANGB
!$omp threadprivate(/GEOMPQ/)
      COMMON /INTAC2/ EI1,EI2,CUX
      COMMON /MAXC  / CMAX(MXGTOT),CMAXA(MXGSH),CMAXB(MXGSH), &
                      CMAXC(MXGSH),CMAXD(MXGSH),ISMLP(MXG2),ISMLQ
!$omp threadprivate(/MAXC/)                      
      COMMON /POPOUT/ LPOPI,LPOPJ,LPOPK,LPOPL
!$omp threadprivate(/POPOUT/)      
      COMMON /SHLLFO/ NGA,LA,EXA(MXGSH),CSA(MXGSH),CPA(MXGSH), &
                      NGB,LB,EXB(MXGSH),CSB(MXGSH),CPB(MXGSH), &
                      NGC,LC,EXC(MXGSH),CSC(MXGSH),CPC(MXGSH), &
                      NGD,LD,EXD(MXGSH),CSD(MXGSH),CPD(MXGSH)
!$omp threadprivate(/SHLLFO/)                      
      COMMON /SHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL, &
                      MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL, &
                      NIJ,IJ,KL,IJKL
!$omp threadprivate(/SHLNOS/)                      
!
!     ==================================================================
!
!JMS  LABELLED COMMON JMSGYH DEFINED FOR COMPUTATIONAL EFFICIENCY.
!JMS  IT IS ONLY USED IN THIS MODULE INT2B AND IN MODULE INT2R.
!
      COMMON /JMSGYH/ SQ(0:1,0:1)
!$omp threadprivate(/JMSGYH/)
      COMMON /KI2 / ACY,ACY2,AQX,AQX2,AQXY,Y03,Y04
!$omp threadprivate(/KI2/)      
      COMMON /KI3 / R00(25),R01(120),R02(156),R03(80),R04(15)
!$omp threadprivate(/KI3/)
      COMMON /FQ04  / FQ(0:4),FQ0(5),FQ1(2,9),FQ2(27),FQ3(20),FQ4(5)
!$omp threadprivate(/FQ04  /)      
      LOGICAL :: LRINT
      COMMON /NLRCF / LRINT
!$omp threadprivate(/NLRCF /)


!
! BEGIN Q LOOP
!
      IKL= 0
      DO K=1,NGC
         X03= EXC(K)
         DO L=1,NGD
            X04= EXD(L)
            X34= X03+X04
            X43= ONE/X34
            Y03= X03*X43
            Y04= ONE-Y03
            Y34= Y03*X04
            R34Y34= R34*Y34
            IF(R34Y34.GT.CUX) CYCLE
            E34= X43* EXP(-R34Y34)
            TST= E34*CMAXC(K)*CMAXD(L)
            IF(TST.LE.EI2) CYCLE
            ISMLQ= 0
            IF(TST.LE.EI1) ISMLQ= 1
!
! CQX = COMPONENT OF CQ ALONG PENULTIMATE X-AXIS
! CQZ = COMPONENT OF CQ ALONG PENULTIMATE Z-AXIS
!
            CQ = RCD*Y04
            CQX= CQ*SING
            CQZ= CQ*COSG
!
! FIND COORDINATES OF Q RELATIVE TO AXES AT A
! QPR IS PERPENDICULAR FROM Q TO AB
!
            AQX= ACX+CQX
            AQX2=AQX*AQX
            AQXY=AQX*ACY
            AQZ= ACZ+CQZ
            QPS= AQX2+ACY2


            SQ(0,0)= E34*CSC(K)*CSD(L)
            SQ(1,0)= E34*CSC(K)*CPD(L)
            SQ(0,1)= E34*CPC(K)*CSD(L)
            SQ(1,1)= E34*CPC(K)*CPD(L)
!
! USE SPECIAL FAST ROUTINE FOR INNER LOOPS FOR 0000 ... 1111
!
!JMS  ZEROING OF THE FQx (x=0..4) ARRAYS OF LABELLED COMMON /FQ04/
!JMS  TAKES PLACE IN THE INTJx (x=1..6) SUBROUTINES
!
!JMS  ZEROING OF THE R0x (x=0..4) ARRAYS OF LABELLED COMMON /KI3 /
!JMS  TOOK PLACE IN THE CALL (FOR IKL=0) OF THE INTKx (x=2..6)
!JMS  SUBROUTINES
!

!
! FROM HERE, MCMURCHIE-DAVIDSON ALGORITHM IS USED.
!
! GENERATE (R] INTEGRALS
!
            IKL= IKL+1
            CALL INTJ2_gpu(LRINT,ISMLP,ISMLQ)
            CALL INTK2_2_gpu(IKL,SQ,X43,FQ0,FQ1, &
                             ACY,AQX,Y03, &
                             R00,R01)

         ENDDO
      ENDDO

      END




      SUBROUTINE genr70_Q3_gpu(JTYPE,P12,R34,ACX,ACZ,RCD,SING,COSG)
      use mx_limits, only: MXGTOT,MXSH,MXGSH,MXG2

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      double precision,parameter :: ONE=1.0D00

      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS, &
                      TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2), &
                      D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2), &
                      NGANGB
!$omp threadprivate(/GEOMPQ/)
      COMMON /INTAC2/ EI1,EI2,CUX
      COMMON /MAXC  / CMAX(MXGTOT),CMAXA(MXGSH),CMAXB(MXGSH), &
                      CMAXC(MXGSH),CMAXD(MXGSH),ISMLP(MXG2),ISMLQ
!$omp threadprivate(/MAXC/)                      
      COMMON /POPOUT/ LPOPI,LPOPJ,LPOPK,LPOPL
!$omp threadprivate(/POPOUT/)      
      COMMON /SHLLFO/ NGA,LA,EXA(MXGSH),CSA(MXGSH),CPA(MXGSH), &
                      NGB,LB,EXB(MXGSH),CSB(MXGSH),CPB(MXGSH), &
                      NGC,LC,EXC(MXGSH),CSC(MXGSH),CPC(MXGSH), &
                      NGD,LD,EXD(MXGSH),CSD(MXGSH),CPD(MXGSH)
!$omp threadprivate(/SHLLFO/)                      
      COMMON /SHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL, &
                      MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL, &
                      NIJ,IJ,KL,IJKL
!$omp threadprivate(/SHLNOS/)                      
!
!     ==================================================================
!
!JMS  LABELLED COMMON JMSGYH DEFINED FOR COMPUTATIONAL EFFICIENCY.
!JMS  IT IS ONLY USED IN THIS MODULE INT2B AND IN MODULE INT2R.
!
      COMMON /JMSGYH/ SQ(0:1,0:1)
!$omp threadprivate(/JMSGYH/)
      COMMON /KI2 / ACY,ACY2,AQX,AQX2,AQXY,Y03,Y04
!$omp threadprivate(/KI2/)      
      COMMON /KI3 / R00(25),R01(120),R02(156),R03(80),R04(15)
!$omp threadprivate(/KI3/)
      COMMON /FQ04  / FQ(0:4),FQ0(5),FQ1(2,9),FQ2(27),FQ3(20),FQ4(5)
!$omp threadprivate(/FQ04  /)      
      LOGICAL :: LRINT
      COMMON /NLRCF / LRINT
!$omp threadprivate(/NLRCF /)


!
! BEGIN Q LOOP
!
      IKL= 0
      DO K=1,NGC
         X03= EXC(K)
         DO L=1,NGD
            X04= EXD(L)
            X34= X03+X04
            X43= ONE/X34
            Y03= X03*X43
            Y04= ONE-Y03
            Y34= Y03*X04
            R34Y34= R34*Y34
            IF(R34Y34.GT.CUX) CYCLE
            E34= X43* EXP(-R34Y34)
            TST= E34*CMAXC(K)*CMAXD(L)
            IF(TST.LE.EI2) CYCLE
            ISMLQ= 0
            IF(TST.LE.EI1) ISMLQ= 1
!
! CQX = COMPONENT OF CQ ALONG PENULTIMATE X-AXIS
! CQZ = COMPONENT OF CQ ALONG PENULTIMATE Z-AXIS
!
            CQ = RCD*Y04
            CQX= CQ*SING
            CQZ= CQ*COSG
!
! FIND COORDINATES OF Q RELATIVE TO AXES AT A
! QPR IS PERPENDICULAR FROM Q TO AB
!
            AQX= ACX+CQX
            AQX2=AQX*AQX
            AQXY=AQX*ACY
            AQZ= ACZ+CQZ
            QPS= AQX2+ACY2

            SQ(0,0)= E34*CSC(K)*CSD(L)
            SQ(1,0)= E34*CSC(K)*CPD(L)
            SQ(0,1)= E34*CPC(K)*CSD(L)
            SQ(1,1)= E34*CPC(K)*CPD(L)
!
! USE SPECIAL FAST ROUTINE FOR INNER LOOPS FOR 0000 ... 1111
!
!JMS  ZEROING OF THE FQx (x=0..4) ARRAYS OF LABELLED COMMON /FQ04/
!JMS  TAKES PLACE IN THE INTJx (x=1..6) SUBROUTINES
!
!JMS  ZEROING OF THE R0x (x=0..4) ARRAYS OF LABELLED COMMON /KI3 /
!JMS  TOOK PLACE IN THE CALL (FOR IKL=0) OF THE INTKx (x=2..6)
!JMS  SUBROUTINES
!

!
! FROM HERE, MCMURCHIE-DAVIDSON ALGORITHM IS USED.
!
! GENERATE (R] INTEGRALS
!
            IKL= IKL+1
            CALL INTJ3_gpu(LRINT,ISMLP,ISMLQ)
            CALL INTK3_2_gpu(IKL,SQ,X43,FQ0,FQ1,FQ2, &
                             ACY,ACY2,AQX,AQX2,AQXY,Y03,Y04, &
                             R00,R01,R02)

         ENDDO
      ENDDO

      END



      SUBROUTINE genr70_Q4_gpu(JTYPE,P12,R34,ACX,ACZ,RCD,SING,COSG)
      use mx_limits, only: MXGTOT,MXSH,MXGSH,MXG2

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      double precision,parameter :: ONE=1.0D00

      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS, &
                      TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2), &
                      D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2), &
                      NGANGB
!$omp threadprivate(/GEOMPQ/)
      COMMON /INTAC2/ EI1,EI2,CUX
      COMMON /MAXC  / CMAX(MXGTOT),CMAXA(MXGSH),CMAXB(MXGSH), &
                      CMAXC(MXGSH),CMAXD(MXGSH),ISMLP(MXG2),ISMLQ
!$omp threadprivate(/MAXC/)                      
      COMMON /POPOUT/ LPOPI,LPOPJ,LPOPK,LPOPL
!$omp threadprivate(/POPOUT/)      
      COMMON /SHLLFO/ NGA,LA,EXA(MXGSH),CSA(MXGSH),CPA(MXGSH), &
                      NGB,LB,EXB(MXGSH),CSB(MXGSH),CPB(MXGSH), &
                      NGC,LC,EXC(MXGSH),CSC(MXGSH),CPC(MXGSH), &
                      NGD,LD,EXD(MXGSH),CSD(MXGSH),CPD(MXGSH)
!$omp threadprivate(/SHLLFO/)                      
      COMMON /SHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL, &
                      MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL, &
                      NIJ,IJ,KL,IJKL
!$omp threadprivate(/SHLNOS/)                      
!
!     ==================================================================
!
!JMS  LABELLED COMMON JMSGYH DEFINED FOR COMPUTATIONAL EFFICIENCY.
!JMS  IT IS ONLY USED IN THIS MODULE INT2B AND IN MODULE INT2R.
!
      COMMON /JMSGYH/ SQ(0:1,0:1)
!$omp threadprivate(/JMSGYH/)
      COMMON /KI2 / ACY,ACY2,AQX,AQX2,AQXY,Y03,Y04
!$omp threadprivate(/KI2/)      
      COMMON /KI3 / R00(25),R01(120),R02(156),R03(80),R04(15)
!$omp threadprivate(/KI3/)
      COMMON /FQ04  / FQ(0:4),FQ0(5),FQ1(2,9),FQ2(27),FQ3(20),FQ4(5)
!$omp threadprivate(/FQ04  /)      
      LOGICAL :: LRINT
      COMMON /NLRCF / LRINT
!$omp threadprivate(/NLRCF /)


!
! BEGIN Q LOOP
!
      IKL= 0
      DO K=1,NGC
         X03= EXC(K)
         DO L=1,NGD
            X04= EXD(L)
            X34= X03+X04
            X43= ONE/X34
            Y03= X03*X43
            Y04= ONE-Y03
            Y34= Y03*X04
            R34Y34= R34*Y34
            IF(R34Y34.GT.CUX) CYCLE
            E34= X43* EXP(-R34Y34)
            TST= E34*CMAXC(K)*CMAXD(L)
            IF(TST.LE.EI2) CYCLE
            ISMLQ= 0
            IF(TST.LE.EI1) ISMLQ= 1
!
! CQX = COMPONENT OF CQ ALONG PENULTIMATE X-AXIS
! CQZ = COMPONENT OF CQ ALONG PENULTIMATE Z-AXIS
!
            CQ = RCD*Y04
            CQX= CQ*SING
            CQZ= CQ*COSG
!
! FIND COORDINATES OF Q RELATIVE TO AXES AT A
! QPR IS PERPENDICULAR FROM Q TO AB
!
            AQX= ACX+CQX
            AQX2=AQX*AQX
            AQXY=AQX*ACY
            AQZ= ACZ+CQZ
            QPS= AQX2+ACY2

            IF(JTYPE.NE.1) THEN
               SQ(0,0)= E34*CSC(K)*CSD(L)
               SQ(1,0)= E34*CSC(K)*CPD(L)
               SQ(0,1)= E34*CPC(K)*CSD(L)
               SQ(1,1)= E34*CPC(K)*CPD(L)
            ELSE
               SQ(0,0)= E34*CSC(K)*CSD(L)
            ENDIF
!
! USE SPECIAL FAST ROUTINE FOR INNER LOOPS FOR 0000 ... 1111
!
!JMS  ZEROING OF THE FQx (x=0..4) ARRAYS OF LABELLED COMMON /FQ04/
!JMS  TAKES PLACE IN THE INTJx (x=1..6) SUBROUTINES
!
!JMS  ZEROING OF THE R0x (x=0..4) ARRAYS OF LABELLED COMMON /KI3 /
!JMS  TOOK PLACE IN THE CALL (FOR IKL=0) OF THE INTKx (x=2..6)
!JMS  SUBROUTINES
!

!
! FROM HERE, MCMURCHIE-DAVIDSON ALGORITHM IS USED.
!
! GENERATE (R] INTEGRALS
!
            IKL= IKL+1
            CALL INTJ4_gpu(LRINT,ISMLP,ISMLQ)
            CALL INTK4_2_gpu(IKL,SQ,X43,FQ0,FQ1,FQ2, &
                             ACY,ACY2,AQX,AQX2,AQXY,Y03, &
                             R00,R01,R02)

         ENDDO
      ENDDO

      END



      SUBROUTINE genr70_Q5_gpu(JTYPE,P12,R34,ACX,ACZ,RCD,SING,COSG)
      use mx_limits, only: MXGTOT,MXSH,MXGSH,MXG2

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      double precision,parameter :: ONE=1.0D00

      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS, &
                      TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2), &
                      D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2), &
                      NGANGB
!$omp threadprivate(/GEOMPQ/)
      COMMON /INTAC2/ EI1,EI2,CUX
      COMMON /MAXC  / CMAX(MXGTOT),CMAXA(MXGSH),CMAXB(MXGSH), &
                      CMAXC(MXGSH),CMAXD(MXGSH),ISMLP(MXG2),ISMLQ
!$omp threadprivate(/MAXC/)                      
      COMMON /POPOUT/ LPOPI,LPOPJ,LPOPK,LPOPL
!$omp threadprivate(/POPOUT/)      
      COMMON /SHLLFO/ NGA,LA,EXA(MXGSH),CSA(MXGSH),CPA(MXGSH), &
                      NGB,LB,EXB(MXGSH),CSB(MXGSH),CPB(MXGSH), &
                      NGC,LC,EXC(MXGSH),CSC(MXGSH),CPC(MXGSH), &
                      NGD,LD,EXD(MXGSH),CSD(MXGSH),CPD(MXGSH)
!$omp threadprivate(/SHLLFO/)                      
      COMMON /SHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL, &
                      MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL, &
                      NIJ,IJ,KL,IJKL
!$omp threadprivate(/SHLNOS/)                      
!
!     ==================================================================
!
!JMS  LABELLED COMMON JMSGYH DEFINED FOR COMPUTATIONAL EFFICIENCY.
!JMS  IT IS ONLY USED IN THIS MODULE INT2B AND IN MODULE INT2R.
!
      COMMON /JMSGYH/ SQ(0:1,0:1)
!$omp threadprivate(/JMSGYH/)
      COMMON /KI2 / ACY,ACY2,AQX,AQX2,AQXY,Y03,Y04
!$omp threadprivate(/KI2/)      
      COMMON /KI3 / R00(25),R01(120),R02(156),R03(80),R04(15)
!$omp threadprivate(/KI3/)
      COMMON /FQ04  / FQ(0:4),FQ0(5),FQ1(2,9),FQ2(27),FQ3(20),FQ4(5)
!$omp threadprivate(/FQ04  /)      
      LOGICAL :: LRINT
      COMMON /NLRCF / LRINT
!$omp threadprivate(/NLRCF /)


!
! BEGIN Q LOOP
!
      IKL= 0
      DO K=1,NGC
         X03= EXC(K)
         DO L=1,NGD
            X04= EXD(L)
            X34= X03+X04
            X43= ONE/X34
            Y03= X03*X43
            Y04= ONE-Y03
            Y34= Y03*X04
            R34Y34= R34*Y34
            IF(R34Y34.GT.CUX) CYCLE
            E34= X43* EXP(-R34Y34)
            TST= E34*CMAXC(K)*CMAXD(L)
            IF(TST.LE.EI2) CYCLE
            ISMLQ= 0
            IF(TST.LE.EI1) ISMLQ= 1
!
! CQX = COMPONENT OF CQ ALONG PENULTIMATE X-AXIS
! CQZ = COMPONENT OF CQ ALONG PENULTIMATE Z-AXIS
!
            CQ = RCD*Y04
            CQX= CQ*SING
            CQZ= CQ*COSG
!
! FIND COORDINATES OF Q RELATIVE TO AXES AT A
! QPR IS PERPENDICULAR FROM Q TO AB
!
            AQX= ACX+CQX
            AQX2=AQX*AQX
            AQXY=AQX*ACY
            AQZ= ACZ+CQZ
            QPS= AQX2+ACY2

            SQ(0,0)= E34*CSC(K)*CSD(L)
            SQ(1,0)= E34*CSC(K)*CPD(L)
            SQ(0,1)= E34*CPC(K)*CSD(L)
            SQ(1,1)= E34*CPC(K)*CPD(L)
!
! USE SPECIAL FAST ROUTINE FOR INNER LOOPS FOR 0000 ... 1111
!
!JMS  ZEROING OF THE FQx (x=0..4) ARRAYS OF LABELLED COMMON /FQ04/
!JMS  TAKES PLACE IN THE INTJx (x=1..6) SUBROUTINES
!
!JMS  ZEROING OF THE R0x (x=0..4) ARRAYS OF LABELLED COMMON /KI3 /
!JMS  TOOK PLACE IN THE CALL (FOR IKL=0) OF THE INTKx (x=2..6)
!JMS  SUBROUTINES
!

!
! FROM HERE, MCMURCHIE-DAVIDSON ALGORITHM IS USED.
!
! GENERATE (R] INTEGRALS
!
            IKL= IKL+1
            CALL INTJ5_gpu(LRINT,ISMLP,ISMLQ)
            CALL INTK5_2_gpu(IKL,SQ,X43,FQ0,FQ1,FQ2,FQ3, &
                             ACY,ACY2,AQX,AQX2,AQXY,Y03,Y04, &
                             R00,R01,R02,R03)

         ENDDO
      ENDDO

      END



      SUBROUTINE genr70_Q6_gpu(JTYPE,P12,R34,ACX,ACZ,RCD,SING,COSG)
      use mx_limits, only: MXGTOT,MXSH,MXGSH,MXG2

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      double precision,parameter :: ONE=1.0D00

      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS, &
                      TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2), &
                      D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2), &
                      NGANGB
!$omp threadprivate(/GEOMPQ/)
      COMMON /INTAC2/ EI1,EI2,CUX
      COMMON /MAXC  / CMAX(MXGTOT),CMAXA(MXGSH),CMAXB(MXGSH), &
                      CMAXC(MXGSH),CMAXD(MXGSH),ISMLP(MXG2),ISMLQ
!$omp threadprivate(/MAXC/)                      
      COMMON /POPOUT/ LPOPI,LPOPJ,LPOPK,LPOPL
!$omp threadprivate(/POPOUT/)      
      COMMON /SHLLFO/ NGA,LA,EXA(MXGSH),CSA(MXGSH),CPA(MXGSH), &
                      NGB,LB,EXB(MXGSH),CSB(MXGSH),CPB(MXGSH), &
                      NGC,LC,EXC(MXGSH),CSC(MXGSH),CPC(MXGSH), &
                      NGD,LD,EXD(MXGSH),CSD(MXGSH),CPD(MXGSH)
!$omp threadprivate(/SHLLFO/)                      
      COMMON /SHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL, &
                      MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL, &
                      NIJ,IJ,KL,IJKL
!$omp threadprivate(/SHLNOS/)                      
!
!     ==================================================================
!
!JMS  LABELLED COMMON JMSGYH DEFINED FOR COMPUTATIONAL EFFICIENCY.
!JMS  IT IS ONLY USED IN THIS MODULE INT2B AND IN MODULE INT2R.
!
      COMMON /JMSGYH/ SQ(0:1,0:1)
!$omp threadprivate(/JMSGYH/)
      COMMON /KI2 / ACY,ACY2,AQX,AQX2,AQXY,Y03,Y04
!$omp threadprivate(/KI2/)      
      COMMON /KI3 / R00(25),R01(120),R02(156),R03(80),R04(15)
!$omp threadprivate(/KI3/)
      COMMON /FQ04  / FQ(0:4),FQ0(5),FQ1(2,9),FQ2(27),FQ3(20),FQ4(5)
!$omp threadprivate(/FQ04  /)      
      LOGICAL :: LRINT
      COMMON /NLRCF / LRINT
!$omp threadprivate(/NLRCF /)


!
! BEGIN Q LOOP
!
      IKL= 0
      DO K=1,NGC
         X03= EXC(K)
         DO L=1,NGD
            X04= EXD(L)
            X34= X03+X04
            X43= ONE/X34
            Y03= X03*X43
            Y04= ONE-Y03
            Y34= Y03*X04
            R34Y34= R34*Y34
            IF(R34Y34.GT.CUX) CYCLE
            E34= X43* EXP(-R34Y34)
            TST= E34*CMAXC(K)*CMAXD(L)
            IF(TST.LE.EI2) CYCLE
            ISMLQ= 0
            IF(TST.LE.EI1) ISMLQ= 1
!
! CQX = COMPONENT OF CQ ALONG PENULTIMATE X-AXIS
! CQZ = COMPONENT OF CQ ALONG PENULTIMATE Z-AXIS
!
            CQ = RCD*Y04
            CQX= CQ*SING
            CQZ= CQ*COSG
!
! FIND COORDINATES OF Q RELATIVE TO AXES AT A
! QPR IS PERPENDICULAR FROM Q TO AB
!
            AQX= ACX+CQX
            AQX2=AQX*AQX
            AQXY=AQX*ACY
            AQZ= ACZ+CQZ
            QPS= AQX2+ACY2

            IF(JTYPE.NE.1) THEN
               SQ(0,0)= E34*CSC(K)*CSD(L)
               SQ(1,0)= E34*CSC(K)*CPD(L)
               SQ(0,1)= E34*CPC(K)*CSD(L)
               SQ(1,1)= E34*CPC(K)*CPD(L)
            ELSE
               SQ(0,0)= E34*CSC(K)*CSD(L)
            ENDIF
!
! USE SPECIAL FAST ROUTINE FOR INNER LOOPS FOR 0000 ... 1111
!
!JMS  ZEROING OF THE FQx (x=0..4) ARRAYS OF LABELLED COMMON /FQ04/
!JMS  TAKES PLACE IN THE INTJx (x=1..6) SUBROUTINES
!
!JMS  ZEROING OF THE R0x (x=0..4) ARRAYS OF LABELLED COMMON /KI3 /
!JMS  TOOK PLACE IN THE CALL (FOR IKL=0) OF THE INTKx (x=2..6)
!JMS  SUBROUTINES
!

!
! FROM HERE, MCMURCHIE-DAVIDSON ALGORITHM IS USED.
!
! GENERATE (R] INTEGRALS
!
            IKL= IKL+1
            CALL INTJ6_gpu(LRINT,ISMLP,ISMLQ)
            CALL INTK6_2_gpu(IKL,SQ,RAB,X43,FQ0,FQ1,FQ2,FQ3,FQ4, &
                             ACY,ACY2,AQX,AQX2,AQXY,Y03,Y04, &
                             R00,R01,R02,R03,R04)
         ENDDO
      ENDDO

      END



!       SUBROUTINE genr70_Q_gpu(JTYPE,P12,R34,ACX,ACZ,RCD,SING,COSG)
!       use mx_limits, only: MXGTOT,MXSH,MXGSH,MXG2

!       IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!       double precision,parameter :: ONE=1.0D00

!       COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS, &
!                       TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2), &
!                       D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2), &
!                       NGANGB
! !$omp threadprivate(/GEOMPQ/)
!       COMMON /INTAC2/ EI1,EI2,CUX
!       COMMON /MAXC  / CMAX(MXGTOT),CMAXA(MXGSH),CMAXB(MXGSH), &
!                       CMAXC(MXGSH),CMAXD(MXGSH),ISMLP(MXG2),ISMLQ
! !$omp threadprivate(/MAXC/)                      
!       COMMON /POPOUT/ LPOPI,LPOPJ,LPOPK,LPOPL
! !$omp threadprivate(/POPOUT/)      
!       COMMON /SHLLFO/ NGA,LA,EXA(MXGSH),CSA(MXGSH),CPA(MXGSH), &
!                       NGB,LB,EXB(MXGSH),CSB(MXGSH),CPB(MXGSH), &
!                       NGC,LC,EXC(MXGSH),CSC(MXGSH),CPC(MXGSH), &
!                       NGD,LD,EXD(MXGSH),CSD(MXGSH),CPD(MXGSH)
! !$omp threadprivate(/SHLLFO/)                      
!       COMMON /SHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL, &
!                       MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL, &
!                       NIJ,IJ,KL,IJKL
! !$omp threadprivate(/SHLNOS/)                      
! !
! !     ==================================================================
! !
! !JMS  LABELLED COMMON JMSGYH DEFINED FOR COMPUTATIONAL EFFICIENCY.
! !JMS  IT IS ONLY USED IN THIS MODULE INT2B AND IN MODULE INT2R.
! !
!       COMMON /JMSGYH/ SQ(0:1,0:1)
! !$omp threadprivate(/JMSGYH/)
!       COMMON /KI2 / ACY,ACY2,AQX,AQX2,AQXY,Y03,Y04
! !$omp threadprivate(/KI2/)      
!       COMMON /KI3 / R00(25),R01(120),R02(156),R03(80),R04(15)
! !$omp threadprivate(/KI3/)
!       COMMON /FQ04  / FQ(0:4),FQ0(5),FQ1(2,9),FQ2(27),FQ3(20),FQ4(5)
! !$omp threadprivate(/FQ04  /)      
!       LOGICAL :: LRINT
!       COMMON /NLRCF / LRINT
! !$omp threadprivate(/NLRCF /)



!       CALL OMP_SET_DYNAMIC(.FALSE.)



! !
! ! BEGIN Q LOOP
! !
!       IKL= 0
!       DO K=1,NGC
!          X03= EXC(K)
!          DO L=1,NGD
!             X04= EXD(L)
!             X34= X03+X04
!             X43= ONE/X34
!             Y03= X03*X43
!             Y04= ONE-Y03
!             Y34= Y03*X04
!             R34Y34= R34*Y34
!             IF(R34Y34.GT.CUX) CYCLE
!             E34= X43* EXP(-R34Y34)
!             TST= E34*CMAXC(K)*CMAXD(L)
!             IF(TST.LE.EI2) CYCLE
!             ISMLQ= 0
!             IF(TST.LE.EI1) ISMLQ= 1
! !
! ! CQX = COMPONENT OF CQ ALONG PENULTIMATE X-AXIS
! ! CQZ = COMPONENT OF CQ ALONG PENULTIMATE Z-AXIS
! !
!             CQ = RCD*Y04
!             CQX= CQ*SING
!             CQZ= CQ*COSG
! !
! ! FIND COORDINATES OF Q RELATIVE TO AXES AT A
! ! QPR IS PERPENDICULAR FROM Q TO AB
! !
!             AQX= ACX+CQX
!             AQX2=AQX*AQX
!             AQXY=AQX*ACY
!             AQZ= ACZ+CQZ
!             QPS= AQX2+ACY2

!             IF(JTYPE.NE.1) THEN
!                SQ(0,0)= E34*CSC(K)*CSD(L)
!                SQ(1,0)= E34*CSC(K)*CPD(L)
!                SQ(0,1)= E34*CPC(K)*CSD(L)
!                SQ(1,1)= E34*CPC(K)*CPD(L)
!             ELSE
!                SQ(0,0)= E34*CSC(K)*CSD(L)
!             ENDIF
! !
! ! USE SPECIAL FAST ROUTINE FOR INNER LOOPS FOR 0000 ... 1111
! !
! !JMS  ZEROING OF THE FQx (x=0..4) ARRAYS OF LABELLED COMMON /FQ04/
! !JMS  TAKES PLACE IN THE INTJx (x=1..6) SUBROUTINES
! !
! !JMS  ZEROING OF THE R0x (x=0..4) ARRAYS OF LABELLED COMMON /KI3 /
! !JMS  TOOK PLACE IN THE CALL (FOR IKL=0) OF THE INTKx (x=2..6)
! !JMS  SUBROUTINES
! !

! !
! ! FROM HERE, MCMURCHIE-DAVIDSON ALGORITHM IS USED.
! !
! ! GENERATE (R] INTEGRALS
! !
!             IKL= IKL+1
!             IF(JTYPE.EQ. 1) THEN
!                CALL INTJ1_gpu(LRINT,ISMLP,ISMLQ)
!                R00(1)= R00(1)+FQ0(1)*SQ(0,0)
!             ELSEIF(JTYPE.EQ. 2) THEN
!                CALL INTJ2_gpu(LRINT,ISMLP,ISMLQ)
!                CALL INTK2_2_gpu(IKL,SQ,X43,FQ0,FQ1, &
!                                 ACY,AQX,Y03, &
!                                 R00,R01)
!             ELSEIF(JTYPE.EQ. 3) THEN
!                CALL INTJ3_gpu(LRINT,ISMLP,ISMLQ)
!                CALL INTK3_2_gpu(IKL,SQ,X43,FQ0,FQ1,FQ2, &
!                                 ACY,ACY2,AQX,AQX2,AQXY,Y03,Y04, &
!                                 R00,R01,R02)
!             ELSEIF(JTYPE.EQ. 4) THEN
!                CALL INTJ4_gpu(LRINT,ISMLP,ISMLQ)
!                CALL INTK4_2_gpu(IKL,SQ,X43,FQ0,FQ1,FQ2, &
!                                 ACY,ACY2,AQX,AQX2,AQXY,Y03, &
!                                 R00,R01,R02)
!             ELSEIF(JTYPE.EQ. 5) THEN
!                CALL INTJ5_gpu(LRINT,ISMLP,ISMLQ)
!                CALL INTK5_2_gpu(IKL,SQ,X43,FQ0,FQ1,FQ2,FQ3, &
!                                 ACY,ACY2,AQX,AQX2,AQXY,Y03,Y04, &
!                                 R00,R01,R02,R03)
!             ELSEIF(JTYPE.EQ. 6) THEN
!                CALL INTJ6_gpu(LRINT,ISMLP,ISMLQ)
!                CALL INTK6_2_gpu(IKL,SQ,RAB,X43,FQ0,FQ1,FQ2,FQ3,FQ4, &
!                                 ACY,ACY2,AQX,AQX2,AQXY,Y03,Y04, &
!                                 R00,R01,R02,R03,R04)
!             ENDIF



!          ENDDO
!       ENDDO

!       END


!*MODULE INT2B   *DECK INTJ1
!>
!>    @brief   SSSS case
!>
!>    @details integration of a SSSS case
!>
      SUBROUTINE INTJ1_gpu(LRINT,ISMLP,ISMLQ)
      use mx_limits, only: MXGTOT,MXGSH,MXG2
      USE lrcdft, ONLY: EMU2

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! GENERATE JTYPE= 1 INTEGRALS
!
      double precision,parameter :: ZER=0.0D00
      double precision,parameter :: ONE=1.0D00
      double precision,parameter :: PI4=0.78539816339744831D00
      integer,parameter :: NTX=4
      integer,parameter :: NPF=450
      integer,parameter :: NGRD=7
      integer,parameter :: NPX=1000
      integer,parameter :: MXQT=16

      COMMON /FMTTBL/ FGRID(0:NTX,0:NPF,0:NGRD),XGRID(0:NTX,0:NPX), &
                      TMAX,RFINC(0:NGRD),RXINC, &
                      RMR(MXQT),TLGM(0:MXQT),NORD
      COMMON /FQ04  / FQ(0:4),FQ0(5),FQ1(18),FQ2(27),FQ3(20),FQ4(5)
!$omp threadprivate(/FQ04  /)
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS, &
                      TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2), &
                      D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2), &
                      NGANGB
!$omp threadprivate(/GEOMPQ/)
      integer :: ISMLP(MXG2),ISMLQ
      LOGICAL :: LRINT



      FQ0(1)= ZER

      DO I=1,NGANGB
         ISML= ISMLQ+ISMLP(I)
         IF(ISML.GE.2) CYCLE
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
         END IF
         XVA=(PQS+QPS)*RHO
         RHO= RHO+RHO
         N=0
         IF(XVA.LE.TMAX) THEN
!
!     FM(T) EVALUATION
!
            TV= XVA*RFINC(N)
            IP= NINT(TV)
            FX=    FGRID(4,IP,N) *TV
            FX=(FX+FGRID(3,IP,N))*TV
            FX=(FX+FGRID(2,IP,N))*TV
            FX=(FX+FGRID(1,IP,N))*TV
            FX= FX+FGRID(0,IP,N)

            FQ(N)= FX
            FQF= FQZ*SQRT(X41)
               FQ(0)= FQ(0)*FQF
         ELSE
            XIN= ONE/XVA
            FQ(0)= FQZ*SQRT(PI4*XIN*X41)
         ENDIF

         FQ0(1)= FQ0(1)+FQ(0)
      ENDDO

      END


!*MODULE INT2B   *DECK INTJ2
!>
!>    @brief   PSSS case
!>
!>    @details integration of a PSSS case
!>
      SUBROUTINE INTJ2_gpu(LRINT,ISMLP,ISMLQ)
      USE lrcdft, ONLY: LCFLAG, EMU, EMU2, LRFILE
      use mx_limits, only: mxgtot,mxgsh,mxg2

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! GENERATE JTYPE= 2 INTEGRALS
!
      double precision,parameter :: ZER=0.0D+00
      double precision,parameter :: PT5=0.5D+00
      double precision,parameter :: ONE=1.0D+00
      double precision,parameter :: PI4=0.78539816339744831D+00
      integer,parameter :: NTX=4
      integer,parameter :: NPF=450
      integer,parameter :: NGRD=7
      integer,parameter :: NPX=1000
      integer,parameter :: MXQT=16

      COMMON /FMTTBL/ FGRID(0:NTX,0:NPF,0:NGRD),XGRID(0:NTX,0:NPX), &
                      TMAX,RFINC(0:NGRD),RXINC, &
                      RMR(MXQT),TLGM(0:MXQT),NORD
      COMMON /FQ04  / FQ(0:4),FQ0(5),FQ1(18),FQ2(27),FQ3(20),FQ4(5)
!$omp threadprivate(/FQ04  /)
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS, &
                      TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2), &
                      D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2), &
                      NGANGB
!$omp threadprivate(/GEOMPQ/)
      integer :: ISMLP(MXG2),ISMLQ
      LOGICAL :: LRINT



      FQ0(1)= ZER
      FQ1(1)= ZER
      FQ1(2)= ZER

      DO I=1,NGANGB
         ISML= ISMLQ+ISMLP(I)
         IF(ISML.GE.2) CYCLE
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
         N=1
         IF(XVA.LE.TMAX) THEN
!
!     FM(T) EVALUATION
!
            TV= XVA*RFINC(0)
            IP= NINT(TV)
            FX=    FGRID(4,IP,0) *TV
            FX=(FX+FGRID(3,IP,0))*TV
            FX=(FX+FGRID(2,IP,0))*TV
            FX=(FX+FGRID(1,IP,0))*TV
            FX= FX+FGRID(0,IP,0)

            FQ(0)=FX

            TV= XVA*RFINC(N)
            IP= NINT(TV)
            FX=    FGRID(4,IP,N) *TV
            FX=(FX+FGRID(3,IP,N))*TV
            FX=(FX+FGRID(2,IP,N))*TV
            FX=(FX+FGRID(1,IP,N))*TV
            FX= FX+FGRID(0,IP,N)

            FQ(N)= FX
            FQF= FQZ*SQRT(X41)
               FQ(0)= FQ(0)*FQF
            FQF= FQF*RHO
               FQ(1)= FQ(1)*FQF
         ELSE
            XIN= ONE/XVA
            FQ(0)= FQZ*SQRT(PI4*XIN*X41)
            ROX= RHO*XIN
            FQF= PT5*ROX
               FQ(1)= FQ(0)*FQF
         ENDIF

         FQ0(1)= FQ0(1)+FQ(0)
         FQ1(1)= FQ1(1)+FQ(1)
         FQ1(2)= FQ1(2)+FQ(1)*PQR

      ENDDO

      END


!*MODULE INT2B   *DECK INTJ3
!>
!>    @brief   PPSS case
!>
!>    @details integration of a PPSS case
!>
      SUBROUTINE INTJ3_gpu(LRINT,ISMLP,ISMLQ)
      USE lrcdft, ONLY: LCFLAG, EMU, EMU2, LRFILE
      use mx_limits, only: mxgtot,mxgsh,mxg2

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! GENERATE JTYPE= 3 INTEGRALS
!
      double precision,parameter :: ZER=0.0D+00
      double precision,parameter :: PT5=0.5D+00
      double precision,parameter :: ONE=1.0D+00
      double precision,parameter :: PI4=0.78539816339744831D+00
      integer,parameter :: NTX=4
      integer,parameter :: NPF=450
      integer,parameter :: NGRD=7
      integer,parameter :: NPX=1000
      integer,parameter :: MXQT=16
      COMMON /FMTTBL/ FGRID(0:NTX,0:NPF,0:NGRD),XGRID(0:NTX,0:NPX), &
                      TMAX,RFINC(0:NGRD),RXINC, &
                      RMR(MXQT),TLGM(0:MXQT),NORD
      COMMON /FQ04  / FQ(0:4),FQ0(5),FQ1(18),FQ2(27),FQ3(20),FQ4(5)
!$omp threadprivate(/FQ04  /)
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS, &
                      TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2), &
                      D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2), &
                      NGANGB
!$omp threadprivate(/GEOMPQ/)
      integer :: ISMLP(MXG2),ISMLQ
      LOGICAL :: LRINT



      FQ0(1)= ZER
      FQ1(1)= ZER
      FQ1(2)= ZER
      FQ2(1)= ZER
      FQ2(2)= ZER
      FQ2(3)= ZER

      DO I=1,NGANGB
         ISML= ISMLQ+ISMLP(I)
         IF(ISML.GE.2) CYCLE
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
!
!     FM(T) EVALUATION
!
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

            FQ(N)= FX
            T2= XVA+XVA
               FQ(2-1)=(T2*FQ(2)+ET)*RMR(2)
               FQ(1-1)=(T2*FQ(1)+ET)*RMR(1)

            FQF= FQZ*SQRT(X41)
               FQ(0)= FQ(0)*FQF
            FQF= FQF*RHO
               FQ(1)= FQ(1)*FQF
            FQF= FQF*RHO
               FQ(2)= FQ(2)*FQF
         ELSE
            XIN= ONE/XVA
            FQ(0)= FQZ*SQRT(PI4*XIN*X41)
            ROX= RHO*XIN
            FQF= PT5*ROX
               FQ(1)= FQ(0)*FQF
            FQF= FQF+ROX
               FQ(2)= FQ(1)*FQF
         ENDIF

         FQ0(1)= FQ0(1)+FQ(0)
         FQ1(1)= FQ1(1)+FQ(1)
         FQ1(2)= FQ1(2)+FQ(1)*PQR
         FQ2(1)= FQ2(1)+FQ(2)
         FQ2(2)= FQ2(2)+FQ(2)*PQR
         FQ2(3)= FQ2(3)+FQ(2)*PQS
      ENDDO

      END


!*MODULE INT2B   *DECK INTJ4
!>
!>    @brief   PSPS case
!>
!>    @details integration of a PSPS case
!>
      SUBROUTINE INTJ4_gpu(LRINT,ISMLP,ISMLQ)
      USE lrcdft, ONLY: LCFLAG, EMU, EMU2, LRFILE
      use mx_limits, only: mxgtot,mxgsh,mxg2

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! GENERATE JTYPE= 4 INTEGRALS
!
      double precision,parameter :: ZER=0.0D+00
      double precision,parameter :: PT5=0.5D+00
      double precision,parameter :: ONE=1.0D+00
      double precision,parameter :: PI4=0.78539816339744831D+00
      integer,parameter :: NTX=4
      integer,parameter :: NPF=450
      integer,parameter :: NGRD=7
      integer,parameter :: NPX=1000
      integer,parameter :: MXQT=16
      COMMON /FMTTBL/ FGRID(0:NTX,0:NPF,0:NGRD),XGRID(0:NTX,0:NPX), & 
                      TMAX,RFINC(0:NGRD),RXINC, & 
                      RMR(MXQT),TLGM(0:MXQT),NORD
      COMMON /FQ04  / FQ(0:4),FQ0(5),FQ1(18),FQ2(27),FQ3(20),FQ4(5)
!$omp threadprivate(/FQ04  /)
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS, &
                      TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2), &
                      D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2), &
                      NGANGB
!$omp threadprivate(/GEOMPQ/)
      LOGICAL :: LRINT
      integer :: ISMLP(MXG2),ISMLQ

      FQ0(1)= ZER
      FQ0(2)= ZER
      DO J=1,6
         FQ1(J)= ZER
      ENDDO
      FQ2(1)= ZER
      FQ2(2)= ZER
      FQ2(3)= ZER

      DO I=1,NGANGB
         ISML= ISMLQ+ISMLP(I)
         IF(ISML.GE.2) CYCLE
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
         N=2
         IF(XVA.LE.TMAX) THEN
!
!     FM(T) EVALUATION
!
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

            FQ(N)= FX
            T2= XVA+XVA
               FQ(2-1)=(T2*FQ(2)+ET)*RMR(2)
               FQ(1-1)=(T2*FQ(1)+ET)*RMR(1)

            FQF= FQZ*SQRT(X41)
               FQ(0)= FQ(0)*FQF
            FQF= FQF*RHO
               FQ(1)= FQ(1)*FQF
            FQF= FQF*RHO
               FQ(2)= FQ(2)*FQF
         ELSE
            XIN= ONE/XVA
            FQ(0)= FQZ*SQRT(PI4*XIN*X41)
            ROX= RHO*XIN
            FQF= PT5*ROX
               FQ(1)= FQ(0)*FQF
            FQF= FQF+ROX
               FQ(2)= FQ(1)*FQF
         ENDIF

         XMD1= TX21(I)
         Y01 = TY01(I)
         DP00= D00P(I)
         TMP1= XMD1*PQR
         TMP2= XMD1*PQS

         FQ0(1)= FQ0(1)+FQ(0)*DP00
         FQ0(2)= FQ0(2)+FQ(0)*Y01
         FQ1(1)= FQ1(1)+FQ(1)*DP00
         FQ1(2)= FQ1(2)+FQ(1)*DP00*PQR
         FQ1(3)= FQ1(3)+FQ(1)*XMD1
         FQ1(4)= FQ1(4)+FQ(1)*TMP1
         FQ1(5)= FQ1(5)+FQ(1)*Y01
         FQ1(6)= FQ1(6)+FQ(1)*Y01 *PQR
         FQ2(1)= FQ2(1)+FQ(2)*XMD1
         FQ2(2)= FQ2(2)+FQ(2)*TMP1
         FQ2(3)= FQ2(3)+FQ(2)*TMP2

      ENDDO

      END

!*MODULE INT2B   *DECK INTJ5
!>
!>    @brief   PPPS case
!>
!>    @details integration of a PPPS case
!>
      SUBROUTINE INTJ5_gpu(LRINT,ISMLP,ISMLQ)
      USE lrcdft, ONLY: LCFLAG, EMU, EMU2, LRFILE
      use mx_limits, only: mxgtot,mxgsh,mxg2

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! GENERATE JTYPE= 5 INTEGRALS
!
      double precision,parameter :: ZER=0.0D+00
      double precision,parameter :: PT5=0.5D+00
      double precision,parameter :: ONE=1.0D+00
      double precision,parameter :: PI4=0.78539816339744831D+00
      integer,parameter :: NTX=4
      integer,parameter :: NPF=450
      integer,parameter :: NGRD=7
      integer,parameter :: NPX=1000
      integer,parameter :: MXQT=16
      COMMON /FMTTBL/ FGRID(0:NTX,0:NPF,0:NGRD),XGRID(0:NTX,0:NPX), &
                      TMAX,RFINC(0:NGRD),RXINC, &
                      RMR(MXQT),TLGM(0:MXQT),NORD
      COMMON /FQ04  / FQ(0:4),FQ0(5),FQ1(2,9),FQ2(3,9),FQ3(20),FQ4(5)
!$omp threadprivate(/FQ04  /)
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS, &
                      TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2), &
                      D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2), &
                      NGANGB
!$omp threadprivate(/GEOMPQ/)
      double precision :: WORK(3,3)
      LOGICAL :: LRINT
      integer :: ISMLP(MXG2),ISMLQ


      FQ0(1)= ZER
      FQ0(2)= ZER
      DO J=1,3
         FQ1(1,J)= ZER
         FQ1(2,J)= ZER

         FQ2(1,J)= ZER
         FQ2(2,J)= ZER
         FQ2(3,J)= ZER
      ENDDO
      FQ3(1)= ZER
      FQ3(2)= ZER
      FQ3(3)= ZER
      FQ3(4)= ZER

      DO I=1,NGANGB
         ISML= ISMLQ+ISMLP(I)
         IF(ISML.GE.2) CYCLE
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
!
!     FM(T) EVALUATION...DOWNWARD RECURSION FOR JTYPE.GE.5
!
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

            FQ(N)= FX
            T2= XVA+XVA
               FQ(3-1)=(T2*FQ(3)+ET)*RMR(3)
               FQ(2-1)=(T2*FQ(2)+ET)*RMR(2)
               FQ(1-1)=(T2*FQ(1)+ET)*RMR(1)

            FQF= FQZ*SQRT(X41)
               FQ(0)= FQ(0)*FQF
            FQF= FQF*RHO
               FQ(1)= FQ(1)*FQF
            FQF= FQF*RHO
               FQ(2)= FQ(2)*FQF
            FQF= FQF*RHO
               FQ(3)= FQ(3)*FQF
         ELSE
            XIN= ONE/XVA
            FQ(0)= FQZ*SQRT(PI4*XIN*X41)
            ROX= RHO*XIN
            FQF= PT5*ROX
               FQ(1)= FQ(0)*FQF
            FQF= FQF+ROX
               FQ(2)= FQ(1)*FQF
            FQF= FQF+ROX
               FQ(3)= FQ(2)*FQF
         ENDIF

         XMD1= TX21(I)
         Y01 = TY01(I)
         DP00= D00P(I)
         WORK(1,1)= DP00
         WORK(1,2)= Y01
         WORK(1,3)= XMD1
         DO J=1,3
            WORK(2,J)= WORK(1,J)*PQR
            WORK(3,J)= WORK(1,J)*PQS
         ENDDO

         FQ0(1)= FQ0(1)+FQ(0)*WORK(1,1)
         FQ0(2)= FQ0(2)+FQ(0)*WORK(1,2)
         DO J=1,3
            FQ1(1,J)= FQ1(1,J)+FQ(1)*WORK(1,J)
            FQ1(2,J)= FQ1(2,J)+FQ(1)*WORK(2,J)

            FQ2(1,J)= FQ2(1,J)+FQ(2)*WORK(1,J)
            FQ2(2,J)= FQ2(2,J)+FQ(2)*WORK(2,J)
            FQ2(3,J)= FQ2(3,J)+FQ(2)*WORK(3,J)
         ENDDO
         FQ3(1)= FQ3(1)+FQ(3)*WORK(1,3)
         FQ3(2)= FQ3(2)+FQ(3)*WORK(2,3)
         FQ3(3)= FQ3(3)+FQ(3)*WORK(3,3)
         FQ3(4)= FQ3(4)+FQ(3)*WORK(3,3)*PQR
      ENDDO

      END

!*MODULE INT2B   *DECK INTJ6
!>
!>    @brief   PPPP case
!>
!>    @details integration of a PPPP case
!>
      SUBROUTINE INTJ6_gpu(LRINT,ISMLP,ISMLQ)
      USE lrcdft, ONLY: LCFLAG, EMU, EMU2, LRFILE
      use mx_limits, only: mxgtot,mxgsh,mxg2

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! GENERATE JTYPE= 6 INTEGRALS
!
      double precision,parameter :: ZER=0.0D00
      double precision,parameter :: PT5=0.5D00
      double precision,parameter :: ONE=1.0D00
      double precision,parameter :: PI4=0.78539816339744831D00
      integer,parameter :: NTX=4
      integer,parameter :: NPF=450
      integer,parameter :: NGRD=7
      integer,parameter :: NPX=1000
      integer,parameter :: MXQT=16
      COMMON /FMTTBL/ FGRID(0:NTX,0:NPF,0:NGRD),XGRID(0:NTX,0:NPX), &
                      TMAX,RFINC(0:NGRD),RXINC, &
                      RMR(MXQT),TLGM(0:MXQT),NORD
      COMMON /FQ04  / FQ(0:4),FQ0(5),FQ1(2,9),FQ2(3,9),FQ3(4,5),FQ4(5)
!$omp threadprivate(/FQ04  /)
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS, &
                      TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2), &
                      D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2), &
                      NGANGB
!$omp threadprivate(/GEOMPQ/)
      double precision :: WORK(4,9)
      LOGICAL :: LRINT
      integer :: ISMLP(MXG2),ISMLQ



      DO J=1,5
         FQ0(J)= ZER
      ENDDO
      DO J=1,8
         FQ1(1,J)= ZER
         FQ1(2,J)= ZER
      ENDDO
         FQ1(1,9)= ZER
      DO J=1,9
         FQ2(1,J)= ZER
         FQ2(2,J)= ZER
         FQ2(3,J)= ZER
      ENDDO
      DO J=1,5
         FQ3(1,J)= ZER
         FQ3(2,J)= ZER
         FQ3(3,J)= ZER
         FQ3(4,J)= ZER
      ENDDO
      DO J=1,5
         FQ4(J)= ZER
      ENDDO

      DO I=1,NGANGB
         ISML= ISMLQ+ISMLP(I)
         IF(ISML.GE.2) CYCLE
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
!
!     FM(T) EVALUATION...DOWNWARD RECURSION FOR JTYPE.GE.5
!
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

            FQ(N)= FX
            T2= XVA+XVA
            DO M=N,1,-1
               FQ(M-1)=(T2*FQ(M)+ET)*RMR(M)
            END DO

            FQF= FQZ*SQRT(X41)
            DO M=0,N
               FQ(M)= FQ(M)*FQF
               FQF= FQF*RHO
            ENDDO
         ELSE
            XIN= ONE/XVA
            FQ(0)= FQZ*SQRT(PI4*XIN*X41)
            ROX= RHO*XIN
            FQF= PT5*ROX
            DO M=1,N
               FQ(M)= FQ(M-1)*FQF
               FQF= FQF+ROX
            ENDDO
         ENDIF

         PQT= PQR*PQS
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
         DO J=1,4
            WORK(2,J)= WORK(1,J)*PQR
            WORK(3,J)= WORK(1,J)*PQS
         ENDDO
         DO J=5,9
            WORK(2,J)= WORK(1,J)*PQR
            WORK(3,J)= WORK(1,J)*PQS
            WORK(4,J)= WORK(1,J)*PQT
         ENDDO

         DO J=1,5
            FQ0(J)= FQ0(J)+FQ(0)*WORK(1,J)
         ENDDO
         DO J=1,8
            FQ1(1,J)= FQ1(1,J)+FQ(1)*WORK(1,J)
            FQ1(2,J)= FQ1(2,J)+FQ(1)*WORK(2,J)
         ENDDO
            FQ1(1,9)= FQ1(1,9)+FQ(1)*WORK(1,9)
         DO J=1,9
            FQ2(1,J)= FQ2(1,J)+FQ(2)*WORK(1,J)
            FQ2(2,J)= FQ2(2,J)+FQ(2)*WORK(2,J)
            FQ2(3,J)= FQ2(3,J)+FQ(2)*WORK(3,J)
         ENDDO
         DO J=1,5
            FQ3(1,J)= FQ3(1,J)+FQ(3)*WORK(1,J+4)
            FQ3(2,J)= FQ3(2,J)+FQ(3)*WORK(2,J+4)
            FQ3(3,J)= FQ3(3,J)+FQ(3)*WORK(3,J+4)
            FQ3(4,J)= FQ3(4,J)+FQ(3)*WORK(4,J+4)
         ENDDO
         FQ4(1)= FQ4(1)+FQ(4)*WORK(1,9)
         FQ4(2)= FQ4(2)+FQ(4)*WORK(2,9)
         FQ4(3)= FQ4(3)+FQ(4)*WORK(3,9)
         FQ4(4)= FQ4(4)+FQ(4)*WORK(4,9)
         FQ4(5)= FQ4(5)+FQ(4)*WORK(4,9)*PQR

      ENDDO

      END




!*MODULE INT2B   *DECK INTK2
!>
!>    @brief   PSSS case
!>
!>    @details integration of a PSSS case
!>
      SUBROUTINE INTK2_1_gpu(IKL,R00,R01)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      double precision,parameter :: ZER=0.0D+00
!
! GENERATE JTYPE= 2 INTEGRALS
!
      double precision :: R00(25),R01(3,40)

      R00(1)= ZER
      R00(2)= ZER
      R01(1,1)= ZER
      R01(2,1)= ZER
      R01(3,1)= ZER

      END



!*MODULE INT2B   *DECK INTK2
!>
!>    @brief   PSSS case
!>
!>    @details integration of a PSSS case
!>
      SUBROUTINE INTK2_2_gpu(IKL,SQ,X43,FQ0,FQ1,ACY,AQX,Y03,R00,R01)
      use mx_limits, only: mxgsh,mxg2

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! GENERATE JTYPE= 2 INTEGRALS
!
     
      double precision :: X43
      double precision :: SQ(4)
      double precision :: FQ0(5),FQ1(2,9)
      double precision :: ACY,AQX,Y03
      double precision :: R00(25),R01(3,40)


      XMD2= X43 *0.5D00
      XMDT= XMD2*SQ(2)

      R00(1)= R00(1)+FQ0(1)*SQ(1)
      R00(2)= R00(2)-FQ0(1)*SQ(2)*Y03

      R01(1,1)= R01(1,1)-FQ1(1,1)*AQX *XMDT
      R01(2,1)= R01(2,1)-FQ1(1,1)*ACY *XMDT
      R01(3,1)= R01(3,1)+FQ1(2,1)     *XMDT


      END


!*MODULE INT2B   *DECK INTK3
!>
!>    @brief   PPSS case
!>
!>    @details integration of a PPSS case
!>
      SUBROUTINE INTK3_1_gpu(IKL,R00,R01,R02)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      double precision,parameter :: ZER=0.0D00
!
! GENERATE JTYPE= 3 INTEGRALS
!
      double precision :: R00(25),R01(3,40),R02(156)


      R00(1)= ZER
      R00(2)= ZER
      R00(3)= ZER
      R00(4)= ZER
      R00(5)= ZER
      DO I= 1, 4
         R01(1,I)= ZER
         R01(2,I)= ZER
         R01(3,I)= ZER
      ENDDO
      R02(1)= ZER
      R02(2)= ZER
      R02(3)= ZER
      R02(4)= ZER
      R02(5)= ZER
      R02(6)= ZER

      END



!*MODULE INT2B   *DECK INTK3
!>
!>    @brief   PPSS case
!>
!>    @details integration of a PPSS case
!>
      SUBROUTINE INTK3_2_gpu(IKL,SQ,X43,FQ0,FQ1,FQ2,ACY,ACY2,AQX,AQX2,AQXY,Y03,Y04,R00,R01,R02)
      use mx_limits, only: mxg2

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! GENERATE JTYPE= 3 INTEGRALS
!
      double precision :: X43
      double precision :: SQ(4)
      double precision :: FQ0(5),FQ1(2,9),FQ2(27)
      double precision :: ACY,ACY2,AQX,AQX2,AQXY,Y03,Y04
      double precision :: R00(25),R01(3,40),R02(156)

!       COMMON /KI3 / R00(25),R01(3,40),R02(156),R03(80),R04(15)
! !$omp threadprivate(/KI3/)

      double precision :: WORK(4)


      XMD2= X43 *0.5D+00
      XMD3= XMD2*SQ(4)
      XMDT= XMD3*XMD2
      WORK(1)= XMD2*SQ(2)
      WORK(2)= XMD2*SQ(3)
      WORK(3)=-XMD3*Y03
      WORK(4)= XMD3*Y04

      R00(1)= R00(1)+FQ0(1)     *SQ(1)
      R00(2)= R00(2)-FQ0(1)*Y03 *SQ(2)
      R00(3)= R00(3)+FQ0(1)*Y04 *SQ(3)
      R00(4)= R00(4)+FQ0(1)*XMD3
      R00(5)= R00(5)-FQ0(1)*Y03 *SQ(4)*Y04

      FQ11=-FQ1(1,1)*AQX
      FQ12=-FQ1(1,1)*ACY
      FQ13= FQ1(2,1)

      DO I= 1, 4
         R01(1,I)= R01(1,I)+FQ11*WORK(I)
         R01(2,I)= R01(2,I)+FQ12*WORK(I)
         R01(3,I)= R01(3,I)+FQ13*WORK(I)
      ENDDO

      R02(1)= R02(1)+(FQ2(1)*AQX2-FQ1(1,1))*XMDT
      R02(2)= R02(2)+(FQ2(1)*ACY2-FQ1(1,1))*XMDT
      R02(3)= R02(3)+(FQ2(3)-FQ1(1,1))     *XMDT
      R02(4)= R02(4)+ FQ2(1)*AQXY          *XMDT
      R02(5)= R02(5)- FQ2(2)*AQX           *XMDT
      R02(6)= R02(6)- FQ2(2)*ACY           *XMDT


      END



!*MODULE INT2B   *DECK INTK4
!>
!>    @brief   PSPS case
!>
!>    @details integration of a PSPS case
!>
      SUBROUTINE INTK4_1_gpu(IKL,R00,R01,R02)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      double precision,parameter :: ZER=0.0D00
!
! GENERATE JTYPE= 4 INTEGRALS
!

      double precision :: R00(25),R01(3,40),R02(156)

      DO I=1,4
         R00(I)= ZER
      ENDDO

      DO I= 1, 4
         DO J=1,3
            R01(J,I)= ZER
         ENDDO
      ENDDO

      DO I=1,6
         R02(I)= ZER
      ENDDO

      END




!*MODULE INT2B   *DECK INTK4
!>
!>    @brief   PSPS case
!>
!>    @details integration of a PSPS case
!>
      SUBROUTINE INTK4_2_gpu(IKL,SQ,X43,FQ0,FQ1,FQ2,ACY,ACY2,AQX,AQX2,AQXY,Y03,R00,R01,R02)
      use mx_limits, only: mxg2

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! GENERATE JTYPE= 4 INTEGRALS
!
      double precision :: X43
      double precision :: SQ(4)
      double precision :: FQ0(5),FQ1(2,9),FQ2(27)
      double precision :: ACY,ACY2,AQX,AQX2,AQXY,Y03
      double precision :: R00(25),R01(3,40),R02(156)

      double precision ::  WORK(4)

      

      XMD2= X43 *0.5D00
      WORK(1)= XMD2*SQ(2)
      WORK(2)=      SQ(1)
      WORK(3)= WORK(1)
      WORK(4)=-Y03 *SQ(2)

      R00(1)= R00(1)+FQ0(1)*WORK(2)
      R00(2)= R00(2)+FQ0(1)*WORK(4)
      R00(3)= R00(3)+FQ0(2)*WORK(2)
      R00(4)= R00(4)+FQ0(2)*WORK(4)

      DO I= 1, 3
         R01(1,I)= R01(1,I)-FQ1(1,I)*AQX *WORK(I)
         R01(2,I)= R01(2,I)-FQ1(1,I)*ACY *WORK(I)
         R01(3,I)= R01(3,I)+FQ1(2,I)     *WORK(I)
      ENDDO

      R01(1,4)= R01(1,4)-FQ1(1,2)*AQX *WORK(4)
      R01(2,4)= R01(2,4)-FQ1(1,2)*ACY *WORK(4)
      R01(3,4)= R01(3,4)+FQ1(2,2)     *WORK(4)

      R02(1)= R02(1)+(FQ2(1)*AQX2-FQ1(1,2))*WORK(1)
      R02(2)= R02(2)+(FQ2(1)*ACY2-FQ1(1,2))*WORK(1)
      R02(3)= R02(3)+(FQ2(3)     -FQ1(1,2))*WORK(1)
      R02(4)= R02(4)+ FQ2(1)*AQXY          *WORK(1)
      R02(5)= R02(5)- FQ2(2)*AQX           *WORK(1)
      R02(6)= R02(6)- FQ2(2)*ACY           *WORK(1)

      END



!*MODULE INT2B   *DECK INTK5
!>
!>    @brief   PPPS case
!>
!>    @details integration of a PPPS case
!>
      SUBROUTINE INTK5_1_gpu(IKL,R00,R01,R02,R03)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      double precision,parameter :: ZER=0.0D00
!
! GENERATE JTYPE= 5 INTEGRALS
!
      double precision :: R00(5,5),R01(3,40),R02(6,26),R03(10,8)

      
      DO I= 1, 2
         R00(1,I)= ZER
         R00(2,I)= ZER
         R00(3,I)= ZER
         R00(4,I)= ZER
         R00(5,I)= ZER
      ENDDO
      DO I= 1,13
         R01(1,I)= ZER
         R01(2,I)= ZER
         R01(3,I)= ZER
      ENDDO
      DO I= 1, 6
         R02(1,I)= ZER
         R02(2,I)= ZER
         R02(3,I)= ZER
         R02(4,I)= ZER
         R02(5,I)= ZER
         R02(6,I)= ZER
      ENDDO
      DO J= 1,10
         R03(J,1)= ZER
      ENDDO

      END


!*MODULE INT2B   *DECK INTK5
!>
!>    @brief   PPPS case
!>
!>    @details integration of a PPPS case
!>
      SUBROUTINE INTK5_2_gpu(IKL,SQ,X43,FQ0,FQ1,FQ2,FQ3,ACY,ACY2,AQX,AQX2,AQXY,Y03,Y04,R00,R01,R02,R03)
      use mx_limits, only: mxgsh,mxg2

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! GENERATE JTYPE= 5 INTEGRALS
!
      double precision,parameter :: F03=3.0D00
      
      double precision :: X43                    
      double precision :: SQ(4)
      double precision :: FQ0(5),FQ1(2,9),FQ2(3,9),FQ3(4,5)
      double precision :: ACY,ACY2,AQX,AQX2,AQXY,Y03,Y04
      double precision :: R00(5,5),R01(3,40),R02(6,26),R03(10,8)

      double precision :: WORK(9),FWK(6,3)

 

      XMD2= X43 *0.5D+00
      XMD3= XMD2*SQ(4)
      XMDT= XMD3*XMD2

      XMDTY=-XMDT*ACY
      XMDTX=-XMDT*AQX
      XMDTXY=XMDT*AQXY

      WORK(1)=      SQ(1)
      WORK(2)=-Y03 *SQ(2)
      WORK(3)= Y04 *SQ(3)
      WORK(4)=-Y03 *SQ(4)*Y04
      WORK(5)= XMD3
      WORK(6)= XMD2*SQ(2)
      WORK(7)= XMD2*SQ(3)
      WORK(8)=-XMD3*Y03
      WORK(9)= XMD3*Y04

      DO I= 1, 2
         DO J= 1, 5
            R00(J,I)= R00(J,I)+FQ0(I)*WORK(J)
         ENDDO
      ENDDO

      DO I= 1, 3
         FWK(1,I)=-FQ1(1,I)*AQX
         FWK(2,I)=-FQ1(1,I)*ACY
         FWK(3,I)= FQ1(2,I)
      ENDDO

      DO I= 1, 4
         R01(1,I)= R01(1,I)+FWK(1,1)*WORK(I+ 5)
         R01(2,I)= R01(2,I)+FWK(2,1)*WORK(I+ 5)
         R01(3,I)= R01(3,I)+FWK(3,1)*WORK(I+ 5)
      ENDDO

      DO I= 5, 8
         R01(1,I)= R01(1,I)+FWK(1,2)*WORK(I+ 1)
         R01(2,I)= R01(2,I)+FWK(2,2)*WORK(I+ 1)
         R01(3,I)= R01(3,I)+FWK(3,2)*WORK(I+ 1)
      ENDDO

      DO I= 9,13
         R01(1,I)= R01(1,I)+FWK(1,3)*WORK(I- 8)
         R01(2,I)= R01(2,I)+FWK(2,3)*WORK(I- 8)
         R01(3,I)= R01(3,I)+FWK(3,3)*WORK(I- 8)
      ENDDO

      DO I= 1, 2
         R02(1,I)= R02(1,I)+(FQ2(1,I)*AQX2-FQ1(1,I))*XMDT
         R02(2,I)= R02(2,I)+(FQ2(1,I)*ACY2-FQ1(1,I))*XMDT
         R02(3,I)= R02(3,I)+(FQ2(3,I)     -FQ1(1,I))*XMDT
         R02(4,I)= R02(4,I)+ FQ2(1,I)               *XMDTXY
         R02(5,I)= R02(5,I)+ FQ2(2,I)               *XMDTX
         R02(6,I)= R02(6,I)+ FQ2(2,I)               *XMDTY
      ENDDO
  
      FWK(1,3)= FQ2(1,3)*AQX2-FQ1(1,3)
      FWK(2,3)= FQ2(1,3)*ACY2-FQ1(1,3)
      FWK(3,3)= FQ2(3,3)     -FQ1(1,3)
      FWK(4,3)= FQ2(1,3)*AQXY
      FWK(5,3)=-FQ2(2,3)*AQX
      FWK(6,3)=-FQ2(2,3)*ACY

      DO I= 3, 6
         DO J= 1, 6
            R02(J,I)= R02(J,I)+FWK(J,3)*WORK(I+ 3)
         ENDDO
      ENDDO

      R03( 1,1)= R03( 1,1)+(FQ3(1,1)*AQX2-FQ2(1,3)*F03)*XMDTX
      R03( 2,1)= R03( 2,1)+(FQ3(1,1)*AQX2-FQ2(1,3)    )*XMDTY
      R03( 3,1)= R03( 3,1)+(FQ3(2,1)*AQX2-FQ2(2,3)    )*XMDT
      R03( 4,1)= R03( 4,1)+(FQ3(1,1)*ACY2-FQ2(1,3)    )*XMDTX
      R03( 5,1)= R03( 5,1)+ FQ3(2,1)                   *XMDTXY
      R03( 6,1)= R03( 6,1)+(FQ3(3,1)     -FQ2(1,3)    )*XMDTX
      R03( 7,1)= R03( 7,1)+(FQ3(1,1)*ACY2-FQ2(1,3)*F03)*XMDTY
      R03( 8,1)= R03( 8,1)+(FQ3(2,1)*ACY2-FQ2(2,3)    )*XMDT
      R03( 9,1)= R03( 9,1)+(FQ3(3,1)     -FQ2(1,3)    )*XMDTY
      R03(10,1)= R03(10,1)+(FQ3(4,1)     -FQ2(2,3)*F03)*XMDT

      END


!*MODULE INT2B   *DECK INTK6
!>
!>    @brief   PPPP case
!>
!>    @details integration of a PPPP case
!>
      SUBROUTINE INTK6_1_gpu(IKL,R00,R01,R02,R03,R04)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      double precision,parameter :: ZER=0.0D00
!
! GENERATE JTYPE= 6 INTEGRALS
!

      double precision :: R00(5,5),R01(3,40),R02(6,26),R03(10,8),R04(15)     


      DO I= 1, 5
         R00(1,I)= ZER
         R00(2,I)= ZER
         R00(3,I)= ZER
         R00(4,I)= ZER
         R00(5,I)= ZER
      ENDDO
      DO I= 1,40
         R01(1,I)= ZER
         R01(2,I)= ZER
         R01(3,I)= ZER
      ENDDO
      DO I= 1,26
         R02(1,I)= ZER
         R02(2,I)= ZER
         R02(3,I)= ZER
         R02(4,I)= ZER
         R02(5,I)= ZER
         R02(6,I)= ZER
      ENDDO
      DO I= 1, 8
         DO J= 1,10
            R03(J,I)= ZER
         ENDDO
      ENDDO
      DO J= 1,15
         R04(J)= ZER
      ENDDO

      END

!*MODULE INT2B   *DECK INTK6
!>
!>    @brief   PPPP case
!>
!>    @details integration of a PPPP case
!>
      SUBROUTINE INTK6_2_gpu(IKL,SQ,RAB,X43,FQ0,FQ1,FQ2,FQ3,FQ4,ACY,ACY2,AQX,AQX2,AQXY,Y03,Y04,R00,R01,R02,R03,R04)
      use mx_limits, only: mxgsh,mxg2

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! GENERATE JTYPE= 6 INTEGRALS
!
      double precision,parameter :: F03=3.0D00
      double precision,parameter :: F06=6.0D00
      
      double precision :: RAB,X43                  
      double precision :: SQ(4)
      double precision :: FQ0(5),FQ1(2,9),FQ2(3,9),FQ3(4,5),FQ4(5)
      double precision :: ACY,ACY2,AQX,AQX2,AQXY,Y03,Y04
      double precision :: R00(5,5),R01(3,40),R02(6,26),R03(10,8),R04(15)

      double precision :: WORK(9),FWK(10),FQW(6,9)

  
      XMD2= X43 *0.5D00
      XMD3= XMD2*SQ(4)
      XMDT= XMD3*XMD2

      XMDTY=-XMDT*ACY
      XMDTX=-XMDT*AQX
      XMDTXY=XMDT*AQXY

      WORK(1)=      SQ(1)
      WORK(2)=-Y03 *SQ(2)
      WORK(3)= Y04 *SQ(3)
      WORK(4)=-Y03 *SQ(4)*Y04
      WORK(5)= XMD3
      WORK(6)= XMD2*SQ(2)
      WORK(7)= XMD2*SQ(3)
      WORK(8)=-XMD3*Y03
      WORK(9)= XMD3*Y04

      DO J= 1, 5
         DO I= 1, 5
            R00(I,J)= R00(I,J)+FQ0(J)*WORK(I)
         ENDDO
      ENDDO

      DO I= 1, 8
         FQW(1,I)=-FQ1(1,I)*AQX
         FQW(2,I)=-FQ1(1,I)*ACY
         FQW(3,I)= FQ1(2,I)
      ENDDO
      TMP1= FQ1(1,5)*RAB +FQ1(1,8)
      TMP2= FQ1(2,5)*RAB +FQ1(2,8)
         FQW(1,9)=-TMP1*AQX
         FQW(2,9)=-TMP1*ACY
         FQW(3,9)= TMP2
      DO I= 1, 4
         R01(1,I)= R01(1,I)+FQW(1,1)*WORK(I+ 5)
         R01(2,I)= R01(2,I)+FQW(2,1)*WORK(I+ 5)
         R01(3,I)= R01(3,I)+FQW(3,1)*WORK(I+ 5)
      ENDDO
      DO I= 5, 8
         R01(1,I)= R01(1,I)+FQW(1,2)*WORK(I+ 1)
         R01(2,I)= R01(2,I)+FQW(2,2)*WORK(I+ 1)
         R01(3,I)= R01(3,I)+FQW(3,2)*WORK(I+ 1)
      ENDDO
      DO I= 9,12
         R01(1,I)= R01(1,I)+FQW(1,3)*WORK(I- 3)
         R01(2,I)= R01(2,I)+FQW(2,3)*WORK(I- 3)
         R01(3,I)= R01(3,I)+FQW(3,3)*WORK(I- 3)
      ENDDO
      DO I=13,16
         R01(1,I)= R01(1,I)+FQW(1,4)*WORK(I- 7)
         R01(2,I)= R01(2,I)+FQW(2,4)*WORK(I- 7)
         R01(3,I)= R01(3,I)+FQW(3,4)*WORK(I- 7)
      ENDDO     
      DO I=17,20
         R01(1,I)= R01(1,I)+FQW(1,5)*WORK(I-11)
         R01(2,I)= R01(2,I)+FQW(2,5)*WORK(I-11)
         R01(3,I)= R01(3,I)+FQW(3,5)*WORK(I-11)
      ENDDO
      DO I=21,25
         R01(1,I)= R01(1,I)+FQW(1,6)*WORK(I-20)
         R01(2,I)= R01(2,I)+FQW(2,6)*WORK(I-20)
         R01(3,I)= R01(3,I)+FQW(3,6)*WORK(I-20)
      ENDDO
      DO I=26,30
         R01(1,I)= R01(1,I)+FQW(1,7)*WORK(I-25)
         R01(2,I)= R01(2,I)+FQW(2,7)*WORK(I-25)
         R01(3,I)= R01(3,I)+FQW(3,7)*WORK(I-25)
      ENDDO
      DO I=31,35
         R01(1,I)= R01(1,I)+FQW(1,8)*WORK(I-30)
         R01(2,I)= R01(2,I)+FQW(2,8)*WORK(I-30)
         R01(3,I)= R01(3,I)+FQW(3,8)*WORK(I-30)
      ENDDO
      DO I=36,40
         R01(1,I)= R01(1,I)+FQW(1,9)*WORK(I-35)
         R01(2,I)= R01(2,I)+FQW(2,9)*WORK(I-35)
         R01(3,I)= R01(3,I)+FQW(3,9)*WORK(I-35)
      ENDDO

      DO I= 1, 5
         R02(1,I)= R02(1,I)+(FQ2(1,I)*AQX2-FQ1(1,I))*XMDT
         R02(2,I)= R02(2,I)+(FQ2(1,I)*ACY2-FQ1(1,I))*XMDT
         R02(3,I)= R02(3,I)+(FQ2(3,I)     -FQ1(1,I))*XMDT
         R02(4,I)= R02(4,I)+ FQ2(1,I)               *XMDTXY
         R02(5,I)= R02(5,I)+ FQ2(2,I)               *XMDTX
         R02(6,I)= R02(6,I)+ FQ2(2,I)               *XMDTY
      ENDDO
      TMP0= FQ1(1,5)*RAB +FQ1(1,8)
      TMP1= FQ2(1,5)*RAB +FQ2(1,8)
      TMP2= FQ2(2,5)*RAB +FQ2(2,8)
      TMP3= FQ2(3,5)*RAB +FQ2(3,8)
         FQW(1,5)= TMP1*AQX2-TMP0
         FQW(2,5)= TMP1*ACY2-TMP0
         FQW(3,5)= TMP3     -TMP0
         FQW(4,5)= TMP1*AQXY
         FQW(5,5)=-TMP2*AQX
         FQW(6,5)=-TMP2*ACY
      DO I= 6, 9
         FQW(1,I)= FQ2(1,I)*AQX2-FQ1(1,I)
         FQW(2,I)= FQ2(1,I)*ACY2-FQ1(1,I)
         FQW(3,I)= FQ2(3,I)     -FQ1(1,I)
         FQW(4,I)= FQ2(1,I)*AQXY
         FQW(5,I)=-FQ2(2,I)*AQX
         FQW(6,I)=-FQ2(2,I)*ACY
      ENDDO

      DO I= 6, 9
         DO J= 1, 6
            R02(J,I)= R02(J,I)+FQW(J,6)*WORK(I)
         ENDDO
      ENDDO

      DO I=10,13
         DO J= 1, 6
            R02(J,I)= R02(J,I)+FQW(J,7)*WORK(I-4)
         ENDDO
      ENDDO

      DO I=14,17
         DO J= 1, 6
            R02(J,I)= R02(J,I)+FQW(J,8)*WORK(I-8)
         ENDDO
      ENDDO
  
      DO I=18,21
         DO J= 1, 6
            R02(J,I)= R02(J,I)+FQW(J,5)*WORK(I-12)
         ENDDO
      ENDDO

      DO I=22,26
         DO J= 1, 6
            R02(J,I)= R02(J,I)+FQW(J,9)*WORK(I-21)
         ENDDO
      ENDDO

      FQ2(1,5)= FQ2(1,5)*RAB +FQ2(1,8)
      FQ2(2,5)= FQ2(2,5)*RAB +FQ2(2,8)
      FQ3(1,1)= FQ3(1,1)*RAB +FQ3(1,4)
      FQ3(2,1)= FQ3(2,1)*RAB +FQ3(2,4)
      FQ3(3,1)= FQ3(3,1)*RAB +FQ3(3,4)
      FQ3(4,1)= FQ3(4,1)*RAB +FQ3(4,4)
      DO I= 1, 4
         R03( 1,I)= R03( 1,I)+(FQ3(1,I)*AQX2-FQ2(1,I+4)*F03)*XMDTX
         R03( 2,I)= R03( 2,I)+(FQ3(1,I)*AQX2-FQ2(1,I+4)    )*XMDTY
         R03( 3,I)= R03( 3,I)+(FQ3(2,I)*AQX2-FQ2(2,I+4)    )*XMDT
         R03( 4,I)= R03( 4,I)+(FQ3(1,I)*ACY2-FQ2(1,I+4)    )*XMDTX
         R03( 5,I)= R03( 5,I)+ FQ3(2,I)                     *XMDTXY
         R03( 6,I)= R03( 6,I)+(FQ3(3,I)     -FQ2(1,I+4)    )*XMDTX
         R03( 7,I)= R03( 7,I)+(FQ3(1,I)*ACY2-FQ2(1,I+4)*F03)*XMDTY
         R03( 8,I)= R03( 8,I)+(FQ3(2,I)*ACY2-FQ2(2,I+4)    )*XMDT
         R03( 9,I)= R03( 9,I)+(FQ3(3,I)     -FQ2(1,I+4)    )*XMDTY
         R03(10,I)= R03(10,I)+(FQ3(4,I)     -FQ2(2,I+4)*F03)*XMDT
      ENDDO
      
      FWK( 1)=-(FQ3(1,5)*AQX2-FQ2(1,9)*F03)*AQX
      FWK( 2)=-(FQ3(1,5)*AQX2-FQ2(1,9)    )*ACY
      FWK( 3)=  FQ3(2,5)*AQX2-FQ2(2,9)
      FWK( 4)=-(FQ3(1,5)*ACY2-FQ2(1,9)    )*AQX
      FWK( 5)=  FQ3(2,5)                   *AQXY
      FWK( 6)=-(FQ3(3,5)     -FQ2(1,9)    )*AQX
      FWK( 7)=-(FQ3(1,5)*ACY2-FQ2(1,9)*F03)*ACY
      FWK( 8)=  FQ3(2,5)*ACY2-FQ2(2,9)
      FWK( 9)=-(FQ3(3,5)     -FQ2(1,9)    )*ACY
      FWK(10)=  FQ3(4,5)     -FQ2(2,9)*F03

      DO I= 5, 8
         DO J= 1,10
            R03(J,I)= R03(J,I)+FWK(J)*WORK(I+ 1)
         ENDDO
      ENDDO

      AQX4= AQX2*AQX2
      ACY4= ACY2*ACY2
      X2Y2= AQX2*ACY2
      Q2C2= AQX2+ACY2
      R04( 1)= R04( 1)+(FQ4(1)*AQX4-FQ3(1,5)*F06*AQX2                   &
                                   +FQ2(1,9)*F03          )*XMDT
      R04( 2)= R04( 2)+(FQ4(1)*AQX2-FQ3(1,5)*F03          )*XMDTXY
      R04( 3)= R04( 3)+(FQ4(2)*AQX2-FQ3(2,5)*F03          )*XMDTX
      R04( 4)= R04( 4)+(FQ4(1)*X2Y2-FQ3(1,5)*Q2C2+FQ2(1,9))*XMDT
      R04( 5)= R04( 5)+(FQ4(2)*AQX2-FQ3(2,5)              )*XMDTY
      R04( 6)= R04( 6)+(FQ4(3)*AQX2-FQ3(1,5)*AQX2-FQ3(3,5)              &
                                                 +FQ2(1,9))*XMDT
      R04( 7)= R04( 7)+(FQ4(1)*ACY2-FQ3(1,5)*F03          )*XMDTXY
      R04( 8)= R04( 8)+(FQ4(2)*ACY2-FQ3(2,5)              )*XMDTX
      R04( 9)= R04( 9)+(FQ4(3)     -FQ3(1,5)              )*XMDTXY
      R04(10)= R04(10)+(FQ4(4)     -FQ3(2,5)*F03          )*XMDTX
      R04(11)= R04(11)+(FQ4(1)*ACY4-FQ3(1,5)*F06*ACY2                   &
                                   +FQ2(1,9)*F03          )*XMDT
      R04(12)= R04(12)+(FQ4(2)*ACY2-FQ3(2,5)*F03          )*XMDTY
      R04(13)= R04(13)+(FQ4(3)*ACY2-FQ3(1,5)*ACY2-FQ3(3,5)              &
                                   +FQ2(1,9)              )*XMDT
      R04(14)= R04(14)+(FQ4(4)     -FQ3(2,5)*F03          )*XMDTY
      R04(15)= R04(15)+(FQ4(5)  -FQ3(3,5)*F06+FQ2(1,9)*F03)*XMDT

      END


!*MODULE INT2B   *DECK SP0S1S
!>
!>    @brief   s,p rotated axis type selection
!>
!>    @details s,p rotated axis type selection
!>
      SUBROUTINE SP0S1S_gpu(JTYPE,IKL)
      use mx_limits, only: MXGTOT,MXGSH,MXG2
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON /JMSGYH/ SQ(4)
!$omp threadprivate(/JMSGYH/)      
      COMMON /FQ04  / FQ(0:4),FQ0(5),FQ1(2,9),FQ2(27),FQ3(20),FQ4(5)
!$omp threadprivate(/FQ04  /)      
      COMMON /KI3 / R00(25),R01(120),R02(156),R03(80),R04(15)
!$omp threadprivate(/KI3/)

! ========================================================
      LOGICAL :: LRINT
      COMMON /NLRCF / LRINT
!$omp threadprivate(/NLRCF /)
      COMMON /MAXC  / CMAX(MXGTOT),CMAXA(MXGSH),CMAXB(MXGSH), &
                      CMAXC(MXGSH),CMAXD(MXGSH),ISMLP(MXG2),ISMLQ
!$omp threadprivate(/MAXC  /)
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS, &
                      TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2), &
                      D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2), &
                      NGANGB
!$omp threadprivate(/GEOMPQ/)
      COMMON /KI2 / ACY,ACY2,AQX,AQX2,AQXY,Y03,Y04
!$omp threadprivate(/KI2 /)
! ========================================================



      CALL OMP_SET_DYNAMIC(.FALSE.)

!
! FROM HERE, MCMURCHIE-DAVIDSON ALGORITHM IS USED.
!
! GENERATE (R] INTEGRALS
!
      IF(JTYPE.EQ. 1) THEN
         CALL INTJ1_gpu(LRINT,ISMLP,ISMLQ)
         R00(1)= R00(1)+FQ0(1)*SQ(1)
      ELSEIF(JTYPE.EQ. 2) THEN
         CALL INTJ2_gpu(LRINT,ISMLP,ISMLQ)
         CALL INTK2_2_gpu(IKL,SQ,X43,FQ0,FQ1, &
                          ACY,AQX,Y03, &
                          R00,R01)
      ELSEIF(JTYPE.EQ. 3) THEN
         CALL INTJ3_gpu(LRINT,ISMLP,ISMLQ)
         CALL INTK3_2_gpu(IKL,SQ,X43,FQ0,FQ1,FQ2, &
                          ACY,ACY2,AQX,AQX2,AQXY,Y03,Y04, &
                          R00,R01,R02)
      ELSEIF(JTYPE.EQ. 4) THEN
         CALL INTJ4_gpu(LRINT,ISMLP,ISMLQ)
         CALL INTK4_2_gpu(IKL,SQ,X43,FQ0,FQ1,FQ2, &
                          ACY,ACY2,AQX,AQX2,AQXY,Y03, &
                          R00,R01,R02)
      ELSEIF(JTYPE.EQ. 5) THEN
         CALL INTJ5_gpu(LRINT,ISMLP,ISMLQ)
         CALL INTK5_2_gpu(IKL,SQ,X43,FQ0,FQ1,FQ2,FQ3, &
                          ACY,ACY2,AQX,AQX2,AQXY,Y03,Y04, &
                          R00,R01,R02,R03)
      ELSEIF(JTYPE.EQ. 6) THEN
         CALL INTJ6_gpu(LRINT,ISMLP,ISMLQ)
         CALL INTK6_2_gpu(IKL,SQ,RAB,X43,FQ0,FQ1,FQ2,FQ3,FQ4, &
                          ACY,ACY2,AQX,AQX2,AQXY,Y03,Y04, &
                          R00,R01,R02,R03,R04)
      ENDIF

      END

!*MODULE INT2B   *DECK MCDV2
!>
!>    @brief   PSSS case
!>
!>    @details integration of a PSSS case
!>
      SUBROUTINE MCDV2_gpu(F,QX,QZ,R00,R01)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      double precision :: F(4,4,4,4)
      double precision :: R00(25),R01(3,40)

      F(1,1,1,1)=           +R00(    1)
      F(2,1,1,1)=+R01( 1, 1)+R00(    2)*QX
      F(3,1,1,1)=+R01( 2, 1)
      F(4,1,1,1)=+R01( 3, 1)+R00(    2)*QZ

      END



!*MODULE INT2B   *DECK MCDV3
!>
!>    @brief   PPSS case
!>
!>    @details integration of a PPSS case
!>
      SUBROUTINE MCDV3_gpu(F,QX,QZ,R00,R01,R02)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      double precision :: F(4,4,4,4)
      double precision :: R00(25),R01(3,40),R02(6,26) 


      F(1,1,1,1)=                      +R00(    1)
      F(2,1,1,1)=            R01( 1, 1)+R00(    2)*QX
      F(3,1,1,1)=            R01( 2, 1)
      F(4,1,1,1)=            R01( 3, 1)+R00(    2)*QZ

      F(1,2,1,1)=            R01( 1, 2)+R00(    3)*QX
      F(2,2,1,1)= R02( 1, 1)           +R00(    4)                      &
                          +(+R01( 1, 3)+R01( 1, 4)+R00(    5)*QX)*QX
      F(3,2,1,1)= R02( 4, 1)+R01( 2, 4)*QX
      F(4,2,1,1)= R02( 5, 1)+R01( 3, 4)*QX                              &
                          +(+R01( 1, 3)+R00(    5)*QX)*QZ

      F(1,3,1,1)=            R01( 2, 2)
      F(2,3,1,1)= R02( 4, 1)+R01( 2, 3)*QX
      F(3,3,1,1)= R02( 2, 1)           +R00(    4)
      F(4,3,1,1)= R02( 6, 1)+R01( 2, 3)*QZ

      F(1,4,1,1)=            R01( 3, 2)+R00(    3)*QZ
      F(2,4,1,1)= R02( 5, 1)+R01( 1, 4)*QZ                              &
                          +(+R01( 3, 3)+R00(    5)*QZ)*QX
      F(3,4,1,1)= R02( 6, 1)+R01( 2, 4)*QZ
      F(4,4,1,1)= R02( 3, 1)           +R00(    4)                      &
                          +(+R01( 3, 3)+R01( 3, 4)+R00(    5)*QZ)*QZ

      END

!*MODULE INT2B   *DECK MCDV4
!>
!>    @brief   PSPS case
!>
!>    @details integration of a PSPS case
!>
      SUBROUTINE MCDV4_gpu(F,QX,QZ,R00,R01,R02)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      double precision :: F(4,4,4,4)
      double precision :: R00(25),R01(3,40),R02(6,26)


      F(1,1,1,1)=                      +R00(    1)
      F(2,1,1,1)=           +R01( 1, 1)+R00(    2)*QX
      F(3,1,1,1)=           +R01( 2, 1)
      F(4,1,1,1)=           +R01( 3, 1)+R00(    2)*QZ

      F(1,1,2,1)=           -R01( 1, 2)
      F(2,1,2,1)=-R02( 1, 1)-R01( 1, 4)*QX
      F(3,1,2,1)=-R02( 4, 1)
      F(4,1,2,1)=-R02( 5, 1)-R01( 1, 4)*QZ

      F(1,1,3,1)=           -R01( 2, 2)
      F(2,1,3,1)=-R02( 4, 1)-R01( 2, 4)*QX
      F(3,1,3,1)=-R02( 2, 1)
      F(4,1,3,1)=-R02( 6, 1)-R01( 2, 4)*QZ

      F(1,1,4,1)=           -R01( 3, 2)+R00(    3)
      F(2,1,4,1)=-R02( 5, 1)-R01( 3, 4)*QX                              &
                            +R01( 1, 3)+R00(    4)*QX
      F(3,1,4,1)=-R02( 6, 1)+R01( 2, 3)
      F(4,1,4,1)=-R02( 3, 1)-R01( 3, 4)*QZ                              &
                            +R01( 3, 3)+R00(    4)*QZ

      END


!*MODULE INT2B   *DECK MCDV5
!>
!>    @brief   PPPS case
!>
!>    @details integration of a PPPS case
!>
      SUBROUTINE MCDV5_gpu(F,QX,QZ,R00,R01,R02,R03)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!jms
!     Simplified calculation of F(I,J,K,L) for cases
!     where:  I = 1..4,  J = 1..4,  K = 1..KX  and  L = 1..LX
!     using auxiliary arrays B1, B2 and B3
!jms
      integer,parameter :: KX=4,LX=1
      double precision :: R00(25),R01(3,40),R02(6,26),R03(10,8)
      double precision :: F(4,4,4,4),S01(13,3),S02(6,6)
      double precision :: B1(  5,KX,LX),B2(4,3,KX,LX),B3(  6,KX,LX)
      integer :: INS(KX,LX),INI(KX,LX),IND( 6,KX,LX),IN6(6)
      DATA INS/ 0, 0, 0, 1/
      DATA INI/ 0, 0, 1, 1/
      DATA IN6/ 1, 4, 5, 2, 6, 3/

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
      DO J=1,13
         DO I=1, 3
            S01(J,I)= R01(I,J)
         ENDDO
      ENDDO

      DO J=1, 6
         DO I=1, 6
            S02(J,I)= R02(I,J)
         ENDDO
      ENDDO
!jms
!     Auxiliary arrays to simplify the formulation of F(I,J,K,L)
!     where:  I = 1..4,  J = 1..4,  K = 1  and  L = 1
!jms
      DO J=1,6
         M= IN6(J)
         B3(  J,1,1)=+R02( M,1)
      ENDDO

      DO J=1,3
         B2(1,J,1,1)=+S01( 1,J)
         B2(2,J,1,1)=+S01( 2,J)
         B2(3,J,1,1)=+S01( 3,J)
         B2(4,J,1,1)=+S01( 4,J)
      ENDDO

      DO I=1,5
         B1(I  ,1,1)=+R00(I)
      ENDDO
!jms
!     Auxiliary arrays to simplify the formulation of F(I,J,K,L)
!     where:  I = 1..4,  J = 1..4,  K = 2  and  L = 1
!jms
      DO J=1,6
         B3(  J,2,1)=-R03( J,1)
      ENDDO

      DO J=1,3
         M= IN6(J)
         B2(1,J,2,1)=-S02( 3,M)
         B2(2,J,2,1)=-S02( 4,M)
         B2(3,J,2,1)=-S02( 5,M)
         B2(4,J,2,1)=-S02( 6,M)
      ENDDO

         J=1
      DO I=1,5
         B1(I  ,2,1)=-S01(I+ 8,J)
      ENDDO
!jms
!     Auxiliary arrays to simplify the formulation of F(I,J,K,L)
!     where:  I = 1..4,  J = 1..4,  K = 3  and  L = 1
!jms
      DO J=1,6
         K= IND(J,3,1)
         B3(  J,3,1)=-R03( K,1)
      ENDDO

      DO J=1,3
         K= IND(J,3,1)
         M= IN6(K)
         B2(1,J,3,1)=-S02( 3,M)
         B2(2,J,3,1)=-S02( 4,M)
         B2(3,J,3,1)=-S02( 5,M)
         B2(4,J,3,1)=-S02( 6,M)
      ENDDO

         J=1
         K= IND(J,3,1)
      DO I=1,5
         B1(I  ,3,1)=-S01(I+ 8,K)
      ENDDO
!jms
!     Auxiliary arrays to simplify the formulation of F(I,J,K,L)
!     where:  I = 1..4,  J = 1..4,  K = 4  and  L = 1
!jms
      DO J=1,6
         K= IND(J,4,1)
         M= IN6(J)
         B3(  J,4,1)=-R03( K,1)+R02( M,2)
      ENDDO

      DO J=1,3
         K= IND(J,4,1)
         M= IN6(K)
         B2(1,J,4,1)=-S02( 3,M)+S01( 5,J)
         B2(2,J,4,1)=-S02( 4,M)+S01( 6,J)
         B2(3,J,4,1)=-S02( 5,M)+S01( 7,J)
         B2(4,J,4,1)=-S02( 6,M)+S01( 8,J)
      ENDDO

         J=1
         K= IND(J,4,1)
      DO I=1,5
         B1(I  ,4,1)=-S01(I+ 8,K)+R00(I+ 5)
      ENDDO

      DO L=1,LX
         DO K=1,KX

            F(1,1,K,L)=+B1(  1,K,L)
            F(2,1,K,L)=+B2(1,1,K,L)+B1(  2,K,L)*QX
            F(3,1,K,L)=+B2(1,2,K,L)
            F(4,1,K,L)=+B2(1,3,K,L)+B1(  2,K,L)*QZ

            F(1,2,K,L)=+B2(2,1,K,L)+B1(  3,K,L)*QX
            F(2,2,K,L)=+B3(  1,K,L)+B1(  5,K,L)                         &
                     +(+B2(3,1,K,L)+B2(4,1,K,L)+B1(  4,K,L)*QX)*QX
            F(3,2,K,L)=+B3(  2,K,L)+B2(4,2,K,L)*QX
            F(4,2,K,L)=+B3(  3,K,L)+B2(4,3,K,L)*QX                      &
                     +(+B2(3,1,K,L)+B1(  4,K,L)*QX)*QZ

            F(1,3,K,L)=+B2(2,2,K,L)
            F(2,3,K,L)=+B3(  2,K,L)+B2(3,2,K,L)*QX
            F(3,3,K,L)=+B3(  4,K,L)+B1(  5,K,L)
            F(4,3,K,L)=+B3(  5,K,L)+B2(3,2,K,L)*QZ

            F(1,4,K,L)=+B2(2,3,K,L)+B1(  3,K,L)*QZ
            F(2,4,K,L)=+B3(  3,K,L)+B2(4,1,K,L)*QZ                      &
                     +(+B2(3,3,K,L)+B1(  4,K,L)*QZ)*QX
            F(3,4,K,L)=+B3(  5,K,L)+B2(4,2,K,L)*QZ
            F(4,4,K,L)=+B3(  6,K,L)+B1(  5,K,L)                         &
                     +(+B2(3,3,K,L)+B2(4,3,K,L)+B1(  4,K,L)*QZ)*QZ

         ENDDO
      ENDDO

      END


!*MODULE INT2B   *DECK MCDV6
!>
!>    @brief   PPPP case
!>
!>    @details integration of a PPPP case
!>
      SUBROUTINE MCDV6_gpu(F,QX,QZ,R00,R01,R02,R03,R04)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      double precision :: F(4,4,4,4)
      double precision :: R00(25),R01(3,40),R02(6,26),R03(10,8),R04(15)
      double precision :: S01(40,3),S02(26,6),S03(8,10)
      LOGICAL :: LSYM06
!jms
!     Simplified calculation of F(I,J,K,L) for cases
!     where:  I = 1..4,  J = 1..4,  K = 1..KX  and  L = 1..LX
!     using auxiliary arrays B1, B2 and B3
!jms
      PARAMETER (KX=4)
      PARAMETER (LX=4)
      DIMENSION       B1(  5,KX,LX),B2(4,3,KX,LX),B3(  6,KX,LX)

      DIMENSION INS(KX,LX),INI(KX,LX),IND( 6,KX,LX)
      DIMENSION IN6(6)
      DATA INS/ 0, 0, 0, 1,  0, 0, 0, 1,  0, 0, 1, 2,  1, 1, 2, 3/
      DATA INI/ 0, 0, 1, 1,  0, 0, 1, 1,  1, 1, 2, 2,  1, 1, 2, 2/
      DATA IN6/ 1, 4, 5, 2, 6, 3/

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
      DO J=1,40
         DO I=1, 3
            S01(J,I)= R01(I,J)
         ENDDO
      ENDDO

      DO J=1,26
         DO I=1, 6
            S02(J,I)= R02(I,J)
         ENDDO
      ENDDO

      DO J=1, 8
         DO I=1,10
            S03(J,I)= R03(I,J)
         ENDDO
      ENDDO

!jms
!     Auxiliary arrays to simplify the formulation of F(I,J,K,L)
!     where:  I = 1..4,  J = 1..4,  K = 1  and  L = 1
!jms
      DO J=1,6
         M= IN6(J)
         B3(  J,1,1)=+R02( M,1)
      ENDDO

      DO J=1,3
         B2(1,J,1,1)=+S01( 1,J)
         B2(2,J,1,1)=+S01( 2,J)
         B2(3,J,1,1)=+S01( 3,J)
         B2(4,J,1,1)=+S01( 4,J)
      ENDDO

      DO I=1,5
         B1(I  ,1,1)=+R00(I)
      ENDDO
!jms
!     Auxiliary arrays to simplify the formulation of F(I,J,K,L)
!     where:  I = 1..4,  J = 1..4,  K = 2  and  L = 1
!jms
      DO J=1,6
         B3(  J,2,1)=-R03( J,2)
      ENDDO

      DO J=1,3
         M= IN6(J)
         B2(1,J,2,1)=-S02( 6,M)
         B2(2,J,2,1)=-S02( 7,M)
         B2(3,J,2,1)=-S02( 8,M)
         B2(4,J,2,1)=-S02( 9,M)
      ENDDO

      J=1
      DO I=1,5
         B1(I  ,2,1)=-S01(I+20,J)
      ENDDO
!jms
!     Auxiliary arrays to simplify the formulation of F(I,J,K,L)
!     where:  I = 1..4,  J = 1..4,  K = 3  and  L = 1
!jms
      DO J=1,6
         K= IND(J,3,1)
         B3(  J,3,1)=-R03( K,2)
      ENDDO

      DO J=1,3
         K= IND(J,3,1)
         M= IN6(K)
         B2(1,J,3,1)=-S02( 6,M)
         B2(2,J,3,1)=-S02( 7,M)
         B2(3,J,3,1)=-S02( 8,M)
         B2(4,J,3,1)=-S02( 9,M)
      ENDDO

      J=1
      K= IND(J,3,1)
      DO I=1,5
         B1(I  ,3,1)=-S01(I+20,K)
      ENDDO
!jms
!     Auxiliary arrays to simplify the formulation of F(I,J,K,L)
!     where:  I = 1..4,  J = 1..4,  K = 4  and  L = 1
!jms
      DO J=1,6
         K= IND(J,4,1)
         M= IN6(J)
         B3(  J,4,1)=-R03( K,2)+R02( M,2)
      ENDDO

      DO J=1,3
         K= IND(J,4,1)
         M= IN6(K)
         B2(1,J,4,1)=-S02( 6,M)+S01( 5,J)
         B2(2,J,4,1)=-S02( 7,M)+S01( 6,J)
         B2(3,J,4,1)=-S02( 8,M)+S01( 7,J)
         B2(4,J,4,1)=-S02( 9,M)+S01( 8,J)
      ENDDO

      J=1
      K= IND(J,4,1)
      DO I=1,5
         B1(I  ,4,1)=-S01(I+20,K)+R00(I+ 5)
      ENDDO
!jms
!     Auxiliary arrays to simplify the formulation of F(I,J,K,L)
!     where:  I = 1..4,  J = 1..4,  K = 1  and  L = 2
!jms
      DO J=1,6
         B3(  J,1,2)=-R03( J,3)
      ENDDO

      DO J=1,3
         M= IN6(J)
         B2(1,J,1,2)=-S02(10,M)
         B2(2,J,1,2)=-S02(11,M)
         B2(3,J,1,2)=-S02(12,M)
         B2(4,J,1,2)=-S02(13,M)
      ENDDO

      J=1
      DO I=1,5
         B1(I  ,1,2)=-S01(I+25,J)
      ENDDO
!jms
!     Auxiliary arrays to simplify the formulation of F(I,J,K,L)
!     where:  I = 1..4,  J = 1..4,  K = 2  and  L = 2
!jms
      DO J=1,6
         M= IN6(J)
         B3(  J,2,2)=+R04(   J)+R02( M,5)
      ENDDO

      DO J=1,3
         B2(1,J,2,2)=+S03( 5,J)+S01(17,J)
         B2(2,J,2,2)=+S03( 6,J)+S01(18,J)
         B2(3,J,2,2)=+S03( 7,J)+S01(19,J)
         B2(4,J,2,2)=+S03( 8,J)+S01(20,J)
      ENDDO

      J=1
      DO I=1,5
         B1(I  ,2,2)=+S02(I+21,J)+R00(I+20)
      ENDDO
!jms
!     Auxiliary arrays to simplify the formulation of F(I,J,K,L)
!     where:  I = 1..4,  J = 1..4,  K = 3  and  L = 2
!jms
      DO J=1,6
         K= IND(J,3,2)
         B3(  J,3,2)=+R04(   K)
      ENDDO

      DO J=1,3
         K= IND(J,3,2)
         B2(1,J,3,2)=+S03( 5,K)
         B2(2,J,3,2)=+S03( 6,K)
         B2(3,J,3,2)=+S03( 7,K)
         B2(4,J,3,2)=+S03( 8,K)
      ENDDO

         J=1
         K= IND(J,3,2)
         M= IN6(K)
      DO I=1,5
         B1(I  ,3,2)=+S02(I+21,M)
      ENDDO
!jms
!     Auxiliary arrays to simplify the formulation of F(I,J,K,L)
!     where:  I = 1..4,  J = 1..4,  K = 4  and  L = 2
!jms
      DO J=1,6
         K= IND(J,4,2)
         B3(  J,4,2)=+R04(   K)-R03( J,4)
      ENDDO

      DO J=1,3
         K= IND(J,4,2)
         M= IN6(J)
         B2(1,J,4,2)=+S03( 5,K)-S02(14,M)
         B2(2,J,4,2)=+S03( 6,K)-S02(15,M)
         B2(3,J,4,2)=+S03( 7,K)-S02(16,M)
         B2(4,J,4,2)=+S03( 8,K)-S02(17,M)
      ENDDO

         J=1
         K= IND(J,4,2)
         M= IN6(K)
      DO I=1,5
         B1(I  ,4,2)=+S02(I+21,M)-S01(I+30,J)
      ENDDO
!jms
!     Auxiliary arrays to simplify the formulation of F(I,J,K,L)
!     where:  I = 1..4,  J = 1..4,  K = 1  and  L = 3
!jms
      DO J=1,6
         K= IND(J,1,3)
         B3(  J,1,3)=-R03( K,3)
      ENDDO

      DO J=1,3
         K= IND(J,1,3)
         M= IN6(K)
         B2(1,J,1,3)=-S02(10,M)
         B2(2,J,1,3)=-S02(11,M)
         B2(3,J,1,3)=-S02(12,M)
         B2(4,J,1,3)=-S02(13,M)
      ENDDO

         J=1
         K= IND(J,1,3)
      DO I=1,5
         B1(I  ,1,3)=-S01(I+25,K)
      ENDDO
!jms
!     Auxiliary arrays to simplify the formulation of F(I,J,K,L)
!     where:  I = 1..4,  J = 1..4,  K = 2  and  L = 3
!
!     next statements commented because F(I,J,2,3)
!     will be obtained at the end from  F(I,J,3,2)
!jms
!     DO 323 J=1,6
! 323 B3(  J,2,3)= B3(  J,3,2)
!     DO 322 J=1,3
!        DO 322 I=1,4
! 322 B2(I,J,2,3)= B2(I,J,3,2)
!        DO 321 I=1,5
! 321 B1(I  ,2,3)= B1(I  ,3,2)
!jms
!     Auxiliary arrays to simplify the formulation of F(I,J,K,L)
!     where:  I = 1..4,  J = 1..4,  K = 3  and  L = 3
!jms
      DO J=1,6
         K= IND(J,3,3)
         M= IN6(J)
         B3(  J,3,3)=+R04(   K)+R02( M,5)
      ENDDO

      DO J=1,3
         K= IND(J,3,3)
         B2(1,J,3,3)=+S03( 5,K)+S01(17,J)
         B2(2,J,3,3)=+S03( 6,K)+S01(18,J)
         B2(3,J,3,3)=+S03( 7,K)+S01(19,J)
         B2(4,J,3,3)=+S03( 8,K)+S01(20,J)
      ENDDO

         J=1
         K= IND(J,3,3)
         M= IN6(K)
      DO I=1,5
         B1(I  ,3,3)=+S02(I+21,M)+R00(I+20)
      ENDDO
!jms
!     Auxiliary arrays to simplify the formulation of F(I,J,K,L)
!     where:  I = 1..4,  J = 1..4,  K = 4  and  L = 3
!jms
      B3(  1,4,3)=+R04(    5)-R03( 2, 4)

      B3(  2,4,3)=+R04(    8)-R03( 4, 4)
      B3(  3,4,3)=+R04(    9)-R03( 5, 4)

      B3(  4,4,3)=+R04(   12)-R03( 7, 4)
      B3(  5,4,3)=+R04(   13)-R03( 8, 4)
      B3(  6,4,3)=+R04(   14)-R03( 9, 4)

      B2(1,1,4,3)=+S03( 5, 5)-S02(14, 4)
      B2(2,1,4,3)=+S03( 6, 5)-S02(15, 4)
      B2(3,1,4,3)=+S03( 7, 5)-S02(16, 4)
      B2(4,1,4,3)=+S03( 8, 5)-S02(17, 4)

      B2(1,2,4,3)=+S03( 5, 8)-S02(14, 2)
      B2(2,2,4,3)=+S03( 6, 8)-S02(15, 2)
      B2(3,2,4,3)=+S03( 7, 8)-S02(16, 2)
      B2(4,2,4,3)=+S03( 8, 8)-S02(17, 2)
      B2(1,3,4,3)=+S03( 5, 9)-S02(14, 6)
      B2(2,3,4,3)=+S03( 6, 9)-S02(15, 6)
      B2(3,3,4,3)=+S03( 7, 9)-S02(16, 6)
      B2(4,3,4,3)=+S03( 8, 9)-S02(17, 6)

      DO I=1,5
         B1(I  ,4,3)=+S02(I+21,6)-S01(I+30,2)
      ENDDO
!jms
!     Auxiliary arrays to simplify the formulation of F(I,J,K,L)
!     where:  I = 1..4,  J = 1..4,  K = 1  and  L = 4
!jms
      DO J=1,6
         K= IND(J,1,4)
         M= IN6(J)
         B3(  J,1,4)=-R03( K,3)+R02( M,3)
      ENDDO

      DO J=1,3
         K= IND(J,1,4)
         M= IN6(K)
         B2(1,J,1,4)=-S02(10,M)+S01( 9,J)
         B2(2,J,1,4)=-S02(11,M)+S01(10,J)
         B2(3,J,1,4)=-S02(12,M)+S01(11,J)
         B2(4,J,1,4)=-S02(13,M)+S01(12,J)
      ENDDO

         J=1
         K= IND(J,1,4)
      DO I=1,5
         B1(I  ,1,4)=-S01(I+25,K)+R00(I+10)
      ENDDO
!jms
!     Auxiliary arrays to simplify the formulation of F(I,J,K,L)
!     where:  I = 1..4,  J = 1..4,  K = 2  and  L = 4
!jms
      DO J=1,6
         K= IND(J,2,4)
         B3(  J,2,4)=+R04(   K)-R03( J,1)
      ENDDO

      DO J=1,3
         K= IND(J,2,4)
         M= IN6(J)
         B2(1,J,2,4)=+S03( 5,K)-S02(18,M)
         B2(2,J,2,4)=+S03( 6,K)-S02(19,M)
         B2(3,J,2,4)=+S03( 7,K)-S02(20,M)
         B2(4,J,2,4)=+S03( 8,K)-S02(21,M)
      ENDDO

         J=1
         K= IND(J,2,4)
         M= IN6(K)
      DO I=1,5
         B1(I  ,2,4)=+S02(I+21,M)-S01(I+35,J)
      ENDDO
!jms
!     Auxiliary arrays to simplify the formulation of F(I,J,K,L)
!     where:  I = 1..4,  J = 1..4,  K = 3  and  L = 4
!jms
      B3(  1,3,4)=+R04(    5)-R03( 2, 1)

      B3(  2,3,4)=+R04(    8)-R03( 4, 1)
      B3(  3,3,4)=+R04(    9)-R03( 5, 1)

      B3(  4,3,4)=+R04(   12)-R03( 7, 1)
      B3(  5,3,4)=+R04(   13)-R03( 8, 1)
      B3(  6,3,4)=+R04(   14)-R03( 9, 1)

      B2(1,1,3,4)=+S03( 5, 5)-S02(18, 4)
      B2(2,1,3,4)=+S03( 6, 5)-S02(19, 4)
      B2(3,1,3,4)=+S03( 7, 5)-S02(20, 4)
      B2(4,1,3,4)=+S03( 8, 5)-S02(21, 4)

      B2(1,2,3,4)=+S03( 5, 8)-S02(18, 2)
      B2(2,2,3,4)=+S03( 6, 8)-S02(19, 2)
      B2(3,2,3,4)=+S03( 7, 8)-S02(20, 2)
      B2(4,2,3,4)=+S03( 8, 8)-S02(21, 2)
      B2(1,3,3,4)=+S03( 5, 9)-S02(18, 6)
      B2(2,3,3,4)=+S03( 6, 9)-S02(19, 6)
      B2(3,3,3,4)=+S03( 7, 9)-S02(20, 6)
      B2(4,3,3,4)=+S03( 8, 9)-S02(21, 6)

      DO I=1,5
         B1(I  ,3,4)=+S02(I+21,6)-S01(I+35,2)
      ENDDO
!jms
!     Auxiliary arrays to simplify the formulation of F(I,J,K,L)
!     where:  I = 1..4,  J = 1..4,  K = 4  and  L = 4
!jms
      B3(  1,4,4)=R04(    6)-S03( 1, 3)-S03( 4, 3)+S02( 4, 1)+S02( 5, 1)

      B3(  2,4,4)=R04(    9)-S03( 1, 5)-S03( 4, 5)+S02( 4, 4)+S02( 5, 4)
      B3(  3,4,4)=R04(   10)-S03( 1, 6)-S03( 4, 6)+S02( 4, 5)+S02( 5, 5)

      B3(  4,4,4)=R04(   13)-S03( 1, 8)-S03( 4, 8)+S02( 4, 2)+S02( 5, 2)
      B3(  5,4,4)=R04(   14)-S03( 1, 9)-S03( 4, 9)+S02( 4, 6)+S02( 5, 6)
      B3(  6,4,4)=R04(   15)-S03( 1,10)-S03( 4,10)+S02( 4, 3)+S02( 5, 3)

      B2(1,1,4,4)=S03( 5, 6)-S02(14, 5)-S02(18, 5)+S01(13, 1)+S01(17, 1)
      B2(2,1,4,4)=S03( 6, 6)-S02(15, 5)-S02(19, 5)+S01(14, 1)+S01(18, 1)
      B2(3,1,4,4)=S03( 7, 6)-S02(16, 5)-S02(20, 5)+S01(15, 1)+S01(19, 1)
      B2(4,1,4,4)=S03( 8, 6)-S02(17, 5)-S02(21, 5)+S01(16, 1)+S01(20, 1)

      B2(1,2,4,4)=S03( 5, 9)-S02(14, 6)-S02(18, 6)+S01(13, 2)+S01(17, 2)
      B2(2,2,4,4)=S03( 6, 9)-S02(15, 6)-S02(19, 6)+S01(14, 2)+S01(18, 2)
      B2(3,2,4,4)=S03( 7, 9)-S02(16, 6)-S02(20, 6)+S01(15, 2)+S01(19, 2)
      B2(4,2,4,4)=S03( 8, 9)-S02(17, 6)-S02(21, 6)+S01(16, 2)+S01(20, 2)
      B2(1,3,4,4)=S03( 5,10)-S02(14, 3)-S02(18, 3)+S01(13, 3)+S01(17, 3)
      B2(2,3,4,4)=S03( 6,10)-S02(15, 3)-S02(19, 3)+S01(14, 3)+S01(18, 3)
      B2(3,3,4,4)=S03( 7,10)-S02(16, 3)-S02(20, 3)+S01(15, 3)+S01(19, 3)
      B2(4,3,4,4)=S03( 8,10)-S02(17, 3)-S02(21, 3)+S01(16, 3)+S01(20, 3)

      DO I=1,5
         B1(I  ,4,4)=+S02(I+21,3)-S01(I+30,3)-S01(I+35,3)               &
                                 +R00(I+15)  +R00(I+20)
      ENDDO

      LSYM06=(KX.EQ.4 .AND. LX.EQ.4)
      DO L=1,LX
         DO K=1,KX

            IF(LSYM06) THEN
               IF(K.EQ.2 .AND. L.EQ.3) CYCLE
            ENDIF

            F(1,1,K,L)=+B1(  1,K,L)
            F(2,1,K,L)=+B2(1,1,K,L)+B1(  2,K,L)*QX
            F(3,1,K,L)=+B2(1,2,K,L)
            F(4,1,K,L)=+B2(1,3,K,L)+B1(  2,K,L)*QZ

            F(1,2,K,L)=+B2(2,1,K,L)+B1(  3,K,L)*QX
            F(2,2,K,L)=+B3(  1,K,L)+B1(  5,K,L)                         &
                     +(+B2(3,1,K,L)+B2(4,1,K,L)+B1(  4,K,L)*QX)*QX
            F(3,2,K,L)=+B3(  2,K,L)+B2(4,2,K,L)*QX
            F(4,2,K,L)=+B3(  3,K,L)+B2(4,3,K,L)*QX                      &
                     +(+B2(3,1,K,L)+B1(  4,K,L)*QX)*QZ

            F(1,3,K,L)=+B2(2,2,K,L)
            F(2,3,K,L)=+B3(  2,K,L)+B2(3,2,K,L)*QX
            F(3,3,K,L)=+B3(  4,K,L)+B1(  5,K,L)
            F(4,3,K,L)=+B3(  5,K,L)+B2(3,2,K,L)*QZ

            F(1,4,K,L)=+B2(2,3,K,L)+B1(  3,K,L)*QZ
            F(2,4,K,L)=+B3(  3,K,L)+B2(4,1,K,L)*QZ                      &
                     +(+B2(3,3,K,L)+B1(  4,K,L)*QZ)*QX
            F(3,4,K,L)=+B3(  5,K,L)+B2(4,2,K,L)*QZ
            F(4,4,K,L)=+B3(  6,K,L)+B1(  5,K,L)                         &
                     +(+B2(3,3,K,L)+B2(4,3,K,L)+B1(  4,K,L)*QZ)*QZ

         ENDDO
      ENDDO

      IF(LSYM06) THEN
         DO J=1,4
            DO I=1,4
               F(I,J,2,3)= F(I,J,3,2)
            ENDDO
         ENDDO
      ENDIF

      END



!*MODULE INT2B   *DECK R30S1S
!>
!>    @brief   rotate up to 256 s,p integrals to space fixed axes
!>
!>    @details rotate up to 256 s,p integrals to space fixed axes
!>             INCOMING AND OUTGOING INTEGRALS IN F, while P(1,1),...
!>             ARE DIRECTION COSINES OF SPACE FIXED AXES WRT AXES AT P
!>
      SUBROUTINE R30S1S_J2_gpu(JTYPE,F,P)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      double precision :: F(0:3,0:3,0:3,0:3),P(3,3)
      double precision :: T(3)

      T(1)= F(1,0,0,0)
      T(2)= F(2,0,0,0)
      T(3)= F(3,0,0,0)
      F(1,0,0,0)= T(1)*P(1,1)+T(2)*P(2,1)+T(3)*P(3,1)
      F(2,0,0,0)= T(1)*P(1,2)+T(2)*P(2,2)+T(3)*P(3,2)
      F(3,0,0,0)= T(1)*P(1,3)+T(2)*P(2,3)+T(3)*P(3,3)

      END


      SUBROUTINE R30S1S_J3_gpu(JTYPE,F,P)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      double precision :: F(0:3,0:3,0:3,0:3),P(3,3)
      double precision :: T(3)

      DO I=0,3
         T(1)= F(I,1,0,0)
         T(2)= F(I,2,0,0)
         T(3)= F(I,3,0,0)
         F(I,1,0,0)= T(1)*P(1,1)+T(2)*P(2,1)+T(3)*P(3,1)
         F(I,2,0,0)= T(1)*P(1,2)+T(2)*P(2,2)+T(3)*P(3,2)
         F(I,3,0,0)= T(1)*P(1,3)+T(2)*P(2,3)+T(3)*P(3,3)
      ENDDO
      DO J=0,3
         T(1)= F(1,J,0,0)
         T(2)= F(2,J,0,0)
         T(3)= F(3,J,0,0)
         F(1,J,0,0)= T(1)*P(1,1)+T(2)*P(2,1)+T(3)*P(3,1)
         F(2,J,0,0)= T(1)*P(1,2)+T(2)*P(2,2)+T(3)*P(3,2)
         F(3,J,0,0)= T(1)*P(1,3)+T(2)*P(2,3)+T(3)*P(3,3)
      ENDDO

      END


      SUBROUTINE R30S1S_J4_gpu(JTYPE,F,P)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      double precision :: F(0:3,0:3,0:3,0:3),P(3,3)
      double precision :: T(3)

      
      DO I=0,3
         T(1)= F(I,0,1,0)
         T(2)= F(I,0,2,0)
         T(3)= F(I,0,3,0)
         F(I,0,1,0)= T(1)*P(1,1)+T(2)*P(2,1)+T(3)*P(3,1)
         F(I,0,2,0)= T(1)*P(1,2)+T(2)*P(2,2)+T(3)*P(3,2)
         F(I,0,3,0)= T(1)*P(1,3)+T(2)*P(2,3)+T(3)*P(3,3)
      ENDDO
      DO K=0,3
         T(1)= F(1,0,K,0)
         T(2)= F(2,0,K,0)
         T(3)= F(3,0,K,0)
         F(1,0,K,0)= T(1)*P(1,1)+T(2)*P(2,1)+T(3)*P(3,1)
         F(2,0,K,0)= T(1)*P(1,2)+T(2)*P(2,2)+T(3)*P(3,2)
         F(3,0,K,0)= T(1)*P(1,3)+T(2)*P(2,3)+T(3)*P(3,3)
      ENDDO

      END

      
      SUBROUTINE R30S1S_J5_gpu(JTYPE,F,P)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      double precision :: F(0:3,0:3,0:3,0:3),P(3,3)
      double precision :: T(3)

      
      DO J=0,3
         DO I=0,3
            T(1)= F(I,J,1,0)
            T(2)= F(I,J,2,0)
            T(3)= F(I,J,3,0)
            F(I,J,1,0)= T(1)*P(1,1)+T(2)*P(2,1)+T(3)*P(3,1)
            F(I,J,2,0)= T(1)*P(1,2)+T(2)*P(2,2)+T(3)*P(3,2)
            F(I,J,3,0)= T(1)*P(1,3)+T(2)*P(2,3)+T(3)*P(3,3)
         ENDDO
      ENDDO
      DO K=0,3
         DO I=0,3
            T(1)= F(I,1,K,0)
            T(2)= F(I,2,K,0)
            T(3)= F(I,3,K,0)
            F(I,1,K,0)= T(1)*P(1,1)+T(2)*P(2,1)+T(3)*P(3,1)
            F(I,2,K,0)= T(1)*P(1,2)+T(2)*P(2,2)+T(3)*P(3,2)
            F(I,3,K,0)= T(1)*P(1,3)+T(2)*P(2,3)+T(3)*P(3,3)
         ENDDO
      ENDDO
      DO K=0,3
         DO J=0,3
            T(1)= F(1,J,K,0)
            T(2)= F(2,J,K,0)
            T(3)= F(3,J,K,0)
            F(1,J,K,0)= T(1)*P(1,1)+T(2)*P(2,1)+T(3)*P(3,1)
            F(2,J,K,0)= T(1)*P(1,2)+T(2)*P(2,2)+T(3)*P(3,2)
            F(3,J,K,0)= T(1)*P(1,3)+T(2)*P(2,3)+T(3)*P(3,3)
         ENDDO
      ENDDO

      END

!*MODULE INT2B   *DECK R30S1S
!>
!>    @brief   rotate up to 256 s,p integrals to space fixed axes
!>
!>    @details rotate up to 256 s,p integrals to space fixed axes
!>             INCOMING AND OUTGOING INTEGRALS IN F, while P(1,1),...
!>             ARE DIRECTION COSINES OF SPACE FIXED AXES WRT AXES AT P
!>
      SUBROUTINE R30S1S_J6_gpu(JTYPE,F,P)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      double precision :: F(0:3,0:3,0:3,0:3),P(3,3)
      double precision :: T(3)

            
      DO K=0,3
         DO J=0,3
            DO I=0,3
               T(1)= F(I,J,K,1)
               T(2)= F(I,J,K,2)
               T(3)= F(I,J,K,3)
               F(I,J,K,1)= T(1)*P(1,1)+T(2)*P(2,1)+T(3)*P(3,1)
               F(I,J,K,2)= T(1)*P(1,2)+T(2)*P(2,2)+T(3)*P(3,2)
               F(I,J,K,3)= T(1)*P(1,3)+T(2)*P(2,3)+T(3)*P(3,3)
            ENDDO
         ENDDO
      ENDDO
      DO L=0,3
         DO J=0,3
            DO I=0,3
               T(1)= F(I,J,1,L)
               T(2)= F(I,J,2,L)
               T(3)= F(I,J,3,L)
               F(I,J,1,L)= T(1)*P(1,1)+T(2)*P(2,1)+T(3)*P(3,1)
               F(I,J,2,L)= T(1)*P(1,2)+T(2)*P(2,2)+T(3)*P(3,2)
               F(I,J,3,L)= T(1)*P(1,3)+T(2)*P(2,3)+T(3)*P(3,3)
            ENDDO
         ENDDO
      ENDDO
      DO L=0,3
         DO K=0,3
            DO I=0,3
               T(1)= F(I,1,K,L)
               T(2)= F(I,2,K,L)
               T(3)= F(I,3,K,L)
               F(I,1,K,L)= T(1)*P(1,1)+T(2)*P(2,1)+T(3)*P(3,1)
               F(I,2,K,L)= T(1)*P(1,2)+T(2)*P(2,2)+T(3)*P(3,2)
               F(I,3,K,L)= T(1)*P(1,3)+T(2)*P(2,3)+T(3)*P(3,3)
            ENDDO
         ENDDO
      ENDDO
      DO L=0,3
         DO K=0,3
            DO J=0,3
               T(1)= F(1,J,K,L)
               T(2)= F(2,J,K,L)
               T(3)= F(3,J,K,L)
               F(1,J,K,L)= T(1)*P(1,1)+T(2)*P(2,1)+T(3)*P(3,1)
               F(2,J,K,L)= T(1)*P(1,2)+T(2)*P(2,2)+T(3)*P(3,2)
               F(3,J,K,L)= T(1)*P(1,3)+T(2)*P(2,3)+T(3)*P(3,3)
            ENDDO
         ENDDO
      ENDDO


      END
