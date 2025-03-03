SUBROUTINE DFZERO (F, B, C, R, RE, AE, IFLAG)
IMPLICIT NONE
REAL(8)  :: A,ACBS,ACMB,AE,AW,B,C,CMB,D1MACH,ER,F,FA,FB,FC,FX,FZ,P,Q,R,RE,RW,T,TOL,Z
INTEGER  :: IC,IFLAG,KOUNT


   ER = 2.0_8 * D1MACH

! Initialize.

   Z = R
   IF (R .LE. MIN(B,C)  .OR.  R .GE. MAX(B,C)) Z = C
   RW = MAX(RE,ER)
   AW = MAX(AE,0._8)
   IC = 0
   T = Z
   FZ = F(T)
   FC = FZ
   T = B
   FB = F(T)
   KOUNT = 2
   IF (SIGN(1.0_8,FZ) .EQ. SIGN(1.0_8,FB)) GO TO 1
   C = Z
   GO TO 2
1  IF (Z .EQ. C) GO TO 2
   T = C
   FC = F(T)
   KOUNT = 3
   IF (SIGN(1.0_8,FZ) .EQ. SIGN(1.0_8,FC)) GO TO 2
   B = Z
   FB = FZ
2  A = C
   FA = FC
   ACBS = ABS(B-C)
   FX = MAX(ABS(FB),ABS(FC))

3  IF (ABS(FC) .GE. ABS(FB)) GO TO 4

! Perform interchange.

   A = B
   FA = FB
   B = C
   FB = FC
   C = A
   FC = FA

4  CMB = 0.5_8*(C-B)
   ACMB = ABS(CMB)
   TOL = RW*ABS(B) + AW

!   Test stopping criterion and function count.

   IF (ACMB .LE. TOL) GO TO 10
   IF (FB .EQ. 0._8) GO TO 11
   IF (KOUNT .GE. 500) GO TO 14

! Calculate new iterate implicitly as B+P/Q, where we arrange
! P .GE. 0.  The implicit form is used to prevent overflow.

   P = (B-A)*FB
   Q = FA - FB
   IF (P .GE. 0._8) GO TO 5
   P = -P
   Q = -Q

! Update A and check for satisfactory reduction in the size of the
! bracketing interval.  If not, perform bisection.

5  A = B
   FA = FB
   IC = IC + 1
   IF (IC .LT. 4) GO TO 6
   IF (8.0_8*ACMB .GE. ACBS) GO TO 8
   IC = 0
   ACBS = ACMB

! Test for too small a change.

6  IF (P .GT. ABS(Q)*TOL) GO TO 7

! Increment by TOLerance.

   B = B + SIGN(TOL,CMB)
   GO TO 9

! Root ought to be between B and (C+B)/2.

7  IF (P .GE. CMB*Q) GO TO 8

! Use secant rule.

   B = B + P/Q
   GO TO 9

! Use bisection (C+B)/2.

8  B = B + CMB

! Have completed computation for new iterate B.

9  T = B
   FB = F(T)
   KOUNT = KOUNT + 1

! Decide whether next step is interpolation or extrapolation.

   IF (SIGN(1.0_8,FB) .NE. SIGN(1.0_8,FC)) GO TO 3
   C = A
   FC = FA
   GO TO 3

! Finished.  Process results for proper setting of IFLAG.

10 IF (SIGN(1.0_8,FB) .EQ. SIGN(1.0_8,FC)) GO TO 13
   IF (ABS(FB) .GT. FX) GO TO 12
   IFLAG = 1
   RETURN
11 IFLAG = 2
   RETURN
12 IFLAG = 3
   RETURN
13 IFLAG = 4
   RETURN
14 IFLAG = 5
   RETURN

END SUBROUTINE DFZERO



SUBROUTINE tWeights(Nlag, QIn, n, output)
IMPLICIT NONE
INTEGER :: wc, wf, i, n
REAL(8) :: Nlag, m
REAL(8) :: QIn(n), output(n), QOut(n)
REAL(8) :: w(ceiling(Nlag)), LagArr(ceiling(Nlag))

wc = ceiling(Nlag)
wf = floor(Nlag)

m = 2._8 / Nlag**2

do i = 1, wf

!   w(i) = m * i - m / 2._8
  w(i) = (2 * m * i - m) / 2._8

enddo

if (wf < wc) then

  w(wc) = (m * Nlag**2 - m * wf**2) / 2._8

endif

w = w / sum(w)

LagArr(:) = 0
do i = 1, n

  LagArr = LagArr + QIn(i) * w
  QOut(i) = LagArr(1)

  LagArr(1:(wc-1)) =  LagArr(2:wc)
  LagArr(wc) = 0

enddo

output = QOut

END SUBROUTINE tWeights
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL(8) FUNCTION FUN_M1(sEn)
IMPLICIT NONE
REAL(8) :: sEn
REAL(8) :: sSt, P, K, m, alpha, Epot, delT
COMMON /coeff1/ sSt, P, K, m, alpha, Epot, delT
FUN_M1 = sEn - sSt + delT * (-P + K * sEn**alpha + Epot * (1 - exp(-sEn / m)))
END FUNCTION FUN_M1


SUBROUTINE MI(x1, x2, n, alpha_Qq_FR, K_Qq_FR, Ce, m_E_FR, dT, output) 
IMPLICIT NONE
INTEGER :: n, i
INTEGER :: iflag
REAL(8) :: sMin, sMax, EPS
REAL(8) :: x1(n), x2(n), output(n)
REAL(8) :: alpha_Qq_FR, K_Qq_FR, Ce, m_E_FR, dT
REAL(8) :: sSt, P, K, m, alpha, Epot, delT

REAL(8), EXTERNAL :: FUN_M1
COMMON /coeff1/ sSt, P, K, m, alpha, Epot, delT


 K=K_Qq_FR
 alpha=alpha_Qq_FR
 delT=dT
 m=m_E_FR
 EPS=1e-9_8
 sSt = 0._8
 
 do i = 1, n

  sMin =0._8
  P = x1(i)
  Epot = Ce*x2(i)
  sMax = min(sSt + P * delT, 10000._8)
  CALL DFZERO(FUN_M1, sMin, sMax, sMin, EPS, EPS, iflag)
  sSt=sMin
  output(i)=K*(sSt)**alpha
  
 enddo

END SUBROUTINE MI

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL(8) FUNCTION FUN_M2(sEn)
IMPLICIT NONE
REAL(8) :: sEn
REAL(8) :: sSt, Smax_U, P, K, Gamma, Beta, Epot, delT
COMMON /coeff2/ sSt, Smax_U, P, K, Gamma, Beta, Epot, delT
FUN_M2 = sEn-sSt-delT*(P-P*(sEn/Smax_U)**Beta-K*(sEn/Smax_U)-Epot*(sEn/Smax_U)*(1+Gamma)/((sEn/Smax_U)+Gamma))
END FUNCTION FUN_M2

SUBROUTINE MII(x1, x2, n, Ce, Beta_Qq_UR, Smax_UR, K_Qb_UR, Beta_E_UR, SiniFr_UR, dT, output) 
IMPLICIT NONE
INTEGER :: n, i
INTEGER :: iflag
REAL(8) :: sMin, sMax, EPS
REAL(8) :: x1(n), x2(n), output(n)
REAL(8) :: Ce, Beta_Qq_UR, Smax_UR, K_Qb_UR, Beta_E_UR, SiniFr_UR, dT
REAL(8) :: sSt, Smax_U, P, K, Gamma, Beta, Epot, delT, Qu, Qb

REAL(8), EXTERNAL :: FUN_M2
COMMON /coeff2/ sSt, Smax_U, P, K, Gamma, Beta, Epot, delT


 K=K_Qb_UR
 Beta=Beta_Qq_UR
 delT=dT
 Gamma=Beta_E_UR
 Smax_U=Smax_UR
 EPS=1e-9_8
 sSt = SiniFr_UR*Smax_U
 
 do i = 1, n

  P = x1(i)
  Epot = Ce*x2(i)
  sMax=min(sSt+P*delT,Smax_U)
  sMin=max(sSt-Epot*delT,0._8)
  CALL DFZERO(FUN_M2, sMin, sMax, sMin, EPS, EPS, iflag)
  sSt=sMin
  Qu=P*(sSt/Smax_U)**Beta
  Qb=K*(sSt/Smax_U)
  output(i)=Qu+Qb
  
 enddo

END SUBROUTINE MII

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL(8) FUNCTION FUN_M3_1(sEn)
IMPLICIT NONE
REAL(8) :: sEn
REAL(8) :: sSt, Smax_U, P, Gamma, Beta, Epot, delT
COMMON /coeff3_1/ sSt, Smax_U, P, Gamma, Beta, Epot, delT
FUN_M3_1 = sEn-sSt-delT*(P-P*(1-((1-(sEn/Smax_U))*(1+Beta)/(1-(sEn/Smax_U)+Beta)))-Epot*(sEn/Smax_U)*(1+Gamma)/((sEn/Smax_U)+Gamma))
END FUNCTION FUN_M3_1

REAL(8) FUNCTION FUN_M3_2(sEn_f)
IMPLICIT NONE
REAL(8) :: sEn_f
REAL(8) :: sSt_f, P_f, K_f, alpha, delT_f
COMMON /coeff3_2/ sSt_f, P_f, K_f, alpha,delT_f
FUN_M3_2 = sEn_f-sSt_f-P_f*delT_f+K_f*(sEn_f**alpha)*delT_f
END FUNCTION FUN_M3_2

SUBROUTINE MIII(x1, x2, n, alpha_Qq_FR, K_Qq_FR, Ce, Smax_UR, Beta_Qq_UR, Beta_E_UR, SiniFr_UR, dT, output) 
IMPLICIT NONE
INTEGER :: n, i
INTEGER :: iflag
REAL(8) :: sMin, sMax, sMin_f, sMax_f, EPS
REAL(8) :: x1(n), x2(n), output(n)
REAL(8) :: alpha_Qq_FR, K_Qq_FR, Ce, Smax_UR, Beta_Qq_UR, Beta_E_UR, SiniFr_UR, dT
REAL(8) :: sSt, Smax_U, P, Gamma, Beta, Epot, delT
REAL(8) :: sSt_f, P_f, K_f, alpha, delT_f

REAL(8), EXTERNAL :: FUN_M3_1, FUN_M3_2
COMMON /coeff3_1/ sSt, Smax_U, P, Gamma, Beta, Epot, delT
COMMON /coeff3_2/ sSt_f, P_f, K_f, alpha, delT_f

 K_f=K_Qq_FR
 alpha=alpha_Qq_FR
 delT_f=dT
 
 Beta=Beta_Qq_UR
 delT=dT
 Gamma=Beta_E_UR
 Smax_U=Smax_UR
 EPS=1e-9_8
 sSt = SiniFr_UR*Smax_U
 sSt_f= 0._8
 
 do i = 1, n

  P = x1(i)
  Epot = Ce*x2(i)
  sMax=min(sSt+P*delT,Smax_U)
  sMin=max(sSt-Epot*delT,0._8)
  CALL DFZERO(FUN_M3_1, sMin, sMax, sMin, EPS, EPS, iflag)
  sSt=sMin
  P_f=P*(1-((1-(sSt/Smax_U))*(1+Beta)/(1-(sSt/Smax_U)+Beta)))
  sMax_f=min(sSt_f+P_f*delT_f,10000._8)
  sMin_f=0._8
  CALL DFZERO(FUN_M3_2, sMin_f, sMax_f, sMin_f, EPS, EPS, iflag )
  sSt_f=sMin_f
  output(i)=K_f*(sSt_f)**alpha
  
 enddo

END SUBROUTINE MIII

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL(8) FUNCTION FUN_M4(sEn)
IMPLICIT NONE
REAL(8) :: sEn
REAL(8) :: sSt, Smax_U, P, Gamma, Beta, Epot, delT
COMMON /coeff4/ sSt, Smax_U, P, Gamma, Beta, Epot, delT
FUN_M4 = sEn-sSt-delT*(P-P*(sEn/Smax_U)**Beta-Epot*(sEn/Smax_U)*(1+Gamma)/((sEn/Smax_U)+Gamma))
END FUNCTION FUN_M4

SUBROUTINE MIV(x1, x2, n, alpha_Qq_FR, K_Qq_FR, Ce, Smax_UR, Beta_Qq_UR, Beta_E_UR, SiniFr_UR, dT, output) 
IMPLICIT NONE
INTEGER :: n, i
INTEGER :: iflag
REAL(8) :: sMin, sMax, sMin_f, sMax_f, EPS
REAL(8) :: x1(n), x2(n), output(n)
REAL(8) :: alpha_Qq_FR, K_Qq_FR, Ce, Smax_UR, Beta_Qq_UR, Beta_E_UR, SiniFr_UR, dT
REAL(8) :: sSt, Smax_U, P, Gamma, Beta, Epot, delT
REAL(8) :: sSt_f, P_f, K_f, alpha, delT_f

REAL(8), EXTERNAL :: FUN_M3_2, FUN_M4
COMMON /coeff3_2/ sSt_f, P_f, K_f, alpha, delT_f
COMMON /coeff4/ sSt, Smax_U, P, Gamma, Beta, Epot, delT

 K_f=K_Qq_FR
 alpha=alpha_Qq_FR
 delT_f=dT
 
 Beta=Beta_Qq_UR
 delT=dT
 Gamma=Beta_E_UR
 Smax_U=Smax_UR
 EPS=1e-9_8
 sSt = SiniFr_UR*Smax_U
 sSt_f=0._8
 
 do i = 1, n

  P = x1(i)
  Epot = Ce*x2(i)
  sMax=min(sSt+P*delT,Smax_U)
  sMin=max(sSt-Epot*delT,0._8)
  CALL DFZERO(FUN_M4, sMin, sMax, sMin, EPS, EPS, iflag)
  sSt=sMin
  P_f=P*(sSt/Smax_U)**Beta
  sMax_f=min(sSt_f+P_f*delT_f,10000._8)
  sMin_f=0._8
  CALL DFZERO(FUN_M3_2, sMin_f, sMax_f, sMin_f, EPS, EPS, iflag )
  sSt_f=sMin_f
  output(i)=K_f*(sSt_f)**alpha
 enddo

END SUBROUTINE MIV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE MV(x1, x2, n, alpha_Qq_FR, K_Qq_FR, Ce, Smax_UR, Beta_Qq_UR, Beta_E_UR, SiniFr_UR, dT, Tlag, output) 
IMPLICIT NONE
INTEGER :: n, i
INTEGER :: iflag
REAL(8) :: sMin, sMax, sMin_f, sMax_f, EPS
REAL(8) :: x1(n), x2(n), output(n)
REAL(8) :: alpha_Qq_FR, K_Qq_FR, Ce, Smax_UR, Beta_Qq_UR, Beta_E_UR, SiniFr_UR, dT, Tlag
REAL(8) :: sSt, Smax_U, P, Gamma, Beta, Epot, delT
REAL(8) :: sSt_f, P_f, K_f, alpha, delT_f
REAL(8) :: P_f_vec(n), P_fl(n)

REAL(8), EXTERNAL :: FUN_M3_2, FUN_M4
COMMON /coeff3_2/ sSt_f, P_f, K_f, alpha, delT_f
COMMON /coeff4/ sSt, Smax_U, P, Gamma, Beta, Epot, delT


 K_f=K_Qq_FR
 alpha=alpha_Qq_FR
 delT_f=dT
 
 Beta=Beta_Qq_UR
 delT=dT
 Gamma=Beta_E_UR
 Smax_U=Smax_UR
 EPS=1e-9_8
 sSt = SiniFr_UR*Smax_U
 sSt_f=0._8
 
 do i = 1, n

  P = x1(i)
  Epot = Ce*x2(i)
  sMax=min(sSt+P*delT,Smax_U)
  sMin=max(sSt-Epot*delT,0._8)
  CALL DFZERO(FUN_M4, sMin, sMax, sMin, EPS, EPS, iflag)
  sSt=sMin
  P_f=P*(sSt/Smax_U)**Beta
  P_f_vec(i)=P_f
  
 end do

  CALL tWeights(Tlag, P_f_vec, n, P_fl)
  
 do i = 1, n

  P_f=P_fl(i)
  sMax_f=min(sSt_f+P_f*delT_f,10000._8)
  sMin_f=0._8
  CALL DFZERO(FUN_M3_2, sMin_f, sMax_f, sMin_f, EPS, EPS, iflag )
  sSt_f=sMin_f
  output(i)=K_f*(sSt_f)**alpha
  
 enddo

END SUBROUTINE MV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL(8) FUNCTION FUN_M6_1(sEn_I)
IMPLICIT NONE
REAL(8) :: sEn_I
REAL(8) :: sSt_I, Smax_I, P_I, Gamma_I, Epot_I, delT_I
COMMON /coeff6_1/ sSt_I, Smax_I, P_I, Gamma_I, Epot_I, delT_I
FUN_M6_1 = sEn_I-sSt_I-delT_I*(P_I-P_I*(1-((1-(sEn_I/Smax_I))*(1+Gamma_I)/(1-(sEn_I/Smax_I)+Gamma_I))) &
& -Epot_I*(sEn_I/Smax_I)*(1+Gamma_I)/((sEn_I/Smax_I)+Gamma_I))
END FUNCTION FUN_M6_1

REAL(8) FUNCTION FUN_M6_2(sEn)
IMPLICIT NONE
REAL(8) :: sEn
REAL(8) :: sSt, Smax_U, P_U, Beta, Gamma, Epot_U, delT
COMMON /coeff6_2/ sSt, Smax_U, P_U, Beta, Gamma, Epot_U, delT
FUN_M6_2 = sEn-sSt-delT*(P_U-P_U*(sEn/Smax_U)**Beta-(Epot_U)*(sEn/Smax_U)*(1+Gamma)/((sEn/Smax_U)+Gamma))
END FUNCTION FUN_M6_2

SUBROUTINE MVI(x1, x2, n, alpha_Qq_FR, Beta_Qq_UR, K_Qq_FR, Ce, Smax_UR, Smax_IR, m_QE_IR, Beta_E_UR, SiniFr_UR, dT, Tlag, output) 
IMPLICIT NONE
INTEGER :: n, i
INTEGER :: iflag
REAL(8) :: sMin, sMax, sMin_f, sMax_f, Min_ir, Max_ir, EPS
REAL(8) :: x1(n), x2(n), output(n)
REAL(8) :: alpha_Qq_FR, Beta_Qq_UR, K_Qq_FR, Ce, Smax_UR, Smax_IR, m_QE_IR, Beta_E_UR, SiniFr_UR, dT, Tlag
REAL(8) :: sSt_I, Smax_I, P_I, Gamma_I, Epot_I, delT_I, E_I
REAL(8) :: sSt, Smax_U, P_U, Beta, Gamma, delT, Epot_U
REAL(8) :: sSt_f, P_f, K_f, alpha, delT_f
REAL(8) :: P_f_vec(n), P_fl(n)

REAL(8), EXTERNAL :: FUN_M3_2, FUN_M6_1, FUN_M6_2
COMMON /coeff3_2/ sSt_f, P_f, K_f, alpha, delT_f
COMMON /coeff6_1/ sSt_I, Smax_I, P_I, Gamma_I, Epot_I, delT_I
COMMON /coeff6_2/ sSt, Smax_U, P_U, Beta, Gamma, Epot_U, delT

 delT_I=dT
 Gamma_I=m_QE_IR
 Smax_I=Smax_IR
 
 K_f=K_Qq_FR
 alpha=alpha_Qq_FR
 delT_f=dT
 
 Beta=Beta_Qq_UR
 delT=dT
 Gamma=Beta_E_UR
 Smax_U=Smax_UR
 EPS=1e-9_8
 
 sSt = SiniFr_UR*Smax_U
 sSt_f=0._8
 sSt_I=0._8
 
 do i = 1, n

  P_I = x1(i)
  Epot_I = Ce*x2(i)
  Max_ir=min(sSt_I+P_I*delT_I,Smax_I)
  Min_ir=0._8
  CALL DFZERO(FUN_M6_1, Min_ir, Max_ir, Min_ir, EPS, EPS, iflag)
  sSt_I=Min_ir
  E_I=Epot_I*(sSt_I/Smax_I)*(1+Gamma_I)/((sSt_I/Smax_I)+Gamma_I)
  P_U=P_I*(1-((1-(sSt_I/Smax_I))*(1+Gamma_I)/(1-(sSt_I/Smax_I)+Gamma_I)))
  sMax=min(sSt+P_U*delT,Smax_U)
  sMin=max(sSt-Epot_I*delT,0._8)
  Epot_U=(Epot_I-E_I)
  CALL DFZERO(FUN_M6_2, sMin, sMax, sMin, EPS, EPS, iflag)
  sSt=sMin
  P_f=P_U*(sSt/Smax_U)**Beta
  P_f_vec(i)=P_f
  
 end do

  CALL tWeights(Tlag, P_f_vec, n, P_fl)
  

 do i = 1, n
  P_f=P_fl(i)
  sMax_f=min(sSt_f+P_f*delT_f,10000._8)
  sMin_f=0._8
  CALL DFZERO(FUN_M3_2, sMin_f, sMax_f, sMin_f, EPS, EPS, iflag )
  sSt_f=sMin_f
  output(i)=K_f*(sSt_f)**alpha
 enddo

END SUBROUTINE MVI

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINE MVII(x1, x2, n, alpha_Qq_FR, K_Qq_FR, Ce, Smax_UR, Beta_Qq_UR, Beta_E_UR, SiniFr_UR, D_R, K_Qq_RR,&
! &dT, Tlag, output, ctrlVar1, ctrlVar2) 
SUBROUTINE MVII(x1, x2, n, alpha_Qq_FR, K_Qq_FR, Ce, Smax_UR, Beta_Qq_UR, Beta_E_UR, SiniFr_UR, D_R, K_Qq_RR,&
&dT, Tlag, output) 
IMPLICIT NONE
INTEGER :: n, i
INTEGER :: iflag
REAL(8) :: sMin, sMax, sMin_f, sMax_f, EPS
REAL(8) :: x1(n), x2(n), output(n)
REAL(8) :: alpha_Qq_FR, K_Qq_FR, Ce, Smax_UR, Beta_Qq_UR, Beta_E_UR, SiniFr_UR, D_R, K_Qq_RR, dT, Tlag
REAL(8) :: sSt, Smax_U, P, Gamma, Beta, Epot, delT
REAL(8) :: sSt_f, P_f, K_f, alpha, delT_f, Qf(n)
REAL(8) :: sSt_r, P_r, K_r, delT_r, Qr(n)
REAL(8) :: P_f_vec(n), P_fl(n)

! REAL(8) :: ctrlVar1(n), ctrlVar2(n)

REAL(8), EXTERNAL :: FUN_M3_2, FUN_M4
COMMON /coeff3_2/ sSt_f, P_f, K_f, alpha, delT_f
COMMON /coeff4/ sSt, Smax_U, P, Gamma, Beta, Epot, delT

 K_f=K_Qq_FR
 alpha=alpha_Qq_FR
 delT_f=dT
 
 Beta=Beta_Qq_UR
 delT=dT
 Gamma=Beta_E_UR
 Smax_U=Smax_UR
 EPS=1e-9_8
 
 K_r=K_Qq_RR
 delT_r=dT
 
 sSt = SiniFr_UR*Smax_U
 sSt_f=0._8
 sSt_r=0._8
 
 do i = 1, n
  
  P_r = D_R*x1(i)
  sSt_r=(P_r*delT_r+sSt_r)/(1+K_r*delT_r)
  sSt_r=min(sSt_r,10000._8)
  sSt_r=max(sSt_r,0._8)
  Qr(i)=K_r* sSt_r
  P=x1(i)-P_r
  Epot = Ce*x2(i)
  sMax=min(sSt+P*delT,Smax_U)
  sMin=max(sSt-Epot*delT,0._8)
  CALL DFZERO(FUN_M4, sMin, sMax, sMin, EPS, EPS, iflag)
  sSt=sMin
  P_f=P*(sSt/Smax_U)**Beta
  P_f_vec(i)=P_f  
  
end do

CALL tWeights(Tlag, P_f_vec, n, P_fl)
  
do i = 1, n

  P_f=P_fl(i)     
!   ctrlVar2(i) = P_f
  sMax_f=min(sSt_f+P_f*delT_f,10000._8) 
  sMin_f=0._8
  CALL DFZERO(FUN_M3_2, sMin_f, sMax_f, sMin_f, EPS, EPS, iflag )
  sSt_f=sMin_f
  Qf(i)=K_f*(sSt_f)**alpha
!   ctrlVar1(i) = Qf(i)
    
enddo

output=Qf+Qr

END SUBROUTINE MVII

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL(8) FUNCTION FUN_M8(sEn_f)
IMPLICIT NONE
REAL(8) :: sEn_f
REAL(8) :: sSt_f, P_f, m_f, K_f, Epot_f, delT_f
COMMON /coeff8/ sSt_f, P_f, m_f, K_f, Epot_f, delT_f
FUN_M8 = sEn_f-sSt_f-delT_f*(P_f-(K_f*sEn_f)-(Epot_f)*(1-exp(-sEn_f/m_f)))
END FUNCTION FUN_M8

SUBROUTINE MVIII(x1, x2, n,  K_Qq_FR, Ce, K_Qq_SR, D_S, m_E_FR, dT, output) 
IMPLICIT NONE
INTEGER :: n, i
INTEGER :: iflag
REAL(8) :: sMin_f, sMax_f, EPS
REAL(8) :: x1(n), x2(n), output(n)
REAL(8) :: K_Qq_FR, Ce, K_Qq_SR, D_S, m_E_FR, dT
REAL(8) :: sSt_f, P_f, m_f, K_f, Epot_f, delT_f, Qf(n)
REAL(8) :: sSt_s, P_s, K_s, delT_s, Qs(n)

REAL(8),EXTERNAL :: FUN_M8
COMMON /coeff8/ sSt_f, P_f, m_f, K_f, Epot_f, delT_f

 K_f=K_Qq_FR
 m_f=m_E_FR
 delT_f=dT
 sSt_f=0._8 
 EPS=1e-9_8
 
 K_s=K_Qq_SR
 delT_s=dT
 sSt_s=0._8

do i = 1, n

  P_s = D_S*x1(i)
  P_f = x1(i)-P_s
  Epot_f = Ce*x2(i)
  sSt_s=(P_s*delT_s+sSt_s)/(1+K_s*delT_s)
  sSt_s=min(sSt_s,10000._8)
  sSt_s=max(sSt_s,0._8)
  Qs(i)=K_s* sSt_s
  
  sMax_f=min(sSt_f+P_f*delT_f,10000._8)
  sMin_f=0._8
  CALL DFZERO(FUN_M8, sMin_f, sMax_f, sMin_f, EPS, EPS, iflag)
  sSt_f=sMin_f
  Qf(i)=K_f*sSt_f
  
end do

output=Qf+Qs

END SUBROUTINE MVIII

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MIX(x1, x2, n, Beta_Qq_UR, K_Qq_FR, Ce, K_Qq_SR, D_S, Smax_UR, Beta_E_UR, SiniFr_UR, dT, output) 
IMPLICIT NONE
INTEGER :: n, i
INTEGER :: iflag
REAL(8) :: sMin, sMax, EPS
REAL(8) :: x1(n), x2(n), output(n)
REAL(8) :: Beta_Qq_UR, K_Qq_FR, Ce, K_Qq_SR, D_S, Smax_UR, Beta_E_UR, SiniFr_UR, dT
REAL(8) :: sSt, Smax_U, P, Gamma, Beta, Epot, delT, Qu
REAL(8) :: sSt_s, P_s, K_s, delT_s, Qs(n)
REAL(8) :: sSt_f, P_f, K_f, delT_f, Qf(n)

REAL(8), EXTERNAL :: FUN_M4
COMMON /coeff4/ sSt, Smax_U, P, Gamma, Beta, Epot, delT

 K_f=K_Qq_FR
 delT_f=dT
 sSt_f=0._8
 
 Beta=Beta_Qq_UR
 delT=dT
 Gamma=Beta_E_UR
 Smax_U=Smax_UR
 EPS=1e-9_8
 sSt = SiniFr_UR*Smax_U
 
 K_s=K_Qq_SR
 delT_s=dT
 sSt_s=0._8

do i = 1, n
  
  P=x1(i)
  Epot = Ce*x2(i)
  sMax=min(sSt+P*delT,Smax_U)
  sMin=max(sSt-Epot*delT,0._8)
  CALL DFZERO(FUN_M4, sMin, sMax, sMin, EPS, EPS, iflag)
  sSt=sMin
  Qu=P*(sSt/Smax_U)**Beta
  
  P_s = D_S*Qu
  P_f = Qu-P_s

  sSt_s=(P_s*delT_s+sSt_s)/(1+K_s*delT_s)
  sSt_s=min(sSt_s,10000._8)
  sSt_s=max(sSt_s,0._8)
  Qs(i)=K_s* sSt_s
  
  sSt_f=(P_f*delT_f+sSt_f)/(1+K_f*delT_f)
  sSt_f=min(sSt_f,10000._8)
  sSt_f=max(sSt_f,0._8)
  Qf(i)=K_f* sSt_f

end do

output=Qf+Qs

END SUBROUTINE MIX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL(8) FUNCTION FUN_M5(sEn)
IMPLICIT NONE
REAL(8) :: sEn
REAL(8) :: sSt, Smax_U, P, Gamma, Epot, delT
COMMON /coeff5/ sSt, Smax_U, P, Gamma, Epot, delT
FUN_M5 = sEn-sSt-delT*(P-P*(sEn/Smax_U)-Epot*(sEn/Smax_U)*(1+Gamma)/((sEn/Smax_U)+Gamma))
END FUNCTION FUN_M5

SUBROUTINE MX(x1, x2, n, K_Qq_FR, Ce, K_Qq_SR, D_S, Smax_UR, Beta_E_UR, SiniFr_UR, dT, Tlag, output) 
IMPLICIT NONE
INTEGER :: n, i
INTEGER :: iflag
REAL(8) :: sMin, sMax, EPS
REAL(8) :: x1(n), x2(n), output(n)
REAL(8) :: K_Qq_FR, Ce, K_Qq_SR, D_S, Smax_UR, Beta_E_UR, SiniFr_UR, dT, Tlag
REAL(8) :: sSt, Smax_U, P, Gamma, Epot, delT, Qu
REAL(8) :: sSt_s, P_s, K_s, delT_s, Qs(n)
REAL(8) :: sSt_f, P_f, K_f, delT_f, Qf(n)
REAL(8) :: P_f_vec(n), P_fl(n)

REAL(8), EXTERNAL :: FUN_M5
COMMON /coeff5/ sSt, Smax_U, P, Gamma, Epot, delT

 K_f=K_Qq_FR
 delT_f=dT
 sSt_f=0._8

 delT=dT
 Gamma=Beta_E_UR
 Smax_U=Smax_UR
 EPS=1e-9_8
 sSt = SiniFr_UR*Smax_U
 
 K_s=K_Qq_SR
 delT_s=dT
 sSt_s=0._8

do i = 1, n
  
  P=x1(i)
  Epot = Ce*x2(i)
  sMax=min(sSt+P*delT,Smax_U)
  sMin=max(sSt-Epot*delT,0._8)
  CALL DFZERO(FUN_M5, sMin, sMax, sMin, EPS, EPS, iflag)
  sSt=sMin
  Qu=P*(sSt/Smax_U)
  
  P_s = D_S*Qu
  P_f = Qu-P_s
  P_f_vec(i)= P_f

  sSt_s=(P_s*delT_s+sSt_s)/(1+K_s*delT_s)
  sSt_s=min(sSt_s,10000._8)
  sSt_s=max(sSt_s,0._8)
  Qs(i)=K_s* sSt_s

end do

  CALL tWeights(Tlag, P_f_vec, n, P_fl)
  
do i = 1, n  

  sSt_f=(P_fl(i)*delT_f+sSt_f)/(1+K_f*delT_f)
  sSt_f=min(sSt_f,10000._8)
  sSt_f=max(sSt_f,0._8)
  Qf(i)=K_f* sSt_f

end do

output=Qf+Qs

END SUBROUTINE MX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE MXI(x1, x2, n, Beta_Qq_UR, K_Qq_FR, Ce, K_Qq_SR, D_S, Smax_UR, Beta_E_UR, SiniFr_UR, dT, Tlag, output) 
IMPLICIT NONE
INTEGER :: n, i
INTEGER :: iflag
REAL(8) :: sMin, sMax, EPS
REAL(8) :: x1(n), x2(n), output(n)
REAL(8) :: Beta_Qq_UR, K_Qq_FR, Ce, K_Qq_SR, D_S, Smax_UR, Beta_E_UR, SiniFr_UR, dT, Tlag
REAL(8) :: sSt, Smax_U, P, Gamma, Beta, Epot, delT, Qu
REAL(8) :: sSt_s, P_s, K_s, delT_s, Qs(n)
REAL(8) :: sSt_f, P_f, K_f, delT_f, Qf(n)
REAL(8) :: P_f_vec(n), P_fl(n)

REAL(8), EXTERNAL :: FUN_M4
COMMON /coeff4/ sSt, Smax_U, P, Gamma, Beta, Epot, delT

 K_f=K_Qq_FR
 delT_f=dT
 sSt_f=0._8
 
 Beta=Beta_Qq_UR
 delT=dT
 Gamma=Beta_E_UR
 Smax_U=Smax_UR
 EPS=1e-9_8
 sSt = SiniFr_UR*Smax_U
 
 K_s=K_Qq_SR
 delT_s=dT
 sSt_s=0._8

do i = 1, n
  
  P=x1(i)
  Epot = Ce*x2(i)
  sMax=min(sSt+P*delT,Smax_U)
  sMin=max(sSt-Epot*delT,0._8)
  CALL DFZERO(FUN_M4, sMin, sMax, sMin, EPS, EPS, iflag)
  sSt=sMin
  Qu=P*(sSt/Smax_U)**Beta
  
  P_s = D_S*Qu
  P_f = Qu-P_s
  P_f_vec(i)= P_f

  sSt_s=(P_s*delT_s+sSt_s)/(1+K_s*delT_s)
  sSt_s=min(sSt_s,10000._8)
  sSt_s=max(sSt_s,0._8)
  Qs(i)=K_s* sSt_s

end do

  CALL tWeights(Tlag, P_f_vec, n, P_fl)
  
do i = 1, n  

  sSt_f=(P_fl(i)*delT_f+sSt_f)/(1+K_f*delT_f)
  sSt_f=min(sSt_f,10000._8)
  sSt_f=max(sSt_f,0._8)
  Qf(i)=K_f* sSt_f

end do

output=Qf+Qs

END SUBROUTINE MXI

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MXII(x1, x2, n, Beta_Qq_UR, K_Qq_FR, Ce, K_Qq_SR, D_S, SiniFr_UR, Smax_UR, & 
& Smax_IR, m_QE_IR, Beta_E_UR, dT, Tlag, output) 
IMPLICIT NONE
INTEGER :: n, i
INTEGER :: iflag
REAL(8) :: sMin, sMax, Min_ir, Max_ir, EPS
REAL(8) :: x1(n), x2(n), output(n)
REAL(8) :: Beta_Qq_UR, K_Qq_FR, Ce, K_Qq_SR, D_S, SiniFr_UR, Smax_UR, Smax_IR, m_QE_IR, Beta_E_UR, dT, Tlag
REAL(8) :: sSt_I, Smax_I, P_I, Gamma_I, Epot_I, delT_I, E_I
REAL(8) :: sSt, Smax_U, P_U, Beta, Gamma, delT, Epot_U, Qu
REAL(8) :: sSt_f, P_f, K_f, delT_f, Qf(n)
REAL(8) :: sSt_s, P_s, K_s, delT_s, Qs(n)
REAL(8) :: P_f_vec(n), P_fl(n)

REAL(8), EXTERNAL :: FUN_M6_1, FUN_M6_2
COMMON /coeff6_1/ sSt_I, Smax_I, P_I, Gamma_I, Epot_I, delT_I
COMMON /coeff6_2/ sSt, Smax_U, P_U, Beta, Gamma, Epot_U, delT

 delT_I=dT
 Gamma_I=m_QE_IR
 Smax_I=Smax_IR
 
 K_f=K_Qq_FR
 delT_f=dT
 
 Beta=Beta_Qq_UR
 delT=dT
 Gamma=Beta_E_UR
 Smax_U=Smax_UR
 EPS=1e-9_8
 
 K_s=K_Qq_SR
 delT_s=dT
 sSt_s=0._8
 
 sSt = SiniFr_UR*Smax_U
 sSt_f=0._8
 sSt_I=0._8
 sSt_s=0._8
 
 do i = 1, n

  P_I = x1(i)
  Epot_I = Ce*x2(i)
  Max_ir=min(sSt_I+P_I*delT_I,Smax_I)
  Min_ir=0._8
  CALL DFZERO(FUN_M6_1, Min_ir, Max_ir, Min_ir, EPS, EPS, iflag)
  sSt_I=Min_ir
  E_I=Epot_I*(sSt_I/Smax_I)*(1+Gamma_I)/((sSt_I/Smax_I)+Gamma_I)
  P_U=P_I*(1-((1-(sSt_I/Smax_I))*(1+Gamma_I)/(1-(sSt_I/Smax_I)+Gamma_I)))
  sMax=min(sSt+P_U*delT,Smax_U)
  sMin=max(sSt-Epot_I*delT,0._8)
  Epot_U=(Epot_I-E_I)
  CALL DFZERO(FUN_M6_2, sMin, sMax, sMin, EPS, EPS, iflag)
  sSt=sMin
  Qu=P_U*(sSt/Smax_U)**Beta
  P_s=D_S*Qu
  sSt_s=(P_s*delT_s+sSt_s)/(1+K_s*delT_s)
  sSt_s=min(sSt_s,10000._8)
  sSt_s=max(sSt_s,0._8)
  Qs(i)=K_s* sSt_s
  P_f=Qu-P_s
  P_f_vec(i)=P_f
 
 end do

  CALL tWeights(Tlag, P_f_vec, n, P_fl)
  
 do i = 1, n
 
  sSt_f=(P_fl(i)*delT_f+sSt_f)/(1+K_f*delT_f)
  sSt_f=min(sSt_f,10000._8)
  sSt_f=max(sSt_f,0._8)
  Qf(i)=K_f* sSt_f
  
 enddo
 
 output=Qs+Qf

END SUBROUTINE MXII
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Constituitive functions

! choice=1: power
! power function
SUBROUTINE power(x,m,y)
IMPLICIT NONE
REAL*8 :: x,m,y
y=x**m
END SUBROUTINE power

! choice=2: rhfun, only in case of interception reservoir, in case of UCR choice=2 is mlc, otherwise choice=2 is linear
! reflected hyperbolic function
SUBROUTINE rhfun(x,m,y)
IMPLICIT NONE
REAL*8 :: x,m,y
y=1-((1-x)*(1+m)/(1-x+m))
END SUBROUTINE rhfun

! choice=3: mkinetic
! monod type kinetics
SUBROUTINE mkinetic(x,m,y)
IMPLICIT NONE
REAL*8 :: x,m,y
y=(x*(1+m))/(x+m)
END SUBROUTINE mkinetic

!Tessier function
SUBROUTINE tfun(x,m,y)
IMPLICIT NONE
REAL*8 :: x,m,y
y=1-exp(-x/m)
END SUBROUTINE tfun

!choice=4: mlc
!modified logistic curve
SUBROUTINE mlc(x,m,lamda,y)
IMPLICIT NONE
REAL*8 :: x,m,lamda,y
y=(1+exp(-m*(1-lamda)))*(exp(-m*x)-1)/(1+exp(-m*(x-lamda)))*(exp(-m)-1)
END SUBROUTINE mlc

!functions
REAL*8 FUNCTION FUN_Mw(sEn_w)
REAL*8 :: sEn_w
REAL*8 :: sSt_w, P_w, m_w, T, Tm_w, Kq_w, delT_w, out_w
COMMON /coeffi/ sSt_w, P_w, m_w, T, Tm_w, Kq_w, delT_w, out_w
call tfun(sEn_w,m_w,out_w)
if (T <= Tm_w) then
FUN_Mw = sEn_w-sSt_w-delT_w*P_w
else if (T > Tm_w) then
FUN_Mw = sEn_w-sSt_w-delT_w*(P_w-Kq_w*(T-Tm_w)*out_w)
end if
END FUNCTION FUN_Mw


REAL*8 FUNCTION FUN_Mi(sEn_i)
REAL*8 :: sEn_i
INTEGER:: choice_i
REAL*8 :: sSt_i, Smax_i, P_i, m_i, Epot_i, delT_i, out_i, out1_i, out2_i
COMMON /coeffi/ sSt_i, Smax_i, P_i, m_i, Epot_i, delT_i, choice_i
call mkinetic((sEn_i/Smax_i),m_i,out_i)
call power((sEn_i/Smax_i),m_i,out1_i)
call rhfun((sEn_i/Smax_i),m_i,out2_i)
if (choice_i==1) then
FUN_Mi = sEn_i-sSt_i-delT_i*(P_i-P_i*out1_i-Epot_i*out_i)
else if (choice_i==2) then
FUN_Mi = sEn_i-sSt_i-delT_i*(P_i-P_i*out2_i-Epot_i*out_i)
end if
END FUNCTION FUN_Mi

REAL*8 FUNCTION FUN_Mu(sEn_u)
REAL*8 :: sEn_u
INTEGER:: choice_u
REAL*8 :: sSt_u, Smax_u, P_u, Beta, mu, Gamma, Epot_u, Kb_u, delT_u, out_u, out1_u, out2_u, out3_u
COMMON /coeffu/ sSt_u, Smax_u, P_u, Beta, mu, Gamma, Epot_u,  Kb_u, delT_u, choice_u
call mkinetic((sEn_u/Smax_u),Gamma,out_u)
call power((sEn_u/Smax_u),Beta,out1_u)
call mkinetic((sEn_u/Smax_u),Beta,out2_u)
call mlc((sEn_u/Smax_u),Beta,mu,out3_u)
if (choice_u==1) then
FUN_Mu = sEn_u-sSt_u-delT_u*(P_u-P_u*out1_u-Kb_u*sEn_u-(Epot_u)*out_u)
else if (choice_u==2) then
FUN_Mu = sEn_u-sSt_u-delT_u*(P_u-P_u*(sEn_u/Smax_u)-Kb_u*sEn_u-(Epot_u)*out_u)
else if (choice_u==3) then
FUN_Mu = sEn_u-sSt_u-delT_u*(P_u-P_u*out2_u-Kb_u*sEn_u-(Epot_u)*out_u)
else if (choice_u==4) then
FUN_Mu = sEn_u-sSt_u-delT_u*(P_u-P_u*out3_u-Kb_u*sEn_u-(Epot_u)*out_u)
end if
END FUNCTION FUN_Mu

REAL*8 FUNCTION FUN_Mf(sEn_f)
REAL*8 :: sEn_f
INTEGER:: choice_f
INTEGER :: res_u, res_s, res_cr
REAL*8 :: sSt_f, P_f, Kq_f, Kb_f, m_f, alpha_f, Epot_f,  D_FD, delT_f, out_f, out1_f 
COMMON /coefff/ sSt_f,P_f,Kq_f,Kb_f,m_f,alpha_f,Epot_f,D_FD,delT_f,choice_f,res_u,res_s,res_cr

 call tfun(sEn_f,m_f,out_f)
 call power(sEn_f,alpha_f,out1_f)
if (res_u==0 .AND. res_cr==0 .AND. res_s==0) then !UR, CR, SR are off
   if (choice_f==1) then
   FUN_Mf = sEn_f - sSt_f - delT_f * (P_f - Kq_f * out1_f -  (Epot_f) * out_f)
   else if (choice_f==2) then
   FUN_Mf = sEn_f - sSt_f - delT_f * (P_f - Kq_f * sEn_f -  (Epot_f) * out_f)
   end if
else if (res_u==0 .AND. res_cr==1 .AND. res_s==0) then !UR, SR are off, CR is on
   if (choice_f==1) then
   FUN_Mf = sEn_f - sSt_f - delT_f * (P_f - Kq_f * out1_f)
   else if (choice_f==2) then
   FUN_Mf = sEn_f - sSt_f - delT_f * (P_f - Kq_f * sEn_f)
   end if
else if (res_u==1 .AND. res_cr==0 .AND. res_s==1) then !UR, SR are on, CR is off
    if (choice_f==1) then
    FUN_Mf = sEn_f -sSt_f -delT_f *(P_f -Kq_f* out1_f -Kb_f* sEn_f)
    else if (choice_f==2) then
    FUN_Mf = sEn_f -sSt_f -delT_f *(P_f -Kq_f* sEn_f -Kb_f* sEn_f)
    end if
else if (res_u==1 .AND. res_cr==1 .AND. res_s==1) then !UR, SR, CR are on
    if (choice_f==1) then
    FUN_Mf = sEn_f -sSt_f -delT_f *(P_f -Kq_f* out1_f -Kb_f* sEn_f)
    else if (choice_f==2) then
    FUN_Mf = sEn_f -sSt_f -delT_f *(P_f -Kq_f* sEn_f -Kb_f* sEn_f)
    end if    
else if (res_u==1 .AND. res_cr==0 .AND. res_s==0) then !UR is on, CR and SR are off
    if (choice_f==1) then
    FUN_Mf = sEn_f-sSt_f-delT_f*(P_f-Kq_f*out1_f)
    else if (choice_f==2) then
    FUN_Mf = sEn_f-sSt_f-delT_f*(P_f-Kq_f*sEn_f)
    end if
else if (res_u==1 .AND. res_cr==1 .AND. res_s==0) then !UR, CR are on,SR is off
    if (choice_f==1) then
    FUN_Mf = sEn_f-sSt_f-delT_f*(P_f-Kq_f*out1_f)
    else if (choice_f==2) then
    FUN_Mf = sEn_f-sSt_f-delT_f*(P_f-Kq_f*sEn_f)
    end if
else if (res_u==0 .AND. res_cr==0 .AND. res_s==1) then !UR, CR are off, SR is on
    if (choice_f==1) then
     FUN_Mf = sEn_f-sSt_f-delT_f*(P_f-Kq_f*out1_f-Kb_f*sEn_f-(1-D_FD)*(Epot_f)*out_f)
     else if (choice_f==2) then
     FUN_Mf = sEn_f-sSt_f-delT_f*(P_f-Kq_f*sEn_f-Kb_f*sEn_f-(1-D_FD)*(Epot_f)*out_f)
     end if
else if (res_u==0 .AND. res_cr==1 .AND. res_s==1) then !UR is off, CR, SR are on
    if (choice_f==1) then
     FUN_Mf = sEn_f-sSt_f-delT_f*(P_f-Kq_f*out1_f-Kb_f*sEn_f)
     else if (choice_f==2) then
     FUN_Mf = sEn_f-sSt_f-delT_f*(P_f-Kq_f*sEn_f-Kb_f*sEn_f)
     end if
end if
END FUNCTION FUN_Mf

REAL*8 FUNCTION FUN_Mucr(sEn_ucr)
REAL*8 :: sEn_ucr
INTEGER:: choice_ucr
REAL*8 :: sSt_cr,sSt_ucr,sEn_scr,P_ucr,S_max_cr,S_evmax_cr,U_max_ucr,Epot_ucr,Beta_cr,Beta_ucr,mu_ucr
REAL*8 :: Kb_ucr,delT_ucr,out_ucr,out1_ucr,out2_ucr,Sm1
COMMON /coeffs/ sSt_cr,sSt_ucr,sEn_scr,P_ucr,S_max_cr,S_evmax_cr,U_max_ucr,Epot_ucr,Beta_cr,Beta_ucr,&
mu_ucr,Kb_ucr,delT_ucr,choice_ucr

 call power(((sEn_ucr/(S_max_cr-sEn_scr))/U_max_ucr),Beta_ucr,out1_ucr)
 call mlc(((sEn_ucr/(S_max_cr-sEn_scr))/U_max_ucr),Beta_ucr,mu_ucr,out2_ucr)
 call mkinetic(((sEn_ucr+(sEn_scr-(S_max_cr-S_evmax_cr)))/S_evmax_cr),Beta_cr,out_ucr)

Sm_1 =S_max_cr-S_evmax_cr
if (choice_ucr==1) then
FUN_Mucr = sEn_ucr-sSt_ucr-delT_ucr*(P_ucr-P_ucr*out1_ucr-Kb_ucr*sEn_ucr-Epot_ucr*out_ucr*(sEn_ucr/(sEn_ucr+(sEn_scr-Sm1))))-&
(sEn_ucr/(S_max_cr-sEn_scr)*(sSt_cr-sSt_ucr))
else if (choice_ucr==2) then 
FUN_Mucr = sEn_ucr-sSt_ucr-delT_ucr*(P_ucr-P_ucr*out2_ucr-Kb_ucr*sEn_ucr-Epot_ucr*out_ucr*(sEn_ucr/(sEn_ucr+(sEn_scr-Sm1))))-&
(sEn_ucr/(S_max_cr-sEn_scr)*(sSt_cr-sSt_ucr))
end if 
END FUNCTION FUN_Mucr 

REAL*8 FUNCTION FUN_Mscr(sEn_scr)
REAL*8 :: sEn_scr
REAL*8 :: sSt_cr,sSt_scr,sEn_ucr,P_scr,S_max_cr,S_evmax_cr,S_min_ucr,Epot_scr,Beta_cr,Beta_scr,KQb_scr,KQd_scr,&
delT_scr,out_scr,out1_scr
REAL*8 :: Sm2, Sm3 
COMMON /coeffs/ sSt_cr,sSt_scr,sEn_ucr,P_scr,S_max_cr,S_evmax_cr,S_min_ucr,Epot_scr,Beta_cr,Beta_scr,KQb_scr,&
KQd_scr,delT_scr

 call power((sEn_scr/(S_max_cr-S_min_ucr)),Beta_scr,out1_scr)
 call mkinetic(((sEn_ucr+(sEn_scr-(S_max_cr-S_evmax_cr)))/S_evmax_cr),Beta_cr,out_scr)

Sm2=sEn_scr-(S_max_cr-S_evmax_cr)
Sm3=(S_max_cr-sEn_scr)/(S_max_cr-sEn_scr-sEn_ucr)
FUN_Mscr=sEn_scr-sSt_scr-delT_scr*Sm3*(P_scr-P_scr*out1_scr-KQb_scr*sEn_scr-KQd_scr*sEn_scr-Epot_scr*out_scr*&
(Sm2/(sEn_ucr+Sm2)))
END FUNCTION FUN_Mscr

REAL*8 FUNCTION FUN_Ms(sEn_s)
REAL*8 :: sEn_s
INTEGER:: choice_s
INTEGER:: res_u_1,res_cr_1
REAL*8 :: sSt_s,Kq_s,m_s,alpha_s,P_s,Epot_s,D_FD1,delT_s,out_s,out1_s
COMMON /coeffs/ sSt_s,Kq_s,m_s,alpha_s,P_s,Epot_s,D_FD1,delT_s,choice_s,res_u_1,res_cr_1

call power(sEn_s,alpha_s,out1_s)
call tfun(sEn_s,m_s,out_s)

if (res_u_1==0 .AND. res_cr_1==0 .AND. choice_s==1) then 
FUN_Ms= sEn_s-sSt_s-delT_s*(P_s-Kq_s*out1_s-D_FD1*(Epot_s)*out_s)
else if (res_u_1==0 .AND. res_cr_1==1 .AND. choice_s==1) then 
FUN_Ms= sEn_s-sSt_s-delT_s*(P_s-Kq_s*out1_s)
else if (res_u_1==0 .AND. res_cr_1==0 .AND. choice_s==2) then 
FUN_Ms= sEn_s-sSt_s-delT_s*(P_s-Kq_s*sEn_s-D_FD1*(Epot_s)*out_s)
else if (res_u_1==0 .AND. res_cr_1==1 .AND. choice_s==2) then 
FUN_Ms= sEn_s-sSt_s-delT_s*(P_s-Kq_s*sEn_s)
else if (res_u_1==1 .AND. res_cr_1==0 .AND. choice_s==1) then 
FUN_Ms= sEn_s-sSt_s-delT_s*(P_s-Kq_s*out1_s)
else if (res_u_1==1 .AND. res_cr_1==1 .AND. choice_s==1) then 
FUN_Ms= sEn_s-sSt_s-delT_s*(P_s-Kq_s*out1_s)
else if (res_u_1==1 .AND. res_cr_1==0 .AND. choice_s==2) then 
FUN_Ms= sEn_s-sSt_s-delT_s*(P_s-Kq_s*sEn_s)
else if (res_u_1==1 .AND. res_cr_1==1 .AND. choice_s==2) then 
FUN_Ms= sEn_s-sSt_s-delT_s*(P_s-Kq_s*sEn_s)
end if

END FUNCTION FUN_Ms


SUBROUTINE combined_tank(x1, x2, x3, n, Ce, Cp_WR, m_Q_WR, Kq_WR, Tp_WR, Tm_WR, K_Qq_FR, K_Qb_FR, &
& m_E_FR, alpha_Qq_FR, Beta_Qq_UR, Smax_UR, Beta_E_UR, SiniFr_UR, K_Qb_UR, mu_Qq_UR, Smax_IR, m_QE_IR, &
& K_Qq_RR, Smax_CR, Umax_uCR, Smin_uCR, Beta_Qq_uCR, mu_Qq_uCR, & 
& K_Qb_uCR, Sevmax_CR, Beta_E_CR, Beta_Qq_sCR, K_Qb_sCR, K_Qd_sCR, K_Qq_SR, m_E_SR, &
& alpha_Qq_SR, P_ED_max, m_P_ED, D_S, D_I, D_F, D_R, D_C, Tlag, dT, Lag_RR, Lag_FR, Lag_SR, &
& option_i, option_u, option_f,  option_s,option_c, w_res, i_res, r_res, u_res, f_res, s_res,c_res, output) 

IMPLICIT NONE
REAL*8,EXTERNAL:: FUN_Mw
REAL*8,EXTERNAL:: FUN_Mi
REAL*8,EXTERNAL:: FUN_Mu
REAL*8,EXTERNAL:: FUN_Mf
REAL*8,EXTERNAL:: FUN_Mucr
REAL*8,EXTERNAL:: FUN_Mscr
REAL*8,EXTERNAL:: FUN_Ms

INTEGER :: n, i
INTEGER :: iflag
REAL*8 :: EPS,chkWB,Sw_start,Si_start,Su_start,Ss_start,Sf_start,Sr_start
REAL*8 :: x1(n), x2(n), x3(n), output(n)
REAL*8 :: Ce, K_Qq_FR, K_Qb_FR,  m_E_FR, alpha_Qq_FR
REAL*8 :: Cp_WR, m_Q_WR, Kq_WR, Tp_WR, Tm_WR
REAL*8 :: Beta_Qq_UR, Smax_UR, Beta_E_UR, SiniFr_UR, K_Qb_UR, mu_Qq_UR
REAL*8 :: Smax_IR, m_QE_IR
REAL*8 :: K_Qq_RR
REAL*8 :: Smax_CR, Umax_uCR, Smin_uCR, Beta_Qq_uCR, mu_Qq_uCR, K_Qb_uCR
REAL*8 :: Sevmax_CR, Beta_E_CR, Beta_Qq_sCR, K_Qb_sCR, K_Qd_sCR
REAL*8 :: K_Qq_SR, m_E_SR, alpha_Qq_SR
REAL*8 :: D_S, D_I, D_F, D_R, D_C, P_ED_max, m_P_ED
REAL*8 :: D_FD,D_FD1
REAL*8 :: Tlag, dT
REAL*8 :: sSt_u, Smax_u, P_u, Beta, mu, Gamma, Epot_u,  Kb_u, delT_u, out_u, out1_u, out2_u, out3_u, max_u, min_u !FUN_Mu
REAL*8 :: sSt_r, P_r, K_r, delT_r
REAL*8 :: sSt_f, P_f, Kq_f, Kb_f, m_f, alpha_f, Epot_f, delT_f, out_f , out1_f , max_f, min_f !FUN_Mf
REAL*8 :: sSt_s,Kq_s,m_s,alpha_s,P_s,Epot_s,delT_s,out_s,out1_s,max_s,min_s !FUN_Ms
REAL*8 :: sSt_i, Smax_i, P_i, m_i, Epot_i, delT_i, out_i, out1_i, out2_i, max_i, min_i !FUN_Mi
REAL*8 :: sSt_w, P_w, T, m_w, delT_w, out_w, Cp_w, Kq_w, Tp_w, Tm_w, max_w, min_w !FUN_Mw
REAL*8 :: sSt_cr,sSt_ucr,P_ucr,S_max_cr,S_evmax_cr,U_max_ucr,Epot_ucr,Beta_cr,Beta_ucr,mu_ucr ! FUN_Mucr
REAL*8 :: Kb_ucr,delT_ucr,out_ucr,out1_ucr,out2_ucr, max_ucr, min_ucr,Sm1 ! FUN_Mucr
REAL*8 :: sSt_scr,P_scr,S_min_ucr,Epot_scr,Beta_scr,KQb_scr,KQd_scr,delT_scr,out_scr,out1_scr, max_scr,min_scr,Sm2,Sm3 !FUN_Mscr
INTEGER :: choice_u, choice_f,choice_s, choice_i, choice_ucr,res_w, res_i, res_r, res_u, res_u_1, res_f, res_s, res_cr, res_cr_1
INTEGER :: lag_r, lag_f, lag_s 
INTEGER :: option_i, option_f, option_u, option_s, option_c
INTEGER :: w_res, i_res, r_res, u_res, f_res, s_res, c_res
INTEGER :: Lag_RR, Lag_FR, Lag_SR
REAL*8 :: sEn_w, sEn_i, sEn_u, sEn_f, sEn_s, sEn_scr, sEn_ucr
REAL*8 :: P,Epot !inputs
REAL*8 :: Q_WD, Q_WR, Sw(n) !Snow Reservoir 
REAL*8 :: Q_ID,Q_IR, Si(n), Ei(n)  !Interception Reservoir
REAL*8 :: Sr(n), P_RD, P_r_vec(n), P_rl(n),Q_RR(n)      !Riparian Reservoir
REAL*8 :: Q_SD, P_ED, P_SD, Q_ED,temp          !Infiltration excess overland flow
REAL*8 :: Su(n), Eu(n), Qq_UR, Qb_UR, Q_UR(n) !Unsaturated Reservoir
REAL*8 :: Sf(n), Ef(n), P_fl_a, P_fl_b, P_f_vec(n), P_fl(n), Qq_FR, Qb_FR, Q_FR(n) !Fast Reservoir
REAL*8 :: Sucr(n), Eucr(n), Qq_UCR, Qb_UCR, Qa_CD, Qb_CD,P_SL_vec(n),max_theta(n) ! UCR
REAL*8 :: Sscr(n), Escr(n), Qq_SCR, Qb_SCR, Qd_SCR, Q_CR ! SCR
REAL*8 :: Ss(n), Es(n), Qq_SR, Q_SR(n), P_SL, P_SR(n) !Slow Reservoir
COMMON /coeffi/ sSt_w, P_w, m_w, T, Tm_w, Kq_w, delT_w, out_w
COMMON /coeffi/ sSt_i,Smax_i,P_i,m_i,Epot_i,delT_i,choice_i
COMMON /coeffu/ sSt_u,Smax_u,P_u,Beta,mu,Gamma,Epot_u,Kb_u,delT_u,choice_u
COMMON /coefff/ sSt_f,P_f,Kq_f,Kb_f,m_f,alpha_f,Epot_f,D_FD,delT_f,choice_f,res_u,res_s,res_cr
COMMON /coeffs/ sSt_s,Kq_s,m_s,alpha_s,P_s,Epot_s,D_FD1,delT_s,choice_s,res_u_1,res_cr_1
COMMON /coeffs/ sSt_ucr,sEn_scr,P_ucr,S_max_cr,S_evmax_cr,U_max_ucr,Epot_ucr,Beta_cr,Beta_ucr,mu_ucr,Kb_ucr,delT_ucr,choice_ucr
COMMON /coeffs/ sSt_cr,sSt_scr,sEn_ucr,P_scr,S_min_ucr,Epot_scr,Beta_scr,KQb_scr,KQd_scr,delT_scr


!fast
Kq_f=K_Qq_FR
Kb_f=K_Qb_FR
m_f=m_E_FR
alpha_f=alpha_Qq_FR
lag_f=Lag_FR

!Unsaturated
Beta=Beta_Qq_UR
Smax_u=Smax_UR
Gamma=Beta_E_UR
Kb_u=K_Qb_UR
mu=mu_Qq_UR

!Interception
Smax_i=Smax_IR
m_i=m_QE_IR

!Snow
Cp_w=Cp_WR   
Tp_w=Tp_WR
Tm_w=Tm_WR
Kq_w=Kq_WR
m_w=m_Q_WR

!Riparian
K_r=K_Qq_RR
lag_r=Lag_RR

!UCR
S_max_cr=Smax_CR
U_max_ucr=Umax_uCR 
S_min_ucr=Smin_uCR
Beta_ucr=Beta_Qq_uCR
mu_ucr=mu_Qq_uCR 
Kb_ucr=K_Qb_uCR

!SCR
S_evmax_cr=Sevmax_CR
Beta_cr=Beta_E_CR
Beta_scr=Beta_Qq_sCR
KQb_scr=K_Qb_sCR
KQd_scr=K_Qd_sCR

!Slow
Kq_s=K_Qq_SR
m_s=m_E_SR
alpha_s=alpha_Qq_SR
lag_s=Lag_SR

sSt_u=SiniFr_UR*Smax_U
sSt_r=0._8
sSt_f=0._8
sSt_i=0._8
sSt_w=0._8
sSt_s=0._8
sSt_ucr=Smin_uCR
sSt_scr=0._8

Su_start=SiniFr_UR*Smax_U
Sr_start=0._8
Sf_start=0._8
Si_start=0._8
Sw_start=0._8
Ss_start=0._8

delT_i=dT
delT_w=dT
delT_f=dT
delT_s=dT
delT_u=dT
delT_r=dT
delT_ucr=dT
delT_scr=dT

choice_i=option_i
choice_u=option_u
choice_s=option_s
choice_f=option_f
choice_ucr=option_c

res_i=i_res
res_w=w_res
res_f=f_res
res_s=s_res
res_u=u_res
res_r=r_res
res_u_1=u_res
res_cr=c_res
res_cr_1=c_res

D_FD=D_F
D_FD1=D_F

EPS=1e-9_8

do i=1,n
 
 P=x1(i)
 Epot=Ce*x2(i)
 T=X3(i)
 if (T <= Tp_w) then
 P_w=P*Cp_w
 else if (T > Tp_w) then
 P_w=0._8
 end if
 if (res_w==1) then
 max_w=min(sSt_w+P_w*delT_w,10000._8)
 min_w=0._8
 CALL DFZERO(FUN_Mw,min_w,max_w,min_w,EPS,EPS,iflag)
 sSt_w=min_w
 Sw(i)=sSt_w
 if (T <= Tm_w) then
 Q_WR=0._8
 else if (T > Tm_w) then
 call tfun(Sw(i),m_w,out_w)
 Q_WR=Kq_w*(T-Tm_w)*out_w
 end if
 else if (res_w==0) then
 Sw(i)=0._8
 P_w=0._8
 Q_WR=0._8
 end if
 
 if (res_w==0) then
 Q_WD=P
 else if (res_w==1 .AND. T <= Tp_w ) then
 Q_WD=0._8
 else if (res_w==1 .AND. T > Tp_w) then
 Q_WD=P
 end if
 
 P_i=D_I*Q_WD 
 Q_ID=Q_WD-P_i 
 if (res_i==1) then
 Epot_i=Epot
 max_i=min(sSt_i+P_i*delT_i,Smax_i)
 min_i=0._8
 CALL DFZERO(FUN_Mi,min_i,max_i,min_i,EPS,EPS,iflag)
 sSt_i=min_i
 Si(i)=sSt_i
 if(choice_i==1) then
 call power((Si(i)/Smax_i),m_i,out1_i)
 Q_IR=P_i*out1_i
 else if (choice_i==2) then
 call rhfun((Si(i)/Smax_i),m_i,out2_i)
 Q_IR=P_i*out2_i
 end if
 call mkinetic((Si(i)/Smax_i),m_i,out_i)
 Ei(i)=Epot_i*out_i
 else if (res_i==0) then
 Si(i)=0._8
 Q_IR=P_i
 Epot_i=0._8
 Ei(i)=0._8
 end if
 
 P_RD= Q_IR+Q_ID+Q_WR
 P_r=D_R*P_RD
 P_r_vec(i)=P_r
 P_SD=P_RD-P_r_vec(i)
 Q_SD=D_S*P_SD 
 P_ED=P_SD-Q_SD
 temp=P_ED-P_ED_max
 Q_ED=max(0._8,temp)
 Q_ED=Q_ED/(1+exp(i*1/m_P_ED))
 P_u=P_ED-Q_ED  
 P_fl_a=Q_ED+Q_SD 
 
 if (lag_r==1) then
 CALL tWeights(Tlag, P_r_vec, n, P_rl)
 else if (lag_r==0) then
 P_rl=P_r_vec
 end if

 if(res_r==1) then
 !sSt_r=(P_r_vec(i)*delT_r+sSt_r)/(1+K_r*delT_r)
 P_r=P_rl(i)
 sSt_r=(P_r*delT_r+sSt_r)/(1+K_r*delT_r)
 sSt_r=min(sSt_r,10000._8)
 sSt_r=max(sSt_r,0._8)
 Sr(i)=sSt_r
 Q_RR(i)=K_r*Sr(i)
 else if(res_r==0) then
 Sr(i)=0._8
 Q_RR(i)=P_r_vec(i)
 end if
 
 if(res_u==1) then
 Epot_u=Epot-Ce*Ei(i)
 max_u=min(sSt_u+P_u*delT_u,Smax_u)
 min_u=max(sSt_u-Epot_u*delT_u,0._8)
 CALL DFZERO(FUN_Mu,min_u,max_u,min_u,EPS,EPS,iflag)
 sSt_u=min_u
 Su(i)=sSt_u
 call mkinetic((Su(i)/Smax_u),Gamma,out_u)
 Eu(i)=Epot_u*out_u
 if (choice_u==1) then
 call power((Su(i)/Smax_u),Beta,out1_u)
 Qq_UR=P_u*out1_u
 else if (choice_u==2) then
 Qq_UR=P_u*Su(i)/Smax_u
 else if (choice_u==3) then
 call mkinetic((Su(i)/Smax_u),Beta,out2_u)
 Qq_UR=P_u*out2_u
 else if (choice_u==4) then
 call mlc((Su(i)/Smax_u),Beta,mu,out3_u)
 Qq_UR=P_u*out3_u
 end if
 Qb_UR=Kb_u*Su(i)
 Q_UR(i)=D_FD*(Qq_UR)+Qb_UR
 else if (res_u==0) then
 Su_start=0._8
 Su(i)=0._8
 Eu(i)=0._8
 Qq_UR=0._8
 Qb_UR=0._8
 Q_UR=0._8
 end if
 P_fl_b=(1._8-D_FD)*Qq_UR
 P_f=P_fl_a + P_fl_b
 P_f_vec(i) = P_f
 
 if (lag_f==1) then
 CALL tWeights(Tlag, P_f_vec, n, P_fl)
 else if (lag_f==0) then
 P_fl=P_f_vec
 end if

 if (res_f==1) then
   Epot_f=Epot-Ce*Ei(i)
   P_f=P_fl(i)
   max_f=min(sSt_f+P_f*delT_f,10000._8)
   min_f=0._8
  CALL DFZERO(FUN_Mf,min_f,max_f,min_f,EPS,EPS,iflag)
  sSt_f=min_f
  Sf(i)=sSt_f
  call tfun(Sf(i),m_f,out_f)
   if (res_u==0 .AND. res_cr==0 .AND. res_s==1) then
   Ef(i)=(1-D_FD)*Epot_f*out_f
   else if (res_u==0 .AND. res_cr==0 .AND. res_s==0) then
   Ef(i)=Epot_f*out_f
   else if (res_u==1) then
   Ef(i)=0._8
   else if (res_cr==1) then
   Ef(i)=0._8
   end if
   if (choice_f==1) then
   call power(Sf(i),alpha_f,out1_f)
   Qq_FR=Kq_f*out1_f
   else if (choice_f==2) then
   Qq_FR=Kq_f*Sf(i)
   end if
   if (res_s==1) then
   Qb_FR=Kb_f*Sf(i)
   else if (res_s==0) then
   Qb_FR=0._8
   end if 
   Q_FR(i)=Qq_FR
 else if (res_f==0) then
  Sf(i)=0._8
  Ef(i)=0._8
  Qq_FR=0._8
  Qb_FR=0._8
  Q_FR(i)=P_f_vec(i)
 end if

 if(res_u==1) then
 P_ucr=Qb_UR+D_FD*Qq_UR
 else if(res_u==0) then
 P_ucr=P_u
 end if

 if(res_cr==1) then
 Epot_ucr=Epot-Ce*Ei(i)-Ce*Eu(i)
 S_min_ucr=min(S_max_cr,S_min_ucr)
 max_ucr=min(sSt_ucr+P_ucr*delT_ucr,S_max_cr-sSt_scr)
 min_ucr=max(sSt_ucr-Epot_ucr*delT_ucr,S_min_ucr)
 CALL DFZERO(FUN_Mucr,min_ucr,max_ucr,min_ucr,EPS,EPS,iflag)
 sSt_ucr=min_ucr
 Sucr(i)=sSt_ucr
 call mkinetic(((Sucr(i)+(Sscr(i-1)-(S_max_cr-S_evmax_cr)))/S_evmax_cr),Beta_cr,out_ucr)
 Eucr(i)=Epot_ucr*out_ucr*(Sucr(i)/(Sucr(i)+(sSt_scr-(S_max_cr-S_evmax_cr))))
 max_theta(i)=(Sucr(i)/(S_max_cr-sSt_scr))
 max_theta(i)=min(max_theta(i),U_max_ucr)
 if (choice_ucr==1) then
 call power((max_theta(i)/U_max_ucr),Beta_ucr,out1_ucr)
 Qq_UCR=P_ucr*out1_ucr
 else if (choice_ucr==2) then
 call mlc((max_theta(i)/U_max_ucr),Beta_ucr,mu_ucr,out2_ucr)
 Qq_UCR=P_ucr*out2_ucr
 end if
 Qb_UCR=Kb_ucr*Sucr(i)
 Qa_CD=D_C*Qq_UCR
 P_scr=Qa_CD+Qb_UCR
 Epot_scr=Epot-Ce*Ei(i)-Ce*Eu(i)-Ce*Eucr(i)
 max_scr=min(sSt_scr+P_scr*delT_scr,S_max_cr-S_min_ucr)
 min_scr=max(sSt_scr-Epot_scr*delT_scr,0._8)
 CALL DFZERO(FUN_Mscr,min_scr,max_scr,min_scr,EPS,EPS,iflag)
 sSt_scr=min_scr
 Sscr(i)=sSt_scr
 S_evmax_cr=min(S_max_cr,S_evmax_cr)
 call mkinetic(((Sucr(i)+(Sscr(i)-(S_max_cr-S_evmax_cr)))/S_evmax_cr),Beta_cr,out_scr)
 Escr(i)=Epot_scr*out_scr*((Sscr(i)-(S_max_cr-S_evmax_cr))/(Sucr(i)+(Sscr(i)-(S_max_cr-S_evmax_cr))))
 call power((Sscr(i)/(S_max_cr-S_min_ucr)),Beta_scr,out1_scr)
 Qq_SCR=P_scr*out1_scr
 Qb_SCR=KQb_scr*Sscr(i)
 Qd_SCR=KQd_scr*Sscr(i)
 Q_CR=Qb_SCR+Qd_SCR
 Qb_CD=(1-D_C)*Qq_UCR
 P_SL=Qb_CD+Qq_SCR
 else if (res_cr==0) then
 Sucr(i)=0
 Eucr(i)=0
 Qq_uCR=0
 Qb_UCR=0
 Sscr(i)=0
 Escr(i)=0
 Qq_SCR=0
 Qb_SCR=0
 Qd_SCR=0
 Q_CR=0
 P_SL=P_ucr
 end if
 
 P_SL_vec(i)=P_SL
 if (lag_s==1) then
 CALL tWeights(Tlag, P_SL_vec, n, P_SR)
 else if (lag_s==0) then
 P_SR=P_SL_vec
 end if

 if(res_s==1) then
 P_s=P_SR(i)+Qb_FR
 if (res_u_1==1 .AND. choice_s==1) then
 sSt_s=(P_s*delT_s+sSt_s)/(1+Kq_s*delT_s)
 sSt_s=min(sSt_s,10000._8)
 sSt_s=max(sSt_s,0._8)
 else 
 max_s=min(sSt_s+P_s*delT_s,10000._8)
 min_s=0._8
 Epot_s=Epot-Ce*Ei(i)
 call DFZERO(FUN_Ms, min_s, max_s, min_s, EPS, EPS, iflag)
 sSt_s=min_s
 end if
 Ss(i)=sSt_s
 if (res_u_1==0 .AND. res_cr_1==0) then
 call tfun(Ss(i),m_s,out_s)
 Es(i)=D_FD*Epot_s*out_s
 else if (res_u_1==1) then
 Es(i)=0._8
 else if (res_cr_1==1) then
 Es(i)=0._8
 endif
 if(choice_s==1) then
 call power(sSt_s,alpha_s,out1_s)
 Qq_SR=Kq_s*out1_s
 else if (choice_s==2) then
 Qq_SR=Kq_s*Ss(i)
 end if
 Q_SR(i)=Qq_SR
 else if (res_s==0) then
 Ss(i)=0._8
 Es(i)=0._8
 Q_SR(i)=P_s
 Qq_SR=0._8
 end if
end do

output=Q_RR+Q_FR+Q_SR+Q_CR


!chkWB=sum(x1)-sum(output)-sum(uloss)-sum(Ei)-sum(Eu)-sum(Ef)-sum(Es)-(Si(n)+Su(n)+Sf(n)+Ss(n)+Sr(n)) &
!& +Si_start+Su_start+Ss_start+Sf_start+Sr_start

!chkWB=sum(x1)-sum(output)-sum(uloss)-sum(Ei)-sum(Eu)-sum(Ef)-sum(Es)-(Si(n)+Su(n)+Sf(n)+Ss(n)+Sr(n)) &
!& +Si_start+Su_start+Ss_start+Sf_start+Sr_start-sum(LagArr)

!print*,"chkWB",chkWB

END SUBROUTINE combined_tank


