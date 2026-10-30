      function redrot (jin,kin,min,beta)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     This function uses eq. (4.1.23) of Edmonds 
c     to calculate the reduced rotation matrix element 
c     d(j,k,m;beta) = <jk|exp(+i*beta*Jy/hbar)|jm>.
c     ----------------------------------------------------------------- 
c     
      parameter (zero = 0.0d0)
      parameter (half = 0.5d0)
      parameter (one  = 1.0d0)
      parameter (two  = 2.0d0)
c
c     half integer angular momenta
c
      rj=dfloat(jin)
      rk=dfloat(kin)
      rm=dfloat(min)
  
      sj = half*nint(two*rj)
      sk = half*nint(two*rk)
      sm = half*nint(two*rm)
c
c     projection ranges
c
      redrot = zero
      if (sk.gt.sj .or. sk.lt.-sj)  return
      if (sm.gt.sj .or. sm.lt.-sj)  return
      if (mod(sj-sk,one) .ne. zero) return
      if (mod(sj-sm,one) .ne. zero) return      
c
c     reflection symmetries
c      
      if (sk+sm .ge. zero) then
        if (sk-sm .ge. zero) then
          tk = sk
          tm = sm
          isign = 0
        else
          tk = sm
          tm = sk
          isign = sk-sm
        endif
      else
        if (sk-sm .ge. zero) then
          tk = -sm
          tm = -sk
          isign = 0
        else
          tk = -sk
          tm = -sm
          isign = sk-sm
        endif
      endif
c
c     evaluation
c
      n = sj-tk
      ia = tk-tm
      ib = tk+tm
      a = ia
      b = ib
      beta2 = half*beta
      cosb2 = cos(beta2)
      sinb2 = sin(beta2)
      cosb = (cosb2-sinb2)*(cosb2+sinb2)
      d1 = pjacob0(n,a,b,cosb)
      d2 = cosb2**ib*sinb2**ia
      d3 = d1*d2
      d4 = d3*d3
      ti = tm
      do i = 1,ia
         ti = ti+one
         d4 = d4*(sj+ti)/(sj-ti+one)
      enddo
      d4 = sqrt(d4)
      redrot = sign(d4,d3)
      if (mod(isign,2) .ne. 0) redrot = -redrot
      redrot=dsqrt((2.d0*dfloat(jin)+1.d0)/2.d0)*redrot
      return
      end

      function pjacob0 (n,a,b,x)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     Jacobi polynomial p(n,a,b;x)
c     Abramowitz and Stegun eq. (22.7.1)
c     ----------------------------------------------------------------- 
c
      parameter (zero = 0.0d0)
      parameter (half = 0.5d0)
      parameter (one  = 1.0d0)
      parameter (two  = 2.0d0)
c
      if (n .eq. 0) then
        fp = one
      else
        f = one
        apa = a+a
        apb = a+b
        amb = a-b
        apbamb = apb*amb
        apbp1 = apb+one
        apbp2 = apb+two
        onek = zero
        twok = zero
        fp = half*(amb+apbp2*x)
        do k = 1,n-1
          onek = onek+one
          twok = twok+two 
          a1 = (twok+two)*(onek+apbp1)*(twok+apb)
          a2 = (twok+apbp1)*apbamb
          a3 = (twok+apb)*(twok+apbp1)*(twok+apbp2)
          a4 = (twok+apa)*(onek+b)*(twok+apbp2)
          fm = f
          f = fp
          fp = ((a2+a3*x)*f-a4*fm)/a1
        enddo
      endif
      pjacob0 = fp
      return
      end

C -----------------------------------------------------------------     MMC0125
      SUBROUTINE ORTPOL(IPAR,LMAX,M,E,W,LD)                             MMC0126
C                                                                       MMC0127
C THIS SUBROUTINE CALCULATES ZEROES AND WEIGHTS OF                      MMC0128
C ORTHOGONAL POLYNOMIALS USING JACOBI MATRIX METHOD.                    MMC0129
C THIS MATRIX IS SYMMETRIC AND TRIDIAGONAL , HENCE                      MMC0130
C EISPACK SUBROUTINE TQL2 IS USED IN ITS DIAGONALISATION.               MMC0131
C                                                                       MMC0132
C IPAR = selects the type of orthogonal polynomials used in quadrature: MMC0133
C      = 0 NORMALIZED ASSOCIATED LEGENDRE FUNCTIONS (FERRERS')          MMC0134
C      = 1 NORMALIZED HERMITE POLYNOMIALS                               MMC0135
C                                                                       MMC0136
C zeroes are stored in area E, weights in W                             MMC0137
C                                                                       MMC0138
      PARAMETER (NGLMAX=200)                                             MMC0139
      IMPLICIT REAL*8(A-H,O-Z)                                          MMC0140
      COMMON /DGBDT1/ EMAX, ISDIG, SMAX                                 MMC0141
      COMMON /PRINT/ IPRZP, IPRZV                                       MMC0142
C                                                                       MMC0143
c     DIMENSION A(NGLMAX),E(NGLMAX),Z(NGLMAX,NGLMAX)                    MMC0144
c     DIMENSION W(NGLMAX)                                               MMC0145
c
      DIMENSION E(1),W(1),A(NGLMAX),Z(NGLMAX,NGLMAX)
      PI=ACOS(-1.0D0)                                                   MMC0146
      NORT=NGLMAX                                                       MMC0147
C                                                                       MMC0148
C                 ***JACOBI MATRIX***                                   MMC0149
C                                                                       MMC0150
C JACOBI MATRIX FOR NORMALIZED ASSOCIATED LEGENDRE FUNCTIONS            MMC0151
C    m         m           m                                            MMC0152
C x P (x) = a P   (x) + b P   (x)                                       MMC0153
C    L         L-1         L+1                                          MMC0154
C where      a=sqrt((L+m)(L-m)/((2L+1)(2L-1)))                          MMC0155
C and        b = sqrt((L-m+1)(L+1+m)/((2L+1)(2L+3)))                    MMC0156
C                                                                       MMC0157
      A(1)=0.0D0                                                        MMC0158
      IF(IPAR.EQ.0) THEN                                                MMC0159
         M1=M+1                                                         MMC0160
         LD=LMAX-M                                                      MMC0161
         E(LD)=0.0D0                                                    MMC0162
         IL=0                                                           MMC0163
         DO 100 L=M1,LMAX                                               MMC0164
            IL=IL+1                                                     MMC0165
            LL=L-1                                                      MMC0166
            L1=L-M1+1                                                   MMC0167
            RLM=DFLOAT(L-M)*(L+M)                                       MMC0168
            RLL=DFLOAT(2*LL+1)*(2*LL+3)                                 MMC0169
            E(IL)=0.0D0                                                 MMC0170
            A(IL+1)=DSQRT(RLM/RLL)                                      MMC0171
  100    CONTINUE                                                       MMC0172
      ELSE                                                              MMC0173
C                                                                       MMC0174
C JACOBI MATRIX FOR NORMALIZED GAUSS-HERMITE POLYNOMIALS                MMC0175
C                                                                       MMC0176
C x H (x) = SQRT((n+1)/2) H   (x) + SQRT(n/2) H   (x)                   MMC0177
C    n                     n+1                 n-1                      MMC0178
C                                                                       MMC0179
         LD=LMAX                                                        MMC0180
         E(LD)=0.0D0                                                    MMC0181
         DO 200 L=1,LD-1                                                MMC0182
            E(L)=0.0D0                                                  MMC0183
            A(L+1)=DSQRT(L/2.0D0)                                       MMC0184
  200    CONTINUE                                                       MMC0185
      END IF                                                            MMC0186
C                                                                       MMC0187
C ZEROES OF ORTHOGONAL POLYNOMIAL ARE EIGEN-VALUES OF JACOBI MATRIX     MMC0188
C                                                                       MMC0189
      DO 210 I=1,NORT                                                   MMC0190
         Z(I,I)=1.0D0                                                   MMC0191
         DO 210 J=I+1,LD                                                MMC0192
            Z(I,J)=0.0D0                                                MMC0193
            Z(J,I)=0.0D0                                                MMC0194
  210 CONTINUE                                                          MMC0195
      CALL TQL2(NORT,LD,E,A,Z,IERR)                                     MMC0196
      IF(IERR.EQ.0) GOTO 215                                            MMC0197
         WRITE(6,1010) IERR                                             MMC0198
         STOP                                                           MMC0199
  215 CONTINUE                                                          MMC0200
C                                                                       MMC0201
C PRINT OUT ZEROES OF POLYNOMIALS                                       MMC0202
C                                                                       MMC0203
      IF(IPAR.EQ.0) THEN                                                MMC0204
        IF(IPRZP.EQ.0) GOTO 218                                         MMC0205
        WRITE(6,1020) M,LMAX                                            MMC0206
        WRITE(6,1030) (E(I),I=1,LD)                                     MMC0207
        IF(IPRZV.EQ.0) GOTO 218                                         MMC0208
        WRITE(6,1033)                                                   MMC0209
C EIGEN VECTORS ARE STORED AS COLUMNS OF Z-MATRIX                       MMC0210
        DO 217 I=1,LD                                                   MMC0211
           ZN=0.0D0                                                     MMC0212
           DO 216 J=1,LD                                                MMC0213
  216      ZN=ZN+Z(J,I)*Z(J,I)                                          MMC0214
           WRITE(6,1030) (Z(J,I),J=1,LD)                                MMC0215
           WRITE(6,1035) ZN                                             MMC0216
  217   CONTINUE                                                        MMC0217
  218 CONTINUE                                                          MMC0218
      ELSE                                                              MMC0219
        IF(IPRZP.EQ.0) GOTO 230                                         MMC0220
        WRITE(6,1040) LMAX                                              MMC0221
        WRITE(6,1030) (E(I),I=1,LD)                                     MMC0222
        IF(IPRZV.EQ.0) GOTO 230                                         MMC0223
        WRITE(6,1033)                                                   MMC0224
        DO 227 I=1,LD                                                   MMC0225
           ZN=0.0D0                                                     MMC0226
           DO 226 J=1,LD                                                MMC0227
  226      ZN=ZN+Z(J,I)*Z(J,I)                                          MMC0228
           WRITE(6,1030) (Z(J,I),J=1,LD)                                MMC0229
           WRITE(6,1035) ZN                                             MMC0230
  227   CONTINUE                                                        MMC0231
  230 CONTINUE                                                          MMC0232
      ENDIF                                                             MMC0233
C                                                                       MMC0234
C WEIGHTS ARE CALCULATED FROM THE FIRST COMPONENT OF EIGEN              MMC0235
C VECTORS                                                               MMC0236
      IF(IPAR.EQ.0) THEN                                                MMC0237
C WEIGHTS FOR ASSOCIATED LEGENDRE FUNCTIONS                             MMC0238
        M2=2*M                                                          MMC0239
        RF1=RFF(M2)                                                     MMC0240
        RF2=RFF2(2*M-1)                                                 MMC0241
        DO 250 IZ=1,LD                                                  MMC0242
           GL1=(2.0D0*M+1.0D0)*(1-E(IZ)*E(IZ))**M                       MMC0243
           GL=RF2*DSQRT(GL1/(2.0D0*RF1))                                MMC0244
           HELP=Z(1,IZ)/GL                                              MMC0245
           W(IZ)=HELP*HELP                                              MMC0246
  250   CONTINUE                                                        MMC0247
        IF(IPRZP.EQ.0) GOTO 260                                         MMC0248
        WRITE(6,1050)                                                   MMC0249
        WRITE(6,1030) (W(I),I=1,LD)                                     MMC0250
  260   CONTINUE                                                        MMC0251
      ELSE                                                              MMC0252
C WEIGHTS FOR HERMITE POLYNOMIALS : Ho(x)=1/sqrt(sqrt(PI))              MMC0253
        DO 350 IZ=1,LD                                                  MMC0254
           W(IZ)=Z(1,IZ)*Z(1,IZ)*DSQRT(PI)                              MMC0255
  350   CONTINUE                                                        MMC0256
        IF(IPRZP.EQ.0) GOTO 360                                         MMC0257
        WRITE(6,1050)                                                   MMC0258
        WRITE(6,1030) (W(I),I=1,LD)                                     MMC0259
  360   CONTINUE                                                        MMC0260
      ENDIF                                                             MMC0261
 1005 FORMAT(//10X,'ORTHOGONAL POLYNOMIALS'//)                          MMC0262
 1010 FORMAT(/1X,'STOP IN ORTPOL :IN TQL2 IERR =',I3/)                  MMC0263
 1020 FORMAT(//1X,'ZEROES OF ASSOCIATED LEGENDRE FUNCTIONS OF ORDER M ='MMC0264
     *,I3,' AND DEGREE L =',I3/)                                        MMC0265
 1030 FORMAT(5D16.8)                                                    MMC0266
 1033 FORMAT(/1X,'EIGEN VECTORS :'/)                                    MMC0267
 1035 FORMAT(1X,'NORM =',D20.8/)                                        MMC0268
 1040 FORMAT(//1X,'ZEROES OF HERMITE POLYNOMIAL OF ORDER N =',I4/)      MMC0269
 1050 FORMAT(/1X,'CORRESPONDING WEIGHTS'/)                              MMC0270
      RETURN                                                            MMC0271
      END                                                               MMC0272
C ----------------------------------------------------------------      MMC0273
      REAL*8 FUNCTION RFF(M)                                            MMC0274
      RFF=1.0D0                                                         MMC0275
      IF(M.LE.1) THEN                                                   MMC0276
        RETURN                                                          MMC0277
      ELSE                                                              MMC0278
        DO 100 I=2,M                                                    MMC0279
           RFF=RFF*I                                                    MMC0280
  100   CONTINUE                                                        MMC0281
      RETURN                                                            MMC0282
      ENDIF                                                             MMC0283
      END                                                               MMC0284
C ----------------------------------------------------------------      MMC0285
      REAL*8 FUNCTION RFF2(M)                                           MMC0286
      RFF2=1.0D0                                                        MMC0287
      IF(M.LE.1) THEN                                                   MMC0288
        RETURN                                                          MMC0289
      ELSE                                                              MMC0290
      DO 100 I=1,M,2                                                    MMC0291
         RFF2=RFF2*I                                                    MMC0292
  100 CONTINUE                                                          MMC0293
      END IF                                                            MMC0294
      RETURN                                                            MMC0295
      END                                                               MMC0296
C -----------------------------------------------------------------
      FUNCTION BRENT(AX,BX,CX,F,TOL,XMIN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL F
      PARAMETER (ITMAX=100,CGOLD=0.3819660D0,ZEPS=1.0D-10)
      A=MIN(AX,CX)
      B=MAX(AX,CX)
      V=BX
      W=V
      X=V
      E=0.0D0
      FX=F(X)
      FV=FX
      FW=FX
      DO 180,ITER=1,ITMAX
	XM=0.5D0*(A+B)
	TOL1=TOL*ABS(X)+ZEPS
	TOL2=2.0D0*TOL1
	IF (ABS(X-XM).LE.(TOL2-0.5D0*(B-A))) GO TO 3
	IF (ABS(E).GT.TOL1) THEN
	  R=(X-W)*(FX-FV)
	  Q=(X-V)*(FX-FW)
	  P=(X-V)*Q-(X-W)*R
	  Q=2.0D0*(Q-R)
	  IF (Q.GT.0.0D0) P=-P
	  Q=ABS(Q)
	  ETEMP=E
	  E=D
	  IF (ABS(P).GE.ABS(0.5D0*Q*ETEMP).OR.P.LE.Q*(A-X).OR.
     -        P.GE.Q*(B-X)) GO TO 1
	  D=P/Q
	  U=X+D
	  IF (U-A.LT.TOL2 .OR. B-U.LT.TOL2) D=SIGN(TOL1,XM-X)
	  GO TO 2
	END IF
1       IF (X.GE.XM) THEN
	  E=A-X
	ELSE
	  E=B-X
	END IF
	D=CGOLD*E
2       IF (ABS(D).GE.TOL1) THEN
	  U=X+D
	ELSE
	  U=X+SIGN(TOL1,D)
	END IF
	FU=F(U)
	IF (FU.LE.FX) THEN
	  IF (U.GE.FX) THEN
	    A=X
	  ELSE
	    B=X
	  END IF
	  V=W
	  FV=FW
	  W=X
	  FW=FX
	  X=U
	  FX=FU
	ELSE
	  IF (U.LT.X) THEN
	    A=U
	  ELSE
	    B=U
	  END IF
	  IF (FU.LE.FW .OR. W.EQ.X) THEN
	    V=W
	    FV=FW
	    W=U
	    FW=FU
	  ELSE IF (FU.LE.FV .OR. V.EQ.X .OR. V.EQ.W) THEN
	    V=U
	    FV=FU
	  END IF
	END IF
180   CONTINUE
      PAUSE 'BRENT EXCEEDED MAXIMUM ITERATIONS.'
3     XMIN=X
      BRENT=FX
      RETURN
      END
c-----------------------------------------------------------------
      SUBROUTINE RBES(N,Z,ZJN,ZJNP,ZYN,ZYNP)
      IMPLICIT REAL*8 (A-H,O-Z)
C***********************************************************************
C***  THIS SUBROUTINE CALCULATES THE REGULAR RICATTI-BESSEL FUNCTION ZJN
C***  AND ITS DERIVATIVE ZJNP AND THE IRREGULAR RICATTI-BESSEL FUNCTION
C***  ZYN AND ITS DERIVATIVE ZYNP OF ORDER N FOR ARGUMENT Z.  THIS IS
C***  A COPY OF LIGHTS ROUTINE ATTRIBUTED TO R.GORDON.
C***********************************************************************
      DATA ZERO,HALF,ONE,TWO,THREE/0.D0,0.5D0,1.D0,2.D0,3.D0/
C***  EVALUATE TRIG FUNCTIONS
      SZ=DSIN(Z)
      CZ=DCOS(Z)
      ZINV=ONE/Z
      AJ=SZ
      BJ=SZ*ZINV-CZ
      AY=-CZ
      BY=-SZ-CZ*ZINV
      IF(N-1)1,2,3
C***  FUNCTIONS OF ZERO ORDER
    1 ZJN=AJ
      ZJNP=CZ
      ZYN=AY
      ZYNP=SZ
      RETURN
C***  FUNCTIONS OF ORDER ONE
    2 ZJN=BJ
      ZJNP=AJ-ZINV*BJ
      ZYN=BY
      ZYNP=AY-ZINV*BY
      RETURN
C***  FUNCTIONS OF ORDER GREATER THAN ONE
    3 XN=N
      DELF=ZINV+ZINV
      FACTOR=DELF+ZINV
C***  TEST TO SEE IF RECURRENCE CAN BE USED
      IF(Z.GT.XN)GO TO 100
C***  IN NON CLASSICAL REGION USE RECURSION RELATION FOR IRREGULAR
C***  SOLUTION AND ASCENDING POWER SERIES FOR REGULAR SOLUTION.
C***  THE COEFFICIENT,C, IN FRONT OF POWER SERIES IS CALCULATED DURING
C***  FORWARD RECURRENCE.
      ZSQ=Z*Z
      TOP=XN+XN
      XJ=THREE
      C=ZSQ/THREE
C***  IRREGULAR SOLUTION
   25 ZYN=BY*FACTOR-AY
      AY=BY
      BY=ZYN
      FACTOR=FACTOR+DELF
      XJ=XJ+TWO
      C=C*(Z/XJ)
      IF(XJ.LT.TOP)GO TO 25
      ZYNP=AY-ZINV*BY*XN
C***  REGULAR SOLUTION
      FACTOR=-HALF*ZSQ
      U=ONE
      TERM=ONE
      DU=ONE
      DTERM=ONE
      D2NP3=TOP+THREE
      DEN=D2NP3
      DFACT=ZERO
   35 DFACT=DFACT+ONE
      TERM=TERM*(FACTOR/(DFACT*DEN))
      U=U+TERM
      DEN=DEN+TWO
      DTERM=DTERM*(FACTOR/(DFACT*DEN))
      DU=DU+DTERM
C***  TEST CONVERGENCE TO SINGLE PRECISION ACCURACY
      IF(ABS(TERM).GT.1.0D-9)GO TO 35
      ZJN=U*C
      ZJNP=(XN+ONE)*ZJN*ZINV-(Z*DU/D2NP3)*C
      RETURN
C***  FOR CLASSICAL CASE USE FORWARD RECURSION FOR REGULAR
C***  AND IRREGULAR SOLUTIONS
  100 CONST=ZINV*XN
      TOP=CONST+CONST
  200 AJ=FACTOR*BJ-AJ
      AY=FACTOR*BY-AY
      FACTOR=FACTOR+DELF
      IF(FACTOR.GT.TOP)GO TO 250
      BJ=FACTOR*AJ-BJ
      BY=FACTOR*AY-BY
      FACTOR=FACTOR+DELF
      IF(FACTOR.LT.TOP)GO TO 200
      ZJN=BJ
      ZJNP=AJ-CONST*BJ
      ZYN=BY
      ZYNP=AY-CONST*BY
      RETURN
  250 ZJN=AJ
      ZJNP=BJ-CONST*AJ
      ZYN=AY
      ZYNP=BY-CONST*AY
      END
      DOUBLE PRECISION FUNCTION EPSLON (X)                              EPS00010
      DOUBLE PRECISION X                                                EPS00020
C                                                                       EPS00030
C     ESTIMATE UNIT ROUNDOFF IN QUANTITIES OF SIZE X.                   EPS00040
C                                                                       EPS00050
      DOUBLE PRECISION A,B,C,EPS                                        EPS00060
C                                                                       EPS00070
C     THIS PROGRAM SHOULD FUNCTION PROPERLY ON ALL SYSTEMS              EPS00080
C     SATISFYING THE FOLLOWING TWO ASSUMPTIONS,                         EPS00090
C        1.  THE BASE USED IN REPRESENTING FLOATING POINT               EPS00100
C            NUMBERS IS NOT A POWER OF THREE.                           EPS00110
C        2.  THE QUANTITY  A  IN STATEMENT 10 IS REPRESENTED TO         EPS00120
C            THE ACCURACY USED IN FLOATING POINT VARIABLES              EPS00130
C            THAT ARE STORED IN MEMORY.                                 EPS00140
C     THE STATEMENT NUMBER 10 AND THE GO TO 10 ARE INTENDED TO          EPS00150
C     FORCE OPTIMIZING COMPILERS TO GENERATE CODE SATISFYING            EPS00160
C     ASSUMPTION 2.                                                     EPS00170
C     UNDER THESE ASSUMPTIONS, IT SHOULD BE TRUE THAT,                  EPS00180
C            A  IS NOT EXACTLY EQUAL TO FOUR-THIRDS,                    EPS00190
C            B  HAS A ZERO FOR ITS LAST BIT OR DIGIT,                   EPS00200
C            C  IS NOT EXACTLY EQUAL TO ONE,                            EPS00210
C            EPS  MEASURES THE SEPARATION OF 1.0 FROM                   EPS00220
C                 THE NEXT LARGER FLOATING POINT NUMBER.                EPS00230
C     THE DEVELOPERS OF EISPACK WOULD APPRECIATE BEING INFORMED         EPS00240
C     ABOUT ANY SYSTEMS WHERE THESE ASSUMPTIONS DO NOT HOLD.            EPS00250
C                                                                       EPS00260
C     THIS VERSION DATED 4/6/83.                                        EPS00270
C                                                                       EPS00280
      A = 4.0D0/3.0D0                                                   EPS00290
   10 B = A - 1.0D0                                                     EPS00300
      C = B + B + B                                                     EPS00310
      EPS = DABS(C-1.0D0)                                               EPS00320
      IF (EPS .EQ. 0.0D0) GO TO 10                                      EPS00330
      EPSLON = EPS*DABS(X)                                              EPS00340
      RETURN                                                            EPS00350
      END                                                               EPS00360
      DOUBLE PRECISION FUNCTION PYTHAG(A,B)                             PYT00010
      DOUBLE PRECISION A,B                                              PYT00020
C                                                                       PYT00030
C     FINDS DSQRT(A**2+B**2) WITHOUT OVERFLOW OR DESTRUCTIVE UNDERFLOW  PYT00040
C                                                                       PYT00050
      DOUBLE PRECISION P,R,S,T,U                                        PYT00060
      P = DMAX1(DABS(A),DABS(B))                                        PYT00070
      IF (P .EQ. 0.0D0) GO TO 20                                        PYT00080
      R = (DMIN1(DABS(A),DABS(B))/P)**2                                 PYT00090
   10 CONTINUE                                                          PYT00100
         T = 4.0D0 + R                                                  PYT00110
         IF (T .EQ. 4.0D0) GO TO 20                                     PYT00120
         S = R/T                                                        PYT00130
         U = 1.0D0 + 2.0D0*S                                            PYT00140
         P = U*P                                                        PYT00150
         R = (S/U)**2 * R                                               PYT00160
      GO TO 10                                                          PYT00170
   20 PYTHAG = P                                                        PYT00180
      RETURN                                                            PYT00190
      END                                                               PYT00200
      SUBROUTINE RS(NM,N,A,W,MATZ,Z,FV1,FV2,IERR)                       RS 00010
C                                                                       RS 00020
      INTEGER N,NM,IERR,MATZ                                            RS 00030
      DOUBLE PRECISION A(NM,N),W(N),Z(NM,N),FV1(N),FV2(N)               RS 00040
C                                                                       RS 00050
C     THIS SUBROUTINE CALLS THE RECOMMENDED SEQUENCE OF                 RS 00060
C     SUBROUTINES FROM THE EIGENSYSTEM SUBROUTINE PACKAGE (EISPACK)     RS 00070
C     TO FIND THE EIGENVALUES AND EIGENVECTORS (IF DESIRED)             RS 00080
C     OF A REAL SYMMETRIC MATRIX.                                       RS 00090
C                                                                       RS 00100
C     ON INPUT                                                          RS 00110
C                                                                       RS 00120
C        NM  MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL    RS 00130
C        ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM            RS 00140
C        DIMENSION STATEMENT.                                           RS 00150
C                                                                       RS 00160
C        N  IS THE ORDER OF THE MATRIX  A.                              RS 00170
C                                                                       RS 00180
C        A  CONTAINS THE REAL SYMMETRIC MATRIX.                         RS 00190
C                                                                       RS 00200
C        MATZ  IS AN INTEGER VARIABLE SET EQUAL TO ZERO IF              RS 00210
C        ONLY EIGENVALUES ARE DESIRED.  OTHERWISE IT IS SET TO          RS 00220
C        ANY NON-ZERO INTEGER FOR BOTH EIGENVALUES AND EIGENVECTORS.    RS 00230
C                                                                       RS 00240
C     ON OUTPUT                                                         RS 00250
C                                                                       RS 00260
C        W  CONTAINS THE EIGENVALUES IN ASCENDING ORDER.                RS 00270
C                                                                       RS 00280
C        Z  CONTAINS THE EIGENVECTORS IF MATZ IS NOT ZERO.              RS 00290
C                                                                       RS 00300
C        IERR  IS AN INTEGER OUTPUT VARIABLE SET EQUAL TO AN ERROR      RS 00310
C           COMPLETION CODE DESCRIBED IN THE DOCUMENTATION FOR TQLRAT   RS 00320
C           AND TQL2.  THE NORMAL COMPLETION CODE IS ZERO.              RS 00330
C                                                                       RS 00340
C        FV1  AND  FV2  ARE TEMPORARY STORAGE ARRAYS.                   RS 00350
C                                                                       RS 00360
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,    RS 00370
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY RS 00380
C                                                                       RS 00390
C     THIS VERSION DATED AUGUST 1983.                                   RS 00400
C                                                                       RS 00410
C     ------------------------------------------------------------------RS 00420
C                                                                       RS 00430
      IF (N .LE. NM) GO TO 10                                           RS 00440
      IERR = 10 * N                                                     RS 00450
      GO TO 50                                                          RS 00460
C                                                                       RS 00470
   10 IF (MATZ .NE. 0) GO TO 20                                         RS 00480
C     .......... FIND EIGENVALUES ONLY ..........                       RS 00490
      CALL  TRED1(NM,N,A,W,FV1,FV2)                                     RS 00500
      CALL  TQLRAT(N,W,FV2,IERR)                                        RS 00510
      GO TO 50                                                          RS 00520
C     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........      RS 00530
   20 CALL  TRED2(NM,N,A,W,FV1,Z)                                       RS 00540
      CALL  TQL2(NM,N,W,FV1,Z,IERR)                                     RS 00550
   50 RETURN                                                            RS 00560
      END                                                               RS 00570
      SUBROUTINE TQL2(NM,N,D,E,Z,IERR)                                  TQL00010
C                                                                       TQL00020
      INTEGER I,J,K,L,M,N,II,L1,L2,NM,MML,IERR                          TQL00030
      DOUBLE PRECISION D(N),E(N),Z(NM,N)                                TQL00040
      DOUBLE PRECISION C,C2,C3,DL1,EL1,F,G,H,P,R,S,S2,TST1,TST2,PYTHAG  TQL00050
C                                                                       TQL00060
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQL2,     TQL00070
C     NUM. MATH. 11, 293-306(1968) BY BOWDLER, MARTIN, REINSCH, AND     TQL00080
C     WILKINSON.                                                        TQL00090
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).   TQL00100
C                                                                       TQL00110
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS            TQL00120
C     OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE QL METHOD.               TQL00130
C     THE EIGENVECTORS OF A FULL SYMMETRIC MATRIX CAN ALSO              TQL00140
C     BE FOUND IF  TRED2  HAS BEEN USED TO REDUCE THIS                  TQL00150
C     FULL MATRIX TO TRIDIAGONAL FORM.                                  TQL00160
C                                                                       TQL00170
C     ON INPUT                                                          TQL00180
C                                                                       TQL00190
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL         TQL00200
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM          TQL00210
C          DIMENSION STATEMENT.                                         TQL00220
C                                                                       TQL00230
C        N IS THE ORDER OF THE MATRIX.                                  TQL00240
C                                                                       TQL00250
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.          TQL00260
C                                                                       TQL00270
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX        TQL00280
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.               TQL00290
C                                                                       TQL00300
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE           TQL00310
C          REDUCTION BY  TRED2, IF PERFORMED.  IF THE EIGENVECTORS      TQL00320
C          OF THE TRIDIAGONAL MATRIX ARE DESIRED, Z MUST CONTAIN        TQL00330
C          THE IDENTITY MATRIX.                                         TQL00340
C                                                                       TQL00350
C      ON OUTPUT                                                        TQL00360
C                                                                       TQL00370
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN          TQL00380
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT          TQL00390
C          UNORDERED FOR INDICES 1,2,...,IERR-1.                        TQL00400
C                                                                       TQL00410
C        E HAS BEEN DESTROYED.                                          TQL00420
C                                                                       TQL00430
C        Z CONTAINS ORTHONORMAL EIGENVECTORS OF THE SYMMETRIC           TQL00440
C          TRIDIAGONAL (OR FULL) MATRIX.  IF AN ERROR EXIT IS MADE,     TQL00450
C          Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED       TQL00460
C          EIGENVALUES.                                                 TQL00470
C                                                                       TQL00480
C        IERR IS SET TO                                                 TQL00490
C          ZERO       FOR NORMAL RETURN,                                TQL00500
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN               TQL00510
C                     DETERMINED AFTER 30 ITERATIONS.                   TQL00520
C                                                                       TQL00530
C     CALLS PYTHAG FOR  DSQRT(A*A + B*B) .                              TQL00540
C                                                                       TQL00550
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,    TQL00560
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY TQL00570
C                                                                       TQL00580
C     THIS VERSION DATED AUGUST 1983.                                   TQL00590
C                                                                       TQL00600
C     ------------------------------------------------------------------TQL00610
C                                                                       TQL00620
      IERR = 0                                                          TQL00630
      IF (N .EQ. 1) GO TO 1001                                          TQL00640
C                                                                       TQL00650
      DO 100 I = 2, N                                                   TQL00660
  100 E(I-1) = E(I)                                                     TQL00670
C                                                                       TQL00680
      F = 0.0D0                                                         TQL00690
      TST1 = 0.0D0                                                      TQL00700
      E(N) = 0.0D0                                                      TQL00710
C                                                                       TQL00720
      DO 240 L = 1, N                                                   TQL00730
         J = 0                                                          TQL00740
         H = DABS(D(L)) + DABS(E(L))                                    TQL00750
         IF (TST1 .LT. H) TST1 = H                                      TQL00760
C     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........         TQL00770
         DO 110 M = L, N                                                TQL00780
            TST2 = TST1 + DABS(E(M))                                    TQL00790
            IF (TST2 .EQ. TST1) GO TO 120                               TQL00800
C     .......... E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT               TQL00810
C                THROUGH THE BOTTOM OF THE LOOP ..........              TQL00820
  110    CONTINUE                                                       TQL00830
C                                                                       TQL00840
  120    IF (M .EQ. L) GO TO 220                                        TQL00850
  130    IF (J .EQ. 30) GO TO 1000                                      TQL00860
         J = J + 1                                                      TQL00870
C     .......... FORM SHIFT ..........                                  TQL00880
         L1 = L + 1                                                     TQL00890
         L2 = L1 + 1                                                    TQL00900
         G = D(L)                                                       TQL00910
         P = (D(L1) - G) / (2.0D0 * E(L))                               TQL00920
         R = PYTHAG(P,1.0D0)                                            TQL00930
         D(L) = E(L) / (P + DSIGN(R,P))                                 TQL00940
         D(L1) = E(L) * (P + DSIGN(R,P))                                TQL00950
         DL1 = D(L1)                                                    TQL00960
         H = G - D(L)                                                   TQL00970
         IF (L2 .GT. N) GO TO 145                                       TQL00980
C                                                                       TQL00990
         DO 140 I = L2, N                                               TQL01000
  140    D(I) = D(I) - H                                                TQL01010
C                                                                       TQL01020
  145    F = F + H                                                      TQL01030
C     .......... QL TRANSFORMATION ..........                           TQL01040
         P = D(M)                                                       TQL01050
         C = 1.0D0                                                      TQL01060
         C2 = C                                                         TQL01070
         EL1 = E(L1)                                                    TQL01080
         S = 0.0D0                                                      TQL01090
         MML = M - L                                                    TQL01100
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........             TQL01110
         DO 200 II = 1, MML                                             TQL01120
            C3 = C2                                                     TQL01130
            C2 = C                                                      TQL01140
            S2 = S                                                      TQL01150
            I = M - II                                                  TQL01160
            G = C * E(I)                                                TQL01170
            H = C * P                                                   TQL01180
            R = PYTHAG(P,E(I))                                          TQL01190
            E(I+1) = S * R                                              TQL01200
            S = E(I) / R                                                TQL01210
            C = P / R                                                   TQL01220
            P = C * D(I) - S * G                                        TQL01230
            D(I+1) = H + S * (C * G + S * D(I))                         TQL01240
C     .......... FORM VECTOR ..........                                 TQL01250
            DO 180 K = 1, N                                             TQL01260
               H = Z(K,I+1)                                             TQL01270
               Z(K,I+1) = S * Z(K,I) + C * H                            TQL01280
               Z(K,I) = C * Z(K,I) - S * H                              TQL01290
  180       CONTINUE                                                    TQL01300
C                                                                       TQL01310
  200    CONTINUE                                                       TQL01320
C                                                                       TQL01330
         P = -S * S2 * C3 * EL1 * E(L) / DL1                            TQL01340
         E(L) = S * P                                                   TQL01350
         D(L) = C * P                                                   TQL01360
         TST2 = TST1 + DABS(E(L))                                       TQL01370
         IF (TST2 .GT. TST1) GO TO 130                                  TQL01380
  220    D(L) = D(L) + F                                                TQL01390
  240 CONTINUE                                                          TQL01400
C     .......... ORDER EIGENVALUES AND EIGENVECTORS ..........          TQL01410
      DO 300 II = 2, N                                                  TQL01420
         I = II - 1                                                     TQL01430
         K = I                                                          TQL01440
         P = D(I)                                                       TQL01450
C                                                                       TQL01460
         DO 260 J = II, N                                               TQL01470
            IF (D(J) .GE. P) GO TO 260                                  TQL01480
            K = J                                                       TQL01490
            P = D(J)                                                    TQL01500
  260    CONTINUE                                                       TQL01510
C                                                                       TQL01520
         IF (K .EQ. I) GO TO 300                                        TQL01530
         D(K) = D(I)                                                    TQL01540
         D(I) = P                                                       TQL01550
C                                                                       TQL01560
         DO 280 J = 1, N                                                TQL01570
            P = Z(J,I)                                                  TQL01580
            Z(J,I) = Z(J,K)                                             TQL01590
            Z(J,K) = P                                                  TQL01600
  280    CONTINUE                                                       TQL01610
C                                                                       TQL01620
  300 CONTINUE                                                          TQL01630
C                                                                       TQL01640
      GO TO 1001                                                        TQL01650
C     .......... SET ERROR -- NO CONVERGENCE TO AN                      TQL01660
C                EIGENVALUE AFTER 30 ITERATIONS ..........              TQL01670
 1000 IERR = L                                                          TQL01680
 1001 RETURN                                                            TQL01690
      END                                                               TQL01700
      SUBROUTINE TQLRAT(N,D,E2,IERR)                                    TQL00010
C                                                                       TQL00020
      INTEGER I,J,L,M,N,II,L1,MML,IERR                                  TQL00030
      DOUBLE PRECISION D(N),E2(N)                                       TQL00040
      DOUBLE PRECISION B,C,F,G,H,P,R,S,T,EPSLON,PYTHAG                  TQL00050
C                                                                       TQL00060
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQLRAT,   TQL00070
C     ALGORITHM 464, COMM. ACM 16, 689(1973) BY REINSCH.                TQL00080
C                                                                       TQL00090
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC              TQL00100
C     TRIDIAGONAL MATRIX BY THE RATIONAL QL METHOD.                     TQL00110
C                                                                       TQL00120
C     ON INPUT                                                          TQL00130
C                                                                       TQL00140
C        N IS THE ORDER OF THE MATRIX.                                  TQL00150
C                                                                       TQL00160
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.          TQL00170
C                                                                       TQL00180
C        E2 CONTAINS THE SQUARES OF THE SUBDIAGONAL ELEMENTS OF THE     TQL00190
C          INPUT MATRIX IN ITS LAST N-1 POSITIONS.  E2(1) IS ARBITRARY. TQL00200
C                                                                       TQL00210
C      ON OUTPUT                                                        TQL00220
C                                                                       TQL00230
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN          TQL00240
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND          TQL00250
C          ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE            TQL00260
C          THE SMALLEST EIGENVALUES.                                    TQL00270
C                                                                       TQL00280
C        E2 HAS BEEN DESTROYED.                                         TQL00290
C                                                                       TQL00300
C        IERR IS SET TO                                                 TQL00310
C          ZERO       FOR NORMAL RETURN,                                TQL00320
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN               TQL00330
C                     DETERMINED AFTER 30 ITERATIONS.                   TQL00340
C                                                                       TQL00350
C     CALLS PYTHAG FOR  DSQRT(A*A + B*B) .                              TQL00360
C                                                                       TQL00370
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,    TQL00380
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY TQL00390
C                                                                       TQL00400
C     THIS VERSION DATED AUGUST 1983.                                   TQL00410
C                                                                       TQL00420
C     ------------------------------------------------------------------TQL00430
C                                                                       TQL00440
      IERR = 0                                                          TQL00450
      IF (N .EQ. 1) GO TO 1001                                          TQL00460
C                                                                       TQL00470
      DO 100 I = 2, N                                                   TQL00480
  100 E2(I-1) = E2(I)                                                   TQL00490
C                                                                       TQL00500
      F = 0.0D0                                                         TQL00510
      T = 0.0D0                                                         TQL00520
      E2(N) = 0.0D0                                                     TQL00530
C                                                                       TQL00540
      DO 290 L = 1, N                                                   TQL00550
         J = 0                                                          TQL00560
         H = DABS(D(L)) + DSQRT(E2(L))                                  TQL00570
         IF (T .GT. H) GO TO 105                                        TQL00580
         T = H                                                          TQL00590
         B = EPSLON(T)                                                  TQL00600
         C = B * B                                                      TQL00610
C     .......... LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT .......... TQL00620
  105    DO 110 M = L, N                                                TQL00630
            IF (E2(M) .LE. C) GO TO 120                                 TQL00640
C     .......... E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT              TQL00650
C                THROUGH THE BOTTOM OF THE LOOP ..........              TQL00660
  110    CONTINUE                                                       TQL00670
C                                                                       TQL00680
  120    IF (M .EQ. L) GO TO 210                                        TQL00690
  130    IF (J .EQ. 30) GO TO 1000                                      TQL00700
         J = J + 1                                                      TQL00710
C     .......... FORM SHIFT ..........                                  TQL00720
         L1 = L + 1                                                     TQL00730
         S = DSQRT(E2(L))                                               TQL00740
         G = D(L)                                                       TQL00750
         P = (D(L1) - G) / (2.0D0 * S)                                  TQL00760
         R = PYTHAG(P,1.0D0)                                            TQL00770
         D(L) = S / (P + DSIGN(R,P))                                    TQL00780
         H = G - D(L)                                                   TQL00790
C                                                                       TQL00800
         DO 140 I = L1, N                                               TQL00810
  140    D(I) = D(I) - H                                                TQL00820
C                                                                       TQL00830
         F = F + H                                                      TQL00840
C     .......... RATIONAL QL TRANSFORMATION ..........                  TQL00850
         G = D(M)                                                       TQL00860
         IF (G .EQ. 0.0D0) G = B                                        TQL00870
         H = G                                                          TQL00880
         S = 0.0D0                                                      TQL00890
         MML = M - L                                                    TQL00900
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........             TQL00910
         DO 200 II = 1, MML                                             TQL00920
            I = M - II                                                  TQL00930
            P = G * H                                                   TQL00940
            R = P + E2(I)                                               TQL00950
            E2(I+1) = S * R                                             TQL00960
            S = E2(I) / R                                               TQL00970
            D(I+1) = H + S * (H + D(I))                                 TQL00980
            G = D(I) - E2(I) / G                                        TQL00990
            IF (G .EQ. 0.0D0) G = B                                     TQL01000
            H = G * P / R                                               TQL01010
  200    CONTINUE                                                       TQL01020
C                                                                       TQL01030
         E2(L) = S * G                                                  TQL01040
         D(L) = H                                                       TQL01050
C     .......... GUARD AGAINST UNDERFLOW IN CONVERGENCE TEST .......... TQL01060
         IF (H .EQ. 0.0D0) GO TO 210                                    TQL01070
         IF (DABS(E2(L)) .LE. DABS(C/H)) GO TO 210                      TQL01080
         E2(L) = H * E2(L)                                              TQL01090
         IF (E2(L) .NE. 0.0D0) GO TO 130                                TQL01100
  210    P = D(L) + F                                                   TQL01110
C     .......... ORDER EIGENVALUES ..........                           TQL01120
         IF (L .EQ. 1) GO TO 250                                        TQL01130
C     .......... FOR I=L STEP -1 UNTIL 2 DO -- ..........               TQL01140
         DO 230 II = 2, L                                               TQL01150
            I = L + 2 - II                                              TQL01160
            IF (P .GE. D(I-1)) GO TO 270                                TQL01170
            D(I) = D(I-1)                                               TQL01180
  230    CONTINUE                                                       TQL01190
C                                                                       TQL01200
  250    I = 1                                                          TQL01210
  270    D(I) = P                                                       TQL01220
  290 CONTINUE                                                          TQL01230
C                                                                       TQL01240
      GO TO 1001                                                        TQL01250
C     .......... SET ERROR -- NO CONVERGENCE TO AN                      TQL01260
C                EIGENVALUE AFTER 30 ITERATIONS ..........              TQL01270
 1000 IERR = L                                                          TQL01280
 1001 RETURN                                                            TQL01290
      END                                                               TQL01300
      SUBROUTINE TRED1(NM,N,A,D,E,E2)                                   TRE00010
C                                                                       TRE00020
      INTEGER I,J,K,L,N,II,NM,JP1                                       TRE00030
      DOUBLE PRECISION A(NM,N),D(N),E(N),E2(N)                          TRE00040
      DOUBLE PRECISION F,G,H,SCALE                                      TRE00050
C                                                                       TRE00060
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED1,    TRE00070
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.   TRE00080
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).   TRE00090
C                                                                       TRE00100
C     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX                   TRE00110
C     TO A SYMMETRIC TRIDIAGONAL MATRIX USING                           TRE00120
C     ORTHOGONAL SIMILARITY TRANSFORMATIONS.                            TRE00130
C                                                                       TRE00140
C     ON INPUT                                                          TRE00150
C                                                                       TRE00160
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL         TRE00170
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM          TRE00180
C          DIMENSION STATEMENT.                                         TRE00190
C                                                                       TRE00200
C        N IS THE ORDER OF THE MATRIX.                                  TRE00210
C                                                                       TRE00220
C        A CONTAINS THE REAL SYMMETRIC INPUT MATRIX.  ONLY THE          TRE00230
C          LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED.               TRE00240
C                                                                       TRE00250
C     ON OUTPUT                                                         TRE00260
C                                                                       TRE00270
C        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL TRANS-             TRE00280
C          FORMATIONS USED IN THE REDUCTION IN ITS STRICT LOWER         TRE00290
C          TRIANGLE.  THE FULL UPPER TRIANGLE OF A IS UNALTERED.        TRE00300
C                                                                       TRE00310
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX.    TRE00320
C                                                                       TRE00330
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL         TRE00340
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO.      TRE00350
C                                                                       TRE00360
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.    TRE00370
C          E2 MAY COINCIDE WITH E IF THE SQUARES ARE NOT NEEDED.        TRE00380
C                                                                       TRE00390
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,    TRE00400
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY TRE00410
C                                                                       TRE00420
C     THIS VERSION DATED AUGUST 1983.                                   TRE00430
C                                                                       TRE00440
C     ------------------------------------------------------------------TRE00450
C                                                                       TRE00460
      DO 100 I = 1, N                                                   TRE00470
         D(I) = A(N,I)                                                  TRE00480
         A(N,I) = A(I,I)                                                TRE00490
  100 CONTINUE                                                          TRE00500
C     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........               TRE00510
      DO 300 II = 1, N                                                  TRE00520
         I = N + 1 - II                                                 TRE00530
         L = I - 1                                                      TRE00540
         H = 0.0D0                                                      TRE00550
         SCALE = 0.0D0                                                  TRE00560
         IF (L .LT. 1) GO TO 130                                        TRE00570
C     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........       TRE00580
         DO 120 K = 1, L                                                TRE00590
  120    SCALE = SCALE + DABS(D(K))                                     TRE00600
C                                                                       TRE00610
         IF (SCALE .NE. 0.0D0) GO TO 140                                TRE00620
C                                                                       TRE00630
         DO 125 J = 1, L                                                TRE00640
            D(J) = A(L,J)                                               TRE00650
            A(L,J) = A(I,J)                                             TRE00660
            A(I,J) = 0.0D0                                              TRE00670
  125    CONTINUE                                                       TRE00680
C                                                                       TRE00690
  130    E(I) = 0.0D0                                                   TRE00700
         E2(I) = 0.0D0                                                  TRE00710
         GO TO 300                                                      TRE00720
C                                                                       TRE00730
  140    DO 150 K = 1, L                                                TRE00740
            D(K) = D(K) / SCALE                                         TRE00750
            H = H + D(K) * D(K)                                         TRE00760
  150    CONTINUE                                                       TRE00770
C                                                                       TRE00780
         E2(I) = SCALE * SCALE * H                                      TRE00790
         F = D(L)                                                       TRE00800
         G = -DSIGN(DSQRT(H),F)                                         TRE00810
         E(I) = SCALE * G                                               TRE00820
         H = H - F * G                                                  TRE00830
         D(L) = F - G                                                   TRE00840
         IF (L .EQ. 1) GO TO 285                                        TRE00850
C     .......... FORM A*U ..........                                    TRE00860
         DO 170 J = 1, L                                                TRE00870
  170    E(J) = 0.0D0                                                   TRE00880
C                                                                       TRE00890
         DO 240 J = 1, L                                                TRE00900
            F = D(J)                                                    TRE00910
            G = E(J) + A(J,J) * F                                       TRE00920
            JP1 = J + 1                                                 TRE00930
            IF (L .LT. JP1) GO TO 220                                   TRE00940
C                                                                       TRE00950
            DO 200 K = JP1, L                                           TRE00960
               G = G + A(K,J) * D(K)                                    TRE00970
               E(K) = E(K) + A(K,J) * F                                 TRE00980
  200       CONTINUE                                                    TRE00990
C                                                                       TRE01000
  220       E(J) = G                                                    TRE01010
  240    CONTINUE                                                       TRE01020
C     .......... FORM P ..........                                      TRE01030
         F = 0.0D0                                                      TRE01040
C                                                                       TRE01050
         DO 245 J = 1, L                                                TRE01060
            E(J) = E(J) / H                                             TRE01070
            F = F + E(J) * D(J)                                         TRE01080
  245    CONTINUE                                                       TRE01090
C                                                                       TRE01100
         H = F / (H + H)                                                TRE01110
C     .......... FORM Q ..........                                      TRE01120
         DO 250 J = 1, L                                                TRE01130
  250    E(J) = E(J) - H * D(J)                                         TRE01140
C     .......... FORM REDUCED A ..........                              TRE01150
         DO 280 J = 1, L                                                TRE01160
            F = D(J)                                                    TRE01170
            G = E(J)                                                    TRE01180
C                                                                       TRE01190
            DO 260 K = J, L                                             TRE01200
  260       A(K,J) = A(K,J) - F * E(K) - G * D(K)                       TRE01210
C                                                                       TRE01220
  280    CONTINUE                                                       TRE01230
C                                                                       TRE01240
  285    DO 290 J = 1, L                                                TRE01250
            F = D(J)                                                    TRE01260
            D(J) = A(L,J)                                               TRE01270
            A(L,J) = A(I,J)                                             TRE01280
            A(I,J) = F * SCALE                                          TRE01290
  290    CONTINUE                                                       TRE01300
C                                                                       TRE01310
  300 CONTINUE                                                          TRE01320
C                                                                       TRE01330
      RETURN                                                            TRE01340
      END                                                               TRE01350
      SUBROUTINE TRED2(NM,N,A,D,E,Z)                                    TRE00010
C                                                                       TRE00020
      INTEGER I,J,K,L,N,II,NM,JP1                                       TRE00030
      DOUBLE PRECISION A(NM,N),D(N),E(N),Z(NM,N)                        TRE00040
      DOUBLE PRECISION F,G,H,HH,SCALE                                   TRE00050
C                                                                       TRE00060
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED2,    TRE00070
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.   TRE00080
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).   TRE00090
C                                                                       TRE00100
C     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX TO A              TRE00110
C     SYMMETRIC TRIDIAGONAL MATRIX USING AND ACCUMULATING               TRE00120
C     ORTHOGONAL SIMILARITY TRANSFORMATIONS.                            TRE00130
C                                                                       TRE00140
C     ON INPUT                                                          TRE00150
C                                                                       TRE00160
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL         TRE00170
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM          TRE00180
C          DIMENSION STATEMENT.                                         TRE00190
C                                                                       TRE00200
C        N IS THE ORDER OF THE MATRIX.                                  TRE00210
C                                                                       TRE00220
C        A CONTAINS THE REAL SYMMETRIC INPUT MATRIX.  ONLY THE          TRE00230
C          LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED.               TRE00240
C                                                                       TRE00250
C     ON OUTPUT                                                         TRE00260
C                                                                       TRE00270
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX.    TRE00280
C                                                                       TRE00290
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL         TRE00300
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO.      TRE00310
C                                                                       TRE00320
C        Z CONTAINS THE ORTHOGONAL TRANSFORMATION MATRIX                TRE00330
C          PRODUCED IN THE REDUCTION.                                   TRE00340
C                                                                       TRE00350
C        A AND Z MAY COINCIDE.  IF DISTINCT, A IS UNALTERED.            TRE00360
C                                                                       TRE00370
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,    TRE00380
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY TRE00390
C                                                                       TRE00400
C     THIS VERSION DATED AUGUST 1983.                                   TRE00410
C                                                                       TRE00420
C     ------------------------------------------------------------------TRE00430
C                                                                       TRE00440
      DO 100 I = 1, N                                                   TRE00450
C                                                                       TRE00460
         DO 80 J = I, N                                                 TRE00470
   80    Z(J,I) = A(J,I)                                                TRE00480
C                                                                       TRE00490
         D(I) = A(N,I)                                                  TRE00500
  100 CONTINUE                                                          TRE00510
C                                                                       TRE00520
      IF (N .EQ. 1) GO TO 510                                           TRE00530
C     .......... FOR I=N STEP -1 UNTIL 2 DO -- ..........               TRE00540
      DO 300 II = 2, N                                                  TRE00550
         I = N + 2 - II                                                 TRE00560
         L = I - 1                                                      TRE00570
         H = 0.0D0                                                      TRE00580
         SCALE = 0.0D0                                                  TRE00590
         IF (L .LT. 2) GO TO 130                                        TRE00600
C     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........       TRE00610
         DO 120 K = 1, L                                                TRE00620
  120    SCALE = SCALE + DABS(D(K))                                     TRE00630
C                                                                       TRE00640
         IF (SCALE .NE. 0.0D0) GO TO 140                                TRE00650
  130    E(I) = D(L)                                                    TRE00660
C                                                                       TRE00670
         DO 135 J = 1, L                                                TRE00680
            D(J) = Z(L,J)                                               TRE00690
            Z(I,J) = 0.0D0                                              TRE00700
            Z(J,I) = 0.0D0                                              TRE00710
  135    CONTINUE                                                       TRE00720
C                                                                       TRE00730
         GO TO 290                                                      TRE00740
C                                                                       TRE00750
  140    DO 150 K = 1, L                                                TRE00760
            D(K) = D(K) / SCALE                                         TRE00770
            H = H + D(K) * D(K)                                         TRE00780
  150    CONTINUE                                                       TRE00790
C                                                                       TRE00800
         F = D(L)                                                       TRE00810
         G = -DSIGN(DSQRT(H),F)                                         TRE00820
         E(I) = SCALE * G                                               TRE00830
         H = H - F * G                                                  TRE00840
         D(L) = F - G                                                   TRE00850
C     .......... FORM A*U ..........                                    TRE00860
         DO 170 J = 1, L                                                TRE00870
  170    E(J) = 0.0D0                                                   TRE00880
C                                                                       TRE00890
         DO 240 J = 1, L                                                TRE00900
            F = D(J)                                                    TRE00910
            Z(J,I) = F                                                  TRE00920
            G = E(J) + Z(J,J) * F                                       TRE00930
            JP1 = J + 1                                                 TRE00940
            IF (L .LT. JP1) GO TO 220                                   TRE00950
C                                                                       TRE00960
            DO 200 K = JP1, L                                           TRE00970
               G = G + Z(K,J) * D(K)                                    TRE00980
               E(K) = E(K) + Z(K,J) * F                                 TRE00990
  200       CONTINUE                                                    TRE01000
C                                                                       TRE01010
  220       E(J) = G                                                    TRE01020
  240    CONTINUE                                                       TRE01030
C     .......... FORM P ..........                                      TRE01040
         F = 0.0D0                                                      TRE01050
C                                                                       TRE01060
         DO 245 J = 1, L                                                TRE01070
            E(J) = E(J) / H                                             TRE01080
            F = F + E(J) * D(J)                                         TRE01090
  245    CONTINUE                                                       TRE01100
C                                                                       TRE01110
         HH = F / (H + H)                                               TRE01120
C     .......... FORM Q ..........                                      TRE01130
         DO 250 J = 1, L                                                TRE01140
  250    E(J) = E(J) - HH * D(J)                                        TRE01150
C     .......... FORM REDUCED A ..........                              TRE01160
         DO 280 J = 1, L                                                TRE01170
            F = D(J)                                                    TRE01180
            G = E(J)                                                    TRE01190
C                                                                       TRE01200
            DO 260 K = J, L                                             TRE01210
  260       Z(K,J) = Z(K,J) - F * E(K) - G * D(K)                       TRE01220
C                                                                       TRE01230
            D(J) = Z(L,J)                                               TRE01240
            Z(I,J) = 0.0D0                                              TRE01250
  280    CONTINUE                                                       TRE01260
C                                                                       TRE01270
  290    D(I) = H                                                       TRE01280
  300 CONTINUE                                                          TRE01290
C     .......... ACCUMULATION OF TRANSFORMATION MATRICES ..........     TRE01300
      DO 500 I = 2, N                                                   TRE01310
         L = I - 1                                                      TRE01320
         Z(N,L) = Z(L,L)                                                TRE01330
         Z(L,L) = 1.0D0                                                 TRE01340
         H = D(I)                                                       TRE01350
         IF (H .EQ. 0.0D0) GO TO 380                                    TRE01360
C                                                                       TRE01370
         DO 330 K = 1, L                                                TRE01380
  330    D(K) = Z(K,I) / H                                              TRE01390
C                                                                       TRE01400
         DO 360 J = 1, L                                                TRE01410
            G = 0.0D0                                                   TRE01420
C                                                                       TRE01430
            DO 340 K = 1, L                                             TRE01440
  340       G = G + Z(K,I) * Z(K,J)                                     TRE01450
C                                                                       TRE01460
            DO 360 K = 1, L                                             TRE01470
               Z(K,J) = Z(K,J) - G * D(K)                               TRE01480
  360    CONTINUE                                                       TRE01490
C                                                                       TRE01500
  380    DO 400 K = 1, L                                                TRE01510
  400    Z(K,I) = 0.0D0                                                 TRE01520
C                                                                       TRE01530
  500 CONTINUE                                                          TRE01540
C                                                                       TRE01550
  510 DO 520 I = 1, N                                                   TRE01560
         D(I) = Z(N,I)                                                  TRE01570
         Z(N,I) = 0.0D0                                                 TRE01580
  520 CONTINUE                                                          TRE01590
C                                                                       TRE01600
      Z(N,N) = 1.0D0                                                    TRE01610
      E(1) = 0.0D0                                                      TRE01620
      RETURN                                                            TRE01630
      END                                                               TRE01640

c-------------------------------------------------------------------------
      SUBROUTINE GAUSSQ(KIND, N, ALPHA, BETA, KPTS, ENDPTS, B, T, W)    
      IMPLICIT REAL*8 (A-H,O-Z) 
C                                                                       
C           THIS SET OF ROUTINES COMPUTES THE NODES X(I) AND WEIGHTS    
C        C(I) FOR GAUSSIAN-TYPE QUADRATURE RULES WITH PRE-ASSIGNED      
C        NODES.  THESE ARE USED WHEN ONE WISHES TO APPROXIMATE          
C                                                                       
C                 INTEGRAL (FROM A TO B)  F(X) W(X) DX                  
C                                                                       
C                              N                                         
C        BY                   SUM C  F(X )                               
C                             I=1  I    I                                
C                                                                        
C        HERE W(X) IS ONE OF SIX POSSIBLE NON-NEGATIVE WEIGHT            
C        FUNCTIONS (LISTED BELOW), AND F(X) IS THE                       
C        FUNCTION TO BE INTEGRATED.  GAUSSIAN QUADRATURE IS PARTICULARLY 
C        USEFUL ON INFINITE INTERVALS (WITH APPROPRIATE WEIGHT           
C        FUNCTIONS), SINCE THEN OTHER TECHNIQUES OFTEN FAIL.             
C                                                                        
C           ASSOCIATED WITH EACH WEIGHT FUNCTION W(X) IS A SET OF        
C        ORTHOGONAL POLYNOMIALS.  THE NODES X(I) ARE JUST THE ZEROES     
C        OF THE PROPER N-TH DEGREE POLYNOMIAL.                           
C                                                                        
C     INPUT PARAMETERS                                                   
C                                                                        
C        KIND     AN INTEGER BETWEEN 1 AND 6 GIVING THE TYPE OF          
C                 QUADRATURE RULE                                        
C                                                                        
C        KIND = 1=  LEGENDRE QUADRATURE, W(X) = 1 ON (-1, 1)             
C        KIND = 2=  CHEBYSHEV QUADRATURE OF THE FIRST KIND               
C                   W(X) = 1/SQRT(1 - X*X) ON (-1, +1)                   
C        KIND = 3=  CHEBYSHEV QUADRATURE OF THE SECOND KIND              
C                   W(X) = SQRT(1 - X*X) ON (-1, 1)                      
C        KIND = 4=  HERMITE QUADRATURE, W(X) = EXP(-X*X) ON              
C                   (-INFINITY, +INFINITY)                               
C        KIND = 5=  JACOBI QUADRATURE, W(X) = (1-X)**ALPHA * (1+X)**     
C                   BETA ON (-1, 1), ALPHA, BETA .GT. -1.                
C                   NOTE= KIND=2 AND 3 ARE A SPECIAL CASE OF THIS.       
C        KIND = 6=  GENERALIZED LAGUERRE QUADRATURE, W(X) = EXP(-X)*     
C                   X**ALPHA ON (0, +INFINITY), ALPHA .GT. -1            
C                                                                        
C        N        THE NUMBER OF POINTS USED FOR THE QUADRATURE RULE      
C        ALPHA    REAL PARAMETER USED ONLY FOR GAUSS-JACOBI AND GAUSS-   
C                 LAGUERRE QUADRATURE (OTHERWISE USE 0.).                
C        BETA     REAL PARAMETER USED ONLY FOR GAUSS-JACOBI QUADRATURE-- 
C                 (OTHERWISE USE 0.).                                    
C        KPTS     (INTEGER) NORMALLY 0, UNLESS THE LEFT OR RIGHT END-    
C                 POINT (OR BOTH) OF THE INTERVAL IS REQUIRED TO BE A    
C                 NODE (THIS IS CALLED GAUSS-RADAU OR GAUSS-LOBATTO      
C                 QUADRATURE).  THEN KPTS IS THE NUMBER OF FIXED         
C                 ENDPOINTS (1 OR 2).                                    
C        ENDPTS   REAL ARRAY OF LENGTH 2.  CONTAINS THE VALUES OF        
C                 ANY FIXED ENDPOINTS, IF KPTS = 1 OR 2.                 
C        B        REAL SCRATCH ARRAY OF LENGTH N                         
C                                                                        
C     OUTPUT PARAMETERS (BOTH ARRAYS OF LENGTH N)                        
C                                                                        
C        T        WILL CONTAIN THE DESIRED NODES X(1),,,X(N)             
C        W        WILL CONTAIN THE DESIRED WEIGHTS C(1),,,C(N)           
C                                                                        
C     SUBROUTINES REQUIRED                                               
C                                                                        
C        GBSLVE, CLASS, AND GBTQL2 ARE PROVIDED. UNDERFLOW MAY SOMETIMES 
C        OCCUR, BUT IT IS HARMLESS IF THE UNDERFLOW INTERRUPTS ARE       
C        TURNED OFF AS THEY ARE ON THIS MACHINE.                         
C                                                                        
C     ACCURACY                                                           
C                                                                        
C        THE ROUTINE WAS TESTED UP TO N = 512 FOR LEGENDRE QUADRATURE,   
C        UP TO N = 136 FOR HERMITE, UP TO N = 68 FOR LAGUERRE, AND UP    
C        TO N = 10 OR 20 IN OTHER CASES.  IN ALL BUT TWO INSTANCES,      
C        COMPARISON WITH TABLES IN REF. 3 SHOWED 12 OR MORE SIGNIFICANT  
C        DIGITS OF ACCURACY.  THE TWO EXCEPTIONS WERE THE WEIGHTS FOR    
C        HERMITE AND LAGUERRE QUADRATURE, WHERE UNDERFLOW CAUSED SOME    
C        VERY SMALL WEIGHTS TO BE SET TO ZERO.  THIS IS, OF COURSE,      
C        COMPLETELY HARMLESS.                                            
C                                                                        
C     METHOD                                                             
C                                                                        
C           THE COEFFICIENTS OF THE THREE-TERM RECURRENCE RELATION       
C        FOR THE CORRESPONDING SET OF ORTHOGONAL POLYNOMIALS ARE         
C        USED TO FORM A SYMMETRIC TRIDIAGONAL MATRIX, WHOSE              
C        EIGENVALUES (DETERMINED BY THE IMPLICIT QL-METHOD WITH          
C        SHIFTS) ARE JUST THE DESIRED NODES.  THE FIRST COMPONENTS OF    
C        THE ORTHONORMALIZED EIGENVECTORS, WHEN PROPERLY SCALED,         
C        YIELD THE WEIGHTS.  THIS TECHNIQUE IS MUCH FASTER THAN USING A  
C        ROOT-FINDER TO LOCATE THE ZEROES OF THE ORTHOGONAL POLYNOMIAL.  
C        FOR FURTHER DETAILS, SEE REF. 1.  REF. 2 CONTAINS DETAILS OF    
C        GAUSS-RADAU AND GAUSS-LOBATTO QUADRATURE ONLY.                  
C                                                                        
C     REFERENCES                                                         
C                                                                        
C        1.  GOLUB, G. H., AND WELSCH, J. H.,  CALCULATION OF GAUSSIAN   
C            QUADRATURE RULES,  MATHEMATICS OF COMPUTATION 23 (APRIL,    
C            1969), PP. 221-230.                                         
C        2.  GOLUB, G. H.,  SOME MODIFIED MATRIX EIGENVALUE PROBLEMS,    
C            SIAM REVIEW 15 (APRIL, 1973), PP. 318-334 (SECTION 7).      
C        3.  STROUD AND SECREST, GAUSSIAN QUADRATURE FORMULAS, PRENTICE- 
C            HALL, ENGLEWOOD CLIFFS, N.J., 1966.                         
C                                                                       
C     ..................................................................
C                                                                       
      REAL*8  MUZERO                                                    
      DIMENSION  B(N),T(N),W(N),ENDPTS(2)                               
C                                                                       
      CALL CLASS (KIND, N, ALPHA, BETA, B, T, MUZERO)                   
C                                                                       
C           THE MATRIX OF COEFFICIENTS IS ASSUMED TO BE SYMMETRIC.      
C           THE ARRAY T CONTAINS THE DIAGONAL ELEMENTS, THE ARRAY       
C           B THE OFF-DIAGONAL ELEMENTS.                                
C           MAKE APPROPRIATE CHANGES IN THE LOWER RIGHT 2 BY 2          
C           SUBMATRIX.                                                  
C                                                                       
      IF (KPTS.EQ.0)  GO TO 100                                         
      IF (KPTS.EQ.2)  GO TO  50                                         
C                                                                       
C           IF KPTS=1, ONLY T(N) MUST BE CHANGED                        
C                                                                       
      T(N) =GBSLVE(ENDPTS(1), N, T, B)*B(N-1)**2 + ENDPTS(1)            
      GO TO 100                                                         
C                                                                       
C           IF KPTS=2, T(N) AND B(N-1) MUST BE RECOMPUTED               
C                                                                       
   50 GAM =GBSLVE(ENDPTS(1), N, T, B)                                   
      T1 = ((ENDPTS(1) - ENDPTS(2))/(GBSLVE(ENDPTS(2), N, T, B) - GAM)) 
      B(N-1) =  SQRT(T1)                                                
      T(N) = ENDPTS(1) + GAM*T1                                         
C                                                                       
C           NOTE THAT THE INDICES OF THE ELEMENTS OF B RUN FROM 1 TO N-1
C           AND THUS THE VALUE OF B(N) IS ARBITRARY.                    
C           NOW COMPUTE THE EIGENVALUES OF THE SYMMETRIC TRIDIAGONAL    
C           MATRIX, WHICH HAS BEEN MODIFIED AS NECESSARY.               
C           THE METHOD USED IS A QL-TYPE METHOD WITH ORIGIN SHIFTING    
C                                                                       
  100 W(1) = 1.0D0                                                      
      DO 105 I = 2, N                                                   
  105    W(I) = 0.0D0                                                   
C                                                                       
      CALL GBTQL2 (N, T, B, W, IERR)                                    
      DO 110 I = 1, N                                                   
  110    W(I) = MUZERO * W(I) * W(I)                                    
C                                                                       
      RETURN                                                            
      END                                                               
C                                                                       
C                                                                       
C                                                                       
      DOUBLE PRECISION FUNCTION GBSLVE(SHIFT, N, A, B)                                   
      IMPLICIT REAL*8 (A-H,O-Z) 
C                                                                       
C       THIS PROCEDURE PERFORMS ELIMINATION TO SOLVE FOR THE            
C       N-TH COMPONENT OF THE SOLUTION DELTA TO THE EQUATION            
C                                                                       
C             (JN - SHIFT*IDENTITY) * DELTA  = EN,                      
C                                                                       
C       WHERE EN IS THE VECTOR OF ALL ZEROES EXCEPT FOR 1 IN            
C       THE N-TH POSITION.                                              
C                                                                       
C       THE MATRIX JN IS SYMMETRIC TRIDIAGONAL, WITH DIAGONAL           
C       ELEMENTS A(I), OFF-DIAGONAL ELEMENTS B(I).  THIS EQUATION       
C       MUST BE SOLVED TO OBTAIN THE APPROPRIATE CHANGES IN THE LOWER   
C       2 BY 2 SUBMATRIX OF COEFFICIENTS FOR ORTHOGONAL POLYNOMIALS.    
C                                                                       
C                                                                       
      DIMENSION  A(N),B(N)                                              
C                                                                       
      ALPHA = A(1) - SHIFT                                              
      NM1 = N - 1                                                       
      DO 10 I = 2, NM1                                                  
   10    ALPHA = A(I) - SHIFT - B(I-1)**2/ALPHA                         
      GBSLVE = 1.0D0  /ALPHA                                            
      RETURN                                                            
      END                                                               
C                                                                       
C                                                                       
C                                                                       
      SUBROUTINE CLASS(KIND, N, ALPHA, BETA, B, A, MUZERO)              
      IMPLICIT REAL*8 (A-H,O-Z) 
C                                                                       
C           THIS PROCEDURE SUPPLIES THE COEFFICIENTS A(J), B(J) OF THE  
C        RECURRENCE RELATION                                            
C                                                                       
C             B P (X) = (X - A ) P   (X) - B   P   (X)                  
C              J J            J   J-1       J-1 J-2                     
C                                                                       
C        FOR THE VARIOUS CLASSICAL (NORMALIZED) ORTHOGONAL POLYNOMIALS, 
C        AND THE ZERO-TH MOMENT                                         
C                                                                       
C             MUZERO = INTEGRAL W(X) DX                                 
C                                                                       
C        OF THE GIVEN POLYNOMIAL   WEIGHT FUNCTION W(X).  SINCE THE     
C        POLYNOMIALS ARE ORTHONORMALIZED, THE TRIDIAGONAL MATRIX IS     
C        GUARANTEED TO BE SYMMETRIC.                                    
C                                                                       
C           THE INPUT PARAMETER ALPHA IS USED ONLY FOR LAGUERRE AND     
C        JACOBI POLYNOMIALS, AND THE PARAMETER BETA IS USED ONLY FOR    
C        JACOBI POLYNOMIALS.  THE LAGUERRE AND JACOBI POLYNOMIALS       
C        REQUIRE THE GAMMA FUNCTION.                                    
C                                                                       
C     ..................................................................
C                                                                       
      DIMENSION  A(N),B(N)                                              
      REAL*8  MUZERO                                                    
      DATA PI / 3.141592653589793D0  /                                  
C                                                                       
      NM1 = N - 1                                                       
      GO TO (10, 20, 30, 40, 50, 60), KIND                              
C                                                                       
C              KIND = 1=  LEGENDRE POLYNOMIALS P(X)                     
C              ON (-1, +1), W(X) = 1.                                   
C                                                                       
   10 MUZERO = 2.0D0                                                    
      DO 11 I = 1, NM1                                                  
         A(I) = 0.0D0                                                   
         ABI = I                                                        
   11    B(I) = ABI/ SQRT(4*ABI*ABI - 1.0D0  )                          
      A(N) = 0.0D0                                                      
      RETURN                                                            
C                                                                       
C              KIND = 2=  CHEBYSHEV POLYNOMIALS OF THE FIRST KIND T(X)  
C              ON (-1, +1), W(X) = 1 / SQRT(1 - X*X)                    
C                                                                       
   20 MUZERO = PI                                                       
      DO 21 I = 1, NM1                                                  
         A(I) = 0.0D0                                                   
   21    B(I) = 0.5D0                                                   
      B(1) =  SQRT(0.5D0  )                                             
      A(N) = 0.0D0                                                      
      RETURN                                                            
C                                                                       
C              KIND = 3=  CHEBYSHEV POLYNOMIALS OF THE SECOND KIND U(X) 
C              ON (-1, +1), W(X) = SQRT(1 - X*X)                        
C                                                                       
   30 MUZERO = PI/2.0D0                                                 
      DO 31 I = 1, NM1                                                  
         A(I) = 0.0D0                                                   
   31    B(I) = 0.5D0                                                   
      A(N) = 0.0D0                                                      
      RETURN                                                            
C                                                                       
C              KIND = 4=  HERMITE POLYNOMIALS H(X) ON (-INFINITY,       
C              +INFINITY), W(X) = EXP(-X**2)                            
C                                                                       
   40 MUZERO =  SQRT(PI)                                                
      DO 41 I = 1, NM1                                                  
         A(I) = 0.0D0                                                   
   41    B(I) =  SQRT(I/2.0D0  )                                        
      A(N) = 0.0D0                                                      
      RETURN                                                            
C                                                                       
C              KIND = 5=  JACOBI POLYNOMIALS P(ALPHA, BETA)(X) ON       
C              (-1, +1), W(X) = (1-X)**ALPHA + (1+X)**BETA, ALPHA AND   
C              BETA GREATER THAN -1                                     
C                                                                       
   50 AB = ALPHA + BETA                                                 
      ABI = 2.0D0   + AB                                                
      MUZERO = 2.0D0   ** (AB + 1.0D0  ) * DGAMMA(ALPHA + 1.0D0  ) * DGA
     VMMA(
     X BETA + 1.0D0  ) / DGAMMA(ABI)                                    
      A(1) = (BETA - ALPHA)/ABI                                         
      B(1) =  SQRT(4.0D0  *(1.0D0   + ALPHA)*(1.0D0   + BETA)/((ABI + 1.
     V0D0  )* 
     1  ABI*ABI))                                                       
      A2B2 = BETA*BETA - ALPHA*ALPHA                                    
      DO 51 I = 2, NM1                                                  
         ABI = 2.0D0  *I + AB                                           
         A(I) = A2B2/((ABI - 2.0D0  )*ABI)                              
   51    B(I) =  SQRT (4.0D0  *I*(I + ALPHA)*(I + BETA)*(I + AB)/       
     1   ((ABI*ABI - 1)*ABI*ABI))                                       
      ABI = 2.0D0  *N + AB                                              
      A(N) = A2B2/((ABI - 2.0D0  )*ABI)                                 
      RETURN                                                            
C                                                                       
C              KIND = 6=  LAGUERRE POLYNOMIALS L(ALPHA)(X) ON           
C              (0, +INFINITY), W(X) = EXP(-X) * X**ALPHA, ALPHA GREATER 
C              THAN -1.                                                 
C                                                                       
   60 MUZERO = DGAMMA(ALPHA + 1.0D0  )                                  
      DO 61 I = 1, NM1                                                  
         A(I) = 2.0D0  *I - 1.0D0   + ALPHA                             
   61    B(I) =  SQRT(I*(I + ALPHA))                                    
      A(N) = 2.0D0  *N - 1 + ALPHA                                      
      RETURN                                                            
      END                                                               
C     ------------------------------------------------------------------
C                                                                       
      SUBROUTINE GBTQL2(N, D, E, Z, IERR)                               
      IMPLICIT REAL*8 (A-H,O-Z) 
C                                                                       
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE IMTQL2,   
C     NUM. MATH. 12, 377-383(1968) BY MARTIN AND WILKINSON,             
C     AS MODIFIED IN NUM. MATH. 15, 450(1970) BY DUBRULLE.              
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 241-248(1971).   
C                                                                       
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND FIRST COMPONENTS OF THE 
C     EIGENVECTORS OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE IMPLICIT QL 
C     METHOD, AND IS ADAPTED FROM THE EISPAK ROUTINE IMTQL2             
C                                                                       
C     ON INPUT=                                                         
C                                                                       
C        N IS THE ORDER OF THE MATRIX;                                  
C                                                                       
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX;          
C                                                                       
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX        
C          IN ITS FIRST N-1 POSITIONS.  E(N) IS ARBITRARY;              
C                                                                       
C        Z CONTAINS THE FIRST ROW OF THE IDENTITY MATRIX.               
C                                                                       
C      ON OUTPUT=                                                       
C                                                                       
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN          
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT          
C          UNORDERED FOR INDICES 1, 2, ..., IERR-1;                     
C                                                                       
C        E HAS BEEN DESTROYED;                                          
C                                                                       
C        Z CONTAINS THE FIRST COMPONENTS OF THE ORTHONORMAL EIGENVECTORS
C          OF THE SYMMETRIC TRIDIAGONAL MATRIX.  IF AN ERROR EXIT IS    
C          MADE, Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED 
C          EIGENVALUES;                                                 
C                                                                       
C        IERR IS SET TO                                                 
C                                                                       
C        IERR IS SET TO                                                 
C          ZERO       FOR NORMAL RETURN,                                
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN               
C                     DETERMINED AFTER 30 ITERATIONS.                   
C                                                                       
C     ------------------------------------------------------------------
C                                                                       
      INTEGER I, J, K, L, M, N, II, MML, IERR                           
      DIMENSION  D(N),E(N),Z(N)                                         
      REAL*8  MACHEP                                                    
C                                                                       
C     ========== MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING     
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.   
C                MACHEP = 16.0D0**(-13) FOR LONG FORM ARITHMETIC        
C                ON S360 ==========                                     
       MACHEP=1.0D-14                                                   
C                                                                       
      IERR = 0                                                          
      IF (N .EQ. 1) GO TO 1001                                          
C                                                                       
      E(N) = 0.0D0                                                      
      DO 240 L = 1, N                                                   
         J = 0                                                          
C     ========== LOOK FOR SMALL SUB-DIAGONAL ELEMENT ==========         
  105    DO 110 M = L, N                                                
            IF (M .EQ. N) GO TO 120                                     
            IF ( ABS(E(M)) .LE. MACHEP * ( ABS(D(M)) +  ABS(D(M+1))))   
     X         GO TO 120                                                
  110    CONTINUE                                                       
C                                                                       
  120    P = D(L)                                                       
         IF (M .EQ. L) GO TO 240                                        
         IF (J .EQ. 30) GO TO 1000                                      
         J = J + 1                                                      
C     ========== FORM SHIFT ==========                                  
         G = (D(L+1) - P) / (2.0D0   * E(L))                            
         R =  SQRT(G*G+1.0D0  )                                         
         G = D(M) - P + E(L) / (G +  SIGN(R, G))                        
         S = 1.0D0                                                      
         C = 1.0D0                                                      
         P = 0.0D0                                                      
         MML = M - L                                                    
C     ========== FOR I=M-1 STEP -1 UNTIL L DO -- ==========             
         DO 200 II = 1, MML                                             
            I = M - II                                                  
            F = S * E(I)                                                
            B = C * E(I)                                                
            IF ( ABS(F) .LT.  ABS(G)) GO TO 150                         
            C = G / F                                                   
            R =  SQRT(C*C+1.0D0  )                                      
            E(I+1) = F * R                                              
            S = 1.0D0   / R                                             
            C = C * S                                                   
            GO TO 160                                                   
  150       S = F / G                                                   
            R =  SQRT(S*S+1.0D0  )                                      
            E(I+1) = G * R                                              
            C = 1.0D0   / R                                             
            S = S * C                                                   
  160       G = D(I+1) - P                                              
            R = (D(I) - G) * S + 2.0D0   * C * B                        
            P = S * R                                                   
            D(I+1) = G + P                                              
            G = C * R - B                                               
C     ========== FORM FIRST COMPONENT OF VECTOR ==========              
            F = Z(I+1)                                                  
            Z(I+1) = S * Z(I) + C * F                                   
            Z(I) = C * Z(I) - S * F                                     
C                                                                       
  200    CONTINUE                                                       
C                                                                       
         D(L) = D(L) - P                                                
         E(L) = G                                                       
         E(M) = 0.0D0                                                   
         GO TO 105                                                      
  240 CONTINUE                                                          
C     ========== ORDER EIGENVALUES AND EIGENVECTORS ==========          
      DO 300 II = 2, N                                                  
         I = II - 1                                                     
         K = I                                                          
         P = D(I)                                                       
C                                                                       
         DO 260 J = II, N                                               
            IF (D(J) .GE. P) GO TO 260                                  
            K = J                                                       
            P = D(J)                                                    
  260    CONTINUE                                                       
C                                                                       
         IF (K .EQ. I) GO TO 300                                        
         D(K) = D(I)                                                    
         D(I) = P                                                       
C                                                                       
         P = Z(I)                                                       
         Z(I) = Z(K)                                                    
         Z(K) = P                                                       
C                                                                       
  300 CONTINUE                                                          
C                                                                       
      GO TO 1001                                                        
C     ========== SET ERROR -- NO CONVERGENCE TO AN                      
C                EIGENVALUE AFTER 30 ITERATIONS ==========              
 1000 IERR = L                                                          
 1001 RETURN                                                            
C     ========== LAST CARD OF GBTQL2 ==========                         
      END                                                               
      DOUBLE PRECISION FUNCTION  DGAMMA(Z)                                               
      IMPLICIT REAL*8 (A-H,O-Z) 
C  THIS IS A PROCEDURE THAT EVALUATES GAMMA(Z) FOR                      
C     0 LT Z LE 3 TO 16 SIGNIFICANT FIGURES                             
C    IT IS BASED ON A CHEBYSHEV-TYPE POLYNOMIAL                         
C   APPROXIMATION GIVEN IN H. WERNER AND R. COLLINGE, MATH. COMPUT.     
C    15 (1961), PP. 195-97.                                             
C   APPROXIMATIONS TO THE GAMMA FUNCTION, ACCURATE UP TO 18 SIGNIFICANT 
C   DIGITS, MAY BE FOUND IN THE PAPER QUOTED ABOVE                      
C                                                                       
C                                                                       
C                                                                       
      DIMENSION  A(18)                                                  
C                                                                       
       A(1)=1.0D0                                                       
       A(2)=.4227843350984678D0                                         
       A(3)=.4118403304263672D0                                         
      A(4)=.0815769192502609D0                                          
      A(5)=.0742490106800904D0                                          
      A(6)=-.0002669810333484D0                                         
      A(7)=.0111540360240344D0                                          
      A(8)=-.0028525821446197D0                                         
      A(9)=.0021036287024598D0                                          
      A(10)=-.0009184843690991D0                                        
      A(11)=.0004874227944768D0                                         
      A(12)=-.0002347204018919D0                                        
      A(13)=.0001115339519666D0                                         
      A(14)=-.0000478747983834D0                                        
      A(15)=.0000175102727179D0                                         
      A(16)=-.0000049203750904D0                                        
      A(17)=.0000009199156407D0                                         
      A(18)=-.0000000839940496D0                                        
C                                                                       
C                                                                       
C                                                                       
      IF(Z.LE.1.0D0  ) GO TO 10                                         
      IF(Z.LE.2.0D0  ) GO TO 20                                         
      T=Z-2.0D0                                                         
      GO TO 30                                                          
10    T=Z                                                               
      GO TO 30                                                          
20    T=Z-1.0D0                                                         
30    P=A(18)                                                           
      DO 40 K1=1,17                                                     
      K=18-K1                                                           
      P=T*P+A(K)                                                        
40    CONTINUE                                                          
C                                                                       
      IF(Z.GT.2.0D0  ) GO TO 50                                         
      IF(Z.GT.1.0D0  ) GO TO 60                                         
      DGAMMA=P/(Z*(Z+1.0D0  ))                                          
      RETURN                                                            
60    DGAMMA=P/Z                                                        
      RETURN                                                            
50    DGAMMA=P                                                          
      RETURN                                                            
      END                                                               


c***********************************************************************
c
c                     function yjm0
c
c     this function returns the value of the spherical
c     harmonic y(j,m) when cos(theta)=x and phi=0.
c
c     input:
c       x      cos(theta) 
c       j,m    order of the spherical harmonic
c
c     output:
c       yjm0   value of Y(j,m|cos(theta),0)
c
c***********************************************************************
      function yjm0(j,m,x)
      implicit real*8 (a-h,o-z)
      pi = dacos(-1.d0)
c-----------------------------------------------------------------------
c
c      -- first calculate the associated legendre function
c
c-----------------------------------------------------------------------
      mabs = abs(m)
      pjm0 = plgndr(j,mabs,x)
c-----------------------------------------------------------------------
c
c     -- calculate the spherical harmonics
c
c-----------------------------------------------------------------------
c
c     -- positive m
c
      factor = (2.*j+1) * factrl(j-mabs) / (4.*pi*factrl(j+mabs) )
      sqfact = sqrt(factor)
      yjm0 = (-1.)**mabs * sqfact * pjm0
c
c     -- negative m
c
      if(m .lt. 0) then
         phase = (-1.)**m
         yjm0  = phase * yjm0
      end if
      yjm0=sqrt(2.d0*pi)*yjm0
      return
      end
c***********************************************************************
c
c                     function dyjm0
c
c     this function returns the value of the derivative of the 
c     spherical harmonic y(j,m) wrt theta when cos(theta)=x 
c     and phi=0.
c
c     input:
c       x      cos(theta) 
c       j,m    order of the spherical harmonic
c
c     output:
c       dyjm0  value of dY(j,m|cos(theta),0)/dtheta
c
c***********************************************************************
      function dyjm0(j,m,x)
      implicit real*8 (a-h,o-z)
      pi = dacos(-1.d0)
      if(j .eq. 0) then
         dyjm0 = 0.
         return
      end if
c-----------------------------------------------------------------------
c
c      -- first calculate the associated legendre functions
c
c-----------------------------------------------------------------------
      if(m .eq. 0) then
         pjmp = plgndr(j,1,x)
         factor = -factrl(j-1) / factrl(j+1)
         pjmm = factor * pjmp
      else if(m .gt. 0) then
         pjmp = plgndr(j,m+1,x)
         pjmm = plgndr(j,m-1,x)
      else if(m .lt. 0) then
         mabs = iabs(m+1)
         pjmp = plgndr(j,mabs,x)
         factor = (-1.)**mabs * factrl(j-mabs) / factrl(j+mabs)
         pjmp = factor * pjmp
         mabs = iabs(m-1)
         pjmm = plgndr(j,mabs,x)
         factor = (-1.)**mabs * factrl(j-mabs) / factrl(j+mabs)
         pjmm = factor * pjmm
      end if
      dpjm0 = .5 * ( pjmp - (j+m)*(j-m+1)*pjmm )
c-----------------------------------------------------------------------
c
c     -- calculate the derivative of the spherical harmonic
c
c-----------------------------------------------------------------------
      factor = (2.*j+1) * factrl(j-m) / (4.*pi*factrl(j+m) )
      sqfact = sqrt(factor)
      dyjm0 = -(-1.)**m * sqfact * dpjm0
      return
      end
c***********************************************************************
c
c                       function plgndr
c
c     This function was taken from the book "Numerical Recipes -
c     the Art of Scientific Computing " by W.H. Press, B.P. Flannery,
c     S.A. Teukolsky and W.T. Vetterling, published by Cambridge Univ.
c     Press, New York, 1986. See p.182.
c
c     Computes the associated Legendre polynomial P_{l,m}(x). Here,
c     m and l are integers satisfying (0 .le. m .le. l), while x 
c     lies in the range (-1 .le. x .le. 1).
c
c***********************************************************************
      function plgndr(l,m,x)
      implicit real*8(a-h,o-z)
      if(m.lt.0.or.m.gt.l.or.abs(x).gt.1.) then
         write(6,*) ' *** fatal error in function plgndr '
         write(6,*) ' *** bad arguements (l,m,x): ', l, m, x
      end if
      pmm=1.
      if(m.gt.0) then
        somx2=sqrt((1.-x)*(1.+x))
        fact=1.
        do 11 i=1,m
          pmm=-pmm*fact*somx2
          fact=fact+2.
11      continue
      endif
      if(l.eq.m) then
        plgndr=pmm
      else
        pmmp1=x*(2*m+1)*pmm
        if(l.eq.m+1) then
          plgndr=pmmp1
        else
          do 12 ll=m+2,l
            pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
            pmm=pmmp1
            pmmp1=pll
12        continue
          plgndr=pll
        endif
      endif
      return
      end
c***********************************************************************
c
c                       function factrl
c
c     This function was taken from the book "Numerical Recipes -
c     the Art of Scientific Computing " by W.H. Press, B.P. Flannery,
c     S.A. Teukolsky and W.T. Vetterling, published by Cambridge Univ.
c     Press, New York, 1986. See p.182.
c
c     Computes the value of N! as a floating point number.
c
c***********************************************************************
      function factrl(n)
      implicit real*8(a-h,o-z)
      dimension a(33)
      data ntop,a(1)/0,1./
      if (n.lt.0) then
        write(6,*) ' *** fatal error in routine factrl '
        write(6,*) ' *** negative factorial n: ', n
        stop
      else if (n.le.ntop) then
        factrl=a(n+1)
      else if (n.le.32) then
        do 11 j=ntop+1,n
          a(j+1)=j*a(j)
11      continue
        ntop=n
        factrl=a(n+1)
      else
	x=n+1.d0
        factrl=exp(gammln(x))
      endif
      return
      end
c***********************************************************************
c
c                       function gammln
c
c     This function was taken from the book "Numerical Recipes -
c     the Art of Scientific Computing " by W.H. Press, B.P. Flannery,
c     S.A. Teukolsky and W.T. Vetterling, published by Cambridge Univ.
c     Press, New York, 1986. See p.182.
c
c     Computes the value ln[gamma(xx)] for xx>0. Full accuracy is
c     obtained for xx>1. For 0<xx<1, the reflection formula (6.1.4)
c     can be used first.
c
c***********************************************************************
      function gammln(xx)
      implicit real*8(a-h,o-z)
      real*8 cof(6),stp,half,one,fpf,x,tmp,ser
      data cof,stp/76.18009173d0,-86.50532033d0,24.01409822d0,
     *    -1.231739516d0,.120858003d-2,-.536382d-5,2.50662827465d0/
      data half,one,fpf/0.5d0,1.0d0,5.5d0/
      x=xx-one
      tmp=x+fpf
      tmp=(x+half)*log(tmp)-tmp
      ser=one
      do 11 j=1,6
        x=x+one
        ser=ser+cof(j)/x
11    continue
      gammln=tmp+log(stp*ser)
      return
      end
c===========================================================
c	implicit real*8 (a-h,o-z)
c	dimension fac(0:100),pl(100)

c	fac(0)=1.D0
c	do i=1,100
c	  fac(i)=factrl(i)
c	end do
c
c	read*,l,m,x
c	call yjm0M(l,m,x,pl,fac)
c
c	m0=iabs(m)
c	do i=1,l-m0+1
c	  j=m0+i-1
c	  write(*,*) pl(i),yjm0(j,m,x),j,m
c	end do
c
c	stop
c	end
c***********************************************************************




c***********************************************************************
c
c                     subroutine yjm0M
c
c     this function returns the value of the spherical
c     harmonic y(l,m) when cos(theta)=x and phi=0 for
c     m <= l <= j
c
c     input:
c       x      cos(theta)
c       j,m    order of the spherical harmonic
c
c     output:
c       y(*)   value of Y(j,m|cos(theta),0)*sqrt(2 pi)
c
c***********************************************************************
c***********************************************************************
      subroutine yjm0M(l,m0,x,y,fac)
      implicit real*8(a-h,o-z)
      dimension y(1),fac(0:1)

      tmp=(-1)**m0
      m=m0
      if(m0.lt.0) then
        m=iabs(m0)
        tmp=1.d0
      end if

      if(m.gt.l.or.abs(x).gt.1.) then
         write(6,*) ' *** fatal error in function yjm0M '
         write(6,*) ' *** bad arguements (l,m,x): ', l, m0, x
      end if

      pmm=1.D0
      if(m.gt.0) then
        somx2=sqrt((1.-x)*(1.+x))
        fact=1.D0
        do 11 i=1,m
          pmm=-pmm*fact*somx2
          fact=fact+2.D0
11      continue
      endif

      y(1)=pmm
      y(2)=x*(2*m+1)*pmm

      N=m-1
      do 12 ll=m+2,l
        y(ll-N)=(x*(2*ll-1)*y(ll-N-1)-(ll+m-1)*y(ll-N-2))/(ll-m)
12    continue

      do 20 j=m,l
c      factor = sqrt( (2.*j+1) * fac(j-m) / (4.*pi*fac(j+m) ) )   
      factor = sqrt( (2.*j+1) * fac(j-m) / (2.*fac(j+m) ) )   ! * sqrt(2 pi)
      y(j-N)=y(j-N)*factor*tmp
20    continue

      return
      end
c***********************************************************************
	SUBROUTINE DJMM(J,MP,M,BETA,R,ID)
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION AF(100)

c when id.eq.1 , d matrix reduced to lengendre ploynomials(m'=0 and m=0)
c when id.eq.2 , d matrix reduced to associate ploynomials(m=0)
c when id.eq.3 , d matrix, all three cases are normalized for theta

	PI=DACOS(-1.D0)
	CALL NFACTO(50,AF,0)
	CALL DLEB(J,MP,M,BETA,R,AF)
	
	IF(ID.EQ.1) THEN
	  IF(MP.NE.0.OR.M.NE.0) STOP 'BAD INPUT OF M'
          R=DSQRT((2.D0*J+1)/2.D0)*R
	END IF

	IF(ID.EQ.2) THEN
	  IF(M.NE.0) STOP 'BAD INPUT OF M'
	
c there is a factor between the normal lengendre and reduced d ,
c it is  dsqrt((j+m)!/(j-m)!)

	  RR=R*DSQRT(AF(J+MP+1)/AF(J-MP+1))
	  FACTOR=DSQRT((2.D0*J+1)*AF(J-MP+1)/(4.D0*PI*AF(J+MP+1)))
	  R=(-1.D0)**MP*FACTOR*RR
	
	  R=DSQRT(2.D0*PI)*R
	END IF

	IF(ID.EQ.3) R=(2.D0*PI)*DSQRT((2.D0*J+1)/(8.D0*PI**2))*R

	RETURN
	END 	
c=================================================================
         SUBROUTINE DLEB(J,MP,M,BETA,R,AF)
         IMPLICIT REAL*8(A-H,O-Z)                                       CHE05720
C***********************************************************************
C        THIS SUBROUTINE CALCULATES THE SMALL D COEF                    CHE05740
C        EQ. 4.13 ON PAGE 52 OF M.E.ROSE                                CHE05750
C***********************************************************************
      DIMENSION AF(1)
         C=DCOS(BETA*0.5D0)
         S=-DSIN(BETA*0.5D0)
      N1=J-MP
      N2=J+M
      N3=MP-M
      NC=2*J+M-MP
      NS=MP-M
      A1=DSQRT(AF(N2+1)*AF(J-M+1)*AF(J+MP+1)*AF(N1+1))
c
c          determine the summation index bound
c
       KST=MAX0(0,M-MP)+1
       KEND=MIN0(N1,N2)+1                       
   10  R=0.D0
       DO 1 KK=KST,KEND
        K=KK-1
        NC1=NC-2*K
         NS1=NS+2*K
          K=KK-1
           X=(-1)**K
          IF(NC1.NE.0) X=X*C**NC1
         IF(NS1.NE.0) X=X*S**NS1
        Y=AF(N1-K+1)*AF(N2-K+1)*AF(N3+KK)*AF(KK)
   1   R=R+X/Y
       R=A1*R
       RETURN
      END
c=====================================================================
      SUBROUTINE NFACTO(N,A,IOP)
      IMPLICIT REAL*8(A-H,O-Z)
C*************************************************************************
C  IF IOP=0: CALCULATES THE NFACTORIAL AN(I)=(N-1)!
C  IF IOP GREATER THAN 0: CALCULATES THE LOG OF FACTORIAL AN(I)=DLOG((N-1)!)
C*************************************************************************
      DIMENSION A(N)
      IF(N.LT.1) STOP 'N IS LESS THAN 1 IN NFACTO'
      IF(IOP.GT.0) GO TO 3
      A(1)=1.D0
      DO 1 I=2,N
      A(I)=DFLOAT(I-1)*A(I-1)
   1  CONTINUE
      RETURN
   3  CONTINUE
      A(1)=0.D0
      DO 2 I=2,N
      A(I)=A(I-1)+DLOG(DFLOAT(I-1))
   2  CONTINUE
      RETURN
      END

c=======================================================================
      SUBROUTINE SINFT(Y,N,NA,SUM)
C   Calculates the sine transform of NA sets of data points
C   stored in array Y(NA,N).  The number of N must be a power of 2.  
C   on exit Y is replaced by its transform.  This program without 
C   changes, also calculates the inverse sine transform, but in this case 
C    the output array should be multiplied by 2/N
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 Y,Y1,Y2,SUM(NA)
      DIMENSION Y(NA,N)

      THETA=3.141592653589793D0/N
      WR=1.D0
      WI=0.D0
      WPR=-2.0D0*DSIN(0.5D0*THETA)**2
      WPI=DSIN(THETA)

      DO I=1,NA
      Y(I,1)=(0.D0,0.D0)
      END DO

      M=N/2
      DO 10 J=1,M
       WTEMP=WR
       WR=WR*WPR-WI*WPI+WR
       WI=WI*WPR+WTEMP*WPI+WI
       DO I=1,NA
       Y1=WI*(Y(I,J+1)+Y(I,N-J+1))
       Y2=0.5D0*(Y(I,J+1)-Y(I,N-J+1))
       Y(I,J+1)=Y1+Y2
       Y(I,N-J+1)=Y1-Y2
       END DO
   10 CONTINUE
      CALL REALFT(Y,M,NA,1)

      DO I=1,NA
      SUM(I)=(0.D0,0.D0)
      Y(I,1)=0.5D0*Y(I,1)
      Y(I,2)=0.D0
      END DO

      DO 20 J=1,N-1,2
       DO I=1,NA
       SUM(I)=SUM(I)+Y(I,J)
       Y(I,J)=Y(I,J+1)
       Y(I,J+1)=SUM(I)
       END DO
   20 CONTINUE

      FAC=SQRT(2.D0/N)
      DO J=1,N
      DO I=1,NA
      Y(I,J)=FAC*Y(I,J)
      END DO
      END DO
     
      RETURN
      END

      SUBROUTINE REALFT(DATA,N,NA,ISIGN)
C   Calculates the Fourier Transform of a set of 2N real-valued
C   data points.  Replaces this data (which is stored in array DATA)
C   by the positive frequency half half of its complex Fourier Transform.  
C   the real-valued first and last components of the complex transform
C   are returned as elements DATA(1) and DATA(2) respectively.  N must be a power of 2.  
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 DATA(NA,2*N),H1R,H1I,H2R,H2I
      THETA=3.141592653589793D0/DFLOAT(N)
      C1=0.5D0
      IF(ISIGN.EQ.1) THEN
      C2=-0.5D0
      CALL FOUR1(DATA,N,NA,1)
      ELSE
      C2=0.5D0
      THETA=-THETA
      ENDIF
      WPR=-2.D0*DSIN(0.5D0*THETA)**2
      WPI=DSIN(THETA)
      WR=1.D0+WPR
      WI=WPI
      N2P3=2*N+3
      DO 10 I=2,N/2
       I1=2*I-1
       I2=I1+1
       I3=N2P3-I2
       I4=I3+1
       WRS=WR
       WIS=WI
       DO J=1,NA
         H1R=C1*(DATA(J,I1)+DATA(J,I3))
         H1I=C1*(DATA(J,I2)-DATA(J,I4))
         H2R=-C2*(DATA(J,I2)+DATA(J,I4))
         H2I=C2*(DATA(J,I1)-DATA(J,I3))
         DATA(J,I1)=H1R+WRS*H2R-WIS*H2I
         DATA(J,I2)=H1I+WRS*H2I+WIS*H2R
         DATA(J,I3)=H1R-WRS*H2R+WIS*H2I
         DATA(J,I4)=-H1I+WRS*H2I+WIS*H2R
       END DO
       WTEMP=WR
       WR=WR*WPR-WI*WPI+WR
       WI=WI*WPR+WTEMP*WPI+WI
   10 CONTINUE
      IF(ISIGN.EQ.1) THEN
      DO I=1,NA
        H1R=DATA(I,1)
        DATA(I,1)=H1R+DATA(I,2)
        DATA(I,2)=H1R-DATA(I,2)
      END DO
      ELSE
      DO I=1,NA
        H1R=DATA(I,1)
        DATA(I,1)=C1*(H1R+DATA(I,2))
        DATA(I,2)=C1*(H1R-DATA(I,2))
        CALL FOUR1(DATA,N,NA,-1)
      END DO
      ENDIF
      RETURN
      END

      SUBROUTINE FOUR1(DATA,NN,NA,ISIGN)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 DATA,TEMPR,TEMPI
      DIMENSION DATA(NA,2*NN)
      N=2*NN
      J=1
      DO 10 I=1,N,2
        IF(J.GT.I) THEN
	  DO K=1,NA
            TEMPR=DATA(K,J)
            TEMPI=DATA(K,J+1)
            DATA(K,J)=DATA(K,I)
            DATA(K,J+1)=DATA(K,I+1)
            DATA(K,I)=TEMPR
            DATA(K,I+1)=TEMPI
	  END DO
        ENDIF
        M=N/2
  1     IF((M.GE.2).AND.(J.GT.M)) THEN
        J=J-M
        M=M/2
        GO TO 1
        ENDIF
        J=J+M
   10 CONTINUE
      MMAX=2
  2   IF(N.GT.MMAX) THEN
      ISTEP=2*MMAX
      THETA=6.28318530717959D0/(ISIGN*MMAX)
      WPR=-2.D0*DSIN(0.5D0*THETA)**2
      WPI=DSIN(THETA)
      WR=1.D0
      WI=0.D0
      DO 30 M=1,MMAX,2
       DO 20 I=M,N,ISTEP
         J=I+MMAX
	 DO K=1,NA
           TEMPR=WR*DATA(K,J)-WI*DATA(K,J+1)
           TEMPI=WR*DATA(K,J+1)+WI*DATA(K,J)
           DATA(K,J)=DATA(K,I)-TEMPR
           DATA(K,J+1)=DATA(K,I+1)-TEMPI
           DATA(K,I)=DATA(K,I)+TEMPR
           DATA(K,I+1)=DATA(K,I+1)+TEMPI
	 END DO
   20  CONTINUE
      WTEMP=WR
      WR=WR*WPR-WI*WPI+WR
      WI=WI*WPR+WTEMP*WPI+WI
   30 CONTINUE
      MMAX=ISTEP
      GO TO 2
      ENDIF
      RETURN
      END
c======================================================================
      SUBROUTINE SINFT1(Y,N)
C   Calculates the sine transform of NC sets of real-valued data points
C   stored in array Y(N,NC).  The number of N must be a power of 2.  
C   on exit Y is replaced by its transform.  This program without 
C   changes, also calculates the inverse sine transform, but in this case 
C    the output array should be multiplied by 2/N
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 Y(N),Y1,Y2,SUM
c      DIMENSION Y(NN,NC)
      THETA=3.141592653589793D0/N
      WR=1.D0
      WI=0.D0
      WPR=-2.0D0*DSIN(0.5D0*THETA)**2
      WPI=DSIN(THETA)
      Y(1)=(0.D0,0.D0)
      M=N/2
      DO 10 J=1,M
       WTEMP=WR
       WR=WR*WPR-WI*WPI+WR
       WI=WI*WPR+WTEMP*WPI+WI
       Y1=WI*(Y(J+1)+Y(N-J+1))
       Y2=0.5D0*(Y(J+1)-Y(N-J+1))
       Y(J+1)=Y1+Y2
       Y(N-J+1)=Y1-Y2
   10 CONTINUE
      CALL REALFT1(Y,M,1)
      SUM=0.D0
      Y(1)=0.5D0*Y(1)
      Y(2)=0.D0
      DO 20 J=1,N-1,2
       SUM=SUM+Y(J)
       Y(J)=Y(J+1)
       Y(J+1)=SUM
   20 CONTINUE

      FAC=SQRT(2.D0/N)
      DO I=1,N
         Y(I)=FAC*Y(I)
      END DO
	
      RETURN
      END
c======================================================================
      SUBROUTINE COSFT1(Y,N,ISIGN)
C   Calculates the sine transform of NC sets of real-valued data points
C   stored in array Y(N,NC).  The number of N must be a power of 2.  
C   on exit Y is replaced by its transform.  This program without 
C   changes, also calculates the inverse sine transform, but in this case 
C    the output array should be multiplied by 2/N
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 Y(N),Y1,Y2,SUM,EVEN,ODD,ENFO,SUMO,SUME
      THETA=3.141592653589793D0/N
      WR=1.D0
      WI=0.D0
      WPR=-2.0D0*DSIN(0.5D0*THETA)**2
      WPI=DSIN(THETA)
      SUM=Y(1)
      M=N/2
      DO 10 J=1,M-1
       WTEMP=WR
       WR=WR*WPR-WI*WPI+WR
       WI=WI*WPR+WTEMP*WPI+WI
       Y1=0.5D0*(Y(J+1)+Y(N-J+1))
       Y2=Y(J+1)-Y(N-J+1)
       Y(J+1)=Y1-WI*Y2
       Y(N-J+1)=Y1+WI*Y2
       SUM=SUM+WR*Y2
   10 CONTINUE
      CALL REALFT1(Y,M,1)
      Y(2)=SUM
      DO 20 J=4,N,2
       SUM=SUM+Y(J)
       Y(J)=SUM
   20 CONTINUE

      IF(ISIGN.EQ.-1) THEN
	EVEN=Y(1)
	ODD=Y(2)
  	DO 30 I=3,N-1,2
	  EVEN=EVEN+Y(I)
          ODD=ODD+Y(I+1)
   30   CONTINUE
	ENFO=2.D0*(EVEN-ODD)
	SUMO=Y(1)-ENFO
	SUME=2.D0*ODD/N-SUMO
        Y(1)=0.5*ENFO
	Y(2)=Y(2)-SUME
	DO 40 I=3,N-1,2
	  Y(I)=Y(I)-SUMO
	  Y(I+1)=Y(I+1)-SUME
   40	END DO
      END IF

c
      FAC=SQRT(2.D0/N)
      DO I=1,N
         Y(I)=FAC*Y(I)
      END DO
	
      RETURN
      END
c=============================================================
      SUBROUTINE REALFT1(DATA,N,ISIGN)
C   Calculates the Fourier Transform of a set of 2N real-valued
C   data points.  Replaces this data (which is stored in array DATA)
C   by the positive frequency half half of its complex Fourier Transform.  
C   the real-valued first and last components of the complex transform
C   are returned as elements DATA(1) and DATA(2) respectively.  N must be a power of 2.  
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 DATA(2*N),H1R,H1I,H2R,H2I
c      DIMENSION DATA(2*N)
      THETA=3.141592653589793D0/DFLOAT(N)
      C1=0.5D0
      IF(ISIGN.EQ.1) THEN
      C2=-0.5D0
      CALL FOUR11(DATA,N,1)
      ELSE
      C2=0.5D0
      THETA=-THETA
      ENDIF
      WPR=-2.D0*DSIN(0.5D0*THETA)**2
      WPI=DSIN(THETA)
      WR=1.D0+WPR
      WI=WPI
      N2P3=2*N+3
      DO 10 I=2,N/2
       I1=2*I-1
       I2=I1+1
       I3=N2P3-I2
       I4=I3+1
       WRS=WR
       WIS=WI
       H1R=C1*(DATA(I1)+DATA(I3))
       H1I=C1*(DATA(I2)-DATA(I4))
       H2R=-C2*(DATA(I2)+DATA(I4))
       H2I=C2*(DATA(I1)-DATA(I3))
       DATA(I1)=H1R+WRS*H2R-WIS*H2I
       DATA(I2)=H1I+WRS*H2I+WIS*H2R
       DATA(I3)=H1R-WRS*H2R+WIS*H2I
       DATA(I4)=-H1I+WRS*H2I+WIS*H2R
       WTEMP=WR
       WR=WR*WPR-WI*WPI+WR
       WI=WI*WPR+WTEMP*WPI+WI
   10 CONTINUE
      IF(ISIGN.EQ.1) THEN
      H1R=DATA(1)
      DATA(1)=H1R+DATA(2)
      DATA(2)=H1R-DATA(2)
      ELSE
      H1R=DATA(1)
      DATA(1)=C1*(H1R+DATA(2))
      DATA(2)=C1*(H1R-DATA(2))
      CALL FOUR11(DATA,N,-1)
      ENDIF
      RETURN
      END
      SUBROUTINE FOUR11(DATA,NN,ISIGN)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 DATA,TEMPR,TEMPI
      DIMENSION DATA(2*NN)
      N=2*NN
      J=1
      DO 10 I=1,N,2
        IF(J.GT.I) THEN
        TEMPR=DATA(J)
        TEMPI=DATA(J+1)
        DATA(J)=DATA(I)
        DATA(J+1)=DATA(I+1)
        DATA(I)=TEMPR
        DATA(I+1)=TEMPI
        ENDIF
        M=N/2
  1     IF((M.GE.2).AND.(J.GT.M)) THEN
        J=J-M
        M=M/2
        GO TO 1
        ENDIF
        J=J+M
   10 CONTINUE
      MMAX=2
  2   IF(N.GT.MMAX) THEN
      ISTEP=2*MMAX
      THETA=6.28318530717959D0/(ISIGN*MMAX)
      WPR=-2.D0*DSIN(0.5D0*THETA)**2
      WPI=DSIN(THETA)
      WR=1.D0
      WI=0.D0
      DO 30 M=1,MMAX,2
       DO 20 I=M,N,ISTEP
         J=I+MMAX
         TEMPR=WR*DATA(J)-WI*DATA(J+1)
         TEMPI=WR*DATA(J+1)+WI*DATA(J)
         DATA(J)=DATA(I)-TEMPR
         DATA(J+1)=DATA(I+1)-TEMPI
         DATA(I)=DATA(I)+TEMPR
         DATA(I+1)=DATA(I+1)+TEMPI
   20  CONTINUE
      WTEMP=WR
      WR=WR*WPR-WI*WPI+WR
      WI=WI*WPR+WTEMP*WPI+WI
   30 CONTINUE
      MMAX=ISTEP
      GO TO 2
      ENDIF
      RETURN
      END

c==========================================================
      function dm(ij,ik,im,beta)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     This function uses eq. (4.1.23) of Edmonds 
c     to calculate the reduced rotation matrix element 
c     d(j,k,m;beta) = <jk|exp(+i*beta*Jy/hbar)|jm>.
c     ----------------------------------------------------------------- 
c     
      parameter (zero = 0.0d0)
      parameter (half = 0.5d0)
      parameter (one  = 1.0d0)
      parameter (two  = 2.0d0)
c
c     half integer angular momenta
c
      sj = half*nint(two*ij)
      sk = half*nint(two*ik)
      sm = half*nint(two*im)
c
c     projection ranges
c
      dm = zero
      if (sk.gt.sj .or. sk.lt.-sj)  return
      if (sm.gt.sj .or. sm.lt.-sj)  return
      if (mod(sj-sk,one) .ne. zero) return
      if (mod(sj-sm,one) .ne. zero) return      
c
c     reflection symmetries
c      
      if (sk+sm .ge. zero) then
        if (sk-sm .ge. zero) then
          tk = sk
          tm = sm
          isign = 0
        else
          tk = sm
          tm = sk
          isign = sk-sm
        endif
      else
        if (sk-sm .ge. zero) then
          tk = -sm
          tm = -sk
          isign = 0
        else
          tk = -sk
          tm = -sm
          isign = sk-sm
        endif
      endif
c
c     evaluation
c
      n = sj-tk
      ia = tk-tm
      ib = tk+tm
      a = ia
      b = ib
      beta2 = half*beta
      cosb2 = cos(beta2)
      sinb2 = sin(beta2)
      cosb = (cosb2-sinb2)*(cosb2+sinb2)
      d1 = pjacob(n,a,b,cosb)
      d2 = cosb2**ib*sinb2**ia
      d3 = d1*d2
      d4 = d3*d3
      ti = tm
      do i = 1,ia
         ti = ti+one
         d4 = d4*(sj+ti)/(sj-ti+one)
      enddo
      d4 = sqrt(d4)
      dm = sign(d4,d3)
      if (mod(isign,2) .ne. 0) dm = -dm
      return
      end

      function pjacob (n,a,b,x)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     Jacobi polynomial p(n,a,b;x)
c     Abramowitz and Stegun eq. (22.7.1)
c     ----------------------------------------------------------------- 
c
      parameter (zero = 0.0d0)
      parameter (half = 0.5d0)
      parameter (one  = 1.0d0)
      parameter (two  = 2.0d0)
c
      if (n .eq. 0) then
        fp = one
      else
        f = one
        apa = a+a
        apb = a+b
        amb = a-b
        apbamb = apb*amb
        apbp1 = apb+one
        apbp2 = apb+two
        onek = zero
        twok = zero
        fp = half*(amb+apbp2*x)
        do k = 1,n-1
          onek = onek+one
          twok = twok+two 
          a1 = (twok+two)*(onek+apbp1)*(twok+apb)
          a2 = (twok+apbp1)*apbamb
          a3 = (twok+apb)*(twok+apbp1)*(twok+apbp2)
          a4 = (twok+apa)*(onek+b)*(twok+apbp2)
          fm = f
          f = fp
          fp = ((a2+a3*x)*f-a4*fm)/a1
        enddo
      endif
      pjacob = fp
      return
      end

