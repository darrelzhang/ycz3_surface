c==========================================================================
c be careful with the real*8 and real*4 data precision
c a safe way is to use real*8 as temporary storage
c and then assign value to the global variables after finish computation
c the real*4 global variables are:
c TMATRU, TMATAU: bound states--PODVR transformation matrix
c TMAT2D: PODVR direct product--2D H eigenvectors transformation matrix
c==========================================================================
        subroutine BOUND2D
        use COMPARAS, only : NRU, NVBRU, NAU, NVBAU, RUmin, RUmax,
     .   AUmin, AUmax, Zinf, R0, RU0, AU0, RU0, A10, A20, B10,
     .   B20, VINF 
        use REDUCED_MASS, only : TMS_3Z, RMS_CZ3 
        use QUADPOINTS, only : RUQ, AUQ
        use TRANSMATS, only : TMAT2D
        use SMALL_ARRAYS, only : EIG2D, Vref2D, HAR2D, TMATAU, TMATRU,
     .    COEFFRU, COEFFAU

        implicit none

        integer NG, N2D, i, j, k, l, m, n, ms2, ms, ns2, ns, ierr
        real*8 FAC, PI, EV, ru, au, vx, vx1, vx2, x, s, CM
        real*8, allocatable :: TRU(:,:), TAU(:,:), TMP(:,:), TMP2(:,:)
        real*8, allocatable :: VPOT_RU(:),RUQ0(:),EIGRU(:),HAR_RU(:,:)
        real*8, allocatable :: VPOT_AU(:),AUQ0(:),EIGAU(:),HAR_AU(:,:)
        real*8, allocatable :: VPOT_2D(:)
        real*8, allocatable :: TMP2D(:,:), TMP2D_2(:,:)
        real*8, allocatable :: WK(:)
        real*8 VV, VVV, ANmax, Xmin, Xmax, Smin, Smax
        real*8, allocatable :: XQ(:), SQ(:), ETX(:), ETS(:)
        real*8, allocatable :: TMATX(:,:), TMATS(:,:)
        real*8, allocatable :: TMP_x(:,:), TMP_s(:,:)
        real*8, allocatable :: VPOT2D_XS(:)
        integer NNX, NNS

        real*4 eng
        real*4, allocatable :: vec(:), WK2(:)
      
        PI=DACOS(-1.d0); EV=27.2116d0; CM=219475d0

        NG=MIN(3*MAX(NRU,NAU), 1000)

        print*, 'NRU=',NRU, 'NAU=',NAU
        print*, 'NVBRU=',NVBRU, 'NVBAU=',NVBAU
        print*, 'NG=',NG

! solve bound state on (RU,AU)
        allocate(RUQ0(NRU),EIGRU(NRU),AUQ0(NAU),EIGAU(NAU),stat=ierr)
        allocate(TRU(NRU,NRU),TAU(NAU,NAU),stat=ierr)
        call QUAD1D(NRU,RUmin,RUmax,RUQ0)
        call QUAD1D(NAU,AUmin,AUmax,AUQ0)

        FAC=DSQRT(2.D0/(NRU+1))
        DO I=1,NRU; DO J=1,NRU
          TRU(J,I)=FAC*DSIN(I*J*PI/(NRU+1))
        END DO; END DO

        FAC=DSQRT(2.D0/(NAU+1))
        DO I=1,NAU; DO J=1,NAU
          TAU(J,I)=FAC*DSIN(I*J*PI/(NAU+1))
        END DO; END DO

        ms=0
        do j=1, NRU; ru=RUQ0(j)
        do i=1, NAU; au=AUQ0(i)
          ms=ms+1
          call POTDRIVA(Zinf,R0,ru,au,A10,A20,B10,B20,vx)
          write(60,*) ru, au, (vx-VINF)*EV
        end do; write(60,"()"); end do
 
! 1D bound states for RU at AU0
        allocate(VPOT_RU(NRU),stat=ierr)
        do i=1, NRU; ru=RUQ0(i)
          call POTDRIVA(Zinf,R0,ru,AU0,A10,A20,B10,B20,vx)
          VPOT_RU(i)=vx-VINF
          write(61,*) ru, VPOT_RU(i)*EV
        end do
        open(100, file='RUPOT.1D',form='unformatted',status='unknown',
     .    iostat=ierr)
        write(100) (VPOT_RU(i),i=1,NRU)
        close(100)
! T(ru,au0)
        allocate(HAR_RU(NRU,NRU), stat=ierr)
        call KMATRIX_RU_FBR(NG,NRU,RUQ0,AU0,HAR_RU,RMS_CZ3,TMS_3Z)
! transfer to sine DVR
        allocate(TMP(NRU,NRU),stat=ierr)
        TMP=0.d0
        do k=1, NRU
        do j=1, NRU
        do i=1, NRU
          TMP(j,k)=TMP(j,k)+TRU(i,j)*HAR_RU(i,k)
        end do
        end do
        end do
        HAR_RU=0.d0
        do k=1, NRU
        do j=1, NRU
        do i=1, NRU
          HAR_RU(i,k)=HAR_RU(i,k)+TRU(j,k)*TMP(i,j)
        end do
        end do
        end do
! H(RU,AU0)
        do i=1, NRU; HAR_RU(i,i)=HAR_RU(i,i)+VPOT_RU(i); end do
        deallocate(TMP)

! ru PODVR
        allocate(RUQ(NVBRU),stat=ierr)
        allocate(TMP(NVBRU,NVBRU),stat=ierr)
        allocate(COEFFRU(NRU,NRU),TMATRU(NVBRU,NVBRU),stat=ierr)
        call VIB1D_P(NRU,NVBRU,RUQ0,RUQ,HAR_RU,EIGRU,COEFFRU,TMP)
        TMATRU=SNGL(TMP)  ! real*8 -> real*4
        print*, 'RU PODVR points'
        do i=1, NVBRU; print*, i, RUQ(i); end do
        print*

        print*, 'RU bound states value'
        do i=1, min(NRU,10)
           print*, i, EIGRU(i)*CM, (EIGRU(i)-EIGRU(1))*CM
        end do
        !do i=1, min(NRU,10); print*, i, EIGRU(i)*EV; end do
        print*

        deallocate(TMP)

! 1D bound states for AU
        allocate(VPOT_AU(NAU),stat=ierr)
        do i=1, NAU; au=AUQ0(i)
          call POTDRIVA(Zinf,R0,RU0,au,A10,A20,B10,B20,vx)
          VPOT_AU(i)=vx-VINF
          write(62,*) au, VPOT_AU(i)*EV
        end do
        open(100, file='AUPOT.1D',form='unformatted',status=
     .    'unknown',iostat=ierr)
        write(100) (VPOT_AU(i),i=1,NAU)
        close(100)

        allocate(HAR_AU(NAU,NAU),stat=ierr)
        call KMATRIX_AU_FBR(NG,NAU,AUQ0,RU0,AU0,HAR_AU,RMS_CZ3,TMS_3Z)
! T(ru0,au) in sine basis
        allocate(TMP(NAU,NAU),stat=ierr)
        TMP=0.d0
        do k=1, NAU
        do j=1, NAU
        do i=1, NAU
          TMP(j,k)=TMP(j,k)+TAU(i,j)*HAR_AU(i,k)
        end do
        end do
        end do
        HAR_AU=0.d0
        do k=1, NAU
        do j=1, NAU
        do i=1, NAU
          HAR_AU(i,k)=HAR_AU(i,k)+TAU(j,k)*TMP(i,j)
        end do
        end do
        end do

! H(RU0,AU)
        do i=1, NAU; HAR_AU(i,i)=HAR_AU(i,i)+VPOT_AU(i); end do

        deallocate(TMP)

! PODVR
        allocate(AUQ(NVBAU),stat=ierr)
        allocate(TMP(NVBAU,NVBAU),stat=ierr)
        allocate(COEFFAU(NAU,NAU),TMATAU(NVBAU,NVBAU),stat=ierr)
        call VIB1D_P(NAU,NVBAU,AUQ0,AUQ,HAR_AU,EIGAU,COEFFAU,TMP)
        TMATAU=SNGL(TMP)  ! real*4
        print*, 'AU PODVR points'
        do i=1, NVBAU; print*, i, AUQ(i); end do
        print*

        print*, 'AU bound states value'
        do i=1, min(NAU,10)
          print*, i, EIGAU(i)*CM, (EIGAU(i)-EIGAU(1))*CM; 
        end do
        !do i=1, min(NAU,10); print*, i, EIGAU(i)*EV; end do
        print*

        deallocate(TMP)

! HAR2D to hold 2D Hamiltonian
        N2D=NRU*NAU
        allocate(HAR2D(N2D,N2D),stat=ierr)
        HAR2D=0.d0
        call HAR2D_RESIDUAL(NRU,NAU,RU0,AU0,RUQ0,AUQ0,HAR2D,RMS_CZ3,TMS_3Z)

        N2D=NVBRU*NVBAU
        allocate(TMP(N2D,N2D),EIG2D(N2D),stat=ierr)
! transfer HAR2D to BS by COEFFRU and COEFFAU
        call SINETOBS(NRU,NAU,NVBRU,NVBAU,COEFFRU,COEFFAU,HAR2D,TMP)

! add two 1D bound state values
        ms=0
        do j=1, NVBRU
        do i=1, NVBAU
          ms=ms+1
          TMP(ms,ms)=TMP(ms,ms)+EIGRU(j)+EIGAU(i)
        end do
        end do

! transfer TMP to PODVR
        call BSTOPODVR(NVBRU,NVBAU,TMATRU,TMATAU,TMP)

! 1D reference potential in PODVR points
        deallocate(VPOT_RU, VPOT_AU)
        allocate(VPOT_RU(NVBRU), stat=ierr)
        do i=1, NVBRU; ru=RUQ(i)
          call POTDRIVA(Zinf,R0,ru,AU0,A10,A20,B10,B20,vx)
          VPOT_RU(i)=vx-VINF
        end do
        allocate(VPOT_AU(NVBAU), stat=ierr)
        do i=1, NVBAU; au=AUQ(i)
          call POTDRIVA(Zinf,R0,RU0,au,A10,A20,B10,B20,vx)
          VPOT_AU(i)=vx-VINF
        end do

        allocate(Vref2D(N2D),stat=ierr)
        ms=0
        do j=1, NVBRU; ru=RUQ(j)
        do i=1, NVBAU; au=AUQ(i)
          ms=ms+1
          call POTDRIVA(Zinf,R0,ru,au,A10,A20,B10,B20,vx)
          Vref2D(ms)=vx-VPOT_AU(i)-VPOT_RU(j)-VINF
        end do; end do
        open(100, file='POT.2D',form='unformatted',status=
     .    'unknown',iostat=ierr)
        write(100) (Vref2D(i),i=1,NVBRU*NVBAU)
        close(100)

        do i=1, NVBRU*NVBAU; TMP(i,i)=TMP(i,i)+Vref2D(i); end do

        deallocate(HAR2D)
        allocate(HAR2D(N2D,N2D),stat=ierr)
        HAR2D=TMP
        print*,'H(ru,au) eigenvalue in PODVR by diagonalization'
        call EV_TEST(N2D,HAR2D,EIG2D)
        print*

! assign to TMAT2D to make change of  the data format
        allocate(TMAT2D(N2D,N2D),stat=ierr)
! need a transpose to follow the usual way of using the ztrans subroutines
        do j=1, N2D; do i=1, N2D
          TMAT2D(i,j)=SNGL(HAR2D(j,i))
        end do; end do
 
        ms=0
        do j=1, NVBRU; ru=RUQ(j)
        do i=1, NVBAU; au=AUQ(i)
          ms=ms+1
          call POTDRIVA(Zinf,R0,ru,au,A10,A20,B10,B20,vx)
          Vref2D(ms)=vx-VINF
!          write(223,*) ms, Vref2D(ms)
        end do; end do
        open(100, file='POT.PODVR.2D',form='unformatted',status=
     .    'unknown',iostat=ierr)
        write(100) (Vref2D(i),i=1,NVBRU*NVBAU)
        close(100)

! transfer to bound state basis
        allocate(TMP2(N2D,N2D),stat=ierr)
        TMP2=0.d0
        do l=1, N2D
        do k=1, NVBRU
        do j=1, NVBRU
        do i=1, NVBAU
          ms=(j-1)*NVBAU+i; ms2=(k-1)*NVBAU+i
          TMP2(ms,l)=TMP2(ms,l)+HAR2D(ms2,l)*TMATRU(j,k)
        end do
        end do
        end do
        end do

        HAR2D=0.d0
        do l=1, N2D
        do k=1, NVBRU
        do j=1, NVBAU
        do i=1, NVBAU
          ms=(k-1)*NVBAU+i; ms2=(k-1)*NVBAU+j
          HAR2D(ms,l)=HAR2D(ms,l)+TMP2(ms2,l)*TMATAU(i,j)
        end do
        end do
        end do
        end do

! check expansion coefficients
        do k=1, N2D
          write(110,*) 'state=', k
        do j=1, NVBRU
        do i=1, NVBAU
          ms=(j-1)*NVBAU+i
          write(110,*) j, i, abs(HAR2D(ms,k))**2.d0
        end do
        end do
          write(110,"()")
        end do

        return

        print*,'H(ru,au) eigenvalue in PODVR by Lanczos'
        allocate(vec(N2D),stat=ierr)
        call BOUND2D_ARPACK(eng, vec, TMP)
       
        allocate(WK2(N2D), stat=ierr)
        WK2=0.d0
        do k=1, NVBRU
        do j=1, NVBRU
        do i=1, NVBAU
          ms=(j-1)*NVBAU+i; ms2=(k-1)*NVBAU+i
          WK2(ms)=WK2(ms)+vec(ms2)*TMATRU(j,k)
        end do
        end do
        end do

        vec=0.d0
        do k=1, NVBRU
        do j=1, NVBAU
        do i=1, NVBAU
          ms=(k-1)*NVBAU+i; ms2=(k-1)*NVBAU+j
          vec(ms)=vec(ms)+WK2(ms2)*TMATAU(i,j)
        end do
        end do
        end do

! check expansion coefficients
        do j=1, NVBRU
        do i=1, NVBAU
          ms=(j-1)*NVBAU+i
          write(115,*) j, i, abs(vec(ms))**2.d0
        end do
        end do

        end subroutine BOUND2D
c================================================================
c this subroutine constructs the 2nd order derivatives of r
c while au=au0
c==================================================================
       subroutine KMATRIX_RU_FBR(NG,NRU,RUQ,AU0,HAR_RU,mus,mux)
       implicit none
       integer NRU
       real*8 RUQ(NRU), HAR_RU(NRU,NRU)
       real*8, allocatable :: TMP_RU(:,:), TMP(:,:), TRUQ(:)

       integer i, j, k, l, ms, ms2, ierr
       real*8 PI, RUmin, RUmax, RUL, DRU
       real*8 ru, mus, mux, FAC, sum, TDRU, AU0
       real*8 A, B, C
       integer NG

       PI=DACOS(-1.d0)
       DRU=RUQ(2)-RUQ(1);RUmin=RUQ(1)-DRU;RUmax=RUQ(NRU)+DRU;RUL=RUmax-RUmin
       A=sin(au0)**2.d0/(2.d0*mux)+cos(au0)**2.d0/(2.d0*mus)
       B=cos(au0)**2.d0/(2.d0*mux)+sin(au0)**2.d0/(2.d0*mus)
       C=1.d0/mux-1.d0/mus

       allocate(TMP_RU(NRU,NRU),TMP(NRU,NRU),stat=ierr)

       allocate(TRUQ(NG), stat=ierr)

       HAR_RU(:,:)=0.d0

       TDRU=RUL/(NG+1)
       do i=1, NG; TRUQ(i)=RUmin+i*TDRU; end do

! second order derivative of ru
       do k=1, NRU; do j=1, NRU
         sum=0.d0
       do i=1, NG
         sum=sum+DSIN(k*PI*i/(NG+1))*DSIN(j*PI*i/(NG+1))
       end do
         TMP_RU(j,k)=sum*2.d0/(NG+1)*(j*PI/RUL)**2.d0*A
       end do; end do

       HAR_RU=TMP_RU

       end subroutine KMATRIX_RU_FBR
c==================================================================
c this subroutine is to calculate the 2nd order derivatives
c of au while ru=ru0
c================================================================== 
       subroutine KMATRIX_AU_FBR(NG,NAU,AUQ,RU0,AU0,HAR_AU,mus,mux)
       implicit none
       integer NAU
       real*8 AUQ(NAU), HAR_AU(NAU,NAU)
       real*8, allocatable :: TMP_AU(:,:), TMP(:,:), TAUQ(:)

       integer i, j, k, l, ms, ms2, ierr
       integer NG
       real*8 PI, AUmin, AUmax, AUL, DAU, A, B, C
       real*8 au, mus, mux, FAC, sum, TDAU, RU0, AU0

       PI=DACOS(-1.d0)
       DAU=AUQ(2)-AUQ(1);AUmin=AUQ(1)-DAU;AUmax=AUQ(NAU)+DAU;AUL=AUmax-AUmin
       C=1.d0/mux-1.d0/mus
       B=cos(au0)**2.d0/(2.d0*mux)+sin(au0)**2.d0/(2.d0*mus)

       allocate(TMP_AU(NAU,NAU),TMP(NAU,NAU),stat=ierr)

       HAR_AU(:,:)=0.d0
       
       allocate(TAUQ(NG), stat=ierr)
       TDAU=AUL/(NG+1)
       do i=1, NG; TAUQ(i)=AUmin+i*TDAU; end do

! second order derivative of au
       do k=1, NAU; do j=1, NAU
         sum=0.d0
       do i=1, NG; au=TAUQ(i)
         sum=sum+DSIN(k*PI*i/(NG+1))*DSIN(j*PI*i/(NG+1))
       end do
         TMP_AU(j,k)=sum*2.d0/(NG+1)*(j*PI/AUL)**2.d0/RU0**2.d0*B
       end do; end do
 
       HAR_AU=TMP_AU

       end subroutine KMATRIX_AU_FBR
c====================================================================
       subroutine VIB1D_P(NR,NVB,RQ,DVR,HAR,EIG,COEFF,TMAT)
       implicit none
       integer NR, NVB, INFO, LWORK, i, ierr, NFOUND
       real*8 RQ(NR), DVR(NVB), EIG(NR) 
       real*8 HAR(NR,NR), TMAT(NVB,NVB), COEFF(NR,NR)
       real*8, allocatable :: WK1(:), WK2(:)
       integer, allocatable :: IFAIL(:)
       real*8 ABSTOL, VL, VU

       ! diagonalize HAR
       ABSTOL=0.d0; LWORK=8*NR; VL=0.d0; VU=0.d0
       allocate(WK1(8*NR),WK2(5*NR),IFAIL(NR),stat=ierr)
       call DSYEVX('V','I','U',NR,HAR,NR,VL,VU,1,NR,
     $    ABSTOL,NFOUND,EIG,HAR,NR,WK1,LWORK,WK2,IFAIL,INFO)

       ! projection to sine basis
       call BSTOSINE(NR,RQ,HAR,COEFF)

       ! PODVR
       call PODVR_B(NR,NVB,HAR,TMAT,RQ,DVR)

       end subroutine VIB1D_P
c===============================================================
       subroutine BSTOSINE(NRU,RUQ,HAR,COEFF)
       implicit none
       integer NRU
       real*8 RUQ(NRU), HAR(NRU,NRU), COEFF(NRU,NRU)
       real*8 FAC,PI, DR, AN
       integer I, J, K

        FAC=SQRT(2.D0/(NRU+1)); PI=ACOS(-1.D0); DR=RUQ(2)-RUQ(1)

        DO I=1,NRU; DO J=1,NRU; AN=0.D0
          DO K=1,NRU
            AN=AN+HAR(K,I)*FAC*SIN(J*K*PI/(NRU+1))
          END DO
          COEFF(J,I)=AN
        END DO; END DO

       end subroutine BSTOSINE
c====================================================================
        SUBROUTINE PODVR_B(NR,NVB,WAVE,TMAT,RQ1,DVR)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION WAVE(NR,1),TMAT(NVB,1),RQ1(1),DVR(1),
     $  WK1(NR),WK2(NR)

c construct X-matrix.
        DO 20 I=1,NVB
        DO 20 J=1,NVB
           TMAT(J,I)=0.D0
        DO K=1,NR
        TMAT(J,I)=TMAT(J,I)+WAVE(K,J)*WAVE(K,I)*RQ1(K)
        END DO
20      CONTINUE

        CALL RS(NVB,NVB,TMAT,DVR,1,TMAT,WK1,WK2,IERR)
        IF(IERR.NE.0) STOP 'PODVR'
        DO I=1,NVB
           IF(TMAT(1,I).LT.0.D0) THEN
           DO J=1,NVB
           TMAT(J,I)=-TMAT(J,I)
           END DO
           END IF
        END DO

        end subroutine PODVR_B
c==========================================================
       subroutine EV_TEST(N,HAR,EIG)
       implicit none
       integer N,NFOUND,LWORK,ierr,INFO,i
       real*8 HAR(N,N),EIG(N),ABSTOL,VL,VU,EV,CM
       real*8, allocatable :: WK1(:), WK2(:)
       integer, allocatable :: IFAIL(:)

       EV=27.2116d0; CM=219475d0
       
       ABSTOL=0.d0; LWORK=8*N; VL=0.d0; VU=0.d0
       allocate(WK1(8*N),WK2(5*N),IFAIL(N),stat=ierr)
        CALL DSYEVX('V','I','U',N,HAR,N,VL,VU,1,N,
     $    ABSTOL,NFOUND,EIG,HAR,N,WK1,LWORK,WK2,IFAIL,INFO)

       !call RS(N,N,HAR,EIG,1,HAR,WK1,WK2,IERR)
       
       do i=1, 10
         !print*, i, EIG(i)*EV, (EIG(i)-EIG(1))*EV
         print*, i, EIG(i)*CM, (EIG(i)-EIG(1))*CM
       end do
       print*

       end subroutine EV_TEST
c=====================================================================
       subroutine BSTOPODVR(NVBRU,NVBAU,TMATRU,TMATAU,TMP)
       implicit none
       integer NVBRU,NVBAU,i,j,k,l,m,n,ms,ms2,ns,ns2,ierr
       real*4 TMATRU(NVBRU,NVBRU),TMATAU(NVBAU,NVBAU)
       real*8 TMP(NVBRU*NVBAU,NVBRU*NVBAU)
       real*8, allocatable :: WK(:,:)

       allocate(WK(NVBRU*NVBAU,NVBRU*NVBAU),stat=ierr)
       WK=0.d0
       do n=1, NVBRU; do m=1, NVBRU
       do l=1, NVBRU; do k=1, NVBRU
       do j=1, NVBAU; do i=1, NVBAU       
         ms=(l-1)*NVBAU+j; ms2=(k-1)*NVBAU+i       
         ns=(n-1)*NVBAU+j; ns2=(m-1)*NVBAU+i
         WK(ns2,ns)=WK(ns2,ns)+TMATRU(k,m)*TMP(ms2,ms)*TMATRU(l,n)
       end do; end do
       end do; end do
       end do; end do

       TMP=0.d0
       do n=1, NVBRU; do m=1, NVBRU
       do l=1, NVBAU; do k=1, NVBAU
       do j=1, NVBAU; do i=1, NVBAU       
         ms=(n-1)*NVBAU+j; ms2=(m-1)*NVBAU+i       
         ns=(n-1)*NVBAU+l; ns2=(m-1)*NVBAU+k
         TMP(ns2,ns)=TMP(ns2,ns)+TMATAU(i,k)*WK(ms2,ms)*TMATAU(j,l)
       end do; end do
       end do; end do
       end do; end do

       end subroutine BSTOPODVR
c======================================================================
       subroutine SINETOBS(NRU,NAU,NVBRU,NVBAU,TRU,TAU,HAR2D,TMP)
       implicit none
       integer NRU,NAU,NVBRU,NVBAU,i,j,k,l,m, n, ms, ms2, ns, ns2,ierr
       real*8 TRU(NRU,NVBRU), TAU(NAU,NVBAU)
       real*8 HAR2D(NRU*NAU,NRU*NAU),TMP(NVBRU*NVBAU,NVBRU*NVBAU)
       real*8, allocatable :: WK(:,:)

       allocate(WK(NVBRU*NAU,NVBRU*NAU),stat=ierr)
       WK=0.d0
       do n=1, NVBRU; do m=1, NVBRU
       do l=1, NRU; do k=1, NRU
       do j=1, NAU; do i=1, NAU       
         ms=(l-1)*NAU+j; ms2=(k-1)*NAU+i       
         ns=(n-1)*NAU+j; ns2=(m-1)*NAU+i
         WK(ns2,ns)=WK(ns2,ns)+TRU(k,m)*HAR2D(ms2,ms)*TRU(l,n)
       end do; end do
       end do; end do
       end do; end do

       TMP=0.d0
       do n=1, NVBRU; do m=1, NVBRU
       do l=1, NVBAU; do k=1, NVBAU
       do j=1, NAU; do i=1, NAU       
         ms=(n-1)*NAU+j; ms2=(m-1)*NAU+i       
         ns=(n-1)*NVBAU+l; ns2=(m-1)*NVBAU+k
         TMP(ns2,ns)=TMP(ns2,ns)+TAU(i,k)*WK(ms2,ms)*TAU(j,l)
       end do; end do
       end do; end do
       end do; end do

       end subroutine SINETOBS
c=========================================================================
        subroutine HAR2D_RESIDUAL(NRU,NAU,RU0,AU0,RUQ,AUQ,HAR2D,mus,mux)
        implicit none
        integer NRU, NAU
        real*8 RUQ(NRU),AUQ(NAU)
        real*8 HAR2D(NRU*NAU,NRU*NAU)
        real*8, allocatable :: TMP_RU(:,:), TMP_AU(:,:)
        real*8, allocatable :: TMP(:,:), TMP2(:,:), TMP3(:,:)
        real*8, allocatable :: T(:,:)
        real*8, allocatable :: TRUQ(:), TAUQ(:)
 
        real*8 RU0, AU0
        real*8 PI, RUmin, RUmax, RUL, AUmin, AUmax, AUL, DRU, DAU
        real*8 ru, au, mus, mux, FAC, FAC1, sum, TDAU, TDRU
        real*8 A, B, C
        integer NG
        integer i, j, k, l, ms, ms2, ierr
 
        PI=DACOS(-1.d0)
        C=1.d0/mux-1.d0/mus
        DRU=RUQ(2)-RUQ(1)
        RUmin=RUQ(1)-DRU;RUmax=RUQ(NRU)+DRU;RUL=RUmax-RUmin
        DAU=AUQ(2)-AUQ(1)
        AUmin=AUQ(1)-DAU;AUmax=AUQ(NAU)+DAU;AUL=AUmax-AUmin
 
! the matrix element K=<lm|T(ru,au)|nk> on sine FBR basis
        allocate(T(NRU*NAU,NRU*NAU),stat=ierr)
        allocate(TMP_RU(NRU,NRU),TMP_AU(NAU,NAU),stat=ierr)
        allocate(TMP(NRU,NRU),TMP2(NAU,NAU),TMP3(NAU,NAU),stat=ierr)
 
        NG=MAX(10*NRU, 10*NAU)
        allocate(TRUQ(NG), TAUQ(NG), stat=ierr)
        HAR2D(:,:)=0.d0
        TDRU=RUL/(NG+1)
        do i=1, NG; TRUQ(i)=RUmin+i*TDRU; end do
        TDAU=AUL/(NG+1)
        do i=1, NG; TAUQ(i)=AUmin+i*TDAU; end do
 
! part I (second order derivative of ru)
        do k=1, NRU; do j=1, NRU
          sum=0.d0
        do i=1, NG
          sum=sum+DSIN(k*PI*i/(NG+1))*DSIN(j*PI*i/(NG+1))
        end do
          TMP_RU(j,k)=sum*2.d0/(NG+1)*(j*PI/RUL)**2.d0
        end do; end do
 
! the factor A(au) 
        do k=1, NAU; do j=1, NAU
          sum=0.d0
        do i=1, NG; au=TAUQ(i)
          FAC=0.5d0*(DCOS(au)**2.d0/mus+DSIN(au)**2.d0/mux)
          sum=sum+FAC*DSIN(k*i*PI/(NG+1))*DSIN(j*i*PI/(NG+1))
        end do
          TMP_AU(j,k)=sum*2.d0/(NG+1)
        end do; end do

!  minus A(au0)
        FAC=0.5d0*(DCOS(AU0)**2.d0/mus+DSIN(AU0)**2.d0/mux)
        do i=1, NAU
          TMP_AU(i,i)=TMP_AU(i,i)-FAC
        end do
       
! direct product
        T(:,:)=0.d0
        do l=1, NRU; do k=1, NRU
        do j=1, NAU; do i=1, NAU
          ms=(l-1)*NAU+j; ms2=(k-1)*NAU+i
          T(ms2,ms)=TMP_RU(k,l)*TMP_AU(i,j)
        end do; end do
        end do; end do

        HAR2D=T  ! the residual of second order term

! part II (first order derivative of ru) 
        do k=1, NRU; do j=1, NRU
          sum=0.d0
        do i=1, NG; ru=TRUQ(i)
          sum=sum+DSIN(k*i*PI/(NG+1))*DCOS(j*i*PI/(NG+1))/ru
        end do
          TMP_RU(j,k)=sum*2.d0/(NG+1)*(j*PI/RUL)
        end do; end do

        do k=1, NAU; do j=1, NAU
          sum=0.d0
        do i=1, NG; au=TAUQ(i)
          FAC=-0.5d0*DCOS(2.d0*au)*C
          sum=sum+DSIN(k*PI*i/(NG+1))*FAC*DSIN(j*PI*i/(NG+1))
        end do
          TMP_AU(j,k)=sum*2.d0/(NG+1)
        end do; end do

        T(:,:)=0.d0
        do l=1, NRU; do k=1, NRU
        do j=1, NAU; do i=1, NAU
          ms=(l-1)*NAU+j; ms2=(k-1)*NAU+i
          T(ms2,ms)=TMP_RU(k,l)*TMP_AU(i,j)
        end do; end do
        end do; end do

        HAR2D=HAR2D+T

! part III  (cross term)
        do k=1, NRU; do j=1, NRU
          sum=0.d0
        do i=1, NG; ru=TRUQ(i)
          sum=sum+DSIN(k*i*PI/(NG+1))*DCOS(j*i*PI/(NG+1))/ru
        end do
          TMP_RU(j,k)=sum*2.d0/(NG+1)*(j*PI/RUL)
        end do; end do

        do k=1, NAU; do j=1, NAU
          sum=0.d0
        do i=1, NG; au=TAUQ(i)
          FAC=0.5d0*dsin(2.d0*au)*(1.d0/mus-1.d0/mux)
          sum=sum+DSIN(k*i*PI/(NG+1))*FAC*DCOS(j*i*PI/(NG+1))
        end do
          TMP_AU(j,k)=sum*2.d0/(NG+1)*(j*PI/AUL)
        end do; end do
       
        T(:,:)=0.d0
        do ms=1, NRU*NAU; l=(ms-1)/NAU+1; j=MOD(ms-1, NAU)+1
        do ms2=1, NRU*NAU; k=(ms2-1)/NAU+1; i=MOD(ms2-1, NAU)+1
          T(ms2,ms)=TMP_RU(k,l)*TMP_AU(i,j)
        end do; end do

        HAR2D=HAR2D+T

100     continue

! part IV (second order derivative of au)
        do k=1, NAU; do j=1, NAU
          sum=0.d0
        do i=1, NG; au=TAUQ(i)
          FAC=DSIN(au)**2.d0/mus+DCOS(au)**2.d0/mux
          sum=sum+DSIN(k*PI*i/(NG+1))*FAC*DSIN(j*PI*i/(NG+1))
        end do
          TMP_AU(j,k)=0.5d0*sum*2.d0/(NG+1)*(j*PI/AUL)**2.d0
        end do; end do

! the factor r**(-2.d0)
        do k=1, NRU; do j=1, NRU
          sum=0.d0
        do i=1, NG; ru=TRUQ(i)
          sum=sum+DSIN(k*i*PI/(NG+1))*DSIN(j*i*PI/(NG+1))/ru**2.d0
        end do
          TMP_RU(j,k)=sum*2.d0/(NG+1)
        end do; end do


        T(:,:)=0.d0
        do ms=1, NRU*NAU; l=(ms-1)/NAU+1; j=MOD(ms-1, NAU)+1
        do ms2=1, NRU*NAU; k=(ms2-1)/NAU+1; i=MOD(ms2-1, NAU)+1
          T(ms2,ms)=TMP_RU(k,l)*TMP_AU(i,j)
        end do; end do

! FAC1 is from KMAT_AU_FBR
        FAC1=(DSIN(au0)**2.d0/mus+DCOS(au0)**2.d0/mux)/RU0**2.d0
        FAC1=0.5d0*FAC1
        ms=0
        do j=1, NRU
        do i=1, NAU
          ms=ms+1
          T(ms,ms)=T(ms,ms)-FAC1*(i*PI/AUL)**2.d0
        end do
        end do
 
        HAR2D=HAR2D+T

! Part V (first order derivative of au)
        do k=1, NAU; do j=1, NAU
          sum=0.d0
        do i=1, NG; au=TAUQ(i)
          FAC=DSIN(au)*DCOS(au)
          sum=sum+DSIN(k*i*PI/(NG+1))*FAC*DCOS(j*i*PI/(NG+1))
        end do
          TMP_AU(j,k)=sum*2.d0/(NG+1)*(j*PI/AUL)*1.5d0*(1.d0/mux-1.d0/mus)
        end do; end do
 
        T(:,:)=0.d0
        do ms=1, NRU*NAU; l=(ms-1)/NAU+1; j=MOD(ms-1, NAU)+1
        do ms2=1, NRU*NAU; k=(ms2-1)/NAU+1; i=MOD(ms2-1, NAU)+1
          T(ms2,ms)=TMP_RU(k,l)*TMP_AU(i,j)
        end do; end do

        HAR2D=HAR2D+T

! grid term
        TMP_AU=0.d0
        do k=1, NAU; do j=1, NAU
          sum=0.d0
        do i=1, NG; au=TAUQ(i)
          A=sin(au)**2.d0/(2.d0*mux)+cos(au)**2.d0/(2.d0*mus)
          B=SIN(au)**2.d0/(2.d0*mus)+COS(au)**2.d0/(2.d0*mux)
          FAC=(2.d0*B-3.d0*A)/4.d0
          sum=sum+DSIN(k*PI*i/(NG+1))*FAC*DSIN(j*PI*i/(NG+1))
        end do
          TMP_AU(j,k)=sum*2.d0/(NG+1)
        end do; end do

        T(:,:)=0.d0
        do ms=1, NRU*NAU; l=(ms-1)/NAU+1; j=MOD(ms-1, NAU)+1
        do ms2=1, NRU*NAU; k=(ms2-1)/NAU+1; i=MOD(ms2-1, NAU)+1
          T(ms2,ms)=TMP_RU(k,l)*TMP_AU(i,j)
        end do; end do

        HAR2D=HAR2D+T

        end subroutine HAR2D_RESIDUAL
c======================================================================================
        subroutine BOUND2D_ARPACK(eng, vec, HAR2D)
        use COMPARAS, only : istat1, N2D
        implicit none

        integer nbasis, nroots, memuse, maxcyc
        real*4 tolvib, eng
        real*4 vec(1)
        real*8 HAR2D(1)
        

        nbasis=N2D
        nroots=istat1*2
        memuse=istat1*4
        maxcyc=4000
        tolvib=1.0d-6

        print*, 'BOUND2D_ARPACK'
        print*, '============================================='
        print*, 'nbasis=', nbasis
        print*, 'nroots=', nroots
        print*, 'memuse=', memuse
        print*, 'maxcyc=', maxcyc
        print*, 'tolvib=', tolvib
        print*, '============================================='
        print*

        rewind(21)
        call SSARPACK(21,nbasis,nroots,memuse,maxcyc,tolvib,istat1,
     .    eng, vec, HAR2D)
        !stop

        end subroutine BOUND2D_ARPACK
c=====================================================================
        subroutine HMAT2D_ARPACK(WF1, WF2, HAR2D)
        use COMPARAS, only : N2D
        use TRANSMATS, only : TMAT2D
        use SMALL_ARRAYS, only : EIG2D

        implicit none

        integer i, j, k,  ms, ns, ierr
        real*4 tkval, WF1(1), WF2(1)
        real*4, allocatable :: WF3(:)
        real*8 AN
        real*8 HAR2D(N2D,1)

        do i=1, N2D
          WF2(i)=WF1(i)
        end do

!        allocate(WF3(N2D),stat=ierr)

!        WF3=0.d0
!        do j=1, N2D
!        do i=1, N2D
!          WF3(j)=WF3(j)+WF2(i)*TMAT2D(i,j)
!        end do
!        end do
!        WF2=WF3

!        WF3=0.d0
!        do j=1, N2D
!        do i=1, N2D
!          WF3(i)=WF3(i)+WF2(j)*TMAT2D(i,j)
!        end do
!        end do
!        WF2=WF3

! H=T+V
!        do j=1, N2D
!          tkVal=EIG2D(j)
!          WF2(j)=WF2(j)*tkVal
!        end do

!        deallocate(WF3, stat=ierr)
        
        do j=1, N2D
          WF2(j)=0.d0
        do i=1, N2D
          WF2(j)=WF2(j)+WF1(i)*HAR2D(i,j)
        end do
        end do

        end subroutine HMAT2D_ARPACK 
c========================================================================
        subroutine HMATDRIV(WF1, WF2, HAR2D)
        implicit none
        real*4 WF1(1), WF2(1)
        real*8 HAR2D(1)

        call HMAT2D_ARPACK(WF1, WF2, HAR2D)

        end subroutine HMATDRIV
c======================================================================
