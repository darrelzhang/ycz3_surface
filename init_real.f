        subroutine WFINIT
        use WFARRAYS
        use ENGARRAYS, only : ESK0, AE0
        use COMPARAS, only : NZ, NZINT, NRB, NVB, N2D, NBASJJ, 
     .    NBASJJ14, NKBZinit, J1init, KBZmin, KBZmax, E1, E2, NE0,
     &    IPOT
        implicit none
        integer maxnodes, ierr, i, j, k, l, ms, ns, ns2
        integer j12Val, j2Val, ik2, kb2Val, kbzVal, ikbz, ibasJJ
        real*8 DE, AN, EV

        complex*16 WGP
        complex*16, allocatable :: W1(:)

	if(IPOT /= 1) return

        EV=27.2116d0
        allocate(ESK0(NE0),AE0(NE0),stat=ierr)
        
        print*, '=================================='
        print*, 'Parameters for initial wavefunction'
        print*
        print*,  'E1=',E1*EV,'E2=',E2*EV

! scattering energy array
        DE=(E2-E1)/(NE0+1)
        do i = 1, NE0
          ESK0(i) = E1 + (i-1)*DE
        end do

!        open(1, file='scatt_energy',status='old')
!        do i=1, NE0
!          read(1,*) ESK0(i)
!          ESK0(i)=ESK0(i)/27.2116d0
!          print*, i, ESK0(i)*27.2116
!        end do
!        close(1)


!        allocate(W1(NZ*NKBZinit), stat=ierr)
!        call GAUSS_PACKET(W1)

        call BS_YCZ3

!        allocate(WA(NZ*NVB*N2D*NBASJJ), stat=ierr)
!        call COMBWF(WA,EIGEN_VEC14,W1)
!
!        AN=0.d0
!        do i=1, NZ*NVB*N2D*NBASJJ
!          AN=AN+abs(WA(i))**2.d0
!        end do
!        print*, 'norm of initial wavefunction WA is:', AN
!
!        deallocate(W1, stat=ierr)

        allocate(WA(NBASJJ*N2D*NVB*NZ), stat=ierr) !##
        call COMBWF(WA,EIGEN_VEC14,W1) !##
	deallocate( WA ); return

!        allocate(WI(NZINT*NRB*N2D*NBASJJ), stat=ierr)
!        WI(:)=0.d0

        end subroutine WFINIT
c=========================================================================
        subroutine GAUSS_PACKET(WP1)
        use TRANSMATS, only : TMATZ
        use QUADPOINTS, only : ZQ
        use REDUCED_MASS, only : TMS_YCZ3
        use ENGARRAYS, only : ESK0, AE0
        use COMPARAS, only : NZ, ESKinit, Zinit, DELTA, DT, NKBZinit, NE0

        implicit none

        integer :: i, j, k, MS, NS, ierr
        real*8 DZ, FAC, Zval, PI, AN
        complex*16 CI, WP1(1)
        complex*16, allocatable :: W2(:), W3(:)

        CI=(0.0d0, 1.d0)
        PI=DACOS(-1.d0)

! gauss_packet
        DZ=ZQ(2)-ZQ(1)
        ESKinit=DSQRT(2.d0*TMS_YCZ3*ESKinit)        
        FAC=DSQRT(DZ)/(PI*DELTA**2)**(0.25d0)

! gaussian wave packet on grids
        allocate(W2(NZ), stat=ierr)
        AN=0.d0
        do 10 i=1,NZ
          Zval=ZQ(i)
          W2(i)=FAC*DEXP(-((Zval-Zinit)/DELTA)**2/2)*CDEXP(-CI*ESKinit*Zval)
          AN=AN+abs(W2(i))**2.d0
          write(216,*) Zval, abs(W2(i))**2.d0
 10     continue

        call CAL_AE0(W2, WP1)

        do i=1, NE0
          write(220,*) ESK0(i)*27.2116, AE0(i)
        end do

! transfer W1 to FBR
        allocate(W3(NZ*NKBZinit),stat=ierr)
        W3(:)=0.d0
        do k=1, NZ
          MS=(k-1)*NKBZinit
        do j=1, NZ
          NS=(j-1)*NKBZinit
        do i=1, NKBZinit
          W3(MS+i)=W3(MS+i)+WP1(NS+i)*TMATZ(j,k)
        end do
        end do
        end do

        do i=1, NKBZinit*NZ
          write(222,*) i, abs(W3(i))**2.d0
        end do

        AN=0.d0
        do i=1, NKBZinit*NZ
          WP1(i)=W3(i)
          AN=AN+abs(WP1(i))**2.d0
        end do
        print*, 'Norm of WP is: ', AN
        
        deallocate(W2, W3, stat=ierr)

        end subroutine GAUSS_PACKET
c=========================================================================
        subroutine CAL_AE0(W2,W1)
        use QUADPOINTS, only : ZQ
        use ENGARRAYS, only : ESK0, AE0
        use REDUCED_MASS, only : TMS_YCZ3
        use COMPARAS, only : NE0, NZ, Zmax, Zinit, DT, Jtot, J1init, 
     .    KBZinit, NKBZinit, NSTEP0

        implicit none

        integer :: maxnodes, i, j, ierr, nshift, jrot, lrot, lrot2
        real*8 :: PI, Zstart, FAC, FAC1, FAC2, DZ, Zval, pVal, 
     .    zz, zp, zu, zup 
        real*8, allocatable :: Znodes(:)
        complex*16 CI, CAE, W1(1), W2(1)
        complex*16, allocatable :: WFZ(:), kmatZ(:), vmatZ(:)

        PI=DACOS(-1.d0)
        DZ=ZQ(2)-ZQ(1)
        CI=(0.d0,1.d0)
        maxnodes=4096

! gauss_packet in a large scope of quadrature
        Zstart=2*Zinit-Zmax
        if (Zstart.lt.0.d0) stop 'Zstart < 0.0'
        do i=1, NZ
          if (ZQ(i).gt.Zstart) goto 1000
        end do
 1000        continue
        nshift=i-1
        if (nshift.le.0) stop 'nshift <= 0'
        Zstart=ZQ(nshift)

        jrot=J1init
        lrot2=Jtot*(Jtot+1)+jrot*(jrot+1)-2*KBZinit*KBZinit
        lrot =0.5d0*(dsqrt(4.0d0*lrot2+1.0d0)-1.0d0)+0.0001d0
        FAC1=(2*DT)*(PI/(maxnodes*DZ))**2/(4*TMS_YCZ3)
        FAC2=(2*DT)*lrot2/(2*TMS_YCZ3)

        allocate(Znodes(maxnodes),stat=ierr)
        allocate(WFZ(maxnodes),kmatZ(maxnodes),vmatZ(maxnodes),stat=ierr)

! gaussian wave-packet on long grids
        WFZ(1)=0.d0
        Znodes(1)=Zstart
        do i=1, maxnodes-1
          Zval=Zstart+(i-1)*DZ
          Znodes(i+1)=Zval
          WFZ(i+1)=0.d0
          if (i+nshift.le.NZ) WFZ(i+1)=W2(i+nshift)
          
          kmatZ(i)=cdexp(CI*FAC1*i**2)
          vmatZ(i)=cdexp(CI*FAC2/(Zval**2))
        end do

c backward propagation
        call SINFT1(WFZ,maxnodes)    ! WFZ on FBR
        do 10 i=1, NSTEP0
          do 20 j=1, maxnodes-1
            WFZ(j+1)=WFZ(j+1)*kmatZ(j)
20          continue
          call SINFT1(WFZ,maxnodes)
          do 30 j=1, maxnodes-1
            WFZ(j+1)=WFZ(j+1)*vmatZ(j)
30          continue
          call SINFT1(WFZ,maxnodes)
          do 40 j=1, maxnodes-1
            WFZ(j+1)=WFZ(j+1)*kmatZ(j)
40          continue
10        continue
        call SINFT1(WFZ,maxnodes)   ! WFZ on DVR

c        do i=1, maxnodes
c          write(217,*) i, abs(WFZ(i))**2.d0
c        end do

c to calculate AE0
        FAC=DSQRT(2.0d0*TMS_YCZ3/PI)*DSQRT(DZ)
        do i=1, NE0
          pVal=DSQRT(2.d0*TMS_YCZ3*ESK0(i))

          CAE=0.d0
          do j=1, maxnodes
            Zval=Znodes(j)
            call RBES(lrot,pVal*Zval,zz,zp,zu,zup)

            CAE=CAE+FAC*(zp+CI*zz)*WFZ(j)/DSQRT(pVal)/2.d0
          end do
          AE0(i)=2*PI*ABS(CAE)**2.d0
        end do

c forward propagation
        call GPFP(nshift, jrot, maxnodes, DZ, Znodes, WFZ, W1)  ! W1 on DVR

        deallocate(Znodes, WFZ, kmatZ, vmatZ, stat=ierr)

        end subroutine CAL_AE0
c=======================================================================
        subroutine GPFP(nshift, jrot, maxnodes, DZ, Znodes, WFZ, FGAUSS)
        use QUADPOINTS, only : ZQ
        use REDUCED_MASS, only : TMS_YCZ3
        use COMPARAS, only : NZ, DT, Jtot, KBZinit, NSTEP0, KBZmin, 
     .    NKBZinit

        implicit none
        integer :: nshift, jrot, maxnodes, i, j, ik, jk, lk, ierr
        real*8 :: DZ, Znodes(maxnodes), PI, FAC, Zval, kbzVal, lamda0, 
     .    lamda1, lamda2
        real*8, allocatable :: vmatCP(:,:), kmatCP(:), WK1(:), WK2(:)
        complex*16 WFZ(maxnodes), vtemp, CI, FGAUSS(NKBZinit, NZ)
        complex*16, allocatable :: WFZ2(:,:), kmatZ(:), vmatZ(:,:,:), 
     .    CWK(:,:), CWK2(:)

        PI=DACOS(-1.d0)
        CI=(0.d0, 1.d0)

        allocate(vmatCP(NKBZinit, NKBZinit), stat=ierr)
        allocate(kmatCP(NKBZinit), stat=ierr)
        allocate(WK1(NKBZinit), WK2(NKBZinit), stat=ierr)
        vmatCP(:,:)=0.d0; kmatCP(:)=0.d0; WK1(:)=0.d0; WK2(:)=0.d0

        allocate(WFZ2(NKBZinit, maxnodes), stat=ierr)
        allocate(kmatZ(maxnodes), stat=ierr)
        allocate(vmatZ(NKBZinit, NKBZinit, maxnodes), stat=ierr)
        allocate(CWK(NKBZinit, maxnodes), stat=ierr)
        allocate(CWK2(NKBZinit), stat=ierr)

        WFZ2(:,:)=0.d0; kmatZ(:)=0.d0; vmatZ(:,:,:)=0.d0; CWK(:,:)=0.d0
        CWK2(:)=0.d0

        FAC=2*DT*(PI/(maxnodes*DZ))**2/(2*TMS_YCZ3)/2

        if (NSTEP0.ne.0) call SINFT1(WFZ, maxnodes)  ! WFZ on FBR


! GP for init K block on FBR if NSTEP0 != 0
! GP for init K block on DVR if NSTEP0 = 0
        do i = 1, maxnodes-1
          WFZ2(KBZinit-KBZmin+1,i+1)=WFZ(i+1)
c          write(218,*) i, abs(WFZ2(KBZinit-KBZmin+1,i+1))**2.d0
        end do
        if (NSTEP0.eq.0) goto 1000

! construct kinetic energy and centrifugal potential
        do 10 i=1, maxnodes-1

          Zval=Znodes(i+1)
          kmatZ(i)=CDEXP(-CI*FAC*i**2)
          
          vmatCP(:,:) = 0.d0
          do 20 ik=1, NKBZinit
            kbzVal=ik+KBZmin-1
            lamda0=Jtot*(Jtot+1)+jrot*(jrot+1)-2*kbzVal**2
            vmatCP(ik,ik)=lamda0/(2*TMS_YCZ3*Zval**2)
20          end do
          do 30 ik=1,NKBZinit-1
            kbzVal=ik+KBZmin-1
            lamda1=DSQRT(Jtot*(Jtot+1.d0)-kbzVal*(kbzVal+1.d0))
            lamda2=DSQRT(jrot*(jrot+1.d0)-kbzVal*(kbzVal+1.d0))
            lamda0=-lamda1*lamda2/(2*TMS_YCZ3*Zval**2)

            if (kbzVal.eq.0) lamda0=DSQRT(2.d0)*lamda0
            vmatCP(ik,ik+1)=lamda0
            vmatCP(ik+1,ik)=lamda0
30          end do

          call RS(NKBZinit,NKBZinit,vmatCP,kmatCP,1,vmatCP,WK1,WK2,ierr)
          
          do 40 ik=1, NKBZinit
            do 40 jk=1, NKBZinit
              vtemp=0.d0
              
              do 50 lk=1, NKBZinit
                vtemp=vtemp+vmatCP(jk,lk)*vmatCP(ik,lk)*
     .            CDEXP(-CI*(2.d0*DT)*kmatCP(lk))
 50              continue
              vmatZ(jk,ik,i)=vtemp
 40          continue

10        continue

! forward popagation
        do 80 i=1, NSTEP0

          do 90 j=1, maxnodes-1
          do 100 ik=1, NKBZinit
            WFZ2(ik,j+1)=WFZ2(ik,j+1)*kmatZ(j)
 100          continue
 90          continue

          call SINFT(WFZ2, maxnodes, NKBZinit, CWK2)

          do 120 j=1, maxnodes-1
            do 130 ik=1, NKBZinit
              vtemp=0.d0
              do 140 jk=1, NKBZinit
                vtemp=vtemp+WFZ2(jk,j+1)*vmatZ(jk,ik,j)
 140              continue
              CWK(ik,j+1)=vtemp
 130            continue
 120          continue

          call SINFT(CWK, maxnodes, NKBZinit, CWK2)

          do 190 j=1, maxnodes-1
          do 200 ik=1, NKBZinit
            WFZ2(ik,j+1)=CWK(ik,j+1)*kmatZ(j)
 200          continue
 190          continue

 80        continue


        call SINFT(WFZ2, maxnodes, NKBZinit, CWK2)   ! WFZ2 on DVR
c        do i=1, maxnodes
c        do ik=1, NKBZinit
c            write(219,*) ik, i, abs(WFZ2(ik, i))**2.d0
c        end do
c        end do

 1000        continue

! FGAUSS on DVR
        FGAUSS(:,:)=0.d0
        do i=1, maxnodes
        do ik=1, NKBZinit
          if (i+nshift.le.NZ) then 
            FGAUSS(ik,i+nshift)=WFZ2(ik,i+1)
c            write(221,*) ik, i, FGAUSS(ik, i+nshift)
          end if
        end do
        end do

        deallocate(WFZ2, CWK, CWK2, vmatCP, kmatCP, vmatZ, kmatZ, 
     &    WK1, WK2, stat=ierr)

        end subroutine GPFP
c=======================================================================
        subroutine BS_YCZ3
        use WFARRAYS, only : EIGEN_VEC14
        use ENGARRAYS, only : EIGEN_ENG14
        use COMPARAS, only : NVB, N2D, NBASJJ14
        implicit none
        integer ierr, ms, i, j, k


        call COUPLE14       ! (j12, ja2, kb2) with j1=J1init

        call VPOT14         ! V(R, AU, A2, B2) with Z=Zinf, RU=RU0, A1=0.0 
                            ! and B1=0.0
        allocate(EIGEN_VEC14(NVB*N2D*NBASJJ14), stat=ierr)
        call BOUND14(EIGEN_ENG14, EIGEN_VEC14)

        print*, 'EIGEN_ENG14=', EIGEN_ENG14, EIGEN_ENG14*27.2116d0

        return

        end subroutine BS_YCZ3
c===================================================================
        subroutine COMBWF(WF15,WF14,FGAUSS)
        use ROTBASIS 
        use COMPARAS, only : NZ, NZINT, NRB, NVB, N2D, NBASJJ, IPOT,
     &    NBASJJ14, NKBZinit, J1init, KBZmin, KBZmax, NVBRU, NVBAU  
     &        ,NA2node, NB2node  ,KB2max ,NKB20
        use ENGARRAYS, only : EIGEN_ENG14
        use TRANSMATS, only : TMAT2D
        use SMALL_ARRAYS, only : TMATRU, TMATAU
	use ENGARRAYS, only : VMAT14
	use SMALL_ARRAYS, only : BROTR, BIA_CZ3, BIC_CZ3

        implicit none
        integer i2, i, j, k, l, ms, ns, ierr  ,Iparity   
	integer j12Val_1, j2Val_1, ik2_1, kb2Val_1 ,ia2, ib2, m2Val, m2Max

        integer j12Val, j2Val, ik2, kb2Val, kbzVal, ikbz, ibasJJ ,jbasJJ
        real*8 AN, dsqrt2     ,krot
!        real*8, allocatable :: WK(:), WK2(:), WK3(:)
!        real*4, allocatable :: WTMP(:,:), WTMP2(:)
       
        complex*8 WF15(NBASJJ,N2D,NVB,NZ)
        complex*16 FGAUSS(NKBZinit,NZ)
        real*8     WF14(NBASJJ14,N2D,NVB)  
        real*8 :: norm

!        real*8 WF2(NBASJJ14,N2D,NVB), WF3(0:KB2max,N2D,NVB)
!        real*8 WK(0:KB2max)

        dsqrt2=dsqrt(2.d0)

	if(IPOT == 1) then
        open(100, file='WF_CH4',form='binary',
     &                 status='UNKNOWN',ACTION='WRITE')
        WRITE(100) EIGEN_ENG14

        norm=0.D0
        ms=0
        do k=1, NVB
        do j=1, N2D
          ms=ms+1
        do i=1, NBASJJ14
	   WRITE(100) WF14(i,j,k)
          norm=norm + abs(WF14(i,j,k))**2.d0
	  if(abs(WF14(i,j,k)).GE.1.D-4) write(1414,*) WF14(i,j,k)
        end do
        end do
        end do
        close(100)

        print*,'WF_CH4 norm=',norm

	end if

!============================
!        sig=0.D0
!         WF3=0.D0
!          do k=1,NVB
!            do j=1,N2D
!              do i=1,NBASJJ14
!                 ik2=IDXJJ14(3,i);kb2Val=KB2(ik2)
!                WF3(kb2Val,j,k)=WF3(Kb2Val,j,k)+DABS(WF14(i,j,k))**2
!                sig=sig+abs(WF14(i,j,k))**2.d0
!	      write(88,*)kb2Val,WF3(kb2Val,j,k)
!              end do
!            end do
!          end do
!         print*,'norm of initial bound state' ,sig
!
!         sig=0.D0
!         WK=0.D0
!         do k=1,NVB
!         do j=1,N2D
!         do ik2=1,NKB20 
!	   kb2Val=KB2(ik2)
!            sig=sig+WF3(kb2Val,j,k)
!          WK(kb2Val)=WK(kb2Val)+WF3(kb2Val,j,k)
!         enddo
!         enddo
!         enddo
!          print*,'norm of initial bound state=' ,sig
!
!         sig=0.D0
!         do ik2=1,NKB20 
!	   kb2Val=KB2(ik2)
!           write(8,*)Kb2Val,WK(Kb2Val)
!           sig=sig+WK(Kb2Val)
!         enddo
!         print*,'norm of initial bound state==' ,sig


        end subroutine COMBWF
c===================================================================

