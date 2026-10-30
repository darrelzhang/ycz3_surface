        subroutine AUBASIS
        use COMPARAS, only : NAU,AUmin,AUmax,NVBAU,NVBRU,RUmin,N2D
        use SMALL_ARRAYS, only : EIG2D,Vref2D
        use QUADPOINTS, only : AUQ, RUQ
        use TRANSMATS, only : TMAT2D
        implicit none
        integer ierr, i, j
        real(8), allocatable :: VPOT_AU(:),AUQ0(:),EIGAU(:),VrefAU(:)
        real*8, allocatable :: TMATAU(:,:)        !###

        allocate(AUQ0(NAU),VPOT_AU(NAU),EIGAU(NAU),stat=ierr)
        allocate(AUQ(NVBAU),VrefAU(NVBAU),stat=ierr)
        allocate(TMATAU(NVBAU,NVBAU),stat=ierr)
	if(NVBAU /= 1) then
          call AUPOT1D(NAU,AUmin,AUmax,AUQ0,VPOT_AU)
	call AUVIB1D(AUQ0,VPOT_AU,EIGAU,TMATAU,VrefAU)
	else
	  NAU=1;  AUQ0(1)=AUmin;  AUQ(1)=AUmin;
	  VPOT_AU(1)=0.d0
	  VrefAU(1)=0.d0
	  EIGAU(1)=0.d0
	  TMATAU(:,:)=1.d0
	end if

!        print*, "Parameters of AU basis"
!        print*, "============================="
!        print*
!        print*, "Eigenenergies"
!        do i=1, NVBAU
!          write(*,*) i, EIGAU(i)
!        end do
!        write(*,*)
!        print*, "PODVR Points and Reference Potential for AU"
!        !open(10,file="aau",status="UNKNOWN")
!        do i=1, NVBAU
!        write(*,'(I5,3F15.10)') i, AUQ(i), VrefAU(i), AUQ(i)*180.d0/dacos(-1.d0)
!          !write(10,*) AUQ(i)
!        end do
!	!close(10)
!        write(*,*)
!        print*, "============================="
!        print*

c        do j=1, NVBAU
c        do i=1, NVBAU
c          write(204,*) i, j, TMATAU(i,j)
c        end do
c        end do

        ! we need to set up 2D related information, while NVBRU is 1
        allocate(RUQ(NVBRU),stat=ierr)
        RUQ(NVBRU)=RUmin
        !open(10,file="rru",status="UNKNOWN")
        do i=1, NVBAU
          !write(10,*) RUQ(i)
        end do
        !close(10)
        allocate(TMAT2D(NVBAU,NVBAU),stat=ierr)
        TMAT2D=TMATAU
        allocate(EIG2D(NVBAU),stat=ierr)
        do i=1, NVBAU; EIG2D(i)=EIGAU(i); end do
        allocate(Vref2D(NVBAU),stat=ierr)
        Vref2D=VrefAU
!	 Vref2D(1)=Vref2D(2)
!	 Vref2D(2)=Vref2D(3)
        end subroutine AUBASIS
c------------------------------------------------------------------
c This file contains subroutines to resolve bound state of umbrella
! vibration mode (au). Codes from Dr. Yang Minghui's design.
c Some minor changes are made for readability improvement.
c------------------------------------------------------------------
        SUBROUTINE AUPOT1D(NR,R1,R2,RQ0,VPOT)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION VPOT(1),RQ0(1)

        DR=(R2-R1)/(NR+1)
        DO I=1,NR
          RQ0(I)=R1+DR*I
        END DO

        call REFPOT4(NR,RQ0,VPOT)
!	 DO i=1,14
!	    VPOT(i)=VPOT(15)
!	 ENDDO 

        open(100, file='AUPOT.1D',form='formatted',status='unknown')
        do i=1,NR
          write(100,*) RQ0(i)*180.d0/dacos(-1.d0), VPOT(i)*27.2116d0
        end do
        close(100)
 
        RETURN
        END
c==================================================================
        subroutine AUVIB1D(AUQ0,VPOT_AU,EIGAU,TMATAU,VrefAU)
        use QUADPOINTS, only :  AUQ
        use REDUCED_MASS, only : TMS_3Z,RMS_CZ3
        use COMPARAS, only : NAU, NVBAU, AM, RUmin

        implicit none

        real*8 WK1(NAU), WK2(NAU) 
        real*8 AUQ0(1),VPOT_AU(1),EIGAU(1),VrefAU(1)
        real*8 TMATAU(NVBAU,NVBAU)        !###
        real*8, allocatable :: HMATAU(:,:)

        integer i, j, IERR
        real*8 :: rtmp

!        allocate(EIGAU(NVBAU),VrefAU(NVBAU),TMATAU(NVBAU,NVBAU),
!     .    AUQ(NVBAU),stat=ierr)
        allocate(HMATAU(NAU, NAU), stat=ierr)

! to express kinetic operator on sine DVR 
        call KCZ3DVR(AUQ0,HMATAU)
! to add reference potential
        do i = 1, NAU
          HMATAU(i,i)=HMATAU(i,i)+VPOT_AU(i)
        end do

! to diagonalize HMATAU
        call RS(NAU,NAU,HMATAU,EIGAU,1,HMATAU,WK1,WK2,IERR)

! to de-phase the eigenvectors
        call WFPHASE(NAU,HMATAU)

! to create PODVR for AU
        call AUDVR(AUQ0,TMATAU,HMATAU)

! to calculate the Vref on PODVR points
        call REFPOT4(NVBAU,AUQ,VrefAU)
!	 VrefAU(1)=VrefAU(2)
!	 VrefAU(2)=VrefAU(3)

        open(100, file='AUPODVR.POT',form='formatted',status='unknown')
!	 open(100, file='AUPODVR.POT',form='unformatted',status='unknown')
!	 write(100) (VrefAU(i),i=1,NVBAU)
        do i=1,NVBAU
          write(100,*) AUQ(i)*180.d0/dacos(-1.d0), VrefAU(i)*27.2116d0
        end do
        close(100)
!	 open(121,file='umbrewf15',form='formatted',status='unknown')
!	 do i=1,NVBAU
!	    write(121,*)(HMATAU(j,i),j=1,NAU)
!	 enddo
!	 close(121)

        end subroutine AUVIB1D
c===========================================================================
        subroutine KCZ3DVR(AUQ0,HMATAU)
        use COMPARAS, only : NAU, AUmin, AUmax, RUmin
        use REDUCED_MASS, only : RMS_CZ3, TMS_3Z
        implicit none

        integer jbas, ibas, inode, i, j, ierr
        real*8 AUQ0(1), HMATAU(NAU,NAU)
        real*8, allocatable :: TMAT(:,:), HAR(:,:)

        real*8 const1, const2, const3, const4, DAU, tmp, f1,
     .    tmp1, tmp2, au


        allocate(TMAT(NAU,NAU),HAR(NAU,NAU),stat=ierr)

! to set up TMAT (sin transformation matrix)
        const1=dsqrt(2.0d0/(NAU+1))
        const2=dacos(-1.d0)/(NAU+1)
        do 10 i = 1, NAU
          do 20 j = 1, i
            TMAT(i,j)=const1*sin(i*j*const2)
            TMAT(j,i)=TMAT(i,j)
 20       continue
 10     continue

! to set up matrix of kinetic operator on sine basis Kij=<j|K|i>
        DAU=AUQ0(2)-AUQ0(1)
        AUmin=AUQ0(1)-DAU
        AUmax=AUQ0(NAU)+DAU
        const1=2.0d0/(AUmax-AUmin)
        const2=DACOS(-1.0d0)/(AUmax-AUmin)
        const3=2.0d0/(NAU+1)
        const4=DACOS(-1.d0)/(NAU+1)

        do 50 jbas=1,NAU
        do 60 ibas=1,NAU
          tmp=0.d0
            
        do 70 inode=1,NAU
          au=AUQ0(inode)

          f1=dsin(ibas*inode*const4)
          tmp1=0.5*(jbas*const2)**2
     .      *(dcos(au)**2/TMS_3Z+dsin(au)**2/RMS_CZ3)
     .      *dsin(jbas*inode*const4)
          tmp2=(jbas*const2)*dsin(au)*dcos(au)
     .      *(1.0d0/TMS_3Z-1.d0/RMS_CZ3)*dcos(jbas*inode*const4)
          tmp=tmp+f1*(tmp1+tmp2)
 70     continue
          tmp=tmp*DAU*const1/RUmin**2
          HMATAU(ibas,jbas)=tmp
 60     continue
 50     continue

c representation transformation (sine FBR -> sine DVR)
        call DGEMM('N','N',NAU,NAU,NAU,
     &                    1.0d+0,HMATAU,  NAU,
     &                           TMAT,  NAU,
     &                    0.0d+0,HAR,  NAU)
        call DGEMM('N','N',NAU,NAU,NAU,
     &                    1.0d+0,TMAT,  NAU,
     &                           HAR,  NAU,
     &                    0.0d+0,HMATAU,  NAU)

        end subroutine KCZ3DVR
c===========================================================================
        subroutine WFPHASE(NAU,HMATAU)
        implicit none
        
        integer ibas, idvr, NAU
        real*8 tmp, HMATAU(NAU,1)

        do 10 ibas=1,NAU
          tmp=0.0d0
          do 20 idvr=1,NAU
            if (abs(HMATAU(idvr,ibas)).gt.abs(tmp)) tmp=HMATAU(idvr,ibas)
 20       continue
          if (tmp.lt.0.0d0) then
            do 40 idvr=1,NAU
              HMATAU(idvr,ibas)=-HMATAU(idvr,ibas)
 40         continue
          end if
 10     continue

        end subroutine WFPHASE
c==========================================================================
        subroutine AUDVR(AUQ0,TMATAU,HMATAU)
        use QUADPOINTS, only :  AUQ
        use COMPARAS, only : NAU, NVBAU

        implicit none
        
        integer I, J, IDVR, IERR
        real*8 AUQ0(NAU),HMATAU(NAU,NAU)
        real*8 TMATAU(NVBAU,NVBAU)  !###
        real*8 TMAT(NVBAU,NVBAU)
        real*8 WK1(NVBAU), WK2(NVBAU), RX

        DO 20 I=1,NVBAU
        DO 20 J=1,NVBAU
          RX=0
          DO 10 IDVR=1,NAU
             RX=RX+AUQ0(IDVR)*HMATAU(IDVR,I)*HMATAU(IDVR,J)
 10       CONTINUE
          TMAT(J,I)=RX
 20     CONTINUE

        CALL RS(NVBAU,NVBAU,TMAT,AUQ,1,TMAT,WK1,WK2,IERR)

        TMATAU=TMAT
        
        end subroutine AUDVR

