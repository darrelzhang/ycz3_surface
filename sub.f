c================================================================
        subroutine ZBASIS
        use COMPARAS, only : NZ, NZINT, Zmin, Zmax
        use QUADPOINTS, only : ZQ
        use TRANSMATS, only : TMATZ, TMATIZ
        use SMALL_ARRAYS, only : ETASY, ETINT, BROTZ
        use REDUCED_MASS, only : TMS_YCZ3
        implicit none
        integer ierr, i, j

        allocate(ZQ(NZ), stat=ierr)
        allocate(ETASY(NZ),ETINT(NZINT),stat=ierr)
        allocate(TMATZ(NZ,NZ),TMATIZ(NZINT,NZINT),stat=ierr)

        call QUAD1D(NZ,Zmin,Zmax,ZQ)
        call SINBASIS(NZ,Zmin,Zmax,TMS_YCZ3,TMATZ,ETASY)
        call SINBASIS(NZINT,Zmin,ZQ(NZINT+1),TMS_YCZ3,TMATIZ,ETINT)

!        print*, "The parameters of Z basis"
!        print*, "====================="
!        print*, "Eigenenergies and Points"
!        print*
!	!open(10,file="zz",status="UNKNOWN")
!        do i=1, NZ
!          write(*,*) i, ETASY(i), ZQ(i)
!	  !write(10,*) ZQ(i)
!        end do
!	!close(10)
!        write(*,*)
!        do i=1, NZINT
!          write(*,*) i, ETINT(i), ZQ(i)
!        end do
!        print*, "====================="
!        print* 

c        do j=1, NZ
c        do i=1, NZ
c          write(200,*) i, j, TMATZ(i,j)
c        end do
c        end do
c        close(200)

c        do j=1, NZINT
c        do i=1, NZINT
c          write(201,*) i, j, TMATIZ(i,j)
c        end do
c        end do
c        close(201)

        allocate(BROTZ(NZ), stat=ierr)
        do i=1, NZ
          BROTZ(i)=1.d0/(2*TMS_YCZ3*ZQ(i)**2.d0)
        end do

        end subroutine ZBASIS
c================================================================
        subroutine R1BASIS
        use COMPARAS, only : NVB, NRB, NR, Rmin, Rmax, IPOT
        use QUADPOINTS, only : RQ, RQI
        use TRANSMATS, only : TMATR, TMATIR
        use SMALL_ARRAYS, only : RQ0, VPOT, EIG0, VrefCH, HAR  
        implicit none
        integer ierr, i, j
        real*8 CM,rtmp

        allocate(RQ0(NR),VPOT(NR), stat=ierr)
        call POT1D(NR,Rmin,Rmax,RQ0,VPOT)

        allocate(EIG0(NRB),VrefCH(NRB+NVB), stat=ierr)
        allocate(HAR(NR,NR),TMATR(NVB,NVB),TMATIR(NRB,NRB),stat=ierr)
        allocate(RQ(NVB), RQI(NRB), stat=ierr)
        if (NVB.eq.1) then
          TMATR(1,1)=1.d0
          TMATIR(1,1)=1.d0
          EIG0(1)=0.d0
          VrefCH(NRB+1)=0.d0
          RQ(1)=RQ0(1)
          RQI(1)=RQ0(1)

          return
        end if
        call VIB1D

        call REFPOT2(NRB,RQI,VrefCH)
        call REFPOT2(NVB,RQ,VrefCH(NRB+1))

!	if(IPOT == 1) then
          open(100,file='R1PODVR.POT',form='formatted',status='unknown')
!	  open(100,file='R1PODVR.POT',form='unformatted',status='unknown')
!	   write(100) (VrefCH(i),i=1,NRB+NVB)
          do i=1,NRB
            write(100,*) RQI(i), VrefCH(i)*27.2116d0
          end do
          do i=1,NVB
            write(100,*) RQ(i), VrefCH(NRB+i)*27.2116d0
          end do
          close(100)
!	end if


!        print*, "The parameters of R1 basis"
!        print*, "=========================="
!        print*, "Eigenenergies"
!        print*
!
!        CM=219475d0
!        do i=1, NRB
!          write(*,*) i, EIG0(i)*CM, (EIG0(i)-EIG0(1))*CM
!        end do
!        write(*,*)
!
!        print*, "Reference Potential"
!        do i=1, NRB+NVB
!          write(*,*) i, VrefCH(i)
!        end do
!        write(*,*)
!
!        print*, "PODVR Points in Interaction Area"
!	!open(10,file="rr",status="UNKNOWN")
!        do i=1, NRB
!          write(*,*) i, RQI(i)
!          !write(10,*) RQI(i)
!        end do
!	!close(10)
!        write(*,*)
!        print*, "PODVR Points in Asymptotic Area"
!	!open(10,file="rr.int",status="UNKNOWN")
!        do i=1, NVB
!          write(*,*) i, RQ(i)
!	  !write(10,*) RQ(i)
!        end do
!	!close(10)
!        write(*,*)
!
!
!        print*, "=========================="
!        write(*,*)
!
!c        do j=1, NVB
!c        do i=1, NVB
!c           write(202,*) i, j, TMATR(i,j)
!c        end do
!c        end do
!c        close(202)
!
!c        do j=1, NRB
!c        do i=1, NRB
!c          write(203,*) i, j, TMATIR(i,j)
!c        end do
!c        end do        
!c        close(203)

        end subroutine R1BASIS
c================================================================
        SUBROUTINE QUAD1D(NZ,Z1,Z2,ZQ)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION ZQ(1)

        DZ=(Z2-Z1)/(NZ+1)
        DO I=1,NZ
          ZQ(I)=Z1+I*DZ
        END DO

        RETURN
        END
c================================================================
        SUBROUTINE SINBASIS(NB1,R1,R2,RMS,TMAT,ET)
        IMPLICIT REAL*8 (A-H,O-Z)
        REAL*4 TMAT(NB1,1)
        DIMENSION ET(1),RQ1(1)

        PI=DACOS(-1.d0)
        RL=R2-R1
        DO I=1,NB1
           ET(I)=(I*PI/RL)**2/2.D0/RMS
        END DO

        FAC=DSQRT(2.D0/(NB1+1))
        DO I=1,NB1
        DO J=1,NB1
           TMAT(J,I)=FAC*DSIN(I*J*PI/(NB1+1))
        END DO
        END DO

        RETURN
        END
c=================================================================
        SUBROUTINE POT1D(NR,R1,R2,RQ0,VPOT)
	use COMPARAS, only : IPOT
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION VPOT(1),RQ0(1)


        DR=(R2-R1)/(NR+1)
        DO I=1,NR
          RQ0(I)=R1+DR*I
        END DO

        if (NR.eq.1) RQ0(1)=R1

        call REFPOT2(NR,RQ0,VPOT)

	!if(IPOT /= 1) return

        open(100, file='R1POT.1D',form='formatted',status='unknown')
        do i=1,NR
          write(100,*) RQ0(i), VPOT(i)*27.2116d0
        end do
        close(100)

        RETURN
        END 
c=================================================================
        SUBROUTINE VIB1D
c------------------------------------------------
c Vref contains the reference for NVB and NRB
c------------------------------------------------
        use QUADPOINTS, only : RQ, RQI
        use TRANSMATS, only : TMATR, TMATIR
        use REDUCED_MASS, only : RMS_YCZ3
        use COMPARAS, only : NR, NVB, NRB
        use SMALL_ARRAYS, only : BROTR, RQ0,VPOT,HAR,EIG0,VrefCH

        implicit real*8 (A-H,O-Z)
        dimension EIGj(NR),HARj(NR,NR)
        dimension DTMAT(NRB*NRB)

        CM=219474.6D0

c solve vibrational energy
        CALL VIBSOL(NR,NR,RQ0,RMS_YCZ3,HAR,VPOT,EIGj)
	 print*, "C--H 1D ref eigenenergy:"
        DO I=1,NRB
          EIG0(I)=EIGj(I)
	  print*,EIG0(I)*CM
        END DO

c transfer HAR to sine basis and save in PARAMETER22
c        FAC=SQRT(2.D0/(NR+1))
c        PI=ACOS(-1.D0)
c        DR=RQ0(2)-RQ0(1)
c        WRITE(51,*) NR,RQ0(1)-DR,RQ0(NR)+DR
c        DO I=1,NRB
c        DO J=1,NR
c           AN=0.D0
c        DO K=1,NR
c           AN=AN+HAR(K,I)*FAC*SIN(J*K*PI/(NR+1))
c        END DO
c           WRITE(51,*) AN      !5151
c        END DO
c        END DO

! PODVR for NRB
        CALL PODVR(NR,NRB,HAR,DTMAT,RQ0,RQI,RMS_YCZ3,EIG0,VrefCH)
        DO 505 I=1,NRB
        DO 505 J=1,NRB
505     TMATIR(J,I)=DTMAT((I-1)*NRB+J)

c generate PODVR for NVB
        IF(NVB.GE.NRB) RETURN
        CALL PODVR(NR,NVB,HAR,DTMAT,RQ0,RQ,RMS_YCZ3,EIG0,VrefCH(NRB+1))
        DO 506 I=1,NVB
        DO 506 J=1,NVB
506     TMATR(J,I)=DTMAT((I-1)*NVB+J)

        allocate(BROTR(NRB+NVB), stat=ierr)
        do i=1, NRB
          BROTR(i)=1.0d0/(2*RMS_YCZ3*RQI(i)**2)
        end do

        do i=1, NVB
          BROTR(NRB+i)=1.d0/(2*RMS_YCZ3*RQ(i)**2)
        end do

        RETURN
        END
c==================================================================
        SUBROUTINE VIBSOL(NQ1max,NQ1,RQ1,TMS,HAR,Vpot,EIG)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION  RQ1(NQ1),HAR(NQ1max,NQ1),Vpot(1),EIG(NQ1),
     $  WK1(NQ1),WK2(NQ1)

c solve eigenvalue problem for one dimensional potential.
        RL=RQ1(NQ1)-RQ1(1)+2.D0*(RQ1(2)-RQ1(1))
        NDVR=NQ1+1
        DO 20 I=1,NQ1
           HAR(I,I)=Vpot(I)+DVRKE1(I,NDVR,RL,TMS)
        DO 20 J=1,I-1
           HAR(J,I)=DVRKE2(I,J,NDVR,RL,TMS)
20        HAR(I,J)=HAR(J,I)


        CALL RS(NQ1max,NQ1,HAR,EIG,1,HAR,WK1,WK2,IERR)

        IF(IERR.NE.0) STOP 'DVRBASIS1'

        DO I=1,NQ1
          AN=0.D0
          DO J=1,NQ1
            IF(ABS(HAR(J,I)) .GT. AN) THEN
              AN=ABS(HAR(J,I))
              ID=J
            END IF
          END DO
          
          IF(HAR(ID,I).LT.0.D0) THEN
          DO J=1,NQ1
            HAR(J,I)=-HAR(J,I)
          END DO
          END IF
        END DO
          
        RETURN
        END
c------------------------------------------------------------------------
        FUNCTION DVRKE1(I,N,RL,TMS)
        IMPLICIT REAL*8 (A-H,O-Z)
        DATA PI/3.14159265358979D0/

        DVRKE1=PI*PI/4.D0/TMS/RL/RL*((2*N*N+1.D0)/3
     $  -1.D0/DSIN(PI*I/N)**2)
        RETURN
        END
C------------------------------------------------------------------------
        FUNCTION DVRKE2(I,J,N,RL,TMS)
        IMPLICIT REAL*8 (A-H,O-Z)
        DATA PI/3.14159265358979D0/

        DVRKE2=PI*PI/4.D0/TMS/RL/RL*(-1)**(I-J)
     $*(1.D0/DSIN(PI*(I-J)/2.D0/N)**2-1.D0/DSIN(PI*(I+J)/2.D0/N)**2)
        RETURN
        END
c======================================================================
        SUBROUTINE PODVR(NR,NVB,WAVE,TMAT,RQ1,DVR,RMS,EVIB,Vref)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION WAVE(NR,1),TMAT(NVB,1),RQ1(1),DVR(1),
     $  WK1(NR),WK2(NR),Vref(1),EVIB(1)

c construct X-matrix.
        DO 20 I=1,NVB
        DO 20 J=1,NVB
           TMAT(J,I)=0.D0
        DO K=1,NR
        TMAT(J,I)=TMAT(J,I)+WAVE(K,J)*WAVE(K,I)*RQ1(K)
        END DO
20        CONTINUE
c20      TMAT(I,J)=TMAT(J,I)

        CALL RS(NVB,NVB,TMAT,DVR,1,TMAT,WK1,WK2,IERR)
        IF(IERR.NE.0) STOP 'PODVR'
        DO I=1,NVB
           IF(TMAT(1,I).LT.0.D0) THEN
           DO J=1,NVB
           TMAT(J,I)=-TMAT(J,I)
           END DO
           END IF
        END DO

c calculate the reference potentail on DVR point
        DO I=1,NVB
          WK2(I)=0.D0
        END DO

        DO N0=1,NVB

        FAC=DSQRT(2.D0/(NR+1))
        PI=DACOS(-1.D0)
        DO I=1,NR
          WK1(I)=0.D0
          DO J=1,NR
            WK1(I)=WK1(I)+WAVE(J,N0)*FAC*DSIN(I*J*PI/(NR+1))
          END DO
        END DO
        RL=RQ1(NR)-RQ1(1)+2.D0*(RQ1(2)-RQ1(1))
        FAC=DSQRT(2.D0/RL)
        R1=RQ1(1)-(RQ1(2)-RQ1(1))

        DO I=1,NVB
          WF=0.D0
          AKE=0.D0
          DO J=1,NR
            WF=WF+WK1(J)*FAC*DSIN(J*PI*(DVR(I)-R1)/RL)
            AKE=AKE+WK1(J)*FAC*DSIN(J*PI*(DVR(I)-R1)/RL)
     $          *(J*PI/RL)**2/(2.D0*RMS)
          END DO
          IF(ABS(WF).GT.WK2(I)) THEN
            WK2(I)=ABS(WF)
            Vref(I)=EVIB(N0)-AKE/WF
          END IF
        END DO

        END DO

        RETURN
        END
c=========================================================================
        subroutine FFLUX
        use COMPARAS, only : NR, NRB, RFLUX, Rmin, Rmax
        use SMALL_ARRAYS, only : WVIBR, HAR
        implicit none
        real*8 PI, FAC, RL, CSIN(NR,NR)
        integer I,J,K, ierr

        allocate(WVIBR(NRB,2), stat=ierr)

c construct VIB_WF matrix in sine basis from sin_DVR basis.
        PI=ACOS(-1.D0)
        FAC=SQRT(2.D0/(NR+1))
        DO 100 J=1,NR
        DO 100 I=1,NR
        CSIN(I,J)=0.D0
        DO 100 K=1,NR
100     CSIN(I,J)=CSIN(I,J)+HAR(K,J)*FAC*SIN(I*K*PI/(NR+1))

c calculate VIB and VIB_der value at NFR point.
        RL=Rmax-Rmin
        RFLUX=RFLUX-Rmin
        FAC=SQRT(2.D0/RL)
        DO K=1,NRB
           WVIBR(K,1)=0.D0
           WVIBR(K,2)=0.D0
           DO J=1,NR
             WVIBR(K,1)=WVIBR(K,1)+CSIN(J,K)*FAC*SIN(J*PI*RFLUX/RL)
             WVIBR(K,2)=WVIBR(K,2)+CSIN(J,K)*FAC*J*PI/RL*COS(J*PI*RFLUX/RL)
           END DO
        END DO

c        do I=1, NRB
c          write(214,*) I, WVIBR(I,1)
c          write(215,*) I, WVIBR(I,2)
c        end do

        end subroutine FFLUX
c=========================================================================
        subroutine WFABS
        use QUADPOINTS, only : ZQ, RQI
        use ENGARRAYS, only : FZABS, FRABS
        use COMPARAS, only : NZ, NRB, ZABS, CZABS, PZABS, RABS, CRABS,
     .    PRABS, NZABS, NRABS, DT, Zmin, Zmax, Rmin, Rmax
        implicit none
        integer I, ierr

        allocate(FZABS(NZ),FRABS(NRB), stat=ierr)

        call ABS_POT(PZABS,NZ,ZQ,CZABS,ZABS,NZABS,FZABS,DT,1)
        call ABS_POT(PRABS,NRB,RQI,CRABS,RABS,NRABS,FRABS,DT,1)

        FRABS(18)=0.997322535522763
        FRABS(19)=0.987453765516319
        FRABS(20)=0.973007533788959
        FRABS(21)=0.954644697538218
        FRABS(22)=0.932688067321339
        FRABS(23)=0.907390159989158
        FRABS(24)=0.879028301893798
        FRABS(25)=0.847968217285209

        print*, "The parameters for absorption potentials"
        print*, "=========================="
        PRINT*,'NZABS=',NZABS,'ZABS=',ZQ(NZABS)
        PRINT*,'NRABS=',NRABS,'RABS=',RQI(NRABS)
        print*
        print*, "=========================="
        print*

c        do I=1, NZ
c          write(212,*) I, ZQ(I), FZABS(I)
c        end do
c
c        do I=1, NRB
c          write(213,*) I, RQI(I), FRABS(I)
c        end do

        end subroutine WFABS
c=========================================================================
        SUBROUTINE ABS_POT(PABS,NQ1,RQ1,CABS,RABS,NABS,FABS,DT,ID)
c absorption potentail.

        IMPLICIT REAL*8 (A-H,O-Z)
        REAL*4 FABS(1)
        DIMENSION RQ1(NQ1)

        IF(ID.EQ.1) THEN
        RL=RQ1(NQ1)+RQ1(2)-RQ1(1)
c        RL=R2-R1
        RABS=RL-RABS
        DO I=1,NQ1
          IF(RQ1(I).GE.RABS) THEN
            FABS(I)=EXP(-CABS*((RQ1(I)-RABS)/(RL-RABS))**PABS*DT)
          ELSE
            NABS=I
            FABS(I)=1.D0
          END IF
c          WRITE(7,*) RQ1(I),FABS(I)
        END DO

        ELSE IF(ID.EQ.-1) THEN
          R1=RQ1(1)-(RQ1(2)-RQ1(1))
          RABS=RABS+R1
        DO I=1,NQ1
          IF(RQ1(I).GE.RABS) THEN
            FABS(I)=1.D0
          ELSE
            FABS(I)=EXP(-CABS*((RABS-RQ1(I))/(RABS-R1))**PABS*DT)
            NABS=I
          END IF
c          WRITE(7,*) RQ1(I),FABS(I)
        END DO

        ELSE
          PRINT*,'WRONG OPTION IN WFABS'
          STOP
        END IF

        RETURN
        END
c=========================================================================
        subroutine CPLXENG
        use QUADPOINTS, only : ZQ, RQ, RQI
        use ENGARRAYS, only : CENT1, CENT2
        use COMPARAS, only : AM, NZ, NZINT, NVB, NRB, N2D, DT
        use SMALL_ARRAYS, only : ETASY, ETINT, EIG0, EIG2D
        implicit none
        integer i, j, k, ms, ierr

        complex*16 CI

        CI=(0.d0, 1.d0)

        allocate(CENT1(NZ*NVB*N2D), CENT2(NZINT*(NRB-NVB)*N2D), stat=ierr)
        ms=0
        do i=1, NZ
          do j=1, NVB
            do k=1, N2D
              ms=ms+1
              CENT1(ms)=cdexp(-CI*(ETASY(i)+EIG0(j)+EIG2D(k))*DT/2.d0)
            end do
          end do
        end do

        ms=0
        do i=1, NZINT
          do j=NVB+1, NRB
            do k=1, N2D
              ms=ms+1
              CENT2(ms)=cdexp(-CI*(ETINT(i)+EIG0(j)+EIG2D(k))*DT/2.d0)
            end do
          end do
        end do

        call BROTCZ3(AM(2),AM(3))

        end subroutine CPLXENG
c=========================================================================
        subroutine PLOT1D(NZ,NVB,NVB2,NA,W1,WK,OPTION)
        implicit real*8 (A-H,O-Z)
        complex*8 W1(1)
        character*2 OPTION
        real*8 WK(1)

        if (OPTION == "ZZ") then
          do j=1,NZ
            ms=(j-1)*NVB*NVB2*NA
            do i=1,NVB*NVB2*NA
              WK(j)=WK(j)+abs(W1(ms+i))**2.d0
            end do
          end do
        end if

        if (OPTION == "R1") then
          do k=1,NZ
            do j=1,NVB
              ms=(k-1)*NVB*NVB2*NA+(j-1)*NVB2*NA
              do i=1,NVB2*NA
                WK(j)=WK(j)+abs(W1(ms+i))**2.d0
              end do
            end do
          end do
        end if

        if (OPTION == "R2") then
          do k=1,NZ*NVB
            do j=1,NVB2
              ms=(k-1)*NVB2*NA+(j-1)*NA
              do i=1,NA
                WK(j)=WK(j)+abs(W1(ms+i))**2.d0
              end do
            end do
          end do
        end if

        end subroutine
c===================================================================
        subroutine PLOT2D(NZ,NVB,NVB2,NA,W1,WK,OPTION)
        implicit real*8 (A-H,O-Z)
        complex*8 W1(1)
        character*4 OPTION
        real*8 WK(1)

        if (OPTION == "ZZR1") then
          do k=1,NZ
            do j=1,NVB
              ns=(k-1)*NVB+j
              ms=(k-1)*NVB*NVB2*NA+(j-1)*NVB2*NA
              do i=1,NVB2*NA
                WK(ns)=WK(ns)+abs(W1(ms+i))**2.d0
              end do
            end do
          end do
        end if

        if (OPTION == "ZZR2") then
          do l=1,NZ
            do k=1,NVB
              do j=1,NVB2
                ns=(l-1)*NVB2+j
                ms=(l-1)*NVB*NVB2*NA+(k-1)*NVB2*NA+(j-1)*NA
                do i=1,NA
                  WK(ns)=WK(ns)+abs(W1(ms+i))**2.d0
                end do
              end do
            end do
          end do
        end if

        if (OPTION == "R1R2") then
          do k=1,NZ
            do j=1,NVB*NVB2
              ms=(k-1)*NVB*NVB2*NA+(j-1)*NA
              do i=1,NA
                WK(j)=WK(j)+abs(W1(ms+i))**2.d0
              end do
            end do
          end do
        end if

        end subroutine
c==================================================================
        function ITOC(I)
        character*2 ITOC
        character*1 CH

        LOW=MOD(I,10)
        IHIGH=MOD(I/10,10)

        ITOC(1:1)=CHAR(48+IHIGH)
        ITOC(2:2)=CHAR(48+LOW)

        RETURN
        END
c=======================================================================
        function ITOC3(I)
        character*3 ITOC3
        character*1 CH

        ITOC3(1:1)=CHAR( 48 + I/100 )
        ITOC3(2:2)=CHAR( 48 + mod(I,100)/10 )
        ITOC3(3:3)=CHAR( 48 + mod(I,10) )

        RETURN
        END
c=======================================================================
        function ITOC4(I)
        character*4 ITOC4
        character*1 CH

        ITOC4(1:1)=CHAR( 48 + I/1000 )
        ITOC4(2:2)=CHAR( 48 + mod(I,1000)/100 )
        ITOC4(3:3)=CHAR( 48 + mod(I,100)/10 )
        ITOC4(4:4)=CHAR( 48 + mod(I,10) )

        RETURN
        END
c=======================================================================
        subroutine BACKUP
        use WFARRAYS, only : WA, WI
        implicit none


        end subroutine BACKUP
c=======================================================================
        subroutine PLOTWAVE(TYPE)
        use WFARRAYS, only : WA, WI
        use QUADPOINTS, only : ZQ, RQ, RQI
        use COMPARAS, only : NZ, NZINT, NVB, NRB, N2D, NBASJJ
        implicit none
        integer ierr, ms, i, j, k, NS0
        character*3 TYPE
        character*10 tmpFilename
        real*8, allocatable :: WK1(:), WK2(:)
        real*8, allocatable :: WK12(:)


        NS0=NZINT*NVB*N2D*NBASJJ
        
        if (TYPE == "ASY") then
          allocate(WK1(NZ-NZINT),WK2(NVB),WK12((NZ-NZINT)*NVB),stat=ierr)
          WK1(:)=0.d0; WK2(:)=0.d0; WK12(:)=0.d0
          call PLOT1D(NZ-NZINT,NVB,N2D,NBASJJ,WA(NS0+1),WK1,'ZZ')
          call PLOT1D(NZ-NZINT,NVB,N2D,NBASJJ,WA(NS0+1),WK2,'R1')
          call PLOT2D(NZ-NZINT,NVB,N2D,NBASJJ,WA(NS0+1),WK12,'ZZR1')

          tmpFilename="ZZ_ASY"
          open(101, FILE = tmpFilename,STATUS="UNKNOWN",ACTION="WRITE", 
     .      FORM="FORMATTED",POSITION="APPEND",IOSTAT = ierr)
          do i=NZINT+1, NZ
            write(101,*) ZQ(i), WK1(i-NZINT)
          end do
          write(101,*)
          close(101)
          tmpFilename="R1_ASY"
          open(101,FILE = tmpFilename,STATUS="UNKNOWN",ACTION="WRITE", 
     .      FORM="FORMATTED",POSITION="APPEND",IOSTAT=ierr)
          do i=1, NVB
            write(101,*) RQ(i), WK2(i)
          end do
          write(101,*)
          close(101)
          tmpFilename="ZZR1_ASY"
          open(101,FILE=tmpFilename,STATUS="UNKNOWN",ACTION = "WRITE", 
     .      FORM="FORMATTED",POSITION = "APPEND",IOSTAT = ierr)
          ms=0
          do j=1, NZ-NZINT
          do i=1, NVB
            ms=ms+1
            write(101,*) ZQ(NZINT+j), RQ(i), WK12(ms)
          end do
            write(101,"()")
          end do
          write(101,*)
          close(101)

        else if (TYPE == "INT") then

          allocate(WK1(NZINT), WK2(NRB), WK12(NZINT*NRB), stat=ierr)
          WK1(:)=0.d0; WK2(:)=0.d0; WK12(:)=0.d0
          call PLOT1D(NZINT,NRB,N2D,NBASJJ,WI,WK1,'ZZ')
          call PLOT1D(NZINT,NRB,N2D,NBASJJ,WI,WK2,'R1')
          call PLOT2D(NZINT,NRB,N2D,NBASJJ,WI,WK12,'ZZR1')

          tmpFilename="ZZ_INT"
          open(101, FILE=tmpFilename,STATUS="UNKNOWN",ACTION="WRITE", 
     .      FORM="FORMATTED",POSITION="APPEND",IOSTAT=ierr)
          do i=1, NZINT
            write(101,*) ZQ(i), WK1(i)
          end do
          write(101,*)
          close(101)
          tmpFilename="R1_INT"
          open(101,FILE=tmpFilename,STATUS="UNKNOWN",ACTION="WRITE", 
     .      FORM="FORMATTED",POSITION="APPEND",IOSTAT = ierr)
          do i=1, NRB
            write(101,*) RQI(i), WK2(i)
          end do
          write(101,*)
          close(101)
          tmpFilename="ZZR1_INT"
          open(101,FILE=tmpFilename,STATUS="UNKNOWN",ACTION="WRITE", 
     .      FORM="FORMATTED",POSITION="APPEND",IOSTAT=ierr)
          ms=0
          do j=1, NZINT
          do i=1, NRB
            ms=ms+1
            write(101,*) ZQ(j), RQI(i), WK12(ms)
          end do
            write(101,"()")
          end do
          write(101,*)
          close(101)

        else
          stop 'wrong option in  plotwave'
        end if

        deallocate(WK1, WK2, WK12, stat=ierr)

        end subroutine PLOTWAVE
c=====================================================================
        subroutine CAL_PROB
        use WFARRAYS, only : PHIER
        use ENGARRAYS, only : ESK0, AE0, EIGEN_ENG14
        use COMPARAS, only : NZINT, N2D, NBASJJ, NE0, DT
        use REDUCED_MASS, only : RMS_YCZ3
        implicit none
        integer NE, i, NFLUXR
        real*8 EV, FXR, FACTID

        EV=27.2116d0
        FACTID=1.d0/RMS_YCZ3*DT*DT
        NFLUXR=NZINT*N2D*NBASJJ

        do NE=1, NE0
          FXR=0.d0
          do i=1, NFLUXR
            FXR=FXR+IMAG(CONJG(PHIER(i,NE,1))*PHIER(i,NE,2))
          end do
          write(500,"(f10.5,4x,f14.8)") ESK0(NE)*EV, FXR*FACTID/AE0(NE)
c          write(500,"(f10.5,4x,f14.8)") (ESK0(NE)+EIGEN_ENG14)*EV, FXR*FACTID/AE0(NE)
        end do
        write(500, "()")

        end subroutine CAL_PROB
c====================================================================
        Subroutine REACTION
        use COMPARAS, only : TTOTAL,DT,TPRINT,NZ,NZINT,N2D,NBASJJ,NE0,
     .    NA1node, NA2node, NB1node, NB2node,NRB,NVB
        use WFARRAYS, only : WA,WI,CWK0,PHIER
        use ENGARRAYS, only : VMASY_SAVE, VMINT_SAVE

        implicit none
        integer i, j, k, ms, NSTEP, NPRINT, PRT_FLAG, ierr, nodeax, nodebx
        real*8 TTOT, AN, AN1, AN2, ANABSZ
        real*8, allocatable :: WK(:)

        NSTEP=INT(TTOTAL/DT)
        NPRINT=INT(TPRINT/DT)
        nodeax=NA1node*NA2node
        nodebx=NB1node*NB2node

! CWK0 is a frequently used temporary array
        allocate(CWK0((NZ+NZINT)*N2D*NBASJJ), stat=ierr)

! PHIER is TD energy dependent wavefunction
        allocate(PHIER(NZINT*N2D*NBASJJ,NE0,2), stat=ierr)

        PHIER(:,:,:)=0.d0

        j=0
        TTOT=0.d0
        PRT_FLAG=1
        ANABSZ=0.d0

        open(100,file='TTIME',action='write',form='formatted',status='unknown')
        write(100,*) TTOTAL
        close(100)

        
        do i=1, NSTEP
          j=j+1
          TTOT=TTOT+DT
          call SCHPROG(TTOT,PRT_FLAG,nodeax,nodebx,VMASY_SAVE,VMINT_SAVE,ANABSZ)
          PRT_FLAG=0
          call REACTIVE_FLUX(TTOT)
          if (j.eq.NPRINT) then 
            PRT_FLAG=1
            j=0
            call CNORM(NZ*NVB*N2D*NBASJJ,WA,AN1)
            call CNORM(NZINT*(NRB-NVB)*N2D*NBASJJ,WI,AN2)
            call CAL_PROB
            write(*,9100) TTOT, AN1, AN2, AN1+AN2
            open(100,file='TTIME',action='read',form='formatted',status='old')
            read(100,*) TTOTAL
            close(100)
            NSTEP=INT(TTOTAL/DT)
          end if
        end do

 9100   format(f10.2,3f15.10)

        call CAL_PROB

        end subroutine REACTION
c==========================================================================
        subroutine BROTCZ3(AMC,AMZ)
        use QUADPOINTS, only :  AUQ, RUQ
        use COMPARAS, only : NVBRU, NVBAU, N2D
        use SMALL_ARRAYS, only : BIA_CZ3, BIC_CZ3

        implicit none
        
        integer ibas, jbas, ms, ierr
        real*8 AMC, AMZ, xcz3, scz3

        print*, AMC, AMZ
        allocate(BIA_CZ3(N2D),BIC_CZ3(N2D),stat=ierr)

        ms=0
        do 10 jbas=1,NVBRU
        do 10 ibas=1,NVBAU
          ms=ms+1
           
          xcz3=(RUQ(jbas)*dsin(AUQ(ibas)))**2
          scz3=(RUQ(jbas)*dcos(AUQ(ibas)))**2
          BIA_CZ3(ms)=1.5d0*AMZ*(xcz3+(2*AMC/(AMC+3*AMZ))*scz3)
          BIC_CZ3(ms)=3.0d0*AMZ*xcz3

          BIA_CZ3(ms)=0.5d0/BIA_CZ3(ms)
          BIC_CZ3(ms)=0.5d0/BIC_CZ3(ms)
 10     continue
        
        end subroutine BROTCZ3
