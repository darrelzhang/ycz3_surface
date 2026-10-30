c=========================================================================
        subroutine READINP
        use COMPARAS
        implicit none
	include 'mpif.h'
        integer I, ierr

        ! read input

	open(100,file='input',status='old')
        read(100,*) IDPOT, PES_TYPE, NPOTS, IPOT
        read(100,*) (ATOM(I), I=1, 3) ! Y:ATOM(2), C:ATOM(3), Z:ATOM(4) 
        read(100,*) NZ, NZINT, Zmin, Zmax   
        read(100,*) NR, NVB, NRB, Rmin, Rmax
        read(100,*) NRU, NVBRU, RUmin, RUmax
        read(100,*) NAU, NVBAU, AUmin, AUmax
        read(100,*) NA1node, NA2node, NB1node, NB2node
        read(100,*) KBZmax, NB3node, Ca0 
        read(100,*) J12max, JS1, JD1, J1max, JS2, JD2, J2max, KB2max 
        read(100,*) Zinf, R0, RU0, AU0, A10, A20, B10, B20  ! RU0 is the equlibrium bond distance of C-Z
        read(100,*) Zinit, ESKinit, DELTA, KBZinit, J1init, istat1
        read(100,*) NE0, E1, E2
        read(100,*) TTOTAL, DT, TPRINT, T_ABSZ, T_ABSR, NSTEP0, RFLUX, VCUT
        read(100,*) ZABS, CZABS, PZABS, RABS, CRABS, PRABS
	read(100,*) nfoldb2
	read(100,*) X00, Y00
	close(100)

        call MPI_COMM_RANK(MPI_COMM_WORLD, IPOT, ierr)
	IPOT=IPOT+1

        call MPI_COMM_SIZE(MPI_COMM_WORLD, NPOTS, ierr)
!	print*,'NPOT= ',NPOTS, 'IPOT= ',IPOT

        end subroutine READINP
c==========================================================================
        subroutine CHKINP
        use COMPARAS
        use REDUCED_MASS
        implicit real*8(a-h,o-z)
        integer :: i, NANODES

        EV=27.2116d0
        CM=219474.6d0

        ESKinit=ESKinit/EV
        E1=E1/EV
        E2=E2/EV
        VCUT=VCUT/CM

        call MASS(1,ATOM(1),AM(1))
        call MASS(1,ATOM(2),AM(2))
        call MASS(1,ATOM(3),AM(3))

        TMS_YCZ3=AM(1)+AM(2)+3*AM(3)
        RMS_YCZ3=AM(1)*(AM(2)+3*AM(3))/(AM(1)+AM(2)+3*AM(3))
        RMS_CZ3=3*AM(2)*AM(3)/(AM(2)+3*AM(3))
        TMS_3Z=3*AM(3)
        print*, 'TMS_YCZ3=', TMS_YCZ3
        print*, 'RMS_YCZ3=', RMS_YCZ3
        print*, 'RMS_CZ3=', RMS_CZ3
        print*, 'TMS_3Z=', TMS_3Z

! check parameters
        if (Zmin.lt.0) stop 'Zmin < 0'
        if (Zmin.gt.Zmax) stop 'Zmin > Zmax'
        if (NZINT.gt.NZ) stop 'NZINT > NZ'
        if (Rmin.lt.0) stop 'Rmin < 0'
        if (Rmin.gt.Rmax) stop 'Rmin > Rmax'
        if (NR < NVB .or. NR < NRB) stop 'NR < NVB or NR < NRB'
        if (NVB.gt.NRB) stop 'NVB > NRB'
        if (NVBRU.gt.NVBRU) stop 'NRU > NVBRU'
        if (RUmin.gt.RUmax) stop 'RUmin > RUmax'
        if (NVBAU.gt.NVBAU) stop 'NAU > NVBAU'
        if (AUmin.gt.AUmax) stop 'AUmin > AUmax'
        if (KBZmin.lt.0) stop 'KBZmin < 0'
        if (Jtot.lt.0) stop 'Jtot < 0'
        if (KBZmin.gt.KBZmax) stop 'KBZmin > KBZmax'
        if (KBZmax.gt.J1max) stop 'KBZmax > J1max'
        if (KB2max.gt.J2max) stop 'KB2max > J2max'
        if (J1max.lt.0) stop 'J1max < 0'
        if (J2max.lt.0) stop 'J2max < 0'
        if (NA1node.lt.J1max+1) stop 'NA1node < J1max+1'
        if (NA2node.lt.J2max+1) stop 'NA2node < J2max+1'
        if (NB1node.lt.2*min(J1max,J2max)+1) 
     .    stop 'NB1node < 2*min(J1max,J2max)+1'
        if (NB2node.lt.2*KB2max+1) stop 'NB2node < 2*KB2max+1'
        if (KBZinit.gt.min(Jtot,J1init)) 
     .    stop'KBZinit< min(Jtotal,J1init)'
        if (KBZinit.lt.KBZmin.or.KBZinit.gt.KBZmax) 
     .    stop 'KBZinit not in (KBZmin,KBZmax)'
        if (J1init.gt.J1max) stop'J1init.gt.J1max'

	ca0=ca0*DACOS(-1.d0)/180.d0

        M1max = min(J1max, J2max)
        istat1=abs(istat1)
        AU0=AU0*DACOS(-1.d0)/180.d0
	       AUmin=AUmin*DACOS(-1.d0)/180.d0
	       AUmax=AUmax*DACOS(-1.d0)/180.d0

        NKBZinit=min(J1init, KBZmax)-KBZmin+1

        if(mod(NB2node,nfoldb2).ne.0)NB2node=nfoldb2*(NB2node/nfoldb2+1)
        NB2node=NB2node/nfoldb2

        NANODES=NA1node*NA2node*NB1node*NB2node
        N2D=NVBRU*NVBAU
!        open(41,file='VMATRIX.ASY',status='unknown',
!     &            form='unformatted',access='direct',recl=NANODES)
!        open(42,file='VMATRIX.INT',status='unknown',
!     &            form='unformatted',access='direct',recl=NANODES)
!        open(43,file='VMATASY.IDX',status='unknown',
!     &            form='formatted')
!        open(44,file='VMATINT.IDX',status='unknown',
!     &            form='formatted')


        end subroutine CHKINP
c========================================================================
        SUBROUTINE MASS(N,ATOM,AMASS)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION AMASS(N)
        CHARACTER*5 ATOM(N)

        DO I=1,N

        IF(ATOM(I).EQ.'H' .or. ATOM(I).EQ.'h') THEN
          AMASS(I)=1837.15D0
        ELSE IF(ATOM(I).EQ.'D' .or. ATOM(I).EQ.'d') THEN
          AMASS(I)=3674.3D0
        ELSE IF(ATOM(I).EQ.'C' .or. ATOM(I).EQ.'c') THEN
          AMASS(I)=21870.D0
        ELSE IF(ATOM(I).EQ.'N' .or. ATOM(I).EQ.'n') THEN
          AMASS(I)=25515.D0
        ELSE IF(ATOM(I).EQ.'O' .or. ATOM(I).EQ.'o') THEN
          AMASS(I)=29160.D0
        ELSE IF(ATOM(I).EQ.'F' .or. ATOM(I).EQ.'f') THEN
          AMASS(I)=34627.87D0
        ELSE IF(ATOM(I).EQ.'R' .or. ATOM(I).EQ.'r') THEN
          AMASS(I)=64701.21D0
        ELSE IF(ATOM(I).EQ.'K' .or. ATOM(I).EQ.'k') THEN
          AMASS(I)=64701.21D0
	ELSE IF(ATOM(I).EQ.'L' .or. ATOM(I).EQ.'l')THEN
          AMASS(I)=34.968852D0*1822.7d0
        ELSE
          PRINT*,'MASS OF', ATOM(I), 'IS NOT DEFINED'
        END IF

        END DO

        RETURN
        END
c====================================================================
        
        
