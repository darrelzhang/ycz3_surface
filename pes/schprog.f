        subroutine SCHPROG(TTOT,PRT_FLAG,nodeax,nodebx,VMA,VMI,ANABSZ)
        use WFARRAYS
        use ENGARRAYS, only : CENT1,CENT2,FZABS,FRABS,VCUT4D_ASY,VCUT4D_INT
        use TRANSMATS
        use ROTBASIS
        use CPBASIS
        use COMPARAS
        use SMALL_ARRAYS, only : BIA_CZ3, BIC_CZ3, BROTZ, BROTR,
     .    LINK1_ASY, LINK2_ASY, LINK1_INT, LINK2_INT

        implicit none
        integer PRT_FLAG, i, j, k, MS, NS, NAAU, INC, NS0, NAVB, NARB, 
     .    NARB2, NQA, NA, NSTART, NGRID_ASY, NGRID_INT, kbzVal, 
     .    idatas, icards, npoten, nodebx, nodeax, IP, NP, ierr

        real*4 VMA(nodebx,nodeax,NCARDS_ASY), VMI(nodebx,nodeax,NCARDS_INT)

        real*4 VCUT0
	real*8 AN, TTOT, ANABSZ
        complex*8 CTMP

        complex*8, allocatable :: WT1(:,:), WT2(:,:), WT3(:,:)

        integer omp_get_max_threads, omp_get_thread_num
        integer itest, jtest

        NP=omp_get_max_threads()
c	print*, 'NP=', NP

        NAVB=N2D*NBASJJ*NVB
        NARB=N2D*NBASJJ*NRB
        NARB2=N2D*NBASJJ*(NRB-NVB)
        NAAU=N2D*NBASJJ
        NS0=NZINT*NAVB
        NGRID_ASY=(NZ-NZINT)*NVB*N2D
        NGRID_INT=NZINT*NRB*N2D
        NQA=NA1node*NA2node*NB1node*NB2node

        allocate(WT1(NQA*2,NP),WT2(NQA*2,NP),WT3(NQA*2,NP),stat=ierr)

! exp(-i*K(R,r,au)*DT/2.d0) for WA and WI
        MS=0
        INC=1
C$OMP PARALLEL DEFAULT(SHARED)
C$OMP&  PRIVATE(i,CTMP)
C$OMP DO
        do i=1, NZ*NVB*N2D
          CTMP=CENT1(i)
          call CSCAL(NBASJJ, CTMP, WA((i-1)*NBASJJ+1), INC)
        end do
C$OMP END DO NOWAIT
C$OMP END PARALLEL

        call CTRANS(NZ,NZ,1,NAVB,WA,CWK0,TMATZ,'BG')

        INC=1
C$OMP PARALLEL DEFAULT(SHARED)
C$OMP&  PRIVATE(i,CTMP)
C$OMP DO
        do i=1, NZINT*(NRB-NVB)*N2D
          CTMP=CENT2(i)
          call CSCAL(NBASJJ, CTMP, WI((i-1)*NBASJJ+1), INC)
        end do
C$OMP END DO NOWAIT
C$OMP END PARALLEL

        call CTRANS(NZINT,NZINT,1,NARB2,WI,CWK0,TMATIZ,'BG')

! L shape ( WA -> WI)
        DO 10 I=NZINT,1,-1
          DO J=NRB,NVB+1,-1   ! shift WI up by NVB
            MS=(I-1)*NARB+(J-1)*NAAU
            NS=(I-1)*NARB2+(J-NVB-1)*NAAU
	    DO K=NAAU,1,-1
	      WI(MS+K)=WI(NS+K)
	    END DO
	  END DO

          DO J=NVB,1,-1
            MS=(I-1)*NARB+(J-1)*NAAU
            NS=(I-1)*NAVB+(J-1)*NAAU
  	    DO K=NAAU,1,-1
	      WI(MS+K)=WA(NS+K)
	    END DO
	  END DO
10	continue


! FBR->DVR (r, au) for WA and WI
        call CTRANS(NVB,NVB,(NZ-NZINT),NAAU,WA(NS0+1),CWK0,TMATR,'BG')
        call CTRANS(N2D,N2D,(NZ-NZINT)*NVB,NBASJJ,WA(NS0+1),CWK0,TMAT2D,'BG')
        call CTRANS(NRB,NRB,NZINT,NAAU,WI,CWK0,TMATIR,'BG')
        call CTRANS(N2D,N2D,NZINT*NRB,NBASJJ,WI,CWK0,TMAT2D,'BG')

! rotational kinetic energy
        call KINEROT(NZ-NZINT,NVB,WA(NS0+1),BROTR(NRB+1))
        call KINEROT(NZINT,NRB,WI, BROTR)

! centrifugal potential
        call CENTRIPOT(NZ-NZINT,NVB,WA(NS0+1), BROTZ(NZINT+1))
        call CENTRIPOT(NZINT,NRB,WI,BROTZ)


! 7D potential energy surface for WA and WI

C$OMP PARALLEL DEFAULT(SHARED)
C$OMP&  PRIVATE(i,icards,idatas,npoten,IP,VCUT0,CTMP)
C$OMP DO
      do i=1, NGRID_ASY
          IP=omp_get_thread_num()+1
          npoten=LINK1_ASY(3,i)
          VCUT0=VCUT4D_ASY(i)
	  icards=LINK1_ASY(1,i)
	  idatas=LINK1_ASY(2,i)

	  if (npoten.eq.0) then

	    CTMP=CDEXP((0.0,-1.0)*DT*VCUT0)
	    call CSCAL(NBASJJ,CTMP,WA(NS0+(i-1)*NBASJJ+1),INC)

	  else

	    call VMATRIX(npoten,LINK2_ASY(1,i),VCUT0,VMA(1,idatas,icards),
     .        WA(NS0+(i-1)*NBASJJ+1),WT1(1,IP),WT2(1,IP),WT3(1,IP))

	  end if
        end do
C$OMP END DO NOWAIT
C$OMP END PARALLEL

	if (PRT_FLAG.ne.0) call PLOTWAVE("ASY")

 5000   continue

C$OMP PARALLEL DEFAULT(SHARED)
C$OMP&  PRIVATE(i,icards,idatas,npoten,IP,VCUT0,CTMP)
C$OMP DO
        do i=1, NGRID_INT
          IP=omp_get_thread_num()+1
	  icards=LINK1_INT(1,i)
	  idatas=LINK1_INT(2,i)
          npoten=LINK1_INT(3,i)
          VCUT0=VCUT4D_INT(i)

	  if (npoten.eq.0) then

	    CTMP=CDEXP((0.0,-1.0)*DT*VCUT0)
	    call CSCAL(NBASJJ,CTMP,WI((i-1)*NBASJJ+1),INC)

	  else


	    call VMATRIX(npoten,LINK2_INT(1,i),VCUT0,VMI(1,idatas,icards),
     .        WI((i-1)*NBASJJ+1),WT1(1,IP),WT2(1,IP),WT3(1,IP))

	  end if
        end do
C$OMP END DO NOWAIT
C$OMP END PARALLEL


	if (PRT_FLAG.ne.0) call PLOTWAVE("INT")

! centrifugal potential
        call CENTRIPOT(NZ-NZINT,NVB,WA(NS0+1),BROTZ(NZINT+1))
        call CENTRIPOT(NZINT,NRB,WI,BROTZ)

! rotational kinetic energy
        call KINEROT(NZ-NZINT,NVB,WA(NS0+1),BROTR(NRB+1))
        call KINEROT(NZINT,NRB,WI,BROTR)

! absorption potential
        if (TTOT.gt.T_ABSZ.and.ANABSZ.gt.2.0d0) then
          ms=0
C$OMP PARALLEL DEFAULT(SHARED)
C$OMP&  PRIVATE(i,j,ms)
C$OMP DO

          do i=NZABS, NZ
            do j=1, NAVB
              ms=(i-1)*NAVB+j
              WA(ms)=WA(ms)*FZABS(i)
            end do
          end do
C$OMP END DO NOWAIT
C$OMP END PARALLEL

	else

	  do i=NZABS,NZ
	    do j=1, NAVB
	      ms=(i-1)*NAVB+j
	      ANABSZ=ANABSZ+abs(WA(ms))**2.d0
	    end do
	  end do

        end if

        if (TTOT.gt.T_ABSR) then
          ms=0
C$OMP PARALLEL DEFAULT(SHARED)
C$OMP&  PRIVATE(i,j,k,ms)
C$OMP DO
         do i=1, NZINT
            do j=NRABS, NRB
              do k=1, N2D*NBASJJ
                ms=(i-1)*NARB+(j-1)*N2D*NBASJJ+k
                WI(ms)=WI(ms)*FRABS(j)
              end do
            end do
          end do
C$OMP END DO NOWAIT
C$OMP END PARALLEL
        end if

! DVR->FBR(r,au) for WA and WI
        call CTRANS(NVB,NVB,(NZ-NZINT),NAAU,WA(NS0+1),CWK0,TMATR,'GB')
        call CTRANS(N2D,N2D,(NZ-NZINT)*NVB,NBASJJ,WA(NS0+1),CWK0,TMAT2D,'GB')

        call CTRANS(NRB,NRB,NZINT,NAAU,WI,CWK0,TMATIR,'GB')
        call CTRANS(N2D,N2D,NZINT*NRB,NBASJJ,WI,CWK0,TMAT2D,'GB')

! L shape (WI -> WA)
        DO 100 I=1,NZINT
          DO J=1,NVB
            MS=(I-1)*NARB+(J-1)*NAAU
            NS=(I-1)*NAVB+(J-1)*NAAU
	    DO K=1,NAAU
	      WA(NS+K)=WI(MS+K)
	    END DO
          END DO

          DO J=NVB+1,NRB   
            MS=(I-1)*NARB+(J-1)*NAAU
            NS=(I-1)*NARB2+(J-NVB-1)*NAAU
	    DO K=1,NAAU
	      WI(NS+K)=WI(MS+K)
	    END DO
          END DO
100     continue

! exp(-i*K(R,r,au)*DT/2.d0) for WA and WI
        call CTRANS(NZ,NZ,1,NAVB,WA,CWK0,TMATZ,'GB')
        INC=1
C$OMP PARALLEL DEFAULT(SHARED)
C$OMP&  PRIVATE(i,CTMP)
C$OMP DO
        do i=1, NZ*NVB*N2D
          CTMP=CENT1(i)
          call CSCAL(NBASJJ, CTMP, WA((i-1)*NBASJJ+1), INC)
        end do
C$OMP END DO NOWAIT
C$OMP END PARALLEL


        call CTRANS(NZINT,NZINT,1,NARB2,WI,CWK0,TMATIZ,'GB')
        INC=1
C$OMP PARALLEL DEFAULT(SHARED)
C$OMP&  PRIVATE(i,CTMP)
C$OMP DO
        do i=1, NZINT*(NRB-NVB)*N2D
          CTMP=CENT2(i)
          call CSCAL(NBASJJ, CTMP, WI((i-1)*NBASJJ+1), INC)
        end do
C$OMP END DO NOWAIT
C$OMP END PARALLEL


        deallocate(WT1,WT2,WT3,stat=ierr)

        end subroutine SCHPROG
c=======================================================================
        subroutine REACTIVE_FLUX(TTOT)
        use TRANSMATS, only : TMATZ, TMATIZ
        use WFARRAYS, only : WA, WI, CWK0, PHIER
        use COMPARAS, only : NZ, NZINT, NVB, NRB, N2D, NBASJJ, NE0
        use ENGARRAYS, only : ESK0, EIGEN_ENG14
        use SMALL_ARRAYS, only : WVIBR
        implicit none
        integer NE, I, J, K, MS, NS, NSC, NFLUXR
        real*8 TTOT, ET
        complex*16 CI, CTMP

        CI=(0.d0, 1.d0)
        
        NFLUXR=NZINT*N2D*NBASJJ

        do NSC=1,2

          call CZERO1D(NFLUXR, CWK0(NFLUXR+1))

C$OMP PARALLEL DEFAULT(SHARED)
C$OMP&  PRIVATE(K,J,I,NS,MS)
C$OMP DO          
          do 20 K=1, NZINT
          do 20 J=NVB+1, NRB
            NS=(K-1)*N2D*NBASJJ
            MS=(K-1)*(NRB-NVB)*N2D*NBASJJ+(J-NVB-1)*N2D*NBASJJ
            do 30 I=1, N2D*NBASJJ
              CWK0(NFLUXR+NS+I)=CWK0(NFLUXR+NS+I)+WI(MS+I)*WVIBR(J,NSC)
 30         continue
 20       continue
C$OMP END DO NOWAIT
C$OMP END PARALLEL

          call CZERO1D(NFLUXR, CWK0)

C$OMP PARALLEL DEFAULT(SHARED)
C$OMP&  PRIVATE(K,J,I,NS,MS)
C$OMP DO          
          
          do 40 K=1, NZINT
            MS=(K-1)*N2D*NBASJJ
          do 40 J=1, NZINT
            NS=(J-1)*N2D*NBASJJ
          do 40 I=1, N2D*NBASJJ
            CWK0(MS+I)=CWK0(MS+I)+CWK0(NFLUXR+NS+I)*TMATIZ(J,K)
 40       continue
C$OMP END DO NOWAIT
C$OMP END PARALLEL

          call CZERO1D(NZ*N2D*NBASJJ, CWK0(NFLUXR+1))

C$OMP PARALLEL DEFAULT(SHARED)
C$OMP&  PRIVATE(K,J,I,NS,MS)
C$OMP DO          

          do 120 K=1, NZ
          do 120 J=1, NVB
            NS=(K-1)*N2D*NBASJJ
            MS=(K-1)*NVB*N2D*NBASJJ+(J-1)*N2D*NBASJJ
            do 130 I=1, N2D*NBASJJ
              CWK0(NFLUXR+NS+I)=CWK0(NFLUXR+NS+I)+WA(MS+I)*WVIBR(J,NSC)
 130         continue
 120       continue
C$OMP END DO NOWAIT
C$OMP END PARALLEL

C$OMP PARALLEL DEFAULT(SHARED)
C$OMP&  PRIVATE(K,J,I,NS,MS)
C$OMP DO          
         do 140 K=1, NZINT
            MS=(K-1)*N2D*NBASJJ
          do 140 J=1, NZ
            NS=(J-1)*N2D*NBASJJ
          do 140 I=1, N2D*NBASJJ
            CWK0(MS+I)=CWK0(MS+I)+CWK0(NFLUXR+NS+I)*TMATZ(J,K)
 140      continue
C$OMP END DO NOWAIT
C$OMP END PARALLEL

C$OMP PARALLEL DEFAULT(SHARED)
C$OMP&  PRIVATE(NE,I,ET,CTMP)
C$OMP DO                    
          do NE=1,NE0
            ET=(ESK0(NE)+EIGEN_ENG14)*TTOT
            CTMP=CDEXP(CI*ET)
            do I=1,NFLUXR
              PHIER(I,NE,NSC)=PHIER(I,NE,NSC)+CWK0(I)*CTMP
            end do
          end do
C$OMP END DO NOWAIT
C$OMP END PARALLEL

        end do

        end subroutine REACTIVE_FLUX
c=======================================================================
        subroutine KINEROT(N1, N2, W1, BROTR)
        use COMPARAS, only : NVB, NRB, N2D, NBASJJ,DT
        use SMALL_ARRAYS, only : BIA_CZ3, BIC_CZ3
        use ENGARRAYS, only :KROTJJ

        integer N1, N2, i, j, k, ms, j12Val, j2Val, kb2Val, ik2

        real*4 krot, BROTR(1)
        complex*8 W1(1), CI

        CI=(0.0, 1.0)

        ms=0
C$OMP PARALLEL DEFAULT(SHARED)
C$OMP&  PRIVATE(i,j,k,l,ms,krot)
C$OMP DO
        do i=1, N1
          do j=1, N2
            do k=1, N2D
              do l=1, NBASJJ
                ms=(i-1)*N2*N2D*NBASJJ+(j-1)*N2D*NBASJJ+(k-1)*NBASJJ+l
                krot=BROTR(j)*KROTJJ(1,l)+BIA_CZ3(k)*KROTJJ(2,l)+(BIC_CZ3(k)-BIA_CZ3(k))*KROTJJ(3,l)
                W1(ms)=W1(ms)*cexp(CI*krot)
              end do
            end do
          end do
        end do
C$OMP END DO NOWAIT
C$OMP END PARALLEL


        end subroutine KINEROT
c=====================================================================
        subroutine CENTRIPOT(N1, N2, W1, BROTZ)
        use COMPARAS, only : NVB, NRB, N2D, NBASJJ, NDIMCP, NDIMCP2,
     .    NCP_COUPLED, DT
        use CPBASIS 
        use ROTBASIS, only : IDXJJ, KB2

        implicit none

        integer N1, N2, i, j, k, l, ms
        complex*8 CI, W1(1)
	real*4 BROTZ(1), kcp

        CI=(0.0, -1.0)*DT/2.0
        
        if (NCP_COUPLED.eq.1) then ! no coupling

C$OMP PARALLEL DEFAULT(SHARED)
C$OMP&  PRIVATE(i,j,k,l,ms,kcp)
C$OMP DO
          do i=1, N1
            do j=1, N2
              do k=1, N2D
                do l=1, NBASJJ
                  ms=(i-1)*N2*N2D*NBASJJ+(j-1)*N2D*NBASJJ+(k-1)*NBASJJ+l
                  kcp=BROTZ(i)*CPVAL(l)
                  W1(ms)=W1(ms)*cexp(CI*kcp)
                end do
              end do
            end do
          end do
C$OMP END DO NOWAIT
C$OMP END PARALLEL

        else
          ms=0
          do i=1, N1
            do j=1, N2
              do k=1, N2D
                ms=ms+1
                call CPTRANS(W1((ms-1)*NBASJJ+1),BROTZ(i))
              end do
            end do
          end do

        end if

        end subroutine CENTRIPOT
c===================================================================
        subroutine CPTRANS(W1,BROTZ)
        use COMPARAS, only : NBASJJ, NDIMCP, NDIMCP2, DT, NCP_COUPLED
        use CPBASIS
        implicit none
        integer*4 ierr, i, j, ibasjj, ms
        real*4 BROTZ
        complex*8 W1(1), ctmp, CI
        complex*8, allocatable :: cworks(:)

        allocate(cworks(NBASJJ), stat=ierr)

        cworks=0.0
        do i=1, NCP_COUPLED
          ctmp=0.0
          do j=1, NCP_COUPLED
            ibasjj=idxCP2(1,j)
            ms=(i-1)*NCP_COUPLED+j
            ctmp=ctmp+W1(ibasjj)*CPMAT(ms)
          end do
          cworks(i)=ctmp*cdexp(-(0.0,1.0)*DT*BROTZ*CPVAL(i)/2.0)
        end do

        do i=1, NCP_COUPLED
          ctmp=0.0
          do j=1, NCP_COUPLED
            ms=(j-1)*NCP_COUPLED+i
            ctmp=ctmp+cworks(j)*CPMAT(ms)
          end do
          ibasjj=idxCP2(1,i)
          W1(i)=ctmp
        end do

        deallocate(cworks, stat=ierr)

        end subroutine CPTRANS

