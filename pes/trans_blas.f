c=======================================================================
	SUBROUTINE CTRANS(NB,NR,N1,N2,W1,W2,TMAT,OPTION)
c ========================================
c W1 --> W1  with W2 a small working space
c ========================================
	IMPLICIT REAL*8 (A-H,O-Z)
	CHARACTER*2 OPTION
	REAL*4 TMAT(NB,1),alpha,beta,W1(1),W2(2000*(NR+NB),1)
	integer omp_get_thread_num

	alpha=1.D0
	beta=0.D0

	N2D=2*N2
	MM=2000 

	IF(OPTION .EQ. 'BG') THEN   ! GRID ---> BASIS
                                    ! sum over the first index.

	  IF(N2.EQ.1) THEN
	     PRINT*,'YOU SHOULD USE CTRANS1'
	     STOP
	  ELSE

c rearrange W1
	    IF(NR.GT.NB) THEN
	      INC=-1
	      DO I=N1,1,-1
	      DO J=NB,1,-1
	        MS=(I-1)*NB*N2D+(J-1)*N2D
	        NS=(I-1)*NR*N2D+(J-1)*N2D
                CALL SCOPY(N2D,W1(MS+1),INC,W1(NS+1),INC)
c	      DO K=N2D,1,-1
c	        W1(NS+K)=W1(MS+K)
c	      END DO
	      END DO
	      END DO
	    END IF
	
	    NP=INT(N2D/MM)
	    IF(MOD(N2D,MM).NE.0) NP=NP+1

	    NBR=MAX(NB,NR)

C$OMP PARALLEL SHARED (NBR,NR,NB,N1,N2D,TMAT,W1,W2,NP,alpha,beta,MM)
C$OMP&  PRIVATE (I,MS,NS,NSTART,NEND,N22,N,IP,NF,J,K)
C$OMP DO
            DO N=1,NP
	      IP=omp_get_thread_num()+1
              NSTART=(N-1)*MM
              NEND=MIN(N*MM,N2D)
              N22=NEND-NSTART

	    DO I=1,N1
	     MS=(I-1)*NBR*N2D+1
	     CALL SGEMM('N','N',N22,NR,NB,alpha,W1(MS+NSTART),N2D,TMAT,
     $       NB,beta,W2(1,IP),MM)
	
	     DO J=1,NR
	       MS=(I-1)*NBR*N2D+1
	       NS=(J-1)*N2D
	       NF=(J-1)*MM
	     DO K=1,N22
	       W1(MS+NSTART+NS+K-1)=W2(NF+K,IP)
	     END DO
	     END DO
	     
	    END DO
	    END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL

            IF(NB.GT.NR) THEN
	      INC=1
              DO I=1,N1
              DO J=1,NR
                MS=(I-1)*NB*N2D+(J-1)*N2D
                NS=(I-1)*NR*N2D+(J-1)*N2D
                CALL SCOPY(N2D,W1(MS+1),INC,W1(NS+1),INC)
c              DO K=1,N2D
c                W1(NS+K)=W1(MS+K)
c              END DO
              END DO
              END DO
            END IF

	  END IF
	ELSE IF (OPTION .EQ. 'GB') THEN    ! BASIS  --->  GRID
          IF(N2.EQ.1) THEN
	    PRINT*,'YOU SHOULD USE CTRANS1'
	    STOP
          ELSE

            NP=INT(N2D/MM)
            IF(MOD(N2D,MM).NE.0) NP=NP+1

	    IF(NB.GT.NR) THEN
              DO I=N1,1,-1
              DO J=NR,1,-1
                MS=(I-1)*NR*N2D+(J-1)*N2D
                NS=(I-1)*NB*N2D+(J-1)*N2D
              DO K=N2D,1,-1
                W1(NS+K)=W1(MS+K)
              END DO
              END DO
              END DO
            END IF

	    NBR=MAX(NB,NR)

C$OMP PARALLEL SHARED (NBR,NR,NB,N1,N2D,TMAT,W1,W2,NP,alpha,beta,MM)
C$OMP&  PRIVATE (I,MS,NS,NSTART,NEND,N22,N,IP,NF,J,K)
C$OMP DO
            DO N=1,NP
	      IP=omp_get_thread_num()+1
              NSTART=(N-1)*MM
              NEND=MIN(N*MM,N2D)
              N22=NEND-NSTART

            DO I=1,N1
             MS=(I-1)*NBR*N2D+1
             CALL SGEMM('N','T',N22,NB,NR,alpha,W1(MS+NSTART),N2D,TMAT,
     $       NB,beta,W2(1,IP),MM)

	     DO J=1,NB
               MS=(I-1)*NBR*N2D+1
               NS=(J-1)*N2D
               NF=(J-1)*MM
             DO K=1,N22
               W1(MS+NSTART+NS+K-1)=W2(NF+K,IP)
             END DO
             END DO

            END DO
	    END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL

            IF(NB.LT.NR) THEN
              DO I=1,N1
              DO J=1,NB
                MS=(I-1)*NR*N2D+(J-1)*N2D
                NS=(I-1)*NB*N2D+(J-1)*N2D
              DO K=1,N2D
                W1(NS+K)=W1(MS+K)
              END DO
              END DO
              END DO
            END IF

          END IF
	ELSE
          STOP 'OPTION IS WRONG IN STRANS'
	END IF

	RETURN
	END 
c=======================================================================
	SUBROUTINE CTRANS0(NR,NB,N1,N2,W1,W2,TMAT,OPTION)
c ==============
c W1 --> W2
c ==============
	IMPLICIT REAL*8 (A-H,O-Z)
	CHARACTER*2 OPTION
	REAL*4 TMAT(NR,1),alpha,beta,W1(1),W2(1)

	alpha=1.D0
	beta=0.D0

	N2D=2*N2
	MM=2000 

	IF(OPTION .EQ. 'BG') THEN   ! GRID ---> BASIS
	  IF(N2.EQ.1) THEN
             CALL CTRANS1(NR,NB,N1,W1,W2,TMAT,OPTION)
	  ELSE
	
	    NP=INT(N2D/MM)
	    IF(MOD(N2D,MM).NE.0) NP=NP+1

C$OMP PARALLEL SHARED (NR,NB,N1,N2D,TMAT,W1,W2,NP,alpha,beta,MM)
C$OMP&  PRIVATE (I,MS,NS,NSTART,NEND,N22,N)
C$OMP DO
            DO N=1,NP
              NSTART=(N-1)*MM
              NEND=MIN(N*MM,N2D)
              N22=NEND-NSTART

	    DO I=1,N1
	     MS=(I-1)*NB*N2D+1
	     NS=(I-1)*NR*N2D+1
	     CALL SGEMM('N','N',N22,NB,NR,alpha,W1(MS+NSTART),N2D,TMAT,
     $       NR,beta,W2(NS+NSTART),N2D)
	     
	    END DO
	    END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL

	  END IF
	ELSE IF (OPTION .EQ. 'GB') THEN    ! BASIS  --->  GRID
          IF(N2.EQ.1) THEN
            CALL CTRANS1(NR,NB,N1,W1,W2,TMAT,OPTION)
          ELSE

            NP=INT(N2D/MM)
            IF(MOD(N2D,MM).NE.0) NP=NP+1

C$OMP PARALLEL SHARED (NR,NB,N1,N2D,TMAT,W1,W2,NP,alpha,beta,MM)
C$OMP&  PRIVATE (I,MS,NS,NSTART,NEND,N22,N)
C$OMP DO
            DO N=1,NP
              NSTART=(N-1)*MM
              NEND=MIN(N*MM,N2D)
              N22=NEND-NSTART

            DO I=1,N1
             MS=(I-1)*NB*N2D+1
             NS=(I-1)*NR*N2D+1
             CALL SGEMM('N','T',N22,NR,NB,alpha,W1(NS+NSTART),N2D,TMAT,
     $       NR,beta,W2(MS+NSTART),N2D)
            END DO
	    END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL

          END IF
	ELSE
          STOP 'OPTION IS WRONG IN STRANS'
	END IF

	RETURN
	END 
c==========================================================================
        SUBROUTINE CTRANS1(NQ2,NB2,NQ1,W1,W2,TMAT,OPTION)
c this routine transfer wavefucntion between basis and grid for
c a wavefunction with external channel NQ1
c Basis --> Grid  ID=1               Grid  --> Basis ID=-1
        IMPLICIT REAL*8 (A-H,O-Z)
        REAL*4 TMAT(NQ2,1)
	CHARACTER*2 OPTION
	COMPLEX*8 W1(1),W2(1)

	IF (OPTION .EQ. 'BG') THEN  ! BASIS  --->  GRID
          DO I=1,NQ1*NQ2
            W2(I)=(0.D0,0.D0)
          END DO

C$OMP PARALLEL IF (NQ1.GT.10) SHARED (NQ1,NQ2,NB2,W1,W2,TMAT)
C$OMP&  PRIVATE (I,MS,NS,J,K)
C$OMP DO

          DO 100 I=1,NQ1
            MS=(I-1)*NQ2
            NS=(I-1)*NB2
          DO 100 J=1,NQ2
          DO 100 K=1,NB2
100       W2(MS+J)=W2(MS+J)+TMAT(K,J)*W1(NS+K)
C$OMP END DO NOWAIT
C$OMP END PARALLEL

	ELSE IF(OPTION .EQ. 'GB') THEN   ! GRID ---> BASIS

	  DO I=1,NQ1*NB2
            W2(I)=(0.D0,0.D0)
          END DO

C$OMP PARALLEL IF (NQ1.GT.10) SHARED (NQ1,NQ2,NB2,W1,W2,TMAT)
C$OMP&  PRIVATE (I,MS,NS,J,K)
C$OMP DO
          DO 200 I=1,NQ1
            MS=(I-1)*NQ2
            NS=(I-1)*NB2
          DO 200 J=1,NQ2
          DO 200 K=1,NB2
200       W2(NS+K)=W2(NS+K)+TMAT(K,J)*W1(MS+J)
C$OMP END DO NOWAIT
C$OMP END PARALLEL

	ELSE
	  STOP 'CTRANS1'
	END IF

        RETURN
        END
c==========================================================================
	SUBROUTINE STRANS(NR,NB,N1,N2,W1,W2,TMAT,OPTION)
	IMPLICIT REAL*8 (A-H,O-Z)
	CHARACTER*2 OPTION
	REAL*4 TMAT(NR,1),W1(1),W2(1),alpha,beta

	alpha=1.D0
	beta=0.D0
	IF(OPTION .EQ. 'BG') THEN   ! GRID ---> BASIS
	  IF(N2.EQ.1) THEN
	    CALL SGEMM('T','N',NB,N1,NR,alpha,TMAT,NR,W1,NR,beta,W2,NB)
	    RETURN
	  ELSE
	    DO I=1,N1
	     MS=(I-1)*NR*N2+1
	     NS=(I-1)*NB*N2+1
	     CALL SGEMM('N','N',N2,NB,NR,alpha,W1(MS),N2,TMAT,NR,beta,
     $       W2(NS),N2)
	    END DO
	  END IF
	ELSE IF (OPTION .EQ. 'GB') THEN    ! BASIS  --->  GRID
          IF(N2.EQ.1) THEN
            CALL SGEMM('N','N',NR,N1,NB,alpha,TMAT,NR,W1,NB,beta,W2,NR)
            RETURN
          ELSE
            DO I=1,N1
             MS=(I-1)*NR*N2+1
             NS=(I-1)*NB*N2+1
             CALL SGEMM('N','T',N2,NR,NB,alpha,W1(NS),N2,TMAT,NR,beta,
     $       W2(MS),N2)
            END DO
          END IF
	ELSE
          STOP 'OPTION IS WRONG IN STRANS'
	END IF

	RETURN
	END 
c==========================================================================
	SUBROUTINE DTRANS(NR,NB,N1,N2,W1,W2,TMAT,OPTION)
	IMPLICIT REAL*8 (A-H,O-Z)
	CHARACTER*2 OPTION
	REAL*8 TMAT(NR,1),W1(1),W2(1),alpha,beta

	alpha=1.D0
	beta=0.D0
	IF(OPTION .EQ. 'BG') THEN   ! GRID ---> BASIS
	  IF(N2.EQ.1) THEN
	    CALL DGEMM('T','N',NB,N1,NR,alpha,TMAT,NR,W1,NR,beta,W2,NB)
	    RETURN
	  ELSE
	    DO I=1,N1
	     MS=(I-1)*NR*N2+1
	     NS=(I-1)*NB*N2+1
	     CALL DGEMM('N','N',N2,NB,NR,alpha,W1(MS),N2,TMAT,NR,beta,
     $       W2(NS),N2)
	    END DO
	  END IF
	ELSE IF (OPTION .EQ. 'GB') THEN ! BASIS  --->  GRID
          IF(N2.EQ.1) THEN
            CALL DGEMM('N','N',NR,N1,NB,alpha,TMAT,NR,W1,NB,beta,W2,NR)
            RETURN
          ELSE
            DO I=1,N1
             MS=(I-1)*NR*N2+1
             NS=(I-1)*NB*N2+1
             CALL DGEMM('N','T',N2,NR,NB,alpha,W1(NS),N2,TMAT,NR,beta,
     $       W2(MS),N2)
            END DO
          END IF
	ELSE
          STOP 'OPTION IS WRONG IN DTRANS'
	END IF

	RETURN
	END 
c==========================================================================
c	function omp_get_thread_num
c
c	implicit none
c
c	integer omp_get_thread_num

c	omp_get_thread_num=0
c
c	end
c==========================================================================
c	function omp_get_max_threads
c	
c	implicit none
c
c	integer omp_get_max_threads
c
c	omp_get_max_threads=1
c
c	end
