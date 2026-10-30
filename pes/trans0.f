
	SUBROUTINE ZTRANS(NB,NR,N1,N2,W1,W2,TMAT,OPTION)
	IMPLICIT REAL*8 (A-H,O-Z)
	PARAMETER (NCmax=400)   ! maximum chunk size for a job
	CHARACTER*2 OPTION

	Mmax=MAX(NB,N1,N2)
	N22=N2*2
	MNRB=MAX(NB,NR)
	IF(N1.GE.4 .AND. N22.LE.NCmax .AND. NR.EQ.NB) THEN
	  CALL DTRANS1(NB,NR,N1,N22,W1,W2,TMAT,OPTION)  ! parallelize loop N1
	ELSE
	  NC=MIN(N22,NCmax)
	  CALL DTRANS2(MNRB,NB,NR,N1,N22,W1,W2,TMAT,OPTION,NC)  ! parallelize loop N2
	END IF
	RETURN
	END
c=======================================================================	

	SUBROUTINE ZTRANSS(NB,NR,N1,N2,W1,W2,TMAT,OPTION,W3,W4)
	IMPLICIT REAL*8 (A-H,O-Z)
	PARAMETER (NCmax=400)   ! maximum chunk size for a job
	CHARACTER*2 OPTION
	real*8,allocatable:: WTEMP(:)
	COMPLEX*16 W1(1),W2(1),W3(1),W4(1)
	Mmax=MAX(NB,N1,N2)
c	allocate(WTEMP(NB*N1*N2*4),stat=ierr)
	N22=N2*2
	MNRB=MAX(NB,NR)
	IF(N1.GE.4 .AND. N22.LE.NCmax .AND. NR.EQ.NB) THEN
	  CALL DTRANS1(NB,NR,N1,N22,W1,W2,TMAT,OPTION)  ! parallelize loop N1
	ELSE
	  NC=MIN(N22,NCmax)
	  CALL DTRANS2(MNRB,NB,NR,N1,N22,W1,W2,TMAT,OPTION,NC)  ! parallelize loop N2
	END IF
c	deallocate(WTEMP)
	DO I=1,NB*N1*N2
	W2(I)=W1(I)
	END DO

	RETURN
	END
c=======================================================================
	SUBROUTINE DTRANS1(NB,NR,N1,N2,W1,W2,TMAT,OPTION)
c
c W1 --> W1 with W2 a small working space by parallelizing N1 loop
c
	IMPLICIT REAL*8 (A-H,O-Z)
	CHARACTER*2 OPTION
	REAL*8 TMAT(NB,1),alpha,beta,W1(1),W2(N2*NB,1)
	integer omp_get_thread_num, omp_get_max_threads

	NP=omp_get_max_threads()    

C$OMP PARALLEL SHARED (NR,NB,N1,N2,TMAT,W1,W2,NP,OPTION)
C$OMP&  PRIVATE (NSTART,NEND,N22,N,IP)
C$OMP DO
        DO N=1,NP
	  IP=omp_get_thread_num()+1
          NSTART=(((N-1.D0)*N1)/NP+1.D-9)
          NEND=((1.D0*N*N1)/NP+1.D-9)
          N22=NEND-NSTART
	  CALL DTRANS0(NR,NB,N22,N2,W1(NSTART*NR*N2+1),W2(1,IP),TMAT,OPTION)
	END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL

	RETURN
	END 
c=======================================================================
	SUBROUTINE DTRANS0(NB,NR,N1,N2,W1,W2,TMAT,OPTION)
c ========================================
c W1 --> W1  with W2 a small working space
c ========================================
	IMPLICIT REAL*8 (A-H,O-Z)
	CHARACTER*2 OPTION
	REAL*8 TMAT(NB,1),alpha,beta,W1(1),W2(NB*N2)

	alpha=1.D0
	beta=0.D0

	IF(OPTION .EQ. 'BG') THEN   ! GRID ---> BASIS
	  DO I=1,N1
	    MS=(I-1)*NB*N2+1
	    CALL DGEMM('N','N',N2,NR,NB,alpha,W1(MS),N2,TMAT,
     $      NB,beta,W2,N2)
	
	    INC=1
	    MS=(I-1)*NR*N2+1
	    DO J=1,NR*N2
	      W1(MS+J-1)=W2(J)
	    END DO
c	    CALL DCOPY(NR*N2,W2,INC,W1(MS),INC)
	  END DO
	ELSE IF (OPTION .EQ. 'GB') THEN    ! BASIS  --->  GRID
          DO I=1,N1
            MS=(I-1)*NR*N2+1
            CALL DGEMM('N','T',N2,NB,NR,alpha,W1(MS),N2,TMAT,
     $      NB,beta,W2,N2)

	    INC=1
	    MS=(I-1)*NB*N2+1
c	    CALL DCOPY(NB*N2,W2,INC,W1(MS),INC)
	    DO J=1,NB*N2
	     W1(MS+J-1)=W2(J)
	    END DO
         END DO
	END IF

	RETURN
	END 
c=======================================================================
	SUBROUTINE DTRANS2(MNRB,NB,NR,N1,N2,W1,W2,TMAT,OPTION,NC)
c ========================================
c W1 --> W1  with W2 a small working space
c ========================================
	IMPLICIT REAL*8 (A-H,O-Z)
	CHARACTER*2 OPTION
	REAL*8 TMAT(NB,1),alpha,beta,W1(1),W2(NC*MNRB,1)
	integer omp_get_thread_num, omp_get_max_threads

	N2D=N2
	alpha=1.D0
	beta=0.D0

	NP=omp_get_max_threads()    ! Minimum number of 
	MM=N2D/NP
	IF(MOD(N2D,NP).NE.0) MM=MM+1
	IF(MM.GT.NC) THEN
          MM=NC
          NP=(N2D/(MM*NP)+1)*NP
        END IF

	IF(OPTION .EQ. 'BG') THEN   ! GRID ---> BASIS
                                    ! sum over the first index.

	  IF(N2.EQ.1) THEN
	     PRINT*,'YOU SHOULD USE ZTRANS1'
	     STOP
	  ELSE

c rearrange W1
	    IF(NR.GT.NB) THEN
	      INC=-1
	      DO I=N1,1,-1
	      DO J=NB,1,-1
	        MS=(I-1)*NB*N2D+(J-1)*N2D
	        NS=(I-1)*NR*N2D+(J-1)*N2D
                CALL DCOPY(N2D,W1(MS+1),INC,W1(NS+1),INC)
	      END DO
	      END DO
	    END IF
	
	    NBR=MAX(NB,NR)
C$OMP PARALLEL SHARED (NBR,NR,NB,N1,N2D,TMAT,W1,W2,NP,alpha,beta,MM)
C$OMP&  PRIVATE (I,MS,NS,NSTART,NEND,N22,N,IP,NF,J,K)
C$OMP DO
            DO N=1,NP
	      IP=omp_get_thread_num()+1
              NSTART=(((N-1.D0)*N2D)/NP+1.D-9)+1
	      NEND=((1.D0*N*N2D)/NP+1.D-9)
              N22=NEND-NSTART+1

	    DO I=1,N1
	     MS=(I-1)*NBR*N2D
	     CALL DGEMM('N','N',N22,NR,NB,alpha,W1(MS+NSTART),N2D,TMAT,
     $       NB,beta,W2(1,IP),MM)
	
	     DO J=1,NR
	       MS=(I-1)*NBR*N2D
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
                CALL DCOPY(N2D,W1(MS+1),INC,W1(NS+1),INC)
              END DO
              END DO
            END IF

	  END IF
	ELSE IF (OPTION .EQ. 'GB') THEN    ! BASIS  --->  GRID
          IF(N2.EQ.1) THEN
	    PRINT*,'YOU SHOULD USE CTRANS1'
	    STOP
          ELSE

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
              NSTART=(((N-1.D0)*N2D)/NP+1.D-9)+1
	      NEND=((1.D0*N*N2D)/NP+1.D-9)
              N22=NEND-NSTART+1

            DO I=1,N1
             MS=(I-1)*NBR*N2D
             CALL DGEMM('N','T',N22,NB,NR,alpha,W1(MS+NSTART),N2D,TMAT,
     $       NB,beta,W2(1,IP),MM)

	     DO J=1,NB
               MS=(I-1)*NBR*N2D
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
