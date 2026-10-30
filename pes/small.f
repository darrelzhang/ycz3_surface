c================================================================
        SUBROUTINE NORCH_K(TIME,NA,NK,IDK,N1,W1,N2,W2)
        IMPLICIT REAL*8 (A-H,O-Z)
	COMPLEX*8 W1(1),W2(1)
	REAL*8 ANORM(NK,2)
	DIMENSION IDK(1)

	IF(TIME.GT.1.D-10) THEN
	  WRITE(*,"('  TIME = ', F10.2)") TIME
	END IF

	IF(NK.EQ.1) THEN
	  ANORM1=0.D0
	  DO K=1,NA*N1
            ANORM1=ANORM1+CONJG(W1(K))*W1(K)
          END DO

	ELSE
	  ANORM=0.D0
  	  DO I=1,N1
	  DO J=1,NK
	    NAk=IDK(J+1)-IDK(J)
	    MS=(I-1)*NA+IDK(J)-1
	  DO K=1,NAk
	    ANORM(J,1)=ANORM(J,1)+CONJG(W1(MS+K))*W1(MS+K)
	  END DO
	  END DO
	  END DO

	  ANORM1=0.D0
	  DO I=1,NK	
	    ANORM1=ANORM1+ANORM(I,1)
	  END DO
	END IF

	NP=NK/5
	IF(MOD(NK,5).NE.0) NP=NP+1

	IF(N2.EQ.0) THEN   ! Norm for one Vector
	   WRITE(*,"(2x,'NORM = ', E14.8)") ANORM1
	   IF(NK.GT.1) THEN
	     WRITE(*,"(2x,'K_NORM =', 5E12.4)") (ANORM(I,1),I=1,MIN(5,NK))
	     DO J=1,NP-1
	       NKleft=MIN(5,NK-5*J)
	       WRITE(*,"(10x, 5E12.4)") (ANORM(J*5+I,1),I=1,NKleft)
	     END DO
	   END IF
	ELSE
          IF(NK.EQ.1) THEN
            ANORM2=0.D0
            DO K=1,NA*N2
              ANORM2=ANORM2+CONJG(W2(K))*W2(K)
            END DO
          ELSE
            DO I=1,N2
            DO J=1,NK
              NAk=IDK(J+1)-IDK(J)
              MS=(I-1)*NA+IDK(J)-1
            DO K=1,NAk
              ANORM(J,2)=ANORM(J,2)+CONJG(W2(MS+K))*W2(MS+K)
            END DO
            END DO
            END DO

            ANORM2=0.D0
            DO I=1,NK     
              ANORM2=ANORM2+ANORM(I,2)
            END DO
          END IF

	  WRITE(*,"(2x,'NORM = ',E14.8,4x,'NORM1 = ',E12.5,4x,'NORM2 = ',E12.5)")
     $    ANORM1+ANORM2, ANORM1, ANORM2

	  IF(NK.GT.1) THEN
             WRITE(*,"(2x,'K_NORM1 =', 5E12.4)") (ANORM(I,1),I=1,MIN(5,NK))
             DO J=1,NP-1
               NKleft=MIN(5,NK-5*J)
               WRITE(*,"(11x, 5E12.4)") (ANORM(J*5+I,1),I=1,NKleft)
	     END DO

	     WRITE(*,"(2x,'K_NORM2 =', 5E12.4)") (ANORM(I,2),I=1,MIN(5,NK))
             DO J=1,NP-1
               NKleft=MIN(5,NK-5*J)
               WRITE(*,"(11x, 5E12.4)") (ANORM(J*5+I,2),I=1,NKleft)
             END DO

             WRITE(*,"(2x,'K_NORM  =', 5E12.4)") (ANORM(I,1)+ANORM(I,2),I=1,MIN(5,NK))
             DO J=1,NP-1
               NKleft=MIN(5,NK-5*J)
               WRITE(*,"(11x, 5E12.4)") (ANORM(J*5+I,1)+ANORM(J*5+I,2),I=1,NKleft)
             END DO
	  END IF
	END IF

	RETURN
	END
c========================================================================
        SUBROUTINE ZZERO1D(N,W0)
        COMPLEX*16 W0(N)

        DO I=1,N
          W0(I)=(0.D0,0.D0)
        END DO
        RETURN
        END
c==================================================================
	SUBROUTINE CZERO1D(N,W0)
	COMPLEX*8 W0(N)

	DO I=1,N
	  W0(I)=(0.D0,0.D0)
	END DO
	RETURN
	END 
c==================================================================
        SUBROUTINE DZERO1D(N,W0)
        REAL*8 W0(N)

        DO I=1,N
          W0(I)=0.D0
        END DO
        RETURN
	END
c==================================================================
        SUBROUTINE SZERO1D(N,W0)
        REAL*4 W0(N)

        DO I=1,N
          W0(I)=0.D0
        END DO
        RETURN
        END
c==================================================================
	SUBROUTINE CNORM(N,W0,AN)
	IMPLICIT REAL*8 (A-H,O-Z)
	COMPLEX*8 W0(1)

	AN=0.D0
	DO I=1,N
	  AN=AN+CONJG(W0(I))*W0(I)
	END DO

	RETURN
	END
c==================================================================
        SUBROUTINE ZNORM(N,W0,AN)
        IMPLICIT REAL*8 (A-H,O-Z)
        COMPLEX*16 W0(1)

        AN=0.D0
        DO I=1,N
          AN=AN+CONJG(W0(I))*W0(I)
        END DO

        RETURN
        END
c==================================================================
        SUBROUTINE SNORM(N,W0,AN)
        IMPLICIT REAL*8 (A-H,O-Z)
        REAL*4 W0(1)

        AN=0.D0
        DO I=1,N
          AN=AN+W0(I)*W0(I)
        END DO

        RETURN
        END
c==================================================================
        SUBROUTINE DNORM(N,W0,AN)
        IMPLICIT REAL*8 (A-H,O-Z)
        REAL*8 W0(1)

        AN=0.D0
        DO I=1,N
          AN=AN+W0(I)*W0(I)
        END DO

        RETURN
        END
c==================================================================
        SUBROUTINE NORCH(N,W0,AN)
        COMPLEX*8 W0(1)
        REAL*8 AN

        AN=0.D0
        DO I=1,N
          AN=AN+CONJG(W0(I))*W0(I)
        END DO

        RETURN
        END
