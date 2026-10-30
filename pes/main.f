c========================================================================
c      A program for time-dependent quantum dynamics calculation
c      of X+YCZ3->XY+CZ3 reaction using the model Hamilton of
c      Juliana Palma and David C. Clary, JCP 112(4) 1859, 2000
c      Iparity = parity * (-1)**Jtotal
c========================================================================
        program X_YCZ3
        use COMPARAS, only : NVBRU
	
	implicit none

	include 'mpif.h'

	integer :: ierr

	call MPI_INIT(ierr)

!        print*, 'READINP'; print*
        call READINP;          

!        print*, 'CHKINP';  print*
        call CHKINP;    

!	call PES_INIT_BM()

!	call MPI_Barrier(MPI_COMM_WORLD, ierr)

!        print*, 'POTZERO'; print*
        call POTZERO;        

!	 STOP
!        print*, 'BASCOUP'; print*
        call BASCOUP;        

!        print*, 'ZBASIS';  print*
        call ZBASIS;    

!        print*, 'R1BASIS'; print*
        call R1BASIS;   

        if (NVBRU.eq.1) then
!          print*, 'AUBASIS'; print*
          call AUBASIS;   
        else
!          print*, 'BOUND2D'; print*
          call BOUND2D
        end if

!        print*, 'A1BASIS'; print*
        call A1BASIS;        

!        print*, 'A2BASIS'; print*
        call A2BASIS;        

!        print*, 'B1BASIS'; print*
        call B1BASIS;        

!        print*, 'B2BASIS'; print*
        call B2BASIS;        

	call B3BASIS; 

!        print*, 'COEFFCG'; print*
        call COEFFCG;   

!	 STOP 'check centri potential'
!        print*, 'WFABS';   print*
        call WFABS;     

!        print*, 'FFLUX';   print*
        call FFLUX;     

!        print*, 'CPLXENG'; print*
        call CPLXENG;   

!        print*, 'POTENTIALS'; print*
         call POTENTIALS;

!        print*, 'WFINIT';  print*
        call WFINIT;    

!        print*, 'START OF REACTION'; print*
!        call REACTION;      
!        print*, 'END OF REACTION'; print*

	call MPI_FINALIZE(ierr)

        end program X_YCZ3
