        subroutine COUPLE14
        use COMPARAS, only : J1init, J2max, J12max, KB2max, NKB2,
     &    NBASJJ, NBASJJ14, parity5, NA2node, NB2node
        use ROTBASIS, only : IDXJJ14, KB2, CGC14

        implicit none
        integer :: j2Val, ik2, kb2Val, j12Val, ibasJJ14, m2Max,
     &    m2Val, ierr
        real*8  :: CG

        allocate(IDXJJ14(3,NBASJJ), stat=ierr)

        if (NA2node.eq.1.and.NB2node.eq.1) then
          NBASJJ14=1
          IDXJJ14(1,1)=0; IDXJJ14(2,1)=0; IDXJJ14(3,1)=1
          CGC14=1.d0
          return
        end if

        NBASJJ14=0
        do 10 j2Val=0, J2max
          do 10 ik2=1, NKB2
            kb2Val=KB2(ik2)
            if (j2Val.lt.abs(kb2Val)) goto 10
            
            do 20 j12Val=abs(J1init-j2Val), min(J12max, J1init+j2Val)
              NBASJJ14=NBASJJ14+1
              IDXJJ14(1,NBASJJ14)=j12Val
              IDXJJ14(2,NBASJJ14)=j2Val
              IDXJJ14(3,NBASJJ14)=ik2
              write(223,*) NBASJJ14, j12Val, j2Val, ik2 ,kb2Val
 20            continue
 10        continue

! CGC14
        allocate(CGC14(-J2max:J2max, NBASJJ14), stat=ierr)
        CGC14(:,:)=0.d0

        do ibasJJ14=1, NBASJJ14
          j12Val=IDXJJ14(1,ibasJJ14)
          j2Val=IDXJJ14(2,ibasJJ14)

          m2Max=min(J1init, j2Val)
          do m2Val=-m2Max,m2Max
            CGC14(m2Val,ibasJJ14)=CG(j2Val,m2Val,j12Val,0,J1init)
     .        *DSQRT((2*j12Val+1.0d0)/(2*J1init+1))
            write(225,*) m2Val, ibasJJ14, CGC14(m2Val, ibasJJ14)
          end do
        end do

        end subroutine COUPLE14
c=====================================================================
        subroutine VPOT14
        use QUADPOINTS, only : RQ, RUQ, AUQ, A2node, B2node 
        use ENGARRAYS, only : VMAT14
        use COMPARAS, only : NVB, NRB, N2D, NVBRU, NVBAU, NA2node, NB2node
        use SMALL_ARRAYS, only : VrefCH, Vref2D

        implicit none
        integer i, j, j2, k, l, ms,ns, ierr
        real*8 R, RU, AU, A2, B2, VRef, vx
        real*8 :: tmp

        allocate(VMAT14(NA2node*NB2node*N2D*NVB), stat=ierr)
        
        ms=0
        do i =1, NVB;  R=RQ(i)
        do j =1, NVBRU;  RU=RUQ(j)
        do j2=1, NVBAU; AU=AUQ(j2)
          VRef=VrefCH(NRB+i)+Vref2D((j-1)*NVBAU+j2)
        do k=1, NA2node;  A2=A2node(k)
        do l=1, NB2node;  B2=B2node(l)  
          ms=ms+1

          call REFPOT5(R,RU,AU,A2,B2,vx)

	   if(A2 > 1.2d0) vx=1.0d0
             VMAT14(ms)=vx-VRef
            !write(224,*) ms, VMAT14(ms)
        end do
        end do
        end do
        end do
        end do

	  tmp=180.D0/dacos(-1.D0)
	 open(100, file='VPOT14',form='formatted',status='unknown')
	   ms=0
	  do i=1, NVB
	  do j=1, N2D
	  do k=1, NA2node
 	  do l=1, NB2node
	     ms=ms+1
	     write(100,"(5f15.8)")RQ(i), AUQ(j)*tmp, A2node(k)*tmp, 
     &                            B2node(l)*tmp, VMAT14(ms)
	  end do
	  end do
	  end do
	  end do
	close(100)

	return
        end subroutine VPOT14
c======================================================================
        subroutine BOUND14(eigen_eng, eigen_vec)
        use COMPARAS, only : istat1, NVB, N2D, NBASJJ14
        implicit none

        integer    :: nbasis, nroots, memuse, maxcyc
	 real*8    :: tolvib, eigen_eng 
        complex*16 :: eigen_vec(1:NBASJJ14*N2D*NVB)   
        

        nbasis=NBASJJ14*N2D*NVB
        nroots=istat1*2
        memuse=istat1*4
        maxcyc=1000000
        tolvib=1.d-10       

!        print*, 'BOUND14'
!        print*, '============================================='
!        print*, 'nbasis=', nbasis
!        print*, 'nroots=', nroots
!        print*, 'memuse=', memuse
!        print*, 'maxcyc=', maxcyc
!        print*, 'tolvib=', tolvib
!        print*, '============================================='
!        print*
!
        rewind(21)
        call ZNARPACK(21,nbasis,nroots,memuse,maxcyc,tolvib,istat1,
     &    eigen_eng, eigen_vec) 
c        stop
c        call READWF14(nbasis,istat1,eigen_eng, eigen_vec)

        end subroutine BOUND14
c=====================================================================
        subroutine HMAT14(WF1, WF2)
        use TRANSMATS, only : TMATR, TMAT2D
        use ENGARRAYS, only : VMAT14
        use COMPARAS, only : NVB, NRB, N2D, NA2node, NB2node, NBASJJ14
        use SMALL_ARRAYS, only : EIG0, EIG2D, BROTR, BIA_CZ3, BIC_CZ3

        implicit none

        integer    :: i, j, k,  ms, ns, ierr, l
        real*8     :: tkval
	complex*16 :: WF1(1:NVB*N2D*NBASJJ14), WF2(1:NVB*N2D*NBASJJ14)   
	complex*16, allocatable :: WF3(:) 

        do i=1, NVB*N2D*NBASJJ14
          WF2(i)=WF1(i)
        end do

        allocate(WF3(NVB*N2D*NBASJJ14*30),stat=ierr)
	if(ierr /= 0) STOP 'WF3 error in HMAT14'

        call ZTRANS(NVB,NVB,1,N2D*NBASJJ14,WF2,WF3,TMATR,'BG')   
!	  WF3(:)=(0.D0,0.D0) 
!	   ms=0; ns=0 
!	  do i=1, NVB
!	  do k=1, NVB
!	  do j=1, N2D
!	     do l=1, NBASJJ14
!	     ms=(i-1)*N2D*NBASJJ14+(j-1)*NBASJJ14+l
!	     ns=(k-1)*N2D*NBASJJ14+(j-1)*NBASJJ14+l
!	     WF3(ms)=WF3(ms)+WF2(ns)*TMATR(i,k)  
!	     enddo
!	  enddo
!	  enddo
!	  enddo

        call ZTRANS(N2D,N2D,NVB,NBASJJ14,WF2,WF3,TMAT2D,'BG')   
!	  WF2(:)=(0.D0,0.D0) 
!	   ms=0; ns=0 
!	  do i=1, NVB
!	  do j=1, N2D
!	  do k=1, N2D
!	     do l=1, NBASJJ14
!	     ms=(i-1)*N2D*NBASJJ14+(j-1)*NBASJJ14+l
!	     ns=(i-1)*N2D*NBASJJ14+(k-1)*NBASJJ14+l
!	     WF2(ms)=WF2(ms)+WF3(ns)*TMAT2D(j,k)  
!	     enddo
!	  enddo
!	  enddo
!	  enddo

        ms=0
        do i=1,NVB
        do j=1,N2D
          ms=ms+1
          call RVMAT14(BROTR(NRB+i),BIA_CZ3(j),BIC_CZ3(j), 
     .      WF2((ms-1)*NBASJJ14+1),WF3,VMAT14((ms-1)*NA2node*NB2node+1))
        end do
        end do

        call ZTRANS(NVB,NVB,1,N2D*NBASJJ14,WF2,WF3,TMATR,'GB')   
!	  WF3(:)=(0.D0,0.D0) 
!	   ms=0; ns=0 
!	  do i=1, NVB
!	  do k=1, NVB
!	  do j=1, N2D
!	     do l=1, NBASJJ14
!	     ms=(i-1)*N2D*NBASJJ14+(j-1)*NBASJJ14+l
!	     ns=(k-1)*N2D*NBASJJ14+(j-1)*NBASJJ14+l
!	     WF3(ms)=WF3(ms)+WF2(ns)*TMATR(k,i)  
!	     enddo
!	  enddo
!	  enddo
!	  enddo
        call ZTRANS(N2D,N2D,NVB,NBASJJ14,WF2,WF3,TMAT2D,'GB')   
!	  WF2(:)=(0.D0,0.D0) 
!	   ms=0; ns=0 
!	  do i=1, NVB
!	  do j=1, N2D
!	  do k=1, N2D
!	     do l=1, NBASJJ14
!	     ms=(i-1)*N2D*NBASJJ14+(j-1)*NBASJJ14+l
!	     ns=(i-1)*N2D*NBASJJ14+(k-1)*NBASJJ14+l
!	     WF2(ms)=WF2(ms)+WF3(ns)*TMAT2D(k,j)  
!	     enddo
!	  enddo
!	  enddo
!	  enddo

! H=T+V
        ms=0
        do i=1, NVB
        do j=1, N2D
           tkVal=EIG0(i)+EIG2D(j)
        do k=1, NBASJJ14
          ms=ms+1
          WF2(ms)=WF2(ms)+WF1(ms)*tkVal
        end do
        end do
        end do

        deallocate(WF3, stat=ierr)

        end subroutine HMAT14
c======================================================================
        subroutine RVMAT14(BROTR,BIA_CZ3,BIC_CZ3,RWF2,RWF1,VMAT)
        use ROTBASIS, only : IDXJJ14, CGC14, DMATA2, CTMATB2, KB2 
        use COMPARAS, only : NA2node, NB2node, NBASJJ14, NKB2,
     &    J1init, J2max

	implicit none

	complex*16 :: RWF1(1), RWF2(1), wf3Val 
	 real*8   :: VMAT(NB2node,NA2node), krot    

	complex*16, allocatable :: WF3(:,:,:), WF4(:,:), WF5(:,:) 
        
	 real*8 BROTR, BIA_CZ3, BIC_CZ3    

        integer i, ia2, ib2, j12Val, j2Val, ik2, kb2Val, m2Val, m2Max, 
     &    ierr

        if (NBASJJ14.eq.1) then
          RWF2(1)=RWF2(1)*VMAT(1,1)
          return
        end if

        do i=1, NBASJJ14
           j12Val=IDXJJ14(1,i)
           j2Val=IDXJJ14(2,i)
           ik2=IDXJJ14(3,i); kb2Val=KB2(ik2)

           krot=j12Val*(j12Val+1)*BROTR+j2Val*(j2Val+1)*BIA_CZ3+
     .          kb2Val*kb2Val*(BIC_CZ3-BIA_CZ3)
           RWF1(i)=RWF2(i)*krot
        end do

c potential value
        allocate(WF3(1:NKB2, 0:J2max, -J2max:J2max), stat=ierr)
        allocate(WF4(NKB2, NA2node), stat=ierr)
        allocate(WF5(NB2node, NA2node), stat=ierr)

        WF3(:,:,:)=(0.D0,0.D0) 
        do  i=1, NBASJJ14
          j12Val=IDXJJ14(1,i)
          j2Val=IDXJJ14(2,i)
          ik2=IDXJJ14(3,i); kb2Val=KB2(ik2)
          m2Max=min(J1init, j2Val)
          
          do m2Val=-m2Max,m2Max
            WF3(ik2, j2Val, m2Val)=WF3(ik2, j2Val, m2Val)+
     &                             RWF2(i)*CGC14(m2Val, i)
	  enddo
	end do


        m2Max=min(J1init, J2max)
        do 20 m2Val=-m2Max, m2Max

          WF4(:,:)=(0.D0,0.D0)   
          do j2Val=0, J2max
          do ik2  =1, NKB2
            kb2Val=KB2(ik2)
          do ia2=1, NA2node
            WF4(ik2,ia2)=WF4(ik2,ia2)+WF3(ik2,j2Val,m2Val)*
     &                                DMATA2(ia2,kb2Val,j2Val,m2Val)
          end do
          end do
          end do

          WF5(:,:)=(0.D0,0.D0) 
          do ia2=1, NA2node
          do ik2=1, NKB2
          do ib2=1, NB2node
	   WF5(ib2,ia2)=WF5(ib2,ia2)+WF4(ik2,ia2)*CTMATB2(ib2,ik2)
          end do
          end do
          end do

          do ia2=1, NA2node
          do ib2=1, NB2node
	   WF5(ib2,ia2)=WF5(ib2,ia2)*VMAT(ib2,ia2)
          end do
          end do

          WF4(:,:)=(0.D0,0.D0) 
          do ia2=1, NA2node
	  do ik2=1, NKB2
          do ib2=1, NB2node
	   WF4(ik2,ia2)=WF4(ik2,ia2)+WF5(ib2,ia2)*CONJG(CTMATB2(ib2,ik2))
          end do
          end do
          end do

          do j2Val=0, J2max
          do ik2=1, NKB2
            wf3Val=(0.D0,0.D0) 
            kb2Val=KB2(ik2)
          do ia2=1, NA2node
             wf3Val=wf3Val+WF4(ik2,ia2)*DMATA2(ia2,kb2Val,j2Val,m2Val)
          end do
             WF3(ik2, j2Val, m2Val)=wf3Val
          end do
          end do

20      end do
	 

        RWF2(1:NBASJJ14)=(0.D0,0.D0)  !###
        do i=1, NBASJJ14
          j12Val=IDXJJ14(1,i)
          j2Val=IDXJJ14(2,i)
          ik2=IDXJJ14(3,i); kb2Val=KB2(ik2)

          m2Max=min(J1init, j2Val)
          
          do m2Val=-m2Max,m2Max
            RWF2(i)=RWF2(i)+WF3(ik2, j2Val, m2Val)*CGC14(m2Val,i)
          end do
        end do

        do i=1, NBASJJ14
          RWF2(i)=RWF2(i)+RWF1(i)
        end do
 
        deallocate(WF3,  WF4,  WF5,  stat=ierr)

        end subroutine RVMAT14
c=========================================================================
        subroutine READWF14(nbasis, istat1, eigen_eng, eigen_vec)
        implicit none
        integer i, nbasis, istat1, iroot, ierr
        real*4 eigen_vec(1), eigen_eng, eigen0, eigenx
        real*4, allocatable :: ev_tmp(:)
        real*8 wvnumb

        wvnumb=219474.6d0
        rewind(21)
        allocate(ev_tmp(nbasis),stat=ierr)
        do 2000 iroot=1,istat1
          if (iroot.ne.istat1) read(21) eigenx, ev_tmp
          if (iroot.eq.istat1) read(21) eigen_eng, eigen_vec

          if (iroot.eq.1) eigen0=eigenx
!          write(*, 9110)iroot,eigenx*wvnumb,(eigenx-eigen0)*wvnumb
 2000   continue
!        write(*,9910)

        deallocate(ev_tmp,stat=ierr)

c        do i=1, nbasis
c          write(226,*) i, abs(eigen_vec(i))**2.d0
c        end do

 9100        format(3x,'Eigenvalues for YCZ3 molecule')
 9110        format(3x,i5,20x,2f20.10)
 9910        format(3x,75('='))

        end subroutine READWF14
c======================================================================
        subroutine HMATDRIV2(WF1, WF2)
        implicit none
        complex*16 WF1(1), WF2(1) 

        call HMAT14(WF1, WF2)

        end subroutine HMATDRIV2
c======================================================================
