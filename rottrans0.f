c====================================================================
c	subroutine VMATRIX(npoten,LINK2,VCUT0,VM,W1)
	subroutine VMATRIX(npoten,LINK2,VCUT0,VM,W1,WT1,WT2,WT3)
        use COMPARAS, only : NBASJJ,NA1node,NA2node,NB1node,NB2node,
     .    KBZmin,KBZmax
        use ROTBASIS, only : IDKBZ
        implicit none
        integer npoten, LINK2(NA1node*NA2node),NA,NQA,ierr
        integer kbzVal
        real*4 VCUT0, VM(NB1node*NB2node,npoten)
        complex*8 W1(NBASJJ)
c        complex*8, allocatable :: WT1(:), WT2(:), WT3(:)
        complex*8 WT1(2*NA1node*NA2node*NB1node*NB2node) 
	complex*8 WT2(2*NA1node*NA2node*NB1node*NB2node)
	complex*8 WT3(2*NA1node*NA2node*NB1node*NB2node)

	NQA=NA1node*NA2node*NB1node*NB2node

c	allocate(WT1(NQA*2), WT2(NQA*2), WT3(NQA*2), stat=ierr)

	WT1(1:2*NQA)=0.0
	WT2(1:2*NQA)=0.0
	WT3(1:2*NQA)=0.0
        do kbzVal=KBZmin, KBZmax
          NA=IDKBZ(2,kbzVal)
          if (kbzVal.eq.0) then
            call ROTTRANS0(npoten,kbzVal,NA,LINK2,VCUT0,W1(IDKBZ(1,kbzVal)),
     .        VM,WT1,WT2,WT3)
          else
c            call ROTTRANS(kbzVal,NA,LINK2,VCUT0,W1,VM)
          end if
        end do

        end subroutine VMATRIX
c====================================================================
! kbzVal=0 case, a basic version
c====================================================================
        subroutine ROTTRANS0(npoten,kbzVal,NA,LINK2,VCUT0,W1,VM,WT1,WT2,WT3)
        use ROTBASIS
        use COMPARAS
        implicit none
        integer LINK2(NA1node*NA2node), npoten
        integer kbzVal, j1Val, j2Val, j12Val, kb2Val, istart, NA, ierr
        integer i, m1max0, m1, ik2, ik1, ia1, ia2, ib1, ib2, j1Min, j2Min

        real*4 VM(NB1node*NB2node,npoten), VCUT0
        complex*8 W1(NA), CI, CTMP, CTMP1, CTMP2
        complex*8 WT1(2*NA1node*NA2node*NB1node*NB2node)
	complex*8 WT2(2*NA1node*NA2node*NB1node*NB2node)
	complex*8 WT3(2*NA1node*NA2node*NB1node*NB2node)

! J12->m1
        call J12M1TRANS(1,kbzVal,NA,W1,WT1)

! J1 -> NA1node
        call J1A1TRANS(1,kbzVal,WT1,WT2)

! J2 -> NA2node
        call J2A2TRANS(1,WT2,WT1)

! (m1, m2) <-> (NB1node, NB2node) and potential action
        call BBVT0(npoten,VCUT0,LINK2,RTMATB1,RTMATB2,VM,WT1,WT2,WT3)

C NA2node -> J2
        call J2A2TRANS(-1,WT2,WT1)

! NA1node -> J1
        call J1A1TRANS(-1,kbzVal,WT1,WT2)

c M1->J12
        call J12M1TRANS(-1,kbzVal,NA,W1,WT1)

        end subroutine ROTTRANS0
c====================================================================
        subroutine BBVT0(npoten,VCUT0,LINK2,RT1,RT2,VM,WT,WT2,WT3)
        use COMPARAS, only : NA1node, NA2node, NB1node, NB2node, M1max,
     .    NKB20, DT

        implicit none
        integer LINK2(NA1node*NA2node),npoten
        real*4 VCUT0,VM(NB2node,NB1node,npoten)
        real*4 RT1(NB1node,2*M1max+1,2), RT2(NB2node,NKB20,2)
        complex*8 WT(NKB20,2*M1max+1,NA1node*NA2node),CTMP0,CTMP,CTMP1,CTMP2

        complex*8 WT2(NKB20,NB1node,2), WT3(NB2node,NB1node)

        integer nodeax, ik1, ib1, ib2, ik2, iax, ipoten, ierr, INC, m1

        nodeax=NA1node*NA2node
        CTMP0=cdexp((0.d0,-1.d0)*VCUT0*DT)
        INC=1
        do iax=1, nodeax
          ipoten=LINK2(iax)

          if (ipoten.eq.0) then

c            call CSCAL(NKB20*(2*M1max+1),CTMP0,WT(1,1,iax),INC)
	     do ik1=1, 2*M1max+1
	     do ik2=1, NKB20
               WT(ik2,ik1,iax)=WT(ik2,ik1,iax)*CTMP0
             end do
             end do

          else

! (m1->b1)
c            do ik1=1, 2*M1max+1
c            do ib1=1, NB1node
c            do ik2=1, NKB20
c              WT2(ik2,ib1)=WT2(ik2,ib1)+WT(ik2,ik1,iax)*RT1(ib1,ik1,1)
c              WT3(ik2,ib1)=WT3(ik2,ib1)+WT(ik2,ik1,iax)*RT1(ib1,ik1,2)
c            end do
c            end do
c            end do


c	     WT2(:,:,:)=0.0
             call SGEMM('N','T',2*NKB20,NB1node,2*M1max+1,
     &                          1.0e+0,WT(1,1,iax),2*NKB20,
     &                          RT1(1,1,1), NB1node,
     &                          0.0e+0,WT2(1,1,1), 2*NKB20)

             call SGEMM('N','T',2*NKB20,NB1node,2*M1max+1,
     &                          1.0e+0,WT(1,1,iax),2*NKB20,
     &                          RT1(1,1,2), NB1node,
     &                          0.0e+0,WT2(1,1,2), 2*NKB20)

! (m2->b2)
            do ib1=1, NB1node
            do ib2=1, NB2node
              CTMP=0.0
              do ik2=1, NKB20
                CTMP=CTMP+RT2(ib2,ik2,1)*WT2(ik2,ib1,1)+RT2(ib2,ik2,2)*WT2(ik2,ib1,2)
              end do
              WT3(ib2,ib1)=CTMP
            end do
            end do
            
! exp(-i*V*dt)
            do ib1=1, NB1node
            do ib2=1, NB2node
              CTMP=CDEXP((0.d0,-1.d0)*DT*VM(ib2,ib1,ipoten))
              WT3(ib2,ib1)=WT3(ib2,ib1)*CTMP
            end do
            end do                

! (b2->m2)
            do ib1=1, NB1node
            do ik2=1, NKB20
              CTMP1=0.0
              CTMP2=0.0
              do ib2=1, NB2node
                CTMP1=CTMP1+RT2(ib2,ik2,1)*WT3(ib2,ib1)
                CTMP2=CTMP2+RT2(ib2,ik2,2)*WT3(ib2,ib1)
              end do
              WT2(ik2,ib1,1)=CTMP1
              WT2(ik2,ib1,2)=CTMP2
            end do
            end do

!(b1->m1)
c            WT(1:NKB20,1:2*M1max+1,iax)=0.0
c            do ik1=1, 2*M1max+1
c            do ib1=1, NB1node
c            do ik2=1, NKB20
c              WT(ik2,ik1,iax)=WT(ik2,ik1,iax)+WT2(ik2,ib1)*RT1(ib1,ik1,1)+WT3(ik2,ib1)*RT1(ib1,ik1,2)
c            end do
c            end do
c            end do

            call SGEMM('N','N',2*NKB20,2*M1max+1,NB1node,
     &        1.0e+0,WT2(1,1,1),NKB20*2,RT1(1,1,1),NB1node,
     &        0.0e+0,WT(1,1,iax),NKB20*2)
            call SGEMM('N','N',2*NKB20,2*M1max+1,NB1node,
     &        1.0e+0,WT2(1,1,2),NKB20*2,RT1(1,1,2),NB1node,
     &        1.0e+0,WT(1,1,iax),NKB20*2)

            end if

        end do 

        end subroutine BBVT0
c====================================================================
        subroutine J12M1TRANS(IFLAG,kbzVal,NA,W1,WT1)
        use ROTBASIS
        use COMPARAS
        implicit none
        integer kbzVal, j1Val, j2Val, j12Val, kb2Val, istart, NA
        integer i, m1max0, m1, ik2, ik1, ia1, ia2, ib1, ib2, ms, j1Min, j2Min

        integer IFLAG

        complex*8 W1(NA)
        complex*8 WT1(0:J1max,0:J2max,-M1max:M1max,1:NKB20)

! J12->m1

        if (IFLAG.eq.1) then

        WT1(0:J1max,0:J2max,-M1max:M1max,1:NKB20)=0.0
        do i=1, NA
          istart=IDKBZ(1,kbzVal)
          j1Val=IDXJJ(1,istart+i-1)
          j2Val=IDXJJ(2,istart+i-1)
          j12Val=IDXJJ(3,istart+i-1)
          ik2=IDXJJ(4,istart+i-1); kb2Val=KB2(ik2)

          m1max0=min(M1max,j1Val,j2Val)

          do m1=-m1max0, m1max0
            WT1(j1Val,j2Val,m1,ik2)=WT1(j1Val,j2Val,m1,ik2)+W1(i)*CGC(m1,istart+i-1)
          end do
        end do

        else

c M1->J12
        W1(1:NA)=0.0
        do i=1, NA
          istart=IDKBZ(1,kbzVal)
          j1Val=IDXJJ(1,istart+i-1)
          j2Val=IDXJJ(2,istart+i-1)
          j12Val=IDXJJ(3,istart+i-1)
          ik2=IDXJJ(4,istart+i-1); kb2Val=KB2(ik2)

          m1max0=min(M1max,j1Val,j2Val)
          do m1=-m1max0, m1max0
            W1(i)=W1(i)+WT1(j1Val,j2Val,m1,ik2)*CGC(m1,istart+i-1)
          end do
        end do

        end if

        end subroutine J12M1TRANS
c====================================================================
        subroutine J1A1TRANS(IFLAG,kbzVal,WT1,WT2)
        use ROTBASIS
        use COMPARAS
        implicit none

        integer kbzVal, j1Val, j2Val, j12Val, kb2Val, istart, NA, ierr
        integer i, m1max0, m1, ik2, ik1, ia1, ia2, ib1, ib2, ms, j1Min, j2Min
        integer IFLAG

        complex*8 WT1(0:J1max,0:J2max,-M1max:M1max,1:NKB20)
        complex*8 WT2(1:NA1node,0:J2max,-M1max:M1max,1:NKB20)

        if (IFLAG.eq.1) then

! J1 -> NA1node
        WT2(1:NA1node,0:J2max,-M1max:M1max,1:NKB20)=0.0
        do ik2=1, NKB20
          kb2Val=KB2(ik2)
          do m1=-M1max, M1max
            j1Min=max(abs(m1), kbzVal)
            j2Min=max(abs(m1), abs(kb2Val))
            do j2Val=j2Min, J2max
              do ia1=1, NA1node
                do j1Val=j1Min, J1max
                  WT2(ia1,j2Val,m1,ik2)=WT2(ia1,j2Val,m1,ik2) +
     .              WT1(j1Val,j2Val,m1,ik2)*DMATA1(ia1,m1,j1Val,kbzVal)
                end do
              end do
            end do
          end do
        end do

        else

c A1->J1
        WT1(1:NA1node,0:J2max,-M1max:M1max,1:NKB20)=0.0
        do ik2=1, NKB20
          kb2Val=KB2(ik2)
          do m1=-M1max, M1max
            j1Min=max(abs(m1), kbzVal)
            j2Min=max(abs(m1), abs(kb2Val))
            do j2Val=j2Min, J2max
              do j1Val=j1Min, J1max
                do ia1=1, NA1node
                  WT1(j1Val,j2Val,m1,ik2)=WT1(j1Val,j2Val,m1,ik2) + 
     .              WT2(ia1,j2Val,m1,ik2)*DMATA1(ia1,m1,j1Val,kbzVal)
                end do
              end do
            end do
          end do
        end do

        end if

        end subroutine J1A1TRANS
c====================================================================
        subroutine J2A2TRANS(IFLAG,WT2,WT1)
        use ROTBASIS
        use COMPARAS
        implicit none

        integer kbzVal, j1Val, j2Val, j12Val, kb2Val, istart, NA, ierr
        integer i, m1max0, m1, ik2, ik1, ia1, ia2, ib1, ib2, ms, j1Min, j2Min
        integer IFLAG
        complex*8 WT1(1:NKB20,-M1max:M1max,1:NA2node,1:NA1node)
        complex*8 WT2(1:NA1node,0:J2max,-M1max:M1max,1:NKB20)
        complex*8 CTMP

        if (IFLAG.eq.1) then

c J2->NA2node
        do ik2=1, NKB20
          kb2Val=KB2(ik2)
          do m1=-M1max, M1max
            j2Min=max(abs(m1), abs(kb2Val))
            do ia2=1, NA2node
              do ia1=1, NA1node
                CTMP=0.d0
                do j2Val=j2Min, J2max
                  CTMP=CTMP+WT2(ia1,j2Val,m1,ik2)*DMATA2(ia2,kb2Val,j2Val,m1)
                end do
                WT1(ik2,m1,ia2,ia1)=CTMP
              end do
            end do
          end do
        end do

        else

c A2->J2
        do ik2=1, NKB20
          kb2Val=KB2(ik2)
          do m1=-M1max, M1max
            j2Min=max(abs(m1), abs(kb2Val))
            do j2Val=j2Min, J2max
              do ia1=1, NA1node
                CTMP=0.d0
                do ia2=1, NA2node
                  CTMP=CTMP+WT1(ik2,m1,ia2,ia1)*DMATA2(ia2,kb2Val,j2Val,m1)
                end do
                WT2(ia1,j2Val,m1,ik2)=CTMP
              end do
            end do
          end do
        end do

        end if

        end subroutine J2A2TRANS
c====================================================================
