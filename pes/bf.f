c=========================================================================
! to calculate rotation coupling basis number
! Jtotal: the total angular momentum number of XYCZ3, "good" quantum number
! J1: the total angular momentum number of YCZ3
! J2: the rotational angular momentum of CZ3
! J12: the orbital rotation number of Y with respect to CZ3
! KBZ: the projection of Jtotal on BF-Z axis, it is also the projection
! of J1 on BF-Z axis
! KB2: the projection of J2 on CZ3 rotation axis
! 0 <= KBZ <= min(Jtotal, J1)
! |J2-J12| <= J1 <= J2+J12 => |J1-J2| <= J12 <= min(J12max, J1+J2)
! after using parity-adpated rotation basis,
! if KBZ = 0, KB2 >= 0 
! if KBZ=0, and KB2 = 0, Iparity*(-1)**(J1+J12+J2+KB2) = 1
! the rotation basis coupling is doing in this way:
! 1) count the KB2 
! 2) set up (J2,KB2) pair
! 3) set up (J1,J2,J12,KB2) pair
! 4) set up TABJJ table
c=========================================================================
	subroutine BASCOUP
	use COMPARAS, only :  KBZmax,  J12max,
     .    JS1, JD1, J1max, JS2, JD2, J2max, KB2max, NKBZ, NKB2, NKB20, 
     .    NBASJJ, NBASCZ3, DT  ,J1init, nfoldb2
	use ROTBASIS, only : IDKBZ, IDXCZ3, IDXJJ, TABJJ, KB2
	use ENGARRAYS, only : KROTJJ
	implicit none

	integer mKB2, kbzVal, kb2Val, mBASCZ3, j1Val, j2Val, j12Val,
     &    ik2, mBASJJ, ibasjj, ierr, iBASCZ3, Iparity, KB2min
	integer, allocatable :: IDX_J1KBZ(:,:)
	integer  :: NBAS_J1KBZ, iBAS_J1KBZ 

	NKBZ=2*KBZmax+1  !###

	mKB2=(2*KB2max+1)
	allocate(KB2(mKB2), stat=ierr)

	KB2min=-KB2max
!	if (KBZmax.gt.0) KB2min=-KB2max    ! ???
	
c for ortho isotopomer, mod(kb2Val, 3) = 0
c for para isotopomer, mod(kb2Val, 3) = 1,2
	NKB2=0
	NKB20=0
	do 10 kb2Val=KB2max, KB2min, -nfoldb2
!	  if(mod(kb2Val,3) .EQ. 0) cycle  !###
	  NKB2=NKB2+1
	  if (kb2Val.ge.0) NKB20=NKB20+1
	
	  KB2(NKB2)=kb2Val
c	  write(300,*) NKB2, KB2(NKB2)
10	continue

c to set up (j2, kb2) pair
	mBASCZ3=NKB2*(J2max+1)
	NBASCZ3=0
	allocate(IDXCZ3(2,mBASCZ3), stat=ierr)
	do 20 j2Val=0, J2max
	  do 30 ik2=1, NKB2
	    kb2Val=KB2(ik2)
            if (abs(kb2Val).gt.j2Val) goto 30
	  
	    NBASCZ3=NBASCZ3+1
	    IDXCZ3(1,NBASCZ3)=j2Val
	    IDXCZ3(2,NBASCZ3)=ik2
c	    write(301,*) NBASCZ3, IDXCZ3(1,NBASCZ3), kb2Val
30	  continue  
20	continue
!================================================
	 allocate(IDX_J1KBZ(2,(2*KBZmax+1)*(J1max+1)),stat=ierr)
	 NBAS_J1KBZ=0; IDX_J1KBZ(:,:)=0
	 do J1val=0, J1max
	 do kbzval=-KBZmax, KBZmax
	    if(ABS(kbzval).GT.J1val) cycle
	    NBAS_J1KBZ=NBAS_J1KBZ+1
	    IDX_J1KBZ(1,NBAS_J1KBZ)=J1val
	    IDX_J1KBZ(2,NBAS_J1KBZ)=kbzval
	 enddo
	 enddo 
!================================================

c to set up (j1, j2, j12, kb2) pair
	mBASJJ=mKB2*(J1max+1)*(J2max+1)*(J12max+1)*NKBZ
	allocate(IDXJJ(5,mBASJJ),stat=ierr)
	NBASJJ=0; IDXJJ(:,:)=0
	 do iBAS_J1KBZ=1, NBAS_J1KBZ
	     j1Val =IDX_J1KBZ(1,iBAS_J1KBZ)
	     kbzVal=IDX_J1KBZ(2,iBAS_J1KBZ)  

	    do  iBASCZ3=1, NBASCZ3
	      j2Val=IDXCZ3(1,iBASCZ3)
	      ik2=IDXCZ3(2,iBASCZ3)
	      kb2Val=KB2(ik2)
	
	      do j12Val=abs(j1Val-j2Val),min((j1Val+j2Val),J12max)

	        NBASJJ=NBASJJ+1
	        IDXJJ(1,NBASJJ)=j1Val
	        IDXJJ(2,NBASJJ)=j2Val
	        IDXJJ(3,NBASJJ)=j12Val
	        IDXJJ(4,NBASJJ)=ik2
	        IDXJJ(5,NBASJJ)=kbzVal

c	        write(302,*) NBASJJ,IDXJJ(1,NBASJJ),IDXJJ(2,NBASJJ),IDXJJ(3,NBASJJ),IDXJJ(4,NBASJJ)

	      end do !j12Val
	    end do !iBASCZ3
	 end do  !iBAS_J1KBZ 

c to set up TABJJ which is used for calculating (Jtotal-J)^2
	allocate(TABJJ(0:J12max,0:J2max,0:J1max,-KB2max:KB2max,-KBZmax:KBZmax)
     &           ,stat=ierr)
	TABJJ(:,:,:,:,:)=0
	
	  do  ibasjj=1, NBASJJ
	    j1Val =IDXJJ(1,ibasjj)
	    j2Val =IDXJJ(2,ibasjj)
	    j12Val=IDXJJ(3,ibasjj)
	    ik2   =IDXJJ(4,ibasjj); kb2Val=KB2(ik2)
	    kbzVal=IDXJJ(5,ibasjj)

	    TABJJ(j12Val, j2Val, j1Val, kb2Val, kbzVal)=ibasjj
	  end do 

! rotation kinetic energy
	allocate(KROTJJ(3,NBASJJ), stat=ierr)
	do ibasjj=1,NBASJJ
	  j2Val=IDXJJ(2,ibasjj); j12Val=IDXJJ(3,ibasjj)
	  ik2=IDXJJ(4,ibasjj); kb2Val=KB2(ik2)
	  
	  KROTJJ(1,ibasjj)=-j12Val*(j12Val+1.0)*DT/2.0
	  KROTJJ(2,ibasjj)=-j2Val*(j2Val+1.0)*DT/2.0
	  KROTJJ(3,ibasjj)=-kb2Val*kb2Val*DT/2.0
	  
	end do

!        print*, "The parameters of j-j coupling"
!        print*, "====================="
!
!        write(*,*)
!	print*, 'NKBZ=',NKBZ,'NKB2=',NKB2,'NKB20=',NKB20
!	print*, 'NBASCZ3=',NBASCZ3,'NBASJJ=',NBASJJ
!        print*, "====================="
!	write(*,*)

	end subroutine BASCOUP
c============================================================================
	subroutine COEFFCG
	use ROTBASIS, only : CGC, IDXJJ, IDKBZ
	use COMPARAS, only : NBASJJ, J1max, KBZmax
	implicit none
	integer kbzVal, ibasjj, j1Val, j2Val, j12Val, m1Val, ierr
	real*8 CG

	allocate(CGC(-J1max:J1max,NBASJJ), stat=ierr)
	CGC(:,:)=0.D0

	  do  ibasjj = 1, NBASJJ
	    j1Val=IDXJJ(1,ibasjj)
	    j2Val=IDXJJ(2,ibasjj)
	    j12Val=IDXJJ(3,ibasjj)

	    do  m1Val=-j1Val, j1Val
	      CGC(m1Val,ibasjj)=CG(j2Val,m1Val,j12Val,0,j1Val)*dsqrt
     &          ((2*j12Val+1.0d0)/(2*j1Val+1))
c	      write(211,*) j2Val, m1Val, j12Val, j1Val, CGC(m1Val,ibasjj)
	     end do
	  end do

	return
	end subroutine COEFFCG
c============================================================================
       subroutine A1BASIS
       use QUADPOINTS, only : A1node, A1weight
       use ROTBASIS, only : DMATA1
       use COMPARAS, only : NA1node, J1max, KBZmax
       implicit none
       integer kind, kpts, i, j, k, l, ms, ierr
       real*8 endpts(2),WK(NA1node)
       real*8 alpha, beta
       real*8, allocatable :: TMP(:,:,:,:)

       allocate(A1node(NA1node), A1weight(NA1node), stat=ierr)
       allocate(DMATA1(NA1node,-J1max:J1max,0:J1max,-KBZmax:KBZmax),stat=ierr)
       allocate(TMP(NA1node,-J1max:J1max,0:J1max,-KBZmax:KBZmax),stat=ierr)

       kind=1; kpts=0; alpha=0.d0; beta=0.d0
       call GAUSSQ(kind,NA1node,alpha,beta,kpts,endpts,WK,A1node,A1weight)
       
       do i = 1, NA1node
         A1node(i)=dacos(A1node(i))
       end do

       call ROTMAT1(NA1node,NA1node,J1max,-KBZmax,KBZmax,-J1max,J1max,
     .   A1node,A1weight,TMP)
	
       DMATA1=TMP

!	print*, "Parameters for A1 basis"
!	print*, "=================================="
!	print*
!	print*, "Points and Weights for A1"
!	!open(10,file="aa1",status="UNKNOWN")
!	do i=1, NA1node
!	  write(*,*) i, A1node(i), A1weight(i)
!	  !write(10,*) A1node(i)
!	end do
!	!close(10)
!	write(*,*)
!	print*, "=================================="

c	do l=-KBZmax,KBZmax
c	do k=0,J1max
c	do j=-J1max,J1max
c	do i=1,NA1node
c	  write(205,*) i, j, k, l, TMP(i,j,k,l)
c	end do
c	end do
c	end do
c	end do


       deallocate(TMP,stat=ierr)

       end subroutine A1BASIS
c==========================================================================
       subroutine A2BASIS
       use QUADPOINTS, only : A2node, A2weight
       use ROTBASIS, only : DMATA2
       use COMPARAS, only : NA2node, J2max
       implicit none
       integer :: kind, kpts, i, j, k, l, ierr
       real*8 WK(NA2node), endpts(2)
       real*8 alpha, beta
       real*8, allocatable :: TMP(:,:,:,:)

       allocate(A2node(NA2node), A2weight(NA2node), stat=ierr)
       allocate(DMATA2(NA2node,-J2max:J2max,0:J2max,-J2max:J2max),stat=ierr)
       allocate(TMP(NA2node,-J2max:J2max,0:J2max,-J2max:J2max),stat=ierr)

       if (NA2node.eq.1) then
         A2node=0.d0
         DMATA2=1.d0
         return
       end if

       kind=1; kpts=0; alpha=0.d0; beta=0.d0
       call GAUSSQ(kind,NA2node,alpha,beta,kpts,endpts,WK,A2node,A2weight)
       
       do i = 1, NA2node
         A2node(i)=dacos(A2node(i))
       end do

       call ROTMAT1(NA2node,NA2node,J2max,-J2max,J2max,-J2max,J2max,
     &                      A2node,A2weight,TMP)

       DMATA2=TMP

!	print*, "Parameters for A2 basis"
!	print*, "=================================="
!	print*
!	print*, "Points and Weights for A2"
!	!open(10,file="aa2",status="UNKNOWN")
!	do i=1, NA2node
!	  write(*,*) i, A2node(i), A2weight(i)
!	  !write(10,*) A2node(i)
!	end do
!	!close(10)
!	write(*,*)
!	print*, "=================================="

!	do l=-J2max,J2max
!	do k=0,J2max
!	do j=-J2max,J2max
!	do i=1,NA2node
!	  write(206,*) i, j, k, l, TMP(i,j,k,l)
!	end do
!	end do
!	end do
!	end do
!	close(206)

       deallocate(TMP, stat=ierr)

       end subroutine A2BASIS
c==========================================================================
       subroutine B1BASIS
       use QUADPOINTS, only : B1node, B1weight
       use ROTBASIS, only : CTMATB1, RTMATB1
       use COMPARAS, only : NB1node, M1max, parity, Jtot
       implicit none
       integer :: nfoldb1, im1, ib1, ierr, i, j
       real*8  :: AN


       allocate(B1node(NB1node), B1weight(NB1node), stat=ierr)
       allocate(CTMATB1(NB1node,-M1max:M1max), stat=ierr)
       allocate(RTMATB1(NB1node,-M1max:M1max,1:2), stat=ierr)

       nfoldb1=1
       call GAUSCHEB(1,nfoldb1,NB1node,B1node,B1weight)

!	print*, "Parameters for B1 basis"
!	print*, "=================================="
!	print*
!	print*, "Points and Weights for B1"
!	!open(10,file="bb1",status="UNKNOWN")
!	do i=1, NB1node
!	  write(*,*) i, B1node(i), B1weight(i)
!	  !write(10,*) B1node(i)
!	end do
!	!close(10)
!	write(*,*)
!	print*, "=================================="
!	print*

	do im1=-M1max, M1max
	   do ib1=1, NB1node
	      CTMATB1(ib1,im1)=cdexp((0.D0,1.D0)*im1*B1node(ib1))*B1weight(ib1)
	      RTMATB1(ib1,im1,1)=dcos(im1*B1node(ib1))*B1weight(ib1) 
	      RTMATB1(ib1,im1,2)=dsin(im1*B1node(ib1))*B1weight(ib1) 
	   enddo
	enddo

	open(100,file='ctmatb1_norm')
	do i=-M1max, M1max
	do im1=-M1max, M1max
	   AN=0.D0 
	   do ib1=1, NB1node
	      AN=AN+CTMATB1(ib1,im1)*CONJG(CTMATB1(ib1,i))
	   enddo
	   write(100,*) im1,i, AN  
	enddo
	enddo
	close(100)

       end subroutine B1BASIS
c==========================================================================
       subroutine B2BASIS
       use QUADPOINTS, only : B2node, B2weight 
       use ROTBASIS, only : CTMATB2, KB2, RTMATB2
       use COMPARAS, only : NKB2, NB2node, nfoldb2
       implicit none
       integer :: ik2, ib2, kb2Val, ierr, i 
	real*8 :: AN

       allocate(B2node(NB2node), B2weight(NB2node), stat=ierr)
       allocate(CTMATB2(NB2node,1:NKB2), stat=ierr)
       allocate(RTMATB2(NB2node,1:NKB2,1:2), stat=ierr)

       if (NB2node.eq.1) then
         B2node=0.D0
	 RTMATB2(:,:,1)=1.D0
	 RTMATB2(:,:,2)=0.D0
         CTMATB2=(1.D0,0.D0)
         return
       end if
        call GAUSCHEB(nfoldb2,nfoldb2,NB2node,B2node,B2weight)  !###


	print*, "Parameters for B2 basis"
	print*, "=================================="
	print*
	print*, "Points and Weights for B2"
	do i=1, NB2node
	  write(*,*) i, B2node(i), B2weight(i)
	end do
	write(*,*)
	print*, "=================================="
	print*

	do ik2=1,NKB2
           kb2Val=KB2(ik2)
	   do  ib2=1,NB2node
             CTMATB2(ib2,ik2)=cdexp((0.D0,1.D0)*kb2Val*B2node(ib2))*B2weight(ib2)
	     RTMATB2(ib2,ik2,1)=dcos(kb2Val*B2node(ib2))*B2weight(ib2)
	     RTMATB2(ib2,ik2,2)=dsin(kb2Val*B2node(ib2))*B2weight(ib2)
	   enddo
	enddo

	open(100,file='ctmatb2_norm')
	do i=1, NKB2
	do ik2=1, NKB2
	   AN=0.D0 
	   do ib2=1, NB2node
	      AN=AN+CTMATB2(ib2,ik2)*CONJG(CTMATB2(ib2,i))
	   enddo
	   write(100,*) ik2,i, AN  
	enddo
	enddo
	close(100)

	return
       end subroutine B2BASIS
c=============================================================================
	subroutine B3BASIS
	use QUADPOINTS, only : B3node, B3weight 
	use ROTBASIS,    only : CTMATB3 
	use COMPARAS,    only : KBZmax, NB3node, ca0 
	implicit none
	integer  :: nfoldb3, ierr
	integer  :: ik3, ib3

	allocate(B3node(1:NB3node),B3weight(1:NB3node),stat=ierr)
	allocate(CTMATB3(1:NB3node,-KBZmax:KBZmax), stat=ierr)

	if(NB3node.eq.1) then
	  B3node(1)=ca0
	  CTMATB3(1,-KBZmax:KBZmax)=(1.D0,0.D0)
	endif

	nfoldb3=1
	call GAUSCHEB(1,nfoldb3,NB3node,B3node,B3weight)

	do ik3=-KBZmax, KBZmax
	do ib3=1, NB3node
	   CTMATB3(ib3,ik3)=cdexp((0.D0,1.D0)*B3node(ib3))*B3weight(ib3) 
	enddo
	enddo

	return
	end subroutine B3BASIS
c=============================================================================
       subroutine GAUSCHEB(kind,nfoldb1,nodeb1,xmatb1,wmatb1)
       implicit none
       integer :: kind,nfoldb1,nodeb1,inodb1
       real*8  :: xmatb1(nodeb1),wmatb1(nodeb1),delta,pi
c
       pi=dacos(-1.0d0)
       if (kind.eq.1) goto 1000
       if (kind.eq.2) goto 2000
       if (kind.eq.3) goto 3000
       stop 'no such choice for GAUSCHEB'
c
c......kind=1, no symmetry [-pi,pi)

 1000  continue
       delta=2*pi/nodeb1
       do 1200 inodb1=1,nodeb1
          xmatb1(inodb1)=(1.0d0*inodb1-0.5d0)*delta
          wmatb1(inodb1)=dsqrt(delta/(2*pi))
 1200  continue
       return
c
c......kind=2, symmetry plane

 2000  continue
       if (nfoldb1.ne.2) stop 'error found for nfoldb1 in GAUSCHEB 2'
       delta=pi/nodeb1
       do 2200 inodb1=1,nodeb1
          xmatb1(inodb1)=(1.0d0*inodb1-0.5d0)*delta
          wmatb1(inodb1)=dsqrt(delta/(2*pi))
 2200  continue
       return
c
c......kind=3, CnV symmetry

 3000  continue
       if (nfoldb1.le.1) stop 'error found for nfoldb1 in GASUCHEB 3'
       delta=2*pi/nodeb1
       do 3200 inodb1=1,nodeb1
          xmatb1(inodb1)=((1.0d0*inodb1-0.5d0)*delta-pi)/nfoldb1
          wmatb1(inodb1)=dsqrt(delta/(2*pi))
 3200  continue
       return
       end
c===========================================================================
!
!  ***************************************************************************
!  *                                                                         *
!  *   Program    :  rotmat1                                                 *
!  *   Function   :  reduced rotation matrix d[j,kz1,kz2](theta)(normalized) *
!  *                                                                         *
!  ***************************************************************************
c
c..this routine calculate reduced rotation matrix 
c  include dsqrt(weight) of gauss-legendre quadrature
c  the value of reduced rotational matrix at gauss-quadrature nodes 
c  are saved in rotmata1(nodea1,kz2min:kz2max,0:ja2max,kz1min:kz1max) 
c
c  maxnod    :  physic dimension of nodea1
c  nodea1    :  the point to be calculated 
c  ja1max    :  the maximum order of angular momentum 
c  kz1min    :  the minimum value of the first index 
c  kz1max    :  the maximum value of the first index
c  kz2min    :  the minimum value of the second index
c  kz2max    :  the minimum value of the second index
c  xnoda1    :  array saving the value of point
c  wnoda1    :  weight of the point (using Gauss-Legendre quadrature) 
c  rotmata1  :  array saving values of function 

       subroutine ROTMAT1(maxnod,nodea1,ja1max,kz1min,kz1max,
     &                    kz2min,kz2max,xnoda1,wnoda1,rotmata1)
       implicit none 
       integer :: memstat,maxnod,nodea1,inoda1,ja1,kz1,kz2,v,vmin,vmax,
     &            ja1max,kz1min,kz1max,kz2min,kz2max

       real*16 :: factor,sinah,cosah,wsqrt,rotmat
       real*8  :: xnoda1(nodea1),wnoda1(nodea1),
     &            rotmata1(maxnod,kz2min:kz2max,0:ja1max,kz1min:kz1max) 

       real*16,dimension(:),allocatable :: factn,fsinah,fcosah
       allocate(factn(0:4*ja1max),fsinah(0:8*ja1max),fcosah(0:8*ja1max), 
     & stat=memstat)
       if (memstat.ne.0) stop 'failed to allocate memory in ROTMAT0'
c
c..calculte N! used in this routine, log(N!) saved in factn(N)

       if (abs(kz1min).gt.ja1max) stop 'error value of kz1min'
       if (abs(kz1max).gt.ja1max) stop 'error value of kz1max'
       if (abs(kz2min).gt.ja1max) stop 'error value of kz2min'
       if (abs(kz2max).gt.ja1max) stop 'error value of kz2max'
       if (kz1min.gt.kz1max) stop 'error value of kz1min or kz1max'
       if (kz2min.gt.kz2max) stop 'error value of kz2min or kz2max'
       if (nodea1.gt.maxnod) stop 'error value of nodea1 or maxnod'
 
       factn(0)=0.0q0
       do 1200 ja1=1,4*ja1max
          factn(ja1)=factn(ja1-1)+qlog(1.0q0*ja1)
 1200  continue 
c
c..now calculate the reduaced rotation matrix 
c
c       call DZERO1D(maxnod*(kz1max-kz1min+1)*(ja1max+1)*
c     &              (kz2max-kz2min+1),rotmata1) 
       rotmata1(:,:,:,:) = 0.d0

       do 4000 inoda1=1,nodea1
          wsqrt=qsqrt(wnoda1(inoda1)*1.0q0)
          sinah=-qsin(xnoda1(inoda1)/2.0q0)
          cosah=+qcos(xnoda1(inoda1)/2.0q0) 

          fsinah(0)=1.0q0 
          fcosah(0)=1.0q0 
          do 2200 ja1=1,8*ja1max
             fsinah(ja1)=fsinah(ja1-1)*sinah
             fcosah(ja1)=fcosah(ja1-1)*cosah
 2200     continue 

          do 3600 kz1=kz1min,kz1max
          do 3600 kz2=kz2min,kz2max
          do 3600 ja1=0,ja1max
             if (ja1.lt.abs(kz1)) goto 3600 
             if (ja1.lt.abs(kz2)) goto 3600 

             vmin=max(0,kz2-kz1) 
             vmax=min(ja1-kz1,ja1+kz2) 
             
             rotmat=0.0q0 
             do 3000 v=vmin,vmax

c               if (ja1-kz1-v.gt.4*ja1max) goto 8000
c               if (ja1+kz2-v.gt.4*ja1max) goto 8000 
c               if (kz1-kz2+v.gt.4*ja1max) goto 8000 
c               if (v        .gt.4*ja1max) goto 8000 

c               if (ja1-kz1-v.lt.0) goto 8000
c               if (ja1+kz2-v.lt.0) goto 8000 
c               if (kz1-kz2+v.lt.0) goto 8000 
c               if (v        .lt.0) goto 8000 

                factor=factn(ja1-kz1-v)+factn(ja1+kz2-v)
     &                +factn(kz1-kz2+v)+factn(v)
                factor=((-1)**v)/qexp(factor) 

c               if (kz1-kz2+2*v.lt.0) goto 8000
c               if (kz1-kz2+2*v.gt.6*ja1max) goto 8000 
c               if (kz2-kz1-2*v+2*ja1.lt.0) goto 8000
c               if (kz2-kz1-2*v+2*ja1.gt.6*ja1max) goto 8000 

                rotmat=rotmat+factor*fsinah(kz1-kz2+2*v)
     &                              *fcosah(kz2-kz1-2*v+2*ja1) 
 3000        continue 
             factor=factn(ja1+kz1)+factn(ja1-kz1)
     &             +factn(ja1+kz2)+factn(ja1-kz2) 
             factor=qexp(0.5q0*factor)
             factor=factor*qsqrt((2.0q0*ja1+1.0q0)/2) 
             rotmata1(inoda1,kz2,ja1,kz1)=wsqrt*factor*rotmat
 3600     continue
 4000  continue 
c
c release memory 
c
       deallocate(factn,fsinah,fcosah,stat=memstat) 
       if (memstat.ne.0)stop 'failed to release memory in ROTMAT0'
       return
       end 
