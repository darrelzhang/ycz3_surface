!  ***************************************************************************
!  *                                                                         *
!  *   Program    :  potdriv                                                 *
!  *   Function   :  Drive routine of potential energy                       *
!  *                                                                         *
!  ***************************************************************************
       subroutine POTDRIV(rs,r1,ru,au,a1,a2,b1,b2,b3,  vx) 
       use COMPARAS, only : VINF, Zinf, Zinit, Zmax
       implicit none 
       real*8 :: rs,r1,ru,au,a1,a2,b1,b2,b3, vx
       real*8 :: rscut, r1_limit
       integer kcoord

!       rscut=2*Zinit-Zmax
!       if (rs.gt.rscut) then
!         rs=Zinf;a1=0.d0;b1=0.d0
!       end if
!	r1_limit=4.0D0
	if (rs .GT. Zinf) then 
	  rs=Zinf
	endif 

	call POTDRIVA(rs,r1,ru,au,a1,a2,b1,b2,b3,vx) 
	
	vx=vx-VINF

       end 
!
!  ***************************************************************************
!  *                                                                         *
!  *   Program    :  potdriva                                                *
!  *   Function   :  Drive routine of potential energy                       *
!  *                                                                         *
!  ***************************************************************************
       subroutine POTDRIVA(rs,r1,ru,au,a1,a2,b1,b2,b3,vx) 
       use COMPARAS, only : PES_TYPE
       implicit none 
       integer :: j ,i
       real*8  :: rs,r1,ru,au,a1,b1,a2,b2,b3,vx,
     &            cart(3,5),cart1(3,5)
c
c......coordinate transform 
c      CARTES give the coordinates for Y(1),C(2),Z1(3),Z2(4),Z3(5)

	 call CARTES15(rs,r1,ru,au,a1,a2,b1,b2,b3,cart)
	 Call ARRANGE(cart,cart1)  
	 cart1(:,:)=cart1(:,:)*0.529177D0

	 call nninter(cart1,vx)       

	 vx=vx/27.2116D0

       end 
!  ***************************************************************************
!  *                                                                         *
!  *   Program    :  arrange                                                 *
!  *   Function   :  re-arrange cartesian coordinates for CH5 potential      *
!  *                                                                         *
!  ***************************************************************************
       subroutine ARRANGE(cart,cart1) 
c       use comdata, only : kcoord
       implicit none 
       integer :: iatom, j, kcoord
       real*8  :: cart(3,5),cart1(3,5)

       kcoord=15
!=========================================================
       call DCOPY1D(3,cart(1,1),cart1(1,1))
       call DCOPY1D(3,cart(1,2),cart1(1,5))
       call DCOPY1D(3,cart(1,3),cart1(1,2))
       call DCOPY1D(3,cart(1,4),cart1(1,3))
       call DCOPY1D(3,cart(1,5),cart1(1,4))
!==========================================================

c......make all components of the first atom is postive

       do 6000 j=1,3
          if (cart1(j,1).lt.0.0d0) then
             do 5800 iatom=1,5
                cart1(j,iatom)=-cart1(j,iatom)
 5800        continue
          endif
 6000  continue
       return
       end 
!
!  ***************************************************************************
!  *                                                                         *
!  *   Program    :  distest                                                 *
!  *   Function   :  print distances                                         *
!  *                                                                         *
!  ***************************************************************************
       subroutine DISTEST(kchanl,natom,cart)
       implicit none 
       integer :: kchanl,natom,i,j,k
       real*8  :: cart(3,natom),dist(1000),rij
c
       do 4000 i=1,natom 
          do 3800 j=1,i
             rij=0 
             do 3000 k=1,3
                rij=rij+(cart(1,i)-cart(1,j))**2
 3000        continue 
             dist(j)=dsqrt(rij) 
 3800     continue 
          write(kchanl,9400)(dist(j),j=1,i) 
 4000  continue 
 9400  format(6f12.5) 
       return
       end 
!
!  ***************************************************************************
!  *                                                                         *
!  *   Program    :  cartest                                                 *
!  *   Function   :  print cartesian                                         *
!  *                                                                         *
!  ***************************************************************************
       subroutine CARTEST(kchanl,natom,cart) 
       implicit none 
       integer :: kchanl,natom,i
       real*8  :: cart(3,natom)
c
       do 4000 i=1,natom
          write(kchanl,9400)cart(1,i),cart(2,i),cart(3,i) 
 4000  continue 
 9400  format(6f12.5) 
       return
       end 
!  ***************************************************************************
!  *                                                                         *
!  *   Program    :  cartes15                                                *
!  *   Function   :  X-YCZ3 system, jacobi to cartesian                      *
!  *                                                                         *
!  ***************************************************************************
       subroutine CARTES15(rs,r1,ru,au,a1,a2,b1,b2,b3, cartes)
       use COMPARAS, only : AM, X00, Y00
       implicit none
       integer :: icoord,iatom, j
       real*8  :: rs,r1,ru,au,a1,b1,a2,b2,b3,cartes(3,5),cent, cent_ch4(3)
c
c......the eular angle a,b,c is defined as
C      the rotation is always from space-fixed to body-fixed frame
C
C      a => rotate (X,Z,z) in a plane, if in a plane, a=0
C      b => rotate Z-axis to z-axis
C      c => rotate X-axis to x-axis
C
C
C      at first, the molecule is placed on the space-fixed frame
c      the rotate the atoms to body-fixed frame, for the case 
c
c      1) CZ3 symmetry axis is body-fixed z-axis and also define 
C         the X-axis of space-fixed. so a=0, 
C      2) Y-CZ3 axis is body-finxed z-axis and also define 
C         the X-axis of space-fixed
C
c......calculate cartesian of CZ3 group, in their body-fixed frame
 
       call CARTESCZ3(AM(2),AM(3),ru,au,cartes(1,2))
c
c......rotate CZ3

       call ROTEULAR(4,cartes(1,2),0.0d0,a2,b2)
c
c......the coordinates of Y atom 

       cartes(1,1)=0.0d0
       cartes(2,1)=0.0d0
       cartes(3,1)=r1
c
c......move the COM of YCZ3 to origin 

       cent=cartes(3,1)*AM(1)/(AM(1)+AM(2)+3*AM(3))
       do 2400 iatom=1,5
          cartes(3,iatom)=cartes(3,iatom)-cent
 2400  continue
c
c......rotate YCZ3

!	 gama=0.D0
!	 gama=30.D0/180.D0*DACOS(-1.D0)
!	 gama=150.D0/180.D0*DACOS(-1.D0)
       call ROTEULAR(5,cartes(1,1),b3,a1,b1)
c
! top 
! atop(1.249,2.160); fcc(1.249,3.594); hcc(2.494,2.875)
! top
!	cent_ch4(1)=1.24899995326996D0/0.529177D0
!	cent_ch4(2)=2.16000008583069D0/0.529177D0
	cent_ch4(1)=X00
	cent_ch4(2)=Y00
! fcc
!	cent_ch4(1)=1.249D0/0.529177D0
!	cent_ch4(2)=3.594D0/0.529177D0
! hcc
!	cent_ch4(1)=2.494D0/0.529177D0
!	cent_ch4(2)=2.875D0/0.529177D0

	cent_ch4(3)=rs+6.097899/0.529177D0
c
       do 3000 icoord=1,3
          do 2800 iatom =1,5
             cartes(icoord,iatom)=cartes(icoord,iatom)+cent_ch4(icoord)
 2800     continue 
 3000  continue 

 9100  format('coordinates'/(3f10.5))
       return
       end
!
!  ***************************************************************************
!  *                                                                         *
!  *   Program    :  roteular                                                *
!  *   Function   :  rotate a vector by eular angle (a,b,c)                  *
!  *                                                                         *
!  ***************************************************************************
       subroutine ROTEULAR(natom,cartes,a,b,c) 
       implicit none
       integer :: natom,iatom,ixyz,jxyz
       real*8  :: a,b,c,tmatrot(3,3),cartes(3,natom),
     &            sina,cosa,sinb,cosb,sinc,cosc,temp(3)
c
       sina=dsin(a); cosa=dcos(a)
       sinb=dsin(b); cosb=dcos(b)
       sinc=dsin(c); cosc=dcos(c)
c
       tmatrot(1,1)= cosa*cosb*cosc-sina*sinc
       tmatrot(2,1)= sina*cosb*cosc+cosa*sinc
       tmatrot(3,1)=-sinb*cosc

       tmatrot(1,2)=-cosa*cosb*sinc-sina*cosc
       tmatrot(2,2)=-sina*cosb*sinc+cosa*cosc
       tmatrot(3,2)= sinb*sinc

       tmatrot(1,3)= cosa*sinb
       tmatrot(2,3)= sina*sinb
       tmatrot(3,3)= cosb
c
       do 3000 iatom=1,natom
          do 2000 ixyz=1,3
             temp(ixyz)=0.0d0
             do 1800 jxyz=1,3
                temp(ixyz)=temp(ixyz)
     &                    +tmatrot(ixyz,jxyz)*cartes(jxyz,iatom)
 1800        continue
 2000     continue
          do 2200 ixyz=1,3
             cartes(ixyz,iatom)=temp(ixyz)
 2200     continue
 3000  continue
       return
       end 
!
!  ***************************************************************************
!  *                                                                         *
!  *   Program    :  cartescz3                                               *
!  *   Function   :  CZ3 system, jacobi to cartesian                         *
!  *                                                                         *
!  ***************************************************************************
       subroutine cartescz3(amc,amz,ru,au,cartes)
       implicit none 
       integer :: iatom,ixyz
       real*8  :: amc,amz,ru,au,cartes(3,4),cent,xu,zu
       real*8  :: sin120,cos120
       data       sin120/+0.866025403784439d0/,
     &            cos120/-0.500000000000000d0/
c
c......the coordinates for CZ3 molecule, C at origin at first 

       xu=ru*dsin(au)
       zu=ru*dcos(au)

       cartes(1,1)=0.0d0
       cartes(2,1)=0.0d0
       cartes(3,1)=0.0d0

       cartes(1,2)=xu
       cartes(2,2)=0.0d0
       cartes(3,2)=zu

       cartes(1,3)=xu*cos120
       cartes(2,3)=xu*sin120
       cartes(3,3)=zu

       cartes(1,4)=+cartes(1,3)
       cartes(2,4)=-cartes(2,3)
       cartes(3,4)=zu
c
c......move the COM of CZ3 to origin, only z-axis coordinates are changed

       cent=3*zu*amz/(amc+3*amz)
       do 1400 iatom=1,4
          cartes(3,iatom)=cartes(3,iatom)-cent
 1400  continue
       return
       end 
c ==========================================================================
c      copy 1-dimen real array
c ==========================================================================
       subroutine dcopy1d(ndim,xdat,ydat)
       implicit none
       integer :: ndim,idim
       real*8  :: xdat(ndim),ydat(ndim)
       do 1200 idim=1,ndim
          ydat(idim)=xdat(idim)
 1200  continue
       return
       end
c========================================================================
