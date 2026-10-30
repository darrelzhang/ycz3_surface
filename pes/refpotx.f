!  ***************************************************************************
!  *                                                                         *
!  *   Program    :  potzero                                                 *
!  *   Function   :  potential energy zero point                             *
!  *                                                                         *
!  ***************************************************************************
       subroutine POTZERO
       use COMPARAS, only : VINF, Zinf, R0, RU0, AU0, A10, A20, B10, B20, ca0
       implicit none 
	real*8 vx ,vm   ,R1,RU,AU
	integer i, j ,k
	 call POTDRIV(Zinf,R0,RU0,AU0,A10,A20,B10,B20,ca0, vx)

       VINF=vx

       print*, 'VINF=', VINF

       end
!
!  ***************************************************************************
!  *                                                                         *
!  *   Program    :  prepot2                                                 *
!  *   Function   :  pre-treatment of reference potential of r1              *
!  *                                                                         *
!  ***************************************************************************
       subroutine PREPOT2
       return
       end
!
!  ***************************************************************************
!  *                                                                         *
!  *   Program    :  refpot2                                                 *
!  *   Function   :  reference potential of r1 (X-Y)                         *
!  *                                                                         *
!  ***************************************************************************
       subroutine REFPOT2(noder1,xnodr1,vnodr1) 
       use COMPARAS, only : NVBRU,VINF,Zinf,R0,RUmin,RU0,AU0,A10,A20,B10,B20,ca0
       implicit none 
       integer :: noder1,inodr1
       real*8  :: xnodr1(noder1),vnodr1(noder1),r1,ru

        
       ru=RU0
       if (NVBRU.eq.1) ru=RUmin
       do 2000 inodr1=1,noder1
          r1=xnodr1(inodr1)
          call POTDRIV(Zinf,r1,ru,AU0,A10,A20,B10,B20,ca0,vnodr1(inodr1))
          print*, inodr1, r1, vnodr1(inodr1)
 2000  continue
       return
       end
!
!  ***************************************************************************
!  *                                                                         *
!  *   Program    :  prepot4                                                 *
!  *   Function   :  pre-treatment of reference potential of au              *
!  *                                                                         *
!  ***************************************************************************
       subroutine PREPOT4
       return
       end
!
!  ***************************************************************************
!  *                                                                         *
!  *   Program    :  refpot4                                                 *
!  *   Function   :  reference potential of au, inversion of CZ3             *
!  *                                                                         *
!  ***************************************************************************
       subroutine REFPOT4(nodeau,xnodau,vnodau) 
       use COMPARAS, only : NVBRU,VINF,Zinf,R0,RUmin,RU0,AU0,A10,A20,B10,B20,ca0
       implicit none 
       integer :: nodeau,inodau, i
       real*8  :: xnodau(nodeau),vnodau(nodeau),au,ru 
       real*8 :: au0r,au0p, ru0p,rhf0p,  rr,zz, vx0p

       ru=RU0
       if (NVBRU.eq.1) ru=RUmin
       do 2000 inodau=1,nodeau
          au=xnodau(inodau)
            call POTDRIV(Zinf,R0,ru,au,A10,A20,B10,B20,ca0,vnodau(inodau))
 2000  continue
       return
       end
!
!  ***************************************************************************
!  *                                                                         *
!  *   Program    :  refpot41                                                *
!  *   Function   :  reference potential of au, inversion of CZ3             *
!  *                                                                         *
!  ***************************************************************************
       subroutine REFPOT41(nodeau,xnodau,vnodau)
       implicit none
       integer :: nodeau
       real*8  :: xnodau(nodeau),vnodau(nodeau)
c
       call REFPOT4 (nodeau,xnodau,vnodau)   ! real ref pot
c      call REFPOT4X(nodeau,xnodau,vnodau)   ! fitted ref pot
       return
       end
!
!  ***************************************************************************
!  *                                                                         *
!  *   Program    :  refpot4x                                                *
!  *   Function   :  reference potential of au, inversion of CZ3             *
!  *                                                                         *
!  ***************************************************************************
       subroutine REFPOT4X(nodeau,xnodau,vnodau)
       implicit none
       integer :: nodeau,inodau
       real*8  :: xnodau(nodeau),vnodau(nodeau),au,v2,v4,pi2
c
c      v2=4.04d-2;     v4=0.1586d0
       v2=-5.046d-2;   v4=0.311d0

       pi2=dacos(-1.0d0)/2
c
       do 2000 inodau=1,nodeau
          au=xnodau(inodau)-pi2
          vnodau(inodau)=v2*au**2+v4*au**4
 2000  continue
       return
       end
!
!  ***************************************************************************
!  *                                                                         *
!  *   Program    :  prepot5                                                 *
!  *   Function   :  pre-treatment of reference potential of ycz3            *
!  *                                                                         *
!  ***************************************************************************
       subroutine PREPOT5
       return
       end
!
!  ***************************************************************************
!  *                                                                         *
!  *   Program    :  refpot5                                                 *
!  *   Function   :  reference potential of YCZ3                             *
!  *                                                                         *
!  ***************************************************************************
       subroutine REFPOT5(r1,ru,au,a2,b2,vx)
       use COMPARAS, only : Zinf, ca0
       implicit none
       real*8  :: rs,r1,ru,au,a1,a2,b1,b2,vx
c
       rs=Zinf
       a1=0.0d0
       b1=0.0d0
c
       call POTDRIV(rs,r1,ru,au,a1,a2,b1,b2,ca0, vx) 
       return
       end 
!==========================================================
