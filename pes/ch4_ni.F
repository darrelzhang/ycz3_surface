	module nnparam
	implicit none
	save 

	integer   :: ndim=25

	end module
! --------------------------------------------------------------------

        subroutine nninter(r0,vx)
          use nnparam
        implicit none
        real*8,intent(in) :: r0(15)
        real*8,intent(out) :: vx
        integer dx
        real*8 :: r(ndim),dv(ndim)
        integer,save :: init_nn=0
        real*8 v1,v2,f1,chmax,zcm
  
           call nnenergy(r0,r,f1,zcm)          
            if(zcm.ge.8.45) then 
           call pot_ch4_part1(r,v1)
              vx=v1
  
          elseif(zcm.gt.8.4) then  
           call pot_ch4_part1(r,v1)
           call pot_ch4_part2(r,v2)      
            vx=v1*f1+v2*(1-f1) 
  
          else
           call pot_ch4_part2(r,v2)
            vx=v2
         endif  
          return
        end subroutine
! --------------------------------------------------------------------
!       subroutine nnenergy(r0,r,f1,chmax)
        subroutine nnenergy(r0,r,f1,zcm)
        use nnparam
         implicit none
        real*8 :: r0(15),r(ndim),r4,d(15)
        integer :: i
        real*8 zcm
        real*8 f1,pi
                
!       xcm= (r0(1)+r0(4)+r0(7)+r0(10)+r0(13)*12.)/16.
!       ycm= (r0(2)+r0(5)+r0(8)+r0(11)+r0(14)*12.)/16.        
!       call PESPBC(r0,xcm,ycm)
        r(1)=sqrt((r0(1)-r0(13))**2+(r0(2)-r0(14))**2+(r0(3)-r0(15))**2)
        r(2)=sqrt((r0(4)-r0(13))**2+(r0(5)-r0(14))**2+(r0(6)-r0(15))**2)
        r(3)=sqrt((r0(7)-r0(13))**2+(r0(8)-r0(14))**2+(r0(9)-r0(15))**2)
          r4=(r0(10)-r0(13))**2+(r0(11)-r0(14))**2+(r0(12)-r0(15))**2
        r(4)=sqrt(r4)
        r(5)=sqrt((r0(1)-r0(4))**2+(r0(2)-r0(5))**2+(r0(3)-r0(6))**2)
        r(6)=sqrt((r0(1)-r0(7))**2+(r0(2)-r0(8))**2+(r0(3)-r0(9))**2) 
        r(7)=sqrt((r0(1)-r0(10))**2+(r0(2)-r0(11))**2+(r0(3)-r0(12))**2)
        r(8)=sqrt((r0(4)-r0(7))**2+(r0(5)-r0(8))**2+(r0(6)-r0(9))**2)
        r(9)=sqrt((r0(4)-r0(10))**2+(r0(5)-r0(11))**2+(r0(6)-r0(12))**2)
        r(10)=sqrt((r0(7)-r0(10))**2+(r0(8)-r0(11))**2+(r0(9)-r0(12))**2)
!*******  Center of mass only for CH4 *******************
!       zcm= (r0(3)+r0(6)+r0(9)+r0(12)+r0(15)*12.)/16.
        zcm=r0(15)
        call dis2surf(r0,d)
        r(11:ndim)=d(1:15)  
!************to choose which fitting functional *********
         pi=4.0*atan(1.0) 
        if(zcm.ge.8.45) then
            f1=1.0
        elseif(zcm.le.8.4) then
            f1=0.0
        else
            f1=0.5*(1-cos(pi*(zcm-8.4)/0.05))
        endif


!************ end **************************************
        return
        end subroutine
!_______________________________________________________
           subroutine dis2surf(r0,d)
          implicit none
         real*8 r0(15),d(15),p(5,3)
         real*8 x1,x2,x3,y1,y2,y3
         real*8 z1,z2,z3
         integer i
         x1=1.247; y1=2.160  ! atop site
         x2=2.494; y2=2.875   ! hcc  site
         x3=1.247; y3=3.594   ! fcc  site
         z1=6.098;z2=4.079;z3=2.036
          do i=1,5
            p(i,1)=r0(i*3-2)
            p(i,2)=r0(i*3-1)
            p(i,3)=r0(i*3)
            d(i*3-2)=sqrt((p(i,1)-x1)**2+(p(i,2)-y1)**2+(p(i,3)-z1)**2)
            d(i*3-1)=sqrt((p(i,1)-x2)**2+(p(i,2)-y2)**2+(p(i,3)-z2)**2)
            d(i*3)=sqrt((p(i,1)-x3)**2+(p(i,2)-y3)**2+(p(i,3)-z3)**2)
          enddo
         return
        end
!!! note: r in ang. and vfit in eV.
!!!   pot_ch4(r0,v)
#DEFINE PESDIR "./"
! --------------------------------------------------------------------
        subroutine pot_ch4_part1(r0,v) ! input bond distances of CH4 in ang.
        implicit none
        real*8,intent(in)  :: r0(25)
        real*8,intent(out) :: v
        real*8 :: r(25),r_all(25,24),v_all(24),w_all(24)
        integer :: ir(25)
        r=r0
        ir=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,
     &       24,25]
        call sortdist_ch4(r,ir)
        call copy24_ch4(r,r_all)
        call weight24_ch4(r_all,w_all)       
        call nnsymfit_ch4_part1(25,24,r_all,v_all)
        v_all=v_all*w_all
        v=sum(v_all)/sum(w_all)    !/27.21138505d0
        return
        end subroutine pot_ch4_part1
! --------------------------------------------------------------------
        subroutine pot_ch4_part2(r0,v) ! input bond distances of CH4 in ang.
        implicit none
        real*8,intent(in)  :: r0(25)
        real*8,intent(out) :: v
        real*8 :: r(25),r_all(25,24),v_all(24),w_all(24)
        integer :: ir(25)
        r=r0
        ir=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,
     &       24,25]
        call sortdist_ch4(r,ir)
        call copy24_ch4(r,r_all)
        call weight24_ch4(r_all,w_all)
        call nnsymfit_ch4_part2(25,24,r_all,v_all)
        v_all=v_all*w_all
        v=sum(v_all)/sum(w_all)     !    /27.21138505d0
        return
        end subroutine pot_ch4_part2
  
! --------------------------------------------------------------------
        subroutine sortdist_ch4(r,ir)
        implicit none
        real*8,intent(inout) :: r(25)
        integer,intent(inout) :: ir(25)
        integer,parameter :: id(5,8)=[-1, 1, 2, 3, 4, 1,-1, 5, 6, 7, 
     &                   2, 5,-1, 8, 9, 3, 6, 8,-1,10, 4, 7, 9,10,-1,
     &                   23, 11, 14, 17, 20, 24, 12, 15, 18, 21, 25, 
     &                   13, 16, 19, 22]
        integer i,j,k,itmp
        real*8 tmp
        do i=2,4
        do j=i+1,5
        if( r(i-1).gt.r(j-1) ) then   !  increasing
          ! exchange Hi with Hj
          do k=1,8
            if(id(i,k).eq.-1 .or. id(j,k).eq.-1) cycle
            tmp=r(id(i,k))
            r(id(i,k))=r(id(j,k))
            r(id(j,k))=tmp
            itmp=ir(id(i,k))
            ir(id(i,k))=ir(id(j,k))
            ir(id(j,k))=itmp
          enddo
        endif
        enddo
        enddo
        return
        end subroutine sortdist_ch4
! --------------------------------------------------------------------
        subroutine copy24_ch4(r,r_all)
        implicit none
        real*8,intent(in) :: r(25)
        real*8,intent(out) :: r_all(25,24)
        r_all(1:25,  1)=r([ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 
     &                     14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25])
        r_all(1:25,  2)=r([ 1, 2, 4, 3, 5, 7, 6, 9, 8, 10, 11, 12, 13, 
     &                     14, 15, 16, 20, 21, 22, 17, 18, 19, 23, 24, 25])
        r_all(1:25,  3)=r([ 1, 3, 2, 4, 6, 5, 7, 8, 10, 9, 11, 12, 13, 
     &                     17, 18, 19, 14, 15, 16, 20, 21, 22, 23, 24, 25])
        r_all(1:25,  4)=r([ 1, 3, 4, 2, 6, 7, 5, 10, 8, 9, 11, 12, 13, 
     &                     17, 18, 19, 20, 21, 22, 14, 15, 16, 23, 24, 25])
        r_all(1:25,  5)=r([ 1, 4, 2, 3, 7, 5, 6, 9, 10, 8, 11, 12, 13, 
     &                     20, 21, 22, 14, 15, 16, 17, 18, 19, 23, 24, 25])
        r_all(1:25,  6)=r([ 1, 4, 3, 2, 7, 6, 5, 10, 9, 8, 11, 12, 13, 
     &                     20, 21, 22, 17, 18 ,19, 14, 15, 16, 23, 24, 25])
        r_all(1:25,  7)=r([ 2, 1, 3, 4, 5, 8, 9, 6, 7, 10, 14, 15, 16, 
     &                     11, 12, 13, 17, 18, 19, 20, 21, 22, 23, 24, 25])
        r_all(1:25,  8)=r([ 2, 1, 4, 3, 5, 9, 8, 7, 6, 10, 14, 15, 16, 
     &                     11, 12, 13, 20, 21, 22, 17, 18, 19, 23, 24, 25])
        r_all(1:25,  9)=r([ 2, 3, 1, 4, 8, 5, 9, 6, 10, 7, 14, 15, 16, 
     &                     17, 18, 19, 11, 12, 13, 20, 21, 22, 23, 24, 25])
        r_all(1:25, 10)=r([ 2, 3, 4, 1, 8, 9, 5, 10, 6, 7, 14, 15, 16, 
     &                     17, 18, 19, 20, 21, 22, 11, 12, 13, 23, 24, 25])
        r_all(1:25, 11)=r([ 2, 4, 1, 3, 9, 5, 8, 7, 10, 6, 14, 15, 16, 
     &                     20, 21, 22, 11, 12, 13, 17, 18, 19, 23, 24, 25])
        r_all(1:25, 12)=r([ 2, 4, 3, 1, 9, 8, 5, 10, 7, 6, 14, 15, 16, 
     &                     20, 21, 22, 17, 18, 19, 11, 12, 13, 23, 24, 25])
        r_all(1:25, 13)=r([ 3, 1, 2, 4, 6, 8, 10, 5, 7, 9, 17, 18, 19, 
     &                     11, 12, 13, 14, 15, 16, 20, 21, 22, 23, 24, 25])
        r_all(1:25, 14)=r([ 3, 1, 4, 2, 6, 10, 8, 7, 5, 9, 17, 18, 19, 
     &                     11, 12, 13, 20, 21, 22, 14, 15, 16, 23, 24, 25])
        r_all(1:25, 15)=r([ 3, 2, 1, 4, 8, 6, 10, 5, 9, 7, 17, 18, 19, 
     &                     14, 15, 16, 11, 12, 13, 20, 21, 22, 23, 24, 25])
        r_all(1:25, 16)=r([ 3, 2, 4, 1, 8, 10, 6, 9, 5, 7, 17, 18, 19, 
     &                     14, 15, 16, 20, 21, 22, 11, 12, 13, 23, 24, 25])
        r_all(1:25, 17)=r([ 3, 4, 1, 2, 10, 6, 8, 7, 9, 5, 17, 18, 19, 
     &                     20, 21, 22, 11, 12, 13, 14, 15, 16, 23, 24, 25])
        r_all(1:25, 18)=r([ 3, 4, 2, 1, 10, 8, 6, 9, 7, 5, 17, 18, 19, 
     &                     20, 21, 22, 14, 15, 16, 11, 12, 13, 23, 24, 25])
        r_all(1:25, 19)=r([ 4, 1, 2, 3, 7, 9, 10, 5, 6, 8, 20, 21, 22, 
     &                     11, 12, 13, 14, 15, 16, 17, 18, 19, 23, 24, 25])
        r_all(1:25, 20)=r([ 4, 1, 3, 2, 7, 10, 9, 6, 5, 8, 20, 21, 22, 
     &                     11, 12, 13, 17, 18, 19, 14, 15, 16, 23, 24, 25])
        r_all(1:25, 21)=r([ 4, 2, 1, 3, 9, 7, 10, 5, 8, 6, 20, 21, 22, 
     &                     14, 15, 16, 11, 12, 13, 17, 18, 19, 23, 24, 25])
        r_all(1:25, 22)=r([ 4, 2, 3, 1, 9, 10, 7, 8, 5, 6, 20, 21, 22, 
     &                     14, 15, 16, 17, 18, 19, 11, 12, 13, 23, 24, 25])
        r_all(1:25, 23)=r([ 4, 3, 1, 2, 10, 7, 9, 6, 8, 5, 20, 21, 22, 
     &                     17, 18, 19, 11, 12, 13, 14, 15, 16, 23, 24, 25])
        r_all(1:25, 24)=r([ 4, 3, 2, 1, 10, 9, 7, 8, 6, 5, 20, 21, 22, 
     &                     17, 18, 19, 14, 15, 16, 11, 12, 13, 23, 24, 25])
        return
        end subroutine copy24_ch4
! --------------------------------------------------------------------
        subroutine weight24_ch4(r_all,w_all)
        implicit none
        real*8,intent(in) :: r_all(25,24)
        real*8,intent(out) :: w_all(24)
        real*8 dw
        integer i
        do i=1,24
          dw=dist(r_all(1:4,1),r_all(1:4,i),4)
          w_all(i)=dexp(-20.d0*dw**2)
        enddo
        return
        contains
  
          function dist(r1,r2,ndim)
          implicit none
          real*8,intent(in) :: r1(ndim),r2(ndim)
          integer,intent(in) :: ndim
          integer i
          real*8 dist
          dist=0.d0
          do i=1,ndim
            dist=dist+(r1(i)-r2(i))**2
          enddo
          dist=dsqrt(dist)
          return
          end function dist
  
        end subroutine weight24_ch4
! --------------------------------------------------------------------
        subroutine nnsymfit_ch4_part1(ndim,ntot,r0,vx)
        implicit none
        integer,intent(in) :: ndim,ntot
        real*8,intent(in) :: r0(ndim,ntot)
        real*8,intent(out) :: vx(ntot)
  
        integer :: i
        real*8 :: dv(ndim,ntot),dvx(ndim),tmp
  
        integer,save :: init=0
        integer,save :: fid
        character*10,save :: fp="file pointer"
        integer,parameter :: s0=25,s1=60,s2=60
        character*99,save :: wfile='nnpara1.in' ! dist in ang., v in eV
        real*8,save :: w1(s0,s1),b1(s1),w2(s1,s2),b2(s2),w3(s2),b3
        real*8,save :: rg(2,s0),vg(2)
  
        if(ndim.ne.s0) stop "ndim .ne. s0"
  
        if (init.eq.0) then
          fid=123456
          open(fid,file=PESDIR//trim(wfile),status='old',action='read')
          read(fid,*) w1,b1,w2,b2,w3,b3,rg,vg
          close(fid)
          init=1
        endif
  
        if(ntot.ge.24) then
          call nsimx(ntot,r0,vx,0,dv,s0,s1,s2,rg,w1,b1,w2,b2,w3,b3,vg)
        else
          do i=1,ntot
          call nsim(r0(1,i),vx(i),0,dvx,s0,s1,s2,rg,w1,b1,w2,b2,w3,b3,vg)
          enddo
        endif
  
        return
        end subroutine

! --------------------------------------------------------------------
        subroutine nnsymfit_ch4_part2(ndim,ntot,r0,vx)
        implicit none
        integer,intent(in) :: ndim,ntot
        real*8,intent(in) :: r0(ndim,ntot)
        real*8,intent(out) :: vx(ntot)
  
        integer :: i
        real*8 :: dv(ndim,ntot),dvx(ndim),tmp
  
        integer,save :: init=0
        integer,save :: fid
        character*10,save :: fp="file pointer"
        integer,parameter :: s0=25,s1=60,s2=60
        character*99,save :: wfile='nnpara2.in' ! dist in ang., v in eV
        real*8,save :: w1(s0,s1),b1(s1),w2(s1,s2),b2(s2),w3(s2),b3
        real*8,save :: rg(2,s0),vg(2)
  
        if(ndim.ne.s0) stop "ndim .ne. s0"
  
        if (init.eq.0) then
          fid=123456
          open(fid,file=PESDIR//trim(wfile),status='old',action='read')
          read(fid,*) w1,b1,w2,b2,w3,b3,rg,vg
          close(fid)
          init=1
        endif
  
        if(ntot.ge.24) then
          call nsimx(ntot,r0,vx,0,dv,s0,s1,s2,rg,w1,b1,w2,b2,w3,b3,vg)
        else
          do i=1,ntot
          call nsim(r0(1,i),vx(i),0,dvx,s0,s1,s2,rg,w1,b1,w2,b2,w3,b3,vg)
          enddo
        endif
  
        return
        end subroutine

! ====================================================================
! ---- simulate a neural network with two hidden layers
! ----      n0-n1-n2-1
! ---- blas routines used in this subroutine: dgemv, ddot
        subroutine nsim(r0,v,idv,dv,n0,n1,n2,rg,w1,b1,w2,b2,w3,b3,vg)
        implicit none
        integer,intent(in) :: n0,n1,n2,idv
        real*8,intent(in)  :: r0(n0),rg(2,n0),vg(2)
        real*8,intent(in)  :: w1(n0,n1),b1(n1)
        real*8,intent(in)  :: w2(n1,n2),b2(n2)
        real*8,intent(in)  :: w3(n2),b3
        real*8,intent(out) :: v,dv(n0)
        integer :: i,j,k
        real*8  :: r(n0),rgg(n0),vgg,ax(n1),bx(n2)
        real*8  :: dvtm,rtmp,rt1(n1),rt2(n2)
        real*8,external :: ddot
        v=0.d0
        r=r0
        vgg=vg(2)-vg(1)
        ! mapminmax [-1,1]
        do i=1,n0
          rgg(i)=rg(2,i)-rg(1,i)
          r(i)=2.d0*(r(i)-rg(1,i))/rgg(i)-1.d0
        end do
        ! 1st layer
        rt1=b1
        call dgemv('t',n0,n1,1.d0,w1,n0,r,1,1.d0,rt1,1)
        ax=dtanh(rt1)
        ! 2nd layer
        rt2=b2
        call dgemv('t',n1,n2,1.d0,w2,n1,ax,1,1.d0,rt2,1)
        bx=dtanh(rt2)
        ! output layer
        v=b3+ddot(n2,w3,1,bx,1)
        !reverse map
        v=vgg*(v+1.d0)/2+vg(1)
        if(idv.ne.1) return
        ! calculate first derivatives, dv(i)=dv/dr(i)
        dv=0.d0
        do i=1,n0
          do k=1,n2
            dvtm=0.d0
            do j=1,n1
              dvtm=dvtm+w2(j,k)*w1(i,j)*(1-ax(j)**2)
            enddo
            dv(i)=dv(i)+w3(k)*dvtm*(1-bx(k)**2)
          enddo
          dv(i)=dv(i)*vgg/rgg(i)
        enddo
        return
        end subroutine nsim
! --------------------------------------------------------------------
! ---- simulate a neural network with two hidden layers
! ----    n0-n1(tansig)-n2(tansig)-1(pureline)
! ---- blas routines used in this subroutine: dgemm, dgemv
! ---- the syntax keeps the same with nsim, except input nt at first
! ---- Chen Jun, 2013/11/26
!
        subroutine nsimx(nt,r,v,idv,dv,n0,n1,n2,rg,w1,b1,w2,b2,w3,b3,vg)
!
! ---- description of parameters
! ---- nt  : total number of points in r,v and dv
! ---- r   : input of nn, r(n0,nt), r(1:n0,i) is the i-th input
! ---- v   : output of nn, v(nt), v(i) is the i-th input
! ---- idv : if idv==0, calc. v only, else calc. v and derivatives dv
! ---- dv  : derivatives of dv/dr, dv(n0,nt), dv(j,i)=d v(i) / d r(j,i)
! ---- n0  : neurons of the input layer, dimension of r
! ---- n1  : neurons of the first hidden layer
! ---- n2  : neurons of the second hidden layer
! ---- rg  : input ranges of the training set, rg(2,n0)
! ---- w1  : weights connect input layer with 1-st hidden layer
! ---- b1  : biases of 1-st hidden layer
! ---- w2  : weights connect 1-st hidden layer with 2-nd hidden layer
! ---- b2  : biases of 2-nd hidden layer
! ---- w3  : weights connect 2-nd hidden layer with output layer
! ---- b3  : biase if output layer
! ---- vg  : output range of the training set, vg(2)
!
        implicit none
        integer,intent(in) :: nt,idv,n0,n1,n2
        real*8,intent(in) :: r(n0,nt)
        real*8,intent(out) :: v(nt),dv(n0,nt)
!
        real*8,intent(in) :: w1(n0,n1),b1(n1),w2(n1,n2),b2(n2),w3(n2),b3
        real*8,intent(in) :: rg(2,n0),vg(2)
!
        real*8 :: x0(n0,nt),x1(nt,n1),x2(nt,n2),tmp,rgg(n0),vgg
        integer :: i,j,k,n
!
        v=0.d0
        dv=0.d0
        do i=1,n0
          rgg(i)=rg(2,i)-rg(1,i)
        enddo
        vgg=vg(2)-vg(1)
!
        x0=r
!
        do i=1,n0
          x0(i,1:nt)=(x0(i,1:nt)-rg(1,i))/rgg(i)*2.d0-1.d0
        enddo
!
        do i=1,n1
         x1(1:nt,i)=b1(i)
        enddo
!
! ---- x1=trans(x0(n0,nt))*w1(n0,n1)+x1(nt,n1)
        call dgemm('t','n',nt,n1,n0,1.d0,x0,n0,w1,n0,1.d0,x1,nt)
!
        x1=dtanh(x1)
!
        do i=1,n2
         x2(1:nt,i)=b2(i)
        enddo
!
! ---- x2=x1(nt,n1)*w2(n1,n2)+x2(nt,n2)
        call dgemm('n','n',nt,n2,n1,1.d0,x1,nt,w2,n1,1.d0,x2,nt)
!
        x2=dtanh(x2)
!
        v(1:nt)=b3
!
! ---- v=x2(nt,n2)*w3(n2)+v(nt)
        call dgemv('n',nt,n2,1.d0,x2,nt,w3,1,1.d0,v,1)
!
        v=(v+1.d0)/2.d0*vgg+vg(1)
!
        if (idv.eq.0) return
! ---- calculate first derivatives
        do n=1,nt
          do i=1,n0
            do k=1,n2
              tmp=0.d0
              do j=1,n1
                tmp=tmp+w2(j,k)*w1(i,j)*(1-x1(n,j)**2)
              enddo
              dv(i,n)=dv(i,n)+w3(k)*tmp*(1-x2(n,k)**2)
            enddo
            dv(i,n)=dv(i,n)*vgg/(rgg(i))
          enddo
        enddo
        return
        end subroutine nsimx
! ====================================================================
