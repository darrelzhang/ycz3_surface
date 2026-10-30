c===============================================================
! seven-dimensional potential surface
        subroutine POTENTIALS
        use QUADPOINTS, only : ZQ,RQ,RQI
        use ENGARRAYS, only : VMASY, VMINT, VMASY_SAVE, VMINT_SAVE, VCUT4D_ASY,
     .    VCUT4D_INT
        use COMPARAS, only : NZ,NZINT,NVB,NRB,N2D,NA1node,NA2node,
     .    NB1node,NB2node,IDPOT,NPOTS,IPOT,PES_TYPE, NCARDS_ASY,NCARDS_INT,VCUT
        use SMALL_ARRAYS, only : IDVSYM, LINK1_ASY, LINK2_ASY,
     .    LINK1_INT, LINK2_INT, VrefCH, Vref2D

        implicit none

	include 'mpif.h'

        integer ierr, NASYNODES, NINTNODES, i, j, k, l, m, n,  ms, ns
        integer NANODES, NRNODES_ASY, NRNODES_INT, NZ1, NZ2
        integer NL, N1, N2, NPL, irnode, icards, idatas

        !real*4 vtmp
        real*4, allocatable :: vsave(:)

        character*2 ITOC
        character*3 ITOC3
        character*4 ITOC4

        NANODES=NA1node*NA2node*NB1node*NB2node
        NRNODES_ASY=(NZ-NZINT)*NVB*N2D
        NRNODES_INT=NZINT*NRB*N2D

        NASYNODES=(NZ-NZINT)*NVB*N2D*NANODES
        NINTNODES=NZINT*NRB*N2D*NANODES

        allocate(IDVSYM(NB1node*NB2node), stat=ierr)
        call SYMMPOT(NB1node,NB2node,IDVSYM)

	if (IDPOT.eq.1) then

	  call POT8D_MPI(IPOT-1, NPOTS)

        else if (IDPOT.eq.2) then  ! do the partial potential calculation

!	  goto 9997

	  print*,"Start Calc ASY Potential"

          print*, 'NZ=', NZ, 'NZINT=', NZINT, 'NPOTS=', NPOTS, 'IPOT=', IPOT
          NZ1=(NZ-NZINT)*(IPOT-1)/NPOTS+1+NZINT
          NZ2=(NZ-NZINT)*IPOT/NPOTS+NZINT
          print*, 'NZ1=', NZ1, "Z1=", ZQ(NZ1), 'NZ2=', NZ2, "Z2=", ZQ(NZ2)

c          allocate(VMASY(NANODES,(NZ2-NZ1+1)*NVB*N2D),stat=ierr)
          allocate(VMASY(NANODES,1),stat=ierr)

!          open(1,file='POT'//itoc(IPOT)//'.ASY',status='UNKNOWN',
          open(1,file='POT'//itoc4(IPOT)//'.ASY',status='UNKNOWN',
     .      	 form='binary')

	  print*,'POT'//itoc4(IPOT)//'.ASY'

          !if (PES_TYPE.eq.1) call POT8D(1, NZ1,NZ2,NVB,ZQ,RQ,VMASY)
          !if (PES_TYPE.eq.2) call POT8D_CLUSTER(NZ1,NZ2,NVB,ZQ,RQ,VMASY)
          call POT8D(1, NZ1,NZ2,NVB,ZQ,RQ,VMASY)

c          do i=1, (NZ2-NZ1+1)*NVB*N2D
cc            do j=1, NANODES
cc              write(1) VMASY(j,i)
cc            end do
c            write(1) VMASY(1:NANODES,i)
c          end do

          close(1)

          deallocate(VMASY, stat=ierr)

9997	  continue	  
	  
	  print*,"Start Calc INT Potential"        
          NZ1=NZINT*(IPOT-1)/NPOTS+1
          NZ2=NZINT*IPOT/NPOTS
          print*, 'NZ1=', NZ1, "Z1=", ZQ(NZ1), 'NZ2=', NZ2, "Z2=", ZQ(NZ2)

c          allocate(VMINT(NANODES,(NZ2-NZ1+1)*NRB*N2D),stat=ierr)
          allocate(VMINT(NANODES,1),stat=ierr)

!          open(2, file='POT'//itoc(IPOT)//'.INT',status='UNKNOWN',
          open(2, file='POT'//itoc4(IPOT)//'.INT',status='UNKNOWN',
     .      	  form='binary')

	  print*,'POT'//itoc4(IPOT)//'.INT'

          !if(PES_TYPE.eq.1) call POT8D(2, NZ1,NZ2,NRB,ZQ,RQI,VMINT)
          !if(PES_TYPE.eq.2) call POT8D_CLUSTER(NZ1,NZ2,NRB,ZQ,RQI,VMINT)
          call POT8D(2, NZ1,NZ2,NRB,ZQ,RQI,VMINT)

c          do i=1, (NZ2-NZ1+1)*NRB*N2D
cc            do j=1, NANODES
cc              write(2) VMINT(j,i)
cc            end do
c            write(2) VMINT(1:NANODES,i)
c          end do

          close(2)

          deallocate(VMINT, stat=ierr)

          print*, 'partial potential surface done'
          call MPI_Barrier(MPI_COMM_WORLD, ierr)
          return

        else if (IDPOT.eq.3) then  ! read potential data 

          allocate(LINK1_ASY(3,(NZ-NZINT)*NVB*N2D), stat=ierr)
          LINK1_ASY(:,:)=0
          allocate(VCUT4D_ASY((NZ-NZINT)*NVB*N2D), stat=ierr)
          allocate(LINK2_ASY(NA1node*NA2node,(NZ-NZINT)*NVB*N2D),stat=ierr)
          LINK2_ASY(:,:)=0

          irnode=0
          icards=0        
          idatas=0
          NPL=0
          allocate(vsave(NB1node*NB2node*NA1node*NA2node),stat=ierr)
          vsave(:)=0.
          do IPOT=1, NPOTS

            open(1, file='POT'//itoc(IPOT)//'.ASY',status='UNKNOWN',form
     .        ='UNFORMATTED')
            NZ1=(NZ-NZINT)*(IPOT-1)/NPOTS+1+NZINT
            print*, 'NZ1=', NZ1
            NZ2=(NZ-NZINT)*IPOT/NPOTS+NZINT
            print*, 'NZ2=', NZ2
            allocate(VMASY(NANODES,(NZ2-NZ1+1)*NVB*N2D),stat=ierr)
            ms=0
            do j=NZ1, NZ2
              do i=1, NVB*N2D
                ms=ms+1
c                do k=1, NANODES
c                  read(1) VMASY(k,ms)
c                end do
                read(1) VMASY(1:NANODES,ms)
              end do
            end do
            close(1)

            call VLINK(41,NZ1,NZ2,NPL,NVB,irnode,icards,idatas,VrefCH(NRB+1),
     .      VMASY,vsave,LINK1_ASY,LINK2_ASY,VCUT,VCUT4D_ASY)

            NPL=irnode

            deallocate(VMASY)

          end do
          if (icards.ne.0) then
            icards=icards+1
            write(41, rec=icards) vsave
          end if

          NCARDS_ASY=icards

          allocate(VMASY_SAVE(NANODES,NCARDS_ASY), stat=ierr)

          call READPOT(41, NANODES, NCARDS_ASY, VMASY_SAVE)

          call POTPLOT(30,NZ-NZINT,NVB,N2D,NA1node*NA2node,NB1node*NB2node,
     &    ZQ(NZINT+1),RQ,VrefCH(NRB+1),Vref2D,NCARDS_ASY,LINK1_ASY,VCUT4D_ASY,
     &    VMASY_SAVE)
 
! save data
          write(43,*) NCARDS_ASY
          do j=1, (NZ-NZINT)*NVB*N2D
            write(43,*) LINK1_ASY(1,j),LINK1_ASY(2,j),LINK1_ASY(3,j)
          end do
          do j=1, (NZ-NZINT)*NVB*N2D
          do i=1, NA1node*NA2node
            write(43,*) LINK2_ASY(i,j)
          end do
          end do
          close(43)

          stop 'IDPOT 3'

        else if (IDPOT.eq.4) then

          allocate(LINK1_INT(3,NZINT*NRB*N2D), stat=ierr)
          LINK1_INT(:,:)=0
          allocate(VCUT4D_INT(NZINT*NRB*N2D), stat=ierr)
          allocate(LINK2_INT(NA1node*NA2node,NZINT*NRB*N2D), stat=ierr)
          LINK2_INT(:,:)=0

          irnode=0
          icards=0        
          idatas=0
          NPL=0
          allocate(vsave(NB1node*NB2node*NA1node*NA2node),stat=ierr)
          vsave(:)=0.

          do IPOT=1, NPOTS

            NZ1=NZINT*(IPOT-1)/NPOTS+1
            print*, 'NZ1=', NZ1
            NZ2=NZINT*IPOT/NPOTS
            print*, 'NZ2=', NZ2
            allocate(VMINT(NANODES,(NZ2-NZ1+1)*NRB*N2D),stat=ierr)
            open(1, file='POT'//itoc(IPOT)//'.INT',status='UNKNOWN',form='UNFORMATTED')
            ns=0
            do j=NZ1, NZ2
              do i=1, NRB*N2D
                ns=ns+1
c                do k=1, NANODES
c                  read(1) VMINT(k,ns)
c                end do
                read(1) VMINT(1:NANODES,ns)
              end do
            end do
            close(1)

            call VLINK(42,NZ1,NZ2,NPL,NRB,irnode,icards,idatas,
     .        VrefCH,VMINT,vsave,LINK1_INT,LINK2_INT,VCUT,VCUT4D_INT)

            NPL=irnode

            deallocate(VMINT)

          end do
          if (icards.ne.0) then
            icards=icards+1
            write(42, rec=icards) vsave
          end if

          NCARDS_INT=icards

          allocate(VMINT_SAVE(NANODES,NCARDS_INT), stat=ierr)
          call READPOT(42, NANODES, NCARDS_INT, VMINT_SAVE)
          call POTPLOT(31,NZINT,NRB,N2D,NA1node*NA2node,NB1node*NB2node,
     &    ZQ,RQI,VrefCH,Vref2D,NCARDS_INT,LINK1_INT,VCUT4D_INT,VMINT_SAVE)

! save data
          write(44,*) NCARDS_INT
          do j=1, NZINT*NRB*N2D
            write(44,*) LINK1_INT(1,j),LINK1_INT(2,j),LINK1_INT(3,j)
          end do
          do j=1, NZINT*NRB*N2D
          do i=1, NA1node*NA2node
            write(44,*) LINK2_INT(i,j)
          end do
          end do
          close(44)

          stop 'IDPOT 4'

        else if (IDPOT.eq.5 .or. IDPOT.eq.50) then ! read post-cutoff potential

          call READ_POST_CUTOFF_POTS

        end if
            
        end subroutine POTENTIALS
c========================================================================
        subroutine POT8D(ID, NZ1,NZ2,NVB,ZQ,RQ,VM)

        use COMPARAS, only : NVBRU,NVBAU,NA1node,NA2node,NB1node,NB2node
        use QUADPOINTS, only : RUQ,AUQ,A1node,A2node,B1node,B2node

        implicit none

        integer NZ1, NZ2, NVB, ID
        integer ms, i, j, j2, k, l, m, n, n2
        real*8 vx, Z, R1, RU, AU, A1, B1, A2, B2
        real*8 ZQ(1), RQ(1) 
        real*4 VM(1)
          
        print*, 'calculating potential values'

        do i=NZ1, NZ2
          print*, 'i=', i
          Z=ZQ(i)
          do j=1, NVB; R1=RQ(j)
          do j2=1, NVBRU; RU=RUQ(j2)
            do k=1, NVBAU; AU=AUQ(k)

              ms=0
              do l=1, NA1node; A1=A1node(l)
                do n=1, NA2node; A2=A2node(n)
                  do m=1, NB1node; B1=B1node(m)
                    do n2=1, NB2node; B2=B2node(n2)
                      call POTDRIV(Z,R1,RU,AU,A1,A2,B1,B2,vx)  
                      ms=ms+1
                      VM(ms)=vx
                    end do
                  end do
                end do
              end do
	      write(ID) VM(1:NA1node*NA2node*NB1node*NB2node)

            end do
          end do
          end do
        end do

        end subroutine POT8D
c=======================================================================
        subroutine IDVCUT1
        use ENGARRAYS, only : VMASY, VMASY_SAVE, VCUT4D_ASY
        use COMPARAS, only : NZ, NZINT, NVB, NRB, N2D, NA1node, NA2node,
     .    NB1node, NB2node, VCUT, NCARDS_ASY
        use SMALL_ARRAYS, only : VrefCH, Vref2D, LINK1_ASY, LINK2_ASY
        use QUADPOINTS, only : ZQ, RQ, RQI

        implicit none
        
        integer ierr, NANODES, i, j

        NANODES=NA1node*NA2node*NB1node*NB2node

        allocate(LINK1_ASY(3,(NZ-NZINT)*NVB*N2D), stat=ierr)
        LINK1_ASY(:,:)=0

        allocate(VCUT4D_ASY((NZ-NZINT)*NVB*N2D), stat=ierr)

        allocate(LINK2_ASY(NA1node*NA2node,(NZ-NZINT)*NVB*N2D),stat=ierr)
        LINK2_ASY(:,:)=0

        call VLINK(15,NCARDS_ASY,NZ-NZINT,NVB,VrefCH(NRB+1),VMASY,
     .    LINK1_ASY,LINK2_ASY,VCUT,VCUT4D_ASY)

        deallocate(VMASY)

        allocate(VMASY_SAVE(NANODES,NCARDS_ASY), stat=ierr)

        call READPOT(15, NANODES, NCARDS_ASY, VMASY_SAVE)

        call POTPLOT(30,NZ-NZINT,NVB,N2D,NA1node*NA2node,NB1node*NB2node,
     &    ZQ(NZINT+1),RQ,VrefCH(NRB+1),Vref2D,NCARDS_ASY,LINK1_ASY,VCUT4D_ASY,
     &    VMASY_SAVE)
 
! save data
        write(13,*) NCARDS_ASY

        do j=1, (NZ-NZINT)*NVB*N2D
          write(13,*) LINK1_ASY(1,j),LINK1_ASY(2,j),LINK1_ASY(3,j)
        end do
        do j=1, (NZ-NZINT)*NVB*N2D
        do i=1, NA1node*NA2node
          write(13,*) LINK2_ASY(i,j)
        end do
        end do

        end subroutine IDVCUT1
c=======================================================================
        subroutine IDVCUT2
        use ENGARRAYS, only : VMINT, VMINT_SAVE, VCUT4D_INT
        use COMPARAS, only : NZ, NZINT, NVB, NRB, N2D, NA1node, NA2node,
     .    NB1node, NB2node, VCUT, NCARDS_INT
        use SMALL_ARRAYS, only : VrefCH, Vref2D, LINK1_INT, LINK2_INT
        use QUADPOINTS, only : ZQ, RQ, RQI

        implicit none
        
        integer ierr, NANODES, i, j

        NANODES=NA1node*NA2node*NB1node*NB2node

        allocate(LINK1_INT(3,NZINT*NRB*N2D), stat=ierr)
        LINK1_INT(:,:)=0

        allocate(VCUT4D_INT(NZINT*NRB*N2D), stat=ierr)
        
        allocate(LINK2_INT(NA1node*NA2node,NZINT*NRB*N2D), stat=ierr)
        LINK2_INT(:,:)=0

        call VLINK(16,NCARDS_INT,NZINT,NRB,VrefCH,VMINT,LINK1_INT,
     .    LINK2_INT,VCUT,VCUT4D_INT)

        deallocate(VMINT)
        allocate(VMINT_SAVE(NANODES,NCARDS_INT), stat=ierr)

        call READPOT(16, NANODES, NCARDS_INT, VMINT_SAVE)

        call POTPLOT(31,NZINT,NRB,N2D,NA1node*NA2node,NB1node*NB2node,
     &    ZQ,RQI,VrefCH,Vref2D,NCARDS_INT,LINK1_INT,VCUT4D_INT,VMINT_SAVE)

! save data
        write(14,*) NCARDS_INT

        do j=1, NZINT*NRB*N2D
          write(14,*) LINK1_INT(1,j),LINK1_INT(2,j),LINK1_INT(3,j)
        end do
        do j=1, NZINT*NRB*N2D
        do i=1, NA1node*NA2node
          write(14,*) LINK2_INT(i,j)
        end do
        end do

        end subroutine IDVCUT2
c=======================================================================
        subroutine SYMMPOT(NB1node,NB2node,IDVSYM)
        implicit none
        integer NB1node,NB2node,ib2,ib1,inodbx
        integer IDVSYM(NB2node,1)

        inodbx=0
        do 2000 ib1=1,NB1node
        do 2000 ib2=1,NB2node
          inodbx=inodbx+1
          IDVSYM(ib2,ib1)=inodbx
 2000   continue

        end subroutine SYMMPOT
c======================================================================
        subroutine VLINK(kchanl,NZ1,NZ2,NPL,NRB,irnode,icards,idatas,VrefCH,VM,
     .    vsave,LINK1,LINK2,VCUT,VCUT4D)

        use COMPARAS, only : N2D, NA1node, NA2node, NB1node, NB2node
        use SMALL_ARRAYS, only : Vref2D

        implicit none

        integer NZ1, NZ2, NRB, NPL 
        real*8 VrefCH(1), VCUT, vzero
        real*4 VM(NB1node*NB2node*NA1node*NA2node,1), VCUT4D(1)
        real*4 vsave(NB1node*NB2node*NA1node*NA2node)

        integer LINK1(3,1), LINK2(NA1node*NA2node,1)

        integer j, k, l, kchanl
        integer irnode, icards, idatas


        do 10 l=NZ1, NZ2
        do 10 k=1, NRB
        do 10 j=1, N2D
          irnode=irnode+1 
          vzero=VrefCH(k)+Vref2D(j)

          VCUT4D(irnode)=VCUT-vzero
           
          call VLINK0(kchanl,icards,idatas,NB1node*NB2node,NA1node*NA2node,
     .      LINK1(1,irnode),LINK2(1,irnode),vzero,VCUT,VM(1,irnode-NPL),vsave)

          write(*,9100) l, k, j, irnode, LINK1(1,irnode), LINK1(2,irnode), LINK1(3,irnode)

10      continue        

9100    format(8i10)


        end subroutine VLINK
c==========================================================================================
       subroutine VLINK0(kchanl,icards,idatas,nodebx,nodeax,
     &   linkv11,linkv11a,vzero,vcut,vanode,vasave)
       implicit none
       integer :: kchanl,icards,idatas,nodebx,nodeax,linkv11(3),
     &            linkv11a(nodeax),nsaves,inodax,inodbx
       real*4  :: vanode(nodebx,nodeax),vasave(nodebx,nodeax)
       real*8  :: vzero,vcut
       real*4  :: vmina

        integer :: i, j
c
c......potenital saved on disk as (nanode,ncard11)
c      potential will be cut off at the level of nodebx
c      it appears in the main program as vnodes(nodebx,nodeax,ncard11)
c      linkv11(1) saves the location of ncard11
c      linkv11(2) saves the location of nodeax
c      linkv11(3) saves the potential length for quadrature (au,r1,rs)

       linkv11(1)=icards+1
       linkv11(2)=idatas+1
c       call IZERO1D(nodeax,linkv11a)
c
c......nsaves is the datas of nodebx block in nanode=nodebx*nodeax
c      idatas is the position saved on DISK
c      icards is the cards saved on DISK

       nsaves=0
       do 3000 inodax=1,nodeax
          vmina=1.0d+30
          do 1200 inodbx=1,nodebx
             vmina=min(vmina,vanode(inodbx,inodax))
 1200     continue

          if (vmina.le.vcut) then
             nsaves=nsaves+1
             idatas=idatas+1
             linkv11a(inodax)=nsaves
             do 1500 inodbx=1,nodebx
             vasave(inodbx,idatas)=min(vanode(inodbx,inodax),vcut)-vzero
 1500        continue
             if (idatas.eq.nodeax) then
                icards=icards+1
                idatas=0
c                do i=1, nodeax
c                do j=1, nodebx
                 write(kchanl,rec=icards) vasave
c                  write(600,*) vasave(j,i)
c                end do
c                end do
             endif
          endif
 3000  continue

       linkv11(3)=nsaves

       end subroutine VLINK0

c===================================================================
        subroutine READPOT(kchanl,NANODES,NCARDS,VM)
        implicit none
        integer kchanl, NANODES, NCARDS
        integer i, j
        real*4 VM(NANODES,1)

        do j=1, NCARDS
          read(kchanl,rec=j) (VM(i,j),i=1,NANODES)
        end do
        
        end subroutine READPOT
c===================================================================
       subroutine POTPLOT(kchanl,nbasrs,nbasr1,nbasau,nodeax,nodebx,
     &                           xmatrs,xmatr1,vmatr1,vmatau,
     &                           nvsave,linkvx,vcutrx,vnodes)
c
c......declare varialbes

       implicit none 
       integer :: kchanl,nbasrs,nbasr1,nbasau,nodeax,nodebx,
     &            nvsave,linkvx(3,nbasau,nbasr1,nbasrs),
     &            inodrs,inodr1,inodau,icards,idatas,npoten
       real*8  :: xmatrs(nbasrs),xmatr1(nbasr1),
     &            vmatr1(nbasr1),vmatau(nbasau), rs,r1,vzero
       real*4  :: vcutrx(nbasau,nbasr1,nbasrs),
     &            vnodes(nodebx,nodeax,nvsave),vmin, vmin1
c
c......find the minimum for (rs,r1)

       rewind(kchanl)

       do 5000 inodrs=1,nbasrs;           rs=xmatrs(inodrs)
          do 4800 inodr1=1,nbasr1;        r1=xmatr1(inodr1)
             vmin=1.0d+30
             do 4600 inodau=1,nbasau
                icards=linkvx(1,inodau,inodr1,inodrs)
                idatas=linkvx(2,inodau,inodr1,inodrs)
                npoten=linkvx(3,inodau,inodr1,inodrs)

                vzero=vmatr1(inodr1)+vmatau(inodau)
             
                vmin1=1.0d+30 
                if (npoten.eq.0) then
                   vmin1=vcutrx(inodau,inodr1,inodrs)
                else
                   call POTPLOT1(npoten*nodebx,vmin1,
     &                           vnodes(1,idatas,icards))
                endif 
                vmin=min(vmin,vmin1+vzero) 
 4600        continue 
             write(kchanl,9100)rs,r1,vmin*27.2116d0
 4800     continue
          write(kchanl,9100)
 5000  continue 
       close(kchanl) 
 9100  format(3f20.10)
       return
       end 
!
!  ***************************************************************************
!  *                                                                         *
!  *   Program    :  potplot1                                                *
!  *   Function   :  find the minimum of potential                           *
!  *                                                                         *
!  ***************************************************************************
       subroutine POTPLOT1(nvnode,vmin,vnodes)
       implicit none 
       integer :: nvnode,ivnode
       real*4  :: vmin,vnodes(nvnode) 
c
       vmin=1.0d+30  
       do 2000 ivnode=1,nvnode
          vmin=min(vmin,vnodes(ivnode)) 
 2000  continue 
       return
       end 
c===================================================================
        subroutine READ_POST_CUTOFF_POTS
        use ENGARRAYS, only : VMASY_SAVE, VMINT_SAVE,
     .    VCUT4D_ASY, VCUT4D_INT
        use COMPARAS, only : NZ, NZINT, NVB, NRB, N2D, NA1node, NA2node,
     .    NB1node, NB2node, VCUT, NCARDS_ASY, NCARDS_INT
        use SMALL_ARRAYS, only : VrefCH, Vref2D, LINK1_ASY, LINK2_ASY,
     .    LINK1_INT, LINK2_INT

        implicit none

        integer l, k, i, j, ierr, irnode, NANODES
        real*8 vzero

! read index information
        read(43,*) NCARDS_ASY
        read(44,*) NCARDS_INT

        NANODES=NA1node*NA2node*NB1node*NB2node
        allocate(LINK1_ASY(3,(NZ-NZINT)*NVB*N2D), stat=ierr)
        LINK1_ASY(:,:)=0

        allocate(VCUT4D_ASY((NZ-NZINT)*NVB*N2D), stat=ierr)

        allocate(LINK2_ASY(NA1node*NA2node,(NZ-NZINT)*NVB*N2D),stat=ierr)
        LINK2_ASY(:,:)=0

        do j=1, (NZ-NZINT)*NVB*N2D
          read(43,*) LINK1_ASY(1,j),LINK1_ASY(2,j),LINK1_ASY(3,j)
        end do
        do j=1, (NZ-NZINT)*NVB*N2D
        do i=1, NA1node*NA2node
          read(43,*) LINK2_ASY(i,j)
        end do
        end do

        allocate(LINK1_INT(3,NZINT*NRB*N2D), stat=ierr)
        LINK1_INT(:,:)=0

        allocate(VCUT4D_INT(NZINT*NRB*N2D), stat=ierr)
        
        allocate(LINK2_INT(NA1node*NA2node,NZINT*NRB*N2D), stat=ierr)
        LINK2_INT(:,:)=0

        do j=1, NZINT*NRB*N2D
          read(44,*) LINK1_INT(1,j),LINK1_INT(2,j),LINK1_INT(3,j)
        end do
        do j=1, NZINT*NRB*N2D
        do i=1, NA1node*NA2node
          read(44,*) LINK2_INT(i,j)
        end do
        end do

! calculate cutoff pots
        irnode=0
        do 10 l=NZINT+1,NZ
        do 10 k=1, NVB
        do 10 j=1, N2D
          irnode=irnode+1 
          vzero=VrefCH(NRB+k)+Vref2D(j)

          VCUT4D_ASY(irnode)=VCUT-vzero
           
10      continue

        irnode=0
        do 20 l=1,NZINT
        do 20 k=1, NRB
        do 20 j=1, N2D
          irnode=irnode+1 
          vzero=VrefCH(k)+Vref2D(j)

          VCUT4D_INT(irnode)=VCUT-vzero
           
20      continue

! read potentials
        allocate(VMASY_SAVE(NANODES,NCARDS_ASY), stat=ierr)
        call READPOT(41, NANODES, NCARDS_ASY, VMASY_SAVE)

        allocate(VMINT_SAVE(NANODES,NCARDS_INT), stat=ierr)
        call READPOT(42, NANODES, NCARDS_INT, VMINT_SAVE)

        end subroutine READ_POST_CUTOFF_POTS
c=========================================================================
        subroutine POT8D_CLUSTER(NZ1,NZ2,NVB,ZQ,RQ,VM)

        use COMPARAS,only:N2D,NVBRU,NVBAU,NA1node,NA2node,NB1node,NB2node,VINF
        use QUADPOINTS,only:RUQ,AUQ,A1node,A2node,B1node,B2node

        implicit none

        integer NZ1, NZ2, NVB, NANODES
        integer ms, ns, i, j, k, k2, l, m, n, n2, n3, n4, n5
        integer ic, jc, kc, lc, mc, nc, n2c
        integer ip, jp, kp, lp, mp, np, n2p
        integer it
        real*8 vx
        real*8 Z, R1, RU, AU, A1, B1, A2, B2
        real*8 Zp, R1p, A1p, B1p, A2p, B2p
        real*8 vvv, Vmin, VVmin
        real*8 ZQ(1), RQ(1) 
        real*4 VM(1)
        real*8 DIST0(15)
        real*8 dists(15, 81)    ! 3 ^ 4 = 81
        real*8 VTMP(81)

        NANODES=NA1node*NA2node*NB1node*NB2node

        VVmin=100.d0

        print*, 'calculating potential values'
        print*, 'NZ1=',NZ1,'NZ2=',NZ2,'NVB=',NVB
        print*
        do 10 i=NZ1, NZ2
          Z=ZQ(i)
        do 20 j=1,NVB
          R1=RQ(j)
        do 25 k=1, NVBRU
          RU=RUQ(k)
        do 30 k2=1, NVBAU
          AU=AUQ(k2)
          write(*,100) i, j, k, k2
        do 40 l=1, NA1node, 3
          lc=min(l+1,NA1node)
          A1=A1node(lc)
        do 50 n=1, NA2node, 3
          nc=min(n+1,NA2node)
          A2=A2node(nc)
        do 60 m=1, NB1node, 3
          mc=min(m+1,NB1node)
          B1=B1node(mc)
        do 70 n2=1, NB2node, 3
          n2c=min(n2+1, NB2node)
          B2=B2node(n2c)

!          call JTOD(Z,R1,RU,AU,A1,A2,B1,B2,dist0)

          ns=0
          do lp=max(lc-1,l), min(lc+1,NA1node)
            A1p=A1node(lp)
          do np=max(nc-1,n), min(nc+1,NA2node)
            A2p=A2node(np)
          do mp=max(mc-1,m), min(mc+1,NB1node)
            B1p=B1node(mp)
          do n2p=max(n2c-1,n2), min(n2c+1,NB2node)
            B2p=B2node(n2p)
            ns=ns+1
            
!            call JTOD(Z,R1,RU,AU,A1p,A2p,B1p,B2p,dists(1,ns))
            
          end do
          end do
          end do
          end do

!          call Energy0(ns,15,dist0,dists,VTMP)

          ns=0
          do lp=max(lc-1,l), min(lc+1,NA1node)
          do np=max(nc-1,n), min(nc+1,NA2node)
          do mp=max(mc-1,m), min(mc+1,NB1node)
          do n2p=max(n2c-1,n2), min(n2c+1,NB2node)
            ns=ns+1
            ms=(i-NZ1)*NVB*N2D*NANODES+(j-1)*N2D*NANODES
     .         +(k-1)*NVBAU*NANODES+(k2-1)*NANODES
     .         +(lp-1)*NA2node*NB1node*NB2node
     .         +(np-1)*NB1node*NB2node+(mp-1)*NB2node
            VM(ms+n2p)=min(VTMP(ns),-40.743d0)-VINF
          end do
          end do
          end do
          end do


 70     continue
 60     continue
 50     continue
 40     continue
 30     continue
 25     continue
 20     continue
 10     continue

100        format(X,4i5)
9100    format(f20.10,X,f20.10)

        end subroutine POT8D_CLUSTER

!============================================================================
	subroutine POT8D_MPI
        use COMPARAS,only:N2D,NVBRU,NZ,NZINT, NVBAU,NA1node,NA2node,
     &    NB1node,NB2node,NB3node,VINF,IPOT,NPOTS, NVB,NRB, IPOT,NPOTS
        use QUADPOINTS,only:RUQ,AUQ,A1node,A2node,B1node,B2node, RQ,RQI,
     &                      ZQ  ,B3node

	implicit none

	include 'mpif.h'

	integer :: myid, numprocs, status(MPI_STATUS_SIZE),Idone,Idone0
	integer :: numsend,sender, nsend, iall,iall0
	integer :: ms,ns, i,j,k, m,n, nLeft
	integer :: i1,i2, i3,i4, i5,i6,i7,i8, i9
	integer :: NGrid,NGrid1,NGrid2, NAnodes, ierr
	real*8  :: ZZ,R1, RU,AU, A1,A2, B1,B2,B3, vx   
	real*4, allocatable :: VM(:)

	character*4 :: itoc4

	NAnodes=NA1node*NA2node*NB1node*NB2node*NB3node
        allocate(VM(1:NAnodes),stat=ierr)

	myid=IPOT-1
	numprocs=NPOTS

	if(myid == 0) then
          print*, 'calculating potential values'
          print*, 'numprocs=', numprocs
	  print*, kind(i1),kind(i2)
	end if
	print*,'myrank=',myid

	NGrid1=NZINT     *NRB*NVBRU*NVBAU
	NGrid2=(NZ-NZINT)*NVB*NVBRU*NVBAU
	NGrid =NGrid1+NGrid2

	if(myid /= 0) then
          open(20000+myid*10+1,file='POT'//itoc4(myid)//'.INT',status='UNKNOWN',
     .           form='binary')
          open(20000+myid*10+2,file='POT'//itoc4(myid)//'.ASY',status='UNKNOWN',
     .           form='binary')
	end if

	call MPI_Barrier(MPI_COMM_WORLD, ierr)


	if(myid == 0) then ! master 

	  print*,'NGrid 1&2 ', NGrid1,NGrid2, NGrid
	  if(NGrid < numprocs-1) then
            print*,'Error : NGrid < numprocs-1',NGrid,numprocs-1
	    stop
	  end if

	  numsend=0
	  do i=1,numprocs-1
	    numsend=numsend+1

	    nsend=numsend
	    if(numsend > NGrid1) nsend=-(numsend-NGrid1)
            call MPI_SEND(nsend,1, MPI_INTEGER, i,numsend,MPI_COMM_WORLD, ierr)
	    write(9000,*) i,nsend,numsend
	  end do


          do 90 i = 1, NGrid1+NGrid2
             call MPI_RECV(Idone, 1, MPI_INTEGER,
     &                  MPI_ANY_SOURCE, MPI_ANY_TAG,
     &                  MPI_COMM_WORLD, status, ierr)
             sender     = status(MPI_SOURCE)

             if (numsend .lt. NGrid1+NGrid2) then        ! send another row
               numsend = numsend+1

               nsend=numsend
               if(numsend.gt.NGrid1) nsend=-(numsend-NGrid1)
               call MPI_SEND(nsend,1,MPI_INTEGER, sender,numsend,MPI_COMM_WORLD, ierr)
	       write(9000,*) sender,nsend,numsend
             else      ! Tell sender that there is no more work
               call MPI_SEND(MPI_BOTTOM,0,MPI_INTEGER, sender,0,MPI_COMM_WORLD, ierr)
             endif
90        end do
	  
	else ! worker

100       call MPI_RECV(iall0,1,MPI_INTEGER, 
     &                  0,MPI_ANY_TAG,MPI_COMM_WORLD, status, ierr)

	  iall=abs(iall0)
           if (status(MPI_TAG) .eq. 0) then
             go to 200
           else
             print*,'myid ',myid, 'iall0 ',iall0, 'tag ',status(MPI_TAG)
	     write(9100+mod(myid,8),*) 'myid ',myid, 'iall0 ',iall0, 'tag ',status(MPI_TAG)

!	     goto 150 ! uncomment this to test task-distribution

	     nleft=iall

!	     i4=mod(nleft-1,NVBAU) + 1; AU=AUQ(i4)
!	     nleft=(iall-i4)/NVBAU + 1
!
!	     i3=mod(nleft-1,NVBRU) + 1; RU=RUQ(i3)
!	     nleft=(iall-i3)/NVBRU + 1
!	     
!	     if(iall0 > 0) then ! INT
!	       i2=mod(iall-1,NRB) + 1; R1=RQI(i2)
!	       nleft=(iall-i2)/NRB + 1
!
!	       i1=nleft; ZZ=ZQ(i1)
!	       
!	       ms=(i1-1)*NRB*NVBRU*NVBAU + (i2-1)*NVBRU*NVBAU
!     &           +(i3-1)*NVBAU + i4
!	     else ! ASY
!	       i2=mod(iall-1,NVB) + 1; R1=RQ(i2)
!	       nleft=(iall-i2)/NVB + 1
!
!	       i1=nleft + NZINT; ZZ=ZQ(i1)
!
!	       ms=(i1-NZINT -1)*NVB*NVBRU*NVBAU + (i2-1)*NVBRU*NVBAU
!     &           +(i3-1)*NVBAU + i4
!	     end if

	     if(iall0 > 0) then ! INT
	       i1=(nleft-1)/(NRB*NVBRU*NVBAU) + 1; ZZ=ZQ(i1)
     	       nleft=mod(nleft-1,NRB*NVBRU*NVBAU) + 1
	       i2=(nleft-1)/(NVBRU*NVBAU) + 1; R1=RQI(i2)

	       ms=(i1 -1)*NRB + i2
	     else ! ASY
	       i1=(nleft-1)/(NVB*NVBRU*NVBAU) + 1 + NZINT; ZZ=ZQ(i1)
     	       nleft=mod(nleft-1,NVB*NVBRU*NVBAU) + 1
	       i2=(nleft-1)/(NVBRU*NVBAU) + 1; R1=RQ(i2)

	       ms=(i1-NZINT -1)*NVB + i2
	     end if
	     nleft=mod(nleft-1,NVBRU*NVBAU) + 1
	     i3=(nleft-1)/NVBAU + 1; RU=RUQ(i3)
	     nleft=mod(nleft-1,NVBAU) + 1
	     i4=nleft; AU=AUQ(i4)

	     ms=(ms-1)*NVBRU*NVBAU + (i3-1)*NVBAU + i4
	     if(ms /= iall) then
               print*,'ms/=iall',ms,iall,iall0
               stop
	     end if

!	     goto 150 ! uncomment this to test task-distribution

	     ms=0
	     do i5=1,NA1node;	A1=A1node(i5)
	     do i6=1,NA2node;	A2=A2node(i6)
	     do i7=1,NB1node;	B1=B1node(i7)
	     do i8=1,NB2node;	B2=B2node(i8)
	     do i9=1,NB3node;	B3=B3node(i9)
	       ms=ms+1
	       call POTDRIV(ZZ,R1,RU,AU, A1,A2,B1,B2,B3, vx)
	       VM(ms)=vx


	     end do; end do
	     end do; end do
	     end do;


	     if(iall0 > 0) then ! INT
	       write(20000+myid*10+1) i1,i2,i3,i4
	       write(20000+myid*10+1) VM(1:NAnodes)
	     else ! ASY
	       write(20000+myid*10+2) i1-NZINT,i2,i3,i4
	       write(20000+myid*10+2) VM(1:NAnodes)
	     end if

150          continue

             call MPI_SEND(Idone, 1, MPI_INTEGER,
     &                     0,myid,MPI_COMM_WORLD, ierr) 
             go to 100
           end if!...if (status(MPI_TAG) .eq. 0)
	end if

200     continue

!@@	print*,'mywork done',myid

	call MPI_Barrier(MPI_COMM_WORLD, ierr)
	if(myid == 0) print*,'all worker done'

        if(myid /= 0) then
          close(20000+myid*10+1)
          close(20000+myid*10+2)
        end if

        call MPI_Barrier(MPI_COMM_WORLD, ierr)

	end subroutine POT8D_MPI
!============================================================================

