c--------------------------------------------------------------------------
        module QUADPOINTS
        implicit none
        save

! quadrature points array for radial coordinates
        real(8), allocatable :: ZQ(:), RQ(:), RQI(:), RUQ(:)

! angle quadrature points arrays
! AU: umbrella vibration angle
! A1: the angle of Z and R
! A2: the angle of R and S
! B1: the azimuth angle of the rotation of YCZ3 acround the vector R
! B2: the azimuth angle of the rotation of CZ3 around vector s
        real(8), allocatable :: AUQ(:), A1node(:), B1node(:), A2node(:),
     .    B2node(:) 

        real(8), allocatable :: A1weight(:), B1weight(:), A2weight(:),
     .    B2weight(:)
	real(8), allocatable :: B2node_1(:) ,B2weight_1(:)

        end module QUADPOINTS
c---------------------------------------------------------------------------
        module TRANSMATS
        implicit none
        save

        real(4), allocatable :: TMATZ(:,:), TMATIZ(:,:), TMATIR(:,:)
        real*8,  allocatable ::   TMATR(:,:)   !###
        real*8, allocatable :: TMAT2D(:,:)    !###

        end module TRANSMATS
c---------------------------------------------------------------------------
        module ROTBASIS
        implicit none
        save

        real(4), allocatable :: DMATA1(:,:,:,:) 
	real*8,  allocatable :: CGC(:,:) ,DMATA2(:,:,:,:)
       
        real*8,  allocatable :: RTMATB1(:,:,:), RTMATB2(:,:,:)    
        complex*16, allocatable :: CTMATB1(:,:), CTMATB2(:,:)
        
! IDXJJ: to save (j1, j2, j12, ik2) pair, kb2=KB2(ik2)
! TABJJ: to save index for (j12,j1,j2,kb2,kbz)
        integer, allocatable :: IDKBZ(:,:), KB2(:), IDXCZ3(:,:), 
     .    IDXJJ(:,:), TABJJ(:,:,:,:,:)

! rotation basis for solving YCZ3 bound states
        integer, allocatable :: IDXJJ14(:,:)
        real*8, allocatable ::  CGC14(:,:)     !###
!	real*8, allocatable :: TTROT(:,:),ROT_ENG(:)

        end module ROTBASIS
c---------------------------------------------------------------------------
        module CPBASIS
        implicit none
        save

        integer, allocatable :: idxCP1(:,:), idxCP2(:,:)

        real(4), allocatable :: CPVAL(:), CPMAT(:)
 
        end module CPBASIS
c---------------------------------------------------------------------------
        module ENGARRAYS
        implicit none
        save

        real*4, allocatable :: VMASY(:,:), VMINT(:,:)

! the two arrays save potential values after cut off
        real*4, allocatable :: VMASY_SAVE(:,:), VMINT_SAVE(:,:)

        real*4, allocatable :: VCUT4D_ASY(:), VCUT4D_INT(:)

        real*8, allocatable :: VMAT14(:)      !###

        real*4, allocatable :: ROT1(:), ROT2(:)

        complex*8, allocatable :: CENT1(:), CENT2(:)

        real*4, allocatable :: KROTJJ(:,:)

        real*4, allocatable :: FZABS(:),FRABS(:)

        real*8, allocatable :: ESK0(:), AE0(:)

        real*8 :: EIGEN_ENG14  !###

        end module ENGARRAYS
c--------------------------------------------------------------------------
        module WFARRAYS
        implicit none
        save
        
        complex*8, allocatable :: WA(:), WI(:), CWK0(:), PHIER(:,:,:)

        real*8,    allocatable :: EIGEN_VEC14(:)        ! wf for YCZ3 

        end module WFARRAYS
c--------------------------------------------------------------------------
        module REDUCED_MASS
        implicit none
        save

        real(8) :: TMS_YCZ3, RMS_YCZ3, TMS_3Z, RMS_CZ3

        end module REDUCED_MASS
c--------------------------------------------------------------------------
        module COMPARAS
        implicit none
        save

! atom's name and mass
        character(5) :: ATOM(3) ! Y:ATOM(1), C:ATOM(2), Z:ATOM(3)
        real*8 AM(3)

! parameters for Z freedom
        integer :: NZ, NZINT
        real(8) :: Zmin, Zmax

! parameters for R freedom
        integer :: NR, NVB, NRB
        real(8) :: Rmin, Rmax

! parameters for RU freedom
        integer :: NRU, NVBRU
        real(8) :: RUmin, RUmax

! parameters for AU freedom
        integer :: NAU, NVBAU
        real(8) :: AUmin, AUmax

        integer :: N2D  ! the product of NVBRU and NVBAU

! parameters for angular quadrature points
        integer :: NA1node, NA2node, NB1node, NB2node

! parameters for rotation numbers
        integer :: Jtot, KBZmin, KBZmax, parity, parity5 
        integer :: J12max, JS1, JD1, J1max, JS2, JD2, J2max, KB2max
	integer :: nfoldb2

        integer :: NKBZ
        integer :: NKB2, NKB20  
        integer :: M1max
        integer :: NBASCZ3
        integer :: NBASJJ  ! total coupled rotation basis
        integer :: NBASJJ14 ! coupled rotation basis for YCZ3

! parameters for centrifugal potential
        integer :: NDIMCP, NDIMCP2
        integer :: NCP_COUPLED

! parameters for equlibrium conformation  
        real(8) :: Zinf, R0, RU0, AU0, A10, A20, B10, B20   !RU0:equlibrium bond distance of C-Z
        real(8) :: VINF=0 ! potential value for equlibrium conformation

! parameters for initial wavefunction
        integer :: KBZinit, J1init, istat1, NSTEP0, NKBZinit
        real(8) :: Zinit, ESKinit, DELTA

! parameters for scattering energy
        integer :: NE0
        real(8) :: E1, E2

! parameters for wavefunction propagation and flux calculation
        real(8) :: TTOTAL, DT, TPRINT, T_ABSZ, T_ABSR, RFLUX, VCUT

! parameters for absorption potential
        real(8) :: ZABS, CZABS, PZABS, RABS, CRABS, PRABS 
        integer :: NZABS, NRABS

! parameters for potential calculation
        integer :: IDPOT        ! 1 = full, 2 = partial, 3 = readin     
        integer :: NPOTS         ! the best way: NPOTS = OMP_NUM_THREADS
        integer :: IPOT          ! the block for calculation
        integer :: PES_TYPE      ! 1: JG; 2: Ab Initio

        integer :: NCARDS_ASY, NCARDS_INT

        end module COMPARAS
c--------------------------------------------------------------------------
        module SMALL_ARRAYS
        implicit none
        save

        real(8), allocatable :: ETASY(:), ETINT(:)
        real(8), allocatable :: VPOT(:),RQ0(:),EIG0(:),HAR(:,:),VrefCH(:),COFFp(:)
        !real(8), allocatable :: VPOT_RU(:),RUQ0(:),EIGRU(:),VrefRU(:),HAR_RU(:,:)
        !real(8), allocatable :: VPOT_AU(:),AUQ0(:),EIGAU(:),VrefAU(:),HAR_AU(:,:)
        real(8), allocatable :: EIG2D(:),Vref2D(:),HAR2D(:,:)
        real(8), allocatable :: COEFFAU(:,:), COEFFRU(:,:)
        real(4), allocatable ::  TMATAU(:,:),TMATRU(:,:)

! the following arrays will be used in the whole calculation
        real*8,  allocatable :: BIA_CZ3(:), BIC_CZ3(:)    !###
        real(4), allocatable :: BROTZ(:)
	real*8,  allocatable :: BROTR(:)  
        real(4), allocatable :: WVIBR(:,:)

        integer, allocatable :: LINK1_ASY(:,:), LINK2_ASY(:,:)
        integer, allocatable :: LINK1_INT(:,:), LINK2_INT(:,:)
        integer, allocatable :: IDVSYM(:)
        
        end module SMALL_ARRAYS
