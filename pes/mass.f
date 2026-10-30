c=========================================================================
!  ***************************************************************************
!  *                                                                         *
!  *   Program    :  ATOMASS                                                 *
!  *   Function   :  output atom mass in au by input element symbol          *
!  *                                                                         *
!  ***************************************************************************
       subroutine ATOMASS(atom,amass)
       implicit none 
       character :: atom*5,atom1*2,atom2*3,element(57)*2
       integer   :: lenatom,ktype,idata,iascii,iatom,imass,jmass,
     &              kmass(4,57)
       real*8    :: amass,amunit,xmass(4,57)
c
c......datas 

      data element/'H ', 'D ', 'T ','HE','LI','BE','B ','C ','N ','O ',
     $             'F ','NE',
     &             'NA','MG','AL','SI','P ','S ','CL','AR','K ','CA',
     &             'SC','TI','V ','CR','MN','FE','CO','NI','CU','ZN',
     &             'GA','GE','AS','SE','BR','KR','RB','SR','Y ','ZR',
     &             'NB','MO','TC','RU','RH','PD','AG','CD','IN','SN',
     &             'SB','TE','I ','XE','HU'/

       data kmass
     $        / 1,  2,  3,  0,     4,  3,  0,  0,     7,  6,  0,  0,
     $          9,  0,  0,  0,    11, 10,  0,  0,    12, 13, 14,  0,
     $         14, 15,  0,  0,    16, 18, 17,  0,    19,  0,  0,  0,
     $         20, 22, 21,  0,    23,  0,  0,  0,    24, 26, 25,  0,
     $         27,  0,  0,  0,    28, 29, 30,  0,    31,  0,  0,  0,
     $         32, 34, 33, 36,    35, 37,  0,  0,    40, 36, 38,  0,
     $         39, 41, 40,  0,    40, 44, 42, 48,    45,  0,  0,  0,
     $         48, 46, 47, 49,    51, 50,  0,  0,    52, 53, 50, 54,
     $         55,  0,  0,  0,    56, 54, 57, 58,    59,  0,  0,  0,
     $         58, 60, 62, 61,    63, 65,  0,  0,    64, 66, 68, 67,
     $         69, 71,  0,  0,    74, 72, 70, 73,    75,  0,  0,  0,
     $         80, 78, 82, 76,    79, 81,  0,  0,    84, 86, 82, 83,
     $         85,  0,  0,  0,    88, 84, 86, 87,    89,  0,  0,  0,
     $         90, 91, 92, 94,    93,  0,  0,  0,    98, 92, 95, 96,
     $          0,  0,  0,  0,   102, 99,100,104,   103,  0,  0,  0,
     $        106,104,105,108,   107,109,  0,  0,   114,110,111,112,
     $        115,113,  0,  0,   118,116,117,119,   121,123,  0,  0,
     $        130,125,126,128,   127,  0,  0,  0,   132,129,131,134,
     $          5,  6,  7,  8/

       data xmass/
C      H-1             H-2             H-3
     $  1.007825037D0,  2.014101787D0,  3.016049286D0, 0.0d0,
C      HE-4            HE-3
     $  4.00260325D0,   3.016029297D0, 0.0d0,           0.0d0,
C      LI-7            LI-6
     $  7.0160045D0,    6.0151232D0,   0.0d0,           0.0d0,
C      BE-9
     $  9.0121825D0,   0.0d0,           0.0d0,           0.0d0,
C       B-11            B-10
     $ 11.0093053D0,   10.0129380D0,   0.0d0,           0.0d0,
C       C-12            C-13            C-14
     $ 12.000000000D0, 13.003354839D0, 14.003241D0,    0.0d0,
C       N-14            N-15
     $ 14.003074008D0, 15.000108978D0, 0.0d0,           0.0d0,
C       O-16            O-18            O-17
     $ 15.99491464D0,  17.99915939D0,  16.9991306D0,   0.0d0,
C       F-19
     $ 18.99840325D0,  0.0d0,           0.0d0,           0.0d0,
C      NE-20           NE-22           NE-21
     $ 19.9924391D0,   21.9913837D0,   20.9938453D0,   0.0d0,
C      NA-23
     $ 22.9897697D0,   0.0d0,           0.0d0,           0.0d0,
C      MG-24           MG-26           MG-25
     $ 23.9850450D0,   25.9825954D0,   24.9858392D0,   0.0d0,
C      AL-27
     $ 26.9815413D0,   0.0d0,           0.0d0,           0.0d0,
C      SI-28           SI-29           SI-30
     $ 27.9769284D0,   28.9764964D0,   29.9737717D0,   0.0d0,
C       P-31
     $ 30.9737634D0,   0.0d0,           0.0d0,           0.0d0,
C       S-32            S-34            S-33            S-36
     $ 31.9720718D0,   33.96786774D0,  32.9714591D0,   35.9670790D0,
C      CL-35           CL-37
     $ 34.968852729D0, 36.965902624D0, 0.0d0,           0.0d0,
C      AR-40           AR-36           AR-38
     $ 39.9623831D0,   35.967545605D0, 37.9627322D0,    0.0d0,
C       K-39            K-41            K-40
     $ 38.9637079D0,   40.9618254D0,   39.9639988D0,   0.0d0,
C      CA-40           CA-44           CA-42           CA-48
     $ 39.9625907D0,   43.9554848D0,   41.9586218D0,   47.952532D0,
C      SC-45
     $ 44.9559136D0,   0.0d0,           0.0d0,           0.0d0,
C      TI-48           TI-46           TI-47           TI-49
     $ 47.9479467D0,   45.9526327D0,   46.9517649D0,   48.9478705D0,
C       V-51            V-50
     $ 50.9439625D0,   49.9471613D0,   0.0d0,           0.0d0,
C      CR-52           CR-53           CR-50           CR-54
     $ 51.9405097D0,   52.9406510D0,   49.9460463D0,   53.9388822D0,
C      MN-55
     $ 54.9380463D0,   0.0d0,           0.0d0,           0.0d0,
C      FE-56           FE-54           FE-57           FE-58
     $ 55.9349393D0,   53.9396121D0,   56.9353957D0,   57.9332778D0,
C      CO-59
     $ 58.9331978D0,   0.0d0,           0.0d0,           0.0d0,
C      NI-58           NI-60           NI-62           NI-61
     $ 57.9353471D0,   59.9307890D0,   61.9283464D0,   60.9310586D0,
C      CU-63           CU-65
     $ 62.9295992D0,   64.9277924D0,   0.0d0,           0.0d0,
C      ZN-64           ZN-66           ZN-68           ZN-67
     $ 63.9291454D0,   65.9260352D0,   67.9248458D0,   66.9271289D0,
C      GA-69           GA-71
     $ 68.9255809D0,   70.9247006D0,   0.0d0,           0.0d0,
C      GE-74           GE-72           GE-70           GE-73
     $ 73.9211788D0,   71.9220800D0,   69.9242498D0,   72.9234639D0,
C      AS-75
     $ 74.9215955D0,   0.0d0,           0.0d0,           0.0d0,
C      SE-80           SE-78           SE-82           SE-76
     $ 79.9165205D0,   77.9173040D0,   81.916709D0,    75.9192066D0,
C      BR-79           BR-81
     $ 78.9183361D0,   80.916290D0,    0.0d0,           0.0d0,
C      KR-84           KR-86           KR-82           KR-83
     $ 83.9115064D0,   85.910614D0,    81.913483D0,    82.914134D0,
C      RB-85
     $ 84.9117d0,     0.0d0,            0.0d0,            0.0d0,
C      Sr-88           Sr-84           Sr-86           Sr-87
     $ 87.9056d0,      83.9134d0,      85.9094d0,      86.9089d0,
C      Y-89
     $ 88.9054d0,     0.0d0,            0.0d0,            0.0d0,
C      Zr-90           Zr-91           Zr-92           Zr-94
     $ 89.9043d0,      90.9053d0,      91.9046d0,      93.9061d0,
C      Nb-93
     $ 92.9060d0,     0.0d0,            0.0d0,            0.0d0,
C      Mo-98           Mo-92           Mo-95           Mo-96
     $ 97.9055d0,      91.9063d0,      94.90584d0,     95.9046d0,
C      Tc-99, longest lived isotope
     $ 98.9063d0,        0.0d0,            0.0d0,            0.0d0,
C      Ru-102          Ru-99           Ru-100          Ru-104
     $ 101.9037d0,     98.9061d0,      99.9030d0,      103.9055d0,
C      Rh-103
     $ 102.9048d0,    0.0d0,            0.0d0,            0.0d0,
C      Pd-106          Pd-104           Pd-105         Pd-108
     $ 105.9032d0,     103.9036d0,      104.9046d0,    107.90389d0,
C      Ag-107          Ag-109
     $ 106.90509d0,    108.9047d0,     0.0d0,            0.0d0,
C      Cd-114          Cd-110           Cd-111         Cd-112
     $ 113.9036d0,     109.9030d0,      110.9042d0,    111.9028d0,
C      In-115          In-113
     $ 114.9041d0,     112.9043d0,     0.0d0,            0.0d0,
C      Sn-118          Sn-116           Sn-117         Sn-119
     $ 117.9018d0,     115.9021d0,      116.9031d0,    118.9034d0,
C      Sb-121          Sb-123
     $ 120.9038d0,     122.9041d0,     0.0d0,            0.0d0,
C      Te-130          Te-125           Te-126         Te-128
     $ 129.9067d0,     124.9044d0,      125.9032d0,    127.9047d0,
C      I-127
     $ 126.9004d0,    0.0d0,            0.0d0,            0.0d0,
C      Xe-132          Xe-129           Xe-131         Xe-134
     $ 131.9042d0,     128.9048d0,      130.9051d0,    133.9054d0,
c      Hu---------large atom mass 
     $ 1.0d+4,         1.0d+5,          1.0D+6,        1.0D+7/
c
c......the ratio between amu and au 

       amunit=1822.7d0 
c
c......split element symbol and atom mass

       atom=adjustl(atom) 
       lenatom=len_trim(atom)
       do 1200 idata=1,lenatom
          iascii=ichar(atom(idata:idata))
          if (iascii.le.57) goto 1300 
 1200  continue 
       ktype=1
       if (lenatom.gt.2) goto 8100 
       atom1=atom(1:lenatom)
       goto 1400 
 1300  continue 
       ktype=2
       if (idata.ne.2.and.idata.ne.3) goto 8100   
       atom1=atom(1:idata-1)
       atom2=atom(idata:lenatom)
c
c......check atom symbol

 1400  continue 
       atom1=adjustl(atom1) 
       lenatom=len_trim(atom1)
       do 1500 idata=1,lenatom
          iascii=ichar(atom1(idata:idata))
          if (iascii.ge.97) atom1(idata:idata)=char(iascii-32)
          if (iascii.ge.65.and.iascii.le.90)  goto 1500  ! A ~ Z
          if (iascii.ge.97.and.iascii.le.122) goto 1500  ! a ~ z
          goto 8100 
 1500  continue 
       if (ktype.eq.1) goto 1700 
c
c......now check atom mass

       atom2=adjustl(atom2)
       lenatom=len_trim(atom2)
       imass=0 
       do 1600 idata=1,lenatom
          iascii=ichar(atom2(idata:idata))
          if (iascii.gt.57) goto 8100 
          iascii=iascii-48 
          imass=imass*10+iascii
 1600  continue 
c
c......which atom 

 1700  continue 
       iatom=0
 2000  continue 
       iatom=iatom+1
       if (iatom.gt.55) then
          write(*,*)'such an element not found in this lib'
          write(7,*)'such an element not found in this lib'
          stop 
       endif 
       if (atom1.ne.element(iatom)) goto 2000 
c
c......now given atomic mass

       if (ktype.eq.1) then 
          amass=xmass(1,iatom) 
          goto 6000 
       else
          do 3000 jmass=1,4
            if (imass.eq.kmass(jmass,iatom)) then 
               amass=xmass(jmass,iatom) 
               goto 6000 
            endif 
 3000     continue 
          write(*,*)'such an element not found in this lib'
          write(7,*)'such an element not found in this lib'
          stop
       endif 
c
c......to au mass

 6000  continue 
       amass=amass*amunit
       return
c
c......error found

 8100  continue
       write(*,*)'illegal character found in',atom
       write(7,*)'illegal character found in',atom
       stop
       end       

