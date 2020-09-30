INTEGER, PARAMETER :: ntr = 1, nconsume = 1
!  INTEGER,PARAMETER :: NI=24,  NJ=24,  NK=24, ngrid=3, maxout=20832,   maxint=15768,   int1=13824
!  INTEGER,PARAMETER :: NI=24,  NJ=48,  NK=24, ngrid=4, maxout=39992,   maxint=31590,   int1=127648
   
!  INTEGER,PARAMETER :: NI=48,  NJ=24,  NK=32, ngrid=4, maxout=13272,   maxint=9360,    int1=8192
!  INTEGER,PARAMETER :: NI=48,  NJ=24,  NK=32, ngrid=5, maxout=46416,   maxint=37448,   int1=32768
!  INTEGER,PARAMETER :: NI=48,  NJ=32,  NK=16, ngrid=4, maxout=36312,   maxint=28080,   int1=24576
!  INTEGER,PARAMETER :: NI=48,  NJ=48,  NK=24, ngrid=4, maxout=76352,   maxint=63180,   int1=55296
!  INTEGER,PARAMETER :: NI=48,  NJ=48,  NK=32, ngrid=4, maxout=99512 ,  maxint =84240,  int1=73728
!  INTEGER,PARAMETER :: NI=48,  NJ=96,  NK=24, ngrid=4, maxout=149072,  maxint=126360,  int1=110592
!  INTEGER,PARAMETER :: NI=48,  NJ=96,  NK=32, ngrid=5, maxout=194472,  maxint=168516,  int1=147456

!  INTEGER,PARAMETER :: NI=96,  NJ=96,  NK=24, ngrid=4, maxout=291092,  maxint=252720,  int1=221184
!  INTEGER,PARAMETER :: NI=96,  NJ=192, NK=24, ngrid=4, maxout=575132,  maxint=505440,  int1=442368
!  INTEGER,PARAMETER :: NI=96,  NJ=192, NK=32, ngrid=5, maxout=750240,  maxint=674064,  int1=589824
!  INTEGER,PARAMETER :: NI=96,  NJ=192, NK=64, ngrid=6, maxout=1449264, maxint=1348164, int1=1179648
!  INTEGER,PARAMETER :: NI=96,  NJ=320, NK=48, ngrid=5, maxout=1823832, maxint=1685160, int1=1474560
!  INTEGER,PARAMETER :: NI=96,  NJ=384, NK=32, ngrid=5, maxout=1491264, maxint=1348128, int1=1179648
!  INTEGER,PARAMETER :: NI=96,  NJ=480, NK=32, ngrid=5, maxout=1861776, maxint=1685160, int1=1474560
!  INTEGER,PARAMETER :: NI=96,  NJ=480, NK=48, ngrid=5, maxout=2729032, maxint=2527740, int1=2211840

!  INTEGER,PARAMETER :: NI=192, NJ=144, NK=32, ngrid=5, maxout=1116288, maxint=1011096, int1=884736
!  INTEGER,PARAMETER :: NI=192, NJ=192, NK=32, ngrid=5, maxout=1482336, maxint=1348128, int1=1179648
!  INTEGER,PARAMETER :: NI=192, NJ=320, NK=48, ngrid=5, maxout=3603852, maxint=3370320, int1=2949120
!  INTEGER,PARAMETER :: NI=192, NJ=384, NK=32, ngrid=5, maxout=2946528, maxint=2696256, int1=2359296	 
!  INTEGER,PARAMETER :: NI=192, NJ=640, NK=32, ngrid=5, maxout=4898784, maxint=4493760, int1=3932160

!  INTEGER,PARAMETER :: NI=384, NJ=384, NK=32, ngrid=4, maxout=5854352, maxint=5391360, int1=4718592
!  INTEGER,PARAMETER :: NI=768, NJ=384, NK=64, ngrid=4, maxout=22553792,maxint=21565440,int1=18874368
!  INTEGER,PARAMETER :: NI=1536,NJ=768, NK=64, ngrid=4, maxout=89804672,maxint=86261760,int1=75497472
       
! =====================================================================================================     
!  INTEGER,PARAMETER :: NI=96,  NJ=96, NK=32,  ngrid=4, maxout=379472,  maxint=336960,  int1=294912  ! good      
!  INTEGER,PARAMETER :: NI=96,  NJ=96, NK=64,  ngrid=4, maxout=732992 , maxint=673920,  int1=589824  ! IN USE leewaves
  INTEGER,PARAMETER :: NI=96,  NJ=96, NK=128, ngrid=4, maxout=1440032, maxint=1347840, int1=1179648 ! good
!  INTEGER,PARAMETER :: NI=96,  NJ=96, NK=256, ngrid=4, maxout=2854112, maxint=2695680, int1=2359296 ! good
!  INTEGER, PARAMETER :: NI=96,  NJ=192,NK=256, ngrid=6, maxout=5642688, maxint=5392656, int1=4718592   
   
   
!  INTEGER,PARAMETER :: NI=192, NJ=48, NK=256, ngrid=4, maxout=2887112, maxint=2695680, int1=2359296   
!  INTEGER,PARAMETER :: NI=192, NJ=96, NK=64,  ngrid=4, maxout=1448432, maxint=1347840, int1=1179648 
!  INTEGER,PARAMETER :: NI=192, NJ=96, NK=256, ngrid=4, maxout=5640272, maxint=5391360, int1=4718592 ! Large residue
!  INTEGER,PARAMETER :: NI=192, NJ=96, NK=256, ngrid=5, maxout=5642288, maxint=5392512, int1=4718592 ! Large residue
!  INTEGER,PARAMETER :: NI=192, NJ=96, NK=256, ngrid=6, maxout=5642688, maxint=5392656, int1=4718592 ! Large residue
      