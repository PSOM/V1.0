subroutine tracerinit(stepl) 
  USE header
  integer :: i,j,k,stepl
!     initializes tracer fields                                         
!     TRACER 1                                                          
!     =========                                                         
      it=1 


IF(stepl==0) then
  Tr = 0d0
ENDIF

IF(stepl==0) then

! if(ABS(user6-1.)<0.1) then
!       do k=1,NK; do j=1,NJ; do i=1,NI;
!         if(j>100) then
!            Tr(it,i,j,k,0) =0d0
!           else
!            Tr(it,i,j,k,0) =1d0
!         endif
!      enddo;enddo;enddo;
! endif


! if(ABS(user6-2.)<0.1) then
!       do k=1,NK; do j=1,NJ; do i=1,NI;
!         if((zc(i,j,k)*DL)>-80.) then
!            Tr(it,i,j,k,0) =0d0
!           else
!            Tr(it,i,j,k,0) =1d0
!         endif
!      enddo;enddo;enddo;
! endif


! if(ABS(user6-3.)<0.1) then
!       do k=1,NK; do j=1,NJ; do i=1,NI;
!         if(j>96.) then
!            Tr(it,i,j,k,0) =0d0
!           else
!            Tr(it,i,j,k,0) =1d0
!         endif
!      enddo;enddo;enddo;
! endif

!   PRINT*,"user6 ?",user6
! if(ABS(user6-6.)<0.1) then
!   PRINT*,"user6 !",user6
!       do k=1,NK; do j=1,NJ; do i=1,NI;
!         if(j<98. .AND. (zc(i,j,k)*DL)<-80.) then
!            Tr(it,i,j,k,0) =1d0
!           else
!            Tr(it,i,j,k,0) =0d0
!         endif
!      enddo;enddo;enddo;
! endif

!   PRINT*,"user6 ?",user6
! if(ABS(user6-7.)<0.1) then
!   PRINT*,"user6 !",user6
!       do k=1,NK; do j=1,NJ; do i=1,NI;
!         if(j<98. .AND. (zc(i,j,k)*DL-zf(i,j,0)*DL)<20.) then
!            Tr(it,i,j,k,0) =1d0
!           else
!            Tr(it,i,j,k,0) =0d0
!         endif
!      enddo;enddo;enddo;
! endif

!   PRINT*,"user6 ?",user6
! if(ABS(user6-8.)<0.1) then
!   PRINT*,"user6 !",user6
!       do k=1,NK; do j=1,NJ; do i=1,NI;
!         if(j>118. .AND. (zc(i,j,k)*DL)>-30.) then
!            Tr(it,i,j,k,0) =1d0
!           else
!            Tr(it,i,j,k,0) =0d0
!         endif
!      enddo;enddo;enddo;
! endif

!   PRINT*,"user6 ?",user6
! if(ABS(user6-9.)<0.1) then
!   PRINT*,"user6 !",user6
!       do k=1,NK; do j=1,NJ; do i=1,NI;
!         if(j<110. .AND. (zc(i,j,k)*DL-zf(i,j,0)*DL)<20.) then
!            Tr(it,i,j,k,0) =1d0
!           else
!            Tr(it,i,j,k,0) =0d0
!         endif
!      enddo;enddo;enddo;
! endif

!   PRINT*,"user6 ?",user6
! if(ABS(user6-10.)<0.1) then
!   PRINT*,"user6 !",user6
!       do k=1,NK; do j=1,NJ; do i=1,NI;
!         if((zc(i,j,k)*DL)<-80.) then
!            Tr(it,i,j,k,0) =1d0
!           else
!            Tr(it,i,j,k,0) =0d0
!         endif
!      enddo;enddo;enddo;
! endif


!   PRINT*,"user6 ?",user6
! if(ABS(user6-11.)<0.1) then
!   PRINT*,"user6 !",user6
        do k=1,NK; do j=1,NJ; do i=1,NI;
          if(((zc(i,j,k)-zf(i,j,0))*DL)<20.) then
             Tr(it,i,j,k,0) =1d0
            else
             Tr(it,i,j,k,0) =0d0
          endif
       enddo;enddo;enddo;
! endif




ENDIF

! if(stepl==50) then

!   if(ABS(user6-4.)<0.1) then
!         do k=1,NK; do j=1,NJ; do i=1,NI;
!           if(j<100 .AND. j>60. .AND. (zc(i,j,k)*DL)<-70. .AND. vor(i,j,k)/MAXVAL(vor(1:NI,60:100,k))>0.1) then
!              Tr(it,i,j,k,0) =1d0
!             else
!              Tr(it,i,j,k,0) =0d0
!           endif
!        enddo;enddo;enddo;
!   endif
! endif


! if(stepl==50) then

!   if(ABS(user6-5.)<0.1) then
!         do k=1,NK; do j=1,NJ; do i=1,NI;
!           if(j<100 .AND. j>80. .AND. (zc(i,j,k)*DL)<-70. .AND. (vor(i,j,k))/MAXVAL(vor(1:NI,80:100,k))<-0.1) then
!              Tr(it,i,j,k,0) =1d0
!             else
!              Tr(it,i,j,k,0) =0d0
!           endif
!        enddo;enddo;enddo;
!   endif
! endif



end subroutine tracerinit
