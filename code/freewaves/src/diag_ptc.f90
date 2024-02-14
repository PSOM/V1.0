subroutine diag_ptc(ptc)

   USE header

   implicit none
   integer i,j,k
   REAL(kind=rc_kind), dimension(0:NI+1,0:NJ+1, 0:NK+1)  :: ptc

   if(.NOT.(use_Shchepetkin)) then ! calculates baroclinic press gradients
       CALL rpevalgrad_Song_short(0)
   else
       CALL rpevalgrad_Sche(0);
   endif
    
   do k=1, NK
      do j=1, NJ
      do i=1, NI
         ptc(i,j,k)=rp(i,j,k)*P1 + h(i,j)*gpr*10.*R0
!        ptc(i,j,k)=rp(i,j,k)*P1 + h(i,j)*gpr*10.*R0 + qpr*P1*p(i,j,k)
      end do
      end do
      ptc(     0,1:NJ,k)=ptc(     1,1:NJ,k)
      ptc(  NI+1,1:NJ,k)=ptc(    NI,1:NJ,k)
      ptc(0:NI+1,   0,k)=ptc(0:NI+1,   1,k)
      ptc(0:NI+1,NJ+1,k)=ptc(0:NI+1,  NJ,k)
   end do
   ptc(0:NI+1,0:NJ+1,   0)=ptc(0:NI+1,0:NJ+1, 1)
   ptc(0:NI+1,0:NJ+1,NK+1)=ptc(0:NI+1,0:NJ+1,NK)

   return
end
