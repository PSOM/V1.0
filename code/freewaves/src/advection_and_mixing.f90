subroutine advection_and_mixing(m,n,dtimel,step) 
#include "cppdefs.h"

   USE header
   INTEGER :: i,j,k,ix
   INTEGER :: n,m,step
   REAL(kind=rc_kind) :: dtimel, res_term
   REAL(kind=rc_kind), dimension(0:NI+1,0:NJ+1,0:NK+1) :: var,var0            ! var=(s,T,u,v,w,Tr(it,:,:,:)
   REAL(kind=rc_kind), dimension(0:NI+1,0:NJ+1,0:NK+1) :: uvarx !, uvarx_advec  ! uvarx is the advective term
   REAL(kind=rc_kind), dimension(1:NI,  1:NJ,  1:NK  ) :: vardif              ! vardif is dissipation term
   REAL(kind=rc_kind), dimension(       0:NJ+1,0:NK+1) :: Tbar,sbar,ubar
   INTEGER :: av_comp(5+ntr)
   ! REAL(kind=rc_kind), dimension(1:NI,1:NJ,1:NK) :: T_h=0.,u_w=0.,v_w=0.
   
   ! ---------------------------------------------------------------
   ! By default, all variables will be subject to the following loop
   av_comp=1
#ifdef rhoonly ! rho is stored in T, s=35 and does not need to be updated.
   av_comp(2)=0
#endif

!==================================================================================================
!                                        Loop for variables 
!==================================================================================================
   do selvar=1,5+ntr
      uvarx = 0.d0; vardif = 0.d0 ! uvarx_advec = 0.d0; ! uvarx_advec is not needed
   
#ifdef implicit
      mat_A = 0.d0; mat_B = 0.d0; mat_C = 0.d0; mat_D = 0.d0; mat_test = 0.d0;
#endif

      if(av_comp(selvar)==1) then    ! skip s/T if rhoonly
         if(selvar==1) then; var=T(:,:,:,m); var0(:,:,:) = T(:,:,:,0); endif
         if(selvar==2) then; var=s(:,:,:,m); var0(:,:,:) = s(:,:,:,0); endif
         if(selvar==3) then; var=u(:,:,:,m); var0(:,:,:) = u(:,:,:,0); endif
         if(selvar==4) then; var=v(:,:,:,m); var0(:,:,:) = v(:,:,:,0); endif
         if(selvar==5) then; var=w(:,:,:,m); var0(:,:,:) = w(:,:,:,0); endif
         do it=1,ntr
            if(selvar==5+it) then; var=Tr(it,:,:,:,m); endif
         enddo

      ! ---------------------------------------------------------------
      ! computation of the advective terms
         CALL advect(var,uvarx)    ! use QUICK scheme
#ifdef extra_output
         if(selvar==1) then; term_advcT(1:NI,1:NJ,1:NK)=uvarx(1:NI,1:NJ,1:NK)*Jacinv(1:NI,1:NJ,1:NK)/TL;    endif
         if(selvar==3) then; term_advcU(1:NI,1:NJ,1:NK)=uvarx(1:NI,1:NJ,1:NK)*Jacinv(1:NI,1:NJ,1:NK)/TL*UL; endif
         if(selvar==4) then; term_advcV(1:NI,1:NJ,1:NK)=uvarx(1:NI,1:NJ,1:NK)*Jacinv(1:NI,1:NJ,1:NK)/TL*UL; endif
         if(selvar==5) then; term_advcW(1:NI,1:NJ,1:NK)=uvarx(1:NI,1:NJ,1:NK)*Jacinv(1:NI,1:NJ,1:NK)/TL*WL; endif
#endif    
      ! ---------------------------------------------------------------
      ! computation of the vertical dissipation terms
         if(selvar<4.5) then                ! vertical mixing ought to be for all tracers
            vardif=0.;
            call mixing_vertical(var,vardif,m,step)
#ifdef extra_output
            if(selvar==1) then; term_dissTv(1:NI,1:NJ,1:NK)=vardif(1:NI,1:NJ,1:NK)*Jacinv(1:NI,1:NJ,1:NK)/TL;    endif
            if(selvar==3) then; term_dissUv(1:NI,1:NJ,1:NK)=vardif(1:NI,1:NJ,1:NK)*Jacinv(1:NI,1:NJ,1:NK)/TL*UL; endif
            if(selvar==4) then; term_dissVv(1:NI,1:NJ,1:NK)=vardif(1:NI,1:NJ,1:NK)*Jacinv(1:NI,1:NJ,1:NK)/TL*UL; endif
            if(selvar==5) then; term_dissWv(1:NI,1:NJ,1:NK)=vardif(1:NI,1:NJ,1:NK)*Jacinv(1:NI,1:NJ,1:NK)/TL*WL; endif
#endif
            uvarx(1:NI,1:NJ,1:NK)=uvarx(1:NI,1:NJ,1:NK)-vardif(1:NI,1:NJ,1:NK)
         endif
               
      ! ---------------------------------------------------------------
      ! computation of the horizontal dissipation terms
         if(selvar<5.5) then      ! horizontal diffusion can be for tracers too. In that case, selvar need not be < 5.5
            vardif=0.;
            if(biharmonic_horizontal) then
               call mixing_horizontal_biharmonic(var,vardif)
            else
               call mixing_horizontal(var,vardif)
            endif
      !      call mixing_isopycnal(var,vardif,10.);PRINT*,"VISCOUS REDI";
      !      call mixing_isopycnal_biharmonic(var,vardif,1.);PRINT*,"BIHARMONIC REDI"
#ifdef extra_output
            if(selvar==1) then; term_dissTh(1:NI,1:NJ,1:NK)=vardif(1:NI,1:NJ,1:NK)*Jacinv(1:NI,1:NJ,1:NK)/TL;    endif
            if(selvar==3) then; term_dissUh(1:NI,1:NJ,1:NK)=vardif(1:NI,1:NJ,1:NK)*Jacinv(1:NI,1:NJ,1:NK)/TL*UL; endif
            if(selvar==4) then; term_dissVh(1:NI,1:NJ,1:NK)=vardif(1:NI,1:NJ,1:NK)*Jacinv(1:NI,1:NJ,1:NK)/TL*UL; endif
            if(selvar==5) then; term_dissWh(1:NI,1:NJ,1:NK)=vardif(1:NI,1:NJ,1:NK)*Jacinv(1:NI,1:NJ,1:NK)/TL*WL; endif
#endif
            uvarx(1:NI,1:NJ,1:NK)=uvarx(1:NI,1:NJ,1:NK)-vardif(1:NI,1:NJ,1:NK)
         endif
      ! uvarx_advec(1:NI,1:NJ,1:NK) = uvarx(1:NI,1:NJ,1:NK) ! uvarx_advec is not needed. Will eliminate it later.

   ! ---------------------------------------------------------------
   !                    Solve equation of motion
   ! ---------------------------------------------------------------

#ifdef implicit
         if (selvar.lt.4.5) then
            do i=1,NI
            do j=1,NJ
               mat_A(i,j,1:NK,selvar)  = dtimel*Jacinv(i,j,1:NK)*mat_A(i,j,1:NK,selvar)
               mat_B(i,j,1:NK,selvar)  = ( (dtimel*Jacinv(i,j,1:NK)*mat_B(i,j,1:NK,selvar)) - 1)
               mat_C(i,j,1:NK-1,selvar)= dtimel*Jacinv(i,j,1:NK-1)*mat_C(i,j,1:NK-1,selvar)
               mat_D(i,j,1:NK,selvar)  = 0.d0 - var0(i,j,1:NK ) - (dtimel*Jacinv(i,j,1:NK)*mat_D(i,j,1:NK,selvar)) + &
               & (dtimel*Jacinv(i,j,1:NK)*uvarx_advec(i,j,1:NK))
               call solve_tridiag(mat_A(i,j,1:NK,selvar),mat_B(i,j,1:NK,selvar),mat_C(i,j,1:NK,selvar),mat_D(i,j,1:NK,selvar),mat_test(i,j,1:NK,selvar),NK)
            enddo
            enddo
         endif

         if(selvar==1) then; T(1:NI,1:NJ,1:NK,n)= mat_test(1:NI,1:NJ,1:NK,selvar); endif
         if(selvar==2) then; s(1:NI,1:NJ,1:NK,n)= mat_test(1:NI,1:NJ,1:NK,selvar); endif
         if(selvar==3) then; cx(1:NI,1:NJ,1:NK) = mat_test(1:NI,1:NJ,1:NK,selvar); endif
         if(selvar==4) then; cy(1:NI,1:NJ,1:NK) = mat_test(1:NI,1:NJ,1:NK,selvar); endif
         !if(selvar==5) then; cz(1:NI,1:NJ,1:NK) = mat_test(1:NI,1:NJ,1:NK,selvar); endif
#else
         ! -------------------  T  -------------------
         if(selvar==1) then
            if (restore_sT) then         
               do j=0,NJ+1
               do k=0,NK+1
                  Tbar(j,k)= 0.d0
                  do i=0,NI+1
                     Tbar(j,k)= Tbar(j,k) + T(i,j,k,0)
                  enddo
                  Tbar(j,k)= Tbar(j,k)/(dble(NI)+2)
               enddo
               enddo
            
               do k=1,NK
               do j=1,NJ
                  res_term= -restore_rate*TL*( Tbar(j,k)-T_restore(j,k) )
#ifdef extra_output
                  term_restT(j,k)= res_term/TL  
#endif
                  do i=1,NI
                     T(i,j,k,n)=T(i,j,k,0) - dtimel* (  Jacinv(i,j,k)*uvarx(i,j,k)-res_term )
                  enddo 
               enddo
               enddo
            else ! not restore T
               T(1:NI,1:NJ,1:NK,n)=T(1:NI,1:NJ,1:NK,0)-dtimel*Jacinv(1:NI,1:NJ,1:NK)*uvarx(1:NI,1:NJ,1:NK)
            endif
         endif ! end selvar==1

         ! -------------------  s  -------------------
         if(selvar==2) then
            write(6,*) "salinity calculated!"
            if (restore_sT) then         
               do j=0,NJ+1
               do k=0,NK+1
                  sbar(j,k)= 0.d0
                  do i=0,NI+1
                     sbar(j,k)= sbar(j,k) + s(i,j,k,0)
                  enddo
                  sbar(j,k)= sbar(j,k)/(dble(NI)+2)
               enddo
               enddo
            
               do k=1,NK
               do j=1,NJ
                  res_term= -restore_rate*TL*( sbar(j,k)-s_restore(j,k) )
                  do i=1,NI
                     s(i,j,k,n)=s(i,j,k,0) - dtimel* (  Jacinv(i,j,k)*uvarx(i,j,k)-res_term )
                  enddo 
               enddo
               enddo   
            else ! not restore s
               s(1:NI,1:NJ,1:NK,n)=s(1:NI,1:NJ,1:NK,0)-dtimel*Jacinv(1:NI,1:NJ,1:NK)*uvarx(1:NI,1:NJ,1:NK)
            endif       
         endif ! end selvar==2

         ! -------------------  u  -------------------
         if(selvar==3) then
            if (restore_u) then
               do j=0,NJ+1
               do k=0,NK+1
                  ubar(j,k)= 0.d0
                  do i=0,NI+1
                     ubar(j,k)= ubar(j,k) + u(i,j,k,0)
                  enddo
                  ubar(j,k)= ubar(j,k)/(dble(NI)+2)
               enddo
               enddo

               do k=1,NK
               do j=1,NJ
                  res_term= -restore_rate*TL*( ubar(j,k)-u_restore(j,k) )
#ifdef extra_output                  
                  term_restU(j,k)= res_term/TL*UL
#endif
                  do i=1,NI
                     cx(i,j,k)=u(i,j,k,0) - dtimel* (  Jacinv(i,j,k)*uvarx(i,j,k)-res_term )
                  enddo 
               enddo
               enddo
            else ! not restore u
               cx(1:NI,1:NJ,1:NK) =u(1:NI,1:NJ,1:NK,0)-dtimel*Jacinv(1:NI,1:NJ,1:NK)*uvarx(1:NI,1:NJ,1:NK)
            endif
         endif ! end selvar
   
         ! -------------------  v  -------------------
         if(selvar==4) then; cy(1:NI,1:NJ,1:NK) =v(1:NI,1:NJ,1:NK,0)-dtimel*Jacinv(1:NI,1:NJ,1:NK)*uvarx(1:NI,1:NJ,1:NK); endif
         ! -------------------  w  -------------------
         if(selvar==5) then; cz(1:NI,1:NJ,1:NK) =w(1:NI,1:NJ,1:NK,0)-dtimel*Jacinv(1:NI,1:NJ,1:NK)*uvarx(1:NI,1:NJ,1:NK); endif
         ! -------------------  tr  -------------------
         do it=1,ntr
            if(selvar==5+it) then; 
               Tr(it,1:NI,1:NJ,1:NK,n)=Tr(it,1:NI,1:NJ,1:NK,0)-dtimel*Jacinv(1:NI,1:NJ,1:NK)*uvarx(1:NI,1:NJ,1:NK); 
            endif
         enddo
#endif
        
      endif !(av_comp(selvar)==1)
   enddo !selvar

   return 
                                                                        
END