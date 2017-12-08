subroutine mixing_isopycnal(var,vardif,rv_kr) 
 ! ---------------------------------------------                     
 ! diffusion along isopycnals.
 ! 06 Feb 2013. JBG.
 !
 USE header
 !                                                                       
 IMPLICIT NONE 

 INTEGER :: i,j,k
 INTEGER :: ip1,im1
 REAL :: rv_kr

 REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1, 0:NK+1) :: var
 REAL(kind=rc_kind) :: vardif(NI,NJ,NK)

 REAL(kind=rc_kind) :: sxfx(0:NI), sxfy(0:NJ),sxfz(0:NK)
 REAL(kind=rc_kind) :: syfx(0:NI), syfy(0:NJ),syfz(0:NK)

 REAL(kind=rc_kind) :: drdxfx(0:NI), drdxfy(0:NJ),drdxfz(0:NK)
 REAL(kind=rc_kind) :: drdyfx(0:NI), drdyfy(0:NJ),drdyfz(0:NK)
 REAL(kind=rc_kind) :: drdzfx(0:NI), drdzfy(0:NJ),drdzfz(0:NK)

 REAL(kind=rc_kind) :: dvdxfx(0:NI), dvdxfy(0:NJ),dvdxfz(0:NK)
 REAL(kind=rc_kind) :: dvdyfx(0:NI), dvdyfy(0:NJ),dvdyfz(0:NK)
 REAL(kind=rc_kind) :: dvdzfx(0:NI), dvdzfy(0:NJ),dvdzfz(0:NK)

 REAL(kind=rc_kind) :: flfx(0:NI), flfy(0:NJ),flfz(0:NK)

 REAL(kind=rc_kind), dimension(0:NI+1,0:NJ+1, 0:NK+1) :: kredi
 REAL(kind=rc_kind) :: fac

 LOGICAL :: lv_smsl=.FALSE.

 kredi(:,:,:)=rv_kr

 fac= 1.0/(UL*LEN)
                                                                        

 !*********************
 ! FLUX IN X-DIRECTION
 !---------------------

 DO k=1,NK    
   DO j=1,NJ

     DO i=1,NI ! For all faces_x (face 0 will be equal to face NI),
       IF(i==NI) then ! For the case of face_x(NI), the neighbour on the right is the cell 1.
         ip1=1
        ELSE
         ip1=i+1
       ENDIF

       !----------------
       ! VARIATIONS IN X

                !    delta_var               /       delta_x 
       drdxfx(i)=( rho(ip1,j,k)-rho(i,j,k) ) * ( 0.5d0*(ux(i,j)+ux(ip1,j)) )
       dvdxfx(i)=( var(ip1,j,k)-var(i,j,k) ) * ( 0.5d0*(ux(i,j)+ux(ip1,j)) )


       !----------------
       ! VARIATIONS IN Y

       if(j==1) then
                  ! [          rho_face_x(i,j+1,k)     -        rho_face_x(i,j,k)       ] / [    vy estimated between faces_x (i,j) and (i,j+1) ]   
         drdyfx(i)=( 0.5*(rho(i,j+1,k)+rho(ip1,j+1,k)) - 0.5*(rho(i,1,k)+rho(ip1,1,k) ) ) * ( ( 0.25d0* (vy(i,j)+vy(ip1,j)+vy(i,j+1)+vy(ip1,j+1)) ))  
         dvdyfx(i)=( 0.5*(var(i,j+1,k)+var(ip1,j+1,k)) - 0.5*(var(i,1,k)+var(ip1,1,k) ) ) * ( ( 0.25d0* (vy(i,j)+vy(ip1,j)+vy(i,j+1)+vy(ip1,j+1)) ))  
        elseif(j==NJ) then 
                  ! [          rho_face_x(i,j,k)     -          rho_face_x(i,j-1,k)       ] / [    vy estimated between faces_x (i,j) and (i,j-1) ]   
         drdyfx(i)=( 0.5*(rho(i,NJ,k)+rho(ip1,NJ,k)) - 0.5*(rho(i,j-1,k)+rho(ip1,j-1,k) ) ) * ( ( 0.25d0* (vy(i,NJ)+vy(ip1,NJ)+vy(i,j-1)+vy(ip1,j-1))) )  
         dvdyfx(i)=( 0.5*(var(i,NJ,k)+var(ip1,NJ,k)) - 0.5*(var(i,j-1,k)+var(ip1,j-1,k) ) ) * ( ( 0.25d0* (vy(i,NJ)+vy(ip1,NJ)+vy(i,j-1)+vy(ip1,j-1))) )  
        else
                  ! [        rho_face_x(i,j+1,k)       -        rho_face_x(i,j-1,k)         ] / [    vy estimated between faces_x (i,j) and (i,j+1)   +  vy estimated between faces_x (i,j) and (i,j-1) ]   
         drdyfx(i)=( 0.5*(rho(i,j+1,k)+rho(ip1,j+1,k)) - 0.5*(rho(i,j-1,k)+rho(ip1,j-1,k) ) ) * 0.25d0 * ( ( 0.25d0* (vy(i,j)+vy(ip1,j)+vy(i,j+1)+vy(ip1,j+1)) + 0.25d0*(vy(i,j)+vy(ip1,j)+vy(i,j-1)+vy(ip1,j-1)) ) )  
         dvdyfx(i)=( 0.5*(var(i,j+1,k)+var(ip1,j+1,k)) - 0.5*(var(i,j-1,k)+var(ip1,j-1,k) ) ) * 0.25d0 * ( ( 0.25d0* (vy(i,j)+vy(ip1,j)+vy(i,j+1)+vy(ip1,j+1)) + 0.25d0*(vy(i,j)+vy(ip1,j)+vy(i,j-1)+vy(ip1,j-1)) ) )  
       endif


       !----------------
       ! VARIATIONS IN Z

       if(k==1) then
                  ! [         rho_face_x(i,j,k+1)      -         rho_face_x(i,j,k)      ] / [    wz estimated between faces_x (i,j,k) and (i,j,k+1)   ]   
         drdzfx(i)=( 0.5*(rho(i,j,k+1)+rho(ip1,j,k+1)) - 0.5*(rho(i,j,k)+rho(ip1,j,k) ) ) * ( ( 0.25d0* (wz(i,j,k)+wz(ip1,j,k)+wz(i,j,k+1)+wz(ip1,j,k+1))  ) )  
         dvdzfx(i)=( 0.5*(var(i,j,k+1)+var(ip1,j,k+1)) - 0.5*(var(i,j,k)+var(ip1,j,k) ) ) * ( ( 0.25d0* (wz(i,j,k)+wz(ip1,j,k)+wz(i,j,k+1)+wz(ip1,j,k+1))  ) )  

        elseif(k==NK) then
                ! [          rho_face_x(i,j,k)     -          rho_face_x(i,j,k-1)     ] / [    wz estimated between faces_x (i,j,k) and (i,j,k-1) ]   
         drdzfx(i)=( 0.5*(rho(i,j,NK)+rho(ip1,j,NK)) - 0.5*(rho(i,j,k-1)+rho(ip1,j,k-1) ) ) * ( ( 0.25d0*(wz(i,j,NK)+wz(ip1,j,NK)+wz(i,j,k-1)+wz(ip1,j,k-1)) ) )  
         dvdzfx(i)=( 0.5*(var(i,j,NK)+var(ip1,j,NK)) - 0.5*(var(i,j,k-1)+var(ip1,j,k-1) ) ) * ( ( 0.25d0*(wz(i,j,NK)+wz(ip1,j,NK)+wz(i,j,k-1)+wz(ip1,j,k-1)) ) )  

        else
                ! [          rho_face_x(i,j,k+1)   -          rho_face_x(i,j,k-1)     ] / [    wz estimated between faces_x (i,j,k) and (i,j,k+1)   +  wz estimated between faces_x (i,j,k) and (i,j,k-1) ]   
         drdzfx(i)=( 0.5*(rho(i,j,k+1)+rho(ip1,j,k+1)) - 0.5*(rho(i,j,k-1)+rho(ip1,j,k-1) ) ) * 0.25d0 * ( ( 0.25d0* (wz(i,j,k)+wz(ip1,j,k)+wz(i,j,k+1)+wz(ip1,j,k+1)) + 0.25d0*(wz(i,j,k)+wz(ip1,j,k)+wz(i,j,k-1)+wz(ip1,j,k-1)) ) )  
         dvdzfx(i)=( 0.5*(var(i,j,k+1)+var(ip1,j,k+1)) - 0.5*(var(i,j,k-1)+var(ip1,j,k-1) ) ) * 0.25d0 * ( ( 0.25d0* (wz(i,j,k)+wz(ip1,j,k)+wz(i,j,k+1)+wz(ip1,j,k+1)) + 0.25d0*(wz(i,j,k)+wz(ip1,j,k)+wz(i,j,k-1)+wz(ip1,j,k-1)) ) )  
       endif

       !------------------------
       ! COMPUTATION OF THE FLUX
     
       if(lv_smsl) then
         flfx(i)=kredi(i,j,k)/(                          drdzfx(i)**2)* ( (             drdzfx(i)**2)*dvdxfx(i) -                   0*dvdyfx(i) - drdxfx(i)*drdzfx(i)*dvdzfx(i) ) 
        else
         flfx(i)=kredi(i,j,k)/(drdxfx(i)**2+drdyfx(i)**2+drdzfx(i)**2)* ( (drdyfx(i)**2+drdzfx(i)**2)*dvdxfx(i) - drdxfx(i)*drdyfx(i)*dvdyfx(i) - drdxfx(i)*drdzfx(i)*dvdzfx(i) ) 
       endif

     ENDDO !_ i

     flfx(0)=flfx(NI)

     !-------------------------------
     ! COMPUTATION OF THE SOURCE TERM

     do i=1,NI
       vardif(i,j,k)=    0   +    fac*Jac(i,j,k)*ux(i,j)* ( flfx(i) - flfx(i-1) )


     ! IF(i==48 .AND. k==18 .AND. j==131) then
     !   PRINT*,vardif(i,j,k),flfx(i),flfx(i-1)
     ! ENDIF 
    enddo

   ENDDO !_ j
  ENDDO !_ k



 !*********************
 ! FLUX IN Y-DIRECTION
 !---------------------

 DO i=1,NI
   IF(i==1) then
     im1=NI;ip1=i+1
    elseif(i==NI) then
     im1=i-1;ip1=1
    else
     im1=i-1;ip1=i+1
   ENDIF
   DO k=1,NK    

     DO j=1,NJ-1 ! For all faces_y


       !----------------
       ! VARIATIONS IN X

                ! [          rho_face_y(i+1,j,k)   -          rho_face_y(i-1,j,k)     ] / [    ux estimated between faces_y (i,j) and (i+1,j)   +  ux estimated between faces_y (i,j) and (i-1,j) ]   
       drdxfy(j)=( 0.5*(rho(ip1,j,k)+rho(ip1,j+1,k)) - 0.5*(rho(im1,j,k)+rho(im1,j+1,k) ) ) * 0.25d0 * ( ( 0.25d0* (ux(i,j)+ux(i,j+1)+ux(ip1,j)+ux(ip1,j+1)) + 0.25d0*(ux(i,j)+ux(i,j+1)+ux(im1,j)+ux(im1,j+1)) ) )  
       dvdxfy(j)=( 0.5*(var(ip1,j,k)+var(ip1,j+1,k)) - 0.5*(var(im1,j,k)+var(im1,j+1,k) ) ) * 0.25d0 * ( ( 0.25d0* (ux(i,j)+ux(i,j+1)+ux(ip1,j)+ux(ip1,j+1)) + 0.25d0*(ux(i,j)+ux(i,j+1)+ux(im1,j)+ux(im1,j+1)) ) )  


       !----------------
       ! VARIATIONS IN Y

                !    delta_var               /       delta_y 
       drdyfy(j)=( rho(i,j+1,k)-rho(i,j,k) ) * ( 0.5d0*(vy(i,j)+vy(i,j+1)) )
       dvdyfy(j)=( var(i,j+1,k)-var(i,j,k) ) * ( 0.5d0*(vy(i,j)+vy(i,j+1)) )


       !----------------
       ! VARIATIONS IN Z


       IF(k==1) then
                  ! [          rho_face_y(i,j,k+1)     -        rho_face_y(i,j,k-1)     ] / [    wz estimated between faces_x (i,j,k) and (i,j,k+1)   +  vy estimated between faces_x (i,j,k) and (i,j,k-1) ]   
         drdzfy(j)=( 0.5*(rho(i,j,k+1)+rho(i,j+1,k+1)) - 0.5*(rho(i,j,k)+rho(i,j+1,k) ) ) * ( ( 0.25d0* (wz(i,j,k)+wz(i,j+1,k)+wz(i,j,k+1)+wz(i,j+1,k+1)) ))  
         dvdzfy(j)=( 0.5*(var(i,j,k+1)+var(i,j+1,k+1)) - 0.5*(var(i,j,k)+var(i,j+1,k) ) ) * ( ( 0.25d0* (wz(i,j,k)+wz(i,j+1,k)+wz(i,j,k+1)+wz(i,j+1,k+1)) ))  

        elseif(k==NK) then

                  ! [        rho_face_y(i,j,k+1)   -           rho_face_y(i,j,k-1)      ] / [    wz estimated between faces_x (i,j,k) and (i,j,k+1)   +  vy estimated between faces_x (i,j,k) and (i,j,k-1) ]   
         drdzfy(j)=( 0.5*(rho(i,j,k)+rho(i,j+1,k)) - 0.5*(rho(i,j,k-1)+rho(i,j+1,k-1) ) ) * ( ( 0.25d0* (wz(i,j,k)+wz(i,j+1,k)+wz(i,j,k-1)+wz(i,j+1,k-1)) ) )  
         dvdzfy(j)=( 0.5*(var(i,j,k)+var(i,j+1,k)) - 0.5*(var(i,j,k-1)+var(i,j+1,k-1) ) ) * ( ( 0.25d0* (wz(i,j,k)+wz(i,j+1,k)+wz(i,j,k-1)+wz(i,j+1,k-1)) ) )  

        else
                    ! [          rho_face_y(i,j,k+1)   -          rho_face_y(i,j,k-1)     ] / [    wz estimated between faces_x (i,j,k) and (i,j,k+1)   +  vy estimated between faces_x (i,j,k) and (i,j,k-1) ]   
         drdzfy(j)=( 0.5*(rho(i,j,k+1)+rho(i,j+1,k+1)) - 0.5*(rho(i,j,k-1)+rho(i,j+1,k-1) ) ) * 0.25d0 * ( ( 0.25d0* (wz(i,j,k)+wz(i,j+1,k)+wz(i,j,k+1)+wz(i,j+1,k+1)) + 0.25d0*(wz(i,j,k)+wz(i,j+1,k)+wz(i,j,k-1)+wz(i,j+1,k-1)) ) )  
         dvdzfy(j)=( 0.5*(var(i,j,k+1)+var(i,j+1,k+1)) - 0.5*(var(i,j,k-1)+var(i,j+1,k-1) ) ) * 0.25d0 * ( ( 0.25d0* (wz(i,j,k)+wz(i,j+1,k)+wz(i,j,k+1)+wz(i,j+1,k+1)) + 0.25d0*(wz(i,j,k)+wz(i,j+1,k)+wz(i,j,k-1)+wz(i,j+1,k-1)) ) )  

      ENDIF

      !------------------------
      ! COMPUTATION OF THE FLUX

      if(lv_smsl) then
        flfy(j)=kredi(i,j,k)/(                          drdzfy(j)**2)* ( -                   0*dvdxfy(j) + (             drdzfy(j)**2)*dvdyfy(j)  - drdyfy(j)*drdzfy(j)*dvdzfy(j) )
       else
        flfy(j)=kredi(i,j,k)/(drdxfy(j)**2+drdyfy(j)**2+drdzfy(j)**2)* ( - drdxfy(j)*drdyfy(j)*dvdxfy(j) + (drdxfy(j)**2+drdzfy(j)**2)*dvdyfy(j)  - drdyfy(j)*drdzfy(j)*dvdzfy(j) )
      endif

    ENDDO !_ j

    flfy(0)=0.;
    flfy(NJ)=0.;


    !-------------------------------
    ! COMPUTATION OF THE SOURCE TERM

    do j=1,NJ
      vardif(i,j,k)= vardif(i,j,k) +    fac*Jac(i,j,k)*vy(i,j)* ( flfy(j) - flfy(j-1) )

      !IF(i==48 .AND. k==18 .AND. j==131) then
      !  PRINT*,vardif(i,j,k),flfy(j),flfy(j-1)
      !ENDIF
    enddo

  ENDDO !_ k
 ENDDO !_ i


 !*********************
 ! FLUX IN Z-DIRECTION
 !---------------------

  DO i=1,NI

    IF(i==1) then
      im1=NI;ip1=i+1
     elseif(i==NI) then
      im1=i-1;ip1=1
     else
      im1=i-1;ip1=i+1
    ENDIF

   DO j=1,NJ    

     DO k=1,NK-1 ! For all faces_k

       !----------------
       ! VARIATIONS IN X

                ! [          rho_face_z(i+1,j,k)   -          rho_face_z(i-1,j,k)     ] / [    ux estimated between faces_z (i,j,k) and (i+1,j,k)   +  ux estimated between faces_z (i,j,k) and (i-1,j,k) ]   
       drdxfz(k)=( 0.5*(rho(ip1,j,k)+rho(ip1,j,k+1)) - 0.5*(rho(im1,j,k)+rho(im1,j,k+1) ) ) * 0.25d0 * ( ( 0.25d0* (ux(i,j)+ux(i,j)+ux(ip1,j)+ux(ip1,j)) + 0.25d0*(ux(i,j)+ux(i,j)+ux(im1,j)+ux(im1,j)) ) )  
       dvdxfz(k)=( 0.5*(var(ip1,j,k)+var(ip1,j,k+1)) - 0.5*(var(im1,j,k)+var(im1,j,k+1) ) ) * 0.25d0 * ( ( 0.25d0* (ux(i,j)+ux(i,j)+ux(ip1,j)+ux(ip1,j)) + 0.25d0*(ux(i,j)+ux(i,j)+ux(im1,j)+ux(im1,j)) ) )  

       !----------------
       ! VARIATIONS IN Y

       IF(j==1) then

                  ! [          rho_face_z(i,j+1,k)   -          rho_face_z(i,j-1,k)     ] / [    vy estimated between faces_z (i,j,k) and (i,j+1,k)   + vy estimated between faces_z (i,j,k) and (i,j-1,k) ]   
         drdyfz(k)=( 0.5*(rho(i,j+1,k)+rho(i,j+1,k+1)) - 0.5*(rho(i,j,k)+rho(i,j,k+1) ) ) * ( ( 0.25d0* (vy(i,j)+vy(i,j)+vy(i,j+1)+vy(i,j+1)) ) )   
         dvdyfz(k)=( 0.5*(var(i,j+1,k)+var(i,j+1,k+1)) - 0.5*(var(i,j,k)+var(i,j,k+1) ) ) * ( ( 0.25d0* (vy(i,j)+vy(i,j)+vy(i,j+1)+vy(i,j+1)) ) )  
      
        elseif(j==NJ) then
      
                  ! [          rho_face_z(i,j+1,k)   -          rho_face_z(i,j-1,k)     ] / [    vy estimated between faces_z (i,j,k) and (i,j+1,k)   + vy estimated between faces_z (i,j,k) and (i,j-1,k) ]   
         drdyfz(k)=( 0.5*(rho(i,j,k)+rho(i,j,k+1)) - 0.5*(rho(i,j-1,k)+rho(i,j-1,k+1) ) ) * ( ( 0.25d0*(vy(i,j)+vy(i,j)+vy(i,j-1)+vy(i,j-1)) ) )  
         dvdyfz(k)=( 0.5*(var(i,j,k)+var(i,j,k+1)) - 0.5*(var(i,j-1,k)+var(i,j-1,k+1) ) ) * ( ( 0.25d0*(vy(i,j)+vy(i,j)+vy(i,j-1)+vy(i,j-1)) ) )  
      
        else
      
                  ! [          rho_face_z(i,j+1,k)   -          rho_face_z(i,j-1,k)     ] / [    vy estimated between faces_z (i,j,k) and (i,j+1,k)   + vy estimated between faces_z (i,j,k) and (i,j-1,k) ]   
         drdyfz(k)=( 0.5*(rho(i,j+1,k)+rho(i,j+1,k+1)) - 0.5*(rho(i,j-1,k)+rho(i,j-1,k+1) ) ) * 0.25d0 * ( ( 0.25d0* (vy(i,j)+vy(i,j)+vy(i,j+1)+vy(i,j+1)) + 0.25d0*(vy(i,j)+vy(i,j)+vy(i,j-1)+vy(i,j-1)) ) )  
         dvdyfz(k)=( 0.5*(var(i,j+1,k)+var(i,j+1,k+1)) - 0.5*(var(i,j-1,k)+var(i,j-1,k+1) ) ) * 0.25d0 * ( ( 0.25d0* (vy(i,j)+vy(i,j)+vy(i,j+1)+vy(i,j+1)) + 0.25d0*(vy(i,j)+vy(i,j)+vy(i,j-1)+vy(i,j-1)) ) )  

       ENDIF

       !----------------
       ! VARIATIONS IN Z

                !    delta_var               /       delta_z 
       drdzfz(k)=( rho(i,j,k+1)-rho(i,j,k) ) * ( 0.5d0*(wz(i,j,k)+wz(i,j,k+1)) )
       dvdzfz(k)=( var(i,j,k+1)-var(i,j,k) ) * ( 0.5d0*(wz(i,j,k)+wz(i,j,k+1)) )


       !------------------------
       ! COMPUTATION OF THE FLUX
   
       if(lv_smsl) then
         flfz(k)=kredi(i,j,k)/(                          drdzfz(k)**2)* ( - drdxfz(k)*drdyfz(k)*dvdxfz(k)  - drdyfz(k)*drdzfz(k)*dvdyfz(k) + (drdxfz(k)**2+drdyfz(k)**2)*dvdzfz(k) )
        else
         flfz(k)=kredi(i,j,k)/(drdxfz(k)**2+drdyfz(k)**2+drdzfz(k)**2)* ( - drdxfz(k)*drdyfz(k)*dvdxfz(k)  - drdyfz(k)*drdzfz(k)*dvdyfz(k) + (drdxfz(k)**2+drdyfz(k)**2)*dvdzfz(k) )
       endif

     ENDDO !_ k

     flfz(0)=0.;
     flfz(NK)=0.;

    !-------------------------------
    ! COMPUTATION OF THE SOURCE TERM

     do k=1,NK
       vardif(i,j,k)= vardif(i,j,k) +    fac*Jac(i,j,k)*wz(i,j,k)* ( flfz(k) - flfz(k-1) )
 
    !  IF(i==48 .AND. k>16.5 .AND. k<18.5 .AND. j==131) then

    !    PRINT*,k,vardif(i,j,k),flfz(k),flfz(k-1)
    !    PRINT*,kredi(i,j,k)/(drdxfz(k)**2+drdyfz(k)**2+drdzfz(k)**2)* ( - drdxfz(k)*drdyfz(k)*dvdxfz(k)                                                                          )
    !    PRINT*,kredi(i,j,k)/(drdxfz(k)**2+drdyfz(k)**2+drdzfz(k)**2)* (                                  - drdyfz(k)*drdzfz(k)*dvdyfz(k)                                         )
    !    PRINT*,kredi(i,j,k)/(drdxfz(k)**2+drdyfz(k)**2+drdzfz(k)**2)* (                                                                  + (drdxfz(k)**2+drdyfz(k)**2)*dvdzfz(k) )
    !    PRINT*,(drdxfz(k)**2+drdyfz(k)**2)
    !    PRINT*,dvdzfz(k)
    !    PRINT*,var(i,j+1,k),var(i,j+1,k+1),var(i,j-1,k),var(i,j-1,k+1) 
    !    PRINT*,k,"-------"
    !  ENDIF
    enddo

   ENDDO !_ j
 ENDDO !_ i

                                                                       
 RETURN 
END subroutine
