subroutine hsolve(h,oldh,hdt,dtime) 
!     --------------                                                    
  USE header, only : NI,NJ,rc_kind
!     modified for periodicew boundaries                                
!     ----------------------------------                                
      integer iter,i,j,kount,l 
      REAL(kind=rc_kind) :: ch(9,NI,NJ),rhs(NI,NJ),                          &
     &     oldh(0:NI+1,0:NJ+1),h(0:NI+1,0:NJ+1),res,                    &
     &     hdt(0:NI+1,0:NJ+1),dti,dtime,tol,maxres,hstar,rlx            
!     &     dti,dtime,tol,maxres,hstar,rlx                              
                                                                        
  tol= 1.d-12 
  dti= 1.d0/dtime 
  kount=0 

!     compute fine grid coefficients and rhs on fine grid               
  call chfine(dtime,ch,rhs) 

  rlx= 1.72 
  do iter=1,100 

    do j=1,NJ 

      i=1 
      hstar= ( ch(2,i,j)*h(i+1,j)      &
          &   +ch(3,i,j)*h(NI,j)       &
          &   +ch(4,i,j)*h(i,j+1)      &
          &   +ch(5,i,j)*h(i,j-1)      &
          &   +ch(6,i,j)*h(i+1,j+1)    &
          &   +ch(7,i,j)*h(NI,j+1)     &
          &   +ch(8,i,j)*h(i+1,j-1)    &
          &   +ch(9,i,j)*h(NI,j-1)     &
          &   - rhs(i,j) )/(-ch(1,i,j))    
      h(i,j)= (1.d0 -rlx)*h(i,j) +rlx*hstar 
      do i=2,NI-1 
        hstar= ( ch(2,i,j)*h(i+1,j)       &
            &   +ch(3,i,j)*h(i-1,j)       &
            &   +ch(4,i,j)*h(i,j+1)       &
            &   +ch(5,i,j)*h(i,j-1)       &
            &   +ch(6,i,j)*h(i+1,j+1)     &
            &   +ch(7,i,j)*h(i-1,j+1)     &
            &   +ch(8,i,j)*h(i+1,j-1)     &
            &   +ch(9,i,j)*h(i-1,j-1)     &
            &   - rhs(i,j) )/(-ch(1,i,j)) 
        h(i,j)= (1.d0 -rlx)*h(i,j) +rlx*hstar 
      enddo ! i 
      i=NI 
      hstar= ( ch(2,i,j)*h(1,j)         &
          &   +ch(3,i,j)*h(i-1,j)       &
          &   +ch(4,i,j)*h(i,j+1)       &
          &   +ch(5,i,j)*h(i,j-1)       &
          &   +ch(6,i,j)*h(1,j+1)       &
          &   +ch(7,i,j)*h(i-1,j+1)     &
          &   +ch(8,i,j)*h(1,j-1)       &
          &   +ch(9,i,j)*h(i-1,j-1)     &
          &   - rhs(i,j) )/(-ch(1,i,j)) 
      h(i,j)= (1.d0 -rlx)*h(i,j) +rlx*hstar 

    enddo ! j

    do l=1,3 
      call hfill(dtime,h) 
    enddo
                                                                          
    maxres= 0.d0 

    do j=1,NJ 
      i=1 
      res= rhs(i,j)-              &
        & (ch(1,i,j)*h(i,j)       &
        & +ch(2,i,j)*h(i+1,j)     &
        & +ch(3,i,j)*h(NI,j)      &
        & +ch(4,i,j)*h(i,j+1)     &
        & +ch(5,i,j)*h(i,j-1)     &
        & +ch(6,i,j)*h(i+1,j+1)   &
        & +ch(7,i,j)*h(NI,j+1)    &
        & +ch(8,i,j)*h(i+1,j-1)   &
        & +ch(9,i,j)*h(NI,j-1))   
      if (abs(res).gt.maxres) then 
        maxres= abs(res) 
      end if 
      do i=2,NI-1    !   res = b - A x'                                                    
        res= rhs(i,j)-              &
          & (ch(1,i,j)*h(i,j)       &
          & +ch(2,i,j)*h(i+1,j)     &
          & +ch(3,i,j)*h(i-1,j)     &
          & +ch(4,i,j)*h(i,j+1)     &
          & +ch(5,i,j)*h(i,j-1)     &
          & +ch(6,i,j)*h(i+1,j+1)   &
          & +ch(7,i,j)*h(i-1,j+1)   &
          & +ch(8,i,j)*h(i+1,j-1)   &
          & +ch(9,i,j)*h(i-1,j-1))  
        if (abs(res).gt.maxres) then 
          maxres= abs(res) 
        end if 
      enddo !i
      i=NI 
      res= rhs(i,j)-              &
        & (ch(1,i,j)*h(i,j)       &
        & +ch(2,i,j)*h(1,j)       &
        & +ch(3,i,j)*h(i-1,j)     &
        & +ch(4,i,j)*h(i,j+1)     &
        & +ch(5,i,j)*h(i,j-1)     &
        & +ch(6,i,j)*h(1,j+1)     &
        & +ch(7,i,j)*h(i-1,j+1)   &
        & +ch(8,i,j)*h(1,j-1)     &
        & +ch(9,i,j)*h(i-1,j-1))  
      if (abs(res).gt.maxres) then 
         maxres= abs(res) 
      end if 
    enddo !j

    !  write(6,*) 'iter,maxres',iter,maxres                          
    if (maxres.gt.3000.d0) then 
      print*, 'hsolve.f90: STOP. res too large, i,j,maxres=',i,j,maxres                                             
      stop 
    end if 
    if (maxres.lt.tol) goto 13 
  enddo !iter 
 
13 continue
                                                                       
  do l=1,1 
     call mprove(h,ch,rhs,dtime) 
  end do 
                                                                        
return 

END                                           
