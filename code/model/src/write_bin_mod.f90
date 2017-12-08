subroutine save3d(Nx,Ny,Nz,var,filename)
   use header, only : rc_kind
   integer :: Nx,Ny,Nz
   character(len=*) :: filename
   REAL(kind=rc_kind) ::  var(Nx,Ny,Nz)
   open(3333,file=trim(filename),form='unformatted',access='stream',status='replace')
   write(3333) real(var,4)
   close(3333)
end subroutine 

subroutine save2d(Nx,Ny,var,filename)
   use header, only : rc_kind
   !Nx,Ny are not necessary in zonal and meridional directions.
   integer :: Nx,Ny
   character(len=*) :: filename
   REAL(kind=rc_kind) ::  var(Nx,Ny)

   open(3333,file=trim(filename),form='unformatted',access='stream',status='replace')
   write(3333) real(var,4)
   close(3333)
end subroutine 

subroutine w_pickup(filename)
   use header, only : h,u,v,w,uf,vf,wf,T,S,ntr,Tr,p,gradhn,hxn,hyn,rc_kind
   character(len=*) :: filename
   open(3333,file=trim(filename),form='unformatted',access='stream',status='replace')
   write(3333) h,u,v,w,uf,vf,wf,T,S,Tr,p,gradhn,hxn,hyn
   close(3333)
   open(3333,file=trim(filename)//'.meta',status='replace')
   write(3333,"(A1,4I3)") 'h',size(h,1),size(h,2),0,0
   write(3333,'(A1,4I3)') 'u',size(u,1),size(u,2),size(u,3),size(u,4)
   write(3333,"(A1,4I3)") 'v',size(v,1),size(v,2),size(v,3),size(v,4)
   write(3333,"(A1,4I3)") 'w',size(w,1),size(w,2),size(w,3),size(w,4)
   write(3333,"(A1,4I3)") 'uf',size(uf,1),size(uf,2),size(uf,3),0
   write(3333,"(A1,4I3)") 'vf',size(vf,1),size(vf,2),size(vf,3),0
   write(3333,"(A1,4I3)") 'wf',size(wf,1),size(wf,2),size(wf,3),0
   write(3333,"(A1,4I3)") 'T',size(T,1),size(T,2),size(T,3),size(T,4)
   write(3333,"(A1,4I3)") 'S',size(S,1),size(S,2),size(S,3),size(S,4)
   write(3333,"(A1,5I3)") 'Tr',ntr,size(S,1),size(S,2),size(S,3),size(S,4)
   write(3333,"(A1,3I3)") 'p',size(p,1),size(p,2),size(p,3)
   write(3333,"(A1,3I3)") 'gradhn',size(gradhn,1),size(gradhn,2),size(gradhn,3)
   write(3333,"(A1,3I3)") 'hxn',size(hxn,1),size(hxn,2),size(hxn,3)
   write(3333,"(A1,3I3)") 'hyn',size(hyn,1),size(hyn,2),size(hyn,3)

!   write(3333,100) 'pcorr',size(pcorr,1),size(pcorr,2),size(pcorr,3)
   close(3333)
end subroutine w_pickup

subroutine r_pickup(step)
   use header, only : h,u,v,w,uf,vf,wf,T,S,Tr,p,gradhn,hxn,hyn,rc_kind,dirout
   integer :: step
   character(len=10) :: stepchar


   open(3333,file=TRIM(dirout)//'op.pickup.'//stepchar(step)//'.bin',form='unformatted',access='stream',status='old')
   read(3333) h,u,v,w,uf,vf,wf,T,S,Tr,p,gradhn,hxn,hyn
   close(3333)
   print*, '# pickup at step '//stepchar(step)
end subroutine r_pickup 


subroutine load_dy(dyM)
   use header, only : NJ,rc_kind
   REAL(kind=rc_kind), intent(out) :: dyM(0:NJ+1)

   call io_error('dyM.data')

   open(3333,file='dyM.data',form='unformatted',access='stream',status='old')
   read(3333) dyM
   close(3333)
   print*, dyM
   print*, '# load dyM.data'
end subroutine load_dy

subroutine io_error(filename)
   use header, only : rc_kind

   character(len=*) filename
   logical :: lexist
   
   inquire(file=trim(filename),exist=lexist)

   if (.not. lexist) then
      print*, "# error: file "//trim(filename)//" does not exsit."
      stop
   endif
end subroutine io_error
