subroutine linerelax(nxm,nym,nzm,cp,p,fn)
!     -----------------------------------------
  USE header, only : NI,NJ,NK,rc_kind
!     z-line Gauss-Seidel relaxation, called from the multigrid solver "mgrid":
!     solve  Grad p = fn  (19-point stencil cp), one tridiagonal solve in k per (i,j) column.
!
!     Multi-colour ordering: the stencil couples column (i,j) to the columns (i+-1,j), (i,j+-1)
!     and (i+-1,j+-1), so the columns are split into colours by the parity of i and of j.
!     Columns of one colour do not depend on each other and are relaxed in parallel (OpenMP);
!     the colours are done one after the other. The result is the same for any number of threads.
!     (The earlier lexicographic sweep, starting at i=NI/2 or a random i, was a data race with
!     OpenMP: threads read rows that other threads were updating.)
!     x is periodic: if nxm is odd, i=nxm (a neighbour of both i=1 and i=nxm-1) gets a third
!     i-colour.
      integer i,j,k,nxm,nym,nzm,info,im1,ip1,ci,cj,nci,ib,ie,is
      REAL(kind=rc_kind) :: p(0:nxm+1,0:nym+1,0:nzm+1),                      &
     &     fn(nxm,nym,nzm),cp(19,nxm,nym,nzm),subd(NK),supd(NK),        &
     &     dia(NK),rh(NK)

nci= 2
if (mod(nxm,2).ne.0) nci= 3
if (nxm.eq.1) nci= 1

do cj=1,2
do ci=1,nci
  if (ci.le.2) then
    ib= ci; ie= nxm; is= 2
    if (nci.eq.3) ie= nxm-1
  else
    ib= nxm; ie= nxm; is= 1
  endif

 !$OMP PARALLEL DO PRIVATE(i,j,im1,ip1,k,rh,subd,supd,dia,info) SHARED(p,cp,fn,nxm,nym,nzm,ib,ie,is,cj)
  do j=cj,nym,2
  do i=ib,ie,is
    im1= i-1
    ip1= i+1
    if (i.eq.1) im1= nxm
    if (i.eq.nxm) ip1= 1
    do k=1,nzm
       rh(k)= -( cp(2,i,j,k)*p(ip1,j,k)   &
       &   +cp(3,i,j,k)*p(i,j+1,k)        &
       &   +cp(4,i,j,k)*p(im1,j,k)        &
       &   +cp(5,i,j,k)*p(i,j-1,k)        &
       &   +cp(8,i,j,k)*p(im1,j+1,k)      &
       &   +cp(9,i,j,k)*p(im1,j-1,k)      &
       &   +cp(10,i,j,k)*p(ip1,j-1,k)     &
       &   +cp(11,i,j,k)*p(ip1,j+1,k)     &
       &   +cp(12,i,j,k)*p(im1,j,k-1)     &
       &   +cp(13,i,j,k)*p(ip1,j,k-1)     &
       &   +cp(14,i,j,k)*p(ip1,j,k+1)     &
       &   +cp(15,i,j,k)*p(im1,j,k+1)     &
       &   +cp(16,i,j,k)*p(i,j-1,k-1)     &
       &   +cp(17,i,j,k)*p(i,j+1,k-1)     &
       &   +cp(18,i,j,k)*p(i,j+1,k+1)     &
       &   +cp(19,i,j,k)*p(i,j-1,k+1)     &
       &   - fn(i,j,k) )
       if (k.ne.1) subd(k)= cp(7,i,j,k)
       if (k.ne.nzm) supd(k)= cp(6,i,j,k)
       dia(k)= cp(1,i,j,k)
     enddo
     rh(1)= rh(1) -cp(7,i,j,1)*p(i,j,0)
     rh(nzm)= rh(nzm) -cp(6,i,j,nzm)*p(i,j,nzm+1)
     call dgtsl(nzm,subd,dia,supd,rh,info)
     if (info.ne.0) then
       write(6,*) 'error in linerelax'
       stop
     endif
     p(i,j,1:nzm)= rh(1:nzm)
  enddo
  enddo
 !$OMP END PARALLEL DO

enddo
enddo

return
END
