      subroutine tracerinit(stepl)
!     --------------------
!     initializes tracer fields 
!     
      USE header
      integer  i,j,k,l,n,stepl
      double precision phydepth,z,trseed
!     
!     Tracer 1 is in mg carbon / m^3

!     if this is a continuation of a run, only trinit is to be initialized
!     Tr 1 is PHY

      phydepth= 30.d0
      trseed=1.2    ! Chl =0.08 mg/l, Chl2C=15

!     TRACER 1
!     =========

!     Initialze PHY
      do j=0,NJ+1
         do i=0,NI+1
            do k=NK,1,-1
               z = -zc(i,j,k)*DL
               if (z.lt.mldepth) Tr(1,i,j,k,0)=trseed  
!               T(1,i,j,k,0)= 10.d0*exp(-z/phydepth)
            end do
         end do
      end do

!     Periodic Boundaries (periodic e-w)
      do n=0,1
         do k=0,NK+1
            do j=0,NJ+1
               do it=1,ntr
                  Tr(it,0,j,k,n)= Tr(it,NI,j,k,n)
                  Tr(it,NI+1,j,k,n)= Tr(it,1,j,k,n)
               end do
            end do
         end do
      end do
      return
      end

