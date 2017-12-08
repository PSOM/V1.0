subroutine tracerinit(stepl) 
  USE header, only : Tr
  integer :: it,stepl
!     initializes tracer fields                                         
!     TRACER 1                                                          
!     =========                                                         
      it=1 
      Tr = 0d0
      Tr(it,:,:,1,0) =1d0
end subroutine tracerinit
