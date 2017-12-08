subroutine error(msg) 
!----------------------------------------------------                   
      USE header
      character :: msg(100)
      print*, msg,'================================================='
      print*, " max(w)",maxval(w)
      print*, " max(u)",maxval(u)
      print*, " max(v)",maxval(v)
      print*, " T(1,NJ,NJ+1,NK,:)",rho(1,NJ:NJ+1,NK)
      print*, " S(1,NJ,NJ+1,NK,:)",T(1,NJ:NJ+1,NK,:)
      print*, " rho(1,NJ,NJ+1,NK)",S(1,NJ:NJ+1,NK,:)
end subroutine error
      
