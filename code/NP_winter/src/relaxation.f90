#include "cppdefs.h"
#ifdef relaxationflag

MODULE relaxation
  USE header
  implicit none

CONTAINS
    ! This subroutine set the coefficient use to restore the boundary conditions and
    ! compensate for the wind curl. It is homogeneous in x, but can vary in y and z
  SUBROUTINE set_coef()
      implicit none
      integer :: i,j,k
      REAL(kind=rc_kind) :: gamma_T, scalefactor(NJ+1), windstressprofile(NJ+1)
      REAL(kind=rc_kind) :: ywindmin,ywindmax,edge
      logical :: exist

      ! specify the restoring timescale (i.e. 60 days)
      gamma_T = 1d0 / (6d0 * 86400d0 * UL/LEN)
        
      ! Creates a continous function throught the domain that scales
      ! the restoring with the wind curl (maximum when wind curl is maximum)

      ! This generates the wind profile used through the domain. It is copy-pasted
      ! from the windstress subroutine
      edge = 0.03    ! tightnness of the padding in the wind stress
      ywindmin = abs(atanh(2.d-2 - 1)/edge/PI) ! tapering region from southern bourndary
      ywindmax = NJ-ywindmin ! tapering region from northern bourndary
      do j=0,NJ+1
         if (j.lt.NJ/2) then
             windstressprofile(j)= 0.5*(tanh(edge*(j-ywindmin)*PI)+1.d0)
         else
             windstressprofile(j)= -0.5*(tanh(edge*(j-ywindmax)*PI)-1.d0)
         end if
      end do
      
      ! Generate the scalefactor, based on the gradient of the windstress profile
      do j = 0,NJ+1
          if (j.le.1) then
              scalefactor(j) = 0
          else if (j.ge.NJ) then
              scalefactor(j) = 0
          else
              scalefactor(j) = ABS(windstressprofile(j+1)-windstressprofile(j-1))
          end if
      enddo
      scalefactor = scalefactor/maxval(scalefactor)

        ! Print the scalefactor function into a file scalefactor.out
        inquire(file=TRIM(dirout)//'scalefactor.out', exist=exist)
        if (exist) then
        else
           open(12, file=TRIM(dirout)//'scalefactor.out', status="new", action="write")
           WRITE(12,*) "j, scalefactor,  windstressprofile"
           
           do j = 1,NJ
               WRITE(12,*) j, ",", scalefactor(j), ",", windstressprofile(j)
           enddo
           close(12)
        end if
        
        ! Compute the relaxation coeficient (product of relaxation timescale and a scalefactor,
        ! which depends on the wind curl)
        do k = 0, NK+1
            do j = 0,NJ+1
            r_T(j,k) = gamma_T * scalefactor(j)
            enddo
        enddo

  END SUBROUTINE set_coef

    ! This subroutine relaxes the temperature and salinity fields to a reference value
  SUBROUTINE sponge(step,n)   
        implicit none
        REAL(kind=rc_kind) :: Tvarbar(0:NJ+1,0:NK+1),svarbar(0:NJ+1,0:NK+1)
        REAL(KIND=rc_kind) :: Tnorth(NK),Tsouth(NK),snorth(NK),ssouth(NK),drhoml,rhocrit,zdepth
        REAL(KIND=rc_kind) :: Tref, sref
        integer :: i,j,k,n,step
        parameter (drhoml=0.01)
        logical :: exist
        character(len=5) :: x1
                
        ! Compute mld 
        do j=1,NJ
            do i=1,NI
                ! Find MLD mld
                rhocrit =rho(i,j,NK)+ drhoml
                do k=NK-1,1,-1
                    if (rho(i,j,k).ge.rhocrit) then
                        mld(i,j)= -DL*( (zc(i,j,k)-zc(i,j,k+1))*(rhocrit -rho(i,j,k+1))/(rho(i,j,k)-rho(i,j,k+1)) + zc(i,j,k+1))
                        goto 51
                    end if
                end do
                ! if the do loop is computed, and the ML base is not found, set it to D
                mld(i,j)= -zf(i,j,0)*DL
51              continue
            end do
        end do        
        
        ! Compute the reference to relax to
        do k=1,NK
            Tsouth(k) = 0.d0
            ssouth(k) = 0.d0
            Tnorth(k) = 0.d0
            snorth(k) = 0.d0
            do i=1,NI
                Tsouth(k)=Tsouth(k) + T(i,75,k,n)
                ssouth(k)=ssouth(k) + s(i,75,k,n)
                Tnorth(k)=Tnorth(k) + T(i,225,k,n)
                snorth(k)=snorth(k) + s(i,225,k,n)
            end do
            Tsouth(k)= Tsouth(k)/dble(NI)
            ssouth(k)= ssouth(k)/dble(NI)
            Tnorth(k)= Tnorth(k)/dble(NI)
            snorth(k)= snorth(k)/dble(NI)
        end do    

        ! calculate the bar{var} in x
        do k = 0,NK+1
            do j = 0,NJ+1
                Tvarbar(j,k) = sum(T(1:NI,j,k,n))/real(NI)
                svarbar(j,k) = sum(s(1:NI,j,k,n))/real(NI)
            enddo
        enddo
 
        do k = 0,NK+1
            do j = 0,NJ+1
                do i = 0,NI+1
                    ! Restore if below the mld
                     !zdepth=-zf(i,j,k)*DL
                     !if (zdepth.gt.mld(i,j)) then
                        ! Restore to Southern boundary
                        if (j.le.NJ/2) then
                            ! compute the reference temperature and salinity based on the background gradient
                            ! calculated between j = 75 and j = 225
                            !Tref = Tsouth(k) - (yc(75)-yc(j))*(Tnorth(k)-Tsouth(k))/(yc(225)-yc(75))
                            !sref = ssouth(k) - (yc(75)-yc(j))*(snorth(k)-ssouth(k))/(yc(225)-yc(75))
                            !T(i,j,k,n) = T(i,j,k,n) - dtf * r_T(j,k) * (Tvarbar(j,k) - Tref)
                            !s(i,j,k,n) = s(i,j,k,n) - dtf * r_T(j,k) * (svarbar(j,k) - sref)
                            T(i,j,k,n) = T(i,j,k,n) - dtf * r_T(j,k) * (Tvarbar(j,k) - Tsouth(k))
                            s(i,j,k,n) = s(i,j,k,n) - dtf * r_T(j,k) * (svarbar(j,k) - ssouth(k))
                        ! Restore to northern boundary
                        else if(j.ge.NJ/2) then
                            !Tref = Tsouth(k) - (yc(225)-yc(j))*(Tnorth(k)-Tsouth(k))/(yc(225)-yc(75))
                            !sref = ssouth(k) - (yc(225)-yc(j))*(snorth(k)-ssouth(k))/(yc(225)-yc(75))
                            !T(i,j,k,n) = T(i,j,k,n) - dtf * r_T(j,k) * (Tvarbar(j,k) - Tref)
                            !s(i,j,k,n) = s(i,j,k,n) - dtf * r_T(j,k) * (svarbar(j,k) - sref)
                            T(i,j,k,n) = T(i,j,k,n) - dtf * r_T(j,k) * (Tvarbar(j,k) - Tnorth(k))
                            s(i,j,k,n) = s(i,j,k,n) - dtf * r_T(j,k) * (svarbar(j,k) - snorth(k))
                        endif
                     !endif                    
                enddo
            enddo
        enddo
        
        ! Print the mld in mld.out
        if (mod(step,100).eq.0) then
            ! Get a string of the file number for name
            write (x1,'(I5.5)') step
            
            !inquire if the file exist
            inquire(file=TRIM(dirout)//'mld'//trim(x1)//'.out', exist=exist)
            if (exist) then
                ! If file exist, text is appended
                open(12, file=TRIM(dirout)//'mld'//trim(x1)//'.out', status="old", position="append", action="write")
            else
                ! If file doesnt exist, file is created
                open(12, file=TRIM(dirout)//'mld'//trim(x1)//'.out', status="new", action="write")
                WRITE(12,*) "i, ","j, ","MLD"
            end if
            do i = 0,NI+1
                do j = 0,NJ+1
                    WRITE(12,*) i, ",", j, ",", mld(i,j)
                end do
            end do
            close(12)
        end if

  END SUBROUTINE sponge

END MODULE relaxation

#endif
