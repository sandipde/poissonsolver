subroutine check_result(m)
use kind_param
use global_var
use poisson_var
implicit none
integer(I4B)::i,j,k,ix,iy,iz,i1,no,m
real(R8B)::factor
allocate (pot_diff(ndim_x*ndim_y*ndim_z))
 pot_diff=0.0
 min_diff=0
 max_diff=0
no=0
 if (choice==1)then
! kappa  =  a0/2.0
   do i=1,ndim_x
     do j=1,ndim_y
       do k=1,ndim_z
         no=no+1
         r=dsqrt(xi(i)*xi(i)+yi(j)*yi(j)+zi(k)*zi(k))
         if (r==0.0) then
           exact_pot(i,j,k)= 2.0d0/dsqrt(pi) * kappa
         else
         exact_pot(i,j,k)=derf(kappa*r)/r
         endif
        
        pot_diff(no)=abs(exact_pot(i,j,k)-hartree_pot(i,j,k))
       enddo
     enddo
  enddo
write(6,'(F6.2,6X,E10.3,6X,E10.3)') kappa,maxval(pot_diff),minval(pot_diff)
!deallocate (pot_diff)
! call del2(ndim_x,ndim_y,ndim_z,delta_x,hartree_pot,exact_pot)
! exact_pot=-exact_pot/4.0/pi
!allocate(pot_diff((ndim_x/2-1)*(ndim_y/2-1)*(ndim_z/2-1)))
!no=0
! do i=ndim_x/2+2,ndim_x-1
! do j=ndim_y/2+2,ndim_y-1
!  do k=ndim_z/2+2,ndim_z-1
!	no=no+1
!        pot_diff(no)=abs(exact_pot(i,j,k)-charge_den(i,j,k))
!   enddo
! enddo
!enddo
! write(6,'(F6.2,6X,E10.3,6X,E10.3)') kappa,maxval(pot_diff),minval(pot_diff)
endif
if (choice==2)then
  do i=1,ndim_x
 do j=1,ndim_y
  do k=1,ndim_z
	no=no+1
	 r=dsqrt(xi(i)*xi(i)+yi(j)*yi(j)+zi(k)*zi(k))
	if (r<=rmax)then
   	exact_pot(i,j,k)=charge*(3.0-r*r/rmax/rmax)*0.5/rmax
	else
    	 exact_pot(i,j,k)=charge/r
	endif

	 pot_diff(no)=abs(exact_pot(i,j,k)-hartree_pot(i,j,k))
  enddo
 enddo
enddo
write(6,'(F6.2,6X,E10.3,6X,E10.3)') rmax,maxval(pot_diff),minval(pot_diff)
endif
if (choice==3)then
deallocate (pot_diff)
 call del2(ndim_x,ndim_y,ndim_z,delta_x,hartree_pot,exact_pot)
 exact_pot=-exact_pot/4.0/pi
allocate(pot_diff((ndim_x/2-1)*(ndim_y/2-1)*(ndim_z/2-1)))
 do i=ndim_x/2+2,ndim_x
 do j=ndim_y/2+2,ndim_y
  do k=ndim_z/2+2,ndim_z
	no=no+1
        pot_diff(no)=abs(exact_pot(i,j,k)-charge_den(i,j,k))
   enddo
 enddo
enddo
 write(6,'(F6.2,6X,E10.3,6X,E10.3)') kappa,maxval(pot_diff),minval(pot_diff)
endif
!allocate (pot_diff(ndim_x*ndim_y*ndim_z))
if (choice==4)then
no=0
!factor=sum(charge_den)!*delta_x*delta_y*delta_z
do i=1,ndim_x
 do j=1,ndim_y
  do k=1,ndim_z
	no=no+1
	! r=dsqrt(xi(i)*xi(i)+yi(j)*yi(j)+zi(k)*zi(k))
	exact_pot(i,j,k)= dexp(-kappa*dsqrt(xi(i)*xi(i)+yi(j)*yi(j)+zi(k)*zi(k))**m)!/factor
	pot_diff(no)=abs(exact_pot(i,j,k)-hartree_pot(i,j,k))
  enddo
 enddo
enddo
write(6,'(6X,I4,6X,F6.2,6X,E10.3,6X,E10.3)')m,kappa,maxval(pot_diff),minval(pot_diff)

endif

return
end subroutine check_result

