subroutine del2(nx,ny,nz,h,V,P)
use kind_param
use global_var
implicit none
integer(I4B)::nx,ny,nz
real(R8B),dimension(nx,ny,nz)::V,P
real(R8B)::h
!real(R8B),allocatable,dimension(:,:,:)::V,P
integer(I4B)::i,j,k
!allocate (V(nx,ny,nz))
!allocate (P(nx,ny,nz))
P=0.0
do i=2,nx-1
 do j=2,ny-1
  do k=2,nz-1
   
   P(i,j,k)=(V(i,j,k+1)+V(i,j,k-1)+V(i,j+1,k)+V(i,j-1,k)+V(i+1,j,k)+V(i-1,j,k)-6.0*V(i,j,k))/(h*h)
  enddo
 enddo
enddo

do i=2,nx
 do j=2,ny
   write(55,'(5F15.5,5X)')xi(i),yi(j),V(i,j,nz/2+1),P(i,j,nz/2+1)
 enddo
 write(55,*)" "
enddo
do i=1,nx
 do j=1,ny
  do k=1,nz
   write(22,*)P(i,j,k)/12.56637061
  enddo
 enddo
enddo
return
end subroutine del2
