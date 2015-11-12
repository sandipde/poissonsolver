!__________________________________________________________________________
!   PROGRAM BY SANDIP DE      26.01.09
!__________________________________________________________________________

subroutine charge_density()
use kind_param
use global_var
use poisson_var
implicit none
integer(I4B)::i,j,k
real(8)::tmp
!kappa  =  (0.944863d0)/ang2bohr
!kappa  =  a0/2.0
gauss_const=dsqrt((kappa*kappa/pi)**3.0)
!write(*,*)"g_const=",gauss_const
do i=1,ndim_x
  do j=1,ndim_y
    do k=1,ndim_z
      charge_den(i,j,k)=gauss_const*dexp(-kappa*kappa*(xi(i)*xi(i)+yi(j)*yi(j)+zi(k)*zi(k)))
    enddo
   enddo
enddo
write(66,*)"========================CHARGE DENSITY=================================="
tmp=kappa*min(ndim_x*delta_x,ndim_y*delta_y,ndim_z*delta_z)
write(66,'(10X,"TYPE :: GAUSSIAN" ,/, 10X,"kappa = ",F10.4,/, 10X,"Kappa*L = ",F10.3)')kappa,tmp
write(66,*)"======================================================================="
open (1,file='charge_den.dat')
do i=1,ndim_x
  do j=1,ndim_y
    write(1,200)xi(i),yi(j),zi(ndim_z),charge_den(i,j,ndim_z/2)
  enddo 
    write(1,*)"  "
enddo
200 Format (4(F15.7,2X))
end subroutine charge_density

!=================================================================================================

subroutine gen_charged_sphere()
use kind_param
use global_var
use poisson_var
implicit none
integer(I4B)::i,j,k
real(R8B)::vol_sphere!r,rmax,
!rmax=5.0
vol_sphere=4.0/3.0*pi*rmax*rmax*rmax
 charge_den=0.0
do i=1,ndim_x
  do j=1,ndim_y
    do k=1,ndim_z
	r=dsqrt(xi(i)*xi(i)+yi(j)*yi(j)+zi(k)*zi(k))
        if (r.le.rmax)then
     !        charge_den(i,j,k)=1.0/(vol_sphere)*4.0/3.0*pi*r**3/(vol_sphere)
        charge_den(i,j,k)=charge/(vol_sphere)
	endif 
    enddo
   enddo
enddo
! charge_den=0.0
open (1,file='charge_den.dat')
do i=1,ndim_x
  do j=1,ndim_y
   do k=1,ndim_z
       write(1,*) dsqrt(xi(i)*xi(i)+yi(j)*yi(j)+zi(k)*zi(k)),charge_den(i,j,k)
   enddo
  enddo
  enddo
write(66,*)"========================CHARGE DENSITY=================================="
write(66,'(10X,"TYPE :: CHARGED SPHERE" ,/, 10X,"RADIUS = ",F10.4,/, 10X,"CHARGE = ",F10.3)')rmax,charge
write(66,*)"======================================================================="
end subroutine gen_charged_sphere

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine exp_density()
use kind_param
use global_var
use poisson_var
implicit none
integer(I4B)::i,j,k
real(R8B)::V(ndim_x,ndim_y,ndim_z),tmp
do i=1,ndim_x
 do j=1,ndim_y
  do k=1,ndim_z
    charge_den(i,j,k)=dexp(-kappa*dsqrt((xi(i)*xi(i)+yi(j)*yi(j)+zi(k)*zi(k))))
  enddo
 enddo
enddo
! charge_den=charge_den/(sum(charge_den))

! call del2(ndim_x,ndim_y,ndim_z,delta_x,V,charge_den)
! charge_den=-charge_den/(4.0*pi)
open (1,file='charge_den.dat')
do i=1,ndim_x
  do j=1,ndim_y
    write(1,200)xi(i),yi(j),charge_den(i,j,ndim_z/2)
  enddo 
    write(1,*)"  "
enddo
200 Format (3(F15.7,2X))
write(66,*)"========================CHARGE DENSITY=================================="
tmp=kappa*min(ndim_x*delta_x,ndim_y*delta_y,ndim_z*delta_z)
write(66,'(10X,"TYPE :: EXPONENTIAL" ,/, 10X,"kappa = ",F10.4,/, 10X,"Kappa*L = ",F10.3)')kappa,tmp
write(66,*)"======================================================================="
return
end subroutine exp_density

subroutine general_exp(m)
use kind_param
use global_var
use poisson_var
implicit none
integer(i4B)::i,j,k,m
do i=1,ndim_x
 do j=1,ndim_y
    do k=1,ndim_z
	r=dsqrt(xi(i)*xi(i)+yi(j)*yi(j)+zi(k)*zi(k))
	charge_den(i,j,k)=(kappa*m*r**(m-2))*(m+1-kappa*m*r**m)*exp(-kappa*r**m)
    enddo
  enddo
enddo
 charge_den=charge_den/(4.0*pi)!(sum(charge_den))
open (1,file='charge_den.dat')
do i=1,ndim_x
  do j=1,ndim_y
    write(1,200)xi(i),yi(j),charge_den(i,j,ndim_z/2)
  enddo 
    write(1,*)"  "
enddo
200 Format (3(F15.7,2X))
return
end subroutine general_exp
