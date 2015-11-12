!__________________________________________________________________________
!   PROGRAM BY SANDIP DE      26.01.09
!__________________________________________________________________________

program poisson
use kind_param
use global_var
use poisson_var
implicit none
integer(I4B)::i,j,k,ix,iy,iz,i1,no,m
!real(r8B)::exact_energy,rmax,r,max_diff,min_diff
!real(r8B),allocatable::exact_pot(:,:,:),pot_diff(:)
 open(66,file='poisson.log')
 call read_input()
 write(66,*)"MIN RAM REQUIRED > ",(3.0*(ndim_x+ndim_y+ndim_z)+35.0*ndim_x*ndim_y*ndim_z+16.0*ndim_y*ndim_z)*8.0/(1024*1024),"MB"
 include 'allocate'
 call poisson_allocate()
allocate(exact_pot(ndim_x,ndim_y,ndim_z))
alpha_conv_L=5.0d0
 call gen_mesh()

 
!write(6,*)"plz enter your choice:: 1= gaussian  2=charged sphere 4= gen_exp density"
 !read(*,*) choice 
  choice=1
  if (choice==1)then
!  write(*,*)"give kappa "
  read(*,*)kappa
  call charge_density()
  exact_energy=dsqrt(2.0d0/pi)*kappa
  endif
  if (choice==2)then
!   write(*,*)"give radius of the sphere(no of points)"
   read(*,*)rmax
   rmax=rmax*min(delta_x,delta_y,delta_z)
    write(6,*)"radius=" ,rmax
!   write(6,*)"charge = ?"
!   read(6,*) charge
	charge=0.1
   call gen_charged_sphere()
   exact_energy=1.0/(pi)/(4.0/3.0*pi*rmax**3)
  endif
  if (choice==3)then
   write(*,*)"give kappa "
  read(*,*)kappa
   call exp_density()
  endif
  if (choice==4)then
  write(*,*)"give kappa "
  read(*,*)kappa
	kappa=5.0
  write(*,*)"give exponent"
  read(*,*)m
  call general_exp(m)
  endif
 call double_grid_trans()
 call fft_grid_trans()
 call gen_hartree_pot()
 call poisson_deallocate()
!write(*,*)sum(charge_den)
write(66,*)"___________________________________OUTPUT_____________________________________"
write(66,'(a35,2x,f15.7)')"Integral rho(r) dr  = ", sum(charge_den)*delta_r
!write(66,'(a35,2x,f15.7)')"Integral Vh(r) dr  = ", sum(hartree_pot)/(Lx*Ly*Lz)!*delta_r
!write(66,'(a35,2x,f15.13)')"Analytical energy ~", exact_energy
write(66,'(a35,2x,f15.13)')"Numerical G-space energy =",hartree_energy
hartree_energy =sum(charge_den*hartree_pot) * delta_r
write(66,'(a35,2x,f15.13)')"Numerical R-space energy  = ",hartree_energy
!write(66,'(a35,2x,f15.7)')"Error(%) =",(exact_energy - hartree_energy)/exact_energy*100.0d0
write(66,*)"____________________________________END_________________________________________"
 
!   call analytic_pot(rmax,exact_pot)
 call check_result(m)
open (5,file="hartree_pot.dat")

100 Format (4F20.15,2x)
  do i=1,ndim_x
  do j=1,ndim_y
       write(5,100) xi(i),yi(j),hartree_pot(i,j,ndim_z/2),exact_pot(i,j,ndim_z/2)
   enddo
   write(5,*)" "
  enddo
! else
!   do i=1,ndim_x
!   do j=1,ndim_y
!    do k=1,ndim_z
!        write(5,*) dsqrt(xi(i)*xi(i)+yi(j)*yi(j)+zi(k)*zi(k)),hartree_pot(i,j,k),hartree_pot(i,j,k)/exact_pot(i,j,k)
!    enddo
!   enddo
!   enddo
! endif
  
  

 !call write_in_vmd_format(hartree_pot,'hartree_vmd.dat')
end program poisson
