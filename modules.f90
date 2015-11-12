module kind_param

! * * * * * * * * SYMBOLIC NAMES OF VARIABLE TYPES * * * * * * * *

! 4-, 8-, 2-, and 1-byte integers
  INTEGER, PARAMETER :: I4B=SELECTED_INT_KIND(8)
  INTEGER, PARAMETER :: I8B=SELECTED_INT_KIND(18)
  INTEGER, PARAMETER :: I2B=SELECTED_INT_KIND(4)
  INTEGER, PARAMETER :: I1B=SELECTED_INT_KIND(2)
  
! 8- and 4-byte reals
  INTEGER, PARAMETER :: R8B=SELECTED_REAL_KIND(10)
  INTEGER, PARAMETER :: R4B=SELECTED_REAL_KIND(4)
  
end module kind_param

module global_var
use kind_param
!************MESH*********************************************************

integer(I8B)::ndim_x,ndim_y,ndim_z		                ! grid dimension
real(R8B),DIMENSION(:),allocatable::xi,yi,zi               ! x,y,z points
real(R8B)::delta_x,delta_y,delta_z                                ! Mesh widths
real(R8B)::x_bound,y_bound,z_bound                          ! Bounderies of the Mesh
real(R8B)::Lx,Ly,Lz,cell_vol

!****************************************************************************

!*************Potential*******************************************************
real(R8B),DIMENSION(:,:,:),allocatable::ext_pot
!****************************************************************************

!***********Hamiltoninan****************************************************
integer(I4B),allocatable::nonzero_counter(:),hamilton_index(:)              !hamilton_index stores the indexes
integer(I4B)::compressed_dim,totall_nonzero,hamilton_dim,max_coeff
real(R8B),allocatable::hamilton_elements(:)                                           !stores the elements 

!_______________________________________________________________________________

!------------------------------------------
! DVDSON
   INTEGER(I4B) :: NMAX, NUMEMAX,LIMMAX, IWRSZ,NUME,   &   
                   IIWSZ ,LIM,ILOW,IHIGH, NIV,         &
                   MBLOCK,MAXITER,NLOOPS,NMV,IERR
   REAL(R8B)    :: CRITE , CRITC, CRITR, ORTHO  

   INTEGER(I4B),ALLOCATABLE :: ISELEC(:), IWORK(:)
REAL(R8B), ALLOCATABLE         :: Eigen_vec(:)
   logical ::HIEND

!_____________________________________________________________________
! POISSON
!___________________________________________________________________

real(R8B),allocatable::charge_den(:,:,:),hartree_pot(:,:,:)
real(R8B)::hart2ev, ev2hart, bohr2ang, ang2bohr, a0
parameter  (hart2ev=27.2116d0,   ev2hart=1.0d0/hart2ev)
parameter  (bohr2ang=0.529177d0, ang2bohr= 1.0d0/bohr2ang)
parameter  (a0=0.529177d0)
end module global_var

!################################################################################################################################

module poisson_var
use kind_param
!__________double grid__________________
real(R8B),DIMENSION(:),allocatable::d_x,d_y,d_z
real(R8B)::dx_bound,dy_bound,dz_bound
!_____________Charge densisity_______________________________
real(R8B)::gauss_const,kappa,pi,twopi,fourpi
parameter (pi=3.141592653589793238462643383d0)
parameter(twopi=2.0d0*pi)
parameter(fourpi=4.0d0*pi)
real(R8B),allocatable::fft_charge_den(:,:,:)
real(R8B),allocatable::double_charge_den(:,:,:)      ! charge density in the double grid
integer(I4B)::original_x_start,original_y_start,original_z_start
!_______fft__________________________________
real(R8B),allocatable,dimension(:)::fft_gx,fft_gy,fft_gz
integer(i4B),allocatable,dimension(:)::tr_x_r2f,tr_y_r2f,tr_z_r2f,tr_x_f2r,tr_y_f2r,tr_z_f2r     !fourier<>real mesh maps
double complex,allocatable::gspace_charge_den(:,:,:),c_dfg_vhartree(:,:,:)
real(R8B),allocatable::r_dfg_vhartree(:,:,:),vhartree (:,:,:)
integer(I4B)::plan,invplan,i0,j0,k0
real(R8B)::delta_r
real(R8B),allocatable:: gsqr (:,:,:)
real(R8B)::tpibyLx, tpibyLy, tpibyLz
real(R8B)::gz2, gyz2, g_sqr, g_sqr_max
real(R8B)::four_alpha_sqr, g_exp_term ,fourpi_by_gsqr,alpha_sqr
real(R8B)::hartree_energy
double complex:: vtermg, vsumg
real(R8B)alpha_x,alpha_y,alpha_z,alpha,alpha_conv_L
real(r8B)::exact_energy,rmax,r,max_diff,min_diff,charge
real(r8B),allocatable::exact_pot(:,:,:),pot_diff(:) 
integer(i4B)::choice
end module poisson_var

!__________________________________________________________________________
!   PROGRAM BY SANDIP DE      26.01.09
!__________________________________________________________________________
