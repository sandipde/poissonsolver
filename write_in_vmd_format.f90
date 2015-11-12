subroutine write_in_vmd_format(function,file_name)
use kind_param
use global_var
implicit none
integer(I4B)::i,j,k
real(R8B)::function(ndim_x,ndim_y,ndim_z)
 character(Len= 20)::file_name
100 Format (3(F8.5,2X))
open (8,file=file_name)
!_____________________________________________________________________________________
! GENERATING THE COMMIN FORMAT FILE FOR VMD
!open(8,file='vmdformat')
write(8,*) "QD"
write(8,*) 1
write(8,100) ndim_x*delta_x,0, 0
write(8,100) 0,ndim_x*delta_x,0
write(8,100) 0, 0,ndim_z*delta_z
write(8,*) 1
write(8,*) "Direct"
write(8,100) ndim_x*delta_x*0.5, ndim_y*delta_y*0.5, ndim_z*delta_z*0.5
write(8,*)""
write(8,*) ndim_x,ndim_y,ndim_z
!_______________________________________________________________________________________

write (8,*)function
return
end subroutine write_in_vmd_format
