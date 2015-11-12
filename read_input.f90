
!__________________________________________________________________________
!   PROGRAM BY SANDIP DE      26.01.09
!__________________________________________________________________________
subroutine read_input()
use kind_param
use global_var
implicit none
open(1,file="mesh.in")
read(1,*);read(1,*) ndim_x,ndim_y,ndim_z
read(1,*);read(1,*) delta_x,delta_y,delta_z
close (1)
return
end subroutine read_input
