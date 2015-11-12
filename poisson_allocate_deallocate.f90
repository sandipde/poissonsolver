subroutine poisson_allocate()
use kind_param
use global_var
use poisson_var
implicit none
!____________ALLOCATIONS_____________________________________________

allocate (d_x(-ndim_x+1:ndim_x))
allocate (d_y(-ndim_y+1:ndim_y))
allocate (d_z(-ndim_z+1:ndim_z))
allocate (double_charge_den(-ndim_x+1:ndim_x,-ndim_y+1:ndim_y,-ndim_z+1:ndim_z))
allocate (gspace_charge_den(ndim_x+1,2*ndim_y,2*ndim_z))
allocate (r_dfg_vhartree  (2*ndim_x, 2*ndim_y, 2*ndim_z))
allocate (c_dfg_vhartree  (ndim_x+1,2*ndim_y,2*ndim_z))
return
end subroutine poisson_allocate



!_____________________________________________________________________________________________


subroutine poisson_deallocate()
use kind_param
use global_var
use poisson_var
implicit none
!____________DEALLOCATIONS_____________________________________________

deallocate (double_charge_den)
deallocate (gspace_charge_den)
deallocate (r_dfg_vhartree  )
deallocate (c_dfg_vhartree )
deallocate (vhartree)
return
end subroutine poisson_deallocate
