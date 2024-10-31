subroutine fft_dir(z,plan)
  use cufft
  use openacc
  use paran
  implicit none

  real*8, dimension(NXP2,NY,NZ) :: z
  integer :: plan
  integer :: i,j,k,ierr
  
  !$acc host_data use_device(z)
  ierr = ierr + cufftExecD2Z(plan,z,z)
  !$acc end host_data
  
  ! In this particular case, single loop 
  ! improves performance because of memory organization
  !$acc parallel loop collapse(3) present(z)
  do k=1,NZ
  do j=1,NY
  do i=1,NXP2
  	z(i,j,k) = z(i,j,k)/NY/NX/NZ
  enddo
  enddo
  enddo
  !$acc end parallel loop

return
end
