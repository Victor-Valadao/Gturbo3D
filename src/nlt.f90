subroutine nlt(ax,ay,az,vx,vy,vz,bx,by,bz)
use paran
use cufft
use openacc
use cudafor
use stream
use plan
use commk
use commphys
implicit none

real*8, dimension(2,NX2P1,NY,NZ) :: ax,ay,az,vx,vy,vz,bx,by,bz
real*8 k2,kx,ky,kz
integer i,j,k

!$acc kernels present(vx,vy,vz,ax,ay,az)
vx=ax
vy=ay
vz=az
!$acc end kernels
!call tt(vx,vy,vz,10)
! Calculate u = curl(b) 
call curl(vx,vy,vz)
!call tt(vx,vy,vz,11)
! Calculate w=-laplacian(b)  (a=-w)
!$acc parallel loop collapse(3) present(ax,ay,az,fkxs,fkys,fkzs)
do k=1,NZ
do j=1,NY
do i=1,NX2P1
 k2=fkxs(i)+fkys(j)+fkzs(k)
 ax(1,i,j,k)=k2*ax(1,i,j,k)
 ax(2,i,j,k)=k2*ax(2,i,j,k)
 ay(1,i,j,k)=k2*ay(1,i,j,k)
 ay(2,i,j,k)=k2*ay(2,i,j,k)
 az(1,i,j,k)=k2*az(1,i,j,k)
 az(2,i,j,k)=k2*az(2,i,j,k)
end do
end do
end do
!$acc end parallel
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!call tt(ax,ay,az,12)
! Go to real space
call fft_inv(vx, plan_inv)
call fft_inv(vy, plan_inv)
call fft_inv(vz, plan_inv)
call fft_inv(ax, plan_inv)
call fft_inv(ay, plan_inv)
call fft_inv(az, plan_inv)

!call tt(vx,vy,vz,13)
!call tt(ax,ay,az,14)

! Calulate the cross product u x w and save in the first 3 entries (v)
call cross(vx,vy,vz,ax,ay,az)
!call tt(vx,vy,vz,15)

! Cross product in the spectral space
call fft_dir(vx, plan_dir)
call fft_dir(vy, plan_dir)
call fft_dir(vz, plan_dir)

!call tt(vx,vy,vz,16)

! Calculate NLT = -laplacian^{-1}(curl(u x w))

!$acc parallel loop collapse(3) present(vx,vy,vz,ax,ay,az,bx,by,bz,fkz)
do k=1,NZ
do j=1,NY
do i=1,NX2P1
 k2=fkz(k)
 ax(1,i,j,k)=vx(1,i,j,k)-omega*k2*bx(2,i,j,k)
 ax(2,i,j,k)=vx(2,i,j,k)+omega*k2*bx(1,i,j,k)
 ay(1,i,j,k)=vy(1,i,j,k)-omega*k2*by(2,i,j,k)
 ay(2,i,j,k)=vy(2,i,j,k)+omega*k2*by(1,i,j,k)
 az(1,i,j,k)=vz(1,i,j,k)-omega*k2*bz(2,i,j,k)
 az(2,i,j,k)=vz(2,i,j,k)+omega*k2*bz(1,i,j,k)
end do
end do
end do
!$acc end parallel

call invcurl(ax,ay,az)

if (truncate.eq.1) then
 call trunc(ax)
 call trunc(ay) 
 call trunc(az)
end if

!call tt(ax,ay,az,17)
!call tt(vx,vy,vz,18)

return
end
! --------------------------------------------------------------