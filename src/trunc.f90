subroutine trunc(z)
use paran
use commk
implicit none

complex*16 z(NX2P1,NY,NZ)
real*8 k1,k2,k3
integer i,j,k

!$acc parallel loop collapse(3) present(z, fkxs, fkys, fkzs)
do k=1,NZ
 do j=1,NY
  do i=1,NX2P1
   k3=fkzs(k)/fkzsmax
   k2=k3+fkys(j)/fkysmax
   k1=k2+fkxs(i)/fkxsmax
   if (k1.gt.rtrunc) then
    z(i,j,k)=(0.d0,0.d0)
   end if
  end do
 end do
end do
!$acc end parallel

! Substract a possible existing mean z
!$acc kernels present(z)
z(1,1,1)=(0.d0,0.d0)
!$acc end kernels

return
end

! --------------------------------------------------------------	
! Sub to control incompressibility
subroutine incompr(ax,ay,az)
use paran
use commk
implicit none

complex*16, dimension(NX2P1,NY,NZ) :: ax,ay,az
complex*16 scra
real*8 k2,cc
integer i,j,k

cc=0.d0

!$acc parallel loop collapse(3) present(ax,ay,az,fkx,fky,fkz) 
do 10 k=1,NZ
do 10 j=1,NY
do 10 i=1,NX2P1

if ((i.eq.1).and.(j.eq.1).and.(k.eq.1)) goto 10

k2=fkxs(i)+fkys(j)+fkzs(k)
scra=fkx(i)*ax(i,j,k)+fky(j)*ay(i,j,k)+fkz(k)*az(i,j,k)
cc=cc+cdabs(scra)**2

ax(i,j,k)=ax(i,j,k)-scra*fkx(i)/k2
ay(i,j,k)=ay(i,j,k)-scra*fky(j)/k2
az(i,j,k)=az(i,j,k)-scra*fkz(k)/k2

10 continue
!$acc end parallel
write(6,*)' Initial compressibility = ',real(cc)
call flush(6)

return
end

! --------------------------------------------------------------	
! Sub to initializate with a small (compressible) noise in high wavenumber
subroutine init_noise(bx,by,bz,seed)
use paran
use commk
use commphys
use commforc
implicit none

real*8, dimension(NXP2,NY,NZ) :: bx,by,bz
integer i,j,k,seed
real :: rann


do k=1,NZ
do j=1,NY
do i=1,NX
 bx(i,j,k)=bx(i,j,k)+sdt*(rann(seed)-0.5)
 by(i,j,k)=by(i,j,k)+sdt*(rann(seed)-0.5)
 bz(i,j,k)=bz(i,j,k)+sdt*(rann(seed)-0.5)
enddo
enddo
enddo

return 
end subroutine


! --------------------------------------------------------------	
! Sub to print memory

subroutine show_mem(fre, tot, step)
use cudafor
implicit none

integer :: ierr,step
integer(kind=cuda_count_kind) :: fre,tot
real :: r

ierr=cudaMemGetInfo( fre, tot )
r=(tot-fre)/2.0**30

write(6,"(A25,I4,A3,F12.9,A4,F6.3)")"Alloc mem (Gb) in step", step, " = ",r," of ",tot/2.0**30
call flush(6)

tot=0
fre=0

end
! --------------------------------------------------------------	
! Sub to print memory

subroutine tt(ax,ay,az,r)
use paran
implicit none

complex*16, dimension(NX2P1,NY,NZ) :: ax,ay,az
integer :: r

!$acc update host(ax,ay,az)
write(6,*)r," sum", sqrt(sum(cdabs(ax)**2+cdabs(ay)**2+cdabs(az)**2))

end