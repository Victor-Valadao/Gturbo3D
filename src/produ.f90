subroutine cross(a1,a2,a3,b1,b2,b3)
use paran
implicit none

real*8, dimension(NXP2,NY,NZ) :: a1,a2,a3,b1,b2,b3
real*8 :: sc1,sc2,sc3
integer i,j,k

!$acc parallel loop collapse(3) present(a1,a2,a3,b1,b2,b3) 
do k=1,NZ
do j=1,NY
do i=1,NXP2
 if(i.le.NX) then
 sc1=a2(i,j,k)*b3(i,j,k)-a3(i,j,k)*b2(i,j,k)
 sc2=a3(i,j,k)*b1(i,j,k)-a1(i,j,k)*b3(i,j,k)
 sc3=a1(i,j,k)*b2(i,j,k)-a2(i,j,k)*b1(i,j,k)
 a1(i,j,k)=sc1
 a2(i,j,k)=sc2
 a3(i,j,k)=sc3
 else
 	 a1(i,j,k)=0.d0
 	 a2(i,j,k)=0.d0
	 a3(i,j,k)=0.d0
 endif
end do
end do
end do
!$acc end parallel

return
end
! --------------------------------------------------------------

subroutine curl(ax,ay,az)
use paran
use commk
implicit none
 
real*8, intent(inout ), dimension(2,NX2P1,NY,NZ) :: ax,ay,az 
real*8, dimension(2) :: scrax,scray,scraz
real*8 :: k2
integer :: i,j,k     

!$acc parallel loop collapse(3) present(ax,ay,az,fkx,fky,fkz) 
do k=1,NZ
do j=1,NY
do i=1,NX2P1
 scrax(1)=-fky(j)*az(2,i,j,k)+fkz(k)*ay(2,i,j,k)
 scrax(2)= fky(j)*az(1,i,j,k)-fkz(k)*ay(1,i,j,k)
 scray(1)=-fkz(k)*ax(2,i,j,k)+fkx(i)*az(2,i,j,k)
 scray(2)= fkz(k)*ax(1,i,j,k)-fkx(i)*az(1,i,j,k)
 scraz(1)=-fkx(i)*ay(2,i,j,k)+fky(j)*ax(2,i,j,k)
 scraz(2)= fkx(i)*ay(1,i,j,k)-fky(j)*ax(1,i,j,k)
 ax(1,i,j,k)=scrax(1)
 ax(2,i,j,k)=scrax(2)
 ay(1,i,j,k)=scray(1)
 ay(2,i,j,k)=scray(2)
 az(1,i,j,k)=scraz(1)
 az(2,i,j,k)=scraz(2)
end do
end do
end do         
!$acc end parallel

return 
end
! --------------------------------------------------------------

subroutine invcurl(ax,ay,az)
use paran
use commk
implicit none
 
real*8, intent(inout ), dimension(2,NX2P1,NY,NZ) :: ax,ay,az 
real*8, dimension(2) :: scrax,scray,scraz
real*8 :: k2,kx,ky,kz
integer :: i,j,k    
 
!$acc parallel loop collapse(3) present(ax,ay,az,fkx,fky,fkz) 
do 30 k=1,NZ
do 30 j=1,NY
do 30 i=1,NX2P1
 if ((i.eq.1).and.(j.eq.1).and.(k.eq.1)) goto 30
 k2=fkxs(i)+fkys(j)+fkzs(k)
 kx=fkx(i)/k2
 ky=fky(j)/k2
 kz=fkz(k)/k2
 scrax(1)=-ky*az(2,i,j,k)+kz*ay(2,i,j,k)
 scrax(2)= ky*az(1,i,j,k)-kz*ay(1,i,j,k)
 scray(1)=-kz*ax(2,i,j,k)+kx*az(2,i,j,k)
 scray(2)= kz*ax(1,i,j,k)-kx*az(1,i,j,k)
 scraz(1)=-kx*ay(2,i,j,k)+ky*ax(2,i,j,k)
 scraz(2)= kx*ay(1,i,j,k)-ky*ax(1,i,j,k)
 ax(1,i,j,k)=scrax(1)
 ax(2,i,j,k)=scrax(2)
 ay(1,i,j,k)=scray(1)
 ay(2,i,j,k)=scray(2)
 az(1,i,j,k)=scraz(1)
 az(2,i,j,k)=scraz(2)
30 continue        
!$acc end parallel

return 
end

