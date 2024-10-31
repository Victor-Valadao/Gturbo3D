subroutine step1(ax,ay,az,cx,cy,cz)
use paran
use commk
use commphys
use commsim
implicit none

complex*16, dimension(NX2P1,NY,NZ) :: ax,ay,az,cx,cy,cz
integer i,j,k
real*8 scra

! Half step of the time evolution
if(split.eq.1) then
	
	!$acc parallel loop collapse(3) present(ax,ay,az,cx,cy,cz,emkx,emky,emkz,dt2)
	do k=1,NZ
	do j=1,NY
	do i=1,NX2P1
	 scra=emkx(i)*emky(j)*emkz(k)
	 cx(i,j,k)=scra*( ax(i,j,k) + dt2*cx(i,j,k) )
	 cy(i,j,k)=scra*( ay(i,j,k) + dt2*cy(i,j,k) )
	 cz(i,j,k)=scra*( az(i,j,k) + dt2*cz(i,j,k) ) 
	end do
	end do
	end do
	!$acc end parallel 
	
	else
	
	!$acc parallel loop collapse(3) present(ax,ay,az,cx,cy,cz,emk,dt2)
	do k=1,NZ
	do j=1,NY
	do i=1,NX2P1
	 cx(i,j,k)=emk(i,j,k)*( ax(i,j,k) + dt2*cx(i,j,k) )
	 cy(i,j,k)=emk(i,j,k)*( ay(i,j,k) + dt2*cy(i,j,k) )
	 cz(i,j,k)=emk(i,j,k)*( az(i,j,k) + dt2*cz(i,j,k) ) 
	end do
	end do
	end do
	!$acc end parallel 
	
end if

!$acc kernels present(cx,cy,cz)
cx(1,1,1)=0.d0
cy(1,1,1)=0.d0
cz(1,1,1)=0.d0
!$acc end kernels

return
end

! ================================================================================
! ================================================================================
! ================================================================================

subroutine step2(ax,ay,az,cx,cy,cz)
use paran
use commk
use commphys
use commsim
implicit none

complex*16, dimension(NX2P1,NY,NZ) :: ax,ay,az,cx,cy,cz
integer i,j,k
real*8 scra

! Last step of the time evolution 
if(split.eq.1) then
	
	!$acc parallel loop collapse(3) present(ax,ay,az,cx,cy,cz,emkx,emky,emkz,dt)
	do k=1,NZ
	do j=1,NY
	do i=1,NX2P1
	 scra=emkx(i)*emky(j)*emkz(k)
	 ax(i,j,k)=scra*( scra*ax(i,j,k) + dt*cx(i,j,k) )
	 ay(i,j,k)=scra*( scra*ay(i,j,k) + dt*cy(i,j,k) )
	 az(i,j,k)=scra*( scra*az(i,j,k) + dt*cz(i,j,k) ) 
	end do
	end do
	end do
	!$acc end parallel 
	
	else
	
	!$acc parallel loop collapse(3) present(ax,ay,az,cx,cy,cz,emk,dt)
	do k=1,NZ
	do j=1,NY
	do i=1,NX2P1
	 ax(i,j,k)=emk(i,j,k)*( emk(i,j,k)*ax(i,j,k) + dt*cx(i,j,k) )
	 ay(i,j,k)=emk(i,j,k)*( emk(i,j,k)*ay(i,j,k) + dt*cy(i,j,k) )
	 az(i,j,k)=emk(i,j,k)*( emk(i,j,k)*az(i,j,k) + dt*cz(i,j,k) ) 
	end do
	end do
	end do
	!$acc end parallel 
	
end if

!$acc kernels present(ax,ay,az)
ax(1,1,1)=0.d0
ay(1,1,1)=0.d0
az(1,1,1)=0.d0
!$acc end kernels

return
end

! ================================================================================
! ================================================================================
! ================================================================================

subroutine step(bx,by,bz,bnlx,bnly,bnlz,ux,uy,uz)
use paran
use commk
use commphys
use commsim
implicit none

complex*16, dimension(NX2P1,NY,NZ) :: bx,by,bz,bnlx,bnly,bnlz,ux,uy,uz
integer i,j,k
real*8 scra

scra=0.d0

!$acc kernels present(bx,by,bz,bnlx,bnly,bnlz)
bnlx=bx
bnly=by  
bnlz=bz
!$acc end kernels

! Calculate the nonlinear part
call nlt(bnlx,bnly,bnlz,ux,uy,uz)

! Half step of the time evolution
if(split.eq.1) then
	
	!$acc parallel loop collapse(3) present(bx,by,bz,bnlx,bnly,bnlz,emkx,emky,emkz,dt2)
	do k=1,NZ
	do j=1,NY
	do i=1,NX2P1
	 scra=emkx(i)*emky(j)*emkz(k)
	 bnlx(i,j,k)=scra*( bx(i,j,k) + dt2*bnlx(i,j,k) )
	 bnly(i,j,k)=scra*( by(i,j,k) + dt2*bnly(i,j,k) )
	 bnlz(i,j,k)=scra*( bz(i,j,k) + dt2*bnlz(i,j,k) ) 
	end do
	end do
	end do
	!$acc end parallel 
	
	else
	
	!$acc parallel loop collapse(3) present(bx,by,bz,bnlx,bnly,bnlz,emk,dt2)
	do k=1,NZ
	do j=1,NY
	do i=1,NX2P1
	 bnlx(i,j,k)=emk(i,j,k)*( bx(i,j,k) + dt2*bnlx(i,j,k) )
	 bnly(i,j,k)=emk(i,j,k)*( by(i,j,k) + dt2*bnly(i,j,k) )
	 bnlz(i,j,k)=emk(i,j,k)*( bz(i,j,k) + dt2*bnlz(i,j,k) ) 
	end do
	end do
	end do
	!$acc end parallel 
	
end if

!$acc kernels present(bnlx,bnly,bnlz)
bnlx(1,1,1)=0.d0
bnly(1,1,1)=0.d0
bnlz(1,1,1)=0.d0
!$acc end kernels

! Calculate nonlinear part in half step
call nlt(bnlx,bnly,bnlz,ux,uy,uz)

! Last step of the time evolution 
if(split.eq.1) then
	
	!$acc parallel loop collapse(3) present(bx,by,bz,bnlx,bnly,bnlz,emkx,emky,emkz,dt)
	do k=1,NZ
	do j=1,NY
	do i=1,NX2P1
	 scra=emkx(i)*emky(j)*emkz(k)
	 bx(i,j,k)=scra*( scra*bx(i,j,k) + dt*bnlx(i,j,k) )
	 by(i,j,k)=scra*( scra*by(i,j,k) + dt*bnly(i,j,k) )
	 bz(i,j,k)=scra*( scra*bz(i,j,k) + dt*bnlz(i,j,k) ) 
	end do
	end do
	end do
	!$acc end parallel 
	
	else
	
	!$acc parallel loop collapse(3) present(bx,by,bz,bnlx,bnly,bnlz,emk,dt)
	do k=1,NZ
	do j=1,NY
	do i=1,NX2P1
	 bx(i,j,k)=emk(i,j,k)*( emk(i,j,k)*bx(i,j,k) + dt*bnlx(i,j,k) )
	 by(i,j,k)=emk(i,j,k)*( emk(i,j,k)*by(i,j,k) + dt*bnly(i,j,k) )
	 bz(i,j,k)=emk(i,j,k)*( emk(i,j,k)*bz(i,j,k) + dt*bnlz(i,j,k) ) 
	end do
	end do
	end do
	!$acc end parallel 
	
end if

!$acc kernels present(bx,by,bz)
bx(1,1,1)=0.d0
by(1,1,1)=0.d0
bz(1,1,1)=0.d0
!$acc end kernels

return
end







