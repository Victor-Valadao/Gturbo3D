subroutine defk
use paran
use commk
use commphys
implicit none

real*8 fk2
real*8 scra
integer i,j,k

! Conjugate variables: jc kc

jc(1)=1
do j=2,NY
 jc(j)=NY+2-j
end do

kc(1)=1
do k=2,NZ
 kc(k)=NZ+2-k
end do

! Wavenumber increments
dkx=2.0d0*3.14159265358979323846d0/xlx
dky=2.0d0*3.14159265358979323846d0/xly
dkz=2.0d0*3.14159265358979323846d0/xlz

! Wavenumbers
do i=1,NX2P1
 fkx(i)=dble(i-1)*dkx
end do
do j=1,NY2P1
 fky(j)=dble(j-1)*dky
end do
do j=NY2P2,NY
 fky(j)=dble(j-NY-1)*dky
end do
do k=1,NZ2P1
 fkz(k)=dble(k-1)*dkz
end do
do k=NZ2P2,NZ
 fkz(k)=dble(k-NZ-1)*dkz
end do

! Wavenumbers square
do i=1,NX2P1
 fkxs(i)=fkx(i)*fkx(i)
end do
do j=1,NY
 fkys(j)=fky(j)*fky(j)
end do
do k=1,NZ
 fkzs(k)=fkz(k)*fkz(k)
end do

! Dissipative terms
emk=0.d0
do k=1,NZ
do j=1,NY
do i=1,NX2P1
 fk2=fkxs(i)+fkys(j)+fkzs(k)
 scra=0.d0
 if (fk2.ne.0.0) scra=nu*fk2**alpnu+mu*fk2**alpmu
 emk(i,j,k)=dexp(-0.5*dt*scra)
 if (fk2.eq.0.0) emk(i,j,k)=1.0
end do
end do
end do

do k=1,NZ
 fk2=fkzs(k)
 scra=nu*fk2**alpnu+mu
 emkz(k)=dexp(-0.5d0*dt*scra)  
enddo

do j=1,NY
 fk2=fkys(j)
 scra=nu*fk2**alpnu
 emky(j)=dexp(-0.5d0*dt*scra)  
enddo

do i=1,NX2P1
 fk2=fkxs(i)
 scra=nu*fk2**alpnu
 emkx(i)=dexp(-0.5d0*dt*scra)  
enddo

fkxsmax=fkxs(NX2P1)
fkysmax=fkys(NY2P1)
fkzsmax=fkzs(NZ2P1)
rtrunc=rtrunc*rtrunc

return
end
! --------------------------------------------------------------	

