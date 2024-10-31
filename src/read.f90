subroutine readf(u1,u2,u3,ifr)
use paran
implicit none

integer ifr,i,j,k
real*8, dimension(NXP2,NY,NZ) :: u1,u2,u3
real*4, allocatable :: w(:,:,:)
character*40 nome

allocate (w(NX,NY,NZ))

write(nome,"('./fields/ux.',i3.3)")ifr
open(unit=90,file=nome,form='unformatted',action='read')
 read(90)w

do k=1,NZ
 do j=1,NY
  do i=1,NX
   u1(i,j,k)=dble(w(i,j,k))
  end do
 end do
end do

close(90)

write(nome,"('./fields/uy.',i3.3)")ifr
open(unit=90,file=nome,form='unformatted',action='read')
 read(90)w

do k=1,NZ
 do j=1,NY
  do i=1,NX
   u2(i,j,k)=dble(w(i,j,k))
  end do
 end do
end do

close(90)

write(nome,"('./fields/uz.',i3.3)")ifr
open(unit=90,file=nome,form='unformatted',action='read')
 read(90)w

do k=1,NZ
 do j=1,NY
  do i=1,NX
   u3(i,j,k)=dble(w(i,j,k))
  end do
 end do
end do

close(90)

do k=1,NZ
 do j=1,NY
  do i=NX+1,NXP2
   u1(i,j,k)=0.d0
   u2(i,j,k)=0.d0
   u3(i,j,k)=0.d0
  end do
 end do
end do

deallocate(w)
return
end
! ---------------------------------------------------------------

subroutine readf2d(z,ifr)
  use paran
  integer i, j,k,ifr
  real*8, parameter :: pi2=2.d0*3.14159265358979d0
  real*8 z(NXP2,NY,NZ)
real*4 w(NX,NY)
character*40 nome

write(nome,"('./w.',i8.8)")ifr
open(unit=90,file=nome,form='unformatted',action='read')
 read(90)w
close(90)

z=0.d0
do k=1,NZ
do j=1,NY
do i=1,NX
 z(i,j,k)=dble(w(i,j))*dcos(pi2*k/NZ)
end do
end do
end do

close(90)

return
end