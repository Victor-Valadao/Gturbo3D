subroutine writef(u1,u2,u3,ifr)
use paran
implicit none

integer ifr,i,j,k
real*8, dimension(NXP2,NY,NZ) :: u1,u2,u3
real*4, allocatable :: w(:,:,:)
character*40 nome
 
allocate (w(NX,NY,NZ))
 
write(nome,"('./fields/ux.',i3.3)")ifr
open(unit=99,file=nome,form='unformatted')   

do k=1,NZ
 do j=1,NY
  do i=1,NX
   w(i,j,k)=real(u1(i,j,k))
  end do
 end do
end do
write(99)w

close(99)
write(nome,"('./fields/uy.',i3.3)")ifr
open(unit=99,file=nome,form='unformatted')

do k=1,NZ
 do j=1,NY
  do i=1,NX
   w(i,j,k)=real(u2(i,j,k))
  end do
 end do
end do
write(99)w

close(99)
write(nome,"('./fields/uz.',i3.3)")ifr
open(unit=99,file=nome,form='unformatted') 

do k=1,NZ
 do j=1,NY
  do i=1,NX
   w(i,j,k)=real(u3(i,j,k))
  end do
 end do
end do
write(99)w

close(99)

deallocate(w)
return
end
! ---------------------------------------------------------------
subroutine writef2(u1,u2,u3,ifr,jfr)
use paran
implicit none

integer ifr,jfr,i,j,k
real*8, dimension(NXP2,NY,NZ) :: u1,u2,u3
real*4, allocatable :: w(:,:,:)
character*40 nome
 
allocate (w(NX,NY,NZ))
 
write(nome,"('./fields/ux.',i3.3,'.',i3.3)")ifr,jfr
open(unit=99,file=nome,form='unformatted')   

do k=1,NZ
 do j=1,NY
  do i=1,NX
   w(i,j,k)=real(u1(i,j,k))
  end do
 end do
end do
write(99)w

close(99)
write(nome,"('./fields/uy.',i3.3,'.',i3.3)")ifr,jfr
open(unit=99,file=nome,form='unformatted')   

do k=1,NZ
 do j=1,NY
  do i=1,NX
   w(i,j,k)=real(u2(i,j,k))
  end do
 end do
end do
write(99)w

close(99)
write(nome,"('./fields/uz.',i3.3,'.',i3.3)")ifr,jfr
open(unit=99,file=nome,form='unformatted') 

do k=1,NZ
 do j=1,NY
  do i=1,NX
   w(i,j,k)=real(u3(i,j,k))
  end do
 end do
end do
write(99)w

close(99)

deallocate(w)
return
end