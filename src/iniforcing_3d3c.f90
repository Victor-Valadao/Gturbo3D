! Forza random su una shell tra k1f e k2f
! Gennaio 2003 - Stefano

subroutine iniforcing
use paran
use commk
use commforc
use commphys
implicit none

real*8 k2,k2n,kk
real*8 alfa
integer i,j,k

! power of the forcing spectrum
alfa=0.0

nf=0
inpe=0.d0
inpz=0.d0

i=1
j=1
do k=2,NZ2

 k2 = fkxs(i) + fkys(j) + fkzs(k)
 kk=sqrt(k2)
 k2n=fkxs(i)/fkxsmax+fkys(j)/fkysmax+fkzs(k)/fkzsmax

 if (k2n.le.rtrunc) then
 if ((k2.le.k2f**2).and.(k2.gt.k1f**2)) then
  nf=nf+1
  if (nf.gt.NFORC) stop 'NFORC too small !'
  ikxf(nf)=i
  ikyf(nf)=j
  ikzf(nf)=k
  ff(nf)=1.0d0/kk
  inpz=inpz+k2*k2*(famp*ff(nf))**2
  inpe=inpe+k2*(famp*ff(nf))**2
 end if
 end if

enddo

do j=2,NY2
do k=1,NZ

 k2 = fkxs(i) + fkys(j) + fkzs(k)
 kk=sqrt(k2)
 k2n=fkxs(i)/fkxsmax+fkys(j)/fkysmax+fkzs(k)/fkzsmax

 if (k2n.le.rtrunc) then
 if ((k2.le.k2f**2).and.(k2.gt.k1f**2)) then
  nf=nf+1
  if (nf.gt.NFORC) stop 'NFORC too small !'
  ikxf(nf)=i
  ikyf(nf)=j
  ikzf(nf)=k
  ff(nf)=1.0d0/kk
  inpz=inpz+k2*k2*(famp*ff(nf))**2
  inpe=inpe+k2*(famp*ff(nf))**2
 end if
 end if

enddo
enddo

j=NY2P1
do k=2,NZ2

 k2 = fkxs(i) + fkys(j) + fkzs(k)
 kk=sqrt(k2)
 k2n=fkxs(i)/fkxsmax+fkys(j)/fkysmax+fkzs(k)/fkzsmax

 if (k2n.le.rtrunc) then
 if ((k2.le.k2f**2).and.(k2.gt.k1f**2)) then
  nf=nf+1
  if (nf.gt.NFORC) stop 'NFORC too small !'
  ikxf(nf)=i
  ikyf(nf)=j
  ikzf(nf)=k
  ff(nf)=1.0d0/kk
  inpz=inpz+k2*k2*(famp*ff(nf))**2
  inpe=inpe+k2*(famp*ff(nf))**2
 end if
 end if

enddo

! -----------------------------------------------------------
do i=2,NX2
do j=1,NY
do k=1,NZ

 k2 = fkxs(i) + fkys(j) + fkzs(k)
 kk=sqrt(k2)
 k2n=fkxs(i)/fkxsmax+fkys(j)/fkysmax+fkzs(k)/fkzsmax

 if (k2n.le.rtrunc) then
 if ((k2.le.k2f**2).and.(k2.gt.k1f**2)) then
  nf=nf+1
  if (nf.gt.NFORC) stop 'NFORC too small !'
  ikxf(nf)=i
  ikyf(nf)=j
  ikzf(nf)=k
  ff(nf)=1.0d0/kk
  inpz=inpz+k2*k2*(famp*ff(nf))**2
  inpe=inpe+k2*(famp*ff(nf))**2
 end if
 end if

end do
end do
end do
! -----------------------------------------------------------

i=NX2P1

j=1
do k=2,NZ2

 k2 = fkxs(i) + fkys(j) + fkzs(k)
 kk=sqrt(k2)
 k2n=fkxs(i)/fkxsmax+fkys(j)/fkysmax+fkzs(k)/fkzsmax

 if (k2n.le.rtrunc) then
 if ((k2.le.k2f**2).and.(k2.gt.k1f**2)) then
  nf=nf+1
  if (nf.gt.NFORC) stop 'NFORC too small !'
  ikxf(nf)=i
  ikyf(nf)=j
  ikzf(nf)=k
  ff(nf)=1.0d0/kk
  inpz=inpz+k2*k2*(famp*ff(nf))**2
  inpe=inpe+k2*(famp*ff(nf))**2
 end if
 end if

enddo

do j=2,NY2
do k=1,NZ

 k2 = fkxs(i) + fkys(j) + fkzs(k)
 kk=sqrt(k2)
 k2n=fkxs(i)/fkxsmax+fkys(j)/fkysmax+fkzs(k)/fkzsmax

 if (k2n.le.rtrunc) then
 if ((k2.le.k2f**2).and.(k2.gt.k1f**2)) then
  nf=nf+1
  if (nf.gt.NFORC) stop 'NFORC too small !'
  ikxf(nf)=i
  ikyf(nf)=j
  ikzf(nf)=k
  ff(nf)=1.0d0/kk
  inpz=inpz+k2*k2*(famp*ff(nf))**2
  inpe=inpe+k2*(famp*ff(nf))**2
 end if
 end if

enddo
enddo

j=NY2P1
do k=2,NZ2

 k2 = fkxs(i) + fkys(j) + fkzs(k)
 kk=sqrt(k2)
 k2n=fkxs(i)/fkxsmax+fkys(j)/fkysmax+fkzs(k)/fkzsmax

 if (k2n.le.rtrunc) then
 if ((k2.le.k2f**2).and.(k2.gt.k1f**2)) then
  nf=nf+1
  if (nf.gt.NFORC) stop 'NFORC too small !'
  ikxf(nf)=i
  ikyf(nf)=j
  ikzf(nf)=k
  ff(nf)=1.0d0/kk
  inpz=inpz+k2*k2*(famp*ff(nf))**2
  inpe=inpe+k2*(famp*ff(nf))**2
 end if
 end if

enddo

write(6,*)' --- 3D 3-component forcing --- '
write(6,*)' # of forced wavenumbers = ',nf
write(6,*)' Mean forcing amplitude  = ',real(famp)
write(6,*)' Forcing range           : ',real(k1f),' < k =<',real(k2f)
write(6,*)' Mean energy input       = ',real(inpe)
write(6,*)' Mean enstrophy input    = ',real(inpz)
write(6,*)' '

return
end
