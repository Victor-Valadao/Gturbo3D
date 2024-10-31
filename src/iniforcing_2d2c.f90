! Forza random su una shell tra k1forc e k2forc
! Gennaio 2003 - Stefano

subroutine iniforcing
use paran
use commk
use commforc
use commphys
implicit none

real*8 k2,k2n,kk
integer i,j,k

nf=0
k=1

i=1
do j=2,NY2
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
  ff(nf)=1.d0/kk
 end if
 end if
enddo

do i=2,NX2
do j=1,NY
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
  ff(nf)=1.d0/kk
 end if
 end if
end do
end do

i=NX2P1
do j=2,NY2
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
  ff(nf)=1.d0/kk
 end if
 end if
enddo

do i=1,nf
 k2=fkxs(ikxf(i))+fkys(ikyf(i))+fkzs(ikzf(i))
 inpe=inpe+k2*(famp*ff(i))**2
 inpz=inpz+k2*k2*(famp*ff(i))**2
end do

write(6,*)' --- 2D 2-component forcing --- '
write(6,*)' # of forced wavenumbers = ',nf
write(6,*)' Mean forcing amplitude  = ',real(famp)
write(6,*)' Forcing range           : ',real(k1f),' < k =<',real(k2f)
write(6,*)' Mean energy input       = ',real(inpe)
write(6,*)' Mean enstrophy input    = ',real(inpz)
write(6,*)' '

call flush(6)

return
end


