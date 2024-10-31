! Forcing random statisticamente isotropo
! incompressible div f = 0
! realizzato scegliendo un vettore complesso 
! di modulo famp*e^(i*theta) con theta random 
! nel piano ortogonale a k
! 1 e l'angolo theta polare di k
! 2 e l'angolo phi azimutale di k  
!       
! Versione NX,NY,NZ
!
subroutine forcing(bx,by,bz,seed)
use commk
use commphys
use commforc
implicit none

complex*16 bx(NX2P1,NY,NZ)
complex*16 by(NX2P1,NY,NZ)
complex*16 bz(NX2P1,NY,NZ)
complex*16 effe,effek,effekc
real*8, parameter :: pi2=2.0d0*3.14159265358979323846d0
real*8 te,ste,cte,k2,kk,kxy,s1,s2,c1,c2
integer i,j,k,ik,seed
real :: rann

!$acc parallel loop seq present(bx,by,bz,jc,kc,ff,fkx,fky,fkz,seed,sdt)
do ik=1,nf

i=ikxf(ik)
j=ikyf(ik)
k=ikzf(ik)

! sceglie la fase random nello spazio complesso 
te=pi2*rann(seed)
effe=sdt*famp*cdexp((0.d0,1.d0)*te)

! sceglie la fase random nel piano ortogonale a k
te=pi2*rann(seed)
cte=dcos(te)
ste=dsin(te)

effek=effe*ff(ik)

! rappresentazione singolare per kx=ky=0
if ((i.eq.1).and.(j.eq.1)) then
 bx(i,j,k)=bx(i,j,k)+effek*cte
 by(i,j,k)=by(i,j,k)+effek*ste
  effekc=conjg(effek)
  bx(i,jc(j),kc(k))=bx(i,jc(j),kc(k))+effekc*cte
  by(i,jc(j),kc(k))=by(i,jc(j),kc(k))+effekc*ste
else
 k2=fkxs(i)+fkys(j)
 kxy=sqrt(k2)
 k2=k2+fkzs(k)
 kk=sqrt(k2)
 c1=fkz(k)/kk
 s1=kxy/kk
 c2=fkx(i)/kxy
 s2=fky(j)/kxy
 bx(i,j,k)=bx(i,j,k)+effek*(cte*c2*c1-ste*s2)
 by(i,j,k)=by(i,j,k)+effek*(cte*s2*c1+ste*c2)
 bz(i,j,k)=bz(i,j,k)-effek*cte*s1
 if ((i.eq.1).or.(i.eq.NX2P1)) then
  effekc=conjg(effek)
  bx(i,jc(j),kc(k))=bx(i,jc(j),kc(k))+effekc*(cte*c2*c1-ste*s2)
  by(i,jc(j),kc(k))=by(i,jc(j),kc(k))+effekc*(cte*s2*c1+ste*c2)
  bz(i,jc(j),kc(k))=bz(i,jc(j),kc(k))-effekc*cte*s1
 end if
end if

end do
!$acc end parallel loop

!call incompr(bx,by,bz,cc)
!write(6,*)cc

return
end

! Subroutine to generate the random seed for the forcing

real function rann(irand)
!$acc routine seq  
implicit none
integer irand
integer mask

mask = '7FFFFFFF'X
irand = iand((69069*irand + 1),mask)
rann = dble(irand) / dble(mask)

return
end
