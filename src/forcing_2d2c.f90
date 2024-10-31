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
use paran
use commk
use commphys
use commforc
implicit none

complex*16, dimension(NX2P1,NY,NZ) :: bx,by,bz
complex*16 effe,effek,effekc
real*8, parameter :: pi2=2.d0*3.14159265358979d0
real*8 te,ste,cte
real*8 s1,s2,c1,c2 
real*8 k2,kk,kxy
integer i,j,k,ik,seed
real :: rann

!$acc parallel loop seq present(bx,by,bz,jc,kc,ff,ikxf,ikyf,ikzf,seed,sdt)
do ik=1,nf
	i=ikxf(ik)
	j=ikyf(ik)
	k=ikzf(ik)

	te=pi2*rann(seed)
	effe=sdt*famp*cdexp((0.d0,1.d0)*te)  
	effek=effe*ff(ik)
	 
	bz(i,j,k)=bz(i,j,k)+effek
	if ((i.eq.1).or.(i.eq.NX2P1)) then
		effekc=conjg(effek)
		bz(i,jc(j),kc(k))=bz(i,jc(j),kc(k))+effekc
	end if
end do
!$acc end parallel loop

! call incompr(bx,by,bz)

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

