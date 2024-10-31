!!!!!!!!!!!!!!!!!!!!!
! Calculate spectra !
!!!!!!!!!!!!!!!!!!!!!

subroutine spectrum(bx,by,bz,ifr)
use paran
use commk
use commphys
use cudafor
implicit none

complex*16, dimension(NX2P1,NY,NZ) :: bx,by,bz
real*8, dimension(NBIN) :: sp,sp2d
real*8 :: kk,scra,fac
integer :: i,j,k,ib,ifr
character*40 :: nome

! Initialize 
sp=0.0
sp2d=0.0

do k=1,NZ
do j=1,NY
do i=1,NX2P1
 fac=1.d0
 if ((i.eq.1).or.(i.eq.NX2P1)) fac=0.5d0
 kk=sqrt(fkxs(i)+fkys(j)+fkzs(k))
 ib=int(kk)+1
 if ((ib.gt.0).and.(ib.lt.NBIN)) then
  scra=cdabs(bx(i,j,k))**2+cdabs(by(i,j,k))**2+cdabs(bz(i,j,k))**2
  sp(ib)=sp(ib)+kk*kk*fac*scra
 end if
end do
end do
end do

do j=1,NY
do i=1,NX2P1
	 fac=1.d0
	 if ((i.eq.1).or.(i.eq.NX2P1)) fac=0.5d0
	 kk=sqrt(fkxs(i)+fkys(j)+fkzs(1))
	 ib=int(kk)+1
	 if ((ib.gt.0).and.(ib.lt.NBIN)) then
		scra=cdabs(bz(i,j,1))**2  ! E_{2D} = (v_x^2+v_y^2)|_{kz=0} =k_\perp^2 |b_z|^2
		sp2d(ib)=sp2d(ib)+kk*kk*fac*scra
	 endif
end do
end do

do ib=1,NBIN
 write(20,99)real(ib-1),real(sp(ib)),real(sp2d(ib))
 99 format(3g)
end do

return
end

!!!!!!!!!!!!!!!!!!!!!
! Calculate fluxes  !
!!!!!!!!!!!!!!!!!!!!!

subroutine fluxes(bx,by,bz,ax,ay,az,ifr)
use paran
use commk
use commphys
implicit none

complex*16, dimension(NX2P1,NY,NZ) :: bx,by,bz,ax,ay,az
real*8, dimension(NBIN) :: sp,spf
real*8 kk,scra,fac
integer i,j,k,ib,ifr
character*40 nome

! Initialize 
sp=0.d0

do k=1,NZ
do j=1,NY
do i=1,NX2P1
 fac=2.d0
 if ((i.eq.1).or.(i.eq.NX2P1)) fac=1.d0
  kk=sqrt(fkxs(i)+fkys(j)+fkzs(k))
  ib=int(kk)
  scra=real(bx(i,j,k)*conjg(ax(i,j,k))+&
      by(i,j,k)*conjg(ay(i,j,k))+bz(i,j,k)*conjg(az(i,j,k)))
  if (ib.lt.NBIN) then
    sp(ib)=sp(ib)+fac*kk*kk*scra
  end if
end do
end do
end do

do i=1,NBIN
	spf(i)=sum(sp(i:NBIN))
end do

do ib=1,NBIN
 write(30,99)real(ib),real(spf(ib))
 99 format(3g)
end do

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate dissipative terms !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! --------------------------------------------------
! Computes energy, enstrophy dissipation
! 
! vertical direction = z (x)
! --------------------------------------------------
subroutine glob(ax,ay,az,ifr)
use paran
use commk
use commphys
implicit none

complex*16, dimension(NX2P1,NY,NZ) :: ax,ay,az
real*8, dimension(9) :: temp
real*8 :: scrax,scray,scraz
real*8 :: k2,fac
integer i,j,k,ifr

call curl(ax,ay,az) ! b = velocity
!$acc update host(ax,ay,az)

temp=0.d0  ! use of temp to store energy global quantities

do k=1,NZ
do j=1,NY
do i=1,NX2P1
 fac=1.0
 if ((i.eq.1).or.(i.eq.NX2P1)) fac=0.5d0
 k2=fkxs(i)+fkys(j)+fkzs(k)
 if (k2.gt.0.d0) then
 scrax=fac*cdabs(ax(i,j,k))**2
 scray=fac*cdabs(ay(i,j,k))**2
 scraz=fac*cdabs(az(i,j,k))**2
 temp(1)=temp(1)+scrax                       ! Energy x component
 temp(2)=temp(2)+scray                       ! Energy y component
 temp(3)=temp(3)+scraz                       ! Energy z component
 temp(4)=temp(4)+2.d0*nu*scrax*k2**alpnu     ! Viscous energy diss x component
 temp(5)=temp(5)+2.d0*nu*scray*k2**alpnu     ! Viscous energy diss y component
 temp(6)=temp(6)+2.d0*nu*scraz*k2**alpnu     ! Viscous energy diss z component
 temp(7)=temp(7)+2.d0*mu*scrax*k2**alpmu     ! Friction energy diss x component
 temp(8)=temp(8)+2.d0*mu*scray*k2**alpmu     ! Friction energy diss y component
 temp(9)=temp(9)+2.d0*mu*scraz*k2**alpmu     ! Friction energy diss z component
 end if
end do
end do
end do

write(18,99)real(t),(real(temp(i)),i=1,9),inpe
call flush(18)

99 format(11g)

write(6,*)"t = ",t,", Eb = ",real(sum(temp(4:9))/inpe),", E = ",real(sum(temp(1:3)))
call flush(6)

! vorticity -------------------------------------------
call curl(ax,ay,az) ! b = vorticity
!$acc update host(ax,ay,az)

temp=0.d0  ! reuse of temp to store vorticity global quantities

do k=1,NZ
do j=1,NY
do i=1,NX2P1
 fac=1.0
 if ((i.eq.1).or.(i.eq.NX2P1)) fac=0.5d0
 k2=fkxs(i)+fkys(j)+fkzs(k)
 if (k2.gt.0.d0) then
 scrax=fac*cdabs(ax(i,j,k))**2
 scray=fac*cdabs(ay(i,j,k))**2
 scraz=fac*cdabs(az(i,j,k))**2
 temp(1)=temp(1)+scrax                       ! Enstrophy x component
 temp(2)=temp(2)+scray                       ! Enstrophy y component
 temp(3)=temp(3)+scraz                       ! Enstrophy z component
 temp(4)=temp(4)+2.d0*nu*scrax*k2**alpnu     ! Viscous enstrophy diss x component
 temp(5)=temp(5)+2.d0*nu*scray*k2**alpnu     ! Viscous enstrophy diss y component
 temp(6)=temp(6)+2.d0*nu*scraz*k2**alpnu     ! Viscous enstrophy diss z component
 temp(7)=temp(7)+2.d0*mu*scrax*k2**alpmu     ! Friction enstrophy diss x component
 temp(8)=temp(8)+2.d0*mu*scray*k2**alpmu     ! Friction enstrophy diss y component
 temp(9)=temp(9)+2.d0*mu*scraz*k2**alpmu     ! Friction enstrophy diss z component
 end if
end do
end do
end do

write(19,98)real(t),(real(temp(i)),i=1,9),inpz
call flush(19)

98 format(11g)

return
end