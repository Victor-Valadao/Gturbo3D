subroutine inieuler(bx,by,bz,ifr)
use paran
use plan
use commsim
use commforc
use commphys
use commk
implicit none

complex*16, dimension(NX2P1,NY,NZ) :: bx,by,bz
integer ifr,status,seed
real :: rann

open(unit=9,file='./params.dat')
read(9,*)xlx,xly,xlz
read(9,*)nu,alpnu,mu,alpmu
read(9,*)dt,omega
read(9,*)npas,iout,imix
read(9,*)rtrunc,truncate
read(9,*)famp,k1f,k2f

split=0

if (alpnu < 1.0005) then 
if (alpnu > 0.9995) then 
	if (alpmu < 0.0005) then 
	if (alpmu >-0.0005) then 
		split=1
	end if
	end if
end if
end if

split=0

write(6,*)' 3D-rotating Navier-Stokes simulation'
write(6,*)' Resolution       = ',NX,NY,NZ
write(6,*)' Rotation rate    = ',real(omega)
write(6,*)' Viscosity        = ',real(nu)
write(6,*)' Viscosity order  = ',real(alpnu)
write(6,*)' Friction         = ',real(mu)
write(6,*)' Friction order   = ',real(alpmu)
write(6,*)' '

if (truncate.eq.1) then
 write(6,*)' Truncation at   = ',real(rtrunc)
else
 write(6,*)' No truncation'
end if

if (split.eq.1) then
 write(6,*)' Splitted linear integrator (only for Ekman-NSE)'
else
 write(6,*)' Unsplitted linear integrator'
end if

write(6,*)' '
write(6,*)' Time step           = ',real(dt)
write(6,*)' Starting Frame      = ',ifr
write(6,*)' Total step number   = ',npas
write(6,*)' Outputs Frames/Diag = ',iout,imix
write(6,*)' ------------------------------------------------'
write(6,*)' '
  
call flush(6)

sdt=sqrt(dt)

! Initialize wavenumber, see defk.f90
call defk

! Initialize forcing, see iniforcing.f90
call iniforcing

! Check previous simulations and load b=velocity
  if (ifr.gt.0) then
    call readf(bx,by,bz,ifr)
  else
  	! Otherwise, start from zero
    bx=0.d0
    by=0.d0
    bz=0.d0
		! Adds random noise, useful for 3d2d simulations
		call init_noise(bx,by,bz,seed)
  endif

return
end
! ----------------------------------------------------------------	
