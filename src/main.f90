use cudafor
use openacc
use cufft
use iso_c_binding , only: C_PTR
use paran
use stream
use plan
use commsim
use commphys
use commk
use commforc
implicit none

complex*16, dimension(NX2P1,NY,NZ) :: bx,by,bz,bnlx,bnly,bnlz,ux,uy,uz
integer :: ipas,ifr,jfr,istat,seed,ierr,err_dir,err_inv
integer :: startTime, stopTime, trate
integer(kind=cuda_count_kind) :: fre,tota
real :: s
character*40 :: nome

call system_clock(startTime)

! Shows the amount of memory consumed in the GPU
call show_mem(fre,tota,0)

! InvFFT Complex to Real double precision
err_inv = err_inv + cufftPlan3D(plan_inv,NZ,NY,NX,CUFFT_Z2D)
istat = cudaStreamCreate(stream1)
err_inv = err_inv + cufftSetStream(plan_inv,stream1)

! DirFFT Real to Complex double precision
err_dir = err_dir + cufftPlan3D(plan_dir,NZ,NY,NX,CUFFT_D2Z)
istat = cudaStreamCreate(stream2)
err_dir = err_dir + cufftSetStream(plan_dir,stream2)

! Reads if the simulation is a continuation of a previous one
open(unit=1,file='./curframe.dat')
 read(1,*)ifr,t
close(1)

! Reads the random seed
write(nome,"('./files/seed.',i3.3)")ifr
open(unit=2,file=nome)
 read(2,*)seed
close(2)

! Initializes wave number, forcing, initial field z (real space)
call inieuler(bx,by,bz,ifr)

! Opening files to save global quantities and useful information
write(nome,"('./files/globalE.',i3.3)")ifr
Open(unit=18,file=nome,form='formatted',status='unknown',access='append')

write(nome,"('./files/globalZ.',i3.3)")ifr
Open(unit=19,file=nome,form='formatted',status='unknown',access='append')

write(nome,"('./files/spectra.',i3.3)")ifr
Open(unit=20,file=nome,form='formatted',status='unknown',access='append')

write(nome,"('./files/fluxes.',i3.3)")ifr
Open(unit=30,file=nome,form='formatted',status='unknown',access='append')


! Shows the amount of memory consumed in the GPU
call show_mem(fre,tota,1)

! ========================================================

! Copying data to the GPU data

if(split.eq.1) then
	!$acc enter data copyin(emkx,emky,emkz)     !<-- calc in defk called in inieuler
	else
	!$acc enter data copyin(emk)                !<-- calc in defk called in inieuler
end if

!$acc enter data copyin(ikxf,ikyf,ikzf,ff)    !<-- calc in iniforcing called in inieuler
!$acc enter data copyin(fkxs,fkys,fkzs)       !<-- calc in defk called in inieuler
!$acc enter data copyin(fkx,fky,fkz,jc,kc)    !<-- calc in defk called in inieuler
!$acc enter data copyin(dt,dt2,sdt)           !<-- calc in inieuler/loaded
!$acc enter data copyin(seed)                 !<-- read in main used in forcing/rann
!$acc enter data copyin(bnlx,bnly,bnlz)       !<-- allocate memory for bnl
!$acc enter data copyin(ux,uy,uz)             !<-- allocate memory for u

!$acc enter data copyin(bx,by,bz)             !<-- pass the initial field to the GPU

! Shows the amount of memory consumed in the GPU
call show_mem(fre,tota,2)

! DirFFT the initial field + b=velocity=>potential vorticity + dealiasing
call fft_dir(bx, plan_dir)
call fft_dir(by, plan_dir)
call fft_dir(bz, plan_dir)

! Remove possible compressible part, useful if init_noise on
call incompr(bx,by,bz)

! Pass velocity -> potential vector
call invcurl(bx,by,bz)     

! De-aliasing
if (truncate.eq.1) then
 call trunc(bx)
 call trunc(by)
 call trunc(bz)
end if

! Computation of initial global quantities
! $acc kernels present(bx,by,bz,bnlx,bnly,bnlz)
! bnlx=bx
! bnly=by  
! bnlz=bz
! $acc end kernels
! $acc update host(bnlx,bnly,bnlz)
! call glob(bnlx,bnly,bnlz,0)

jfr=1
! =============================================================

!  Time evolution <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
do ipas=1,npas	
	
	! Put B on BNL
	!$acc kernels present(bx,by,bz,bnlx,bnly,bnlz)
	bnlx=bx
	bnly=by  
	bnlz=bz
	!$acc end kernels
	
	! Calculate the nonlinear part and put on BNL
	call nlt(bnlx,bnly,bnlz,ux,uy,uz,bx,by,bz)

	! Write Spectrum and Flux in the time t (not t+dt)
	if (mod(ipas,imix).eq.0) then
	 !$acc update host(bx,by,bz,bnlx,bnly,bnlz)
	 call fluxes(bx,by,bz,bnlx,bnly,bnlz,ipas)
	 call spectrum(bx,by,bz,ipas)
	end if
	
	! Half RK-step
	call step1(bx,by,bz,bnlx,bnly,bnlz)
	
	! Calculate nonlinear part in half step put on BNL
	call nlt(bnlx,bnly,bnlz,ux,uy,uz,bx,by,bz)
	! Last RK-step
	call step2(bx,by,bz,bnlx,bnly,bnlz)
	
	! Add forcing contribution
	call forcing(bx,by,bz,seed)
	t=t+dt
	
	! Global quantities and save field in t+dt
	if (mod(ipas,imix).eq.0) then
	 !$acc kernels present(bx,by,bz,bnlx,bnly,bnlz)
	 bnlx=bx
	 bnly=by  
	 bnlz=bz
	 !$acc end kernels
	 !$acc update host(bnlx,bnly,bnlz)
	 call glob(bnlx,bnly,bnlz,ipas)  ! bnl=b used to calc glob because inside  
	end if
	
  if ((mod(ipas,iout).eq.0).and.(ipas.ne.npas)) then  ! Save velocity each iout
	  !$acc kernels present(bx,by,bz,bnlx,bnly,bnlz)
	  bnlx=bx
	  bnly=by  
	  bnlz=bz
	  !$acc end kernels
	  
	  ! Calculate velocity=curl of potential vector
	  call curl(bnlx,bnly,bnlz)
		! FFt to save in the physical space
		call fft_inv(bnlx, plan_inv)
		call fft_inv(bnly, plan_inv)
		call fft_inv(bnlz, plan_inv)
		!$acc update host(bnlx,bnly,bnlz,seed) 
  	
  	call writef2(bnlx,bnly,bnlz,ifr,jfr)        ! Save
  	jfr=jfr+1
  	
  	write(6,*)"Saved v field, seed and curframe at timestep = ",ipas
  	call flush(6)
  	
  end if
end do 

! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
ifr=ifr+1

!$acc kernels present(bx,by,bz,bnlx,bnly,bnlz)
bnlx=bx
bnly=by  
bnlz=bz
!$acc end kernels

! Calculate v=curl(b)
call curl(bnlx,bnly,bnlz)

call fft_inv(bnlx, plan_inv)
call fft_inv(bnly, plan_inv)
call fft_inv(bnlz, plan_inv)  ! FFt to save in the physical space
!$acc update host(bnlx,bnly,bnlz,seed) 

call writef(bnlx,bnly,bnlz,ifr)        ! Save

! Export seed for continuation to seed.****
write(nome,"('./files/seed.',i3.3)")ifr
open(unit=2,file=nome)
write(2,*)seed
close(2)

! Export data to curframe.dat to continue, if you want
open(unit=1,file='./curframe.dat')
 write(1,*)ifr,t
close(1)
  
close(18)
close(19)
close(20)
close(30)

! Exit GPU data Freeing the memory

if(split.eq.1) then
	!$acc exit data delete(emkx,emky,emkz) finalize
	else
	!$acc exit data delete(emk) finalize
end if

!$acc exit data delete(ikxf,ikyf,ikzf,ff) finalize
!$acc exit data delete(fkxs,fkys,fkzs) finalize
!$acc exit data delete(fkx,fky,fkz,jc,kc) finalize
!$acc exit data delete(dt,dt2,sdt) finalize
!$acc exit data delete(seed) finalize
!$acc exit data delete(bnlx,bnly,bnlz) finalize
!$acc exit data delete(ux,uy,uz) finalize

!$acc exit data delete(bx,by,bz) finalize

err_inv = err_inv + cufftDestroy(plan_inv)
err_dir = err_dir + cufftDestroy(plan_dir)
istat = cudaStreamDestroy(stream1)
istat = cudaStreamDestroy(stream2)

call system_clock(stopTime,trate)
s=(real(stopTime-startTime)/(real(npas)*real(trate,8)))
write(6,*),"Mean time per iteration      = ", s
call flush(6)

end