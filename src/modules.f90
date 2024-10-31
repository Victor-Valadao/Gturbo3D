! Sim parameters
module paran
 implicit none
 integer, parameter :: NX=1024, NY=1024, NZ=16
 integer, parameter :: NXP2=NX+2
 integer, parameter :: NX2=NX/2, NX2P1=NX2+1, NX2P2=NX2+2
 integer, parameter :: NY2=NY/2, NY2P1=NY2+1, NY2P2=NY2+2
 integer, parameter :: NZ2=NZ/2, NZ2P1=NZ2+1, NZ2P2=NZ2+2
 integer, parameter :: NBIN=NX/2
end module paran

! FFT parameters and plans
module stream
  use cufft
  use cudafor
  implicit none
  integer(kind=cuda_stream_kind) :: stream1
  integer(kind=cuda_stream_kind) :: stream2
end module stream

module plan
	implicit none
	integer plan_dir
	integer plan_inv
end module plan

! physical parameters
module commphys
	implicit none
	real*8 :: xlx,xly,xlz,dkx,dky,dkz
	real*8 :: nu,mu,alpnu,alpmu
	real*8 :: dt,dt2,sdt,t
	real*8 :: inpz,inpe,omega
end module commphys

! forcing parameters
module commforc
	implicit none
	integer, parameter :: NFORC = 2048
	real*8, dimension(NFORC) :: ff
	real*8 :: famp,k1f,k2f
	integer, dimension(NFORC) :: ikxf, ikyf, ikzf
	integer :: nf
end module commforc

! diagnostics parameters
module commsim
	implicit none
	integer :: npas,iout,imix,split
end module commsim

! wavenumber parameters
module commk
	use paran
	implicit none
	real*8, dimension(NX2P1,NY,NZ) :: emk
	real*8, dimension(NX2P1) :: emkx
	real*8, dimension(NY) :: emky
	real*8, dimension(NZ) :: emkz
  real*8, dimension(NX2P1) :: fkxs, fkx
	real*8, dimension(NY) :: fkys, fky
	real*8, dimension(NZ) :: fkzs, fkz
	integer, dimension(NY) :: jc
	integer, dimension(NZ) :: kc
	real*8 :: rtrunc,fkxsmax,fkysmax,fkzsmax
	integer :: truncate
end module commk