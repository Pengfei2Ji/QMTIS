!!****m* ABINIT/interfaces_41_xc_lowlevel
!! NAME
!! interfaces_41_xc_lowlevel
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/41_xc_lowlevel
!!
!! COPYRIGHT
!! Copyright (C) 2010-2016 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!! THIS FILE IS GENERATED AUTOMATICALLY BY abilint.
!! To do that: config/scripts/abilint . .
!!
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module interfaces_41_xc_lowlevel

 implicit none

interface
 subroutine check_kxc(ixc,optdriver)
  implicit none
  integer, intent(in) :: ixc
  integer, intent(in) :: optdriver
 end subroutine check_kxc
end interface

interface
 subroutine drivexc(exc,ixc,npts,nspden,order,rho_updn,vxcrho,ndvxc,ngr2,nd2vxc,nvxcgrho,&  
  &  dvxc,d2vxc,grho2_updn,vxcgrho,el_temp,exexch,fxcT,&  
  &  lrho_updn,vxclrho,tau_updn,vxctau,xc_tb09_c)  !Optional arguments
  use defs_basis
  implicit none
  integer,intent(in),optional :: exexch
  integer,intent(in) :: ixc
  integer,intent(in) :: nd2vxc
  integer,intent(in) :: ndvxc
  integer,intent(in) :: ngr2
  integer,intent(in) :: npts
  integer,intent(in) :: nspden
  integer,intent(in) :: nvxcgrho
  integer,intent(in) :: order
  real(dp),intent(in),optional :: el_temp
  real(dp),intent(in),optional :: xc_tb09_c
  real(dp),intent(out),optional :: d2vxc(npts,nd2vxc)
  real(dp),intent(out),optional :: dvxc(npts,ndvxc)
  real(dp),intent(out) :: exc(npts)
  real(dp),intent(out),optional :: fxcT(:)
  real(dp),intent(in),optional :: grho2_updn(npts,ngr2)
  real(dp),intent(in),optional :: lrho_updn(npts,nspden)
  real(dp),intent(in) :: rho_updn(npts,nspden)
  real(dp),intent(in),optional :: tau_updn(npts,nspden)
  real(dp),intent(out),optional :: vxcgrho(npts,nvxcgrho)
  real(dp),intent(out),optional :: vxclrho(npts,nspden)
  real(dp),intent(out) :: vxcrho(npts,nspden)
  real(dp),intent(out),optional :: vxctau(npts,nspden)
 end subroutine drivexc
end interface

interface
 subroutine drivexc_main(exc,ixc,mgga,ndvxc,nd2vxc,ngr2,npts,nspden,nvxcgrho,order,rho,vxcrho,xclevel,&  
  &  dvxc,d2vxc,el_temp,exexch,fxcT,grho2,lrho,tau,vxcgrho,vxclrho,vxctau,xc_tb09_c) ! Optional arguments
  use defs_basis
  implicit none
  integer,intent(in),optional :: exexch
  integer,intent(in) :: ixc
  integer,intent(in) :: mgga
  integer,intent(in) :: nd2vxc
  integer,intent(in) :: ndvxc
  integer,intent(in) :: ngr2
  integer,intent(in) :: npts
  integer,intent(in) :: nspden
  integer,intent(in) :: nvxcgrho
  integer,intent(in) :: order
  integer,intent(in) :: xclevel
  real(dp),intent(in),optional :: el_temp
  real(dp),intent(in),optional :: xc_tb09_c
  real(dp),intent(out),optional :: d2vxc(npts,nd2vxc)
  real(dp),intent(out),optional :: dvxc(npts,ndvxc)
  real(dp),intent(out) :: exc(npts)
  real(dp),intent(out),optional :: fxcT(:)
  real(dp),intent(in),optional :: grho2(npts,ngr2)
  real(dp),intent(in),optional :: lrho(npts,nspden*mgga)
  real(dp),intent(in) :: rho(npts,nspden)
  real(dp),intent(in),optional :: tau(npts,nspden*mgga)
  real(dp),intent(out),optional :: vxcgrho(npts,nvxcgrho)
  real(dp),intent(out),optional :: vxclrho(npts,nspden*mgga)
  real(dp),intent(out) :: vxcrho(npts,nspden)
  real(dp),intent(out),optional :: vxctau(npts,nspden*mgga)
 end subroutine drivexc_main
end interface

interface
 subroutine echo_xc_name (ixc)
  implicit none
  integer, intent(in) :: ixc
 end subroutine echo_xc_name
end interface

interface
 subroutine invcb(rhoarr,rspts,npts)
  use defs_basis
  implicit none
  integer,intent(in) :: npts
  real(dp),intent(in) :: rhoarr(npts)
  real(dp),intent(out) :: rspts(npts)
 end subroutine invcb
end interface

interface
 subroutine mkdenpos(iwarn,nfft,nspden,option,rhonow,xc_denpos)
  use defs_basis
  implicit none
  integer,intent(inout) :: iwarn
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: option
  real(dp),intent(in) :: xc_denpos
  real(dp),intent(inout) :: rhonow(nfft,nspden)
 end subroutine mkdenpos
end interface

interface
 subroutine size_dvxc(ixc,ndvxc,ngr2,nd2vxc,nspden,nvxcdgr,order,add_tfw)
  implicit none
  integer, intent(in) :: ixc
  integer, intent(out) :: nd2vxc
  integer, intent(out) :: ndvxc
  integer, intent(out) :: ngr2
  integer, intent(in) :: nspden
  integer, intent(out) :: nvxcdgr
  integer, intent(in) :: order
  logical, intent(in), optional :: add_tfw
 end subroutine size_dvxc
end interface

interface
 subroutine xchcth(dvxcdgr,exci,grho2_updn,ixc,npts,nspden,order,rho_updn,vxci)
  use defs_basis
  implicit none
  integer,intent(in) :: ixc
  integer,intent(in) :: npts
  integer,intent(in) :: nspden
  integer,intent(in) :: order
  real(dp),intent(out) :: dvxcdgr(npts,2)
  real(dp),intent(out) :: exci(npts)
  real(dp),intent(in) :: grho2_updn(npts,2*nspden-1)
  real(dp),intent(in) :: rho_updn(npts,nspden)
  real(dp),intent(out) :: vxci(npts,nspden)
 end subroutine xchcth
end interface

interface
 subroutine xchelu(exc,npt,order,rspts,vxc,dvxc)  ! dvxc is optional
  use defs_basis
  implicit none
  integer,intent(in) :: npt
  integer,intent(in) :: order
  real(dp),intent(out),optional :: dvxc(npt)
  real(dp),intent(out) :: exc(npt)
  real(dp),intent(in) :: rspts(npt)
  real(dp),intent(out) :: vxc(npt)
 end subroutine xchelu
end interface

interface
 subroutine xciit(exc,fxc,npt,order,rspts,temp,vxc,&  
  &  dvxc)!Optional argument
  use defs_basis
  implicit none
  integer,intent(in) :: npt
  integer,intent(in) :: order
  real(dp),intent(in) :: temp
  real(dp),intent(out),optional :: dvxc(npt)
  real(dp),intent(out) :: exc(npt)
  real(dp),intent(out) :: fxc(npt)
  real(dp),intent(in) :: rspts(npt)
  real(dp),intent(out) :: vxc(npt)
 end subroutine xciit
end interface

interface
 subroutine xclb(grho2_updn,npts,nspden,rho_updn,vxci)
  use defs_basis
  implicit none
  integer,intent(in) :: npts
  integer,intent(in) :: nspden
  real(dp),intent(in) :: grho2_updn(npts,2*nspden-1)
  real(dp),intent(in) :: rho_updn(npts,nspden)
  real(dp),intent(inout) :: vxci(npts,nspden)
 end subroutine xclb
end interface

interface
 subroutine xcmult (depsxc,nfft,ngrad,nspden,nspgrad,rhonow)
  use defs_basis
  implicit none
  integer,intent(in) :: nfft
  integer,intent(in) :: ngrad
  integer,intent(in) :: nspden
  integer,intent(in) :: nspgrad
  real(dp),intent(in) :: depsxc(nfft,nspgrad)
  real(dp),intent(inout) :: rhonow(nfft,nspden,ngrad*ngrad)
 end subroutine xcmult
end interface

interface
 subroutine xcpbe(exci,npts,nspden,option,order,rho_updn,vxci,ndvxci,ngr2,nd2vxci,&  !Mandatory Arguments
  &  d2vxci,dvxcdgr,dvxci,exexch,grho2_updn)                          !Optional Arguments
  use defs_basis
  implicit none
  integer,intent(in),optional :: exexch
  integer,intent(in) :: nd2vxci
  integer,intent(in) :: ndvxci
  integer,intent(in) :: ngr2
  integer,intent(in) :: npts
  integer,intent(in) :: nspden
  integer,intent(in) :: option
  integer,intent(in) :: order
  real(dp),intent(out),optional :: d2vxci(npts,nd2vxci)
  real(dp),intent(out),optional :: dvxcdgr(npts,3)
  real(dp),intent(out),optional :: dvxci(npts,ndvxci)
  real(dp),intent(out) :: exci(npts)
  real(dp),intent(in),optional :: grho2_updn(npts,ngr2)
  real(dp),intent(in) :: rho_updn(npts,nspden)
  real(dp),intent(out) :: vxci(npts,nspden)
 end subroutine xcpbe
end interface

interface
 subroutine xcpositron(fnxc,grhoe2,ixcpositron,ngr,npt,posdensity0_limit,rhoer,rhopr,vxce,vxcegr,vxcp,&  
  &  dvxce,dvxcp) ! optional arguments
  use defs_basis
  implicit none
  integer,intent(in) :: ixcpositron
  integer,intent(in) :: ngr
  integer,intent(in) :: npt
  logical,intent(in) :: posdensity0_limit
  real(dp),intent(out),optional :: dvxce(npt)
  real(dp),intent(out),optional :: dvxcp(npt)
  real(dp),intent(out) :: fnxc(npt)
  real(dp),intent(in) :: grhoe2(ngr)
  real(dp),intent(in) :: rhoer(npt)
  real(dp),intent(in) :: rhopr(npt)
  real(dp),intent(out) :: vxce(npt)
  real(dp),intent(out) :: vxcegr(ngr)
  real(dp),intent(out) :: vxcp(npt)
 end subroutine xcpositron
end interface

interface
 subroutine xcpzca(exc,npt,order,rhor,rspts,vxc,&  !Mandatory arguments
  &  dvxc)                            !Optional arguments
  use defs_basis
  implicit none
  integer,intent(in) :: npt
  integer,intent(in) :: order
  real(dp),intent(out),optional :: dvxc(npt)
  real(dp),intent(out) :: exc(npt)
  real(dp),intent(in) :: rhor(npt)
  real(dp),intent(in) :: rspts(npt)
  real(dp),intent(out) :: vxc(npt)
 end subroutine xcpzca
end interface

interface
 subroutine xcspol(exc,npts,nspden,order,rspts,vxc,zeta,ndvxc,&  !Mandatory arguments
  &  dvxc)                            !Optional arguments
  use defs_basis
  implicit none
  integer,intent(in) :: ndvxc
  integer,intent(in) :: npts
  integer,intent(in) :: nspden
  integer,intent(in) :: order
  real(dp),intent(out),optional :: dvxc(npts,ndvxc)
  real(dp),intent(out) :: exc(npts)
  real(dp),intent(in) :: rspts(npts)
  real(dp),intent(out) :: vxc(npts,nspden)
  real(dp),intent(in) :: zeta(npts)
 end subroutine xcspol
end interface

interface
 subroutine xctetr(exc,npt,order,rhor,rspts,vxc,&  !Mandatory arguments
  &  d2vxc,dvxc)                    !Optional arguments
  use defs_basis
  implicit none
  integer,intent(in) :: npt
  integer,intent(in) :: order
  real(dp),intent(out),optional :: d2vxc(npt)
  real(dp),intent(out),optional :: dvxc(npt)
  real(dp),intent(out) :: exc(npt)
  real(dp),intent(in) :: rhor(npt)
  real(dp),intent(in) :: rspts(npt)
  real(dp),intent(out) :: vxc(npt)
 end subroutine xctetr
end interface

interface
 subroutine xctfw(temp,exci,fxci,usefxc,rho_updn,vxci,npts,nspden,dvxcdgr,ndvxcdgr,grho2_updn,ngr2)
  use defs_basis
  implicit none
  integer,intent(in) :: ndvxcdgr
  integer,intent(in) :: ngr2
  integer,intent(in) :: npts
  integer,intent(in) :: nspden
  integer,intent(in) :: usefxc
  real(dp),intent(in) :: temp
  real(dp),intent(inout) :: dvxcdgr(npts,ndvxcdgr)
  real(dp),intent(inout) :: exci(npts)
  real(dp),intent(inout) :: fxci(npts*usefxc)
  real(dp),intent(in) :: grho2_updn(npts,ngr2)
  real(dp),intent(in) :: rho_updn(npts,nspden)
  real(dp),intent(inout) :: vxci(npts,nspden)
 end subroutine xctfw
end interface

interface
 subroutine xcwign(exc,npt,order,rspts,vxc,&  !Mandatory arguments
  &  dvxc)                           !Optional arguments
  use defs_basis
  implicit none
  integer,intent(in) :: npt
  integer,intent(in) :: order
  real(dp),intent(out),optional :: dvxc(npt)
  real(dp),intent(out) :: exc(npt)
  real(dp),intent(in) :: rspts(npt)
  real(dp),intent(out) :: vxc(npt)
 end subroutine xcwign
end interface

interface
 subroutine xcxalp(exc,npt,order,rspts,vxc, dvxc)  ! dvxc is optional
  use defs_basis
  implicit none
  integer,intent(in) :: npt
  integer,intent(in) :: order
  real(dp),intent(out),optional :: dvxc(npt)
  real(dp),intent(out) :: exc(npt)
  real(dp),intent(in) :: rspts(npt)
  real(dp),intent(out) :: vxc(npt)
 end subroutine xcxalp
end interface

end module interfaces_41_xc_lowlevel
!!***
