!!****m* ABINIT/interfaces_56_xc
!! NAME
!! interfaces_56_xc
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/56_xc
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

module interfaces_56_xc

 implicit none

interface
 subroutine calc_smeared_density(rhor,kappa_strategy,rhotilder,nfftf,ngfftf,npw,&  
  &  gvec,gprimd,ucvol,mpi_enreg,paral_kgb,kappa_in)
  use defs_basis
  use defs_abitypes
  implicit none
  integer, intent(in) :: kappa_strategy
  integer, intent(in) :: nfftf
  integer, intent(in) :: npw
  integer, intent(in) :: paral_kgb
  real(dp), intent(in), optional :: kappa_in
  type(mpi_type), intent(inout) :: mpi_enreg
  real(dp), intent(in) :: ucvol
  integer,intent(in) :: ngfftf(18)
  real(dp), intent(in) :: gprimd(3,3)
  integer,intent(in) :: gvec(3,npw)
  real(dp), intent(inout) :: rhor(nfftf)
  real(dp), intent(out) :: rhotilder(nfftf)
 end subroutine calc_smeared_density
end interface

interface
 subroutine dfpt_mkvxc(cplex,ixc,kxc,mpi_enreg,nfft,ngfft,nhat1,nhat1dim,nhat1gr,nhat1grdim,&  
  &  nkxc,nspden,n3xccc,option,paral_kgb,qphon,rhor1,rprimd,usexcnhat,vxc1,xccc3d1)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: ixc
  integer,intent(in) :: n3xccc
  integer,intent(in) :: nfft
  integer,intent(in) :: nhat1dim
  integer,intent(in) :: nhat1grdim
  integer,intent(in) :: nkxc
  integer,intent(in) :: nspden
  integer,intent(in) :: option
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: usexcnhat
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: kxc(nfft,nkxc)
  real(dp),intent(in) :: nhat1(cplex*nfft,nspden*nhat1dim)
  real(dp),intent(in) :: nhat1gr(cplex*nfft,nspden,3*nhat1grdim)
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(in),target :: rhor1(cplex*nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: vxc1(cplex*nfft,nspden)
  real(dp),intent(in) :: xccc3d1(cplex*n3xccc)
 end subroutine dfpt_mkvxc
end interface

interface
 subroutine dfpt_mkvxc_noncoll(cplex,ixc,kxc,mpi_enreg,nfft,ngfft,nhat1,nhat1dim,nhat1gr,nhat1grdim,&  
  &  nkxc,nkxc_cur,nspden,n3xccc,optnc,option,optxc,paral_kgb,qphon,rhor,rhor1,rprimd,usexcnhat,vxc1,xccc3d1)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: ixc
  integer,intent(in) :: n3xccc
  integer,intent(in) :: nfft
  integer,intent(in) :: nhat1dim
  integer,intent(in) :: nhat1grdim
  integer,intent(in) :: nkxc
  integer,intent(in) :: nkxc_cur
  integer,intent(in) :: nspden
  integer,intent(in) :: option
  integer,intent(in) :: optnc
  integer,intent(in) :: optxc
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: usexcnhat
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: kxc(nfft,nkxc)
  real(dp),intent(in) :: nhat1(cplex*nfft,nspden*nhat1dim)
  real(dp),intent(in) :: nhat1gr(cplex*nfft,nspden,3*nhat1grdim)
  real(dp),intent(in),target :: qphon(3)
  real(dp),intent(in) :: rhor(cplex*nfft,nspden)
  real(dp),intent(in),target :: rhor1(cplex*nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: vxc1(cplex*nfft,nspden)
  real(dp),intent(in) :: xccc3d1(cplex*n3xccc)
 end subroutine dfpt_mkvxc_noncoll
end interface

interface
 subroutine dfpt_mkvxcgga(cplex,gprimd,kxc,mpi_enreg,nfft,ngfft,&  
  &  nhat1,nhat1dim,nhat1gr,nhat1grdim,nkxc,&  
  &  nspden,paral_kgb,qphon,rhor1tmp,usexcnhat,vxc1)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: nfft
  integer,intent(in) :: nhat1dim
  integer,intent(in) :: nhat1grdim
  integer,intent(in) :: nkxc
  integer,intent(in) :: nspden
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: usexcnhat
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: kxc(nfft,nkxc)
  real(dp),intent(in) :: nhat1(cplex*nfft,2*nhat1dim)
  real(dp),intent(in) :: nhat1gr(cplex*nfft,2,3*nhat1grdim)
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(in),target :: rhor1tmp(cplex*nfft,2)
  real(dp),intent(out) :: vxc1(cplex*nfft,nspden)
 end subroutine dfpt_mkvxcgga
end interface

interface
 subroutine gammapositron(gamma,grhocore2,grhoe2,igamma,ngr,npt,rhocore,rhoer,rhopr,usecore)
  use defs_basis
  implicit none
  integer,intent(in) :: igamma
  integer,intent(in) :: ngr
  integer,intent(in) :: npt
  integer,intent(in) :: usecore
  real(dp),intent(out) :: gamma(npt,2)
  real(dp),intent(in) :: grhocore2(ngr*usecore)
  real(dp),intent(in) :: grhoe2(ngr)
  real(dp),intent(in) :: rhocore(npt*usecore)
  real(dp),intent(in) :: rhoer(npt)
  real(dp),intent(in) :: rhopr(npt)
 end subroutine gammapositron
end interface

interface
 subroutine gammapositron_fft(electronpositron,gamma,gprimd,igamma,mpi_enreg,&  
  &  n3xccc,nfft,ngfft,rhor_e,rhor_p,xccc3d)
  use defs_basis
  use defs_abitypes
  use m_electronpositron
  implicit none
  integer,intent(in) :: igamma
  integer,intent(in) :: n3xccc
  integer,intent(in) :: nfft
  type(electronpositron_type),pointer :: electronpositron
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(out) :: gamma(nfft,2)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: rhor_e(nfft)
  real(dp),intent(in) :: rhor_p(nfft)
  real(dp),intent(in) :: xccc3d(n3xccc)
 end subroutine gammapositron_fft
end interface

interface
 subroutine hartre(cplex,gmet,gsqcut,izero,mpi_enreg,nfft,ngfft,paral_kgb,qphon,rhog,vhartr,&  
  &  divgq0) ! Optional argument
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: izero
  integer,intent(in) :: nfft
  integer,intent(in) :: paral_kgb
  real(dp),intent(in),optional :: divgq0
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(in) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(out) :: vhartr(cplex*nfft)
 end subroutine hartre
end interface

interface
 subroutine hybrid_corr(dtset,ixc,nkxc,mpi_enreg,nfft,ngfft,nspden,rhor,rprimd,hybrid_mixing,vxc,enxc)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: ixc
  integer,intent(in) :: nfft
  integer,intent(in) :: nkxc
  integer,intent(in) :: nspden
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(inout) :: enxc
  real(dp),intent(in) :: hybrid_mixing
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: rhor(nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(inout) :: vxc(nfft,nspden)
 end subroutine hybrid_corr
end interface

interface
 subroutine mkcore(corstr,dyfrx2,grxc,mpi_enreg,natom,nfft,nspden,ntypat,n1,n1xccc,&  
  &  n2,n3,option,rprimd,typat,ucvol,vxc,xcccrc,xccc1d,xccc3d,xred)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: n1
  integer,intent(in) :: n1xccc
  integer,intent(in) :: n2
  integer,intent(in) :: n3
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  integer,intent(in) :: option
  type(mpi_type),intent(in) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  real(dp),intent(out) :: corstr(6)
  real(dp),intent(out) :: dyfrx2(3,3,natom)
  real(dp),intent(inout) :: grxc(3,natom)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: vxc(nfft,nspden)
  real(dp),intent(in) :: xccc1d(n1xccc,6,ntypat)
  real(dp),intent(inout) :: xccc3d(nfft)
  real(dp),intent(in) :: xcccrc(ntypat)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine mkcore
end interface

interface
 subroutine phase(ngfft,ph)
  use defs_basis
  implicit none
  integer,intent(in) :: ngfft
  real(dp),intent(out) :: ph(2*ngfft)
 end subroutine phase
end interface

interface
 subroutine rhohxc(dtset,enxc,gsqcut,izero,kxc,mpi_enreg,nfft,ngfft,&  
  &  nhat,nhatdim,nhatgr,nhatgrdim,nkxc,nk3xc,nspden,n3xccc,option,&  
  &  rhog,rhor,rprimd,strsxc,usexcnhat,vhartr,vxc,vxcavg,xccc3d,&  
  &  k3xc,electronpositron,taug,taur,vxctau,exc_vdw_out,add_tfw) ! optional arguments
  use defs_basis
  use defs_abitypes
  use m_electronpositron
  implicit none
  integer,intent(in) :: izero
  integer,intent(in) :: n3xccc
  integer,intent(in) :: nfft
  integer,intent(in) :: nhatdim
  integer,intent(in) :: nhatgrdim
  integer,intent(in) :: nk3xc
  integer,intent(in) :: nkxc
  integer,intent(in) :: nspden
  integer,intent(in) :: option
  integer,intent(in) :: usexcnhat
  logical,intent(in),optional :: add_tfw
  type(dataset_type),intent(in) :: dtset
  type(electronpositron_type),pointer,optional :: electronpositron
  real(dp),intent(out) :: enxc
  real(dp),intent(out),optional :: exc_vdw_out
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(in) :: mpi_enreg
  real(dp),intent(out) :: vxcavg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(out),optional :: k3xc(1:nfft,1:nk3xc)
  real(dp),intent(out) :: kxc(nfft,nkxc)
  real(dp),intent(in) :: nhat(nfft,nspden*nhatdim)
  real(dp),intent(in) :: nhatgr(nfft,nspden,3*nhatgrdim)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(in),target :: rhor(nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: strsxc(6)
  real(dp),intent(in),optional :: taug(:,:)
  real(dp),intent(in),optional :: taur(:,:)
  real(dp),intent(out) :: vhartr(nfft)
  real(dp),intent(out) :: vxc(nfft,nspden)
  real(dp),intent(out),optional :: vxctau(:,:,:)
  real(dp),intent(in) :: xccc3d(n3xccc)
 end subroutine rhohxc
end interface

interface
 subroutine rhohxcpositron(electronpositron,gprimd,kxcapn,mpi_enreg,nfft,ngfft,nhat,nkxc,nspden,n3xccc,&  
  &  paral_kgb,rhor,strsxc,ucvol,usexcnhat,usepaw,vhartr,vxcapn,vxcavg,xccc3d,xc_denpos)
  use defs_basis
  use defs_abitypes
  use m_electronpositron
  implicit none
  integer,intent(in) :: n3xccc
  integer,intent(in) :: nfft
  integer,intent(in) :: nkxc
  integer,intent(in) :: nspden
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: usepaw
  integer,intent(in) :: usexcnhat
  type(electronpositron_type),pointer :: electronpositron
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  real(dp),intent(out) :: vxcavg
  real(dp),intent(in) :: xc_denpos
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(out) :: kxcapn(nfft,nkxc)
  real(dp),intent(in) :: nhat(nfft,nspden*usepaw)
  real(dp),intent(in) :: rhor(nfft,nspden)
  real(dp),intent(out) :: strsxc(6)
  real(dp),intent(out) :: vhartr(nfft)
  real(dp),intent(out) :: vxcapn(nfft,nspden)
  real(dp),intent(in) :: xccc3d(n3xccc)
 end subroutine rhohxcpositron
end interface

interface
 subroutine xcden (cplex,gprimd,ishift,mpi_enreg,nfft,ngfft,ngrad,nspden,paral_kgb,qphon,rhor,rhonow,&  !Mandatory arguments
  &  lrhonow)              !Optional arguments
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: ishift
  integer,intent(in) :: nfft
  integer,intent(in) :: ngrad
  integer,intent(in) :: nspden
  integer,intent(in) :: paral_kgb
  type(mpi_type),intent(in) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(out),optional :: lrhonow(cplex*nfft,nspden)
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(out) :: rhonow(cplex*nfft,nspden,ngrad*ngrad)
  real(dp),intent(in) :: rhor(cplex*nfft,nspden)
 end subroutine xcden
end interface

interface
 subroutine xchybrid_ncpp_cc(dtset,enxc,mpi_enreg,nfft,ngfft,n3xccc,rhor,rprimd,strsxc,vxc,vxcavg,xccc3d)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: n3xccc
  integer,intent(in) :: nfft
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(out) :: enxc
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(out) :: vxcavg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: rhor(nfft,dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: strsxc(6)
  real(dp),intent(out) :: vxc(nfft,dtset%nspden)
  real(dp),intent(in) :: xccc3d(n3xccc)
 end subroutine xchybrid_ncpp_cc
end interface

interface
 subroutine xcpot (cplex,depsxc,gprimd,ishift,mgga,mpi_enreg,nfft,ngfft,ngrad,nspden,&  
  &  nspgrad,paral_kgb,qphon,rhonow,vxc,&  
  &  vxctau) ! optional argument
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: ishift
  integer,intent(in) :: mgga
  integer,intent(in) :: nfft
  integer,intent(in) :: ngrad
  integer,intent(in) :: nspden
  integer,intent(in) :: nspgrad
  integer,intent(in) :: paral_kgb
  type(mpi_type),intent(in) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: depsxc(cplex*nfft,nspgrad)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(in) :: rhonow(cplex*nfft,nspden,ngrad*ngrad)
  real(dp),intent(inout) :: vxc(cplex*nfft,nspden)
  real(dp),intent(inout),optional :: vxctau(cplex*nfft,nspden,4)
 end subroutine xcpot
end interface

end module interfaces_56_xc
!!***
