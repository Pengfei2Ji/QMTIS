!!****m* ABINIT/interfaces_62_poisson
!! NAME
!! interfaces_62_poisson
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/62_poisson
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

module interfaces_62_poisson

 implicit none

interface
 subroutine psolver_hartree(enhartr, hgrid, icoulomb, me, mpi_comm, nfft, ngfft, nproc,&  
  &  nscforder, nspden, rhor, vhartr, usewvl)
  use defs_basis
  implicit none
  integer, intent(in) :: icoulomb
  integer, intent(in) :: me
  integer, intent(in) :: mpi_comm
  integer, intent(in) :: nfft
  integer, intent(in) :: nproc
  integer, intent(in) :: nscforder
  integer, intent(in) :: nspden
  integer, intent(in) :: usewvl
  real(dp), intent(out) :: enhartr
  integer, intent(in) :: ngfft(3)
  real(dp),intent(in) :: hgrid(3)
  real(dp),intent(in) :: rhor(nfft,nspden)
  real(dp),intent(out) :: vhartr(nfft)
 end subroutine psolver_hartree
end interface

interface
 subroutine psolver_kernel(hgrid, iaction,  icoulomb,&  
  &  iproc, kernel, mpi_comm, ngfft, nproc, nscforder)
  use defs_basis
  use defs_wvltypes
  implicit none
  integer,intent(in) :: iaction
  integer,intent(in) :: icoulomb
  integer,intent(in) :: iproc
  integer,intent(in) :: mpi_comm
  integer,intent(in) :: nproc
  integer,intent(in) :: nscforder
  type(coulomb_operator),intent(inout) :: kernel
  integer, intent(in) :: ngfft(3)
  real(dp),intent(in) :: hgrid(3)
 end subroutine psolver_kernel
end interface

interface
 subroutine psolver_rhohxc(enhartr, enxc, envxc, icoulomb, ixc,&  
  &  mpi_enreg, nfft, ngfft, nhat,nhatdim,&  
  &  nscforder, nspden, n3xccc, rhor, rprimd,&  
  &  usexcnhat,usepaw,usewvl,vhartr, vxc, vxcavg, wvl,wvl_den,wvl_e,&  
  &  xccc3d,xclevel,xc_denpos)
  use defs_basis
  use defs_abitypes
  use defs_wvltypes
  implicit none
  integer, intent(in) :: icoulomb
  integer, intent(in) :: ixc
  integer, intent(in) :: n3xccc
  integer, intent(in) :: nfft
  integer, intent(in) :: nhatdim
  integer, intent(in) :: nscforder
  integer, intent(in) :: nspden
  integer,intent(in) :: usepaw
  integer, intent(in) :: usewvl
  integer,intent(in) :: usexcnhat
  integer,intent(in) :: xclevel
  real(dp), intent(out) :: enhartr
  real(dp), intent(out) :: envxc
  real(dp), intent(out) :: enxc
  type(mpi_type), intent(inout) :: mpi_enreg
  real(dp), intent(out) :: vxcavg
  type(wvl_internal_type), intent(in) :: wvl
  type(wvl_denspot_type), intent(inout) :: wvl_den
  type(wvl_energy_terms), intent(inout) :: wvl_e
  real(dp), intent(in) :: xc_denpos
  integer, intent(in) :: ngfft(18)
  real(dp),intent(in) :: nhat(nfft,nspden*nhatdim)
  real(dp),intent(inout) :: rhor(nfft, nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: vhartr(nfft)
  real(dp),intent(out) :: vxc(nfft, nspden)
  real(dp),intent(in) :: xccc3d(n3xccc)
 end subroutine psolver_rhohxc
end interface

end module interfaces_62_poisson
!!***
