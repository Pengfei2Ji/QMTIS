!!****m* ABINIT/interfaces_62_wvl_wfs
!! NAME
!! interfaces_62_wvl_wfs
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/62_wvl_wfs
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

module interfaces_62_wvl_wfs

 implicit none

interface
 subroutine wvl_hpsitopsi(cprj,dtset, energies, istep, mcprj,mpi_enreg,&  
  &  residm, wvl,xcart)
  use defs_basis
  use m_energies
  use defs_abitypes
  use m_pawcprj
  use defs_wvltypes
  implicit none
  integer, intent(in) :: istep
  integer, intent(in) :: mcprj
  type(dataset_type), intent(in) :: dtset
  type(energies_type), intent(inout) :: energies
  type(mpi_type), intent(in) :: mpi_enreg
  real(dp), intent(inout) :: residm
  type(wvl_data), intent(inout) :: wvl
  type(pawcprj_type),dimension(dtset%natom,mcprj),intent(inout) :: cprj
  real(dp), intent(in) :: xcart(3, dtset%natom)
 end subroutine wvl_hpsitopsi
end interface

interface
 subroutine wvl_nl_gradient(grnl, mpi_enreg, natom, rprimd, wvl, xcart)
  use defs_basis
  use defs_abitypes
  use defs_wvltypes
  implicit none
  integer, intent(in) :: natom
  type(mpi_type),intent(inout) :: mpi_enreg
  type(wvl_data),intent(inout) :: wvl
  real(dp),intent(inout) :: grnl(3,natom)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: xcart(3,natom)
 end subroutine wvl_nl_gradient
end interface

interface
 subroutine wvl_psitohpsi(alphamix,eexctX, eexcu, ehart, ekin_sum, epot_sum, eproj_sum, eSIC_DC,&  
  &  itrp, iter, iscf, me, natom, nfft, nproc, nspden, rpnrm, scf,&  
  &  vexcu, wvl, wvlbigdft, xcart, xcstr,vtrial,vxc)
  use defs_basis
  use defs_wvltypes
  implicit none
  integer, intent(in) :: iscf
  integer, intent(in) :: iter
  integer, intent(in) :: itrp
  integer, intent(in) :: me
  integer, intent(in) :: natom
  integer, intent(in) :: nfft
  integer, intent(in) :: nproc
  integer, intent(in) :: nspden
  real(dp), intent(in) :: alphamix
  real(dp), intent(inout) :: eSIC_DC
  real(dp), intent(inout) :: eexctX
  real(dp), intent(inout) :: eexcu
  real(dp), intent(inout) :: ehart
  real(dp), intent(inout) :: ekin_sum
  real(dp), intent(inout) :: epot_sum
  real(dp), intent(inout) :: eproj_sum
  real(dp), intent(out) :: rpnrm
  logical, intent(in) :: scf
  real(dp), intent(inout) :: vexcu
  type(wvl_data), intent(inout) :: wvl
  logical, intent(in) :: wvlbigdft
  real(dp),intent(out), optional :: vtrial(nfft,nspden)
  real(dp),intent(out), optional :: vxc(nfft,nspden)
  real(dp), intent(inout) :: xcart(3, natom)
  real(dp), dimension(6), intent(out) :: xcstr
 end subroutine wvl_psitohpsi
end interface

interface
 subroutine wvl_read(dtset, hdr0, hdr, mpi_enreg, option, rprimd, wff, wfs, wvl, xred)
  use defs_basis
  use defs_abitypes
  use m_wffile
  use defs_wvltypes
  implicit none
  integer, intent(in) :: option
  type(dataset_type), intent(in) :: dtset
  type(hdr_type), intent(in) :: hdr
  type(hdr_type), intent(in) :: hdr0
  type(mpi_type), intent(in) :: mpi_enreg
  type(wffile_type),intent(in) :: wff
  type(wvl_wf_type), intent(inout) :: wfs
  type(wvl_internal_type), intent(in) :: wvl
  real(dp), intent(in) :: rprimd(3, 3)
  real(dp), intent(in) :: xred(3, dtset%natom)
 end subroutine wvl_read
end interface

interface
 subroutine wvl_write(dtset, eigen, mpi_enreg, option, rprimd, wff, wfs, wvl, xred)
  use defs_basis
  use defs_abitypes
  use m_wffile
  use defs_wvltypes
  implicit none
  integer, intent(in) :: option
  type(dataset_type), intent(in) :: dtset
  type(mpi_type), intent(in) :: mpi_enreg
  type(wffile_type),intent(in) :: wff
  type(wvl_wf_type), intent(in) :: wfs
  type(wvl_internal_type), intent(in) :: wvl
  real(dp), intent(in), target :: eigen(dtset%mband)
  real(dp), intent(in) :: rprimd(3, 3)
  real(dp), intent(in) :: xred(3, dtset%natom)
 end subroutine wvl_write
end interface

interface
 subroutine wvl_tail_corrections(dtset, energies, etotal, mpi_enreg, psps, wvl, xcart)
  use defs_basis
  use m_energies
  use defs_abitypes
  use defs_datatypes
  use defs_wvltypes
  implicit none
  type(dataset_type),intent(in) :: dtset
  type(energies_type),intent(inout) :: energies
  real(dp),intent(out) :: etotal
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  type(wvl_data),intent(inout) :: wvl
  real(dp),intent(in) :: xcart(3,dtset%natom)
 end subroutine wvl_tail_corrections
end interface

end module interfaces_62_wvl_wfs
!!***
