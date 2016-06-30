!!****m* ABINIT/interfaces_43_wvl_wrappers
!! NAME
!! interfaces_43_wvl_wrappers
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/43_wvl_wrappers
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

module interfaces_43_wvl_wrappers

 implicit none

interface
 subroutine nullify_wvl_data(wvl)
  use defs_wvltypes
  implicit none
  type(wvl_data) , intent(inout) :: wvl
 end subroutine nullify_wvl_data
end interface

interface
 subroutine paw2wvl(pawtab,proj,wvl)
  use m_pawtab
  use defs_wvltypes
  implicit none
  type(wvl_projectors_type),intent(inout) :: proj
  type(wvl_internal_type), intent(inout) :: wvl
  type(pawtab_type),intent(in) :: pawtab(:)
 end subroutine paw2wvl
end interface

interface
 subroutine paw2wvl_ij(option,paw_ij,wvl)
  use m_paw_ij
  use defs_wvltypes
  implicit none
  integer,intent(in) :: option
  type(wvl_internal_type), intent(inout) :: wvl
  type(paw_ij_type),intent(in) :: paw_ij(:)
 end subroutine paw2wvl_ij
end interface

interface
 subroutine wvl_cprjreorder(wvl,atm_indx)
  use defs_wvltypes
  implicit none
  type(wvl_internal_type),intent(inout),target :: wvl
  integer,intent(in) :: atm_indx(:)
 end subroutine wvl_cprjreorder
end interface

interface
 subroutine wvl_denspot_free(den)
  use defs_wvltypes
  implicit none
  type(wvl_denspot_type), intent(inout) :: den
 end subroutine wvl_denspot_free
end interface

interface
 subroutine wvl_denspot_set(den,gth_params,ixc,natom,nsppol,rprimd,wvl,&  
  &  wvl_crmult,wvl_frmult,wvl_mpi_comm,xred)
  use defs_basis
  use defs_datatypes
  use defs_wvltypes
  implicit none
  integer,intent(in) :: ixc
  integer,intent(in) :: natom
  integer,intent(in) :: nsppol
  integer,intent(in) :: wvl_mpi_comm
  type(wvl_denspot_type), intent(out) :: den
  type(pseudopotential_gth_type),intent(in) :: gth_params
  type(wvl_internal_type),intent(in) :: wvl
  real(dp), intent(in) :: wvl_crmult
  real(dp), intent(in) :: wvl_frmult
  real(dp), intent(in) :: rprimd(3, 3)
  real(dp), intent(inout) :: xred(3,natom)
 end subroutine wvl_denspot_set
end interface

interface
 subroutine wvl_descr_atoms_set(acell, icoulomb, natom, ntypat, typat, wvl)
  use defs_basis
  use defs_wvltypes
  implicit none
  integer, intent(in) :: icoulomb
  integer, intent(in) :: natom
  integer, intent(in) :: ntypat
  type(wvl_internal_type), intent(inout) :: wvl
  real(dp), intent(in) :: acell(3)
  integer, intent(in) :: typat(natom)
 end subroutine wvl_descr_atoms_set
end interface

interface
 subroutine wvl_descr_atoms_set_sym(wvl, efield, irrzon, nsppol, nsym, phnons,&  
  &  symafm, symrel, tnons, tolsym)
  use defs_basis
  use defs_wvltypes
  implicit none
  integer, intent(in) :: nsppol
  integer, intent(in) :: nsym
  real(dp), intent(in) :: tolsym
  type(wvl_internal_type), intent(inout) :: wvl
  integer, target, intent(in) :: irrzon(:,:,:)
  real(dp), intent(in) :: efield(3)
  real(dp), target, intent(in) :: phnons(:,:,:)
  integer, intent(in) :: symafm(3,3,nsym)
  integer, intent(in) :: symrel(3,3,nsym)
  real(dp), intent(in) :: tnons(3,nsym)
 end subroutine wvl_descr_atoms_set_sym
end interface

interface
 subroutine wvl_descr_free(wvl)
  use defs_wvltypes
  implicit none
  type(wvl_internal_type), intent(inout) :: wvl
 end subroutine wvl_descr_free
end interface

interface
 subroutine wvl_descr_psp_set(filoccup, nsppol, psps, spinat, wvl)
  use defs_basis
  use defs_datatypes
  use defs_wvltypes
  implicit none
  integer, intent(in) :: nsppol
  character(len = *), intent(in) :: filoccup
  type(pseudopotential_type), intent(in) :: psps
  type(wvl_internal_type), intent(inout) :: wvl
  real(dp),intent(in) :: spinat(:,:)
 end subroutine wvl_descr_psp_set
end interface

interface
 subroutine wvl_descr_psp_fill(gth_params, ipsp, ixc, nelpsp, nzatom, pspunit)
  use defs_datatypes
  implicit none
  integer, intent(in) :: ipsp
  integer, intent(in) :: ixc
  integer, intent(in) :: nelpsp
  integer, intent(in) :: nzatom
  integer, intent(in) :: pspunit
  type(pseudopotential_gth_type), intent(inout) :: gth_params
 end subroutine wvl_descr_psp_fill
end interface

interface
 subroutine wvl_paw_free(wvl)
  use defs_wvltypes
  implicit none
  type(wvl_internal_type),intent(inout) :: wvl
 end subroutine wvl_paw_free
end interface

interface
 subroutine wvl_projectors_free(proj)
  use defs_wvltypes
  implicit none
  type(wvl_projectors_type),intent(inout) :: proj
 end subroutine wvl_projectors_free
end interface

interface
 subroutine wvl_projectors_set(me, natom, proj, psps, rprimd, wfs, wvl, wvl_frmult, xred)
  use defs_basis
  use defs_datatypes
  use defs_wvltypes
  implicit none
  integer, intent(in) :: me
  integer, intent(in) :: natom
  type(wvl_projectors_type),intent(inout) :: proj
  type(pseudopotential_type),intent(in) :: psps
  type(wvl_wf_type),intent(in) :: wfs
  type(wvl_internal_type), intent(in) :: wvl
  real(dp), intent(in) :: wvl_frmult
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine wvl_projectors_set
end interface

interface
 subroutine wvl_setBoxGeometry(prtvol, radii, rprimd, xred, wvl, wvl_crmult, wvl_frmult)
  use defs_basis
  use defs_wvltypes
  implicit none
  integer,intent(in) :: prtvol
  type(wvl_internal_type), intent(inout) :: wvl
  real(dp), intent(in) :: wvl_crmult
  real(dp), intent(in) :: wvl_frmult
  real(dp),intent(in) :: radii(:,:)
  real(dp),intent(inout) :: rprimd(3,3)
  real(dp),intent(inout) :: xred(:,:)
 end subroutine wvl_setBoxGeometry
end interface

interface
 subroutine wvl_setngfft(me_wvl, mgfft, nfft, ngfft, nproc_wvl, n1i, n2i, n3i,n3d)
  implicit none
  integer, intent(in) :: me_wvl
  integer, intent(out) :: mgfft
  integer, intent(in) :: n1i
  integer, intent(in) :: n2i
  integer, intent(in) :: n3d
  integer, intent(in) :: n3i
  integer, intent(out) :: nfft
  integer, intent(in) :: nproc_wvl
  integer, intent(out) :: ngfft(18)
 end subroutine wvl_setngfft
end interface

interface
 subroutine wvl_wfs_free(wfs)
  use defs_wvltypes
  implicit none
  type(wvl_wf_type),intent(inout) :: wfs
 end subroutine wvl_wfs_free
end interface

interface
 subroutine wvl_wfs_lr_copy(wfs, wvl)
  use defs_wvltypes
  implicit none
  type(wvl_wf_type), intent(inout) :: wfs
  type(wvl_internal_type), intent(in) :: wvl
 end subroutine wvl_wfs_lr_copy
end interface

interface
 subroutine wvl_wfs_set(alphadiis, spinmagntarget, kpt, me, natom, nband, nkpt, nproc, nspinor,&  
  &  nsppol, nwfshist, occ, psps, rprimd, wfs, wtk, wvl, wvl_crmult, wvl_frmult, xred)
  use defs_basis
  use defs_datatypes
  use defs_wvltypes
  implicit none
  integer, intent(in) :: me
  integer, intent(in) :: natom
  integer, intent(in) :: nband
  integer, intent(in) :: nkpt
  integer, intent(in) :: nproc
  integer, intent(in) :: nspinor
  integer, intent(in) :: nsppol
  integer, intent(in) :: nwfshist
  real(dp), intent(in) :: alphadiis
  type(pseudopotential_type),intent(in) :: psps
  real(dp), intent(in) :: spinmagntarget
  type(wvl_wf_type),intent(out) :: wfs
  type(wvl_internal_type), intent(in) :: wvl
  real(dp), intent(in) :: wvl_crmult
  real(dp), intent(in) :: wvl_frmult
  real(dp), intent(in) :: kpt(3,nkpt)
  real(dp), intent(in) :: occ(:)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp), intent(in) :: wtk(nkpt)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine wvl_wfs_set
end interface

interface
 subroutine derfcf(derfc_yy,yy)
  use defs_basis
  implicit none
  real(dp),intent(out) :: derfc_yy
  real(dp),intent(in) :: yy
 end subroutine derfcf
end interface

interface
 subroutine derf_ab(derf_yy,yy)
  use defs_basis
  implicit none
  real(dp),intent(out) :: derf_yy
  real(dp),intent(in) :: yy
 end subroutine derf_ab
end interface

end module interfaces_43_wvl_wrappers
!!***
