!!****m* ABINIT/interfaces_68_dmft
!! NAME
!! interfaces_68_dmft
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/68_dmft
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

module interfaces_68_dmft

 implicit none

interface
 subroutine compute_levels(cryst_struc,energy_level,hdc,pawang,paw_dmft,nondiag)
  use m_oper
  use m_pawang
  use m_paw_dmft
  use m_crystal
  implicit none
  type(crystal_t),intent(in) :: cryst_struc
  type(oper_type),intent(inout) :: energy_level
  type(oper_type), intent(in) :: hdc
  logical, optional, intent(out) :: nondiag
  type(paw_dmft_type), intent(in) :: paw_dmft
  type(pawang_type), intent(in) :: pawang
 end subroutine compute_levels
end interface

interface
 subroutine datafordmft(cryst_struc,cprj,dimcprj,dtset,eigen,fermie,&  
  &  lda_occup,mband,mband_cprj,mkmem,mpi_enreg,nkpt,my_nspinor,nsppol,occ,&  
  &  paw_dmft,paw_ij,pawang,pawtab,psps,usecprj,unpaw,nbandkss)
  use m_pawtab
  use m_oper
  use m_pawang
  use m_paw_dmft
  use defs_abitypes
  use m_pawcprj
  use m_crystal
  use defs_basis
  use defs_datatypes
  use m_paw_ij
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: mband_cprj
  integer,intent(in) :: mkmem
  integer,intent(in) :: my_nspinor
  integer, optional, intent(in) :: nbandkss
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  integer,intent(in) :: unpaw
  integer,intent(in) :: usecprj
  type(crystal_t),intent(in) :: cryst_struc
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: fermie
  type(oper_type), intent(inout) :: lda_occup
  type(mpi_type),intent(inout) :: mpi_enreg
  type(paw_dmft_type), intent(inout) :: paw_dmft
  type(pawang_type),intent(in) :: pawang
  type(pseudopotential_type),intent(in) :: psps
  type(pawcprj_type), intent(in) :: cprj(cryst_struc%natom,my_nspinor*mband_cprj*mkmem*nsppol*usecprj)
  integer, intent(in) :: dimcprj(cryst_struc%natom)
  real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
  real(dp),intent(in) :: occ(mband*nkpt*nsppol)
  type(paw_ij_type),intent(in) :: paw_ij(:)
  type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)
 end subroutine datafordmft
end interface

interface
 subroutine dmft_solve(cryst_struc,istep,lda_occup,paw_dmft,pawang,pawtab,pawprtvol)
  use m_pawtab
  use m_oper
  use m_pawang
  use m_paw_dmft
  use m_crystal
  implicit none
  integer, intent(in) :: istep
  integer, intent(in) :: pawprtvol
  type(crystal_t),intent(in) :: cryst_struc
  type(oper_type),intent(in) :: lda_occup
  type(paw_dmft_type), intent(inout) :: paw_dmft
  type(pawang_type), intent(in) :: pawang
  type(pawtab_type),intent(in) :: pawtab(cryst_struc%ntypat)
 end subroutine dmft_solve
end interface

interface
 subroutine dyson(green,paw_dmft,self,weiss,opt_weissself)
  use m_self
  use m_paw_dmft
  use m_green
  implicit none
  integer,intent(in) :: opt_weissself
  type(green_type),intent(inout) :: green
  type(paw_dmft_type), intent(inout) :: paw_dmft
  type(self_type),intent(inout) :: self
  type(green_type),intent(inout) :: weiss
 end subroutine dyson
end interface

interface
 subroutine fermi_green(cryst_struc,green,paw_dmft,pawang,&  
  &  self)
  use m_self
  use m_pawang
  use m_paw_dmft
  use m_green
  use m_crystal
  implicit none
  type(crystal_t),intent(in) :: cryst_struc
  type(green_type),intent(inout) :: green
  type(paw_dmft_type), intent(inout) :: paw_dmft
  type(pawang_type),intent(in) :: pawang
  type(self_type), intent(inout) :: self
 end subroutine fermi_green
end interface

interface
 subroutine hubbard_one(cryst_struc,green,hu,paw_dmft,pawang,pawprtvol,hdc,weiss)
  use m_oper
  use m_pawang
  use m_paw_dmft
  use m_hu
  use m_green
  use m_crystal
  implicit none
  integer, intent(in) :: pawprtvol
  type(crystal_t),intent(in) :: cryst_struc
  type(green_type), intent(inout) :: green
  type(oper_type), intent(inout) :: hdc
  type(paw_dmft_type), intent(in) :: paw_dmft
  type(pawang_type), intent(in) :: pawang
  type(green_type), intent(inout) :: weiss
  type(hu_type), intent(inout) :: hu(cryst_struc%ntypat)
 end subroutine hubbard_one
end interface

interface
 subroutine hybridization_asymptotic_coefficient(cryst_struc,paw_dmft,pawang,hybri_coeff)
  use m_matlu
  use m_pawang
  use m_paw_dmft
  use m_crystal
  implicit none
  type(crystal_t),intent(in) :: cryst_struc
  type(paw_dmft_type), intent(in) :: paw_dmft
  type(pawang_type), intent(in) :: pawang
  type(matlu_type), intent(inout) :: hybri_coeff(paw_dmft%natom)
 end subroutine hybridization_asymptotic_coefficient
end interface

interface
 subroutine impurity_solve(cryst_struc,green,hu,paw_dmft,&  
  &  pawang,pawtab,self_old,self_new,weiss,pawprtvol)
  use m_pawang
  use m_paw_dmft
  use m_self
  use m_hu
  use m_green
  use m_crystal
  use m_pawtab
  implicit none
  integer, intent(in) :: pawprtvol
  type(crystal_t),intent(in) :: cryst_struc
  type(green_type), intent(inout) :: green
  type(paw_dmft_type), intent(inout) :: paw_dmft
  type(pawang_type), intent(in) :: pawang
  type(self_type), intent(inout) :: self_new
  type(self_type), intent(inout) :: self_old
  type(green_type), intent(inout) :: weiss
  type(hu_type),intent(inout) :: hu(cryst_struc%ntypat)
  type(pawtab_type),intent(in) :: pawtab(cryst_struc%ntypat)
 end subroutine impurity_solve
end interface

interface
 subroutine ldau_self(cryst_struc,green,paw_dmft,pawtab,self,opt_ldau,prtopt)
  use m_self
  use m_pawtab
  use m_paw_dmft
  use m_green
  use m_crystal
  implicit none
  integer, intent(in) :: opt_ldau
  integer, intent(in) :: prtopt
  type(crystal_t),intent(in) :: cryst_struc
  type(green_type), intent(in) :: green
  type(paw_dmft_type), intent(in) :: paw_dmft
  type(self_type), intent(inout) :: self
  type(pawtab_type),intent(in) :: pawtab(cryst_struc%ntypat)
 end subroutine ldau_self
end interface

interface
 subroutine local_ks_green(green,paw_dmft,prtopt)
  use m_paw_dmft
  use m_green
  implicit none
  integer, intent(in) :: prtopt
  type(green_type), intent(in) :: green
  type(paw_dmft_type), intent(in) :: paw_dmft
 end subroutine local_ks_green
end interface

interface
 subroutine newton(cryst_struc,green,paw_dmft,pawang,self,&  
  &  x_input,x_precision,max_iter,f_precision,ierr_hh,opt_noninter,opt_algo)
  use m_pawang
  use m_paw_dmft
  use m_self
  use m_green
  use m_crystal
  use defs_basis
  implicit none
  integer,intent(out) :: ierr_hh
  integer,intent(in) :: max_iter
  integer,intent(in) :: opt_noninter
  type(crystal_t),intent(in) :: cryst_struc
  real(dp),intent(inout) :: f_precision
  type(green_type),intent(inout) :: green
  real(dp),intent(in), optional :: opt_algo
  type(paw_dmft_type), intent(inout) :: paw_dmft
  type(pawang_type),intent(in) :: pawang
  type(self_type), intent(inout) :: self
  real(dp),intent(inout) :: x_input
  real(dp),intent(inout) :: x_precision
 end subroutine newton
end interface

interface
 subroutine psichi_renormalization(cryst_struc,paw_dmft,pawang,opt)
  use m_pawang
  use m_paw_dmft
  use m_crystal
  implicit none
  integer, optional, intent(in) :: opt
  type(crystal_t),intent(in) :: cryst_struc
  type(paw_dmft_type), intent(inout) :: paw_dmft
  type(pawang_type), intent(in) :: pawang
 end subroutine psichi_renormalization
end interface

interface
 subroutine qmc_prep_ctqmc(cryst_struc,green,self,hu,paw_dmft,pawang,pawprtvol,weiss)
  use m_pawang
  use m_paw_dmft
  use m_self
  use m_hu
  use m_green
  use m_crystal
  implicit none
  integer, intent(in) :: pawprtvol
  type(crystal_t),intent(in) :: cryst_struc
  type(green_type), intent(inout) :: green
  type(paw_dmft_type), intent(inout) :: paw_dmft
  type(pawang_type), intent(in) :: pawang
  type(self_type), intent(in) :: self
  type(green_type), intent(inout) :: weiss
  type(hu_type), intent(in) :: hu(cryst_struc%ntypat)
 end subroutine qmc_prep_ctqmc
end interface

interface
 subroutine spectral_function(cryst_struc,green,hu,paw_dmft,&  
  &  pawang,pawtab,self_old,prtopt)
  use m_pawang
  use m_paw_dmft
  use m_self
  use m_hu
  use m_green
  use m_crystal
  use m_pawtab
  implicit none
  integer, intent(in) :: prtopt
  type(crystal_t),intent(in) :: cryst_struc
  type(green_type), intent(in) :: green
  type(paw_dmft_type), intent(inout) :: paw_dmft
  type(pawang_type), intent(in) :: pawang
  type(self_type), intent(inout) :: self_old
  type(hu_type),intent(inout) :: hu(cryst_struc%ntypat)
  type(pawtab_type),intent(in) :: pawtab(cryst_struc%ntypat)
 end subroutine spectral_function
end interface

interface
 subroutine testcode_ctqmc(dmftqmc_l,fw1_nd,fw1,gtmp_nd,gw_tmp_nd,levels_ctqmc,hybri_limit,&  
  &  nflavor,opt,temp,testrot,testcode,umod)
  use defs_basis
  implicit none
  integer, intent(in) :: dmftqmc_l
  integer, intent(in) :: nflavor
  integer, intent(in) :: opt
  integer, intent(in) :: testcode
  integer, intent(in) :: testrot
  real(dp), intent(in) :: temp
  complex(dpc), intent(out) :: fw1(:,:)
  complex(dpc), intent(out) :: fw1_nd(:,:,:)
  real(dp),  intent(inout) :: gtmp_nd(:,:,:)
  complex(dpc), intent(inout) :: gw_tmp_nd(:,:,:)
  complex(dpc),  intent(inout) :: hybri_limit(:,:)
  real(dp),  intent(inout) :: levels_ctqmc(:)
  real(dp), intent(out) :: umod(2,2)
 end subroutine testcode_ctqmc
end interface

end module interfaces_68_dmft
!!***
