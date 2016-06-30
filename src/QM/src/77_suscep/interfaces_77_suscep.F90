!!****m* ABINIT/interfaces_77_suscep
!! NAME
!! interfaces_77_suscep
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/77_suscep
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

module interfaces_77_suscep

 implicit none

interface
 subroutine suscep_stat(atindx,atindx1,cg,cprj,dielar,dimcprj,doccde,&  
  &  eigen,gbound_diel,gprimd,irrzondiel,istwfk,kg,&  
  &  kg_diel,lmax_diel,&  
  &  mband,mcg,mcprj,mgfftdiel,mkmem,mpi_enreg,mpw,natom,nband,&  
  &  neglect_pawhat,nfftdiel,ngfftdiel,nkpt,npwarr,&  
  &  npwdiel,nspden,nspinor,nsppol,nsym,ntypat,occ,occopt,&  
  &  pawang,pawtab,phnonsdiel,ph1ddiel,rprimd,&  
  &  susmat,symafm,symrel,tnons,typat,ucvol,unpaw,usecprj,usepaw,usetimerev,&  
  &  wtk,ylmdiel)
  use defs_basis
  use defs_abitypes
  use m_pawcprj
  use m_pawang
  use m_pawtab
  implicit none
  integer,intent(in) :: lmax_diel
  integer,intent(in) :: mband
  integer,intent(in) :: mcg
  integer,intent(in) :: mcprj
  integer,intent(in) :: mgfftdiel
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: natom
  integer,intent(in) :: neglect_pawhat
  integer,intent(in) :: nfftdiel
  integer,intent(in) :: nkpt
  integer,intent(in) :: npwdiel
  integer,intent(in) :: nspden
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(in) :: occopt
  integer,intent(in) :: unpaw
  integer,intent(in) :: usecprj
  integer,intent(in) :: usepaw
  integer,intent(in) :: usetimerev
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfftdiel(18)
  integer,intent(in) :: atindx(natom)
  integer,intent(in) :: atindx1(natom)
  real(dp),intent(in) :: cg(2,mcg)
  type(pawcprj_type) :: cprj(natom,mcprj*usecprj)
  real(dp),intent(in) :: dielar(7)
  integer,intent(in) :: dimcprj(natom*usepaw)
  real(dp),intent(in) :: doccde(mband*nkpt*nsppol)
  real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
  integer,intent(in) :: gbound_diel(2*mgfftdiel+8,2)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: irrzondiel(nfftdiel**(1-1/nsym),2,(nspden/nsppol)-3*(nspden/4))
  integer,intent(in) :: istwfk(nkpt)
  integer,intent(in) :: kg(3,mpw*mkmem)
  integer,intent(in) :: kg_diel(3,npwdiel)
  integer,intent(in) :: nband(nkpt*nsppol)
  integer,intent(in) :: npwarr(nkpt)
  real(dp),intent(in) :: occ(mband*nkpt*nsppol)
  type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)
  real(dp),intent(in) :: ph1ddiel(2,(3*(2*mgfftdiel+1)*natom)*usepaw)
  real(dp),intent(in) :: phnonsdiel(2,nfftdiel**(1-1/nsym),(nspden/nsppol)-3*(nspden/4))
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: susmat(2,npwdiel,nspden,npwdiel,nspden)
  integer,intent(in) :: symafm(nsym)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
  integer,intent(in) :: typat(ntypat)
  real(dp),intent(in) :: wtk(nkpt)
  real(dp),intent(in) :: ylmdiel(npwdiel,lmax_diel**2)
 end subroutine suscep_stat
end interface

interface
 subroutine susk(atindx,bdtot_index,cg_mpi,cprj_k,doccde,drhode,eigen,extrap,gbound,&  
  &  gbound_diel,gylmg_diel,icg_mpi,ikpt,isp,istwf_k,kg_diel,kg_k_mpi,&  
  &  lmax_diel,mband,mcg,mgfftdiel,mpi_enreg,&  
  &  natom,nband_k,ndiel4,ndiel5,ndiel6,neglect_pawhat,ngfftdiel,nkpt,&  
  &  npwdiel,npw_k_mpi,nspden,nspden_eff,nspinor,nsppol,ntypat,occ,occopt,occ_deavg,&  
  &  pawang,pawtab,ph3d_diel,rhoextrap,sumdocc,&  
  &  susmat,typat,ucvol,usepaw,wtk)
  use defs_basis
  use defs_abitypes
  use m_pawcprj
  use m_pawang
  use m_pawtab
  implicit none
  integer,intent(in) :: bdtot_index
  integer,intent(in) :: extrap
  integer,intent(in),target :: icg_mpi
  integer,intent(in) :: ikpt
  integer,intent(in) :: isp
  integer,intent(in) :: istwf_k
  integer,intent(in) :: lmax_diel
  integer,intent(in) :: mband
  integer,intent(in) :: mcg
  integer,intent(in) :: mgfftdiel
  integer,intent(in) :: natom
  integer,intent(in) :: nband_k
  integer,intent(in) :: ndiel4
  integer,intent(in) :: ndiel5
  integer,intent(in) :: ndiel6
  integer,intent(in) :: neglect_pawhat
  integer,intent(in) :: nkpt
  integer,intent(in),target :: npw_k_mpi
  integer,intent(in) :: npwdiel
  integer,intent(in) :: nspden
  integer,intent(in) :: nspden_eff
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: ntypat
  integer,intent(in) :: occopt
  integer,intent(in) :: usepaw
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  real(dp),intent(inout) :: sumdocc
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfftdiel(18)
  integer,intent(in) :: atindx(natom)
  real(dp),intent(in),target :: cg_mpi(2,mcg)
  type(pawcprj_type) :: cprj_k(natom,nspinor*nband_k*usepaw)
  real(dp),intent(in) :: doccde(mband*nkpt*nsppol)
  real(dp),intent(inout) :: drhode(2,npwdiel,nspden_eff)
  real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
  integer,intent(inout) :: gbound(2*mgfftdiel+8,2)
  integer,intent(in) :: gbound_diel(2*mgfftdiel+8,2)
  real(dp),intent(in) :: gylmg_diel(npwdiel,lmax_diel**2,ntypat*usepaw)
  integer,intent(in) :: kg_diel(3,npwdiel)
  integer,intent(in),target :: kg_k_mpi(3,npw_k_mpi)
  real(dp),intent(in) :: occ(mband*nkpt*nsppol)
  real(dp),intent(in) :: occ_deavg(mband)
  type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)
  real(dp),intent(in) :: ph3d_diel(2,npwdiel,natom*usepaw)
  real(dp),intent(inout) :: rhoextrap(ndiel4,ndiel5,ndiel6,nspinor)
  real(dp),intent(inout) :: susmat(2,npwdiel,nspden,npwdiel,nspden)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: wtk(nkpt)
 end subroutine susk
end interface

interface
 subroutine suskmm(atindx,bdtot_index,cg,cprj_k,doccde,drhode,eigen,extrap,gbound,&  
  &  gbound_diel,gylmg_diel,icg,ikpt,isp,istwf_k,kg_diel,kg_k,&  
  &  lmax_diel,mband,mcg,mgfftdiel,mpi_enreg,&  
  &  natom,nband_k,ndiel4,ndiel5,ndiel6,neglect_pawhat,ngfftdiel,nkpt,&  
  &  npwdiel,npw_k,nspden,nspden_eff,nspinor,nsppol,ntypat,occ,occopt,occ_deavg,paral_kgb,&  
  &  pawang,pawtab,ph3d_diel,rhoextrap,sumdocc,&  
  &  susmat,typat,ucvol,usepaw,wtk)
  use defs_basis
  use defs_abitypes
  use m_pawcprj
  use m_pawang
  use m_pawtab
  implicit none
  integer,intent(in) :: bdtot_index
  integer,intent(in) :: extrap
  integer,intent(in) :: icg
  integer,intent(in) :: ikpt
  integer,intent(in) :: isp
  integer,intent(in) :: istwf_k
  integer,intent(in) :: lmax_diel
  integer,intent(in) :: mband
  integer,intent(in) :: mcg
  integer,intent(in) :: mgfftdiel
  integer,intent(in) :: natom
  integer,intent(in) :: nband_k
  integer,intent(in) :: ndiel4
  integer,intent(in) :: ndiel5
  integer,intent(in) :: ndiel6
  integer,intent(in) :: neglect_pawhat
  integer,intent(in) :: nkpt
  integer,intent(in) :: npw_k
  integer,intent(in) :: npwdiel
  integer,intent(in) :: nspden
  integer,intent(in) :: nspden_eff
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: ntypat
  integer,intent(in) :: occopt
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: usepaw
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  real(dp),intent(inout) :: sumdocc
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfftdiel(18)
  integer,intent(in) :: atindx(natom)
  real(dp),intent(in) :: cg(2,mcg)
  type(pawcprj_type) :: cprj_k(natom,nspinor*nband_k*usepaw)
  real(dp),intent(in) :: doccde(mband*nkpt*nsppol)
  real(dp),intent(inout) :: drhode(2,npwdiel,nspden_eff)
  real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
  integer,intent(in) :: gbound(2*mgfftdiel+8,2)
  integer,intent(in) :: gbound_diel(2*mgfftdiel+8,2)
  real(dp),intent(in) :: gylmg_diel(npwdiel,lmax_diel**2,ntypat*usepaw)
  integer,intent(in) :: kg_diel(3,npwdiel)
  integer,intent(in) :: kg_k(3,npw_k)
  real(dp),intent(in) :: occ(mband*nkpt*nsppol)
  real(dp),intent(in) :: occ_deavg(mband)
  type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)
  real(dp),intent(in) :: ph3d_diel(2,npwdiel,natom*usepaw)
  real(dp),intent(inout) :: rhoextrap(ndiel4,ndiel5,ndiel6,nspinor)
  real(dp),intent(inout) :: susmat(2,npwdiel,nspden,npwdiel,nspden)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: wtk(nkpt)
 end subroutine suskmm
end interface

end module interfaces_77_suscep
!!***
