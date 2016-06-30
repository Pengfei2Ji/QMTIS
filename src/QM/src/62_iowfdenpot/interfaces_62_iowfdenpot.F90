!!****m* ABINIT/interfaces_62_iowfdenpot
!! NAME
!! interfaces_62_iowfdenpot
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/62_iowfdenpot
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

module interfaces_62_iowfdenpot

 implicit none

interface
 subroutine WffReadEigK(eigen,formeig,headform,ikpt,isppol,mband,mpi_enreg,nband,tim_rwwf,wff)
  use defs_basis
  use defs_abitypes
  use m_wffile
  implicit none
  integer,intent(in) :: formeig
  integer,intent(in) :: headform
  integer,intent(in) :: ikpt
  integer,intent(in) :: isppol
  integer,intent(in) :: mband
  integer,intent(in) :: nband
  integer,intent(in) :: tim_rwwf
  type(mpi_type),intent(inout) :: mpi_enreg
  type(wffile_type),intent(inout) :: wff
  real(dp),intent(out) :: eigen((2*mband)**formeig*mband)
 end subroutine WffReadEigK
end interface

interface
 subroutine WffReadSkipK(formeig,headform,ikpt,isppol,mpi_enreg,wff)
  use defs_abitypes
  use m_wffile
  implicit none
  integer,intent(in) :: formeig
  integer,intent(in) :: headform
  integer,intent(in) :: ikpt
  integer,intent(in) :: isppol
  type(mpi_type),intent(inout) :: mpi_enreg
  type(wffile_type),intent(inout) :: wff
 end subroutine WffReadSkipK
end interface

interface
 subroutine initwf(cg,eig_k,formeig,headform,icg,ikpt,ikptsp_old,&  
  &  spin,mcg,mpi_enreg,nband_k,nkpt,npw,nspinor,occ_k,wff1)
  use defs_basis
  use defs_abitypes
  use m_wffile
  implicit none
  integer,intent(in) :: formeig
  integer,intent(in) :: headform
  integer,intent(in) :: icg
  integer,intent(in) :: ikpt
  integer,intent(inout) :: ikptsp_old
  integer,intent(in) :: mcg
  integer,intent(in) :: nband_k
  integer,intent(in) :: nkpt
  integer,intent(in) :: npw
  integer,intent(in) :: nspinor
  integer,intent(in) :: spin
  type(mpi_type),intent(inout) :: mpi_enreg
  type(wffile_type),intent(inout) :: wff1
  real(dp),intent(inout) :: cg(2,mcg)
  real(dp),intent(inout) :: eig_k((2*nband_k)**formeig*nband_k)
  real(dp),intent(inout) :: occ_k(nband_k)
 end subroutine initwf
end interface

interface
 subroutine out1dm(fnameabo_app_1dm,mpi_enreg,natom,nfft,ngfft,nspden,ntypat,&  
  &  rhor,rprimd,typat,ucvol,vtrial,xred,znucl)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  character(len=fnlen),intent(in) :: fnameabo_app_1dm
  type(mpi_type),intent(in) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: rhor(nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: vtrial(nfft,nspden)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: znucl(ntypat)
 end subroutine out1dm
end interface

interface
 subroutine outwant(dtset,eig,cg,kg,npwarr,mband,mcg,nkpt,nsppol,mkmem,mpw,prtwant)
  use defs_basis
  use defs_abitypes
  implicit none
  integer :: mband
  integer :: mcg
  integer :: mkmem
  integer :: mpw
  integer :: nkpt
  integer :: nsppol
  integer :: prtwant
  type(dataset_type),intent(in) :: dtset
  real(dp) :: cg(2,mcg)
  real(dp) :: eig(mband*nkpt*nsppol)
  integer :: kg(3,mpw*mkmem)
  integer :: npwarr(nkpt)
 end subroutine outwant
end interface

interface
 subroutine randac(debug,headform1,ikptsp_prev,ikpt,isppol,nband,nkpt,nsppol,wffinp)
  use m_wffile
  implicit none
  integer,intent(in) :: debug
  integer,intent(in) :: headform1
  integer,intent(in) :: ikpt
  integer,intent(inout) :: ikptsp_prev
  integer,intent(in) :: isppol
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  type(wffile_type),intent(inout) :: wffinp
  integer,intent(in) :: nband(nkpt*nsppol)
 end subroutine randac
end interface

interface
 subroutine rdnpw(ikpt,isppol,nband_k,npw_k,nspinor,option,unitfile)
  implicit none
  integer,intent(in) :: ikpt
  integer,intent(in) :: isppol
  integer,intent(inout) :: nband_k
  integer,intent(inout) :: npw_k
  integer,intent(inout) :: nspinor
  integer,intent(in) :: option
  integer,intent(in) :: unitfile
 end subroutine rdnpw
end interface

end module interfaces_62_iowfdenpot
!!***
