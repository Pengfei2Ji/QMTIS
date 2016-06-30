!!****m* ABINIT/interfaces_56_io_mpi
!! NAME
!! interfaces_56_io_mpi
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/56_io_mpi
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

module interfaces_56_io_mpi

 implicit none

interface
 subroutine hdr_vs_dtset(Hdr,Dtset)
  use defs_abitypes
  implicit none
  type(dataset_type),intent(in) :: Dtset
  type(hdr_type),intent(in) :: Hdr
 end subroutine hdr_vs_dtset
end interface

interface
 subroutine outxfhist(ab_xfh,natom,option,wff2,ios)
  use m_abimover
  use m_wffile
  implicit none
  integer          ,intent(out) :: ios
  integer          ,intent(in) :: natom
  integer          ,intent(in) :: option
  type(ab_xfh_type),intent(inout) :: ab_xfh
  type(wffile_type),intent(inout) :: wff2
 end subroutine outxfhist
end interface

interface
 subroutine rwwf(cg,eigen,formeig,headform,icg,ikpt,isppol,kg_k,mband,mcg,mpi_enreg,&  
  &  nband,nband_disk,npw,nspinor,occ,option,optkg,tim_rwwf,wff)
  use defs_basis
  use defs_abitypes
  use m_wffile
  implicit none
  integer,intent(in) :: formeig
  integer,intent(in) :: headform
  integer,intent(in) :: icg
  integer,intent(in) :: ikpt
  integer,intent(in) :: isppol
  integer,intent(in) :: mband
  integer,intent(in) :: mcg
  integer,intent(in) :: nband
  integer,intent(inout) :: nband_disk
  integer,intent(in) :: npw
  integer,intent(in) :: nspinor
  integer,intent(in) :: option
  integer,intent(in) :: optkg
  integer,intent(in) :: tim_rwwf
  type(mpi_type), intent(inout) :: mpi_enreg
  type(wffile_type),intent(inout) :: wff
  real(dp),intent(inout),target :: cg(2,mcg)
  real(dp),intent(inout),target :: eigen((2*mband)**formeig*mband)
  integer,intent(inout),target :: kg_k(3,optkg*npw)
  real(dp),intent(inout),target :: occ(mband)
 end subroutine rwwf
end interface

interface
 subroutine readwf(cg,eigen,formeig,headform,icg,ikpt,isppol,kg_k,mband,mcg,mpi_enreg,&  
  &  nband,nband_disk,npw,nspinor,occ,option,optkg,wff)
  use defs_basis
  use defs_abitypes
  use m_wffile
  implicit none
  integer,intent(in) :: formeig
  integer,intent(in) :: headform
  integer,intent(in) :: icg
  integer,intent(in) :: ikpt
  integer,intent(in) :: isppol
  integer,intent(in) :: mband
  integer,intent(in) :: mcg
  integer,intent(in) :: nband
  integer,intent(inout) :: nband_disk
  integer,intent(in) :: npw
  integer,intent(in) :: nspinor
  integer,intent(in) :: option
  integer,intent(in) :: optkg
  type(mpi_type),intent(inout) :: mpi_enreg
  type(wffile_type),intent(inout) :: wff
  real(dp),intent(inout),target :: cg(2,mcg)
  real(dp),intent(inout),target :: eigen((2*mband)**formeig*mband)
  integer,intent(inout),target :: kg_k(3,optkg*npw)
  real(dp),intent(inout),target :: occ(mband)
 end subroutine readwf
end interface

interface
 subroutine writewf(cg,eigen,formeig,icg,ikpt,isppol,kg_k,mband,mcg,mpi_enreg,&  
  &  nband,nband_disk,npw,nspinor,occ,option,optkg,wff)
  use defs_basis
  use defs_abitypes
  use m_wffile
  implicit none
  integer,intent(in) :: formeig
  integer,intent(in) :: icg
  integer,intent(in) :: ikpt
  integer,intent(in) :: isppol
  integer,intent(in) :: mband
  integer,intent(in) :: mcg
  integer,intent(in) :: nband
  integer,intent(in) :: nband_disk
  integer,intent(in) :: npw
  integer,intent(in) :: nspinor
  integer,intent(in) :: option
  integer,intent(in) :: optkg
  type(mpi_type),intent(inout) :: mpi_enreg
  type(wffile_type),intent(inout) :: wff
  real(dp),intent(in),target :: cg(2,mcg)
  real(dp),intent(in),target :: eigen((2*mband)**formeig*mband)
  integer,intent(in),target :: kg_k(3,optkg*npw)
  real(dp),intent(in),target :: occ(mband)
 end subroutine writewf
end interface

end module interfaces_56_io_mpi
!!***
