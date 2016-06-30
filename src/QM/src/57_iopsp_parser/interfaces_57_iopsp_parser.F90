!!****m* ABINIT/interfaces_57_iopsp_parser
!! NAME
!! interfaces_57_iopsp_parser
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/57_iopsp_parser
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

module interfaces_57_iopsp_parser

 implicit none

interface
 subroutine inpspheads(filnam,npsp,pspheads,ecut_tmp)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: npsp
  real(dp),intent(inout) :: ecut_tmp(3,2,10)
  character(len=fnlen), intent(in) :: filnam(npsp)
  type(pspheader_type),intent(inout) :: pspheads(npsp)
 end subroutine inpspheads
end interface

interface
 subroutine psml_die(str)
  implicit none
  character(len=*), intent(in) :: str
 end subroutine psml_die
end interface

interface
 subroutine pawpsxml2ab( psxml, pspheads,option )
  use defs_datatypes
  use m_pawxmlps
  implicit none
  integer ,intent(in) :: option
  type(pspheader_type),intent(inout) :: pspheads
  type(paw_setup_t),intent(in) :: psxml
 end subroutine pawpsxml2ab
end interface

interface
 subroutine upfheader2abi (filpsp, znucl, zion, pspxc, lmax_, n1xccc, nproj_l, nprojso_l)
  use defs_basis
  implicit none
  integer, intent(out) :: lmax_
  integer, intent(inout) :: n1xccc
  integer, intent(out) :: pspxc
  character(len=fnlen), intent(in) :: filpsp
  real(dp), intent(out) :: zion
  real(dp), intent(out) :: znucl
  integer, intent(out) :: nproj_l(0:3)
  integer, intent(out) :: nprojso_l(1:3)
 end subroutine upfheader2abi
end interface

interface
 subroutine upfxc2abi(dft, pspxc)
  implicit none
  integer, intent(out) :: pspxc
  character(len=20), intent(in) :: dft
 end subroutine upfxc2abi
end interface

end module interfaces_57_iopsp_parser
!!***
