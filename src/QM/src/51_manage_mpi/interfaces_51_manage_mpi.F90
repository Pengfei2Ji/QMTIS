!!****m* ABINIT/interfaces_51_manage_mpi
!! NAME
!! interfaces_51_manage_mpi
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/51_manage_mpi
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

module interfaces_51_manage_mpi

 implicit none

interface
 subroutine distrb2(mband,nband,nkpt,nproc,nsppol,mpi_enreg)
  use defs_abitypes
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: nkpt
  integer,intent(in) :: nproc
  integer,intent(in) :: nsppol
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: nband(nkpt*nsppol)
 end subroutine distrb2
end interface

interface
 subroutine distrb2_hf(nbandhf,nkpthf, nproc, nsppol, mpi_enreg)
  use defs_abitypes
  implicit none
  integer,intent(in) :: nbandhf
  integer,intent(in) :: nkpthf
  integer,intent(in) :: nproc
  integer,intent(in) :: nsppol
  type(mpi_type),intent(inout) :: mpi_enreg
 end subroutine distrb2_hf
end interface

interface
 subroutine get_npert_rbz(dtset,nband_rbz,nkpt_rbz,npert)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(out) :: npert
  type(dataset_type),intent(in) :: dtset
  integer,pointer :: nkpt_rbz(:)
  real(dp),pointer :: nband_rbz(:,:)
 end subroutine get_npert_rbz
end interface

interface
 subroutine herald(code_name,code_version,iout)
  implicit none
  integer,intent(in) :: iout
  character(len=*),intent(in) :: code_name
  character(len=*),intent(in) :: code_version
 end subroutine herald
end interface

interface
 subroutine initmpi_atom(dtset,mpi_enreg)
  use defs_abitypes
  implicit none
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
 end subroutine initmpi_atom
end interface

interface
 subroutine initmpi_band(mpi_enreg,nband,nkpt,nsppol)
  use defs_abitypes
  implicit none
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: nband(nkpt*nsppol)
 end subroutine initmpi_band
end interface

interface
 subroutine initmpi_grid(mpi_enreg)
  use defs_abitypes
  implicit none
  type(mpi_type),intent(inout) :: mpi_enreg
 end subroutine initmpi_grid
end interface

interface
 subroutine initmpi_img(dtset,mpi_enreg,option)
  use defs_abitypes
  implicit none
  integer,intent(in) :: option
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
 end subroutine initmpi_img
end interface

interface
 subroutine initmpi_pert(dtset,mpi_enreg)
  use defs_abitypes
  implicit none
  type(dataset_type),intent(inout) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
 end subroutine initmpi_pert
end interface

interface
 subroutine initmpi_seq(mpi_enreg)
  use defs_abitypes
  implicit none
  type(mpi_type),intent(out) :: mpi_enreg
 end subroutine initmpi_seq
end interface

interface
 subroutine initmpi_world(mpi_enreg,nproc)
  use defs_abitypes
  implicit none
  integer, intent(in) :: nproc
  type(mpi_type),intent(inout) :: mpi_enreg
 end subroutine initmpi_world
end interface

interface
 subroutine pspheads_comm(npsp,pspheads,test_paw)
  use defs_datatypes
  implicit none
  integer,intent(in) :: npsp
  integer,intent(inout) :: test_paw
  type(pspheader_type),intent(inout) :: pspheads(npsp)
 end subroutine pspheads_comm
end interface

end module interfaces_51_manage_mpi
!!***
