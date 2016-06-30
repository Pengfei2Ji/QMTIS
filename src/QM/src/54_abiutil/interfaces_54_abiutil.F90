!!****m* ABINIT/interfaces_54_abiutil
!! NAME
!! interfaces_54_abiutil
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/54_abiutil
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

module interfaces_54_abiutil

 implicit none

interface
 subroutine calcdensph(gmet,mpi_enreg,natom,nfft,ngfft,nspden,ntypat,nunit,ratsph,rhor,rprimd,typat,ucvol,xred,&  
  &  intgden,atgridpts)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  integer,intent(in) :: nunit
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  integer,intent(inout),optional :: atgridpts(nfft)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(out),optional :: intgden(nspden,natom)
  real(dp),intent(in) :: ratsph(ntypat)
  real(dp),intent(in) :: rhor(nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine calcdensph
end interface

interface
 subroutine find_getdtset(dtsets,getvalue,getname,idtset,iget,miximage,mxnimage,ndtset_alloc)
  use defs_basis
  use defs_abitypes
  implicit none
  integer, intent(in) :: getvalue
  integer, intent(in) :: idtset
  integer, intent(out) :: iget
  integer, intent(in) :: mxnimage
  integer, intent(in) :: ndtset_alloc
  character(len=*),intent(in) :: getname
  type(dataset_type),intent(in) :: dtsets(0:ndtset_alloc)
  real(dp), intent(out) :: miximage(mxnimage,mxnimage)
 end subroutine find_getdtset
end interface

interface
 subroutine mkfilename(filnam,filnam_out,get,idtset,ird,jdtset_,ndtset,stringfil,stringvar,will_read)
  use defs_basis
  implicit none
  integer,intent(in) :: get
  integer,intent(in) :: idtset
  integer,intent(in) :: ird
  integer,intent(in) :: ndtset
  integer,intent(out) :: will_read
  character(len=fnlen),intent(out) :: filnam_out
  character(len=*),intent(in) :: stringfil
  character(len=*),intent(in) :: stringvar
  character(len=fnlen),intent(in) :: filnam(5)
  integer,intent(in) :: jdtset_(0:ndtset)
 end subroutine mkfilename
end interface

end module interfaces_54_abiutil
!!***
