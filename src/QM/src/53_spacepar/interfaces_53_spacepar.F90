!!****m* ABINIT/interfaces_53_spacepar
!! NAME
!! interfaces_53_spacepar
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/53_spacepar
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

module interfaces_53_spacepar

 implicit none

interface
 subroutine dotprod_vn(cplex,dens,dotr,doti,nfft,nfftot,nspden,option,pot,&  
  &  ucvol,mpi_comm_sphgrid)
  use defs_basis
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in),optional :: mpi_comm_sphgrid
  integer,intent(in) :: nfft
  integer,intent(in) :: nfftot
  integer,intent(in) :: nspden
  integer,intent(in) :: option
  real(dp),intent(out) :: doti
  real(dp),intent(out) :: dotr
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: dens(cplex*nfft,nspden)
  real(dp),intent(in) :: pot(cplex*nfft,nspden)
 end subroutine dotprod_vn
end interface

interface
 subroutine meanvalue_g(ar,diag,filter,istwf_k,mpi_enreg,npw,nspinor,vect,vect1,use_ndo,ar_im)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: filter
  integer,intent(in) :: istwf_k
  integer,intent(in) :: npw
  integer,intent(in) :: nspinor
  integer,intent(in) :: use_ndo
  real(dp),intent(out) :: ar
  real(dp),intent(out),optional :: ar_im
  type(mpi_type),intent(in) :: mpi_enreg
  real(dp),intent(in) :: diag(npw)
  real(dp),intent(in) :: vect(2,npw*nspinor)
  real(dp),intent(in) :: vect1(2,npw*nspinor)
 end subroutine meanvalue_g
end interface

interface
 subroutine multipoles_fftr(arraysp,dipole,mpi_enreg,nfft,ngfft,nspden,rprimd,neworigin)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  type(mpi_type),intent(in) :: mpi_enreg
  integer,intent(in) :: ngfft(3)
  real(dp),intent(in) :: arraysp(nfft,nspden)
  real(dp),intent(out) :: dipole(3,nspden)
  real(dp),intent(in) :: neworigin(3)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine multipoles_fftr
end interface

interface
 subroutine multipoles_out(arraysp,mpi_enreg,natom,nfft,ngfft,nspden,&  
  &  ntypat,rprimd,typat,ucvol,xred,ziontypat)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  type(mpi_type),intent(in) :: mpi_enreg
  real(dp), intent(in) :: ucvol
  integer,intent(in) :: ngfft(3)
  real(dp),intent(in) :: arraysp(nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  integer, intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: ziontypat(ntypat)
 end subroutine multipoles_out
end interface

end module interfaces_53_spacepar
!!***
