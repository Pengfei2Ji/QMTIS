!!****m* ABINIT/interfaces_20_datashare
!! NAME
!! interfaces_20_datashare
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/20_datashare
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

module interfaces_20_datashare

 implicit none

interface
 subroutine vdw_dftd3_data(vdw_dftd3_r0,vdw_dftd3_c6,index_c6,vdw_dftd3_cni,index_cni,&  
  &  vdw_dftd3_cnj,index_cnj)
  use defs_basis
  implicit none
  integer,intent(out) :: index_c6(254)
  integer,intent(out) :: index_cni(27884)
  integer,intent(out) :: index_cnj(13171)
  real(dp),intent(out) :: vdw_dftd3_c6(32385)
  real(dp),intent(out) :: vdw_dftd3_cni(27884)
  real(dp),intent(out) :: vdw_dftd3_cnj(13171)
  real(dp),intent(out) :: vdw_dftd3_r0(4465)
 end subroutine vdw_dftd3_data
end interface

end module interfaces_20_datashare
!!***
