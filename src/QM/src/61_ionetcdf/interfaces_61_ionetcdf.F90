!!****m* ABINIT/interfaces_61_ionetcdf
!! NAME
!! interfaces_61_ionetcdf
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/61_ionetcdf
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

module interfaces_61_ionetcdf

 implicit none

interface
 subroutine write_eig(eigen,filename,kptns,mband,nband,nkpt,nsppol)
  use defs_basis
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  character(len=fnlen),intent(in) :: filename
  real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
  real(dp),intent(in) :: kptns(3,nkpt)
  integer,intent(in) :: nband(nkpt*nsppol)
 end subroutine write_eig
end interface

interface
 subroutine wrt_moldyn_netcdf(amass,dtset,itime,option,moldyn_file,mpi_enreg,&  
  &  results_gs,rprimd,unpos,vel,xcart,xred)
  use defs_basis
  use defs_abitypes
  use m_results_gs
  implicit none
  integer,intent(in) :: itime
  integer,intent(in) :: option
  integer,intent(in) :: unpos
  type(dataset_type),intent(in) :: dtset
  character(fnlen),intent(in) :: moldyn_file
  type(mpi_type),intent(in) :: mpi_enreg
  type(results_gs_type),intent(in) :: results_gs
  real(dp),intent(in) :: amass(dtset%natom)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in),target :: vel(3,dtset%natom)
  real(dp),intent(in) :: xcart(3,dtset%natom)
  real(dp),intent(in) :: xred(3,dtset%natom)
 end subroutine wrt_moldyn_netcdf
end interface

end module interfaces_61_ionetcdf
!!***
