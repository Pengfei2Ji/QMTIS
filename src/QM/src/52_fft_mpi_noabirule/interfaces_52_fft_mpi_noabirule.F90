!!****m* ABINIT/interfaces_52_fft_mpi_noabirule
!! NAME
!! interfaces_52_fft_mpi_noabirule
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/52_fft_mpi_noabirule
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

module interfaces_52_fft_mpi_noabirule

 implicit none

interface
 subroutine bound(dsqmax,dsqmin,gbound,gmet,kpt,ngfft,plane)
  use defs_basis
  implicit none
  integer,intent(out) :: plane
  real(dp),intent(out) :: dsqmax
  real(dp),intent(out) :: dsqmin
  integer,intent(out) :: gbound(3)
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: kpt(3)
 end subroutine bound
end interface

interface
 subroutine getng(boxcutmin,ecut,gmet,kpt,me_fft,mgfft,nfft,ngfft,nproc_fft,nsym,paral_fft,symrel,&  
  &  ngfftc,use_gpu_cuda,unit) ! optional
  use defs_basis
  implicit none
  integer,intent(in) :: me_fft
  integer,intent(out) :: mgfft
  integer,intent(out) :: nfft
  integer,intent(in) :: nproc_fft
  integer,intent(in) :: nsym
  integer,intent(in) :: paral_fft
  integer,optional,intent(in) :: unit
  integer,optional,intent(in) :: use_gpu_cuda
  real(dp),intent(in) :: boxcutmin
  real(dp),intent(in) :: ecut
  integer,intent(inout) :: ngfft(18)
  integer,intent(in),optional :: ngfftc(3)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: kpt(3)
  integer,intent(in) :: symrel(3,3,nsym)
 end subroutine getng
end interface

interface
 subroutine indirect_parallel_Fourier(index,left,mpi_enreg,ngleft,ngright,nleft,nright,paral_kgb,right,sizeindex)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: nleft
  integer,intent(in) :: nright
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: sizeindex
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: ngleft(18)
  integer,intent(in) :: ngright(18)
  integer,intent(in) :: index(sizeindex)
  real(dp),intent(inout) :: left(2,nleft)
  real(dp),intent(in) :: right(2,nright)
 end subroutine indirect_parallel_Fourier
end interface

interface
 subroutine sphereboundary(gbound,istwf_k,kg_k,mgfft,npw)
  implicit none
  integer,intent(in) :: istwf_k
  integer,intent(in) :: mgfft
  integer,intent(in) :: npw
  integer,intent(out) :: gbound(2*mgfft+8,2)
  integer,intent(in) :: kg_k(3,npw)
 end subroutine sphereboundary
end interface

end module interfaces_52_fft_mpi_noabirule
!!***
