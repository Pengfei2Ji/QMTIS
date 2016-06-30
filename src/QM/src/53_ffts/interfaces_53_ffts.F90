!!****m* ABINIT/interfaces_53_ffts
!! NAME
!! interfaces_53_ffts
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/53_ffts
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

module interfaces_53_ffts

 implicit none

interface
 subroutine ccfft(ngfft,isign,n1,n2,n3,n4,n5,n6,ndat,option,work1,work2,comm_fft)
  use defs_basis
  implicit none
  integer,intent(in) :: comm_fft
  integer,intent(in) :: isign
  integer,intent(in) :: n1
  integer,intent(in) :: n2
  integer,intent(in) :: n3
  integer,intent(in) :: n4
  integer,intent(in) :: n5
  integer,intent(in) :: n6
  integer,intent(in) :: ndat
  integer,intent(in) :: option
  integer,intent(in) :: ngfft(18)
  real(dp),intent(inout) :: work1(2,n4*n5*n6*ndat)
  real(dp),intent(inout) :: work2(2,n4*n5*n6*ndat)
 end subroutine ccfft
end interface

interface
 subroutine fftpac(ispden,mpi_enreg,nspden,n1,n2,n3,nd1,nd2,nd3,ngfft,aa,bb,option)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: ispden
  integer,intent(in) :: n1
  integer,intent(in) :: n2
  integer,intent(in) :: n3
  integer,intent(in) :: nd1
  integer,intent(in) :: nd2
  integer,intent(in) :: nd3
  integer,intent(in) :: nspden
  integer,intent(in) :: option
  type(mpi_type),intent(in) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(inout) :: aa(n1*n2*n3/ngfft(10),nspden)
  real(dp),intent(inout) :: bb(nd1,nd2,nd3)
 end subroutine fftpac
end interface

interface
 subroutine fourdp(cplex,fofg,fofr,isign,mpi_enreg,nfft,ngfft,paral_kgb,tim_fourdp)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: isign
  integer,intent(in) :: nfft
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: tim_fourdp
  type(mpi_type),intent(in) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(inout) :: fofg(2,nfft)
  real(dp),intent(inout) :: fofr(cplex*nfft)
 end subroutine fourdp
end interface

interface
 subroutine fourdp_6d(cplex,matrix,isign,MPI_enreg,nfft,ngfft,paral_kgb,tim_fourdp)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: isign
  integer,intent(in) :: nfft
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: tim_fourdp
  type(mpi_type),intent(inout) :: MPI_enreg
  integer,intent(in) :: ngfft(18)
  complex(gwpc),intent(inout) :: matrix(nfft,nfft)
 end subroutine fourdp_6d
end interface

interface
 subroutine fourwf(cplex,denpot,fofgin,fofgout,fofr,gboundin,gboundout,istwf_k,&  
  &  kg_kin,kg_kout,mgfft,mpi_enreg,ndat,ngfft,npwin,npwout,n4,n5,n6,option,&  
  &  paral_kgb,tim_fourwf,weight_r,weight_i,&  
  &  use_gpu_cuda,use_ndo,fofginb) ! Optional arguments
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: istwf_k
  integer,intent(in) :: mgfft
  integer,intent(in) :: n4
  integer,intent(in) :: n5
  integer,intent(in) :: n6
  integer,intent(in) :: ndat
  integer,intent(in) :: npwin
  integer,intent(in) :: npwout
  integer,intent(in) :: option
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: tim_fourwf
  integer,intent(in),optional :: use_gpu_cuda
  integer,intent(in),optional :: use_ndo
  type(mpi_type),intent(in) :: mpi_enreg
  real(dp),intent(in) :: weight_i
  real(dp),intent(in) :: weight_r
  integer,intent(in) :: ngfft(18)
  real(dp),intent(inout) :: denpot(cplex*n4,n5,n6)
  real(dp),intent(inout) :: fofgin(2,npwin*ndat)
  real(dp),intent(inout),optional :: fofginb(:,:)
  real(dp),intent(out) :: fofgout(2,npwout*ndat)
  real(dp),intent(inout) :: fofr(2,n4,n5,n6*ndat)
  integer,intent(in) :: gboundin(2*mgfft+8,2)
  integer,intent(in) :: gboundout(2*mgfft+8,2)
  integer,intent(in) :: kg_kin(3,npwin)
  integer,intent(in) :: kg_kout(3,npwout)
 end subroutine fourwf
end interface

interface
 subroutine kgindex(indpw_k,kg_k,mask,mpi_enreg,ngfft,npw_k)
  use defs_abitypes
  implicit none
  integer,intent(in) :: npw_k
  type(mpi_type),intent(in) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  integer,intent(out) :: indpw_k(npw_k)
  integer,intent(in) :: kg_k(3,npw_k)
  logical,intent(out) :: mask(npw_k)
 end subroutine kgindex
end interface

interface
 subroutine zerosym(array,cplex,n1,n2,n3,&  
  &  ig1,ig2,ig3,comm_fft,distribfft) ! Optional arguments
  use defs_basis
  use m_distribfft
  implicit none
  integer,optional,intent(in) :: comm_fft
  integer,intent(in) :: cplex
  integer,optional,intent(in) :: ig1
  integer,optional,intent(in) :: ig2
  integer,optional,intent(in) :: ig3
  integer,intent(in) :: n1
  integer,intent(in) :: n2
  integer,intent(in) :: n3
  type(distribfft_type),intent(in),target,optional :: distribfft
  real(dp),intent(inout) :: array(cplex,n1*n2*n3)
 end subroutine zerosym
end interface

end module interfaces_53_ffts
!!***
