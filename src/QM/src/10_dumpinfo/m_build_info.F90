!{src2tex{textfont=tt}}
!!****m* ABINIT/m_build_info
!! NAME
!!  m_build_info
!!
!! FUNCTION
!!  This module contains information about this particular version of ABINIT
!!  and its build parameters (useful for debugging).
!!
!! COPYRIGHT
!!  Copyright (C) 2005-2016 ABINIT group (Yann Pouillon, Matteo Giantomassi)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module m_build_info

implicit none

! Try to prevent problems with the length of the argument 
! (should not exceed the line length of 132 character).
! So, always start the string to be replaced at the beginning of a line.

! Parameters set-up by Autoconf
character(len=8),parameter :: abinit_version = &
&   "8.0.6"
character(len=*),parameter :: build_target   = &
&   "x86_64_linux_gnu4.8"

! More info on current version
character(len=*),parameter :: version_major = &
&   "8"
character(len=*),parameter :: version_minor = &
&   "0"
character(len=*),parameter :: version_micro = &
&   "6"
character(len=*),parameter :: version_build = &
&   "20160629"

! Info on compilers. Try to prevent problems with the length of the argument 
! (should not exceed the line length of 132 character).
character(len=*),parameter :: cc_info    = &
&   "gnu4.8"
character(len=*),parameter :: cxx_info   = &
&   "gnu4.8"
character(len=*),parameter :: fc_info    = &
&   "gnu4.8"
character(len=*),parameter :: cc_flags   = &
&   " -g -O2 -mtune=native -march=native  "
character(len=*),parameter :: cxx_flags  = &
&   " -g -O2 -mtune=native -march=native  "
character(len=*),parameter :: fc_flags   = &
&   " -g -ffree-line-length-none  "
character(len=*),parameter :: fc_ldflags = &
&   "    "

! Info on optimizations
character(len=*),parameter :: enable_debug    = &
&   "basic"
character(len=*),parameter :: enable_optim = &
&   "standard"
character(len=*),parameter :: cpu_info = &
&   "unknown_unknown"

! Info on MPI
character(len=*),parameter :: enable_mpi       = &
&   "no"
character(len=*),parameter :: enable_mpi_io    = &
&   "no"

! Info on openMP
character(len=*),parameter :: enable_openmp = &
&   "no"

! Info on GPU
character(len=*),parameter :: enable_gpu    = &
&   "no"

! Info on connectors / fallbacks
character(len=*),parameter :: enable_connectors = &
&   "yes"
character(len=*),parameter :: enable_fallbacks  = &
&   "yes"
character(len=*),parameter :: dft_flavor        = &
&   "none"
character(len=*),parameter :: fft_flavor        = &
&   "none"
character(len=*),parameter :: linalg_flavor     = &
&   "netlib-fallback"
character(len=*),parameter :: math_flavor       = &
&   "none"
character(len=*),parameter :: timer_flavor      = &
&   "abinit"
character(len=*),parameter :: trio_flavor       = &
&   "none"

! Info on experimental features
character(len=*),parameter :: enable_bindings = &
&   "@enable_bindings@"
character(len=*),parameter :: enable_exports = &
&   "no"
character(len=*),parameter :: enable_gw_dpc = &
&   "no"

#if defined HAVE_BZR_BRANCH
! Info on Bazaar branch (if applies)
character(len=*),parameter :: bzr_branch = &
&   ""
character(len=*),parameter :: bzr_revno  = &
&   ""
character(len=*),parameter :: bzr_clean  = &
&   ""
#endif

contains  !===========================================================
!!***

!!****f* ABINIT/m_build_info/dump_config
!! NAME
!!  dump_config
!!
!! FUNCTION
!!  Reports a printout of the information stored in m_build_info,
!!  useful for error messages and debugging.
!!
!! INPUTS
!!  my_unit= Fortran unit number
!!
!! OUTPUT
!!  Only printing
!!
!! PARENTS
!!  Will be filled automatically by the parent script
!!
!! CHILDREN
!!  Will be filled automatically by the parent script
!!
!! SOURCE

subroutine dump_config(my_unit)

implicit none

!Arguments ------------------------------------
integer,intent(in) :: my_unit

!Local variables-------------------------------

! *********************************************************************

! TODO: things that might be added through preprocessing options, e.g.
! date and time of compilation 

write(my_unit,*)
write(my_unit,'(1x,a)') repeat('+',78)
write(my_unit,*)
write(my_unit,'(a)' )' === Build Information === '
write(my_unit,'(2a)')'  Version       : ',trim(abinit_version)
write(my_unit,'(2a)')'  Build target  : ',trim(build_target)
write(my_unit,'(2a)')'  Build date    : ',trim(version_build)
write(my_unit,*)
write(my_unit,'(a)' )' === Compiler Suite === '
write(my_unit,'(2a)')'  C compiler       : ',trim(cc_info)
write(my_unit,'(2a)')'  C++ compiler     : ',trim(cxx_info)
write(my_unit,'(2a)')'  Fortran compiler : ',trim(fc_info)
write(my_unit,'(2a)')'  CFLAGS           : ',trim(cc_flags)
write(my_unit,'(2a)')'  CXXFLAGS         : ',trim(cxx_flags)
write(my_unit,'(2a)')'  FCFLAGS          : ',trim(fc_flags)
write(my_unit,'(2a)')'  FC_LDFLAGS       : ',trim(fc_ldflags)
write(my_unit,*)
write(my_unit,'(a) ')' === Optimizations === '
write(my_unit,'(2a)')'  Debug level        : ',trim(enable_debug)
write(my_unit,'(2a)')'  Optimization level : ',trim(enable_optim)
write(my_unit,'(2a)')'  Architecture       : ',trim(cpu_info)
write(my_unit,*)
write(my_unit,'(a) ')' === Multicore === '
write(my_unit,'(2a)')'  Parallel build : ',trim(enable_mpi)
write(my_unit,'(2a)')'  Parallel I/O   : ',trim(enable_mpi_io)
write(my_unit,'(2a)')'  openMP support : ',trim(enable_openmp)
write(my_unit,'(2a)')'  GPU support    : ',trim(enable_gpu)
write(my_unit,*)
write(my_unit,'(a) ')' === Connectors / Fallbacks === '
write(my_unit,'(2a)')'  Connectors on : ',trim(enable_connectors)
write(my_unit,'(2a)')'  Fallbacks on  : ',trim(enable_fallbacks)
write(my_unit,'(2a)')'  DFT flavor    : ',trim(dft_flavor)
write(my_unit,'(2a)')'  FFT flavor    : ',trim(fft_flavor)
write(my_unit,'(2a)')'  LINALG flavor : ',trim(linalg_flavor)
write(my_unit,'(2a)')'  MATH flavor   : ',trim(math_flavor)
write(my_unit,'(2a)')'  TIMER flavor  : ',trim(timer_flavor)
write(my_unit,'(2a)')'  TRIO flavor   : ',trim(trio_flavor)
write(my_unit,*)
write(my_unit,'(a)' )' === Experimental features === '
write(my_unit,'(2a)')'  Bindings            : ',trim(enable_bindings)
write(my_unit,'(2a)')'  Exports             : ',trim(enable_exports)
write(my_unit,'(2a)')'  GW double-precision : ',trim(enable_gw_dpc)
write(my_unit,*)
#if defined HAVE_BZR_BRANCH
write(my_unit,'(a)' )' === Bazaar branch information === '
write(my_unit,'(2a)')'  Branch ID : ',trim(bzr_branch)
write(my_unit,'(2a)')'  Revision  : ',trim(bzr_revno)
write(my_unit,'(2a)')'  Committed : ',trim(bzr_clean)
write(my_unit,*)
#endif
write(my_unit,'(1x,a)') repeat('+',78)
write(my_unit,*)

end subroutine dump_config
!!***

end module m_build_info
!!***
