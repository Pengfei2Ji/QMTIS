!!****m* ABINIT/interfaces_32_util
!! NAME
!! interfaces_32_util
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/32_util
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

module interfaces_32_util

 implicit none

interface
 subroutine acrossb(a,b,c)
  use defs_basis
  implicit none
  real(dp),intent(in) :: a(3)
  real(dp),intent(in) :: b(3)
  real(dp),intent(out) :: c(3)
 end subroutine acrossb
end interface

interface
 subroutine appdig(integ,string,strinn)
  implicit none
  integer,intent(in) :: integ
  character(len=*),intent(in) :: string
  character(len=*),intent(out) :: strinn
 end subroutine appdig
end interface

interface
 subroutine besjm(arg,besjx,cosx,nn,nx,sinx,xx)
  use defs_basis
  implicit none
  integer,intent(in) :: nn
  integer,intent(in) :: nx
  real(dp),intent(in) :: arg
  real(dp),intent(out) :: besjx(nx)
  real(dp),intent(in) :: cosx(nx)
  real(dp),intent(in) :: sinx(nx)
  real(dp),intent(in) :: xx(nx)
 end subroutine besjm
end interface

interface
 subroutine blow_pawuj(mat,nj,matt)
  use defs_basis
  implicit none
  integer,intent(in) :: nj
  real(dp),intent(in) :: mat(nj,nj)
  real(dp),intent(out) :: matt(nj+1,nj+1)
 end subroutine blow_pawuj
end interface

interface
 function clp(x)
  use defs_basis
  implicit none
  real(dp) :: clp
  real(dp),intent(in) :: x
 end function clp
end interface

interface
 subroutine create_mlms2jmj(lcor,mlmstwojmj)
  use defs_basis
  implicit none
  integer,intent(in) :: lcor
  complex(dpc),intent(out) :: mlmstwojmj(2*(2*lcor+1),2*(2*lcor+1))
 end subroutine create_mlms2jmj
end interface

interface
 subroutine create_slm2ylm(lcor,slmtwoylm)
  use defs_basis
  implicit none
  integer,intent(in) :: lcor
  complex(dpc),intent(out) :: slmtwoylm(2*lcor+1,2*lcor+1)
 end subroutine create_slm2ylm
end interface

interface
 subroutine ctrap(imax,ff,hh,ans)
  use defs_basis
  implicit none
  integer,intent(in) :: imax
  real(dp),intent(out) :: ans
  real(dp),intent(in) :: hh
  real(dp),intent(in) :: ff(imax)
 end subroutine ctrap
end interface

interface
 subroutine ddb_chkname(nmfond,nmxpct,nmxpct2)
  implicit none
  character(len=*),intent(in) :: nmfond
  character(len=*),intent(in) :: nmxpct
  character(len=*),intent(in),optional :: nmxpct2
 end subroutine ddb_chkname
end interface

interface
 subroutine fappnd(filapp,filnam,iapp,&  
  &  suff) ! optional argument
  use defs_basis
  implicit none
  integer,intent(in) :: iapp
  character(len=fnlen),intent(out) :: filapp
  character(len=fnlen),intent(in) :: filnam
  character(len=3),optional,intent(in) :: suff
 end subroutine fappnd
end interface

interface
 subroutine findmin(dedv_1,dedv_2,dedv_predict,&  
  &  d2edv2_1,d2edv2_2,d2edv2_predict,&  
  &  etotal_1,etotal_2,etotal_predict,&  
  &  lambda_1,lambda_2,lambda_predict,status)
  use defs_basis
  implicit none
  integer,intent(out) :: status
  real(dp),intent(out) :: d2edv2_1
  real(dp),intent(out) :: d2edv2_2
  real(dp),intent(out) :: d2edv2_predict
  real(dp),intent(in) :: dedv_1
  real(dp),intent(in) :: dedv_2
  real(dp),intent(out) :: dedv_predict
  real(dp),intent(in) :: etotal_1
  real(dp),intent(in) :: etotal_2
  real(dp),intent(out) :: etotal_predict
  real(dp),intent(in) :: lambda_1
  real(dp),intent(in) :: lambda_2
  real(dp),intent(out) :: lambda_predict
 end subroutine findmin
end interface

interface
 subroutine hermit(chmin,chmout,ierr,ndim)
  use defs_basis
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: ndim
  real(dp),intent(inout) :: chmin(ndim*ndim+ndim)
  real(dp),intent(inout) :: chmout(ndim*ndim+ndim)
 end subroutine hermit
end interface

interface
 subroutine inupper(string)
  implicit none
  character(len=*),intent(inout) :: string
 end subroutine inupper
end interface

interface
 subroutine isfile(filnam,status)
  use defs_basis
  implicit none
  character(len=fnlen),intent(inout) :: filnam
  character(len=3),intent(in) :: status
 end subroutine isfile
end interface

interface
 subroutine kramerskronig(nomega,omega,eps,method,only_check)
  use defs_basis
  implicit none
  integer,intent(in) :: method
  integer,intent(in) :: nomega
  integer,intent(in) :: only_check
  complex,intent(inout) :: eps(nomega)
  real(dp),intent(in) :: omega(nomega)
 end subroutine kramerskronig
end interface

interface
 subroutine littlegroup_q(nsym,qpt,symq,symrec,symafm,timrev,prtvol,use_sym)
  use defs_basis
  implicit none
  integer,intent(in) :: nsym
  integer,intent(in),optional :: prtvol
  integer,intent(out) :: timrev
  integer,intent(in),optional :: use_sym
  real(dp),intent(in) :: qpt(3)
  integer,intent(in) :: symafm(nsym)
  integer,intent(out) :: symq(4,2,nsym)
  integer,intent(in) :: symrec(3,3,nsym)
 end subroutine littlegroup_q
end interface

interface
 subroutine mati3det(mm,det)
  implicit none
  integer,intent(out) :: det
  integer,intent(in) :: mm(3,3)
 end subroutine mati3det
end interface

interface
 subroutine mati3inv(mm,mit)
  implicit none
  integer,intent(out) :: mit(3,3)
  integer,intent(in) :: mm(3,3)
 end subroutine mati3inv
end interface

interface
 subroutine matr3eigval(eigval,matr)
  use defs_basis
  implicit none
  real(dp),intent(out) :: eigval(3)
  real(dp),intent(in) :: matr(3,3)
 end subroutine matr3eigval
end interface

interface
 subroutine matr3inv(aa,ait)
  use defs_basis
  implicit none
  real(dp),intent(in) :: aa(3,3)
  real(dp),intent(out) :: ait(3,3)
 end subroutine matr3inv
end interface

interface
 subroutine matrginv(a,lda,n)
  use defs_basis
  implicit none
  integer,intent(in) :: lda
  integer,intent(in) :: n
  real(dp),intent(inout) :: a(lda,n)
 end subroutine matrginv
end interface

interface
 subroutine mknormpath(nbounds,bounds,gmet,ndiv_small,ndiv,npt_tot,path)
  use defs_basis
  implicit none
  integer,intent(in) :: nbounds
  integer,intent(in) :: ndiv_small
  integer,intent(inout) :: npt_tot
  real(dp),intent(in) :: bounds(3,nbounds)
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(inout) :: ndiv(nbounds-1)
  real(dp),intent(out),optional :: path(3,npt_tot)
 end subroutine mknormpath
end interface

interface
 subroutine nderiv(hh,yy,zz,ndim,norder)
  use defs_basis
  implicit none
  integer,intent(in) :: ndim
  integer,intent(in) :: norder
  real(dp),intent(in) :: hh
  real(dp),intent(in) :: yy(ndim)
  real(dp),intent(out) :: zz(ndim)
 end subroutine nderiv
end interface

interface
 subroutine overlap_g(doti,dotr,mpw,npw_k1,npw_k2,nspinor,pwind_k,vect1,vect2)
  use defs_basis
  implicit none
  integer,intent(in) :: mpw
  integer,intent(in) :: npw_k1
  integer,intent(in) :: npw_k2
  integer,intent(in) :: nspinor
  real(dp),intent(out) :: doti
  real(dp),intent(out) :: dotr
  integer,intent(in) :: pwind_k(mpw)
  real(dp),intent(in) :: vect1(1:2,0:mpw*nspinor)
  real(dp),intent(in) :: vect2(1:2,0:mpw*nspinor)
 end subroutine overlap_g
end interface

interface
 function permutations(nn,kk)
  use defs_basis
  implicit none
  integer,intent(in) :: kk
  integer,intent(in) :: nn
  real(dp) :: permutations
 end function permutations
end interface

interface
 function phim(costheta,sintheta,mm)
  use defs_basis
  implicit none
  integer,intent(in) :: mm
  real(dp),intent(in) :: costheta
  real(dp) :: phim
  real(dp),intent(in) :: sintheta
 end function phim
end interface

interface
 function proc_distrb_cycle(distrb,ikpt,iband1,iband2,isppol,me) 
  implicit none
  integer,intent(in) :: iband1
  integer,intent(in) :: iband2
  integer,intent(in) :: ikpt
  integer,intent(in) :: isppol
  integer,intent(in) :: me
  logical :: proc_distrb_cycle
  integer,allocatable,intent(in) :: distrb(:,:,:)
 end function proc_distrb_cycle
end interface

interface
 subroutine radsintr(funr,funq,mqgrid,mrgrid,qgrid,rgrid,yq1,yqn)
  use defs_basis
  implicit none
  integer , intent(in) :: mqgrid
  integer , intent(in) :: mrgrid
  real(dp), intent(out) :: yq1
  real(dp), intent(out) :: yqn
  real(dp), intent(out) :: funq(mqgrid)
  real(dp), intent(in) :: funr(mrgrid)
  real(dp), intent(in) :: qgrid(mqgrid)
  real(dp), intent(in) :: rgrid(mrgrid)
 end subroutine radsintr
end interface

interface
 function radsmear(r, rsph, rsm)
  use defs_basis
  implicit none
  real(dp), intent(in) :: r
  real(dp) :: radsmear
  real(dp), intent(in) :: rsm
  real(dp), intent(in) :: rsph
 end function radsmear
end interface

interface
 subroutine rotmat(xaxis,zaxis,inversion_flag,umat)
  use defs_basis
  implicit none
  integer,intent(out) :: inversion_flag
  real(dp),intent(out) :: umat(3,3)
  real(dp),intent(in) :: xaxis(3)
  real(dp),intent(in) :: zaxis(3)
 end subroutine rotmat
end interface

interface
 subroutine sbf8(nm,xx,sb_out)
  use defs_basis
  implicit none
  integer,intent(in) :: nm
  real(dp),intent(in) :: xx
  real(dp),intent(out) :: sb_out(nm)
 end subroutine sbf8
end interface

interface
 function set_istwfk(kpoint) result(istwfk)
  use defs_basis
  implicit none
  integer :: istwfk
  real(dp),intent(in) :: kpoint(3)
 end function set_istwfk
end interface

interface
 subroutine smatrix(cg,cgq,cg1_k,ddkflag,dtm_k,icg,icg1,itrs,job,maxbd,&  
  &  mcg_k,mcg_q,mcg1_k,minbd,mpw,nband_occ,npw_k1,npw_k2,nspinor,&  
  &  pwind_k,pwnsfac_k,sflag_k,shiftbd,smat_inv,smat_k,smat_k_paw,usepaw)
  use defs_basis
  implicit none
  integer,intent(in) :: ddkflag
  integer,intent(in) :: icg
  integer,intent(in) :: icg1
  integer,intent(in) :: itrs
  integer,intent(in) :: job
  integer,intent(in) :: maxbd
  integer,intent(in) :: mcg1_k
  integer,intent(in) :: mcg_k
  integer,intent(in) :: mcg_q
  integer,intent(in) :: minbd
  integer,intent(in) :: mpw
  integer,intent(in) :: nband_occ
  integer,intent(in) :: npw_k1
  integer,intent(in) :: npw_k2
  integer,intent(in) :: nspinor
  integer,intent(in) :: shiftbd
  integer,intent(in) :: usepaw
  real(dp),intent(in) :: cg(2,mcg_k)
  real(dp),intent(out) :: cg1_k(2,mcg1_k)
  real(dp),intent(in) :: cgq(2,mcg_q)
  real(dp),intent(out) :: dtm_k(2)
  integer,intent(in) :: pwind_k(mpw)
  real(dp),intent(in) :: pwnsfac_k(4,mpw)
  integer,intent(inout) :: sflag_k(nband_occ)
  real(dp),intent(out) :: smat_inv(2,nband_occ,nband_occ)
  real(dp),intent(inout) :: smat_k(2,nband_occ,nband_occ)
  real(dp),intent(in) :: smat_k_paw(2,usepaw*nband_occ,usepaw*nband_occ)
 end subroutine smatrix
end interface

interface
 subroutine status(counter,filstat,istat,level,routine)
  implicit none
  integer,intent(in) :: counter
  integer,intent(in) :: istat
  integer,intent(in) :: level
  character(len=*),intent(in) :: filstat
  character(len=*),intent(in) :: routine
 end subroutine status
end interface


end module interfaces_32_util
!!***
