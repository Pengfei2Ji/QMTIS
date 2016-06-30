!!****m* ABINIT/interfaces_66_wfs
!! NAME
!! interfaces_66_wfs
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/66_wfs
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

module interfaces_66_wfs

 implicit none

interface
 subroutine chebfi(cg,dtset,eig,enl,gs_hamk,gsc,kinpw,mpi_enreg,nband,npw,nspinor,prtvol,resid)
  use defs_basis
  use defs_abitypes
  use m_hamiltonian
  implicit none
  integer,intent(in) :: nband
  integer,intent(in) :: npw
  integer,intent(in) :: nspinor
  integer,intent(in) :: prtvol
  type(dataset_type),intent(in) :: dtset
  type(gs_hamiltonian_type),intent(inout) :: gs_hamk
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(inout), target :: cg(2,npw*nspinor*nband)
  real(dp),intent(out) :: eig(nband)
  real(dp),intent(out) :: enl(nband*(1-gs_hamk%usepaw))
  real(dp),intent(inout), target :: gsc(2,npw*nspinor*nband)
  real(dp),intent(in) :: kinpw(npw)
  real(dp),intent(out) :: resid(nband)
 end subroutine chebfi
end interface

interface
 function cheb_poly(x, n, a, b) result(y)
  use defs_basis
  implicit none
  integer, intent(in) :: n
  real(dp), intent(in) :: a
  real(dp), intent(in) :: b
  real(dp), intent(in) :: x
  real(dp) :: y
 end function cheb_poly
end interface

interface
 function cheb_oracle(x, a, b, tol, nmax) result(n)
  use defs_basis
  implicit none
  integer :: n
  integer :: nmax
  real(dp), intent(in) :: a
  real(dp), intent(in) :: b
  real(dp) :: tol
  real(dp), intent(in) :: x
 end function cheb_oracle
end interface

interface
 subroutine fock_getghc(cwavef,cwaveprj,ghc,gs_ham,mpi_enreg)
  use defs_basis
  use defs_abitypes
  use m_pawcprj
  use m_hamiltonian
  implicit none
  type(gs_hamiltonian_type),target,intent(inout) :: gs_ham
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(inout) :: cwavef(2,gs_ham%npw_k)
  type(pawcprj_type),intent(in) :: cwaveprj(:,:)
  real(dp),intent(inout) :: ghc(2,gs_ham%npw_k)
 end subroutine fock_getghc
end interface

interface
 subroutine fxphas(cg,gsc,icg,igsc,istwfk,mcg,mgsc,mpi_enreg,nband_k,npw_k,useoverlap)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: icg
  integer,intent(in) :: igsc
  integer,intent(in) :: istwfk
  integer,intent(in) :: mcg
  integer,intent(in) :: mgsc
  integer,intent(in) :: nband_k
  integer,intent(in) :: npw_k
  integer,intent(in) :: useoverlap
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(inout) :: cg(2,mcg)
  real(dp),intent(inout) :: gsc(2,mgsc*useoverlap)
 end subroutine fxphas
end interface

interface
 subroutine fxphas_seq(cg,gsc,icg,igsc,istwfk,mcg,mgsc,nband_k,npw_k,useoverlap)
  use defs_basis
  implicit none
  integer,intent(in) :: icg
  integer,intent(in) :: igsc
  integer,intent(in) :: istwfk
  integer,intent(in) :: mcg
  integer,intent(in) :: mgsc
  integer,intent(in) :: nband_k
  integer,intent(in) :: npw_k
  integer,intent(in) :: useoverlap
  real(dp),intent(inout) :: cg(2,mcg)
  real(dp),intent(inout) :: gsc(2,mgsc*useoverlap)
 end subroutine fxphas_seq
end interface

interface
 subroutine getdc1(cgq,cprjq,dcwavef,dcwaveprj,ibgq,icgq,istwfk,mcgq,mcprjq,&  
  &  mpi_enreg,natom,nband,npw1,nspinor,optcprj,s1cwave0)
  use defs_basis
  use defs_abitypes
  use m_pawcprj
  implicit none
  integer,intent(in) :: ibgq
  integer,intent(in) :: icgq
  integer,intent(in) :: istwfk
  integer,intent(in) :: mcgq
  integer,intent(in) :: mcprjq
  integer,intent(in) :: natom
  integer,intent(in) :: nband
  integer,intent(in) :: npw1
  integer,intent(in) :: nspinor
  integer,intent(in) :: optcprj
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: cgq(2,mcgq)
  type(pawcprj_type),intent(in) :: cprjq(natom,mcprjq)
  real(dp),intent(out) :: dcwavef(2,npw1*nspinor)
  type(pawcprj_type),intent(inout) :: dcwaveprj(natom,nspinor*optcprj)
  real(dp),intent(in) :: s1cwave0(2,npw1*nspinor)
 end subroutine getdc1
end interface

interface
 subroutine getgh1c(berryopt,cwave,cwaveprj,gh1c,grad_berry,gs1c,gs_hamkq,&  
  &  gvnl1,idir,ipert,lambda,mpi_enreg,optlocal,optnl,opt_gvnl1,&  
  &  rf_hamkq,sij_opt,tim_getgh1c,usevnl)
  use defs_basis
  use defs_abitypes
  use m_pawcprj
  use m_hamiltonian
  implicit none
  integer,intent(in) :: berryopt
  integer,intent(in) :: idir
  integer,intent(in) :: ipert
  integer,intent(in) :: opt_gvnl1
  integer,intent(in) :: optlocal
  integer,intent(in) :: optnl
  integer,intent(in) :: sij_opt
  integer,intent(in) :: tim_getgh1c
  integer,intent(in) :: usevnl
  type(gs_hamiltonian_type),intent(inout),target :: gs_hamkq
  real(dp),intent(in) :: lambda
  type(mpi_type),intent(inout) :: mpi_enreg
  type(rf_hamiltonian_type),intent(inout),target :: rf_hamkq
  real(dp),intent(inout) :: cwave(2,gs_hamkq%npw_k*gs_hamkq%nspinor)
  type(pawcprj_type),intent(inout),target :: cwaveprj(:,:)
  real(dp),intent(out) :: gh1c(2,gs_hamkq%npw_kp*gs_hamkq%nspinor)
  real(dp),intent(in) :: grad_berry(:,:)
  real(dp),intent(out) :: gs1c(2,gs_hamkq%npw_kp*gs_hamkq%nspinor)
  real(dp),intent(inout),target :: gvnl1(2,gs_hamkq%npw_kp*gs_hamkq%nspinor)
 end subroutine getgh1c
end interface

interface
 subroutine rf_transgrid_and_pack(isppol,nspden,usepaw,cplex,nfftf,nfft,ngfft,nvloc,&  
  &  pawfgr,mpi_enreg,vtrial,vtrial1,vlocal,vlocal1)
  use m_pawfgr
  use defs_abitypes
  use defs_basis
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: isppol
  integer,intent(in) :: nfft
  integer,intent(in) :: nfftf
  integer,intent(in) :: nspden
  integer,intent(in) :: nvloc
  integer,intent(in) :: usepaw
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawfgr_type),intent(in) :: pawfgr
  integer,intent(in) :: ngfft(18)
  real(dp),intent(out) :: vlocal(ngfft(4),ngfft(5),ngfft(6),nvloc)
  real(dp),intent(out) :: vlocal1(cplex*ngfft(4),ngfft(5),ngfft(6),nvloc)
  real(dp),intent(in),target :: vtrial(nfftf,nspden)
  real(dp),intent(inout),target :: vtrial1(cplex*nfftf,nspden)
 end subroutine rf_transgrid_and_pack
end interface

interface
 subroutine getgh1c_setup(gs_hamkq,rf_hamkq,dtset,psps,kpoint,kpq,idir,ipert,&  ! In
  &  natom,rmet,gprimd,gmet,istwf_k,npw_k,npw1_k,&  ! In
  &  useylmgr1,kg_k,ylm_k,kg1_k,ylm1_k,ylmgr1_k,&  ! In
  &  dkinpw,nkpg,nkpg1,kpg_k,kpg1_k,kinpw1,ffnlk,ffnl1,ph3d,ph3d1,&  ! Out
  &  ddkinpw,dkinpw2,rf_hamk_dir2)                                          ! Optional
  use m_hamiltonian
  use defs_abitypes
  use defs_datatypes
  use defs_basis
  implicit none
  integer,intent(in) :: idir
  integer,intent(in) :: ipert
  integer,intent(in) :: istwf_k
  integer,intent(in) :: natom
  integer,intent(out) :: nkpg
  integer,intent(out) :: nkpg1
  integer,intent(in) :: npw1_k
  integer,intent(in) :: npw_k
  integer,intent(in) :: useylmgr1
  type(dataset_type),intent(in) :: dtset
  type(gs_hamiltonian_type),intent(inout) :: gs_hamkq
  type(pseudopotential_type),intent(in) :: psps
  type(rf_hamiltonian_type),intent(inout),optional :: rf_hamk_dir2
  type(rf_hamiltonian_type),intent(inout) :: rf_hamkq
  real(dp),allocatable,intent(out),optional :: ddkinpw(:)
  real(dp),allocatable,intent(out) :: dkinpw(:)
  real(dp),allocatable,intent(out),optional :: dkinpw2(:)
  real(dp),allocatable,intent(out) :: ffnl1(:,:,:,:)
  real(dp),allocatable,intent(out) :: ffnlk(:,:,:,:)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: kg1_k(3,npw1_k)
  integer,intent(in) :: kg_k(3,npw_k)
  real(dp),allocatable,intent(out) :: kinpw1(:)
  real(dp),allocatable,intent(out) :: kpg1_k(:,:)
  real(dp),allocatable,intent(out) :: kpg_k(:,:)
  real(dp),intent(in) :: kpoint(3)
  real(dp),intent(in) :: kpq(3)
  real(dp),allocatable,intent(out) :: ph3d(:,:,:)
  real(dp),allocatable,intent(out) :: ph3d1(:,:,:)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: ylm1_k(npw1_k,psps%mpsang*psps%mpsang*psps%useylm)
  real(dp),intent(in) :: ylm_k(npw_k,psps%mpsang*psps%mpsang*psps%useylm)
  real(dp),intent(in) :: ylmgr1_k(npw1_k,3+6*((ipert-natom)/10), &
  &         psps%mpsang*psps%mpsang*psps%useylm*useylmgr1)
 end subroutine getgh1c_setup
end interface

interface
 subroutine getgh2c(cwavef,cwaveprj,gh2c,gs2c,gs_hamkq,gvnl2,idir,ipert,lambda,&  
  &  mpi_enreg,optlocal,optnl,rf_hamkq,sij_opt,tim_getgh2c,usevnl)
  use defs_basis
  use defs_abitypes
  use m_pawcprj
  use m_hamiltonian
  implicit none
  integer,intent(in) :: idir
  integer,intent(in) :: ipert
  integer,intent(in) :: optlocal
  integer,intent(in) :: optnl
  integer,intent(in) :: sij_opt
  integer,intent(in) :: tim_getgh2c
  integer,intent(in) :: usevnl
  type(gs_hamiltonian_type),intent(inout),target :: gs_hamkq
  real(dp),intent(in) :: lambda
  type(mpi_type),intent(inout) :: mpi_enreg
  type(rf_hamiltonian_type),intent(inout),target :: rf_hamkq
  real(dp),intent(inout) :: cwavef(:,:)
  type(pawcprj_type),intent(inout),target :: cwaveprj(:,:)
  real(dp),intent(out) :: gh2c(:,:)
  real(dp),intent(out) :: gs2c(:,:)
  real(dp),intent(inout),target :: gvnl2(:,:)
 end subroutine getgh2c
end interface

interface
 subroutine getghc(cpopt,cwavef,cwaveprj,ghc,gsc,gs_ham,gvnlc,lambda,mpi_enreg,ndat,&  
  &  prtvol,sij_opt,tim_getghc,type_calc,&  
  &  kg_fft_k,kg_fft_kp,select_k) ! optional arguments
  use defs_basis
  use defs_abitypes
  use m_pawcprj
  use m_hamiltonian
  implicit none
  integer,intent(in) :: cpopt
  integer,intent(in) :: ndat
  integer,intent(in) :: prtvol
  integer,intent(in),optional :: select_k
  integer,intent(in) :: sij_opt
  integer,intent(in) :: tim_getghc
  integer,intent(in) :: type_calc
  type(gs_hamiltonian_type),intent(inout),target :: gs_ham
  real(dp),intent(in) :: lambda
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in),optional,target :: kg_fft_k(:,:)
  integer,intent(in),optional,target :: kg_fft_kp(:,:)
  real(dp),intent(inout) :: cwavef(:,:)
  type(pawcprj_type),intent(inout),target :: cwaveprj(:,:)
  real(dp),intent(inout) :: ghc(:,:)
  real(dp),intent(out),target :: gsc(:,:)
  real(dp),intent(inout) :: gvnlc(:,:)
 end subroutine getghc
end interface

interface
 subroutine getghc_mGGA(cwavef,ghc_mGGA,gbound_k,gprimd,istwf_k,kg_k,kpt,mgfft,mpi_enreg,&  
  &  ndat,ngfft,npw_k,nvloc,n4,n5,n6,my_nspinor,vxctaulocal,use_gpu_cuda)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: istwf_k
  integer,intent(in) :: mgfft
  integer,intent(in) :: my_nspinor
  integer,intent(in) :: n4
  integer,intent(in) :: n5
  integer,intent(in) :: n6
  integer,intent(in) :: ndat
  integer,intent(in) :: npw_k
  integer,intent(in) :: nvloc
  integer,intent(in) :: use_gpu_cuda
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(inout) :: cwavef(2,npw_k*my_nspinor*ndat)
  integer,intent(in) :: gbound_k(2*mgfft+4)
  real(dp),intent(inout) :: ghc_mGGA(2,npw_k*my_nspinor*ndat)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: kg_k(3,npw_k)
  real(dp),intent(in) :: kpt(3)
  real(dp),intent(inout) :: vxctaulocal(n4,n5,n6,nvloc,4)
 end subroutine getghc_mGGA
end interface

interface
 subroutine getghcnd(cwavef,ghcnd,gs_ham,my_nspinor,ndat)
  use defs_basis
  use m_hamiltonian
  implicit none
  integer,intent(in) :: my_nspinor
  integer,intent(in) :: ndat
  type(gs_hamiltonian_type),intent(in),target :: gs_ham
  real(dp),intent(in) :: cwavef(2,gs_ham%npw_k*my_nspinor*ndat)
  real(dp),intent(out) :: ghcnd(2,gs_ham%npw_k*my_nspinor*ndat)
 end subroutine getghcnd
end interface

interface
 subroutine getgsc(cg,cprj,gs_ham,gsc,ibg,icg,igsc,ikpt,isppol,&  
  &  mcg,mcprj,mgsc,mpi_enreg,natom,nband,npw_k,nspinor,select_k)
  use defs_basis
  use defs_abitypes
  use m_pawcprj
  use m_hamiltonian
  implicit none
  integer,intent(in) :: ibg
  integer,intent(in) :: icg
  integer,intent(in) :: igsc
  integer,intent(in) :: ikpt
  integer,intent(in) :: isppol
  integer,intent(in) :: mcg
  integer,intent(in) :: mcprj
  integer,intent(in) :: mgsc
  integer,intent(in) :: natom
  integer,intent(in) :: nband
  integer,intent(in) :: npw_k
  integer,intent(in) :: nspinor
  integer,intent(in),optional :: select_k
  type(gs_hamiltonian_type),intent(inout),target :: gs_ham
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: cg(2,mcg)
  type(pawcprj_type),intent(in) :: cprj(natom,mcprj)
  real(dp),intent(out) :: gsc(2,mgsc)
 end subroutine getgsc
end interface

interface
 subroutine prep_bandfft_tabs(gs_hamk,ikpt,mkmem,mpi_enreg)
  use m_hamiltonian
  use defs_abitypes
  implicit none
  integer,intent(in) :: ikpt
  integer,intent(in) :: mkmem
  type(gs_hamiltonian_type),intent(inout) :: gs_hamk
  type(mpi_type),intent(inout) :: mpi_enreg
 end subroutine prep_bandfft_tabs
end interface

interface
 subroutine prep_fourwf(rhoaug,blocksize,cwavef,wfraug,iblock,istwf_k,mgfft,&  
  &  mpi_enreg,nband_k,ndat,ngfft,npw_k,n4,n5,n6,occ_k,option_fourwf,ucvol,wtk,&  
  &  bandfft_kpt_tab,use_gpu_cuda) ! Optional arguments
  use m_bandfft_kpt
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: blocksize
  integer,intent(in) :: iblock
  integer,intent(in) :: istwf_k
  integer,intent(in) :: mgfft
  integer,intent(in) :: n4
  integer,intent(in) :: n5
  integer,intent(in) :: n6
  integer,intent(in) :: nband_k
  integer,intent(in) :: ndat
  integer,intent(in) :: npw_k
  integer,intent(in) :: option_fourwf
  integer,intent(in),optional :: use_gpu_cuda
  type(bandfft_kpt_type),optional,target,intent(in) :: bandfft_kpt_tab
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: wtk
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: cwavef(2,npw_k*blocksize)
  real(dp),intent(in) :: occ_k(nband_k)
  real(dp),intent(out) :: rhoaug(n4,n5,n6)
  real(dp),target,intent(inout) :: wfraug(2,n4,n5,n6*ndat)
 end subroutine prep_fourwf
end interface

interface
 subroutine prep_getghc(cwavef,gs_hamk,gvnlc,gwavef,swavef,lambda,blocksize,&  
  &  mpi_enreg,prtvol,sij_opt,cpopt,cwaveprj,&  
  &  already_transposed) ! optional argument
  use defs_basis
  use defs_abitypes
  use m_pawcprj
  use m_hamiltonian
  implicit none
  integer,intent(in) :: blocksize
  integer,intent(in) :: cpopt
  integer,intent(in) :: prtvol
  integer,intent(in) :: sij_opt
  logical, intent(in),optional :: already_transposed
  type(gs_hamiltonian_type),intent(inout) :: gs_hamk
  real(dp),intent(in) :: lambda
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(inout) :: cwavef(:,:)
  type(pawcprj_type), intent(inout) :: cwaveprj(:,:)
  real(dp),intent(inout) :: gvnlc(:,:)
  real(dp),intent(inout) :: gwavef(:,:)
  real(dp),intent(inout) :: swavef(:,:)
 end subroutine prep_getghc
end interface

interface
 subroutine prep_index_wavef_bandpp(nproc_band,bandpp,&  
  nspinor,ndatarecv,&  
  recvcounts,rdispls,&  
  index_wavef_band)
  implicit none
  integer,intent(in) :: bandpp
  integer,intent(in) :: ndatarecv
  integer,intent(in) :: nproc_band
  integer,intent(in) :: nspinor
  integer,allocatable,intent(out) :: index_wavef_band(:)
  integer,intent(in) :: rdispls(nproc_band)
  integer,intent(in) :: recvcounts(nproc_band)
 end subroutine prep_index_wavef_bandpp
end interface

interface
 subroutine prep_nonlop(choice,cpopt,cwaveprj,enlout_block,hamk,idir,lambdablock,&  
  &  blocksize,mpi_enreg,nnlout,paw_opt,signs,gsc,&  
  &  tim_nonlop,cwavef,gvnlc,already_transposed)
  use defs_basis
  use defs_abitypes
  use m_pawcprj
  use m_hamiltonian
  implicit none
  integer,intent(in) :: blocksize
  integer,intent(in) :: choice
  integer,intent(in) :: cpopt
  integer,intent(in) :: idir
  integer,intent(in) :: nnlout
  integer,intent(in) :: paw_opt
  integer,intent(in) :: signs
  integer :: tim_nonlop
  logical,optional,intent(in) :: already_transposed
  type(gs_hamiltonian_type),intent(in) :: hamk
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(inout) :: cwavef(:,:)
  type(pawcprj_type),intent(inout) :: cwaveprj(:,:)
  real(dp),intent(out) :: enlout_block(nnlout*blocksize)
  real(dp),intent(out) :: gsc(:,:)
  real(dp),intent(out) :: gvnlc(:,:)
  real(dp),intent(in) :: lambdablock(blocksize)
 end subroutine prep_nonlop
end interface

interface
 subroutine prep_sort_wavef_spin(nproc_band,nspinor,ndatarecv,recvcounts,rdispls,index_wavef)
  implicit none
  integer,intent(in) :: ndatarecv
  integer,intent(in) :: nproc_band
  integer,intent(in) :: nspinor
  integer,allocatable,intent(out) :: index_wavef(:)
  integer,intent(in) :: rdispls(nproc_band)
  integer,intent(in) :: recvcounts(nproc_band)
 end subroutine prep_sort_wavef_spin
end interface

interface
 subroutine prep_wavef_sym_do(mpi_enreg,bandpp,nspinor,&  
  &  ndatarecv,&  
  &  ndatarecv_tot,ndatasend_sym,tab_proc,&  
  &  cwavef_alltoall,&  
  &  sendcounts_sym,sdispls_sym,&  
  &  recvcounts_sym,rdispls_sym,&  
  &  ewavef_alltoall_sym,&  
  &  index_wavef_send)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: bandpp
  integer,intent(in) :: ndatarecv
  integer,intent(in) :: ndatarecv_tot
  integer,intent(in) :: ndatasend_sym
  integer,intent(in) :: nspinor
  type(mpi_type),intent(in) :: mpi_enreg
  integer,allocatable,intent(out) :: index_wavef_send(:)
  integer,intent(in) :: rdispls_sym(:)
  integer,intent(in) :: recvcounts_sym(:)
  integer,intent(in) :: sdispls_sym(:)
  integer,intent(in) :: sendcounts_sym(:)
  integer,intent(in) :: tab_proc(:)
  real(dp),intent(inout) :: cwavef_alltoall(2,ndatarecv*nspinor*bandpp)
  real(dp),pointer :: ewavef_alltoall_sym(:,:)
 end subroutine prep_wavef_sym_do
end interface

interface
 subroutine prep_wavef_sym_undo(mpi_enreg,bandpp,nspinor,&  
  &  ndatarecv,&  
  &  ndatarecv_tot,ndatasend_sym,idatarecv0,&  
  &  gwavef_alltoall,&  
  &  sendcounts_sym,sdispls_sym,&  
  &  recvcounts_sym,rdispls_sym,&  
  &  gwavef_alltoall_sym,&  
  &  index_wavef_send)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: bandpp
  integer,intent(in) :: idatarecv0
  integer,intent(in) :: ndatarecv
  integer,intent(in) :: ndatarecv_tot
  integer,intent(in) :: ndatasend_sym
  integer,intent(in) :: nspinor
  type(mpi_type),intent(in) :: mpi_enreg
  integer,intent(in) :: index_wavef_send(:)
  integer,intent(in) :: rdispls_sym(:)
  integer,intent(in) :: recvcounts_sym(:)
  integer,intent(in) :: sdispls_sym(:)
  integer,intent(in) :: sendcounts_sym(:)
  real(dp),intent(inout) :: gwavef_alltoall(2,ndatarecv*nspinor*bandpp)
  real(dp),intent(inout) :: gwavef_alltoall_sym(:,:)
 end subroutine prep_wavef_sym_undo
end interface

interface
 subroutine pw_orthon(icg,igsc,istwf_k,mcg,mgsc,nelem,nvec,ortalgo,ovl_vecnm,useoverlap,vecnm,me_g0,comm)
  use defs_basis
  implicit none
  integer,intent(in) :: comm
  integer,intent(in) :: icg
  integer,intent(in) :: igsc
  integer,intent(in) :: istwf_k
  integer,intent(in) :: mcg
  integer,intent(in) :: me_g0
  integer,intent(in) :: mgsc
  integer,intent(in) :: nelem
  integer,intent(in) :: nvec
  integer,intent(in) :: ortalgo
  integer,intent(in) :: useoverlap
  real(dp),intent(inout) :: ovl_vecnm(2,mgsc*useoverlap)
  real(dp),intent(inout) :: vecnm(2,mcg)
 end subroutine pw_orthon
end interface

interface
 subroutine rayleigh_ritz_subdiago(cg,ghc,gsc,gvnlc,eig,istwf_k,mpi_enreg,nband,npw,nspinor,usepaw)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: istwf_k
  integer,intent(in) :: nband
  integer,intent(in) :: npw
  integer,intent(in) :: nspinor
  integer,intent(in) :: usepaw
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(inout) :: cg(2,npw*nspinor*nband)
  real(dp),intent(out) :: eig(nband)
  real(dp),intent(inout) :: ghc(2,npw*nspinor*nband)
  real(dp),intent(inout) :: gsc(2,npw*nspinor*nband)
  real(dp),intent(inout) :: gvnlc(2,npw*nspinor*nband)
 end subroutine rayleigh_ritz_subdiago
end interface

interface
 subroutine rayleigh_ritz_distributed(cg,ghc,gsc,gvnlc,eig,istwf_k,mpi_enreg,nband,npw,nspinor,usepaw)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: istwf_k
  integer,intent(in) :: nband
  integer,intent(in) :: npw
  integer,intent(in) :: nspinor
  integer,intent(in) :: usepaw
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(inout) :: cg(2,npw*nspinor*nband)
  real(dp),intent(out) :: eig(nband)
  real(dp),intent(inout) :: ghc(2,npw*nspinor*nband)
  real(dp),intent(inout) :: gsc(2,npw*nspinor*nband)
  real(dp),intent(inout) :: gvnlc(2,npw*nspinor*nband)
 end subroutine rayleigh_ritz_distributed
end interface

interface
 subroutine from_mat_to_block_cyclic(full_mat, vectsize, nband, block_cyclic_mat, buffsize, blocksize, iproc, nprocs)
  use defs_basis
  implicit none
  integer, intent(in) :: blocksize
  integer, intent(in) :: buffsize
  integer, intent(in) :: iproc
  integer, intent(in) :: nband
  integer, intent(in) :: nprocs
  integer, intent(in) :: vectsize
  real(dp), intent(inout) :: block_cyclic_mat(2, vectsize*buffsize)
  real(dp), intent(in) :: full_mat(2, vectsize*nband)
 end subroutine from_mat_to_block_cyclic
end interface

interface
 subroutine from_block_cyclic_to_mat(full_mat, vectsize, nband, block_cyclic_mat, buffsize, blocksize, iproc, nprocs)
  use defs_basis
  implicit none
  integer, intent(in) :: blocksize
  integer, intent(in) :: buffsize
  integer, intent(in) :: iproc
  integer, intent(in) :: nband
  integer, intent(in) :: nprocs
  integer, intent(in) :: vectsize
  real(dp), intent(in) :: block_cyclic_mat(2, vectsize*buffsize)
  real(dp), intent(inout) :: full_mat(2, vectsize*nband)
 end subroutine from_block_cyclic_to_mat
end interface

interface
 subroutine pack_matrix(mat_in, mat_out, N, cplx)
  use defs_basis
  implicit none
  integer, intent(in) :: N
  integer, intent(in) :: cplx
  real(dp), intent(in) :: mat_in(cplx, N*N)
  real(dp), intent(out) :: mat_out(cplx*N*(N+1)/2)
 end subroutine pack_matrix
end interface

interface
 subroutine subdiago(cg,eig_k,evec,gsc,icg,igsc,istwf_k,&  
  &  mcg,mgsc,nband_k,npw_k,nspinor,paral_kgb,&  
  &  subham,subovl,use_subovl,usepaw,me_g0)
  use defs_basis
  implicit none
  integer,intent(in) :: icg
  integer,intent(in) :: igsc
  integer,intent(in) :: istwf_k
  integer,intent(in) :: mcg
  integer,intent(in) :: me_g0
  integer,intent(in) :: mgsc
  integer,intent(in) :: nband_k
  integer,intent(in) :: npw_k
  integer,intent(in) :: nspinor
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: use_subovl
  integer,intent(in) :: usepaw
  real(dp),intent(inout) :: cg(2,mcg)
  real(dp),intent(out) :: eig_k(nband_k)
  real(dp),intent(out) :: evec(2*nband_k,nband_k)
  real(dp),intent(inout) :: gsc(2,mgsc)
  real(dp),intent(inout) :: subham(nband_k*(nband_k+1))
  real(dp),intent(inout) :: subovl(nband_k*(nband_k+1)*use_subovl)
 end subroutine subdiago
end interface

interface
 subroutine wfconv(ceksp2,cg1,cg2,debug,ecut1,ecut2,ecut2_eff,&  
  &  eig_k1,eig_k2,exchn2n3d,formeig,gmet1,gmet2,icg1,icg2,&  
  &  ikpt1,ikpt10,ikpt2,indkk,inplace,isppol2,istwfk1,istwfk2,&  
  &  kg1,kg2,kptns1,kptns2,mband1,mband2,mcg1,mcg2,mpi_enreg1,mpi_enreg2,&  
  &  mpw1,mpw2,nbd1,nbd2,ngfft1,ngfft2,nkpt1,nkpt2,npw1,npw2,nspinor1,nspinor2,&  
  &  nsym,occ_k1,occ_k2,optorth,randalg,restart,rprimd2,sppoldbl,symrel,tnons)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: ceksp2
  integer,intent(in) :: debug
  integer,intent(in) :: exchn2n3d
  integer,intent(in) :: formeig
  integer,intent(in) :: icg1
  integer,intent(in) :: icg2
  integer,intent(in) :: ikpt1
  integer,intent(inout) :: ikpt10
  integer,intent(in) :: ikpt2
  integer,intent(in) :: inplace
  integer,intent(in) :: isppol2
  integer,intent(in) :: mband1
  integer,intent(in) :: mband2
  integer,intent(in) :: mcg1
  integer,intent(in) :: mcg2
  integer,intent(in) :: mpw1
  integer,intent(in) :: mpw2
  integer,intent(in) :: nbd1
  integer,intent(in) :: nbd2
  integer,intent(in) :: nkpt1
  integer,intent(in) :: nkpt2
  integer,intent(inout) :: npw1
  integer,intent(inout) :: npw2
  integer,intent(in) :: nspinor1
  integer,intent(in) :: nspinor2
  integer,intent(in) :: nsym
  integer,intent(in) :: optorth
  integer,intent(in) :: randalg
  integer,intent(in) :: restart
  integer,intent(in) :: sppoldbl
  real(dp),intent(in) :: ecut1
  real(dp),intent(in) :: ecut2
  real(dp),intent(in) :: ecut2_eff
  type(mpi_type),intent(inout) :: mpi_enreg1
  type(mpi_type),intent(inout) :: mpi_enreg2
  integer,intent(in) :: ngfft1(18)
  integer,intent(in) :: ngfft2(18)
  real(dp),intent(inout) :: cg1(2,mcg1)
  real(dp),intent(inout) :: cg2(2,mcg2)
  real(dp),intent(inout) :: eig_k1(mband1*(2*mband1)**formeig)
  real(dp),intent(inout) :: eig_k2(mband2*(2*mband2)**formeig)
  real(dp),intent(in) :: gmet1(3,3)
  real(dp),intent(in) :: gmet2(3,3)
  integer,intent(in) :: indkk(nkpt2*sppoldbl,6)
  integer,intent(in) :: istwfk1(nkpt1)
  integer,intent(in) :: istwfk2(nkpt2)
  integer,intent(inout) :: kg1(3,mpw1)
  integer,intent(inout) :: kg2(3,mpw2)
  real(dp),intent(in) :: kptns1(3,nkpt1)
  real(dp),intent(in) :: kptns2(3,nkpt2)
  real(dp),intent(inout) :: occ_k1(mband1)
  real(dp),intent(inout) :: occ_k2(mband2)
  real(dp),intent(in) :: rprimd2(3,3)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
 end subroutine wfconv
end interface

end module interfaces_66_wfs
!!***
