!!****m* ABINIT/interfaces_65_paw
!! NAME
!! interfaces_65_paw
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/65_paw
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

module interfaces_65_paw

 implicit none

interface
 subroutine atomden(MPI_enreg,natom,ntypat,typat,ngrid,r_vec_grid,rho,a,b,c,atom_pos,&  
  &  natomgr,natomgrmax,atomrgrid,density,prtvol,calctype)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: natomgrmax
  integer,intent(in) :: ngrid
  integer,intent(in) :: ntypat
  integer,intent(in) :: prtvol
  type(mpi_type),intent(in) :: MPI_enreg
  character(len=7),intent(in) :: calctype
  real(dp),intent(in) :: a(3)
  real(dp),intent(in) :: atom_pos(3,natom)
  real(dp),intent(in) :: atomrgrid(natomgrmax,ntypat)
  real(dp),intent(in) :: b(3)
  real(dp),intent(in) :: c(3)
  real(dp),intent(in) :: density(natomgrmax,ntypat)
  integer,intent(in) :: natomgr(ntypat)
  real(dp),intent(in) :: r_vec_grid(3,ngrid)
  real(dp),intent(inout) :: rho(ngrid)
  integer,intent(in) :: typat(natom)
 end subroutine atomden
end interface

interface
 subroutine calc_ubare(itypatcor,lpawu,pawang,pawrad,pawtab,rmax)
  use m_pawtab
  use m_pawang
  use m_pawrad
  use defs_basis
  implicit none
  integer, intent(in) :: itypatcor
  integer, intent(in) :: lpawu
  type(pawang_type),intent(in) :: pawang
  type(pawrad_type),intent(in) :: pawrad
  type(pawtab_type),target,intent(in) :: pawtab
  real(dp), optional, intent(in) :: rmax
 end subroutine calc_ubare
end interface

interface
 subroutine chkpawovlp(natom,ntypat,pawovlp,pawtab,rmet,typat,xred)
  use defs_basis
  use m_pawtab
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  real(dp) :: pawovlp
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(in) :: rmet(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine chkpawovlp
end interface

interface
 subroutine convert_notation(mu4,eps_alpha,eps_beta,eps_gamma,eps_delta)
  implicit none
  integer,intent(in) :: eps_alpha
  integer,intent(in) :: eps_beta
  integer,optional,intent(in) :: eps_delta
  integer,optional,intent(in) :: eps_gamma
  integer,intent(inout) :: mu4(4)
 end subroutine convert_notation
end interface

interface
 subroutine denfgr(atindx1,gmet,spaceComm_in,my_natom,natom,nattyp,ngfft,nhat,nspinor,nsppol,nspden,ntypat,&  
  &  pawfgr,pawrad,pawrhoij,pawtab,prtvol,rhor,rhor_paw,rhor_n_one,rhor_nt_one,rprimd,typat,ucvol,xred,&  
  &  abs_n_tilde_nt_diff,znucl,mpi_atmtab,comm_atom) ! Optional arguments
  use defs_basis
  use m_pawtab
  use m_pawrad
  use m_pawrhoij
  use m_pawfgr
  implicit none
  integer,optional,intent(in) :: comm_atom
  integer,intent(in) :: my_natom
  integer,intent(in) :: natom
  integer,intent(in) :: nspden
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: ntypat
  integer,intent(in) :: prtvol
  integer,intent(in) :: spaceComm_in
  type(pawfgr_type),intent(in) :: pawfgr
  real(dp),intent(in) :: ucvol
  integer,optional,target,intent(in) :: mpi_atmtab(:)
  integer,intent(in) :: ngfft(18)
  real(dp),optional,intent(out) :: abs_n_tilde_nt_diff(nspden)
  integer,intent(in) :: atindx1(natom)
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: nattyp(ntypat)
  real(dp),intent(in) :: nhat(pawfgr%nfft,nspden)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawrhoij_type),intent(in) :: pawrhoij(my_natom)
  type(pawtab_type),target,intent(in) :: pawtab(ntypat)
  real(dp),intent(in) :: rhor(pawfgr%nfft,nspden)
  real(dp),intent(out) :: rhor_n_one(pawfgr%nfft,nspden)
  real(dp),intent(out) :: rhor_nt_one(pawfgr%nfft,nspden)
  real(dp),intent(out) :: rhor_paw(pawfgr%nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(inout) :: xred(3,natom)
  real(dp),optional,intent(in) :: znucl(ntypat)
 end subroutine denfgr
end interface

interface
 subroutine dsdr_k_paw(cprj_k,cprj_kb,dsdr,dtefield,kdir,kfor,mband,natom,ncpgr,typat)
  use defs_basis
  use m_efield
  use m_pawcprj
  implicit none
  integer,intent(in) :: kdir
  integer,intent(in) :: kfor
  integer,intent(in) :: mband
  integer,intent(in) :: natom
  integer,intent(in) :: ncpgr
  type(efield_type),intent(in) :: dtefield
  type(pawcprj_type),intent(in) :: cprj_k(natom,dtefield%nspinor*mband)
  type(pawcprj_type),intent(in) :: cprj_kb(natom,dtefield%nspinor*mband)
  real(dp),intent(inout) :: dsdr(2,natom,ncpgr,dtefield%nband_occ,dtefield%nband_occ)
  integer,intent(in) :: typat(natom)
 end subroutine dsdr_k_paw
end interface

interface
 subroutine expibi(calc_expibi,dkvecs,natom,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  real(dp),intent(out) :: calc_expibi(2,natom,3)
  real(dp),intent(in) :: dkvecs(3,3)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine expibi
end interface

interface
 subroutine fourier_interpol(cplex,nspden,optin,optout,nfft_in,ngfft_in,nfft_out,ngfft_out,&  
  &  paral_kgb,MPI_enreg,rhor_in,rhor_out,rhog_in,rhog_out)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: nfft_in
  integer,intent(in) :: nfft_out
  integer,intent(in) :: nspden
  integer,intent(in) :: optin
  integer,intent(in) :: optout
  integer,intent(in) :: paral_kgb
  type(mpi_type),intent(inout) :: MPI_enreg
  integer,intent(in) :: ngfft_in(18)
  integer,intent(in) :: ngfft_out(18)
  real(dp),intent(inout) :: rhog_in(2,nfft_in)
  real(dp),intent(out) :: rhog_out(2,nfft_out)
  real(dp),intent(inout) :: rhor_in(cplex*nfft_in,nspden)
  real(dp),intent(out) :: rhor_out(cplex*nfft_out,nspden)
 end subroutine fourier_interpol
end interface

interface
 subroutine initrhoij(cplex,lexexch,lpawu,my_natom,natom,&  
  &  nspden,nspinor,nsppol,ntypat,pawrhoij,pawspnorb,pawtab,spinat,typat,&  
  &  ngrhoij,nlmnmix,use_rhoij_,use_rhoijres,&  ! optional arguments
  &  mpi_atmtab,comm_atom) ! optional arguments (parallelism)
  use m_pawtab
  use m_pawrhoij
  use defs_basis
  implicit none
  integer,intent(in),optional :: comm_atom
  integer,intent(in) :: cplex
  integer,intent(in) :: my_natom
  integer,intent(in) :: natom
  integer,intent(in),optional :: ngrhoij
  integer,intent(in),optional :: nlmnmix
  integer,intent(in) :: nspden
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: ntypat
  integer,intent(in) :: pawspnorb
  integer,intent(in),optional :: use_rhoij_
  integer,intent(in),optional :: use_rhoijres
  integer,optional,target,intent(in) :: mpi_atmtab(:)
  integer,intent(in) :: lexexch(ntypat)
  integer,intent(in) :: lpawu(ntypat)
  type(pawrhoij_type),intent(inout) :: pawrhoij(my_natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(in) :: spinat(3,natom)
  integer,intent(in) :: typat(natom)
 end subroutine initrhoij
end interface

interface
 subroutine int_ang(ang_phipphj,mpsang)
  use defs_basis
  implicit none
  integer,intent(in) :: mpsang
  real(dp),intent(out) :: ang_phipphj(mpsang**2,mpsang**2,8)
 end subroutine int_ang
end interface

interface
 subroutine linear_optics_paw(filnam,filnam_out,mpi_enreg_seq)
  use defs_basis
  use defs_abitypes
  implicit none
  character(len=fnlen),intent(in) :: filnam
  character(len=fnlen),intent(in) :: filnam_out
  type(mpi_type),intent(inout) :: mpi_enreg_seq
 end subroutine linear_optics_paw
end interface

interface
 subroutine make_efg_onsite(efg,my_natom,natom,nsym,ntypat,paw_an,pawang,pawrhoij,pawrad,pawtab,&  
  &  rprimd,symrel,tnons,xred,&  
  &  mpi_atmtab,comm_atom) ! optional arguments (parallelism)
  use m_pawrad
  use m_pawang
  use m_pawrhoij
  use m_paw_an
  use m_pawtab
  use defs_basis
  implicit none
  integer,optional,intent(in) :: comm_atom
  integer,intent(in) :: my_natom
  integer,intent(in) :: natom
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  type(pawang_type),intent(in) :: pawang
  integer,optional,target,intent(in) :: mpi_atmtab(:)
  real(dp),intent(out) :: efg(3,3,natom)
  type(paw_an_type),intent(in) :: paw_an(my_natom)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawrhoij_type),intent(in) :: pawrhoij(my_natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine make_efg_onsite
end interface

interface
 subroutine make_fc_paw(fc,my_natom,natom,nspden,ntypat,pawrhoij,pawrad,pawtab,&  
  &  mpi_atmtab,comm_atom) ! optional arguments (parallelism)
  use defs_basis
  use m_pawrad
  use m_pawrhoij
  use m_pawtab
  implicit none
  integer,optional,intent(in) :: comm_atom
  integer,intent(in) :: my_natom
  integer,intent(in) :: natom
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  integer,optional,target,intent(in) :: mpi_atmtab(:)
  real(dp),intent(out) :: fc(nspden,natom)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawrhoij_type),intent(in) :: pawrhoij(my_natom)
  type(pawtab_type),target,intent(in) :: pawtab(ntypat)
 end subroutine make_fc_paw
end interface

interface
 subroutine nhatgrid(atindx1,gmet,my_natom,natom,nattyp,ngfft,ntypat,&  
  &  optcut,optgr0,optgr1,optgr2,optrad,pawfgrtab,pawtab,rprimd,typat,ucvol,xred,&  
  &  mpi_atmtab,comm_atom,comm_fft,distribfft,typord) ! optional arguments (parallelism)
  use m_distribfft
  use defs_basis
  use m_pawfgrtab
  use m_pawtab
  implicit none
  integer,optional,intent(in) :: comm_atom
  integer,optional,intent(in) :: comm_fft
  integer,intent(in) :: my_natom
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  integer,intent(in) :: optcut
  integer,intent(in) :: optgr0
  integer,intent(in) :: optgr1
  integer,intent(in) :: optgr2
  integer,intent(in) :: optrad
  integer,optional,intent(in) :: typord
  type(distribfft_type),optional,target,intent(in) :: distribfft
  real(dp),intent(in) :: ucvol
  integer,optional,target,intent(in) :: mpi_atmtab(:)
  integer,intent(in) :: ngfft(18)
  integer,intent(in),target :: atindx1(natom)
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in),target :: nattyp(ntypat)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(my_natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine nhatgrid
end interface

interface
 subroutine optics_paw(atindx1,cg,cprj,dimcprj,dtfil,dtset,eigen0,gprimd,hdr,kg,&  
  &  mband,mcg,mcprj,mkmem,mpi_enreg,mpsang,mpw,natom,nkpt,npwarr,nsppol,&  
  &  pawrad,pawtab)
  use defs_basis
  use defs_abitypes
  use m_pawcprj
  use m_pawrad
  use m_pawtab
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: mcg
  integer,intent(in) :: mcprj
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpw
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  type(hdr_type),intent(inout) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: atindx1(natom)
  real(dp),intent(inout) :: cg(2,mcg)
  type(pawcprj_type),target,intent(inout) :: cprj(natom,mcprj)
  integer,intent(in) :: dimcprj(natom)
  real(dp),intent(in) :: eigen0(mband*nkpt*nsppol)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: kg(3,mpw*mkmem)
  integer,intent(in) :: npwarr(nkpt)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(pawtab_type),target,intent(in) :: pawtab(dtset%ntypat)
 end subroutine optics_paw
end interface

interface
 subroutine optics_paw_core(atindx1,cprj,dimcprj,dtfil,dtset,eigen0,filpsp,hdr,&  
  &  mband,mcprj,mkmem,mpi_enreg,mpsang,natom,nkpt,nsppol,pawrad,pawtab)
  use defs_basis
  use defs_abitypes
  use m_pawcprj
  use m_pawrad
  use m_pawtab
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: mcprj
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpsang
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  type(hdr_type),intent(inout) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: atindx1(natom)
  type(pawcprj_type),target,intent(inout) :: cprj(natom,mcprj)
  integer,intent(in) :: dimcprj(natom)
  real(dp),intent(in) :: eigen0(mband*nkpt*nsppol)
  character(len=fnlen),intent(in) :: filpsp(dtset%ntypat)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(pawtab_type),target,intent(in) :: pawtab(dtset%ntypat)
 end subroutine optics_paw_core
end interface

interface
 subroutine partial_dos_fractions_paw(cprj,dimcprj,dos_fractions,dos_fractions_m,&  
  &  dos_fractions_paw1,dos_fractions_pawt1,&  
  &  dtset,fatbands_flag,mbesslang,mcprj,mkmem,&  
  &  mpi_enreg,prtdosm,ndosfraction,paw_dos_flag,pawrad,pawtab)
  use defs_basis
  use defs_abitypes
  use m_pawcprj
  use m_pawrad
  use m_pawtab
  implicit none
  integer,intent(in) :: fatbands_flag
  integer,intent(in) :: mbesslang
  integer,intent(in) :: mcprj
  integer,intent(in) :: mkmem
  integer,intent(in) :: ndosfraction
  integer,intent(in) :: paw_dos_flag
  integer,intent(in) :: prtdosm
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawcprj_type) :: cprj(dtset%natom,mcprj)
  integer,intent(in) :: dimcprj(dtset%natom)
  real(dp),intent(inout) :: dos_fractions(dtset%nkpt,dtset%mband,dtset%nsppol,ndosfraction)
  real(dp),intent(inout) :: dos_fractions_m(dtset%nkpt,dtset%mband, &
  &         dtset%nsppol,ndosfraction*mbesslang*min(max(prtdosm,fatbands_flag),1))
  real(dp),intent(out) :: dos_fractions_paw1(dtset%nkpt,dtset%mband, &
  &         dtset%nsppol,ndosfraction*paw_dos_flag)
  real(dp),intent(out) :: dos_fractions_pawt1(dtset%nkpt,dtset%mband, &
  &         dtset%nsppol,ndosfraction*paw_dos_flag)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(pawtab_type),target,intent(in) :: pawtab(dtset%ntypat)
 end subroutine partial_dos_fractions_paw
end interface

interface
 subroutine paw_mknewh0(my_natom,nsppol,nspden,nfftf,pawspnorb,pawprtvol,Cryst,&  
  &  Pawtab,Paw_an,Paw_ij,Pawang,Pawfgrtab,vxc,vxc_val,vtrial,&  
  &  mpi_atmtab,comm_atom) ! optional arguments (parallelism)
  use defs_basis
  use m_paw_ij
  use m_pawang
  use m_pawfgrtab
  use m_paw_an
  use m_crystal
  use m_pawtab
  implicit none
  integer,optional,intent(in) :: comm_atom
  integer,intent(in) :: my_natom
  integer,intent(in) :: nfftf
  integer,intent(in) :: nspden
  integer,intent(in) :: nsppol
  integer,intent(in) :: pawprtvol
  integer,intent(in) :: pawspnorb
  type(crystal_t),intent(in) :: Cryst
  type(pawang_type),intent(in) :: Pawang
  integer,optional,target,intent(in) :: mpi_atmtab(:)
  type(paw_an_type),intent(in) :: Paw_an(my_natom)
  type(paw_ij_type),intent(inout) :: Paw_ij(my_natom)
  type(pawfgrtab_type),intent(inout) :: Pawfgrtab(my_natom)
  type(pawtab_type),target,intent(in) :: Pawtab(Cryst%ntypat)
  real(dp),intent(in) :: vtrial(nfftf,nspden)
  real(dp),intent(in) :: vxc(nfftf,nspden)
  real(dp),intent(in) :: vxc_val(nfftf,nspden)
 end subroutine paw_mknewh0
end interface

interface
 subroutine paw_symcprj(ik_bz,nspinor,nband_k,Cryst,Kmesh,Pawtab,Pawang,Cprj_bz)
  use m_bz_mesh
  use m_pawcprj
  use m_pawang
  use m_pawtab
  use m_crystal
  implicit none
  integer,intent(in) :: ik_bz
  integer,intent(in) :: nband_k
  integer,intent(in) :: nspinor
  type(crystal_t),intent(in) :: Cryst
  type(kmesh_t),intent(in) :: Kmesh
  type(pawang_type),intent(in) :: Pawang
  type(pawcprj_type),intent(inout) :: Cprj_bz(Cryst%natom,nspinor*nband_k)
  type(pawtab_type),target,intent(in) :: Pawtab(Cryst%ntypat)
 end subroutine paw_symcprj
end interface

interface
 subroutine paw_symcprj_op(ik_bz,nspinor,nband_k,Cryst,Kmesh,Pawtab,Pawang,in_Cprj,out_Cprj)
  use m_bz_mesh
  use m_pawcprj
  use m_pawang
  use m_pawtab
  use m_crystal
  implicit none
  integer,intent(in) :: ik_bz
  integer,intent(in) :: nband_k
  integer,intent(in) :: nspinor
  type(crystal_t),intent(in) :: Cryst
  type(kmesh_t),intent(in) :: Kmesh
  type(pawang_type),intent(in) :: Pawang
  type(pawtab_type),target,intent(in) :: Pawtab(Cryst%ntypat)
  type(pawcprj_type),intent(in) :: in_Cprj(Cryst%natom,nspinor*nband_k)
  type(pawcprj_type),intent(inout) :: out_Cprj(Cryst%natom,nspinor*nband_k)
 end subroutine paw_symcprj_op
end interface

interface
 subroutine pawaccrhoij(atindx,cplex,cwaveprj,cwaveprj1,ipert,isppol,my_natom,natom,&  
  &  nspinor,occ_k,option,pawrhoij,usetimerev,wtk_k,occ_k_2,&  
  &  comm_atom,mpi_atmtab ) ! optional (parallelism)
  use defs_basis
  use m_pawcprj
  use m_pawrhoij
  implicit none
  integer,optional,intent(in) :: comm_atom
  integer,intent(in) :: cplex
  integer,intent(in) :: ipert
  integer,intent(in) :: isppol
  integer,intent(in) :: my_natom
  integer,intent(in) :: natom
  integer,intent(in) :: nspinor
  integer,intent(in) :: option
  real(dp),intent(in) :: occ_k
  real(dp),optional,intent(in) :: occ_k_2
  logical,intent(in) :: usetimerev
  real(dp),intent(in) :: wtk_k
  integer,optional,target,intent(in) :: mpi_atmtab(:)
  integer,intent(in) :: atindx(natom)
  type(pawcprj_type),intent(in) :: cwaveprj(natom,nspinor)
  type(pawcprj_type),intent(in) :: cwaveprj1(natom,nspinor)
  type(pawrhoij_type),intent(inout) :: pawrhoij(my_natom)
 end subroutine pawaccrhoij
end interface

interface
 subroutine pawdenpot(compch_sph,epaw,epawdc,ipert,ixc,&  
  &  my_natom,natom,nspden,ntypat,nucdipmom,nzlmopt,option,paw_an,paw_an0,&  
  &  paw_ij,pawang,pawprtvol,pawrad,pawrhoij,pawspnorb,pawtab,pawxcdev,spnorbscl,xclevel,xc_denpos,ucvol,znucl,&  
  &  electronpositron,mpi_atmtab,comm_atom,vpotzero,hyb_mixing,hyb_mixing_sr) ! optional arguments
  use m_pawrad
  use m_paw_ij
  use m_pawang
  use m_pawrhoij
  use m_paw_an
  use m_pawtab
  use defs_basis
  use m_electronpositron
  implicit none
  integer,optional,intent(in) :: comm_atom
  integer,intent(in) :: ipert
  integer,intent(in) :: ixc
  integer,intent(in) :: my_natom
  integer,intent(in) :: natom
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  integer,intent(in) :: nzlmopt
  integer,intent(in) :: option
  integer,intent(in) :: pawprtvol
  integer,intent(in) :: pawspnorb
  integer,intent(in) :: pawxcdev
  integer,intent(in) :: xclevel
  real(dp),intent(out) :: compch_sph
  type(electronpositron_type),pointer,optional :: electronpositron
  real(dp),intent(out) :: epaw
  real(dp),intent(out) :: epawdc
  real(dp),intent(in),optional :: hyb_mixing
  real(dp),intent(in),optional :: hyb_mixing_sr
  type(pawang_type),intent(in) :: pawang
  real(dp), intent(in) :: spnorbscl
  real(dp), intent(in) :: ucvol
  real(dp), intent(in) :: xc_denpos
  integer,optional,target,intent(in) :: mpi_atmtab(:)
  real(dp),intent(in) :: nucdipmom(3,my_natom)
  type(paw_an_type),intent(inout) :: paw_an(my_natom)
  type(paw_an_type), intent(in) :: paw_an0(my_natom)
  type(paw_ij_type),intent(inout) :: paw_ij(my_natom)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawrhoij_type),intent(inout) :: pawrhoij(my_natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(out),optional :: vpotzero(2)
  real(dp),intent(in) :: znucl(ntypat)
 end subroutine pawdenpot
end interface

interface
 subroutine pawdensities(compch_sph,cplex,iatom,lmselectin,lmselectout,lm_size,nhat1,nspden,nzlmopt,&  
  &  opt_compch,opt_dens,opt_l,opt_print,pawang,pawprtvol,pawrad,pawrhoij,pawtab,rho1,trho1,&  
  &  one_over_rad2) ! optional
  use defs_basis
  use m_pawang
  use m_pawrad
  use m_pawrhoij
  use m_pawtab
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: iatom
  integer,intent(in) :: lm_size
  integer,intent(in) :: nspden
  integer,intent(in) :: nzlmopt
  integer,intent(in) :: opt_compch
  integer,intent(in) :: opt_dens
  integer,intent(in) :: opt_l
  integer,intent(in) :: opt_print
  integer,intent(in) :: pawprtvol
  real(dp),intent(inout) :: compch_sph
  type(pawang_type),intent(in) :: pawang
  type(pawrad_type),intent(in) :: pawrad
  type(pawrhoij_type),intent(in) :: pawrhoij
  type(pawtab_type),intent(in) :: pawtab
  logical,intent(in) :: lmselectin(lm_size)
  logical,intent(inout) :: lmselectout(lm_size)
  real(dp),intent(out) :: nhat1(cplex*pawtab%mesh_size,lm_size,nspden*(1-((opt_dens+1)/2)))
  real(dp),intent(in),target,optional :: one_over_rad2(pawtab%mesh_size)
  real(dp),intent(out) :: rho1(cplex*pawtab%mesh_size,lm_size,nspden)
  real(dp),intent(out) :: trho1(cplex*pawtab%mesh_size,lm_size,nspden*(1-(opt_dens/2)))
 end subroutine pawdensities
end interface

interface
 subroutine pawdfptenergy(delta_energy,ipert1,ipert2,ixc,my_natom,natom,ntypat,nzlmopt_a,nzlmopt_b,&  
  &  paw_an0,paw_an1,paw_ij1,pawang,pawprtvol,pawrad,pawrhoij_a,pawrhoij_b,&  
  &  pawtab,pawxcdev,xclevel,&  
  &  mpi_atmtab,comm_atom) ! optional arguments (parallelism)
  use m_pawrad
  use m_paw_ij
  use m_pawang
  use m_pawrhoij
  use m_paw_an
  use m_pawtab
  use defs_basis
  implicit none
  integer,optional,intent(in) :: comm_atom
  integer,intent(in) :: ipert1
  integer,intent(in) :: ipert2
  integer,intent(in) :: ixc
  integer,intent(in) :: my_natom
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  integer,intent(in) :: nzlmopt_a
  integer,intent(in) :: nzlmopt_b
  integer,intent(in) :: pawprtvol
  integer,intent(in) :: pawxcdev
  integer,intent(in) :: xclevel
  type(pawang_type),intent(in) :: pawang
  integer,optional,target,intent(in) :: mpi_atmtab(:)
  real(dp),intent(out) :: delta_energy(2)
  type(paw_an_type),intent(in) :: paw_an0(my_natom)
  type(paw_an_type),intent(inout) :: paw_an1(my_natom)
  type(paw_ij_type),intent(inout) :: paw_ij1(my_natom)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawrhoij_type),intent(in) :: pawrhoij_a(my_natom)
  type(pawrhoij_type),intent(in) :: pawrhoij_b(my_natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
 end subroutine pawdfptenergy
end interface

interface
 subroutine pawgrnl(atindx1,dimnhat,dyfrnl,dyfr_cplex,eltfrnl,grnl,gsqcut,mgfft,my_natom,natom,&  
  &  nattyp,nfft,ngfft,nhat,nlstr,nspden,nsym,ntypat,optgr,optgr2,optstr,optstr2,&  
  &  pawang,pawfgrtab,pawrhoij,pawtab,ph1d,psps,qphon,rprimd,symrec,typat,vtrial,vxc,xred,&  
  &  mpi_atmtab,comm_atom,comm_fft,mpi_comm_grid,me_g0,paral_kgb,distribfft) ! optional arguments (parallelism)
  use m_distribfft
  use m_pawang
  use m_pawrhoij
  use m_pawfgrtab
  use m_pawtab
  use defs_basis
  use defs_datatypes
  implicit none
  integer,optional,intent(in) :: comm_atom
  integer,optional,intent(in) :: comm_fft
  integer,intent(in) :: dimnhat
  integer,intent(in) :: dyfr_cplex
  integer,optional,intent(in) :: me_g0
  integer,intent(in) :: mgfft
  integer,optional,intent(in) :: mpi_comm_grid
  integer,intent(in) :: my_natom
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(in) :: optgr
  integer,intent(in) :: optgr2
  integer,intent(in) :: optstr
  integer,intent(in) :: optstr2
  integer,optional,intent(in) :: paral_kgb
  type(distribfft_type),optional,target,intent(in) :: distribfft
  real(dp),intent(in) :: gsqcut
  type(pawang_type),intent(in) :: pawang
  type(pseudopotential_type),intent(in) :: psps
  integer,optional,target,intent(in) :: mpi_atmtab(:)
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: atindx1(natom)
  real(dp),intent(inout) :: dyfrnl(dyfr_cplex,3,3,natom,natom*optgr2)
  real(dp),intent(inout) :: eltfrnl(6+3*natom,6)
  real(dp),intent(inout) :: grnl(3*natom*optgr)
  integer,intent(in) :: nattyp(ntypat)
  real(dp),intent(in) :: nhat(nfft,dimnhat)
  real(dp),intent(inout) :: nlstr(6*optstr)
  type(pawfgrtab_type),target,intent(inout) :: pawfgrtab(:)
  type(pawrhoij_type),target,intent(inout) :: pawrhoij(:)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symrec(3,3,nsym)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in),target :: vtrial(nfft,nspden)
  real(dp),intent(in) :: vxc(nfft,nspden)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine pawgrnl
end interface

interface
 subroutine pawgylmg(gprimd,gylmg,kg,kpg,kpt,lmax,nkpg,npw,ntypat,pawtab,ylm)
  use defs_basis
  use m_pawtab
  implicit none
  integer,intent(in) :: lmax
  integer,intent(in) :: nkpg
  integer,intent(in) :: npw
  integer,intent(in) :: ntypat
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(out) :: gylmg(npw,lmax**2,ntypat)
  integer,intent(in) :: kg(3,npw)
  real(dp),intent(in) :: kpg(npw,nkpg)
  real(dp),intent(in) :: kpt(3)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(in) :: ylm(npw,lmax**2)
 end subroutine pawgylmg
end interface

interface
 subroutine pawinit(gnt_option,gsqcut_eff,hyb_range,lcutdens,lmix,mpsang,nphi,nsym,ntheta,&  
  &  pawang,pawrad,pawspnorb,pawtab,pawxcdev,xclevel,usepotzero)
  use defs_basis
  use m_pawang
  use m_pawrad
  use m_pawtab
  implicit none
  integer,intent(in) :: gnt_option
  integer,intent(in) :: lcutdens
  integer,intent(in) :: lmix
  integer,intent(in) :: mpsang
  integer,intent(in) :: nphi
  integer,intent(in) :: nsym
  integer,intent(in) :: ntheta
  integer,intent(in) :: pawspnorb
  integer,intent(in) :: pawxcdev
  integer,intent(in) :: usepotzero
  integer,intent(in) :: xclevel
  real(dp),intent(in) :: gsqcut_eff
  real(dp),intent(in) :: hyb_range
  type(pawang_type),intent(inout) :: pawang
  type(pawrad_type),intent(in) :: pawrad(:)
  type(pawtab_type),target,intent(inout) :: pawtab(:)
 end subroutine pawinit
end interface

interface
 subroutine paw_gencond(Dtset,gnt_option,mode,call_pawinit) 
  use defs_abitypes
  implicit none
  integer,intent(in) :: gnt_option
  type(dataset_type),intent(in) :: Dtset
  logical,intent(out) :: call_pawinit
  character(len=*),intent(in) :: mode
 end subroutine paw_gencond
end interface

interface
 subroutine pawmkaewf(Dtset,crystal,ebands,my_natom,mpw,mband,mcg,mcprj,nkpt,mkmem,nsppol,nband,&  
  &  istwfk,npwarr,kpt,ngfftf,kg,dimcprj,Pawfgrtab,Pawrad,Pawtab,&  
  &  Hdr,Dtfil,cg,Cprj,MPI_enreg,ierr,pseudo_norms,set_k,set_band ,&  
  &  mpi_atmtab,comm_atom) ! Optional arguments
  use m_pawtab
  use m_pawrad
  use defs_abitypes
  use m_pawcprj
  use m_pawfgrtab
  use m_crystal
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in),optional :: comm_atom
  integer,intent(out) :: ierr
  integer,intent(in) :: mband
  integer,intent(in) :: mcg
  integer,intent(in) :: mcprj
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: my_natom
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  integer,intent(in),optional :: set_band
  integer,intent(in),optional :: set_k
  type(datafiles_type),intent(in) :: Dtfil
  type(dataset_type),intent(in) :: Dtset
  type(hdr_type),intent(inout) :: Hdr
  type(mpi_type),intent(inout) :: MPI_enreg
  type(crystal_t),intent(in) :: crystal
  type(ebands_t),intent(in) :: ebands
  integer,optional,target,intent(in) :: mpi_atmtab(:)
  integer,intent(in) :: ngfftf(18)
  type(pawcprj_type),intent(in) :: Cprj(crystal%natom,mcprj)
  type(pawfgrtab_type),intent(in) :: Pawfgrtab(my_natom)
  type(pawrad_type),intent(in) :: Pawrad(crystal%ntypat)
  type(pawtab_type),intent(in) :: Pawtab(crystal%ntypat)
  real(dp),intent(in) :: cg(2,mcg)
  integer,intent(in) :: dimcprj(crystal%natom)
  integer,intent(in) :: istwfk(nkpt)
  integer,intent(in) :: kg(3,mpw*mkmem)
  real(dp),intent(in) :: kpt(3,nkpt)
  integer,intent(in) :: nband(nkpt*nsppol)
  integer,intent(in) :: npwarr(nkpt)
  real(dp),optional,intent(out) :: pseudo_norms(nsppol,nkpt,mband)
 end subroutine pawmkaewf
end interface

interface
 subroutine pawmknhat(compch_fft,cplex,ider,idir,ipert,izero,gprimd,&  
  &  my_natom,natom,nfft,ngfft,nhatgrdim,nspden,ntypat,pawang,pawfgrtab,&  
  &  pawgrnhat,pawnhat,pawrhoij,pawrhoij0,pawtab,qphon,rprimd,ucvol,usewvl,xred,&  
  &  mpi_atmtab,comm_atom,comm_fft,mpi_comm_wvl,me_g0,paral_kgb,distribfft) ! optional arguments
  use m_distribfft
  use m_pawang
  use m_pawrhoij
  use m_pawfgrtab
  use m_pawtab
  use defs_basis
  implicit none
  integer,optional,intent(in) :: comm_atom
  integer,optional,intent(in) :: comm_fft
  integer,intent(in) :: cplex
  integer,intent(in) :: ider
  integer,intent(in) :: idir
  integer,intent(in) :: ipert
  integer,intent(in) :: izero
  integer,optional,intent(in) :: me_g0
  integer,optional,intent(in) :: mpi_comm_wvl
  integer,intent(in) :: my_natom
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nhatgrdim
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  integer,optional,intent(in) :: paral_kgb
  integer,intent(in) :: usewvl
  real(dp),intent(out) :: compch_fft
  type(distribfft_type),optional,intent(in),target :: distribfft
  type(pawang_type),intent(in) :: pawang
  real(dp),intent(in) :: ucvol
  integer,optional,target,intent(in) :: mpi_atmtab(:)
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: gprimd(3,3)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(my_natom)
  real(dp),intent(out) :: pawgrnhat(cplex*nfft,nspden,3*nhatgrdim)
  real(dp),intent(inout) :: pawnhat(cplex*nfft,nspden)
  type(pawrhoij_type),intent(in) :: pawrhoij(my_natom)
  type(pawrhoij_type),intent(in) :: pawrhoij0(my_natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine pawmknhat
end interface

interface
 subroutine pawmknhat_psipsi(cprj1,cprj2,ider,izero,my_natom,natom,nfft,ngfft,nhat12_grdim,&  
  &  nspinor,ntypat,pawang,pawfgrtab,grnhat12,nhat12,pawtab,&  
  &  gprimd,qphon,xred,mpi_atmtab,comm_atom,comm_fft,me_g0,paral_kgb,distribfft) ! optional arguments
  use m_distribfft
  use m_pawang
  use m_pawcprj
  use m_pawfgrtab
  use defs_basis
  use m_pawtab
  implicit none
  integer,optional,intent(in) :: comm_atom
  integer,optional,intent(in) :: comm_fft
  integer,intent(in) :: ider
  integer,intent(in) :: izero
  integer,optional,intent(in) :: me_g0
  integer,intent(in) :: my_natom
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nhat12_grdim
  integer,intent(in) :: nspinor
  integer,intent(in) :: ntypat
  integer,optional,intent(in) :: paral_kgb
  type(distribfft_type),optional,intent(in),target :: distribfft
  type(pawang_type),intent(in) :: pawang
  integer,optional,target,intent(in) :: mpi_atmtab(:)
  integer,intent(in) :: ngfft(18)
  type(pawcprj_type),intent(in) :: cprj1(natom,nspinor)
  type(pawcprj_type),intent(in) :: cprj2(natom,nspinor)
  real(dp),optional, intent(in) :: gprimd(3,3)
  real(dp),intent(out) :: grnhat12(2,nfft,nspinor**2,3*nhat12_grdim)
  real(dp),intent(out) :: nhat12(2,nfft,nspinor**2)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(my_natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),optional, intent(in) :: qphon(3)
  real(dp),optional, intent(in) :: xred(3,natom)
 end subroutine pawmknhat_psipsi
end interface

interface
 subroutine pawmkrho(compch_fft,cplex,gprimd,idir,indsym,ipert,mpi_enreg,&  
  &  my_natom,natom,nspden,nsym,ntypat,paral_kgb,pawang,pawfgr,pawfgrtab,pawprtvol,&  
  &  pawrhoij,pawrhoij_unsym,&  
  &  pawtab,qphon,rhopsg,rhopsr,rhor,rprimd,symafm,symrec,typat,ucvol,usewvl,xred,&  
  &  pawang_sym,pawnhat,pawrhoij0,rhog) ! optional arguments
  use m_pawtab
  use m_pawang
  use m_pawrhoij
  use defs_abitypes
  use m_pawfgrtab
  use m_pawfgr
  use defs_basis
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: idir
  integer,intent(in) :: ipert
  integer,intent(in) :: my_natom
  integer,intent(in) :: natom
  integer,intent(in) :: nspden
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: pawprtvol
  integer,intent(in) :: usewvl
  real(dp),intent(out) :: compch_fft
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pawang_type),intent(in),optional :: pawang_sym
  type(pawfgr_type),intent(in) :: pawfgr
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: indsym(4,nsym,natom)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(my_natom)
  real(dp),intent(inout),target,optional :: pawnhat(cplex*pawfgr%nfft,nspden)
  type(pawrhoij_type),intent(inout),target :: pawrhoij(:)
  type(pawrhoij_type),intent(in),target,optional :: pawrhoij0(my_natom)
  type(pawrhoij_type),intent(inout) :: pawrhoij_unsym(:)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(out),optional :: rhog(2,pawfgr%nfft)
  real(dp),intent(inout) :: rhopsg(2,pawfgr%nfftc)
  real(dp),intent(inout) :: rhopsr(cplex*pawfgr%nfftc,nspden)
  real(dp),intent(inout) :: rhor(cplex*pawfgr%nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symafm(nsym)
  integer,intent(in) :: symrec(3,3,nsym)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine pawmkrho
end interface

interface
 subroutine pawmkrhoij(atindx,atindx1,cprj,dimcprj,istwfk,kptopt,mband,mband_cprj,mcprj,mkmem,mpi_enreg,&  
  &  natom,nband,nkpt,nspinor,nsppol,occ,paral_kgb,paw_dmft,&  
  &  pawprtvol,pawrhoij,unpaw,usewvl,wtk)
  use defs_basis
  use defs_abitypes
  use m_pawcprj
  use m_pawrhoij
  use m_paw_dmft
  implicit none
  integer,intent(in) :: kptopt
  integer,intent(in) :: mband
  integer,intent(in) :: mband_cprj
  integer,intent(in) :: mcprj
  integer,intent(in) :: mkmem
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: pawprtvol
  integer,intent(in) :: unpaw
  integer,intent(in) :: usewvl
  type(mpi_type),intent(inout) :: mpi_enreg
  type(paw_dmft_type),intent(in) :: paw_dmft
  integer,intent(in) :: atindx(natom)
  integer,intent(in) :: atindx1(natom)
  type(pawcprj_type),target,intent(in) :: cprj(natom,mcprj)
  integer,intent(in) :: dimcprj(natom)
  integer,intent(in) :: istwfk(nkpt)
  integer,intent(in) :: nband(nkpt*nsppol)
  real(dp),intent(in) :: occ(mband*nkpt*nsppol)
  type(pawrhoij_type),intent(inout),target :: pawrhoij(:)
  real(dp),intent(in) :: wtk(nkpt)
 end subroutine pawmkrhoij
end interface

interface
 subroutine pawnabla_init(mpsang,ntypat,pawrad,pawtab)
  use m_pawtab
  use m_pawrad
  implicit none
  integer,intent(in) :: mpsang
  integer,intent(in) :: ntypat
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawtab_type),target,intent(inout) :: pawtab(ntypat)
 end subroutine pawnabla_init
end interface

interface
 subroutine pawnhatfr(ider,idir,ipert,my_natom,natom,nspden,ntypat,&  
  &  pawang,pawfgrtab,pawrhoij,pawtab,rprimd,&  
  &  mpi_atmtab,comm_atom) ! optional arguments (parallelism)
  use m_pawtab
  use m_pawang
  use m_pawfgrtab
  use m_pawrhoij
  use defs_basis
  implicit none
  integer,optional,intent(in) :: comm_atom
  integer,intent(in) :: ider
  integer,intent(in) :: idir
  integer,intent(in) :: ipert
  integer,intent(in) :: my_natom
  integer,intent(in) :: natom
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  type(pawang_type),intent(in) :: pawang
  integer,optional,target,intent(in) :: mpi_atmtab(:)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(my_natom)
  type(pawrhoij_type),intent(in) :: pawrhoij(my_natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine pawnhatfr
end interface

interface
 subroutine pawpolev(my_natom,natom,ntypat,pawrhoij,pawtab,pelev,&  
  &  comm_atom) ! optional argument (parallelism)
  use m_pawtab
  use m_pawrhoij
  use defs_basis
  implicit none
  integer,optional,intent(in) :: comm_atom
  integer,intent(in) :: my_natom
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  type(pawrhoij_type),intent(in) :: pawrhoij(my_natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(out) :: pelev(3)
 end subroutine pawpolev
end interface

interface
 subroutine pawprt(dtset,my_natom,paw_ij,pawrhoij,pawtab,&  
  &  electronpositron,&  ! optional argument
  &  mpi_atmtab,comm_atom) ! optional arguments (parallelism)
  use m_pawtab
  use m_paw_ij
  use defs_abitypes
  use m_electronpositron
  use m_pawrhoij
  implicit none
  integer,optional,intent(in) :: comm_atom
  integer,intent(in) :: my_natom
  type(dataset_type),intent(in) :: dtset
  type(electronpositron_type),pointer,optional :: electronpositron
  integer,optional,target,intent(in) :: mpi_atmtab(:)
  type(paw_ij_type),target,intent(inout) :: paw_ij(my_natom)
  type(pawrhoij_type),target,intent(inout) :: pawrhoij(my_natom)
  type(pawtab_type),target,intent(in) :: pawtab(dtset%ntypat)
 end subroutine pawprt
end interface

interface
 subroutine pawpuxinit(dmatpuopt,exchmix,f4of2_sla,f6of2_sla,jpawu,llexexch,llpawu,&  
  &  ntypat,pawang,pawprtvol,pawrad,pawtab,upawu,use_dmft,useexexch,usepawu,ucrpa)
  use defs_basis
  use m_pawang
  use m_pawrad
  use m_pawtab
  implicit none
  integer,intent(in) :: dmatpuopt
  integer,intent(in) :: ntypat
  integer,intent(in) :: pawprtvol
  integer,optional, intent(in) :: ucrpa
  integer,intent(in) :: use_dmft
  integer,intent(in) :: useexexch
  integer,intent(in) :: usepawu
  real(dp),intent(in) :: exchmix
  type(pawang_type), intent(in) :: pawang
  real(dp),intent(in) :: f4of2_sla(ntypat)
  real(dp),intent(in) :: f6of2_sla(ntypat)
  real(dp),intent(in) :: jpawu(ntypat)
  integer,intent(in) :: llexexch(ntypat)
  integer,intent(in) :: llpawu(ntypat)
  type(pawrad_type),intent(inout) :: pawrad(ntypat)
  type(pawtab_type),target,intent(inout) :: pawtab(ntypat)
  real(dp),intent(in) :: upawu(ntypat)
 end subroutine pawpuxinit
end interface

interface
 subroutine pawsushat(atindx,cprj_k,gbound_diel,gylmg_diel,iband1,iband2,ispinor1,ispinor2,istwf_k,kg_diel,&  
  &  lmax_diel,mgfftdiel,natom,nband,ndiel4,ndiel5,ndiel6,&  
  &  ngfftdiel,npwdiel,nspinor,ntypat,optreal,&  
  &  pawang,pawtab,ph3d_diel,typat,wfprod,wfraug,&  
  &  mpi_atmtab,comm_atom,comm_fft,me_g0,paral_kgb,distribfft) ! optional arguments (parallelism)
  use defs_basis
  use m_pawang
  use m_pawcprj
  use m_distribfft
  use m_pawtab
  implicit none
  integer,optional,intent(in) :: comm_atom
  integer,optional,intent(in) :: comm_fft
  integer,intent(in) :: iband1
  integer,intent(in) :: iband2
  integer,intent(in) :: ispinor1
  integer,intent(in) :: ispinor2
  integer,intent(in) :: istwf_k
  integer,intent(in) :: lmax_diel
  integer,optional,intent(in) :: me_g0
  integer,intent(in) :: mgfftdiel
  integer,intent(in) :: natom
  integer,intent(in) :: nband
  integer,intent(in) :: ndiel4
  integer,intent(in) :: ndiel5
  integer,intent(in) :: ndiel6
  integer,intent(in) :: npwdiel
  integer,intent(in) :: nspinor
  integer,intent(in) :: ntypat
  integer,intent(in) :: optreal
  integer,optional,intent(in) :: paral_kgb
  type(distribfft_type),optional,intent(in),target :: distribfft
  type(pawang_type),intent(in) :: pawang
  integer,optional,target,intent(in) :: mpi_atmtab(:)
  integer,intent(in) :: ngfftdiel(18)
  integer,intent(in) :: atindx(natom)
  type(pawcprj_type),intent(in) :: cprj_k(natom,nspinor*nband)
  integer,intent(in) :: gbound_diel(2*mgfftdiel+8,2)
  real(dp),intent(in) :: gylmg_diel(npwdiel,lmax_diel**2,ntypat)
  integer,intent(in) :: kg_diel(3,npwdiel)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(in) :: ph3d_diel(2,npwdiel,natom)
  integer,intent(in) :: typat(natom)
  real(dp),intent(inout) :: wfprod(2,npwdiel*(1-optreal))
  real(dp),intent(inout) :: wfraug(2,ndiel4,ndiel5,ndiel6*optreal)
 end subroutine pawsushat
end interface

interface
 subroutine pawuenergy(iatom,eldaumdc,eldaumdcdc,pawprtvol,pawtab,paw_ij,&  
  &  dmft_dc,e_ee,e_dc,e_dcdc,u_dmft,j_dmft) ! optional arguments (DMFT)
  use defs_basis
  use m_paw_ij
  use m_pawtab
  implicit none
  integer,optional,intent(in) :: dmft_dc
  integer,intent(in) :: iatom
  integer,intent(in) :: pawprtvol
  real(dp),optional,intent(inout) :: e_dc
  real(dp),optional,intent(inout) :: e_dcdc
  real(dp),optional,intent(inout) :: e_ee
  real(dp),intent(inout) :: eldaumdc
  real(dp),intent(inout) :: eldaumdcdc
  real(dp),optional,intent(in) :: j_dmft
  type(paw_ij_type),intent(in) :: paw_ij
  type(pawtab_type),intent(in) :: pawtab
  real(dp),optional,intent(in) :: u_dmft
 end subroutine pawuenergy
end interface

interface
 subroutine pawuj_det(dtpawuj,ndtpawuj,ujdet_filename,ures)
  use defs_basis
  use defs_abitypes
  implicit none
  integer :: ndtpawuj
  character(len=*),intent(in) :: ujdet_filename
  real(dp),intent(out) :: ures
  type(macro_uj_type),intent(in) :: dtpawuj(0:ndtpawuj)
 end subroutine pawuj_det
end interface

interface
 subroutine pawuj_ini(dtpawuj,ndtset)
  use defs_abitypes
  implicit none
  integer :: ndtset
  type(macro_uj_type),intent(inout) :: dtpawuj(0:ndtset)
 end subroutine pawuj_ini
end interface

interface
 subroutine pawuj_red(dtset,dtpawuj,fatvshift,my_natom,natom,ntypat,paw_ij,pawrad,pawtab,ndtpawuj,&  
  &  mpi_atmtab,comm_atom) ! optional arguments (parallelism)
  use defs_basis
  use m_paw_ij
  use defs_abitypes
  use m_pawrad
  use m_pawtab
  implicit none
  integer,optional,intent(in) :: comm_atom
  integer,intent(in) :: my_natom
  integer,intent(in) :: natom
  integer,intent(in) :: ndtpawuj
  integer,intent(in) :: ntypat
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: fatvshift
  integer,optional,target,intent(in) :: mpi_atmtab(:)
  type(macro_uj_type),intent(inout) :: dtpawuj(0:ndtpawuj)
  type(paw_ij_type),intent(in) :: paw_ij(my_natom)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
 end subroutine pawuj_red
end interface

interface
 subroutine linvmat(inmat,oumat,nat,nam,option,gam,prtvol)
  use defs_basis
  implicit none
  integer,intent(in) :: nat
  integer,intent(in),optional :: option
  integer,intent(in),optional :: prtvol
  real(dp),intent(in) :: gam
  character(len=500),intent(in) :: nam
  real(dp),intent(in) :: inmat(nat,nat)
  real(dp),intent(inout) :: oumat(:,:)
 end subroutine linvmat
end interface

interface
 subroutine lprtmat(commnt,chan,prtvol,mmat,nat)
  use defs_basis
  implicit none
  integer,intent(in) :: chan
  integer,intent(in) :: nat
  integer,intent(in) :: prtvol
  character(len=500),intent(in) :: commnt
  real(dp),intent(in) :: mmat(nat,nat)
 end subroutine lprtmat
end interface

interface
 subroutine lcalcu(magv,natom,rprimd,xred,chi,chi0,pawujat,ures,prtvol,gam,opt)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in),optional :: opt
  integer,intent(in),optional :: pawujat
  integer,intent(in),optional :: prtvol
  real(dp),intent(in),optional :: gam
  real(dp),intent(out) :: ures
  real(dp),intent(in) :: chi(natom)
  real(dp),intent(in) :: chi0(natom)
  integer,intent(in) :: magv(natom)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine lcalcu
end interface

interface
 subroutine pawuj_free(dtpawuj)
  use defs_abitypes
  implicit none
  type(macro_uj_type),intent(inout) :: dtpawuj
 end subroutine pawuj_free
end interface

interface
 subroutine pawxenergy(eexex,pawprtvol,pawrhoij,pawtab)
  use defs_basis
  use m_pawrhoij
  use m_pawtab
  implicit none
  integer,intent(in) :: pawprtvol
  real(dp),intent(inout) :: eexex
  type(pawrhoij_type),intent(in) :: pawrhoij
  type(pawtab_type),intent(in) :: pawtab
 end subroutine pawxenergy
end interface

interface
 subroutine prtfatbands(dos_fractions_m,dtset,fildata,fermie,eigen,&  
  &  mbesslang,m_dos_flag,ndosfraction,pawfatbnd,pawtab)
  use defs_basis
  use defs_abitypes
  use m_pawtab
  implicit none
  integer,intent(in) :: m_dos_flag
  integer,intent(in) :: mbesslang
  integer,intent(in) :: ndosfraction
  integer,intent(in) :: pawfatbnd
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: fermie
  character(len=fnlen),intent(in) :: fildata
  real(dp),intent(in) :: dos_fractions_m(dtset%nkpt,dtset%mband,dtset%nsppol,ndosfraction*mbesslang)
  real(dp),intent(in) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)
 end subroutine prtfatbands
end interface

interface
 subroutine qijb_kk(calc_qijb,dkvecs,expibi,gprimd,lmn2max,natom,ntypat,&  
  &  pawang,pawrad,pawtab,typat)
  use defs_basis
  use m_pawang
  use m_pawrad
  use m_pawtab
  implicit none
  integer,intent(in) :: lmn2max
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  type(pawang_type),intent(in) :: pawang
  real(dp),intent(out) :: calc_qijb(2,lmn2max,natom,3)
  real(dp),intent(in) :: dkvecs(3,3)
  real(dp),intent(in) :: expibi(2,natom,3)
  real(dp),intent(in) :: gprimd(3,3)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  integer,intent(in) :: typat(natom)
 end subroutine qijb_kk
end interface

interface
 subroutine read_atomden(MPI_enreg,natom,nspden,ntypat,pawfgr,&  
  &  rhor_paw,typat,rprimd,xred,prtvol,file_prefix)
  use m_pawfgr
  use defs_abitypes
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  integer,intent(in) :: prtvol
  type(mpi_type),intent(in) :: MPI_enreg
  character(len=7), intent(in) :: file_prefix
  type(pawfgr_type),intent(in) :: pawfgr
  real(dp),intent(inout) :: rhor_paw(pawfgr%nfft,nspden)
  real(dp), intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
  real(dp), intent(in) :: xred(3,natom)
 end subroutine read_atomden
end interface

interface
 subroutine setnoccmmp(compute_dmat,dimdmat,dmatpawu,dmatudiag,impose_dmat,indsym,my_natom,natom,&  
  &  natpawu,nspinor,nsppol,nsym,ntypat,paw_ij,pawang,pawprtvol,pawrhoij,pawtab,&  
  &  spinat,symafm,typat,useexexch,usepawu,&  
  &  mpi_atmtab,comm_atom) ! optional arguments (parallelism)
  use defs_basis
  use m_paw_ij
  use m_pawang
  use m_pawrhoij
  use m_pawtab
  implicit none
  integer,optional,intent(in) :: comm_atom
  integer,intent(in) :: compute_dmat
  integer,intent(in) :: dimdmat
  integer,intent(in) :: dmatudiag
  integer,intent(in) :: impose_dmat
  integer,intent(in) :: my_natom
  integer,intent(in) :: natom
  integer,intent(in) :: natpawu
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(in) :: pawprtvol
  integer,intent(in) :: useexexch
  integer,intent(in) :: usepawu
  type(pawang_type),intent(in) :: pawang
  integer,optional,target,intent(in) :: mpi_atmtab(:)
  real(dp),intent(in) :: dmatpawu(dimdmat,dimdmat,nspinor*nsppol,natpawu*impose_dmat)
  integer,intent(in) :: indsym(4,nsym,natom)
  type(paw_ij_type),intent(inout) :: paw_ij(my_natom)
  type(pawrhoij_type),intent(in) :: pawrhoij(my_natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(in) :: spinat(3,natom)
  integer,intent(in) :: symafm(nsym)
  integer,intent(in) :: typat(natom)
 end subroutine setnoccmmp
end interface

interface
 subroutine setrhoijpbe0(dtset,initialized,istep,istep_mix,&  
  &  mpi_comm_read,my_natom,natom,ntypat,pawrhoij,pawtab,typat,&  
  &  mpi_atmtab,comm_atom) ! optional arguments (parallelism)
  use m_pawtab
  use defs_abitypes
  use m_pawrhoij
  implicit none
  integer,optional,intent(in) :: comm_atom
  integer,intent(in) :: initialized
  integer,intent(in) :: istep
  integer,intent(inout) :: istep_mix
  integer,intent(in) :: mpi_comm_read
  integer,intent(in) :: my_natom
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  type(dataset_type),intent(in) :: dtset
  integer,optional,target,intent(in) :: mpi_atmtab(:)
  type(pawrhoij_type),intent(inout) :: pawrhoij(my_natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  integer,intent(in) :: typat(natom)
 end subroutine setrhoijpbe0
end interface

interface
 subroutine setsymrhoij(gprimd,lmax,nsym,pawprtvol,rprimd,sym,zarot)
  use defs_basis
  implicit none
  integer,intent(in) :: lmax
  integer,intent(in) :: nsym
  integer,intent(in) :: pawprtvol
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: sym(3,3,nsym)
  real(dp),intent(out) :: zarot(2*lmax+1,2*lmax+1,lmax+1,nsym)
 end subroutine setsymrhoij
end interface

interface
 subroutine simple_j_dia(jdia,natom,nfft,pawfgrtab)
  use defs_basis
  use m_pawfgrtab
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  real(dp),intent(out) :: jdia(3,3,nfft)
  type(pawfgrtab_type),intent(in) :: pawfgrtab(natom)
 end subroutine simple_j_dia
end interface

interface
 subroutine smatrix_k_paw(cprj_k,cprj_kb,dtefield,kdir,kfor,mband,natom,smat_k_paw,typat)
  use defs_basis
  use m_efield
  use m_pawcprj
  implicit none
  integer,intent(in) :: kdir
  integer,intent(in) :: kfor
  integer,intent(in) :: mband
  integer,intent(in) :: natom
  type(efield_type),intent(in) :: dtefield
  type(pawcprj_type),intent(in) :: cprj_k(natom,dtefield%nspinor*mband)
  type(pawcprj_type),intent(in) :: cprj_kb(natom,dtefield%nspinor*mband)
  real(dp),intent(out) :: smat_k_paw(2,dtefield%nband_occ,dtefield%nband_occ)
  integer,intent(in) :: typat(natom)
 end subroutine smatrix_k_paw
end interface

interface
 subroutine smatrix_pawinit(atindx1,cm2,cprj,ikpt1,ikpt2,isppol,&  
  &  g1,gprimd,kpt,mband,mbandw,mkmem,mpi_enreg,&  
  &  natom,nband,nkpt,nspinor,nsppol,ntypat,pawang,pawrad,pawtab,rprimd,&  
  &  seed_name,typat,xred)
  use m_pawrad
  use m_pawang
  use defs_abitypes
  use m_pawcprj
  use m_pawtab
  use defs_basis
  implicit none
  integer,intent(in) :: ikpt1
  integer,intent(in) :: ikpt2
  integer,intent(in) :: isppol
  integer,intent(in) :: mband
  integer,intent(in) :: mbandw
  integer,intent(in) :: mkmem
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: ntypat
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  character(len=fnlen) :: seed_name
  integer,intent(in) :: g1(3)
  integer,intent(in) :: atindx1(natom)
  real(dp),intent(inout) :: cm2(2,mbandw,mbandw)
  type(pawcprj_type) :: cprj(natom,nspinor*mband*mkmem*nsppol)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: kpt(3,nkpt)
  integer,intent(in) :: nband(nsppol*nkpt)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine smatrix_pawinit
end interface

interface
 subroutine spline_paw_fncs(dphi,dtphi,nnl,npts,pawrad,pawtab,points,phi,tphi)
  use defs_basis
  use m_pawrad
  use m_pawtab
  implicit none
  integer,intent(in) :: nnl
  integer,intent(in) :: npts
  type(pawrad_type),intent(in) :: pawrad
  type(pawtab_type),intent(in) :: pawtab
  real(dp),intent(out) :: dphi(npts,nnl)
  real(dp),intent(out) :: dtphi(npts,nnl)
  real(dp),intent(out) :: phi(npts,nnl)
  real(dp),intent(in) :: points(npts)
  real(dp),intent(out) :: tphi(npts,nnl)
 end subroutine spline_paw_fncs
end interface

interface
 subroutine transgrid(cplex,mpi_enreg,nspden,optgrid,optin,optout,paral_kgb,pawfgr,rhog,rhogf,rhor,rhorf)
  use m_pawfgr
  use defs_abitypes
  use defs_basis
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: nspden
  integer,intent(in) :: optgrid
  integer,intent(in) :: optin
  integer,intent(in) :: optout
  integer,intent(in) :: paral_kgb
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawfgr_type),intent(in) :: pawfgr
  real(dp),intent(inout) :: rhog(2,pawfgr%nfftc)
  real(dp),intent(inout) :: rhogf(2,pawfgr%nfft)
  real(dp),intent(inout) :: rhor(cplex*pawfgr%nfftc,nspden)
  real(dp),intent(inout) :: rhorf(cplex*pawfgr%nfft,nspden)
 end subroutine transgrid
end interface

interface
 subroutine wvl_nhatgrid(atindx1,geocode,h,i3s,natom,natom_tot,&  
  &  nattyp,ntypat,n1,n1i,n2,n2i,n3,n3pi,optcut,optgr0,optgr1,optgr2,optrad,&  
  &  pawfgrtab,pawtab,psppar,rprimd,shift,xred)
  use defs_basis
  use m_pawfgrtab
  use m_pawtab
  implicit none
  integer,intent(in) :: i3s
  integer,intent(in) :: n1
  integer,intent(in) :: n1i
  integer,intent(in) :: n2
  integer,intent(in) :: n2i
  integer,intent(in) :: n3
  integer,intent(in) :: n3pi
  integer,intent(in) :: natom
  integer,intent(in) :: natom_tot
  integer,intent(in) :: ntypat
  integer,intent(in) :: optcut
  integer,intent(in) :: optgr0
  integer,intent(in) :: optgr1
  integer,intent(in) :: optgr2
  integer,intent(in) :: optrad
  integer,intent(in) :: shift
  character(1),intent(in) :: geocode
  integer,intent(in) :: atindx1(natom)
  real(dp),intent(in) :: h(3)
  integer,intent(in) :: nattyp(ntypat)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(in) :: psppar(0:4,0:6,ntypat)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(inout) :: xred(3,natom)
 end subroutine wvl_nhatgrid
end interface

end module interfaces_65_paw
!!***
