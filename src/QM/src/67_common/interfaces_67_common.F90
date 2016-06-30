!!****m* ABINIT/interfaces_67_common
!! NAME
!! interfaces_67_common
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/67_common
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

module interfaces_67_common

 implicit none

interface
 subroutine berryphase(atindx1,bdberry,cg,gprimd,istwfk,kberry,kg,kpt_,&  
  &  kptopt,kptrlatt,mband,mcg,&  
  &  mkmem,mpw,natom,nattyp,nband,nberry,npwarr,nspinor,nsppol,ntypat,&  
  &  nkpt_,rprimd,ucvol,xred,zion)
  use defs_basis
  implicit none
  integer,intent(in) :: kptopt
  integer,intent(in) :: mband
  integer,intent(in) :: mcg
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: natom
  integer,intent(in) :: nberry
  integer,intent(in) :: nkpt_
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: ntypat
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: bdberry(4)
  integer,intent(in) :: kptrlatt(3,3)
  integer,intent(in) :: atindx1(natom)
  real(dp),intent(in) :: cg(2,mcg)
  real(dp),intent(in) :: gprimd(1:3,1:3)
  integer,intent(in) :: istwfk(nkpt_)
  integer,intent(in) :: kberry(3,nberry)
  integer,intent(in) :: kg(3,mpw*mkmem)
  real(dp),intent(in) :: kpt_(3,nkpt_)
  integer,intent(in) :: nattyp(ntypat)
  integer,intent(in) :: nband(nkpt_*nsppol)
  integer,intent(in) :: npwarr(nkpt_)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: zion(ntypat)
 end subroutine berryphase
end interface

interface
 subroutine berryphase_new(atindx1,cg,cprj,dtefield,dtfil,dtset,psps,&  
  &  gprimd,hdr,indlmn,kg,lmnmax,mband,mcg,mcprj,&  
  &  mkmem,mpi_enreg,mpw,my_natom,natom,npwarr,nsppol,ntypat,&  
  &  nkpt,option,pawrhoij,pawtab,pel,pelev,pion,ptot,red_ptot,pwind,&  !!REC
  &  pwind_alloc,pwnsfac,&  
  &  rprimd,typat,ucvol,unit_out,usecprj,usepaw,xred,zion)
  use m_pawrhoij
  use defs_abitypes
  use m_pawcprj
  use m_pawtab
  use defs_basis
  use m_efield
  use defs_datatypes
  implicit none
  integer, intent(in) :: lmnmax
  integer, intent(in) :: mband
  integer, intent(in) :: mcg
  integer, intent(in) :: mcprj
  integer, intent(in) :: mkmem
  integer, intent(in) :: mpw
  integer, intent(in) :: my_natom
  integer, intent(in) :: natom
  integer, intent(in) :: nkpt
  integer, intent(in) :: nsppol
  integer, intent(in) :: ntypat
  integer, intent(in) :: option
  integer, intent(in) :: pwind_alloc
  integer, intent(in) :: unit_out
  integer, intent(in) :: usecprj
  integer, intent(in) :: usepaw
  type(efield_type), intent(inout) :: dtefield
  type(datafiles_type), intent(in) :: dtfil
  type(dataset_type), intent(in) :: dtset
  type(hdr_type), intent(inout) :: hdr
  type(mpi_type), intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  real(dp), intent(in) :: ucvol
  integer, intent(in) :: atindx1(natom)
  real(dp), intent(in) :: cg(2,mcg)
  type(pawcprj_type),intent(in) :: cprj(natom,mcprj*usecprj)
  real(dp), intent(in) :: gprimd(3,3)
  integer, intent(in) :: indlmn(6,lmnmax,ntypat)
  integer, intent(in) :: kg(3,mpw*mkmem)
  integer, intent(in) :: npwarr(nkpt)
  type(pawrhoij_type), intent(in) :: pawrhoij(my_natom*usepaw)
  type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)
  real(dp), intent(out) :: pel(3)
  real(dp), intent(out) :: pelev(3)
  real(dp), intent(out) :: pion(3)
  real(dp), intent(out) :: ptot(3)
  integer, intent(in) :: pwind(pwind_alloc,2,3)
  real(dp), intent(in) :: pwnsfac(2,pwind_alloc)
  real(dp), intent(out) :: red_ptot(3)
  real(dp), intent(in) :: rprimd(3,3)
  integer, intent(in) :: typat(natom)
  real(dp), intent(inout) :: xred(3,natom)
  real(dp), intent(in) :: zion(ntypat)
 end subroutine berryphase_new
end interface

interface
 subroutine calc_efg(mpi_enreg,my_natom,natom,nfft,ngfft,nspden,nsym,ntypat,paral_kgb,&  
  &  paw_an,pawang,pawrad,pawrhoij,pawtab,&  
  &  ptcharge,prtefg,quadmom,rhor,rprimd,symrel,tnons,typat,ucvol,usepaw,xred,zion,&  
  &  mpi_atmtab,comm_atom) ! optional arguments (parallelism)
  use m_pawrad
  use m_pawang
  use m_pawrhoij
  use defs_abitypes
  use m_paw_an
  use defs_basis
  use m_pawtab
  implicit none
  integer,optional,intent(in) :: comm_atom
  integer,intent(in) :: my_natom
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: prtefg
  integer,intent(in) :: usepaw
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  real(dp),intent(in) :: ucvol
  integer,optional,target,intent(in) :: mpi_atmtab(:)
  integer,intent(in) :: ngfft(18)
  type(paw_an_type),intent(in) :: paw_an(my_natom)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawrhoij_type),intent(in) :: pawrhoij(my_natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(in) :: ptcharge(ntypat)
  real(dp),intent(in) :: quadmom(ntypat)
  real(dp),intent(in) :: rhor(nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
  integer,intent(in) :: typat(natom)
  real(dp),intent(inout) :: xred(3,natom)
  real(dp),intent(in) :: zion(ntypat)
 end subroutine calc_efg
end interface

interface
 subroutine calc_fc(my_natom,natom,nspden,ntypat,pawrad,pawrhoij,pawtab,typat,usepaw,&  
  &  mpi_atmtab,comm_atom) ! optional arguments (parallelism)
  use m_pawtab
  use m_pawrad
  use m_pawrhoij
  implicit none
  integer,optional,intent(in) :: comm_atom
  integer,intent(in) :: my_natom
  integer,intent(in) :: natom
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  integer,intent(in) :: usepaw
  integer,optional,target,intent(in) :: mpi_atmtab(:)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawrhoij_type),intent(in) :: pawrhoij(my_natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  integer,intent(in) :: typat(natom)
 end subroutine calc_fc
end interface

interface
 subroutine cgq_builder(berryflag,cg,cgq,dtefield,dtset,ikpt,ikpt_loc,isppol,mcg,mcgq,&  
  &  me_distrb,mkgq,mpi_enreg,my_nspinor,nband_k,nproc_distrb,&  
  &  npwarr,pwnsfac,pwnsfacq,pwind_alloc,spaceComm_distrb)
  use defs_basis
  use m_efield
  use defs_abitypes
  implicit none
  integer,intent(in) :: ikpt
  integer,intent(in) :: ikpt_loc
  integer,intent(in) :: isppol
  integer,intent(in) :: mcg
  integer,intent(in) :: mcgq
  integer,intent(in) :: me_distrb
  integer,intent(in) :: mkgq
  integer,intent(in) :: my_nspinor
  integer,intent(in) :: nband_k
  integer,intent(in) :: nproc_distrb
  integer,intent(in) :: pwind_alloc
  integer,intent(in) :: spaceComm_distrb
  logical,intent(in) :: berryflag
  type(efield_type), intent(inout) :: dtefield
  type(dataset_type), intent(in) :: dtset
  type(mpi_type), intent(inout) :: mpi_enreg
  real(dp),intent(in) :: cg(2,mcg)
  real(dp),intent(out) :: cgq(2,mcgq)
  integer,intent(in) :: npwarr(dtset%nkpt)
  real(dp),intent(in) :: pwnsfac(2,pwind_alloc)
  real(dp),intent(out) :: pwnsfacq(2,mkgq)
 end subroutine cgq_builder
end interface

interface
 subroutine cgwf(berryopt,cg,cgq,chkexit,cpus,dphase_k,dtefield,&  
  &  filnam_ds1,gsc,gs_hamk,icg,igsc,ikpt,inonsc,&  
  &  isppol,mband,mcg,mcgq,mgsc,mkgq,mpi_enreg,&  
  &  mpw,nband,nbdblock,nkpt,nline,npw,npwarr,&  
  &  nspinor,nsppol,ortalg,prtvol,pwind,&  
  &  pwind_alloc,pwnsfac,pwnsfacq,quit,resid,subham,subovl,&  
  &  subvnl,tolrde,tolwfr,use_subovl,wfoptalg,zshift)
  use defs_basis
  use m_efield
  use defs_abitypes
  use m_hamiltonian
  implicit none
  integer,intent(in) :: berryopt
  integer,intent(in) :: chkexit
  integer,intent(in) :: icg
  integer,intent(in) :: igsc
  integer,intent(in) :: ikpt
  integer,intent(in) :: inonsc
  integer,intent(in) :: isppol
  integer,intent(in) :: mband
  integer,intent(in) :: mcg
  integer,intent(in) :: mcgq
  integer,intent(in) :: mgsc
  integer,intent(in) :: mkgq
  integer,intent(in) :: mpw
  integer,intent(in) :: nband
  integer,intent(in) :: nbdblock
  integer,intent(in) :: nkpt
  integer,intent(in) :: nline
  integer,intent(in) :: npw
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: ortalg
  integer,intent(in) :: prtvol
  integer,intent(in) :: pwind_alloc
  integer,intent(in) :: quit
  integer,intent(in) :: use_subovl
  integer,intent(in) :: wfoptalg
  real(dp),intent(in) :: cpus
  type(efield_type),intent(inout) :: dtefield
  character(len=*),intent(in) :: filnam_ds1
  type(gs_hamiltonian_type),intent(inout) :: gs_hamk
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: tolrde
  real(dp),intent(in) :: tolwfr
  real(dp),intent(inout) :: cg(2,mcg)
  real(dp),intent(in) :: cgq(2,mcgq)
  real(dp),intent(inout) :: dphase_k(3)
  real(dp),intent(inout) :: gsc(2,mgsc)
  integer,intent(in) :: npwarr(nkpt)
  integer,intent(in) :: pwind(pwind_alloc,2,3)
  real(dp),intent(in) :: pwnsfac(2,pwind_alloc)
  real(dp),intent(in) :: pwnsfacq(2,mkgq)
  real(dp),intent(out) :: resid(nband)
  real(dp),intent(out) :: subham(nband*(nband+1))
  real(dp),intent(out) :: subovl(nband*(nband+1)*use_subovl)
  real(dp),intent(out) :: subvnl(nband*(nband+1)*(1-gs_hamk%usepaw))
  real(dp),intent(in) :: zshift(nband)
 end subroutine cgwf
end interface

interface
 subroutine clnup1(acell,dtset,eigen,fermie,&  
  &  fnameabo_dos,fnameabo_eig,fred,&  
  &  mpi_enreg,nfft,ngfft,occ,prtfor,&  
  &  resid,rhor,rprimd,vxcavg,xred)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: nfft
  integer,intent(in) :: prtfor
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: fermie
  character(len=*),intent(in) :: fnameabo_dos
  character(len=*),intent(in) :: fnameabo_eig
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: vxcavg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp),intent(in) :: fred(3,dtset%natom)
  real(dp),intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp),intent(in) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp),intent(in) :: rhor(nfft,dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: xred(3,dtset%natom)
 end subroutine clnup1
end interface

interface
 subroutine clnup2(n1xccc,fred,gresid,grewtn,grvdw,grxc,iscf,natom,ngrvdw,&  
  &  prtfor,prtstr,prtvol,start,strten,synlgr,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: iscf
  integer,intent(in) :: n1xccc
  integer,intent(in) :: natom
  integer,intent(in) :: ngrvdw
  integer,intent(in) :: prtfor
  integer,intent(in) :: prtstr
  integer,intent(in) :: prtvol
  real(dp),intent(in) :: fred(3,natom)
  real(dp),intent(in) :: gresid(3,natom)
  real(dp),intent(in) :: grewtn(3,natom)
  real(dp),intent(in) :: grvdw(3,ngrvdw)
  real(dp),intent(in) :: grxc(3,natom)
  real(dp),intent(in) :: start(3,natom)
  real(dp),intent(in) :: strten(6)
  real(dp),intent(in) :: synlgr(3,natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine clnup2
end interface

interface
 subroutine conducti_nc(filnam,filnam_out,mpi_enreg)
  use defs_basis
  use defs_abitypes
  implicit none
  character(len=fnlen) :: filnam
  character(len=fnlen) :: filnam_out
  type(mpi_type),intent(inout) :: mpi_enreg
 end subroutine conducti_nc
end interface

interface
 subroutine conducti_paw(filnam,filnam_out,mpi_enreg)
  use defs_basis
  use defs_abitypes
  implicit none
  character(len=fnlen) :: filnam
  character(len=fnlen) :: filnam_out
  type(mpi_type),intent(inout) :: mpi_enreg
 end subroutine conducti_paw
end interface

interface
 subroutine conducti_paw_core(filnam,filnam_out,mpi_enreg)
  use defs_basis
  use defs_abitypes
  implicit none
  character(len=fnlen) :: filnam
  character(len=fnlen) :: filnam_out
  type(mpi_type),intent(inout) :: mpi_enreg
 end subroutine conducti_paw_core
end interface

interface
 subroutine constrf(diffor,fcart,forold,fred,iatfix,ionmov,maxfor,natom,&  
  &  nconeq,prtvol,rprimd,wtatcon,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: ionmov
  integer,intent(in) :: natom
  integer,intent(in) :: nconeq
  integer,intent(in) :: prtvol
  real(dp),intent(out) :: diffor
  real(dp),intent(out) :: maxfor
  real(dp),intent(inout) :: fcart(3,natom)
  real(dp),intent(inout) :: forold(3,natom)
  real(dp),intent(inout) :: fred(3,natom)
  integer,intent(in) :: iatfix(3,natom)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: wtatcon(3,natom,nconeq)
  real(dp),intent(inout) :: xred(3,natom)
 end subroutine constrf
end interface

interface
 subroutine dens_in_sph(cmax,cg,gmet,istwfk,kg_k,natom,ngfft,mpi_enreg,npw_k,&  
  &  paral_kgb,ph1d,rmax,ucvol)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: istwfk
  integer,intent(in) :: natom
  integer,intent(in) :: npw_k
  integer,intent(in) :: paral_kgb
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  real(dp),intent(inout) :: cg(2,npw_k)
  real(dp),intent(out) :: cmax(natom)
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: kg_k(3,npw_k)
  real(dp),intent(in) :: ph1d(2,(2*ngfft(1)+1+2*ngfft(2)+1+2*ngfft(3)+1)*natom)
  real(dp),intent(in) :: rmax(natom)
 end subroutine dens_in_sph
end interface

interface
 subroutine dielmt(dielinv,gmet,kg_diel,&  
  &  npwdiel,nspden,occopt,prtvol,susmat)
  use defs_basis
  implicit none
  integer,intent(in) :: npwdiel
  integer,intent(in) :: nspden
  integer,intent(in) :: occopt
  integer,intent(in) :: prtvol
  real(dp),intent(out) :: dielinv(2,npwdiel,nspden,npwdiel,nspden)
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: kg_diel(3,npwdiel)
  real(dp),intent(in) :: susmat(2,npwdiel,nspden,npwdiel,nspden)
 end subroutine dielmt
end interface

interface
 subroutine dieltcel(dielinv,gmet,kg_diel,kxc,&  
  &  nfft,ngfft,nkxc,npwdiel,nspden,occopt,option,paral_kgb,prtvol,susmat)
  use defs_basis
  implicit none
  integer,intent(in) :: nfft
  integer,intent(in) :: nkxc
  integer,intent(in) :: npwdiel
  integer,intent(in) :: nspden
  integer,intent(in) :: occopt
  integer,intent(in) :: option
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: prtvol
  integer,intent(in) :: ngfft(18)
  real(dp),intent(out) :: dielinv(2,npwdiel,nspden,npwdiel,nspden)
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: kg_diel(3,npwdiel)
  real(dp),intent(in) :: kxc(nfft,nkxc)
  real(dp),intent(in) :: susmat(2,npwdiel,nspden,npwdiel,nspden)
 end subroutine dieltcel
end interface

interface
 subroutine emispec_paw(filnam,filnam_out,mpi_enreg)
  use defs_basis
  use defs_abitypes
  implicit none
  character(len=fnlen) :: filnam
  character(len=fnlen) :: filnam_out
  type(mpi_type),intent(inout) :: mpi_enreg
 end subroutine emispec_paw
end interface

interface
 subroutine energy(cg,compch_fft,dtset,electronpositron,&  
  &  energies,eigen,etotal,gsqcut,indsym,irrzon,kg,mcg,mpi_enreg,my_natom,nfftf,ngfftf,nhat,&  
  &  nhatgr,nhatgrdim,npwarr,n3xccc,occ,optene,paw_dmft,paw_ij,pawang,pawfgr,&  
  &  pawfgrtab,pawrhoij,pawtab,phnons,ph1d,psps,resid,rhog,rhor,rprimd,strsxc,symrec,&  
  &  taug,taur,usexcnhat,vhartr,vtrial,vpsp,vxc,vxctau,wfs,wvl,wvl_den,wvl_e,xccc3d,xred,ylm,&  
  &  add_tfw) ! optional argument
  use defs_wvltypes
  use m_pawrhoij
  use m_pawtab
  use m_paw_ij
  use m_pawang
  use m_paw_dmft
  use m_energies
  use defs_abitypes
  use m_pawfgrtab
  use m_pawfgr
  use defs_basis
  use defs_datatypes
  use m_electronpositron
  implicit none
  integer,intent(in) :: mcg
  integer,intent(in) :: my_natom
  integer,intent(in) :: n3xccc
  integer,intent(in) :: nfftf
  integer,intent(in) :: nhatgrdim
  integer,intent(in) :: optene
  integer,intent(in) :: usexcnhat
  logical,intent(in),optional :: add_tfw
  real(dp),intent(out) :: compch_fft
  type(dataset_type),intent(in) :: dtset
  type(electronpositron_type),pointer :: electronpositron
  type(energies_type),intent(inout) :: energies
  real(dp),intent(out) :: etotal
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  type(paw_dmft_type), intent(inout) :: paw_dmft
  type(pawang_type),intent(in) :: pawang
  type(pawfgr_type),intent(in) :: pawfgr
  type(pseudopotential_type),intent(in) :: psps
  type(wvl_wf_type),intent(inout) :: wfs
  type(wvl_internal_type), intent(in) :: wvl
  type(wvl_denspot_type), intent(inout) :: wvl_den
  type(wvl_energy_terms),intent(inout) :: wvl_e
  integer, intent(in) :: ngfftf(18)
  real(dp), intent(in) :: cg(2,mcg)
  real(dp), intent(in) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
  integer, intent(in) :: indsym(4,dtset%nsym,dtset%natom)
  integer :: irrzon(dtset%nfft**(1-1/dtset%nsym),2,(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  integer :: kg(3,dtset%mpw*dtset%mkmem)
  real(dp), intent(inout) :: nhat(nfftf,dtset%nspden*psps%usepaw)
  real(dp),intent(in) :: nhatgr(nfftf,dtset%nspden,3*nhatgrdim)
  integer, intent(in) :: npwarr(dtset%nkpt)
  real(dp), intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  type(paw_ij_type), intent(in) :: paw_ij(my_natom*psps%usepaw)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(my_natom)
  type(pawrhoij_type),target,intent(inout) :: pawrhoij(my_natom*psps%usepaw)
  type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)
  real(dp), intent(in) :: ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom)
  real(dp), intent(in) :: phnons(2,dtset%nfft**(1-1/dtset%nsym), &
  &         (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  real(dp), intent(inout) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp), intent(inout) :: rhog(2,nfftf)
  real(dp), intent(inout) :: rhor(nfftf,dtset%nspden)
  real(dp), intent(in) :: rprimd(3,3)
  real(dp), intent(out) :: strsxc(6)
  integer, intent(in) :: symrec(3,3,dtset%nsym)
  real(dp), intent(inout) :: taug(2,nfftf*dtset%usekden)
  real(dp), intent(inout) :: taur(nfftf,dtset%nspden*dtset%usekden)
  real(dp), intent(out) :: vhartr(nfftf)
  real(dp), intent(in) :: vpsp(nfftf)
  real(dp), intent(out) :: vtrial(nfftf,dtset%nspden)
  real(dp), intent(out) :: vxc(nfftf,dtset%nspden)
  real(dp),intent(out) :: vxctau(nfftf,dtset%nspden,4*dtset%usekden)
  real(dp), intent(in) :: xccc3d(n3xccc)
  real(dp), intent(in) :: xred(3,dtset%natom)
  real(dp), intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 end subroutine energy
end interface

interface
 subroutine erlxconv(etotals,iexit,ihist,itime,mxhist,ntime,tolmxde)
  use defs_basis
  implicit none
  integer,intent(inout) :: iexit
  integer,intent(in) :: ihist
  integer,intent(in) :: itime
  integer,intent(in) :: mxhist
  integer,intent(in) :: ntime
  real(dp), intent(in) :: tolmxde
  real(dp),intent(in) :: etotals(mxhist)
 end subroutine erlxconv
end interface

interface
 subroutine etheta(bcut,chc,detovc,detovd,dhc,dhd,efield_dot,e0,e1,&  
  &  hel,nkpt,nstr,sdeg,theta)
  use defs_basis
  implicit none
  integer,intent(in) :: nkpt
  real(dp),intent(in) :: chc
  real(dp),intent(in) :: dhc
  real(dp),intent(in) :: dhd
  real(dp),intent(out) :: e0
  real(dp),intent(out) :: e1
  real(dp),intent(in) :: sdeg
  real(dp),intent(in) :: theta
  integer,intent(in) :: hel(2,3)
  integer,intent(in) :: nstr(3)
  real(dp),intent(in) :: bcut(2,3)
  real(dp),intent(in) :: detovc(2,2,3)
  real(dp),intent(in) :: detovd(2,2,3)
  real(dp),intent(in) :: efield_dot(3)
 end subroutine etheta
end interface

interface
 subroutine etotfor(atindx1,deltae,diffor,dtefield,dtset,&  
  &  elast,electronpositron,energies,&  
  &  etotal,favg,fcart,forold,fred,gresid,grewtn,grhf,grnl,grvdw,&  
  &  grxc,gsqcut,indsym,kxc,maxfor,mgfft,mpi_enreg,my_natom,nattyp,&  
  &  nfft,ngfft,ngrvdw,nhat,nkxc,ntypat,nvresid,n1xccc,n3xccc,optene,optforces,optres,&  
  &  pawang,pawfgrtab,pawrhoij,pawtab,ph1d,red_ptot,psps,rhog,rhor,rprimd,symrec,synlgr,&  
  &  usepaw,vhartr,vpsp,vxc,wvl,wvl_den,xccc3d,xred)
  use defs_wvltypes
  use m_pawang
  use m_pawrhoij
  use m_energies
  use defs_abitypes
  use m_pawfgrtab
  use m_pawtab
  use defs_basis
  use m_efield
  use defs_datatypes
  use m_electronpositron
  implicit none
  integer,intent(in) :: mgfft
  integer,intent(in) :: my_natom
  integer,intent(in) :: n1xccc
  integer,intent(in) :: n3xccc
  integer,intent(in) :: nfft
  integer,intent(in) :: ngrvdw
  integer,intent(in) :: nkxc
  integer,intent(in) :: ntypat
  integer,intent(in) :: optene
  integer,intent(in) :: optforces
  integer,intent(in) :: optres
  integer,intent(in) :: usepaw
  real(dp),intent(out) :: deltae
  real(dp),intent(out) :: diffor
  type(efield_type),intent(in) :: dtefield
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(inout) :: elast
  type(electronpositron_type),pointer :: electronpositron
  type(energies_type),intent(inout) :: energies
  real(dp),intent(out) :: etotal
  real(dp),intent(in) :: gsqcut
  real(dp),intent(out) :: maxfor
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pseudopotential_type),intent(in) :: psps
  type(wvl_internal_type), intent(in) :: wvl
  type(wvl_denspot_type), intent(inout) :: wvl_den
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: atindx1(dtset%natom)
  real(dp),intent(out) :: favg(3)
  real(dp),intent(inout) :: fcart(3,dtset%natom)
  real(dp),intent(inout) :: forold(3,dtset%natom)
  real(dp),intent(out) :: fred(3,dtset%natom)
  real(dp),intent(out) :: gresid(3,dtset%natom)
  real(dp),intent(in) :: grewtn(3,dtset%natom)
  real(dp),intent(out) :: grhf(3,dtset%natom)
  real(dp),intent(inout) :: grnl(3*dtset%natom)
  real(dp),intent(in) :: grvdw(3,ngrvdw)
  real(dp),intent(inout) :: grxc(3,dtset%natom)
  integer,intent(in) :: indsym(4,dtset%nsym,dtset%natom)
  real(dp),intent(in) :: kxc(nfft,nkxc)
  integer,intent(in) :: nattyp(ntypat)
  real(dp),intent(inout) :: nhat(nfft,dtset%nspden*psps%usepaw)
  real(dp),intent(inout) :: nvresid(nfft,dtset%nspden)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(my_natom*psps%usepaw)
  type(pawrhoij_type),intent(inout) :: pawrhoij(my_natom*psps%usepaw)
  type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*dtset%natom)
  real(dp),intent(in) :: red_ptot(3)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(in) :: rhor(nfft,dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symrec(3,3,dtset%nsym)
  real(dp),intent(out) :: synlgr(3,dtset%natom)
  real(dp),intent(in) :: vhartr(nfft)
  real(dp),intent(in) :: vpsp(nfft)
  real(dp),intent(in) :: vxc(nfft,dtset%nspden)
  real(dp),intent(in) :: xccc3d(n3xccc)
  real(dp),intent(inout) :: xred(3,dtset%natom)
 end subroutine etotfor
end interface

interface
 subroutine evdw_wannier(csix,corrvdw,origmwan,natom,nsppol,orignwan,tdocc_wan,vdw_nfrag,&  
  &  vdw_supercell,vdw_typfrag,vdw_xc,rprimd,wann_centres,wann_spreads,xcart)    
  use defs_basis
  implicit none
  integer , intent(in) :: natom
  integer , intent(in) :: nsppol
  integer , intent(in) :: origmwan
  integer , intent(in) :: vdw_nfrag
  integer , intent(in) :: vdw_xc
  real(dp), intent(out) :: corrvdw
  integer , intent(in) :: vdw_supercell(3)
  real(dp), intent(out) :: csix(origmwan,origmwan,nsppol,nsppol)
  integer , intent(in) :: orignwan(nsppol)
  real(dp), intent(in) :: rprimd(3,3)
  real(dpc), intent(in) :: tdocc_wan(origmwan,nsppol)
  integer , intent(in) :: vdw_typfrag(natom)
  real(dp), intent(in) :: wann_centres(3,origmwan,nsppol)
  real(dp), intent(in) :: wann_spreads(origmwan,nsppol)
  real(dp), intent(in) :: xcart(3,natom)
 end subroutine evdw_wannier
end interface

interface
 subroutine getFu(sn,sl,rn,rl,occn,occl,fu) ! sn-->spread(n), sl-->spread(l), rn --> rc(n), rl --> rc(l) 
  use defs_basis
  implicit none
  real(dp),intent(out) :: fu
  real(dp),intent(in) :: occl
  real(dp),intent(in) :: occn
  real(dp),intent(in) :: rl
  real(dp),intent(in) :: rn
  real(dp),intent(in) :: sl
  real(dp),intent(in) :: sn
 end subroutine getFu
end interface

interface
 subroutine order_wannier(mwan,natom,nwan,nsppol,ord,vdw_typfrag,wanncent,xcart)
  use defs_basis
  implicit none
  integer, intent(in) :: mwan
  integer, intent(in) :: natom
  integer, intent(in) :: nsppol
  integer, intent(in) :: nwan(nsppol)
  integer, intent(inout) :: ord(mwan,nsppol)
  integer, intent(in) :: vdw_typfrag(natom)
  real(dp),intent(in) :: wanncent(3,mwan,nsppol)
  real(dp),intent(in) :: xcart(3,natom)
 end subroutine order_wannier
end interface

interface
 subroutine ovlp_wann(mwan,nwan,nsppol,ord,wanncent,wannspr,xi)
  use defs_basis
  implicit none
  integer, intent(in) :: mwan
  integer, intent(in) :: nsppol
  integer, intent(in) :: nwan(nsppol)
  integer, intent(in) :: ord(mwan,nsppol)
  real(dp),intent(in) :: wanncent(3,mwan,nsppol)
  real(dp),intent(in) :: wannspr(mwan,nsppol)
  real(dp), intent(out) :: xi(mwan,nsppol)
 end subroutine ovlp_wann
end interface

interface
 subroutine vv10limit(sn,sl,rn,rl,fu) ! sn-->spread(n), sl-->spread(l), rn --> rc(n), rl --> rc(l) 
  use defs_basis
  implicit none
  real(dp),intent(out) :: fu
  real(dp),intent(in) :: rl
  real(dp),intent(in) :: rn
  real(dp),intent(in) :: sl
  real(dp),intent(in) :: sn
 end subroutine vv10limit
end interface

interface
 subroutine amalgam(amagr,ngr,nsppol,nw,mwan,ord,nwan,vdw_nfrag,wanncent,wannspr) 
  use defs_basis
  implicit none
  integer,intent(in) :: mwan
  integer,intent(out) :: ngr
  integer,intent(in) :: nsppol
  integer,intent(in) :: vdw_nfrag
  integer,intent(out) :: amagr(mwan,nsppol,mwan/2)
  integer,intent(out) :: nw(nsppol,mwan/2)
  integer,intent(in) :: nwan(nsppol)
  integer,intent(in) :: ord(mwan,nsppol)
  real(dp),intent(in) :: wanncent(3,mwan,nsppol)
  real(dp),intent(in) :: wannspr(mwan,nsppol)
 end subroutine amalgam
end interface

interface
 subroutine extraprho(atindx,atindx1,cg,dtset,gmet,gprimd,gsqcut,istep,&  
  &  kg,mcg,mgfft,mpi_enreg,mqgrid,my_natom,nattyp,nfft,ngfft,npwarr,ntypat,pawrhoij,&  
  &  pawtab,ph1d,psps,qgrid,rhor,rprimd,scf_history,ucvol,usepaw,&  
  &  xred_new,xred_old,ylm,zion,znucl)
  use m_pawrhoij
  use m_scf_history
  use defs_abitypes
  use m_pawtab
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: istep
  integer,intent(in) :: mcg
  integer,intent(in) :: mgfft
  integer,intent(in) :: mqgrid
  integer,intent(in) :: my_natom
  integer,intent(in) :: nfft
  integer,intent(in) :: ntypat
  integer,intent(in) :: usepaw
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  type(scf_history_type),intent(inout) :: scf_history
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(in) :: atindx1(dtset%natom)
  real(dp), intent(inout) :: cg(2,mcg)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
  integer,intent(in) :: nattyp(ntypat)
  integer,intent(in) :: npwarr(dtset%nkpt)
  type(pawrhoij_type),intent(inout) :: pawrhoij(my_natom*usepaw)
  type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*dtset%natom)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(inout) :: rhor(nfft,dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(inout) :: xred_new(3,dtset%natom)
  real(dp),intent(in) :: xred_old(3,dtset%natom)
  real(dp),intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
  real(dp),intent(in) :: zion(ntypat)
  real(dp),intent(in) :: znucl(ntypat)
 end subroutine extraprho
end interface

interface
 subroutine extrapwf(atindx,atindx1,cg,dtset,istep,kg,mcg,mgfft,mpi_enreg,&  
  &  nattyp,ngfft,npwarr,ntypat,pawtab,psps,rprimd,scf_history,usepaw,xred_old,ylm)
  use m_scf_history
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use m_pawtab
  implicit none
  integer,intent(in) :: istep
  integer,intent(in) :: mcg
  integer,intent(in) :: mgfft
  integer,intent(in) :: ntypat
  integer,intent(in) :: usepaw
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  type(scf_history_type),intent(inout) :: scf_history
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(in) :: atindx1(dtset%natom)
  real(dp), intent(inout) :: cg(2,mcg)
  integer,intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
  integer,intent(in) :: nattyp(ntypat)
  integer,intent(in) :: npwarr(dtset%nkpt)
  type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: xred_old(3,dtset%natom)
  real(dp),intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 end subroutine extrapwf
end interface

interface
 subroutine fconv(fcart,iatfix,iexit,itime,natom,ntime,optcell,strfact,strtarget,strten,tolmxf)
  use defs_basis
  implicit none
  integer,intent(inout) :: iexit
  integer,intent(in) :: itime
  integer,intent(in) :: natom
  integer,intent(in) :: ntime
  integer,intent(in) :: optcell
  real(dp),intent(in) :: strfact
  real(dp),intent(in) :: tolmxf
  real(dp),intent(in) :: fcart(3,natom)
  integer,intent(in) :: iatfix(3,natom)
  real(dp),intent(in) :: strtarget(6)
  real(dp),intent(in) :: strten(6)
 end subroutine fconv
end interface

interface
 subroutine fixsym(iatfix,indsym,natom,nsym)
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nsym
  integer,intent(in) :: iatfix(3,natom)
  integer,intent(in) :: indsym(4,nsym,natom)
 end subroutine fixsym
end interface

interface
 subroutine forces(atindx1,diffor,dtefield,dtset,favg,fcart,forold,fred,gresid,grewtn,&  
  &  grhf,grnl,grvdw,grxc,gsqcut,indsym,&  
  &  maxfor,mgfft,mpi_enreg,n1xccc,n3xccc,&  
  &  nattyp,nfft,ngfft,ngrvdw,ntypat,&  
  &  pawtab,ph1d,psps,rhog,rhor,rprimd,symrec,synlgr,&  
  &  vresid,vxc,wvl,wvl_den,xred,&  
  &  electronpositron) ! optional argument
  use defs_abitypes
  use m_electronpositron
  use m_pawtab
  use defs_basis
  use m_efield
  use defs_datatypes
  use defs_wvltypes
  implicit none
  integer,intent(in) :: mgfft
  integer,intent(in) :: n1xccc
  integer,intent(in) :: n3xccc
  integer,intent(in) :: nfft
  integer,intent(in) :: ngrvdw
  integer,intent(in) :: ntypat
  real(dp),intent(out) :: diffor
  type(efield_type),intent(in) :: dtefield
  type(dataset_type),intent(in) :: dtset
  type(electronpositron_type),pointer,optional :: electronpositron
  real(dp),intent(in) :: gsqcut
  real(dp),intent(out) :: maxfor
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  type(wvl_internal_type), intent(in) :: wvl
  type(wvl_denspot_type), intent(inout) :: wvl_den
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: atindx1(dtset%natom)
  real(dp),intent(out) :: favg(3)
  real(dp),intent(inout) :: fcart(3,dtset%natom)
  real(dp),intent(inout) :: forold(3,dtset%natom)
  real(dp),intent(out) :: fred(3,dtset%natom)
  real(dp),intent(out) :: gresid(3,dtset%natom)
  real(dp),intent(in) :: grewtn(3,dtset%natom)
  real(dp),intent(out) :: grhf(3,dtset%natom)
  real(dp),intent(in) :: grnl(3*dtset%natom)
  real(dp),intent(in) :: grvdw(3,ngrvdw)
  real(dp),intent(inout) :: grxc(3,dtset%natom)
  integer,intent(in) :: indsym(4,dtset%nsym,dtset%natom)
  integer,intent(in) :: nattyp(ntypat)
  type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*dtset%natom)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(in) :: rhor(nfft,dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symrec(3,3,dtset%nsym)
  real(dp),intent(out) :: synlgr(3,dtset%natom)
  real(dp),intent(inout) :: vresid(nfft,dtset%nspden)
  real(dp),intent(in) :: vxc(nfft,dtset%nspden)
  real(dp),intent(inout) :: xred(3,dtset%natom)
 end subroutine forces
end interface

interface
 subroutine forstr(atindx1,cg,cprj,diffor,dtefield,dtset,eigen,electronpositron,energies,favg,fcart,&  
  &  forold,fred,gresid,grewtn,grhf,grvdw,grxc,gsqcut,indsym,&  
  &  kg,kxc,maxfor,mcg,mcprj,mgfftf,mpi_enreg,my_natom,n3xccc,nattyp,&  
  &  nfftf,ngfftf,ngrvdw,nhat,nkxc,npwarr,&  
  &  ntypat,nvresid,occ,optfor,optres,paw_ij,pawang,pawfgr,&  
  &  pawfgrtab,pawrhoij,pawtab,ph1d,ph1df,psps,rhog,rhor,rprimd,stress_needed,&  
  &  strsxc,strten,symrec,synlgr,ucvol,usecprj,vhartr,vpsp,&  
  &  vxc,wvl,wvl_den,xccc3d,xred,ylm,ylmgr,qvpotzero)
  use defs_wvltypes
  use m_pawtab
  use m_paw_ij
  use m_pawang
  use m_pawrhoij
  use m_energies
  use defs_abitypes
  use m_pawcprj
  use m_pawfgrtab
  use m_pawfgr
  use defs_basis
  use m_efield
  use defs_datatypes
  use m_electronpositron
  implicit none
  integer,intent(in) :: mcg
  integer,intent(in) :: mcprj
  integer,intent(in) :: mgfftf
  integer,intent(in) :: my_natom
  integer,intent(in) :: n3xccc
  integer,intent(in) :: nfftf
  integer,intent(in) :: ngrvdw
  integer,intent(in) :: nkxc
  integer,intent(in) :: ntypat
  integer,intent(in) :: optfor
  integer,intent(in) :: optres
  integer,intent(in) :: stress_needed
  integer,intent(in) :: usecprj
  real(dp),intent(inout) :: diffor
  type(efield_type),intent(in) :: dtefield
  type(dataset_type),intent(in) :: dtset
  type(electronpositron_type),pointer :: electronpositron
  type(energies_type),intent(in) :: energies
  real(dp),intent(in) :: gsqcut
  real(dp),intent(inout) :: maxfor
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pawfgr_type),intent(in) :: pawfgr
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: qvpotzero
  real(dp),intent(in) :: ucvol
  type(wvl_internal_type), intent(in) :: wvl
  type(wvl_denspot_type), intent(inout) :: wvl_den
  integer,intent(in) :: ngfftf(18)
  integer,intent(in) :: atindx1(dtset%natom)
  real(dp),intent(in) :: cg(2,mcg)
  type(pawcprj_type),intent(in) :: cprj(dtset%natom,mcprj*usecprj)
  real(dp),intent(in) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp),intent(out) :: favg(3)
  real(dp),intent(inout) :: fcart(3,dtset%natom)
  real(dp),intent(inout) :: forold(3,dtset%natom)
  real(dp),intent(inout) :: fred(3,dtset%natom)
  real(dp),intent(inout) :: gresid(3,dtset%natom)
  real(dp),intent(in) :: grewtn(3,dtset%natom)
  real(dp),intent(inout) :: grhf(3,dtset%natom)
  real(dp),intent(in) :: grvdw(3,ngrvdw)
  real(dp),intent(inout) :: grxc(3,dtset%natom)
  integer,intent(in) :: indsym(4,dtset%nsym,dtset%natom)
  integer,intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
  real(dp),intent(in) :: kxc(dtset%nfft,nkxc)
  integer,intent(in) :: nattyp(ntypat)
  real(dp),intent(inout) :: nhat(nfftf,dtset%nspden*psps%usepaw)
  integer,intent(in) :: npwarr(dtset%nkpt)
  real(dp),intent(inout) :: nvresid(nfftf,dtset%nspden)
  real(dp),intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  type(paw_ij_type),intent(in) :: paw_ij(my_natom*psps%usepaw)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(my_natom*psps%usepaw)
  type(pawrhoij_type),intent(inout) :: pawrhoij(my_natom*psps%usepaw)
  type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)
  real(dp),intent(in) :: ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom)
  real(dp),intent(in) :: ph1df(2,3*(2*mgfftf+1)*dtset%natom)
  real(dp),intent(in) :: rhog(2,nfftf)
  real(dp),intent(inout) :: rhor(nfftf,dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: strsxc(6)
  real(dp),intent(inout) :: strten(6)
  integer,intent(in) :: symrec(3,3,dtset%nsym)
  real(dp),intent(inout) :: synlgr(3,dtset%natom)
  real(dp),intent(in) :: vhartr(nfftf)
  real(dp),intent(in) :: vpsp(nfftf)
  real(dp),intent(in) :: vxc(nfftf,dtset%nspden)
  real(dp),intent(inout) :: xccc3d(n3xccc)
  real(dp),intent(inout) :: xred(3,dtset%natom)
  real(dp),intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
  real(dp),intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 end subroutine forstr
end interface

interface
 subroutine forstrnps(cg,cprj,ecut,ecutsm,effmass,eigen,electronpositron,&  
  &  grnl,istwfk,kg,kinstr,npsstr,kpt,mband,mcg,mcprj,mgfft,mkmem,mpi_enreg,mpsang,&  
  &  mpw,my_natom,natom,nband,nfft,ngfft,nkpt,nloalg,npwarr,nspden,nspinor,nsppol,nsym,&  
  &  ntypat,nucdipmom,occ,optfor,paw_ij,pawtab,ph1d,psps,rprimd,&  
  &  stress_needed,symrec,typat,usecprj,use_gpu_cuda,wtk,xred,ylm,ylmgr)
  use m_paw_ij
  use defs_abitypes
  use m_pawcprj
  use m_pawtab
  use defs_basis
  use defs_datatypes
  use m_electronpositron
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: mcg
  integer,intent(in) :: mcprj
  integer,intent(in) :: mgfft
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpw
  integer,intent(in) :: my_natom
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspden
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(in) :: optfor
  integer,intent(in) :: stress_needed
  integer,intent(in) :: use_gpu_cuda
  integer,intent(in) :: usecprj
  real(dp),intent(in) :: ecut
  real(dp),intent(in) :: ecutsm
  real(dp),intent(in) :: effmass
  type(electronpositron_type),pointer :: electronpositron
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: nloalg(3)
  real(dp),intent(in) :: cg(2,mcg)
  type(pawcprj_type),intent(in) :: cprj(natom,mcprj*usecprj)
  real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
  real(dp),intent(out) :: grnl(3*natom*optfor)
  integer,intent(in) :: istwfk(nkpt)
  integer,intent(in) :: kg(3,mpw*mkmem)
  real(dp),intent(out) :: kinstr(6)
  real(dp),intent(in) :: kpt(3,nkpt)
  integer,intent(in) :: nband(nkpt*nsppol)
  real(dp),intent(out) :: npsstr(6)
  integer,intent(in) :: npwarr(nkpt)
  real(dp),intent(in) :: nucdipmom(3,my_natom)
  real(dp),intent(in) :: occ(mband*nkpt*nsppol)
  type(paw_ij_type),intent(in) :: paw_ij(my_natom*psps%usepaw)
  type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symrec(3,3,nsym)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: wtk(nkpt)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: ylm(mpw*mkmem,mpsang*mpsang*psps%useylm)
  real(dp),intent(in) :: ylmgr(mpw*mkmem,3,mpsang*mpsang*psps%useylm)
 end subroutine forstrnps
end interface

interface
 subroutine fresid(dtset,gresid,mpi_enreg,nfft,ngfft,ntypat,option,&  
  &  pawtab,rhor,rprimd,ucvol,work,xred_new,xred_old,znucl)
  use defs_basis
  use defs_abitypes
  use m_pawtab
  implicit none
  integer,intent(in) :: nfft
  integer,intent(in) :: ntypat
  integer,intent(in) :: option
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  real(dp),intent(out) :: gresid(3,dtset%natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat*dtset%usepaw)
  real(dp),intent(in) :: rhor(nfft,dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(inout) :: work(nfft,dtset%nspden)
  real(dp),intent(in) :: xred_new(3,dtset%natom)
  real(dp),intent(in) :: xred_old(3,dtset%natom)
  real(dp),intent(in) :: znucl(ntypat)
 end subroutine fresid
end interface

interface
 subroutine fresidrsp(atindx1,dtset,gmet,gprimd,gresid,gsqcut,mgfft,mpi_enreg,mqgrid,nattyp,nfft,&  
  &  ngfft,ntypat,psps,pawtab,ph1d,qgrid,ucvol,usepaw,vresid,zion,znucl)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use m_pawtab
  implicit none
  integer,intent(in) :: mgfft
  integer,intent(in) :: mqgrid
  integer,intent(in) :: nfft
  integer,intent(in) :: ntypat
  integer,intent(in) :: usepaw
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: atindx1(dtset%natom)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(out) :: gresid(3,dtset%natom)
  integer,intent(in) :: nattyp(ntypat)
  type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*dtset%natom)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(in) :: vresid(nfft,dtset%nspden)
  real(dp),intent(in) :: zion(ntypat)
  real(dp),intent(in) :: znucl(ntypat)
 end subroutine fresidrsp
end interface

interface
 subroutine getcgqphase(dtset, timrev, cg,  mcg,  cgq, mcgq, mpi_enreg,&  
  &  nkpt_rbz, npwarr, npwar1, phasecg)
  use defs_basis
  use defs_abitypes
  implicit none
  integer, intent(in) :: mcg
  integer, intent(in) :: mcgq
  integer, intent(in) :: nkpt_rbz
  integer, intent(in) :: timrev
  type(dataset_type), intent(in) :: dtset
  type(mpi_type), intent(in) :: mpi_enreg
  real(dp), intent(in) :: cg(2,mcg)
  real(dp), intent(in) :: cgq(2,mcgq)
  integer, intent(in) :: npwar1(nkpt_rbz)
  integer, intent(in) :: npwarr(nkpt_rbz)
  real(dp),intent(out) :: phasecg(2, dtset%mband*dtset%mband*nkpt_rbz*dtset%nsppol)
 end subroutine getcgqphase
end interface

interface
 subroutine hartre1(cplex,gmet,gsqcut,nfft,ngfft,paral_kgb,qphon,rhog,vhartr,ehvalues,rcut,ucvol)
  use defs_basis
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: nfft
  integer,intent(in) :: paral_kgb
  real(dp),intent(in) :: gsqcut
  real(dp),intent(in) :: rcut
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  real(dp),intent(out) :: ehvalues(3)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(out) :: vhartr(cplex*nfft)
 end subroutine hartre1
end interface

interface
 subroutine init_e_field_vars(dtefield,dtset,gmet,gprimd,kg,&  
  &  mpi_enreg,npwarr,occ,pawang,pawrad,pawtab,psps,&  
  &  pwind,pwind_alloc,pwnsfac,rprimd,symrec,xred)
  use m_pawrad
  use m_pawang
  use defs_abitypes
  use m_pawtab
  use defs_basis
  use m_efield
  use defs_datatypes
  implicit none
  integer,intent(out) :: pwind_alloc
  type(efield_type),intent(inout) :: dtefield
  type(dataset_type),intent(inout) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pseudopotential_type),intent(in) :: psps
  integer,pointer :: pwind(:,:,:)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
  integer,intent(in) :: npwarr(dtset%nkpt)
  real(dp),intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat*psps%usepaw)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*psps%usepaw)
  real(dp),pointer :: pwnsfac(:,:)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symrec(3,3,dtset%nsym)
  real(dp),intent(in) :: xred(3,dtset%natom)
 end subroutine init_e_field_vars
end interface

interface
 subroutine initberry(dtefield,dtset,gmet,gprimd,kg,mband,&  
  &  mkmem,mpi_enreg,mpw,natom,nkpt,npwarr,nsppol,&  
  &  nsym,ntypat,occ,pawang,pawrad,pawtab,psps,&  
  &  pwind,pwind_alloc,pwnsfac,&  
  &  rprimd,symrec,typat,usepaw,xred)
  use m_pawrad
  use m_pawang
  use defs_abitypes
  use m_pawtab
  use defs_basis
  use m_efield
  use defs_datatypes
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(out) :: pwind_alloc
  integer,intent(in) :: usepaw
  type(efield_type),intent(inout) :: dtefield
  type(dataset_type),intent(inout) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pseudopotential_type),intent(in) :: psps
  integer,pointer :: pwind(:,:,:)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: kg(3,mpw*mkmem)
  integer,intent(in) :: npwarr(nkpt)
  real(dp),intent(in) :: occ(mband*nkpt*nsppol)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),pointer :: pwnsfac(:,:)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symrec(3,3,nsym)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine initberry
end interface

interface
 subroutine initmv(cgindex,dtset,gmet,kg,kneigh,kg_neigh,kptindex,&  
  &  kpt3,mband,mkmem,mpi_enreg,mpw,nband,nkpt2,&  
  &  nkpt3,nneigh,npwarr,nsppol,occ,pwind)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: nkpt2
  integer,intent(in) :: nkpt3
  integer,intent(in) :: nneigh
  integer,intent(in) :: nsppol
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(out) :: cgindex(nkpt2,nsppol)
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: kg(3,mpw*mkmem)
  integer,intent(in) :: kg_neigh(30,nkpt2,3)
  integer,intent(in) :: kneigh(30,nkpt2)
  real(dp),intent(in) :: kpt3(3,nkpt3)
  integer,intent(in) :: kptindex(2,nkpt3)
  integer,intent(in) :: nband(nkpt2*nsppol)
  integer,intent(in) :: npwarr(nkpt2)
  real(dp),intent(in) :: occ(mband*nkpt2*nsppol)
  integer,intent(out) :: pwind(mpw,nneigh,mkmem)
 end subroutine initmv
end interface

interface
 subroutine initro(atindx,densty,gmet,gsqcut,izero,mgfft,mpi_enreg,mqgrid,natom,nattyp,&  
  &  nfft,ngfft,nspden,ntypat,paral_kgb,psps,pawtab,ph1d,qgrid,rhog,rhor,spinat,ucvol,usepaw,zion,znucl)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use m_pawtab
  implicit none
  integer,intent(in) :: izero
  integer,intent(in) :: mgfft
  integer,intent(in) :: mqgrid
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: usepaw
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(in) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: atindx(natom)
  real(dp),intent(in) :: densty(ntypat,4)
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: nattyp(ntypat)
  type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(out) :: rhog(2,nfft)
  real(dp),intent(out) :: rhor(nfft,nspden)
  real(dp),intent(in) :: spinat(3,natom)
  real(dp),intent(in) :: zion(ntypat)
  real(dp),intent(in) :: znucl(ntypat)
 end subroutine initro
end interface

interface
 subroutine ionion_realSpace(dtset, eew, grewtn, rprimd, xred, zion)
  use defs_basis
  use defs_abitypes
  implicit none
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(out) :: eew
  real(dp),intent(out) :: grewtn(3,dtset%natom)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: xred(3,dtset%natom)
  real(dp),intent(in) :: zion(dtset%ntypat)
 end subroutine ionion_realSpace
end interface

interface
 subroutine ionion_surface(dtset, eew, grewtn, me, nproc, rprimd, wvl, wvl_den, xred)
  use defs_basis
  use defs_abitypes
  use defs_wvltypes
  implicit none
  integer, intent(in) :: me
  integer, intent(in) :: nproc
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(out) :: eew
  type(wvl_internal_type), intent(in) :: wvl
  type(wvl_denspot_type), intent(inout) :: wvl_den
  real(dp),intent(out) :: grewtn(3,dtset%natom)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: xred(3,dtset%natom)
 end subroutine ionion_surface
end interface

interface
 subroutine jellium(gmet,gsqcut,mpi_enreg,nfft,ngfft,nspden,&  
  &  option,paral_kgb,slabwsrad,rhog,rhor,rprimd,vjell,slabzstart,slabzend)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: option
  integer,intent(in) :: paral_kgb
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: slabwsrad
  real(dp),intent(in) :: slabzend
  real(dp),intent(in) :: slabzstart
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(inout) :: rhog(2,nfft)
  real(dp),intent(inout) :: rhor(nfft,min(option,nspden))
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: vjell(nfft)
 end subroutine jellium
end interface

interface
 subroutine ks_ddiago(Diago_ctl,nband_k,nfftc,mgfftc,ngfftc,natom,&  
  &  typat,nfftf,nspinor,nspden,nsppol,Pawtab,Pawfgr,Paw_ij,&  
  &  Psps,rprimd,vtrial,xred,onband_diago,eig_ene,eig_vec,Cprj_k,comm,ierr,&  
  &  Electronpositron) ! Optional arguments
  use m_pawtab
  use m_paw_ij
  use m_pawcprj
  use m_pawfgr
  use m_hamiltonian
  use defs_basis
  use defs_datatypes
  use m_electronpositron
  implicit none
  integer,intent(in) :: comm
  integer,intent(out) :: ierr
  integer,intent(in) :: mgfftc
  integer,intent(in) :: natom
  integer,intent(in) :: nband_k
  integer,intent(in) :: nfftc
  integer,intent(in) :: nfftf
  integer,intent(in) :: nspden
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(out) :: onband_diago
  type(ddiago_ctl_type),intent(in) :: Diago_ctl
  type(electronpositron_type),optional,pointer :: Electronpositron
  type(pawfgr_type),intent(in) :: Pawfgr
  type(pseudopotential_type),intent(in) :: Psps
  integer,intent(in) :: ngfftc(18)
  type(pawcprj_type),pointer :: Cprj_k(:,:)
  type(paw_ij_type),intent(in) :: Paw_ij(natom*Psps%usepaw)
  type(pawtab_type),intent(in) :: Pawtab(Psps%ntypat*Psps%usepaw)
  real(dp),pointer :: eig_ene(:)
  real(dp),pointer :: eig_vec(:,:,:)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(inout) :: vtrial(nfftf,nspden)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine ks_ddiago
end interface

interface
 subroutine linemin(bcut,chc,costh,detovc,detovd,dhc,dhd,dphase_aux1,&  
  &  efield_dot,iline,nkpt,nstr,hel,phase_end,phase_init,sdeg,sinth,thetam)
  use defs_basis
  implicit none
  integer,intent(in) :: iline
  integer,intent(in) :: nkpt
  real(dp),intent(in) :: chc
  real(dp),intent(out) :: costh
  real(dp),intent(in) :: dhc
  real(dp),intent(in) :: dhd
  real(dp),intent(in) :: sdeg
  real(dp),intent(out) :: sinth
  real(dp),intent(out) :: thetam
  integer,intent(out) :: hel(2,3)
  integer,intent(in) :: nstr(3)
  real(dp),intent(out) :: bcut(2,3)
  real(dp),intent(in) :: detovc(2,2,3)
  real(dp),intent(in) :: detovd(2,2,3)
  real(dp),intent(inout) :: dphase_aux1(3)
  real(dp),intent(in) :: efield_dot(3)
  real(dp),intent(out) :: phase_end(3)
  real(dp),intent(inout) :: phase_init(3)
 end subroutine linemin
end interface

interface
 subroutine mag_constr(natom,spinat,nspden,magconon,magcon_lambda,rprimd,&  
  mpi_enreg,nfft,ngfft,ntypat,ratsph,rhor,&  
  typat,Vmagconstr,xred)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: magconon
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  real(dp),intent(in) :: magcon_lambda
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(out) :: Vmagconstr(nfft,nspden)
  real(dp),intent(in) :: ratsph(ntypat)
  real(dp),intent(in) :: rhor(nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: spinat(3,natom)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine mag_constr
end interface

interface
 subroutine mag_constr_e(magconon,magcon_lambda,mpi_enreg,natom,nfft,ngfft,nspden,ntypat,ratsph,rhor,rprimd,spinat,typat,xred)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: magconon
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  real(dp),intent(in) :: magcon_lambda
  type(mpi_type),intent(inout) :: mpi_enreg
  integer, intent(in) :: ngfft(18)
  real(dp),intent(in) :: ratsph(ntypat)
  real(dp),intent(in) :: rhor(nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: spinat(3,natom)
  integer, intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine mag_constr_e
end interface

interface
 subroutine make_efg_el(efg,mpi_enreg,natom,nfft,ngfft,nspden,nsym,paral_kgb,rhor,rprimd,symrel,tnons,xred)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: nsym
  integer,intent(in) :: paral_kgb
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(out) :: efg(3,3,natom)
  real(dp),intent(in) :: rhor(nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine make_efg_el
end interface

interface
 subroutine make_efg_ion(efg,natom,nsym,ntypat,rprimd,symrel,tnons,typat,ucvol,xred,zion)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  real(dp) :: ucvol
  real(dp),intent(out) :: efg(3,3,natom)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
  integer,intent(in) :: typat(natom)
  real(dp),intent(inout) :: xred(3,natom)
  real(dp),intent(in) :: zion(ntypat)
 end subroutine make_efg_ion
end interface

interface
 subroutine make_grad_berry(cg,cgq,cprj_k,detovc,dimlmn,dimlmn_srt,direc,dtefield,grad_berry,&  
  &  gs_hamk,iband,icg,ikpt,isppol,mband,mcg,mcgq,mkgq,mpi_enreg,mpw,natom,nkpt,&  
  &  npw,npwarr,nspinor,nsppol,pwind,pwind_alloc,pwnsfac,pwnsfacq)
  use defs_basis
  use defs_abitypes
  use m_efield
  use m_pawcprj
  use m_hamiltonian
  implicit none
  integer,intent(in) :: iband
  integer,intent(in) :: icg
  integer,intent(in) :: ikpt
  integer,intent(in) :: isppol
  integer,intent(in) :: mband
  integer,intent(in) :: mcg
  integer,intent(in) :: mcgq
  integer,intent(in) :: mkgq
  integer,intent(in) :: mpw
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt
  integer,intent(in) :: npw
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: pwind_alloc
  type(efield_type),intent(inout) :: dtefield
  type(gs_hamiltonian_type),intent(in) :: gs_hamk
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: cg(2,mcg)
  real(dp),intent(in) :: cgq(2,mcgq)
  type(pawcprj_type),intent(in) :: cprj_k(natom,dtefield%nband_occ*gs_hamk%usepaw*dtefield%nspinor)
  real(dp),intent(out) :: detovc(2,2,3)
  integer,intent(in) :: dimlmn(natom)
  integer,intent(in) :: dimlmn_srt(natom)
  real(dp),intent(inout) :: direc(2,npw*nspinor)
  real(dp),intent(out) :: grad_berry(2,npw*nspinor)
  integer,intent(in) :: npwarr(nkpt)
  integer,intent(in) :: pwind(pwind_alloc,2,3)
  real(dp),intent(in) :: pwnsfac(2,pwind_alloc)
  real(dp),intent(in) :: pwnsfacq(2,mkgq)
 end subroutine make_grad_berry
end interface

interface
 subroutine mkcore_inner(corfra,core_mesh,dyfrx2,&  
  &  grxc1,grxc2,grxc3,ifftsph,msz,natom,ncmax,nfft,nfgd,nfgd_r0,nspden,n3xccc,option,pawtab,&  
  &  rmet,rprimd,rr,rred,rshpm1,strdia,ucvol,vxc,xccc3d)
  use defs_basis
  use m_pawrad
  use m_pawtab
  implicit none
  integer , intent(in) :: msz
  integer , intent(in) :: n3xccc
  integer , intent(in) :: natom
  integer , intent(in) :: ncmax
  integer , intent(in) :: nfft
  integer , intent(in) :: nfgd
  integer , intent(in) :: nfgd_r0
  integer , intent(in) :: nspden
  integer , intent(in) :: option
  type(pawrad_type),intent(in) :: core_mesh
  real(dp), intent(out) :: grxc1
  real(dp), intent(out) :: grxc2
  real(dp), intent(out) :: grxc3
  type(pawtab_type),intent(in) :: pawtab
  real(dp),intent(in) :: rshpm1
  real(dp),intent(inout) :: strdia
  real(dp), intent(in) :: ucvol
  real(dp),intent(inout) :: corfra(3,3)
  real(dp),intent(inout) :: dyfrx2(3,3,natom)
  integer , intent(in) :: ifftsph(ncmax)
  real(dp), intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp), intent(in) :: rr(ncmax)
  real(dp),pointer :: rred(:,:)
  real(dp),intent(in) :: vxc(nfft,nspden)
  real(dp),intent(inout) :: xccc3d(n3xccc)
 end subroutine mkcore_inner
end interface

interface
 subroutine mkcore_paw(atindx1,corstr,dyfrx2,&  
  &  grxc,icoulomb,natom,mpi_enreg,&  
  &  nattyp,nfft,ngfft,nspden,ntypat,n3xccc,option,pawrad,pawtab,&  
  &  psppar,rprimd,&  
  &  ucvol,vxc,xccc3d,xred)
  use defs_basis
  use defs_abitypes
  use m_pawrad
  use m_pawtab
  implicit none
  integer,intent(in) :: icoulomb
  integer,intent(in) :: n3xccc
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  integer,intent(in) :: option
  type(mpi_type),intent(in) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: atindx1(natom)
  real(dp),intent(out) :: corstr(6)
  real(dp),intent(out) :: dyfrx2(3,3,natom)
  real(dp),intent(out) :: grxc(3,natom)
  integer,intent(in) :: nattyp(ntypat)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(in) :: psppar(0:4,0:6,ntypat)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: vxc(nfft,nspden)
  real(dp),intent(inout) :: xccc3d(nfft)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine mkcore_paw
end interface

interface
 subroutine mkcore_wvl(atindx1,corstr,dyfrx2,geocode,grxc,h,natom,&  
  &  nattyp,nfft,nscatterarr,nspden,ntypat,n1,n1i,n2,n2i,n3,n3pi,&  
  &  n3xccc,option,pawrad,pawtab,psppar,rprimd,ucvol,&  
  &  vxc,xccc3d,xred,mpi_comm_wvl)
  use defs_basis
  use m_pawrad
  use m_pawtab
  implicit none
  integer,intent(in),optional :: mpi_comm_wvl
  integer,intent(in) :: n1
  integer,intent(in) :: n1i
  integer,intent(in) :: n2
  integer,intent(in) :: n2i
  integer,intent(in) :: n3
  integer,intent(in) :: n3pi
  integer,intent(in) :: n3xccc
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  integer,intent(in) :: option
  character(1),intent(in) :: geocode
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: nscatterarr(4)
  integer,intent(in) :: atindx1(natom)
  real(dp),intent(out) :: corstr(6)
  real(dp),intent(out) :: dyfrx2(3,3,natom)
  real(dp),intent(out) :: grxc(3,natom)
  real(dp),intent(in) :: h(3)
  integer,intent(in) :: nattyp(ntypat)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(in) :: psppar(0:4,0:6,ntypat)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: vxc(nfft,nspden)
  real(dp),intent(out) :: xccc3d(n3xccc)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine mkcore_wvl
end interface

interface
 subroutine mkgrid_fft(ffti3_local,fftn3_distrib,gridcart,nfft,ngfft,rprimd)
  use defs_basis
  implicit none
  integer, intent(in) :: nfft
  integer, dimension(*), intent(in) :: ffti3_local
  integer, dimension(*), intent(in) :: fftn3_distrib
  integer,intent(in) :: ngfft(18)
  real(dp), dimension(3,nfft), intent(out) :: gridcart
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine mkgrid_fft
end interface

interface
 subroutine mklocl(dtset, dyfrlo,eei,gmet,gprimd,grtn,gsqcut,lpsstr,mgfft,&  
  &  mpi_enreg,natom,nattyp,nfft,ngfft,nspden,ntypat,option,pawtab,ph1d,psps,qprtrb,&  
  &  rhog,rhor,rprimd,ucvol,vprtrb,vpsp,wvl,wvl_den,xred)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use defs_wvltypes
  use m_pawtab
  implicit none
  integer,intent(in) :: mgfft
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  integer,intent(in) :: option
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: eei
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: ucvol
  type(wvl_internal_type), intent(in) :: wvl
  type(wvl_denspot_type), intent(inout) :: wvl_den
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: qprtrb(3)
  real(dp),intent(out) :: dyfrlo(3,3,natom)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(out) :: grtn(3,natom)
  real(dp),intent(out) :: lpsstr(6)
  integer,intent(in) :: nattyp(ntypat)
  type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(in) :: rhor(nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: vprtrb(2)
  real(dp),intent(inout) :: vpsp(nfft)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine mklocl
end interface

interface
 subroutine mklocl_realspace(grtn,icoulomb,mpi_enreg, natom, nattyp, nfft, ngfft,&  
  &  nscforder,nspden,ntypat,option,paral_kgb,pawtab,psps, rhog, rhor,&  
  &  rprimd, typat,ucvol,usewvl,vpsp, xred)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use m_pawtab
  implicit none
  integer,intent(in) :: icoulomb
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nscforder
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  integer,intent(in) :: option
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: usewvl
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  real(dp),intent(out) :: grtn(3,natom)
  integer,intent(in) :: nattyp(ntypat)
  type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(in) :: rhor(nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(out) :: vpsp(nfft)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine mklocl_realspace
end interface

interface
 subroutine createIonicPotential_new(fftn3_distrib,ffti3_local,geocode,iproc,&  
  &  nproc,nat,ntypes,iatype,psppar,nelpsp,rxyz,gridcart,&  
  &  hxh,hyh,hzh,n1i,n2i,n3d,n3i,kernel,pot_ion,spaceworld,pawtab,usepaw)
  use m_pawtab
  use defs_wvltypes
  implicit none
  integer, intent(in) :: iproc
  integer, intent(in) :: n1i
  integer, intent(in) :: n2i
  integer, intent(in) :: n3d
  integer, intent(in) :: n3i
  integer, intent(in) :: nat
  integer, intent(in) :: nproc
  integer, intent(in) :: ntypes
  integer, intent(in) :: spaceworld
  integer, intent(in) :: usepaw
  character(len=1), intent(in) :: geocode
  real(kind=8), intent(in) :: hxh
  real(kind=8), intent(in) :: hyh
  real(kind=8), intent(in) :: hzh
  type(coulomb_operator), intent(in) :: kernel
  integer, dimension(*), intent(in) :: ffti3_local
  integer, dimension(*), intent(in) :: fftn3_distrib
  real(kind=8), dimension(*), intent(inout) :: pot_ion
  real(kind=8), dimension(3,n1i*n2i*n3d), intent(in) :: gridcart
  integer, dimension(nat), intent(in) :: iatype
  integer, dimension(ntypes), intent(in) :: nelpsp
  type(pawtab_type),intent(in) :: pawtab(ntypes*usepaw)
  real(kind=8), dimension(0:4,0:6,ntypes), intent(in) :: psppar
  real(kind=8), dimension(3,nat), intent(in) :: rxyz
 end subroutine createIonicPotential_new
end interface

interface
 subroutine ind_positions_(periodic,i,n,j,go)
  implicit none
  integer, intent(in) :: i
  integer, intent(out) :: j
  integer, intent(in) :: n
  logical, intent(out) :: go
  logical, intent(in) :: periodic
 end subroutine ind_positions_
end interface

interface
 subroutine local_forces_new(fftn3_distrib,ffti3_local,&  
  geocode,iproc,ntypes,nat,iatype,rxyz,gridcart,psppar,nelpsp,hxh,hyh,hzh,&  
  n1,n2,n3,n3d,rho,pot,floc)
  implicit none
  integer, intent(in) :: iproc
  integer, intent(in) :: n1
  integer, intent(in) :: n2
  integer, intent(in) :: n3
  integer, intent(in) :: n3d
  integer, intent(in) :: nat
  integer, intent(in) :: ntypes
  character(len=1), intent(in) :: geocode
  real(kind=8), intent(in) :: hxh
  real(kind=8), intent(in) :: hyh
  real(kind=8), intent(in) :: hzh
  integer, dimension(*), intent(in) :: ffti3_local
  integer, dimension(*), intent(in) :: fftn3_distrib
  real(kind=8), dimension(*), intent(in) :: pot
  real(kind=8), dimension(*), intent(in) :: rho
  real(kind=8), dimension(3,nat), intent(out) :: floc
  real(kind=8), dimension(3,n1*n2*n3d), intent(in) :: gridcart
  integer, dimension(nat), intent(in) :: iatype
  integer, dimension(ntypes), intent(in) :: nelpsp
  real(kind=8), dimension(0:4,0:6,ntypes), intent(in) :: psppar
  real(kind=8), dimension(3,nat), intent(in) :: rxyz
 end subroutine local_forces_new
end interface

interface
 subroutine mklocl_recipspace(dyfrlo,eei,gmet,gprimd,grtn,gsqcut,lpsstr,mgfft,&  
  &  mpi_enreg,mqgrid,natom,nattyp,nfft,ngfft,ntypat,option,paral_kgb,ph1d,qgrid,qprtrb,&  
  &  rhog,ucvol,vlspl,vprtrb,vpsp)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: mgfft
  integer,intent(in) :: mqgrid
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: ntypat
  integer,intent(in) :: option
  integer,intent(in) :: paral_kgb
  real(dp),intent(in) :: eei
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: qprtrb(3)
  real(dp),intent(out) :: dyfrlo(3,3,natom)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(out) :: grtn(3,natom)
  real(dp),intent(out) :: lpsstr(6)
  integer,intent(in) :: nattyp(ntypat)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(in) :: vlspl(mqgrid,2,ntypat)
  real(dp),intent(in) :: vprtrb(2)
  real(dp),intent(inout) :: vpsp(nfft)
 end subroutine mklocl_recipspace
end interface

interface
 subroutine mklocl_wavelets(efield, grtn, mpi_enreg, natom, nfft,&  
  &  nspden, option, rprimd, vpsp, wvl_den, wvl_descr, xcart)
  use defs_basis
  use defs_abitypes
  use defs_wvltypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: option
  type(mpi_type),intent(inout) :: mpi_enreg
  type(wvl_denspot_type), intent(inout) :: wvl_den
  type(wvl_internal_type), intent(in) :: wvl_descr
  real(dp),intent(in) :: efield(3)
  real(dp),intent(inout) :: grtn(3,natom)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp), intent(inout) :: vpsp(nfft)
  real(dp),intent(inout) :: xcart(3,natom)
 end subroutine mklocl_wavelets
end interface

interface
 subroutine mkresi(cg,eig_k,gs_hamk,icg,ikpt,isppol,mcg,mpi_enreg,nband,prtvol,resid_k)
  use defs_basis
  use defs_abitypes
  use m_hamiltonian
  implicit none
  integer,intent(in) :: icg
  integer,intent(in) :: ikpt
  integer,intent(in) :: isppol
  integer,intent(in) :: mcg
  integer,intent(in) :: nband
  integer,intent(in) :: prtvol
  type(gs_hamiltonian_type),intent(inout) :: gs_hamk
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: cg(2,mcg)
  real(dp),intent(out) :: eig_k(nband)
  real(dp),intent(out) :: resid_k(nband)
 end subroutine mkresi
end interface

interface
 subroutine mkrho(cg,dtset,gprimd,irrzon,kg,mcg,mpi_enreg,npwarr,occ,paw_dmft,phnons,&  
  &  rhog,rhor,rprimd,tim_mkrho,ucvol,wvl_den,wvl_wfs,&  
  &  option) !optional
  use defs_basis
  use defs_abitypes
  use defs_wvltypes
  use m_paw_dmft
  implicit none
  integer,intent(in) :: mcg
  integer,intent(in),optional :: option
  integer,intent(in) :: tim_mkrho
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  type(paw_dmft_type), intent(in) :: paw_dmft
  real(dp),intent(in) :: ucvol
  type(wvl_denspot_type), intent(inout) :: wvl_den
  type(wvl_wf_type),intent(inout) :: wvl_wfs
  real(dp), intent(in) :: cg(2,mcg)
  real(dp), intent(in) :: gprimd(3,3)
  integer, intent(in) :: irrzon(dtset%nfft**(1-1/dtset%nsym),2, &
  &         (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
  integer, intent(in) :: npwarr(dtset%nkpt)
  real(dp), intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp), intent(in) :: phnons(2,(dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3))**(1-1/dtset%nsym), &
  &         (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  real(dp), intent(out) :: rhog(2,dtset%nfft)
  real(dp), intent(out) :: rhor(dtset%nfft,dtset%nspden)
  real(dp), intent(in) :: rprimd(3,3)
 end subroutine mkrho
end interface

interface
 subroutine mksubham(cg,ghc,gsc,gvnlc,iblock,icg,igsc,istwf_k,&  
  &  isubh,isubo,mcg,mgsc,nband_k,nbdblock,npw_k,&  
  &  nspinor,subham,subovl,subvnl,use_subovl,use_vnl,me_g0)
  use defs_basis
  implicit none
  integer,intent(in) :: iblock
  integer,intent(in) :: icg
  integer,intent(in) :: igsc
  integer,intent(in) :: istwf_k
  integer,intent(inout) :: isubh
  integer,intent(inout) :: isubo
  integer,intent(in) :: mcg
  integer,intent(in) :: me_g0
  integer,intent(in) :: mgsc
  integer,intent(in) :: nband_k
  integer,intent(in) :: nbdblock
  integer,intent(in) :: npw_k
  integer,intent(in) :: nspinor
  integer,intent(in) :: use_subovl
  integer,intent(in) :: use_vnl
  real(dp),intent(in) :: cg(2,mcg)
  real(dp),intent(inout) :: ghc(2,npw_k*nspinor)
  real(dp),intent(in) :: gsc(2,mgsc)
  real(dp),intent(inout) :: gvnlc(2,npw_k*nspinor)
  real(dp),intent(inout) :: subham(nband_k*(nband_k+1))
  real(dp),intent(inout) :: subovl(nband_k*(nband_k+1)*use_subovl)
  real(dp),intent(inout) :: subvnl(nband_k*(nband_k+1)*use_vnl)
 end subroutine mksubham
end interface

interface
 subroutine mlwfovlp(atindx1,cg,cprj,dtset,dtfil,eigen,gprimd,hdr,kg,&  
  &  mband,mcg,mcprj,mgfftc,mkmem,mpi_enreg,mpw,natom,&  
  &  nattyp,nfft,ngfft,nkpt,npwarr,nsppol,ntypat,occ,&  
  &  pawang,pawrad,pawtab,prtvol,psps,rprimd,ucvol,xred)
  use m_pawrad
  use m_pawang
  use defs_abitypes
  use m_pawcprj
  use m_pawtab
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: mcg
  integer,intent(in) :: mcprj
  integer,intent(in) :: mgfftc
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  integer,intent(in) :: ntypat
  integer,intent(in) :: prtvol
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  type(hdr_type),intent(in) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: ucvol
  integer :: ngfft(18)
  integer,intent(in) :: atindx1(natom)
  real(dp),intent(in) :: cg(2,mcg)
  type(pawcprj_type) :: cprj(natom,mcprj)
  real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
  real(dp),intent(in) :: gprimd(3,3)
  integer :: kg(3,mpw*mkmem)
  integer :: nattyp(ntypat)
  integer :: npwarr(nkpt)
  real(dp),intent(in) :: occ(mband*nkpt*nsppol)
  type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
  type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine mlwfovlp
end interface

interface
 subroutine mlwfovlp_proj(A_matrix,band_in,cg,cprj,dtset,gprimd,just_augmentation,kg,&  
  &  lproj,max_num_bands,mband,mkmem,mpi_enreg,mpw,mwan,natom,nattyp,&  
  &  nkpt,npwarr,nspinor,&  
  &  nsppol,ntypat,num_bands,nwan,pawtab,proj_l,proj_m,proj_radial,&  
  &  proj_site,proj_x,proj_z,proj_zona,psps,spin,ucvol)
  use defs_basis
  use defs_abitypes
  use m_pawcprj
  use defs_datatypes
  use m_pawtab
  implicit none
  integer,intent(in) :: lproj
  integer,intent(in) :: max_num_bands
  integer,intent(in) :: mband
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: mwan
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: ntypat
  integer,intent(in) :: spin
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  real(dp) :: ucvol
  complex(dpc),intent(out) :: A_matrix(max_num_bands,mwan,nkpt,nsppol)
  logical,intent(in) :: band_in(mband,nsppol)
  real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  type(pawcprj_type) :: cprj(natom,nspinor*mband*mkmem*nsppol)
  real(dp),intent(in) :: gprimd(3,3)
  logical,intent(in) :: just_augmentation(mwan,nsppol)
  integer,intent(in) :: kg(3,mpw*mkmem)
  integer :: nattyp(ntypat)
  integer,intent(in) :: npwarr(nkpt)
  integer,intent(in) :: num_bands(nsppol)
  integer,intent(in) :: nwan(nsppol)
  type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)
  integer,intent(in) :: proj_l(mband,nsppol)
  integer,intent(in) :: proj_m(mband,nsppol)
  integer,intent(inout) :: proj_radial(mband,nsppol)
  real(dp),intent(in) :: proj_site(3,mband,nsppol)
  real(dp),intent(in) :: proj_x(3,mband,nsppol)
  real(dp),intent(in) :: proj_z(3,mband,nsppol)
  real(dp),intent(in) :: proj_zona(mband,nsppol)
 end subroutine mlwfovlp_proj
end interface

interface
 subroutine mlwfovlp_projpaw(A_paw,band_in,cprj,just_augmentation,max_num_bands,mband,mkmem,&  
  &  mwan,natom,nband,nkpt,&  
  &  nspinor,nsppol,ntypat,nwan,pawrad,pawtab,&  
  &  proj_l,proj_m,proj_radial,proj_site,proj_x,proj_z,proj_zona,psps,&  
  &  rprimd,spin,typat,xred)
  use defs_basis
  use m_pawcprj
  use m_pawrad
  use defs_datatypes
  use m_pawtab
  implicit none
  integer,intent(in) :: max_num_bands
  integer,intent(in) :: mband
  integer,intent(in) :: mkmem
  integer,intent(in) :: mwan
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: ntypat
  integer,intent(in) :: spin
  type(pseudopotential_type),intent(in) :: psps
  complex(dpc),intent(out) :: A_paw(max_num_bands,mwan,nkpt,nsppol)
  logical,intent(in) :: band_in(mband,nsppol)
  type(pawcprj_type) :: cprj(natom,nspinor*mband*mkmem*nsppol)
  logical,intent(in) :: just_augmentation(mwan,nsppol)
  integer,intent(in) :: nband(nsppol*nkpt)
  integer,intent(in) :: nwan(nsppol)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  integer,intent(in) :: proj_l(mband,nsppol)
  integer,intent(in) :: proj_m(mband,nsppol)
  integer,intent(in) :: proj_radial(mband,nsppol)
  real(dp),intent(in) :: proj_site(3,mband,nsppol)
  real(dp),intent(in) :: proj_x(3,mband,nsppol)
  real(dp),intent(in) :: proj_z(3,mband,nsppol)
  real(dp),intent(in) :: proj_zona(mband,nsppol)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine mlwfovlp_projpaw
end interface

interface
 subroutine mlwfovlp_pw(cg,cm1,g1,iwav,kg,mband,mkmem,mpi_enreg,mpw,nfft,ngfft,nkpt,nntot,&  
  &  npwarr,nspinor,nsppol,ovikp,seed_name,spin)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: nfft
  integer,intent(in) :: nkpt
  integer,intent(in) :: nntot
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: spin
  type(mpi_type),intent(in) :: mpi_enreg
  character(len=fnlen) :: seed_name
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  real(dp),intent(out) :: cm1(2,mband,mband,nntot,nkpt,nsppol)
  integer,intent(in) :: g1(3,nkpt,nntot)
  integer,intent(in) :: iwav(mband,nkpt,nsppol)
  integer,intent(in) :: kg(3,mpw*mkmem)
  integer,intent(in) :: npwarr(nkpt)
  integer,intent(in) :: ovikp(nkpt,nntot)
 end subroutine mlwfovlp_pw
end interface

interface
 subroutine mlwfovlp_radial(alpha,lmax,lmax2,radial,rvalue,xx)
  use defs_basis
  implicit none
  integer,intent(in) :: lmax
  integer,intent(in) :: lmax2
  integer,intent(in) :: rvalue
  real(dp),intent(in) :: alpha
  real(dp),intent(in) :: xx
  real(dp),intent(out) :: radial(lmax2)
 end subroutine mlwfovlp_radial
end interface

interface
 subroutine mlwfovlp_seedname(fname_w90,filew90_win,filew90_wout,filew90_amn,&  
  &  filew90_ramn,filew90_mmn,filew90_eig,nsppol,seed_name)
  use defs_basis
  implicit none
  integer,intent(in) :: nsppol
  character(len=fnlen),intent(in) :: fname_w90
  character(len=fnlen),intent(out) :: filew90_amn(nsppol)
  character(len=fnlen),intent(out) :: filew90_eig(nsppol)
  character(len=fnlen),intent(out) :: filew90_mmn(nsppol)
  character(len=fnlen),intent(out) :: filew90_ramn(nsppol)
  character(len=fnlen),intent(out) :: filew90_win(nsppol)
  character(len=fnlen),intent(out) :: filew90_wout(nsppol)
  character(len=fnlen),intent(out) :: seed_name(nsppol)
 end subroutine mlwfovlp_seedname
end interface

interface
 subroutine mlwfovlp_setup(atom_symbols,band_in,dtset,filew90_win,gamma_only,&  
  &  g1,lwanniersetup,mband,natom,nband_inc,nkpt,&  
  &  nntot,num_bands,num_nnmax,nsppol,nwan,ovikp,&  
  &  proj_l,proj_m,proj_radial,proj_site,proj_x,proj_z,proj_zona,&  
  &  real_lattice,recip_lattice,rprimd,seed_name,spin,spinors,xcart,xred)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: lwanniersetup
  integer,intent(in) :: mband
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt
  integer,intent(out) :: nntot
  integer,intent(in) :: nsppol
  integer,intent(in) :: num_nnmax
  integer,intent(in) :: spin
  type(dataset_type),intent(in) :: dtset
  logical,intent(in) :: gamma_only
  logical,intent(in) :: spinors
  character(len=3),intent(out) :: atom_symbols(natom)
  logical,intent(out) :: band_in(mband,nsppol)
  character(len=fnlen),intent(in) :: filew90_win(nsppol)
  integer,intent(out) :: g1(3,nkpt,num_nnmax)
  integer,intent(out) :: nband_inc(nsppol)
  integer,intent(out) :: num_bands(nsppol)
  integer,intent(out) :: nwan(nsppol)
  integer,intent(out) :: ovikp(nkpt,num_nnmax)
  integer,intent(out) :: proj_l(mband,nsppol)
  integer,intent(out) :: proj_m(mband,nsppol)
  integer,intent(out) :: proj_radial(mband,nsppol)
  real(dp),intent(out) :: proj_site(3,mband,nsppol)
  real(dp),intent(out) :: proj_x(3,mband,nsppol)
  real(dp),intent(out) :: proj_z(3,mband,nsppol)
  real(dp),intent(out) :: proj_zona(mband,nsppol)
  real(dp),intent(in) :: real_lattice(3,3)
  real(dp),intent(in) :: recip_lattice(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  character(len=fnlen),intent(in) :: seed_name(nsppol)
  real(dp),intent(out) :: xcart(3,natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine mlwfovlp_setup
end interface

interface
 subroutine mlwfovlp_ylmfac(ylmc_fac,lmax,lmax2,mband,nwan,proj_l,proj_m,proj_x,proj_z)
  use defs_basis
  implicit none
  integer, intent(in) :: lmax
  integer, intent(in) :: lmax2
  integer, intent(in) :: mband
  integer, intent(in) :: nwan
  integer,intent(in) :: proj_l(mband)
  integer,intent(in) :: proj_m(mband)
  real(dp),intent(in) :: proj_x(3,mband)
  real(dp),intent(in) :: proj_z(3,mband)
  complex(dp),intent(out) :: ylmc_fac(lmax2,nwan)
 end subroutine mlwfovlp_ylmfac
end interface

interface
 subroutine mlwfovlp_ylmfar(ylmr_fac,lmax,lmax2,mband,nwan,proj_l,proj_m,proj_x,proj_z)
  use defs_basis
  implicit none
  integer, intent(in) :: lmax
  integer, intent(in) :: lmax2
  integer, intent(in) :: mband
  integer, intent(in) :: nwan
  integer,intent(in) :: proj_l(mband)
  integer,intent(in) :: proj_m(mband)
  real(dp),intent(in) :: proj_x(3,mband)
  real(dp),intent(in) :: proj_z(3,mband)
  real(dp),intent(out) :: ylmr_fac(lmax2,nwan)
 end subroutine mlwfovlp_ylmfar
end interface

interface
 subroutine moddiel(cplex,dielar,mpi_enreg,nfft,ngfft,nspden,optreal,optres,paral_kgb,qphon,rprimd,vresid,vrespc)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: optreal
  integer,intent(in) :: optres
  integer,intent(in) :: paral_kgb
  type(mpi_type),intent(in) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: dielar(7)
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: vresid(cplex*nfft,nspden)
  real(dp),intent(out) :: vrespc(cplex*nfft,nspden)
 end subroutine moddiel
end interface

interface
 subroutine msig(fcti,npti,xi,filnam_out_sig)
  use defs_basis
  implicit none
  integer,intent(in) :: npti
  character(len=fnlen),intent(in) :: filnam_out_sig
  real(dp),intent(in) :: fcti(npti)
  real(dp),intent(in) :: xi(npti)
 end subroutine msig
end interface

interface
 subroutine newkpt(ceksp2,cg,debug,ecut1,ecut2,ecut2_eff,eigen,exchn2n3d,fill,&  
  &  formeig,gmet1,gmet2,headform1,indkk,iout,ireadwf,&  
  &  istwfk1,istwfk2,kg2,kptns1,kptns2,mband2,mcg,mkmem1,mkmem2,&  
  &  mpi_enreg1,mpi_enreg2,mpw1,mpw2,my_nkpt2,nband1,nband2,&  
  &  ngfft1,ngfft2,nkpt1,nkpt2,npwarr1,npwarr2,nspinor1,nspinor2,&  
  &  nsppol1,nsppol2,nsym,occ,optorth,prtvol,randalg,restart,rprimd,&  
  &  sppoldbl,symrel,tnons,unkg2,wffinp,wffout)
  use defs_basis
  use defs_abitypes
  use m_wffile
  implicit none
  integer,intent(in) :: ceksp2
  integer,intent(in) :: debug
  integer,intent(in) :: exchn2n3d
  integer,intent(in) :: fill
  integer,intent(in) :: formeig
  integer,intent(in) :: headform1
  integer,intent(in) :: iout
  integer,intent(in) :: ireadwf
  integer,intent(in) :: mband2
  integer,intent(in) :: mcg
  integer,intent(in) :: mkmem1
  integer,intent(in) :: mkmem2
  integer,intent(in) :: mpw1
  integer,intent(in) :: mpw2
  integer,intent(in) :: my_nkpt2
  integer,intent(in) :: nkpt1
  integer,intent(in) :: nkpt2
  integer,intent(in) :: nspinor1
  integer,intent(in) :: nspinor2
  integer,intent(in) :: nsppol1
  integer,intent(in) :: nsppol2
  integer,intent(in) :: nsym
  integer,intent(in) :: optorth
  integer,intent(in) :: prtvol
  integer,intent(in) :: randalg
  integer,intent(in) :: restart
  integer,intent(in) :: sppoldbl
  integer,intent(in) :: unkg2
  real(dp),intent(in) :: ecut1
  real(dp),intent(in) :: ecut2
  real(dp),intent(in) :: ecut2_eff
  type(mpi_type),intent(inout) :: mpi_enreg1
  type(mpi_type),intent(inout) :: mpi_enreg2
  type(wffile_type),intent(inout) :: wffinp
  type(wffile_type),intent(inout) :: wffout
  integer,intent(in) :: ngfft1(18)
  integer,intent(in) :: ngfft2(18)
  real(dp),intent(inout) :: cg(2,mcg)
  real(dp),intent(inout) :: eigen(mband2*(2*mband2)**formeig*nkpt2*nsppol2)
  real(dp),intent(in) :: gmet1(3,3)
  real(dp),intent(in) :: gmet2(3,3)
  integer,intent(in) :: indkk(nkpt2*sppoldbl,6)
  integer,intent(in) :: istwfk1(nkpt1)
  integer,intent(in) :: istwfk2(nkpt2)
  integer,intent(in) :: kg2(3,mpw2*mkmem2)
  real(dp),intent(in) :: kptns1(3,nkpt1)
  real(dp),intent(in) :: kptns2(3,nkpt2)
  integer,intent(in) :: nband1(nkpt1*nsppol1)
  integer,intent(in) :: nband2(nkpt2*nsppol2)
  integer,intent(in) :: npwarr1(nkpt1)
  integer,intent(in) :: npwarr2(nkpt2)
  real(dp),intent(inout) :: occ(mband2*nkpt2*nsppol2)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
 end subroutine newkpt
end interface

interface
 subroutine nres2vres(dtset,gsqcut,izero,kxc,mpi_enreg,my_natom,nfft,ngfft,nhat,&  
  &  nkxc,nresid,n3xccc,optnc,optxc,pawang,pawfgrtab,pawrhoij,pawtab,&  
  &  rhor,rprimd,usepaw,vresid,xccc3d,xred)
  use m_pawang
  use m_pawrhoij
  use defs_abitypes
  use m_pawfgrtab
  use m_pawtab
  use defs_basis
  implicit none
  integer,intent(in) :: izero
  integer,intent(in) :: my_natom
  integer,intent(in) :: n3xccc
  integer,intent(in) :: nfft
  integer,intent(in) :: nkxc
  integer,intent(in) :: optnc
  integer,intent(in) :: optxc
  integer,intent(in) :: usepaw
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: kxc(nfft,nkxc)
  real(dp),intent(inout) :: nhat(nfft,dtset%nspden*usepaw)
  real(dp),intent(in) :: nresid(nfft,dtset%nspden)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(my_natom*usepaw)
  type(pawrhoij_type),intent(in) :: pawrhoij(my_natom*usepaw)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*usepaw)
  real(dp),intent(in) :: rhor(nfft,dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: vresid(nfft,dtset%nspden)
  real(dp),intent(in) :: xccc3d(n3xccc)
  real(dp),intent(in) :: xred(3,dtset%natom)
 end subroutine nres2vres
end interface

interface
 subroutine odamix(deltae,dtset,elast,energies,etotal,&  
  &  gprimd,gsqcut,kxc,mpi_enreg,my_natom,nfft,ngfft,nhat,&  
  &  nkxc,ntypat,nvresid,n3xccc,optres,paw_ij,&  
  &  paw_an,pawang,pawfgrtab,pawrad,pawrhoij,pawtab,&  
  &  red_ptot,psps,rhog,rhor,rprimd,strsxc,ucvol,usepaw,&  
  &  usexcnhat,vhartr,vpsp,vtrial,vxc,vxcavg,xccc3d,xred,&  
  &  taug,taur,vxctau,add_tfw) ! optional arguments
  use m_pawrad
  use m_paw_ij
  use m_pawang
  use m_pawrhoij
  use m_energies
  use defs_abitypes
  use m_pawfgrtab
  use m_paw_an
  use m_pawtab
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: my_natom
  integer,intent(in) :: n3xccc
  integer,intent(in) :: nfft
  integer,intent(in) :: nkxc
  integer,intent(in) :: ntypat
  integer,intent(in) :: optres
  integer,intent(in) :: usepaw
  integer,intent(in) :: usexcnhat
  logical,intent(in),optional :: add_tfw
  real(dp),intent(out) :: deltae
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(inout) :: elast
  type(energies_type),intent(inout) :: energies
  real(dp),intent(out) :: etotal
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: ucvol
  real(dp),intent(out) :: vxcavg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(inout) :: kxc(nfft,nkxc)
  real(dp),intent(inout) :: nhat(nfft,dtset%nspden*usepaw)
  real(dp),intent(inout) :: nvresid(nfft,dtset%nspden)
  type(paw_an_type),intent(inout) :: paw_an(my_natom)
  type(paw_ij_type),intent(inout) :: paw_ij(my_natom)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(my_natom)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawrhoij_type),intent(inout) :: pawrhoij(my_natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(in) :: red_ptot(3)
  real(dp),intent(inout) :: rhog(2,nfft)
  real(dp),intent(inout) :: rhor(nfft,dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: strsxc(6)
  real(dp),intent(in),optional :: taug(2,nfft*dtset%usekden)
  real(dp),intent(in),optional :: taur(nfft,dtset%nspden*dtset%usekden)
  real(dp),intent(inout) :: vhartr(nfft)
  real(dp),intent(in) :: vpsp(nfft)
  real(dp),intent(inout) :: vtrial(nfft,dtset%nspden)
  real(dp),intent(inout) :: vxc(nfft,dtset%nspden)
  real(dp),intent(inout),optional :: vxctau(nfft,dtset%nspden,4*dtset%usekden)
  real(dp),intent(inout) :: xccc3d(n3xccc)
  real(dp),intent(in) :: xred(3,dtset%natom)
 end subroutine odamix
end interface

interface
 subroutine optics_vloc(cg,dtfil,dtset,eigen0,gprimd,hdr,kg,mband,mcg,mkmem,mpi_enreg,mpw,&  
  &  nkpt,npwarr,nsppol)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: mcg
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  type(hdr_type),intent(inout) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(inout) :: cg(2,mcg)
  real(dp),intent(in) :: eigen0(mband*nkpt*nsppol)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: kg(3,mpw*mkmem)
  integer,intent(in) :: npwarr(nkpt)
 end subroutine optics_vloc
end interface

interface
 subroutine partial_dos_fractions(cg,dos_fractions,dos_fractions_m,&  
  &  dtset,hdr,mbesslang,mcg,mpi_enreg,m_dos_flag,ndosfraction,partial_dos)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: m_dos_flag
  integer,intent(in) :: mbesslang
  integer,intent(in) :: mcg
  integer,intent(in) :: ndosfraction
  integer,intent(in) :: partial_dos
  type(dataset_type),intent(in) :: dtset
  type(hdr_type),intent(inout) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(inout) :: cg(2,mcg)
  real(dp),intent(out) :: dos_fractions(dtset%nkpt,dtset%mband,dtset%nsppol,ndosfraction)
  real(dp),intent(out) :: dos_fractions_m(dtset%nkpt,dtset%mband, &
  &         dtset%nsppol,ndosfraction*mbesslang*m_dos_flag)
 end subroutine partial_dos_fractions
end interface

interface
 subroutine posdoppler(cg,cprj,Crystal,dimcprj,dtfil,dtset,electronpositron,&  
  &  filpsp,kg,mcg,mcprj,mpi_enreg,my_natom,&  
  &  n3xccc,nfft,ngfft,nhat,npwarr,occ,pawang,pawrad,&  
  &  pawrhoij,pawtab,rhor,xccc3d)
  use m_pawtab
  use m_pawrad
  use m_pawang
  use m_pawrhoij
  use defs_abitypes
  use m_pawcprj
  use m_crystal
  use defs_basis
  use m_electronpositron
  implicit none
  integer,intent(in) :: mcg
  integer,intent(in) :: mcprj
  integer,intent(in) :: my_natom
  integer,intent(in) :: n3xccc
  integer,intent(in) :: nfft
  type(crystal_t) :: Crystal
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  type(electronpositron_type),pointer :: electronpositron
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  integer,intent(in) :: ngfft(18)
  real(dp),intent(inout),target :: cg(2,mcg)
  type(pawcprj_type),target :: cprj(dtset%natom,mcprj)
  integer,intent(in) :: dimcprj(dtset%natom)
  character(len=fnlen),intent(in) :: filpsp(dtset%ntypat)
  integer,intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
  real(dp),intent(in) :: nhat(nfft,dtset%nspden*dtset%usepaw)
  integer,intent(in) :: npwarr(dtset%nkpt)
  real(dp),intent(in),target :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat*dtset%usepaw)
  type(pawrhoij_type),intent(in),target :: pawrhoij(my_natom*dtset%usepaw)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*dtset%usepaw)
  real(dp),intent(in),target :: rhor(nfft,dtset%nspden)
  real(dp),intent(in) :: xccc3d(n3xccc)
 end subroutine posdoppler
end interface

interface
 subroutine poslifetime(dtset,electronpositron,gprimd,my_natom,mpi_enreg,n3xccc,nfft,ngfft,nhat,&  
  &  option,pawang,pawrad,pawrhoij,pawtab,rate,rate_paw,rhor,ucvol,xccc3d,&  
  &  rhor_dop_el,pawrhoij_dop_el,pawrhoij_ep) ! optional arguments
  use m_pawrad
  use m_pawang
  use m_pawrhoij
  use defs_abitypes
  use m_pawtab
  use defs_basis
  use m_electronpositron
  implicit none
  integer,intent(in) :: my_natom
  integer,intent(in) :: n3xccc
  integer,intent(in) :: nfft
  integer,intent(in) :: option
  type(dataset_type), intent(in) :: dtset
  type(electronpositron_type),pointer :: electronpositron
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type), intent(in) :: pawang
  real(dp),intent(out) :: rate
  real(dp),intent(out) :: rate_paw
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: nhat(nfft,dtset%nspden*dtset%usepaw)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat*dtset%usepaw)
  type(pawrhoij_type),intent(in) :: pawrhoij(my_natom*dtset%usepaw)
  type(pawrhoij_type),optional,intent(in) :: pawrhoij_dop_el(my_natom*dtset%usepaw)
  type(pawrhoij_type),optional,target,intent(in) :: pawrhoij_ep(my_natom*dtset%usepaw)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*dtset%usepaw)
  real(dp),intent(in),target :: rhor(nfft,dtset%nspden)
  real(dp),optional,intent(in) :: rhor_dop_el(nfft)
  real(dp),intent(in) :: xccc3d(n3xccc)
 end subroutine poslifetime
end interface

interface
 subroutine posratecore(dtset,electronpositron,iatom,my_natom,mesh_sizej,mpi_enreg,&  
  &  option,pawang,pawrad,pawrhoij,pawrhoij_ep,&  
  &  pawtab,rate,rhocorej)
  use m_pawrad
  use m_pawang
  use m_pawrhoij
  use defs_abitypes
  use defs_basis
  use m_pawtab
  use m_electronpositron
  implicit none
  integer,intent(in) :: iatom
  integer,intent(in) :: mesh_sizej
  integer,intent(in) :: my_natom
  integer,intent(in) :: option
  type(dataset_type), intent(in) :: dtset
  type(electronpositron_type),pointer :: electronpositron
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type), intent(in) :: pawang
  real(dp),intent(out) :: rate
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat*dtset%usepaw)
  type(pawrhoij_type),intent(in) :: pawrhoij(my_natom*dtset%usepaw)
  type(pawrhoij_type),intent(in),target :: pawrhoij_ep(my_natom*dtset%usepaw)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*dtset%usepaw)
  real(dp),intent(in) :: rhocorej(mesh_sizej)
 end subroutine posratecore
end interface

interface
 subroutine prtefield(dtset,dtefield,iunit,rprimd)
  use m_efield
  use defs_abitypes
  use defs_basis
  implicit none
  integer :: iunit
  type(efield_type),intent(in) :: dtefield
  type(dataset_type),intent(inout) :: dtset
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine prtefield
end interface

interface
 subroutine prteigrs(eigen,enunit,fermie,fname_eig,iout,iscf,kptns,kptopt,mband,nband,&  
  &  nkpt,nnsclo_now,nsppol,occ,occopt,option,prteig,prtvol,resid,tolwfr,vxcavg,wtk)
  use defs_basis
  implicit none
  integer,intent(in) :: enunit
  integer,intent(in) :: iout
  integer,intent(in) :: iscf
  integer,intent(in) :: kptopt
  integer,intent(in) :: mband
  integer,intent(in) :: nkpt
  integer,intent(in) :: nnsclo_now
  integer,intent(in) :: nsppol
  integer,intent(in) :: occopt
  integer,intent(in) :: option
  integer,intent(in) :: prteig
  integer,intent(in) :: prtvol
  real(dp),intent(in) :: fermie
  character(len=*),intent(in) :: fname_eig
  real(dp),intent(in) :: tolwfr
  real(dp),intent(in) :: vxcavg
  real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
  real(dp),intent(in) :: kptns(3,nkpt)
  integer,intent(in) :: nband(nkpt*nsppol)
  real(dp),intent(in) :: occ(mband*nkpt*nsppol)
  real(dp),intent(in) :: resid(mband*nkpt*nsppol)
  real(dp),intent(in) :: wtk(nkpt)
 end subroutine prteigrs
end interface

interface
 subroutine prtene(dtset,energies,iout,usepaw)
  use m_energies
  use defs_abitypes
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: usepaw
  type(dataset_type),intent(in) :: dtset
  type(energies_type),intent(in) :: energies
 end subroutine prtene
end interface

interface
 subroutine prtimg(dynimage,imagealgo_str,imgmov,iout,mpi_enreg,nimage,nimage_tot,&  
  &  prt_all_images,prtvolimg,resimg)
  use m_results_img
  use defs_abitypes
  implicit none
  integer,intent(in) :: imgmov
  integer,intent(in) :: iout
  integer,intent(in) :: nimage
  integer,intent(in) :: nimage_tot
  integer,intent(in) :: prtvolimg
  character(len=60),intent(in) :: imagealgo_str
  type(mpi_type),intent(inout) :: mpi_enreg
  logical,intent(in) :: prt_all_images
  integer,intent(in) :: dynimage(nimage_tot)
  type(results_img_type),target,intent(inout) :: resimg(nimage)
 end subroutine prtimg
end interface

interface
 subroutine prtrhomxmn(iout,mpi_enreg,nfft,ngfft,nspden,option,rhor,optrhor,ucvol)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: option
  integer,intent(in),optional :: optrhor
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in),optional :: ucvol
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: rhor(nfft,nspden)
 end subroutine prtrhomxmn
end interface

interface
 subroutine prtxf(fred,iatfix,iout,iwfrc,natom,rprimd,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: iwfrc
  integer,intent(in) :: natom
  real(dp),intent(in) :: fred(3,natom)
  integer,intent(in) :: iatfix(3,natom)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine prtxf
end interface

interface
 subroutine prtxvf(fcart,fred,iatfix,iout,natom,prtvel,vel,xcart,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: natom
  integer,intent(in) :: prtvel
  real(dp),intent(in) :: fcart(3,natom)
  real(dp),intent(in) :: fred(3,natom)
  integer,intent(in) :: iatfix(3,natom)
  real(dp),intent(in) :: vel(3,natom)
  real(dp),intent(in) :: xcart(3,natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine prtxvf
end interface

interface
 subroutine rhotov(dtset,energies,gprimd,gsqcut,istep,kxc,mpi_enreg,nfft,ngfft,&  
  &  nhat,nhatgr,nhatgrdim,nkxc,vresidnew,n3xccc,optene,optres,optxc,&  
  &  rhog,rhor,rprimd,strsxc,ucvol,usepaw,usexcnhat,&  
  &  vhartr,vnew_mean,vpsp,vres_mean,vres2,vtrial,vxcavg,vxc,wvl,xccc3d,xred,&  
  &  electronpositron,taug,taur,vxctau,add_tfw) ! optional arguments
  use defs_basis
  use m_energies
  use defs_abitypes
  use m_electronpositron
  use defs_wvltypes
  implicit none
  integer,intent(in) :: istep
  integer,intent(in) :: n3xccc
  integer,intent(in) :: nfft
  integer,intent(in) :: nhatgrdim
  integer,intent(in) :: nkxc
  integer,intent(in) :: optene
  integer,intent(in) :: optres
  integer,intent(in) :: optxc
  integer,intent(in) :: usepaw
  integer,intent(in) :: usexcnhat
  logical,intent(in),optional :: add_tfw
  type(dataset_type),intent(inout) :: dtset
  type(electronpositron_type),pointer,optional :: electronpositron
  type(energies_type),intent(inout) :: energies
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  real(dp),intent(out) :: vres2
  real(dp),intent(out) :: vxcavg
  type(wvl_data), intent(inout) :: wvl
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(out) :: kxc(nfft,nkxc)
  real(dp),intent(in) :: nhat(nfft,dtset%nspden*usepaw)
  real(dp),intent(in) :: nhatgr(nfft,dtset%nspden,3*nhatgrdim)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(inout) :: rhor(nfft,dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: strsxc(6)
  real(dp),intent(in),optional :: taug(2,nfft*dtset%usekden)
  real(dp),intent(in),optional :: taur(nfft,dtset%nspden*dtset%usekden)
  real(dp),intent(inout) :: vhartr(nfft)
  real(dp),intent(out) :: vnew_mean(dtset%nspden)
  real(dp),intent(inout) :: vpsp(nfft)
  real(dp),intent(out) :: vres_mean(dtset%nspden)
  real(dp),intent(out) :: vresidnew(nfft,dtset%nspden)
  real(dp),intent(inout) :: vtrial(nfft,dtset%nspden)
  real(dp),intent(inout) :: vxc(nfft,dtset%nspden)
  real(dp),intent(out),optional :: vxctau(nfft,dtset%nspden,4*dtset%usekden)
  real(dp),intent(inout) :: xccc3d(n3xccc)
  real(dp),intent(inout) :: xred(3,dtset%natom)
 end subroutine rhotov
end interface

interface
 subroutine scprqt(choice,cpus,deltae,diffor,dtset,&  
  &  eigen,etotal,favg,fcart,fermie,fname_eig,filnam1,initGS,&  
  &  iscf,istep,kpt,maxfor,moved_atm_inside,mpi_enreg,&  
  &  nband,nkpt,nstep,occ,optres,&  
  &  prtfor,prtxml,quit,res2,resid,residm,response,tollist,usepaw,&  
  &  vxcavg,wtk,xred,conv_retcode,&  
  &  electronpositron) ! optional argument)
  use defs_basis
  use defs_abitypes
  use m_electronpositron
  implicit none
  integer,intent(in) :: choice
  integer,intent(out) :: conv_retcode
  integer,intent(in) :: initGS
  integer,intent(in) :: iscf
  integer,intent(in) :: istep
  integer,intent(in) :: moved_atm_inside
  integer,intent(in) :: nkpt
  integer,intent(in) :: nstep
  integer,intent(in) :: optres
  integer,intent(in) :: prtfor
  integer,intent(in) :: prtxml
  integer,intent(out) :: quit
  integer,intent(in) :: response
  integer,intent(in) :: usepaw
  real(dp),intent(in) :: cpus
  real(dp),intent(in) :: deltae
  real(dp),intent(in) :: diffor
  type(dataset_type),intent(in) :: dtset
  type(electronpositron_type),pointer,optional :: electronpositron
  real(dp),intent(in) :: etotal
  real(dp),intent(in) :: fermie
  character(len=fnlen),intent(in) :: filnam1
  character(len=fnlen),intent(in) :: fname_eig
  real(dp),intent(in) :: maxfor
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: res2
  real(dp),intent(in) :: residm
  real(dp),intent(in) :: vxcavg
  real(dp),intent(in) :: eigen(dtset%mband*nkpt*dtset%nsppol)
  real(dp),intent(in) :: favg(3)
  real(dp),intent(in) :: fcart(3,dtset%natom)
  real(dp),intent(in) :: kpt(3,nkpt)
  integer,intent(in) :: nband(nkpt*dtset%nsppol)
  real(dp),intent(in) :: occ(dtset%mband*nkpt*dtset%nsppol)
  real(dp),intent(in) :: resid(dtset%mband*nkpt*dtset%nsppol)
  real(dp),intent(in) :: tollist(12)
  real(dp),intent(in) :: wtk(nkpt)
  real(dp),intent(in) :: xred(3,dtset%natom)
 end subroutine scprqt
end interface

interface
 subroutine setup1(acell,amass,amu,bantot,dtset,ecut_eff,ecutc_eff,gmet,&  
  &  gprimd,gsqcut_eff,gsqcutc_eff,natom,ngfft,ngfftc,nkpt,nsppol,&  
  &  response,rmet,rprim,rprimd,ucvol,usepaw)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(out) :: bantot
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  integer,intent(in) :: response
  integer,intent(in) :: usepaw
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: ecut_eff
  real(dp),intent(in) :: ecutc_eff
  real(dp),intent(out) :: gsqcut_eff
  real(dp),intent(out) :: gsqcutc_eff
  real(dp),intent(out) :: ucvol
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: ngfftc(18)
  real(dp),intent(in) :: acell(3)
  real(dp),intent(out) :: amass(natom)
  real(dp),intent(in) :: amu(dtset%ntypat)
  real(dp),intent(out) :: gmet(3,3)
  real(dp),intent(out) :: gprimd(3,3)
  real(dp),intent(out) :: rmet(3,3)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(out) :: rprimd(3,3)
 end subroutine setup1
end interface

interface
 subroutine setup_positron(atindx,atindx1,cg,cprj,dtefield,dtfil,dtset,ecore,eigen,etotal,electronpositron,&  
  &  energies,forces_needed,fred,gmet,gprimd,grewtn,grvdw,gsqcut,hdr,ifirst_gs,indsym,istep,istep_mix,kg,&  
  &  kxc,maxfor,mcg,mcprj,mgfft,mpi_enreg,my_natom,n3xccc,nattyp,nfft,ngfft,ngrvdw,nhat,nkxc,npwarr,nvresid,occ,optres,&  
  &  paw_ij,pawang,pawfgr,pawfgrtab,pawrhoij,pawtab,ph1d,ph1dc,psps,rhog,rhor,&  
  &  rprimd,stress_needed,strsxc,symrec,ucvol,usecprj,vhartr,vpsp,vxc,&  
  &  xccc3d,xred,ylm,ylmgr)
  use m_pawtab
  use m_paw_ij
  use m_pawang
  use m_pawrhoij
  use m_energies
  use defs_abitypes
  use m_pawcprj
  use m_pawfgrtab
  use m_pawfgr
  use defs_basis
  use m_efield
  use defs_datatypes
  use m_electronpositron
  implicit none
  integer,intent(in) :: forces_needed
  integer,intent(in) :: ifirst_gs
  integer,intent(in) :: istep
  integer,intent(inout) :: istep_mix
  integer,intent(in) :: mcg
  integer,intent(in) :: mcprj
  integer,intent(in) :: mgfft
  integer,intent(in) :: my_natom
  integer,intent(in) :: n3xccc
  integer,intent(in) :: nfft
  integer,intent(in) :: ngrvdw
  integer,intent(in) :: nkxc
  integer,intent(in) :: optres
  integer,intent(in) :: stress_needed
  integer,intent(in) :: usecprj
  type(efield_type),intent(in) :: dtefield
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: ecore
  type(electronpositron_type),pointer :: electronpositron
  type(energies_type),intent(inout) :: energies
  real(dp),intent(in) :: etotal
  real(dp),intent(in) :: gsqcut
  type(hdr_type),intent(inout) :: hdr
  real(dp),intent(in) :: maxfor
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pawfgr_type),intent(in) :: pawfgr
  type(pseudopotential_type), intent(in) :: psps
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(in) :: atindx1(dtset%natom)
  real(dp),intent(inout) :: cg(2,mcg)
  type(pawcprj_type),intent(inout) :: cprj(dtset%natom,mcprj*usecprj)
  real(dp),intent(inout) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp),intent(inout) :: fred(3,dtset%natom)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: grewtn(3,dtset%natom)
  real(dp),intent(in) :: grvdw(3,ngrvdw)
  integer,intent(in) :: indsym(4,dtset%nsym,dtset%natom)
  integer,intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
  real(dp),intent(in) :: kxc(nfft,nkxc)
  integer,intent(in) :: nattyp(dtset%natom)
  real(dp),intent(inout) :: nhat(nfft,dtset%nspden*dtset%usepaw)
  integer,intent(in) :: npwarr(dtset%nkpt)
  real(dp),intent(inout) :: nvresid(nfft,dtset%nspden)
  real(dp),intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  type(paw_ij_type),intent(in) :: paw_ij(my_natom*dtset%usepaw)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(my_natom*dtset%usepaw)
  type(pawrhoij_type),intent(inout) :: pawrhoij(my_natom*dtset%usepaw)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*dtset%usepaw)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*dtset%natom)
  real(dp),intent(in) :: ph1dc(2,(3*(2*dtset%mgfft+1)*dtset%natom)*dtset%usepaw)
  real(dp),intent(inout) :: rhog(2,nfft)
  real(dp),intent(inout) :: rhor(nfft,dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: strsxc(6)
  integer,intent(in) :: symrec(3,3,dtset%nsym)
  real(dp),intent(in) :: vhartr(nfft)
  real(dp),intent(in) :: vpsp(nfft)
  real(dp),intent(in) :: vxc(nfft,dtset%nspden)
  real(dp),intent(inout) :: xccc3d(n3xccc)
  real(dp),intent(inout) :: xred(3,dtset%natom)
  real(dp),intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
  real(dp),intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 end subroutine setup_positron
end interface

interface
 subroutine setvtr(atindx1,dtset,energies,gmet,gprimd,grewtn,grvdw,gsqcut,&  
  &  istep,kxc,mgfft,moved_atm_inside,moved_rhor,mpi_enreg,&  
  &  nattyp,nfft,ngfft,ngrvdw,nhat,nhatgr,nhatgrdim,nkxc,ntypat,n1xccc,n3xccc,&  
  &  optene,pawrad,pawtab,ph1d,psps,rhog,rhor,rmet,rprimd,strsxc,&  
  &  ucvol,usexcnhat,vhartr,vpsp,vtrial,vxc,vxcavg,wvl,xccc3d,xred,&  
  &  electronpositron,taug,taur,vxctau,add_tfw) ! optionals arguments
  use m_pawrad
  use m_energies
  use defs_abitypes
  use m_electronpositron
  use m_pawtab
  use defs_basis
  use defs_datatypes
  use defs_wvltypes
  implicit none
  integer,intent(in) :: istep
  integer,intent(in) :: mgfft
  integer,intent(inout) :: moved_atm_inside
  integer,intent(inout) :: moved_rhor
  integer,intent(in) :: n1xccc
  integer,intent(in) :: n3xccc
  integer,intent(in) :: nfft
  integer,intent(in) :: ngrvdw
  integer,intent(in) :: nhatgrdim
  integer,intent(in) :: nkxc
  integer,intent(in) :: ntypat
  integer,intent(in) :: optene
  integer,intent(in) :: usexcnhat
  logical,intent(in),optional :: add_tfw
  type(dataset_type),intent(inout) :: dtset
  type(electronpositron_type),pointer,optional :: electronpositron
  type(energies_type),intent(inout) :: energies
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: ucvol
  real(dp),intent(out) :: vxcavg
  type(wvl_data), intent(inout) :: wvl
  integer, intent(in) :: ngfft(18)
  integer, intent(in) :: atindx1(dtset%natom)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(out) :: grewtn(3,dtset%natom)
  real(dp),intent(out) :: grvdw(3,ngrvdw)
  real(dp),intent(out) :: kxc(nfft,nkxc)
  integer, intent(in) :: nattyp(ntypat)
  real(dp),intent(in) :: nhat(nfft,dtset%nspden*psps%usepaw)
  real(dp),intent(in) :: nhatgr(:,:,:)
  type(pawrad_type),intent(in) :: pawrad(ntypat*dtset%usepaw)
  type(pawtab_type),intent(in) :: pawtab(ntypat*dtset%usepaw)
  real(dp),intent(inout) :: ph1d(2,3*(2*mgfft+1)*dtset%natom)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(inout) :: rhor(nfft,dtset%nspden)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: strsxc(6)
  real(dp),intent(in),optional :: taug(2,nfft*dtset%usekden)
  real(dp),intent(inout),optional :: taur(nfft,dtset%nspden*dtset%usekden)
  real(dp),intent(inout) :: vhartr(nfft)
  real(dp),intent(inout) :: vpsp(nfft)
  real(dp),intent(inout) :: vtrial(nfft,dtset%nspden)
  real(dp),intent(inout) :: vxc(nfft,dtset%nspden)
  real(dp),intent(out),optional :: vxctau(nfft,dtset%nspden,4*dtset%usekden)
  real(dp),intent(inout) :: xccc3d(n3xccc)
  real(dp),intent(in) :: xred(3,dtset%natom)
 end subroutine setvtr
end interface

interface
 subroutine spin_current(cg,dtfil,dtset,gprimd,hdr,kg,mcg,mpi_enreg,psps)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: mcg
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  type(hdr_type),intent(inout) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: cg(2,mcg)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
 end subroutine spin_current
end interface

interface
 subroutine stress(atindx1,berryopt,dtefield,eei,efield,ehart,eii,gsqcut,ixc,kinstr,&  
  &  mgfft,mpi_enreg,mqgrid,n1xccc,n3xccc,natom,nattyp,&  
  &  nfft,ngfft,nlstr,nspden,nsym,ntypat,paral_kgb,psps,pawtab,ph1d,&  
  &  prtvol,qgrid,red_efieldbar,rhog,rprimd,strten,strsxc,symrec,&  
  &  typat,usepaw,vdw_tol,vdw_tol_3bt,vdw_xc,&  
  &  vlspl,vxc,xccc1d,xccc3d,xcccrc,xred,zion,znucl,qvpotzero,&  
  &  electronpositron) ! optional argument
  use defs_abitypes
  use m_pawtab
  use defs_basis
  use m_efield
  use defs_datatypes
  use m_electronpositron
  implicit none
  integer,intent(in) :: berryopt
  integer,intent(in) :: ixc
  integer,intent(in) :: mgfft
  integer,intent(in) :: mqgrid
  integer,intent(in) :: n1xccc
  integer,intent(in) :: n3xccc
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: prtvol
  integer,intent(in) :: usepaw
  integer,intent(in) :: vdw_xc
  type(efield_type),intent(in) :: dtefield
  real(dp),intent(in) :: eei
  real(dp),intent(in) :: ehart
  real(dp),intent(in) :: eii
  type(electronpositron_type),pointer,optional :: electronpositron
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: qvpotzero
  real(dp),intent(in) :: vdw_tol
  real(dp),intent(in) :: vdw_tol_3bt
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: atindx1(natom)
  real(dp),intent(in) :: efield(3)
  real(dp),intent(in) :: kinstr(6)
  integer,intent(in) :: nattyp(ntypat)
  real(dp),intent(in) :: nlstr(6)
  type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(in) :: red_efieldbar(3)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: strsxc(6)
  real(dp),intent(out) :: strten(6)
  integer,intent(in) :: symrec(3,3,nsym)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: vlspl(mqgrid,2,ntypat)
  real(dp),intent(in) :: vxc(nfft,nspden)
  real(dp),intent(in) :: xccc1d(n1xccc*(1-usepaw),6,ntypat)
  real(dp),intent(inout) :: xccc3d(n3xccc)
  real(dp),intent(in) :: xcccrc(ntypat)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: zion(ntypat)
  real(dp),intent(in) :: znucl(ntypat)
 end subroutine stress
end interface

interface
 subroutine strhar(ehart,gprimd,gsqcut,harstr,mpi_enreg,nfft,ngfft,rhog,ucvol,&  
  &  rhog2) ! optional argument
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: nfft
  real(dp),intent(in) :: ehart
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(in) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(out) :: harstr(6)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(in),optional :: rhog2(2,nfft)
 end subroutine strhar
end interface

interface
 subroutine sygrad(fred,natom,dedt,nsym,symrec,indsym)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nsym
  real(dp),intent(in) :: dedt(3,natom)
  real(dp),intent(out) :: fred(3,natom)
  integer,intent(in) :: indsym(4,nsym,natom)
  integer,intent(in) :: symrec(3,3,nsym)
 end subroutine sygrad
end interface

interface
 subroutine symrhg(cplex,gprimd,irrzon,mpi_enreg,nfft,nfftot,ngfft,nspden,nsppol,nsym,paral_kgb,&  
  &  phnons,rhog,rhor,rprimd,symafm,symrel)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: nfft
  integer,intent(in) :: nfftot
  integer,intent(in) :: nspden
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym
  integer,intent(in) :: paral_kgb
  type(mpi_type),intent(in) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: irrzon(nfftot**(1-1/nsym),2,(nspden/nsppol)-3*(nspden/4))
  real(dp),intent(in) :: phnons(2,nfftot**(1-1/nsym),(nspden/nsppol)-3*(nspden/4))
  real(dp),intent(out) :: rhog(2,nfft)
  real(dp),intent(inout) :: rhor(cplex*nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symafm(nsym)
  integer,intent(in) :: symrel(3,3,nsym)
 end subroutine symrhg
end interface

interface
 subroutine testsusmat(compute,dielop,dielstrt,dtset,istep)
  use defs_abitypes
  implicit none
  integer,intent(in) :: dielop
  integer,intent(in) :: dielstrt
  integer,intent(in) :: istep
  logical,intent(out) :: compute
  type(dataset_type),intent(in) :: dtset
 end subroutine testsusmat
end interface

interface
 subroutine uderiv(bdberry,cg,gprimd,hdr,istwfk,kberry,kg,kpt_,kptopt,kptrlatt,&  
  &  mband,mcg,mkmem,mpi_enreg,mpw,natom,nband,nberry,npwarr,nspinor,nsppol,nkpt_,&  
  &  unddk,fnameabo_1wf)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: kptopt
  integer,intent(in) :: mband
  integer,intent(in) :: mcg
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: natom
  integer,intent(in) :: nberry
  integer,intent(in) :: nkpt_
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: unddk
  character(len=fnlen),intent(in) :: fnameabo_1wf
  type(hdr_type),intent(inout) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: bdberry(4)
  integer,intent(in) :: kberry(3,20)
  integer,intent(in) :: kptrlatt(3,3)
  real(dp),intent(in) :: cg(2,mcg)
  real(dp),intent(in) :: gprimd(1:3,1:3)
  integer,intent(in) :: istwfk(nkpt_)
  integer,intent(in) :: kg(3,mpw*mkmem)
  real(dp),intent(in) :: kpt_(3,nkpt_)
  integer,intent(in) :: nband(nkpt_*nsppol)
  integer,intent(in) :: npwarr(nkpt_)
 end subroutine uderiv
end interface

interface
 subroutine update_e_field_vars(atindx,atindx1,cg,dimcprj,dtefield,dtfil,dtset,&  
  &  efield_old_cart,gmet,gprimd,hdr,idir,kg,mcg,&  
  &  mkmem,mpi_enreg,mpw,my_natom,natom,nattyp,ngfft,nkpt,npwarr,ntypat,&  
  &  pawrhoij,pawtab,pel_cg,pelev,pion,psps,ptot,ptot_cart,pwind,&  
  &  pwind_alloc,pwnsfac,red_efield2,red_efield2_old,red_ptot,rmet,rprimd,&  
  &  scfcv_level,scfcv_quit,scfcv_step,ucvol,unit_out,&  
  &  usepaw,xred,ylm,ylmgr)
  use m_pawrhoij
  use defs_abitypes
  use m_pawtab
  use defs_basis
  use m_efield
  use defs_datatypes
  implicit none
  integer, intent(in) :: idir
  integer, intent(in) :: mcg
  integer, intent(in) :: mkmem
  integer, intent(in) :: mpw
  integer, intent(in) :: my_natom
  integer, intent(in) :: natom
  integer, intent(in) :: nkpt
  integer, intent(in) :: ntypat
  integer, intent(in) :: pwind_alloc
  integer, intent(in) :: scfcv_level
  integer, intent(in) :: scfcv_quit
  integer, intent(in) :: scfcv_step
  integer, intent(in) :: unit_out
  integer, intent(in) :: usepaw
  type(efield_type), intent(inout) :: dtefield
  type(datafiles_type), intent(in) :: dtfil
  type(dataset_type), intent(inout) :: dtset
  type(hdr_type), intent(inout) :: hdr
  type(mpi_type), intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  real(dp), intent(in) :: ucvol
  integer, intent(in) :: ngfft(18)
  integer, intent(in) :: atindx(natom)
  integer, intent(in) :: atindx1(natom)
  real(dp), intent(in) :: cg(2,mcg)
  integer, intent(in) :: dimcprj(usepaw*natom)
  real(dp), intent(inout) :: efield_old_cart(3)
  real(dp), intent(in) :: gmet(3,3)
  real(dp), intent(in) :: gprimd(3,3)
  integer, intent(in) :: kg(3,mpw*mkmem)
  integer, intent(in) :: nattyp(ntypat)
  integer, intent(in) :: npwarr(nkpt)
  type(pawrhoij_type), intent(in) :: pawrhoij(my_natom*usepaw)
  type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)
  real(dp), intent(out) :: pel_cg(3)
  real(dp), intent(out) :: pelev(3)
  real(dp), intent(out) :: pion(3)
  real(dp), intent(out) :: ptot(3)
  real(dp), intent(inout) :: ptot_cart(3)
  integer, intent(in) :: pwind(pwind_alloc,2,3)
  real(dp), intent(in) :: pwnsfac(2,pwind_alloc)
  real(dp), intent(inout) :: red_efield2(3)
  real(dp), intent(inout) :: red_efield2_old(3)
  real(dp), intent(out) :: red_ptot(3)
  real(dp), intent(in) :: rmet(3,3)
  real(dp), intent(in) :: rprimd(3,3)
  real(dp), intent(inout) :: xred(3,natom)
  real(dp), intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
  real(dp), intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 end subroutine update_e_field_vars
end interface

interface
 subroutine vdw_dftd2(e_vdw_dftd2,ixc,natom,ntypat,prtvol,typat,rprimd,vdw_tol,xred,znucl,&  
  &  dyn_vdw_dftd2,elt_vdw_dftd2,fred_vdw_dftd2,str_vdw_dftd2,qphon) ! Optionals
  use defs_basis
  implicit none
  integer,intent(in) :: ixc
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  integer,intent(in) :: prtvol
  real(dp),intent(out) :: e_vdw_dftd2
  real(dp),intent(in) :: vdw_tol
  real(dp),intent(out),optional :: dyn_vdw_dftd2(2,3,natom,3,natom)
  real(dp),intent(out),optional :: elt_vdw_dftd2(6+3*natom,6)
  real(dp),intent(out),optional :: fred_vdw_dftd2(3,natom)
  real(dp),intent(in),optional :: qphon(3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out),optional :: str_vdw_dftd2(6)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: znucl(ntypat)
 end subroutine vdw_dftd2
end interface

interface
 subroutine vdw_dftd3(e_vdw_dftd3,ixc,natom,ntypat,prtvol,typat,rprimd,vdw_xc,&  
  &  vdw_tol,vdw_tol_3bt,xred,znucl,dyn_vdw_dftd3,elt_vdw_dftd3,&  
  &  fred_vdw_dftd3,str_vdw_dftd3,qphon)
  use defs_basis
  implicit none
  integer,intent(in) :: ixc
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  integer,intent(in) :: prtvol
  integer,intent(in) :: vdw_xc
  real(dp),intent(out) :: e_vdw_dftd3
  real(dp),intent(in) :: vdw_tol
  real(dp),intent(in) :: vdw_tol_3bt
  real(dp),intent(out),optional :: dyn_vdw_dftd3(2,3,natom,3,natom)
  real(dp),intent(out),optional :: elt_vdw_dftd3(6+3*natom,6)
  real(dp),intent(out),optional :: fred_vdw_dftd3(3,natom)
  real(dp),intent(in),optional :: qphon(3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out),optional :: str_vdw_dftd3(6)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: znucl(ntypat)
 end subroutine vdw_dftd3
end interface

interface
 subroutine vso_realspace_local(dtset,hdr,position_op,psps,vso_realspace)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  type(dataset_type),intent(in) :: dtset
  type(hdr_type),intent(inout) :: hdr
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: position_op(3,dtset%ngfft(1),dtset%ngfft(2),dtset%ngfft(3))
  real(dp),intent(out) :: vso_realspace(2,dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3), &
  &         dtset%nspinor,dtset%nspinor,3)
 end subroutine vso_realspace_local
end interface

interface
 subroutine vtorhotf(dtfil,dtset,ek,enl,entropy,fermie,gprimd,grnl,&  
  &  irrzon,mpi_enreg,natom,nfft,nspden,nsppol,nsym,phnons,rhog,rhor,rprimd,ucvol,vtrial)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(out) :: ek
  real(dp),intent(out) :: enl
  real(dp),intent(out) :: entropy
  real(dp),intent(out) :: fermie
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(out) :: grnl(3*natom)
  integer,intent(in) :: irrzon((dtset%ngfft(1)*dtset%ngfft(1)*dtset%ngfft(1))**(1-1/nsym), &
  &         2,(nspden/nsppol)-3*(nspden/4))
  real(dp),intent(in) :: phnons(2,(dtset%ngfft(1)*dtset%ngfft(1)*dtset%ngfft(1))**(1-1/nsym), &
  &         (nspden/nsppol)-3*(nspden/4))
  real(dp),intent(inout) :: rhog(2,nfft)
  real(dp),intent(inout) :: rhor(nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: vtrial(nfft,nspden)
 end subroutine vtorhotf
end interface

interface
 function zfermim12(xx)
  use defs_basis
  implicit none
  real(dp), intent(in) :: xx
  real(dp) :: zfermim12
 end function zfermim12
end interface

interface
 function zfermi12(xx)
  use defs_basis
  implicit none
  real(dp), intent(in) :: xx
  real(dp) :: zfermi12
 end function zfermi12
end interface

interface
 function zfermi1(xx)
  use defs_basis
  implicit none
  real(dp), intent(in) :: xx
  real(dp) :: zfermi1
 end function zfermi1
end interface

interface
 function zfermi32(xx)
  use defs_basis
  implicit none
  real(dp), intent(in) :: xx
  real(dp) :: zfermi32
 end function zfermi32
end interface

interface
 function zfermi2(xx)
  use defs_basis
  implicit none
  real(dp), intent(in) :: xx
  real(dp) :: zfermi2
 end function zfermi2
end interface

interface
 function zfermi52(xx)
  use defs_basis
  implicit none
  real(dp), intent(in) :: xx
  real(dp) :: zfermi52
 end function zfermi52
end interface

interface
 function zfermi3(xx)
  use defs_basis
  implicit none
  real(dp), intent(in) :: xx
  real(dp) :: zfermi3
 end function zfermi3
end interface

interface
 function ifermim12(ff)
  use defs_basis
  implicit none
  real(dp), intent(in) :: ff
  real(dp) :: ifermim12
 end function ifermim12
end interface

interface
 function ifermi12(ff)
  use defs_basis
  implicit none
  real(dp), intent(in) :: ff
  real(dp) :: ifermi12
 end function ifermi12
end interface

interface
 function ifermi32(ff)
  use defs_basis
  implicit none
  real(dp), intent(in) :: ff
  real(dp) :: ifermi32
 end function ifermi32
end interface

interface
 function ifermi52(ff)
  use defs_basis
  implicit none
  real(dp), intent(in) :: ff
  real(dp) :: ifermi52
 end function ifermi52
end interface

interface
 function fp12a1 (x)
  use defs_basis
  implicit none
  real(dp) :: fp12a1
  real(dp),intent(in) :: x
 end function fp12a1
end interface

interface
 function fp32a1 (x)
  use defs_basis
  implicit none
  real(dp) :: fp32a1
  real(dp),intent(in) :: x
 end function fp32a1
end interface

interface
 function xp12a1 (y)
  use defs_basis
  implicit none
  real(dp) :: xp12a1
  real(dp),intent(in) :: y
 end function xp12a1
end interface

interface
 function fm12a1 (x)
  use defs_basis
  implicit none
  real(dp) :: fm12a1
  real(dp),intent(in) :: x
 end function fm12a1
end interface

interface
 subroutine fm12a1t (cktf,rtnewt,tsmear,vtrial,rhor_middx,rhor_mid,nfft)
  use defs_basis
  implicit none
  integer,intent(in) :: nfft
  real(dp),intent(in) :: cktf
  real(dp),intent(in) :: rtnewt
  real(dp),intent(in) :: tsmear
  real(dp),intent(out) :: rhor_mid(nfft)
  real(dp),intent(out) :: rhor_middx(nfft)
  real(dp),intent(in) :: vtrial(nfft)
 end subroutine fm12a1t
end interface

interface
 subroutine waveformat(cg,cg_disk,cg_index,cg_new,dk,ii,ikpt,&  
  &  ikpt_,isgn,isppol,jj,jkpt,jkpt_,kg_kpt,kpt,kg_jl,maxband,mband,mcg,mcg_disk,&  
  &  minband,mkmem,mpw,nkpt,nkpt_,npwarr,nsppol,nspinor,shift_g_2,tr)
  use defs_basis
  implicit none
  integer,intent(in) :: ii
  integer,intent(in) :: ikpt
  integer,intent(in) :: ikpt_
  integer,intent(in) :: isgn
  integer,intent(in) :: isppol
  integer,intent(in) :: jj
  integer,intent(in) :: jkpt
  integer,intent(in) :: jkpt_
  integer,intent(in) :: maxband
  integer,intent(in) :: mband
  integer,intent(in) :: mcg
  integer,intent(in) :: mcg_disk
  integer,intent(in) :: minband
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: nkpt
  integer,intent(in) :: nkpt_
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  real(dp),intent(in) :: cg(2,mcg)
  real(dp),intent(in) :: cg_disk(2,mcg_disk,2)
  integer,intent(in) :: cg_index(mband,nkpt_,nsppol)
  real(dp),intent(out) :: cg_new(2,mpw,maxband)
  real(dp),intent(in) :: dk(3)
  integer,intent(in) :: kg_jl(3,mpw,2)
  integer,intent(in) :: kg_kpt(3,mpw*nspinor,nkpt_)
  real(dp),intent(in) :: kpt(3,nkpt)
  integer,intent(in) :: npwarr(nkpt_)
  logical,intent(in) :: shift_g_2(nkpt,nkpt)
  real(dp),intent(in) :: tr(2)
 end subroutine waveformat
end interface

interface
 subroutine wvl_initro(&  
  &  atindx1,geocode,h,me,&  
  &  natom,nattyp,nfft,nspden,ntypat,&  
  &  n1,n1i,n2,n2i,n3,&  
  &  pawrad,pawtab,psppar,&  
  &  rhor,rprimd,spinat,wvl_den,xc_denpos,xred,zion)
  use defs_basis
  use m_pawrad
  use defs_wvltypes
  use m_pawtab
  implicit none
  integer,intent(in) :: me
  integer,intent(in) :: n1
  integer,intent(in) :: n1i
  integer,intent(in) :: n2
  integer,intent(in) :: n2i
  integer,intent(in) :: n3
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  character(1),intent(in) :: geocode
  type(wvl_denspot_type), intent(inout) :: wvl_den
  real(dp),intent(in) :: xc_denpos
  integer,intent(in) :: atindx1(natom)
  real(dp),intent(in) :: h(3)
  integer,intent(in) :: nattyp(ntypat)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(in) :: psppar(0:4,0:6,ntypat)
  real(dp),intent(inout) :: rhor(nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: spinat(3,natom)
  real(dp),intent(inout) :: xred(3,natom)
  real(dp),intent(in) :: zion(ntypat)
 end subroutine wvl_initro
end interface

interface
 subroutine wvl_mkrho(dtset, irrzon, mpi_enreg, phnons, rhor, wvl_wfs, wvl_den)
  use defs_basis
  use defs_abitypes
  use defs_wvltypes
  implicit none
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(in) :: mpi_enreg
  type(wvl_denspot_type), intent(inout) :: wvl_den
  type(wvl_wf_type),intent(inout) :: wvl_wfs
  integer, target, intent(in) :: irrzon(dtset%nfft**(1-1/dtset%nsym), &
  &         2,(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  real(dp), target, intent(in) :: phnons(2,(dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3))**(1-1/dtset%nsym), &
  &         (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  real(dp),intent(inout) :: rhor(dtset%nfft,dtset%nspden)
 end subroutine wvl_mkrho
end interface

end module interfaces_67_common
!!***
