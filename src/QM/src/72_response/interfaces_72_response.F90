!!****m* ABINIT/interfaces_72_response
!! NAME
!! interfaces_72_response
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/72_response
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

module interfaces_72_response

 implicit none

interface
 subroutine corrmetalwf1(cgq,cprjq,cwavef,cwave1,cwaveprj,cwaveprj1,edocc,eig1,fermie1,ghc,iband,&  
  &  ibgq,icgq,istwf_k,mcgq,mcprjq,mpi_enreg,natom,nband,npw1,nspinor,occ,rocceig,timcount,&  
  &  usepaw,wf_corrected)
  use defs_basis
  use defs_abitypes
  use m_pawcprj
  implicit none
  integer,intent(in) :: iband
  integer,intent(in) :: ibgq
  integer,intent(in) :: icgq
  integer,intent(in) :: istwf_k
  integer,intent(in) :: mcgq
  integer,intent(in) :: mcprjq
  integer,intent(in) :: natom
  integer,intent(in) :: nband
  integer,intent(in) :: npw1
  integer,intent(in) :: nspinor
  integer,intent(in) :: timcount
  integer,intent(in) :: usepaw
  integer,intent(out) :: wf_corrected
  real(dp),intent(in) :: fermie1
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: cgq(2,mcgq)
  type(pawcprj_type),intent(in) :: cprjq(natom,mcprjq)
  real(dp),intent(out) :: cwave1(2,npw1*nspinor)
  real(dp),intent(in) :: cwavef(2,npw1*nspinor)
  type(pawcprj_type),intent(in) :: cwaveprj(natom,nspinor*usepaw)
  type(pawcprj_type),intent(inout) :: cwaveprj1(natom,nspinor*usepaw)
  real(dp),intent(out) :: edocc(nband)
  real(dp),intent(in) :: eig1(2*nband**2)
  real(dp),intent(in) :: ghc(2,npw1*nspinor)
  real(dp),intent(in) :: occ(nband)
  real(dp),intent(in) :: rocceig(nband,nband)
 end subroutine corrmetalwf1
end interface

interface
 subroutine d2dk_init(cg,d2dk,dtset,eig0_k,eig1_k,gs_hamkq,icg,idir,ikpt,ipert,isppol,mpi_enreg,mpw,&  
  nband_k,nsppol,rf_hamkq,rf_hamk_dir2,occ_k,rocceig,wffddk,ddk_f)
  use defs_abitypes
  use m_wffile
  use m_hamiltonian
  use defs_basis
  use m_d2dk
  use m_wfk
  implicit none
  integer,intent(in) :: icg
  integer,intent(in) :: idir
  integer,intent(in) :: ikpt
  integer,intent(in) :: ipert
  integer,intent(in) :: isppol
  integer,intent(in) :: mpw
  integer,intent(in) :: nband_k
  integer,intent(in) :: nsppol
  type(d2dk_t),intent(inout) :: d2dk
  type(dataset_type),intent(in) :: dtset
  type(gs_hamiltonian_type),intent(inout) :: gs_hamkq
  type(mpi_type),intent(inout) :: mpi_enreg
  type(rf_hamiltonian_type),intent(inout),target :: rf_hamk_dir2
  type(rf_hamiltonian_type),intent(inout),target :: rf_hamkq
  real(dp),intent(in) :: cg(2,mpw*gs_hamkq%nspinor*dtset%mband*nsppol)
  type(wfk_t),intent(inout) :: ddk_f(2)
  real(dp),intent(in) :: eig0_k(dtset%mband)
  real(dp),intent(inout) :: eig1_k(2*dtset%mband**2)
  real(dp),intent(in) :: occ_k(nband_k)
  real(dp),intent(in) :: rocceig(nband_k,nband_k)
  type(wffile_type),intent(inout) :: wffddk(2)
 end subroutine d2dk_init
end interface

interface
 subroutine d2frnl(becfrnl,cg,dtfil,dtset,dyfrnl,dyfr_cplex,dyfr_nondiag,efmasdeg,efmasfr,eigen,eltfrnl,&  
  &  gsqcut,has_allddk,indsym,kg,mgfftf,mpi_enreg,mpsang,my_natom,natom,nfftf,ngfft,ngfftf,npwarr,&  
  &  occ,paw_ij,pawang,pawbec,pawfgrtab,pawpiezo,pawrad,pawrhoij,pawtab,ph1d,ph1df,piezofrnl,psps,&  
  &  rprimd,rfphon,rfstrs,symrec,vtrial,vxc,xred,ylm,ylmgr)
  use m_pawrad
  use m_paw_ij
  use m_pawang
  use m_pawrhoij
  use defs_abitypes
  use m_efmas_defs
  use m_pawtab
  use defs_basis
  use defs_datatypes
  use m_pawfgrtab
  implicit none
  integer,intent(in) :: dyfr_cplex
  integer,intent(in) :: dyfr_nondiag
  integer,intent(in) :: mgfftf
  integer,intent(in) :: mpsang
  integer,intent(in) :: my_natom
  integer,intent(in) :: natom
  integer,intent(in) :: nfftf
  integer,intent(in) :: pawbec
  integer,intent(in) :: pawpiezo
  integer,intent(in) :: rfphon
  integer,intent(in) :: rfstrs
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: gsqcut
  logical,intent(inout) :: has_allddk
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pseudopotential_type),intent(in) :: psps
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: ngfftf(18)
  real(dp),intent(out) :: becfrnl(3,natom,3*pawbec)
  real(dp),intent(in) :: cg(2,dtset%mpw*dtset%nspinor*dtset%mband*dtset%mkmem*dtset%nsppol)
  real(dp),intent(out) :: dyfrnl(dyfr_cplex,3,3,natom,1+(natom-1)*dyfr_nondiag)
  type(efmasdeg_type),allocatable,intent(out) :: efmasdeg(:)
  type(efmasfr_type),allocatable,intent(out) :: efmasfr(:,:)
  real(dp),intent(in) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp),intent(out) :: eltfrnl(6+3*natom,6)
  integer,intent(in) :: indsym(4,dtset%nsym,natom)
  integer,intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
  integer,intent(in) :: npwarr(dtset%nkpt)
  real(dp),intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  type(paw_ij_type),intent(in) :: paw_ij(my_natom)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(my_natom*psps%usepaw)
  type(pawrad_type),intent(in) :: pawrad(psps%ntypat)
  type(pawrhoij_type),intent(inout),target :: pawrhoij(my_natom*psps%usepaw)
  type(pawtab_type),intent(in) :: pawtab(psps%ntypat)
  real(dp),intent(in) :: ph1d(2,3*(2*dtset%mgfft+1)*natom)
  real(dp),intent(in) :: ph1df(2,3*(2*mgfftf+1)*natom)
  real(dp),intent(out) :: piezofrnl(6,3*pawpiezo)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symrec(3,3,dtset%nsym)
  real(dp),intent(in),target :: vtrial(nfftf,dtset%nspden)
  real(dp),intent(in) :: vxc(nfftf,dtset%nspden)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: ylm(dtset%mpw*dtset%mkmem,mpsang*mpsang*psps%useylm)
  real(dp),intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,9,mpsang*mpsang*psps%useylm)
 end subroutine d2frnl
end interface

interface
 subroutine d2kindstr2(cwavef,ecut,ecutsm,effmass,ekinout,gmet,gprimd,&  
  &  istwfk,kg_k,kpt,npw,nspinor)
  use defs_basis
  implicit none
  integer,intent(in) :: istwfk
  integer,intent(in) :: npw
  integer,intent(in) :: nspinor
  real(dp),intent(in) :: ecut
  real(dp),intent(in) :: ecutsm
  real(dp),intent(in) :: effmass
  real(dp),intent(in) :: cwavef(2,npw*nspinor)
  real(dp),intent(inout) :: ekinout(36)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: kg_k(3,npw)
  real(dp),intent(in) :: kpt(3)
 end subroutine d2kindstr2
end interface

interface
 subroutine ddb_io_out (dscrpt,filnam,matom,mband,&  
  &  mkpt,msym,mtypat,unddb,vrsddb,&  
  &  acell,amu,dilatmx,ecut,ecutsm,intxc,iscf,ixc,kpt,kptnrm,&  
  &  natom,nband,ngfft,nkpt,nspden,nspinor,nsppol,nsym,ntypat,occ,occopt,&  
  &  pawecutdg,rprim,dfpt_sciss,spinat,symafm,symrel,tnons,tolwfr,tphysel,tsmear,&  
  &  typat,usepaw,wtk,xred,zion,znucl)
  use defs_basis
  implicit none
  integer,intent(in) :: intxc
  integer,intent(in) :: iscf
  integer,intent(in) :: ixc
  integer,intent(in) :: matom
  integer,intent(in) :: mband
  integer,intent(in) :: mkpt
  integer,intent(in) :: msym
  integer,intent(in) :: mtypat
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspden
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(in) :: occopt
  integer,intent(in) :: unddb
  integer,intent(in) :: usepaw
  integer,intent(in) :: vrsddb
  real(dp),intent(in) :: dfpt_sciss
  real(dp),intent(in) :: dilatmx
  character(len=fnlen),intent(in) :: dscrpt
  real(dp),intent(in) :: ecut
  real(dp),intent(in) :: ecutsm
  character(len=fnlen),intent(in) :: filnam
  real(dp),intent(in) :: kptnrm
  real(dp),intent(in) :: pawecutdg
  real(dp),intent(in) :: tolwfr
  real(dp),intent(in) :: tphysel
  real(dp),intent(in) :: tsmear
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: amu(mtypat)
  real(dp),intent(in) :: kpt(3,mkpt)
  integer,intent(in) :: nband(mkpt)
  real(dp),intent(in) :: occ(mband*mkpt)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(in) :: spinat(3,matom)
  integer,intent(in) :: symafm(msym)
  integer,intent(in) :: symrel(3,3,msym)
  real(dp),intent(in) :: tnons(3,msym)
  integer,intent(in) :: typat(matom)
  real(dp),intent(in) :: wtk(mkpt)
  real(dp),intent(in) :: xred(3,matom)
  real(dp),intent(in) :: zion(mtypat)
  real(dp),intent(in) :: znucl(mtypat)
 end subroutine ddb_io_out
end interface

interface
 subroutine ddb_write_blok(ddb,iblok,choice,mband,mpert,msize,nkpt,nunit,&  
  &  blkval2,kpt) !optional
  use defs_basis
  use m_ddb
  implicit none
  integer,intent(in) :: choice
  integer,intent(in) :: iblok
  integer,intent(in) :: mband
  integer,intent(in) :: mpert
  integer,intent(in) :: msize
  integer,intent(in) :: nkpt
  integer,intent(in) :: nunit
  type(ddb_type),intent(in) :: ddb
  real(dp),intent(in),optional :: blkval2(2,msize,mband,nkpt)
  real(dp),intent(in),optional :: kpt(3,nkpt)
 end subroutine ddb_write_blok
end interface

interface
 subroutine dfpt_accrho(counter,cplex,cwave0,cwave1,cwavef,cwaveprj0,cwaveprj1,&  
  &  eloc0_k,filstat,gs_hamkq,iband,idir,ipert,isppol,kptopt,&  
  &  mpi_enreg,natom,nband_k,ncpgr,npw_k,npw1_k,nspinor,occ_k,&  
  &  option,pawrhoij1,prtvol,rhoaug1,tim_fourwf,wf_corrected,&  
  &  wtk_k,comm_atom,mpi_atmtab)
  use defs_basis
  use defs_abitypes
  use m_pawcprj
  use m_pawrhoij
  use m_hamiltonian
  implicit none
  integer,optional,intent(in) :: comm_atom
  integer,intent(in) :: counter
  integer,intent(in) :: cplex
  integer,intent(in) :: iband
  integer,intent(in) :: idir
  integer,intent(in) :: ipert
  integer,intent(in) :: isppol
  integer,intent(in) :: kptopt
  integer,intent(in) :: natom
  integer,intent(in) :: nband_k
  integer,intent(in) :: ncpgr
  integer,intent(in) :: npw1_k
  integer,intent(in) :: npw_k
  integer,intent(in) :: nspinor
  integer,intent(in) :: option
  integer,intent(in) :: prtvol
  integer,intent(in) :: tim_fourwf
  integer,intent(in) :: wf_corrected
  real(dp),intent(out) :: eloc0_k
  character(len=fnlen),intent(in) :: filstat
  type(gs_hamiltonian_type),intent(inout),target :: gs_hamkq
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: wtk_k
  integer,optional,target,intent(in) :: mpi_atmtab(:)
  real(dp),intent(in),target :: cwave0(2,npw_k*nspinor)
  real(dp),intent(in),target :: cwave1(2,npw1_k*nspinor)
  real(dp),intent(in),target :: cwavef(2,npw1_k*nspinor)
  type(pawcprj_type),intent(in) :: cwaveprj0(natom,nspinor*gs_hamkq%usecprj)
  type(pawcprj_type),intent(in) :: cwaveprj1(natom,nspinor*gs_hamkq%usepaw)
  real(dp),intent(in) :: occ_k(nband_k)
  type(pawrhoij_type),intent(inout) :: pawrhoij1(:)
  real(dp),intent(inout) :: rhoaug1(cplex*gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6,gs_hamkq%nvloc)
 end subroutine dfpt_accrho
end interface

interface
 subroutine dfpt_cgwf(band,berryopt,cgq,cwavef,cwave0,cwaveprj,cwaveprj0,d2dk,dcwavef,&  
  &  eig0nk,eig0_kq,eig1_k,ghc,gh1c_n,grad_berry,gsc,gscq,&  
  &  gs_hamkq,gvnlc,gvnl1,icgq,idir,ipert,igscq,&  
  &  mcgq,mgscq,mpi_enreg,mpw1,natom,nband,nbdbuf,nline,npw,npw1,nspinor,&  
  &  opt_gvnl1,prtvol,quit,resid,rf_hamkq,dfpt_sciss,tolrde,tolwfr,&  
  &  usedcwavef,wfoptalg,nlines_done)
  use defs_basis
  use m_hamiltonian
  use defs_abitypes
  use m_pawcprj
  use m_d2dk
  implicit none
  integer,intent(in) :: band
  integer,intent(in) :: berryopt
  integer,intent(in) :: icgq
  integer,intent(in) :: idir
  integer,intent(in) :: igscq
  integer,intent(in) :: ipert
  integer,intent(in) :: mcgq
  integer,intent(in) :: mgscq
  integer,intent(in) :: mpw1
  integer,intent(in) :: natom
  integer,intent(in) :: nband
  integer,intent(in) :: nbdbuf
  integer,intent(in) :: nline
  integer,intent(inout) :: nlines_done
  integer,intent(in) :: npw
  integer,intent(in) :: npw1
  integer,intent(in) :: nspinor
  integer,intent(in) :: opt_gvnl1
  integer,intent(in) :: prtvol
  integer,intent(in) :: quit
  integer,intent(in) :: usedcwavef
  integer,intent(in) :: wfoptalg
  type(d2dk_t), intent(in) :: d2dk
  real(dp),intent(in) :: dfpt_sciss
  real(dp),intent(in) :: eig0nk
  type(gs_hamiltonian_type),intent(inout) :: gs_hamkq
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(out) :: resid
  type(rf_hamiltonian_type),intent(inout) :: rf_hamkq
  real(dp),intent(in) :: tolrde
  real(dp),intent(in) :: tolwfr
  real(dp),intent(in) :: cgq(2,mcgq)
  real(dp),intent(inout) :: cwave0(2,npw*nspinor)
  real(dp),intent(inout) :: cwavef(2,npw1*nspinor)
  type(pawcprj_type),intent(inout) :: cwaveprj(natom,nspinor)
  type(pawcprj_type),intent(inout) :: cwaveprj0(natom,nspinor*gs_hamkq%usecprj)
  real(dp),intent(inout) :: dcwavef(2,npw1*nspinor*((usedcwavef+1)/2))
  real(dp),intent(in) :: eig0_kq(nband)
  real(dp),intent(inout) :: eig1_k(2*nband**2)
  real(dp),intent(out) :: gh1c_n(2,npw1*nspinor)
  real(dp),intent(out) :: ghc(2,npw1*nspinor)
  real(dp),intent(in) :: grad_berry(2,mpw1*nspinor,nband)
  real(dp),intent(out) :: gsc(2,npw1*nspinor*gs_hamkq%usepaw)
  real(dp),intent(in) :: gscq(2,mgscq)
  real(dp),intent(inout) :: gvnl1(2,npw1*nspinor)
  real(dp),intent(inout) :: gvnlc(2,npw1*nspinor)
 end subroutine dfpt_cgwf
end interface

interface
 subroutine dfpt_dyfro(atindx1,dyfrnl,dyfrlo,dyfrwf,dyfrxc,dyfr_cplex,dyfr_nondiag,&  
  &  gmet,gprimd,gsqcut,indsym,mgfft,mpi_enreg,mqgrid,natom,nattyp,&  
  &  nfft,ngfft,nspden,nsym,ntypat,n1xccc,n3xccc,paral_kgb,psps,pawtab,ph1d,qgrid,&  
  &  qphon,rhog,rprimd,symq,symrec,typat,ucvol,usepaw,vlspl,vxc,&  
  &  xcccrc,xccc1d,xccc3d,xred)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use m_pawtab
  implicit none
  integer,intent(in) :: dyfr_cplex
  integer,intent(in) :: dyfr_nondiag
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
  integer,intent(in) :: usepaw
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: atindx1(natom)
  real(dp),intent(out) :: dyfrlo(3,3,natom)
  real(dp),intent(inout) :: dyfrnl(dyfr_cplex,3,3,natom,1+(natom-1)*dyfr_nondiag)
  real(dp),intent(out) :: dyfrwf(dyfr_cplex,3,3,natom,1+(natom-1)*dyfr_nondiag)
  real(dp),intent(inout) :: dyfrxc(3,3,natom)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: indsym(4,nsym,natom)
  integer,intent(in) :: nattyp(ntypat)
  type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symq(4,2,nsym)
  integer,intent(in) :: symrec(3,3,nsym)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: vlspl(mqgrid,2,ntypat)
  real(dp),intent(in) :: vxc(nfft,nspden)
  real(dp),intent(in) :: xccc1d(n1xccc,6,ntypat)
  real(dp),intent(inout) :: xccc3d(n3xccc)
  real(dp),intent(in) :: xcccrc(ntypat)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine dfpt_dyfro
end interface

interface
 subroutine dfpt_dyout(becfrnl,berryopt,blkflg,carflg,ddboun,ddkfil,dyew,dyfrlo,dyfrnl,&  
  &  dyfrx1,dyfrx2,dyfr_cplex,dyfr_nondiag,dyvdw,d2cart,d2cart_bbb,&  
  &  d2eig0,d2k0,d2lo,d2loc0,d2matr,d2nl,d2nl0,d2nl1,d2ovl,d2vn,&  
  &  eltcore,elteew,eltfrhar,eltfrkin,eltfrloc,eltfrnl,eltfrxc,eltvdw,&  
  &  has_full_piezo,has_allddk,iout,mband,mpert,natom,ntypat,&  
  &  outd2,pawbec,pawpiezo,piezofrnl,prtbbb,prtvol,qphon,qzero,typat,rfdir,&  
  &  rfpert,rfphon,rfstrs,usepaw,usevdw,zion)
  use defs_basis
  implicit none
  integer,intent(in) :: berryopt
  integer,intent(in) :: ddboun
  integer,intent(in) :: dyfr_cplex
  integer,intent(in) :: dyfr_nondiag
  integer,intent(in) :: iout
  integer,intent(in) :: mband
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  integer,intent(in) :: outd2
  integer,intent(in) :: pawbec
  integer,intent(in) :: pawpiezo
  integer,intent(in) :: prtbbb
  integer,intent(in) :: prtvol
  integer,intent(in) :: qzero
  integer, intent(in) :: rfphon
  integer, intent(in) :: rfstrs
  integer, intent(in) :: usepaw
  integer, intent(in) :: usevdw
  logical,intent(in) :: has_allddk
  logical,intent(in) :: has_full_piezo
  integer,intent(in) :: ddkfil(3)
  integer,intent(in) :: rfdir(3)
  real(dp),intent(in) :: becfrnl(3,natom,3*pawbec)
  integer,intent(in) :: blkflg(3,mpert,3,mpert)
  integer,intent(in) :: carflg(3,mpert,3,mpert)
  real(dp),intent(in) :: d2cart(2,3,mpert,3,mpert)
  real(dp),intent(inout) :: d2cart_bbb(2,3,3,mpert,mband,mband*prtbbb)
  real(dp),intent(in) :: d2eig0(2,3,mpert,3,mpert)
  real(dp),intent(in) :: d2k0(2,3,mpert,3,mpert)
  real(dp),intent(in) :: d2lo(2,3,mpert,3,mpert)
  real(dp),intent(in) :: d2loc0(2,3,mpert,3,mpert)
  real(dp),intent(in) :: d2matr(2,3,mpert,3,mpert)
  real(dp),intent(in) :: d2nl(2,3,mpert,3,mpert)
  real(dp),intent(in) :: d2nl0(2,3,mpert,3,mpert)
  real(dp),intent(in) :: d2nl1(2,3,mpert,3,mpert)
  real(dp),intent(in) :: d2ovl(2,3,mpert,3,mpert*usepaw)
  real(dp),intent(in) :: d2vn(2,3,mpert,3,mpert)
  real(dp),intent(in) :: dyew(2,3,natom,3,natom)
  real(dp),intent(in) :: dyfrlo(3,3,natom)
  real(dp),intent(in) :: dyfrnl(dyfr_cplex,3,3,natom,1+(natom-1)*dyfr_nondiag)
  real(dp),intent(in) :: dyfrx1(2,3,natom,3,natom)
  real(dp),intent(in) :: dyfrx2(3,3,natom)
  real(dp),intent(in) :: dyvdw(2,3,natom,3,natom*usevdw)
  real(dp),intent(in) :: eltcore(6,6)
  real(dp),intent(in) :: elteew(6+3*natom,6)
  real(dp),intent(in) :: eltfrhar(6,6)
  real(dp),intent(in) :: eltfrkin(6,6)
  real(dp),intent(in) :: eltfrloc(6+3*natom,6)
  real(dp),intent(in) :: eltfrnl(6+3*natom,6)
  real(dp),intent(in) :: eltfrxc(6+3*natom,6)
  real(dp),intent(in) :: eltvdw(6+3*natom,6*usevdw)
  real(dp),intent(in) :: piezofrnl(6,3*pawpiezo)
  real(dp),intent(in) :: qphon(3)
  integer,intent(in) :: rfpert(mpert)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: zion(ntypat)
 end subroutine dfpt_dyout
end interface

interface
 subroutine dfpt_dyxc1(atindx,blkflgfrx1,dyfrx1,gmet,gsqcut,ixc,kxc,mgfft,mpert,mpi_enreg,mqgrid,&  
  &  natom,nfft,ngfft,nkxc,nspden,ntypat,n1xccc,paral_kgb,psps,pawtab,&  
  &  ph1d,qgrid,qphon,rfdir,rfpert,rprimd,timrev,typat,ucvol,usepaw,xcccrc,xccc1d,xred,rhor)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use m_pawtab
  implicit none
  integer,intent(in) :: ixc
  integer,intent(in) :: mgfft
  integer,intent(in) :: mpert
  integer,intent(in) :: mqgrid
  integer,intent(in) :: n1xccc
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nkxc
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: timrev
  integer,intent(in) :: usepaw
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: rfdir(3)
  integer,intent(in) :: atindx(natom)
  integer,intent(out) :: blkflgfrx1(3,natom,3,natom)
  real(dp),intent(out) :: dyfrx1(2,3,natom,3,natom)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: kxc(nfft,nkxc)
  type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(in) :: qphon(3)
  integer,intent(in) :: rfpert(mpert)
  real(dp),optional,intent(in) :: rhor(nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xccc1d(n1xccc,6,ntypat)
  real(dp),intent(in) :: xcccrc(ntypat)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine dfpt_dyxc1
end interface

interface
 subroutine dfpt_eltfrhar(eltfrhar,rprimd,gsqcut,mpi_enreg,nfft,ngfft,rhog)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: nfft
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(in) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(out) :: eltfrhar(6,6)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine dfpt_eltfrhar
end interface

interface
 subroutine dfpt_eltfrkin(cg,eltfrkin,ecut,ecutsm,effmass,&  
  &  istwfk,kg,kptns,mband,mgfft,mkmem,mpi_enreg,&  
  &  mpw,nband,nkpt,ngfft,npwarr,nspinor,nsppol,occ,rprimd,wtk)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: mgfft
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  real(dp),intent(in) :: ecut
  real(dp),intent(in) :: ecutsm
  real(dp),intent(in) :: effmass
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  real(dp),intent(out) :: eltfrkin(6,6)
  integer,intent(in) :: istwfk(nkpt)
  integer,intent(in) :: kg(3,mpw*mkmem)
  real(dp),intent(in) :: kptns(3,nkpt)
  integer,intent(in) :: nband(nkpt*nsppol)
  integer,intent(in) :: npwarr(nkpt)
  real(dp),intent(in) :: occ(mband*nkpt*nsppol)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: wtk(nkpt)
 end subroutine dfpt_eltfrkin
end interface

interface
 subroutine dfpt_eltfrloc(atindx,eltfrloc,gmet,gprimd,gsqcut,mgfft,&  
  &  mpi_enreg,mqgrid,natom,nattyp,nfft,ngfft,ntypat,ph1d,qgrid,rhog,vlspl)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: mgfft
  integer,intent(in) :: mqgrid
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: ntypat
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(in) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: atindx(natom)
  real(dp),intent(out) :: eltfrloc(6+3*natom,6)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: nattyp(ntypat)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(in) :: vlspl(mqgrid,2,ntypat)
 end subroutine dfpt_eltfrloc
end interface

interface
 subroutine dfpt_eltfrxc(atindx,dtset,eltfrxc,enxc,gsqcut,kxc,mpi_enreg,mgfft,&  
  &  nattyp,nfft,ngfft,ngfftf,nhat,nkxc,n3xccc,pawtab,ph1d,psps,rhor,rprimd,&  
  &  usexcnhat,vxc,xccc3d,xred)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use m_pawtab
  implicit none
  integer,intent(in) :: mgfft
  integer,intent(in) :: n3xccc
  integer,intent(in) :: nfft
  integer,intent(in) :: nkxc
  integer,intent(in) :: usexcnhat
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: enxc
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(inout) :: psps
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: ngfftf(18)
  integer,intent(in) :: atindx(dtset%natom)
  real(dp),intent(out) :: eltfrxc(6+3*dtset%natom,6)
  real(dp),intent(inout) :: kxc(nfft,nkxc)
  integer,intent(in) :: nattyp(dtset%ntypat)
  real(dp),intent(in) :: nhat(nfft,dtset%nspden*psps%usepaw)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*dtset%usepaw)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*dtset%natom)
  real(dp),intent(in),target :: rhor(nfft,dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: vxc(nfft,dtset%nspden)
  real(dp),intent(in) :: xccc3d(n3xccc)
  real(dp),intent(in) :: xred(3,dtset%natom)
 end subroutine dfpt_eltfrxc
end interface

interface
 subroutine dfpt_etot(berryopt,deltae,eberry,edocc,eeig0,eew,efrhar,efrkin,efrloc,&  
  &  efrnl,efrx1,efrx2,ehart1,ek0,ek1,eii,elast,eloc0,elpsp1,&  
  &  enl0,enl1,epaw1,etotal,evar,evdw,exc1,ipert,natom,optene)
  use defs_basis
  implicit none
  integer,intent(in) :: berryopt
  integer,intent(in) :: ipert
  integer,intent(in) :: natom
  integer,intent(in) :: optene
  real(dp),intent(out) :: deltae
  real(dp),intent(in) :: eberry
  real(dp),intent(in) :: edocc
  real(dp),intent(in) :: eeig0
  real(dp),intent(in) :: eew
  real(dp),intent(in) :: efrhar
  real(dp),intent(in) :: efrkin
  real(dp),intent(in) :: efrloc
  real(dp),intent(in) :: efrnl
  real(dp),intent(in) :: efrx1
  real(dp),intent(in) :: efrx2
  real(dp),intent(in) :: ehart1
  real(dp),intent(in) :: eii
  real(dp),intent(in) :: ek0
  real(dp),intent(in) :: ek1
  real(dp),intent(inout) :: elast
  real(dp),intent(in) :: eloc0
  real(dp),intent(in) :: elpsp1
  real(dp),intent(in) :: enl0
  real(dp),intent(in) :: enl1
  real(dp),intent(in) :: epaw1
  real(dp),intent(out) :: etotal
  real(dp),intent(out) :: evar
  real(dp),intent(in) :: evdw
  real(dp),intent(in) :: exc1
 end subroutine dfpt_etot
end interface

interface
 subroutine dfpt_ewald(dyew,gmet,my_natom,natom,qphon,rmet,sumg0,typat,ucvol,xred,zion,&  
  &  mpi_atmtab,comm_atom ) ! optional arguments (parallelism))
  use defs_basis
  implicit none
  integer,optional,intent(in) :: comm_atom
  integer,intent(in) :: my_natom
  integer,intent(in) :: natom
  integer,intent(in) :: sumg0
  real(dp),intent(in) :: ucvol
  integer,optional,target,intent(in) :: mpi_atmtab(:)
  real(dp),intent(out) :: dyew(2,3,natom,3,natom)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(in) :: rmet(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: zion(*)
 end subroutine dfpt_ewald
end interface

interface
 subroutine dfpt_gatherdy(becfrnl,berryopt,blkflg,carflg,dyew,dyfrwf,dyfrx1,&  
  &  dyfr_cplex,dyfr_nondiag,dyvdw,d2bbb,d2cart,d2cart_bbb,d2matr,d2nfr,&  
  &  eltcore,elteew,eltfrhar,eltfrkin,eltfrloc,eltfrnl,eltfrxc,eltvdw,&  
  &  gprimd,mband,mpert,natom,ntypat,outd2,pawbec,pawpiezo,piezofrnl,prtbbb,&  
  &  rfasr,rfpert,rprimd,typat,ucvol,usevdw,zion)
  use defs_basis
  implicit none
  integer,intent(in) :: berryopt
  integer,intent(in) :: dyfr_cplex
  integer,intent(in) :: dyfr_nondiag
  integer,intent(in) :: mband
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  integer,intent(in) :: outd2
  integer,intent(in) :: pawbec
  integer,intent(in) :: pawpiezo
  integer,intent(in) :: prtbbb
  integer,intent(in) :: rfasr
  integer,intent(in) :: usevdw
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: becfrnl(3,natom,3*pawbec)
  integer,intent(inout) :: blkflg(3,mpert,3,mpert)
  integer,intent(out) :: carflg(3,mpert,3,mpert)
  real(dp),intent(in) :: d2bbb(2,3,3,mpert,mband,mband*prtbbb)
  real(dp),intent(out) :: d2cart(2,3,mpert,3,mpert)
  real(dp),intent(out) :: d2cart_bbb(2,3,3,mpert,mband,mband*prtbbb)
  real(dp),intent(out) :: d2matr(2,3,mpert,3,mpert)
  real(dp),intent(in) :: d2nfr(2,3,mpert,3,mpert)
  real(dp),intent(in) :: dyew(2,3,natom,3,natom)
  real(dp),intent(in) :: dyfrwf(dyfr_cplex,3,3,natom,1+(natom-1)*dyfr_nondiag)
  real(dp),intent(in) :: dyfrx1(2,3,natom,3,natom)
  real(dp),intent(in) :: dyvdw(2,3,natom,3,natom*usevdw)
  real(dp),intent(in) :: eltcore(6,6)
  real(dp),intent(in) :: elteew(6+3*natom,6)
  real(dp),intent(in) :: eltfrhar(6,6)
  real(dp),intent(in) :: eltfrkin(6,6)
  real(dp),intent(in) :: eltfrloc(6+3*natom,6)
  real(dp),intent(in) :: eltfrnl(6+3*natom,6)
  real(dp),intent(in) :: eltfrxc(6+3*natom,6)
  real(dp),intent(in) :: eltvdw(6+3*natom,6*usevdw)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: piezofrnl(6,3*pawpiezo)
  integer,intent(in) :: rfpert(mpert)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: zion(ntypat)
 end subroutine dfpt_gatherdy
end interface

interface
 subroutine dfpt_mkcore(cplex,idir,ipert,natom,ntypat,n1,n1xccc,&  
  &  n2,n3,qphon,rprimd,typat,ucvol,xcccrc,xccc1d,xccc3d1,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: idir
  integer,intent(in) :: ipert
  integer,intent(in) :: n1
  integer,intent(in) :: n1xccc
  integer,intent(in) :: n2
  integer,intent(in) :: n3
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xccc1d(n1xccc,6,ntypat)
  real(dp),intent(out) :: xccc3d1(cplex*n1*n2*n3)
  real(dp),intent(in) :: xcccrc(ntypat)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine dfpt_mkcore
end interface

interface
 subroutine dfpt_mkrho(cg,cg1,cplex,gprimd,irrzon,istwfk_rbz,&  
  &  kg,kg1,mband,mgfft,mkmem,mk1mem,mpi_enreg,mpw,mpw1,nband_rbz,&  
  &  nfft,ngfft,nkpt_rbz,npwarr,npwar1,nspden,nspinor,nsppol,nsym,&  
  &  occ_rbz,paral_kgb,phnons,rhog1,rhor1,rprimd,symafm,symrel,ucvol,wtk_rbz)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: mband
  integer,intent(in) :: mgfft
  integer,intent(in) :: mk1mem
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: mpw1
  integer,intent(in) :: nfft
  integer,intent(in) :: nkpt_rbz
  integer,intent(in) :: nspden
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym
  integer,intent(in) :: paral_kgb
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  real(dp),intent(in) :: cg1(2,mpw1*nspinor*mband*mk1mem*nsppol)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: irrzon(nfft**(1-1/nsym),2,(nspden/nsppol)-3*(nspden/4))
  integer,intent(in) :: istwfk_rbz(nkpt_rbz)
  integer,intent(in) :: kg(3,mpw*mkmem)
  integer,intent(in) :: kg1(3,mpw1*mk1mem)
  integer,intent(in) :: nband_rbz(nkpt_rbz*nsppol)
  integer,intent(in) :: npwar1(nkpt_rbz)
  integer,intent(in) :: npwarr(nkpt_rbz)
  real(dp),intent(in) :: occ_rbz(mband*nkpt_rbz*nsppol)
  real(dp),intent(in) :: phnons(2,nfft**(1-1/nsym),(nspden/nsppol)-3*(nspden/4))
  real(dp),intent(out) :: rhog1(2,nfft)
  real(dp),intent(out) :: rhor1(cplex*nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symafm(nsym)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: wtk_rbz(nkpt_rbz)
 end subroutine dfpt_mkrho
end interface

interface
 subroutine dfpt_mkvxcstr(cplex,idir,ipert,kxc,mpi_enreg,natom,nfft,ngfft,nhat,nhat1,&  
  &  nkxc,nspden,n3xccc,option,paral_kgb,qphon,rhor,rhor1,rprimd,usepaw,usexcnhat,vxc1,xccc3d1)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: idir
  integer,intent(in) :: ipert
  integer,intent(in) :: n3xccc
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nkxc
  integer,intent(in) :: nspden
  integer,intent(in) :: option
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: usepaw
  integer,intent(in) :: usexcnhat
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: kxc(nfft,nkxc)
  real(dp),target,intent(in) :: nhat(nfft,nspden)
  real(dp),target,intent(in) :: nhat1(cplex*nfft,nspden)
  real(dp),intent(in) :: qphon(3)
  real(dp),target,intent(in) :: rhor(nfft,nspden)
  real(dp),target,intent(in) :: rhor1(cplex*nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: vxc1(cplex*nfft,nspden)
  real(dp),intent(in) :: xccc3d1(cplex*n3xccc)
 end subroutine dfpt_mkvxcstr
end interface

interface
 subroutine dfpt_mkvxcstrgga(cplex,dgprimdds,gprimd,kxc,mpi_enreg,nfft,ngfft,&  
  &  nkxc,nspden,paral_kgb,qphon,rhor1tmp,vxc1)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: nfft
  integer,intent(in) :: nkxc
  integer,intent(in) :: nspden
  integer,intent(in) :: paral_kgb
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: dgprimdds(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: kxc(nfft,nkxc)
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(in) :: rhor1tmp(cplex*nfft,2)
  real(dp),intent(out) :: vxc1(cplex*nfft,nspden)
 end subroutine dfpt_mkvxcstrgga
end interface

interface
 subroutine dfpt_newvtr(cplex,dbl_nnsclo,dielar,dtset,etotal,ffttomix,&  
  &  initialized,iscf,ispmix,istep,mix,mixtofft,&  
  &  mpi_enreg,my_natom,nfft,nfftmix,ngfft,ngfftmix,npawmix,pawrhoij,&  
  &  qphon,rhor,rprimd,usepaw,vresid,vtrial)
  use defs_basis
  use defs_abitypes
  use m_ab7_mixing
  use m_pawrhoij
  implicit none
  integer,intent(in) :: cplex
  integer,intent(inout) :: dbl_nnsclo
  integer,intent(in) :: initialized
  integer,intent(in) :: iscf
  integer,intent(in) :: ispmix
  integer,intent(in) :: istep
  integer,intent(in) :: my_natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nfftmix
  integer,intent(in) :: npawmix
  integer,intent(in) :: usepaw
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: etotal
  type(ab7_mixing_object), intent(inout) :: mix
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: ngfftmix(18)
  real(dp),intent(in) :: dielar(7)
  integer,intent(in) :: ffttomix(nfft*(1-nfftmix/nfft))
  integer,intent(in) :: mixtofft(nfftmix*(1-nfftmix/nfft))
  type(pawrhoij_type),intent(inout) :: pawrhoij(my_natom*usepaw)
  real(dp),intent(in) :: qphon(3)
  real(dp), intent(in), target :: rhor(cplex*nfft,dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(inout) :: vresid(cplex*nfft,dtset%nspden)
  real(dp),intent(inout) :: vtrial(cplex*nfft,dtset%nspden)
 end subroutine dfpt_newvtr
end interface

interface
 subroutine dfpt_nselt(blkflg,cg,cg1,cplex,&  
  &  d2bbb,d2lo,d2nl,ecut,ecutsm,effmass,&  
  &  gmet,gprimd,gsqcut,idir,&  
  &  ipert,istwfk_rbz,kg,kg1,kpt_rbz,kxc,mband,mgfft,&  
  &  mkmem,mk1mem,mpert,mpi_enreg,mpsang,mpw,mpw1,&  
  &  natom,nband_rbz,nfft,ngfft,&  
  &  nkpt_rbz,nkxc,nloalg,npwarr,npwar1,nspden,nspinor,nsppol,&  
  &  nsym1,ntypat,occ_rbz,&  
  &  paral_kgb,ph1d,prtbbb,psps,qphon,rhog,&  
  &  rhor,rhor1,rmet,rprimd,symrc1,typat,ucvol,&  
  &  wtk_rbz,&  
  &  xred,ylm,ylm1,ylmgr,ylmgr1)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: idir
  integer,intent(in) :: ipert
  integer,intent(in) :: mband
  integer,intent(in) :: mgfft
  integer,intent(in) :: mk1mem
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpert
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpw
  integer,intent(in) :: mpw1
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nkpt_rbz
  integer,intent(in) :: nkxc
  integer,intent(in) :: nspden
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym1
  integer,intent(in) :: ntypat
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: prtbbb
  real(dp),intent(in) :: ecut
  real(dp),intent(in) :: ecutsm
  real(dp),intent(in) :: effmass
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: nloalg(3)
  integer,intent(inout) :: blkflg(3,mpert,3,mpert)
  real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  real(dp),intent(in) :: cg1(2,mpw1*nspinor*mband*mk1mem*nsppol)
  real(dp),intent(out) :: d2bbb(2,3,3,mpert,mband,mband*prtbbb)
  real(dp),intent(inout) :: d2lo(2,3,mpert,3,mpert)
  real(dp),intent(inout) :: d2nl(2,3,mpert,3,mpert)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: istwfk_rbz(nkpt_rbz)
  integer,intent(in) :: kg(3,mpw*mkmem)
  integer,intent(in) :: kg1(3,mpw1*mk1mem)
  real(dp),intent(in) :: kpt_rbz(3,nkpt_rbz)
  real(dp),intent(in) :: kxc(nfft,nkxc)
  integer,intent(in) :: nband_rbz(nkpt_rbz*nsppol)
  integer,intent(in) :: npwar1(nkpt_rbz)
  integer,intent(in) :: npwarr(nkpt_rbz)
  real(dp),intent(in) :: occ_rbz(mband*nkpt_rbz*nsppol)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(in) :: rhor(nfft,nspden)
  real(dp),intent(in) :: rhor1(cplex*nfft,nspden)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symrc1(3,3,nsym1)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: wtk_rbz(nkpt_rbz)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: ylm(mpw*mkmem,mpsang*mpsang*psps%useylm)
  real(dp),intent(in) :: ylm1(mpw1*mk1mem,mpsang*mpsang*psps%useylm)
  real(dp),intent(in) :: ylmgr(mpw*mkmem,3,mpsang*mpsang*psps%useylm)
  real(dp),intent(in) :: ylmgr1(mpw1*mk1mem,3,mpsang*mpsang*psps%useylm)
 end subroutine dfpt_nselt
end interface

interface
 subroutine dfpt_nstdy(atindx,blkflg,cg,cg1,cplex,dtfil,dtset,d2bbb,d2lo,d2nl,eigen0,eigen1,&  
  &  gmet,gsqcut,idir,indkpt1,indsy1,ipert,istwfk_rbz,kg,kg1,kpt_rbz,kxc,mkmem,mk1mem,&  
  &  mpert,mpi_enreg,mpw,mpw1,nattyp,nband_rbz,nfft,ngfft,nkpt,nkpt_rbz,nkxc,&  
  &  npwarr,npwar1,nspden,nsppol,nsym1,occ_rbz,ph1d,psps,rhor1,rmet,rprimd,&  
  &  symrc1,ucvol,wtk_rbz,xred,ylm,ylm1,rhor)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: idir
  integer,intent(in) :: ipert
  integer,intent(in) :: mk1mem
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpert
  integer,intent(in) :: mpw
  integer,intent(in) :: mpw1
  integer,intent(in) :: nfft
  integer,intent(in) :: nkpt
  integer,intent(in) :: nkpt_rbz
  integer,intent(in) :: nkxc
  integer,intent(in) :: nspden
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym1
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(inout) :: blkflg(3,mpert,3,mpert)
  real(dp),intent(in) :: cg(2,mpw*dtset%nspinor*dtset%mband*mkmem*nsppol)
  real(dp),intent(in) :: cg1(2,mpw1*dtset%nspinor*dtset%mband*mk1mem*nsppol)
  real(dp),intent(inout) :: d2bbb(2,3,3,mpert,dtset%mband,dtset%mband*dtset%prtbbb)
  real(dp),intent(inout) :: d2lo(2,3,mpert,3,mpert)
  real(dp),intent(inout) :: d2nl(2,3,mpert,3,mpert)
  real(dp),intent(in) :: eigen0(dtset%mband*nkpt_rbz*nsppol)
  real(dp),intent(in) :: eigen1(2*dtset%mband*dtset%mband*nkpt_rbz*nsppol)
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: indkpt1(nkpt_rbz)
  integer,intent(in) :: indsy1(4,nsym1,dtset%natom)
  integer,intent(in) :: istwfk_rbz(nkpt_rbz)
  integer,intent(in) :: kg(3,mpw*mkmem)
  integer,intent(in) :: kg1(3,mpw1*mk1mem)
  real(dp),intent(in) :: kpt_rbz(3,nkpt_rbz)
  real(dp),intent(in) :: kxc(nfft,nkxc)
  integer,intent(in) :: nattyp(dtset%ntypat)
  integer,intent(in) :: nband_rbz(nkpt_rbz*nsppol)
  integer,intent(in) :: npwar1(nkpt_rbz)
  integer,intent(in) :: npwarr(nkpt_rbz)
  real(dp),intent(in) :: occ_rbz(dtset%mband*nkpt_rbz*nsppol)
  real(dp),intent(in) :: ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom)
  real(dp),optional,intent(in) :: rhor(nfft,nspden)
  real(dp),intent(in) :: rhor1(cplex*nfft,nspden)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symrc1(3,3,nsym1)
  real(dp),intent(in) :: wtk_rbz(nkpt_rbz)
  real(dp),intent(in) :: xred(3,dtset%natom)
  real(dp),intent(in) :: ylm(mpw*mkmem,psps%mpsang*psps%mpsang*psps%useylm)
  real(dp),intent(in) :: ylm1(mpw1*mk1mem,psps%mpsang*psps%mpsang*psps%useylm)
 end subroutine dfpt_nstdy
end interface

interface
 subroutine dfpt_nsteltwf(cg,cg1,d2nl_k,ecut,ecutsm,effmass,gs_hamk,icg,icg1,ikpt,isppol,&  
  &  istwf_k,kg_k,kg1_k,kpoint,mband,mkmem,mk1mem,mpert,mpi_enreg,mpw,mpw1,natom,nband_k,&  
  &  npw_k,npw1_k,nspinor,nsppol,ntypat,occ_k,psps,rmet,wtk_k,ylm,ylmgr)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use m_hamiltonian
  implicit none
  integer,intent(in) :: icg
  integer,intent(in) :: icg1
  integer,intent(in) :: ikpt
  integer,intent(in) :: isppol
  integer,intent(in) :: istwf_k
  integer,intent(in) :: mband
  integer,intent(in) :: mk1mem
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpert
  integer,intent(in) :: mpw
  integer,intent(in) :: mpw1
  integer,intent(in) :: natom
  integer,intent(inout) :: nband_k
  integer,intent(inout) :: npw1_k
  integer,intent(inout) :: npw_k
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: ntypat
  real(dp),intent(in) :: ecut
  real(dp),intent(in) :: ecutsm
  real(dp),intent(in) :: effmass
  type(gs_hamiltonian_type) :: gs_hamk
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: wtk_k
  real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  real(dp),intent(in) :: cg1(2,mpw1*nspinor*mband*mk1mem*nsppol)
  real(dp),intent(inout) :: d2nl_k(2,3,mpert)
  integer,intent(in) :: kg1_k(3,npw1_k)
  integer,intent(in) :: kg_k(3,npw_k)
  real(dp),intent(in) :: kpoint(3)
  real(dp),intent(in) :: occ_k(nband_k)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: ylm(npw_k,psps%mpsang*psps%mpsang)
  real(dp),intent(in) :: ylmgr(npw_k,3,psps%mpsang*psps%mpsang)
 end subroutine dfpt_nsteltwf
end interface

interface
 subroutine dfpt_nstpaw(blkflg,cg,cgq,cg1,cplex,cprj,cprjq,docckqde,doccde_rbz,dtfil,dtset,d2lo,d2nl,d2ovl,&  
  &  eigenq,eigen0,eigen1,eovl1,gmet,gprimd,gsqcut,idir,indkpt1,indsy1,ipert,irrzon1,istwfk_rbz,&  
  &  kg,kg1,kpt_rbz,kxc,mgfftf,mkmem,mkqmem,mk1mem,&  
  &  mpert,mpi_enreg,mpw,mpw1,nattyp,nband_rbz,ncpgr,nfftf,ngfftf,nhat,nhat1,&  
  &  nkpt,nkpt_rbz,nkxc,npwarr,npwar1,nspden,nspinor,nsppol,nsym1,n3xccc,occkq,occ_rbz,&  
  &  paw_an,paw_an1,paw_ij,paw_ij1,pawang,pawang1,pawfgr,pawfgrtab,pawrad,pawrhoij,&  
  &  pawrhoij1,pawtab,phnons1,ph1d,ph1df,psps,rhog,rhor,rhor1,rmet,rprimd,symaf1,symrc1,symrl1,&  
  &  ucvol,usecprj,usepaw,usexcnhat,useylmgr1,vhartr1,vpsp1,vtrial,vtrial1,vxc,&  
  &  wtk_rbz,xccc3d1,xred,ylm,ylm1,ylmgr1)
  use m_pawtab
  use m_pawrad
  use m_paw_ij
  use m_pawang
  use m_pawrhoij
  use defs_abitypes
  use m_pawcprj
  use m_pawfgrtab
  use m_paw_an
  use m_pawfgr
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: idir
  integer,intent(in) :: ipert
  integer,intent(in) :: mgfftf
  integer,intent(in) :: mk1mem
  integer,intent(in) :: mkmem
  integer,intent(in) :: mkqmem
  integer,intent(in) :: mpert
  integer,intent(in) :: mpw
  integer,intent(in) :: mpw1
  integer,intent(in) :: n3xccc
  integer,intent(in) :: ncpgr
  integer,intent(in) :: nfftf
  integer,intent(in) :: nkpt
  integer,intent(in) :: nkpt_rbz
  integer,intent(in) :: nkxc
  integer,intent(in) :: nspden
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym1
  integer,intent(in) :: usecprj
  integer,intent(in) :: usepaw
  integer,intent(in) :: usexcnhat
  integer,intent(in) :: useylmgr1
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(out) :: eovl1
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pawang_type),intent(in) :: pawang1
  type(pawfgr_type),intent(in) :: pawfgr
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfftf(18)
  integer,intent(inout) :: blkflg(3,mpert,3,mpert)
  real(dp),intent(in) :: cg(2,mpw*nspinor*dtset%mband*mkmem*nsppol)
  real(dp),intent(in) :: cg1(2,mpw1*nspinor*dtset%mband*mk1mem*nsppol)
  real(dp),intent(in) :: cgq(2,mpw1*nspinor*dtset%mband*mkqmem*nsppol)
  type(pawcprj_type),intent(in) :: cprj(dtset%natom,nspinor*dtset%mband*mkmem*nsppol*usecprj)
  type(pawcprj_type),intent(in) :: cprjq(dtset%natom,nspinor*dtset%mband*mkqmem*nsppol*usecprj)
  real(dp),intent(inout) :: d2lo(2,3,mpert,3,mpert)
  real(dp),intent(inout) :: d2nl(2,3,mpert,3,mpert)
  real(dp),intent(inout) :: d2ovl(2,3,mpert,3,mpert*usepaw)
  real(dp),intent(in) :: doccde_rbz(dtset%mband*nkpt_rbz*nsppol)
  real(dp),intent(in) :: docckqde(dtset%mband*nkpt_rbz*nsppol)
  real(dp),intent(in) :: eigen0(dtset%mband*nkpt_rbz*nsppol)
  real(dp),intent(in) :: eigen1(2*dtset%mband*dtset%mband*nkpt_rbz*nsppol)
  real(dp),intent(in) :: eigenq(dtset%mband*nkpt_rbz*nsppol)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: indkpt1(nkpt_rbz)
  integer,intent(in) :: indsy1(4,nsym1,dtset%natom)
  integer,intent(in) :: irrzon1(dtset%nfft**(1-1/nsym1),2,(nspden/nsppol)-3*(nspden/4))
  integer,intent(in) :: istwfk_rbz(nkpt_rbz)
  integer,intent(in) :: kg(3,mpw*mkmem)
  integer,intent(in) :: kg1(3,mpw1*mk1mem)
  real(dp),intent(in) :: kpt_rbz(3,nkpt_rbz)
  real(dp),intent(in) :: kxc(nfftf,nkxc)
  integer,intent(in) :: nattyp(dtset%ntypat)
  integer,intent(in) :: nband_rbz(nkpt_rbz*nsppol)
  real(dp),intent(in) :: nhat(nfftf,nspden)
  real(dp),intent(in) :: nhat1(cplex*nfftf,nspden*usepaw)
  integer,intent(in) :: npwar1(nkpt_rbz)
  integer,intent(in) :: npwarr(nkpt_rbz)
  real(dp),intent(in) :: occ_rbz(dtset%mband*nkpt_rbz*nsppol)
  real(dp),intent(in) :: occkq(dtset%mband*nkpt_rbz*dtset%nsppol)
  type(paw_an_type),intent(in) :: paw_an(:)
  type(paw_an_type),intent(inout) :: paw_an1(:)
  type(paw_ij_type),intent(in) :: paw_ij(:)
  type(paw_ij_type),intent(inout) :: paw_ij1(:)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(:)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat*usepaw)
  type(pawrhoij_type),intent(in) :: pawrhoij(:)
  type(pawrhoij_type),intent(in) :: pawrhoij1(:)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*usepaw)
  real(dp),intent(in) :: ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom)
  real(dp),intent(in) :: ph1df(2,3*(2*mgfftf+1)*dtset%natom)
  real(dp),intent(in) :: phnons1(2,dtset%nfft**(1-1/nsym1),(nspden/nsppol)-3*(nspden/4))
  real(dp),intent(in) :: rhog(2,nfftf)
  real(dp),intent(in) :: rhor(cplex*nfftf,nspden)
  real(dp),intent(in) :: rhor1(cplex*nfftf,nspden)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symaf1(nsym1)
  integer,intent(in) :: symrc1(3,3,nsym1)
  integer,intent(in) :: symrl1(3,3,nsym1)
  real(dp),intent(in) :: vhartr1(cplex*nfftf)
  real(dp),target,intent(in) :: vpsp1(cplex*nfftf)
  real(dp),target,intent(in) :: vtrial(nfftf,nspden)
  real(dp),intent(in) :: vtrial1(cplex*nfftf,nspden)
  real(dp),intent(in) :: vxc(nfftf,nspden)
  real(dp),intent(in) :: wtk_rbz(nkpt_rbz)
  real(dp),target,intent(in) :: xccc3d1(cplex*n3xccc)
  real(dp),intent(in) :: xred(3,dtset%natom)
  real(dp),intent(in) :: ylm(mpw*mkmem,psps%mpsang*psps%mpsang*psps%useylm)
  real(dp),intent(in) :: ylm1(mpw1*mk1mem,psps%mpsang*psps%mpsang*psps%useylm)
  real(dp),intent(in) :: ylmgr1(mpw1*mk1mem,3,psps%mpsang*psps%mpsang*psps%useylm*useylmgr1)
 end subroutine dfpt_nstpaw
end interface

interface
 subroutine dfpt_nstwf(cg,cg1,ddkfil,dtset,d2bbb_k,d2nl_k,eig_k,eig1_k,gs_hamkq,&  
  &  icg,icg1,idir,ikpt,ipert,isppol,istwf_k,kg_k,kg1_k,kpt,kpq,&  
  &  mkmem,mk1mem,mpert,mpi_enreg,mpw,mpw1,nband_k,npw_k,npw1_k,nsppol,&  
  &  occ_k,psps,rmet,wffddk,ddks,wtk_k,ylm,ylm1)
  use defs_abitypes
  use m_wffile
  use m_hamiltonian
  use defs_basis
  use defs_datatypes
  use m_wfk
  implicit none
  integer,intent(in) :: icg
  integer,intent(in) :: icg1
  integer,intent(in) :: idir
  integer,intent(in) :: ikpt
  integer,intent(in) :: ipert
  integer,intent(in) :: isppol
  integer,intent(in) :: istwf_k
  integer,intent(in) :: mk1mem
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpert
  integer,intent(in) :: mpw
  integer,intent(in) :: mpw1
  integer,intent(inout) :: nband_k
  integer,intent(inout) :: npw1_k
  integer,intent(inout) :: npw_k
  integer,intent(in) :: nsppol
  type(dataset_type),intent(in) :: dtset
  type(gs_hamiltonian_type),intent(inout) :: gs_hamkq
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: wtk_k
  integer,intent(in) :: ddkfil(3)
  real(dp),intent(in) :: cg(2,mpw*dtset%nspinor*dtset%mband*mkmem*nsppol)
  real(dp),intent(in) :: cg1(2,mpw1*dtset%nspinor*dtset%mband*mk1mem*nsppol)
  real(dp),intent(out) :: d2bbb_k(2,3,dtset%mband,dtset%mband*dtset%prtbbb)
  real(dp),intent(inout) :: d2nl_k(2,3,mpert)
  type(wfk_t),intent(inout) :: ddks(3)
  real(dp),intent(inout) :: eig1_k(2*nsppol*dtset%mband**2)
  real(dp),intent(in) :: eig_k(dtset%mband*nsppol)
  integer,intent(in) :: kg1_k(3,npw1_k)
  integer,intent(in) :: kg_k(3,npw_k)
  real(dp),intent(in) :: kpq(3)
  real(dp),intent(in) :: kpt(3)
  real(dp),intent(in) :: occ_k(nband_k)
  real(dp),intent(in) :: rmet(3,3)
  type(wffile_type),intent(inout) :: wffddk(3)
  real(dp),intent(in) :: ylm(npw_k,psps%mpsang*psps%mpsang*psps%useylm)
  real(dp),intent(in) :: ylm1(npw1_k,psps%mpsang*psps%mpsang*psps%useylm)
 end subroutine dfpt_nstwf
end interface

interface
 subroutine dfpt_phfrq(amu,displ,d2cart,eigval,eigvec,indsym,&  
  &  mpert,msym,natom,nsym,ntypat,phfrq,qphnrm,qphon,rprimd,&  
  &  symdynmat,symrel,symafm,typat,ucvol)
  use defs_basis
  implicit none
  integer,intent(in) :: mpert
  integer,intent(in) :: msym
  integer,intent(in) :: natom
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(in) :: symdynmat
  real(dp),intent(in) :: qphnrm
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: amu(ntypat)
  real(dp),intent(in) :: d2cart(2,3,mpert,3,mpert)
  real(dp),intent(out) :: displ(2*3*natom*3*natom)
  real(dp),intent(out) :: eigval(3*natom)
  real(dp),intent(out) :: eigvec(2*3*natom*3*natom)
  integer,intent(in) :: indsym(4,msym*natom)
  real(dp),intent(out) :: phfrq(3*natom)
  real(dp),intent(inout) :: qphon(3)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symafm(nsym)
  integer,intent(in) :: symrel(3,3,nsym)
  integer,intent(in) :: typat(natom)
 end subroutine dfpt_phfrq
end interface

interface
 subroutine dfpt_prtene(berryopt,eberry,edocc,eeig0,eew,efrhar,efrkin,efrloc,efrnl,efrx1,efrx2,&  
  &  ehart01,ehart1,eii,ek0,ek1,eloc0,elpsp1,enl0,enl1,eovl1,epaw1,evdw,exc1,iout,ipert,natom,&  
  &  usepaw,usevdw)
  use defs_basis
  implicit none
  integer,intent(in) :: berryopt
  integer,intent(in) :: iout
  integer,intent(in) :: ipert
  integer,intent(in) :: natom
  integer,intent(in) :: usepaw
  integer,intent(in) :: usevdw
  real(dp),intent(in) :: eberry
  real(dp),intent(in) :: edocc
  real(dp),intent(in) :: eeig0
  real(dp),intent(in) :: eew
  real(dp),intent(in) :: efrhar
  real(dp),intent(in) :: efrkin
  real(dp),intent(in) :: efrloc
  real(dp),intent(in) :: efrnl
  real(dp),intent(in) :: efrx1
  real(dp),intent(in) :: efrx2
  real(dp),intent(in) :: ehart01
  real(dp),intent(in) :: ehart1
  real(dp),intent(in) :: eii
  real(dp),intent(in) :: ek0
  real(dp),intent(in) :: ek1
  real(dp),intent(in) :: eloc0
  real(dp),intent(in) :: elpsp1
  real(dp),intent(in) :: enl0
  real(dp),intent(in) :: enl1
  real(dp),intent(in) :: eovl1
  real(dp),intent(in) :: epaw1
  real(dp),intent(in) :: evdw
  real(dp),intent(in) :: exc1
 end subroutine dfpt_prtene
end interface

interface
 subroutine dfpt_prtph(displ,eivec,enunit,iout,natom,phfrq,qphnrm,qphon)
  use defs_basis
  implicit none
  integer,intent(in) :: eivec
  integer,intent(in) :: enunit
  integer,intent(in) :: iout
  integer,intent(in) :: natom
  real(dp),intent(in) :: qphnrm
  real(dp),intent(in) :: displ(2,3*natom,3*natom)
  real(dp),intent(in) :: phfrq(3*natom)
  real(dp),intent(in) :: qphon(3)
 end subroutine dfpt_prtph
end interface

interface
 subroutine dfpt_rhofermi(cg,cgq,cplex,cprj,cprjq,&  
  &  doccde_rbz,docckqde,dtfil,dtset,eigenq,eigen0,eigen1,fe1fixed,gmet,gprimd,idir,&  
  &  indsy1,ipert,irrzon1,istwfk_rbz,kg,kg1,kpt_rbz,mband,mkmem,mkqmem,mk1mem,mpi_enreg,&  
  &  mpw,mpw1,my_natom,natom,nband_rbz,ncpgr,nfftf,ngfftf,nhatfermi,nkpt_rbz,npwarr,npwar1,nspden,&  
  &  nsppol,nsym1,occkq,occ_rbz,paw_ij,pawang,pawang1,pawfgr,pawfgrtab,pawrad,pawrhoijfermi,pawtab,&  
  &  phnons1,ph1d,prtvol,psps,rhorfermi,rmet,rprimd,symaf1,symrc1,symrl1,&  
  &  ucvol,usecprj,useylmgr1,vtrial,vxc,wtk_rbz,xred,ylm,ylm1,ylmgr1)
  use m_pawtab
  use m_pawrad
  use m_paw_ij
  use m_pawang
  use m_pawrhoij
  use defs_abitypes
  use m_pawcprj
  use m_pawfgrtab
  use m_pawfgr
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: idir
  integer,intent(in) :: ipert
  integer,intent(in) :: mband
  integer,intent(in) :: mk1mem
  integer,intent(in) :: mkmem
  integer,intent(in) :: mkqmem
  integer,intent(in) :: mpw
  integer,intent(in) :: mpw1
  integer,intent(in) :: my_natom
  integer,intent(in) :: natom
  integer,intent(in) :: ncpgr
  integer,intent(in) :: nfftf
  integer,intent(in) :: nkpt_rbz
  integer,intent(in) :: nspden
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym1
  integer,intent(in) :: prtvol
  integer,intent(in) :: usecprj
  integer,intent(in) :: useylmgr1
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(out) :: fe1fixed
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pawang_type),intent(in) :: pawang1
  type(pawfgr_type),intent(in) :: pawfgr
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfftf(18)
  real(dp),intent(in) :: cg(2,mpw*dtset%nspinor*mband*mkmem*nsppol)
  real(dp),intent(in) :: cgq(2,mpw1*dtset%nspinor*mband*mkqmem*nsppol)
  type(pawcprj_type),intent(in) :: cprj(natom,dtset%nspinor*mband*mkmem *nsppol*usecprj)
  type(pawcprj_type),intent(in) :: cprjq(natom,dtset%nspinor*mband*mkqmem*nsppol*usecprj)
  real(dp),intent(in) :: doccde_rbz(mband*nkpt_rbz*nsppol)
  real(dp),intent(in) :: docckqde(mband*nkpt_rbz*nsppol)
  real(dp),intent(in) :: eigen0(mband*nkpt_rbz*nsppol)
  real(dp),intent(out) :: eigen1(2*mband*mband*nkpt_rbz*nsppol)
  real(dp),intent(in) :: eigenq(mband*nkpt_rbz*nsppol)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: indsy1(4,nsym1,natom)
  integer,intent(in) :: irrzon1(dtset%nfft**(1-1/nsym1),2,(nspden/nsppol)-3*(nspden/4))
  integer,intent(in) :: istwfk_rbz(nkpt_rbz)
  integer,intent(in) :: kg(3,mpw*mkmem)
  integer,intent(in) :: kg1(3,mpw1*mk1mem)
  real(dp),intent(in) :: kpt_rbz(3,nkpt_rbz)
  integer,intent(in) :: nband_rbz(nkpt_rbz*nsppol)
  real(dp),intent(out) :: nhatfermi(:,:)
  integer,intent(in) :: npwar1(nkpt_rbz)
  integer,intent(in) :: npwarr(nkpt_rbz)
  real(dp),intent(in) :: occ_rbz(mband*nkpt_rbz*nsppol)
  real(dp),intent(in) :: occkq(mband*nkpt_rbz*nsppol)
  type(paw_ij_type),intent(in) :: paw_ij(my_natom*psps%usepaw)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(my_natom*psps%usepaw)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat*psps%usepaw)
  type(pawrhoij_type),target,intent(inout) :: pawrhoijfermi(my_natom*psps%usepaw)
  type(pawtab_type), intent(in) :: pawtab(dtset%ntypat*psps%usepaw)
  real(dp),intent(in) :: ph1d(2,3*(2*dtset%mgfft+1)*natom)
  real(dp),intent(in) :: phnons1(2,dtset%nfft**(1-1/nsym1),(nspden/nsppol)-3*(nspden/4))
  real(dp),intent(out) :: rhorfermi(cplex*nfftf,nspden)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symaf1(nsym1)
  integer,intent(in) :: symrc1(3,3,nsym1)
  integer,intent(in) :: symrl1(3,3,nsym1)
  real(dp),intent(in) :: vtrial(nfftf,nspden)
  real(dp),intent(in) :: vxc(nfftf,nspden)
  real(dp),intent(in) :: wtk_rbz(nkpt_rbz)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: ylm(mpw*mkmem,psps%mpsang*psps%mpsang*psps%useylm)
  real(dp),intent(in) :: ylm1(mpw1*mk1mem,psps%mpsang*psps%mpsang*psps%useylm)
  real(dp),intent(in) :: ylmgr1(mpw1*mk1mem,3,psps%mpsang*psps%mpsang*psps%useylm*useylmgr1)
 end subroutine dfpt_rhofermi
end interface

interface
 subroutine dfpt_rhotov(cplex,ehart01,ehart1,elpsp1,exc1,gmet,gprimd,gsqcut,idir,ipert,&  
  &  ixc,kxc,mpi_enreg,natom,nfft,ngfft,nhat,nhat1,nhat1gr,nhat1grdim,nkxc,nspden,n3xccc,&  
  &  optene,optres,paral_kgb,qphon,rhog,rhog1,rhor,rhor1,rprimd,ucvol,&  
  &  usepaw,usexcnhat,vhartr1,vpsp1,vresid1,vres2,vtrial1,vxc1,xccc3d1)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: idir
  integer,intent(in) :: ipert
  integer,intent(in) :: ixc
  integer,intent(in) :: n3xccc
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nhat1grdim
  integer,intent(in) :: nkxc
  integer,intent(in) :: nspden
  integer,intent(in) :: optene
  integer,intent(in) :: optres
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: usepaw
  integer,intent(in) :: usexcnhat
  real(dp),intent(inout) :: ehart01
  real(dp),intent(inout) :: ehart1
  real(dp),intent(inout) :: elpsp1
  real(dp),intent(inout) :: exc1
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  real(dp),intent(out) :: vres2
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: kxc(nfft,nkxc)
  real(dp),intent(in) :: nhat(nfft,nspden)
  real(dp),intent(in) :: nhat1(cplex*nfft,nspden)
  real(dp),intent(in) :: nhat1gr(cplex*nfft,nspden,3*nhat1grdim)
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(in) :: rhog1(2,nfft)
  real(dp),target,intent(in) :: rhor(nfft,nspden)
  real(dp),target,intent(in) :: rhor1(cplex*nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),target,intent(out) :: vhartr1(:)
  real(dp),intent(in) :: vpsp1(cplex*nfft)
  real(dp),intent(out) :: vresid1(cplex*nfft,nspden)
  real(dp),intent(inout) :: vtrial1(cplex*nfft,nspden)
  real(dp),target,intent(out) :: vxc1(:,:)
  real(dp),intent(in) :: xccc3d1(cplex*n3xccc)
 end subroutine dfpt_rhotov
end interface

interface
 subroutine dfpt_scfcv(atindx,blkflg,cg,cgq,cg1,cg1_active,cplex,cprj,cprjq,cpus,&  
  &  dielt,dim_eig2rf,doccde_rbz,docckqde,dtfil,dtset,&  
  &  d2bbb,d2lo,d2nl,d2ovl,eberry,edocc,eeig0,eew,efrhar,efrkin,efrloc,efrnl,efrx1,efrx2,&  
  &  ehart01,ehart1,eigenq,eigen0,eigen1,eii,ek0,ek1,eloc0,elpsp1,&  
  &  enl0,enl1,eovl1,epaw1,etotal,evdw,exc1,fermie,gh0c1_set,gh1c_set,hdr,idir,indkpt1,&  
  &  indsy1,initialized,ipert,irrzon1,istwfk_rbz,&  
  &  kg,kg1,kpt_rbz,kxc,mgfftf,mkmem,mkqmem,mk1mem,&  
  &  mpert,mpi_enreg,mpw,mpw1,my_natom,nattyp,nband_rbz,ncpgr,&  
  &  nfftf,ngfftf,nhat,nkpt,nkpt_rbz,nkxc,npwarr,npwar1,nspden,&  
  &  nsym1,n3xccc,occkq,occ_rbz,&  
  &  paw_an,paw_ij,pawang,pawang1,pawfgr,pawfgrtab,pawrad,pawrhoij,pawrhoij1,pawtab,&  
  &  pertcase,phnons1,ph1d,ph1df,&  
  &  prtbbb,psps,qphon,resid,residm,rhog,rhog1,&  
  &  rhor,rhor1,rprimd,symaf1,symrc1,symrl1,&  
  &  usecprj,useylmgr,useylmgr1,wffddk,ddk_f,vpsp1,vtrial,vxc,&  
  &  wtk_rbz,xccc3d1,xred,ylm,ylm1,ylmgr,ylmgr1,zeff,conv_retcode)
  use m_pawtab
  use m_pawrad
  use m_paw_ij
  use m_pawang
  use m_pawrhoij
  use defs_abitypes
  use m_pawcprj
  use m_pawfgrtab
  use m_paw_an
  use m_pawfgr
  use defs_basis
  use m_wfk
  use defs_datatypes
  use m_wffile
  implicit none
  integer,intent(out) :: conv_retcode
  integer,intent(in) :: cplex
  integer,intent(in) :: dim_eig2rf
  integer,intent(in) :: idir
  integer,intent(inout) :: initialized
  integer,intent(in) :: ipert
  integer,intent(in) :: mgfftf
  integer,intent(in) :: mk1mem
  integer,intent(in) :: mkmem
  integer,intent(in) :: mkqmem
  integer,intent(in) :: mpert
  integer,intent(in) :: mpw
  integer,intent(in) :: mpw1
  integer,intent(in) :: my_natom
  integer,intent(in) :: n3xccc
  integer,intent(in) :: ncpgr
  integer,intent(in) :: nfftf
  integer,intent(in) :: nkpt
  integer,intent(in) :: nkpt_rbz
  integer,intent(in) :: nkxc
  integer,intent(in) :: nspden
  integer,intent(in) :: nsym1
  integer,intent(in) :: pertcase
  integer,intent(in) :: prtbbb
  integer,intent(in) :: usecprj
  integer,intent(in) :: useylmgr
  integer,intent(in) :: useylmgr1
  real(dp),intent(in) :: cpus
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(out) :: eberry
  real(dp),intent(out) :: edocc
  real(dp),intent(out) :: eeig0
  real(dp),intent(in) :: eew
  real(dp),intent(in) :: efrhar
  real(dp),intent(in) :: efrkin
  real(dp),intent(in) :: efrloc
  real(dp),intent(in) :: efrnl
  real(dp),intent(in) :: efrx1
  real(dp),intent(in) :: efrx2
  real(dp),intent(out) :: ehart01
  real(dp),intent(out) :: ehart1
  real(dp),intent(in) :: eii
  real(dp),intent(out) :: ek0
  real(dp),intent(out) :: ek1
  real(dp),intent(out) :: eloc0
  real(dp),intent(out) :: elpsp1
  real(dp),intent(out) :: enl0
  real(dp),intent(out) :: enl1
  real(dp),intent(out) :: eovl1
  real(dp),intent(out) :: epaw1
  real(dp),intent(out) :: etotal
  real(dp),intent(out) :: evdw
  real(dp),intent(out) :: exc1
  real(dp),intent(inout) :: fermie
  type(hdr_type),intent(inout) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pawang_type),intent(in) :: pawang1
  type(pawfgr_type),intent(in) :: pawfgr
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(out) :: residm
  integer,intent(in) :: ngfftf(18)
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(out) :: blkflg(3,mpert,3,mpert)
  real(dp),intent(in) :: cg(2,mpw*dtset%nspinor*dtset%mband*mkmem*dtset%nsppol)
  real(dp),intent(inout) :: cg1(2,mpw1*dtset%nspinor*dtset%mband*mk1mem*dtset%nsppol)
  real(dp),intent(out) :: cg1_active(2,mpw1*dtset%nspinor*dtset%mband*mk1mem*dtset%nsppol*dim_eig2rf)
  real(dp),intent(in) :: cgq(2,mpw1*dtset%nspinor*dtset%mband*mkqmem*dtset%nsppol)
  type(pawcprj_type),intent(in) :: cprj(dtset%natom,dtset%nspinor*dtset%mband*mkmem*dtset%nsppol*usecprj)
  type(pawcprj_type),intent(in) :: cprjq(dtset%natom,dtset%nspinor*dtset%mband*mkqmem*dtset%nsppol*usecprj)
  real(dp),intent(out) :: d2bbb(2,3,3,mpert,dtset%mband,dtset%mband*prtbbb)
  real(dp),intent(out) :: d2lo(2,3,mpert,3,mpert)
  real(dp),intent(out) :: d2nl(2,3,mpert,3,mpert)
  real(dp),intent(out) :: d2ovl(2,3,mpert,3,mpert*psps%usepaw)
  type(wfk_t),intent(inout) :: ddk_f(2)
  real(dp),intent(in) :: dielt(3,3)
  real(dp),intent(in) :: doccde_rbz(dtset%mband*nkpt_rbz*dtset%nsppol)
  real(dp),intent(in) :: docckqde(dtset%mband*nkpt_rbz*dtset%nsppol)
  real(dp),intent(in) :: eigen0(dtset%mband*nkpt_rbz*dtset%nsppol)
  real(dp),intent(out) :: eigen1(2*dtset%mband*dtset%mband*nkpt_rbz*dtset%nsppol)
  real(dp),intent(in) :: eigenq(dtset%mband*nkpt_rbz*dtset%nsppol)
  real(dp),intent(out) :: gh0c1_set(2,mpw1*dtset%nspinor*dtset%mband*mk1mem*dtset%nsppol*dim_eig2rf)
  real(dp),intent(out) :: gh1c_set(2,mpw1*dtset%nspinor*dtset%mband*mk1mem*dtset%nsppol*dim_eig2rf)
  integer,intent(in) :: indkpt1(nkpt_rbz)
  integer,intent(in) :: indsy1(4,nsym1,dtset%natom)
  integer,intent(in) :: irrzon1(dtset%nfft**(1-1/nsym1),2,(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  integer,intent(in) :: istwfk_rbz(nkpt_rbz)
  integer,intent(in) :: kg(3,mpw*mkmem)
  integer,intent(in) :: kg1(3,mpw1*mk1mem)
  real(dp),intent(in) :: kpt_rbz(3,nkpt_rbz)
  real(dp),intent(in) :: kxc(nfftf,nkxc)
  integer,intent(in) :: nattyp(psps%ntypat)
  integer,intent(in) :: nband_rbz(nkpt_rbz*dtset%nsppol)
  real(dp),intent(in) :: nhat(dtset%nfft,dtset%nspden)
  integer,intent(in) :: npwar1(nkpt_rbz)
  integer,intent(in) :: npwarr(nkpt_rbz)
  real(dp),intent(in) :: occ_rbz(dtset%mband*nkpt_rbz*dtset%nsppol)
  real(dp),intent(in) :: occkq(dtset%mband*nkpt_rbz*dtset%nsppol)
  type(paw_an_type),intent(in) :: paw_an(my_natom*psps%usepaw)
  type(paw_ij_type),intent(in) :: paw_ij(my_natom*psps%usepaw)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(my_natom*psps%usepaw)
  type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
  type(pawrhoij_type),intent(in) :: pawrhoij(my_natom*psps%usepaw)
  type(pawrhoij_type),intent(inout) :: pawrhoij1(my_natom*psps%usepaw)
  type(pawtab_type), intent(in) :: pawtab(psps%ntypat*psps%usepaw)
  real(dp),intent(in) :: ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom)
  real(dp),intent(in) :: ph1df(2,3*(2*mgfftf+1)*dtset%natom)
  real(dp),intent(in) :: phnons1(2,dtset%nfft**(1-1/nsym1),(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(out) :: resid(dtset%mband*nkpt_rbz*nspden)
  real(dp),intent(in) :: rhog(2,nfftf)
  real(dp),intent(inout) :: rhog1(2,nfftf)
  real(dp),intent(in) :: rhor(nfftf,nspden)
  real(dp),intent(inout) :: rhor1(cplex*nfftf,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symaf1(nsym1)
  integer,intent(in) :: symrc1(3,3,nsym1)
  integer,intent(in) :: symrl1(3,3,nsym1)
  real(dp),intent(in) :: vpsp1(cplex*nfftf)
  real(dp),target,intent(in) :: vtrial(nfftf,nspden)
  real(dp),intent(in) :: vxc(nfftf,nspden)
  type(wffile_type),intent(inout) :: wffddk(2)
  real(dp),intent(in) :: wtk_rbz(nkpt_rbz)
  real(dp),intent(in) :: xccc3d1(cplex*n3xccc)
  real(dp),intent(inout) :: xred(3,dtset%natom)
  real(dp),intent(in) :: ylm(mpw*mkmem,psps%mpsang*psps%mpsang*psps%useylm)
  real(dp),intent(in) :: ylm1(mpw1*mk1mem,psps%mpsang*psps%mpsang*psps%useylm)
  real(dp),intent(in) :: ylmgr(mpw*mkmem,3,psps%mpsang*psps%mpsang*psps%useylm*useylmgr)
  real(dp),intent(in) :: ylmgr1(mpw1*mk1mem,3+6*((ipert-dtset%natom)/10), &
  &         psps%mpsang*psps%mpsang*psps%useylm*useylmgr1)
  real(dp),intent(in) :: zeff(3,3,dtset%natom)
 end subroutine dfpt_scfcv
end interface

interface
 subroutine dfpt_sydy(cplex,dyfrow,indsym,natom,nondiag,nsym,qphon,sdyfro,symq,symrec)
  use defs_basis
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: natom
  integer,intent(in) :: nondiag
  integer,intent(in) :: nsym
  real(dp),intent(in) :: dyfrow(cplex,3,3,natom,1+(natom-1)*nondiag)
  integer,intent(in) :: indsym(4,nsym,natom)
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(out) :: sdyfro(cplex,3,3,natom,1+(natom-1)*nondiag)
  integer,intent(in) :: symq(4,2,nsym)
  integer,intent(in) :: symrec(3,3,nsym)
 end subroutine dfpt_sydy
end interface

interface
 subroutine dfpt_sygra(natom,desym,deunsy,indsym,ipert,nsym,qpt,symrec)
  use defs_basis
  implicit none
  integer,intent(in) :: ipert
  integer,intent(in) :: natom
  integer,intent(in) :: nsym
  real(dp),intent(out) :: desym(2,3,natom)
  real(dp),intent(in) :: deunsy(2,3,natom)
  integer,intent(in) :: indsym(4,nsym,natom)
  real(dp),intent(in) :: qpt(3)
  integer,intent(in) :: symrec(3,3,nsym)
 end subroutine dfpt_sygra
end interface

interface
 subroutine dfpt_symph(iout,acell,eigvec,indsym,natom,nsym,phfrq,rprim,symrel)
  use defs_basis
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: natom
  integer,intent(in) :: nsym
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: eigvec(2*3*natom*3*natom)
  integer,intent(in) :: indsym(4,nsym,natom)
  real(dp),intent(in) :: phfrq(3*natom)
  real(dp),intent(in) :: rprim(3,3)
  integer,intent(in) :: symrel(3,3,nsym)
 end subroutine dfpt_symph
end interface

interface
 subroutine dfpt_vlocal(atindx,cplex,gmet,gsqcut,idir,ipert,&  
  &  mpi_enreg,mqgrid,natom,nattyp,nfft,ngfft,&  
  &  ntypat,n1,n2,n3,paral_kgb,ph1d,qgrid,qphon,ucvol,vlspl,vpsp1,xred)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: idir
  integer,intent(in) :: ipert
  integer,intent(in) :: mqgrid
  integer,intent(in) :: n1
  integer,intent(in) :: n2
  integer,intent(in) :: n3
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: ntypat
  integer,intent(in) :: paral_kgb
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: atindx(natom)
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: nattyp(ntypat)
  real(dp),intent(in) :: ph1d(2,(2*n1+1+2*n2+1+2*n3+1)*natom)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(in) :: vlspl(mqgrid,2,ntypat)
  real(dp),intent(out) :: vpsp1(cplex*nfft)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine dfpt_vlocal
end interface

interface
 subroutine dfpt_vtorho(cg,cgq,cg1,cg1_active,cplex,cprj,cprjq,cprj1,dbl_nnsclo,&  
  &  dim_eig2rf,doccde_rbz,docckqde,dtefield,dtfil,dtset,&  
  &  edocc,eeig0,eigenq,eigen0,eigen1,ek0,ek1,eloc0,enl0,enl1,&  
  &  fermie1,gh0c1_set,gh1c_set,gmet,gprimd,idir,indsy1,&  
  &  ipert,irrzon1,istwfk_rbz,kg,kg1,kpt_rbz,mband,&  
  &  mkmem,mkqmem,mk1mem,mpi_enreg,mpw,mpw1,my_natom,&  
  &  natom,nband_rbz,ncpgr,nfftf,nhat1,nkpt_rbz,npwarr,npwar1,nres2,nspden,&  
  &  nsppol,nsym1,ntypat,nvresid1,occkq,occ_rbz,optres,&  
  &  paw_ij,paw_ij1,pawang,pawang1,pawfgr,pawfgrtab,pawrhoij,pawrhoij1,pawtab,&  
  &  phnons1,ph1d,prtvol,psps,pwindall,qmat,resid,residm,rhog1,rhor1,rmet,rprimd,symaf1,symrc1,symrl1,ucvol,&  
  &  usecprj,useylmgr1,wffddk,ddk_f,vtrial,vtrial1,wtk_rbz,xred,ylm,ylm1,ylmgr1,cg1_out)
  use m_pawtab
  use m_paw_ij
  use m_pawang
  use m_wfk
  use m_pawrhoij
  use defs_abitypes
  use m_pawcprj
  use m_pawfgrtab
  use m_pawfgr
  use defs_basis
  use m_efield
  use defs_datatypes
  use m_wffile
  implicit none
  integer,optional,intent(in) :: cg1_out
  integer,intent(in) :: cplex
  integer,intent(in) :: dbl_nnsclo
  integer,intent(in) :: dim_eig2rf
  integer,intent(in) :: idir
  integer,intent(in) :: ipert
  integer,intent(in) :: mband
  integer,intent(in) :: mk1mem
  integer,intent(in) :: mkmem
  integer,intent(in) :: mkqmem
  integer,intent(in) :: mpw
  integer,intent(in) :: mpw1
  integer,intent(in) :: my_natom
  integer,intent(in) :: natom
  integer,intent(in) :: ncpgr
  integer,intent(in) :: nfftf
  integer,intent(in) :: nkpt_rbz
  integer,intent(in) :: nspden
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym1
  integer,intent(in) :: ntypat
  integer,intent(in) :: optres
  integer,intent(in) :: prtvol
  integer,intent(in) :: usecprj
  integer,intent(in) :: useylmgr1
  type(efield_type),intent(in) :: dtefield
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(out) :: edocc
  real(dp),intent(out) :: eeig0
  real(dp),intent(out) :: ek0
  real(dp),intent(out) :: ek1
  real(dp),intent(out) :: eloc0
  real(dp),intent(out) :: enl0
  real(dp),intent(out) :: enl1
  real(dp),intent(in) :: fermie1
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(out) :: nres2
  type(pawang_type),intent(in) :: pawang
  type(pawang_type),intent(in) :: pawang1
  type(pawfgr_type),intent(in) :: pawfgr
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(out) :: residm
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: cg(2,mpw*dtset%nspinor*mband*mkmem*nsppol)
  real(dp),intent(inout) :: cg1(2,mpw1*dtset%nspinor*mband*mk1mem*nsppol)
  real(dp),intent(inout) :: cg1_active(2,mpw1*dtset%nspinor*mband*mk1mem*nsppol*dim_eig2rf)
  real(dp),intent(in) :: cgq(2,mpw1*dtset%nspinor*mband*mkqmem*nsppol)
  type(pawcprj_type),intent(in) :: cprj(natom,dtset%nspinor*mband*mkmem *nsppol*usecprj)
  type(pawcprj_type),intent(inout) :: cprj1(natom,dtset%nspinor*mband*mk1mem*nsppol*usecprj)
  type(pawcprj_type),intent(in) :: cprjq(natom,dtset%nspinor*mband*mkqmem*nsppol*usecprj)
  type(wfk_t),intent(inout) :: ddk_f(2)
  real(dp),intent(in) :: doccde_rbz(mband*nkpt_rbz*nsppol)
  real(dp),intent(in) :: docckqde(mband*nkpt_rbz*nsppol)
  real(dp),intent(in) :: eigen0(mband*nkpt_rbz*nsppol)
  real(dp),intent(out) :: eigen1(2*mband*mband*nkpt_rbz*nsppol)
  real(dp),intent(in) :: eigenq(mband*nkpt_rbz*nsppol)
  real(dp),intent(out) :: gh0c1_set(2,mpw1*dtset%nspinor*mband*mk1mem*nsppol*dim_eig2rf)
  real(dp),intent(out) :: gh1c_set(2,mpw1*dtset%nspinor*mband*mk1mem*nsppol*dim_eig2rf)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: indsy1(4,nsym1,natom)
  integer,intent(in) :: irrzon1(dtset%nfft**(1-1/nsym1),2,(nspden/nsppol)-3*(nspden/4))
  integer,intent(in) :: istwfk_rbz(nkpt_rbz)
  integer,intent(in) :: kg(3,mpw*mkmem)
  integer,intent(in) :: kg1(3,mpw1*mk1mem)
  real(dp),intent(in) :: kpt_rbz(3,nkpt_rbz)
  integer,intent(in) :: nband_rbz(nkpt_rbz*nsppol)
  real(dp), intent(out) :: nhat1(cplex*nfftf,dtset%nspden*psps%usepaw)
  integer,intent(in) :: npwar1(nkpt_rbz,2)
  integer,intent(in) :: npwarr(nkpt_rbz,2)
  real(dp),intent(inout) :: nvresid1(cplex*nfftf,nspden)
  real(dp),intent(in) :: occ_rbz(mband*nkpt_rbz*nsppol)
  real(dp),intent(in) :: occkq(mband*nkpt_rbz*nsppol)
  type(paw_ij_type),intent(in) :: paw_ij(my_natom*psps%usepaw)
  type(paw_ij_type),intent(in) :: paw_ij1(my_natom*psps%usepaw)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(my_natom*psps%usepaw)
  type(pawrhoij_type),intent(in) :: pawrhoij(my_natom*psps%usepaw)
  type(pawrhoij_type),target,intent(inout) :: pawrhoij1(my_natom*psps%usepaw)
  type(pawtab_type), intent(in) :: pawtab(ntypat*psps%usepaw)
  real(dp),intent(in) :: ph1d(2,3*(2*dtset%mgfft+1)*natom)
  real(dp),intent(in) :: phnons1(2,dtset%nfft**(1-1/nsym1),(nspden/nsppol)-3*(nspden/4))
  integer,intent(in) :: pwindall(max(mpw,mpw1)*mkmem,8,3)
  real(dp),intent(in) :: qmat(2,dtefield%nband_occ,dtefield%nband_occ,nkpt_rbz,2,3)
  real(dp),intent(out) :: resid(mband*nkpt_rbz*nsppol)
  real(dp),intent(out) :: rhog1(2,nfftf)
  real(dp),intent(inout) :: rhor1(cplex*nfftf,nspden)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symaf1(nsym1)
  integer,intent(in) :: symrc1(3,3,nsym1)
  integer,intent(in) :: symrl1(3,3,nsym1)
  real(dp),intent(in),target :: vtrial(nfftf,nspden)
  real(dp),intent(inout),target :: vtrial1(cplex*nfftf,nspden)
  type(wffile_type),intent(inout) :: wffddk(2)
  real(dp),intent(in) :: wtk_rbz(nkpt_rbz)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: ylm(mpw*mkmem,psps%mpsang*psps%mpsang*psps%useylm)
  real(dp),intent(in) :: ylm1(mpw1*mk1mem,psps%mpsang*psps%mpsang*psps%useylm)
  real(dp),intent(in) :: ylmgr1(mpw1*mk1mem,3+6*((ipert-natom)/10), &
  &         psps%mpsang*psps%mpsang*psps%useylm*useylmgr1)
 end subroutine dfpt_vtorho
end interface

interface
 subroutine dfpt_vtowfk(cg,cgq,cg1,cg1_active,cplex,cprj,cprjq,cprj1,&  
  &  dim_eig2rf,dtfil,dtset,&  
  &  edocc_k,eeig0_k,eig0_k,eig0_kq,eig1_k,&  
  &  ek0_k,ek1_k,eloc0_k,enl0_k,enl1_k,&  
  &  fermie1,gh0c1_set,gh1c_set,grad_berry,gs_hamkq,&  
  &  ibg,ibgq,ibg1,icg,icgq,icg1,idir,ikpt,ipert,&  
  &  isppol,mband,mcgq,mcprjq,mkmem,mk1mem,&  
  &  mpi_enreg,mpw,mpw1,natom,nband_k,ncpgr,&  
  &  nnsclo_now,npw_k,npw1_k,nspinor,nsppol,&  
  &  n4,n5,n6,occ_k,pawrhoij1,prtvol,psps,resid_k,rf_hamkq,rf_hamk_dir2,rhoaug1,rocceig,&  
  &  wffddk,ddk_f,wtk_k,nlines_done,cg1_out)
  use m_pawrhoij
  use defs_abitypes
  use m_pawcprj
  use m_wffile
  use m_hamiltonian
  use defs_basis
  use defs_datatypes
  use m_wfk
  implicit none
  integer,optional,intent(in) :: cg1_out
  integer,intent(in) :: cplex
  integer,intent(in) :: dim_eig2rf
  integer,intent(in) :: ibg
  integer,intent(in) :: ibg1
  integer,intent(in) :: ibgq
  integer,intent(in) :: icg
  integer,intent(in) :: icg1
  integer,intent(in) :: icgq
  integer,intent(in) :: idir
  integer,intent(in) :: ikpt
  integer,intent(in) :: ipert
  integer,intent(in) :: isppol
  integer,intent(in) :: mband
  integer,intent(in) :: mcgq
  integer,intent(in) :: mcprjq
  integer,intent(in) :: mk1mem
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: mpw1
  integer,intent(in) :: n4
  integer,intent(in) :: n5
  integer,intent(in) :: n6
  integer,intent(in) :: natom
  integer,intent(in) :: nband_k
  integer,intent(in) :: ncpgr
  integer,intent(inout) :: nlines_done
  integer,intent(in) :: nnsclo_now
  integer,intent(in) :: npw1_k
  integer,intent(in) :: npw_k
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: prtvol
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: fermie1
  type(gs_hamiltonian_type),intent(inout) :: gs_hamkq
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  type(rf_hamiltonian_type),intent(inout) :: rf_hamk_dir2
  type(rf_hamiltonian_type),intent(inout) :: rf_hamkq
  real(dp),intent(in) :: wtk_k
  real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  real(dp),intent(inout) :: cg1(2,mpw1*nspinor*mband*mk1mem*nsppol)
  real(dp),intent(inout) :: cg1_active(2,mpw1*nspinor*mband*mk1mem*nsppol*dim_eig2rf)
  real(dp),intent(in) :: cgq(2,mcgq)
  type(pawcprj_type),intent(in) :: cprj(natom,nspinor*mband*mkmem*nsppol*gs_hamkq%usecprj)
  type(pawcprj_type),intent(inout) :: cprj1(natom,nspinor*mband*mk1mem*nsppol*gs_hamkq%usecprj)
  type(pawcprj_type),intent(in) :: cprjq(natom,mcprjq)
  type(wfk_t),intent(inout) :: ddk_f(2)
  real(dp),intent(inout) :: edocc_k(nband_k)
  real(dp),intent(inout) :: eeig0_k(nband_k)
  real(dp),intent(in) :: eig0_k(nband_k)
  real(dp),intent(in) :: eig0_kq(nband_k)
  real(dp),intent(inout) :: eig1_k(2*nband_k**2)
  real(dp),intent(out) :: ek0_k(nband_k)
  real(dp),intent(inout) :: ek1_k(nband_k)
  real(dp),intent(out) :: eloc0_k(nband_k)
  real(dp),intent(out) :: enl0_k(nband_k)
  real(dp),intent(out) :: enl1_k(nband_k)
  real(dp),intent(inout) :: gh0c1_set(2,mpw1*nspinor*mband*mk1mem*nsppol*dim_eig2rf)
  real(dp),intent(inout) :: gh1c_set(2,mpw1*nspinor*mband*mk1mem*nsppol*dim_eig2rf)
  real(dp),intent(in) :: grad_berry(2,mpw1*nspinor,nband_k)
  real(dp),intent(in) :: occ_k(nband_k)
  type(pawrhoij_type),intent(inout) :: pawrhoij1(natom*gs_hamkq%usepaw)
  real(dp),intent(out) :: resid_k(nband_k)
  real(dp),intent(inout) :: rhoaug1(cplex*n4,n5,n6,gs_hamkq%nvloc)
  real(dp),intent(in) :: rocceig(nband_k,nband_k)
  type(wffile_type),intent(inout) :: wffddk(2)
 end subroutine dfpt_vtowfk
end interface

interface
 subroutine dfpt_wfkfermi(cg,cgq,cplex,cprj,cprjq,&  
  &  dtfil,eig0_k,eig1_k,fe1fixed_k,fe1norm_k,gs_hamkq,&  
  &  ibg,ibgq,icg,icgq,idir,ikpt,ipert,isppol,&  
  &  kptopt,mband,mcgq,mcprjq,mkmem,mpi_enreg,mpw,nband_k,ncpgr,&  
  &  npw_k,npw1_k,nspinor,nsppol,occ_k,pawrhoijfermi,prtvol,&  
  &  rf_hamkq,rhoaug,rocceig,wtk_k)
  use defs_basis
  use defs_abitypes
  use m_pawcprj
  use m_pawrhoij
  use m_hamiltonian
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: ibg
  integer,intent(in) :: ibgq
  integer,intent(in) :: icg
  integer,intent(in) :: icgq
  integer,intent(in) :: idir
  integer,intent(in) :: ikpt
  integer,intent(in) :: ipert
  integer,intent(in) :: isppol
  integer,intent(in) :: kptopt
  integer,intent(in) :: mband
  integer,intent(in) :: mcgq
  integer,intent(in) :: mcprjq
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(inout) :: nband_k
  integer,intent(in) :: ncpgr
  integer,intent(in) :: npw1_k
  integer,intent(inout) :: npw_k
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: prtvol
  type(datafiles_type),intent(in) :: dtfil
  type(gs_hamiltonian_type),intent(inout) :: gs_hamkq
  type(mpi_type),intent(inout) :: mpi_enreg
  type(rf_hamiltonian_type),intent(inout) :: rf_hamkq
  real(dp),intent(in) :: wtk_k
  real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  real(dp),intent(in) :: cgq(2,mcgq)
  type(pawcprj_type),intent(in) :: cprj(gs_hamkq%natom,nspinor*mband*mkmem*nsppol*gs_hamkq%usecprj)
  type(pawcprj_type),intent(in) :: cprjq(gs_hamkq%natom,mcprjq)
  real(dp),intent(in) :: eig0_k(nband_k)
  real(dp),intent(inout) :: eig1_k(2*nband_k**2)
  real(dp),intent(out) :: fe1fixed_k(nband_k)
  real(dp),intent(out) :: fe1norm_k(nband_k)
  real(dp),intent(in) :: occ_k(nband_k)
  type(pawrhoij_type),intent(inout) :: pawrhoijfermi(gs_hamkq%natom*gs_hamkq%usepaw)
  real(dp),intent(inout) :: rhoaug(cplex*gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6)
  real(dp),intent(in) :: rocceig(nband_k,nband_k)
 end subroutine dfpt_wfkfermi
end interface

interface
 subroutine dfptff_bec(cg,cg1,dtefield,natom,d2lo,idirpert,ipert,mband,mkmem,&  
  &  mpw,mpw1,mpert,nkpt,npwarr,npwar1,nsppol,nspinor,pwindall,qmat,rprimd)
  use defs_basis
  use m_efield
  implicit none
  integer,intent(in) :: idirpert
  integer,intent(in) :: ipert
  integer,intent(in) :: mband
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpert
  integer,intent(in) :: mpw
  integer,intent(in) :: mpw1
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  type(efield_type),intent(in) :: dtefield
  real(dp),intent(in) :: cg(2,mpw*mband*mkmem*nspinor*nsppol)
  real(dp),intent(in) :: cg1(2,mpw1*mband*mkmem*nspinor*nsppol)
  real(dp),intent(inout) :: d2lo(2,3,mpert,3,mpert)
  integer,intent(in) :: npwar1(nkpt)
  integer,intent(in) :: npwarr(nkpt)
  integer,intent(in) :: pwindall(max(mpw,mpw1)*mkmem,8,3)
  real(dp),intent(in) :: qmat(2,dtefield%nband_occ,dtefield%nband_occ,nkpt,2,3)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine dfptff_bec
end interface

interface
 subroutine dfptff_die(cg,cg1,dtefield,d2lo,idirpert,ipert,mband,mkmem,&  
  &  mpw,mpw1,mpert,nkpt,npwarr,npwar1,nsppol,nspinor,pwindall,qmat,rprimd)
  use defs_basis
  use m_efield
  implicit none
  integer,intent(in) :: idirpert
  integer,intent(in) :: ipert
  integer,intent(in) :: mband
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpert
  integer,intent(in) :: mpw
  integer,intent(in) :: mpw1
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  type(efield_type),intent(in) :: dtefield
  real(dp),intent(in) :: cg(2,mpw*mband*mkmem*nspinor*nsppol)
  real(dp),intent(in) :: cg1(2,mpw1*mband*mkmem*nspinor*nsppol)
  real(dp),intent(inout) :: d2lo(2,3,mpert,3,mpert)
  integer,intent(in) :: npwar1(nkpt)
  integer,intent(in) :: npwarr(nkpt)
  integer,intent(in) :: pwindall(max(mpw,mpw1)*mkmem,8,3)
  real(dp),intent(in) :: qmat(2,dtefield%nband_occ,dtefield%nband_occ,nkpt,2,3)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine dfptff_die
end interface

interface
 subroutine dfptff_ebp(cg,cg1,dtefield,eberry,mband,mkmem,&  
  &  mpw,mpw1,nkpt,npwarr,npwar1,nsppol,nspinor,pwindall,qmat)
  use defs_basis
  use m_efield
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: mpw1
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  type(efield_type),intent(in) :: dtefield
  real(dp),intent(out) :: eberry
  real(dp),intent(in) :: cg(2,mpw*mband*mkmem*nspinor*nsppol)
  real(dp),intent(in) :: cg1(2,mpw1*mband*mkmem*nspinor*nsppol)
  integer,intent(in) :: npwar1(nkpt)
  integer,intent(in) :: npwarr(nkpt)
  integer,intent(in) :: pwindall(max(mpw,mpw1)*mkmem,8,3)
  real(dp),intent(in) :: qmat(2,dtefield%nband_occ,dtefield%nband_occ,nkpt,2,3)
 end subroutine dfptff_ebp
end interface

interface
 subroutine dfptff_edie(cg,cg1,dtefield,eberry,idir_efield,mband,mkmem,&  
  &  mpw,mpw1,nkpt,npwarr,npwar1,nsppol,nspinor,pwindall,qmat,rprimd)
  use defs_basis
  use m_efield
  implicit none
  integer,intent(in) :: idir_efield
  integer,intent(in) :: mband
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: mpw1
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  type(efield_type),intent(in) :: dtefield
  real(dp),intent(out) :: eberry
  real(dp),intent(in) :: cg(2,mpw*mband*mkmem*nspinor*nsppol)
  real(dp),intent(in) :: cg1(2,mpw1*mband*mkmem*nspinor*nsppol)
  integer,intent(in) :: npwar1(nkpt)
  integer,intent(in) :: npwarr(nkpt)
  integer,intent(in) :: pwindall(max(mpw,mpw1)*mkmem,8,3)
  real(dp),intent(in) :: qmat(2,dtefield%nband_occ,dtefield%nband_occ,nkpt,2,3)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine dfptff_edie
end interface

interface
 subroutine dfptff_gbefd(cg,cg1,dtefield,grad_berry,idir_efield,ikpt,isppol,mband,mpw,mpw1,mkmem,mk1mem,nkpt,&  
  &  npwarr,npwar1,nspinor,nsppol,qmat,pwindall,rprimd)
  use defs_basis
  use m_efield
  implicit none
  integer,intent(in) :: idir_efield
  integer,intent(in) :: ikpt
  integer,intent(in) :: isppol
  integer,intent(in) :: mband
  integer,intent(in) :: mk1mem
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: mpw1
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  type(efield_type),intent(in) :: dtefield
  real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  real(dp),intent(in) :: cg1(2,mpw1*nspinor*mband*mk1mem*nsppol)
  real(dp),intent(out) :: grad_berry(2,mpw1,dtefield%nband_occ)
  integer,intent(in) :: npwar1(nkpt)
  integer,intent(in) :: npwarr(nkpt)
  integer,intent(in) :: pwindall(max(mpw,mpw1)*mkmem,8,3)
  real(dp),intent(in) :: qmat(2,dtefield%nband_occ,dtefield%nband_occ,nkpt,2,3)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine dfptff_gbefd
end interface

interface
 subroutine dfptff_gradberry(cg,cg1,dtefield,grad_berry,ikpt,isppol,mband,mpw,mpw1,mkmem,mk1mem,nkpt,&  
  &  npwarr,npwar1,nspinor,nsppol,qmat,pwindall)
  use defs_basis
  use m_efield
  implicit none
  integer,intent(in) :: ikpt
  integer,intent(in) :: isppol
  integer,intent(in) :: mband
  integer,intent(in) :: mk1mem
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: mpw1
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  type(efield_type),intent(in) :: dtefield
  real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  real(dp),intent(in) :: cg1(2,mpw1*nspinor*mband*mk1mem*nsppol)
  real(dp),intent(out) :: grad_berry(2,mpw1,dtefield%nband_occ)
  integer,intent(in) :: npwar1(nkpt)
  integer,intent(in) :: npwarr(nkpt)
  integer,intent(in) :: pwindall(max(mpw,mpw1)*mkmem,8,3)
  real(dp),intent(in) :: qmat(2,dtefield%nband_occ,dtefield%nband_occ,nkpt,2,3)
 end subroutine dfptff_gradberry
end interface

interface
 subroutine  dfptff_initberry(dtefield,dtset,gmet,kg,kg1,mband,mkmem,mpi_enreg,&  
  &  mpw,mpw1,nkpt,npwarr,npwar1,nsppol,occ,pwindall,rprimd)
  use defs_basis
  use m_efield
  use defs_abitypes
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: mpw1
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  type(efield_type),intent(inout) :: dtefield
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: kg(3,mpw*mkmem)
  integer,intent(in) :: kg1(3,mpw1*mkmem)
  integer,intent(in) :: npwar1(nkpt)
  integer,intent(in) :: npwarr(nkpt)
  real(dp),intent(in) :: occ(mband*nkpt*nsppol)
  integer,intent(out) :: pwindall(max(mpw,mpw1)*mkmem,8,3)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine dfptff_initberry
end interface

interface
 subroutine dfptnl_doutput(blkflg,d3,mband,mpert,nkpt,natom,ntypat,unddb)
  use defs_basis
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt
  integer,intent(in) :: ntypat
  integer,intent(in) :: unddb
  integer,intent(in) :: blkflg(3,mpert,3,mpert,3,mpert)
  real(dp),intent(in) :: d3(2,3,mpert,3,mpert,3,mpert)
 end subroutine dfptnl_doutput
end interface

interface
 subroutine dfptnl_resp(cg,cg1,cg3,cplex,dtfil,dtset,d3lo,&  
  &  i1dir,i2dir,i3dir,i1pert,i2pert,i3pert,&  
  &  kg,mband,mgfft,mkmem,mk1mem,&  
  &  mpert,mpi_enreg,mpsang,mpw,natom,nfft,nkpt,nspden,nspinor,nsppol,&  
  &  npwarr,occ,ph1d,psps,rprimd,vtrial1,xred,ylm)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: i1dir
  integer,intent(in) :: i1pert
  integer,intent(in) :: i2dir
  integer,intent(in) :: i2pert
  integer,intent(in) :: i3dir
  integer,intent(in) :: i3pert
  integer,intent(in) :: mband
  integer,intent(in) :: mgfft
  integer,intent(in) :: mk1mem
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpert
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpw
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspden
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  real(dp),intent(in) :: cg1(2,mpw*nspinor*mband*mk1mem*nsppol)
  real(dp),intent(in) :: cg3(2,mpw*nspinor*mband*mk1mem*nsppol)
  real(dp),intent(inout) :: d3lo(2,3,mpert,3,mpert,3,mpert)
  integer,intent(in) :: kg(3,mpw*mkmem)
  integer,intent(in) :: npwarr(nkpt)
  real(dp),intent(in) :: occ(mband*nkpt*nsppol)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(inout) :: vtrial1(cplex*nfft,nspden)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: ylm(mpw*mkmem,mpsang*mpsang*psps%useylm)
 end subroutine dfptnl_resp
end interface

interface
 subroutine eig2stern(occ,bdeigrf,clflg,cg1_pert,dim_eig2nkq,dim_eig2rf,eigen0,eigenq,&  
  &  eigen1,eig2nkq,elph2_imagden,esmear,gh0c1_pert,gh1c_pert,ieig2rf,istwfk_pert,&  
  &  mband,mk1mem,mpert,npert,mpi_enreg,mpw1,nkpt_rbz,npwar1,nspinor,nsppol,smdelta,&  
  &  dtset,eigbrd,eigenq_fine,hdr_fine,hdr0) 
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: bdeigrf
  integer,intent(in) :: dim_eig2nkq
  integer,intent(in) :: dim_eig2rf
  integer,intent(in) :: ieig2rf
  integer,intent(in) :: mband
  integer,intent(in) :: mk1mem
  integer,intent(in) :: mpert
  integer,intent(in) :: mpw1
  integer,intent(in) :: nkpt_rbz
  integer,intent(in) :: npert
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: smdelta
  type(dataset_type), intent(in) :: dtset
  real(dp),intent(in) :: elph2_imagden
  real(dp),intent(in) :: esmear
  type(hdr_type),intent(in),optional :: hdr0
  type(hdr_type),intent(in),optional :: hdr_fine
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: cg1_pert(2,mpw1*nspinor*mband*mk1mem*nsppol*dim_eig2rf,3,mpert)
  integer,intent(in) :: clflg(3,mpert)
  real(dp),intent(out) :: eig2nkq(2,mband*nsppol,nkpt_rbz,3,npert,3,npert*dim_eig2nkq)
  real(dp),intent(out),optional :: eigbrd(2,mband*nsppol,nkpt_rbz,3,npert,3,npert)
  real(dp),intent(inout) :: eigen0(nkpt_rbz*mband*nsppol)
  real(dp),intent(in) :: eigen1(nkpt_rbz*2*nsppol*mband**2,3,mpert)
  real(dp),intent(inout) :: eigenq(nkpt_rbz*mband*nsppol)
  real(dp),intent(in),pointer,optional :: eigenq_fine(:,:,:)
  real(dp),intent(in) :: gh0c1_pert(2,mpw1*nspinor*mband*mk1mem*nsppol*dim_eig2rf,3,mpert)
  real(dp),intent(in) :: gh1c_pert(2,mpw1*nspinor*mband*mk1mem*nsppol*dim_eig2rf,3,mpert)
  integer,intent(in) :: istwfk_pert(nkpt_rbz,3,mpert)
  integer,intent(in) :: npwar1(nkpt_rbz,mpert)
  real(dp), intent(in) :: occ(mband*nkpt_rbz*nsppol)
 end subroutine eig2stern
end interface

interface
 subroutine eig2tot(dtfil,xred,psps,pawtab,natom,bdeigrf,clflg,dim_eig2nkq,eigen0,eigenq,eigen1,eig2nkq,&  
  &  elph2_imagden,esmear,ieig2rf,mband,mpert,npert,mpi_enreg,doccde,&  
  &  nkpt_rbz,nsppol,smdelta,rprimd,dtset,occ_rbz,hdr0,eigbrd,eigenq_fine,hdr_fine) 
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use m_pawtab
  implicit none
  integer,intent(in) :: bdeigrf
  integer,intent(in) :: dim_eig2nkq
  integer,intent(in) :: ieig2rf
  integer,intent(in) :: mband
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt_rbz
  integer,intent(in) :: npert
  integer,intent(in) :: nsppol
  integer,intent(in) :: smdelta
  type(datafiles_type), intent(in) :: dtfil
  type(dataset_type), intent(in) :: dtset
  real(dp),intent(in) :: elph2_imagden
  real(dp),intent(in) :: esmear
  type(hdr_type),intent(in) :: hdr0
  type(hdr_type),intent(in),optional :: hdr_fine
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type), intent(inout) :: psps
  integer,intent(in) :: clflg(3,mpert)
  real(dp),intent(in) :: doccde(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp),intent(inout) :: eig2nkq(2,mband*nsppol,nkpt_rbz,3,npert,3,npert*dim_eig2nkq)
  real(dp),intent(inout),optional :: eigbrd(2,mband*nsppol,nkpt_rbz,3,npert,3,npert)
  real(dp),intent(in) :: eigen0(nkpt_rbz*mband*nsppol)
  real(dp),intent(in) :: eigen1(nkpt_rbz*2*nsppol*mband**2,3,mpert)
  real(dp),intent(in) :: eigenq(nkpt_rbz*mband*nsppol)
  real(dp),intent(in),pointer,optional :: eigenq_fine(:,:,:)
  real(dp),intent(in) :: occ_rbz(mband*nkpt_rbz*nsppol)
  type(pawtab_type), intent(inout) :: pawtab(psps%ntypat*psps%usepaw)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine eig2tot
end interface

interface
 subroutine eigen_meandege(eigen0,eigenresp,eigenresp_mean,mband,nband,nkpt,nsppol,option)
  use defs_basis
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  integer,intent(in) :: option
  real(dp),intent(in) :: eigen0(mband*nkpt*nsppol)
  real(dp),intent(in) :: eigenresp((3-option)*mband**(3-option)*nkpt*nsppol)
  real(dp),intent(out) :: eigenresp_mean(mband*nkpt*nsppol)
  integer,intent(in) :: nband(nkpt*nsppol)
 end subroutine eigen_meandege
end interface

interface
 subroutine elph2_fanddw(dim_eig2nkq,displ,eig2nkq,eigen_corr,gprimd,mband,natom,nkpt,nsppol,option,phfrq,prtvol)
  use defs_basis
  implicit none
  integer,intent(in) :: dim_eig2nkq
  integer,intent(in) :: mband
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  integer,intent(in) :: option
  integer,intent(in) :: prtvol
  real(dp),intent(in) :: displ(2*3*natom*3*natom)
  real(dp),intent(in) :: eig2nkq(2,mband*nsppol,nkpt,3,natom,3,natom*dim_eig2nkq)
  real(dp),intent(out) :: eigen_corr(mband*nkpt*nsppol)
  real(dp) :: gprimd(3,3)
  real(dp),intent(in) :: phfrq(3*natom)
 end subroutine elph2_fanddw
end interface

interface
 subroutine elt_ewald(elteew,gmet,gprimd,my_natom,natom,ntypat,rmet,rprimd,&  
  &  typat,ucvol,xred,zion,&  
  &  mpi_atmtab,comm_atom) ! optional arguments (parallelism)
  use defs_basis
  implicit none
  integer,optional,intent(in) :: comm_atom
  integer,intent(in) :: my_natom
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  real(dp),intent(in) :: ucvol
  integer,optional,target,intent(in) :: mpi_atmtab(:)
  real(dp),intent(out) :: elteew(6+3*natom,6)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: zion(ntypat)
 end subroutine elt_ewald
end interface

interface
 subroutine eltxccore(eltfrxc,is2_in,my_natom,natom,nfft,ntypat,&  
  &  n1,n1xccc,n2,n3,rprimd,typat,ucvol,vxc_core,vxc10_core,vxc1is_core,&  
  &  xcccrc,xccc1d,xred,&  
  &  mpi_atmtab,comm_atom) ! optional arguments (parallelism)
  use defs_basis
  implicit none
  integer,optional,intent(in) :: comm_atom
  integer,intent(in) :: is2_in
  integer,intent(in) :: my_natom
  integer,intent(in) :: n1
  integer,intent(in) :: n1xccc
  integer,intent(in) :: n2
  integer,intent(in) :: n3
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: ntypat
  real(dp),intent(in) :: ucvol
  integer,optional,target,intent(in) :: mpi_atmtab(:)
  real(dp),intent(inout) :: eltfrxc(6+3*natom,6)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: vxc10_core(nfft)
  real(dp),intent(in) :: vxc1is_core(nfft)
  real(dp),intent(in) :: vxc_core(nfft)
  real(dp),intent(in) :: xccc1d(n1xccc,6,ntypat)
  real(dp),intent(in) :: xcccrc(ntypat)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine eltxccore
end interface

interface
 subroutine gaugetransfo(cg_k,cwavef,cwavef_d,eig_k,eig1_k,iband,nband_k,&  
  &  mband,npw_k,npw1_k,nspinor,nsppol,occ_k)
  use defs_basis
  implicit none
  integer,intent(in) :: iband
  integer,intent(in) :: mband
  integer,intent(in) :: nband_k
  integer,intent(in) :: npw1_k
  integer,intent(in) :: npw_k
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  real(dp),intent(in) :: cg_k(2,npw_k*nspinor*nband_k)
  real(dp),intent(in) :: cwavef(2,npw1_k*nspinor)
  real(dp),intent(out) :: cwavef_d(2,npw1_k*nspinor)
  real(dp),intent(in) :: eig1_k(2*nsppol*mband**2)
  real(dp),intent(in) :: eig_k(mband*nsppol)
  real(dp),intent(in) :: occ_k(nband_k)
 end subroutine gaugetransfo
end interface

interface
 subroutine getshell(gmet,kneigh,kg_neigh,kptindex,kptopt,kptrlatt,kpt2,&  
  &  kpt3,mkmem,mkmem_max,mpi_enreg,mvwtk,&  
  &  nkpt2,nkpt3,nneigh,nshiftk,rmet,rprimd,shiftk,wtk2)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: kptopt
  integer,intent(in) :: mkmem
  integer,intent(out) :: mkmem_max
  integer,intent(in) :: nkpt2
  integer,intent(in) :: nkpt3
  integer,intent(out) :: nneigh
  integer,intent(inout) :: nshiftk
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(inout) :: kptrlatt(3,3)
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(out) :: kg_neigh(30,nkpt2,3)
  integer,intent(out) :: kneigh(30,nkpt2)
  real(dp),intent(in) :: kpt2(3,nkpt2)
  real(dp),intent(out) :: kpt3(3,nkpt3)
  integer,intent(out) :: kptindex(2,nkpt3)
  real(dp),intent(out) :: mvwtk(30,nkpt2)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: shiftk(3,nshiftk)
  real(dp),intent(in) :: wtk2(nkpt2)
 end subroutine getshell
end interface

interface
 subroutine hartrestr(gmet,gprimd,gsqcut,idir,ipert,mpi_enreg,natom,nfft,ngfft,&  
  &  paral_kgb,rhog,vhartr1)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: idir
  integer,intent(in) :: ipert
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: paral_kgb
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(in) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(out) :: vhartr1(nfft)
 end subroutine hartrestr
end interface

interface
 subroutine inpgkk(bantot1,eigen1,filegkk,hdr1)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: bantot1
  character(len=fnlen),intent(in) :: filegkk
  type(hdr_type), intent(out) :: hdr1
  real(dp),intent(out) :: eigen1(bantot1)
 end subroutine inpgkk
end interface

interface
 subroutine kptfine_av(center,qptrlatt,kpt_fine,nkpt_fine,kpt_fine_sub,&  
  &  nkpt_sub,wgt_sub)
  use defs_basis
  implicit none
  integer,intent(in) :: nkpt_fine
  integer,intent(out) :: nkpt_sub
  integer,pointer :: kpt_fine_sub(:)
  integer,intent(in) :: qptrlatt(3,3)
  real(dp),intent(in) :: center(3)
  real(dp),intent(in) :: kpt_fine(3,nkpt_fine)
  real(dp),pointer :: wgt_sub(:)
 end subroutine kptfine_av
end interface

interface
 subroutine newfermie1(cplex,fermie1,fe1fixed,ipert,istep,ixc,my_natom,natom,nfft,nfftot,&  
  &  nhatfermi,nspden,ntypat,occopt,paw_an,paw_an1,paw_ij1,pawang,pawnzlm,pawrad,&  
  &  pawrhoij1,pawrhoijfermi,pawtab,pawxcdev,prtvol,rhorfermi,&  
  &  ucvol,usepaw,usexcnhat,vtrial1,vxc1,xclevel,&  
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
  integer,intent(in) :: cplex
  integer,intent(in) :: ipert
  integer,intent(in) :: istep
  integer,intent(in) :: ixc
  integer,intent(in) :: my_natom
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nfftot
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  integer,intent(in) :: occopt
  integer,intent(in) :: pawnzlm
  integer,intent(in) :: pawxcdev
  integer,intent(in) :: prtvol
  integer,intent(in) :: usepaw
  integer,intent(in) :: usexcnhat
  integer,intent(in) :: xclevel
  real(dp),intent(in) :: fe1fixed
  real(dp),intent(inout) :: fermie1
  type(pawang_type),intent(in) :: pawang
  real(dp),intent(in) :: ucvol
  integer,optional,target,intent(in) :: mpi_atmtab(:)
  real(dp),intent(in) :: nhatfermi(:,:)
  type(paw_an_type),intent(in) :: paw_an(my_natom*usepaw)
  type(paw_an_type),intent(inout) :: paw_an1(my_natom*usepaw)
  type(paw_ij_type),intent(inout) :: paw_ij1(my_natom*usepaw)
  type(pawrad_type),intent(in) :: pawrad(ntypat*usepaw)
  type(pawrhoij_type),intent(in) :: pawrhoij1(my_natom*usepaw)
  type(pawrhoij_type),intent(in) :: pawrhoijfermi(my_natom*usepaw)
  type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)
  real(dp),intent(in) :: rhorfermi(nfft,nspden)
  real(dp),intent(in) :: vtrial1(cplex*nfft,nspden)
  real(dp),intent(in) :: vxc1(:,:)
 end subroutine newfermie1
end interface

interface
 subroutine outbsd(bdeigrf,dtset,eig2nkq,mpert,nkpt_rbz,unitout)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: bdeigrf
  integer,intent(in) :: mpert
  integer,intent(in) :: nkpt_rbz
  integer,intent(in) :: unitout
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: eig2nkq(2,dtset%mband*dtset%nsppol,nkpt_rbz,3,mpert,3,mpert)
 end subroutine outbsd
end interface

interface
 subroutine outddbnc (dtfil, dtset, hdr0, psps, natom, mpert, rprimd, xred, qpt, d2matr, blkflg)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  type(datafiles_type), intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  type(hdr_type),intent(in) :: hdr0
  type(pseudopotential_type),intent(in) :: psps
  integer,intent(in) :: blkflg(3,mpert,3,mpert)
  real(dp),intent(in) :: d2matr(2,3,mpert,3,mpert)
  real(dp),intent(in) :: qpt(3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine outddbnc
end interface

interface
 subroutine outgkk(bantot0,bantot1,outfile,eigen0,eigen1,hdr0,hdr1,mpi_enreg,phasecg)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: bantot0
  integer,intent(in) :: bantot1
  type(hdr_type),intent(inout) :: hdr0
  type(hdr_type),intent(inout) :: hdr1
  type(mpi_type),intent(inout) :: mpi_enreg
  character(len=fnlen),intent(in) :: outfile
  real(dp),intent(in) :: eigen0(bantot0)
  real(dp),intent(in) :: eigen1(2*bantot1)
  real(dp),intent(in) :: phasecg(2,bantot1)
 end subroutine outgkk
end interface

interface
 subroutine qmatrix(cg,dtefield,qmat,mpw,mpw1,mkmem,mband,npwarr,nkpt,nspinor,nsppol,pwindall)
  use defs_basis
  use m_efield
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: mpw1
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  type(efield_type),intent(in) :: dtefield
  real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  integer,intent(in) :: npwarr(nkpt)
  integer,intent(in) :: pwindall(max(mpw,mpw1)*mkmem,8,3)
  real(dp),intent(out) :: qmat(2,dtefield%nband_occ,dtefield%nband_occ,nkpt,2,3)
 end subroutine qmatrix
end interface

interface
 subroutine redgr (frin,frredgr,mpi_enreg,nfft,ngfft,paral_kgb)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: nfft
  integer,intent(in) :: paral_kgb
  type(mpi_type) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: frin(nfft)
  real(dp),intent(out) :: frredgr(nfft,3)
 end subroutine redgr
end interface

interface
 subroutine rotate_rho(cplex, itirev, mpi_enreg, nfft, ngfft, nspden,&  
  &  rhor1, rhog1_eq, rhor1_eq, symrel1, tnon)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: itirev
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  type(mpi_type),intent(in) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  integer, intent(in) :: symrel1(3,3)
  real(dp),intent(out) :: rhog1_eq(2,nfft)
  real(dp),intent(inout) :: rhor1(cplex*nfft,nspden)
  real(dp),intent(out) :: rhor1_eq(cplex*nfft,nspden)
  real(dp),intent(in) :: tnon(3)
 end subroutine rotate_rho
end interface

interface
 subroutine smeared_delta(eigen0,eigenq,esmear,mband,smdelta,smdfunc)
  use defs_basis
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: smdelta
  real(dp),intent(in) :: esmear
  real(dp),intent(in) :: eigen0(mband)
  real(dp),intent(in) :: eigenq(mband)
  real(dp),intent(out) :: smdfunc(mband,mband)
 end subroutine smeared_delta
end interface

interface
 subroutine sytens(indsym,mpert,natom,nsym,rfpert,symrec,symrel)
  implicit none
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  integer,intent(in) :: nsym
  integer,intent(in) :: indsym(4,nsym,natom)
  integer,intent(inout) :: rfpert(3,mpert,3,mpert,3,mpert)
  integer,intent(in) :: symrec(3,3,nsym)
  integer,intent(in) :: symrel(3,3,nsym)
 end subroutine sytens
end interface

interface
 subroutine vlocalstr(gmet,gprimd,gsqcut,istr,mgfft,mpi_enreg,&  
  &  mqgrid,natom,nattyp,nfft,ngfft,ntypat,paral_kgb,ph1d,qgrid,&  
  &  ucvol,vlspl,vpsp1)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: istr
  integer,intent(in) :: mgfft
  integer,intent(in) :: mqgrid
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: ntypat
  integer,intent(in) :: paral_kgb
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: nattyp(ntypat)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(in) :: vlspl(mqgrid,2,ntypat)
  real(dp),intent(out) :: vpsp1(nfft)
 end subroutine vlocalstr
end interface

interface
 subroutine wrtloctens(blkflg,d2bbb,d2nl,mband,mpert,natom,prtbbb,rprimd,usepaw)
  use defs_basis
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  integer,intent(in) :: prtbbb
  integer,intent(in) :: usepaw
  integer,intent(in) :: blkflg(3,mpert,3,mpert)
  real(dp),intent(inout) :: d2bbb(2,3,3,mpert,mband,mband*prtbbb)
  real(dp),intent(inout) :: d2nl(2,3,mpert,3,mpert)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine wrtloctens
end interface

end module interfaces_72_response
!!***
