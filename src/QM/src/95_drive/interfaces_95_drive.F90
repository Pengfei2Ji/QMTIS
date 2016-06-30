!!****m* ABINIT/interfaces_95_drive
!! NAME
!! interfaces_95_drive
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/95_drive
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

module interfaces_95_drive

 implicit none

interface
 subroutine bethe_salpeter(acell,codvsn,Dtfil,Dtset,Pawang,Pawrad,Pawtab,Psps,rprim,xred)
  use m_pawrad
  use m_pawang
  use defs_abitypes
  use m_pawtab
  use defs_basis
  use defs_datatypes
  implicit none
  type(datafiles_type),intent(inout) :: Dtfil
  type(dataset_type),intent(inout) :: Dtset
  type(pawang_type),intent(inout) :: Pawang
  type(pseudopotential_type),intent(inout) :: Psps
  character(len=6),intent(in) :: codvsn
  type(pawrad_type),intent(inout) :: Pawrad(Psps%ntypat*Psps%usepaw)
  type(pawtab_type),intent(inout) :: Pawtab(Psps%ntypat*Psps%usepaw)
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(in) :: xred(3,Dtset%natom)
 end subroutine bethe_salpeter
end interface

interface
 subroutine dfpt_looppert(atindx,blkflg,codvsn,cpus,dim_eigbrd,dim_eig2nkq,doccde,&  
  &  ddkfil,dtfil,dtset,dyew,dyfrlo,dyfrnl,dyfrx1,dyfrx2,dyvdw,&  
  &  dyfr_cplex,dyfr_nondiag,d2bbb,d2lo,d2nl,d2ovl,efmasdeg,efmasfr,eigbrd,eig2nkq,&  
  &  eltcore,elteew,eltfrhar,eltfrkin,eltfrloc,eltfrnl,eltfrxc,eltvdw,&  
  &  etotal,fermie,iexit,indsym,kxc,&  
  &  mkmem,mkqmem,mk1mem,mpert,mpi_enreg,my_natom,nattyp,&  
  &  nfftf,nhat,nkpt,nkxc,nspden,nsym,occ,&  
  &  paw_an,paw_ij,pawang,pawfgr,pawfgrtab,pawrad,pawrhoij,pawtab,&  
  &  pertsy,prtbbb,psps,rfpert,rhog,rhor,symq,symrec,timrev,&  
  &  usecprj,usevdw,vtrial,vxc,vxcavg,xred,clflg,occ_rbz_pert,eigen0_pert,eigenq_pert,&  
  &  eigen1_pert,nkpt_rbz,eigenq_fine,hdr_fine,hdr0)
  use m_pawtab
  use m_pawrad
  use m_paw_ij
  use m_pawang
  use m_pawrhoij
  use defs_abitypes
  use m_efmas_defs
  use m_paw_an
  use m_pawfgr
  use defs_basis
  use defs_datatypes
  use m_pawfgrtab
  implicit none
  integer, intent(in) :: dim_eig2nkq
  integer, intent(in) :: dim_eigbrd
  integer, intent(in) :: dyfr_cplex
  integer, intent(in) :: dyfr_nondiag
  integer, intent(out) :: iexit
  integer, intent(in) :: mk1mem
  integer, intent(in) :: mkmem
  integer, intent(in) :: mkqmem
  integer, intent(in) :: mpert
  integer, intent(inout) :: my_natom
  integer, intent(in) :: nfftf
  integer, intent(in) :: nkpt
  integer, intent(out) :: nkpt_rbz
  integer, intent(in) :: nkxc
  integer, intent(in) :: nspden
  integer, intent(in) :: nsym
  integer, intent(in) :: prtbbb
  integer, intent(in) :: timrev
  integer, intent(in) :: usecprj
  integer, intent(in) :: usevdw
  character(len=6), intent(in) :: codvsn
  real(dp), intent(in) :: cpus
  type(datafiles_type), intent(in) :: dtfil
  type(dataset_type), intent(in), target :: dtset
  real(dp), intent(inout) :: etotal
  real(dp), intent(inout) :: fermie
  type(hdr_type),intent(out) :: hdr0
  type(hdr_type),intent(out) :: hdr_fine
  type(mpi_type), intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pawfgr_type),intent(in) :: pawfgr
  type(pseudopotential_type), intent(inout) :: psps
  real(dp), intent(in) :: vxcavg
  integer, intent(out) :: ddkfil(3)
  integer, intent(in) :: atindx(dtset%natom)
  integer, intent(inout) :: blkflg(3,mpert,3,mpert)
  integer, intent(out) :: clflg(3,mpert)
  real(dp), intent(inout) :: d2bbb(2,3,3,mpert,dtset%mband,dtset%mband*prtbbb)
  real(dp), intent(inout) :: d2lo(2,3,mpert,3,mpert)
  real(dp), intent(inout) :: d2nl(2,3,mpert,3,mpert)
  real(dp), intent(inout) :: d2ovl(2,3,mpert,3,mpert*psps%usepaw)
  real(dp), intent(in) :: doccde(dtset%mband*nkpt*dtset%nsppol)
  real(dp), intent(in) :: dyew(2,3,dtset%natom,3,dtset%natom)
  real(dp), intent(in) :: dyfrlo(3,3,dtset%natom)
  real(dp), intent(in) :: dyfrnl(dyfr_cplex,3,3,dtset%natom,1+(dtset%natom-1)*dyfr_nondiag)
  real(dp), intent(in) :: dyfrx1(2,3,dtset%natom,3,dtset%natom)
  real(dp), intent(in) :: dyfrx2(3,3,dtset%natom)
  real(dp), intent(in) :: dyvdw(2,3,dtset%natom,3,dtset%natom*usevdw)
  type(efmasdeg_type),allocatable,intent(in) :: efmasdeg(:)
  type(efmasfr_type),allocatable,intent(in) :: efmasfr(:,:)
  real(dp), intent(out) :: eig2nkq(2,dtset%mband*dtset%nsppol,nkpt, &
  &         3,dtset%natom,3,dtset%natom*dim_eig2nkq)
  real(dp), intent(out) :: eigbrd(2,dtset%mband*dtset%nsppol,nkpt, &
  &         3,dtset%natom,3,dtset%natom*dim_eigbrd)
  real(dp),intent(out) :: eigen0_pert(:)
  real(dp),pointer :: eigen1_pert(:,:,:)
  real(dp),pointer :: eigenq_fine(:,:,:)
  real(dp),intent(out) :: eigenq_pert(:)
  real(dp), intent(in) :: eltcore(6,6)
  real(dp), intent(in) :: elteew(6+3*dtset%natom,6)
  real(dp), intent(in) :: eltfrhar(6,6)
  real(dp), intent(in) :: eltfrkin(6,6)
  real(dp), intent(in) :: eltfrloc(6+3*dtset%natom,6)
  real(dp), intent(in) :: eltfrnl(6+3*dtset%natom,6)
  real(dp), intent(in) :: eltfrxc(6+3*dtset%natom,6)
  real(dp), intent(in) :: eltvdw(6+3*dtset%natom,6*usevdw)
  integer, intent(in) :: indsym(4,nsym,dtset%natom)
  real(dp), intent(in) :: kxc(nfftf,nkxc)
  integer, intent(in) :: nattyp(dtset%ntypat)
  real(dp), intent(in) :: nhat(nfftf,nspden)
  real(dp), intent(in) :: occ(dtset%mband*nkpt*dtset%nsppol)
  real(dp),intent(out) :: occ_rbz_pert(:)
  type(paw_an_type),allocatable,target,intent(inout) :: paw_an(:)
  type(paw_ij_type),allocatable,target,intent(inout) :: paw_ij(:)
  type(pawfgrtab_type),allocatable,target,intent(inout) :: pawfgrtab(:)
  type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
  type(pawrhoij_type),allocatable,target,intent(inout) :: pawrhoij(:)
  type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)
  integer, intent(in) :: pertsy(3,mpert)
  integer, intent(in) :: rfpert(mpert)
  real(dp), intent(in) :: rhog(2,nfftf)
  real(dp), intent(in) :: rhor(nfftf,nspden)
  integer, intent(in) :: symq(4,2,nsym)
  integer, intent(in) :: symrec(3,3,nsym)
  real(dp), intent(in) :: vtrial(nfftf,nspden)
  real(dp), intent(in) :: vxc(nfftf,nspden)
  real(dp), intent(inout) :: xred(3,dtset%natom)
 end subroutine dfpt_looppert
end interface

interface
 subroutine dfptnl_loop(blkflg,cg,cgindex,dtfil,dtset,d3lo,&  
  &  etotal,gmet,gprimd,gsqcut,&  
  &  hdr,kg,kneigh,kg_neigh,kptindex,kpt3,kxc,k3xc,mband,mgfft,mkmem,mkmem_max,mk1mem,&  
  &  mpert,mpi_enreg,mpw,mvwtk,natom,nfft,nkpt,nkpt3,nkxc,nk3xc,nneigh,nspinor,nsppol,&  
  &  npwarr,occ,psps,pwind,&  
  &  rfpert,rprimd,ucvol,xred)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: mgfft
  integer,intent(in) :: mk1mem
  integer,intent(in) :: mkmem
  integer,intent(in) :: mkmem_max
  integer,intent(in) :: mpert
  integer,intent(in) :: mpw
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nk3xc
  integer,intent(in) :: nkpt
  integer,intent(in) :: nkpt3
  integer,intent(in) :: nkxc
  integer,intent(in) :: nneigh
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(inout) :: dtset
  real(dp),intent(inout) :: etotal
  real(dp),intent(in) :: gsqcut
  type(hdr_type),intent(inout) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: ucvol
  integer,intent(inout) :: blkflg(3,mpert,3,mpert,3,mpert)
  real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  integer,intent(in) :: cgindex(nkpt,nsppol)
  real(dp),intent(inout) :: d3lo(2,3,mpert,3,mpert,3,mpert)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: k3xc(nfft,nk3xc)
  integer,intent(in) :: kg(3,mk1mem*mpw)
  integer,intent(in) :: kg_neigh(30,nkpt,3)
  integer,intent(in) :: kneigh(30,nkpt)
  real(dp),intent(in) :: kpt3(3,nkpt3)
  integer,intent(in) :: kptindex(2,nkpt3)
  real(dp),intent(in) :: kxc(nfft,nkxc)
  real(dp),intent(in) :: mvwtk(30,nkpt)
  integer,intent(in) :: npwarr(nkpt)
  real(dp),intent(inout) :: occ(mband*nkpt*nsppol)
  integer,intent(in) :: pwind(mpw,nneigh,mkmem)
  integer,intent(in) :: rfpert(3,mpert,3,mpert,3,mpert)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine dfptnl_loop
end interface

interface
 subroutine driver(codvsn,cpui,dtsets,filnam,filstat,&  
  &  mpi_enregs,ndtset,ndtset_alloc,npsp,pspheads,results_out)
  use defs_basis
  use m_results_out
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: ndtset
  integer,intent(in) :: ndtset_alloc
  integer,intent(in) :: npsp
  character(len=6),intent(in) :: codvsn
  real(dp),intent(in) :: cpui
  character(len=fnlen),intent(in) :: filstat
  type(dataset_type),intent(inout) :: dtsets(0:ndtset_alloc)
  character(len=fnlen),intent(in) :: filnam(5)
  type(mpi_type),intent(inout) :: mpi_enregs(0:ndtset_alloc)
  type(pspheader_type),intent(in) :: pspheads(npsp)
  type(results_out_type),target,intent(inout) :: results_out(0:ndtset_alloc)
 end subroutine driver
end interface

interface
 subroutine dtfil_init(dtfil,dtset,filnam,filstat,idtset,jdtset_,mpi_enreg,ndtset,&  
  &  image_index) ! optional argument
  use defs_basis
  use defs_abitypes
  implicit none
  integer, intent(in) :: idtset
  integer, optional, intent(in) :: image_index
  integer, intent(in) :: ndtset
  type(datafiles_type),intent(inout) :: dtfil
  type(dataset_type),intent(in) :: dtset
  character(len=fnlen),intent(in) :: filstat
  type(mpi_type),intent(in) :: mpi_enreg
  character(len=fnlen),intent(in) :: filnam(5)
  integer :: jdtset_(0:ndtset)
 end subroutine dtfil_init
end interface

interface
 subroutine dtfil_init_img(dtfil,dtset,dtsets,idtset,jdtset,ndtset,ndtset_alloc)
  use defs_abitypes
  implicit none
  integer, intent(in) :: idtset
  integer, intent(in) :: ndtset
  integer, intent(in) :: ndtset_alloc
  type(datafiles_type),intent(out) :: dtfil
  type(dataset_type),intent(in) :: dtset
  type(dataset_type),intent(in) :: dtsets(0:ndtset_alloc)
  integer :: jdtset(0:ndtset)
 end subroutine dtfil_init_img
end interface

interface
 subroutine dtfil_init_time(dtfil,iapp)
  use defs_abitypes
  implicit none
  integer, intent(in) :: iapp
  type(datafiles_type),intent(inout) :: dtfil
 end subroutine dtfil_init_time
end interface

interface
 subroutine eph(acell,codvsn,dtfil,dtset,pawang,pawrad,pawtab,psps,rprim,xred)
  use m_pawrad
  use m_pawang
  use defs_abitypes
  use m_pawtab
  use defs_basis
  use defs_datatypes
  implicit none
  character(len=6),intent(in) :: codvsn
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(inout) :: dtset
  type(pawang_type),intent(inout) :: pawang
  type(pseudopotential_type),intent(inout) :: psps
  real(dp),intent(in) :: acell(3)
  type(pawrad_type),intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
  type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(in) :: xred(3,dtset%natom)
 end subroutine eph
end interface

interface
 subroutine gstate(args_gs,acell,codvsn,cpui,dtfil,dtset,iexit,initialized,&  
  &  mpi_enreg,npwtot,occ,pawang,pawrad,pawtab,&  
  &  psps,results_gs,rprim,scf_history,vel,vel_cell,wvl,xred)
  use m_pawrad
  use m_pawang
  use m_scf_history
  use m_args_gs
  use defs_abitypes
  use defs_basis
  use m_pawtab
  use m_results_gs
  use defs_datatypes
  use defs_wvltypes
  implicit none
  integer,intent(inout) :: iexit
  integer,intent(inout) :: initialized
  type(args_gs_type),intent(in) :: args_gs
  character(len=6),intent(in) :: codvsn
  real(dp),intent(in) :: cpui
  type(datafiles_type),intent(inout) :: dtfil
  type(dataset_type),intent(inout) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(inout) :: pawang
  type(pseudopotential_type),intent(inout) :: psps
  type(results_gs_type),intent(inout) :: results_gs
  type(scf_history_type),target,intent(inout) :: scf_history
  type(wvl_data),intent(inout) :: wvl
  real(dp),intent(inout) :: acell(3)
  integer,intent(out) :: npwtot(dtset%nkpt)
  real(dp),intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  type(pawrad_type),intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
  type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)
  real(dp),intent(inout) :: rprim(3,3)
  real(dp),intent(inout) :: vel(3,dtset%natom)
  real(dp),intent(inout) :: vel_cell(3,3)
  real(dp),intent(inout) :: xred(3,dtset%natom)
 end subroutine gstate
end interface

interface
 subroutine gstateimg(acell_img,amu_img,codvsn,cpui,dtfil,dtset,etotal_img,fcart_img,&  
  &  fred_img,iexit,mixalch_img,mpi_enreg,nimage,npwtot,occ_img,&  
  &  pawang,pawrad,pawtab,psps,&  
  &  rprim_img,strten_img,vel_cell_img,vel_img,wvl,xred_img,&  
  &  filnam,filstat,idtset,jdtset,ndtset) ! optional arguments
  use m_pawrad
  use m_pawang
  use defs_abitypes
  use m_pawtab
  use defs_basis
  use defs_datatypes
  use defs_wvltypes
  implicit none
  integer,optional,intent(in) :: idtset
  integer,intent(inout) :: iexit
  integer,optional,intent(in) :: ndtset
  integer,intent(in) :: nimage
  character(len=6),intent(in) :: codvsn
  real(dp),intent(in) :: cpui
  type(datafiles_type),target,intent(inout) :: dtfil
  type(dataset_type),target,intent(inout) :: dtset
  character(len=fnlen),optional,intent(in) :: filstat
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(inout) :: pawang
  type(pseudopotential_type),intent(inout) :: psps
  type(wvl_data),intent(inout) :: wvl
  integer,optional,intent(in) :: jdtset(:)
  real(dp),intent(inout) :: acell_img(3,nimage)
  real(dp),intent(inout) :: amu_img(dtset%ntypat,nimage)
  real(dp), intent(out) :: etotal_img(nimage)
  real(dp), intent(out) :: fcart_img(3,dtset%natom,nimage)
  character(len=fnlen),optional,intent(in) :: filnam(:)
  real(dp), intent(out) :: fred_img(3,dtset%natom,nimage)
  real(dp),intent(inout) :: mixalch_img(dtset%npspalch,dtset%ntypalch,nimage)
  integer,intent(out) :: npwtot(dtset%nkpt)
  real(dp),intent(inout) :: occ_img(dtset%mband*dtset%nkpt*dtset%nsppol,nimage)
  type(pawrad_type),intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
  type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)
  real(dp),intent(inout) :: rprim_img(3,3,nimage)
  real(dp), intent(out) :: strten_img(6,nimage)
  real(dp),intent(inout) :: vel_cell_img(3,3,nimage)
  real(dp),intent(inout) :: vel_img(3,dtset%natom,nimage)
  real(dp),intent(inout) :: xred_img(3,dtset%natom,nimage)
 end subroutine gstateimg
end interface

interface
 subroutine gwls_sternheimer(acell_img,amu_img,codvsn,cpui,dtfil,dtset,etotal_img,fcart_img,&  
  &  fred_img,iexit,mixalch_img,mpi_enreg,nimage,npwtot,occ_img,&  
  &  pawang,pawrad,pawtab,psps,rprim_img,strten_img,vel_cell_img,vel_img,xred_img,&  
  &  filnam,filstat,idtset,jdtset,ndtset) ! optional arguments
  use m_pawrad
  use m_pawang
  use defs_abitypes
  use m_pawtab
  use defs_basis
  use defs_datatypes
  implicit none
  integer,optional,intent(in) :: idtset
  integer,intent(inout) :: iexit
  integer,optional,intent(in) :: ndtset
  integer,intent(in) :: nimage
  character(len=6),intent(in) :: codvsn
  real(dp),intent(in) :: cpui
  type(datafiles_type),target,intent(inout) :: dtfil
  type(dataset_type),target,intent(inout) :: dtset
  character(len=fnlen),optional,intent(in) :: filstat
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(inout) :: pawang
  type(pseudopotential_type),intent(inout) :: psps
  integer,optional,intent(in) :: jdtset(:)
  real(dp),intent(inout) :: acell_img(3,nimage)
  real(dp),intent(inout) :: amu_img(dtset%ntypat,nimage)
  real(dp), intent(out) :: etotal_img(nimage)
  real(dp), intent(out) :: fcart_img(3,dtset%natom,nimage)
  character(len=fnlen),optional,intent(in) :: filnam(:)
  real(dp), intent(out) :: fred_img(3,dtset%natom,nimage)
  real(dp),intent(inout) :: mixalch_img(dtset%npspalch,dtset%ntypalch,nimage)
  integer,intent(out) :: npwtot(dtset%nkpt)
  real(dp),intent(inout) :: occ_img(dtset%mband*dtset%nkpt*dtset%nsppol,nimage)
  type(pawrad_type),intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
  type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)
  real(dp),intent(inout) :: rprim_img(3,3,nimage)
  real(dp), intent(out) :: strten_img(6,nimage)
  real(dp),intent(inout) :: vel_cell_img(3,3,nimage)
  real(dp),intent(inout) :: vel_img(3,dtset%natom,nimage)
  real(dp),intent(inout) :: xred_img(3,dtset%natom,nimage)
 end subroutine gwls_sternheimer
end interface

interface
 subroutine iofn1(filnam,filstat,comm)
  use defs_basis
  implicit none
  integer,intent(in) :: comm
  character(len=fnlen), intent(out) :: filstat
  character(len=fnlen), intent(out) :: filnam(5)
 end subroutine iofn1
end interface

interface
 subroutine mover(scfcv_args,ab_xfh,acell,amass,dtfil,&  
  &  electronpositron,rhog,rhor,rprimd,vel,vel_cell,xred,xred_old)
  use m_scfcv
  use defs_basis
  use m_abimover
  use defs_abitypes
  use m_electronpositron
  implicit none
  type(ab_xfh_type),intent(inout) :: ab_xfh
  type(datafiles_type),intent(inout),target :: dtfil
  type(electronpositron_type),pointer :: electronpositron
  type(scfcv_t),intent(inout) :: scfcv_args
  real(dp),intent(inout) :: acell(3)
  real(dp), intent(in),target :: amass(:)
  real(dp), pointer :: rhog(:,:)
  real(dp), pointer :: rhor(:,:)
  real(dp), intent(inout) :: rprimd(3,3)
  real(dp), intent(inout) :: vel(3,scfcv_args%dtset%natom)
  real(dp), intent(inout) :: vel_cell(3,3)
  real(dp), intent(inout) :: xred(3,scfcv_args%dtset%natom)
  real(dp), intent(inout) :: xred_old(3,scfcv_args%dtset%natom)
 end subroutine mover
end interface

interface
 subroutine nonlinear(codvsn,dtfil,dtset,etotal,iexit,&  
  &  mband,mgfft,mkmem,mpi_enreg,mpw,natom,nfft,nkpt,npwtot,nspden,&  
  &  nspinor,nsppol,nsym,occ,pawrad,pawtab,psps,xred)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use m_pawrad
  use m_pawtab
  implicit none
  integer,intent(in) :: iexit
  integer,intent(in) :: mband
  integer,intent(in) :: mgfft
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspden
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym
  character(len=6),intent(in) :: codvsn
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(inout) :: dtset
  real(dp),intent(inout) :: etotal
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(inout) :: psps
  integer,intent(out) :: npwtot(nkpt)
  real(dp),intent(inout) :: occ(mband*nkpt*nsppol)
  type(pawrad_type),intent(inout) :: pawrad(psps%ntypat,psps%usepaw)
  type(pawtab_type),intent(inout) :: pawtab(psps%ntypat,psps%usepaw)
  real(dp),intent(inout) :: xred(3,natom)
 end subroutine nonlinear
end interface

interface
 subroutine pawuj_drive(scfcv,&  
  &  dtset,electronpositron,rhog,rhor,rprimd,&  
  &  xred,xred_old)
  use m_scfcv
  use defs_basis
  use defs_abitypes
  use m_electronpositron
  implicit none
  type(dataset_type),intent(inout) :: dtset
  type(electronpositron_type),pointer :: electronpositron
  type(scfcv_t), intent(inout) :: scfcv
  real(dp), pointer :: rhog(:,:)
  real(dp), pointer :: rhor(:,:)
  real(dp), intent(inout) :: rprimd(3,3)
  real(dp), intent(inout) :: xred(3,dtset%natom)
  real(dp), intent(inout) :: xred_old(3,dtset%natom)
 end subroutine pawuj_drive
end interface

interface
 subroutine respfn(codvsn,cpui,dtfil,dtset,etotal,iexit,&  
  &  mkmems,mpi_enreg,npwtot,occ,pawang,pawrad,pawtab,psps,results_respfn,xred)
  use m_pawrad
  use m_pawang
  use defs_abitypes
  use m_pawtab
  use defs_basis
  use m_results_respfn
  use defs_datatypes
  implicit none
  integer,intent(inout) :: iexit
  character(len=6),intent(in) :: codvsn
  real(dp),intent(in) :: cpui
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(inout) :: etotal
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(inout) :: pawang
  type(pseudopotential_type),intent(inout) :: psps
  type(results_respfn_type),intent(inout) :: results_respfn
  integer,intent(in) :: mkmems(3)
  integer,intent(inout) :: npwtot(dtset%nkpt)
  real(dp),intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  type(pawrad_type),intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
  type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)
  real(dp),intent(inout) :: xred(3,dtset%natom)
 end subroutine respfn
end interface

interface
 subroutine screening(acell,codvsn,Dtfil,Dtset,Pawang,Pawrad,Pawtab,Psps,rprim)
  use m_pawrad
  use m_pawang
  use defs_abitypes
  use m_pawtab
  use defs_basis
  use defs_datatypes
  implicit none
  type(datafiles_type),intent(in) :: Dtfil
  type(dataset_type),intent(inout) :: Dtset
  type(pawang_type),intent(inout) :: Pawang
  type(pseudopotential_type),intent(inout) :: Psps
  character(len=6),intent(in) :: codvsn
  type(pawrad_type),intent(inout) :: Pawrad(Psps%ntypat*Dtset%usepaw)
  type(pawtab_type),intent(inout) :: Pawtab(Psps%ntypat*Dtset%usepaw)
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: rprim(3,3)
 end subroutine screening
end interface

interface
 subroutine sigma(acell,codvsn,Dtfil,Dtset,Pawang,Pawrad,Pawtab,Psps,rprim,converged)
  use m_pawrad
  use m_pawang
  use defs_abitypes
  use m_pawtab
  use defs_basis
  use defs_datatypes
  implicit none
  type(datafiles_type),intent(in) :: Dtfil
  type(dataset_type),intent(inout) :: Dtset
  type(pawang_type),intent(inout) :: Pawang
  type(pseudopotential_type),intent(inout) :: Psps
  character(len=6),intent(in) :: codvsn
  logical,intent(out) :: converged
  type(pawrad_type),intent(inout) :: Pawrad(Psps%ntypat*Psps%usepaw)
  type(pawtab_type),intent(inout) :: Pawtab(Psps%ntypat*Psps%usepaw)
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: rprim(3,3)
 end subroutine sigma
end interface

interface
 subroutine testfi(builtintest,etotal,filstat,fred,natom,strten,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: builtintest
  integer,intent(in) :: natom
  real(dp),intent(in) :: etotal
  character(len=fnlen),intent(in) :: filstat
  real(dp),intent(in) :: fred(3,natom)
  real(dp),intent(in) :: strten(6)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine testfi
end interface

interface
 subroutine timana(mpi_enreg,natom,nband,ndtset,nfft,nkpt,npwtot,nsppol,timopt)
  use defs_abitypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: ndtset
  integer,intent(in) :: nfft
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  integer,intent(in) :: timopt
  type(mpi_type),intent(in) :: mpi_enreg
  integer,intent(in) :: nband(nkpt*nsppol)
  integer,intent(in) :: npwtot(nkpt)
 end subroutine timana
end interface

interface
 subroutine wfk_analyze(acell,codvsn,dtfil,dtset,pawang,pawrad,pawtab,psps,rprim,xred)
  use m_pawrad
  use m_pawang
  use defs_abitypes
  use m_pawtab
  use defs_basis
  use defs_datatypes
  implicit none
  character(len=6),intent(in) :: codvsn
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(inout) :: dtset
  type(pawang_type),intent(inout) :: pawang
  type(pseudopotential_type),intent(inout) :: psps
  real(dp),intent(in) :: acell(3)
  type(pawrad_type),intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
  type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(in) :: xred(3,dtset%natom)
 end subroutine wfk_analyze
end interface

end module interfaces_95_drive
!!***
