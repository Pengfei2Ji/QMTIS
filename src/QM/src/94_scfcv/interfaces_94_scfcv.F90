!!****m* ABINIT/interfaces_94_scfcv
!! NAME
!! interfaces_94_scfcv
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/94_scfcv
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

module interfaces_94_scfcv

 implicit none

interface
 subroutine afterscfloop(atindx,atindx1,cg,computed_forces,cprj,cpus,&  
  &  deltae,diffor,dtefield,dtfil,dtset,eigen,electronpositron,elfr,&  
  &  energies,etotal,favg,fcart,forold,fred,gresid,grewtn,grhf,grhor,grvdw,&  
  &  grxc,gsqcut,hdr,indsym,irrzon,istep,kg,kxc,lrhor,maxfor,mcg,mcprj,mgfftf,&  
  &  moved_atm_inside,mpi_enreg,my_natom,n3xccc,nattyp,nfftf,ngfft,ngfftf,ngrvdw,nhat,&  
  &  nkxc,npwarr,nvresid,occ,optres,paw_an,paw_ij,pawang,pawfgr,&  
  &  pawfgrtab,pawrhoij,pawtab,pel,pel_cg,ph1d,ph1df,phnons,pion,prtfor,prtxml,&  
  &  psps,pwind,pwind_alloc,pwnsfac,res2,resid,residm,results_gs,&  
  &  rhog,rhor,rprimd,stress_needed,strsxc,strten,symrec,synlgr,taug,&  
  &  taur,tollist,usecprj,vhartr,vpsp,vtrial,vxc,vxcavg,wvl,&  
  &  xccc3d,xred,ylm,ylmgr,qvpotzero,conv_retcode)
  use defs_wvltypes
  use m_results_gs
  use m_pawtab
  use m_paw_ij
  use m_pawang
  use m_pawrhoij
  use m_energies
  use defs_abitypes
  use m_pawcprj
  use m_pawfgrtab
  use m_paw_an
  use m_pawfgr
  use defs_basis
  use m_efield
  use defs_datatypes
  use m_electronpositron
  implicit none
  integer,intent(inout) :: computed_forces
  integer,intent(out) :: conv_retcode
  integer,intent(in) :: istep
  integer,intent(in) :: mcg
  integer,intent(in) :: mcprj
  integer,intent(in) :: mgfftf
  integer,intent(in) :: moved_atm_inside
  integer,intent(in) :: my_natom
  integer,intent(in) :: n3xccc
  integer,intent(in) :: nfftf
  integer,intent(in) :: ngrvdw
  integer,intent(in) :: nkxc
  integer,intent(in) :: optres
  integer,intent(in) :: prtfor
  integer,intent(in) :: prtxml
  integer,intent(in) :: pwind_alloc
  integer,intent(in) :: stress_needed
  integer,intent(in) :: usecprj
  real(dp),intent(in) :: cpus
  real(dp),intent(in) :: deltae
  real(dp),intent(inout) :: diffor
  type(efield_type),intent(inout) :: dtefield
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(inout) :: dtset
  type(electronpositron_type),pointer :: electronpositron
  type(energies_type),intent(inout) :: energies
  real(dp),intent(inout) :: etotal
  real(dp),intent(in) :: gsqcut
  type(hdr_type),intent(inout) :: hdr
  real(dp),intent(inout) :: maxfor
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pawfgr_type),intent(in) :: pawfgr
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: qvpotzero
  real(dp),intent(in) :: res2
  real(dp),intent(in) :: residm
  type(results_gs_type),intent(inout) :: results_gs
  real(dp),intent(inout) :: vxcavg
  type(wvl_data),intent(inout) :: wvl
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: ngfftf(18)
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(in) :: atindx1(dtset%natom)
  real(dp),intent(inout) :: cg(2,mcg)
  type(pawcprj_type),intent(in) :: cprj(dtset%natom,mcprj*usecprj)
  real(dp),intent(inout) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp),pointer :: elfr(:,:)
  real(dp),intent(inout) :: favg(3)
  real(dp),intent(inout) :: fcart(3,dtset%natom)
  real(dp),intent(inout) :: forold(3,dtset%natom)
  real(dp),intent(inout) :: fred(3,dtset%natom)
  real(dp),intent(inout) :: gresid(3,dtset%natom)
  real(dp),intent(in) :: grewtn(3,dtset%natom)
  real(dp),intent(inout) :: grhf(3,dtset%natom)
  real(dp),pointer :: grhor(:,:,:)
  real(dp),intent(in) :: grvdw(3,ngrvdw)
  real(dp),intent(inout) :: grxc(3,dtset%natom)
  integer,intent(in) :: indsym(4,dtset%nsym,dtset%natom)
  integer,intent(in) :: irrzon(dtset%nfft**(1-1/dtset%nsym),2,(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  integer,intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
  real(dp),intent(inout) :: kxc(nfftf,nkxc)
  real(dp),pointer :: lrhor(:,:)
  integer,intent(in) :: nattyp(dtset%ntypat)
  real(dp),intent(inout) :: nhat(nfftf,dtset%nspden*psps%usepaw)
  integer,intent(in) :: npwarr(dtset%nkpt)
  real(dp),intent(inout) :: nvresid(nfftf,dtset%nspden)
  real(dp),intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  type(paw_an_type),intent(inout) :: paw_an(my_natom*psps%usepaw)
  type(paw_ij_type),intent(inout) :: paw_ij(my_natom*psps%usepaw)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(my_natom*psps%usepaw)
  type(pawrhoij_type),intent(inout) :: pawrhoij(my_natom*psps%usepaw)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*psps%usepaw)
  real(dp),intent(inout) :: pel(3)
  real(dp),intent(inout) :: pel_cg(3)
  real(dp),intent(inout) :: ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom)
  real(dp),intent(inout) :: ph1df(2,3*(2*mgfftf+1)*dtset%natom)
  real(dp),intent(in) :: phnons(2,dtset%nfft**(1-1/dtset%nsym), &
  &         (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  real(dp),intent(inout) :: pion(3)
  integer,intent(in) :: pwind(pwind_alloc,2,3)
  real(dp),intent(in) :: pwnsfac(2,pwind_alloc)
  real(dp),intent(in) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp),intent(inout) :: rhog(2,nfftf)
  real(dp),intent(inout) :: rhor(nfftf,dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(inout) :: strsxc(6)
  real(dp),intent(inout) :: strten(6)
  integer,intent(in) :: symrec(3,3,dtset%nsym)
  real(dp),intent(inout) :: synlgr(3,dtset%natom)
  real(dp),pointer :: taug(:,:)
  real(dp),pointer :: taur(:,:)
  real(dp),intent(in) :: tollist(12)
  real(dp),intent(inout) :: vhartr(nfftf)
  real(dp),intent(in) :: vpsp(nfftf)
  real(dp),intent(inout) :: vtrial(nfftf,dtset%nspden)
  real(dp),intent(inout) :: vxc(nfftf,dtset%nspden)
  real(dp),intent(inout) :: xccc3d(n3xccc)
  real(dp),intent(inout) :: xred(3,dtset%natom)
  real(dp),intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
  real(dp),intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 end subroutine afterscfloop
end interface

interface
 subroutine elpolariz(atindx1,cg,cprj,dtefield,dtfil,dtset,etotal,enefield,gprimd,hdr,&  
  &  kg,mband,mcg,mcprj,mkmem,mpi_enreg,mpw,my_natom,natom,nattyp,nkpt,&  
  &  npwarr,nsppol,ntypat,pawrhoij,pawtab,&  
  &  pel,pel_cg,pelev,pion,psps,pwind,pwind_alloc,&  
  &  pwnsfac,rprimd,ucvol,usecprj,xred)
  use m_pawrhoij
  use defs_abitypes
  use m_pawcprj
  use m_pawtab
  use defs_basis
  use m_efield
  use defs_datatypes
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: mcg
  integer,intent(in) :: mcprj
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: my_natom
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  integer,intent(in) :: ntypat
  integer,intent(in) :: pwind_alloc
  integer,intent(in) :: usecprj
  type(efield_type),intent(inout) :: dtefield
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(inout) :: dtset
  real(dp),intent(inout) :: enefield
  real(dp),intent(inout) :: etotal
  type(hdr_type),intent(inout) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: atindx1(natom)
  real(dp),intent(in) :: cg(2,mcg)
  type(pawcprj_type),intent(in) :: cprj(natom,mcprj*usecprj)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: kg(3,mpw*mkmem)
  integer,intent(in) :: nattyp(ntypat)
  integer,intent(in) :: npwarr(nkpt)
  type(pawrhoij_type), intent(in) :: pawrhoij(my_natom*psps%usepaw)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*psps%usepaw)
  real(dp),intent(inout) :: pel(3)
  real(dp),intent(in) :: pel_cg(3)
  real(dp),intent(inout) :: pelev(3)
  real(dp),intent(inout) :: pion(3)
  integer,intent(in) :: pwind(pwind_alloc,2,3)
  real(dp),intent(in) :: pwnsfac(2,pwind_alloc)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(inout) :: xred(3,natom)
 end subroutine elpolariz
end interface

interface
 subroutine outscfcv(atindx1,cg,compch_fft,compch_sph,cprj,dimcprj,dmatpawu,dtfil,dtset,&  
  &  ecut,eigen,electronpositron,elfr,etotal,fermie,gmet,gprimd,grhor,hdr,kg,&  
  &  lrhor,mband,mcg,mcprj,mgfftc,mkmem,mpi_enreg,mpsang,mpw,my_natom,natom,&  
  &  nattyp,nfft,ngfft,nhat,nkpt,npwarr,nspden,nsppol,nsym,ntypat,n3xccc,occ,&  
  &  paw_dmft,pawang,pawfgr,pawfgrtab,pawrad,pawrhoij,pawtab,paw_an,paw_ij,&  
  &  prtvol,psps,results_gs,rhor,rprimd,&  
  &  taur,ucvol,usecprj,vhartr,vpsp,vtrial,vxc,wvl_den,xccc3d,xred)
  use m_pawtab
  use defs_wvltypes
  use m_pawrhoij
  use m_pawrad
  use m_paw_ij
  use m_pawang
  use m_paw_dmft
  use defs_abitypes
  use m_pawcprj
  use m_pawfgrtab
  use m_paw_an
  use m_pawfgr
  use defs_basis
  use m_results_gs
  use defs_datatypes
  use m_electronpositron
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: mcg
  integer,intent(in) :: mcprj
  integer,intent(in) :: mgfftc
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpw
  integer,intent(in) :: my_natom
  integer,intent(in) :: n3xccc
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspden
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(in) :: prtvol
  integer,intent(in) :: usecprj
  real(dp),intent(in) :: compch_fft
  real(dp),intent(in) :: compch_sph
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: ecut
  type(electronpositron_type),pointer :: electronpositron
  real(dp),intent(inout) :: etotal
  real(dp),intent(in) :: fermie
  type(hdr_type),intent(inout) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  type(paw_dmft_type), intent(inout) :: paw_dmft
  type(pawang_type),intent(in) :: pawang
  type(pawfgr_type),intent(in) :: pawfgr
  type(pseudopotential_type),intent(in) :: psps
  type(results_gs_type),intent(in) :: results_gs
  real(dp),intent(in) :: ucvol
  type(wvl_denspot_type), intent(in) :: wvl_den
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: atindx1(natom)
  real(dp),intent(inout) :: cg(2,mcg)
  type(pawcprj_type),intent(inout) :: cprj(natom,mcprj*usecprj)
  integer,intent(in) :: dimcprj(natom*usecprj)
  real(dp),intent(in) :: dmatpawu(:,:,:,:)
  real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
  real(dp),pointer :: elfr(:,:)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),pointer :: grhor(:,:,:)
  integer,intent(in) :: kg(3,mpw*mkmem)
  real(dp),pointer :: lrhor(:,:)
  integer,intent(in) :: nattyp(ntypat)
  real(dp),intent(inout) :: nhat(nfft,nspden*psps%usepaw)
  integer,intent(in) :: npwarr(nkpt)
  real(dp),intent(in) :: occ(mband*nkpt*nsppol)
  type(paw_an_type),intent(inout) :: paw_an(my_natom*psps%usepaw)
  type(paw_ij_type),intent(inout) :: paw_ij(my_natom*psps%usepaw)
  type(pawfgrtab_type),intent(in) :: pawfgrtab(my_natom*psps%usepaw)
  type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
  type(pawrhoij_type),target,intent(inout) :: pawrhoij(my_natom*psps%usepaw)
  type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)
  real(dp),intent(inout) :: rhor(nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),pointer :: taur(:,:)
  real(dp),intent(in) :: vhartr(nfft)
  real(dp),intent(in) :: vpsp(nfft)
  real(dp),intent(inout) :: vtrial(nfft,nspden)
  real(dp),intent(inout) :: vxc(nfft,nspden)
  real(dp),intent(in) :: xccc3d(n3xccc)
  real(dp),intent(inout) :: xred(3,natom)
 end subroutine outscfcv
end interface

interface
 subroutine scfcv(atindx,atindx1,cg,cpus,dmatpawu,dtefield,dtfil,dtpawuj,&  
  &  dtset,ecore,eigen,electronpositron,fatvshift,hdr,indsym,&  
  &  initialized,irrzon,kg,mcg,mpi_enreg,my_natom,nattyp,ndtpawuj,nfftf,npwarr,occ,&  
  &  paw_dmft,pawang,pawfgr,pawrad,pawrhoij,pawtab,phnons,psps,pwind,&  
  &  pwind_alloc,pwnsfac,rec_set,resid,results_gs,rhog,rhor,rprimd,&  
  &  scf_history,symrec,taug,taur,wffnew,wvl,xred,xred_old,ylm,ylmgr,conv_retcode)
  use m_pawtab
  use defs_wvltypes
  use defs_rectypes
  use m_results_gs
  use m_pawrad
  use m_pawang
  use m_paw_dmft
  use m_scf_history
  use defs_abitypes
  use m_pawrhoij
  use m_pawfgr
  use defs_basis
  use m_efield
  use defs_datatypes
  use m_electronpositron
  use m_wffile
  implicit none
  integer,intent(out) :: conv_retcode
  integer,intent(inout) :: initialized
  integer,intent(in) :: mcg
  integer,intent(in) :: my_natom
  integer,intent(in) :: ndtpawuj
  integer,intent(inout) :: nfftf
  integer,intent(in) :: pwind_alloc
  real(dp),intent(in) :: cpus
  type(efield_type),intent(inout) :: dtefield
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(inout) :: dtset
  real(dp),intent(in) :: ecore
  type(electronpositron_type),pointer :: electronpositron
  real(dp),intent(in) :: fatvshift
  type(hdr_type),intent(inout) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  type(paw_dmft_type), intent(inout) :: paw_dmft
  type(pawang_type),intent(in) :: pawang
  type(pawfgr_type),intent(inout) :: pawfgr
  type(pseudopotential_type),intent(in) :: psps
  type(recursion_type),intent(inout) :: rec_set
  type(results_gs_type),intent(inout) :: results_gs
  type(scf_history_type),intent(inout) :: scf_history
  type(wffile_type),intent(inout) :: wffnew
  type(wvl_data),intent(inout) :: wvl
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(in) :: atindx1(dtset%natom)
  real(dp), intent(inout) :: cg(2,mcg)
  real(dp), intent(inout) :: dmatpawu(:,:,:,:)
  type(macro_uj_type),intent(inout) :: dtpawuj(0:ndtpawuj)
  real(dp), intent(inout) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
  integer,intent(in) :: indsym(4,dtset%nsym,dtset%natom)
  integer, intent(in) :: irrzon(dtset%nfft**(1-1/dtset%nsym),2, &
  &         (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
  integer, intent(in) :: nattyp(psps%ntypat)
  integer, intent(in) :: npwarr(dtset%nkpt)
  real(dp), intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  type(pawrad_type), intent(in) :: pawrad(psps%ntypat*psps%usepaw)
  type(pawrhoij_type), intent(inout) :: pawrhoij(my_natom*psps%usepaw)
  type(pawtab_type), intent(in) :: pawtab(psps%ntypat*psps%usepaw)
  real(dp), intent(in) :: phnons(2,dtset%nfft**(1-1/dtset%nsym), &
  &         (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  integer, intent(in) :: pwind(pwind_alloc,2,3)
  real(dp), intent(in) :: pwnsfac(2,pwind_alloc)
  real(dp), intent(inout) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp), pointer :: rhog(:,:)
  real(dp), pointer :: rhor(:,:)
  real(dp), intent(in) :: rprimd(3,3)
  integer, intent(in) :: symrec(3,3,dtset%nsym)
  real(dp), pointer :: taug(:,:)
  real(dp), pointer :: taur(:,:)
  real(dp), intent(inout) :: xred(3,dtset%natom)
  real(dp), intent(inout) :: xred_old(3,dtset%natom)
  real(dp), intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
  real(dp), intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 end subroutine scfcv
end interface

end module interfaces_94_scfcv
!!***
