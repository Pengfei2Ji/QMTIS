!!****m* ABINIT/interfaces_70_gw
!! NAME
!! interfaces_70_gw
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/70_gw
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

module interfaces_70_gw

 implicit none

interface
 subroutine calc_rpa_functional(gwrpacorr,iqcalc,iq,Ep,Pvc,Qmesh,Dtfil,gmet,chi0,spaceComm,ec_rpa)
  use m_vcoul
  use m_bz_mesh
  use defs_abitypes
  use m_gwdefs
  use defs_basis
  implicit none
  integer,intent(in) :: gwrpacorr
  integer,intent(in) :: iq
  integer,intent(in) :: iqcalc
  integer,intent(in) :: spaceComm
  type(datafiles_type),intent(in) :: Dtfil
  type(em1params_t),intent(in) :: Ep
  type(vcoul_t),intent(in) :: Pvc
  type(kmesh_t),intent(in) :: Qmesh
  complex(gwpc),intent(inout) :: chi0(Ep%npwe,Ep%npwe,Ep%nomega)
  real(dp),intent(inout) :: ec_rpa(gwrpacorr)
  real(dp),intent(in) :: gmet(3,3)
 end subroutine calc_rpa_functional
end interface

interface
 subroutine calc_sig_ppm_comp(npwc,nomega,rhotwgp,botsq,otq,omegame0i_io,zcut,theta_mu_minus_e0i,ket,ppmodel,npwx,npwc1,npwc2)
  use defs_basis
  implicit none
  integer,intent(in) :: nomega
  integer,intent(in) :: npwc
  integer,intent(in) :: npwc1
  integer,intent(in) :: npwc2
  integer,intent(in) :: npwx
  integer,intent(in) :: ppmodel
  real(dp),intent(in) :: omegame0i_io
  real(dp),intent(in) :: theta_mu_minus_e0i
  real(dp),intent(in) :: zcut
  complex(gwpc),intent(in) :: botsq(npwc,npwc1)
  complex(gwpc),intent(inout) :: ket(npwc,nomega)
  complex(gwpc),intent(in) :: otq(npwc,npwc2)
  complex(gwpc),intent(in) :: rhotwgp(npwx)
 end subroutine calc_sig_ppm_comp
end interface

interface
 subroutine calc_sigc_cd(npwc,npwx,nspinor,nomega,nomegae,nomegaer,nomegaei,rhotwgp,&  
  &  omega,epsm1q,omegame0i,theta_mu_minus_e0i,ket,plasmafreq,npoles_missing,calc_poles,method)
  use defs_basis
  implicit none
  integer, intent(in), optional :: method
  integer,intent(in) :: nomega
  integer,intent(in) :: nomegae
  integer,intent(in) :: nomegaei
  integer,intent(in) :: nomegaer
  integer,intent(inout) :: npoles_missing
  integer,intent(in) :: npwc
  integer,intent(in) :: npwx
  integer,intent(in) :: nspinor
  real(dp),intent(in) :: plasmafreq
  real(dp),intent(in) :: theta_mu_minus_e0i
  logical, intent(in), optional :: calc_poles(nomega)
  complex(gwpc),intent(in) :: epsm1q(npwc,npwc,nomegae)
  complex(gwpc),intent(inout) :: ket(nspinor*npwc,nomega)
  complex(dpc),intent(in) :: omega(nomegae)
  real(dp),intent(in) :: omegame0i(nomega)
  complex(gwpc),intent(in) :: rhotwgp(npwx*nspinor)
 end subroutine calc_sigc_cd
end interface

interface
 subroutine calc_sigc_me(sigmak_ibz,ikcalc,nomega_sigc,minbnd,maxbnd,&  
  &  Dtset,Cryst,QP_BSt,Sigp,Sr,Er,Gsph_Max,Gsph_c,Vcp,Kmesh,Qmesh,Ltg_k,&  
  &  PPm,Pawtab,Pawang,Paw_pwff,Pawfgrtab,Paw_onsite,Psps,Wfd,Wfdf,allQP_sym,&  
  &  gwc_ngfft,rho_ngfft,rho_nfftot,rhor,use_aerhor,aepaw_rhor,sigcme_tmp)
  use m_vcoul
  use m_pawtab
  use m_sigma
  use m_esymm
  use defs_basis
  use m_bz_mesh
  use m_pawang
  use m_pawfgrtab
  use m_wfd
  use m_paw_pwaves_lmn
  use defs_abitypes
  use m_ppmodel
  use m_crystal
  use m_gsphere
  use m_pawpwij
  use defs_datatypes
  use m_gwdefs
  use m_screening
  implicit none
  integer,intent(in) :: ikcalc
  integer,intent(in) :: maxbnd
  integer,intent(in) :: minbnd
  integer,intent(in) :: nomega_sigc
  integer,intent(in) :: rho_nfftot
  integer,intent(in) :: sigmak_ibz
  integer,intent(in) :: use_aerhor
  type(crystal_t),intent(in) :: Cryst
  type(dataset_type),intent(in) :: Dtset
  type(epsilonm1_results),intent(inout) :: Er
  type(gsphere_t),intent(in) :: Gsph_Max
  type(gsphere_t),intent(in) :: Gsph_c
  type(kmesh_t),intent(in) :: Kmesh
  type(littlegroup_t),intent(in) :: Ltg_k
  type(ppmodel_t),intent(inout) :: PPm
  type(pseudopotential_type),intent(in) :: Psps
  type(ebands_t),target,intent(in) :: QP_BSt
  type(kmesh_t),intent(in) :: Qmesh
  type(sigparams_t),target,intent(in) :: Sigp
  type(sigma_t),intent(in) :: Sr
  type(vcoul_t),intent(in) :: Vcp
  type(wfd_t),target,intent(inout) :: Wfd
  type(wfd_t),target,intent(inout) :: Wfdf
  type(pawang_type),intent(in) :: pawang
  integer,intent(in) :: gwc_ngfft(18)
  integer,intent(in) :: rho_ngfft(18)
  type(paw_pwaves_lmn_t),intent(in) :: Paw_onsite(Cryst%natom*Psps%usepaw)
  type(pawpwff_t),intent(in) :: Paw_pwff(Psps%ntypat*Psps%usepaw)
  type(pawfgrtab_type),intent(inout) :: Pawfgrtab(Cryst%natom*Psps%usepaw)
  type(pawtab_type),intent(in) :: Pawtab(Psps%ntypat)
  real(dp),intent(in) :: aepaw_rhor(rho_nfftot,Wfd%nspden*use_aerhor)
  type(esymm_t),target,intent(in) :: allQP_sym(Wfd%nkibz,Wfd%nsppol)
  real(dp),intent(in) :: rhor(rho_nfftot,Wfd%nspden)
  complex(dpc),intent(out) :: sigcme_tmp(nomega_sigc,minbnd:maxbnd, &
  &         minbnd:maxbnd,Wfd%nsppol*Sigp%nsig_ab)
 end subroutine calc_sigc_me
end interface

interface
 subroutine calc_sigc_pole_cd(npwc,npwx,nspinor,ncoeff,nomega,nomegae,nomegaer,nomegaei,rhotwgp,&  
  &  omega,epsm1q,omegame0i,theta_mu_minus_e0i,ket,npoles_missing,calc_poles)
  use defs_basis
  implicit none
  integer,intent(in) :: ncoeff
  integer,intent(in) :: nomega
  integer,intent(in) :: nomegae
  integer,intent(in) :: nomegaei
  integer,intent(in) :: nomegaer
  integer,intent(inout) :: npoles_missing
  integer,intent(in) :: npwc
  integer,intent(in) :: npwx
  integer,intent(in) :: nspinor
  real(dp),intent(in) :: theta_mu_minus_e0i
  logical, intent(in), optional :: calc_poles(nomega)
  real(gwp) :: epsm1q(npwc,npwc,ncoeff)
  complex(gwpc),intent(inout) :: ket(nspinor*npwc,nomega)
  complex(dpc),intent(in) :: omega(nomegae)
  real(dp),intent(in) :: omegame0i(nomega)
  complex(gwpc),intent(in) :: rhotwgp(npwx*nspinor)
 end subroutine calc_sigc_pole_cd
end interface

interface
 subroutine calc_sigx_me(sigmak_ibz,ikcalc,minbnd,maxbnd,Cryst,QP_BSt,Sigp,Sr,Gsph_x,Vcp,Kmesh,Qmesh,&  
  &  Ltg_k,Pawtab,Pawang,Paw_pwff,Pawfgrtab,Paw_onsite,Psps,Wfd,Wfdf,allQP_sym,gwx_ngfft,ngfftf,&  
  &  prtvol,pawcross)
  use m_vcoul
  use m_pawtab
  use m_sigma
  use m_esymm
  use m_bz_mesh
  use m_pawang
  use m_wfd
  use m_paw_pwaves_lmn
  use m_pawfgrtab
  use m_crystal
  use m_gsphere
  use m_pawpwij
  use defs_datatypes
  use m_gwdefs
  implicit none
  integer,intent(in) :: ikcalc
  integer,intent(in) :: maxbnd
  integer,intent(in) :: minbnd
  integer,intent(in) :: pawcross
  integer,intent(in) :: prtvol
  integer,intent(in) :: sigmak_ibz
  type(crystal_t),intent(in) :: Cryst
  type(gsphere_t),intent(in) :: Gsph_x
  type(kmesh_t),intent(in) :: Kmesh
  type(littlegroup_t),intent(in) :: Ltg_k
  type(pawang_type),intent(in) :: Pawang
  type(pseudopotential_type),intent(in) :: Psps
  type(ebands_t),target,intent(in) :: QP_BSt
  type(kmesh_t),intent(in) :: Qmesh
  type(sigparams_t),target,intent(in) :: Sigp
  type(sigma_t),intent(inout) :: Sr
  type(vcoul_t),intent(in) :: Vcp
  type(wfd_t),target,intent(inout) :: Wfd
  type(wfd_t),target,intent(inout) :: Wfdf
  integer,intent(in) :: gwx_ngfft(18)
  integer,intent(in) :: ngfftf(18)
  type(paw_pwaves_lmn_t),intent(in) :: Paw_onsite(Cryst%natom*Psps%usepaw)
  type(pawpwff_t),intent(in) :: Paw_pwff(Psps%ntypat*Psps%usepaw)
  type(pawfgrtab_type),intent(inout) :: Pawfgrtab(Cryst%natom*Psps%usepaw)
  type(pawtab_type),intent(in) :: Pawtab(Psps%ntypat)
  type(esymm_t),target,intent(in) :: allQP_sym(Wfd%nkibz,Wfd%nsppol)
 end subroutine calc_sigx_me
end interface

interface
 subroutine calc_ucrpa(itypatcor,cryst,Kmesh,lpawu,M1_q_m,Qmesh,npwe,&  
  &  npw,nsym,rhot1_q_m,nomega,omegamin,omegamax,bandinf,bandsup,optimisation,ucvol,Wfd,fname)
  use m_wfd
  use m_bz_mesh
  use defs_basis
  use m_crystal
  implicit none
  integer, intent(in) :: bandinf
  integer, intent(in) :: bandsup
  integer, intent(in) :: itypatcor
  integer, intent(in) :: lpawu
  integer, intent(in) :: nomega
  integer, intent(in) :: npw
  integer, intent(in) :: npwe
  integer, intent(in) :: nsym
  type(crystal_t),intent(in) :: Cryst
  type(kmesh_t),intent(in) :: Kmesh
  type(kmesh_t),intent(in) :: Qmesh
  type(wfd_t),intent(inout) :: Wfd
  character(len=fnlen), intent(in) :: fname
  real(dp), intent(in) :: omegamax
  real(dp), intent(in) :: omegamin
  character(len=*), intent(in) :: optimisation
  real(dp), intent(in) :: ucvol
  complex(dpc), intent(in) :: M1_q_m(cryst%nattyp(itypatcor),Wfd%nspinor, &
  &         Wfd%nspinor,2*lpawu+1,2*lpawu+1,npw,Qmesh%nibz)
  complex(dpc), intent(in) :: rhot1_q_m(cryst%nattyp(itypatcor), &
  &         Wfd%nspinor,Wfd%nspinor,2*lpawu+1,2*lpawu+1,npw,Qmesh%nibz)
 end subroutine calc_ucrpa
end interface

interface
 subroutine calc_vhxc_me(Wfd,Mflags,Mels,Cryst,Dtset,gsqcutf_eff,nfftf,ngfftf,&  
  &  vtrial,vhartr,vxc,Psps,Pawtab,Paw_an,Pawang,Pawfgrtab,Paw_ij,dijexc_core,&  
  &  rhor,rhog,usexcnhat,nhat,nhatgr,nhatgrdim,kstab,&  
  &  taug,taur) ! optional arguments
  use m_melemts
  use m_pawtab
  use m_paw_ij
  use m_pawang
  use m_wfd
  use defs_abitypes
  use m_pawfgrtab
  use m_paw_an
  use m_crystal
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: nfftf
  integer,intent(in) :: nhatgrdim
  integer,intent(in) :: usexcnhat
  type(crystal_t),intent(in) :: Cryst
  type(dataset_type),intent(in) :: Dtset
  type(melements_t),intent(out) :: Mels
  type(melflags_t),intent(in) :: Mflags
  type(pawang_type),intent(in) :: Pawang
  type(pseudopotential_type),intent(in) :: Psps
  type(wfd_t),target,intent(inout) :: Wfd
  real(dp),intent(in) :: gsqcutf_eff
  integer,intent(in) :: ngfftf(18)
  type(paw_an_type),intent(in) :: Paw_an(Cryst%natom)
  type(paw_ij_type),intent(inout) :: Paw_ij(Cryst%natom)
  type(pawfgrtab_type),intent(inout) :: Pawfgrtab(Cryst%natom)
  type(pawtab_type),intent(in) :: Pawtab(Cryst%ntypat*Wfd%usepaw)
  real(dp),intent(in) :: dijexc_core(:,:,:)
  integer,intent(in) :: kstab(2,Wfd%nkibz,Wfd%nsppol)
  real(dp),intent(in) :: nhat(nfftf,Wfd%nspden*Wfd%usepaw)
  real(dp),intent(in) :: nhatgr(nfftf,Wfd%nspden,3*nhatgrdim)
  real(dp),intent(in) :: rhog(2,nfftf)
  real(dp),intent(in) :: rhor(nfftf,Wfd%nspden)
  real(dp),intent(in),optional :: taug(2,nfftf*Dtset%usekden)
  real(dp),intent(in),optional :: taur(nfftf,Wfd%nspden*Dtset%usekden)
  real(dp),intent(in) :: vhartr(nfftf)
  real(dp),intent(in) :: vtrial(nfftf,Wfd%nspden)
  real(dp),intent(in) :: vxc(nfftf,Wfd%nspden)
 end subroutine calc_vhxc_me
end interface

interface
 subroutine update_cprj(natom,nkibz,nbnds,nsppol,nspinor,m_lda_to_qp,dimlmn,Cprj_ibz)
  use defs_basis
  use m_pawcprj
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nbnds
  integer,intent(in) :: nkibz
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  type(pawcprj_type),intent(inout) :: Cprj_ibz(natom,nspinor*nbnds*nkibz*nsppol)
  integer,intent(in) :: dimlmn(natom)
  complex(dpc),intent(in) :: m_lda_to_qp(nbnds,nbnds,nkibz,nsppol)
 end subroutine update_cprj
end interface

interface
 subroutine cchi0(use_tr,Dtset,Cryst,qpoint,Ep,Psps,Kmesh,QP_BSt,Gsph_epsG0,&  
  &  Pawtab,Pawang,Paw_pwff,Pawfgrtab,Paw_onsite,nbvw,ngfft_gw,nfftot_gw,ngfftf,nfftf_tot,&  
  &  chi0,ktabr,ktabrf,Ltg_q,chi0_sumrule,Wfd,Wfdf)
  use m_pawtab
  use m_bz_mesh
  use m_pawang
  use m_wfd
  use m_paw_pwaves_lmn
  use defs_abitypes
  use m_pawfgrtab
  use m_gsphere
  use m_crystal
  use defs_basis
  use m_pawpwij
  use defs_datatypes
  use m_gwdefs
  implicit none
  integer,intent(in) :: nbvw
  integer,intent(in) :: nfftf_tot
  integer,intent(in) :: nfftot_gw
  type(crystal_t),intent(in) :: Cryst
  type(dataset_type),intent(in) :: Dtset
  type(em1params_t),intent(in) :: Ep
  type(gsphere_t),intent(in) :: Gsph_epsG0
  type(kmesh_t),intent(in) :: Kmesh
  type(littlegroup_t),intent(in) :: Ltg_q
  type(pawang_type),intent(in) :: Pawang
  type(pseudopotential_type),intent(in) :: Psps
  type(ebands_t),target,intent(in) :: QP_BSt
  type(wfd_t),target,intent(inout) :: Wfd
  type(wfd_t),target,intent(inout) :: Wfdf
  logical,intent(in) :: use_tr
  integer,intent(in) :: ngfft_gw(18)
  integer,intent(in) :: ngfftf(18)
  type(paw_pwaves_lmn_t),intent(in) :: Paw_onsite(Cryst%natom)
  type(pawpwff_t),intent(in) :: Paw_pwff(Psps%ntypat*Psps%usepaw)
  type(pawfgrtab_type),intent(inout) :: Pawfgrtab(Cryst%natom*Psps%usepaw)
  type(pawtab_type),intent(in) :: Pawtab(Psps%ntypat*Psps%usepaw)
  complex(gwpc),intent(out) :: chi0(Ep%npwe*Ep%nI,Ep%npwe*Ep%nJ,Ep%nomega)
  real(dp),intent(out) :: chi0_sumrule(Ep%npwe)
  integer,intent(in) :: ktabr(nfftot_gw,Kmesh%nbz)
  integer,intent(in) :: ktabrf(nfftf_tot*Dtset%pawcross,Kmesh%nbz)
  real(dp),intent(in) :: qpoint(3)
 end subroutine cchi0
end interface

interface
 subroutine cchi0q0(use_tr,Dtset,Cryst,Ep,Psps,Kmesh,QP_BSt,KS_BSt,Gsph_epsG0,&  
  &  Pawang,Pawrad,Pawtab,Paw_ij,Paw_pwff,Pawfgrtab,Paw_onsite,ktabr,ktabrf,nbvw,ngfft_gw,&  
  &  nfftot_gw,ngfftf,nfftf_tot,chi0,chi0_head,chi0_lwing,chi0_uwing,Ltg_q,chi0_sumrule,Wfd,Wfdf)
  use m_pawtab
  use m_pawrad
  use defs_basis
  use m_bz_mesh
  use m_paw_ij
  use m_pawang
  use m_wfd
  use m_paw_pwaves_lmn
  use defs_abitypes
  use m_pawfgrtab
  use m_crystal
  use m_gsphere
  use m_pawpwij
  use defs_datatypes
  use m_gwdefs
  implicit none
  integer,intent(in) :: nbvw
  integer,intent(in) :: nfftf_tot
  integer,intent(in) :: nfftot_gw
  type(crystal_t),intent(in) :: Cryst
  type(dataset_type),intent(in) :: Dtset
  type(em1params_t),intent(in) :: Ep
  type(gsphere_t),intent(in) :: Gsph_epsG0
  type(ebands_t),target,intent(in) :: KS_BSt
  type(kmesh_t),intent(in) :: Kmesh
  type(littlegroup_t),intent(in) :: Ltg_q
  type(pawang_type),intent(in) :: Pawang
  type(pseudopotential_type),intent(in) :: Psps
  type(ebands_t),target,intent(in) :: QP_BSt
  type(wfd_t),target,intent(inout) :: Wfd
  type(wfd_t),target,intent(inout) :: Wfdf
  logical,intent(in) :: use_tr
  integer,intent(in) :: ngfft_gw(18)
  integer,intent(in) :: ngfftf(18)
  type(paw_ij_type),intent(in) :: Paw_ij(Cryst%natom*Psps%usepaw)
  type(paw_pwaves_lmn_t),intent(in) :: Paw_onsite(Cryst%natom)
  type(pawpwff_t),intent(in) :: Paw_pwff(Psps%ntypat*Psps%usepaw)
  type(pawfgrtab_type),intent(inout) :: Pawfgrtab(Cryst%natom*Psps%usepaw)
  type(pawrad_type),intent(in) :: Pawrad(Psps%ntypat*Psps%usepaw)
  type(pawtab_type),intent(in) :: Pawtab(Psps%ntypat*Psps%usepaw)
  complex(gwpc),intent(out) :: chi0(Ep%npwe,Ep%npwe,Ep%nomega)
  complex(dpc),intent(out) :: chi0_head(3,3,Ep%nomega)
  complex(dpc),intent(out) :: chi0_lwing(Ep%npwe*Ep%nI,Ep%nomega,3)
  real(dp),intent(out) :: chi0_sumrule(Ep%npwe)
  complex(dpc),intent(out) :: chi0_uwing(Ep%npwe*Ep%nJ,Ep%nomega,3)
  integer,intent(in) :: ktabr(nfftot_gw,Kmesh%nbz)
  integer,intent(in) :: ktabrf(nfftf_tot*Dtset%pawcross,Kmesh%nbz)
 end subroutine cchi0q0
end interface

interface
 subroutine chi0q0_intraband(Wfd,Cryst,Ep,Psps,BSt,Gsph_epsG0,Pawang,Pawrad,Pawtab,Paw_ij,Paw_pwff,use_tr,usepawu,&  
  &  ngfft_gw,chi0,chi0_head,chi0_lwing,chi0_uwing)
  use m_pawtab
  use m_pawrad
  use defs_basis
  use m_paw_ij
  use m_pawang
  use m_wfd
  use m_crystal
  use m_gsphere
  use m_pawpwij
  use defs_datatypes
  use m_gwdefs
  implicit none
  integer,intent(in) :: usepawu
  type(ebands_t),intent(in) :: BSt
  type(crystal_t),intent(in) :: Cryst
  type(em1params_t),intent(in) :: Ep
  type(gsphere_t),intent(in) :: Gsph_epsG0
  type(pawang_type),intent(in) :: Pawang
  type(pseudopotential_type),intent(in) :: Psps
  type(wfd_t),target,intent(inout) :: Wfd
  logical,intent(in) :: use_tr
  integer,intent(in) :: ngfft_gw(18)
  type(paw_ij_type),intent(in) :: Paw_ij(Cryst%natom*Psps%usepaw)
  type(pawpwff_t),intent(in) :: Paw_pwff(Psps%ntypat*Psps%usepaw)
  type(pawrad_type),intent(in) :: Pawrad(Psps%ntypat*Psps%usepaw)
  type(pawtab_type),intent(in) :: Pawtab(Psps%ntypat*Psps%usepaw)
  complex(gwpc),intent(out) :: chi0(Ep%npwe*Ep%nI,Ep%npwe*Ep%nJ,Ep%nomega)
  complex(dpc),intent(out) :: chi0_head(3,3,Ep%nomega)
  complex(dpc),intent(out) :: chi0_lwing(Ep%npwe*Ep%nI,Ep%nomega,3)
  complex(dpc),intent(out) :: chi0_uwing(Ep%npwe*Ep%nJ,Ep%nomega,3)
 end subroutine chi0q0_intraband
end interface

interface
 subroutine cohsex_me(sigmak_ibz,ikcalc,nomega_sigc,minbnd,maxbnd,Cryst,QP_BSt,Sigp,Sr,Er,Gsph_c,Vcp,&  
  &  Kmesh,Qmesh,Ltg_k,Pawtab,Pawang,Paw_pwff,Psps,Wfd,allQP_sym,gwc_ngfft,iomode,prtvol,sigcme_tmp)
  use m_vcoul
  use m_pawtab
  use m_sigma
  use m_esymm
  use defs_basis
  use m_bz_mesh
  use m_pawang
  use m_wfd
  use m_crystal
  use m_gsphere
  use m_pawpwij
  use defs_datatypes
  use m_gwdefs
  use m_screening
  implicit none
  integer,intent(in) :: ikcalc
  integer,intent(in) :: iomode
  integer,intent(in) :: maxbnd
  integer,intent(in) :: minbnd
  integer,intent(in) :: nomega_sigc
  integer,intent(in) :: prtvol
  integer,intent(in) :: sigmak_ibz
  type(crystal_t),intent(in) :: Cryst
  type(epsilonm1_results),intent(inout) :: Er
  type(gsphere_t),intent(in) :: Gsph_c
  type(kmesh_t),intent(in) :: Kmesh
  type(littlegroup_t),intent(in) :: Ltg_k
  type(pseudopotential_type),intent(in) :: Psps
  type(ebands_t),target,intent(in) :: QP_BSt
  type(kmesh_t),intent(in) :: Qmesh
  type(sigparams_t),target,intent(in) :: Sigp
  type(sigma_t),intent(in) :: Sr
  type(vcoul_t),intent(in) :: Vcp
  type(wfd_t),target,intent(inout) :: Wfd
  type(pawang_type),intent(in) :: pawang
  integer,intent(in) :: gwc_ngfft(18)
  type(pawpwff_t),intent(in) :: Paw_pwff(Psps%ntypat*Psps%usepaw)
  type(pawtab_type),intent(in) :: Pawtab(Psps%ntypat)
  type(esymm_t),target,intent(in) :: allQP_sym(Wfd%nkibz,Wfd%nsppol)
  complex(dpc),intent(out) :: sigcme_tmp(nomega_sigc,minbnd:maxbnd, &
  &         minbnd:maxbnd,Wfd%nsppol*Sigp%nsig_ab)
 end subroutine cohsex_me
end interface

interface
 subroutine check_zarot(npwvec,Cryst,ngfft,gvec,psps,pawang,grottb,grottbm1)
  use m_pawang
  use defs_datatypes
  use m_crystal
  implicit none
  integer,intent(in) :: npwvec
  type(crystal_t),intent(in) :: Cryst
  type(pawang_type),intent(in) :: pawang
  type(pseudopotential_type),intent(in) :: psps
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: grottb(npwvec,Cryst%timrev,Cryst%nsym)
  integer,intent(in) :: grottbm1(npwvec,Cryst%timrev,Cryst%nsym)
  integer,intent(in) :: gvec(3,npwvec)
 end subroutine check_zarot
end interface

interface
 subroutine paw_check_symcprj(Wfd,ik_bz,band,spin,sym_mode,Cryst,Kmesh,Psps,Pawtab,Pawang,Cprj_bz)
  use m_bz_mesh
  use m_pawang
  use m_wfd
  use m_pawcprj
  use m_crystal
  use m_pawtab
  use defs_datatypes
  implicit none
  integer,intent(in) :: band
  integer,intent(in) :: ik_bz
  integer,intent(in) :: spin
  integer,intent(in) :: sym_mode
  type(crystal_t),intent(in) :: Cryst
  type(kmesh_t),intent(in) :: Kmesh
  type(pawang_type),intent(in) :: Pawang
  type(pseudopotential_type),intent(in) :: Psps
  type(wfd_t),intent(inout) :: Wfd
  type(pawcprj_type),intent(out) :: Cprj_bz(Cryst%natom,Wfd%nspinor)
  type(pawtab_type),intent(in) :: Pawtab(Psps%ntypat*Psps%usepaw)
 end subroutine paw_check_symcprj
end interface

interface
 function dotproductqrc(r,c,b1,b2,b3)
  use defs_basis
  implicit none
  complex(gwpc) :: dotproductqrc
  real(dp),intent(in) :: b1(3)
  real(dp),intent(in) :: b2(3)
  real(dp),intent(in) :: b3(3)
  complex(gwpc),intent(in) :: c(3)
  real(dp),intent(in) :: r(3)
 end function dotproductqrc
end interface

interface
 subroutine fsumrule(nomega,omega,eps,omegaplasma,method)
  use defs_basis
  implicit none
  integer,intent(in) :: method
  integer,intent(in) :: nomega
  real(dp),intent(in) :: omegaplasma
  real(dp),intent(in) :: eps(nomega)
  real(dp),intent(in) :: omega(nomega)
 end subroutine fsumrule
end interface

interface
 subroutine make_transitions(Wfd,chi0alg,nbnds,nbvw,nsppol,symchi,timrev,TOL_DELTA_OCC,&  
  &  max_rest,min_rest,my_max_rest,my_min_rest,Kmesh,Ltg_q,gw_energy,occ,qpoint,bbp_ks_distrb)
  use m_wfd
  use defs_basis
  use m_bz_mesh
  implicit none
  integer,intent(in) :: chi0alg
  integer,intent(in) :: nbnds
  integer,intent(in) :: nbvw
  integer,intent(in) :: nsppol
  integer,intent(in) :: symchi
  integer,intent(in) :: timrev
  type(kmesh_t),intent(in) :: Kmesh
  type(littlegroup_t),intent(in) :: Ltg_q
  real(dp),intent(in) :: TOL_DELTA_OCC
  type(wfd_t),intent(in) :: Wfd
  real(dp),intent(out) :: max_rest
  real(dp),intent(out) :: min_rest
  real(dp),intent(out) :: my_max_rest
  real(dp),intent(out) :: my_min_rest
  integer,intent(in) :: bbp_ks_distrb(Wfd%mband,Wfd%mband,Kmesh%nbz,Wfd%nsppol)
  real(dp),intent(in) :: gw_energy(nbnds,Kmesh%nibz,nsppol)
  real(dp),intent(in) :: occ(nbnds,Kmesh%nibz,nsppol)
  real(dp),intent(in) :: qpoint(3)
 end subroutine make_transitions
end interface

interface
 subroutine sigma_distribute_bks(Wfd,Kmesh,Ltg_kgw,Qmesh,nsppol,can_symmetrize,kptgw,mg0,my_nbks,proc_distrb,got,bks_mask,global)
  use m_wfd
  use m_bz_mesh
  use defs_basis
  implicit none
  integer,intent(out) :: my_nbks
  integer,intent(in) :: nsppol
  type(kmesh_t),intent(in) :: Kmesh
  type(littlegroup_t),intent(in) :: Ltg_kgw
  type(kmesh_t),intent(in) :: Qmesh
  type(wfd_t),intent(inout) :: Wfd
  logical,optional,intent(in) :: global
  integer,intent(in) :: mg0(3)
  logical,optional,intent(in) :: bks_mask(Wfd%mband,Kmesh%nbz,nsppol)
  logical,intent(in) :: can_symmetrize(Wfd%nsppol)
  integer,optional,intent(inout) :: got(Wfd%nproc)
  real(dp),intent(in) :: kptgw(3)
  integer,intent(out) :: proc_distrb(Wfd%mband,Kmesh%nbz,nsppol)
 end subroutine sigma_distribute_bks
end interface

interface
 subroutine chi0_bbp_mask(Ep,use_tr,QP_BSt,mband,ikmq_ibz,ik_ibz,spin,spin_fact,bbp_mask)
  use defs_basis
  use defs_datatypes
  use m_gwdefs
  implicit none
  integer,intent(in) :: ik_ibz
  integer,intent(in) :: ikmq_ibz
  integer,intent(in) :: mband
  integer,intent(in) :: spin
  type(em1params_t),intent(in) :: Ep
  type(ebands_t),target,intent(in) :: QP_BSt
  real(dp),intent(in) :: spin_fact
  logical,intent(in) :: use_tr
  logical,intent(out) :: bbp_mask(mband,mband)
 end subroutine chi0_bbp_mask
end interface

interface
 subroutine completechi0_deltapart(ik_bz,qzero,symchi,npwe,npwvec,nomega,nspinor,&  
  &  nfftot,ngfft,igfft0,Gsph_FFT,Ltg_q,green_enhigh_w,wfwfg,chi0)
  use m_gsphere
  use defs_basis
  use m_bz_mesh
  implicit none
  integer,intent(in) :: ik_bz
  integer,intent(in) :: nfftot
  integer,intent(in) :: nomega
  integer,intent(in) :: npwe
  integer,intent(in) :: npwvec
  integer,intent(in) :: nspinor
  integer,intent(in) :: symchi
  type(gsphere_t),intent(in) :: Gsph_FFT
  type(littlegroup_t),intent(in) :: Ltg_q
  logical,intent(in) :: qzero
  integer,intent(in) :: ngfft(18)
  complex(gwpc),intent(inout) :: chi0(npwe,npwe,nomega)
  complex(dpc),intent(in) :: green_enhigh_w(nomega)
  integer,intent(in) :: igfft0(npwvec)
  complex(gwpc),intent(in) :: wfwfg(nfftot*nspinor**2)
 end subroutine completechi0_deltapart
end interface

interface
 subroutine output_chi0sumrule(qeq0,iq,npwe,omegaplasma,chi0sumrule,epsm1_w0,vc_sqrt)
  use defs_basis
  implicit none
  integer,intent(in) :: iq
  integer,intent(in) :: npwe
  real(dp),intent(in) :: omegaplasma
  logical,intent(in) :: qeq0
  real(dp),intent(inout) :: chi0sumrule(npwe)
  complex(gwpc),intent(in) :: epsm1_w0(npwe,npwe)
  complex(gwpc),intent(in) :: vc_sqrt(npwe)
 end subroutine output_chi0sumrule
end interface

interface
 subroutine accumulate_chi0sumrule(ik_bz,symchi,npwe,factor,delta_ene,&  
  &  Ltg_q,Gsph_epsG0,npwepG0,rhotwg,chi0sumrule)
  use defs_basis
  use m_gsphere
  use m_bz_mesh
  implicit none
  integer,intent(in) :: ik_bz
  integer,intent(in) :: npwe
  integer,intent(in) :: npwepG0
  integer,intent(in) :: symchi
  type(gsphere_t),intent(in) :: Gsph_epsG0
  type(littlegroup_t),target,intent(in) :: Ltg_q
  real(dp),intent(in) :: delta_ene
  real(dp),intent(in) :: factor
  real(dp),intent(inout) :: chi0sumrule(npwe)
  complex(gwpc),intent(in) :: rhotwg(npwepG0)
 end subroutine accumulate_chi0sumrule
end interface

interface
 subroutine mlwfovlp_qp(cg,Cprj_BZ,dtset,dtfil,eigen,mband,mcg,mcprj,mkmem,mpw,natom,&  
  &  nkpt,npwarr,nspden,nsppol,ntypat,Hdr,Pawtab,rprimd,MPI_enreg)
  use defs_basis
  use defs_abitypes
  use m_pawcprj
  use m_pawtab
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: mcg
  integer,intent(in) :: mcprj
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspden
  integer,intent(in) :: nsppol
  integer,intent(in) :: ntypat
  type(hdr_type),intent(in) :: Hdr
  type(mpi_type),intent(inout) :: MPI_enreg
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  type(pawcprj_type),target,intent(inout) :: Cprj_BZ(natom,mcprj)
  type(pawtab_type),intent(in) :: Pawtab(ntypat*Dtset%usepaw)
  real(dp),intent(inout) :: cg(2,mcg)
  real(dp),intent(inout) :: eigen(mband*nkpt*nsppol)
  integer,intent(in) :: npwarr(nkpt)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine mlwfovlp_qp
end interface

interface
 subroutine calc_coh(nspinor,nsig_ab,nfftot,ngfft,npwc,gvec,wfg2_jk,epsm1q_o,vc_sqrt,i_sz,iqibz,same_band,sigcohme)
  use defs_basis
  implicit none
  integer,intent(in) :: iqibz
  integer,intent(in) :: nfftot
  integer,intent(in) :: npwc
  integer,intent(in) :: nsig_ab
  integer,intent(in) :: nspinor
  real(dp),intent(in) :: i_sz
  logical,intent(in) :: same_band
  integer,intent(in) :: ngfft(18)
  complex(gwpc),intent(in) :: epsm1q_o(npwc,npwc)
  integer,intent(in) :: gvec(3,npwc)
  complex(gwpc),intent(out) :: sigcohme(nsig_ab)
  complex(gwpc),intent(in) :: vc_sqrt(npwc)
  complex(gwpc),intent(in) :: wfg2_jk(nsig_ab*nfftot)
 end subroutine calc_coh
end interface

interface
 subroutine calc_coh_comp(iqibz,i_sz,same_band,nspinor,nsig_ab,ediff,npwc,gvec,&  
  &  ngfft,nfftot,wfg2_jk,vc_sqrt,botsq,otq,sigcohme)
  use defs_basis
  implicit none
  integer,intent(in) :: iqibz
  integer,intent(in) :: nfftot
  integer,intent(in) :: npwc
  integer,intent(in) :: nsig_ab
  integer,intent(in) :: nspinor
  real(dp),intent(in) :: ediff
  real(dp),intent(in) :: i_sz
  logical,intent(in) :: same_band
  integer,intent(in) :: ngfft(18)
  complex(gwpc),intent(in) :: botsq(npwc,npwc)
  integer,intent(in) :: gvec(3,npwc)
  complex(gwpc),intent(in) :: otq(npwc,npwc)
  complex(gwpc),intent(out) :: sigcohme(nsig_ab)
  complex(gwpc),intent(in) :: vc_sqrt(npwc)
  complex(gwpc),intent(in) :: wfg2_jk(nsig_ab*nfftot)
 end subroutine calc_coh_comp
end interface

interface
 subroutine paw_qpscgw(Wfd,nscf,nfftf,ngfftf,Dtset,Cryst,Kmesh,Psps,QP_BSt,&  
  &  Pawang,Pawrad,Pawtab,Pawfgrtab,prev_Pawrhoij,&  
  &  QP_pawrhoij,QP_paw_ij,QP_paw_an,QP_energies,qp_nhat,nhatgrdim,qp_nhatgr,qp_compch_sph,qp_compch_fft)
  use m_pawrad
  use defs_basis
  use m_bz_mesh
  use m_paw_ij
  use m_pawang
  use m_pawrhoij
  use m_wfd
  use m_energies
  use defs_abitypes
  use m_pawfgrtab
  use m_paw_an
  use m_crystal
  use m_pawtab
  use defs_datatypes
  implicit none
  integer,intent(in) :: nfftf
  integer,intent(in) :: nhatgrdim
  integer,intent(in) :: nscf
  type(crystal_t),intent(in) :: Cryst
  type(dataset_type),intent(in) :: Dtset
  type(kmesh_t),intent(in) :: Kmesh
  type(pawang_type),intent(in) :: Pawang
  type(pseudopotential_type),intent(in) :: Psps
  type(ebands_t),intent(in) :: QP_BSt
  type(energies_type),intent(inout) :: QP_energies
  type(wfd_t),intent(inout) :: Wfd
  real(dp),intent(out) :: qp_compch_fft
  real(dp),intent(out) :: qp_compch_sph
  integer,intent(in) :: ngfftf(18)
  type(pawfgrtab_type),intent(inout) :: Pawfgrtab(Cryst%natom)
  type(pawrad_type),intent(in) :: Pawrad(Psps%ntypat*Psps%usepaw)
  type(pawtab_type),intent(in) :: Pawtab(Psps%ntypat*Psps%usepaw)
  type(paw_an_type),intent(inout) :: QP_paw_an(Cryst%natom)
  type(paw_ij_type),intent(out) :: QP_paw_ij(Cryst%natom)
  type(pawrhoij_type),intent(out) :: QP_pawrhoij(Cryst%natom)
  type(pawrhoij_type),intent(inout) :: prev_Pawrhoij(Cryst%natom)
  real(dp),intent(out) :: qp_nhat(nfftf,Dtset%nspden)
  real(dp),intent(out) :: qp_nhatgr(nfftf,Dtset%nspden,3*nhatgrdim)
 end subroutine paw_qpscgw
end interface

interface
 subroutine prep_calc_ucrpa(sigmak_ibz,ikcalc,itypatcor,minbnd,maxbnd,Cryst,QP_BSt,Sigp,Gsph_x,Vcp,Kmesh,Qmesh,lpawu,&  
  &  M1_q_m,Pawtab,Pawang,Paw_pwff,Pawfgrtab,Paw_onsite,&  
  &  Psps,Wfd,Wfdf,allQP_sym,gwx_ngfft,ngfftf,&  
  &  prtvol,pawcross,rhot1_q_m)
  use m_vcoul
  use m_pawtab
  use defs_basis
  use m_esymm
  use m_bz_mesh
  use m_pawang
  use m_wfd
  use m_paw_pwaves_lmn
  use m_pawfgrtab
  use m_crystal
  use m_gsphere
  use m_pawpwij
  use defs_datatypes
  use m_gwdefs
  implicit none
  integer,intent(in) :: ikcalc
  integer,intent(in) :: itypatcor
  integer,intent(in) :: lpawu
  integer,intent(in) :: maxbnd
  integer,intent(in) :: minbnd
  integer,intent(in) :: pawcross
  integer,intent(in) :: prtvol
  integer,intent(in) :: sigmak_ibz
  type(crystal_t),intent(in) :: Cryst
  type(gsphere_t),intent(in) :: Gsph_x
  type(kmesh_t),intent(in) :: Kmesh
  type(pawang_type),intent(in) :: Pawang
  type(pseudopotential_type),intent(in) :: Psps
  type(ebands_t),target,intent(in) :: QP_BSt
  type(kmesh_t),intent(in) :: Qmesh
  type(sigparams_t),target,intent(in) :: Sigp
  type(vcoul_t),intent(in) :: Vcp
  type(wfd_t),target,intent(inout) :: Wfd
  type(wfd_t),target,intent(inout) :: Wfdf
  integer,intent(in) :: gwx_ngfft(18)
  integer,intent(in) :: ngfftf(18)
  complex(dpc), intent(out) :: M1_q_m(cryst%nattyp(itypatcor),Wfd%nspinor, &
  &         Wfd%nspinor,2*lpawu+1,2*lpawu+1,sigp%npwx,Qmesh%nibz)
  type(paw_pwaves_lmn_t),intent(in) :: Paw_onsite(Cryst%natom)
  type(pawpwff_t),intent(in) :: Paw_pwff(Psps%ntypat*Psps%usepaw)
  type(pawfgrtab_type),intent(inout) :: Pawfgrtab(Cryst%natom*Psps%usepaw)
  type(pawtab_type),intent(in) :: Pawtab(Psps%ntypat)
  type(esymm_t),target,intent(in) :: allQP_sym(Wfd%nkibz,Wfd%nsppol)
  complex(dpc), intent(out) :: rhot1_q_m(cryst%nattyp(itypatcor), &
  &         Wfd%nspinor,Wfd%nspinor,2*lpawu+1,2*lpawu+1,sigp%npwx,Qmesh%nibz)
 end subroutine prep_calc_ucrpa
end interface

interface
 subroutine random_stopping_power(iqibz,npvel,pvelmax,Ep,Gsph_epsG0,Qmesh,Vcp,Cryst,Dtfil,epsm1,rspower)
  use m_vcoul
  use m_bz_mesh
  use defs_abitypes
  use m_gsphere
  use m_crystal
  use defs_basis
  use m_gwdefs
  implicit none
  integer,intent(in) :: iqibz
  integer,intent(in) :: npvel
  type(crystal_t),intent(in) :: Cryst
  type(datafiles_type),intent(in) :: Dtfil
  type(em1params_t),intent(in) :: Ep
  type(gsphere_t),intent(in) :: Gsph_epsG0
  type(kmesh_t),intent(in) :: Qmesh
  type(vcoul_t),intent(in) :: Vcp
  complex(gwpc),intent(in) :: epsm1(Ep%npwe,Ep%npwe,Ep%nomega)
  real(dp),intent(in) :: pvelmax(3)
  real(dp),intent(inout) :: rspower(npvel)
 end subroutine random_stopping_power
end interface

interface
 subroutine read_plowannier(cryst,bandinf,bandsup,coeffW_BZ,itypatcor,Kmesh,lcor,luwindow,nspinor,nsppol,pawang,prtvol,ucrpa_bands)
  use defs_basis
  use m_pawang
  use m_bz_mesh
  use m_crystal
  implicit none
  integer, intent(out) :: bandinf
  integer, intent(out) :: bandsup
  integer, intent(out) :: itypatcor
  integer, intent(out) :: lcor
  integer, intent(in) :: nspinor
  integer, intent(in) :: nsppol
  integer, intent(in) :: prtvol
  type(crystal_t),intent(in) :: Cryst
  type(kmesh_t),intent(in) :: Kmesh
  type(pawang_type),intent(in) :: Pawang
  logical, intent(inout) :: luwindow
  integer, intent(in) :: ucrpa_bands(2)
  complex(dpc), allocatable, intent(inout) :: coeffW_BZ(:,:,:,:,:,:)
 end subroutine read_plowannier
end interface

interface
 subroutine setup_screening(codvsn,acell,rprim,ngfftf,wfk_fname,dtfil,Dtset,Psps,Pawtab,&  
  &  ngfft_gw,Hdr_wfk,Hdr_out,Cryst,Kmesh,Qmesh,KS_BSt,Ltg_q,Gsph_epsG0,Gsph_wfn,Vcp,Ep,comm)
  use m_vcoul
  use m_pawtab
  use m_bz_mesh
  use defs_abitypes
  use m_gsphere
  use m_crystal
  use defs_basis
  use defs_datatypes
  use m_gwdefs
  implicit none
  integer,intent(in) :: comm
  type(crystal_t),intent(out) :: Cryst
  type(dataset_type),intent(inout) :: Dtset
  type(em1params_t),intent(out) :: Ep
  type(gsphere_t),intent(out) :: Gsph_epsG0
  type(gsphere_t),intent(out) :: Gsph_wfn
  type(hdr_type),intent(out) :: Hdr_out
  type(hdr_type),intent(out) :: Hdr_wfk
  type(ebands_t),intent(out) :: KS_BSt
  type(kmesh_t),intent(out) :: Kmesh
  type(pseudopotential_type),intent(in) :: Psps
  type(kmesh_t),intent(out) :: Qmesh
  type(vcoul_t),intent(out) :: Vcp
  character(len=6),intent(in) :: codvsn
  type(datafiles_type),intent(in) :: dtfil
  character(len=fnlen),intent(in) :: wfk_fname
  integer,intent(out) :: ngfft_gw(18)
  integer,intent(in) :: ngfftf(18)
  type(littlegroup_t),pointer :: Ltg_q(:)
  type(pawtab_type),intent(in) :: Pawtab(Psps%ntypat*Dtset%usepaw)
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: rprim(3,3)
 end subroutine setup_screening
end interface

interface
 subroutine chi0_bksmask(Dtset,Ep,Kmesh,nbvw,nbcw,my_rank,nprocs,bks_mask,keep_ur,ierr)
  use m_bz_mesh
  use defs_abitypes
  use m_gwdefs
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: my_rank
  integer,intent(in) :: nbcw
  integer,intent(in) :: nbvw
  integer,intent(in) :: nprocs
  type(dataset_type),intent(in) :: Dtset
  type(em1params_t),intent(in) :: Ep
  type(kmesh_t),intent(in) :: Kmesh
  logical,intent(out) :: bks_mask(Ep%nbnds,Kmesh%nibz,Dtset%nsppol)
  logical,intent(out) :: keep_ur(Ep%nbnds,Kmesh%nibz,Dtset%nsppol)
 end subroutine chi0_bksmask
end interface

interface
 subroutine setup_sigma(codvsn,wfk_fname,acell,rprim,ngfftf,Dtset,Dtfil,Psps,Pawtab,&  
  &  gwx_ngfft,gwc_ngfft,Hdr_wfk,Hdr_out,Cryst,Kmesh,Qmesh,KS_BSt,Gsph_Max,Gsph_x,Gsph_c,Vcp,Er,Sigp,comm)
  use m_vcoul
  use m_pawtab
  use m_bz_mesh
  use defs_abitypes
  use m_gsphere
  use m_crystal
  use defs_basis
  use defs_datatypes
  use m_gwdefs
  use m_screening
  implicit none
  integer,intent(in) :: comm
  type(crystal_t),intent(out) :: Cryst
  type(datafiles_type),intent(in) :: Dtfil
  type(dataset_type),intent(inout) :: Dtset
  type(epsilonm1_results),intent(out) :: Er
  type(gsphere_t),intent(out) :: Gsph_Max
  type(gsphere_t),intent(out) :: Gsph_c
  type(gsphere_t),intent(out) :: Gsph_x
  type(hdr_type),intent(out) :: Hdr_out
  type(hdr_type),intent(out) :: Hdr_wfk
  type(ebands_t),intent(out) :: KS_BSt
  type(kmesh_t),intent(out) :: Kmesh
  type(pseudopotential_type),intent(in) :: Psps
  type(kmesh_t),intent(out) :: Qmesh
  type(sigparams_t),intent(out) :: Sigp
  type(vcoul_t),intent(out) :: Vcp
  character(len=6),intent(in) :: codvsn
  character(len=*),intent(in) :: wfk_fname
  integer,intent(out) :: gwc_ngfft(18)
  integer,intent(out) :: gwx_ngfft(18)
  integer,intent(in) :: ngfftf(18)
  type(pawtab_type),intent(in) :: Pawtab(Psps%ntypat*Dtset%usepaw)
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: rprim(3,3)
 end subroutine setup_sigma
end interface

interface
 subroutine sigma_tables(Sigp,Kmesh,Bnd_sym)
  use m_bz_mesh
  use m_esymm
  use m_gwdefs
  implicit none
  type(kmesh_t),intent(in) :: Kmesh
  type(sigparams_t),intent(inout) :: Sigp
  type(esymm_t),optional,intent(in) :: Bnd_sym(Kmesh%nibz,Sigp%nsppol)
 end subroutine sigma_tables
end interface

interface
 subroutine sigma_bksmask(Dtset,Sigp,Kmesh,my_rank,nprocs,my_spins,bks_mask,keep_ur,ierr)
  use m_bz_mesh
  use defs_abitypes
  use m_gwdefs
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: my_rank
  integer,intent(in) :: nprocs
  type(dataset_type),intent(in) :: Dtset
  type(kmesh_t),intent(in) :: Kmesh
  type(sigparams_t),intent(in) :: Sigp
  integer,allocatable,intent(out) :: my_spins(:)
  logical,intent(out) :: bks_mask(Sigp%nbnds,Kmesh%nibz,Sigp%nsppol)
  logical,intent(out) :: keep_ur(Sigp%nbnds,Kmesh%nibz,Sigp%nsppol)
 end subroutine sigma_bksmask
end interface

end module interfaces_70_gw
!!***
