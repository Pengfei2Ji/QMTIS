!!****m* ABINIT/interfaces_68_rsprc
!! NAME
!! interfaces_68_rsprc
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/68_rsprc
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

module interfaces_68_rsprc

 implicit none

interface
 subroutine moddiel_csrb(dielar,dtset,gprimd,mpi_enreg,rdiemac,rhor_in)
  use defs_basis
  use defs_abitypes
  implicit none
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: dielar(7)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(out) :: rdiemac(dtset%nfft,dtset%nspden)
  real(dp),intent(in) :: rhor_in(dtset%nfft,dtset%nspden)
 end subroutine moddiel_csrb
end interface

interface
 subroutine newrho(atindx,dbl_nnsclo,dielar,dielinv,dielstrt,dtn_pc,dtset,etotal,fcart,ffttomix,&  
  &  gmet,grhf,gsqcut,initialized,ispmix,istep,kg_diel,kxc,mgfft,mix,mixtofft,&  
  &  moved_atm_inside,mpi_enreg,my_natom,nattyp,nfft,&  
  &  nfftmix,nfftmix_per_nfft,ngfft,ngfftmix,nkxc,npawmix,npwdiel,&  
  &  nresid,ntypat,n1xccc,pawrhoij,pawtab,&  
  &  ph1d,psps,rhog,rhor,rprimd,susmat,usepaw,vtrial,wvl,wvl_den,xred)
  use m_pawrhoij
  use defs_abitypes
  use m_pawtab
  use defs_basis
  use defs_datatypes
  use m_ab7_mixing
  use defs_wvltypes
  implicit none
  integer,intent(inout) :: dbl_nnsclo
  integer,intent(in) :: dielstrt
  integer,intent(in) :: initialized
  integer,intent(in) :: ispmix
  integer,intent(in) :: istep
  integer,intent(in) :: mgfft
  integer,intent(in) :: moved_atm_inside
  integer,intent(in) :: my_natom
  integer,intent(in) :: n1xccc
  integer,intent(in) :: nfft
  integer,intent(in) :: nfftmix
  integer,intent(in) :: nfftmix_per_nfft
  integer,intent(in) :: nkxc
  integer,intent(in) :: npawmix
  integer,intent(in) :: npwdiel
  integer,intent(in) :: ntypat
  integer,intent(in) :: usepaw
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: etotal
  real(dp),intent(in) :: gsqcut
  type(ab7_mixing_object), intent(inout) :: mix
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  type(wvl_internal_type), intent(in) :: wvl
  type(wvl_denspot_type), intent(inout) :: wvl_den
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: ngfftmix(18)
  integer,intent(in) :: atindx(dtset%natom)
  real(dp),intent(in) :: dielar(7)
  real(dp),intent(inout) :: dielinv(2,npwdiel,dtset%nspden,npwdiel,dtset%nspden)
  real(dp),intent(inout), target :: dtn_pc(3,dtset%natom)
  real(dp),intent(in) :: fcart(3,dtset%natom)
  integer,intent(in) :: ffttomix(nfft*(nfftmix_per_nfft))
  real(dp), intent(inout) :: gmet(3,3)
  real(dp),intent(in) :: grhf(3,dtset%natom)
  integer,intent(in) :: kg_diel(3,npwdiel)
  real(dp),intent(inout) :: kxc(nfft,nkxc)
  integer,intent(in) :: mixtofft(nfftmix*nfftmix_per_nfft)
  integer,intent(in) :: nattyp(ntypat)
  real(dp),intent(inout) :: nresid(nfft,dtset%nspden)
  type(pawrhoij_type),intent(inout) :: pawrhoij(my_natom*psps%usepaw)
  type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)
  real(dp),intent(inout) :: ph1d(2,3*(2*mgfft+1)*dtset%natom)
  real(dp),intent(inout) :: rhog(2,nfft)
  real(dp),intent(inout) :: rhor(nfft,dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: susmat(2,npwdiel,dtset%nspden,npwdiel,dtset%nspden)
  real(dp),intent(in), target :: vtrial(nfft,dtset%nspden)
  real(dp), intent(inout), target :: xred(3,dtset%natom)
 end subroutine newrho
end interface

interface
 subroutine newvtr(atindx,dbl_nnsclo,dielar,dielinv,dielstrt,&  
  &  dtn_pc,dtset,etotal,fcart,ffttomix,&  
  &  gmet,grhf,gsqcut,&  
  &  initialized,ispmix,&  
  &  istep,&  
  &  kg_diel,kxc,mgfft,mix,mixtofft,&  
  &  moved_atm_inside,mpi_enreg,my_natom,nattyp,nfft,nfftmix,&  
  &  ngfft,ngfftmix,nkxc,npawmix,npwdiel,&  
  &  nstep,ntypat,n1xccc,&  
  &  pawrhoij,&  
  &  ph1d,&  
  &  psps,rhor,rprimd,susmat,usepaw,&  
  &  vhartr,vnew_mean,vpsp,vresid,&  
  &  vtrial,vxc,xred,&  
  &  nfftf,&  
  &  pawtab,&  
  &  rhog,&  
  &  wvl)
  use m_pawrhoij
  use defs_abitypes
  use m_pawtab
  use defs_basis
  use defs_datatypes
  use m_ab7_mixing
  use defs_wvltypes
  implicit none
  integer,intent(inout) :: dbl_nnsclo
  integer,intent(in) :: dielstrt
  integer,intent(in) :: initialized
  integer,intent(in) :: ispmix
  integer,intent(in) :: istep
  integer,intent(in) :: mgfft
  integer,intent(in) :: moved_atm_inside
  integer,intent(in) :: my_natom
  integer,intent(in) :: n1xccc
  integer,intent(in) :: nfft
  integer,intent(in) :: nfftf
  integer,intent(in) :: nfftmix
  integer,intent(in) :: nkxc
  integer,intent(in) :: npawmix
  integer,intent(in) :: npwdiel
  integer,intent(in) :: nstep
  integer,intent(in) :: ntypat
  integer,intent(in) :: usepaw
  type(dataset_type),intent(inout) :: dtset
  real(dp),intent(in) :: etotal
  real(dp),intent(in) :: gsqcut
  type(ab7_mixing_object),intent(inout) :: mix
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  type(wvl_data), intent(inout) :: wvl
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: ngfftmix(18)
  integer,intent(in) :: atindx(dtset%natom)
  real(dp),intent(in) :: dielar(7)
  real(dp),intent(inout) :: dielinv(2,npwdiel,dtset%nspden,npwdiel,dtset%nspden)
  real(dp),intent(inout), target :: dtn_pc(3,dtset%natom)
  real(dp),intent(in) :: fcart(3,dtset%natom)
  integer,intent(in) :: ffttomix(nfft*(1-nfftmix/nfft))
  real(dp),intent(inout) :: gmet(3,3)
  real(dp),intent(in) :: grhf(3,dtset%natom)
  integer,intent(in) :: kg_diel(3,npwdiel)
  real(dp),intent(inout) :: kxc(nfft,nkxc)
  integer,intent(in) :: mixtofft(nfftmix*(1-nfftmix/nfft))
  integer,intent(in) :: nattyp(ntypat)
  type(pawrhoij_type),intent(inout) :: pawrhoij(my_natom*usepaw)
  type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)
  real(dp),intent(inout) :: ph1d(2,3*(2*mgfft+1)*dtset%natom)
  real(dp),intent(inout) :: rhog(2,nfftf)
  real(dp),intent(inout), target :: rhor(nfft,dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: susmat(2,npwdiel,dtset%nspden,npwdiel,dtset%nspden)
  real(dp),intent(in) :: vhartr(nfft)
  real(dp),intent(in) :: vnew_mean(dtset%nspden)
  real(dp),intent(inout) :: vpsp(nfft)
  real(dp),intent(inout) :: vresid(nfft,dtset%nspden)
  real(dp),intent(inout) :: vtrial(nfft,dtset%nspden)
  real(dp),intent(in) :: vxc(nfft,dtset%nspden)
  real(dp),intent(inout), target :: xred(3,dtset%natom)
 end subroutine newvtr
end interface

interface
 subroutine prcref(atindx,dielar,dielinv,&  
  &  dielstrt,dtn_pc,dtset,etotal,fcart,ffttomix,gmet,gsqcut,&  
  &  istep,kg_diel,kxc,&  
  &  mgfft,moved_atm_inside,mpi_enreg,my_natom,&  
  &  nattyp,nfft,nfftprc,ngfft,ngfftprc,nkxc,npawmix,npwdiel,ntypat,n1xccc,&  
  &  optreal,optres,pawrhoij,pawtab,ph1d,psps,rhog,rhoijrespc,rhor,rprimd,&  
  &  susmat,vhartr,vpsp,vresid,vrespc,vxc,wvl,wvl_den,xred)
  use m_pawrhoij
  use defs_abitypes
  use m_pawtab
  use defs_basis
  use defs_datatypes
  use defs_wvltypes
  implicit none
  integer,intent(in) :: dielstrt
  integer,intent(in) :: istep
  integer,intent(in) :: mgfft
  integer,intent(in) :: moved_atm_inside
  integer,intent(in) :: my_natom
  integer,intent(in) :: n1xccc
  integer,intent(in) :: nfft
  integer,intent(in) :: nfftprc
  integer,intent(in) :: nkxc
  integer,intent(in) :: npawmix
  integer,intent(in) :: npwdiel
  integer,intent(in) :: ntypat
  integer,intent(in) :: optreal
  integer,intent(in) :: optres
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: etotal
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  type(wvl_internal_type), intent(in) :: wvl
  type(wvl_denspot_type), intent(inout) :: wvl_den
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: ngfftprc(18)
  integer,intent(in) :: atindx(dtset%natom)
  real(dp),intent(in) :: dielar(7)
  real(dp),intent(inout) :: dielinv(2,npwdiel,dtset%nspden,npwdiel,dtset%nspden)
  real(dp),intent(out) :: dtn_pc(3,dtset%natom)
  real(dp),intent(in) :: fcart(3,dtset%natom)
  integer,intent(in) :: ffttomix(nfft*(1-nfftprc/nfft))
  real(dp),intent(inout) :: gmet(3,3)
  integer,intent(in) :: kg_diel(3,npwdiel)
  real(dp),intent(inout) :: kxc(nfft,nkxc)
  integer,intent(in) :: nattyp(ntypat)
  type(pawrhoij_type),intent(inout) :: pawrhoij(my_natom*psps%usepaw)
  type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)
  real(dp),intent(inout) :: ph1d(2,3*(2*mgfft+1)*dtset%natom)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(out) :: rhoijrespc(npawmix)
  real(dp),intent(in) :: rhor(nfft,dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: susmat(2,npwdiel,dtset%nspden,npwdiel,dtset%nspden)
  real(dp),intent(in) :: vhartr(nfft)
  real(dp),intent(inout) :: vpsp(nfft)
  real(dp),intent(in) :: vresid(nfftprc*optreal,dtset%nspden)
  real(dp),intent(out) :: vrespc(nfftprc*optreal,dtset%nspden)
  real(dp),intent(in) :: vxc(nfft,dtset%nspden)
  real(dp),intent(inout) :: xred(3,dtset%natom)
 end subroutine prcref
end interface

interface
 subroutine prcref_PMA(atindx,dielar,dielinv,&  
  &  dielstrt,dtn_pc,dtset,fcart,ffttomix,gmet,gsqcut,&  
  &  istep,kg_diel,kxc,&  
  &  mgfft,moved_atm_inside,mpi_enreg,my_natom,&  
  &  nattyp,nfft,nfftprc,ngfft,ngfftprc,nkxc,npawmix,npwdiel,ntypat,n1xccc,&  
  &  optreal,optres,pawrhoij,ph1d,psps,rhog, rhoijrespc,rhor,rprimd,&  
  &  susmat,vhartr,vpsp,vresid,vrespc,vxc,xred,&  
  &  etotal,pawtab,wvl)
  use m_pawrhoij
  use defs_abitypes
  use m_pawtab
  use defs_basis
  use defs_datatypes
  use defs_wvltypes
  implicit none
  integer,intent(in) :: dielstrt
  integer,intent(in) :: istep
  integer,intent(in) :: mgfft
  integer,intent(in) :: moved_atm_inside
  integer,intent(in) :: my_natom
  integer,intent(in) :: n1xccc
  integer,intent(in) :: nfft
  integer,intent(in) :: nfftprc
  integer,intent(in) :: nkxc
  integer,intent(in) :: npawmix
  integer,intent(in) :: npwdiel
  integer,intent(in) :: ntypat
  integer,intent(in) :: optreal
  integer,intent(in) :: optres
  type(dataset_type),intent(inout) :: dtset
  real(dp),intent(in) :: etotal
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  type(wvl_data), intent(inout) :: wvl
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: ngfftprc(18)
  integer,intent(in) :: atindx(dtset%natom)
  real(dp),intent(in) :: dielar(7)
  real(dp),intent(inout) :: dielinv(2,npwdiel,dtset%nspden,npwdiel,dtset%nspden)
  real(dp),intent(out) :: dtn_pc(3,dtset%natom)
  real(dp),intent(in) :: fcart(3,dtset%natom)
  integer,intent(in) :: ffttomix(nfft*(1-nfftprc/nfft))
  real(dp),intent(inout) :: gmet(3,3)
  integer,intent(in) :: kg_diel(3,npwdiel)
  real(dp),intent(inout) :: kxc(nfft,nkxc)
  integer,intent(in) :: nattyp(ntypat)
  type(pawrhoij_type),intent(inout) :: pawrhoij(my_natom*psps%usepaw)
  type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)
  real(dp),intent(inout) :: ph1d(2,3*(2*mgfft+1)*dtset%natom)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(out) :: rhoijrespc(npawmix)
  real(dp),intent(in) :: rhor(nfft,dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: susmat(2,npwdiel,dtset%nspden,npwdiel,dtset%nspden)
  real(dp),intent(in) :: vhartr(nfft)
  real(dp),intent(inout) :: vpsp(nfft)
  real(dp),intent(in) :: vresid(nfftprc*optreal,dtset%nspden)
  real(dp),intent(out) :: vrespc(nfftprc*optreal,dtset%nspden)
  real(dp),intent(in) :: vxc(nfft,dtset%nspden)
  real(dp),intent(inout) :: xred(3,dtset%natom)
 end subroutine prcref_PMA
end interface

interface
 subroutine prcrskerker1(dtset,mpi_enreg,nfft,nspden,ngfft,dielar,etotal,gprimd,vresid,vrespc,base)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  type(dataset_type),intent(in) :: dtset
  real(dp) :: etotal
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: base(nfft)
  real(dp),intent(in) :: dielar(7)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: vresid(nfft,nspden)
  real(dp),intent(out) :: vrespc(nfft,nspden)
 end subroutine prcrskerker1
end interface

interface
 subroutine prcrskerker2(dtset,nfft,nspden,ngfft,dielar,gprimd,rprimd,vresid,vrespc,natom,xred,mpi_enreg,ucvol)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(in) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: dielar(7)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: vresid(nfft,nspden)
  real(dp),intent(out) :: vrespc(nfft,nspden)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine prcrskerker2
end interface

interface
 subroutine wvl_prcref(dielar,iprcel,my_natom,nfftprc,npawmix,nspden,pawrhoij,&  
  &  rhoijrespc,usepaw,vresid,vrespc)
  use defs_basis
  use m_pawrhoij
  implicit none
  integer , intent(in) :: iprcel
  integer , intent(in) :: my_natom
  integer , intent(in) :: nfftprc
  integer , intent(in) :: npawmix
  integer , intent(in) :: nspden
  integer , intent(in) :: usepaw
  real(dp), intent(in) :: dielar(7)
  type(pawrhoij_type),intent(inout) :: pawrhoij(my_natom*usepaw)
  real(dp),intent(out) :: rhoijrespc(npawmix)
  real(dp), intent(in) :: vresid(nfftprc,nspden)
  real(dp),intent(out) :: vrespc(nfftprc,nspden)
 end subroutine wvl_prcref
end interface

end module interfaces_68_rsprc
!!***
