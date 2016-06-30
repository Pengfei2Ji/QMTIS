!!****m* ABINIT/interfaces_41_geometry
!! NAME
!! interfaces_41_geometry
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/41_geometry
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

module interfaces_41_geometry

 implicit none

interface
 subroutine bldgrp(msym,nogen,nsym,symafm,symrel,tnons)
  use defs_basis
  implicit none
  integer,intent(in) :: msym
  integer,intent(inout) :: nogen
  integer,intent(in) :: nsym
  integer,intent(inout) :: symafm(msym)
  integer,intent(inout) :: symrel(3,3,msym)
  real(dp),intent(inout) :: tnons(3,msym)
 end subroutine bldgrp
end interface

interface
 subroutine bonds_lgth_angles(coordn,fnameabo_app_geo,natom,ntypat,rprimd,typat,xred,znucl)
  use defs_basis
  implicit none
  integer,intent(in) :: coordn
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  character(len=*),intent(in) :: fnameabo_app_geo
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(inout) :: xred(3,natom)
  real(dp),intent(in) :: znucl(ntypat)
 end subroutine bonds_lgth_angles
end interface

interface
 subroutine chiscwrt(chi_org,disv_org,nat_org,sdisv_org,smult_org,nsh_org,chi_sc,&  
  &  disv_sc,nat_sc,smult_sc,nsh_sc,opt,prtvol) 
  use defs_basis
  implicit none
  integer,intent(in) :: nat_org
  integer,intent(in) :: nat_sc
  integer,intent(in) :: nsh_org
  integer,intent(in) :: nsh_sc
  integer,intent(in),optional :: opt
  integer,intent(in),optional :: prtvol
  real(dp),intent(in) :: chi_org(nat_org)
  real(dp),intent(out) :: chi_sc(nat_sc)
  real(dp),intent(in) :: disv_org(nat_org)
  real(dp),intent(in) :: disv_sc(nat_sc)
  real(dp),intent(in) :: sdisv_org(nsh_org)
  integer,intent(in) :: smult_org(nsh_org)
  integer,intent(in) :: smult_sc(nsh_sc)
 end subroutine chiscwrt
end interface

interface
 subroutine chkdilatmx(dilatmx,rprimd,rprimd_orig,dilatmx_errmsg)
  use defs_basis
  implicit none
  real(dp),intent(in) :: dilatmx
  character(len=500),intent(out) :: dilatmx_errmsg
  real(dp),intent(inout) :: rprimd(3,3)
  real(dp),intent(in) :: rprimd_orig(3,3)
 end subroutine chkdilatmx
end interface

interface
 subroutine chkgrp(nsym,symafm,symrel,ierr)
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: nsym
  integer,intent(in) :: symafm(nsym)
  integer,intent(in) :: symrel(3,3,nsym)
 end subroutine chkgrp
end interface

interface
 subroutine sg_multable(nsym,symafm,symrel,tnons,tnons_tol,ierr,multable,toinv)
  use defs_basis
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: nsym
  real(dp),intent(in) :: tnons_tol
  integer,optional,intent(out) :: multable(4,nsym,nsym)
  integer,intent(in) :: symafm(nsym)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
  integer,optional,intent(out) :: toinv(4,nsym)
 end subroutine sg_multable
end interface

interface
 subroutine chkorthsy(gprimd,iexit,nsym,rmet,rprimd,symrel)
  use defs_basis
  implicit none
  integer,intent(inout) :: iexit
  integer,intent(in) :: nsym
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symrel(3,3,nsym)
 end subroutine chkorthsy
end interface

interface
 subroutine chkprimit(chkprim,multi,nsym,symafm,symrel)
  implicit none
  integer,intent(in) :: chkprim
  integer,intent(out) :: multi
  integer,intent(in) :: nsym
  integer,intent(in) :: symafm(nsym)
  integer,intent(in) :: symrel(3,3,nsym)
 end subroutine chkprimit
end interface

interface
 subroutine chkrprimd(acell,rprim,rprimd,iout)
  use defs_basis
  implicit none
  integer,intent(in) :: iout
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine chkrprimd
end interface

interface
 function dbeta(cosbeta,ll,mp,mm)
  use defs_basis
  implicit none
  integer,intent(in) :: ll
  integer,intent(in) :: mm
  integer,intent(in) :: mp
  real(dp),intent(in) :: cosbeta
  real(dp) :: dbeta
 end function dbeta
end interface

interface
 function dist2(v1,v2,rprimd,option)
  use defs_basis
  implicit none
  integer,intent(in),optional :: option
  real(dp) :: dist2
  real(dp),intent(in),optional :: rprimd(3,3)
  real(dp),intent(in) :: v1(3)
  real(dp),intent(in) :: v2(3)
 end function dist2
end interface

interface
 subroutine fcart2fred(fcart,fred,rprimd,natom)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  real(dp),intent(in) :: fcart(3,natom)
  real(dp),intent(out) :: fred(3,natom)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine fcart2fred
end interface

interface
 subroutine fillcell(natom,natrd,nsym,nucdipmom,spinat,symafm,symrel,tnons,tolsym,typat,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: natrd
  integer,intent(in) :: nsym
  real(dp),intent(in) :: tolsym
  real(dp),intent(inout) :: nucdipmom(3,natom)
  real(dp),intent(inout) :: spinat(3,natom)
  integer,intent(in) :: symafm(nsym)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
  integer,intent(inout) :: typat(natom)
  real(dp),intent(inout) :: xred(3,natom)
 end subroutine fillcell
end interface

interface
 subroutine fred2fcart(favg,fcart,fred,gprimd,jellslab,natom)
  use defs_basis
  implicit none
  integer,intent(in) :: jellslab
  integer,intent(in) :: natom
  real(dp),intent(out) :: favg(3)
  real(dp),intent(out) :: fcart(3,natom)
  real(dp),intent(in) :: fred(3,natom)
  real(dp),intent(in) :: gprimd(3,3)
 end subroutine fred2fcart
end interface

interface
 subroutine gensymshub(genafm,spgroup,spgroupma,shubnikov)
  use defs_basis
  implicit none
  integer,intent(out) :: shubnikov
  integer,intent(in) :: spgroup
  integer,intent(in) :: spgroupma
  real(dp),intent(out) :: genafm(3)
 end subroutine gensymshub
end interface

interface
 subroutine gensymshub4(genafm,msym,nsym,symafm,symrel,tnons)
  use defs_basis
  implicit none
  integer,intent(in) :: msym
  integer,intent(inout) :: nsym
  real(dp),intent(in) :: genafm(3)
  integer,intent(inout) :: symafm(msym)
  integer,intent(inout) :: symrel(3,3,msym)
  real(dp),intent(inout) :: tnons(3,msym)
 end subroutine gensymshub4
end interface

interface
 subroutine gensymspgr(brvltt,msym,nsym,shubnikov,spgaxor,spgorig,&  
  &  spgroup,spgroupma,symafm,symrel,tnons)
  use defs_basis
  implicit none
  integer,intent(inout) :: brvltt
  integer,intent(in) :: msym
  integer,intent(out) :: nsym
  integer,intent(in) :: shubnikov
  integer,intent(in) :: spgaxor
  integer,intent(in) :: spgorig
  integer,intent(in) :: spgroup
  integer,intent(in) :: spgroupma
  integer,intent(inout) :: symafm(msym)
  integer,intent(inout) :: symrel(3,3,msym)
  real(dp),intent(inout) :: tnons(3,msym)
 end subroutine gensymspgr
end interface

interface
 subroutine getptgroupma(ptgroup,ptgroupha,ptgroupma)
  implicit none
  integer,intent(out) :: ptgroupma
  character(len=5),intent(in) :: ptgroup
  character(len=5),intent(in) :: ptgroupha
 end subroutine getptgroupma
end interface

interface
 subroutine getspinrot(rprimd,spinrot,symrel_conv)
  use defs_basis
  implicit none
  integer,intent(in) :: symrel_conv(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: spinrot(4)
 end subroutine getspinrot
end interface

interface
 subroutine holocell(cell_base,enforce,foundc,iholohedry,tolsym)
  use defs_basis
  implicit none
  integer,intent(in) :: enforce
  integer,intent(out) :: foundc
  integer,intent(in) :: iholohedry
  real(dp),intent(in) :: tolsym
  real(dp),intent(inout) :: cell_base(3,3)
 end subroutine holocell
end interface

interface
 subroutine ioniondist(natom,rprimd,xred,inm,option,varlist,magv,atp,prtvol)
  use defs_basis
  implicit none
  integer,intent(in),optional :: atp
  integer,intent(in) :: natom
  integer,intent(in) :: option
  integer,intent(in),optional :: prtvol
  real(dp),intent(out) :: inm(natom,natom)
  integer,intent(in),optional :: magv(natom)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in),optional :: varlist(natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine ioniondist
end interface

interface
 subroutine irreducible_set_pert(indsym,mpert,natom,nsym,pertsy,rfdir,rfpert,symq,symrec,symrel)
  implicit none
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  integer,intent(in) :: nsym
  integer,intent(in) :: rfdir(3)
  integer,intent(in) :: indsym(4,nsym,natom)
  integer,intent(out) :: pertsy(3,mpert)
  integer,intent(in) :: rfpert(mpert)
  integer,intent(in) :: symq(4,2,nsym)
  integer,intent(in) :: symrec(3,3,nsym)
  integer,intent(in) :: symrel(3,3,nsym)
 end subroutine irreducible_set_pert
end interface

interface
 subroutine littlegroup_pert(gprimd,idir,indsym,iout,ipert,natom,nsym,nsym1,&  
  &  rfmeth,symafm,symaf1,symq,symrec,symrel,symrl1,syuse,tnons,tnons1)
  use defs_basis
  implicit none
  integer,intent(in) :: idir
  integer,intent(in) :: iout
  integer,intent(in) :: ipert
  integer,intent(in) :: natom
  integer,intent(in) :: nsym
  integer,intent(out) :: nsym1
  integer,intent(in) :: rfmeth
  integer,intent(in) :: syuse
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: indsym(4,nsym,natom)
  integer,intent(out) :: symaf1(nsym)
  integer,intent(in) :: symafm(nsym)
  integer,intent(in) :: symq(4,2,nsym)
  integer,intent(in) :: symrec(3,3,nsym)
  integer,intent(in) :: symrel(3,3,nsym)
  integer,intent(out) :: symrl1(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
  real(dp),intent(out) :: tnons1(3,nsym)
 end subroutine littlegroup_pert
end interface

interface
 subroutine metric(gmet,gprimd,iout,rmet,rprimd,ucvol)
  use defs_basis
  implicit none
  integer,intent(in) :: iout
  real(dp),intent(out) :: ucvol
  real(dp),intent(out) :: gmet(3,3)
  real(dp),intent(out) :: gprimd(3,3)
  real(dp),intent(out) :: rmet(3,3)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine metric
end interface

interface
 subroutine mkeuler(rot,cosbeta,cosalp,sinalp,cosgam,singam,isn)
  use defs_basis
  implicit none
  integer,intent(out) :: isn
  real(dp),intent(out) :: cosalp
  real(dp),intent(out) :: cosbeta
  real(dp),intent(out) :: cosgam
  real(dp),intent(out) :: sinalp
  real(dp),intent(out) :: singam
  real(dp),intent(in) :: rot(3,3)
 end subroutine mkeuler
end interface

interface
 subroutine mkradim(acell,rprim,rprimd)
  use defs_basis
  implicit none
  real(dp),intent(out) :: acell(3)
  real(dp),intent(out) :: rprim(3,3)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine mkradim
end interface

interface
 subroutine mkrdim(acell,rprim,rprimd)
  use defs_basis
  implicit none
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(out) :: rprimd(3,3)
 end subroutine mkrdim
end interface

interface
 subroutine mksupercell(xred_org,magv_org,rprimd_org,nat_org,nat_sc,xred_sc,magv_sc,rprimd_sc,ext,prtvol) 
  use defs_basis
  implicit none
  integer,intent(in) :: nat_org
  integer,intent(in) :: nat_sc
  integer,intent(in),optional :: prtvol
  integer,intent(in) :: ext(3)
  integer,intent(in),optional :: magv_org(nat_org)
  real(dp),intent(out) :: magv_sc(nat_sc)
  real(dp),intent(in) :: rprimd_org(3,3)
  real(dp),intent(out) :: rprimd_sc(3,3)
  real(dp),intent(in) :: xred_org(3,nat_org)
  real(dp),intent(out) :: xred_sc(3,nat_sc)
 end subroutine mksupercell
end interface

interface
 subroutine polcart(red_ptot,pel,pel_cart,pelev,pion,pion_cart,polunit,&  
  &  ptot_cart,rprimd,ucvol,unit_out,usepaw)
  use defs_basis
  implicit none
  integer,intent(in) :: polunit
  integer,intent(in) :: unit_out
  integer,intent(in) :: usepaw
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: pel(3)
  real(dp),intent(out) :: pel_cart(3)
  real(dp),intent(in) :: pelev(3)
  real(dp),intent(in) :: pion(3)
  real(dp),intent(out) :: pion_cart(3)
  real(dp),intent(out) :: ptot_cart(3)
  real(dp),intent(in) :: red_ptot(3)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine polcart
end interface

interface
 subroutine prt_cif(brvltt, ciffname, natom, nsym, ntypat, rprimd,&  
  &  spgaxor, spgroup, spgorig, symrel, tnon, typat, xred, znucl)
  use defs_basis
  implicit none
  integer, intent(in) :: brvltt
  integer,intent(in) :: natom
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer, intent(in) :: spgaxor
  integer, intent(in) :: spgorig
  integer, intent(in) :: spgroup
  character(len=*), intent(in) :: ciffname
  real(dp), intent(in) :: rprimd(3,3)
  integer, intent(in) :: symrel(3,3,nsym)
  real(dp), intent(in) :: tnon(3,nsym)
  integer, intent(in) :: typat(natom)
  real(dp), intent(in) :: xred(3,natom)
  real(dp), intent(in) :: znucl(ntypat)
 end subroutine prt_cif
end interface

interface
 subroutine symrel2string(symrel1, tnon, string)
  use defs_basis
  implicit none
  character(len=80), intent(out) :: string
  integer, intent(in) :: symrel1(3,3)
  real(dp), intent(in) :: tnon(3)
 end subroutine symrel2string
end interface

interface
 subroutine prtposcar(fcart, fnameradix, natom, ntypat, rprimd, typat, ucvol, xred, znucl)
  use defs_basis
  implicit none
  integer, intent(in) :: natom
  integer, intent(in) :: ntypat
  character(len=fnlen), intent(in) :: fnameradix
  real(dp), intent(in) :: ucvol
  real(dp), intent(in) :: fcart(3,natom)
  real(dp), intent(in) :: rprimd(3,3)
  integer, intent(in) :: typat(natom)
  real(dp), intent(in) :: xred(3,natom)
  real(dp), intent(in) :: znucl(ntypat)
 end subroutine prtposcar
end interface

interface
 subroutine prtspgroup(bravais,genafm,iout,jdtset,ptgroupma,spgroup)
  use defs_basis
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: jdtset
  integer,intent(in) :: ptgroupma
  integer,intent(in) :: spgroup
  integer,intent(in) :: bravais(11)
  real(dp),intent(inout) :: genafm(3)
 end subroutine prtspgroup
end interface

interface
 subroutine ptgmadata(ptgroupma,ptgrpmasb)
  implicit none
  integer,intent(in) :: ptgroupma
  character(len=10),intent(out) :: ptgrpmasb
 end subroutine ptgmadata
end interface

interface
 subroutine randomcellpos(natom,npsp,ntypat,random_atpos,ratsph,rprim,rprimd,typat,xred,znucl,acell)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: npsp
  integer,intent(in) :: ntypat
  integer,intent(in) :: random_atpos
  real(dp), intent(inout) :: acell(3)
  real(dp),intent(in) :: ratsph(ntypat)
  real(dp), intent(inout) :: rprim(3,3)
  real(dp), intent(inout) :: rprimd(3,3)
  integer, intent(in) :: typat(natom)
  real(dp), intent(inout) :: xred(3,natom)
  real(dp), intent(in) :: znucl(npsp)
 end subroutine randomcellpos
end interface

interface
 subroutine remove_inversion(nsym,symrel,tnons,nsym_out,symrel_out,tnons_out,pinv)
  use defs_basis
  implicit none
  integer,intent(in) :: nsym
  integer,intent(out) :: nsym_out
  integer,intent(out) :: pinv
  integer,pointer :: symrel_out(:,:,:)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
  real(dp),pointer :: tnons_out(:,:)
 end subroutine remove_inversion
end interface

interface
 subroutine shellstruct(xred,rprimd,natom,magv,distv,smult,sdisv,nsh,atp,prtvol) 
  use defs_basis
  implicit none
  integer,intent(in),optional :: atp
  integer,intent(in) :: natom
  integer,intent(out) :: nsh
  integer,intent(in),optional :: prtvol
  real(dp),intent(out) :: distv(natom)
  integer,intent(in),optional :: magv(natom)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: sdisv(natom)
  integer,intent(out) :: smult(natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine shellstruct
end interface

interface
 subroutine smallprim(metmin,minim,rprimd)
  use defs_basis
  implicit none
  real(dp),intent(out) :: metmin(3,3)
  real(dp),intent(out) :: minim(3,3)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine smallprim
end interface

interface
 subroutine spgdata(brvsb,intsb,intsbl,ptintsb,ptschsb,&  
  &  schsb,spgaxor,spgroup,sporder,spgorig)
  implicit none
  integer,intent(in) :: spgaxor
  integer,intent(in) :: spgorig
  integer,intent(in) :: spgroup
  integer,intent(out) :: sporder
  character(len=1),intent(out) :: brvsb
  character(len=15),intent(out) :: intsb
  character(len=35),intent(out) :: intsbl
  character(len=15),intent(out) :: ptintsb
  character(len=15),intent(out) :: ptschsb
  character(len=15),intent(out) :: schsb
 end subroutine spgdata
end interface

interface
 subroutine strconv(frac,gprimd,cart)
  use defs_basis
  implicit none
  real(dp),intent(inout) :: cart(6)
  real(dp),intent(in) :: frac(6)
  real(dp),intent(in) :: gprimd(3,3)
 end subroutine strconv
end interface

interface
 subroutine stresssym(gprimd,nsym,stress,sym)
  use defs_basis
  implicit none
  integer,intent(in) :: nsym
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(inout) :: stress(6)
  integer,intent(in) :: sym(3,3,nsym)
 end subroutine stresssym
end interface

interface
 subroutine symanal(bravais,chkprim,genafm,msym,nsym,ptgroupma,rprimd,spgroup,symafm,symrel,tnons,tolsym)
  use defs_basis
  implicit none
  integer,intent(in) :: chkprim
  integer,intent(in) :: msym
  integer,intent(in) :: nsym
  integer,intent(out) :: ptgroupma
  integer,intent(out) :: spgroup
  real(dp),intent(in) :: tolsym
  integer,intent(out) :: bravais(11)
  real(dp),intent(out) :: genafm(3)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symafm(msym)
  integer,intent(in) :: symrel(3,3,msym)
  real(dp),intent(inout) :: tnons(3,msym)
 end subroutine symanal
end interface

interface
 subroutine symatm(indsym,natom,nsym,symrec,tnons,tolsym,typat,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nsym
  real(dp), intent(in) :: tolsym
  integer,intent(out) :: indsym(4,nsym,natom)
  integer,intent(in) :: symrec(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine symatm
end interface

interface
 subroutine symaxes(center,iholohedry,isym,isymrelconv,label,ordersym,tnons_order,trialt,type_axis)
  use defs_basis
  implicit none
  integer,intent(in) :: center
  integer,intent(in) :: iholohedry
  integer,intent(in) :: isym
  integer,intent(in) :: ordersym
  integer,intent(in) :: tnons_order
  integer,intent(out) :: type_axis
  character(len=128),intent(out) :: label
  integer,intent(in) :: isymrelconv(3,3)
  real(dp),intent(in) :: trialt(3)
 end subroutine symaxes
end interface

interface
 subroutine symbrav(bravais,msym,nsym,ptgroup,rprimd,symrel,tolsym,axis)
  use defs_basis
  implicit none
  integer,intent(in) :: msym
  integer,intent(in) :: nsym
  character(len=5),intent(out) :: ptgroup
  real(dp),intent(in) :: tolsym
  integer,optional,intent(out) :: axis(3)
  integer,intent(out) :: bravais(11)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symrel(3,3,msym)
 end subroutine symbrav
end interface

interface
 subroutine symcharac(center, determinant, iholohedry, isym, label, symrel, tnons, type_axis)
  use defs_basis
  implicit none
  integer, intent(in) :: center
  integer, intent(in) :: determinant
  integer, intent(in) :: iholohedry
  integer, intent(in) :: isym
  integer, intent(out) :: type_axis
  character(len = 128), intent(out) :: label
  integer,intent(in) :: symrel(3,3)
  real(dp),intent(in) :: tnons(3)
 end subroutine symcharac
end interface

interface
 subroutine symchk(difmin,eatom,natom,tratom,transl,trtypat,typat,xred)
  use defs_basis
  implicit none
  integer,intent(out) :: eatom
  integer,intent(in) :: natom
  integer,intent(in) :: trtypat
  integer,intent(out) :: transl(3)
  real(dp),intent(out) :: difmin(3)
  real(dp),intent(in) :: tratom(3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine symchk
end interface

interface
 subroutine symdet(determinant,nsym,sym)
  implicit none
  integer,intent(in) :: nsym
  integer,intent(out) :: determinant(nsym)
  integer,intent(in) :: sym(3,3,nsym)
 end subroutine symdet
end interface

interface
 subroutine symfind(berryopt,efield,gprimd,jellslab,msym,natom,noncoll,nptsym,nsym,&  
  &  ptsymrel,spinat,symafm,symrel,tnons,tolsym,typat,use_inversion,xred,&  
  &  nucdipmom)
  use defs_basis
  implicit none
  integer,intent(in) :: berryopt
  integer,intent(in) :: jellslab
  integer,intent(in) :: msym
  integer,intent(in) :: natom
  integer,intent(in) :: noncoll
  integer,intent(in) :: nptsym
  integer,intent(out) :: nsym
  integer,intent(in) :: use_inversion
  real(dp),intent(in) :: tolsym
  real(dp),intent(in) :: efield(3)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),optional :: nucdipmom(3,natom)
  integer,intent(in) :: ptsymrel(3,3,msym)
  real(dp),intent(in) :: spinat(3,natom)
  integer,intent(inout) :: symafm(msym)
  integer,intent(inout) :: symrel(3,3,msym)
  real(dp),intent(inout) :: tnons(3,msym)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine symfind
end interface

interface
 subroutine symkpt(chksymbreak,gmet,indkpt1,iout,kptns,nkpt,nkpt1,nsym1,&  
  &  symrc1,timrev,wtk,wtk_folded)
  use defs_basis
  implicit none
  integer,intent(in) :: chksymbreak
  integer,intent(in) :: iout
  integer,intent(in) :: nkpt
  integer,intent(out) :: nkpt1
  integer,intent(in) :: nsym1
  integer,intent(in) :: timrev
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(inout) :: indkpt1(nkpt)
  real(dp),intent(in) :: kptns(3,nkpt)
  integer,intent(in) :: symrc1(3,3,nsym1)
  real(dp),intent(in) :: wtk(nkpt)
  real(dp),intent(out) :: wtk_folded(nkpt)
 end subroutine symkpt
end interface

interface
 subroutine symlatt(bravais,msym,nptsym,ptsymrel,rprimd,tolsym)
  use defs_basis
  implicit none
  integer,intent(in) :: msym
  integer,intent(out) :: nptsym
  real(dp),intent(in) :: tolsym
  integer,intent(out) :: bravais(11)
  integer,intent(out) :: ptsymrel(3,3,msym)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine symlatt
end interface

interface
 subroutine symlist_bcc(additional_info,nsym,n_axes,spgroup)
  implicit none
  integer,intent(in) :: additional_info
  integer,intent(in) :: nsym
  integer,intent(out) :: spgroup
  integer,intent(in) :: n_axes(31)
 end subroutine symlist_bcc
end interface

interface
 subroutine symlist_fcc(nsym,n_axes,spgroup)
  implicit none
  integer,intent(in) :: nsym
  integer,intent(out) :: spgroup
  integer,intent(in) :: n_axes(31)
 end subroutine symlist_fcc
end interface

interface
 subroutine symlist_others(brvltt,nsym,n_axes,spgroup)
  implicit none
  integer,intent(in) :: brvltt
  integer,intent(in) :: nsym
  integer,intent(out) :: spgroup
  integer,intent(in) :: n_axes(31)
 end subroutine symlist_others
end interface

interface
 subroutine symlist_prim(additional_info,nsym,n_axes,spgroup)
  implicit none
  integer,intent(in) :: additional_info
  integer,intent(in) :: nsym
  integer,intent(out) :: spgroup
  integer,intent(in) :: n_axes(31)
 end subroutine symlist_prim
end interface

interface
 subroutine symmetrize_rprimd(bravais,nsym,rprimd,symrel,tolsym)
  use defs_basis
  implicit none
  integer,intent(in) :: nsym
  real(dp) :: tolsym
  integer,intent(in) :: bravais(11)
  real(dp),intent(inout) :: rprimd(3,3)
  integer,intent(in) :: symrel(3,3,nsym)
 end subroutine symmetrize_rprimd
end interface

interface
 subroutine symmetrize_xred(indsym,natom,nsym,symrel,tnons,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nsym
  integer,intent(in) :: indsym(4,nsym,natom)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
  real(dp),intent(inout) :: xred(3,natom)
 end subroutine symmetrize_xred
end interface

interface
 subroutine symplanes(center,iholohedry,isym,isymrelconv,itnonsconv,label,type_axis)
  use defs_basis
  implicit none
  integer,intent(in) :: center
  integer,intent(in) :: iholohedry
  integer,intent(in) :: isym
  integer,intent(out) :: type_axis
  character(len = 128), intent(out) :: label
  integer,intent(in) :: isymrelconv(3,3)
  real(dp),intent(in) :: itnonsconv(3)
 end subroutine symplanes
end interface

interface
 subroutine symptgroup(iholohedry,nsym,ptgroup,symrel)
  implicit none
  integer,intent(out) :: iholohedry
  integer,intent(in) :: nsym
  character(len=5),intent(out) :: ptgroup
  integer,intent(in) :: symrel(3,3,nsym)
 end subroutine symptgroup
end interface

interface
 subroutine symredcart(aprim,bprim,symcart,symred)
  use defs_basis
  implicit none
  integer,intent(in) :: symred(3,3)
  real(dp),intent(in) :: aprim(3,3)
  real(dp),intent(in) :: bprim(3,3)
  real(dp),intent(out) :: symcart(3,3)
 end subroutine symredcart
end interface

interface
 subroutine symrelrot(nsym,rprimd,rprimd_new,symrel,tolsym)
  use defs_basis
  implicit none
  integer,intent(in) :: nsym
  real(dp),intent(in) :: tolsym
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: rprimd_new(3,3)
  integer,intent(inout) :: symrel(3,3,nsym)
 end subroutine symrelrot
end interface

interface
 subroutine symsgcube(msym,nsym,shubnikov,spgaxor,spgorig,spgroup,&  
  &  spgroupma,symafm,symrel,tnons)
  use defs_basis
  implicit none
  integer,intent(in) :: msym
  integer,intent(inout) :: nsym
  integer,intent(in) :: shubnikov
  integer,intent(in) :: spgaxor
  integer,intent(in) :: spgorig
  integer,intent(in) :: spgroup
  integer,intent(in) :: spgroupma
  integer,intent(inout) :: symafm(msym)
  integer,intent(inout) :: symrel(3,3,msym)
  real(dp),intent(out) :: tnons(3,msym)
 end subroutine symsgcube
end interface

interface
 subroutine symsghexa(brvltt,msym,nsym,shubnikov,spgaxor,spgorig,spgroup,&  
  &  spgroupma,symafm,symrel,tnons)
  use defs_basis
  implicit none
  integer,intent(inout) :: brvltt
  integer,intent(in) :: msym
  integer,intent(in) :: nsym
  integer,intent(in) :: shubnikov
  integer,intent(in) :: spgaxor
  integer,intent(in) :: spgorig
  integer,intent(in) :: spgroup
  integer,intent(in) :: spgroupma
  integer,intent(inout) :: symafm(msym)
  integer,intent(inout) :: symrel(3,3,msym)
  real(dp),intent(inout) :: tnons(3,msym)
 end subroutine symsghexa
end interface

interface
 subroutine symsgmono(brvltt,msym,nsym,shubnikov,spgaxor,spgorig,spgroup,&  
  &  spgroupma,symafm,symrel,tnons)
  use defs_basis
  implicit none
  integer,intent(inout) :: brvltt
  integer,intent(in) :: msym
  integer,intent(in) :: nsym
  integer,intent(in) :: shubnikov
  integer,intent(in) :: spgaxor
  integer,intent(in) :: spgorig
  integer,intent(in) :: spgroup
  integer,intent(in) :: spgroupma
  integer,intent(inout) :: symafm(msym)
  integer,intent(inout) :: symrel(3,3,msym)
  real(dp),intent(inout) :: tnons(3,msym)
 end subroutine symsgmono
end interface

interface
 subroutine symsgortho(msym,nsym,shubnikov,spgaxor,spgorig,spgroup,&  
  &  spgroupma,symafm,symrel,tnons)
  use defs_basis
  implicit none
  integer,intent(in) :: msym
  integer,intent(in) :: nsym
  integer,intent(in) :: shubnikov
  integer,intent(in) :: spgaxor
  integer,intent(in) :: spgorig
  integer,intent(in) :: spgroup
  integer,intent(in) :: spgroupma
  integer,intent(inout) :: symafm(msym)
  integer,intent(inout) :: symrel(3,3,msym)
  real(dp),intent(inout) :: tnons(3,msym)
 end subroutine symsgortho
end interface

interface
 subroutine symsgtetra(msym,nsym,shubnikov,spgaxor,spgorig,spgroup,&  
  &  spgroupma,symafm,symrel,tnons)
  use defs_basis
  implicit none
  integer,intent(in) :: msym
  integer,intent(in) :: nsym
  integer,intent(in) :: shubnikov
  integer,intent(in) :: spgaxor
  integer,intent(in) :: spgorig
  integer,intent(in) :: spgroup
  integer,intent(in) :: spgroupma
  integer,intent(inout) :: symafm(msym)
  integer,intent(inout) :: symrel(3,3,msym)
  real(dp),intent(inout) :: tnons(3,msym)
 end subroutine symsgtetra
end interface

interface
 subroutine symspgr(bravais,nsym,spgroup,symrel,tnons,tolsym)
  use defs_basis
  implicit none
  integer,intent(in) :: nsym
  integer,intent(out) :: spgroup
  real(dp),intent(in) :: tolsym
  integer,intent(in) :: bravais(11)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(inout) :: tnons(3,nsym)
 end subroutine symspgr
end interface

interface
 subroutine xcart2xred(natom,rprimd,xcart,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: xcart(3,natom)
  real(dp),intent(out) :: xred(3,natom)
 end subroutine xcart2xred
end interface

interface
 subroutine xred2xcart(natom,rprimd,xcart,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: xcart(3,natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine xred2xcart
end interface

end module interfaces_41_geometry
!!***
