!!****m* ABINIT/interfaces_45_geomoptim
!! NAME
!! interfaces_45_geomoptim
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/45_geomoptim
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

module interfaces_45_geomoptim

 implicit none

interface
 subroutine calc_b_matrix(deloc,natom,rprimd,xcart,b_matrix)
  use defs_basis
  use m_abimover
  implicit none
  integer,intent(in) :: natom
  type(delocint),intent(inout) :: deloc
  real(dp),intent(out) :: b_matrix(deloc%ninternal,3*natom)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: xcart(3,natom)
 end subroutine calc_b_matrix
end interface

interface
 subroutine dbond_length_d1(r1,r2,bb)
  use defs_basis
  implicit none
  real(dp),intent(out) :: bb(3)
  real(dp),intent(in) :: r1(3)
  real(dp),intent(in) :: r2(3)
 end subroutine dbond_length_d1
end interface

interface
 subroutine dang_d1(r1,r2,r3,bb)
  use defs_basis
  implicit none
  real(dp),intent(out) :: bb(3)
  real(dp),intent(in) :: r1(3)
  real(dp),intent(in) :: r2(3)
  real(dp),intent(in) :: r3(3)
 end subroutine dang_d1
end interface

interface
 subroutine dang_d2(r1,r2,r3,bb)
  use defs_basis
  implicit none
  real(dp),intent(out) :: bb(3)
  real(dp),intent(in) :: r1(3)
  real(dp),intent(in) :: r2(3)
  real(dp),intent(in) :: r3(3)
 end subroutine dang_d2
end interface

interface
 subroutine ddihedral_d1(r1,r2,r3,r4,bb)
  use defs_basis
  implicit none
  real(dp),intent(out) :: bb(3)
  real(dp),intent(in) :: r1(3)
  real(dp),intent(in) :: r2(3)
  real(dp),intent(in) :: r3(3)
  real(dp),intent(in) :: r4(3)
 end subroutine ddihedral_d1
end interface

interface
 subroutine ddihedral_d2(r1,r2,r3,r4,bb)
  use defs_basis
  implicit none
  real(dp),intent(out) :: bb(3)
  real(dp),intent(in) :: r1(3)
  real(dp),intent(in) :: r2(3)
  real(dp),intent(in) :: r3(3)
  real(dp),intent(in) :: r4(3)
 end subroutine ddihedral_d2
end interface

interface
 subroutine deloc2xcart(deloc,natom,rprimd,xcart,deloc_int,btinv,u_matrix)
  use defs_basis
  use m_abimover
  implicit none
  integer,intent(in) :: natom
  type(delocint),intent(inout) :: deloc
  real(dp),intent(out) :: btinv(3*(natom-1),3*natom)
  real(dp),intent(in) :: deloc_int(3*(natom-1))
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(inout) :: u_matrix(deloc%ninternal,3*(natom-1))
  real(dp),intent(inout) :: xcart(3,natom)
 end subroutine deloc2xcart
end interface

interface
 subroutine fred2fdeloc(btinv,deloc_force,fred,natom,gprimd)
  use defs_basis
  implicit none
  integer, intent(in) :: natom
  real(dp),intent(in) :: btinv(3*(natom-1),3*natom)
  real(dp),intent(out) :: deloc_force(3*(natom-1))
  real(dp),intent(in) :: fred(3,natom)
  real(dp),intent(in) :: gprimd(3,3)
 end subroutine fred2fdeloc
end interface

interface
 subroutine isotemp(amass,dtion,ekin,iatfix,ktemp,mttk_vars,natom,nnos,qmass,vel)
  use defs_basis
  use m_abimover
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nnos
  real(dp),intent(in) :: dtion
  real(dp),intent(out) :: ekin
  real(dp),intent(in) :: ktemp
  type(mttk_type) :: mttk_vars
  integer,intent(in) :: iatfix(:,:)
  real(dp),intent(in) :: amass(natom)
  real(dp),intent(in) :: qmass(:)
  real(dp),intent(inout) :: vel(3,natom)
 end subroutine isotemp
end interface

interface
 subroutine isopress(amass,bmass,dtion,ekin,iatfix,ktemp,mttk_vars,natom,nnos,qmass,&  
  &  strten,strtarget,ucvol,vel,vlogv)
  use defs_basis
  use m_abimover
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nnos
  real(dp),intent(in) :: bmass
  real(dp),intent(in) :: dtion
  real(dp),intent(out) :: ekin
  real(dp),intent(in) :: ktemp
  type(mttk_type) :: mttk_vars
  real(dp),intent(in) :: ucvol
  real(dp),intent(inout) :: vlogv
  integer,intent(in) :: iatfix(:,:)
  real(dp),intent(in) :: amass(natom)
  real(dp),intent(in) :: qmass(:)
  real(dp),intent(in) :: strtarget(6)
  real(dp),intent(in) :: strten(6)
  real(dp),intent(inout) :: vel(3,natom)
 end subroutine isopress
end interface

interface
 subroutine isostress(amass,bmass,dtion,ekin,iatfix,ktemp,mttk_vars,natom,nnos,&  
  &  qmass,strten,strtarget,ucvol,vel)
  use defs_basis
  use m_abimover
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nnos
  real(dp),intent(in) :: bmass
  real(dp),intent(in) :: dtion
  real(dp),intent(out) :: ekin
  real(dp),intent(in) :: ktemp
  type(mttk_type) :: mttk_vars
  real(dp),intent(in) :: ucvol
  integer, intent(in) :: iatfix(:,:)
  real(dp),intent(in) :: amass(natom)
  real(dp),intent(in) :: qmass(:)
  real(dp),intent(in) :: strtarget(6)
  real(dp),intent(in) :: strten(6)
  real(dp),intent(inout) :: vel(3,natom)
 end subroutine isostress
end interface

interface
 subroutine matpointsym(iatom,mat3,natom,nsym,rprimd,symrel,tnons,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: iatom
  integer,intent(in) :: natom
  integer,intent(in) :: nsym
  real(dp),intent(inout) :: mat3(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine matpointsym
end interface

interface
 subroutine pimd_langevin_npt(etotal,forces,itimimage,natom,pimd_param,prtvolimg,&  
  &  rprimd,rprimd_next,rprimd_prev,stressin,trotter,vel,vel_cell,&  
  &  volume,xred,xred_next,xred_prev)
  use m_pimd
  use defs_basis
  implicit none
  integer,intent(in) :: itimimage
  integer,intent(in) :: natom
  integer,intent(in) :: prtvolimg
  integer,intent(in) :: trotter
  type(pimd_type),intent(in) :: pimd_param
  real(dp),intent(in) :: volume
  real(dp),intent(in) :: etotal(trotter)
  real(dp),intent(inout) :: forces(3,natom,trotter)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: rprimd_next(3,3)
  real(dp),intent(in) :: rprimd_prev(3,3)
  real(dp),intent(in) :: stressin(3,3,trotter)
  real(dp),intent(inout) :: vel(3,natom,trotter)
  real(dp),intent(inout) :: vel_cell(3,3,trotter)
  real(dp),intent(in),target :: xred(3,natom,trotter)
  real(dp),intent(out) :: xred_next(3,natom,trotter)
  real(dp),intent(in),target :: xred_prev(3,natom,trotter)
 end subroutine pimd_langevin_npt
end interface

interface
 subroutine pimd_langevin_nvt(etotal,forces,itimimage,natom,pimd_param,prtvolimg,&  
  &  rprimd,stressin,trotter,vel,volume,xred,xred_next,xred_prev)
  use m_pimd
  use defs_basis
  implicit none
  integer,intent(in) :: itimimage
  integer,intent(in) :: natom
  integer,intent(in) :: prtvolimg
  integer,intent(in) :: trotter
  type(pimd_type),intent(in) :: pimd_param
  real(dp),intent(in) :: volume
  real(dp),intent(in) :: etotal(trotter)
  real(dp),intent(inout) :: forces(3,natom,trotter)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: stressin(3,3,trotter)
  real(dp),intent(inout) :: vel(3,natom,trotter)
  real(dp),intent(in),target :: xred(3,natom,trotter)
  real(dp),intent(out) :: xred_next(3,natom,trotter)
  real(dp),intent(in),target :: xred_prev(3,natom,trotter)
 end subroutine pimd_langevin_nvt
end interface

interface
 subroutine pimd_nosehoover_npt(etotal,forces,itimimage,natom,pimd_param,prtvolimg,&  
  &  rprimd,rprimd_next,rprimd_prev,stressin,trotter,vel,vel_cell,&  
  &  volume,xred,xred_next,xred_prev)
  use m_pimd
  use defs_basis
  implicit none
  integer,intent(in) :: itimimage
  integer,intent(in) :: natom
  integer,intent(in) :: prtvolimg
  integer,intent(in) :: trotter
  type(pimd_type),intent(in) :: pimd_param
  real(dp),intent(in) :: volume
  real(dp),intent(in) :: etotal(trotter)
  real(dp),intent(inout) :: forces(3,natom,trotter)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: rprimd_next(3,3)
  real(dp),intent(in) :: rprimd_prev(3,3)
  real(dp),intent(in) :: stressin(3,3,trotter)
  real(dp),intent(inout) :: vel(3,natom,trotter)
  real(dp),intent(inout) :: vel_cell(3,3,trotter)
  real(dp),intent(in),target :: xred(3,natom,trotter)
  real(dp),intent(out) :: xred_next(3,natom,trotter)
  real(dp),intent(in),target :: xred_prev(3,natom,trotter)
 end subroutine pimd_nosehoover_npt
end interface

interface
 subroutine pimd_nosehoover_nvt(etotal,forces,itimimage,natom,pimd_param,prtvolimg,&  
  &  rprimd,stressin,trotter,vel,volume,xred,xred_next,xred_prev)
  use m_pimd
  use defs_basis
  implicit none
  integer,intent(in) :: itimimage
  integer,intent(in) :: natom
  integer,intent(in) :: prtvolimg
  integer,intent(in) :: trotter
  type(pimd_type),intent(inout) :: pimd_param
  real(dp),intent(in) :: volume
  real(dp),intent(in) :: etotal(trotter)
  real(dp),intent(inout) :: forces(3,natom,trotter)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: stressin(3,3,trotter)
  real(dp),intent(inout) :: vel(3,natom,trotter)
  real(dp),intent(in),target :: xred(3,natom,trotter)
  real(dp),intent(out) :: xred_next(3,natom,trotter)
  real(dp),intent(in),target :: xred_prev(3,natom,trotter)
 end subroutine pimd_nosehoover_nvt
end interface

interface
 subroutine prec_simple(ab_mover,forstr,hist,icycle,itime,iexit)
  use m_abimover
  use m_abihist
  implicit none
  integer,intent(in) :: icycle
  integer,intent(in) :: iexit
  integer,intent(in) :: itime
  type(abimover),intent(in) :: ab_mover
  type(abiforstr),intent(inout) :: forstr
  type(abihist),intent(inout) :: hist
 end subroutine prec_simple
end interface

interface
 subroutine pred_bfgs(ab_mover,ab_xfh,forstr,hist,ionmov,itime,zDEBUG,iexit)
  use m_abimover
  use m_abihist
  implicit none
  integer,intent(in) :: iexit
  integer,intent(in) :: ionmov
  integer,intent(in) :: itime
  type(abimover),intent(in) :: ab_mover
  type(ab_xfh_type),intent(inout) :: ab_xfh
  type(abiforstr),intent(in) :: forstr
  type(abihist),intent(inout) :: hist
  logical,intent(in) :: zDEBUG
 end subroutine pred_bfgs
end interface

interface
 subroutine pred_delocint(ab_mover,ab_xfh,forstr,hist,ionmov,itime,zDEBUG,iexit)
  use m_abimover
  use m_abihist
  implicit none
  integer,intent(in) :: iexit
  integer,intent(in) :: ionmov
  integer,intent(in) :: itime
  type(abimover),intent(inout) :: ab_mover
  type(ab_xfh_type),intent(inout) :: ab_xfh
  type(abiforstr),intent(in) :: forstr
  type(abihist),intent(inout) :: hist
  logical,intent(in) :: zDEBUG
 end subroutine pred_delocint
end interface

interface
 subroutine pred_diisrelax(ab_mover,hist,itime,ntime,zDEBUG,iexit)
  use m_abimover
  use m_abihist
  implicit none
  integer,intent(in) :: iexit
  integer,intent(in) :: itime
  integer,intent(in) :: ntime
  type(abimover),intent(in) :: ab_mover
  type(abihist),intent(inout) :: hist
  logical,intent(in) :: zDEBUG
 end subroutine pred_diisrelax
end interface

interface
 subroutine pred_isokinetic(ab_mover,hist,itime,ntime,zDEBUG,iexit)
  use m_abimover
  use m_abihist
  implicit none
  integer,intent(in) :: iexit
  integer,intent(in) :: itime
  integer,intent(in) :: ntime
  type(abimover),intent(in) :: ab_mover
  type(abihist),intent(inout) :: hist
  logical,intent(in) :: zDEBUG
 end subroutine pred_isokinetic
end interface

interface
 subroutine pred_isothermal(ab_mover,hist,itime,mttk_vars,ntime,zDEBUG,iexit)
  use m_abimover
  use m_abihist
  implicit none
  integer,intent(in) :: iexit
  integer,intent(in) :: itime
  integer,intent(in) :: ntime
  type(abimover),intent(in) :: ab_mover
  type(abihist),intent(inout) :: hist
  type(mttk_type),intent(inout) :: mttk_vars
  logical,intent(in) :: zDEBUG
 end subroutine pred_isothermal
end interface

interface
 subroutine pred_langevin(ab_mover,hist,icycle,itime,ncycle,ntime,zDEBUG,iexit,skipcycle)
  use m_abimover
  use m_abihist
  implicit none
  integer,intent(in) :: icycle
  integer,intent(in) :: iexit
  integer,intent(in) :: itime
  integer,intent(inout) :: ncycle
  integer,intent(in) :: ntime
  type(abimover),intent(in) :: ab_mover
  type(abihist),intent(inout) :: hist
  logical,intent(out) :: skipcycle
  logical,intent(in) :: zDEBUG
 end subroutine pred_langevin
end interface

interface
 subroutine pred_moldyn(ab_mover,hist,icycle,itime,ncycle,ntime,zDEBUG,iexit)
  use m_abimover
  use m_abihist
  implicit none
  integer,intent(in) :: icycle
  integer,intent(in) :: iexit
  integer,intent(in) :: itime
  integer,intent(inout) :: ncycle
  integer,intent(in) :: ntime
  type(abimover),intent(in) :: ab_mover
  type(abihist),intent(inout) :: hist
  logical,intent(in) :: zDEBUG
 end subroutine pred_moldyn
end interface

interface
 function fdtion(ab_mover,hist,itime)
  use defs_basis
  use m_abimover
  use m_abihist
  implicit none
  integer,intent(in) :: itime
  type(abimover),intent(in) :: ab_mover
  real(dp) :: fdtion
  type(abihist),intent(in) :: hist
 end function fdtion
end interface

interface
 subroutine pred_nose(ab_mover,hist,itime,ntime,zDEBUG,iexit)
  use m_abimover
  use m_abihist
  implicit none
  integer,intent(in) :: iexit
  integer,intent(in) :: itime
  integer,intent(in) :: ntime
  type(abimover),intent(in) :: ab_mover
  type(abihist),intent(inout) :: hist
  logical,intent(in) :: zDEBUG
 end subroutine pred_nose
end interface

interface
 subroutine pred_simple(ab_mover,hist,iexit)
  use m_abimover
  use m_abihist
  implicit none
  integer,intent(in) :: iexit
  type(abimover),intent(in) :: ab_mover
  type(abihist),intent(inout) :: hist
 end subroutine pred_simple
end interface

interface
 subroutine pred_srkna14(ab_mover,hist,icycle,zDEBUG,iexit,skipcycle)
  use m_abimover
  use m_abihist
  implicit none
  integer,intent(in) :: icycle
  integer,intent(in) :: iexit
  type(abimover),intent(in) :: ab_mover
  type(abihist),intent(inout) :: hist
  logical,intent(out) :: skipcycle
  logical,intent(in) :: zDEBUG
 end subroutine pred_srkna14
end interface

interface
 subroutine pred_steepdesc(ab_mover,forstr,hist,itime,zDEBUG,iexit)
  use m_abimover
  use m_abihist
  implicit none
  integer,intent(in) :: iexit
  integer,intent(in) :: itime
  type(abimover),intent(in) :: ab_mover
  type(abiforstr),intent(in) :: forstr
  type(abihist),intent(inout) :: hist
  logical,intent(in) :: zDEBUG
 end subroutine pred_steepdesc
end interface

interface
 subroutine pred_verlet(ab_mover,hist,ionmov,itime,ntime,zDEBUG,iexit)
  use m_abimover
  use m_abihist
  implicit none
  integer,intent(in) :: iexit
  integer,intent(in) :: ionmov
  integer,intent(in) :: itime
  integer,intent(in) :: ntime
  type(abimover),intent(in) :: ab_mover
  type(abihist),intent(inout) :: hist
  logical,intent(in) :: zDEBUG
 end subroutine pred_verlet
end interface

interface
 subroutine predict_copy(itimimage,list_dynimage,ndynimage,nimage,ntimimage,results_img)
  use m_results_img
  implicit none
  integer,intent(in) :: itimimage
  integer,intent(in) :: ndynimage
  integer,intent(in) :: nimage
  integer,intent(in) :: ntimimage
  integer,intent(in) :: list_dynimage(ndynimage)
  type(results_img_type),intent(inout) :: results_img(nimage,ntimimage)
 end subroutine predict_copy
end interface

interface
 subroutine predict_neb(itimimage,list_dynimage,mep_param,mpi_enreg,natom,&  
  &  ndynimage,nimage,nimage_tot,ntimimage,results_img)
  use m_results_img
  use m_mep
  use defs_abitypes
  implicit none
  integer,intent(in) :: itimimage
  integer,intent(in) :: natom
  integer,intent(in) :: ndynimage
  integer,intent(in) :: nimage
  integer,intent(in) :: nimage_tot
  integer,intent(in) :: ntimimage
  type(mep_type),intent(inout) :: mep_param
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: list_dynimage(ndynimage)
  type(results_img_type),intent(inout) :: results_img(nimage,ntimimage)
 end subroutine predict_neb
end interface

interface
 subroutine predict_pimd(imgmov,itimimage,mpi_enreg,natom,nimage,nimage_tot,&  
  &  ntimimage,pimd_param,prtvolimg,results_img)
  use m_pimd
  use m_results_img
  use defs_abitypes
  implicit none
  integer,intent(in) :: imgmov
  integer,intent(in) :: itimimage
  integer,intent(in) :: natom
  integer,intent(in) :: nimage
  integer,intent(in) :: nimage_tot
  integer,intent(in) :: ntimimage
  integer,intent(in) :: prtvolimg
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pimd_type),intent(inout) :: pimd_param
  type(results_img_type) :: results_img(nimage,ntimimage)
 end subroutine predict_pimd
end interface

interface
 subroutine predict_steepest(itimimage,list_dynimage,mep_param,natom,ndynimage,nimage,&  
  &  ntimimage,results_img)
  use m_results_img
  use m_mep
  implicit none
  integer,intent(in) :: itimimage
  integer,intent(in) :: natom
  integer,intent(in) :: ndynimage
  integer,intent(in) :: nimage
  integer,intent(in) :: ntimimage
  type(mep_type),intent(inout) :: mep_param
  integer,intent(in) :: list_dynimage(ndynimage)
  type(results_img_type) :: results_img(nimage,ntimimage)
 end subroutine predict_steepest
end interface

interface
 subroutine predict_string(itimimage,list_dynimage,mep_param,mpi_enreg,natom,&  
  &  ndynimage,nimage,nimage_tot,ntimimage,results_img)
  use m_results_img
  use m_mep
  use defs_abitypes
  implicit none
  integer,intent(in) :: itimimage
  integer,intent(in) :: natom
  integer,intent(in) :: ndynimage
  integer,intent(in) :: nimage
  integer,intent(in) :: nimage_tot
  integer,intent(in) :: ntimimage
  type(mep_type),intent(inout) :: mep_param
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: list_dynimage(ndynimage)
  type(results_img_type) :: results_img(nimage,ntimimage)
 end subroutine predict_string
end interface

interface
 subroutine predictimg(deltae,imagealgo_str,imgmov,itimimage,list_dynimage,&  
  &  ga_param,mep_param,mpi_enreg,natom,ndynimage,nimage,nimage_tot,ntimimage,&  
  &  pimd_param,prtvolimg,results_img)
  use m_results_img
  use m_ga
  use m_mep
  use m_pimd
  use defs_abitypes
  use defs_basis
  implicit none
  integer,intent(in) :: imgmov
  integer,intent(in) :: itimimage
  integer,intent(in) :: natom
  integer,intent(in) :: ndynimage
  integer,intent(in) :: nimage
  integer,intent(in) :: nimage_tot
  integer,intent(in) :: ntimimage
  integer,intent(in) :: prtvolimg
  real(dp),intent(in) :: deltae
  type(ga_type),intent(inout) :: ga_param
  character(len=60),intent(in) :: imagealgo_str
  type(mep_type),intent(inout) :: mep_param
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pimd_type),intent(inout) :: pimd_param
  integer,intent(in) :: list_dynimage(ndynimage)
  type(results_img_type) :: results_img(nimage,ntimimage)
 end subroutine predictimg
end interface

interface
 subroutine prtxfase(ab_mover,hist,iout,pos)
  use m_abimover
  use m_abihist
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: pos
  type(abimover),intent(in) :: ab_mover
  type(abihist),intent(in) :: hist
 end subroutine prtxfase
end interface

interface
 subroutine strainsym(nsym,rprimd0,rprimd,rprimd_symm,symrel)
  use defs_basis
  implicit none
  integer,intent(in) :: nsym
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: rprimd0(3,3)
  real(dp),intent(out) :: rprimd_symm(3,3)
  integer,intent(in) :: symrel(3,3,nsym)
 end subroutine strainsym
end interface

interface
 subroutine xcart2deloc(deloc,natom,rprimd,xcart,&  
  &  bt_inv_matrix,u_matrix,deloc_int,prim_int)
  use defs_basis
  use m_abimover
  implicit none
  integer,intent(in) :: natom
  type(delocint),intent(inout) :: deloc
  real(dp),intent(out) :: bt_inv_matrix(3*(natom-1),3*natom)
  real(dp),intent(out) :: deloc_int(3*(natom-1))
  real(dp),intent(out) :: prim_int(deloc%ninternal)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(inout) :: u_matrix(deloc%ninternal,3*(natom-1))
  real(dp),intent(in) :: xcart(3,natom)
 end subroutine xcart2deloc
end interface

interface
 subroutine calc_btinv_matrix(b_matrix,natom,ninternal,bt_inv_matrix,u_matrix)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: ninternal
  real(dp),intent(in) :: b_matrix(ninternal,3*natom)
  real(dp),intent(out) :: bt_inv_matrix(3*(natom-1),3*natom)
  real(dp),intent(inout) :: u_matrix(ninternal,3*(natom-1))
 end subroutine calc_btinv_matrix
end interface

interface
 subroutine align_u_matrices(natom,ninternal,u_matrix,u_matrix_old,s_matrix,f_eigs)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: ninternal
  real(dp),intent(inout) :: f_eigs(3*natom)
  real(dp),intent(inout) :: s_matrix(3*natom,3*natom)
  real(dp),intent(inout) :: u_matrix(ninternal,3*(natom-1))
  real(dp),intent(in) :: u_matrix_old(ninternal,3*(natom-1))
 end subroutine align_u_matrices
end interface

interface
 subroutine xfh_recover_deloc(ab_xfh,ab_mover,acell,acell0,cycl_main,&  
  &  fred,hessin,ndim,rprim,rprimd0,strten,ucvol,ucvol0,vin,vin_prev,&  
  &  vout,vout_prev,xred,deloc,deloc_int,deloc_force,btinv,gprimd,prim_int,&  
  &  u_matrix)
  use defs_basis
  use m_abimover
  implicit none
  integer,intent(out) :: cycl_main
  integer,intent(in) :: ndim
  type(abimover),intent(in) :: ab_mover
  type(ab_xfh_type),intent(inout) :: ab_xfh
  type(delocint),intent(inout) :: deloc
  real(dp),intent(inout) :: ucvol
  real(dp),intent(inout) :: ucvol0
  real(dp),intent(inout) :: acell(3)
  real(dp),intent(in) :: acell0(3)
  real(dp),intent(inout) :: btinv(3*(ab_mover%natom-1),3*ab_mover%natom)
  real(dp),intent(inout) :: deloc_force(3*(ab_mover%natom-1))
  real(dp),intent(inout) :: deloc_int(3*(ab_mover%natom-1))
  real(dp),intent(inout) :: fred(3,ab_mover%natom)
  real(dp),intent(inout) :: gprimd(3,3)
  real(dp),intent(inout) :: hessin(:,:)
  real(dp),intent(inout) :: prim_int(:)
  real(dp),intent(inout) :: rprim(3,3)
  real(dp),intent(inout) :: rprimd0(3,3)
  real(dp),intent(inout) :: strten(6)
  real(dp),intent(inout) :: u_matrix(:,:)
  real(dp),intent(inout) :: vin(:)
  real(dp),intent(inout) :: vin_prev(:)
  real(dp),intent(inout) :: vout(:)
  real(dp),intent(inout) :: vout_prev(:)
  real(dp),intent(inout) :: xred(3,ab_mover%natom)
 end subroutine xfh_recover_deloc
end interface

interface
 subroutine xfh_recover_new(ab_xfh,ab_mover,acell,acell0,cycl_main,fred,&  
  &  hessin,ndim,rprim,rprimd0,strten,ucvol,ucvol0,vin,vin_prev,vout,&  
  &  vout_prev,xred)
  use defs_basis
  use m_abimover
  implicit none
  integer,intent(out) :: cycl_main
  integer,intent(in) :: ndim
  type(abimover),intent(in) :: ab_mover
  type(ab_xfh_type),intent(inout) :: ab_xfh
  real(dp),intent(inout) :: ucvol
  real(dp),intent(inout) :: ucvol0
  real(dp),intent(inout) :: acell(3)
  real(dp),intent(in) :: acell0(3)
  real(dp),intent(inout) :: fred(3,ab_mover%natom)
  real(dp),intent(inout) :: hessin(:,:)
  real(dp),intent(inout) :: rprim(3,3)
  real(dp),intent(inout) :: rprimd0(3,3)
  real(dp),intent(inout) :: strten(6)
  real(dp),intent(inout) :: vin(:)
  real(dp),intent(inout) :: vin_prev(:)
  real(dp),intent(inout) :: vout(:)
  real(dp),intent(inout) :: vout_prev(:)
  real(dp),intent(inout) :: xred(3,ab_mover%natom)
 end subroutine xfh_recover_new
end interface

interface
 subroutine xfh_update(ab_xfh,acell,fred_corrected,natom,rprim,strten,xred)
  use defs_basis
  use m_abimover
  implicit none
  integer,intent(in) :: natom
  type(ab_xfh_type),intent(inout) :: ab_xfh
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: fred_corrected(3,natom)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(in) :: strten(6)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine xfh_update
end interface

interface
 subroutine xfpack_f2vout(fred,natom,ndim,optcell,strtarget,strten,ucvol,vout)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: ndim
  integer,intent(in) :: optcell
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: fred(3,natom)
  real(dp),intent(in) :: strtarget(6)
  real(dp),intent(in) :: strten(6)
  real(dp),intent(out) :: vout(ndim)
 end subroutine xfpack_f2vout
end interface

interface
 subroutine xfpack_vin2x(acell,acell0,natom,ndim,nsym,optcell,&  
  &  rprim,rprimd0,symrel,ucvol,ucvol0,vin,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: ndim
  integer,intent(in) :: nsym
  integer,intent(in) :: optcell
  real(dp),intent(out) :: ucvol
  real(dp),intent(in) :: ucvol0
  real(dp),intent(inout) :: acell(3)
  real(dp),intent(in) :: acell0(3)
  real(dp),intent(inout) :: rprim(3,3)
  real(dp),intent(in) :: rprimd0(3,3)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: vin(ndim)
  real(dp),intent(out) :: xred(3,natom)
 end subroutine xfpack_vin2x
end interface

interface
 subroutine xfpack_x2vin(acell,acell0,natom,ndim,nsym,optcell,&  
  &  rprim,rprimd0,symrel,ucvol,ucvol0,vin,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: ndim
  integer,intent(in) :: nsym
  integer,intent(in) :: optcell
  real(dp),intent(inout) :: ucvol
  real(dp),intent(in) :: ucvol0
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: acell0(3)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(in) :: rprimd0(3,3)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(out) :: vin(ndim)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine xfpack_x2vin
end interface

end module interfaces_45_geomoptim
!!***
