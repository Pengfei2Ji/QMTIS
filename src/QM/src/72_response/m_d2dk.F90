!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_d2dk
!! NAME
!! m_d2dk
!!
!! FUNCTION
!! This module defines structures and provides procedures used to compute the 2nd order Sternheimer 
!! equation with respect to the wave vector perturbation.
!!
!! COPYRIGHT
!!  Copyright (C) 2015-2016 ABINIT group (LB,MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_d2dk

 use defs_basis
 use m_profiling_abi

 implicit none

 private
!!***

!----------------------------------------------------------------------

!!****t* m_d2dk/d2dk_t
!! NAME
!!  d2dk_t
!!
!! FUNCTION
!!  Datatype gathering all needed data
!!
!! SOURCE

 type,public :: d2dk_t

!scalars
  integer :: ndir ! number of directions to consider
  integer :: nband_k ! number of bands
  integer :: size_wf ! number of coeffs in a wavefunction

!arrays
  integer :: idirs(2) ! directions of the perturbations (ndir=1 : idirs(1)=idirs(2) , ndir=2 : idirs(1)/=idirs(2))

  real(dp),allocatable :: RHS_Stern(:,:)
   ! Right-hand side of the 2nd order Sternheimer equation, for every bands.
   ! Namely, for a band "n" :
   ! |(RHS_Stern)_n> = (H^(2)-epsilon_n S^(2)) |u^(0)_n> + 2(H^(1)-epsilon_n S^(1))|u^(1)_n>
   !                       - 2 sum_m ( lambda_mn^(1) ( S^(1) |u^(0)_m> + S^(0) |u^(1)_m> )
   ! ( /!\ : in the sum, the index m runs over occupied bands )
   ! where :
   !  - epsilon_n = eigenvalue of the GS Hamiltonian
   !  - lambda_mn^(1) = <u^(0)_m| (H^(1)-(epsilon_n+epsilon_m)/2 S^(1) |u^(0)_n> (1st order Lagrange multiplier)
   !  - X^(2) = d^2X/(dk_dir1 dk_dir2)
   !  - 2X^(1)Y^(1) = dX/dk_dir1 dY/dk_dir2 + dX/dk_dir2 dY/dk_dir1
   ! Note : P_c^* will be apply to |(RHS_Stern)_n> in dfpt_cgwf.F90
   ! **
   ! Computed in "d2dk_init"

  real(dp),allocatable :: amn(:,:)
   ! Scalar needed for the "orhtonormalization condition", see "dcwavef(:,:)"
   ! Namely :
   ! A_mn = <u^(0)_m| S^(2) |u^(0)_n> + 2 <u^(1)_m| S^(0) |u^(1)_n>
   !    + 2 <u^(1)_m| S^(1) |u^(0)_n> + 2 <u^(0)_m| S^(1) |u^(1)_n>
   ! **
   ! Computed in "d2dk_init", stored only for testing purposes

  real(dp),allocatable :: dcwavef(:,:)
   ! Vector needed to enforce the "orhtonormalization condition" on 2nd order wave functions.
   ! Namely :
   ! |dcwavef_n> = -1/2 sum_m A_mn |u^(0)_m>
   ! **
   ! Computed in "d2dk_init"

  real(dp),allocatable :: lambda_mn(:,:)
   ! 2nd order Lagrange multiplier.
   ! Namely :
   ! lambda_mn =  <u^(0)_m| H^(2) |u^(0)_n> + 2 <u^(1)_m| H^(0) |u^(1)_n>
   !          + 2 <u^(1)_m| H^(1) |u^(0)_n> + 2 <u^(0)_m| H^(1) |u^(1)_n>
   !          - A_mn * (epsilon_m + epsilon_n) / 2 
   ! **
   ! Computed in "d2dk_init"

 end type d2dk_t

 public :: d2dk_accumulate_bands
 public :: d2dk_apply_hamiltonian
 public :: d2dk_destroy

!!***

!----------------------------------------------------------------------

CONTAINS  !==============================================================================
!!***

!----------------------------------------------------------------------

!!****f* m_d2dk/d2dk_accumulate_bands
!! NAME
!! d2dk_init
!!
!! FUNCTION
!! Compute the scalar product < vi | v1j > and add it to d2dk%lambda_mn
!! If necessary, compute < vi | v2j > and add it to d2dk%amn.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine d2dk_accumulate_bands(d2dk,choice,gs_hamkq,mpi_enreg,iband,jband,print_info,vi,v1j,v2j)

 use defs_basis
 use defs_abitypes
 use m_hamiltonian
 use m_cgtools

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'd2dk_accumulate_bands'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: choice,iband,jband,print_info
 type(d2dk_t),intent(inout) :: d2dk
 type(gs_hamiltonian_type),intent(in) :: gs_hamkq
 type(MPI_type),intent(in) :: mpi_enreg

!arrays
 real(dp),intent(in) :: vi(2,d2dk%size_wf),v1j(2,d2dk%size_wf),v2j(2,d2dk%size_wf)
 
!Local variables ---------------------------------------
!scalars
 integer :: nband_k,size_wf
 real(dp) :: dotr,doti,factor
 character(len=500) :: msg
 character(len=22) :: bra_i,ket_j,op1,op2

! *************************************************************************

 nband_k = d2dk%nband_k
 size_wf = d2dk%size_wf
 factor = one
 if(d2dk%ndir==1 .and. choice /= 3) factor = two ! in order to not compute same terms twice

 call dotprod_g(dotr,doti,gs_hamkq%istwf_k,size_wf,2,vi,v1j,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)

 if(print_info/=0) then
   select case (choice)
    case(1)
      write(bra_i,'(a,i1,a)') '    < du/dk_dir2(',iband,') | '
      write(ket_j,'(a,i1,a)') ' | du/dk_dir1(',jband,') >    '
      write(op1,'(a)') '        H^(0)         '
      write(op2,'(a)') '        S^(0)         '
    case(2)
      write(bra_i,'(a,i1,a)') '    < du/dk_dir2(',iband,') | '
      write(ket_j,'(a,i1,a)') ' |   u^(0)   (',jband,') >    '
      write(op1,'(a)') '       dH/dk_dir1     '
      write(op2,'(a)') '       dS/dk_dir1     '
    case(3)
      write(bra_i,'(a,i1,a)') '    <   u^(0)   (',iband,') | '
      write(ket_j,'(a,i1,a)') ' |   u^(0)   (',jband,') >    '
      write(op1,'(a)') 'd^2H/(dk_dir1 dk_dir2)'
      write(op2,'(a)') 'd^2S/(dk_dir1 dk_dir2)'
    case(4)
      write(bra_i,'(a,i1,a)') '    <   u^(0)   (',iband,') | '
      write(ket_j,'(a,i1,a)') ' | du/dk_dir1(',jband,') >    '
      write(op1,'(a)') '       dH/dk_dir2     '
      write(op2,'(a)') '       dS/dk_dir2     '
   end select
   write(msg,'(3a,2(a,es22.13E3))') bra_i,op1,ket_j,' = ',dotr,',',doti
   call wrtout(std_out,msg)
 end if

 d2dk%lambda_mn(:,iband+(jband-1)*nband_k) = factor*(/dotr,doti/) + d2dk%lambda_mn(:,iband+(jband-1)*nband_k)

 if (choice == 1 .or. gs_hamkq%usepaw==1) then
   call dotprod_g(dotr,doti,gs_hamkq%istwf_k,size_wf,2,vi,v2j,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)

   if(print_info/=0) then
     write(msg,'(3a,2(a,es22.13E3))') bra_i,op2,ket_j,' = ',dotr,',',doti
     call wrtout(std_out,msg)
   end if

   d2dk%amn(:,iband+(jband-1)*nband_k) = factor*(/dotr,doti/) + d2dk%amn(:,iband+(jband-1)*nband_k)

 end if ! end choice

end subroutine d2dk_accumulate_bands
!!***

!----------------------------------------------------------------------

!!****f* m_d2dk/d2dk_apply_hamiltonian
!! NAME
!! d2dk_apply_hamiltonian
!!
!! FUNCTION
!! Apply the KS Hamiltonian (or derivative) to an input wave function.
!! If asked, it also does some checks.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!      d2dk_init
!!
!! CHILDREN
!!
!! SOURCE

subroutine d2dk_apply_hamiltonian(cg_jband,d2dk,eig0_jband,eig1_k_jband,gs_hamkq,idir,ipert,&
                                  mpi_enreg,print_info,prtvol,rf_hamk_idir,&
                                  work1,work2,work3)

 use defs_basis
 use defs_abitypes
 use m_hamiltonian
 use m_cgtools

 use m_pawcprj,  only : pawcprj_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'd2dk_apply_hamiltonian'
 use interfaces_14_hidewrite
 use interfaces_66_wfs
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: print_info,prtvol,idir,ipert
 real(dp),intent(in) :: eig0_jband
 type(d2dk_t),intent(in) :: d2dk
 type(gs_hamiltonian_type),intent(inout) :: gs_hamkq
 type(rf_hamiltonian_type),intent(inout),target :: rf_hamk_idir
 type(MPI_type),intent(inout) :: mpi_enreg

!arrays
 real(dp),intent(inout) :: cg_jband(2,d2dk%size_wf*print_info)
 real(dp),intent(in) :: eig1_k_jband(2)
 real(dp),intent(inout) :: work1(2,d2dk%size_wf),work2(2,d2dk%size_wf),work3(2,d2dk%size_wf)

!Local variables ---------------------------------------
!scalars
 integer,parameter :: tim_getghc=1,tim_getgh1c=1,tim_getgh2c=1
! GETGH1C/GETGH2C OPTIONS :
 integer,parameter :: berryopt=0 ! no berry phase
 integer,parameter :: optlocal=0 ! no local contribution (k pertubation)
 integer,parameter :: optnl=2 ! non local terms must be computed here
 integer,parameter :: usevnl=0 ! non-local potential is not computed yet, to be computed here
 integer,parameter :: opt_gvnl1=0 ! gvnl1 is used as ouput only
! ****************
 integer :: cpopt,sij_opt,natom,nband_k,size_wf
 real(dp) :: dotr,doti,dotr2,doti2,tol_test
 character(len=500) :: msg
 
!arrays
 type(pawcprj_type),allocatable :: cwaveprj(:,:)
 real(dp) :: dum1(1,1)
 real(dp),allocatable :: gvnl1(:,:),gvnlc(:,:)
 
! *********************************************************************

 natom = gs_hamkq%natom
 nband_k = d2dk%nband_k
 size_wf = d2dk%size_wf
 sij_opt=1;if (gs_hamkq%usepaw==0) sij_opt=0
 tol_test=tol12; if (gs_hamkq%usepaw==1) tol_test=tol8

 if (ipert == 0) then
   cpopt=-1 ! nonlop option : <p_lmn|in> (and derivatives) are computed here (and not saved)
   ABI_ALLOCATE(gvnlc,(2,size_wf))

!  Test if < u^(0) | H^(0) | u^(0) > = eps_0
   if(print_info/=0) then

     call getghc(cpopt,cg_jband,cwaveprj,work2,work3,gs_hamkq,gvnlc,zero,mpi_enreg,1,prtvol,&
                 sij_opt,tim_getghc,0,select_k=KPRIME_H_KPRIME)

     call dotprod_g(dotr,doti,gs_hamkq%istwf_k,size_wf,2,cg_jband,work2,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
!     write(msg,'(2(a,es22.13E3))') 'D2DK TEST GETGHC (before) dotr = ',dotr,  '    doti = ',doti
!     call wrtout(std_out,msg)
!     write(msg,'(a,es22.13E3)') 'D2DK TEST GETGHC Eig = ',eig0_jband
!     call wrtout(std_out,msg)
     dotr = dotr - eig0_jband
     dotr = sqrt(dotr**2+doti**2)
     if (dotr > tol_test) then
       write(msg,'(a,es22.13E3)') 'D2DK TEST GETGHC : NOT PASSED dotr = ',dotr
       call wrtout(std_out,msg)
     end if

   end if

   call getghc(cpopt,work1,cwaveprj,work2,work3,gs_hamkq,gvnlc,zero,&
               mpi_enreg,1,prtvol,sij_opt,tim_getghc,0,select_k=KPRIME_H_KPRIME)
   ABI_DEALLOCATE(gvnlc)

 else if (ipert == natom+1) then

!  Test if < u^(0) | ( H^(1) - eps^(0) S^(1) ) | u^(0) > = eps^(1)
   if(print_info/=0) then

     call getgh1c(berryopt,cg_jband,cwaveprj,work2,dum1,work3,gs_hamkq,gvnl1,idir,ipert,zero,&
                mpi_enreg,optlocal,optnl,opt_gvnl1,rf_hamk_idir,sij_opt,tim_getgh1c,usevnl)

     call dotprod_g(dotr,doti,gs_hamkq%istwf_k,size_wf,2,cg_jband,work2,mpi_enreg%me_g0, mpi_enreg%comm_spinorfft)
     if (gs_hamkq%usepaw==1) then
       call dotprod_g(dotr2,doti2,gs_hamkq%istwf_k,size_wf,2,cg_jband,work3,mpi_enreg%me_g0, mpi_enreg%comm_spinorfft)
       dotr = dotr - eig0_jband*dotr2
       doti = doti - eig0_jband*doti2
     end if
!     write(msg,'(2(a,es22.13E3))') 'D2DK TEST GETGH1 (before) : dotr = ',dotr,' doti = ',doti
!     call wrtout(std_out,msg)
     dotr = dotr - eig1_k_jband(1)
     doti = doti - eig1_k_jband(2)
!     write(msg,'(2(a,es22.13E3))') 'D2DK TEST GETGH1 : Eig1 = ',eig1_k_jband(1),',',eig1_k_jband(2)
!     call wrtout(std_out,msg)
     dotr = sqrt(dotr**2+doti**2)
     if (dotr > tol_test) then
       write(msg,'(a,es22.13E3)') 'D2DK TEST GETGH1 (1) : NOT PASSED dotr = ',dotr
       call wrtout(std_out,msg)
     end if
   end if

   call getgh1c(berryopt,work1,cwaveprj,work2,dum1,work3,gs_hamkq,gvnl1,idir,ipert,zero,&
                mpi_enreg,optlocal,optnl,opt_gvnl1,rf_hamk_idir,sij_opt,tim_getgh1c,usevnl)

 else if (ipert == natom+10 .or. ipert == natom+11) then
!    Compute  : d^2H/(dk_dir1 dk_dir2)|u^(0)>  (in work2)
!    and      : d^2S/(dk_dir1 dk_dir2)|u^(0)>  (in work3)
     call getgh2c(work1,cwaveprj,work2,work3,gs_hamkq,dum1,idir,ipert,eig0_jband,&
                   mpi_enreg,optlocal,optnl,rf_hamk_idir,sij_opt,tim_getgh2c,usevnl)
 
 end if

end subroutine d2dk_apply_hamiltonian
!!***

!----------------------------------------------------------------------

!!****f* m_d2dk/d2dk_destroy
!! NAME
!! d2dk_init
!!
!! FUNCTION
!!  Free all allocated arrays in a d2dk_t object
!!
!! INPUTS
!! d2dk : the d2dk_t object to free
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine d2dk_destroy(d2dk)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'd2dk_destroy'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 type(d2dk_t),intent(inout) :: d2dk
!arrays

!Local variables ---------------------------------------
!scalars

! *************************************************************************

 if (allocated(d2dk%RHS_Stern)) then
   ABI_DEALLOCATE(d2dk%RHS_Stern)
 end if
 if (allocated(d2dk%dcwavef)) then
   ABI_DEALLOCATE(d2dk%dcwavef)
 end if
 if (allocated(d2dk%amn)) then
   ABI_DEALLOCATE(d2dk%amn)
 end if
 if (allocated(d2dk%lambda_mn)) then
   ABI_DEALLOCATE(d2dk%lambda_mn)
 end if

end subroutine d2dk_destroy
!!***

END MODULE m_d2dk
