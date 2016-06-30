!{\src2tex{textfont=tt}}
!!****f* ABINIT/d2dk_init
!!
!! NAME
!! d2dk_init
!!
!! FUNCTION
!! Compute terms needed for the 2nd order Sternheimer equation with respect k perturbation.
!! All terms are stored in a d2dk_t object.
!!
!! COPYRIGHT
!! Copyright (C) 2015-2016 ABINIT group (LB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cg(2,mpw*nspinor*mband*nsppol)=planewave coefficients of wavefunctions at k
!!  d2dk : the object we want to initialize
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  eig0_k(mband*nsppol)=GS eigenvalues at k (hartree)
!!  eig1_k(2*mband*mband*nsppol)=1st-order eigenvalues at k,q (hartree)
!!  gs_hamkq <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k+q
!!  icg=shift to be applied on the location of data in the array cg
!!  idir=direction of the perturbation
!!  ikpt=number of the k-point
!!  ipert=type of the perturbation
!!  isppol=index of current spin component
!!  mpi_enreg=informations about MPI parallelization
!!  mpw=maximum dimensioned size of npw or wfs at k
!!  nband_k=number of bands at this k point for that spin polarization
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  rf_hamkq <type(rf_hamiltonian_type)>=all data for the 1st-order Hamiltonian at k,q
!!  rf_hamk_dir2 <type(rf_hamiltonian_type)>= (used only when ipert=natom+11, so q=0)
!!    same as rf_hamkq, but the direction of the perturbation is different
!!  occ_k(nband_k)=occupation number for each band (usually 2) for each k.
!!  rocceig(nband_k,nband_k)= (occ_kq(m)-occ_k(n))/(eig0_kq(m)-eig0_k(n))
!!  wffddk=struct info for wf ddk file.
!!  ddk<wfk_t>=struct info for DDK file.
!!
!! OUTPUT
!!  d2dk%RHS_Stern
!!
!! NOTES
!!
!! PARENTS
!!      dfpt_vtowfk
!!
!! CHILDREN
!!      cg_zaxpy,cg_zscal,d2dk_accumulate_bands,d2dk_apply_hamiltonian
!!      dotprod_g,initmpi_seq,wffreaddatarec,wfk_read_bks,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine d2dk_init(cg,d2dk,dtset,eig0_k,eig1_k,gs_hamkq,icg,idir,ikpt,ipert,isppol,mpi_enreg,mpw,&
                     nband_k,nsppol,rf_hamkq,rf_hamk_dir2,occ_k,rocceig,wffddk,ddk_f)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi
 use m_errors
 use m_wffile
 use m_wfk
 use m_hamiltonian
 use m_cgtools
 use m_d2dk

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'd2dk_init'
 use interfaces_14_hidewrite
 use interfaces_51_manage_mpi
!End of the abilint section

 implicit none

! *************************************************************************
!Arguments -------------------------------
!scalars

 integer,intent(in) :: icg,idir,ipert,isppol,ikpt
 integer,intent(in) :: mpw,nband_k,nsppol
 type(dataset_type),intent(in) :: dtset
 type(gs_hamiltonian_type),intent(inout) :: gs_hamkq
 type(rf_hamiltonian_type),intent(inout),target :: rf_hamkq,rf_hamk_dir2
 type(MPI_type),intent(inout) :: mpi_enreg

!arrays
 real(dp),intent(in) :: cg(2,mpw*gs_hamkq%nspinor*dtset%mband*nsppol)
 real(dp),intent(in) :: eig0_k(dtset%mband)
 real(dp),intent(inout) :: eig1_k(2*dtset%mband**2) ! Here eig1_k contains 2nd order eigenvalues...
 real(dp),intent(in) :: occ_k(nband_k),rocceig(nband_k,nband_k)
 type(d2dk_t),intent(inout) :: d2dk
 type(wffile_type),intent(inout) :: wffddk(2)
 type(wfk_t),intent(inout) :: ddk_f(2)

!Local variables-------------------------------
!scalars
 integer,parameter :: berryopt=0,tim_getghc=1,tim_getgh1c=1,tim_getgh2c=1,level=19
 integer :: factor,iband,idir1,idir2,ierr
 integer :: igs,ipert1,ipw,ispinor,jband,kdir1
 integer :: master,me,my_comm_atom,natom,print_info
 integer :: size_wf,spaceworld,shift_band1,shift_band2
 integer :: shift_dir1_lambda,shift_dir1,shift_dir2,shift_jband_lambda
 real(dp) :: doti,dotr,invocc,tol_final
 logical :: paral_atom
 character(len=500) :: msg
 type(MPI_type) :: mpi_enreg_seq
!arrays
 integer,pointer :: my_atmtab(:)
 real(dp) :: lambda_ij(2),eig1_k_jband(2)
 real(dp),allocatable :: cg_jband(:,:),dsusdu(:,:),dudk(:,:)
 real(dp),allocatable :: eig1_k_stored(:),eig1_k_tmp(:)
 real(dp),allocatable :: work1(:,:),work2(:,:),work3(:,:)
 type(rf_hamiltonian_type),pointer :: rf_hamk_idir

! *********************************************************************

 DBG_ENTER("COLL")

#ifdef DEV_MG_WFK
 ABI_UNUSED((/wffddk(1)%unwff/))
#endif

 size_wf=gs_hamkq%npw_k*gs_hamkq%nspinor
 natom = gs_hamkq%natom
 print_info = 0
 if (dtset%prtvol == -level) print_info = 1 ! also active a lot of tests
!Define some atrributes of the d2dk object
 d2dk%nband_k = nband_k
 d2dk%size_wf = size_wf

 if(ipert<natom+10.or.ipert>natom+11) then
   write(msg,'(a)') 'ipert must be equal to natom+10 or natom+11 for d2dk calculations.'
   MSG_BUG(msg)
 end if

!Set up parallelism
 master=0;me=mpi_enreg%me_kpt
 spaceworld=mpi_enreg%comm_cell
 paral_atom=(mpi_enreg%my_natom/=natom)
 my_comm_atom=mpi_enreg%comm_atom
! my_natom=mpi_enreg%my_natom
 my_atmtab=>mpi_enreg%my_atmtab

!Fake MPI data to be used in sequential calls to parallel routines
 call initmpi_seq(mpi_enreg_seq)
 mpi_enreg_seq%my_natom=natom

!Define the two directions of the perturbation
 if (ipert==natom+10) then ! diagonal term
   d2dk%ndir=1
   d2dk%idirs(1)=idir ; d2dk%idirs(2)=idir
 else ! case ipert=natom+11 ! off-diagonal term
   d2dk%ndir=2
   if(idir==1) then
     d2dk%idirs(1)=2 ; d2dk%idirs(2)=3
   else if (idir==2) then
     d2dk%idirs(1)=1 ; d2dk%idirs(2)=3
   else if (idir==3) then
     d2dk%idirs(1)=1 ; d2dk%idirs(2)=2
   end if
 end if

! Allocate work space for one band
 ABI_ALLOCATE(work1,(2,size_wf))

! **************************************************************************************************
! Get info from ddk files
! **************************************************************************************************

! "eig1_k_stored" contains dLambda_{nm}/dk_dir every bands n and m and ndir (=1 or 2) directions
 ABI_ALLOCATE(eig1_k_stored,(2*d2dk%ndir*nband_k**2))

! "dudk" contains du/dk_dir for every bands and ndir (=1 or 2) directions
 ABI_STAT_ALLOCATE(dudk,(2,d2dk%ndir*nband_k*size_wf), ierr)
 ABI_CHECK(ierr==0, "out of memory in m_d2dk : dudk")
 
 ABI_ALLOCATE(eig1_k_tmp,(2*nband_k))
 
 if(print_info/=0) then
   write(msg,'(4(a,i2))') 'D2DK_INIT : ipert-natom = ',ipert-natom,' , idir = ',idir,&
   ' , ikpt = ',ikpt,' , isppol = ',isppol
   call wrtout(std_out,msg,'COLL')
 end if

 do kdir1=1,d2dk%ndir
   idir1=d2dk%idirs(kdir1)
   do iband=1,nband_k
#ifndef DEV_MG_WFK
     call WffReadDataRec(eig1_k_tmp,ierr,2*nband_k,wffddk(kdir1))
     call WffReadDataRec(work1,ierr,2,size_wf,wffddk(kdir1))
#else
     call wfk_read_bks(ddk_f(kdir1), iband, ikpt, isppol, xmpio_single, cg_bks=work1,eig1_bks=eig1_k_tmp)
#endif
!    Filter the wavefunctions for large modified kinetic energy
!    The GS wavefunctions should already be non-zero
     do ispinor=1,gs_hamkq%nspinor
       igs=(ispinor-1)*gs_hamkq%npw_k
       do ipw=1+igs,gs_hamkq%npw_k+igs
         if(gs_hamkq%kinpw_kp(ipw-igs)>huge(zero)*1.d-11)then
           work1(1,ipw)=zero
           work1(2,ipw)=zero
         end if
       end do
     end do
!    Copy work1 in "dudk"
     shift_band1=(iband-1)*size_wf
     shift_dir1=(kdir1-1)*nband_k*size_wf
     dudk(:,1+shift_band1+shift_dir1:size_wf+shift_band1+shift_dir1)=work1(:,:)
!    Copy eig1_k_tmp in "eig1_k_stored"
     shift_band1=(iband-1)*2*nband_k
     shift_dir1=2*(kdir1-1)*nband_k**2
     eig1_k_stored(1+shift_band1+shift_dir1:2*nband_k+shift_band1+shift_dir1)=eig1_k_tmp(:)
   end do
 end do

 ABI_DEALLOCATE(eig1_k_tmp)

! **************************************************************************************************
 ipert1 = gs_hamkq%natom+1 ! For H^(1) and S^(1)
! **************************************************************************************************

! **************************************************************************************************
! COMPUTATION OF "dsusdu", A PART OF "A_mn" AND A PART OF "Lambda_mn" (see defs in m_d2dk)
! **************************************************************************************************

! "dsusdu" contains dS/dk_dir |u_band> + S|du_band/dk_dir> for every bands and ndir (=1 or 2) directions 
 ABI_STAT_ALLOCATE(dsusdu,(2,d2dk%ndir*nband_k*size_wf), ierr)
 ABI_CHECK(ierr==0, "out of memory in d2dk_init : dsusdu")
 dsusdu=zero

 ABI_ALLOCATE(d2dk%amn,(2,nband_k**2))
 d2dk%amn=zero

 ABI_ALLOCATE(d2dk%lambda_mn,(2,nband_k**2))
 d2dk%lambda_mn(:,:)=zero

! Allocate work spaces for one band
 ABI_ALLOCATE(work2,(2,size_wf))
 ABI_ALLOCATE(work3,(2,size_wf))
 ABI_ALLOCATE(cg_jband,(2,size_wf*print_info))

 factor=1
 if(ipert==natom+10) factor=2 ! in order to not compute same terms twice

 do kdir1=1,d2dk%ndir
   idir1=d2dk%idirs(kdir1)
   if(ipert==natom+10) then
     shift_dir1=0
     shift_dir2=0
   else
     shift_dir1=(kdir1-1)*nband_k*size_wf
     shift_dir2=(2-kdir1)*nband_k*size_wf
   end if
   if (kdir1==1) rf_hamk_idir => rf_hamkq
   if (kdir1==2) rf_hamk_idir => rf_hamk_dir2

!  LOOP OVER BANDS
   do jband=1,nband_k ! = band n
     shift_band1=(jband-1)*size_wf

     if (abs(occ_k(jband))>tol8) then

       if(print_info/=0) cg_jband(:,:) = cg(:,1+shift_band1+icg:size_wf+shift_band1+icg)

!      Extract first order wavefunction for jband (in work1)
       work1(:,:)=dudk(:,1+shift_band1+shift_dir1:size_wf+shift_band1+shift_dir1)

!      Compute H^(0) | du/dk_dir1 > (in work2) and S^(0) | du/dk_dir1 > (in work3)
       call d2dk_apply_hamiltonian(cg_jband,d2dk,eig0_k(jband),eig1_k_stored,gs_hamkq,idir1,0,&
       mpi_enreg,print_info,dtset%prtvol,rf_hamk_idir,work1,work2,work3)

       if (gs_hamkq%usepaw==0) work3(:,:)=work1(:,:) ! Store | du/dk_dir1 > in work3

!      Copy infos in dsusdu
       dsusdu(:,1+shift_band1+shift_dir1:size_wf+shift_band1+shift_dir1)=work3(:,:)&
       +dsusdu(:,1+shift_band1+shift_dir1:size_wf+shift_band1+shift_dir1)

!      For every occupied iband, we compute :
!      < du/dk_dir2(iband) | H^(0) | du/dk_dir1(jband) > and add it to lambda_mn
!      < du/dk_dir2(iband) | S^(0) | du/dk_dir1(jband) > and add it to amn
       do iband=1,d2dk%nband_k  ! = band m
         if (abs(occ_k(iband))>tol8) then
           shift_band2=(iband-1)*size_wf
           work1(:,:)=dudk(:,1+shift_band2+shift_dir2:size_wf+shift_band2+shift_dir2)
           call d2dk_accumulate_bands(d2dk,1,gs_hamkq,mpi_enreg,iband,jband,print_info,&
           work1,work2,work3)
         end if
       end do

       if (print_info/=0) then
         eig1_k_jband(1) = eig1_k_stored(jband  +(jband-1)*(2*nband_k+1)+2*(kdir1-1)*nband_k**2) ! dir1
         eig1_k_jband(2) = eig1_k_stored(jband+1+(jband-1)*(2*nband_k+1)+2*(kdir1-1)*nband_k**2) ! dir1
       end if

!      Extract GS wavefunction for jband
       work1 = cg(:,1+shift_band1+icg:size_wf+shift_band1+icg)

!      Compute dH/dk_dir1 | u^(0) > (in work2) and dS/dk_dir1 | u^(0) > (in work3)
       call d2dk_apply_hamiltonian(cg_jband,d2dk,eig0_k(jband),eig1_k_jband,gs_hamkq,idir1,ipert1,&
       mpi_enreg,print_info,dtset%prtvol,rf_hamk_idir,work1,work2,work3)

!      Copy infos in dsusdu
       if (gs_hamkq%usepaw==1) then
         dsusdu(:,1+shift_band1+shift_dir1:size_wf+shift_band1+shift_dir1)=work3(:,:)&
         +dsusdu(:,1+shift_band1+shift_dir1:size_wf+shift_band1+shift_dir1)
       end if

!      For every occupied iband, we compute :
!      < du/dk_dir2(iband) | dH/dk_dir1 | u^(0)(jband) > and add it to lambda_mn
!      < du/dk_dir2(iband) | dS/dk_dir1 | u^(0)(jband) > and add it to amn
       do iband=1,d2dk%nband_k  ! = band m
         if (abs(occ_k(iband))>tol8) then
           shift_band2=(iband-1)*size_wf
           work1(:,:)=dudk(:,1+shift_band2+shift_dir2:size_wf+shift_band2+shift_dir2)
           call d2dk_accumulate_bands(d2dk,2,gs_hamkq,mpi_enreg,iband,jband,print_info,&
           work1,work2,work3)
         end if
       end do

     end if ! empty band test
   end do ! jband
 end do ! idir1

! **************************************************************************************************
! COMPUTATION OF "RHS_Stern", THE LAST PART OF "A_mn" AND A PART OF "Lambda_mn"
! **************************************************************************************************

 ABI_STAT_ALLOCATE(d2dk%RHS_Stern,(2,nband_k*size_wf), ierr)
 ABI_CHECK(ierr==0, "out of memory in m_d2dk : RHS_Stern")
 d2dk%RHS_Stern(:,:)=zero

 do jband=1,nband_k
   if (abs(occ_k(jband))>tol8) then
     shift_band1=(jband-1)*size_wf

!    Extract GS wavefunction (in work1)
     work1(:,:)=cg(:,1+shift_band1+icg:size_wf+shift_band1+icg)

!    Compute  : d^2H/(dk_dir1 dk_dir2)|u^(0)>  (in work2)
!    and      : d^2S/(dk_dir1 dk_dir2)|u^(0)>  (in work3)
     call d2dk_apply_hamiltonian(cg_jband,d2dk,eig0_k(jband),eig1_k_jband,gs_hamkq,idir,ipert,&
     mpi_enreg,print_info,dtset%prtvol,rf_hamkq,work1,work2,work3)

!    For every occupied iband, we compute :
!    < u^(0)(iband) | d^2H/(dk_dir1 dk_dir2) | u^(0)(jband) > and add it to lambda_mn
!    < u^(0)(iband) | d^2S/(dk_dir1 dk_dir2) | u^(0)(jband) > and add it to amn
     do iband=1,d2dk%nband_k  ! = band m
       if (abs(occ_k(iband))>tol8) then
         shift_band2=(iband-1)*size_wf
         work1(:,:)=cg(:,1+shift_band2+icg:size_wf+shift_band2+icg)
         call d2dk_accumulate_bands(d2dk,3,gs_hamkq,mpi_enreg,iband,jband,print_info,&
         work1,work2,work3)
       end if
     end do

!    Add d^2H/(dk_dir1 dk_dir2)|u^(0)> to RHS_Stern :
     if (gs_hamkq%usepaw==1) work2(:,:)=work2(:,:)-eig0_k(jband)*work3(:,:) ! if PAW : we add H^(2)-eps^(0) S^(2)
     work1(:,:)=d2dk%RHS_Stern(:,1+shift_band1:size_wf+shift_band1)
     call cg_zaxpy(size_wf,(/one,zero/),work2,work1)
     d2dk%RHS_Stern(:,1+shift_band1:size_wf+shift_band1)=work1(:,:)

     do kdir1=1,d2dk%ndir
       shift_dir1_lambda=(kdir1-1)*2*nband_k**2
       if(ipert==natom+10) then
         idir2=d2dk%idirs(kdir1)
         shift_dir1=0
         shift_dir2=0
         rf_hamk_idir => rf_hamkq
       else
         idir2=d2dk%idirs(2-kdir1+1)
         shift_dir1=(kdir1-1)*nband_k*size_wf
         shift_dir2=(2-kdir1)*nband_k*size_wf
         if (kdir1==1) rf_hamk_idir => rf_hamk_dir2 ! dir2
         if (kdir1==2) rf_hamk_idir => rf_hamkq ! dir1
       end if

       if(print_info/=0) then
         cg_jband(:,:) = cg(:,1+shift_band1+icg:size_wf+shift_band1+icg)
         if(ipert == natom+10) then
           eig1_k_jband(1) = eig1_k_stored(jband  +(jband-1)*(2*nband_k+1)) ! dir2 = dir1
           eig1_k_jband(2) = eig1_k_stored(jband+1+(jband-1)*(2*nband_k+1)) ! dir2 = dir1
         else
           eig1_k_jband(1) = eig1_k_stored(jband  +(jband-1)*(2*nband_k+1)+2*(2-kdir1)*nband_k**2) ! dir2
           eig1_k_jband(2) = eig1_k_stored(jband+1+(jband-1)*(2*nband_k+1)+2*(2-kdir1)*nband_k**2) ! dir2
         end if
       end if

!      Extract first order wavefunction : | du/dk_dir1 > (in work1)
       work1(:,:)=dudk(:,1+shift_band1+shift_dir1:size_wf+shift_band1+shift_dir1)

!      Compute dH/dk_dir2 | du/dk_dir1 > (in work2) and dS/dk_dir2 | du/dk_dir1 > (in work3)
       call d2dk_apply_hamiltonian(cg_jband,d2dk,eig0_k(jband),eig1_k_jband,gs_hamkq,idir2,ipert1,&
       mpi_enreg,print_info,dtset%prtvol,rf_hamk_idir,work1,work2,work3)

!      For every occupied iband, we compute :
!      < u^(0)(iband) | dH/dk_dir2 | du/dk_dir1(jband) > and add it to lambda_mn
!      < u^(0)(iband) | dS/dk_dir2 | du/dk_dir1(jband) > and add it to amn
       do iband=1,d2dk%nband_k  ! = band m
         if (abs(occ_k(iband))>tol8) then
           shift_band2=(iband-1)*size_wf
           work1(:,:)=cg(:,1+shift_band2+icg:size_wf+shift_band2+icg)
           call d2dk_accumulate_bands(d2dk,4,gs_hamkq,mpi_enreg,iband,jband,print_info,&
           work1,work2,work3)
         end if
       end do

!      Add dH/dk_dir2 | du/dk_dir1 > to RHS_Stern :
       if (gs_hamkq%usepaw==1) work2(:,:)=work2(:,:)-eig0_k(jband)*work3(:,:) ! if PAW : we add H^(1)-eps^(0) S^(1)
       work1(:,:)=d2dk%RHS_Stern(:,1+shift_band1:size_wf+shift_band1)
       call cg_zaxpy(size_wf,(/factor*one,zero/),work2,work1)
       d2dk%RHS_Stern(:,1+shift_band1:size_wf+shift_band1)=work1(:,:)

!      Compute : -factor * sum_iband ( Lambda^(1)_{iband,jband} * dsusdu_{iband} )
       shift_jband_lambda=(jband-1)*2*nband_k
       do iband=1,nband_k
         if (abs(occ_k(iband))>tol8) then ! if empty band, nothing to do

!          Extract lambda_ij(iband,jband) for dir1
           lambda_ij(1)=eig1_k_stored(2*iband-1+shift_jband_lambda+shift_dir1_lambda)
           lambda_ij(2)=eig1_k_stored(2*iband  +shift_jband_lambda+shift_dir1_lambda)

!          Extract dsusdu for iband and dir2 (in work2)
           shift_band2=(iband-1)*size_wf
           work2(:,:)=dsusdu(:,1+shift_band2+shift_dir2:size_wf+shift_band2+shift_dir2)

!          Compute Lambda_{iband,jband} * dsusdu_{iband} (in work2)
           call cg_zscal(size_wf,lambda_ij,work2)

!          Add it to RHS_Stern
           work1(:,:)=d2dk%RHS_Stern(:,1+shift_band1:size_wf+shift_band1)
           call cg_zaxpy(size_wf,(/-factor*one,zero/),work2,work1) ! do not forget the minus sign!
           d2dk%RHS_Stern(:,1+shift_band1:size_wf+shift_band1)=work1(:,:)

         end if ! empty band test
       end do ! iband

     end do ! kdir1

   end if ! empty band test
 end do ! jband

 ABI_DEALLOCATE(cg_jband)
 ABI_DEALLOCATE(dudk)
 ABI_DEALLOCATE(dsusdu)
 ABI_DEALLOCATE(eig1_k_stored)

! Compute the part of 2nd order wavefunction that belongs to the space of empty bands
 do jband=1,nband_k
   shift_band1=(jband-1)*size_wf
   if (abs(occ_k(jband))>tol8) then
     invocc = one/occ_k(jband)
     work1(:,:)=d2dk%RHS_Stern(:,1+shift_band1:size_wf+shift_band1)
     do iband=1,nband_k
       if (iband /= jband) then
         if(print_info/=0) then
           if (abs(occ_k(iband) - occ_k(jband)) > tol12 .and. occ_k(iband) > tol8) then
             write(msg,'(a,i2,a,i2)') 'D2DK TEST ACTIVE SPACE : jband = ',jband,' iband = ',iband
             call wrtout(std_out,msg)
             write(msg,'(a)') 'ERROR : occ_k(iband) /= occ_k(jband) (and both are >0)'
             call wrtout(std_out,msg)
           end if
           if (abs(eig0_k(iband) - eig0_k(jband)) < tol8 ) then
             write(msg,'(a,i2,a,i2)') 'D2DK TEST ACTIVE SPACE : jband = ',jband,' iband = ',iband
             call wrtout(std_out,msg)
             write(msg,'(a,es22.13e3)') 'WARNING : DEGENERATE BANDS  Eig(jband) = Eig(jband) = ',eig0_k(jband)
             call wrtout(std_out,msg)
           end if
           if ( (eig0_k(iband) - eig0_k(jband) < -tol12) .and. (jband < iband) ) then
             write(msg,'(a,i2,a,i2)') 'D2DK TEST ACTIVE SPACE : jband = ',jband,' iband = ',iband
             call wrtout(std_out,msg)
             write(msg,'(a)') 'ERROR : Eig(jband) < Eig(iband) with jband < iband'
             call wrtout(std_out,msg)
             write(msg,'(a,es22.13e3)') 'Eig(jband) = ',eig0_k(jband)
             call wrtout(std_out,msg)
             write(msg,'(a,es22.13e3)') 'Eig(iband) = ',eig0_k(iband)
             call wrtout(std_out,msg)
           end if
         end if ! end tests
         if ( abs(occ_k(iband))<tol8 ) then ! for empty bands only
           shift_band2=(iband-1)*size_wf
           work2(:,:)=cg(:,1+shift_band2+icg:size_wf+shift_band2+icg)
           call dotprod_g(dotr,doti,gs_hamkq%istwf_k,size_wf,2,work2,work1,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
!          Store it in a_mn
!          /!\ There is a factor "-2" to simplify the use of amn in the following.
!          /!\ Occupied and empty bands will be treated in a same way.
           d2dk%amn(:,iband+(jband-1)*nband_k)=-two*rocceig(iband,jband)*invocc*(/dotr,doti/)&
           +d2dk%amn(:,iband+(jband-1)*nband_k)
         end if ! empty band test
       end if ! iband \= jband
     end do ! iband
   end if  ! empty band test
 end do ! jband

! **************************************************************************************************
!  COMPUTATION OF "dcwavef" AND "Lambda_mn" FROM "A_mn"
! **************************************************************************************************

 ABI_STAT_ALLOCATE(d2dk%dcwavef,(2,nband_k*size_wf), ierr)
 ABI_CHECK(ierr==0, "out of memory in m_d2dk : dcwavef")
 d2dk%dcwavef=zero

 do jband=1,nband_k
   shift_band1=(jband-1)*size_wf
   if (abs(occ_k(jband))>tol8) then
     do iband=1,nband_k
       shift_band2=(iband-1)*size_wf

!      Extract GS wavefunction for iband
       work1(:,:)=cg(:,1+shift_band2+icg:size_wf+shift_band2+icg)

       work2(:,:)=d2dk%dcwavef(:,1+shift_band1:size_wf+shift_band1)
       call cg_zaxpy(size_wf,-half*d2dk%amn(:,iband+(jband-1)*nband_k),work1,work2)
       d2dk%dcwavef(:,1+shift_band1:size_wf+shift_band1)=work2(:,:)

       if (abs(occ_k(iband))>tol8 .and. abs(occ_k(jband))>tol8) then
         d2dk%lambda_mn(:,iband+(jband-1)*nband_k) = -half*(eig0_k(iband)+eig0_k(jband))*d2dk%amn(:,iband+(jband-1)*nband_k)&
         + d2dk%lambda_mn(:,iband+(jband-1)*nband_k)

         eig1_k(2*iband-1+(jband-1)*2*nband_k) = d2dk%lambda_mn(1,iband+(jband-1)*nband_k)
         eig1_k(2*iband  +(jband-1)*2*nband_k) = d2dk%lambda_mn(2,iband+(jband-1)*nband_k)

       end if ! empty band test
     end do ! iband
   end if ! empty band test
 end do ! jband

! **************************************************************************************************
!  FINAL TEST
! **************************************************************************************************

 tol_final = tol6
 if(print_info/=0) then
   do jband=1,nband_k
     if (abs(occ_k(jband))>tol8) then
       shift_band1=(jband-1)*size_wf
       work1(:,:)=d2dk%RHS_Stern(:,1+shift_band1:size_wf+shift_band1)
       work2(:,:)=cg(:,1+shift_band1+icg:size_wf+shift_band1+icg)
       call dotprod_g(dotr,doti,gs_hamkq%istwf_k,size_wf,2,work1,work2,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
       dotr = dotr - d2dk%lambda_mn(1,jband+(jband-1)*nband_k)
       doti = doti - d2dk%lambda_mn(2,jband+(jband-1)*nband_k)
       dotr = sqrt(dotr**2+doti**2)
       if (dotr > tol_final) then
         write(msg,'(a,i2,a,es22.13E3)') 'D2DK TEST FINAL iband = ',jband,' : NOT PASSED dotr = ',dotr
         call wrtout(std_out,msg)
       else
         write(msg,'(a,i2,a,es22.13E3,a,es7.1E2)') &
         'D2DK TEST FINAL iband = ',jband,' : OK. |test| = ',dotr,' < ',tol_final
         call wrtout(std_out,msg)
       end if
     end if
   end do
 end if

! **************************************************************************************************
!  JOB FINISHED
! **************************************************************************************************

! Deallocations of arrays
 if (print_info==0) ABI_DEALLOCATE(d2dk%amn)
 ABI_DEALLOCATE(work1)
 ABI_DEALLOCATE(work2)
 ABI_DEALLOCATE(work3)

! call timab(566,2,tsec)

 DBG_EXIT("COLL")

end subroutine d2dk_init
!!***
