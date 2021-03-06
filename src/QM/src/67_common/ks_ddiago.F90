!{\src2tex{textfont=tt}}
!!****f* ABINIT/ks_ddiago
!! NAME
!! ks_ddiago
!!
!! FUNCTION
!!  This routine performs the direct diagonalization of the Kohn-Sham Hamiltonian
!!  for a given k-point and spin. The routine drives the following operations:
!!
!!    1) Re-computing <G|H|G_prim> matrix elements for all (G, G_prim).
!!       starting from the knowledge of the local potential on the real-space FFT mesh.
!!
!!    2) Diagonalizing H in the plane-wave basis.
!!
!!  It is called in outkss.F90 during the generation of the KSS file
!!  needed for a GW post-treatment. Since many-body calculations usually
!!  require a large number of eigenstates eigen-functions, a direct
!!  diagonalization of the Hamiltonian might reveal more stable than iterative
!!  techniques that might be problematic when several high energy states are required.
!!  The main drawback of the direct diagonalization is the bad scaling with the size
!!  of the basis set (npw**3) and the large memory requirements.
!!  At present, only norm-conserving pseudopotentials are implemented.
!!
!! COPYRIGHT
!! Copyright (C) 2000-2016 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  kpoint(3)
!!  prtvol=Integer Flags  defining verbosity level
!!  ecut=cut-off energy for plane wave basis sphere (Ha)
!!  mgfftc=maximum size of 1D FFTs (coarse mesh).
!!  natom=number of atoms in cell.
!!  nfftf=(effective) number of FFT grid points in the dense FFT mesh (for this processor)
!!         (nfftf=nfft for norm-conserving potential runs)
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nspden=number of density components
!!  Pawtab(Psps%ntypat*Psps%usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  Pawfgr<pawfgr_type>=fine grid parameters and related data
!!  Paw_ij(natom) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  Psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rprimd(3,3)=dimensional primitive translations for real space (bohr)
!!  vtrial(nfftf,nspden)=the trial potential
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  comm=MPI communicator.
!!  [Electronpositron] <electronpositron_type>=quantities for the electron-positron annihilation.
!!  nfftc=Number of points in the coarse FFT mesh.
!!  ngfftc(18)=Info about 3D FFT for the coarse mesh, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  Diago_ctl<ddiago_ctl_type>=Datatype storing variables and options controlling the direct diagonalization.
!!
!! OUTPUT
!!  ierr=Status error.
!!  onband_diago
!!
!! SIDE EFFECTS
!!  eig_ene(:)=Pointer used for allocating and storing the eigenvalues (hartree)
!!    input: pointer to NULL
!!     output: eig_ene(onband_diago)=The calculatated eigenvalues in ascending order.
!!
!!  eig_vec(:,:,:)=Pointer used for allocating and holding the wave functions at this k-point and spin.
!!    input: pointer to NULL
!!    output: eig_vec(2,npw_k*nspinor,onband_diago)=The calculated eigenvectors.
!!
!!  Cprj_k(natom,nspinor*onband_diago) PAW only===
!!   input: pointer to NULL
!!   output: Projected eigenstates <Proj_i|Cnk> from output eigenstates.
!!
!! NOTES
!! * The routine can be time consuming (in particular when computing <G1|H|G2> elements for all (G1,G2)).
!!   So, it is recommended to call it once per run.
!!
!! * The routine RE-compute all Hamiltonian terms. So it is equivalent to an additional electronic SC cycle.
!!   (This has no effect is convergence was reach. If not, eigenvalues/vectors may differs from the conjugate gradient ones)
!!
!! NOTES
!!  Please, do NOT pass Dtset% to this routine. Either use a local variable properly initialized
!!  or add the additional variable to ddiago_ctl_type and change the creation method accordingly.
!!  ks_ddiago is designed such that it is possible to diagonalize the Hamiltonian at an arbitrary k-point
!!  or spin (not efficient but easy to code). Therefore ks_ddiago is useful non only for
!!  the KSS generation but also for testing more advanced iterative algorithms as well as interpolation techniques.
!!
!! PARENTS
!!      m_shirley,outkss
!!
!! CHILDREN
!!      destroy_hamiltonian,destroy_mpi_enreg,fftpac,getcprj,getghc
!!      init_distribfft_seq,init_hamiltonian,initmpi_seq,initylmg,kpgsph
!!      load_k_hamiltonian,load_spin_hamiltonian,metric,mkffnl,mkkin,mkkpg
!!      pawcprj_alloc,pawcprj_free,pawcprj_reorder,transgrid,wrtout,xheev
!!      xheevx,xhegv,xhegvx,xmpi_barrier
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine ks_ddiago(Diago_ctl,nband_k,nfftc,mgfftc,ngfftc,natom,&
& typat,nfftf,nspinor,nspden,nsppol,Pawtab,Pawfgr,Paw_ij,&
& Psps,rprimd,vtrial,xred,onband_diago,eig_ene,eig_vec,Cprj_k,comm,ierr,&
& Electronpositron) ! Optional arguments

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_profiling_abi
 use m_xmpi
 use m_errors
 use m_hamiltonian

 use m_abilasi,           only : xheev, xhegv, xheevx, xhegvx
 use m_electronpositron,  only : electronpositron_type
 use m_fftcore,           only : kpgsph
 use m_mpinfo,            only : destroy_mpi_enreg
 use m_pawtab,            only : pawtab_type
 use m_paw_ij,            only : paw_ij_type
 use m_pawcprj,           only : pawcprj_type, pawcprj_alloc, pawcprj_free, &
&                                pawcprj_reorder
 use m_pawfgr,            only : pawfgr_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ks_ddiago'
 use interfaces_14_hidewrite
 use interfaces_41_geometry
 use interfaces_51_manage_mpi
 use interfaces_53_ffts
 use interfaces_56_recipspace
 use interfaces_65_paw
 use interfaces_66_nonlocal
 use interfaces_66_wfs
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mgfftc,natom,comm,nband_k
 integer,intent(in) :: nfftf,nsppol,nspden,nspinor,nfftc
 integer,intent(out) :: ierr,onband_diago
 type(pseudopotential_type),intent(in) :: Psps
 type(pawfgr_type),intent(in) :: Pawfgr
 type(ddiago_ctl_type),intent(in) :: Diago_ctl
!arrays
 integer,intent(in) :: typat(natom)
 integer,intent(in) :: ngfftc(18)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(inout) :: vtrial(nfftf,nspden)
 real(dp),intent(in) :: xred(3,natom)
 real(dp),pointer :: eig_ene(:),eig_vec(:,:,:)
 type(pawcprj_type),pointer :: Cprj_k(:,:)
 type(pawtab_type),intent(in) :: Pawtab(Psps%ntypat*Psps%usepaw)
 type(paw_ij_type),intent(in) :: Paw_ij(natom*Psps%usepaw)
 type(electronpositron_type),optional,pointer :: Electronpositron

!Local variables-------------------------------
!scalars
 integer,parameter :: mkmem_=1,tim_getghc=4,paral_kgb=0
 integer :: cprj_choice,cpopt,dimffnl,ib,ider,idir,isppol,npw_k
 integer :: ikg,master,istwf_k,exchn2n3d,prtvol
 integer :: jj,n1,n2,n3,n4,n5,n6,negv,nkpg,nprocs,npw_k_test
 integer :: my_rank,optder
 integer :: ispden,ndat,type_calc,sij_opt,igsp2,cplex_ghg
 integer :: iband,iorder_cprj,ibs1,ibs2
 real(dp) :: cfact,ucvol,ecutsm,effmass,lambda,size_mat,ecut
 logical :: do_full_diago
 character(len=50) :: jobz,range
 character(len=80) :: frmt1,frmt2
 character(len=10) :: stag(2)
 character(len=500) :: msg
 type(MPI_type) :: MPI_enreg_seq
 type(gs_hamiltonian_type) :: gs_hamk
!arrays
 integer :: nloalg(3)
 integer,allocatable :: kg_k(:,:)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3),kptns_(3,1)
 real(dp) :: kpoint(3),ylmgr_dum(1,1,1)
 real(dp) :: rhodum(1)
 real(dp),allocatable :: ph3d(:,:,:),pwave(:,:)
 real(dp),allocatable :: ffnl(:,:,:,:),kinpw(:),kpg_k(:,:)
 real(dp),allocatable :: vlocal(:,:,:,:),ylm_k(:,:),dum_ylm_gr_k(:,:,:),vlocal_tmp(:,:,:)
 real(dp),allocatable :: ghc(:,:),gvnlc(:,:),gsc(:,:)
 real(dp),allocatable :: ghg_mat(:,:,:),gtg_mat(:,:,:)
 real(dp),allocatable :: cgrvtrial(:,:)
 real(dp),pointer :: cwavef(:,:)
 type(pawcprj_type),allocatable :: Cwaveprj(:,:)

! *********************************************************************

 DBG_ENTER("COLL")

 nprocs  = xmpi_comm_size(comm)
 my_rank = xmpi_comm_rank(comm)
 master=0

 if (nprocs>1) then
   MSG_WARNING(" ks_ddiago not supported in parallel. Running in sequential.")
 end if

 call initmpi_seq(MPI_enreg_seq) ! Fake MPI_type for sequential part.
 call init_distribfft_seq(MPI_enreg_seq%distribfft,'c',ngfftc(2),ngfftc(3),'all')
 if (Pawfgr%usefinegrid/=0) then
   call init_distribfft_seq(MPI_enreg_seq%distribfft,'f',Pawfgr%ngfft(2),pawfgr%ngfft(3),'all')
 end if

 isppol  = Diago_ctl%isppol
 kpoint  = Diago_ctl%kpoint
 istwf_k = Diago_ctl%istwf_k
!% nband_k = Diago_ctl%nband_k
 npw_k   = Diago_ctl%npw_k
 nloalg  = Diago_ctl%nloalg
 ecut    = Diago_ctl%ecut
 ecutsm  = Diago_ctl%ecutsm
 effmass = Diago_ctl%effmass
 prtvol  = Diago_ctl%prtvol

 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 if (nsppol==1) stag=(/'          ','          '/)
 if (nsppol==2) stag=(/'SPIN UP:  ','SPIN DOWN:'/)
!
!The coarse FFT mesh.
 n1=ngfftc(1); n2=ngfftc(2); n3=ngfftc(3)
 n4=ngfftc(4); n5=ngfftc(5); n6=ngfftc(6)
!
!====================
!=== Check input ====
!====================
 ierr=0
!
!* istwfk must be 1 for each k-point
 if (istwf_k/=1) then
   write(msg,'(7a)')&
&   ' istwfk/=1 not allowed:',ch10,&
&   ' States output not programmed for time-reversal symmetry.',ch10,&
&   ' Action : change istwfk in input file (put it to 1 for all kpt).',ch10,&
&   ' Program does not stop but _KSS file will not be created...'
   MSG_WARNING(msg)
   ierr=ierr+1
 end if

 if (paral_kgb/=0) then
   write(msg,'(3a)')&
&   ' paral_kgb/=0 not allowed:',ch10,&
&   ' Program does not stop but _KSS file will not be created...'
   MSG_WARNING(msg)
   ierr=ierr+1
 end if

 if (ierr/=0) RETURN ! Houston we have a problem!
!
!Initialize the Hamiltonian datatype on the coarse FFT mesh.
 if (PRESENT(Electronpositron)) then
   call init_hamiltonian(gs_hamk,Psps,pawtab,nspinor,nspden,natom,typat,xred,nfftc,mgfftc,ngfftc,&
&   rprimd,nloalg,usecprj=0,Electronpositron=Electronpositron)
 else
   call init_hamiltonian(gs_hamk,Psps,pawtab,nspinor,nspden,natom,typat,xred,nfftc,mgfftc,ngfftc,&
&   rprimd,nloalg,usecprj=0)
 end if

!Check on the number of stored bands.
 if (nband_k==-1.or.nband_k>=npw_k*nspinor) then
   onband_diago=npw_k*nspinor
   write(msg,'(4a)')ch10,&
&   ' Since the number of bands to be computed was (-1) or',ch10,&
&   ' too large, it has been set to the max. value npw_k*nspinor. '
   call wrtout(std_out,msg,'COLL')
 else
   onband_diago=nband_k
 end if

!do_full_diago = (onband_diago==npw_k*nspinor)
 do_full_diago = Diago_ctl%do_full_diago

 if (do_full_diago) then
   write(msg,'(6a)')ch10,&
&   ' Since the number of bands to be computed',ch10,&
&   ' is equal to the nb of G-vectors found for this k-pt,',ch10,&
&   ' the program will perform complete diagonalizations.'
 else
   write(msg,'(6a)')ch10,&
&   ' Since the number of bands to be computed',ch10,&
&   ' is less than the number of G-vectors found,',ch10,&
&   ' the program will perform partial diagonalizations.'
 end if
 if (prtvol>0) then
   call wrtout(std_out,msg,'COLL')
 end if
!
!* Set up local potential vlocal with proper dimensioning, from vtrial.
!* Select spin component of interest if nspden<=2 as nvloc==1, for nspden==4, nvloc==4
!* option=2: vtrial(n1*n2*n3,ispden) --> vlocal(nd1,nd2,nd3) real case

!nvloc=1; if (nspden==4) nvloc=4
 ABI_MALLOC(vlocal,(n4,n5,n6,gs_hamk%nvloc))

 if (nspden/=4)then
   if (Psps%usepaw==0.or.Pawfgr%usefinegrid==0) then
     call fftpac(isppol,MPI_enreg_seq,nspden,n1,n2,n3,n4,n5,n6,ngfftc,vtrial,vlocal,2)
   else ! Move from fine to coarse FFT mesh (PAW)
     ABI_MALLOC(cgrvtrial,(nfftc,nspden))
     call transgrid(1,MPI_enreg_seq,nspden,-1,0,0,paral_kgb,Pawfgr,rhodum,rhodum,cgrvtrial,vtrial)
     call fftpac(isppol,MPI_enreg_seq,nspden,n1,n2,n3,n4,n5,n6,ngfftc,cgrvtrial,vlocal,2)
     ABI_FREE(cgrvtrial)
   end if
 else
   ABI_MALLOC(vlocal_tmp,(n4,n5,n6))
   if (Psps%usepaw==0.or.Pawfgr%usefinegrid==0) then
     do ispden=1,nspden
       call fftpac(ispden,MPI_enreg_seq,nspden,n1,n2,n3,n4,n5,n6,ngfftc,vtrial,vlocal_tmp,2)
       vlocal(:,:,:,ispden)=vlocal_tmp(:,:,:)
     end do
   else ! Move from fine to coarse FFT mesh (PAW)
     ABI_MALLOC(cgrvtrial,(nfftc,nspden))
     call transgrid(1,MPI_enreg_seq,nspden,-1,0,0,paral_kgb,Pawfgr,rhodum,rhodum,cgrvtrial,vtrial)
     do ispden=1,nspden
       call fftpac(ispden,MPI_enreg_seq,nspden,n1,n2,n3,n4,n5,n6,ngfftc,cgrvtrial,vlocal_tmp,2)
       vlocal(:,:,:,ispden)=vlocal_tmp(:,:,:)
     end do
     ABI_FREE(cgrvtrial)
   end if
   ABI_FREE(vlocal_tmp)
 end if

!Continue to initialize the Hamiltonian (spin-dependent part)
 call load_spin_hamiltonian(gs_hamk,isppol,paw_ij=Paw_ij,vlocal=vlocal)
!
!* Calculate G-vectors, for this k-point.
!* Count the number of planewaves as a check.
 exchn2n3d=0; ikg=0
 ABI_MALLOC(kg_k,(3,npw_k))

 call kpgsph(ecut,exchn2n3d,gmet,ikg,0,istwf_k,kg_k,kpoint,0,MPI_enreg_seq,0,npw_k_test)
 ABI_CHECK(npw_k_test==npw_k,"npw_k_test/=npw_k")

 call kpgsph(ecut,exchn2n3d,gmet,ikg,0,istwf_k,kg_k,kpoint,mkmem_,MPI_enreg_seq,npw_k,npw_k_test)

!========================
!==== Kinetic energy ====
!========================
 ABI_MALLOC(kinpw,(npw_k))
! call mkkin(ecut,ecutsm,effmass,gmet,kg_k,kinpw,kpoint,npw_k)
 call mkkin(ecut,ecutsm,effmass,gmet,kg_k,kinpw,kpoint,npw_k,0,0)
!
!================================
!==== Non-local form factors ====
!================================
 ABI_MALLOC(ylm_k,(npw_k,Psps%mpsang**2*Psps%useylm))

 if (Psps%useylm==1) then
   optder=0
   ABI_MALLOC(dum_ylm_gr_k,(npw_k,3+6*(optder/2),Psps%mpsang**2))
   kptns_(:,1) = kpoint

!  Here mband is not used if paral_compil_kpt=0
   call initylmg(gprimd,kg_k,kptns_,mkmem_,MPI_enreg_seq,Psps%mpsang,npw_k,(/nband_k/),1,&
&   (/npw_k/),1,optder,rprimd,ylm_k,dum_ylm_gr_k)

   ABI_FREE(dum_ylm_gr_k)
 end if

!Compute (k+G) vectors (only if useylm=1)
 nkpg=3*nloalg(3)
 ABI_MALLOC(kpg_k,(npw_k,nkpg))
 if (nkpg>0) then
   call mkkpg(kg_k,kpg_k,kpoint,nkpg,npw_k)
 end if

!Compute nonlocal form factors ffnl at all (k+G):
 idir=0; ider=0; dimffnl=1+ider ! Now the derivative is not needed anymore.
 ABI_MALLOC(ffnl,(npw_k,dimffnl,Psps%lmnmax,Psps%ntypat))

 call mkffnl(Psps%dimekb,dimffnl,Psps%ekb,ffnl,Psps%ffspl,gmet,gprimd,ider,idir,Psps%indlmn,&
& kg_k,kpg_k,kpoint,Psps%lmnmax,Psps%lnmax,Psps%mpsang,Psps%mqgrid_ff,nkpg,npw_k,&
& Psps%ntypat,Psps%pspso,Psps%qgrid_ff,rmet,Psps%usepaw,Psps%useylm,ylm_k,ylmgr_dum)

 ABI_FREE(ylm_k)

!Load k-dependent part in the Hamiltonian datastructure
 ABI_MALLOC(ph3d,(2,npw_k,gs_hamk%matblk))
 call load_k_hamiltonian(gs_hamk,kpt_k=kpoint,istwf_k=istwf_k,npw_k=npw_k,kinpw_k=kinpw,&
& kg_k=kg_k,kpg_k=kpg_k,ffnl_k=ffnl,ph3d_k=ph3d,&
& compute_ph3d=.true.,compute_gbound=.true.)

!Prepare the call to getghc.
 ndat=1; lambda=zero; type_calc=0         ! For applying the whole Hamiltonian
 sij_opt=0; if (Psps%usepaw==1) sij_opt=1 ! For PAW, <k+G|1+S|k+G"> is also needed.

 cpopt=-1    ! If cpopt=-1, <p_lmn|in> (and derivatives) are computed here (and not saved)
 if (Psps%usepaw==1.and..FALSE.) then ! TODO Calculate <p_lmn|k+G>.
   cpopt = 0  ! <p_lmn|in> are computed here and saved
 end if

 ABI_MALLOC(ghc  ,(2,npw_k*nspinor*ndat))
 ABI_MALLOC(gvnlc,(2,npw_k*nspinor*ndat))
 ABI_MALLOC(gsc  ,(2,npw_k*nspinor*ndat*(sij_opt+1)/2))

 cplex_ghg=2
 size_mat = cplex_ghg*(npw_k*nspinor)**2*dp*b2Mb
!; if (Psps%usepaw==1) size_mat=two*size_mat
!write(msg,'(a,f0.3,a)')" Memory required by the Hamiltonian matrix: ",size_mat," [Mb]."
!call wrtout(std_out,msg,"COLL")

 write(msg,'(a,f0.3,a)')" Out-of-memory in ghg_mat. Memory required by the Hamiltonian matrix: ",size_mat," [Mb]."
 ABI_STAT_MALLOC(ghg_mat,(cplex_ghg,npw_k*nspinor,npw_k*nspinor), ierr)
 ABI_CHECK(ierr==0, msg)

 write(msg,'(a,f0.3,a)')" Out-of-memory in gtg_mat. Memory required by the PAW overlap operator: ",size_mat," [Mb]."
 ABI_STAT_MALLOC(gtg_mat,(cplex_ghg,npw_k*nspinor,npw_k*nspinor*Psps%usepaw), ierr)
 ABI_CHECK(ierr==0, msg)

 ABI_DT_MALLOC(Cwaveprj,(natom,nspinor*(1+cpopt)*gs_hamk%usepaw))
 if (cpopt==0) then  ! Cwaveprj is ordered, see nonlop_ylm.
   call pawcprj_alloc(Cwaveprj,0,gs_hamk%dimcprj)
 end if

 ABI_MALLOC(pwave,(2,npw_k*nspinor))
 pwave=zero ! Initialize plane-wave array:

 if (prtvol>0) then
   call wrtout(std_out,' Calculating <G|H|G''> elements','PERS')
 end if

 do igsp2=1,npw_k*nspinor ! Loop over the |beta,G''> component.

   pwave(1,igsp2)=one      ! Get <:|H|beta,G''> and <:|T_{PAW}|beta,G''>

   call getghc(cpopt,pwave,Cwaveprj,ghc,gsc,gs_hamk,gvnlc,lambda,MPI_enreg_seq,ndat,&
&   prtvol,sij_opt,tim_getghc,type_calc)

!  Fill the upper triangle.
   ghg_mat(:,1:igsp2,igsp2) = ghc(:,1:igsp2)
   if (Psps%usepaw==1) gtg_mat(:,1:igsp2,igsp2) = gsc(:,1:igsp2)

   pwave(1,igsp2)=zero ! Reset the |G,beta> component that has been treated.
 end do

!Free workspace memory allocated so far.
 ABI_FREE(pwave)
 ABI_FREE(kinpw)
 ABI_FREE(vlocal)
 ABI_FREE(ghc)
 ABI_FREE(gvnlc)
 ABI_FREE(gsc)

 if (Psps%usepaw==1.and.cpopt==0) then
   call pawcprj_free(Cwaveprj)
 end if
 ABI_DT_FREE(Cwaveprj)
!
!===========================================
!=== Diagonalization of <G|H|G''> matrix ===
!===========================================
!
!*** Allocate the pointers ***
 ABI_MALLOC(eig_ene,(onband_diago))
 ABI_MALLOC(eig_vec,(cplex_ghg,npw_k*nspinor,onband_diago))

 jobz =Diago_ctl%jobz  !jobz="Vectors"

 if (do_full_diago) then ! * Complete diagonalization

   write(msg,'(2a,3es16.8,3x,3a,i5)')ch10,&
&   ' Begin complete diagonalization for kpt= ',kpoint(:),stag(isppol),ch10,&
&   ' - Size of mat.=',npw_k*nspinor
   if (prtvol>0) then
     call wrtout(std_out,msg,'PERS')
   end if

   if (Psps%usepaw==0) then
     call xheev(  jobz,"Upper",cplex_ghg,npw_k*nspinor,ghg_mat,eig_ene)
   else
     call xhegv(1,jobz,"Upper",cplex_ghg,npw_k*nspinor,ghg_mat,gtg_mat,eig_ene)
   end if

   eig_vec(:,:,:)=  ghg_mat

 else ! * Partial diagonalization

   range=Diago_ctl%range !range="Irange"

   write(msg,'(2a,3es16.8,3a,i5,a,i5)')ch10,&
&   ' Begin partial diagonalization for kpt= ',kpoint,stag(isppol),ch10,&
&   ' - Size of mat.=',npw_k*nspinor,' - # bnds=',onband_diago
   if (prtvol>0) then
     call wrtout(std_out,msg,'PERS')
   end if

   if (Psps%usepaw==0) then
     call xheevx(jobz,range,"Upper",cplex_ghg,npw_k*nspinor,ghg_mat,zero,zero,&
&     1,onband_diago,-tol8,negv,eig_ene,eig_vec,npw_k*nspinor)
   else
     call xhegvx(1,jobz,range,"Upper",cplex_ghg,npw_k*nspinor,ghg_mat,gtg_mat,zero,zero,&
&     1,onband_diago,-tol8,negv,eig_ene,eig_vec,npw_k*nspinor)
   end if

 end if

 ABI_FREE(ghg_mat)
 ABI_FREE(gtg_mat)
!
!========================================================
!==== Calculate <Proj_i|Cnk> from output eigenstates ====
!========================================================
 if (Psps%usepaw==1) then

   iorder_cprj=1 !  Ordered (order does change wrt input file); will be changed later
   ABI_DT_MALLOC(Cprj_k,(natom,nspinor*onband_diago))
   call pawcprj_alloc(Cprj_k,0,gs_hamk%dimcprj)

   idir=0; cprj_choice=1  ! Only projected wave functions.

   do iband=1,onband_diago
     ibs1=nspinor*(iband-1)+1
     ibs2=ibs1; if (nspinor==2) ibs2=ibs2+1
     cwavef => eig_vec(1:2,1:npw_k,iband)

     call getcprj(cprj_choice,0,cwavef,Cprj_k(:,ibs1:ibs2),&
&     gs_hamk%ffnl_k,idir,gs_hamk%indlmn,gs_hamk%istwf_k,gs_hamk%kg_k,&
&     gs_hamk%kpg_k,gs_hamk%kpt_k,gs_hamk%lmnmax,gs_hamk%mgfft,MPI_enreg_seq,&
&     gs_hamk%natom,gs_hamk%nattyp,gs_hamk%ngfft,gs_hamk%nloalg,gs_hamk%npw_k,gs_hamk%nspinor,&
&     gs_hamk%ntypat,gs_hamk%phkxred,gs_hamk%ph1d,gs_hamk%ph3d_k,gs_hamk%ucvol,gs_hamk%useylm)
   end do

!  Reorder the cprj (order is now the same as in input file)
   call pawcprj_reorder(Cprj_k,gs_hamk%atindx1)

!  deallocate(cwavef)
 end if !usepaw==1

 if (prtvol>0) then ! Write eigenvalues.
   cfact=Ha_eV ; frmt1='(8x,9(1x,f7.2))' ; frmt2='(8x,9(1x,f7.2))'
   write(msg,'(a,3es16.8,3x,a)')' Eigenvalues in eV for kpt= ',kpoint,stag(isppol)
   call wrtout(std_out,msg,'COLL')

   write(msg,frmt1)(eig_ene(ib)*cfact,ib=1,MIN(9,onband_diago))
   call wrtout(std_out,msg,'COLL')
   if (onband_diago>9) then
     do jj=10,onband_diago,9
       write(msg,frmt2) (eig_ene(ib)*cfact,ib=jj,MIN(jj+8,onband_diago))
       call wrtout(std_out,msg,'COLL')
     end do
   end if
 end if

 ABI_FREE(kpg_k)
 ABI_FREE(kg_k)
 ABI_FREE(ph3d)
 ABI_FREE(ffnl)

 call destroy_mpi_enreg(MPI_enreg_seq)
 call destroy_hamiltonian(gs_hamk)

 call xmpi_barrier(comm)

 DBG_EXIT("COLL")

end subroutine ks_ddiago
!!***
