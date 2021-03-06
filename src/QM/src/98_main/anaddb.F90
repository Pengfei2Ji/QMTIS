!{\src2tex{textfont=tt}}
!!****p* ABINIT/anaddb
!! NAME
!! anaddb
!!
!! FUNCTION
!! Main routine for analysis of the interatomic force constants and associated properties.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2016 ABINIT group (XG,DCA,JCC,CL,XW)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  (main routine)
!!
!! OUTPUT
!!  (main routine)
!!
!! PARENTS
!!
!! CHILDREN
!!      abi_io_redirect,abimem_init,abinit_doctor,anaddb_dtset_free,anaddb_init
!!      asria_calc,asria_corr,asrprs,crystal_free,ddb_diel,ddb_elast,ddb_free
!!      ddb_from_file,ddb_getdims,ddb_internalstr,ddb_piezo,dfpt_phfrq
!!      dfpt_prtph,dfpt_symph,elast_ncwrite,electrooptic,elphon,flush_unit
!!      gtblk9,gtdyn9,harmonic_thermo,herald,ifc_free,ifc_init,ifc_outphbtrap
!!      ifc_print,instrng,int2char4,inupper,invars9,isfile,mkherm,mkphbs
!!      mkphdos,outvars_anaddb,phdos_free,phdos_ncwrite,phdos_print
!!      phdos_print_debye,ramansus,relaxpol,thmeig,timein,wrtout,xmpi_bcast
!!      xmpi_end,xmpi_init,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

program anaddb

 use defs_basis
 use m_build_info
 use m_xmpi
 use m_xomp
 use m_profiling_abi
 use m_errors
 use m_ifc
 use m_ddb
 use m_phonons        
 use iso_c_binding
 use m_nctk
#ifdef HAVE_TRIO_NETCDF
 use netcdf
#endif

 use m_dfpt_io,        only : elast_ncwrite
 use m_io_tools,       only : open_file, flush_unit, num_opened_units, show_units
 use m_fstrings,       only : int2char4, itoa, sjoin, strcat
 use m_numeric_tools,  only : mkherm
 use m_time ,          only : asctime
 use m_anaddb_dataset, only : anaddb_dataset_type, anaddb_dtset_free, outvars_anaddb, invars9
 use m_crystal,        only : crystal_t, crystal_free
 use m_crystal_io,     only : crystal_ncwrite
 use m_dynmat,         only : asria_calc,asria_corr, chneu9, asrprs, gtdyn9

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'anaddb'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_42_parser
 use interfaces_51_manage_mpi
 use interfaces_72_response
 use interfaces_77_ddb
!End of the abilint section

 implicit none

!Arguments -----------------------------------

!Local variables-------------------------------
! Set array dimensions
 integer :: msym !  msym =maximum number of symmetry elements in space group
!Define input and output unit numbers (some are defined in defs_basis -all should be there ...):
 integer,parameter :: ddbun=2,master=0 ! FIXME: these should not be reserved unit numbers!
 integer :: dimekb,dims,comm,iatom,iblok,iblok_stress,idir,ii,index
 integer :: ierr,iphl2,lenstr,lmnmax,mband,mtyp,mpert,msize,natom,nblok,nblok2
 integer :: nkpt,nph2l,nsym,ntypat,option,rftyp,usepaw,nproc,my_rank
 logical :: iam_master
 integer :: rfelfd(4),rfphon(4),rfstrs(4),ngqpt_coarse(3)
 integer,allocatable :: d2flg(:)
 real(dp) :: etotal,tcpu,tcpui,twall,twalli
 real(dp),target :: dielt(3,3) 
 real(dp) :: compl(6,6),compl_clamped(6,6),compl_stress(6,6)
 real(dp) :: dielt_rlx(3,3),elast(6,6),elast_clamped(6,6),elast_stress(6,6)
 real(dp) :: epsinf(3,3),red_ptot(3),pel(3)
 real(dp) :: piezo(6,3),qphnrm(3),qphon(3,3),strten(6),tsec(2)
 real(dp),allocatable :: d2asr(:,:,:,:,:),d2cart(:,:),dchide(:,:,:)
 real(dp),allocatable :: dchidt(:,:,:,:),displ(:),eigval(:,:)
 real(dp),allocatable :: eigvec(:,:,:,:,:),fact_oscstr(:,:,:),instrain(:,:)
 real(dp),allocatable :: fred(:,:),lst(:),phfrq(:)
 real(dp),allocatable :: rsus(:,:,:)
 real(dp),allocatable :: singular(:),uinvers(:,:), vtinvers(:,:)
 real(dp),target,allocatable :: zeff(:,:,:)
 real(dp),allocatable :: d2asr_res(:,:,:,:,:)
 character(len=10) :: procstr
 character(len=24) :: codename
 character(len=24) :: start_datetime
 character(len=strlen) :: string
 character(len=fnlen) :: filnam(7),elph_base_name,tmpfilename,phdos_fname,ec_fname
 character(len=500) :: message
 type(anaddb_dataset_type) :: inp
 type(phonon_dos_type) :: Phdos
 type(ifc_type) :: Ifc,Ifc_coarse
 type(ddb_type) :: ddb
 type(crystal_t) :: Crystal
#ifdef HAVE_TRIO_NETCDF
 integer :: phdos_ncid, ana_ncid, ec_ncid, ncerr
 integer :: na_dir_varid,na_phmodes_varid
#endif

!******************************************************************

 ! Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm=xmpi_world)

 ! Initialize MPI
 call xmpi_init()

 ! MPI variables
 comm = xmpi_world
 nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 iam_master = (my_rank == master)

!Initialize memory profiling if it is activated !if a full abimem.mocc report is desired, 
!set the argument of abimem_init to "2" instead of "0"
!note that abimem.mocc files can easily be multiple GB in size so don't use this option normally
#ifdef HAVE_MEM_PROFILING
 call abimem_init(0)
#endif

 ! Initialisation of the timing
 call timein(tcpui,twalli)

 if (iam_master) then
   codename='ANADDB'//repeat(' ',18)
   call herald(codename,abinit_version,std_out)
 end if

 start_datetime = asctime()

 ! Initialise the code : write heading, and read names of files.
 if (iam_master) then
   call anaddb_init(filnam)
 end if
 call xmpi_bcast (filnam, master, comm, ierr)

 ! make log file for non-master procs
 if (.not. iam_master) then
   call int2char4(my_rank, procstr)
   ABI_CHECK((procstr(1:1)/='#'),'Bug: string length too short!')
   tmpfilename = trim(filnam(2)) // "_LOG_P" // trim(procstr)
   if (open_file(tmpfilename, message, unit=std_out, form="formatted", action="write") /= 0) then
     MSG_ERROR(message)
   end if
 end if 

!******************************************************************

 call timein(tcpu,twall)
 write(message, '(a,f11.3,a,f11.3,a)' )'-begin at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
 call wrtout(std_out,message,'COLL')

 ! Must read natom from the DDB before being able to allocate some arrays needed for invars9
 call ddb_getdims(dimekb,filnam(3),lmnmax,mband,mtyp,msym,natom,nblok,nkpt,ntypat,ddbun,usepaw,DDB_VERSION,comm)

 mpert=natom+6
 msize=3*mpert*3*mpert; if (mtyp==3) msize=msize*3*mpert

 ! Read the input file, and store the information in a long string of characters
 ! strlen from defs_basis module
 option=1
 if (iam_master) then
   call instrng (filnam(1),lenstr,option,strlen,string)

   !To make case-insensitive, map characters to upper case:
   call inupper(string(1:lenstr))
 end if

 call xmpi_bcast(string,master, comm, ierr)
 call xmpi_bcast(lenstr,master, comm, ierr)

 ! Read the inputs
 call invars9(inp,lenstr,natom,string)

 ! Open output files and ab_out (might change its name if needed)
 ! MJV 1/2010 : now output file is open, but filnam(2) continues unmodified
 ! so the other output files are overwritten instead of accumulating.
 if (iam_master) then
   tmpfilename = filnam(2)
   call isfile(tmpfilename,'new')
   if (open_file(tmpfilename,message,unit=ab_out,form='formatted',status='new') /= 0) then
     MSG_ERROR(message)
   end if
   rewind (unit=ab_out)
   call herald(codename,abinit_version,ab_out)

   ! Echo the inputs to console and main output file
   call outvars_anaddb(inp,std_out)
   call outvars_anaddb(inp,ab_out)
 else
   ab_out = dev_null
 end if

 nph2l=inp%nph2l
 ABI_ALLOCATE(lst,(nph2l))

!******************************************************************

 ! Read the DDB information, also perform some checks, and symmetrize partially the DDB
 write(message, '(a,a)' )' read the DDB information and perform some checks',ch10
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')

 call timein(tcpu,twall)
 write(message, '(a,f11.3,a,f11.3,a,a)' )'-begin at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec',ch10
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')

 ABI_ALLOCATE(instrain,(3*natom,6))
 ABI_ALLOCATE(d2cart,(2,msize))

 call ddb_from_file(ddb,filnam(3),inp%brav,natom,inp%natifc,inp%atifc,Crystal,comm)
 nsym = Crystal%nsym

 ABI_ALLOCATE(displ,(2*3*natom*3*natom))
 ABI_ALLOCATE(eigval,(3,natom))
 ABI_ALLOCATE(eigvec,(2,3,natom,3,natom))
 ABI_ALLOCATE(phfrq,(3*natom))
 ABI_ALLOCATE(zeff,(3,3,natom))

 ! Open the netcdf file that will contain the anaddb results
 if (iam_master) then
#ifdef HAVE_TRIO_NETCDF
   NCF_CHECK_MSG(nctk_open_create(ana_ncid, "anaddb.nc", xmpi_comm_self), "Creating anaddb.nc")
   NCF_CHECK(nctk_def_basedims(ana_ncid))
   NCF_CHECK(nctk_defnwrite_ivars(ana_ncid, ["anaddb_version"], [1]))
   NCF_CHECK(crystal_ncwrite(crystal, ana_ncid))
#endif
 end if

!**********************************************************************
!**********************************************************************

 ! Acoustic Sum Rule
 ! In case the interatomic forces are not calculated, the
 ! ASR-correction (d2asr) has to be determined here from the Dynamical matrix at Gamma.

 !acorr = ddb_make_asrq0corr(ddb, asr, rftyp, xcart) result(acorr)
 !call asrq0corr_free(acorr)

 ABI_ALLOCATE(d2asr,(2,3,natom,3,natom))
 d2asr = zero

 ! Pre allocate array used if asr in [3,4]
 dims=3*natom*(3*natom-1)/2
 ABI_CALLOC(uinvers,(1:dims,1:dims))
 ABI_CALLOC(vtinvers,(1:dims,1:dims))
 ABI_CALLOC(singular,(1:dims))

 if (inp%ifcflag==0 .or. inp%instrflag/=0 .or. inp%elaflag/=0) then
   ! Find the Gamma block in the DDB (no need for E-field entries)
   qphon(:,1)=zero
   qphnrm(1)=zero
   rfphon(1:2)=1
   rfelfd(:)=0
   rfstrs(:)=0
   rftyp=inp%rfmeth

   call gtblk9(ddb,iblok,qphon,qphnrm,rfphon,rfelfd,rfstrs,rftyp)

   d2asr = zero
   if (iblok /=0) then
     select case (inp%asr)
     case (0)
       continue 

     case (1,2)
       call asria_calc(inp%asr,d2asr,ddb%val(:,:,iblok),ddb%mpert,ddb%natom)

     case (3,4)
       ! Rotational invariance for 1D and 0D systems
       call asrprs(inp%asr,1,3,uinvers,vtinvers,singular,ddb%val(:,:,iblok),ddb%mpert,ddb%natom,Crystal%xcart)

     case (5)
       ! d2cart is a temp variable here
       d2cart = ddb%val(:,:,iblok)
       ! calculate diagonal correction
       call asria_calc(2,d2asr,d2cart,ddb%mpert,ddb%natom)
       ! apply diagonal correction
       call asria_corr(2,d2asr,d2cart,ddb%mpert,ddb%natom)
       ! hermitianize
       call mkherm(d2cart,3*ddb%mpert)
       ! remove remaining ASR rupture due to Hermitianization
       ABI_ALLOCATE(d2asr_res,(2,3,ddb%natom,3,ddb%natom))
       call asria_calc(inp%asr,d2asr_res,d2cart,ddb%mpert,ddb%natom)
       ! full correction is sum of both
       d2asr = d2asr + d2asr_res
       ABI_DEALLOCATE(d2asr_res)

     case default
       write(message,'(a,i0)')"Wrong value for asr: ",inp%asr
       MSG_ERROR(message)
     end select
   end if
 end if

 ! Get Dielectric Tensor and Effective Charges
 ! (initialized to one_3D and zero if the derivatives are not available in the DDB file)
 iblok = ddb_get_dielt_zeff(ddb,crystal,inp%rfmeth,inp%chneut,inp%selectz,dielt,zeff)

 if (my_rank == master) then
#ifdef HAVE_TRIO_NETCDF
   ! TODO: Cartesian or reduced?
   ncerr = nctk_def_arrays(ana_ncid, [&
   nctkarr_t('emacro_cart', "dp", 'number_of_cartesian_directions, number_of_cartesian_directions'),&
   nctkarr_t('becs_cart', "dp", "number_of_cartesian_directions, number_of_cartesian_directions, number_of_atoms")],&
   defmode=.True.)
   NCF_CHECK(ncerr)

   NCF_CHECK(nctk_set_datamode(ana_ncid))
   NCF_CHECK(nf90_put_var(ana_ncid, nctk_idname(ana_ncid, 'emacro_cart'), dielt))
   NCF_CHECK(nf90_put_var(ana_ncid, nctk_idname(ana_ncid, 'becs_cart'), zeff))
#endif
 end if

 ! Structural response at fixed polarization
 if (inp%polflag == 1) then
   ABI_ALLOCATE(d2flg,(msize))

   if(iblok/=0)then
     ! Save the second-order derivatives
     d2cart(1:2,1:msize) = ddb%val(1:2,1:msize,iblok)
     d2flg(1:msize) = ddb%flg(1:msize,iblok)

   else 
     ! the gamma blok has not been found
     if (inp%relaxat==0 .and. inp%relaxstr==0) then
       ! The gamma blok is not needed
       d2cart(1:2,1:msize)=zero
       d2flg(1:msize)=1
     else 
       ! There is a problem !
       write(message, '(7a)' )&
&       'The dynamical matrix at Gamma is needed, in order to perform ',ch10,&
&       "relaxation at constant polarisation (Na Sai's method)",ch10,&
&       'However, this was not found in the DDB.',ch10,&
&       'Action: complete your DDB with the dynamical matrix at Gamma.'
       MSG_ERROR(message)
     end if
   end if ! iblok not found

   ! Extract the block with the total energy
   if (ddb_get_etotal(ddb,etotal) == 0) then
     MSG_ERROR("DDB file does not contain GS etotal")
   end if

   ! Extract the block with the gradients
   ABI_ALLOCATE(fred,(3,natom))
   qphon(:,:) = zero
   qphnrm(:) = zero
   rfphon(:) = 0
   rfstrs(:) = 0
   rftyp = 4
   rfelfd(:) = 2
   if (inp%relaxat == 1) rfphon(:) = 1
   if (inp%relaxstr == 1) rfstrs(:) = 3

   call gtblk9(ddb,iblok,qphon,qphnrm,rfphon,rfelfd,rfstrs,rftyp)

   if (inp%relaxat == 1) then
     index = 0
     do iatom = 1, natom
       do idir = 1, 3
         index = index + 1
         fred(idir,iatom) = ddb%val(1,index,iblok)
       end do
     end do
   end if

   pel(1:3) = ddb%val(1,3*natom+4:3*natom+6,iblok)

   if (inp%relaxstr == 1) then
     index = 3*natom + 6
     do ii = 1, 6
       index = index + 1
       strten(ii) = ddb%val(1,index,iblok)
     end do
   end if

   ! when called from here red_ptot is not set! So set it to zero
   red_ptot(:)=zero

   call relaxpol(Crystal,d2flg,d2cart,etotal,fred,inp%iatfix,&
&   ab_out,inp%istrfix,mpert,msize,inp%natfix,natom,&
&   inp%nstrfix,pel,red_ptot,inp%relaxat,inp%relaxstr,&
&   strten,inp%targetpol,usepaw)

   ABI_DEALLOCATE(fred)
   ABI_DEALLOCATE(d2flg)
 end if

!***************************************************************************

 ! Compute non-linear optical susceptibilities and
 ! First-order change in the linear dielectric susceptibility induced by an atomic displacement
 if (inp%nlflag > 0) then
   ABI_ALLOCATE(dchide,(3,3,3))
   ABI_ALLOCATE(dchidt,(natom,3,3,3))

   if (ddb_get_dchidet(ddb,inp%ramansr,dchide,dchidt) == 0) then
     message = "Cannot find block corresponding to non-linear optical susceptibilities in DDB file"
     MSG_ERROR(message)
   end if
 end if ! nlflag

!**********************************************************************
! Interatomic Forces Calculation
!**********************************************************************
 if (inp%ifcflag ==1) then
   ! ifc to be calculated for interpolation
   write(message, '(a,a,(80a),a,a,a,a)' ) ch10,('=',ii=1,80),ch10,ch10,&
&   ' Calculation of the interatomic forces ',ch10
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

   call timein(tcpu,twall)
   write(message, '(a,f11.3,a,f11.3,a)' )'-begin at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

   if (inp%qrefine > 1) then
     ! Gaal-Nagy's algorithm in PRB <b>73</b> 014117.

     ! Build the IFCs using the coarse q-mesh.
     ngqpt_coarse(1:3) = inp%ngqpt(1:3)/inp%qrefine
     call ifc_init(Ifc_coarse,Crystal,ddb,&
&     inp%brav,inp%asr,inp%symdynmat,inp%dipdip,inp%rfmeth,ngqpt_coarse,inp%nqshft,inp%q1shft,dielt,zeff,&
&     inp%nsphere,inp%rifcsph,inp%prtsrlr,inp%enunit,prtfreq=.True.)

     ! And now use the coarse q-mesh to fill the entries in dynmat(q) 
     ! on the dense q-mesh that cannot be obtained from the DDB file.
     call ifc_init(Ifc,Crystal,ddb,&
&     inp%brav,inp%asr,inp%symdynmat,inp%dipdip,inp%rfmeth,inp%ngqpt(1:3),inp%nqshft,inp%q1shft,dielt,zeff,&
&     inp%nsphere,inp%rifcsph,inp%prtsrlr,inp%enunit,prtfreq=.True.,Ifc_coarse=Ifc_coarse)
     call ifc_free(Ifc_coarse)

   else
     call ifc_init(Ifc,Crystal,ddb,&
&     inp%brav,inp%asr,inp%symdynmat,inp%dipdip,inp%rfmeth,inp%ngqpt(1:3),inp%nqshft,inp%q1shft,dielt,zeff,&
&     inp%nsphere,inp%rifcsph,inp%prtsrlr,inp%enunit,prtfreq=.True.)
   end if

   !Print analysis of the real-space interatomic force constants
   if(inp%ifcout/=0)then
#ifdef HAVE_TRIO_NETCDF
     call ifc_print(Ifc,dielt,zeff,inp%ifcana,inp%atifc,inp%ifcout,inp%prt_ifc,ncid=ana_ncid)
#else
     call ifc_print(Ifc,dielt,zeff,inp%ifcana,inp%atifc,inp%ifcout,inp%prt_ifc)
#endif
   end if
 end if

!**********************************************************************
!**********************************************************************

!Short-Range/Long-Range decomposition of the phonon frequencies
!if (inp%prtsrlr == 1) then
!call wrtout(std_out,' anaddb    : start of the SR/LR decomposition ','COLL')
!end if

!**********************************************************************

!Electron-phonon section
 if (inp%elphflag == 1) then
   call elphon(inp,Crystal,Ifc,filnam)
 end if

!**********************************************************************

!Phonon density of states calculation, Start if interatomic forces have been calculated
 if (inp%ifcflag==1 .and. any(inp%prtdos==[1,2])) then
   write(message,'(a,(80a),4a)')ch10,('=',ii=1,80),ch10,ch10,&
&   ' Calculation of phonon density of states ',ch10
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')

   call mkphdos(Phdos,Crystal,Ifc, inp%prtdos,inp%dosdeltae,inp%dossmear, inp%ng2qpt, inp%q2shft)

   phdos_fname = TRIM(filnam(2))//"_PHDOS"
   call phdos_print(Phdos,phdos_fname)
   call phdos_print_debye(Phdos, Crystal%ucvol)

   if (iam_master) then
#ifdef HAVE_TRIO_NETCDF
     ncerr = nctk_open_create(phdos_ncid, trim(phdos_fname)//".nc", xmpi_comm_self)
     NCF_CHECK_MSG(ncerr, "Creating PHDOS.nc file")
     NCF_CHECK(crystal_ncwrite(Crystal, phdos_ncid))
     call phdos_ncwrite(Phdos, phdos_ncid) 
     NCF_CHECK(nf90_close(phdos_ncid))
#endif
   end if

   call phdos_free(Phdos)
 end if

 if (iam_master.and.inp%ifcflag==1 .and. inp%outboltztrap==1) then
   call ifc_outphbtrap(Ifc,Crystal,inp%ng2qpt,1,inp%q2shft,filnam(2))
 end if

 ! Phonon density of states and thermodynamical properties calculation
 ! Start if interatomic forces and thermal flags are on
 if (inp%ifcflag==1 .and. any(inp%thmflag==[1,2])) then

   write(message, '(a,(80a),a,a,a,a,a,a,a,a)' ) ch10,('=',ii=1,80),ch10,ch10,&
&   ' Calculation of phonon density of states, ',ch10,&
&   '    thermodynamical properties, ',ch10,&
&   '    and Debye-Waller factors.',ch10
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')

   call timein(tcpu,twall)
   write(message, '(a,f11.3,a,f11.3,a)' )'-begin at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

   if (inp%thmflag==1) then
     call harmonic_thermo(Ifc,Crystal,ddb%amu,inp,ab_out,filnam(2),tcpui,twalli,comm)

   else if (inp%thmflag==2) then
     write(message, '(a,(80a),a,a,a,a)' ) ch10,('=',ii=1,80),ch10,&
&     ch10,' Entering thm9 routine with thmflag=2 ',ch10
     call wrtout(std_out,message,'COLL')
     call wrtout(ab_out,message,'COLL')

     call harmonic_thermo(Ifc,Crystal,ddb%amu,inp,ab_out,filnam(2),tcpui,twalli,comm,&
&     thmflag=inp%thmflag)
   end if
 end if

!**********************************************************************

 ! Now treat the first list of vectors (without non-analyticities)
 call mkphbs(Ifc,Crystal,inp,ddb,d2asr,filnam(2),singular,tcpui,twalli,uinvers,vtinvers,zeff,comm)

!***********************************************************************

!Test thmeig
! MG: FIXME
! =====================================================================================================
! COULD SOMEONE PLEASE CLEAN THE INTERFACE OF THMEIG? WHY DO WE HAVE SO MANY VARIABLES WITH INTENT(OUT)
! =====================================================================================================
! real(dp),intent(out) :: ucvol !new
! integer,intent(inout) :: natom,nkpt,nsym,ntypat,occopt,nblok2 !new in ==> inout
! integer,intent(out) :: symrel(3,3,msym)
! integer,intent(out) :: indsym(4,nsym,natom),symrec(3,3,msym) !new
! integer,intent(inout) :: typat(natom),atifc(natom)! new in ==> inout
! real(dp),intent(out) :: zion(ntypat),tnons(3,msym),gmet(3,3) !new
! real(dp),intent(out) :: gprim(3,3),rmet(3,3),xcart(3,natom) !new
! real(dp),intent(inout) :: acell(3),amu(ntypat),rprim(3,3),xred(3,natom)! new in ==> inout
 if (inp%thmflag>=3 .and. inp%thmflag<=8) then

!  Obtain the number of bloks contained in this file.
   call ddb_getdims(dimekb,filnam(5),lmnmax,mband,mtyp,msym,natom,nblok2,nkpt,ntypat,ddbun,usepaw,DDB_VERSION,comm)

   !write(std_out,*)'Entering thmeig: '
   elph_base_name=trim(filnam(2))//"_ep"
   call thmeig(inp%a2fsmear,ddb%acell,ddb%amu,inp,d2asr,&
&   elph_base_name,mband,mpert,msize,natom,nkpt,inp%ntemper,&
&   ntypat,ddb%rprim,inp%telphint,inp%temperinc,&
&   inp%tempermin,inp%thmflag,Crystal%typat,Crystal%xred,&
&   ddb,ddbun,dimekb,filnam(5),ab_out,& !new
&  lmnmax,msym,nblok2,Crystal%nsym,ddb%occopt,Crystal%symrel,Crystal%tnons,usepaw,Crystal%zion,& !new
&  Crystal%symrec,inp%natifc,Crystal%gmet,ddb%gprim,Crystal%indsym,Crystal%rmet,inp%atifc,& !new
&  Crystal%ucvol,Crystal%xcart,comm) !new
 end if

!**********************************************************************

!Now treat the second list of vectors (only at the Gamma point,
!but can include non-analyticities), as well as the frequency-dependent dielectric tensor

 if (inp%nlflag > 0) then
   ABI_ALLOCATE(rsus,(3*natom,3,3))
 end if
 ABI_ALLOCATE(fact_oscstr,(2,3,3*natom))

 if (nph2l/=0 .or. inp%dieflag==1) then

   write(message, '(a,(80a),a,a,a,a)' ) ch10,('=',ii=1,80),ch10,&
&   ch10,' Treat the second list of vectors ',ch10
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

   call timein(tcpu,twall)
   write(message, '(a,f11.3,a,f11.3,a)' )&
&   '-begin at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

   ! Before examining every direction or the dielectric tensor, generates the dynamical matrix at gamma
   qphon(:,1)=zero
   qphnrm(1)=zero

   ! Generation of the dynamical matrix in cartesian coordinates
   if (inp%ifcflag==1) then

     ! Get d2cart using the interatomic forces and the
     ! long-range coulomb interaction through Ewald summation
     call gtdyn9(ddb%acell,Ifc%atmfrc,dielt,inp%dipdip,&
&     Ifc%dyewq0,d2cart,Crystal%gmet,ddb%gprim,mpert,natom,&
&     Ifc%nrpt,qphnrm(1),qphon,Crystal%rmet,ddb%rprim,Ifc%rpt,&
&     Ifc%trans,Crystal%ucvol,Ifc%wghatm,Crystal%xred,zeff)

   else if (inp%ifcflag==0) then

     ! Look after the information in the DDB
     rfphon(1:2)=1
     rfelfd(1:2)=2
     rfstrs(1:2)=0
     rftyp=inp%rfmeth
     call gtblk9(ddb,iblok,qphon,qphnrm,rfphon,rfelfd,rfstrs,rftyp)

     ! Copy the dynamical matrix in d2cart
     d2cart(:,1:msize)=ddb%val(:,:,iblok)

     ! Eventually impose the acoustic sum rule
     call asria_corr(inp%asr,d2asr,d2cart,mpert,natom)
   end if ! end of the generation of the dynamical matrix at gamma.

   if (nph2l/=0) then

     if (my_rank == master) then
#ifdef HAVE_TRIO_NETCDF
       NCF_CHECK(nctk_def_basedims(ana_ncid, defmode=.True.))

       ncerr = nctk_def_dims(ana_ncid, [&
       nctkdim_t("number_of_non_analytical_directions", nph2l), nctkdim_t('number_of_phonon_modes', 3*natom)],&
       defmode=.True.)
       NCF_CHECK(ncerr)

       ncerr = nctk_def_arrays(ana_ncid, [&
       nctkarr_t('non_analytical_directions', "dp", "number_of_cartesian_directions, number_of_non_analytical_directions"),&
       nctkarr_t('non_analytical_phonon_modes', "dp", "number_of_phonon_modes, number_of_non_analytical_directions")],&
       defmode=.True.)
       NCF_CHECK(ncerr)

       NCF_CHECK(nctk_set_datamode(ana_ncid))

       NCF_CHECK(nf90_inq_varid(ana_ncid, "non_analytical_directions", na_dir_varid))
       NCF_CHECK(nf90_put_var(ana_ncid,na_dir_varid,inp%qph2l))
#endif
     end if

     ! Examine every wavevector of this list
     do iphl2=1,nph2l

       ! Initialisation of the phonon wavevector
       qphon(:,1)=inp%qph2l(:,iphl2)
       qphnrm(1)=inp%qnrml2(iphl2)

       ! Calculation of the eigenvectors and eigenvalues of the dynamical matrix
       call dfpt_phfrq(ddb%amu,displ,d2cart,eigval,eigvec,Crystal%indsym,&
&       mpert,msym,natom,nsym,ntypat,phfrq,qphnrm(1),qphon,Crystal%rprimd,inp%symdynmat,&
&       Crystal%symrel,Crystal%symafm,Crystal%typat,Crystal%ucvol)

       ! Write the phonon frequencies
       call dfpt_prtph(displ,inp%eivec,inp%enunit,ab_out,natom,phfrq,qphnrm(1),qphon)

       if (my_rank == master) then
#ifdef HAVE_TRIO_NETCDF
         NCF_CHECK(nf90_inq_varid(ana_ncid, "non_analytical_phonon_modes", na_phmodes_varid))
         NCF_CHECK(nf90_put_var(ana_ncid,na_phmodes_varid,phfrq*Ha_eV,start=[1, iphl2], count=[3*natom, 1]))
#endif
       end if

       ! Determine the symmetries of the phonon modes at Gamma
       if (sum(abs(qphon(:,1)))<DDB_QTOL) then
         call dfpt_symph(ab_out,ddb%acell,eigvec,Crystal%indsym,natom,nsym,phfrq,ddb%rprim,Crystal%symrel)
       end if

       ! Write Raman susceptibilities
       if (inp%nlflag == 1) then
         call ramansus(d2cart,dchide,dchidt,displ,mpert,natom,phfrq,qphon,qphnrm(1),rsus,Crystal%ucvol)
       end if

       ! Prepare the evaluation of the Lyddane-Sachs-Teller relation
       if(inp%dieflag==1 .and. natom>1)then
         lst(iphl2)=zero
         ! The fourth mode should have positive frequency, otherwise,
         ! there is an instability, and the LST relationship should not be evaluated
         if(phfrq(4)>tol6)then
           do ii=4,3*natom
             lst(iphl2)=lst(iphl2)+2*log(phfrq(ii))
           end do
         end if
       end if

     end do ! iphl2
   end if ! nph2l/=0

   ! The frequency-dependent dielectric tensor (and oscillator strength).
   if (inp%dieflag==1)then
     write(message, '(a,a,a,a,a,a)' )&
&     ' the frequency-dependent dielectric tensor (and also once more',ch10,&
&     ' the phonons at gamma - without non-analytic part )',ch10,ch10,&
&     ' The frequency-dependent dielectric tensor'
     call wrtout(std_out,message,'COLL')

     ! Initialisation of the phonon wavevector
     qphon(:,1)=zero
     qphnrm(1)=zero

     ! Calculation of the eigenvectors and eigenvalues of the dynamical matrix
     call dfpt_phfrq(ddb%amu,displ,d2cart,eigval,eigvec,Crystal%indsym,&
&     mpert,msym,natom,nsym,ntypat,phfrq,qphnrm(1),qphon,&
&     Crystal%rprimd,inp%symdynmat,Crystal%symrel,Crystal%symafm,Crystal%typat,Crystal%ucvol)

     ! Write the phonon frequencies (not to ab_out, however)
     call dfpt_prtph(displ,0,inp%enunit,-1,natom,phfrq,qphnrm(1),qphon)

     ! Evaluation of the oscillator strengths and frequency-dependent dielectric tensor.
     call ddb_diel(Crystal,ddb%amu,inp,dielt_rlx,displ,d2cart,epsinf,fact_oscstr,&
&     ab_out,lst,mpert,natom,nph2l,phfrq)

     ! write(std_out,*)'after ddb_diel, dielt_rlx(:,:)=',dielt_rlx(:,:)
   end if

   ! If the electronic dielectric tensor only is needed...
   if (inp%dieflag==2.or.inp%dieflag==3.or. inp%dieflag==4) then
!    Everything is already in place...
     call ddb_diel(Crystal,ddb%amu,inp,dielt_rlx,displ,d2cart,epsinf,fact_oscstr,&
&     ab_out,lst,mpert,natom,nph2l,phfrq)
   end if

 end if ! End the condition of either nph2l/=0  or  dieflag==1

!**********************************************************************

!In case nph2l was equal to 0, the electronic dielectric tensor has to be computed independently.

 if (inp%dieflag==2 .and. inp%nph2l==0) then
   call wrtout(std_out,'nph2l=0, so compute the electronic dielectric tensor independently','COLL')
   ! Look after the second derivative matrix at gamma in the DDB
   ! Note that the information on the dielectric tensor is completely
   ! independent of the interatomic force constant calculation
   qphon(:,1)=zero
   qphnrm(1)=zero
   rfphon(1:2)=0
   rfelfd(1:2)=2
   rfstrs(1:2)=0
   rftyp=inp%rfmeth

   call gtblk9(ddb,iblok,qphon,qphnrm,rfphon,rfelfd,rfstrs,rftyp)

   d2cart(:,1:msize)=ddb%val(:,:,iblok)

   ! Print the electronic dielectric tensor
   call ddb_diel(Crystal,ddb%amu,inp,dielt_rlx,displ,d2cart,epsinf,fact_oscstr,ab_out,lst,mpert,natom,nph2l,phfrq)
 end if

!**********************************************************************

 ! Compute the electrooptic tensor
 if (inp%nlflag == 1) then
   ! In case dieflag = 2, recompute phonon frequencies and eigenvectors without non-analyticity
   if (inp%dieflag == 2) then
     qphon(:,1)=zero
     qphnrm(1)=zero
     call dfpt_phfrq(ddb%amu,displ,d2cart,eigval,eigvec,Crystal%indsym,&
&     mpert,msym,natom,nsym,ntypat,phfrq,qphnrm(1),qphon,&
&     Crystal%rprimd,inp%symdynmat,Crystal%symrel,Crystal%symafm,Crystal%typat,Crystal%ucvol)
   end if

   rsus = zero
   call ramansus(d2cart,dchide,dchidt,displ,mpert,natom,phfrq(1),qphon,qphnrm(1),rsus,Crystal%ucvol)

   call electrooptic(dchide,inp%dieflag,epsinf,fact_oscstr,natom,phfrq,inp%prtmbm,rsus,Crystal%ucvol)
 end if ! condition on nlflag

 ABI_DEALLOCATE(fact_oscstr)
 if (inp%nlflag > 0)  then
   ABI_DEALLOCATE(dchide)
   ABI_DEALLOCATE(dchidt)
   ABI_DEALLOCATE(rsus)
 end if

!**********************************************************************

! Here treating the internal strain tensors at Gamma point
 if (inp%instrflag/=0) then

   write(message, '(a,a,(80a),a,a,a,a)') ch10,('=',ii=1,80),ch10,ch10,&
&   ' Calculation of the internal-strain  tensor',ch10
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

   call timein(tcpu,twall)
   write(message,'(a,f11.3,a,f11.3,a)')'-begin at tcpu',tcpu-tcpui,'   and twall',twall-twalli,'sec'
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

   if (inp%instrflag==1) then
     call wrtout(std_out,'instrflag=1, so extract the internal strain constant from the 2DTE','COLL')

     ! looking after the no. of blok that contains the internal strain tensor
     qphon(:,1)=zero
     qphnrm(1)=zero
     rfphon(1:2)=0
     rfelfd(1:2)=0
     rfstrs(1:2)=3
     rftyp=inp%rfmeth
     call gtblk9(ddb,iblok,qphon,qphnrm,rfphon,rfelfd,rfstrs,rftyp)
     ! then print the internal stain tensor
     call ddb_internalstr(inp%asr,ddb%val,d2asr,iblok,instrain,ab_out,mpert,natom,ddb%nblok)
   end if
 end if !end the part for internal strain

!**********************************************************************

!here treating the elastic tensors at Gamma Point
 if (inp%elaflag/=0) then
   write(message, '(a,a,(80a),a,a,a,a,a,a)') ch10,('=',ii=1,80),ch10,ch10,&
&   ' Calculation of the elastic and compliances tensor (Voigt notation)',ch10
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

   call timein(tcpu,twall)
   write(message,'(a,f11.3,a,f11.3,a)')'-begin at tcpu',tcpu-tcpui,'   and twall',twall-twalli,'sec'
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

   if (any(inp%elaflag == [1,2,3,4,5])) then
     call wrtout(std_out,'so extract the elastic constant from the 2DTE','COLL')

     ! look after the blok no. that contains the stress tensor
     qphon(:,1)=zero
     qphnrm(1)=zero
     rfphon(1:2)=0
     rfelfd(1:2)=0
     rfstrs(1:2)=0
     rftyp=4

     call gtblk9(ddb,iblok,qphon,qphnrm,rfphon,rfelfd,rfstrs,rftyp)
     iblok_stress=iblok

     ! look after the blok no.iblok that contains the elastic tensor
     qphon(:,1)=zero
     qphnrm(1)=zero
     rfphon(1:2)=0
     rfelfd(1:2)=0
     rfstrs(1:2)=3

     ! for both diagonal and shear parts
     rftyp=inp%rfmeth
     call gtblk9(ddb,iblok,qphon,qphnrm,rfphon,rfelfd,rfstrs,rftyp)

     ! print the elastic tensor
     call ddb_elast(inp,ddb%val,compl,compl_clamped,compl_stress,d2asr,&
&     elast,elast_clamped,elast_stress,iblok,iblok_stress,&
&     instrain,ab_out,mpert,natom,ddb%nblok,Crystal%ucvol)
     ec_fname = TRIM(filnam(2))//"_EC.nc"
#ifdef HAVE_TRIO_NETCDF
     if (iam_master) then
       ncerr = nctk_open_create(ec_ncid, ec_fname, xmpi_comm_self) 
       NCF_CHECK_MSG(ncerr, "Creating EC.nc file")
       NCF_CHECK(crystal_ncwrite(Crystal, ec_ncid))
       call elast_ncwrite(compl,compl_clamped,compl_stress,elast,elast_clamped,elast_stress,ec_ncid)
       NCF_CHECK(nf90_close(ec_ncid))
     end if
#endif
   end if
 end if !ending the part for elastic tensors

!**********************************************************************

!here treating the piezoelectric tensor at Gamma Point
 if (inp%piezoflag/=0 .or. inp%dieflag==4 .or. inp%elaflag==4) then
   write(message, '(a,a,(80a),a,a,a,a,a)') ch10,('=',ii=1,80),ch10,ch10,&
&   ' Calculation of the tensor related to piezoelectric effetc',ch10,&
&   '  (Elastic indices in Voigt notation)',ch10
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

   call timein(tcpu,twall)
   write(message,'(a,f11.3,a,f11.3,a)')'-begin at tcpu',tcpu-tcpui,'   and twall',twall-twalli,'sec'
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

   if (any(inp%piezoflag == [1,2,3,4,5,6,7]) .or. inp%dieflag==4 .or.inp%elaflag==4) then
     call wrtout(std_out,'extract the piezoelectric constant from the 2DTE','COLL')

     ! looking for the gamma point block
     qphon(:,1)=zero
     qphnrm(1)=zero
     rfphon(1:2)=0
     rfelfd(1:2)=0
     rfstrs(1:2)=3
     ! for both diagonal and shear parts
     rftyp=inp%rfmeth

     call gtblk9(ddb,iblok,qphon,qphnrm,rfphon,rfelfd,rfstrs,rftyp)

     ! then print out the piezoelectric constants
     call ddb_piezo(inp,ddb%val,dielt_rlx,elast,iblok,instrain,ab_out,mpert,natom,ddb%nblok,piezo,Crystal%ucvol)
   end if
 end if

!**********************************************************************

 ABI_DEALLOCATE(displ)
 ABI_DEALLOCATE(d2asr)
 ABI_DEALLOCATE(d2cart)
 ABI_DEALLOCATE(eigval)
 ABI_DEALLOCATE(eigvec)
 ABI_DEALLOCATE(lst)
 ABI_DEALLOCATE(phfrq)
 ABI_DEALLOCATE(zeff)
 ABI_DEALLOCATE(instrain)
 ABI_DEALLOCATE(uinvers)
 ABI_DEALLOCATE(vtinvers)
 ABI_DEALLOCATE(singular)

 call anaddb_dtset_free(inp)
 call ddb_free(ddb)
 call crystal_free(Crystal)
 call ifc_free(Ifc)

 ! Close files
 if (iam_master) then
#ifdef HAVE_TRIO_NETCDF
   NCF_CHECK(nf90_close(ana_ncid))
#endif
 end if

 call timein(tcpu,twall)
 tsec(1)=tcpu-tcpui
 tsec(2)=twall-twalli

 write(message, '(a,i4,a,f13.1,a,f13.1)' )' Proc.',my_rank,' individual time (sec): cpu=',tsec(1),'  wall=',tsec(2)
 call wrtout(std_out,message,"COLL")

 if (iam_master) then
   write(ab_out, '(a,a,a,i4,a,f13.1,a,f13.1)' )'-',ch10,&
&   '- Proc.',my_rank,' individual time (sec): cpu=',tsec(1),'  wall=',tsec(2)
 end if

 call xmpi_sum(tsec,comm,ierr)

 write(message, '(a,(80a),a,a,a,f11.3,a,f11.3,a,a,a,a)' ) ch10,&
& ('=',ii=1,80),ch10,ch10,&
& '+Total cpu time',tsec(1),&
& '  and wall time',tsec(2),' sec',ch10,ch10,&
& ' anaddb : the run completed succesfully.'
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')

 if (iam_master) then
   ! Write YAML document with the final summary.
   ! we use this doc to test whether the calculation is completed.
   write(std_out,"(a)")"--- !FinalSummary"
   write(std_out,"(a)")"program: anaddb"
   write(std_out,"(2a)")"version: ",trim(abinit_version)
   write(std_out,"(2a)")"start_datetime: ",start_datetime
   write(std_out,"(2a)")"end_datetime: ",asctime()
   write(std_out,"(a,f13.1)")"overall_cpu_time: ",tsec(1)
   write(std_out,"(a,f13.1)")"overall_wall_time: ",tsec(2)
   write(std_out,"(a,i0)")"mpi_procs: ",xmpi_comm_size(xmpi_world)
   write(std_out,"(a,i0)")"omp_threads: ",xomp_get_num_threads(open_parallel=.True.)
   !write(std_out,"(a,i0)")"num_warnings: ",nwarning
   !write(std_out,"(a,i0)")"num_comments: ",ncomment
   write(std_out,"(a)")"..."
   call flush_unit(std_out)
 end if

!Write information on file about the memory before ending mpi module, if memory profiling is enabled
 call abinit_doctor(filnam(2))

 call flush_unit(ab_out)
 call flush_unit(std_out)

 if (iam_master) close(ab_out)

 call xmpi_end()

 end program anaddb
!!***
