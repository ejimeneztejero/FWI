!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!	This file contains the modules:
!!	(mod_parfile, mod_data_arrays)
!!	Author:	Clara Estela Jimenez Tejero. 
!!	License: GNU General Public License v3.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module mod_parfile

!! INFORMATION ABOUT: WAS PARAMETERS, MCS PARAMETERS, MODEL PARAMETERS, FWI PARAMETERS

!! OJO, se puede hacer multishooting con OBS y con Shots de Streamer?

implicit none

  !!unidades de lectura archivos (ojo con el num de OBS, Shots en MCS e iteraciones)
  integer, parameter :: unit_was=1000,unit_mcs=5000
  integer, parameter :: unit_misfit=20,unit_iter=30,unit_FWI=50
 
  INTEGER, parameter :: size_su_header = 60
  INTEGER(4), parameter :: byte_fldr = 9
  INTEGER(4), parameter :: byte_tracl = 1
  INTEGER(4), parameter :: byte_tracr = 5
  INTEGER(4), parameter :: byte_offset = 37
  INTEGER(4), parameter :: byte_scalco = 71
  INTEGER(4), parameter :: byte_sx = 73
  INTEGER(4), parameter :: byte_sy = 77
  
  INTEGER(4) :: byte_shotnumber_was
  INTEGER(4) :: byte_shotnumber_mcs

  integer :: inv_WAS,inv_MCS
  integer :: inv_source
  integer :: NumSS_was,NumSS_mcs,SS_was,SS_mcs
  integer :: NumShots_was,NumShots_mcs
  integer :: NumOBS,NumChannels	!!fuentes reales y OBS reales
  integer :: NumSou_MCS,NumSou_WAS,NumRec_MCS,NumRec_WAS
  integer :: seed_option
  integer :: niter
  integer :: nmodel
  integer :: nxmodel,nymodel 
  integer :: nxmodel_was 
  integer :: nxmodel_mcs

  integer :: nt_was,nt_mcs
  integer :: NumFreq,nFWI
  integer :: read_vp    !!ojo, un modelo de velocidad para was, y otro para mcs? o uno comun?

  integer :: endianness_data,endianness_machine
  integer :: save_txt,save_xyz,save_z,save_su
  integer :: vp_z,vp_xyz,vp_su
  integer :: offset_header
  integer :: maxbytes

  integer :: misfit_type,BP_type
  integer :: method,typef
  integer :: added_grid

  real :: offset_unit
  real :: drec
  real :: dshots_was,dshots_mcs
  real :: near_offset,dmodel
  real :: dt_was,dt_mcs
  real :: shot_depth_was,shot_depth_mcs
  real :: streamer_depth
  real :: streamer_length
  real :: added_space_model,water_velocity
  real :: tfin_was,tfin_mcs
  real :: water_depth
  real :: step_freq,freq_ricker
  real :: f_init,f_final
  real :: vpi,vpf
  real :: length_model  !!ojo
  real :: length_model_was  !!ojo
  real :: length_model_mcs  !!ojo

  character(len=500) :: folder_input_was,folder_input_mcs,folder_input_model
  character(len=500) :: folder_output
  character(len=500) :: folder_VEL,folder_GRAD,folder_DATA,folder_SOU
  character(len=500) :: su_file_FWI
  character(len=500) :: was_data,mcs_data
  character(len=500) :: nav_obs !! navegacion (posiciones OBS-WAS)
  character(len=500) :: nav_shot_was,nav_shot_mcs
  character(len=500) :: vp_file
  character(len=500) :: par_file

  character(len=500), allocatable :: su_file_was(:)
  character(len=500), allocatable :: su_file_mcs(:)

  contains

  subroutine read_parfile(rank,numtasks)

  implicit none
  include 'mpif.h'

  integer :: numtasks,rank,ierr,errcode,status(MPI_STATUS_SIZE)

  ! Input related variables
  character(len=200) :: buffer,label
  integer :: ps,icount,i,ifile,interval,nlines
  integer, parameter :: fh = 10
  integer :: ios = 0
  integer :: line = 0

  character(len=500) :: command0,command
  character(len=50) :: Str,access,form,num_split
  character(len=500) :: file_name,file_name2

  logical :: su_exist,nav_exist,vp_exist,data_exist
  logical :: input_exist, output_exist

  call MPI_COMM_SIZE(MPI_COMM_WORLD,numtasks,ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)

  ierr=0;

  access = 'STREAM'
  form = 'UNFORMATTED'

  icount = iargc()  ! The number does not includes the executable name, so if user passed one argument, 1 is returned.
  if ( icount.eq.1 ) then
	call getarg(1, par_file)	! The file name of the executable.
	if(rank.eq.0)write(*,*)'name par_file: ',trim(adjustl(par_file))
	file_name = trim(adjustl(par_file))
	INQUIRE(file=file_name,EXIST=su_exist)
	if(.NOT. su_exist)     then
  		if(rank.eq.0)write(*,*)'ERROR: Par_file named: ', trim(adjustl(par_file)),' does not exist'
        	if(rank.eq.0)call ascii_art(2)
		call MPI_barrier(MPI_COMM_WORLD,ierr)
		call MPI_Abort(MPI_COMM_WORLD, errcode, ierr)
		stop
  	endif
  endif

  su_file_FWI = 'null'

  nav_shot_was = 'null'
  nav_shot_mcs = 'null'
  nav_obs = 'null'
  vp_file = 'null'

  byte_shotnumber_was= byte_tracr
  byte_shotnumber_mcs= byte_fldr
  endianness_data=1;endianness_machine=0;
  offset_header=0;offset_unit=1;

  inv_source=0;
  inv_WAS=0;inv_MCS=0;
  seed_option=0;
  added_space_model=0!!en metros

  water_velocity=1500.
  dmodel=0;
  shot_depth_was=0;
  shot_depth_mcs=0;
  streamer_depth=0;
  streamer_length=0;
  f_init=1.;f_final=10.;
  save_txt=0;save_z=0;save_xyz=0;save_su=1;
  vp_z=0;vp_xyz=0;vp_su=0;
  misfit_type=1!!por defecto L2

  method=1
  typef=3  
  misfit_type=1
  BP_type=1

  tfin_was=-100.;
  tfin_mcs=-100.;
  nt_was=-100;dt_was=-100;
  nt_mcs=-100;dt_mcs=-100;
  NumOBS=0;NumChannels=0;NumShots_was=0;NumShots_mcs=0;
  NumSou_WAS=0;NumSou_MCS=0;
  NumRec_WAS=0;NumRec_MCS=0;
  NumSS_was=0;NumSS_mcs=0;
  SS_was=0;SS_mcs=0;
  read_vp=0
  vpi=1500.;vpf=4000;
  nmodel=1
  nxmodel=1
  nymodel=1
  nxmodel_was=1
  nxmodel_mcs=1
  near_offset=0  
  freq_ricker=10

  open(fh, file=par_file)

  ! ios is negative if an end of record condition is encountered or if
  ! an endfile condition was detected. It is positive if an error was
  ! detected.  ios is zero otherwise.

  do while (ios == 0)

     read(fh, '(A)', iostat=ios) buffer
     if (ios == 0) then
        line = line + 1

	! Verificar si la línea es un comentario
        if (buffer(1:1) .ne. '#') then

        ! Find the first instance of whitespace.  Split label and data.
        ps = scan(buffer,' ')
        label = buffer(1:ps)
        buffer = buffer(ps+1:)

        select case (label)

        case ('inv_source:')
           read(buffer, *, iostat=ios) inv_source
        case ('fi:')
           read(buffer, *, iostat=ios) f_init
        case ('ff:')
           read(buffer, *, iostat=ios) f_final
        case ('step_freq:')
           read(buffer, *, iostat=ios) step_freq
        case ('nFWI:')
           read(buffer, *, iostat=ios) nFWI
        case ('freq_ricker:')
           read(buffer, *, iostat=ios) freq_ricker
        case ('endianness_data:')
           read(buffer, *, iostat=ios) endianness_data
        case ('endianness_machine:')
           read(buffer, *, iostat=ios) endianness_machine
        case ('save_su:')
           read(buffer, *, iostat=ios) save_su
        case ('save_z:')
           read(buffer, *, iostat=ios) save_z
        case ('save_xyz:')
           read(buffer, *, iostat=ios) save_xyz
        case ('method:')
           read(buffer, *, iostat=ios) method
        case ('typef:')
           read(buffer, *, iostat=ios) typef
        case ('misfit_type:')
           read(buffer, *, iostat=ios) misfit_type
        case ('BP_type:')
           read(buffer, *, iostat=ios) BP_type
        case ('vp_file:')
           read(buffer, *, iostat=ios) vp_file
        case ('vp_su:')
           read(buffer, *, iostat=ios) vp_su
        case ('vp_z:')
           read(buffer, *, iostat=ios) vp_z
        case ('vp_xyz:')
           read(buffer, *, iostat=ios) vp_xyz
        case ('vpi:')
           read(buffer, *, iostat=ios) vpi
        case ('vpf:')
           read(buffer, *, iostat=ios) vpf
        case ('output_folder:')
           read(buffer, *, iostat=ios) folder_output
        case ('nxmodel:')
           read(buffer, *, iostat=ios) nmodel
        case ('nymodel:')
           read(buffer, *, iostat=ios) nymodel
        case ('dmodel:')
           read(buffer, *, iostat=ios) dmodel
        case ('water_velocity:')
           read(buffer, *, iostat=ios) water_velocity
        case ('seed_option:')
           read(buffer, *, iostat=ios) seed_option 
        case ('input_folder_model:')
	    read(buffer, *, iostat=ios) folder_input_model
        case ('inv_WAS:')
           read(buffer, *, iostat=ios) inv_WAS
        case ('inv_MCS:')
           read(buffer, *, iostat=ios) inv_MCS

!        case default
!           if(rank.eq.0)print *, 'Skipping invalid label at line', line

        end select

        end if  ! Fin del if para verificar comentarios
     end if
  end do

close(fh)


if(inv_WAS.ne.0)	then

ios=0
open(fh, file=par_file)

do while (ios == 0)
read(fh, '(A)', iostat=ios) buffer
if (ios == 0) then
        line = line + 1
	! Verificar si la línea es un comentario
        if (buffer(1:1) .ne. '#') then
        ! Find the first instance of whitespace.  Split label and data.
        ps = scan(buffer,' ')
        label = buffer(1:ps)
        buffer = buffer(ps+1:)

        select case (label)

	        case ('byte_shotnumber_was:')
	           read(buffer, *, iostat=ios) byte_shotnumber_was
	        case ('was_file_list:')
	           read(buffer, *, iostat=ios) was_data
	        case ('nav_shot_was:')
	           read(buffer, *, iostat=ios) nav_shot_was
	        case ('nav_obs:')
	           read(buffer, *, iostat=ios) nav_obs
	        case ('input_folder_was:')
	           read(buffer, *, iostat=ios) folder_input_was
	        case ('dt_was:')
		read(buffer, *, iostat=ios) dt_was
		case ('nt_was:')
	           read(buffer, *, iostat=ios) nt_was
	        case ('NumShots_was:')
	           read(buffer, *, iostat=ios) NumShots_was
	        case ('NumOBS:')
	           read(buffer, *, iostat=ios) NumOBS
	        case ('NumSS_was:')
	            read(buffer, *, iostat=ios) NumSS_was
	        case ('shot_depth_was:')
	           read(buffer, *, iostat=ios) shot_depth_was

end select
end if  ! Fin del if para verificar comentarios
end if
end do
close(fh)

endif


if(inv_MCS.ne.0)	then

ios=0
open(fh, file=par_file)

do while (ios == 0)

read(fh, '(A)', iostat=ios) buffer
if (ios == 0) then
        line = line + 1

	! Verificar si la línea es un comentario
        if (buffer(1:1) .ne. '#') then

        ! Find the first instance of whitespace.  Split label and data.
        ps = scan(buffer,' ')
        label = buffer(1:ps)
        buffer = buffer(ps+1:)

        select case (label)

	        case ('byte_shotnumber_mcs:')
	            read(buffer, *, iostat=ios) byte_shotnumber_mcs
	        case ('mcs_file_list:')
	           read(buffer, *, iostat=ios) mcs_data
	        case ('nav_shot_mcs:')
	           read(buffer, *, iostat=ios) nav_shot_mcs
	        case ('offset_header_mcs:')
	            read(buffer, *, iostat=ios) offset_header
	        case ('offset_unit_mcs:')
		    read(buffer, *, iostat=ios) offset_unit
        	case ('input_folder_mcs:')
        	    read(buffer, *, iostat=ios) folder_input_mcs
        	case ('dt_mcs:')
        	   read(buffer, *, iostat=ios) dt_mcs
        	case ('nt_mcs:')
        	   read(buffer, *, iostat=ios) nt_mcs
        	case ('NumShots_mcs:')
        	   read(buffer, *, iostat=ios) NumShots_mcs
        	case ('near_offset:')
        	   read(buffer, *, iostat=ios) near_offset
        	case ('NumRec:')
        	   read(buffer, *, iostat=ios) NumChannels
        	case ('NumSS_mcs:')
        	    read(buffer, *, iostat=ios) NumSS_mcs
        	case ('drec:')
        	   read(buffer, *, iostat=ios) drec
        	case ('shot_depth_mcs:')
        	    read(buffer, *, iostat=ios) shot_depth_mcs
        	case ('streamer_depth:')
        	   read(buffer, *, iostat=ios) streamer_depth

	case default

!if(rank.eq.0)print *, 'WARNING in file ',trim(adjustl(par_file)),': skipping invalid label at line', line

end select

end if  ! Fin del if para verificar comentarios
end if
end do
close(fh)

end if


if(inv_WAS.ne.0)folder_input_was = trim(adjustl(folder_input_was)) // '/'
if(inv_MCS.ne.0)folder_input_mcs = trim(adjustl(folder_input_mcs)) // '/'
folder_input_model = trim(adjustl(folder_input_model)) // '/'

folder_output = trim(adjustl(folder_output)) // '/'
folder_VEL=trim(adjustl(folder_output)) // "MODEL/"
folder_GRAD=trim(adjustl(folder_output)) // "GRADIENT/"
folder_DATA=trim(adjustl(folder_output)) // "DATA/"
folder_SOU=trim(adjustl(folder_output)) // "SOURCE/"

call MPI_barrier(MPI_COMM_WORLD,ierr)
if(rank.eq.0)	then
	command="mkdir " // trim(adjustl(folder_output))
	call system(command)
	command="mkdir " // trim(adjustl(folder_DATA))
	call system(command)
	command="mkdir " // trim(adjustl(folder_VEL))
	call system(command)
	command="mkdir " // trim(adjustl(folder_GRAD))
	call system(command)
	command="mkdir " // trim(adjustl(folder_SOU))
	call system(command)
	command="cp " // trim(adjustl(par_file)) // " " // trim(adjustl(folder_output))
	call system(command)
endif
call MPI_barrier(MPI_COMM_WORLD,ierr)

if(rank.eq.0)call warnings_errors(numtasks)

nxmodel=nmodel
length_model=(nxmodel-1)*dmodel

NumRec_WAS=NumShots_was	!! IN WAS DATA: RECIPROCITY. SHOTS ACT AS RECEIVERS:
NumSou_WAS=NumOBS	!! sources in WAS are OBSs

NumSou_MCS=NumShots_mcs	!! source in MCS are shotgathers
NumRec_MCS=NumChannels	!! Receivers are channels

if(NumSS_was.ne.0)SS_was=1
if(NumSS_mcs.ne.0)SS_mcs=1
if(NumSS_was.eq.0)NumSS_was=NumSou_WAS
if(NumSS_mcs.eq.0)NumSS_mcs=NumSou_MCS

tfin_was=(nt_was-1)*dt_was
tfin_mcs=(nt_mcs-1)*dt_mcs

streamer_length=drec*(NumRec_MCS-1)
added_grid=floor(added_space_model/dmodel)

NumFreq=1+(f_final-f_init)/step_freq

niter=NumFreq*nFWI

file_name=trim(adjustl(folder_output)) // "f1_f2_rIter_Iter.dat"
open(unit_iter,file=file_name,status='unknown')
if(rank.eq.0)write(unit_iter,*)"# f1, f2, relative iteration, iteration"

file_name=trim(adjustl(folder_output)) // 'Misfit_iter.dat'
open(unit_misfit,file=file_name,status='unknown')
if(rank.eq.0)write(unit_misfit,*)"# iteration, misfit"

if(inv_WAS.ne.0) then
    allocate(su_file_was(NumOBS))
    file_name= trim(adjustl(folder_input_was)) // trim(adjustl(was_data))
    open(unit=10,file=file_name,status='old')
    do i=1,NumOBS
        read(10,*)file_name2
        su_file_was(i) = trim(adjustl(file_name2))
    enddo
    close(10)
endif

if(inv_MCS.ne.0) then
    allocate(su_file_mcs(NumShots_mcs))
    file_name= trim(adjustl(folder_input_mcs)) // trim(adjustl(mcs_data))
    open(unit=10,file=file_name,status='old')
    do i=1,NumShots_mcs
        read(10,*)file_name2
        su_file_mcs(i) = trim(adjustl(file_name2))
    enddo
    close(10)
endif

maxbytes=1900000000
	
if(save_z.ne.0.or.save_xyz.ne.0)save_txt=1

file_name = trim(adjustl(folder_input_model)) // trim(adjustl(vp_file))
INQUIRE(FILE=file_name, EXIST=vp_exist)
if(vp_exist)read_vp=1

if(rank.eq.0)	then

	write(*,*)
	write(*,*)'*******************************************'
	write(*,*)'RUNNING PROGRAM FOR NEXT SET OF PARAMETERS:'
	write(*,*)'*******************************************'
        write(*,*)'Output folder: ',adjustl(trim(folder_output))
	write(*,*)'dmodel (m): ',dmodel
	write(*,*)'nxmmodel: ',nxmodel
	write(*,*)'nymodel: ',nymodel
	write(*,*)'method: ',method
	write(*,*)'typef: ',typef
	write(*,*)'misfit_type: ',misfit_type
	write(*,*)'BP_type: ',BP_type
	if(read_vp.eq.1)    then
        	write(*,*)'Input folder model: ',adjustl(trim(folder_input_model))
        	write(*,*)'Vp file: ',adjustl(trim(vp_file))
	endif
        if(read_vp.eq.0)write(*,*)'Vp is a gradient, between:',vpi,' and ',vpf
        if(inv_WAS.ne.0)    then
        	write(*,*)'Inversion parameters for WAS: '
        	write(*,*)'Input WAS folder: ',adjustl(trim(folder_input_was))
        	write(*,*)'NumSS_was: ',NumSS_was
        	write(*,*)'NumOBS: ',NumOBS
        	write(*,*)'NumShots WAS: ',NumShots_was
        	write(*,*)'name of navigation obs file: ',adjustl(trim(nav_obs))
        	write(*,*)'name of navigation shots WAS file: ',adjustl(trim(nav_shot_was))
        	write(*,*)'name of file with the name of the WAS-SU-data files: ',adjustl(trim(was_data))
	        write(*,*)'dt_was: ',dt_was
		write(*,*)'nt_was: ',nt_was
        endif
        if(inv_MCS.ne.0) then
        	write(*,*)'Inversion parameters for MCS: '
        	write(*,*)'Input MCS folder: ',adjustl(trim(folder_input_mcs))
        	write(*,*)'NumSS_mcs: ',NumSS_mcs
        	write(*,*)'NumChannels: ',NumChannels
        	write(*,*)'NumShots MCS: ',NumShots_mcs
        	write(*,*)'**name of navigation shots MCS file: ',adjustl(trim(nav_shot_mcs))
        	write(*,*)'**name of file with the name of the MCS (shotgathers) SU-data files: ',adjustl(trim(mcs_data))
	        write(*,*)'dt_mcs: ',dt_mcs
		write(*,*)'nt_mcs: ',nt_mcs
        endif

endif

end subroutine read_parfile

subroutine warnings_errors(numtasks)

implicit none
integer :: rank,numtasks

  ! Input related variables
  character(len=200) :: buffer,label
  integer :: ps,icount,i,ifile,interval,nlines
  integer, parameter :: fh = 10
  integer :: ios = 0
  integer :: line = 0
  integer :: error = 0
  character(len=500) :: command0,command
  character(len=50) :: Str,access,form,num_split
  character(len=500) :: file_name,file_name2

  logical :: su_exist,nav_exist,vp_exist,data_exist
  logical :: input_exist, output_exist

file_name = trim(adjustl(folder_input_model)) // trim(adjustl(vp_file))
INQUIRE(FILE=file_name, EXIST=vp_exist)
if(.NOT. vp_exist)	then
	write(*,*)"WARNING: vp_file not given or not found"
	write(*,*)"Therefore the initial Vp model will be set to a gradient with values, vpi, vpf: ", vpi, vpf
endif

if(endianness_data.ne.0.and.endianness_data.ne.1)	then
	write(*,*)'ERROR: endianness_data should be set to 0 (little endian)&
	 or 1 (big endian) in ',trim(adjustl(par_file))
	error=1
endif

if(endianness_machine.ne.0.and.endianness_machine.ne.1)	then
	write(*,*)'ERROR: endianness_machine should be set to 0 (little endian)&
	 or 1 (big endian) in ',trim(adjustl(par_file))
	error=1
endif

if(inv_was.ne.0)    then

    file_name = trim(adjustl(folder_input_was)) // trim(adjustl(was_data))
    INQUIRE(FILE=file_name, EXIST=data_exist)
    if(.NOT. data_exist)     then
        write(*,*)'ERROR: "file_data_list: " file not found'
        error=1
    else
        !!	OBS DATA FILE
        file_name= trim(adjustl(folder_input_was)) // trim(adjustl(was_data))
        open(unit=10,file=file_name,status='old')
        nlines=0
        do
            read(10,*, END=10)
            nlines = nlines + 1
        enddo
        10 close (10)

        if(nlines.lt.NumOBS)    then
            write(*,*)'ERROR: ', abs(nlines-NumOBS),' OBS data files are missing in was_file_list'
            error=1
        endif
    endif

    file_name = trim(adjustl(folder_input_was)) // trim(adjustl(nav_shot_was))
    INQUIRE(FILE=file_name, EXIST=nav_exist)
    if(.NOT. nav_exist)	then
        write(*,*)'ERROR: nav_shot_was not found'
        error=1
    endif
    file_name = trim(adjustl(folder_input_was)) // trim(adjustl(nav_obs))
    INQUIRE(FILE=file_name, EXIST=nav_exist)
    if(.NOT. nav_exist)	then
        write(*,*)'ERROR: nav_obs not found'
        error=1
    endif

if(dt_was.eq.0) then
	write(*,*) 'ERROR: Please give a value to dt (seconds) in ',trim(adjustl(par_file))
	error=1
endif

if(nt_was.eq.0) then
	write(*,*)'ERROR: Please give a value to nt in ',trim(adjustl(par_file))
	error=1
endif


endif


if(inv_MCS.ne.0)    then

    file_name = trim(adjustl(folder_input_mcs)) // trim(adjustl(mcs_data))
    INQUIRE(FILE=file_name, EXIST=data_exist)
    if(.NOT. data_exist)     then
        write(*,*)'ERROR: "file_data_list: " file not found'
        error=1
    else
        file_name= trim(adjustl(folder_input_mcs)) // trim(adjustl(mcs_data))
        open(unit=10,file=file_name,status='old')
        nlines=0
        do
            read(10,*, END=11)
            nlines = nlines + 1
        enddo
        11 close (10)

        if(nlines.lt.NumShots_mcs)    then
            write(*,*)'ERROR: ', abs(nlines-NumShots_mcs),' MCS data files are missing in mcs_file_list'
            error=1
        endif
    endif

    file_name = trim(adjustl(folder_input_mcs)) // trim(adjustl(nav_shot_mcs))
    INQUIRE(FILE=file_name, EXIST=nav_exist)
    if(.NOT. nav_exist)    then
        write(*,*)'ERROR: nav_shot_mcs not found: ',file_name
        error=1
    endif

    if(dt_mcs.eq.0) then
	write(*,*) 'ERROR: Please give a value to dt (seconds) in ',trim(adjustl(par_file))
	error=1
    endif

    if(nt_mcs.eq.0) then
	write(*,*)'ERROR: Please give a value to nt in ',trim(adjustl(par_file))
	error=1
    endif

endif

if(inv_WAS.ne.0.and.shot_depth_was.eq.0) then
	write(*,*)'ERROR: Please give a value to shot_depth_was (meters) in ',trim(adjustl(par_file))
	error=1
endif

if(inv_MCS.ne.0.and.shot_depth_mcs.eq.0) then
    write(*,*)'ERROR: Please give a value to shot_depth_mcs (meters) in ',trim(adjustl(par_file))
    error=1
endif

!if(inv_MCS.ne.0.and.inv_WAS.ne.0.and.NumOBS+NumShots_mcs.ne.numtasks) then
!	write(*,*)"Num CPU must be equal to num NumShots_mcs+NumOBS: ", NumShots_mcs+NumOBS
!	error=1
!endif
!if(inv_MCS.eq.0.and.inv_WAS.ne.0.and.NumOBS.ne.numtasks) then
!	write(*,*)"Num CPU must be equal to NumOBS: ", NumOBS
!	error=1
!endif
!if(inv_MCS.ne.0.and.inv_WAS.eq.0.and.NumShots_mcs.ne.numtasks) then
!	write(*,*)"Num CPU must be equal to NumShots_mcs: ", NumShots_mcs
!	error=1
!endif

if(error.eq.1) then
        call ascii_art(2)
        stop
endif

end subroutine warnings_errors

end module mod_parfile

module mod_data_arrays

implicit none

!!      WAS DATA
        integer :: shotID_was_1, shotID_was_n
        integer, allocatable :: shotID_nav_was(:),shotID_su_was(:,:),shotID_was_(:)
        integer, allocatable :: pos_bat_grid_was(:) !!position shot, no dimensions (points)
        integer, allocatable :: pos_shot_grid_was(:) !!position shot, no dimensions (points)
        real, allocatable :: pos_bat_was(:) !!position shot
        real, allocatable :: pos_shot_was(:) !!position shot

        integer, allocatable :: pos_xobs_grid(:) !!position obs x
        integer, allocatable :: pos_zobs_grid(:) !!position obs z
        real, allocatable :: pos_xobs(:) !!position obs x
        real, allocatable :: pos_zobs(:) !!position obs z
        real, allocatable    :: bat_model_was(:)
        integer, allocatable    :: bat_model_grid_was(:)
        integer, allocatable :: nxSou_WAS(:),nySou_WAS(:)
        integer, allocatable :: nxRec_WAS(:,:),nyRec_WAS(:,:)
        integer, allocatable :: sizeof_was(:)

!!      MCS DATA
        integer :: shotID_mcs_1, shotID_mcs_n
        integer, allocatable :: shotID_nav_mcs(:),shotID_su_mcs(:),shotID_mcs_(:)
        integer, allocatable :: pos_bat_grid_mcs(:) !!position shot, no dimensions (points)
        integer, allocatable :: pos_shot_grid_mcs(:) !!position shot, no dimensions (points)
        real, allocatable :: pos_bat_mcs(:) !!position shot
        real, allocatable :: pos_shot_mcs(:) !!position shot
        integer, allocatable :: pos_rec_grid(:,:) !!position channels in mcs x
        real, allocatable :: pos_rec(:,:) !!position mcs x
        real, allocatable :: MCS_raw_Source(:,:),MCS_raw_Data(:,:)
        real, allocatable    :: bat_model_mcs(:)
        real, allocatable    :: offset_su(:,:)
        integer, allocatable    :: bat_model_grid_mcs(:)
        integer, allocatable :: nxSou_MCS(:),nySou_MCS(:)
        integer, allocatable :: nxRec_MCS(:,:),nyRec_MCS(:,:)

!!      GENERAL
        real, allocatable :: Model_ini(:,:),lambda(:,:),rho(:,:),Freq(:)

contains

subroutine allocate_data_arrays()
use mod_parfile

implicit none

allocate(Freq(NumFreq))
Freq=0.

    if(inv_WAS.ne.0)    then

	allocate(sizeof_was(NumShots_was))
        allocate(shotID_nav_was(NumShots_was))
        allocate(shotID_su_was(NumShots_was,NumOBS),shotID_was_(NumShots_was))
        allocate(pos_shot_was(NumShots_was),pos_shot_grid_was(NumShots_was))
        allocate(pos_bat_was(NumShots_was),pos_bat_grid_was(NumShots_was))

        shotID_was_1=0;shotID_was_n=0;
        shotID_nav_was=0;
        shotID_su_was=0;shotID_was_=0;
        pos_shot_was=0;pos_shot_grid_was=0;
        pos_bat_was=0;pos_bat_grid_was=0;

        allocate(pos_xobs(NumOBS),pos_zobs(NumOBS))
        allocate(pos_xobs_grid(NumOBS),pos_zobs_grid(NumOBS))
        pos_xobs=0;pos_xobs_grid=0;
        pos_zobs=0;pos_zobs_grid=0;

!!      Geometry en WAS
        allocate(nxSou_WAS(NumSou_WAS),nySou_WAS(NumSou_WAS))
        allocate(nxRec_WAS(NumRec_WAS,NumSou_WAS),nyRec_WAS(NumRec_WAS,NumSou_WAS))
        nxSou_WAS=0;nySou_WAS=0;nxRec_WAS=0;nyRec_WAS=0;

    endif

        if(inv_MCS.ne.0)    then

            allocate(shotID_nav_mcs(NumShots_mcs))
            allocate(shotID_su_mcs(NumShots_mcs),shotID_mcs_(NumShots_mcs))
            allocate(pos_shot_mcs(NumShots_mcs),pos_shot_grid_mcs(NumShots_mcs))
            allocate(pos_bat_mcs(NumShots_mcs),pos_bat_grid_mcs(NumShots_mcs))

            shotID_mcs_1=0;shotID_mcs_n=0;
            shotID_nav_mcs=0;
            shotID_su_mcs=0;shotID_mcs_=0;
            pos_shot_mcs=0;pos_shot_grid_mcs=0;
            pos_bat_mcs=0;pos_bat_grid_mcs=0;

            allocate(offset_su(NumChannels,NumShots_mcs))
            allocate(pos_rec(NumChannels,NumShots_mcs))
            allocate(pos_rec_grid(NumChannels,NumShots_mcs))
            offset_su=0;pos_rec=0;pos_rec_grid=0;

!!		Geometry en MCS
            allocate(nxSou_MCS(NumSou_MCS),nySou_MCS(NumSou_MCS))
            allocate(nxRec_MCS(NumRec_MCS,NumSou_MCS),nyRec_MCS(NumRec_MCS,NumSou_MCS))
            allocate(MCS_raw_Source(nt_mcs,NumSou_MCS),MCS_raw_Data(nt_mcs,NumRec_MCS))
            nxSou_MCS=0;nySou_MCS=0;nxRec_MCS=0;nyRec_MCS=0;
            MCS_raw_Source=0.;MCS_raw_Data=0.;

        endif

end subroutine allocate_data_arrays

subroutine allocate_model_arrays_was(rank)
use mod_parfile

    implicit none
	integer	:: rank

    allocate(bat_model_was(nxmodel),bat_model_grid_was(nxmodel))
    bat_model_was=0.;bat_model_grid_was=0;

end subroutine allocate_model_arrays_was

subroutine allocate_model_arrays_mcs(rank)
use mod_parfile

    implicit none
	integer	:: rank

    allocate(bat_model_mcs(nxmodel),bat_model_grid_mcs(nxmodel))
    bat_model_mcs=0.;bat_model_grid_mcs=0;

end subroutine allocate_model_arrays_mcs

subroutine allocate_model_arrays(rank)
use mod_parfile

    implicit none
	integer	:: rank

	allocate(Model_ini(nymodel,nxmodel),lambda(nymodel,nxmodel),rho(nymodel,nxmodel))
        Model_ini=0.;lambda=0.;rho=0.;

end subroutine allocate_model_arrays

subroutine deallocate_arrays()
use mod_parfile

        if(inv_WAS.ne.0)    then
            deallocate(shotID_nav_was)
            deallocate(shotID_su_was,shotID_was_)
            deallocate(pos_shot_was,pos_shot_grid_was)
            deallocate(pos_bat_was,pos_bat_grid_was)
            deallocate(nxSou_WAS,nySou_WAS,nxRec_WAS,nyRec_WAS)
            deallocate(bat_model_was,bat_model_grid_was)
	    deallocate(sizeof_was)
        endif

        if(inv_MCS.ne.0)    then
            deallocate(shotID_nav_mcs)
            deallocate(shotID_su_mcs,shotID_mcs_)
            deallocate(pos_shot_mcs,pos_shot_grid_mcs)
            deallocate(offset_su)
            deallocate(pos_bat_mcs,pos_bat_grid_mcs)
            deallocate(nxSou_MCS,nySou_MCS,nxRec_MCS,nyRec_MCS)
            deallocate(MCS_raw_Source,MCS_raw_Data)
            deallocate(bat_model_mcs,bat_model_grid_mcs)
        endif

        deallocate(Model_ini,lambda,rho,Freq)

end subroutine deallocate_arrays


end module mod_data_arrays
