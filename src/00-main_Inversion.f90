program Inversion_Main

use mod_parfile
use mod_data_arrays
implicit none
include 'mpif.h'

integer :: numtasks, rank, ierr, status(MPI_STATUS_SIZE)
integer :: i,j,k
integer :: itimes, ntimes
integer :: ifreq,iter
real :: f1,f2
character(len=500) :: Strtmp,file_name
call MPI_INIT(ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,numtasks,ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)

!!!!!-------------------
!!!!! Getting input data parameters
!!!!!-------------------

if(rank.eq.0)    then
        write(*,*)
        write(*,*)'*******************************'
        write(*,*)'READING PARAMETERS FROM PARFILE'
        write(*,*)'*******************************'

endif

call read_parfile(rank,numtasks)

if(rank.eq.0)    then
        write(*,*)
        write(*,*)'*** PARAMETER CHECK OK ***'
endif

call allocate_data_arrays()

if(rank.eq.0)    then
        write(*,*)
        write(*,*)'*******************************'
        write(*,*)'READING AND CHECKING INPUT FILES '
        write(*,*)'*******************************'

endif

if(inv_WAS.ne.0) call get_was_data(rank)
if(inv_MCS.ne.0) call get_mcs_data(rank)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call allocate_model_arrays(rank)    !!from was and/or mcs

if(rank.eq.0)   then
        write(*,*)
        write(*,*)'*********************************'
        write(*,*)'GET INITIAL MODEL '
        write(*,*)'*********************************'
endif

call get_initial_models(rank)   

!!!!!-------------------
!!!!! Starting inversion
!!!!!-------------------

if(rank.eq.0)    then
        write(*,*)
        write(*,*)'*******************************'
        write(*,*)'STARTING INVERSION'
        write(*,*)'*******************************'
endif

	iter=1
	do ifreq=2,NumFreq		!siempre un 2?, repasar, mejor que empiece en 1

		call get_f1f2(ifreq,f1,f2)

		if(rank.eq.0)	then	
			write(*,*) 
			if(typef.eq.1)write(*,*) 'Lowpass at FREQ = ',f2
			if(typef.eq.3)write(*,*) 'Bandpass between FREQ = ',f1,' and ',f2
		endif

		call FWI(iter,f1,f2)

	enddo

    if(inv_WAS.ne.0) call close_was_files()
    if(inv_MCS.ne.0) call close_mcs_files()

	call deallocate_arrays()

	write(*,*)'end rank',rank

	call MPI_FINALIZE(ierr)

	!!!!!-----------------
	!!!!! Ending inversion
	!!!!!-----------------	

end program Inversion_Main
