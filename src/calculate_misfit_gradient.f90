subroutine calculate_misfit_gradient(iter,step,f1,f2,dt,nt,tfin,S0,Misfit,Grad)
use mod_parfile
use mod_data_arrays

implicit none
include 'mpif.h'

integer :: numtasks,rank
integer :: ierr,status(MPI_STATUS_SIZE),TAG

integer :: iter,step
integer :: itimes,icount,nsamples,ntimes
integer :: i,j,k,itype,kk
integer :: isou,iSS
integer :: NumRec,NumSou
integer :: store_snap,nt_Store
integer :: nt,nt_data

integer :: nSSwas
integer :: NumRec_SS,NumSou_SS
integer :: SourceNum

real :: S0(nt_mcs)
real :: tfin,dt,dt_data
real :: f1,f2,temp
real :: Misfit,Misfit_,Misfit_thread
real :: Grad(nymodel,nxmodel)

integer, allocatable :: nxSou(:),nxRec(:,:)
integer, allocatable :: nySou(:),nyRec(:,:)
integer, allocatable :: nxSou_SS(:),nxRec_SS(:)
integer, allocatable :: nySou_SS(:),nyRec_SS(:)
real, allocatable :: Field_Data(:,:),Field_Source(:,:)
real, allocatable :: Fields(:,:),Data_target(:,:),Data_synth(:,:)
real, allocatable :: Data_source(:,:),AdjSource(:,:)
real, allocatable :: Coding(:)
real, allocatable :: Vpmodel(:,:)
real, allocatable :: Grad_(:,:)
real, allocatable :: Grad_thread(:,:)

character(len=1000) :: file_name_shot,file_name
character(len=50) :: Str_iter,Str_iter2,Str_f1,Str_f2

call MPI_COMM_SIZE(MPI_COMM_WORLD,numtasks,ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
TAG=0

itype=0
store_snap=10

allocate(Vpmodel(nymodel,nxmodel))
allocate(Grad_(nymodel,nxmodel))
allocate(Grad_thread(nymodel,nxmodel))

!!  SS
nSSwas=ceiling(1.*NumOBS/NumSS_was)     !!numero de shots para este SS
allocate(Coding(nSSwas))
Coding=1
!call coding_fun(nSSwas,Coding)

nsamples=NumSS_was+NumSS_mcs

ntimes=nsamples
if(numtasks.gt.1)    then
        if(nsamples.gt.numtasks)        then
                ntimes=ceiling(1.*nsamples/numtasks)
        endif
        if(nsamples.le.numtasks)        then
               ntimes=1
        endif
endif

do itimes=1,ntimes

icount=(itimes-1)*numtasks+rank+1        !!aquí sucede la paralelización

if(rank.eq.0.and.ntimes.gt.1)   then
        write(*,*)
        write(*,*)"oooooooooooooooooooooooooooooooooooooooooooooooooooooo"
        write(*,*)'ROUND ',itimes, 'OUT OF', ntimes
        write(*,*)'icount ', icount
        write(*,*)"oooooooooooooooooooooooooooooooooooooooooooooooooooooo"
        write(*,*)
endif

if(rank.eq.0)    then
        write(*,*)
        write(*,*)'READING FIELD DATA ...'
endif

    Misfit_=0.
    Grad_=0.
    nt=1

    itype=0
    iSS=1
    nt_data=1
    dt_data=1
    tfin=1
    NumRec=1
    NumSou=1
    NumSou_SS=1
    NumRec_SS=1

if(icount.le.nsamples)  then

if(inv_WAS.ne.0.and.icount.le.NumSS_was)    then ! WAS data

    itype=1
    iSS=icount
    nt_data=nt_was
    dt_data=dt_was
    tfin=tfin_was
    NumRec=NumRec_WAS
    NumSou=NumSou_WAS
!!  SS para WAS
    NumSou_SS=ceiling(1.*NumOBS/NumSS_was)     !!numero de shots para este SS
    NumRec_SS=NumRec

endif

if(inv_MCS.ne.0.and.icount.gt.NumSS_was)    then ! MCS data

    itype=2
    iSS=icount-NumSS_was
    nt_data=nt_mcs
    dt_data=dt_mcs
    tfin=tfin_mcs
    NumRec=NumRec_MCS
    NumSou=NumSou_MCS
    NumSou_SS=1
    NumRec_SS=NumRec

endif

endif!nsamples

allocate(nxSou(NumSou),nySou(NumSou))
allocate(nxRec(NumRec,NumSou),nyRec(NumRec,NumSou))
allocate(nxSou_SS(NumSou_SS),nySou_SS(NumSou_SS))
allocate(nxRec_SS(NumRec_SS),nyRec_SS(NumRec_SS))
allocate(Field_Data(nt_data,NumRec_SS))
allocate(Field_Source(nt_data,NumSou_SS))
nxSou=0;nySou=0;
nxRec=0;nyRec=0;
nxSou_SS=0;nySou_SS=0;
nxRec_SS=0;nyRec_SS=0;
Field_Data=0;

if(icount.le.nsamples)  then

if(inv_WAS.ne.0.and.icount.le.NumSS_was)    then ! WAS data

    call get_obs_data(rank,iSS,NumSou_SS,NumRec_SS,Coding,Field_data)    !! raw_data SU-WAS files
    !call Ricker_was(rank,iSS,NumSou_SS,Coding,S0)
    Field_Source(1:nt_data,1)=S0(1:nt_data);

    nxSou=nxSou_WAS
    nySou=nySou_WAS
    nxRec=nxRec_WAS
    nyRec=nyRec_WAS

    do isou=1,NumSou_SS
        nxSou_SS(isou)=nxSou(SourceNum(iSS,isou,NumSS_was))
        nySou_SS(isou)=nySou(SourceNum(iSS,isou,NumSS_was))
        nxRec_SS=nxRec(:,isou);nyRec_SS=nyRec(:,isou);
    enddo

endif

if(inv_MCS.ne.0.and.icount.gt.NumSS_was)        then ! MCS data

    call get_shot_data(iSS,Field_Data)  !! raw_data SU-MCS files
    Field_Source(1:nt_data,1)=S0(1:nt_data);

    nxSou=nxSou_MCS
    nySou=nySou_MCS
    nxRec=nxRec_MCS
    nyRec=nyRec_MCS

    nxSou_SS=nxSou(iSS)
    nySou_SS=nySou(iSS)
    nxRec_SS=nxRec(:,iSS)
    nyRec_SS=nyRec(:,iSS)

endif


!!!!!-----------------------------------
!!!!!------ Frequency filter of source and data
!!!!!-----------------------------------

  call time_filter(Field_Data,nt_data,dt_data,NumRec_SS,typef,f1,f2)
  call time_filter(Field_Source,nt_data,dt_data,NumSou_SS,typef,f1,f2)

!!!!!------ New time sampling

dt=dmodel/sqrt(maxval(lambda/rho))/2.
nt=1+ceiling(tfin/dt)

temp=nt/store_snap
nt_Store=ceiling(temp)
nt=nt_Store*store_snap

if(rank.eq.0)write(*,*)
if(rank.eq.0)write(*,*)'calc1 - tfin,nt,dt: ',tfin,nt,dt

endif !nsamples

if(icount.eq.1)	then
	 write(111,*)"icount 1 ",rank,NumSou_SS,NumRec_SS,dmodel
endif
if(icount.eq.2)then
	 write(112,*)"icount 2 " 
endif

write(*,*)"icount ",icount 


!!!!!------ Allocate Synthetic trace matrix, field snapshot and time vector 

allocate(Data_source(nt,NumSou_SS))
allocate(Data_target(nt,NumRec_SS))
allocate(Data_synth(nt,NumRec_SS))
allocate(AdjSource(nt,NumRec_SS))
allocate(Fields(nymodel,nxmodel*(nt_Store+1)))
Data_Source=0
Data_target=0
Data_synth=0
AdjSource=0
Fields=0

! write(110+icount,*)rank,NumSou_SS,NumRec_SS,nxSou_SS,nySou_SS,nxRec_SS,nyRec_SS,&
! dmodel,dt,nymodel,nxmodel,store_snap,nt_Store,nt

if(icount.le.nsamples)  then

        if(rank.eq.0)write(*,*)"3 *KK",Field_Source(10,1)

        !!!!!------ Interpolation and Filtering Seismic Data with filter_function1D     
        call interpolation_akima(Field_Data,dt_data,nt_data,Data_target,NumRec_SS,dt,nt)

        !!!!!------ Filtering source with filter_function1D
        call interpolation_akima(Field_Source,dt_data,nt_data,Data_source,NumSou_SS,dt,nt)

        !!!!!---------------------          
        !!!!!------ Forward solver
        !!!!!---------------------

        if(rank.eq.0)write(*,*)"4 *KK"

        Vpmodel=sqrt(lambda/1000.)

        if(rank.eq.0)write(*,*)"5 *KK"

        if(rank.eq.0)write(*,*)
        if(rank.eq.0)write(*,*)"00 Forward solver ..."

        call solver_acoustic_SS(rank,1,NumSou_SS,NumRec_SS,nxSou_SS,nySou_SS,nxRec_SS,nyRec_SS,&
        Data_source,AdjSource,dmodel,dt,nymodel,nxmodel,Vpmodel,&
        store_snap,nt_Store,nt,Grad_thread,Fields,Data_synth)

        if(rank.eq.0)write(*,*)"00 After forward solver ..."

	if(step.eq.1)	then
        !!!!!----- Write new synthetic data
        if(inv_WAS.ne.0.and.icount.le.NumSS_was)    then
                if(rank.eq.0)write(*,*)
		if(rank.eq.0)write(*,*)"-----------------------------------------------------"
                if(rank.eq.0)write(*,*)'... writting OBS synth data in: ',trim(folder_output)
		if(rank.eq.0)write(*,*)"-----------------------------------------------------"
                call write_new_data(rank,itype,iter,nt,dt,iSS,NumSou_SS,NumRec_SS,Data_synth)   !! mejorar cabeceras
        endif
        if(inv_MCS.ne.0.and.icount.gt.NumSS_was.and.icount.le.nsamples)    then
                if(rank.eq.0)write(*,*)
		if(rank.eq.0)write(*,*)"------------------------------------------------------------"
                if(rank.eq.0)write(*,*)'... writting shotgather synth data in: ',trim(folder_output)
		if(rank.eq.0)write(*,*)"------------------------------------------------------------"
                call write_new_data(rank,itype,iter,nt,dt,iSS,NumSou_SS,NumRec_SS,Data_synth)   !! mejorar cabeceras
        endif
	endif

!!!!!----------------------
!!!!!------ Adjoint problem
!!!!!----------------------

        if(rank.eq.0) write(*,*)
        if(rank.eq.0) write(*,*)"Misfit ..."

        Misfit_thread=0.
        call Misfit_function(rank,f1,Data_synth,Data_target,NumRec_SS,nt,dt,Misfit_thread,AdjSource)

!!        do j=1,nt
!!               write(77,*) dt*(j-1),Data_target(j,1),Data_synth(j,1)
!!        enddo

endif

call MPI_REDUCE(Misfit_thread,Misfit_,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(Misfit_,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)

if(rank.eq.0.and.ntimes.gt.1)write(*,*) 'Partial Misfit Test =',Misfit_

Misfit=Misfit+Misfit_

if(icount.le.nsamples.and.step.eq.1)  then

        if(rank.eq.0) write(*,*)
        if(rank.eq.0) write(*,*)"Propagation of residuals ..."
        if(rank.eq.0) write(*,*)"and Gradient calculation ..."

        Vpmodel=sqrt(lambda/1000.)
        Grad_thread=0.

        call solver_acoustic_SS(rank,2,NumSou_SS,NumRec_SS,nxSou_SS,nySou_SS,nxRec_SS,nyRec_SS,&
        Data_source,AdjSource,dmodel,dt,nymodel,nxmodel,Vpmodel,&
        store_snap,nt_Store,nt,Grad_thread,Fields,Data_synth)

!        do j=1,nymodel
!               write(88,*) j,Grad_thread(j,nxmodel/2)
!        enddo

endif


!!!!!-------------------------      
!!!!!------ End Forward solver
!!!!!-------------------------

if(step.eq.1)  then

!!!!!----------------------
!!!!!------ Adjoint problem
!!!!!----------------------

!!! Grad calculado por cada thread, se suma entre todos en el rank 0
!!! y el resultado se envia a todos los rank

	call MPI_REDUCE(Grad_thread,Grad_,nxmodel*nymodel,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(Grad_,nxmodel*nymodel,MPI_REAL,0,MPI_COMM_WORLD,ierr)

	Grad_=Grad_/(numtasks)
	Grad=Grad+Grad_

	do i=1,nxmodel
	        if(inv_MCS.ne.0)kk=bat_model_grid_mcs(i)
	        if(inv_MCS.eq.0)kk=bat_model_grid_was(i)
	        do j=1,kk
	                Grad(j,i)=0.    !!water is not inverted
	        enddo
	enddo

endif

deallocate(nxSou,nySou,nxRec,nyRec)
deallocate(nxSou_SS,nySou_SS,nxRec_SS,nyRec_SS)
deallocate(Field_Data,Field_Source)
deallocate(AdjSource,Data_target,Data_source,Data_synth,Fields)

enddo   !itimes 

deallocate(Vpmodel)
deallocate(Grad_thread)
deallocate(Coding)

end subroutine calculate_misfit_gradient
