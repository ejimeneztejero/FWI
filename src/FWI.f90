subroutine FWI(iter,f1,f2)

use mod_parfile
use mod_data_arrays

implicit none
include 'mpif.h'

integer :: numtasks,rank
integer :: ierr,status(MPI_STATUS_SIZE),TAG

integer :: iFWI
integer :: nsamples
integer :: i,j,k,iter
integer :: step
integer :: loopind,looptotal,looptest,ind
integer :: ntold,irep
integer :: nt,nt_data

real :: tfin
real :: f1,f2
real :: dt
real :: Misfit
real :: Misfit_test0,Misfit_test1,Misfit_test!,Misfit_test_
real :: nut,nut0,tmp
real :: Beta,EndStep
real :: S0(nt_mcs)

real, allocatable :: Vpmodel(:,:)
real, allocatable :: var_orig(:,:)
real, allocatable :: Grad(:,:)
real, allocatable :: Search(:,:)
real, allocatable :: OldGrad(:,:)

character(len=1000) :: file_name_shot,file_name
character(len=50) :: Str_iter,Str_iter2,Str_f1,Str_f2,name_iter

call MPI_COMM_SIZE(MPI_COMM_WORLD,numtasks,ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)

TAG=0

allocate(Vpmodel(nymodel,nxmodel))
allocate(var_orig(nymodel,nxmodel))
allocate(Grad(nymodel,nxmodel))
allocate(Search(nymodel,nxmodel))
allocate(OldGrad(nymodel,nxmodel))

iFWI=1
EndStep=0.
OldGrad=0
var_orig=lambda

nsamples=NumSS_was+NumSS_mcs
if(rank.eq.0)	then
if(numtasks.eq.1.and.nsamples.gt.1) then
	write(*,*)
	write(*,*)'NO PARALELIZATION'
endif

if(numtasks.gt.nsamples) then
	write(*,*)
	write(*,*)'WARNING: DO NOT WASTE RESOURCES, BETTER CHOOSE YOUR CPU NUMBER LESS OR EQUAL THAN: ',nsamples
endif
endif


call MPI_barrier(MPI_COMM_WORLD,ierr)

22	do while (iFWI .le. nFWI)


if(rank.eq.0)   then
        write(*,*)
        write(*,*)"`'-.,_,.-'``'-.,_,.='``'-.,_,.-'``'-.,_,.='````'-.,_,.-'``'-.,_,.='``'-.             "
        write(*,*)"`'-.,_,.-'``'-.,_,.='``'-.,_,.-'``'-.,_,.='````'-.,_,.-'``'-.,_,.='``'-.             "
        write(*,*)"`'-.,_,.-'``'-.,_,.='``'-.,_,.-'``'-.,_,.='````'-.,_,.-'``'-.,_,.='``'-.             "
        write(*,*)
        write(*,*) ' ....... iteration number .......'
        write(*,*)		iFWI	
        write(*,*)
        write(*,*)"`'-.,_,.-'``'-.,_,.='``'-.,_,.-'``'-.,_,.='````'-.,_,.-'``'-.,_,.='``'-.             "
        write(*,*)"`'-.,_,.-'``'-.,_,.='``'-.,_,.-'``'-.,_,.='````'-.,_,.-'``'-.,_,.='``'-.             "
        write(*,*)"`'-.,_,.-'``'-.,_,.='``'-.,_,.-'``'-.,_,.='````'-.,_,.-'``'-.,_,.='``'-.             "
        write(*,*)
endif

iFWI=iFWI+1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!     MAIN LOOP Paralelization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(rank.eq.0)write(*,*)"-------------------------------"
if(rank.eq.0)write(*,*) 'START Adjoint Method'
if(rank.eq.0)write(*,*)"-------------------------------"

Grad=0.
Misfit=0.
step=1

if(inv_source.eq.1.and.iFWI-1.eq.1)     then
        Vpmodel=sqrt(lambda/1000.)
        call source_estimation(iFWI-1,f1,f2,nxmodel,nymodel,nt_mcs,NumRec_MCS,dmodel,dt_mcs,&
        nxSou_mcs(1),nySou_mcs(1),nxRec_mcs(:,1),nyRec_mcs(:,1),Vpmodel,S0)     !ojo
endif
        
if(inv_source.eq.0)     then
        call Ricker_mcs(rank,1,S0)
endif

!!        call time_filter(S0,nt_mcs,dt_mcs,1,typef,f1,f2)
!        do j=1,nt_mcs
!	        write(40+iFWI,*) dt_mcs*(j-1),S0(j)    
!        enddo
!        close(40+iFWI)

call calculate_misfit_gradient(iter,step,f1,f2,dt,nt,tfin,S0,Misfit,Grad)

if(rank.eq.0)write(*,*)
if(rank.eq.0)write(*,*)"*************************"
if(rank.eq.0)write(*,*)"TOTAL MISFIT:", Misfit
if(rank.eq.0)write(*,*)"*************************"
if(rank.eq.0)write(*,*)

if(rank.eq.0)write(*,*)"-------------------------------"
if(rank.eq.0)write(*,*) 'END Adjoint Method'
if(rank.eq.0)write(*,*)"-------------------------------"
if(rank.eq.0)write(*,*)

!!!!!---------------------
!!!!!------ End of Adjoint problem and Model gradient computation
!!!!!---------------------

call MPI_barrier(MPI_COMM_WORLD,ierr)
	

if(rank.eq.0)	then
write(*,*)
write(*,*)"-------------------------------------------------------------"
write(*,*)"---- Determine Search direction as a function of the gradient"
write(*,*)"-------------------------------------------------------------"
write(*,*)
endif
		    				
if (method .eq. 1) then !NLCG

	if (iFWI-1 .eq. 1) then
		Search=-Grad
	else
		call Beta_calculous(Grad,OldGrad,nxmodel,nymodel,Beta) !!!!!------ Non-linear Conjugate Gradient (Polak-Ribiere) 
 		if (Beta .gt. 0.) then
		 	Search=-Grad+Beta*Search
		 else
			Search=-Grad
			!!!write(*,*) '..... Gradient Reinitialization .....'
		endif
	endif
		
else !!!!!----- Steepest descent
		Search=-Grad
endif
				
OldGrad=Grad	!former iteration

!!!!!-----------------------------------------------------------------
!!!!!------ End determine Search direction as function of the gradient
!!!!!-----------------------------------------------------------------
	
		    
	!!!!!------  Obtain a test step that achieves the percentual minimum variation that we allow
	!!!!!		 to the model      
	            
        if (EndStep.eq.0.) then
		nut=1./100./maxval(abs(Search/log(var_orig)))
	else
		nut=EndStep
	endif

!!	nut0=nut

	if(rank.eq.0)	then
		write(*,*)
		write(*,*) 'Start Optimization with nut= ',nut
		write(*,*) maxval(abs(Search/log(var_orig))),maxval(abs(Search)),maxval(abs(log(var_orig)))
		write(*,*)		
	endif	


if(rank.eq.0)	then
write(*,*)
write(*,*)"---------------------------------------------"
write(*,*)"------ Line Search Optimisation -------------"
write(*,*)"---------------------------------------------"
write(*,*)
endif

	Misfit_test0=Misfit
	Misfit_test1=Misfit
	loopind=0
	looptest=0
	looptotal=0

        call MPI_barrier(MPI_COMM_WORLD,ierr)

20	do while (looptest .eq. 0)
			
		loopind=loopind+1
		looptotal=looptotal+1
			
		!!!!!------  Update test model
		   
		lambda=exp(Search*nut)*var_orig
	
		!!!!!----- keep in memory the two test misfit values
			
		if (loopind .ge. 2) then
			
			Misfit_test0=Misfit_test1
			Misfit_test1=Misfit_test
				
		endif
			
		!!!!!------ New time sampling calculation (CFL conditions)

		ntold=nt
		dt=dmodel/sqrt(maxval(lambda/rho))/2.
		nt=1+ceiling(tfin/dt)

		if(rank.eq.0)write(*,*) 'FWI: nut = ',nut,'nt: ',nt, 'tfin: ',tfin,' dt: ',dt

		Misfit_test=0.
		step=2
		call calculate_misfit_gradient(iter,step,f1,f2,dt,nt,tfin,S0,Misfit_test,Grad)

	        call MPI_barrier(MPI_COMM_WORLD,ierr)

	!!!!!-------  End Stop frequency inversion condition -----!!!!!

		if(rank.eq.0)write(*,*)
		if(rank.eq.0)write(*,*) 'Total Misfit Test=',Misfit_test
		if(rank.eq.0)write(*,*)

	if (looptotal.gt.6) then
						
		EndStep=0.
		Search=0.
		irep=irep+1

		if(rank.eq.0) write(*,*)"repite num ",irep

		if(irep.ge.2) then
			if(rank.eq.0)write(*,*)'freq estres'
			iFWI=nFWI+1	!para evitar el loop interminable de iteraciones
		else 
			if(rank.eq.0)write(*,*)'REPITE ITERACION'
			iFWI=iFWI-1
		endif
		goto 22
	endif
				
	if ((Misfit_test-Misfit) .ge. 0.) then
		if (loopind .eq. 1) then
			if(rank.eq.0)write(*,*)'nut too big, it gets divided by 2'
			if(rank.eq.0)write(*,*)'Start loop again'
			nut=nut/2
			loopind=0
			goto 20	
		endif		
	endif
					
	if ((Misfit_test-Misfit_test1) .le. 0.) then
		nut=nut*2
	else
		if(rank.eq.0)then
			write(*,*)'out of loop'
		endif
		looptest=1	
	endif
	
	enddo!20 do while

if(rank.eq.0)	then
write(*,*)
write(*,*)"---------------------------------------------"
write(*,*)"------ END Optimisation -------------"
write(*,*)"---------------------------------------------"
write(*,*)
endif


if(rank.eq.0)	then
write(*,*)
write(*,*)"---------------------------------------------"
write(*,*)"------ Solve quadratic interpolation---------"
write(*,*)"---------------------------------------------"
write(*,*)
endif

	if (loopind .eq. 2) then    
		call polyfit(nut,Misfit,Misfit_test1,Misfit_test,loopind,EndStep)
	else
  		call polyfit(nut,Misfit_test0,Misfit_test1,Misfit_test,loopind,EndStep)
	endif
		
!	if(rank.eq.0)write(*,*)		
!	if(rank.eq.0)write(*,*) 'EndStep =',EndStep

	!!!---------------------------------------------
	!!!------ End Line Search Optimisation algorithm
	!!!---------------------------------------------
		
	!!!!!------ Update Model
								
	lambda=exp(Search*EndStep)*var_orig
	var_orig=lambda
        Vpmodel=sqrt(lambda/1000.)

	!!!!!----- Save iteration result
	call MPI_barrier(MPI_COMM_WORLD,ierr)
        if(rank.eq.0)   then

		write(*,*)
		write(*,*)"---------------------------------------------"
		write(*,*)"     UPDATE MODELS IN: ",trim(folder_output)
		write(*,*)"---------------------------------------------"
		write(*,*)
	
		!!!!!----- Write frequencies and relative and global iteration
		write(unit_iter,*) f1,f2,iFWI-1,iter

		!!!!!----- Write Misfit at each (global) iteration
		write(unit_misfit,*) iter,Misfit

		write(Str_iter,*) iter
		write(Str_iter2,*) iFWI-1
		write(Str_f1, "(F3.1)") f1
		write(Str_f2, "(F3.1)") f2

		call file_iter(iter,name_iter)

		!!!!!----- Write new Grad model at each (global) iteration
		file_name=trim(folder_GRAD) // trim("grad") // trim(name_iter) // trim(adjustl(Str_iter))
		call write_model(nxmodel,nymodel,file_name,Grad)

		!!!!!----- Write new Vp model at each (global) iteration
		file_name=trim(folder_VEL) // trim("vp") // trim(name_iter) // trim(adjustl(Str_iter))

		call write_model(nxmodel,nymodel,file_name,Vpmodel)

	endif
	call MPI_barrier(MPI_COMM_WORLD,ierr)

	iter=iter+1

enddo   !iFWI

deallocate(Vpmodel)
deallocate(var_orig)
deallocate(Grad)
deallocate(Search)
deallocate(OldGrad)

!!!!-----------------------
!!!!------ End of Inversion
!!!!-----------------------	

end subroutine FWI
