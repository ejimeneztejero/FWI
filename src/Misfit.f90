subroutine Misfit_function(rank,freq,Data_Synth,Data_trace,NumRec,nt,dt,Misfit,AdjSource)

implicit none

	integer :: TID,option,iref,ioff,NumRec,nt,i,j,k,inorm,rank
	integer :: offset1,offset2,indfile,grad
	real ::  dt,Misfit,sigma_av,maxT,maxS,freq
	real :: Data_trace(nt,NumRec),Data_Synth(nt,NumRec),AdjSource(nt,NumRec)
	integer :: EAj(nt),FAj(nt)
	real, allocatable :: time(:)
        real, allocatable :: Data_S(:,:),Data_T(:,:),DData_S(:,:)
        real, allocatable :: Data_(:),DData_(:)
	logical	:: isnotCero1,isnotCero2

	allocate(time(nt))
	allocate(Data_S(nt,NumRec),Data_T(nt,NumRec),DData_S(nt,NumRec))
	allocate(Data_(nt),DData_(nt))
	Data_S=0;Data_T=0.;

	do k=1,nt
		time(k)=(k-1)*dt
	enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!! TRAZAS NORMALIZADAS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	inorm=1	!! trazas normalizadas
!	inorm=0	!! trazas SIN normalizar

	do j=1,NumRec

		call is_not_cero(nt,Data_trace(:,j),isnotCero1)
		call is_not_cero(nt,Data_Synth(:,j),isnotCero2)
	
		if(isnotCero1)	then

			maxS=1
			maxT=1

			if(inorm.eq.1)	then

				maxS=maxval(abs(Data_Synth(:,j)))
				maxT=maxval(abs(Data_trace(:,j)))

				if(maxS.ne.0.and.maxT.ne.0)	then	
					Data_S(:,j)=Data_Synth(:,j)/maxS
					Data_T(:,j)=Data_trace(:,j)/maxT

				endif
			endif

!			if(option.eq.3)	then
!				Data_(:)=Data_S(:,j)
!				call derivative(nt,time,Data_,DData_) 
!				DData_S(:,j)=DData_(:)
!			endif

		endif

	enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	Misfit=0.
	AdjSource=0.

	call Misfit_option(freq,Data_S,Data_T,nt,dt,NumRec,Misfit,AdjSource)

	deallocate(time)
	deallocate(Data_S,Data_T,DData_S)
	deallocate(Data_,DData_)

return
end


subroutine Misfit_option(freq,Data_synth,Data_real,nt,dt,NumRec,Misfit,AdjSource)


use mod_parfile, only: misfit_type

        implicit none

        integer :: j,k,l,nt,NumRec

        real :: dt,time_period,Misfit,freq
        real :: Data_synth(nt,NumRec)
        real :: Data_real(nt,NumRec)
        real :: AdjSource(nt,NumRec)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!     L2
	if(misfit_type.eq.1) then

                call Misfit_L2(Data_synth,Data_real,NumRec,nt,dt,Misfit,AdjSource)

        endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!     phase, evnelope, Sn
        if(misfit_type.eq.2.or.misfit_type.eq.3.or.misfit_type.eq.4)    then

                call  Misfit_phase_env_Sn(Data_synth,Data_real,nt,dt,NumRec,Misfit,AdjSource)

        endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!      CCTT

	if(misfit_type.eq.5)    then
                time_period=1./freq
                call Misfit_CCTT(time_period,Data_synth,Data_real,nt,dt,NumRec,Misfit,AdjSource)
        endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end subroutine Misfit_option
