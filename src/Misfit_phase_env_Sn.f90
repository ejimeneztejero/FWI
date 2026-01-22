subroutine Misfit_phase_env_Sn(Data_synth,Data_real,nt,dt,NumRec,Misfit,AdjSource)

use mod_parfile, only: misfit_type

	implicit none

	integer :: j,k,l,ntt,nttS,nttR,nt,ll,sec1,NumRec

	real :: dt,Misfit,sgn
	real :: Data_synth(nt,NumRec)
	real :: Data_real(nt,NumRec)
	real :: phase_synth(nt,NumRec),phase_real(nt,NumRec)
	real :: env_synth(nt,NumRec),env_real(nt,NumRec)
	real :: Sn_synth(nt,NumRec),Sn_real(nt,NumRec)
	real :: AdjSource(nt,NumRec)

	real, allocatable :: S(:),T(:),DS(:)
	real, allocatable :: L1(:),Sh(:),S_ad(:)
	real, allocatable :: envS(:),envT(:)
	real, allocatable :: phaseS(:),phaseT(:)
	real, allocatable :: SnS(:),SnT(:)
	real, allocatable :: DQS(:),DQT(:)
	real, allocatable :: DQ_synth(:,:),DQ_real(:,:)

	real, parameter :: delta=0.1

!!!	phase, evnelope, Sn

		nttR=0;nttS=0;

		allocate(S(nt),T(nt),DS(nt))
		allocate(L1(nt),Sh(nt),S_ad(nt))
		allocate(envS(nt),envT(nt))
		allocate(phaseS(nt),phaseT(nt))
		allocate(SnS(nt),SnT(nt))

		allocate(DQS(nt),DQT(nt))
		allocate(DQ_synth(nt,NumRec),DQ_real(nt,NumRec))

		S=0;T=0;DS=0;
		envS=0;envT=0;
		phaseS=0;phaseT=0;
		SnS=0;SnT=0;

		DQS=0;DQT=0;
		DQ_synth=0;DQ_real=0;

		do j=1,NumRec
		
		        do k=1,nt
				S(k)=Data_synth(k,j)
				T(k)=Data_real(k,j)
			enddo

			S=S/maxval(abs(S))
			T=T/maxval(abs(T))

!!			calculo sintético
			call first_arrivals(S,nt,nttS)
			call sub1_cycle(1,j,nt,dt,nttS,S,envS,phaseS,SnS,DQS)
               
!!			calculo real
			call first_arrivals(T,nt,nttR)
			call sub1_cycle(2,j,nt,dt,nttR,T,envT,phaseT,SnT,DQT)

!!!			diferencia de fase
			do k=1,nt

				phase_synth(k,j)=phaseS(k)
				phase_real(k,j)=phaseT(k)

				Sn_synth(k,j)=SnS(k)
				Sn_real(k,j)=SnT(k)

				env_synth(k,j)=envS(k)
				env_real(k,j)=envT(k)

				DQ_synth(k,j)=DQS(k)
				DQ_real(k,j)=DQT(k)

			enddo

		enddo	!!j
	
		if(misfit_type.eq.2)	then	!!phase triangular

			do j=1,NumRec
			        do k=1,nt

					L1(k)=phase_synth(k,j)-phase_real(k,j)
	        	                Misfit=Misfit+0.5*dt*L1(k)**2

		                        Sh(k)=env_synth(k,j)*DQ_synth(k,j)
        		                !if(Sh(k).ne.0)sgn=Sh(k)/abs(Sh(k))
	                   	        S_ad(k)=L1(k)*abs(Sh(k))/(env_synth(k,j)+delta)**2.

				enddo
			enddo

		endif

		if(misfit_type.eq.3)	then	!!env

			do j=1,NumRec
			        do k=1,nt

					L1(k)=env_synth(k,j)-env_real(k,j)
	        	                Misfit=Misfit+0.5*dt*L1(k)**2
					S_ad(k)=-L1(k)*Data_synth(k,j)/(env_synth(k,j)+delta)

				enddo
			enddo

		endif

		if(misfit_type.eq.4)	then	!!señal normalizada

			do j=1,NumRec
			        do k=1,nt

					L1(k)=Sn_synth(k,j)-Sn_real(k,j)
	        	                Misfit=Misfit+0.5*dt*L1(k)**2
					S_ad(k)=-L1(k)*DQ_synth(k,j)*DQ_synth(k,j)/(env_synth(k,j)+delta)

				enddo
			enddo

		endif


		do j=1,NumRec
		        do k=1,nt
              			AdjSource(k,j)=-S_ad(nt-k+1)
       			enddo
		enddo


	deallocate(S,T,DS)
	deallocate(L1,Sh,S_ad)
	deallocate(envS,envT)
	deallocate(phaseS,phaseT)
	deallocate(SnS,SnT)
	deallocate(DQS,DQT,DQ_synth,DQ_real)

end subroutine Misfit_phase_env_Sn
