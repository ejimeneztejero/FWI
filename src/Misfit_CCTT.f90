subroutine Misfit_CCTT(time_period,Data_Synth,Data_trace,nt,dt,NumRec,Misfit,AdjSource)
implicit none
	
	integer :: nt,NumRec,i,j,k,kfin,delta,indfile,jj
	integer :: ntau,itau,m,imax

	real :: dt,maxS,maxT,Tp,Iu22,tau0,u22,t1,t2,time_period,Misfit
	real :: vtau(NumRec)	
	real :: Data_trace(nt,NumRec),Data_Synth(nt,NumRec)
        real :: AdjSource(nt,NumRec)

        real, allocatable :: Data_S(:),Data_T(:),DData_S(:)
	real, allocatable :: time(:)
	real, allocatable :: tau(:),Int_0(:),C_C(:)

	allocate(Data_S(nt),DData_S(nt),Data_T(nt))
	allocate(time(nt))

	vtau=0.

	do k=1,nt
		time(k)=(k-1)*dt
	enddo

	Tp=time(nt)	!periodo de integracion (-Tp,Tp)

	ntau=2*(nt-1)+1

	allocate(tau(ntau),Int_0(ntau),C_C(ntau))

        do itau=1,ntau
	        tau(itau)=-Tp+(2.*Tp/(ntau-1))*(itau-1)
	enddo

        delta=ceiling(time_period/dt)	!!rango de búsqueda

	tau0=0
	
	do j=1,NumRec
	
		Data_S(:)=Data_Synth(:,j)
                call derivative(nt,time,Data_S,DData_S)

		Data_T(:)=Data_trace(:,j)

!!!!		CALCULO DE LA FUNCION CROSSCORRELATION 

		call CC_function(j,tau,Data_S,Data_T,DData_S,nt,dt,ntau,C_C,Int_0)

!!!!		BUSQUEDA TAU (cero de Int0 y el maximo de la funcion C_C) 
		
		call tau_search(j,nt,ntau,delta,tau,C_C,Int_0,tau0,imax) !1+2

		vtau(j)=tau0

!!!		write(*,*)j,vtau(j)

                Misfit=Misfit+0.5*tau0**2.d0

!!!!            ADJOINT SOURCE - ESTRATEGIA FISCHNER

                Iu22=0.
                u22=0.

                do k=1,nt
	                Iu22=Iu22+dt*DData_S(k)**2.d0
                enddo

                u22=sqrt(Iu22)

                do k=1,nt
	                i=nt-k+1
        	        AdjSource(k,j)=-tau0*DData_S(i)/u22**2.d0
                enddo


        enddo!j, NumRec
        
	deallocate(Data_S,DData_S,Data_T)
	deallocate(time)
	deallocate(tau,Int_0,C_C)

return
end
