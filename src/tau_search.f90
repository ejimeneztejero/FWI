subroutine tau_search(j,nt,ntau,delta,tau,C_C,Int_0,tau0,imax)
implicit none

	integer :: j,nt,ntau,indfile
	integer :: it,itau,imax,delta,delta_tau
	real :: tau0,ymax,difm,tau_old,xmax,dif
	real :: tau(ntau),C_C(ntau),Int_0(ntau)

		imax=0
		ymax=-1e6
		it=0
		tau_old=0
		difm=100000.

		delta_tau=delta/2

	        do itau=nt-delta_tau,nt+delta_tau

			xmax=(C_C(itau)+C_C(itau+1))/2

			if(Int_0(itau)*Int_0(itau+1).lt.0.and.xmax.gt.0)then
				imax=imax+1!num de maximos
				
				if(xmax.gt.ymax)then 
					if(j.eq.1)	then
						!it=itau
						!ymax=xmax
						if(tau_old.eq.0)then
							it=itau
							ymax=xmax
							tau_old=tau(itau)
						endif
						if(tau_old.ne.0.and.abs(tau(itau)).lt.abs(tau_old))then
							it=itau
							ymax=xmax
							tau_old=tau(itau)	
						endif
					endif
					if(j.ne.1)	then
						dif=abs(tau(itau)-tau0)	!!criterio 2: la minima diferencia entre los tiempos de la actual j a la anterior (j-1)
						if(dif.lt.difm)then
							it=itau
							ymax=xmax
							difm=dif
						endif
					endif


				endif
			endif

		enddo

!!!!!!!!	calculo tau0
	
	if(imax.ne.0)tau0=tau(it)
	

return
end
