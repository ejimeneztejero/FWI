subroutine f1f2correlation(iplus,ttau,nt,dt,k1,k2,kk)
implicit none

	integer :: iplus,nt,k1,k2,kk,nttau
	real :: ttau,dt

	nttau=ceiling(abs(ttau)/dt)+1 !puntos en el segmento ttau con resolution dt

!!!!!!!!	t+T

	if(iplus.eq.1)then

		if(ttau.gt.0)then						
			k1=nttau	
			kk=1-nttau	! negativo

		elseif(ttau.lt.0)then			
			k2=nt-nttau+1	
			kk=-1+nttau	! positivo
		endif

!!!!!!!!	t-T

	else

		if(ttau.lt.0)then						
			k1=nttau
			kk=1-nttau

		elseif(ttau.gt.0)then			
			k2=nt-nttau+1
			kk=-1+nttau

		endif

	endif

return
end
