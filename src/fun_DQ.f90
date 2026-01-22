subroutine fun_DQ(nt,nl_1,nl_l,nl_i,nl_f,S_n,DQ)

implicit none

!!	CÁLCULO DQ

	integer :: nt,signo_q,k,nl_1,nl_l,nl_i,nl_f
	real :: dif_x1,S_n(nt),DQ(nt)
	
!!	signo coseno
	signo_q=0
!	do k=nl_i,nl_f !!son estos límites a partir del segundo maximo hasta el penúltimo, para que la señal empiece en cero
	do k=nl_1,nl_l-1

!		dif_x1=S2(k)-S2(k-1)
!		if(dif_x1.ge.0)signo_q=1
!		if(dif_x1.lt.0)signo_q=-1
!		phase(k)=triang(k)/(-signo_q)

		dif_x1=S_n(k+1)-S_n(k)
		if(dif_x1.ge.0)signo_q=1
		if(dif_x1.lt.0)signo_q=-1
		DQ(k)=-signo_q*sqrt(abs(1.-S_n(k)**2.))
	
	enddo
	
end subroutine fun_DQ
