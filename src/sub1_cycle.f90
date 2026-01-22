subroutine sub1_cycle(num,jrec,nt,dt,ntt,S,env,triang,S_n,DQ)

implicit none

	integer :: num,nt,ntt,k,l,m,j,jrec
	real :: dt,S(nt),DQ(nt)
	real :: env(nt),triang(nt),S_n(nt)
	real, allocatable :: time(:)
	
	env=0.;triang=0.;S_n=0.

	allocate(time(nt))

	do k=1,nt
		time(k)=(k-1)*dt
	enddo

!!	CÁLCULO señal triangular directamente, S_n
	call sub2_cycle(num,jrec,nt,ntt,dt,time,S,triang,S_n,env,DQ)

!!	CÁLCULO ENVELOPE
!	env=abs(S)
!	do k=1,nt
!		if(abs(S_n(k)).gt.0)env(k)=abs(S(k)/S_n(k))
!	enddo

	deallocate(time)

return
end
