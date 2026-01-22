!!	calculo el primer y el último cero
subroutine fun_ceros(nt,nl_1,nl_l,nl_i,nl_f,S,cero_i,cero_f)
implicit none

	integer :: k,nt,nl_1,nl_l,nl_i,nl_f,cero_i,cero_f
	real :: S(nt)

	cero_i=1
	cero_f=nt
	
	do k=nl_1,nl_i-1
		if(S(k)*S(k+1).lt.0)	then
			cero_i=k+1
			goto 132
		endif
	enddo
	132 continue

	do k=nl_f,nl_l
		if(S(nl_l-k+nl_f)*S(nl_l-k+nl_f-1).lt.0)	then
			cero_f=k
			goto 243
		endif
	enddo
	243 continue
	
!	write(*,*)'n1,n2,nl-1,nl',nl_1*dt,nl_i*dt,nl_f*dt,nl_l*dt
!	write(*,*)'FA,cero i,cero f',ntt*dt,cero_i*dt,cero_f*dt

return
end
