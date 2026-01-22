subroutine F_A(S,nt,FA)
implicit none

	integer :: nt,k,FA
	real :: umbral,S(nt)		!!! tiene que estar normalizada. 	
	
	umbral=0.001
	do k=3,nt-2

		if(abs(S(k)).ge.umbral)	then
		if(abs(S(k+1)).ge.umbral.and.abs(S(k-1)).ge.umbral)	then
		if(abs(S(k+2)).ge.umbral.and.abs(S(k-2)).ge.umbral)	then
			FA=k
			goto 22
		endif
		endif
		endif	
		
	enddo
	22 continue

	
return
end

