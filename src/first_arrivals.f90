subroutine first_arrivals(Data_trace,nt,ntt)
implicit none		

	integer :: nt,FAj
	integer :: j,k,ntt,l
	real :: Data_trace(nt)
        real, allocatable :: T(:) 

	allocate(T(nt))

!!!	PRIMERAS LLEGADAS

	T=Data_trace/maxval(abs(Data_trace(:)))

	call F_A(T,nt,ntt)!FAT, el tiempo donde la señal comienza a no ser cero

	deallocate(T)


return
end
