subroutine time_filter(Data_trace,nt,dt,Num,typef,f1,fc)
implicit none

	integer :: nt,Num,nt2,zero_pad,ind,typef
	real :: dt,fe,var,fc,f1
	real :: Data_trace(nt,Num)
 	real,allocatable :: tmp(:)
	integer :: j,i,k
	logical :: isnotCero     	! Variable para almacenar si el vector es cero

	zero_pad=24	
	nt2=zero_pad*nt
	allocate(tmp(nt2))
	fe=1/dt

	do i=1,Num

		tmp=0.
		tmp((zero_pad-1)*nt+1:nt2)=Data_trace(:,i)

		call is_not_cero(nt,Data_trace(:,i),isnotCero)

		if(isnotCero)	then

!			do j=1,nt
!				write(55,*)j,Data_trace(j,i)
!			enddo
			call filters_function1D(tmp,nt2,fe,1,typef,f1,fc)
			Data_trace(:,i)=tmp(nt2-nt+1:nt2)


		endif

	enddo

	
	deallocate(tmp)

return
end
