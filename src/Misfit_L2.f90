subroutine Misfit_L2(Data_Synth,Data_trace,NumRec,nt,dt,Misfit,AdjSource)
implicit none
	
	integer :: nt,NumRec,j,k
	real :: dt ,Misfit
	real :: Data_trace(nt,NumRec),Data_Synth(nt,NumRec),AdjSource(nt,NumRec)

	do j=1,NumRec
		do k=1,nt
			Misfit=Misfit+0.5*dt*(Data_Synth(k,j)-Data_trace(k,j))**2
			AdjSource(k,j)=-(Data_Synth(nt-k+1,j)-Data_trace(nt-k+1,j))
		enddo
!                write(22,*) j,dt,maxval(abs(Data_trace(:,j))),maxval(abs(Data_Synth(:,j))),Misfit
	enddo
!	write(*,*)"Misfit L2, adj source",AdjSource(1,1)

return	
end
