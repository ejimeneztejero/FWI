subroutine CC_function(jrec,tau,Data_S,Data_T,DData_S,nt,dt,ntau,C_C,Int_0)
implicit none

	integer :: NumRec,nt,ntau,jrec
	integer :: itau,k1,k2,kk,iplus,k
	real :: dt,CC,Int0
	real :: Data_T(nt),Data_S(nt),DData_S(nt)
	real :: C_C(ntau),Int_0(ntau),tau(ntau)

        do itau=1,ntau			

		k1=1
		k2=nt	!no tocar
		kk=0
		iplus=1

		call f1f2correlation(iplus,tau(itau),nt,dt,k1,k2,kk)!no tocar

		CC=0.
		Int0=0.

		do k=k1,k2
			if((k+kk).ge.1.and.(k+kk).le.nt)	then
				CC=CC+dt*Data_T(k)*Data_S(k+kk)
				Int0=Int0+dt*Data_T(k)*DData_S(k+kk)
			endif
		enddo

		C_C(itau)=CC
		Int_0(itau)=Int0

!!!!!		write(555,*)jrec,itau,CC,Int0

	enddo!tau

return
end
