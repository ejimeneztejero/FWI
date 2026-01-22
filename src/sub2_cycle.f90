subroutine sub2_cycle(num,jrec,nt,ntt,dt,time,S,triang,S_n,env,DQ)
!!	CÁLCULO señal triangular directamente

implicit none

	real, parameter :: pi = 3.14159265

	integer nt,ntt,k,l,entra,jrec,num
	integer nl,nl_1,nl_l,nl_i,nl_f,cero_i,cero_f

	real :: dt,time(nt),S(nt),triang(nt),S_n(nt),env(nt),DQ(nt)
	real, allocatable :: DS(:),DDS(:),xc(:),yc(:)
	real, allocatable :: xl(:),yl(:),tml(:),triang2(:)
	real, allocatable :: ye(:),yl_e(:),env2(:) !!dimensiones en l
	
	allocate(DS(nt),DDS(nt),xc(nt),yc(nt),ye(nt))
	DS=0.;DDS=0.;xc=0.;yc=0.;ye=0;
	
	env=abs(S)

!!!	primera y segunda derivada del valor absoluto de la señal

	call derivative(nt,time,S,DS) 
	call derivative(nt,time,DS,DDS)

	l=0;entra=0
	do k=ntt,nt-10
	
!!!!		if(num.eq.1.and.jrec.eq.20)write(111,*)(k-1)*dt,S(k),DS(k),DDS(k)

		if(DS(k)*DS(k+1).lt.0.or.(DS(k)*DS(k+1).eq.0.and.DS(k)*DS(k+2).lt.0.)) then !!punto critico, maximo o minimo

		entra=entra+1

		!minimo
                if(DDS(k).gt.0) then

			l=l+1		

			xc(l)=time(k)	
			yc(l)=pi/2.
			ye(l)=abs(S(k))

		endif
		
		!máximo
                if(DDS(k).lt.0) then

			l=l+1

			xc(l)=time(k)	
			yc(l)=-pi/2.
			ye(l)=abs(S(k))

		endif
				
		endif
					
	enddo


	if(l.eq.0.or.entra.eq.0) then
		write(*,*)'No hay datos en jrec: ',jrec
	endif

	if(l.ne.0)	then

		allocate(xl(l),yl(l),yl_e(l))
		do k=1,l
			xl(k)=xc(k)
			yl(k)=yc(k)
			yl_e(k)=ye(k)
		enddo

		nl_1=xl(1)/dt+1 
		nl_l=xl(l)/dt+1
	
		nl_i=xl(2)/dt+1	 ! segundo max or min
		nl_f=xl(l-1)/dt+1! penúltimo max or min
	
		nl=nl_l-nl_1+1 

		allocate(tml(nl),triang2(nl),env2(nl))
		triang2=0.
		do k=1,nl
			tml(k)=time(k+nl_1-1)
		enddo

		call interpola_lineal1D(l,xl,yl,nl,tml,triang2)
		call interpola_lineal1D(l,xl,yl_e,nl,tml,env2)


		do k=nl_1,nl_l
			env(k)=env2(k-nl_1+1) !envelope
			triang(k)=triang2(k-nl_1+1) !Triangular signal
			S_n(k)=-sin(triang(k))	!normalized signal
		enddo

		cero_i=0;cero_f=0
	
!!		Para que las señales empiece en cero		
		call fun_ceros(nt,nl_1,nl_l,nl_i,nl_f,S_n,cero_i,cero_f)	

		if(cero_i.gt.1.and.cero_f.lt.nt) then
			triang(1:cero_i)=0.
			triang(cero_f:nt)=0.
			S_n(1:cero_i)=0.
			S_n(cero_f:nt)=0.	
		endif

!!      	CÁLCULO DQ, cuadratura?
        	call fun_DQ(nt,nl_1,nl_l,nl_i,nl_f,S_n,DQ)

		deallocate(xl,yl,tml,triang2)!entre nL1 y nL

	endif

	deallocate(DS,DDS,xc,yc)
			
return
end
