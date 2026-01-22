subroutine interpola_lineal1D(l,xc,yc,nl,tml,S_e2)

implicit none

        integer i,j,k,l,nl
	real :: xc(l),yc(l),tml(nl),S_e2(nl),m

	do i=1,l
        do j=1,nl

	m=(yc(i+1)-yc(i))/(xc(i+1)-xc(i))

	if (tml(j).ge.xc(i).and.tml(j).lt.xc(i+1)) then
               	S_e2(j)=yc(i)+(tml(j)-xc(i))*m
	end if

	enddo
	enddo

        do j=1,nl
		S_e2(j)=S_e2(j)+0.001
	enddo

return
end

