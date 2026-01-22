	subroutine space_deriv2x(acc,coef,nyt,nxt,p,deriv2_px)

	integer nyt,nxt,acc
	integer i,i1,i2,j,j1,j2,k
	real p(nyt,nxt),deriv2_px(nyt,nxt)
	real coef(acc+1)

	j1=1+acc
	j2=nxt-acc

	do j=j1,j2
                deriv2_px(:,j)=0.
                do i=1,acc+1
                        deriv2_px(:,j)=deriv2_px(:,j)+coef(i)*(p(:,j+(i-1))+p(:,j-(i-1)) )
                enddo
	enddo

        do j=1,j1-1
                deriv2_px(:,j)=deriv2_px(:,acc+1)
        enddo
	do j=j2+1,nxt
                deriv2_px(:,j)=deriv2_px(:,nxt-acc)
        enddo

	return
	end
	
	subroutine space_deriv2y(acc,coef,nyt,nxt,p,deriv2_py)

	integer nyt,nxt,acc
	integer i,i1,i2,j,j1,j2,k
	real p(nyt,nxt),deriv2_py(nyt,nxt)
	real coef(acc+1)

	j1=1+acc
	j2=nyt-acc

	do k=1,nxt
		do j=j1,j2
			deriv2_py(j,k)=0.
			do i=1,acc+1
				deriv2_py(j,k)=deriv2_py(j,k)+(p(j+(i-1),k)+p(j-(i-1),k))*coef(i)
			enddo
		enddo
	enddo

	do k=1,nxt
		do j=1,j1-1
			deriv2_py(j,k)=deriv2_py(acc+1,k)
		enddo

		do j=j2+1,nyt
			deriv2_py(j,k)=deriv2_py(nyt-acc,k)
		enddo
	enddo

	return
	end

	subroutine space_derivx(acc,coef,nyt,nxt,p,deriv_px)

	integer nyt,nxt,acc
	integer i,i1,i2,j,j1,j2
	real p(nyt,nxt),deriv_px(nyt,nxt)
	real coef(acc+1)

	j1=1+acc
	j2=nxt-acc

	do j=j1,j2

		do k=1,nyt
			deriv_px(k,j)=0.
		enddo

		do i=1,acc+1
		do k=1,nyt
			deriv_px(k,j)=deriv_px(k,j)+(p(k,j+(i-1))-p(k,j-(i-1)))*coef(i)
		enddo
		enddo

	enddo

	do j=1,j1-1
	do k=1,nyt
		deriv_px(k,j)=deriv_px(k,acc+1)
	enddo
	enddo

	do j=j2+1,nxt
	do k=1,nyt
		deriv_px(k,j)=deriv_px(k,nxt-acc)
	enddo
	enddo				    

	return
	end
		
	subroutine space_derivy(acc,coef,nyt,nxt,p,deriv_py)

	integer nyt,nxt,acc
	integer i,i1,i2,j,j1,j2
	real p(nyt,nxt),deriv_py(nyt,nxt)
	real coef(acc+1)
	
	j1=1+acc
	j2=nyt-acc

	do k=1,nxt
		do j=j1,j2
			deriv_py(j,k)=0.
			do i=1,acc+1
				deriv_py(j,k)=deriv_py(j,k)+(p(j+(i-1),k)-p(j-(i-1),k))*coef(i)
			enddo
		enddo
	enddo

	do k=1,nxt
		do j=1,j1-1
			deriv_py(j,k)=deriv_py(acc+1,k)
		enddo
	enddo

	do k=1,nxt
		do j=j2+1,nyt
			deriv_py(j,k)=deriv_py(nyt-acc,k)
		enddo
	enddo
	
	return
	end
