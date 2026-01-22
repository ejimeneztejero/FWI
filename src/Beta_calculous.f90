subroutine Beta_calculous(Grad,OldGrad,nx,ny,Beta)
implicit none

	integer :: nx,ny,i
	real :: Grad(ny,nx) ,OldGrad(ny,nx) ,Beta ,Num, Den
	real,allocatable :: Vec(:),VecOld(:)
	
	allocate(Vec(nx),VecOld(nx))
	
	Num=0.
	Den=0.
	
	do i=1,ny
		
		Vec=Grad(i,1:nx)
		VecOld=OldGrad(i,1:nx)
		Num=Num+dot_product(Vec,(Vec-VecOld))		
		Den=Den+dot_product(VecOld,VecOld)
		
	enddo
	
	Beta=Num/Den
	
return
end
