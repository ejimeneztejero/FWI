subroutine polyfit(nut,Misfit, Misfit_test1, Misfit_test,loopind,EndStep)
implicit none

	include 'mpif.h'
	integer :: numtasks, rank, ierr, status(MPI_STATUS_SIZE)
	integer :: loopind ,j ,i ,k ,length
	real :: Misfit ,Misfit_test1 ,Misfit_test, nut, EndStep
	double precision :: vec_nut(3) , mean, std, Vec_Misfit(3) ,MinMisfit
	double precision :: Coef(3) ,A(3,3) ,y(3)
	integer,allocatable :: INDX(:)
	double precision, allocatable :: Nut_int(:),Poly_int(:)

	call MPI_COMM_SIZE(MPI_COMM_WORLD,numtasks,ierr)
	call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)

	length=2001
	
	allocate(INDX(3),Nut_int(length),Poly_int(length))
	
	if (loopind .eq. 2) then
		
		vec_nut(1)=0.
		vec_nut(2)=0.5*nut
		vec_nut(3)=nut
		
		mean=sum(vec_nut)/3.
		std=(1./2.*sum((vec_nut-mean)**2))**0.5
		vec_nut=(vec_nut-mean)/std
		
		Vec_Misfit(1)=Misfit
		Vec_Misfit(2)=Misfit_test1
		Vec_Misfit(3)=Misfit_test
		
		do j=1,3
			do i=1,3
				
				A(j,i)=vec_nut(j)**(i-1)
			enddo 
		enddo
				
		do i=1,3
			y(i)=Vec_Misfit(i)
		enddo
					
		call LEGS(A,3,y,Coef,INDX)
					
		do j=1,length
			Nut_int(j)=(j-1)*nut/(length-1)
		enddo
						
		Poly_int=((Nut_int-mean)/std)**2*Coef(3)+((Nut_int-mean)/std)*Coef(2)+Coef(1)
						
		MinMisfit=minval(Poly_int)
						
		do j=1,length
			if (MinMisfit .eq. Poly_int(j)) then
				EndStep=Nut_int(j)
			endif
		enddo
		
	elseif (loopind .gt. 2) then	
		
		vec_nut(1)=0.25*nut
		vec_nut(2)=0.5*nut
		vec_nut(3)=nut
		
 		mean=sum(vec_nut)/3.
 		std=(1./2.*sum((vec_nut-mean)**2))**0.5
		vec_nut=(vec_nut-mean)/std
		
		Vec_Misfit(1)=Misfit
		Vec_Misfit(2)=Misfit_test1
		Vec_Misfit(3)=Misfit_test
		
		do j=1,3
			do i=1,3
				
				A(j,i)=vec_nut(j)**(i-1)
				
			enddo 
		enddo
		
		do i=1,3
			y(i)=Vec_Misfit(i)
		enddo
					
		call LEGS(A,3,y,Coef,INDX)
					
		do j=1,length
			Nut_int(j)=0.25*nut+(j-1)*(0.75*nut)/(length-1)
		enddo
		
		Poly_int=((Nut_int-mean)/std)**2*Coef(3)+((Nut_int-mean)/std)*Coef(2)+Coef(1)
		
		MinMisfit=minval(Poly_int)
		
		
		do j=1,length
			if (MinMisfit .eq. Poly_int(j)) then
				EndStep=Nut_int(j)
			endif
		enddo
		
	end if
	if(rank.eq.0)write(*,*) 'Expected misfit ='
     	if(rank.eq.0)write(*,*) MinMisfit

return	
end 

	SUBROUTINE LEGS (A,N,B,X,INDX)
	!
	! Subroutine to solve the equation A(N,N)*X(N) = B(N) with the
	! partial-pivoting Gaussian elimination scheme.
	! Copyright (c) Tao Pang 2001.
	!
	IMPLICIT NONE
	INTEGER, INTENT (IN) :: N
	INTEGER :: I,J
	INTEGER, INTENT (OUT), DIMENSION (N) :: INDX
	double precision, INTENT (INOUT), DIMENSION (N,N) :: A
	double precision, INTENT (INOUT), DIMENSION (N) :: B
	double precision,INTENT (OUT), DIMENSION (N) :: X
	!
	CALL ELGS (A,N,INDX)
	!
	DO I = 1, N-1
		DO J = I+1, N
			B(INDX(J)) = B(INDX(J))-A(INDX(J),I)*B(INDX(I))
		END DO
	END DO
	!
	X(N) = B(INDX(N))/A(INDX(N),N)
	DO I = N-1, 1, -1
		X(I) = B(INDX(I))
		DO J = I+1, N
			X(I) = X(I)-A(INDX(I),J)*X(J)
		END DO
		X(I) =  X(I)/A(INDX(I),I)
	END DO
	!
	END SUBROUTINE LEGS
	!
	SUBROUTINE ELGS (A,N,INDX)
	!
	! Subroutine to perform the partial-pivoting Gaussian elimination.
	! A(N,N) is the original matrix in the input and transformed matrix
	! plus the pivoting element ratios below the diagonal in the output.
	! INDX(N) records the pivoting order.  Copyright (c) Tao Pang 2001.
	!
	IMPLICIT NONE
	INTEGER, INTENT (IN) :: N
	INTEGER :: I,J,K,ITMP
	INTEGER, INTENT (OUT), DIMENSION (N) :: INDX
	double precision :: C1,PI,PI1,PJ
	double precision, INTENT (INOUT), DIMENSION (N,N) :: A
	double precision, DIMENSION (N) :: C
	!
	! Initialize the index
	!
	DO I = 1, N
	INDX(I) = I
	END DO
	!
	! Find the rescaling factors, one from each row
	!
	DO I = 1, N
	C1= 0.0
	DO J = 1, N
	C1 = max(C1,ABS(A(I,J)))
	END DO
	C(I) = C1
	END DO
	!
	! Search the pivoting (largest) element from each column
	!
	DO J = 1, N-1
	PI1 = 0.0
	DO I = J, N
	PI = ABS(A(INDX(I),J))/C(INDX(I))
	IF (PI.GT.PI1) THEN
	PI1 = PI
	K   = I
	ENDIF
	END DO
	!
	! Interchange the rows via INDX(N) to record pivoting order
	!
	ITMP    = INDX(J)
	INDX(J) = INDX(K)
	INDX(K) = ITMP
	DO I = J+1, N
	PJ  = A(INDX(I),J)/A(INDX(J),J)
	!
	! Record pivoting ratios below the diagonal
	!
	A(INDX(I),J) = PJ
	!
	! Modify other elements accordingly
	!
	DO K = J+1, N
	A(INDX(I),K) = A(INDX(I),K)-PJ*A(INDX(J),K)
	END DO
	END DO
	END DO
	!
	END SUBROUTINE ELGS
