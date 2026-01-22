    !*********************************************************************/
    !*********************************************************************/
    !*********************************************************************/


!        call fourier_1D(fs(:,is),sou_Fourier,w,nt,dt,zero_pad,nt_pad)

    subroutine fourier_1D( trace, trace_Fourier, w, nt, dt, zero_pad,nt_pad )

      implicit none

      ! The variables which are passed to the function.
      integer, intent(in)      :: nt
      integer, intent(in)      :: nt_pad,zero_pad
      real, intent(in)     :: dt
      real, intent(in)     :: trace(nt)
      real, intent(out)    :: w(nt_pad)
      complex, intent(out) :: trace_Fourier(nt_pad)

      ! The variables which are generated inside the function.
      integer :: k, lensav, lenwrk, iaux
      real, allocatable :: wsave(:), work(:)
      real :: dw
      integer :: stal
      integer :: ind_middle

      !---------------------------------------------------------------------/
      !
      if ( mod(nt_pad,2) /= 0 ) then
        ind_middle = (nt_pad-1)/2
      else
	ind_middle = nt_pad/2
      end if

      !---------------------------------------------------------------------/
      !
      dw = 1./( dt*dble(nt_pad) )
      w(1) = 0.
      do k=2,ind_middle+1
        w(k)          =  dble(k-1)*dw
        w(nt_pad-k+2) = -dble(k-1)*dw
      end do

      !---------------------------------------------------------------------/
      !
!      lensav = 2*nt_pad + ceiling(dlog(dble(nt_pad))/dlog(2.0d0)) + 4	!dani
      lensav=2*nt_pad+ceiling(log(nt_pad/1.))+4	!yo
      lenwrk = 2*nt_pad

      allocate(wsave(lensav), work(lenwrk))

      call cfft1i(nt_pad,wsave,lensav,iaux)
      !---------------------------------------------------------------------/
      !
      trace_Fourier(1:nt) = dcmplx( trace )
      if ( nt_pad /= nt ) trace_Fourier(nt+1:nt_pad) = 0.

      call cfft1f(nt_pad,1,trace_Fourier,nt_pad,wsave,lensav,work,lenwrk,iaux)

      !---------------------------------------------------------------------/
      deallocate( wsave, work)

    end subroutine fourier_1D

    !*********************************************************************/
    !*********************************************************************/
    !*********************************************************************/
  
  subroutine inv_fourier_1D( trace_Fourier, trace, nt, nt_pad )
      implicit none

      ! The variables which are passed to the function.
      integer, intent(in)     :: nt
      integer, intent(in)     :: nt_pad
      complex, intent(in) :: trace_Fourier(nt_pad)
      real, intent(out)   :: trace(nt)

      ! The variables which are generated inside the function.
      integer :: k
      integer :: lensav, lenwrk, iaux
      real, allocatable :: wsave(:), work(:)
      integer :: stal
      complex :: trace_Fourier_aux(nt_pad)

      !---------------------------------------------------------------------/
      !
      lensav=2*nt_pad+ceiling(log(nt_pad/1.))+4	!yo
!      lensav = 2*nt_pad + ceiling(dlog(dble(nt_pad))/dlog(2.0d0)) + 4	!dani
      lenwrk = 2*nt_pad

      allocate( wsave(lensav), work(lenwrk), stat=stal ); if ( stal/=0 ) stop ' ***** AE-1: inv_fourier_trans_1D ***** '

      call cfft1i(nt_pad,wsave,lensav,iaux)

      !---------------------------------------------------------------------/
      trace_Fourier_aux = trace_Fourier
      call cfft1b(nt_pad,1,trace_Fourier_aux,nt_pad,wsave,lensav,work,lenwrk,iaux)

      trace = real(trace_Fourier_aux(1:nt))/nt_pad

      !---------------------------------------------------------------------/
      deallocate( wsave, work, stat=stal ); if ( stal/=0 ) stop ' ***** AE-1: inv_fourier_trans_1D ***** '

    end subroutine inv_fourier_1D
