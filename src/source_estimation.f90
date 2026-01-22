subroutine source_estimation(iter,f1,f2,nx,ny,nt,nr,dmodel,dt,&
nxs,nys,nxr,nyr,model,sst)

implicit none

integer :: rank,i,j,k,ir,iw,ier,ixm,iter
integer :: nx,ny,nt,nr,nt_pad
integer :: it1,it2,if1,if2
integer :: nxs,nys
integer :: nxr(nr),nyr(nr)
integer	:: lensav, lenwrk

real	:: xmid,t_wb,zwb,z_wb(nx)
real	:: dt,t1,t2,dt1,dt2,tdir
real	:: f1,f2,fc,df
real	:: dmodel
real	:: eps
real	:: model(ny,nx)
real	:: fdata(nt,nr)
real	:: sst(nt)

real, allocatable :: wsave(:), work(:)

real, allocatable :: delta(:),s0(:)
real, allocatable :: Gs(:,:)
real, allocatable :: Go(:,:)

complex, allocatable :: Gws(:,:),Gwo(:,:)
complex, allocatable :: kk(:),bb(:)
complex, allocatable :: ss(:)

allocate(Go(nt,nr))
allocate(Gs(nt,nr))
allocate(delta(nt))
allocate(s0(nt))
S0=0.;Gs=0.;delta=0;
delta(1)=1/dt;

call get_shot_data(1,Go)  !! raw_data SU-MCS files

if(iter.eq.1)S0=delta
!if(iter.eq.1)call Ricker_mcs(rank,1,S0)
!!if(iter.gt.1)s0=sst

!!write(54,*)iter,delta(:)
!!write(55,*)iter,S0(:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

nt_pad=nt*2

allocate(Gws(nt_pad,nr),Gwo(nt_pad,nr))
allocate(kk(nt_pad),bb(nt_pad),ss(nt_pad))
Gws=0.;Gwo=0.;
kk=(0.0,0.0);bb=(0.0,0.0);ss=(0.0,0.0);

lensav = 2*nt_pad + ceiling(log(real(nt_pad))) + 4
lenwrk = 2*nt_pad

allocate(wsave(lensav),work(lenwrk))
call cfft1i(nt_pad, wsave, lensav, ier)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

df=1./(nt*dt)
if1=1+floor(f1/df);if2=1+floor(f2/df)

fc=(f1+f2)/2.
dt2=3./fc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Generacion de shots sinteticos con impulso delta:
!! Tiempo maximo para propagar solo en agua timew(:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do i = 1, nx
z_wb(i) = (ny-1)*dmodel   ! default
   do j = 1, ny
      if (model(j,i) > 1500.0) then
         z_wb(i) = (j-1)*dmodel
         exit
      endif
   enddo
enddo

call solver_forward_simple(nr,nxs,nys,nxr,nyr,s0,dmodel,dt,ny,nx,model,nt,Gs)

!! Ventana para G0 y Gs =0 en (1,it2):
do ir=1,nr

	xmid = 0.5*(nxr(ir) + nxs) * dmodel
	ixm  = int(xmid/dmodel) + 1
	zwb = z_wb(ixm)

	tdir = abs(nxr(ir)-nxs) * dmodel / 1500.0
	t_wb = 2.0 * sqrt( (abs(nxr(ir)-nxs)*dmodel*0.5)**2 + zwb**2 ) / 1500.0

!!	t2 = min( tdir + dt2, t_wb - 0.5/fc )
	t2 = t_wb - 0.5/fc
	it2 = min(nt, 1 + int(t2/dt))

	Gs(it2+1:nt,ir) = 0.0
	Go(it2+1:nt,ir) = 0.0

!!	write(66,*)ir,it2,t2

enddo

Gs=Gs/maxval(abs(Gs))
Go=Go/maxval(abs(Go))

if(iter.eq.1)open(unit=12,file="Gs_Go_01.dat",status='unknown')
if(iter.eq.2)open(unit=12,file="Gs_Go_02.dat",status='unknown')
do j=1,nt
       write(12,*) dt*(j-1),Gs(j,1),Go(j,1)
enddo
close(12)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! FFT para G0 y Gs:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Gwo = (0.0,0.0);Gws = (0.0,0.0)

do ir = 1, nr

	Gwo(1:nt,ir) = Go(1:nt,ir)
	call cfft1f(nt_pad,1, Gwo(:,ir), nt_pad,wsave, lensav, work, lenwrk, ier)

	Gws(1:nt,ir) = Gs(1:nt,ir)
	call cfft1f(nt_pad,1, Gws(:,ir), nt_pad,wsave, lensav, work, lenwrk, ier)

enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! From here it depends on BandPass (f1,f2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculo de kk(w), bb(w) y s(w), dentro de banda (if1,if2):
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do iw=if1,if2

	kk(iw)=0
	bb(iw)=0
	
	do ir=1,nr
		kk(iw)=kk(iw)+conjg(Gws(iw,ir))*Gws(iw,ir)
		bb(iw)=bb(iw)+conjg(Gws(iw,ir))*Gwo(iw,ir)
	enddo

enddo

eps=0.001*maxval(real(kk(if1:if2)))
do iw=if1,if2
	ss(iw)=bb(iw)/(kk(iw)+eps)
enddo

!do iw = if1, if2
!   ss(nt_pad - iw + 2) = conjg(ss(iw))
!enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! IFFT para s(t)=IFFT(s(w)) :
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call cfft1b(nt_pad,1,ss,nt_pad,wsave,lensav,work,lenwrk,ier)

sst(:) = real(ss(1:nt))/real(nt_pad)

!open(unit=12,file="source_estimation.dat",status='unknown')
!do j=1,nt
!	write(12,*) dt*(j-1),sst(j)
!enddo
!close(12)

deallocate(kk,bb,ss)
deallocate(Go)
deallocate(Gs)
deallocate(Gws,Gwo)
deallocate(delta)

deallocate(wsave,work)

end subroutine source_estimation
