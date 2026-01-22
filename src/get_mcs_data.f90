!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!	This file contains subroutines for extracting input data features
!!	Author:	Clara Estela Jimenez Tejero. 
!!	License: GNU General Public License v3.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_mcs_data(rank)

USE mod_parfile
USE mod_data_arrays

implicit none

integer :: i,j,k,dif
integer :: rank,ishot
character(len=500) :: file_name
logical :: file_exists

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(rank.eq.0)    then
        write(*,*)
        write(*,*)'*************************************'
        write(*,*)'CHECK MCS DATA'
        write(*,*)'*************************************'
endif

!! FOR SHOTS-MCS
call check_nav_shot_mcs(rank) !perhaps ony rank=0. Intentar que solo haya warnings y errores
call read_nav_shot_mcs(rank)

call open_shot_mcs(rank)     !!opens unit file for MCS data
call get_pos_channels(rank)    !!ojo, repasar

!if(rank.eq.0)   then
!do i=1,NumShots_mcs
!        if(shotID_su_mcs(i).ne.shotID_nav_mcs(i))   then
!                write(*,*)'ERROR: shotID must be coincident in su and nav file'
!		write(*,*)i,shotID_su_mcs(i),shotID_nav_mcs(i)
!                call ascii_art(2)
!                stop
!        endif
!enddo
!endif

if(rank.eq.0)write(*,*)'FOR MCS: nx,length: ',nxmodel_mcs,length_model_mcs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(rank.eq.0)    then
        write(*,*)
        write(*,*)'*************************************'
        write(*,*)'GEOMETRY FOR MCS DATA'
        write(*,*)'*************************************'
endif

!! GEOMETRY FOR INVERSION OF MCS DATA (OBSs ACT AS SOURCES, SHOTS AS RECEIVERS)
call get_geometry_mcs(rank) !! valores para nxSou_MCS,nySou_MCS,nxRec_MCS,nyRec_MCS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call allocate_model_arrays_mcs(rank)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(rank.eq.0)	then
        write(*,*)
        write(*,*)'*********************************'
      	write(*,*)'CALCULATION OF BATHYMETRY MODEL '
        write(*,*)'*********************************'
endif

call get_bathymetry_model_mcs(rank)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(rank.eq.0)	then
        write(*,*)
        write(*,*)'**************************'
        write(*,*)'CREATING SOURCE FUNCTION'
        write(*,*)'**************************'
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END subroutine get_mcs_data

subroutine read_nav_shot_mcs(rank)	

USE mod_parfile
USE mod_data_arrays

implicit none

integer :: i,j,k,rank
integer :: iline,itr,ERR,nlines
integer :: shotID
real :: xs,ys,bat,xs1,xsn
character (len=500) :: file_name

!!	LECTURA DATOS SHOTGATHERS
file_name= trim(adjustl(folder_input_mcs)) // trim(adjustl(nav_shot_mcs))
open(unit=10,file=file_name,status='old')

do i=1,NumShots_mcs !!number of lines nav_shot

    shotID=0;xs=0;ys=0.;bat=0.;
    read(10,*)shotID,xs,bat

    shotID_nav_mcs(i)=shotID

    pos_shot_mcs(i) = added_space_model+abs(xs)
    pos_shot_grid_mcs(i)=1+ceiling(pos_shot_mcs(i)/dmodel)

    pos_bat_mcs(i) = abs(bat)
    pos_bat_grid_mcs(i) = 1+ceiling(pos_bat_mcs(i)/dmodel)

!    write(111,*)pos_shot_mcs(i),pos_shot_grid_mcs(i),pos_bat_mcs(i),pos_bat_grid_mcs(i) 

    if(i.eq.1)xs1=xs
    if(i.eq.NumShots_mcs)xsn=xs

enddo
CLOSE(10)

   dshots_mcs=abs(xsn-xs1)/(NumShots_mcs-1)
   length_model_mcs=xs+2*added_space_model
   nxmodel_mcs=1+ceiling(length_model_mcs/dmodel)

if(rank.eq.0)    then
   write(*,*)"dshots_mcs: ",dshots_mcs
   if(length_model_mcs.gt.length_model)write(*,*)"WARNING: length_model_mcs>length_model"
   if(nxmodel_mcs.gt.nxmodel)write(*,*)"WARNING: nxmodel_mcs>nxmodel"
endif

end subroutine  read_nav_shot_mcs

subroutine get_geometry_mcs(rank)

USE mod_parfile
USE mod_data_arrays

implicit none
integer :: i,j,k,ishot,irec,rank

do ishot=1,NumShots_MCS

	nxSou_MCS(ishot)=pos_shot_grid_mcs(ishot)
	nySou_MCS(ishot)=1+ceiling(shot_depth_mcs/dmodel)

	do irec=1,NumRec_MCS
		nxRec_MCS(irec,ishot)=pos_rec_grid(irec,ishot)
		nyRec_MCS(irec,ishot)=1+ceiling(streamer_depth/dmodel)
	enddo

!	do k=1,NumRec_MCS
!	        write(50,*)pos_shot_mcs(ishot),pos_rec(k,ishot)
!	        write(51,*)nxSou_MCS(ishot),nySou_MCS(ishot),nxRec_MCS(k,ishot),nyRec_MCS(k,ishot)
!	enddo

enddo

end subroutine get_geometry_mcs


SUBROUTINE get_shot_data(ishot,raw_data)

USE mod_parfile
USE mod_data_arrays

IMPLICIT NONE

INTEGER :: nh
INTEGER :: i,j,k,ishot
INTEGER(4) :: pos_r
real    :: raw_data(nt_mcs,NumRec_MCS)

REAL(4), ALLOCATABLE :: Data_tmp(:,:),sudata(:,:)

nh = size_su_header

allocate(Data_tmp(nt_mcs,NumRec_MCS))
allocate(sudata(nt_mcs+nh,NumRec_MCS))
Data_tmp=0.;sudata=0.;
raw_Data=0.

!!! Loop trace by trace
do j=1,NumRec_MCS
        pos_r = 1 + (j-1)*(nh+nt_mcs)*4
        READ(unit_mcs+ishot, pos=pos_r)  sudata(1:nt_mcs+nh,j)
enddo

do j = 1,NumRec_MCS
        do i = 1, nt_mcs
                Data_tmp(i,j) = sudata(i+nh,j)
        enddo
enddo

raw_Data = Data_tmp

!!!     write(*,*)"SALE DE READ_OBS"

deallocate(Data_tmp)
deallocate(sudata)

END SUBROUTINE get_shot_data


subroutine check_nav_shot_mcs(rank)

use mod_parfile
USE mod_data_arrays

implicit none

        integer :: i,nlines,nums,rank
        integer :: shoti,shotf,nx_user
        real    :: xi,xf,xs,bat
        CHARACTER(len=50) :: access,form,Str
        CHARACTER(len=500) :: file_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!     NAVIGATION FILE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        file_name= trim(adjustl(folder_input_mcs)) // trim(adjustl(nav_shot_mcs))
        open(unit=10,file=file_name,status='old')
        nlines=0
        do
            read(10,*, END=10)
            nlines = nlines + 1
        enddo
        10 close (10)

        open(unit=10,file=file_name,status='old')
        do i=1,nlines

                read(10,*)nums,xs,bat

                if(i.eq.1)      then
                        shoti=nums
                        xi=xs
                endif

                if(i.eq.nlines) then
                        shotf=nums
                        xf=xs
                endif

!                shotID_nav_mcs(i)=nums

        enddo
        close(10)

        if(rank.eq.0)   then
                write(*,*)'Navigation file in mcs contains shots from ',shoti,' to ',shotf
                write(*,*)'NumShots total:', shotf-shoti+1,', Num Lines Nav_file:',nlines
        endif

        if(nlines.lt.NumShots_mcs)  then

                if(rank.eq.0)then
                        write(*,*)'ERROR: navigation, some shots are missing in nav_shot_mcs: ',nlines,NumShots_mcs
                        call ascii_art(2)
                        stop
                endif
        endif
        if(nlines.gt.NumShots_mcs)  then

                if(rank.eq.0)then
                        write(*,*)'WARNING: navigation contains more shots than wanted to be read: ',nlines,NumShots_mcs
                endif
        endif

        nx_user=1+ceiling(xf/dmodel)

        if(rank.eq.0)   then
                write(*,*)'nxmodel for MCS should be >=', nx_user
        endif

        if(nx_user.gt.nmodel)    then

            if(rank.eq.0)then
                write(*,*)
                write(*,*)"ERROR: fix you nxmodel. It should be equal to:"
                write(*,*)"1+distance(shotN-shot1)/dmodel = ", nx_user
                write(*,*)", where: distance(shotN-shot1) = ",xf-xi," and dmodel= ",dmodel
                call ascii_art(2)
                stop
            endif

    endif

end subroutine check_nav_shot_mcs

SUBROUTINE open_shot_mcs(rank)

USE mod_parfile
USE mod_data_arrays

IMPLICIT NONE
integer     :: i,j,k,ishot,rank
CHARACTER(len=50) :: OBS_num
CHARACTER(len=500) :: file_name

CHARACTER(len=50) :: access,form

access = 'stream'
form = 'unformatted'

do ishot=1,NumShots_mcs

file_name = trim(adjustl(folder_input_mcs)) // su_file_mcs(ishot)

if(endianness_machine.eq.0)     then
        if(endianness_data.eq.0)open(unit_mcs+ishot,FILE=file_name,ACCESS=access,FORM=form,STATUS='old')
        if(endianness_data.eq.1)open(unit_mcs+ishot,FILE=file_name,ACCESS=access,FORM=form,CONVERT='BIG_ENDIAN',STATUS='old')
endif

if(endianness_machine.eq.1)     then
        if(endianness_data.eq.0)open(unit_mcs+ishot,FILE=file_name,ACCESS=access,FORM=form,CONVERT='LITTLE_ENDIAN',STATUS='old')
        if(endianness_data.eq.1)open(unit_mcs+ishot,FILE=file_name,ACCESS=access,FORM=form,STATUS='old')
endif

enddo

END SUBROUTINE open_shot_mcs


subroutine Ricker_mcs(rank,iSS,raw_source)

use mod_parfile
use mod_data_arrays

implicit none

    integer :: i,j,k,rank,ishot,iSS,NumSou_SS
    integer :: SourceNum
    real    :: raw_source(nt_mcs)
    real :: pi,f0,wc,t0,tmp,S0(nt_mcs),t(nt_mcs)
    real :: n

	!!!!!------ Using a Ricker wavelet as impulse

        pi=3.14159265
        f0=freq_ricker
        wc=f0*2*pi
        t0=1.5/f0
        tmp=10./f0

        do i=1,nt_mcs
                t(i)=(i-1)*dt_mcs
        enddo

	S0=(1.-2.*(pi*f0*(t-t0))**2)*exp(-1.*(pi*f0*(t-t0))**2)

        do i=1,nt_mcs
                if(t(i).ge.tmp)then
                        S0(i)=0.
                endif
        enddo

        raw_source=0
        do i=1,nt_mcs
	        raw_source(i)=S0(i)
        enddo

end subroutine Ricker_mcs

subroutine get_bathymetry_model_mcs(rank)

USE mod_parfile
USE mod_data_arrays

implicit none

        integer :: i,j,ERR,rank

        real, allocatable :: xmodel(:),xmodel_(:),xbat_(:)
        character (len=500) :: file_name

        allocate(xmodel(nxmodel))
        allocate(xmodel_(NumShots_mcs+2),xbat_(NumShots_mcs+2))

        do i=2,NumShots_mcs+1
                xbat_(i)=pos_bat_mcs(i-1)
                xmodel_(i)=pos_shot_mcs(i-1)
        enddo

        xmodel_(1)=0;
        xmodel_(NumShots_mcs+2)=length_model;!!x

        xbat_(1)=xbat_(2);
        xbat_(NumShots_mcs+2)=xbat_(NumShots_mcs+1);!!y

!       do i=1,NumShots_mcs+2
!		    write(17,*)i,xmodel_(i),xbat_(i)
!	    enddo

        do i=1,nxmodel
            bat_model_mcs(i)=0
            xmodel(i)=(i-1)*dmodel
 	    !write(18,*)i,xmodel(i)!bat_model_mcs(i)
        enddo

        ERR=0
        call INTRPL(NumShots_mcs+2,xmodel_,xbat_,nxmodel,xmodel,bat_model_mcs,ERR)

        do i=1,nxmodel
                bat_model_grid_mcs(i)=1+ceiling(bat_model_mcs(i)/dmodel)
        enddo

        file_name=trim(adjustl(folder_output))// "bathymetry_model_meters_mcs.txt"
        open(12,file=file_name,status='unknown')
        do i=1,nxmodel
                write(12,*)xmodel(i),-bat_model_mcs(i)
        enddo
        close(12)

        deallocate(xmodel)
        deallocate(xmodel_,xbat_)

end subroutine get_bathymetry_model_mcs


subroutine get_pos_channels(rank)

USE mod_parfile
USE mod_data_arrays

implicit none
include 'mpif.h'

    integer :: numtasks,rank,ierr,errcode,status(MPI_STATUS_SIZE)
    integer :: i,j,k,ishot,irec
    INTEGER :: nh

    INTEGER*2 :: pos_scalco,scalco

    INTEGER*4 :: shotID,pos_read
    INTEGER*4 :: pos_byte_trace

    INTEGER*4 :: pos_offset,pos_sx,pos_sy
    INTEGER*4 :: sx,sy,offset

    call MPI_COMM_SIZE(MPI_COMM_WORLD,numtasks,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)

    ierr=0;

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    .SU FILE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! TO READ ShotID:

nh=size_su_header ! Size su header = 60*4 bytes
pos_read=byte_shotnumber_mcs
pos_offset = byte_offset

!------------------------------------------------------------------------

do ishot=1,NumShots_MCS

	READ(unit_mcs+ishot,pos=pos_read) shotID
	
	shotID_su_mcs(ishot)=shotID

        if(offset_header.eq.1)    then

	!    Loop trace by trace
	    do irec=1,NumRec_MCS

            pos_byte_trace=(irec-1)*(nh+nt_mcs)*4

            READ(unit_mcs+ishot,pos=pos_byte_trace+pos_offset) offset

            if(offset.eq.0)    then
                write(*,*)'ERROR: offset in header is 0. &
                set offset_header: 0 in file: ',trim(adjustl(par_file)),&
                ' or include geometry in SU-MCS file'
                call ascii_art(2)
                call MPI_Abort(MPI_COMM_WORLD, errcode, ierr)
                stop
            endif

            if(offset_unit.gt.0) then
                offset_su(irec,ishot)=abs(offset)*offset_unit
            endif

            if(offset_unit.lt.0) then
                offset_su(irec,ishot)=abs(offset)/abs(offset_unit)
            endif

  	    pos_rec(irec,ishot)=pos_shot_mcs(ishot)-offset_su(irec,ishot)
	    pos_rec_grid(irec,ishot)=1+ceiling(pos_rec(irec,ishot)/dmodel)

	    enddo   !irec

        endif !offset_header

        if(offset_header.eq.0)    then

	!    Loop trace by trace
	    do irec=1,NumRec_MCS

		pos_rec(irec,ishot)=pos_shot_mcs(ishot)-near_offset-(irec-1)*drec
!		write(112,*)ishot,irec,pos_rec(irec,ishot),pos_shot_mcs(ishot),near_offset,(irec-1)*drec
		pos_rec_grid(irec,ishot)=1+ceiling(pos_rec(irec,ishot)/dmodel)

	    enddo   !irec

        endif !offset_header


enddo	 !ishot

end subroutine get_pos_channels


function rec_ss(irec,isou,NumSou_SS)
USE mod_parfile
USE mod_data_arrays

implicit none
integer	:: isou,irec,idsou,NumSou_SS
integer :: rec_ss
real	:: dsou

	dsou=(NumSS_mcs)*dshots_mcs   !!space between shots in a SS
	idsou=floor(dsou)/drec

 	rec_ss=irec+(isou-1)*idsou

return
end function


function RecNum(irec_ss,isou,NumSou_SS)
USE mod_parfile
USE mod_data_arrays

implicit none
integer	:: isou,irec_ss,idsou,NumSou_SS
integer :: RecNum
real	:: dsou

	dsou=(NumSS_mcs)*dshots_mcs   !!space between shots in a SS
	idsou=floor(dsou)/drec

 	RecNum=irec_ss-(isou-1)*idsou

return
end function
