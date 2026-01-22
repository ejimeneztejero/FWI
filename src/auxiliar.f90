!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!	This file contains axiliary subroutines, mostly for reading/writting data
!!	Author:	Clara Estela Jimenez Tejero. 
!!	License: GNU General Public License v3.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE get_su_data(ni,ns,unit_su,vp_data)

USE mod_parfile
USE mod_data_arrays

IMPLICIT NONE

INTEGER :: nh
INTEGER :: i,j,k,ni,ns,unit_su
INTEGER(4) :: pos_r
REAL    :: vp_data(ns,ni)
REAL(4), ALLOCATABLE :: Data_tmp(:,:),sudata(:,:)

nh = size_su_header

allocate(Data_tmp(ns,ni))
allocate(sudata(ns+nh,ni))
Data_tmp=0.;sudata=0.;

do i=1,ni
        pos_r = 1 + (i-1)*(nh+ns)*4
        READ(unit_su, pos=pos_r)  sudata(1:ns+nh,i)
enddo

do i=1,ni
        do j=1,ns
                Data_tmp(j,i) = sudata(j+nh,i)
        enddo
enddo

vp_data=Data_tmp

!do i=1,ni
!do j=1,ns
!	write(111,*)i,j,vp_data(i,j)
!enddo
!enddo

deallocate(Data_tmp)
deallocate(sudata)

END SUBROUTINE get_su_data



subroutine file_iter(iter,name_iter)
USE mod_parfile
implicit none
integer :: iter
character(len=50) :: name_iter

if(niter.lt.10) name_iter=trim("_iter")
if(niter.lt.100)        then            
	if(iter.lt.10) name_iter=trim("_iter0")
	if(iter.ge.10) name_iter=trim("_iter")
endif
if(niter.ge.100.and.niter.lt.1000)      then            
	if(iter.lt.10) name_iter=trim("_iter00")
	if(iter.ge.10.and.iter.lt.100) name_iter=trim("_iter0")
	if(iter.ge.100) name_iter=trim("_iter")
endif

end subroutine file_iter

function SourceNum(iSS,k,NumSS)	!ojo, repasar que pasa en rounds>1
USE mod_parfile

integer :: SourceNum
integer :: iSS		!rank
integer :: k		!cada fuente
integer :: NumSS	!num ss

SourceNum=iSS+NumSS*(k-1)	!esta mal
!!tmp=floor(FirstSou+1.*(ShotNumThread-1)+THREADS*(i-1))	!OJO!!! para distribuir

! iOBS=SourceNum(iSS,k,NumSS_was)
! SourceNum(iSS,NumSou_SS,NumSS_mcs))

RETURN
end function


subroutine coding_fun(NumSou_SS,Coding)
use mod_parfile

implicit none
integer :: i,nn,NumSou_SS
real    :: n,Coding(NumSou_SS)
integer, allocatable:: seed(:)

Coding=1
nn=NumSou_SS
allocate(seed(nn))

if(NumSou_SS.gt.1) then

	if(seed_option.eq.5)    then
!	if(seed_option.eq.0)    then
		call random_seed(size=nn)
		seed = 123456789    ! putting arbitrary seed to all elements
		call random_seed(put=seed)
	endif

	do i=1,NumSou_SS

		call random_number(n)

		if (n .ge. 0.5) then

                        Coding(i)=1.0

                else

                        Coding(i)=-(1.0)!ojo, para promediar y evitar cross-talk noise, signo -

                endif

                write(*,*)"seed option, coding: ",seed_option,i,Coding(i)

        enddo

!       if(seed_option.eq.0)deallocate(seed)

endif

deallocate(seed)

end subroutine coding_fun

subroutine get_initial_models(rank)

  USE mod_parfile
  USE mod_data_arrays

  implicit none

  integer :: i,j,jwater,rank
  integer, allocatable :: bat_model_grid(:)
  real	:: xi,xj
  real, allocatable :: mi(:,:)
  character(len=500) :: file_name
  character(len=1000) :: file_init
  CHARACTER(len=50) :: access,form

  access = 'stream'
  form = 'unformatted'

  allocate(bat_model_grid(nxmodel))
  if(inv_was.eq.1)bat_model_grid=bat_model_grid_was
  if(inv_mcs.eq.1)bat_model_grid=bat_model_grid_mcs

!!! Vp
    Model_ini=water_velocity

    if(read_vp.eq.1)    then    !!ojo, modelo inicial deberia ya contener el agua

!!	ARREGLAR

        allocate(mi(nymodel,nxmodel))
        file_name=trim(folder_input_model) // trim(vp_file)

	if(vp_z.eq.1)	then
	        open(unit=10,file=file_name,status='old')
	        do j=1,nymodel
			read(10,*) (mi(j,i),i=1,nxmodel)
		enddo
	endif
	if(vp_xyz.eq.1)	then
	        open(unit=10,file=file_name,status='old')
                do i=1,nxmodel
                        do j=1,nymodel
                                read(10,*) xi,xj,mi(j,i)
                        enddo
                enddo
	endif

	if(vp_su.eq.1)	then
	        open(10,FILE=file_name,ACCESS=access,FORM=form,CONVERT='BIG_ENDIAN',STATUS='old')
		call get_su_data(nxmodel,nymodel,10,mi)
	endif

	Model_ini=mi

        deallocate(mi)
	close(10)

!     Si no hay modelo de velocidad entonces impone un gradiente
    else

	do i=1,nxmodel
		jwater=bat_model_grid(i)
!		write(76,*)i,jwater,bat_model_grid(i)
                do j=jwater,nymodel
                        Model_ini(j,i)=vpi+(j-jwater)*(vpf-vpi)/(nymodel-jwater)
                enddo
        enddo

    endif

!!  Density rho
    rho=1000.

!!  Bulk modulus, K:lambda
    lambda=rho*(Model_ini)**2.

    if(rank.eq.0)   then
	    file_init=trim(folder_VEL) // trim("initial_model")
	    call write_model(nxmodel,nymodel,file_init,Model_ini)
    endif

  deallocate(bat_model_grid)

end subroutine get_initial_models

subroutine get_f1f2(ifreq,f1,f2)!!repasar
use mod_parfile
use mod_data_arrays

implicit none

integer :: i,ifreq
real    :: f1,f2

do i=1,NumFreq
    Freq(i)=f_init+(i-1)*step_freq
enddo

if(BP_type.eq.1)   then    !!bandpass increases
    f1=f_init
    f2=Freq(ifreq)
endif

if(BP_type.eq.2)   then    !!same band pass
    f2=Freq(ifreq)
    f1=f2-step_freq
endif

if(BP_type.eq.3)   then

    if(ifreq.eq.2) then
        f1=f_init
        f2=Freq(2)
        misfit_type=3
    endif
    if(ifreq.gt.2) then
        f1=f_init
        f2=Freq(ifreq-1)
        misfit_type=1
    endif

endif

end subroutine get_f1f2

subroutine offset_sou_rec( off_sr, is, ir )
USE mod_parfile
USE mod_data_arrays
	
implicit none
integer, intent(in) :: is, ir
real 	:: ds,dr,off_sr

ds=pos_shot_was(ir)-pos_xobs(is)
dr=shot_depth_was-pos_zobs(is)

off_sr = sqrt( ds**2+ dr**2)

!write(11,*)is,ir,ds,dr,off_sr
!write(12,*)is,ir,pos_shot_was(ir),pos_xobs(is),pos_zobs(is)

end subroutine offset_sou_rec


subroutine is_not_cero(nt,data1d,isnotCero)
implicit none

    integer :: i,nt
    real :: data1d(nt)
    logical :: isnotCero

    isnotCero=.FALSE.

    ! Verificar si todas las componentes son cero
    do i = 1, nt
        if (data1d(i) /= 0.0) then
            isnotCero = .TRUE.
            exit ! Puedes salir del bucle en cuanto encuentres un valor no cero
        end if
    end do

end subroutine is_not_cero

subroutine write_model(ni,nj,file_name,model)	!!paralelizar

USE mod_parfile
USE mod_data_arrays

implicit none
integer :: i,j,k,ni,nj
real	:: xi,xj
real	:: model(nj,ni)
character(len=1000) :: Str,command
character(len=1000) :: file_name,file_name2
character(len=1000) :: file_vp_temp,file_vp_bin,file_vp_su


	if(save_z.ne.0)	then

		file_name2= trim(adjustl(file_name)) // ".z"

		open(unit=12,file=file_name2,status='unknown')

!                write(777,*)"nx.ny: ",ni,nj

		do j=1,nj
			write(12,'(20000(e12.5,2x))') (model(j,i),i=1,ni)
		enddo

		close(12)

	endif


	if(save_xyz.ne.0)       then

		file_name2= trim(adjustl(file_name)) // ".xyz"
		open(unit=12,file=file_name2,status='unknown')

	        do i=1,ni
			xi=(i-1)*dmodel
        	        do j=1,nj
				xj=(j-1)*dmodel
        	                write(12,*) xi,xj,model(j,i)
        	        enddo
        	enddo

		close(12)

	endif   


	if(save_su.ne.0)       then

                file_vp_temp=trim(adjustl(file_name)) // ".temp"

                open(unit=12,file=file_vp_temp,status='unknown')
                do i=1,ni
                do j=1,nj
                        write(12,*) model(j,i)
                enddo
                enddo
                close(12)

                file_vp_bin= trim(adjustl(file_name)) // ".bin"
                open(unit=14,file=file_vp_bin,status='unknown')

                file_vp_su= trim(adjustl(file_name)) // ".su"
                open(unit=16,file=file_vp_su,status='unknown')

                !write(*,*)"vp_temp: ",trim(adjustl(file_vp_temp))
                !write(*,*)"vp_bin: ",trim(adjustl(file_vp_bin))
                !write(*,*)"vp_su: ",trim(adjustl(file_vp_su))

                write(Str,*) nj

                command= "a2b < " // trim(adjustl(file_vp_temp)) // " n1=1 > " // trim(adjustl(file_vp_bin))
                !write(*,*)"command: ",trim(adjustl(command))
                call system(command)

                command= "suaddhead < " // trim(adjustl(file_vp_bin)) // &
                " ns=" // trim(adjustl(Str)) // " > " // trim(adjustl(file_vp_su))
                !write(*,*)"command: ",trim(adjustl(command))
                call system(command)

                close(14)
                close(16)

!                command="rm " // trim(adjustl(file_vp_bin)) // " " // trim(adjustl(file_vp_temp))
!                call system(command)

	endif


end subroutine write_model


subroutine write_source(rank,itype,nsou,nt,dt,file_name,model)	!!paralelizar

USE mod_parfile
USE mod_data_arrays

implicit none
integer :: i,j,k,rank
integer :: nsou,nt,nt_Data,itype
real	:: dt,dt_Data
real	:: model(nt,nsou)
real, allocatable :: model2(:,:)
character(len=1000) :: file_name,file_name2

	if(itype.eq.1)	then
		nt_Data=nt_was
		dt_Data=dt_was
		file_name2= trim(adjustl(file_name)) // "_was.dat"
	endif

	if(itype.eq.2)	then
		nt_Data=nt_mcs
		dt_Data=dt_mcs
		file_name2= trim(adjustl(file_name)) // "_mcs.dat"
	endif

	allocate(model2(nt_Data,nsou))

!!        !!!!!------ Filtering source with filter_function1D
        call interpolation_akima(model,dt,nt,model2,nsou,dt_Data,nt_Data)

	open(unit=12,file=file_name2,status='unknown')
	do j=1,nt_Data
		write(12,*) (model2(j,i),i=1,nsou)
	enddo
	close(12)

end subroutine write_source

SUBROUTINE write_new_data(rank,itype,iter,nt,dt,iSS,NumSou_SS,NumRec_SS,SG)
USE mod_parfile
USE mod_data_arrays

IMPLICIT NONE
integer :: i,j,k,iSS
integer :: num,inum,itype,unit_temp
integer	:: iter,rank
integer :: nh,nt,nt2,dt_head
integer :: NumSou_SS,NumRec_SS
integer(4) :: pos_r	
real :: dt,dt2
real :: SG(nt,NumRec_SS)
REAL(4), ALLOCATABLE :: sudata(:,:),shot(:,:),shot_new(:,:)
CHARACTER(len=50) :: access,form,iter_num,data_num,name_iter
CHARACTER(len=500) :: file_name,file_name2,file_temp
CHARACTER(len=500) :: command,command1,command2
CHARACTER(len=500) :: Str,Str2,Str3,Str_FWI

nh = size_su_header
unit_temp=12

ALLOCATE(shot(nt,NumRec_SS))
shot=0.;
access = 'stream'
form = 'unformatted'

write(iter_num,*) iter
write(data_num,*) iSS

call file_iter(iter,name_iter)

if(itype.eq.1)	then
	nt2=nt_was
	dt2=dt_was
	if(NumSou_SS.eq.1)then
	file_name2 = 'OBS' // trim(adjustl(data_num)) //  trim(adjustl(name_iter)) &
	// trim(adjustl(iter_num)) // '.su' 
	endif

	if(NumSou_SS.gt.1)then
	file_name2 = 'sOBS' // trim(adjustl(data_num)) //  trim(adjustl(name_iter)) // &
	trim(adjustl(iter_num)) // '.su' 
	endif

endif

if(itype.eq.2)	then
	nt2=nt_mcs
	dt2=dt_mcs
	file_name2 = 'shot' // trim(adjustl(data_num)) //  trim(adjustl(name_iter)) // &
	trim(adjustl(iter_num)) // '.su' 
endif

ALLOCATE(shot_new(nt2,NumRec_SS))
shot_new=0.

!! new data
do j=1,NumRec_SS
	do i=1,nt
		shot(i,j)=SG(i,j)
	enddo
enddo

!! interpolar a valores originales de field data (dt_Data,nt_Data)
call interpolation_akima(shot,dt,nt,shot_new,NumRec_SS,dt2,nt2)

!! write new data into new file
file_temp=trim(adjustl(folder_output)) // "temp" //  trim(adjustl(data_num)) // ".su"
open(unit_temp,FILE=file_temp,ACCESS=access,FORM=form,STATUS='unknown')

file_name = trim(adjustl(folder_DATA)) // trim(adjustl(file_name2))
num=unit_FWI+iSS
open(num,FILE=file_name,ACCESS=access,FORM=form,CONVERT='BIG_ENDIAN',STATUS='unknown')

do j=1,NumRec_SS
        pos_r=1+(j-1)*(nt2)*4
        WRITE(unit_temp, pos=pos_r+nh*4) shot_new(1:nt2,j) !! real traces
enddo

dt_head=dt2*1e6
write(Str,*) nt2
write(Str2,*) dt_head
write(Str3,*) NumRec_SS

command1= "suaddhead ns=" // trim(adjustl(Str))  // " < " // trim(adjustl(file_temp))
command2= "| sushw key=dt a=" // trim(adjustl(Str2)) // " > " // trim(adjustl(file_name))
command=trim(adjustl(command1)) // " " // trim(adjustl(command2))
!write(*,*) "command: ",trim(adjustl(command))
call system(command)

command="rm " // trim(adjustl(file_temp))
call system(command)

deallocate(shot)
close(unit_temp)
close(num)

end subroutine write_new_data

subroutine close_was_files()

USE mod_parfile
USE mod_data_arrays

implicit none
integer	:: k,iOBS

do iOBS=1,NumSou_WAS
!	write(*,*)"close obs, unit: ",unit_was+iOBS
	close(unit_was+iOBS)
enddo

end subroutine close_was_files

subroutine close_mcs_files()

USE mod_parfile
USE mod_data_arrays

implicit none
integer    :: k,ishot

do ishot=1,NumSou_MCS
!    write(*,*)"close obs, unit: ",unit_was+iOBS
    close(unit_mcs+ishot)
enddo

end subroutine close_mcs_files

subroutine get_FWI_file_unit(iter,iSS,num)

USE mod_parfile
implicit none
integer :: num,iSS
integer :: iter

	num=unit_FWI+10*iSS+iter

end subroutine get_FWI_file_unit


SUBROUTINE ascii_art(i)

USE mod_parfile

implicit none


integer :: i

if(i.eq.1)	then

write(*,*)
write(*,*)
write(*,*)                                                                      
write(*,*)"		Release (2025)			"
write(*,*)"		Author: Clara Estela Jiménez Tejero	"
write(*,*)"		email: ejimenez@icm.csic.es 		"
write(*,*)"		Barcelona Center for Subsurface Imaging "
write(*,*)"		Instituto de Ciencias Marinas (ICM-CSIC)"
write(*,*)
write(*,*)
write(*,*)"	                     |				"
write(*,*)"	                     |				"
write(*,*)"	            |        |				"
write(*,*)"	          |-|-|      |				"
write(*,*)"	            |        |				"
write(*,*)"	            | {O}    |				"
write(*,*)"	            '--|     |				"
write(*,*)"	              .|]_   |				"
write(*,*)"	        _.-=.' |     |				"	
write(*,*)"	       |    |  |]_   |				"
write(*,*)"	       |_.-='  |   __|__				"
write(*,*)"	        _.-='  |\   /|\				"
write(*,*)"	       |    |  |-'-'-'-'-.				"
write(*,*)"	       |_.-='  '========='				"
write(*,*)"	            `   |     |				"
write(*,*)"	             `. |    / \				"
write(*,*)"	               ||   /   \____.--=''''==--.._		"
write(*,*)"	               ||_.'--=='    |__  __  __  _.'	"
write(*,*)"	               ||  |    |    |\ ||  ||  || |                        ___	"	
write(*,*)"	  ____         ||__|____|____| \||__||__||_/    __________________/|   |	"
write(*,*)"	 |    |______  |===.---. .---.========''''=-._ |     |     |     / |   |	"
write(*,*)"	 |    ||     |\| |||   | |   |      '===' ||  \|_____|_____|____/__|___|	"
write(*,*)"	 |-.._||_____|_\___'---' '---'______....---===''======//=//////========|	"
write(*,*)"	 |--------------\------------------/-----------------//-//////---------/	"
write(*,*)"	 |               \                /                 // //////         /	"
write(*,*)"	 |                \______________/                 // //////         /	"
write(*,*)"	 |                                        _____===//=//////=========/	"
write(*,*)"	 |=================================================================/		"
write(*,*)"	  -----------------------------------------------------------------		"
write(*,*)"	``'-.,_,.-'``'-.,_,.='``'-.,_,.-'``'-.,_,.='````'-.,_,.-'``'-.,_,.='``'-.		"
write(*,*)"	``'-.,_,.-'``'-.,_,.='``'-.,_,.-'``'-.,_,.='````'-.,_,.-'``'-.,_,.='``'-.		"
write(*,*)"	``'-.,_,.-'``'-.,_,.='``'-.,_,.-'``'-.,_,.='````'-.,_,.-'``'-.,_,.='``'-.		"
write(*,*)
write(*,*)	
write(*,*)"				         ______		"
write(*,*)"				        /      \ 	"	
write(*,*)"				       /        \ 	"
write(*,*)"				       |        |	"
write(*,*)"				    )  o        o   (	"
write(*,*)"				   (    \      /    )	"
write(*,*)"				  _ \___/||||||\___/ _	"
write(*,*)"				   \____/ |||| \____/ `	"
write(*,*)"				   ,-.___/ || \__,-._	"
write(*,*)"				  /    ___/  \__	"
write(*,*)"				     _/         `--	"
write(*,*)
write(*,*)
write(*,*)"	``'-.,_,.-'``'-.,_,.='``'-.,_,.-'``'-.,_,.='````'-.,_,.-'``'-.,_,.='``'-.		"
write(*,*)"	``'-.,_,.-'``'-.,_,.='``'-.,_,.-'``'-.,_,.='````'-.,_,.-'``'-.,_,.='``'-.		"
write(*,*)"	``'-.,_,.-'``'-.,_,.='``'-.,_,.-'``'-.,_,.='````'-.,_,.-'``'-.,_,.='``'-.		"
write(*,*)
write(*,*)


endif

if(i.eq.2)	then

write(*,*)
write(*,*)
write(*,*)'Fix your inputs please'
write(*,*)'                                |     |'
write(*,*)'                                \\_V_//'
write(*,*)'                                \/=|=\/'
write(*,*)'                                 [=v=]'
write(*,*)'                               __\___/_____'
write(*,*)'                              /..[  _____  ]'
write(*,*)'                             /_  [ [  M /] ]'
write(*,*)'                            /../.[ [ M /@] ]'
write(*,*)'                           <-->[_[ [M /@/] ]'
write(*,*)'                          /../ [.[ [ /@/ ] ]'
write(*,*)'     _________________]\ /__/  [_[ [/@/ C] ]'
write(*,*)'    <_________________>>0---]  [=\ \@/ C / /'
write(*,*)'       ___      ___   ]/000o   /__\ \ C / /'
write(*,*)'          \    /              /....\ \_/ /'
write(*,*)'       ....\||/....           [___/=\___/'
write(*,*)'      .    .  .    .          [...] [...]'
write(*,*)'     .      ..      .         [___/ \___]'
write(*,*)'     .    0 .. 0    .         <---> <--->'
write(*,*)'  /\/\.    .  .    ./\/\      [..]   [..]'
write(*,*)' / / / .../|  |\... \ \ \    _[__]   [__]_'
write(*,*)'/ / /       \/       \ \ \  [____>   <____]'
write(*,*)
write(*,*)

endif

if(i.eq.3)	then

	write(*,*)
	write(*,*)
	write(*,*)
	write(*,*)
	write(*,*)"	 _______  ___   __    _  ___   _______  __   __  _______  ______  	"
	write(*,*)"	|       ||   | |  |  | ||   | |       ||  | |  ||       ||      | 	"
	write(*,*)"	|    ___||   | |   |_| ||   | |  _____||  |_|  ||    ___||  _    |	"
	write(*,*)"	|   |___ |   | |       ||   | | |_____ |       ||   |___ | | |   |	"
	write(*,*)"	|    ___||   | |  _    ||   | |_____  ||       ||    ___|| |_|   |	"
	write(*,*)"	|   |    |   | | | |   ||   |  _____| ||   _   ||   |___ |       |	"
	write(*,*)"	|___|    |___| |_|  |__||___| |_______||__| |__||_______||______|	"
	write(*,*)
	write(*,*)

endif

if(i.eq.4)	then

write(*,*)""
write(*,*)""
write(*,*)"*** WARNING: numtasks do not need to be greater than variable NumOBS"
write(*,*)"***          You are wasting energy. Please, next time you run a job, set in you MPI execution line, numtasks=",NumOBS
write(*,*)""

endif

END SUBROUTINE ascii_art
