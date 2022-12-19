!=======================================================================================================
!
!          Subroutine to initialize netcdf files for REcoM diagnostics
!
!=======================================================================================================
subroutine REcoM_init_output(yearnew, yearold, cyearnew)
  use recom_ocean_settings
  use recom_config
  use REcoM_declarations
  implicit none
  
#include "netcdf.inc" 

  integer, INTENT(IN)       :: yearnew, yearold
  character(4), INTENT(IN)  :: cyearnew
  integer                   :: status, ncid, dimid_rec, j
  integer                   :: dimid_2d, dimid_3d, dimids(2)
  integer                   :: time_varid, iter_varid, benthos_varid(4), hplus_varid
  character(2000)           :: longname, filename
  
  if(yearnew==yearold) return
  if (recom_mype/=0) return
  
!--  Initializing output files for snap-shots  -------------------------------------------------

  filename=trim(recom_ResultPath)//recom_runid//'.'//cyearnew//'.bio.nc'
  
  ! create a file
  status = nf_create(filename, nf_clobber, ncid)
  if (status.ne.nf_noerr) call handle_err(status)
  
  ! Define the dimensions
  status = nf_def_dim(ncid, 'nodes_2d', recom_nod2d, dimid_2d)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_dim(ncid, 'nodes_3d', recom_nod3d, dimid_3d)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_dim(ncid, 'time', NF_UNLIMITED, dimid_rec)
  if (status .ne. nf_noerr) call handle_err(status)
  ! Define the time and iteration variables
  status = nf_def_var(ncid, 'time', NF_DOUBLE, 1, dimid_rec, time_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'iter', NF_INT, 1, dimid_rec, iter_varid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! Define the netCDF variables for 2D fields.
  ! In Fortran, the unlimited dimension must come
  ! last on the list of dimids.
  dimids(1) = dimid_2d
  dimids(2) = dimid_rec
  
  ! Benthos
  status = nf_def_var(ncid, 'BenN', NF_DOUBLE, 2, dimids, benthos_varid(1))
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'BenC', NF_DOUBLE, 2, dimids, benthos_varid(2))
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'BenSi', NF_DOUBLE, 2, dimids, benthos_varid(3))
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'BenCalc', NF_DOUBLE, 2, dimids, benthos_varid(4))
  if (status .ne. nf_noerr) call handle_err(status)
  
  status = nf_def_var(ncid, 'Hplus', NF_DOUBLE, 2, dimids, Hplus_varid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! Define the netCDF variables for 3D fields
  
  !dimids(1) = dimid_3d
  !dimids(2) = dimid_rec
  
  ! Benthos
  longname='Benthos Nitrogen'
  status = nf_put_att_text(ncid, benthos_varid(1), 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, benthos_varid(1), 'units', 7, 'mmol/m3')
  if (status .ne. nf_noerr) call handle_err(status)

  longname='Benthos Carbon'
  status = nf_put_att_text(ncid, benthos_varid(2), 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, benthos_varid(2), 'units', 7, 'mmol/m3')
  if (status .ne. nf_noerr) call handle_err(status)

  longname='Benthos Silicate'
  status = nf_put_att_text(ncid, benthos_varid(3), 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, benthos_varid(3), 'units', 7, 'mmol/m3')
  if (status .ne. nf_noerr) call handle_err(status)

  longname='Benthos Calcite'
  status = nf_put_att_text(ncid, benthos_varid(4), 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, benthos_varid(4), 'units', 7, 'mmol/m3')
  if (status .ne. nf_noerr) call handle_err(status)
  
  longname='Conc. of H-plus ions in the surface water'
  status = nf_put_att_text(ncid, Hplus_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, Hplus_varid, 'units', 6, 'mol/kg')
  if (status .ne. nf_noerr) call handle_err(status)

  status = nf_enddef(ncid)
  if (status .ne. nf_noerr) call handle_err(status)

  status=nf_close(ncid)
  if (status .ne. nf_noerr) call handle_err(status)

end subroutine REcoM_init_output

!=======================================================================================================
!
!          Subroutine that writes output (benthos, CO2, Fe and diagnostics) from REcoM
!
!=======================================================================================================
!subroutine REcoM_write_snapshots(cyearnew,recom_istep)
!  use recom_ocean_settings
!  use REcoM_GloVar
 ! use REcoM_LocVar
 ! use REcoM_variables
! use recom_config
!  use REcoM_declarations
 ! implicit none
  
!#include "netcdf.inc"

!  character(4), INTENT(IN)  :: cyearnew
!  integer, INTENT(IN)       :: recom_istep

!  integer                   :: time_varid, iter_varid, benthos_varid(4), Hplus_varid
!  integer                   :: status, ncid, j
!  integer                   :: start(2), count(2), n3
!  real(kind=8)              :: sec_in_year
!  character(100)            :: filename
!  real(kind=8), allocatable :: aux2(:) !, aux3(:)
  
!  allocate(aux2(recom_nod2D))
!!  allocate(aux3(recom_nod3D))

!  if (recom_mype==0) then
!    sec_in_year=recom_dt*recom_istep
    
!    ! open files
!     filename=trim(recom_ResultPath)//recom_runid//'.'//cyearnew//'.bio.nc'
!     status = nf_open(filename, nf_write, ncid)
!     if (status .ne. nf_noerr) call handle_err(status)
     
     ! inquire variable id
!     status=nf_inq_varid(ncid, 'time', time_varid)
!     if (status .ne. nf_noerr) call handle_err(status)
!     status=nf_inq_varid(ncid, 'iter', iter_varid)
!     if (status .ne. nf_noerr) call handle_err(status)
     
!     status=nf_inq_varid(ncid, 'BenN', benthos_varid(1))
!     if (status .ne. nf_noerr) call handle_err(status)
!     status=nf_inq_varid(ncid, 'BenC', benthos_varid(2))
!     if (status .ne. nf_noerr) call handle_err(status)
!     status=nf_inq_varid(ncid, 'BenSi', benthos_varid(3))
!     if (status .ne. nf_noerr) call handle_err(status)
!     status=nf_inq_varid(ncid, 'BenCalc', benthos_varid(4))
!     if (status .ne. nf_noerr) call handle_err(status)
     
!     status=nf_inq_varid(ncid, 'Hplus', Hplus_varid)
!     if (status .ne. nf_noerr) call handle_err(status)

! write variables

     ! time and iteration
!     status=nf_put_vara_double(ncid, time_varid, recom_save_count_restart, 1, sec_in_year)
!     if (status .ne. nf_noerr) call handle_err(status)
!     status=nf_put_vara_int(ncid, iter_varid, recom_save_count_restart, 1, recom_istep)
!     if (status .ne. nf_noerr) call handle_err(status)
     
!  end if    !! recom_mype==0 
  
!! 2d fields
!  ! Benthos
!  call broadcast2D(Benthos(1,:),aux2)
!  if(recom_mype==0) then            
!    start=(/1,recom_save_count_restart/)
!    count=(/recom_nod2d, 1/)
!    status=nf_put_vara_double(ncid, benthos_varid(1), start, count, aux2) 
!    if (status .ne. nf_noerr) call handle_err(status)
!  end if
    
!  call broadcast2D(Benthos(2,:),aux2)
!  if(recom_mype==0) then            
!    start=(/1,recom_save_count_restart/)
!    count=(/recom_nod2d, 1/)
!    status=nf_put_vara_double(ncid, benthos_varid(2), start, count, aux2) 
!    if (status .ne. nf_noerr) call handle_err(status)
!  end if
    
!  call broadcast2D(Benthos(3,:),aux2)
!  if(recom_mype==0) then            
!    start=(/1,recom_save_count_restart/)
!    count=(/recom_nod2d, 1/)
!    status=nf_put_vara_double(ncid, benthos_varid(3), start, count, aux2) 
!    if (status .ne. nf_noerr) call handle_err(status)
!  end if
  
!  call broadcast2D(Benthos(4,:),aux2)
!  if(recom_mype==0) then            
!    start=(/1,recom_save_count_restart/)
!    count=(/recom_nod2d, 1/)
!    status=nf_put_vara_double(ncid, benthos_varid(4), start, count, aux2) 
!    if (status .ne. nf_noerr) call handle_err(status)
!  end if
  
!  call broadcast2D(GloHplus,aux2)
!  if(recom_mype==0) then            
!    start=(/1,recom_save_count_restart/)
!    count=(/recom_nod2d, 1/)
!    status=nf_put_vara_double(ncid, Hplus_varid, start, count, aux2) 
!    if (status .ne. nf_noerr) call handle_err(status)
!  end if

!  deallocate(aux2)
!!  deallocate(aux2,aux3)
  
!  if(recom_mype==0) then
!    status=nf_close(ncid)
!    if (status .ne. nf_noerr) call handle_err(status)
!  endif
         
!end subroutine REcoM_write_snapshots

!====================================================================
! Fortran binary (unformatted) restarts:
! ====================================================================
subroutine save_recom_binary(tracer, cyearnew)
  use recom_ocean_settings
  use recom_config
  implicit none
   character(4), INTENT(IN)  :: cyearnew
   integer                   :: fileID
   integer                   :: n, phys_num
   real(kind=8), allocatable :: array_3d(:)
   character(2000)            :: filename
   real(kind=8), dimension(recom_todim_nod3d,recom_num_tracer), INTENT(IN) :: tracer
   
   allocate(array_3d(recom_nod3D))
   ! CN: add variable to determine how many physical tracers exist to determine                      
   ! the index of the first bgc tracer                                                               
   phys_num = recom_num_tracer - bgc_num
   !
   ! Write restart
   if (recom_mype==0) then
      write(*,*) 'Writing recom restart'
      fileID=125
      filename = trim(recom_ResultPath)//'recom_restart'//cyearnew//'.bin'
      open(fileID,file=filename, form='unformatted')
   end if
   ! CN: make flexible depending on how many physical tracers exist                                  
   do n=phys_num+1, recom_num_tracer
!   do n=3, recom_num_tracer
      call broadcast3D(tracer(:,n), array_3d)
      if (recom_mype==0) then
         write(fileID) array_3d
      end if
   end do
   if (recom_mype==0) then
      close(fileID)
   end if
   deallocate(array_3d)
  end subroutine save_recom_binary


