! Benthos is declared and initialized.
! BGC tracers are initialized
! This routine is called from 'oce_input' and is therefore only called when fesom is doing a restart run.
!========================================================================================================
  
subroutine recom_init(tracer, yearnew, cyearnew, yearold, cyearold)
  use REcoM_declarations
  use REcoM_GloVar
  use REcoM_locVar
  use recom_config
  use recom_ocean_settings
  implicit none 
  
  integer, INTENT(IN) :: yearnew, yearold
  character(4), INTENT(IN) :: cyearnew, cyearold

#include "netcdf.inc"
  real(kind=8), dimension(recom_todim_nod3d,recom_num_tracer), INTENT(INOUT) :: tracer
 
  integer                   :: status, ncid, dimid_rec, nrec, varid
  integer                   :: tra_varid(bgc_num), benthos_varid(4), hplus_varid
  integer                   :: din_varid, dsi_varid, dfe_varid, do2_varid, alk_varid, dic_varid, zoo2n_varid,zoo2c_varid
  Integer                   :: istart(2), icount(2), size2D, size3D     ! Number of 2D and 3D nodes on current processor
  Integer                   :: i, j
  character(2000)                    :: filename
  character(2)                       :: trind
  character(4)                       :: tr_name
  real(kind=8), allocatable          :: aux2(:), aux3(:)
  real(kind=8), allocatable          :: ncdata(:)
  character(2000)                    :: CO2filename
  character(20)                      :: CO2vari
  integer                            :: CO2start, CO2count

  if (.not. use_REcoM) return
  
  allocate(aux2(recom_nod2D),aux3(recom_nod3D))
  size2D = recom_myDim_nod2D+recom_eDim_nod2D
  size3D = recom_myDim_nod3D+recom_eDim_nod3D
  
!--------------------------------------------------------------------------------------------------------
! Allocating and initializing atm dep of iron and arrays for diagnostics
! When felimit is not used, arrays are set to zero
  
  allocate(GloFeDust(size2D))
  allocate(AtmFeInput(size2D))
  GloFeDust       = 0.d0
  AtmFeInput      = 0.d0

  if (recom_mype==0) then
    if (UseDustClimAlbani) then
      write(*,*),'Monthly dust climatology is used (from Albani)'
    else
      if (yearnew .LT. 1979) then
        write(*,*),'Monthly dust climatology (Mahowald) is used as interannually varying files do not exist for current year'
      else
        if (UseDustClim .eq. .true.) then
          write(*,*),'Monthly dust climatology is used (from Mahowald)'
        else
       	  write(*,*),'Interannually varying monthly dust is used (from Mahowald)'
        end if ! DustClim
      end if   ! year < 1979
    end if     ! Albani 
  end if       ! recom_mype=0

! Also need to be allocated when nitrogen sources and sinks not are used
  allocate(GloNDust(size2D))
  allocate(AtmNInput(size2D))
  GloNDust         = 0.d0
  AtmNInput        = 0.d0

  allocate(GloHplus(size2D))
  allocate(GloPCO2surf(size2D))
  allocate(GlodPCO2surf(size2D))
  allocate(GloCO2flux(size2D))
#ifdef __cpl_echam 
  allocate(GloCO2flux_noicemask(size2D))
#endif
  allocate(GloCO2flux_seaicemask(size2D))
  allocate(GloO2flux_seaicemask(size2D))
  allocate(PistonVelocity(size2D))
  allocate(alphaCO2(size2D))



  GloHplus         = exp(-8.d0 * log(10.d0))          !=10**(-8)
  GloPCO2surf      = 0.d0
  GlodPCO2surf      = 0.d0
  GloCO2flux       = 0.d0
#ifdef __cpl_echam 
  GloCO2flux_noicemask       = 0.d0
#endif
  GloCO2flux_seaicemask = 0.0d0 
  GloO2flux_seaicemask = 0.0d0 
  PistonVelocity     = 0.d0 ! Piston velocity                                                    
  alphaCO2           = 0.d0 ! CO2 solubility                                                     

  AtmCO2   = 0.d0
  Hplus    = 0.d0
  pco2surf = 0.d0
  dflux    = 0.d0
  co2flux_seaicemask = 0.d0
  o2flux_seaicemask = 0.d0
  dpco2surf= 0.d0
  co2      = 0.d0                  

  if (Diagnostics) then
    allocate(diags2D(12,size2D))          
    diags2D(:,:)      = 0.d0
    allocate(diags3D(size3D,diags3d_num))
    diags3D(:,:)      = 0.d0
  end if  

  allocate(CO23D(size3D))             
  CO23D(:)          = 0.d0

  allocate(pH3D(size3D))          
  pH3D(:)           = 0.d0

  allocate(pCO23D(size3D))         
  pCO23D(:)         = 0.d0

  allocate(HCO33D(size3D))         
  HCO33D(:)         = 0.d0

  allocate(CO33D(size3D))          
  CO33D(:)          = 0.d0

  allocate(OmegaC3D(size3D))       
  OmegaC3D(:)       = 0.d0

  allocate(kspc3D(size3D))         
  kspc3D(:)         = 0.d0

  allocate(rhoSW3D(size3D))          
  rhoSW3D(:)        = 0.d0

  allocate(Nutlim_phy3D(size3D))    
  Nutlim_phy3D(:)   = 0.d0

  allocate(Nutlim_dia3D(size3D))     
  Nutlim_dia3D(:)   = 0.d0

  allocate(Nutlim_cocco3D(size3D))   
  Nutlim_cocco3D(:) = 0.d0

  allocate(Tlim_arr3D(size3D))       
  Tlim_arr3D(:)     = 0.d0

  allocate(Tlim_cocco3D(size3D))     
  Tlim_cocco3D(:)   = 0.d0

  allocate(Llim_phy3D(size3D))       
  Llim_phy3D(:)     = 0.d0

  allocate(Llim_dia3D(size3D))      
  Llim_dia3D(:)     = 0.d0

  allocate(Llim_cocco3D(size3D))   
  Llim_cocco3D(:)   = 0.d0

  allocate(CO2lim_phy3D(size3D))  
  CO2lim_phy3D(:)   = 0.d0

  allocate(CO2lim_dia3D(size3D))  
  CO2lim_dia3D(:)   = 0.d0

  allocate(CO2lim_cocco3D(size3D))   
  CO2lim_cocco3D(:) = 0.d0

  allocate(PR_phy3D(size3D))      
  PR_phy3D(:)       = 0.d0

  allocate(PR_dia3D(size3D))        
  PR_dia3D(:)       = 0.d0

  allocate(PR_cocco3D(size3D))     
  PR_cocco3D(:)     = 0.d0

  allocate(Cal_Tlim3D(size3D))     
  Cal_Tlim3D(:)     = 0.d0

  allocate(Cal_CO2lim3D(size3D))    
  Cal_CO2lim3D(:)   = 0.d0

  allocate(Cal_Nlim3D(size3D))    
  Cal_Nlim3D(:)     = 0.d0

  allocate(Cal_pure3D(size3D))        
  Cal_pure3D(:)     = 0.d0

  allocate(PAR3D(size3D))
  PAR3D(:)          = 0.d0

  allocate(DenitBen(size2D)) 
  DenitBen(:)       = 0.d0

  allocate(Benthos(4,Size2D))



  allocate(index_recom_tracer(bgc_num)) 
  do j=1, bgc_num
    write(trind,'(i2.2)') j
    tr_name='tr'//trind
    do i=1, recom_num_tracer ! all tracers, incl T and S
      if(recom_prog_tracer_name(i) == tr_name) index_recom_tracer(j)=i
    end do
  end do

!---------------------------------------------------------------------------------------------------------
! Initialization of tracers

  if ((REcoM_restart .eq. .true.) .and. (recom_binary_init .eq. .false.)) then
    ! open file
    filename=trim(recom_ResultPath)//recom_runid//'.'//cyearold//'.oce.nc'
    status = nf_open(filename, nf_nowrite, ncid)
    if (status .ne. nf_noerr) call handle_err(status)
    
    ! inquire variable id
    do j=1,bgc_num
      write(trind,'(i2.2)') j
      status=nf_inq_varid(ncid, 'tr'//trind, tra_varid(j))
      if (status .ne. nf_noerr) call handle_err(status)	
    end do
    
    ! read variables
  
    ! which record to read
    if(recom_restartflag=='last') then
      status = nf_inq_dimid(ncid, 'time', dimid_rec)
      if(status .ne. nf_noerr) call handle_err(status)
      status = nf_inq_dimlen(ncid, dimid_rec, nrec)
      if(status .ne. nf_noerr) call handle_err(status)
    else
      read(recom_restartflag,'(i4)') nrec
    end if
    
    ! 3d fields
    istart=(/1,nrec/)
    icount=(/recom_nod3d, 1/)
    do j=1,bgc_num
      status=nf_get_vara_double(ncid, tra_varid(j), istart, icount, aux3) 
      if (status .ne. nf_noerr) call handle_err(status)
      tracer(:,index_recom_tracer(j))=aux3(recom_myList_nod3D)
    end do
        
    status=nf_close(ncid)
    if (status .ne. nf_noerr) call handle_err(status)
    if(recom_mype==0) write(*,*),'Tracers have been initialized as restart from netcdf files'
    
!-------------------------------------------------------------------------------------------------
! Initialization of benthos from file -.bio.diag.nc

    ! open file
    filename=trim(recom_ResultPath)//recom_runid//'.'//cyearold//'.bio.nc'
    status = nf_open(filename, nf_nowrite, ncid)
    if (status .ne. nf_noerr) call handle_err(status)    
    
    ! inquire variable id
    status=nf_inq_varid(ncid, 'BenN', benthos_varid(1))
	if (status .ne. nf_noerr) call handle_err(status)
	    
	status=nf_inq_varid(ncid, 'BenC', benthos_varid(2))
    if (status .ne. nf_noerr) call handle_err(status)
         
    status=nf_inq_varid(ncid, 'BenSi', benthos_varid(3))
    if (status .ne. nf_noerr) call handle_err(status)
         
    status=nf_inq_varid(ncid, 'BenCalc', benthos_varid(4))
    if (status .ne. nf_noerr) call handle_err(status)

    status=nf_inq_varid(ncid, 'Hplus', hplus_varid)
    if (status .ne. nf_noerr) call handle_err(status)
    
    ! read variables
  
    ! which record to read
    if(recom_restartflag=='last') then
      status = nf_inq_dimid(ncid, 'time', dimid_rec)
      if(status .ne. nf_noerr) call handle_err(status)
      status = nf_inq_dimlen(ncid, dimid_rec, nrec)
      if(status .ne. nf_noerr) call handle_err(status)
    else
      read(recom_restartflag,'(i4)') nrec
    end if
  
    ! 2d fields
    istart=(/1,nrec/)
    icount=(/recom_nod2d, 1/)
    do j=1,4
          status=nf_get_vara_double(ncid, benthos_varid(j), istart, icount, aux2) 
          if (status .ne. nf_noerr) call handle_err(status)
          Benthos(j,:)=aux2(recom_myList_nod2D)
    end do

    status=nf_get_vara_double(ncid, hplus_varid, istart, icount, aux2) 
    if (status .ne. nf_noerr) call handle_err(status)
    GloHplus(:)=aux2(recom_myList_nod2D)
    
    status=nf_close(ncid)
    if (status .ne. nf_noerr) call handle_err(status)
    if(recom_mype==0) write(*,*),'Tracers have been initialized for REcoM using values from last run'



  else if ((REcoM_restart .eq. .true.) .and. (recom_binary_init .eq. .true.)) then


    call read_recom_binary(tracer, cyearold)
    if(recom_mype==0) write(*,*),'Tracers have been initialized as restart from binary file'

  else

    ! Initializing REcoM to do spin-up
    ! DIN, DIC, alk, Si and Fe are initialized from databases (not as spin-up)
    
    ! DIN is initialized
    if(recom_mype==0) write(*,*),'Reading REcoM2 DIN for restart'
    filename=trim(REcoMDataPath)//'InitDIN.nc'
    status = nf_open(filename, nf_nowrite, ncid)
    if (status .ne. nf_noerr) call handle_err(status)

    ! inquire variable id
    status=nf_inq_varid(ncid, 'DIN', din_varid)
    if (status .ne. nf_noerr) call handle_err(status)
    
    ! read variables
    istart=(/1,1/)
    icount=(/recom_nod3d, 1/)
    status=nf_get_vara_double(ncid, din_varid, istart, icount, aux3) 
    if (status .ne. nf_noerr) call handle_err(status)
    tracer(:,index_recom_tracer(idin))=aux3(recom_myList_nod3D)

    status=nf_close(ncid)

    ! DIC is initialized
    if (REcoM_PI .eq. .false.) then ! Use present day DIC for spin-up


      if(recom_mype==0) write(*,*),'Reading REcoM2 DIC for restart'
      filename=trim(REcoMDataPath)//'InitDIC.nc'
      status = nf_open(filename, nf_nowrite, ncid)
      if (status .ne. nf_noerr) call handle_err(status)

      ! inquire variable id
      status=nf_inq_varid(ncid, 'DIC', dic_varid)
      if (status .ne. nf_noerr) call handle_err(status)
    
      ! read variables
      istart=(/1,1/)
      icount=(/recom_nod3d, 1/)
      status=nf_get_vara_double(ncid, dic_varid, istart, icount, aux3) 
      if (status .ne. nf_noerr) call handle_err(status)
      tracer(:,index_recom_tracer(idic))=aux3(recom_myList_nod3D)

      status=nf_close(ncid)

    else ! Initialize DIC with pre-industrial concentrations

      if(recom_mype==0) write(*,*),'Reading REcoM2 preindustrial DIC for restart'
      filename=trim(REcoMDataPath)//'InitDIC_PI.nc'
      status = nf_open(filename, nf_nowrite, ncid)
      if (status .ne. nf_noerr) call handle_err(status)

      ! inquire variable id
      status=nf_inq_varid(ncid, 'DIC', dic_varid)
      if (status .ne. nf_noerr) call handle_err(status)
    
      ! read variables
      istart=(/1,1/)
      icount=(/recom_nod3d, 1/)
      status=nf_get_vara_double(ncid, dic_varid, istart, icount, aux3) 
      if (status .ne. nf_noerr) call handle_err(status)
      tracer(:,index_recom_tracer(idic))=aux3(recom_myList_nod3D)

      status=nf_close(ncid)


  end if ! Use pre-industrial dic for spin-up

    ! Alkalinity is initialized
    if(recom_mype==0) write(*,*),'Reading REcoM2 alkalinity for restart'
    filename=trim(REcoMDataPath)//'InitAlk.nc'

    status = nf_open(filename, nf_nowrite, ncid)
    if (status .ne. nf_noerr) call handle_err(status)

    status=nf_inq_varid(ncid, 'Alk', alk_varid)
    if (status .ne. nf_noerr) call handle_err(status)

    ! read variables
    istart=(/1,1/)
    icount=(/recom_nod3d, 1/)
    status=nf_get_vara_double(ncid, alk_varid, istart, icount, aux3) 
    if (status .ne. nf_noerr) call handle_err(status)
    tracer(:,index_recom_tracer(ialk))=aux3(recom_myList_nod3D)

    status=nf_close(ncid)

    ! DSi is initialized
    if(recom_mype==0) write(*,*),'Reading REcoM2 silicon for restart'
    filename=trim(REcoMDataPath)//'InitDSi.nc'

    status = nf_open(filename, nf_nowrite, ncid)
    if (status .ne. nf_noerr) call handle_err(status)
    
    status=nf_inq_varid(ncid, 'DSi', dsi_varid)
    if (status .ne. nf_noerr) call handle_err(status)

    ! read variables
    istart=(/1,1/)
    icount=(/recom_nod3d, 1/)
    status=nf_get_vara_double(ncid, dsi_varid, istart, icount, aux3) 
    if (status .ne. nf_noerr) call handle_err(status)
    tracer(:,index_recom_tracer(isi))=aux3(recom_myList_nod3D)

    status=nf_close(ncid)

    ! DFe is initialized
    if(recom_mype==0) write(*,*),'Reading REcoM2 iron for restart'
    filename=trim(REcoMDataPath)//'InitDFe.nc'

    status = nf_open(filename, nf_nowrite, ncid)
    if (status .ne. nf_noerr) call handle_err(status)
    
    status=nf_inq_varid(ncid, 'DFe', dfe_varid)
    if (status .ne. nf_noerr) call handle_err(status)	
    
    ! read variables
    istart=(/1,1/)
    icount=(/recom_nod3d, 1/)
    status=nf_get_vara_double(ncid, dfe_varid, istart, icount, aux3) 
    if (status .ne. nf_noerr) call handle_err(status)
    tracer(:,index_recom_tracer(ife))=aux3(recom_myList_nod3D) 

    status=nf_close(ncid)

    ! oxygen is initialized
    if(recom_mype==0) write(*,*),'Reading REcoM2 oxygen for restart'
    filename=trim(REcoMDataPath)//'InitO2.nc' ! micromoles_per_kilogram 

    status = nf_open(filename, nf_nowrite, ncid)
    if (status .ne. nf_noerr) call handle_err(status)
    
    status=nf_inq_varid(ncid, 'DO2', do2_varid)
    if (status .ne. nf_noerr) call handle_err(status)

    ! read variables
    istart=(/1,1/)
    icount=(/recom_nod3d, 1/)
    status=nf_get_vara_double(ncid, do2_varid, istart, icount, aux3) 
    if (status .ne. nf_noerr) call handle_err(status)
    tracer(:,index_recom_tracer(ioxy))=aux3(recom_myList_nod3D)

    status=nf_close(ncid)
  
    if (REcoM_Second_Zoo .and. zoo2_initial_field) then                                                                    
     !Zoo2N initialized                                    
    
     if(recom_mype==0) write(*,*),'Reading REcoM2 zoo2n for restart'
     filename=trim(REcoMDataPath)//'Init_zoo2n.nc'

     status = nf_open(filename, nf_nowrite, ncid)
     if (status .ne. nf_noerr) call handle_err(status)

     status=nf_inq_varid(ncid, 'zoo2n', zoo2n_varid)
     if (status .ne. nf_noerr) call handle_err(status)

     ! read variables                                                                                                                                                                                        
     istart=(/1,1/)
     icount=(/recom_nod3d, 1/)
     status=nf_get_vara_double(ncid, zoo2n_varid, istart, icount, aux3)
     if (status .ne. nf_noerr) call handle_err(status)
     tracer(:,index_recom_tracer(izoo2n))=aux3(recom_myList_nod3D)

     status=nf_close(ncid)

     ! Zoo2C is initialized                                                                                            

     if(recom_mype==0) write(*,*),'Reading REcoM2 zoo2c for restart'
     filename=trim(REcoMDataPath)//'Init_zoo2c.nc'

     status = nf_open(filename, nf_nowrite, ncid)
     if (status .ne. nf_noerr) call handle_err(status)

     status=nf_inq_varid(ncid, 'zoo2c', zoo2c_varid)
     if (status .ne. nf_noerr) call handle_err(status)

     ! read variables                                                                                              

     istart=(/1,1/)
     icount=(/recom_nod3d, 1/)
     status=nf_get_vara_double(ncid, zoo2c_varid, istart, icount, aux3)
     if (status .ne. nf_noerr) call handle_err(status)
     tracer(:,index_recom_tracer(izoo2c))=aux3(recom_myList_nod3D)

     status=nf_close(ncid)

    endif



! change to fesom tracer counter, not hard coded counters:
    tracer(:,index_recom_tracer(iphyn))  = tiny            ! tracer 6  = PhyN
    tracer(:,index_recom_tracer(iphyc))  = tiny * Redfield ! tracer 7  = PhyC
    tracer(:,index_recom_tracer(ipchl))  = tiny * 1.56d0   ! tracer 8  = PhyChl
    tracer(:,index_recom_tracer(idetn))  = tiny            ! tracer 9  = DetN
    tracer(:,index_recom_tracer(idetc)) = tiny            ! tracer 10 = DetC
    tracer(:,index_recom_tracer(ihetn)) = tiny            ! tracer 11 = HetN
    tracer(:,index_recom_tracer(ihetc)) = tiny * Redfield ! tracer 12 = HetC
    tracer(:,index_recom_tracer(idon)) = tiny            ! tracer 13 = DON
    tracer(:,index_recom_tracer(idoc)) = tiny            ! tracer 14 = DOC
    tracer(:,index_recom_tracer(idian)) = tiny            ! tracer 15 = DiaN
    tracer(:,index_recom_tracer(idiac)) = tiny * Redfield ! tracer 16 = DiaC
    tracer(:,index_recom_tracer(idchl)) = tiny * 1.56d0   ! tracer 17 = DiaCh
    tracer(:,index_recom_tracer(idiasi)) = tiny            ! tracer 18 = DiaSi
    tracer(:,index_recom_tracer(idetsi)) = tiny            ! tracer 19 = DetSi
    tracer(:,index_recom_tracer(icocn)) = tiny             ! tracer 29 = CoccoN 
    tracer(:,index_recom_tracer(icocc)) = tiny * Redfield  ! tracer 30 = CoccoC
    tracer (:,index_recom_tracer(icchl)) = tiny * 1.56d0   ! tracer 31 = CoccoChl 
!#ifdef REcoM_calcification
    tracer(:,index_recom_tracer(iphycal)) = tiny * Redfield  ! tracer 22 = PhyCalc
    tracer(:,index_recom_tracer(idetcal)) = tiny             ! tracer 23 = DetCalc
!#endif
 if (REcoM_Second_Zoo) then
   if (REcoM_Second_Zoo .and. zoo2_initial_field) then
     tracer(:,index_recom_tracer(idetz2n)) = tiny            ! tracer 26 = DetZ2N
     tracer(:,index_recom_tracer(idetz2c)) = tiny            ! tracer 27 = DetZ2C
     tracer(:,index_recom_tracer(idetz2si)) = tiny            ! tracer 28 = DetZ2Si
     tracer(:,index_recom_tracer(idetz2calc)) = tiny            ! tracer 29 = DetZ2Calc
     else
     tracer(:,index_recom_tracer(izoo2n)) = tiny            ! tracer 24 = Zoo2N                                      
     tracer(:,index_recom_tracer(izoo2c)) = tiny            ! tracer 25 = Zoo2C 
     tracer(:,index_recom_tracer(idetz2n)) = tiny            ! tracer 26 = DetZ2N                              
     tracer(:,index_recom_tracer(idetz2c)) = tiny            ! tracer 27 = DetZ2C                                    
     tracer(:,index_recom_tracer(idetz2si)) = tiny            ! tracer 28 = DetZ2Si                            
     tracer(:,index_recom_tracer(idetz2calc)) = tiny            ! tracer 29 = DetZ2Calc 
   endif
  endif
!---------------------------------------------------------------------------------------------------------
! Initialization of benthos
! Benthic layer consists of Benthos(1) = N, Benthos(2)=C, Benthos(3)=Si, Benthos(4)=calc
  
    do j = 1, size2D
      do i = 1,4
        Benthos(i,j) = tiny
      end do
    end do

    if(recom_mype==0) write(*,*),'Tracers have been initialized as spinup from WOA/glodap netcdf files'
    
  end if ! .not. REcoM_restart

!===============================================================================
! File for alkalinity restoring is opened, surface field is saved
!===============================================================================

! this is now used for surface distribution of all surface fluxes (alk,dic,fe,N,o2), hence independent of alkalinity restoring  
  allocate(recom_sfc_force(recom_ToDim_nod2d,5))

  if (restore_alkalinity) then
    allocate(Alk_surf(size2D))

    ! Alkalinity is initialized
    filename=trim(REcoMDataPath)//'InitAlk.nc'

    status = nf_open(filename, nf_nowrite, ncid)
    if (status .ne. nf_noerr) call handle_err(status)

    status=nf_inq_varid(ncid, 'Alk', alk_varid)
    if (status .ne. nf_noerr) call handle_err(status)

    ! read variables
    istart=(/1,1/)
    icount=(/recom_nod2d, 1/)
    status=nf_get_vara_double(ncid, alk_varid, istart, icount, aux2) 
    if (status .ne. nf_noerr) call handle_err(status)
    Alk_surf = aux2(recom_myList_nod2D)

    status=nf_close(ncid)
    if (status .ne. nf_noerr) call handle_err(status)

    if(recom_mype==0) write(*,*),'Alkalinity restoring = true. Field read in.'

  end if 

!-------------------------------------------------------------------------------

    call RiverineNutrients(Size2D)
 
!-------------------------------------------------------------------------------
! Initializing output netcdf files

    call REcoM_init_output(yearnew, yearold, cyearnew)
  
  deallocate(aux2,aux3)
        
end subroutine recom_init

!===============================================================================

! ====================================================================
  subroutine read_recom_binary(tracer, cyearold)
  use recom_ocean_settings
  use recom_config
  implicit none
   character(4), INTENT(IN) :: cyearold
   real(kind=8), dimension(recom_todim_nod3d,recom_num_tracer), INTENT(INOUT) :: tracer
   integer                   :: fileID
   integer                   :: n, phys_num
   real(kind=8), allocatable :: array_3d(:)
   character(2000)           :: filename
   
   allocate(array_3d(recom_nod3D))
! add variable to determine how many physical tracers exist to determine                         
! the index of the first bgc tracer                                                                  
   phys_num = recom_num_tracer - bgc_num

!   ! Read restart
   if (recom_mype==0) then
     write(*,*) 'Reading restart'
   end if
      fileID=125
      filename = trim(recom_ResultPath)//'recom_restart'//cyearold//'.bin'
      open(fileID,file=filename, form='unformatted')

! make flexible depending on how many physical tracers exist                                     
   do n=1, recom_num_tracer-phys_num
   read(fileID) array_3d
   tracer(:,n+phys_num)=array_3d(recom_myList_nod3D)
   end do

   close(fileID)
   deallocate(array_3d)

  end subroutine read_recom_binary

!===============================================================================
! Riverine input of DIN, DON and DOC is provided in the netcdf file 'RiverineInput.nc'
! (Mayorga et al. (2010)). It entails fluxes (mmol/day) from ~6000 rivers.
! The position of these rivers on the current grid is found in this subroutine, 
! and the flux is added to 2D fields, which have zeros at all other locations.
!------------------------------------------------------------------------------- 

  subroutine RiverineNutrients(Size2D)
  use MPI
  use recom_ocean_settings
  use REcoM_GloVar
  use recom_config
  implicit none

#include "netcdf.inc"
  
  real(kind=8), allocatable :: ncdata(:), lon_fesom(:), lat_fesom(:), mindist(:)
  real(kind=8)              :: lon_orig_loc, lat_orig_loc, lon_fesom_loc, lat_fesom_loc
  real(kind=8)              :: in(2), out(2)
  integer	            :: status, ncid, Size2D, mypeout
  integer                   :: Riverstart(2), Rivercount(2)
  integer                   :: i, row
  character(2000)            :: filename
  character(5)              :: lonvari, latvari, DINvari, DONvari, DOCvari, DSivari

  if (NitrogenSS) then
    filename = trim(REcoMDataPath)//'RiverineInput.nc'
    lonvari  = 'lon'
    latvari  = 'lat'
    DINvari  = 'DIN'
    DONvari  = 'DON'
    DOCvari  = 'DOC'
    DSivari  = 'DSi'

  ! open file
    status=nf_open(filename, nf_nowrite, ncid)
    if (status.ne.nf_noerr)then
      print*,'ERROR: CANNOT READ riverine FILE CORRECTLY !!!!!'
      print*,'Error in opening netcdf file '//filename
      stop
    endif

    allocate(RiverineLonOrig(5917))
    allocate(RiverineLatOrig(5917))
    allocate(RiverineDINOrig(5917))
    allocate(RiverineDONOrig(5917))
    allocate(RiverineDOCOrig(5917))
    allocate(RiverineDSiOrig(5917))
    RiverineLonOrig=0.d0
    RiverineLatOrig=0.d0
    RiverineDINOrig=0.d0
    RiverineDONOrig=0.d0
    RiverineDOCOrig=0.d0
    RiverineDSiOrig=0.d0

    allocate(ncdata(5917))
  
    Riverstart = (/1,1/)
    Rivercount= (/5917,1/)
 
    status=nf_inq_varid(ncid, 'lon', lonvari)
    if (status .ne. nf_noerr) call handle_err(status)
    status=nf_get_vara_double(ncid,lonvari,Riverstart,Rivercount,ncdata)
    RiverineLonOrig(:)=ncdata(:)

    status=nf_inq_varid(ncid, 'lat', latvari)
    if (status .ne. nf_noerr) call handle_err(status)
    status=nf_get_vara_double(ncid,latvari,Riverstart,Rivercount,ncdata)
    RiverineLatOrig(:)=ncdata(:)

    status=nf_inq_varid(ncid, 'DIN', DINvari)
    if (status .ne. nf_noerr) call handle_err(status)
    status=nf_get_vara_double(ncid,DINvari,Riverstart,Rivercount,ncdata)
    RiverineDINOrig(:)=ncdata(:)

    status=nf_inq_varid(ncid, 'DON', DONvari)
    if (status .ne. nf_noerr) call handle_err(status)
    status=nf_get_vara_double(ncid,DONvari,Riverstart,Rivercount,ncdata)
    RiverineDONOrig(:)=ncdata(:)

    status=nf_inq_varid(ncid, 'DOC', DOCvari)
    if (status .ne. nf_noerr) call handle_err(status)
    status=nf_get_vara_double(ncid,DOCvari,Riverstart,Rivercount,ncdata)
    RiverineDOCOrig(:)=ncdata(:)

    status=nf_inq_varid(ncid, 'DSi', DOCvari)
    if (status .ne. nf_noerr) call handle_err(status)
    status=nf_get_vara_double(ncid,DSivari,Riverstart,Rivercount,ncdata)
    RiverineDSiOrig(:)=ncdata(:)

    deallocate(ncdata)
    status=nf_close(ncid)
    if (status .ne. nf_noerr) call handle_err(status)

    ! Finding nearest ocean node
    allocate(RiverDIN2D(size2D))
    allocate(RiverDON2D(size2D))
    allocate(RiverDOC2D(size2D))
    allocate(RiverDSi2D(size2D))
    RiverDIN2D = 0.d0
    RiverDON2D = 0.d0
    RiverDOC2D = 0.d0
    RiverDSi2D = 0.d0
  
    allocate(lon_fesom(recom_mydim_nod2D))
    allocate(lat_fesom(recom_mydim_nod2D))
    if (recom_rotated_grid) then
      do i = 1, recom_myDim_nod2d
        call r2g(lon_fesom_loc, lat_fesom_loc, recom_coord_nod2d(1, i), recom_coord_nod2d(2, i)) !
        lon_fesom(i) = lon_fesom_loc/recom_rad
        lat_fesom(i) = lat_fesom_loc/recom_rad
      end do
    else
      lon_fesom = recom_coord_nod2d(1, 1:Size2D)/recom_rad ! fesom's longitudes on the current processor
      lat_fesom = recom_coord_nod2d(2, 1:Size2D)/recom_rad ! fesom's latitudes on the current processor
    end if

    allocate(Mindist(5917))
    do i =1,5917
      lon_orig_loc = RiverineLonOrig(i)
      lat_orig_loc = RiverineLatOrig(i)
      MinDist(i) = minval(sqrt((lon_fesom - lon_orig_loc)**2+(lat_fesom - lat_orig_loc)**2))
      in = (/MinDist(i),real(recom_mype, 8)/)
    
      call MPI_Barrier(MPI_COMM_RECOM,MPIerr_recom)
      call MPI_AllREDUCE(in, out, 1, MPI_2DOUBLE_PRECISION, MPI_MINLOC, &
         MPI_COMM_RECOM, MPIerr_recom)
    
      mypeout = int(out(2))
      if (recom_mype==mypeout) then
         row = minloc((lon_fesom - lon_orig_loc)**2+(lat_fesom - lat_orig_loc)**2,1)
         RiverDIN2D(row)= RiverDIN2D(row) + RiverineDINOrig(i)
         RiverDON2D(row)= RiverDON2D(row) + RiverineDONOrig(i)
         RiverDOC2D(row)= RiverDOC2D(row) + RiverineDOCOrig(i)
         RiverDSi2D(row)= RiverDSi2D(row) + 0.7*RiverineDSiOrig(i)
      end if
    end do

    deallocate(Mindist)

    deallocate(RiverineLonOrig)
    deallocate(RiverineLatOrig)
    deallocate(RiverineDINOrig)
    deallocate(RiverineDONOrig)
    deallocate(RiverineDOCOrig)
    deallocate(RiverineDSiOrig)

    if(recom_mype==0) write(*,*),'Finished riverine input'
  else ! if not NitrogenSS
    allocate(RiverDIN2D(size2D))
    allocate(RiverDON2D(size2D))
    allocate(RiverDOC2D(size2D))
    allocate(RiverDSi2D(size2D))
    RiverDIN2D = 0.d0
    RiverDON2D = 0.d0
    RiverDOC2D = 0.d0
    RiverDSi2D = 0.d0
  end if
  end subroutine RiverineNutrients

