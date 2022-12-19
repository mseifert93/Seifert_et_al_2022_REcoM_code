!> \file recom.f90
!! \BRIEF 
!> ===============================================================================
!! Here REcoM is called
!! ===============================================================================
subroutine call_REcoM(shortwave                                        &
            , u_wind, v_wind                                           &
            , slp                                                      &
#ifdef __cpl_echam
            , co2c_oce, w10w_oce                                       &
#endif
            , a_ice, tracer                                            &
            , month, yearnew, yearold, cyearnew, daynew                &
            , ndpyr, day_in_month, timeold                             &
            , aux_todim_nod3d, aux_num_tracer, recom_istep, mocsy_restart, mocsy_step_per_day )   

  use recom_ocean_settings
  use REcoM_declarations
  use REcoM_LocVar
  use REcoM_GloVar
  use recom_config

  implicit none

  integer, INTENT(IN) :: month, yearnew, yearold, daynew, ndpyr, day_in_month
! interface needs these too, doesn't recognize that they are part of the module:
! but of course need a different name then the one loaded by the module
  integer, INTENT(IN) :: aux_todim_nod3d, aux_num_tracer, recom_istep, mocsy_step_per_day    
  logical :: mocsy_restart                                                            
  real(kind=8), INTENT(IN) :: timeold
  character(4), INTENT(IN) :: cyearnew
  real(kind=8), dimension(:), INTENT(IN) :: shortwave, u_wind, v_wind, slp
#ifdef __cpl_echam 
  real(kind=8), dimension(:), INTENT(IN) :: w10w_oce, co2c_oce
#endif
  real(kind=8), dimension(:), INTENT(IN) :: a_ice
  real(kind=8), dimension(aux_todim_nod3d,aux_num_tracer), INTENT(INOUT) :: tracer

  integer         :: i, row, nn, inds(recom_max_num_layers),n,idiags, phys_num, aux
  real(kind=8)    :: z(recom_max_num_layers), C(recom_max_num_layers,bgc_num)
  real(kind=8)    :: Temp(recom_max_num_layers,1),Sali,Sali_depth(recom_max_num_layers,1)        
  real(kind=8)    :: lat, lon, SW, Loc_slp,vol, PAR(recom_max_num_layers), CO2_watercolumn(recom_max_num_layers) 
  real(kind=8)    :: pH_watercolumn(recom_max_num_layers), pCO2_watercolumn(recom_max_num_layers), HCO3_watercolumn(recom_max_num_layers)
  real(kind=8)    :: CO3_watercolumn(recom_max_num_layers)                                               
  real(kind=8)    :: OmegaC_watercolumn(recom_max_num_layers), kspc_watercolumn(recom_max_num_layers), rhoSW_watercolumn(recom_max_num_layers)
  real(kind=8)    :: Nutlim_phy(recom_max_num_layers), Nutlim_dia(recom_max_num_layers), Nutlim_cocco(recom_max_num_layers) 
  real(kind=8)    :: Tlim_arr(recom_max_num_layers), Tlim_cocco(recom_max_num_layers)                             
  real(kind=8)    :: Llim_phy(recom_max_num_layers), Llim_dia(recom_max_num_layers), Llim_cocco(recom_max_num_layers) 
  real(kind=8)    :: CO2lim_phy(recom_max_num_layers), CO2lim_dia(recom_max_num_layers), CO2lim_cocco(recom_max_num_layers)
  real(kind=8)    :: PR_phy(recom_max_num_layers), PR_dia(recom_max_num_layers), PR_cocco(recom_max_num_layers)         
  real(kind=8)    :: Cal_Tlim(recom_max_num_layers), Cal_CO2lim(recom_max_num_layers), Cal_Nlim(recom_max_num_layers)  
  real(kind=8)    :: Cal_pure(recom_max_num_layers)                                                               
  
  ! ---- use_REcom = .false. : Disable recom  ----
  if (.not. use_REcoM) return

  ! add variable to determine how many physical tracers exist to determine
  ! the index of the first bgc tracer (to make age tracer work with REcoM)                           
  phys_num = recom_num_tracer - bgc_num

  !> Initialize REcoM diagnostics
  call REcoM_init_output(yearnew, yearold, cyearnew)

  !> Restore surface alkalinity
  call alk_restore(tracer)  

  !> Read atmospheric deposition fields
  !! a) iron dust 
  !! b) co2 
  !! c) aeolian nitrogen 
  call Atm_input(day_in_month, timeold, yearnew, cyearnew, month, recom_istep)

  ! ---- Each surface node  ----
  do row = 1,recom_myDim_nod2d

     !> the visible part (350nm-700nm for REcoM)                                                
     SW                = parFrac*shortwave(row)
     !> assume no shortwave penetration under the ice (cavities are not taken into account)
     SW                = SW * (1.d0 - a_ice(row))
     !> sea ice concentration
     Loc_ice_conc      = a_ice(row)
     !> sea level pressure
     Loc_slp           = slp(row)
!     if(recom_mype==1) write(*,*), 'SW1=',SW
!     if(recom_mype==1) write(*,*), 'slp=',slp

     if (recom_rotated_grid) then
!     if(recom_mype==1) write(*,*), 'in rotated grid'
        !> rotated grid 
        !! change model grid to geographical coordinates (rad)
        call r2g(lon, lat, recom_coord_nod2d(1, row), recom_coord_nod2d(2, row))
     else
        lon            = recom_coord_nod2d(1, row)
        lat            = recom_coord_nod2d(2, row)
     end if

     !> number of vertical layers in current water column
     nn                = recom_num_layers_below_nod2D(row)+1

     !> Nodes are placed in vertical 
     !! columns below the surface nodes
     
     !> get 3d node indices
     inds(1:nn)        = recom_nod3D_below_nod2D(1:nn,row)
     !> potential temperature in current water column                                         
     Temp(1:nn,1)      = tracer(inds(1:nn),1)
     ! surface salinity
     Sali              = tracer(inds(1),2)
     ! salinity for the whole water column                  
     Sali_depth(1:nn,1)= tracer(inds(1:nn),2)                      
     ! CO2 for the whole watercolumn                           
     CO2_watercolumn(1:nn)    = CO23D(inds(1:nn))              
     pH_watercolumn(1:nn)     = pH3D(inds(1:nn))                   
     pCO2_watercolumn(1:nn)   = pCO23D(inds(1:nn))             
     HCO3_watercolumn(1:nn)   = HCO33D(inds(1:nn))               
     CO3_watercolumn(1:nn)    = CO33D(inds(1:nn))              
     OmegaC_watercolumn(1:nn) = OmegaC3D(inds(1:nn))            
     kspc_watercolumn(1:nn)   = kspc3D(inds(1:nn))               
     rhoSW_watercolumn(1:nn)  = rhoSW3D(inds(1:nn))             
     ! Additional output to analyse limitations                
     Nutlim_phy(1:nn)         = Nutlim_phy3D(inds(1:nn))        
     Nutlim_dia(1:nn)         = Nutlim_dia3D(inds(1:nn))      
     Nutlim_cocco(1:nn)       = Nutlim_cocco3D(inds(1:nn))    
     Tlim_arr(1:nn)           = Tlim_arr3D(inds(1:nn))        
     Tlim_cocco(1:nn)         = Tlim_cocco3D(inds(1:nn))     
     Llim_phy(1:nn)           = Llim_phy3D(inds(1:nn))           
     Llim_dia(1:nn)           = Llim_dia3D(inds(1:nn))                  
     Llim_cocco(1:nn)         = Llim_cocco3D(inds(1:nn))          
     CO2lim_phy(1:nn)         = CO2lim_phy3D(inds(1:nn))            
     CO2lim_dia(1:nn)         = CO2lim_dia3D(inds(1:nn))           
     CO2lim_cocco(1:nn)       = CO2lim_cocco3D(inds(1:nn))          
     PR_phy(1:nn)             = PR_phy3D(inds(1:nn))                 
     PR_dia(1:nn)             = PR_dia3D(inds(1:nn))                
     PR_cocco(1:nn)           = PR_cocco3D(inds(1:nn))                  
     Cal_Tlim(1:nn)           = Cal_Tlim3D(inds(1:nn))                  
     Cal_CO2lim(1:nn)         = Cal_CO2lim3D(inds(1:nn))         
     Cal_Nlim(1:nn)           = Cal_Nlim3D(inds(1:nn))             
     Cal_pure(1:nn)           = Cal_pure3D(inds(1:nn))          

     !> C is set equal to the concentration of tracer in current water column
     C(1:nn,1:bgc_num) = tracer(inds(1:nn), phys_num+1:bgc_num+phys_num) 

!     C(1:nn,1:bgc_num) = tracer(inds(1:nn), 3:bgc_num+2)                 

     !> Depth of the nodes in the column are calculated (negative)   
     z(1:nn)           = recom_coord_nod3D(3, inds(1:nn))

     !> The PAR in the current water column is initialized
     PAR(1:nn)       = 0.d0
     LocBenthos(1:4)   = Benthos(1:4,row)

     !> Local conc of [H+]-ions from last time time step. 
     !! Stored in LocVar, used as first guess for H+ conc.in subroutine CO2flux
     Hplus             = GloHplus(row)

     !> Wind from FESOM is temporarily stored in module LocVar                          
#ifdef __cpl_echam 
     ULoc              = w10w_oce(row) !* 0.25 !Sensitivity OG
#else
     ULoc              = sqrt(u_wind(row)**2+v_wind(row)**2)            
#endif

     !> a_ice(row): Sea ice concentration in the local node
     FeDust            = GloFeDust(row) * dust_sol * (1 - a_ice(row))  ! Iron
     NDust             = GloNDust(row) * (1 - a_ice(row))              ! Nitrogen

     !> co2 concentration from the atmosphere
#ifdef __cpl_echam 

     aux = 28.970/44.0095*1000000 

     LocAtmCO2         = co2c_oce(row) * aux !* 28.970 / 44.0095 * 1000000 
#else
     LocAtmCO2         = AtmCO2(month)
#endif

     !> Riverine input of DIN, DON and DOC is provided in the netcdf file 'RiverineInput.nc'
     !! (Mayorga et al. (2010)). It entails fluxes (mmol/day) from ~6000 rivers.
     !! The position of these rivers on the current grid is found in this subroutine, 
     !! and the flux is added to 2D fields, which have zeros at all other locations.

     vol=0.d0
     do n=1, nn
        vol=vol+sum(recom_voltetra(recom_nod_in_elem3D(n)%addresses)/4.d0)
     end do
     
     LocRiverDIN       = RiverDIN2D(row)/vol 
     LocRiverDON       = RiverDON2D(row)/vol
     LocRiverDOC       = RiverDOC2D(row)/vol
     LocRiverDSi       = RiverDSi2D(row)/vol

     allocate(Diags3Dloc(Nn,diags3d_num)) 
     Diags3Dloc(:,:) = 0.d0

     !> call REcoM_Forcing
     call REcoM_Forcing(lon, lat, z(1:nn)                              & ! lon, lat in rad, z is negative
           , nn                                                        & ! total number of vertical nodes 
           , C(1:nn,1:bgc_num)                                         & ! tracer fields excluding temperature and salinity 
           , SW                                                        & ! Shortwave radiation on the ocean surface [W/m2]
           , Loc_slp                                                   & ! sea level pressure 
           , recom_dt                                                  & ! [s]ocean timestep -> recom timestep
           , Temp                                                      & ! Temperature in water column [degrees C] 
           , Sali                                                      & ! Sea surface salinity
           , Sali_depth                                                & ! Salinity for the whole depth
           , CO2_watercolumn(1:nn)                                     & ! CO2 for the whole watercolumn
           , pH_watercolumn(1:nn)                                      & ! pH for the whole watercolumn
           , pCO2_watercolumn(1:nn)                                    & ! pCO2 for the whole watercolumn
           , HCO3_watercolumn(1:nn)                                    & ! HCO3 for the whole watercolumn
           , CO3_watercolumn(1:nn)                                     & ! CO3 for the whole watercolumn
           , OmegaC_watercolumn(1:nn)                                  & ! OmegaC for the whole watercolumn
           , kspc_watercolumn(1:nn)                                    & ! stoichiometric solubility product for calcite [mol^2/kg^2]
           , rhoSW_watercolumn(1:nn)                                   & ! in-situ density of seawater [mol/m^3]
           , Nutlim_phy(1:nn)                                          & ! nutrient limitation of small phytoplankton
           , Nutlim_dia(1:nn)                                          & ! nutrient limitation of diatoms
           , Nutlim_cocco(1:nn)                                        & ! nutrient limitiation of coccolithophores
           , Tlim_arr(1:nn)                                            & ! temperature limitation according to the Arrhenius function
           , Tlim_cocco(1:nn)                                          & ! temperature limitation of coccolithophores
           , Llim_phy(1:nn)                                            & ! light limitation of small phytoplankton
           , Llim_dia(1:nn)                                            & ! light limitation of diatoms
           , Llim_cocco(1:nn)                                          & ! light limitation of coccolithophores
           , CO2lim_phy(1:nn)                                          & ! CO2 limitation of small phytoplankton
           , CO2lim_dia(1:nn)                                          & ! CO2 limitation of diatoms
           , CO2lim_cocco(1:nn)                                        & ! CO2 limitation of coccolithophores
           , PR_phy(1:nn)                                              & ! Photosynthesis rate of small phytoplankton
           , PR_dia(1:nn)                                              & ! Photosynthesis rate of diatoms
           , PR_cocco(1:nn)                                            & ! Photosynthesis rate of coccolithophores
           , Cal_Tlim(1:nn)                                            & ! Temperature dependence of calcification
           , Cal_CO2lim(1:nn)                                          & ! CO2 dependence of calcification
           , Cal_Nlim(1:nn)                                            & ! Nitrate dependence of calcification
           , Cal_pure(1:nn)                                            & ! PIC only dependent on PICPOCmax,  CoccoC, T, N, CO2
           , PAR(1:nn)                                                 & ! Photosynthetically available radiation
           , recom_istep                                               & ! Timestep in REcoM used for carbonate system calculation
           , mocsy_restart                                             & ! Marks the very first year of a run
           , mocsy_step_per_day                                        & ! 
           , daynew                                                    & ! updated day (for cos of angle of incidence -> average PAR )
           , ndpyr)                                                      ! number of days in yearnew (for cos of angle of incidence -> average PAR)

     !> REcoM is called and changes concentrations of tracer                               
     tracer(inds(1:nn),phys_num+1:recom_num_tracer) = C(1:nn,1:bgc_num)

     !> Local variables that have been changed during the time-step are stored so they can be saved
     Benthos(1:4,row)    = LocBenthos(1:4)                                ! Updating Benthos values
     Diags2D(:,row)      = LocDiags2D(:)                                  ! Updating diagnostics
     GloPCO2surf(row)    = pco2surf(1)
     GloCO2flux(row)     = dflux(1)
     GloCO2flux_seaicemask(row)  = co2flux_seaicemask(1)                  !  [mmol/m2/s]
#ifdef __cpl_echam 
     GloCO2flux_noicemask(row)   = co2flux(1)*44.0095*0.001 ! in mocsy without sea-ice mol/m2/s -> kg/m2/s
#endif
     GloO2flux_seaicemask(row)   = o2flux_seaicemask(1)                  !  [mmol/m2/s]
     GlodPCO2surf(row)           = dpco2surf(1)
     PistonVelocity(row)         = kw660(1) ! write Piston velocity                                
     alphaCO2(row)               = K0(1)  ! write CO2 solubility                                   

     GloHplus(row)              = hplus
     AtmFeInput(row)            = FeDust
     AtmNInput(row)             = NDust 
     DenitBen(row)              = LocDenit
     PAR3D(inds(1:nn))          = PAR(1:nn)
     CO23D(inds(1:nn))          = CO2_watercolumn(1:nn)       
     pH3D(inds(1:nn))           = pH_watercolumn(1:nn)        
     pCO23D(inds(1:nn))         = pCO2_watercolumn(1:nn)        
     HCO33D(inds(1:nn))         = HCO3_watercolumn(1:nn)         
     CO33D(inds(1:nn))          = CO3_watercolumn(1:nn)         
     OmegaC3D(inds(1:nn))       = OmegaC_watercolumn(1:nn)     
     kspc3D(inds(1:nn))         = kspc_watercolumn(1:nn)     
     rhoSW3D(inds(1:nn))        = rhoSW_watercolumn(1:nn)       
     Nutlim_phy3D(inds(1:nn))   = Nutlim_phy(1:nn)                
     Nutlim_dia3D(inds(1:nn))   = Nutlim_dia(1:nn)            
     Nutlim_cocco3D(inds(1:nn)) = Nutlim_cocco(1:nn)          
     Tlim_arr3D(inds(1:nn))     = Tlim_arr(1:nn)       
     Tlim_cocco3D(inds(1:nn))   = Tlim_cocco(1:nn)               
     Llim_phy3D(inds(1:nn))     = Llim_phy(1:nn)           
     Llim_dia3D(inds(1:nn))     = Llim_dia(1:nn)            
     Llim_cocco3D(inds(1:nn))   = Llim_cocco(1:nn)          
     CO2lim_phy3D(inds(1:nn))   = CO2lim_phy(1:nn)          
     CO2lim_dia3D(inds(1:nn))   = CO2lim_dia(1:nn)          
     CO2lim_cocco3D(inds(1:nn)) = CO2lim_cocco(1:nn)      
     PR_phy3D(inds(1:nn))       = PR_phy(1:nn)             
     PR_dia3D(inds(1:nn))       = PR_dia(1:nn)           
     PR_cocco3D(inds(1:nn))     = PR_cocco(1:nn)            
     Cal_Tlim3D(inds(1:nn))     = Cal_Tlim(1:nn)               
     Cal_CO2lim3D(inds(1:nn))   = Cal_CO2lim(1:nn)         
     Cal_Nlim3D(inds(1:nn))     = Cal_Nlim(1:nn)          
     Cal_pure3D(inds(1:nn))     = Cal_pure(1:nn)             
     do idiags = 1,diags3d_num
      Diags3D(inds(1:nn),idiags) = Diags3Dloc(1:nn,idiags) ! 1=NPPnano, 2=NPPdia, 3=NPPcocco, 20 diagnostics see further below recom_sms
     end do    
    
    deallocate(Diags3Dloc)
     
  end do

  call surface_fluxes

  !> recom_num_tracer is number of all tracers incl T and S   
  !! communicate to other processors
  do i=phys_num+1, recom_num_tracer 
    call com_3D(tracer(:, i))  
  end do

end subroutine call_REcoM
!===============================================================================
! REcoM_Forcing
!===============================================================================
subroutine REcoM_Forcing(Lonr, Latr, zNodes                            &
            , Nn                                                       &
            , state                                                    &
            , SurfSW                                                   &
            , Loc_slp                                                  &
            , dt_s                                                     &
            , Temp                                                     &
            , Sali                                                     &
            , Sali_depth                                               &  
            , CO2_watercolumn                                          &   
            , pH_watercolumn                                           &    
            , pCO2_watercolumn                                         &       
            , HCO3_watercolumn                                         &   
            , CO3_watercolumn                                          & 
            , OmegaC_watercolumn                                       &  
            , kspc_watercolumn                                         &  
            , rhoSW_watercolumn                                        &  
            , Nutlim_phy                                               &
            , Nutlim_dia                                               &   
            , Nutlim_cocco                                             &      
            , Tlim_arr                                                 &     
            , Tlim_cocco                                               &   
            , Llim_phy                                                 &     
            , Llim_dia                                                 &    
            , Llim_cocco                                               &    
            , CO2lim_phy                                               &    
            , CO2lim_dia                                               &    
            , CO2lim_cocco                                             &    
            , PR_phy                                                   &     
            , PR_dia                                                   &   
            , PR_cocco                                                 &     
            , Cal_Tlim                                                 &   
            , Cal_CO2lim                                               &   
            , Cal_Nlim                                                 &   
            , Cal_pure                                                 &    
            , PAR                                                      &
            , recom_istep                                              &    
            , mocsy_restart                                            &   
            , mocsy_step_per_day                                       &   
            , daynew                                                   &
            , ndpyr)
  use recom_ocean_settings
  use REcoM_declarations
  use REcoM_LocVar
  use recom_config
  use GASX
  Implicit none

  integer, INTENT(IN) :: daynew, ndpyr, recom_istep, mocsy_step_per_day
  logical             :: mocsy_restart

! From FESOM
  Real(kind=8)                       :: Latr
  Real(kind=8)                       :: Lonr
  Integer                            :: Nn		     ! Total number of nodes
  Real(kind=8),dimension(Nn)	     :: zNodes		     ! Depth of nodes
  Real(kind=8),dimension(Nn,bgc_num) :: state                ! ChlA conc in phytoplankton [mg/m3]	
  Real(kind=8)                       :: SurfSW               ! [W/m2] ShortWave radiation at surface
  Real(kind=8)                       :: Loc_slp              ! [Pa] se-level pressure
  Real(kind=8)			     :: dt_s                 ! Size of time steps [s]
  Real(kind=8),dimension(Nn)         :: Temp                 ! [degrees C] Ocean temperature
  Real(kind=8),dimension(Nn)         :: Sali_depth           ! Salinity for the entire water column
  real(kind=8),dimension(Nn)         :: PAR
! Subroutine Cobeta
  Real(kind=8)                       :: cosAI                ! cos of angle of incidence
! Subroutine Depth
  Real(kind=8),dimension(Nn+1)       :: zF                   ! [m] Depth of fluxes
  Real(kind=8),dimension(Nn)         :: dzF                  ! [m] Thickness of control volumes  
  Real(kind=8),dimension(Nn)         :: recipdzF             ! [1/m]
  Real(kind=8),dimension(Nn+1,5)     :: SinkVel              ! [m/day]
  Real(kind=8),dimension(Nn-1)       :: thick                ! [m] Vertical distance between two nodes = Thickness 
  Real(kind=8),dimension(Nn-1)         :: recipthick         ! [1/m] reciprocal of thick
! Subroutine CO2Flux /mocsy
  Real(kind=8)                       :: REcoM_DIC(1)            ! [mmol/m3] Conc of DIC in the surface water, used to calculate CO2 flux
  Real(kind=8)                       :: REcoM_Alk(1)            ! [mmol/m3] Conc of Alk in the surface water, used to calculate CO2 flux
  Real(kind=8)                       :: REcoM_Si(1)             ! [mol/m3] Conc of Si in the surface water, used to calculate CO2 flux
  Real(kind=8)                       :: REcoM_Phos(1)           ! [mol/m3] Conc of Phos in the surface water, used to calculate the CO2 flux
  Real(kind=8)                       :: Sali(1)                 ! Salinity of current surface layer
  Real(kind=8)                       :: Latd(1)                 ! latitude in degree
  Real(kind=8)                       :: Lond(1)                 ! longitude in degree
  Real(kind=8)                       :: REcoM_T(1)                 ! temperature again, for mocsy minimum defined as -2
  Real(kind=8)                       :: REcoM_S(1)                 ! temperature again, for mocsy minimum defined as 21
! atm pressure, now read in as forcing!!
  Real(kind=8)                       :: Patm(1)                 ! atmospheric pressure [atm]
! Subroutine o2flux /mocsy 
 Real(kind=8)                       :: ppo(1)                 ! atmospheric pressure, divided by 1 atm 
 Real(kind=8)                       :: REcoM_O2(1)            ! [mmol/m3] Conc of O2 in the surface water, used to calculate O2 flux
  
! Subroutine REcoM_sms
  Real(kind=8),dimension(Nn,bgc_num) :: sms                  ! matrix that entail changes in tracer concentrations
!Diagnostics
  Integer                            :: idiags
  Integer                            :: diags2d_num
! For entire watercolumn carbonate chemistry
  Real(kind=8),dimension(Nn)         :: CO2_watercolumn  
  Real(kind=8),dimension(Nn)         :: pH_watercolumn     
  Real(kind=8),dimension(Nn)         :: pCO2_watercolumn  
  Real(kind=8),dimension(Nn)         :: HCO3_watercolumn  
  Real(kind=8),dimension(Nn)         :: CO3_watercolumn   
  Real(kind=8),dimension(Nn)         :: OmegaC_watercolumn  
  Real(kind=8),dimension(Nn)         :: kspc_watercolumn   
  Real(kind=8),dimension(Nn)         :: rhoSW_watercolumn
! New output to analyse limitations
  Real(kind=8),dimension(Nn)         :: Nutlim_phy      
  Real(kind=8),dimension(Nn)         :: Nutlim_dia     
  Real(kind=8),dimension(Nn)         :: Nutlim_cocco    
  Real(kind=8),dimension(Nn)         :: Tlim_arr        
  Real(kind=8),dimension(Nn)         :: Tlim_cocco     
  Real(kind=8),dimension(Nn)         :: Llim_phy       
  Real(kind=8),dimension(Nn)         :: Llim_dia       
  Real(kind=8),dimension(Nn)         :: Llim_cocco     
  Real(kind=8),dimension(Nn)         :: CO2lim_phy      
  Real(kind=8),dimension(Nn)         :: CO2lim_dia      
  Real(kind=8),dimension(Nn)         :: CO2lim_cocco    
  Real(kind=8),dimension(Nn)         :: PR_phy           
  Real(kind=8),dimension(Nn)         :: PR_dia          
  Real(kind=8),dimension(Nn)         :: PR_cocco        
  Real(kind=8),dimension(Nn)         :: Cal_Tlim       
  Real(kind=8),dimension(Nn)         :: Cal_CO2lim     
  Real(kind=8),dimension(Nn)         :: Cal_Nlim       
  Real(kind=8),dimension(Nn)         :: Cal_pure      

  !> compute cos(AngleOfIncidence)
  !! it is used to average PAR
  call Cobeta(cosAI,Latr,daynew,ndpyr)

  !>  Depth and thickness computations
  call Depth_calculations(zNodes,Nn,zF,dzF,recipdzF,SinkVel,thick,recipthick)

!-------------------------------------------------------------------------------
! Atmospheric input of CO2
!  REcoM_DIC = max(tiny,state(one,idic))
!  REcoM_Alk = max(tiny,state(one,ialk))
!  call CO2Flux(REcoM_DIC,REcoM_Alk,Temp(one),Sali)
!  dflux     = dflux * (1.d0 - Loc_ice_conc)

  !>----- new:mocsy -----! 
  !! convert from mmol/m3 to mol/m3 
  !! surface layer

  REcoM_DIC  = max(tiny*1e-3,state(one,idic)*1e-3)
  REcoM_Alk  = max(tiny*1e-3,state(one,ialk)*1e-3)
  REcoM_Si   = max(tiny*1e-3,state(one,isi)*1e-3)
  REcoM_Phos = max(tiny*1e-3,state(one,idin)*1e-3) /16 ! convert N to P with Redfield ratio
  REcoM_T    = max(2.d0, Temp(1)) ! minimum set to 2 degC: K1/K2 Lueker valid between 2degC-35degC and 19-43psu
  REcoM_T    = min(REcoM_T, 40.d0) ! maximum set to 40 degC: K1/K2 Lueker valid between 2degC-35degC and 19-43psu
  REcoM_S    = max(21.d0, Sali(1)) ! minimum set to 21: K1/K2 Lueker valid between 2degC-35degC and 19-43psu, else causes trouble in regions with S between 19 and 21 and ice conc above 97%
  REcoM_S    = min(REcoM_S, 43.d0) ! maximum set to 43: K1/K2 Lueker valid between 2degC-35degC and 19-43psu, else causes trouble 

  Patm       = Loc_slp/Pa2atm ! convert from Pa to atm. 
  Lond       = Lonr/recom_rad ! convert from rad to degree
  Latd       = Latr/recom_rad ! convert from rad to degree
! first calculate piston velocity kw660, which is an input to the flxco2 calculation:
! pistonvel already scaled for ice-free area:
  call pistonvel(ULoc, Loc_ice_conc, Nmocsy, kw660)

  call flxco2(co2flux, co2ex, dpco2,                                                    &
                  ph, pco2surf, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD, rhoSW, p, tempis, K0, &
                  REcoM_T, REcoM_S, REcoM_Alk, REcoM_DIC, REcoM_Si, REcoM_Phos, kw660, LocAtmCO2, Patm, thick(One), Nmocsy, Lond,Latd,     &
                  optCON='mol/m3',optT='Tpot   ',optP='m ',optB='u74',optK1K2='l  ',optKf='dg',optGAS='Pinsitu',optS='Sprc')
! changed optK1K2='l  ' to 'm10'
  if((co2flux(1)>1.e10) .or. (co2flux(1)<-1.e10)) then
      print*,'ERROR: co2 flux !'
      print*,'pco2surf: ',pco2surf
      print*,'co2: ',co2
      print*,'rhoSW: ', rhoSW
      print*,'temp: ',REcoM_T
      print*,'tempis: ',tempis
      print*,'K0: ',K0
      print*,'REcoM_S: ', REcoM_S
      print*, 'REcoM_Alk: ', REcom_Alk
      print*, 'REcoM_DIC: ', REcoM_DIC
      print*, 'REcoM_Si: ', REcoM_Si
      print*, 'REcoM_Phos: ', REcoM_Phos
      print*, 'kw660: ',kw660
      print*, 'LocAtmCO2: ', LocAtmCO2
      print*, 'Patm: ', Patm
      print*, 'dzF(One): ',dzF(One) 
      print*, 'Nmocsy: ', Nmocsy
      print*, 'Lond: ', Lond
      print*, 'Latd: ', Latd   
      print*, 'ULoc: ', ULoc
      print*, 'Loc_ice_conc: ', Loc_ice_conc
      stop
    endif

! use ice-free area and also convert from mol/m2/s to mmol/m2/d
!   if(recom_mype==1) write(*,*), 'co2flux (mol/m2/s) =',co2flux
! ice-fraction is already considered in piston-velocity, so don't apply it here
   dflux     = co2flux *1.e3 *SecondsPerDay  !* (1.d0 - Loc_ice_conc)

!   if(recom_mype==1) write(*,*), 'dflux (mmol/m2/d) =',dflux
   co2flux_seaicemask = co2flux * 1.e3 !  [mmol/m2/s]  * (1.d0 - Loc_ice_conc)
   dpco2surf = dpco2 ![uatm]

! then oxygen
   ppo = Loc_slp/Pa2atm !1 !slp divided by 1 atm
   REcoM_O2 = max(tiny*1e-3,state(one,ioxy)*1e-3) ! convert from mmol/m3 to mol/m3 for mocsy

   call  o2flux(REcoM_T, REcoM_S, kw660, ppo, REcoM_O2, Nmocsy, o2ex)

  !    **********************************************************************
  !    Purpose: Compute time rate of change of O2 in the surface layer due to air-sea gas exchange [mol/(m^3 *s)]                                                         
  !    Input:
  !      T       model surface temperature (deg C)
  !      S       model surface salinity (permil)
  !      kw660   gas transfer velocity at a Schmidt number of 660, accounting
  !              for sea ice fraction (m/s)
  !      ppo     surface pressure divided by 1 atm.
  !      o2      surface ocean O2 concentration (mol/m^3) 
  !      dz1     thickness of surface grid box (m) - taken out (jh)
  !    Output:
  !      o2ex    time rate of change of oxygen in the surface layer due
  !              to air-sea exchange (mol/m^2/s) (not per m3 -jh)       

   o2flux_seaicemask = o2ex * 1.e3 ! back to mmol here [mmol/m2/s] 


!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Source-Minus-Sinks
  call REcoM_sms(cosAI                                                 & ! cos of angle of incidence
        , Nn                                                           & ! Total number of nodes in the vertical
        , state                                                        & ! tracer fields excluding temperature and salinity
        , dzF                                                          & ! [m] Thickness of control volumes
        , recipdzF                                                     & ! [1/m]
        , thick                                                        & ! [m] Vertical distance between nodes = Thickness
        , recipthick                                                   & ! [1/m] reciprocal of thick
        , SurfSW                                                       & ! [W/m2] ShortWave radiation at surface
        , dt_s                                                         & ! [s] timestep
        , sms                                                          & ! Source-Minus-Sinks term
        , Temp                                                         & ! [degrees C] Ocean temperature
        , Sali_depth                                                   & 
        , CO2_watercolumn                                              & ! [mol/m3]
        , pH_watercolumn                                               & ! on total scale
        , pCO2_watercolumn                                             & ! [uatm]
        , HCO3_watercolumn                                             & ! [mol/m3]
        , CO3_watercolumn                                              & ! [mol/m3]
        , OmegaC_watercolumn                                           & ! calcite saturation state
        , kspc_watercolumn                                             & ! stoichiometric solubility product [mol^2/kg^2]
        , rhoSW_watercolumn                                            & ! in-situ density of seawater [kg/m3]
        , Nutlim_phy                                                   & ! nutrient limitation of small phytoplankton
        , Nutlim_dia                                                   & ! nutrient limitation of diatoms
        , Nutlim_cocco                                                 & ! nutrient limitation of coccolithophores
        , Tlim_arr                                                     & ! temperature limitation according to Arrhenius function
        , Tlim_cocco                                                   & ! temperature limitation of coccolithophores
        , Llim_phy                                                     & ! light limitation of small phytoplankton
        , Llim_dia                                                     & ! light limitation of diatoms
        , Llim_cocco                                                   & ! light limitation of coccolithophores
        , CO2lim_phy                                                   & ! CO2 limitation of small phytoplankton
        , CO2lim_dia                                                   & ! CO2 limitation of diatoms
        , CO2lim_cocco                                                 & ! CO2 limitation of coccolithophores
        , PR_phy                                                       & ! Photosynthesis rate of small phytoplankton
        , PR_dia                                                       & ! Photosynthesis rate of diatoms
        , PR_cocco                                                     & ! Photosynthesis rate of coccolithophores
        , Cal_Tlim                                                     & ! Temperature dependence of calcification
        , Cal_CO2lim                                                   & ! CO2 dependence of calcification
        , Cal_Nlim                                                     & ! Nitrate dependence of calcification
        , Cal_pure                                                     & ! PIC only dependent on PICPOCmax, CoccoC, T, N, CO2
        , SinkVel                                                      &
        , zF                                                           & ! [m] Depth of fluxes
        , PAR                                                          &
        , Latd                                                         & ! [degree] latitude 
        , Lond                                                         & ! [degree] longitude
        , recom_istep                                                  & 
        , mocsy_restart                                                & 
        , mocsy_step_per_day                                           & 
        , Loc_slp                                                      & ! [Pa] sea-level pressure
        , daynew                                                       & 
        , Latr)                                                           ! [radyan] latitude 

  state = max(tiny,state + sms)

  state(:,ipchl)  = max(tiny_chl,state(:,ipchl))
  state(:,iphyn)  = max(tiny_N,  state(:,iphyn))
  state(:,iphyc)  = max(tiny_C,  state(:,iphyc))
  state(:,idchl)  = max(tiny_chl,state(:,idchl))
  state(:,idian)  = max(tiny_N_d,state(:,idian))
  state(:,idiac)  = max(tiny_C_d,state(:,idiac))
  state(:,idiasi) = max(tiny_Si, state(:,idiasi))
  state(:,icchl)  = max(tiny_chl,state(:,icchl))                       
  state(:,icocn)  = max(tiny_N_c,state(:,icocn))                     
  state(:,icocc)  = max(tiny_C_c,state(:,icocc))                    
  
!-------------------------------------------------------------------------------
! Diagnostics
  diags2d_num = 12                                         
  if (Diagnostics) then
	do idiags = one,diags2d_num
	  LocDiags2D(idiags) = sum(diags3Dloc(:,idiags) * dzF)
	end do
  end if

end subroutine REcoM_Forcing

!===============================================================================
! Subroutine REcoM_sms
!===============================================================================
subroutine REcoM_sms(cosAI                                             &
            , Nn                                                       &
            , state                                                    &
            , dzF                                                      &
            , recipdzF                                                 &
            , thick                                                    &
            , recipthick                                               &
            , SurfSW                                                   &
            , dt_s                                                     &
            , sms                                                      &
            , Temp                                                     &
            , Sali_depth                                               & 
            , CO2_watercolumn                                          & 
            , pH_watercolumn                                           &
            , pCO2_watercolumn                                         & 
            , HCO3_watercolumn                                         & 
            , CO3_watercolumn                                          & 
            , OmegaC_watercolumn                                       & 
            , kspc_watercolumn                                         & 
            , rhoSW_watercolumn                                        &
            , Nutlim_phy                                               & 
            , Nutlim_dia                                               & 
            , Nutlim_cocco                                             &
            , Tlim_arr                                                 & 
            , Tlim_cocco                                               & 
            , Llim_phy                                                 & 
            , Llim_dia                                                 & 
            , Llim_cocco                                               & 
            , CO2lim_phy                                               & 
            , CO2lim_dia                                               & 
            , CO2lim_cocco                                             & 
            , PR_phy                                                   &
            , PR_dia                                                   & 
            , PR_cocco                                                 & 
            , Cal_Tlim                                                 & 
            , Cal_CO2lim                                               & 
            , Cal_Nlim                                                 & 
            , Cal_pure                                                 & 
            , SinkVel                                                  &
            , zF                                                       &
            , PAR                                                      &
            , Latd                                                     &
            , Lond                                                     & 
            , recom_istep                                              & 
            , mocsy_restart                                            & 
            , mocsy_step_per_day                                       & 
            , Loc_slp                                                  & 
            , daynew                                                   &                           
            , Latr)                                                      

  Use REcoM_declarations
  use REcoM_LocVar
  use recom_config
  use mvars                                                           
  use recom_ocean_settings                                            

  Implicit none

  Real(kind=8),                      intent(in)    :: cosAI                ! cos of angle of incidence
  Integer,                           intent(in)    :: Nn                   ! Total number of nodes in the vertical
  Real(kind=8),dimension(Nn,bgc_num),intent(inout) :: state                ! ChlA conc in phytoplankton [mg/m3]
  Real(kind=8),dimension(Nn),        intent(in)    :: dzF                  ! [m] Thickness of control volumes
  Real(kind=8),dimension(Nn),        intent(in)    :: recipdzF             ! [1/m] 
  Real(kind=8),dimension(Nn-1),      intent(in)    :: thick                ! [m] Vertical distance between nodes = Thickness
  Real(kind=8),dimension(Nn-1),      intent(in)    :: recipthick           ! [1/m] reciprocal of thick
  Real(kind=8),                      intent(in)    :: SurfSW               ! [W/m2] ShortWave radiation at surface
  Real(kind=8),                      intent(in)    :: Loc_slp              ! [Pa] sea-level pressure
  Real(kind=8),                      intent(in)    :: dt_s                 ! [s] Size of time steps
  Real(kind=8),dimension(Nn,bgc_num),intent(inout) :: sms                  ! Source-Minus-Sinks term
  Real(kind=8),dimension(Nn)        ,intent(in)    :: Temp                 ! [degrees C] Ocean temperature
  Real(kind=8),dimension(Nn)        ,intent(in)    :: Sali_depth           ! Salinity for the entire water column
  Real(kind=8),dimension(Nn+1,5)    ,intent(in)    :: SinkVel              ! 
  Real(kind=8),dimension(Nn+1)      ,intent(in)    :: zF                   ! [m] Depth of fluxes
  Real(kind=8),dimension(Nn),intent(inout)         :: PAR
  Real(kind=8),                      intent(in)    :: Latd(1)              ! latitude in degree
  Real(kind=8),                      intent(in)    :: Lond(1)              ! longitude in degree 
  Real(kind=8)                                     :: dt                   ! Size of time steps [day]
  Real(kind=8),dimension(Nn)                       :: Sink
  Real(kind=8)                                     :: dt_sink              ! Size of local time step
  Real(kind=8)                                     :: Fc                   ! Flux of labile C into sediment, used for denitrification calculation [umolC/cm2/s]
  Real(kind=8)                                     :: recip_hetN_plus      ! MB's addition to heterotrophic respiration
  Real(kind=8)                                     :: recip_res_het        ! [day] Reciprocal of respiration by heterotrophs and mortality (loss to detritus)
!  Real(kind=8)                                     :: Denit                ! Fraction of org. matter degradation occurring due to denitrification, dimensionless
  Integer                                          :: k,step,ii
  Integer, INTENT(IN)                              :: daynew 
  Integer, INTENT(IN)                              :: recom_istep         
  logical                                          :: mocsy_restart       
  Integer, INTENT(IN)                              :: mocsy_step_per_day 
  Real(kind=8),intent(in)                          :: Latr
  Real(kind=8)                                     :: REcoM_T_depth(1)        ! temperature for the whole water column for mocsy minimum defined as -2
  Real(kind=8)                                     :: REcoM_S_depth(1)       
  Real(kind=8)                                     :: REcoM_DIC_depth(1)   
  Real(kind=8)                                     :: REcoM_Alk_depth(1)    
  Real(kind=8)                                     :: REcoM_Si_depth(1)     
  Real(kind=8)                                     :: REcoM_Phos_depth(1)   
  Real(kind=8),dimension(Nn),intent(inout)         :: CO2_watercolumn   
  Real(kind=8),dimension(Nn),intent(inout)         :: pH_watercolumn    
  Real(kind=8),dimension(Nn),intent(inout)         :: pCO2_watercolumn   
  Real(kind=8),dimension(Nn),intent(inout)         :: HCO3_watercolumn  
  Real(kind=8),dimension(Nn),intent(inout)         :: CO3_watercolumn    
  Real(kind=8),dimension(Nn),intent(inout)         :: OmegaC_watercolumn 
  Real(kind=8),dimension(Nn),intent(inout)         :: kspc_watercolumn   
  Real(kind=8),dimension(Nn),intent(inout)         :: rhoSW_watercolumn  
  Real(kind=8),dimension(Nn),intent(inout)         :: Nutlim_phy      
  Real(kind=8),dimension(Nn),intent(inout)         :: Nutlim_dia      
  Real(kind=8),dimension(Nn),intent(inout)         :: Nutlim_cocco   
  Real(kind=8),dimension(Nn),intent(inout)         :: Tlim_arr        
  Real(kind=8),dimension(Nn),intent(inout)         :: Tlim_cocco       
  Real(kind=8),dimension(Nn),intent(inout)         :: Llim_phy          
  Real(kind=8),dimension(Nn),intent(inout)         :: Llim_dia        
  Real(kind=8),dimension(Nn),intent(inout)         :: Llim_cocco      
  Real(kind=8),dimension(Nn),intent(inout)         :: CO2lim_phy      
  Real(kind=8),dimension(Nn),intent(inout)         :: CO2lim_dia     
  Real(kind=8),dimension(Nn),intent(inout)         :: CO2lim_cocco   
  Real(kind=8),dimension(Nn),intent(inout)         :: PR_phy           
  Real(kind=8),dimension(Nn),intent(inout)         :: PR_dia          
  Real(kind=8),dimension(Nn),intent(inout)         :: PR_cocco        
  Real(kind=8),dimension(Nn),intent(inout)         :: Cal_Tlim        
  Real(kind=8),dimension(Nn),intent(inout)         :: Cal_CO2lim      
  Real(kind=8),dimension(Nn),intent(inout)         :: Cal_Nlim       
  Real(kind=8),dimension(Nn),intent(inout)         :: Cal_pure       
  Real(kind=8)                                     :: Patm_depth(1)   

  Real(kind=8)                                     :: & 
    DIN,     & ! [mmol/m3] Dissolved Inorganic Nitrogen 	
    DIC,     & ! [mmol/m3] Dissolved Inorganic Carbon
    Alk,     & !? [mmol/m3] Total Alkalinity
    PhyN,    & ! [mmol/m3] Intracellular conc of Nitrogen in small phytoplankton
    PhyC,    & ! [mmol/m3] Intracellular conc of Carbon in small phytoplankton
    Chl,     & ! [mg/m3] Current intracellular ChlA conc.
    DetN,    & ! [mmol/m3] Conc of N in Detritus
    DetC,    & ! [mmol/m3] Conc of C in Detritus
    HetN,    & ! [mmol/m3] Conc of N in heterotrophs
    HetC,    & ! [mmol/m3] Conc of C in heterotrophs
    DON,     & ! [mmol/m3] Dissolved organic N in the water
    EOC,     & !? [mmol/m3] Extracellular Organic C conc
    DiaN,    &
    DiaC,    &
    DiaChl,  &
    DiaSi,   &
    CoccoN,  & 
    CoccoC,  & 
    CoccoChl,& 
    DetSi,   &
    Si,      &
    Fe,      &
    PhyCalc, &
    DetCalc, &
    O2,      &
!    if (REcoM_Second_Zoo .eq. true) then
     Zoo2N,    &
     Zoo2C,    &
     DetZ2N,   &
     DetZ2C,   &
     DetZ2Si,  &
     DetZ2Calc,&
!    endif 
     FreeFe
     	
  sms = zero
  tiny_N   = tiny_chl/chl2N_max
  tiny_N_d = tiny_chl/chl2N_max_d
  tiny_N_c = tiny_chl/chl2N_max_c
  tiny_C   = tiny_N  /NCmax    
  tiny_C_d = tiny_N_d/NCmax_d
  tiny_C_c = tiny_N_c/NCmax_c   
  tiny_Si  = tiny_C_d/SiCmax


  recip_res_het = 1./res_het

  Patm_depth       = Loc_slp/Pa2atm             ! convert from Pa to atm.
  
!-------------------------------------------------------------------------------
! Size of REcoM time steps are calculated [day]
!-------------------------------------------------------------------------------

  rTref =  real(one)/recom_Tref
  
  dt	 =	dt_s/SecondsPerDay ! Size of fysics time step [day]
  dt	 =	dt/real(biostep)   ! Size of REcoM time step [day]

!-------------------------------------------------------------------------------
!Main time loop starts
  do step  = one,biostep

    kdzUpper	= 0.d0	            ! Upper light attenuation of top cell is set to zero
    do k  = one,Nn
      do ii = one,bgc_num
        if (abs(sms(k,ii)) .le. tiny) sms(k,ii) = zero
      end do
    end do 

!-------------------------------------------------------------------------------
! Main vertical loop starts
  do k = one,Nn

    DIN    = max(tiny,state(k,idin)   + sms(k,idin  )) ! Avoids division by zero
    DIC    = max(tiny,state(k,idic)   + sms(k,idic  )) ! and updates Conc between
    ALK    = max(tiny,state(k,ialk)   + sms(k,ialk  )) ! local steps in REcoM when
    PhyN   = max(tiny_N,state(k,iphyn)  + sms(k,iphyn )) ! biostep > 1
    PhyC   = max(tiny_C,state(k,iphyc)  + sms(k,iphyc ))
    Chl    = max(tiny_chl,state(k,ipchl)  + sms(k,ipchl ))
    DetN   = max(tiny,state(k,idetn)  + sms(k,idetn ))
    DetC   = max(tiny,state(k,idetc)  + sms(k,idetc ))
    HetN   = max(tiny,state(k,ihetn)  + sms(k,ihetn ))
    HetC   = max(tiny,state(k,ihetc)  + sms(k,ihetc ))
    if (REcoM_Second_Zoo) then 
     Zoo2N  = max(tiny,state(k,izoo2n)  + sms(k,izoo2n ))
     Zoo2C  = max(tiny,state(k,izoo2c)  + sms(k,izoo2c ))
     DetZ2N = max(tiny,state(k,idetz2n)  + sms(k,idetz2n ))
     DetZ2C = max(tiny,state(k,idetz2c)  + sms(k,idetz2c ))
     DetZ2Si = max(tiny,state(k,idetz2si)  + sms(k,idetz2si ))
     DetZ2Calc = max(tiny,state(k,idetz2calc)  + sms(k,idetz2calc ))
    endif
    DON    = max(tiny,state(k,idon)   + sms(k,idon  ))
    EOC    = max(tiny,state(k,idoc)   + sms(k,idoc  ))
    DiaN   = max(tiny_N,state(k,idian)  + sms(k,idian ))
    DiaC   = max(tiny_C,state(k,idiac)  + sms(k,idiac ))
    DiaChl = max(tiny_chl,state(k,idchl)  + sms(k,idchl ))
    DiaSi  = max(tiny_si,state(k,idiasi) + sms(k,idiasi))
    DetSi  = max(tiny,state(k,idetsi) + sms(k,idetsi))
    Si     = max(tiny,state(k,isi)    + sms(k,isi   ))
    CoccoN = max(tiny,state(k,icocn) + sms(k,icocn  ))
    CoccoC = max(tiny,state(k,icocc) + sms(k,icocc  ))
    CoccoChl = max(tiny,state(k,icchl) + sms(k,icchl))
    Fe     = max(tiny,state(k,ife)    + sms(k,ife   ))
    O2     = max(tiny,state(k,ioxy)    + sms(k,ioxy   ))
    FreeFe = zero
! For Mocsy
    REcoM_T_depth = max(2.d0, Temp(k))   ! minimum set to 2 degC: K1/K2 Lueker valid between 2degC-35degC and 19-43psu
    REcoM_T_depth = min(REcoM_T_depth, 40.d0)! maximum set to 40 degC: K1/K2 Lueker valid between 2degC-35degC and 19-43psu
    REcoM_S_depth = max(21.d0, Sali_depth(k))! minimum set to 21: K1/K2 Lueker valid between 2degC-35degC and 19-43psu, else causes trouble in regions with S between 19 and 21 and ice conc above 97%                      
    REcoM_S_depth    = min(REcoM_S_depth, 43.d0)  ! maximum set to 43: K1/K2 Lueker valid between 2degC-35degC and 19-43psu, else causes trouble   
    REcoM_DIC_depth  = max(tiny*1e-3,state(k,idic)*1e-3   + sms(k,idic  )*1e-3)    
    REcoM_Alk_depth  = max(tiny*1e-3,state(k,ialk)*1e-3   + sms(k,ialk  )*1e-3)     
    REcoM_Si_depth   = max(tiny*1e-3,state(k,isi)*1e-3    + sms(k,isi   )*1e-3)   
    REcoM_Phos_depth = max(tiny*1e-3,state(k,idin)*1e-3   + sms(k,idin  )*1e-3) /16  ! convert N to P with Redfield
!#ifdef REcoM_calcification		
    PhyCalc= max(tiny,state(k,iphycal)+ sms(k,iphycal))
    DetCalc= max(tiny,state(k,idetcal)+ sms(k,idetcal))
!#endif

    quota          =  PhyN / PhyC
    recipquota     =  real(one) / quota
    Chl2C          =  Chl  / PhyC
    Chl2N          =  Chl  / PhyN
    CHL2C_plast    =  Chl2C * (quota/(quota - NCmin))
    
    quota_dia      =  DiaN / DiaC
    recipQuota_dia =  real(one)/quota_dia
    Chl2C_dia      =  DiaChl / DiaC
    Chl2N_dia      =  DiaChl / DiaN
    CHL2C_plast_dia = Chl2C_dia * (quota_dia/(quota_dia - NCmin_d))
    qSiC           =  DiaSi / DiaC
    qSiN           =  DiaSi / DiaN

    quota_cocco    = CoccoN / CoccoC         
    recipQuota_cocco = real(one)/quota_cocco 
    Chl2C_cocco    = CoccoChl / CoccoC     
    Chl2N_cocco    = CoccoChl / CoccoN  
    CHL2C_plast_cocco = Chl2C_cocco * (quota_cocco/(quota_cocco - NCmin_c))

    recipQZoo      =  HetC / HetN
    recip_hetN_plus = 1. / (hetN + tiny_het) ! MB's addition for more stable zoo respiration
    if (REcoM_Second_Zoo) then  
      recipQZoo2     =  Zoo2C / Zoo2N
    endif
    if (Grazing_detritus) then
      recipDet = DetC / DetN
      recipDet2 = DetZ2C / DetZ2N
    end if

!-------------------------------------------------------------------------------
! Temperature dependence of rates
!-------------------------------------------------------------------------------
    rTloc          = real(one)/(Temp(k) + C2K)
    arrFunc        = exp( -Ae * ( rTloc - rTref))
    CoccoTFunc     = max(0.1419d0 * Temp(k)**0.8151d0,tiny)        ! (function from Fielding 2013; is based on observational GR, but range fits best to our arrFunc; they use T in degree celsius, so C2K is not needed)
    if (REcoM_Second_Zoo) then 
      arrFuncZoo2   = exp(t1_zoo2/t2_zoo2 - t1_zoo2*rTloc)/(1 + exp(t3_zoo2/t4_zoo2 - t3_zoo2*rTloc))
    endif
!-Silicate temperature dependence
    reminSiT = max(0.023d0 * 2.6d0**((Temp(k)-10.)/10.),reminSi)

    Tlim_arr(k)    = arrFunc                              
    Tlim_cocco(k)  = CoccoTFunc                         

!-------------------------------------------------------------------------------
! Light                                                                               
!-------------------------------------------------------------------------------
! Has to be calculated here already to use the 1%PAR depth.

    if (k==1) then
       PARave    = max(tiny,SurfSW)
       PAR(k)    = PARave
       chl_upper = (Chl + DiaChl + CoccoChl)                                                                                                                                                                              

    else
       chl_lower = Chl + DiaChl + CoccoChl              
       Chlave    = (chl_upper+chl_lower)*0.5

       kappa          =  k_w + a_chl * (Chlave)
       kappastar      =  kappa / cosAI
       kdzLower       =  kdzUpper + kappastar * thick(k)
       Lowerlight     =  SurfSW * exp(-kdzLower)
       Lowerlight     =  max(tiny,Lowerlight) 
       PARave         =  Lowerlight
       PAR(k)         =  PARave
       chl_upper      =  chl_lower
       kdzUpper       =  kdzLower
    end if

!-------------------------------------------------------------------------------
! Depth component of Mocsy (see http://ocmip5.ipsl.jussieu.fr/mocsy/pyth.html)                                                                 
!-------------------------------------------------------------------------------

! Calculate the carbonate system for the very first time step of the first year of the run
    if (recom_istep==1) then
       call vars_sprac(ph_depth, pco2_depth, fco2_depth, co2_depth, hco3_depth, co3_depth, OmegaA_depth, OmegaC_depth, kspc_depth, BetaD_depth, rhoSW_depth, p_depth, tempis_depth, & 
                            REcoM_T_depth, REcoM_S_depth, REcoM_Alk_depth, REcoM_DIC_depth, REcoM_Si_depth, REcoM_Phos_depth, Patm_depth, zF(k), Latd, Nmocsy,                      & 
                            optCON='mol/m3', optT='Tpot   ', optP='m ', optB='u74', optK1K2='l  ', optKf='dg', optGAS='Pinsitu', optS='Sprc')                         
       CO2_watercolumn(k)    = co2_depth(1)
       pH_watercolumn(k)     = ph_depth(1)
       pCO2_watercolumn(k)   = pco2_depth(1)
       HCO3_watercolumn(k)   = hco3_depth(1)
       CO3_watercolumn(k)    = co3_depth(1)  
       OmegaC_watercolumn(k) = OmegaC_depth(1)
       kspc_watercolumn(k)   = kspc_depth(1)  
       rhoSW_watercolumn(k)  = rhoSW_depth(1)
    endif 


! Calculate carbonate system every 7 days for depths < 1%PAR, and every 30 days for the depths below.
    logfile_outfreq_7  = mocsy_step_per_day*7                    ! mocsy_step_per_day seems to be defined in the fesom namelist.
    logfile_outfreq_30 = mocsy_step_per_day*30

    if (PARave > 0.01*SurfSW .and. mod(recom_istep,logfile_outfreq_7)==0) then
       call vars_sprac(ph_depth, pco2_depth, fco2_depth, co2_depth, hco3_depth, co3_depth, OmegaA_depth, OmegaC_depth, kspc_depth, BetaD_depth, rhoSW_depth, p_depth, tempis_depth, & 
                            REcoM_T_depth, REcoM_S_depth, REcoM_Alk_depth, REcoM_DIC_depth, REcoM_Si_depth, REcoM_Phos_depth, Patm_depth, zF(k), Latd, Nmocsy,                      & 
                            optCON='mol/m3', optT='Tpot   ', optP='m ', optB='u74', optK1K2='l  ', optKf='dg', optGAS='Pinsitu', optS='Sprc')                          
       CO2_watercolumn(k)    = co2_depth(1)
       pH_watercolumn(k)     = ph_depth(1)
       pCO2_watercolumn(k)   = pco2_depth(1)
       HCO3_watercolumn(k)   = hco3_depth(1)
       CO3_watercolumn(k)    = co3_depth(1)  
       OmegaC_watercolumn(k) = OmegaC_depth(1)
       kspc_watercolumn(k)   = kspc_depth(1)  
       rhoSW_watercolumn(k)  = rhoSW_depth(1) 
    elseif (PARave < 0.01*SurfSW .and. mod(recom_istep,logfile_outfreq_30)==0) then
       call vars_sprac(ph_depth, pco2_depth, fco2_depth, co2_depth, hco3_depth, co3_depth, OmegaA_depth, OmegaC_depth, kspc_depth, BetaD_depth, rhoSW_depth, p_depth, tempis_depth, & 
                           REcoM_T_depth, REcoM_S_depth, REcoM_Alk_depth, REcoM_DIC_depth, REcoM_Si_depth, REcoM_Phos_depth, Patm_depth, zF(k), Latd, Nmocsy,                       & 
                           optCON='mol/m3', optT='Tpot   ', optP='m ', optB='u74', optK1K2='l  ', optKf='dg', optGAS='Pinsitu', optS='Sprc')                          
       CO2_watercolumn(k)    = co2_depth(1)
       pH_watercolumn(k)     = ph_depth(1)
       pCO2_watercolumn(k)   = pco2_depth(1)
       HCO3_watercolumn(k)   = hco3_depth(1)
       CO3_watercolumn(k)    = co3_depth(1)  
       OmegaC_watercolumn(k) = OmegaC_depth(1) 
       kspc_watercolumn(k)   = kspc_depth(1) 
       rhoSW_watercolumn(k)  = rhoSW_depth(1)
    endif
    

!-------------------------------------------------------------------------------
! CO2 dependence of rates
!-------------------------------------------------------------------------------
! Convert pH to proton concentration
   h_depth(1) = 10.**(-ph_depth(1))

! Conversion factor between [mol/m3] (model) and [umol/kg] (function): (1000 * 1000) / 1024
   Cunits = 976.5625
   ! Conversion not needed for [H], because in model and function derived from pH and therefore in [mol/L]

! Small phytoplankton
   PhyCO2 = 1.162e+00 * HCO3_watercolumn(k) * Cunits / (4.888e+01 + HCO3_watercolumn(k) * Cunits) - exp(-2.255e-01 * CO2_watercolumn(k) * Cunits) - 1.023e+07 * 10.**(-pH_watercolumn(k))
   PhyCO2 = min(PhyCO2,3.0)
   if((co2flux(1)>1.e10) .or. (co2flux(1)<-1.e10)) then
     write(*,*), 'PhyCO2 =',PhyCO2
   endif 
   CO2lim_phy(k)   = PhyCO2 

! Diatoms
   DiaCO2 = 1.040e+00 * HCO3_watercolumn(k) * Cunits / (2.890e+01 + HCO3_watercolumn(k) * Cunits) - exp(-8.778e-01 * CO2_watercolumn(k) * Cunits) - 2.640e+06 * 10.**(-pH_watercolumn(k))
   DiaCO2 = min(DiaCO2,3.0) 
   if((co2flux(1)>1.e10) .or. (co2flux(1)<-1.e10)) then
      write(*,*), 'DiaCO2 =',DiaCO2
   endif
   CO2lim_dia(k)   = DiaCO2 

! Coccolithophores
   CoccoCO2 = 1.109e+00 * HCO3_watercolumn(k) * Cunits / (3.767e+01 + HCO3_watercolumn(k) * Cunits) - exp(-3.912e-01 * CO2_watercolumn(k) * Cunits) - 9.450e+06 * 10.**(-pH_watercolumn(k))
   CoccoCO2 = min(CoccoCO2,3.0) 
   if((co2flux(1)>1.e10) .or. (co2flux(1)<-1.e10)) then
      write(*,*), 'CoccoCO2 =',CoccoCO2
   endif
   CO2lim_cocco(k) = CoccoCO2


!------------------------------------------------------------------------------
! Calcite dissolution dependent on OmegaC 
!------------------------------------------------------------------------------

   if (OmegaC_diss) then
    Ca                 = (0.02128d0/40.078d0) * Sali_depth(k)/1.80655d0 ! Calcium ion concentration [mol/kg], function from varsolver.f90
    CO3_sat            = (kspc_watercolumn(k) / Ca) * rhoSW_watercolumn(k) ! Saturated carbonate ion concentration, converted to [mol/m3]
    calc_diss          = calc_diss_omegac * max(zero,(1-(CO3_watercolumn(k)/CO3_sat)))**(calc_diss_exp) ! Dissolution rate scaled by carbonate ratio, after Aumont et al. 2015
   else
    calc_diss      = calc_diss_rate * SinkVel(k,ivdet) /20.d0 ! Dissolution rate of CaCO3 scaled by the sinking velocity at the current depth
   endif
    calc_diss2     = calc_diss_rate2  ! Dissolution rate of CaCO3 scaled by the sinking velocity at the current depth  seczoo

   calc_diss_ben   = calc_diss_rate * SinkVel(k,ivdet) /20.0

!-------------------------------------------------------------------------------
! Photosynthesis section, light parameters and rates
!-------------------------------------------------------------------------------
! Small phytoplankton

    qlimitFac     = recom_limiter(NMinSlope,NCmin,quota)
      feLimitFac  = Fe/(k_Fe + Fe)
      qlimitFac   = min(qlimitFac,feLimitFac)
    Nutlim_phy(k) = qlimitFac                             

    pMax          = P_cm * qlimitFac * arrFunc

!-------------------------------------------------------------------------------
! Diatoms

    qlimitFac     = recom_limiter(NMinSlope,NCmin_d,quota_dia)
    qlimitFacTmp  = recom_limiter(SiMinSlope,SiCmin,qSiC)
    qlimitFac     = min(qLimitFac,qlimitFacTmp)
      feLimitFac  = Fe/(k_Fe_d + Fe)
      qlimitFac   = min(qlimitFac,feLimitFac)
    Nutlim_dia(k) = qlimitFac                           

    pMax_dia      = P_cm_d * qlimitFac * arrFunc

!------------------------------------------------------------------------------
! Coccolithophores

    qlimitFac     = recom_limiter(NMinSlope,NCmin_c,quota_cocco) 
      feLimitFac  = Fe/(k_Fe_c + Fe)                            
      qlimitFac   = min(qlimitFac,feLimitFac)                  
    Nutlim_cocco(k)=qlimitFac                                   

    pMax_cocco    = P_cm_c * qlimitFac * CoccoTFunc    

!-------------------------------------------------------------------------------
! Small phytoplankton photosynthesis rate

    if ( pMax .lt. tiny .OR. PARave /= PARave                  &
         .OR. CHL2C /= CHL2C) then
      Cphot       = zero
      Llim_phy(k) = zero                                            
    else
      if (CO2lim) then
        Cphot     = pMax*( real(one) &
                   - exp(-alfa * Chl2C * PARave / pMax)) * PhyCO2    
      else
        Cphot     = pMax*( real(one) &
                   - exp(-alfa * Chl2C * PARave / pMax))
      endif
      Llim_phy(k) = real(one) - exp(-alfa * Chl2C * PARave / pMax)  
    end if
    if ( Cphot .lt. tiny) Cphot = zero

    PR_phy(k)     = Cphot                                         
    
 !------------------------------------------------------------------------------
 ! Diatom photosynthesis rate

    if ( pMax_dia .lt. tiny .OR. PARave /= PARave               &
         .OR. CHL2C_dia /= CHL2C_dia) then
      Cphot_dia   = zero
      Llim_dia(k) = zero                                          
    else
      if (CO2lim) then
        Cphot_dia = pMax_dia * (real(one) &
           	- exp( -alfa_d * Chl2C_dia * PARave / pMax_dia)) * DiaCO2 
      else
        Cphot_dia = pMax_dia * (real(one) &
                - exp( -alfa_d * Chl2C_dia * PARave / pMax_dia))
      endif
      Llim_dia(k) = real(one)-exp(-alfa_d*Chl2C_dia*PARave/pMax_dia) 
    end if
    if (Cphot_dia .lt. tiny) Cphot_dia = zero

    PR_dia(k)     = Cphot_dia                                      

!--------------------------------------------------------------------------------
! Coccolithophore photosynthesis rate

    if ( pMax_cocco .lt. tiny .OR. Parave /= Parave             &    
         .OR. CHL2C_cocco /= CHL2C_cocco) then                 
       Cphot_cocco = zero                                     
       Llim_cocco(k)=zero                                      
    else                                                        
       if (CO2lim) then
         Cphot_cocco = pMax_cocco * (real(one) &              
                 - exp( -alfa_c * Chl2C_cocco * PARave / pMax_cocco)) * CoccoCO2 
       else
         Cphot_cocco = pMax_cocco * (real(one) &                 
                 - exp( -alfa_c * Chl2C_cocco * PARave / pMax_cocco)) 
       endif
       Llim_cocco(k)=real(one)-exp(-alfa_c*Chl2C_cocco*PARave/pMax_cocco)
    end if                                                        
    if (Cphot_cocco .lt. tiny) Cphot_cocco = zero                

    PR_cocco(k)    = Cphot_cocco                                

!-------------------------------------------------------------------------------- 
! chlorophyll degradation

    KOchl = deg_Chl
    KOchl_dia = deg_Chl_d
    KOchl_cocco = deg_Chl_c                       
        
    if (use_photodamage) then
!    Phyto chla loss                                                                                                                                                      
! adding a minimum value for photodamage 
         if (pMax .lt. tiny .OR. PARave /= PARave                 &
             .OR. CHL2C_plast /= CHL2C_plast) then
           KOchl = deg_Chl*0.1d0
        else
           KOchl = deg_Chl*(1-exp(-alfa * CHL2C_plast * PARave / pMax))
           KOchl = max((deg_Chl*0.1d0), KOchl)
        end if

!    Diatoms chla loss                                                                                                                                                     
        if (pMax_dia .lt. tiny .OR. PARave /= PARave             &
                 .OR. CHL2C_plast_dia /= CHL2C_plast_dia) then
           KOchl_dia = deg_Chl_d*0.1d0
        else
           KOchl_dia = deg_Chl_d * ( 1 -                         &
                exp( -alfa_d * CHL2C_plast_dia * PARave / pMax_dia ))
           KOchl_dia = max((deg_Chl_d*0.1d0), KOchl_dia)
        end if

!    Coccolithophores chla loss 
        if (pMax_cocco .lt. tiny .OR. PARave /= Parave           &  
                 .OR. CHL2C_plast_cocco /= CHL2C_plast_cocco) then   
           KOchl_cocco = deg_Chl_c*0.1d0                          
        else                                                      
           KOchl_cocco = deg_Chl_c * ( 1 -                       &   
                exp( -alfa_c * CHL2C_plast_cocco * PARave / pMax_cocco ))
           KOchl_cocco = max((deg_Chl_c*0.1d0), KOchl_cocco)   
        end if

    if (KOchl /= KOchl) then
       print*,' KOchl is ', KOchl
       print*,' deg_Chl is ', deg_Chl
       print*,' alfa is ', alfa
       print*,' CHL2C is ', CHL2C_plast
       print*,' PARave is ', PARave
       print*,' pMax is ', pMax
       stop
    end if
    if (KOchl_dia /= KOchl_dia) then
       print*,' KOchl_dia is ', KOchl_dia
       print*,' deg_Chl_d is ', deg_Chl_d
       print*,' alfa_d is ', alfa_d
       print*,' CHL2C_d is ', CHL2C_plast_dia
       print*,' PARave is ', PARave
       print*,' pMax_d is ', pMax_dia
       stop
    end if
    if (KOchl_cocco /= KOchl_cocco) then                         
       print*,' KOchl_cocco is ', KOchl_cocco                  
       print*,' deg_Chl_c is ', deg_Chl_c                       
       print*,' alfa_c is ', alfa_c                             
       print*,' CHL2C_c is ', CHL2C_plast_cocco                   
       print*,' PARave is ', PARave                              
       print*,' pMax_c is ', pMax_cocco                            
       stop                                                        
    end if                                                    

   end if ! photodamage  
    
!-------------------------------------------------------------------------------
! Assimilation section
!-------------------------------------------------------------------------------
! Compute assimilation from Geider et al 1998

    V_cm           = V_cm_fact
    limitFacN      = recom_limiter(NMaxSlope,quota,NCmax)
    N_assim        = V_cm * pMax * NCuptakeRatio &                ! [mmol N / (mmol C * day)]
                      * limitFacN * (DIN/(DIN + k_din))

    V_cm           = V_cm_fact_d
    limitFacN_dia  = recom_limiter(NMaxSlope,quota_dia,NCmax_d)
    N_assim_dia    = V_cm * pMax_dia * NCUptakeRatio_d &
                      * limitFacN_dia * DIN/(DIN + k_din_d)

    V_cm           = V_cm_fact_c                                    
    limitFacN_cocco= recom_limiter(NMaxSlope,quota_cocco,NCmax_c)   
    N_assim_cocco  = V_cm * pMax_cocco * NCUptakeRatio_c &         
                      * limitFacN_cocco * DIN/(DIN + k_din_c)      

    limitFacSi     = recom_limiter(SiMaxSlope,qSiC,SiCmax)  &
                      * limitFacN_dia
    Si_assim       = V_cm * P_cm_d * arrFunc * SiCUptakeRatio &
                      * limitFacSi * Si/(Si + k_si)

!-------------------------------------------------------------------------------
! Iron chemistry
 
      freeFe      = iron_chemistry(Fe,totalligand,ligandStabConst)

!-------------------------------------------------------------------------------
! Chlorophyll synthesis

    chlSynth       = zero
    if (PARave .ge. tiny .AND. PARave .eq. PARave) then
       chlSynth = N_assim * Chl2N_max                    &
          * min( real(one),Cphot/(alfa*Chl2C*PARave))
    end if
    ChlSynth_dia   = zero
    if (PARave .ge. tiny .AND. PARave .eq. PARave) then
       ChlSynth_dia = N_assim_dia                        &
          * Chl2N_max_d * min(real(one),                 &
            Cphot_dia /(alfa_d * Chl2C_dia * PARave))
    end if
    ChlSynth_cocco    = zero                                        
    if (PARave .ge. tiny .AND. PARave .eq. PARave) then            
       ChlSynth_cocco = N_assim_cocco                    &         
          * Chl2N_max_c * min(real(one),                 &         
            Cphot_cocco /(alfa_c * Chl2C_cocco * PARave))          
    end if                                                        

!-------------------------------------------------------------------------------
! Phytoplankton respiration rate

    phyRespRate       = res_phy * limitFacN + biosynth * N_assim
    phyRespRate_dia   = res_phy_d * limitFacN_dia +          &
            biosynth * N_assim_dia + biosynthSi * Si_assim
    phyRespRate_cocco = res_phy_c * limitFacN_cocco +        &    
            biosynth * N_assim_cocco                              
!-------------------------------------------------------------------------------
! Zooplankton grazing on small phytoplankton and diatoms and coccolithophores!

     if (REcoM_Grazing_Variable_Preference) then
       if (Grazing_detritus) then
         if (Graz_pref_new) then
          varpzPhy      = pzPhy * PhyN /(pzPhy*PhyN + pzDia*DiaN + pzCocco*CoccoN + PzDet*DetN + pzDetZ2*DetZ2N)    
          varpzDia      = pzDia * DiaN /(pzPhy*PhyN + pzDia*DiaN + pzCocco*CoccoN + PzDet*DetN + pzDetZ2*DetZ2N)    
          varpzCocco    = pzCocco * CoccoN /(pzPhy*PhyN + pzDia*DiaN + pzCocco*CoccoN + PzDet*DetN + pzDetZ2*DetZ2N) 
          varpzDet      = pzDet * DetN /(pzPhy*PhyN + pzDia*DiaN + pzCocco*CoccoN + PzDet*DetN + pzDetZ2*DetZ2N)   
          varpzDetZ2    = pzDetZ2 * DetN /(pzPhy*PhyN + pzDia*DiaN + pzCocco*CoccoN + PzDet*DetN + pzDetZ2*DetZ2N)  
         else
          DiaNsq        = DiaN * DiaN
          varpzDia      = pzDia * DiaNsq /(sDiaNsq + DiaNsq)
          PhyNsq        = PhyN * PhyN
          varpzPhy      = pzPhy * PhyNsq /(sPhyNsq + PhyNsq)
          CoccoNsq      = CoccoN * CoccoN                             
          varpzCocco    = pzCocco * CoccoNsq /(sCoccoNsq + CoccoNsq)  
          DetNsq        = DetN * DetN
          varpzDet      = pzDet * DetNsq /(sDetNsq + DetNsq)
          DetZ2Nsq        = DetZ2N * DetZ2N
          varpzDetZ2      = pzDetZ2 * DetZ2Nsq /(sDetZ2Nsq + DetZ2Nsq)
         endif
         fDiaN         = varpzDia * DiaN
         fPhyN         = varpzPhy * PhyN
         fCoccoN       = varpzCocco * CoccoN                      
         fDetN         = varpzDet * DetN
         fDetZ2N       = varpzDetZ2 * DetZ2N 
       else
        if (Graz_pref_new) then
         varpzPhy      = pzPhy * PhyN /(pzPhy*PhyN + pzDia*DiaN + pzCocco*CoccoN)     
         varpzDia      = pzDia * DiaN /(pzPhy*PhyN + pzDia*DiaN + pzCocco*CoccoN)       
         varpzCocco    = pzCocco * CoccoN /(pzPhy*PhyN + pzDia*DiaN + pzCocco*CoccoN)  
        else
         DiaNsq        = DiaN * DiaN
         varpzDia      = pzDia * DiaNsq /(sDiaNsq + DiaNsq)
         PhyNsq        = PhyN * PhyN
         varpzPhy      = pzPhy * PhyNsq /(sPhyNsq + PhyNsq)
         CoccoNsq      = CoccoN * CoccoN                           
         varpzCocco    = pzCocco * CoccoNsq /(sCoccoNsq + CoccoNsq)  
        end if
         fDiaN         = varpzDia * DiaN
         fPhyN         = varpzPhy * PhyN
         fCoccoN       = varpzCocco * CoccoN                        
       end if
     else
       fDiaN         = pzDia * DiaN
       fPhyN         = pzPhy * PhyN
       fCoccoN       = pzCocco * CoccoN                          
       if (Grazing_detritus) then
        fDetN        = pzDet * DetN
        fDetZ2N      = pzDetZ2 * DetZ2N
       end if
     end if
     if (Grazing_detritus) then
       food            = fPhyN + fDiaN + fCoccoN + fDetN + fDetZ2N   
       foodsq          = food * food
       grazingFlux     = (Graz_max * foodsq)/(epsilon + foodsq) * HetN * arrFunc
       grazingFlux_phy = grazingFlux * fphyN / food
       grazingFlux_Dia = grazingFlux * fDiaN / food
       grazingFlux_Cocco = grazingFlux * fCoccoN / food           
       grazingFlux_Det = grazingFlux * fDetN / food
       grazingFlux_DetZ2 = grazingFlux * fDetZ2N / food
     else
       food            = fPhyN + fDiaN + fCoccoN                    
       foodsq          = food * food
       grazingFlux     = (Graz_max * foodsq)/(epsilon + foodsq) * HetN * arrFunc
       grazingFlux_phy = grazingFlux * fphyN / food
       grazingFlux_Dia = grazingFlux * fDiaN / food
       grazingFlux_Cocco = grazingFlux * fCoccoN / food         
     endif


     if (REcoM_Second_Zoo) then
!-------------------------------------------------------------------------------
! Second Zooplankton grazing on small phytoplankton, diatoms and heterotrophs and coccolithophores
        
      if (REcoM_Grazing_Variable_Preference) then
        if (Grazing_detritus) then
         if (Graz_pref_new) then
          varpzDia2      = pzDia2 * DiaN /(pzPhy2 * PhyN + PzDia2 * DiaN + pzCocco2 * CoccoN + pzHet * HetN + pzDet2 * DetN + pzDetZ22 * DetZ2N)     
          varpzPhy2      = pzPhy2 * PhyN /(pzPhy2 * PhyN + PzDia2 * DiaN + pzCocco2 * CoccoN + pzHet * HetN + pzDet2 * DetN + pzDetZ22 * DetZ2N)     
          varpzCocco2    = pzCocco2 * CoccoN /(pzPhy2 * PhyN + PzDia2 * DiaN + pzCocco2 * CoccoN + pzHet * HetN + pzDet2 * DetN + pzDetZ22 * DetZ2N) 
          varpzHet       = pzHet * HetN /(pzPhy2 * PhyN + PzDia2 * DiaN + pzCocco2 * CoccoN + pzHet * HetN + pzDet2 * DetN + pzDetZ22 * DetZ2N)      
          varpzDet2      = pzDet2 * DetN /(pzPhy2 * PhyN + PzDia2 * DiaN + pzCocco2 * CoccoN + pzHet * HetN + pzDet2 * DetN + pzDetZ22 * DetZ2N)      
          varpzDetZ22    = pzDetZ22 * DetZ2N /(pzPhy2 * PhyN + PzDia2 * DiaN + pzCocco2 * CoccoN + pzHet * HetN + pzDet2 * DetN + pzDetZ22 * DetZ2N)  
         else
          DiaNsq2        = DiaN * DiaN
          varpzDia2      = pzDia2 * DiaNsq2 /(sDiaNsq2 + DiaNsq2)
          fDiaN2         = varpzDia2 * DiaN
          PhyNsq2        = PhyN * PhyN
          varpzPhy2      = pzPhy2 * PhyNsq2 /(sPhyNsq2 + PhyNsq2)
          fPhyN2         = varpzPhy2 * PhyN
          varpzCocco2    = pzCocco2 * CoccoNsq2 /(sCoccoNsq2 + CoccoNsq2) 
          fCoccoN2       = varpzCocco2 * CoccoN                        
          HetNsq         = HetN * HetN
          varpzHet       = pzHet * HetNsq /(sHetNsq + HetNsq)
          fHetN          = varpzHet * HetN
          DetNsq         = DetN * DetN
          varpzDet2      = pzDet2 * DetNsq /(sDetNsq2 + DetNsq)
          DetZ2Nsq       = DetZ2N * DetZ2N
          varpzDetZ22    = pzDetZ22 * DetZ2Nsq /(sDetZ2Nsq2 + DetZ2Nsq)
         end if
          fDiaN2         = varpzDia2 * DiaN
          fPhyN2         = varpzPhy2 * PhyN
          fCoccoN2       = varpzCocco2 * CoccoN                       
          fHetN          = varpzHet * HetN
          fDetN2         = varpzDet2 * DetN
          fDetZ2N2       = varpzDetZ22 * DetZ2N
        else
          if (Graz_pref_new) then
           varpzDia2      = pzDia2 * DiaN /(pzPhy2 * PhyN + PzDia2 * DiaN + pzCocco2 * CoccoN + pzHet * HetN)      
           varpzPhy2      = pzPhy2 * PhyN /(pzPhy2 * PhyN + PzDia2 * DiaN + pzCocco2 * CoccoN + pzHet * HetN)      
           varpzCocco2    = pzCocco2 * CoccoN /(pzPhy2 * PhyN + PzDia2 * DiaN + pzCocco2 * CoccoN + pzHet * HetN)  
           varpzHet       = pzHet * HetN /(pzPhy2 * PhyN + PzDia2 * DiaN + pzCocco2 * CoccoN + pzHet * HetN)      
          else
           DiaNsq2        = DiaN * DiaN
           varpzDia2      = pzDia2 * DiaNsq2 /(sDiaNsq2 + DiaNsq2)
           fDiaN2         = varpzDia2 * DiaN
           PhyNsq2        = PhyN * PhyN
           varpzPhy2      = pzPhy2 * PhyNsq2 /(sPhyNsq2 + PhyNsq2)
           fPhyN2         = varpzPhy2 * PhyN
           CoccoNsq2      = CoccoN * CoccoN                            
           varpzCocco2    = pzCocco2 * CoccoNsq2 /(sCoccoNsq2 + CoccoNsq2)
           fCoccoN2       = varpzCocco2 * CoccoN                         
           HetNsq         = HetN * HetN
           varpzHet       = pzHet * HetNsq /(sHetNsq + HetNsq)
          end if
          fDiaN2         = varpzDia2 * DiaN
          fPhyN2         = varpzPhy2 * PhyN
          fCoccoN2       = varpzCocco2 * CoccoN                      
          fHetN          = varpzHet * HetN
        end if
      else
       fDiaN2         = pzDia2 * DiaN
       fPhyN2         = pzPhy2 * PhyN
       fCoccoN2       = pzCocco2 * CoccoN                             
       fHetN          = pzHet * HetN
       if (Grazing_detritus) then
        fDetN2        = pzDet2 * DetN
        fDetZ2N2      = pzDetZ22 * DetZ2N
       end if
      end if
     
      if (Grazing_detritus) then
        food2             = fPhyN2 + fDiaN2 + fCoccoN2 + fHetN + fDetN2 + fDetZ2N2   
        foodsq2           = food2 * food2
        grazingFlux2     = (Graz_max2 * foodsq2)/(epsilon2 + foodsq2) * Zoo2N * arrFuncZoo2
        grazingFlux_phy2  = grazingFlux2 * fphyN2 / food2
        grazingFlux_Dia2  = grazingFlux2 * fDiaN2 / food2
        grazingFlux_Cocco2 = grazingFlux2 * fCoccoN2 / food2                  
        grazingFlux_het2  = grazingFlux2 * fHetN / food2
        grazingFlux_Det2  = grazingFlux2 * fDetN2 / food2
        grazingFlux_DetZ22  = grazingFlux2 * fDetZ2N2 / food2

        grazingFluxcarbonzoo2 = (grazingFlux_phy2 * recipQuota * grazEff2) &
                          + (grazingFlux_Dia2 * recipQuota_Dia * grazEff2) &
                          + (grazingFlux_Cocco2 * recipQuota_Cocco * grazEff2) &    
                          + (grazingFlux_het2 * recipQZoo * grazEff2)      &
                          + (grazingFlux_Det2 * recipDet * grazEff2)       &
                          + (grazingFlux_DetZ22 *recipDet2 * grazEff2)    
      else
        food2             = fPhyN2 + fDiaN2 + fCoccoN2 + fHetN                  
        foodsq2           = food2 * food2
        grazingFlux2      = (Graz_max2 * foodsq2)/(epsilon2 + foodsq2) * Zoo2N * arrFuncZoo2
        grazingFlux_phy2  = grazingFlux2 * fphyN2 / food2
        grazingFlux_Dia2  = grazingFlux2 * fDiaN2 / food2
        grazingFlux_Cocco2 = grazingFlux2 * fCoccoN2 / food2                     
        grazingFlux_het2  = grazingFlux2 * fHetN / food2

        grazingFluxcarbonzoo2 = (grazingFlux_phy2 * recipQuota * grazEff2) &
                          + (grazingFlux_Dia2 * recipQuota_Dia * grazEff2) &
                          + (grazingFlux_Cocco2 * recipQuota_Cocco * grazEff2) &  
                          + (grazingFlux_het2 * recipQZoo * grazEff2)
      end if
     end if
!-------------------------------------------------------------------------------
! Heterotrophic respiration is assumed to drive zooplankton back to Redfield C:N
! if their C:N becomes higher than Redfield
    if (het_resp_noredfield) then
      HetRespFlux    = res_het *  arrFunc * HetC
    else
      if (HetRespFlux_plus) then
       HetRespFlux    = recip_res_het * arrFunc * (hetC    * recip_hetN_plus - redfield) * HetC
      else
! default computation scheme 
       HetRespFlux    = recip_res_het * arrFunc * (recipQZoo    - redfield) * HetC  
      endif
      HetRespFlux    = max(zero,HetRespFlux)
    endif
!-------------------------------------------------------------------------------
! Quadratic zooplanton mortality
    hetLossFlux    = loss_het * HetN * HetN

 if (REcoM_Second_Zoo) then
!-------------------------------------------------------------------------------
! if (REcoM_Second_Zoo .eq. .true.) then
! Second zooplankton respiration is assumed to drive zooplankton back to Redfield C:N
! if their C:N becomes higher than Redfield
  call krill_resp(daynew,Latr)
   if((grazingFluxcarbonzoo2/Zoo2C) <= 0.1)then
    res_zoo2_f = 0.1*(grazingFluxcarbonzoo2/Zoo2C*100)
     else
    res_zoo2_f = 1.
   end if  
   recip_res_zoo22 = res_zoo2*(1.+res_zoo2_f + res_zoo2_a)
   Zoo2RespFlux    = recip_res_zoo22 * Zoo2C                   
!-------------------------------------------------------------------------------
! Quadratic second zooplanton mortality
    Zoo2LossFlux    = loss_zoo2 * zoo2N * zoo2N

    if(zoo2_fecal_loss) then
    Zoo2fecalloss_n   = fecal_rate_n * grazingFlux2
    Zoo2fecalloss_c   = fecal_rate_c * grazingFluxcarbonzoo2
     else
    Zoo2fecalloss_n   = 0.0
    Zoo2fecalloss_c   = 0.0
    end if
 end if
!-------------------------------------------------------------------------------
! Phytoplankton and detritus aggregation
 if (diatom_mucus) then
  qlimitFac     = recom_limiter(NMinSlope,NCmin_d,quota_dia)
  qlimitFacTmp  = recom_limiter(SiMinSlope,SiCmin,qSiC)
  qlimitFac     = min(qLimitFac,qlimitFacTmp)
  feLimitFac  = Fe/(k_Fe_d + Fe)
  qlimitFac   = min(qlimitFac,feLimitFac)
  if (REcoM_Second_Zoo) then 
    aggregationrate = agg_PD * DetN + agg_PD * DetZ2N  &
                     + agg_PP * PhyN + agg_PP * CoccoN &       
                     + agg_PP * (1 - qlimitFac) * DiaN
  else
    aggregationrate = agg_PD * DetN + agg_PP * PhyN &
                    + agg_PP * CoccoN &                          
                    + agg_PP * (1 - qlimitFac) * DiaN
  endif
 else
  if (REcoM_Second_Zoo) then
    aggregationrate = agg_PD * DetN + agg_PD * DetZ2N &
                     + agg_PP * PhyN + agg_PP * CoccoN + agg_PP * DiaN 
  else
    aggregationrate = agg_PD * DetN + agg_PP * PhyN &
                    + agg_PP * CoccoN + agg_PP * DiaN        
  endif
 endif
!-------------------------------------------------------------------------------
! Terms required for the formation and dissolution of CaCO3
!#ifdef REcoM_calcification
    if (Temp(k) < 10.6) then                                                 ! PICPOC definition after Krumhardt et al. 2017
      PICPOCtemp = 0.104d0 * Temp(k) - 0.108d0
       else
      PICPOCtemp = 1
    end if
    PICPOCtemp   = max(tiny,PICPOCtemp)
    Cal_Tlim(k)  = PICPOCtemp                                     

    PICPOCCO2     = 1.102e+00 * HCO3_watercolumn(k) * Cunits / (4.238e+01 + HCO3_watercolumn(k) * Cunits) - exp(-7.079e-01 * CO2_watercolumn(k) * Cunits) - 1.343e+07 * 10.**(-pH_watercolumn(k))
    PICPOCCO2     = min(PICPOCCO2,3.0)                                    
    Cal_CO2lim(k) = PICPOCCO2                                      

    PICPOCN     = -0.31 * (DIN/(DIN + k_din_c)) + 1.31
    PICPOCN     = max(tiny,PICPOCN)
    Cal_Nlim(k) = PICPOCN                                          
    
    if (CO2lim) then
      calcification = 1.d0 * Cphot_cocco * CoccoC * PICPOCtemp * PICPOCN * PICPOCCO2 
    else
      calcification = 1.d0 * Cphot_cocco * CoccoC * PICPOCtemp * PICPOCN
    endif
    Cal_pure(k)   = calcification                                   
  
    calc_loss_agg = aggregationRate * PhyCalc                        

  if(REcoM_Second_Zoo)  then
     calc_loss_gra  =  grazingFlux_Cocco   &                            
                     * recipQuota_Cocco/(CoccoC + tiny)    * PhyCalc   
     calc_loss_gra2 = grazingFlux_Cocco2           &                   
                     * recipQuota_Cocco/(CoccoC + tiny)    * PhyCalc    
   else
    calc_loss_gra = grazingFlux_Cocco           &                        
                    * recipQuota_Cocco/(CoccoC + tiny)    * PhyCalc    
  endif  
!#endif
		
!-------------------------------------------------------------------------------
! Sources minus sinks are calculated
!-------------------------------------------------------------------------------
! DIN
    sms(k,idin)      = (                       &
      - N_assim                      * PhyC    &
      - N_assim_Dia                  * DiaC    &
      - N_assim_Cocco                * CoccoC  &    
      + rho_N * arrFunc              * DON     &
      + LocRiverDIN                            &
                                             ) * dt + sms(k,idin)  
!-------------------------------------------------------------------------------
! DIC
  if(REcoM_Second_Zoo)  then
    sms(k,idic)      = (                       &
     - Cphot                         * PhyC    &
     + phyRespRate                   * PhyC    &
     - Cphot_Dia                     * DiaC    &
     + phyRespRate_Dia               * DiaC    &
     - Cphot_Cocco                   * CoccoC  &   
     + phyRespRate_Cocco             * CoccoC  &       
     + rho_C1 * arrFunc              * EOC     &
     + HetRespFlux                             &                      
     + Zoo2RespFlux                            &
!#ifdef REcoM_calcification                     
     + calc_diss                     * DetCalc &
     + calc_loss_gra * calc_diss_guts          &
     + calc_loss_gra2 * calc_diss_guts         &
     + calc_diss2                     * DetZ2Calc &
     - calcification                           &
!#endif
                                             ) * dt + sms(k,idic)
  else
    sms(k,idic)      = (                       &
     - Cphot                         * PhyC    &
     + phyRespRate                   * PhyC    &
     - Cphot_Dia                     * DiaC    &
     + phyRespRate_Dia               * DiaC    &
     - Cphot_Cocco                   * CoccoC  &  
     + phyRespRate_Cocco             * CoccoC  &    
     + rho_C1 * arrFunc              * EOC     &
     + HetRespFlux                             & 
!#ifdef REcoM_calcification                     
     + calc_diss                     * DetCalc &
     + calc_loss_gra * calc_diss_guts          &
     - calcification                           &
!#endif
                                             ) * dt + sms(k,idic)
  endif

!-------------------------------------------------------------------------------
! Alkalinity (Assumes that N:P follows a constant Redfield ratio
! N_assimC: 1.0625 = 1/16 + 1
   if (REcoM_Second_Zoo) then
    sms(k,ialk)      = (                       &
      + 1.0625 * N_assim             * PhyC    &
      + 1.0625 * N_assim_Dia         * DiaC    &
      + 1.0625 * N_assim_Cocco       * CoccoC  &        
      - 1.0625 * rho_N * arrFunc     * DON     &
!#ifdef REcoM_calcification                                                                                                                                                    
                        
      + 2.d0 * calc_diss             * DetCalc &
      + 2.d0 * calc_loss_gra * calc_diss_guts  &
      + 2.d0 * calc_loss_gra2 * calc_diss_guts &
      + 2.d0 * calc_diss2             * DetZ2Calc &
      - 2.d0 * calcification                   &
!#endif                                                                                                                                                                        
                         
                                             ) * dt + sms(k,ialk)
   else
    sms(k,ialk)      = (                       &
      + 1.0625 * N_assim             * PhyC    &
      + 1.0625 * N_assim_Dia         * DiaC    &
      + 1.0625 * N_assim_Cocco       * CoccoC  &      
      - 1.0625 * rho_N * arrFunc     * DON     &
!#ifdef REcoM_calcification
      + 2.d0 * calc_diss             * DetCalc &
      + 2.d0 * calc_loss_gra * calc_diss_guts  &
      - 2.d0 * calcification                   &
!#endif 
                                             ) * dt + sms(k,ialk) 
   endif
!-------------------------------------------------------------------------------
! Phytoplankton N
   if (REcoM_Second_Zoo) then
     sms(k,iphyn)      = (                      &
      + N_assim                      * PhyC    &
      - lossN * limitFacN            * PhyN    &
      - aggregationRate              * phyN    &
      - grazingFlux_phy                        &
      - grazingFlux_phy2                       &     
  
                                             ) * dt + sms(k,iphyn)
   else
    sms(k,iphyn)      = (                      &
      + N_assim                      * PhyC    &
      - lossN * limitFacN            * PhyN    &
      - aggregationRate              * phyN    &
      - grazingFlux_phy                        &
                                              ) * dt + sms(k,iphyn)
   endif
!-------------------------------------------------------------------------------
! Phytoplankton C
   if (REcoM_Second_Zoo) then

    sms(k,iphyc)      = (                      &
      + Cphot                        * PhyC    &
      - lossC * limitFacN            * PhyC    &
      - phyRespRate                  * PhyC    &
      - aggregationRate              * PhyC    &
      - grazingFlux_phy * recipQuota           &
      - grazingFlux_phy2 * recipQuota          &
                                             ) * dt + sms(k,iphyc)
   else
    sms(k,iphyc)      = (                      &
      + Cphot                        * PhyC    &
      - lossC * limitFacN            * PhyC    &
      - phyRespRate                  * PhyC    &
      - aggregationRate              * PhyC    &
      - grazingFlux_phy * recipQuota           &
                                             ) * dt + sms(k,iphyc)
   endif
!-------------------------------------------------------------------------------
! Phytoplankton ChlA
   if (REcoM_Second_Zoo) then
   
    sms(k,ipchl)       = (                       &
     	+ chlSynth                     * phyC    &
     	- KOchl                        * Chl     &
     	- aggregationRate              * Chl     &
     	- grazingFlux_phy * Chl2N                &
        - grazingFlux_phy2 * Chl2N               & 
                                               ) * dt + sms(k,ipchl)
   else
    sms(k,ipchl)       = (                       &
        + chlSynth                     * phyC    &
        - KOchl                        * Chl     &
        - aggregationRate              * Chl     &
        - grazingFlux_phy * Chl2N                &
                                               ) * dt + sms(k,ipchl)

   endif
!-------------------------------------------------------------------------------
! Detritus N
   if (Grazing_detritus) then
    sms(k,idetn)       = (                       &
	+ grazingFlux_phy                        &
        - grazingFlux_phy * grazEff              &
        + grazingFlux_dia                        &
        - grazingFlux_dia * grazEff              &
        + grazingFlux_Cocco                      &          
        - grazingFlux_Cocco * grazEff            &       
        - grazingFlux_Det * grazEff              & !!!Sloppy feeding is  thought because of grazing flux multiplied with grazeff 
        - grazingFlux_Det2 * grazEff             &
        + aggregationRate              * PhyN    &
        + aggregationRate              * DiaN    &
        + aggregationRate              * CoccoN  &   
        + hetLossFlux                            &
        - reminN * arrFunc             * DetN    &
                                               ) * dt + sms(k,idetn)
   else
    sms(k,idetn)       = (                       &
        + grazingFlux_phy                        &
        + grazingFlux_dia                        &
        + grazingFlux_Cocco                      &        
        - grazingFlux * grazEff                  &
        + aggregationRate              * PhyN    &
        + aggregationRate              * DiaN    &
        + aggregationRate              * CoccoN  &      
        + hetLossFlux                            &
        - reminN * arrFunc             * DetN    &
                                               ) * dt + sms(k,idetn)
   end if   
!-------------------------------------------------------------------------------
! Detritus C
   if (Grazing_detritus) then
    sms(k,idetc)       = (                            &
        + grazingFlux_phy * recipQuota                 &
        - grazingFlux_phy * recipQuota * grazEff       &
        + grazingFlux_Dia * recipQuota_Dia             &
        - grazingFlux_Dia * recipQuota_Dia * grazEff   &
        + grazingFlux_Cocco * recipQuota_Cocco         &      
        - grazingFlux_Cocco * recipQuota_Cocco * grazEff &   
        - grazingFlux_Det * recipDet * grazEff         &
        - grazingFlux_Det2 * recipDet2 * grazEff       &     
        + aggregationRate              * phyC          &
        + aggregationRate              * DiaC          &
        + aggregationRate              * CoccoC        &    
        + hetLossFlux * recipQZoo                      &
        - reminC * arrFunc             * DetC          &
                                             )   * dt + sms(k,idetc)
   else
    sms(k,idetc)       = (                             &
        + grazingFlux_phy * recipQuota                 &
        - grazingFlux_phy * recipQuota * grazEff       &
        + grazingFlux_Dia * recipQuota_Dia             &
        - grazingFlux_Dia * recipQuota_Dia * grazEff   &
        + grazingFlux_Cocco * recipQuota_Cocco         &   
        - grazingFlux_Cocco * recipQuota_Cocco * grazEff &     
        + aggregationRate              * phyC          &
        + aggregationRate              * DiaC          &
        + aggregationRate              * CoccoC        &     
        + hetLossFlux * recipQZoo                      &
        - reminC * arrFunc             * DetC          &
                                             )   * dt + sms(k,idetc)
   end if
!-------------------------------------------------------------------------------
! Heterotrophic N
   if (REcoM_Second_Zoo) then
    sms(k,ihetn)       = (                       &
    	+ grazingFlux * grazEff                  &
        - grazingFlux_het2                       &
     	- hetLossFlux                            &
     	- lossN_z                      * HetN    &
                                               ) * dt + sms(k,ihetn)
   else
    sms(k,ihetn)       = (                       &
        + grazingFlux * grazEff                  &
        - hetLossFlux                            &
        - lossN_z                      * HetN    &
                                               ) * dt + sms(k,ihetn)
   endif
!-------------------------------------------------------------------------------
! Heterotrophic C
   if (REcoM_Second_Zoo) then
    if (Grazing_detritus) then
     sms(k,ihetc)      = (                            &
        + grazingFlux_phy * recipQuota * grazEff     &
        + grazingFlux_Dia * recipQuota_Dia * grazEff &
        + grazingFlux_Cocco*recipQuota_Cocco*grazEff &    
        + grazingFlux_Det * recipDet * grazEff       &
        + grazingFlux_DetZ2 * recipDet2 * grazEff    &
        - hetLossFlux * recipQZoo                    &
        - grazingFlux_het2 * recipQZoo               &
        - lossC_z                      * HetC        &
        - hetRespFlux                                &
                                                ) * dt + sms(k,ihetc)
    else
     sms(k,ihetc)      = (                            &
        + grazingFlux_phy * recipQuota * grazEff     &
        + grazingFlux_Dia * recipQuota_Dia * grazEff &
        + grazingFlux_Cocco*recipQuota_Cocco*grazEff &       
        - hetLossFlux * recipQZoo                    &
        - grazingFlux_het2 * recipQZoo               &
        - lossC_z                      * HetC        &
        - hetRespFlux                                &
                                                ) * dt + sms(k,ihetc)

    end  if
   else
    sms(k,ihetc)      = (                              &
        + grazingFlux_phy * recipQuota * grazEff     &
        + grazingFlux_Dia * recipQuota_Dia * grazEff &
        + grazingFlux_Cocco*recipQuota_Cocco*grazEff &    
        - hetLossFlux * recipQZoo                    &
        - lossC_z                      * HetC        &
        - hetRespFlux                                &
                                                ) * dt + sms(k,ihetc)
   endif
!-------------------------------------------------------------------------------

   if (REcoM_Second_Zoo) then
 ! Second Zooplankton N                                                                                              
      sms(k,izoo2n)       = (                        &
         + grazingFlux2 * grazEff2                  &
         - Zoo2LossFlux                             &
         - lossN_z2                      * Zoo2N    &
         - Zoo2fecalloss_n                            & 
                                               ) * dt + sms(k,izoo2n)
 !-------------------------------------------------------------------------------                       
  ! Second Zooplankton C                                                                                 
                        
     if (Grazing_detritus) then
      sms(k,izoo2c)      = (                             &
         + grazingFlux_phy2 * recipQuota * grazEff2     &
         + grazingFlux_Dia2 * recipQuota_Dia * grazEff2 &
         + grazingFlux_Cocco2*recipQuota_Cocco*grazEff2 &   
         + grazingFlux_het2 * recipQZoo * grazEff2      &
         + grazingFlux_Det2 * recipDet * grazEff2       &
         + grazingFlux_DetZ22 * recipDet2 * grazEff2    &
         - zoo2LossFlux * recipQZoo2                    &
         - lossC_z2                      * Zoo2C        &
         - Zoo2RespFlux                                 &
         - Zoo2fecalloss_c                              &
                                                ) * dt + sms(k,izoo2c)
     else
      sms(k,izoo2c)      = (                             &
         + grazingFlux_phy2 * recipQuota * grazEff2     &
         + grazingFlux_Dia2 * recipQuota_Dia * grazEff2 &
         + grazingFlux_het2 * recipQZoo * grazEff2      &
         + grazingFlux_Cocco2*recipQuota_Cocco*grazEff2 &  
         - zoo2LossFlux * recipQZoo2                    &
         - lossC_z2                      * Zoo2C        &
         - Zoo2RespFlux                                 &
         - Zoo2fecalloss_c                              &
                                                ) * dt + sms(k,izoo2c)
     end if   
 !---------------------------------------------------------------------------------
  ! Second Zooplankton Detritus N
     if (Grazing_detritus) then
      sms(k,idetz2n)       = (                       &
         + grazingFlux_phy2                       &
         - grazingFlux_phy2 * grazEff2            &         
         + grazingFlux_dia2                       &
         - grazingFlux_dia2 * grazEff2            & 
         + grazingFlux_Cocco2                     &     
         - grazingFlux_Cocco2 * grazEff2          &      
         + grazingFlux_het2                       &
         - grazingFlux_het2 * grazEff2            &
         - grazingFlux_DetZ2 * grazEff2           &
         - grazingFlux_DetZ22 * grazEff2          &   
         + Zoo2LossFlux                           &
         + Zoo2fecalloss_n                          &
         - reminN * arrFunc             * DetZ2N  &
                                               ) * dt + sms(k,idetz2n)
     else
      sms(k,idetz2n)       = (                       &
         + grazingFlux_phy2                       &
         + grazingFlux_dia2                       &
         + grazingFlux_Cocco2                     &      
         + grazingFlux_het2                       &
         - grazingFlux2 * grazEff2                &
         + Zoo2LossFlux                           &
         + Zoo2fecalloss_n                          &
         - reminN * arrFunc             * DetZ2N  &
                                               ) * dt + sms(k,idetz2n)
     end if
 !---------------------------------------------------------------------------------                                                                                            
                                                           
  ! Second Zooplankton Detritus C
     if (Grazing_detritus) then
      sms(k,idetz2c)       = (                             &
        + grazingFlux_phy2 * recipQuota                &
        - grazingFlux_phy2 * recipQuota * grazEff2     &
        + grazingFlux_Dia2 * recipQuota_Dia            &
        - grazingFlux_Dia2 * recipQuota_Dia * grazEff2 &
        + grazingFlux_Cocco2 * recipQuota_Cocco        &   
        - grazingFlux_Cocco2*recipQuota_Cocco*grazEff2 &  
        + grazingFlux_het2 * recipQZoo                 &
        - grazingFlux_het2 * recipQZoo * grazEff2      &
        - grazingFlux_DetZ2 * recipDet * grazEff2      &
        - grazingFlux_DetZ22 * recipDet2 * grazEff2    &
        + Zoo2LossFlux * recipQZoo2                    &
        + Zoo2fecalloss_c                   & 
        - reminC * arrFunc             * DetZ2C          &
                                             )   * dt + sms(k,idetz2c)
     else
      sms(k,idetz2c)       = (                             &
        + grazingFlux_phy2 * recipQuota                &
        - grazingFlux_phy2 * recipQuota * grazEff2     &
        + grazingFlux_Dia2 * recipQuota_Dia            &
        - grazingFlux_Dia2 * recipQuota_Dia * grazEff2 &
        + grazingFlux_Cocco2 * recipQuota_Cocco        &  
        - grazingFlux_Cocco2*recipQuota_Cocco*grazEff2 &  
        + grazingFlux_het2 * recipQZoo                 &
        - grazingFlux_het2 * recipQZoo * grazEff2      &
        + Zoo2LossFlux * recipQZoo2                    &
        + Zoo2fecalloss_c                    &
        - reminC * arrFunc             * DetZ2C          &
                                             )   * dt + sms(k,idetz2c)
     end if

 !-------------------------------------------------------------------------------------
  !Second Zooplankton  Detritus Si                                                                                                                                             
                                                           
     sms(k,idetz2si)     = (                           &
         + grazingFlux_dia2 * qSiN                    &
         - reminSiT                        * DetZ2Si &
                                             ) * dt + sms(k,idetz2si)

 !-------------------------------------------------------------------------------                                                                                              
 ! Detritus calcite   
     sms(k,idetz2calc)   = (               &
       + calc_loss_gra2                  &
       - calc_loss_gra2 * calc_diss_guts &
       - calc_diss2     * DetZ2Calc       &
                                           ) * dt + sms(k,idetz2calc)
    
   endif                                                                               
!-------------------------------------------------------------------------------
! DON (Extracellular organic N)
   if (REcoM_Second_Zoo) then
    sms(k,idon)      = (                        &
      + lossN * limitFacN              * phyN   &
      + lossN_d * limitFacN_Dia        * DiaN   &
      + lossN_c * limitFacN_Cocco      * CoccoN &       
      + reminN * arrFunc               * DetN   &
      + reminN * arrFunc               * DetZ2N &
      + lossN_z                        * HetN   &
      + lossN_z2                       * Zoo2N  &
      - rho_N * arrFunc                * DON    &
      + LocRiverDON                             &
                                             ) * dt + sms(k,idon)
   else
    sms(k,idon)      = (                        &
      + lossN * limitFacN              * phyN   &
      + lossN_d * limitFacN_Dia        * DiaN   &
      + lossN_c * limitFacN_Cocco      * CoccoN &      
      + reminN * arrFunc               * DetN   &
      + lossN_z                        * HetN   &
      - rho_N * arrFunc                * DON    &
      + LocRiverDON                             &
                                              ) * dt + sms(k,idon)
   endif
!-------------------------------------------------------------------------------
! EOC
   if (REcoM_Second_Zoo) then
    sms(k,idoc)       = (                       &
      + lossC * limitFacN              * phyC   &
      + lossC_d * limitFacN_dia        * DiaC   &
      + lossC_c * limitFacN_cocco      * CoccoC &       
      + reminC * arrFunc               * DetC   &
      + reminC * arrFunc               * DetZ2C &
      + lossC_z                        * HetC   &
      + lossC_z2                       * Zoo2C  &
      - rho_c1 * arrFunc               * EOC    &
      + LocRiverDOC                             &
                                              ) * dt + sms(k,idoc)	
   else   
    sms(k,idoc)       = (                       &
      + lossC * limitFacN              * phyC   &
      + lossC_d * limitFacN_dia        * DiaC   &
      + lossC_c * limitFacN_cocco      * CoccoC &      
      + reminC * arrFunc               * DetC   &
      + lossC_z                        * HetC   &
      - rho_c1 * arrFunc               * EOC    &
      + LocRiverDOC                             &
                                              ) * dt + sms(k,idoc)
   endif 	
!-------------------------------------------------------------------------------
! Diatom N
   if (REcoM_Second_Zoo) then
    sms(k,idian)      = (                      &
      + N_assim_dia                    * DiaC  &
      - lossN_d * limitFacN_dia        * DiaN  &
      - aggregationRate                * DiaN  &
      - grazingFlux_Dia                        &
      - grazingFlux_Dia2                       &
                                             ) * dt + sms(k,idian)
   else
    sms(k,idian)      = (                      &
      + N_assim_dia                    * DiaC  &
      - lossN_d * limitFacN_dia        * DiaN  &
      - aggregationRate                * DiaN  &
      - grazingFlux_Dia                        &
                                                 ) * dt + sms(k,idian)
   endif
!-------------------------------------------------------------------------------
! Diatom C
   if (REcoM_Second_Zoo) then
    sms(k,idiac)      = (                      &
      + Cphot_dia                      * DiaC  &
      - lossC_d * limitFacN_dia        * DiaC  &
      - phyRespRate_dia                * DiaC  &
      - aggregationRate                * DiaC  &
      - grazingFlux_dia * recipQuota_dia       &
      - grazingFlux_dia2 * recipQuota_dia      &
     	                                     ) * dt + sms(k,idiac)
   else
    sms(k,idiac)      = (                      &
      + Cphot_dia                      * DiaC  &
      - lossC_d * limitFacN_dia        * DiaC  &
      - phyRespRate_dia                * DiaC  &
      - aggregationRate                * DiaC  &
      - grazingFlux_dia * recipQuota_dia       &
                                             ) * dt + sms(k,idiac)
   endif
!-------------------------------------------------------------------------------
! Diatom Chl
   if (REcoM_Second_Zoo) then
    sms(k,idchl)      = (                       &
      + chlSynth_dia                   * DiaC   &
      - KOchl_dia                      * DiaChl &
      - aggregationRate                * DiaChl &
      - grazingFlux_dia * Chl2N_dia             &
      - grazingFlux_dia2 * Chl2N_dia            &                 
                                             ) * dt + sms(k,idchl)
   else
    sms(k,idchl)      = (                       &
      + chlSynth_dia                   * DiaC   &
      - KOchl_dia                      * DiaChl &
      - aggregationRate                * DiaChl &
      - grazingFlux_dia * Chl2N_dia             &
                                             ) * dt + sms(k,idchl)
   endif 
   
!-------------------------------------------------------------------------------
! Diatom Si
   if (REcoM_Second_Zoo) then
    sms(k,idiasi)        = (                    &
      + Si_assim                        * DiaC  &
      - lossN_d * limitFacN_dia         * DiaSi &
      - aggregationRate                 * DiaSi &
      - grazingFlux_dia * qSiN                  &
      - grazingFlux_dia2 * qSiN                 &
                                             ) * dt + sms(k,idiasi)
   else
    sms(k,idiasi)        = (                    &
      + Si_assim                        * DiaC  &
      - lossN_d * limitFacN_dia         * DiaSi &
      - aggregationRate                 * DiaSi &
      - grazingFlux_dia * qSiN                  &
                                                 ) * dt + sms(k,idiasi)  
   endif

!-------------------------------------------------------------------------------
! Coccolithophore N
   if (REcoM_Second_Zoo) then                                     
    sms(k,icocn)      = (                      &                
      + N_assim_cocco                 * CoccoC &               
      - lossN_c * limitFacN_cocco     * CoccoN &              
      - aggregationRate               * CoccoN &                
      - grazingFlux_Cocco                      &             
      - grazingFlux_Cocco2                     &            
                                             ) * dt + sms(k,icocn) 
   else                                                      
    sms(k,icocn)      = (                      &               
      + N_assim_cocco                 * CoccoC &             
      - lossN_c * limitFacN_cocco     * CoccoN &               
      - aggregationRate               * CoccoN &             
      - grazingFlux_Cocco                      &              
                                                 ) * dt + sms(k,icocn) 
   endif                                                      
!-------------------------------------------------------------------------------
! Coccolithophore C
    if (REcoM_Second_Zoo) then                              
    sms(k,icocc)      = (                      &                    
      + Cphot_cocco                   * CoccoC &                
      - lossC_c * limitFacN_cocco     * CoccoC &                   
      - phyRespRate_cocco             * CoccoC &                   
      - aggregationRate               * CoccoC &                   
      - grazingFlux_cocco * recipQuota_cocco   &                 
      - grazingFlux_Cocco2* recipQuota_cocco   &                  
                                             ) * dt + sms(k,icocc)  
   else                                                          
    sms(k,icocc)      = (                      &               
      + Cphot_cocco                   * CoccoC &             
      - lossC_c * limitFacN_cocco     * CoccoC &                
      - phyRespRate_cocco             * CoccoC &              
      - aggregationRate               * CoccoC &            
      - grazingFlux_cocco * recipQuota_cocco   &               
                                             ) * dt + sms(k,icocc)
   endif

!-------------------------------------------------------------------------------
! Coccolithophore Chl 
    if (REcoM_Second_Zoo) then                                  
    sms(k,icchl)      = (                         &          
      + ChlSynth_cocco                 * CoccoC   &           
      - KOchl_cocco                    * CoccoChl &         
      - aggregationRate                * CoccoChl &          
      - grazingFlux_cocco * Chl2N_cocco           &           
      - grazingFlux_Cocco2* Chl2N_cocco           &           
                                             ) * dt + sms(k,icchl)   
   else                                                         
    sms(k,icchl)      = (                         &                
      + chlSynth_cocco                 * CoccoC   &            
      - KOchl_cocco                    * CoccoChl &             
      - aggregationRate                * CoccoChl &           
      - grazingFlux_cocco * Chl2N_cocco           &          
                                             ) * dt + sms(k,icchl) 
   endif
!-------------------------------------------------------------------------------
! Detritus Si
   if (REcoM_Second_Zoo) then
    sms(k,idetsi)     = (                       &
      + aggregationRate                 * DiaSi &
      + lossN_d * limitFacN_dia         * DiaSi &
      + grazingFlux_dia * qSiN                  &
      + grazingFlux_dia2 * qSiN                 &
      - reminSiT                        * DetSi &
                                             ) * dt + sms(k,idetsi)
   else
    sms(k,idetsi)     = (                       &
      + aggregationRate                 * DiaSi &
      + lossN_d * limitFacN_dia         * DiaSi &
      + grazingFlux_dia * qSiN                  &
      - reminSiT                        * DetSi &
                                             ) * dt + sms(k,idetsi)
   endif
!-------------------------------------------------------------------------------
! Silicate
   if (REcoM_Second_Zoo) then
    sms(k,isi)        = (                         &
      - Si_assim                        * DiaC    &
      + reminSiT                        * DetSi   &
      + reminSiT                        * DetZ2Si &
      + LocRiverDSi                              &

                                             ) * dt + sms(k,isi)
   else 
    sms(k,isi)        = (                         &
      - Si_assim                        * DiaC    &
      + reminSiT                        * DetSi   &
      + LocRiverDSi                              &

                                             ) * dt + sms(k,isi)
   endif
!-------------------------------------------------------------------------------
! Fe
   if (REcoM_Second_Zoo) then
    if (use_Fe2N) then
         sms(k,ife) = ( Fe2N * (                  &
            - N_assim                 * PhyC      &
            - N_assim_dia             * DiaC      &
            - N_assim_cocco           * CoccoC    &         
            + lossN*limitFacN         * PhyN      &
            + lossN_d*limitFacN_dia   * DiaN      &
            + lossN_c*limitFacN_cocco * CoccoN    &          
            + reminN * arrFunc        * DetN      &
            + reminN * arrFunc        * DetZ2N    &
            + lossN_z                 * HetN      &
            + lossN_z2                * Zoo2N     &         
                                              )   &
            - kScavFe                 * DetC * FreeFe &
            - kScavFe                 * DetZ2C * FreeFe &
                          ) * dt           + sms(k,ife)
    else
          sms(k,ife)      = ( Fe2C *(          &
          -  Cphot                  * PhyC     &
          -  Cphot_dia              * DiaC     &
          -  Cphot_cocco            * CoccoC   &           
          +  phyRespRate            * PhyC     &
          +  phyRespRate_dia        * DiaC     &
          +  phyRespRate_cocco      * CoccoC   &        
          +  lossC*limitFacN        * phyC     &
          +  lossC_d*limitFacN_dia  * diaC     &
          +  lossC_c*limitFacN_cocco* CoccoC   &             
          +  reminC * arrFunc       * detC     &
          +  reminC * arrFunc       * DetZ2C   &
          +  lossC_z                * hetC     &
          +  hetRespFlux                       &
          +  lossC_z2               * Zoo2C    &
          +  zoo2RespFlux                      &
                                          )    &
          -  kScavFe                * DetC * FreeFe   & 
          -  kScavFe                * DetZ2C * FreeFe & 
                                            ) * dt + sms(k,ife)
    end if
   else
    if (use_Fe2N) then
         sms(k,ife) = ( Fe2N * (                  &
            - N_assim                 * PhyC      &
            - N_assim_dia             * DiaC      &
            - N_assim_cocco           * CoccoC    &        
            + lossN*limitFacN         * PhyN      &
            + lossN_d*limitFacN_dia   * DiaN      &
            + lossN_c*limitFacN_cocco * CoccoN    &          
            + reminN * arrFunc        * DetN      &
            + lossN_z                 * HetN      &
                                              )   &
            - kScavFe                 * DetC * FreeFe &
                          ) * dt           + sms(k,ife)
    else
          sms(k,ife)      = ( Fe2C *(          &
          -  Cphot                  * PhyC     &
          -  Cphot_dia              * DiaC     &
          -  Cphot_cocco            * CoccoC   &          
          +  phyRespRate            * PhyC     &
          +  phyRespRate_dia        * DiaC     &
          +  phyRespRate_cocco      * CoccoC   &         
          +  lossC*limitFacN        * phyC     &
          +  lossC_d*limitFacN_dia  * diaC     &
          +  lossC_c*limitFacN_cocco* CoccoC   &           
          +  reminC * arrFunc       * detC     &
          +  lossC_z                * hetC     &
          +  hetRespFlux                       &
                                          )    &
          -  kScavFe                * DetC * FreeFe   &
                                            ) * dt + sms(k,ife)
    end if
   endif

!-------------------------------------------------------------------------------
! Coccolithophore calcite
!#ifdef REcoM_calcification
   if (REcoM_Second_Zoo) then
    sms(k,iphycal)    = (                   &
      + calcification                       &
      - lossC_c * limitFacN_cocco * phyCalc &               
      - phyRespRate_cocco         * phyCalc &                 
      - calc_loss_agg                       &
      - calc_loss_gra                       &
      - calc_loss_gra2                      &

                                            ) * dt + sms(k,iphycal)
   else
    sms(k,iphycal)    = (                   &
      + calcification                       &
      - lossC_c * limitFacN_cocco * phyCalc &              
      - phyRespRate_cocco         * phyCalc &                 
      - calc_loss_agg                       &
      - calc_loss_gra                       &
                                            ) * dt + sms(k,iphycal)
   endif
!-------------------------------------------------------------------------------
! Detritus calcite
    sms(k,idetcal)   = (                    &
      + lossC_c * limitFacN_cocco * phyCalc &               
      + phyRespRate_cocco         * phyCalc &               
      + calc_loss_agg                       &
      + calc_loss_gra                       &
      - calc_loss_gra * calc_diss_guts      &
      - calc_diss     * DetCalc             &
                                           ) * dt + sms(k,idetcal)
!#endif
!-------------------------------------------------------------------------------
! Oxygen
   if (REcoM_Second_Zoo) then
    sms(k,ioxy)   = (               &
      + Cphot               * phyC  &
      - phyRespRate         * phyC  &
      + Cphot_dia           * diaC  &
      - phyRespRate_dia     * diaC  &
      + Cphot_cocco         * CoccoC&      
      - phyRespRate_cocco   * CoccoC&      
      - rho_C1  * arrFunc   * EOC   &
      - hetRespFlux                 &
      - Zoo2RespFlux                 &

                                        )*redO2C * dt + sms(k,ioxy)     
   
   else
    sms(k,ioxy)   = (               &
      + Cphot               * phyC  &
      - phyRespRate         * phyC  &
      + Cphot_dia           * diaC  &
      - phyRespRate_dia     * diaC  &
      + Cphot_cocco         * CoccoC&       
      - phyRespRate_cocco   * CoccoC&       
      - rho_C1  * arrFunc   * EOC   &
      - hetRespFlux                 &
                                      )*redO2C * dt + sms(k,ioxy)
   endif

!-------------------------------------------------------------------------------
! Diagnostics: Averaged rates
	
	recipbiostep    = 1.d0/real(biostep)

!*** Net primary production [mmol C /(m3 * day)]
	Diags3Dloc(k,1) = Diags3Dloc(k,1) + (   &
     	+ Cphot                   * PhyC  &
     	- PhyRespRate             * PhyC  &
     	) * recipbiostep

	Diags3Dloc(k,2) = Diags3Dloc(k,2) + (   &
     	+ Cphot_dia               * DiaC  &
     	- PhyRespRate_dia         * DiaC  &
     	) * recipbiostep

        Diags3Dloc(k,21) = Diags3Dloc(k,21) + (   &               
        + Cphot_cocco             * CoccoC  &                  
        - PhyRespRate_cocco       * CoccoC  &                   
        ) * recipbiostep


!*** Gross primary production [mmol C /(m3 * day)]
	Diags3Dloc(k,3) = Diags3Dloc(k,3) + (   &
     	+ Cphot                   * PhyC  &
     	) * recipbiostep

	Diags3Dloc(k,4) = Diags3Dloc(k,4) + (   &
     	+ Cphot_dia               * DiaC  &
     	) * recipbiostep

        Diags3Dloc(k,22) = Diags3Dloc(k,22) + (   &             
       + Cphot_cocco             * CoccoC  &
      ) * recipbiostep

!*** Net N-assimilation [mmol N/(m3 * day)]
	Diags3Dloc(k,5) = Diags3Dloc(k,5) + (   & 
     	+ N_assim                 * PhyC  &
     	- lossN * limitFacN       * PhyN  &
     	) * recipbiostep

	Diags3Dloc(k,6) = Diags3Dloc(k,6) + (   & 
     	+ N_assim_dia             * DiaC  &
     	- lossN * limitFacN_dia   * DiaN  &
     	) * recipbiostep

        Diags3Dloc(k,23) = Diags3Dloc(k,23) + (   &                 
        + N_assim_cocco           * CoccoC  &
        - lossN * limitFacN_cocco * CoccoN  &
        ) * recipbiostep

	Diags3Dloc(k,7) = Diags3Dloc(k,7) + (   &
     	+ N_assim                 * PhyC  &
     	) * recipbiostep

	Diags3Dloc(k,8) = Diags3Dloc(k,8) + (   & 
     	+ N_assim_dia             * DiaC  &
     	) * recipbiostep

        Diags3Dloc(k,24) = Diags3Dloc(k,24) + ( &                 
        + N_assim_cocco            * CoccoC &             
        ) * recipbiostep                                

!*** Total grazing of first zooplankton (with graz_eff, i.e. what reaches ZOO)
        Diags3Dloc(k,9) = Diags3Dloc(k,9) + (   &                             
        + grazingFlux_phy * recipQuota * grazEff     &
        + grazingFlux_Dia * recipQuota_Dia * grazEff &
        + grazingFlux_Cocco * recipQuota_Cocco * grazEff &     
        ) * recipbiostep                                 

!*** Grazing on small Phytoplankton by First Zooplankton (without grazeff, i.e. loss term for PHY)
        Diags3Dloc(k,10) = Diags3Dloc(k,10)+ (   &   
        + grazingFlux_phy * recipQuota           &
        ) * recipbiostep                                 

!*** Grazing on diatoms by First Zooplankton (without grazeff, i.e. loss term for DIA)
        Diags3Dloc(k,11) = Diags3Dloc(k,11) + (  &                                                 
        + grazingFlux_dia * recipQuota_dia       & 
        ) * recipbiostep  

!*** Grazing on cocclithophores by First Zooplankton (without grazeff, i.e. loss term for COCCO)
        Diags3Dloc(k,25) = Diags3Dloc(k,25) +(   &            
        + grazingFlux_Cocco * recipQuota_cocco   &          
        ) * recipbiostep                                      

!*** zooplankton1 respiration
        Diags3Dloc(k,12) = Diags3Dloc(k,12) + (   &
        + HetRespFlux                             &
        ) * recipbiostep

!*** calc_diss
        Diags3Dloc(k,13) = Diags3Dloc(k,13) + (   &
        + calc_diss * DetCalc                     &
        ) * recipbiostep

!***    aggregation by  small phytoplankton                                                                      
        Diags3Dloc(k,14) = Diags3Dloc(k,14) + (   &   
        + aggregationrate * PhyC                  &
        ) * recipbiostep

!***    aggregation by  diatoms                                                                                  
        Diags3Dloc(k,15) = Diags3Dloc(k,15) + (   &
        + aggregationrate * DiaC                  &
        ) * recipbiostep  

!***    aggregation by coccolithophores
        Diags3Dloc(k,26) = Diags3Dloc(k,26) + (   &          
        + aggregationrate * CoccoC                &            
        ) * recipbiostep                                      

!*** excrection of DOC by phytoplankton
        Diags3Dloc(k,16) = Diags3Dloc(k,16) + (   &
        + lossC * limitFacN              * phyC   &
        ) * recipbiostep  
  
!*** excrection of DOC by diatoms
        Diags3Dloc(k,17) = Diags3Dloc(k,17) + (   &
        + lossC_d * limitFacN_dia        * DiaC   &
        ) * recipbiostep  

!*** excretion of DOC by coccolithophores
        Diags3Dloc(k,27) = Diags3Dloc(k,27) + (   &       
        + lossC_c * limitFacN_cocco      * CoccoC &            
        ) * recipbiostep                                    

!*** calcification
        Diags3Dloc(k,18) = Diags3Dloc(k,18) + (   &
        + calcification                           &
        ) * recipbiostep  

! phy respiration
	Diags3Dloc(k,19) = Diags3Dloc(k,19) + (   &
     	+ PhyRespRate             * PhyC          &
     	) * recipbiostep

! dia respiration
	Diags3Dloc(k,20) = Diags3Dloc(k,20) + (   &
     	+ PhyRespRate_dia         * DiaC          &
     	) * recipbiostep

! cocco resipration
        Diags3Dloc(k,28) = Diags3Dloc(k,28) + (   &            
        + PhyRespRate_cocco       * CoccoC        &          
        ) * recipbiostep                                    


  end do ! Main vertikal loop ends
  
!-------------------------------------------------------------------------------
! Remineralization from the sediments into the bottom layer
!		kLoc = Nn      
!*** DIN ***		
    decayBenthos(1) = decayRateBenN * LocBenthos(1)
    LocBenthos(1)      = LocBenthos(1)   - decaybenthos(1) * dt ! / depth of benthos
    sms(Nn,idin)    = sms(Nn,idin) + decayBenthos(1) * dt  &
                      * recipdzF(Nn)

!*** DIC ***
    decayBenthos(2) = decayRateBenC * LocBenthos(2)
    LocBenthos(2)      = LocBenthos(2)   - decaybenthos(2) * dt ! / depth of benthos
    sms(Nn,idic)    = sms(Nn,idic) + decayBenthos(2) * dt  &
                      * recipdzF(Nn)

!*** Si ***
    decayBenthos(3) = decayRateBenSi * LocBenthos(3)
    LocBenthos(3)      = LocBenthos(3)   - decaybenthos(3) * dt ! / depth of benthos
    sms(Nn,isi)     = sms(Nn,isi)  + decayBenthos(3) * dt  &
                      * recipdzF(Nn)

!*** Calc: DIC, Alk ***
!#ifdef REcoM_calcification
    decayBenthos(4) = calc_diss_ben * LocBenthos(4)  
    LocBenthos(4)      = LocBenthos(4)   - decayBenthos(4) * dt ! / depth of benthos
    sms(Nn,idic)    = sms(Nn,idic) + decayBenthos(4) * dt  &
                      * recipdzF(Nn)
    sms(Nn,ialk)    = sms(Nn,ialk) + decayBenthos(4) * dt  &
                      * recipdzF(Nn) * 2.d0
!#endif
!*** DFe ***
  if(use_Fe2N) then 
    Ironflux          = decayRateBenN * LocBenthos(1) * Fe2N_benthos
  else
    Ironflux          = decayRateBenC * LocBenthos(2) * Fe2C_benthos
  end if
  sms(Nn,ife)       = sms(Nn,ife) + Ironflux * recipdzF(Nn) * dt
  

!*** O2 ***
   sms(Nn,ioxy)    = sms(Nn,ioxy) - decayBenthos(2) * redO2C *dt  &
                      * recipdzF(Nn)

  end do ! Main time loop ends

!-------------------------------------------------------------------------------
! Sinking

  dt_sink = dt * real(biostep)	
!  dt = dt_s / SecondsPerDay ! Physics time step converted to [day] as sinking velocity is in [m/day]

  if (VPhy .gt. 0.1) then
    call REcoM_sinking(dt,Nn,SinkVel(:,ivphy),dzF,recipDzF,state(:,iphyn),sink,zF)
    sms(:,iphyn) = sms(:,iphyn) + sink

    call REcoM_sinking(dt,Nn,SinkVel(:,ivphy),dzF,recipDzF,state(:,iphyc),sink,zF)
    sms(:,iphyc) = sms(:,iphyc) + sink	

    call REcoM_sinking(dt,Nn,SinkVel(:,ivphy),dzF,recipDzF,state(:,ipchl),sink,zF)
    sms(:,ipchl) = sms(:,ipchl) + sink

!#ifdef REcoM_calcification
  
  end if		

  if (VDia .gt. 0.1) then
    call REcoM_sinking(dt,Nn,SinkVel(:,ivdia),dzF,recipDzF,state(:,idian),sink,zF)
    sms(:,idian) = sms(:,idian) + sink

    call REcoM_sinking(dt,Nn,SinkVel(:,ivdia),dzF,recipDzF,state(:,idiac),sink,zF)
    sms(:,idiac) = sms(:,idiac) + sink	
    
    call REcoM_sinking(dt,Nn,SinkVel(:,ivdia),dzF,recipDzF,state(:,idchl),sink,zF)
    sms(:,idchl) = sms(:,idchl) + sink
  
    call REcoM_sinking(dt,Nn,SinkVel(:,ivdia),dzF,recipDzF,state(:,idiasi),sink,zF)
    sms(:,idiasi) = sms(:,idiasi) + sink
  end if		


  if (VCocco .gt. 0.1) then                                                          
    call REcoM_sinking(dt,Nn,SinkVel(:,ivcoc),dzF,recipDzF,state(:,icocn),sink,zF)   
    sms(:,icocn) = sms(:,icocn) + sink                                              

    call REcoM_sinking(dt,Nn,SinkVel(:,ivcoc),dzF,recipDzF,state(:,icocc),sink,zF)    
    sms(:,icocc) = sms(:,icocc) + sink                                               

    call REcoM_sinking(dt,Nn,SinkVel(:,ivcoc),dzF,recipDzF,state(:,icchl),sink,zF)   
    sms(:,icchl) = sms(:,icchl) + sink                                               

!#ifdef REcoM_calcification                                                                                                                                                                                

    call REcoM_sinking(dt,Nn,SinkVel(:,ivcoc),dzF,recipDzF,state(:,iphycal),sink,zF)  
    sms(:,iphycal) = sms(:,iphycal) + sink                                           

!#endif                                                      
  end if                                                                        


  if (VDet .gt. 0.1) then
    call REcoM_sinking(dt,Nn,SinkVel(:,ivdet),dzF,recipDzF,state(:,idetn),sink,zF)
    sms(:,idetn) = sms(:,idetn) + sink

    call REcoM_sinking(dt,Nn,SinkVel(:,ivdet),dzF,recipDzF,state(:,idetc),sink,zF)
    sms(:,idetc) = sms(:,idetc) + sink	

    call REcoM_sinking(dt,Nn,SinkVel(:,ivdet),dzF,recipDzF,state(:,idetsi),sink,zF)
    sms(:,idetsi) = sms(:,idetsi) + sink

!#ifdef REcoM_calcification

    call REcoM_sinking(dt,Nn,SinkVel(:,ivdet),dzF,recipDzF,state(:,idetcal),sink,zF)
    sms(:,idetcal) = sms(:,idetcal) + sink
!#endif
   if (REcoM_Second_Zoo) then
    call REcoM_sinking2(dt,Nn,SinkVel(:,ivdetsc),dzF,recipDzF,state(:,idetz2n),sink,zF)
    sms(:,idetz2n) = sms(:,idetz2n) + sink

    call REcoM_sinking2(dt,Nn,SinkVel(:,ivdetsc),dzF,recipDzF,state(:,idetz2c),sink,zF)
    sms(:,idetz2c) = sms(:,idetz2c) + sink

    call REcoM_sinking2(dt,Nn,SinkVel(:,ivdetsc),dzF,recipDzF,state(:,idetz2si),sink,zF)
    sms(:,idetz2si) = sms(:,idetz2si) + sink

    call REcoM_sinking2(dt,Nn,SinkVel(:,ivdetsc),dzF,recipDzF,state(:,idetz2calc),sink,zF)
    sms(:,idetz2calc) = sms(:,idetz2calc) + sink
   end if
  end if		

!-------------------------------------------------------------------------------
! Benthic layer: Detritus and phytoplankton sink into the benthic layer and are 
! lost from the water column
! (but remineralized and re-released in dissolved form later
 if (REcoM_Second_Zoo) then
  if (allow_var_sinking) then
    Vben_det = Vdet_a * abs(zF(Nn)) + VDet  
    Vben_det_seczoo = VDet_zoo2
    Vben_phy = Vdet_a * abs(zF(Nn)) + VPhy
    Vben_dia = Vdet_a * abs(zF(Nn)) + VDia
    Vben_coc = Vdet_a * abs(zF(Nn)) + VCocco                                         
  else
    Vben_det = VDet
    Vben_det_seczoo = VDet_zoo2
    Vben_phy = VPhy
    Vben_dia = VDia
    Vben_coc = VCocco                                                  
  endif
 else
  if (allow_var_sinking) then
    Vben_det = Vdet_a * abs(zF(Nn)) + VDet
    Vben_phy = Vdet_a * abs(zF(Nn)) + VPhy
    Vben_dia = Vdet_a * abs(zF(Nn)) + VDia
    Vben_coc = Vdet_a * abs(zF(Nn)) + VCocco                                
  else
    Vben_det = VDet
    Vben_phy = VPhy
    Vben_dia = VDia
    Vben_coc = VCocco                                                  
  endif
 endif

!*** Det *** (index: 1=N,2=C,3=Si,4=calc)
    wFluxDet(1)    = Vben_det * state(Nn,idetn)
    sms(Nn,idetn)  = sms(Nn,idetn)    - wFluxDet(1) * dt  &
                     * recipDzF(Nn)
    LocBenthos(1)  = LocBenthos(1)    + wFluxDet(1) * dt ! / thickness of benthos layer

    wFluxDet(2)    = Vben_det * state(Nn,idetc)
    Fc             = 0.1d0 * wFluxDet(2)                                     ! Conversion of detC from [mmolC/m2/day] => [umol/cm2/day]
    LocDenit = zero

    sms(Nn,idetc)  = sms(Nn,idetc)    - wFluxDet(2) * dt  &
                     * recipDzF(Nn)
    LocBenthos(2)  = LocBenthos(2)    + wFluxDet(2) * dt ! / thickness of benthos layer

    wFluxDet(3)    = Vben_det * state(Nn,idetsi)
    sms(Nn,idetsi) = sms(Nn,idetsi)   - wFluxDet(3) * dt  &
                     * recipDzF(Nn)
    LocBenthos(3)  = LocBenthos(3)    + wFluxDet(3) * dt ! / thickness of benthos layer
	  
!#ifdef REcoM_calcification 
    wFluxDet(4)     = Vben_det * state(Nn,idetcal)
    sms(Nn,idetcal) = sms(Nn,idetcal) - wFluxDet(4) * dt &
                      * recipDzF(Nn) 
    LocBenthos(4)   = LocBenthos(4)   + wFluxDet(4) * dt ! / thickness of benthos layer
!#endif

!---------------------------------

!*** Phy *** (index: 1=N,2=C)
    wFluxPhy(1)     = Vben_phy * state(Nn,iphyn)
    sms(Nn,iphyn)   = sms(Nn,iphyn)   - wFluxPhy(1) * dt  &
                      * recipDzF(Nn)
    LocBenthos(1)   = LocBenthos(1)   + wFluxPhy(1) * dt ! / thickness of benthos layer

    wFluxPhy(2)     = Vben_phy * state(Nn,iphyc)
    sms(Nn,iphyc)   = sms(Nn,iphyc)   - wFluxPhy(2) * dt  &
                      * recipDzF(Nn)
    LocBenthos(2)   = LocBenthos(2)   + wFluxPhy(2) * dt ! / thickness of benthos layer

    wFluxPhy(3)     = Vben_phy * state(Nn,ipchl)                            
    sms(Nn,ipchl)   = sms(Nn,ipchl)   - wFluxPhy(3) * dt  &              
     	              * recipDzF(Nn)
!------------------------------------
if (REcoM_Second_Zoo) then
    wFluxDet(5)    = Vben_det_seczoo * state(Nn,idetz2n)
    sms(Nn,idetz2n)  = sms(Nn,idetz2n)    - wFluxDet(5) * dt  &
                     * recipDzF(Nn)
    LocBenthos(1)  = LocBenthos(1)    + wFluxDet(5) * dt ! / thickness of benthos layer  

    wFluxDet(6)    = Vben_det_seczoo * state(Nn,idetz2c)
    sms(Nn,idetz2c)  = sms(Nn,idetz2c)    - wFluxDet(6) * dt  &
                     * recipDzF(Nn)
    LocBenthos(2)  = LocBenthos(2)    + wFluxDet(6) * dt ! / thickness of benthos layer  

    wFluxDet(7)    = Vben_det_seczoo * state(Nn,idetz2si)
    sms(Nn,idetz2si)  = sms(Nn,idetz2si)    - wFluxDet(7) * dt  &
                        * recipDzF(Nn)
    LocBenthos(3)  = LocBenthos(3)    + wFluxDet(7) * dt ! / thickness of benthos layer  

    wFluxDet(8)    = Vben_det_seczoo * state(Nn,idetz2calc)
    sms(Nn,idetz2calc)  = sms(Nn,idetz2calc)    - wFluxDet(8) * dt  &
                        * recipDzF(Nn)
    LocBenthos(4)  = LocBenthos(4)    + wFluxDet(8) * dt ! / thickness of benthos layer  
endif
!------------------------------------

!*** Dia *** (index: 1=N,2=C,3=Si)
    wFluxDia(1)     = Vben_dia * state(Nn,idian)
    sms(Nn,idian)   = sms(Nn,idian)   - wFluxDia(1) * dt  &
                      * recipDzF(Nn)
    LocBenthos(1)   = LocBenthos(1)   + wFluxDia(1) * dt ! / thickness of benthos layer
 
    wFluxDia(2)     = Vben_dia * state(Nn,idiac)
    sms(Nn,idiac)   = sms(Nn,idiac)   - wFluxDia(2) * dt  &
                      * recipDzF(Nn)
    LocBenthos(2)   = LocBenthos(2)   + wFluxDia(2) * dt ! / thickness of benthos layer

    wFluxDia(3)     = Vben_dia * state(Nn,idiasi)
    sms(Nn,idiasi)  = sms(Nn,idiasi)  - wFluxDia(3) * dt  &
                      * recipDzF(Nn)
    LocBenthos(3)   = LocBenthos(3)   + wFluxDia(3) * dt ! / thickness of benthos layer
    
    wFluxDia(4)     = Vben_dia * state(Nn,idchl)
    sms(Nn,idchl)   = sms(Nn,idchl)   - wFluxDia(4) * dt  &
     	              * recipDzF(Nn)
!-------------------------------------

!*** Cocco *** (index: 1=N,2=C)
                                                                                                                                                                                                                                              
    wFluxCocco(1)   = Vben_coc * state(Nn,icocn)                                           
    sms(Nn,icocn)   = sms(Nn,icocn)   - wFluxCocco(1) * dt  &                                  
                      * recipDzF(Nn)                                                           
    LocBenthos(1)   = LocBenthos(1)   + wFluxCocco(1) * dt ! / thickness of benthos layer      

    wFluxCocco(2)   = Vben_coc * state(Nn,icocc)                                                 
    sms(Nn,icocc)   = sms(Nn,icocc)   - wFluxCocco(2) * dt  &                                
                      * recipDzF(Nn)                                                          
    LocBenthos(2)   = LocBenthos(2)   + wFluxCocco(2) * dt ! / thickness of benthos layer      
!#ifdef REcoM_calcification                                                                    
    wFluxCocco(3)     = Vben_coc * state(Nn,iphycal)                                            
    sms(Nn,iphycal) = sms(Nn,iphycal) - wFluxCocco(3) * dt  &                                     
                      * recipDzF(Nn)                                                             
    LocBenthos(3)   = LocBenthos(3)   + wFluxCocco(3) * dt ! / thickness of benthos layer     
!#endif                                                                                        

    wFluxCocco(4)   = Vben_coc * state(Nn,icchl)                                          
    sms(Nn,icchl)   = sms(Nn,icchl)   - wFluxCocco(4) * dt  &                         
                      * recipDzF(Nn)                                                  

end subroutine REcoM_sms

!-------------------------------------------------------------------------------
! Function for calculating limiter
!-------------------------------------------------------------------------------

function recom_limiter(slope,qa,qb)
  use recom_config
  Implicit None
  Real(kind=8) :: recom_limiter
  Real(kind=8) :: slope, qa, qb
  Real(kind=8) :: dq
	
  dq = qa - qb
  if (REcoM_Geider_limiter) then
    recom_limiter = max(min( -slope*dq, 1.d0),0.d0)
  else
    recom_limiter = 1.d0 - exp( -slope*( abs(dq)-dq )**2)
  endif
  return
  end

!-------------------------------------------------------------------------------
! Function for iron chemistry
!-------------------------------------------------------------------------------

function iron_chemistry(Fe, totalLigand, ligandStabConst)
  implicit none

  Real(kind=8) :: iron_chemistry
  Real(kind=8) :: Fe, totalLigand, ligandStabConst ! Input
  Real(kind=8) :: FreeFe                          ! Output
  Real(kind=8) :: ligand,FeL,a,b,c,discrim

! Abbrevations
  a = ligandstabConst
  b = ligandstabConst * (Fe - totalLigand) + 1.d0
  c = -totalLigand
  discrim = b*b - 4.d0 * a * c
	
  if (a .ne. 0.d0 .and. discrim .ge. 0.d0) then
    ligand = ( -b + sqrt(discrim) ) / (2.d0 * a)
    FeL    = totalLigand - ligand
    freeFe = Fe - FeL
  else ! No free iron
    freeFe = 0.d0
  end if

  iron_chemistry = freeFe

  return
  end

!===============================================================================
! Subroutine calculating sinking of Detritus, diatoms and possibly phy and coccos
!===============================================================================
subroutine REcoM_sinking(dt,Nn,wF,dzF,recipDzF,C,sink,zF)
  use recom_config
!  use REcoM_constants
  implicit none
  
  Real(kind=8)                 :: dt
  Integer                      :: Nn
  Integer		               :: k,km1,km2,kp1
  Real(kind=8),dimension(Nn)   :: C
  Real(kind=8),dimension(Nn+1) :: wF            ! Vertical velocity for fluxes
  Real(kind=8)                 :: wLoc,wM,wP
  Real(kind=8)                 :: Rjp,Rj,Rjm
  Real(kind=8),dimension(Nn)   :: dzF,recipdzF
  Real(kind=8)                 :: cfl, d0, d1, thetaP, thetaM, psiP, psiM
  Real(kind=8)                 :: onesixth	= 	1.d0/6.d0
  Real(kind=8)                 :: wflux
  Real(kind=8)                 :: wfluxkp1
  Real(kind=8),dimension(Nn)   :: sink	
  Real(kind=8),dimension(Nn+1) :: zF                   ! [m] Depth of fluxes
  Real(8),dimension(Nn+1)	   :: new_wF
  
  if (allow_var_sinking) then
    do k=one,Nn+one
      new_wF(k)= Vdet_a * abs(zF(k)) + wF(k) ! Applies to phy, dia and det. and coccos
    enddo
  else
    new_wF = wF
  endif
  
  wfluxkp1   = 0.d0
  wflux      = 0.d0	
	
  do k=Nn,2,-1
    km1      = max(1,k-1)
    km2      = max(1,k-2)
    kp1      = min(Nn,k+1)
    wLoc     = -new_wF(k)
    wP       = wLoc + abs(wLoc)
    wM       = wLoc - abs(wLoc)
    Rjp	   = C(k)	-	C(kp1)
    Rj       = C(km1)-	C(k)
    Rjm	   = C(km2)-	C(km1)	
    cfl      = abs(wLoc * dt * recipDzF(k))
    d0       = (2.d0 - cfl)*(1.d0 - cfl)*onesixth
    d1       = (1.d0 - cfl*cfl)*onesixth	
    thetaP   = Rjm/(1.d-20+Rj)
    psiP     = d0 + d1*thetaP
    psiP     = max(0.d0, min(min(1.d0,psiP), &
                 (1.d0-cfl)/(1.d-20+cfl)*thetaP))
    thetaM   = Rjp/(1.d-20 + Rj)	
    psiM     = d0 + d1*thetaM
    psiM     = max(0.d0, min(min(1.d0,psiM), &
		(1.d0-cfl)/(1.d-20-cfl)*thetaM))
    wflux    = ( 0.5*wP*(C(k)  + psiM * Rj)+ &
		0.5*wM*(C(km1)+ psiP * Rj))
    sink(k)  =	-(wflux-wfluxkp1)*recipDzF(k)*dt
    wfluxkp1 =	wflux
  enddo
  k          = 1
  wflux      = 0
  sink(k)    = -(wflux-wfluxkp1)*recipDzF(k)*dt
	
end subroutine REcoM_sinking

!===============================================================================                                                                                                                            
! Subroutine for Second Detritus sinking                                                                                                                                                                
!=============================================================================== 
subroutine REcoM_sinking2(dt,Nn,wF,dzF,recipDzF,C,sink,zF)
  use recom_config
!  use REcoM_constants
  implicit none

  Real(kind=8)                 :: dt
  Integer                      :: Nn
  Integer                              :: k,km1,km2,kp1
  Real(kind=8),dimension(Nn)   :: C
  Real(kind=8),dimension(Nn+1) :: wF            ! Vertical velocity for fluxes                                                                                                                              
  Real(kind=8)                 :: wLoc,wM,wP
  Real(kind=8)                 :: Rjp,Rj,Rjm
  Real(kind=8),dimension(Nn)   :: dzF,recipdzF
  Real(kind=8)                 :: cfl, d0, d1, thetaP, thetaM, psiP, psiM
  Real(kind=8)                 :: onesixth      =       1.d0/6.d0
  Real(kind=8)                 :: wflux
  Real(kind=8)                 :: wfluxkp1
  Real(kind=8),dimension(Nn)   :: sink  
  Real(kind=8),dimension(Nn+1) :: zF                   ! [m] Depth of fluxes                                                                                                                             !   
  Real(8),dimension(Nn+1)          :: new_wF
!
  new_wF = wF

  wfluxkp1   = 0.d0
  wflux      = 0.d0

  do k=Nn,2,-1
    km1      = max(1,k-1)
    km2      = max(1,k-2)
    kp1      = min(Nn,k+1)
    wLoc     = -new_wF(k)
    wP       = wLoc + abs(wLoc)
    wM       = wLoc - abs(wLoc)
    Rjp    = C(k)       -       C(kp1)
    Rj       = C(km1)-  C(k)
    Rjm    = C(km2)-    C(km1)
    cfl      = abs(wLoc * dt * recipDzF(k))
    d0       = (2.d0 - cfl)*(1.d0 - cfl)*onesixth
    d1       = (1.d0 - cfl*cfl)*onesixth
    thetaP   = Rjm/(1.d-20+Rj)
    psiP     = d0 + d1*thetaP
    psiP     = max(0.d0, min(min(1.d0,psiP), &
                 (1.d0-cfl)/(1.d-20+cfl)*thetaP))
    thetaM   = Rjp/(1.d-20 + Rj)
    psiM     = d0 + d1*thetaM
    psiM     = max(0.d0, min(min(1.d0,psiM), &
                (1.d0-cfl)/(1.d-20-cfl)*thetaM))
    wflux    = ( 0.5*wP*(C(k)  + psiM * Rj)+ &
                0.5*wM*(C(km1)+ psiP * Rj))
    sink(k)  =  -(wflux-wfluxkp1)*recipDzF(k)*dt
    wfluxkp1 =  wflux
  enddo
  k          = 1

  wflux      = 0
  sink(k)    = -(wflux-wfluxkp1)*recipDzF(k)*dt

end subroutine REcoM_sinking2


!===============================================================================
! Subroutine to restore surface alkalinity
!===============================================================================
subroutine alk_restore(tracer)
  use recom_ocean_settings
  use REcoM_GloVar
!  use REcoM_constants
  use recom_config
  implicit none

  integer                   :: row, m, elnodes(3), elnodes2(3), col, q 
  real(kind=8)              :: aux, auxf, entries(3)
  real(kind=8), allocatable :: relax_alk(:)
  real(kind=8), dimension(recom_todim_nod3d,recom_num_tracer), INTENT(IN) :: tracer

  if (.not. restore_alkalinity) return

  allocate(relax_alk(recom_ToDim_nod2d))
  relax_alk       = 0.d0
  recom_sfc_force = 0.d0

  do row=1,recom_ToDim_nod2d
    m              = recom_nod3D_below_nod2D(1,row) 
    aux            = Alk_surf(row) - tracer(m,index_recom_tracer(ialk)) ! tracer 5 is alkalinity if only recom adds extra tracers
    relax_alk(row) = restore_alk_surf * aux
  end do


  do row=1,recom_myDim_elem2d
    elnodes2 = recom_elem2D_nodes(:,row)
    elnodes  = recom_nod3D_below_nod2D(1,elnodes2)
    auxf     = recom_voltriangle(row)/12.d0

    entries(1:3)  = auxf * relax_alk(elnodes2)

    do q=1,3
        col=elnodes2(q)
        recom_sfc_force(col,1)=recom_sfc_force(col,1)+sum(entries(1:3))+entries(q)
    end do

  end do

  deallocate(relax_alk)

end subroutine alk_restore

!===============================================================================
! Subroutine to distribute surface fluxes of din, fe, dic
!===============================================================================
subroutine surface_fluxes
  use recom_ocean_settings
  use REcoM_GloVar
!  use REcoM_constants
  use recom_config
  implicit none

  integer                   :: row, m, elnodes(3), elnodes2(3), col, q 
  real(kind=8)              :: aux, entries(12)


  do row=1,recom_myDim_elem2d
    elnodes2 = recom_elem2D_nodes(:,row)
    elnodes  = recom_nod3D_below_nod2D(1,elnodes2)
    aux     = recom_voltriangle(row)/12.d0

    entries(1:3)    = aux * AtmNInput(elnodes2) 
    entries(4:6)    = aux * AtmFeInput(elnodes2)
    entries(7:9)    = aux * GloCO2flux_seaicemask(elnodes2)
    entries(10:12)  = aux * GloO2flux_seaicemask(elnodes2)
!if(recom_mype==0) then
!print*,' GloCO2flux_seaicemask ', GloCO2flux_seaicemask(elnodes2)
!print*,' GloO2flux_seaicemask ', GloO2flux_seaicemask(elnodes2)
!endif

    do q=1,3
        col=elnodes2(q)
        recom_sfc_force(col,2)=recom_sfc_force(col,2)+sum(entries(1:3))  +entries(q)     ! din
        recom_sfc_force(col,3)=recom_sfc_force(col,3)+sum(entries(4:6))  +entries(q+3)   ! fe
        recom_sfc_force(col,4)=recom_sfc_force(col,4)+sum(entries(7:9))  +entries(q+6)   ! dic
        recom_sfc_force(col,5)=recom_sfc_force(col,5)+sum(entries(10:12))+entries(q+9)   ! oxy
    end do

  end do
 
end subroutine surface_fluxes
