!===============================================================================
! Subroutine for calculating flux-depth and thickness of control volumes
!===============================================================================
subroutine Depth_calculations(zNodes,Nn,zF,dzF,recipdzF,wF,thick,recipthick)
  use recom_config
!  use REcoM_constants
  use recom_ocean_settings
  implicit none

! Input
  Integer,                     intent(in)  :: Nn	     ! Total number of nodes
  Real(kind=8),dimension(Nn),  intent(in)  :: zNodes	     ! Depth of nodes
! Output
  real(kind=8),dimension(Nn+1),intent(out) :: zF             ! [m] Depth of fluxes
  real(kind=8),dimension(Nn),  intent(out) :: dzF            ! [m] Thickness of control volumes  
  real(kind=8),dimension(Nn),  intent(out) :: recipdzF       ! [1/m] 
  real(kind=8),dimension(Nn-1),intent(out) :: thick          ! [m] Distance between two nodes = thickness
  real(kind=8),dimension(Nn+1),intent(out) :: recipthick     ! [1/m] reciprocal of thickness
  real(kind=8),dimension(Nn+1,5),intent(out) :: wF           ! [m/day] Velocities of fluxes at the border of the control volumes
  real(kind=8),dimension(Nn,5)               :: wNodes       ! [m/day] Vertical velocities at each node
  Integer                                  :: k              ! Index for depth                     
!-------------------------------------------------------------------------------
! Sinking velocities
  wNodes(:,ivphy) = VPhy  
  wNodes(:,ivdia) = VDia
  wNodes(:,ivcoc) = VCocco          
  wNodes(:,ivdet) = VDet  
  wNodes(:,ivdetsc) = VDet_zoo2
  
  wNodes(1,:)     = 0.d0
  wF(1,:)         = 0.d0
  wF(Nn+1,:)      = 0.d0
!----------------------------------------------------
! calculate thickness of vertical layers
  do k = one,Nn-1
    thick(k)	=	zNodes(k) - zNodes(k+1)    ! depth is defined negative downwards, hence thickness is always positive
    recipthick(k)	=	1.d0/thick(k)
  end do
  
!-------------------------------------------------------------------------------
! Depth of the borders of the control volumes are calculated

  zF(1)		  = zNodes(1)
  zF(Nn+1)	  = zNodes(Nn)
if (REcoM_Second_Zoo) then
  do k=one,Nn-one
    zF(k+1)	  = zNodes(k)+0.5*(zNodes(k+1)-zNodes(k))
    wF(k+1,ivphy) = 0.5 * (wNodes(k,ivphy) + wNodes(k+1,ivphy)) 
    wF(k+1,ivdia) = 0.5 * (wNodes(k,ivdia) + wNodes(k+1,ivdia))
    wF(k+1,ivcoc) = 0.5 * (wNodes(k,ivcoc) + wNodes(k+1,ivcoc))         
    wF(k+1,ivdet) = 0.5 * (wNodes(k,ivdet) + wNodes(k+1,ivdet))
    wF(k+1,ivdetsc) = 0.5 * (wNodes(k,ivdetsc) + wNodes(k+1,ivdetsc))
  enddo
else
   do k=one,Nn-one
     zF(k+1)       = zNodes(k)+0.5*(zNodes(k+1)-zNodes(k))
     wF(k+1,ivphy) = 0.5 * (wNodes(k,ivphy) + wNodes(k+1,ivphy))
     wF(k+1,ivdia) = 0.5 * (wNodes(k,ivdia) + wNodes(k+1,ivdia))
     wF(k+1,ivcoc) = 0.5 * (wNodes(k,ivcoc) + wNodes(k+1,ivcoc))   
     wF(k+1,ivdet) = 0.5 * (wNodes(k,ivdet) + wNodes(k+1,ivdet))
  enddo
endif
!-------------------------------------------------------------------------------
! Thickness of control volumes are calculated (distance between fluxes)

  do k = one,Nn
    dzF(k)	=	zF(k)-zF(k+1)
    recipDzF(k)	=	1.d0/dzF(k)
  end do
!-----------------------------------------
!! check thicknesses:
!  if (recom_mype==1) write(*,*),'thick: ',thick
!  if (recom_mype==1) write(*,*),'dzF: ',dzF
  


end subroutine Depth_calculations

!===============================================================================
! Subroutine for calculating cos(AngleOfIncidence)
!===============================================================================
subroutine Cobeta(cosAngleOfIncidence,Latr, daynew, ndpyr)
  use recom_ocean_settings
  Implicit none

  integer, INTENT(IN) :: daynew, ndpyr
! Values from FESOM
  Real(kind=8),intent(in)    :: Latr
	
! Declarations
  Real(kind=8)               :: yearfrac              ! The fraction of the year that has passed [0 1]
  Real(kind=8)               :: yDay                  ! yearfrac in radians [0 2*pi]
  Real(kind=8)               :: declination   = 0.d0  ! Declination of the sun at present lat and time
  Real(kind=8)               :: CosAngleNoon  = 0.d0  ! Cos(Angle of Incidence) at Noon ?
  Real(kind=8),intent(out)   :: cosAngleOfIncidence   ! Cos(Angle of incidence)
! Constants
  Real(kind=8)		     :: nWater        = 1.33

! find day (****NOTE for year starting in winter*****)  
! Paltridge, G. W. and C. M. R. Platt, Radiative Processes in Meteorology and Climatology, Developments in Atmospheric Sciences, vol. 5, Elsevier Scientific Publishing Company, Amsterdam, Oxford, New York, 1976, ISBN 0-444-41444-4.
  yearfrac    = mod(real(daynew),real(ndpyr))/real(ndpyr)
  yDay        = 2 * fesom_pi * yearfrac
  declination = 0.006918           &
     	- 0.399912 * cos(yDay)     &
     	+ 0.070257 * sin(yDay)     &
     	- 0.006758 * cos(2 * yDay) &
     	+ 0.000907 * sin(2 * yDay) &
     	- 0.002697 * cos(3 * yDay) &
     	+ 0.001480 * sin(3 * yDay) 
  cosAngleNoon   = sin(LatR) * sin(declination) &
     				 + cos(LatR) * cos(declination)
  cosAngleOfIncidence &
                  = sqrt(1-(1-cosAngleNoon**2)/nWater**2)

end subroutine Cobeta

!===============================================================================
! Subroutine calculating the CO2-flux from the atmosphere to the ocean
!===============================================================================
!subroutine CO2Flux(DIC,TAlk,Temp,Sali)
!  use REcoM_Constants
!  use REcoM_LocVar
!  Implicit none
	
!  Real(kind=8),intent(in) :: DIC,TAlk
!  Real(kind=8),intent(in) :: Temp, Sali
!!  Real(kind=8),intent(out):: dflux
!! Local variables
!  Real(kind=8)            :: piston_vel, t1
!  Real(kind=8)            :: schmidt_num
!  Real(kind=8)            :: hplus2,co2star,co2starair,dco2star
!  Real(kind=8)            :: pCO2
!  Real(kind=8)            :: scl

  
!!-------------------------------------------------------------------------------
!! Piston velocity
!! Adapted from the OCMIP program, updated and extended by Christoph Voelker
!! This is the piston vel at Schmidt num = 660 cm/h
!  t1 = 0.5246d0 + 1.6256e-2 * Temp + 4.9946e-4 * Temp * Temp
!  piston_vel  = 2.5d0 * t1 + 0.3d0 * Uloc * Uloc
	
!  Schmidt_num = 2073.1d0 - 125.62d0 * Temp + 3.6276d0 * Temp * Temp &
!              - 4.3219e-2 * Temp * Temp * Temp
!! Piston velocity in [m/day] and at current Schmidt number 
!  piston_vel  = piston_vel * sqrt(660.d0/Schmidt_num) * 0.24d0		

!  call REcoM_cal_constants(Temp,Sali)

!! Calculating concentrations for borate
!! Total borate: Uppstrom (1974)	
!  scl = Sali/1.80655              ! Chlorinity
!  bt  = 0.000232d0 * scl/10.811d0 ! Total borate (unit: [mol/kg]) 
  
!! DIC is converted from [mmol/m3] to [mol/kg] (Molality)
!! Atm. CO2 is converted from [uatm] to [atm]
!  dic_molal   = DIC   * permil
!  talk_molal  = TAlk  * permil
!  pCO2        = LocAtmCO2 * permeg
  
!!  ! Iterating to find the pH (or rather the [H+] concentration) for given DIC and TALK
!!  if use_mocsy_solver
!!   H = solve_at_general(talk_molal, dic_molal, bt,                         & 
!!                       pt,     sit,                                        &
!!                       St, Ft,                                             &
!!                       K0, K1, K2, Kb, Kw, Ks, Kf, K1p, K2p, K3p, Ksi )
!!  else
!   call REcoM_iterate_ph
!!  end if 	
!  hplus = 10**(-8.1d0)
	
!! Calculate [CO2*] as defined in DOE Methods Handbook 1994, Ver. 2, ORNL/CDIAC-74, 
!! Dickson and Goyet, eds. (Ch. 2 p 10, Eg A.49)
!  hplus2     = hplus * hplus
!  co2star    = hplus2 + k1 * hplus + k1 * k2
!  co2star    = dic_molal * hplus2 /co2star
!  co2starair = pCO2 * ff
!  dco2star   = co2starair - co2star
	
!  pco2surf   = co2star/ff
	
!  dco2star   = dco2star/permil ! Unit conversion [mol/kg] => [mmol/m3]
!  pco2surf   = pco2surf/permeg ! Unit conversion [atm] => [uatm]
	
!  dflux      = piston_vel * dco2star

!end subroutine CO2Flux

!===============================================================================
! Subroutine that calculates constants for CO2-flux subroutine
!===============================================================================
!subroutine REcoM_cal_constants(Temp,Sali)
!  use REcoM_LocVar
!  Implicit none

!! Input
!  Real(kind=8) :: Temp,Sali
!! Auxiliary variables
!  Real(kind=8) :: tk,tk100,tk1002,invtk,dlogtk,s2,sqrts,s15
!  Real(kind=8) :: kb1,kb2,kw1 

!! Calculation of auxiliary variables
!  tk     = 273.15d0 + Temp
!  tk100  = tk/100.d0
!  tk1002 = tk100 * tk100
!  invtk  = 1.d0/tk
!  dlogtk = log(tk)
!  s2     = Sali * Sali
!  sqrts  = sqrt(Sali)
!  s15    = exp(1.5d0 * log (Sali)) ! Sali^1.5

!! ff = k0 * (1-pH2O)*correction for non-ideality
!! Weiss & Price (1980, Mar. Chem., 8, 347-359; Eq 13 with table 6 values)
!  ff = exp(-162.8301d0 + 218.2968d0 / tk100         &
!     + 90.9241d0 * log(tk100) - 1.47696d0 * tk1002  &
!     + Sali * (.025695d0 - .025225 * tk100          &
!     + 0.0049867 * tk1002))
     
!! k1 = [H][HCO3]/[H2CO3]
!! k2 = [H][CO3]/[HCO3]
!! Millero p. 644 (1995) using Mehrbach et al. data on seawater
!  k1 = 3670.7d0 * invtk - 62.008d0 + 9.7944 * dlogtk &
!     - 0.0118 * Sali + 0.000116 * s2
!  k1 = exp(-k1 * log(10.d0))

!  k2 = 1394.7d0 * invtk + 4.777 - 0.0184d0 * Sali + 0.0000118 * s2
!  k2 = exp(-k2 * log(10.d0))
	
!! kb = [H][BO2]/[HBO2]
!! Millero p. 669 (1995) using data from Dickson 1990
!  kb1 = -8966.9d0 - 2890.53d0 * sqrts - 77.942d0 * sali &
!      + 1.728d0 * s15 - 0.0996d0 * s2 
!  kb2 = -24.4344d0 - 25.085d0 * sqrts - 0.2474d0 * sali 
!  kb  = kb1 * invtk + 148.0248d0 + 137.1942d0 * sqrts & 
!      + 1.62142 * sali + kb2 * dlogtk + 0.053105d0 * sqrts * tk
!  kb  = exp(kb) 
	
!! kw = [H][OH]
!! Millero p. 670 (1995) using composite data
!  kw1 = 118.67d0 * invtk - 5.977d0 + 1.0495 * dlogtk
!  kw  = -13847.26d0 * invtk + 148.9652d0 - 23.6521d0 * dlogtk &
!      + kw1 * sqrts - 0.01615 * Sali
!  kw  = exp(kw)
		
!end subroutine

!!===============================================================================
!! Iteration subroutine for calculating ph for the CO2-flux
!!===============================================================================
!subroutine REcoM_iterate_ph
!  use REcoM_LocVar
!  use REcoM_constants
!  Implicit none
	
!  Real(kind=8) :: xguess, xsol ! Values for ph before and after iteration
!  Real(kind=8) :: FL,DF,FH,xl,xh,swap,dxold,dx,F,tmp
!  Integer      :: Maxit = 100, Niter  ! Maximum number of iterations and actual number of iterations
	
!  xguess = Hplus
	
!  call REcoM_talk_difference(x1,FL,DF)
!  call REcoM_talk_difference(x2,FH,DF)	
!  if (fl .lt. zero) then
!    xl   = x1
!    xh   = x2
!  else
!    xh   = x1
!    xl   = x2
!    swap = fl
!    fl   = fh
!    fh   = swap
!  end if   
!  xsol  = xguess
!  dxold = abs(x2-x1) 
!  dx    = dxold
	
!  call REcoM_talk_difference(xsol,f,df)
	
!  DO NITER=1,Maxit
!    IF(((xsol-XH)*DF-F)*((xsol-XL)*DF-F) &
!      .GE. 0.0 .OR. ABS(2.0*F) .GT. ABS(DXOLD*DF)) THEN
!      DXOLD=DX
!      DX=0.5*(XH-XL)
!      xsol=XL+DX
!      IF(XL .EQ. xsol) then
!        hplus = xsol
!        RETURN
!      end if
!    ELSE
!      DXOLD=DX
!      DX=F/DF
!      tmp=xsol
!      xsol=xsol-DX
!      IF(tmp .EQ. xsol)then
!        hplus = xsol
!        RETURN
!      end if
!    END IF
!    IF(ABS(DX) .LT. XACC) then
!      hplus = xsol
!      RETURN
!    end if	
!    CALL recom_talk_difference(xsol,F,DF)
!    IF(F .LT. 0.0) THEN
!      XL=xsol
!      FL=F
!    ELSE
!      XH=xsol
!      FH=F
!    END IF
!  END DO
	
!end subroutine REcoM_iterate_ph
	
!--------------------------------------------------------------------------------	
!subroutine REcoM_talk_difference(x,fn,df)
!  use REcoM_locVar
!  implicit none
	
!  Real(kind=8) :: x,fn,df
!  Real(kind=8) :: x2,x3,k12,b,b2,db,rb,rx
	
!! This routine expresses TA as a function of DIC, htotal and constants.
!! It also calculates the derivative of this function with respect to 
!! htotal. It is used in the iterative solution for htotal. In the call
!! "x" is the input value for htotal, "fn" is the calculated value for TA
!! and "df" is the value for dTA/dhtotal

!  x2  = x * x
!  x3  = x2 * x
!  k12 = k1 * k2
!  b   = x2 + k1 * x + k12
!  b2  = b * b
!  db  = 2.d0 * x + k1
!  rb  = 0.d0
!  if (b .ne. 0.d0) rb = 1.d0/b
!  rx  = 0.d0
!  if (x .ne. 0.d0) rx = 1.d0/x
	
!! function: fn = HCO3 + 2 CO3 + borate + OH - Hfree - TA
!  fn  = k1 * x * dic_molal * rb     &
!      + 2.d0 * dic_molal * k12 * rb &
!      + bt * kb /(kb + x)           &
!      + kw * rx                     &
!      - x - talk_molal              
     
!! df = dfn/dx
!	df  = k1 * dic_molal * rb - k1 * x * dic_molal * db * rb * rb &
!        - 2.d0 * dic_molal * k12 * db * rb * rb                   &       
!        - bt * kb/(kb + x)**2 - kw * rx * rx - 1.d0
     
!end subroutine REcoM_talk_difference

!================================================================================
! Subroutine controlling and calculating atm. dep. of Fe
!================================================================================
subroutine Atm_input(day_in_month, timeold, yearnew, cyearnew, month, recom_istep)
  use recom_ocean_settings
!  use REcoM_constants
  use REcoM_GloVar
  use REcoM_locVar
  use recom_config
  implicit none
  
#include "netcdf.inc"
 
  integer, INTENT(IN) :: day_in_month, month, yearnew, recom_istep
  character(4), INTENT(IN) :: cyearnew
  real(kind=8), INTENT(IN) :: timeold
 
  real(kind=8), allocatable :: ncdata(:)
  integer	                :: status, ncid, varid
  character(2000)            :: Ironname, CO2filename, DustNfilename
  character(20)             :: Ironvari, CO2vari, Nvari
  !character(10)             :: runoff_data_source='CORE2'
  integer, dimension(2)     :: istart, icount
  integer                   :: CO2start, CO2count
  integer                   :: firstyearofcurrentCO2cycle, totnumyear
  character(4)              :: currentCO2year_char
  
!-Checking if files need to be opened---------------------------------------------

  if(day_in_month==1 .and. mod(int(timeold),int(SecondsPerDay))==0) then ! file is opened and read every month
  
  if(UseFeDust) then
  
!-Testing if files exist for the year in question---------------------------------
  if (UseDustClimAlbani) then
    Ironname = trim(REcoMDataPath)//'DustClimMonthlyAlbani.nc'
    Ironvari     = 'DustClim'
  else  
    if((yearnew .LT. 1979) .OR.(yearnew .GT. 2008)) then
      Ironname = trim(REcoMDataPath)//'DustClimMonthlyMahowald.nc'
      Ironvari     = 'DustClim'
    else
      if (UseDustClim) then
        Ironname = trim(REcoMDataPath)//'DustClimMonthlyMahowald.nc'
        Ironvari     = 'DustClim'
      else
        Ironname = trim(REcoMDataPath)//'DustMonthlyMahowald.nc'
        Ironvari     = 'Dust'//cyearnew
      end if
    end if 
  end if    

!-Opening Iron--------------------------------------------------------------------
    
    ! open file
    status=nf_open(Ironname, nf_nowrite, ncid)
    if (status.ne.nf_noerr)then
      print*,'ERROR: CANNOT READ fe FILE CORRECTLY !!!!!'
      print*,'Error in opening netcdf file '//Ironname
      stop
    endif
    
    ! data
    allocate(ncdata(recom_nod2D))
    status=nf_inq_varid(ncid, Ironvari, varid)
    istart = (/1,month/)
    icount= (/recom_nod2D,1/)
    status=nf_get_vara_double(ncid,varid,istart,icount,ncdata)
    GloFeDust(:)=ncdata(recom_myList_nod2D)
    deallocate(ncdata)
    
    status=nf_close(ncid)
  endif ! FeDust
  endif ! time to open
    
!-Iron finished--------------------------------------------------------------------

!-Opening CO2----------------------------------------------------------------------
  if (recom_istep == 1) then ! The year has changed
    if (constant_CO2) then
      AtmCO2(:) = CO2_for_spinup
!     if (recom_mype==1) write(*,*),'in atm_input: Atm CO2=',AtmCO2

    else  

     CO2filename = trim(REcoMDataPath)//'MonthlyAtmCO2.nc'

     totnumyear                 = lastyearoffesomcycle-firstyearoffesomcycle+1
     firstyearofcurrentCO2cycle = lastyearoffesomcycle-numofCO2cycles*totnumyear+(currentCO2cycle-1)*totnumyear
    
     currentCO2year = firstyearofcurrentCO2cycle + (yearnew-firstyearoffesomcycle)+1
     if(recom_mype==0) write(*,*),currentCO2year, firstyearofcurrentCO2cycle, yearnew, firstyearoffesomcycle
     write(currentCO2year_char,'(i4)') currentCO2year
     CO2vari     = 'AtmCO2_'//currentCO2year_char
    
     ! open file
     status=nf_open(CO2filename, nf_nowrite, ncid)
     if (status.ne.nf_noerr)then
        print*,'ERROR: CANNOT READ CO2 FILE CORRECTLY !!!!!'
        print*,'Error in opening netcdf file '//CO2filename
        stop
     endif    
	
	! data
     allocate(ncdata(12))
     status=nf_inq_varid(ncid, CO2vari, varid)
     CO2start = 1
     CO2count = 12
     status=nf_get_vara_double(ncid,varid,CO2start,CO2count,ncdata)
     AtmCO2(:)=ncdata(:)
     deallocate(ncdata)
     if (recom_mype==0) write(*,*),'Current carbon year=',currentCO2year
     if (recom_mype==0) write(*,*),'Atm CO2=',AtmCO2(1)

    status=nf_close(ncid)
     
   endif  ! constant CO2 or changing
    
 ! Aeolian nitrogen deposition
   if (useAeolianN) then

     DustNfilename = trim(REcoMDataPath)//'AeolianNitrogenDep.nc'
     if (yearnew .lt. 2010) then
        Nvari         = 'NDep'//cyearnew
     else
        Nvari = 'NDep2009'
     endif
     ! open file
!     write(*,*),DustNfilename
     status=nf_open(DustNfilename, nf_nowrite, ncid)
     if (status.ne.nf_noerr)then
!        print*,'ERROR: CANNOT READ nitrogen dust FILE CORRECTLY !!!!!'
        if(recom_mype==0) write(*,*),'Error in opening netcdf file '//DustNfilename
        stop
     endif
     ! data
     allocate(ncdata(recom_nod2D))
     status=nf_inq_varid(ncid, Nvari, varid)
     istart = (/1,1/)
     icount= (/recom_nod2D,1/)
     status=nf_get_vara_double(ncid,varid,istart,icount,ncdata)
     GloNDust(:)=ncdata(recom_myList_nod2D)
     deallocate(ncdata)
   else
     GloNDust(:) = 0.d0 ! no aeolian N input 
   end if
  else
    return
  end if
  
end subroutine Atm_input


!================================================================================
! Calculating second zooplankton respiration rates
!================================================================================
 subroutine krill_resp(daynew,Latr)
   use recom_ocean_settings
   use REcoM_declarations
!  use REcoM_constants
   use REcoM_LocVar
   implicit none
 
   integer, INTENT(IN) :: daynew  
   ! Values from FESOM                                                                                                 
   Real(kind=8),intent(in)    :: Latr
 
  if ((Latr .lt. 0).and.(daynew .le. 105)) then
       res_zoo2_a = 0.d0
   else if((Latr .lt. 0).and.(105 .le. daynew).and.(daynew .le. 150)) then
       res_zoo2_a = (-1./90.)*daynew +7./6.
   else if((Latr .lt. 0).and.(150 .lt. daynew).and.(daynew .lt. 250)) then
       res_zoo2_a = -0.5
   else if((Latr .lt. 0).and.(250 .le. daynew).and.(daynew .le. 295)) then
       res_zoo2_a = (1/90.)*daynew - 59./18.
   else if((Latr .lt. 0).and.(daynew .gt. 295)) then
       res_zoo2_a = 0.d0
  end if

!!For Northern Hemisphere
  if ((Latr .ge. 0).and.(daynew .le. 65)) then
       res_zoo2_a = -0.5
   else if((Latr .ge. 0).and.(285 .le. daynew).and.(daynew .le. 330)) then
       res_zoo2_a = (-1./90.)*daynew +57./18.
   else if((Latr .ge. 0).and.(330 .lt. daynew)) then
       res_zoo2_a = -0.5
   else if((Latr .ge. 0).and.(65 .le. daynew).and.(daynew .le. 110)) then
       res_zoo2_a = (1/90.)*daynew - 22./18.
   else if((Latr .ge. 0).and.(110 .lt. daynew).and.(daynew .lt. 285)) then
       res_zoo2_a = 0.d0
  end if
 
 end subroutine krill_resp
