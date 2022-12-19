module recom_ocean_settings
!         use o_DATA_TYPES ! we wanted to avoid using fesom modules in recom, but I haven't found
         ! a way to avoid this one, as the same definition of a data type here as in fesom would 
         ! NOT be considered to be the same data type by the compiler. Either have to use a module 
         ! as done here or use 'SEQUENCE' but this would cause other trouble.
!        use MPI
        ! MPI_COMM_WORLD wird benutzt, sollte MPI_COMM_RECOM oder so sein...
        ! in FESOM ist der richtige comm MPI_COMM_FESOM

        ! Somewhere in fesom you need:
        !use g_parfe, only :  mype, mydim_nod2D, todim_nod2D, todim_nod3D, mydim_elem2D,&
        !                     mylist_nod2D, edim_nod2D, mydim_nod3D, &
        !                     edim_nod3D, mylist_nod3D, MPIerr, MPI_COMM_FESOM
        !use o_mesh, only : nod2D, coord_nod2D, nod3D, coord_nod3D, nod_in_elem3D, &
        !                   num_layers_below_nod2D, nod3D_below_nod2D, max_num_layers 
        !use o_array, only : prog_tracer_name, tracer
        !use o_param, only : pi, num_tracer, rad
        !use o_elements, only : voltetra, voltriangle, elem2D_nodes 
        !! istep needs to go out here, else it is never updated
        !use g_config, only : ResultPath, runid, dt, istep, save_count_restart, &
        !                     restartflag, MeshPath, rotated_grid
        !call recom_init_ocean_settings(mype, mype_nod2D, todim_nod2D, todim_nod3D, &
        !                mydim_elem2d, mylist_nod2D, edim_nod2D, &
        !                mydim_nod3D, edim_nod3D, mylist_nod3D, &
        !                nod2D, coord_nod2D, nod3D, coord_nod3D, &
        !                nod_in_elem3D, num_layers_below_nod2D, &
        !                nod3D_below_nod2D, max_num_layers, &
        !                prog_tracer_name, num_tracer, pi, rad, &
        !                voltetra, voltriangle, elem2D_nodes, &
        !                ResultPath, runid, dt, istep, &
        !                save_count_restart, restartflag, MeshPath, &
        !                rotated_grid, rotate_matrix, MPI_COMM_FESOM, MPIerr)


        ! calls to recom functions changed:

        !use g_forcing_arrays, only : shortwave, u_wind, v_wind 
        !use i_array, only: a_ice
        !use o_array, only: tracer
        !use g_clock, only: month, yearnew, yearold, cyearnew, daynew, ndpyr,&
        !                   day_in_month, timeold
        !call call_REcoM(shortwave, u_wind, v_wind, a_ice, tracer, month,&
        !                yearnew, yearold, cyearnew, daynew, ndpyr, &
        !                day_in_month, timeold)

        !use o_array, only: tracer
        !use g_clock, only: cyearnew
        !call save_recom_binary(tracer, cyearnew)

        !use g_clock, only: cyearnew
        !call REcoM_write_snapshots(cyearnew)

        !use g_clock, only: yearnew, yearold, cyearnew
        !call REcoM_init_output(yearnew, yearold, cyearnew)

        !use o_array, only: tracer
        !use g_clock, only: yearnew, yearold, cyearnew, cyearold
        !call recom_init(tracer, yearnew, cyearnew, yearold, cyearold)

        !use o_array, only: tracer
        !use g_clock, only: cyearold
        !call read_recom_binary(tracer, cyearold)

        !use g_clock, only: daynew, ndpyr
        !call Cobeta(cosAngleOfIncidence,Latr, daynew, ndpyr)

        !use g_clock, only: daynew, ndpyr
        !call REcoM_Forcing(Latr,zNodes,Nn,state,SurfSW,dt_s,Temp,Sali,PAR,daynew,&
        !                 ndpyr)
   
        !use g_clock, only: day_in_month, timeold, yearnew, cyearnew, month
        !call Atm_input(day_in_month, timeold, yearnew, cyearnew, month)

!        see comment at top of routine, has to be loaded with a module to be considered the same 
        type addresstype
                integer                                :: nmb
                integer(KIND=4), dimension(:), pointer :: addresses
        end type addresstype

        save

        integer :: recom_mype, recom_mydim_nod2D, recom_todim_nod2D, &
                   recom_mydim_elem2d, recom_edim_nod2D, recom_mydim_nod3D, &
                   recom_edim_nod3D, recom_todim_nod3D, recom_mydim_elem3D

        integer :: recom_nod2D, recom_nod3D, recom_max_num_layers 

        integer :: MPI_COMM_RECOM, MPIerr_recom

        integer, allocatable, dimension(:) :: recom_mylist_nod2D 
        integer, allocatable, dimension(:) :: recom_mylist_nod3D

        integer, allocatable, dimension(:)           :: recom_num_layers_below_nod2D
        integer(KIND=4), allocatable, dimension(:,:) :: recom_nod3D_below_nod2D
        integer(KIND=4), allocatable, dimension(:,:) :: recom_elem3D_nodes

        integer :: recom_num_tracer
        integer :: recom_save_count_restart
        integer(KIND=4), allocatable, dimension(:,:) :: recom_elem2D_nodes

        real(kind=8), allocatable, dimension(:,:)    :: recom_coord_nod2D, recom_coord_nod3D
        real(kind=8) :: fesom_pi, recom_rad
        real(kind=8)                    :: recom_dt
        real(kind=8), allocatable, dimension(:)         :: recom_voltetra, recom_voltriangle
        real(kind=8)        :: recom_rotate_matrix(3,3)

        character(4), allocatable, dimension(:)         :: recom_prog_tracer_name
        character(2000)                :: recom_ResultPath, recom_MeshPath
        character(5)                    :: recom_runid
        character*4                     :: recom_restartflag

        logical                         :: recom_rotated_grid

        type(addresstype), allocatable, dimension(:) :: recom_nod_in_elem3D

        contains

          subroutine r2g(lon, lat, rlon, rlat)
                ! Convert the rotated coordinates to geographical coordinates  
                ! lon, lat          :: [radian] geographical coordinates
                ! rlon, rlat        :: [radian] rotated coordinates
                !
                implicit none
                real(kind=8), intent(out)      :: lon, lat
                real(kind=8), intent(in)       :: rlon, rlat
                real(kind=8)                   :: xr, yr, zr, xg, yg, zg
                !
                ! Rotated Cartesian coordinates:
                xr=cos(rlat)*cos(rlon)
                yr=cos(rlat)*sin(rlon)
                zr=sin(rlat)

                ! Geographical Cartesian coordinates:
                xg=recom_rotate_matrix(1,1)*xr + recom_rotate_matrix(2,1)*yr + recom_rotate_matrix(3,1)*zr
                yg=recom_rotate_matrix(1,2)*xr + recom_rotate_matrix(2,2)*yr + recom_rotate_matrix(3,2)*zr  
                zg=recom_rotate_matrix(1,3)*xr + recom_rotate_matrix(2,3)*yr + recom_rotate_matrix(3,3)*zr  

                ! Geographical coordinates:
                lat=asin(zg)
                if(yg==0. .and. xg==0.) then
                        lon=0.0     ! exactly at the poles
                else
                        lon=atan2(yg,xg)
                end if
        end subroutine r2g



        subroutine recom_init_ocean_settings(in_mype, in_mydim_nod2D, in_todim_nod2D, &
                        in_todim_nod3D, in_elem3D_nodes, in_mydim_elem3D, &
                        in_mydim_elem2d, in_mylist_nod2D, in_edim_nod2D, &
                        in_mydim_nod3D, in_edim_nod3D, in_mylist_nod3D, &
                        in_nod2D, in_coord_nod2D, in_nod3D, in_coord_nod3D, &
                        in_num_layers_below_nod2D, & !in_nod_in_elem3D, & 
                        in_nod3D_below_nod2D, in_max_num_layers, &
                        in_prog_tracer_name, in_num_tracer, in_pi, in_rad, &
                        in_voltetra, in_voltriangle, in_elem2D_nodes, &
                        in_ResultPath, in_runid, in_dt, &
                        in_save_count_restart, in_restartflag, in_MeshPath, &
                        in_rotated_grid, in_rotate_matrix, in_MPI_COMM_RECOM, &
                        in_MPIerr)
                 use MPI
                 use recom_config

                INTEGER, INTENT(IN) :: in_mype, in_mydim_nod2D, in_todim_nod2D, &
                        in_mydim_elem2d, in_edim_nod2D, in_mydim_elem3D, &
                        in_mydim_nod3D, in_edim_nod3D, in_todim_nod3D, &
                        in_nod2D, in_nod3D, &
                        in_max_num_layers, in_num_tracer,&
                        in_save_count_restart, &
                        in_MPI_COMM_RECOM,  in_MPIerr
                integer, dimension(:), INTENT(IN) :: in_myList_nod2D
                integer, dimension(:), INTENT(IN) :: in_myList_nod3D

                integer(KIND=4), dimension(:,:) :: in_elem2D_nodes
                integer, dimension(:) :: in_num_layers_below_nod2D
                integer(KIND=4), dimension(:,:) :: in_nod3D_below_nod2D
                integer(KIND=4), dimension(:,:) :: in_elem3D_nodes

                REAL(KIND=8), INTENT(IN) :: in_pi, in_rad, in_dt

                REAL(KIND=8), DIMENSION(:), INTENT(IN) :: in_voltetra, &
                        in_voltriangle

                real(kind=8), dimension(:,:)    :: in_coord_nod2D, &
                                                   in_coord_nod3D
                real(kind=8)        :: in_rotate_matrix(3,3)

                CHARACTER(4), DIMENSION(:), INTENT(IN) :: in_prog_tracer_name
                character(2000)                :: in_ResultPath, in_MeshPath
                character(5)                   :: in_runid
                character*4                     :: in_restartflag

                logical                         :: in_rotated_grid

!                type(addresstype), dimension(:) :: in_nod_in_elem3D

                integer                            :: j,k,tet(4)
                integer, allocatable, dimension(:) :: ind


               if (.not. use_REcoM) return

                recom_mype=in_mype
                recom_mydim_nod2D=in_mydim_nod2D
                recom_todim_nod2D=in_todim_nod2D
                recom_todim_nod3D=in_todim_nod3D
                recom_mydim_elem2d=in_mydim_elem2d
                recom_mydim_elem3d=in_mydim_elem3d
                recom_mylist_nod2D=in_mylist_nod2D
                recom_edim_nod2D=in_edim_nod2D
                recom_mydim_nod3D=in_mydim_nod3D
                recom_edim_nod3D=in_edim_nod3D
                recom_mylist_nod3D=in_mylist_nod3D
                recom_MPIerr=in_MPIerr
                recom_nod2D=in_nod2D
                recom_nod3D=in_nod3D
                recom_max_num_layers=in_max_num_layers
                recom_num_tracer=in_num_tracer
                fesom_pi=in_pi
                recom_rad=in_rad
                recom_ResultPath=in_ResultPath
                recom_runid=in_runid
                recom_dt=in_dt
                recom_save_count_restart=in_save_count_restart
                recom_restartflag=in_restartflag
                recom_MeshPath=in_MeshPath
                recom_rotated_grid=in_rotated_grid
                recom_rotate_matrix=in_rotate_matrix
                recom_MPI_COMM_RECOM=in_MPI_COMM_RECOM

!                allocate(recom_mylist_nod2D(recom_mydim_nod2D+recom_edim_nod2D))
                recom_mylist_nod2D = in_mylist_nod2D
!                ALLOCATE(prog_tracer_name(num_tracer))
                recom_prog_tracer_name=in_prog_tracer_name
!                ALLOCATE(voltetra(SIZE(in_voltetra)))
                recom_voltetra=in_voltetra
!                ALLOCATE(voltriangle(SIZE(in_voltriangle)))
                recom_voltriangle=in_voltriangle
!                allocate(elem2D_nodes(3, myDim_elem2D))
                recom_elem2D_nodes=in_elem2D_nodes
!                allocate(coord_nod2D(2, myDim_nod2D+eDim_nod2D))
                recom_coord_nod2D=in_coord_nod2D
!                allocate(num_layers_below_nod2D(myDim_nod2D+eDim_nod2D))
                recom_num_layers_below_nod2D=in_num_layers_below_nod2D
!                allocate(nod3D_below_nod2D(max_num_layers,myDim_nod2D+eDim_nod2D))
                recom_nod3D_below_nod2D=in_nod3D_below_nod2D
!                allocate(coord_nod3D(3,myDim_nod3D+eDim_nod3D)) 
                recom_coord_nod3D=in_coord_nod3D
!                allocate(recom_nod_in_elem3D(recom_myDim_nod3D))
!                recom_nod_in_elem3D=in_nod_in_elem3D
                recom_elem3D_nodes = in_elem3D_nodes

! now have to repeat the definition of nod_in_elem3D because user-derived data
! type cannot be passed to a subroutine (easily)

                ! Construction of nod_in_elem3D
                allocate(ind(recom_myDim_nod3D+recom_eDim_nod3D))
                ind=0
                do j=1,recom_myDim_elem3D
                   tet=recom_elem3D_nodes(:,j)
                   ind(tet)=ind(tet)+1
                end do
                allocate(recom_nod_in_elem3D(recom_myDim_nod3D))
                recom_nod_in_elem3D%nmb=ind(1:recom_myDim_nod3D)
                do j=1,recom_myDim_nod3D   
                   allocate(recom_nod_in_elem3D(j)%addresses(ind(j)))
                end do
                ind=0
                do j=1,recom_myDim_elem3D   
                   tet=recom_elem3D_nodes(:,j)
                   ind(tet)=ind(tet)+1
                   do k=1,4
                      if(tet(k)<=recom_myDim_nod3D) then                         
                         recom_nod_in_elem3D(tet(k))%addresses(ind(tet(k)))=j     
                      end if
                   end do
                 end do

                 deallocate(ind)
  ! the list of elements is ordered, and no sorting is needed

        end subroutine recom_init_ocean_settings

end module recom_ocean_settings

