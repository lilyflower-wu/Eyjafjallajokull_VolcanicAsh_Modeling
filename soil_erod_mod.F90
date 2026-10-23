!===============================================================================
!===============================================================================
module soil_erod_mod
  use shr_kind_mod,     only: r8 => shr_kind_r8, cl => shr_kind_cl
  use cam_logfile,      only: iulog
  use spmd_utils,       only: masterproc
  use cam_abortutils,   only: endrun

  implicit none
  private

  public :: soil_erod_init
  public :: soil_erodibility
  public :: soil_erod_fact

#if ( defined DUSTTYPES )
  public ::  illite_frac_accum
  public ::  illite_frac_aitken
  public ::  illite_frac_coarse
  public ::  kaolinite_frac_accum
  public ::  kaolinite_frac_aitken
  public ::  kaolinite_frac_coarse
  public ::  mont_frac_accum
  public ::  mont_frac_aitken
  public ::  mont_frac_coarse
  public ::  hematite_frac_accum
  public ::  hematite_frac_aitken
  public ::  hematite_frac_coarse
  public ::  quartz_frac_accum
  public ::  quartz_frac_aitken
  public ::  quartz_frac_coarse
  public ::  calcite_frac_accum
  public ::  calcite_frac_aitken
  public ::  calcite_frac_coarse
  public ::  feldspar_frac_accum
  public ::  feldspar_frac_aitken
  public ::  feldspar_frac_coarse
  public ::  gypsum_frac_accum
  public ::  gypsum_frac_aitken
  public ::  gypsum_frac_coarse
!--++Longlei--
  public ::  illite_clay
  public ::  kaolinite_clay
  public ::  kaolinite_silt
  public ::  mont_clay
  public ::  hematite_clay
  public ::  hematite_silt
  public ::  quartz_clay
  public ::  quartz_silt
  public ::  calcite_clay
  public ::  calcite_silt
  public ::  feldspar_clay
  public ::  feldspar_silt
  public ::  gypsum_silt
#if ( defined IRONPROC )
  public ::  comb_accum_sol
  public ::  comb_aitken_sol
  public ::  comb_coarse_sol
  public ::  comb_accum_insol
  public ::  comb_aitken_insol
  public ::  comb_coarse_insol
#endif
#endif

  real(r8), allocatable ::  soil_erodibility(:,:)  ! soil erodibility factor
#if ( defined DUSTTYPES )
  real(r8), allocatable ::  illite_frac_accum(:,:)
  real(r8), allocatable ::  illite_frac_aitken(:,:)
  real(r8), allocatable ::  illite_frac_coarse(:,:)
  real(r8), allocatable ::  kaolinite_frac_accum(:,:)
  real(r8), allocatable ::  kaolinite_frac_aitken(:,:)
  real(r8), allocatable ::  kaolinite_frac_coarse(:,:)
  real(r8), allocatable ::  mont_frac_accum(:,:)
  real(r8), allocatable ::  mont_frac_aitken(:,:)
  real(r8), allocatable ::  mont_frac_coarse(:,:)
  real(r8), allocatable ::  hematite_frac_accum(:,:)
  real(r8), allocatable ::  hematite_frac_aitken(:,:)
  real(r8), allocatable ::  hematite_frac_coarse(:,:)
  real(r8), allocatable ::  quartz_frac_accum(:,:)
  real(r8), allocatable ::  quartz_frac_aitken(:,:)
  real(r8), allocatable ::  quartz_frac_coarse(:,:)
  real(r8), allocatable ::  calcite_frac_accum(:,:)
  real(r8), allocatable ::  calcite_frac_aitken(:,:)
  real(r8), allocatable ::  calcite_frac_coarse(:,:)
  real(r8), allocatable ::  feldspar_frac_accum(:,:)
  real(r8), allocatable ::  feldspar_frac_aitken(:,:)
  real(r8), allocatable ::  feldspar_frac_coarse(:,:)
  real(r8), allocatable ::  gypsum_frac_accum(:,:)
  real(r8), allocatable ::  gypsum_frac_aitken(:,:)
  real(r8), allocatable ::  gypsum_frac_coarse(:,:)
!--++Longlei--
  real(r8), allocatable ::  illite_clay(:,:)
  real(r8), allocatable ::  kaolinite_clay(:,:)
  real(r8), allocatable ::  kaolinite_silt(:,:)
  real(r8), allocatable ::  mont_clay(:,:)
  real(r8), allocatable ::  hematite_clay(:,:)
  real(r8), allocatable ::  hematite_silt(:,:)
  real(r8), allocatable ::  quartz_clay(:,:)
  real(r8), allocatable ::  quartz_silt(:,:)
  real(r8), allocatable ::  calcite_clay(:,:)
  real(r8), allocatable ::  calcite_silt(:,:)
  real(r8), allocatable ::  feldspar_clay(:,:)
  real(r8), allocatable ::  feldspar_silt(:,:)
  real(r8), allocatable ::  gypsum_silt(:,:)
#if ( defined IRONPROC )
  real(r8), allocatable ::  comb_accum_sol(:,:)
  real(r8), allocatable ::  comb_aitken_sol(:,:)
  real(r8), allocatable ::  comb_coarse_sol(:,:)
  real(r8), allocatable ::  comb_accum_insol(:,:)
  real(r8), allocatable ::  comb_aitken_insol(:,:)
  real(r8), allocatable ::  comb_coarse_insol(:,:)
  !real(r8), pointer, dimension(:) :: curr_data_times => null()
  !real(r8), pointer, dimension(:) :: next_data_times => null()
  integer :: curr_data_times 
  real(r8) :: next_data_times 

#endif
#endif
  real(r8) :: soil_erod_fact                       ! tuning parameter for dust emissions

contains

  !=============================================================================
  !=============================================================================
  subroutine soil_erod_init( dust_emis_fact, soil_erod_file, comb_emis_file, comb_file_list )
    use interpolate_data, only: lininterp_init, lininterp, lininterp_finish, interp_type
    use ppgrid,           only: begchunk, endchunk, pcols
    use mo_constants,     only: pi, d2r
    use pio,              only: file_desc_t,pio_inq_dimid,pio_inq_dimlen,pio_get_var,pio_inq_varid, PIO_NOWRITE
    use phys_grid,        only: get_ncols_p, get_rlat_all_p, get_rlon_all_p
    use cam_pio_utils,    only: cam_pio_openfile
    use ioFileMod,        only: getfil
    ! Longlei Li
    use cam_control_mod, only : restart_run
    !use datetime_mod, only: datetime

    real(r8),         intent(in) :: dust_emis_fact
    character(len=*), intent(in) :: soil_erod_file
    character(len=*), intent(in) :: comb_emis_file
    character(len=*), intent(in) :: comb_file_list 

    real(r8), allocatable :: soil_erodibility_in(:,:)  ! temporary input array
    real(r8), allocatable :: dst_lons(:)
    real(r8), allocatable :: dst_lats(:)
#if ( defined DUSTTYPES )
    real(r8), allocatable :: illite_frac_in_accum(:,:)
    real(r8), allocatable :: illite_frac_in_aitken(:,:)
    real(r8), allocatable :: illite_frac_in_coarse(:,:)
    real(r8), allocatable :: kaolinite_frac_in_accum(:,:)
    real(r8), allocatable :: kaolinite_frac_in_aitken(:,:)
    real(r8), allocatable :: kaolinite_frac_in_coarse(:,:)
    real(r8), allocatable :: mont_frac_in_accum(:,:)
    real(r8), allocatable :: mont_frac_in_aitken(:,:)
    real(r8), allocatable :: mont_frac_in_coarse(:,:)
    real(r8), allocatable :: hematite_frac_in_accum(:,:)
    real(r8), allocatable :: hematite_frac_in_aitken(:,:)
    real(r8), allocatable :: hematite_frac_in_coarse(:,:)
    real(r8), allocatable :: quartz_frac_in_accum(:,:)
    real(r8), allocatable :: quartz_frac_in_aitken(:,:)
    real(r8), allocatable :: quartz_frac_in_coarse(:,:)
    real(r8), allocatable :: calcite_frac_in_accum(:,:)
    real(r8), allocatable :: calcite_frac_in_aitken(:,:)
    real(r8), allocatable :: calcite_frac_in_coarse(:,:)
    real(r8), allocatable :: feldspar_frac_in_accum(:,:)
    real(r8), allocatable :: feldspar_frac_in_aitken(:,:)
    real(r8), allocatable :: feldspar_frac_in_coarse(:,:)
    real(r8), allocatable :: gypsum_frac_in_accum(:,:)
    real(r8), allocatable :: gypsum_frac_in_aitken(:,:)
    real(r8), allocatable :: gypsum_frac_in_coarse(:,:)
    real(r8), allocatable ::  illite_in_clay(:,:)
    real(r8), allocatable ::  kaolinite_in_clay(:,:)
    real(r8), allocatable ::  kaolinite_in_silt(:,:)
    real(r8), allocatable ::  mont_in_clay(:,:)
    real(r8), allocatable ::  hematite_in_clay(:,:)
    real(r8), allocatable ::  hematite_in_silt(:,:)
    real(r8), allocatable ::  quartz_in_clay(:,:)
    real(r8), allocatable ::  quartz_in_silt(:,:)
    real(r8), allocatable ::  calcite_in_clay(:,:)
    real(r8), allocatable ::  calcite_in_silt(:,:)
    real(r8), allocatable ::  feldspar_in_clay(:,:)
    real(r8), allocatable ::  feldspar_in_silt(:,:)
    real(r8), allocatable ::  gypsum_in_silt(:,:)
#if ( defined IRONPROC )
    real(r8), allocatable ::  comb_in_accum_sol(:,:)
    real(r8), allocatable ::  comb_in_aitken_sol(:,:)
    real(r8), allocatable ::  comb_in_coarse_sol(:,:)
    real(r8), allocatable ::  comb_in_accum_insol(:,:)
    real(r8), allocatable ::  comb_in_aitken_insol(:,:)
    real(r8), allocatable ::  comb_in_coarse_insol(:,:)
    real(r8), allocatable ::  comb_lons(:)
    real(r8), allocatable ::  comb_lats(:)

    integer :: nlatc, nlonc
    type(file_desc_t) :: ncid2
    type(interp_type) :: lon_wgtsc, lat_wgtsc
#endif
#endif
    character(len=cl)     :: infile
    integer :: did, vid, nlat, nlon
    type(file_desc_t) :: ncid

    type(interp_type) :: lon_wgts, lat_wgts
    real(r8) :: to_lats(pcols), to_lons(pcols)
    integer :: c, ncols, ierr
    real(r8), parameter :: zero=0._r8, twopi=2._r8*pi

    soil_erod_fact = dust_emis_fact

    ! Summary to log file
    if (masterproc) then
       write(iulog,*) 'soil_erod_mod: comb iron emission dataset: ', trim(comb_emis_file)
       !write(iulog,*) 'soil_erod_mod: comb iron emission dataset filename list: ', trim(comb_file_list)
       write(iulog,*) 'soil_erod_mod: soil erodibility dataset: ', trim(soil_erod_file)
       write(iulog,*) 'soil_erod_mod: soil_erod_fact = ', soil_erod_fact
    end if

    ! for soil erodibility in mobilization, apply inside CAM instead of lsm.
    ! read in soil erodibility factors, similar to Zender's boundary conditions

    ! Get file name.  
    call getfil(soil_erod_file, infile, 0)
    call cam_pio_openfile (ncid, trim(infile), PIO_NOWRITE)

    ! Get input data resolution.
    ierr = pio_inq_dimid( ncid, 'lon', did )
    ierr = pio_inq_dimlen( ncid, did, nlon )

    ierr = pio_inq_dimid( ncid, 'lat', did )
    ierr = pio_inq_dimlen( ncid, did, nlat )

    allocate(dst_lons(nlon))
    allocate(dst_lats(nlat))
    allocate(soil_erodibility_in(nlon,nlat))

#if ( defined DUSTTYPES )
    allocate(illite_frac_in_accum(nlon,nlat))
    allocate(illite_frac_in_aitken(nlon,nlat))
    allocate(illite_frac_in_coarse(nlon,nlat))

    allocate(kaolinite_frac_in_accum(nlon,nlat))
    allocate(kaolinite_frac_in_aitken(nlon,nlat))
    allocate(kaolinite_frac_in_coarse(nlon,nlat))

    allocate(mont_frac_in_accum(nlon,nlat))
    allocate(mont_frac_in_aitken(nlon,nlat))
    allocate(mont_frac_in_coarse(nlon,nlat))

    allocate(hematite_frac_in_accum(nlon,nlat))
    allocate(hematite_frac_in_aitken(nlon,nlat))
    allocate(hematite_frac_in_coarse(nlon,nlat))

    allocate(quartz_frac_in_accum(nlon,nlat))
    allocate(quartz_frac_in_aitken(nlon,nlat))
    allocate(quartz_frac_in_coarse(nlon,nlat))

    allocate(calcite_frac_in_accum(nlon,nlat))
    allocate(calcite_frac_in_aitken(nlon,nlat))
    allocate(calcite_frac_in_coarse(nlon,nlat))

    allocate(feldspar_frac_in_accum(nlon,nlat))
    allocate(feldspar_frac_in_aitken(nlon,nlat))
    allocate(feldspar_frac_in_coarse(nlon,nlat))

    allocate(gypsum_frac_in_accum(nlon,nlat))
    allocate(gypsum_frac_in_aitken(nlon,nlat))
    allocate(gypsum_frac_in_coarse(nlon,nlat))

    allocate(illite_in_clay(nlon,nlat))
    allocate(kaolinite_in_clay(nlon,nlat))
    allocate(kaolinite_in_silt(nlon,nlat))
    allocate(mont_in_clay(nlon,nlat))
    allocate(hematite_in_clay(nlon,nlat))
    allocate(hematite_in_silt(nlon,nlat))
    allocate(quartz_in_clay(nlon,nlat))
    allocate(quartz_in_silt(nlon,nlat))
    allocate(calcite_in_clay(nlon,nlat))
    allocate(calcite_in_silt(nlon,nlat))
    allocate(feldspar_in_clay(nlon,nlat))
    allocate(feldspar_in_silt(nlon,nlat))
    allocate(gypsum_in_silt(nlon,nlat))
#endif

    ierr = pio_inq_varid( ncid, 'lon', vid )
    ierr = pio_get_var( ncid, vid, dst_lons  )

    ierr = pio_inq_varid( ncid, 'lat', vid )
    ierr = pio_get_var( ncid, vid, dst_lats  )

    ierr = pio_inq_varid( ncid, 'mbl_bsn_fct_geo', vid )
    ierr = pio_get_var( ncid, vid, soil_erodibility_in )

#if ( defined DUSTTYPES )
    ierr = pio_inq_varid( ncid, 'FracI', vid )            !FracIa
    ierr = pio_get_var( ncid, vid, illite_in_clay )

    ierr = pio_inq_varid( ncid, 'FracK1', vid )            !FracKa
    ierr = pio_get_var( ncid, vid, kaolinite_in_clay )
    ierr = pio_inq_varid( ncid, 'FracK2', vid )            !FracKa
    ierr = pio_get_var( ncid, vid, kaolinite_in_silt )

    ierr = pio_inq_varid( ncid, 'FracS', vid )            !FracSa
    ierr = pio_get_var( ncid, vid, mont_in_clay )

    ierr = pio_inq_varid( ncid, 'FracH1', vid )            !FracHa
    ierr = pio_get_var( ncid, vid, hematite_in_clay )
    ierr = pio_inq_varid( ncid, 'FracH2', vid )            !FracHc
    ierr = pio_get_var( ncid, vid, hematite_in_silt )

    ierr = pio_inq_varid( ncid, 'FracQ1', vid )            !FracQa
    ierr = pio_get_var( ncid, vid, quartz_in_clay )
    ierr = pio_inq_varid( ncid, 'FracQ2', vid )            !FracQc
    ierr = pio_get_var( ncid, vid, quartz_in_silt )

    ierr = pio_inq_varid( ncid, 'FracC1', vid )            !FracCa
    ierr = pio_get_var( ncid, vid, calcite_in_clay )
    ierr = pio_inq_varid( ncid, 'FracC2', vid )            !FracCc
    ierr = pio_get_var( ncid, vid, calcite_in_silt )

    ierr = pio_inq_varid( ncid, 'FracF1', vid )            !FracFa
    ierr = pio_get_var( ncid, vid, feldspar_in_clay )
    ierr = pio_inq_varid( ncid, 'FracF2', vid )            !FracFa
    ierr = pio_get_var( ncid, vid, feldspar_in_silt )

    ierr = pio_inq_varid( ncid, 'FracG', vid )            !FracFc
    ierr = pio_get_var( ncid, vid, gypsum_in_silt )

#if ( defined IRONPROC )
    call getfil(comb_emis_file, infile, 0)

    call cam_pio_openfile (ncid, trim(infile), PIO_NOWRITE)

    ! Get input data resolution.
    ierr = pio_inq_dimid( ncid, 'lon', did )
    ierr = pio_inq_dimlen( ncid, did, nlonc )

    ierr = pio_inq_dimid( ncid, 'lat', did )
    ierr = pio_inq_dimlen( ncid, did, nlatc )

    allocate(comb_lons(nlonc))
    allocate(comb_lats(nlatc))

    allocate(comb_in_accum_sol(nlonc,nlatc))
    allocate(comb_in_aitken_sol(nlonc,nlatc))
    allocate(comb_in_coarse_sol(nlonc,nlatc))
    allocate(comb_in_accum_insol(nlonc,nlatc))
    allocate(comb_in_aitken_insol(nlonc,nlatc))
    allocate(comb_in_coarse_insol(nlonc,nlatc))

    ierr = pio_inq_varid( ncid, 'lon', vid )
    ierr = pio_get_var( ncid, vid, comb_lons  )
    ierr = pio_inq_varid( ncid, 'lat', vid )
    ierr = pio_get_var( ncid, vid, comb_lats  )

    ierr = pio_inq_varid( ncid, 'Fine_InsolFe_Comb', vid )       !Combustion fine insoluble   
    ierr = pio_get_var( ncid, vid, comb_in_accum_insol )
    ierr = pio_inq_varid( ncid, 'Coar_InsolFe_Comb', vid )       !Combustion coarse insoluble
    ierr = pio_get_var( ncid, vid, comb_in_coarse_insol )
    ierr = pio_inq_varid( ncid, 'Fine_SolFe_Comb', vid )         !Combustion fine soluble    
    ierr = pio_get_var( ncid, vid, comb_in_accum_sol )
    ierr = pio_inq_varid( ncid, 'Coar_SolFe_Comb', vid )         !Combustion coarse soluble   
    ierr = pio_get_var( ncid, vid, comb_in_coarse_sol )
!---new---

! partitioning of fine mode combustion fe emissions to aiken and accumulation modes
    comb_in_aitken_sol = comb_in_accum_sol *0.1_r8
    comb_in_accum_sol  = comb_in_accum_sol *0.9_r8
    comb_in_aitken_insol = comb_in_accum_insol *0.1_r8
    comb_in_accum_insol  = comb_in_accum_insol *0.9_r8
#endif
#endif
    !-----------------------------------------------------------------------
    !     	... convert to radians and setup regridding
    !-----------------------------------------------------------------------
    dst_lats(:) = d2r * dst_lats(:)
    dst_lons(:) = d2r * dst_lons(:)

#if ( defined IRONPROC )
    comb_lats(:) = d2r * comb_lats(:)
    comb_lons(:) = d2r * comb_lons(:)
#endif

    allocate( soil_erodibility(pcols,begchunk:endchunk), stat=ierr )
#if ( defined DUSTTYPES )
    allocate( illite_clay(pcols,begchunk:endchunk) )
    allocate( kaolinite_clay(pcols,begchunk:endchunk) )
    allocate( kaolinite_silt(pcols,begchunk:endchunk) )
    allocate( mont_clay(pcols,begchunk:endchunk) )
    allocate( hematite_clay(pcols,begchunk:endchunk) )
    allocate( hematite_silt(pcols,begchunk:endchunk) )
    allocate( quartz_clay(pcols,begchunk:endchunk) )
    allocate( quartz_silt(pcols,begchunk:endchunk) )
    allocate( calcite_clay(pcols,begchunk:endchunk) )
    allocate( calcite_silt(pcols,begchunk:endchunk) )
    allocate( feldspar_clay(pcols,begchunk:endchunk) )
    allocate( feldspar_silt(pcols,begchunk:endchunk) )
    allocate( gypsum_silt(pcols,begchunk:endchunk) )
#if ( defined IRONPROC )
    allocate( comb_accum_sol(pcols,begchunk:endchunk) )
    allocate( comb_aitken_sol(pcols,begchunk:endchunk) )
    allocate( comb_coarse_sol(pcols,begchunk:endchunk) )
    allocate( comb_accum_insol(pcols,begchunk:endchunk) )
    allocate( comb_aitken_insol(pcols,begchunk:endchunk) )
    allocate( comb_coarse_insol(pcols,begchunk:endchunk) )
#endif
#endif

    if( ierr /= 0 ) then
       write(iulog,*) 'soil_erod_init: failed to allocate soil_erodibility_in, ierr = ',ierr
       call endrun('soil_erod_init: failed to allocate soil_erodibility_in')
    end if

    !-----------------------------------------------------------------------
    !     	... regrid ..
    !-----------------------------------------------------------------------
    do c=begchunk,endchunk
       ncols = get_ncols_p(c)
       call get_rlat_all_p(c, pcols, to_lats)
       call get_rlon_all_p(c, pcols, to_lons)

       call lininterp_init(dst_lons, nlon, to_lons, ncols, 2, lon_wgts, zero, twopi)
       call lininterp_init(dst_lats, nlat, to_lats, ncols, 1, lat_wgts)

       call lininterp(soil_erodibility_in(:,:), nlon,nlat , soil_erodibility(:,c), ncols, lon_wgts,lat_wgts)

#if ( defined DUSTTYPES )
       call lininterp(illite_in_clay(:,:),    nlon, nlat, illite_clay(:,c),    ncols, lon_wgts, lat_wgts)
       call lininterp(kaolinite_in_clay(:,:), nlon, nlat, kaolinite_clay(:,c), ncols, lon_wgts, lat_wgts)
       call lininterp(kaolinite_in_silt(:,:), nlon, nlat, kaolinite_silt(:,c), ncols, lon_wgts, lat_wgts)
       call lininterp(mont_in_clay(:,:),      nlon, nlat, mont_clay(:,c),      ncols, lon_wgts, lat_wgts)
       call lininterp(hematite_in_clay(:,:),  nlon, nlat, hematite_clay(:,c),  ncols, lon_wgts, lat_wgts)
       call lininterp(hematite_in_silt(:,:),  nlon, nlat, hematite_silt(:,c),  ncols, lon_wgts, lat_wgts)
       call lininterp(quartz_in_clay(:,:),    nlon, nlat, quartz_clay(:,c),    ncols, lon_wgts, lat_wgts)
       call lininterp(quartz_in_silt(:,:),    nlon, nlat, quartz_silt(:,c),    ncols, lon_wgts, lat_wgts)
       call lininterp(calcite_in_clay(:,:),   nlon, nlat, calcite_clay(:,c),   ncols, lon_wgts, lat_wgts)
       call lininterp(calcite_in_silt(:,:),   nlon, nlat, calcite_silt(:,c),   ncols, lon_wgts, lat_wgts)
       call lininterp(feldspar_in_clay(:,:),  nlon, nlat, feldspar_clay(:,c),  ncols, lon_wgts, lat_wgts)
       call lininterp(feldspar_in_silt(:,:),  nlon, nlat, feldspar_silt(:,c),  ncols, lon_wgts, lat_wgts)
       call lininterp(gypsum_in_silt(:,:),    nlon, nlat, gypsum_silt(:,c),    ncols, lon_wgts, lat_wgts)
#if ( defined IRONPROC )
       call lininterp_init(comb_lons, nlonc, to_lons, ncols, 2, lon_wgtsc, zero, twopi)
       call lininterp_init(comb_lats, nlatc, to_lats, ncols, 1, lat_wgtsc)

       call lininterp(comb_in_accum_sol(:,:) , nlon, nlat, comb_accum_sol(:,c) , ncols, lon_wgts, lat_wgts)
       call lininterp(comb_in_aitken_sol(:,:), nlon, nlat, comb_aitken_sol(:,c), ncols, lon_wgts, lat_wgts)
       call lininterp(comb_in_coarse_sol(:,:), nlon, nlat, comb_coarse_sol(:,c), ncols, lon_wgts, lat_wgts)
       call lininterp(comb_in_accum_insol(:,:) , nlon, nlat, comb_accum_insol(:,c) , ncols, lon_wgts, lat_wgts)
       call lininterp(comb_in_aitken_insol(:,:), nlon, nlat, comb_aitken_insol(:,c), ncols, lon_wgts, lat_wgts)
       call lininterp(comb_in_coarse_insol(:,:), nlon, nlat, comb_coarse_insol(:,c), ncols, lon_wgts, lat_wgts)
#endif
#endif
       call lininterp_finish(lat_wgts)
       call lininterp_finish(lon_wgts)
    end do
    deallocate( soil_erodibility_in, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'soil_erod_init: failed to deallocate soil_erodibility_in, ierr = ',ierr
       call endrun('soil_erod_init: failed to deallocate soil_erodibility_in')
    end if
#if ( defined DUSTTYPES )
    deallocate( illite_in_clay, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate illite_in_clay, ierr = ',ierr
       call endrun
    end if

    deallocate( kaolinite_in_clay, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate kaolinite_in_clay, ierr = ',ierr
       call endrun
    end if

    deallocate( kaolinite_in_silt, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate kaolinite_in_silt, ierr = ',ierr
       call endrun
    end if

    deallocate( mont_in_clay, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate mont_in_clay, ierr = ',ierr
       call endrun
    end if

    deallocate( hematite_in_clay, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate hematite_in_clay, ierr = ',ierr
       call endrun
    end if

    deallocate( hematite_in_silt, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate hematite_in_silt, ierr = ',ierr
       call endrun
    end if

    deallocate( quartz_in_clay, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate quartz_in_clay, ierr = ',ierr
       call endrun
    end if
    deallocate( quartz_in_silt, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate quartz_in_silt, ierr = ',ierr
       call endrun
    end if

    deallocate( calcite_in_clay, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate calcite_in_clay, ierr = ',ierr
       call endrun
    end if

    deallocate( calcite_in_silt, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate calcite_in_silt, ierr = ',ierr
       call endrun
    end if

    deallocate( feldspar_in_clay, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate feldspar_in_clay, ierr = ',ierr
       call endrun
    end if

    deallocate( feldspar_in_silt, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate feldspar_in_silt, ierr = ',ierr
       call endrun
    end if

    deallocate( gypsum_in_silt, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate gypsum_in_silt, ierr = ',ierr
       call endrun
    end if
#if ( defined IRONPROC )
    deallocate( comb_in_accum_sol, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'combfe_initialize: failed to deallocate comb_in_accum_sol, ierr = ',ierr
       call endrun
    end if
    deallocate( comb_in_aitken_sol, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'combfe_initialize: failed to deallocate comb_in_aitken_sol, ierr = ',ierr
       call endrun
    end if
    deallocate( comb_in_coarse_sol, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'combfe_initialize: failed to deallocate comb_in_coarse_sol, ierr = ',ierr
       call endrun
    end if
    deallocate( comb_in_accum_insol, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'combfe_initialize: failed to deallocate comb_in_accum_insol, ierr = ',ierr
       call endrun
    end if
    deallocate( comb_in_aitken_insol, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'combfe_initialize: failed to deallocate comb_in_aitken_insol, ierr = ',ierr
       call endrun
    end if
    deallocate( comb_in_coarse_insol, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'combfe_initialize: failed to deallocate comb_in_coarse_insol, ierr = ',ierr
       call endrun
    end if 
   deallocate( comb_lats )
    deallocate( comb_lons )
#endif
#endif

    deallocate( dst_lats )
    deallocate( dst_lons )

  end  subroutine soil_erod_init

end module soil_erod_mod
