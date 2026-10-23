!===============================================================================
! Dust for Modal Aerosol Model
!===============================================================================
module dust_model 
  use shr_kind_mod,     only: r8 => shr_kind_r8, cl => shr_kind_cl
  use spmd_utils,       only: masterproc
  use cam_abortutils,   only: endrun
  use modal_aero_data,  only: ntot_amode, ndst=>nDust

  implicit none
  private

  public :: dust_names
  public :: dust_nbin
  public :: dust_nnum
#if ( defined DUSTTYPES )
  public :: mine_nbin 
#if ( defined IRONPROC )
  public :: iron_nbin 
#endif
#endif
  public :: dust_indices
  public :: dust_emis
  public :: dust_readnl
  public :: dust_init
  public :: dust_active

  integer, protected :: dust_nbin != 2
  integer, protected :: dust_nnum != 2
#if ( defined DUSTTYPES )
  integer, protected :: mine_nbin
#if ( defined IRONPROC )
  integer, protected :: iron_nbin
#endif
#endif

#if  ( defined MODAL_AERO_3MODE )
  character(len=6), protected, allocatable :: dust_names(:)
#elif ( defined MODAL_AERO_4MODE )
#if ( defined DUSTTYPES )
  character(len=7), protected, allocatable :: dust_names(:)
#else
  character(len=6), protected, allocatable :: dust_names(:)
#endif
#elif ( defined MODAL_AERO_7MODE )
  character(len=6), protected, allocatable :: dust_names(:)
#endif
  real(r8), allocatable :: dust_dmt_grd(:)
  real(r8), allocatable :: dust_emis_sclfctr(:)
#if ( defined DUSTTYPES )
  real(r8), allocatable :: brit_emis_sclfctr(:)
#if ( defined IRONPROC )
  real(r8), allocatable :: feRR_emis_sclfctr(:)
  real(r8), allocatable :: feMe_emis_sclfctr(:)
  real(r8), allocatable :: feSl_emis_sclfctr(:)
#endif
#endif

  integer , protected, allocatable :: dust_indices(:)
  real(r8), allocatable :: dust_dmt_vwr(:)
  real(r8), allocatable :: dust_stk_crc(:)

  real(r8)          :: dust_emis_fact = -1.e36_r8        ! tuning parameter for dust emissions
  character(len=cl) :: soil_erod_file = 'soil_erod_file' ! full pathname for soil erodibility dataset
  character(len=cl) :: comb_emis_file = 'comb_emis_file' ! full pathname for combustion dataset 
  character(len=cl) :: comb_file_list = 'comb_file_list'

  logical :: dust_active = .false.

 contains

  !=============================================================================
  ! reads dust namelist options
  !=============================================================================
  subroutine dust_readnl(nlfile)

    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use mpishorthand
    use cam_logfile,       only: iulog

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'dust_readnl'

    namelist /dust_nl/ dust_emis_fact, soil_erod_file, comb_emis_file, comb_file_list

    !-----------------------------------------------------------------------------
    ! Read namelist
    if (masterproc) then
       unitn = getunit()
write(iulog,*)'longlei, nlfile', trim(nlfile)
!write( iulog, '(1x,a)' ) &
!              nlfile
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'dust_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, dust_nl, iostat=ierr)
write(iulog,*)'longlei, dust_emis_fact',dust_emis_fact
write(iulog,*)'longlei, comb_emis_file',comb_emis_file
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

#ifdef SPMD
    ! Broadcast namelist variables
    call mpibcast(dust_emis_fact, 1,                   mpir8,   0, mpicom)
    call mpibcast(soil_erod_file, len(soil_erod_file), mpichar, 0, mpicom)
    call mpibcast(comb_emis_file, len(comb_emis_file), mpichar, 0, mpicom)
    call mpibcast(comb_file_list, len(comb_file_list), mpichar, 0, mpicom)
#endif

  end subroutine dust_readnl

  !=============================================================================
  !=============================================================================
  subroutine dust_init()
    use soil_erod_mod, only: soil_erod_init
    use constituents,  only: cnst_get_ind
    use rad_constituents, only: rad_cnst_get_info
    use dust_common,   only: dust_set_params
    use cam_logfile,       only: iulog

    integer :: l, m, mm, ndx, nspec, n
    character(len=32) :: spec_name
    integer, parameter :: mymodes(7) = (/ 2, 1, 3, 4, 5, 6, 7 /) ! tricky order ...

    !dust_nbin = ndst ! not this number any more -Longlei
    !dust_nnum = ndst

#if  ( defined MODAL_AERO_3MODE )
    dust_nbin = 2
    dust_nnum = 2
#elif ( defined MODAL_AERO_4MODE )
#if ( defined DUSTTYPES )
#if ( defined IRONPROC )
    dust_nbin = 48 !42
    dust_nnum = 3
    mine_nbin = 33
    iron_nbin = 15
#else
    dust_nbin = 24
    dust_nnum = 3
    mine_nbin = 39
#endif
#else
    dust_nbin = 3
    dust_nnum = 3
#endif
#elif ( defined MODAL_AERO_7MODE )
    dust_nbin = 2
    dust_nnum = 2
#endif

    !allocate( dust_names(2*ndst) )  ! demension changed -Longlei
    !allocate( dust_indices(2*ndst) )
    !allocate( dust_dmt_grd(ndst+1) )
    !allocate( dust_emis_sclfctr(ndst) )
    !allocate( dust_dmt_vwr(ndst) )
    !allocate( dust_stk_crc(ndst) )

    allocate( dust_names(dust_nbin+dust_nnum) )
    allocate( dust_indices(dust_nbin+dust_nnum) )
    allocate( dust_dmt_grd(dust_nnum+1) )
    allocate( dust_emis_sclfctr(dust_nbin) ) 
#if ( defined DUSTTYPES )
    allocate( dust_dmt_vwr(dust_nnum) )
    allocate( dust_stk_crc(dust_nnum) )
    allocate( brit_emis_sclfctr(mine_nbin) )
#if ( defined IRONPROC )
    allocate( feRR_emis_sclfctr(iron_nbin) )
    allocate( feMe_emis_sclfctr(iron_nbin) )
    allocate( feSl_emis_sclfctr(iron_nbin) )
#endif
#else
    allocate( dust_dmt_vwr(dust_nbin) )
    allocate( dust_stk_crc(dust_nbin) )
#endif

    if ( ntot_amode == 3 ) then
       dust_dmt_grd(:) = (/ 0.1e-6_r8, 1.0e-6_r8, 10.0e-6_r8/)
       dust_emis_sclfctr(:) = (/ 0.011_r8,0.989_r8 /)
    elseif ( ntot_amode == 4 ) then
       dust_dmt_grd(:) = (/ 0.01e-6_r8, 0.1e-6_r8, 1.0e-6_r8, 10.0e-6_r8 /) ! Aitken dust
#if ( defined DUSTTYPES )
! based on Scanza et al 2015 and Hamilton et al 2019 we remove hematite from silt
! Journet maps adds chlorite, mica, and vermiculite and named it as kaolinite in the silt 
! Journet maps adds Feldspar in both clay and silt fraction

     brit_emis_sclfctr(:) = &
   (/ 1.000_r8, 1.000_r8, 1.000_r8, 0.000_r8, 1.000_r8, 0.000_r8, 1.000_r8, 0.000_r8, 1.000_r8, 0.000_r8, 1.000_r8, 0.000_r8, 0.000_r8,& ! aitken
      1.000_r8, 1.000_r8, 1.000_r8, 0.000_r8, 1.000_r8, 0.000_r8, 1.000_r8, 0.000_r8, 1.000_r8, 0.000_r8, 1.000_r8, 0.000_r8, 0.000_r8,& ! accum
      0.695_r8, 0.695_r8, 0.695_r8, 0.305_r8, 0.695_r8, 0.305_r8, 0.695_r8, 0.000_r8, 0.695_r8, 0.305_r8, 0.695_r8, 0.305_r8, 0.305_r8/) ! coarse
      ! ^illi^   ^mont^   ^kaol_c^  ^kaol_s^   ^fel_c^  ^fel_s^   ^hem_c^   ^hem_s    ^qua_c    ^qua_s    ^cal_c    ^cal_s   ^ gypsu

#if ( defined IRONPROC )
       dust_emis_sclfctr(:) = &
     (/ 0.001_r8, 0.001_r8, 0.001_r8, 0.001_r8, 0.001_r8, 0.001_r8, 0.001_r8, 0.001_r8, 0.001_r8, 0.001_r8, 0.001_r8, 0.001_r8, 0.001_r8, 0.001_r8, 0.001_r8, 0.001_r8, &
        0.011_r8, 0.011_r8, 0.011_r8, 0.011_r8, 0.011_r8, 0.011_r8, 0.011_r8, 0.011_r8, 0.011_r8, 0.011_r8, 0.011_r8, 0.011_r8, 0.011_r8, 0.011_r8, 0.011_r8, 0.011_r8, &
        0.989_r8, 0.989_r8, 0.989_r8, 0.989_r8, 0.989_r8, 0.989_r8, 0.989_r8, 0.989_r8, 0.989_r8, 0.989_r8, 0.989_r8, 0.989_r8, 0.989_r8, 0.989_r8, 0.989_r8, 0.989_r8 /)

       dust_names(:) = (/ 'dst1_a2', 'dst2_a2', 'dst3_a2', 'dst4_a2', 'dst5_a2', 'dst6_a2', 'dst7_a2', 'dst8_a2', &
                          'fets_a2', 'fetm_a2', 'fess_a2', 'fesm_a2', 'fetc_a2', 'fesc_a2', 'fetf_a2', 'fesf_a2', &
                          'dst1_a1', 'dst2_a1', 'dst3_a1', 'dst4_a1', 'dst5_a1', 'dst6_a1', 'dst7_a1', 'dst8_a1', &
                          'fets_a1', 'fetm_a1', 'fess_a1', 'fesm_a1', 'fetc_a1', 'fesc_a1', 'fetf_a1', 'fesf_a1', &
                          'dst1_a3', 'dst2_a3', 'dst3_a3', 'dst4_a3', 'dst5_a3', 'dst6_a3', 'dst7_a3', 'dst8_a3', &
                          'fets_a3', 'fetm_a3', 'fess_a3', 'fesm_a3', 'fetc_a3', 'fesc_a3', 'fetf_a3', 'fesf_a3', &
                          'num_a2 ', 'num_a1 ', 'num_a3 ' /)

! currently assume that aitken = accumulation. Sum feRR+feMed+feSlow = 1.0
! RR now equals 10% (down from 25%) as per Conway et al. 2018 

       feRR_emis_sclfctr(:) = &
                       (/ 0.003_r8, 0.005_r8, 0.004_r8, 0.003_r8, 0.000_r8, &   ! aitken
                          0.003_r8, 0.005_r8, 0.004_r8, 0.003_r8, 0.000_r8, &   ! accum
                          0.003_r8, 0.005_r8, 0.004_r8, 0.003_r8, 0.000_r8 /)   ! coarse
                          !^illite^   ^mont^    ^kaol^   ^felds^   ^hemat^

       feMe_emis_sclfctr(:) = &
                       (/ 0.997_r8, 0.995_r8, 0.039_r8, 0.027_r8, 0.000_r8, &   ! aitken
                          0.997_r8, 0.995_r8, 0.039_r8, 0.027_r8, 0.000_r8, &   ! accum
                          0.997_r8, 0.995_r8, 0.039_r8, 0.027_r8, 0.000_r8 /)   ! coarse
                          !^illite^   ^mont^    ^kaol^   ^felds^   ^hemat^

       feSl_emis_sclfctr(:) = &
                       (/ 0.000_r8, 0.000_r8, 0.957_r8, 0.970_r8, 1.000_r8, &   ! aitken
                          0.000_r8, 0.000_r8, 0.957_r8, 0.970_r8, 1.000_r8, &   ! accum
                          0.000_r8, 0.000_r8, 0.957_r8, 0.970_r8, 1.000_r8 /)   ! coarse
                          !^illite^   ^mont^    ^kaol^   ^felds^   ^hemat^
#else
       dust_emis_sclfctr(:) = &
                              (/ 0.001_r8, 0.001_r8, 0.001_r8, 0.001_r8, 0.001_r8, 0.001_r8, 0.001_r8, 0.001_r8, &
                                 0.011_r8, 0.011_r8, 0.011_r8, 0.011_r8, 0.011_r8, 0.011_r8, 0.011_r8, 0.011_r8, &
                                 0.989_r8, 0.989_r8, 0.989_r8, 0.989_r8, 0.989_r8, 0.989_r8, 0.989_r8, 0.989_r8 /)

       dust_names(:) = (/ 'dst1_a2', 'dst2_a2', 'dst3_a2', 'dst4_a2', 'dst5_a2', 'dst6_a2', 'dst7_a2', 'dst8_a2', &
                          'dst1_a1', 'dst2_a1', 'dst3_a1', 'dst4_a1', 'dst5_a1', 'dst6_a1', 'dst7_a1', 'dst8_a1', &
                          'dst1_a3', 'dst2_a3', 'dst3_a3', 'dst4_a3', 'dst5_a3', 'dst6_a3', 'dst7_a3', 'dst8_a3', &
                          'num_a2 ', 'num_a1 ', 'num_a3 ' /)
#endif
#else
       dust_emis_sclfctr(:) = (/ 1.65E-05_r8, 0.011_r8, 0.989_r8 /) ! Aitken dust
#endif
    else if( ntot_amode == 7 ) then
       dust_dmt_grd(:) = (/ 0.1e-6_r8, 2.0e-6_r8, 10.0e-6_r8/)
       dust_emis_sclfctr(:) = (/ 0.13_r8, 0.87_r8 /)
    endif

  !  ndx = 0
  !  do mm = 1, ntot_amode
  !     m = mymodes(mm)
  !     call rad_cnst_get_info(0, m, nspec=nspec)
  !     do l = 1, nspec
  !        call rad_cnst_get_info(0, m, l, spec_name=spec_name )
  !        if (spec_name(:3) == 'dst') then
  !           ndx=ndx+1
  !           dust_names(ndx) = spec_name
  !           dust_names(ndst+ndx) = 'num_'//spec_name(5:)
  !           call cnst_get_ind(dust_names(     ndx), dust_indices(     ndx))
  !           call cnst_get_ind(dust_names(ndst+ndx), dust_indices(ndst+ndx))
  !        endif
  !     enddo
  !  enddo

    do n = 1, dust_nbin
               !write( iulog, '(1x,a)' ) &
               !      'dust_names=',dust_names(n)
       call cnst_get_ind(dust_names(n), dust_indices(n))
    end do
    do n = 1, dust_nnum
       call cnst_get_ind(dust_names(dust_nbin+n), dust_indices(dust_nbin+n))
    enddo

    dust_active = any(dust_indices(:) > 0)
    if (.not.dust_active) return
   
    call  soil_erod_init( dust_emis_fact, soil_erod_file, comb_emis_file, comb_file_list )
#if ( defined DUSTTYPES )
    call dust_set_params( dust_nnum, dust_dmt_grd, dust_dmt_vwr, dust_stk_crc )
#else
    call dust_set_params( dust_nbin, dust_dmt_grd, dust_dmt_vwr, dust_stk_crc )
#endif

    !write(iulog,*) 'marks 01 in dust_init of dust_model'
  end subroutine dust_init

  !===============================================================================
  !===============================================================================
  subroutine dust_emis( ncol, lchnk, dust_flux_in, cflx, soil_erod )
    use soil_erod_mod, only : soil_erod_fact
    use soil_erod_mod, only : soil_erodibility
#if ( defined DUSTTYPES )
    ! apply brittle theory online -Longlei
    use soil_erod_mod, only :  illite_clay,    kaolinite_clay,   kaolinite_silt, mont_clay
    use soil_erod_mod, only :  hematite_clay,  hematite_silt,    quartz_clay,   quartz_silt
    use soil_erod_mod, only :  calcite_clay,   calcite_silt,     feldspar_clay, feldspar_silt, gypsum_silt
#endif
#if ( defined IRONPROC )
    use soil_erod_mod, only :  comb_accum_sol, comb_aitken_sol, comb_coarse_sol
    use soil_erod_mod, only :  comb_accum_insol, comb_aitken_insol, comb_coarse_insol
#endif
    use mo_constants,  only : dust_density
    use physconst,     only : pi
    use cam_history,   only : outfld
    use ppgrid,        only : pcols
    use cam_logfile,       only: iulog

  ! args
    integer,  intent(in)    :: ncol, lchnk
    real(r8), intent(in)    :: dust_flux_in(:,:)
    real(r8), intent(inout) :: cflx(:,:)
    real(r8), intent(out)   :: soil_erod(:)

  ! local vars
    integer :: i, m, idst, inum
    real(r8) :: x_mton
    real(r8),parameter :: soil_erod_threshold = 0.1_r8
#if ( defined DUSTTYPES )
    real(r8) ::  frac_I_accum(pcols), frac_I_aitken(pcols), frac_I_coarse(pcols)
    real(r8) ::  frac_K_accum(pcols), frac_K_aitken(pcols), frac_K_coarse(pcols)
    real(r8) ::  frac_M_accum(pcols), frac_M_aitken(pcols), frac_M_coarse(pcols)
    real(r8) ::  frac_H_accum(pcols), frac_H_aitken(pcols), frac_H_coarse(pcols)
    real(r8) ::  frac_Q_accum(pcols), frac_Q_aitken(pcols), frac_Q_coarse(pcols)
    real(r8) ::  frac_C_accum(pcols), frac_C_aitken(pcols), frac_C_coarse(pcols)
    real(r8) ::  frac_F_accum(pcols), frac_F_aitken(pcols), frac_F_coarse(pcols)
    real(r8) ::  frac_G_accum(pcols), frac_G_aitken(pcols), frac_G_coarse(pcols)
    real(r8) ::  dust_flux_sum(pcols)
#endif

#if ( defined IRONPROC )
    real(r8) ::  frac_Fe_accum(pcols),   frac_Fe_aitken(pcols),   frac_Fe_coarse(pcols)
    real(r8) ::  frac_Fets_accum(pcols), frac_Fets_aitken(pcols), frac_Fets_coarse(pcols)
    real(r8) ::  frac_Fetm_accum(pcols), frac_Fetm_aitken(pcols), frac_Fetm_coarse(pcols)
    real(r8) ::  frac_Fess_accum(pcols), frac_Fess_aitken(pcols), frac_Fess_coarse(pcols)
    real(r8) ::  frac_Fesm_accum(pcols), frac_Fesm_aitken(pcols), frac_Fesm_coarse(pcols)
    real(r8) ::  frac_Fetc_accum(pcols), frac_Fetc_aitken(pcols), frac_Fetc_coarse(pcols)
    real(r8) ::  frac_Fesc_accum(pcols), frac_Fesc_aitken(pcols), frac_Fesc_coarse(pcols)
#endif

#if ( defined DUSTTYPES )
!brittle frag theory
    do i = 1, ncol
       frac_I_aitken(i) = illite_clay(i,lchnk)    * brit_emis_sclfctr(1)  
       frac_I_accum(i)  = illite_clay(i,lchnk)    * brit_emis_sclfctr(14) 
       frac_I_coarse(i) = illite_clay(i,lchnk)    * brit_emis_sclfctr(27) 

       frac_M_aitken(i) = mont_clay(i,lchnk)      * brit_emis_sclfctr(2)  
       frac_M_accum(i)  = mont_clay(i,lchnk)      * brit_emis_sclfctr(15) 
       frac_M_coarse(i) = mont_clay(i,lchnk)      * brit_emis_sclfctr(28) 

       frac_K_aitken(i) = kaolinite_clay(i,lchnk) * brit_emis_sclfctr(3)  + kaolinite_silt(i,lchnk) * brit_emis_sclfctr(4)   
       frac_K_accum(i)  = kaolinite_clay(i,lchnk) * brit_emis_sclfctr(16) + kaolinite_silt(i,lchnk) * brit_emis_sclfctr(17)
       frac_K_coarse(i) = kaolinite_clay(i,lchnk) * brit_emis_sclfctr(29) + kaolinite_silt(i,lchnk) * brit_emis_sclfctr(30)

       frac_F_aitken(i) = feldspar_clay(i,lchnk)  * brit_emis_sclfctr(5)  + feldspar_silt(i,lchnk)  * brit_emis_sclfctr(6)  
       frac_F_accum(i)  = feldspar_clay(i,lchnk)  * brit_emis_sclfctr(18) + feldspar_silt(i,lchnk)  * brit_emis_sclfctr(19) 
       frac_F_coarse(i) = feldspar_clay(i,lchnk)  * brit_emis_sclfctr(31) + feldspar_silt(i,lchnk)  * brit_emis_sclfctr(32) 

       frac_H_aitken(i) = hematite_clay(i,lchnk)  * brit_emis_sclfctr(7)  + hematite_silt(i,lchnk)  * brit_emis_sclfctr(8)
       frac_H_accum(i)  = hematite_clay(i,lchnk)  * brit_emis_sclfctr(20) + hematite_silt(i,lchnk)  * brit_emis_sclfctr(21)
       frac_H_coarse(i) = hematite_clay(i,lchnk)  * brit_emis_sclfctr(33) + hematite_silt(i,lchnk)  * brit_emis_sclfctr(34)

       frac_Q_aitken(i) = quartz_clay(i,lchnk)    * brit_emis_sclfctr(9)  + quartz_silt(i,lchnk)    * brit_emis_sclfctr(10)
       frac_Q_accum(i)  = quartz_clay(i,lchnk)    * brit_emis_sclfctr(22) + quartz_silt(i,lchnk)    * brit_emis_sclfctr(23)
       frac_Q_coarse(i) = quartz_clay(i,lchnk)    * brit_emis_sclfctr(35) + quartz_silt(i,lchnk)    * brit_emis_sclfctr(36)

       frac_C_aitken(i) = calcite_clay(i,lchnk)   * brit_emis_sclfctr(11) + calcite_silt(i,lchnk)   * brit_emis_sclfctr(12)
       frac_C_accum(i)  = calcite_clay(i,lchnk)   * brit_emis_sclfctr(24) + calcite_silt(i,lchnk)   * brit_emis_sclfctr(25)
       frac_C_coarse(i) = calcite_clay(i,lchnk)   * brit_emis_sclfctr(37) + calcite_silt(i,lchnk)   * brit_emis_sclfctr(38)

       frac_G_aitken(i) = gypsum_silt(i,lchnk)    * brit_emis_sclfctr(13)
       frac_G_accum(i)  = gypsum_silt(i,lchnk)    * brit_emis_sclfctr(26)
       frac_G_coarse(i) = gypsum_silt(i,lchnk)    * brit_emis_sclfctr(39)

#if ( defined IRONPROC )
! iron dust mass    
       frac_Fets_aitken(i)=(frac_I_aitken(i) *0.0400_r8*feSl_emis_sclfctr(1)) + &
                           (frac_M_aitken(i) *0.1100_r8*feSl_emis_sclfctr(2)) + &
                           (frac_K_aitken(i) *0.0024_r8*feSl_emis_sclfctr(3)) + &
                           (frac_F_aitken(i) *0.0034_r8*feSl_emis_sclfctr(4)) + &
                           (frac_H_aitken(i) *0.5750_r8*feSl_emis_sclfctr(5))
       frac_Fets_accum(i) =(frac_I_accum (i) *0.0400_r8*feSl_emis_sclfctr(6)) + &
                           (frac_M_accum (i) *0.1100_r8*feSl_emis_sclfctr(7)) + &
                           (frac_K_accum (i) *0.0024_r8*feSl_emis_sclfctr(8)) + &
                           (frac_F_accum (i) *0.0034_r8*feSl_emis_sclfctr(9)) + &
                           (frac_H_accum (i) *0.5750_r8*feSl_emis_sclfctr(10))
       frac_Fets_coarse(i)=(frac_I_coarse(i) *0.0400_r8*feSl_emis_sclfctr(11))+ &
                           (frac_M_coarse(i) *0.1100_r8*feSl_emis_sclfctr(12))+ &
                           (frac_K_coarse(i) *0.0024_r8*feSl_emis_sclfctr(13))+ &
                           (frac_F_coarse(i) *0.0034_r8*feSl_emis_sclfctr(14))+ &
                           (frac_H_coarse(i) *0.5750_r8*feSl_emis_sclfctr(15))

       frac_Fetm_aitken(i)=(frac_I_aitken(i) *0.0400_r8*feMe_emis_sclfctr(1)) + &
                           (frac_M_aitken(i) *0.1100_r8*feMe_emis_sclfctr(2)) + &
                           (frac_K_aitken(i) *0.0024_r8*feMe_emis_sclfctr(3)) + &
                           (frac_F_aitken(i) *0.0034_r8*feMe_emis_sclfctr(4)) + &
                           (frac_H_aitken(i) *0.5750_r8*feMe_emis_sclfctr(5))
       frac_Fetm_accum(i) =(frac_I_accum (i) *0.0400_r8*feMe_emis_sclfctr(6)) + &
                           (frac_M_accum (i) *0.1100_r8*feMe_emis_sclfctr(7)) + &
                           (frac_K_accum (i) *0.0024_r8*feMe_emis_sclfctr(8)) + &
                           (frac_F_accum (i) *0.0034_r8*feMe_emis_sclfctr(9)) + &
                           (frac_H_accum (i) *0.5750_r8*feMe_emis_sclfctr(10))
       frac_Fetm_coarse(i)=(frac_I_coarse(i) *0.0400_r8*feMe_emis_sclfctr(11))+ &
                           (frac_M_coarse(i) *0.1100_r8*feMe_emis_sclfctr(12))+ &
                           (frac_K_coarse(i) *0.0024_r8*feMe_emis_sclfctr(13))+ &
                           (frac_F_coarse(i) *0.0034_r8*feMe_emis_sclfctr(14))+ &
                           (frac_H_coarse(i) *0.5750_r8*feMe_emis_sclfctr(15))

       frac_Fess_accum(i) =0.0_r8
       frac_Fess_aitken(i)=0.0_r8
       frac_Fess_coarse(i)=0.0_r8

       frac_Fesm_aitken(i)=(frac_I_aitken(i) *0.0400_r8*feRR_emis_sclfctr(1)) + &
                           (frac_M_aitken(i) *0.1100_r8*feRR_emis_sclfctr(2)) + &
                           (frac_K_aitken(i) *0.0024_r8*feRR_emis_sclfctr(3)) + &
                           (frac_F_aitken(i) *0.0034_r8*feRR_emis_sclfctr(4)) + &
                           (frac_H_aitken(i) *0.5750_r8*feRR_emis_sclfctr(5))
       frac_Fesm_accum(i) =(frac_I_accum (i) *0.0400_r8*feRR_emis_sclfctr(6)) + &
                           (frac_M_accum (i) *0.1100_r8*feRR_emis_sclfctr(7)) + &
                           (frac_K_accum (i) *0.0024_r8*feRR_emis_sclfctr(8)) + &
                           (frac_F_accum (i) *0.0034_r8*feRR_emis_sclfctr(9)) + &
                           (frac_H_accum (i) *0.5750_r8*feRR_emis_sclfctr(10))
       frac_Fesm_coarse(i)=(frac_I_coarse(i) *0.0400_r8*feRR_emis_sclfctr(11))+ &
                           (frac_M_coarse(i) *0.1100_r8*feRR_emis_sclfctr(12))+ &
                           (frac_K_coarse(i) *0.0024_r8*feRR_emis_sclfctr(13))+ &
                           (frac_F_coarse(i) *0.0034_r8*feRR_emis_sclfctr(14))+ &
                           (frac_H_coarse(i) *0.5750_r8*feRR_emis_sclfctr(15))

! to conserve dust mass must remove the iron mass tracer content
       frac_I_aitken(i) = frac_I_aitken(i) *(1.0_r8-0.0400_r8)
       frac_I_accum(i)  = frac_I_accum (i) *(1.0_r8-0.0400_r8)
       frac_I_coarse(i) = frac_I_coarse(i) *(1.0_r8-0.0400_r8)
  
       frac_M_aitken(i) = frac_M_aitken(i) *(1.0_r8-0.1100_r8)
       frac_M_accum(i)  = frac_M_accum (i) *(1.0_r8-0.1100_r8)
       frac_M_coarse(i) = frac_M_coarse(i) *(1.0_r8-0.1100_r8)

       frac_K_aitken(i) = frac_K_aitken(i) *(1.0_r8-0.0024_r8)
       frac_K_accum(i)  = frac_K_accum (i) *(1.0_r8-0.0024_r8)
       frac_K_coarse(i) = frac_K_coarse(i) *(1.0_r8-0.0024_r8)

       frac_F_aitken(i) = frac_F_aitken(i) *(1.0_r8-0.0034_r8)
       frac_F_accum(i)  = frac_F_accum (i) *(1.0_r8-0.0034_r8)
       frac_F_coarse(i) = frac_F_coarse(i) *(1.0_r8-0.0034_r8)

       frac_H_aitken(i) = frac_H_aitken(i) *(1.0_r8-0.5750_r8)
       frac_H_accum(i)  = frac_H_accum (i) *(1.0_r8-0.5750_r8)
       frac_H_coarse(i) = frac_H_coarse(i) *(1.0_r8-0.5750_r8)


!combustion iron
     !frac_Fetc_accum(i)  = comb_accum_insol(i,lchnk) 
     !frac_Fetc_aitken(i) = comb_aitken_insol(i,lchnk)
     !frac_Fetc_coarse(i) = comb_coarse_insol(i,lchnk)

     !frac_Fesc_accum(i)  = comb_accum_sol(i,lchnk) 
     !frac_Fesc_aitken(i) = comb_aitken_sol(i,lchnk)
     !frac_Fesc_coarse(i) = comb_coarse_sol(i,lchnk)
 
!updated to read in from file - dsh
     frac_Fetc_accum(i)  = 0.0_r8
     frac_Fetc_aitken(i) = 0.0_r8
     frac_Fetc_coarse(i) = 0.0_r8

     frac_Fesc_accum(i)  = 0.0_r8
     frac_Fesc_aitken(i) = 0.0_r8
     frac_Fesc_coarse(i) = 0.0_r8
#endif
    end do
    call outfld('Illite_accum',frac_I_accum(:),pcols,lchnk)
    call outfld('Illite_aitken',frac_I_aitken(:),pcols,lchnk)
    call outfld('Illite_coarse',frac_I_coarse(:),pcols,lchnk)

    call outfld('Kaolinite_accum',frac_K_accum(:),pcols,lchnk)
    call outfld('Kaolinite_aitken',frac_K_aitken(:),pcols,lchnk)
    call outfld('Kaolinite_coarse',frac_K_coarse(:),pcols,lchnk)

    call outfld('Mont_accum',frac_M_accum(:),pcols,lchnk)
    call outfld('Mont_aitken',frac_M_aitken(:),pcols,lchnk)
    call outfld('Mont_coarse',frac_M_coarse(:),pcols,lchnk)

    call outfld('Hematite_accum',frac_H_accum(:),pcols,lchnk)
    call outfld('Hematite_aitken',frac_H_aitken(:),pcols,lchnk)
    call outfld('Hematite_coarse',frac_H_coarse(:),pcols,lchnk)

    call outfld('Quartz_accum',frac_Q_accum(:),pcols,lchnk)
    call outfld('Quartz_aitken',frac_Q_aitken(:),pcols,lchnk)
    call outfld('Quartz_coarse',frac_Q_coarse(:),pcols,lchnk)

    call outfld('Calcite_accum',frac_C_accum(:),pcols,lchnk)
    call outfld('Calcite_aitken',frac_C_aitken(:),pcols,lchnk)
    call outfld('Calcite_coarse',frac_C_coarse(:),pcols,lchnk)

    call outfld('Feldspar_accum',frac_F_accum(:),pcols,lchnk)
    call outfld('Feldspar_aitken',frac_F_aitken(:),pcols,lchnk)
    call outfld('Feldspar_coarse',frac_F_coarse(:),pcols,lchnk)

    call outfld('Gypsum_accum',frac_G_accum(:),pcols,lchnk)
    call outfld('Gypsum_aitken',frac_G_aitken(:),pcols,lchnk)
    call outfld('Gypsum_coarse',frac_G_coarse(:),pcols,lchnk)

#if ( defined IRONPROC )
    call outfld('Iron_accum',frac_Fe_accum(:),pcols,lchnk)
    call outfld('Iron_aitken',frac_Fe_aitken(:),pcols,lchnk)
    call outfld('Iron_coarse',frac_Fe_coarse(:),pcols,lchnk)

    call outfld('fetcSF_a1' ,frac_Fetc_accum(:),pcols,lchnk)
    call outfld('fetcSF_a2',frac_Fetc_aitken(:),pcols,lchnk)
    call outfld('fetcSF_a3',frac_Fetc_coarse(:),pcols,lchnk)
    call outfld('fescSF_a1' ,frac_Fesc_accum(:),pcols,lchnk)
    call outfld('fescSF_a2',frac_Fesc_aitken(:),pcols,lchnk)
    call outfld('fescSF_a3',frac_Fesc_coarse(:),pcols,lchnk)
#endif

#endif

    ! set dust emissions
#if ( defined DUSTTYPES )
#if ( defined IRONPROC )
    col_loop: do i =1,ncol
       dust_flux_sum(i) = sum( -dust_flux_in(i,:) )

       !soil_erod(i) = soil_erodibility( i, lchnk )
       soil_erod(i) = 1._r8 ! Kok scheme

       if( soil_erod(i) .lt. soil_erod_threshold ) soil_erod(i) = 0._r8

       ! rebin and adjust dust emissons..
       ! Aitken mode minerals
       ! moving tracers 1 and 2 into 3
       idst = dust_indices(1)
       cflx(i,idst) = 0.0_r8 !sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(1)*frac_I_aitken(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(2)
       cflx(i,idst) = 0.0_r8 !sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(2)*frac_K_aitken(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(3)
       cflx(i,idst) = 0.0_r8!sum( -dust_flux_in(i,:) )*(dust_emis_sclfctr(1)+dust_emis_sclfctr(2)+dust_emis_sclfctr(3))*(frac_I_aitken(i)+frac_K_aitken(i)+frac_M_aitken(i))*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(4)
       cflx(i,idst) = 0.0_r8 !sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(4)*frac_H_aitken(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(5)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(4)*frac_H_aitken(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(6)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(6)*frac_C_aitken(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(7)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(7)*frac_F_aitken(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(8)
       !cflx(i,idst) = sum( -dust_flux_in(i,:) )*(dust_emis_sclfctr(1)*frac_I_aitken(i)+dust_emis_sclfctr(2)*frac_K_aitken(i)+dust_emis_sclfctr(3)*frac_M_aitken(i)+dust_emis_sclfctr(5)*frac_Q_aitken(i)+dust_emis_sclfctr(8)*frac_G_aitken(i))*soil_erod(i)/soil_erod_fact*1.15_r8
       cflx(i,idst) = ((sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(1)*frac_I_aitken(i)*soil_erod(i)/soil_erod_fact*1.15_r8) + &
                       (sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(2)*frac_K_aitken(i)*soil_erod(i)/soil_erod_fact*1.15_r8) + &
                       (sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(3)*frac_M_aitken(i)*soil_erod(i)/soil_erod_fact*1.15_r8) + &
                       (sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(5)*frac_Q_aitken(i)*soil_erod(i)/soil_erod_fact*1.15_r8) + &
                       (sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(8)*frac_G_aitken(i)*soil_erod(i)/soil_erod_fact*1.15_r8))
       idst = dust_indices(9)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(9)*frac_Fets_aitken(i)*soil_erod(i)/soil_erod_fact*1.15_r8 
       idst = dust_indices(10)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(10)*frac_Fetm_aitken(i)*soil_erod(i)/soil_erod_fact*1.15_r8 
       idst = dust_indices(11)
       cflx(i,idst) = 0.0_r8
       idst = dust_indices(12)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(12)*frac_Fesm_aitken(i)*soil_erod(i)/soil_erod_fact*1.15_r8 
       idst = dust_indices(13)
       cflx(i,idst) = 0.0_r8 
       idst = dust_indices(14)
       cflx(i,idst) = 0.0_r8 
       idst = dust_indices(15)
       cflx(i,idst) = 0.0_r8 
       idst = dust_indices(16)
       cflx(i,idst) = 0.0_r8

       ! Accumulation mode minerals
       idst = dust_indices(17)
       cflx(i,idst) = 0.0_r8 !sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(17)*frac_I_accum(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(18)
       cflx(i,idst) = 0.0_r8 !sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(18)*frac_K_accum(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(19)
       cflx(i,idst) = 0.0_r8 !sum( -dust_flux_in(i,:) )*(dust_emis_sclfctr(17)+dust_emis_sclfctr(18)+dust_emis_sclfctr(19))*(frac_I_accum(i)+frac_K_accum(i)+frac_M_accum(i))*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(20)
       cflx(i,idst) = 0.0_r8 !sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(20)*frac_H_accum(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(21)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(20)*frac_H_accum(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(22)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(22)*frac_C_accum(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(23)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(23)*frac_F_accum(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(24)
       !cflx(i,idst) = sum( -dust_flux_in(i,:) )*(dust_emis_sclfctr(17)*frac_I_accum(i)+dust_emis_sclfctr(18)*frac_K_accum(i)+dust_emis_sclfctr(19)*frac_M_accum(i)+dust_emis_sclfctr(21)*frac_Q_accum(i)+dust_emis_sclfctr(24)*frac_G_accum(i))*soil_erod(i)/soil_erod_fact*1.15_r8
       cflx(i,idst) = ((sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(17)*frac_I_accum(i)*soil_erod(i)/soil_erod_fact*1.15_r8) + &
                       (sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(18)*frac_K_accum(i)*soil_erod(i)/soil_erod_fact*1.15_r8) + &
                       (sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(19)*frac_M_accum(i)*soil_erod(i)/soil_erod_fact*1.15_r8) + &
                       (sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(21)*frac_Q_accum(i)*soil_erod(i)/soil_erod_fact*1.15_r8) + &
                       (sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(24)*frac_G_accum(i)*soil_erod(i)/soil_erod_fact*1.15_r8))       
       idst = dust_indices(25)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(25)*frac_Fets_accum(i)*soil_erod(i)/soil_erod_fact*1.15_r8 
       idst = dust_indices(26)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(26)*frac_Fetm_accum(i)*soil_erod(i)/soil_erod_fact*1.15_r8 
       idst = dust_indices(27)
       cflx(i,idst) = 0.0_r8
       idst = dust_indices(28)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(28)*frac_Fesm_accum(i)*soil_erod(i)/soil_erod_fact*1.15_r8 
       idst = dust_indices(29)
       cflx(i,idst) = 0.0_r8 
       idst = dust_indices(30)
       cflx(i,idst) = 0.0_r8 
       idst = dust_indices(31)
       cflx(i,idst) = 0.0_r8
       idst = dust_indices(32)
       cflx(i,idst) = 0.0_r8

       ! Coarse mode minerals
       idst = dust_indices(33)
       cflx(i,idst) = 0.0_r8 !sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(33)*frac_I_coarse(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(34)
       cflx(i,idst) = 0.0_r8 !sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(34)*frac_K_coarse(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(35)
       cflx(i,idst) = 0.0_r8 !sum( -dust_flux_in(i,:) )*(dust_emis_sclfctr(33)+dust_emis_sclfctr(34)+dust_emis_sclfctr(35))*(frac_I_coarse(i)+frac_K_coarse(i)+frac_M_coarse(i))*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(36)
       cflx(i,idst) = 0.0_r8 !sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(36)*frac_H_coarse(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(37)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(36)*frac_H_coarse(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(38)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(38)*frac_C_coarse(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(39)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(39)*frac_F_coarse(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(40)
       !cflx(i,idst) = sum( -dust_flux_in(i,:) )*(dust_emis_sclfctr(33)*frac_I_coarse(i)+dust_emis_sclfctr(34)*frac_K_coarse(i)+dust_emis_sclfctr(35)*frac_M_coarse(i)+dust_emis_sclfctr(37)*frac_Q_coarse(i)+dust_emis_sclfctr(40)*frac_G_coarse(i))*soil_erod(i)/soil_erod_fact*1.15_r8
       cflx(i,idst) = ((sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(33)*frac_I_coarse(i)*soil_erod(i)/soil_erod_fact*1.15_r8) + &
                       (sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(34)*frac_K_coarse(i)*soil_erod(i)/soil_erod_fact*1.15_r8) + &
                       (sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(35)*frac_M_coarse(i)*soil_erod(i)/soil_erod_fact*1.15_r8) + &
                       (sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(37)*frac_Q_coarse(i)*soil_erod(i)/soil_erod_fact*1.15_r8) + &
                       (sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(40)*frac_G_coarse(i)*soil_erod(i)/soil_erod_fact*1.15_r8))
       idst = dust_indices(41)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(41)*frac_Fets_coarse(i)*soil_erod(i)/soil_erod_fact*1.15_r8 
       idst = dust_indices(42)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(42)*frac_Fetm_coarse(i)*soil_erod(i)/soil_erod_fact*1.15_r8 
       idst = dust_indices(43)
       cflx(i,idst) = 0.0_r8
       idst = dust_indices(44)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(44)*frac_Fesm_coarse(i)*soil_erod(i)/soil_erod_fact*1.15_r8 
       idst = dust_indices(45)
       cflx(i,idst) = 0.0_r8 
       idst = dust_indices(46)
       cflx(i,idst) = 0.0_r8 
       idst = dust_indices(47)
       cflx(i,idst) = 0.0_r8 
       idst = dust_indices(48)
       cflx(i,idst) = 0.0_r8

       x_mton = 6._r8 / (pi * dust_density * (dust_dmt_vwr(1)**3._r8))
       inum = dust_indices(1+dust_nbin)
       cflx(i,inum) =(cflx(i,dust_indices(1)) + cflx(i,dust_indices(2))  + cflx(i,dust_indices(3))  + cflx(i,dust_indices(4)) + &
                      cflx(i,dust_indices(5)) + cflx(i,dust_indices(6))  + cflx(i,dust_indices(7))  + cflx(i,dust_indices(8)) +  &
                      cflx(i,dust_indices(9)) + cflx(i,dust_indices(10)) + cflx(i,dust_indices(11)) + cflx(i,dust_indices(12))   &
                      )*x_mton

       x_mton = 6._r8 / (pi * dust_density * (dust_dmt_vwr(2)**3._r8))
       inum = dust_indices(2+dust_nbin)
       cflx(i,inum) =(cflx(i,dust_indices(17)) + cflx(i,dust_indices(18)) + cflx(i,dust_indices(19)) + cflx(i,dust_indices(20)) + &
                      cflx(i,dust_indices(21))+ cflx(i,dust_indices(22)) + cflx(i,dust_indices(23)) + cflx(i,dust_indices(24))  + &
                      cflx(i,dust_indices(25))+ cflx(i,dust_indices(26)) + cflx(i,dust_indices(27)) + cflx(i,dust_indices(28))   &
                       )*x_mton

       x_mton = 6._r8 / (pi * dust_density * (dust_dmt_vwr(3)**3._r8))
       inum = dust_indices(3+dust_nbin)
       cflx(i,inum) = (cflx(i,dust_indices(33))+ cflx(i,dust_indices(34)) + cflx(i,dust_indices(35)) + cflx(i,dust_indices(36)) + &
                       cflx(i,dust_indices(37))+ cflx(i,dust_indices(38)) + cflx(i,dust_indices(39)) + cflx(i,dust_indices(40)) +  &
                       cflx(i,dust_indices(41))+ cflx(i,dust_indices(42)) + cflx(i,dust_indices(43)) + cflx(i,dust_indices(44))   &
                        )*x_mton
    end do col_loop
#else
! do i change this ?????????????????????????????????????????????????
    col_loop: do i =1,ncol
       dust_flux_sum(i) = sum( -dust_flux_in(i,:) )

       !soil_erod(i) = soil_erodibility( i, lchnk )
       soil_erod(i) = 1._r8 ! Kok scheme

       if( soil_erod(i) .lt. soil_erod_threshold ) soil_erod(i) = 0._r8

       ! rebin and adjust dust emissons..
       ! Aitken mode minerals
       idst = dust_indices(1)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(1)*frac_I_aitken(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(2)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(2)*frac_K_aitken(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(3)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(3)*frac_M_aitken(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(4)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(4)*frac_H_aitken(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(5)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(5)*frac_Q_aitken(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(6)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(6)*frac_C_aitken(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(7)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(7)*frac_F_aitken(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(8)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(8)*frac_G_aitken(i)*soil_erod(i)/soil_erod_fact*1.15_r8

       ! Accumulation mode minerals
       idst = dust_indices(9)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(9)*frac_I_accum(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(10)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(10)*frac_K_accum(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(11)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(11)*frac_M_accum(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(12)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(12)*frac_H_accum(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(13)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(13)*frac_Q_accum(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(14)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(14)*frac_C_accum(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(15)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(15)*frac_F_accum(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(16)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(16)*frac_G_accum(i)*soil_erod(i)/soil_erod_fact*1.15_r8

       ! Coarse mode minerals
       idst = dust_indices(17)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(17)*frac_I_coarse(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(18)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(18)*frac_K_coarse(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(19)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(19)*frac_M_coarse(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(20)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(20)*frac_H_coarse(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(21)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(21)*frac_Q_coarse(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(22)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(22)*frac_C_coarse(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(23)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(23)*frac_F_coarse(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(24)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*dust_emis_sclfctr(24)*frac_G_coarse(i)*soil_erod(i)/soil_erod_fact*1.15_r8

       x_mton = 6._r8 / (pi * dust_density * (dust_dmt_vwr(1)**3._r8))
       inum = dust_indices(1+dust_nbin)
       cflx(i,inum) = (cflx(i,dust_indices(1)) + cflx(i,dust_indices(2)) + cflx(i,dust_indices(3)) + cflx(i,dust_indices(4)) + &
                       cflx(i,dust_indices(5)) + cflx(i,dust_indices(6)) + cflx(i,dust_indices(7)) + cflx(i,dust_indices(8))   &
                       )*x_mton

       x_mton = 6._r8 / (pi * dust_density * (dust_dmt_vwr(2)**3._r8))
       inum = dust_indices(2+dust_nbin)
       cflx(i,inum) = (cflx(i,dust_indices(9)) + cflx(i,dust_indices(10)) + cflx(i,dust_indices(11)) + cflx(i,dust_indices(12)) + &
                       cflx(i,dust_indices(13))+ cflx(i,dust_indices(14)) + cflx(i,dust_indices(15)) + cflx(i,dust_indices(16))   &
                       )*x_mton

       x_mton = 6._r8 / (pi * dust_density * (dust_dmt_vwr(3)**3._r8))
       inum = dust_indices(3+dust_nbin)
       cflx(i,inum) = (cflx(i,dust_indices(17))+ cflx(i,dust_indices(18)) + cflx(i,dust_indices(19)) + cflx(i,dust_indices(20)) + &
                       cflx(i,dust_indices(21))+ cflx(i,dust_indices(22)) + cflx(i,dust_indices(23)) + cflx(i,dust_indices(24))   &
                        )*x_mton
    end do col_loop
#endif
#else
    col_loop: do i =1,ncol

       soil_erod(i) = soil_erodibility( i, lchnk )

       if( soil_erod(i) .lt. soil_erod_threshold ) soil_erod(i) = 0._r8

       ! rebin and adjust dust emissons..
       do m = 1,dust_nbin

          idst = dust_indices(m)

          cflx(i,idst) = sum( -dust_flux_in(i,:) ) &
               * dust_emis_sclfctr(m)*soil_erod(i)/soil_erod_fact*1.15_r8

          x_mton = 6._r8 / (pi * dust_density * (dust_dmt_vwr(m)**3._r8))                

          inum = dust_indices(m+dust_nbin)

          cflx(i,inum) = cflx(i,idst)*x_mton

       enddo

    end do col_loop
#endif

#if ( defined DUSTTYPES )
  call outfld('dust_flux',dust_flux_sum(:),pcols,lchnk)
#endif
    !write(iulog,*) 'Longlei marks end in dust_model'
  end subroutine dust_emis

end module dust_model
