!-------------------------------------------------------------------
! Manages reading and interpolation of prescribed aerosol deposition 
! fluxes.  These are the deposition fluxes sent to the surface.
!
! Created by: Francis Vitt
!-------------------------------------------------------------------
module aerodep_flx

  use shr_kind_mod,     only : r8 => shr_kind_r8
  use cam_abortutils,   only : endrun
  use spmd_utils,       only : masterproc
  use tracer_data,      only : trfld, trfile
  use cam_logfile,      only : iulog
  use ppgrid,           only : pcols, pver, begchunk, endchunk

  implicit none
  private
  save 

  type(trfld), pointer :: fields(:)
  type(trfile)         :: file

  public :: aerodep_flx_init
  public :: aerodep_flx_adv
  public :: aerodep_flx_readnl
  public :: aerodep_flx_prescribed

  logical :: has_aerodep_flx = .false.
  integer, parameter, public :: N_BULK = 14
#if ( defined DUSTTYPES )
#if ( defined IRONPROC )
  integer, parameter, public :: N_MODAL = 206 !142 !126
#else
  integer, parameter, public :: N_MODAL = 78
#endif
#else
  integer, parameter, public :: N_MODAL = 22
#endif
  integer :: number_flds

  character(len=256) :: filename = 'NONE'
  character(len=256) :: filelist = ' '
  character(len=256) :: datapath = ' '
  character(len=32)  :: datatype = 'SERIAL'
  logical            :: rmv_file = .false.
  integer            :: cycle_yr = 0
  integer            :: fixed_ymd = 0
  integer            :: fixed_tod = 0
  character(len=32)  :: specifier(N_MODAL) = ' '

  ! for bulk aerosol fluxes

  character(len=12), parameter :: bulk_names(N_BULK) = (/ &
       'BCDEPWET    ', 'BCPHODRY    ', 'BCPHIDRY    ',  &
       'OCDEPWET    ', 'OCPHODRY    ', 'OCPHIDRY    ',  &
       'DSTX01DD    ', 'DSTX02DD    ', 'DSTX03DD    ', 'DSTX04DD    ', &
       'DSTX01WD    ', 'DSTX02WD    ', 'DSTX03WD    ', 'DSTX04WD    ' /)

  integer :: index_bulk_map(N_BULK)

  integer :: ibcphiwet,ibcphidry,ibcphodry
  integer :: iocphiwet,iocphidry,iocphodry

  integer :: idstdry1,idstdry2,idstdry3,idstdry4
  integer :: idstwet1,idstwet2,idstwet3,idstwet4

  ! for modal aerosol fluxes
#if ( defined DUSTTYPES )
#if ( defined IRONPROC )
  character(len=12), parameter :: modal_names(N_MODAL) = (/ &
       'bc_a1DDF    ', 'bc_c1DDF    ', 'pom_a1DDF   ', 'pom_c1DDF   ',  &
       'soa_a1DDF   ', 'soa_c1DDF   ', 'soa_a2DDF   ', 'soa_c2DDF   ',  &
       'dst1_a1DDF  ', 'dst1_c1DDF  ', 'dst1_a3DDF  ', 'dst1_c3DDF  ',  &
       'dst2_a1DDF  ', 'dst2_c1DDF  ', 'dst2_a3DDF  ', 'dst2_c3DDF  ',  &
       'dst3_a1DDF  ', 'dst3_c1DDF  ', 'dst3_a3DDF  ', 'dst3_c3DDF  ',  &
       'dst4_a1DDF  ', 'dst4_c1DDF  ', 'dst4_a3DDF  ', 'dst4_c3DDF  ',  &
       'dst5_a1DDF  ', 'dst5_c1DDF  ', 'dst5_a3DDF  ', 'dst5_c3DDF  ',  &
       'dst6_a1DDF  ', 'dst6_c1DDF  ', 'dst6_a3DDF  ', 'dst6_c3DDF  ',  &
       'dst7_a1DDF  ', 'dst7_c1DDF  ', 'dst7_a3DDF  ', 'dst7_c3DDF  ',  &
       'dst8_a1DDF  ', 'dst8_c1DDF  ', 'dst8_a3DDF  ', 'dst8_c3DDF  ',  &
       'dst1_a2DDF  ', 'dst1_c2DDF  ', &
       'dst2_a2DDF  ', 'dst2_c2DDF  ', &
       'dst3_a2DDF  ', 'dst3_c2DDF  ', &
       'dst4_a2DDF  ', 'dst4_c2DDF  ', &
       'dst5_a2DDF  ', 'dst5_c2DDF  ', &
       'dst6_a2DDF  ', 'dst6_c2DDF  ', &
       'dst7_a2DDF  ', 'dst7_c2DDF  ', &
       'dst8_a2DDF  ', 'dst8_c2DDF  ', &
       'bc_a1SFWET  ', 'bc_c1SFWET  ', 'pom_a1SFWET ', 'pom_c1SFWET ',  &
       'soa_a1SFWET ', 'soa_c1SFWET ', 'dst1_a1SFWET', 'dst1_c1SFWET',  &
       'dst1_a3SFWET', 'dst1_c3SFWET', 'dst2_a1SFWET', 'dst2_c1SFWET',  &
       'dst2_a3SFWET', 'dst2_c3SFWET', 'dst3_a1SFWET', 'dst3_c1SFWET',  &
       'dst3_a3SFWET', 'dst3_c3SFWET', 'dst4_a1SFWET', 'dst4_c1SFWET',  &
       'dst4_a3SFWET', 'dst4_c3SFWET', 'dst5_a1SFWET', 'dst5_c1SFWET',  &
       'dst5_a3SFWET', 'dst5_c3SFWET', 'dst6_a1SFWET', 'dst6_c1SFWET',  &
       'dst6_a3SFWET', 'dst6_c3SFWET', 'dst7_a1SFWET', 'dst7_c1SFWET',  &
       'dst7_a3SFWET', 'dst7_c3SFWET', 'dst8_a1SFWET', 'dst8_c1SFWET',  &
       'dst8_a3SFWET', 'dst8_c3SFWET', &
       'dst1_a2SFWET', 'dst1_c2SFWET', &
       'dst2_a2SFWET', 'dst2_c2SFWET', &
       'dst3_a2SFWET', 'dst3_c2SFWET', &
       'dst4_a2SFWET', 'dst4_c2SFWET', &
       'dst5_a2SFWET', 'dst5_c2SFWET', &
       'dst6_a2SFWET', 'dst6_c2SFWET', &
       'dst7_a2SFWET', 'dst7_c2SFWET', &
       'dst8_a2SFWET', 'dst8_c2SFWET', &
       'fets_a1DDF  ', 'fets_c1DDF  ', 'fets_a3DDF  ', 'fets_c3DDF  ', &
       'fetm_a1DDF  ', 'fetm_c1DDF  ', 'fetm_a3DDF  ', 'fetm_c3DDF  ', &
       'fess_a1DDF  ', 'fess_c1DDF  ', 'fess_a3DDF  ', 'fess_c3DDF  ', &
       'fesm_a1DDF  ', 'fesm_c1DDF  ', 'fesm_a3DDF  ', 'fesm_c3DDF  ', &
       'fetc_a1DDF  ', 'fetc_c1DDF  ', 'fetc_a3DDF  ', 'fetc_c3DDF  ', &
       'fesc_a1DDF  ', 'fesc_c1DDF  ', 'fesc_a3DDF  ', 'fesc_c3DDF  ', &
       'fetf_a1DDF  ', 'fetf_c1DDF  ', 'fetf_a3DDF  ', 'fetf_c3DDF  ', &
       'fesf_a1DDF  ', 'fesf_c1DDF  ', 'fesf_a3DDF  ', 'fesf_c3DDF  ', &
       'fets_a2DDF  ', 'fets_c2DDF  ', & 
       'fetm_a2DDF  ', 'fetm_c2DDF  ', &
       'fess_a2DDF  ', 'fess_c2DDF  ', &
       'fesm_a2DDF  ', 'fesm_c2DDF  ', &
       'fetc_a2DDF  ', 'fetc_c2DDF  ', &
       'fesc_a2DDF  ', 'fesc_c2DDF  ', &
       'fetf_a2DDF  ', 'fetf_c2DDF  ', &
       'fesf_a2DDF  ', 'fesf_c2DDF  ', &
       'fets_a1SFWET', 'fets_c1SFWET', 'fets_a3SFWET', 'fets_c3SFWET', &
       'fetm_a1SFWET', 'fetm_c1SFWET', 'fetm_a3SFWET', 'fetm_c3SFWET', &
       'fess_a1SFWET', 'fess_c1SFWET', 'fess_a3SFWET', 'fess_c3SFWET', &
       'fesm_a1SFWET', 'fesm_c1SFWET', 'fesm_a3SFWET', 'fesm_c3SFWET', &
       'fetc_a1SFWET', 'fetc_c1SFWET', 'fetc_a3SFWET', 'fetc_c3SFWET', &
       'fesc_a1SFWET', 'fesc_c1SFWET', 'fesc_a3SFWET', 'fesc_c3SFWET', &
       'fetf_a1SFWET', 'fetf_c1SFWET', 'fetf_a3SFWET', 'fetf_c3SFWET', &
       'fesf_a1SFWET', 'fesf_c1SFWET', 'fesf_a3SFWET', 'fesf_c3SFWET', &
       'fets_a2SFWET', 'fets_c2SFWET', &
       'fetm_a2SFWET', 'fetm_c2SFWET', &
       'fess_a2SFWET', 'fess_c2SFWET', &
       'fesm_a2SFWET', 'fesm_c2SFWET', &
       'fetc_a2SFWET', 'fetc_c2SFWET', &
       'fesc_a2SFWET', 'fesc_c2SFWET', &
       'fetf_a2SFWET', 'fetf_c2SFWET', &
       'fesf_a2SFWET', 'fesf_c2SFWET' /)
#else
  character(len=12), parameter :: modal_names(N_MODAL) = (/ &
       'bc_a1DDF    ', 'bc_c1DDF    ', 'pom_a1DDF   ', 'pom_c1DDF   ',  &
       'soa_a1DDF   ', 'soa_c1DDF   ', 'soa_a2DDF   ', 'soa_c2DDF   ',  &
       'dst1_a1DDF  ', 'dst1_c1DDF  ', 'dst1_a3DDF  ', 'dst1_c3DDF  ',  &
       'dst2_a1DDF  ', 'dst2_c1DDF  ', 'dst2_a3DDF  ', 'dst2_c3DDF  ',  &
       'dst3_a1DDF  ', 'dst3_c1DDF  ', 'dst3_a3DDF  ', 'dst3_c3DDF  ',  &
       'dst4_a1DDF  ', 'dst4_c1DDF  ', 'dst4_a3DDF  ', 'dst4_c3DDF  ',  &
       'dst5_a1DDF  ', 'dst5_c1DDF  ', 'dst5_a3DDF  ', 'dst5_c3DDF  ',  &
       'dst6_a1DDF  ', 'dst6_c1DDF  ', 'dst6_a3DDF  ', 'dst6_c3DDF  ',  &
       'dst7_a1DDF  ', 'dst7_c1DDF  ', 'dst7_a3DDF  ', 'dst7_c3DDF  ',  &
       'dst8_a1DDF  ', 'dst8_c1DDF  ', 'dst8_a3DDF  ', 'dst8_c3DDF  ',  &
       'bc_a1SFWET  ', 'bc_c1SFWET  ', 'pom_a1SFWET ', 'pom_c1SFWET ',  &
       'soa_a1SFWET ', 'soa_c1SFWET ', 'dst1_a1SFWET', 'dst1_c1SFWET',  &
       'dst1_a3SFWET', 'dst1_c3SFWET', 'dst2_a1SFWET', 'dst2_c1SFWET',  &
       'dst2_a3SFWET', 'dst2_c3SFWET', 'dst3_a1SFWET', 'dst3_c1SFWET',  &
       'dst3_a3SFWET', 'dst3_c3SFWET', 'dst4_a1SFWET', 'dst4_c1SFWET',  &
       'dst4_a3SFWET', 'dst4_c3SFWET', 'dst5_a1SFWET', 'dst5_c1SFWET',  &
       'dst5_a3SFWET', 'dst5_c3SFWET', 'dst6_a1SFWET', 'dst6_c1SFWET',  &
       'dst6_a3SFWET', 'dst6_c3SFWET', 'dst7_a1SFWET', 'dst7_c1SFWET',  &
       'dst7_a3SFWET', 'dst7_c3SFWET', 'dst8_a1SFWET', 'dst8_c1SFWET',  &
       'dst8_a3SFWET', 'dst8_c3SFWET'/)
#endif
#else
  character(len=12), parameter :: modal_names(N_MODAL) = (/ &
       'bc_a1DDF    ', 'bc_c1DDF    ', 'pom_a1DDF   ', 'pom_c1DDF   ',  &
       'soa_a1DDF   ', 'soa_c1DDF   ', 'soa_a2DDF   ', 'soa_c2DDF   ',  &
       'dst_a1DDF   ', 'dst_c1DDF   ', 'dst_a3DDF   ', 'dst_c3DDF   ',  &
       'bc_a1SFWET  ', 'bc_c1SFWET  ', 'pom_a1SFWET ', 'pom_c1SFWET ',  &
       'soa_a1SFWET ', 'soa_c1SFWET ', 'dst_a1SFWET ', 'dst_c1SFWET ',  &
       'dst_a3SFWET ', 'dst_c3SFWET ' /)
#endif

  integer :: index_modal_map(N_MODAL)

#if ( defined DUSTTYPES )
#if ( defined IRONPROC )
  integer, parameter :: idx_bc1 = 1
  integer, parameter :: idx_pom1 = 2
  integer, parameter :: idx_soa1 = 3
  integer, parameter :: idx_soa2 = 4
  integer, parameter :: idx_dst11 = 5
  integer, parameter :: idx_dst21 = 6
  integer, parameter :: idx_dst31 = 7
  integer, parameter :: idx_dst41 = 8
  integer, parameter :: idx_dst51 = 9
  integer, parameter :: idx_dst61 = 10
  integer, parameter :: idx_dst71 = 11
  integer, parameter :: idx_dst81 = 12
  integer, parameter :: idx_dst13 = 13
  integer, parameter :: idx_dst23 = 14
  integer, parameter :: idx_dst33 = 15
  integer, parameter :: idx_dst43 = 16
  integer, parameter :: idx_dst53 = 17
  integer, parameter :: idx_dst63 = 18
  integer, parameter :: idx_dst73 = 19
  integer, parameter :: idx_dst83 = 20
  integer, parameter :: idx_dst12 = 21
  integer, parameter :: idx_dst22 = 22
  integer, parameter :: idx_dst32 = 23
  integer, parameter :: idx_dst42 = 24
  integer, parameter :: idx_dst52 = 25
  integer, parameter :: idx_dst62 = 26
  integer, parameter :: idx_dst72 = 27
  integer, parameter :: idx_dst82 = 28
  integer, parameter :: idx_ncl3  = 29
  integer, parameter :: idx_so43  = 30
  integer, parameter :: idx_fets1 = 31
  integer, parameter :: idx_fetm1 = 32
  integer, parameter :: idx_fess1 = 33
  integer, parameter :: idx_fesm1 = 34
  integer, parameter :: idx_fetc1 = 35
  integer, parameter :: idx_fesc1 = 36
  integer, parameter :: idx_fetf1 = 37
  integer, parameter :: idx_fesf1 = 38
  integer, parameter :: idx_fets3 = 39
  integer, parameter :: idx_fetm3 = 40
  integer, parameter :: idx_fess3 = 41
  integer, parameter :: idx_fesm3 = 42
  integer, parameter :: idx_fesc3 = 43
  integer, parameter :: idx_fetc3 = 44
  integer, parameter :: idx_fetf3 = 45
  integer, parameter :: idx_fesf3 = 46
  integer, parameter :: idx_fetm2 = 47
  integer, parameter :: idx_fets2 = 48
  integer, parameter :: idx_fess2 = 49
  integer, parameter :: idx_fesm2 = 50
  integer, parameter :: idx_fetc2 = 51
  integer, parameter :: idx_fesc2 = 52
  integer, parameter :: idx_fetf2 = 53
  integer, parameter :: idx_fesf2 = 54
#else
  integer, parameter :: idx_bc1 = 1
  integer, parameter :: idx_pom1 = 2
  integer, parameter :: idx_soa1 = 3
  integer, parameter :: idx_soa2 = 4
  integer, parameter :: idx_dst11 = 5
  integer, parameter :: idx_dst21 = 6
  integer, parameter :: idx_dst31 = 7
  integer, parameter :: idx_dst41 = 8
  integer, parameter :: idx_dst51 = 9
  integer, parameter :: idx_dst61 = 10
  integer, parameter :: idx_dst71 = 11
  integer, parameter :: idx_dst81 = 12
  integer, parameter :: idx_dst13 = 13
  integer, parameter :: idx_dst23 = 14
  integer, parameter :: idx_dst33 = 15
  integer, parameter :: idx_dst43 = 16
  integer, parameter :: idx_dst53 = 17
  integer, parameter :: idx_dst63 = 18
  integer, parameter :: idx_dst73 = 19
  integer, parameter :: idx_dst83 = 20
  integer, parameter :: idx_ncl3 = 21
  integer, parameter :: idx_so43 = 22
#endif
#else
  integer, parameter :: idx_bc1 = 1
  integer, parameter :: idx_pom1 = 2
  integer, parameter :: idx_soa1 = 3
  integer, parameter :: idx_soa2 = 4
  integer, parameter :: idx_dst1 = 5
  integer, parameter :: idx_dst3 = 6
  integer, parameter :: idx_ncl3 = 7
  integer, parameter :: idx_so43 = 8
#endif

#if ( defined DUSTTYPES )
#if ( defined IRONPROC )
  integer, parameter :: nmodal_idxs = 54 !34
#else
  integer, parameter :: nmodal_idxs = 22
#endif
#else
  integer, parameter :: nmodal_idxs = 8
#endif

  integer :: idx_bc1_dryis = -1
  integer :: idx_bc1_drycw = -1
  integer :: idx_pom1_dryis = -1
  integer :: idx_pom1_drycw = -1
  integer :: idx_soa1_dryis = -1
  integer :: idx_soa1_drycw = -1
  integer :: idx_soa2_dryis = -1
  integer :: idx_soa2_drycw = -1
#if ( defined DUSTTYPES )
#if ( defined IRONPROC )
  integer :: idx_dst11_dryis = -1
  integer :: idx_dst11_drycw = -1
  integer :: idx_dst21_dryis = -1
  integer :: idx_dst21_drycw = -1
  integer :: idx_dst31_dryis = -1
  integer :: idx_dst31_drycw = -1
  integer :: idx_dst41_dryis = -1
  integer :: idx_dst41_drycw = -1
  integer :: idx_dst51_dryis = -1
  integer :: idx_dst51_drycw = -1
  integer :: idx_dst61_dryis = -1
  integer :: idx_dst61_drycw = -1
  integer :: idx_dst71_dryis = -1
  integer :: idx_dst71_drycw = -1
  integer :: idx_dst81_dryis = -1
  integer :: idx_dst81_drycw = -1
  integer :: idx_fets1_dryis = -1
  integer :: idx_fets1_drycw = -1
  integer :: idx_fetm1_dryis = -1
  integer :: idx_fetm1_drycw = -1
  integer :: idx_fess1_dryis = -1
  integer :: idx_fess1_drycw = -1
  integer :: idx_fesm1_dryis = -1
  integer :: idx_fesm1_drycw = -1
  integer :: idx_fetc1_dryis = -1
  integer :: idx_fetc1_drycw = -1
  integer :: idx_fesc1_dryis = -1
  integer :: idx_fesc1_drycw = -1
  integer :: idx_fetf1_dryis = -1
  integer :: idx_fetf1_drycw = -1
  integer :: idx_fesf1_dryis = -1
  integer :: idx_fesf1_drycw = -1

  integer :: idx_dst12_dryis = -1
  integer :: idx_dst12_drycw = -1
  integer :: idx_dst22_dryis = -1
  integer :: idx_dst22_drycw = -1
  integer :: idx_dst32_dryis = -1
  integer :: idx_dst32_drycw = -1
  integer :: idx_dst42_dryis = -1
  integer :: idx_dst42_drycw = -1
  integer :: idx_dst52_dryis = -1
  integer :: idx_dst52_drycw = -1
  integer :: idx_dst62_dryis = -1
  integer :: idx_dst62_drycw = -1
  integer :: idx_dst72_dryis = -1
  integer :: idx_dst72_drycw = -1
  integer :: idx_dst82_dryis = -1
  integer :: idx_dst82_drycw = -1
  integer :: idx_fets2_dryis = -1
  integer :: idx_fets2_drycw = -1
  integer :: idx_fetm2_dryis = -1
  integer :: idx_fetm2_drycw = -1
  integer :: idx_fess2_dryis = -1
  integer :: idx_fess2_drycw = -1
  integer :: idx_fesm2_dryis = -1
  integer :: idx_fesm2_drycw = -1
  integer :: idx_fetc2_dryis = -1
  integer :: idx_fetc2_drycw = -1
  integer :: idx_fesc2_dryis = -1
  integer :: idx_fesc2_drycw = -1
  integer :: idx_fetf2_dryis = -1
  integer :: idx_fetf2_drycw = -1
  integer :: idx_fesf2_dryis = -1
  integer :: idx_fesf2_drycw = -1

  integer :: idx_dst13_dryis = -1
  integer :: idx_dst13_drycw = -1
  integer :: idx_dst23_dryis = -1
  integer :: idx_dst23_drycw = -1
  integer :: idx_dst33_dryis = -1
  integer :: idx_dst33_drycw = -1
  integer :: idx_dst43_dryis = -1
  integer :: idx_dst43_drycw = -1
  integer :: idx_dst53_dryis = -1
  integer :: idx_dst53_drycw = -1
  integer :: idx_dst63_dryis = -1
  integer :: idx_dst63_drycw = -1
  integer :: idx_dst73_dryis = -1
  integer :: idx_dst73_drycw = -1
  integer :: idx_dst83_dryis = -1
  integer :: idx_dst83_drycw = -1
  integer :: idx_fets3_dryis = -1
  integer :: idx_fets3_drycw = -1
  integer :: idx_fetm3_dryis = -1
  integer :: idx_fetm3_drycw = -1
  integer :: idx_fess3_dryis = -1
  integer :: idx_fess3_drycw = -1
  integer :: idx_fesm3_dryis = -1
  integer :: idx_fesm3_drycw = -1
  integer :: idx_fetc3_dryis = -1
  integer :: idx_fetc3_drycw = -1
  integer :: idx_fesc3_dryis = -1
  integer :: idx_fesc3_drycw = -1
  integer :: idx_fetf3_dryis = -1
  integer :: idx_fetf3_drycw = -1
  integer :: idx_fesf3_dryis = -1
  integer :: idx_fesf3_drycw = -1
#else
  integer :: idx_dst11_dryis = -1
  integer :: idx_dst11_drycw = -1
  integer :: idx_dst21_dryis = -1
  integer :: idx_dst21_drycw = -1
  integer :: idx_dst31_dryis = -1
  integer :: idx_dst31_drycw = -1
  integer :: idx_dst41_dryis = -1
  integer :: idx_dst41_drycw = -1
  integer :: idx_dst51_dryis = -1
  integer :: idx_dst51_drycw = -1
  integer :: idx_dst61_dryis = -1
  integer :: idx_dst61_drycw = -1
  integer :: idx_dst71_dryis = -1
  integer :: idx_dst71_drycw = -1
  integer :: idx_dst81_dryis = -1
  integer :: idx_dst81_drycw = -1

  integer :: idx_dst13_dryis = -1
  integer :: idx_dst13_drycw = -1
  integer :: idx_dst23_dryis = -1
  integer :: idx_dst23_drycw = -1
  integer :: idx_dst33_dryis = -1
  integer :: idx_dst33_drycw = -1
  integer :: idx_dst43_dryis = -1
  integer :: idx_dst43_drycw = -1
  integer :: idx_dst53_dryis = -1
  integer :: idx_dst53_drycw = -1
  integer :: idx_dst63_dryis = -1
  integer :: idx_dst63_drycw = -1
  integer :: idx_dst73_dryis = -1
  integer :: idx_dst73_drycw = -1
  integer :: idx_dst83_dryis = -1
  integer :: idx_dst83_drycw = -1
#endif
#else  
  integer :: idx_dst1_dryis = -1
  integer :: idx_dst1_drycw = -1
  integer :: idx_dst3_dryis = -1
  integer :: idx_dst3_drycw = -1
#endif

  integer :: idx_bc1_wetis = -1
  integer :: idx_bc1_wetcw = -1
  integer :: idx_pom1_wetis = -1
  integer :: idx_pom1_wetcw = -1
  integer :: idx_soa1_wetis = -1
  integer :: idx_soa1_wetcw = -1
#if ( defined DUSTTYPES )
#if ( defined IRONPROC )
  integer :: idx_dst11_wetis = -1
  integer :: idx_dst11_wetcw = -1
  integer :: idx_dst21_wetis = -1
  integer :: idx_dst21_wetcw = -1
  integer :: idx_dst31_wetis = -1
  integer :: idx_dst31_wetcw = -1
  integer :: idx_dst41_wetis = -1
  integer :: idx_dst41_wetcw = -1
  integer :: idx_dst51_wetis = -1
  integer :: idx_dst51_wetcw = -1
  integer :: idx_dst61_wetis = -1
  integer :: idx_dst61_wetcw = -1
  integer :: idx_dst71_wetis = -1
  integer :: idx_dst71_wetcw = -1
  integer :: idx_dst81_wetis = -1
  integer :: idx_dst81_wetcw = -1
  integer :: idx_fets1_wetis = -1
  integer :: idx_fets1_wetcw = -1
  integer :: idx_fetm1_wetis = -1
  integer :: idx_fetm1_wetcw = -1
  integer :: idx_fess1_wetis = -1
  integer :: idx_fess1_wetcw = -1
  integer :: idx_fesm1_wetis = -1
  integer :: idx_fesm1_wetcw = -1
  integer :: idx_fetc1_wetis = -1
  integer :: idx_fetc1_wetcw = -1
  integer :: idx_fesc1_wetis = -1
  integer :: idx_fesc1_wetcw = -1
  integer :: idx_fetf1_wetis = -1
  integer :: idx_fetf1_wetcw = -1
  integer :: idx_fesf1_wetis = -1
  integer :: idx_fesf1_wetcw = -1

  integer :: idx_dst12_wetis = -1
  integer :: idx_dst12_wetcw = -1
  integer :: idx_dst22_wetis = -1
  integer :: idx_dst22_wetcw = -1
  integer :: idx_dst32_wetis = -1
  integer :: idx_dst32_wetcw = -1
  integer :: idx_dst42_wetis = -1
  integer :: idx_dst42_wetcw = -1
  integer :: idx_dst52_wetis = -1
  integer :: idx_dst52_wetcw = -1
  integer :: idx_dst62_wetis = -1
  integer :: idx_dst62_wetcw = -1
  integer :: idx_dst72_wetis = -1
  integer :: idx_dst72_wetcw = -1
  integer :: idx_dst82_wetis = -1
  integer :: idx_dst82_wetcw = -1
  integer :: idx_fets2_wetis = -1
  integer :: idx_fets2_wetcw = -1
  integer :: idx_fetm2_wetis = -1
  integer :: idx_fetm2_wetcw = -1
  integer :: idx_fess2_wetis = -1
  integer :: idx_fess2_wetcw = -1
  integer :: idx_fesm2_wetis = -1
  integer :: idx_fesm2_wetcw = -1
  integer :: idx_fetc2_wetis = -1
  integer :: idx_fetc2_wetcw = -1
  integer :: idx_fesc2_wetis = -1
  integer :: idx_fesc2_wetcw = -1
  integer :: idx_fetf2_wetis = -1
  integer :: idx_fetf2_wetcw = -1
  integer :: idx_fesf2_wetis = -1
  integer :: idx_fesf2_wetcw = -1

  integer :: idx_dst13_wetis = -1
  integer :: idx_dst13_wetcw = -1
  integer :: idx_dst23_wetis = -1
  integer :: idx_dst23_wetcw = -1
  integer :: idx_dst33_wetis = -1
  integer :: idx_dst33_wetcw = -1
  integer :: idx_dst43_wetis = -1
  integer :: idx_dst43_wetcw = -1
  integer :: idx_dst53_wetis = -1
  integer :: idx_dst53_wetcw = -1
  integer :: idx_dst63_wetis = -1
  integer :: idx_dst63_wetcw = -1
  integer :: idx_dst73_wetis = -1
  integer :: idx_dst73_wetcw = -1
  integer :: idx_dst83_wetis = -1
  integer :: idx_dst83_wetcw = -1
  integer :: idx_fets3_wetis = -1
  integer :: idx_fets3_wetcw = -1
  integer :: idx_fetm3_wetis = -1
  integer :: idx_fetm3_wetcw = -1
  integer :: idx_fess3_wetis = -1
  integer :: idx_fess3_wetcw = -1
  integer :: idx_fesm3_wetis = -1
  integer :: idx_fesm3_wetcw = -1
  integer :: idx_fetc3_wetis = -1
  integer :: idx_fetc3_wetcw = -1
  integer :: idx_fesc3_wetis = -1
  integer :: idx_fesc3_wetcw = -1
  integer :: idx_fetf3_wetis = -1
  integer :: idx_fetf3_wetcw = -1
  integer :: idx_fesf3_wetis = -1
  integer :: idx_fesf3_wetcw = -1
#else
  integer :: idx_dst11_wetis = -1
  integer :: idx_dst11_wetcw = -1
  integer :: idx_dst21_wetis = -1
  integer :: idx_dst21_wetcw = -1
  integer :: idx_dst31_wetis = -1
  integer :: idx_dst31_wetcw = -1
  integer :: idx_dst41_wetis = -1
  integer :: idx_dst41_wetcw = -1
  integer :: idx_dst51_wetis = -1
  integer :: idx_dst51_wetcw = -1
  integer :: idx_dst61_wetis = -1
  integer :: idx_dst61_wetcw = -1
  integer :: idx_dst71_wetis = -1
  integer :: idx_dst71_wetcw = -1
  integer :: idx_dst81_wetis = -1
  integer :: idx_dst81_wetcw = -1

  integer :: idx_dst13_wetis = -1
  integer :: idx_dst13_wetcw = -1
  integer :: idx_dst23_wetis = -1
  integer :: idx_dst23_wetcw = -1
  integer :: idx_dst33_wetis = -1
  integer :: idx_dst33_wetcw = -1
  integer :: idx_dst43_wetis = -1
  integer :: idx_dst43_wetcw = -1
  integer :: idx_dst53_wetis = -1
  integer :: idx_dst53_wetcw = -1
  integer :: idx_dst63_wetis = -1
  integer :: idx_dst63_wetcw = -1
  integer :: idx_dst73_wetis = -1
  integer :: idx_dst73_wetcw = -1
  integer :: idx_dst83_wetis = -1
  integer :: idx_dst83_wetcw = -1
#endif
#else  
  integer :: idx_dst1_wetis = -1
  integer :: idx_dst1_wetcw = -1
  integer :: idx_dst3_wetis = -1
  integer :: idx_dst3_wetcw = -1
#endif
  logical :: modal_fluxes = .false.

contains

!-------------------------------------------------------------------
! parses the list of dep fluxes specified in aerodep_flx_specifier namelist
! variable and sets up index variables
!-------------------------------------------------------------------
  subroutine aerodep_flx_init()

    use tracer_data, only : trcdata_init
    use cam_history, only : addfld, horiz_only
    use physics_buffer, only : physics_buffer_desc
    use modal_aero_deposition, only : modal_aero_deposition_init

    implicit none

    integer :: ndx, istat, i

    if ( has_aerodep_flx ) then
       if ( masterproc ) then
          write(iulog,*) 'aero dep fluxes are prescribed in :'//trim(filename)
       endif
    else
       return
    endif

    allocate(file%in_pbuf(size(specifier)))
    file%in_pbuf(:) = .false.
    call trcdata_init( specifier, filename, filelist, datapath, fields, file, &
                       rmv_file, cycle_yr, fixed_ymd, fixed_tod, datatype)

    number_flds = 0
    if (associated(fields)) number_flds = size( fields )

    if( number_flds < 1 ) then
       has_aerodep_flx = .false.
       if (masterproc) then
          write(iulog,*) 'aerodep_flx_init: no aerosol deposition fluxes have been specified'
       endif
       return
    end if

    index_bulk_map(:) = -1
    index_modal_map(:) = -1

    do i = 1,number_flds

       ndx = get_ndx( fields(i)%fldnam, bulk_names )
       if (ndx >0) then
          index_bulk_map(ndx) = i
       else
          ndx = get_ndx( fields(i)%fldnam, modal_names )
          if (ndx >0) then
             index_modal_map(ndx) = i
          endif
       endif
       if (ndx>0) then
          call addfld(trim(fields(i)%fldnam)//'_D', horiz_only, 'A',fields(i)%units, 'prescribed aero dep' )
       else
          call endrun('aerodep_flx_init: aerosol flux name not recognized: '//trim(fields(i)%fldnam))
       endif
    enddo

    modal_fluxes = any(index_modal_map(:)>0)

    if (modal_fluxes) then

#if (defined DUSTTYPES)
#if (defined IRONPROC)
       idx_bc1_dryis  = index_modal_map(1)
       idx_bc1_drycw  = index_modal_map(2)
       idx_pom1_dryis = index_modal_map(3)
       idx_pom1_drycw = index_modal_map(4)
       idx_soa1_dryis = index_modal_map(5)
       idx_soa1_drycw = index_modal_map(6)
       idx_soa2_dryis = index_modal_map(7)
       idx_soa2_drycw = index_modal_map(8)
       idx_dst11_dryis = index_modal_map(9)
       idx_dst11_drycw = index_modal_map(10)
       idx_dst21_dryis = index_modal_map(11)
       idx_dst21_drycw = index_modal_map(12)
       idx_dst31_dryis = index_modal_map(13)
       idx_dst31_drycw = index_modal_map(14)
       idx_dst41_dryis = index_modal_map(15)
       idx_dst41_drycw = index_modal_map(16)
       idx_dst51_dryis = index_modal_map(17)
       idx_dst51_drycw = index_modal_map(18)
       idx_dst61_dryis = index_modal_map(19)
       idx_dst61_drycw = index_modal_map(20)
       idx_dst71_dryis = index_modal_map(21)
       idx_dst71_drycw = index_modal_map(22)
       idx_dst81_dryis = index_modal_map(23)
       idx_dst81_drycw = index_modal_map(24)

       idx_dst13_dryis = index_modal_map(25)
       idx_dst13_drycw = index_modal_map(26)
       idx_dst23_dryis = index_modal_map(27)
       idx_dst23_drycw = index_modal_map(28)
       idx_dst33_dryis = index_modal_map(29)
       idx_dst33_drycw = index_modal_map(30)
       idx_dst43_dryis = index_modal_map(31)
       idx_dst43_drycw = index_modal_map(32)
       idx_dst53_dryis = index_modal_map(33)
       idx_dst53_drycw = index_modal_map(34)
       idx_dst63_dryis = index_modal_map(35)
       idx_dst63_drycw = index_modal_map(36)
       idx_dst73_dryis = index_modal_map(37)
       idx_dst73_drycw = index_modal_map(38)
       idx_dst83_dryis = index_modal_map(39)
       idx_dst83_drycw = index_modal_map(40)

       idx_dst12_dryis = index_modal_map(41)
       idx_dst12_drycw = index_modal_map(42)
       idx_dst22_dryis = index_modal_map(43)
       idx_dst22_drycw = index_modal_map(44)
       idx_dst32_dryis = index_modal_map(45)
       idx_dst32_drycw = index_modal_map(46)
       idx_dst42_dryis = index_modal_map(47)
       idx_dst42_drycw = index_modal_map(48)
       idx_dst52_dryis = index_modal_map(49)
       idx_dst52_drycw = index_modal_map(50)
       idx_dst62_dryis = index_modal_map(51)
       idx_dst62_drycw = index_modal_map(52)
       idx_dst72_dryis = index_modal_map(53)
       idx_dst72_drycw = index_modal_map(54)
       idx_dst82_dryis = index_modal_map(55)
       idx_dst82_drycw = index_modal_map(56)

       idx_bc1_wetis   = index_modal_map(57)
       idx_bc1_wetcw   = index_modal_map(58)
       idx_pom1_wetis  = index_modal_map(59)
       idx_pom1_wetcw  = index_modal_map(60)
       idx_soa1_wetis  = index_modal_map(61)
       idx_soa1_wetcw  = index_modal_map(62)
       idx_dst11_wetis = index_modal_map(63)
       idx_dst11_wetcw = index_modal_map(64)
       idx_dst21_wetis = index_modal_map(65)
       idx_dst21_wetcw = index_modal_map(66)
       idx_dst31_wetis = index_modal_map(67)
       idx_dst31_wetcw = index_modal_map(68)
       idx_dst41_wetis = index_modal_map(69)
       idx_dst41_wetcw = index_modal_map(70)
       idx_dst51_wetis = index_modal_map(71)
       idx_dst51_wetcw = index_modal_map(72)
       idx_dst61_wetis = index_modal_map(73)
       idx_dst61_wetcw = index_modal_map(74)
       idx_dst71_wetis = index_modal_map(75)
       idx_dst71_wetcw = index_modal_map(76)
       idx_dst81_wetis = index_modal_map(77)
       idx_dst81_wetcw = index_modal_map(78)

       idx_dst13_wetis = index_modal_map(79)
       idx_dst13_wetcw = index_modal_map(80)
       idx_dst23_wetis = index_modal_map(81)
       idx_dst23_wetcw = index_modal_map(82)
       idx_dst33_wetis = index_modal_map(83)
       idx_dst33_wetcw = index_modal_map(84)
       idx_dst43_wetis = index_modal_map(85)
       idx_dst43_wetcw = index_modal_map(86)
       idx_dst53_wetis = index_modal_map(87)
       idx_dst53_wetcw = index_modal_map(88)
       idx_dst63_wetis = index_modal_map(89)
       idx_dst63_wetcw = index_modal_map(90)
       idx_dst73_wetis = index_modal_map(91)
       idx_dst73_wetcw = index_modal_map(92)
       idx_dst83_wetis = index_modal_map(93)
       idx_dst83_wetcw = index_modal_map(94)

       idx_dst12_wetis = index_modal_map(95)
       idx_dst12_wetcw = index_modal_map(96)
       idx_dst22_wetis = index_modal_map(97)
       idx_dst22_wetcw = index_modal_map(98)
       idx_dst32_wetis = index_modal_map(99)
       idx_dst32_wetcw = index_modal_map(100)
       idx_dst42_wetis = index_modal_map(101)
       idx_dst42_wetcw = index_modal_map(102)
       idx_dst52_wetis = index_modal_map(103)
       idx_dst52_wetcw = index_modal_map(104)
       idx_dst62_wetis = index_modal_map(105)
       idx_dst62_wetcw = index_modal_map(106)
       idx_dst72_wetis = index_modal_map(107)
       idx_dst72_wetcw = index_modal_map(108)
       idx_dst82_wetis = index_modal_map(109)
       idx_dst82_wetcw = index_modal_map(110)

       idx_fets1_dryis = index_modal_map(111)
       idx_fets1_drycw = index_modal_map(112)
       idx_fetm1_dryis = index_modal_map(113)
       idx_fetm1_drycw = index_modal_map(114)
       idx_fess1_dryis = index_modal_map(115)
       idx_fess1_drycw = index_modal_map(116)
       idx_fesm1_dryis = index_modal_map(117)
       idx_fesm1_drycw = index_modal_map(118)
       idx_fetc1_dryis = index_modal_map(119)
       idx_fetc1_drycw = index_modal_map(120)
       idx_fesc1_dryis = index_modal_map(121)
       idx_fesc1_drycw = index_modal_map(122)
       idx_fetf1_dryis = index_modal_map(123)
       idx_fetf1_drycw = index_modal_map(124)
       idx_fesf1_dryis = index_modal_map(125)
       idx_fesf1_drycw = index_modal_map(126)

       idx_fets3_dryis = index_modal_map(127)
       idx_fets3_drycw = index_modal_map(128)
       idx_fetm3_dryis = index_modal_map(129)
       idx_fetm3_drycw = index_modal_map(130)
       idx_fess3_dryis = index_modal_map(131)
       idx_fess3_drycw = index_modal_map(132)
       idx_fesm3_dryis = index_modal_map(133)
       idx_fesm3_drycw = index_modal_map(134)
       idx_fetc3_dryis = index_modal_map(135)
       idx_fetc3_drycw = index_modal_map(136)
       idx_fesc3_dryis = index_modal_map(137)
       idx_fesc3_drycw = index_modal_map(138)
       idx_fetf3_dryis = index_modal_map(139)
       idx_fetf3_drycw = index_modal_map(140)
       idx_fesf3_dryis = index_modal_map(141)
       idx_fesf3_drycw = index_modal_map(142)

       idx_fets2_dryis = index_modal_map(143)
       idx_fets2_drycw = index_modal_map(144)
       idx_fetm2_dryis = index_modal_map(145)
       idx_fetm2_drycw = index_modal_map(146)
       idx_fess2_dryis = index_modal_map(147)
       idx_fess2_drycw = index_modal_map(148)
       idx_fesm2_dryis = index_modal_map(149)
       idx_fesm2_drycw = index_modal_map(150)
       idx_fetc2_dryis = index_modal_map(151)
       idx_fetc2_drycw = index_modal_map(152)
       idx_fesc2_dryis = index_modal_map(153)
       idx_fesc2_drycw = index_modal_map(154)
       idx_fetf2_dryis = index_modal_map(155)
       idx_fetf2_drycw = index_modal_map(156)
       idx_fesf2_dryis = index_modal_map(157)
       idx_fesf2_drycw = index_modal_map(158)

       idx_fets1_wetis = index_modal_map(159)
       idx_fets1_wetcw = index_modal_map(160)
       idx_fetm1_wetis = index_modal_map(161)
       idx_fetm1_wetcw = index_modal_map(162)
       idx_fess1_wetis = index_modal_map(163)
       idx_fess1_wetcw = index_modal_map(164)
       idx_fesm1_wetis = index_modal_map(165)
       idx_fesm1_wetcw = index_modal_map(166)
       idx_fetc1_wetis = index_modal_map(167)
       idx_fetc1_wetcw = index_modal_map(168)
       idx_fesc1_wetis = index_modal_map(169)
       idx_fesc1_wetcw = index_modal_map(170)
       idx_fetf1_wetis = index_modal_map(171)
       idx_fetf1_wetcw = index_modal_map(172)
       idx_fesf1_wetis = index_modal_map(173)
       idx_fesf1_wetcw = index_modal_map(174)

       idx_fets3_wetis = index_modal_map(175)
       idx_fets3_wetcw = index_modal_map(176)
       idx_fetm3_wetis = index_modal_map(177)
       idx_fetm3_wetcw = index_modal_map(178)
       idx_fess3_wetis = index_modal_map(179)
       idx_fess3_wetcw = index_modal_map(180)
       idx_fesm3_wetis = index_modal_map(181)
       idx_fesm3_wetcw = index_modal_map(182)
       idx_fetc3_wetis = index_modal_map(183)
       idx_fetc3_wetcw = index_modal_map(184)
       idx_fesc3_wetis = index_modal_map(185)
       idx_fesc3_wetcw = index_modal_map(186)
       idx_fetf3_wetis = index_modal_map(187)
       idx_fetf3_wetcw = index_modal_map(188)
       idx_fesf3_wetis = index_modal_map(189)
       idx_fesf3_wetcw = index_modal_map(190)

       idx_fets2_wetis = index_modal_map(191)
       idx_fets2_wetcw = index_modal_map(192)
       idx_fetm2_wetis = index_modal_map(193)
       idx_fetm2_wetcw = index_modal_map(194)
       idx_fess2_wetis = index_modal_map(195)
       idx_fess2_wetcw = index_modal_map(196)
       idx_fesm2_wetis = index_modal_map(197)
       idx_fesm2_wetcw = index_modal_map(198)
       idx_fetc2_wetis = index_modal_map(199)
       idx_fetc2_wetcw = index_modal_map(200)
       idx_fesc2_wetis = index_modal_map(201)
       idx_fesc2_wetcw = index_modal_map(202)
       idx_fetf2_wetis = index_modal_map(203)
       idx_fetf2_wetcw = index_modal_map(204)
       idx_fesf2_wetis = index_modal_map(205)
       idx_fesf2_wetcw = index_modal_map(206)

       call modal_aero_deposition_init( &
                 bcphi_indices=(/idx_bc1/), &
                 ocphi_indices=(/idx_pom1,idx_soa1/), &
                 ocpho_indices=(/idx_soa2/), &
                 ufine_dust_indices=(/idx_dst12,idx_dst22,idx_dst32,idx_dst42,idx_dst52,idx_dst62,idx_dst72,idx_dst82,&
                                      idx_fets2,idx_fetm2,idx_fess2,idx_fesm2,idx_fetc2,idx_fesc2,idx_fetf2,idx_fesf2/),&
                 fine_dust_indices=(/idx_dst11,idx_dst21,idx_dst31,idx_dst41,idx_dst51,idx_dst61,idx_dst71,idx_dst81,&
                                     idx_fets1,idx_fetm1,idx_fess1,idx_fesm1,idx_fetc1,idx_fesc1,idx_fetf1,idx_fesf1/),&
                 crse_dust_indices=(/idx_dst31,idx_dst23,idx_dst33,idx_dst43,idx_dst53,idx_dst63,idx_dst73,idx_dst83,&
                                     idx_fets3,idx_fetm3,idx_fess3,idx_fesm3,idx_fetc3,idx_fesc3,idx_fetf3,idx_fesf3/) )
#else
       idx_bc1_dryis  = index_modal_map(1)
       idx_bc1_drycw  = index_modal_map(2)
       idx_pom1_dryis = index_modal_map(3)
       idx_pom1_drycw = index_modal_map(4)
       idx_soa1_dryis = index_modal_map(5)
       idx_soa1_drycw = index_modal_map(6)
       idx_soa2_dryis = index_modal_map(7)
       idx_soa2_drycw = index_modal_map(8)
       idx_dst11_dryis = index_modal_map(9)
       idx_dst11_drycw = index_modal_map(10)
       idx_dst21_dryis = index_modal_map(11)
       idx_dst21_drycw = index_modal_map(12)
       idx_dst31_dryis = index_modal_map(13)
       idx_dst31_drycw = index_modal_map(14)
       idx_dst41_dryis = index_modal_map(15)
       idx_dst41_drycw = index_modal_map(16)
       idx_dst51_dryis = index_modal_map(17)
       idx_dst51_drycw = index_modal_map(18)
       idx_dst61_dryis = index_modal_map(19)
       idx_dst61_drycw = index_modal_map(20)
       idx_dst71_dryis = index_modal_map(21)
       idx_dst71_drycw = index_modal_map(22)
       idx_dst81_dryis = index_modal_map(23)
       idx_dst81_drycw = index_modal_map(24)

       idx_dst13_dryis = index_modal_map(25)
       idx_dst13_drycw = index_modal_map(26)
       idx_dst23_dryis = index_modal_map(27)
       idx_dst23_drycw = index_modal_map(28)
       idx_dst33_dryis = index_modal_map(29)
       idx_dst33_drycw = index_modal_map(30)
       idx_dst43_dryis = index_modal_map(31)
       idx_dst43_drycw = index_modal_map(32)
       idx_dst53_dryis = index_modal_map(33)
       idx_dst53_drycw = index_modal_map(34)
       idx_dst63_dryis = index_modal_map(35)
       idx_dst63_drycw = index_modal_map(36)
       idx_dst73_dryis = index_modal_map(37)
       idx_dst73_drycw = index_modal_map(38)
       idx_dst83_dryis = index_modal_map(39)
       idx_dst83_drycw = index_modal_map(40)

       idx_bc1_wetis  = index_modal_map(41)
       idx_bc1_wetcw  = index_modal_map(42)
       idx_pom1_wetis = index_modal_map(43)
       idx_pom1_wetcw = index_modal_map(44)
       idx_soa1_wetis = index_modal_map(45)
       idx_soa1_wetcw = index_modal_map(46)
       idx_dst11_wetis = index_modal_map(47)
       idx_dst11_wetcw = index_modal_map(48)
       idx_dst21_wetis = index_modal_map(49)
       idx_dst21_wetcw = index_modal_map(50)
       idx_dst31_wetis = index_modal_map(51)
       idx_dst31_wetcw = index_modal_map(52)
       idx_dst41_wetis = index_modal_map(53)
       idx_dst41_wetcw = index_modal_map(54)
       idx_dst51_wetis = index_modal_map(55)
       idx_dst51_wetcw = index_modal_map(56)
       idx_dst61_wetis = index_modal_map(57)
       idx_dst61_wetcw = index_modal_map(58)
       idx_dst71_wetis = index_modal_map(59)
       idx_dst71_wetcw = index_modal_map(60)
       idx_dst81_wetis = index_modal_map(61)
       idx_dst81_wetcw = index_modal_map(62)

       idx_dst13_wetis = index_modal_map(63)
       idx_dst13_wetcw = index_modal_map(64)
       idx_dst23_wetis = index_modal_map(65)
       idx_dst23_wetcw = index_modal_map(66)
       idx_dst33_wetis = index_modal_map(67)
       idx_dst33_wetcw = index_modal_map(68)
       idx_dst43_wetis = index_modal_map(69)
       idx_dst43_wetcw = index_modal_map(70)
       idx_dst53_wetis = index_modal_map(71)
       idx_dst53_wetcw = index_modal_map(72)
       idx_dst63_wetis = index_modal_map(73)
       idx_dst63_wetcw = index_modal_map(74)
       idx_dst73_wetis = index_modal_map(75)
       idx_dst73_wetcw = index_modal_map(76)
       idx_dst83_wetis = index_modal_map(77)
       idx_dst83_wetcw = index_modal_map(78)

       call modal_aero_deposition_init( &
                 bcphi_indices=(/idx_bc1/), &
                 ocphi_indices=(/idx_pom1,idx_soa1/), &
                 ocpho_indices=(/idx_soa2/), &
                 fine_dust_indices=(/idx_dst11,idx_dst21,idx_dst31,idx_dst41,idx_dst51,idx_dst61,idx_dst71,idx_dst81/),&
                 crse_dust_indices=(/idx_dst31,idx_dst23,idx_dst33,idx_dst43,idx_dst53,idx_dst63,idx_dst73,idx_dst83/) )
#endif
#else 
       idx_bc1_dryis  = index_modal_map(1)
       idx_bc1_drycw  = index_modal_map(2)
       idx_pom1_dryis = index_modal_map(3)
       idx_pom1_drycw = index_modal_map(4)
       idx_soa1_dryis = index_modal_map(5)
       idx_soa1_drycw = index_modal_map(6)
       idx_soa2_dryis = index_modal_map(7)
       idx_soa2_drycw = index_modal_map(8)
       idx_dst1_dryis = index_modal_map(9)
       idx_dst1_drycw = index_modal_map(10)
       idx_dst3_dryis = index_modal_map(11)
       idx_dst3_drycw = index_modal_map(12)

       idx_bc1_wetis  = index_modal_map(13)
       idx_bc1_wetcw  = index_modal_map(14)
       idx_pom1_wetis = index_modal_map(15)
       idx_pom1_wetcw = index_modal_map(16)
       idx_soa1_wetis = index_modal_map(17)
       idx_soa1_wetcw = index_modal_map(18)
       idx_dst1_wetis = index_modal_map(19)
       idx_dst1_wetcw = index_modal_map(20)
       idx_dst3_wetis = index_modal_map(21)
       idx_dst3_wetcw = index_modal_map(22)

       call modal_aero_deposition_init( bcphi_indices=(/idx_bc1/), &
                                        ocphi_indices=(/idx_pom1,idx_soa1/), &
                                        ocpho_indices=(/idx_soa2/), &
                                        fine_dust_indices=(/idx_dst1/),&
                                        crse_dust_indices=(/idx_dst3/) )
#endif
    else

       ibcphiwet = index_bulk_map(1)
       ibcphodry = index_bulk_map(2)
       ibcphidry = index_bulk_map(3)
       iocphiwet = index_bulk_map(4)
       iocphodry = index_bulk_map(5)
       iocphidry = index_bulk_map(6)
       idstdry1  = index_bulk_map(7)
       idstdry2  = index_bulk_map(8)
       idstdry3  = index_bulk_map(9)
       idstdry4  = index_bulk_map(10)
       idstwet1  = index_bulk_map(11)
       idstwet2  = index_bulk_map(12)
       idstwet3  = index_bulk_map(13)
       idstwet4  = index_bulk_map(14)

    endif

  end subroutine aerodep_flx_init
!-------------------------------------------------------------------
! sets namelist options
!-------------------------------------------------------------------
subroutine aerodep_flx_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'aerodep_flx_readnl'

   character(len=32)  :: aerodep_flx_specifier(N_MODAL)
   character(len=256) :: aerodep_flx_file
   character(len=256) :: aerodep_flx_filelist
   character(len=256) :: aerodep_flx_datapath
   character(len=32)  :: aerodep_flx_type
   logical            :: aerodep_flx_rmfile
   integer            :: aerodep_flx_cycle_yr
   integer            :: aerodep_flx_fixed_ymd
   integer            :: aerodep_flx_fixed_tod

   namelist /aerodep_flx_nl/ &
      aerodep_flx_specifier, &
      aerodep_flx_file,      &
      aerodep_flx_filelist,  &
      aerodep_flx_datapath,  &
      aerodep_flx_type,      &
      aerodep_flx_rmfile,    &
      aerodep_flx_cycle_yr,  &
      aerodep_flx_fixed_ymd, &
      aerodep_flx_fixed_tod      
   !-----------------------------------------------------------------------------

   ! Initialize namelist variables from local module variables.
   aerodep_flx_specifier= specifier
   aerodep_flx_file     = filename
   aerodep_flx_filelist = filelist
   aerodep_flx_datapath = datapath
   aerodep_flx_type     = datatype
   aerodep_flx_rmfile   = rmv_file
   aerodep_flx_cycle_yr = cycle_yr
   aerodep_flx_fixed_ymd= fixed_ymd
   aerodep_flx_fixed_tod= fixed_tod

   ! Read namelist
   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'aerodep_flx_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, aerodep_flx_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(aerodep_flx_specifier,len(aerodep_flx_specifier(1))*N_MODAL,     mpichar, 0, mpicom)
   call mpibcast(aerodep_flx_file,     len(aerodep_flx_file),     mpichar, 0, mpicom)
   call mpibcast(aerodep_flx_filelist, len(aerodep_flx_filelist), mpichar, 0, mpicom)
   call mpibcast(aerodep_flx_datapath, len(aerodep_flx_datapath), mpichar, 0, mpicom)
   call mpibcast(aerodep_flx_type,     len(aerodep_flx_type),     mpichar, 0, mpicom)
   call mpibcast(aerodep_flx_rmfile,   1, mpilog,  0, mpicom)
   call mpibcast(aerodep_flx_cycle_yr, 1, mpiint,  0, mpicom)
   call mpibcast(aerodep_flx_fixed_ymd,1, mpiint,  0, mpicom)
   call mpibcast(aerodep_flx_fixed_tod,1, mpiint,  0, mpicom)
#endif

   ! Update module variables with user settings.
   specifier  = aerodep_flx_specifier
   filename   = aerodep_flx_file
   filelist   = aerodep_flx_filelist
   datapath   = aerodep_flx_datapath
   datatype   = aerodep_flx_type
   rmv_file   = aerodep_flx_rmfile
   cycle_yr   = aerodep_flx_cycle_yr
   fixed_ymd  = aerodep_flx_fixed_ymd
   fixed_tod  = aerodep_flx_fixed_tod

   ! Turn on prescribed volcanics if user has specified an input dataset.
   if (len_trim(filename) > 0 .and. filename.ne.'NONE' ) has_aerodep_flx = .true.

end subroutine aerodep_flx_readnl

!-------------------------------------------------------------------
! sets the aerosol deposition fluxes in the cam_out structure 
! to be sent to the surface models
!-------------------------------------------------------------------
  subroutine aerodep_flx_set( cam_out, ncol, lchnk )
    use camsrfexch,       only : cam_out_t     

    type(cam_out_t),     intent(inout) :: cam_out
    integer,             intent(in)    :: ncol, lchnk
    
    if( .not. has_aerodep_flx ) return
    
    if (modal_fluxes) then
       call set_modal_fluxes( cam_out, ncol, lchnk )
    else
       call set_bulk_fluxes( cam_out, ncol, lchnk )
    endif

  end subroutine aerodep_flx_set

!-------------------------------------------------------------------
! advances the prescribed fluxes to the current time step
!-------------------------------------------------------------------
  subroutine aerodep_flx_adv( state, pbuf2d, cam_out )

    use tracer_data,      only : advance_trcdata
    use physics_types,    only : physics_state
    use camsrfexch,       only : cam_out_t
    use physics_buffer, only : physics_buffer_desc

    implicit none

    type(physics_state), intent(in)    :: state(begchunk:endchunk)                 
    type(cam_out_t),     intent(inout) :: cam_out(begchunk:endchunk)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    integer :: c, ncol
    
    if( .not. has_aerodep_flx ) return

    call advance_trcdata( fields, file, state, pbuf2d  )

!$OMP PARALLEL DO PRIVATE (C, NCOL)
    do c = begchunk, endchunk
       ncol = state(c)%ncol
       call aerodep_flx_set( cam_out(c), ncol, c )
    enddo

  end subroutine aerodep_flx_adv

!-------------------------------------------------------------------
! returns true if aerosol dep fluxes are prescribed from dataset
!-------------------------------------------------------------------
  function aerodep_flx_prescribed()
    logical :: aerodep_flx_prescribed
    aerodep_flx_prescribed = has_aerodep_flx
  endfunction aerodep_flx_prescribed

! private methods
!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine set_bulk_fluxes( cam_out, ncol, lchnk )
    use camsrfexch,            only : cam_out_t     

    ! Arguments
    type(cam_out_t), intent(inout) :: cam_out
    integer,         intent(in)    :: ncol, lchnk

    call set_fluxes( cam_out%bcphiwet, ibcphiwet, ncol, lchnk )
    call set_fluxes( cam_out%bcphidry, ibcphidry, ncol, lchnk )
    call set_fluxes( cam_out%bcphodry, ibcphodry, ncol, lchnk )

    call set_fluxes( cam_out%ocphiwet, iocphiwet, ncol, lchnk )
    call set_fluxes( cam_out%ocphidry, iocphidry, ncol, lchnk )
    call set_fluxes( cam_out%ocphodry, iocphodry, ncol, lchnk )

    call set_fluxes( cam_out%dstdry1, idstdry1, ncol, lchnk )
    call set_fluxes( cam_out%dstdry2, idstdry2, ncol, lchnk )
    call set_fluxes( cam_out%dstdry3, idstdry3, ncol, lchnk )
    call set_fluxes( cam_out%dstdry4, idstdry4, ncol, lchnk )

    call set_fluxes( cam_out%dstwet1, idstwet1, ncol, lchnk )
    call set_fluxes( cam_out%dstwet2, idstwet2, ncol, lchnk )
    call set_fluxes( cam_out%dstwet3, idstwet3, ncol, lchnk )
    call set_fluxes( cam_out%dstwet4, idstwet4, ncol, lchnk )

  end subroutine set_bulk_fluxes

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine set_modal_fluxes( cam_out, ncol, lchnk )
    use camsrfexch,            only : cam_out_t     
    use modal_aero_deposition, only : set_srf_drydep, set_srf_wetdep

    ! Arguments
    type(cam_out_t), intent(inout) :: cam_out
    integer,         intent(in)    :: ncol, lchnk

    ! local vars
    integer :: i
    real(r8) :: aerdepdryis(pcols,nmodal_idxs)
    real(r8) :: aerdepdrycw(pcols,nmodal_idxs)
    real(r8) :: aerdepwetis(pcols,nmodal_idxs)
    real(r8) :: aerdepwetcw(pcols,nmodal_idxs)

    ! bin the fluxes as using modal_aero_deposition...

    aerdepdryis(:,:) = 0._r8
    aerdepdrycw(:,:) = 0._r8
    aerdepwetis(:,:) = 0._r8
    aerdepwetcw(:,:) = 0._r8

    call set_fluxes( aerdepwetis(:ncol,idx_bc1 ), idx_bc1_wetis , ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_bc1 ), idx_bc1_wetcw , ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_pom1), idx_pom1_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_pom1), idx_pom1_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_soa1), idx_soa1_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_soa1), idx_soa1_wetcw, ncol, lchnk )
#if ( defined DUSTTYPES )
    call set_fluxes( aerdepwetis(:ncol,idx_dst11), idx_dst11_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst11), idx_dst11_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst21), idx_dst21_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst21), idx_dst21_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst31), idx_dst31_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst31), idx_dst31_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst41), idx_dst41_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst41), idx_dst41_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst51), idx_dst51_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst51), idx_dst51_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst61), idx_dst61_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst61), idx_dst61_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst71), idx_dst71_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst71), idx_dst71_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst81), idx_dst81_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst81), idx_dst81_wetcw, ncol, lchnk )

    call set_fluxes( aerdepwetis(:ncol,idx_dst13), idx_dst13_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst13), idx_dst13_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst23), idx_dst23_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst23), idx_dst23_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst33), idx_dst33_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst33), idx_dst33_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst43), idx_dst43_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst43), idx_dst43_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst53), idx_dst53_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst53), idx_dst53_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst63), idx_dst63_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst63), idx_dst63_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst73), idx_dst73_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst73), idx_dst73_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst83), idx_dst83_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst83), idx_dst83_wetcw, ncol, lchnk )

#if ( defined IRONPROC )
    call set_fluxes( aerdepwetis(:ncol,idx_fets1), idx_fets1_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_fets1), idx_fets1_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_fetm1), idx_fetm1_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_fetm1), idx_fetm1_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_fess1), idx_fess1_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_fess1), idx_fess1_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_fesm1), idx_fesm1_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_fesm1), idx_fesm1_wetcw, ncol, lchnk )

    call set_fluxes( aerdepwetis(:ncol,idx_fets3), idx_fets3_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_fets3), idx_fets3_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_fetm3), idx_fetm3_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_fetm3), idx_fetm3_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_fess3), idx_fess3_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_fess3), idx_fess3_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_fesm3), idx_fesm3_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_fesm3), idx_fesm3_wetcw, ncol, lchnk )
#endif
#else   
    call set_fluxes( aerdepwetis(:ncol,idx_dst1), idx_dst1_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst1), idx_dst1_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst3), idx_dst3_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst3), idx_dst3_wetcw, ncol, lchnk )
#endif

    call set_fluxes( aerdepdryis(:ncol,idx_bc1 ), idx_bc1_dryis , ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_bc1 ), idx_bc1_drycw , ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_pom1), idx_pom1_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_pom1), idx_pom1_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_soa1), idx_soa1_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_soa1), idx_soa1_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_soa2), idx_soa2_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_soa2), idx_soa2_drycw, ncol, lchnk )
#if ( defined DUSTTYPES )
    call set_fluxes( aerdepdryis(:ncol,idx_dst11), idx_dst11_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst11), idx_dst11_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst21), idx_dst21_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst21), idx_dst21_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst31), idx_dst31_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst31), idx_dst31_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst41), idx_dst41_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst41), idx_dst41_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst51), idx_dst51_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst51), idx_dst51_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst61), idx_dst61_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst61), idx_dst61_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst71), idx_dst71_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst71), idx_dst71_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst81), idx_dst81_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst81), idx_dst81_drycw, ncol, lchnk )

    call set_fluxes( aerdepdryis(:ncol,idx_dst13), idx_dst13_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst13), idx_dst13_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst23), idx_dst23_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst23), idx_dst23_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst33), idx_dst33_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst33), idx_dst33_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst43), idx_dst43_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst43), idx_dst43_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst53), idx_dst53_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst53), idx_dst53_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst63), idx_dst63_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst63), idx_dst63_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst73), idx_dst73_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst73), idx_dst73_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst83), idx_dst83_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst83), idx_dst83_drycw, ncol, lchnk )

#if ( defined IRONPROC )
    call set_fluxes( aerdepdryis(:ncol,idx_fets1), idx_fets1_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_fets1), idx_fets1_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_fetm1), idx_fetm1_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_fetm1), idx_fetm1_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_fess1), idx_fess1_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_fess1), idx_fess1_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_fesm1), idx_fesm1_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_fesm1), idx_fesm1_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_fetc1), idx_fetc1_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_fetc1), idx_fetc1_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_fesc1), idx_fesc1_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_fesc1), idx_fesc1_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_fetf1), idx_fetf1_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_fetf1), idx_fetf1_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_fesf1), idx_fesf1_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_fesf1), idx_fesf1_drycw, ncol, lchnk )

    call set_fluxes( aerdepdryis(:ncol,idx_fets3), idx_fets3_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_fets3), idx_fets3_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_fetm3), idx_fetm3_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_fetm3), idx_fetm3_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_fess3), idx_fess3_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_fess3), idx_fess3_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_fesm3), idx_fesm3_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_fesm3), idx_fesm3_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_fetc3), idx_fetc3_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_fetc3), idx_fetc3_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_fesc3), idx_fesc3_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_fesc3), idx_fesc3_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_fetf3), idx_fetf3_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_fetf3), idx_fetf3_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_fesf3), idx_fesf3_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_fesf3), idx_fesf3_drycw, ncol, lchnk )
#endif
#else   
    call set_fluxes( aerdepdryis(:ncol,idx_dst1), idx_dst1_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst1), idx_dst1_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst3), idx_dst3_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst3), idx_dst3_drycw, ncol, lchnk )
#endif

#if (defined IRONPROC)
    call set_srf_drydep(lchnk, aerdepdryis, aerdepdrycw, cam_out)
    call set_srf_wetdep(lchnk, aerdepwetis, aerdepwetcw, cam_out)
#else
    call set_srf_drydep(aerdepdryis, aerdepdrycw, cam_out)
    call set_srf_wetdep(aerdepwetis, aerdepwetcw, cam_out)
#endif
  end  subroutine set_modal_fluxes

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine set_fluxes( fluxes, fld_indx, ncol, lchnk )
    use cam_history,  only : outfld

    real(r8), intent(inout) :: fluxes(:)
    integer,  intent(in)    :: fld_indx, ncol, lchnk

    integer :: i

    if (fld_indx<1) return

    do i = 1,ncol
       ! modal aero wet dep history fields are negative
       fluxes(i) = fields(fld_indx)%data(i,1,lchnk)
    enddo

    call outfld(trim(fields(fld_indx)%fldnam)//'_D', fluxes(:ncol), ncol, lchnk )

  endsubroutine set_fluxes

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  integer function get_ndx( name, list )

    implicit none
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: list(:)

    integer :: i
    integer :: maxnum

    maxnum = size(list)

    get_ndx = -1
    do i = 1, maxnum
      if ( trim(name) == trim(list(i)) ) then
        get_ndx = i
        return
      endif
    enddo

  end function get_ndx

end module aerodep_flx
