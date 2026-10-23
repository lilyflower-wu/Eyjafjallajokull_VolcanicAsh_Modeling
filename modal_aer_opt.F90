module modal_aer_opt

! parameterizes aerosol coefficients using chebychev polynomial
! parameterize aerosol radiative properties in terms of
! surface mode wet radius and wet refractive index

! Ghan and Zaveri, JGR 2007.

! uses Wiscombe's (1979) mie scattering code


use shr_kind_mod,      only: r8 => shr_kind_r8, shr_kind_cl
use ppgrid,            only: pcols, pver, pverp
use constituents,      only: pcnst
use spmd_utils,        only: masterproc
use ref_pres,          only: top_lev => clim_modal_aero_top_lev
use physconst,         only: rhoh2o, rga, rair
use radconstants,      only: nswbands, nlwbands, idx_sw_diag, idx_uv_diag, idx_nir_diag
use rad_constituents,  only: n_diag, rad_cnst_get_call_list, rad_cnst_get_info, rad_cnst_get_aer_mmr, &
                             rad_cnst_get_aer_props, rad_cnst_get_mode_props
use physics_types,     only: physics_state

use physics_buffer, only : pbuf_get_index,physics_buffer_desc, pbuf_get_field
use pio,               only: file_desc_t, var_desc_t, pio_inq_dimlen, pio_inq_dimid, pio_inq_varid, &
                             pio_get_var, pio_nowrite, pio_closefile
use cam_pio_utils,     only: cam_pio_openfile
use cam_history,       only: addfld, add_default, outfld, horiz_only
use cam_history_support, only: fillvalue
use cam_logfile,       only: iulog
use perf_mod,          only: t_startf, t_stopf
use cam_abortutils,    only: endrun

use modal_aero_wateruptake, only: modal_aero_wateruptake_dr
use modal_aero_calcsize,    only: modal_aero_calcsize_diag

implicit none
private
save

public :: modal_aer_opt_readnl, modal_aer_opt_init, modal_aero_sw, modal_aero_lw


character(len=*), parameter :: unset_str = 'UNSET'

! Namelist variables:
character(shr_kind_cl)      :: modal_optics_file = unset_str   ! full pathname for modal optics dataset
character(shr_kind_cl)      :: water_refindex_file = unset_str ! full pathname for water refractive index dataset

! Dimension sizes in coefficient arrays used to parameterize aerosol radiative properties
! in terms of refractive index and wet radius
integer, parameter :: ncoef=5, prefr=7, prefi=10

real(r8) :: xrmin, xrmax

! refractive index for water read in read_water_refindex
complex(r8) :: crefwsw(nswbands) ! complex refractive index for water visible
complex(r8) :: crefwlw(nlwbands) ! complex refractive index for water infrared

! physics buffer indices
integer :: dgnumwet_idx = -1
integer :: qaerwat_idx  = -1

character(len=4) :: diag(0:n_diag) = (/'    ','_d1 ','_d2 ','_d3 ','_d4 ','_d5 ', &
                                       '_d6 ','_d7 ','_d8 ','_d9 ','_d10'/)

!===============================================================================
CONTAINS
!===============================================================================

subroutine modal_aer_opt_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'modal_aer_opt_readnl'

   namelist /modal_aer_opt_nl/ water_refindex_file
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'modal_aer_opt_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, modal_aer_opt_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   call mpibcast(water_refindex_file, len(water_refindex_file), mpichar, 0, mpicom)
#endif


end subroutine modal_aer_opt_readnl

!===============================================================================

subroutine modal_aer_opt_init()

   use ioFileMod,        only: getfil
   use phys_control,     only: phys_getopts

   ! Local variables

   integer  :: i, m
   real(r8) :: rmmin, rmmax       ! min, max aerosol surface mode radius treated (m)
   character(len=256) :: locfile
   
   logical           :: history_amwg            ! output the variables used by the AMWG diag package
   logical           :: history_aero_optics     ! output aerosol optics diagnostics
   logical           :: history_dust            ! output dust diagnostics

   logical :: call_list(0:n_diag)
   integer :: ilist, nmodes, m_ncoef, m_prefr, m_prefi
   integer :: errcode

   character(len=*), parameter :: routine='modal_aer_opt_init'
   character(len=10) :: fldname
   character(len=128) :: lngname

   !----------------------------------------------------------------------------

   rmmin = 0.01e-6_r8
   rmmax = 25.e-6_r8
   xrmin = log(rmmin)
   xrmax = log(rmmax)

   ! Check that dimension sizes in the coefficient arrays used to
   ! parameterize aerosol radiative properties are consistent between this
   ! module and the mode physprop files.
   call rad_cnst_get_call_list(call_list)
   do ilist = 0, n_diag
      if (call_list(ilist)) then
         call rad_cnst_get_info(ilist, nmodes=nmodes)
         do m = 1, nmodes
            call rad_cnst_get_mode_props(ilist, m, ncoef=m_ncoef, prefr=m_prefr, prefi=m_prefi)
            if (m_ncoef /= ncoef .or. m_prefr /= prefr .or. m_prefi /= prefi) then
               write(iulog,*) routine//': ERROR - file and module values do not match:'
               write(iulog,*) '   ncoef:', ncoef, m_ncoef
               write(iulog,*) '   prefr:', prefr, m_prefr
               write(iulog,*) '   prefi:', prefi, m_prefi
               call endrun(routine//': ERROR - file and module values do not match')
            end if
         end do
      end if
   end do

   ! Initialize physics buffer indices for dgnumwet and qaerwat.  Note the implicit assumption
   ! that the loops over modes in the optics calculations will use the values for dgnumwet and qaerwat
   ! that are set in the aerosol_wet_intr code.
   dgnumwet_idx = pbuf_get_index('DGNUMWET',errcode)
   if (errcode < 0) then
      call endrun(routine//' ERROR: cannot find physics buffer field DGNUMWET')
   end if
   qaerwat_idx  = pbuf_get_index('QAERWAT',errcode)
   if (errcode < 0) then
      call endrun(routine//' ERROR: cannot find physics buffer field QAERWAT')
   end if

   call getfil(water_refindex_file, locfile)
   call read_water_refindex(locfile)
   if (masterproc) write(iulog,*) "modal_aer_opt_init: read water refractive index file:", trim(locfile)

   call phys_getopts(history_amwg_out        = history_amwg, &
                     history_aero_optics_out = history_aero_optics, &
                     history_dust_out        = history_dust )

   ! Add diagnostic fields to history output.

   call addfld ('EXTINCT',    (/ 'lev' /), 'A','/m','Aerosol extinction 550 nm, day only',                   &
                flag_xyfill=.true.)
   !++ Longlei Li
   call addfld ('EXTDUST',    (/ 'lev' /), 'A','/m','Dust extinction 550 nm, day only',                   &
                flag_xyfill=.true.)
   !-- Longlei Li

   call addfld ('EXTINCTUV',  (/ 'lev' /), 'A','/m','Aerosol extinction 350 nm, day only',                   &
                flag_xyfill=.true.)
   call addfld ('EXTINCTNIR', (/ 'lev' /), 'A','/m','Aerosol extinction 1020 nm, day only',                  &
                flag_xyfill=.true.)
   call addfld ('ABSORB',     (/ 'lev' /), 'A','/m','Aerosol absorption, day only',                          &
                flag_xyfill=.true.)
   call addfld ('AODVIS',     horiz_only,  'A','  ','Aerosol optical depth 550 nm, day only',                &
                flag_xyfill=.true.)
   call addfld ('AODVISst',   horiz_only,  'A','  ','Stratospheric aerosol optical depth 550 nm, day only',  &
                flag_xyfill=.true.)
   call addfld ('AODUV',      horiz_only,  'A','  ','Aerosol optical depth 350 nm, day only',                &
                flag_xyfill=.true.)
   call addfld ('AODUVst',    horiz_only,  'A','  ','Stratospheric aerosol optical depth 350 nm, day only',  &
                flag_xyfill=.true.)
   call addfld ('AODNIR',     horiz_only,  'A','  ','Aerosol optical depth 1020 nm, day only',               &
                flag_xyfill=.true.)
   call addfld ('AODNIRst',   horiz_only,  'A','  ','Stratospheric aerosol optical depth 1020 nm, day only', &
                flag_xyfill=.true.)
   call addfld ('AODABS',     horiz_only,  'A','  ','Aerosol absorption optical depth 550 nm, day only',     &
                flag_xyfill=.true.)
   call addfld ('AODxASYM',   horiz_only,  'A','  ','Aerosol optical depth 550 * asymmetry factor, day only',&
                flag_xyfill=.true.)
   call addfld ('EXTxASYM',   (/ 'lev' /), 'A','  ','extinction 550 nm * asymmetry factor, day only',        &
                flag_xyfill=.true.)

   call addfld ('EXTINCTdn',    (/ 'lev' /), 'A','/m','Aerosol extinction 550 nm, day night',                 &
                  flag_xyfill=.true.)

   !++ Longlei Li
   call addfld ('EXTDUSTdn',    (/ 'lev' /), 'A','/m','Dust extinction 550 nm, day night',                 &
                  flag_xyfill=.true.)
   !-- Longlei Li

   call addfld ('EXTINCTUVdn',  (/ 'lev' /), 'A','/m','Aerosol extinction 350 nm, day night',                 &
                  flag_xyfill=.true.)
   call addfld ('EXTINCTNIRdn', (/ 'lev' /), 'A','/m','Aerosol extinction 1020 nm, day night',                &
                  flag_xyfill=.true.)
   call addfld ('ABSORBdn',     (/ 'lev' /), 'A','/m','Aerosol absorption, day night',                        &
                  flag_xyfill=.true.)
   call addfld ('AODVISdn',     horiz_only,  'A','  ','Aerosol optical depth 550 nm, day night',              &
                  flag_xyfill=.true.)
   call addfld ('AODVISstdn',   horiz_only,  'A','  ','Stratospheric aerosol optical depth 550 nm, day night',&
                  flag_xyfill=.true.)
   call addfld ('AODUVdn',      horiz_only,  'A','  ','Aerosol optical depth 350 nm, day night',              &
                  flag_xyfill=.true.)
   call addfld ('AODUVstdn',    horiz_only,  'A','  ','Stratospheric aerosol optical depth 350 nm, day night',&
                  flag_xyfill=.true.)
   call addfld ('AODNIRdn',     horiz_only,  'A','  ','Aerosol optical depth 1020 nm, day night',             &
                  flag_xyfill=.true.)
   call addfld ('AODNIRstdn',   horiz_only,  'A','  ','Stratospheric aerosol optical depth 1020 nm, day night',&
                 flag_xyfill=.true.)
   call addfld ('AODABSdn',     horiz_only,  'A','  ','Aerosol absorption optical depth 550 nm, day night',   &
                  flag_xyfill=.true.)
   call addfld ('AODxASYMdn',   horiz_only,  'A','  ','Aerosol optical depth 550 * asymmetry factor, day night',&
                flag_xyfill=.true.)
   call addfld ('EXTxASYMdn',   (/ 'lev' /), 'A','  ','extinction 550 * asymmetry factor, day night',        &
                flag_xyfill=.true.)

   call rad_cnst_get_info(0, nmodes=nmodes)

   do m = 1, nmodes

      write(fldname,'(a,i1)') 'BURDEN', m
      write(lngname,'(a,i1)') 'Aerosol burden, day only, mode ', m
      call addfld (fldname, horiz_only, 'A', 'kg/m2', lngname, flag_xyfill=.true.)
      if (m>3 .and. history_aero_optics) then
         call add_default (fldname, 1, ' ')
      endif

      write(fldname,'(a,i1)') 'AODMODE', m
      write(lngname,'(a,i1)') 'Aerosol optical depth, day only, 550 nm mode ', m
      call addfld (fldname, horiz_only, 'A', '  ', lngname, flag_xyfill=.true.)
      if (m>3 .and. history_aero_optics) then
         call add_default (fldname, 1, ' ')
      endif

      write(fldname,'(a,i1)') 'AODDUST', m
      write(lngname,'(a,i1,a)') 'Aerosol optical depth, day only, 550 nm mode ',m,' from dust'
      call addfld (fldname, horiz_only, 'A', '  ', lngname, flag_xyfill=.true.)
      if (m>3 .and. history_aero_optics) then
         call add_default (fldname, 1, ' ')
      endif

      !++ Longlei Li ++
      write(fldname,'(a,i1)') 'AODASH', m
      write(lngname,'(a,i1,a)') 'Aerosol optical depth, day only, 550 nm mode ',m,' from ash'
      call addfld (fldname, horiz_only, 'A', '  ', lngname, flag_xyfill=.true.)
      if (m>3 .and. history_aero_optics) then
         call add_default (fldname, 1, ' ')
      endif
      !-- Longlei Li --

      write(fldname,'(a,i1)') 'BURDENdn', m
      write(lngname,'(a,i1)') 'Aerosol burden, day night, mode ', m
      call addfld (fldname, horiz_only, 'A', 'kg/m2', lngname, flag_xyfill=.true.)
      if (m>3 .and. history_aero_optics) then
         call add_default (fldname, 1, ' ')
      endif

      write(fldname,'(a,i1)') 'AODdnMODE', m
      write(lngname,'(a,i1)') 'Aerosol optical depth 550 nm, day night, mode ', m
      call addfld (fldname, horiz_only, 'A', '  ', lngname, flag_xyfill=.true.)
      if (m>3 .and. history_aero_optics) then
         call add_default (fldname, 1, ' ')
      endif

      write(fldname,'(a,i1)') 'AODdnDUST', m
      write(lngname,'(a,i1,a)') 'Aerosol optical depth 550 nm, day night, mode ',m,' from dust'
      call addfld (fldname, horiz_only, 'A', '  ', lngname, flag_xyfill=.true.)
      if (m>3 .and. history_aero_optics) then
         call add_default (fldname, 1, ' ')
      endif

      !++ Longlei Li
      write(fldname,'(a,i1)') 'EXTDUST', m
      write(lngname,'(a,i1,a)') 'Dust extinction 550 nm, day only, mode ',m,' from dust'
      call addfld (fldname, horiz_only, 'A', '  ', lngname, flag_xyfill=.true.)
      if (m>3 .and. history_aero_optics) then
         call add_default (fldname, 1, ' ')
      endif
      !-- Longlei Li

   enddo

   call addfld ('AODDUST',       horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from dust, day only',         &
        flag_xyfill=.true.)
   call addfld ('AODDUSTS1',       horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from dust1, day only',         &
        flag_xyfill=.true.)
   call addfld ('AODDUSTS2',       horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from dust2, day only',         &
        flag_xyfill=.true.)
   call addfld ('AODDUSTS3',       horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from dust3, day only',         &
        flag_xyfill=.true.)
   call addfld ('AODDUSTS4',       horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from dust4, day only',         &
        flag_xyfill=.true.)
   call addfld ('AODDUSTS5',       horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from dust5, day only',         &
        flag_xyfill=.true.)
   call addfld ('AODDUSTS6',       horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from dust6, day only',         &
        flag_xyfill=.true.)
   call addfld ('AODDUSTS7',       horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from dust7, day only',         &
        flag_xyfill=.true.)
   call addfld ('AODDUSTS8',       horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from dust8, day only',         &
        flag_xyfill=.true.)
   !++ Longlei Li ++
   call addfld ('AODASH',        horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from ash, day only',          &
        flag_xyfill=.true.)
   !-- Longlei Li --
   call addfld ('AODSO4',        horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from SO4, day only',          &
        flag_xyfill=.true.)
   call addfld ('AODPOM',        horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from POM, day only',          &
        flag_xyfill=.true.)
   call addfld ('AODSOA',        horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from SOA, day only',          &
        flag_xyfill=.true.)
   call addfld ('AODBC',         horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from BC, day only',           &
        flag_xyfill=.true.)
   call addfld ('AODSS',         horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from seasalt, day only',      &
        flag_xyfill=.true.)
   call addfld ('AODABSBC',      horiz_only, 'A','  ',    'Aerosol absorption optical depth 550 nm from BC, day only',&
        flag_xyfill=.true.)
   call addfld ('BURDENDUST',    horiz_only, 'A','kg/m2', 'Dust aerosol burden, day only'        ,                    &
        flag_xyfill=.true.)
   !++ Longlei Li ++
   call addfld ('BURDENASH',     horiz_only, 'A','kg/m2', 'Ash aerosol burden, day only'         ,                    &
        flag_xyfill=.true.)
   !-- Longlei Li --
   call addfld ('BURDENSO4',     horiz_only, 'A','kg/m2', 'Sulfate aerosol burden, day only'     ,                    &
        flag_xyfill=.true.)
   call addfld ('BURDENPOM',     horiz_only, 'A','kg/m2', 'POM aerosol burden, day only'         ,                    &
        flag_xyfill=.true.)
   call addfld ('BURDENSOA',     horiz_only, 'A','kg/m2', 'SOA aerosol burden, day only'         ,                    &
        flag_xyfill=.true.)
   call addfld ('BURDENBC',      horiz_only, 'A','kg/m2', 'Black carbon aerosol burden, day only',                    &
        flag_xyfill=.true.)
   call addfld ('BURDENSEASALT', horiz_only, 'A','kg/m2', 'Seasalt aerosol burden, day only'     ,                    &
        flag_xyfill=.true.)
   call addfld ('SSAVIS',        horiz_only, 'A','  ',    'Aerosol single-scatter albedo, day only',                  &
        flag_xyfill=.true.)

   call addfld ('AODDUSTdn',       horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from dust, day night',         &
        flag_xyfill=.true.)
   call addfld ('AODDUSTS1dn',       horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from dust1, day night',         &
        flag_xyfill=.true.)
   call addfld ('AODDUSTS2dn',       horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from dust2, day night',         &
        flag_xyfill=.true.)
   call addfld ('AODDUSTS3dn',       horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from dust3, day night',         &
        flag_xyfill=.true.)
   call addfld ('AODDUSTS4dn',       horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from dust4, day night',         &
        flag_xyfill=.true.)
   call addfld ('AODDUSTS5dn',       horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from dust5, day night',         &
        flag_xyfill=.true.)
   call addfld ('AODDUSTS6dn',       horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from dust6, day night',         &
        flag_xyfill=.true.)
   call addfld ('AODDUSTS7dn',       horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from dust7, day night',         &
        flag_xyfill=.true.)
   call addfld ('AODDUSTS8dn',       horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from dust8, day night',         &
        flag_xyfill=.true.)
   call addfld ('AODSO4dn',        horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from SO4, day night',          &
        flag_xyfill=.true.)
   call addfld ('AODPOMdn',        horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from POM, day night',          &
        flag_xyfill=.true.)
   call addfld ('AODSOAdn',        horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from SOA, day night',          &
        flag_xyfill=.true.)
   call addfld ('AODBCdn',         horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from BC, day night',           &
        flag_xyfill=.true.)
   call addfld ('AODSSdn',         horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from seasalt, day night',      &
        flag_xyfill=.true.)
   call addfld ('AODABSBCdn',      horiz_only, 'A','  ',    'Aerosol absorption optical depth 550 nm from BC, day night',&
        flag_xyfill=.true.)
   call addfld ('BURDENDUSTdn',    horiz_only, 'A','kg/m2', 'Dust aerosol burden, day night'        ,                    &
        flag_xyfill=.true.)
   call addfld ('BURDENSO4dn',     horiz_only, 'A','kg/m2', 'Sulfate aerosol burden, day night'     ,                    &
        flag_xyfill=.true.)
   call addfld ('BURDENPOMdn',     horiz_only, 'A','kg/m2', 'POM aerosol burden, day night'         ,                    &
        flag_xyfill=.true.)
   call addfld ('BURDENSOAdn',     horiz_only, 'A','kg/m2', 'SOA aerosol burden, day night'         ,                    &
        flag_xyfill=.true.)
   call addfld ('BURDENBCdn',      horiz_only, 'A','kg/m2', 'Black carbon aerosol burden, day night',                    &
        flag_xyfill=.true.)
   call addfld ('BURDENSEASALTdn', horiz_only, 'A','kg/m2', 'Seasalt aerosol burden, day night'     ,                    &
        flag_xyfill=.true.)
   call addfld ('SSAVISdn',        horiz_only, 'A','  ',    'Aerosol single-scatter albedo, day night',                  &
        flag_xyfill=.true.)

 
   if (history_amwg) then 
      call add_default ('AODDUST1'     , 1, ' ')
      call add_default ('AODDUST3'     , 1, ' ')
      call add_default ('AODDUST'      , 1, ' ')
      !++ Longlei Li ++
      call add_default ('AODASH'       , 1, ' ')
      !-- Longlei Li --
      call add_default ('AODVIS'       , 1, ' ')
   end if

   if (history_dust) then 
      call add_default ('AODDUST1'     , 1, ' ')
      call add_default ('AODDUST2'     , 1, ' ')
      call add_default ('AODDUST3'     , 1, ' ')
   end if

   if (history_aero_optics) then 
      call add_default ('AODDUST1'     , 1, ' ')
      call add_default ('AODDUST3'     , 1, ' ')
      call add_default ('ABSORB'       , 1, ' ')
      call add_default ('AODMODE1'     , 1, ' ')
      call add_default ('AODMODE2'     , 1, ' ')
      call add_default ('AODMODE3'     , 1, ' ')
      call add_default ('AODVIS'       , 1, ' ')
      call add_default ('AODUV'        , 1, ' ')
      call add_default ('AODNIR'       , 1, ' ')
      call add_default ('AODABS'       , 1, ' ')
      call add_default ('AODABSBC'     , 1, ' ')
      call add_default ('AODDUST'      , 1, ' ')
      !++ Longlei Li ++
      call add_default ('AODASH'       , 1, ' ')
      !-- Longlei Li --
      call add_default ('AODSO4'       , 1, ' ')
      call add_default ('AODPOM'       , 1, ' ')
      call add_default ('AODSOA'       , 1, ' ')
      call add_default ('AODBC'        , 1, ' ')
      call add_default ('AODSS'        , 1, ' ')
      call add_default ('BURDEN1'      , 1, ' ')
      call add_default ('BURDEN2'      , 1, ' ')
      call add_default ('BURDEN3'      , 1, ' ')
      call add_default ('BURDENDUST'   , 1, ' ')
      !++ Longlei Li ++
      call add_default ('BURDENASH'    , 1, ' ')
      !-- Longlei Li --
      call add_default ('BURDENSO4'    , 1, ' ')
      call add_default ('BURDENPOM'    , 1, ' ')
      call add_default ('BURDENSOA'    , 1, ' ')
      call add_default ('BURDENBC'     , 1, ' ')
      call add_default ('BURDENSEASALT', 1, ' ')
      call add_default ('SSAVIS'       , 1, ' ')
      call add_default ('EXTINCT'      , 1, ' ')
      call add_default ('AODxASYM'     , 1, ' ')
      call add_default ('EXTxASYM'     , 1, ' ')
      !++ Longlei Li
      call add_default ('EXTDUST'      , 1, ' ')
      !-- Longlei Li
     
      call add_default ('AODdnDUST1'     , 1, ' ')
      call add_default ('AODdnDUST3'     , 1, ' ')
      call add_default ('ABSORBdn'       , 1, ' ')
      call add_default ('AODdnMODE1'     , 1, ' ')
      call add_default ('AODdnMODE2'     , 1, ' ')
      call add_default ('AODdnMODE3'     , 1, ' ')
      call add_default ('AODVISdn'       , 1, ' ')
      call add_default ('AODUVdn'        , 1, ' ')
      call add_default ('AODNIRdn'       , 1, ' ')
      call add_default ('AODABSdn'       , 1, ' ')
      call add_default ('AODABSBCdn'     , 1, ' ')
      call add_default ('AODDUSTdn'      , 1, ' ')
      call add_default ('AODSO4dn'       , 1, ' ')
      call add_default ('AODPOMdn'       , 1, ' ')
      call add_default ('AODSOAdn'       , 1, ' ')
      call add_default ('AODBCdn'        , 1, ' ')
      call add_default ('AODSSdn'        , 1, ' ')
      call add_default ('BURDENdn1'      , 1, ' ')
      call add_default ('BURDENdn2'      , 1, ' ')
      call add_default ('BURDENdn3'      , 1, ' ')
      call add_default ('BURDENDUSTdn'   , 1, ' ')
      call add_default ('BURDENSO4dn'    , 1, ' ')
      call add_default ('BURDENPOMdn'    , 1, ' ')
      call add_default ('BURDENSOAdn'    , 1, ' ')
      call add_default ('BURDENBCdn'     , 1, ' ')
      call add_default ('BURDENSEASALTdn', 1, ' ')
      call add_default ('SSAVISdn'       , 1, ' ')
      call add_default ('EXTINCTdn'      , 1, ' ')
      !++ Longlei Li
      call add_default ('EXTDUSTdn'      , 1, ' ')
      !-- Longlei Li

      call add_default ('AODxASYMdn'     , 1, ' ')
      call add_default ('EXTxASYMdn'     , 1, ' ')
  end if

   do ilist = 1, n_diag
      if (call_list(ilist)) then
         
         call addfld ('EXTINCT'//diag(ilist),  (/ 'lev' /), 'A','/m', &
              'Aerosol extinction', flag_xyfill=.true.)

         !++ Longlei Li
         call addfld ('EXTDUST'//diag(ilist),  (/ 'lev' /), 'A','/m', &
              'Dust extinction', flag_xyfill=.true.)
         !-- Longlei Li

         call addfld ('ABSORB'//diag(ilist),   (/ 'lev' /), 'A','/m', &
              'Aerosol absorption', flag_xyfill=.true.)
         call addfld ('AODVIS'//diag(ilist),   horiz_only,  'A','  ', &
              'Aerosol optical depth 550 nm', flag_xyfill=.true.)
         call addfld ('AODVISst'//diag(ilist), horiz_only,  'A','  ', &
              'Stratospheric aerosol optical depth 550 nm', flag_xyfill=.true.)
         call addfld ('AODABS'//diag(ilist),   horiz_only,  'A','  ', &
              'Aerosol absorption optical depth 550 nm', flag_xyfill=.true.)

         call addfld ('EXTINCTdn'//diag(ilist),    (/ 'lev' /), 'A','/m',&
              'Aerosol extinction 550 nm, day night', flag_xyfill=.true.)
         !++ Longlei Li
         call addfld ('EXTDUSTdn'//diag(ilist),    (/ 'lev' /), 'A','/m',&
              'Aerosol extinction 550 nm, day night', flag_xyfill=.true.)
         !-- Longlei Li

         call addfld ('ABSORBdn'//diag(ilist),     (/ 'lev' /), 'A','/m',&
              'Aerosol absorption, day night',  flag_xyfill=.true.)
         call addfld ('AODVISdn'//diag(ilist),     horiz_only,  'A','  ',&
              'Aerosol optical depth 550 nm, day night',  flag_xyfill=.true.)
         call addfld ('AODVISstdn'//diag(ilist),   horiz_only,  'A','  ',&
              'Stratospheric aerosol optical depth 550 nm, day night', flag_xyfill=.true.)
         call addfld ('AODABSdn'//diag(ilist),     horiz_only,  'A','  ',&
              'Aerosol absorption optical depth 550 nm, day night',  flag_xyfill=.true.)
         call addfld ('EXTxASYMdn'//diag(ilist),   (/ 'lev' /), 'A','  ',&
              'extinction 550 * asymmetry factor, day night',  flag_xyfill=.true.)
         call addfld ('EXTxASYM'//diag(ilist),   (/ 'lev' /), 'A','  ',&
              'extinction 550 nm * asymmetry factor, day only',   flag_xyfill=.true.)

         if (history_aero_optics) then
            call add_default ('EXTINCT'//diag(ilist), 1, ' ')
            !++ Longlei Li
            call add_default ('EXTDUST'//diag(ilist), 1, ' ')
            !-- Longlei Li
            call add_default ('ABSORB'//diag(ilist),  1, ' ')
            call add_default ('AODVIS'//diag(ilist),  1, ' ')
            call add_default ('AODVISst'//diag(ilist),  1, ' ')
            call add_default ('AODABS'//diag(ilist),  1, ' ')
         end if

      end if
   end do

end subroutine modal_aer_opt_init

!===============================================================================

subroutine modal_aero_sw(list_idx, state, pbuf, nnite, idxnite, &
                         tauxar, wa, ga, fa)

   ! calculates aerosol sw radiative properties
   
   use tropopause, only : tropopause_findChemTrop

   integer,             intent(in) :: list_idx       ! index of the climate or a diagnostic list
   type(physics_state), intent(in), target :: state          ! state variables
   
   type(physics_buffer_desc), pointer :: pbuf(:)
   integer,             intent(in) :: nnite          ! number of night columns
   integer,             intent(in) :: idxnite(nnite) ! local column indices of night columns

   real(r8), intent(out) :: tauxar(pcols,0:pver,nswbands) ! layer extinction optical depth
   real(r8), intent(out) :: wa(pcols,0:pver,nswbands)     ! layer single-scatter albedo
   real(r8), intent(out) :: ga(pcols,0:pver,nswbands)     ! asymmetry factor
   real(r8), intent(out) :: fa(pcols,0:pver,nswbands)     ! forward scattered fraction

   ! Local variables
   integer :: i, ifld, isw, k, l, m, nc, ns
   integer :: lchnk                    ! chunk id
   integer :: ncol                     ! number of active columns in the chunk
   integer :: nmodes
   integer :: nspec
   integer :: troplevchem(pcols)       ! Chemical tropopause level
   integer :: istat

   real(r8) :: mass(pcols,pver)        ! layer mass
   real(r8) :: air_density(pcols,pver) ! (kg/m3)

   real(r8),    pointer :: specmmr(:,:)        ! species mass mixing ratio
   real(r8)             :: specdens            ! species density (kg/m3)
   complex(r8), pointer :: specrefindex(:)     ! species refractive index
   character*32         :: spectype            ! species type
   real(r8)             :: hygro_aer           ! 

   real(r8), pointer :: dgnumwet(:,:)     ! number mode wet diameter
   real(r8), pointer :: qaerwat(:,:)      ! aerosol water (g/g)

   real(r8), pointer :: dgnumdry_m(:,:,:) ! number mode dry diameter for all modes
   real(r8), pointer :: dgnumwet_m(:,:,:) ! number mode wet diameter for all modes
   real(r8), pointer :: qaerwat_m(:,:,:)  ! aerosol water (g/g) for all modes
   real(r8), pointer :: wetdens_m(:,:,:)  ! 
   real(r8), pointer :: hygro_m(:,:,:)  ! 
   real(r8), pointer :: dryvol_m(:,:,:)  ! 
   real(r8), pointer :: dryrad_m(:,:,:)  ! 
   real(r8), pointer :: drymass_m(:,:,:)  ! 
   real(r8), pointer :: so4dryvol_m(:,:,:)  ! 
   real(r8), pointer :: naer_m(:,:,:)  ! 

   real(r8) :: sigma_logr_aer         ! geometric standard deviation of number distribution
   real(r8) :: radsurf(pcols,pver)    ! aerosol surface mode radius
   real(r8) :: logradsurf(pcols,pver) ! log(aerosol surface mode radius)
   real(r8) :: cheb(ncoef,pcols,pver)

   real(r8)    :: refr(pcols)     ! real part of refractive index
   real(r8)    :: refi(pcols)     ! imaginary part of refractive index
   complex(r8) :: crefin(pcols)   ! complex refractive index
   real(r8), pointer :: refrtabsw(:,:) ! table of real refractive indices for aerosols
   real(r8), pointer :: refitabsw(:,:) ! table of imag refractive indices for aerosols
   real(r8), pointer :: extpsw(:,:,:,:) ! specific extinction
   real(r8), pointer :: abspsw(:,:,:,:) ! specific absorption
   real(r8), pointer :: asmpsw(:,:,:,:) ! asymmetry factor

   real(r8) :: vol(pcols)      ! volume concentration of aerosol specie (m3/kg)
   real(r8) :: dryvol(pcols)   ! volume concentration of aerosol mode (m3/kg)
   real(r8) :: watervol(pcols) ! volume concentration of water in each mode (m3/kg)
   real(r8) :: wetvol(pcols)   ! volume concentration of wet mode (m3/kg)

   integer  :: itab(pcols), jtab(pcols)
   real(r8) :: ttab(pcols), utab(pcols)
   real(r8) :: cext(pcols,ncoef), cabs(pcols,ncoef), casm(pcols,ncoef)
   real(r8) :: pext(pcols)     ! parameterized specific extinction (m2/kg)
   real(r8) :: specpext(pcols) ! specific extinction (m2/kg)
   real(r8) :: dopaer(pcols)   ! aerosol optical depth in layer
   real(r8) :: pabs(pcols)     ! parameterized specific absorption (m2/kg)
   real(r8) :: pasm(pcols)     ! parameterized asymmetry factor
   real(r8) :: palb(pcols)     ! parameterized single scattering albedo

   ! Diagnostics
   real(r8) :: extinct(pcols,pver)
   ! ++ Longlei Li
   real(r8) :: extdust(pcols,pver)
   real(r8) :: extdust1(pcols,pver),extdust2(pcols,pver),extdust3(pcols,pver),extdust4(pcols,pver),&
               extdust5(pcols,pver),extdust6(pcols,pver),extdust7(pcols,pver),extdust8(pcols,pver)
   real(r8) :: dustextmode(pcols,pver)     ! mode dust extinction 

   real(r8) :: ashvol(pcols)
   real(r8) :: ashaodmode(pcols)
   ! -- Longlei Li

   real(r8) :: extinctnir(pcols,pver)
   real(r8) :: extinctuv(pcols,pver)
   real(r8) :: absorb(pcols,pver)
   real(r8) :: aodvis(pcols)               ! extinction optical depth
   real(r8) :: aodvisst(pcols)             ! stratospheric extinction optical depth
   real(r8) :: aodabs(pcols)               ! absorption optical depth
   real(r8) :: asymvis(pcols)              ! asymmetry factor * optical depth
   real(r8) :: asymext(pcols,pver)         ! asymmetry factor * extinction

   real(r8) :: aodabsbc(pcols)             ! absorption optical depth of BC

   real(r8) :: ssavis(pcols)
   real(r8) :: dustvol(pcols)              ! volume concentration of dust in aerosol mode (m3/kg)
   real(r8) :: dustvol1(pcols),dustvol2(pcols),dustvol3(pcols),dustvol4(pcols),&
               dustvol5(pcols),dustvol6(pcols),dustvol7(pcols),dustvol8(pcols)

   ! ++ Longei Li ++         

   real(r8) :: fevol(pcols)              ! volume concentration of iron in aerosol mode (m3/kg)
   real(r8) :: fetmvol(pcols),fesmvol(pcols),fetsvol(pcols),fessvol(pcols),&
               fetfvol(pcols),fesfvol(pcols),fetcvol(pcols),fescvol(pcols)
   ! -- Longlei Li -- 

   real(r8) :: burden(pcols)
   real(r8) :: burdendust(pcols), burdenso4(pcols), burdenbc(pcols), &
               burdenpom(pcols), burdensoa(pcols), burdenseasalt(pcols)
! separate dust bins
   real(r8) :: burdendust1(pcols),burdendust2(pcols),burdendust3(pcols), &
               burdendust4(pcols),burdendust5(pcols),burdendust6(pcols), &
               burdendust7(pcols),burdendust8(pcols)

   !++ Longlei Li ++
   real(r8) :: burdenash(pcols), burdenfe(pcols)
   real(r8) :: burdenfetm(pcols),burdenfesm(pcols),burdenfets(pcols), &
               burdenfess(pcols),burdenfetf(pcols),burdenfesf(pcols), &
               burdenfetc(pcols),burdenfesc(pcols)
   !-- Longlei Li --

   real(r8) :: aodmode(pcols)
   real(r8) :: dustaodmode(pcols)          ! dust aod in aerosol mode

   real(r8) :: specrefr, specrefi
   real(r8) :: scatdust(pcols), scatso4(pcols), scatbc(pcols), &
               scatpom(pcols), scatsoa(pcols), scatseasalt(pcols)
   real(r8) :: scatdust1(pcols),scatdust2(pcols),scatdust3(pcols), &
               scatdust4(pcols),scatdust5(pcols),scatdust6(pcols), &
               scatdust7(pcols),scatdust8(pcols)
   !++ Longlei Li ++
   real(r8) :: scatash(pcols), scatfe(pcols)
   real(r8) :: scatfetm(pcols),scatfesm(pcols),scatfets(pcols), &
               scatfess(pcols),scatfetf(pcols),scatfesf(pcols), &
               scatfetc(pcols),scatfesc(pcols)
   !-- Longlei Li --

   real(r8) :: absdust(pcols), absso4(pcols), absbc(pcols), &
               abspom(pcols), abssoa(pcols), absseasalt(pcols)
   real(r8) :: absdust1(pcols),absdust2(pcols),absdust3(pcols), &
               absdust4(pcols),absdust5(pcols),absdust6(pcols), &
               absdust7(pcols),absdust8(pcols)
   !++ Longlei Li ++
   real(r8) :: absash(pcols), absfe(pcols)
   real(r8) :: absfetm(pcols),absfesm(pcols),absfets(pcols), &
               absfess(pcols),absfetf(pcols),absfesf(pcols), &
               absfetc(pcols),absfesc(pcols)
   !-- Longlei Li --

   real(r8) :: hygrodust(pcols), hygroso4(pcols), hygrobc(pcols), &
               hygropom(pcols), hygrosoa(pcols), hygroseasalt(pcols)
   real(r8) :: hygrodust1(pcols),hygrodust2(pcols),hygrodust3(pcols), &
               hygrodust4(pcols),hygrodust5(pcols),hygrodust6(pcols), &
               hygrodust7(pcols),hygrodust8(pcols)

   !++ Longlei Li ++
   real(r8) :: hygroash(pcols), hygrofe(pcols)
   real(r8) :: hygrofetm(pcols),hygrofesm(pcols),hygrofets(pcols), &
               hygrofess(pcols),hygrofetf(pcols),hygrofesf(pcols), &
               hygrofetc(pcols),hygrofesc(pcols)
   !-- Longlei Li --

   real(r8) :: scath2o, absh2o, sumscat, sumabs, sumhygro
   real(r8) :: aodc                        ! aod of component
   !++ Longlei Li ++ 
   real(r8) :: aodcash
   !-- Longlei Li --   

   real(r8) :: aodc1,aodc2,aodc3,aodc4,aodc5,aodc6,aodc7,aodc8                        ! aod of component
   real(r8) :: eaod                        ! enhanced aod 

   ! total species AOD
   real(r8) :: dustaod(pcols), so4aod(pcols), bcaod(pcols), &
               pomaod(pcols), soaaod(pcols), seasaltaod(pcols)
   !real(r8) :: dustaod1(pcols),dustaod2(pcols),dustaod3(pcols), &
   !            dustaod4(pcols),dustaod5(pcols),dustaod6(pcols), &
   !            dustaod7(pcols),dustaod8(pcols)

   ! ++ Longlei Li ++
   real(r8) :: ashaod(pcols)
   ! -- Longlei Li --


   logical :: savaervis ! true if visible wavelength (0.55 micron)
   logical :: savaernir ! true if near ir wavelength (~0.88 micron)
   logical :: savaeruv  ! true if uv wavelength (~0.35 micron)

   real(r8) :: aoduv(pcols)               ! extinction optical depth in uv
   real(r8) :: aoduvst(pcols)             ! stratospheric extinction optical depth in uv
   real(r8) :: aodnir(pcols)              ! extinction optical depth in nir
   real(r8) :: aodnirst(pcols)            ! stratospheric extinction optical depth in nir


   character(len=32) :: outname

   ! debug output
   integer, parameter :: nerrmax_dopaer=1000
   integer  :: nerr_dopaer = 0
   real(r8) :: volf            ! volume fraction of insoluble aerosol
   character(len=*), parameter :: subname = 'modal_aero_sw'
   !----------------------------------------------------------------------------
    !write(iulog,*) 'Longlei marks begin in modal_aer_opt'
   lchnk = state%lchnk
   ncol  = state%ncol

   ! initialize output variables
   tauxar(:ncol,:,:) = 0._r8
   wa(:ncol,:,:)     = 0._r8
   ga(:ncol,:,:)     = 0._r8
   fa(:ncol,:,:)     = 0._r8

   ! zero'th layer does not contain aerosol
   tauxar(1:ncol,0,:)  = 0._r8
   wa(1:ncol,0,:)      = 0.925_r8
   ga(1:ncol,0,:)      = 0.850_r8
   fa(1:ncol,0,:)      = 0.7225_r8

   mass(:ncol,:)        = state%pdeldry(:ncol,:)*rga
   air_density(:ncol,:) = state%pmid(:ncol,:)/(rair*state%t(:ncol,:))

   ! diagnostics for visible band summed over modes
   extinct(1:ncol,:)     = 0.0_r8
   ! ++ Longlei Li
   extdust(1:ncol,:)     = 0.0_r8
   extdust1(1:ncol,:)     = 0.0_r8
   extdust2(1:ncol,:)     = 0.0_r8
   extdust3(1:ncol,:)     = 0.0_r8
   extdust4(1:ncol,:)     = 0.0_r8
   extdust5(1:ncol,:)     = 0.0_r8
   extdust6(1:ncol,:)     = 0.0_r8
   extdust7(1:ncol,:)     = 0.0_r8
   extdust8(1:ncol,:)     = 0.0_r8
   ! -- Longlei Li
   absorb(1:ncol,:)      = 0.0_r8
   aodvis(1:ncol)        = 0.0_r8
   aodvisst(1:ncol)      = 0.0_r8
   aodabs(1:ncol)        = 0.0_r8
   burdendust(:ncol)     = 0.0_r8
   ! ++ Longlei Li ++ 
   burdenash(:ncol)      = 0.0_r8
   ! -- Longlei Li --
   burdenso4(:ncol)      = 0.0_r8
   burdenpom(:ncol)      = 0.0_r8
   burdensoa(:ncol)      = 0.0_r8
   burdenbc(:ncol)       = 0.0_r8
   burdenseasalt(:ncol)  = 0.0_r8
   ssavis(1:ncol)        = 0.0_r8
   asymvis(1:ncol)       = 0.0_r8
   asymext(1:ncol,:)     = 0.0_r8
   burdendust1(:ncol)     = 0.0_r8
   burdendust2(:ncol)     = 0.0_r8
   burdendust3(:ncol)     = 0.0_r8
   burdendust4(:ncol)     = 0.0_r8
   burdendust5(:ncol)     = 0.0_r8
   burdendust6(:ncol)     = 0.0_r8
   burdendust7(:ncol)     = 0.0_r8
   burdendust8(:ncol)     = 0.0_r8

   aodabsbc(:ncol)       = 0.0_r8 
   dustaod(:ncol)        = 0.0_r8
   ! ++ Longlei Li ++ 
   ashaod(:ncol)        = 0.0_r8
   ! -- Longlei Li --

   so4aod(:ncol)         = 0.0_r8
   pomaod(:ncol)         = 0.0_r8
   soaaod(:ncol)         = 0.0_r8
   bcaod(:ncol)          = 0.0_r8
   seasaltaod(:ncol)     = 0.0_r8
   !dustaod1(:ncol)     = 0.0_r8
   !dustaod2(:ncol)     = 0.0_r8
   !dustaod3(:ncol)     = 0.0_r8
   !dustaod4(:ncol)     = 0.0_r8
   !dustaod5(:ncol)     = 0.0_r8
   !dustaod6(:ncol)     = 0.0_r8
   !dustaod7(:ncol)     = 0.0_r8
   !dustaod8(:ncol)     = 0.0_r8

   ! diags for other bands
   extinctuv(1:ncol,:)   = 0.0_r8
   extinctnir(1:ncol,:)  = 0.0_r8
   aoduv(:ncol)          = 0.0_r8
   aodnir(:ncol)         = 0.0_r8
   aoduvst(:ncol)        = 0.0_r8
   aodnirst(:ncol)       = 0.0_r8
   call tropopause_findChemTrop(state, troplevchem)

   ! loop over all aerosol modes
   call rad_cnst_get_info(list_idx, nmodes=nmodes)

   if (list_idx == 0) then
      ! water uptake and wet radius for the climate list has already been calculated
      call pbuf_get_field(pbuf, dgnumwet_idx, dgnumwet_m)
      call pbuf_get_field(pbuf, qaerwat_idx,  qaerwat_m)
   else
      ! If doing a diagnostic calculation then need to calculate the wet radius
      ! and water uptake for the diagnostic modes
      allocate(dgnumdry_m(pcols,pver,nmodes),  dgnumwet_m(pcols,pver,nmodes), &
               qaerwat_m(pcols,pver,nmodes),   wetdens_m(pcols,pver,nmodes), &
               hygro_m(pcols,pver,nmodes),     dryvol_m(pcols,pver,nmodes), &
               dryrad_m(pcols,pver,nmodes),    drymass_m(pcols,pver,nmodes),  &
               so4dryvol_m(pcols,pver,nmodes), naer_m(pcols,pver,nmodes),     stat=istat)
      if (istat > 0) then
         call endrun('modal_aero_sw: allocation FAILURE: arrays for diagnostic calcs')
      end if
      call modal_aero_calcsize_diag(state, pbuf, list_idx, dgnumdry_m, hygro_m, &
                                    dryvol_m, dryrad_m, drymass_m, so4dryvol_m, naer_m)  
      call modal_aero_wateruptake_dr(state, pbuf, list_idx, dgnumdry_m, dgnumwet_m, &
                                     qaerwat_m, wetdens_m,  hygro_m, dryvol_m, dryrad_m, &
                                     drymass_m, so4dryvol_m, naer_m)
   endif

   do m = 1, nmodes

      ! diagnostics for visible band for each mode
      burden(:ncol)       = 0._r8
      aodmode(1:ncol)     = 0.0_r8
      dustaodmode(1:ncol) = 0.0_r8
      ! ++ Longlei Li
      dustextmode(1:ncol,:) = 0.0_r8
      ! -- Longlei Li

      dgnumwet => dgnumwet_m(:,:,m)
      qaerwat  => qaerwat_m(:,:,m)

      ! get mode properties
      call rad_cnst_get_mode_props(list_idx, m, sigmag=sigma_logr_aer, refrtabsw=refrtabsw , &
         refitabsw=refitabsw, extpsw=extpsw, abspsw=abspsw, asmpsw=asmpsw)

      ! get mode info
      call rad_cnst_get_info(list_idx, m, nspec=nspec)

      ! calc size parameter for all columns
      call modal_size_parameters(ncol, sigma_logr_aer, dgnumwet, radsurf, logradsurf, cheb)

      do isw = 1, nswbands
         savaervis = (isw .eq. idx_sw_diag)
         savaeruv  = (isw .eq. idx_uv_diag)
         savaernir = (isw .eq. idx_nir_diag)

         do k = top_lev, pver

            ! form bulk refractive index
            crefin(:ncol) = (0._r8, 0._r8)
            dryvol(:ncol) = 0._r8
            dustvol(:ncol) = 0._r8

            !++ Longlei Li ++
            ashvol(:ncol)  = 0._r8
            !-- Longlei Li --

            scatdust(:ncol)     = 0._r8
            absdust(:ncol)      = 0._r8
            hygrodust(:ncol)    = 0._r8
            scatso4(:ncol)      = 0._r8
            absso4(:ncol)       = 0._r8
            hygroso4(:ncol)     = 0._r8
            scatbc(:ncol)       = 0._r8
            absbc(:ncol)        = 0._r8
            hygrobc(:ncol)      = 0._r8
            scatpom(:ncol)      = 0._r8
            abspom(:ncol)       = 0._r8
            hygropom(:ncol)     = 0._r8
            scatsoa(:ncol)      = 0._r8
            abssoa(:ncol)       = 0._r8
            hygrosoa(:ncol)     = 0._r8
            scatseasalt(:ncol)  = 0._r8
            absseasalt(:ncol)   = 0._r8
            hygroseasalt(:ncol) = 0._r8

            !++ Longlei Li ++
            hygroash(:ncol)     = 0._r8
            scatash(:ncol)      = 0._r8
            absash(:ncol)       = 0._r8
            !-- Longlei Li --

            dustvol1(:ncol) = 0._r8
            dustvol2(:ncol) = 0._r8
            dustvol3(:ncol) = 0._r8
            dustvol4(:ncol) = 0._r8
            dustvol5(:ncol) = 0._r8
            dustvol6(:ncol) = 0._r8
            dustvol7(:ncol) = 0._r8
            dustvol8(:ncol) = 0._r8

            hygrodust1(:ncol)    = 0._r8
            scatdust1(:ncol)     = 0._r8
            absdust1(:ncol)      = 0._r8
            hygrodust2(:ncol)    = 0._r8
            scatdust2(:ncol)     = 0._r8
            absdust2(:ncol)      = 0._r8
            hygrodust3(:ncol)    = 0._r8
            scatdust3(:ncol)     = 0._r8
            absdust3(:ncol)      = 0._r8
            hygrodust4(:ncol)    = 0._r8
            scatdust4(:ncol)     = 0._r8
            absdust4(:ncol)      = 0._r8
            hygrodust5(:ncol)    = 0._r8
            scatdust5(:ncol)     = 0._r8
            absdust5(:ncol)      = 0._r8
            hygrodust6(:ncol)    = 0._r8
            scatdust6(:ncol)     = 0._r8
            absdust6(:ncol)      = 0._r8
            hygrodust7(:ncol)    = 0._r8
            scatdust7(:ncol)     = 0._r8
            absdust7(:ncol)      = 0._r8
            hygrodust8(:ncol)    = 0._r8
            scatdust8(:ncol)     = 0._r8
            absdust8(:ncol)      = 0._r8

            ! aerosol species loop
            do l = 1, nspec
               call rad_cnst_get_aer_mmr(list_idx, m, l, 'a', state, pbuf, specmmr)
               call rad_cnst_get_aer_props(list_idx, m, l, density_aer=specdens, &
                                           refindex_aer_sw=specrefindex, spectype=spectype, &
                                           hygro_aer=hygro_aer)
!++Longlei Li
                !write( iulog, '(2a,3(1x,i12),2x,i12,a)' ) &
                !     'subr initaermodes_setspecptrs lptr_', &
                !     ' - m, l, lmassa, cnst_name = ', m, l, lmassa,pcnst, cnst_name(lmassa)
                !write( iulog, '(2a,2x,5(1pe14.5))' ) &
                !     'For high AOD, Longlei', &
                !     ' specdens = ', specdens

                !if (k == pver) then
                !    write( iulog, '(2a,2x,a)' ) &
                !         'For high AOD, Longlei', &
                !         ' speci name = ', spectype
                !end if
!--Longlei Li

!++Longlei, force sulfate to be background level 
               !if (m == 3 .and. trim(spectype) == 'sulfate') then
               !if ( (m == 3 .or. m == 1) .and. l == 1) then
                  !specmmr(1:ncol,top_lev:pver)  = 1.0e-36_r8
                  !qaerwat(1:ncol,top_lev:pver)  = 1.0e-36_r8
               !end if
!--

               do i = 1, ncol
                  vol(i)      = specmmr(i,k)/specdens
                  dryvol(i)   = dryvol(i) + vol(i)
                  crefin(i)   = crefin(i) + vol(i)*specrefindex(isw)
!++Longlei Li
                !if (k == pver) then
                !    write( iulog, '(2a,2x,5(1pe14.5))' ) &
                !         'For high AOD, Longlei', &
                !         ' specmmr(i,k) = , vol(i) = ', specmmr(i,k), vol(i)
                !end if
!--Longlei Li
               end do

               ! compute some diagnostics for visible band only
               if (savaervis) then

                  specrefr = real(specrefindex(isw))
                  specrefi = aimag(specrefindex(isw))

                  do i = 1, ncol
                     burden(i) = burden(i) + specmmr(i,k)*mass(i,k)
                  end do

!                  if (trim(spectype) == 'dust') then
!                     do i = 1, ncol
!                        burdendust(i) = burdendust(i) + specmmr(i,k)*mass(i,k)
!                        dustvol(i)    = vol(i)
!                        scatdust(i)   = vol(i)*specrefr
!                        absdust(i)    = -vol(i)*specrefi
!                        hygrodust(i)  = vol(i)*hygro_aer
!                     end do
!                  end if

#if ( defined DUSTTYPES )
                  if (trim(spectype) == 'dust-type1') then
                     do i = 1, ncol
                        !++ Longlei Li, ash now (no longer dust) ++
                        burdenash(i) = burdenash(i) + specmmr(i,k)*mass(i,k)
                        ashvol(i)    = ashvol(i) + vol(i)
                        scatash(i)   = scatash(i) + vol(i)*specrefr
                        absash(i)    = absash(i) - vol(i)*specrefi
                        hygroash(i)  = hygroash(i) + vol(i)*hygro_aer
                        burdendust1(i) = burdendust1(i) + specmmr(i,k)*mass(i,k)
                        dustvol1(i)    = dustvol1(i) + vol(i)
                        scatdust1(i)   = scatdust1(i) + vol(i)*specrefr
                        absdust1(i)    = absdust1(i) - vol(i)*specrefi
                        hygrodust1(i)  = hygrodust1(i) + vol(i)*hygro_aer
                        !-- Longlei Li --
                     end do
                  end if
                  if (trim(spectype) == 'dust-type2') then
                     do i = 1, ncol
                        !++ Longlei Li, ash now (no longer dust) ++
                        burdenash(i) = burdenash(i) + specmmr(i,k)*mass(i,k)
                        ashvol(i)    = ashvol(i) + vol(i)
                        scatash(i)   = scatash(i) + vol(i)*specrefr
                        absash(i)    = absash(i) - vol(i)*specrefi
                        hygroash(i)  = hygroash(i) + vol(i)*hygro_aer
                        burdendust2(i) = burdendust2(i) + specmmr(i,k)*mass(i,k)
                        dustvol2(i)    = dustvol2(i) + vol(i)
                        scatdust2(i)   = scatdust2(i) + vol(i)*specrefr
                        absdust2(i)    = absdust2(i) - vol(i)*specrefi
                        hygrodust2(i)  = hygrodust2(i) + vol(i)*hygro_aer
                        !-- Longlei Li --
                     end do
                  end if
                  if (trim(spectype) == 'dust-type3') then
                     do i = 1, ncol
                        !++ Longlei Li, ash now (no longer dust) ++
                        burdenash(i) = burdenash(i) + specmmr(i,k)*mass(i,k)
                        ashvol(i)    = ashvol(i) + vol(i)
                        scatash(i)   = scatash(i) + vol(i)*specrefr
                        absash(i)    = absash(i) - vol(i)*specrefi
                        hygroash(i)  = hygroash(i) + vol(i)*hygro_aer
                        burdendust3(i) = burdendust3(i) + specmmr(i,k)*mass(i,k)
                        dustvol3(i)    = dustvol3(i) + vol(i)
                        scatdust3(i)   = scatdust3(i) + vol(i)*specrefr
                        absdust3(i)    = absdust3(i) - vol(i)*specrefi
                        hygrodust3(i)  = hygrodust3(i) + vol(i)*hygro_aer
                        !-- Longlei Li --
                     end do
                  end if
                  if (trim(spectype) == 'dust-type4') then
                     do i = 1, ncol
                        !++ Longlei Li, ash now (no longer dust) ++
                        burdenash(i) = burdenash(i) + specmmr(i,k)*mass(i,k)
                        ashvol(i)    = ashvol(i) + vol(i)
                        scatash(i)   = scatash(i) + vol(i)*specrefr
                        absash(i)    = absash(i) - vol(i)*specrefi
                        hygroash(i)  = hygroash(i) + vol(i)*hygro_aer
                        burdendust4(i) = burdendust4(i) + specmmr(i,k)*mass(i,k)
                        dustvol4(i)    = dustvol4(i) + vol(i)
                        scatdust4(i)   = scatdust4(i) + vol(i)*specrefr
                        absdust4(i)    = absdust4(i) - vol(i)*specrefi
                        hygrodust4(i)  = hygrodust4(i) + vol(i)*hygro_aer
                        !-- Longlei Li --
                     end do
                  end if
                  if (trim(spectype) == 'dust-type5') then
                     do i = 1, ncol
                        burdendust(i) = burdendust(i) + specmmr(i,k)*mass(i,k)
                        dustvol(i)    = dustvol(i) + vol(i)
                        scatdust(i)   = scatdust(i) + vol(i)*specrefr
                        absdust(i)    = absdust(i) - vol(i)*specrefi
                        hygrodust(i)  = hygrodust(i) + vol(i)*hygro_aer
                        burdendust5(i) = burdendust5(i) + specmmr(i,k)*mass(i,k)
                        dustvol5(i)    = dustvol5(i) + vol(i)
                        scatdust5(i)   = scatdust5(i) + vol(i)*specrefr
                        absdust5(i)    = absdust5(i) - vol(i)*specrefi
                        hygrodust5(i)  = hygrodust5(i) + vol(i)*hygro_aer
                     end do
                  end if
                  if (trim(spectype) == 'dust-type6') then
                     do i = 1, ncol
                        burdendust(i) = burdendust(i) + specmmr(i,k)*mass(i,k)
                        dustvol(i)    = dustvol(i) + vol(i)
                        scatdust(i)   = scatdust(i) + vol(i)*specrefr
                        absdust(i)    = absdust(i) - vol(i)*specrefi
                        hygrodust(i)  = hygrodust(i) + vol(i)*hygro_aer
                        burdendust6(i) = burdendust6(i) + specmmr(i,k)*mass(i,k)
                        dustvol6(i)    = dustvol6(i) + vol(i)
                        scatdust6(i)   = scatdust6(i) + vol(i)*specrefr
                        absdust6(i)    = absdust6(i) - vol(i)*specrefi
                        hygrodust6(i)  = hygrodust6(i) + vol(i)*hygro_aer
                     end do
                  end if
                  if (trim(spectype) == 'dust-type7') then
                     do i = 1, ncol
                        burdendust(i) = burdendust(i) + specmmr(i,k)*mass(i,k)
                        dustvol(i)    = dustvol(i) + vol(i)
                        scatdust(i)   = scatdust(i) + vol(i)*specrefr
                        absdust(i)    = absdust(i) - vol(i)*specrefi
                        hygrodust(i)  = hygrodust(i) + vol(i)*hygro_aer
                        burdendust7(i) = burdendust7(i) + specmmr(i,k)*mass(i,k)
                        dustvol7(i)    = dustvol7(i) + vol(i)
                        scatdust7(i)   = scatdust7(i) + vol(i)*specrefr
                        absdust7(i)    = absdust7(i) - vol(i)*specrefi
                        hygrodust7(i)  = hygrodust7(i) + vol(i)*hygro_aer
                     end do
                  end if
                  if (trim(spectype) == 'dust-type8') then
                     do i = 1, ncol
                        burdendust(i) = burdendust(i) + specmmr(i,k)*mass(i,k)
                        dustvol(i)    = dustvol(i) + vol(i)
                        scatdust(i)   = scatdust(i) + vol(i)*specrefr
                        absdust(i)    = absdust(i) - vol(i)*specrefi
                        hygrodust(i)  = hygrodust(i) + vol(i)*hygro_aer
                        burdendust8(i) = burdendust8(i) + specmmr(i,k)*mass(i,k)
                        dustvol8(i)    = dustvol8(i) + vol(i)
                        scatdust8(i)   = scatdust8(i) + vol(i)*specrefr
                        absdust8(i)    = absdust8(i) - vol(i)*specrefi
                        hygrodust8(i)  = hygrodust8(i) + vol(i)*hygro_aer
                     end do
                  end if
! ++ Longlei Li ++
#if ( defined IRONPROC )
                  if (trim(spectype) == 'fetm') then
                     do i = 1, ncol
                        burdenfe(i) = burdenfe(i) + specmmr(i,k)*mass(i,k)
                        fevol(i)    = fevol(i) + vol(i)
                        scatfe(i)   = scatfe(i) + vol(i)*specrefr
                        absfe(i)    = absfe(i) - vol(i)*specrefi
                        hygrofe(i)  = hygrofe(i) + vol(i)*hygro_aer

                        burdenfetm(i) = burdenfetm(i) + specmmr(i,k)*mass(i,k)
                        fetmvol(i)    = fetmvol(i) + vol(i)
                        scatfetm(i)   = scatfetm(i) + vol(i)*specrefr
                        absfetm(i)    = absfetm(i) - vol(i)*specrefi
                        hygrofetm(i)  = hygrofetm(i) + vol(i)*hygro_aer
                     end do
                  end if
                  if (trim(spectype) == 'fesm') then
                     do i = 1, ncol
                        burdenfe(i) = burdenfe(i) + specmmr(i,k)*mass(i,k)
                        fevol(i)    = fevol(i) + vol(i)
                        scatfe(i)   = scatfe(i) + vol(i)*specrefr
                        absfe(i)    = absfe(i) - vol(i)*specrefi
                        hygrofe(i)  = hygrofe(i) + vol(i)*hygro_aer

                        burdenfesm(i) = burdenfesm(i) + specmmr(i,k)*mass(i,k)
                        fesmvol(i)    = fesmvol(i) + vol(i)
                        scatfesm(i)   = scatfesm(i) + vol(i)*specrefr
                        absfesm(i)    = absfesm(i) - vol(i)*specrefi
                        hygrofesm(i)  = hygrofesm(i) + vol(i)*hygro_aer
                     end do
                  end if
                  if (trim(spectype) == 'fets') then
                     do i = 1, ncol
                        burdenfe(i) = burdenfe(i) + specmmr(i,k)*mass(i,k)
                        fevol(i)    = fevol(i) + vol(i)
                        scatfe(i)   = scatfe(i) + vol(i)*specrefr
                        absfe(i)    = absfe(i) - vol(i)*specrefi
                        hygrofe(i)  = hygrofe(i) + vol(i)*hygro_aer

                        burdenfets(i) = burdenfets(i) + specmmr(i,k)*mass(i,k)
                        fetsvol(i)    = fetsvol(i) + vol(i)
                        scatfets(i)   = scatfets(i) + vol(i)*specrefr
                        absfets(i)    = absfets(i) - vol(i)*specrefi
                        hygrofets(i)  = hygrofets(i) + vol(i)*hygro_aer
                     end do
                  end if
                  if (trim(spectype) == 'fess') then
                     do i = 1, ncol
                        burdenfe(i) = burdenfe(i) + specmmr(i,k)*mass(i,k)
                        fevol(i)    = fevol(i) + vol(i)
                        scatfe(i)   = scatfe(i) + vol(i)*specrefr
                        absfe(i)    = absfe(i) - vol(i)*specrefi
                        hygrofe(i)  = hygrofe(i) + vol(i)*hygro_aer

                        burdenfess(i) = burdenfess(i) + specmmr(i,k)*mass(i,k)
                        fessvol(i)    = fessvol(i) + vol(i)
                        scatfess(i)   = scatfess(i) + vol(i)*specrefr
                        absfess(i)    = absfess(i) - vol(i)*specrefi
                        hygrofess(i)  = hygrofess(i) + vol(i)*hygro_aer
                     end do
                  end if
                  if (trim(spectype) == 'fetf') then
                     do i = 1, ncol
                        burdenfe(i) = burdenfe(i) + specmmr(i,k)*mass(i,k)
                        fevol(i)    = fevol(i) + vol(i)
                        scatfe(i)   = scatfe(i) + vol(i)*specrefr
                        absfe(i)    = absfe(i) - vol(i)*specrefi
                        hygrofe(i)  = hygrofe(i) + vol(i)*hygro_aer

                        burdenfetf(i) = burdenfetf(i) + specmmr(i,k)*mass(i,k)
                        fetfvol(i)    = fetfvol(i) + vol(i)
                        scatfetf(i)   = scatfetf(i) + vol(i)*specrefr
                        absfetf(i)    = absfetf(i) - vol(i)*specrefi
                        hygrofetf(i)  = hygrofetf(i) + vol(i)*hygro_aer
                     end do
                  end if
                  if (trim(spectype) == 'fesf') then
                     do i = 1, ncol
                        burdenfe(i) = burdenfe(i) + specmmr(i,k)*mass(i,k)
                        fevol(i)    = fevol(i) + vol(i)
                        scatfe(i)   = scatfe(i) + vol(i)*specrefr
                        absfe(i)    = absfe(i) - vol(i)*specrefi
                        hygrofe(i)  = hygrofe(i) + vol(i)*hygro_aer

                        burdenfesf(i) = burdenfesf(i) + specmmr(i,k)*mass(i,k)
                        fesfvol(i)    = fesfvol(i) + vol(i)
                        scatfesf(i)   = scatfesf(i) + vol(i)*specrefr
                        absfesf(i)    = absfesf(i) - vol(i)*specrefi
                        hygrofesf(i)  = hygrofesf(i) + vol(i)*hygro_aer
                     end do
                  end if
                  if (trim(spectype) == 'fetc') then
                     do i = 1, ncol
                        burdenfe(i) = burdenfe(i) + specmmr(i,k)*mass(i,k)
                        fevol(i)    = fevol(i) + vol(i)
                        scatfe(i)   = scatfe(i) + vol(i)*specrefr
                        absfe(i)    = absfe(i) - vol(i)*specrefi
                        hygrofe(i)  = hygrofe(i) + vol(i)*hygro_aer

                        burdenfetc(i) = burdenfetc(i) + specmmr(i,k)*mass(i,k)
                        fetcvol(i)    = fetcvol(i) + vol(i)
                        scatfetc(i)   = scatfetc(i) + vol(i)*specrefr
                        absfetc(i)    = absfetc(i) - vol(i)*specrefi
                        hygrofetc(i)  = hygrofetc(i) + vol(i)*hygro_aer
                     end do
                  end if
                  if (trim(spectype) == 'fesc') then
                     do i = 1, ncol
                        burdenfe(i) = burdenfe(i) + specmmr(i,k)*mass(i,k)
                        fevol(i)    = fevol(i) + vol(i)
                        scatfe(i)   = scatfe(i) + vol(i)*specrefr
                        absfe(i)    = absfe(i) - vol(i)*specrefi
                        hygrofe(i)  = hygrofe(i) + vol(i)*hygro_aer

                        burdenfesc(i) = burdenfesc(i) + specmmr(i,k)*mass(i,k)
                        fescvol(i)    = fescvol(i) + vol(i)
                        scatfesc(i)   = scatfesc(i) + vol(i)*specrefr
                        absfesc(i)    = absfesc(i) - vol(i)*specrefi
                        hygrofesc(i)  = hygrofesc(i) + vol(i)*hygro_aer
                     end do
                  end if
#endif
! -- Longlei Li --
#else
                  if (trim(spectype) == 'dust') then
                     do i = 1, ncol
                        burdendust(i) = burdendust(i) + specmmr(i,k)*mass(i,k)
                        dustvol(i)    = vol(i)
                        scatdust(i)   = vol(i)*specrefr
                        absdust(i)    = -vol(i)*specrefi
                        hygrodust(i)  = vol(i)*hygro_aer
                     end do
                  end if
#endif
                  if (trim(spectype) == 'sulfate') then
                     do i = 1, ncol
                        burdenso4(i) = burdenso4(i) + specmmr(i,k)*mass(i,k)
                        scatso4(i)   = vol(i)*specrefr
                        absso4(i)    = -vol(i)*specrefi
                        hygroso4(i)  = vol(i)*hygro_aer
                     end do
                  end if
                  if (trim(spectype) == 'black-c') then
                     do i = 1, ncol
                        burdenbc(i) = burdenbc(i) + specmmr(i,k)*mass(i,k)
                        scatbc(i)   = vol(i)*specrefr
                        absbc(i)    = -vol(i)*specrefi
                        hygrobc(i)  = vol(i)*hygro_aer
                   end do
                  end if
                  if (trim(spectype) == 'p-organic') then
                     do i = 1, ncol
                        burdenpom(i) = burdenpom(i) + specmmr(i,k)*mass(i,k)
                        scatpom(i)   = vol(i)*specrefr
                        abspom(i)    = -vol(i)*specrefi
                        hygropom(i)  = vol(i)*hygro_aer
                      end do
                  end if
                  if (trim(spectype) == 's-organic') then
                     do i = 1, ncol
                        burdensoa(i) = burdensoa(i) + specmmr(i,k)*mass(i,k)
                        scatsoa(i)   = vol(i)*specrefr
                        abssoa(i)    = -vol(i)*specrefi
                        hygrosoa(i)  = vol(i)*hygro_aer
                     end do
                  end if
                  if (trim(spectype) == 'seasalt') then
                     do i = 1, ncol
                        burdenseasalt(i) = burdenseasalt(i) + specmmr(i,k)*mass(i,k)
                        scatseasalt(i)   = vol(i)*specrefr
                        absseasalt(i)    = -vol(i)*specrefi
                        hygroseasalt(i)  = vol(i)*hygro_aer
                      end do
                  end if

               end if
            end do ! species loop

            do i = 1, ncol
               watervol(i) = qaerwat(i,k)/rhoh2o
               wetvol(i) = watervol(i) + dryvol(i)
               if (watervol(i) < 0._r8) then
                  if (abs(watervol(i)) .gt. 1.e-1_r8*wetvol(i)) then
                     write(iulog,'(a,2e10.2,a)') 'watervol,wetvol=', &
                        watervol(i), wetvol(i), ' in '//subname
                  end if
                  watervol(i) = 0._r8
                  wetvol(i) = dryvol(i)
               end if

               ! volume mixing
               crefin(i) = crefin(i) + watervol(i)*crefwsw(isw)
               crefin(i) = crefin(i)/max(wetvol(i),1.e-60_r8)
               refr(i)   = real(crefin(i))
               refi(i)   = abs(aimag(crefin(i)))

                 ! ! Longlei tests
                 ! write(iulog,*) 'wetvol=', wetvol(i), ' dryvol=', dryvol(i), ' watervol=', watervol(i)
                 ! write(iulog,*) 'crefin=',crefin(i)
                 ! ! test ends
            end do

            ! call t_startf('binterp')

            ! interpolate coefficients linear in refractive index
            ! first call calcs itab,jtab,ttab,utab
            itab(:ncol) = 0
            call binterp(extpsw(:,:,:,isw), ncol, ncoef, prefr, prefi, &
                         refr, refi, refrtabsw(:,isw), refitabsw(:,isw), &
                         itab, jtab, ttab, utab, cext)
            call binterp(abspsw(:,:,:,isw), ncol, ncoef, prefr, prefi, &
                         refr, refi, refrtabsw(:,isw), refitabsw(:,isw), &
                         itab, jtab, ttab, utab, cabs)
            call binterp(asmpsw(:,:,:,isw), ncol, ncoef, prefr, prefi, &
                         refr, refi, refrtabsw(:,isw), refitabsw(:,isw), &
                         itab, jtab, ttab, utab, casm)

            ! call t_stopf('binterp')

            ! parameterized optical properties
            do i=1,ncol

               if (logradsurf(i,k) .le. xrmax) then
                  pext(i) = 0.5_r8*cext(i,1)
                  do nc = 2, ncoef
                     pext(i) = pext(i) + cheb(nc,i,k)*cext(i,nc)
                  enddo
                  pext(i) = exp(pext(i))
               else
                  pext(i) = 1.5_r8/(radsurf(i,k)*rhoh2o) ! geometric optics
               endif

               ! convert from m2/kg water to m2/kg aerosol
               specpext(i) = pext(i)
               pext(i) = pext(i)*wetvol(i)*rhoh2o
               pabs(i) = 0.5_r8*cabs(i,1)
               pasm(i) = 0.5_r8*casm(i,1)
               do nc = 2, ncoef
                  pabs(i) = pabs(i) + cheb(nc,i,k)*cabs(i,nc)
                  pasm(i) = pasm(i) + cheb(nc,i,k)*casm(i,nc)
               enddo
               pabs(i) = pabs(i)*wetvol(i)*rhoh2o
               pabs(i) = max(0._r8,pabs(i))
               pabs(i) = min(pext(i),pabs(i))

               palb(i) = 1._r8-pabs(i)/max(pext(i),1.e-40_r8)
               palb(i) = 1._r8-pabs(i)/max(pext(i),1.e-40_r8)

               dopaer(i) = pext(i)*mass(i,k)
!++Longlei Li
                !if (k == pver) then
                !    write( iulog, '(2a,2x,6(1pe14.5))' ) &
                !         'Longlei  ', &
                !         'dryvol(i) =, wetvol(i)=, pext(i)=, mass(i,k)=, dopaer(i)= ', &
                !           dryvol(i),   wetvol(i),  pext(i),  mass(i,k),  dopaer(i)
                !end if
!--Longlei Li
            end do

            if (savaeruv) then
               do i = 1, ncol
                 extinctuv(i,k) = extinctuv(i,k) + dopaer(i)*air_density(i,k)/mass(i,k)
                 aoduv(i) = aoduv(i) + dopaer(i)
                  if (k.le.troplevchem(i)) then
                    aoduvst(i) = aoduvst(i) + dopaer(i)
                  end if
               end do
            end if

            if (savaernir) then
               do i = 1, ncol
                  extinctnir(i,k) = extinctnir(i,k) + dopaer(i)*air_density(i,k)/mass(i,k)
                  aodnir(i) = aodnir(i) + dopaer(i)
                  if (k.le.troplevchem(i)) then
                    aodnirst(i) = aodnirst(i) + dopaer(i)
                  end if
               end do
            endif

            ! Save aerosol optical depth at longest visible wavelength
            ! sum over layers
            if (savaervis) then
               ! aerosol extinction (/m)
               do i = 1, ncol
                  extinct(i,k) = extinct(i,k) + dopaer(i)*air_density(i,k)/mass(i,k)
                  absorb(i,k)  = absorb(i,k) + pabs(i)*air_density(i,k)
                  aodvis(i)    = aodvis(i) + dopaer(i)
                  aodabs(i)    = aodabs(i) + pabs(i)*mass(i,k)
                  aodmode(i)   = aodmode(i) + dopaer(i)
                  ssavis(i)    = ssavis(i) + dopaer(i)*palb(i)
                  asymvis(i)    = asymvis(i) + dopaer(i)*pasm(i)
                  asymext(i,k)  = asymext(i,k) + dopaer(i)*pasm(i)*air_density(i,k)/mass(i,k)
                  if (k.le.troplevchem(i)) then
                    aodvisst(i) = aodvisst(i) + dopaer(i)
                  end if

                  if (wetvol(i) > 1.e-40_r8) then

                     dustaodmode(i) = dustaodmode(i) + dopaer(i)*dustvol(i)/wetvol(i)

                     ! ++ Longlei Li
                     dustextmode(i,k) = dustextmode(i,k) + dopaer(i)*air_density(i,k)/mass(i,k)*dustvol(i)/wetvol(i)
                     ! -- Longlei Li

                     ! partition optical depth into contributions from each constituent
                     ! assume contribution is proportional to refractive index X volume

                     scath2o        = watervol(i)*real(crefwsw(isw))
		     absh2o         = -watervol(i)*aimag(crefwsw(isw))
                     ! Longlei Li: ash added
		     sumscat        = scatso4(i) + scatpom(i) + scatsoa(i) + scatbc(i) + &
                                      scatdust(i) + scatseasalt(i) + scath2o + scatfe(i) + scatash(i) 
		     sumabs         = absso4(i) + abspom(i) + abssoa(i) + absbc(i) + &
                                      absdust(i) + absseasalt(i) + absh2o + absfe(i) + absash(i)
                     sumhygro       = hygroso4(i) + hygropom(i) + hygrosoa(i) + hygrobc(i) + &
                                      hygrodust(i) + hygroseasalt(i) + hygrofe(i) + hygroash(i)

                     scatdust(i)    = (scatdust(i) + scath2o*hygrodust(i)/sumhygro)/sumscat
                     scatdust1(i)    = (scatdust1(i) + scath2o*hygrodust1(i)/sumhygro)/sumscat
                     scatdust2(i)    = (scatdust2(i) + scath2o*hygrodust2(i)/sumhygro)/sumscat
                     scatdust3(i)    = (scatdust3(i) + scath2o*hygrodust3(i)/sumhygro)/sumscat
                     scatdust4(i)    = (scatdust4(i) + scath2o*hygrodust4(i)/sumhygro)/sumscat
                     scatdust5(i)    = (scatdust5(i) + scath2o*hygrodust5(i)/sumhygro)/sumscat
                     scatdust6(i)    = (scatdust6(i) + scath2o*hygrodust6(i)/sumhygro)/sumscat
                     scatdust7(i)    = (scatdust7(i) + scath2o*hygrodust7(i)/sumhygro)/sumscat
                     scatdust8(i)    = (scatdust8(i) + scath2o*hygrodust8(i)/sumhygro)/sumscat
                     absdust(i)     = (absdust(i) + absh2o*hygrodust(i)/sumhygro)/sumabs
                     absdust1(i)     = (absdust1(i) + absh2o*hygrodust1(i)/sumhygro)/sumabs
                     absdust2(i)     = (absdust2(i) + absh2o*hygrodust2(i)/sumhygro)/sumabs
                     absdust3(i)     = (absdust3(i) + absh2o*hygrodust3(i)/sumhygro)/sumabs
                     absdust4(i)     = (absdust4(i) + absh2o*hygrodust4(i)/sumhygro)/sumabs
                     absdust5(i)     = (absdust5(i) + absh2o*hygrodust5(i)/sumhygro)/sumabs
                     absdust6(i)     = (absdust6(i) + absh2o*hygrodust6(i)/sumhygro)/sumabs
                     absdust7(i)     = (absdust7(i) + absh2o*hygrodust7(i)/sumhygro)/sumabs
                     absdust8(i)     = (absdust8(i) + absh2o*hygrodust8(i)/sumhygro)/sumabs

                     !++ Longlei Li ++
                     scatash(i)    = (scatash(i) + scath2o*hygroash(i)/sumhygro)/sumscat
                     !-- Longlei Li --
                     scatso4(i)     = (scatso4(i) + scath2o*hygroso4(i)/sumhygro)/sumscat
                     absso4(i)      = (absso4(i) + absh2o*hygroso4(i)/sumhygro)/sumabs

                     scatpom(i)     = (scatpom(i) + scath2o*hygropom(i)/sumhygro)/sumscat
                     abspom(i)      = (abspom(i) + absh2o*hygropom(i)/sumhygro)/sumabs

                     scatsoa(i)     = (scatsoa(i) + scath2o*hygrosoa(i)/sumhygro)/sumscat
                     abssoa(i)      = (abssoa(i) + absh2o*hygrosoa(i)/sumhygro)/sumabs

                     scatbc(i)      = (scatbc(i) + scath2o*hygrobc(i)/sumhygro)/sumscat
                     absbc(i)       = (absbc(i) + absh2o*hygrobc(i)/sumhygro)/sumabs

                     scatseasalt(i) = (scatseasalt(i) + scath2o*hygroseasalt(i)/sumhygro)/sumscat
                     absseasalt(i)  = (absseasalt(i) + absh2o*hygroseasalt(i)/sumhygro)/sumabs
                     
                     aodabsbc(i)    = aodabsbc(i) + absbc(i)*dopaer(i)*(1.0_r8-palb(i))

                     aodc           = (absdust(i)*(1.0_r8 - palb(i)) + palb(i)*scatdust(i))*dopaer(i)
                     ! ++Longlei Li++
                     aodcash        = (absash(i)*(1.0_r8 - palb(i)) + palb(i)*scatash(i))*dopaer(i)
                     ! ++Longlei Li--
                     aodc1           = (absdust1(i)*(1.0_r8 - palb(i)) + palb(i)*scatdust1(i))*dopaer(i)
                     aodc2           = (absdust2(i)*(1.0_r8 - palb(i)) + palb(i)*scatdust2(i))*dopaer(i)
                     aodc3           = (absdust3(i)*(1.0_r8 - palb(i)) + palb(i)*scatdust3(i))*dopaer(i)
                     aodc4           = (absdust4(i)*(1.0_r8 - palb(i)) + palb(i)*scatdust4(i))*dopaer(i)
                     aodc5           = (absdust5(i)*(1.0_r8 - palb(i)) + palb(i)*scatdust5(i))*dopaer(i)
                     aodc6           = (absdust6(i)*(1.0_r8 - palb(i)) + palb(i)*scatdust6(i))*dopaer(i)
                     aodc7           = (absdust7(i)*(1.0_r8 - palb(i)) + palb(i)*scatdust7(i))*dopaer(i)
                     aodc8           = (absdust8(i)*(1.0_r8 - palb(i)) + palb(i)*scatdust8(i))*dopaer(i)
                     !dustaod(i)     = dustaod(i) + aodc
                     ! ++Longlei Li
                     ! Account for dust asphericity (communication with JFK, Kok et al. 2017):
                     ! the extinction enhancement by 16% and 28% for the first two modes and the coarse mode, respectively
                     if (m == 1 .or. m == 2) then
                       dustaod(i)     = dustaod(i) + (aodc * 1.16_r8)
                       !dustaod1(i)     = dustaod1(i) + (aodc1 * 1.16_r8)
                       !dustaod2(i)     = dustaod2(i) + (aodc2 * 1.16_r8)
                       !dustaod3(i)     = dustaod3(i) + (aodc3 * 1.16_r8)
                       !dustaod4(i)     = dustaod4(i) + (aodc4 * 1.16_r8)
                       !dustaod5(i)     = dustaod5(i) + (aodc5 * 1.16_r8)
                       !dustaod6(i)     = dustaod6(i) + (aodc6 * 1.16_r8)
                       !dustaod7(i)     = dustaod7(i) + (aodc7 * 1.16_r8)
                       !dustaod8(i)     = dustaod8(i) + (aodc8 * 1.16_r8)

                       ashaod(i)       = ashaod(i) + (aodcash * 1.16_r8)
                       !ashaod1(i)      = ashaod1(i) + (aodc1 * 1.16_r8)
                       !ashaod2(i)      = ashaod2(i) + (aodc2 * 1.16_r8)
                       !ashaod3(i)      = ashaod3(i) + (aodc3 * 1.16_r8)
                       !ashaod4(i)      = ashaod4(i) + (aodc4 * 1.16_r8)

                       extdust(i,k)   = extdust(i,k) + aodc* 1.16_r8*air_density(i,k)/mass(i,k)
                       extdust1(i,k)   = extdust1(i,k) + aodc1* 1.16_r8*air_density(i,k)/mass(i,k)
                       extdust2(i,k)   = extdust2(i,k) + aodc2* 1.16_r8*air_density(i,k)/mass(i,k)
                       extdust3(i,k)   = extdust3(i,k) + aodc3* 1.16_r8*air_density(i,k)/mass(i,k)
                       extdust4(i,k)   = extdust4(i,k) + aodc4* 1.16_r8*air_density(i,k)/mass(i,k)
                       extdust5(i,k)   = extdust5(i,k) + aodc5* 1.16_r8*air_density(i,k)/mass(i,k)
                       extdust6(i,k)   = extdust6(i,k) + aodc6* 1.16_r8*air_density(i,k)/mass(i,k)
                       extdust7(i,k)   = extdust7(i,k) + aodc7* 1.16_r8*air_density(i,k)/mass(i,k)
                       extdust8(i,k)   = extdust8(i,k) + aodc8* 1.16_r8*air_density(i,k)/mass(i,k)
                     else
                       dustaod(i)     = dustaod(i) + (aodc * 1.28_r8)
                       extdust(i,k)   = extdust(i,k) + aodc* 1.28_r8*air_density(i,k)/mass(i,k)
                       !dustaod1(i)     = dustaod1(i) + (aodc1 * 1.28_r8)
                       !extdust1(i,k)   = extdust1(i,k) + aodc1* 1.28_r8*air_density(i,k)/mass(i,k)
                       !dustaod2(i)     = dustaod2(i) + (aodc2 * 1.28_r8)
                       !extdust2(i,k)   = extdust2(i,k) + aodc2* 1.28_r8*air_density(i,k)/mass(i,k)
                       !dustaod3(i)     = dustaod3(i) + (aodc3 * 1.28_r8)
                       !extdust3(i,k)   = extdust3(i,k) + aodc3* 1.28_r8*air_density(i,k)/mass(i,k)
                       !dustaod4(i)     = dustaod4(i) + (aodc4 * 1.28_r8)
                       !extdust4(i,k)   = extdust4(i,k) + aodc4* 1.28_r8*air_density(i,k)/mass(i,k)
                       !dustaod5(i)     = dustaod5(i) + (aodc5 * 1.28_r8)
                       !extdust5(i,k)   = extdust5(i,k) + aodc5* 1.28_r8*air_density(i,k)/mass(i,k)
                       !dustaod6(i)     = dustaod6(i) + (aodc6 * 1.28_r8)
                       !extdust6(i,k)   = extdust6(i,k) + aodc6* 1.28_r8*air_density(i,k)/mass(i,k)
                       !dustaod7(i)     = dustaod7(i) + (aodc7 * 1.28_r8)
                       !extdust7(i,k)   = extdust7(i,k) + aodc7* 1.28_r8*air_density(i,k)/mass(i,k)
                       !dustaod8(i)     = dustaod8(i) + (aodc8 * 1.28_r8)
                       !extdust8(i,k)   = extdust8(i,k) + aodc8* 1.28_r8*air_density(i,k)/mass(i,k)

                       ashaod(i)     = ashaod(i) + (aodcash * 1.28_r8)
                       !ashaod1(i)     = dustaod1(i) + (aodc1 * 1.28_r8)
                       !ashaod2(i)     = dustaod2(i) + (aodc2 * 1.28_r8)
                       !ashaod3(i)     = dustaod3(i) + (aodc3 * 1.28_r8)
                       !ashaod4(i)     = dustaod4(i) + (aodc4 * 1.28_r8)
                     end if
                     ! --Longlei Li
   
                     aodc           = (absso4(i)*(1.0_r8 - palb(i)) + palb(i)*scatso4(i))*dopaer(i)
                     so4aod(i)      = so4aod(i) + aodc

                     aodc           = (abspom(i)*(1.0_r8 - palb(i)) + palb(i)*scatpom(i))*dopaer(i)
                     pomaod(i)      = pomaod(i) + aodc

                     aodc           = (abssoa(i)*(1.0_r8 - palb(i)) + palb(i)*scatsoa(i))*dopaer(i)
                     soaaod(i)      = soaaod(i) + aodc

                     aodc           = (absbc(i)*(1.0_r8 - palb(i)) + palb(i)*scatbc(i))*dopaer(i)
                     bcaod(i)       = bcaod(i) + aodc

                     aodc           = (absseasalt(i)*(1.0_r8 - palb(i)) + palb(i)*scatseasalt(i))*dopaer(i)
                     seasaltaod(i)  = seasaltaod(i) + aodc

                  endif

               end do
            endif

            do i = 1, ncol

               if ((dopaer(i) <= -1.e-10_r8) .or. (dopaer(i) >= 30._r8)) then

                  if (dopaer(i) <= -1.e-10_r8) then
                     write(iulog,*) "ERROR: Negative aerosol optical depth &
                          &in this layer."
                  else
                     write(iulog,*) "WARNING: Aerosol optical depth is &
                          &unreasonably high in this layer at shortwave."
                  end if

                  write(iulog,*) 'dopaer(', i, ',', k, ',', m, ',', lchnk, ')=', dopaer(i)
                  ! write(iulog,*) 'itab,jtab,ttab,utab=',itab(i),jtab(i),ttab(i),utab(i)
                  write(iulog,*) 'k=', k, ' pext=', pext(i), ' specext=', specpext(i)
                  write(iulog,*) 'wetvol=', wetvol(i), ' dryvol=', dryvol(i), ' watervol=', watervol(i)
                  ! write(iulog,*) 'cext=',(cext(i,l),l=1,ncoef)
                  ! write(iulog,*) 'crefin=',crefin(i)
                  write(iulog,'(a,i12,a)') 'nspec=,  specname=  ', nspec, spectype
                  ! write(iulog,*) 'cheb=', (cheb(nc,m,i,k),nc=2,ncoef)
                  do l = 1, nspec
                     call rad_cnst_get_aer_mmr(list_idx, m, l, 'a', state, pbuf, specmmr)
                     call rad_cnst_get_aer_props(list_idx, m, l, density_aer=specdens, &
                                                 refindex_aer_sw=specrefindex)
                     volf = specmmr(i,k)/specdens
                     write(iulog,*) 'l=', l, 'vol(l)=', volf
                     write(iulog,*) 'isw=', isw, 'specrefindex(isw)=', specrefindex(isw)
                     write(iulog,*) 'specdens=', specdens
                  end do

                  nerr_dopaer = nerr_dopaer + 1
!                  if (nerr_dopaer >= nerrmax_dopaer) then
                  if (dopaer(i) < -1.e-10_r8) then
                  !if (nerr_dopaer >= nerrmax_dopaer .and. dopaer(i) < -1.e-10_r8) then
                     write(iulog,*) '*** halting in '//subname//' after nerr_dopaer =', nerr_dopaer
                     call endrun('exit from '//subname)
                  end if

               end if
            end do

            do i=1,ncol
               !tauxar(i,k,isw) = tauxar(i,k,isw) + dopaer(i)
               !wa(i,k,isw)     = wa(i,k,isw)     + dopaer(i)*palb(i)
               !ga(i,k,isw)     = ga(i,k,isw)     + dopaer(i)*palb(i)*pasm(i)
               !fa(i,k,isw)     = fa(i,k,isw)     + dopaer(i)*palb(i)*pasm(i)*pasm(i)
               ! ++Longlei Li
               if (savaervis .and. wetvol(i) > 1.e-40_r8) then

                  if (m == 1 .or. m == 2) then
                     eaod           = ((absdust(i)*(1.0_r8 - palb(i)) + palb(i)*scatdust(i))*0.16_r8+1.00_r8)*dopaer(i)
                  else if (m == 3) then
                     eaod           = ((absdust(i)*(1.0_r8 - palb(i)) + palb(i)*scatdust(i))*0.28_r8+1.00_r8)*dopaer(i)
                  else
                     eaod           = dopaer(i)
                  end if

                  tauxar(i,k,isw) = tauxar(i,k,isw) + eaod
                  wa(i,k,isw)     = wa(i,k,isw)     + eaod*palb(i)
                  ga(i,k,isw)     = ga(i,k,isw)     + eaod*palb(i)*pasm(i)
                  fa(i,k,isw)     = fa(i,k,isw)     + eaod*palb(i)*pasm(i)*pasm(i)
               else
                  tauxar(i,k,isw) = tauxar(i,k,isw) + dopaer(i)
                  wa(i,k,isw)     = wa(i,k,isw)     + dopaer(i)*palb(i)
                  ga(i,k,isw)     = ga(i,k,isw)     + dopaer(i)*palb(i)*pasm(i)
                  fa(i,k,isw)     = fa(i,k,isw)     + dopaer(i)*palb(i)*pasm(i)*pasm(i)
               end if
               ! --Longlei Li
            end do

         end do ! pver

      end do ! sw bands

      ! mode diagnostics
      ! The diagnostics are currently only output for the climate list.  Code mods will
      ! be necessary to provide output for the rad_diag lists.
      if (list_idx == 0) then

         write(outname,'(a,i1)') 'BURDENdn', m
         call outfld(trim(outname), burden, pcols, lchnk)

         write(outname,'(a,i1)') 'AODdnMODE', m
         call outfld(trim(outname), aodmode, pcols, lchnk)

         write(outname,'(a,i1)') 'AODdnDUST', m
         call outfld(trim(outname), dustaodmode, pcols, lchnk)
         
         do i = 1, nnite
            burden(idxnite(i))  = fillvalue
            aodmode(idxnite(i)) = fillvalue
            dustaodmode(idxnite(i)) = fillvalue
         end do

         write(outname,'(a,i1)') 'BURDEN', m
         call outfld(trim(outname), burden, pcols, lchnk)

         write(outname,'(a,i1)') 'AODMODE', m
         call outfld(trim(outname), aodmode, pcols, lchnk)

         write(outname,'(a,i1)') 'AODDUST', m
         call outfld(trim(outname), dustaodmode, pcols, lchnk)

         !++ Longlei Li
         write(outname,'(a,i1)') 'EXTDUST', m
         call outfld(trim(outname), dustextmode, pcols, lchnk)
         !-- Longlei Li

      end if

   end do ! nmodes

   if (list_idx > 0) then
      deallocate(dgnumdry_m)
      deallocate(dgnumwet_m)
      deallocate(qaerwat_m)
      deallocate(wetdens_m)
      deallocate(hygro_m)
      deallocate(dryvol_m)
      deallocate(dryrad_m)
      deallocate(drymass_m)
      deallocate(so4dryvol_m)
      deallocate(naer_m)
   end if

   ! Output visible band diagnostics for quantities summed over the modes
   ! These fields are put out for diagnostic lists as well as the climate list.

   call outfld('EXTINCTdn'//diag(list_idx),  extinct, pcols, lchnk)
   !++ Longlei Li
   call outfld('EXTDUSTdn'//diag(list_idx),  extinct, pcols, lchnk)
   !-- Longlei Li

   call outfld('ABSORBdn'//diag(list_idx),   absorb,  pcols, lchnk)
   call outfld('AODVISdn'//diag(list_idx),   aodvis,  pcols, lchnk)
   call outfld('AODABSdn'//diag(list_idx),   aodabs,  pcols, lchnk)
   call outfld('AODVISstdn'//diag(list_idx), aodvisst,pcols, lchnk)
   call outfld('EXTxASYMdn'//diag(list_idx), asymext, pcols, lchnk)
   
   do i = 1, nnite
      extinct(idxnite(i),:) = fillvalue
      absorb(idxnite(i),:)  = fillvalue
      aodvis(idxnite(i))    = fillvalue
      aodabs(idxnite(i))    = fillvalue
      aodvisst(idxnite(i))  = fillvalue
      asymext(idxnite(i),:) = fillvalue
   end do

   call outfld('EXTINCT'//diag(list_idx),  extinct, pcols, lchnk)
   !++ Longlei Li
   call outfld('EXTDUST'//diag(list_idx),  extinct, pcols, lchnk)
   !-- Longlei Li
   call outfld('ABSORB'//diag(list_idx),   absorb,  pcols, lchnk)
   call outfld('AODVIS'//diag(list_idx),   aodvis,  pcols, lchnk)
   call outfld('AODABS'//diag(list_idx),   aodabs,  pcols, lchnk)
   call outfld('AODVISst'//diag(list_idx), aodvisst,pcols, lchnk)
   call outfld('EXTxASYM'//diag(list_idx), asymext, pcols, lchnk)

   ! These diagnostics are output only for climate list
   if (list_idx == 0) then
      do i = 1, ncol
         if (aodvis(i) > 1.e-10_r8) then
            ssavis(i) = ssavis(i)/aodvis(i)
         else
            ssavis(i) = 0.925_r8
         endif
      end do
      
      call outfld('SSAVISdn',        ssavis,        pcols, lchnk)
      call outfld('AODxASYMdn',      asymvis,       pcols, lchnk)

      call outfld('EXTINCTUVdn',     extinctuv,     pcols, lchnk)
      call outfld('EXTINCTNIRdn',    extinctnir,    pcols, lchnk)
      call outfld('AODUVdn',         aoduv,         pcols, lchnk)
      call outfld('AODNIRdn',        aodnir,        pcols, lchnk)
      call outfld('AODUVstdn',       aoduvst,       pcols, lchnk)
      call outfld('AODNIRstdn',      aodnirst,      pcols, lchnk)

      call outfld('BURDENDUSTdn',    burdendust,    pcols, lchnk)
      !call outfld('BURDENDUST1dn',    burdendust1,    pcols, lchnk)
      !call outfld('BURDENDUST2dn',    burdendust2,    pcols, lchnk)
      !call outfld('BURDENDUST3dn',    burdendust3,    pcols, lchnk)
      !call outfld('BURDENDUST4dn',    burdendust4,    pcols, lchnk)
      !call outfld('BURDENDUST5dn',    burdendust5,    pcols, lchnk)
      !call outfld('BURDENDUST6dn',    burdendust6,    pcols, lchnk)
      !call outfld('BURDENDUST7dn',    burdendust7,    pcols, lchnk)
      !call outfld('BURDENDUST8dn',    burdendust8,    pcols, lchnk)
      call outfld('BURDENSO4dn' ,    burdenso4,     pcols, lchnk)
      call outfld('BURDENPOMdn' ,    burdenpom,     pcols, lchnk)
      call outfld('BURDENSOAdn' ,    burdensoa,     pcols, lchnk)
      call outfld('BURDENBCdn'  ,    burdenbc,      pcols, lchnk)
      call outfld('BURDENSEASALTdn', burdenseasalt, pcols, lchnk)

      call outfld('AODABSBCdn',      aodabsbc,      pcols, lchnk)

      call outfld('AODDUSTdn',       dustaod,       pcols, lchnk)
      !call outfld('AODDUSTS1dn',       dustaod1,       pcols, lchnk)
      !call outfld('AODDUSTS2dn',       dustaod2,       pcols, lchnk)
      !call outfld('AODDUSTS3dn',       dustaod3,       pcols, lchnk)
      !call outfld('AODDUSTS4dn',       dustaod4,       pcols, lchnk)
      !call outfld('AODDUSTS5dn',       dustaod5,       pcols, lchnk)
      !call outfld('AODDUSTS6dn',       dustaod6,       pcols, lchnk)
      !call outfld('AODDUSTS7dn',       dustaod7,       pcols, lchnk)
      !call outfld('AODDUSTS8dn',       dustaod8,       pcols, lchnk)
      call outfld('AODSO4dn',        so4aod,        pcols, lchnk)
      call outfld('AODPOMdn',        pomaod,        pcols, lchnk)
      call outfld('AODSOAdn',        soaaod,        pcols, lchnk)
      call outfld('AODBCdn',         bcaod,         pcols, lchnk)
      call outfld('AODSSdn',         seasaltaod,    pcols, lchnk)


      do i = 1, nnite
         ssavis(idxnite(i))     = fillvalue
         asymvis(idxnite(i))    = fillvalue

         aoduv(idxnite(i))      = fillvalue
         aodnir(idxnite(i))     = fillvalue
         aoduvst(idxnite(i))    = fillvalue
         aodnirst(idxnite(i))   = fillvalue
         extinctuv(idxnite(i),:)  = fillvalue
         extinctnir(idxnite(i),:) = fillvalue

         burdendust(idxnite(i)) = fillvalue
         ! ++ Longlei Li ++
         burdenash(idxnite(i)) = fillvalue
         ! -- Longlei Li --
         burdendust1(idxnite(i)) = fillvalue
         burdendust2(idxnite(i)) = fillvalue
         burdendust3(idxnite(i)) = fillvalue
         burdendust4(idxnite(i)) = fillvalue
         burdendust5(idxnite(i)) = fillvalue
         burdendust6(idxnite(i)) = fillvalue
         burdendust7(idxnite(i)) = fillvalue
         burdendust8(idxnite(i)) = fillvalue
         burdenso4(idxnite(i))  = fillvalue
         burdenpom(idxnite(i))  = fillvalue
         burdensoa(idxnite(i))  = fillvalue
         burdenbc(idxnite(i))   = fillvalue
         burdenseasalt(idxnite(i)) = fillvalue

         aodabsbc(idxnite(i))   = fillvalue

         dustaod(idxnite(i))    = fillvalue
         ! ++ Longlei Li ++
         ashaod(idxnite(i))    = fillvalue
         ! -- Longlei Li --

         !dustaod1(idxnite(i))    = fillvalue
         !dustaod2(idxnite(i))    = fillvalue
         !dustaod3(idxnite(i))    = fillvalue
         !dustaod4(idxnite(i))    = fillvalue
         !dustaod5(idxnite(i))    = fillvalue
         !dustaod6(idxnite(i))    = fillvalue
         !dustaod7(idxnite(i))    = fillvalue
         !dustaod8(idxnite(i))    = fillvalue
         so4aod(idxnite(i))     = fillvalue
         pomaod(idxnite(i))     = fillvalue
         soaaod(idxnite(i))     = fillvalue
         bcaod(idxnite(i))      = fillvalue
         seasaltaod(idxnite(i)) = fillvalue
       end do

      call outfld('SSAVIS',        ssavis,        pcols, lchnk)
      call outfld('AODxASYM',      asymvis,       pcols, lchnk)

      call outfld('EXTINCTUV',     extinctuv,     pcols, lchnk)
      call outfld('EXTINCTNIR',    extinctnir,    pcols, lchnk)
      call outfld('AODUV',         aoduv,         pcols, lchnk)
      call outfld('AODNIR',        aodnir,        pcols, lchnk)
      call outfld('AODUVst',       aoduvst,       pcols, lchnk)
      call outfld('AODNIRst',      aodnirst,      pcols, lchnk)

      call outfld('BURDENDUST',    burdendust,    pcols, lchnk)
      !call outfld('BURDENDUST1',    burdendust1,    pcols, lchnk)
      !call outfld('BURDENDUST2',    burdendust2,    pcols, lchnk)
      !call outfld('BURDENDUST3',    burdendust3,    pcols, lchnk)
      !call outfld('BURDENDUST4',    burdendust4,    pcols, lchnk)
      !call outfld('BURDENDUST5',    burdendust5,    pcols, lchnk)
      !call outfld('BURDENDUST6',    burdendust6,    pcols, lchnk)
      !call outfld('BURDENDUST7',    burdendust7,    pcols, lchnk)
      !call outfld('BURDENDUST8',    burdendust8,    pcols, lchnk)
      call outfld('BURDENSO4' ,    burdenso4,     pcols, lchnk)
      call outfld('BURDENPOM' ,    burdenpom,     pcols, lchnk)
      call outfld('BURDENSOA' ,    burdensoa,     pcols, lchnk)
      call outfld('BURDENBC'  ,    burdenbc,      pcols, lchnk)
      call outfld('BURDENSEASALT', burdenseasalt, pcols, lchnk)

      ! ++ Longlei Li ++ 
      call outfld('BURDENASH',    burdenash,    pcols, lchnk)
      ! -- Longlei Li --

      call outfld('AODABSBC',      aodabsbc,      pcols, lchnk)

      call outfld('AODDUST',       dustaod,       pcols, lchnk)

      ! ++ Longlei Li ++ 
      call outfld('AODASH',       ashaod,       pcols, lchnk)
      ! -- Longlei Li --

      !call outfld('AODDUSTS1',       dustaod1,       pcols, lchnk)
      !call outfld('AODDUSTS2',       dustaod2,       pcols, lchnk)
      !call outfld('AODDUSTS3',       dustaod3,       pcols, lchnk)
      !call outfld('AODDUSTS4',       dustaod4,       pcols, lchnk)
      !call outfld('AODDUSTS5',       dustaod5,       pcols, lchnk)
      !call outfld('AODDUSTS6',       dustaod6,       pcols, lchnk)
      !call outfld('AODDUSTS7',       dustaod7,       pcols, lchnk)
      !call outfld('AODDUSTS8',       dustaod8,       pcols, lchnk)
      call outfld('AODSO4',        so4aod,        pcols, lchnk)
      call outfld('AODPOM',        pomaod,        pcols, lchnk)
      call outfld('AODSOA',        soaaod,        pcols, lchnk)
      call outfld('AODBC',         bcaod,         pcols, lchnk)
      call outfld('AODSS',         seasaltaod,    pcols, lchnk)
   end if

    !write(iulog,*) 'Longlei marks end in modal_aer_opt'
end subroutine modal_aero_sw

!===============================================================================

subroutine modal_aero_lw(list_idx, state, pbuf, tauxar)

   ! calculates aerosol lw radiative properties

   integer,             intent(in)  :: list_idx ! index of the climate or a diagnostic list
   type(physics_state), intent(in), target :: state    ! state variables
   
   type(physics_buffer_desc), pointer :: pbuf(:)

   real(r8), intent(out) :: tauxar(pcols,pver,nlwbands) ! layer absorption optical depth

   ! Local variables
   integer :: i, ifld, ilw, k, l, m, nc, ns
   integer :: lchnk                    ! chunk id
   integer :: ncol                     ! number of active columns in the chunk
   integer :: nmodes
   integer :: nspec
   integer :: istat

   real(r8), pointer :: dgnumwet(:,:)  ! wet number mode diameter (m)
   real(r8), pointer :: qaerwat(:,:)   ! aerosol water (g/g)

   real(r8), pointer :: dgnumdry_m(:,:,:) ! number mode dry diameter for all modes
   real(r8), pointer :: dgnumwet_m(:,:,:) ! number mode wet diameter for all modes
   real(r8), pointer :: qaerwat_m(:,:,:)  ! aerosol water (g/g) for all modes
   real(r8), pointer :: wetdens_m(:,:,:)  ! 
   real(r8), pointer :: hygro_m(:,:,:)  !
   real(r8), pointer :: dryvol_m(:,:,:)  !
   real(r8), pointer :: dryrad_m(:,:,:)  !
   real(r8), pointer :: drymass_m(:,:,:)  !
   real(r8), pointer :: so4dryvol_m(:,:,:)  !
   real(r8), pointer :: naer_m(:,:,:)  !

   real(r8) :: sigma_logr_aer          ! geometric standard deviation of number distribution
   real(r8) :: alnsg_amode             ! log of geometric standard deviation of number distribution
   real(r8) :: xrad(pcols)
   real(r8) :: cheby(ncoef,pcols,pver)  ! chebychef polynomials

   real(r8) :: mass(pcols,pver) ! layer mass

   real(r8),    pointer :: specmmr(:,:)        ! species mass mixing ratio
   real(r8)             :: specdens            ! species density (kg/m3)
   complex(r8), pointer :: specrefindex(:)     ! species refractive index
!++Longlei Li
   character*32         :: spectype            ! species type
!--

   real(r8) :: vol(pcols)       ! volume concentration of aerosol specie (m3/kg)
   real(r8) :: dryvol(pcols)    ! volume concentration of aerosol mode (m3/kg)
   real(r8) :: wetvol(pcols)    ! volume concentration of wet mode (m3/kg)
   real(r8) :: watervol(pcols)  ! volume concentration of water in each mode (m3/kg)
   real(r8) :: refr(pcols)      ! real part of refractive index
   real(r8) :: refi(pcols)      ! imaginary part of refractive index
   complex(r8) :: crefin(pcols) ! complex refractive index
   real(r8), pointer :: refrtablw(:,:) ! table of real refractive indices for aerosols
   real(r8), pointer :: refitablw(:,:) ! table of imag refractive indices for aerosols
   real(r8), pointer :: absplw(:,:,:,:) ! specific absorption

   integer  :: itab(pcols), jtab(pcols)
   real(r8) :: ttab(pcols), utab(pcols)
   real(r8) :: cabs(pcols,ncoef)
   real(r8) :: pabs(pcols)      ! parameterized specific absorption (m2/kg)
   real(r8) :: dopaer(pcols)    ! aerosol optical depth in layer

   integer, parameter :: nerrmax_dopaer=1000
   integer  :: nerr_dopaer = 0
   real(r8) :: volf             ! volume fraction of insoluble aerosol

   character(len=*), parameter :: subname = 'modal_aero_lw'
   !----------------------------------------------------------------------------

   lchnk = state%lchnk
   ncol  = state%ncol

   ! initialize output variables
   tauxar(:ncol,:,:) = 0._r8

   ! dry mass in each cell
   mass(:ncol,:) = state%pdeldry(:ncol,:)*rga

   ! loop over all aerosol modes
   call rad_cnst_get_info(list_idx, nmodes=nmodes)

   if (list_idx == 0) then
      ! water uptake and wet radius for the climate list has already been calculated
      call pbuf_get_field(pbuf, dgnumwet_idx, dgnumwet_m)
      call pbuf_get_field(pbuf, qaerwat_idx,  qaerwat_m)
   else
      ! If doing a diagnostic calculation then need to calculate the wet radius
      ! and water uptake for the diagnostic modes
      allocate(dgnumdry_m(pcols,pver,nmodes),  dgnumwet_m(pcols,pver,nmodes), &
               qaerwat_m(pcols,pver,nmodes),   wetdens_m(pcols,pver,nmodes), &
               hygro_m(pcols,pver,nmodes),     dryvol_m(pcols,pver,nmodes), &
               dryrad_m(pcols,pver,nmodes),    drymass_m(pcols,pver,nmodes),  &
               so4dryvol_m(pcols,pver,nmodes), naer_m(pcols,pver,nmodes),     stat=istat)

      if (istat > 0) then
         call endrun('modal_aero_lw: allocation FAILURE: arrays for diagnostic calcs')
      end if
      call modal_aero_calcsize_diag(state, pbuf, list_idx, dgnumdry_m, hygro_m, &
                                    dryvol_m, dryrad_m, drymass_m, so4dryvol_m, naer_m)  
      call modal_aero_wateruptake_dr(state, pbuf, list_idx, dgnumdry_m, dgnumwet_m, &
                                     qaerwat_m, wetdens_m,  hygro_m, dryvol_m, dryrad_m, &
                                     drymass_m, so4dryvol_m, naer_m)
   endif

   do m = 1, nmodes

      dgnumwet => dgnumwet_m(:,:,m)
      qaerwat  => qaerwat_m(:,:,m)

      ! get mode properties
      call rad_cnst_get_mode_props(list_idx, m, sigmag=sigma_logr_aer, refrtablw=refrtablw , &
         refitablw=refitablw, absplw=absplw)

      ! get mode info
      call rad_cnst_get_info(list_idx, m, nspec=nspec)

      ! calc size parameter for all columns
      ! this is the same calculation that's done in modal_size_parameters, but there
      ! some intermediate results are saved and the chebyshev polynomials are stored
      ! in a array with different index order.  Could be unified.
      do k = top_lev, pver
         do i = 1, ncol
            alnsg_amode = log( sigma_logr_aer )
            ! convert from number diameter to surface area
            xrad(i) = log(0.5_r8*dgnumwet(i,k)) + 2.0_r8*alnsg_amode*alnsg_amode
            ! normalize size parameter
            xrad(i) = max(xrad(i), xrmin)
            xrad(i) = min(xrad(i), xrmax)
            xrad(i) = (2*xrad(i)-xrmax-xrmin)/(xrmax-xrmin)
            ! chebyshev polynomials
            cheby(1,i,k) = 1.0_r8
            cheby(2,i,k) = xrad(i)
            do nc = 3, ncoef
               cheby(nc,i,k) = 2.0_r8*xrad(i)*cheby(nc-1,i,k)-cheby(nc-2,i,k)
            end do
         end do
      end do

      do ilw = 1, nlwbands

         do k = top_lev, pver

            ! form bulk refractive index. Use volume mixing for infrared
            crefin(:ncol) = (0._r8, 0._r8)
            dryvol(:ncol) = 0._r8

            ! aerosol species loop
            do l = 1, nspec
               call rad_cnst_get_aer_mmr(list_idx, m, l, 'a', state, pbuf, specmmr)
               call rad_cnst_get_aer_props(list_idx, m, l, density_aer=specdens, &
                                           refindex_aer_lw=specrefindex, spectype=spectype)
                                                             ! Longlei called spectype here

!++Longlei
               !if ( (m == 3 .or. m == 1) .and. l == 1) then
                  !specmmr(1:ncol,top_lev:pver)  = 1.0e-36_r8
                  !qaerwat(1:ncol,top_lev:pver)  = 1.0e-36_r8
               !end if
!--
               do i = 1, ncol
                  vol(i)    = specmmr(i,k)/specdens
                  dryvol(i) = dryvol(i) + vol(i)
                  crefin(i) = crefin(i) + vol(i)*specrefindex(ilw)
               end do
            end do

            do i = 1, ncol
               watervol(i) = qaerwat(i,k)/rhoh2o
               wetvol(i)   = watervol(i) + dryvol(i)
               if (watervol(i) < 0.0_r8) then
                  if (abs(watervol(i)) .gt. 1.e-1_r8*wetvol(i)) then
                     write(iulog,*) 'watervol,wetvol,dryvol=',watervol(i),wetvol(i),dryvol(i),' in '//subname
                  end if
                  watervol(i) = 0._r8
                  wetvol(i)   = dryvol(i)
               end if

               crefin(i) = crefin(i) + watervol(i)*crefwlw(ilw)
               if (wetvol(i) > 1.e-40_r8) crefin(i) = crefin(i)/wetvol(i)
               refr(i) = real(crefin(i))
               refi(i) = aimag(crefin(i))
            end do

            ! interpolate coefficients linear in refractive index
            ! first call calcs itab,jtab,ttab,utab
            itab(:ncol) = 0
            call binterp(absplw(:,:,:,ilw), ncol, ncoef, prefr, prefi, &
                         refr, refi, refrtablw(:,ilw), refitablw(:,ilw), &
                         itab, jtab, ttab, utab, cabs)

            ! parameterized optical properties
            do i = 1, ncol
               pabs(i) = 0.5_r8*cabs(i,1)
               do nc = 2, ncoef
                  pabs(i) = pabs(i) + cheby(nc,i,k)*cabs(i,nc)
               end do
               pabs(i)   = pabs(i)*wetvol(i)*rhoh2o
               pabs(i)   = max(0._r8,pabs(i))
               dopaer(i) = pabs(i)*mass(i,k)
!++Longlei Li
               !if (dopaer(i) >= 20._r8) then
               !      write( iulog, '(2a,2x,1(1pe14.5), a, i12)' ) &
               !           'For high AOD, Longlei', &
               !           ' mass(i,k) = ,  pabs(i)= , wetvol', mass(i,k), pabs(i), wetvol(i)
               !end if
!--Longlei Li
            end do

            do i = 1, ncol

               if ((dopaer(i) <= -1.e-10_r8) .or. (dopaer(i) >= 20._r8)) then

                  if (dopaer(i) <= -1.e-10_r8) then
                     write(iulog,*) "ERROR: Negative aerosol optical depth &
                          &in this layer."
                  else
                     write(iulog,*) "WARNING: Aerosol optical depth is &
                          &unreasonably high in this layer at longwave."
                  end if

                  write(iulog,*) 'dopaer(',i,',',k,',',m,',',lchnk,')=', dopaer(i)
                  write(iulog,*) 'k=',k,' pabs=', pabs(i)
                  write(iulog,*) 'wetvol=',wetvol(i),' dryvol=',dryvol(i),     &
                     ' watervol=',watervol(i)
                  write(iulog,*) 'cabs=', (cabs(i,l),l=1,ncoef)
                  write(iulog,*) 'crefin=', crefin(i)
                  write(iulog,*) 'nspec=', nspec
                  do l = 1,nspec
                     call rad_cnst_get_aer_mmr(list_idx, m, l, 'a', state, pbuf, specmmr)
                     call rad_cnst_get_aer_props(list_idx, m, l, density_aer=specdens, &
                                                 refindex_aer_lw=specrefindex, spectype=spectype)
                                                                  ! Longlei called spectype here
                     volf = specmmr(i,k)/specdens
!++Longlei Li
                     write( iulog, '(2a,2x,1(1pe14.5), a, i12)' ) &
                          'For high AOD, Longlei', &
                          ' specmmr(i,k) = ,  spectype=, m= ', specmmr(i,k), spectype, m
!--Longlei Li
                     write(iulog,*) 'l=',l,'vol(l)=',volf
                     write(iulog,*) 'ilw=',ilw,' specrefindex(ilw)=',specrefindex(ilw)
                     write(iulog,*) 'specdens=',specdens
                  end do

                  nerr_dopaer = nerr_dopaer + 1
                  if (nerr_dopaer >= nerrmax_dopaer .or. dopaer(i) < -1.e-10_r8) then
                     write(iulog,*) '*** halting in '//subname//' after nerr_dopaer =', nerr_dopaer
                     call endrun()
                  end if

               end if
            end do

            do i = 1, ncol
               tauxar(i,k,ilw) = tauxar(i,k,ilw) + dopaer(i)
            end do

         end do ! k = top_lev, pver

      end do  ! nlwbands

   end do ! m = 1, nmodes

   if (list_idx > 0) then
      deallocate(dgnumdry_m)
      deallocate(dgnumwet_m)
      deallocate(qaerwat_m)
      deallocate(wetdens_m)
      deallocate(hygro_m)
      deallocate(dryvol_m)
      deallocate(dryrad_m)
      deallocate(drymass_m)
      deallocate(so4dryvol_m)
      deallocate(naer_m)
   end if

end subroutine modal_aero_lw

!===============================================================================
! Private routines
!===============================================================================

subroutine read_water_refindex(infilename)

   ! read water refractive index file and set module data

   character*(*), intent(in) :: infilename   ! modal optics filename

   ! Local variables

   integer            :: i, ierr
   type(file_desc_t)  :: ncid              ! pio file handle
   integer            :: did               ! dimension ids
   integer            :: dimlen            ! dimension lengths
   type(var_desc_t)   :: vid               ! variable ids
   real(r8) :: refrwsw(nswbands), refiwsw(nswbands) ! real, imaginary ref index for water visible
   real(r8) :: refrwlw(nlwbands), refiwlw(nlwbands) ! real, imaginary ref index for water infrared
   !----------------------------------------------------------------------------

   ! open file
   call cam_pio_openfile(ncid, infilename, PIO_NOWRITE)

   ! inquire dimensions.  Check that file values match parameter values.

   ierr = pio_inq_dimid(ncid, 'lw_band', did)
   ierr = pio_inq_dimlen(ncid, did, dimlen)
   if (dimlen .ne. nlwbands) then
      write(iulog,*) 'lw_band len=', dimlen, ' from ', infilename, ' ne nlwbands=', nlwbands
      call endrun('read_modal_optics: bad lw_band value')
   endif

   ierr = pio_inq_dimid(ncid, 'sw_band', did)
   ierr = pio_inq_dimlen(ncid, did, dimlen)
   if (dimlen .ne. nswbands) then
      write(iulog,*) 'sw_band len=', dimlen, ' from ', infilename, ' ne nswbands=', nswbands
      call endrun('read_modal_optics: bad sw_band value')
   endif

   ! read variables
   ierr = pio_inq_varid(ncid, 'refindex_real_water_sw', vid)
   ierr = pio_get_var(ncid, vid, refrwsw)

   ierr = pio_inq_varid(ncid, 'refindex_im_water_sw', vid)
   ierr = pio_get_var(ncid, vid, refiwsw)

   ierr = pio_inq_varid(ncid, 'refindex_real_water_lw', vid)
   ierr = pio_get_var(ncid, vid, refrwlw)

   ierr = pio_inq_varid(ncid, 'refindex_im_water_lw', vid)
   ierr = pio_get_var(ncid, vid, refiwlw)

   ! set complex representation of refractive indices as module data
   do i = 1, nswbands
      crefwsw(i)  = cmplx(refrwsw(i), abs(refiwsw(i)),kind=r8)
   end do
   do i = 1, nlwbands
      crefwlw(i)  = cmplx(refrwlw(i), abs(refiwlw(i)),kind=r8)
   end do

   call pio_closefile(ncid)

end subroutine read_water_refindex

!===============================================================================

subroutine modal_size_parameters(ncol, sigma_logr_aer, dgnumwet, radsurf, logradsurf, cheb)

   integer,  intent(in)  :: ncol
   real(r8), intent(in)  :: sigma_logr_aer  ! geometric standard deviation of number distribution
   real(r8), intent(in)  :: dgnumwet(:,:)   ! aerosol wet number mode diameter (m)
   real(r8), intent(out) :: radsurf(:,:)    ! aerosol surface mode radius
   real(r8), intent(out) :: logradsurf(:,:) ! log(aerosol surface mode radius)
   real(r8), intent(out) :: cheb(:,:,:)

   integer  :: i, k, nc
   real(r8) :: alnsg_amode
   real(r8) :: explnsigma
   real(r8) :: xrad(pcols) ! normalized aerosol radius
   !-------------------------------------------------------------------------------

   alnsg_amode = log(sigma_logr_aer)
   explnsigma = exp(2.0_r8*alnsg_amode*alnsg_amode)

   do k = top_lev, pver
      do i = 1, ncol
         ! convert from number mode diameter to surface area
         radsurf(i,k) = 0.5_r8*dgnumwet(i,k)*explnsigma
         logradsurf(i,k) = log(radsurf(i,k))
         ! normalize size parameter
         xrad(i) = max(logradsurf(i,k),xrmin)
         xrad(i) = min(xrad(i),xrmax)
         xrad(i) = (2._r8*xrad(i)-xrmax-xrmin)/(xrmax-xrmin)
         ! chebyshev polynomials
         cheb(1,i,k) = 1._r8
         cheb(2,i,k) = xrad(i)
         do nc = 3, ncoef
            cheb(nc,i,k) = 2._r8*xrad(i)*cheb(nc-1,i,k)-cheb(nc-2,i,k)
         end do
      end do
   end do

end subroutine modal_size_parameters

!===============================================================================

    subroutine binterp(table,ncol,km,im,jm,x,y,xtab,ytab,ix,jy,t,u,out)

        !     bilinear interpolation of table
        !
        implicit none
        integer im,jm,km,ncol
        real(r8) table(km,im,jm),xtab(im),ytab(jm),out(pcols,km)
        integer i,ix(pcols),ip1,j,jy(pcols),jp1,k,ic,ip1m(pcols),jp1m(pcols),ixc,jyc
        real(r8) x(pcols),dx,t(pcols),y(pcols),dy,u(pcols),tu(pcols),tuc(pcols),tcu(pcols),tcuc(pcols)
        
        if(ix(1).gt.0) go to 30
        if(im.gt.1)then
            do ic=1,ncol
                do i=1,im
                    if(x(ic).lt.xtab(i))go to 10
                enddo
                10 ix(ic)=max0(i-1,1)
                ip1=min(ix(ic)+1,im)
                dx=(xtab(ip1)-xtab(ix(ic)))
                if(abs(dx).gt.1.e-20_r8)then
                    t(ic)=(x(ic)-xtab(ix(ic)))/dx
                else
                    t(ic)=0._r8
                endif
            end do
        else
            ix(:ncol)=1
            t(:ncol)=0._r8
        endif
        if(jm.gt.1)then
            do ic=1,ncol
                do j=1,jm
                    if(y(ic).lt.ytab(j))go to 20
                enddo
                20 jy(ic)=max0(j-1,1)
                jp1=min(jy(ic)+1,jm)
                dy=(ytab(jp1)-ytab(jy(ic)))
                if(abs(dy).gt.1.e-20_r8)then
                    u(ic)=(y(ic)-ytab(jy(ic)))/dy
                else
                    u(ic)=0._r8
                endif
            end do
        else
            jy(:ncol)=1
            u(:ncol)=0._r8
        endif
        30 continue
        do ic=1,ncol
            tu(ic)=t(ic)*u(ic)
            tuc(ic)=t(ic)-tu(ic)
            tcuc(ic)=1._r8-tuc(ic)-u(ic)
            tcu(ic)=u(ic)-tu(ic)
            jp1m(ic)=min(jy(ic)+1,jm)
            ip1m(ic)=min(ix(ic)+1,im)
        enddo
        do ic=1,ncol
            jyc=jy(ic)
            ixc=ix(ic)
            jp1=jp1m(ic)
            ip1=ip1m(ic)
            do k=1,km
                out(ic,k) = tcuc(ic) * table(k,ixc,jyc) + tuc(ic) * table(k,ip1,jyc) + &
                              tu(ic) * table(k,ip1,jp1) + tcu(ic) * table(k,ixc,jp1)
            end do
        end do
        return
    end subroutine binterp

end module modal_aer_opt
