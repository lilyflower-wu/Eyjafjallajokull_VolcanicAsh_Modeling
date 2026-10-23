module modal_aero_deposition

!------------------------------------------------------------------------------------------------
! Purpose:
!
! Partition the contributions from modal components of wet and dry 
! deposition at the surface into the fields passed to the coupler.
!
! *** N.B. *** Currently only a simple scheme for the 3-mode version
!              of MAM has been implemented.
!
! Revision history:
! Feb 2009  M. Flanner, B. Eaton   Original version for trop_mam3.
! Jul 2011  F Vitt -- made avaliable to be used in a prescribed modal aerosol mode (no prognostic MAM)
! Mar 2012  F Vitt -- made changes for to prevent abort when 7-mode aeroslol model is used
!                     some of the needed consituents do not exist in 7-mode so bin_fluxes will be false
! May 2014  F Vitt -- included contributions from MAM4 aerosols and added soa_a2 to the ocphiwet fluxes
!------------------------------------------------------------------------------------------------

use shr_kind_mod,     only: r8 => shr_kind_r8
use camsrfexch,       only: cam_out_t     
use constituents,     only: cnst_get_ind, pcnst
use cam_abortutils,   only: endrun
use rad_constituents, only: rad_cnst_get_info
use cam_logfile,      only: iulog

implicit none
private
save

public :: &
   modal_aero_deposition_init, &
   set_srf_drydep,             &
   set_srf_wetdep

! Private module data

logical :: initialized = .false.
integer :: bcphi_ndx( pcnst ) = -1
integer :: bcpho_ndx( pcnst ) = -1
integer :: ocphi_ndx( pcnst ) = -1
integer :: ocpho_ndx( pcnst ) = -1
integer :: crse_dust_ndx( pcnst ) = -1
integer :: fine_dust_ndx( pcnst ) = -1
integer :: ufine_dust_ndx( pcnst ) = -1
integer :: bcphi_cnt = 0
integer :: ocphi_cnt = 0
integer :: bcpho_cnt = 0
integer :: ocpho_cnt = 0
integer :: crse_dust_cnt = 0
integer :: fine_dust_cnt = 0
#if (defined IRONPROC)
integer :: ufine_dust_cnt = 0
integer :: ufine_dufe_cnt = 0
integer :: fine_dufe_cnt = 0
integer :: crse_dufe_cnt = 0
integer :: ufine_anfe_cnt = 0
integer :: fine_anfe_cnt = 0
integer :: crse_anfe_cnt = 0
integer :: ufine_bbfe_cnt = 0
integer :: fine_bbfe_cnt = 0
integer :: crse_bbfe_cnt = 0
integer :: ufine_dusfe_cnt = 0
integer :: fine_dusfe_cnt = 0
integer :: crse_dusfe_cnt = 0
integer :: ufine_ansfe_cnt = 0
integer :: fine_ansfe_cnt = 0
integer :: crse_ansfe_cnt = 0
integer :: ufine_bbsfe_cnt = 0
integer :: fine_bbsfe_cnt = 0
integer :: crse_bbsfe_cnt = 0
#endif
!==============================================================================
contains
!==============================================================================
#if (defined IRONPROC)
subroutine modal_aero_deposition_init( bcphi_indices, bcpho_indices, ocphi_indices, &
                                ocpho_indices, ufine_dust_indices, fine_dust_indices, crse_dust_indices )
#else
subroutine modal_aero_deposition_init( bcphi_indices, bcpho_indices, ocphi_indices, &
                                ocpho_indices, fine_dust_indices, crse_dust_indices )
#endif

  ! set aerosol indices for re-mapping surface deposition fluxes:
  ! *_a1 = accumulation mode
  ! *_a2 = aitken mode
  ! *_a3 = coarse mode
     use cam_logfile,       only: iulog 
  ! can be initialized with user specified indices
  ! if called from aerodep_flx module (for prescribed modal aerosol fluxes) then these indices are specified
  integer, optional, intent(in) :: bcphi_indices(:)     ! hydrophilic black carbon
  integer, optional, intent(in) :: bcpho_indices(:)     ! hydrophobic black carbon
  integer, optional, intent(in) :: ocphi_indices(:)     ! hydrophilic organic carbon
  integer, optional, intent(in) :: ocpho_indices(:)     ! hydrophobic organic carbon 
  integer, optional, intent(in) :: fine_dust_indices(:) ! fine dust
  integer, optional, intent(in) :: crse_dust_indices(:) ! coarse dust
  ! local vars
  integer :: i, pcnt, scnt, scnt1, scnt2, scnt3, scnt4, scnt5, scnt6
#if (defined IRONPROC)
  integer, optional, intent(in) :: ufine_dust_indices(:) ! ufine (Aitken) dust
  integer :: scnt7, scnt8, scnt9, scnt10, scnt11, scnt12, scnt13, scnt14
  integer :: upcnt, uscnt, uscnt1, uscnt2, uscnt3, uscnt4, uscnt5, uscnt6
  integer :: uscnt7, uscnt8, uscnt9, uscnt10, uscnt11, uscnt12, uscnt13, uscnt14
  character(len=16), parameter :: ufine_dust_modes(1) = (/ 'aitken          '/)
#endif
  character(len=16), parameter :: fine_dust_modes(2) =  (/ 'accum           ', 'fine_dust       '/)
  character(len=16), parameter :: crse_dust_modes(2) =  (/ 'coarse          ', 'coarse_dust     '/)
  character(len=16), parameter :: hydrophilic_carbon_modes(1) = (/'accum           '/)
  character(len=16), parameter :: hydrophobic_carbon_modes(3) = (/'aitken          ',  'coarse          ', 'primary_carbon  '/)

  ! if already initialized abort the run
  if (initialized) then
     call endrun('modal_aero_deposition is already initialized')
  endif

  if (present(bcphi_indices)) then
     bcphi_cnt = size(bcphi_indices)
     bcphi_ndx(1:bcphi_cnt) = bcphi_indices (1:bcphi_cnt)
  else
     call get_indices( type='black-c', modes=hydrophilic_carbon_modes, indices=bcphi_ndx, count=bcphi_cnt )
  endif
  if (present(bcpho_indices)) then
     bcpho_cnt = size(bcpho_indices)
     bcpho_ndx(1:bcpho_cnt) = bcpho_indices (1:bcpho_cnt)
  else
     call get_indices( type='black-c', modes=hydrophobic_carbon_modes, indices=bcpho_ndx, count=bcpho_cnt )
  endif

  if (present(ocphi_indices)) then
     ocphi_cnt = size(ocphi_indices)
     ocphi_ndx(1:ocphi_cnt) = ocphi_indices (1:ocphi_cnt)
  else
     call get_indices( type='s-organic', modes=hydrophilic_carbon_modes, indices=ocphi_ndx, count=pcnt )
     call get_indices( type='p-organic', modes=hydrophilic_carbon_modes, indices=ocphi_ndx(pcnt+1:), count=scnt )
     ocphi_cnt = pcnt+scnt
  endif
  if (present(ocpho_indices)) then
     ocpho_cnt = size(ocpho_indices)
     ocpho_ndx(1:ocpho_cnt) = ocpho_indices (1:ocpho_cnt)
  else
     call get_indices( type='s-organic', modes=hydrophobic_carbon_modes, indices=ocpho_ndx, count=pcnt )
     call get_indices( type='p-organic', modes=hydrophobic_carbon_modes, indices=ocpho_ndx(pcnt+1:), count=scnt )
     ocpho_cnt = pcnt+scnt
  endif

  if (present(fine_dust_indices)) then
     fine_dust_cnt = size(fine_dust_indices)
     fine_dust_ndx(1:fine_dust_cnt) = fine_dust_indices(1:fine_dust_cnt)
  else
!-lli
#if ( defined DUSTTYPES )
     call get_indices( type='dust-type1', modes=ufine_dust_modes, indices=ufine_dust_ndx, count=upcnt )
     call get_indices( type='dust-type2', modes=ufine_dust_modes, indices=ufine_dust_ndx(pcnt+1:), count=uscnt )
     call get_indices( type='dust-type3', modes=ufine_dust_modes, indices=ufine_dust_ndx(pcnt+2:), count=uscnt1 )
     call get_indices( type='dust-type4', modes=ufine_dust_modes, indices=ufine_dust_ndx(pcnt+3:), count=uscnt2 )
     call get_indices( type='dust-type5', modes=ufine_dust_modes, indices=ufine_dust_ndx(pcnt+4:), count=uscnt3 )
     call get_indices( type='dust-type6', modes=ufine_dust_modes, indices=ufine_dust_ndx(pcnt+5:), count=uscnt4 )
     call get_indices( type='dust-type7', modes=ufine_dust_modes, indices=ufine_dust_ndx(pcnt+6:), count=uscnt5 )
     call get_indices( type='dust-type8', modes=ufine_dust_modes, indices=ufine_dust_ndx(pcnt+7:), count=uscnt6 )
     call get_indices( type='dust-type1', modes=fine_dust_modes, indices=fine_dust_ndx, count=pcnt )
     call get_indices( type='dust-type2', modes=fine_dust_modes, indices=fine_dust_ndx(pcnt+1:), count=scnt )
     call get_indices( type='dust-type3', modes=fine_dust_modes, indices=fine_dust_ndx(pcnt+2:), count=scnt1 )
     call get_indices( type='dust-type4', modes=fine_dust_modes, indices=fine_dust_ndx(pcnt+3:), count=scnt2 )
     call get_indices( type='dust-type5', modes=fine_dust_modes, indices=fine_dust_ndx(pcnt+4:), count=scnt3 )
     call get_indices( type='dust-type6', modes=fine_dust_modes, indices=fine_dust_ndx(pcnt+5:), count=scnt4 )
     call get_indices( type='dust-type7', modes=fine_dust_modes, indices=fine_dust_ndx(pcnt+6:), count=scnt5 )
     call get_indices( type='dust-type8', modes=fine_dust_modes, indices=fine_dust_ndx(pcnt+7:), count=scnt6 )
#if (defined IRONPROC)
     call get_indices( type='fets      ', modes=ufine_dust_modes, indices=ufine_dust_ndx(pcnt+8:), count=uscnt7 )
     call get_indices( type='fetm      ', modes=ufine_dust_modes, indices=ufine_dust_ndx(pcnt+9:), count=uscnt8 )
     call get_indices( type='fess      ', modes=ufine_dust_modes, indices=ufine_dust_ndx(pcnt+10:),count=uscnt9 )
     call get_indices( type='fesm      ', modes=ufine_dust_modes, indices=ufine_dust_ndx(pcnt+11:),count=uscnt10 )
     call get_indices( type='fetc      ', modes=ufine_dust_modes, indices=ufine_dust_ndx(pcnt+12:),count=uscnt11 )
     call get_indices( type='fesc      ', modes=ufine_dust_modes, indices=ufine_dust_ndx(pcnt+13:),count=uscnt12 )
     call get_indices( type='fetf      ', modes=ufine_dust_modes, indices=ufine_dust_ndx(pcnt+14:),count=uscnt13 )
     call get_indices( type='fesf      ', modes=ufine_dust_modes, indices=ufine_dust_ndx(pcnt+15:),count=uscnt14 )
     call get_indices( type='fets      ', modes=fine_dust_modes, indices=fine_dust_ndx(pcnt+8:), count=scnt7 )
     call get_indices( type='fetm      ', modes=fine_dust_modes, indices=fine_dust_ndx(pcnt+9:), count=scnt8 )
     call get_indices( type='fess      ', modes=fine_dust_modes, indices=fine_dust_ndx(pcnt+10:),count=scnt9 )
     call get_indices( type='fesm      ', modes=fine_dust_modes, indices=fine_dust_ndx(pcnt+11:),count=scnt10 )
     call get_indices( type='fetc      ', modes=fine_dust_modes, indices=fine_dust_ndx(pcnt+12:),count=scnt11 )
     call get_indices( type='fesc      ', modes=fine_dust_modes, indices=fine_dust_ndx(pcnt+13:),count=scnt12 )
     call get_indices( type='fetf      ', modes=fine_dust_modes, indices=fine_dust_ndx(pcnt+14:),count=scnt13 )
     call get_indices( type='fesf      ', modes=fine_dust_modes, indices=fine_dust_ndx(pcnt+15:),count=scnt14 )
     fine_dust_cnt  = pcnt+scnt+scnt1+scnt2+scnt3+scnt4+scnt5+scnt6+scnt7+scnt8+scnt9+scnt10
     ufine_dust_cnt = upcnt+uscnt+uscnt1+uscnt2+uscnt3+uscnt4+uscnt5+uscnt6+uscnt7+uscnt8+uscnt9+uscnt10
     fine_dufe_cnt  = scnt7+scnt8+scnt9+scnt10
     ufine_dufe_cnt = uscnt7+uscnt8+uscnt9+uscnt10
     fine_anfe_cnt  = scnt11+scnt12
     ufine_anfe_cnt = uscnt11+uscnt12
     fine_bbfe_cnt  = scnt13+scnt14
     ufine_bbfe_cnt = uscnt13+uscnt14
     fine_dusfe_cnt = scnt9+scnt10
     ufine_dusfe_cnt= uscnt9+uscnt10
     fine_ansfe_cnt = scnt12
     ufine_ansfe_cnt= uscnt12
     fine_bbsfe_cnt = scnt14
     ufine_bbsfe_cnt= uscnt14
     !write(iulog,*) fine_dust_cnt, 'dsh fine_dust_cnt'
     !write(iulog,*) ufine_dust_cnt, 'dsh ufine_udust_cnt'
     !write(iulog,*) fine_dufe_cnt, 'dsh fine_dufe_cnt,'
     !write(iulog,*) ufine_dufe_cnt, 'dsh fine_dufe_cnt'
     !write(iulog,*) fine_anfe_cnt, 'dshfine_anfe_cnt'
     !write(iulog,*) ufine_anfe_cnt, 'dsh ufine_anfe_cnt'
     !write(iulog,*) fine_bbfe_cnt, 'dsh fine_bbfe_cnt'
     !write(iulog,*) ufine_bbfe_cnt, 'dsh ufine_bbfe_cnt'
     !write(iulog,*) fine_dusfe_cnt, 'dsh fine_dusfe_cnt'
     !write(iulog,*) ufine_dusfe_cnt, 'dsh ufine_dusfe_cnt'
     !write(iulog,*) fine_ansfe_cnt, 'dsh fine_ansfe_cnt'
     !write(iulog,*) ufine_ansfe_cnt, 'dsh ufine_ansfe_cnt'
     !write(iulog,*) fine_bbsfe_cnt, 'dsh fine_bbsfe_cnt'
     !write(iulog,*) ufine_bbsfe_cnt, 'dsh ufine_bbsfe_cnt'
#else
     fine_dust_cnt = pcnt+scnt+scnt1+scnt2+scnt3+scnt4+scnt5+scnt6
#endif
!write(iulog,'(a,i12)') 'fine_dust_cnt=', fine_dust_cnt
#else
     call get_indices( type='dust', modes=fine_dust_modes, indices=fine_dust_ndx, count=fine_dust_cnt )
#endif
  endif
  if (present(crse_dust_indices)) then
     crse_dust_cnt = size(crse_dust_indices)
     crse_dust_ndx(1:crse_dust_cnt) = crse_dust_indices(1:crse_dust_cnt)
  else
#if ( defined DUSTTYPES )
     call get_indices( type='dust-type1', modes=crse_dust_modes, indices=crse_dust_ndx, count=pcnt )
     call get_indices( type='dust-type2', modes=crse_dust_modes, indices=crse_dust_ndx(pcnt+1:), count=scnt )
     call get_indices( type='dust-type3', modes=crse_dust_modes, indices=crse_dust_ndx(pcnt+2:), count=scnt1 )
     call get_indices( type='dust-type4', modes=crse_dust_modes, indices=crse_dust_ndx(pcnt+3:), count=scnt2 )
     call get_indices( type='dust-type5', modes=crse_dust_modes, indices=crse_dust_ndx(pcnt+4:), count=scnt3 )
     call get_indices( type='dust-type6', modes=crse_dust_modes, indices=crse_dust_ndx(pcnt+5:), count=scnt4 )
     call get_indices( type='dust-type7', modes=crse_dust_modes, indices=crse_dust_ndx(pcnt+6:), count=scnt5 )
     call get_indices( type='dust-type8', modes=crse_dust_modes, indices=crse_dust_ndx(pcnt+7:), count=scnt6 )
#if (defined IRONPROC)
     call get_indices( type='fets      ', modes=crse_dust_modes, indices=crse_dust_ndx(pcnt+8:), count=scnt7 )
     call get_indices( type='fetm      ', modes=crse_dust_modes, indices=crse_dust_ndx(pcnt+9:), count=scnt8 )
     call get_indices( type='fess      ', modes=crse_dust_modes, indices=crse_dust_ndx(pcnt+10:),count=scnt9 )
     call get_indices( type='fesm      ', modes=crse_dust_modes, indices=crse_dust_ndx(pcnt+11:),count=scnt10 )
     call get_indices( type='fetc      ', modes=crse_dust_modes, indices=crse_dust_ndx(pcnt+12:),count=scnt11 )
     call get_indices( type='fesc      ', modes=crse_dust_modes, indices=crse_dust_ndx(pcnt+13:),count=scnt12 )
     call get_indices( type='fetf      ', modes=crse_dust_modes, indices=crse_dust_ndx(pcnt+14:),count=scnt13 )
     call get_indices( type='fesf      ', modes=crse_dust_modes, indices=crse_dust_ndx(pcnt+15:),count=scnt14 )
     crse_dust_cnt = pcnt+scnt+scnt1+scnt2+scnt3+scnt4+scnt5+scnt6+scnt7+scnt8+scnt9+scnt10
     crse_dufe_cnt = scnt7+scnt8+scnt9+scnt10
     crse_anfe_cnt = scnt11+scnt12
     crse_bbfe_cnt = scnt13+scnt14
     crse_dusfe_cnt= scnt9+scnt10
     crse_ansfe_cnt= scnt12
     crse_bbsfe_cnt= scnt14
#else
     crse_dust_cnt = pcnt+scnt+scnt1+scnt2+scnt3+scnt4+scnt5+scnt6
#endif
!write(iulog,'(a,i12)') 'crse_dust_cnt=', crse_dust_cnt
#else
     call get_indices( type='dust', modes=crse_dust_modes, indices=crse_dust_ndx, count=crse_dust_cnt )
#endif
  endif

  initialized = .true.

end subroutine modal_aero_deposition_init

!==============================================================================
#if (defined IRONPROC)
subroutine set_srf_wetdep(lchnk, aerdepwetis, aerdepwetcw, cam_out)
  use ppgrid,         only: pcols
  use constituents,   only: cnst_name
  use modal_aero_data,only: cnst_name_cw
  use cam_history,    only: outfld
  integer,  intent(in) :: lchnk ! chunk id
#else
subroutine set_srf_wetdep(aerdepwetis, aerdepwetcw, cam_out)
#endif

! Set surface wet deposition fluxes passed to coupler.

   ! Arguments:
   real(r8), intent(in) :: aerdepwetis(:,:)  ! aerosol wet deposition (interstitial)
   real(r8), intent(in) :: aerdepwetcw(:,:)  ! aerosol wet deposition (cloud water)
   type(cam_out_t), intent(inout) :: cam_out     ! cam export state

   ! Local variables:
   integer :: i, ispec, idx
   integer :: ncol                      ! number of columns

   real(r8) :: bcphiwet_sum, ocphiwet_sum
#if ( defined IRONPROC )
! total arrays
  real(r8) :: dutotdep(pcols)
  real(r8) :: fetotdep(pcols)
  real(r8) :: fesoldep(pcols)
  real(r8) :: fesbydep(pcols)
  real(r8) :: fedutotdep(pcols)
  real(r8) :: fedusoldep(pcols)
  real(r8) :: fedusbydep(pcols)
  real(r8) :: feantotdep(pcols)
  real(r8) :: feansoldep(pcols)
  real(r8) :: feansbydep(pcols)
  real(r8) :: febbtotdep(pcols)
  real(r8) :: febbsoldep(pcols)
  real(r8) :: febbsbydep(pcols)
#endif

   !----------------------------------------------------------------------------

  if (.not.initialized) call endrun('set_srf_wetdep: modal_aero_deposition has not been initialized')

   ncol = cam_out%ncol

   cam_out%bcphiwet(:) = 0._r8
   cam_out%ocphiwet(:) = 0._r8

#if ( defined IRONPROC )
   dutotdep(:) = 0._r8
   fetotdep(:) = 0._r8
   fesoldep(:) = 0._r8
   fesbydep(:) = 0._r8
   fedutotdep(:) = 0._r8
   fedusoldep(:) = 0._r8
   fedusbydep(:) = 0._r8
   feantotdep(:) = 0._r8
   feansoldep(:) = 0._r8
   feansbydep(:) = 0._r8
   febbtotdep(:) = 0._r8
   febbsoldep(:) = 0._r8
   febbsbydep(:) = 0._r8
#endif

   ! derive cam_out variables from deposition fluxes
   !  note: wet deposition fluxes are negative into surface, 
   !        dry deposition fluxes are positive into surface.
   !        srf models want positive definite fluxes.
   do i = 1, ncol

      ! black carbon fluxes
      do ispec=1,bcphi_cnt
         cam_out%bcphiwet(i) = cam_out%bcphiwet(i) &
                             - (aerdepwetis(i,bcphi_ndx(ispec))+aerdepwetcw(i,bcphi_ndx(ispec)))
      enddo
      do ispec=1,bcpho_cnt
         cam_out%bcphiwet(i) = cam_out%bcphiwet(i) &
                             - (aerdepwetis(i,bcpho_ndx(ispec))+aerdepwetcw(i,bcpho_ndx(ispec)))
      enddo

      ! organic carbon fluxes
      do ispec=1,ocphi_cnt
         cam_out%ocphiwet(i) = cam_out%ocphiwet(i) &
                             - (aerdepwetis(i,ocphi_ndx(ispec))+aerdepwetcw(i,ocphi_ndx(ispec)))
      enddo
      do ispec=1,ocpho_cnt
         cam_out%ocphiwet(i) = cam_out%ocphiwet(i) &
                             - (aerdepwetis(i,ocpho_ndx(ispec))+aerdepwetcw(i,ocpho_ndx(ispec)))
      enddo

      ! dust fluxes
      cam_out%dstwet1(i) = 0._r8
      cam_out%dstwet2(i) = 0._r8
      cam_out%dstwet3(i) = 0._r8
      cam_out%dstwet4(i) = 0._r8

#if ( defined IRONPROC )
      ! bulk bin1 (fine) dust deposition equals Aitken and accumulation mode deposition:
      do ispec=1,ufine_dust_cnt
         cam_out%dstwet1(i) = cam_out%dstwet1(i) &
                            -(aerdepwetis(i,ufine_dust_ndx(ispec))+aerdepwetcw(i,ufine_dust_ndx(ispec)))
         !write(iulog,*) ispec,'ufine_dust_cnt', ufine_dust_cnt, cnst_name(ufine_dust_ndx(ispec)), cnst_name_cw(ufine_dust_ndx(ispec))
      enddo
      do ispec=1,fine_dust_cnt
         cam_out%dstwet1(i) = cam_out%dstwet1(i) &
                            -(aerdepwetis(i,fine_dust_ndx(ispec))+aerdepwetcw(i,fine_dust_ndx(ispec)))
         !write(iulog,*) ispec, 'fine_dust_cnt', fine_dust_cnt, cnst_name(fine_dust_ndx(ispec)), cnst_name_cw(fine_dust_ndx(ispec))
      enddo

      !  Assign all coarse-mode dust to bulk size bin 3:
      do ispec=1,crse_dust_cnt
         cam_out%dstwet3(i) = cam_out%dstwet3(i) &
                            -(aerdepwetis(i,crse_dust_ndx(ispec))+aerdepwetcw(i,crse_dust_ndx(ispec)))
         !write(iulog,*) ispec, 'crse_dust_cnt', crse_dust_cnt, cnst_name(crse_dust_ndx(ispec)), cnst_name_cw(crse_dust_ndx(ispec))
      enddo

      !  total iron dust
      do ispec=9,9+ufine_dufe_cnt-1
         fedutotdep(i) = fedutotdep(i) &
                       -(aerdepwetis(i,ufine_dust_ndx(ispec))+aerdepwetcw(i,ufine_dust_ndx(ispec)))
         !write(iulog,*) ispec, 'ufine_dufe_cnt', ufine_dufe_cnt, cnst_name(ufine_dust_ndx(ispec)), cnst_name_cw(ufine_dust_ndx(ispec))
      enddo
      do ispec=9,9+fine_dufe_cnt-1
         fedutotdep(i) = fedutotdep(i) &
                       -(aerdepwetis(i,fine_dust_ndx(ispec))+aerdepwetcw(i,fine_dust_ndx(ispec)))
         !write(iulog,*) ispec, 'fine_dufe_cnt', fine_dufe_cnt, cnst_name(fine_dust_ndx(ispec)), cnst_name_cw(fine_dust_ndx(ispec))
      enddo
      do ispec=9,9+crse_dufe_cnt-1
         fedutotdep(i) = fedutotdep(i) &
                       -(aerdepwetis(i,crse_dust_ndx(ispec))+aerdepwetcw(i,crse_dust_ndx(ispec)))
         !write(iulog,*) ispec, 'crse_dufe_cnt', crse_dufe_cnt, cnst_name(crse_dust_ndx(ispec)), cnst_name_cw(crse_dust_ndx(ispec))
      enddo

      !  solub iron dust
      do ispec=11,11+ufine_dusfe_cnt-1
         fedusoldep(i) = fedusoldep(i) &
                       -(aerdepwetis(i,ufine_dust_ndx(ispec))+aerdepwetcw(i,ufine_dust_ndx(ispec)))
         !write(iulog,*) ispec, 'ufine_dusfe_cnt', ufine_dusfe_cnt, cnst_name(ufine_dust_ndx(ispec)), cnst_name_cw(ufine_dust_ndx(ispec))
      enddo
      do ispec=11,11+fine_dusfe_cnt-1
         fedusoldep(i) = fedusoldep(i) &
                       -(aerdepwetis(i,fine_dust_ndx(ispec))+aerdepwetcw(i,fine_dust_ndx(ispec)))
         !write(iulog,*) ispec, 'fine_dusfe_cnt', fine_dusfe_cnt, cnst_name(fine_dust_ndx(ispec)), cnst_name_cw(fine_dust_ndx(ispec))
      enddo
      do ispec=11,11+crse_dusfe_cnt-1
         fedusoldep(i) = fedusoldep(i) &
                       -(aerdepwetis(i,crse_dust_ndx(ispec))+aerdepwetcw(i,crse_dust_ndx(ispec)))
         !write(iulog,*) ispec, 'crse_dusfe_cnt', crse_dusfe_cnt, cnst_name(crse_dust_ndx(ispec)), cnst_name_cw(crse_dust_ndx(ispec))
      enddo

      ! total iron anth
      do ispec=13,13+ufine_anfe_cnt-1
         feantotdep(i) = feantotdep(i) &
                       -(aerdepwetis(i,ufine_dust_ndx(ispec))+aerdepwetcw(i,ufine_dust_ndx(ispec)))
         !write(iulog,*) ispec, 'ufine_anfe_cnt', ufine_anfe_cnt, cnst_name(ufine_dust_ndx(ispec)), cnst_name_cw(ufine_dust_ndx(ispec))
      enddo
      do ispec=13,13+fine_anfe_cnt-1
         feantotdep(i) = feantotdep(i) &
                       -(aerdepwetis(i,fine_dust_ndx(ispec))+aerdepwetcw(i,fine_dust_ndx(ispec)))
         !write(iulog,*) ispec, 'fine_anfe_cnt', fine_anfe_cnt, cnst_name(fine_dust_ndx(ispec)), cnst_name_cw(fine_dust_ndx(ispec))
      enddo
      do ispec=13,13+crse_anfe_cnt-1
         feantotdep(i) = feantotdep(i) &
                       -(aerdepwetis(i,crse_dust_ndx(ispec))+aerdepwetcw(i,crse_dust_ndx(ispec)))
         !write(iulog,*) ispec, 'crse_anfe_cnt', crse_anfe_cnt, cnst_name(crse_dust_ndx(ispec)), cnst_name_cw(crse_dust_ndx(ispec))
      enddo

      ! solub iron anth
      do ispec=14,14+ufine_ansfe_cnt-1
         feansoldep(i) = feansoldep(i) &
                       -(aerdepwetis(i,ufine_dust_ndx(ispec))+aerdepwetcw(i,ufine_dust_ndx(ispec)))
         !write(iulog,*) ispec, 'ufine_ansfe_cnt', ufine_ansfe_cnt, cnst_name(ufine_dust_ndx(ispec)), cnst_name_cw(ufine_dust_ndx(ispec))
      enddo
      do ispec=14,14+fine_ansfe_cnt-1
         feansoldep(i) = feansoldep(i) &
                       -(aerdepwetis(i,fine_dust_ndx(ispec))+aerdepwetcw(i,fine_dust_ndx(ispec)))
         !write(iulog,*) ispec, 'fine_ansfe_cnt', fine_ansfe_cnt, cnst_name(fine_dust_ndx(ispec)), cnst_name_cw(fine_dust_ndx(ispec))
      enddo
      do ispec=14,14+crse_ansfe_cnt-1
         feansoldep(i) = feansoldep(i) &
                       -(aerdepwetis(i,crse_dust_ndx(ispec))+aerdepwetcw(i,crse_dust_ndx(ispec)))
         !write(iulog,*) ispec, 'crse_ansfe_cnt', crse_ansfe_cnt, cnst_name(crse_dust_ndx(ispec)), cnst_name_cw(crse_dust_ndx(ispec))
      enddo

      ! total iron fire 
      do ispec=15,15+ufine_bbfe_cnt-1
         febbtotdep(i) = febbtotdep(i) &
                       -(aerdepwetis(i,ufine_dust_ndx(ispec))+aerdepwetcw(i,ufine_dust_ndx(ispec)))
         !write(iulog,*) ispec, 'ufine_bbfe_cnt', ufine_bbfe_cnt, cnst_name(ufine_dust_ndx(ispec)), cnst_name_cw(ufine_dust_ndx(ispec))
      enddo
      do ispec=15,15+fine_bbfe_cnt-1
         febbtotdep(i) = febbtotdep(i) &
                       -(aerdepwetis(i,fine_dust_ndx(ispec))+aerdepwetcw(i,fine_dust_ndx(ispec)))
         !write(iulog,*) ispec, 'fine_bbfe_cnt', fine_bbfe_cnt, cnst_name(fine_dust_ndx(ispec)), cnst_name_cw(fine_dust_ndx(ispec))
      enddo
      do ispec=15,15+crse_bbfe_cnt-1
         febbtotdep(i) = febbtotdep(i) &
                       -(aerdepwetis(i,crse_dust_ndx(ispec))+aerdepwetcw(i,crse_dust_ndx(ispec)))
         !write(iulog,*) ispec, 'crse_bbfe_cnt', crse_bbfe_cnt, cnst_name(crse_dust_ndx(ispec)), cnst_name_cw(crse_dust_ndx(ispec))
      enddo

      ! solub iron fire 
      do ispec=16,16+ufine_bbsfe_cnt-1
         febbsoldep(i) = febbsoldep(i) &
                       -(aerdepwetis(i,ufine_dust_ndx(ispec))+aerdepwetcw(i,ufine_dust_ndx(ispec)))
         !write(iulog,*) ispec, 'ufine_bbsfe_cnt', ufine_bbsfe_cnt, cnst_name(ufine_dust_ndx(ispec)), cnst_name_cw(ufine_dust_ndx(ispec))
      enddo
      do ispec=16,16+fine_bbsfe_cnt-1
         febbsoldep(i) = febbsoldep(i) &
                       -(aerdepwetis(i,fine_dust_ndx(ispec))+aerdepwetcw(i,fine_dust_ndx(ispec)))
         !write(iulog,*) ispec, 'fine_bbsfe_cnt', fine_bbsfe_cnt, cnst_name(fine_dust_ndx(ispec)), cnst_name_cw(fine_dust_ndx(ispec))
      enddo
      do ispec=16,16+crse_bbsfe_cnt-1
         febbsoldep(i) = febbsoldep(i) &
                       -(aerdepwetis(i,crse_dust_ndx(ispec))+aerdepwetcw(i,crse_dust_ndx(ispec)))
         !write(iulog,*) ispec, 'crse_bbsfe_cnt', crse_bbsfe_cnt, cnst_name(crse_dust_ndx(ispec)), cnst_name_cw(crse_dust_ndx(ispec))
      enddo

      dutotdep(i) = cam_out%dstwet1(i) + cam_out%dstwet3(i)
      fetotdep(i) = fedutotdep(i) + feantotdep(i) + febbtotdep(i)
      fesoldep(i) = fedusoldep(i) + feansoldep(i) + febbsoldep(i)

      if (fesoldep(i) .ge. 1.e-30_r8) then   ! only give value if wet dep is present and occuring, otherwise 0%
       ! bound %solubility from 0.001-100%
       fesbydep(i)   = max(min((fesoldep(i)  /fetotdep(i)  *100._r8),100._r8),0.001_r8)
       fedusbydep(i) = max(min((fedusoldep(i)/fedutotdep(i)*100._r8),100._r8),0.001_r8)
       feansbydep(i) = max(min((feansoldep(i)/feantotdep(i)*100._r8),100._r8),0.001_r8)
       febbsbydep(i) = max(min((febbsoldep(i)/febbtotdep(i)*100._r8),100._r8),0.001_r8)
      else
       fesbydep(i)=0.0_r8
       fedusbydep(i)=0.0_r8
       feansbydep(i)=0.0_r8
       febbsbydep(i)=0.0_r8
      endif
#else
      ! bulk bin1 (fine) dust deposition equals accumulation mode deposition:
      do ispec=1,fine_dust_cnt
         cam_out%dstwet1(i) = cam_out%dstwet1(i) &
                            -(aerdepwetis(i,fine_dust_ndx(ispec))+aerdepwetcw(i,fine_dust_ndx(ispec)))
      enddo

      !  Assign all coarse-mode dust to bulk size bin 3:
      do ispec=1,crse_dust_cnt
         cam_out%dstwet3(i) = cam_out%dstwet3(i) &
                            -(aerdepwetis(i,crse_dust_ndx(ispec))+aerdepwetcw(i,crse_dust_ndx(ispec)))
      enddo
#endif
      ! in rare cases, integrated deposition tendency is upward
      if (cam_out%bcphiwet(i) .lt. 0._r8) cam_out%bcphiwet(i) = 0._r8
      if (cam_out%ocphiwet(i) .lt. 0._r8) cam_out%ocphiwet(i) = 0._r8
      if (cam_out%dstwet1(i)  .lt. 0._r8) cam_out%dstwet1(i)  = 0._r8
      if (cam_out%dstwet3(i)  .lt. 0._r8) cam_out%dstwet3(i)  = 0._r8
   enddo

#if ( defined IRONPROC )
      call outfld( 'DUTOTWET',   dutotdep,   pcols, lchnk)
      call outfld( 'FETOTWET',   fetotdep,   pcols, lchnk)
      call outfld( 'FESOLWET',   fesoldep,   pcols, lchnk)
      call outfld( 'FESBYWET',   fesbydep,   pcols, lchnk)
      call outfld( 'FEDUTOTWET', fedutotdep, pcols, lchnk)
      call outfld( 'FEDUSOLWET', fedusoldep, pcols, lchnk)
      call outfld( 'FEDUSBYWET', fedusbydep, pcols, lchnk)
      call outfld( 'FEANTOTWET', feantotdep, pcols, lchnk)
      call outfld( 'FEANSOLWET', feansoldep, pcols, lchnk)
      call outfld( 'FEANSBYWET', feansbydep, pcols, lchnk)
      call outfld( 'FEBBTOTWET', febbtotdep, pcols, lchnk)
      call outfld( 'FEBBSOLWET', febbsoldep, pcols, lchnk)
      call outfld( 'FEBBSBYWET', febbsbydep, pcols, lchnk)
#endif

end subroutine set_srf_wetdep

!==============================================================================

#if (defined IRONPROC)
subroutine set_srf_drydep(lchnk, aerdepdryis, aerdepdrycw, cam_out)
  use cam_history,    only: outfld
  use constituents,   only: cnst_name
  use modal_aero_data,only: cnst_name_cw
  use ppgrid,         only: pcols
  integer,  intent(in) :: lchnk ! chunk id
#else
subroutine set_srf_drydep(aerdepdryis, aerdepdrycw, cam_out)
#endif
! Set surface dry deposition fluxes passed to coupler.
   
   ! Arguments:
   real(r8), intent(in) :: aerdepdryis(:,:)  ! aerosol dry deposition (interstitial)
   real(r8), intent(in) :: aerdepdrycw(:,:)  ! aerosol dry deposition (cloud water)
   type(cam_out_t), intent(inout) :: cam_out     ! cam export state

   ! Local variables:
   integer :: i, ispec, idx
   integer :: ncol                      ! number of columns
   real(r8):: bcphidry_sum, ocphidry_sum, ocphodry_sum
#if ( defined IRONPROC )
! total arrays
  real(r8) :: dutotdep(pcols)
  real(r8) :: fetotdep(pcols)
  real(r8) :: fesoldep(pcols)
  real(r8) :: fesbydep(pcols)
  real(r8) :: fedutotdep(pcols)
  real(r8) :: fedusoldep(pcols)
  real(r8) :: fedusbydep(pcols)
  real(r8) :: feantotdep(pcols)
  real(r8) :: feansoldep(pcols)
  real(r8) :: feansbydep(pcols)
  real(r8) :: febbtotdep(pcols)
  real(r8) :: febbsoldep(pcols)
  real(r8) :: febbsbydep(pcols)
#endif
   !----------------------------------------------------------------------------

   if (.not.initialized) call endrun('set_srf_drydep: modal_aero_deposition has not been initialized')

   ncol = cam_out%ncol

   cam_out%bcphidry(:) = 0._r8
   cam_out%bcphodry(:) = 0._r8
   cam_out%ocphidry(:) = 0._r8
   cam_out%ocphodry(:) = 0._r8
#if ( defined IRONPROC )
   dutotdep(:) = 0._r8
   fetotdep(:) = 0._r8
   fesoldep(:) = 0._r8
   fesbydep(:) = 0._r8
   fedutotdep(:) = 0._r8
   fedusoldep(:) = 0._r8
   fedusbydep(:) = 0._r8
   feantotdep(:) = 0._r8
   feansoldep(:) = 0._r8
   feansbydep(:) = 0._r8
   febbtotdep(:) = 0._r8
   febbsoldep(:) = 0._r8
   febbsbydep(:) = 0._r8
#endif

   ! derive cam_out variables from deposition fluxes
   !  note: wet deposition fluxes are negative into surface, 
   !        dry deposition fluxes are positive into surface.
   !        srf models want positive definite fluxes.
   do i = 1, ncol

      ! black carbon fluxes
      do ispec=1,bcphi_cnt
         cam_out%bcphidry(i) = cam_out%bcphidry(i) &
                             + (aerdepdryis(i,bcphi_ndx(ispec))+aerdepdrycw(i,bcphi_ndx(ispec)))
      enddo
      do ispec=1,bcpho_cnt
         cam_out%bcphodry(i) = cam_out%bcphodry(i) &
                             + (aerdepdryis(i,bcpho_ndx(ispec))+aerdepdrycw(i,bcpho_ndx(ispec)))
      enddo

      ! organic carbon fluxes
      do ispec=1,ocphi_cnt
         cam_out%ocphidry(i) = cam_out%ocphidry(i) &
                             + (aerdepdryis(i,ocphi_ndx(ispec))+aerdepdrycw(i,ocphi_ndx(ispec)))
      enddo
      do ispec=1,ocpho_cnt
         cam_out%ocphodry(i) = cam_out%ocphodry(i) &
                             + (aerdepdryis(i,ocpho_ndx(ispec))+aerdepdrycw(i,ocpho_ndx(ispec)))
      enddo

      ! dust fluxes
      cam_out%dstdry1(i) = 0._r8
      cam_out%dstdry2(i) = 0._r8
      cam_out%dstdry3(i) = 0._r8
      cam_out%dstdry4(i) = 0._r8
#if ( defined IRONPROC )
      ! bulk bin1 (fine) dust deposition equals accumulation mode deposition:
      do ispec=1,ufine_dust_cnt
         cam_out%dstdry1(i) = cam_out%dstdry1(i) &
                            +(aerdepdryis(i,ufine_dust_ndx(ispec))+aerdepdrycw(i,ufine_dust_ndx(ispec)))
         !write(iulog,*) ispec,'ufine_dust_cnt', ufine_dust_cnt, cnst_name(ufine_dust_ndx(ispec)), cnst_name_cw(ufine_dust_ndx(ispec))
      enddo
      do ispec=1,fine_dust_cnt
         cam_out%dstdry1(i) = cam_out%dstdry1(i) &
                            +(aerdepdryis(i,fine_dust_ndx(ispec))+aerdepdrycw(i,fine_dust_ndx(ispec)))
         !write(iulog,*) ispec, 'fine_dust_cnt', fine_dust_cnt, cnst_name(fine_dust_ndx(ispec)), cnst_name_cw(fine_dust_ndx(ispec))
      enddo

      !  Assign all coarse-mode dust to bulk size bin 3:
      do ispec=1,crse_dust_cnt
         cam_out%dstdry3(i) = cam_out%dstdry3(i) &
                            +(aerdepdryis(i,crse_dust_ndx(ispec))+aerdepdrycw(i,crse_dust_ndx(ispec)))
         !write(iulog,*) ispec, 'crse_dust_cnt', crse_dust_cnt, cnst_name(crse_dust_ndx(ispec)), cnst_name_cw(crse_dust_ndx(ispec))
      enddo

      !  total iron dust
      do ispec=9,9+ufine_dufe_cnt-1
         fedutotdep(i) = fedutotdep(i) &
                       +(aerdepdryis(i,ufine_dust_ndx(ispec))+aerdepdrycw(i,ufine_dust_ndx(ispec)))
         !write(iulog,*) ispec, 'ufine_dufe_cnt', ufine_dufe_cnt, cnst_name(ufine_dust_ndx(ispec)), cnst_name_cw(ufine_dust_ndx(ispec))
      enddo
      do ispec=9,9+fine_dufe_cnt-1
         fedutotdep(i) = fedutotdep(i) &
                       +(aerdepdryis(i,fine_dust_ndx(ispec))+aerdepdrycw(i,fine_dust_ndx(ispec)))
         !write(iulog,*) ispec, 'fine_dufe_cnt', fine_dufe_cnt, cnst_name(fine_dust_ndx(ispec)), cnst_name_cw(fine_dust_ndx(ispec))
      enddo
      do ispec=9,9+crse_dufe_cnt-1
         fedutotdep(i) = fedutotdep(i) &
                       +(aerdepdryis(i,crse_dust_ndx(ispec))+aerdepdrycw(i,crse_dust_ndx(ispec)))
         !write(iulog,*) ispec, 'crse_dufe_cnt', crse_dufe_cnt, cnst_name(crse_dust_ndx(ispec)), cnst_name_cw(crse_dust_ndx(ispec))
      enddo

      !  solub iron dust
      do ispec=11,11+ufine_dusfe_cnt-1
         fedusoldep(i) = fedusoldep(i) &
                       +(aerdepdryis(i,ufine_dust_ndx(ispec))+aerdepdrycw(i,ufine_dust_ndx(ispec)))
         !write(iulog,*) ispec, 'ufine_dusfe_cnt', ufine_dusfe_cnt, cnst_name(ufine_dust_ndx(ispec)), cnst_name_cw(ufine_dust_ndx(ispec))
      enddo
      do ispec=11,11+fine_dusfe_cnt-1
         fedusoldep(i) = fedusoldep(i) &
                       +(aerdepdryis(i,fine_dust_ndx(ispec))+aerdepdrycw(i,fine_dust_ndx(ispec)))
         !write(iulog,*) ispec, 'fine_dusfe_cnt', fine_dusfe_cnt, cnst_name(fine_dust_ndx(ispec)), cnst_name_cw(fine_dust_ndx(ispec))
      enddo
      do ispec=11,11+crse_dusfe_cnt-1
         fedusoldep(i) = fedusoldep(i) &
                       +(aerdepdryis(i,crse_dust_ndx(ispec))+aerdepdrycw(i,crse_dust_ndx(ispec)))
         !write(iulog,*) ispec, 'crse_dusfe_cnt', crse_dusfe_cnt, cnst_name(crse_dust_ndx(ispec)), cnst_name_cw(crse_dust_ndx(ispec))
      enddo

      ! total iron anth
      do ispec=13,13+ufine_anfe_cnt-1
         feantotdep(i) = feantotdep(i) &
                       +(aerdepdryis(i,ufine_dust_ndx(ispec))+aerdepdrycw(i,ufine_dust_ndx(ispec)))
         !write(iulog,*) ispec, 'ufine_anfe_cnt', ufine_anfe_cnt, cnst_name(ufine_dust_ndx(ispec)), cnst_name_cw(ufine_dust_ndx(ispec))
      enddo
      do ispec=13,13+fine_anfe_cnt-1
         feantotdep(i) = feantotdep(i) &
                       +(aerdepdryis(i,fine_dust_ndx(ispec))+aerdepdrycw(i,fine_dust_ndx(ispec)))
         !write(iulog,*) ispec, 'fine_anfe_cnt', fine_anfe_cnt, cnst_name(fine_dust_ndx(ispec)), cnst_name_cw(fine_dust_ndx(ispec))
      enddo
      do ispec=13,13+crse_anfe_cnt-1
         feantotdep(i) = feantotdep(i) &
                       +(aerdepdryis(i,crse_dust_ndx(ispec))+aerdepdrycw(i,crse_dust_ndx(ispec)))
         !write(iulog,*) ispec, 'crse_anfe_cnt', crse_anfe_cnt, cnst_name(crse_dust_ndx(ispec)), cnst_name_cw(crse_dust_ndx(ispec))
      enddo

      ! solub iron anth
      do ispec=14,14+ufine_ansfe_cnt-1
         feansoldep(i) = feansoldep(i) &
                       +(aerdepdryis(i,ufine_dust_ndx(ispec))+aerdepdrycw(i,ufine_dust_ndx(ispec)))
         !write(iulog,*) ispec, 'ufine_ansfe_cnt', ufine_ansfe_cnt, cnst_name(ufine_dust_ndx(ispec)), cnst_name_cw(ufine_dust_ndx(ispec))
      enddo
      do ispec=14,14+fine_ansfe_cnt-1
         feansoldep(i) = feansoldep(i) &
                       +(aerdepdryis(i,fine_dust_ndx(ispec))+aerdepdrycw(i,fine_dust_ndx(ispec)))
         !write(iulog,*) ispec, 'fine_ansfe_cnt', fine_ansfe_cnt, cnst_name(fine_dust_ndx(ispec)), cnst_name_cw(fine_dust_ndx(ispec))
      enddo
      do ispec=14,14+crse_ansfe_cnt-1
         feansoldep(i) = feansoldep(i) &
                       +(aerdepdryis(i,crse_dust_ndx(ispec))+aerdepdrycw(i,crse_dust_ndx(ispec)))
         !write(iulog,*) ispec, 'crse_ansfe_cnt', crse_ansfe_cnt, cnst_name(crse_dust_ndx(ispec)), cnst_name_cw(crse_dust_ndx(ispec))
      enddo

      ! total iron fire 
      do ispec=15,15+ufine_bbfe_cnt-1
         febbtotdep(i) = febbtotdep(i) &
                       +(aerdepdryis(i,ufine_dust_ndx(ispec))+aerdepdrycw(i,ufine_dust_ndx(ispec)))
         !write(iulog,*) ispec, 'ufine_bbfe_cnt', ufine_bbfe_cnt, cnst_name(ufine_dust_ndx(ispec)), cnst_name_cw(ufine_dust_ndx(ispec))
      enddo
      do ispec=15,15+fine_bbfe_cnt-1
         febbtotdep(i) = febbtotdep(i) &
                       +(aerdepdryis(i,fine_dust_ndx(ispec))+aerdepdrycw(i,fine_dust_ndx(ispec)))
         !write(iulog,*) ispec, 'fine_bbfe_cnt', fine_bbfe_cnt, cnst_name(fine_dust_ndx(ispec)), cnst_name_cw(fine_dust_ndx(ispec))
      enddo
      do ispec=15,15+crse_bbfe_cnt-1
         febbtotdep(i) = febbtotdep(i) &
                       +(aerdepdryis(i,crse_dust_ndx(ispec))+aerdepdrycw(i,crse_dust_ndx(ispec)))
         !write(iulog,*) ispec, 'crse_bbfe_cnt', crse_bbfe_cnt, cnst_name(crse_dust_ndx(ispec)), cnst_name_cw(crse_dust_ndx(ispec))
      enddo

      ! solub iron fire 
      do ispec=16,16+ufine_bbsfe_cnt-1
         febbsoldep(i) = febbsoldep(i) &
                       +(aerdepdryis(i,ufine_dust_ndx(ispec))+aerdepdrycw(i,ufine_dust_ndx(ispec)))
         !write(iulog,*) ispec, 'ufine_bbsfe_cnt', ufine_bbsfe_cnt, cnst_name(ufine_dust_ndx(ispec)), cnst_name_cw(ufine_dust_ndx(ispec))
      enddo
      do ispec=16,16+fine_bbsfe_cnt-1
         febbsoldep(i) = febbsoldep(i) &
                       +(aerdepdryis(i,fine_dust_ndx(ispec))+aerdepdrycw(i,fine_dust_ndx(ispec)))
         !write(iulog,*) ispec, 'fine_bbsfe_cnt', fine_bbsfe_cnt, cnst_name(fine_dust_ndx(ispec)), cnst_name_cw(fine_dust_ndx(ispec))
      enddo
      do ispec=16,16+crse_bbsfe_cnt-1
         febbsoldep(i) = febbsoldep(i) &
                       +(aerdepdryis(i,crse_dust_ndx(ispec))+aerdepdrycw(i,crse_dust_ndx(ispec)))
         !write(iulog,*) ispec, 'crse_bbsfe_cnt', crse_bbsfe_cnt, cnst_name(crse_dust_ndx(ispec)), cnst_name_cw(crse_dust_ndx(ispec))
      enddo

      dutotdep(i) = cam_out%dstdry1(i) + cam_out%dstdry3(i)
      fetotdep(i) = fedutotdep(i) + feantotdep(i) + febbtotdep(i)
      fesoldep(i) = fedusoldep(i) + feansoldep(i) + febbsoldep(i)

      ! bound %solubility from 0.001-100%
      fesbydep(i)   = max(min((fesoldep(i)  /fetotdep(i)  *100._r8),100._r8),0.001_r8)
      fedusbydep(i) = max(min((fedusoldep(i)/fedutotdep(i)*100._r8),100._r8),0.001_r8)
      feansbydep(i) = max(min((feansoldep(i)/feantotdep(i)*100._r8),100._r8),0.001_r8)
      febbsbydep(i) = max(min((febbsoldep(i)/febbtotdep(i)*100._r8),100._r8),0.001_r8)
#else
      ! bulk bin1 (fine) dust deposition equals accumulation mode deposition:
      do ispec=1,fine_dust_cnt
         cam_out%dstdry1(i) = cam_out%dstdry1(i) &
                            + (aerdepdryis(i,fine_dust_ndx(ispec))+aerdepdrycw(i,fine_dust_ndx(ispec)))
      enddo
      !  Assign all coarse-mode dust to bulk size bin 3:
      do ispec=1,crse_dust_cnt
         cam_out%dstdry3(i) = cam_out%dstdry3(i) &
                            + (aerdepdryis(i,crse_dust_ndx(ispec))+aerdepdrycw(i,crse_dust_ndx(ispec)))
      enddo
#endif
      ! in rare cases, integrated deposition tendency is upward
      if (cam_out%bcphidry(i) .lt. 0._r8) cam_out%bcphidry(i) = 0._r8
      if (cam_out%bcphodry(i) .lt. 0._r8) cam_out%bcphodry(i) = 0._r8
      if (cam_out%ocphidry(i) .lt. 0._r8) cam_out%ocphidry(i) = 0._r8
      if (cam_out%ocphodry(i) .lt. 0._r8) cam_out%ocphodry(i) = 0._r8
      if (cam_out%dstdry1(i)  .lt. 0._r8) cam_out%dstdry1(i)  = 0._r8
      if (cam_out%dstdry3(i)  .lt. 0._r8) cam_out%dstdry3(i)  = 0._r8
   enddo

#if ( defined IRONPROC )
      call outfld( 'DUTOTDRY',   dutotdep,   pcols, lchnk)
      call outfld( 'FETOTDRY',   fetotdep,   pcols, lchnk)
      call outfld( 'FESOLDRY',   fesoldep,   pcols, lchnk)
      call outfld( 'FESBYDRY',   fesbydep,   pcols, lchnk)
      call outfld( 'FEDUTOTDRY', fedutotdep, pcols, lchnk)
      call outfld( 'FEDUSOLDRY', fedusoldep, pcols, lchnk)
      call outfld( 'FEDUSBYDRY', fedusbydep, pcols, lchnk)
      call outfld( 'FEANTOTDRY', feantotdep, pcols, lchnk)
      call outfld( 'FEANSOLDRY', feansoldep, pcols, lchnk)
      call outfld( 'FEANSBYDRY', feansbydep, pcols, lchnk)
      call outfld( 'FEBBTOTDRY', febbtotdep, pcols, lchnk)
      call outfld( 'FEBBSOLDRY', febbsoldep, pcols, lchnk)
      call outfld( 'FEBBSBYDRY', febbsbydep, pcols, lchnk)
#endif

end subroutine set_srf_drydep

!==============================================================================
subroutine get_indices( type, modes, indices, count )
use cam_logfile,     only: iulog
  character(len=*), intent(in) :: type
  character(len=*), intent(in) :: modes(:)
  integer, intent(out) :: indices(:)
  integer, intent(out) :: count

  integer :: l, n, ndx, nmodes, nspec
  character(len=32) :: spec_type, spec_name, mode_type

  call rad_cnst_get_info(0, nmodes=nmodes)

  count = 0
  indices(:) = -1

  if (nmodes==7) return ! historically turned off for mam7

  do n = 1, nmodes

     call rad_cnst_get_info(0, n, mode_type=mode_type, nspec=nspec)

     if ( any(modes==trim(mode_type)) ) then

        do l = 1,nspec

           call rad_cnst_get_info(0, n, l, spec_type=spec_type, spec_name=spec_name)
           call cnst_get_ind(spec_name, ndx, abort=.false.)
           if (ndx>0) then
              if (trim(spec_type) == trim(type)) then
                 count = count+1
                 indices(count) = ndx
              endif
           endif
        enddo

     endif

  enddo

end subroutine get_indices
!==============================================================================

end module modal_aero_deposition
