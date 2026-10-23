!----------------------------------------------------------------------------------
! Modal aerosol implementation
!
! Based on the sox_cldaero_mod but implemeneted for iron chemistry
! Update of CAM4 scheme in Scanza et al 2018 
!
! Douglas Hamilton 2018
!----------------------------------------------------------------------------------
module fe_cldaero_mod

  use shr_kind_mod,    only : r8 => shr_kind_r8
  use cam_abortutils,  only : endrun
  use ppgrid,          only : pcols, pver
  use mo_chem_utls,    only : get_spc_ndx
  use cldaero_mod,     only : cldaero_conc_t, cldaero_allocate, cldaero_deallocate
  use modal_aero_data, only : ntot_amode, lptr2_soa_cw_amode,lptr2_soa_a_amode, lptr2_pom_cw_amode,nsoa
  use modal_aero_data, only : numptrcw_amode, lptr_dust6_a_amode, lptr_so4_a_amode
  use modal_aero_data, only : lptr_fess_a_amode, lptr_fesm_a_amode, lptr_fets_a_amode, lptr_fetm_a_amode
  use modal_aero_data, only : lptr_fess_cw_amode, lptr_fesm_cw_amode, lptr_fets_cw_amode, lptr_fetm_cw_amode
  use modal_aero_data, only : lptr_fesc_a_amode, lptr_fesc_cw_amode, lptr_fetc_a_amode, lptr_fetc_cw_amode
  use modal_aero_data, only : lptr_fesf_a_amode, lptr_fesf_cw_amode, lptr_fetf_a_amode, lptr_fetf_cw_amode
  use modal_aero_data, only : cnst_name_cw, specmw_so4_amode
  use cam_history,     only : outfld, horiz_only
  use cam_history,     only : addfld, add_default
  use chem_mods,       only : adv_mass
  use physconst,       only : gravit
  use phys_control,    only : phys_getopts
  use cldaero_mod,     only : cldaero_uptakerate
  use chem_mods,       only : gas_pcnst
  use cam_logfile,     only : iulog
  use physics_types,   only : physics_state

  implicit none
  private

  public :: fe_cldaero_init
  public :: fe_cldaero_update

  integer :: id_fetsa

  real(r8), parameter :: small_value = 1.e-35_r8

contains

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------

  subroutine fe_cldaero_init

    integer :: l, m
    logical :: history_aerosol      ! Output the MAM aerosol tendencies

    id_fetsa = get_spc_ndx('fets_a1')

    write(iulog,*)'spc_ndx fe', id_fetsa

    if (id_fetsa <1 ) then
      call endrun('fe_cldaero_init:MAM mech does not include necessary iron species' &
                  //' -- should not invoke fe_cldaero_mod ')
    endif
  
    call phys_getopts( history_aerosol_out        = history_aerosol   )
    !
    !   add to history
    !
    call addfld ('fesssrc_c1', (/ 'lev' /), 'A', 'kg/kg/s', 'accum  fess in-cloud source')
    call addfld ('fesssrc_c2', (/ 'lev' /), 'A', 'kg/kg/s', 'Aitken fess in-cloud source')
    call addfld ('fesssrc_c3', (/ 'lev' /), 'A', 'kg/kg/s', 'coarse fess in-cloud source')
    call addfld ('fesmsrc_c1', (/ 'lev' /), 'A', 'kg/kg/s', 'accum  fesm in-cloud source')
    call addfld ('fesmsrc_c2', (/ 'lev' /), 'A', 'kg/kg/s', 'Aitken fesm in-cloud source')
    call addfld ('fesmsrc_c3', (/ 'lev' /), 'A', 'kg/kg/s', 'coarse fesm in-cloud source')
    call addfld ('fescsrc_c1', (/ 'lev' /), 'A', 'kg/kg/s', 'accum  fesc in-cloud source')
    call addfld ('fescsrc_c2', (/ 'lev' /), 'A', 'kg/kg/s', 'Aitken fesc in-cloud source')
    call addfld ('fescsrc_c3', (/ 'lev' /), 'A', 'kg/kg/s', 'coarse fesc in-cloud source')
    call addfld ('fesfsrc_c1', (/ 'lev' /), 'A', 'kg/kg/s', 'accum  fesf in-cloud source')
    call addfld ('fesfsrc_c2', (/ 'lev' /), 'A', 'kg/kg/s', 'Aitken fesf in-cloud source')
    call addfld ('fesfsrc_c3', (/ 'lev' /), 'A', 'kg/kg/s', 'coarse fesf in-cloud source')
    call addfld ('fesswrk_c1', horiz_only,    'A', 'kg/m2/s', 'accum  fesc in-cloud column source')
    call addfld ('fesswrk_c2', horiz_only,    'A', 'kg/m2/s', 'Aitken fesc in-cloud column source')
    call addfld ('fesswrk_c3', horiz_only,    'A', 'kg/m2/s', 'coarse fesc in-cloud column source')
    call addfld ('fesmwrk_c1', horiz_only,    'A', 'kg/m2/s', 'accum  fesc in-cloud column source')
    call addfld ('fesmwrk_c2', horiz_only,    'A', 'kg/m2/s', 'Aitken fesc in-cloud column source')
    call addfld ('fesmwrk_c3', horiz_only,    'A', 'kg/m2/s', 'coarse fesc in-cloud column source')
    call addfld ('fescwrk_c1', horiz_only,    'A', 'kg/m2/s', 'accum  fesc in-cloud column source')
    call addfld ('fescwrk_c2', horiz_only,    'A', 'kg/m2/s', 'Aitken fesc in-cloud column source')
    call addfld ('fescwrk_c3', horiz_only,    'A', 'kg/m2/s', 'coarse fesc in-cloud column source')
    call addfld ('fesfwrk_c1', horiz_only,    'A', 'kg/m2/s', 'accum  fesf in-cloud column source')
    call addfld ('fesfwrk_c2', horiz_only,    'A', 'kg/m2/s', 'Aitken fesf in-cloud column source')
    call addfld ('fesfwrk_c3', horiz_only,    'A', 'kg/m2/s', 'coarse fesf in-cloud column source')
    call addfld ('fesssrc_a1', (/ 'lev' /), 'A', 'kg/kg/s', 'accum  fess interstitial source')
    call addfld ('fesssrc_a2', (/ 'lev' /), 'A', 'kg/kg/s', 'Aitken fess interstitial source')
    call addfld ('fesssrc_a3', (/ 'lev' /), 'A', 'kg/kg/s', 'coarse fess interstitial source')
    call addfld ('fesmsrc_a1', (/ 'lev' /), 'A', 'kg/kg/s', 'accum  fesm interstitial source')
    call addfld ('fesmsrc_a2', (/ 'lev' /), 'A', 'kg/kg/s', 'Aitken fesm interstitial source')
    call addfld ('fesmsrc_a3', (/ 'lev' /), 'A', 'kg/kg/s', 'coarse fesm interstitial source')
    call addfld ('fescsrc_a1', (/ 'lev' /), 'A', 'kg/kg/s', 'accum  fesc interstitial source')
    call addfld ('fescsrc_a2', (/ 'lev' /), 'A', 'kg/kg/s', 'Aitken fesc interstitial source')
    call addfld ('fescsrc_a3', (/ 'lev' /), 'A', 'kg/kg/s', 'coarse fesc interstitial source')
    call addfld ('fesfsrc_a1', (/ 'lev' /), 'A', 'kg/kg/s', 'accum  fesf interstitial source')
    call addfld ('fesfsrc_a2', (/ 'lev' /), 'A', 'kg/kg/s', 'Aitken fesf interstitial source')
    call addfld ('fesfsrc_a3', (/ 'lev' /), 'A', 'kg/kg/s', 'coarse fesf interstitial source')
    call addfld ('fesswrk_a1', horiz_only,    'A', 'kg/m2/s', 'accum  fess interstitial column source')
    call addfld ('fesswrk_a2', horiz_only,    'A', 'kg/m2/s', 'Aitken fess interstitial column source')
    call addfld ('fesswrk_a3', horiz_only,    'A', 'kg/m2/s', 'coarse fess interstitial column source')
    call addfld ('fesmwrk_a1', horiz_only,    'A', 'kg/m2/s', 'accum  fesm interstitial column source')
    call addfld ('fesmwrk_a2', horiz_only,    'A', 'kg/m2/s', 'Aitken fesm interstitial column source')
    call addfld ('fesmwrk_a3', horiz_only,    'A', 'kg/m2/s', 'coarse fesm interstitial column source')
    call addfld ('fescwrk_a1', horiz_only,    'A', 'kg/m2/s', 'accum  fesc interstitial column source')
    call addfld ('fescwrk_a2', horiz_only,    'A', 'kg/m2/s', 'Aitken fesc interstitial column source')
    call addfld ('fescwrk_a3', horiz_only,    'A', 'kg/m2/s', 'coarse fesc interstitial column source')
    call addfld ('fesfwrk_a1', horiz_only,    'A', 'kg/m2/s', 'accum  fesf interstitial column source')
    call addfld ('fesfwrk_a2', horiz_only,    'A', 'kg/m2/s', 'Aitken fesf interstitial column source')
    call addfld ('fesfwrk_a3', horiz_only,    'A', 'kg/m2/s', 'coarse fesf interstitial column source')

    call addfld ('myPH_a1',(/ 'lev' /), 'A','pH','my pH')
    call addfld ('myPH_a2',(/ 'lev' /), 'A','pH','my pH')
    call addfld ('myPH_a3',(/ 'lev' /), 'A','pH','my pH')
    call addfld ('MOXASRC',(/ 'lev' /), 'A','umol/L','oxalate source')
    call addfld ('MOXASOA',(/ 'lev' /), 'A','mol/mol','oxalate source')
    call addfld ('mbar',(/ 'lev' /), 'A','amu','mean wet atmospheric mass')
    call addfld ('pdel',(/ 'lev' /), 'A','Pa' ,'pressure difference')

    if ( history_aerosol ) then
     call add_default ('fesssrc_c1',1,' ')
     call add_default ('fesssrc_c2',1,' ')
     call add_default ('fesssrc_c3',1,' ')
     call add_default ('fesmsrc_c1',1,' ')
     call add_default ('fesmsrc_c2',1,' ')
     call add_default ('fesmsrc_c3',1,' ')
     call add_default ('fescsrc_c1',1,' ')
     call add_default ('fescsrc_c2',1,' ')
     call add_default ('fescsrc_c3',1,' ')
     call add_default ('fesfsrc_c1',1,' ')
     call add_default ('fesfsrc_c2',1,' ')
     call add_default ('fesfsrc_c3',1,' ')
     call add_default ('fesswrk_c1',1,' ')
     call add_default ('fesswrk_c2',1,' ')
     call add_default ('fesswrk_c3',1,' ')
     call add_default ('fesmwrk_c1',1,' ')
     call add_default ('fesmwrk_c2',1,' ')
     call add_default ('fesmwrk_c3',1,' ')
     call add_default ('fescwrk_c1',1,' ')
     call add_default ('fescwrk_c2',1,' ')
     call add_default ('fescwrk_c3',1,' ')
     call add_default ('fesfwrk_c1',1,' ')
     call add_default ('fesfwrk_c2',1,' ')
     call add_default ('fesfwrk_c3',1,' ')
     call add_default ('fesssrc_a1',1,' ')
     call add_default ('fesssrc_a2',1,' ')
     call add_default ('fesssrc_a3',1,' ')
     call add_default ('fesmsrc_a1',1,' ')
     call add_default ('fesmsrc_a2',1,' ')
     call add_default ('fesmsrc_a3',1,' ')
     call add_default ('fescsrc_a1',1,' ')
     call add_default ('fescsrc_a2',1,' ')
     call add_default ('fescsrc_a3',1,' ')
     call add_default ('fesfsrc_a1',1,' ')
     call add_default ('fesfsrc_a2',1,' ')
     call add_default ('fesfsrc_a3',1,' ')
     call add_default ('fesswrk_a1',1,' ')
     call add_default ('fesswrk_a2',1,' ')
     call add_default ('fesswrk_a3',1,' ')
     call add_default ('fesmwrk_a1',1,' ')
     call add_default ('fesmwrk_a2',1,' ')
     call add_default ('fesmwrk_a3',1,' ')
     call add_default ('fescwrk_a1',1,' ')
     call add_default ('fescwrk_a2',1,' ')
     call add_default ('fescwrk_a3',1,' ')
     call add_default ('fesfwrk_a1',1,' ')
     call add_default ('fesfwrk_a2',1,' ')
     call add_default ('fesfwrk_a3',1,' ')
     call add_default ('myPH_a1', 1, ' ')
     call add_default ('myPH_a2', 1, ' ')
     call add_default ('myPH_a3', 1, ' ')
     call add_default ('MOXASRC', 1, ' ')
     call add_default ('MOXASOA', 1, ' ')
     call add_default ('mbar', 1, ' ')
     call add_default ('pdel', 1, ' ')
    endif

  end subroutine fe_cldaero_init

!----------------------------------------------------------------------------------
! Update the mixing ratios
!----------------------------------------------------------------------------------

  subroutine fe_cldaero_update( &
       ncol, lchnk, loffset, dtime, mbar, pdel, press, tfld, cldnum, cldfrc, cfact, xlwc, &
       xso4c, qcw, qin )

    ! args 

    integer,  intent(in) :: ncol
    integer,  intent(in) :: lchnk ! chunk id
    integer,  intent(in) :: loffset

    real(r8), intent(in) :: dtime ! time step (sec)

    real(r8), intent(in) :: mbar(:,:) ! mean wet atmospheric mass ( amu )
    real(r8), intent(in) :: pdel(:,:) 
    real(r8), intent(in) :: press(:,:)
    real(r8), intent(in) :: tfld(:,:) ! temp

    real(r8), intent(in) :: cldnum(:,:)
    real(r8), intent(in) :: cldfrc(:,:)
    real(r8), intent(in) :: cfact(:,:)
    real(r8), intent(in) :: xlwc(:,:)

    real(r8), intent(in) :: xso4c(:,:)
    real(r8), intent(inout) :: qcw(:,:,:) ! cloud-borne aerosol (vmr)
    real(r8), intent(inout) :: qin(:,:,:) ! xported species ( vmr )


    ! local vars ...

    integer :: l, n, m
    integer :: i,k
    real(r8) :: xl

    ! dsh add for iron
    
    real(r8), dimension(pcols,pver,ntot_amode) :: dqdt_fets
    real(r8), dimension(pcols,pver,ntot_amode) :: dqdt_fetm
    real(r8), dimension(pcols,pver,ntot_amode) :: dqdt_fetc
    real(r8), dimension(pcols,pver,ntot_amode) :: dqdt_fetf
    real(r8), dimension(pcols,pver,ntot_amode) :: dqdt_fess
    real(r8), dimension(pcols,pver,ntot_amode) :: dqdt_fesm
    real(r8), dimension(pcols,pver,ntot_amode) :: dqdt_fesc
    real(r8), dimension(pcols,pver,ntot_amode) :: dqdt_fesf
    real(r8), dimension(pcols,pver,ntot_amode) :: myPH
    real(r8), dimension(pcols,pver) :: xoc
    real(r8), dimension(pcols,pver) :: xsoa
    real(r8), dimension(pcols,pver) :: xocsoa
    real(r8), dimension(pcols,pver) :: molarc2o4
    real(r8), dimension(pcols,ntot_amode) :: fesswrk_c
    real(r8), dimension(pcols,ntot_amode) :: fesmwrk_c
    real(r8), dimension(pcols,ntot_amode) :: fescwrk_c
    real(r8), dimension(pcols,ntot_amode) :: fesfwrk_c
    real(r8), dimension(pcols,ntot_amode) :: fesswrk_a
    real(r8), dimension(pcols,ntot_amode) :: fesmwrk_a
    real(r8), dimension(pcols,ntot_amode) :: fescwrk_a
    real(r8), dimension(pcols,ntot_amode) :: fesfwrk_a
    real(r8)  k_fast, k_slow, A_fast, A_slow, m_fast, m_slow, R_fast, R_slow
    real(r8)  aH1f, aH1c, aH3, Ca, S4 
  
    ! make sure dqdt is zero initially, for budgets
    dqdt_fess(:,:,:) = 0.0_r8
    dqdt_fesm(:,:,:) = 0.0_r8
    dqdt_fets(:,:,:) = 0.0_r8
    dqdt_fetm(:,:,:) = 0.0_r8
    dqdt_fesc(:,:,:) = 0.0_r8
    dqdt_fetc(:,:,:) = 0.0_r8
    dqdt_fesf(:,:,:) = 0.0_r8
    dqdt_fetf(:,:,:) = 0.0_r8
    xsoa(:,:)        = 0.0_r8   !secondary organic carbon - carried as C(mw=12) in model
    molarc2o4(:,:)   = 0.0_r8   !Scanza 2018 calc from [SO4] or [OC+SOA], updated to be [SOA] only here
    fesswrk_c(:,:)   = 0.0_r8   !Column total tendancy - cloudbourne
    fesmwrk_c(:,:)   = 0.0_r8   !Column total tendancy - cloudbourne
    fescwrk_c(:,:)   = 0.0_r8   !Column total tendancy - cloudbourne
    fesfwrk_c(:,:)   = 0.0_r8   !Column total tendancy - cloudbourne

    do k = 1, pver
      do i = 1, ncol
        do n = 1, ntot_amode
           l = lptr2_soa_cw_amode(n,nsoa) - loffset
           if ( l > 0 ) then 
            xsoa(i,k) = xsoa(i,k) + qcw(i,k,l)
           endif
           l = lptr2_soa_a_amode(n,nsoa) - loffset
           if ( l > 0 ) then
            xsoa(i,k) = xsoa(i,k) + qin(i,k,l)
           endif
        end do
      end do
    end do

    lev_loop: do k = 1,pver
       col_loop: do i = 1,ncol
          cloud: if (cldfrc(i,k) >= 1.0e-5_r8) then
             xl = xlwc(i,k) 

             IF (XL .ge. 1.e-8_r8) THEN !! WHEN CLOUD IS PRESENTED

                molarc2o4(i,k) = 150._r8*(xsoa(i,k)/1.4e-8_r8)           !dsh updated for oxalate SOA only umol/L      

                k_slow = min(9.532e-9*molarc2o4(i,k)+2.975e-8,1.727e-7)  ! min value ensures that the rate is not calc  
                k_fast = min(2.3e-7*molarc2o4(i,k)+4.8e-7,3.96e-6)       ! for a conc of C2O4 > 15

                ! compute tendencies for Fe aerosol-in-cloud-water
                do n = 1, ntot_amode
                  
                  l = lptr_fets_cw_amode(n) - loffset
                  if ( l > 0 ) then
                   if ((qcw(i,k,l)) .gt. 0._r8) then
                     dqdt_fess(i,k,n) = min(k_slow*(qcw(i,k,l)),(qcw(i,k,l))/dtime) 
                     dqdt_fets(i,k,n) = -dqdt_fess(i,k,n) 
                   else
                     dqdt_fess(i,k,n) = 0._r8
                     dqdt_fets(i,k,n) = 0._r8
                   endif                      
                  endif             ! units (mol/mol)Fe/s         

                  l = lptr_fetm_cw_amode(n) - loffset
                  if ( l > 0 ) then
                   if ((qcw(i,k,l)) .gt. 0._r8) then
                     dqdt_fesm(i,k,n) = min(k_fast*(qcw(i,k,l)),(qcw(i,k,l))/dtime)  
                     dqdt_fetm(i,k,n) = -dqdt_fesm(i,k,n)
                   else
                    dqdt_fesm(i,k,n) = 0._r8
                    dqdt_fetm(i,k,n) = 0._r8
                   endif                      
                  endif             ! units (mol/mol)Fe/s         
                 
                  l = lptr_fetc_cw_amode(n) - loffset
                  if ( l > 0 ) then
                   if ((qcw(i,k,l)) .gt. 0._r8) then
                     dqdt_fesc(i,k,n) = min(k_fast*(qcw(i,k,l)),(qcw(i,k,l))/dtime) 
                     dqdt_fetc(i,k,n) = -dqdt_fesc(i,k,n)
                   else
                    dqdt_fesc(i,k,n) = 0._r8
                    dqdt_fetc(i,k,n) = 0._r8
                   endif
                  endif             ! units (mol/mol)Fe/s                 

                  l = lptr_fetf_cw_amode(n) - loffset
                  if ( l > 0 ) then
                   if ((qcw(i,k,l)) .gt. 0._r8) then
                     dqdt_fesf(i,k,n) = min(k_fast*(qcw(i,k,l)),(qcw(i,k,l))/dtime)
                     dqdt_fetf(i,k,n) = -dqdt_fesf(i,k,n)
                   else
                    dqdt_fesf(i,k,n) = 0._r8
                    dqdt_fetf(i,k,n) = 0._r8
                   endif
                  endif             ! units (mol/mol)Fe/s   
                end do

             ENDIF !! WHEN CLOUD IS PRESENTED
          endif cloud
       enddo col_loop
    enddo lev_loop
    
!   Tendency due to aqueous chemistry 

   do n = 1, ntot_amode
      m = lptr_fetc_cw_amode(n)
      l = m - loffset
      if (l > 0) then
       do k=1,pver
        do i=1,ncol
         fesswrk_c(i,n) = fesswrk_c(i,n) + dqdt_fess(i,k,n) * adv_mass(l)/mbar(i,k)*pdel(i,k)/gravit
         fesmwrk_c(i,n) = fesmwrk_c(i,n) + dqdt_fesm(i,k,n) * adv_mass(l)/mbar(i,k)*pdel(i,k)/gravit
         fescwrk_c(i,n) = fescwrk_c(i,n) + dqdt_fesc(i,k,n) * adv_mass(l)/mbar(i,k)*pdel(i,k)/gravit
         fesfwrk_c(i,n) = fesfwrk_c(i,n) + dqdt_fesf(i,k,n) * adv_mass(l)/mbar(i,k)*pdel(i,k)/gravit
        end do
       end do
      endif 
   end do

! output fields
    call outfld('fesssrc_c1', dqdt_fess(1,1,1), pcols, lchnk)
    call outfld('fesssrc_c2', dqdt_fess(1,1,2), pcols, lchnk)
    call outfld('fesssrc_c3', dqdt_fess(1,1,3), pcols, lchnk)
    call outfld('fesmsrc_c1', dqdt_fesm(1,1,1), pcols, lchnk)
    call outfld('fesmsrc_c2', dqdt_fesm(1,1,2), pcols, lchnk)
    call outfld('fesmsrc_c3', dqdt_fesm(1,1,3), pcols, lchnk)
    call outfld('fescsrc_c1', dqdt_fesc(1,1,1), pcols, lchnk)
    call outfld('fescsrc_c2', dqdt_fesc(1,1,2), pcols, lchnk)
    call outfld('fescsrc_c3', dqdt_fesc(1,1,3), pcols, lchnk)
    call outfld('fesfsrc_c1', dqdt_fesf(1,1,1), pcols, lchnk)
    call outfld('fesfsrc_c2', dqdt_fesf(1,1,2), pcols, lchnk)
    call outfld('fesfsrc_c3', dqdt_fesf(1,1,3), pcols, lchnk)
    call outfld('fesswrk_c1', fesswrk_c(1,1), pcols, lchnk)
    call outfld('fesswrk_c2', fesswrk_c(1,2), pcols, lchnk)
    call outfld('fesswrk_c3', fesswrk_c(1,3), pcols, lchnk)
    call outfld('fesmwrk_c1', fesmwrk_c(1,1), pcols, lchnk)
    call outfld('fesmwrk_c2', fesmwrk_c(1,2), pcols, lchnk)
    call outfld('fesmwrk_c3', fesmwrk_c(1,3), pcols, lchnk)
    call outfld('fescwrk_c1', fescwrk_c(1,1), pcols, lchnk)
    call outfld('fescwrk_c2', fescwrk_c(1,2), pcols, lchnk)
    call outfld('fescwrk_c3', fescwrk_c(1,3), pcols, lchnk)
    call outfld('fesfwrk_c1', fesfwrk_c(1,1), pcols, lchnk)
    call outfld('fesfwrk_c2', fesfwrk_c(1,2), pcols, lchnk)
    call outfld('fesfwrk_c3', fesfwrk_c(1,3), pcols, lchnk)
    call outfld('MOXASRC'   , molarc2o4       , pcols, lchnk)
    call outfld('MOXASOA'   , xsoa            , pcols, lchnk)


    !==============================================================
    ! ... Update the in-cloud mixing ratios
    !==============================================================
    do k = 1,pver
       do n = 1, ntot_amode

         !write(iulog,*)'starting update'

         l = lptr_fess_cw_amode(n) - loffset
           if ( l > 0 ) then
             qcw(:ncol,k,l) = qcw(:ncol,k,l) + dqdt_fess(:ncol,k,n)*dtime
             qcw(:,k,l) = MAX(qcw(:,k,l), small_value )
           endif  

         l = lptr_fesm_cw_amode(n) - loffset
           if ( l > 0 ) then
             qcw(:ncol,k,l) = qcw(:ncol,k,l) + dqdt_fesm(:ncol,k,n)*dtime
             qcw(:,k,l) = MAX(qcw(:,k,l), small_value )
           endif  

          l = lptr_fets_cw_amode(n) - loffset
           if ( l > 0 ) then
             qcw(:ncol,k,l) = qcw(:ncol,k,l) + dqdt_fets(:ncol,k,n)*dtime
             qcw(:,k,l) = MAX(qcw(:,k,l), small_value )
           endif  

          l = lptr_fetm_cw_amode(n) - loffset
           if ( l > 0 ) then
             qcw(:ncol,k,l) = qcw(:ncol,k,l) + dqdt_fetm(:ncol,k,n)*dtime
             qcw(:,k,l) = MAX(qcw(:,k,l), small_value )
           endif  

          l = lptr_fetc_cw_amode(n) - loffset
           if ( l > 0 ) then
             qcw(:ncol,k,l) = qcw(:ncol,k,l) + dqdt_fetc(:ncol,k,n)*dtime
             qcw(:,k,l) = MAX(qcw(:,k,l), small_value )
           endif

          l = lptr_fesc_cw_amode(n) - loffset
           if ( l > 0 ) then
             qcw(:ncol,k,l) = qcw(:ncol,k,l) + dqdt_fesc(:ncol,k,n)*dtime
             qcw(:,k,l) = MAX(qcw(:,k,l), small_value )
           endif

          l = lptr_fetf_cw_amode(n) - loffset
           if ( l > 0 ) then
             qcw(:ncol,k,l) = qcw(:ncol,k,l) + dqdt_fetf(:ncol,k,n)*dtime
             qcw(:,k,l) = MAX(qcw(:,k,l), small_value )
           endif

          l = lptr_fesf_cw_amode(n) - loffset
           if ( l > 0 ) then
             qcw(:ncol,k,l) = qcw(:ncol,k,l) + dqdt_fesf(:ncol,k,n)*dtime
             qcw(:,k,l) = MAX(qcw(:,k,l), small_value )
           endif
       end do
     end do

    ! reset dqdt to zero initially, for budgets
    ! set Ca, S4, myPH to zero initially, for budgets
    dqdt_fess(:,:,:) = 0.0_r8
    dqdt_fesm(:,:,:) = 0.0_r8
    dqdt_fets(:,:,:) = 0.0_r8
    dqdt_fetm(:,:,:) = 0.0_r8
    dqdt_fesc(:,:,:) = 0.0_r8
    dqdt_fetc(:,:,:) = 0.0_r8
    dqdt_fesf(:,:,:) = 0.0_r8
    dqdt_fetf(:,:,:) = 0.0_r8
    Ca = 0.0_r8
    S4 = 0.0_r8
    myPH(:,:,:) = 0.0_r8
    fesswrk_a(:,:)   = 0.0_r8   !Column total tendancy - interstitial
    fesmwrk_a(:,:)   = 0.0_r8   !Column total tendancy - interstitial
    fescwrk_a(:,:)   = 0.0_r8   !Column total tendancy - interstitial
    fesfwrk_a(:,:)   = 0.0_r8   !Column total tendancy - interstitial
    
    !! calculate dqdt for acid processing of interstitial aerosol

    A_fast = 90._r8         ! specific surface area m2/g (=illite)
    A_slow = 100._r8        ! specific surface area m2/g (=hematite)
    m_fast = 0.39_r8        ! empirical reaction order
    m_slow = 0.5_r8         ! empirical reaction order
    aH1f   = 0.1_r8         ! pH = 1 fine mode low limit
    aH1c   = 0.01_r8        ! pH = 2 coarse mode low limit
    aH3    = 3.162e-8_r8    ! pH = 7.5 fine and coarse mode high limit

    do k = 1,pver
     do i  = 1,ncol 
       k_fast = 1.3e-11_r8*exp(6.7e3_r8*(1._r8/298._r8 - 1._r8/(tfld(i,k))))
       k_slow = 1.8e-11_r8*exp(9.2e3_r8*(1._r8/298._r8 - 1._r8/(tfld(i,k))))
       do n = 1, ntot_amode

          l = lptr_dust6_a_amode(n) - loffset
          if ( l > 0 ) then      
            Ca = qin(i,k,l)
          endif
          l = lptr_so4_a_amode(n) - loffset
          if ( l > 0 ) then
            S4 = qin(i,k,l)
            if (S4 .ge. Ca .and. n .lt. 3) then              ! Fine (Aitken and Accum) modes
              R_fast =k_fast*A_fast*389.34*aH1f**m_fast
              R_slow =k_slow*A_slow*159.59*aH1f**m_slow
              myPH(i,k,n) = 1.0_r8
            else if (S4 .ge. Ca .and. n .gt. 2) then         ! Coarse mode
              R_fast =k_fast*A_fast*389.34*aH1c**m_fast
              R_slow =k_slow*A_slow*159.59*aH1c**m_slow
              myPH(i,k,n) = 2.0_r8
            else
              R_fast =k_fast*A_fast*389.34*aH3**m_fast
              R_slow =k_slow*A_slow*159.59*aH3**m_slow
              myPH(i,k,n) = 7.5_r8
            endif
         endif   
       
         l = lptr_fets_a_amode(n) - loffset
         if ( l > 0 ) then
            if ((qin(i,k,l)) .gt. 0._r8) then
               dqdt_fess(i,k,n) = min(R_slow*(qin(i,k,l)),(qin(i,k,l))/dtime) 
               dqdt_fets(i,k,n) = -dqdt_fess(i,k,n)
            else
               dqdt_fess(i,k,n) = 0._r8
               dqdt_fets(i,k,n) = 0._r8
            endif
          endif             ! units (mol/mol)Fe/s         

         l = lptr_fetm_a_amode(n) - loffset
         if ( l > 0 ) then
            if ((qin(i,k,l)) .gt. 0._r8) then
               dqdt_fesm(i,k,n) = min(R_fast*(qin(i,k,l)),(qin(i,k,l))/dtime) 
               dqdt_fetm(i,k,n) = -dqdt_fesm(i,k,n)
            else
               dqdt_fesm(i,k,n) = 0._r8
               dqdt_fetm(i,k,n) = 0._r8
            endif
            endif             ! units (mol/mol)Fe/s         
   
         l = lptr_fetc_a_amode(n) - loffset
         if ( l > 0 ) then
            if ((qin(i,k,l)) .gt. 0._r8) then
               dqdt_fesc(i,k,n) = min(R_fast*(qin(i,k,l)),(qin(i,k,l))/dtime)
               dqdt_fetc(i,k,n) = -dqdt_fesc(i,k,n)
            else
               dqdt_fesc(i,k,n) = 0._r8
               dqdt_fetc(i,k,n) = 0._r8
            endif
            endif             ! units (mol/mol)Fe/s         

         l = lptr_fetf_a_amode(n) - loffset
         if ( l > 0 ) then
            if ((qin(i,k,l)) .gt. 0._r8) then
               dqdt_fesf(i,k,n) = min(R_fast*(qin(i,k,l)),(qin(i,k,l))/dtime)
               dqdt_fetf(i,k,n) = -dqdt_fesf(i,k,n)
            else
               dqdt_fesf(i,k,n) = 0._r8
               dqdt_fetf(i,k,n) = 0._r8
            endif
            endif             ! units (mol/mol)Fe/s   

        end do !modes
      end do ! ncol
    end do  ! plev

!   Tendency due to aqueous chemistry 

   do n = 1, ntot_amode
      m = lptr_fetc_a_amode(n)
      l = m - loffset
      if (l > 0) then
       do k=1,pver
        do i=1,ncol
         fesswrk_a(i,n) = fesswrk_a(i,n) + dqdt_fess(i,k,n) * adv_mass(l)/mbar(i,k)*pdel(i,k)/gravit
         fesmwrk_a(i,n) = fesmwrk_a(i,n) + dqdt_fesm(i,k,n) * adv_mass(l)/mbar(i,k)*pdel(i,k)/gravit
         fescwrk_a(i,n) = fescwrk_a(i,n) + dqdt_fesc(i,k,n) * adv_mass(l)/mbar(i,k)*pdel(i,k)/gravit
         fesfwrk_a(i,n) = fesfwrk_a(i,n) + dqdt_fesf(i,k,n) * adv_mass(l)/mbar(i,k)*pdel(i,k)/gravit
        end do
       end do
      endif 
   end do

! output fields
    call outfld('fesssrc_a1', dqdt_fess(1,1,1), pcols, lchnk)
    call outfld('fesssrc_a2', dqdt_fess(1,1,2), pcols, lchnk)
    call outfld('fesssrc_a3', dqdt_fess(1,1,3), pcols, lchnk)
    call outfld('fesmsrc_a1', dqdt_fesm(1,1,1), pcols, lchnk)
    call outfld('fesmsrc_a2', dqdt_fesm(1,1,2), pcols, lchnk)
    call outfld('fesmsrc_a3', dqdt_fesm(1,1,3), pcols, lchnk)
    call outfld('fescsrc_a1', dqdt_fesc(1,1,1), pcols, lchnk)
    call outfld('fescsrc_a2', dqdt_fesc(1,1,2), pcols, lchnk)
    call outfld('fescsrc_a3', dqdt_fesc(1,1,3), pcols, lchnk)
    call outfld('fesfsrc_a1', dqdt_fesf(1,1,1), pcols, lchnk)
    call outfld('fesfsrc_a2', dqdt_fesf(1,1,2), pcols, lchnk)
    call outfld('fesfsrc_a3', dqdt_fesf(1,1,3), pcols, lchnk)
    call outfld('fesswrk_a1', fesswrk_a(1,1), pcols, lchnk)
    call outfld('fesswrk_a2', fesswrk_a(1,2), pcols, lchnk)
    call outfld('fesswrk_a3', fesswrk_a(1,3), pcols, lchnk)
    call outfld('fesmwrk_a1', fesmwrk_a(1,1), pcols, lchnk)
    call outfld('fesmwrk_a2', fesmwrk_a(1,2), pcols, lchnk)
    call outfld('fesmwrk_a3', fesmwrk_a(1,3), pcols, lchnk)
    call outfld('fescwrk_a1', fescwrk_a(1,1), pcols, lchnk)
    call outfld('fescwrk_a2', fescwrk_a(1,2), pcols, lchnk)
    call outfld('fescwrk_a3', fescwrk_a(1,3), pcols, lchnk)
    call outfld('fesfwrk_a1', fesfwrk_a(1,1), pcols, lchnk)
    call outfld('fesfwrk_a2', fesfwrk_a(1,2), pcols, lchnk)
    call outfld('fesfwrk_a3', fesfwrk_a(1,3), pcols, lchnk)
    call outfld('myPH_a1'   , myPH(1,1,1)     , pcols, lchnk)
    call outfld('myPH_a2'   , myPH(1,1,2)     , pcols, lchnk)
    call outfld('myPH_a3'   , myPH(1,1,3)     , pcols, lchnk)
    call outfld('mbar'      , mbar     , pcols, lchnk)
    call outfld('pdel'      , pdel     , pcols, lchnk)

    !==============================================================
    ! ... Update the interstitial mixing ratios
    !==============================================================
    do k = 1,pver
       do n = 1, ntot_amode

         l = lptr_fess_a_amode(n) - loffset
           if ( l > 0 ) then
             qin(:ncol,k,l) = qin(:ncol,k,l) + dqdt_fess(:ncol,k,n)*dtime
             qin(:,k,l) = MAX(qin(:,k,l), small_value )
           endif

         l = lptr_fesm_a_amode(n) - loffset
           if ( l > 0 ) then
             qin(:ncol,k,l) = qin(:ncol,k,l) + dqdt_fesm(:ncol,k,n)*dtime
             qin(:,k,l) = MAX(qin(:,k,l), small_value )
           endif

          l = lptr_fets_a_amode(n) - loffset
           if ( l > 0 ) then
             qin(:ncol,k,l) = qin(:ncol,k,l) + dqdt_fets(:ncol,k,n)*dtime
             qin(:,k,l) = MAX(qin(:,k,l), small_value )
           endif

          l = lptr_fetm_a_amode(n) - loffset
           if ( l > 0 ) then
             qin(:ncol,k,l) = qin(:ncol,k,l) + dqdt_fetm(:ncol,k,n)*dtime
             qin(:,k,l) = MAX(qin(:,k,l), small_value )
           endif

          l = lptr_fetc_a_amode(n) - loffset
           if ( l > 0 ) then
             qin(:ncol,k,l) = qin(:ncol,k,l) + dqdt_fetc(:ncol,k,n)*dtime
             qin(:,k,l) = MAX(qin(:,k,l), small_value )
           endif

          l = lptr_fesc_a_amode(n) - loffset
           if ( l > 0 ) then
             qin(:ncol,k,l) = qin(:ncol,k,l) + dqdt_fesc(:ncol,k,n)*dtime
             qin(:,k,l) = MAX(qin(:,k,l), small_value )
           endif

          l = lptr_fetf_a_amode(n) - loffset
           if ( l > 0 ) then
             qin(:ncol,k,l) = qin(:ncol,k,l) + dqdt_fetf(:ncol,k,n)*dtime
             qin(:,k,l) = MAX(qin(:,k,l), small_value )
           endif

          l = lptr_fesf_a_amode(n) - loffset
           if ( l > 0 ) then
             qin(:ncol,k,l) = qin(:ncol,k,l) + dqdt_fesf(:ncol,k,n)*dtime
             qin(:,k,l) = MAX(qin(:,k,l), small_value )
           endif
       end do
     end do

  end subroutine fe_cldaero_update

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------

end module fe_cldaero_mod
