! radiation_ifs_rrtmgp.F90 - Interface to IFS implementation of RRTMGP-NN, 
! which can use either RRTMGP k-distributions loaded from a netCDF file,
! or faster NN emulators of these distributions
!
! (C) Copyright 2015- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
! Author:  Peter U kkonen
! Email:   peterukk@gmail.com
!
! Modifications
!   2020-01-30  P. Ukkonen  Initial version
!   2022-02-01  P. Ukkonen  Reworked to incur minimal changes to existing ecRAD code

module radiation_ifs_rrtmgp

  implicit none
  private  
  public  :: setup_gas_optics_ifs_rrtmgp, gas_optics_ifs_rrtmgp, set_gas_units

  ! List of character for case conversion, need this for changing an entry
  ! in a string array
  character(len=26), parameter :: LOWER_CASE_CHARS = 'abcdefghijklmnopqrstuvwxyz'
  character(len=26), parameter :: UPPER_CASE_CHARS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

contains
  
  !---------------------------------------------------------------------
  ! Setup the IFS implementation of RRTMGP gas absorption model
  subroutine setup_gas_optics_ifs_rrtmgp(config)

    use parkind1,         only : jprb
    use yomhook,          only : lhook, dr_hook, jphook
    use radiation_io,     only : nulerr, radiation_abort    
    use radiation_config, only : config_type, IGasModelRRTMGP_NN,IGasModelRRTMGP, ISolverSPARTACUS
    use radiation_spectral_definition, only &
         &  : SolarReferenceTemperature, TerrestrialReferenceTemperature
    use mo_load_coefficients,     only:   load_and_init
    use mo_gas_concentrations,    only :  ty_gas_concs
    use mo_gas_optics_rrtmgp,     only :  ty_gas_optics_rrtmgp

    type(config_type),      intent(inout), target    :: config

    integer :: irep ! For implied do
    integer :: idx_hcfc22
    logical :: use_neural_nets =.false.

    integer, parameter :: RRTMGP_GPOINT_REORDERING_LW_256(256)  = (/ & 
        & 225, 241, 226, 242, 227, 193, 228, 113, 114, 115, 97, 116, 194, 98, 81, 117, 229,  & 
        & 243, 129, 82, 83, 118, 99, 84, 195, 85, 230, 86, 87, 100, 65, 119, 231, 177, 244, 88,  & 
        & 66, 196, 130, 101, 232, 67, 102, 245, 233, 120, 89, 178, 197, 68, 131, 234, 103, 235,  & 
        & 246, 198, 236, 90, 69, 179, 237, 33, 104, 91, 132, 238, 121, 34, 199, 92, 247, 239, 35,  & 
        & 70, 93, 180, 133, 36, 240, 94, 105, 200, 122, 17, 95, 209, 37, 71, 123, 181, 248, 96,  & 
        & 134, 124, 106, 38, 145, 201, 18, 107, 182, 72, 49, 135, 108, 125, 39, 146, 109, 210,  & 
        & 202, 249, 161, 19, 110, 203, 183, 147, 50, 126, 204, 111, 73, 136, 148, 112, 51, 162,  & 
        & 40, 205, 20, 250, 1, 127, 52, 128, 251, 163, 149, 211, 184, 74, 206, 252, 21, 53, 164,  & 
        & 75, 253, 41, 137, 150, 54, 254, 76, 22, 207, 165, 2, 255, 212, 77, 256, 208, 185, 42,  & 
        & 55, 78, 151, 166, 3, 23, 138, 43, 213, 139, 4, 79, 44, 56, 167, 186, 214, 140, 152, 45,  & 
        & 5, 24, 80, 187, 141, 215, 46, 168, 6, 188, 57, 47, 142, 48, 153, 189, 216, 7, 25, 190,  & 
        & 169, 143, 58, 154, 155, 191, 144, 59, 156, 8, 170, 157, 192, 26, 171, 217, 60, 158, & 
        & 172, 61, 27, 159, 173, 160, 62, 174, 9, 63, 28, 175, 218, 64, 176, 219, 29, 220, 10, & 
        & 11, 221, 30, 12, 222, 13, 223, 14, 224, 31, 15, 32, 16 /)
    integer, parameter :: RRTMGP_GPOINT_REORDERING_SW_224(224)  = (/ & 
        & 97, 81, 65, 113, 129, 49, 130, 82, 114, 131, 98, 115, 50, 66, 145, 116, 132, 146, & 
        & 17, 117, 147, 51, 148, 1, 83, 118, 149, 133, 67, 150, 119, 161, 52, 162, 18, 151, 163, & 
        & 120, 2, 164, 99, 165, 134, 68, 152, 166, 53, 167, 168, 169, 177, 170, 19, 171, 178, & 
        & 172, 173, 179, 121, 69, 3, 180, 54, 174, 181, 84, 135, 182, 183, 184, 20, 185, 186, & 
        & 187, 188, 189, 190, 191, 192, 175, 153, 100, 193, 70, 122, 55, 194, 176, 136, 21, 123, & 
        & 33, 71, 195, 124, 4, 154, 56, 101, 22, 125, 155, 85, 72, 137, 126, 196, 34, 23, 102, & 
        & 127, 156, 57, 128, 86, 5, 73, 138, 197, 24, 103, 157, 139, 209, 58, 87, 35, 210, 198, & 
        & 74, 140, 59, 211, 75, 104, 199, 60, 25, 212, 6, 141, 200, 76, 158, 88, 36, 201, 202, & 
        & 203, 204, 205, 206, 207, 208, 213, 61, 77, 142, 214, 215, 216, 217, 37, 218, 219, 220, & 
        & 221, 222, 224, 223, 26, 105, 7, 27, 62, 78, 143, 28, 159, 38, 89, 29, 144, 63, 106, 30, & 
        & 160, 79, 39, 8, 64, 31, 107, 32, 90, 108, 80, 91, 40, 109, 92, 93, 110, 9, 94, 111, 41, & 
        & 112, 95, 96, 10, 42, 11, 43, 12, 44, 45, 13, 46, 47, 48, 14, 15, 16 /)
    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_ifs_rrtmgp:setup_gas_optics_ifs_rrtmgp',0,hook_handle)

    ! HCFC22 is called CFC22 in RRTMGP look up tables; we can change the name here
    idx_hcfc22 = string_loc_in_array('hcfc22',config%rrtmgp_gas_names)
    if (idx_hcfc22 /= -1) config%rrtmgp_gas_names(idx_hcfc22) = 'cfc22'

    ! Load k-distributions
    call load_and_init(config%k_dist_lw, trim(config%rrtmgp_gas_optics_file_name_lw), config%rrtmgp_gas_names)
    call load_and_init(config%k_dist_sw, trim(config%rrtmgp_gas_optics_file_name_sw), config%rrtmgp_gas_names)

    ! Cloud and aerosol properties can only be defined per band
    if (config%i_gas_model_sw == IGasModelRRTMGP_NN .or. config%i_gas_model_sw == IGasModelRRTMGP) then
      config%do_cloud_aerosol_per_sw_g_point = .false.
      config%n_g_sw = config%k_dist_sw%get_ngpt()

      ! Store band positions if using generalized cloud or aerosol
      ! Read from the RRTMGP k-distribution derived type
      call config%gas_optics_sw%spectral_def%allocate_bands_only( SolarReferenceTemperature, &
          & config%k_dist_sw%band_lims_wvn(1,:), config%k_dist_sw%band_lims_wvn(2,:))

      config%i_band_from_g_sw =  config%k_dist_sw%get_gpoint_bands() 
      config%n_bands_sw = config%k_dist_sw%get_nband()
      allocate(config%i_band_from_reordered_g_sw(config%n_g_sw))
      allocate(config%i_g_from_reordered_g_sw(config%n_g_sw))

      ! if (config%i_solver_sw == ISolverSpartacus) then
      !      ! SPARTACUS requires g points ordered in approximately
      !      ! increasing order of optical depth
      !      config%i_g_from_reordered_g_sw = RRTMGP_GPOINT_REORDERING_SW_224
      ! else
          ! Implied-do for no reordering
        config%i_g_from_reordered_g_sw = (/ (irep, irep=1,config%n_g_sw) /)
        ! end if

        config%i_band_from_reordered_g_sw &
              & = config%i_band_from_g_sw(config%i_g_from_reordered_g_sw)
    end if


    if (config%i_gas_model_lw == IGasModelRRTMGP_NN .or. config%i_gas_model_lw == IGasModelRRTMGP) then
      config%do_cloud_aerosol_per_lw_g_point = .false.
      config%n_g_lw = config%k_dist_lw%get_ngpt()

      call config%gas_optics_lw%spectral_def%allocate_bands_only( TerrestrialReferenceTemperature, &
      & config%k_dist_lw%band_lims_wvn(1,:), config%k_dist_lw%band_lims_wvn(2,:))

      config%i_band_from_g_lw =  config%k_dist_lw%get_gpoint_bands() 
      config%n_bands_lw = config%k_dist_lw%get_nband()
      allocate(config%i_band_from_reordered_g_lw(config%n_g_lw))
      allocate(config%i_g_from_reordered_g_lw(config%n_g_lw))

      if (config%i_solver_lw == ISolverSpartacus) then
        ! SPARTACUS requires g points ordered in approximately
        ! increasing order of optical depth
        if (config%n_g_lw == 256) then
          config%i_g_from_reordered_g_lw = RRTMGP_GPOINT_REORDERING_LW_256
        else 
          write(nulerr,'(a)') '*** Error in setup_gas_optics_ifs_rrtmgp: SPARTACUS reordering missing for this LW k-distribution'
          call radiation_abort()
        end if
      else
        ! Implied-do for no reordering
        config%i_g_from_reordered_g_lw = (/ (irep, irep=1,config%n_g_lw) /)
      end if

      config%i_band_from_reordered_g_lw &
            & = config%i_band_from_g_lw(config%i_g_from_reordered_g_lw)
    end if


    if (config%i_gas_model_sw == IGasModelRRTMGP_NN) then
        call config%rrtmgp_neural_nets(1) % load_netcdf(trim(config%rrtmgp_neural_net_sw_tau))
        call config%rrtmgp_neural_nets(2) % load_netcdf(trim(config%rrtmgp_neural_net_sw_ray))
        if (config%n_g_sw /= size(config%rrtmgp_neural_nets(1)%coeffs_output_mean)) then 
          write(nulerr,'(a,i0,a,i0)') '*** Error in setup_gas_optics_ifs_rrtmgp: NN n_g_sw of ', &
              & config%n_g_sw, ' doesnt match the LUT ng:', size(config%rrtmgp_neural_nets(1)%coeffs_output_mean)
          call radiation_abort()
        end if 
    end if 
    if (config%i_gas_model_lw == IGasModelRRTMGP_NN) then
        call config%rrtmgp_neural_nets(3) % load_netcdf(trim(config%rrtmgp_neural_net_lw))
        ! call config%rrtmgp_neural_nets(3) % load_netcdf(trim(config%rrtmgp_neural_net_lw_tau))
        ! call config%rrtmgp_neural_nets(4) % load_netcdf(trim(config%rrtmgp_neural_net_lw_pfrac))
        if (.not. config%n_g_lw == size(config%rrtmgp_neural_nets(3)%coeffs_output_mean)/2 .or. &
        config%n_g_lw == size(config%rrtmgp_neural_nets(3)%coeffs_output_mean)) then 
         write(nulerr,'(a,i0,a,i0)') '*** Error in setup_gas_optics_ifs_rrtmgp: NN ng_lw of ', &
              & config%n_g_lw, ' doesnt match the LUT ng:', size(config%rrtmgp_neural_nets(3)%coeffs_output_mean)
        call radiation_abort()
      end if 
    end if

    if (config%i_gas_model_sw == IGasModelRRTMGP_NN .or. config%i_gas_model_lw == IGasModelRRTMGP_NN ) then 
      use_neural_nets = .true.
    end if

    ! Scale the spectral solar source function with user-provided solar irradiance
    ! REMOVED - now done within the gas optics call
    ! call stop_on_err(config%k_dist_sw%set_tsi(solar_irradiance ))

    ! The i_spec_* variables are used solely for storing spectral
    ! data, and this can either be by band or by g-point
    if (config%do_save_spectral_flux) then
      if (config%do_save_gpoint_flux) then
        config%n_spec_sw = config%n_g_sw
        config%n_spec_lw = config%n_g_lw
        config%i_spec_from_reordered_g_sw => config%i_g_from_reordered_g_sw
        config%i_spec_from_reordered_g_lw => config%i_g_from_reordered_g_lw
      else
        config%n_spec_sw = config%n_bands_sw
        config%n_spec_lw = config%n_bands_lw
        config%i_spec_from_reordered_g_sw => config%i_band_from_reordered_g_sw
        config%i_spec_from_reordered_g_lw => config%i_band_from_reordered_g_lw
      end if
    else
      config%n_spec_sw = 0
      config%n_spec_lw = 0
      nullify(config%i_spec_from_reordered_g_sw)
      nullify(config%i_spec_from_reordered_g_lw)
    end if

    if (lhook) call dr_hook('radiation_ifs_rrtmgp:setup_gas_optics_ifs_rrtmgp',1,hook_handle)

  end subroutine setup_gas_optics_ifs_rrtmgp

  !---------------------------------------------------------------------
  ! Scale gas mixing ratios according to required units
  subroutine set_gas_units(gas)

  use radiation_gas,           only : gas_type, IVolumeMixingRatio
  type(gas_type),    intent(inout) :: gas

  call gas%set_units(IVolumeMixingRatio)

  end subroutine set_gas_units


  !---------------------------------------------------------------------
  ! Compute gas optical depths, shortwave scattering, Planck function
  ! and incoming shortwave radiation at top-of-atmosphere
  subroutine gas_optics_ifs_rrtmgp(ncol,nlev,istartcol,iendcol, &
        &  config, single_level, thermodynamics, gas, & 
        &  od_lw, od_sw, ssa_sw, lw_albedo, planck_hl, lw_emission, &
        &  incoming_sw)

    use parkind1,                 only : jprb, jpim
    use yomhook,                  only : lhook, dr_hook, jphook
    use radiation_config,         only : config_type, ISolverSpartacus, IGasModelRRTMGP, IGasModelRRTMGP_NN
    use radiation_thermodynamics, only : thermodynamics_type
    use radiation_single_level,   only : single_level_type
    use radiation_gas
    use radiation_io,     only : nulerr, radiation_abort     
    use mo_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp
    use mo_gas_concentrations,    only: ty_gas_concs
    use radiation_gas_constants
    use mo_rte_util_array,     only:  any_vals_outside


    integer, intent(in) :: ncol               ! number of columns
    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type), intent(in) :: config
    type(single_level_type),  intent(in) :: single_level
    type(thermodynamics_type),intent(in) :: thermodynamics
    type(gas_type),           intent(in) :: gas

    ! Longwave albedo of the surface
    real(jprb), dimension(config%n_g_lw,istartcol:iendcol), intent(in)  :: lw_albedo

    ! Gaseous layer optical depth in longwave and shortwave, and
    ! shortwave single scattering albedo (i.e. fraction of extinction
    ! due to Rayleigh scattering) at each g-point
    real(jprb), dimension(config%n_g_lw,nlev,istartcol:iendcol), intent(out) :: od_lw
    real(jprb), dimension(config%n_g_sw,nlev,istartcol:iendcol), intent(out) :: od_sw, ssa_sw

    ! The Planck function (emitted flux from a black body) at half
    ! levels at each longwave g-point
    real(jprb), dimension(config%n_g_lw,nlev+1,istartcol:iendcol), intent(out)   :: planck_hl
    ! Planck function for the surface (W m-2)
    real(jprb), dimension(config%n_g_lw,istartcol:iendcol),      intent(out)     :: lw_emission

    ! The incoming shortwave flux into a plane perpendicular to the
    ! incoming radiation at top-of-atmosphere in each of the shortwave
    ! g-points
    real(jprb), dimension(config%n_g_sw,istartcol:iendcol),  intent(out)         :: incoming_sw

    type(ty_gas_concs)    :: gas_rrtmgp

    ! RRTMGP sources are provided in W/m2-str; factor of pi converts to flux units
    real(jprb), parameter :: pi = acos(-1._jprb)

    real(jprb), dimension(:,:), allocatable :: pressure_fl, temperature_fl ! (nlev, ncol_loc)
    real(jprb), dimension(:,:), allocatable :: pressure_hl, temperature_hl ! (nlev+1, ncol_loc)

    real(jprb) :: solar_src_tmp(config%n_g_sw), norm

    integer :: ret, jcol, jlev, igpt, ncol_loc, jgas
    logical :: is_volume_mixing_ratio
    logical :: use_neural_nets =.false.

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_ifs_rrtmgp:gas_optics',0,hook_handle)
   
    if (config%i_gas_model_sw == IGasModelRRTMGP_NN .or. config%i_gas_model_lw == IGasModelRRTMGP_NN ) then 
      use_neural_nets = .true.
    end if

    ncol_loc = iendcol-istartcol+1

    allocate(pressure_hl(nlev+1, ncol_loc), temperature_hl(nlev+1, ncol_loc))
    allocate(pressure_fl(nlev, ncol_loc),   temperature_fl(nlev, ncol_loc))

    pressure_hl     = transpose(thermodynamics%pressure_hl(istartcol:iendcol,:))
    temperature_hl  = transpose(thermodynamics%temperature_hl(istartcol:iendcol,:))

    pressure_fl &
        &  = 0.5_jprb * transpose(thermodynamics%pressure_hl(istartcol:iendcol,1:nlev) &
        &               +thermodynamics%pressure_hl(istartcol:iendcol,2:nlev+1))

    temperature_fl &
        &  = 0.5_jprb * transpose(thermodynamics%temperature_hl(istartcol:iendcol,1:nlev) &
        &               +thermodynamics%temperature_hl(istartcol:iendcol,2:nlev+1))   

    ! RRTMGP doesn't work with pressures that are extremely low, need to cap it to the minimum here
    pressure_fl(1,:)    = config%k_dist_sw%get_press_min() + epsilon(config%k_dist_sw%get_press_min())

    ! temperature also needs to be in a certain range (160 - 355 K)
    temperature_hl = max(config%k_dist_lw%get_temp_min(), temperature_hl)
    temperature_hl = min(config%k_dist_lw%get_temp_max(), temperature_hl)
    temperature_fl = max(config%k_dist_lw%get_temp_min(), temperature_fl)
    temperature_fl = min(config%k_dist_lw%get_temp_max(), temperature_fl)

    ! Check that the gas concentrations are stored in volume mixing
    ! ratio with no scaling; if not, return a vector of scalings
    call gas%assert_units(IVolumeMixingRatio, scale_factor=1.0_jprb, &
         &                istatus=is_volume_mixing_ratio)
    if (.not. is_volume_mixing_ratio) then
      write(nulerr,'(a)') '*** RRTMGP requires volume mixing ratios, mixing it with another gas optics not supported'
      call radiation_abort()
    !   call gas%get_scaling(IVolumeMixingRatio, concentration_scaling)
    ! else
    !   concentration_scaling = 1.0_jprb
    end if  

    ! Initialize RRTMGP gas derived type
    call stop_on_err(gas_rrtmgp%init(config%rrtmgp_gas_names))

    ! Write nitrogen with a constant value - here we assume it's missing from ecRAD gas type 
    ! (but that 'N2' is already in config%rrtmgp_gas_names)
    ! Bit iffy to hardcode like this, but avoids polluting gas%mixing_ratio with more gases that are not used
    call stop_on_err(gas_rrtmgp%set_vmr("n2", 0.781000018_jprb))

    do jgas = 1,NMaxGases
      if (gas%is_present(jgas)) then  
        ! Write the gas, transposing because RRTMGP-NN has columns outermost
        ! Concentrations must already be in volume mixing ratios!
        if (GasLowerCaseName(jgas)=="hcfc22") then
          ! HCFC22 has different name in RRTMGP
          call stop_on_err(gas_rrtmgp%set_vmr('cfc22',transpose(gas%mixing_ratio(istartcol:iendcol,:,jgas))))
        else 
          call stop_on_err(gas_rrtmgp%set_vmr(GasLowerCaseName(jgas), &
              & transpose(gas%mixing_ratio(istartcol:iendcol,:,jgas))))
        end if
      end if
    end do

    ! do jgas = 1, size(gas_rrtmgp%concs)
    !   print *, "max of rrtmgp gas ", gas_rrtmgp%gas_name(jgas), ":", maxval(gas_rrtmgp%concs(jgas)%conc)
    ! end do

    if (config%do_sw .and. (config%i_gas_model_sw == IGasModelRRTMGP .or. &
        &   config%i_gas_model_sw == IGasModelRRTMGP_NN)) then

      if (use_neural_nets) then
        call stop_on_err(config%k_dist_sw%gas_optics_ext_ecrad( &
            &   ncol_loc, nlev, config%n_g_sw, &
            &   pressure_fl, &
            &   pressure_hl, &
            &   temperature_fl, &
            &   gas_rrtmgp, &
            &   od_sw(:,:,istartcol:iendcol), ssa_sw(:,:,istartcol:iendcol), & 
                                  ! incoming_sw, & ! <--- incoming_sw was previously an output, 
                                  ! RRTMGP filling it across (ng,ncol) from the internal solar_source(ng) 
            &   neural_nets = config%rrtmgp_neural_nets(1:2))) 
      else 
        call stop_on_err(config%k_dist_sw%gas_optics_ext_ecrad( &
            &   ncol_loc, nlev, config%n_g_sw, &
            &   pressure_fl, &
            &   pressure_hl, &
            &   temperature_fl, &
            &   gas_rrtmgp, &
            &   od_sw(:,:,istartcol:iendcol), ssa_sw(:,:,istartcol:iendcol)))!, incoming_sw)) 
      end if

      ! NEW: Scale the RRTMGP solar source function here, to avoid having to call k_dist_sw%set_solar_src()
      ! in setup_radiation, because that would require a new input in setup_radiation (single_level%solar_irradiance)
      ! set_solar_src():
      norm = 1._jprb/sum(config%k_dist_sw%solar_source(:))
      solar_src_tmp(:) = config%k_dist_sw%solar_source(:) * single_level%solar_irradiance * norm
      ! This bit was previously done inside the gas optics code:
      do jcol = istartcol,iendcol
        do igpt = 1,config%n_g_sw
          incoming_sw(igpt,jcol) = solar_src_tmp(igpt)
        end do
      end do

      od_sw(:,:,istartcol:iendcol) = max(config%min_gas_od_sw, od_sw(:,:,istartcol:iendcol))

    end if

    if (config%do_lw .and. (config%i_gas_model_lw == IGasModelRRTMGP .or. &
        &   config%i_gas_model_lw == IGasModelRRTMGP_NN)) then
        if (use_neural_nets) then
          call stop_on_err(config%k_dist_lw%gas_optics_int_ecrad( &
              &   ncol_loc, nlev, config%n_g_lw, &
              &   pressure_fl, &
              &   pressure_hl, &
              &   temperature_fl, &
              &   temperature_hl, &
              &   single_level%skin_temperature(istartcol:iendcol), &
              &   gas_rrtmgp, &
              &   od_lw(:,:,istartcol:iendcol), lw_emission(:,istartcol:iendcol), planck_hl(:,:,istartcol:iendcol), &
              &   neural_nets = config%rrtmgp_neural_nets(3:))) 
        else 
          call stop_on_err(config%k_dist_lw%gas_optics_int_ecrad( &
              &   ncol_loc, nlev, config%n_g_lw, &
              &   pressure_fl, &
              &   pressure_hl, &
              &   temperature_fl, &
              &   temperature_hl, &
              &   single_level%skin_temperature(istartcol:iendcol), &
              &   gas_rrtmgp, &
              &   od_lw(:,:,istartcol:iendcol), lw_emission(:,istartcol:iendcol), planck_hl(:,:,istartcol:iendcol)))
        end if
        ! lw_emission at this point is actually the planck function of
        ! the surface
        lw_emission = lw_emission * (1.0_jprb - lw_albedo)
        ! RRTMGP sources are provided in W/m2-str; factor of pi converts to flux units
        planck_hl = planck_hl * pi 
        ! source%lay_source = source%lay_source * pi 
        lw_emission = lw_emission * pi
  
        if (config%i_solver_lw == ISolverSpartacus) then
          !    if (.true.) then
          ! We need to rearrange the gas optics info in memory: reordering
          ! the g points in order of approximately increasing optical
          ! depth (for efficient 3D processing on only the regions of the
          ! spectrum that are optically thin for gases) and reorder in
          ! pressure since the the functions above treat pressure
          ! decreasing with increasing index.  Note that the output gas
          ! arrays have dimensions in a different order to the inputs,
          ! so there is some inefficiency here.
          do jcol = istartcol,iendcol
            lw_emission(:,jcol) = lw_emission(config%i_g_from_reordered_g_lw,jcol)
            do jlev = 1, nlev
              od_lw(:,jlev,jcol) = od_lw(config%i_g_from_reordered_g_lw,jlev,jcol)
              planck_hl(:,jlev,jcol) = planck_hl(config%i_g_from_reordered_g_lw,jlev,jcol)
            end do
            planck_hl(:,nlev+1,jcol) = planck_hl(config%i_g_from_reordered_g_lw,nlev+1,jcol)
          end do
        end if 

        od_lw(:,:,istartcol:iendcol) = max(config%min_gas_od_lw, od_lw(:,:,istartcol:iendcol))

    end if

    ! if(any_vals_outside(od_lw(:,:,istartcol:iendcol), 0.0_jprb, 10000.0_jprb)) then
    !   ! print *, " OD LW OUTSIDE", "MAX", maxval(od_lw(:,:,istartcol:iendcol)), "MIN", minval(od_lw(:,:,istartcol:iendcol))
    !   stop
    ! end if

    deallocate(pressure_fl, temperature_fl)
    deallocate(pressure_hl, temperature_hl)


   if (lhook) call dr_hook('radiation_ifs_rrtmgp:gas_optics',1,hook_handle)
    
  end subroutine gas_optics_ifs_rrtmgp


  pure function string_loc_in_array(string, array)
    character(len=*),               intent(in) :: string
    character(len=*), dimension(:), intent(in) :: array
    integer                                    :: string_loc_in_array

    integer :: i
    character(len=len_trim(string)) :: lc_string

    string_loc_in_array = -1
    lc_string = lower_case(trim(string))
    do i = 1, size(array)
      if(lc_string == lower_case(trim(array(i)))) then
        string_loc_in_array = i
        exit
      end if
    end do
  end function string_loc_in_array

  pure function lower_case( input_string ) result( output_string )
    character(len=*), intent(in) :: input_string
    character(len=len(input_string)) :: output_string
    integer :: i, n

    ! Copy input string
    output_string = input_string

    ! Convert case character by character
    do i = 1, len(output_string)
      n = index(UPPER_CASE_CHARS, output_string(i:i))
      if ( n /= 0 ) output_string(i:i) = LOWER_CASE_CHARS(n:n)
    end do
  end function


  subroutine stop_on_err(error_msg)
    use iso_c_binding
    use radiation_io,   only : nulerr, radiation_abort

    character(len=*), intent(in) :: error_msg
  
    if(error_msg /= "") then
        write(nulerr,'(a,a)') &
           &  '*** radiation_ifs_rrtmgp error: ', trim(error_msg)
      call radiation_abort()
    end if
  end subroutine stop_on_err

  
  end module radiation_ifs_rrtmgp
  
