  ! radiation_spartacus_sw.F90 - SPARTACUS shortwave solver
!
! (C) Copyright 2014- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
!
! Modifications
!   2017-04-11  R. Hogan  Receive albedos at g-points
!   2017-04-22  R. Hogan  Store surface fluxes at all g-points
!   2017-06-30  R. Hogan  Reformulate to use total_albedo_direct not total_source
!   2017-07-03  R. Hogan  Explicit calculation of encroachment
!   2017-10-23  R. Hogan  Renamed single-character variables
!   2018-02-20  R. Hogan  Corrected "computed" encroachment
!   2018-03-09  R. Hogan  Security on computed encroachment, transmittance and reflectance
!   2018-08-29  R. Hogan  Reformulate horizontal migration distances in step_migrations
!   2018-09-02  R. Hogan  Bug fix in x_direct in step_migrations
!   2018-09-03  R. Hogan  Security via min_cloud_effective_size
!   2018-09-04  R. Hogan  Use encroachment_scaling and read encroachment edge_length from upper layer
!   2018-09-13  R. Hogan  Added "Fractal" encroachment option
!   2018-09-14  R. Hogan  Added "Zero" encroachment option
!   2018-10-08  R. Hogan  Call calc_region_properties
!   2018-10-15  R. Hogan  Added call to fast_expm_exchange instead of expm
!   2019-01-12  R. Hogan  Use inv_inhom_effective_size if allocated
!   2019-02-10  R. Hogan  Renamed "encroachment" to "entrapment"
!   2022-09-01  P. Ukkonen  Optimizations for much better performance with ECCKD, including: 
!                           batching section 3 computations, faster kernels, ng can be defined at compile time


module radiation_spartacus_sw

  public

! Allow size of inner dimension (number of g-points) to be known at compile time if NG_SW is defined
#ifdef NG_SW
  integer, parameter, private :: ng = NG_SW
#else
#define ng ng_sw_in
#endif

contains

  ! Small routine for scaling cloud optical depth in the cloudy
  ! regions
#include "radiation_optical_depth_scaling.h"

  ! This module contains just one exported subroutine, the shortwave
  ! solver using the Speedy Algorithm for Radiative Transfer through
  ! Cloud Sides (SPARTACUS), which can represent 3D effects using a
  ! matrix form of the two-stream equations.
  !
  ! Sections:
  !   1: Prepare general variables and arrays
  !   2: Prepare column-specific variables and arrays
  !   3: First loop over layers
  !     3.1: Layer-specific initialization
  !     3.2: Compute gamma variables
  !       3.2a: Clear-sky case
  !       3.2b: Cloudy case
  !     3.3: Compute reflection, transmission and sources
  !       3.3a: g-points with 3D effects
  !       3.3b: g-points without 3D effects
  !   4: Compute total albedos
  !     4.1: Adding method
  !     4.2: Overlap and entrapment
  !   5: Compute fluxes
  subroutine solver_spartacus_sw(ng_sw_in, nlev,istartcol,iendcol, &
       &  config, single_level, thermodynamics, cloud, & 
       &  od, ssa, g, od_cloud, ssa_cloud, g_cloud, &
       &  albedo_direct, albedo_diffuse, incoming_sw, &
       &  flux)

    use parkind1, only           : jprb
    use yomhook,  only           : lhook, dr_hook, jphook
    use radiation_io, only             : nulout
    use radiation_config, only         : config_type, IPdfShapeGamma, &
         &  IEntrapmentZero, IEntrapmentEdgeOnly, IEntrapmentExplicit, &
         &  IEntrapmentExplicitNonFractal, IEntrapmentMaximum
    use radiation_single_level, only   : single_level_type
    use radiation_thermodynamics, only : thermodynamics_type
    use radiation_cloud, only          : cloud_type
    use radiation_regions, only        : calc_region_properties
    use radiation_overlap, only        : calc_overlap_matrices_nocol, calc_overlap_matrices_nocol_dp
    use radiation_flux, only           : flux_type, &
         &                               indexed_sum, add_indexed_sum
    use radiation_matrix
    use radiation_two_stream, only     : calc_two_stream_gammas_sw, &
         & calc_ref_trans_sw, calc_frac_scattered_diffuse_sw, &
         &  SwDiffusivity
    use radiation_constants, only      : Pi, GasConstantDryAir, AccelDueToGravity

    implicit none

    ! Inputs
    integer, intent(in) :: ng_sw_in           ! number of g-points
    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type),        intent(in) :: config
    type(single_level_type),  intent(in) :: single_level
    type(thermodynamics_type),intent(in) :: thermodynamics
    type(cloud_type),         intent(in) :: cloud

    ! Gas and aerosol optical depth, single-scattering albedo and
    ! asymmetry factor at each shortwave g-point
    real(jprb), intent(in), dimension(ng,nlev,istartcol:iendcol) :: &
         &  od, ssa, g

    ! Cloud and precipitation optical depth, single-scattering albedo and
    ! asymmetry factor in each shortwave band
    real(jprb), intent(in), dimension(config%n_bands_sw,nlev,istartcol:iendcol)   :: &
         &  od_cloud, ssa_cloud, g_cloud

    ! Direct and diffuse surface albedos, and the incoming shortwave
    ! flux into a plane perpendicular to the incoming radiation at
    ! top-of-atmosphere in each of the shortwave g points
    real(jprb), intent(in), dimension(ng,istartcol:iendcol) :: &
         &  albedo_direct, albedo_diffuse, incoming_sw

    ! Output
    type(flux_type), intent(inout):: flux
    ! Local variables
    ! integer :: nreg !, ng
    integer, parameter :: nreg = 3
    integer :: nregactive ! =1 in clear layer, =nreg in a cloudy layer
    integer :: jcol, jlev, jg, jreg, iband, jreg2,jreg3, j1, j2
#ifdef EXPLICIT_EDGE_ENTRAPMENT
    integer :: jreg4
#endif
    integer :: ng3D ! Number of g-points with small enough gas optical
                    ! depth that 3D effects need to be represented

    ! Ratio of gas constant for dry air to acceleration due to gravity
    real(jprb), parameter :: R_over_g = GasConstantDryAir / AccelDueToGravity

    real(jprb) :: mu0, one_over_mu0, tan_sza, transfer_scaling

    ! The tangent of the effective zenith angle for diffuse radiation
    ! that is appropriate for 3D transport
    real(jprb), parameter :: tan_diffuse_angle_3d = Pi * 0.5_jprb

    ! The minimum cosine of solar zenith angle to consider for 3D
    ! effects, equivalent to one solar radius (0.2615 degrees)
    real(jprb), parameter :: min_mu0_3d = 0.004625_jprb

    ! Optical depth, single scattering albedo and asymmetry factor in
    ! each region and (except for asymmetry) at each g-point
!     real(jprb), dimension(ng, nreg) &
!          &  :: od_region, ssa_region
    real(jprb), dimension(ng, nlev) :: od_region_clear  

    ! Scattering optical depths of gases and clouds
    real(jprb) :: scat_od, scat_od_cloud

    ! The area fractions of each region
    real(jprb) :: region_fracs(1:nreg,nlev,istartcol:iendcol)

    ! The scaling used for the optical depth in the cloudy regions
    real(jprb) :: od_scaling(2:nreg,nlev,istartcol:iendcol)

    ! The length of the interface between regions jreg and jreg+1 per
    ! unit area of gridbox, equal to L_diff^ab in Hogan and Shonk
    ! (2013). When the first index is 3, this is the length of the
    ! interface between regions 3 and 1. This is actually the
    ! effective length oriented to a photon with random azimuth angle,
    ! so is the true edge length divided by pi.
    real(jprb) :: edge_length(3,nlev)

    ! Element i,j gives the rate of 3D transfer of diffuse/direct
    ! radiation from region i to region j, multiplied by the thickness
    ! of the layer in m
    real(jprb) :: transfer_rate_diffuse(nreg,nreg)
    ! real(jprb) :: transfer_rate_direct(nreg,nreg)

    ! Directional overlap matrices defined at all layer interfaces
    ! including top-of-atmosphere and the surface
    real(jprb), dimension(nreg,nreg,nlev+1) :: u_matrix, v_matrix

    ! Two-stream variables

    real(jprb), dimension(ng, nlev) &
          &  :: gamma1_clear, gamma2_clear, gamma3_clear
          
    ! Matrix Gamma multiplied by the layer thickness z1, so units
    ! metres.  After calling expm, this array contains the matrix
    ! exponential of the original.
    real(jprb), dimension(:,:,:,:), allocatable, target  &
          &    ::   Gamma_z1

    ! Diffuse reflection and transmission matrices of each layer
    real(jprb), dimension(ng, nreg, &
         &  nreg, nlev) :: reflectance, transmittance

    ! Clear-sky diffuse reflection and transmission matrices of each
    ! layer
    real(jprb), dimension(ng, nlev) :: ref_clear, trans_clear

    ! Matrices translating the direct flux entering the layer from
    ! above to the reflected radiation exiting upwards (ref_dir) and
    ! the scattered radiation exiting downwards (trans_dir_diff),
    ! along with the direct unscattered transmission matrix
    ! (trans_dir_dir).
    real(jprb), dimension(ng, nreg, nreg, nlev) &
         &  :: ref_dir, trans_dir_diff, trans_dir_dir
    ! ...clear-sky equivalents
    real(jprb), dimension(ng, nlev) &
         &  :: ref_dir_clear, trans_dir_diff_clear, trans_dir_dir_clear

    ! The fluxes downwelling from the bottom of the layer due to
    ! scattering by the direct beam within the layer
    real(jprb), dimension(ng, nreg) :: source_dn
    ! ...clear-sky equivalent
    real(jprb), dimension(ng) :: source_dn_clear

    ! The fluxes upwelling just above the base of a layer due to
    ! reflection of the direct downwelling beam; this is just used as
    ! a temporary variable
    real(jprb), dimension(ng, nreg) :: total_source

    ! Direct downwelling flux below and above an interface between
    ! layers into a plane perpendicular to the direction of the sun
    real(jprb), dimension(ng, nreg) &
         &  :: direct_dn_below, direct_dn_above
    ! ...clear-sky equivalent (no distinction between "above/below")
    real(jprb), dimension(ng) :: direct_dn_clear

    ! Total albedo of the atmosphere/surface just above a layer
    ! interface with respect to downwelling diffuse and direct
    ! radiation at that interface, where level index = 1 corresponds
    ! to the top-of-atmosphere
    real(jprb), dimension(ng, nreg, &
         &  nreg, nlev+1) :: total_albedo, total_albedo_direct
    ! ...clear-sky equivalent
    real(jprb), dimension(ng, nlev+1) &
         &  :: total_albedo_clear, total_albedo_clear_direct

    ! As total_albedo, but just below a layer interface
    real(jprb), dimension(ng, nreg, nreg) &
         &  :: total_albedo_below, total_albedo_below_direct

    ! Temporary array for applying adding method with entrapment to
    ! albedo matrices
    real(jprb), dimension(ng, nreg, nreg) &
         &  :: tmp, tmp2

    ! Horizontal migration distance (m) of reflected light
    real(jprb), dimension(ng, nreg) &
         &  :: x_diffuse, x_direct
    ! Temporary variables when applying overlap rules
    real(jprb), dimension(ng, nreg) &
         &  :: x_diffuse_above, x_direct_above
         
#define USE_FAST_EXPM_EXCHANGE 1
! #ifdef USE_FAST_EXPM_EXCHANGE
!     real(jprb), dimension(ng, nreg-1, nreg-1) :: entrapment
! #else
    real(jprb), dimension(ng, 2, nreg, nreg, nreg) :: entrapment
! #endif
    ! The following is used to store matrices of the form I-A*B that
    ! are used on the denominator of some expressions
    real(jprb) :: denominator(ng,nreg,nreg)

    ! Clear-sky equivalent, but actually its reciprocal to replace
    ! some divisions by multiplications
    real(jprb) :: inv_denom_scalar

    ! Final step in working out how much transport between regions
    ! above occurs
    real(jprb), dimension(ng) :: fractal_fac_dir
    real(jprb), dimension(ng) :: fractal_fac_dif

    ! Inverse of cloud effective size (m^-1)
    real(jprb) :: inv_effective_size

    ! Layer depth (m)
    real(jprb) :: dz, layer_depth(nlev)

    ! Upwelling and downwelling fluxes above/below layer interfaces
    real(jprb), dimension(ng, nreg) &
         &  :: flux_up_above, flux_dn_above, flux_dn_below
    ! Clear-sky upwelling and downwelling fluxes (which don't
    ! distinguish between whether they are above/below a layer
    ! interface)
    real(jprb), dimension(ng) :: flux_up_clear, flux_dn_clear

    ! Index of top-most cloudy layer, or nlev+1 if no cloud
    integer :: i_cloud_top

    ! Keep a count of the number of calls to the two ways of computing
    ! reflectance/transmittance matrices
    integer :: n_calls_expm, n_calls_meador_weaver

    ! Identify clear-sky layers, with pseudo layers for outer space
    ! and below the ground, both treated as single-region clear skies
    logical :: is_clear_sky_layer(0:nlev+1)

    ! Used in computing rates of lateral radiation transfer
    real(jprb), parameter :: four_over_pi = 4.0_jprb / Pi

    ! Maximum entrapment coefficient
    real(jprb) :: max_entr

    ! Faster broadband flux computation 
    real(jprb) :: sums_up, sums_dn, sums_dn_dir, sums_up_clear, sums_dn_clear, sums_dn_dir_clear

    real(jphook) :: hook_handle

    !TEMP
    integer :: bad_col

#ifdef PARKIND1_SINGLE
    real(jprb), parameter :: coeff = 1000.0_jprb * 100.0_jprb * epsilon(1.0_jprb)
#else 
    real(jprb), parameter :: coeff = 100.0_jprb * epsilon(1.0_jprb)
#endif
    ! New for optimized SPARTACUS
    ! VARIABLES FOR BATCHED COMPUTATIONS
    logical :: is_cloudy_layer(1:nlev), are_clouds_below
    real(jprb), dimension(:,:,:), allocatable                        &
         &    ::   transfer_rate_dir, transfer_rate_dif,             &
         &         od_region_cld, ssa_region_cld, g_region_cld,      &  
         &         gamma1_cld, gamma2_cld, gamma3_cld    
    integer :: jtop, jbot, nlev_cld, nlev_cld_limit

    if (lhook) call dr_hook('radiation_spartacus_sw:solver_spartacus_sw',0,hook_handle)

    ! --------------------------------------------------------
    ! Section 1: Prepare general variables and arrays
    ! --------------------------------------------------------

    ! Copy array dimensions to local variables for convenience
    ! nreg = config%nregions
    ! ng   = ng_sw

    ! BLOCKING OF SECTION 3
    ! For better performance the expm computations in section 3 are batched across several 
    ! cloudy layers. It may be useful to cap how many adjacent cloudy levels are batched 
    ! (too large and cache use suffers)
    ! Since the size of Gamma is ng*nlevs*9*9, and ng changes with the gas optics model (32 for ECCKD),
    ! some expression which depends on ng and working precision (and ideally cache size) should be used
    ! nlev_cld_limit = 1536 / (ng* sizeof(1.0_jprb))  != 12 if ng=32 and wp=sp (sizeof=4)
    nlev_cld_limit = max(nint(800.0_jprb / (ng* sizeof(1.0_jprb))),1)  != 2 for RRTMG, 6 for ECCKD
    ! for optimal performance user can hand-tune nlev_cld_limit

    ! Reset count of number of calls to the two ways to compute
    ! reflectance/transmission matrices
    n_calls_expm          = 0
    n_calls_meador_weaver = 0

    ! Compute the wavelength-independent region fractions and
    ! optical-depth scalings
    call calc_region_properties(nlev,nreg,istartcol,iendcol, &
         &  config%i_cloud_pdf_shape == IPdfShapeGamma, &
         &  cloud%fraction, cloud%fractional_std, region_fracs, &
         &  od_scaling, config%cloud_fraction_threshold)

    if (config%iverbose >= 3) then
         write(nulout,'(a)',advance='no') '  Processing columns'
    end if

    ! Main loop over columns
    do jcol = istartcol, iendcol

      ! --------------------------------------------------------
      ! Section 2: Prepare column-specific variables and arrays
      ! --------------------------------------------------------
 
      if (config%iverbose >= 3) then
      write(nulout,'(a)',advance='no') '.'
      end if
 
      ! Copy local cosine of the solar zenith angle
      mu0 = single_level%cos_sza(jcol)
 
      ! Skip profile if sun is too low in the sky
      if (mu0 < 1.0e-10_jprb) then
        flux%sw_dn(jcol,:) = 0.0_jprb
        flux%sw_up(jcol,:) = 0.0_jprb
        if (allocated(flux%sw_dn_direct)) then
          flux%sw_dn_direct(jcol,:) = 0.0_jprb
        end if
        if (config%do_clear) then
          flux%sw_dn_clear(jcol,:) = 0.0_jprb
          flux%sw_up_clear(jcol,:) = 0.0_jprb
          if (allocated(flux%sw_dn_direct_clear)) then
            flux%sw_dn_direct_clear(jcol,:) = 0.0_jprb
          end if
        end if

        if (config%do_save_spectral_flux) then
          flux%sw_dn_band(:,jcol,:) = 0.0_jprb
          flux%sw_up_band(:,jcol,:) = 0.0_jprb
          if (allocated(flux%sw_dn_direct_band)) then
            flux%sw_dn_direct_band(:,jcol,:) = 0.0_jprb
          end if
          if (config%do_clear) then
            flux%sw_dn_clear_band(:,jcol,:) = 0.0_jprb
            flux%sw_up_clear_band(:,jcol,:) = 0.0_jprb
            if (allocated(flux%sw_dn_direct_clear_band)) then
              flux%sw_dn_direct_clear_band(:,jcol,:) = 0.0_jprb
            end if
          end if
        end if

        flux%sw_dn_diffuse_surf_g(:,jcol) = 0.0_jprb
        flux%sw_dn_direct_surf_g(:,jcol)  = 0.0_jprb
        if (config%do_clear) then
          flux%sw_dn_diffuse_surf_clear_g(:,jcol) = 0.0_jprb
          flux%sw_dn_direct_surf_clear_g(:,jcol)  = 0.0_jprb
        end if

        cycle
      end if ! sun is below the horizon

      ! At this point mu0 >= 1.0e-10

      ! Used to compute rate of attenuation of direct solar beam
      ! through the atmosphere; for stability we limit this to one
      ! solar radius above the horizon.  Elsewhere, mu0 is used to
      ! scale the incoming fluxes so is not changed.
      one_over_mu0 = 1.0_jprb / max(min_mu0_3d, mu0)

      ! The rate at which direct radiation enters cloud sides is
      ! proportional to the tangent of the solar zenith angle
      if (one_over_mu0 > 1.0_jprb) then
        tan_sza = sqrt(one_over_mu0*one_over_mu0 - 1.0_jprb &
             &         + config%overhead_sun_factor)
      else
        ! Just in case we get mu0 > 1...
        tan_sza = sqrt(config%overhead_sun_factor)
      end if

      ! Define which layers contain cloud; assume that
      ! cloud%crop_cloud_fraction has already been called
      is_clear_sky_layer = .true.
      is_cloudy_layer = .false.
      i_cloud_top = nlev+1
      do jlev = nlev,1,-1
        if (cloud%fraction(jcol,jlev) > 0.0_jprb) then
          is_clear_sky_layer(jlev) = .false.
          i_cloud_top = jlev
        end if


        layer_depth(jlev) = R_over_g &
             &  * (thermodynamics%pressure_hl(jcol,jlev+1) &
             &     - thermodynamics%pressure_hl(jcol,jlev)) &
             &  * (thermodynamics%temperature_hl(jcol,jlev) &
             &     + thermodynamics%temperature_hl(jcol,jlev+1)) &
             &  / (thermodynamics%pressure_hl(jcol,jlev) &
             &     + thermodynamics%pressure_hl(jcol,jlev+1))

        if (config%do_3d_effects .and. &
            &  allocated(cloud%inv_cloud_effective_size) .and. &
            &  .not. (nreg == 2 .and. cloud%fraction(jcol,jlev) &
            &  > 1.0-config%cloud_fraction_threshold)) then
          if (cloud%inv_cloud_effective_size(jcol,jlev) > 0.0_jprb) then

            if (.not. is_clear_sky_layer(jlev)) is_cloudy_layer(jlev) = .true. 
            ! Cloudy in 3D sense, fulfilling above criteria

            ! Compute cloud edge length per unit area of gridbox
            ! from rearranging Hogan & Shonk (2013) Eq. 45, but
            ! adding a factor of (1-frac) so that a region that
            ! fully occupies the gridbox (frac=1) has an edge of
            ! zero. We can therefore use the fraction of clear sky,
            ! region_fracs(1,jlev,jcol) for convenience instead. The
            ! pi on the denominator means that this is actually edge
            ! length with respect to a light ray with a random
            ! azimuthal direction.
            edge_length(1,jlev) = four_over_pi &
                &  * region_fracs(1,jlev,jcol)*(1.0_jprb-region_fracs(1,jlev,jcol)) &
                &  * min(cloud%inv_cloud_effective_size(jcol,jlev), &
                &        1.0_jprb / config%min_cloud_effective_size)
            if (nreg > 2) then
              ! The corresponding edge length between the two cloudy
              ! regions is computed in the same way but treating the
              ! optically denser of the two regions (region 3) as
              ! the cloud; note that the fraction of this region may
              ! be less than that of the optically less dense of the
              ! two regions (region 2).  For increased flexibility,
              ! the user may specify the effective size of
              ! inhomogeneities separately from the cloud effective
              ! size.
              if (allocated(cloud%inv_inhom_effective_size)) then
                edge_length(2,jlev) = four_over_pi &
                    &  * region_fracs(3,jlev,jcol)*(1.0_jprb-region_fracs(3,jlev,jcol)) &
                    &  * min(cloud%inv_inhom_effective_size(jcol,jlev), &
                    &        1.0_jprb / config%min_cloud_effective_size)
              else
                edge_length(2,jlev) = four_over_pi &
                    &  * region_fracs(3,jlev,jcol)*(1.0_jprb-region_fracs(3,jlev,jcol)) &
                    &  * min(cloud%inv_cloud_effective_size(jcol,jlev), &
                    &        1.0_jprb / config%min_cloud_effective_size)
              end if

              ! In the case of three regions, some of the cloud
              ! boundary may go directly to the "thick" region
              if (config%clear_to_thick_fraction > 0.0_jprb) then
                edge_length(3,jlev) = config%clear_to_thick_fraction &
                    &  * min(edge_length(1,jlev), edge_length(2,jlev))
                edge_length(1,jlev) = edge_length(1,jlev) - edge_length(3,jlev)
                edge_length(2,jlev) = edge_length(2,jlev) - edge_length(3,jlev)
              else 
                edge_length(3,jlev) = 0.0_jprb
              end if
            end if
          end if

        end if
        
      end do ! jlev
 
      ! call calc_overlap_matrices_nocol(nlev,nreg,is_clear_sky_layer, &
      ! &  region_fracs(:,:,jcol), cloud%overlap_param(jcol,:), &
      ! &  v_matrix, u_matrix=u_matrix, decorrelation_scaling=config%cloud_inhom_decorr_scaling, &
      ! &  cloud_fraction_threshold=config%cloud_fraction_threshold, &
      ! &  use_beta_overlap=config%use_beta_overlap, &
      ! &  cloud_cover=flux%cloud_cover_sw(jcol))
      call calc_overlap_matrices_nocol_dp(nlev, nreg, is_clear_sky_layer, &
      &  region_fracs(:,:,jcol), cloud%overlap_param(jcol,:), &
      &  v_matrix, u_matrix=u_matrix, decorrelation_scaling=config%cloud_inhom_decorr_scaling, &
      &  cloud_fraction_threshold=config%cloud_fraction_threshold, &
      &  use_beta_overlap=config%use_beta_overlap, &
      &  cloud_cover=flux%cloud_cover_sw(jcol))

      ! Horizontal migration distances of reflected radiation at the
      ! surface are zero
      x_diffuse = 0.0_jprb
      x_direct  = 0.0_jprb
 
      ! --------------------------------------------------------
      ! Section 3: First loop over layers
      ! --------------------------------------------------------
      ! In this section the reflectance, transmittance and sources
      ! are computed for each layer
 
      ! ------------------------ TEST BATCHING FOR SECTION 3 ------------------------------------------------------
      ! Improve performance by batching together first all clear-sky layers, and then adjacent cloudy layers
 
      ! As a result the vectorized dimension in expm becomes larger, but if it becomes too large performance
      ! will suffer when the array doesn't fit in fast cache
      ! To make things simpler, 3D computations are done for all g-points in a cloudy layer, 
      ! which in the shortwave does not waste many computations because ng3D ~ ng
      ! To avoid overflow in expm the optical depth for clear-sky region must then be capped 

      ! 1. Reftrans computations for clear-sky-region : these are done for all layers, also for cloudy layers
      od_region_clear = od(:,:,jcol)

      call calc_ref_trans_sw(ng,nlev, &
          &  mu0, od_region_clear, ssa(:,:,jcol), g(:,:,jcol), &
          &  ref_clear, trans_clear, ref_dir_clear, trans_dir_diff_clear, trans_dir_dir_clear, &
          &  gamma1_clear, gamma2_clear, gamma3_clear)

      ! for 3D computations, cap the optical depth of clear-sky region to a threshold
      od_region_clear = min(od_region_clear,config%max_gas_od_3D)
 
      ! Initialization 
      trans_dir_dir = 0.0_jprb
      reflectance   = 0.0_jprb
      transmittance = 0.0_jprb
      ref_dir       = 0.0_jprb
      trans_dir_diff= 0.0_jprb
 
      trans_dir_dir (:,1,1,:) = trans_dir_dir_clear 
      reflectance   (:,1,1,:) = ref_clear
      transmittance (:,1,1,:) = trans_clear
      ref_dir       (:,1,1,:) = ref_dir_clear
      trans_dir_diff(:,1,1,:) = trans_dir_diff_clear
 
      ! 2.  Cloudy-sky computations
 
      are_clouds_below = .false.
 
      !jtop = findloc(is_cloudy_layer(1:nlev), .true.,dim=1) ! returns index of first cloudy layer
      ! findloc not working for some compilers. manual implementation:
      jtop = 0
      do jlev = 1, nlev
        if (is_cloudy_layer(jlev)) then
          jtop = jlev
          exit
        end if
      end do
 
      if (jtop>0) are_clouds_below = .true.
 
      do while (are_clouds_below)
        ! Find the bottom of this cloudy layer
        ! we could already be at the lowest level, then jtop=jbot=nlev and there's only this one cloudy layer to compute
        if (jtop == nlev) then
          jbot = nlev
        else
          ! otherwise, find the bottom by starting from the top and break if there's a clear-sky layer:
          ! jbot is above this layer, or the lowest level if clouds reach the surface
          do jlev = jtop+1,nlev
            if (.not. is_cloudy_layer(jlev)) then
              jbot = jlev-1 ! bottom found, exit loop
              exit
            else if (jlev==nlev) then
              jbot = nlev
            end if
          end do
        end if

        nlev_cld = jbot - jtop + 1
          
        ! Limit the number of layers here so that Gamma1 is a suitable
        ! size for better cache performance
        ! So instead of Gamma_z1(ng,jtop:jbot,        9,9),
        !               Gamma_z1(ng,jtop:jbot_block,  9,9)
        ! where jbot_block is set to something between jtop+1 and jbot
        if (nlev_cld > nlev_cld_limit) then
          jbot = jtop + nlev_cld_limit - 1
          nlev_cld = jbot - jtop + 1
        end if

        ! Allocations
        allocate(transfer_rate_dir(jtop:jbot, nreg, nreg), transfer_rate_dif(jtop:jbot, nreg, nreg))
        allocate(od_region_cld(ng,jtop:jbot,2:nreg), ssa_region_cld(ng,jtop:jbot,2:nreg))
        allocate(g_region_cld (ng,jtop:jbot,2:nreg))
        allocate(gamma1_cld(ng,jtop:jbot,2:nreg), gamma2_cld(ng,jtop:jbot,2:nreg), &
            & gamma3_cld(ng,jtop:jbot,2:nreg))
        allocate(Gamma_z1(ng, jtop:jbot,3*nreg,3*nreg))
        
        ! Array-wise assignments
        transfer_rate_dir   = 0.0_jprb
        transfer_rate_dif   = 0.0_jprb
        Gamma_z1           = 0.0_jprb

        ! --------------------------------------------------------------
        ! ------------------ START CLOUDY COMPUTATIONS -----------------
        ! --------------------------------------------------------------

        ! --------------------------------------------------------
        ! Section 3.2: Compute transfer rates and gamma variables
        ! --------------------------------------------------------
        do jreg = 1, nreg-1
          do jlev = jtop, jbot
            ! Depth of current layer
            dz = layer_depth(jlev)

            ! Compute lateral transfer rates from region jreg to
            ! jreg+1 following Hogan & Shonk (2013) Eq. 47, but
            ! multiplied by dz because the transfer rate is
            ! vertically integrated across the depth of the layer
            if (region_fracs(jreg,jlev,jcol) > epsilon(1.0_jprb)) then
              transfer_rate_dir(jlev,jreg,jreg+1) = dz &
              & * edge_length(jreg,jlev) * tan_sza / region_fracs(jreg,jlev,jcol)
              transfer_rate_dif(jlev,jreg,jreg+1) = dz &
              & * edge_length(jreg,jlev) * tan_diffuse_angle_3d / region_fracs(jreg,jlev,jcol)
            end if
            ! Compute transfer rates from region jreg+1 to
            ! jreg
            if (region_fracs(jreg+1,jlev,jcol) > epsilon(1.0_jprb)) then
              transfer_rate_dir(jlev,jreg+1,jreg) = dz &
              & * edge_length(jreg,jlev) * tan_sza / region_fracs(jreg+1,jlev,jcol)
              transfer_rate_dif(jlev,jreg+1,jreg) = dz &
              & * edge_length(jreg,jlev) * tan_diffuse_angle_3d / region_fracs(jreg+1,jlev,jcol)
            end if
          end do
        end do

        do jlev = jtop, jbot
          dz = layer_depth(jlev)
          ! Compute transfer rates directly between regions 1 and 3
          if (edge_length(3,jlev) > 0.0_jprb) then
            if (region_fracs(1,jlev,jcol) > epsilon(1.0_jprb)) then
              transfer_rate_dir(jlev,1,3) = dz &
              &  * edge_length(3,jlev) * tan_sza / region_fracs(1,jlev,jcol)
              transfer_rate_dif(jlev,1,3) = dz &
              &  * edge_length(3,jlev) * tan_diffuse_angle_3d / region_fracs(1,jlev,jcol)
            end if
            if (region_fracs(3,jlev,jcol) > epsilon(1.0_jprb)) then
              transfer_rate_dir(jlev,3,1) = dz &
              &  * edge_length(3,jlev) * tan_sza / region_fracs(3,jlev,jcol)
              transfer_rate_dif(jlev,3,1) = dz &
              &  * edge_length(3,jlev) * tan_diffuse_angle_3d / region_fracs(3,jlev,jcol)
            end if
          end if
        end do

        ! Don't allow the transfer rate out of a region to be
        ! equivalent to a loss of exp(-10) through the layer
        transfer_rate_dir = min(transfer_rate_dir, config%max_3d_transfer_rate)
        transfer_rate_dif = min(transfer_rate_dif, config%max_3d_transfer_rate)

        ! Compute scattering properties of the regions at each
        ! g-point, mapping from the cloud properties
        ! defined in each band.
        ! Loop over cloudy regions
        do jreg = 2,nreg
          do jlev = jtop, jbot
            do jg = 1,ng
              ! Mapping from g-point to band
              if (ng == config%n_bands_sw) then ! help vectorization for ECCKD
                iband = jg
              else
                iband = config%i_band_from_reordered_g_sw(jg)
              end if
              
              ! Scattering optical depth of clear-sky region
              scat_od = od(jg,jlev,jcol)*ssa(jg,jlev,jcol)
    
              scat_od_cloud = od_cloud(iband,jlev,jcol) &
                  &  * ssa_cloud(iband,jlev,jcol)*od_scaling(jreg,jlev,jcol)
              ! Add scaled cloud optical depth to clear-sky value
              od_region_cld(jg,jlev,jreg) = od(jg,jlev,jcol) &
                  &  + od_cloud(iband,jlev,jcol)*od_scaling(jreg,jlev,jcol)
              ! Compute single-scattering albedo and asymmetry
              ! factor of gas-cloud combination
              ssa_region_cld(jg,jlev,jreg) = (scat_od+scat_od_cloud) / od_region_cld(jg,jlev,jreg)
              g_region_cld(jg,jlev,jreg)   = (scat_od*g(jg,jlev,jcol) &
                  &  + scat_od_cloud * g_cloud(iband,jlev,jcol)) &
                  &  / (scat_od + scat_od_cloud)
    
              ! Apply maximum cloud optical depth for stability in the
              ! 3D case
              od_region_cld(jg,jlev,jreg) = min(od_region_cld(jg,jlev,jreg),config%max_cloud_od_sw)
            end do   
          end do 

          ! Calculate two-stream variables gamma1-gamma3 for cloudy regions
          call calc_two_stream_gammas_sw(ng*nlev_cld, &
          &  mu0, ssa_region_cld(:,:,jreg), g_region_cld(:,:,jreg), &
          &  gamma1_cld(:,:,jreg), gamma2_cld(:,:,jreg), gamma3_cld(:,:,jreg))
        end do 

        ! --------------------------------------------------------------
        ! Section 3.3: Compute reflection, transmission and emission
        ! --------------------------------------------------------------
                            
        ! 3D effects need to be represented in (ng*nlevs) square matrices each of dimension 3*nreg by 3*nreg,
        ! computing the matrix exponential, then computing the various transmission/reflectance matrices from that

        ! Write diagonal elements of Gamma_z1
        ! Clear-sky region
        jreg = 1
        call write_gamma_diag(ng,nlev_cld, nreg, jreg, od_region_clear(:,jtop:jbot), &
        &   gamma1_clear(:,jtop:jbot), gamma2_clear(:,jtop:jbot), gamma3_clear(:,jtop:jbot), &
        &   ssa(:,jtop:jbot,jcol), one_over_mu0, Gamma_z1)
        ! Cloudy regions
        do jreg = 2, nreg
          call write_gamma_diag(ng,nlev_cld, nreg, jreg, od_region_cld(:,:,jreg), &
          &   gamma1_cld(:,:,jreg), gamma2_cld(:,:,jreg), gamma3_cld(:,:,jreg), &
          &   ssa_region_cld(:,:,jreg), one_over_mu0, Gamma_z1)
        end do

        do jreg = 1,nreg-1
          do jlev = jtop, jbot 
            do jg = 1, ng
              ! Write the elements of -Gamma1*z1 concerned with 3D
              ! transport
              Gamma_z1(jg,jlev,jreg,jreg) = Gamma_z1(jg,jlev,jreg,jreg) &
                  &  + transfer_rate_dif(jlev,jreg,jreg+1)
              Gamma_z1(jg,jlev,jreg+1,jreg+1) = Gamma_z1(jg,jlev,jreg+1,jreg+1) &
                  &  + transfer_rate_dif(jlev,jreg+1,jreg)
              Gamma_z1(jg,jlev,jreg+1,jreg) = -transfer_rate_dif(jlev,jreg,jreg+1)
              Gamma_z1(jg,jlev,jreg,jreg+1) = -transfer_rate_dif(jlev,jreg+1,jreg)
              ! Write the elements of +Gamma0*z1 concerned with 3D
              ! transport
              Gamma_z1(jg,jlev,jreg+2*nreg,jreg+2*nreg) &
                  &  = Gamma_z1(jg,jlev,jreg+2*nreg,jreg+2*nreg) &
                  &  - transfer_rate_dir(jlev,jreg,jreg+1)
              Gamma_z1(jg,jlev,jreg+2*nreg+1,jreg+2*nreg+1) &
                  &  = Gamma_z1(jg,jlev,jreg+2*nreg+1,jreg+2*nreg+1) &
                  &  - transfer_rate_dir(jlev,jreg+1,jreg)
              Gamma_z1(jg,jlev,jreg+2*nreg+1,jreg+2*nreg) = transfer_rate_dir(jlev,jreg,jreg+1)
              Gamma_z1(jg,jlev,jreg+2*nreg,jreg+2*nreg+1) = transfer_rate_dir(jlev,jreg+1,jreg)
            end do
          end do
        end do

        ! Possible flow between regions a and c
        do jlev = jtop, jbot
          if (edge_length(3,jlev) > 0.0_jprb) then
            do jg = 1, ng
              ! Diffuse transport
              Gamma_z1(jg,jlev,1,1) = Gamma_z1(jg,jlev,1,1) + transfer_rate_dif(jlev,1,3)
              Gamma_z1(jg,jlev,3,3) = Gamma_z1(jg,jlev,3,3) + transfer_rate_dif(jlev,3,1)
              Gamma_z1(jg,jlev,3,1) = -transfer_rate_dif(jlev,1,3)
              Gamma_z1(jg,jlev,1,3) = -transfer_rate_dif(jlev,3,1)
              ! Direct transport
              Gamma_z1(jg,jlev,1+2*nreg,1+2*nreg) = Gamma_z1(jg,jlev,1+2*nreg,1+2*nreg) &
                  &  - transfer_rate_dir(jlev,1,3)
              Gamma_z1(jg,jlev,3+2*nreg,3+2*nreg) = Gamma_z1(jg,jlev,3+2*nreg,3+2*nreg) &
                  &  - transfer_rate_dir(jlev,3,1)
              Gamma_z1(jg,jlev,3+2*nreg,1+2*nreg) = transfer_rate_dir(jlev,1,3)
              Gamma_z1(jg,jlev,1+2*nreg,3+2*nreg) = transfer_rate_dir(jlev,3,1)
            end do
          end if
        end do 

        ! Copy Gamma1*z1
        Gamma_z1(:,:,nreg+1:nreg+nreg,nreg+1:nreg+nreg) = -Gamma_z1(:,:,1:nreg,1:nreg)
        ! Copy Gamma2*z1
        Gamma_z1(:,:,1:nreg,nreg+1:nreg+nreg) = -Gamma_z1(:,:,nreg+1:nreg+nreg,1:nreg)

        ! Compute the matrix exponential of Gamma_z1, returning the result in-place
        ng3D = ng*nlev_cld

        ! Additional security on elements fed to matrix exponential
        ! in single precision
        if (jprb <= 4) then
          Gamma_z1 = min(Gamma_z1, 18.0_jprb)
          ! Gamma_z1 = max(-20.0_jprb, min(Gamma_z1, 20.0_jprb))
        end if

        call expm_sw(ng3D, ng, nlev_cld, Gamma_z1)

        ! Update count of expm calls
        n_calls_expm = n_calls_expm + ng3D

        ! Following computations are not able to be batched by collapsing ng and nlev because of
        ! different dimension order - trans_dir_dir(ng,nreg,nreg,nlev) and Gamma(ng,nlev,nreg,nreg)
        associate(gamma_z11=>tmp, gamma_z21 => total_albedo_below, &
              gamma_z22 => total_albedo_below_direct, gamma_z23 => total_albedo_below_direct)
        do jlev = jtop, jbot   
          ! Direct transmission matrix
          trans_dir_dir(:,:,:,jlev) = min(1.0_jprb,max(0.0_jprb, &
              &  Gamma_z1(:,jlev,2*nreg+1:3*nreg,2*nreg+1:3*nreg)))
          ! Diffuse reflectance matrix; security on negative values necessary 
          ! occasionally for very low cloud fraction and very high in-cloud optical depth
          gamma_z11 = Gamma_z1(:,jlev,1:nreg,1:nreg)
          call solve_mat_3_sw(ng, gamma_z11, &
              & Gamma_z1(:,jlev,1:nreg,nreg+1:2*nreg), reflectance(:,:,:,jlev))
          reflectance(:,:,:,jlev) = min(1.0_jprb,max(0.0_jprb,-reflectance(:,:,:,jlev)))
          ! Diffuse transmission matrix
          gamma_z21 = Gamma_z1(:,jlev,nreg+1:2*nreg,1:nreg)
          gamma_z22 = Gamma_z1(:,jlev,nreg+1:2*nreg,nreg+1:2*nreg)
          do j2 = 1,3
            do j1 = 1,3
              do jg = 1, ng
                ! T = gamma_z21 * R
                transmittance(jg,j1,j2,jlev) = gamma_z21(jg,j1,1)*reflectance(jg,1,j2,jlev) + &
                      & gamma_z21(jg,j1,2)*reflectance(jg,2,j2,jlev) + &
                      & gamma_z21(jg,j1,3)*reflectance(jg,3,j2,jlev) 
                transmittance(jg,j1,j2,jlev) = transmittance(jg,j1,j2,jlev) + gamma_z22(jg,j1,j2)
                transmittance(jg,j1,j2,jlev) = max(0.0_jprb, transmittance(jg,j1,j2,jlev))
                transmittance(jg,j1,j2,jlev) = min(1.0_jprb - reflectance(jg,j1,j2,jlev), &
                  & transmittance(jg,j1,j2,jlev))
              end do
            end do
          end do

          ! Transfer matrix between downward direct and upward diffuse
          call solve_mat_3_sw(ng, gamma_z11, &
              & Gamma_z1(:,jlev,1:nreg,2*nreg+1:3*nreg), ref_dir(:,:,:,jlev))
          ref_dir(:,:,:,jlev) = min(1.0_jprb,max(0.0_jprb,-ref_dir(:,:,:,jlev)))
          ! Transfer matrix between downward direct and downward diffuse in layer interface below.  
          ! Include correction for trans_dir_diff out of plausible bounds (note that Meador &
          ! Weaver has the same correction in radiation_two_stream.F90 - this is not just an expm thing)
          gamma_z23 = Gamma_z1(:,jlev,nreg+1:2*nreg,2*nreg+1:3*nreg)
          do j2 = 1,3
            do j1 = 1,3
              do jg = 1, ng
                trans_dir_diff(jg,j1,j2,jlev) = gamma_z21(jg,j1,1)*ref_dir(jg,1,j2,jlev) + &
                    & gamma_z21(jg,j1,2)*ref_dir(jg,2,j2,jlev) + &
                    & gamma_z21(jg,j1,3)*ref_dir(jg,3,j2,jlev)
                trans_dir_diff(jg,j1,j2,jlev) = trans_dir_diff(jg,j1,j2,jlev) + gamma_z23(jg,j1,j2)
                trans_dir_diff(jg,j1,j2,jlev) = max(0.0_jprb, trans_dir_diff(jg,j1,j2,jlev))
                trans_dir_diff(jg,j1,j2,jlev) = min(mu0, trans_dir_diff(jg,j1,j2,jlev))
              end do 
            end do
          end do
        end do    
        end associate

        ! --------------------------------------------------------------
        ! ------------------- END CLOUDY COMPUTATIONS ------------------
        ! --------------------------------------------------------------

        ! Deallocations
        deallocate(transfer_rate_dir, transfer_rate_dif)
        deallocate(od_region_cld, ssa_region_cld, g_region_cld)
        deallocate(gamma1_cld, gamma2_cld, gamma3_cld)
        deallocate(Gamma_z1)

        ! CLOUD WHILE CONSTRUCT: update clouds_below and find new jtop if true
        if (jbot==nlev) are_clouds_below=.false. ! because surface reached

        ! does another cloudy layer exist?
        if (any(is_cloudy_layer(jbot+1:nlev))) then
          ! find the cloud top 
          do jlev = jbot+1, nlev  
            if (is_cloudy_layer(jlev)) then
              jtop = jlev ! top found, exit loop
              exit
            end if
          end do
        else  ! no further cloudy regions below
          are_clouds_below=.false.
        end if
     end do

      ! --------------------------------------------------------
      ! Section 4: Compute total albedos
      ! --------------------------------------------------------

     total_albedo(:,:,:,nlev+1)        = 0.0_jprb
     total_albedo_direct(:,:,:,nlev+1) = 0.0_jprb

     if (config%do_clear) then
          total_albedo_clear(:,:)        = 0.0_jprb
          total_albedo_clear_direct(:,:) = 0.0_jprb
     end if

     ! Calculate the upwelling radiation scattered from the direct
     ! beam incident on the surface, and copy the surface albedo
     ! into total_albedo
     do jreg = 1,nreg
          total_albedo(:,jreg,jreg,nlev+1)        = albedo_diffuse(:,jcol)
          total_albedo_direct(:,jreg,jreg,nlev+1) = mu0 * albedo_direct(:,jcol)
     end do

     if (config%do_clear) then
          ! Surface albedo is the same
          total_albedo_clear(:,nlev+1) = total_albedo(1:ng,1,1,nlev+1)
          total_albedo_clear_direct(:,nlev+1) = total_albedo_direct(1:ng,1,1,nlev+1)
     end if

     ! Work back up through the atmosphere computing the total albedo
     ! of the atmosphere below that point using the adding method
     do jlev = nlev,1,-1

        ! --------------------------------------------------------
        ! Section 4.1: Adding method
        ! --------------------------------------------------------
        if (config%do_clear) then
          ! Use adding method for clear-sky arrays; note that there
          ! is no need to consider "above" and "below" quantities
          ! since with no cloud overlap to worry about, these are
          ! the same
          do jg = 1, ng
            inv_denom_scalar = 1.0_jprb &
                &  / (1.0_jprb - total_albedo_clear(jg,jlev+1)*ref_clear(jg,jlev))
            total_albedo_clear(jg,jlev) = ref_clear(jg,jlev) &
                &  + trans_clear(jg,jlev)*trans_clear(jg,jlev)*total_albedo_clear(jg,jlev+1) &
                &  * inv_denom_scalar
            total_albedo_clear_direct(jg,jlev) = ref_dir_clear(jg,jlev) &
                &  + (trans_dir_dir_clear(jg,jlev) * total_albedo_clear_direct(jg,jlev+1) &
                &    +trans_dir_diff_clear(jg,jlev) * total_albedo_clear(jg,jlev+1)) &
                &  * trans_clear(jg,jlev) * inv_denom_scalar
          end do 
        end if
 
        if (is_clear_sky_layer(jlev)) then
          ! Clear-sky layer: use scalar adding method
          total_albedo_below = 0.0_jprb
          total_albedo_below_direct = 0.0_jprb
          do jg = 1, ng
            inv_denom_scalar = 1.0_jprb &
                &  / (1.0_jprb - total_albedo(jg,1,1,jlev+1)*ref_clear(jg,jlev))
            total_albedo_below(jg,1,1) = ref_clear(jg,jlev) &
                &  + trans_clear(jg,jlev)  * trans_clear(jg,jlev) &
                &  * total_albedo(jg,1,1,jlev+1) * inv_denom_scalar
            total_albedo_below_direct(jg,1,1) = ref_dir_clear(jg,jlev) &
                &  + (trans_dir_dir_clear(jg,jlev)*total_albedo_direct(jg,1,1,jlev+1) &
                &    +trans_dir_diff_clear(jg,jlev)*total_albedo(jg,1,1,jlev+1)) &
                &  * trans_clear(jg,jlev) * inv_denom_scalar
          end do 
        else 
          ! Cloudy layer: use matrix adding method
          call identity_minus_mat_x_mat_3_sw(ng, &
                &  total_albedo(:,:,:,jlev+1), reflectance(:,:,:,jlev), denominator)

          ! Let's write the terms explicitly for clarity as well
          ! as performance (avoids array temporaries being created)
          !  B=C*D, then solve AX=B, finally Y = Z*X
          associate(B1=>total_albedo_below, B2=>total_albedo_below_direct, X=>tmp)

          call mat_x_mat_3_sw(ng,&                                ! 1. B = C*D
              & total_albedo(:,:,:,jlev+1), transmittance(:,:,:,jlev), B1)
          call solve_mat_3_sw(ng,denominator, B1, X)              ! 2. solve AX=B
          call mat_x_mat_3_sw(ng, transmittance(:,:,:,jlev), X, & ! 3. Y = Z*X ...
              & total_albedo_below)
          total_albedo_below = total_albedo_below + reflectance(:,:,:,jlev) ! + R

          ! Direct
          call mat_x_mat_3_sw(ng, &                               ! 1. B = C1*D1 ...
              &  total_albedo_direct(:,:,:,jlev+1), trans_dir_dir(:,:,:,jlev), X)
          call mat_x_mat_3_sw(ng, &                               ! + C2+D3
              & total_albedo(:,:,:,jlev+1), trans_dir_diff(:,:,:,jlev), B2)
          B2 = X + B2
          call solve_mat_3_sw(ng,denominator,B2, X)               ! 2. solve AX=B
          call mat_x_mat_3_sw(ng,transmittance(:,:,:,jlev),X,B2)  ! Z*X
          total_albedo_below_direct = ref_dir(:,:,:,jlev) + B2    ! 3. Y = R + Z*X

          end associate
        end if
        ! --------------------------------------------------------
        ! Section 4.2: Overlap and entrapment
        ! --------------------------------------------------------
        
#ifndef PRINT_ENTRAPMENT_DATA
        if ((config%i_3d_sw_entrapment == IEntrapmentExplicitNonFractal &
             &  .or. config%i_3d_sw_entrapment == IEntrapmentExplicit) &
             &  .and. jlev >= i_cloud_top) then
#else
        if (config%i_3d_sw_entrapment == IEntrapmentExplicitNonFractal &
             &  .or. config%i_3d_sw_entrapment == IEntrapmentExplicit) then
#endif
          !  "Explicit entrapment": we have the horizontal migration
          !  distances just above the base of the layer, and need to
          !  step them to just below the top of the same layer
          call step_migrations(ng, nreg, cloud%fraction(jcol,jlev), &
               & layer_depth(jlev), tan_diffuse_angle_3d, tan_sza, &
               &  reflectance(:,:,:,jlev), transmittance(:,:,:,jlev), &
               &  ref_dir(:,:,:,jlev), trans_dir_dir(:,:,:,jlev), &
               &  trans_dir_diff(:,:,:,jlev), total_albedo(:,:,:,jlev+1), &
               &  total_albedo_direct(:,:,:,jlev+1), &
               &  x_diffuse, x_direct)

#ifdef PRINT_ENTRAPMENT_DATA
          ! Write out for later analysis: these are the entrapment
          ! statistics at the top of layer "jlev"
          ! Note that number of scattering events is now not computed,
          ! so print "1.0"
          write(101,'(i4,i4,9e14.6)') jcol, jlev, &
                &  x_direct(1,1:3), x_diffuse(1,1:3), 1.0_jprb,1.0_jprb,1.0_jprb
#endif

        end if

        ! Account for cloud overlap when converting albedo and source
        ! below a layer interface to the equivalent values just above
        if (is_clear_sky_layer(jlev) .and. is_clear_sky_layer(jlev-1)) then
          ! If both layers are cloud free, this is trivial...
          total_albedo(:,:,:,jlev) = 0.0_jprb
          total_albedo(:,1,1,jlev) = total_albedo_below(:,1,1)
          total_albedo_direct(:,:,:,jlev) = 0.0_jprb
          total_albedo_direct(:,1,1,jlev) = total_albedo_below_direct(:,1,1)

        else if (config%i_3d_sw_entrapment == IEntrapmentMaximum &
             &  .or. is_clear_sky_layer(jlev-1)) then
          ! "Maximum entrapment": use the overlap matrices u_matrix and v_matrix
          ! (this is the original SPARTACUS method)
          tmp = mat_x_singlemat_3_sw(ng,&
              & total_albedo_below,v_matrix(:,:,jlev))
          total_albedo(:,:,:,jlev) = singlemat_x_mat_3_sw(ng,&
              &  u_matrix(:,:,jlev), tmp)
          tmp = mat_x_singlemat_3_sw(ng,&
              & total_albedo_below_direct, v_matrix(:,:,jlev))
          total_albedo_direct(:,:,:,jlev) = singlemat_x_mat_3_sw(ng,&
              &  u_matrix(:,:,jlev), tmp)
        else if (config%i_3d_sw_entrapment == IEntrapmentZero) then
          ! "Zero entrapment": even radiation transported
          ! laterally between regions in the layers below is
          ! reflected back up into the same region. First diffuse
          ! radiation:
          total_albedo(:,:,:,jlev) = 0.0_jprb
          do jreg = 1,nreg    ! Target layer (jlev-1)
            do jreg2 = 1,nreg ! Current layer (jlev)
              total_albedo(:,jreg,jreg,jlev) = total_albedo(:,jreg,jreg,jlev) &
                   &  + sum(total_albedo_below(:,:,jreg2),2) &
                   &  * v_matrix(jreg2,jreg,jlev)
            end do
          end do
          ! ...then direct radiation:
          total_albedo_direct(:,:,:,jlev) = 0.0_jprb
          do jreg = 1,nreg    ! Target layer (jlev-1)
            do jreg2 = 1,nreg ! Current layer (jlev)
              total_albedo_direct(:,jreg,jreg,jlev) = total_albedo_direct(:,jreg,jreg,jlev) &
                   &  + sum(total_albedo_below_direct(:,:,jreg2),2) &
                   &  * v_matrix(jreg2,jreg,jlev)
            end do
          end do

        else
          ! Controlled entrapment

#ifdef EXPLICIT_EDGE_ENTRAPMENT
          ! If "EXPLICIT_EDGE_ENTRAPMENT" is defined then we use the
          ! explicit entrapment approach for both horizontal transport
          ! within regions, and horizontal transport between regions
          ! (otherwise, horizontal transport between regions is
          ! automatically treated using maximum entrapment). This is
          ! experimental, which is why it is not a run-time option.

          if (config%i_3d_sw_entrapment == IEntrapmentEdgeOnly) then
#endif
          ! Add the contribution from off-diagonal elements of the
          ! albedo matrix in the lower layer, i.e. radiation that
          ! flows between regions...

          ! First diffuse radiation:
          tmp = total_albedo_below
          do jreg = 1,nreg
            tmp(:,jreg,jreg) = 0.0_jprb
          end do 
          ! total_albedo(:,:,:,jlev) = singlemat_x_mat(ng,ng,nreg,&
          !     &  u_matrix(:,:,jlev), & ! temporary array  created for next argument
          !     &  mat_x_singlemat(ng,ng,nreg,tmp,&  
          !     &  v_matrix(:,:,jlev))) 
          tmp2 = mat_x_singlemat_3_sw(ng,tmp,v_matrix(:,:,jlev))
          total_albedo(:,:,:,jlev) = singlemat_x_mat_3_sw(ng,u_matrix(:,:,jlev), tmp2)
          ! ...then direct radiation:
          tmp = total_albedo_below_direct
          do jreg = 1,nreg
            tmp(:,jreg,jreg) = 0.0_jprb
          end do
          tmp2 = mat_x_singlemat_3_sw(ng,tmp,v_matrix(:,:,jlev))
          total_albedo_direct(:,:,:,jlev) = singlemat_x_mat_3_sw(ng,&
            &  u_matrix(:,:,jlev), tmp2)

#ifdef EXPLICIT_EDGE_ENTRAPMENT
end if
#endif
          ! Now the contribution from the diagonals of the albedo
          ! matrix in the lower layer
          if (config%i_3d_sw_entrapment == IEntrapmentEdgeOnly &
               &  .or. (.not. config%do_3d_effects)) then
            ! "Edge-only entrapment": the operation we perform is
            ! essentially diag(total_albedo) += matmul(transpose(v_matrix),
            ! diag(total_albedo_below)).
            do jreg = 1,nreg
              do jreg2 = 1,nreg
                total_albedo(:,jreg,jreg,jlev) &
                     &  = total_albedo(:,jreg,jreg,jlev) &
                     &  + total_albedo_below(:,jreg2,jreg2) &
                     &  * v_matrix(jreg2,jreg,jlev)
                total_albedo_direct(:,jreg,jreg,jlev) &
                     &  = total_albedo_direct(:,jreg,jreg,jlev) &
                     &  + total_albedo_below_direct(:,jreg2,jreg2) &
                     &  * v_matrix(jreg2,jreg,jlev)
              end do
            end do

          else
            associate(albedo_part=>entrapment)
            ! "Explicit entrapment"
! #ifndef USE_FAST_EXPM_EXCHANGE
            entrapment = 0.0_jprb
! #endif
            do jreg2 = 1,nreg
              ! Loop through each region in the lower layer. For one
              ! of the regions in the lower layer, we are imagining it
              ! to be divided into "nreg" subregions that map onto the
              ! regions in the upper layer. The rate of exchange
              ! between these subregions is computed via a coupled
              ! differential equation written in terms of a singular
              ! exchange matrix (there are only terms for the exchange
              ! between subregions, but no propagation effects since
              ! we already know the albedo of this region). This is
              ! solved using the matrix-exponential method.

              ! Use the following array for the transfer of either
              ! diffuse or direct radiation (despite the name), per
              ! unit horizontal distance travelled
              transfer_rate_diffuse = 0.0_jprb

              ! As we need to reference the layer above the interface,
              ! don't do the following on the highest layer
              if (jlev > 1) then

                ! Given a horizontal migration distance, there is
                ! still uncertainty about how much entrapment occurs
                ! associated with how one assumes cloud boundaries
                ! line up in adjacent layers. "overhang_factor"
                ! can be varied between 0.0 (the boundaries line up to
                ! the greatest extent possible given the overlap
                ! parameter) and 1.0 (the boundaries line up to the
                ! minimum extent possible); here this is used to
                ! produce a scaling factor for the transfer rate.
                transfer_scaling = 1.0_jprb - (1.0_jprb - config%overhang_factor) & 
                     &  * cloud%overlap_param(jcol,jlev-1) &
                     &  * min(region_fracs(jreg2,jlev,jcol),region_fracs(jreg2,jlev-1,jcol)) &
                     &  / max(config%cloud_fraction_threshold, region_fracs(jreg2,jlev,jcol))

                do jreg = 1, nreg-1
                  ! Compute lateral transfer rates from region jreg to
                  ! jreg+1 as before, but without length scale which
                  ! is wavelength dependent.

                  ! Recall that overlap indexing is
                  ! u_matrix(upper_region, lower_region, level, column).
                  transfer_rate_diffuse(jreg,jreg+1) = transfer_scaling &
                       &  * edge_length(jreg,jlev-1) / max(u_matrix(jreg,jreg2,jlev),1.0e-5_jprb)
                  ! Compute transfer rates from region jreg+1 to jreg
                  transfer_rate_diffuse(jreg+1,jreg) = transfer_scaling &
                       &  * edge_length(jreg,jlev-1) / max(u_matrix(jreg+1,jreg2,jlev),1.0e-5_jprb)
                end do
              
                ! Compute transfer rates directly between regions 1
                ! and 3 (not used below)
                if (edge_length(3,jlev) > 0.0_jprb) then
                  transfer_rate_diffuse(1,3) = transfer_scaling &
                       &  * edge_length(3,jlev-1) / max(u_matrix(1,jreg2,jlev),1.0e-5_jprb)
                  transfer_rate_diffuse(3,1) = transfer_scaling &
                       &  * edge_length(3,jlev-1) / max(u_matrix(3,jreg2,jlev),1.0e-5_jprb)
                end if

#ifdef USE_FAST_EXPM_EXCHANGE
                ! If we use fast_expm_exchange then ensure "a" and "d"
                ! are not equal as this causes problems in the
                ! eigenvector computation
                ! print *, abs(transfer_rate_diffuse(1,2)-transfer_rate_diffuse(3,2))
                if (abs(transfer_rate_diffuse(1,2)-transfer_rate_diffuse(3,2)) &
                    ! & < 100.0_jprb * epsilon(1.0_jprb) &
                      & < coeff &
                      &   * (transfer_rate_diffuse(1,2)+transfer_rate_diffuse(3,2))) then
                  transfer_rate_diffuse(1,2) = transfer_rate_diffuse(1,2) &
                        &  * (1.0_jprb - 400.0_jprb*epsilon(1.0_jprb)) * transfer_rate_diffuse(1,2)
                end if
#endif

              end if

              ! Compute matrix of exchange coefficients
              inv_effective_size = min(cloud%inv_cloud_effective_size(jcol,jlev-1), &
                   &                   1.0_jprb/config%min_cloud_effective_size)

              if (config%i_3d_sw_entrapment == IEntrapmentExplicit) then
                fractal_fac_dif = 1.0_jprb / sqrt(max(1.0_jprb, 2.5_jprb*x_diffuse(:,jreg2) &
                      &                                         * inv_effective_size))
                fractal_fac_dir = 1.0_jprb / sqrt(max(1.0_jprb, 2.5_jprb*x_direct(:,jreg2) &
                     &                                         * inv_effective_size))
              else
                fractal_fac_dif = 1.0_jprb
                fractal_fac_dir = 1.0_jprb
              end if

#ifdef USE_FAST_EXPM_EXCHANGE
!             ! only need e21, e12, e32, e23 for fastexpm and e11, e22 for the minvals, 
              ! but e11 is always just zero and e22 is just -e12
              ! jreg = 1
              ! e11 = 0 - 0 = 0
              do jg = 1, ng
                entrapment(jg,1,jreg2,2,1) = transfer_rate_diffuse(1,2)*x_diffuse(jg,jreg2) &
                    & * fractal_fac_dif(jg)
                entrapment(jg,2,jreg2,2,1) = transfer_rate_diffuse(1,2)*x_direct(jg,jreg2) & 
                    & * fractal_fac_dir(jg)
                entrapment(jg,1,jreg2,1,2) = transfer_rate_diffuse(2,1)*x_diffuse(jg,jreg2) &
                    & * fractal_fac_dif(jg)
                entrapment(jg,2,jreg2,1,2) = transfer_rate_diffuse(2,1)*x_direct(jg,jreg2) &
                    & * fractal_fac_dir(jg)
                ! e22 = 0 - e12
                ! jreg = 2jg
                ! e22 = e22 - e32 = e22 - 0 = -e12 (used in security below)
                entrapment(jg,1,jreg2,3,2) = transfer_rate_diffuse(2,3)*x_diffuse(jg,jreg2) &
                    & * fractal_fac_dif(jg)
                entrapment(jg,2,jreg2,3,2) = transfer_rate_diffuse(2,3)*x_direct(jg,jreg2) & 
                    & * fractal_fac_dir(jg)
                entrapment(jg,1,jreg2,2,3) = transfer_rate_diffuse(3,2)*x_diffuse(jg,jreg2) & 
                    & * fractal_fac_dif(jg)
                entrapment(jg,2,jreg2,2,3) = transfer_rate_diffuse(3,2)*x_direct(jg,jreg2) &
                    & * fractal_fac_dir(jg)
              end do 
              ! e33 : not needed ( equals -e23)
#else
              do jreg = 1,nreg-1
                ! Diffuse transport down and up with one random
                ! scattering event
                do jg = 1, ng
                  entrapment(jg,1,jreg2,jreg,jreg)   = entrapment(jg,1,jreg2,jreg,jreg) &
                      & - entrapment(jg,1,jreg2,jreg+1,jreg)
                  entrapment(jg,1,jreg2,jreg+1,jreg) = transfer_rate_diffuse(jreg,jreg+1) &
                      &  * x_diffuse(jg,jreg2) * fractal_fac_dif(jg)
                  entrapment(jg,1,jreg2,jreg,jreg+1) = transfer_rate_diffuse(jreg+1,jreg) &
                      &  * x_diffuse(jg,jreg2) * fractal_fac_dif(jg)
                  entrapment(jg,1,jreg2,jreg+1,jreg+1) = entrapment(jg,1,jreg2,jreg+1,jreg+1) &
                      & - entrapment(jg,1,jreg2,jreg,jreg+1)
                  entrapment(jg,2,jreg2,jreg,jreg)     = entrapment(jg,2,jreg2,jreg,jreg)     &
                      & - entrapment(jg,2,jreg2,jreg+1,jreg)
                  entrapment(jg,2,jreg2,jreg+1,jreg) = transfer_rate_diffuse(jreg,jreg+1) &
                      &  * x_direct(jg,jreg2) * fractal_fac_dif(jg)
                  entrapment(jg,2,jreg2,jreg,jreg+1) = transfer_rate_diffuse(jreg+1,jreg) &
                      &  * x_direct(jg,jreg2) * fractal_fac_dif(jg)
                  entrapment(jg,2,jreg2,jreg+1,jreg+1) = entrapment(jg,2,jreg2,jreg+1,jreg+1) &
                      & - entrapment(jg,2,jreg2,jreg,jreg+1)
                end do 
              end do
#endif
              ! If rate of exchange is excessive the expm can throw a
              ! floating point exception, even if it tends towards a
              ! trival limit, so we cap the maximum input to expm by
              ! scaling down if necessary
              do jg = 1,ng
                ! max_entr = -min(entrapment(jg,1,1),entrapment(jg,2,2),entrapment(jg,3,3))
                ! max_entr = -min(entrapment(jg,1,1),entrapment(jg,2,2)) 
                !             ^first argument is zero, second is -e12
                if (entrapment(jg,1,jreg2,1,2) > config%max_cloud_od_sw) then
                  ! Scale down all inputs for this g point
                  !entrapment(jg,1,:,:) = entrapment(jg,1,:,:) *(config%max_cloud_od_sw/max_entr)
                  entrapment(jg,1,jreg2,:,:) = entrapment(jg,1,jreg2,:,:) &
                      & * (config%max_cloud_od_sw/entrapment(jg,1,jreg2,1,2))
                end if
                if (entrapment(jg,2,jreg2,1,2) > config%max_cloud_od_sw) then
                  ! Scale down all inputs for this g point
                  entrapment(jg,2,jreg2,:,:) = entrapment(jg,2,jreg2,:,:) &
                    & * (config%max_cloud_od_sw/entrapment(jg,2,jreg2,1,2))
                end if
              end do     

            end do   
#ifdef USE_FAST_EXPM_EXCHANGE
            ! Since the matrix to be exponentiated has a simple
            ! structure we may use a faster method described in the
            ! appendix of Hogan et al. (GMD 2018)
            call fast_expm_exchange_3(2*nreg*ng, entrapment)

! #define DEBUG_FAST_EXPM 1
#ifdef DEBUG_FAST_EXPM
            if ((maxval(albedo_part)>1.5_jprb)) then 
              bad_col = jcol
              print *, "1 WARNING maxval alb>1.5 after fastexp", maxval(albedo_part)
              !return 
              if ((maxval(albedo_part)>100.0_jprb)) then
                if ((maxval(albedo_part)>1000.0_jprb)) then
                  bad_col = 1000+jcol
                  print *, "ERROR maxval alb>1000 after fastexp (FPE overflow)", maxval(albedo_part)
                  return
                else 
                  bad_col = 100+jcol
                  print *, "ERROR maxval alb>100 after fastexp (div by zero)", maxval(albedo_part)
                  return
                end if 
              end if 
              albedo_part = min(albedo_part, 1.0_jprb)
            end if 
#else 
            albedo_part = min(albedo_part, 1.0_jprb)
#endif
#else
            ! Use matrix exponential to compute rate of exchange
            call expm(2*nreg*ng, 2*nreg*ng, nreg, entrapment, IMatrixPatternDense)
            n_calls_expm = n_calls_expm + ng
#endif
            ! Scale to get the contribution to the diffuse albedo & increment diffuse albedo
            do jreg = 1,nreg
              do jreg3 = 1, nreg
                do jreg2 = 1,nreg
                  total_albedo(:,jreg3,jreg,jlev) = total_albedo(:,jreg3,jreg,jlev) &
                      & + albedo_part(:,1,jreg2,jreg3,jreg) &
                      & * v_matrix(jreg2,jreg,jlev) * total_albedo_below(:,jreg2,jreg2)
                  total_albedo_direct(:,jreg3,jreg,jlev) =total_albedo_direct(:,jreg3,jreg,jlev) + &
                      & + albedo_part(:,2,jreg2,jreg3,jreg) &
                      & * v_matrix(jreg2,jreg,jlev) * total_albedo_below_direct(:,jreg2,jreg2)
                end do
              end do 
            end do   
            end associate
          end if
        end if

        if ((config%i_3d_sw_entrapment == IEntrapmentExplicitNonFractal &
             &  .or. config%i_3d_sw_entrapment == IEntrapmentExplicit) &
             &  .and. .not. (is_clear_sky_layer(jlev) .and. is_clear_sky_layer(jlev-1))) then
          ! Horizontal migration distances are averaged when
          ! applying overlap rules, so equation is
          ! x_above=matmul(transpose(v_matrix),x_below)

          ! We do this into temporary arrays...
          x_direct_above = 0.0_jprb
          x_diffuse_above = 0.0_jprb
          
          nregactive = nreg
          if (is_clear_sky_layer(jlev)) then
            nregactive = 1
          end if

          do jreg = 1,nreg          ! Target layer (jlev-1)
            do jreg2 = 1,nregactive ! Current layer (jlev)
              x_direct_above(:,jreg) = x_direct_above(:,jreg) &
                   &  + x_direct(:,jreg2) * v_matrix(jreg2,jreg,jlev)
              x_diffuse_above(:,jreg) = x_diffuse_above(:,jreg) &
                   &  + x_diffuse(:,jreg2) * v_matrix(jreg2,jreg,jlev)
            end do
          end do
          
          !... then copy out of the temporary arrays
          x_direct = x_direct_above
          x_diffuse = x_diffuse_above
        end if

      end do ! Reverse loop over levels

      ! --------------------------------------------------------
      ! Section 5: Compute fluxes
      ! --------------------------------------------------------

      ! Top-of-atmosphere fluxes into the regions of the top-most
      ! layer, zero since we assume no diffuse downwelling
      flux_dn_below = 0.0_jprb
      ! Direct downwelling flux (into a plane perpendicular to the
      ! sun) entering the top of each region in the top-most layer
      do jreg = 1,nreg
        direct_dn_below(:,jreg) = incoming_sw(:,jcol)*region_fracs(jreg,1,jcol)
      end do
      ! We're using flux_up_above as a container; actually its
      ! interpretation at top of atmosphere here is just 'below' the
      ! TOA interface, so using the regions of the first model layer
      flux_up_above = mat_x_vec_3_sw(ng,total_albedo_direct(:,:,:,1),direct_dn_below)

      if (config%do_clear) then
        flux_dn_clear    = 0.0_jprb
        direct_dn_clear  = incoming_sw(:,jcol)
        flux_up_clear    = direct_dn_clear*total_albedo_clear_direct(:,1)
      end if

      ! Store the TOA broadband fluxes
      flux%sw_up(jcol,1) = sum(sum(flux_up_above,1))
      flux%sw_dn(jcol,1) = mu0 * sum(direct_dn_clear(:))
      if (allocated(flux%sw_dn_direct)) then
        flux%sw_dn_direct(jcol,1) = flux%sw_dn(jcol,1)
      end if
      if (config%do_clear) then
        flux%sw_up_clear(jcol,1) = sum(flux_up_clear)
        flux%sw_dn_clear(jcol,1) = flux%sw_dn(jcol,1)
        if (allocated(flux%sw_dn_direct_clear)) then
          flux%sw_dn_direct_clear(jcol,1) = flux%sw_dn_clear(jcol,1)
        end if
      end if

      ! Save the spectral fluxes if required
      if (config%do_save_spectral_flux) then
        call indexed_sum(sum(flux_up_above(:,:),2), &
             &           config%i_spec_from_reordered_g_sw, &
             &           flux%sw_up_band(:,jcol,1))
        call indexed_sum(sum(direct_dn_below(:,:),2), &
             &           config%i_spec_from_reordered_g_sw, &
             &           flux%sw_dn_band(:,jcol,1))
        flux%sw_dn_band(:,jcol,1) = mu0 * flux%sw_dn_band(:,jcol,1)
        if (allocated(flux%sw_dn_direct_band)) then
          flux%sw_dn_direct_band(:,jcol,1) = flux%sw_dn_band(:,jcol,1)
        end if
        if (config%do_clear) then
          flux%sw_dn_clear_band(:,jcol,1) = flux%sw_dn_band(:,jcol,1)
          call indexed_sum(flux_up_clear, &
               &           config%i_spec_from_reordered_g_sw, &
               &           flux%sw_up_clear_band(:,jcol,1))
          if (allocated(flux%sw_dn_direct_clear_band)) then
            flux%sw_dn_direct_clear_band(:,jcol,1) &
                 &   = flux%sw_dn_clear_band(:,jcol,1)
          end if
        end if
      end if

      ! Final loop back down through the atmosphere to compute fluxes
      do jlev = 1,nlev

#ifdef PRINT_ENTRAPMENT_DATA
        if (config%i_3d_sw_entrapment == IEntrapmentExplicitNonFractal &
             &  .or. config%i_3d_sw_entrapment == IEntrapmentExplicit) then
          ! Save downwelling direct and diffuse fluxes at the top of
          ! layer "jlev" in each of the regions of layer "jlev" 
          write(102,'(i4,i4,6e14.6)') jcol, jlev, direct_dn_below(1,1:3), flux_dn_below(1,1:3)
        end if
#endif

        if (config%do_clear) then
          ! Compute the solar downwelling "source" at the base of the
          ! layer due to scattering of the direct beam within it
          do jg = 1,ng
          source_dn_clear(jg) = trans_dir_diff_clear(jg,jlev)*direct_dn_clear(jg)
          ! Compute direct downwelling flux in each region at base of
          ! current layer
          direct_dn_clear(jg) = trans_dir_dir_clear(jg,jlev)*direct_dn_clear(jg)
          end do 
        end if

        ! source_dn(:,:) = mat_x_vec_3_sw(ng,trans_dir_diff(:,:,:,jlev),direct_dn_below, &
        !     &  is_clear_sky_layer(jlev))
        ! direct_dn_above = mat_x_vec_3_sw(ng,trans_dir_dir(:,:,:,jlev),direct_dn_below, &
        !     &  is_clear_sky_layer(jlev))
        associate(A1=>trans_dir_diff(:,:,:,jlev), A2=>trans_dir_dir(:,:,:,jlev), &
            & b=>direct_dn_below)
        if (is_clear_sky_layer(jlev)) then
          source_dn             = 0.0_jprb
          source_dn(:,1)        = A1(:,1,1)*b(:,1)
          direct_dn_above       = 0.0_jprb
          direct_dn_above(:,1)  = A2(:,1,1)*b(:,1)
        else
          do jg = 1, ng
            source_dn(jg,1) = A1(jg,1,1)*b(jg,1) + A1(jg,1,2)*b(jg,2) + A1(jg,1,3)*b(jg,3)
            source_dn(jg,2) = A1(jg,2,1)*b(jg,1) + A1(jg,2,2)*b(jg,2) + A1(jg,2,3)*b(jg,3)
            source_dn(jg,3) = A1(jg,3,1)*b(jg,1) + A1(jg,3,2)*b(jg,2) + A1(jg,3,3)*b(jg,3)
            direct_dn_above(jg,1) = A2(jg,1,1)*b(jg,1) + A2(jg,1,2)*b(jg,2) + A2(jg,1,3)*b(jg,3) 
            direct_dn_above(jg,2) = A2(jg,2,1)*b(jg,1) + A2(jg,2,2)*b(jg,2) + A2(jg,2,3)*b(jg,3) 
            direct_dn_above(jg,3) = A2(jg,3,1)*b(jg,1) + A2(jg,3,2)*b(jg,2) + A2(jg,3,3)*b(jg,3) 
          end do
        end if
        end associate

        if (config%do_clear) then
          ! Scalar operations for clear-sky fluxes
          do jg = 1, ng
            flux_dn_clear(jg) = (trans_clear(jg,jlev)*flux_dn_clear(jg) &
                &  + ref_clear(jg,jlev)*total_albedo_clear_direct(jg,jlev+1)*direct_dn_clear(jg) &
                &  + source_dn_clear(jg)) &
                &  / (1.0_jprb - ref_clear(jg,jlev)*total_albedo_clear(jg,jlev+1))
            flux_up_clear(jg) = total_albedo_clear_direct(jg,jlev+1)*direct_dn_clear(jg) &
                &  + total_albedo_clear(jg,jlev+1)*flux_dn_clear(jg)
          end do 
        end if

        if (is_clear_sky_layer(jlev)) then
          ! Scalar operations for clear-sky layer
          flux_dn_above = 0.0_jprb
          flux_up_above = 0.0_jprb
          do jg = 1,ng
            flux_dn_above(jg,1) = (trans_clear(jg,jlev)*flux_dn_below(jg,1) &
                &  + ref_clear(jg,jlev)*total_albedo_direct(jg,1,1,jlev+1)*direct_dn_above(jg,1) &
                &  + source_dn(jg,1)) &
                &  / (1.0_jprb - ref_clear(jg,jlev)*total_albedo(jg,1,1,jlev+1))     
            flux_up_above(jg,1) = total_albedo_direct(jg,1,1,jlev+1)*direct_dn_above(jg,1) &
                &  + total_albedo(jg,1,1,jlev+1)*flux_dn_above(jg,1)
          end do 
        else
          ! Matrix operations for cloudy layer
          call identity_minus_mat_x_mat_3_sw(ng,reflectance(:,:,:,jlev), &
              &  total_albedo(:,:,:,jlev+1), denominator)

          total_source = mat_x_vec_3_sw(ng,total_albedo_direct(:,:,:,jlev+1),direct_dn_above)

          flux_dn_above = solve_vec_3_sw(ng,denominator, &
              &  mat_x_vec_3_sw(ng,transmittance(:,:,:,jlev),flux_dn_below) &
              &  + mat_x_vec_3_sw(ng,reflectance(:,:,:,jlev), total_source(:,:)) &
              &  + source_dn(:,:))

          flux_up_above = mat_x_vec_3_sw(ng,total_albedo(:,:,:,jlev+1), &
              &  flux_dn_above) + total_source(:,:)
        end if

        ! Account for overlap rules in translating fluxes just above
        ! a layer interface to the values just below
        if (is_clear_sky_layer(jlev) .and. is_clear_sky_layer(jlev+1)) then
          ! Regions in current layer map directly on to regions in
          ! layer below
          flux_dn_below = flux_dn_above
          direct_dn_below = direct_dn_above
        else
          ! Apply downward overlap matrix to compute direct
          ! downwelling flux entering the top of each region in the
          ! layer below
          flux_dn_below = singlemat_x_vec_sw(ng,v_matrix(:,:,jlev+1), &
               &  flux_dn_above)
          direct_dn_below = singlemat_x_vec_sw(ng,v_matrix(:,:,jlev+1), &
               &  direct_dn_above)
        end if

        ! Compute and store the broadband fluxes
        sums_up = 0.0_jprb; sums_dn = 0.0_jprb; sums_dn_dir = 0.0_jprb
        if (is_clear_sky_layer(jlev)) then 
#ifdef __NEC__
          !NEC$ shortloop
#else
          !$omp simd reduction(+:sums_up, sums_dn, sums_dn_dir)
#endif
          do jg = 1, ng  
            sums_up = sums_up + flux_up_above(jg,1) 
            sums_dn = sums_dn + flux_dn_above(jg,1) 
            sums_dn_dir = sums_dn_dir + direct_dn_above(jg,1)
          end do
        else
#ifdef __NEC__
          !NEC$ shortloop
#else
          !$omp simd reduction(+:sums_up, sums_dn, sums_dn_dir)
#endif
          do jg = 1, ng  
            sums_up = sums_up + flux_up_above(jg,1) + flux_up_above(jg,2) + flux_up_above(jg,3)
            sums_dn = sums_dn + flux_dn_above(jg,1) + flux_dn_above(jg,2) + flux_dn_above(jg,3)
            sums_dn_dir = sums_dn_dir + direct_dn_above(jg,1) + direct_dn_above(jg,2) + direct_dn_above(jg,3)
          end do
        end if 
        flux%sw_up(jcol,jlev+1) = sums_up 
        flux%sw_dn(jcol,jlev+1) = mu0*sums_dn_dir + sums_dn
        if (allocated(flux%sw_dn_direct)) flux%sw_dn_direct(jcol,jlev+1) = mu0*sums_dn_dir
        if (config%do_clear) then
          sums_up_clear = 0.0_jprb; sums_dn_clear = 0.0_jprb; sums_dn_dir_clear = 0.0_jprb
#ifdef __NEC__
          !NEC$ shortloop
#else
          !$omp simd reduction(+:sums_up_clear, sums_dn_clear, sums_dn_dir_clear)
#endif
          do jg = 1, ng  
            sums_up_clear = sums_up_clear + flux_up_clear(jg)
            sums_dn_clear = sums_dn_clear + flux_dn_clear(jg)
            sums_dn_dir_clear = sums_dn_dir_clear + direct_dn_clear(jg)
          end do
          flux%sw_up_clear(jcol,jlev+1) = sums_up_clear 
          flux%sw_dn_clear(jcol,jlev+1) = mu0*sums_dn_dir_clear + sums_dn_clear
          if (allocated(flux%sw_dn_direct_clear)) then
            flux%sw_dn_direct_clear(jcol,jlev+1) =  mu0*sums_dn_dir_clear
          end if
        end if

        ! Save the spectral fluxes if required
        if (config%do_save_spectral_flux) then
          ! Down
          call indexed_sum(sum(direct_dn_above,2), &
               &           config%i_spec_from_reordered_g_sw, &
               &           flux%sw_dn_band(:,jcol,jlev+1))
          flux%sw_dn_band(:,jcol,jlev+1) = mu0 * flux%sw_dn_band(:,jcol,jlev+1)

          if (allocated(flux%sw_dn_direct_band)) then
            flux%sw_dn_direct_band(:,jcol,jlev+1) &
                 &   = flux%sw_dn_band(:,jcol,jlev+1)
          end if
          if (config%do_clear) then
            call indexed_sum(direct_dn_clear, &
                 &           config%i_spec_from_reordered_g_sw, &
                 &           flux%sw_dn_clear_band(:,jcol,jlev+1))
            flux%sw_dn_clear_band(:,jcol,jlev+1) = mu0 &
                 &   * flux%sw_dn_clear_band(:,jcol,jlev+1)
            if (allocated(flux%sw_dn_direct_clear_band)) then
              flux%sw_dn_direct_clear_band(:,jcol,jlev+1) &
                   &  = flux%sw_dn_clear_band(:,jcol,jlev+1)
            end if
          end if
          ! Up
          call indexed_sum(sum(flux_up_above,2), &
               &           config%i_spec_from_reordered_g_sw, &
               &           flux%sw_up_band(:,jcol,jlev+1))
          call add_indexed_sum(sum(flux_dn_above,2), &
               &           config%i_spec_from_reordered_g_sw, &
               &           flux%sw_dn_band(:,jcol,jlev+1))
          if (config%do_clear) then
            call indexed_sum(flux_up_clear, &
                 &           config%i_spec_from_reordered_g_sw, &
                 &           flux%sw_up_clear_band(:,jcol,jlev+1))
            call add_indexed_sum(flux_dn_clear, &
                 &           config%i_spec_from_reordered_g_sw, &
                 &           flux%sw_dn_clear_band(:,jcol,jlev+1))
          end if
        end if

      end do ! Final loop over levels

      ! Store surface spectral fluxes, if required (after the end of
      ! the final loop over levels, the current values of these arrays
      ! will be the surface values)
      flux%sw_dn_diffuse_surf_g(:,jcol) = sum(flux_dn_above,2)
      flux%sw_dn_direct_surf_g(:,jcol)  = mu0 * sum(direct_dn_above,2)
      if (config%do_clear) then
        flux%sw_dn_diffuse_surf_clear_g(:,jcol) = flux_dn_clear
        flux%sw_dn_direct_surf_clear_g(:,jcol)  = mu0 * direct_dn_clear
      end if
      
    end do ! Loop over columns
    
    if (config%iverbose >= 3) then
      write(nulout,*)
    end if

    ! Report number of calls to each method of solving single-layer
    ! two-stream equations
    if (config%iverbose >= 4) then
      write(nulout,'(a,i0)') '  Matrix-exponential calls: ', n_calls_expm
      write(nulout,'(a,i0)') '  Meador-Weaver calls: ', n_calls_meador_weaver
    end if

    if (lhook) call dr_hook('radiation_spartacus_sw:solver_spartacus_sw',1,hook_handle)

  end subroutine solver_spartacus_sw


  ! Step the horizontal migration distances from the base of a layer
  ! to the top, accounting for the extra distance travelled within the
  ! layer
  pure subroutine step_migrations(ng_sw_in, nreg, cloud_frac, &
     &  layer_depth, tan_diffuse_angle_3d, tan_sza, &
     &  reflectance, transmittance, ref_dir, trans_dir_dir, &
     &  trans_dir_diff, total_albedo_diff, total_albedo_dir, &
     &  x_diffuse, x_direct)
  
  use parkind1, only : jprb

  implicit none

  ! Inputs

  ! Number of g points and regions
  integer, intent(in) :: ng_sw_in, nreg
  ! Cloud fraction
  real(jprb), intent(in) :: cloud_frac
  ! Layer depth (m), tangent of diffuse zenith angle and tangent of
  ! solar zenith angle
  real(jprb), intent(in) :: layer_depth, tan_diffuse_angle_3d, tan_sza
  ! Reflectance and transmittance to diffuse downwelling radiation
  real(jprb), intent(in), dimension(ng, nreg, nreg) :: reflectance, transmittance
  ! Reflectance and transmittance to direct downwelling radiation
  real(jprb), intent(in), dimension(ng, nreg, nreg) :: ref_dir, trans_dir_dir
  ! Transmittance involving direct entering a layer from the top and
  ! diffuse leaving from the bottom
  real(jprb), intent(in), dimension(ng, nreg, nreg) :: trans_dir_diff

  ! Total albedo of direct and diffuse radiation of the atmosphere
  ! below the layer in question
  real(jprb), intent(in), dimension(ng, nreg, nreg) &
       &  :: total_albedo_diff, total_albedo_dir

  ! Inputs/outputs

  ! Horizontal migration distance (m) of reflected light
  real(jprb), intent(inout), dimension(ng, nreg) :: x_diffuse, x_direct

  ! Local variables

  ! Top albedo, i.e. the albedo of the top of the layer assuming no
  ! lateral transport
  real(jprb) :: top_albedo

  ! Multiple-scattering amplitude enhancement
  real(jprb) :: ms_enhancement

  ! Multiple-scattering distance enhancement
  real(jprb) :: x_enhancement


  real(jprb) :: x_layer_diffuse, x_layer_direct, one_min_ref_x_albedo
  integer :: jreg, istartreg, iendreg, jg

  istartreg = 1
  iendreg   = nreg

  if (cloud_frac <= 0.0_jprb) then
    ! Clear-sky layer: don't waste time on cloudy regions
    iendreg = 1
  else if (cloud_frac >= 1.0_jprb) then
    ! Overcast layer: don't waste time on clear region
    istartreg = 2
  end if

  ! This is the mean horizontal distance travelled by diffuse
  ! radiation that travels from the top of a layer to the centre and
  ! is then scattered back up and out
  x_layer_diffuse = layer_depth * tan_diffuse_angle_3d/sqrt(2.0_jprb) 

  ! This is the mean horizontal distance travelled by direct
  ! radiation that travels from the top of a layer to the centre and
  ! is then scattered back up and out
  x_layer_direct  = layer_depth * sqrt(tan_sza*tan_sza &
       &                             + tan_diffuse_angle_3d*tan_diffuse_angle_3d) * 0.5_jprb

  do jreg = istartreg,iendreg
    do jg = 1, ng
    ! Geometric series enhancement due to multiple scattering: the
    ! amplitude enhancement is equal to the limit of
    ! T*[1+RA+(RA)^2+(RA)^3+...]

      one_min_ref_x_albedo = 1.0_jprb - reflectance(jg,jreg,jreg)*total_albedo_diff(jg,jreg,jreg)
      ms_enhancement = transmittance(jg,jreg,jreg) / one_min_ref_x_albedo

      ! ...and the distance enhancement is approximately equal to the
      ! limit of T*[1+sqrt(2)*RA+sqrt(3)*(RA)^2+sqrt(4)*(RA)^3+...]

      ! x_enhancement = one_min_ref_x_albedo**(-1.5_jprb)
      ! pow() is very slow, luckily we can exploit that X**-1.5 = 1 / X**(1+0.5) = 1 / (X**1 * X**0.5) 
      x_enhancement = 1 / (one_min_ref_x_albedo * sqrt(one_min_ref_x_albedo))

      ! Horizontal migration of direct downwelling radiation
      top_albedo = max(1.0e-8_jprb, ref_dir(jg,jreg,jreg) + ms_enhancement &
          &  * (trans_dir_diff(jg,jreg,jreg)*total_albedo_diff(jg,jreg,jreg) &
          &     +trans_dir_dir(jg,jreg,jreg)*total_albedo_dir(jg,jreg,jreg)))
      ! The following is approximate and has been found to
      ! occasionally go negative
      x_direct(jg,jreg) = max(0.0_jprb, x_layer_direct &
          &  + ((trans_dir_diff(jg,jreg,jreg)*total_albedo_diff(jg,jreg,jreg)*x_enhancement &
          &      +trans_dir_dir(jg,jreg,jreg)*total_albedo_dir(jg,jreg,jreg)*(x_enhancement-1.0_jprb)) &
          &     *(x_diffuse(jg,jreg)+x_layer_diffuse) &
          &    +trans_dir_dir(jg,jreg,jreg)*total_albedo_dir(jg,jreg,jreg) &
          &     *(x_direct(jg,jreg)+x_layer_direct)) &
          &    * transmittance(jg,jreg,jreg) / top_albedo)

      ! Horizontal migration of diffuse downwelling radiation
      top_albedo = max(1.0e-8_jprb, reflectance(jg,jreg,jreg) &
          &  + ms_enhancement*transmittance(jg,jreg,jreg)*total_albedo_diff(jg,jreg,jreg))
      x_diffuse(jg,jreg) = x_layer_diffuse + x_enhancement*total_albedo_diff(jg,jreg,jreg) &
          &  *(transmittance(jg,jreg,jreg)*transmittance(jg,jreg,jreg)) &
          &  * (x_diffuse(jg,jreg) + x_layer_diffuse) / top_albedo
    end do
  end do
  if (iendreg < nreg) then
    x_diffuse(:,iendreg+1:nreg)      = 0.0_jprb
    x_direct(:,iendreg+1:nreg)       = 0.0_jprb
  else if (istartreg == 2) then
    x_diffuse(:,1)      = 0.0_jprb
    x_direct(:,1)       = 0.0_jprb
  end if

end subroutine step_migrations

pure subroutine write_gamma_diag(ng_sw_in, nlev, nreg, jreg, od_region, &
  &   gamma1, gamma2, gamma3, ssa, one_over_mu0, Gamma_z1)

  use parkind1, only : jprb

  ! Inputs
  ! Number of g points times levels to process (collapsed dimension)
  integer, intent(in) :: ng_sw_in, nlev
  ! Number of regions and region index 
  integer, intent(in) :: nreg, jreg

  real(jprb), intent(in), dimension(ng*nlev) :: od_region, gamma1, gamma2, gamma3, ssa
  real(jprb), intent(in) :: one_over_mu0

  real(jprb), intent(inout), dimension(ng*nlev, 3*nreg, 3*nreg) :: Gamma_z1

  integer :: jg

  do jg = 1, ng*nlev
    ! Write the diagonal elements of -Gamma1*z1
    Gamma_z1(jg,jreg,jreg) = od_region(jg)*gamma1(jg)
    ! Write the diagonal elements of +Gamma2*z1
    Gamma_z1(jg,jreg+nreg,jreg) = od_region(jg)*gamma2(jg)
    ! Write the diagonal elements of -Gamma3*z1
    Gamma_z1(jg,jreg,jreg+2*nreg) &
    &  = -od_region(jg)*ssa(jg) * gamma3(jg)

    ! Write the diagonal elements of +Gamma4*z1
    Gamma_z1(jg,jreg+nreg,jreg+2*nreg) &
    &  = od_region(jg)*ssa(jg) * (1.0_jprb - gamma3(jg))

    ! Write the diagonal elements of +Gamma0*z1
    Gamma_z1(jg,jreg+2*nreg,jreg+2*nreg) = -od_region(jg)*one_over_mu0
  end do
    
end subroutine write_gamma_diag

end module radiation_spartacus_sw