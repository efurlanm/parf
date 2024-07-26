! radiation_tripleclouds_sw.F90 - Shortwave "Tripleclouds" solver
!
! (C) Copyright 2016- ECMWF.
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
!   2017-10-23  R. Hogan  Renamed single-character variables
!   2018-10-08  R. Hogan  Call calc_region_properties
!   2019-01-02  R. Hogan  Fixed problem of do_save_spectral_flux .and. .not. do_sw_direct
!   2020-09-18  R. Hogan  Replaced some array expressions with loops for speed
!   2022-09-01  P. Ukkonen  Optimizations for much better performance with ECCKD, including: 
!                           batching section 3 computations, faster kernels, ng can be defined at compile time

module radiation_tripleclouds_sw

  public

contains
  ! Provides elemental function "delta_eddington"
#include "radiation_delta_eddington.h"

  ! Small routine for scaling cloud optical depth in the cloudy
  ! regions
#include "radiation_optical_depth_scaling.h"

  !---------------------------------------------------------------------
  ! This module contains just one subroutine, the shortwave
  ! "Tripleclouds" solver in which cloud inhomogeneity is treated by
  ! dividing each model level into three regions, one clear and two
  ! cloudy (with differing optical depth). This approach was described
  ! by Shonk and Hogan (2008).
  subroutine solver_tripleclouds_sw(ng_sw_in, nlev,istartcol,iendcol, &
       &  config, single_level, cloud, & 
       &  od, ssa, g, od_cloud, ssa_cloud, g_cloud, &
       &  albedo_direct, albedo_diffuse, incoming_sw, &
       &  flux)

    use parkind1, only           : jprb
    use yomhook,  only           : lhook, dr_hook, jphook

!    use radiation_io, only             : nulout
    use radiation_config, only         : config_type, IPdfShapeGamma
    use radiation_single_level, only   : single_level_type
    use radiation_cloud, only          : cloud_type
    use radiation_regions, only        : calc_region_properties
    use radiation_overlap, only        : calc_overlap_matrices_nocol
    use radiation_flux, only           : flux_type, &
         &                               indexed_sum, add_indexed_sum
    use radiation_matrix, only         : singlemat_x_vec, singlemat_x_vec_sw
    use radiation_two_stream, only     : calc_ref_trans_sw
#ifdef USE_TIMING
    ! Timing library
    use gptl,                  only: gptlstart, gptlstop
#endif
    implicit none

! Allow size of inner dimension (number of g-points) to be known at compile time if NG_SW is defined
#ifdef NG_SW
integer, parameter :: ng = NG_SW
#else
#define ng ng_sw_in
#endif

    ! Inputs
    integer, intent(in) :: ng_sw_in           ! number of g-points
    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type),        intent(in) :: config
    type(single_level_type),  intent(in) :: single_level
    type(cloud_type),         intent(in) :: cloud

    ! Gas and aerosol optical depth, single-scattering albedo and
    ! asymmetry factor at each shortwave g-point
    real(jprb), intent(in), dimension(ng,nlev,istartcol:iendcol) :: &
         &  od, ssa, g

    ! Cloud and precipitation optical depth, single-scattering albedo and
    ! asymmetry factor in each shortwave band
    real(jprb), intent(in), dimension(config%n_bands_sw,nlev,istartcol:iendcol) :: &
         &  od_cloud, ssa_cloud, g_cloud
    ! Optical depth, single scattering albedo and asymmetry factor in
    ! each g-point including gas, aerosol and clouds
    ! real(jprb), dimension(ng) :: od_total, ssa_total, g_total

    ! Now allocatable due to batching computations for cloudy regions
    real(jprb), dimension(:,:,:), allocatable, target &
         &  :: od_tot_batch, ssa_tot_batch, g_tot_batch

    ! Direct and diffuse surface albedos, and the incoming shortwave
    ! flux into a plane perpendicular to the incoming radiation at
    ! top-of-atmosphere in each of the shortwave g points
    real(jprb), intent(in), dimension(ng,istartcol:iendcol) :: &
         &  albedo_direct, albedo_diffuse, incoming_sw

    ! Output
    type(flux_type), intent(inout):: flux

    ! Local constants
    integer, parameter :: nregions = 3

    ! In a clear-sky layer this will be 1, otherwise equal to nregions
    ! integer :: nreg

    ! Local variables
    
    ! The area fractions of each region
    real(jprb) :: region_fracs(1:nregions,nlev,istartcol:iendcol)

    ! The scaling used for the optical depth in the cloudy regions
    real(jprb) :: od_scaling(2:nregions,nlev,istartcol:iendcol)

    ! Directional overlap matrices defined at all layer interfaces
    ! including top-of-atmosphere and the surface
    real(jprb), dimension(nregions,nregions,nlev+1) ::  v_matrix !, u_matrix

    ! Two-stream variables
    ! Diffuse reflection and transmission matrices of each layer
    real(jprb), dimension(ng, 2:nregions, nlev), target &
         &  :: reflectance, transmittance

    ! Terms translating the direct flux entering the layer from above
    ! to the reflected radiation exiting upwards (ref_dir) and the
    ! scattered radiation exiting downwards (trans_dir_diff), along with the
    ! direct unscattered transmission matrix (trans_dir_dir).
    real(jprb), dimension(ng, 2:nregions, nlev), target &
         &  :: ref_dir, trans_dir_diff, trans_dir_dir

    real(jprb), dimension(ng, nlev) &
    	 &  :: reflectance_clear, transmittance_clear, &
         & ref_dir_clear, trans_dir_diff_clear, trans_dir_dir_clear

    ! Total albedo of the atmosphere/surface just above a layer
    ! interface with respect to downwelling diffuse and direct
    ! (respectively) radiation at that interface, where level index =
    ! 1 corresponds to the top-of-atmosphere
    real(jprb), dimension(ng, nregions, nlev+1) &
         &  :: total_albedo, total_albedo_direct

    ! ...equivalent values for clear-skies
    real(jprb), dimension(ng, nlev+1) &
         &  :: total_albedo_clear, total_albedo_clear_direct

    ! Total albedo of the atmosphere just below a layer interface
    real(jprb), dimension(ng, nregions) &
         &  :: total_albedo_below, total_albedo_below_direct

    ! Direct downwelling flux below and above an interface between
    ! layers into a plane perpendicular to the direction of the sun
    real(jprb), dimension(ng, nregions) &
         &  :: direct_dn!, direct_dn_below
    ! Diffuse equivalents
    real(jprb), dimension(ng, nregions) &
         &  :: flux_dn, flux_up!, flux_dn_below

    ! ...clear-sky equivalent (no distinction between "above/below")
    real(jprb), dimension(ng) &
         &  :: direct_dn_clear, flux_dn_clear, flux_up_clear

    ! Clear-sky equivalent, but actually its reciprocal to replace
    ! some divisions by multiplications
    real(jprb), dimension(ng, nregions) :: inv_denom

    ! Identify clear-sky layers, with pseudo layers for outer space
    ! and below the ground, both treated as single-region clear skies
    logical :: is_clear_sky_layer(0:nlev+1), are_clouds_below

    ! Scattering optical depth of gas+aerosol and of cloud
    real(jprb) :: scat_od, scat_od_cloud

    real(jprb) :: mu0

    ! Faster broadband flux computation 
    real(jprb) :: sums_up, sums_dn, sums_dn_dir, sums_up_clear, sums_dn_clear, sums_dn_dir_clear

    integer :: jcol, jlev, jg, jreg, jband, jreg2, jtop, jbot, nlev_cloud
#ifdef USE_TIMING
    integer :: ret
#endif
    ! integer :: nlev_block, nblocks, istartlev, iendlev, jb
    real(jphook) :: hook_handle
  

    if (lhook) call dr_hook('radiation_tripleclouds_sw:solver_tripleclouds_sw',0,hook_handle)

    ! --------------------------------------------------------
    ! Section 1: Prepare general variables and arrays
    ! --------------------------------------------------------
    ! Initialize so compiler doesn't complain
    jbot = 0

    ! Compute the wavelength-independent region fractions and
    ! optical-depth scalings
    call calc_region_properties(nlev,nregions,istartcol,iendcol, &
         &  config%i_cloud_pdf_shape == IPdfShapeGamma, &
         &  cloud%fraction, cloud%fractional_std, region_fracs, &
         &  od_scaling, config%cloud_fraction_threshold)

    ! Main loop over columns
    do jcol = istartcol, iendcol

      ! Define which layers contain cloud; assume that
      ! cloud%crop_cloud_fraction has already been called
      is_clear_sky_layer = .true.
      do jlev = 1,nlev
        if (cloud%fraction(jcol,jlev) > 0.0_jprb) then
          is_clear_sky_layer(jlev) = .false.
        end if
      end do

      ! Compute wavelength-independent overlap matrix v_matrix
      call calc_overlap_matrices_nocol(nlev,nregions, is_clear_sky_layer, &
      &  region_fracs(:,:,jcol), cloud%overlap_param(jcol,:), &
      &  v_matrix, decorrelation_scaling=config%cloud_inhom_decorr_scaling, &
      &  cloud_fraction_threshold=config%cloud_fraction_threshold, &
      &  use_beta_overlap=config%use_beta_overlap, &
      &  cloud_cover=flux%cloud_cover_sw(jcol))
      ! --------------------------------------------------------
      ! Section 2: Prepare column-specific variables and arrays
      ! --------------------------------------------------------

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


      ! --------------------------------------------------------
      ! Section 3: Loop over layers to compute reflectance and transmittance
      ! --------------------------------------------------------
      ! In this section the reflectance, transmittance and sources
      ! are computed for each layer

      ! ...clear-sky equivalents
      call calc_ref_trans_sw(ng,nlev, &
          &  mu0, od(:,:,jcol), ssa(:,:,jcol), g(:,:,jcol), &
          &  reflectance_clear, transmittance_clear, &
          &  ref_dir_clear, trans_dir_diff_clear, &
          &  trans_dir_dir_clear )  

      ! Cloudy-sky equivalents
      ! Optimized computations for cloudy sky: batch together consecutively cloudy layers
      ! For this we need the cloud boundaries i.e. top and bottom indices
 
      are_clouds_below = .false.
      !  jtop = findloc(is_clear_sky_layer(1:nlev), .false.,dim=1) ! returns index of first cloudy layer

      ! findloc not working for nvfortran. manual implementation:
      jtop = 0
      do jlev = 1, nlev
          if (.not.is_clear_sky_layer(jlev)) then
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
            if (is_clear_sky_layer(jlev)) then
              jbot = jlev-1 ! bottom found, exit loop
              exit
            else if (jlev==nlev) then
              jbot = nlev
            end if
          end do
        end if

        nlev_cloud = jbot - jtop + 1
        allocate(od_tot_batch(ng,2,jtop:jbot))
        allocate(ssa_tot_batch(ng,2,jtop:jbot))
        allocate(g_tot_batch(ng,2,jtop:jbot))

        do jlev = jtop, jbot
          do jreg = 2, nregions
            if (ng == config%n_bands_sw) then 
              ! ECCKD, simpler loop indices (can vectorize better)
              do jg = 1,ng
                ! Mapping from g-point to band
                scat_od = od(jg,jlev,jcol)*ssa(jg,jlev,jcol)
                scat_od_cloud = od_cloud(jg,jlev,jcol) &
                    &  * ssa_cloud(jg,jlev,jcol) * od_scaling(jreg,jlev,jcol)
                ! Add scaled cloud optical depth to clear-sky value
                od_tot_batch(jg,jreg-1,jlev) = od(jg,jlev,jcol) &
                    &  + od_cloud(jg,jlev,jcol)*od_scaling(jreg,jlev,jcol)
                ! Compute single-scattering albedo and asymmetry
                ! factor of gas-cloud combination
                ssa_tot_batch(jg,jreg-1,jlev) = (scat_od+scat_od_cloud) &
                    &  / od_tot_batch(jg,jreg-1,jlev)
                g_tot_batch(jg,jreg-1,jlev) = (scat_od*g(jg,jlev,jcol) &
                    &         + scat_od_cloud * g_cloud(jg,jlev,jcol)) &
                    &      / (scat_od + scat_od_cloud)
              end do
            else
              do jg = 1,ng
                ! Mapping from g-point to band
                jband = config%i_band_from_reordered_g_sw(jg)
                scat_od = od(jg,jlev,jcol)*ssa(jg,jlev,jcol)
                scat_od_cloud = od_cloud(jband,jlev,jcol) &
                    &  * ssa_cloud(jband,jlev,jcol) * od_scaling(jreg,jlev,jcol)
                ! Add scaled cloud optical depth to clear-sky value
                od_tot_batch(jg,jreg-1,jlev) = od(jg,jlev,jcol) &
                    &  + od_cloud(jband,jlev,jcol)*od_scaling(jreg,jlev,jcol)
                ! Compute single-scattering albedo and asymmetry
                ! factor of gas-cloud combination
                ssa_tot_batch(jg,jreg-1,jlev) = (scat_od+scat_od_cloud) &
                    &  / od_tot_batch(jg,jreg-1,jlev)
                g_tot_batch(jg,jreg-1,jlev) = (scat_od*g(jg,jlev,jcol) &
                    &         + scat_od_cloud * g_cloud(jband,jlev,jcol)) &
                    &      / (scat_od + scat_od_cloud)
              end do
            end if
          end do
        end do

        call calc_ref_trans_sw(ng,2*nlev_cloud, &
           &  mu0, od_tot_batch, ssa_tot_batch, g_tot_batch, &
           &  reflectance(:,:,jtop:jbot), transmittance(:,:,jtop:jbot), &
           &  ref_dir(:,:,jtop:jbot), trans_dir_diff(:,:,jtop:jbot), trans_dir_dir(:,:,jtop:jbot))

        deallocate(od_tot_batch,ssa_tot_batch,g_tot_batch)
 
        if (jbot== nlev) are_clouds_below=.false.

        ! does another cloudy layer exist?
        if (any(.not. is_clear_sky_layer(jbot+1:nlev))) then
             ! find the cloud top 
             do jlev = jbot+1, nlev  
                  if (.not. is_clear_sky_layer(jlev)) then
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

      ! total_albedo(:,:,:) = 0.0_jprb
      ! total_albedo_direct(:,:,:) = 0.0_jprb
      total_albedo(:,:,nlev+1) = 0.0_jprb
      total_albedo_direct(:,:,nlev+1) = 0.0_jprb

      do jg = 1,ng
        ! Copy surface albedo in clear-sky region
        total_albedo(jg,1,nlev+1) = albedo_diffuse(jg,jcol)

        ! If direct albedo is available, use it; otherwise copy from
        ! diffuse albedo
        total_albedo_direct(jg,1,nlev+1) &
          &  = mu0 * albedo_direct(jg,jcol)
      end do

      ! If there is cloud in the lowest layer then we need the albedos
      ! underneath
      if (.not. is_clear_sky_layer(nlev)) then
        do jreg = 2,nregions
          total_albedo(:,jreg,nlev+1)        = total_albedo(:,1,nlev+1)
          total_albedo_direct(:,jreg,nlev+1) = total_albedo_direct(:,1,nlev+1)
        end do
      end if

      if (config%do_clear) then
        total_albedo_clear(:,nlev+1)        = total_albedo(:,1,nlev+1)
        total_albedo_clear_direct(:,nlev+1) = total_albedo_direct(:,1,nlev+1)
      end if

      ! Work up from the surface computing the total albedo of the
      ! atmosphere below that point using the adding method
      do jlev = nlev,1,-1

        total_albedo_below        = 0.0_jprb
        total_albedo_below_direct = 0.0_jprb

        if (config%do_clear) then
          ! For clear-skies there is no need to consider "above" and
          ! "below" quantities since with no cloud overlap to worry
          ! about, these are the same
          do jg = 1,ng
            inv_denom(jg,1) = 1.0_jprb &
                 &  / (1.0_jprb - total_albedo_clear(jg,jlev+1)*reflectance_clear(jg,jlev))
  
            total_albedo_clear(jg,jlev) = reflectance_clear(jg,jlev) &
                 &  + transmittance_clear(jg,jlev) * transmittance_clear(jg,jlev) &
                 &  * total_albedo_clear(jg,jlev+1) * inv_denom(jg,1)
  
            total_albedo_clear_direct(jg,jlev) = ref_dir_clear(jg,jlev) &
                 &  + (trans_dir_dir_clear(jg,jlev)*total_albedo_clear_direct(jg,jlev+1) &
                 &     +trans_dir_diff_clear(jg,jlev)*total_albedo_clear(jg,jlev+1)) &
                 &  * transmittance_clear(jg,jlev) * inv_denom(jg,1)
             end do
        end if

        do jg = 1,ng
          inv_denom(jg,1) = 1.0_jprb &
               &  / (1.0_jprb - total_albedo(jg,1,jlev+1)*reflectance_clear(jg,jlev))

          total_albedo_below(jg,1) = reflectance_clear(jg,jlev) &
               &  + transmittance_clear(jg,jlev)  * transmittance_clear(jg,jlev) &
               &  * total_albedo(jg,1,jlev+1) * inv_denom(jg,1)

          total_albedo_below_direct(jg,1) = ref_dir_clear(jg,jlev) &
               &  + (trans_dir_dir_clear(jg,jlev)*total_albedo_direct(jg,1,jlev+1) &
               &     +trans_dir_diff_clear(jg,jlev)*total_albedo(jg,1,jlev+1)) &
               &  * transmittance_clear(jg,jlev) * inv_denom(jg,1)
        end do

        if (.not. is_clear_sky_layer(jlev)) then
          inv_denom(:,2:) = 1.0_jprb / (1.0_jprb - total_albedo(:,2:,jlev+1)*reflectance(:,2:,jlev))
          total_albedo_below(:,2:) = reflectance(:,2:,jlev) &
               &  + transmittance(:,2:,jlev)  * transmittance(:,2:,jlev) &
               &  * total_albedo(:,2:,jlev+1) * inv_denom(:,2:)
          total_albedo_below_direct(:,2:) = ref_dir(:,2:,jlev) &
               &  + (trans_dir_dir(:,2:,jlev)*total_albedo_direct(:,2:,jlev+1) &
               &     +trans_dir_diff(:,2:,jlev)*total_albedo(:,2:,jlev+1)) &
               &  * transmittance(:,2:,jlev) * inv_denom(:,2:)
        end if

        ! Account for cloud overlap when converting albedo below a
        ! layer interface to the equivalent values just above
        if (is_clear_sky_layer(jlev) .and. is_clear_sky_layer(jlev-1)) then
          total_albedo(:,:,jlev)        = total_albedo_below(:,:)
          total_albedo_direct(:,:,jlev) = total_albedo_below_direct(:,:)
        else
          ! Use overlap matrix and exclude "anomalous" horizontal
          ! transport described by Shonk & Hogan (2008).  Therefore,
          ! the operation we perform is essentially diag(total_albedo)
          ! = matmul(transpose(v_matrix)), diag(total_albedo_below)).
          total_albedo(:,:,jlev) = 0.0_jprb
          total_albedo_direct(:,:,jlev) = 0.0_jprb
          do jreg = 1,nregions
            do jreg2 = 1,nregions
              total_albedo(:,jreg,jlev) &
                   &  = total_albedo(:,jreg,jlev) &
                   &  + total_albedo_below(:,jreg2) &
                   &  * v_matrix(jreg2,jreg,jlev)
              total_albedo_direct(:,jreg,jlev) &
                   &  = total_albedo_direct(:,jreg,jlev) &
                   &  + total_albedo_below_direct(:,jreg2) &
                   &  * v_matrix(jreg2,jreg,jlev)
            end do
          end do

        end if

      end do ! Reverse loop over levels

      ! --------------------------------------------------------
      ! Section 5: Compute fluxes
      ! --------------------------------------------------------

      ! Top-of-atmosphere fluxes into the regions of the top-most
      ! layer, zero since we assume no diffuse downwelling
      flux_dn = 0.0_jprb
      ! Direct downwelling flux (into a plane perpendicular to the
      ! sun) entering the top of each region in the top-most layer
      do jreg = 1,nregions
        direct_dn(:,jreg) = incoming_sw(:,jcol)*region_fracs(jreg,1,jcol)
      end do
      flux_up = direct_dn*total_albedo_direct(:,:,1)

      if (config%do_clear) then
        flux_dn_clear = 0.0_jprb
        direct_dn_clear(:) = incoming_sw(:,jcol)
        flux_up_clear = direct_dn_clear*total_albedo_clear_direct(:,1)
      end if

      ! Store the TOA broadband fluxes
      flux%sw_up(jcol,1) = sum(sum(flux_up,1))
      flux%sw_dn(jcol,1) = mu0 * sum(sum(direct_dn,1))
      if (allocated(flux%sw_dn_direct)) then
        flux%sw_dn_direct(jcol,1) = flux%sw_dn(jcol,1)
      end if
      if (config%do_clear) then
        flux%sw_up_clear(jcol,1) = sum(flux_up_clear)
        flux%sw_dn_clear(jcol,1) = mu0 * sum(direct_dn_clear)
        if (allocated(flux%sw_dn_direct_clear)) then
          flux%sw_dn_direct_clear(jcol,1) = flux%sw_dn_clear(jcol,1)
        end if
      end if

      ! Save the spectral fluxes if required; some redundancy here as
      ! the TOA downwelling flux is the same in clear and cloudy skies
      if (config%do_save_spectral_flux) then
        call indexed_sum(sum(flux_up,2), &
             &           config%i_spec_from_reordered_g_sw, &
             &           flux%sw_up_band(:,jcol,1))
        call indexed_sum(sum(direct_dn,2), &
             &           config%i_spec_from_reordered_g_sw, &
             &           flux%sw_dn_band(:,jcol,1))
        flux%sw_dn_band(:,jcol,1) = &
             &  mu0 * flux%sw_dn_band(:,jcol,1)
        if (allocated(flux%sw_dn_direct_band)) then
          flux%sw_dn_direct_band(:,jcol,1) = flux%sw_dn_band(:,jcol,1)
        end if
        call add_indexed_sum(sum(flux_dn,2), &
             &           config%i_spec_from_reordered_g_sw, &
             &           flux%sw_dn_band(:,jcol,1))
        if (config%do_clear) then
          call indexed_sum(flux_up_clear, &
               &           config%i_spec_from_reordered_g_sw, &
               &           flux%sw_up_clear_band(:,jcol,1))
          call indexed_sum(direct_dn_clear, &
               &           config%i_spec_from_reordered_g_sw, &
               &           flux%sw_dn_clear_band(:,jcol,1))
          flux%sw_dn_clear_band(:,jcol,1) &
               &  = mu0 * flux%sw_dn_clear_band(:,jcol,1)
          if (allocated(flux%sw_dn_direct_clear_band)) then
            flux%sw_dn_direct_clear_band(:,jcol,1) &
                 &  = flux%sw_dn_clear_band(:,jcol,1)
          end if
          call add_indexed_sum(flux_dn_clear, &
               &           config%i_spec_from_reordered_g_sw, &
               &           flux%sw_dn_clear_band(:,jcol,1))
        end if
      end if

      ! Final loop back down through the atmosphere to compute fluxes
      do jlev = 1,nlev
        if (config%do_clear) then
          do jg = 1,ng
            flux_dn_clear(jg) = (transmittance_clear(jg,jlev)*flux_dn_clear(jg) + direct_dn_clear(jg) &
               &  * (trans_dir_dir_clear(jg,jlev)*total_albedo_clear_direct(jg,jlev+1)*reflectance_clear(jg,jlev) &
               &     + trans_dir_diff_clear(jg,jlev) )) &
               &  / (1.0_jprb - reflectance_clear(jg,jlev)*total_albedo_clear(jg,jlev+1))
            direct_dn_clear(jg) = trans_dir_dir_clear(jg,jlev)*direct_dn_clear(jg)
            flux_up_clear(jg) = direct_dn_clear(jg)*total_albedo_clear_direct(jg,jlev+1) &
               &        +   flux_dn_clear(jg)*total_albedo_clear(jg,jlev+1)
          end do
        end if

        ! Fluxes for the clear sky region
        do jg = 1,ng
          flux_dn(jg,1) = (transmittance_clear(jg,jlev)*flux_dn(jg,1) + direct_dn(jg,1) &
               &  * (trans_dir_dir_clear(jg,jlev)*total_albedo_direct(jg,1,jlev+1)*reflectance_clear(jg,jlev) &
               &     + trans_dir_diff_clear(jg,jlev) )) &
               &  / (1.0_jprb - reflectance_clear(jg,jlev)*total_albedo(jg,1,jlev+1))
          direct_dn(jg,1) = trans_dir_dir_clear(jg,jlev)*direct_dn(jg,1)
          flux_up(jg,1) = direct_dn(jg,1)*total_albedo_direct(jg,1,jlev+1) &
               &  +        flux_dn(jg,1)*total_albedo(jg,1,jlev+1)
        end do
        ! Fluxes for cloudy regions, if they exist
        if (.not. is_clear_sky_layer(jlev)) then
          flux_dn(:,2:) = (transmittance(:,2:,jlev)*flux_dn(:,2:) + direct_dn(:,2:) &
               &  * (trans_dir_dir(:,2:,jlev)*total_albedo_direct(:,2:,jlev+1)*reflectance(:,2:,jlev) &
               &     + trans_dir_diff(:,2:,jlev) )) &
               &  / (1.0_jprb - reflectance(:,2:,jlev)*total_albedo(:,2:,jlev+1))
          direct_dn(:,2:) = trans_dir_dir(:,2:,jlev)*direct_dn(:,2:)
          flux_up(:,2:) = direct_dn(:,2:)*total_albedo_direct(:,2:,jlev+1) &
               &  +   flux_dn(:,2:)*total_albedo(:,2:,jlev+1)
        end if

        if (.not. (is_clear_sky_layer(jlev) &
             &    .and. is_clear_sky_layer(jlev+1))) then
          ! Moved from above:
          ! The zero initialization is quite slow and only necessary if either the current and above layer 
          ! aren't clear-sky because we can simply avoid adding zeroes in the broadband reduction
          if (is_clear_sky_layer(jlev)) then
            flux_dn(:,2:)  = 0.0_jprb
            flux_up(:,2:)  = 0.0_jprb
            direct_dn(:,2:)= 0.0_jprb
          end if 
          ! Account for overlap rules in translating fluxes just above
          ! a layer interface to the values just below
          flux_dn = singlemat_x_vec_sw(ng, &
              &  v_matrix(:,:,jlev+1), flux_dn)
          direct_dn = singlemat_x_vec_sw(ng, &
              &  v_matrix(:,:,jlev+1), direct_dn)
        end if ! Otherwise the fluxes in each region are the same so nothing to do

        ! Compute and store the broadband fluxes
        if (is_clear_sky_layer(jlev) .and. is_clear_sky_layer(jlev+1)) then
          sums_up = 0.0_jprb; sums_dn = 0.0_jprb; sums_dn_dir = 0.0_jprb
#ifdef __NEC__
          !NEC$ shortloop
#else
          !$omp simd reduction(+:sums_up, sums_dn, sums_dn_dir)
#endif
          do jg = 1, ng  
            sums_up = sums_up + flux_up(jg,1) 
            sums_dn = sums_dn + flux_dn(jg,1) 
            sums_dn_dir = sums_dn_dir + direct_dn(jg,1)
          end do
        else
          sums_up = 0.0_jprb; sums_dn = 0.0_jprb; sums_dn_dir = 0.0_jprb
#ifdef __NEC__
          !NEC$ shortloop
#else
          !$omp simd reduction(+:sums_up, sums_dn, sums_dn_dir)
#endif
          do jg = 1, ng  
            sums_up = sums_up + flux_up(jg,1) + flux_up(jg,2) + flux_up(jg,3)
            sums_dn = sums_dn + flux_dn(jg,1) + flux_dn(jg,2) + flux_dn(jg,3)
            sums_dn_dir = sums_dn_dir + direct_dn(jg,1) + direct_dn(jg,2) + direct_dn(jg,3)
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
          call indexed_sum(sum(flux_up,2), &
               &           config%i_spec_from_reordered_g_sw, &
               &           flux%sw_up_band(:,jcol,jlev+1))
          call indexed_sum(sum(direct_dn,2), &
               &           config%i_spec_from_reordered_g_sw, &
               &           flux%sw_dn_band(:,jcol,jlev+1))
          flux%sw_dn_band(:,jcol,jlev+1) = &
               &  mu0 * flux%sw_dn_band(:,jcol,jlev+1)
          if (allocated(flux%sw_dn_direct_band)) then
            flux%sw_dn_direct_band(:,jcol,jlev+1) &
                 &  = flux%sw_dn_band(:,jcol,jlev+1)
          end if
          call add_indexed_sum(sum(flux_dn,2), &
               &           config%i_spec_from_reordered_g_sw, &
               &           flux%sw_dn_band(:,jcol,jlev+1))
          if (config%do_clear) then
            call indexed_sum(flux_up_clear, &
                 &           config%i_spec_from_reordered_g_sw, &
                 &           flux%sw_up_clear_band(:,jcol,jlev+1))
            call indexed_sum(direct_dn_clear, &
                 &           config%i_spec_from_reordered_g_sw, &
                 &           flux%sw_dn_clear_band(:,jcol,jlev+1))
            flux%sw_dn_clear_band(:,jcol,jlev+1) = &
                 &  mu0 * flux%sw_dn_clear_band(:,jcol,jlev+1)
            if (allocated(flux%sw_dn_direct_clear_band)) then
              flux%sw_dn_direct_clear_band(:,jcol,jlev+1) &
                   &  = flux%sw_dn_clear_band(:,jcol,jlev+1)
            end if
            call add_indexed_sum(flux_dn_clear, &
                 &           config%i_spec_from_reordered_g_sw, &
                 &           flux%sw_dn_clear_band(:,jcol,jlev+1))
          end if
        end if
      end do ! Final loop over levels
      
      ! Store surface spectral fluxes, if required (after the end of
      ! the final loop over levels, the current values of these arrays
      ! will be the surface values)
      flux%sw_dn_diffuse_surf_g(:,jcol) = sum(flux_dn,2)
      flux%sw_dn_direct_surf_g(:,jcol)  = mu0 * sum(direct_dn,2)
      if (config%do_clear) then
        flux%sw_dn_diffuse_surf_clear_g(:,jcol) = flux_dn_clear
        flux%sw_dn_direct_surf_clear_g(:,jcol)  = mu0 * direct_dn_clear
      end if

    end do ! Loop over columns

    if (lhook) call dr_hook('radiation_tripleclouds_sw:solver_tripleclouds_sw',1,hook_handle)

  end subroutine solver_tripleclouds_sw

end module radiation_tripleclouds_sw
