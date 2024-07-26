! radiation_tripleclouds_lw.F90 - Longwave "Tripleclouds" solver
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
!   2017-04-28  R. Hogan  Receive emission/albedo rather than planck/emissivity
!   2017-04-22  R. Hogan  Store surface fluxes at all g-points
!   2017-10-23  R. Hogan  Renamed single-character variables
!   2018-10-08  R. Hogan  Call calc_region_properties
!   2020-09-18  R. Hogan  Replaced some array expressions with loops
!   2020-09-19  R. Hogan  Implement the cloud-only-scattering optimization
!   2022-09-01  P. Ukkonen  Optimizations for much better performance with ECCKD, including: 
!                           batching section 3 computations, faster kernels, ng can be defined at compile time

module radiation_tripleclouds_lw

  public

contains
  ! Small routine for scaling cloud optical depth in the cloudy
  ! regions
#include "radiation_optical_depth_scaling.h"

  !---------------------------------------------------------------------
  ! This module contains just one subroutine, the longwave
  ! "Tripleclouds" solver in which cloud inhomogeneity is treated by
  ! dividing each model level into three regions, one clear and two
  ! cloudy (with differing optical depth). This approach was described
  ! by Shonk and Hogan (2008).
  subroutine solver_tripleclouds_lw(ng_lw_in, nlev,istartcol,iendcol, &
       &  config, cloud, & 
       &  od, ssa, g, od_cloud, ssa_cloud, g_cloud, planck_hl, &
       &  emission, albedo, &
       &  flux)

    use parkind1, only           : jprb
    use yomhook,  only           : lhook, dr_hook, jphook

!    use radiation_io, only             : nulout
    use radiation_config, only         : config_type, IPdfShapeGamma
    use radiation_cloud, only          : cloud_type
    use radiation_regions, only        : calc_region_properties
    use radiation_overlap, only        : calc_overlap_matrices_nocol
    use radiation_flux, only           : flux_type, indexed_sum
    use radiation_matrix, only         : singlemat_x_vec, singlemat_x_vec_lw
    use radiation_two_stream, only     : calc_ref_trans_lw, &
         &                               calc_no_scattering_transmittance_lw
    use radiation_adding_ica_lw, only  : adding_ica_lw, calc_fluxes_no_scattering_lw
    use radiation_lw_derivatives, only : calc_lw_derivatives_region
    implicit none

! Allow size of inner dimension (number of g-points) to be known at compile time if NG_LW is defined
#ifdef NG_LW
  integer, parameter :: ng = NG_LW
#else
#define ng ng_lw_in
#endif

    ! Inputs
    integer, intent(in) :: ng_lw_in           ! number of g-points
    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type),        intent(in) :: config
    type(cloud_type),         intent(in) :: cloud

    ! Gas and aerosol optical depth of each layer at each longwave
    ! g-point
    real(jprb), intent(in), dimension(ng,nlev,istartcol:iendcol) :: od

    ! Gas and aerosol single-scattering albedo and asymmetry factor,
    ! only if longwave scattering by aerosols is to be represented
    real(jprb), intent(in), &
         &  dimension(config%n_g_lw_if_scattering,nlev,istartcol:iendcol) :: ssa, g

    ! Cloud and precipitation optical depth of each layer in each
    ! longwave band
    real(jprb), intent(in) :: od_cloud(config%n_bands_lw,nlev,istartcol:iendcol)

    ! Cloud and precipitation single-scattering albedo and asymmetry
    ! factor, only if longwave scattering by clouds is to be
    ! represented
    real(jprb), intent(in), dimension(config%n_bands_lw_if_scattering, &
         &                            nlev,istartcol:iendcol) :: ssa_cloud, g_cloud

    ! Planck function (emitted flux from a black body) at half levels
    ! and at the surface at each longwave g-point
    real(jprb), intent(in), dimension(ng,nlev+1,istartcol:iendcol) :: planck_hl

    ! Emission (Planck*emissivity) and albedo (1-emissivity) at the
    ! surface at each longwave g-point
    real(jprb), intent(in), dimension(ng, istartcol:iendcol) :: emission, albedo

    ! Local variables
    integer, parameter :: nregions = 3

    ! Optical depth, single scattering albedo and asymmetry factor in
    ! each g-point including gas, aerosol and clouds
    ! real(jprb), dimension(ng) :: od_total, ssa_total, g_total
    real(jprb), dimension(ng, 2:nregions) :: od_total, ssa_total, g_total, &
          &                                  planck_hl_bot, planck_hl_top

    ! Modified optical depth after Tripleclouds scaling to represent
    ! cloud inhomogeneity
    real(jprb), dimension(ng) :: od_cloud_new

    ! Output
    type(flux_type), intent(inout):: flux
    
    ! In a clear-sky layer this will be 1, otherwise equal to nregions
    ! integer :: nreg

    ! The area fractions of each region
    real(jprb) :: region_fracs(1:nregions,nlev,istartcol:iendcol)

    ! The scaling used for the optical depth in the cloudy regions
    real(jprb) :: od_scaling(2:nregions,nlev,istartcol:iendcol)

    ! Directional overlap matrices defined at all layer interfaces
    ! including top-of-atmosphere and the surface
    real(jprb), dimension(nregions,nregions,nlev+1) :: u_matrix, v_matrix

    ! Two-stream variables
    real(jprb), dimension(ng) :: gamma1, gamma2

    ! Diffuse reflection and transmission matrices of each layer
    ! real(jprb), dimension(ng, nregions, nlev) :: reflectance, transmittance

    ! Since we already have separate clear-sky variables, we only have to
    ! allocate the cloudy regions here, except for transmittance which needs
    ! all regions because it's used in in calc_lw_derivatives_region
    real(jprb), dimension(ng, 2:nregions, nlev) :: reflectance
    real(jprb), dimension(ng, nregions, nlev) ::  transmittance

    ! Emission by a layer into the upwelling or downwelling diffuse
    ! streams, cloudy regions
    real(jprb), dimension(ng, 2:nregions, nlev) &
         &  :: source_up, source_dn

    ! Clear-sky reflectance and transmittance
    real(jprb), dimension(ng, nlev) &
         &  :: ref_clear, trans_clear

    ! ...clear-sky equivalent emissions
    real(jprb), dimension(ng, nlev) &
         &  :: source_up_clear, source_dn_clear

    ! Total albedo of the atmosphere/surface just above a layer
    ! interface with respect to downwelling diffuse radiation at that
    ! interface, where level index = 1 corresponds to the
    ! top-of-atmosphere
    real(jprb), dimension(ng,nregions,nlev+1) :: total_albedo 
    ! Upwelling radiation just above a layer interface due to emission
    ! below that interface, where level index = 1 corresponds to the
    ! top-of-atmosphere
    real(jprb), dimension(ng,nregions,nlev+1) :: total_source
    ! Total albedo and source of the atmosphere just below a layer interface
    real(jprb), dimension(ng, nregions) &
         &  :: total_albedo_below, total_source_below

    ! Downwelling flux below and above an interface between
    ! layers into a plane perpendicular to the direction of the sun
    real(jprb), dimension(ng, nregions) &
         &  :: flux_dn, flux_up!, flux_dn_below

    ! Faster broadband flux computation 
    real(jprb) :: sums_up, sums_dn

    ! ...clear-sky equivalent (no distinction between "above/below")
    real(jprb), dimension(ng, nlev+1) &
         &  :: flux_dn_clear, flux_up_clear

    ! Clear-sky equivalent, but actually its reciprocal to replace
    ! some divisions by multiplications
    real(jprb), dimension(ng, nregions) :: inv_denom

    ! Identify clear-sky layers, with pseudo layers for outer space
    ! and below the ground, both treated as single-region clear skies
    logical :: is_clear_sky_layer(0:nlev+1)

    ! Index of the highest cloudy layer
    integer :: i_cloud_top

    integer :: jcol, jlev, jg, jreg, jreg2, ret

    real(jprb), parameter :: eps = epsilon(1.0_jprb)
    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_tripleclouds_lw:solver_tripleclouds_lw',0,hook_handle)

    ! --------------------------------------------------------
    ! Section 1: Prepare general variables and arrays
    ! --------------------------------------------------------
    ! Copy array dimensions to local variables for convenience
    ! ng   = ng_lw

    ! Compute the wavelength-independent region fractions and
    ! optical-depth scalings
    call calc_region_properties(nlev,nregions,istartcol,iendcol, &
         &  config%i_cloud_pdf_shape == IPdfShapeGamma, &
         &  cloud%fraction, cloud%fractional_std, region_fracs, &
         &  od_scaling, config%cloud_fraction_threshold)

    ! Main loop over columns
    do jcol = istartcol, iendcol

      ! --------------------------------------------------------
      ! Section 2: Prepare column-specific variables and arrays
      ! --------------------------------------------------------

      ! Define which layers contain cloud; assume that
      ! cloud%crop_cloud_fraction has already been called
      is_clear_sky_layer = .true.
      i_cloud_top = nlev+1
      do jlev = 1,nlev
        if (cloud%fraction(jcol,jlev) > 0.0_jprb) then
          is_clear_sky_layer(jlev) = .false.
          ! Get index to the first cloudy layer from the top
          if (i_cloud_top > jlev) then
            i_cloud_top = jlev
          end if
        end if
      end do
      if (config%do_lw_aerosol_scattering) then
        ! This is actually the first layer in which we need to
        ! consider scattering
        i_cloud_top = 1
      end if

      ! Compute wavelength-independent overlap matrices u_matrix and v_matrix 
      call calc_overlap_matrices_nocol(nlev,nregions, is_clear_sky_layer,&
      &  region_fracs(:,:,jcol), cloud%overlap_param(jcol,:), &
      &  v_matrix, u_matrix=u_matrix, decorrelation_scaling=config%cloud_inhom_decorr_scaling, &
      &  cloud_fraction_threshold=config%cloud_fraction_threshold, &
      &  use_beta_overlap=config%use_beta_overlap, &
      &  cloud_cover=flux%cloud_cover_lw(jcol))

      ! --------------------------------------------------------
      ! Section 3: Clear-sky calculation
      ! --------------------------------------------------------
      if (.not. config%do_lw_aerosol_scattering) then
        ! No scattering in clear-sky flux calculation; note that here
        ! the first two dimensions of the input arrays are unpacked
        ! into vectors inside the routine  
        call calc_no_scattering_transmittance_lw(ng*nlev, od(:,:,jcol), &
        &  planck_hl(:,1:nlev,jcol), planck_hl(:,2:nlev+1,jcol), &
        &  trans_clear, source_up_clear, source_dn_clear)

        ! ref_clear = 0.0_jprb
        ! Simple down-then-up method to compute fluxes
        call calc_fluxes_no_scattering_lw(ng, nlev, &
             &  trans_clear, source_up_clear, source_dn_clear, &
             &  emission(:,jcol), albedo(:,jcol), &
             &  flux_up_clear, flux_dn_clear)

      else
        ! Scattering in clear-sky flux calculation
        call calc_ref_trans_lw(ng*nlev, &
            &  od(:,:,jcol), ssa(:,:,jcol), g(:,:,jcol), &
            &  planck_hl(:,1:nlev,jcol), planck_hl(:,2:nlev+1,jcol), &
            &  ref_clear, trans_clear, &
            &  source_up_clear, source_dn_clear)
        ! Use adding method to compute fluxes
        call adding_ica_lw(ng, nlev, &
             &  ref_clear, trans_clear, source_up_clear, source_dn_clear, &
             &  emission(:,jcol), albedo(:,jcol), &
             &  flux_up_clear, flux_dn_clear)
      end if

      if (config%do_clear) then
        ! Sum over g-points to compute broadband fluxes
        do jlev = 1,nlev+1
          sums_up = 0.0_jprb; sums_dn = 0.0_jprb
#ifdef __NEC__
          !NEC$ shortloop
#else
          !$omp simd reduction(+:sums_up, sums_dn)
#endif
          do jg = 1, ng  
            sums_up = sums_up + flux_up_clear(jg,jlev)
            sums_dn = sums_dn + flux_dn_clear(jg,jlev)
          end do
          flux%lw_up_clear(jcol,jlev) = sums_up
          flux%lw_dn_clear(jcol,jlev) = sums_dn
       end do

        ! Store surface spectral downwelling fluxes / TOA upwelling
        do jg = 1,ng
          flux%lw_dn_surf_clear_g(jg,jcol) = flux_dn_clear(jg,nlev+1)
          flux%lw_up_toa_clear_g (jg,jcol) = flux_up_clear(jg,1)
        end do
        ! Save the spectral fluxes if required
        if (config%do_save_spectral_flux) then
          do jlev = 1,nlev+1
            call indexed_sum(flux_up_clear(:,jlev), &
                 &           config%i_spec_from_reordered_g_lw, &
                 &           flux%lw_up_clear_band(:,jcol,jlev))
            call indexed_sum(flux_dn_clear(:,jlev), &
                 &           config%i_spec_from_reordered_g_lw, &
                 &           flux%lw_dn_clear_band(:,jcol,jlev))
          end do
        end if
      end if

      ! --------------------------------------------------------
      ! Section 4: Loop over cloudy layers to compute reflectance and transmittance
      ! --------------------------------------------------------
      ! In this section the reflectance, transmittance and sources
      ! are computed for each layer
      ! G-points are batched (collapsed) with cloudy regions (ng*2) for bigger vector length (helps ECCKD)
      
      ! Firstly, ensure clear-sky transmittance is valid for whole
      ! depth of the atmosphere, because even above cloud it is used
      ! by the LW derivatives

      transmittance(:,1,:) = trans_clear(:,:)
      ! Dummy values in cloudy regions above cloud top
      if (i_cloud_top > 0) then
        transmittance(:,2:,1:min(i_cloud_top,nlev)) = 1.0_jprb
      end if

      ! Set values of cloudy-regions in clear-sky layers to zero
      do jlev = i_cloud_top,nlev ! Start at cloud top and work down
        if (is_clear_sky_layer(jlev)) then
          transmittance(:,2:,jlev) = 1.0_jprb
          ! Only transmittance needs zeroes in clear-sky regions, as it's used by
          ! lw_derivatives - cloudy regions are otherwise not used in clear-sky layers
          ! reflectance(:,:,jlev)   = 0.0_jprb
          ! source_up(:,:,jlev)     = 0.0_jprb
          ! source_dn(:,:,jlev)     = 0.0_jprb
        else  

          do jreg = 2,nregions
            ! For batched computations, planck functions need to expanded into regions
            planck_hl_top(:,jreg) = planck_hl(:,jlev,jcol) 
            planck_hl_bot(:,jreg) = planck_hl(:,jlev+1,jcol) 

            ! Cloudy sky
            ! Add scaled cloud optical depth to clear-sky value
            if (config%n_g_lw == config%n_bands_lw) then ! ensures vectorization for ECCKD model
              do jg = 1, ng
                od_cloud_new(jg) = od_cloud(jg,jlev,jcol) * od_scaling(jreg,jlev,jcol)
                od_total(jg,jreg) = od(jg,jlev,jcol) + od_cloud_new(jg)
              end do
            else 
              od_cloud_new = od_cloud(config%i_band_from_reordered_g_lw,jlev,jcol) &
                  &       * od_scaling(jreg,jlev,jcol)
              od_total(:,jreg) = od(:,jlev,jcol) + od_cloud_new
            end if
            if (config%do_lw_cloud_scattering) then
              if (config%do_lw_aerosol_scattering) then
                do jg = 1, ng
                  ssa_total(jg,jreg) = (ssa(jg,jlev,jcol)*od(jg,jlev,jcol) &
                        &     + ssa_cloud(config%i_band_from_reordered_g_lw(jg),jlev,jcol) &
                        &     *  od_cloud_new(jg)) & 
                        &     / max(od_total(jg,jreg), eps)
                  ! ^when od_total is zero, numerator is also zero, so eps can be set to anything
                  g_total = (g(jg,jlev,jcol)*ssa(jg,jlev,jcol)*od(jg,jlev,jcol) &
                        &     +   g_cloud(config%i_band_from_reordered_g_lw(jg),jlev,jcol) &
                        &     * ssa_cloud(config%i_band_from_reordered_g_lw(jg),jlev,jcol) &
                        &     *  od_cloud_new(jg)) &
                        &     / max(ssa_total(jg,jreg)*od_total(jg,jreg),eps)
                  ! ^same principle
                end do
              else
                if (config%n_g_lw == config%n_bands_lw) then ! ensures vectorization for ECCKD model
                  do jg = 1, ng
                    ssa_total(jg,jreg) = ssa_cloud(jg,jlev,jcol) &
                        &     * od_cloud_new(jg) / max(od_total(jg,jreg),eps)
                    g_total(jg,jreg) = g_cloud(jg,jlev,jcol) &
                        &     * ssa_cloud(jg,jlev,jcol) &
                        &     * od_cloud_new(jg) / max(ssa_total(jg,jreg)*od_total(jg,jreg),eps)
                  end do
                else
                  do jg = 1, ng
                    ssa_total(jg,jreg)  = ssa_cloud(config%i_band_from_reordered_g_lw(jg),jlev,jcol) &
                        &     * od_cloud_new(jg) / max(od_total(jg,jreg),eps)
                    g_total(jg,jreg) = g_cloud(config%i_band_from_reordered_g_lw(jg),jlev,jcol) &
                        &     * ssa_cloud(config%i_band_from_reordered_g_lw(jg),jlev,jcol) &
                        &     * od_cloud_new(jg) / max(ssa_total(jg,jreg)*od_total(jg,jreg),eps)
                  end do
                end if

              end if ! do_lw_aerosol_scattering
            end if ! do_lw_cloud_scattering
  
          end do

          if (config%do_lw_cloud_scattering) then
            call calc_ref_trans_lw(ng*2, &
            &  od_total(:,2:nregions), ssa_total(:,2:nregions), g_total(:,2:nregions), &
            &  planck_hl_top(:,2:nregions), planck_hl_bot(:,2:nregions), &
            &  reflectance(:,2:nregions,jlev), transmittance(:,2:nregions,jlev), &
            &  source_up(:,2:nregions,jlev), source_dn(:,2:nregions,jlev))
          else
            ! No-scattering case: use simpler functions for
            ! transmission and emission
            call calc_no_scattering_transmittance_lw(ng*2, &
                &  od_total(:,2:nregions), &
                &  planck_hl_top(:,2:nregions), planck_hl_bot(:,2:nregions), &
                &  transmittance(:,2:nregions,jlev), &
                &  source_up(:,2:nregions,jlev), source_dn(:,2:nregions,jlev))

            reflectance(:,2:nregions, jlev) = 0.0_jprb
          end if

          ! Emission is scaled by the size of each region
          do jreg = 2,nregions
            source_up(:,jreg,jlev) = region_fracs(jreg,jlev,jcol) * source_up(:,jreg,jlev)
            source_dn(:,jreg,jlev) = region_fracs(jreg,jlev,jcol) * source_dn(:,jreg,jlev)
          end do
          source_up_clear(:,jlev) = region_fracs(1,jlev,jcol) * source_up_clear(:,jlev)
          source_dn_clear(:,jlev) = region_fracs(1,jlev,jcol) * source_dn_clear(:,jlev)
        
        end if
      end do

      ! --------------------------------------------------------
      ! Section 5: Compute total sources and albedos at each half level
      ! --------------------------------------------------------

      total_albedo(:,:,i_cloud_top:nlev+1) = 0.0_jprb

      ! Calculate the upwelling radiation emitted by the surface, and
      ! copy the surface albedo into total_albedo 
      do jreg = 1,nregions
        do jg = 1,ng
          ! region_fracs(jreg,nlev,jcol) is the fraction of each region in the
          ! lowest model level
          total_source(jg,jreg,nlev+1) = region_fracs(jreg,nlev,jcol)*emission(jg,jcol)
          total_albedo(jg,jreg,nlev+1) = albedo(jg,jcol)
        end do
      end do

      ! Work up from the surface computing the total albedo of the
      ! atmosphere and the total upwelling due to emission below each
      ! level below using the adding method
      do jlev = nlev,i_cloud_top,-1
        total_albedo_below        = 0.0_jprb
        total_source_below        = 0.0_jprb

        ! Clear-sky region
        if (.not. config%do_lw_aerosol_scattering) then ! ref_clear is zero (inv_denom is one)
          total_albedo_below(:,1) = transmittance(:,1,jlev)*transmittance(:,1,jlev)*total_albedo(:,1,jlev+1)
          total_source_below(:,1) = source_up_clear(:,jlev) &
                &  + transmittance(:,1,jlev)*(total_source(:,1,jlev+1) &
                &  + total_albedo(:,1,jlev+1)*source_dn_clear(:,jlev))
        else 
          inv_denom(:,1) = 1.0_jprb &
                &  / (1.0_jprb - total_albedo(:,1,jlev+1)*ref_clear(:,jlev))
          total_albedo_below(:,1) = ref_clear(:,jlev) &
                &  + transmittance(:,1,jlev)*transmittance(:,1,jlev)*total_albedo(:,1,jlev+1) &
                &  * inv_denom(:,1)
          total_source_below(:,1) = source_up_clear(:,jlev) &
                &  + transmittance(:,1,jlev)*(total_source(:,1,jlev+1) &
                &  + total_albedo(:,1,jlev+1)*source_dn_clear(:,jlev)) &
                &  * inv_denom(:,1)
        end if 

        if (.not. is_clear_sky_layer(jlev)) then
          inv_denom(:,2:) = 1.0_jprb / (1.0_jprb - total_albedo(:,2:,jlev+1)*reflectance(:,2:,jlev))
          total_albedo_below(:,2:) = reflectance(:,2:,jlev) &
               &  + transmittance(:,2:,jlev)*transmittance(:,2:,jlev)*total_albedo(:,2:,jlev+1) &
               &  * inv_denom(:,2:)
          total_source_below(:,2:) = source_up(:,2:,jlev) &
               &  + transmittance(:,2:,jlev)*(total_source(:,2:,jlev+1) &
               &  + total_albedo(:,2:,jlev+1)*source_dn(:,2:,jlev)) &
               &  * inv_denom(:,2:)
        end if

        ! Account for cloud overlap when converting albedo below a
        ! layer interface to the equivalent values just above
        if (is_clear_sky_layer(jlev) .and. is_clear_sky_layer(jlev-1)) then
          total_albedo(:,:,jlev) = total_albedo_below(:,:)
          total_source(:,:,jlev) = total_source_below(:,:)
        else
          total_source(:,:,jlev) = singlemat_x_vec_lw(ng,&
                &  u_matrix(:,:,jlev), total_source_below)  
          ! Use overlap matrix and exclude "anomalous" horizontal
          ! transport described by Shonk & Hogan (2008).  Therefore,
          ! the operation we perform is essentially diag(total_albedo)
          ! = matmul(transpose(v_matrix), diag(total_albedo_below)).
          do jreg = 1,nregions
            do jreg2 = 1,nregions
              total_albedo(:,jreg,jlev) &
                   &  = total_albedo(:,jreg,jlev) &
                   &  + total_albedo_below(:,jreg2) &
                   &  * v_matrix(jreg2,jreg,jlev)
            end do
          end do
        end if
        
      end do ! Reverse loop over levels

      ! Section 6: Compute fluxes up to top-of-atmosphere
      ! --------------------------------------------------------
      ! Compute the fluxes just above the highest cloud

      flux_up(:,1) = total_source(:,1,i_cloud_top) &
           &  + total_albedo(:,1,i_cloud_top)*flux_dn_clear(:,i_cloud_top)
      flux_up(:,2:) = 0.0_jprb
      flux%lw_up(jcol,i_cloud_top) = sum(flux_up(:,1))
      if (config%do_save_spectral_flux) then
        call indexed_sum(flux_up(:,1), &
             &           config%i_spec_from_reordered_g_lw, &
             &           flux%lw_up_band(:,jcol,i_cloud_top))
      end if
      ! do jlev = i_cloud_top-1,1,-1
      do jlev = i_cloud_top,1,-1
        if (jlev /= i_cloud_top) then
          sums_up = 0.0_jprb
#ifdef __NEC__
          !NEC$ shortloop
#else
          !$omp simd reduction(+:sums_up)
#endif
          do jg = 1, ng
            flux_up(jg,1) = trans_clear(jg,jlev)*flux_up(jg,1) + source_up_clear(jg,jlev)
            sums_up       = sums_up + flux_up(jg,1)
          end do 
          flux%lw_up(jcol,jlev) = sums_up
        end if 
        if (config%do_save_spectral_flux) then
          if (jlev /= i_cloud_top) then
            call indexed_sum(flux_up(:,1), &
                &           config%i_spec_from_reordered_g_lw, &
                &           flux%lw_up_band(:,jcol,jlev))
          end if 
    ! --------------------------------------------------------
    ! Section 7: Copy over downwelling fluxes above cloud top
    ! --------------------------------------------------------
          if (config%do_clear) then
            ! Clear-sky fluxes have already been averaged: use these
            flux%lw_dn_band(:,jcol,jlev) = flux%lw_dn_clear_band(:,jcol,jlev)
          else 
            call indexed_sum(flux_dn_clear(:,jlev), &
                 &           config%i_spec_from_reordered_g_lw, &
                 &           flux%lw_dn_band(:,jcol,jlev))
          end if 
        end if
        ! Down-welling flux
        if (config%do_clear) then
          ! Clear-sky fluxes have already been averaged: use these
          flux%lw_dn(jcol,jlev) = flux%lw_dn_clear(jcol,jlev)
        else
          flux%lw_dn(jcol,:) = sum(flux_dn_clear(:,jlev))
        end if
      end do
      flux%lw_up_toa_g(:,jcol) = sum(flux_up,2)

      ! --------------------------------------------------------
      ! Section 8: Compute fluxes down to surface
      ! --------------------------------------------------------

      ! Copy over downwelling spectral fluxes at top of first
      ! scattering layer, using overlap matrix to translate to the
      ! regions of the first layer of cloud
      do jreg = 1,nregions
        flux_dn(:,jreg)  = v_matrix(jreg,1,i_cloud_top)*flux_dn_clear(:,i_cloud_top)
      end do

      ! Final loop back down through the atmosphere to compute fluxes
      do jlev = i_cloud_top,nlev

        ! Clear-sky region
        if (.not. config%do_lw_aerosol_scattering) then ! ref_clear is zero
          flux_dn(:,1) = (transmittance(:,1,jlev)*flux_dn(:,1) + source_dn_clear(:,jlev) )
          flux_up(:,1) = total_source(:,1,jlev+1) + flux_dn(:,1)*total_albedo(:,1,jlev+1)
        else 
          flux_dn(:,1) = (transmittance(:,1,jlev)*flux_dn(:,1) &
                &  + ref_clear(:,jlev)*total_source(:,1,jlev+1) + source_dn_clear(:,jlev) ) &
                &  / (1.0_jprb - ref_clear(:,jlev)*total_albedo(:,1,jlev+1))
          flux_up(:,1) = total_source(:,1,jlev+1) + flux_dn(:,1)*total_albedo(:,1,jlev+1)
        end if 

        if (is_clear_sky_layer(jlev)) then
          ! flux_dn(:,2:)  = 0.0_jprb
          ! flux_up(:,2:)  = 0.0_jprb
        else
          flux_dn(:,2:) = (transmittance(:,2:,jlev)*flux_dn(:,2:) &
               &     + reflectance(:,2:,jlev)*total_source(:,2:,jlev+1) + source_dn(:,2:,jlev) ) &
               &  / (1.0_jprb - reflectance(:,2:,jlev)*total_albedo(:,2:,jlev+1))
          flux_up(:,2:) = total_source(:,2:,jlev+1) + flux_dn(:,2:)*total_albedo(:,2:,jlev+1)
        end if

        if (.not. (is_clear_sky_layer(jlev) &
             &    .and. is_clear_sky_layer(jlev+1))) then
            ! moved from above
            if (is_clear_sky_layer(jlev)) then
              flux_up(:,2:)  = 0.0_jprb
              flux_dn(:,2:)  = 0.0_jprb
            end if
            ! Account for overlap rules in translating fluxes just above
            ! a layer interface to the values just below
            flux_dn = singlemat_x_vec_lw(ng, v_matrix(:,:,jlev+1), flux_dn) 
        end if ! Otherwise the fluxes in each region are the same so nothing to do

        ! Store the broadband fluxes
        sums_up = 0.0_jprb; sums_dn = 0.0_jprb
        if (is_clear_sky_layer(jlev) .and. is_clear_sky_layer(jlev+1)) then
#ifdef __NEC__
          !NEC$ shortloop
#else
          !$omp simd reduction(+:sums_up, sums_dn)
#endif
          do jg = 1, ng  
            sums_up = sums_up + flux_up(jg,1)
            sums_dn = sums_dn + flux_dn(jg,1)
          end do
        else
#ifdef __NEC__
          !NEC$ shortloop
#else
          !$omp simd reduction(+:sums_up, sums_dn)
#endif
          do jg = 1, ng  
            sums_up = sums_up + flux_up(jg,1) + flux_up(jg,2) + flux_up(jg,3)
            sums_dn = sums_dn + flux_dn(jg,1) + flux_dn(jg,2) + flux_dn(jg,3)
          end do
        end if
        flux%lw_up(jcol,jlev+1) = sums_up
        flux%lw_dn(jcol,jlev+1) = sums_dn

        ! Save the spectral fluxes if required
        if (config%do_save_spectral_flux) then
          call indexed_sum(sum(flux_up,2), &
               &           config%i_spec_from_reordered_g_lw, &
               &           flux%lw_up_band(:,jcol,jlev+1))
          call indexed_sum(sum(flux_dn,2), &
               &           config%i_spec_from_reordered_g_lw, &
               &           flux%lw_dn_band(:,jcol,jlev+1))
         end if
      end do ! Final loop over levels
      
      ! Store surface spectral downwelling fluxes, which at this point
      ! are at the surface
      flux%lw_dn_surf_g(:,jcol) = sum(flux_dn,2)

      ! Compute the longwave derivatives needed by Hogan and Bozzo
      ! (2015) approximate radiation update scheme
      if (config%do_lw_derivatives) then
        ! Note that at this point flux_up contains the spectral
        ! fluxes into the regions of the lowest layer; we sum over
        ! regions first to provide a simple spectral flux upwelling
        ! from the surface
        call calc_lw_derivatives_region(ng, nlev, nregions, jcol, transmittance, &
             &  u_matrix(:,:,:), sum(flux_up,2), flux%lw_derivatives)
      end if

    end do ! Loop over columns

    if (lhook) call dr_hook('radiation_tripleclouds_lw:solver_tripleclouds_lw',1,hook_handle)

  end subroutine solver_tripleclouds_lw

end module radiation_tripleclouds_lw
