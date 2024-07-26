! radiation_spartacus_lw.F90 - SPARTACUS longwave solver
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
!   2017-04-11  R. Hogan  Receive emission/albedo rather than planck/emissivity
!   2017-04-22  R. Hogan  Store surface fluxes at all g-points
!   2018-09-03  R. Hogan  Security via min_cloud_effective_size
!   2018-10-08  R. Hogan  Call calc_region_properties
!   2019-01-12  R. Hogan  Use inv_inhom_effective_size if allocated
!   2022-09-01  P. Ukkonen  Optimizations for much better performance with ECCKD, including: 
!                           batching section 3 computations, faster kernels, ng can be defined at compile time
module radiation_spartacus_lw

  public

! Allow size of inner dimension (number of g-points) to be known at compile time if NG_LW is defined
#ifdef NG_LW
  integer, parameter, private :: ng = NG_LW
#else
#define ng ng_lw_in
#endif

contains

  ! Small routine for scaling cloud optical depth in the cloudy
  ! regions
#include "radiation_optical_depth_scaling.h"

  ! This module contains just one subroutine, the longwave solver
  ! using the Speedy Algorithm for Radiative Transfer through Cloud
  ! Sides (SPARTACUS), which can represent 3D effects using a matrix
  ! form of the two-stream equations.
  !
  ! Sections:
  !   1: Prepare general variables and arrays
  !   2: Prepare column-specific variables and arrays
  !   3: First loop over layers
  !     3.1: Layer-specific initialization
  !     3.2: Compute gamma variables
  !       3.2a: Clear-sky case
  !       3.2b: Cloudy case
  !     3.3: Compute reflection, transmission and emission
  !       3.3a: g-points with 3D effects
  !       3.3b: g-points without 3D effects
  !   4: Compute total sources and albedos
  !   5: Compute fluxes
  subroutine solver_spartacus_lw(ng_lw_in, nlev,istartcol,iendcol, &
       &  config, thermodynamics, cloud, & 
       &  od, ssa, g, od_cloud, ssa_cloud, g_cloud, planck_hl, &
       &  emission, albedo, &
       &  flux)

    use parkind1,                 only : jprb
    use yomhook,                  only : lhook, dr_hook, jphook

    use radiation_io,             only : nulout
    use radiation_config,         only : config_type, IPdfShapeGamma
    use radiation_thermodynamics, only : thermodynamics_type
    use radiation_cloud,          only : cloud_type
    use radiation_regions,        only : calc_region_properties
    use radiation_overlap,        only : calc_overlap_matrices_nocol
    use radiation_flux,           only : flux_type, indexed_sum
    use radiation_matrix
    use radiation_two_stream,     only : calc_ref_trans_lw, & 
        & calc_no_scattering_transmittance_lw, LwDiffusivityWP, calc_two_stream_gammas_lw
    use radiation_lw_derivatives, only : calc_lw_derivatives_matrix
    use radiation_constants,      only : Pi, GasConstantDryAir, AccelDueToGravity
    implicit none

    ! Inputs
    integer, intent(in) :: ng_lw_in           ! number of g-points
    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type), intent(in)        :: config
    type(thermodynamics_type),intent(in) :: thermodynamics
    type(cloud_type),   intent(in)       :: cloud

    ! Gas and aerosol optical depth of each layer at each longwave
    ! g-point
    real(jprb), intent(in), dimension(ng,nlev,istartcol:iendcol) :: od

    ! Gas and aerosol single-scattering albedo and asymmetry factor,
    ! only if longwave scattering by aerosols is to be represented
    real(jprb), intent(in), target, &
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
    real(jprb), intent(in), dimension(ng,nlev+1,istartcol:iendcol) &
         &  :: planck_hl

    ! Emission (Planck*emissivity) and albedo (1-emissivity) at the
    ! surface at each longwave g-point
    real(jprb), intent(in), dimension(ng, istartcol:iendcol) &
         &  :: emission, albedo

    ! Output
    type(flux_type),          intent(inout):: flux

    ! Local variables
    integer, parameter :: nreg = 3
    integer :: jcol, jlev, jg, jj, jreg, jband, jreg2
    integer :: ng3D, ng3D_old! Number of g-points with small enough gas optical
                    ! depth that 3D effects need to be represented

    ! Ratio of gas constant for dry air to acceleration due to gravity
    real(jprb), parameter &
         &  :: R_over_g = GasConstantDryAir / AccelDueToGravity

    ! Used in computing rates of lateral radiation transfer
    real(jprb), parameter :: four_over_pi = 4.0_jprb / Pi

    ! Used in computing rates of lateral radiation transfersolver_spartacus_lw
    real(jprb), parameter :: tan_diffuse_angle_3d = Pi * 0.5_jprb
    ! Incorrect value of tand(53) used by Hogan & Shonk (2013)
    !    real(jprb), parameter :: tan_diffuse_angle_3d = 1.32704482162041

    ! Optical depth, single scattering albedo and asymmetry factor in
    ! each region at each g-point
    ! real(jprb), dimension(ng, 2:nreg,nlev) &
    !      &  :: od_region_cld, ssa_region_cld, g_region_cld
    real(jprb), dimension(ng,nlev) :: od_region_clear

    real(jprb), dimension(:,:), contiguous, pointer   :: ssa_clear, g_clear
    real(jprb), dimension(ng,nlev), target :: zero_array

    ! Scattering optical depths of gases and clouds
    real(jprb) :: scat_od(ng)
    real(jprb) :: scat_od_cloud

    ! The area fractions of each region
    real(jprb) :: region_fracs(1:nreg,nlev,istartcol:iendcol)

    ! The scaling used for the optical depth in the cloudy regions
    real(jprb) :: od_scaling(2:nreg,nlev,istartcol:iendcol)

    ! The length of the interface between regions per unit area of
    ! gridbox, equal to L_diff^ab in Hogan and Shonk (2013). This is
    ! actually the effective length oriented to a photon with random
    ! azimuth angle. The three values are the edge length between
    ! regions a-b, b-c and a-c.
    real(jprb) :: edge_length(3,nlev)

    ! Element i,j gives the rate of 3D transfer of diffuse
    ! radiation from region i to region j, multiplied by the thickness
    ! of the layer in m
    ! real(jprb) :: transfer_rate(nlev,nreg,nreg)

    ! Directional overlap matrices defined at all layer interfaces
    ! including top-of-atmosphere and the surface
    real(jprb), dimension(nreg,nreg,nlev+1) :: u_matrix, v_matrix

    ! Two-stream variables
    ! real(jprb), dimension(ng, nreg-1, nlev) :: gamma1_cld, gamma2_cld 
    real(jprb), dimension(ng, nlev) :: gamma1_clear, gamma2_clear

    ! Matrix Gamma multiplied by the layer thickness z1, so units
    ! metres.  After calling expm, this array contains the matrix
    ! exponential of the original.
    real(jprb), allocatable :: Gamma_z1(:,:,:)

    ! Diffuse reflection and transmission matrices of each layer
    real(jprb), dimension(ng, nreg, nreg, nlev) :: reflectance, transmittance

    ! Clear-sky diffuse reflection and transmission matrices of each
    ! layer
    real(jprb), dimension(ng, nlev) :: ref_clear, trans_clear

    ! If the Planck term at the top of a layer in each region is the
    ! vector b0 then the following is vector [-b0; b0]*dz
    real(jprb), dimension(:,:), allocatable :: planck_top

    ! The difference between the Planck terms at the bottom and top of
    ! a layer, doubled as with planck_top; in terms of b' in the
    ! paper, planck_diff is [-b'; b']*dz*dz
    real(jprb), dimension(:,:), allocatable :: planck_diff

    ! Parts of the particular solution associated with the
    ! inhomogeneous ODE. In terms of quantities from the paper,
    ! solution0 is [c0;d0] and solution_diff is [c';d']*dz.
    ! real(jprb), dimension(:,:), allocatable :: solution0, solution_diff

    ! Used for computing the Planck emission per layer
    real(jprb), dimension(:,:), allocatable :: tmp_vectors

    ! The fluxes upwelling from the top of a layer (source_up) and
    ! downwelling from the bottom of the layer (source_dn) due to
    ! emission
    real(jprb), dimension(ng, nreg, nlev)   :: source_up
    real(jprb), dimension(ng, nreg, nlev)   :: source_dn
    ! ...clear-sky equivalents
    real(jprb), dimension(ng, nlev)         :: source_up_clear, source_dn_clear

    ! Upwelling radiation just above a layer interface due to emission
    ! below that interface, where level index = 1 corresponds to the
    ! top-of-atmosphere
    real(jprb), dimension(ng, nreg, nlev+1) :: total_source
    ! ...clear-sky equivalent
    real(jprb), dimension(ng, nlev+1) :: total_source_clear

    ! As total_source, but just below a layer interface
    real(jprb), dimension(ng, nreg)   :: total_source_below

    ! Total albedo of the atmosphere/surface just above a layer
    ! interface with respect to downwelling diffuse radiation at that
    ! interface, where level index = 1 corresponds to the
    ! top-of-atmosphere
    real(jprb), dimension(ng, nreg, nreg, nlev+1) :: total_albedo
    ! ...clear-sky equivalent
    real(jprb), dimension(ng, nlev+1) :: total_albedo_clear

    ! As total_albedo, but just below a layer interface
    real(jprb), dimension(ng, nreg, nreg) :: total_albedo_below
    ! Temporary variables for albedo computations
    real(jprb), dimension(ng, nreg, nreg) :: X
    real(jprb), dimension(ng, nreg)       :: b

    ! The following is used to store matrices of the form I-A*B that
    ! are used on the denominator of some expressions
    real(jprb), dimension(ng, nreg, nreg) :: denominator
    ! Clear-sky equivalent, but actually its reciprocal to replace
    ! some divisions by multiplications
    real(jprb), dimension(ng) :: inv_denom_scalar

    ! Layer depth (m)
    real(jprb) :: dz

    ! Upwelling and downwelling fluxes above/below layer interfaces
    real(jprb), dimension(ng, nreg) &
         &  :: flux_up_above, flux_dn_above, flux_dn_below
    ! Clear-sky upwelling and downwelling fluxes (which don't
    ! distinguish between whether they are above/below a layer
    ! interface)
    real(jprb), dimension(ng) :: flux_up_clear, flux_dn_clear

    ! Faster broadband flux computation 
    real(jprb) :: sums_up, sums_dn

    ! Parameterization of internal flux distribution in a cloud that
    ! can lead to the flux just about to exit a cloud being different
    ! from the mean in-cloud value.  "lateral_od" is the typical
    ! absorption optical depth through the clormse of 7.34W/m2 for the boundary fluxes and 0.175K/d for the heating rate.ud in a horizontal
    ! direction.
    real(jprb) :: lateral_od, sqrt_1_minus_ssa, side_emiss_thick
    real(jprb), dimension(ng) ::  side_emiss

    ! In the optically thin limit, the flux near the edge is greater
    ! than that in the interior
    real(jprb), parameter :: side_emiss_thin = 1.4107

    real(jprb) :: aspect_ratio

    ! Keep a count of the number of calls to the two ways of computing
    ! reflectance/transmittance matrices
    integer :: n_calls_expm, n_calls_meador_weaver

    ! Identify clear-sky layers, with pseudo layers for outer space
    ! and below the ground, both treated as single-region clear skies
    logical :: is_clear_sky_layer(0:nlev+1)

    ! Layer depth (m)
    real(jprb) :: layer_depth(nlev)

    ! CLOUDY BATCHING
    logical :: is_cloudy_layer(1:nlev), are_clouds_below
    real(jprb), dimension(:,:,:), allocatable ::  gamma1_cld, gamma2_cld, &
        & od_region_cld, ssa_region_cld, g_region_cld, LU, transfer_rate, &
        & gamma_z11, gamma_z12, gamma_z21, gamma_z22
    real(jprb), dimension(:,:), allocatable :: &
       & source_up_tmp, source_dn_tmp
    integer, dimension(:,:), allocatable :: inds
    integer :: jtop, jbot, nlev_cld, j1, j2, ng3D_tot, ng_batch_limit

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_spartacus_lw:solver_spartacus_lw',0,hook_handle)

    ! --------------------------------------------------------
    ! Section 1: Prepare general variables and arrays
    ! --------------------------------------------------------

    ! For better performance the expm computations in section 3 are batched across several 
    ! cloudy layers. It may be useful to cap how many g-points (across adjacent layers) 
    ! are batched - too large and cache use may suffer
    ng_batch_limit = max(nint(4*300.0_jprb / (sizeof(1.0_jprb))),1) ! 300 for SP, 150 for DP
    ! for optimal performance user can hand-tune ng_batch_limit but for LW and ecCKD 
    ! ng_batch_limit had little impact on performance on ECMWF machine 

    ! Copy array dimensions to local variables for convenience
    ! nreg = config%nregions

    ! Reset count of number of calls to the two ways to compute
    ! reflectance/transmission matrices
    n_calls_expm = 0
    n_calls_meador_weaver = 0

    ! Initialize dz to avoid compiler warning
    dz = 1.0_jprb

    zero_array = 0.0_jprb
    ng3D_old = 0

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

      ! Define which layers contain cloud; assume that
      ! cloud%crop_cloud_fraction has already been called
      is_clear_sky_layer = .true.
      is_cloudy_layer = .false.
      do jlev = nlev,1,-1
        if (cloud%fraction(jcol,jlev) > 0.0_jprb) then
          is_clear_sky_layer(jlev) = .false.
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

            if (.not. is_clear_sky_layer(jlev)) then
                is_cloudy_layer(jlev) = .true. 
              ! Cloudy in 3D sense, fulfilling all above criteria
            end if

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


      ! Compute wavelength-independent overlap matrices u_matrix and v_matrix
      call calc_overlap_matrices_nocol(nlev,nreg, is_clear_sky_layer, &
      &  region_fracs(:,:,jcol), cloud%overlap_param(jcol,:), &
      &  v_matrix, u_matrix=u_matrix, decorrelation_scaling=config%cloud_inhom_decorr_scaling, &
      &  cloud_fraction_threshold=config%cloud_fraction_threshold, &
      &  use_beta_overlap=config%use_beta_overlap, &
      &  cloud_cover=flux%cloud_cover_lw(jcol))

      ! --------------------------------------------------------
      ! Section 3
      ! --------------------------------------------------------
      ! In this section the reflectance, transmittance and sources
      ! are computed for each layer

      ! Set optical depth of clear-sky region (region 1) to the
      ! gas/aerosol optical depth
      od_region_clear = od(:,:,jcol)

      if (config%do_lw_aerosol_scattering) then
        ssa_clear => ssa(:,:,jcol)
        g_clear   => g(:,:,jcol)
      else
        ssa_clear => zero_array
        g_clear   => zero_array
      end if

      ! Compute reflectance, transmittance and upward and downward
      ! sources for clear skies, using the Meador-Weaver formulas
      ! for reflectance and transmittance, and equivalent solutions
      ! to the coupled ODEs for the sources.

      if (config%do_lw_aerosol_scattering) then
        call calc_ref_trans_lw(ng*nlev, &
            &  od_region_clear, ssa_clear, g_clear, &
            &  planck_hl(:,1:nlev,jcol), planck_hl(:,2:nlev+1,jcol), &
            &  ref_clear, trans_clear, &
            &  source_up_clear, source_dn_clear, &
            &  gamma1_clear, gamma2_clear)
      else

        call calc_no_scattering_transmittance_lw(ng*nlev, &
        &  od(:,:, jcol), &
        &  planck_hl(:,1:nlev,jcol), planck_hl(:,2:nlev+1,jcol), & 
        &  trans_clear, &
        &  source_up_clear, source_dn_clear)

        gamma1_clear  = 1.65999997_jprb
        gamma2_clear  = 0.0_jprb
        ref_clear     = 0.0_jprb
      end if

      ! Initialization 
      reflectance   = 0.0_jprb
      transmittance = 0.0_jprb
      source_up     = 0.0_jprb
      source_dn     = 0.0_jprb

      ! only the non-3D g-points of the clear-sky region should be set to clear-sky values, 
      ! but we will later overwrite the 3D g-points
      reflectance   (:,1,1,:) = ref_clear
      transmittance (:,1,1,:) = trans_clear

      do jlev = 1, nlev 
        source_up(:,1,jlev) = region_fracs(1,jlev,jcol)*source_up_clear(:,jlev)
        source_dn(:,1,jlev) = region_fracs(1,jlev,jcol)*source_dn_clear(:,jlev)
      end do


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
          ! do jlev = jtop,nlev
            if (.not. is_cloudy_layer(jlev)) then
              jbot = jlev-1 ! bottom found, exit loop
              exit
            else if (jlev==nlev) then
              jbot = nlev
            end if
          end do
        end if

        ! nlev_cld = jbot - jtop + 1
        allocate(inds(2,jtop:jbot))
        
        ng3D = ng
        ng3D_tot = 0
        j2 = 0
        ! The 3D g-points across adjacent cloudy layers are gathered into one array like this:
        ! Each layer has a varying amount of 3D g-points, ng3D, but because the g-points have
        ! been ordered (roughly) by optical depth, it's always the first 1:ng3D g-points in a given 
        ! layer that 3D computations are performed for 
        ! Then we just need to keep track of the layer-wise index pairs that give the start and end 
        ! indices of a given layer's 3D g-points in the "gather array"
        ! Find indices
        do jlev = jtop, jbot
          j1 = j2 + 1
          do jg = 1,ng
            if (od_region_clear(jg,jlev) > config%max_gas_od_3D) then
              ng3D = jg-1
              exit
            end if
          end do
          j2 = j1 + ng3D -1
          inds(1,jlev) = j1 
          inds(2,jlev) = j2
          ng3D_tot = ng3D_tot + ng3D
          ! print *, "jlev ", jlev, "inds", inds(:,jlev), "ng3D for this lev", ng3D
          if (ng3D_tot > ng_batch_limit) then
            ! print *, ng3D_tot, ".. gpoints exceeding batch lim, setting jbot to ", jlev, "previously", jbot
            jbot = jlev
            exit
          end if
        end do

        nlev_cld = jbot - jtop + 1

        allocate(od_region_cld(ng,2:nreg,jtop:jbot), ssa_region_cld(ng,2:nreg,jtop:jbot))
        allocate(g_region_cld (ng,2:nreg,jtop:jbot))
        allocate(gamma1_cld(ng,2:nreg, jtop:jbot), gamma2_cld(ng,2:nreg, jtop:jbot))

        ! od_region_cld = 0.0_jprb
        ssa_region_cld = 0.0_jprb ! No Rayleigh scattering in longwave by default
        g_region_cld = 0.0_jprb

        ! Loop over cloudy levels
        do jlev = jtop, jbot
          if (config%do_lw_cloud_scattering) then
            ! Scattering optical depth of clear-sky region
            scat_od = od_region_clear(:,jlev)*ssa_clear(:,jlev)
          end if
          ! Loop over cloudy regions
          do jreg = 2,nreg
            ! Loop over g-points
            do jg = 1,ng
              if (ng == config%n_bands_lw) then 
                ! encourage vectorization for ECCKD with simpler loop indices
                jband = jg
              else
                jband = config%i_band_from_reordered_g_lw(jg)
              end if

              ! Add scaled cloud optical depth to clear-sky value
              od_region_cld(jg,jreg,jlev) = od_region_clear(jg,jlev) &
                  &  + od_cloud(jband,jlev,jcol)*od_scaling(jreg,jlev,jcol)

              if (config%do_lw_cloud_scattering) then
                ! Compute single-scattering albedo of gas-cloud
                ! combination
                scat_od_cloud = od_cloud(jband,jlev,jcol) &
                      &  * ssa_cloud(jband,jlev,jcol)*od_scaling(jreg,jlev,jcol)
                ssa_region_cld(jg,jreg,jlev) = (scat_od(jg)+scat_od_cloud) &
                      &  / od_region_cld(jg,jreg,jlev)

                ! Compute asymmetry factor of gas-cloud
                ! combination
                g_region_cld(jg,jreg,jlev) = &
                      &  (scat_od(jg)*g_clear(jg,jlev) + scat_od_cloud*g_cloud(jband,jlev,jcol)) &
                      &  / max((scat_od(jg) + scat_od_cloud), tiny(g_region_cld(jg,jreg,jlev)))
                ! equal to 0 if scat_od + scat_od_cloud  = 0
              end if
              ! else ssa_region is already set to zero

              ! Apply maximum cloud optical depth for stability in the 3D case
              od_region_cld(jg,jreg,jlev) = min(od_region_cld(jg,jreg,jlev), config%max_cloud_od_lw)
            end do
          end do
        end do

        ! Calculate two-stream variables gamma1 and gamma2 of
        ! all g-points and cloudy regions at once
        call calc_two_stream_gammas_lw(ng*(nreg-1)*nlev_cld, &
        &  ssa_region_cld, g_region_cld, &
        &  gamma1_cld(:,:,jtop:jbot), gamma2_cld(:,:,jtop:jbot))

        do jlev = jtop, jbot
          ! Get the indices corresponding to 3D g-points at this level
          j1 = inds(1,jlev)
          j2 = inds(2,jlev)

          ng3D = j2 - j1 +1

          if (ng3D < ng) then
            ! Some of the g points are to be treated using the
            ! conventional plane-parallel method

            ! Compute reflectance, transmittance and upward and downward
            ! sources, for each cloudy region, using the Meador-Weaver
            ! formulas for reflectance and transmittance, and equivalent
            ! solutions to the coupled ODEs for the sources.
            do jreg = 2, nreg
              call calc_ref_trans_lw(ng-ng3D, &
                  &  od_region_cld(ng3D+1:ng,jreg,jlev), &
                  &  ssa_region_cld(ng3D+1:ng,jreg,jlev), g_region_cld(ng3D+1:ng,jreg,jlev), &
                  ! array temps are created for the next two arguments
                  &  region_fracs(jreg,jlev,jcol)*planck_hl(ng3D+1:ng,jlev,jcol), &
                  &  region_fracs(jreg,jlev,jcol)*planck_hl(ng3D+1:ng,jlev+1,jcol), &
                  &  reflectance(ng3D+1:ng,jreg,jreg,jlev), &
                  &  transmittance(ng3D+1:ng,jreg,jreg,jlev), &
                  &  source_up(ng3D+1:ng,jreg,jlev), &
                  &  source_dn(ng3D+1:ng,jreg,jlev))
            end do
            n_calls_meador_weaver &
                &  = n_calls_meador_weaver + (ng-ng3D)*(nreg-1)
          end if

        end do

        if (ng3D_tot>0) then
          allocate(Gamma_z1(ng3D_tot, 2*nreg, 2*nreg), LU(ng3D_tot, 2*nreg, 2*nreg))
          allocate(planck_diff(ng3D_tot, 2*nreg), planck_top(ng3D_tot, 2*nreg))
          
          allocate(tmp_vectors(ng3D_tot, nreg))
          allocate(source_up_tmp(ng3D_tot,nreg), source_dn_tmp(ng3D_tot,nreg))

          allocate(gamma_z11(ng3D_tot, nreg, nreg), gamma_z12(ng3D_tot, nreg, nreg))
          allocate(gamma_z21(ng3D_tot, nreg, nreg), gamma_z22(ng3D_tot, nreg, nreg))

          allocate(transfer_rate(nreg, nreg, jtop:jbot))
          transfer_rate = 0.0_jprb
          Gamma_z1= 0.0_jprb
          planck_top = 0.0_jprb
          planck_diff = 0.0_jprb

          ! --------------------------------------------------------
          ! Compute transfer rates
          ! --------------------------------------------------------
          do jlev = jtop, jbot
            ! Depth of the current layer
            dz = layer_depth(jlev)

            do jreg = 1, nreg-1
              ! Compute lateral transfer rate from region jreg to
              ! jreg+1 following Hogan & Shonk (2013) Eq. 47, but
              ! multiplied by dz because the transfer rate is
              ! vertically integrated across the depth of the layer
              if (region_fracs(jreg,jlev,jcol) > epsilon(1.0_jprb)) then
                transfer_rate(jreg,jreg+1,jlev) = dz &
                      &  * edge_length(jreg,jlev) &
                      &  * tan_diffuse_angle_3d / region_fracs(jreg,jlev,jcol)
              end if
              ! Compute transfer rate from region jreg+1 to jreg
              if (region_fracs(jreg+1,jlev,jcol) > epsilon(1.0_jprb)) then
                transfer_rate(jreg+1,jreg,jlev) = dz &
                      &  * edge_length(jreg,jlev) &
                      &  * tan_diffuse_angle_3d / region_fracs(jreg+1,jlev,jcol)
              end if
            end do
            ! now set ()

            ! Compute transfer rates directly between regions 1 and 3
            if (edge_length(3,jlev) > 0.0_jprb) then
              if (region_fracs(1,jlev,jcol) > epsilon(1.0_jprb)) then
                transfer_rate(1,3,jlev) = dz &
                    &  * edge_length(3,jlev) &
                    &  * tan_diffuse_angle_3d / region_fracs(1,jlev,jcol)
              end if
              if (region_fracs(3,jlev,jcol) > epsilon(1.0_jprb)) then
                transfer_rate(3,1,jlev) = dz &
                    &  * edge_length(3,jlev) &
                    &  * tan_diffuse_angle_3d / region_fracs(3,jlev,jcol)
              end if
            end if
          end do

          ! Don't allow the transfer rate out of a region to be
          ! equivalent to a loss of exp(-10) through the layer
          transfer_rate = min(transfer_rate, config%max_3d_transfer_rate)

          do jlev = jtop, jbot 
            ! Get the indices corresponding to 3D g-points at this level
            j1 = inds(1,jlev)
            j2 = inds(2,jlev)
            ng3D = j2 - j1 +1

            ! Depth of the current layer
            dz = layer_depth(jlev)
            ! First clear-sky region
            jreg = 1
            jg = 1
            do jj = j1,j2 
              ! Write the diagonal elements of -Gamma1*z1
              Gamma_z1(jj,jreg,jreg)      = od_region_clear(jg,jlev)*gamma1_clear(jg,jlev)
              ! Write the diagonal elements of +Gamma2*z1
              if (config%do_lw_aerosol_scattering)  then 
                Gamma_z1(jj,jreg+nreg,jreg) = od_region_clear(jg,jlev)*gamma2_clear(jg,jlev)
                ! else zero because gamma2_clear is zero
              end if 
              ! Write the vectors corresponding to the inhomogeneous
              ! parts of the matrix ODE
              planck_top(jj,nreg+jreg) = od_region_clear(jg,jlev) &
                    &  *(1.0_jprb-ssa_clear(jg,jlev))*region_fracs(jreg,jlev,jcol) &
                    &  *planck_hl(jg,jlev,jcol)*LwDiffusivityWP
              planck_top(jj,jreg) = -planck_top(jj,nreg+jreg)
              planck_diff(jj,nreg+jreg) = od_region_clear(jg,jlev) &
                    &  * (1.0_jprb-ssa_clear(jg,jlev))*region_fracs(jreg,jlev,jcol) &
                    &  * (planck_hl(jg,jlev+1,jcol) &
                    &  -planck_hl(jg,jlev,jcol))*LwDiffusivityWP
              planck_diff(jj,jreg) = -planck_diff(jj,nreg+jreg)
              jg = jg + 1
            end do 
            ! Cloudy regions
            do jreg = 2,nreg
              jg = 1
              do jj = j1,j2 
                ! Write the diagonal elements of -Gamma1*z1
                Gamma_z1(jj,jreg,jreg)      = od_region_cld(jg,jreg,jlev)*gamma1_cld(jg,jreg,jlev)
                ! Write the diagonal elements of +Gamma2*z1
                Gamma_z1(jj,jreg+nreg,jreg) = od_region_cld(jg,jreg,jlev)*gamma2_cld(jg,jreg,jlev)

                ! Write the vectors corresponding to the inhomogeneous
                ! parts of the matrix ODE
                planck_top(jj,nreg+jreg) = od_region_cld(jg,jreg,jlev) &
                    &  *(1.0_jprb-ssa_region_cld(jg,jreg,jlev))*region_fracs(jreg,jlev,jcol) &
                    &  *planck_hl(jg,jlev,jcol)*LwDiffusivityWP
                planck_top(jj,jreg) = -planck_top(jj,nreg+jreg)
                planck_diff(jj,nreg+jreg) = od_region_cld(jg,jreg,jlev) &
                    &  * (1.0_jprb-ssa_region_cld(jg,jreg,jlev))*region_fracs(jreg,jlev,jcol) &
                    &  * (planck_hl(jg,jlev+1,jcol) &
                    &  -planck_hl(jg,jlev,jcol))*LwDiffusivityWP
                planck_diff(jj,jreg) = -planck_diff(jj,nreg+jreg)
                jg = jg + 1
              end do 
            end do

            ! Parameterization for the effective emissivity of the side
            ! of the cloud
            if (config%do_lw_side_emissivity &
              & .and. region_fracs(1,jlev,jcol) > 0.0_jprb .and. region_fracs(2,jlev,jcol) > 0.0_jprb &
              & .and. config%do_3d_effects &
              & .and. cloud%inv_cloud_effective_size(jcol,jlev) > 0.0_jprb) then
              aspect_ratio = 1.0_jprb / (min(cloud%inv_cloud_effective_size(jcol,jlev), &
                  &                         1.0_jprb / config%min_cloud_effective_size) &
                  &                     * region_fracs(1,jlev,jcol) * dz)
                do jg = 1, ng3D
                  lateral_od = (aspect_ratio / (nreg-1.0_jprb)) &
                      &  * sum(od_region_cld(jg,2:nreg,jlev)*(1.0_jprb-ssa_region_cld(jg,2:nreg,jlev)))
                  sqrt_1_minus_ssa = sqrt(1.0_jprb - ssa_region_cld(jg,2,jlev))
                  side_emiss_thick = 2.0_jprb * sqrt_1_minus_ssa / &
                      & (sqrt_1_minus_ssa &
                      &  + sqrt(1.0_jprb-ssa_region_cld(jg,2,jlev)*g_region_cld(jg,2,jlev)))
                  side_emiss(jg) = (side_emiss_thin - side_emiss_thick) &
                      &  / (lateral_od + 1.0_jprb)   + side_emiss_thick
                end do
            else
              side_emiss(1:ng3D) = 1.0_jprb
            end if


            do jreg = 1,nreg-1
              ! Add the terms assocated with 3D transport to Gamma1*z1.
              ! First the rate from region jreg to jreg+1
              Gamma_z1(j1:j2,jreg,jreg) = Gamma_z1(j1:j2,jreg,jreg) + transfer_rate(jreg,jreg+1,jlev)
              Gamma_z1(j1:j2,jreg+1,jreg) = -transfer_rate(jreg,jreg+1,jlev)

              ! Next the rate from region jreg+1 to jreg
              if (jreg > 1) then
                ! Flow between one cloudy region and another
                Gamma_z1(j1:j2,jreg+1,jreg+1) = Gamma_z1(j1:j2,jreg+1,jreg+1) &
                    &  + transfer_rate(jreg+1,jreg,jlev)
                Gamma_z1(j1:j2,jreg,jreg+1) = -transfer_rate(jreg+1,jreg,jlev)
              else
                ! Only for the lateral transfer between cloud and clear
                ! skies do we account for the effective emissivity of
                ! the side of the cloud
                jg = 1
                do jj = j1,j2 
                  Gamma_z1(jj,jreg+1,jreg+1) = Gamma_z1(jj,jreg+1,jreg+1) &
                      &  + side_emiss(jg) * transfer_rate(jreg+1,jreg,jlev)
                  Gamma_z1(jj,jreg,jreg+1) = -side_emiss(jg)*transfer_rate(jreg+1,jreg,jlev)
                  jg = jg + 1
                end do
              end if
            end do

            ! Possible flow between regions a and c
            if (edge_length(3,jlev) > 0.0_jprb) then
              Gamma_z1(j1:j2,1,1) = Gamma_z1(j1:j2,1,1) + transfer_rate(1,3,jlev)
              Gamma_z1(j1:j2,3,1) = -transfer_rate(1,3,jlev)
              Gamma_z1(j1:j2,3,3) = Gamma_z1(j1:j2,3,3) + side_emiss(1:ng3D) * transfer_rate(3,1,jlev)
              Gamma_z1(j1:j2,1,3) = -side_emiss(1:ng3D)*transfer_rate(3,1,jlev)
            end if

            ! Copy Gamma1*z1
            Gamma_z1(j1:j2,nreg+1:2*nreg,nreg+1:2*nreg) = -Gamma_z1(j1:j2,1:nreg,1:nreg)
            ! Copy Gamma2*z1
            Gamma_z1(j1:j2,1:nreg,nreg+1:2*nreg)  = -Gamma_z1(j1:j2,nreg+1:2*nreg,1:nreg)
          end do

          ! Compute the parts of the particular solution
          ! solution_diff(:,1:2*nreg) = solve_vec(ng3D_tot,ng3D_tot,2*nreg,Gamma_z1,planck_diff)
          ! solution_diff(:,1:2*nreg) = - solution_diff(:,1:2*nreg) 
          ! solution0(:,1:2*nreg) = solve_vec(ng3D_tot,ng3D,2*nreg,Gamma_z1, &
          !      &  solution_diff-planck_top)
          ! We can skip a redundant LU factorization by explicitly calling the subroutines
          ! call lu_factorization(ng3D_tot, ng3D_tot, 2*nreg, Gamma_z1, LU)
          call lu_factorization_lw(ng3D_tot, Gamma_z1, LU)

          associate(solution_diff=>planck_diff, solution0=>planck_top)

          call lu_substitution_lw_inplace(ng3D_tot, LU, solution_diff) ! overwrite planck_diff with solution
          solution_diff = - solution_diff

          planck_top = solution_diff - planck_top
          call lu_substitution_lw_inplace(ng3D_tot, LU, solution0) ! overwrite planck_top with solution

          ! Additional security on elements fed to matrix exponential
          ! in single precision
          if (jprb <= 4) then
            Gamma_z1 = min(Gamma_z1, 22.0_jprb)
          end if

          call expm_lw(ng3D_tot, Gamma_z1)

          ! Compute sources, reflectance and transmittance
          ! Avoid temporary arrays being made for these noncontiguous sections of Gamma, 
          ! which are used as input arguments more than once
          gamma_z11 = Gamma_z1(:,1:nreg,1:nreg) 
          gamma_z12 = Gamma_z1(:,1:nreg,nreg+1:2*nreg)
          gamma_z21 = Gamma_z1(:,nreg+1:2*nreg,1:nreg)
          gamma_z22 = Gamma_z1(:,nreg+1:2*nreg,nreg+1:2*nreg)

          ! Upward and downward sources due to emission within the layer
          tmp_vectors = mat_x_vec_3(ng3D_tot,  gamma_z12, solution0(:,nreg+1:2*nreg))
          do jreg = 1, nreg
            do jg = 1, ng3D_tot
              tmp_vectors(jg,jreg) =  solution0(jg,jreg) + solution_diff(jg,jreg) - tmp_vectors(jg,jreg)
            end do
          end do
          source_up_tmp  = solve_vec_3_ng(ng3D_tot,gamma_z11, tmp_vectors)

          do jreg = 1, nreg
            do jg = 1, ng3D_tot
              source_up_tmp(jg,jreg) =  solution0(jg,jreg) - source_up_tmp(jg,jreg) 
              tmp_vectors(jg,jreg)  = source_up_tmp(jg,jreg) - solution0(jg,jreg)
            end do
          end do
          source_dn_tmp = mat_x_vec_3(ng3D_tot,  gamma_z21, tmp_vectors) &
              &  + solution0(:,nreg+1:2*nreg) &
              &  - mat_x_vec_3(ng3D_tot, gamma_z22, solution0(:,nreg+1:2*nreg)) &
              &  + solution_diff(:,nreg+1:2*nreg)

          end associate ! solution0, solution_diff

          associate(reflectance_tmp=>gamma_z12, transmittance_tmp=>gamma_z11)

          ! Diffuse reflectance matrix
          ! solve_mat_inplace(n,A,B) : Solve AX=B, ovewriting B with X
          call solve_mat_3_inplace(ng3D_tot, gamma_z11, gamma_z12)
          reflectance_tmp = -reflectance_tmp

          ! Diffuse transmission matrix
          ! call  mat_x_mat_3(ng3D_tot, gamma_z21, reflectance_tmp, transmittance_tmp)
          ! transmittance_tmp = transmittance_tmp + gamma_z22
          do j2 = 1,3
            do j1 = 1,3
              transmittance_tmp(:,j1,j2) = gamma_z21(:,j1,1)*reflectance_tmp(:,1,j2) &
                  & + gamma_z21(:,j1,2)*reflectance_tmp(:,2,j2) &
                  & + gamma_z21(:,j1,3)*reflectance_tmp(:,3,j2) 
              transmittance_tmp(:,j1,j2) = transmittance_tmp(:,j1,j2) + gamma_z22(:,j1,j2)
            end do
          end do

          ! Write results to main arrays
          do jlev = jtop, jbot 
            j1 = inds(1,jlev)
            j2 = inds(2,jlev)
            ng3D = j2 - j1 +1

            do jreg = 1, nreg
              do jreg2 = 1, nreg
                transmittance(1:ng3D,jreg2,jreg,jlev) = transmittance_tmp(j1:j2,jreg2,jreg)
                reflectance(1:ng3D,jreg2,jreg,jlev) = reflectance_tmp(j1:j2,jreg2,jreg)
              end do
              source_up(1:ng3D,jreg,jlev) = source_up_tmp(j1:j2,jreg)
              source_dn(1:ng3D,jreg,jlev) = source_dn_tmp(j1:j2,jreg)
            end do
          end do

          end associate
          ! --------------------------------------------------------------
          ! ------------------- END CLOUDY COMPUTATIONS ------------------
          ! -------------------------------------------------------------- 
          ! Deallocations
          deallocate(transfer_rate, source_up_tmp, source_dn_tmp)
          deallocate(gamma_z11, gamma_z12, gamma_z21, gamma_z22)
          deallocate(tmp_vectors) !, solution_diff, solution0)
          deallocate(planck_diff, planck_top)
          deallocate(Gamma_z1, LU)

        end if

        deallocate(inds)
        deallocate(od_region_cld, ssa_region_cld, g_region_cld)
        deallocate(gamma1_cld, gamma2_cld)

        ! CLOUD WHILE CONSTRUCT: update clouds_below and find new jtop if true
        if (jbot==nlev) then 
          are_clouds_below=.false. ! surface reached
        else 
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
        end if
     end do

      ! --------------------------------------------------------
      ! Section 4: Compute total sources and albedos
      ! --------------------------------------------------------
      ! total_albedo(:,:,:,:) = 0.0_jprb
      ! total_source(:,:,:) = 0.0_jprb
      total_albedo(:,:,:,nlev+1) = 0.0_jprb

      if (config%do_clear) then
        total_albedo_clear(:,:) = 0.0_jprb
        total_source_clear(:,:) = 0.0_jprb
      end if

      ! Calculate the upwelling radiation emitted by the surface, and
      ! copy the surface albedo into total_albedo 
      do jreg = 1,nreg
        do jg = 1,ng
          ! region_fracs(jreg,nlev,jcol) is the fraction of each
          ! region in the lowest model level
          total_source(jg,jreg,nlev+1) = region_fracs(jreg,nlev,jcol)*emission(jg,jcol)
          total_albedo(jg,jreg,jreg,nlev+1) = albedo(jg,jcol)
        end do
      end do
      ! Equivalent surface values for computing clear-sky fluxes 
      if (config%do_clear) then
        do jg = 1,ng
          total_source_clear(jg,nlev+1) = emission(jg,jcol)
        end do
        ! In the case of surface albedo there is no dependence on
        ! cloud fraction so we can copy the all-sky value
        total_albedo_clear(1:ng,nlev+1) = total_albedo(1:ng,1,1,nlev+1)
      end if


      ! Loop back up through the atmosphere computing the total albedo
      ! and the total upwelling due to emission below each level
      do jlev = nlev,1,-1
        if (config%do_clear) then
          ! Use adding method for clear-sky arrays; note that there
          ! is no need to consider "above" and "below" quantities
          ! since with no cloud overlap to worry about, these are
          ! the same
          if (config%do_lw_aerosol_scattering)  then 
            inv_denom_scalar(:) = 1.0_jprb &
                &  / (1.0_jprb - total_albedo_clear(:,jlev+1)*ref_clear(:,jlev))
            total_albedo_clear(:,jlev) = ref_clear(:,jlev) &
                &  + trans_clear(:,jlev)*trans_clear(:,jlev)*total_albedo_clear(:,jlev+1) &
                &  * inv_denom_scalar(:)
            total_source_clear(:,jlev) = source_up_clear(:,jlev) &
                &  + trans_clear(:,jlev)*(total_source_clear(:,jlev+1) &
                &  + total_albedo_clear(:,jlev+1)*source_dn_clear(:,jlev)) &
                &  * inv_denom_scalar(:)
          else  ! ref_clear = 0 --> inv_denom_scalar = 1
            total_albedo_clear(:,jlev) = &
                &  trans_clear(:,jlev)*trans_clear(:,jlev)*total_albedo_clear(:,jlev+1)
            total_source_clear(:,jlev) = source_up_clear(:,jlev) &
                &  + trans_clear(:,jlev)*(total_source_clear(:,jlev+1) &
                &  + total_albedo_clear(:,jlev+1)*source_dn_clear(:,jlev))
          end if 
        end if

        if (is_clear_sky_layer(jlev)) then
          ! Clear-sky layer: use scalar adding method
          total_albedo_below = 0.0_jprb
          total_source_below = 0.0_jprb
          if (config%do_lw_aerosol_scattering)  then 
            inv_denom_scalar(:) = 1.0_jprb &
                &  / (1.0_jprb - total_albedo(:,1,1,jlev+1)*reflectance(:,1,1,jlev))
            total_albedo_below(:,1,1) = reflectance(:,1,1,jlev) &
                &  + transmittance(:,1,1,jlev)  * transmittance(:,1,1,jlev) &
                &  * total_albedo(:,1,1,jlev+1) * inv_denom_scalar(:)
            total_source_below(:,1) = source_up(:,1,jlev) &
                &  + transmittance(:,1,1,jlev)*(total_source(:,1,jlev+1) &
                &  + total_albedo(:,1,1,jlev+1)*source_dn(:,1,jlev)) &
                &  * inv_denom_scalar(:)
          else 
            total_albedo_below(:,1,1) = reflectance(:,1,1,jlev) &
                &  + transmittance(:,1,1,jlev)  * transmittance(:,1,1,jlev) &
                &  * total_albedo(:,1,1,jlev+1)
            total_source_below(:,1) = source_up(:,1,jlev) &
                &  + transmittance(:,1,1,jlev)*(total_source(:,1,jlev+1) &
                &  + total_albedo(:,1,1,jlev+1)*source_dn(:,1,jlev))
          end if 

        else if (config%do_3d_effects .or. &
             &   config%do_3d_lw_multilayer_effects) then
          ! Cloudy layer: use matrix adding method
          call identity_minus_mat_x_mat_3_lw(ng, total_albedo(:,:,:,jlev+1), &
                &  reflectance(:,:,:,jlev), denominator)
          ! Let's write the terms explicitly (avoids array temporaries)
          associate(B1=>total_albedo_below)
            ! B=C*D, then solve AX=B, finally Y = Z*X
            call mat_x_mat_3_lw(ng, &                             	! 1. B = C*D
                  &  total_albedo(:,:,:,jlev+1), transmittance(:,:,:,jlev), B1) 
            call solve_mat_3_lw(ng,denominator, B1, X)            	! 2. solve AX=B
            call mat_x_mat_3_lw(ng, transmittance(:,:,:,jlev), X, & ! 3. Y = Z*X ...
                  & total_albedo_below)
            total_albedo_below = total_albedo_below + reflectance(:,:,:,jlev) ! + R
          end associate
          ! Source
          b = mat_x_vec_3_lw(ng,total_albedo(:,:,:,jlev+1),source_dn(:,:,jlev))
          b = b + total_source(:,:,jlev+1)
          call solve_vec_3_lw_inplace(ng, denominator, b) ! Solve Ax=b, overwriting b with x
          total_source_below = source_up(:,:,jlev)  & 
              & + mat_x_vec_3_lw(ng,transmittance(:,:,:,jlev), b)
      
        else 
          ! Cloudy layer for which reflectance, transmittance and
          ! total_albedo matrices are diagonal
          total_albedo_below = 0.0_jprb
          total_source_below = 0.0_jprb
          do jreg = 1,nreg
            inv_denom_scalar(:) = 1.0_jprb / (1.0_jprb &
                 &  - total_albedo(:,jreg,jreg,jlev+1)*reflectance(:,jreg,jreg,jlev))
            total_albedo_below(:,jreg,jreg) = reflectance(:,jreg,jreg,jlev) &
                 &  + transmittance(:,jreg,jreg,jlev)*transmittance(:,jreg,jreg,jlev) &
                 &  * total_albedo(:,jreg,jreg,jlev+1) &
                 &  * inv_denom_scalar(:)
            total_source_below(:,jreg) = source_up(:,jreg,jlev) &
                 &  + transmittance(:,jreg,jreg,jlev)*(total_source(:,jreg,jlev+1) &
                 &  + total_albedo(:,jreg,jreg,jlev+1)*source_dn(:,jreg,jlev)) &
                 &  * inv_denom_scalar(:)
          end do

        end if

        ! Account for cloud overlap when converting albedo and
        ! source below a layer interface to the equivalent values
        ! just above
        if (is_clear_sky_layer(jlev) .and. is_clear_sky_layer(jlev-1)) then
          ! If both layers are cloud free, this is trivial...
          total_albedo(:,:,:,jlev) = 0.0_jprb
          total_albedo(:,1,1,jlev) = total_albedo_below(:,1,1)
          total_source(:,:,jlev) = 0.0_jprb
          total_source(:,1,jlev) = total_source_below(:,1)
        else 
          total_source(:,:,jlev) = singlemat_x_vec_lw(ng,&
               &  u_matrix(:,:,jlev), total_source_below)

!          if (config%do_3d_effects .or. config%do_3d_lw_multilayer_effects) then
          if (config%do_3d_lw_multilayer_effects) then
            ! Use the overlap matrices u_matrix and v_matrix
            total_albedo(:,:,:,jlev) = singlemat_x_mat(ng,ng,nreg,&
                 &  u_matrix(:,:,jlev), &
                 &  mat_x_singlemat(ng,ng,nreg,total_albedo_below,&
                 &  v_matrix(:,:,jlev)))
          else
            total_albedo(:,:,:,jlev) = 0.0_jprb
            ! "total_albedo" is diagonal and we wish to exclude
            ! anomalous horizontal transport described by Shonk &
            ! Hogan (2008).  Therefore, the operation we perform is
            ! essentially diag(total_albedo) = matmul(transpose(v_matrix),
            ! diag(total_albedo_below)).
            do jreg = 1,nreg
              do jreg2 = 1,nreg
                total_albedo(:,jreg,jreg,jlev) &
                     &  = total_albedo(:,jreg,jreg,jlev) &
                     &  + total_albedo_below(:,jreg2,jreg2) &
                     &  * v_matrix(jreg2,jreg,jlev)
              end do
            end do
          end if
        end if

      end do ! Reverse loop over levels

      ! --------------------------------------------------------
      ! Section 5: Compute fluxes
      ! --------------------------------------------------------

      ! Top-of-atmosphere fluxes into the regions of the top-most
      ! layer: zero since we assume no extraterrestrial longwave
      ! radiation (the shortwave scheme accounts for solar radiation
      ! in the "longwave" part of the spectrum)
      flux_dn_below = 0.0_jprb
      flux%lw_dn(jcol,1) = 0.0_jprb
      if (config%do_clear) then
        flux_dn_clear = 0.0_jprb
        flux%lw_dn_clear(jcol,1) = 0.0_jprb
      end if

      ! Store the outgoing longwave radiation at top-of-atmosphere
      flux%lw_up(jcol,1) = sum(sum(total_source(:,:,1),1))
      if (config%do_clear) then
        flux%lw_up_clear(jcol,1) = sum(total_source_clear(:,1))
      end if

      if (config%do_save_spectral_flux) then
        call indexed_sum(sum(total_source(:,:,1),2), &
             &           config%i_spec_from_reordered_g_lw, &
             &           flux%lw_up_band(:,jcol,1))
        flux%lw_dn_band(:,jcol,1) = 0.0_jprb
        if (config%do_clear) then
          call indexed_sum(total_source_clear(:,1), &
               &           config%i_spec_from_reordered_g_lw, &
               &           flux%lw_up_clear_band(:,jcol,1))
          flux%lw_dn_clear_band(:,jcol,1) = 0.0_jprb
        end if
      end if

      ! Final loop back down through the atmosphere to compute fluxes
      do jlev = 1,nlev
        if (config%do_clear) then
          ! Scalar operations for clear-sky fluxes
          if (config%do_lw_aerosol_scattering)  then 
            flux_dn_clear(:) = (trans_clear(:,jlev)*flux_dn_clear(:) &
                + ref_clear(:,jlev)*total_source_clear(:,jlev+1) &
                + source_dn_clear(:,jlev)) &
                / (1.0_jprb - ref_clear(:,jlev)*total_albedo_clear(:,jlev+1))
          else ! ref_clear = 0
            flux_dn_clear(:) = trans_clear(:,jlev)*flux_dn_clear(:) &
              + source_dn_clear(:,jlev)
          end if 
          flux_up_clear(:) = total_source_clear(:,jlev+1) &
               + total_albedo_clear(:,jlev+1)*flux_dn_clear(:)
        end if

        if (is_clear_sky_layer(jlev)) then
          ! Scalar operations for clear-sky layer
          flux_dn_above = 0.0_jprb
          flux_up_above = 0.0_jprb

          flux_dn_above(:,1) = (transmittance(:,1,1,jlev)*flux_dn_below(:,1) &
               &  + reflectance(:,1,1,jlev)*total_source(:,1,jlev+1) &
               &  + source_dn(:,1,jlev)) &
               &  / (1.0_jprb - reflectance(:,1,1,jlev)*total_albedo(:,1,1,jlev+1))
          flux_up_above(:,1) = total_source(:,1,jlev+1) &
               &  + total_albedo(:,1,1,jlev+1)*flux_dn_above(:,1)
        else if (config%do_3d_effects .or. config%do_3d_lw_multilayer_effects) then
          ! Matrix operations for cloudy layer
          call identity_minus_mat_x_mat_3_lw(ng,reflectance(:,:,:,jlev), &
                &  total_albedo(:,:,:,jlev+1), denominator)
          flux_dn_above = solve_vec_3_lw(ng,denominator, &
              &  mat_x_vec_3_lw(ng,transmittance(:,:,:,jlev),flux_dn_below) &
              &  + mat_x_vec_3_lw(ng,reflectance(:,:,:,jlev), &
              &  total_source(:,:,jlev+1)) &
              &  + source_dn(:,:,jlev))
          flux_up_above = mat_x_vec_3_lw(ng,total_albedo(:,:,:,jlev+1), &
              &  flux_dn_above) + total_source(:,:,jlev+1)
        else
          do jreg = 1,nreg
            ! Scalar operations for all regions, requiring that
            ! reflectance, transmittance and total_albedo are diagonal
            flux_dn_above(:,jreg) = (transmittance(:,jreg,jreg,jlev)*flux_dn_below(:,jreg) &
                 &  + reflectance(:,jreg,jreg,jlev)*total_source(:,jreg,jlev+1) &
                 &  + source_dn(:,jreg,jlev)) &
                 &  / (1.0_jprb - reflectance(:,jreg,jreg,jlev) &
                 &              * total_albedo(:,jreg,jreg,jlev+1))
            flux_up_above(:,jreg) = total_source(:,jreg,jlev+1) &
                 &  + total_albedo(:,jreg,jreg,jlev+1)*flux_dn_above(:,jreg)
          end do
        end if

        ! Account for overlap rules in translating fluxes just above
        ! a layer interface to the values just below
        if (is_clear_sky_layer(jlev) .and. is_clear_sky_layer(jlev+1)) then
          flux_dn_below = flux_dn_above
        else
          flux_dn_below = singlemat_x_vec_lw(ng,v_matrix(:,:,jlev+1), &
               &    flux_dn_above)
        end if

        ! Store the broadband fluxes
        sums_up = 0.0_jprb; sums_dn = 0.0_jprb
        if (is_clear_sky_layer(jlev)) then
#ifdef __NEC__
          !NEC$ shortloop
#else
          !$omp simd reduction(+:sums_up, sums_dn)
#endif
          do jg = 1, ng  
            sums_up = sums_up + flux_up_above(jg,1) 
            sums_dn = sums_dn + flux_dn_above(jg,1)
          end do
        else 
#ifdef __NEC__
          !NEC$ shortloop
#else
          !$omp simd reduction(+:sums_up, sums_dn)
#endif
          do jg = 1, ng  
            sums_up = sums_up + flux_up_above(jg,1) + flux_up_above(jg,2) + flux_up_above(jg,3)
            sums_dn = sums_dn + flux_dn_above(jg,1) + flux_dn_above(jg,2) + flux_dn_above(jg,3)
          end do
        end if
        flux%lw_up(jcol,jlev+1) = sums_up
        flux%lw_dn(jcol,jlev+1) = sums_dn

        if (config%do_clear) then
          sums_up = 0.0_jprb; sums_dn = 0.0_jprb
#ifdef __NEC__
          !NEC$ shortloop
#else
          !$omp simd reduction(+:sums_up, sums_dn)
#endif
          do jg = 1, ng  
            sums_up = sums_up + flux_up_clear(jg)
            sums_dn = sums_dn + flux_dn_clear(jg)
          end do
          flux%lw_up_clear(jcol,jlev+1) = sums_up
          flux%lw_dn_clear(jcol,jlev+1) = sums_dn
        end if

        if (config%do_save_spectral_flux) then
          call indexed_sum(sum(flux_up_above,2), &
               &           config%i_spec_from_reordered_g_lw, &
               &           flux%lw_up_band(:,jcol,jlev+1))
          call indexed_sum(sum(flux_dn_above,2), &
               &           config%i_spec_from_reordered_g_lw, &
               &           flux%lw_dn_band(:,jcol,jlev+1))
          if (config%do_clear) then
            call indexed_sum(flux_up_clear, &
                 &           config%i_spec_from_reordered_g_lw, &
                 &           flux%lw_up_clear_band(:,jcol,jlev+1))
            call indexed_sum(flux_dn_clear, &
                 &           config%i_spec_from_reordered_g_lw, &
                 &           flux%lw_dn_clear_band(:,jcol,jlev+1))
          end if
        end if
        
      end do ! Final loop over levels

      ! Store surface spectral downwelling fluxes, which at this point
      ! are at the surface
      flux%lw_dn_surf_g(:,jcol) = sum(flux_dn_above,2)
      if (config%do_clear) then
        flux%lw_dn_surf_clear_g(:,jcol) = flux_dn_clear
      end if

      ! Compute the longwave derivatives needed by Hogan and Bozzo
      ! (2015) approximate radiation update scheme
      if (config%do_lw_derivatives) then
        ! Note that at this point flux_up_above contains the spectral
        ! fluxes into the regions of the lowest layer; we sum over
        ! regions first to provide a simple spectral flux upwelling
        ! from the surface
        call calc_lw_derivatives_matrix(ng, nlev, nreg, jcol, transmittance, &
             &  u_matrix(:,:,:), sum(flux_up_above,2), flux%lw_derivatives)
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

    if (lhook) call dr_hook('radiation_spartacus_lw:solver_spartacus_lw',1,hook_handle)

  end subroutine solver_spartacus_lw

end module radiation_spartacus_lw