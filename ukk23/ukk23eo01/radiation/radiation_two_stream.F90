! radiation_two_stream.F90 - Compute two-stream coefficients
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
!   2017-05-04  P Dueben/R Hogan  Use JPRD where double precision essential
!   2017-07-12  R Hogan  Optimized LW coeffs in low optical depth case
!   2017-07-26  R Hogan  Added calc_frac_scattered_diffuse_sw routine
!   2017-10-23  R Hogan  Renamed single-character variables
!   2021-02-19  R Hogan  Security for shortwave singularity
!   2022-11-22  P Ukkonen/R Hogan  Single precision uses no double precision

#include "ecrad_config.h"

module radiation_two_stream

  use parkind1, only : jprb, jprd

  implicit none
  public

  ! Elsasser's factor: the effective factor by which the zenith
  ! optical depth needs to be multiplied to account for longwave
  ! transmission at all angles through the atmosphere.  Alternatively
  ! think of acos(1/lw_diffusivity) to be the effective zenith angle
  ! of longwave radiation.
  real(jprd), parameter :: LwDiffusivity   = 1.66_jprd
  real(jprb), parameter :: LwDiffusivityWP = 1.66_jprb ! Working precision version

  ! Shortwave diffusivity factor assumes hemispheric isotropy, assumed
  ! by Zdunkowski's scheme and most others; note that for efficiency
  ! this parameter is not used in the calculation of the gamma values,
  ! but is used in the SPARTACUS solver.
  real(jprb), parameter :: SwDiffusivity = 2.00_jprb

  ! Make minimum k value depend on precision, allowing to avoid JPRD altogether
#ifdef PARKIND1_SINGLE
  ! real(jprb), parameter :: KMin = 1.e4_jprb * epsilon(1._jprb)
  real(jprb), parameter :: KMinSw = 1.e-4_jprb
  real(jprb), parameter :: KMinLw = 1.e-4_jprb
#else
  real(jprb), parameter :: KMinSw = 1.e-12_jprb
  real(jprb), parameter :: KMinLw = 1.e-12_jprb
#endif
  ! The routines in this module can be called millions of times, so
  ! calling Dr Hook for each one may be a significant overhead.
  ! Uncomment the following to turn Dr Hook on.
!#define DO_DR_HOOK_TWO_STREAM

  ! Private parameters used in longwave calculations only
  real(jprb), parameter :: Half= 0.5_jprb
  real(jprb), parameter :: One = 1.0_jprb
  real(jprb), parameter :: Two = 2.0_jprb
  ! For optical depths lower than this value, use the linear
  ! approximation
  real(jprb), parameter :: OdThresholdLw = 1.0e-3_jprb

  private :: OdThresholdLw, Half, One, Two, KMinSw, KMinLw

! Allow size of inner dimension (number of g-points) to be known at compile time if NG_SW is defined
#ifdef NG_SW
  integer, parameter, private :: ng_sw = NG_SW
#else
#define ng_sw ng_sw_in
#endif

contains

  !---------------------------------------------------------------------
  ! Calculate the two-stream coefficients gamma1 and gamma2 for the
  ! longwave
  subroutine calc_two_stream_gammas_lw(ng, ssa, g, &
       &                               gamma1, gamma2)

#ifdef DO_DR_HOOK_TWO_STREAM
    use yomhook, only : lhook, dr_hook, jphook
#endif

    integer, intent(in) :: ng
    ! Sngle scattering albedo and asymmetry factor:
    real(jprb), intent(in),  dimension(ng) :: ssa, g
    real(jprb), intent(out), dimension(ng) :: gamma1, gamma2

    real(jprb) :: factor

    integer    :: jg

#ifdef DO_DR_HOOK_TWO_STREAM
    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_two_stream:calc_two_stream_gammas_lw',0,hook_handle)
#endif
! Added for DWD (2020)
!NEC$ shortloop
    do jg = 1, ng
      ! Fu et al. (1997), Eq 2.9 and 2.10:
      !      gamma1(jg) = LwDiffusivity * (1.0_jprb - 0.5_jprb*ssa(jg) &
      !           &                    * (1.0_jprb + g(jg)))
      !      gamma2(jg) = LwDiffusivity * 0.5_jprb * ssa(jg) &
      !           &                    * (1.0_jprb - g(jg))
      ! Reduce number of multiplications
      factor = (LwDiffusivity * 0.5_jprb) * ssa(jg)
      gamma1(jg) = LwDiffusivity - factor*(1.0_jprb + g(jg))
      gamma2(jg) = factor * (1.0_jprb - g(jg))
    end do

#ifdef DO_DR_HOOK_TWO_STREAM
    if (lhook) call dr_hook('radiation_two_stream:calc_two_stream_gammas_lw',1,hook_handle)
#endif

  end subroutine calc_two_stream_gammas_lw


  !---------------------------------------------------------------------
  ! Calculate the two-stream coefficients gamma1-gamma4 in the
  ! shortwave
  subroutine calc_two_stream_gammas_sw(ng, mu0, ssa, g, &
       &                               gamma1, gamma2, gamma3)

#ifdef DO_DR_HOOK_TWO_STREAM
    use yomhook, only : lhook, dr_hook, jphook
#endif

    integer, intent(in) :: ng
    ! Cosine of solar zenith angle, single scattering albedo and
    ! asymmetry factor:
    real(jprb), intent(in)                :: mu0
    real(jprb), intent(in),  dimension(ng) :: ssa, g
    real(jprb), intent(out), dimension(ng) :: gamma1, gamma2, gamma3

    real(jprb) :: factor

    integer    :: jg

#ifdef DO_DR_HOOK_TWO_STREAM
    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_two_stream:calc_two_stream_gammas_sw',0,hook_handle)
#endif

    ! Zdunkowski "PIFM" (Zdunkowski et al., 1980; Contributions to
    ! Atmospheric Physics 53, 147-66)
! Added for DWD (2020)
!NEC$ shortloop
    do jg = 1, ng
      !      gamma1(jg) = 2.0_jprb  - ssa(jg) * (1.25_jprb + 0.75_jprb*g(jg))
      !      gamma2(jg) = 0.75_jprb *(ssa(jg) * (1.0_jprb - g(jg)))
      !      gamma3(jg) = 0.5_jprb  - (0.75_jprb*mu0)*g(jg)
      ! Optimized version:
      factor = 0.75_jprb*g(jg)
      gamma1(jg) = 2.0_jprb  - ssa(jg) * (1.25_jprb + factor)
      gamma2(jg) = ssa(jg) * (0.75_jprb - factor)
      gamma3(jg) = 0.5_jprb  - mu0*factor
    end do

#ifdef DO_DR_HOOK_TWO_STREAM
    if (lhook) call dr_hook('radiation_two_stream:calc_two_stream_gammas_sw',1,hook_handle)
#endif

  end subroutine calc_two_stream_gammas_sw


  !---------------------------------------------------------------------
  ! Compute the longwave reflectance and transmittance to diffuse
  ! radiation using the Meador & Weaver formulas, as well as the
  ! upward flux at the top and the downward flux at the base of the
  ! layer due to emission from within the layer assuming a linear
  ! variation of Planck function within the layer.
  ! Optimized version: single precision, allow exp() and other expensive parts to 
  ! vectorize by placing conditional for small optical depths in a separate loop
  subroutine calc_ref_trans_lw(ng, &
       &    od, ssa, g, &
       &    planck_top, planck_bot, &
       &    reflectance, transmittance, source_up, source_dn, &
       &    gamma1_out, gamma2_out)

#ifdef DO_DR_HOOK_TWO_STREAM
    use yomhook, only : lhook, dr_hook, jphook
#endif
    integer, intent(in) :: ng

    ! Optical depth and single scattering albedo
    real(jprb), intent(in), dimension(ng) :: od

    ! Sngle scattering albedo and asymmetry factor:
    real(jprb), intent(in),  dimension(ng) :: ssa, g

        ! The Planck terms (functions of temperature) at the top and
    ! bottom of the layer
    real(jprb), intent(in), dimension(ng) :: planck_top, planck_bot

    ! The diffuse reflectance and transmittance, i.e. the fraction of
    ! diffuse radiation incident on a layer from either top or bottom
    ! that is reflected back or transmitted through
    real(jprb), intent(out), dimension(ng) :: reflectance, transmittance

    ! The upward emission at the top of the layer and the downward
    ! emission at its base, due to emission from within the layer
    real(jprb), intent(out), dimension(ng) :: source_up, source_dn

    real(jprb), intent(out), optional, target, dimension(ng) :: gamma1_out, gamma2_out

    integer  :: jg

    ! The two transfer coefficients from the two-stream
    ! differential equations
    real(jprb), dimension(ng), target             :: gamma1_loc, gamma2_loc
    real(jprb), dimension(:), contiguous, pointer :: gamma1, gamma2
    !real(jprb), dimension(ng)                     :: k_exponent, exponential
    real(jprb) :: reftrans_factor
    real(jprb) :: exponential2 ! = exp(-2*k_exponent*od)

    real(jprb) :: coeff, coeff_up_top, coeff_up_bot, coeff_dn_top, coeff_dn_bot, factor

#ifdef DO_DR_HOOK_TWO_STREAM
    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_two_stream:calc_ref_trans_lw',0,hook_handle)
#endif

    associate (exponential=>transmittance, k_exponent=>source_up)

    if (present(gamma1_out)) then
      gamma1 => gamma1_out
      gamma2 => gamma2_out
    else
      gamma1 => gamma1_loc
      gamma2 => gamma2_loc
    end if

    do jg = 1, ng
      factor = (LwDiffusivity * 0.5_jprb) * ssa(jg)
      gamma1(jg) = LwDiffusivity - factor*(1.0_jprb + g(jg))
      gamma2(jg) = factor * (1.0_jprb - g(jg))
      k_exponent(jg) = sqrt(max((gamma1(jg) - gamma2(jg)) * (gamma1(jg) + gamma2(jg)), &
                KMinLw)) ! Eq 18 of Meador & Weaver (1980)
    end do

    exponential = exp(-k_exponent*od)

    do jg = 1, ng
      exponential2 = exponential(jg)*exponential(jg)
      reftrans_factor = 1.0 / (k_exponent(jg)  + gamma1(jg) + (k_exponent(jg) - gamma1(jg))*exponential2)
      ! Meador & Weaver (1980) Eq. 25
      reflectance(jg) = gamma2(jg) * (1.0_jprb - exponential2) * reftrans_factor
      ! Meador & Weaver (1980) Eq. 26
      transmittance(jg) = 2.0_jprb * k_exponent(jg) * exponential(jg) * reftrans_factor
      !
      ! Toon et al. (JGR 1989) Eqs 26-27
      !
      !   if (od(jg) > 1.0e-8_jprb) then
      coeff = (planck_bot(jg)-planck_top(jg)) / (od(jg)*(gamma1(jg)+gamma2(jg)))
      coeff_up_top  =  coeff + planck_top(jg)
      coeff_up_bot  =  coeff + planck_bot(jg)
      coeff_dn_top  = -coeff + planck_top(jg)
      coeff_dn_bot  = -coeff + planck_bot(jg)
      source_up(jg) = (coeff_up_top - reflectance(jg) * coeff_dn_top - transmittance(jg) * coeff_up_bot)
      source_dn(jg) = (coeff_dn_bot - reflectance(jg) * coeff_up_bot - transmittance(jg) * coeff_dn_top)
    !   else
    !     source_up(jg) = 0._jprb
    !     source_dn(jg) = 0._jprb
    !   end if
    end do

    do jg = 1, ng
      if (od(jg) < OdThresholdLw) then
        if (od(jg) < 1.0e-8_jprb) then 
          source_up(jg) = 0._jprb
          source_dn(jg) = 0._jprb
        else 
          ! In between e-3..e-8, must use Robin's original code, otherwise not secure 
          ! for instance when LW aerosol scattering is turned on 
          source_up(jg) = (One - reflectance(jg) - transmittance(jg)) &
          &       * Half * (planck_top(jg) + planck_bot(jg))
          source_dn(jg) = source_up(jg)
        end if 
      end if 
    end do

    end associate
    
#ifdef DO_DR_HOOK_TWO_STREAM
    if (lhook) call dr_hook('radiation_two_stream:calc_ref_trans_lw',1,hook_handle)
#endif
  
  end subroutine calc_ref_trans_lw

  
  !---------------------------------------------------------------------
  ! Compute the longwave transmittance to diffuse radiation in the
  ! no-scattering case, as well as the upward flux at the top and the
  ! downward flux at the base of the layer due to emission from within
  ! the layer assuming a linear variation of Planck function within
  ! the layer.
  ! Optimized version
  subroutine calc_no_scattering_transmittance_lw(ng, &
    &    od, planck_top, planck_bot, transmittance, source_up, source_dn)

#ifdef DO_DR_HOOK_TWO_STREAM
 use yomhook, only : lhook, dr_hook, jphook
#endif
    integer, intent(in) :: ng

    ! Optical depth and single scattering albedo
    real(jprb), intent(in), dimension(ng) :: od

    ! The Planck terms (functions of temperature) at half levels and full levels
    real(jprb), intent(in), dimension(ng)  :: planck_top, planck_bot

    ! The diffuse transmittance, i.e. the fraction of diffuse
    ! radiation incident on a layer from either top or bottom that is
    ! reflected back or transmitted through
    real(jprb), intent(out), dimension(ng) :: transmittance

    ! The upward emission at the top of the layer and the downward
    ! emission at its base, due to emission from within the layer
    real(jprb), intent(out), dimension(ng) :: source_up, source_dn

    ! Local variables
    real(jprb) :: coeff, coeff_up_top, coeff_up_bot, coeff_dn_top, coeff_dn_bot !, planck_mean
    real(jprb), dimension(ng) :: source_up_low_od
    integer :: jg
    
#ifdef DO_DR_HOOK_TWO_STREAM
    real(jphook) :: hook_handle
    if (lhook) call dr_hook('radiation_two_stream:calc_no_scattering_transmittance_lw',0,hook_handle)
#endif

    associate(od_loc=>source_up)

      od_loc = LwDiffusivity*od
#ifndef DWD_TWO_STREAM_OPTIMIZATIONS
      transmittance = exp(-od_loc)
#endif
      do jg = 1, ng
#ifdef DWD_TWO_STREAM_OPTIMIZATIONS
          transmittance(jg) = exp(-od_loc(jg))
#endif
          source_up_low_od(jg) = od_loc(jg)*0.5_jprb*(planck_top(jg)+planck_bot(jg))
        ! Compute upward and downward emission assuming the Planck
        ! function to vary linearly with optical depth within the layer
        ! (e.g. Wiscombe , JQSRT 1976).
          od_loc(jg) = (planck_bot(jg)-planck_top(jg)) / od_loc(jg)
          coeff_up_top  =  od_loc(jg) + planck_top(jg)
          coeff_up_bot  =  od_loc(jg) + planck_bot(jg)
          coeff_dn_top  = -od_loc(jg) + planck_top(jg)
          coeff_dn_bot  = -od_loc(jg) + planck_bot(jg)
          source_up(jg) =  coeff_up_top - transmittance(jg) * coeff_up_bot
          source_dn(jg) =  coeff_dn_bot - transmittance(jg) * coeff_dn_top
      end do
    end associate
      
    do jg = 1, ng
      if (od(jg) < 1.0e-3) then
        ! Linear limit at low optical depth
        ! coeff = LwDiffusivity*od(jg)
        ! transmittance(jg) = 1.0_jprb - coeff
        ! source_up(jg) = coeff * 0.5_jprb * (planck_top(jg)+planck_bot(jg))
        source_up(jg) = source_up_low_od(jg)
        source_dn(jg) = source_up(jg)
      end if
    end do

#ifdef DO_DR_HOOK_TWO_STREAM
 if (lhook) call dr_hook('radiation_two_stream:calc_no_scattering_transmittance_lw',1,hook_handle)
#endif

  end subroutine calc_no_scattering_transmittance_lw

  ! Upward and downward emission at layer edges assuming a linear variation of 
  ! Planck function within the layer, using Clough et al., 1992, doi:10.1029/92JD01419, eq 13
  ! with Pade approximant to avoid conditional for low tau (eq 15) - considerably faster 
  ! but introduces errors which may be considerable
  subroutine calc_no_scattering_transmittance_lw_pade(ng, &
    &    od, planck_top, planck_bot, transmittance, source_up, source_dn)

#ifdef DO_DR_HOOK_TWO_STREAM
 use yomhook, only : lhook, dr_hook, jphook
#endif

    integer, intent(in) :: ng

    ! Optical depth and single scattering albedo
    real(jprb), intent(in), dimension(ng) :: od

    ! The Planck terms (functions of temperature) at half levels and full levels
    real(jprb), intent(in), dimension(ng)  :: planck_top, planck_bot

    ! The diffuse transmittance, i.e. the fraction of diffuse
    ! radiation incident on a layer from either top or bottom that is
    ! reflected back or transmitted through
    real(jprb), intent(out), dimension(ng) :: transmittance

    ! The upward emission at the top of the layer and the downward
    ! emission at its base, due to emission from within the layer
    real(jprb), intent(out), dimension(ng) :: source_up, source_dn

    real(jprb) :: planck_fl, coeff
    integer :: jg

#ifdef DO_DR_HOOK_TWO_STREAM
    real(jphook) :: hook_handle
    if (lhook) call dr_hook('radiation_two_stream:calc_no_scattering_transmittance_lw_pade',0,hook_handle)
#endif

    associate(od_loc=>source_up)

      od_loc = LwDiffusivity*od
      transmittance = exp(-od_loc)

      do jg = 1, ng
        ! Alternative source computation using a Pade approximant for the linear-in-tau solution
        ! This method requires no conditional but introduces some error
        ! See Clough et al., 1992, doi:10.1029/92JD01419, Eq 15
        coeff = 0.2_jprb * od_loc(jg)
        planck_fl = 0.5_jprb * (planck_top(jg) + planck_bot(jg))
        source_dn(jg)  = (1.0_jprb-transmittance(jg)) * (planck_fl + coeff*planck_bot(jg)) / (1 + coeff)          
        source_up(jg)  = (1.0_jprb-transmittance(jg)) * (planck_fl + coeff*planck_top(jg)) / (1 + coeff)
      end do 
  
    end associate

#ifdef DO_DR_HOOK_TWO_STREAM
 if (lhook) call dr_hook('radiation_two_stream:calc_no_scattering_transmittance_lw_pade',1,hook_handle)
#endif

  end subroutine calc_no_scattering_transmittance_lw_pade

  !---------------------------------------------------------------------
  ! Compute the shortwave reflectance and transmittance to diffuse
  ! radiation using the Meador & Weaver formulas, as well as the
  ! "direct" reflection and transmission, which really means the rate
  ! of transfer of direct solar radiation (into a plane perpendicular
  ! to the direct beam) into diffuse upward and downward streams at
  ! the top and bottom of the layer, respectively.  Finally,
  ! trans_dir_dir is the transmittance of the atmosphere to direct
  ! radiation with no scattering.
  ! Faster version with all variables in working precision (jprb), inlined
  ! computations of gamma, and exponentials kept out of loops to facilitate vectorization
  ! Adapted from RTE (Radiative Transfer for Energetics) code (Robert Pincus)
  ! Optimizations by PU and RH
  subroutine calc_ref_trans_sw(ng_sw_in, nlev, mu0, od, ssa, &
       &      asymmetry, ref_diff, trans_diff, &
       &      ref_dir, trans_dir_diff, trans_dir_dir, &
       &      gamma1_out, gamma2_out, gamma3_out)

#ifdef DO_DR_HOOK_TWO_STREAM
    use yomhook, only : lhook, dr_hook, jphook
#endif
    integer, intent(in) :: ng_sw_in, nlev

    ! Cosine of solar zenith angle
    real(jprb), intent(in) :: mu0

    ! Optical depth, single scattering albedo and asymmetry factor
    real(jprb), intent(in), dimension(ng_sw*nlev) :: od, ssa, asymmetry

    ! The direct reflectance and transmittance, i.e. the fraction of
    ! incoming direct solar radiation incident at the top of a layer
    ! that is either reflected back (ref_dir) or scattered but
    ! transmitted through the layer to the base (trans_dir_diff)
    real(jprb), intent(out), dimension(ng_sw*nlev) :: ref_dir, trans_dir_diff

    ! The diffuse reflectance and transmittance, i.e. the fraction of
    ! diffuse radiation incident on a layer from either top or bottom
    ! that is reflected back or transmitted through
    real(jprb), intent(out), dimension(ng_sw*nlev) :: ref_diff, trans_diff

    ! Transmittance of the direct been with no scattering
    real(jprb), intent(out), dimension(ng_sw*nlev) :: trans_dir_dir

    ! The transfer coefficients from the two-stream differential
    ! equations
#ifndef DWD_TWO_STREAM_OPTIMIZATIONS
    real(jprb), intent(out), optional, target, dimension(ng_sw*nlev) :: gamma1_out, gamma2_out, gamma3_out
    real(jprb), dimension(ng_sw*nlev), target :: gamma1_loc, gamma2_loc, gamma3_loc
    real(jprb), dimension(:), contiguous, pointer :: gamma1, gamma2, gamma3
    real(jprb), dimension(ng_sw*nlev) :: gamma4
    real(jprb), dimension(ng_sw*nlev) :: k_exponent !alpha1, alpha2, k_exponent
    real(jprb), dimension(ng_sw*nlev) :: exponential  ! = exp(-k_exponent*od)
#else
    real(jprb) :: gamma1, gamma2, gamma3, gamma4 
    real(jprb) :: alpha1, alpha2, k_exponent
    real(jprb) :: exponential ! = exp(-k_exponent*od)
#endif
    
    real(jprb) :: reftrans_factor, factor
    real(jprb) :: exponential2 ! = exp(-2*k_exponent*od)
    real(jprb) :: k_mu0, k_gamma3, k_gamma4
    real(jprb) :: k_2_exponential, one_minus_kmu0_sqr
    integer    :: jg

#ifdef DO_DR_HOOK_TWO_STREAM
    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_two_stream:calc_ref_trans_sw',0,hook_handle)
#endif

#ifndef DWD_TWO_STREAM_OPTIMIZATIONS
    if (present(gamma1_out)) then
      gamma1 => gamma1_out
      gamma2 => gamma2_out
      gamma3 => gamma3_out
    else
      gamma1 => gamma1_loc
      gamma2 => gamma2_loc
      gamma3 => gamma3_loc
    end if

    ! Transmission of unscattered direct solar radiation through the
    ! layer
    ! trans_dir_dir(:) = exp(-od(:) * (1.0_jprb/mu0))

    trans_dir_dir = -od(:) * (1.0_jprb/mu0)
    trans_dir_dir = max(trans_dir_dir, -1000.0_jprb)
    ! print *, "min max tr", minval(trans_dir_dir), maxval(trans_dir_dir)
    ! GCC 9.3 strange error: intermediate values of ~ -8000 cause a FPE when vectorizing exp(), but not 
    ! in non-vectorized loop, nor with larger negative values!

    trans_dir_dir = exp(trans_dir_dir)

    associate(alpha1=>trans_dir_diff, alpha2=>ref_dir)

    do jg = 1, ng_sw*nlev

      factor = 0.75_jprb*asymmetry(jg)

      gamma1(jg) = 2.0_jprb  - ssa(jg) * (1.25_jprb + factor)
      gamma2(jg) = ssa(jg) * (0.75_jprb - factor)
      gamma3(jg) = 0.5_jprb  - mu0*factor
      gamma4(jg) = 1.0_jprb - gamma3(jg)

      alpha1(jg) = gamma1(jg)*gamma4(jg) + gamma2(jg)*gamma3(jg) ! Eq. 16
      alpha2(jg) = gamma1(jg)*gamma3(jg) + gamma2(jg)*gamma4(jg) ! Eq. 17

      ! The following computation is very slow on ifort+Ryzen when inside a loop!
#if !defined(__INTEL_COMPILER)
      k_exponent(jg) = sqrt(max((gamma1(jg) - gamma2(jg)) * (gamma1(jg) + gamma2(jg)), &
            &       KMinSw)) ! Eq 18
    end do  
#else 
    end do    
    k_exponent = sqrt(max((gamma1- gamma2) * (gamma1 + gamma2), KMinSw)) ! Eq 18
#endif

    exponential(:) = exp(-k_exponent(:)*od(:))

    do jg = 1, ng_sw*nlev
      k_mu0 = k_exponent(jg)*mu0
      one_minus_kmu0_sqr = 1.0_jprb - k_mu0*k_mu0
      k_gamma3 = k_exponent(jg)*gamma3(jg)
      k_gamma4 = k_exponent(jg)*gamma4(jg)
      exponential2 = exponential(jg)*exponential(jg)
      k_2_exponential = 2.0_jprb * k_exponent(jg) * exponential(jg)
      reftrans_factor = 1.0_jprb / (k_exponent(jg) + gamma1(jg) + (k_exponent(jg) - gamma1(jg))*exponential2)
        
      ! Meador & Weaver (1980) Eq. 25
      ref_diff(jg) = gamma2(jg) * (1.0_jprb - exponential2) * reftrans_factor
      !ref_diff(jg)       = max(0.0_jprb, min(ref_diff(jg)), 1.0_jprb)

      ! Meador & Weaver (1980) Eq. 26, with security (which is
      ! sometimes needed, but apparently not on ref_diff)
      trans_diff(jg) = max(0.0_jprb, min(k_2_exponential * reftrans_factor, 1.0_jprb-ref_diff(jg)))
        
      ! Here we need mu0 even though it wasn't in Meador and Weaver
      ! because we are assuming the incoming direct flux is defined to
      ! be the flux into a plane perpendicular to the direction of the
      ! sun, not into a horizontal plane
      reftrans_factor = mu0 * ssa(jg) * reftrans_factor &
            &  / merge(one_minus_kmu0_sqr, epsilon(1.0_jprb), abs(one_minus_kmu0_sqr) > epsilon(1.0_jprb))
      
      ! Meador & Weaver (1980) Eq. 14, multiplying top & bottom by
      ! exp(-k_exponent*od) in case of very high optical depths
      ref_dir(jg) = reftrans_factor &
           &  * ( (1.0_jprb - k_mu0) * (alpha2(jg) + k_gamma3) &
           &     -(1.0_jprb + k_mu0) * (alpha2(jg) - k_gamma3)*exponential2 &
           &     -k_2_exponential*(gamma3(jg) - alpha2(jg)*mu0)*trans_dir_dir(jg) )
        
      ! Meador & Weaver (1980) Eq. 15, multiplying top & bottom by
      ! exp(-k_exponent*od), minus the 1*exp(-od/mu0) term
      ! representing direct unscattered transmittance.
      trans_dir_diff(jg) = reftrans_factor * &
           & ( k_2_exponential*(gamma4(jg) + alpha1(jg)*mu0) &
           & - trans_dir_dir(jg) &
           & * ( (1.0_jprb + k_mu0) * (alpha1(jg) + k_gamma4) &
           &    -(1.0_jprb - k_mu0) * (alpha1(jg) - k_gamma4) * exponential2) )

      ! Final check that ref_dir + trans_dir_diff <= 1
      ref_dir(jg)        = max(0.0_jprb, min(ref_dir(jg), mu0*(1.0_jprb-trans_dir_dir(jg))))
      trans_dir_diff(jg) = max(0.0_jprb, min(trans_dir_diff(jg), mu0*(1.0_jprb-trans_dir_dir(jg))-ref_dir(jg)))

    end do

    end associate

#else
    ! GPU-capable and vector-optimized version for ICON
    do jg = 1, ng_sw*nlev

      trans_dir_dir(jg) = max(-max(od(jg) * (1.0_jprb/mu0),0.0_jprb),-1000.0_jprb)
      trans_dir_dir(jg) = exp(trans_dir_dir(jg))

      ! Zdunkowski "PIFM" (Zdunkowski et al., 1980; Contributions to
      ! Atmospheric Physics 53, 147-66)
      factor = 0.75_jprb*asymmetry(jg)

      gamma1 = 2.0_jprb  - ssa(jg) * (1.25_jprb + factor)
      gamma2 = ssa(jg) * (0.75_jprb - factor)
      gamma3 = 0.5_jprb  - mu0*factor
      gamma4 = 1.0_jprb - gamma3

      alpha1 = gamma1*gamma4 + gamma2*gamma3 ! Eq. 16
      alpha2 = gamma1*gamma3 + gamma2*gamma4 ! Eq. 17
#ifdef PARKIND1_SINGLE
      k_exponent = sqrt(max((gamma1 - gamma2) * (gamma1 + gamma2), 1.0e-6_jprb))  ! Eq 18
#else
      k_exponent = sqrt(max((gamma1 - gamma2) * (gamma1 + gamma2), 1.0e-12_jprb)) ! Eq 18
#endif

      exponential = exp(-k_exponent*od(jg))

      k_mu0 = k_exponent*mu0
      one_minus_kmu0_sqr = 1.0_jprb - k_mu0*k_mu0
      k_gamma3 = k_exponent*gamma3
      k_gamma4 = k_exponent*gamma4
      exponential2 = exponential*exponential
      k_2_exponential = 2.0_jprb * k_exponent * exponential
      reftrans_factor = 1.0_jprb / (k_exponent + gamma1 + (k_exponent - gamma1)*exponential2)
        
      ! Meador & Weaver (1980) Eq. 25
      ref_diff(jg) = gamma2 * (1.0_jprb - exponential2) * reftrans_factor
        
      ! Meador & Weaver (1980) Eq. 26
      trans_diff(jg) = k_2_exponential * reftrans_factor
        
      ! Here we need mu0 even though it wasn't in Meador and Weaver
      ! because we are assuming the incoming direct flux is defined to
      ! be the flux into a plane perpendicular to the direction of the
      ! sun, not into a horizontal plane
      reftrans_factor = mu0 * ssa(jg) * reftrans_factor &
            &  / merge(one_minus_kmu0_sqr, epsilon(1.0_jprb), abs(one_minus_kmu0_sqr) > epsilon(1.0_jprb))
      
      ! Meador & Weaver (1980) Eq. 14, multiplying top & bottom by
      ! exp(-k_exponent*od) in case of very high optical depths
      ref_dir(jg) = reftrans_factor &
           &  * ( (1.0_jprb - k_mu0) * (alpha2 + k_gamma3) &
           &     -(1.0_jprb + k_mu0) * (alpha2 - k_gamma3)*exponential2 &
           &     -k_2_exponential*(gamma3 - alpha2*mu0)*trans_dir_dir(jg) )
        
      ! Meador & Weaver (1980) Eq. 15, multiplying top & bottom by
      ! exp(-k_exponent*od), minus the 1*exp(-od/mu0) term
      ! representing direct unscattered transmittance.
      trans_dir_diff(jg) = reftrans_factor * ( k_2_exponential*(gamma4 + alpha1*mu0) &
           & - trans_dir_dir(jg) &
           & * ( (1.0_jprb + k_mu0) * (alpha1 + k_gamma4) &
           &    -(1.0_jprb - k_mu0) * (alpha1 - k_gamma4) * exponential2) )

      ! Final check that ref_dir + trans_dir_diff <= 1
      ref_dir(jg)        = max(0.0_jprb, min(ref_dir(jg), mu0*(1.0_jprb-trans_dir_dir(jg))))
      trans_dir_diff(jg) = max(0.0_jprb, min(trans_dir_diff(jg), mu0*(1.0_jprb-trans_dir_dir(jg))-ref_dir(jg)))

    end do
#endif

#ifdef DO_DR_HOOK_TWO_STREAM
    if (lhook) call dr_hook('radiation_two_stream:calc_ref_trans_sw',1,hook_handle)
#endif
 
  end subroutine calc_ref_trans_sw
  
  !---------------------------------------------------------------------
  ! Compute the fraction of shortwave transmitted diffuse radiation
  ! that is scattered during its transmission, used to compute
  ! entrapment in SPARTACUS
  subroutine calc_frac_scattered_diffuse_sw(ng, od, &
       &      gamma1, gamma2, frac_scat_diffuse)
    
#ifdef DO_DR_HOOK_TWO_STREAM
    use yomhook, only : lhook, dr_hook, jphook
#endif

    integer, intent(in) :: ng

    ! Optical depth
    real(jprb), intent(in), dimension(ng) :: od

    ! The first two transfer coefficients from the two-stream
    ! differentiatial equations (computed by calc_two_stream_gammas)
    real(jprb), intent(in), dimension(ng) :: gamma1, gamma2

    ! The fraction of shortwave transmitted diffuse radiation that is
    ! scattered during its transmission
    real(jprb), intent(out), dimension(ng) :: frac_scat_diffuse

    real(jprd) :: k_exponent, reftrans_factor
    real(jprd) :: exponential  ! = exp(-k_exponent*od)
    real(jprd) :: exponential2 ! = exp(-2*k_exponent*od)
    real(jprd) :: k_2_exponential
    integer    :: jg

#ifdef DO_DR_HOOK_TWO_STREAM
    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_two_stream:calc_frac_scattered_diffuse_sw',0,hook_handle)
#endif

! Added for DWD (2020)
!NEC$ shortloop
    do jg = 1, ng
      ! Note that if the minimum value is reduced (e.g. to 1.0e-24)
      ! then noise starts to appear as a function of solar zenith
      ! angle
      k_exponent = sqrt(max((gamma1(jg) - gamma2(jg)) * (gamma1(jg) + gamma2(jg)), &
           &       1.0e-12_jprd)) ! Eq 18
      exponential = exp(-k_exponent*od(jg))
      exponential2 = exponential*exponential
      k_2_exponential = 2.0_jprd * k_exponent * exponential
        
      reftrans_factor = 1.0_jprd / (k_exponent + gamma1(jg) + (k_exponent - gamma1(jg))*exponential2)
        
      ! Meador & Weaver (1980) Eq. 26.
      ! Until 1.1.8, used LwDiffusivity instead of 2.0, although the
      ! effect is very small
      !      frac_scat_diffuse(jg) = 1.0_jprb - min(1.0_jprb,exp(-LwDiffusivity*od(jg)) &
      !           &  / max(1.0e-8_jprb, k_2_exponential * reftrans_factor))
      frac_scat_diffuse(jg) = 1.0_jprb &
           &  - min(1.0_jprb,exp(-2.0_jprb*od(jg)) &
           &  / max(1.0e-8_jprb, k_2_exponential * reftrans_factor))
    end do
    
#ifdef DO_DR_HOOK_TWO_STREAM
    if (lhook) call dr_hook('radiation_two_stream:calc_frac_scattered_diffuse_sw',1,hook_handle)
#endif
 
  end subroutine calc_frac_scattered_diffuse_sw

end module radiation_two_stream
