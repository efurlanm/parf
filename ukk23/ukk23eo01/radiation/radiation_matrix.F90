! radiation_matrix.F90 - SPARTACUS matrix operations
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
!   2018-10-15  R. Hogan    Added fast_expm_exchange_[23]
!   2020-12-15  P. Ukkonen  Optimized kernels
!
! This module provides the neccessary mathematical functions for the
! SPARTACUS radiation scheme: matrix multiplication, matrix solvers
! and matrix exponentiation, but (a) multiple matrices are operated on
! at once with array access indended to facilitate vectorization, and
! (b) optimization for 2x2 and 3x3 matrices.  There is probably
! considerable scope for further optimization. Note that this module
! is not used by the McICA solver.

module radiation_matrix

  use parkind1, only : jprb, jprd
#ifdef USE_TIMING
    ! Timing library
    use gptl,                  only: gptlstart, gptlstop
#endif

  implicit none
  public

  ! Codes to describe sparseness pattern, where the SHORTWAVE
  ! pattern is of the form:
  ! (x x x)
  ! (x x x)
  ! (0 0 x)
  ! where each element may itself be a square matrix.  
  integer, parameter :: IMatrixPatternDense     = 0
  integer, parameter :: IMatrixPatternShortwave = 1
#ifdef USE_TIMING
  integer, private :: ret
#endif

! Allow size of inner dimension (number of g-points) to be known at compile time for some routines
#ifdef NG_SW
    integer, parameter, private :: ng_sw = NG_SW
    integer, parameter, private :: ng_sww = 6*NG_SW
#else
#define ng_sw ng_sw_in
#define ng_sww ng_sw_in
#endif
#ifdef NG_LW
    integer, parameter, private :: ng_lw = NG_LW
#else
#define ng_lw ng_lw_in
#endif

  private :: solve_vec_2, solve_vec_3, solve_mat_2, solve_mat_n, &
       &     diag_mat_right_divide_3

  interface fast_expm_exchange
    module procedure fast_expm_exchange_2, fast_expm_exchange_3
  end interface fast_expm_exchange

contains

  ! --- MATRIX-VECTOR MULTIPLICATION ---

  !---------------------------------------------------------------------
  ! Treat A as n m-by-m square matrices (with the n dimension varying
  ! fastest) and b as n m-element vectors, and perform matrix-vector
  ! multiplications on first iend pairs
  function mat_x_vec(n,iend,m,A,b,do_top_left_only_in)

    use yomhook, only : lhook, dr_hook, jphook

    integer,    intent(in)                   :: n, m, iend
    real(jprb), intent(in), dimension(:,:,:) :: A
    real(jprb), intent(in), dimension(:,:)   :: b
    logical,    intent(in), optional         :: do_top_left_only_in
    real(jprb),             dimension(iend,m):: mat_x_vec

    integer :: j1, j2
    logical :: do_top_left_only

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_matrix:mat_x_vec',0,hook_handle)

    if (present(do_top_left_only_in)) then
      do_top_left_only = do_top_left_only_in
    else
      do_top_left_only = .false.
    end if

    ! Array-wise assignment
    mat_x_vec = 0.0_jprb

    if (do_top_left_only) then
      mat_x_vec(1:iend,1) = A(1:iend,1,1)*b(1:iend,1)
    else
      do j1 = 1,m
        do j2 = 1,m
          mat_x_vec(1:iend,j1) = mat_x_vec(1:iend,j1) &
               &               + A(1:iend,j1,j2)*b(1:iend,j2)
        end do
      end do
    end if

    if (lhook) call dr_hook('radiation_matrix:mat_x_vec',1,hook_handle)

  end function mat_x_vec

  pure function mat_x_vec_3(ng,A,b)

    integer,    intent(in)                    :: ng
    real(jprb), intent(in), dimension(ng,3,3) :: A
    real(jprb), intent(in), dimension(ng,3)   :: b
    real(jprb),             dimension(ng,3):: mat_x_vec_3

    integer :: jg

    do jg = 1, ng   
      ! both inner and outer loop of the matrix loops j1 and j2 unrolled
      ! inner loop:            j2=1               j2=2                  j2=3 
      mat_x_vec_3(jg,1) = A(jg,1,1)*b(jg,1) + A(jg,1,2)*b(jg,2) + A(jg,1,3)*b(jg,3) ! j1=1
      mat_x_vec_3(jg,2) = A(jg,2,1)*b(jg,1) + A(jg,2,2)*b(jg,2) + A(jg,2,3)*b(jg,3) ! j1=2
      mat_x_vec_3(jg,3) = A(jg,3,1)*b(jg,1) + A(jg,3,2)*b(jg,2) + A(jg,3,3)*b(jg,3) ! j1=3
    end do

  end function mat_x_vec_3

  pure function mat_x_vec_3_lw(ng_lw_in,A,b)

    integer,    intent(in)                       :: ng_lw_in
    real(jprb), intent(in), dimension(ng_lw,3,3) :: A
    real(jprb), intent(in), dimension(ng_lw,3)   :: b
    real(jprb),             dimension(ng_lw,3):: mat_x_vec_3_lw

    integer :: jg

    do jg = 1, ng_lw   
      ! both inner and outer loop of the matrix loops j1 and j2 unrolled
      ! inner loop:            j2=1               j2=2                  j2=3 
      mat_x_vec_3_lw(jg,1) = A(jg,1,1)*b(jg,1) + A(jg,1,2)*b(jg,2) + A(jg,1,3)*b(jg,3) ! j1=1
      mat_x_vec_3_lw(jg,2) = A(jg,2,1)*b(jg,1) + A(jg,2,2)*b(jg,2) + A(jg,2,3)*b(jg,3) ! j1=2
      mat_x_vec_3_lw(jg,3) = A(jg,3,1)*b(jg,1) + A(jg,3,2)*b(jg,2) + A(jg,3,3)*b(jg,3) ! j1=3
    end do

  end function mat_x_vec_3_lw

  pure function mat_x_vec_3_sw(ng_sw_in,A,b,do_top_left_only_in)

    integer,    intent(in)                        :: ng_sw_in
    real(jprb), intent(in), dimension(ng_sw,3,3)  :: A
    real(jprb), intent(in), dimension(ng_sw,3)    :: b
    real(jprb),             dimension(ng_sw,3)    :: mat_x_vec_3_sw
    logical,    intent(in), optional              :: do_top_left_only_in
    logical :: do_top_left_only
    integer :: jg

    if (present(do_top_left_only_in)) then
      do_top_left_only = do_top_left_only_in
    else
      do_top_left_only = .false.
    end if

    if (do_top_left_only) then
      mat_x_vec_3_sw      = 0.0_jprb
      mat_x_vec_3_sw(:,1) = A(:,1,1)*b(:,1)
    else
      do jg = 1, ng_sw   
        ! inner loop:            j2=1               j2=2                  j2=3 
        mat_x_vec_3_sw(jg,1) = A(jg,1,1)*b(jg,1) + A(jg,1,2)*b(jg,2) + A(jg,1,3)*b(jg,3) ! j1=1
        mat_x_vec_3_sw(jg,2) = A(jg,2,1)*b(jg,1) + A(jg,2,2)*b(jg,2) + A(jg,2,3)*b(jg,3) ! j1=2
        mat_x_vec_3_sw(jg,3) = A(jg,3,1)*b(jg,1) + A(jg,3,2)*b(jg,2) + A(jg,3,3)*b(jg,3) ! j1=3
      end do
    end if
  end function mat_x_vec_3_sw

  !---------------------------------------------------------------------
  ! Treat A as an m-by-m square matrix and b as n m-element vectors
  ! (with the n dimension varying fastest), and perform matrix-vector
  ! multiplications on first iend pairs
  function singlemat_x_vec(n,iend,m,A,b)

    use yomhook, only : lhook, dr_hook, jphook

    integer,    intent(in)                    :: n, m, iend
    real(jprb), intent(in), dimension(m,m)    :: A
    real(jprb), intent(in), dimension(:,:)    :: b
    real(jprb),             dimension(iend,m) :: singlemat_x_vec

    integer    :: j1, j2
    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_matrix:single_mat_x_vec',0,hook_handle)

    ! Array-wise assignment
    singlemat_x_vec = 0.0_jprb

    do j1 = 1,m
      do j2 = 1,m
        singlemat_x_vec(1:iend,j1) = singlemat_x_vec(1:iend,j1) &
             &                    + A(j1,j2)*b(1:iend,j2)
      end do
    end do

    if (lhook) call dr_hook('radiation_matrix:single_mat_x_vec',1,hook_handle)

  end function singlemat_x_vec

  pure function singlemat_x_vec_sw(ng_sw_in,A,b)
    integer,    intent(in)                      :: ng_sw_in
    real(jprb), intent(in), dimension(3,3)      :: A
    real(jprb), intent(in), dimension(ng_sw,3)  :: b
    real(jprb),             dimension(ng_sw,3)  :: singlemat_x_vec_sw
    integer    :: jg

    do jg = 1, ng_sw 
      ! both inner and outer loop of the matrix loops j1 and j2 unrolled
      ! inner loop:                   j2=1             j2=2             j2=3 
      singlemat_x_vec_sw(jg,1) = A(1,1)*b(jg,1) + A(1,2)*b(jg,2) + A(1,3)*b(jg,3) ! j1=1
      singlemat_x_vec_sw(jg,2) = A(2,1)*b(jg,1) + A(2,2)*b(jg,2) + A(2,3)*b(jg,3) ! j1=2
      singlemat_x_vec_sw(jg,3) = A(3,1)*b(jg,1) + A(3,2)*b(jg,2) + A(3,3)*b(jg,3) ! j1=3
    end do
  
  end function singlemat_x_vec_sw

  pure function singlemat_x_vec_lw(ng_lw_in,A,b)
    integer,    intent(in)                      :: ng_lw_in
    real(jprb), intent(in), dimension(3,3)      :: A
    real(jprb), intent(in), dimension(ng_lw,3)  :: b
    real(jprb),             dimension(ng_lw,3)  :: singlemat_x_vec_lw
    integer    :: jg

    do jg = 1, ng_lw
      ! both inner and outer loop of the matrix loops j1 and j2 unrolled
      ! inner loop:                   j2=1             j2=2             j2=3 
      singlemat_x_vec_lw(jg,1) = A(1,1)*b(jg,1) + A(1,2)*b(jg,2) + A(1,3)*b(jg,3) ! j1=1
      singlemat_x_vec_lw(jg,2) = A(2,1)*b(jg,1) + A(2,2)*b(jg,2) + A(2,3)*b(jg,3) ! j1=2
      singlemat_x_vec_lw(jg,3) = A(3,1)*b(jg,1) + A(3,2)*b(jg,2) + A(3,3)*b(jg,3) ! j1=3
    end do
  
  end function singlemat_x_vec_lw
  
  ! --- SQUARE MATRIX-MATRIX MULTIPLICATION ---

  !---------------------------------------------------------------------
  ! Treat A and B each as n m-by-m square matrices (with the n
  ! dimension varying fastest) and perform matrix multiplications on
  ! all n matrix pairs
  function mat_x_mat(n,iend,m,A,B,i_matrix_pattern)

    use yomhook, only : lhook, dr_hook, jphook

    integer,    intent(in)                      :: n, m, iend
    integer,    intent(in), optional            :: i_matrix_pattern
    real(jprb), intent(in), dimension(:,:,:)    :: A, B

    real(jprb),             dimension(iend,m,m) :: mat_x_mat
    integer    :: j1, j2, j3
    integer    :: mblock, m2block
    integer    :: i_actual_matrix_pattern
    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_matrix:mat_x_mat',0,hook_handle)

    if (present(i_matrix_pattern)) then
      i_actual_matrix_pattern = i_matrix_pattern
    else
      i_actual_matrix_pattern = IMatrixPatternDense
    end if

    ! Array-wise assignment
    mat_x_mat = 0.0_jprb

    if (i_actual_matrix_pattern == IMatrixPatternShortwave) then
      ! Matrix has a sparsity pattern
      !     (C D E)
      ! A = (F G H)
      !     (0 0 I)
      mblock = m/3
      m2block = 2*mblock 
      ! Do the top-left (C, D, F, G)
      do j2 = 1,m2block
        do j1 = 1,m2block
          do j3 = 1,m2block
            mat_x_mat(1:iend,j1,j2) = mat_x_mat(1:iend,j1,j2) &
                 &                  + A(1:iend,j1,j3)*B(1:iend,j3,j2)
          end do
        end do
      end do
      do j2 = m2block+1,m
        ! Do the top-right (E & H)
        do j1 = 1,m2block
          do j3 = 1,m
            mat_x_mat(1:iend,j1,j2) = mat_x_mat(1:iend,j1,j2) &
                 &                  + A(1:iend,j1,j3)*B(1:iend,j3,j2)
          end do
        end do
        ! Do the bottom-right (I)
        do j1 = m2block+1,m
          do j3 = m2block+1,m
            mat_x_mat(1:iend,j1,j2) = mat_x_mat(1:iend,j1,j2) &
                 &                  + A(1:iend,j1,j3)*B(1:iend,j3,j2)
          end do
        end do
      end do
    else
      ! Ordinary dense matrix
      do j2 = 1,m
        do j1 = 1,m
          do j3 = 1,m
            mat_x_mat(1:iend,j1,j2) = mat_x_mat(1:iend,j1,j2) &
                 &                  + A(1:iend,j1,j3)*B(1:iend,j3,j2)
          end do
        end do
      end do
    end if

    if (lhook) call dr_hook('radiation_matrix:mat_x_mat',1,hook_handle)

  end function mat_x_mat

  pure subroutine mat_x_mat_3(ng,A,B,C)

    integer,    intent(in)                          :: ng
    real(jprb), intent(in),     dimension(ng,3,3)  :: A, B
    !dir$ assume_aligned A:64, B:64
    real(jprb), intent(inout),  dimension(ng,3,3)  :: C
    !dir$ assume_aligned C:64

    integer    :: j1, j2

    do j2 = 1,3
      do j1 = 1,3
        C(:,j1,j2) = A(:,j1,1)*B(:,1,j2) + A(:,j1,2)*B(:,2,j2) + A(:,j1,3)*B(:,3,j2) 
      end do
    end do

  end subroutine mat_x_mat_3

  pure subroutine mat_x_mat_3_sw(ng_sw_in,A,B,C)

    integer,    intent(in)                            :: ng_sw_in
    real(jprb), intent(in),     dimension(ng_sw,3,3)  :: A, B
    !dir$ assume_aligned A:64, B:64
    real(jprb), intent(inout),  dimension(ng_sw,3,3)  :: C
    !dir$ assume_aligned C:64
    integer    :: j1, j2!, i,j,k

    do j2 = 1,3
      do j1 = 1,3
        C(:,j1,j2) = A(:,j1,1)*B(:,1,j2) + A(:,j1,2)*B(:,2,j2) + A(:,j1,3)*B(:,3,j2) 
      end do
    end do

  end subroutine mat_x_mat_3_sw

  pure subroutine mat_x_mat_3_lw(ng_lw_in,A,B,C)

    integer,    intent(in)                            :: ng_lw_in
    real(jprb), intent(in),     dimension(ng_lw,3,3)  :: A, B
    real(jprb), intent(inout),  dimension(ng_lw,3,3)  :: C
    integer    :: j1, j2

    do j2 = 1,3
      do j1 = 1,3
        C(:,j1,j2) = A(:,j1,1)*B(:,1,j2) + A(:,j1,2)*B(:,2,j2) + A(:,j1,3)*B(:,3,j2) 
      end do
    end do

  end subroutine mat_x_mat_3_lw

  pure function mat_x_mat_6(ng,A,B)

    integer,    intent(in)                    :: ng
    real(jprb), intent(in), dimension(ng,6,6) :: A, B
    real(jprb),             dimension(ng,6,6) :: mat_x_mat_6
    integer    :: j1, j2

    do j2 = 1,6
      do j1 = 1,6
        mat_x_mat_6(:,j1,j2) =  A(:,j1,1)*B(:,1,j2) + A(:,j1,2)*B(:,2,j2) + A(:,j1,3)*B(:,3,j2) &
            &                 + A(:,j1,4)*B(:,4,j2) + A(:,j1,5)*B(:,5,j2) + A(:,j1,6)*B(:,6,j2) 
      end do
    end do

  end function mat_x_mat_6

  !---------------------------------------------------------------------
  ! Treat A as an m-by-m matrix and B as n m-by-m square matrices
  ! (with the n dimension varying fastest) and perform matrix
  ! multiplications on the first iend matrix pairs
  function singlemat_x_mat(n,iend,m,A,B)

    use yomhook, only : lhook, dr_hook, jphook

    integer,    intent(in)                      :: n, m, iend
    real(jprb), intent(in), dimension(m,m)      :: A
    real(jprb), intent(in), dimension(:,:,:)    :: B
    real(jprb),             dimension(iend,m,m) :: singlemat_x_mat

    integer    :: j1, j2, j3
    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_matrix:singlemat_x_mat',0,hook_handle)

    ! Array-wise assignment
    singlemat_x_mat = 0.0_jprb

    do j2 = 1,m
      do j1 = 1,m
        do j3 = 1,m
          singlemat_x_mat(1:iend,j1,j2) = singlemat_x_mat(1:iend,j1,j2) &
               &                        + A(j1,j3)*B(1:iend,j3,j2)
        end do
      end do
    end do

    if (lhook) call dr_hook('radiation_matrix:singlemat_x_mat',1,hook_handle)

  end function singlemat_x_mat

  pure function singlemat_x_mat_3_sw(ng_sw_in,A,B)

    integer,    intent(in)                        :: ng_sw_in
    real(jprb), intent(in), dimension(3,3)        :: A
    real(jprb), intent(in), dimension(ng_sw,3,3)  :: B
    real(jprb),             dimension(ng_sw,3,3)  :: singlemat_x_mat_3_sw

    integer    :: j2, j3, jg

    ! unroll loops, just one SIMD write per matrix element
    do j2 = 1,3
      do jg = 1, ng_sw  ! inner loop:     j3=1                j3=2                j3=3
        singlemat_x_mat_3_sw(jg,1,j2) = A(1,1)*B(jg,1,j2) + A(1,2)*B(jg,2,j2) + A(1,3)*B(jg,3,j2) ! j1=1
        singlemat_x_mat_3_sw(jg,2,j2) = A(2,1)*B(jg,1,j2) + A(2,2)*B(jg,2,j2) + A(2,3)*B(jg,3,j2) ! j1=2
        singlemat_x_mat_3_sw(jg,3,j2) = A(3,1)*B(jg,1,j2) + A(3,2)*B(jg,2,j2) + A(3,3)*B(jg,3,j2) ! j1=3
      end do
    end do

  end function singlemat_x_mat_3_sw

  !---------------------------------------------------------------------
  ! Treat B as an m-by-m matrix and A as n m-by-m square matrices
  ! (with the n dimension varying fastest) and perform matrix
  ! multiplications on the first iend matrix pairs
  function mat_x_singlemat(n,iend,m,A,B)

    use yomhook, only : lhook, dr_hook, jphook

    integer,    intent(in)                      :: n, m, iend
    real(jprb), intent(in), dimension(:,:,:)    :: A
    real(jprb), intent(in), dimension(m,m)      :: B

    real(jprb),             dimension(iend,m,m) :: mat_x_singlemat
    integer    :: j1, j2, j3
    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_matrix:mat_x_singlemat',0,hook_handle)

    ! Array-wise assignment
    mat_x_singlemat = 0.0_jprb

    do j2 = 1,m
      do j1 = 1,m
        do j3 = 1,m
          mat_x_singlemat(1:iend,j1,j2) = mat_x_singlemat(1:iend,j1,j2) &
               &                        + A(1:iend,j1,j3)*B(j3,j2)
        end do
      end do
    end do

    if (lhook) call dr_hook('radiation_matrix:mat_x_singlemat',1,hook_handle)

  end function mat_x_singlemat

  pure function mat_x_singlemat_3_sw(ng_sw_in, A,B)

    integer,    intent(in)                        :: ng_sw_in
    real(jprb), intent(in), dimension(ng_sw,3,3)  :: A
    real(jprb), intent(in), dimension(3,3)        :: B

    real(jprb),             dimension(ng_sw,3,3)  :: mat_x_singlemat_3_sw
    integer    :: j1, j2, j3, jg

    do j2 = 1,3
      do jg = 1,ng_sw  ! inner loop:      j3=1                j3=2                j3=3
        mat_x_singlemat_3_sw(jg,1,j2) = A(jg,1,1)*B(1,j2) + A(jg,1,2)*B(2,j2) +  A(jg,1,3)*B(3,j2) ! j1 = 1
        mat_x_singlemat_3_sw(jg,2,j2) = A(jg,2,1)*B(1,j2) + A(jg,2,2)*B(2,j2) +  A(jg,2,3)*B(3,j2) ! j1 = 2
        mat_x_singlemat_3_sw(jg,3,j2) = A(jg,3,1)*B(1,j2) + A(jg,3,2)*B(2,j2) +  A(jg,3,3)*B(3,j2) ! j1 = 3
      end do
    end do

  end function mat_x_singlemat_3_sw

  !---------------------------------------------------------------------
  ! Compute I-A*B where I is the identity matrix and A & B are n
  ! m-by-m square matrices
  function identity_minus_mat_x_mat(n,iend,m,A,B,i_matrix_pattern)

    use yomhook, only : lhook, dr_hook, jphook

    integer,    intent(in)                   :: n, m, iend
    integer,    intent(in), optional         :: i_matrix_pattern
    real(jprb), intent(in), dimension(:,:,:) :: A, B
    real(jprb),             dimension(iend,m,m) :: identity_minus_mat_x_mat

    integer    :: j
    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_matrix:identity_mat_x_mat',0,hook_handle)

    if (present(i_matrix_pattern)) then
      identity_minus_mat_x_mat = mat_x_mat(n,iend,m,A,B,i_matrix_pattern)
    else
      identity_minus_mat_x_mat = mat_x_mat(n,iend,m,A,B)
    end if

    identity_minus_mat_x_mat = - identity_minus_mat_x_mat
    do j = 1,m
      identity_minus_mat_x_mat(1:iend,j,j) &
           &     = 1.0_jprb + identity_minus_mat_x_mat(1:iend,j,j)
    end do

    if (lhook) call dr_hook('radiation_matrix:identity_mat_x_mat',1,hook_handle)

  end function identity_minus_mat_x_mat

  pure subroutine identity_minus_mat_x_mat_3_sw(ng_sw_in,A,B,C)

    integer,    intent(in)                            :: ng_sw_in
    real(jprb), intent(in),     dimension(ng_sw,3,3)  :: A, B
    real(jprb), intent(inout),  dimension(ng_sw,3,3)  :: C

    integer    :: j1,j2

    do j2 = 1,3
      do j1 = 1,3
        C(:,j1,j2) = A(:,j1,1)*B(:,1,j2) + A(:,j1,2)*B(:,2,j2) + A(:,j1,3)*B(:,3,j2) 
      end do
    end do

    C = - C
    do j1 = 1,3
      C(:,j1,j1) = 1.0_jprb + C(:,j1,j1)
    end do

  end subroutine identity_minus_mat_x_mat_3_sw

  subroutine identity_minus_mat_x_mat_3_lw(ng_lw_in,A,B,C)

    integer,    intent(in)                            :: ng_lw_in
    real(jprb), intent(in),     dimension(ng_lw,3,3)  :: A, B
    real(jprb), intent(inout),  dimension(ng_lw,3,3)  :: C
    integer    :: j1,j2

    do j2 = 1,3
      do j1 = 1,3
        C(:,j1,j2) = A(:,j1,1)*B(:,1,j2) + A(:,j1,2)*B(:,2,j2) + A(:,j1,3)*B(:,3,j2) 
      end do
    end do

    C = - C
    do j1 = 1,3
      C(:,j1,j1) = 1.0_jprb + C(:,j1,j1)
    end do

  end subroutine identity_minus_mat_x_mat_3_lw

  !---------------------------------------------------------------------
  ! Replacement for matmul in the case that the first matrix is sparse
  function sparse_x_dense(sparse, dense)

    real(jprb), intent(in) :: sparse(:,:), dense(:,:)
    real(jprb) :: sparse_x_dense(size(sparse,1),size(dense,2))

    integer :: j1, j2, j3 ! Loop indices
    integer :: n1, n2, n3 ! Array sizes

    n1 = size(sparse,1)
    n2 = size(sparse,2)
    n3 = size(dense,2)
    
    sparse_x_dense = 0.0_jprb
    do j2 = 1,n2
      do j1 = 1,n1
        if (sparse(j1,j2) /= 0.0_jprb) then
          sparse_x_dense(j1,:) = sparse_x_dense(j1,:) + sparse(j1,j2)*dense(j2,:)
        end if
      end do
    end do
    
  end function sparse_x_dense

  ! --- REPEATEDLY SQUARE A MATRIX ---

  !---------------------------------------------------------------------
  ! Square m-by-m matrix "A" nrepeat times. A will be corrupted by
  ! this function.
  function repeated_square(m,A,nrepeat,i_matrix_pattern)
    integer,    intent(in)           :: m, nrepeat
    real(jprb), intent(inout)        :: A(m,m)
    integer,    intent(in), optional :: i_matrix_pattern
    real(jprb)                       :: repeated_square(m,m)

    integer :: j1, j2, j3, j4
    integer :: mblock, m2block
    integer :: i_actual_matrix_pattern

    if (present(i_matrix_pattern)) then
      i_actual_matrix_pattern = i_matrix_pattern
    else
      i_actual_matrix_pattern = IMatrixPatternDense
    end if

    if (i_actual_matrix_pattern == IMatrixPatternShortwave) then
      ! Matrix has a sparsity pattern
      !     (C D E)
      ! A = (F G H)
      !     (0 0 I)
      mblock = m/3
      m2block = 2*mblock
      do j4 = 1,nrepeat
        repeated_square = 0.0_jprb
        ! Do the top-left (C, D, F & G)
        do j2 = 1,m2block
          do j1 = 1,m2block
            do j3 = 1,m2block
              repeated_square(j1,j2) = repeated_square(j1,j2) &
                   &                 + A(j1,j3)*A(j3,j2)
            end do
          end do
        end do
        do j2 = m2block+1, m
          ! Do the top-right (E & H)
          do j1 = 1,m2block
            do j3 = 1,m
              repeated_square(j1,j2) = repeated_square(j1,j2) &
                   &                 + A(j1,j3)*A(j3,j2)
            end do
          end do
          ! Do the bottom-right (I)
          do j1 = m2block+1, m
            do j3 = m2block+1,m
              repeated_square(j1,j2) = repeated_square(j1,j2) &
                   &                 + A(j1,j3)*A(j3,j2)
            end do
          end do
        end do
        if (j4 < nrepeat) then
          A = repeated_square
        end if
      end do
    else
      ! Ordinary dense matrix
      do j4 = 1,nrepeat
        repeated_square = 0.0_jprb
        do j2 = 1,m
          do j1 = 1,m
            do j3 = 1,m
              repeated_square(j1,j2) = repeated_square(j1,j2) &
                   &                 + A(j1,j3)*A(j3,j2)
            end do
          end do
        end do
        if (j4 < nrepeat) then
          A = repeated_square
        end if
      end do
    end if

  end function repeated_square

  pure subroutine repeated_square_6(nrepeat,A,B)
    integer,    intent(in)           :: nrepeat
    real(jprb), intent(inout)        :: A(6,6)
    real(jprb), intent(out)          :: B(6,6)

    integer :: j1, j2, j3, j4

    ! Ordinary dense matrix
    do j4 = 1,nrepeat
      do j2 = 1,6
        do j1 = 1,6
          B(j1,j2) =  A(j1,1)*A(1,j2) + A(j1,2)*A(2,j2) + &
               &      A(j1,3)*A(3,j2) + A(j1,4)*A(4,j2) + &
               &      A(j1,5)*A(5,j2) + A(j1,6)*A(6,j2) 
        end do
      end do
      if (j4 < nrepeat) then
        A = B
      end if
    end do

  end subroutine repeated_square_6

  pure subroutine repeated_square_sw_9(nrepeat,A,B)
    integer,    intent(in)            :: nrepeat
    real(jprb), intent(inout)         :: A(9,9)
    real(jprb), intent(out)           :: B(9,9)

    integer :: j1, j2, j3

    do j3 = 1,nrepeat
        B = 0.0_jprb
        ! Do the top-left (C, D, F & G)
        do j2 = 1,6  
          do j1 = 1,6 
            B(j1,j2) = A(j1,1)*A(1,j2) + A(j1,2)*A(2,j2) &
            &          + A(j1,3)*A(3,j2) + A(j1,4)*A(4,j2) &
            &          + A(j1,5)*A(5,j2) + A(j1,6)*A(6,j2)
          end do
        end do
        do j2 = 7,9
          ! Do the top-right (E & H)
          do j1 =  1,6
            B(j1,j2) = A(j1,1)*A(1,j2) + A(j1,2)*A(2,j2) &
              &       + A(j1,3)*A(3,j2) + A(j1,4)*A(4,j2) + A(j1,5)*A(5,j2) &
              &       + A(j1,6)*A(6,j2) + A(j1,7)*A(7,j2) + A(j1,8)*A(8,j2) &
              &       + A(j1,9)*A(9,j2)
          end do
          ! Do the bottom-right (I)
          do j1 = 7,9
            B(j1,j2) =  A(j1,7)*A(7,j2) + A(j1,8)*A(8,j2) + A(j1,9)*A(9,j2)  
          end do
        end do
  
        if (j3 < nrepeat) then
            A = B
        end if
    end do
  end subroutine repeated_square_sw_9


  ! --- SOLVE LINEAR EQUATIONS ---

  !---------------------------------------------------------------------
  ! Solve Ax=b to obtain x.  Version optimized for 2x2 matrices using
  ! Cramer's method: "A" contains n 2x2 matrices and "b" contains n
  ! 2-element vectors; returns A^-1 b.
  pure subroutine solve_vec_2(n,iend,A,b,x)

    integer,    intent(in)  :: n, iend
    real(jprb), intent(in)  :: A(:,:,:)
    real(jprb), intent(in)  :: b(:,:)
    real(jprb), intent(out) :: x(:,:)

    real(jprb) :: inv_det(iend)

    inv_det = 1.0_jprb / (  A(1:iend,1,1)*A(1:iend,2,2) &
         &                - A(1:iend,1,2)*A(1:iend,2,1))

    x(1:iend,1) = inv_det*(A(1:iend,2,2)*b(1:iend,1)-A(1:iend,1,2)*b(1:iend,2))
    x(1:iend,2) = inv_det*(A(1:iend,1,1)*b(1:iend,2)-A(1:iend,2,1)*b(1:iend,1))

  end subroutine solve_vec_2


  !---------------------------------------------------------------------
  ! Solve AX=B to obtain X, i.e. the matrix right-hand-side version of
  ! solve_vec_2, with A, X and B all containing n 2x2 matrices;
  ! returns A^-1 B using Cramer's method.
  pure subroutine solve_mat_2(n,iend,A,B,X)
    integer,    intent(in)  :: n, iend
    real(jprb), intent(in)  :: A(:,:,:)
    real(jprb), intent(in)  :: B(:,:,:)
    real(jprb), intent(out) :: X(:,:,:)

    real(jprb) :: inv_det(iend)

    inv_det = 1.0_jprb / (  A(1:iend,1,1)*A(1:iend,2,2) &
         &                - A(1:iend,1,2)*A(1:iend,2,1))

    X(1:iend,1,1) = inv_det*( A(1:iend,2,2)*B(1:iend,1,1) &
         &                   -A(1:iend,1,2)*B(1:iend,2,1))
    X(1:iend,2,1) = inv_det*( A(1:iend,1,1)*B(1:iend,2,1) &
         &                   -A(1:iend,2,1)*B(1:iend,1,1))
    X(1:iend,1,2) = inv_det*( A(1:iend,2,2)*B(1:iend,1,2) &
         &                   -A(1:iend,1,2)*B(1:iend,2,2))
    X(1:iend,2,2) = inv_det*( A(1:iend,1,1)*B(1:iend,2,2) &
         &                   -A(1:iend,2,1)*B(1:iend,1,2))

  end subroutine solve_mat_2


  !---------------------------------------------------------------------
  ! Solve Ax=b optimized for 3x3 matrices, using LU
  ! factorization and substitution without pivoting.
  pure subroutine solve_vec_3(n,iend,A,b,x)
    integer,    intent(in)  :: n, iend
    real(jprb), intent(in)  :: A(:,:,:)
    real(jprb), intent(in)  :: b(:,:)
    real(jprb), intent(out) :: x(:,:)

    real(jprb), dimension(iend) :: L21, L31, L32
    real(jprb), dimension(iend) :: U22, U23, U33
    real(jprb), dimension(iend) :: y2, y3

    ! Some compilers unfortunately don't support assocate
    !    associate (U11 => A(:,1,1), U12 => A(:,1,2), U13 => A(1,3), &
    !         y1 => b(:,1), x1 => solve_vec3(:,1), &
    !         x2 => solve_vec3(:,2), x3 => solve_vec3(:,3))

    ! LU decomposition:
    !     ( 1        )   (U11 U12 U13)
    ! A = (L21  1    ) * (    U22 U23)
    !     (L31 L32  1)   (        U33)
    L21 = A(1:iend,2,1) / A(1:iend,1,1)
    L31 = A(1:iend,3,1) / A(1:iend,1,1)
    U22 = A(1:iend,2,2) - L21*A(1:iend,1,2)
    U23 = A(1:iend,2,3) - L21*A(1:iend,1,3)
    L32 =(A(1:iend,3,2) - L31*A(1:iend,1,2)) / U22
    U33 = A(1:iend,3,3) - L31*A(1:iend,1,3) - L32*U23

    ! Solve Ly = b by forward substitution
    y2 = b(1:iend,2) - L21*b(1:iend,1)
    y3 = b(1:iend,3) - L31*b(1:iend,1) - L32*y2

    ! Solve Ux = y by back substitution
    x(1:iend,3) = y3/U33
    x(1:iend,2) = (y2 - U23*x(1:iend,3)) / U22
    x(1:iend,1) = (b(1:iend,1) - A(1:iend,1,2)*x(1:iend,2) &
         &         - A(1:iend,1,3)*x(1:iend,3)) / A(1:iend,1,1)
    !    end associate

  end subroutine solve_vec_3

  pure function solve_vec_3_sw(ng_sw_in,A,b)
    integer,    intent(in)  :: ng_sw_in
    real(jprb), intent(in)  :: A(ng_sw,3,3)
    real(jprb), intent(in)  :: b(ng_sw,3)
    real(jprb) :: solve_vec_3_sw(ng_sw,3)

    real(jprb), dimension(ng_sw) :: L21, L31, L32
    real(jprb), dimension(ng_sw) :: U22, U23, U33
    real(jprb), dimension(ng_sw) :: y2, y3

    ! LU decomposition:
    !     ( 1        )   (U11 U12 U13)
    ! A = (L21  1    ) * (    U22 U23)
    !     (L31 L32  1)   (        U33)
    L21 = A(:,2,1) / A(:,1,1)
    L31 = A(:,3,1) / A(:,1,1)
    U22 = A(:,2,2) - L21*A(:,1,2)
    U23 = A(:,2,3) - L21*A(:,1,3)
    L32 =(A(:,3,2) - L31*A(:,1,2)) / U22
    U33 = A(:,3,3) - L31*A(:,1,3) - L32*U23

    ! Solve Ly = b by forward substitution
    y2 = b(:,2) - L21*b(:,1)
    y3 = b(:,3) - L31*b(:,1) - L32*y2

    ! Solve Ux = y by back substitution
    solve_vec_3_sw(:,3) = y3/U33
    solve_vec_3_sw(:,2) = (y2 - U23*solve_vec_3_sw(:,3)) / U22
    solve_vec_3_sw(:,1) = (b(:,1) - A(:,1,2)*solve_vec_3_sw(:,2) &
         &         - A(:,1,3)*solve_vec_3_sw(:,3)) / A(:,1,1)

  end function solve_vec_3_sw

  ! Like solve_vec_3 but overwrite b with x
  pure subroutine solve_vec_3_lw_inplace(ng_lw_in,A,b)
    integer,    intent(in)  :: ng_lw_in
    real(jprb), intent(in)  :: A(ng_lw,3,3)
    real(jprb), intent(inout)  :: b(ng_lw,3)
    !real(jprb), intent(out) :: x(:,:)

    real(jprb), dimension(ng_lw) :: L21, L31, L32
    real(jprb), dimension(ng_lw) :: U22, U23, U33
    real(jprb), dimension(ng_lw) :: y2, y3

    ! LU decomposition:
    !     ( 1        )   (U11 U12 U13)
    ! A = (L21  1    ) * (    U22 U23)
    !     (L31 L32  1)   (        U33)
    L21 = A(:,2,1) / A(:,1,1)
    L31 = A(:,3,1) / A(:,1,1)
    U22 = A(:,2,2) - L21*A(:,1,2)
    U23 = A(:,2,3) - L21*A(:,1,3)
    L32 =(A(:,3,2) - L31*A(:,1,2)) / U22
    U33 = A(:,3,3) - L31*A(:,1,3) - L32*U23

    ! Solve Ly = b by forward substitution
    y2 = b(:,2) - L21*b(:,1)
    y3 = b(:,3) - L31*b(:,1) - L32*y2

    ! Solve Ux = y by back substitution
    b(:,3) = y3/U33
    b(:,2) = (y2 - U23*b(:,3)) / U22
    b(:,1) = (b(:,1) - A(:,1,2)*b(:,2) &
         &         - A(:,1,3)*b(:,3)) / A(:,1,1)

  end subroutine solve_vec_3_lw_inplace

  pure function solve_vec_3_lw(ng_lw_in,A,b)
    integer,    intent(in)  :: ng_lw_in
    real(jprb), intent(in)  :: A(ng_lw,3,3)
    real(jprb), intent(in)  :: b(ng_lw,3)
    real(jprb) :: solve_vec_3_lw(ng_lw,3)

    real(jprb), dimension(ng_lw) :: L21, L31, L32
    real(jprb), dimension(ng_lw) :: U22, U23, U33
    real(jprb), dimension(ng_lw) :: y2, y3

    ! LU decomposition:
    !     ( 1        )   (U11 U12 U13)
    ! A = (L21  1    ) * (    U22 U23)
    !     (L31 L32  1)   (        U33)
    L21 = A(:,2,1) / A(:,1,1)
    L31 = A(:,3,1) / A(:,1,1)
    U22 = A(:,2,2) - L21*A(:,1,2)
    U23 = A(:,2,3) - L21*A(:,1,3)
    L32 =(A(:,3,2) - L31*A(:,1,2)) / U22
    U33 = A(:,3,3) - L31*A(:,1,3) - L32*U23

    ! Solve Ly = b by forward substitution
    y2 = b(:,2) - L21*b(:,1)
    y3 = b(:,3) - L31*b(:,1) - L32*y2

    ! Solve Ux = y by back substitution
    solve_vec_3_lw(:,3) = y3/U33
    solve_vec_3_lw(:,2) = (y2 - U23*solve_vec_3_lw(:,3)) / U22
    solve_vec_3_lw(:,1) = (b(:,1) - A(:,1,2)*solve_vec_3_lw(:,2) &
         &         - A(:,1,3)*solve_vec_3_lw(:,3)) / A(:,1,1)

  end function solve_vec_3_lw


  pure function solve_vec_3_ng(ng,A,b)
    integer,    intent(in)  :: ng
    real(jprb), intent(in)  :: A(ng,3,3)
    real(jprb), intent(in)  :: b(ng,3)
    real(jprb) :: solve_vec_3_ng(ng,3)

    real(jprb), dimension(ng) :: L21, L31, L32
    real(jprb), dimension(ng) :: U22, U23, U33
    real(jprb), dimension(ng) :: y2, y3

    ! LU decomposition:
    !     ( 1        )   (U11 U12 U13)
    ! A = (L21  1    ) * (    U22 U23)
    !     (L31 L32  1)   (        U33)
    L21 = A(:,2,1) / A(:,1,1)
    L31 = A(:,3,1) / A(:,1,1)
    U22 = A(:,2,2) - L21*A(:,1,2)
    U23 = A(:,2,3) - L21*A(:,1,3)
    L32 =(A(:,3,2) - L31*A(:,1,2)) / U22
    U33 = A(:,3,3) - L31*A(:,1,3) - L32*U23

    ! Solve Ly = b by forward substitution
    y2 = b(:,2) - L21*b(:,1)
    y3 = b(:,3) - L31*b(:,1) - L32*y2

    ! Solve Ux = y by back substitution
    solve_vec_3_ng(:,3) = y3/U33
    solve_vec_3_ng(:,2) = (y2 - U23*solve_vec_3_ng(:,3)) / U22
    solve_vec_3_ng(:,1) = (b(:,1) - A(:,1,2)*solve_vec_3_ng(:,2) &
         &         - A(:,1,3)*solve_vec_3_ng(:,3)) / A(:,1,1)

  end function solve_vec_3_ng

  !---------------------------------------------------------------------
  ! Return X = B A^-1 = (A^-T B)^T optimized for 3x3 matrices, where B
  ! is a diagonal matrix, using LU factorization and substitution with
  ! no pivoting.
  pure subroutine diag_mat_right_divide_3(n,iend,A,B,X)
    integer,    intent(in)  :: n, iend
    real(jprb), intent(in)  :: A(iend,3,3)
    real(jprb), intent(in)  :: B(iend,3)
    real(jprb), intent(out) :: X(n,3,3)

    real(jprb), dimension(iend) :: L21, L31, L32
    real(jprb), dimension(iend) :: U22, U23, U33
    real(jprb), dimension(iend) :: y2, y3

    !    associate (U11 => A(:,1,1), U12 => A(:,1,2), U13 => A(1,3))
    ! LU decomposition of the *transpose* of A:
    !       ( 1        )   (U11 U12 U13)
    ! A^T = (L21  1    ) * (    U22 U23)
    !       (L31 L32  1)   (        U33)
    L21 = A(1:iend,1,2) / A(1:iend,1,1)
    L31 = A(1:iend,1,3) / A(1:iend,1,1)
    U22 = A(1:iend,2,2) - L21*A(1:iend,2,1)
    U23 = A(1:iend,3,2) - L21*A(1:iend,3,1)
    L32 =(A(1:iend,2,3) - L31*A(1:iend,2,1)) / U22
    U33 = A(1:iend,3,3) - L31*A(1:iend,3,1) - L32*U23

    ! Solve X(1,:) = A^-T ( B(1) )
    !                     (  0   )
    !                     (  0   )
    ! Solve Ly = B(:,:,j) by forward substitution
    ! y1 = B(:,1)
    y2 = - L21*B(1:iend,1)
    y3 = - L31*B(1:iend,1) - L32*y2
    ! Solve UX(:,:,j) = y by back substitution
    X(1:iend,1,3) = y3 / U33
    X(1:iend,1,2) = (y2 - U23*X(1:iend,1,3)) / U22
    X(1:iend,1,1) = (B(1:iend,1) - A(1:iend,2,1)*X(1:iend,1,2) &
         &          - A(1:iend,3,1)*X(1:iend,1,3)) / A(1:iend,1,1)

    ! Solve X(2,:) = A^-T (  0   )
    !                     ( B(2) )
    !                     (  0   )
    ! Solve Ly = B(:,:,j) by forward substitution
    ! y1 = 0
    ! y2 = B(1:iend,2)
    y3 = - L32*B(1:iend,2)
    ! Solve UX(:,:,j) = y by back substitution
    X(1:iend,2,3) = y3 / U33
    X(1:iend,2,2) = (B(1:iend,2) - U23*X(1:iend,2,3)) / U22
    X(1:iend,2,1) = (-A(1:iend,2,1)*X(1:iend,2,2) &
         &           -A(1:iend,3,1)*X(1:iend,2,3)) / A(1:iend,1,1)

    ! Solve X(3,:) = A^-T (  0   )
    !                     (  0   )
    !                     ( B(3) )
    ! Solve Ly = B(:,:,j) by forward substitution
    ! y1 = 0
    ! y2 = 0
    ! y3 = B(1:iend,3)
    ! Solve UX(:,:,j) = y by back substitution
    X(1:iend,3,3) = B(1:iend,3) / U33
    X(1:iend,3,2) = -U23*X(1:iend,3,3) / U22
    X(1:iend,3,1) = (-A(1:iend,2,1)*X(1:iend,3,2) &
         &          - A(1:iend,3,1)*X(1:iend,3,3)) / A(1:iend,1,1)

  end subroutine diag_mat_right_divide_3


  !---------------------------------------------------------------------
  ! Treat A as n m-by-m matrices and return the LU factorization of A
  ! compressed into a single matrice (with L below the diagonal and U
  ! on and above the diagonal; the diagonal elements of L are 1). No
  ! pivoting is performed.
  pure subroutine lu_factorization(n, iend, m, A, LU)
    integer,    intent(in)  :: n, m, iend
    real(jprb), intent(in)  :: A(:,:,:)
    real(jprb), intent(out) :: LU(iend,m,m)

    real(jprb) :: s(iend)
    integer    :: j1, j2, j3

    ! This routine is adapted from an in-place one, so we first copy
    ! the input into the output.
    LU(1:iend,1:m,1:m) = A(1:iend,1:m,1:m)

    do j2 = 1, m
      do j1 = 1, j2-1
        s = LU(1:iend,j1,j2)
        do j3 = 1, j1-1
          s = s - LU(1:iend,j1,j3) * LU(1:iend,j3,j2)
        end do
        LU(1:iend,j1,j2) = s
      end do
      do j1 = j2, m
        s = LU(1:iend,j1,j2)
        do j3 = 1, j2-1
          s = s - LU(1:iend,j1,j3) * LU(1:iend,j3,j2)
        end do
        LU(1:iend,j1,j2) = s
      end do
      if (j2 /= m) then
        s = 1.0_jprb / LU(1:iend,j2,j2)
        do j1 = j2+1, m
          LU(1:iend,j1,j2) = LU(1:iend,j1,j2) * s
        end do
      end if
    end do

  end subroutine lu_factorization

  pure subroutine lu_factorization_lw(n, A, LU)
    integer,    intent(in)  :: n
    real(jprb), intent(in)  :: A(n,6,6)
    real(jprb), intent(out) :: LU(n,6,6)
    real(jprb) :: s(n)
    integer    :: j1, j2, j3

    ! This routine is adapted from an in-place one, so we first copy
    ! the input into the output.
    LU = A

    do j2 = 1, 6
      do j1 = 1, j2-1
        do j3 = 1, j1-1
          LU(:,j1,j2) = LU(:,j1,j2) - LU(:,j1,j3) * LU(:,j3,j2)
        end do
      end do
      do j1 = j2, 6
        do j3 = 1, j2-1
          LU(:,j1,j2) = LU(:,j1,j2) - LU(:,j1,j3) * LU(:,j3,j2)
        end do
      end do
      if (j2 /= 6) then
        s = 1.0_jprb / LU(:,j2,j2)
        do j1 = j2+1, 6
          LU(:,j1,j2) = LU(:,j1,j2) * s
        end do
      end if
    end do

  end subroutine lu_factorization_lw

  !---------------------------------------------------------------------
  ! Treat LU as an LU-factorization of an original matrix A, and
  ! return x where Ax=b. LU consists of n m-by-m matrices and b as n
  ! m-element vectors.
  pure subroutine lu_substitution(n,iend,m,LU,b,x)
    ! CHECK: dimensions should be ":"?
    integer,    intent(in) :: n, m, iend
    real(jprb), intent(in) :: LU(iend,m,m)
    real(jprb), intent(in) :: b(:,:)
    real(jprb), intent(out):: x(iend,m)

    integer :: j1, j2

    x(1:iend,1:m) = b(1:iend,1:m)

    ! First solve Ly=b
    do j2 = 2, m
      do j1 = 1, j2-1
        x(1:iend,j2) = x(1:iend,j2) - x(1:iend,j1)*LU(1:iend,j2,j1)
      end do
    end do
    ! Now solve Ux=y
    do j2 = m, 1, -1
      do j1 = j2+1, m
        x(1:iend,j2) = x(1:iend,j2) - x(1:iend,j1)*LU(1:iend,j2,j1)
      end do
      x(1:iend,j2) = x(1:iend,j2) / LU(1:iend,j2,j2)
    end do

  end subroutine lu_substitution

  ! Optimized version, overwrite input vector b with output vector x
  pure subroutine lu_substitution_lw_inplace(n,LU,b) ! (n,LU,b,x)
    integer,    intent(in) :: n
    real(jprb), intent(in) :: LU(n,6,6)
    real(jprb), intent(inout) :: b(n,6)
    ! real(jprb), intent(out):: x(n,6)
    integer :: j1, j2
    ! x = b

    ! First solve Ly=b
    ! do j2 = 2, 6
    !   do j1 = 1, j2-1
    !     b(:,j2) = b(:,j2) - b(:,j1)*LU(:,j2,j1)
    !   end do
    ! end do
    ! j2=2 => j1=1,1
    b(:,2) = b(:,2) - b(:,1)*LU(:,2,1)
    ! j2=3 => j1=1,2
    b(:,3) = b(:,3) - b(:,1)*LU(:,3,1) - b(:,2)*LU(:,3,2) 
    ! j2=4 => j1=1,3
    b(:,4) = b(:,4) - b(:,1)*LU(:,4,1) - b(:,2)*LU(:,4,2) - b(:,3)*LU(:,4,3) 
    ! j2=5 => j1=1,4
    b(:,5) = b(:,5) - b(:,1)*LU(:,5,1) - b(:,2)*LU(:,5,2) - b(:,3)*LU(:,5,3) - b(:,4)*LU(:,5,4) 
    ! j2=6 => j1=1,5
    b(:,6) = b(:,6) - b(:,1)*LU(:,6,1) - b(:,2)*LU(:,6,2) - b(:,3)*LU(:,6,3) - b(:,4)*LU(:,6,4) &
        & - b(:,5)*LU(:,6,5)
    ! Now solve Ux=y
    ! do j2 = 6, 1, -1
    !   do j1 = j2+1, 6
    !     b(:,j2) = b(:,j2) - b(:,j1)*LU(:,j2,j1)
    !   end do
    !   b(:,j2) = b(:,j2) / LU(:,j2,j2)
    ! end do
    ! j2=6 => j1=7,6
    b(:,6) = b(:,6) / LU(:,6,6)
    ! j2=5 => j1=6,6
    b(:,5) = (b(:,5) - b(:,6)*LU(:,5,6)) / LU(:,5,5)
    ! j2=4 => j1=5,6
    b(:,4) = ( b(:,4) - b(:,5)*LU(:,4,5) - b(:,6)*LU(:,4,6)) / LU(:,4,4)
    ! j2=3 => j1=4,6
    b(:,3) = (b(:,3) - b(:,4)*LU(:,3,4) - b(:,5)*LU(:,3,5) - b(:,6)*LU(:,3,6)) / LU(:,3,3)
    ! j2=2 => j1=3,6
    b(:,2) = ( b(:,2) - b(:,3)*LU(:,2,3) - b(:,4)*LU(:,2,4) - b(:,5)*LU(:,2,5) - b(:,6)*LU(:,2,6))  &
        &   / LU(:,2,2)
    ! j2=1 => j1=2,6
    b(:,1) = ( b(:,1) - b(:,2)*LU(:,1,2) - b(:,3)*LU(:,1,3) - b(:,4)*LU(:,1,4) - b(:,5)*LU(:,1,5) &
        &   - b(:,6)*LU(:,1,6) ) / LU(:,1,1)

  end subroutine lu_substitution_lw_inplace


  !---------------------------------------------------------------------
  ! Return matrix X where AX=B. LU, A, X, B all consist of n m-by-m
  ! matrices.
  pure subroutine solve_mat_n(n,iend,m,A,B,X)
    integer,    intent(in) :: n, m, iend
    real(jprb), intent(in) :: A(:,:,:)
    real(jprb), intent(in) :: B(:,:,:)
    real(jprb), intent(out):: X(iend,m,m)

    real(jprb) :: LU(iend,m,m)

    integer :: j

    call lu_factorization(n,iend,m,A,LU)

    do j = 1, m
      call lu_substitution(n,iend,m,LU,B(1:,1:m,j),X(1:iend,1:m,j))
!      call lu_substitution(n,iend,m,LU,B(1:n,1:m,j),X(1:iend,1:m,j))
    end do

  end subroutine solve_mat_n

  ! Solve AX=B, where A, X and B consist of ng m-by-m matrices
  ! To minimize memory use, overwrite B with X, and A is corrupted
  pure subroutine solve_mat_inplace(n,m,A,B)
    integer,    intent(in) :: n, m
    real(jprb), intent(inout) :: A(n,m,m)
    real(jprb), intent(inout) :: B(n,m,m)

    real(jprb) :: s(n)

    integer :: j1, j2, j3

    ! LU factorization
    do j2 = 1, m
      do j1 = 1, j2-1
        s = A(:,j1,j2)
        do j3 = 1, j1-1
          s = s - A(:,j1,j3) * A(:,j3,j2)
        end do
        A(:,j1,j2) = s
      end do
      do j1 = j2, m
        s = A(:,j1,j2)
        do j3 = 1, j2-1
          s = s - A(:,j1,j3) * A(:,j3,j2)
        end do
        A(:,j1,j2) = s
      end do
      if (j2 /= m) then
        s = 1.0_jprb / A(:,j2,j2)
        do j1 = j2+1, m
          A(:,j1,j2) = A(:,j1,j2) * s
        end do
      end if
    end do

    ! LU substitution
    do j3 = 1, m
      ! First solve Ly=b
      do j2 = 2, m
        do j1 = 1, j2-1
          B(:,j2,j3) = B(:,j2,j3) - B(:,j1,j3)*A(:,j2,j1)
        end do
      end do
      ! Now solve Ux=y  (UB=y)
      do j2 = m, 1, -1
        do j1 = j2+1, m
          B(:,j2,j3) = B(:,j2,j3) - B(:,j1,j3)*A(:,j2,j1)
        end do
        B(:,j2,j3) = B(:,j2,j3) / A(:,j2,j2)
      end do
    end do

  end subroutine solve_mat_inplace

  pure subroutine solve_mat_3_inplace(ng,A,B)
    integer,    intent(in)  :: ng
    real(jprb), intent(in)  :: A(ng,3,3)
    real(jprb), intent(inout)  :: B(ng,3,3)
    ! real(jprb), intent(out) :: X(ng_lw,3,3)

    real(jprb), dimension(ng) :: L21, L31, L32
    real(jprb), dimension(ng) :: U22, U23, U33
    real(jprb), dimension(ng) :: y2, y3

    integer :: j
    associate(X=>B)

      !    associate (U11 => A(:,1,1), U12 => A(:,1,2), U13 => A(1,3))
      ! LU decomposition:
      !     ( 1        )   (U11 U12 U13)
      ! A = (L21  1    ) * (    U22 U23)
      !     (L31 L32  1)   (        U33)
      L21 = A(:,2,1) / A(:,1,1)
      L31 = A(:,3,1) / A(:,1,1)
      U22 = A(:,2,2) - L21*A(:,1,2)
      U23 = A(:,2,3) - L21*A(:,1,3)
      L32 =(A(:,3,2) - L31*A(:,1,2)) / U22
      U33 = A(:,3,3) - L31*A(:,1,3) - L32*U23

      do j = 1,3
        ! Solve Ly = B(:,:,j) by forward substitution
        ! y1 = B(:,1,j)
        y2 = B(:,2,j) - L21*B(:,1,j)
        y3 = B(:,3,j) - L31*B(:,1,j) - L32*y2
        ! Solve UX(:,:,j) = y by back substitution
        X(:,3,j) = y3 / U33
        X(:,2,j) = (y2 - U23*X(:,3,j)) / U22
        X(:,1,j) = (B(:,1,j) - A(:,1,2)*X(:,2,j) &
            &          - A(:,1,3)*X(:,3,j)) / A(:,1,1)
      end do
    end associate 
  end subroutine solve_mat_3_inplace

    !---------------------------------------------------------------------
  ! Solve AX=B optimized for 3x3 matrices, using LU factorization and
  ! substitution with no pivoting.
  pure subroutine solve_mat_3(n,iend,A,B,X)
    integer,    intent(in)  :: n, iend
    real(jprb), intent(in)  :: A(:,:,:)
    real(jprb), intent(in)  :: B(:,:,:)
    real(jprb), intent(out) :: X(:,:,:)

    real(jprb), dimension(iend) :: L21, L31, L32
    real(jprb), dimension(iend) :: U22, U23, U33
    real(jprb), dimension(iend) :: y2, y3

    integer :: j

    !    associate (U11 => A(:,1,1), U12 => A(:,1,2), U13 => A(1,3))
    ! LU decomposition:
    !     ( 1        )   (U11 U12 U13)
    ! A = (L21  1    ) * (    U22 U23)
    !     (L31 L32  1)   (        U33)
    L21 = A(1:iend,2,1) / A(1:iend,1,1)
    L31 = A(1:iend,3,1) / A(1:iend,1,1)
    U22 = A(1:iend,2,2) - L21*A(1:iend,1,2)
    U23 = A(1:iend,2,3) - L21*A(1:iend,1,3)
    L32 =(A(1:iend,3,2) - L31*A(1:iend,1,2)) / U22
    U33 = A(1:iend,3,3) - L31*A(1:iend,1,3) - L32*U23

    do j = 1,3
      ! Solve Ly = B(:,:,j) by forward substitution
      ! y1 = B(:,1,j)
      y2 = B(1:iend,2,j) - L21*B(1:iend,1,j)
      y3 = B(1:iend,3,j) - L31*B(1:iend,1,j) - L32*y2
      ! Solve UX(:,:,j) = y by back substitution
      X(1:iend,3,j) = y3 / U33
      X(1:iend,2,j) = (y2 - U23*X(1:iend,3,j)) / U22
      X(1:iend,1,j) = (B(1:iend,1,j) - A(1:iend,1,2)*X(1:iend,2,j) &
           &          - A(1:iend,1,3)*X(1:iend,3,j)) / A(1:iend,1,1)
    end do

  end subroutine solve_mat_3

  pure subroutine solve_mat_3_sw(ng_sw_in,A,B,X)
    integer,    intent(in)  :: ng_sw_in
    real(jprb), intent(in)  :: A(ng_sw,3,3)
    real(jprb), intent(in)  :: B(ng_sw,3,3)
    real(jprb), intent(out) :: X(ng_sw,3,3)

    real(jprb), dimension(ng_sw) :: L21, L31, L32
    real(jprb), dimension(ng_sw) :: U22, U23, U33
    real(jprb), dimension(ng_sw) :: one_over_A11
    real(jprb) :: y2, y3
    integer :: j, jg

    ! LU decomposition:
    !     ( 1        )   (U11 U12 U13)
    ! A = (L21  1    ) * (    U22 U23)
    !     (L31 L32  1)   (        U33)
    do jg = 1, ng_sw
      one_over_A11(jg) = 1 / A(jg,1,1)
      L21(jg) = A(jg,2,1) * one_over_A11(jg) !/ A(jg,1,1)
      L31(jg) = A(jg,3,1) * one_over_A11(jg) !/ A(jg,1,1)
      U22(jg) = A(jg,2,2) - L21(jg)*A(jg,1,2)
      U23(jg) = A(jg,2,3) - L21(jg)*A(jg,1,3)
      L32(jg) =(A(jg,3,2) - L31(jg)*A(jg,1,2)) / U22(jg)
      U33(jg) = A(jg,3,3) - L31(jg)*A(jg,1,3) - L32(jg)*U23(jg)
    end do 
    do j = 1,3
      do jg = 1, ng_sw
        ! Solve Ly = B(:,:,j) by forward substitution
        ! y1 = B(jg,1,j)
        y2 = B(jg,2,j) - L21(jg)*B(jg,1,j)
        y3 = B(jg,3,j) - L31(jg)*B(jg,1,j) - L32(jg)*y2
        ! Solve UX(:,:,j) = y by back substitution
        X(jg,3,j) = y3 / U33(jg)
        X(jg,2,j) = (y2 - U23(jg)*X(jg,3,j)) / U22(jg)
        X(jg,1,j) = (B(jg,1,j) - A(jg,1,2)*X(jg,2,j) &
            &          - A(jg,1,3)*X(jg,3,j)) * one_over_A11(jg) !/ A(:,1,1)
      end do 
    end do

  end subroutine solve_mat_3_sw

  pure subroutine solve_mat_3_lw(ng_lw_in,A,B,X)
    integer,    intent(in)  :: ng_lw_in
    real(jprb), intent(in)  :: A(ng_lw,3,3)
    real(jprb), intent(in)  :: B(ng_lw,3,3)
    real(jprb), intent(out) :: X(ng_lw,3,3)

    real(jprb), dimension(ng_lw) :: L21, L31, L32
    real(jprb), dimension(ng_lw) :: U22, U23, U33
    real(jprb), dimension(ng_lw) :: y2, y3, one_over_A11

    integer :: j

    ! LU decomposition:
    !     ( 1        )   (U11 U12 U13)
    ! A = (L21  1    ) * (    U22 U23)
    !     (L31 L32  1)   (        U33)
    one_over_A11 = 1 / A(:,1,1)
    L21 = A(:,2,1) * one_over_A11 !/ A(:,1,1)
    L31 = A(:,3,1) * one_over_A11 !/ A(:,1,1)
    U22 = A(:,2,2) - L21*A(:,1,2)
    U23 = A(:,2,3) - L21*A(:,1,3)
    L32 =(A(:,3,2) - L31*A(:,1,2)) / U22
    U33 = A(:,3,3) - L31*A(:,1,3) - L32*U23

    do j = 1,3
      ! Solve Ly = B(:,:,j) by forward substitution
      ! y1 = B(:,1,j)
      y2 = B(:,2,j) - L21*B(:,1,j)
      y3 = B(:,3,j) - L31*B(:,1,j) - L32*y2
      ! Solve UX(:,:,j) = y by back substitution
      X(:,3,j) = y3 / U33
      X(:,2,j) = (y2 - U23*X(:,3,j)) / U22
      X(:,1,j) = (B(:,1,j) - A(:,1,2)*X(:,2,j) &
           &          - A(:,1,3)*X(:,3,j)) * one_over_A11 !/ A(:,1,1)
    end do

  end subroutine solve_mat_3_lw

  !---------------------------------------------------------------------
  ! Solve Ax=b, where A consists of n m-by-m matrices and x and b
  ! consist of n m-element vectors. For m=2 or m=3, this function
  ! calls optimized versions, otherwise it uses general LU
  ! decomposition without pivoting.
  function solve_vec(n,iend,m,A,b)

    use yomhook, only : lhook, dr_hook, jphook

    integer,    intent(in)  :: n, m, iend
    real(jprb), intent(in)  :: A(:,:,:)
    real(jprb), intent(in)  :: b(:,:)

    real(jprb)              :: solve_vec(iend,m)
    real(jprb)              :: LU(iend,m,m)
    real(jphook)            :: hook_handle

    if (lhook) call dr_hook('radiation_matrix:solve_vec',0,hook_handle)

    if (m == 2) then
      call solve_vec_2(n,iend,A,b,solve_vec)
    elseif (m == 3) then
      call solve_vec_3(n,iend,A,b,solve_vec)
    else
      call lu_factorization(n,iend,m,A,LU)
      call lu_substitution(n,iend,m,LU,b,solve_vec)
    end if

    if (lhook) call dr_hook('radiation_matrix:solve_vec',1,hook_handle)

  end function solve_vec

  !---------------------------------------------------------------------
  ! Solve AX=B, where A, X and B consist of n m-by-m matrices. For m=2
  ! or m=3, this function calls optimized versions, otherwise it uses
  ! general LU decomposition without pivoting.
  function solve_mat(n,iend,m,A,B)

    use yomhook, only : lhook, dr_hook, jphook

    integer,    intent(in)  :: n, m, iend
    real(jprb), intent(in)  :: A(:,:,:)
    real(jprb), intent(in)  :: B(:,:,:)

    real(jprb)              :: solve_mat(iend,m,m)
    real(jphook)            :: hook_handle

    if (lhook) call dr_hook('radiation_matrix:solve_mat',0,hook_handle)

    if (m == 2) then
      call solve_mat_2(n,iend,A,B,solve_mat)
    elseif (m == 3) then
      call solve_mat_3(n,iend,A,B,solve_mat)
    else
      call solve_mat_n(n,iend,m,A,B,solve_mat)
    end if

    if (lhook) call dr_hook('radiation_matrix:solve_mat',1,hook_handle)

  end function solve_mat

  pure subroutine mat_square_sw(ng_sw_in,nlev_b,A,C)
    ! Square 9x9 matrices, assuming SW pattern with sparsity
    !     (C D E)
    ! A = (F G H)
    !     (0 0 I)
    integer,    intent(in)                                :: ng_sw_in, nlev_b
    real(jprb), intent(in),  dimension(ng_sw*nlev_b,9,9)  :: A
    !dir$ assume_aligned A:64
    real(jprb), intent(out), dimension(ng_sw*nlev_b,9,9)  :: C
    !dir$ assume_aligned C:64
    integer    :: j1, j2, j3
    
    ! Do the top-left (C, D, F, G)
    do j2 = 1,6  
      do j1 = 1,6 
        ! do j3 = 1,6 
        !   C(:,j1,j2) = C(:,j1,j2) + A(:,j1,j3)*A(:,j3,j2)
        ! end do
        ! flatten last loop, only one write SIMD instruction
        C(:,j1,j2) = A(:,j1,1)*A(:,1,j2) + A(:,j1,2)*A(:,2,j2) &
        &          + A(:,j1,3)*A(:,3,j2) + A(:,j1,4)*A(:,4,j2) &
        &          + A(:,j1,5)*A(:,5,j2) + A(:,j1,6)*A(:,6,j2)
        ! C(:,j1,j2) = sum(A(:,j1,1:6)*A(:,1:6,j2),2)
      end do
    end do

    do j2 = 7,9
      C(:,:,j2) = 0.0_jprb
      ! Do the top-right (E & H)
      do j1 =  1,6
        do j3 = 1,9
          C(:,j1,j2) = C(:,j1,j2) + A(:,j1,j3)*A(:,j3,j2)
        end do
      end do

      ! Do the bottom-right (I)
      do j1 = 7,9
        do j3 = 7,9
          C(:,j1,j2) =  C(:,j1,j2) + A(:,j1,j3)*A(:,j3,j2)
        end do
      end do
    end do

  end subroutine mat_square_sw

  pure subroutine mat_x_mat_sw(ng_sw_in,nlev_b,A,B,C)

    integer,    intent(in)                    :: ng_sw_in, nlev_b
    real(jprb), intent(in), dimension(ng_sw*nlev_b,9,9) :: A, B
    !dir$ assume_aligned A:64,B:64
    real(jprb), intent(out),dimension(ng_sw*nlev_b,9,9) :: C
    !dir$ assume_aligned C:64
    integer    :: j1, j2, j3
    !     (C    D    E)
    ! A = (F=   G    H)
    !     (0    0    I)

    ! Do the top-left (C, D, F, G)
    do j2 = 1,6
      do j1 = 1,6
        C(:,j1,j2) =  A(:,j1,1)*B(:,1,j2) + A(:,j1,2)*B(:,2,j2) + A(:,j1,3)*B(:,3,j2) &
             &      + A(:,j1,4)*B(:,4,j2) + A(:,j1,5)*B(:,5,j2) + A(:,j1,6)*B(:,6,j2)
      end do
    end do
    do j2 = 7,9
      ! Do the top-right (E & H)
      do j1 = 1,6
        C(:,j1,j2) =  A(:,j1,1)*B(:,1,j2) + A(:,j1,2)*B(:,2,j2) + A(:,j1,3)*B(:,3,j2) &
             &      + A(:,j1,4)*B(:,4,j2) + A(:,j1,5)*B(:,5,j2) + A(:,j1,6)*B(:,6,j2) &
             &      + A(:,j1,7)*B(:,7,j2) + A(:,j1,8)*B(:,8,j2) + A(:,j1,9)*B(:,9,j2)
      end do
      ! Do the bottom-right (I)
      do j1 = 7,9
        C(:,j1,j2) = A(:,j1,7)*B(:,7,j2) + A(:,j1,8)*B(:,8,j2) + A(:,j1,9)*B(:,9,j2)
      end do
    end do
    ! Lower left corner with zeros
    C(:,7:9,1:6) = 0.0_jprb

  end subroutine mat_x_mat_sw

  pure subroutine mat_x_mat_sw_repeats(ng_sw_in, nlev_b, A, B, C)
    integer,    intent(in)                              :: ng_sw_in, nlev_b
    real(jprb), intent(in), dimension(ng_sw*nlev_b,9,9) :: A, B
    !dir$ assume_aligned A:64,B:64
    real(jprb), intent(out),dimension(ng_sw*nlev_b,9,9) :: C
    !dir$ assume_aligned C:64
    integer    :: j1, j2, j22
    ! Input matrices have pattern
    !     (C    D     E)
    ! A = (F=-D G=-C  H)
    !     (0    0     I)
    ! As a result, output matrices has pattern
    !     (C    D    E)
    ! C = (F=D  G=C  H)
    !     (0    0    I)
    do j2 = 1,3
      j22 = j2 + 6
      do j1 = 1,6 
        ! Do the top-left (C, F)
        C(:,j1,j2) = A(:,j1,1)*B(:,1,j2) + A(:,j1,2)*B(:,2,j2) + A(:,j1,3)*B(:,3,j2) &
        &          + A(:,j1,4)*B(:,4,j2) + A(:,j1,6)*B(:,6,j2)
        ! Do the top-right (E & H)
        C(:,j1,j22) = A(:,j1,1)*B(:,1,j22) + A(:,j1,2)*B(:,2,j22) + A(:,j1,3)*B(:,3,j22) &
        &      + A(:,j1,4)*B(:,4,j22) + A(:,j1,5)*B(:,5,j22) + A(:,j1,6)*B(:,6,j22) &
        &      + A(:,j1,7)*B(:,7,j22) + A(:,j1,8)*B(:,8,j22) + A(:,j1,9)*B(:,9,j22)
      end do
      do j1 = 7,9  ! Do the bottom-right (I)
        C(:,j1,j22) = A(:,j1,7)*B(:,7,j22) + A(:,j1,8)*B(:,8,j22) + A(:,j1,9)*B(:,9,j22)
      end do
    end do
    C(:,1:3,4:6) = C(:,4:6,1:3) ! D = F
    C(:,4:6,4:6) = C(:,1:3,1:3) ! G = C
    C(:,7:9,1:6) = 0.0_jprb     ! Lower left corner
  
  end subroutine mat_x_mat_sw_repeats

  pure subroutine mat_x_mat_sw_repeatss(n, A, B, C)
    integer,    intent(in)                              :: n
    real(jprb), intent(in), dimension(n,9,9) :: A, B
    !dir$ assume_aligned A:64,B:64
    real(jprb), intent(out),dimension(n,9,9) :: C
    !dir$ assume_aligned C:64
    integer    :: j1, j2, j22
    ! Input matrices have pattern
    !     (C    D     E)
    ! A = (F=-D G=-C  H)
    !     (0    0     I)
    ! As a result, output matrices has pattern
    !     (C    D    E)
    ! C = (F=D  G=C  H)
    !     (0    0    I)
    do j2 = 1,3
      j22 = j2 + 6
      do j1 = 1,6 
        ! Do the top-left (C, F)
        C(:,j1,j2) = A(:,j1,1)*B(:,1,j2) + A(:,j1,2)*B(:,2,j2) + A(:,j1,3)*B(:,3,j2) &
        &          + A(:,j1,4)*B(:,4,j2) + A(:,j1,6)*B(:,6,j2)
        ! Do the top-right (E & H)
        C(:,j1,j22) = A(:,j1,1)*B(:,1,j22) + A(:,j1,2)*B(:,2,j22) + A(:,j1,3)*B(:,3,j22) &
        &      + A(:,j1,4)*B(:,4,j22) + A(:,j1,5)*B(:,5,j22) + A(:,j1,6)*B(:,6,j22) &
        &      + A(:,j1,7)*B(:,7,j22) + A(:,j1,8)*B(:,8,j22) + A(:,j1,9)*B(:,9,j22)
      end do
      do j1 = 7,9  ! Do the bottom-right (I)
        C(:,j1,j22) = A(:,j1,7)*B(:,7,j22) + A(:,j1,8)*B(:,8,j22) + A(:,j1,9)*B(:,9,j22)
      end do
    end do
    C(:,1:3,4:6) = C(:,4:6,1:3) ! D = F
    C(:,4:6,4:6) = C(:,1:3,1:3) ! G = C
    C(:,7:9,1:6) = 0.0_jprb     ! Lower left corner
  
  end subroutine mat_x_mat_sw_repeatss

  !---------------------------------------------------------------------
  ! Solve AX=B, where A, X and B consist of ng m-by-m matrices
  ! Overwrite B with X. A is corrupted
  pure subroutine solve_mat_sw(ng_sw_in,nlev_b,A,B)
    integer,    intent(in)    :: ng_sw_in, nlev_b
    real(jprb), intent(inout) :: A(ng_sw*nlev_b,9,9) ! A=LU is corrupted
    !dir$ assume_aligned A:64
    real(jprb), intent(inout) :: B(ng_sw*nlev_b,9,9) ! X = B, both input and output
    !dir$ assume_aligned B:64
    ! Local variables
    integer :: j1, j2, j3
    integer, parameter :: m = 9
    integer, parameter :: mblock = 3  ! m/3   
    integer, parameter :: m2block = 6 ! 2*mblock 
    real(jprb) :: diagdiv(ng_sw*nlev_b,9)
    !     (C   D  E)
    ! A = (-D -C  H)
    !     (0   0  I)
 
    ! factorization of A into LU
    ! First do columns 1-6, for which only rows 1-6 have non-negative entries
    do j2 = 1, m2block
      do j1 = 1, j2-1
        do j3 = 1, j1-1
          A(:,j1,j2) = A(:,j1,j2)- A(:,j1,j3) * A(:,j3,j2)
        end do
      end do
      do j1 = j2, m2block
        do j3 = 1, j2-1
          A(:,j1,j2) = A(:,j1,j2) - A(:,j1,j3) * A(:,j3,j2)
        end do
      end do
      ! s = 1.0_jprb / A(:,j2,j2)
      diagdiv(:,j2) =  1.0_jprb / A(:,j2,j2)
      do j1 = j2+1, m2block
        A(:,j1,j2) = A(:,j1,j2) * diagdiv(:,j2) ! * s
      end do
    end do
    ! Remaining columns
    do j2 = m2block+1, m
      do j1 = 1, j2-1
        do j3 = 1, j1-1
          A(:,j1,j2) = A(:,j1,j2) - A(:,j1,j3) * A(:,j3,j2)
        end do
      end do
      do j1 = j2, m
        do j3 = 1, j2-1
          A(:,j1,j2)= A(:,j1,j2) - A(:,j1,j3) * A(:,j3,j2)
        end do
      end do
      diagdiv(:,j2) =  1.0_jprb / A(:,j2,j2)
      if (j2 /= m) then
        ! s = 1.0_jprb / A(:,j2,j2)
        do j1 = j2+1, m
          A(:,j1,j2) = A(:,j1,j2) * diagdiv(:,j2) ! * s
        end do
      end if
    end do
    !---------------------------------------------------------------------
    ! Treat LU as an LU-factorization of an original matrix A, and
    ! return x where Ax=b. LU consists of n m-by-m matrices and b as n
    ! m-element vectors.
    ! Here B is both input b and output x, and A has been LU factorized, combining L and U
    ! into one matrix where the diagonal is the diagonal of U (L has ones in the diagonal)

    ! A and B both have following structure:
    !     (C   D  E)
    !     (F   G  H)
    !     (0   0  I)
    ! Separate j3 (columns) into two regions to avoid redundant operations with zero
    do j3 = 1,m2block ! in this region B(:,7:9),A(:,7:9) are 0
      ! First solve Ly=b
      ! do j2 = 2, m2block
      !   do j1 = 1, j2-1
      !     B(:,j2,j3) = B(:,j2,j3) - B(:,j1,j3)*A(:,j2,j1)
      !   end do
      !   ! No division because diagonal of L is unity
      ! end do
      B(:,2,j3) = B(:,2,j3) - B(:,1,j3)*A(:,2,1)
      B(:,3,j3) = B(:,3,j3) - B(:,1,j3)*A(:,3,1) - B(:,2,j3)*A(:,3,2) 
      B(:,4,j3) = B(:,4,j3) - B(:,1,j3)*A(:,4,1) - B(:,2,j3)*A(:,4,2) - B(:,3,j3)*A(:,4,3) 
      B(:,5,j3) = B(:,5,j3) - B(:,1,j3)*A(:,5,1) - B(:,2,j3)*A(:,5,2) - B(:,3,j3)*A(:,5,3) - B(:,4,j3)*A(:,5,4)
      B(:,6,j3) = B(:,6,j3) - B(:,1,j3)*A(:,6,1) - B(:,2,j3)*A(:,6,2) - B(:,3,j3)*A(:,6,3) &
        & - B(:,4,j3)*A(:,6,4) - B(:,5,j3)*A(:,6,5) 
      ! Now solve Ux=y
      ! do j2 = m2block, 1, -1
      !   do j1 = j2+1, m2block
      !     B(:,j2,j3) = ( B(:,j2,j3) - B(:,j1,j3)*A(:,j2,j1) )
      !   end do
      !   B(:,j2,j3) = B(:,j2,j3) * diagdiv(:,j2) ! / A(:,j2,j2) ! Divide by diagonal of A=U
      ! end do
      B(:,6,j3) = B(:,6,j3) * diagdiv(:,6)
      B(:,5,j3) = B(:,5,j3) - B(:,6,j3)*A(:,5,6)
      B(:,5,j3) = B(:,5,j3) * diagdiv(:,5)
      B(:,4,j3) = B(:,4,j3) - B(:,5,j3)*A(:,4,5) - B(:,6,j3)*A(:,4,6)
      B(:,4,j3) = B(:,4,j3) * diagdiv(:,4)
      B(:,3,j3) = B(:,3,j3) - B(:,4,j3)*A(:,3,4) - B(:,5,j3)*A(:,3,5) - B(:,6,j3)*A(:,3,6)
      B(:,3,j3) = B(:,3,j3) * diagdiv(:,3)
      B(:,2,j3) = B(:,2,j3) - B(:,3,j3)*A(:,2,3) - B(:,4,j3)*A(:,2,4) - B(:,5,j3)*A(:,2,5) - B(:,6,j3)*A(:,2,6)
      B(:,2,j3) = B(:,2,j3) * diagdiv(:,2)
      B(:,1,j3) = B(:,1,j3) - B(:,2,j3)*A(:,1,2) - B(:,3,j3)*A(:,1,3) - B(:,4,j3)*A(:,1,4) &
        & - B(:,5,j3)*A(:,1,5) - B(:,6,j3)*A(:,1,6) 
      B(:,1,j3) = B(:,1,j3) * diagdiv(:,1) 
    end do

    do j3 = m2block+1,m ! columns 7-9: here B has nonzero values for all rows, but A doesn't 
      ! First solve Ly=b
      ! do j2 = 2, m2block 
      !   do j1 = 1, j2-1
      !     B(:,j2,j3) = B(:,j2,j3) - B(:,j1,j3)*A(:,j2,j1)
      !   end do
      !   ! No division because diagonal of L is unity
      ! end do
      B(:,2,j3) = B(:,2,j3) - B(:,1,j3)*A(:,2,1)
      B(:,3,j3) = B(:,3,j3) - B(:,1,j3)*A(:,3,1) - B(:,2,j3)*A(:,3,2) 
      B(:,4,j3) = B(:,4,j3) - B(:,1,j3)*A(:,4,1) - B(:,2,j3)*A(:,4,2) - B(:,3,j3)*A(:,4,3) 
      B(:,5,j3) = B(:,5,j3) - B(:,1,j3)*A(:,5,1) - B(:,2,j3)*A(:,5,2) - B(:,3,j3)*A(:,5,3) - B(:,4,j3)*A(:,5,4)
      B(:,6,j3) = B(:,6,j3) - B(:,1,j3)*A(:,6,1) - B(:,2,j3)*A(:,6,2) - B(:,3,j3)*A(:,6,3) &
        & - B(:,4,j3)*A(:,6,4) - B(:,5,j3)*A(:,6,5) 
      ! When j2 = 7, the A terms are all 0, because A(7,1:6)=0
      ! when j2 = 8, only the last j1 has nonzero A
      ! When j2 = 9, two last j1 are nonzero
      B(:,8,j3) = B(:,8,j3) - B(:,7,j3)*A(:,8,7)   ! j2 = 8
      B(:,9,j3) = B(:,9,j3) - B(:,8,j3)*A(:,9,8) - B(:,7,j3)*A(:,9,7) ! j2 = 9
      ! Now solve Ux=y
      do j2 = m, 1, -1
        do j1 = j2+1, m
          B(:,j2,j3) = B(:,j2,j3) - B(:,j1,j3)*A(:,j2,j1)
        end do
        B(:,j2,j3) = B(:,j2,j3) * diagdiv(:,j2) ! / A(:,j2,j2) ! Divide by diagonal of A=U
      end do
    end do
  end subroutine solve_mat_sw

  ! --- MATRIX EXPONENTIATION ---

  !---------------------------------------------------------------------
  ! Perform matrix exponential of n m-by-m matrices stored in A (where
  ! index n varies fastest) using the Higham scaling and squaring
  ! method. The result is placed in A. This routine is intended for
  ! speed so is accurate only to single precision.  For simplicity and
  ! to aid vectorization, the Pade approximant of order 7 is used for
  ! all input matrices, perhaps leading to a few too many
  ! multiplications for matrices with a small norm.
  subroutine expm(n,iend,m,A,i_matrix_pattern)

    use yomhook, only : lhook, dr_hook, jphook

    integer,    intent(in)      :: n, m, iend
    real(jprb), intent(inout)   :: A(n,m,m)
    integer,    intent(in)      :: i_matrix_pattern

    real(jprb), parameter :: theta(3) = (/4.258730016922831e-01_jprb, &
         &                                1.880152677804762e+00_jprb, &
         &                                3.925724783138660e+00_jprb/) 
    real(jprb), parameter :: c(8) = (/17297280.0_jprb, 8648640.0_jprb, &
         &                1995840.0_jprb, 277200.0_jprb, 25200.0_jprb, &
         &                1512.0_jprb, 56.0_jprb, 1.0_jprb/)

    real(jprb), dimension(iend,m,m) :: A2, A4, A6
    real(jprb), dimension(iend,m,m) :: U, V

    real(jprb) :: normA(iend), sum_column(iend)

    integer    :: j1, j2, j3
    real(jprb) :: frac(iend)
    integer    :: expo(iend)
    real(jprb) :: scaling(iend)

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_matrix:expm',0,hook_handle)

    normA = 0.0_jprb

    ! Compute the 1-norms of A
    do j3 = 1,m
      sum_column(:) = 0.0_jprb
      do j2 = 1,m
        do j1 = 1,iend
          sum_column(j1) = sum_column(j1) + abs(A(j1,j2,j3))
        end do
      end do
      do j1 = 1,iend
        if (sum_column(j1) > normA(j1)) then
          normA(j1) = sum_column(j1)
        end if
      end do
    end do

    frac = fraction(normA/theta(3))
    expo = exponent(normA/theta(3))
    where (frac == 0.5_jprb)
      expo = expo - 1
    end where

    where (expo < 0)
      expo = 0
    end where

    ! Scale the input matrices by a power of 2
    scaling = 2.0_jprb**(-expo)
    do j3 = 1,m
      do j2 = 1,m
        A(1:iend,j2,j3) = A(1:iend,j2,j3) * scaling
      end do
    end do
    ! Pade approximant of degree 7
    A2 = mat_x_mat(n,iend,m,A, A, i_matrix_pattern)
    A4 = mat_x_mat(n,iend,m,A2,A2,i_matrix_pattern)
    A6 = mat_x_mat(n,iend,m,A2,A4,i_matrix_pattern)

    V = c(8)*A6 + c(6)*A4 + c(4)*A2
    do j3 = 1,m
      V(:,j3,j3) = V(:,j3,j3) + c(2)
    end do
    U = mat_x_mat(n,iend,m,A,V,i_matrix_pattern)
    V = c(7)*A6 + c(5)*A4 + c(3)*A2
    ! Add a multiple of the identity matrix
    do j3 = 1,m
      V(:,j3,j3) = V(:,j3,j3) + c(1)
    end do

    V = V-U
    U = 2.0_jprb*U
    A(1:iend,1:m,1:m) = solve_mat(n,iend,m,V,U)

    ! Add the identity matrix
    do j3 = 1,m
      A(1:iend,j3,j3) = A(1:iend,j3,j3) + 1.0_jprb
    end do

    ! Loop through the matrices
    do j1 = 1,iend
      if (expo(j1) > 0) then
        ! Square matrix j1 expo(j1) times          
        A(j1,:,:) = repeated_square(m,A(j1,:,:),expo(j1),i_matrix_pattern)
      end if
    end do

    if (lhook) call dr_hook('radiation_matrix:expm',1,hook_handle)

  end subroutine expm

  subroutine expm_lw(n,A)

    use yomhook, only : lhook, dr_hook, jphook

    integer,    intent(in)      :: n
    real(jprb), intent(inout)   :: A(n,6,6)

    real(jprb), parameter :: theta(3) = (/4.258730016922831e-01_jprb, &
         &                                1.880152677804762e+00_jprb, &
         &                                3.925724783138660e+00_jprb/) 
    real(jprb), parameter :: c(8) = (/17297280.0_jprb, 8648640.0_jprb, &
         &                1995840.0_jprb, 277200.0_jprb, 25200.0_jprb, &
         &                1512.0_jprb, 56.0_jprb, 1.0_jprb/)

    real(jprb), dimension(n,6,6) :: A2, A4, A6
    real(jprb), dimension(n,6,6) :: U, V
    real(jprb), dimension(6,6)   :: temp_in, temp_out

    real(jprb) :: normA(n), sum_column(n)
    ! real(jprb) :: frac(n)
    ! real(jprb) :: scaling(n)
    integer    :: j1, j2, j3, minexpo, nrepeat, jn
    integer    :: expo(n)

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_matrix:expm',0,hook_handle)

    normA = 0.0_jprb

    ! Compute the 1-norms of A
    do j3 = 1,6
      sum_column(:) = 0.0_jprb
      do j2 = 1,6
        do j1 = 1,n
          sum_column(j1) = sum_column(j1) + abs(A(j1,j2,j3))
        end do
      end do
      do j1 = 1,n
        if (sum_column(j1) > normA(j1)) then
          normA(j1) = sum_column(j1)
        end if
      end do
    end do


    associate(frac=>normA, scaling=>normA, normdiv=>sum_column)

      normdiv = normA/theta(3)
      frac = fraction(normdiv)
      expo = exponent(normdiv)
      where (frac == 0.5_jprb)
        expo = expo - 1
      end where
      expo = max(expo,0)
      minexpo = minval(expo)

      ! Scale the input matrices by a power of 2
      scaling = 2.0_jprb**(-expo)
      do j3 = 1,6
        do j2 = 1,6
          A(:,j2,j3) = A(:,j2,j3) * scaling
        end do
      end do

    end associate

    ! Pade approximant of degree 7
    A2 = mat_x_mat_6(n,A, A)
    A4 = mat_x_mat_6(n,A2,A2)
    A6 = mat_x_mat_6(n,A2,A4)

    V = c(8)*A6 + c(6)*A4 + c(4)*A2
    do j3 = 1,6
      V(:,j3,j3) = V(:,j3,j3) + c(2)
    end do
    U = mat_x_mat_6(n,A,V)
    V = c(7)*A6 + c(5)*A4 + c(3)*A2

    ! Add a multiple of the identity matrix
    do j3 = 1,6
      V(:,j3,j3) = V(:,j3,j3) + c(1)
    end do

    V = V-U
    A = 2.0_jprb*U

    call solve_mat_inplace(n,6,V,A)

    ! Add the identity matrix
    do j3 = 1,6
      A(:,j3,j3) = A(:,j3,j3) + 1.0_jprb
    end do

    ! To improve efficiency, square all matrices with the minimum expo first, and then square individual matrices as needed
    do j1 = 1,minexpo
      V = mat_x_mat_6(n,A,A) 
      A = V
    end do
    
    expo = expo - minexpo

    do j1 = 1,n
      if (expo(j1) > 0) then
        nrepeat = expo(j1)
        ! print *, "jg ", j1, "/", n, " expo", expo(j1)
        !Square matrix nrepeat times 
        temp_in =  A(j1,:,:)

        ! call repeated_square_6(nrepeat,temp_in,temp_out)
        do jn = 1,nrepeat
          do j2 = 1,6
            do j3 = 1,6
              temp_out(j3,j2) =  temp_in(j3,1)*temp_in(1,j2) + temp_in(j3,2)*temp_in(2,j2) + &
                   &      temp_in(j3,3)*temp_in(3,j2) + temp_in(j3,4)*temp_in(4,j2) + &
                   &      temp_in(j3,5)*temp_in(5,j2) + temp_in(j3,6)*temp_in(6,j2) 
            end do
          end do
          if (jn < nrepeat) then
            temp_in = temp_out
          end if
        end do
        
        A(j1,:,:) = temp_out
      end if
    end do

    if (lhook) call dr_hook('radiation_matrix:expm',1,hook_handle)

  end subroutine expm_lw

  !---------------------------------------------------------------------
  ! Like expm, but optimized for the shortwave, which has
  ! a special matrix structure with zeros and repeated elements. 
  ! Further assumes nreg = 3  =>  m = 9, and computations are done for all 
  ! N elements, where N (the inner dimension) is the number of individual matrices
  ! For performance reasons, N is increased by batching several cloudy layers in the
  ! higher level code: N = nlev_b*ng_sw, where ng_sw is the number of g-points and nlev_b
  ! is the number of batched layers (this is capped at 6 for ecCKD, i.e. nlev_b=1..6)
  subroutine expm_sw(N,ng_sw_in,nlev_b,A)

    use yomhook, only : lhook, dr_hook, jphook

    integer,    intent(in)      :: N, ng_sw_in, nlev_b
    real(jprb), intent(inout)   :: A(ng_sw*nlev_b,9,9)
    !dir$ assume_aligned A:64
    real(jprb), parameter :: theta(3) = (/4.258730016922831e-01_jprb, &
         &                                1.880152677804762e+00_jprb, &
         &                                3.925724783138660e+00_jprb/) 
    real(jprb), parameter :: c(8) = (/17297280.0_jprb, 8648640.0_jprb, &
         &                1995840.0_jprb, 277200.0_jprb, 25200.0_jprb, &
         &                1512.0_jprb, 56.0_jprb, 1.0_jprb/)
    real(jprb), dimension(ng_sw*nlev_b,9,9) :: A2, A4, A6
    real(jprb), dimension(ng_sw*nlev_b,9,9) :: V
    real(jprb), dimension(9,9)    :: temp_in, temp_out
    real(jprb), dimension(ng_sw*nlev_b) :: normA, sum_column
    integer    :: j1, j2, j3, jg, minexpo, nrepeat, jS, jE
    integer    :: expo(ng_sw*nlev_b)
    real(jphook) :: hook_handle
#ifdef USE_TIMING
    integer :: ret
#endif

    if (lhook) call dr_hook('radiation_matrix:expm_sw',0,hook_handle)
#ifdef USE_TIMING
    ret =  gptlstart('expm_sw')
#endif
    ! Compute the 1-norms of A
    normA = 0.0_jprb 
    do j3 = 1,9
      sum_column(:) = 0.0_jprb
      do j2 = 1,9
        do j1 = 1,N
          sum_column(j1) = sum_column(j1) + abs(A(j1,j2,j3))
        end do
      end do
      normA = max(normA,sum_column)
    end do

    associate(frac=>normA, scaling=>normA, normdiv=>sum_column, V2=>A6, U=>A2)

      normdiv = normA/theta(3)
      frac = fraction(normdiv)
      expo = exponent(normdiv)
      where (frac == 0.5_jprb)
        expo = expo - 1
      end where
      expo = max(expo,0)
      minexpo = minval(expo)
      
      ! Scale the input matrices by a power of 2
      scaling = 2.0_jprb**(-expo)
      do j3 = 1,9
        do j2 = 1,9
          if (.not.((j2>6) .and. (j3<7))) then
            A(:,j2,j3) = A(:,j2,j3) * scaling
          end if 
        end do
      end do

      ! Pade approximant of degree 7
      ! Input and output matrices have zeroes in the lower left corner AND repeated elements
      ! call mat_square_sw_repeats(N,A,A2) 
      ! call mat_square_sw_repeats(N,A2,A4)
      ! call mat_x_mat_sw_repeats(N,A2,A4,A6) 
      call mat_x_mat_sw_repeats(ng_sw,nlev_b,A,A,A2)
      call mat_x_mat_sw_repeats(ng_sw,nlev_b,A2,A2,A4)
      call mat_x_mat_sw_repeats(ng_sw,nlev_b,A2,A4,A6) 
      do j3 = 1,9
        do j2 = 1,9
          if (.not.((j2>6) .and. (j3<7))) then
            V(:,j2,j3)  = c(8)*A6(:,j2,j3) + c(6)*A4(:,j2,j3) + c(4)*A2(:,j2,j3)
            V2(:,j2,j3) = c(7)*A6(:,j2,j3) + c(5)*A4(:,j2,j3) + c(3)*A2(:,j2,j3)
          end if 
        end do
        ! Add a multiple of the identity matrix
        V(:,j3,j3)  = V(:,j3,j3)  + c(2)
        V2(:,j3,j3) = V2(:,j3,j3) + c(1)
      end do

      ! call mat_x_mat_sw(N,A,V,U) ! matrices have zeroes in the lower left corner, but no repeats
      call mat_x_mat_sw(ng_sw,nlev_b,A,V,U) ! matrices have zeroes in the lower left corner, but no repeats

      do j3 = 1,9
        do j2 = 1,9
          if (.not.((j2>6) .and. (j3<7))) then
            V2(:,j2,j3) = V2(:,j2,j3) - U(:,j2,j3) 
            A(:,j2,j3)  = 2.0_jprb*U(:,j2,j3)
          end if 
        end do 
      end do 
 
      call solve_mat_sw(ng_sw,nlev_b,V2,A)

    end associate 

    ! Add the identity matrix
    do j3 = 1,9
      A(:,j3,j3) = A(:,j3,j3) + 1.0_jprb
    end do

    ! Loop through the matrices    
    ! Improve efficiency by squaring all matrices with the minimum expo first
    do j1 = 1,minexpo
      call mat_square_sw(ng_sw,nlev_b,A,V) ! Matrices have zeroes in the corner but no repeated elements
      A = V
    end do

    expo = expo - minexpo

    ! Now square individual matrices as needed, depending on the expo of that matrix
    ! if (any(expo>0)) then 
    !   do jg = 1,ng
    !     if (expo(jg) > 0) then
    !       ! print *, "jg ", jg, "/", ng, " expo", expo(jg)
    !       nrepeat = expo(jg)
    !       !Square matrix nrepeat times 
    !       temp_in =  A(jg,:,:)
    !       call repeated_square_sw_9(nrepeat,temp_in,temp_out)
    !       A(jg,:,:) = temp_out
    !     end if
    !   end do
    ! end if

    ! Improve efficiency by squaring a group of matrices, based on array indexing
    ! consecutive matrices (in first dimension of A) that still need to be squared
    do while (any(expo>0))
      ! Find consecutive indices that should be squared
      jS=0; jE=0
      do jg = 1,N
        if (expo(jg)>0 .and. jS==0) jS = jg 
        if (jS/=0 .and. expo(jg)<=0) then
            jE = jg-1
            exit
        end if  
      end do
      if (jg==N+1) jE=N 

      if((jE-jS+1)<12) then ! not enough small matrices to make array operations worth it
        do jg = jS,jE
          nrepeat = expo(jg)
          !Square individual matrix nrepeat times 
          temp_in =  A(jg,:,:)
          call repeated_square_sw_9(nrepeat,temp_in,temp_out)
          A(jg,:,:) = temp_out
          expo(jg) = expo(jg)-nrepeat
        end do 
      else ! Square group of matrices
        do j2 = 1,6  
          do j1 = 1,6 
            V(jS:jE,j1,j2) = A(jS:jE,j1,1)*A(jS:jE,1,j2) + A(jS:jE,j1,2)*A(jS:jE,2,j2) &
                &          + A(jS:jE,j1,3)*A(jS:jE,3,j2) + A(jS:jE,j1,4)*A(jS:jE,4,j2) &
                &          + A(jS:jE,j1,5)*A(jS:jE,5,j2) + A(jS:jE,j1,6)*A(jS:jE,6,j2)
          end do
        end do
        do j2 = 7,9
          ! Do the top-right (E & H)
          do j1 =  1,6
            V(jS:jE,j1,j2) = A(jS:jE,j1,1)*A(jS:jE,1,j2) + A(jS:jE,j1,2)*A(jS:jE,2,j2) &
            &   + A(jS:jE,j1,3)*A(jS:jE,3,j2) + A(jS:jE,j1,4)*A(jS:jE,4,j2) + A(jS:jE,j1,5)*A(jS:jE,5,j2) &
            &   + A(jS:jE,j1,6)*A(jS:jE,6,j2) + A(jS:jE,j1,7)*A(jS:jE,7,j2) + A(jS:jE,j1,8)*A(jS:jE,8,j2) &
            &   + A(jS:jE,j1,9)*A(jS:jE,9,j2)
          end do
          ! Do the bottom-right (I)
          do j1 = 7,9
            V(jS:jE,j1,j2) =  A(jS:jE,j1,7)*A(jS:jE,7,j2) &
                &           + A(jS:jE,j1,8)*A(jS:jE,8,j2) + A(jS:jE,j1,9)*A(jS:jE,9,j2)
          end do
        end do
        A(jS:jE,:,:) = V(jS:jE,:,:)
        expo(jS:jE) = expo(jS:jE) - 1
      end if
    end do 
#ifdef USE_TIMING
    ret =  gptlstop('expm_sw')
#endif
    if (lhook) call dr_hook('radiation_matrix:expm_sw',1,hook_handle)

  end subroutine expm_sw


  !---------------------------------------------------------------------
  ! Return the matrix exponential of n 2x2 matrices representing
  ! conservative exchange between SPARTACUS regions, where the
  ! matrices have the structure
  !   (-a   b)
  !   ( a  -b)
  ! and a and b are assumed to be positive or zero.  The solution uses
  ! Putzer's algorithm - see the appendix of Hogan et al. (GMD 2018)
  subroutine fast_expm_exchange_2(n,iend,a,b,R)

    use yomhook, only : lhook, dr_hook, jphook

    integer,                      intent(in)  :: n, iend
    real(jprb), dimension(n),     intent(in)  :: a, b
    real(jprb), dimension(n,2,2), intent(out) :: R

    real(jprb), dimension(iend) :: factor

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_matrix:fast_expm_exchange_2',0,hook_handle)

    ! Security to ensure that if a==b==0 then the identity matrix is returned
    factor = (1.0_jprb - exp(-(a(1:iend)+b(1:iend))))/max(1.0e-12_jprb,a(1:iend)+b(1:iend))

    R(1:iend,1,1) = 1.0_jprb - factor*a(1:iend)
    R(1:iend,2,1) = factor*a(1:iend)
    R(1:iend,1,2) = factor*b(1:iend)
    R(1:iend,2,2) = 1.0_jprb - factor*b(1:iend)

    if (lhook) call dr_hook('radiation_matrix:fast_expm_exchange_2',1,hook_handle)

  end subroutine fast_expm_exchange_2

  !---------------------------------------------------------------------
  ! Return the matrix exponential of n 3x3 matrices representing
  ! conservative exchange between SPARTACUS regions, where the
  ! matrices have the structure
  !   (-a   b   0)
  !   ( a -b-c  d)
  !   ( 0   c  -d)
  ! and a-d are assumed to be positive or zero.  The solution uses the
  ! diagonalization method and is a slight generalization of the
  ! solution provided in the appendix of Hogan et al. (GMD 2018),
  ! which assumed c==d.
  ! subroutine fast_expm_exchange_3(n,iend,a,b,c,d,R)

  !   use yomhook, only : lhook, dr_hook, jphook

  !   real(jprb), parameter :: my_epsilon = 1.0e-12_jprb

  !   integer,                      intent(in)  :: n, iend
  !   real(jprb), dimension(n),     intent(in)  :: a, b, c, d
  !   real(jprb), dimension(n,3,3), intent(out) :: R

  !   ! Eigenvectors
  !   real(jprb), dimension(iend,3,3) :: V

  !   ! Non-zero Eigenvalues
  !   real(jprb), dimension(iend) :: lambda1, lambda2

  !   ! Diagonal matrix of the exponential of the eigenvalues
  !   real(jprb), dimension(iend,3) :: diag

  !   ! Result of diag right-divided by V
  !   real(jprb), dimension(iend,3,3) :: diag_rdivide_V

  !   ! Intermediate arrays
  !   real(jprb), dimension(iend) :: tmp1, tmp2

  !   integer :: j1, j2

  !   real(jphook) :: hook_handle

  !   if (lhook) call dr_hook('radiation_matrix:fast_expm_exchange_3',0,hook_handle)

  !   ! Eigenvalues
  !   tmp1 = 0.5_jprb * (a(1:iend)+b(1:iend)+c(1:iend)+d(1:iend))
  !   tmp2 = sqrt(tmp1*tmp1 - (a(1:iend)*c(1:iend) + a(1:iend)*d(1:iend) + b(1:iend)*d(1:iend)))
  !   lambda1 = -tmp1 + tmp2
  !   lambda2 = -tmp1 - tmp2

  !   ! Eigenvectors, with securities such taht if a--d are all zero
  !   ! then V is non-singular and the identity matrix is returned in R;
  !   ! note that lambdaX is typically negative so we need a
  !   ! sign-preserving security
  !   V(1:iend,1,1) = max(my_epsilon, b(1:iend)) &
  !        &  / sign(max(my_epsilon, abs(a(1:iend) + lambda1)), a(1:iend) + lambda1)
  !   V(1:iend,1,2) = b(1:iend) &
  !        &  / sign(max(my_epsilon, abs(a(1:iend) + lambda2)), a(1:iend) + lambda2)
  !   V(1:iend,1,3) = b(1:iend) / max(my_epsilon, a(1:iend))
  !   V(1:iend,2,:) = 1.0_jprb
  !   V(1:iend,3,1) = c(1:iend) &
  !        &  / sign(max(my_epsilon, abs(d(1:iend) + lambda1)), d(1:iend) + lambda1)
  !   V(1:iend,3,2) = c(1:iend) &
  !        &  / sign(max(my_epsilon, abs(d(1:iend) + lambda2)), d(1:iend) + lambda2)
  !   V(1:iend,3,3) = max(my_epsilon, c(1:iend)) / max(my_epsilon, d(1:iend))
    
  !   diag(:,1) = exp(lambda1)
  !   diag(:,2) = exp(lambda2)
  !   diag(:,3) = 1.0_jprb

  !   ! Compute diag_rdivide_V = diag * V^-1
  !   call diag_mat_right_divide_3(iend,iend,V,diag,diag_rdivide_V)

  !   ! Compute V * diag_rdivide_V
  !   do j1 = 1,3
  !     do j2 = 1,3
  !       R(1:iend,j2,j1) = V(1:iend,j2,1)*diag_rdivide_V(1:iend,1,j1) &
  !            &          + V(1:iend,j2,2)*diag_rdivide_V(1:iend,2,j1) &
  !            &          + V(1:iend,j2,3)*diag_rdivide_V(1:iend,3,j1)
  !     end do
  !   end do

  !   if (lhook) call dr_hook('radiation_matrix:fast_expm_exchange_3',1,hook_handle)

  ! end subroutine fast_expm_exchange_3

  subroutine fast_expm_exchange_3(ng_sw_in,R)

    use yomhook, only : lhook, dr_hook, jphook

! #ifdef PARKIND1_SINGLE
!     real(jprb), parameter :: my_epsilon = 10.0_jprb * epsilon(1.0_jprb)
! #else
    real(jprb), parameter :: my_epsilon = 1.0e-12_jprb
! #endif

    integer,                     intent(in)  :: ng_sw_in
    ! real(jprb), dimension(ng_sww),intent(in)  :: a, b, c, d
    real(jprb), dimension(ng_sww,3,3), intent(inout) :: R

    ! Eigenvectors
    real(jprb), dimension(ng_sww,3,3) :: V

    ! Non-zero Eigenvalues
    real(jprb), dimension(ng_sww,2) :: lambda

    ! Diagonal matrix of the exponential of the eigenvalues
    real(jprb), dimension(ng_sww,2) :: diag

    ! Result of diag right-divided by V
    real(jprb), dimension(ng_sww,3,3) :: X

    ! Intermediate arrays
    real(jprb), dimension(ng_sww) :: tmp1, tmp2
    real(jprb) :: L21, L31, L32, U22_inv, U33_inv, U23, y2, V11_inv
    integer :: j1, j2, jg
    real(jphook) :: hook_handle
! #define DEBUG_FAST_EXPM 1
#ifdef DEBUG_FAST_EXPM
    ! TEMPORARY
    real(jprb), dimension(ng_sww) :: denom
#endif
#ifdef __INTEL_COMPILER
    integer :: ones(ng_sww), minusones(ng_sww), signs(ng_sww)
    ones = 1 
    minusones = -1
#endif


    if (lhook) call dr_hook('radiation_matrix:fast_expm_exchange_3',0,hook_handle)
    associate(a=>R(:,2,1),b=>R(:,1,2),c=>R(:,3,2),d=>R(:,2,3))
      ! Eigenvalues
      tmp1 = 0.5_jprb * (a(:)+b(:)+c(:)+d(:))
      tmp2 = sqrt(max(0.0_jprb, tmp1*tmp1 - (a(:)*c(:) + a(:)*d(:) + b(:)*d(:))))
      ! The eigenvalues must not be the same or the LU decomposition
      ! fails; this can occur occasionally in single precision, which we
      ! avoid by limiting the minimum value of tmp2
      tmp2 = max(tmp2, epsilon(1.0_jprb) * tmp1)

      lambda(:,1) = -tmp1 + tmp2
      lambda(:,2) = -tmp1 - tmp2

      ! Eigenvectors, with securities such taht if a--d are all zero
      ! then V is non-singular and the identity matrix is returned in R;
      ! note that lambdaX is typically negative so we need a
      ! sign-preserving security
      tmp1 = a(:) + lambda(:,1) 
#ifdef __INTEL_COMPILER
      ! faster on ifort
      signs = merge(ones, minusones, tmp1>=0.0_jprb)        
      V(:,1,1) = max(my_epsilon, b(:)) / (signs*max(my_epsilon, abs(tmp1)))
#else 
      V(:,1,1) = max(my_epsilon, b(:)) / sign(max(my_epsilon, abs(tmp1)), tmp1)
#endif 
#ifdef DEBUG_FAST_EXPM
      denom =  1.0_jprb / V(:,1,1)
      if (any(denom==0.0_jprb)) then 
        R = 101.0_jprb
        print *, "1 divide by zero in fastexpm, returning"
        return 
      end if 
#endif
      tmp1 = a(:) + lambda(:,2)
#ifdef __INTEL_COMPILER
      signs = merge(ones, minusones, tmp1>=0.0_jprb)      
      V(:,1,2) = b(:) / (signs*max(my_epsilon, abs(tmp1))) 
#else
      V(:,1,2) = b(:) / sign(max(my_epsilon, abs(tmp1)), tmp1)
#endif 
      V(:,1,3) = b(:) / max(my_epsilon, a(:))
      V(:,2,:) = 1.0_jprb
      tmp1 = d(:) + lambda(:,1)
#ifdef __INTEL_COMPILER
      signs = merge(ones, minusones, tmp1>=0.0_jprb)
      V(:,3,1) = c(:) / (signs*max(my_epsilon, abs(tmp1)))
#else 
      V(:,3,1) = c(:) / sign(max(my_epsilon, abs(tmp1)), tmp1)
#endif 
      tmp1 = d(:) + lambda(:,2)
#ifdef __INTEL_COMPILER
      signs = merge(ones, minusones, tmp1>=0.0_jprb)
      V(:,3,2) = c(:) / (signs*max(my_epsilon, abs(tmp1)))
#else 
      V(:,3,2) = c(:) / sign(max(my_epsilon, abs(tmp1)), tmp1)
#endif 
      V(:,3,3) = max(my_epsilon, c(:)) / max(my_epsilon, d(:))
    end associate
#ifdef DEBUG_FAST_EXPM
    if (maxval(lambda)>20000.0_jprb) then 
      print *, "argument of exp was large, may cause FPE. max min:", maxval(lambda), minval(lambda)
      R = 1001.0_jprb
      print *, "returning from radiation_matrix and saving inputs"
      return
    end if
#endif

    diag = exp(lambda)
    ! diag(:,1) = exp(lambda(:,1))
    ! diag(:,2) = exp(lambda(:,2))
    ! diag(:,3) = 1.0_jprb

    ! ------ Compute X = diag * V^-1 ---------
    !  call diag_mat_right_divide_3(n,V,diag,X)

    ! LU decomposition of the *transpose* of V:
    !       ( 1        )   (U11 U12 U13)
    ! V^T = (L21  1    ) * (    U22 U23)
    !       (L31 L32  1)   (        U33)
    do jg = 1, ng_sww
      V11_inv = 1.0_jprb / V(jg,1,1)

      L21 = V(jg,1,2) *V11_inv ! / V(jg,1,1)
      L31 = V(jg,1,3) *V11_inv ! / V(jg,1,1)
      !U22 = V(jg,2,2) - L21*V(jg,2,1)
#ifdef DEBUG_FAST_EXPM
      denom(jg) = (V(jg,2,2) - L21*V(jg,2,1))
      if (denom(jg)==0.0_jprb) then 
        R = 101.0_jprb
        print *, "2 divide by zero in fastexpm, returning"
        return 
      end if  
#endif
      U22_inv = 1.0_jprb / (V(jg,2,2) - L21*V(jg,2,1))
      U23 = V(jg,3,2) - L21*V(jg,3,1)
      L32 =(V(jg,2,3) - L31*V(jg,2,1)) *U22_inv ! / U22
      ! U33 = V(jg,3,3) - L31*V(jg,3,1) - L32*U23
#ifdef DEBUG_FAST_EXPM
      denom(jg) = ( V(jg,3,3) - L31*V(jg,3,1) - L32*U23 )
      if (denom(jg)==0.0_jprb) then 
        R = 101.0_jprb
        print *, "3 divide by zero in fastexpm, returning"
        return 
      end if 
#endif
      U33_inv = 1.0_jprb / ( V(jg,3,3) - L31*V(jg,3,1) - L32*U23 )

      ! Solve X(1,:) = V^-T ( diag(1) )
      !                     (  0   )
      !                     (  0   )
      ! Solve Ly = diag(:,:,j) by forward substitution
      ! y1 = diag(jg,1)
      y2 = - L21*diag(jg,1)
      ! y3 = - L31*diag(jg,1) - L32*y2
      ! Solve UX = y by back substitution
      X(jg,1,3) = (-L31*diag(jg,1) - L32*y2) *U33_inv ! y3 / U33
      X(jg,1,2) = (y2 - U23*X(jg,1,3)) *U22_inv ! / U22
      X(jg,1,1) = (diag(jg,1) - V(jg,2,1)*X(jg,1,2) &
          &          - V(jg,3,1)*X(jg,1,3)) *V11_inv ! / V(jg,1,1)

      ! Solve X(2,jg) = V^-T (  0   )
      !                     ( diag(2) )
      !                     (  0   )
      ! Solve Ly = diag(:,:,j) by forward substitution
      ! y1 = 0
      ! y2 = diag(jg,2)
      ! y3 = - L32*diag(jg,2)
      ! Solve UX = y by back substitution
      X(jg,2,3) = (-L32*diag(jg,2)) *U33_inv ! y3 / U33
      X(jg,2,2) = (diag(jg,2) - U23*X(jg,2,3)) *U22_inv ! / U22
      X(jg,2,1) = (-V(jg,2,1)*X(jg,2,2) &
          &           -V(jg,3,1)*X(jg,2,3)) * V11_inv ! / V(jg,1,1)

      ! Solve X(3,jg) = V^-T (  0   )
      !                     (  0   )
      !                     ( diag(3) )
      ! Solve Ly = diag(:,:,j) by forward substitution
      ! y1 = 0
      ! y2 = 0
      ! y3 = diag(jg,3)
      ! Solve UX = y by back substitution
      X(jg,3,3) = U33_inv ! diag(jg,3) / U33
      X(jg,3,2) = -U23*X(jg,3,3) *U22_inv ! / U22
      X(jg,3,1) = (-V(jg,2,1)*X(jg,3,2) &
          &          - V(jg,3,1)*X(jg,3,3)) *V11_inv ! / V(jg,1,1)
    end do 

    ! Compute V * X
    do j1 = 1,3
      do j2 = 1,3
        R(:,j2,j1) = V(:,j2,1)*X(:,1,j1) &
            &     + V(:,j2,2)*X(:,2,j1) &
            &     + V(:,j2,3)*X(:,3,j1)
      end do
    end do

  if (lhook) call dr_hook('radiation_matrix:fast_expm_exchange_3',1,hook_handle)

  end subroutine fast_expm_exchange_3
!  generic :: fast_expm_exchange => fast_expm_exchange_2, fast_expm_exchange_3


end module radiation_matrix
