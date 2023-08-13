MODULE utilities
  IMPLICIT NONE

  INTEGER, PARAMETER, PRIVATE :: qsort_limit = 10

  INTEGER, PARAMETER, PRIVATE :: rnd_n = 624
  INTEGER, PARAMETER, PRIVATE :: rnd_m = 397
  INTEGER, PARAMETER, PRIVATE :: rnd_mata = -1727483681
  INTEGER, PARAMETER, PRIVATE :: rnd_umask = -2147483648
  INTEGER, PARAMETER, PRIVATE :: rnd_lmask = 2147483647
  INTEGER, PARAMETER, PRIVATE :: rnd_tmaskb = -1658038656
  INTEGER, PARAMETER, PRIVATE :: rnd_tmaskc = -272236544

  CHARACTER(LEN=3), PARAMETER :: af = "(A)"

  INTEGER, PRIVATE :: rnd_mt(0:rnd_n-1)
  INTEGER, PRIVATE :: rnd_mag01(0:1)
  DATA rnd_mag01 /0, rnd_mata/
  INTEGER, PRIVATE :: rnd_mti = rnd_n + 1
  INTEGER, PRIVATE :: available_handle = 10

CONTAINS
  SUBROUTINE qsort_real(vector, permutation, missing)
    REAL :: vector(:)
    INTEGER, INTENT(OUT) :: permutation(:)
    INTEGER, OPTIONAL, INTENT(OUT) :: missing
    INTEGER :: i
    LOGICAL :: ismissing ! C function

    IF (PRESENT(missing)) missing = 0
    DO i = 1, UBOUND(vector, 1)
      IF (ismissing(vector(i))) THEN
        permutation(i) = 0
        IF (PRESENT(missing)) missing = missing + 1
      ELSE
        permutation(i) = i
      END IF
    END DO
    CALL qsort_real_inner(vector, permutation, 1, UBOUND(vector, 1))
  END SUBROUTINE qsort_real

  RECURSIVE SUBROUTINE qsort_real_inner(vector, permutation, ileft, iright)
    REAL :: vector(:)
    INTEGER, INTENT(OUT) :: permutation(:)
    INTEGER :: left, right, ileft, iright
    INTEGER :: l, r, ll, rr, t, i
    REAL :: pivot, el, maxreal
    left = ileft
    right = iright
    maxreal = HUGE(0.0)

    DO WHILE (right - left .GE. qsort_limit)
      ! skip all missing values at the end
      DO WHILE (right .GT. left .AND. permutation(right) .EQ. 0)
        right = right - 1
      END DO

      ! vector boundaries: left | e | l | r | z | right
      r = right
      rr = right
      ll = left
      l = left

      ! choose pivot (the last element)

      IF (permutation(right).EQ.0) THEN
        pivot = maxreal
      ELSE
        pivot = vector(permutation(right))
      END IF

      ! partition; push equals to the side
      DO
        ! search left
        DO WHILE (l .LE. r)
          IF (permutation(l).EQ.0) THEN
            el = maxreal
          ELSE
            el = vector(permutation(l))
          END IF
          IF (el .GT. pivot) THEN
            EXIT
          ELSE IF (el .EQ. pivot) THEN
            t = permutation(l)
            permutation(l) = permutation(ll)
            permutation(ll) = t
            ll = ll + 1
          END IF
          l = l + 1
        END DO

        ! search right
        DO WHILE (l .LE. r)
          IF (permutation(r).EQ.0) THEN
            el = maxreal
          ELSE
            el = vector(permutation(r))
          END IF
          IF (el .LT. pivot) THEN
            EXIT
          ELSE IF (el .EQ. pivot) THEN
            t = permutation(r)
            permutation(r) = permutation(rr)
            permutation(rr) = t
            rr = rr - 1
          END IF
          r = r - 1
        END DO

        IF (l .GT. r) EXIT

        t = permutation(l)
        permutation(l) = permutation(r)
        permutation(r) = t
        l = l + 1
        r = r - 1
      END DO

      DO i = 0, MIN(ll - left, l - ll) - 1
        t = permutation(left + i)
        permutation(left + i) = permutation(l - i - 1)
        permutation(l - i - 1) = t
      END DO
      DO i = 0, MIN(right - rr, rr - r) - 1
        t = permutation(right - i)
        permutation(right - i) = permutation(r + i + 1)
        permutation(r + i + 1) = t
      END DO

      ! recurse on the smaller partition, then iterate on the larger
      ! to avoid getting stack overflow on degenerate cases
      IF (rr - r .GE. l - ll) THEN
        IF (l - ll .GT. 1) THEN
          CALL qsort_real_inner(vector, permutation, left, left + l - ll - 1)
        END IF
        left = r + 1
        right = right
      ELSE
        IF (rr - r .GT. 1) THEN
          CALL qsort_real_inner(vector, permutation, r + 1, right)
        END IF
        right = left + l - ll - 1
      END IF
    END DO

    ! sort the rest using insertion sort
    DO r = left + 1, right ! take a value
      IF (permutation(r) .NE. 0) THEN
        DO l = left, r - 1 ! find a place to put it
          IF (permutation(l) .EQ. 0) EXIT
          IF (vector(permutation(l)) .GT. vector(permutation(r))) EXIT
        END DO
        t = permutation(r)
        permutation(l+1:r) = permutation(l:r-1)
        permutation(l) = t
      END IF
    END DO
  END SUBROUTINE qsort_real_inner

  SUBROUTINE shuffle(source, dest)
    INTEGER :: source(:)
    INTEGER, INTENT(OUT) :: dest(:)
    INTEGER :: i, j, tmp, bound

    dest = source
    bound = UBOUND(source, 1)
    DO i = bound, 2, -1
      j = rnd_int(i - 1)
      tmp = dest(j)
      dest(j) = dest(i)
      dest(i) = tmp
    END DO
  END SUBROUTINE shuffle

  FUNCTION median(vector)
    REAL :: vector(:)
    REAL :: median
    INTEGER :: permutation(UBOUND(vector, 1))
    INTEGER :: missing
    CALL qsort_real(vector, permutation, missing)
    median = vector(permutation((UBOUND(vector, 1)+1-missing)/2))
  END FUNCTION median

  FUNCTION most_frequent(vector, cat_count, freqs)
    INTEGER :: vector(:)
    INTEGER :: cat_count
    INTEGER, POINTER, OPTIONAL :: freqs(:)
    INTEGER :: most_frequent
    INTEGER, POINTER :: cat_frequencies(:)
    INTEGER :: i
    INTEGER :: most_freq_loc(1)

    IF (PRESENT(freqs)) THEN
      cat_frequencies => freqs
    ELSE
      ALLOCATE (cat_frequencies(cat_count))
    END IF
    cat_frequencies = 0
    DO i = 1, UBOUND(vector, 1)
      IF (vector(i).NE.0) THEN
        cat_frequencies(vector(i)) = cat_frequencies(vector(i)) + 1
      END IF
    END DO
    most_freq_loc = MAXLOC(cat_frequencies)
    most_frequent = most_freq_loc(1)
    
    IF (.NOT.PRESENT(freqs)) THEN
      DEALLOCATE (cat_frequencies)
    END IF
  END FUNCTION most_frequent

  FUNCTION lcase(string)
    CHARACTER(*), INTENT(IN) :: string
    INTEGER :: ls, lc
    CHARACTER(LEN=LEN(string)) :: lcase
    ls = LEN(string)
    lcase = string
    DO lc = 1, ls
      IF (LGE(lcase(lc:lc),'A').AND.LLE(lcase(lc:lc),'Z')) THEN
        lcase(lc:lc) = CHAR(ICHAR(lcase(lc:lc)) + 32)
      END IF
    END DO
  END FUNCTION lcase

  SUBROUTINE squeeze_num(str, x)
    CHARACTER(LEN=*) :: str
    REAL :: x
    INTEGER :: pos, pose, strlen
    IF (x.EQ.0) THEN
      str = "0"
      RETURN
    END IF
    WRITE(str, *) x
    strlen = LEN(str)
    str = ADJUSTR(str)
    pose = SCAN(str, 'E', .TRUE.)
    IF (pose.EQ.0) pose = strlen + 1
    pos = SCAN(str(1:pose-1),'123456789.',.TRUE.)
    IF (str(pos:pos).EQ.'.') pos = pos - 1
    IF (pose.GT.strlen) THEN
      str = ADJUSTL(str(1:pos))
    ELSE
      str = ADJUSTL(str(1:pos) // str(pose:))
    END IF
  END SUBROUTINE squeeze_num

  FUNCTION rnd_int(limit)
    INTEGER :: rnd_int
    INTEGER :: limit
    rnd_int = FLOOR(rnd() * limit + 1)
    IF (rnd_int.GT.limit) rnd_int = limit
  END FUNCTION rnd_int

  FUNCTION rnd_norm()
    REAL :: rnd_norm
    rnd_norm = SQRT(-2 * LOG(rnd())) * COS(6.283185 * rnd())
  END FUNCTION rnd_norm

  FUNCTION rnd()
    REAL :: rnd
    INTEGER :: y

    y = rnd_integer()
    IF (y.LT.0) THEN
      rnd = (DBLE(y) + 2.0d0**32) / (2.0d0**32 - 1.0d0)
    ELSE
      rnd = DBLE(y) / (2.0d0**32 - 1.0d0)
    END IF
  END FUNCTION rnd

  FUNCTION rnd_integer() RESULT (y)
    INTEGER :: kk, y

    IF (rnd_mti.GE.rnd_n) THEN
      IF (rnd_mti.EQ.(rnd_n + 1)) THEN
        CALL seed_rnd(4357) ! TODO seed with a random value
      END IF
      DO kk = 0, rnd_n - rnd_m - 1
        y = IOR(IAND(rnd_mt(kk), rnd_umask), IAND(rnd_mt(kk + 1), rnd_lmask))
        rnd_mt(kk)=IEOR(IEOR(rnd_mt(kk + rnd_m), ISHFT(y, -1)), &
          & rnd_mag01(IAND(y, 1)))
      END DO
      DO kk = rnd_n - rnd_m, rnd_n - 2
        y = IOR(IAND(rnd_mt(kk), rnd_umask), IAND(rnd_mt(kk + 1), rnd_lmask))
        rnd_mt(kk)=IEOR(IEOR(rnd_mt(kk + rnd_m - rnd_n), ISHFT(y, -1)), &
          & rnd_mag01(IAND(y, 1)))
      END DO
      y = IOR(IAND(rnd_mt(rnd_n - 1), rnd_umask), IAND(rnd_mt(0), rnd_lmask))
      rnd_mt(rnd_n - 1)=IEOR(IEOR(rnd_mt(rnd_m - 1), ISHFT(y, -1)), &
        & rnd_mag01(IAND(y, 1)))
      rnd_mti = 0
    END IF
    y = rnd_mt(rnd_mti)
    rnd_mti = rnd_mti + 1
    y = IEOR(y, ISHFT(y, -11))
    y = IEOR(y, IAND(ISHFT(y, 7), rnd_tmaskb))
    y = IEOR(y, IAND(ISHFT(y, 15), rnd_tmaskc))
    y = IEOR(y, ISHFT(y, -18))
  END FUNCTION rnd_integer

  SUBROUTINE seed_rnd(seed)
    INTEGER :: seed
    rnd_mt(0) = seed
    DO rnd_mti = 1, rnd_n - 1
      ! TODO is IAND() necessary?
      rnd_mt(rnd_mti) = 69069 * rnd_mt(rnd_mti - 1)
    END DO
  END SUBROUTINE seed_rnd

  FUNCTION open_file(filename, for_output, err) RESULT (handle)
    INTEGER :: handle
    CHARACTER(LEN=*) :: filename
    LOGICAL :: for_output
    LOGICAL, INTENT(OUT) :: err
    
    err = .FALSE.
    IF (TRIM(filename).NE."-") THEN
      handle = available_handle
      IF (for_output) THEN
        OPEN(handle, ACTION="WRITE", POSITION="REWIND", STATUS="REPLACE", &
          & ERR=9000, FILE=filename)
      ELSE
        OPEN(handle, ACTION="READ", POSITION="REWIND", STATUS="OLD", &
          & ERR=9000, FILE=filename)
      END IF
      available_handle = available_handle + 1
    ELSE
      IF (for_output) THEN
        handle = 6
      ELSE
        handle = 5
      END IF
    END IF
    RETURN

    9000 CONTINUE
    WRITE(0, af) "Error accessing file " // TRIM(filename)
    err = .TRUE.
  END FUNCTION open_file

  SUBROUTINE close_file(handle)
    INTEGER :: handle
    IF (handle.GE.10) CLOSE(handle, ERR=9000)
    9000 CONTINUE ! Ignore a close error
  END SUBROUTINE close_file

  PURE FUNCTION erfcc(x)
    REAL, INTENT(IN) :: x
    REAL :: t, z, erfcc
    z = ABS(x)/1.41421356
    t = 1./(1.+0.5*z)
    erfcc = t*EXP(-z*z-1.26551223+t*(1.00002368+t*(.37409196+t* &
      & (.09678418+t*(-.18628806+t*(.27886807+t*(-1.13520398+t* &
      & (1.48851587+t*(-.82215223+t*.17087277)))))))))
    erfcc = erfcc/2
    IF (x.LT.0.) erfcc = 2.-erfcc
  END FUNCTION erfcc
END MODULE utilities
