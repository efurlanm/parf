MODULE bitvectors
! Usage:
!
!   USE bitvectors
!   TYPE (bitvector), POINTER :: flags(:)
!
!   flags => new_bitvector(60)
!   CALL setbit(flags, 10, .TRUE.)
!   WRITE(*, *) getbit(flags, 10)
!   CALL free_bitvector(flags)

  USE utilities
  IMPLICIT NONE

  INTEGER, PARAMETER :: bits_per_integer = BIT_SIZE(1)

  TYPE bitvector
    INTEGER, POINTER :: bits(:)
    INTEGER :: width
  END TYPE bitvector
  TYPE bitvector_p
    TYPE (bitvector), POINTER :: p
  END TYPE bitvector_p

CONTAINS
  FUNCTION new_bitvector(vector_size) RESULT (bvptr)
    INTEGER :: vector_size
    TYPE (bitvector), POINTER :: bvptr
    INTEGER :: newsize
    newsize = (vector_size + bits_per_integer - 1) / bits_per_integer
    ALLOCATE (bvptr, bvptr%bits(newsize))
    bvptr%bits = 0
    bvptr%width = 1
  END FUNCTION new_bitvector

  FUNCTION new_bitvector2(size1, size2) RESULT (bvptr)
    INTEGER :: size1, size2
    TYPE (bitvector), POINTER :: bvptr
    INTEGER :: newsize
    newsize = (size1 * size2 + bits_per_integer - 1) / bits_per_integer
    ALLOCATE (bvptr, bvptr%bits(newsize))
    bvptr%bits = 0
    bvptr%width = size2
  END FUNCTION new_bitvector2

  FUNCTION clone_bitvector(oldbvptr) RESULT (bvptr)
    TYPE (bitvector), POINTER :: oldbvptr, bvptr
    ALLOCATE (bvptr, bvptr%bits(UBOUND(oldbvptr%bits, 1)))
    bvptr%width = oldbvptr%width
    bvptr%bits = oldbvptr%bits
  END FUNCTION clone_bitvector

  SUBROUTINE free_bitvector(bvptr)
    TYPE (bitvector), POINTER :: bvptr
    IF (ASSOCIATED(bvptr)) THEN
      IF (ASSOCIATED(bvptr%bits)) THEN
        DEALLOCATE (bvptr%bits)
      END IF
      DEALLOCATE (bvptr)
    END IF
  END SUBROUTINE free_bitvector

  SUBROUTINE setbit(bvptr, bit, value)
    TYPE (bitvector), POINTER :: bvptr
    INTEGER :: bit
    LOGICAL :: value
    INTEGER :: i
    INTEGER :: realbit
    realbit = bit - 1
    i = realbit / bits_per_integer + 1
    IF (value) THEN
      bvptr%bits(i) = IBSET(bvptr%bits(i), MODULO(realbit, bits_per_integer))
    ELSE
      bvptr%bits(i) = IBCLR(bvptr%bits(i), MODULO(realbit, bits_per_integer))
    END IF
  END SUBROUTINE setbit

  SUBROUTINE setbit2(bvptr, dim1, dim2, value)
    TYPE (bitvector), POINTER :: bvptr
    INTEGER :: dim1, dim2
    LOGICAL :: value
    CALL setbit(bvptr, (dim1 - 1) * bvptr%width + dim2, value)
  END SUBROUTINE setbit2

  FUNCTION getbit(bvptr, bit)
    TYPE (bitvector), POINTER :: bvptr
    INTEGER :: bit
    LOGICAL :: getbit
    INTEGER :: realbit
    realbit = bit - 1
    getbit = BTEST(bvptr%bits(realbit / bits_per_integer + 1), &
      & MODULO(realbit, bits_per_integer))
  END FUNCTION getbit

  FUNCTION getbit2(bvptr, dim1, dim2)
    TYPE (bitvector), POINTER :: bvptr
    INTEGER :: dim1, dim2
    LOGICAL :: getbit2
    getbit2 = getbit(bvptr, (dim1 - 1)* bvptr%width + dim2)
  END FUNCTION getbit2

  SUBROUTINE increment_bitvector(bvptr)
    TYPE (bitvector), POINTER :: bvptr
    INTEGER :: i
    DO i = 1, UBOUND(bvptr%bits, 1)
      bvptr%bits(i) = bvptr%bits(i) + 1
      IF (bvptr%bits(i).NE.0) RETURN
    END DO
    ! Should never reach here
  END SUBROUTINE

  SUBROUTINE randomise_bitvector(bvptr, bvptr_used)
    TYPE (bitvector), POINTER :: bvptr, bvptr_used
    INTEGER :: i

    DO i = 1, UBOUND(bvptr%bits, 1)
      bvptr%bits(i) = IAND(rnd_integer(), NOT(bvptr_used%bits(i)))
    END DO
  END SUBROUTINE randomise_bitvector

  FUNCTION is_and_bitvector(bvptr, bvptr2)
    TYPE (bitvector), POINTER :: bvptr, bvptr2
    INTEGER :: i
    LOGICAL :: is_and_bitvector

    DO i = 1, UBOUND(bvptr%bits, 1)
      IF (IAND(bvptr%bits(i), bvptr2%bits(i)).NE.0) THEN
        is_and_bitvector = .TRUE.
        RETURN
      END IF
    END DO
    is_and_bitvector = .FALSE.
  END FUNCTION is_and_bitvector

  SUBROUTINE or_bitvector(bvptr, bvptr2)
    TYPE (bitvector), POINTER :: bvptr, bvptr2
    INTEGER :: i

    DO i = 1, UBOUND(bvptr%bits, 1)
      bvptr%bits(i) = IOR(bvptr%bits(i), bvptr2%bits(i))
    END DO
  END SUBROUTINE or_bitvector

  SUBROUTINE eor_bitvector(bvptr, bvptr2)
    TYPE (bitvector), POINTER :: bvptr, bvptr2
    INTEGER :: i

    DO i = 1, UBOUND(bvptr%bits, 1)
      bvptr%bits(i) = IEOR(bvptr%bits(i), bvptr2%bits(i))
    END DO
  END SUBROUTINE eor_bitvector

  FUNCTION not_bitvector(bvptr, length)
    TYPE (bitvector), POINTER :: bvptr, not_bitvector
    INTEGER :: length
    INTEGER :: i, to_clear, top

    not_bitvector => new_bitvector(length)
    top = UBOUND(bvptr%bits, 1)
    DO i = 1, top
      not_bitvector%bits(i) = NOT(bvptr%bits(i))
    END DO
    to_clear = bits_per_integer - MODULO(length, bits_per_integer)
    not_bitvector%bits(top) = ISHFT(ISHFT(not_bitvector%bits(top), &
      & to_clear), -to_clear)
  END FUNCTION not_bitvector
END MODULE bitvectors
