MODULE parallel
  INTEGER :: par_rank
  INTEGER :: par_processes
  LOGICAL :: par_front

  INTEGER*8 :: slowness
CONTAINS
  SUBROUTINE par_init()
    par_rank = 0
    par_processes = 1
    par_front = .TRUE.
  END SUBROUTINE par_init

  SUBROUTINE par_finalize()
  END SUBROUTINE par_finalize

  SUBROUTINE par_bcast_char(message)
    CHARACTER(LEN=*) :: message
  END SUBROUTINE par_bcast_char

  SUBROUTINE par_bcast_int(scalar)
    INTEGER :: scalar
  END SUBROUTINE par_bcast_int

  SUBROUTINE par_bcast_int_vector(vector)
    INTEGER, POINTER :: vector(:)
  END SUBROUTINE par_bcast_int_vector

  SUBROUTINE par_bcast_bool_vector(vector)
    LOGICAL :: vector(:)
  END SUBROUTINE par_bcast_bool_vector

  SUBROUTINE par_bcast_float(scalar)
    REAL :: scalar
  END SUBROUTINE par_bcast_float

  SUBROUTINE par_bcast_bool(scalar)
    LOGICAL :: scalar
  END SUBROUTINE par_bcast_bool

  SUBROUTINE par_notify(process)
    INTEGER :: process
  END SUBROUTINE par_notify

  SUBROUTINE par_wait()
  END SUBROUTINE par_wait

  SUBROUTINE par_recv_int(scalar)
    INTEGER :: scalar
  END SUBROUTINE par_recv_int

  SUBROUTINE par_send_int(scalar, process)
    INTEGER :: scalar
    INTEGER :: process
  END SUBROUTINE par_send_int

  SUBROUTINE par_recv_int_vector(vector)
    INTEGER :: vector(:)
  END SUBROUTINE par_recv_int_vector

  SUBROUTINE par_send_int_vector(vector, process)
    INTEGER :: vector(:)
    INTEGER :: process
  END SUBROUTINE par_send_int_vector

  SUBROUTINE par_recv_real_vector(vector)
    REAL :: vector(:)
  END SUBROUTINE par_recv_real_vector

  SUBROUTINE par_send_real_vector(vector, process)
    REAL :: vector(:)
    INTEGER :: process
  END SUBROUTINE par_send_real_vector

  SUBROUTINE par_max_int(val, pos)
    INTEGER :: val, pos
  END SUBROUTINE par_max_int

  SUBROUTINE par_sum_real(scalar)
    REAL :: scalar
  END SUBROUTINE par_sum_real

  SUBROUTINE par_sum_real_vector(vector)
    REAL :: vector(:)
  END SUBROUTINE par_sum_real_vector

  SUBROUTINE par_sum_real_vector_at_one(source_vector, target_vector, proc)
    REAL :: source_vector(:), target_vector(:)
    INTEGER :: proc

    target_vector = source_vector
  END SUBROUTINE par_sum_real_vector_at_one

  SUBROUTINE par_sum_real_matrix(matrix)
    REAL :: matrix(:, :)
  END SUBROUTINE par_sum_real_matrix

  SUBROUTINE par_sum_int_vector(vector)
    INTEGER :: vector(:)
  END SUBROUTINE par_sum_int_vector

  SUBROUTINE par_gather_real(whole, part, stripes, offsets)
    REAL, POINTER :: whole(:), part(:)
    INTEGER, POINTER :: stripes(:), offsets(:)
    whole = part
  END SUBROUTINE par_gather_real

  SUBROUTINE par_gather_slowness()
    slowness_table = 1.0
  END SUBROUTINE par_gather_slowness

  SUBROUTINE par_get_stripe(n, lower, upper, stripes, offsets)
    INTEGER :: n
    INTEGER, INTENT(OUT) :: lower, upper
    INTEGER, POINTER, OPTIONAL :: stripes(:), offsets(:)

    lower = 1
    upper = n
    IF (PRESENT(offsets)) THEN
      ALLOCATE (stripes(0:0), offsets(0:0))
      stripes(0) = n
      offsets(0) = 0
    END IF
  END SUBROUTINE par_get_stripe
END MODULE parallel
