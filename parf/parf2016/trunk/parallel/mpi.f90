MODULE parallel
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  INTEGER :: par_rank
  INTEGER :: par_processes
  LOGICAL :: par_front

  REAL, POINTER :: slowness_table(:)
CONTAINS
  SUBROUTINE par_init()
    INTEGER :: mpi_err
    REAL :: slowness
    REAL :: benchmark ! C function

    CALL MPI_Init(mpi_err)
    IF (.TRUE.) THEN
      CALL MPI_Comm_rank(MPI_COMM_WORLD, par_rank, mpi_err)
      CALL par_handle_error(mpi_err)
      CALL MPI_Comm_size(MPI_COMM_WORLD, par_processes, mpi_err)
      CALL par_handle_error(mpi_err)
      par_front = par_rank.EQ.0
    END IF
    
    ALLOCATE (slowness_table(0:par_processes-1))
    IF (par_processes.GT.1) THEN
      slowness = benchmark()
      CALL MPI_Allgather(slowness, 1, MPI_REAL, &
        & slowness_table, 1, MPI_REAL, MPI_COMM_WORLD, mpi_err)
      CALL par_handle_error(mpi_err)
    ELSE
      slowness_table = 1.0
    END IF
  END SUBROUTINE par_init

  SUBROUTINE par_finalize()
    INTEGER :: mpi_err

    DEALLOCATE (slowness_table)

    CALL MPI_Finalize(mpi_err)
    CALL par_handle_error(mpi_err)
  END SUBROUTINE par_finalize

!
! The bcast, send and recv, which in MPI use MPI_datatype, are split
! in separate functions for each data type, as to keep compatibility
! with possible other communication protocols/interfaces
!

  SUBROUTINE par_bcast_char(message)
    CHARACTER(*) :: message
    INTEGER :: mpi_err
    CALL MPI_BCast(message, LEN(message), MPI_CHARACTER, 0, &
      & MPI_COMM_WORLD, mpi_err)    
    CALL par_handle_error(mpi_err)
  END SUBROUTINE par_bcast_char

  SUBROUTINE par_bcast_int(scalar)
    INTEGER :: scalar
    INTEGER :: mpi_err
    CALL MPI_BCast(scalar, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)    
    CALL par_handle_error(mpi_err)
  END SUBROUTINE par_bcast_int

  SUBROUTINE par_bcast_int_from_proc(scalar, proc)
    INTEGER :: scalar, proc
    INTEGER :: mpi_err
    CALL MPI_BCast(scalar, 1, MPI_INTEGER, proc, MPI_COMM_WORLD, mpi_err)    
    CALL par_handle_error(mpi_err)
  END SUBROUTINE par_bcast_int_from_proc

  SUBROUTINE par_bcast_int_vector(vector)
    INTEGER :: vector(:)
    INTEGER :: mpi_err
    CALL MPI_BCast(vector(1), UBOUND(vector), MPI_INTEGER, 0, &
      & MPI_COMM_WORLD, mpi_err)    
    CALL par_handle_error(mpi_err)
  END SUBROUTINE par_bcast_int_vector

  SUBROUTINE par_bcast_bool_vector(vector)
    LOGICAL :: vector(:)
    INTEGER :: mpi_err
    CALL MPI_BCast(vector(1), UBOUND(vector), MPI_LOGICAL, 0, &
      & MPI_COMM_WORLD, mpi_err)    
    CALL par_handle_error(mpi_err)
  END SUBROUTINE par_bcast_bool_vector

  SUBROUTINE par_bcast_float(scalar)
    REAL :: scalar
    INTEGER :: mpi_err
    CALL MPI_BCast(scalar, 1, MPI_REAL, 0, MPI_COMM_WORLD, mpi_err)    
    CALL par_handle_error(mpi_err)
  END SUBROUTINE par_bcast_float

  SUBROUTINE par_bcast_bool(boolean)
    LOGICAL :: boolean
    INTEGER :: mpi_err
    CALL MPI_BCast(boolean, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpi_err)    
    CALL par_handle_error(mpi_err)
  END SUBROUTINE par_bcast_bool

  SUBROUTINE par_send_int(scalar, processor)
    INTEGER :: processor
    INTEGER :: mpi_err
    INTEGER :: scalar

    CALL MPI_Send(scalar, 1, MPI_INTEGER, processor, processor, &
      & MPI_COMM_WORLD, mpi_err)
    CALL par_handle_error(mpi_err)
  END SUBROUTINE par_send_int

  SUBROUTINE par_recv_int(scalar)
    INTEGER :: mpi_err
    INTEGER :: scalar

    CALL MPI_Recv(scalar, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_ANY_TAG, &
      & MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
    CALL par_handle_error(mpi_err)
  END SUBROUTINE par_recv_int

  SUBROUTINE par_notify(processor)
    INTEGER :: processor

    CALL par_send_int(processor, processor)
  END SUBROUTINE par_notify

  SUBROUTINE par_wait()
    INTEGER :: dummy

    CALL par_recv_int(dummy)
  END SUBROUTINE par_wait

  SUBROUTINE par_send_int_vector(vector, processor)
    INTEGER :: vector(:)
    INTEGER :: processor
    INTEGER :: mpi_err

    CALL MPI_Send(vector, UBOUND(vector, 1), MPI_INTEGER, processor, &
      & processor, MPI_COMM_WORLD, mpi_err)
    CALL par_handle_error(mpi_err)
  END SUBROUTINE par_send_int_vector

  SUBROUTINE par_recv_int_vector(vector)
    INTEGER :: vector(:)
    INTEGER :: mpi_err

    CALL MPI_Recv(vector, UBOUND(vector, 1), MPI_INTEGER, MPI_ANY_SOURCE, &
      & MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
    CALL par_handle_error(mpi_err)
  END SUBROUTINE par_recv_int_vector

  SUBROUTINE par_send_real_vector(vector, processor)
    REAL :: vector(:)
    INTEGER :: processor
    INTEGER :: mpi_err

    CALL MPI_Send(vector, UBOUND(vector, 1), MPI_REAL, processor, par_rank, &
      & MPI_COMM_WORLD, mpi_err)
    CALL par_handle_error(mpi_err)
  END SUBROUTINE par_send_real_vector

  SUBROUTINE par_recv_real_vector(vector)
    REAL, POINTER :: vector(:)
    INTEGER :: mpi_err

    CALL MPI_Recv(vector, UBOUND(vector, 1), MPI_REAL, MPI_ANY_SOURCE, &
      & MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
    CALL par_handle_error(mpi_err)
  END SUBROUTINE par_recv_real_vector
!
! Reduction.
!
  SUBROUTINE par_max_int(val, pos)
    INTEGER :: val, pos, instr(2), outstr(2), mpi_err

    instr(1) = val
    instr(2) = pos
    CALL MPI_Allreduce(instr, outstr, 1, MPI_2INTEGER, &
      & MPI_MAXLOC, MPI_COMM_WORLD, mpi_err)
    CALL par_handle_error(mpi_err)
    val = outstr(1)
    pos = outstr(2)
  END SUBROUTINE par_max_int

  SUBROUTINE par_sum_real(scalar)
    INTEGER :: mpi_err
    REAL :: scalar
    REAL :: sum_result

    CALL MPI_Allreduce(scalar, sum_result, 1, MPI_REAL, &
      & MPI_SUM, MPI_COMM_WORLD, mpi_err)
    CALL par_handle_error(mpi_err)
    scalar = sum_result
  END SUBROUTINE par_sum_real

  SUBROUTINE par_sum_real_vector(vector)
    REAL :: vector(:)
    INTEGER :: mpi_err
    REAL, POINTER :: sum_result(:)

    ALLOCATE (sum_result(UBOUND(vector, 1)))
    CALL MPI_Allreduce(vector, sum_result, UBOUND(vector, 1), MPI_REAL, &
      & MPI_SUM, MPI_COMM_WORLD, mpi_err)
    CALL par_handle_error(mpi_err)
    vector = sum_result
    DEALLOCATE (sum_result)
  END SUBROUTINE par_sum_real_vector

  SUBROUTINE par_sum_real_vector_at_one(source_vector, target_vector, process)
    REAL :: source_vector(:), target_vector(:)
    INTEGER :: process
    INTEGER :: mpi_err

    CALL MPI_Reduce(source_vector, target_vector, UBOUND(source_vector, 1), &
      & MPI_REAL, MPI_SUM, process, MPI_COMM_WORLD, mpi_err)
    CALL par_handle_error(mpi_err)
  END SUBROUTINE par_sum_real_vector_at_one

  SUBROUTINE par_sum_real_matrix(matrix)
    INTEGER :: mpi_err
    REAL :: matrix(:, :)
    REAL, POINTER :: sum_result(:, :)

    ALLOCATE (sum_result(UBOUND(matrix, 1), UBOUND(matrix, 2)))
    CALL MPI_Allreduce(matrix, sum_result, &
      & UBOUND(matrix, 1) * UBOUND(matrix, 2), MPI_REAL, &
      & MPI_SUM, MPI_COMM_WORLD, mpi_err)
    CALL par_handle_error(mpi_err)
    matrix = sum_result
    DEALLOCATE (sum_result)
  END SUBROUTINE par_sum_real_matrix

  SUBROUTINE par_sum_int_vector(vector)
    INTEGER :: mpi_err
    INTEGER :: vector(:)
    INTEGER, POINTER :: sum_result(:)

    ALLOCATE (sum_result(LBOUND(vector, 1):UBOUND(vector, 1)))
    CALL MPI_Allreduce(vector, sum_result, UBOUND(vector, 1), MPI_INTEGER, &
      & MPI_SUM, MPI_COMM_WORLD, mpi_err)
    CALL par_handle_error(mpi_err)
    vector = sum_result
    DEALLOCATE (sum_result)
  END SUBROUTINE par_sum_int_vector

!
! Gathering
!
  SUBROUTINE par_gather_real(whole, part, stripes, offsets)
    REAL, POINTER :: whole(:), part(:)
    INTEGER, POINTER :: stripes(:), offsets(:)
    INTEGER :: mpi_err

    CALL MPI_Allgatherv(part, stripes(par_rank), MPI_REAL, &
      & whole, stripes, offsets, MPI_REAL, MPI_COMM_WORLD, mpi_err)
    CALL par_handle_error(mpi_err)
  END SUBROUTINE par_gather_real

!
! Striping
!
  FUNCTION par_get_stripes(n) RESULT (units)
    INTEGER :: n
    INTEGER, POINTER :: units(:)
    REAL, POINTER :: time(:)
    REAL :: total_power
    INTEGER :: i, proc, pos(1), extra_units, step

    ALLOCATE (time(0:par_processes-1), units(0:par_processes-1))

    total_power = SUM(1 / slowness_table)
    DO i = 0, par_processes - 1
      units(i) = NINT(n / slowness_table(i) / total_power)
    END DO
    time = units * slowness_table

    ! adjust for rounding errors
    extra_units = SUM(units) - n
    step = -SIGN(1, extra_units)
    DO i = 1, ABS(extra_units)
      IF (extra_units.GT.0) THEN
        pos = MAXLOC(time)
      ELSE
        pos = MINLOC(time + slowness_table)
      END IF
      proc = pos(1) - 1
      units(proc) = units(proc) + step
      time(proc) = time(proc) + step * slowness_table(proc)
    END DO
    DEALLOCATE (time)
  END FUNCTION par_get_stripes

  SUBROUTINE par_get_stripe(n, lower, upper, stripes, offsets)
    INTEGER :: n
    INTEGER, INTENT(OUT) :: lower, upper
    INTEGER, POINTER :: stripe_widths(:)
    INTEGER, POINTER, OPTIONAL :: stripes(:), offsets(:)
    INTEGER, POINTER :: stripe_bounds(:)
    INTEGER :: i

    stripe_widths => par_get_stripes(n)
    ALLOCATE (stripe_bounds(0:par_processes-1))
    IF (PRESENT(offsets)) THEN
      stripes => stripe_widths
      offsets => stripe_bounds
    END IF

    stripe_bounds(0) = 0
    DO i = 1, par_processes - 1
      stripe_bounds(i) = stripe_widths(i - 1) + stripe_bounds(i - 1) 
    END DO

    lower = stripe_bounds(par_rank) + 1
    IF (par_rank.EQ.par_processes - 1) THEN
      upper = n
    ELSE
      upper = stripe_bounds(par_rank + 1)
    END IF

    IF (.NOT.PRESENT(offsets)) THEN
      DEALLOCATE (stripe_bounds)
      DEALLOCATE (stripe_widths)
    END IF
  END SUBROUTINE par_get_stripe

  SUBROUTINE par_handle_error(mpi_err)
    INTEGER :: mpi_err

    IF (mpi_err.NE.0) THEN
      WRITE(0, "(A36, I4)") "WARNING: MPI Communication Error on ", par_rank
    END IF
  END SUBROUTINE par_handle_error
END MODULE parallel
