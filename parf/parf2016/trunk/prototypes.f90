MODULE prototypes
  USE options
  USE instancesets
  USE utilities
  IMPLICIT NONE

  TYPE prototype
    INTEGER :: num_close
    TYPE (prototype_contvar), POINTER :: contvars(:)
    TYPE (prototype_catvar), POINTER :: catvars(:)
  END TYPE prototype
  
  TYPE prototype_contvar
    REAL :: vmed, vhigh, vlow
    REAL :: med, high, low
  END TYPE prototype_contvar

  TYPE prototype_catvar
    INTEGER, POINTER :: freq(:)
    INTEGER :: class
  END TYPE prototype_catvar

  TYPE (prototype), POINTER :: prots(:,:)
CONTAINS
  SUBROUTINE free_prototype(prot)
    TYPE (prototype) :: prot
    INTEGER :: i

    DO i = 1, UBOUND(prot%catvars, 1)
      IF (prot%catvars(i)%class.NE.0) DEALLOCATE (prot%catvars(i)%freq)
    END DO

    DEALLOCATE (prot%catvars, prot%contvars)
  END SUBROUTINE free_prototype

  SUBROUTINE free_prototypes(dd)
    TYPE (datadescription), POINTER :: dd
    INTEGER :: i, j

    IF (opts%num_prot.EQ.0) RETURN
    DO i = 1, dd%attributes(opts%class_attribute_num)%cat_count
      DO j = 1, opts%num_prot
        IF (prots(i, j)%num_close.EQ.0) EXIT
        CALL free_prototype(prots(i, j))
      END DO
    END DO
    DEALLOCATE (prots)
    NULLIFY (protoset%estimated_class)
    IF (ASSOCIATED(protoset)) CALL free_instanceset(protoset)
  END SUBROUTINE free_prototypes

  SUBROUTINE calculate_prototypes()
    TYPE (prototype), POINTER :: prot
    LOGICAL, POINTER :: unseen(:)
    INTEGER :: num_close, max_close, max_close_loc
    TYPE (attribute), POINTER :: attrib, class_attr
    INTEGER :: c, p, i, j, ii, attr, proto_count
    INTEGER, POINTER :: near_index(:), near_perm(:), prox_copy(:)
    REAL :: dt
    LOGICAL :: here
    
    class_attr => trainset%dd%attributes(opts%class_attribute_num)
    ALLOCATE (unseen(UBOUND(trainset%classes, 1)))
    unseen = .TRUE.
    IF (par_front) THEN
      ALLOCATE (prots(class_attr%cat_count, opts%num_prot), &
        & prox_copy(opts%prox_for_prot))
      prots(:,:)%num_close = 0
      proto_count = 0
    END IF

    DO c = 1, class_attr%cat_count
      unseen = .TRUE.
      DO p = 1, opts%num_prot
        IF (p.NE.1) THEN
          ! on each subsequent pass, update and distribute unseen
          unseen(max_close_loc) = .FALSE.
          unseen(prox_copy) = .FALSE.
          CALL par_bcast_bool_vector(unseen)
        END IF
        ! find the unseen instance with the most yet unseen neighbours
        ! classified as c
        max_close = 0
        max_close_loc = 0
        DO i = trainset%lower, trainset%upper
          IF (unseen(i)) THEN
            num_close = 0
            DO j = 1, opts%prox_for_prot
              ii = trainset%prox_index(i, j)
              IF (unseen(ii).AND.trainset%classes(ii).EQ.c) THEN
                num_close = num_close + 1
              END IF
            END DO
            IF (num_close.GT.max_close) THEN
              max_close = num_close
              max_close_loc = i
            END IF
          END IF
        END DO
        i = max_close_loc
        CALL par_max_int(max_close, max_close_loc)
        IF (max_close.EQ.0) EXIT ! no more prototypes for this class?
        here = max_close_loc.EQ.i

        ! get the proximity vector for the "friendliest" instance
        IF (par_front) THEN
          IF (here) THEN
            prox_copy = trainset%prox_index(max_close_loc, &
              & 1:opts%prox_for_prot)
          ELSE
            CALL par_recv_int_vector(prox_copy)
          END IF
        ELSE IF (here) THEN
          CALL par_send_int_vector(trainset%prox_index(max_close_loc, &
            & 1:opts%prox_for_prot), 0)
        END IF

        ! calculate the prototype data
        IF (par_front) THEN
          proto_count = proto_count + 1
          prot => prots(c, p)
          prot%num_close = max_close
          ALLOCATE (prot%contvars(trainset%dd%contvar_count), &
            & prot%catvars(trainset%dd%catvar_count))
          prot%catvars(:)%class = 0

          ! identify the closest same-class instances among the unseen
          ALLOCATE (near_index(max_close), near_perm(max_close))
          j = 0
          DO i = 1, opts%prox_for_prot
            ii = prox_copy(i)
            IF (unseen(ii).AND.trainset%classes(ii).EQ.c) THEN
              j = j + 1
              near_index(j) = ii
            END IF
          END DO

          DO i = 1, UBOUND(trainset%dd%usedvars, 1)
            attrib => trainset%dd%attributes(trainset%dd%usedvars(i))
            attr = attrib%mapping
            SELECT CASE (attrib%cat_count)
              CASE (ignorevar)
                ! nothing
              CASE (contvar)
                ! continuous variable
                dt = trainset%contvar_95_percent(attr) &
                  & - trainset%contvar_5_percent(attr)
                CALL qsort_real(trainset%contvars(near_index, attr), near_perm)
                near_perm = near_index(near_perm)
                ii = near_perm(NINT((max_close - 1) / 4.) + 1)
                prot%contvars(attr)%vlow = trainset%contvars(ii, attr)
                prot%contvars(attr)%low = (trainset%contvars(ii, attr) &
                  & - trainset%contvar_5_percent(attr)) / dt
                ii = near_perm(NINT((max_close - 1) * 3 / 4.) + 1)
                prot%contvars(attr)%vhigh = trainset%contvars(ii, attr)
                prot%contvars(attr)%high = (trainset%contvars(ii, attr) &
                  & - trainset%contvar_5_percent(attr)) / dt
                ii = near_perm(NINT((max_close - 1) / 2.) + 1)
                prot%contvars(attr)%vmed = trainset%contvars(ii, attr)
                prot%contvars(attr)%med = (trainset%contvars(ii, attr) &
                  & - trainset%contvar_5_percent(attr)) / dt
              CASE DEFAULT
                ! categorical attribute
                ALLOCATE (prot%catvars(attr)%freq(attrib%cat_count))
                prot%catvars(attr)%class = most_frequent(trainset%catvars( &
                  & prox_copy, attr), attrib%cat_count, prot%catvars(attr)%freq)
            END SELECT
          END DO

          DEALLOCATE (near_index, near_perm)
        END IF
      END DO
    END DO

    DEALLOCATE (unseen)
    IF (par_front) THEN
      CALL prototype_instanceset(proto_count)
      protoset%relation_name = TRIM(trainset%relation_name) // "-proto"
      DEALLOCATE (prox_copy)
    END IF
  END SUBROUTINE calculate_prototypes

  SUBROUTINE prototype_instanceset(proto_count)
    INTEGER :: proto_count
    TYPE (prototype), POINTER :: prot
    INTEGER :: c, p, a, n, attribute_count, used_count
    TYPE (attribute), POINTER :: class_attr

    protoset => new_instanceset(protoset_type)
    protoset%dd => trainset%dd
    attribute_count = UBOUND(trainset%dd%attributes, 1)
    used_count = UBOUND(trainset%dd%usedvars, 1)
    class_attr => trainset%dd%attributes(opts%class_attribute_num)

    ALLOCATE (protoset%contvars(proto_count, trainset%dd%contvar_count), &
      & protoset%catvars(proto_count, trainset%dd%catvar_count))
    protoset%missing_data => new_bitvector2(proto_count, attribute_count)
    protoset%missing_data%bits = -1 ! all true - all missing
    protoset%estimated_class => protoset%catvars(:, class_attr%mapping)
    protoset%classes => protoset%catvars(:, class_attr%mapping)
    n = 0
    DO c = 1, class_attr%cat_count
      DO p = 1, opts%num_prot
        prot => prots(c, p)
        IF (prot%num_close.EQ.0) EXIT
        n = n + 1
        protoset%contvars(n, :) = prot%contvars(:)%vmed
        protoset%catvars(n, :) = prot%catvars(:)%class
        DO a = 1, used_count
          CALL setbit2(protoset%missing_data, n, trainset%dd%usedvars(a), &
            & .FALSE.)
        END DO
        protoset%catvars(n, class_attr%mapping) = c
        CALL setbit2(protoset%missing_data, n, opts%class_attribute_num, &
          & .FALSE.)
      END DO
    END DO
  END SUBROUTINE prototype_instanceset

  SUBROUTINE print_prototype_analysis()
    INTEGER :: c, p, i, j, attr
    CHARACTER(LEN=20) :: tmp
    CHARACTER(LEN=linelen) :: str
    TYPE (attribute), POINTER :: attrib, class_attr
    TYPE (prototype), POINTER :: prot
    CHARACTER :: star
    INTEGER :: handle
    LOGICAL :: err

    handle = open_file(opts%prototype_analysis, .TRUE., err)
    IF (err) RETURN

    class_attr => trainset%dd%attributes(opts%class_attribute_num)
    DO c = 1, class_attr%cat_count
      WRITE(handle, *) "Class: ", &
        & TRIM(trainset%dd%categories(class_attr%cat_start + c - 1)%name)
      DO p = 1, opts%num_prot
        prot => prots(c, p)
        IF (prot%num_close.EQ.0) EXIT

        CALL squeeze_num(tmp, REAL(p))
        str = "    Prototype #" // TRIM(tmp)
        CALL squeeze_num(tmp, REAL(prot%num_close))
        str = TRIM(str) // " (" // TRIM(tmp) // ")"
        WRITE(handle, *) TRIM(str)

        DO i = 1, UBOUND(trainset%dd%usedvars, 1)
          attrib => trainset%dd%attributes(trainset%dd%usedvars(i))
          attr = attrib%mapping
          SELECT CASE (attrib%cat_count)
            CASE (ignorevar)
              ! nothing
            CASE (contvar)
              ! continuous variable
              WRITE(handle, *) "        ", TRIM(attrib%name)
              WRITE(handle, "(15X, 3F10.2)") prot%contvars(attr)%vlow, &
                & prot%contvars(attr)%vmed, prot%contvars(attr)%vhigh
              WRITE(handle, "(15X, 3F10.2)") prot%contvars(attr)%low, &
                & prot%contvars(attr)%med, prot%contvars(attr)%high
            CASE DEFAULT
              ! categorical variable
              WRITE(handle, *) "        ", TRIM(attrib%name)
              DO j = 1, attrib%cat_count
                WRITE(tmp, *) prot%catvars(attr)%freq(j)
                IF (j.EQ.prot%catvars(attr)%class) THEN
                  star = "*"
                ELSE
                  star = " "
                END IF
                WRITE(handle, *) "              ", star, " ", &
                  & TRIM(trainset%dd%categories(attrib%cat_start &
                  & + j - 1)%name), " (", TRIM(ADJUSTL(tmp)), ")"
              END DO
          END SELECT
        END DO
      END DO
    END DO
    CALL close_file(handle)
  END SUBROUTINE print_prototype_analysis
END MODULE prototypes
