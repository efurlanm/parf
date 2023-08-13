MODULE bootstraps
  USE utilities
  USE options
  USE instancesets
  IMPLICIT NONE

  TYPE bootstrap
    ! For each instance, how many times it was taken into the current bootstrap
    INTEGER, POINTER :: instances(:)

    ! For each class, how many times it was taken into the current bootstrap
    INTEGER, POINTER :: classes(:)

    ! The indices of the bootstrapped samples
    INTEGER, POINTER :: idx(:)

    REAL, POINTER :: weight_sums_by_class(:)
    REAL, POINTER :: weight_sums_by_instance(:)
  END TYPE bootstrap

  ! For each instance, how many trees do not have it in its bootstrap
  INTEGER, POINTER :: out_of_bag(:)
CONTAINS
  SUBROUTINE init_bootstraps(isptr)
    TYPE (instanceset), POINTER :: isptr
    ALLOCATE (out_of_bag(UBOUND(isptr%catvars, 1)))
    out_of_bag = 0
  END SUBROUTINE init_bootstraps

  FUNCTION new_bootstrap(isptr) RESULT (bs)
    TYPE (instanceset), POINTER :: isptr
    TYPE (bootstrap), POINTER :: bs
    INTEGER :: instance_count, class_count, class_catvar, tagged_count
    INTEGER :: i, instance, num_bootstrapped, class, class_start
    REAL, ALLOCATABLE :: weights(:)

    instance_count = UBOUND(isptr%catvars, 1)
    class_count = isptr%dd%attributes(opts%class_attribute_num)%cat_count
    class_start = isptr%dd%attributes(opts%class_attribute_num)%cat_start
    class_catvar = isptr%dd%attributes(opts%class_attribute_num)%mapping

    ALLOCATE (bs, bs%instances(instance_count), bs%classes(class_count), &
      & bs%weight_sums_by_instance(instance_count), &
      & bs%weight_sums_by_class(class_count), weights(class_count))
    bs%instances = 0
    bs%classes = 0
    num_bootstrapped = 0
    tagged_count = UBOUND(isptr%tagged_instances, 1)
    DO i = 1, tagged_count
      ! Only consider the tagged instances for bootstraps
      instance = isptr%tagged_instances(rnd_int(tagged_count))
      class = isptr%catvars(instance, class_catvar)
      bs%instances(instance) = bs%instances(instance) + 1
      bs%classes(class) = bs%classes(class) + 1
    END DO

    DO i = 1, instance_count
      IF (bs%instances(i).EQ.0) THEN
        out_of_bag(i) = out_of_bag(i) + 1
      ELSE
        num_bootstrapped = num_bootstrapped + 1
      END IF
    END DO

    ALLOCATE (bs%idx(num_bootstrapped))
    instance = 0
    DO i = 1, instance_count
      IF (bs%instances(i).NE.0) THEN
        instance = instance + 1
        bs%idx(instance) = i
      END IF
    END DO

    weights = isptr%dd%categories(class_start:class_start+class_count-1)%weight
    DO i = 1, instance_count
      IF (isptr%catvars(i, class_catvar).EQ.missing_cat) CYCLE !GTDEBUG
      bs%weight_sums_by_instance(i) = &
        & bs%instances(i) * weights(isptr%catvars(i, class_catvar))
    END DO
    bs%weight_sums_by_class = bs%classes * weights
    DEALLOCATE (weights)
  END FUNCTION new_bootstrap

  SUBROUTINE free_bootstrap(bs)
    TYPE (bootstrap), POINTER :: bs
    IF (ASSOCIATED(bs)) THEN
      IF (ASSOCIATED(bs%classes)) DEALLOCATE (bs%classes)
      IF (ASSOCIATED(bs%instances)) DEALLOCATE (bs%instances)
      IF (ASSOCIATED(bs%idx)) DEALLOCATE (bs%idx)
      IF (ASSOCIATED(bs%weight_sums_by_class)) &
        & DEALLOCATE (bs%weight_sums_by_class)
      IF (ASSOCIATED(bs%weight_sums_by_instance)) &
        & DEALLOCATE (bs%weight_sums_by_instance)
      DEALLOCATE (bs)
    END IF
  END SUBROUTINE free_bootstrap

  SUBROUTINE finish_bootstraps()
    IF (ASSOCIATED(out_of_bag)) DEALLOCATE (out_of_bag)
  END SUBROUTINE finish_bootstraps
END MODULE bootstraps
