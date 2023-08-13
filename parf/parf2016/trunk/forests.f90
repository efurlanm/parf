MODULE forests
  USE instancesets
  USE bootstraps
  USE trees
  USE bitvectors
  USE importances
  USE parallel
  IMPLICIT NONE

  TYPE forest
    TYPE (tree), POINTER :: trees(:)
    REAL, POINTER :: dgini(:) ! average dgini for the forest
  END TYPE forest

CONTAINS
  FUNCTION new_forest(isptr) RESULT (rfptr)
    TYPE (instanceset), POINTER :: isptr
    TYPE (attribute), POINTER :: class
    TYPE (forest), POINTER :: rfptr
    TYPE (bootstrap), POINTER :: bsptr
    TYPE (node), POINTER :: leaf
    TYPE (tree), POINTER :: thetree
    INTEGER :: i, j, l, o, instance_count, attribute_count, pos(1)
    REAL :: avg_dgini, weight
    INTEGER, POINTER :: leaf_count(:) ! number of instances in each leaf
    INTEGER, POINTER :: node_tags(:)
    INTEGER :: class_attr
    REAL, POINTER :: node_weights(:)
    TYPE (bitvector), POINTER :: used_attrs
    REAL :: total_right, case_right
    INTEGER, POINTER :: seeds(:)
    INTEGER :: lower, upper

    instance_count = UBOUND(isptr%catvars, 1)
    attribute_count = UBOUND(isptr%dd%attributes, 1)
    class => isptr%dd%attributes(opts%class_attribute_num)
    class_attr = class%mapping

    CALL par_get_stripe(opts%num_trees, lower, upper)
    ALLOCATE (rfptr, rfptr%trees(lower:upper), rfptr%dgini(attribute_count), &
      & seeds(0:opts%num_trees))
    DO i = 0, opts%num_trees
      seeds(i) = rnd_integer()
    END DO
    ! TODO reallocation might not be necessary; investigate
    IF (.NOT.ASSOCIATED(isptr%estimated_class)) THEN
      ALLOCATE (isptr%estimated_class(instance_count), &
        & isptr%votes(instance_count, class%cat_count), &
        & isptr%leaf_id(instance_count, lower:upper))
    END IF
    isptr%estimated_class = 0
    isptr%votes = 0
! FIXME what is this NULLIFY line for?
    NULLIFY (isptr%outlying)
    rfptr%dgini = 0

    IF (opts%calc_importances) THEN
      used_attrs => new_bitvector(attribute_count)
    END IF

    DO i = lower, upper
      CALL seed_rnd(seeds(i))
      IF (opts%verboseall) THEN
        IF (par_processes.EQ.1) THEN
          WRITE(6, "(A14, I6)") "        Tree #", i
        ELSE
          WRITE(6, "(A14, I6, A4, I5)") "        Tree #", i, " on ", par_rank
        END IF
      END IF
      bsptr => new_bootstrap(isptr)
      thetree => rfptr%trees(i)
      thetree%oob => new_bitvector(instance_count)
      CALL build_tree(isptr, bsptr, thetree)

      ! check out-of-bag estimates, while we're dealing with a bootstrap
      ALLOCATE (thetree%dgini(attribute_count), &
        & thetree%leaf_bounds(0:thetree%terminal_bound), &
        & thetree%leaf_index(instance_count), &
        & thetree%leaf_pop(thetree%terminal_bound), &
        & leaf_count(thetree%terminal_bound))

      IF (opts%calc_importances) THEN
        imp%oob_count = instance_count - UBOUND(bsptr%idx, 1)
        ALLOCATE (imp%oob_index(imp%oob_count), &
          & imp%rnd_oob_index(imp%oob_count), &
          & node_tags(thetree%terminal_bound), &
          & node_weights(thetree%terminal_bound))
        total_right = 0
        used_attrs%bits = 0
        CALL fill_node_data(thetree%p, node_tags, node_weights, used_attrs)
      END IF

      thetree%leaf_pop = 0
      leaf_count = 0
      thetree%dgini = 0
      o = 0
      DO j = 1, instance_count
        ! collect fast-importance data (gini difference)
        leaf => classify_on_tree(isptr, j, thetree%p, thetree%dgini)

        l = leaf%id
        isptr%leaf_id(j, i) = l
        thetree%leaf_pop(l) = thetree%leaf_pop(l) + bsptr%instances(j)
        leaf_count(l) = leaf_count(l) + 1

        ! see if it's out-of-bag
        IF (bsptr%instances(j).EQ.0) THEN
          CALL setbit(thetree%oob, j, .TRUE.)

          weight = leaf%total_weight / leaf%num_instances
          isptr%votes(j, leaf%attribute) = &
            & isptr%votes(j, leaf%attribute) + weight

          ! collect the importances data
          IF (opts%calc_importances) THEN
            o = o + 1
            imp%oob_index(o) = j

            IF (leaf%attribute.EQ.isptr%catvars(j, class_attr)) THEN
              total_right = total_right + weight
              IF (opts%calc_case_importances) THEN
                imp%total_weights(j) = imp%total_weights(j) &
                  & + weight / imp%oob_count
              END IF
            END IF
          END IF
        END IF
      END DO

      IF (opts%calc_importances) THEN
        ! classify the set with randomly shuffled columns, see importances
        DO o = 1, attribute_count
          IF (getbit(used_attrs, o)) THEN
            CALL shuffle(imp%oob_index, imp%rnd_oob_index)
            case_right = 0
            DO j = 1, imp%oob_count
              leaf => classify_permuted_on_tree(isptr, j, thetree%p, o, &
                & imp%oob_index, imp%rnd_oob_index)
              l = imp%oob_index(j)
              IF (leaf%attribute.EQ.isptr%catvars(l, class_attr)) THEN
                weight = leaf%total_weight / leaf%num_instances
                case_right = case_right + weight
                IF (opts%calc_case_importances) THEN
                  imp%case_weights(l, o) = imp%case_weights(l, o) &
                    & + weight / imp%oob_count
                END IF
              END IF
            END DO
            imp%average_imp(o) = imp%average_imp(o) &
              & + (total_right - case_right) / imp%oob_count
            imp%squared_imp(o) = imp%squared_imp(o) &
              & + ((total_right - case_right) / imp%oob_count) ** 2
          ELSE IF (opts%calc_case_importances) THEN
            DO j = 1, instance_count
              l = isptr%leaf_id(j, i)
              IF (node_tags(l).EQ.isptr%catvars(j, class_attr)) THEN
                imp%case_weights(j, o) = imp%case_weights(j, o) &
                  & + node_weights(l) / imp%oob_count
              END IF
            END DO
          END IF
        END DO
      END IF
      ! Sort instances by leaf
      rfptr%trees(i)%dgini = rfptr%trees(i)%dgini / instance_count
      rfptr%dgini = rfptr%dgini + rfptr%trees(i)%dgini
      thetree%leaf_bounds(0) = 1
      DO j = 1, thetree%terminal_bound
        thetree%leaf_bounds(j) = leaf_count(j) + thetree%leaf_bounds(j - 1)
      END DO
      leaf_count = 0
      DO j = 1, instance_count
        l = isptr%leaf_id(j, i)
        thetree%leaf_index(thetree%leaf_bounds(l - 1) + leaf_count(l)) = j
        leaf_count(l) = leaf_count(l) + 1
      END DO

      ! Cleanup
      IF (opts%calc_importances) THEN
        DEALLOCATE (node_tags, node_weights, imp%oob_index, &
          & imp%rnd_oob_index)
      END IF
      DEALLOCATE (leaf_count)
      CALL free_bootstrap(bsptr)
    END DO

    ! Communicate the results
    CALL par_sum_int_vector(out_of_bag)
    CALL par_sum_real_vector(rfptr%dgini)
    CALL par_sum_real_matrix(isptr%votes)
    IF (opts%calc_importances) THEN
      CALL par_sum_real_vector(imp%average_imp)
      CALL par_sum_real_vector(imp%squared_imp)
      IF (opts%calc_case_importances) THEN
        CALL par_sum_real_vector(imp%total_weights)
        CALL par_sum_real_matrix(imp%case_weights)
      END IF
    END IF
    
    ! TODO Could be parallelized further
    DO j = 1, instance_count
      IF (out_of_bag(j).NE.0) THEN
        isptr%votes(j, :) = isptr%votes(j, :) / out_of_bag(j)
        pos = MAXLOC(isptr%votes(j, :))
        isptr%estimated_class(j) = pos(1)
      END IF
    END DO

    IF (opts%calc_importances) THEN
      CALL free_bitvector(used_attrs)
    END IF

    avg_dgini = SUM(rfptr%dgini) / attribute_count
    rfptr%dgini = rfptr%dgini / avg_dgini

    CALL seed_rnd(seeds(0))
    DEALLOCATE (seeds)
  END FUNCTION new_forest

  SUBROUTINE calc_training_error(trainset)
    TYPE (instanceset), POINTER :: trainset
    INTEGER :: i, oob_count, errorcount, class, instance_count, class_count
    REAL :: error, kappa_value
    LOGICAL :: never_oob

    IF (.NOT.par_front) RETURN
    oob_count = 0
    errorcount = 0
    instance_count = UBOUND(trainset%classes, 1)
    class = trainset%dd%attributes(opts%class_attribute_num)%mapping
    DO i = 1, instance_count
      never_oob = out_of_bag(i).EQ.0
      CALL setbit2(trainset%missing_data, i, opts%class_attribute_num, &
        & never_oob)
      IF (.NOT.never_oob.AND.trainset%catvars(i, class).NE.missing_cat) THEN
        oob_count = oob_count + 1
        IF (trainset%catvars(i, class).NE.trainset%estimated_class(i)) THEN
          errorcount = errorcount + 1
        END IF
      END IF
    END DO
    IF (opts%summary) THEN
      error = REAL(errorcount) / oob_count
      kappa_value = kappa(error)
      WRITE(6, "(A32, 1X, F6.2, A4, I8, 1X, A7, 1X, F6.4, A2)") &
        & "Trainset classification error is", &
        & REAL(errorcount) * 100 / oob_count, "% of", oob_count, "(kappa:", &
        & kappa_value, ")"
    END IF
  END SUBROUTINE calc_training_error

  FUNCTION kappa(error)
    REAL :: error, kappa
    REAL :: expected_miss
    INTEGER :: largest_class(1)

    largest_class = MAXVAL(trainset%class_populations)
    expected_miss = 1 - REAL(largest_class(1)) / UBOUND(trainset%classes, 1)
    kappa = (expected_miss - error) / expected_miss
  END FUNCTION kappa

  SUBROUTINE fix_num_prox(train_count)
    INTEGER :: train_count

    IF (opts%num_prox.LT.0) THEN
      opts%num_prox = -opts%num_prox * train_count / 100
    END IF
    IF (opts%num_prox.GT.train_count) THEN
      opts%num_prox = train_count
    END IF
    IF (opts%prox_for_prot.LT.0) THEN
      opts%prox_for_prot = -opts%prox_for_prot * train_count / 100
    END IF
    IF (opts%prox_for_prot.GT.opts%num_prox) THEN
      opts%prox_for_prot = opts%num_prox
      IF (par_front) &
        & WRITE(0, *) "Warning: Prototype search radius (-yp) cannot exceed&
          & the number of proximate instances to be considered (-p)."
    END IF
  END SUBROUTINE

  SUBROUTINE calculate_proximities(rfptr, isptr)
    TYPE (instanceset), POINTER :: isptr
    TYPE (forest), POINTER :: rfptr
    TYPE (tree), POINTER :: thetree
    INTEGER :: instance_count
    INTEGER :: i, k, l, alloc_err, proc

    TYPE (tree), POINTER :: trees(:)
    INTEGER, POINTER :: leaf_id(:, :)
    INTEGER :: main_proc, num_trees, pos(1), oobsize

    IF (opts%num_prox.EQ.0) RETURN
      
    NULLIFY (leaf_id, trees)
    instance_count = UBOUND(isptr%estimated_class, 1)
    CALL par_get_stripe(instance_count, isptr%lower, isptr%upper, &
      & isptr%par_stripes, isptr%par_offsets)
    ALLOCATE (isptr%prox(isptr%lower:isptr%upper, opts%num_prox), &
      & isptr%prox_index(isptr%lower:isptr%upper, opts%num_prox), &
      & STAT=alloc_err)
    ! TODO parallelise: reduce alloc_err
    IF (alloc_err.NE.0) THEN
      WRITE(0, *) "ERROR: Not enough memory for proximity calculation"
      WRITE(0, *) "       Proximity calculation skipped"
      WRITE(0, *) "Suggestion: Use -p parameter to lessen memory requirements"
      opts%num_prox = 0
      opts%num_prot = 0
      opts%num_scale = 0
      opts%fill_passes = 1
      SELECT CASE (isptr%set_type)
        CASE (trainset_type)
          opts%train_outliers = ""
        CASE (testset_type)
          opts%test_outliers = ""
      END SELECT
      RETURN
    END IF

    ! do we need to gather data to a node?
    IF (par_processes.EQ.1) THEN
      ! if we're running on only one process, then all this is unnecessary
      alloc_err = 1
    ELSE
      ! yes; find the fastest node
      pos = MAXLOC(isptr%par_stripes)
      main_proc = pos(1) - 1
      ! see if we can put all the relevant data there
      IF (par_rank.EQ.main_proc) THEN
        ALLOCATE (leaf_id(instance_count, opts%num_trees), &
          trees(opts%num_trees), STAT=alloc_err)
        IF (alloc_err.EQ.0) THEN
          oobsize = UBOUND(rfptr%trees(LBOUND(rfptr%trees, 1))%oob%bits, 1)
          ! clean the trees for clean deallocation in case of error
          DO i = 1, opts%num_trees
            thetree => trees(i)
            NULLIFY (thetree%leaf_bounds, thetree%leaf_pop, thetree%dgini, &
              & thetree%leaf_index, thetree%oob, thetree%p)
          END DO
        END IF

        i = 0
        DO proc = 0, par_processes - 1
          IF (proc.EQ.main_proc) THEN
            i = UBOUND(rfptr%trees, 1)
            CYCLE
          END IF
          CALL par_send_int(alloc_err, proc)
          IF (alloc_err.EQ.0) THEN
            CALL par_recv_int(num_trees)
            DO k = 1, num_trees
              i = i + 1
              CALL par_recv_int(l)
              thetree => trees(i)
              ALLOCATE (thetree%leaf_bounds(0:l), thetree%leaf_pop(l), &
                & thetree%leaf_index(instance_count), &
                & thetree%oob, thetree%oob%bits(oobsize), STAT=alloc_err)
              CALL par_send_int(alloc_err, proc)
              IF (alloc_err.EQ.0) THEN
                CALL par_recv_int_vector(leaf_id(:, i))
                CALL par_recv_int_vector(thetree%leaf_bounds)
                CALL par_recv_int_vector(thetree%leaf_pop)
                CALL par_recv_int_vector(thetree%leaf_index)
                CALL par_recv_int_vector(thetree%oob%bits)
              ELSE
                EXIT
              END IF
            END DO
          END IF
        END DO

        IF (alloc_err.EQ.0) THEN
          DO i = LBOUND(rfptr%trees, 1), UBOUND(rfptr%trees, 1)
            thetree => trees(i)
            l = UBOUND(rfptr%trees(i)%leaf_bounds, 1)
            ALLOCATE (thetree%leaf_bounds(0:l), thetree%leaf_pop(l), &
              & thetree%leaf_index(instance_count), &
              & thetree%oob, thetree%oob%bits(oobsize), STAT=alloc_err)
            IF (alloc_err.EQ.0) THEN
              leaf_id(:, i) = isptr%leaf_id(:, i)
              trees(i)%leaf_bounds = rfptr%trees(i)%leaf_bounds
              trees(i)%leaf_pop = rfptr%trees(i)%leaf_pop
              trees(i)%leaf_index = rfptr%trees(i)%leaf_index
              trees(i)%oob%bits = rfptr%trees(i)%oob%bits
            END IF
          END DO
        END IF
      ELSE
        CALL par_recv_int(alloc_err)
        IF (alloc_err.EQ.0) THEN
          CALL par_send_int(UBOUND(rfptr%trees, 1) - LBOUND(rfptr%trees, 1) + 1, &
            & main_proc)
          DO i = LBOUND(rfptr%trees, 1), UBOUND(rfptr%trees, 1)
            thetree => rfptr%trees(i)
            CALL par_send_int(UBOUND(rfptr%trees(i)%leaf_bounds, 1), main_proc)
            CALL par_recv_int(alloc_err)
            IF (alloc_err.EQ.0) THEN
              CALL par_send_int_vector(isptr%leaf_id(:, i), main_proc)
              CALL par_send_int_vector(thetree%leaf_bounds, main_proc)
              CALL par_send_int_vector(thetree%leaf_pop, main_proc)
              CALL par_send_int_vector(thetree%leaf_index, main_proc)
              CALL par_send_int_vector(thetree%oob%bits, main_proc)
            ELSE
              EXIT
            END IF
          END DO
        END IF
      END IF
    END IF

    IF (alloc_err.EQ.0) THEN
      IF (opts%verbose.AND.par_front) &
        & WRITE(6, af) "Performing centralised proximity calculation"
      CALL calculate_proximities_internal(isptr, trees, leaf_id, &
        & main_proc)
    END IF

    IF (par_rank.EQ.main_proc) THEN
      ! deallocate
      IF (ASSOCIATED(trees)) THEN
        DO i = 1, opts%num_trees
          CALL free_tree(trees(i))
        END DO
        DEALLOCATE (trees)
      END IF
      IF (ASSOCIATED(leaf_id)) DEALLOCATE (leaf_id)
    END IF

    IF (alloc_err.EQ.0) RETURN

    IF (opts%verbose.AND.par_front.AND.par_processes.GT.1) &
      & WRITE(6, af) "Performing distributed proximity calculation"
    ! failed to collect the berries; do the distributed approach
    CALL calculate_proximities_internal(isptr, rfptr%trees, &
      & isptr%leaf_id, -1)
  END SUBROUTINE calculate_proximities

  SUBROUTINE calculate_proximities_internal(isptr, trees, leaf_id, &
    & main_proc)
    TYPE (instanceset), POINTER :: isptr
    TYPE (tree), POINTER :: trees(:)
    INTEGER, POINTER :: leaf_id(:, :)
    INTEGER :: main_proc
    INTEGER :: j, proc, inst, i, k, l, instance, class
    INTEGER :: train_count, instance_count
    REAL :: prox_squares
    TYPE (tree), POINTER :: thetree
    REAL, POINTER :: curprox(:) ! proximities for the current instance
    REAL, POINTER :: sumcurprox(:) ! final proximities for the current
    INTEGER, POINTER :: prox_perm(:)
    REAL, POINTER :: outlying(:)
    LOGICAL :: is_oob, calc_outliers
    TYPE (attribute), POINTER :: class_attr
    INTEGER :: class_offset, class_mapping

    is_oob = .TRUE.
    SELECT CASE (isptr%set_type)
      CASE (trainset_type)
        calc_outliers = LEN_TRIM(opts%train_outliers).NE.0
      CASE (testset_type)
        calc_outliers = LEN_TRIM(opts%test_outliers).NE.0
      CASE DEFAULT
        calc_outliers = .FALSE.
    END SELECT
    IF (calc_outliers) ALLOCATE(outlying(isptr%lower:isptr%upper))

    class_attr => trainset%dd%attributes(opts%class_attribute_num)
    class_offset = class_attr%cat_start - 1
    class_mapping = class_attr%mapping
    train_count = UBOUND(trainset%classes, 1)
    instance_count = UBOUND(isptr%estimated_class, 1)
    ALLOCATE (curprox(train_count), sumcurprox(train_count))
    IF (opts%num_prox.NE.train_count.OR.opts%num_prot.NE.0) &
      & ALLOCATE (prox_perm(train_count))
    j = 0

    DO proc = 0, par_processes - 1
      DO inst = 1, isptr%par_stripes(proc)
        j = j + 1

        IF (main_proc.EQ.-1.OR.par_rank.EQ.main_proc) THEN
          curprox = 0
          DO i = LBOUND(trees, 1), UBOUND(trees, 1)
            thetree => trees(i)
            l = leaf_id(j, i)
            IF (isptr%set_type.EQ.trainset_type) is_oob = getbit(thetree%oob, j)
            IF (is_oob) THEN
              ! out-of-bag (or not in trainset)
              DO k = thetree%leaf_bounds(l - 1), thetree%leaf_bounds(l) - 1
                instance = thetree%leaf_index(k)
                IF (.NOT.getbit(thetree%oob, instance)) THEN
                  curprox(instance) = curprox(instance) &
                    & + trainset%dd%categories(class_offset &
                    & + trainset%classes(instance))%weight &
                    & / thetree%leaf_pop(leaf_id(instance, i))
                END IF
              END DO
            ELSE
              ! in bootstrap (only in trainset)
              DO k = thetree%leaf_bounds(l - 1), thetree%leaf_bounds(l) - 1
                instance = thetree%leaf_index(k)
                IF (getbit(thetree%oob, instance)) THEN
                  curprox(instance) = curprox(instance) &
                    & + trainset%dd%categories(class_offset &
                    & + trainset%classes(j))%weight / thetree%leaf_pop(l)
                END IF
              END DO
            END IF
          END DO
        END IF

        ! Collect the total row proximity values into sumcurprox
        IF (main_proc.EQ.-1) THEN
          CALL par_sum_real_vector_at_one(curprox, sumcurprox, proc)
        ELSE IF (main_proc.EQ.proc) THEN
          sumcurprox = curprox
        ELSE IF (par_rank.EQ.main_proc) THEN
          CALL par_send_real_vector(curprox, proc)
        ELSE IF (par_rank.EQ.proc) THEN
          CALL par_recv_real_vector(sumcurprox)
        END IF

        IF (proc.EQ.par_rank) THEN
          IF (calc_outliers) THEN
            prox_squares = 0
            class = isptr%classes(j)
            DO k = 1, train_count
              IF (trainset%classes(k).EQ.class) THEN
                prox_squares = prox_squares + sumcurprox(k) * sumcurprox(k)
              END IF
            END DO
            IF (prox_squares.EQ.0) prox_squares = 1
            ! this was instance_count/prox_squares, but now
            ! we divide it in print_outliers,
            ! and the dividend is class population
            outlying(j) = prox_squares
          END IF

          IF (opts%num_prox.EQ.train_count.AND.opts%num_prot.EQ.0) THEN
            isptr%prox_index(j, :) = (/(k, k = 1, train_count)/)
            isptr%prox(j, :) = sumcurprox
          ELSE
            CALL qsort_real(-sumcurprox, prox_perm)
            isptr%prox_index(j, :) = prox_perm
            isptr%prox(j, :) = sumcurprox(prox_perm)
          END IF
        END IF
      END DO
    END DO

    IF (opts%num_prox.NE.train_count.OR.opts%num_prot.NE.0) &
      & DEALLOCATE (prox_perm)
    DEALLOCATE (curprox, sumcurprox)

    IF (calc_outliers) THEN
      IF (par_front) ALLOCATE (isptr%outlying(instance_count))
      CALL par_gather_real(isptr%outlying, outlying, &
        & isptr%par_stripes, isptr%par_offsets)
      DEALLOCATE (outlying)
    END IF
  END SUBROUTINE calculate_proximities_internal

  RECURSIVE SUBROUTINE free_forest(rfptr)
    TYPE (forest), POINTER :: rfptr
    INTEGER :: i
    IF (ASSOCIATED(rfptr)) THEN
      IF (ASSOCIATED(rfptr%trees)) THEN
        DO i = LBOUND(rfptr%trees, 1), UBOUND(rfptr%trees, 1)
          CALL free_tree(rfptr%trees(i))
        END DO
        DEALLOCATE (rfptr%trees)
      END IF
      IF (ASSOCIATED(rfptr%dgini)) DEALLOCATE (rfptr%dgini)
      DEALLOCATE (rfptr)
    END IF
  END SUBROUTINE free_forest

  FUNCTION classify_on_forest(isptr, instance, rfptr) &
    & RESULT (voted_class)
    INTEGER :: voted_class
    TYPE (instanceset), POINTER :: isptr
    INTEGER :: instance
    INTEGER :: pos(1)
    TYPE (forest), POINTER :: rfptr
    TYPE (node), POINTER :: leaf
    INTEGER :: i, class

    DO i = LBOUND(rfptr%trees, 1), UBOUND(rfptr%trees, 1)
      leaf => classify_on_tree(isptr, instance, rfptr%trees(i)%p)
      class = leaf%attribute
      isptr%votes(instance, class) = isptr%votes(instance, class) &
        & + leaf%total_weight / leaf%num_instances
      isptr%leaf_id(instance, i) = leaf%id
    END DO
    CALL setbit2(isptr%missing_data, instance, opts%class_attribute_num, &
      & .FALSE.)
    
    CALL par_sum_real_vector(isptr%votes(instance, :))

    pos = MAXLOC(isptr%votes(instance, :))
    voted_class = pos(1)
  END FUNCTION classify_on_forest

  SUBROUTINE print_forest(rfptr, datadesc) ! TODO parallelise
    TYPE (forest), POINTER :: rfptr
    TYPE (datadescription), POINTER :: datadesc
    INTEGER :: handle, i
    LOGICAL :: err

    handle = open_file(opts%dump_forest, .TRUE., err)
    IF (err) RETURN

    IF (par_processes.GT.1) THEN
      IF (par_front) &
        WRITE(0, af) "Error: Forest dump not supported when distributed"
      RETURN
    END IF

    DO i = 1, UBOUND(rfptr%trees, 1)
      CALL print_tree(rfptr%trees(i)%p, datadesc, handle)
    END DO
    CALL close_file(handle)
  END SUBROUTINE print_forest

  SUBROUTINE save_forest(rfptr)
    TYPE (forest), POINTER :: rfptr
    TYPE (datadescription), POINTER :: datadesc
    INTEGER :: handle, i
    LOGICAL :: err
    CHARACTER(LEN=20) :: tmp, fmtstr

    datadesc => trainset%dd
    IF (par_front) THEN
      i = LEN_TRIM(opts%save_forest)
      IF (i.GT.7) THEN
        IF (opts%save_forest(i-6:i).EQ.".forest") &
          & opts%save_forest = opts%save_forest(1:i-7)
      END IF
      handle = open_file(TRIM(opts%save_forest) // ".forest", .TRUE., err)
    END IF
    CALL par_bcast_bool(err)
    IF (err) RETURN

    IF (par_front) THEN
      ! Data description
      WRITE(handle, *) opts%num_trees
      WRITE(handle, *) UBOUND(datadesc%attributes, 1), &
        & UBOUND(datadesc%categories, 1), &
        & opts%class_attribute_num, UBOUND(trainset%catvars, 1), &
        & UBOUND(datadesc%usedvars, 1)
      DO i = 1, UBOUND(datadesc%attributes, 1)
        WRITE(handle, *) datadesc%attributes(i)%name
        WRITE(handle, *) datadesc%attributes(i)%cat_count, &
          & datadesc%attributes(i)%cat_start, &
          & datadesc%attributes(i)%mapping
      END DO
      WRITE(handle, *) datadesc%usedvars
      DO i = 1, UBOUND(datadesc%categories, 1)
        WRITE(handle, *) datadesc%categories(i)%name
        WRITE(handle, *) datadesc%categories(i)%attribute, &
          & datadesc%categories(i)%weight
      END DO
      WRITE(handle, *) trainset%classes
      WRITE(handle, *) trainset%class_populations

      CALL close_file(handle)
    END IF

    i = FLOOR(LOG10(REAL(opts%num_trees))) + 1
    WRITE(fmtstr, "(A2, I1, A1, I1, A1)") "(I", i, ".", i, ")"
    DO i = LBOUND(rfptr%trees, 1), UBOUND(rfptr%trees, 1)
      WRITE(tmp, fmtstr) i
      err = .NOT.save_tree(rfptr%trees(i), TRIM(opts%save_forest) // "." &
        & // TRIM(ADJUSTL(tmp)) // ".tree", datadesc, trainset%leaf_id(:, i))
      IF (err) RETURN
    END DO
  END SUBROUTINE save_forest

  SUBROUTINE load_forest(rfptr, datadesc)
    TYPE (forest), POINTER :: rfptr
    TYPE (datadescription), POINTER :: datadesc
    INTEGER :: handle, i, lower, upper
    LOGICAL :: err
    INTEGER :: attribute_count, category_count, instance_count, used_count
    CHARACTER(LEN=20) :: tmp, fmtstr

    IF (par_front) THEN
      i = LEN_TRIM(opts%load_forest)
      IF (i.GT.7) THEN
        IF (opts%load_forest(i-6:i).EQ.".forest") &
          & opts%load_forest = opts%load_forest(1:i-7)
      END IF
      handle = open_file(TRIM(opts%load_forest) // ".forest", .FALSE., err)
    END IF
    CALL par_bcast_bool(err)
    IF (err) RETURN

    IF (par_front) THEN
      ! Data description
      READ(handle, *) opts%num_trees
      READ(handle, *) attribute_count, category_count, &
        & opts%class_attribute_num, instance_count, used_count
    END IF
    CALL par_bcast_int(opts%num_trees)
    CALL par_bcast_int(attribute_count)
    CALL par_bcast_int(category_count)
    CALL par_bcast_int(opts%class_attribute_num)
    CALL par_bcast_int(instance_count)
    CALL par_bcast_int(used_count)

    datadesc => new_datadescription()
    CALL par_get_stripe(opts%num_trees, lower, upper)
    ALLOCATE (datadesc%attributes(attribute_count), &
      & datadesc%categories(category_count), &
      & rfptr, rfptr%trees(lower:upper), &
      & trainset%estimated_class(instance_count), &
      & trainset%leaf_id(instance_count, lower:upper), &
      & datadesc%usedvars(used_count))
    NULLIFY (rfptr%dgini)

    datadesc%contvar_count = 0
    datadesc%catvar_count = 0
    DO i = 1, attribute_count
      IF (par_front) THEN
        READ(handle, *) datadesc%attributes(i)%name
        READ(handle, *) datadesc%attributes(i)%cat_count, &
          & datadesc%attributes(i)%cat_start, &
          & datadesc%attributes(i)%mapping
      END IF
      CALL par_bcast_char(datadesc%attributes(i)%name)
      CALL par_bcast_int(datadesc%attributes(i)%cat_count)
      CALL par_bcast_int(datadesc%attributes(i)%cat_start)
      CALL par_bcast_int(datadesc%attributes(i)%mapping)
      SELECT CASE (datadesc%attributes(i)%cat_count)
        CASE (ignorevar)
        CASE (contvar)
          datadesc%contvar_count = datadesc%contvar_count + 1
        CASE DEFAULT
          datadesc%catvar_count = datadesc%catvar_count + 1
      END SELECT
    END DO
    IF (par_front) READ(handle, *) datadesc%usedvars
    CALL par_bcast_int_vector(datadesc%usedvars)

    DO i = 1, category_count
      IF (par_front) THEN
        READ(handle, *) datadesc%categories(i)%name
        READ(handle, *) datadesc%categories(i)%attribute, &
          & datadesc%categories(i)%weight
      END IF
      CALL par_bcast_char(datadesc%categories(i)%name)
      CALL par_bcast_int(datadesc%categories(i)%attribute)
      CALL par_bcast_float(datadesc%categories(i)%weight)
    END DO

    ALLOCATE (trainset%class_populations(0:datadesc%attributes( &
      & opts%class_attribute_num)%cat_count))
    IF (par_front) THEN
      READ(handle, *) trainset%estimated_class
      READ(handle, *) trainset%class_populations
      CALL close_file(handle)
    END IF
    CALL par_bcast_int_vector(trainset%estimated_class)
    CALL par_bcast_int_vector(trainset%class_populations)
    trainset%dd => datadesc

    i = FLOOR(LOG10(REAL(opts%num_trees))) + 1
    WRITE(fmtstr, "(A2, I1, A1, I1, A1)") "(I", i, ".", i, ")"
    DO i = lower, upper
      WRITE(tmp, fmtstr) i
      err = .NOT.load_tree(rfptr%trees(i), TRIM(opts%load_forest) // "." &
        & // TRIM(ADJUSTL(tmp)) // ".tree", datadesc, instance_count, &
        & trainset%leaf_id(:, i))
      IF (err) RETURN
    END DO
  END SUBROUTINE load_forest

  SUBROUTINE classify_instanceset(isptr, rfptr)
    TYPE (instanceset), POINTER :: isptr
    TYPE (forest), POINTER :: rfptr
    TYPE (attribute), POINTER :: class
    INTEGER :: handle
    LOGICAL :: err, do_output
    INTEGER :: errorcount, class_attr, class_base, i, instance_count
    INTEGER :: voted_class, tagged_class
    CHARACTER(LEN=10) :: labelstr
    REAL :: error

    do_output = LEN_TRIM(opts%test_results).NE.0.AND.par_front
    IF (do_output) THEN
      handle = open_file(opts%test_results, .TRUE., err)
      IF (err) RETURN
    
      IF (isptr%tagged) THEN
        WRITE(handle, "(A10, 1X, A13, 17X, A9, 21X, A7)") &
        & "Row", "Classified as", "Tagged as", "Certain"
      ELSE
        WRITE(handle, "(A10, 1X, A13, 17X, A7)") &
        & "Row", "Classified as", "Certain"
      END IF
    END IF

    errorcount = 0
    class => isptr%dd%attributes(opts%class_attribute_num)
    class_attr = class%mapping
    class_base = class%cat_start - 1
    instance_count = UBOUND(isptr%catvars, 1)
    ALLOCATE (isptr%votes(instance_count, class%cat_count), &
      & isptr%estimated_class(instance_count), &
      & isptr%leaf_id(instance_count, opts%num_trees))
    isptr%votes = 0
    DO i = 1, instance_count
      voted_class = classify_on_forest(isptr, i, rfptr)
      isptr%estimated_class(i) = voted_class
      labelstr = label_for(isptr, i)
      IF (isptr%tagged) THEN
        tagged_class = isptr%catvars(i, class_attr)
        IF (tagged_class.NE.missing_cat &
          & .AND.tagged_class.NE.voted_class) THEN
          errorcount = errorcount + 1
          IF (do_output) &
            & WRITE(handle, "(A10, 1X, A30, A30, F6.2, A1)") labelstr, &
            & isptr%dd%categories(class_base + voted_class)%name, &
            & isptr%dd%categories(class_base + tagged_class)%name, &
            & isptr%votes(i, voted_class) * 100 / SUM(isptr%votes(i, :)), "%"
        END IF
      ELSE
        IF (do_output) &
          & WRITE(handle, "(A10, 1X, A30, F6.2, A1)") labelstr, &
          & isptr%dd%categories(class_base + voted_class)%name, &
          & isptr%votes(i, voted_class) * 100 / SUM(isptr%votes(i, :)), "%"
      END IF
    END DO
    IF (isptr%tagged.AND.opts%summary.AND.par_front) THEN
      error = REAL(errorcount) / instance_count
      WRITE(6, "(A32, 1X, F6.2, A4, I8, 1X, A7, 1X, F6.4, A2)") &
        & "Testset classification error is", &
        & error * 100, "% of", instance_count, &
        & "(kappa:", kappa(error), ")"
    END IF
    
    IF (do_output) CALL close_file(handle)
  END SUBROUTINE classify_instanceset

  SUBROUTINE print_interaction(rfptr, dd) ! TODO parallelise
    TYPE (forest), POINTER :: rfptr
    TYPE (datadescription), POINTER :: dd
    INTEGER :: handle, attribute_count, used_count, t, attr
    INTEGER :: i, j, ii, jj, iu, ju
    LOGICAL :: err
    CHARACTER(LEN=20) :: fmtstr
    INTEGER, POINTER :: ranks(:, :), attr_perm(:)
    REAL, POINTER :: effect(:, :), hist(:, :)
    REAL :: rcor, teffect

    handle = open_graphics(opts%interaction, err, "interactions")
    IF (err) RETURN

    IF (opts%do_graphics) THEN
      WRITE(graphics_handle, af) "set xtics 1"
      WRITE(graphics_handle, af) "set grid x"
      WRITE(graphics_handle, af) "set ytics 1"
      WRITE(graphics_handle, af) "set grid y"
    END IF
    attribute_count = UBOUND(dd%attributes, 1)
    used_count = UBOUND(dd%usedvars, 1)
    CALL squeeze_num(fmtstr, attribute_count + 0.5)
    CALL plot_graphics("splot [0.5:" // TRIM(fmtstr) // "][0.5:" &
      & // TRIM(fmtstr) // "] '" // TRIM(opts%interaction) &
      & // "' title 'Interactions' with lines");

    WRITE(fmtstr, *) attribute_count
    fmtstr = "(" // TRIM(ADJUSTL(fmtstr)) // "I8)"
    ALLOCATE (ranks(used_count, opts%num_trees), &
      & attr_perm(used_count), hist(used_count, used_count), &
      & effect(attribute_count, attribute_count))

    ranks = 0 ! ranks(attribute, tree) are attribute ranks for each tree
    hist = 0  ! hist(attribute, rank) is frequency of a rank in in an attribute
    DO t = 1, opts%num_trees
      CALL qsort_real(rfptr%trees(t)%dgini(dd%usedvars), attr_perm)
      j = 0
      DO i = used_count, 1, -1
        attr = attr_perm(i)
        ! FIXME following line probably dgini LE 0
        IF (ABS(rfptr%trees(t)%dgini(dd%usedvars(attr))).LT.8.232D-11) THEN
          EXIT
        END IF

        j = j + 1
        ranks(attr, t) = j

        hist(attr, i) = hist(attr, i) + 1
      END DO
    END DO
    hist = hist / opts%num_trees

    effect = 0
    DO i = 1, used_count
      iu = dd%usedvars(i)
      DO j = 1, used_count
        IF (i.EQ.j) CYCLE
        ju = dd%usedvars(j)
        teffect = 0
        DO ii = 1, used_count
          DO jj = 1, used_count
            teffect = teffect + IABS(ii - jj) * hist(i, jj) * hist(j, ii)
          END DO
        END DO
        rcor = SUM(hist(i, :) * hist(j, :))
        effect(iu, ju) = (SUM(IABS(ranks(i, :) - ranks(j, :))) &
          & / opts%num_trees - teffect / (1 - rcor)) * 100
      END DO
    END DO

    IF (opts%do_graphics) THEN
      DO i = 1, used_count
        ii = dd%usedvars(i)
        IF (i.NE.1) WRITE(handle, *)
        DO j = 1, used_count
          jj = dd%usedvars(j)
          WRITE(handle, "(I6, I6, F10.2)") ii, jj, effect(ii, jj)
        END DO
      END DO
    ELSE
      WRITE(handle, fmtstr) NINT(effect)
    END IF

    DEALLOCATE (ranks, attr_perm, effect, hist)

    CALL close_graphics(handle)
  END SUBROUTINE print_interaction

  SUBROUTINE print_outliers(isptr)
    TYPE (instanceset), POINTER :: isptr
    INTEGER, POINTER :: class_populations(:), class_index(:), perm(:)
    TYPE (attribute), POINTER :: class_attr
    INTEGER :: i, j, class
    INTEGER :: instance_count, handle
    REAL :: med, dev
    LOGICAL :: err
    CHARACTER(LEN=20) :: fmtstr
    CHARACTER(LEN=optlen) :: filename

    IF (isptr%set_type.EQ.testset_type) THEN
      fmtstr = "outliers-test"
      filename = opts%test_outliers
    ELSE
      fmtstr = "outliers-train"
      filename = opts%train_outliers
    END IF

    handle = open_graphics(filename, err, fmtstr)
    IF (err) RETURN
    
    instance_count = UBOUND(isptr%classes, 1)
    class_attr => isptr%dd%attributes(opts%class_attribute_num)
    ALLOCATE (class_populations(class_attr%cat_count), perm(instance_count))

    ! count the instances for each class
    class_populations = 0
    DO i = 1, instance_count
      class = isptr%classes(i)
      class_populations(class) = class_populations(class) + 1
    END DO

    WRITE(handle, af), "    Median  Deviation Tag"
    ! calculate the outliers for each class
    DO class = 1, UBOUND(class_populations, 1)
      ALLOCATE (class_index(class_populations(class)))

      j = 0
      DO i = 1, instance_count
        IF (isptr%classes(i).EQ.class) THEN
          j = j + 1
          class_index(j) = i
        END IF
      END DO

! FIXME moved the division here, and using class_population to normalise
      isptr%outlying(class_index) = &
        & class_populations(class) / isptr%outlying(class_index)
      
!     med = median(isptr%outlying(class_index))
! FIXME i put in average instead of median
med = SUM(isptr%outlying(class_index)) / j

      DO i = 1, j
        ! FIXME why cut at 5med?
        dev = dev + MIN(ABS(isptr%outlying(class_index(i)) - med), 5 * med)
      END DO
      dev = dev / j

      WRITE(handle, '(E10.4, 1X, E10.4, 1X, A30)'), med, dev, &
        & isptr%dd%categories(class_attr%cat_start + class - 1)%name

      DO i = 1, j
        ! FIXME the original had MIN(...., 20.0) - don't see why, dropped it
        isptr%outlying(class_index(i)) = &
          & (isptr%outlying(class_index(i)) - med) / dev
      END DO

      DEALLOCATE (class_index)
    END DO

    CALL qsort_real(isptr%outlying, perm)

    WRITE(handle, af), "  Row#    Outlier Tagged as"
    DO i = instance_count, 1, -1
      IF (isptr%outlying(perm(i)).LT.opts%outlier_cutoff) CYCLE
      WRITE(handle, '(I6, 1X, F10.4, 1X, A30)'), &
        & perm(i), isptr%outlying(perm(i)), &
        & isptr%dd%categories(class_attr%cat_start + &
        & isptr%classes(perm(i)) - 1)%name
    END DO

    DEALLOCATE (class_populations, perm)
    CALL close_graphics(handle)
  END SUBROUTINE print_outliers

  SUBROUTINE calc_scaling()
    INTEGER :: x, i, s, div
    REAL, POINTER :: y(:), avg_instance_prox(:), u_wh(:), ev(:, :), u(:)
    REAL, POINTER :: y_test(:), avg_instance_prox_test(:)
    REAL, POINTER :: u_test(:), ev_test(:, :)
    REAL, POINTER :: y_proto(:), avg_instance_prox_proto(:)
    REAL, POINTER :: u_proto(:), ev_proto(:, :)
    REAL, POINTER :: bl(:), dl(:)
    REAL :: y2, eu, ru, ra, sa, avg_prox, ynorm
    REAL :: avg_prox_test, avg_prox_proto
    INTEGER :: train_count, test_count, proto_count
    TYPE (attribute), POINTER :: class
    LOGICAL :: do_testset, do_protoset
    REAL :: max_crit, crit, last_crit

    class => trainset%dd%attributes(opts%class_attribute_num)

    do_testset = LEN_TRIM(opts%test_scaling).NE.0
    do_protoset = LEN_TRIM(opts%proto_scaling).NE.0

    train_count = UBOUND(trainset%classes, 1)
    ALLOCATE (y(trainset%lower:trainset%upper), &
      & u(trainset%lower:trainset%upper), &
      & avg_instance_prox(trainset%lower:trainset%upper), &
      & ev(trainset%lower:trainset%upper, opts%num_scale), &
      & bl(opts%num_scale), dl(opts%num_scale), &
      & trainset%scaling(trainset%lower:trainset%upper, opts%num_scale), &
      & u_wh(train_count))
    trainset%scaling = 0

    DO i = trainset%lower, trainset%upper
      avg_instance_prox(i) = SUM(trainset%prox(i, :)) / opts%num_prox
    END DO
    avg_prox = SUM(avg_instance_prox) / train_count
    CALL par_sum_real(avg_prox)

    IF (do_testset) THEN
      test_count = UBOUND(testset%classes, 1)
      ALLOCATE (y_test(testset%lower:testset%upper), &
        & u_test(testset%lower:testset%upper), &
        & avg_instance_prox_test(testset%lower:testset%upper), &
        & ev_test(testset%lower:testset%upper, opts%num_scale), &
        & testset%scaling(testset%lower:testset%upper, opts%num_scale))

      DO i = testset%lower, testset%upper
        avg_instance_prox_test(i) = SUM(testset%prox(i, :)) / opts%num_prox
      END DO
      avg_prox_test = SUM(avg_instance_prox_test) / test_count
      CALL par_sum_real(avg_prox_test)
    END IF

    IF (do_protoset) THEN
      proto_count = UBOUND(protoset%classes, 1)
      ALLOCATE (y_proto(protoset%lower:protoset%upper), &
        & u_proto(protoset%lower:protoset%upper), &
        & avg_instance_prox_proto(protoset%lower:protoset%upper), &
        & ev_proto(protoset%lower:protoset%upper, opts%num_scale), &
        & protoset%scaling(protoset%lower:protoset%upper, opts%num_scale))

      DO i = protoset%lower, protoset%upper
        avg_instance_prox_proto(i) = SUM(protoset%prox(i, :)) / opts%num_prox
      END DO
      avg_prox_proto = SUM(avg_instance_prox_proto) / proto_count
      CALL par_sum_real(avg_prox_proto)
    END IF

    outer_loop: DO s = 1, opts%num_scale
      IF (opts%verbose) WRITE(6, "(A14, I6)") "  Coordinate #", s
      ! y = (-1, 1, -1, 1...)
      y = (/(1 - 2 * MOD(i, 2), i = trainset%lower, trainset%upper)/)
      IF (do_testset) THEN
        y_test = (/(1 - 2 * MOD(i, 2), i = testset%lower, testset%upper)/)
      END IF
      IF (do_protoset) THEN
        y_proto = (/(1 - 2 * MOD(i, 2), i = protoset%lower, protoset%upper)/)
      END IF

      x = 1
      max_crit = 0
      last_crit = HUGE(max_crit)
      DO
        y2 = SUM(y * y)
        CALL par_sum_real(y2)
        y2 = SQRT(y2)
        u = y / y2
        CALL par_gather_real(u_wh, u, &
          & trainset%par_stripes, trainset%par_offsets)
        IF (do_testset) u_test = y_test / y2
        IF (do_protoset) u_proto = y_proto / y2
        eu = SUM(u_wh)
        ru = SUM(avg_instance_prox * u)
        CALL par_sum_real(ru)
        DO i = trainset%lower, trainset%upper
          y(i) = (SUM(trainset%prox(i, :) * u_wh(trainset%prox_index(i, :))) &
            & - (avg_instance_prox(i) - avg_prox) * eu - ru) / 2
        END DO
        DO i = 1, s - 1
          bl(i) = SUM(ev(:, i) * u)
          CALL par_sum_real(bl(i))
          y = y - bl(i) * dl(i) * ev(:, i)
        END DO

        IF (do_testset) THEN
          DO i = testset%lower, testset%upper
            y_test(i) = (SUM(testset%prox(i, :) * &
              & u_wh(testset%prox_index(i, :))) &
              & - (avg_instance_prox_test(i) - avg_prox_test) * eu - ru) / 2
          END DO
          DO i = 1, s - 1
            y_test = y_test - bl(i) * dl(i) * ev_test(:, i)
          END DO
        END IF

        IF (do_protoset) THEN
          DO i = protoset%lower, protoset%upper
            y_proto(i) = (SUM(protoset%prox(i, :) * &
              & u_wh(protoset%prox_index(i, :))) &
              & - (avg_instance_prox_proto(i) - avg_prox_proto) * eu - ru) / 2
          END DO
          DO i = 1, s - 1
            y_proto = y_proto - bl(i) * dl(i) * ev_proto(:, i)
          END DO
        END IF

        ra = SUM(y * u)
        CALL par_sum_real(ra)
        ynorm = SUM((y - ra * u) ** 2)
        CALL par_sum_real(ynorm)
        sa = ABS(ra)
        IF (ynorm.LT.sa * 1.0e-7) THEN
          dl(s) = ra
          trainset%scaling(trainset%lower:trainset%upper, s) = u * SQRT(sa)
          ev(:, s) = u
          IF (do_testset) THEN
            testset%scaling(testset%lower:testset%upper, s) = u_test * SQRT(sa)
            ev_test(:, s) = u_test
          END IF
          IF (do_protoset) THEN
            protoset%scaling(protoset%lower:protoset%upper, s) &
              & = u_proto * SQRT(sa)
            ev_proto(:, s) = u_proto
          END IF
          EXIT
        END IF
        crit = -LOG(sa * 1.0e-7 / ynorm)
        IF (max_crit.LT.crit) THEN
          max_crit = crit
        ELSE IF (crit.GT.max_crit) THEN
          crit = max_crit
        END IF
        IF (opts%verbose) THEN
          WRITE(6, "(I6, 1X, F8.2, A1)") x, (1 - crit / max_crit) * 100, "%"
        END IF
        x = x + 1
        IF (last_crit.LT.crit) THEN
          div = div + 1
          IF (div.GT.opts%scaling_divergence) THEN
            IF (par_front) THEN
              WRITE(0, af) "ERROR: The scaling calculation diverges"
              WRITE(0, af) "Suggestion: Allow for scaling calculation&
                & divergence (-sd)"
              WRITE(0, af) "            Use wider proximity matrix, if&
                & shortened (-p)"
              WRITE(0, af) "            Increase the number of trees (-n)"
            END IF
            opts%train_scaling = ""
            IF (do_testset) opts%test_scaling = ""
            IF (do_protoset) opts%proto_scaling = ""
            EXIT outer_loop
          END IF
        ELSE
          div = 0
        END IF
        last_crit = crit
      END DO
    END DO outer_loop

    IF (do_testset) THEN
      DEALLOCATE (y_test, u_test, avg_instance_prox_test, ev_test)
    END IF
    IF (do_protoset) THEN
      DEALLOCATE (y_proto, u_proto, avg_instance_prox_proto, ev_proto)
    END IF
    DEALLOCATE (y, u, u_wh, avg_instance_prox, ev, bl, dl)
  END SUBROUTINE calc_scaling

  SUBROUTINE print_scaling(isptr)
    TYPE (instanceset), POINTER :: isptr
    CHARACTER(LEN=optlen) :: filename
    INTEGER :: s
    TYPE (attribute), POINTER :: class
    CHARACTER(LEN=20) :: fmtstr
    CHARACTER(LEN=linelen) :: plot
    CHARACTER :: width
    INTEGER :: handle
    LOGICAL :: err

    class => isptr%dd%attributes(opts%class_attribute_num)

    IF (par_front) THEN
      SELECT CASE (isptr%set_type)
        CASE (trainset_type)
          fmtstr = "scaling-train"
          filename = opts%train_scaling
        CASE (testset_type)
          fmtstr = "scaling-test"
          filename = opts%test_scaling
        CASE (protoset_type)
          fmtstr = "scaling-proto"
          filename = opts%proto_scaling
      END SELECT

      handle = open_graphics(filename, err, fmtstr)
      IF (err) RETURN

      IF (opts%do_graphics) THEN
        IF (opts%num_scale.LT.3) THEN
          plot = "plot"
        ELSE
          plot = "splot"
        END IF

        SELECT CASE (isptr%set_type)
          CASE (trainset_type)
            IF (LEN_TRIM(opts%test_scaling).EQ.0) THEN
              width = "1"
            ELSE
              width = "0"
            END IF
          CASE (testset_type)
            width = "1"
          CASE (protoset_type)
            width = "2"
        END SELECT

        DO s = 1, class%cat_count
          IF (s.NE.1) plot = TRIM(plot) // ","
          plot = TRIM(plot) // " '" // TRIM(filename) // "' using 2:3"
          IF (opts%num_scale.GT.2) plot = TRIM(plot) // ":4"
          IF (.NOT.graphics_together) THEN
            WRITE(fmtstr, *) s - 1
            plot = TRIM(plot) // " index " // TRIM(ADJUSTL(fmtstr))
          END IF
          WRITE(fmtstr, *) s
          plot = TRIM(plot) // " pt " // TRIM(ADJUSTL(fmtstr)) &
            & // " ps " // width // " title '" &
            & // TRIM(isptr%dd%categories(class%cat_start + s - 1)%name) // "'"
        END DO
        CALL plot_graphics(plot)
      END IF

      WRITE(fmtstr, *) opts%num_scale
      fmtstr = "(I6, 1X, " // TRIM(ADJUSTL(fmtstr)) // "F10.3)"
    END IF

    IF (opts%do_graphics) THEN
      IF (par_front) WRITE(handle, af) "#  Row Coordinates..."
      DO s = 1, class%cat_count
        IF (par_front) THEN
          IF (.NOT.graphics_together.AND.s.NE.1) THEN
            WRITE(handle, *)
            WRITE(handle, *)
          ELSE IF (s.NE.1) THEN
            WRITE(handle, af) "end"
          END IF
          WRITE(handle, af) "# " &
            & // TRIM(isptr%dd%categories(class%cat_start + s - 1)%name)
        END IF
        CALL print_scaling_rows(handle, fmtstr, isptr, s)
      END DO
    ELSE
      IF (par_front) WRITE(handle, af) "  Row# Coordinates..."
      CALL print_scaling_rows(handle, fmtstr, isptr)
    END IF

    IF (par_front) CALL close_graphics(handle)
  END SUBROUTINE print_scaling

  SUBROUTINE print_scaling_rows(handle, fmtstr, isptr, s)
    INTEGER :: handle
    CHARACTER(LEN=*) :: fmtstr
    TYPE (instanceset), POINTER :: isptr
    INTEGER, OPTIONAL :: s
    LOGICAL :: print_all, none_written
    REAL, POINTER :: scaling_row(:)
    INTEGER :: i, row, proc, class

    print_all = .NOT.PRESENT(s)
    IF (print_all) THEN
      class = -1
    ELSE
      class = s
    END IF

    IF (par_front) THEN
      ! front node
      none_written = .TRUE.
      DO i = isptr%lower, isptr%upper
        IF (print_all.OR.isptr%classes(i).EQ.class) THEN
          WRITE(handle, fmtstr) i, isptr%scaling(i, :)
          none_written = .FALSE.
        END IF
      END DO
      IF (none_written.AND..NOT.print_all) WRITE(handle, af) "$"
      
      ! other nodes
      IF (par_processes.GT.1) THEN
        none_written = .TRUE.
        ALLOCATE (scaling_row(opts%num_scale))
        i = isptr%upper
        DO proc = 1, par_processes - 1
          CALL par_notify(proc)
          DO row = 1, isptr%par_stripes(proc)
            i = i + 1
            IF (print_all.OR.isptr%classes(i).EQ.class) THEN
              CALL par_recv_real_vector(scaling_row)
              WRITE(handle, fmtstr) i, scaling_row
              none_written = .FALSE.
            END IF
          END DO
        END DO
        DEALLOCATE (scaling_row)
        IF (none_written.AND..NOT.print_all) WRITE(handle, af) "$"
      END IF
    ELSE
      CALL par_wait()
      DO i = isptr%lower, isptr%upper
        IF (print_all.OR.isptr%classes(i).EQ.class) THEN
          CALL par_send_real_vector(isptr%scaling(i, :), 0)
        END IF
      END DO
    END IF
  END SUBROUTINE print_scaling_rows
END MODULE forests
