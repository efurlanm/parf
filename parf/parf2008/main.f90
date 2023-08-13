PROGRAM random_forest
  USE options
  USE instancesets
  USE bootstraps
  USE utilities
  USE forests
  USE prototypes
  USE parallel
  IMPLICIT NONE

  TYPE (datadescription), POINTER :: datadesc
  TYPE (forest), POINTER :: rfptr
  INTEGER :: fill_pass
  LOGICAL :: last_pass

! MPI ...
  CALL par_init()

  IF (parse_options()) THEN
    CALL init_graphics()
    NULLIFY (testset, trainset, protoset, datadesc, rfptr)
    trainset => new_instanceset(trainset_type)
    IF (LEN_TRIM(opts%load_forest).GT.0) THEN
      IF (opts%verbose) WRITE(6, af) "Loading forest"
      CALL load_forest(rfptr, datadesc)
      trainset%classes => trainset%estimated_class
      CALL fix_num_prox(UBOUND(trainset%estimated_class, 1))
      IF (LEN_TRIM(opts%testset).GT.0) THEN
        IF (opts%verbose) WRITE(6, af) "Loading test set"
        testset => new_instanceset(testset_type)
        testset%dd => datadesc
        IF (.NOT.parse_arff(testset, opts%testset)) GO TO 9999

        IF (opts%verbose) WRITE(6, af) "Classifying testing set"
        CALL classify_instanceset(testset, rfptr)
        testset%classes => testset%estimated_class
      END IF

      IF (opts%last_prox_required) THEN
        IF (opts%verbose) WRITE(6, af) "Calculating proximities"
        CALL calculate_proximities(rfptr, trainset)
        IF (opts%calc_test_prox) THEN
          CALL calculate_proximities(rfptr, testset)
        END IF
      END IF
    ELSE
      IF (opts%verbose .AND. par_processes.EQ.1) THEN
        WRITE(6, af) "Loading training set"
      ELSE IF (opts%verbose.AND.par_front) THEN
        WRITE(6, af) "Loading and distributing training set"
      END IF
      IF (.NOT.parse_arff(trainset, opts%trainset)) GO TO 9999
      CALL fix_num_prox(UBOUND(trainset%catvars, 1))
      datadesc => trainset%dd
      trainset%classes => trainset%catvars(:, &
        & datadesc%attributes(opts%class_attribute_num)%mapping)
      IF (LEN_TRIM(opts%testset).GT.0) THEN
        IF (opts%verbose) WRITE(6, af) "Loading test set"
        testset => new_instanceset(testset_type)
        testset%dd => datadesc
        IF (.NOT.parse_arff(testset, opts%testset)) GO TO 9999
      END IF
      IF (opts%verbose) THEN
        WRITE(6, "(A26, I6)") "Number of training cases: ", &
          & UBOUND(trainset%catvars, 1)
        WRITE(6, "(A26, I6)") "Number of attributes:     ", &
          & UBOUND(datadesc%attributes, 1)
      END IF

      ! Prelude
      IF (opts%verbose) WRITE(6, af) "Counting classes"
      CALL count_classes(trainset)
      IF (opts%fill_passes.NE.0) THEN
        IF (opts%verbose) WRITE(6, af) "Calculating rough fills"
        CALL calculate_rough_fills(trainset)
        IF (opts%verbose) WRITE(6, af) "Filling missing values"
        CALL fill_missing_rough(trainset)
      END IF
      CALL allocate_importance_arrays(trainset)
      CALL init_bootstraps(trainset)

      ! Variations
      DO
        CALL get_num_split_variables(datadesc)
        IF (opts%verbose) THEN
          WRITE(6, "(A26, I6)") "Number of used attributes:", &
            & UBOUND(datadesc%usedvars, 1)
          WRITE(6, "(A26, I6)") "Attributes to split on:   ", &
            & opts%split_variables
        END IF
        CALL zero_importance_arrays()
        fill_pass = 1
        DO WHILE (fill_pass.LE.MAX(1, opts%fill_passes)) ! at least 1 pass
          ! this is a while loop and not a for loop,
          ! to allow early exit in case proximities can't be calculated

          IF (opts%verbose.AND.opts%fill_passes.GT.1) &
            & WRITE(6, "(A6, I2)") "Pass #", fill_pass
          IF (opts%verbose) WRITE(6, af) "Sorting and ranking"
          CALL sort_and_rank(trainset, fill_pass.GT.1)

          IF (opts%verbose) WRITE(6, af) "Growing forest"
          rfptr => new_forest(trainset)

          last_pass = fill_pass.GE.opts%fill_passes &
            & .AND.opts%redo_with_important_vars.EQ.0 &
            & .AND.opts%redo_with_significant_vars.EQ.0

          IF (last_pass) THEN
            CALL calc_training_error(trainset)
            IF (LEN_TRIM(opts%testset).NE.0) THEN
              IF (opts%verbose) WRITE(6, af) "Classifying testing set"
              CALL classify_instanceset(testset, rfptr)
              testset%classes => testset%estimated_class
            END IF
          END IF

          IF (fill_pass.LT.opts%fill_passes.OR.opts%last_prox_required) THEN
            IF (opts%verbose) WRITE(6, af) "Calculating proximities"
            CALL calculate_proximities(rfptr, trainset)
            IF (opts%calc_test_prox.AND.last_pass) THEN
              ! test set proximities only on the very last pass
              CALL calculate_proximities(rfptr, testset)
            END IF
          END IF

          IF (fill_pass.NE.MAX(1, opts%fill_passes)) THEN ! each pass but last
            IF (opts%verbose) WRITE(6, af) "Filling missing values"
            CALL fill_missing_by_prox(trainset)
            CALL free_forest(rfptr)
          END IF

          fill_pass = fill_pass + 1
        END DO

        ! redo with most important variables?
        CALL finalize_importance_arrays(trainset)
        IF (opts%redo_with_important_vars.NE.0) THEN
          opts%redo_with_important_vars = 0 ! redo just once
        ELSE
          EXIT
        END IF
      END DO

      IF (LEN_TRIM(opts%save_forest).NE.0) THEN
        IF (opts%verbose) WRITE(6, af) "Saving forest"
        CALL save_forest(rfptr)
      END IF
    END IF

    ! Finale

    IF (opts%num_prot.NE.0) THEN
      IF (opts%verbose) WRITE(6, af) "Calculating prototypes"
      CALL calculate_prototypes()
    END IF
    IF (opts%num_scale.NE.0) THEN
      IF (opts%verbose) WRITE(6, af) "Calculating scaling coordinates"
      IF (LEN_TRIM(opts%proto_scaling).NE.0) THEN
        CALL classify_instanceset(protoset, rfptr)
        CALL calculate_proximities(rfptr, protoset)
      END IF
      CALL calc_scaling()
    END IF

    IF (par_front) THEN
      IF (LEN_TRIM(opts%train_votes).NE.0) THEN
        IF (opts%verbose) WRITE(6, af) "Printing training set votes"
        CALL print_votes(trainset, opts%train_votes)
      END IF
      IF (LEN_TRIM(opts%train_confusion).NE.0 &
          & .OR.LEN_TRIM(opts%positive_category).NE.0) THEN
        IF (opts%verbose) WRITE(6, af) "Processing training set confusion matrix"
        CALL process_confusion_matrix(trainset, opts%train_confusion)
      END IF
      IF (LEN_TRIM(opts%fast_importances).NE.0) THEN
        IF (opts%verbose) WRITE(6, af) "Printing fast variable importances"
        CALL print_fast_importances(rfptr%dgini, datadesc)
      END IF
      IF (LEN_TRIM(opts%importances).NE.0) THEN
        IF (opts%verbose) WRITE(6, af) "Printing variable importances"
        CALL print_importances(datadesc)
      END IF
      IF (LEN_TRIM(opts%case_importances).NE.0) THEN
        IF (opts%verbose) &
          & WRITE(6, af) "Printing case-by-case variable importances"
        CALL print_case_importances(trainset)
      END IF
      IF (LEN_TRIM(opts%interaction).NE.0) THEN
        IF (opts%verbose) WRITE(6, af) "Printing variable interaction"
        CALL print_interaction(rfptr, datadesc)
      END IF
      IF (LEN_TRIM(opts%prototype_analysis).NE.0) THEN
        IF (opts%verbose) WRITE(6, af) "Printing prototype analysis"
        CALL print_prototype_analysis()
      END IF
      IF (LEN_TRIM(opts%prototypes).NE.0) THEN
        IF (opts%verbose) WRITE(6, af) "Printing prototypes"
        CALL print_arff(opts%prototypes, protoset)
      END IF
      IF (LEN_TRIM(opts%train_outliers).NE.0) THEN
        IF (opts%verbose) WRITE(6, af) "Printing training set outliers"
        CALL print_outliers(trainset)
      END IF
      IF (LEN_TRIM(opts%test_outliers).NE.0) THEN
        IF (opts%verbose) WRITE(6, af) "Printing test set outliers"
        CALL print_outliers(testset)
      END IF
      IF (LEN_TRIM(opts%test_votes).NE.0) THEN
        IF (opts%verbose) WRITE(6, af) "Printing test set votes"
        CALL print_votes(testset, opts%test_votes)
      END IF
      IF (LEN_TRIM(opts%test_arff).NE.0) THEN
        IF (opts%verbose) WRITE(6, af) "Printing test set ARFF"
        CALL print_arff(opts%test_arff, testset)
      END IF
      IF (LEN_TRIM(opts%train_test_arff).NE.0) THEN
        IF (opts%verbose) WRITE(6, af) "Printing train+test set ARFF"
        CALL print_arff(opts%train_test_arff, trainset, testset)
      END IF
      IF (LEN_TRIM(opts%test_confusion).NE.0 &
          & .OR.(LEN_TRIM(opts%positive_category).NE.0 &
          & .AND.LEN_TRIM(opts%testset).NE.0)) THEN
        IF (opts%verbose) WRITE(6, af) "Processing test set confusion matrix"
        CALL process_confusion_matrix(testset, opts%test_confusion)
      END IF
    END IF
    IF (LEN_TRIM(opts%train_scaling).NE.0) THEN
      IF (opts%verbose) WRITE(6, af) "Printing training set scaling coordinates"
      CALL print_scaling(trainset)
    END IF
    IF (LEN_TRIM(opts%test_scaling).NE.0) THEN
      IF (opts%verbose) WRITE(6, af) "Printing test set scaling coordinates"
      CALL print_scaling(testset)
    END IF
    IF (LEN_TRIM(opts%proto_scaling).NE.0) THEN
      IF (opts%verbose) WRITE(6, af) "Printing prototype scaling coordinates"
      CALL print_scaling(protoset)
    END IF
    IF (LEN_TRIM(opts%dump_forest).NE.0) THEN
      IF (opts%verbose) WRITE(6, af) "Printing forest"
      CALL print_forest(rfptr, datadesc)
    END IF

    9999 CONTINUE
    IF (par_front) CALL free_prototypes(datadesc)
    CALL free_importance_arrays()
    CALL free_forest(rfptr)
    CALL free_instanceset(trainset)
    CALL free_instanceset(testset)
    CALL free_datadescription(datadesc)
    CALL finish_bootstraps()
    CALL finish_graphics()
  END IF

! MPI ...
  IF (opts%verbose.AND.par_processes.GT.1) WRITE(6, af) "Finalizing"
  CALL par_finalize()
  IF (opts%verbose.AND.par_processes.GT.1) WRITE(6, af) "Finished"
END PROGRAM random_forest
