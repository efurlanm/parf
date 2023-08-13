MODULE importances
  USE instancesets
  USE graphics

  IMPLICIT NONE

  TYPE importances_details
    ! global
    REAL, POINTER :: total_weights(:), case_weights(:,:)
    REAL, POINTER :: average_imp(:), squared_imp(:)
    REAL, POINTER :: z_score(:), significance(:)
    INTEGER, POINTER :: attr_perm(:)
    ! reusable
    INTEGER :: oob_count
    INTEGER, POINTER :: oob_index(:), rnd_oob_index(:)
  END TYPE importances_details
  TYPE (importances_details) :: imp

  REAL :: minreal
CONTAINS
  SUBROUTINE print_fast_importances(dgini, dd)
    REAL, POINTER :: dgini(:)
    TYPE (datadescription), POINTER :: dd
    INTEGER :: handle
    LOGICAL :: err
    INTEGER, POINTER :: variable_importance_chart(:)
    INTEGER :: i, j, attribute_count
    CHARACTER(LEN=20) :: tmp

    handle = open_graphics(opts%fast_importances, err, "fast_importances")
    IF (err) RETURN

    IF (opts%do_graphics) THEN
      WRITE(graphics_handle, af) "set xtics 1"
      WRITE(graphics_handle, af) "set xlabel 'Attribute'"
      WRITE(graphics_handle, af) "set ylabel 'Importance'"
      WRITE(graphics_handle, af) "set boxwidth 0.9 relative"
    END IF

    attribute_count = UBOUND(dd%attributes, 1)
    CALL squeeze_num(tmp, attribute_count+0.5)
    CALL plot_graphics("plot [0.5:"// TRIM(tmp) //"] '" &
      & // TRIM(opts%fast_importances) &
      & // "' using 1:2 title 'Fast Importances' with boxes fs solid")

    ALLOCATE (variable_importance_chart(attribute_count))
    IF (opts%do_graphics) THEN
      variable_importance_chart = (/(i, i = attribute_count, 1, -1)/)
      WRITE(handle, af), "# Class    dGINI"
    ELSE
      CALL qsort_real(dgini, variable_importance_chart)
      WRITE(handle, af), " dGINI Tag"
    END IF
    DO i = attribute_count, 1, -1
      j = variable_importance_chart(i)
      IF (j.EQ.0) CYCLE
      IF (opts%do_graphics) THEN
        WRITE(handle, '(3X, I4, 1X, F8.4)'), j, dgini(j)
      ELSE
        IF (dgini(j).EQ.0) EXIT
        WRITE(handle, '(F6.2, 1X, A30)'), dgini(j), dd%attributes(j)%name
      END IF
    END DO
    CALL close_graphics(handle)
    DEALLOCATE (variable_importance_chart)
  END SUBROUTINE print_fast_importances

  SUBROUTINE allocate_importance_arrays(isptr)
    TYPE (instanceset), POINTER :: isptr
    INTEGER :: instance_count, attribute_count
    
    IF (opts%calc_importances) THEN
      instance_count = UBOUND(isptr%catvars, 1)
      attribute_count = UBOUND(isptr%dd%attributes, 1)

      ALLOCATE (imp%average_imp(attribute_count), &
        & imp%squared_imp(attribute_count))

      IF (opts%calc_case_importances) THEN
        ALLOCATE (imp%total_weights(instance_count), &
          & imp%case_weights(instance_count, attribute_count))
      END IF

      NULLIFY (imp%z_score, imp%significance, imp%attr_perm)
    END IF
  END SUBROUTINE allocate_importance_arrays

  SUBROUTINE zero_importance_arrays()
    IF (opts%calc_importances) THEN
      imp%average_imp = 0
      imp%squared_imp = 0
      IF (opts%calc_case_importances) THEN
        imp%total_weights = 0
        imp%case_weights = 0
      END IF
    END IF
  END SUBROUTINE zero_importance_arrays

  SUBROUTINE finalize_importance_arrays(isptr)
    TYPE (instanceset), POINTER :: isptr
    INTEGER :: attribute_count, i, j, c
    REAL :: temp, std_error_sq

    minreal = -HUGE(minreal)
    IF (opts%calc_importances) THEN
      attribute_count = UBOUND(imp%average_imp, 1)
      ALLOCATE (imp%z_score(attribute_count), &
        & imp%significance(attribute_count), imp%attr_perm(attribute_count))

      DO i = 1, attribute_count
        temp = imp%average_imp(i) / opts%num_trees ! X/N
        std_error_sq = (imp%squared_imp(i) / opts%num_trees - temp ** 2) &
          & / opts%num_trees
        IF (std_error_sq.GT.0) THEN
          imp%z_score(i) = temp / SQRT(std_error_sq)
          imp%significance(i) = erfcc(imp%z_score(i))
        ELSE
          imp%z_score(i) = minreal
          imp%significance(i) = 1
        END IF
      END DO
      CALL qsort_real(imp%z_score, imp%attr_perm)
      IF (opts%redo_with_significant_vars.GT.0) THEN
        opts%redo_with_important_vars = 0

        DO WHILE (imp%significance(imp%attr_perm( &
            & attribute_count - opts%redo_with_important_vars)) &
            & .LE.opts%redo_with_significant_vars)
          opts%redo_with_important_vars = opts%redo_with_important_vars + 1
          IF (opts%redo_with_important_vars.GE.attribute_count) EXIT
        END DO
      END IF
      opts%redo_with_significant_vars = 0

      ! select the most important variables, if needed
      IF (opts%redo_with_important_vars.GT.0.AND. &
        & UBOUND(isptr%dd%usedvars, 1).GT.opts%redo_with_important_vars) THEN
        DEALLOCATE (isptr%dd%usedvars)
        ALLOCATE (isptr%dd%usedvars(opts%redo_with_important_vars))
        j = 0
        c = 0
        DO i = attribute_count, &
            & attribute_count - opts%redo_with_important_vars + 1, -1
          j = j + 1
          isptr%dd%usedvars(j) = imp%attr_perm(i)
          IF (isptr%dd%attributes(imp%attr_perm(i))%cat_count.EQ.contvar) &
            & c = c + 1
        END DO
        DEALLOCATE (isptr%dd%usedconts)
        ALLOCATE (isptr%dd%usedconts(c))
        c = 0
        DO i = attribute_count, &
            & attribute_count - opts%redo_with_important_vars + 1, -1
          IF (isptr%dd%attributes(imp%attr_perm(i))%cat_count.NE.contvar) CYCLE
          c = c + 1
          isptr%dd%usedconts(c) = imp%attr_perm(i)
        END DO
        CALL index_usedconts(isptr%dd)
        opts%split_variables = opts%split_variables_important
      END IF
    END IF
  END SUBROUTINE finalize_importance_arrays

  SUBROUTINE free_importance_arrays()
    IF (opts%calc_importances) THEN
      IF (opts%calc_case_importances) THEN
        DEALLOCATE (imp%total_weights, imp%case_weights)
      END IF
      DEALLOCATE (imp%average_imp, imp%squared_imp)
      IF (ASSOCIATED(imp%z_score)) DEALLOCATE (imp%z_score)
      IF (ASSOCIATED(imp%significance)) DEALLOCATE (imp%significance)
      IF (ASSOCIATED(imp%attr_perm)) DEALLOCATE (imp%attr_perm)
    END IF
  END SUBROUTINE free_importance_arrays

  SUBROUTINE print_importances(dd)
    TYPE (datadescription), POINTER :: dd
    INTEGER :: i, j, handle, attribute_count
    LOGICAL :: err
    CHARACTER(LEN=20) :: tmp

    handle = open_graphics(opts%importances, err, "importances")
    IF (err) RETURN

    attribute_count = UBOUND(dd%attributes, 1)

    IF (opts%do_graphics) THEN
      WRITE(graphics_handle, af) "set xlabel 'Attribute'"
      WRITE(graphics_handle, af) "set xtics 1"
      WRITE(graphics_handle, af) "set ylabel 'Importance'"
      WRITE(graphics_handle, af) "set y2label 'Significance'"
      WRITE(graphics_handle, af) "set log y2"
      WRITE(graphics_handle, af) "set y2tics 0,0.1"

      CALL squeeze_num(tmp, attribute_count+0.5)
      CALL plot_graphics("plot [0.5:"// TRIM(tmp) //"] '" &
        & // TRIM(opts%importances) &
        & // "' using 1:2:($2/$3) title 'Importance' axis x1y1 with&
        & errorbars, '" // TRIM(opts%importances) // "' using 1:4 axis x1y2&
        & title 'Significance' with points ps 2 pt 5, '" &
        & // TRIM(opts%importances) // "' using 1:3:(0.9) title 'Z-Score'&
        & with boxes")

      WRITE(handle, af) "# Class    Imp   Z-Sc Significan"
      CALL print_importances_graphics(attribute_count, handle)
      IF (graphics_together) THEN ! twice more
        WRITE(handle, af) "end"
        CALL print_importances_graphics(attribute_count, handle)
        WRITE(handle, af) "end"
        CALL print_importances_graphics(attribute_count, handle)
      END IF
    ELSE
      WRITE(handle, af) "  Imp   Z-Sc Significan Tag"
      DO i = attribute_count, 1, -1
        j = imp%attr_perm(i)
        IF (imp%z_score(j).EQ.minreal) EXIT
        WRITE(handle, '(F6.2, 1X, F6.2, 1X, F10.6, 1X, A30)'), &
          & 100 * imp%average_imp(j) / opts%num_trees, &
          & imp%z_score(j), imp%significance(j), dd%attributes(j)%name
      END DO
    END IF
    CALL close_graphics(handle)
  END SUBROUTINE print_importances

  SUBROUTINE print_importances_graphics(attribute_count, handle)
    INTEGER :: attribute_count, handle, j

    DO j = 1, attribute_count
      IF (imp%z_score(j).EQ.minreal) CYCLE
      WRITE(handle, '(1X, I6, 1X, F6.2, 1X, F6.2, 1X, F10.6)'), &
        & j, 100 * imp%average_imp(j) / opts%num_trees, &
        & imp%z_score(j), imp%significance(j)
    END DO
  END SUBROUTINE print_importances_graphics

  SUBROUTINE print_case_importances(isptr)
    TYPE (instanceset), POINTER :: isptr
    INTEGER :: i, handle, attribute_count
    LOGICAL :: err
    CHARACTER(LEN=20) :: fmtstr
    CHARACTER(LEN=10) :: labelstr

    handle = open_graphics(opts%case_importances, err, 'case_importances')
    IF (err) RETURN

    IF (opts%do_graphics) THEN
      WRITE(graphics_handle, af) "set xtics 1"
      WRITE(graphics_handle, af) "set grid x"
      WRITE(graphics_handle, af) "set pm3d"
      labelstr = "0"
    ELSE
      labelstr = "10"
    END IF
    attribute_count = UBOUND(imp%case_weights, 2)
    CALL squeeze_num(fmtstr, attribute_count+0.5)
    CALL plot_graphics("splot [0.5:" // TRIM(fmtstr) // "] '" &
      & // TRIM(opts%case_importances) &
      & // "' matrix title 'Case-by-Case Importance' with pm3d");

    WRITE(fmtstr, *) attribute_count
    fmtstr = "(A" // TRIM(labelstr) // ", " // TRIM(ADJUSTL(fmtstr)) // "F10.2)"
    labelstr = ""
    DO i = 1, UBOUND(imp%case_weights, 1)
      IF (.NOT.opts%do_graphics) THEN
        labelstr = label_for(isptr, i)
      END IF
      WRITE(handle, fmtstr), labelstr, &
        & 100 * (imp%total_weights(i) - imp%case_weights(i,:)) / opts%num_trees
    END DO
    CALL close_graphics(handle)
  END SUBROUTINE print_case_importances
END MODULE importances
