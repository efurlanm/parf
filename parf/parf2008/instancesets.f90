MODULE instancesets
  USE bitvectors
  USE options
  USE utilities
  USE parallel
  IMPLICIT NONE

  INTEGER, PARAMETER :: ignorevar = 0, contvar = 1, catvar = 2 ! and more
  INTEGER, PARAMETER :: parse_ok = 0, parse_comment = 1, &
    & parse_error = 2, parse_end = 3
  INTEGER, PARAMETER :: missing_cat = 0
  REAL :: missing_cont
  DATA missing_cont/B'11111111111111111111111111111111'/

  TYPE attribute
    CHARACTER(LEN=tokenlen) :: name
    INTEGER :: cat_count ! number of categories, or contvar or ignorevar
    INTEGER :: cat_start ! start index in categories, for catvars
    INTEGER :: mapping   ! index in the appropriate *vars array
    INTEGER :: usedcont  ! order in usedconts
  END TYPE attribute

  TYPE category
    CHARACTER(LEN=tokenlen) :: name
    INTEGER :: attribute
    REAL :: weight
  END TYPE

  TYPE instanceset
    INTEGER :: set_type
    CHARACTER(LEN=tokenlen) :: relation_name
    REAL, POINTER :: contvars(:,:)
    INTEGER, POINTER :: catvars(:,:)
    CHARACTER(LEN=tokenlen), POINTER :: ignorevars(:,:)
    TYPE (bitvector), POINTER :: missing_data
    INTEGER, POINTER :: missing_count(:)
    INTEGER, POINTER :: contvar_perms(:,:)
    INTEGER, POINTER :: contvar_ranks(:,:)
    REAL, POINTER :: contvar_5_percent(:)
    REAL, POINTER :: contvar_95_percent(:)
    INTEGER, POINTER :: class_populations(:)
    INTEGER, POINTER :: tagged_instances(:)
    REAL, POINTER :: prox(:, :)
    INTEGER, POINTER :: prox_index(:, :)
    REAL, POINTER :: outlying(:)
    INTEGER, POINTER :: leaf_id(:, :)
    REAL, POINTER :: scaling(:, :)
    TYPE (bitvector), POINTER :: constant_attr
    TYPE (datadescription), POINTER :: dd
    LOGICAL :: tagged
    INTEGER, POINTER :: estimated_class(:)
    REAL, POINTER :: votes(:, :)
    INTEGER, POINTER :: classes(:) ! not allocated, not deallocated
    INTEGER, POINTER :: par_stripes(:), par_offsets(:) ! instance distribution
    INTEGER :: lower, upper
  END TYPE instanceset

  TYPE datadescription
    TYPE (attribute), POINTER :: attributes(:)
    TYPE (category), POINTER :: categories(:)
    REAL, POINTER :: contfills(:)
    INTEGER, POINTER :: catfills(:)
    INTEGER, POINTER :: usedvars(:), usedconts(:)
    INTEGER :: contvar_count, catvar_count
  END TYPE datadescription

  TYPE (instanceset), POINTER :: trainset, testset, protoset
  TYPE (attribute), POINTER :: label
  INTEGER :: labelattr
CONTAINS
  FUNCTION new_datadescription() RESULT (ddptr)
    TYPE (datadescription), POINTER :: ddptr
    ALLOCATE (ddptr)
    NULLIFY (ddptr%attributes, ddptr%categories, &
      & ddptr%contfills, ddptr%catfills, ddptr%usedvars, ddptr%usedconts)
  END FUNCTION new_datadescription

  FUNCTION new_instanceset(set_type) RESULT (isptr)
    INTEGER :: set_type
    TYPE (instanceset), POINTER :: isptr

    ALLOCATE (isptr)
    NULLIFY (isptr%contvars, isptr%catvars, isptr%ignorevars, &
      & isptr%missing_data, isptr%contvar_perms, isptr%contvar_ranks, &
      & isptr%contvar_5_percent, isptr%contvar_95_percent, &
      & isptr%class_populations, isptr%tagged_instances, isptr%dd, &
      & isptr%constant_attr, isptr%missing_count, isptr%estimated_class, &
      & isptr%votes, isptr%prox, isptr%prox_index, isptr%outlying, &
      & isptr%leaf_id, isptr%scaling, isptr%par_stripes, isptr%par_offsets)
    isptr%set_type = set_type
  END FUNCTION new_instanceset

  SUBROUTINE free_datadescription(ddptr)
    TYPE (datadescription), POINTER :: ddptr
    IF (ASSOCIATED(ddptr)) THEN
      IF (ASSOCIATED(ddptr%categories)) DEALLOCATE (ddptr%categories)
      IF (ASSOCIATED(ddptr%attributes)) DEALLOCATE (ddptr%attributes)
      IF (ASSOCIATED(ddptr%contfills)) DEALLOCATE (ddptr%contfills)
      IF (ASSOCIATED(ddptr%catfills)) DEALLOCATE (ddptr%catfills)
      IF (ASSOCIATED(ddptr%usedvars)) DEALLOCATE (ddptr%usedvars)
      IF (ASSOCIATED(ddptr%usedconts)) DEALLOCATE (ddptr%usedconts)
      DEALLOCATE (ddptr)
    END IF
  END SUBROUTINE free_datadescription

  SUBROUTINE free_instanceset(isptr)
    TYPE (instanceset), POINTER :: isptr
    IF (ASSOCIATED(isptr)) THEN
      IF (ASSOCIATED(isptr%ignorevars)) DEALLOCATE (isptr%ignorevars)
      IF (ASSOCIATED(isptr%catvars)) DEALLOCATE (isptr%catvars)
      IF (ASSOCIATED(isptr%contvars)) DEALLOCATE (isptr%contvars)
      IF (ASSOCIATED(isptr%contvar_perms)) DEALLOCATE (isptr%contvar_perms)
      IF (ASSOCIATED(isptr%contvar_ranks)) DEALLOCATE (isptr%contvar_ranks)
      IF (ASSOCIATED(isptr%missing_count)) DEALLOCATE (isptr%missing_count)
      IF (ASSOCIATED(isptr%estimated_class)) DEALLOCATE (isptr%estimated_class)
      IF (ASSOCIATED(isptr%votes)) DEALLOCATE (isptr%votes)
      IF (ASSOCIATED(isptr%prox)) DEALLOCATE (isptr%prox)
      IF (ASSOCIATED(isptr%prox_index)) DEALLOCATE (isptr%prox_index)
      IF (ASSOCIATED(isptr%outlying)) DEALLOCATE (isptr%outlying)
      IF (ASSOCIATED(isptr%leaf_id)) DEALLOCATE (isptr%leaf_id)
      IF (ASSOCIATED(isptr%scaling)) DEALLOCATE (isptr%scaling)
      IF (ASSOCIATED(isptr%par_stripes)) DEALLOCATE (isptr%par_stripes)
      IF (ASSOCIATED(isptr%par_offsets)) DEALLOCATE (isptr%par_offsets)
      IF (ASSOCIATED(isptr%tagged_instances)) &
        & DEALLOCATE (isptr%tagged_instances)
      IF (ASSOCIATED(isptr%contvar_5_percent)) &
        & DEALLOCATE (isptr%contvar_5_percent)
      IF (ASSOCIATED(isptr%contvar_95_percent)) &
        & DEALLOCATE (isptr%contvar_95_percent)
      IF (ASSOCIATED(isptr%class_populations)) &
        & DEALLOCATE (isptr%class_populations)
      IF (ASSOCIATED(isptr%missing_data)) &
        & CALL free_bitvector(isptr%missing_data)
      IF (ASSOCIATED(isptr%constant_attr)) &
        & CALL free_bitvector(isptr%constant_attr)
      DEALLOCATE (isptr)
    END IF
  END SUBROUTINE free_instanceset

  FUNCTION skip_comma_and_spaces(string, pos) RESULT (ok)
    CHARACTER(*) :: string
    INTEGER :: pos, length
    LOGICAL :: ok
    LOGICAL :: comma_ok
    comma_ok = .TRUE.
    ok = .TRUE.
    length = LEN_TRIM(string)
    IF (pos.GT.length) THEN
      ok = .TRUE.
      string = ""
    ELSE
      DO
        ! Only one comma
        IF (string(pos:pos).EQ.",".AND.comma_ok) THEN
          comma_ok = .FALSE.
        ! End on a non-space
        ELSE IF (string(pos:pos).NE." ") THEN
          EXIT
        ELSE IF (pos.EQ.length) THEN
          ok = .FALSE.
          EXIT
        END IF
        pos = pos + 1
      END DO
      string = ADJUSTL(string(pos:))
    END IF
  END FUNCTION skip_comma_and_spaces

  FUNCTION parse_token(line, pos)
    CHARACTER(LEN=linelen) :: line
    INTEGER, INTENT(OUT) :: pos
    CHARACTER(LEN=tokenlen) :: parse_token
    CHARACTER :: delimiter

    ! Quoted?
    IF (line(1:1).EQ."'".OR.line(1:1).EQ.'"') THEN
      delimiter = line(1:1)
      line = line(2:)
      pos = SCAN(line, delimiter)
      parse_token = line(1:pos-1)
      pos = pos + 1
    
    ! Simple token.
    ELSE
      pos = SCAN(line, " ,")
      IF (pos.EQ.0) THEN
        parse_token = line
      ELSE
        parse_token = line(1:pos-1)
      END IF
    END IF
  END FUNCTION parse_token

  FUNCTION parse_line(handle, isptr, instance_no, line, pass, line_no, &
    & skip_attr, count_missing) RESULT (error)
    TYPE (instanceset), POINTER :: isptr
    INTEGER, INTENT(OUT) :: instance_no
    CHARACTER(LEN=*) :: line
    INTEGER :: handle, pass, line_no, skip_attr
    LOGICAL :: count_missing
    INTEGER :: error
    INTEGER :: i, j, pos, attribute_count
    LOGICAL :: eof, err
    CHARACTER(LEN=linelen) :: tmpline
    CHARACTER(LEN=tokenlen) :: cat_name

    tmpline = line
    tmpline = ADJUSTL(tmpline)

    ! Comment or a blank line? Ignore it.
    IF (tmpline(1:1).EQ."%".OR.LEN_TRIM(tmpline).EQ.0) THEN
      error = parse_comment
      RETURN
    END IF

    ! Real data?
    instance_no = instance_no + 1
    IF (pass.EQ.2) THEN
      error = parse_error ! We failed unless we reach the end successfully
      attribute_count = UBOUND(isptr%dd%attributes, 1)
      IF (opts%class_attribute_num.EQ.class_new) THEN
        isptr%catvars(instance_no, UBOUND(isptr%catvars, 2)) = 1
        attribute_count = attribute_count - 1
      END IF
      attribute_loop: DO i = 1, attribute_count
        IF (i.NE.skip_attr) THEN
          DO WHILE (tmpline(1:1).EQ.'&')
            CALL read_line(handle, tmpline, line, line_no, eof, err)
            IF (eof.OR.err) GO TO 7999
          END DO
          cat_name = parse_token(tmpline, pos)
          IF (pos.EQ.0) RETURN
          IF (.NOT.skip_comma_and_spaces(tmpline, pos)) GO TO 7999
        END IF

        ! Missing value?
        IF (i.EQ.skip_attr.OR.cat_name.EQ."?") THEN
          IF (count_missing) THEN
            isptr%missing_count(i) = isptr%missing_count(i) + 1
          END IF
          CALL setbit2(isptr%missing_data, instance_no, i, .TRUE.)

          SELECT CASE (isptr%dd%attributes(i)%cat_count)
            CASE (ignorevar)
              isptr%ignorevars(instance_no, isptr%dd%attributes(i)%mapping) &
                & = cat_name
            CASE (contvar)
              isptr%contvars(instance_no, isptr%dd%attributes(i)%mapping) &
                & = missing_cont
            CASE DEFAULT
              isptr%catvars(instance_no, isptr%dd%attributes(i)%mapping) &
                & = missing_cat
          END SELECT
        ! Real value.
        ELSE
          SELECT CASE (isptr%dd%attributes(i)%cat_count)
            CASE (ignorevar)
              isptr%ignorevars(instance_no, isptr%dd%attributes(i)%mapping) &
                & = cat_name
            CASE (contvar)
              READ (cat_name, *, ERR=7999) &
                & isptr%contvars(instance_no, isptr%dd%attributes(i)%mapping)
            CASE DEFAULT
              ! Recognise the category name
              DO j = isptr%dd%attributes(i)%cat_start, &
                & isptr%dd%attributes(i)%cat_start + &
                & isptr%dd%attributes(i)%cat_count - 1
                IF (isptr%dd%categories(j)%name.EQ.cat_name) THEN
                  isptr%catvars(instance_no, isptr%dd%attributes(i)%mapping) = &
                    & j - isptr%dd%attributes(i)%cat_start + 1
                  CYCLE attribute_loop
                END IF
              END DO
              WRITE(0, af) "Invalid category " // TRIM(cat_name) &
                & // " for attribute " // TRIM(isptr%dd%attributes(i)%name)
              RETURN
          END SELECT
        END IF
      END DO attribute_loop

      ! Success only if we reached the end of the line without incidents
      IF (LEN_TRIM(tmpline(pos:)).LE.0) error = parse_ok
    ELSE
      error = parse_ok
      DO
        pos = LEN_TRIM(tmpline)
        IF (tmpline(pos:pos).EQ."&") THEN
          CALL read_line(handle, tmpline, line, line_no, eof, err)
        ELSE
          EXIT
        END IF
      END DO
    END IF
    7999 CONTINUE
  END FUNCTION parse_line

  SUBROUTINE read_line(handle, line, origline, line_no, eof, err)
    CHARACTER(LEN=*) :: line, origline
    INTEGER :: handle, line_no
    LOGICAL :: eof, err
    INTEGER :: pos

    eof = .FALSE.
    err = .FALSE.
    DO
      IF (par_front) THEN
        line_no = line_no + 1
        READ(handle, "(A)", END=8010, ERR=8030) line
        origline = line
        DO
          pos = SCAN(line, CHAR(9))
          IF (pos.EQ.0) EXIT
          line(pos:pos) = " "
        END DO
        line = ADJUSTL(line)
      END IF

      GOTO 8020
8010  eof = .TRUE.
      line="%"

! Now distribute the clear text to all processes. The broadcast is in
! human readable form, and the comment sign (%) at line(1:1) is eof
! indication

8020  CONTINUE
      ! Comment? Skip it.
      IF (eof.OR.(line(1:1).NE."%".AND.LEN_TRIM(line).NE.0)) THEN
        CALL par_bcast_char(line)
        IF ((.NOT.par_front).AND.(line(1:1).EQ."%")) eof = .TRUE.
        RETURN
      END IF
    END DO
    RETURN

8030 err = .TRUE.
  END SUBROUTINE read_line

  FUNCTION parse_arff(isptr, filename)
    TYPE (instanceset), POINTER :: isptr
    CHARACTER(*) :: filename
    LOGICAL :: parse_arff
    CHARACTER(LEN=linelen) :: origline, line
    INTEGER :: line_no, attribute_count, cat_count, instance_count, cat_start
    INTEGER :: contvar_count, catvar_count, ignorevar_count
    INTEGER :: pos, pos2, pass, i, attr, attr2
    INTEGER :: handle
    CHARACTER(LEN=tokenlen) :: cat_name
    CHARACTER :: delimiter
    LOGICAL :: statement_attribute, statement_ignored
    LOGICAL :: eof, err
    INTEGER :: parse_line_result, constant_attrs
    LOGICAL :: loaddd ! load data description
    INTEGER :: skip_attr
    LOGICAL, ALLOCATABLE :: usedvars(:)

    parse_arff = .FALSE.
    loaddd = .NOT.ASSOCIATED(isptr%dd)
    IF (loaddd) isptr%dd => new_datadescription()

    IF (par_front) THEN
      handle = open_file(filename, .FALSE., err)
      IF (err) GO TO 9000
    END IF

    ! Two passes. In the first pass, get the structure.
    ! Then allocate arrays, and fill them in the second pass.
    pass_loop: DO pass = 1, 2
      line(1:1) = "!"  ! prevent an exceptional situation in parse_line
      IF (pass.EQ.2) THEN
        ! Start of pass 2.
        if (par_front) REWIND(handle)
        IF (loaddd) THEN
          IF (opts%class_attribute_num.EQ.class_new) THEN
            instance_count = instance_count * (1 + opts%new_class_quantity)
            attribute_count = attribute_count + 1
            catvar_count = catvar_count + 1
            cat_count = cat_count + 2
            isptr%tagged = .TRUE.
          END IF
          isptr%dd%contvar_count = contvar_count
          isptr%dd%catvar_count = catvar_count
          ALLOCATE (isptr%dd%attributes(attribute_count),&
            & isptr%dd%categories(cat_count),&
            & usedvars(attribute_count), isptr%missing_count(attribute_count))
          isptr%constant_attr => new_bitvector(attribute_count)
          usedvars = opts%used_default
          isptr%missing_count = 0
          skip_attr = attribute_count + 1 ! i.e. don't skip.
          IF (opts%split_node_size.LT.0) THEN
            opts%split_node_size = -opts%split_node_size * instance_count / 100
          END IF
          IF (opts%split_node_size.GT.instance_count) THEN
            opts%split_node_size = instance_count
          END IF
        ELSE
          ! TODO check compatibility
          i = attribute_count
          attribute_count = UBOUND(isptr%dd%attributes, 1)
          IF (i + 1.EQ.attribute_count) THEN
            skip_attr = opts%class_attribute_num ! skip the class attribute
            catvar_count = catvar_count + 1
            isptr%tagged = .FALSE.
          ELSE IF (i.EQ.attribute_count) THEN
            skip_attr = attribute_count + 1 ! i.e. don't skip.
            isptr%tagged = .TRUE.
          ELSE
            GO TO 9070 ! incompatible attribute schema sizes
          END IF
        END IF
        isptr%missing_data => new_bitvector2(instance_count, &
          & attribute_count)
        ALLOCATE (isptr%contvars(instance_count, contvar_count),&
          & isptr%catvars(instance_count, catvar_count),&
          & isptr%ignorevars(instance_count, ignorevar_count))
        IF (loaddd.AND.opts%class_attribute_num.EQ.class_new) THEN
          isptr%dd%attributes(attribute_count)%name = "class"
          isptr%dd%attributes(attribute_count)%cat_count = 2
          isptr%dd%attributes(attribute_count)%cat_start = cat_count - 1
          isptr%dd%attributes(attribute_count)%mapping = catvar_count
          isptr%dd%categories(cat_count - 1)%attribute = attribute_count
          isptr%dd%categories(cat_count)%attribute = attribute_count
          isptr%dd%categories(cat_count - 1)%name = "original"
          isptr%dd%categories(cat_count)%name = "constructed"
          isptr%dd%categories(cat_count - 1)%weight = 1
          isptr%dd%categories(cat_count)%weight = 1
        END IF
      END IF

      line_no = 0
      attribute_count = 0
      cat_count = 0
      instance_count = 0
      contvar_count = 0
      catvar_count = 0
      ignorevar_count = 0
      DO
        CALL read_line(handle, line, origline, line_no, eof, err)
        IF (eof) GO TO 9010
        IF (err) GO TO 9030
        statement_attribute = lcase(line(1:10)).EQ."@attribute"
        statement_ignored = lcase(line(1:8)).EQ."@ignored"

        ! Data block? Read data 'till EOF.
        IF (lcase(line(1:5)).EQ."@data") THEN
          DO
            CALL read_line(handle, line, origline, line_no, eof, err)

            IF (eof) GO TO 1000
            IF (err) GO TO 9030
            parse_line_result = parse_line(handle, isptr, instance_count, &
              & line, pass, line_no, skip_attr, pass.NE.1.AND.loaddd)
            SELECT CASE (parse_line_result)
              CASE (parse_ok)
              CASE (parse_comment)
              CASE (parse_error)
                GO TO 9020
            END SELECT
          END DO
          1000 CYCLE pass_loop

        ! Attribute specification?
        ELSE IF (statement_attribute.OR.statement_ignored) THEN
          attribute_count = attribute_count + 1
          line = ADJUSTL(line(SCAN(line, ' ')+1:))

          ! Get name
          cat_name = parse_token(line, pos)
          IF (pass.EQ.2) THEN
            IF (loaddd) THEN
              isptr%dd%attributes(attribute_count)%name = cat_name
              isptr%dd%attributes(attribute_count)%cat_count = 0
            ELSE
              ! TODO check compatibility
            END IF
          ELSE IF (loaddd.AND. &
            & lcase(TRIM(cat_name)).EQ.TRIM(opts%class_attribute)) THEN
            opts%class_attribute_num = attribute_count
          END IF
          line = ADJUSTL(line(pos+1:))

          ! Catvar?
          IF (line(1:1).EQ."{") THEN
            cat_start = cat_count + 1
            IF (pass.EQ.2) THEN
              IF (loaddd) THEN
                isptr%dd%attributes(attribute_count)%cat_start = cat_start
              ELSE
                ! TODO check compatibility
              END IF
            END IF
            pos = 2
            IF (.NOT.skip_comma_and_spaces(line, pos)) GO TO 9020
            DO WHILE (line(1:1).NE."}")
              DO WHILE (line(1:1).EQ."&")
                CALL read_line(handle, line, origline, line_no, eof, err)
                IF (eof) GO TO 9010
                IF (err) GO TO 9030
              END DO
              IF (line(1:1).EQ."'".OR.line(1:1).EQ.'"') THEN
                delimiter = line(1:1)
                line = line(2:)
                pos = SCAN(line, delimiter)
                cat_name = line(1:pos-1)
                pos = pos + 1
              ELSE
                pos = SCAN(line, " ,}(")
                cat_name = line(1:pos-1)
              END IF
              cat_count = cat_count + 1
              IF (pos.EQ.0) GO TO 9020
              IF (pass.EQ.2) THEN
                IF (loaddd) THEN
                  isptr%dd%categories(cat_count)%name = cat_name
                  isptr%dd%categories(cat_count)%attribute = attribute_count
                  isptr%dd%attributes(attribute_count)%cat_count = &
                    & isptr%dd%attributes(attribute_count)%cat_count + 1
                ELSE
                  ! TODO check compatibility
                END IF
              END IF
              IF (.NOT.skip_comma_and_spaces(line, pos)) GO TO 9020
              IF (line(1:1).EQ.'(') THEN ! category weight
                pos = SCAN(line, ')')
                IF (pos.EQ.0) GO TO 9020
                cat_name = line(2:pos-1)
                pos = pos + 1
                IF (.NOT.skip_comma_and_spaces(line, pos)) GO TO 9020
                IF (pass.EQ.2) THEN
                  IF (loaddd) THEN
                    READ(cat_name, *), isptr%dd%categories(cat_count)%weight
                  ELSE
                    ! TODO check compatibility
                  END IF
                END IF
              ELSE IF (pass.EQ.2) THEN
                IF (loaddd) THEN
                  isptr%dd%categories(cat_count)%weight = 1
                ELSE
                  ! TODO check compatibility
                END IF
              END IF
            END DO

            IF (cat_start.EQ.cat_count) THEN
              ! only 1 category - ignore this attribute
              cat_count = cat_count - 1
              ignorevar_count = ignorevar_count + 1
              IF (pass.EQ.2) THEN
                IF (loaddd) THEN
                  isptr%dd%attributes(attribute_count)%cat_count = ignorevar
                  isptr%dd%attributes(attribute_count)%mapping = ignorevar_count
                ELSE
                  ! TODO check compatibility
                END IF
              END IF
            ELSE
              catvar_count = catvar_count + 1
              IF (pass.EQ.2) THEN
                IF (loaddd) THEN
                  isptr%dd%attributes(attribute_count)%mapping = catvar_count
                ELSE
                  ! TODO check compatibility
                END IF
              END IF
            END IF

          ! Contvar?
          ELSE IF (lcase(line(1:7)).EQ."numeric" &
            & .OR.lcase(line(1:7)).EQ."integer" & ! We'll treat ints as conts
            & .OR.lcase(line(1:4)).EQ."real") THEN
            contvar_count = contvar_count + 1
            IF (pass.EQ.2) THEN
              IF (loaddd) THEN
                isptr%dd%attributes(attribute_count)%mapping = contvar_count
                isptr%dd%attributes(attribute_count)%cat_count = contvar
              ELSE
                ! TODO check compatibility
              END IF
            END IF

          ! Something else (presumably a string)
          ELSE
            ignorevar_count = ignorevar_count + 1
            IF (pass.EQ.2) THEN
              IF (loaddd) THEN
                isptr%dd%attributes(attribute_count)%mapping = ignorevar_count
                isptr%dd%attributes(attribute_count)%cat_count = ignorevar
              ELSE
                ! TODO check compatibility
              END IF
            END IF
          END IF

          ! Ignore the @ignored
          IF (pass.EQ.2) THEN
            IF (loaddd) THEN
              usedvars(attribute_count) = usedvars(attribute_count).AND.&
                & statement_attribute
            ELSE
              ! TODO check compatibility
            END IF
          END IF
        ELSE IF (pass.EQ.1) THEN
          CYCLE
        ELSE IF (lcase(line(1:9)).EQ."@relation") THEN
          pos = 1
          line = ADJUSTL(line(11:))
          isptr%relation_name = parse_token(line, pos)
        ELSE
          GO TO 9020
        END IF
      END DO
    END DO pass_loop

    IF (loaddd) THEN
      ! See which attribute is the classification attribute
      IF (opts%class_attribute_num.EQ.class_named) THEN
        opts%class_attribute_num = get_attribute_by_name(isptr, &
          & opts%class_attribute)
        IF (opts%class_attribute_num.EQ.0) GO TO 9040
      ELSE IF (opts%class_attribute_num.LT.class_named) THEN
        CALL generate_class(isptr, instance_count)
        opts%class_attribute_num = UBOUND(isptr%dd%attributes, 1)
      END IF
      cat_start = isptr%dd%attributes(opts%class_attribute_num)%cat_start
      cat_count = isptr%dd%attributes(opts%class_attribute_num)%cat_count
      IF (cat_count.LT.2) GO TO 9040

      ! Is there weight override?
      IF (LEN_TRIM(opts%class_weights).NE.0) THEN
        READ (opts%class_weights, *) &
          & isptr%dd%categories(cat_start:cat_start+cat_count-1)%weight
      END IF

      ! What is the positive category?
      IF (LEN_TRIM(opts%positive_category).NE.0) THEN
        opts%positive_category_num = &
          & get_category_by_name(isptr, opts%positive_category)
      END IF

      ! See which attribute is the label
      IF (LEN_TRIM(opts%label).NE.0) THEN
        labelattr = get_attribute_by_name(isptr, opts%label)
        IF (labelattr.EQ.0) GO TO 9080
        label => isptr%dd%attributes(labelattr)
      ELSE
        labelattr = 0
      END IF

      ! See which attributes are not used
      IF (LEN_TRIM(opts%unusedvars).NE.0) THEN
        line = TRIM(opts%unusedvars) // ','
        DO WHILE (LEN_TRIM(line).NE.0)
          pos = SCAN(line, ',')
          IF (pos.EQ.0) pos = LEN_TRIM(line) + 1
          cat_name = line(1:pos-1)
          attr = get_attribute_by_name(isptr, cat_name)
          IF (attr.EQ.0) THEN
            pos2 = SCAN(line, '-')
            cat_name = line(1:pos2-1)
            attr = get_attribute_by_name(isptr, cat_name)
            IF (attr.EQ.0) GO TO 9050
            cat_name = line(pos2+1:pos-1)
            attr2 = get_attribute_by_name(isptr, cat_name)
            IF (attr2.EQ.0) GO TO 9050
            usedvars(attr:attr2) = .FALSE.
          ELSE
            usedvars(attr) = .FALSE.
          END IF
          line = line(pos+1:)
        END DO
      END IF
      ! and which are:
      IF (LEN_TRIM(opts%usedvars).NE.0) THEN
        line = TRIM(opts%usedvars) // ','
        DO WHILE (LEN_TRIM(line).NE.0)
          pos = SCAN(line, ',')
          IF (pos.EQ.0) pos = LEN_TRIM(line) + 1
          cat_name = line(1:pos-1)
          attr = get_attribute_by_name(isptr, cat_name)
          IF (attr.EQ.0) THEN
            pos2 = SCAN(line, '-')
            cat_name = line(1:pos2-1)
            attr = get_attribute_by_name(isptr, cat_name)
            IF (attr.EQ.0) GO TO 9050
            cat_name = line(pos2+1:pos-1)
            attr2 = get_attribute_by_name(isptr, cat_name)
            IF (attr2.EQ.0) GO TO 9050
            usedvars(attr:attr2) = .TRUE.
          ELSE
            usedvars(attr) = .TRUE.
          END IF
          line = line(pos+1:)
        END DO
      END IF
      ! Ignorevars can't be used at all
      usedvars = usedvars.AND.(isptr%dd%attributes(:)%cat_count.NE.ignorevar)
      usedvars(opts%class_attribute_num) = .FALSE.

      ! See which attributes are constant (all present, or all missing)
      constant_attrs = 0
      DO i = 1, attribute_count
        IF (isptr%missing_count(i).EQ.0 &
          & .OR. isptr%missing_count(i).EQ.instance_count) THEN
          CALL setbit(isptr%constant_attr, i, .TRUE.)
          IF (isptr%missing_count(i).EQ.instance_count) usedvars(i) = .FALSE.
          constant_attrs = constant_attrs + 1
        END IF
      END DO
      IF (constant_attrs.EQ.attribute_count) opts%fill_passes = 0

      attr = 0
      pos = 0
      DO i = 1, attribute_count
        ! count used vars
        IF (.NOT.usedvars(i)) CYCLE
        attr = attr + 1
        ! count used contvars
        IF (isptr%dd%attributes(i)%cat_count.EQ.contvar) pos = pos + 1
      END DO
      IF (attr.EQ.0) GO TO 9060
      ALLOCATE (isptr%dd%usedvars(attr), isptr%dd%usedconts(pos))
      attr = 0
      pos = 0
      DO i = 1, attribute_count
        IF (usedvars(i)) THEN
          attr = attr + 1
          isptr%dd%usedvars(attr) = i
          IF (isptr%dd%attributes(i)%cat_count.EQ.contvar) THEN
            pos = pos + 1
            isptr%dd%usedconts(pos) = i
          END IF
        END IF
      END DO
      DEALLOCATE (usedvars)
      CALL index_usedconts(isptr%dd)
    END IF
    
    ! Success
    parse_arff = .TRUE.
    GO TO 9999

    ! Error situations
    9000 WRITE(*, *) "Error reading ", TRIM(filename), ": cannot open file"
    GO TO 9990

    9010 WRITE(*, *) "Error reading ", TRIM(filename), ": no data"
    GO TO 9990

    9020 WRITE(*, *) "Error reading ", TRIM(filename), &
      & ": syntax error on line", line_no
    WRITE (*, *) TRIM(origline)
    GO TO 9990

    9030 WRITE(*, *) "Error reading ", TRIM(filename), " on line", line_no
    GO TO 9990

    9040 WRITE(*, *) "Error: Invalid argument to -c"
    GO TO 9990

    9050 WRITE(*, *) "Error: Invalid argument to -u: ", TRIM(cat_name)
    GO TO 9990

    9060 WRITE(*, *) "Error: Must use at least one attribute"
    GO TO 9990

    9070 WRITE(*, *) "Error: The attribute schema are incompatible"
    GO TO 9990

    9080 WRITE(*, *) "Error: Invalid argument to -ri: ", TRIM(cat_name)
    GO TO 9990

    9990 CALL free_instanceset(isptr)

    ! Finally, close the file
    9999 IF (par_front) CALL close_file(handle)

  END FUNCTION parse_arff

  SUBROUTINE index_usedconts(datadesc)
    TYPE (datadescription), POINTER :: datadesc
    INTEGER :: i

    DO i = 1, UBOUND(datadesc%usedconts, 1)
      datadesc%attributes(datadesc%usedconts(i))%usedcont = i
    END DO
  END SUBROUTINE index_usedconts

  SUBROUTINE generate_class(isptr, original_count)
    TYPE (instanceset), POINTER :: isptr
    INTEGER :: i, j
    INTEGER :: original_count, instance_count, attribute_count

    instance_count = UBOUND(isptr%catvars, 1)
    attribute_count = UBOUND(isptr%catvars, 2)
    IF (opts%class_attribute_num.EQ.class_new) THEN
      ! Generate the second class
      opts%class_attribute_num = UBOUND(isptr%catvars, 2)
      DO i = original_count + 1, instance_count
        DO j = 1, UBOUND(isptr%catvars, 2)
          IF (j.EQ.attribute_count) THEN ! generated attribute is always last
            isptr%catvars(i, j) = 2 ! "constructed"
          ELSE
            isptr%catvars(i, j) = isptr%catvars(rnd_int(original_count), j)
          END IF
        END DO
        DO j = 1, UBOUND(isptr%contvars, 2)
          isptr%contvars(i, j) = isptr%contvars(rnd_int(original_count), j)
        END DO
        ! FIXME This should probably be cut, because it just slows things down.
        ! Also, perhaps filter out the unused vars?
        DO j = 1, UBOUND(isptr%ignorevars, 2)
          isptr%ignorevars(i, j) = isptr%ignorevars(rnd_int(instance_count), j)
        END DO
      END DO
    ELSE IF (opts%class_attribute_num.EQ.class_last) THEN
      opts%class_attribute_num = attribute_count
    END IF
  END SUBROUTINE generate_class

  PURE FUNCTION escape_token(string)
    CHARACTER(*), INTENT(IN) :: string
    CHARACTER(LEN=LEN(string)+2) :: escape_token

    ! Do we need single quotes?
    IF (INDEX(TRIM(string), '"').NE.0) THEN
      escape_token = "'" // TRIM(string) // "'"

    ! Do we need quotes at all?
    ELSE IF (SCAN(TRIM(string), " '").NE.0) THEN
      escape_token = '"' // TRIM(string) // '"'

    ! A simple token
    ELSE
      escape_token = TRIM(string)
    END IF
  END FUNCTION escape_token
  
  SUBROUTINE print_datadescription(dd, relname, handle)
    TYPE (datadescription), POINTER :: dd
    CHARACTER(LEN=*) :: relname
    INTEGER :: handle, i, j
    CHARACTER(LEN=1024) line, val

    ! Relation name
    WRITE(handle, af) "@relation " // TRIM(relname)
    WRITE(handle, *)

    ! Attribute descriptors
    DO i = 1, UBOUND(dd%attributes, 1)
      SELECT CASE (dd%attributes(i)%cat_count)
        CASE (ignorevar)
          line = "string"
        CASE (contvar)
          line = "numeric"
        CASE DEFAULT
          line = ""
          DO j = dd%attributes(i)%cat_start, &
            & dd%attributes(i)%cat_start &
            & + dd%attributes(i)%cat_count - 1
            line = TRIM(line) // ", " // &
              & escape_token(dd%categories(j)%name)
            IF (dd%categories(j)%weight.NE.1) THEN
              CALL squeeze_num(val, dd%categories(j)%weight)
              line = TRIM(line) // " (" // TRIM(val) // ")"
            END IF
          END DO
          line = "{ " // TRIM(line(3:)) // " }"
      END SELECT
      WRITE(handle, af) "@attribute " &
        & // TRIM(escape_token(dd%attributes(i)%name)) // " " // TRIM(line)
    END DO
    WRITE(handle, *)
  END SUBROUTINE print_datadescription

  SUBROUTINE print_instances(isptr, handle)
    TYPE (instanceset), POINTER :: isptr
    INTEGER :: handle

    CALL print_datadescription(isptr%dd, isptr%relation_name, handle)
    WRITE(handle, af) "@data"
    CALL print_data(isptr, handle)
  END SUBROUTINE print_instances

  SUBROUTINE print_data(isptr, handle)
    TYPE (instanceset), POINTER :: isptr
    INTEGER :: handle
    INTEGER :: i, j
    CHARACTER(LEN=1024) line, val
    INTEGER :: cat_tmp

    DO i = 1, UBOUND(isptr%catvars, 1)
      line = ""
      DO j = 1, UBOUND(isptr%dd%attributes, 1)
        IF (getbit2(isptr%missing_data, i, j)) THEN
          line = TRIM(line) // ", ?"
          CYCLE
        END IF
        SELECT CASE (isptr%dd%attributes(j)%cat_count)
          CASE (ignorevar)
            val = escape_token(isptr%ignorevars(i, &
              & isptr%dd%attributes(j)%mapping))
          CASE (contvar)
            CALL squeeze_num(val, &
              & isptr%contvars(i, isptr%dd%attributes(j)%mapping))
          CASE DEFAULT
            IF (j.EQ.opts%class_attribute_num) THEN
              cat_tmp = isptr%classes(i)
            ELSE
              cat_tmp = isptr%catvars(i, isptr%dd%attributes(j)%mapping)
            END IF
            val = escape_token(isptr%dd%categories( &
              & isptr%dd%attributes(j)%cat_start+cat_tmp-1)%name)
        END SELECT
        line = TRIM(line) // ", " // TRIM(val)
      END DO
      WRITE(handle, af) TRIM(line(3:))
    END DO
  END SUBROUTINE print_data

  SUBROUTINE print_arff(filename, isptr, isptr2)
    TYPE (instanceset), POINTER :: isptr
    TYPE (instanceset), POINTER, OPTIONAL :: isptr2
    CHARACTER(LEN=*) :: filename
    INTEGER :: handle
    LOGICAL :: err
    
    handle = open_file(filename, .TRUE., err)
    IF (err) RETURN

    CALL print_instances(isptr, handle)
    IF (PRESENT(isptr2)) THEN
      WRITE(handle, af) "%-------"
      CALL print_data(isptr2, handle)
    END IF

    CALL close_file(handle)
  END SUBROUTINE print_arff

  SUBROUTINE calculate_rough_fills(isptr)
    TYPE (instanceset), POINTER :: isptr
    INTEGER :: i, a, attr

    ! Allocate arrays
    ALLOCATE (isptr%dd%contfills(UBOUND(isptr%contvars, 2)), &
      & isptr%dd%catfills(UBOUND(isptr%catvars, 2)))

    ! Calculate rough fills
    DO i = 1, UBOUND(isptr%dd%usedvars, 1)
      a = isptr%dd%usedvars(i)
      IF (getbit(isptr%constant_attr, a)) CYCLE
      attr = isptr%dd%attributes(a)%mapping
      SELECT CASE (isptr%dd%attributes(a)%cat_count)
        CASE (ignorevar)
        CASE (contvar)
          ! For contvars, median value
          isptr%dd%contfills(attr) = median(isptr%contvars(:,attr))
        CASE DEFAULT
          ! For catvars, the most frequent category
          isptr%dd%catfills(attr) = &
            & most_frequent(isptr%catvars(:,attr), &
            & isptr%dd%attributes(a)%cat_count)
      END SELECT
    END DO
  END SUBROUTINE calculate_rough_fills

  SUBROUTINE fill_missing_rough(isptr)
    TYPE (instanceset), POINTER :: isptr
    INTEGER :: i, j, a, attr, num_instances
    LOGICAL :: ismissing ! C function

    num_instances = UBOUND(isptr%catvars, 1)
    DO i = 1, UBOUND(isptr%dd%usedvars, 1)
      a = isptr%dd%usedvars(i)
      IF (getbit(isptr%constant_attr, a)) CYCLE
      attr = isptr%dd%attributes(a)%mapping
      SELECT CASE (isptr%dd%attributes(a)%cat_count)
        CASE (ignorevar)
        CASE (contvar)
          DO j = 1, num_instances
            IF (ismissing(isptr%contvars(j, attr))) THEN
              isptr%contvars(j, attr) = isptr%dd%contfills(attr)
            END IF
          END DO
        CASE DEFAULT
          DO j = 1, num_instances
            IF ((isptr%catvars(j, attr).EQ.missing_cat)) THEN
              isptr%catvars(j, attr) = isptr%dd%catfills(attr)
            END IF
          END DO
      END SELECT
    END DO
  END SUBROUTINE fill_missing_rough

  SUBROUTINE fill_missing_by_prox(isptr)
    TYPE (instanceset), POINTER :: isptr
    INTEGER :: i, j, k, l, a, cat, attr, num_instances, cat_count
    REAL :: weighted_sum, total_prox
    INTEGER, POINTER :: votes(:)
    ! FIXME

PRINT*,"GTDEBUG ",UBOUND(isptr%prox_index), num_instances, opts%num_prox
    num_instances = UBOUND(isptr%catvars, 1)
    DO i = 1, UBOUND(isptr%dd%usedvars, 1)
      a = isptr%dd%usedvars(i)
      IF (getbit(isptr%constant_attr, a)) CYCLE
      attr = isptr%dd%attributes(a)%mapping
      cat_count = isptr%dd%attributes(i)%cat_count
      SELECT CASE (cat_count)
        CASE (ignorevar)
        CASE (contvar)
          ! Continuous variable:
          ! Make a proximity-weighted average of non-missing values
          DO j = isptr%lower, isptr%upper
            IF (getbit2(isptr%missing_data, j, attr)) THEN
              weighted_sum = 0
              total_prox = 0
              DO k = 1, opts%num_prox
                l = isptr%prox_index(j, k)
                IF (.NOT.getbit2(isptr%missing_data, l, attr)) THEN
PRINT*, j, k, isptr%prox(j, k)
                  weighted_sum = weighted_sum + isptr%prox(j, k) &
                    & * isptr%contvars(l, attr)
                  total_prox = total_prox + isptr%prox(j, k)
                END IF
              END DO
              IF (total_prox.GT.0) THEN
                isptr%contvars(j, attr) = weighted_sum / total_prox
              END IF
            END IF
          END DO
        CASE DEFAULT
          ! Categorical variable:
          ! Sum proximities by category; pick highest proximity total
          ALLOCATE (votes(cat_count))
          DO j = isptr%lower, isptr%upper
            votes = 0
            IF (getbit2(isptr%missing_data, j, attr)) THEN
              DO k = 1, opts%num_prox
                l = isptr%prox_index(j, k)
                IF (.NOT.getbit2(isptr%missing_data, l, attr)) THEN
                  cat = isptr%catvars(l, attr)
                  votes(cat) = votes(cat) + isptr%prox(j, k)
                END IF
              END DO
              l = 1
              DO k = 2, cat_count
                IF (votes(k).GE.votes(l)) l = k
              END DO
              isptr%catvars(j, attr) = l
            END IF
          END DO
          DEALLOCATE (votes)
      END SELECT
    END DO
  END SUBROUTINE fill_missing_by_prox

  ! Called after filling missing values, so no ?'s here
  SUBROUTINE sort_and_rank(isptr, skip_constants)
    TYPE (instanceset), POINTER :: isptr
    LOGICAL :: skip_constants
    INTEGER :: i, j, attr, instance_count, contvar_count
    INTEGER :: percent5, percent95
    INTEGER :: last_idx, idx, temp_idx

    instance_count = UBOUND(isptr%catvars, 1)
    contvar_count = UBOUND(isptr%contvars, 2)
    percent5 = NINT((instance_count - 1)*0.05)+1
    percent95 = NINT((instance_count - 1)*0.95)+1

    ! Allocate arrays, but only the first time
    IF (.NOT.ASSOCIATED(isptr%contvar_ranks)) THEN
      ALLOCATE (isptr%contvar_perms(instance_count, contvar_count), &
        & isptr%contvar_ranks(instance_count, contvar_count))
    END IF
    IF (opts%num_prot .GT. 0) THEN ! Only needed for the prototypes
      IF (ASSOCIATED(isptr%contvar_5_percent)) &
        & DEALLOCATE (isptr%contvar_5_percent)
      IF (ASSOCIATED(isptr%contvar_95_percent)) &
        & DEALLOCATE (isptr%contvar_95_percent)
      ALLOCATE (isptr%contvar_5_percent(contvar_count), &
        & isptr%contvar_95_percent(contvar_count))
    END IF

    DO i = 1, UBOUND(isptr%dd%usedvars, 1)
      j = isptr%dd%usedvars(i)

      ! sort constant attributes only once
      IF (skip_constants .AND. getbit(isptr%constant_attr, j)) CYCLE

      attr = isptr%dd%attributes(j)%mapping
      IF (isptr%dd%attributes(j)%cat_count .EQ. contvar) THEN

        ! Create a permutation vector
        CALL qsort_real(isptr%contvars(:,attr), isptr%contvar_perms(:,attr))

        IF (opts%num_prot .GT. 0) THEN
          ! Get the 5% and 95% values
          isptr%contvar_5_percent(attr) = &
            & isptr%contvars(isptr%contvar_perms(percent5, attr), attr)
          isptr%contvar_95_percent(attr) = &
            & isptr%contvars(isptr%contvar_perms(percent95, attr), attr)
        END IF

        ! Rank the instances
        idx = isptr%contvar_perms(1, attr)
        IF (idx.GT.0) THEN ! only if not all values are missing
          isptr%contvar_ranks(idx, attr) = 1
          DO j = 2, instance_count
            temp_idx = isptr%contvar_perms(j, attr)
            ! FIXME this should be idx.EQ.0 EXIT ?
            IF (temp_idx.EQ.0) CYCLE ! Skip missing values altogether
            last_idx = idx
            idx = temp_idx
            IF (isptr%contvars(idx, attr) &
              & .NE. isptr%contvars(last_idx, attr)) THEN
              isptr%contvar_ranks(idx, attr) = &
                & isptr%contvar_ranks(last_idx, attr) + 1
            ELSE
              isptr%contvar_ranks(idx, attr) = &
                & isptr%contvar_ranks(last_idx, attr)
            END IF
          END DO
        ELSE ! if all values are missing, all are rank 1
          DO j = 1, instance_count
            isptr%contvar_ranks(j, attr) = 1
          END DO
        END IF
      END IF
    END DO
  END SUBROUTINE sort_and_rank

  SUBROUTINE count_classes(isptr)
    TYPE (instanceset), POINTER :: isptr
    INTEGER :: i, c, instance_count, class_count, class_catvar
    instance_count = UBOUND(isptr%catvars, 1)
    class_count = isptr%dd%attributes(opts%class_attribute_num)%cat_count
    class_catvar = isptr%dd%attributes(opts%class_attribute_num)%mapping
    ALLOCATE (isptr%class_populations(0:class_count))
    isptr%class_populations = 0
    DO i = 1, instance_count
      c = isptr%catvars(i, class_catvar)
      isptr%class_populations(c) = isptr%class_populations(c) + 1
    END DO
    ALLOCATE (isptr%tagged_instances(instance_count - &
      & isptr%class_populations(0)))
    c = 0
    DO i = 1, instance_count
      IF (isptr%catvars(i, class_catvar).NE.missing_cat) THEN
        c = c + 1
        isptr%tagged_instances(c) = i
      END IF
    END DO
  END SUBROUTINE count_classes

  FUNCTION get_category_by_name(isptr, str)
    INTEGER :: get_category_by_name
    TYPE (instanceset), POINTER :: isptr
    CHARACTER(LEN=*), INTENT(IN) :: str
    TYPE (attribute), POINTER :: class
    INTEGER :: attribute_count, i

    class => isptr%dd%attributes(opts%class_attribute_num)
    DO i = 1, class%cat_count
      IF (TRIM(isptr%dd%categories(i + class%cat_start - 1)%name) &
          & .EQ.TRIM(ADJUSTL(str))) THEN
        get_category_by_name = i
        RETURN
      END IF
    END DO
    
    ! The name was not recognised, see if it is numeric
    READ(str, *, ERR=9200) get_category_by_name
    IF ((get_category_by_name.GE.1).OR. &
      & (get_category_by_name.LE.class%cat_count)) RETURN

    ! Error
    9200 get_category_by_name = 0
  END FUNCTION get_category_by_name

  FUNCTION get_attribute_by_name(isptr, str)
    INTEGER :: get_attribute_by_name
    TYPE (instanceset), POINTER :: isptr
    CHARACTER(LEN=*), INTENT(IN) :: str
    INTEGER :: attribute_count, i

    attribute_count = UBOUND(isptr%dd%attributes, 1)
    DO i = 1, attribute_count
      IF (TRIM(isptr%dd%attributes(i)%name).EQ.TRIM(ADJUSTL(str))) THEN
        get_attribute_by_name = i
        RETURN
      END IF
    END DO
    
    ! The name was not recognised, see if it is numeric
    READ(str, *, ERR=9210) get_attribute_by_name
    IF ((get_attribute_by_name.GE.1).AND. &
      & (get_attribute_by_name.LE.attribute_count)) RETURN

    ! Error
    9210 get_attribute_by_name = 0
  END FUNCTION get_attribute_by_name

  SUBROUTINE get_num_split_variables(dd)
    TYPE (datadescription), POINTER :: dd
    INTEGER :: num_used

    num_used = UBOUND(dd%usedvars, 1)
    IF (opts%split_variables.LT.0) THEN
      opts%split_variables = NINT(SQRT(REAL(num_used)))
    ELSE IF (opts%split_variables.GT.num_used) THEN
      opts%split_variables = num_used
    END IF
  END SUBROUTINE

  FUNCTION label_for(isptr, i)
    TYPE (instanceset), POINTER :: isptr
    INTEGER :: i, attr
    CHARACTER(LEN=10) :: label_for

    IF (labelattr.NE.0) THEN
      IF (getbit2(isptr%missing_data, i, labelattr)) THEN
        label_for = "?"
      ELSE
        SELECT CASE (label%cat_count)
          CASE (ignorevar)
            WRITE (label_for, "(A10)"), isptr%ignorevars(i, label%mapping)
          CASE (contvar)
            WRITE (label_for, "(F10.2)"), isptr%contvars(i, label%mapping)
          CASE DEFAULT
            WRITE (label_for, "(A10)"), isptr%dd%categories( &
              & label%cat_start + &
              & isptr%catvars(i, label%mapping) - 1)%name
        END SELECT
      END IF
    ELSE
      WRITE(label_for, "(I10)"), i
    END IF
  END FUNCTION label_for

  SUBROUTINE print_votes(isptr, filename)
    TYPE (instanceset), POINTER :: isptr
    CHARACTER(LEN=*) :: filename
    INTEGER :: i
    CHARACTER(LEN=20) :: fmtstr
    CHARACTER(LEN=10) :: labelstr
    INTEGER :: handle
    LOGICAL :: err

    handle = open_file(filename, .TRUE., err)
    IF (err) RETURN

    WRITE(fmtstr, *) UBOUND(isptr%votes, 2)
    fmtstr = "(A10, 1X, " // TRIM(ADJUSTL(fmtstr)) // "F10.4)"

    DO i = 1, UBOUND(isptr%catvars, 1)
      labelstr = label_for(isptr, i)
      WRITE(handle, fmtstr), labelstr, isptr%votes(i, :)
    END DO

    CALL close_file(handle)
  END SUBROUTINE

  SUBROUTINE process_confusion_matrix(isptr, filename)
    TYPE (instanceset), POINTER :: isptr
    CHARACTER(LEN=*) :: filename
    INTEGER, POINTER :: confusion(:, :)
    INTEGER :: instance_count, i
    TYPE (attribute), POINTER :: class
    INTEGER, POINTER :: counter
    CHARACTER(LEN=20) :: fmtstr
    CHARACTER(LEN=5) :: namestr
    INTEGER :: handle
    INTEGER :: tp, fp, tn, fn
    LOGICAL :: err

    instance_count = UBOUND(isptr%catvars, 1)
    class => isptr%dd%attributes(opts%class_attribute_num)
    ALLOCATE (confusion(0:class%cat_count, 0:class%cat_count))
    confusion = 0

    DO i = 1, instance_count
      counter => confusion(isptr%catvars(i, class%mapping), &
        & isptr%estimated_class(i))
      counter = counter + 1
    END DO

    IF (opts%positive_category_num.NE.0) THEN
      tp = confusion(opts%positive_category_num, opts%positive_category_num)
      fp = SUM(confusion(1:, opts%positive_category_num)) - tp
      tn = SUM((/(confusion(i, i), i = 1, class%cat_count)/)) - tp
      fn = SUM(confusion(1:, 1:)) - (fp + tp + tn)
      IF (isptr%set_type.EQ.trainset_type) fmtstr = "Trainset"
      IF (isptr%set_type.EQ.testset_type) fmtstr = " Testset"
      WRITE(6, "(A8, 1X, A10, F8.4, A9, F8.4, A14, F8.4)") fmtstr, &
        & "Precision:", REAL(tp) / (tp + fp), &
        & "; Recall:", REAL(tp) / (tp + fn), &
        & "; Specificity:", REAL(tn) / (tn + fp)
    END IF

    IF (LEN_TRIM(filename).GT.0) THEN

      handle = open_file(filename, .TRUE., err)
      IF (.NOT.err) THEN

        WRITE(fmtstr, *) class%cat_count + 1
        fmtstr = "(A6, " // TRIM(ADJUSTL(fmtstr)) // "(1X, A5))"

        WRITE(handle, fmtstr) "Tag\Cl", "NotCl", &
          & (ADJUSTR(isptr%dd%categories(i)%name(1:5)), &
          & i = class%cat_start, class%cat_start + class%cat_count - 1)

        WRITE(fmtstr, *) class%cat_count + 1
        fmtstr = "(1X, A5, " // TRIM(ADJUSTL(fmtstr)) // "I6)"

        DO i = 0, class%cat_count
          IF (i.EQ.0) THEN
            namestr = "NoTag"
          ELSE
            namestr = isptr%dd%categories(class%cat_start + i - 1)%name
          END IF
          WRITE(handle, fmtstr), namestr, confusion(i, :)
        END DO
      END IF
    END IF
    DEALLOCATE (confusion)

    CALL close_file(handle)
  END SUBROUTINE process_confusion_matrix
END MODULE instancesets
