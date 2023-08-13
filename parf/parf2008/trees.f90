MODULE trees
! This is a binary tree implementation. If both %yes and %no are associated,
! then it is a branch; %attribute is the attribute to branch on. Otherwise,
! neither %yes nor %no are associated, and the structure represents a leaf
! node, where the %attribute gives the class. If the attribute is a
! categorical one, %catvar_classes is associated into a bitvector. If the
! attribute is continuous, %contvar_max gives the cutoff point - the highest
! value to be accepted into %yes, and excluded from %no.

  USE bitvectors
  USE instancesets
  USE bootstraps
  USE utilities
  IMPLICIT NONE

  TYPE node
    TYPE (node), POINTER :: yes, no
    INTEGER :: attribute
    REAL :: total_weight
    INTEGER :: num_instances
    REAL :: contvar_max
    TYPE (bitvector), POINTER :: catvar_classes
    INTEGER :: id
    REAL :: dgini ! TODO do we need to save this?
  END TYPE node
  TYPE tree
    TYPE (node), POINTER :: p
    REAL, POINTER :: dgini(:) ! average dgini for this tree
    INTEGER, POINTER :: leaf_pop(:) ! population of each leaf
    INTEGER, POINTER :: leaf_bounds(:) ! indices in leaf_index for each leaf
    INTEGER, POINTER :: leaf_index(:) ! instance leaf block permutation table
    INTEGER :: terminal_bound
    TYPE (bitvector), POINTER :: oob
  END TYPE tree
CONTAINS
  FUNCTION new_node() RESULT (nodeptr)
    TYPE (node), POINTER :: nodeptr
    ALLOCATE (nodeptr)
    NULLIFY (nodeptr%yes, nodeptr%no, nodeptr%catvar_classes)
  END FUNCTION new_node

  RECURSIVE SUBROUTINE free_node(nodeptr)
    TYPE (node), POINTER :: nodeptr
    IF (ASSOCIATED(nodeptr%yes)) THEN
      CALL free_node(nodeptr%yes)
      CALL free_node(nodeptr%no)
      CALL free_bitvector(nodeptr%catvar_classes)
    END IF
    DEALLOCATE (nodeptr)
  END SUBROUTINE free_node

  SUBROUTINE free_tree(thetree)
    TYPE (tree) :: thetree
    IF (ASSOCIATED(thetree%p)) CALL free_node(thetree%p)
    IF (ASSOCIATED(thetree%dgini)) DEALLOCATE (thetree%dgini)
    IF (ASSOCIATED(thetree%leaf_pop)) DEALLOCATE (thetree%leaf_pop)
    IF (ASSOCIATED(thetree%leaf_index)) DEALLOCATE (thetree%leaf_index)
    IF (ASSOCIATED(thetree%leaf_bounds)) DEALLOCATE (thetree%leaf_bounds)
    IF (ASSOCIATED(thetree%oob)) CALL free_bitvector(thetree%oob)
  END SUBROUTINE free_tree

  SUBROUTINE build_tree(isptr, bs, treeptr)
    TYPE (instanceset), POINTER :: isptr
    TYPE (bootstrap), POINTER :: bs
    TYPE (tree) :: treeptr
    INTEGER :: nonterminal_bound, i, j, k, a, s
    TYPE (attribute), POINTER :: attr
    TYPE (bitvector_p), POINTER :: used_cats(:)
    INTEGER, POINTER :: sorted_idx(:, :)
    INTEGER :: attribute_count, instance_count
    LOGICAL, POINTER :: absent_from_idx(:)

    NULLIFY (treeptr%dgini, treeptr%leaf_pop, treeptr%leaf_index, &
      & treeptr%leaf_bounds)
    attribute_count = UBOUND(isptr%dd%attributes, 1)
    instance_count = UBOUND(isptr%classes, 1)
    ALLOCATE (used_cats(UBOUND(isptr%catvars, 2)), &
      & sorted_idx(UBOUND(bs%idx, 1), UBOUND(isptr%dd%usedconts, 1)), &
      & absent_from_idx(UBOUND(isptr%classes, 1)))

    DO i = 1, attribute_count
      attr => isptr%dd%attributes(i)
      IF (attr%cat_count.GE.catvar) THEN
        used_cats(attr%mapping)%p => new_bitvector(attr%cat_count)
      END IF
    END DO

    absent_from_idx = .TRUE.
    absent_from_idx(bs%idx) = .FALSE.

    DO i = 1, UBOUND(isptr%dd%usedconts, 1)
      a = isptr%dd%attributes(isptr%dd%usedconts(i))%mapping
      k = 1
      DO j = 1, instance_count
        s = isptr%contvar_perms(j, a)
        IF (s.EQ.0) EXIT ! missing data at the end
        IF (absent_from_idx(s)) CYCLE
        sorted_idx(k, i) = s
        k = k + 1
      END DO
      sorted_idx(k:, i) = 0
    END DO
    treeptr%terminal_bound = 0
    nonterminal_bound = 1
    treeptr%p => build_tree_internal(isptr, bs, bs%idx, &
      & treeptr%terminal_bound, nonterminal_bound, used_cats, 1, sorted_idx, &
      & absent_from_idx)
    DO i = 1, UBOUND(used_cats, 1)
      CALL free_bitvector(used_cats(i)%p)
    END DO
    DEALLOCATE (used_cats, sorted_idx, absent_from_idx)
  END SUBROUTINE build_tree

  RECURSIVE FUNCTION build_tree_internal(isptr, bs, idx, &
    & terminal_id, nonterminal_id, used_cats, sort_offset, sorted_idx, &
    & absent_from_idx) RESULT (nodeptr)
    TYPE (instanceset), POINTER :: isptr
    TYPE (bootstrap), POINTER :: bs
    INTEGER :: terminal_id, nonterminal_id
    TYPE (bitvector_p), POINTER :: used_cats(:)
    INTEGER :: sort_offset
    INTEGER, POINTER :: sorted_idx(:, :)
    LOGICAL, POINTER :: absent_from_idx(:)
    INTEGER, POINTER :: idx(:), yes_idx(:), no_idx(:), temp_sort(:)
    INTEGER :: subset_count, yes_count, no_count, cat_count
    INTEGER :: i, j, k, s, sort_end
    LOGICAL, POINTER :: yes_set(:)
    TYPE (node), POINTER :: nodeptr
    TYPE (bitvector), POINTER :: used_temp1, used_temp2
    LOGICAL :: is_catvar, is_big_cat
    TYPE (attribute), POINTER :: curr_attr
    LOGICAL :: ismissing ! C function

    nodeptr => new_node()
    subset_count = UBOUND(idx, 1)

    ! Non-terminal?
    IF (find_best_split(nodeptr, isptr, bs, idx, used_cats, sorted_idx, &
      & sort_offset)) THEN
      nonterminal_id = nonterminal_id - 1
      nodeptr%id = nonterminal_id
      ! Split the sets
      ALLOCATE (yes_set(subset_count))
      yes_count = 0
      DO i = 1, subset_count
        yes_set(i) = satisfies(nodeptr, isptr, idx(i))
        IF (yes_set(i)) yes_count = yes_count + 1
      END DO

      ALLOCATE (yes_idx(yes_count), no_idx(subset_count - yes_count))
      yes_count = 0
      no_count = 0
      DO i = 1, subset_count
        IF (yes_set(i)) THEN
          yes_count = yes_count + 1
          yes_idx(yes_count) = idx(i)
        ELSE
          no_count = no_count + 1
          no_idx(no_count) = idx(i)
        END IF
      END DO
      DEALLOCATE (yes_set)

      ! Prettify - take care of which categories have already been decided upon
      curr_attr => isptr%dd%attributes(nodeptr%attribute)
      cat_count = curr_attr%cat_count
      is_catvar = cat_count.GE.catvar
      IF (is_catvar) THEN
        is_big_cat = curr_attr%cat_count.GE.opts%big_catvar_cat_count
        IF (is_big_cat) THEN
          used_temp1 => used_cats(curr_attr%mapping)%p
          used_cats(curr_attr%mapping)%p => &
            & not_bitvector(nodeptr%catvar_classes, cat_count)
        ELSE
          used_temp1 => clone_bitvector(used_cats(curr_attr%mapping)%p)
          CALL or_bitvector(used_cats(curr_attr%mapping)%p, &
            & nodeptr%catvar_classes)
        END IF
      ELSE
        is_big_cat = .FALSE.
      END IF
      
      ! rearrange the sorted index
      absent_from_idx(no_idx) = .TRUE.
      sort_end = sort_offset + UBOUND(idx, 1) - 1
      ALLOCATE (temp_sort(1:(sort_end - sort_offset + 1)))
      DO i = 1, UBOUND(sorted_idx, 2)
        k = 0
        ! first the yesses
        DO j = sort_offset, sort_end
          s = sorted_idx(j, i)
          IF (s.EQ.0) THEN ! missing values?
            temp_sort(k+1:yes_count-1) = 0
            EXIT ! we're done here
          END IF
          IF (.NOT.absent_from_idx(s)) THEN
            k = k + 1
            temp_sort(k) = s
          END IF
        END DO
        ! then the nos
        DO j = sort_offset, sort_end
          s = sorted_idx(j, i)
          IF (s.EQ.0) THEN ! missing values?
            temp_sort(k+1:) = 0
            EXIT ! we're done here
          END IF
          IF (absent_from_idx(s)) THEN
            k = k + 1
            temp_sort(k) = s
          END IF
        END DO
        sorted_idx(sort_offset:sort_end, i) = temp_sort
      END DO
      DEALLOCATE (temp_sort)

      ! "Yes" subtree
      nodeptr%yes => build_tree_internal(isptr, bs, yes_idx, &
        & terminal_id, nonterminal_id, used_cats, sort_offset, &
        & sorted_idx, absent_from_idx)
      IF (is_big_cat) THEN
        used_temp2 => not_bitvector(used_temp1, cat_count)
        CALL eor_bitvector(used_cats(curr_attr%mapping)%p, used_temp2)
        CALL free_bitvector(used_temp2)
      END IF
      ! "No" subtree
      absent_from_idx(yes_idx) = .TRUE.
      absent_from_idx(no_idx) = .FALSE.
      nodeptr%no => build_tree_internal(isptr, bs, no_idx, &
        & terminal_id, nonterminal_id, used_cats, sort_offset + yes_count, &
        & sorted_idx, absent_from_idx)
      IF (is_catvar) THEN
        CALL free_bitvector(used_cats(curr_attr%mapping)%p)
        used_cats(curr_attr%mapping)%p => used_temp1
      END IF
      DEALLOCATE (yes_idx, no_idx)
    ELSE
      terminal_id = terminal_id + 1
      nodeptr%id = terminal_id
      nodeptr%attribute = most_frequent_class(isptr, idx)
    END IF
  END FUNCTION build_tree_internal

  FUNCTION satisfies(nodeptr, isptr, instance)
    LOGICAL :: satisfies
    TYPE (node), POINTER :: nodeptr
    TYPE (instanceset), POINTER :: isptr
    INTEGER :: instance
    INTEGER :: attr

    attr = isptr%dd%attributes(nodeptr%attribute)%mapping
    IF (isptr%dd%attributes(nodeptr%attribute)%cat_count.EQ.contvar) THEN
      satisfies = isptr%contvars(instance, attr).LE.nodeptr%contvar_max
    ELSE
      satisfies = getbit(nodeptr%catvar_classes, isptr%catvars(instance, attr))
    END IF
  END FUNCTION satisfies

  FUNCTION find_best_split(nodeptr, isptr, bs, idx, used_cats, sorted_idx, &
    & sort_offset) RESULT (nonterminal)
    LOGICAL :: nonterminal
    TYPE (node), POINTER :: nodeptr
    TYPE (instanceset), POINTER :: isptr
    TYPE (bootstrap) :: bs
    INTEGER, POINTER :: idx(:)
    INTEGER, POINTER :: sorted_idx(:, :)
    INTEGER :: sort_offset
    TYPE (bitvector_p), POINTER :: used_cats(:)
    INTEGER :: a, i, j, k, iters
    INTEGER :: instance, class, cat, attr
    REAL, ALLOCATABLE :: weights_per_cat(:, :), weight_by_class(:)
    REAL, ALLOCATABLE :: weight_yes_by_class(:), weight_no_by_class(:)
    TYPE (attribute), POINTER :: curr_attr, class_attr
    REAL :: weight, crit, maxcrit, gini_node
    REAL :: weight_yes_n, weight_yes_d, weight_no_n, weight_no_d
    REAL :: gini_n, gini_d, split_limit
    TYPE (bitvector), POINTER :: cats
    LOGICAL :: small_cat

    NULLIFY (nodeptr%catvar_classes, nodeptr%yes, nodeptr%no)

    class_attr => isptr%dd%attributes(opts%class_attribute_num)
    ALLOCATE (weight_by_class(class_attr%cat_count))
    weight_by_class = 0
    gini_n = 0
    gini_d = 0
    nodeptr%num_instances = 0
    DO k = 1, UBOUND(idx, 1)
      instance = idx(k)
      weight = bs%weight_sums_by_instance(instance)
      class = isptr%catvars(instance, class_attr%mapping)
      weight_by_class(class) = weight_by_class(class) + weight
      nodeptr%num_instances = nodeptr%num_instances + bs%instances(instance)
    END DO
    DO class = 1, class_attr%cat_count
      weight = weight_by_class(class)
      gini_n = gini_n + weight * weight
      gini_d = gini_d + weight
    END DO
    nodeptr%total_weight = gini_d

    IF (UBOUND(idx, 1).LT.opts%split_node_size) THEN
      ! Leaf node
      nonterminal = .FALSE.
    ELSE
      ! Find best split
      nodeptr%attribute = 0
      gini_node = gini_n / gini_d
      split_limit = gini_d * opts%split_ratio_limit
      maxcrit = gini_node

      DO i = 1, opts%split_variables
        j = rnd_int(UBOUND(isptr%dd%usedvars, 1))
        a = isptr%dd%usedvars(j)
        curr_attr => isptr%dd%attributes(a)
        attr = curr_attr%mapping
        IF (curr_attr%cat_count.EQ.contvar) THEN
          ! For continuous variables
          ALLOCATE (weight_yes_by_class(class_attr%cat_count), &
            & weight_no_by_class(class_attr%cat_count))
          weight_yes_n = 0
          weight_yes_d = 0
          weight_no_n = gini_n
          weight_no_d = gini_d
          weight_yes_by_class = 0
          weight_no_by_class = weight_by_class
          DO k = sort_offset, sort_offset + UBOUND(idx, 1) - 2
            instance = sorted_idx(k, curr_attr%usedcont)
            ! missing are at the bottom
            IF (sorted_idx(k + 1, curr_attr%usedcont).EQ.0) EXIT
            class = isptr%catvars(instance, class_attr%mapping)
            IF (class.EQ.missing_cat) CYCLE
            weight = bs%weight_sums_by_instance(instance)
            weight_yes_n = weight_yes_n + weight * (weight + &
              & 2 * weight_yes_by_class(class))
            weight_no_n = weight_no_n + weight * (weight - &
              & 2 * weight_no_by_class(class))
            weight_yes_d = weight_yes_d + weight
            weight_no_d = weight_no_d - weight
            weight_yes_by_class(class) = weight_yes_by_class(class) + weight
            weight_no_by_class(class) = weight_no_by_class(class) - weight
            IF ((isptr%contvar_ranks(instance, attr).LT. &
                & isptr%contvar_ranks(sorted_idx(k + 1, curr_attr%usedcont), &
                & attr)).AND.(MIN(weight_yes_d, weight_no_d).GT.split_limit)) &
                & THEN
              crit = weight_yes_n / weight_yes_d + weight_no_n / weight_no_d
              IF (crit.GT.maxcrit) THEN
                maxcrit = crit
                nodeptr%attribute = a
                ! The cutoff maximum is between this instance and the next one
                nodeptr%contvar_max = (isptr%contvars(instance, attr) &
                    & + isptr%contvars(sorted_idx(k + 1, curr_attr%usedcont), &
                    & attr)) / 2
              END IF
            END IF
          END DO
          DEALLOCATE (weight_yes_by_class, weight_no_by_class)
        ELSE
          ! For categorical variables (note: ignorevars are not in usedvars)
          ALLOCATE (weight_yes_by_class(class_attr%cat_count), &
            & weights_per_cat(class_attr%cat_count, curr_attr%cat_count))
          weights_per_cat = 0

          DO k = 1, UBOUND(idx, 1)
            instance = idx(k)
            cat = isptr%catvars(instance, attr)
            class = isptr%catvars(instance, class_attr%mapping)
            IF (cat.EQ.missing_cat.OR.class.EQ.missing_cat) CYCLE
            weights_per_cat(class, cat) = weights_per_cat(class, cat) + &
              & bs%weight_sums_by_instance(instance)
          END DO

          cats => new_bitvector(curr_attr%cat_count)
          small_cat = curr_attr%cat_count.LT.opts%big_catvar_cat_count
          IF (.NOT.small_cat) iters = 0
          DO
            IF (small_cat) THEN
              CALL increment_bitvector(cats)
              IF (getbit(cats, curr_attr%cat_count)) EXIT
              IF (is_and_bitvector(cats, used_cats(attr)%p)) CYCLE
            ELSE
              iters = iters + 1
              IF (iters.GT.opts%big_catvar_iterations) EXIT
              CALL randomise_bitvector(cats, used_cats(attr)%p)
            END IF
            weight_yes_by_class = 0

            DO cat = 1, curr_attr%cat_count
              IF (getbit(cats, cat)) THEN
                DO class = 1, class_attr%cat_count
                  weight_yes_by_class(class) = weight_yes_by_class(class) &
                    & + weights_per_cat(class, cat)
                END DO
              END IF
            END DO
            weight_yes_n = 0
            weight_yes_d = 0
            weight_no_n = 0
            DO class = 1, class_attr%cat_count
              weight = weight_yes_by_class(class)
              weight_yes_n = weight_yes_n + weight * weight
              weight_yes_d = weight_yes_d + weight
              weight = weight_by_class(class) - weight
              weight_no_n = weight_no_n + weight * weight
            END DO
            weight_no_d = (gini_d - weight_yes_d)
            IF (MIN(weight_yes_d, weight_no_d).EQ.0) CYCLE ! Non-splits
            crit = weight_yes_n / weight_yes_d + &
              & weight_no_n / weight_no_d
            IF ((crit.GT.maxcrit).AND. &
              & (MIN(weight_yes_d, weight_no_d).GT.split_limit)) THEN
              IF (ASSOCIATED(nodeptr%catvar_classes)) THEN
                CALL free_bitvector(nodeptr%catvar_classes)
              END IF
              nodeptr%catvar_classes => clone_bitvector(cats)
              maxcrit = crit
              nodeptr%attribute = a
            END IF
          END DO
          CALL free_bitvector(cats)
          DEALLOCATE (weights_per_cat, weight_yes_by_class)
        END IF
      END DO
      nonterminal = nodeptr%attribute.NE.0
      nodeptr%dgini = maxcrit - gini_node
      DEALLOCATE (weight_by_class)
    END IF
  END FUNCTION find_best_split

  FUNCTION most_frequent_class(isptr, idx)
    INTEGER :: most_frequent_class
    TYPE (instanceset), POINTER :: isptr
    INTEGER, POINTER :: idx(:)
    TYPE (attribute), POINTER :: class_attr
    class_attr => isptr%dd%attributes(opts%class_attribute_num)
    most_frequent_class = most_frequent(isptr%catvars(idx, &
      & class_attr%mapping), class_attr%cat_count)
  END FUNCTION most_frequent_class

  FUNCTION save_tree(thetree, filename, datadesc, leaf_id)
    TYPE (tree) :: thetree
    CHARACTER(LEN=*) :: filename
    TYPE (datadescription), POINTER :: datadesc
    INTEGER :: leaf_id(:)
    LOGICAL :: save_tree
    INTEGER :: handle
    LOGICAL :: err
    
    save_tree = .FALSE.
    handle = open_file(filename, .TRUE., err)
    IF (err) RETURN
    save_tree = .TRUE.
    
    WRITE(handle, *) thetree%terminal_bound
    WRITE(handle, *) thetree%leaf_bounds
    WRITE(handle, *) thetree%leaf_index
    WRITE(handle, *) thetree%leaf_pop
    WRITE(handle, *) thetree%oob%bits
    WRITE(handle, *) leaf_id

    CALL save_tree_internal(thetree%p, datadesc, handle)

    CALL close_file(handle)
  END FUNCTION save_tree

  FUNCTION load_tree(thetree, filename, datadesc, instance_count, leaf_id)
    TYPE (tree) :: thetree
    CHARACTER(LEN=*) :: filename
    TYPE (datadescription), POINTER :: datadesc
    INTEGER :: instance_count
    INTEGER, INTENT(OUT) :: leaf_id(:)
    LOGICAL :: load_tree
    INTEGER :: handle
    LOGICAL :: err
    
    load_tree = .FALSE.
    handle = open_file(filename, .FALSE., err)
    IF (err) RETURN
    load_tree = .TRUE.

    READ(handle, *) thetree%terminal_bound
    ALLOCATE (thetree%leaf_bounds(0:thetree%terminal_bound), &
      & thetree%leaf_index(instance_count), &
      & thetree%leaf_pop(thetree%terminal_bound))
    NULLIFY (thetree%dgini)
    READ(handle, *) thetree%leaf_bounds
    READ(handle, *) thetree%leaf_index
    READ(handle, *) thetree%leaf_pop
    thetree%oob => new_bitvector(instance_count)
    READ(handle, *) thetree%oob%bits
    READ(handle, *) leaf_id

    thetree%p => load_tree_internal(datadesc, handle)

    CALL close_file(handle)
  END FUNCTION load_tree

  RECURSIVE SUBROUTINE save_tree_internal(nodeptr, datadesc, handle)
    TYPE (node), POINTER :: nodeptr
    TYPE (datadescription), POINTER :: datadesc
    INTEGER :: handle
    
    IF (ASSOCIATED(nodeptr%yes)) THEN
      ! Branch node
      WRITE(handle, "(A1, 1X, I6, 1X, I6, 1X, F10.4, 1X, I6)") &
        & "b", nodeptr%id, nodeptr%attribute, nodeptr%total_weight, &
        & nodeptr%num_instances
      IF (datadesc%attributes(nodeptr%attribute)%cat_count.NE.contvar) THEN
        ! Catvar
        WRITE(handle, "((I12))") nodeptr%catvar_classes%bits
      ELSE
        ! Contvar
        WRITE(handle, "(F10.4)") nodeptr%contvar_max
      END IF
      CALL save_tree_internal(nodeptr%yes, datadesc, handle)
      CALL save_tree_internal(nodeptr%no, datadesc, handle)
    ELSE
      WRITE(handle, "(A1, 1X, I6, 1X, I6, 1X, F10.4, 1X, I6)") &
        & "l", nodeptr%id, nodeptr%attribute, nodeptr%total_weight, &
        & nodeptr%num_instances
    END IF
  END SUBROUTINE save_tree_internal

  RECURSIVE FUNCTION load_tree_internal(datadesc, handle) RESULT (nodeptr)
    TYPE (node), POINTER :: nodeptr
    TYPE (datadescription), POINTER :: datadesc
    INTEGER :: handle
    CHARACTER :: branch_or_leaf
    
    nodeptr => new_node()
    READ(handle, *) branch_or_leaf, nodeptr%id, nodeptr%attribute, &
      & nodeptr%total_weight, nodeptr%num_instances
    IF (branch_or_leaf.EQ."b") THEN
      ! Branch node
      IF (datadesc%attributes(nodeptr%attribute)%cat_count.NE.contvar) THEN
        ! Catvar
        nodeptr%catvar_classes => new_bitvector(datadesc%attributes( &
          & nodeptr%attribute)%cat_count)
        READ(handle, *) nodeptr%catvar_classes%bits
      ELSE
        ! Contvar
        READ(handle, *) nodeptr%contvar_max
      END IF
      nodeptr%yes => load_tree_internal(datadesc, handle)
      nodeptr%no => load_tree_internal(datadesc, handle)
    END IF
  END FUNCTION load_tree_internal

  SUBROUTINE print_tree(nodeptr, datadesc, handle)
    TYPE (node), POINTER :: nodeptr
    TYPE (datadescription), POINTER :: datadesc
    INTEGER :: handle

    CALL print_tree_internal(nodeptr, datadesc, 0, "R", handle)
  END SUBROUTINE print_tree

  RECURSIVE SUBROUTINE print_tree_internal(nodeptr, datadesc, level, head, &
    & handle)
    TYPE (node), POINTER :: nodeptr
    TYPE (datadescription), POINTER :: datadesc
    INTEGER :: level
    CHARACTER :: head
    INTEGER :: handle
    CHARACTER(LEN=linelen) :: str, wtstr
    INTEGER :: cat_base, i

    CALL squeeze_num(wtstr, nodeptr%total_weight)
    IF (ASSOCIATED(nodeptr%yes)) THEN
      ! Branch node
      IF (datadesc%attributes(nodeptr%attribute)%cat_count.NE.contvar) THEN
        ! Catvar
        str = ""
        cat_base = datadesc%attributes(nodeptr%attribute)%cat_start - 1
        DO i = 1, datadesc%attributes(nodeptr%attribute)%cat_count
          IF (getbit(nodeptr%catvar_classes, i)) THEN
            str = TRIM(str) // ", " // &
              & TRIM(datadesc%categories(cat_base + i)%name)
          END IF
        END DO
        str = "{ " // TRIM(str(3:)) // " }"
      ELSE
        CALL squeeze_num(str, nodeptr%contvar_max)
        str = TRIM(str) // " or less"
      END IF
      WRITE(handle, af) REPEAT("  ", level) // head // " (" // TRIM(wtstr) &
        & // "): " // TRIM(datadesc%attributes(nodeptr%attribute)%name) &
        & // " is " // TRIM(str) // "?"
      CALL print_tree_internal(nodeptr%yes, datadesc, level + 1, "Y", handle)
      CALL print_tree_internal(nodeptr%no, datadesc, level + 1, "N", handle)
    ELSE
      ! Leaf node
      WRITE(handle, af) REPEAT("  ", level) // head // " (" // TRIM(wtstr) &
        & // "): " // TRIM(datadesc%attributes(opts%class_attribute_num)%name) &
        & // " is [" // TRIM(datadesc%categories(datadesc%attributes( &
        & opts%class_attribute_num)%cat_start + nodeptr%attribute - 1)%name) &
        & // "]."
    END IF
  END SUBROUTINE print_tree_internal

  RECURSIVE FUNCTION classify_on_tree(isptr, instance, nodeptr, dgini) &
    & RESULT (leaf)
    TYPE (instanceset), POINTER :: isptr
    INTEGER :: instance
    TYPE (node), POINTER :: nodeptr
    REAL, POINTER, OPTIONAL :: dgini(:)
    TYPE (node), POINTER :: leaf
    TYPE (attribute), POINTER :: curr_attr
    TYPE (node), POINTER :: branch
    LOGICAL :: yes

    IF (ASSOCIATED(nodeptr%yes)) THEN
      ! Branch node
      curr_attr => isptr%dd%attributes(nodeptr%attribute)
      IF (curr_attr%cat_count.EQ.contvar) THEN
        ! Continuous variable
        yes = isptr%contvars(instance, curr_attr%mapping) &
          & .LE.nodeptr%contvar_max
      ELSE
        ! Categorical variable
        yes = getbit(nodeptr%catvar_classes, &
          & isptr%catvars(instance, curr_attr%mapping))
      END IF
      IF (yes) THEN
        branch => nodeptr%yes
      ELSE
        branch => nodeptr%no
      END IF
      IF (PRESENT(dgini)) THEN
        dgini(nodeptr%attribute) = &
          & dgini(nodeptr%attribute) + nodeptr%dgini
        leaf => classify_on_tree(isptr, instance, branch, dgini)
      ELSE
        leaf => classify_on_tree(isptr, instance, branch)
      END IF
    ELSE
      ! Leaf node
      leaf => nodeptr
    END IF
  END FUNCTION classify_on_tree

  RECURSIVE FUNCTION classify_permuted_on_tree(isptr, instance, nodeptr, &
    & perm_attr, straight, permutation) RESULT (leaf)
    TYPE (instanceset), POINTER :: isptr
    INTEGER :: instance
    TYPE (node), POINTER :: nodeptr
    INTEGER :: perm_attr
    INTEGER, POINTER :: straight(:), permutation(:)
    TYPE (node), POINTER :: leaf
    TYPE (attribute), POINTER :: curr_attr
    TYPE (node), POINTER :: branch
    LOGICAL :: yes
    INTEGER :: real_instance

    IF (ASSOCIATED(nodeptr%yes)) THEN
      ! Branch node
      IF (nodeptr%attribute.EQ.perm_attr) THEN
        real_instance = permutation(instance)
      ELSE
        real_instance = straight(instance)
      END IF

      curr_attr => isptr%dd%attributes(nodeptr%attribute)
      IF (curr_attr%cat_count.EQ.contvar) THEN
        ! Continuous variable
        yes = isptr%contvars(real_instance, curr_attr%mapping) &
          & .LE.nodeptr%contvar_max
      ELSE
        ! Categorical variable
        yes = getbit(nodeptr%catvar_classes, &
          & isptr%catvars(real_instance, curr_attr%mapping))
      END IF
      IF (yes) THEN
        branch => nodeptr%yes
      ELSE
        branch => nodeptr%no
      END IF
      leaf => classify_permuted_on_tree(isptr, instance, branch, &
        & perm_attr, straight, permutation)
    ELSE
      ! Leaf node
      leaf => nodeptr
    END IF
  END FUNCTION classify_permuted_on_tree

  RECURSIVE SUBROUTINE fill_node_data(nodeptr, tags, weights, used_attrs)
    TYPE (node), POINTER :: nodeptr
    INTEGER, POINTER :: tags(:)
    REAL, POINTER :: weights(:)
    TYPE (bitvector), POINTER :: used_attrs

    IF (ASSOCIATED(nodeptr%yes)) THEN
      ! Branch node
      CALL setbit(used_attrs, nodeptr%attribute, .TRUE.)
      CALL fill_node_data(nodeptr%yes, tags, weights, used_attrs)
      CALL fill_node_data(nodeptr%no, tags, weights, used_attrs)
    ELSE
      ! Leaf node
      tags(nodeptr%id) = nodeptr%attribute
      weights(nodeptr%id) = nodeptr%total_weight / nodeptr%num_instances
    END IF
  END SUBROUTINE fill_node_data
END MODULE trees
