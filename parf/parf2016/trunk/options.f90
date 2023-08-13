MODULE options
  USE utilities
  USE parallel
  IMPLICIT NONE

  INTEGER, PARAMETER :: optlen = 256
  INTEGER, PARAMETER :: class_last = -1, class_new = -2, class_named = 0
  INTEGER, PARAMETER :: tokenlen = 32, linelen = 1024

  INTEGER, PARAMETER :: trainset_type = 0, testset_type = 1, protoset_type = 2

  TYPE optionlist
    CHARACTER(LEN=optlen) :: trainset, testset
    INTEGER :: fill_passes
    INTEGER :: redo_with_important_vars
    REAL :: redo_with_significant_vars
    CHARACTER(LEN=optlen) :: class_attribute, positive_category
    INTEGER :: class_attribute_num, positive_category_num
    INTEGER :: split_variables, split_variables_important
    INTEGER :: split_node_size
    REAL :: split_ratio_limit, new_class_quantity
    REAL :: outlier_cutoff
    CHARACTER(LEN=optlen) :: usedvars, unusedvars
    LOGICAL :: used_default
    INTEGER :: big_catvar_cat_count
    INTEGER :: big_catvar_iterations 
    INTEGER :: num_trees, scaling_divergence, labelvar
    INTEGER :: num_prox, num_scale, num_prot, prox_for_prot
    LOGICAL :: verbose, verboseall, summary
    CHARACTER(LEN=optlen) :: dump_forest, save_forest, load_forest
    CHARACTER(LEN=optlen) :: fast_importances, case_importances
    CHARACTER(LEN=optlen) :: importances, interaction
    CHARACTER(LEN=optlen) :: train_outliers, test_outliers
    CHARACTER(LEN=optlen) :: train_scaling, test_scaling, proto_scaling
    CHARACTER(LEN=optlen) :: prototypes, prototype_analysis
    CHARACTER(LEN=optlen) :: train_confusion, test_confusion
    CHARACTER(LEN=optlen) :: train_votes, test_votes, test_results
    CHARACTER(LEN=optlen) :: train_test_arff, test_arff
    CHARACTER(LEN=optlen) :: class_weights, label
    CHARACTER(LEN=optlen) :: graphics, graphics_term
    LOGICAL :: calc_importances, calc_case_importances
    LOGICAL :: last_prox_required, do_graphics, calc_test_prox
  END TYPE optionlist
  INTEGER :: rnd_seed
  INTEGER, PRIVATE :: argc

  TYPE (optionlist) :: opts
CONTAINS
  FUNCTION parse_options()
    INTEGER :: IARGC
    CHARACTER(LEN=optlen) :: arg
    LOGICAL :: parse_options
    LOGICAL :: show_rng_seed
    INTEGER :: i, pos
    INTEGER :: timeseed ! C function

!
! Here options are parsed.
! All PARF options are parsed on the main processor (par_front = par_rank.EQ.0)
! and then distributed to all other participating processors
! 
    IF (.NOT.par_front) GO TO 9999   ! Just do the broadcast reception

    parse_options = .FALSE.
    show_rng_seed = .FALSE.

    opts%do_graphics = .FALSE.
    opts%calc_test_prox = .FALSE.
    opts%summary = .TRUE.
    opts%graphics = ""
    opts%graphics_term = ""
    opts%fill_passes = 1 ! 0 none, 1 rough, n>1 full with n passes
    opts%redo_with_important_vars = 0
    opts%redo_with_important_vars = 0
    opts%class_attribute_num = -1 ! -1 last, -2 new
    opts%new_class_quantity = 1.0
    opts%split_variables = -1 ! -1 means sqrt(non-ignored-vars)
    opts%split_variables_important = -1 ! -1 means sqrt(non-ignored-vars)
    opts%split_node_size = 2 ! i.e. split every node you can
    opts%split_ratio_limit = 0 ! i.e. any split is okay
    opts%usedvars = ""
    opts%class_attribute = ""
    opts%positive_category = ""
    opts%unusedvars = ""
    opts%label = ""
    opts%test_results = ""
    opts%test_arff = ""
    opts%train_test_arff = ""
    opts%big_catvar_cat_count = 12
    opts%big_catvar_iterations = 0 ! 2 ** (big_catvar_cat_count - 1)
    opts%trainset = ""
    opts%class_weights = ""
    opts%testset = ""
    opts%num_trees = MAX(par_processes, 100)
    opts%dump_forest = ""
    opts%save_forest = ""
    opts%load_forest = ""
    opts%fast_importances = ""
    opts%importances = ""
    opts%interaction = ""
    opts%case_importances = ""
    opts%train_outliers = ""
    opts%test_outliers = ""
    opts%train_confusion = ""
    opts%test_confusion = ""
    opts%train_votes = ""
    opts%test_votes = ""
    opts%verbose = .FALSE.
    opts%verboseall = .FALSE.
    opts%num_prox = 0
    opts%prox_for_prot = 0
    opts%num_scale = 0
    opts%scaling_divergence = 10
    opts%num_prot = 0
    opts%train_scaling = ""
    opts%test_scaling = ""
    opts%proto_scaling = ""
    opts%prototypes = ""
    opts%prototype_analysis = ""
    opts%outlier_cutoff = -HUGE(1.0)
    rnd_seed = 0

    argc = IARGC()
    IF (argc.EQ.0) THEN
      parse_options = .FALSE.
      IF (par_front) CALL print_optionhelp()
      RETURN
    END IF
    i = 1
    DO WHILE (i.LE.argc)
      CALL getarg(i, arg)
      SELECT CASE (TRIM(arg))
        CASE ("-t") ! trainset file
          i = i + 1
          IF (nextarg(i, arg)) opts%trainset = arg
        CASE ("-tv") ! trainset votes file
          i = i + 1
          IF (nextarg(i, arg, "-")) opts%train_votes = arg
        CASE ("-tc") ! trainset confusion matrix file
          i = i + 1
          IF (nextarg(i, arg, "-")) opts%train_confusion = arg
        CASE ("-a") ! testset file
          i = i + 1
          IF (nextarg(i, arg)) opts%testset = arg
        CASE ("-av") ! testset votes file
          i = i + 1
          IF (nextarg(i, arg, "-")) opts%test_votes = arg
        CASE ("-ac") ! testset confusion matrix file
          i = i + 1
          IF (nextarg(i, arg, "-")) opts%test_confusion = arg
        CASE ("-ar") ! testset classification results file
          i = i + 1
          IF (nextarg(i, arg, "-")) opts%test_results = arg
        CASE ("-aa") ! testset arff output file
          i = i + 1
          IF (nextarg(i, arg, "-")) opts%test_arff = arg
        CASE ("-ta") ! trainset+testset arff output file
          i = i + 1
          IF (nextarg(i, arg, "-")) opts%train_test_arff = arg
        CASE ("-u") ! (un)used vars
          i = i + 1
          IF (nextarg(i, arg)) opts%usedvars = arg
        CASE ("-uu") ! (un)used vars
          i = i + 1
          IF (nextarg(i, arg)) opts%unusedvars = arg
        CASE ("-ri") ! row identifier (=label)
          i = i + 1
          IF (nextarg(i, arg)) opts%label = arg
        CASE ("-g") ! graphics
          i = i + 1
          IF (nextarg(i, arg, "-")) THEN
            opts%graphics = arg
            opts%do_graphics = .TRUE.
          END IF
        CASE ("-gt") ! graphics terminal type
          i = i + 1
          IF (nextarg(i, arg)) opts%graphics_term = arg
        CASE ("-n") ! number of trees to grow
          i = i + 1
          IF (nextarg(i, arg)) READ(arg, *) opts%num_trees
        CASE ("-cp") ! positive category
          i = i + 1
          IF (nextarg(i, arg)) READ(arg, *) opts%positive_category
        CASE ("-cq") ! new class quantity
          i = i + 1
          IF (nextarg(i, arg, "1")) THEN
            pos = INDEX(arg, "%")
            IF (pos.NE.0) arg = arg(1:pos-1)
            READ(arg, *) opts%new_class_quantity
            IF (pos.NE.0) &
              & opts%new_class_quantity = opts%new_class_quantity / 100
          END IF
        CASE ("-c") ! class attribute: NEW or LAST or name
          i = i + 1
          IF (nextarg(i, arg)) THEN
            opts%class_attribute = lcase(arg)
            SELECT CASE (TRIM(arg))
              CASE ("last")
                opts%class_attribute_num = class_last
              CASE ("new")
                opts%class_attribute_num = class_new
              CASE DEFAULT
                opts%class_attribute_num = class_named
            END SELECT
          END IF
        CASE ("-p") ! number of proximate instances to take into consideration
          i = i + 1
          IF (nextarg(i, arg, "100%")) THEN
            pos = INDEX(arg, "%")
            IF (pos.NE.0) arg = arg(1:pos-1)
            READ(arg, *) opts%num_prox
            IF (pos.NE.0) opts%num_prox = -opts%num_prox ! % are negative
          END IF
        CASE ("-w") ! class weight override
          i = i + 1
          IF (nextarg(i, arg)) opts%class_weights = arg
        CASE ("-st") ! training set scaling
          i = i + 1
          IF (nextarg(i, arg, "-")) opts%train_scaling = arg
        CASE ("-sa") ! test set scaling
          i = i + 1
          IF (nextarg(i, arg, "-")) opts%test_scaling = arg
          opts%calc_test_prox = .TRUE.
        CASE ("-sy") ! prototype set scaling
          i = i + 1
          IF (nextarg(i, arg, "-")) opts%proto_scaling = arg
        CASE ("-s") ! number of scaling coordinates
          i = i + 1
          IF (nextarg(i, arg, "2")) READ(arg, *) opts%num_scale
        CASE ("-sd") ! scaling divergence limit
          i = i + 1
          IF (nextarg(i, arg, "0")) READ(arg, *) opts%scaling_divergence
        CASE ("-y") ! prototypes
          i = i + 1
          IF (nextarg(i, arg, "-")) opts%prototypes = arg
        CASE ("-ya") ! prototype analysis
          i = i + 1
          IF (nextarg(i, arg, "-")) opts%prototype_analysis = arg
        CASE ("-yn") ! number of prototypes
          i = i + 1
          IF (nextarg(i, arg)) READ(arg, *) opts%num_prot
        CASE ("-yp") ! number of proximate instances to look for prototype
          i = i + 1
          IF (nextarg(i, arg, "100%")) THEN
            pos = INDEX(arg, "%")
            IF (pos.NE.0) arg = arg(1:pos-1)
            READ(arg, *) opts%prox_for_prot
            IF (pos.NE.0) opts%prox_for_prot = -opts%prox_for_prot
          END IF
        CASE ("-b") ! big categorical limit
          i = i + 1
          IF (nextarg(i, arg)) READ(arg, *) opts%big_catvar_cat_count
        CASE ("-bi") ! big categorical iterations
          i = i + 1
          IF (nextarg(i, arg)) READ(arg, *) opts%big_catvar_iterations
        CASE ("-f") ! fill mode: 0 none, 1 rough, 2+ passes
          i = i + 1
          IF (nextarg(i, arg)) READ(arg, *) opts%fill_passes
        CASE ("-v") ! number of most important variables to redo
          i = i + 1
          IF (nextarg(i, arg)) READ(arg, *) opts%redo_with_important_vars
        CASE ("-vs") ! the maximum significance of variables to redo
          i = i + 1
          IF (nextarg(i, arg)) READ(arg, *) opts%redo_with_significant_vars
        CASE ("-m") ! number of variables to split on
          i = i + 1
          IF (nextarg(i, arg)) READ(arg, *) opts%split_variables
        CASE ("-mv") ! number of variables to split on in most important var pass
          i = i + 1
          IF (nextarg(i, arg)) READ(arg, *) opts%split_variables_important
        CASE ("-xs") ! minimum node size to split
          i = i + 1
          IF (nextarg(i, arg)) THEN
            pos = INDEX(arg, "%")
            IF (pos.NE.0) arg = arg(1:pos-1)
            READ(arg, *) opts%split_node_size
            IF (pos.NE.0) opts%split_node_size = -opts%split_node_size
          END IF
        CASE ("-xr") ! split cutoff ratio limit
          i = i + 1
          IF (nextarg(i, arg)) THEN
            pos = INDEX(arg, "%")
            IF (pos.NE.0) arg = arg(1:pos-1)
            READ(arg, *) opts%split_ratio_limit
            IF (pos.NE.0) opts%split_ratio_limit = opts%split_ratio_limit / 100
          END IF
        CASE ("-r") ! random number generator seed
          i = i + 1
          IF (nextarg(i, arg, "show")) THEN
            IF (TRIM(arg).NE."show") THEN
              READ(arg, *) rnd_seed
            ELSE
              show_rng_seed = .TRUE.
            END IF
          END IF
        CASE ("-h", "--help") ! help
          parse_options = .FALSE.
          IF (par_front) CALL print_optionhelp()
          RETURN
        CASE ("-fs") ! save forest
          i = i + 1
          IF (nextarg(i, arg)) opts%save_forest = arg
        CASE ("-fl") ! load forest (use an old one)
          i = i + 1
          IF (nextarg(i, arg)) opts%load_forest = arg
        CASE ("-fd") ! text dump of the forest
          i = i + 1
          IF (nextarg(i, arg, "-")) opts%dump_forest = arg
        CASE ("-i") ! importances
          i = i + 1
          IF (nextarg(i, arg, "-")) opts%importances = arg
        CASE ("-ic") ! case importances
          i = i + 1
          IF (nextarg(i, arg, "-")) opts%case_importances = arg
        CASE ("-ii") ! interaction
          i = i + 1
          IF (nextarg(i, arg, "-")) opts%interaction = arg
        CASE ("-if") ! fast importances
          i = i + 1
          IF (nextarg(i, arg, "-")) opts%fast_importances = arg
        CASE ("-ot") ! train outliers
          i = i + 1
          IF (nextarg(i, arg, "-")) opts%train_outliers = arg
        CASE ("-oa") ! test outliers
          i = i + 1
          IF (nextarg(i, arg, "-")) opts%test_outliers = arg
          opts%calc_test_prox = .TRUE.
        CASE ("-o") ! outlier cutoff
          i = i + 1
          IF (nextarg(i, arg, "0")) READ(arg, *) opts%outlier_cutoff
        CASE ("--verbose")
          opts%verbose = par_front
          opts%verboseall = .TRUE.
        CASE DEFAULT
          ! commented out for MPI (MPI_INIT does not work well with GETARG)
          ! IF (par_front) GO TO 9000
      END SELECT
      i = i + 1
    END DO
   
    ! dependencies
    IF (LEN_TRIM(opts%testset).EQ.0) THEN
      opts%test_outliers = ""
      opts%test_scaling = ""
      opts%test_confusion = ""
      opts%test_votes = ""
      opts%test_results = ""
      opts%test_arff = ""
      opts%train_test_arff = ""
      opts%calc_test_prox = .FALSE.
    END IF
    opts%calc_case_importances = LEN_TRIM(opts%case_importances).NE.0
    opts%calc_importances = opts%calc_case_importances &
      & .OR.opts%redo_with_important_vars.NE.0 &
      & .OR.opts%redo_with_significant_vars.NE.0 &
      & .OR.LEN_TRIM(opts%importances).NE.0.OR.LEN_TRIM(opts%interaction).NE.0
    IF ((LEN_TRIM(opts%train_scaling).NE.0 &
      & .OR.LEN_TRIM(opts%test_scaling).NE.0 &
      & .OR.LEN_TRIM(opts%proto_scaling).NE.0) &
      & .AND.opts%num_scale.EQ.0) &
        & opts%num_scale = 2
    IF (LEN_TRIM(opts%train_scaling).EQ.0 &
      & .AND.LEN_TRIM(opts%test_scaling).EQ.0 &
      & .AND.LEN_TRIM(opts%proto_scaling).EQ.0.AND.opts%num_scale.NE.0) THEN
      IF (LEN_TRIM(opts%testset).EQ.0) THEN
        opts%train_scaling = "-"
      ELSE
        opts%test_scaling = "-"
        opts%calc_test_prox = .TRUE.
      END IF
    END IF
    IF ((LEN_TRIM(opts%prototypes).NE.0 &
      & .OR.LEN_TRIM(opts%prototype_analysis).NE.0 &
      & .OR.LEN_TRIM(opts%proto_scaling).NE.0).AND.opts%num_prot.EQ.0) &
        & opts%num_prot = 1
    IF (opts%num_prot.GT.0.AND.LEN_TRIM(opts%trainset).EQ.0) THEN
      WRITE(0, af) "Warning: Can't calculate proximities from a saved forest."
      opts%num_prot = 0
      opts%prototypes = ""
      opts%prototype_analysis = ""
      opts%proto_scaling = ""
    END IF
    IF (LEN_TRIM(opts%prototypes).EQ.0 &
      & .AND.LEN_TRIM(opts%prototype_analysis).EQ.0 &
      & .AND.LEN_TRIM(opts%proto_scaling).EQ.0.AND.opts%num_prot.NE.0) &
        & opts%prototype_analysis = "-"
    opts%last_prox_required = LEN_TRIM(opts%train_outliers).NE.0 &
      & .OR.LEN_TRIM(opts%test_outliers).NE.0 &
      & .OR.opts%num_scale.GT.0.OR.opts%num_prot.GT.0
    IF ((opts%last_prox_required.OR.opts%fill_passes.GT.1) &
      & .AND.opts%num_prox.EQ.0) THEN
      opts%num_prox = -100
    END IF
    IF (opts%num_prot.NE.0.AND.opts%prox_for_prot.EQ.0) THEN
      opts%prox_for_prot = -50
    END IF
    IF (LEN_TRIM(opts%testset).EQ.0) THEN
      opts%test_votes = ""
    END IF
    IF (LEN_TRIM(opts%graphics_term).NE.0.AND.LEN_TRIM(opts%graphics).EQ.0) THEN
      opts%graphics = "-"
      opts%do_graphics = .TRUE.
    END IF
    IF (TRIM(opts%graphics).EQ."-") THEN
      opts%verbose = .FALSE.
      opts%summary = .FALSE.
    END IF

    ! Both -u and -uu?
    opts%used_default = (LEN_TRIM(opts%unusedvars).NE.0).OR. &
      & (LEN_TRIM(opts%usedvars).EQ.0)
    IF (opts%big_catvar_iterations.EQ.0) THEN
      opts%big_catvar_iterations = 2 ** (opts%big_catvar_cat_count - 1)
    END IF
    GO TO 9999

    ! Errors
    9000 WRITE(0, af) "Syntax error: Unknown option: " // TRIM(arg)
    GO TO 9999

    9999 CONTINUE

    ! All processors do the following for parse_options()
    parse_options = .TRUE.

    CALL broadcast_options()

    ! if seed is not specified, make it truly random
    IF (par_front) THEN
      IF (rnd_seed.EQ.0) rnd_seed = timeseed()
      IF (opts%verbose.OR.show_rng_seed) &
        & WRITE(6, "(A6, I12)") "Seed: ", rnd_seed
    END IF
    CALL par_bcast_int(rnd_seed)
    CALL seed_rnd(rnd_seed)
  END FUNCTION parse_options

  FUNCTION nextarg(pos, arg, default_value)
    INTEGER :: pos
    CHARACTER(LEN=optlen), INTENT(OUT) :: arg
    CHARACTER(LEN=*), OPTIONAL :: default_value
    CHARACTER(LEN=optlen) :: tmp
    LOGICAL :: nextarg
    nextarg = pos .LE. argc
    IF (nextarg) THEN
      CALL getarg(pos, tmp)
      IF (tmp(1:1).EQ."-".AND.LEN_TRIM(tmp).NE.1 &
        & .AND.PRESENT(default_value)) THEN
        arg = default_value
        pos = pos - 1
      ELSE
        arg = tmp
      END IF
    ELSE IF (PRESENT(default_value)) THEN
      arg = default_value
      nextarg = .TRUE.
    ELSE
      WRITE(0, af) "Syntax error: Missing argument"
    END IF
  END FUNCTION

  SUBROUTINE print_optionhelp()
    WRITE(*, af) "PARF (C) 2005 Rudjer Boskovic Institute"
    WRITE(*, af) "Goran Topic, Tomislav Smuc; algorithm by Leo Breiman and&
      & Adele Cutler"
    WRITE(*, af) "Licensed under GNU GPL 2.0"
    WRITE(*, *)
    WRITE(*, af) "Usage: rf [OPTION...]"
    WRITE(*, af) "-h | --help   show this message"
    WRITE(*, af) "-t file       file to use as training set"
    WRITE(*, af) "-a file       file to analyse and classify"
    WRITE(*, af) "-tv [file]    training set votes output file"
    WRITE(*, af) "-tc [file]    training set confusion matrix output file"
    WRITE(*, af) "-av [file]    test set votes output file"
    WRITE(*, af) "-ac [file]    test set confusion matrix output file"
    WRITE(*, af) "-ar [file]    test set classification results output file"
    WRITE(*, af) "-aa [file]    test set ARFF output file"
    WRITE(*, af) "-ta [file]    train + test set ARFF output file"
    WRITE(*, af) "-c class      the class attribute, or NEW, or LAST (default)"
    WRITE(*, af) "-cq [n[%]]    quantity of generated class instances (only&
      & with -c NEW)"
    WRITE(*, af) "-cp category  positive category"
    WRITE(*, af) "-n trees      the number of trees to grow"
    WRITE(*, af) "-f n          the fill method: 0=none, 1=rough,&
      & 2+=# of passes"
    WRITE(*, af) "-v n          redo the forest with n most important variables"
    WRITE(*, af) "-vs n         redo the forest with variables more&
      & significant than n"
    WRITE(*, af) "-p [n[%]]     number of proximate cases to take into&
      & consideration"
    WRITE(*, af) "-st [file]    training set scaling output file"
    WRITE(*, af) "-sa [file]    test set scaling output file"
    WRITE(*, af) "-sy [file]    prototype scaling output file"
    WRITE(*, af) "-s n          number of scaling coordinates (default 2)"
    WRITE(*, af) "-sd [n]       max allowed divergent iterations in scaling&
      & calculation"
    WRITE(*, af) "-y [file]     prototypes output file"
    WRITE(*, af) "-ya [file]    detailed prototypes analysis output file"
    WRITE(*, af) "-yn n         number of prototypes per class (if available)"
    WRITE(*, af) "-yp n         number of proximates to look at for prototypes"
    WRITE(*, af) "-m mvar       the number of variables to split on"
    WRITE(*, af) "-mv mvar      the number of split variables with most&
      & important variables"
    WRITE(*, af) "-xs size[%]   the minimum node size to split"
    WRITE(*, af) "-xr ratio[%]  the split cutoff ratio limit (0-1)"
    WRITE(*, af) "-b categories the number of categories that use the fast&
      & split algorithm"
    WRITE(*, af) "-bi iters     the number of fast split iterations per&
      & category"
    WRITE(*, af) "-u(u) var,... comma-separated list of used or unused&
      & attributes"
    WRITE(*, af) "-ri var       row identifier variable"
    WRITE(*, af) "-r seed       random number generator seed"
    WRITE(*, af) "-fd [file]    dump the forest as a text"
    WRITE(*, af) "-fs file      save the forest"
    WRITE(*, af) "-fl file      load the forest"
    WRITE(*, af) "-i [file]     variable importances output file"
    WRITE(*, af) "-ic [file]    case-by-case variable importances output file"
    WRITE(*, af) "-ii [file]    variable interaction output file"
    WRITE(*, af) "-if [file]    fast variable importances output file"
    WRITE(*, af) "-ot [file]    training set outlier measure output file" 
    WRITE(*, af) "-oa [file]    testing set outlier measure output file" 
    WRITE(*, af) "-w wt,...     class weight override"
    WRITE(*, af) "-g [file]     generate gnuplot graphics script"
    WRITE(*, af) "-gt type      gnuplot graphics terminal type"
    WRITE(*, af) "--verbose     print what is done"
  END SUBROUTINE

!
! The bcast_options is a specific PARF function which broadcasts options as a
! string of broadcasts for each element of the PARF options structure.
!

  SUBROUTINE broadcast_options()
    CALL par_bcast_char(opts%trainset)
    CALL par_bcast_char(opts%testset)
    CALL par_bcast_int(opts%fill_passes)
    CALL par_bcast_int(opts%redo_with_important_vars)
    CALL par_bcast_float(opts%redo_with_significant_vars)
    CALL par_bcast_char(opts%class_attribute)
    CALL par_bcast_char(opts%positive_category)
    CALL par_bcast_int(opts%class_attribute_num)
    CALL par_bcast_float(opts%new_class_quantity)
    CALL par_bcast_int(opts%split_variables)
    CALL par_bcast_int(opts%split_variables_important)
    CALL par_bcast_int(opts%split_node_size)
    CALL par_bcast_float(opts%split_ratio_limit)
    CALL par_bcast_float(opts%outlier_cutoff)
    CALL par_bcast_char(opts%usedvars)
    CALL par_bcast_char(opts%unusedvars)
    CALL par_bcast_char(opts%label)
    CALL par_bcast_bool(opts%used_default)
    CALL par_bcast_int(opts%big_catvar_cat_count)
    CALL par_bcast_int(opts%big_catvar_iterations)
    CALL par_bcast_int(opts%num_trees)
    CALL par_bcast_int(opts%scaling_divergence)
    CALL par_bcast_int(opts%num_prox)
    CALL par_bcast_int(opts%num_scale)
    CALL par_bcast_int(opts%num_prot)
    CALL par_bcast_int(opts%prox_for_prot)
    CALL par_bcast_bool(opts%verboseall)
    CALL par_bcast_bool(opts%summary)
    CALL par_bcast_char(opts%dump_forest)
    CALL par_bcast_char(opts%save_forest)
    CALL par_bcast_char(opts%load_forest)
    CALL par_bcast_char(opts%fast_importances)
    CALL par_bcast_char(opts%case_importances)
    CALL par_bcast_char(opts%importances)
    CALL par_bcast_char(opts%interaction)
    CALL par_bcast_char(opts%train_outliers)
    CALL par_bcast_char(opts%test_outliers)
    CALL par_bcast_char(opts%train_scaling)
    CALL par_bcast_char(opts%test_scaling)
    CALL par_bcast_char(opts%proto_scaling)
    CALL par_bcast_char(opts%prototypes)
    CALL par_bcast_char(opts%prototype_analysis)
    CALL par_bcast_char(opts%train_confusion)
    CALL par_bcast_char(opts%test_confusion)
    CALL par_bcast_char(opts%train_votes)
    CALL par_bcast_char(opts%test_votes)
    CALL par_bcast_char(opts%test_results)
    CALL par_bcast_char(opts%train_test_arff)
    CALL par_bcast_char(opts%test_arff)
    CALL par_bcast_char(opts%class_weights)
    CALL par_bcast_char(opts%graphics)
    CALL par_bcast_char(opts%graphics_term)
    CALL par_bcast_bool(opts%calc_importances)
    CALL par_bcast_bool(opts%calc_case_importances)
    CALL par_bcast_bool(opts%last_prox_required)
    CALL par_bcast_bool(opts%do_graphics)
    CALL par_bcast_bool(opts%calc_test_prox)
    CALL par_bcast_int(rnd_seed)
  END SUBROUTINE broadcast_options
END MODULE options
