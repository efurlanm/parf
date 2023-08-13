MODULE graphics
  USE options
  USE utilities

  IMPLICIT NONE

  LOGICAL :: graphics_together ! can be global, as it is never nested
  INTEGER :: graphics_handle
  LOGICAL :: graphics_has_term
  INTEGER :: graphics_num = 0
  CHARACTER(LEN=linelen) :: plot_options

CONTAINS
  SUBROUTINE init_graphics()
    INTEGER :: pos
    LOGICAL :: err

    IF (opts%do_graphics) THEN
      graphics_handle = open_file(opts%graphics, .TRUE., err)
      IF (err) THEN
        opts%do_graphics = .FALSE.
        opts%graphics = ""
        RETURN
      END IF
      graphics_has_term = LEN_TRIM(opts%graphics_term).NE.0
      IF (graphics_has_term) THEN
        WRITE(graphics_handle, af) "set term " // TRIM(opts%graphics_term)
        IF (opts%graphics_term(1:10).EQ."postscript") THEN
          IF (opts%graphics_term(12:14).EQ."eps") THEN
            opts%graphics_term = "eps"
          ELSE
            opts%graphics_term = "ps"
          END IF
        ELSE
          pos = SCAN(opts%graphics_term, " ")
          IF (pos.NE.0) opts%graphics_term = opts%graphics_term(1:pos)
        END IF
      END IF
    END IF
  END SUBROUTINE init_graphics

  SUBROUTINE finish_graphics()
    IF (opts%do_graphics) THEN
      IF (.NOT.graphics_has_term) THEN
        WRITE (graphics_handle, af) "pause -1 'Press Enter'"
      END IF
      CALL close_file(graphics_handle)
    END IF
  END SUBROUTINE finish_graphics

  FUNCTION open_graphics(filename, err, default_outfile) RESULT (handle)
    INTEGER :: handle
    CHARACTER(LEN=*), INTENT(OUT) :: filename
    CHARACTER(LEN=*) :: default_outfile
    LOGICAL, INTENT(OUT) :: err
    CHARACTER(LEN=optlen) :: real_outfile
    INTEGER :: pos
    
    IF (filename.NE."-") THEN
      pos = SCAN(filename, ".", .TRUE.)
      IF (pos.EQ.0) THEN
        real_outfile = filename
      ELSE
        real_outfile = filename(1:pos-1)
      END IF
    ELSE
      real_outfile = default_outfile
    END IF
    
    IF (opts%do_graphics) THEN
      graphics_together = opts%graphics.EQ.filename.OR.TRIM(filename).EQ."-"
      WRITE(graphics_handle, af) "reset"
      IF (graphics_has_term) THEN
        WRITE(graphics_handle, af) "set output '" // TRIM(real_outfile) &
          & // "." // TRIM(opts%graphics_term) // "'"
      ELSE
        graphics_num = graphics_num + 1
        WRITE(graphics_handle, "(A17, I2)") "set terminal x11 ", graphics_num
      END IF
    ELSE
      graphics_together = .FALSE.
    END IF
    IF (graphics_together) THEN
      handle = graphics_handle
      filename = "-"
      err = .FALSE.
    ELSE
      handle = open_file(filename, .TRUE., err)
    END IF
  END FUNCTION open_graphics

  SUBROUTINE plot_graphics(plot)
    CHARACTER(LEN=*) :: plot

    IF (.NOT.opts%do_graphics) RETURN
    IF (graphics_together) THEN
      WRITE(graphics_handle, af) TRIM(plot)
    ELSE
      plot_options = plot
    END IF
  END SUBROUTINE plot_graphics

  SUBROUTINE close_graphics(handle)
    INTEGER :: handle

    IF (.NOT.opts%do_graphics) RETURN
    IF (graphics_together) THEN
      WRITE (handle, af) "end"
    ELSE
      CALL close_file(handle)
      WRITE (graphics_handle, af) TRIM(plot_options)
    END IF
  END SUBROUTINE close_graphics
END MODULE
