      INTEGER FUNCTION mpir_iargc()
        mpir_iargc = iargc()
      END

      SUBROUTINE mpir_getarg(i, s)
        INTEGER i
        CHARACTER*(*) s
        CALL getarg(i,s)
      END
