c-----------------------------------------------------------------------
c     file defs.f
c     local defintions for most computers.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. defs_mod.
c-----------------------------------------------------------------------
      module defs_mod
      implicit none
      
      integer, parameter ::
     $     r4=selected_real_kind(6,37),
     $     r8=selected_real_kind(15,307),
     $     r16=selected_real_kind(33,4931),
     $     i4=selected_int_kind(9),
     $     i8=selected_int_kind(18)

      real(r8), parameter :: pi = 3.14159265358979323846_r8
      real(r8), parameter :: twopi = 6.28318530717958647692_r8
      real(r8), parameter :: alog10 = 2.3025850929940459_r8
      real(r8), parameter :: rtod = 180.0_r8 / pi
      real(r8), parameter :: zero = 0.0_r8
      real(r8), parameter :: one = 1.0_r8
      real(r8), parameter :: two = 2.0_r8
      real(r8), parameter :: half = 0.5_r8
      integer, parameter :: bin_unit = 10
      integer, parameter :: out_unit = 11

      contains
c-----------------------------------------------------------------------
c     subprogram 1. timer.
c     handles machine-dependent timing statistics.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     delcarations.
c-----------------------------------------------------------------------
      SUBROUTINE timer(mode,unit,op_cpuseconds,op_wallseconds)
      
      INTEGER, INTENT(IN) :: mode,unit
      REAL(r4), INTENT(OUT), OPTIONAL :: op_cpuseconds,op_wallseconds

      INTEGER(8), SAVE :: count_rate, wall_start
      REAL(r4), SAVE :: start
      REAL(r4) :: seconds
      INTEGER(8) :: hrs,mins,secs, wall_seconds, count_max
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(1x,a,1p,e10.3,a)
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      IF(mode == 0)THEN
         CALL CPU_TIME(start)
         CALL SYSTEM_CLOCK(wall_start)
      ELSE
         ! report cpu time
         CALL CPU_time(seconds)
         seconds=seconds-start
         secs = int(seconds)
         hrs = secs/(60*60)
         mins = (secs-hrs*60*60)/60
         secs = secs-hrs*60*60-mins*60
         IF(PRESENT(op_cpuseconds))THEN
            ! simply provide the time to the caller
            op_cpuseconds = seconds
         ELSE
            ! write the time to terminal and file
            IF(hrs>0)THEN
               WRITE(*,'(1x,a,i3,a,i2,a,i2,a)') "Total cpu time = ",
     $            hrs," hours, ",mins," minutes, ",secs," seconds"
            ELSEIF(mins>0)THEN
               WRITE(*,'(1x,a,i2,a,i2,a)') "Total cpu time = ",
     $            mins," minutes, ",secs," seconds"
            ELSEIF(secs>0)THEN
               WRITE(*,'(1x,a,i2,a)') "Total cpu time = ",secs,
     $            " seconds"
            ENDIF
            WRITE(unit,10) "Total cpu time = ",seconds," seconds"
         ENDIF
         ! report wall time
         CALL SYSTEM_CLOCK(wall_seconds, count_rate, count_max)
         seconds=(wall_seconds-wall_start)/REAL(count_rate, 8)
         secs = int(seconds)
         hrs = secs/(60*60)
         mins = (secs-hrs*60*60)/60
         secs = secs-hrs*60*60-mins*60
         IF(PRESENT(op_wallseconds))THEN
            op_wallseconds = seconds
         ELSE
            IF(hrs>0)THEN
               WRITE(*,'(1x,a,i3,a,i2,a,i2,a)') "Total wall time = ",
     $            hrs," hours, ",mins," minutes, ",secs," seconds"
            ELSEIF(mins>0)THEN
               WRITE(*,'(1x,a,i2,a,i2,a)') "Total wall time = ",
     $            mins," minutes, ",secs," seconds"
            ELSEIF(secs>0)THEN
               WRITE(*,'(1x,a,i2,a)') "Total wall time = ",secs,
     $            " seconds"
            ENDIF
            WRITE(unit,10) "Total wall time = ",seconds," seconds"
         ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE timer
c-----------------------------------------------------------------------
c     subprogram 2. bin_open.
c     opens a binary input or output file.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE bin_open(unit,name,stat,pos,convert_type)

      CHARACTER(*), INTENT(IN) :: name,stat,pos,convert_type
      INTEGER, INTENT(IN) :: unit
c-----------------------------------------------------------------------
c     open file.
c-----------------------------------------------------------------------
      SELECT CASE(convert_type)
      CASE("none")
         OPEN(UNIT=unit,FILE=name,STATUS=stat,POSITION=pos,
     $        FORM="UNFORMATTED")
      CASE("big")
         OPEN(UNIT=unit,FILE=name,STATUS=stat,POSITION=pos,
     $        FORM="UNFORMATTED",CONVERT="BIG_ENDIAN")
      CASE("little")
         OPEN(UNIT=unit,FILE=name,STATUS=stat,POSITION=pos,
     $        FORM="UNFORMATTED",CONVERT="LITTLE_ENDIAN")
      CASE DEFAULT
         CALL program_stop
     $        ("Cannot recognize convert_type = "//TRIM(convert_type))
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE bin_open
c-----------------------------------------------------------------------
c     subprogram 3. bin_close.
c     close a binary input or output file.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE bin_close(unit)

      INTEGER, INTENT(IN) :: unit
c-----------------------------------------------------------------------
c     work.
c-----------------------------------------------------------------------
      CLOSE(UNIT=unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE bin_close
c-----------------------------------------------------------------------
c     subprogram 4. ascii_open.
c     opens a ascii input or output file.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE ascii_open(unit,name,stat)

      CHARACTER(*), INTENT(IN) :: name,stat
      INTEGER, INTENT(IN) :: unit
c-----------------------------------------------------------------------
c     open file.
c-----------------------------------------------------------------------
      OPEN(UNIT=unit,FILE=name,STATUS=stat)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ascii_open
c-----------------------------------------------------------------------
c     subprogram 5. ascii_close.
c     close a ascii input or output file.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE ascii_close(unit)

      INTEGER, INTENT(IN) :: unit
c-----------------------------------------------------------------------
c     work.
c-----------------------------------------------------------------------
      CLOSE(UNIT=unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ascii_close
c-----------------------------------------------------------------------
c     subprogram 6. program_stop.
c     terminates program with message, calls timer, closes output file.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE program_stop(message)

      CHARACTER(*), INTENT(IN) :: message
c-----------------------------------------------------------------------
c     write completion message.
c-----------------------------------------------------------------------
      CALL timer(1,out_unit)
      CALL ascii_close(out_unit)
      WRITE(*,'(1x,2a)') 'PROGRAM STOP => ', TRIM(message)
c-----------------------------------------------------------------------
c     write completion message.
c-----------------------------------------------------------------------
      STOP
      END SUBROUTINE program_stop

      end module defs_mod