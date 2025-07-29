      PROGRAM TESTVAC
      USE vacuum_mod
      USE vglobal_mod
      IMPLICIT NONE

      !INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15, 307)
      INTEGER, PARAMETER :: mthvac = 900
      INTEGER :: mtheta = 200
      LOGICAL :: complex_flag = .TRUE.
      REAL(r8) :: kernelsignin = 1.0_r8
      LOGICAL :: walflag = .FALSE.
      LOGICAL :: farwalflag = .TRUE.
      INTEGER, PARAMETER :: mpert = 20


      COMPLEX(r8), DIMENSION(mpert, mpert) :: wv
      REAL(r8), DIMENSION(2*(mthvac+5),mpert*2) :: grriio
      REAL(r8), DIMENSION(mthvac+5,4) :: xzptso
      
      ! REAL(r8) :: xs, zs, xt, zt, xtp, ztp
      ! INTEGER :: n

      xs = 2.1629189218489446_r8
      zs = 0.49099087247571593_r8
      xt = 1.4309251319336729_r8
      zt = -1.0909471147721690_r8
      xtp = -0.31581804571944494_r8
      ztp = -0.13298631469337879_r8
      n = 1


      CALL defglo(mthvac)
      CALL green
!       CALL mscvac(wv,mpert,mtheta,mthvac,complex_flag,kernelsignin,
!      $                   walflag, farwalflag, grriio, xzptso)

c-----------------------------------------------------------------------
c     Write a bunch of outputs:
c     bval  
c     aval
c     aval0 
c-----------------------------------------------------------------------

      WRITE(*,*) "bval = ", bval
      WRITE(*,*) "aval = ", aval
      WRITE(*,*) "aval0 = ", aval0


      END PROGRAM TESTVAC