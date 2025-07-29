! First export:
!     export FFLAGS='-O0 -fbacktrace -g -fcheck=all -fPIC'
!     export FC=gfortran
! run using: make && gfortran TESTVAC.f *.o -std=legacy -framework Accelerate -o testvac && ./testvac


      MODULE TESTVAC_MOD
      IMPLICIT NONE
      CONTAINS

c-----------------------------------------------------------------------
c     Test 1: Test the green function from vacuum_ma.f
c-----------------------------------------------------------------------
      SUBROUTINE test_green

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
      CALL defglo(mthvac)

      xs = 2.1629189218489446_r8
      zs = 0.49099087247571593_r8
      xt = 1.4309251319336729_r8
      zt = -1.0909471147721690_r8
      xtp = -0.31581804571944494_r8
      ztp = -0.13298631469337879_r8
      n = 1

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
      
      END SUBROUTINE test_green
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c     Test 2: test fourier and fourier inverse from vacuum_vac.f
c-----------------------------------------------------------------------
      SUBROUTINE test_fourier
      USE vglobal_mod
      USE vacuum_mod
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: gij, gil, cs, gll
      REAL(r8), DIMENSION(:), ALLOCATABLE :: x
      INTEGER :: mtheta=200, mthvac=900, mpert=34
      INTEGER :: i, j
      kernelsign=1
      ntsin0=mtheta+1
      nths0=mthvac
      nfm=mpert
      mtot=mpert
      call global_alloc(nths0,nfm,mtot,ntsin0)
      CALL defglo(900)

      lmin(1) = 1
      lmax(1) = nfm

      WRITE(*,*) nths, nths2, nfm, nfm2
      WRITE(*,*) "lmin(1)=", lmin(1), " lmax(1)=", lmax(1)

      nths = 905
      nths2 = 1810
      nfm = 34
      nfm2 = 68

      ALLOCATE(gij(nths,nths), gil(nths2,nfm2), cs(nths,nfm))
      ALLOCATE(gll(nfm,nfm))
      ALLOCATE(x(nths))
      
      gij = 0.0_r8
      gil = 0.0_r8
      cs = 0.0_r8
      gll = 0.0_r8

      DO i=1, nths
         x(i) = 6.28/REAL(i,r8)
      ENDDO


      DO i=1, nths
         DO j=1, nths
            gij(i,j) = SIN(x(i) + x(j))
         ENDDO

         DO j=1, nfm
            cs(i,j) = SIN(j*x(i))**2
         ENDDO
      ENDDO

      CALL fouran(gij,gil,cs,0,34)

      WRITE(*,*) 'gil(1:5,1:5) matrix:'
      DO i = 1, 5
         WRITE(*,'(5(F12.6,2X))') (gil(109+i,33+j), j=1,5)
      END DO

      mth = nths
      dth = twopi / mth

      CALL foranv(gil,gll,cs,0,34)

      WRITE(*,*) 'gll(1:5,1:5) matrix:'
      DO i = 1, 5
            WRITE(*,'(5(F12.6,2X))') (gll(i,j), j=1,5)
      END DO

      END SUBROUTINE test_fourier
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c     Test 3: test mscvac from vacuum_vac.f
c-----------------------------------------------------------------------
      SUBROUTINE test_mscvac
      USE vacuum_mod
      USE vglobal_mod
      IMPLICIT NONE

      !INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15, 307)
      INTEGER, PARAMETER :: mthvac = 900
      INTEGER :: mtheta = 200
      INTEGER :: i, j
      LOGICAL :: complex_flag = .TRUE.
      REAL(r8) :: kernelsignin = 1.0_r8
      LOGICAL :: walflag = .FALSE.
      LOGICAL :: farwalflag = .TRUE.
      INTEGER, PARAMETER :: mpert = 34
      COMPLEX(r8), DIMENSION(mpert, mpert) :: wv
      REAL(r8), DIMENSION(2*(mthvac+5),mpert*2) :: grriio
      REAL(r8), DIMENSION(mthvac+5,4) :: xzptso
      
      CALL mscvac(wv,mpert,mtheta,mthvac,complex_flag,kernelsignin,
     $                   walflag, farwalflag, grriio, xzptso)

      WRITE(*,*) 'wv(1:5,1:5) matrix:'
      DO i = 1, 5
         WRITE(*,'(5("(",F10.6,",",F10.6,")",2X))') (wv(i,j), j=1,5)
      END DO

      END SUBROUTINE test_mscvac
c-----------------------------------------------------------------------

      END MODULE TESTVAC_MOD

c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c     Main program to run the tests
c-----------------------------------------------------------------------
      
      PROGRAM TESTVAC

      USE TESTVAC_MOD

      WRITE(*,*) "Running tests for vacuum module."
      WRITE(*,*) "-----------------------------------"
      WRITE(*,*) "Test 1: Green function"
      CALL test_green
      WRITE(*,*) "-----------------------------------"
      WRITE(*,*) "Test 2: Fourier and Inverse Fourier"
      CALL test_fourier
      WRITE(*,*) "-----------------------------------"
      WRITE(*,*) "Test 3: MSCVAC"
      CALL test_mscvac

      WRITE(*,*) "-----------------------------------"
      WRITE(*,*) "All tests completed successfully."

      END PROGRAM TESTVAC
