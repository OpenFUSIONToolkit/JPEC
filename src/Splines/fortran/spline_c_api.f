c-----------------------------------------------------------------------
c     file spline_c_api.f
c     an interface to the spline libraries for c variables.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. spline_c_api_mod.
c.    1. spline_c_create
c     2. spline_c_destroy
c     3. spline_c_setup
c     4. spline_c_fit
c     5. spline_c_eval
c     6. spline_c_eval_deriv
c     7. spline_c_eval_deriv2
c     8. spline_c_eval_deriv3
c     9. cspline_c_create
c     10. cspline_c_destroy
c     11. cspline_c_setup
c     12. cspline_c_fit
c     13. cspline_c_eval
c     14. cspline_c_eval_deriv
c     15. cspline_c_eval_deriv2
c     16. cspline_c_eval_deriv3
c     17. bicube_c_create
c     18. bicube_c_destroy
c     19. bicube_c_setup
c     20. bicube_c_fit
c     21. bicube_c_eval
c     22. bicube_c_eval_deriv
c     23. bicube_c_eval_deriv2
c-----------------------------------------------------------------------
c     subprogram 0. spline_c_api_mod
c     module declarations.
c-----------------------------------------------------------------------
      module spline_c_api_mod
      use iso_c_binding
      use spline_mod
      use cspline_mod
      use bicube_mod
      implicit none
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      type, bind(C) :: spline_handle
         type(c_ptr) :: obj
      end type spline_handle

      logical :: debug = .false.

      contains

c-----------------------------------------------------------------------
c     Cubic Spline API
c     This section includes the C API for spline.f.
c-----------------------------------------------------------------------



c-----------------------------------------------------------------------
c     subprogram 1. spline_c_create
c     allocates a spline object
c-----------------------------------------------------------------------
      subroutine spline_c_create(mx, nqty, handle) bind(C)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      integer(c_int), value :: mx, nqty
      type(spline_handle), intent(out) :: handle
      type(spline_type), pointer :: spl

      allocate(spl)
      call spline_alloc(spl, mx, nqty)
      handle%obj = c_loc(spl)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      return
      end subroutine spline_c_create
c-----------------------------------------------------------------------
c     subprogram 2. spline_c_destroy
c     deallocates a spline object
c-----------------------------------------------------------------------
      subroutine spline_c_destroy(handle) bind(C)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      type(spline_handle), value :: handle
      type(spline_type), pointer :: spl
      
      call c_f_pointer(handle%obj, spl)
      if (associated(spl)) then
         call spline_dealloc(spl)
         deallocate(spl)
      else
         print *, "spline_c_destroy: handle is not associated "
     $       // "with a valid spline object." 
      end if
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      return
      end subroutine spline_c_destroy
c-----------------------------------------------------------------------
c     subprogram 3. spline_c_setup
c     sets up the spline object with data.
c-----------------------------------------------------------------------
      subroutine spline_c_setup(handle, xs, fs) bind(C)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      type(spline_handle), value :: handle
      type(c_ptr), value :: xs, fs

      type(spline_type), pointer :: spl
      real(c_double), pointer :: x(:), f(:, :)
      integer(i8) :: mx, nqty
      integer(i8) :: i
c-----------------------------------------------------------------------
c     work.
c-----------------------------------------------------------------------
      call c_f_pointer(handle%obj, spl)
      if (.not. associated(spl)) then
         print *, "spline_c_setup: handle is not associated "
     $       // "with a valid spline object."
         return
      end if

      mx = spl%mx
      nqty = spl%nqty

      call c_f_pointer(xs, x, [mx+1])
      call c_f_pointer(fs, f, [mx+1, nqty])

      ! do i = 0, mx
      !    spl%xs(i) = x(i+1)  ! Fortran is 1-based, C is 0-based
      !    spl%fs(i, 1:nqty) = f(i+1, 1:nqty)
      ! end do

      spl%xs = x
      spl%fs = f

      if (debug) then
         print *, "spline_c_setup: setting up spline with "
     $       // "mx = ", mx, " and nqty = ", nqty
         print *, "xs = ", x(1:mx+1)
         print *, "fs = "
     $       // "(", nqty, " quantities):"
         do i = 1, nqty
            print *, "  fs(:,", i, ") = ", f(:, i)
         end do
      end if
c------------------------------------------------------------------------
c     terminate.
c------------------------------------------------------------------------
      return 
      end subroutine spline_c_setup
c-----------------------------------------------------------------------
c     subprogram 4. spline_c_fit
c     fits the spline to the data.
c-----------------------------------------------------------------------
      subroutine spline_c_fit(handle, endmode) bind(C)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      type(spline_handle), value :: handle
      integer(c_int), value :: endmode
      type(spline_type), pointer :: spl

      character(10) :: endmode_str
c-----------------------------------------------------------------------
c     work.
c-----------------------------------------------------------------------
      call c_f_pointer(handle%obj, spl)
      if (.not. associated(spl)) then
            print *, "spline_c_fit: handle is not associated "
     $       // "with a valid spline object."
            return
      end if

c      select case(endmode)
c      case(1)
c         endmode_str = "natural"
c      case(2)
c         endmode_str = "periodic"
c      case(3)
c         endmode_str = "extrap"
c      case(4)
c         endmode_str = "not-a-knot"
c      end select
      
c      if (debug) then
c         print *, "spline_c_fit: fitting spline with endmode = "
c     $       // TRIM(endmode_str)
c      end if
      call spline_fit(spl, endmode)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      return
      end subroutine spline_c_fit
c-----------------------------------------------------------------------
c     subprogram 5. spline_c_eval
c     evaluates the spline at a given point.
c-----------------------------------------------------------------------
      subroutine spline_c_eval(handle, x, f, ix_op) bind(C)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      type(spline_handle), value :: handle
      real(c_double), value :: x
      type(c_ptr), value :: f
      integer(c_int), intent(in), value :: ix_op

      integer(i4) :: ix ! index of x position in the spline
      real(c_double), pointer :: fi(:)
      type(spline_type), pointer :: spl
c-----------------------------------------------------------------------
c     work.
c-----------------------------------------------------------------------
      call c_f_pointer(handle%obj, spl)
      if (.not. associated(spl)) then
         print *, "spline_c_eval: handle is not associated "
     $       // "with a valid spline object."
         return
      end if

      call c_f_pointer(f, fi, [spl%nqty])

      if (ix_op>0) then
         ix = ix_op
      else
         ix = 0
      end if

      call spline_eval_external(spl, x, ix, fi)
c-----------------------------------------------------------------------
c     copy results back to the C pointer.
c-----------------------------------------------------------------------
      call c_f_pointer(f, fi, [spl%nqty])
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      return
      end subroutine spline_c_eval
c-----------------------------------------------------------------------
c     subprogram 6. spline_c_eval_deriv
c     evaluates the spline and its first derivative at a given point.
c-----------------------------------------------------------------------
      subroutine spline_c_eval_deriv(handle, x, f, f1, ix_op) bind(C)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      type(spline_handle), value :: handle
      real(c_double), value :: x
      type(c_ptr), value :: f, f1
      integer(c_int), intent(in), value :: ix_op

      integer(i4) :: ix ! index of x position in the spline
      real(c_double), pointer :: fi(:), f1i(:)
      type(spline_type), pointer :: spl
c-----------------------------------------------------------------------
c     work.
c-----------------------------------------------------------------------
      call c_f_pointer(handle%obj, spl)
      if (.not. associated(spl)) then
         print *, "spline_c_eval: handle is not associated "
     $       // "with a valid spline object."
         return
      end if

      call c_f_pointer(f, fi, [spl%nqty])
      call c_f_pointer(f1, f1i, [spl%nqty])

      if (ix_op>0) then
         ix = ix_op
      else
         ix = 0
      end if

      call spline_eval_external(spl, x, ix, fi, f1i)
c-----------------------------------------------------------------------
c     copy results back to the C pointer.
c-----------------------------------------------------------------------
      call c_f_pointer(f, fi, [spl%nqty])
      call c_f_pointer(f1, f1i, [spl%nqty])

c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      return
      end subroutine spline_c_eval_deriv
c-----------------------------------------------------------------------
c     subprogram 7. spline_c_eval_deriv_2
c     evaluates the spline and its first and\
c     second derivatives at a given point.
c-----------------------------------------------------------------------
      subroutine spline_c_eval_deriv_2(handle, x, f, f1,
     $     f2, ix_op) bind(C)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      type(spline_handle), value :: handle
      real(c_double), value :: x
      type(c_ptr), value :: f, f1, f2
      integer(c_int), intent(in), value :: ix_op

      integer(i4) :: ix ! index of x position in the spline
      real(c_double), pointer :: fi(:), f1i(:), f2i(:)
      type(spline_type), pointer :: spl
c-----------------------------------------------------------------------
c     work.
c-----------------------------------------------------------------------
      call c_f_pointer(handle%obj, spl)
      if (.not. associated(spl)) then
         print *, "spline_c_eval: handle is not associated "
     $       // "with a valid spline object."
         return
      end if

      call c_f_pointer(f, fi, [spl%nqty])
      call c_f_pointer(f1, f1i, [spl%nqty])
      call c_f_pointer(f2, f2i, [spl%nqty])

      if (ix_op>0) then
         ix = ix_op
      else
         ix = 0
      end if

      call spline_eval_external(spl, x, ix, fi, f1i, f2i)
c-----------------------------------------------------------------------
c     copy results back to the C pointer.
c-----------------------------------------------------------------------
      call c_f_pointer(f, fi, [spl%nqty])
      call c_f_pointer(f1, f1i, [spl%nqty])
      call c_f_pointer(f2, f2i, [spl%nqty])

c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      return
      end subroutine spline_c_eval_deriv_2
c-----------------------------------------------------------------------
c     subprogram 8. spline_c_eval_deriv_3
c     evaluates the spline and its first derivative at a given point.
c-----------------------------------------------------------------------
      subroutine spline_c_eval_deriv_3(handle, x, f, f1,
     $    f2, f3, ix_op) bind(C)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      type(spline_handle), value :: handle
      real(c_double), value :: x
      type(c_ptr), value :: f, f1, f2, f3
      integer(c_int), intent(in), value :: ix_op

      integer(i4) :: ix ! index of x position in the spline
      real(c_double), pointer :: fi(:), f1i(:), f2i(:), f3i(:)
      type(spline_type), pointer :: spl
c-----------------------------------------------------------------------
c     work.
c-----------------------------------------------------------------------
      call c_f_pointer(handle%obj, spl)
      if (.not. associated(spl)) then
         print *, "spline_c_eval: handle is not associated "
     $       // "with a valid spline object."
         return
      end if

      call c_f_pointer(f, fi, [spl%nqty])
      call c_f_pointer(f1, f1i, [spl%nqty])
      call c_f_pointer(f2, f2i, [spl%nqty])
      call c_f_pointer(f3, f3i, [spl%nqty])

      if (ix_op>0) then
         ix = ix_op
      else
         ix = 0
      end if

      call spline_eval_external(spl, x, ix, fi, f1i, f2i, f3i)
c-----------------------------------------------------------------------
c     copy results back to the C pointer.
c-----------------------------------------------------------------------
      call c_f_pointer(f, fi, [spl%nqty])
      call c_f_pointer(f1, f1i, [spl%nqty])
      call c_f_pointer(f2, f2i, [spl%nqty])
      call c_f_pointer(f3, f3i, [spl%nqty])

c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      return
      end subroutine spline_c_eval_deriv_3


c-----------------------------------------------------------------------
c     Complex Cubic Spline API
c     This section includes the C API for cspline.f.
c-----------------------------------------------------------------------


c-----------------------------------------------------------------------
c     subprogram 9. cspline_c_create
c     allocates a spline object
c-----------------------------------------------------------------------
      subroutine cspline_c_create(mx, nqty, handle) bind(C)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      integer(c_int), value :: mx, nqty
      type(spline_handle), intent(out) :: handle
      type(cspline_type), pointer :: cspl

      allocate(cspl)
      call cspline_alloc(cspl, mx, nqty)
      handle%obj = c_loc(cspl)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      return
      end subroutine cspline_c_create
c-----------------------------------------------------------------------
c     subprogram 10. cspline_c_destroy
c     deallocates a spline object
c-----------------------------------------------------------------------
      subroutine cspline_c_destroy(handle) bind(C)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      type(spline_handle), value :: handle
      type(cspline_type), pointer :: cspl

      call c_f_pointer(handle%obj, cspl)
      if (associated(cspl)) then
         call cspline_dealloc(cspl)
         deallocate(cspl)
      else
         print *, "cspline_c_destroy: handle is not associated "
     $       // "with a valid cspline object." 
      end if
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      return
      end subroutine cspline_c_destroy
c-----------------------------------------------------------------------
c     subprogram 11. cspline_c_setup
c     sets up the spline object with data.
c-----------------------------------------------------------------------
      subroutine cspline_c_setup(handle, xs, fs) bind(C)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      type(spline_handle), value :: handle
      type(c_ptr), value :: xs, fs

      type(cspline_type), pointer :: cspl
      real(c_double), pointer :: x(:)
      complex(c_double), pointer :: f(:, :)
      integer(i8) :: mx, nqty
      integer(i8) :: i
c-----------------------------------------------------------------------
c     work.
c-----------------------------------------------------------------------
      call c_f_pointer(handle%obj, cspl)
      if (.not. associated(cspl)) then
         print *, "cspline_c_setup: handle is not associated "
     $       // "with a valid cspline object."
         return
      end if

      mx = cspl%mx
      nqty = cspl%nqty

      call c_f_pointer(xs, x, [mx+1])
      call c_f_pointer(fs, f, [mx+1, nqty])

      cspl%xs = x
      cspl%fs = f

      if (debug) then
         print *, "cspline_c_setup: setting up spline with "
     $       // "mx = ", mx, " and nqty = ", nqty
         print *, "xs = ", x(1:mx+1)
         print *, "fs = "
     $       // "(", nqty, " quantities):"
         do i = 1, nqty
            print *, "  fs(:,", i, ") = ", f(:, i)
         end do
      end if
c------------------------------------------------------------------------
c     terminate.
c------------------------------------------------------------------------
      return 
      end subroutine cspline_c_setup
c-----------------------------------------------------------------------
c     subprogram 12. spline_c_fit
c     fits the spline to the data.
c-----------------------------------------------------------------------
      subroutine cspline_c_fit(handle, endmode) bind(C)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      type(spline_handle), value :: handle
      integer(c_int), value :: endmode
      type(cspline_type), pointer :: cspl

      character(10) :: endmode_str
c-----------------------------------------------------------------------
c     work.
c-----------------------------------------------------------------------
      
      call c_f_pointer(handle%obj, cspl)
      if (.not. associated(cspl)) then
            print *, "cspline_c_fit: handle is not associated "
     $       // "with a valid cspline object."
            return
      end if

c      select case(endmode)
c      case(1)
c         endmode_str = "natural"
c      case(2)
c         endmode_str = "periodic"
c      case(3)
c         endmode_str = "extrap"
c      case(4)
c         endmode_str = "not-a-knot"
c      end select
c      if (debug) then
c         print *, "cspline_c_fit: fitting spline with endmode = "
c     $       // TRIM(endmode_str)
c      end if
      call cspline_fit(cspl, endmode)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      return
      end subroutine cspline_c_fit
c-----------------------------------------------------------------------
c     subprogram 13. cspline_c_eval
c     evaluates the spline at a given point.
c-----------------------------------------------------------------------
      subroutine cspline_c_eval(handle, x, f, ix_op) bind(C)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      type(spline_handle), value :: handle
      real(c_double), value :: x
      type(c_ptr), value :: f
      integer(c_int), intent(in), value :: ix_op

      integer(i4) :: ix ! index of x position in the spline
      complex(c_double), pointer :: fi(:)
      type(cspline_type), pointer :: cspl
c-----------------------------------------------------------------------
c     work.
c-----------------------------------------------------------------------
      call c_f_pointer(handle%obj, cspl)
      if (.not. associated(cspl)) then
         print *, "cspline_c_eval: handle is not associated "
     $       // "with a valid cspline object."
         return
      end if

      call c_f_pointer(f, fi, [cspl%nqty])

      if (ix_op>0) then
         ix = ix_op
      else
         ix = 0
      end if

      call cspline_eval_external(cspl, x, ix, fi)
c-----------------------------------------------------------------------
c     copy results back to the C pointer.
c-----------------------------------------------------------------------
      call c_f_pointer(f, fi, [cspl%nqty])
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      return
      end subroutine cspline_c_eval
c-----------------------------------------------------------------------
c     subprogram 14. cspline_c_eval_deriv
c     evaluates the spline and its first derivative at a given point.
c-----------------------------------------------------------------------
      subroutine cspline_c_eval_deriv(handle, x, f, f1, ix_op) bind(C)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      type(spline_handle), value :: handle
      real(c_double), value :: x
      type(c_ptr), value :: f, f1
      integer(c_int), intent(in), value :: ix_op

      integer(i4) :: ix ! index of x position in the spline
      complex(c_double), pointer :: fi(:), f1i(:)
      type(cspline_type), pointer :: cspl
c-----------------------------------------------------------------------
c     work.
c-----------------------------------------------------------------------
      call c_f_pointer(handle%obj, cspl)
      if (.not. associated(cspl)) then
         print *, "cspline_c_eval: handle is not associated "
     $       // "with a valid cspline object."
         return
      end if

      call c_f_pointer(f, fi, [cspl%nqty])
      call c_f_pointer(f1, f1i, [cspl%nqty])

      if (ix_op>0) then
         ix = ix_op
      else
         ix = 0
      end if

      call cspline_eval_external(cspl, x, ix, fi, f1i)
c-----------------------------------------------------------------------
c     copy results back to the C pointer.
c-----------------------------------------------------------------------
      call c_f_pointer(f, fi, [cspl%nqty])
      call c_f_pointer(f1, f1i, [cspl%nqty])

c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      return
      end subroutine cspline_c_eval_deriv
c-----------------------------------------------------------------------
c     subprogram 15. cspline_c_eval_deriv_2
c     evaluates the spline and its first and\
c     second derivatives at a given point.
c-----------------------------------------------------------------------
      subroutine cspline_c_eval_deriv_2(handle, x, f, f1,
     $     f2, ix_op) bind(C)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      type(spline_handle), value :: handle
      real(c_double), value :: x
      type(c_ptr), value :: f, f1, f2
      integer(c_int), intent(in), value :: ix_op

      integer(i4) :: ix ! index of x position in the spline
      complex(c_double), pointer :: fi(:), f1i(:), f2i(:)
      type(cspline_type), pointer :: cspl
c-----------------------------------------------------------------------
c     work.
c-----------------------------------------------------------------------
      call c_f_pointer(handle%obj, cspl)
      if (.not. associated(cspl)) then
         print *, "cspline_c_eval: handle is not associated "
     $       // "with a valid cspline object."
         return
      end if

      call c_f_pointer(f, fi, [cspl%nqty])
      call c_f_pointer(f1, f1i, [cspl%nqty])
      call c_f_pointer(f2, f2i, [cspl%nqty])

      if (ix_op>0) then
         ix = ix_op
      else
         ix = 0
      end if

      call cspline_eval_external(cspl, x, ix, fi, f1i, f2i)
c-----------------------------------------------------------------------
c     copy results back to the C pointer.
c-----------------------------------------------------------------------
      call c_f_pointer(f, fi, [cspl%nqty])
      call c_f_pointer(f1, f1i, [cspl%nqty])
      call c_f_pointer(f2, f2i, [cspl%nqty])

c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      return
      end subroutine cspline_c_eval_deriv_2
c-----------------------------------------------------------------------
c     subprogram 16. cspline_c_eval_deriv_3
c     evaluates the spline and its first derivative at a given point.
c-----------------------------------------------------------------------
      subroutine cspline_c_eval_deriv_3(handle, x, f, f1,
     $    f2, f3, ix_op) bind(C)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      type(spline_handle), value :: handle
      real(c_double), value :: x
      type(c_ptr), value :: f, f1, f2, f3
      integer(c_int), intent(in), value :: ix_op

      integer(i4) :: ix ! index of x position in the spline
      complex(c_double), pointer :: fi(:), f1i(:), f2i(:), f3i(:)
      type(cspline_type), pointer :: cspl
c-----------------------------------------------------------------------
c     work.
c-----------------------------------------------------------------------
      call c_f_pointer(handle%obj, cspl)
      if (.not. associated(cspl)) then
         print *, "cspline_c_eval: handle is not associated "
     $       // "with a valid cspline object."
         return
      end if

      call c_f_pointer(f, fi, [cspl%nqty])
      call c_f_pointer(f1, f1i, [cspl%nqty])
      call c_f_pointer(f2, f2i, [cspl%nqty])
      call c_f_pointer(f3, f3i, [cspl%nqty])

      if (ix_op>0) then
         ix = ix_op
      else
         ix = 0
      end if

      call cspline_eval_external(cspl, x, ix, fi, f1i, f2i, f3i)
c-----------------------------------------------------------------------
c     copy results back to the C pointer.
c-----------------------------------------------------------------------
      call c_f_pointer(f, fi, [cspl%nqty])
      call c_f_pointer(f1, f1i, [cspl%nqty])
      call c_f_pointer(f2, f2i, [cspl%nqty])
      call c_f_pointer(f3, f3i, [cspl%nqty])

c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      return
      end subroutine cspline_c_eval_deriv_3

c-----------------------------------------------------------------------
c     Bicubic Spline API
c     This section includes the C API for bicube.f.
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c     subprogram 17. bicube_c_create
c     allocates a bicubic spline object
c-----------------------------------------------------------------------
      subroutine bicube_c_create(mx, my, nqty, handle) bind(C)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      integer(c_int), value :: mx, my, nqty
      type(spline_handle), intent(out) :: handle
      type(bicube_type), pointer :: bicube

      allocate(bicube)
      call bicube_alloc(bicube, mx, my, nqty)
      handle%obj = c_loc(bicube)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      return
      end subroutine bicube_c_create
c-----------------------------------------------------------------------
c     subprogram 18. bicube_c_destroy
c     deallocates a bicubic spline object
c-----------------------------------------------------------------------
      subroutine bicube_c_destroy(handle) bind(C)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      type(spline_handle), value :: handle
      type(bicube_type), pointer :: bicube

      call c_f_pointer(handle%obj, bicube)
      if (associated(bicube)) then
         call bicube_dealloc(bicube)
         deallocate(bicube)
      else
         print *, "bicube_c_destroy: handle is not associated "
     $       // "with a valid bicubic spline object."
      end if
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      return
      end subroutine bicube_c_destroy
c-----------------------------------------------------------------------
c     subprogram 19. bicube_c_setup
c     sets up the bicubic spline object with data.
c-----------------------------------------------------------------------
      subroutine bicube_c_setup(handle, xs, ys, fs) bind(C)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      type(spline_handle), value :: handle
      type(c_ptr), value :: xs, ys, fs

      type(bicube_type), pointer :: bicube
      real(c_double), pointer :: x(:), y(:), f(:, :, :)
      integer(i8) :: mx, my, nqty
      integer(i8) :: i
c-----------------------------------------------------------------------
c     work.
c-----------------------------------------------------------------------
      call c_f_pointer(handle%obj, bicube)
      if (.not. associated(bicube)) then
         print *, "bicube_c_setup: handle is not associated "
     $       // "with a valid bicubic spline object."
         return
      end if

      mx = bicube%mx
      my = bicube%my
      nqty = bicube%nqty

      call c_f_pointer(xs, x, [mx+1])
      call c_f_pointer(ys, y, [my+1])
      call c_f_pointer(fs, f, [mx+1, my+1, nqty])

      ! do i = 0, mx
      !    spl%xs(i) = x(i+1)  ! Fortran is 1-based, C is 0-based
      !    spl%fs(i, 1:nqty) = f(i+1, 1:nqty)
      ! end do

      bicube%xs = x
      bicube%ys = y
      bicube%fs = f

      if (debug) then
         print *, "bicube_c_setup: setting up bicubic spline with "
     $       // "mx = ", mx, ", my = ", my, " and nqty = ", nqty
         print *, "xs = ", x(1:mx+1)
         print *, "ys = ", y(1:my+1)
         print *, "fs = "
     $       // "(", nqty, " quantities):"
         do i = 1, nqty
            print *, "  fs(:,:,", i, ") = ", f(:,:,i)
         end do
      end if
c------------------------------------------------------------------------
c     terminate.
c------------------------------------------------------------------------
      return
      end subroutine bicube_c_setup
c-----------------------------------------------------------------------
c     subprogram 20. bicube_c_fit
c     fits the bicubic spline to the data.
c-----------------------------------------------------------------------
      subroutine bicube_c_fit(handle, endmode1, endmode2) bind(C)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      type(spline_handle), value :: handle
      integer(c_int), value :: endmode1, endmode2
      type(bicube_type), pointer :: bicube

      character(10) :: endmode1_str, endmode2_str
c-----------------------------------------------------------------------
c     work.
c-----------------------------------------------------------------------
      call c_f_pointer(handle%obj, bicube)
      if (.not. associated(bicube)) then
            print *, "bicube_c_fit: handle is not associated "
     $       // "with a valid bicubic spline object."
            return
      end if

c      select case(endmode1)
c      case(1)
c         endmode1_str = "natural"
c      case(2)
c         endmode1_str = "periodic"
c      case(3)
c         endmode1_str = "extrap"
c      case(4)
c         endmode1_str = "not-a-knot"
c      end select
c      select case(endmode2)
c      case(1)
c         endmode2_str = "natural"
c      case(2)
c         endmode2_str = "periodic"
c      case(3)
c         endmode2_str = "extrap"
c      case(4)
c         endmode2_str = "not-a-knot"
c      end select

c      if (debug) then
c         print *, "bicube_c_fit: fitting spline with endmode1 = "
c     $       // TRIM(endmode1_str) // " and endmode2 = "
c     $       // TRIM(endmode2_str)
c      end if
      call bicube_fit(bicube, endmode1, endmode2)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      return
      end subroutine bicube_c_fit
c-----------------------------------------------------------------------
c     subprogram 21. bicube_c_eval
c     evaluates the bicubic spline at a given point.
c-----------------------------------------------------------------------
      subroutine bicube_c_eval(handle, x, y, f, ix_op, iy_op) bind(C)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      type(spline_handle), value :: handle
      real(c_double), value :: x, y
      type(c_ptr), value :: f
      integer(c_int), intent(in), value :: ix_op, iy_op

      integer(i4) :: ix ! index of x position in the spline
      integer(i4) :: iy ! index of y position in the spline
      real(c_double), pointer :: fi(:)
      real(r8), pointer :: fix(:), fiy(:), fxy(:), fxx(:), fyy(:)
      type(bicube_type), pointer :: bicube
c-----------------------------------------------------------------------
c     work.
c-----------------------------------------------------------------------
      call c_f_pointer(handle%obj, bicube)
      if (.not. associated(bicube)) then
         print *, "bicube_c_eval: handle is not associated "
     $       // "with a valid bicubic spline object."
         return
      end if

      call c_f_pointer(f, fi, [bicube%nqty])
      allocate(fix(bicube%nqty), fiy(bicube%nqty))
      
      if (ix_op>0) then
         ix = ix_op
      else
         ix = 0
      end if
      
      if (iy_op>0) then
         iy = iy_op
      else
         iy = 0
      end if

      call bicube_eval_external(bicube, x, y, 0, ix, iy, fi, fix, fiy,
     $                                                    fxy, fxx, fyy)
      deallocate(fix, fiy)

c-----------------------------------------------------------------------
c     copy results back to the C pointer.
c-----------------------------------------------------------------------
      call c_f_pointer(f, fi, [bicube%nqty])
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      return
      end subroutine bicube_c_eval
c-----------------------------------------------------------------------
c     subprogram 22. bicube_c_eval_deriv
c     evaluates the bicube and its first derivative at a given point.
c-----------------------------------------------------------------------
      subroutine bicube_c_eval_deriv(handle, x, y,
     $                                f, f1x, f1y, ix_op, iy_op) bind(C)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      type(spline_handle), value :: handle
      real(c_double), value :: x, y
      type(c_ptr), value :: f, f1x, f1y
      integer(c_int), intent(in), value :: ix_op, iy_op

      integer(i4) :: ix ! index of x position in the spline
      integer(i4) :: iy ! index of y position in the spline
      real(c_double), pointer :: fi(:), f1xi(:), f1yi(:)
      real(r8), pointer :: fxy(:), fxx(:), fyy(:)
      type(bicube_type), pointer :: bicube
c-----------------------------------------------------------------------
c     work.
c-----------------------------------------------------------------------
      call c_f_pointer(handle%obj, bicube)
      if (.not. associated(bicube)) then
         print *, "bicube_c_eval_deriv: handle is not associated "
     $       // "with a valid bicubic spline object."
         return
      end if

      call c_f_pointer(f, fi, [bicube%nqty])
      call c_f_pointer(f1x, f1xi, [bicube%nqty])
      call c_f_pointer(f1y, f1yi, [bicube%nqty])
      allocate(fxy(bicube%nqty), fxx(bicube%nqty), fyy(bicube%nqty))

      if (ix_op>0) then
         ix = ix_op
      else
         ix = 0
      end if

      if (iy_op>0) then
         iy = iy_op
      else
         iy = 0
      end if

      call bicube_eval_external(bicube, x, y, 1, ix, iy, fi, f1xi, f1yi,
     $                                                    fxy, fxx, fyy)
      deallocate(fxy, fxx, fyy)
c-----------------------------------------------------------------------
c     copy results back to the C pointer.
c-----------------------------------------------------------------------
      call c_f_pointer(f, fi, [bicube%nqty])
      call c_f_pointer(f1x, f1xi, [bicube%nqty])
      call c_f_pointer(f1y, f1yi, [bicube%nqty])

c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      return
      end subroutine bicube_c_eval_deriv
c-----------------------------------------------------------------------
c     subprogram 23. bicube_c_eval_deriv2
c     evaluates the bicube and its derivatives at a given point.
c-----------------------------------------------------------------------
      subroutine bicube_c_eval_deriv2(handle, x, y, f, f1x, f1y,
     $                           f2xx, f2xy, f2yy, ix_op, iy_op) bind(C)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      type(spline_handle), value :: handle
      real(c_double), value :: x, y
      type(c_ptr), value :: f, f1x, f1y, f2xx, f2xy, f2yy
      integer(c_int), intent(in), value :: ix_op, iy_op

      integer(i4) :: ix ! index of x position in the spline
      integer(i4) :: iy ! index of y position in the spline
      real(c_double), pointer :: fi(:), f1xi(:), f1yi(:)
      real(c_double), pointer :: f2xxi(:), f2xyi(:), f2yyi(:)
      type(bicube_type), pointer :: bicube
c-----------------------------------------------------------------------
c     work.
c-----------------------------------------------------------------------
      call c_f_pointer(handle%obj, bicube)
      if (.not. associated(bicube)) then
         print *, "bicube_c_eval_deriv2: handle is not associated "
     $       // "with a valid bicubic spline object."
         return
      end if

      call c_f_pointer(f, fi, [bicube%nqty])
      call c_f_pointer(f1x, f1xi, [bicube%nqty])
      call c_f_pointer(f1y, f1yi, [bicube%nqty])
      call c_f_pointer(f2xx, f2xxi, [bicube%nqty])
      call c_f_pointer(f2xy, f2xyi, [bicube%nqty])
      call c_f_pointer(f2yy, f2yyi, [bicube%nqty])

      if (ix_op>0) then
         ix = ix_op
      else
         ix = 0
      end if

      if (iy_op>0) then
         iy = iy_op
      else
         iy = 0
      end if

      call bicube_eval_external(bicube, x, y, 2, ix, iy, fi, f1xi, f1yi,
     $                                              f2xxi, f2xyi, f2yyi)
c-----------------------------------------------------------------------
c     copy results back to the C pointer.
c-----------------------------------------------------------------------
      call c_f_pointer(f, fi, [bicube%nqty])
      call c_f_pointer(f1x, f1xi, [bicube%nqty])
      call c_f_pointer(f1y, f1yi, [bicube%nqty])
      call c_f_pointer(f2xx, f2xxi, [bicube%nqty])
      call c_f_pointer(f2xy, f2xyi, [bicube%nqty])
      call c_f_pointer(f2yy, f2yyi, [bicube%nqty])

c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      return
      end subroutine bicube_c_eval_deriv2

c-----------------------------------------------------------------------
      end module spline_c_api_mod
