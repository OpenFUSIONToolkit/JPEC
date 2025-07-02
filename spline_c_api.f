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
c-----------------------------------------------------------------------
c     subprogram 0. spline_c_api_mod
c     module declarations.
c-----------------------------------------------------------------------
      module spline_c_api_mod
      use iso_c_binding
      use spline_mod
      use cspline_mod
      implicit none
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      type, bind(C) :: spline_handle
         type(c_ptr) :: obj
      end type spline_handle

      logical :: debug = .true.

      contains
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
      integer(i4) :: i
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

      print *, "test"

      select case(endmode)
      case(1)
         endmode_str = "natural"
      case(2)
         endmode_str = "periodic"
      case(3)
         endmode_str = "extrap"
      case(4)
         endmode_str = "not-a-knot"
      end select
      if (debug) then
         print *, "spline_c_fit: fitting spline with endmode = "
     $       // TRIM(endmode_str)
      end if
      call spline_fit(spl, TRIM(endmode_str))
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
      integer(c_int), intent(in), optional :: ix_op

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

      if (present(ix_op)) then
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
      integer(c_int), intent(in), optional :: ix_op

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

      if (present(ix_op)) then
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
      integer(c_int), intent(in), optional :: ix_op

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

      if (present(ix_op)) then
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
      integer(c_int), intent(in), optional :: ix_op

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

      if (present(ix_op)) then
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
      integer(i4) :: i
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

      print *, "test"

      select case(endmode)
      case(1)
         endmode_str = "natural"
      case(2)
         endmode_str = "periodic"
      case(3)
         endmode_str = "extrap"
      case(4)
         endmode_str = "not-a-knot"
      end select
      if (debug) then
         print *, "cspline_c_fit: fitting spline with endmode = "
     $       // TRIM(endmode_str)
      end if
      call cspline_fit(cspl, TRIM(endmode_str))
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
      integer(c_int), intent(in), optional :: ix_op

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

      if (present(ix_op)) then
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
      integer(c_int), intent(in), optional :: ix_op

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

      if (present(ix_op)) then
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
      integer(c_int), intent(in), optional :: ix_op

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

      if (present(ix_op)) then
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
      integer(c_int), intent(in), optional :: ix_op

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

      if (present(ix_op)) then
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
      end module spline_c_api_mod
