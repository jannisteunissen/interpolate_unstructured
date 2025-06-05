program test_trace_field
  use iso_fortran_env, only: error_unit
  use m_interp_unstructured

  implicit none
  integer, parameter :: dp = kind(0.0d0)
  type(iu_grid_t)    :: ug
  integer            :: n, i_vx, i_vy

  integer  :: max_steps, n_steps
  real(dp) :: r_start(2)
  real(dp) :: min_dx, max_dx
  real(dp) :: rtol, atol, y(3)
  logical  :: reverse

  call iu_read_grid('test_data/triangle.binda', ug)
  call iu_add_point_data(ug, 'vx', i_vx)
  call iu_add_point_data(ug, 'vy', i_vy)

  do n = 1, ug%n_points
     ug%point_data(n, i_vx) = -ug%points(2, n)
     ug%point_data(n, i_vy) = ug%points(1, n)
  end do

  call iu_write_vtk(ug, 'test_data/test_trace_field.vtu')

  r_start     = [1.500_dp, 0.0_dp]
  rtol        = 1e-3_dp
  atol        = 1e-3_dp
  min_dx      = 1e-5_dp
  max_dx      = 1.0e-1_dp
  reverse     = .false.
  max_steps   = 100

  call iu_integrate_along_field(ug, dummy_func, 2, r_start, [i_vx, i_vy], &
       min_dx, max_dx, max_steps, rtol, atol, reverse, y, n_steps)

  print *, "Solution (x, y, length):", y
  print *, "n_steps:", n_steps

contains

  real(dp) function dummy_func(ndim, r, field)
    integer, intent(in) :: ndim
    real(dp), intent(in) :: r(ndim)
    real(dp), intent(in) :: field(ndim)
    dummy_func = 1.0_dp
  end function dummy_func

end program test_trace_field
