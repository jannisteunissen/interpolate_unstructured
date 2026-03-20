program test_trace_field
  use iso_fortran_env, only: error_unit
  use m_interp_unstructured

  implicit none
  integer, parameter :: dp = kind(0.0d0)
  type(iu_grid_t)    :: ug
  integer            :: n, i_vx, i_vy

  integer  :: max_steps, n_steps, ndim, nvar
  real(dp) :: min_dx, max_dx
  real(dp) :: rtol, atol
  logical  :: reverse, axisymmetric
  real(dp), allocatable ::  y(:, :), y_field(:, :)

  call iu_read_grid('test_data/triangle.binda', ug)
  call iu_add_point_data(ug, 'vx', i_vx)
  call iu_add_point_data(ug, 'vy', i_vy)

  do n = 1, ug%n_points
     ug%point_data(n, i_vx) = -ug%points(2, n)
     ug%point_data(n, i_vy) = ug%points(1, n)
  end do

  call iu_write_vtk(ug, 'test_data/test_trace_field.vtu')

  ndim      = 2
  rtol      = 1e-3_dp
  atol      = 1e-3_dp
  min_dx    = 1e-5_dp
  max_dx    = 1.0e-1_dp
  reverse   = .false.
  axisymmetric = .false.
  max_steps = 100
  ndim      = 2
  nvar      = 1

  allocate(y(ndim+nvar, max_steps))
  allocate(y_field(ndim, max_steps))

  y(1:ndim, 1) = [1.500_dp, 0.0_dp]

  ! 1.5 * pi / 2, so that final solution is about zero
  y(ndim+1:, 1) = -0.75_dp * acos(-1.0_dp)

  call iu_integrate_along_field(ug, ndim, nvar, sub_int, [i_vx, i_vy], &
       min_dx, max_dx, max_steps, rtol, atol, reverse, &
       y, y_field, n_steps, axisymmetric)

  if (n_steps > max_steps) error stop "Boundary not reached"

  print *, "Solution (x, y, length):", y(:, n_steps)
  print *, "n_steps:", n_steps

contains

  subroutine sub_int(ndim, nvar, field, y, dy_var)
    integer, intent(in)   :: ndim
    integer, intent(in)   :: nvar
    real(dp), intent(in)  :: field(ndim)
    real(dp), intent(in)  :: y(ndim+nvar)
    real(dp), intent(out) :: dy_var(nvar)
    dy_var(:) = 1.0_dp
  end subroutine sub_int

end program test_trace_field
