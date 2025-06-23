program test_trace_field
  use iso_fortran_env, only: error_unit
  use m_interp_unstructured

  implicit none
  integer, parameter :: dp = kind(0.0d0)
  type(iu_grid_t)    :: ug
  integer            :: n, i_vx, i_vy

  integer  :: max_steps, n_steps, ndim, nvar
  real(dp) :: r_start(2)
  real(dp) :: min_dx, max_dx
  real(dp) :: rtol, atol
  logical  :: reverse
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
  r_start   = [1.500_dp, 0.0_dp]
  rtol      = 1e-3_dp
  atol      = 1e-3_dp
  min_dx    = 1e-5_dp
  max_dx    = 1.0e-1_dp
  reverse   = .false.
  max_steps = 100
  ndim      = 2
  nvar      = 1

  allocate(y(ndim+nvar, max_steps))
  allocate(y_field(ndim, max_steps))

  call iu_integrate_along_field(ug, ndim, sub_int, r_start, [i_vx, i_vy], &
       min_dx, max_dx, max_steps, rtol, atol, reverse, &
       nvar, y, y_field, n_steps)

  if (n_steps > max_steps) error stop "Boundary not reached"

  print *, "Solution (x, y, length):", y(:, n_steps)
  print *, "n_steps:", n_steps

contains

  subroutine sub_int(ndim, r, field, nvar, integrand)
    integer, intent(in)   :: ndim
    real(dp), intent(in)  :: r(ndim)
    real(dp), intent(in)  :: field(ndim)
    integer, intent(in)   :: nvar
    real(dp), intent(out) :: integrand(nvar)
    integrand(:) = 1.0_dp
  end subroutine sub_int

end program test_trace_field
