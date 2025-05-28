program test_trace_field
  use iso_fortran_env, only: error_unit
  use m_interp_unstructured

  implicit none
  integer, parameter :: dp = kind(0.0d0)
  type(iu_grid_t)    :: ug
  integer            :: n, i_vx, i_vy

  integer, parameter :: max_points = 1000
  integer            :: n_points
  real(dp)           :: r_start(2)
  real(dp)           :: points(2, max_points)
  real(dp)           :: fields(2, max_points)
  real(dp)           :: min_dx, max_dx, abs_tol

  call iu_read_grid('test_data/triangle.binda', ug)
  call iu_add_point_data(ug, 'vx', i_vx)
  call iu_add_point_data(ug, 'vy', i_vy)

  do n = 1, ug%n_points
     ug%point_data(n, i_vx) = -ug%points(2, n)
     ug%point_data(n, i_vy) = ug%points(1, n)
  end do

  call iu_write_vtk(ug, 'test_data/test_trace_field.vtu')

  r_start = [1.5_dp, 0.0_dp]
  abs_tol = 1e-3_dp
  min_dx  = 1e-4_dp
  max_dx  = 1.0e-1_dp

  call iu_trace_field(ug, 2, r_start, [i_vx, i_vy], max_points, n_points, &
       points, fields, min_dx, max_dx, abs_tol)

  if (n_points > max_points) error stop "Boundary not reached"

  call write_field_trace(n_points, points, fields, &
       'test_data/test_trace_field.3D')

contains

  subroutine write_field_trace(n_points, points, fields, filename)
    integer, intent(in)          :: n_points
    real(dp), intent(in)         :: points(2, n_points)
    real(dp), intent(in)         :: fields(2, n_points)
    character(len=*), intent(in) :: filename
    integer                      :: my_unit, n

    open(newunit=my_unit, file=trim(filename), action="write")
    write(my_unit, *) "x y z field_norm"
    do n = 1, n_points
       write(my_unit, *) points(:, n), 0.0_dp, norm2(fields(:, n))
    end do
    close(my_unit)

    print *, "Wrote ", trim(filename)
  end subroutine write_field_trace

end program test_trace_field
