program test_triangle
  use iso_fortran_env, only: error_unit
  use m_interp_unstructured

  implicit none
  integer, parameter :: dp = kind(0.0d0)
  character(len=100) :: fname
  type(iu_grid_t)    :: ug

  fname = 'test_data/tetra'
  call iu_read_grid(trim(fname), 1, ['Polynomial'], ug)

  call test_interpolation(ug, 1000)

contains

  subroutine test_interpolation(ug, n_samples)
    type(iu_grid_t), intent(inout) :: ug
    integer, intent(in)            :: n_samples

    integer               :: n
    real(dp), allocatable :: r_samples(:, :), res(:)
    integer, allocatable  :: i_cell(:)
    real(dp)              :: rmin(3), rmax(3), difference
    real(dp), parameter   :: threshold = 1e-14_dp

    allocate(r_samples(3, n_samples), res(n_samples), i_cell(n_samples))
    rmin = minval(ug%points, dim=2)
    rmax = maxval(ug%points, dim=2)

    do n = 1, n_samples
       call random_number(r_samples(:, n))
       r_samples(:, n) = rmin + r_samples(:, n) * (rmax - rmin)
    end do

    ! Set solution at points
    do n = 1, ug%n_points
       ug%values(n, 1) = solution(ug%points(:, n))
    end do

    i_cell(:) = 0

    do n = 1, n_samples
       call iu_interpolate_at(ug, r_samples(:, n), 1, res(n), i_cell(n))
       difference = abs(res(n) - solution(r_samples(:, n)))

       if (difference > threshold) then
          write(error_unit, *) "Difference of ", difference, &
               " at ", r_samples(:, n)
          error stop "Difference exceeds threshold"
       end if
    end do

  end subroutine test_interpolation

  real(dp) function solution(r)
    real(dp), intent(in) :: r(3)
    solution = sum(r) + 1
  end function solution

end program test_triangle
