program benchmark
  use m_interp_unstructured

  implicit none
  integer, parameter :: dp = kind(0.0d0)
  character(len=100) :: fname
  type(iu_grid_t)    :: ug

  fname = 'test_data/triangle'
  call iu_read_grid(trim(fname), 1, ['Polynomial'], ug)

  call run_benchmark(ug, 1000*1000)

contains

  subroutine run_benchmark(ug, n_samples)
    use iso_fortran_env, only: int64
    type(iu_grid_t), intent(inout) :: ug
    integer, intent(in)            :: n_samples
    real(dp), allocatable          :: r_samples(:, :), velocity(:, :), res(:)
    integer, allocatable           :: i_cell(:)
    real(dp)                       :: rmin(3), rmax(3), domain_size(3)
    real(dp)                       :: cpu_time, dt
    integer                        :: n
    integer(int64)                 :: t_start, t_end, count_rate

    allocate(r_samples(3, n_samples), res(n_samples), i_cell(n_samples))
    allocate(velocity(3, n_samples))

    rmin = minval(ug%points, dim=2)
    rmax = maxval(ug%points, dim=2)
    domain_size = rmax - rmin

    ! Ensure there is some extra space for points to move
    rmin = rmin + 0.1_dp * domain_size
    rmax = rmax - 0.1_dp * domain_size

    do n = 1, n_samples
       call random_number(r_samples(:, n))
       r_samples(:, n) = rmin + r_samples(:, n) * (rmax - rmin)

       call random_number(velocity(:, n))
    end do

    i_cell(:) = 0

    call system_clock(t_start, count_rate)
    do n = 1, n_samples
       call iu_interpolate_at(ug, r_samples(:, n), 1, res(n), i_cell(n))
    end do
    call system_clock(t_end, count_rate)

    cpu_time = (t_end-t_start) / real(count_rate, dp)
    write(*, "(A,I0,A,E10.3,A)") " Wall-clock for ", n_samples, &
       " samples: ", cpu_time, " seconds"

    dt = 0.01_dp * minval(domain_size)
    r_samples = r_samples + dt * velocity

    call system_clock(t_start, count_rate)
    do n = 1, n_samples
       call iu_interpolate_at(ug, r_samples(:, n), 1, res(n), i_cell(n))
    end do
    call system_clock(t_end, count_rate)

    cpu_time = (t_end-t_start) / real(count_rate, dp)
    write(*, "(A,I0,A,E10.3,A)") " Wall-clock for ", n_samples, &
       " samples: ", cpu_time, " seconds"

  end subroutine run_benchmark

end program benchmark
