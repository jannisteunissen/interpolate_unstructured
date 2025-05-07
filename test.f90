program interpolate_unstructured
  use m_interp_unstructured

  implicit none
  integer, parameter :: dp = kind(0.0d0)
  character(len=100) :: fname
  type(iu_grid_t)    :: ug

  fname = 'Electric_Potential_needle'
  call iu_read_grid(trim(fname), 1, ['Color'], ug)

  call test_tracking(ug)

contains

  subroutine test_tracking(ug)
    type(iu_grid_t), intent(inout) :: ug
    real(dp)                       :: r1(3), res
    integer                        :: i_cell

    r1 = [0.0e-3_dp, 5.0_dp, 0.0_dp]
    i_cell = 0
    call iu_interpolate_at(ug, r1, 1, res, i_cell)

    print *, "value at ", r1, " is ", res, " i_cell is ", i_cell
  end subroutine test_tracking

end program
