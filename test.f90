program interpolate_unstructured
  use kdtree2_module

  implicit none
  integer, parameter    :: dp = kind(0.0d0)
  character(len=100)    :: fname

  type ugrid_t
     integer :: n_cells
     integer :: n_points
     integer :: n_points_per_cell
     integer :: n_faces_per_cell
     real(dp), allocatable :: points(:, :)
     integer, allocatable  :: cells(:, :)
     real(dp), allocatable :: cell_coords(:, :, :)
     integer, allocatable  :: neighbors(:, :)
     real(dp), allocatable :: cell_face_normals(:, :, :)
     real(dp), allocatable :: values(:, :)

     ! kd-tree for searching a nearby cell
     type(kdtree2), pointer :: tree => null()
  end type ugrid_t

  integer       :: n, k
  real(dp)      :: vec(3), normal_cell(3), normal_vec(3)
  type(ugrid_t) :: ug

  fname = 'Electric_Potential_needle'

  call read_ugrid(trim(fname), 1, ['Color'], ug)
  call set_cell_coords(ug)
  call set_normal_vectors(ug)

  print *, shape(ug%points)
  print *, shape(ug%cells)
  print *, shape(ug%neighbors)
  print *, shape(ug%values)

  call test_tracking(ug)

contains

  ! Create kd-tree to efficiently find a (nearby) cell for a new search at a
  ! new location
  subroutine build_kdtree(ug)
    type(ugrid_t), intent(inout) :: ug
    real(dp), allocatable        :: cell_centers(:, :)

    allocate(cell_centers(3, ug%n_cells))
    do n = 1, ug%n_cells
       cell_centers(:, n) = sum(ug%cell_coords(:, :, n), dim=2) / ug%n_points_per_cell
    end do

    ug%tree => kdtree2_create(cell_centers, sort=.false., rearrange=.false.)
  end subroutine build_kdtree

  ! Find a nearby cell at a new location
  subroutine find_nearby_cell_kdtree(ug, r, i_cell)
    type(ugrid_t), intent(in) :: ug
    real(dp), intent(in)      :: r(3)
    integer, intent(out)      :: i_cell
    type(kdtree2_result)      :: res(1)

    if (.not. associated(ug%tree)) error stop "Build tree first"

    call kdtree2_n_nearest(ug%tree, r, 1, res)
    i_cell = res(1)%idx
  end subroutine find_nearby_cell_kdtree

  subroutine read_ugrid(basename, n_variables, variable_names, ug)
    character(len=*), intent(in) :: basename
    integer, intent(in)          :: n_variables
    character(len=*), intent(in) :: variable_names(n_variables)
    type(ugrid_t), intent(out)   :: ug

    call read_array_float64_2d(trim(basename) // '_points.bin', ug%points)
    call read_array_int32_2d(trim(basename) // '_cells.bin', ug%cells)
    call read_array_int32_2d(trim(basename) // '_neighbors.bin', ug%neighbors)

    ! Store information about mesh
    ug%n_cells = size(ug%cells, 2)
    ug%n_points_per_cell = size(ug%cells, 1)
    ug%n_faces_per_cell = size(ug%cells, 1) ! TODO: only in 2d
    ug%n_points = size(ug%points, 2)

    ! Convert to 1-based indexing
    ug%cells = ug%cells + 1
    ug%neighbors = ug%neighbors + 1

    allocate(ug%values(ug%n_points, n_variables))

    ! Read data on points
    do n = 1, n_variables
       call read_array_float64_1d(trim(basename) // '_point_data_' // &
            trim(variable_names(n)) // '.bin', ug%n_points, ug%values(:, n))
    end do

  end subroutine read_ugrid

  ! Store points of each cell
  subroutine set_cell_coords(ug)
    type(ugrid_t), intent(inout) :: ug
    integer                      :: n, k

    allocate(ug%cell_coords(3, ug%n_points_per_cell, ug%n_cells))

    do n = 1, ug%n_cells
       do k = 1, ug%n_points_per_cell
          ug%cell_coords(:, k, n) = ug%points(:, ug%cells(k, n))
       end do
    end do
  end subroutine set_cell_coords

  ! Compute normal vectors to cell faces
  subroutine set_normal_vectors(ug)
    type(ugrid_t), intent(inout) :: ug
    integer                      :: n, k

    allocate(ug%cell_face_normals(3, ug%n_points_per_cell, ug%n_cells))

    do n = 1, ug%n_cells
       ! Compute vector normal to the cell
       normal_cell = cross_product_3d(&
            ug%cell_coords(:, 2, n) - ug%cell_coords(:, 1, n), &
            ug%cell_coords(:, 3, n) - ug%cell_coords(:, 2, n))

       do k = 1, ug%n_points_per_cell
          if (k < ug%n_points_per_cell) then
             vec = ug%cell_coords(:, k+1, n) - ug%cell_coords(:, k, n)
          else
             vec = ug%cell_coords(:, 1, n) - ug%cell_coords(:, k, n)
          end if

          ! TODO: order depends on arrangement of points
          normal_vec = cross_product_3d(vec, normal_cell)
          ug%cell_face_normals(:, k, n) = normal_vec / norm2(normal_vec)
       end do
    end do
  end subroutine set_normal_vectors

  subroutine test_tracking(ug)
    type(ugrid_t), intent(inout) :: ug
    real(dp)                     :: r0(3), r1(3)
    integer                      :: i0, i1

    i0 = 1
    r0 = sum(ug%cell_coords(:, :, i0), dim=2) / ug%n_points_per_cell
    r1 = [0.0_dp, 5.0_dp, 0.0_dp]
    call get_new_cell_index(ug, r0, r1, i0, i1)

    print *, "Point is inside cell: ", i1

    call build_kdtree(ug)

    i1 = 0
    call find_nearby_cell_kdtree(ug, r1, i1)

    print *, "Nearest cell from kd-tree: ", i1
  end subroutine test_tracking

  ! Helper function to compute the cross product of two 3D vectors
  function cross_product_3d(a, b) result(cross)
    real(dp), intent(in) :: a(3), b(3)
    real(dp)             :: cross(3)

    cross(1) = a(2) * b(3) - a(3) * b(2)
    cross(2) = a(3) * b(1) - a(1) * b(3)
    cross(3) = a(1) * b(2) - a(2) * b(1)
  end function cross_product_3d

  ! Determine the index of the cell at r1, given that r0 was in cell ic0
  subroutine get_new_cell_index(ug, r0, r1, ic0, ic1)
    type(ugrid_t), intent(in) :: ug
    real(dp), intent(in)      :: r0(3) ! Old position
    real(dp), intent(in)      :: r1(3) ! New position
    integer, intent(in)       :: ic0   ! Old cell index
    integer, intent(out)      :: ic1   ! New cell index

    real(dp) :: distance_left, path_unit_vec(3), r_p(3)
    real(dp) :: face_distance
    integer  :: i_cell, i_face

    distance_left = norm2(r1 - r0)
    path_unit_vec = (r1 - r0)/distance_left
    r_p           = r0
    i_cell        = ic0

    do
       call get_cell_intersection(ug, path_unit_vec, r_p, i_cell, &
            face_distance, i_face)

       distance_left = distance_left - face_distance

       if (distance_left > 0) then
          i_cell = ug%neighbors(i_face, i_cell)
          if (i_cell < 1) return ! Boundary cell
       else
          ! Done
          ic1 = i_cell
          exit
       end if
    end do
  end subroutine get_new_cell_index

  ! Given a cell, a direction, and a start position, determine through which
  ! cell face the path goes as well as the distance to the intersection point
  subroutine get_cell_intersection(ug, path_unit_vec, r_p, i_cell, &
       face_distance, i_face)
    type(ugrid_t), intent(in) :: ug
    real(dp), intent(in)      :: path_unit_vec(3)
    real(dp), intent(inout)   :: r_p(3)
    integer, intent(in)       :: i_cell
    real(dp), intent(out)     :: face_distance
    integer, intent(out)      :: i_face
    real(dp)                  :: path_dot_n, r_face(3), face_normal(3), dist
    real(dp), parameter       :: huge_distance = 1e100_dp

    face_distance = huge_distance
    i_face = -1

    do k = 1, ug%n_faces_per_cell
       face_normal = ug%cell_face_normals(:, k, i_cell)
       path_dot_n = dot_product(path_unit_vec, face_normal)

       ! Only consider faces whose normal points towards path
       if (path_dot_n > 0) then
          ! TODO: point on the cell face, generalize for 3d
          r_face = ug%cell_coords(:, k, i_cell)

          dist = dot_product(r_face - r_p, face_normal) / path_dot_n

          if (dist < face_distance) then
             face_distance = dist
             i_face = k
          end if
       end if
    end do

    r_p = r_p + face_distance * path_unit_vec

  end subroutine get_cell_intersection

  subroutine read_array_float64_1d(filename, array_size, array)
    use iso_fortran_env, only: error_unit
    character(len=*), intent(in) :: filename
    integer, intent(in)          :: array_size
    real(dp), intent(out)        :: array(array_size)

    integer :: my_unit, ndim, data_shape(1)
    character(len=64) :: dtype_str

    ! Open the binary file
    open(newunit=my_unit, file=filename, form='unformatted', &
         access='stream', status='old')

    read(my_unit) dtype_str
    if (dtype_str /= 'float64') then
       write(error_unit, *) ' Found dtype: ', trim(dtype_str)
       error stop 'Wrong dtype, expected float64'
    end if

    read(my_unit) ndim
    if (ndim /= 1) then
       write(error_unit, *) ' Found ndim: ', ndim
       error stop 'Wrong dimension, expected 1'
    end if

    read(my_unit) data_shape

    if (data_shape(1) /= array_size) then
       write(error_unit, *) ' Found data_shape: ', data_shape(1), &
            ', expected: ', array_size
       error stop 'Wrong data_shape'
    end if

    read(my_unit) array
    close(my_unit)

  end subroutine read_array_float64_1d

  subroutine read_array_float64_2d(filename, array)
    use iso_fortran_env, only: error_unit
    character(len=*), intent(in) :: filename
    real(dp), allocatable        :: array(:, :)

    integer :: my_unit, ndim, data_shape(2)
    character(len=64) :: dtype_str

    ! Open the binary file
    open(newunit=my_unit, file=filename, form='unformatted', &
         access='stream', status='old')

    read(my_unit) dtype_str
    if (dtype_str /= 'float64') then
       write(error_unit, *) ' Found dtype: ', trim(dtype_str)
       error stop 'Wrong dtype, expected float64'
    end if

    read(my_unit) ndim
    if (ndim /= 2) then
       write(error_unit, *) ' Found ndim: ', ndim
       error stop 'Wrong dimension, expected 2'
    end if

    read(my_unit) data_shape
    allocate(array(data_shape(2), data_shape(1)))
    read(my_unit) array

    close(my_unit)

  end subroutine read_array_float64_2d

  subroutine read_array_int32_2d(filename, array)
    use iso_fortran_env, only: error_unit
    character(len=*), intent(in) :: filename
    integer, allocatable         :: array(:, :)

    integer :: my_unit, ndim, data_shape(2)
    character(len=64) :: dtype_str

    ! Open the binary file
    open(newunit=my_unit, file=filename, form='unformatted', &
         access='stream', status='old')

    read(my_unit) dtype_str
    if (dtype_str /= 'int32') then
       write(error_unit, *) ' Found dtype: ', trim(dtype_str)
       error stop 'Wrong dtype, expected int32'
    end if

    read(my_unit) ndim
    if (ndim /= 2) then
       write(error_unit, *) ' Found ndim: ', ndim
       error stop 'Wrong dimension, expected 2'
    end if

    read(my_unit) data_shape
    allocate(array(data_shape(2), data_shape(1)))
    read(my_unit) array

    close(my_unit)

  end subroutine read_array_int32_2d

  subroutine interpolate_triangle(p1, p2, p3, v1, v2, v3, point, value)
    real(dp), intent(in) :: p1(3), p2(3), p3(3)  ! vertices of the triangle
    real(dp), intent(in) :: v1, v2, v3            ! values at vertices
    real(dp), intent(in) :: point(3)               ! point for interpolation
    real(dp), intent(out) :: value                  ! interpolated value

    real(dp) :: area, area1, area2, area3
    real(dp) :: lambda1, lambda2, lambda3

    ! Compute areas to find barycentric coordinates
    area = 0.5 * norm2(cross_product_3d(p2 - p1, p3 - p1))

    area1 = 0.5 * norm2(cross_product_3d(point - p2, point - p3))
    area2 = 0.5 * norm2(cross_product_3d(point - p3, point - p1))
    area3 = 0.5 * norm2(cross_product_3d(point - p1, point - p2))

    ! Calculate barycentric coordinates
    lambda1 = area1 / area
    lambda2 = area2 / area
    lambda3 = area3 / area

    ! Perform interpolation
    value = lambda1 * v1 + lambda2 * v2 + lambda3 * v3
  end subroutine interpolate_triangle

end program
