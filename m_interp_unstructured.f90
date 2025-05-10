module m_interp_unstructured
  use iso_fortran_env, only: error_unit
  use kdtree2_module

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)
  integer, parameter :: sp = kind(0.0e0)

  integer, parameter :: ug_triangle = 1
  integer, parameter :: ug_quad = 2
  integer, parameter :: ug_tetra = 3

  integer, parameter :: ug_max_points_per_cell = 4

  real(dp), parameter :: tiny_distance = 1e-100_dp

  type iu_grid_t
     integer :: n_cells
     integer :: n_points
     integer :: n_points_per_cell
     integer :: n_faces_per_cell
     integer :: cell_type
     real(dp), allocatable :: points(:, :)
     integer, allocatable  :: cells(:, :)
     real(dp), allocatable :: cell_points(:, :, :)
     integer, allocatable  :: neighbors(:, :)
     real(dp), allocatable :: cell_face_normals(:, :, :)
     real(dp), allocatable :: values(:, :)

     ! kd-tree for searching a nearby cell
     type(kdtree2) :: tree
  end type iu_grid_t

  ! Public types
  public :: iu_grid_t

  ! Public methods
  public :: iu_read_grid
  public :: iu_interpolate_at

contains

  ! Create kd-tree to efficiently find a (nearby) cell for a new search at a
  ! new location
  subroutine build_kdtree(ug)
    type(iu_grid_t), intent(inout) :: ug
    real(dp), allocatable          :: cell_centers(:, :)
    integer                        :: n

    allocate(cell_centers(3, ug%n_cells))
    do n = 1, ug%n_cells
       cell_centers(:, n) = sum(ug%cell_points(:, :, n), dim=2) / &
            ug%n_points_per_cell
    end do

    ug%tree = kdtree2_create(cell_centers, sort=.false., rearrange=.false.)
  end subroutine build_kdtree

  ! Find a nearby cell at a new location
  integer function find_nearby_cell_kdtree(ug, r) result(i_cell)
    type(iu_grid_t), intent(in) :: ug
    real(dp), intent(in)        :: r(3)
    type(kdtree2_result)        :: res(1)

    if (ug%tree%n == 0) error stop "Build tree first"

    call kdtree2_n_nearest(ug%tree, r, 1, res)
    i_cell = res(1)%idx
  end function find_nearby_cell_kdtree

  subroutine iu_read_grid(basename, n_variables, variable_names, ug)
    character(len=*), intent(in) :: basename
    integer, intent(in)          :: n_variables
    character(len=*), intent(in) :: variable_names(n_variables)
    type(iu_grid_t), intent(out) :: ug
    integer                      :: n, my_unit
    character(len=64)            :: cell_type

    open(newunit=my_unit, file=trim(basename)//'_cell_type.txt', status='old')
    read(my_unit, *) cell_type
    close(my_unit)

    select case (cell_type)
    case ("triangle")
       ug%cell_type = ug_triangle
    case ("quad")
       ug%cell_type = ug_quad
    case ("tetra")
       ug%cell_type = ug_tetra
    case default
       write(error_unit, *) "Cell type '", trim(cell_type), "' not supported"
       error stop "Unsupported cell type"
    end select

    call read_array_float64_2d(trim(basename)//'_points.bin', ug%points)
    call read_array_int32_2d(trim(basename)//'_cells.bin', ug%cells)
    call read_array_int32_2d(trim(basename)//'_neighbors.bin', ug%neighbors)

    ! Store information about mesh
    ug%n_cells = size(ug%cells, 2)
    ug%n_points_per_cell = size(ug%cells, 1)
    ug%n_faces_per_cell = size(ug%cells, 1) ! works for triangle, quad, tetra
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

    ! Store points for each cell
    call set_cell_points(ug)

    ! Store normal vectors for each cell face
    call set_face_normal_vectors(ug)

    ! Build kd-tree for efficient lookup of new points
    call build_kdtree(ug)

  end subroutine iu_read_grid

  ! Store points of each cell
  subroutine set_cell_points(ug)
    type(iu_grid_t), intent(inout) :: ug
    integer                        :: n, k

    allocate(ug%cell_points(3, ug%n_points_per_cell, ug%n_cells))

    do n = 1, ug%n_cells
       do k = 1, ug%n_points_per_cell
          ug%cell_points(:, k, n) = ug%points(:, ug%cells(k, n))
       end do
    end do
  end subroutine set_cell_points

  ! Compute normal vectors to cell faces
  subroutine set_face_normal_vectors(ug)
    type(iu_grid_t), intent(inout) :: ug
    integer                        :: n, k, k1, k2
    real(dp)                       :: vec(3), center(3)
    real(dp)                       :: normal_cell(3), normal_face(3)

    allocate(ug%cell_face_normals(3, ug%n_points_per_cell, ug%n_cells))

    select case (ug%cell_type)
    case (ug_triangle, ug_quad)
       do n = 1, ug%n_cells
          center = sum(ug%cell_points(:, :, n), dim=2) / ug%n_points_per_cell

          ! Compute vector normal to the cell, assuming it to be flat
          normal_cell = cross_product(&
               ug%cell_points(:, 2, n) - ug%cell_points(:, 1, n), &
               ug%cell_points(:, 3, n) - ug%cell_points(:, 2, n))

          do k = 1, ug%n_points_per_cell
             k1 = mod(k, ug%n_points_per_cell) + 1

             vec = ug%cell_points(:, k1, n) - ug%cell_points(:, k, n)
             normal_face = cross_product(vec, normal_cell)

             ! Ensure that face is pointing outwards
             vec = ug%cell_points(:, k, n) - center
             if (dot_product(vec, normal_face) < 0) normal_face = -normal_face

             ug%cell_face_normals(:, k, n) = normal_face / norm2(normal_face)
          end do
       end do
    case (ug_tetra)
       do n = 1, ug%n_cells
          center = sum(ug%cell_points(:, :, n), dim=2) / ug%n_points_per_cell

          do k = 1, ug%n_points_per_cell
             k1 = mod(k, ug%n_points_per_cell) + 1
             k2 = mod(k+1, ug%n_points_per_cell) + 1

             normal_face = cross_product(&
                  ug%cell_points(:, k1, n) - ug%cell_points(:, k, n), &
                  ug%cell_points(:, k2, n) - ug%cell_points(:, k1, n))

             ! Ensure that face is pointing outwards
             vec = ug%cell_points(:, k, n) - center
             if (dot_product(vec, normal_face) < 0) normal_face = -normal_face

             ug%cell_face_normals(:, k, n) = normal_face / norm2(normal_face)
          end do
       end do
    case default
       error stop 'Not implemented'
    end select

  end subroutine set_face_normal_vectors

  ! Find cell containing r, optionally using a nearby cell as a starting point
  integer function find_containing_cell(ug, r, i_guess) result(i_cell)
    type(iu_grid_t), intent(in) :: ug
    real(dp), intent(in)        :: r(3)
    integer, intent(inout)      :: i_guess
    real(dp)                    :: r0(3)

    if (i_guess < 1) then
       i_guess = find_nearby_cell_kdtree(ug, r)
    end if

    ! Start from center of cell i_guess
    r0 = sum(ug%cell_points(:, :, i_guess), dim=2) / ug%n_points_per_cell

    call get_cell_through_neighbors(ug, r0, r, i_guess, i_cell)

    if (i_cell < 1) then
       write(error_unit, *) "ERROR: cannot find cell containing", r
       error stop "Point is probably outside domain"
    end if
  end function find_containing_cell

  ! Interpolate variable with index i_variable at location r
  subroutine iu_interpolate_at(ug, r, i_variable, var_at_r, i_cell)
    type(iu_grid_t), intent(in) :: ug
    real(dp), intent(in)        :: r(3)         ! Interpolate at this point
    integer, intent(in)         :: i_variable   ! Index of variable
    real(dp), intent(out)       :: var_at_r     ! Result at r
    ! On input: guess for nearby cell. Less than 1 means not set.
    ! On output: cell containing the point r.
    integer, intent(inout)      :: i_cell
    real(dp)                    :: values(ug_max_points_per_cell)

    if (i_cell > ug%n_cells) error stop "i_cell > ug%n_cells"

    i_cell = find_containing_cell(ug, r, i_cell)
    values(1:ug%n_points_per_cell) = ug%values(ug%cells(:, i_cell), i_variable)

    select case (ug%cell_type)
    case (ug_triangle)
       call interpolate_triangle(ug%cell_points(:, :, i_cell), &
            values(1:ug%n_points_per_cell), r, var_at_r)
    case (ug_quad)
       call interpolate_quad(ug%cell_points(:, :, i_cell), &
            values(1:ug%n_points_per_cell), r, var_at_r)
    case (ug_tetra)
       call interpolate_tetrahedron(ug%cell_points(:, :, i_cell), &
            values(1:ug%n_points_per_cell), r, var_at_r)
    case default
       error stop "Not implemented"
    end select
  end subroutine iu_interpolate_at

  subroutine interpolate_triangle(points, values, r, res)
    real(dp), intent(in)  :: points(3, 3) ! vertices of the triangle
    real(dp), intent(in)  :: values(3)    ! values at vertices
    real(dp), intent(in)  :: r(3)         ! point for interpolation
    real(dp), intent(out) :: res          ! interpolated value
    real(dp)              :: area_full, areas(3)

    ! Compute areas to find barycentric coordinates. The factors 0.5 cancel in
    ! the final computation and have been left out.
    area_full = norm2(cross_product(points(:, 2) - points(:, 1), &
         points(:, 3) - points(:, 1)))
    areas(1) = norm2(cross_product(r - points(:, 2), r - points(:, 3)))
    areas(2) = norm2(cross_product(r - points(:, 3), r - points(:, 1)))
    areas(3) = norm2(cross_product(r - points(:, 1), r - points(:, 2)))

    ! Perform interpolation
    res = sum(areas * values) / area_full
  end subroutine interpolate_triangle

  ! Based on https://www.cdsimpson.net/2014/10/barycentric-coordinates.html
  ! https://stackoverflow.com/questions/38545520/barycentric-coordinates-of-a-tetrahedron
  subroutine interpolate_tetrahedron(points, values, r, res)
    real(dp), intent(in)  :: points(3, 4) ! vertices of the triangle
    real(dp), intent(in)  :: values(4)    ! values at vertices
    real(dp), intent(in)  :: r(3)         ! point for interpolation
    real(dp), intent(out) :: res          ! interpolated value
    real(dp)              :: weights(4), full_weight
    real(dp), dimension(3) :: v1r, v2r, v12, v13, v14, v23, v24

    v1r = r - points(:, 1)
    v2r = r - points(:, 2)
    v12 = points(:, 2) - points(:, 1)
    v13 = points(:, 3) - points(:, 1)
    v14 = points(:, 4) - points(:, 1)
    v23 = points(:, 3) - points(:, 2)
    v24 = points(:, 4) - points(:, 2)

    weights(1) = scalar_triple_product(v2r, v24, v23)
    weights(2) = scalar_triple_product(v1r, v13, v14)
    weights(3) = scalar_triple_product(v1r, v14, v12)
    weights(4) = scalar_triple_product(v1r, v12, v13)
    full_weight = scalar_triple_product(v12, v13, v14)

    ! Perform interpolation
    res = sum(weights * values) / full_weight

  end subroutine interpolate_tetrahedron

  ! Interpolate a quad element with four points at arbitary locations. The
  ! implementation is based on:
  ! https://www.reedbeta.com/blog/quadrilateral-interpolation-part-2/
  subroutine interpolate_quad(points, values, r, res)
    real(dp), intent(in)  :: points(3, 4) ! vertices of the triangle
    real(dp), intent(in)  :: values(4)    ! values at vertices
    real(dp), intent(in)  :: r(3)         ! point for interpolation
    real(dp), intent(out) :: res          ! interpolated value

    real(dp) :: q(3), b1(3), b2(3), b3(3), denom(3)
    real(dp) :: A, B, C, discrim, coeff(2), tmp(2)
    integer  :: max_dim(1)
    real(dp), parameter :: tiny_value = 1e-20_dp

    ! Below, points(:,1) is the first neighbor of points(:,1), and points(:,4)
    ! is the second neighbor, whereas points(:,3) is the diagonal neighbor
    q = r - points(:, 1)
    b1 = points(:, 2) - points(:, 1) ! First neighbor
    b2 = points(:, 4) - points(:, 1) ! Second neighbor
    b3 = points(:, 1) - points(:, 2) - points(:, 4) + points(:, 3)

    ! Set up quadratic formula
    A = cross_product_z(b2, b3)
    B = cross_product_z(b3, q) - cross_product_z(b1, b2)
    C = cross_product_z(b1, q)
    discrim = B**2 - 4 * A * C

    ! Solve quadratic equation
    if (abs(A) < tiny_value) then
       coeff(2) = -C/B
    else
       coeff(2) = 0.5_dp * (-B - sqrt(discrim)) / A
    end if

    denom = b1 + coeff(2) * b3

    ! The nominator and denominator should be parallel, so we can use a single
    ! coordinate for the division below
    max_dim = maxloc(abs(denom))

    associate (dim => max_dim(1))
      coeff(1) = (q(dim) - b2(dim) * coeff(2)) / denom(dim)
    end associate

    ! Perform bilinear interpolation using the found coefficients
    tmp(1) = values(1) * (1 - coeff(1)) + values(2) * coeff(1)
    tmp(2) = values(4) * (1 - coeff(1)) + values(3) * coeff(1)
    res = tmp(1) * (1 - coeff(2)) + tmp(2) * coeff(2)

  end subroutine interpolate_quad

  ! Helper function to compute the cross product of two 3D vectors
  pure function cross_product(a, b) result(cross)
    real(dp), intent(in) :: a(3), b(3)
    real(dp)             :: cross(3)

    cross(1) = a(2) * b(3) - a(3) * b(2)
    cross(2) = a(3) * b(1) - a(1) * b(3)
    cross(3) = a(1) * b(2) - a(2) * b(1)
  end function cross_product

  ! Return z-coordinate of cross product, useful if a(3) and b(3) are zero
  pure real(dp) function cross_product_z(a, b)
    real(dp), intent(in) :: a(3), b(3)
    cross_product_z = a(1) * b(2) - a(2) * b(1)
  end function cross_product_z

  pure real(dp) function scalar_triple_product(a, b, c) result(stp)
    real(dp), intent(in) :: a(3), b(3), c(3)
    stp = dot_product(a, cross_product(b, c))
  end function scalar_triple_product

  ! Determine the index of the cell at r1, given that r0 was in cell ic0
  subroutine get_cell_through_neighbors(ug, r0, r1, ic0, ic1)
    type(iu_grid_t), intent(in) :: ug
    real(dp), intent(in)        :: r0(3) ! Old position
    real(dp), intent(in)        :: r1(3) ! New position
    integer, intent(in)         :: ic0   ! Old cell index
    integer, intent(out)        :: ic1   ! New cell index

    real(dp) :: distance_left, path_unit_vec(3), r_p(3)
    real(dp) :: face_distance
    integer  :: i_face

    distance_left = norm2(r1 - r0)

    if (distance_left < tiny_distance) then
       ic1 = ic0
       return
    end if

    path_unit_vec = (r1 - r0) / distance_left
    r_p           = r0          ! Current position
    ic1           = ic0         ! Current cell

    do
       call get_cell_intersection(ug, path_unit_vec, r_p, ic1, &
            face_distance, i_face)

       distance_left = distance_left - face_distance

       if (distance_left > 0) then
          ic1 = ug%neighbors(i_face, ic1)
          if (ic1 < 1) exit     ! Exit when moving outside domain
       else
          exit                  ! Done
       end if
    end do

  end subroutine get_cell_through_neighbors

  ! Given a cell, a direction, and a start position, determine through which
  ! cell face the path goes as well as the distance to the intersection point
  subroutine get_cell_intersection(ug, path_unit_vec, r_p, i_cell, &
       face_distance, i_face)
    type(iu_grid_t), intent(in) :: ug
    real(dp), intent(in)        :: path_unit_vec(3)
    real(dp), intent(inout)     :: r_p(3)
    integer, intent(in)         :: i_cell
    real(dp), intent(out)       :: face_distance
    integer, intent(out)        :: i_face
    real(dp)                    :: path_dot_n, r_face(3), face_normal(3), dist
    real(dp), parameter         :: huge_distance = 1e100_dp
    integer                     :: k

    face_distance = huge_distance
    i_face = -1

    do k = 1, ug%n_faces_per_cell
       face_normal = ug%cell_face_normals(:, k, i_cell)
       path_dot_n = dot_product(path_unit_vec, face_normal)

       ! Only consider faces whose normal points towards path_unit_vec
       if (path_dot_n > 0) then
          ! Point on the cell face, works for triangle, quad, tetra
          r_face = ug%cell_points(:, k, i_cell)

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
    real(sp), allocatable :: sp_array(:)

    ! Open the binary file
    open(newunit=my_unit, file=filename, form='unformatted', &
         access='stream', status='old')

    read(my_unit) dtype_str

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

    select case (dtype_str)
    case ('float64')
       read(my_unit) array
    case ('float32')
       allocate(sp_array(array_size))
       read(my_unit) sp_array
       array = sp_array
    case default
       write(error_unit, *) ' Found dtype: ', trim(dtype_str)
       error stop 'Wrong dtype, expected float64 or float32'
    end select

    close(my_unit)

  end subroutine read_array_float64_1d

  subroutine read_array_float64_2d(filename, array)
    use iso_fortran_env, only: error_unit
    character(len=*), intent(in) :: filename
    real(dp), allocatable        :: array(:, :)

    integer :: my_unit, ndim, data_shape(2)
    character(len=64) :: dtype_str
    real(sp), allocatable :: sp_array(:, :)

    ! Open the binary file
    open(newunit=my_unit, file=filename, form='unformatted', &
         access='stream', status='old')

    read(my_unit) dtype_str

    read(my_unit) ndim
    if (ndim /= 2) then
       write(error_unit, *) ' Found ndim: ', ndim
       error stop 'Wrong dimension, expected 2'
    end if

    read(my_unit) data_shape
    allocate(array(data_shape(2), data_shape(1)))

    select case (dtype_str)
    case ('float64')
       read(my_unit) array
    case ('float32')
       allocate(sp_array(data_shape(2), data_shape(1)))
       read(my_unit) sp_array
       array = sp_array
    case default
       write(error_unit, *) ' Found dtype: ', trim(dtype_str)
       error stop 'Wrong dtype, expected float64 or float32'
    end select

    close(my_unit)

  end subroutine read_array_float64_2d

  subroutine read_array_int32_2d(filename, array)
    use iso_fortran_env, only: error_unit, int64
    character(len=*), intent(in) :: filename
    integer, allocatable         :: array(:, :)

    integer :: my_unit, ndim, data_shape(2)
    character(len=64) :: dtype_str
    integer(int64), allocatable :: int64_array(:, :)

    ! Open the binary file
    open(newunit=my_unit, file=filename, form='unformatted', &
         access='stream', status='old')

    read(my_unit) dtype_str

    read(my_unit) ndim
    if (ndim /= 2) then
       write(error_unit, *) ' Found ndim: ', ndim
       error stop 'Wrong dimension, expected 2'
    end if

    read(my_unit) data_shape
    allocate(array(data_shape(2), data_shape(1)))

    select case (dtype_str)
    case ('int64')
       allocate(int64_array(data_shape(2), data_shape(1)))
       read(my_unit) int64_array
       array = int(int64_array)
    case ('int32')
       read(my_unit) array
    case default
       write(error_unit, *) ' Found dtype: ', trim(dtype_str)
       error stop 'Wrong dtype, expected float64 or float32'
    end select

    close(my_unit)

  end subroutine read_array_int32_2d

end module m_interp_unstructured
