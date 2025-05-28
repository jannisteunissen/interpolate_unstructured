module m_interp_unstructured
  use iso_fortran_env, only: error_unit, int64
  use kdtree2_module
  use m_binda
  use m_vtk

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)
  integer, parameter :: sp = kind(0.0e0)

  integer, parameter :: iu_triangle = 1
  integer, parameter :: iu_quad = 2
  integer, parameter :: iu_tetra = 3
  integer, parameter :: iu_ndim_cell_type(3) = [2, 2, 3]

  integer, parameter :: ug_max_points_per_cell = 4

  real(dp), parameter :: tiny_distance = 1e-100_dp

  ! Type for storing an unstructured grid with additional data for efficient
  ! interpolation
  type iu_grid_t
     integer :: n_cells
     integer :: n_points
     integer :: n_points_per_cell
     integer :: n_faces_per_cell
     integer :: cell_type
     integer :: n_point_data
     integer :: n_cell_data

     real(dp) :: rmin(3)
     real(dp) :: rmax(3)

     real(dp), allocatable :: points(:, :)
     integer, allocatable  :: cells(:, :)
     real(dp), allocatable :: cell_points(:, :, :)
     integer, allocatable  :: neighbors(:, :)
     real(dp), allocatable :: cell_face_normals(:, :, :)
     real(dp), allocatable :: point_data(:, :)
     real(dp), allocatable :: cell_data(:, :)

     ! Names of the point data variables
     character(len=128), allocatable :: point_data_names(:)

     ! Names of the cell data variables
     character(len=128), allocatable :: cell_data_names(:)

     ! kd-tree for searching a nearby cell
     type(kdtree2) :: tree
  end type iu_grid_t

  ! Public types
  public :: iu_grid_t
  public :: iu_triangle, iu_quad, iu_tetra
  public :: iu_ndim_cell_type

  ! Public methods
  public :: iu_read_grid
  public :: iu_get_point_data_index
  public :: iu_get_cell_data_index
  public :: iu_add_cell_data
  public :: iu_add_point_data
  public :: iu_get_cell_center
  public :: iu_get_cell
  public :: iu_get_cell_scalar_at
  public :: iu_trace_field
  public :: iu_interpolate_at
  public :: iu_interpolate_scalar_at
  public :: iu_write_vtk

contains

  ! Find index of point data variable, -1 if not present
  subroutine iu_get_point_data_index(ug, name, ix)
    type(iu_grid_t), intent(in)  :: ug
    character(len=*), intent(in) :: name
    integer, intent(out)         :: ix

    do ix = 1, ug%n_point_data
       if (ug%point_data_names(ix) == name) exit
    end do

    if (ix == ug%n_point_data + 1) ix = -1
  end subroutine iu_get_point_data_index

  !> Extend cell data with one variable
  subroutine iu_add_cell_data(ug, name, i_var)
    type(iu_grid_t), intent(inout)  :: ug
    character(len=*), intent(in)    :: name
    integer, intent(out)            :: i_var
    real(dp), allocatable           :: old_data(:, :)
    character(len=128), allocatable :: old_names(:)

    call move_alloc(ug%cell_data, old_data)
    call move_alloc(ug%cell_data_names, old_names)
    allocate(ug%cell_data(ug%n_cells, ug%n_cell_data+1))
    allocate(ug%cell_data_names(ug%n_cell_data+1))

    ug%cell_data(:, 1:ug%n_cell_data) = old_data
    ug%cell_data_names(1:ug%n_cell_data) = old_names

    ug%n_cell_data = ug%n_cell_data + 1
    ug%cell_data_names(ug%n_cell_data) = name
    i_var = ug%n_cell_data
  end subroutine iu_add_cell_data

  !> Extend point data with one variable
  subroutine iu_add_point_data(ug, name, i_var)
    type(iu_grid_t), intent(inout)  :: ug
    character(len=*), intent(in)    :: name
    integer, intent(out)            :: i_var
    real(dp), allocatable           :: old_data(:, :)
    character(len=128), allocatable :: old_names(:)

    call move_alloc(ug%point_data, old_data)
    call move_alloc(ug%point_data_names, old_names)
    allocate(ug%point_data(ug%n_points, ug%n_point_data+1))
    allocate(ug%point_data_names(ug%n_point_data+1))

    ug%point_data(:, 1:ug%n_point_data) = old_data
    ug%point_data_names(1:ug%n_point_data) = old_names

    ug%n_point_data = ug%n_point_data + 1
    ug%point_data_names(ug%n_point_data) = name
    i_var = ug%n_point_data
  end subroutine iu_add_point_data

  ! Find index of cell data variable, -1 if not present
  subroutine iu_get_cell_data_index(ug, name, ix)
    type(iu_grid_t), intent(in)  :: ug
    character(len=*), intent(in) :: name
    integer, intent(out)         :: ix

    do ix = 1, ug%n_cell_data
       if (ug%cell_data_names(ix) == name) exit
    end do

    if (ix == ug%n_cell_data + 1) ix = -1
  end subroutine iu_get_cell_data_index

  ! Create kd-tree to efficiently find a (nearby) cell for a new search at a
  ! new location
  subroutine build_kdtree(ug)
    type(iu_grid_t), intent(inout) :: ug
    real(dp), allocatable          :: cell_centers(:, :)
    integer                        :: n

    allocate(cell_centers(3, ug%n_cells))
    do n = 1, ug%n_cells
       cell_centers(:, n) = iu_get_cell_center(ug, n)
    end do

    ug%tree = kdtree2_create(cell_centers, sort=.false., rearrange=.false.)
  end subroutine build_kdtree

  pure function iu_get_cell_center(ug, i_cell) result(center)
    type(iu_grid_t), intent(in) :: ug
    integer, intent(in)         :: i_cell
    real(dp)                    :: center(3)
    center = sum(ug%cell_points(:, :, i_cell), dim=2) / ug%n_points_per_cell
  end function iu_get_cell_center

  ! Find a nearby cell at a new location
  integer function find_nearby_cell_kdtree(ug, r) result(i_cell)
    type(iu_grid_t), intent(in) :: ug
    real(dp), intent(in)        :: r(3)
    type(kdtree2_result)        :: res(1)

    if (ug%tree%n == 0) error stop "Build tree first"

    call kdtree2_n_nearest(ug%tree, r, 1, res)
    i_cell = res(1)%idx

    if (i_cell < 1 .or. i_cell > ug%n_cells) then
       write(error_unit, *) "r = ", r
       write(error_unit, *) "i_cell = ", i_cell
       error stop "Error in kdtree search"
    end if

  end function find_nearby_cell_kdtree

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
    case (iu_triangle, iu_quad)
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
    case (iu_tetra)
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
       error stop "Not implemented"
    end select

  end subroutine set_face_normal_vectors

  ! Find cell containing r. If the argument i_cell > 1, use cell i_cell as a
  ! starting point.
  subroutine iu_get_cell(ug, r, i_cell)
    type(iu_grid_t), intent(in) :: ug
    real(dp), intent(in)        :: r(3)
    !> Input: guess for i_cell or value < 1. Output: The cell containing r
    integer, intent(inout)      :: i_cell
    integer                     :: i_start
    real(dp)                    :: r0(3)

    if (i_cell < 1) then
       i_start = find_nearby_cell_kdtree(ug, r)
    else
       i_start = i_cell
    end if

    ! Start from center of cell i_cell
    r0 = sum(ug%cell_points(:, :, i_start), dim=2) / ug%n_points_per_cell

    call get_cell_through_neighbors(ug, r0, r, i_start, i_cell)
  end subroutine iu_get_cell

  ! Get the value of cell data at the cell that contains r
  subroutine iu_get_cell_scalar_at(ug, r, i_var, i_guess, res)
    type(iu_grid_t), intent(in) :: ug
    real(dp), intent(in)        :: r(3)  ! Location
    integer, intent(in)         :: i_var ! Index of cell data variable
    ! On input: guess for nearby cell. Less than 1 means not set.
    ! On output: cell containing the point r.
    integer, intent(inout)      :: i_guess
    real(dp), intent(inout)     :: res ! Result

    call iu_get_cell(ug, r, i_guess)
    if (i_guess >= 1) res = ug%cell_data(i_guess, i_var)
  end subroutine iu_get_cell_scalar_at

  ! Interpolate scalar at location r
  subroutine iu_interpolate_scalar_at(ug, r, i_var, var_at_r, i_cell)
    type(iu_grid_t), intent(in) :: ug
    real(dp), intent(in)        :: r(3)     ! Interpolate at this point
    integer, intent(in)         :: i_var    ! Index of variables
    real(dp), intent(out)       :: var_at_r ! Result at r
    ! On input: guess for nearby cell. Less than 1 means not set.
    ! On output: cell containing the point r.
    integer, intent(inout)      :: i_cell
    real(dp)                    :: tmp(1)

    call iu_interpolate_at(ug, r, 1, [i_var], tmp, i_cell)
    var_at_r = tmp(1)
  end subroutine iu_interpolate_scalar_at

  ! Interpolate multiple variables at location r
  subroutine iu_interpolate_at(ug, r, n_vars, i_vars, var_at_r, i_cell)
    type(iu_grid_t), intent(in) :: ug
    real(dp), intent(in)        :: r(3)             ! Interpolate at this point
    integer, intent(in)         :: n_vars           ! Number of variables to interpolate
    integer, intent(in)         :: i_vars(n_vars)   ! Indices of variables
    real(dp), intent(inout)     :: var_at_r(n_vars) ! Result at r
    ! On input: guess for nearby cell. Less than 1 means not set.
    ! On output: cell containing the point r.
    integer, intent(inout)      :: i_cell
    integer                     :: n
    real(dp)                    :: point_data(ug_max_points_per_cell, n_vars)

    if (i_cell > ug%n_cells) error stop "i_cell > ug%n_cells"

    call iu_get_cell(ug, r, i_cell)

    ! Exit if cell not found
    if (i_cell <= 0) return

    do n = 1, n_vars
       point_data(1:ug%n_points_per_cell, n) = &
            ug%point_data(ug%cells(:, i_cell), i_vars(n))
    end do

    select case (ug%cell_type)
    case (iu_triangle)
       call interpolate_triangle(n_vars, ug%cell_points(:, :, i_cell), &
            point_data, r, var_at_r)
    case (iu_quad)
       call interpolate_quad(n_vars, ug%cell_points(:, :, i_cell), &
            point_data, r, var_at_r)
    case (iu_tetra)
       call interpolate_tetrahedron(n_vars, ug%cell_points(:, :, i_cell), &
            point_data, r, var_at_r)
    case default
       error stop "Not implemented"
    end select

  end subroutine iu_interpolate_at

  subroutine interpolate_triangle(n_vars, points, point_data, r, res)
    integer, intent(in)   :: n_vars                ! Number of variables
    real(dp), intent(in)  :: points(3, 3)          ! vertices of the triangle
    ! point_data at vertices
    real(dp), intent(in)  :: point_data(ug_max_points_per_cell, n_vars)
    real(dp), intent(in)  :: r(3)                  ! point for interpolation
    real(dp), intent(out) :: res(n_vars)           ! interpolated value
    real(dp)              :: inv_area, areas(3)
    integer               :: n

    ! Compute areas to find barycentric coordinates. The factors 0.5 cancel in
    ! the final computation and have been left out.
    inv_area = 1/norm2(cross_product(points(:, 2) - points(:, 1), &
         points(:, 3) - points(:, 1)))
    areas(1) = norm2(cross_product(r - points(:, 2), r - points(:, 3)))
    areas(2) = norm2(cross_product(r - points(:, 3), r - points(:, 1)))
    areas(3) = norm2(cross_product(r - points(:, 1), r - points(:, 2)))

    ! Perform interpolation
    do n = 1, n_vars
       res(n) = sum(point_data(1:3, n) * areas) * inv_area
    end do
  end subroutine interpolate_triangle

  ! Based on https://www.cdsimpson.net/2014/10/barycentric-coordinates.html
  ! https://stackoverflow.com/questions/38545520/barycentric-coordinates-of-a-tetrahedron
  subroutine interpolate_tetrahedron(n_vars, points, point_data, r, res)
    integer, intent(in)    :: n_vars       ! Number of variables
    real(dp), intent(in)   :: points(3, 4) ! vertices of the quad
    ! point_data at vertices
    real(dp), intent(in)   :: point_data(ug_max_points_per_cell, n_vars)
    real(dp), intent(in)   :: r(3)         ! point for interpolation
    real(dp), intent(out)  :: res(n_vars)  ! interpolated value
    real(dp)               :: weights(4), inv_weight
    real(dp), dimension(3) :: v1r, v2r, v12, v13, v14, v23, v24
    integer                :: n

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
    inv_weight = 1/scalar_triple_product(v12, v13, v14)

    ! Perform interpolation
    do n = 1, n_vars
       res(n) = sum(point_data(1:4, n) * weights) * inv_weight
    end do

  end subroutine interpolate_tetrahedron

  ! Interpolate a quad element with four points at arbitary locations. The
  ! implementation is based on:
  ! https://www.reedbeta.com/blog/quadrilateral-interpolation-part-2/
  subroutine interpolate_quad(n_vars, points, point_data, r, res)
    integer, intent(in)   :: n_vars       ! Number of variables
    real(dp), intent(in)  :: points(3, 4) ! vertices of the triangle
    ! point_data at vertices
    real(dp), intent(in)  :: point_data(ug_max_points_per_cell, n_vars)
    real(dp), intent(in)  :: r(3)         ! point for interpolation
    real(dp), intent(out) :: res(n_vars)  ! interpolated value

    real(dp) :: q(3), b1(3), b2(3), b3(3), denom(3)
    real(dp) :: A, B, C, discrim, coeff(2), tmp(2)
    integer  :: max_dim(1), n
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
    do n = 1, n_vars
       tmp(1) = point_data(1, n) * (1 - coeff(1)) + point_data(2, n) * coeff(1)
       tmp(2) = point_data(4, n) * (1 - coeff(1)) + point_data(3, n) * coeff(1)
       res(n) = tmp(1) * (1 - coeff(2)) + tmp(2) * coeff(2)
    end do

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

  subroutine iu_read_grid(filename, ug, coord_scale_factor)
    character(len=*), intent(in) :: filename
    type(iu_grid_t), intent(out) :: ug
    ! Scale coordinates by this factor
    real(dp), intent(in), optional :: coord_scale_factor
    type(binda_t)                :: bfile
    integer                      :: ix, i_point_data, i_cell_data

    call binda_open_file(filename, bfile)
    call binda_read_header(bfile)

    call binda_get_index(bfile, "cells", ix)
    if (ix == -1) error stop "cells not found in binda file"
    call binda_read_alloc_int32_2d(bfile, ix, ug%cells)

    select case (bfile%metadata(ix))
    case ("triangle")
       ug%cell_type = iu_triangle
    case ("quad")
       ug%cell_type = iu_quad
    case ("tetra")
       ug%cell_type = iu_tetra
    case default
       write(error_unit, *) "Cell type '", bfile%metadata(ix), &
            "' not supported"
       error stop "Unsupported cell type"
    end select

    call binda_get_index(bfile, "points", ix)
    if (ix == -1) error stop "points not found in binda file"
    call binda_read_alloc_float64_2d(bfile, ix, ug%points)

    call binda_get_index(bfile, "cell_neighbors", ix)
    if (ix == -1) error stop "cell_neighbors not found in binda file"
    call binda_read_alloc_int32_2d(bfile, ix, ug%neighbors)

    if (present(coord_scale_factor)) then
       ug%points(:, :) = ug%points(:, :) * coord_scale_factor
    end if

    ! Store information about mesh
    ug%n_cells = size(ug%cells, 2)
    ug%n_points_per_cell = size(ug%cells, 1)
    ug%n_faces_per_cell = size(ug%cells, 1) ! works for triangle, quad, tetra
    ug%n_points = size(ug%points, 2)
    ug%rmin = minval(ug%points, dim=2)
    ug%rmax = maxval(ug%points, dim=2)

    ! Convert to 1-based indexing
    ug%cells = ug%cells + 1
    ug%neighbors = ug%neighbors + 1

    ! Allocate storage for point data
    ug%n_point_data = count(bfile%name == "point_data")
    allocate(ug%point_data(ug%n_points, ug%n_point_data))
    allocate(ug%point_data_names(ug%n_point_data))

    ! Allocate storage for cell data
    ug%n_cell_data = count(bfile%name == "cell_data")
    allocate(ug%cell_data(ug%n_cells, ug%n_cell_data))
    allocate(ug%cell_data_names(ug%n_cell_data))

    ! Read point and cell data
    i_point_data = 0
    i_cell_data = 0

    do ix = 1, bfile%n_entries
       if (bfile%name(ix) == "point_data") then
          i_point_data = i_point_data + 1
          ug%point_data_names(i_point_data) = bfile%metadata(ix)
          call binda_read_float64_1d(bfile, ix, &
               ug%n_points, ug%point_data(:, i_point_data))
       else if (bfile%name(ix) == "cell_data") then
          i_cell_data = i_cell_data + 1
          ug%cell_data_names(i_cell_data) = bfile%metadata(ix)
          call binda_read_float64_1d(bfile, ix, &
               ug%n_cells, ug%cell_data(:, i_cell_data))
       end if
    end do

    call binda_close_file(bfile)

    ! Store points for each cell
    call set_cell_points(ug)

    ! Store normal vectors for each cell face
    call set_face_normal_vectors(ug)

    ! Build kd-tree for efficient lookup of new points
    call build_kdtree(ug)

  end subroutine iu_read_grid

  !> Write unstructured grid back to VTK file
  subroutine iu_write_vtk(ug, filename)
    type(iu_grid_t), intent(in)  :: ug
    character(len=*), intent(in) :: filename !< Filename for the vtk file
    integer                      :: n
    integer, allocatable         :: offsets(:), cell_types(:), connectivity(:)
    type(vtk_t)                  :: vtkf

    allocate(offsets(ug%n_cells))
    allocate(cell_types(ug%n_cells))
    allocate(connectivity(ug%n_cells*ug%n_points_per_cell))

    select case (ug%cell_type)
    case (iu_triangle)
       cell_types(:) = 5
    case (iu_quad)
       cell_types(:) = 9
    case (iu_tetra)
       cell_types(:) = 10
    case default
       error stop "Unsupported ug%cell_type"
    end select

    do n = 1, ug%n_cells
       offsets(n) = n * ug%n_points_per_cell
    end do

    ! Revert back to zero-based indexing
    connectivity(:) = pack(ug%cells-1, .true.)

    call vtk_ini_xml(vtkf, trim(filename), 'UnstructuredGrid')
    call vtk_dat_xml(vtkf, "UnstructuredGrid", .true.)
    call vtk_unstr_geo_xml(vtkf, ug%n_points, ug%n_cells, ug%points)
    call vtk_unstr_con_xml(vtkf, connectivity, offsets, cell_types, ug%n_cells)

    call vtk_dat_xml(vtkf, "CellData", .true.)
    do n = 1, ug%n_cell_data
       call vtk_var_r8_xml(vtkf, trim(ug%cell_data_names(n)), &
            ug%cell_data(:, n), ug%n_cells)
    end do
    call vtk_dat_xml(vtkf, "CellData", .false.)

    call vtk_dat_xml(vtkf, "PointData", .true.)
    do n = 1, ug%n_point_data
       call vtk_var_r8_xml(vtkf, trim(ug%point_data_names(n)), &
            ug%point_data(:, n), ug%n_points)
    end do
    call vtk_dat_xml(vtkf, "PointData", .false.)

    call vtk_unstr_geo_xml_close(vtkf)
    call vtk_dat_xml(vtkf, "UnstructuredGrid", .false.)
    call vtk_end_xml(vtkf)
  end subroutine iu_write_vtk

  !> Trace the field given by the point data in i_field(:) until a boundary is reached
  subroutine iu_trace_field(ug, ndim, r_start, i_field, max_points, n_points, &
       points, fields, min_dx, max_dx, abs_tol)
    type(iu_grid_t), intent(in) :: ug
    integer, intent(in)         :: ndim !< Problem dimension
    real(dp), intent(in)        :: r_start(ndim) !< Start location
    integer, intent(in)         :: i_field(ndim) !< Field to trace (stored as point data)
    integer, intent(in)         :: max_points !< Maximum number of points along path
    integer, intent(out)        :: n_points !< Number of points along path
    real(dp), intent(inout)     :: points(ndim, max_points) !< Location at each point
    real(dp), intent(inout)     :: fields(ndim, max_points) !< Field at each point
    real(dp), intent(in)        :: min_dx !< Min distance per step
    real(dp), intent(in)        :: max_dx !< Max distance per step
    real(dp), intent(in)        :: abs_tol !< Absolute tolerance in r per step

    integer  :: i_cell
    real(dp) :: r0(3), rq(3), r1(3)
    real(dp) :: field_0(ndim), field_1(ndim)
    real(dp) :: unitvec_0(ndim), unitvec_1(ndim)
    real(dp) :: dx, dx_factor, error_estimate

    if (max_dx < min_dx) error stop "max_dx < min_dx"
    if (max_points < 1) error stop "max_points < 1"
    if (abs_tol <= 0.0_dp) error stop "abs_tol <= 0.0_dp"

    ! Set unused coordinates to zero
    r0(ndim+1:) = 0.0_dp
    rq(ndim+1:) = 0.0_dp
    r1(ndim+1:) = 0.0_dp

    r0(1:ndim) = r_start
    dx         = max_dx
    n_points   = 1
    i_cell     = 0

    ! Get field at initial point
    call iu_interpolate_at(ug, r0, ndim, i_field, field_0, i_cell)
    if (i_cell <= 0) return

    points(:, n_points) = r0(1:ndim)
    fields(:, n_points) = field_0

    do while (n_points <= max_points)
       ! Field for forward Euler step
       call iu_interpolate_at(ug, r0, ndim, i_field, field_0, i_cell)
       unitvec_0 = field_0 / norm2(field_0)

       ! Get field for Heun's method, while reducing step size if a domain
       ! boundary is reached
       do
          rq(1:ndim) = r0(1:ndim) + dx * unitvec_0
          call get_field_unitvec(ndim, rq, unitvec_1, i_cell)

          call iu_interpolate_at(ug, rq, ndim, i_field, field_1, i_cell)
          unitvec_1 = field_1 / norm2(field_1)

          ! Check if still inside domain
          if (i_cell > 0) then
             exit               ! Found field
          else if (dx >= 2 * min_dx) then
             dx = 0.5_dp * dx   ! Reduce dx
          else
             return             ! Reached boundary with small dx
          end if
       end do

       ! New position with Heun's method
       r1(1:ndim) = r0(1:ndim) + 0.5_dp * dx * (unitvec_0 + unitvec_1)

       ! Estimate error in coordinates
       error_estimate = norm2(r1(1:ndim) - rq(1:ndim))

       if (error_estimate < abs_tol) then
          ! Step is accepted
          n_points = n_points + 1
          points(:, n_points) = r1(1:ndim)
          fields(:, n_points) = 0.5_dp * (field_0 + field_1)
          r0(1:ndim) = r1(1:ndim)
       end if

       ! Adjust dx to have approximately an error of 0.5 * abs_tol per step
       dx_factor = min(2.0_dp, (0.5_dp*abs_tol/error_estimate)**(1/3.0_dp))
       dx = max(min_dx, min(max_dx, dx * dx_factor))
    end do

    ! If this is reached, the loop did not exit at a boundary
    n_points = max_points + 1

  contains

    subroutine get_field_unitvec(ndim, r, unitvec, i_cell)
      integer, intent(in)    :: ndim
      real(dp), intent(in)   :: r(ndim)
      real(dp), intent(out)  :: unitvec(ndim)
      integer, intent(inout) :: i_cell

      call iu_interpolate_at(ug, r, ndim, i_field, unitvec, i_cell)
      unitvec = unitvec / norm2(unitvec)
    end subroutine get_field_unitvec

  end subroutine iu_trace_field

end module m_interp_unstructured
