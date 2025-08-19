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
     integer :: n_icell_data

     real(dp) :: rmin(3)
     real(dp) :: rmax(3)

     real(dp), allocatable :: points(:, :)
     logical, allocatable  :: point_is_at_boundary(:)
     integer, allocatable  :: cells(:, :)
     real(dp), allocatable :: cell_points(:, :, :)
     real(dp), allocatable :: cell_volume(:)
     integer, allocatable  :: neighbors(:, :)
     real(dp), allocatable :: cell_face_normals(:, :, :)
     real(dp), allocatable :: point_data(:, :)
     real(dp), allocatable :: cell_data(:, :)
     integer, allocatable  :: icell_data(:, :)

     ! Names of the point data variables
     character(len=128), allocatable :: point_data_names(:)

     ! Names of the cell data variables
     character(len=128), allocatable :: cell_data_names(:)

     ! Names of the integer cell data variables
     character(len=128), allocatable :: icell_data_names(:)

     ! kd-tree for searching a nearby cell
     type(kdtree2) :: tree
  end type iu_grid_t

  interface
     subroutine integrate_sub_t(ndim, r, field, nvar, integrand)
       import
       integer, intent(in)   :: ndim
       real(dp), intent(in)  :: r(ndim)
       real(dp), intent(in)  :: field(ndim)
       integer, intent(in)   :: nvar
       real(dp), intent(out) :: integrand(nvar)
     end subroutine integrate_sub_t
  end interface

  ! Public types
  public :: iu_grid_t
  public :: iu_triangle, iu_quad, iu_tetra
  public :: iu_ndim_cell_type

  ! Public methods
  public :: iu_read_grid
  public :: iu_get_point_data_index
  public :: iu_get_cell_data_index
  public :: iu_get_icell_data_index
  public :: iu_add_cell_data
  public :: iu_add_icell_data
  public :: iu_add_point_data
  public :: iu_get_cell_center
  public :: iu_get_cell
  public :: iu_get_cell_scalar_at
  public :: iu_get_icell_scalar_at
  public :: iu_integrate_along_field
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

  !> Extend cell data with one variable
  subroutine iu_add_icell_data(ug, name, i_var)
    type(iu_grid_t), intent(inout)  :: ug
    character(len=*), intent(in)    :: name
    integer, intent(out)            :: i_var
    integer, allocatable            :: old_data(:, :)
    character(len=128), allocatable :: old_names(:)

    call move_alloc(ug%icell_data, old_data)
    call move_alloc(ug%icell_data_names, old_names)
    allocate(ug%icell_data(ug%n_cells, ug%n_icell_data+1))
    allocate(ug%icell_data_names(ug%n_icell_data+1))

    ug%icell_data(:, 1:ug%n_icell_data) = old_data
    ug%icell_data_names(1:ug%n_icell_data) = old_names

    ug%n_icell_data = ug%n_icell_data + 1
    ug%icell_data_names(ug%n_icell_data) = name
    i_var = ug%n_icell_data
  end subroutine iu_add_icell_data

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

  ! Find index of cell data variable, -1 if not present
  subroutine iu_get_icell_data_index(ug, name, ix)
    type(iu_grid_t), intent(in)  :: ug
    character(len=*), intent(in) :: name
    integer, intent(out)         :: ix

    do ix = 1, ug%n_icell_data
       if (ug%icell_data_names(ix) == name) exit
    end do

    if (ix == ug%n_icell_data + 1) ix = -1
  end subroutine iu_get_icell_data_index

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

  ! Compute normal vectors to cell faces. Also store which points lie at the
  ! domain boundary.
  subroutine set_face_normal_vectors(ug)
    type(iu_grid_t), intent(inout) :: ug
    integer                        :: n, k, k1, k2
    real(dp)                       :: vec(3), center(3)
    real(dp)                       :: normal_cell(3), normal_face(3)

    allocate(ug%cell_face_normals(3, ug%n_points_per_cell, ug%n_cells))
    allocate(ug%point_is_at_boundary(ug%n_points))
    ug%point_is_at_boundary(:) = .false.

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

             if (ug%neighbors(k, n) < 1) then
                ug%point_is_at_boundary(ug%cells([k, k1], n)) = .true.
             end if
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

             if (ug%neighbors(k, n) < 1) then
                ug%point_is_at_boundary(ug%cells([k, k1, k2], n)) = .true.
             end if
          end do
       end do
    case default
       error stop "Not implemented"
    end select

  end subroutine set_face_normal_vectors

  subroutine set_cell_volumes(ug)
    type(iu_grid_t), intent(inout) :: ug
    integer                        :: n
    real(dp)                       :: v12(3), v13(3), v14(3)
    real(dp)                       :: area1, area2
    real(dp), parameter            :: sixth = 1/6.0_dp

    allocate(ug%cell_volume(ug%n_cells))

    select case (ug%cell_type)
    case (iu_triangle)
       do n = 1, ug%n_cells
          associate (p => ug%cell_points(:, :, n))
            ug%cell_volume(n) = 0.5_dp * norm2(&
                 cross_product(p(:, 2) - p(:, 1), p(:, 3) - p(:, 1)))
          end associate
       end do
    case (iu_quad)
       do n = 1, ug%n_cells
          associate (p => ug%cell_points(:, :, n))
            ! Split into two triangles: (p1, p2, p3) and (p1, p3, p4)
            area1 = 0.5_dp * norm2(cross_product(p(:, 2) - p(:, 1), &
                 p(:, 3) - p(:, 1)))
            area2 = 0.5_dp * norm2(cross_product(p(:, 3) - p(:, 1), &
                 p(:, 4) - p(:, 1)))
            ug%cell_volume(n) = area1 + area2
          end associate
       end do
    case (iu_tetra)
       do n = 1, ug%n_cells
          associate (p => ug%cell_points(:, :, n))
            v12 = p(:, 2) - p(:, 1)
            v13 = p(:, 3) - p(:, 1)
            v14 = p(:, 4) - p(:, 1)
            ug%cell_volume(n) = sixth * scalar_triple_product(v12, v13, v14)
          end associate
       end do
    end select
  end subroutine set_cell_volumes

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

  ! Get the value of integer cell data at the cell that contains r
  subroutine iu_get_icell_scalar_at(ug, r, i_var, i_guess, res)
    type(iu_grid_t), intent(in) :: ug
    real(dp), intent(in)        :: r(3)  ! Location
    integer, intent(in)         :: i_var ! Index of cell data variable
    ! On input: guess for nearby cell. Less than 1 means not set.
    ! On output: cell containing the point r.
    integer, intent(inout)      :: i_guess
    integer, intent(inout)      :: res    ! Result

    call iu_get_cell(ug, r, i_guess)
    if (i_guess >= 1) res = ug%icell_data(i_guess, i_var)
  end subroutine iu_get_icell_scalar_at

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
            point_data, ug%cell_volume(i_cell), r, var_at_r)
    case (iu_quad)
       call interpolate_quad(n_vars, ug%cell_points(:, :, i_cell), &
            point_data, r, var_at_r)
    case (iu_tetra)
       call interpolate_tetrahedron(n_vars, ug%cell_points(:, :, i_cell), &
            point_data, ug%cell_volume(i_cell), r, var_at_r)
    case default
       error stop "Not implemented"
    end select

  end subroutine iu_interpolate_at

  subroutine interpolate_triangle(n_vars, points, point_data, area, r, res)
    integer, intent(in)   :: n_vars                ! Number of variables
    real(dp), intent(in)  :: points(3, 3)          ! vertices of the triangle
    ! point_data at vertices
    real(dp), intent(in)  :: point_data(ug_max_points_per_cell, n_vars)
    real(dp), intent(in)  :: r(3)                  ! point for interpolation
    real(dp), intent(in)  :: area                  ! Area of the cell
    real(dp), intent(out) :: res(n_vars)           ! interpolated value
    real(dp)              :: inv_area, areas(3)
    integer               :: n

    ! Compute areas to find barycentric coordinates. The factors 0.5 cancel in
    ! the final computation and have been left out.
    areas(1) = 0.5_dp * norm2(cross_product(r - points(:, 2), r - points(:, 3)))
    areas(2) = 0.5_dp * norm2(cross_product(r - points(:, 3), r - points(:, 1)))
    areas(3) = 0.5_dp * norm2(cross_product(r - points(:, 1), r - points(:, 2)))

    ! Perform interpolation
    inv_area = 1/area
    do n = 1, n_vars
       res(n) = sum(point_data(1:3, n) * areas) * inv_area
    end do
  end subroutine interpolate_triangle

  ! Based on https://www.cdsimpson.net/2014/10/barycentric-coordinates.html
  ! https://stackoverflow.com/questions/38545520/barycentric-coordinates-of-a-tetrahedron
  subroutine interpolate_tetrahedron(n_vars, points, point_data, vol, r, res)
    integer, intent(in)    :: n_vars       ! Number of variables
    real(dp), intent(in)   :: points(3, 4) ! Vertices of the tetrahedron
    ! point_data at vertices
    real(dp), intent(in)   :: point_data(ug_max_points_per_cell, n_vars)
    real(dp), intent(in)   :: vol          ! Volume of the cell
    real(dp), intent(in)   :: r(3)         ! Point for interpolation
    real(dp), intent(out)  :: res(n_vars)  ! Interpolated value
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

    ! Perform interpolation
    inv_weight = 1/(6 * vol)
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
    integer                      :: ix, i_point_data, i_cell_data, i_icell_data

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

    ! Allocate storage for integer cell data
    ug%n_icell_data = count(bfile%name == "icell_data")
    allocate(ug%icell_data(ug%n_cells, ug%n_icell_data))
    allocate(ug%icell_data_names(ug%n_icell_data))

    ! Read point and cell data
    i_point_data = 0
    i_cell_data = 0
    i_icell_data = 0

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
       else if (bfile%name(ix) == "icell_data") then
          i_icell_data = i_icell_data + 1
          ug%icell_data_names(i_icell_data) = bfile%metadata(ix)
          call binda_read_int32_1d(bfile, ix, &
               ug%n_cells, ug%icell_data(:, i_icell_data))
       end if
    end do

    call binda_close_file(bfile)

    ! Store points for each cell
    call set_cell_points(ug)

    ! Store normal vectors for each cell face
    call set_face_normal_vectors(ug)

    ! Store volume of each cell
    call set_cell_volumes(ug)

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
    do n = 1, ug%n_icell_data
       call vtk_var_i4_xml(vtkf, trim(ug%icell_data_names(n)), &
            ug%icell_data(:, n), ug%n_cells)
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

  !> Integrate a function along the field given by the point data in
  !> i_field(:) until a boundary is reached
  subroutine iu_integrate_along_field(ug, ndim, sub_int, r_start, i_field, &
       min_dx, max_dx, max_steps, rtol, atol, reverse, nvar, y, y_field, &
       n_steps, axisymmetric, i_icell_mask, mask_value, boundary_material)
    type(iu_grid_t), intent(in)   :: ug
    procedure(integrate_sub_t)    :: sub_int
    integer, intent(in)           :: ndim !< Number of spatial dimensions
    real(dp), intent(in)          :: r_start(ndim) !< Start location
    integer, intent(in)           :: i_field(ndim) !< Field to trace (stored as point data)
    real(dp), intent(in)          :: min_dx    !< Min distance per step
    real(dp), intent(in)          :: max_dx    !< Max distance per step
    integer, intent(in)           :: max_steps !< Max number of steps
    real(dp), intent(in)          :: rtol    !< Relative tolerance
    real(dp), intent(in)          :: atol    !< Absolute tolerance
    logical, intent(in)           :: reverse !< Go in minus field direction
    integer, intent(in)           :: nvar !< Number of variables to integrate
    !> Solution curve. y(:, i) contains position (ndim) and the nvar
    !> solution variables
    real(dp), intent(inout)       :: y(ndim+nvar, max_steps)
    !> Field along the solution curve
    real(dp), intent(inout)       :: y_field(ndim, max_steps)
    !> Number of steps taken plus one. Equal to max_steps+1 if the boundary
    !> was not reached within max_steps.
    integer, intent(out)          :: n_steps
    !> Whether the domain is axisymmetric
    logical, intent(in)           :: axisymmetric
    !> Use integer cell mask (should be non-negative)
    integer, intent(in), optional :: i_icell_mask
    !< Integrate only in region where mask has this value
    integer, intent(in), optional :: mask_value
    !> What kind of boundary was reached. -1 indicates a physical boundary,
    !> otherwise the cell mask of the boundary stored.
    integer, intent(out), optional :: boundary_material

    real(dp), parameter :: safety_fac = 0.8_dp
    real(dp), parameter :: inv_24 = 1/24.0_dp
    real(dp), parameter :: inv_9 = 1/9.0_dp
    real(dp), parameter :: min_radius = 1e-12_dp

    integer  :: i_cell, i_cell_prev, iteration, last_rejected
    real(dp) :: y_new(ndim+nvar), y_2nd(ndim+nvar)
    real(dp) :: r0(3), r(3), field(ndim)
    real(dp) :: err, scales(ndim+nvar), dx, dx_factor
    real(dp) :: k(ndim+nvar, 4), max_growth
    logical  :: invalid_position

    if (max_dx < min_dx) error stop "max_dx < min_dx"
    if (max_steps < 1) error stop "max_steps < 1"
    if (present(i_icell_mask) .neqv. present(mask_value)) &
         error stop "present(i_icell_mask) .neqv. present(mask_value)"

    ! Set unused coordinates to zero
    r0(ndim+1:) = 0.0_dp
    r(ndim+1:) = 0.0_dp

    ! Initial solution
    n_steps = 1
    y(1:ndim, n_steps) = r_start
    y(ndim+1:ndim+nvar, n_steps) = 0.0_dp

    ! Initialization
    r0(1:ndim)       = r_start
    invalid_position = .false.
    last_rejected    = -100
    dx               = max_dx
    i_cell           = 0

    ! Get field at initial point
    call iu_interpolate_at(ug, r0, ndim, i_field, field, i_cell)

    ! Exit if initial cell is not valid
    if (.not. cell_is_valid(i_cell, mask_value, i_icell_mask)) then
       ! Optionally set boundary type and then exit
       if (present(boundary_material)) then
          if (i_cell <= 0) then
             boundary_material = -1
          else
             boundary_material = ug%icell_data(i_cell, i_icell_mask)
          end if
       end if

       return
    end if

    ! Initialization
    y_field(:, n_steps) = field
    i_cell_prev = i_cell

    ! The code below implements the Bogacki–Shampine Runge-Kutta method, which
    ! is third order accurate. Higher order might not be beneficial since we
    ! are linearly interpolating fields.
    do iteration = 1, huge(1)-1

       ! Handle cases in which the previous position was invalid
       if (invalid_position) then
          if (dx >= 2 * min_dx) then
             ! Reduce step size
             last_rejected = iteration - 1
             dx = 0.1_dp * dx
          else
             ! Optionally set boundary type and then exit
             if (present(boundary_material)) then
                if (i_cell <= 0) then
                   boundary_material = -1
                else
                   boundary_material = ug%icell_data(i_cell, i_icell_mask)
                end if
             end if
             return
          end if
       end if

       invalid_position = .false.
       i_cell = i_cell_prev
       r0(1:ndim) = y(1:ndim, n_steps)   ! Current position

       ! First sub-step, re-uses field from last step
       field = y_field(:, n_steps)

       k(1:ndim, 1) = get_unitvec(ndim, field, reverse)
       call sub_int(ndim, r0(1:ndim), field, nvar, k(ndim+1:ndim+nvar, 1))

       ! Second sub-step
       r(1:ndim) = r0(1:ndim) + 0.5_dp * dx * k(1:ndim, 1)
       if (axisymmetric) r(1) = max(r(1), min_radius)

       call iu_interpolate_at(ug, r, ndim, i_field, field, i_cell)
       if (.not. cell_is_valid(i_cell, mask_value, i_icell_mask)) then
          invalid_position = .true.
          cycle
       end if

       k(1:ndim, 2) = get_unitvec(ndim, field, reverse)
       call sub_int(ndim, r(1:ndim), field, nvar, k(ndim+1:ndim+nvar, 2))

       ! Third sub-step
       r(1:ndim) = r0(1:ndim) + 0.75_dp * dx * k(1:ndim, 2)
       if (axisymmetric) r(1) = max(r(1), min_radius)

       call iu_interpolate_at(ug, r, ndim, i_field, field, i_cell)
       if (.not. cell_is_valid(i_cell, mask_value, i_icell_mask)) then
          invalid_position = .true.
          cycle
       end if

       k(1:ndim, 3) = get_unitvec(ndim, field, reverse)
       call sub_int(ndim, r(1:ndim), field, nvar, k(ndim+1:ndim+nvar, 3))

       ! Update with third-order scheme
       y_new = y(:, n_steps) + dx * inv_9 * &
            (2 * k(:, 1) + 3 * k(:, 2) + 4 * k(:, 3))

       ! Fourth sub-step
       r(1:ndim) = y_new(1:ndim)
       if (axisymmetric) r(1) = max(r(1), min_radius)

       call iu_interpolate_at(ug, r, ndim, i_field, field, i_cell)
       if (.not. cell_is_valid(i_cell, mask_value, i_icell_mask)) then
          invalid_position = .true.
          cycle
       end if

       k(1:ndim, 4) = get_unitvec(ndim, field, reverse)
       call sub_int(ndim, r(1:ndim), field, nvar, k(ndim+1:ndim+nvar, 4))

       ! Estimate with second-order scheme
       y_2nd = y(:, n_steps) + dx * inv_24 * &
            (7 * k(:, 1) + 6 * k(:, 2) + 8 * k(:, 3) + 3 * k(:, 4))

       scales = atol + max(abs(y_new), abs(y_2nd)) * rtol
       err = sqrt(sum(((y_new-y_2nd)/scales)**2) / 3)

       if (err <= 1 .or. dx < 2 * min_dx) then
          ! Step is accepted, advance solution y
          n_steps = n_steps + 1
          if (n_steps > max_steps) return

          y(:, n_steps) = y_new
          if (axisymmetric) y(1, n_steps) = max(y(1, n_steps), min_radius)
          y_field(:, n_steps) = field
          i_cell_prev = i_cell
       else
          last_rejected = iteration
       end if

       if (last_rejected > iteration - 2) then
          ! If steps were recently rejected, do not grow dx
          max_growth = 1.0_dp
       else
          ! Limit increase, to prevent increasing dx at boundaries
          max_growth = 2.0_dp
       end if

       ! Adjust dx to for an error of safety_fac * abs_tol per step
       dx_factor = min(max_growth, safety_fac * (1/err)**(1/3.0_dp))
       dx = max(min_dx, min(max_dx, dx * dx_factor))
    end do

  contains

    function get_unitvec(ndim, field, reverse) result(unitvec)
      integer, intent(in)  :: ndim
      real(dp), intent(in) :: field(ndim)
      logical, intent(in)  :: reverse
      real(dp)             :: unitvec(ndim)

      unitvec = field / norm2(field)
      if (reverse) unitvec = -unitvec
    end function get_unitvec

    logical function cell_is_valid(ix, mask_value, i_icell_mask)
      integer, intent(in)           :: ix
      integer, intent(in), optional :: mask_value
      integer, intent(in), optional :: i_icell_mask

      if (ix <= 0) then
         cell_is_valid = .false.
      else if (present(mask_value)) then
         cell_is_valid = (ug%icell_data(ix, i_icell_mask) == mask_value)
      else
         cell_is_valid = .true.
      end if
    end function cell_is_valid

  end subroutine iu_integrate_along_field

end module m_interp_unstructured
