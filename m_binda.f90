! Module for reading BInary N-dimensional DAta (binda)
module m_binda
  use iso_fortran_env, only: error_unit, int64

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)
  integer, parameter :: sp = kind(0.0e0)

  ! Type for reading binary input files
  type binda_t
     integer                         :: file_unit = -1
     integer                         :: n_entries
     integer(int64)                  :: total_header_size
     character(len=128), allocatable :: name(:)
     character(len=128), allocatable :: dtype(:)
     character(len=128), allocatable :: metadata(:)
     integer(int64), allocatable     :: ndim(:)
     integer(int64), allocatable     :: dshape(:, :)
     integer(int64), allocatable     :: offset(:)
  end type binda_t

  ! Public types
  public :: binda_t

  ! Public methods
  public :: binda_open_file
  public :: binda_close_file
  public :: binda_get_index
  public :: binda_read_header
  public :: binda_read_float64_1d
  public :: binda_read_int32_1d
  public :: binda_read_alloc_float64_2d
  public :: binda_read_alloc_int32_2d

contains

  subroutine binda_open_file(filename, bfile)
    character(len=*), intent(in) :: filename
    type(binda_t), intent(inout) :: bfile

    open(newunit=bfile%file_unit, file=filename, status="old", &
         form="unformatted", access="stream")
  end subroutine binda_open_file

  subroutine binda_close_file(bfile)
    type(binda_t), intent(inout) :: bfile
    close(bfile%file_unit)
  end subroutine binda_close_file

  subroutine binda_read_header(bfile)
    type(binda_t), intent(inout) :: bfile
    integer                      :: i
    character(len=8)             :: identifier
    integer(int64)               :: n_entries_int64

    read(bfile%file_unit) identifier
    if (identifier /= "BINDA") error stop "Wrong file format"

    ! Read the number of entries
    read(bfile%file_unit) n_entries_int64

    ! Convert to default integer
    bfile%n_entries = int(n_entries_int64)

    ! Read total header size
    read(bfile%file_unit) bfile%total_header_size

    allocate(bfile%name(bfile%n_entries))
    allocate(bfile%dtype(bfile%n_entries))
    allocate(bfile%metadata(bfile%n_entries))
    allocate(bfile%ndim(bfile%n_entries))
    allocate(bfile%dshape(8, bfile%n_entries))
    allocate(bfile%offset(bfile%n_entries))

    ! Loop over each entry
    do i = 1, bfile%n_entries
       ! Read entry metadata
       read(bfile%file_unit) bfile%name(i)
       read(bfile%file_unit) bfile%dtype(i)
       read(bfile%file_unit) bfile%metadata(i)
       read(bfile%file_unit) bfile%ndim(i)
       read(bfile%file_unit) bfile%dshape(:, i)
       read(bfile%file_unit) bfile%offset(i)
    end do

  end subroutine binda_read_header

  subroutine binda_read_alloc_int32_2d(bfile, ix, array)
    type(binda_t), intent(in)           :: bfile
    integer, intent(in)                 :: ix
    integer, allocatable, intent(inout) :: array(:, :)
    integer(int64), allocatable         :: int64_array(:, :)

    if (bfile%ndim(ix) /= 2) then
       write(error_unit, *) "Found ndim = ", bfile%ndim(ix), " expecting 2"
       error stop "Invalid ndim"
    end if

    select case (bfile%dtype(ix))
    case ("int32")
       allocate(array(bfile%dshape(2, ix), bfile%dshape(1, ix)))
       read(bfile%file_unit, pos=bfile%offset(ix)+1) array
    case ("int64")
       allocate(int64_array(bfile%dshape(2, ix), bfile%dshape(1, ix)))
       read(bfile%file_unit, pos=bfile%offset(ix)+1) int64_array
       array = int(int64_array)
    case default
       write(error_unit, *) " Found dtype: ", trim(bfile%dtype(ix))
       error stop "Unsupported data type"
    end select
  end subroutine binda_read_alloc_int32_2d

  subroutine binda_read_alloc_float64_2d(bfile, ix, array)
    type(binda_t), intent(in)            :: bfile
    integer, intent(in)                  :: ix
    real(dp), allocatable, intent(inout) :: array(:, :)
    real(sp), allocatable                :: float32_array(:, :)

    if (bfile%ndim(ix) /= 2) error stop "Invalid ndim"

    select case (bfile%dtype(ix))
    case ("float64")
       allocate(array(bfile%dshape(2, ix), bfile%dshape(1, ix)))
       read(bfile%file_unit, pos=bfile%offset(ix)+1) array
    case ("float32")
       allocate(float32_array(bfile%dshape(2, ix), bfile%dshape(1, ix)))
       read(bfile%file_unit, pos=bfile%offset(ix)+1) float32_array
       array = float32_array
    case default
       write(error_unit, *) " Found dtype: ", trim(bfile%dtype(ix))
       error stop "Unsupported data type"
    end select
  end subroutine binda_read_alloc_float64_2d

  subroutine binda_read_float64_1d(bfile, ix, array_size, array)
    type(binda_t), intent(in) :: bfile
    integer, intent(in)       :: ix
    integer, intent(in)       :: array_size
    real(dp), intent(inout)   :: array(array_size)
    real(sp), allocatable     :: float32_array(:)

    if (bfile%ndim(ix) /= 1) error stop "Invalid ndim"
    if (bfile%dshape(1, ix) /= array_size) error stop "Invalid array_size"

    select case (bfile%dtype(ix))
    case ("float64")
       read(bfile%file_unit, pos=bfile%offset(ix)+1) array
    case ("float32")
       allocate(float32_array(array_size))
       read(bfile%file_unit, pos=bfile%offset(ix)+1) float32_array
       array = float32_array
    case default
       write(error_unit, *) " Found dtype: ", trim(bfile%dtype(ix))
       error stop "Unsupported data type"
    end select
  end subroutine binda_read_float64_1d

  subroutine binda_read_int32_1d(bfile, ix, array_size, array)
    type(binda_t), intent(in)   :: bfile
    integer, intent(in)         :: ix
    integer, intent(in)         :: array_size
    integer, intent(inout)      :: array(array_size)
    integer(int64), allocatable :: int64_array(:)

    if (bfile%ndim(ix) /= 1) error stop "Invalid ndim"
    if (bfile%dshape(1, ix) /= array_size) error stop "Invalid array_size"

    select case (bfile%dtype(ix))
    case ("int32")
       read(bfile%file_unit, pos=bfile%offset(ix)+1) array
    case ("int64")
       allocate(int64_array(array_size))
       read(bfile%file_unit, pos=bfile%offset(ix)+1) int64_array
       array = int(int64_array)
    case default
       write(error_unit, *) " Found dtype: ", trim(bfile%dtype(ix))
       error stop "Unsupported data type"
    end select
  end subroutine binda_read_int32_1d

  ! Return first index where name matches, and -1 if there is no match
  subroutine binda_get_index(bfile, name, ix)
    type(binda_t), intent(in)    :: bfile
    character(len=*), intent(in) :: name
    integer, intent(out)         :: ix

    do ix = 1, bfile%n_entries
       if (bfile%name(ix) == name) exit
    end do

    if (ix == bfile%n_entries + 1) ix = -1
  end subroutine binda_get_index

end module m_binda
