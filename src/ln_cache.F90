!!--------------------------------------------------------------------------!
!! libNEGF: a general library for Non-Equilibrium Green's functions.        !
!! Copyright (C) 2012                                                       !
!!                                                                          !
!! This file is part of libNEGF: a library for                              !
!! Non Equilibrium Green's Function calculation                             !
!!                                                                          !
!! Developers: Alessandro Pecchia, Gabriele Penazzi                         !
!! Former Conctributors: Luca Latessa, Aldo Di Carlo                        !
!!                                                                          !
!! libNEGF is free software: you can redistribute it and/or modify          !
!! it under the terms of the GNU Lesse General Public License as published  !
!! by the Free Software Foundation, either version 3 of the License, or     !
!! (at your option) any later version.                                      !
!!                                                                          !
!!  You should have received a copy of the GNU Lesser General Public        !
!!  License along with libNEGF.  If not, see                                !
!!  <http://www.gnu.org/licenses/>.                                         !
!!--------------------------------------------------------------------------!

module ln_cache

  use mat_def, only: z_DNS, assignment(=)
  use ln_precision
  use globals
  use outmatrix, only: outmat_c, inmat_c

  implicit none
  private

  public :: TMatrixCache
  public :: TMatrixCacheMem, TMatrixCacheDisk
  public :: TMatrixCacheDummy
  public :: TMatLabel

  type TMatLabel
    integer :: kpoint = -1
    integer :: energy_point = -1
    integer :: spin = -1
    integer :: row_block = -1
    integer :: col_block = -1
  end type TMatLabel

  interface operator (==)
    module procedure labeleq
  end interface

  type, abstract :: TMatrixCache
  contains
    procedure(abst_add_matrix_to_cache), deferred :: add
    procedure(abst_retrieve_matrix_from_cache), deferred :: retrieve
    procedure(abst_destroy_matrix_cache), deferred :: destroy
    procedure(abst_is_cached), deferred :: is_cached
  end type TMatrixCache

  abstract interface
    subroutine abst_add_matrix_to_cache(this, matrix, label)
      import :: TMatrixCache
      import :: z_DNS
      import :: TMatLabel
      class(TMatrixCache) :: this
      type(z_DNS):: matrix
      type(TMatLabel) :: label
    end subroutine

    subroutine abst_retrieve_matrix_from_cache(this, matrix, label)
      import :: TMatrixCache
      import :: z_DNS
      import :: TMatLabel
      class(TMatrixCache) :: this
      type(z_DNS) :: matrix
      type(TMatLabel) :: label
    end subroutine

    function abst_is_cached(this, label) result(val)
      import :: TMatrixCache
      import :: TMatLabel
      !> The bound class
      class(TMatrixCache) :: this
      !> matrix label identifier
      type(TMatLabel) :: label
      logical :: val
    end function

    subroutine abst_destroy_matrix_cache(this)
      import :: TMatrixCache
      class(TMatrixCache) :: this
    end subroutine
  end interface

  !> Linked list to store matrix in memory
  type :: TMatrixCacheEntry
    type(TMatrixCacheEntry), pointer :: next => null()
    type(z_DNS), allocatable :: matrix
    type(TMatLabel) :: mat_label
  end type

  !> Mem derived class
  type, extends(TMatrixCache) :: TMatrixCacheMem
    type(TMatrixCacheEntry), pointer :: first => null()
  contains
    procedure :: add => mem_add
    procedure :: retrieve => mem_retrieve
    procedure :: destroy => mem_destroy
    procedure :: is_cached => mem_is_cached
  end type

  !> Disk derived class
  type, extends(TMatrixCache) :: TMatrixCacheDisk
    character(len=LST) :: scratch_path
  contains
    procedure :: add => disk_add
    procedure :: retrieve => disk_retrieve
    procedure :: destroy => disk_destroy
    procedure :: is_cached => disk_is_cached
  end type

  !> Dummy cache when no caching is performed
  type, extends(TMatrixCache) :: TMatrixCacheDummy
  contains
    procedure :: add => dummy_add
    procedure :: retrieve => dummy_retrieve
    procedure :: destroy => dummy_destroy
    procedure :: is_cached => dummy_is_cached
  end type

contains

  function labeleq(label1, label2) result(var)
     type(TMatLabel), intent(in) :: label1, label2
     logical :: var
     var = .false.
     if (label1%kpoint        == label2%kpoint        .and. &
       & label1%spin          == label2%spin          .and. &
       & label1%row_block     == label2%row_block     .and. &
       & label1%col_block     == label2%col_block     .and. &
       & label1%energy_point  == label2%energy_point ) then
          var = .true.
     end if
  end function labeleq


  !! Definitions for TMatrixCacheMem

  subroutine mem_add(this, matrix, label)
    class(TMatrixCacheMem) :: this
    type(z_DNS):: matrix
    type(TMatLabel) :: label

    type(TMatrixCacheEntry), pointer :: p

    if (.not. associated(this%first)) then
      allocate (this%first)
      p => this%first
    else
      ! Always add at the beginning, because it is faster.
      allocate(p)
      p%next => this%first
      this%first => p
    end if
    p%matrix = matrix
    p%mat_label = label
  end subroutine

  subroutine mem_retrieve(this, matrix, label)
    class(TMatrixCacheMem) :: this
    type(z_DNS) :: matrix
    type(TMatLabel) :: label

    type(TMatrixCacheEntry), pointer :: p

    if (.not. associated(this%first)) then
      error stop "No entry in matrix cache"
    else
      p => this%first
    end if

    do while (associated(p))
      if (p%mat_label .eq. label) then
        matrix = p%matrix
        return
      end if
      p => p%next
    end do

    error stop "Cannot retrieve matrix"

  end subroutine

  function mem_is_cached(this, label) result(val)
    class(TMatrixCacheMem) :: this
    type(TMatLabel) :: label
    logical :: val

    type(TMatrixCacheEntry), pointer :: p

    if (.not. associated(this%first)) then
      val = .false.
      return
    else
      p => this%first
    end if

    do while (associated(p))
      if (p%mat_label .eq. label) then
        val = .true.
        return
      end if
      p => p%next
    end do

    val = .false.

  end function

  subroutine mem_destroy(this)
    class(TMatrixCacheMem) :: this

    type(TMatrixCacheEntry), pointer :: p, previous

    if (.not. associated(this%first)) then
      return
    end if

    p => this%first
    do while (associated(p))
      previous => p
      p => p%next
      deallocate (previous)
    end do

  end subroutine

  !! Definitions for TMatrixCacheDisk

  subroutine disk_add(this, matrix, label)
    class(TMatrixCacheDisk) :: this
    type(z_DNS):: matrix
    type(TMatLabel) :: label

    real(kind=dp) :: dens
    character(2) :: ofcont
    character(1) :: ofspin
    character(10) :: ofkpnt
    character(10) :: ofpnt
    character(LST) :: filename
    integer :: file_unit

    call disk_indices_to_filename(filename, label)
    open (newunit=file_unit, file=trim(this%scratch_path)//filename, form='UNFORMATTED')
    call outmat_c(file_unit, .false., matrix%val, matrix%nrow, matrix%ncol)
    close (file_unit)

  end subroutine

  subroutine disk_retrieve(this, matrix, label)
    class(TMatrixCacheDisk) :: this
    !> Input matrix to be retrieved: It needs to be already allocated with
    !> the correct number of rows and columns.
    type(z_DNS) :: matrix
    type(TMatLabel) :: label

    real(kind=dp) :: dens
    character(LST) :: filename
    logical :: file_exists
    integer :: file_unit

    call disk_indices_to_filename(filename, label)
    inquire (file=trim(this%scratch_path)//filename, EXIST=file_exists)

    if (.not. file_exists) then
      error stop "Cannot retrieve matrix from disk: file not found"
    end if

    open (newunit=file_unit, file=trim(this%scratch_path)//filename, form='UNFORMATTED')
    call inmat_c(file_unit, .false., matrix%val, matrix%nrow, matrix%ncol)
    close (file_unit)

  end subroutine

  function disk_is_cached(this, label) result(val)
    class(TMatrixCacheDisk) :: this
    type(TMatLabel) :: label
    logical :: val

    character(LST) :: filename
    logical :: file_exists
    integer :: file_unit

    call disk_indices_to_filename(filename, label)
    inquire (file=trim(this%scratch_path)//filename, EXIST=file_exists)

    if (file_exists) then
      val = .true.
    else
      val = .false.
    end if

  end function

  subroutine disk_destroy(this)
    class(TMatrixCacheDisk) :: this
    ! Empty. Defined only for interface.
  end subroutine

  subroutine disk_indices_to_filename(filename, label)
    character(LST), intent(out) :: filename
    type(TMatLabel) :: label

    character(2) :: ofrow
    character(2) :: ofcol
    character(1) :: ofspin
    character(10) :: ofkpnt
    character(10) :: ofEpnt

    write (ofrow, '(i2.2)') label%row_block
    write (ofcol, '(i2.2)') label%col_block

    if (label%kpoint .le. 99) write (ofkpnt, '(i2.2)') label%kpoint
    if (label%kpoint .gt. 99) write (ofkpnt, '(i3.3)') label%kpoint
    if (label%kpoint .gt. 999) write (ofkpnt, '(i4.4)') label%kpoint
    if (label%kpoint .gt. 9999) stop 'ERROR: too many k-points (> 9999)'
    if (label%energy_point .le. 999) write (ofEpnt, '(i3.3)') label%energy_point
    if (label%energy_point .gt. 999) write (ofEpnt, '(i4.4)') label%energy_point
    if (label%energy_point .gt. 9999) write (ofEpnt, '(i5.5)') label%energy_point
    if (label%energy_point .gt. 99999) stop 'ERROR: too many contour points (> 99999)'
    if (label%spin .eq. 1) ofspin = 'u'
    if (label%spin .eq. 2) ofspin = 'd'
    filename = 'Mat'//ofspin//ofrow//'_'//ofcol//'_'//trim(ofkpnt)//'_'//trim(ofEpnt)//'.dat'

  end subroutine

  !! Definitions for TMatrixCacheDummy

  subroutine dummy_add(this, matrix, label)
    class(TMatrixCacheDummy) :: this
    type(z_DNS):: matrix
    type(TMatLabel) :: label

    ! Dummy operation

  end subroutine

  subroutine dummy_retrieve(this, matrix, label)
    class(TMatrixCacheDummy) :: this
    !> Input matrix to be retrieved: It needs to be already allocated with
    !> the correct number of rows and columns.
    type(z_DNS) :: matrix
    type(TMatLabel) :: label

    ! Dummy operation

  end subroutine

  function dummy_is_cached(this, label) result(val)
    class(TMatrixCacheDummy) :: this
    type(TMatLabel) :: label
    logical :: val

    val = .false.

  end function

  subroutine dummy_destroy(this)
    class(TMatrixCacheDummy) :: this
  end subroutine
end module
