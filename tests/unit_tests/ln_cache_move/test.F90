!!--------------------------------------------------------------------------!
!! libNEGF: a general library for Non-Equilibrium Green's functions.        !
!! Copyright (C) 2012 - 2026                                                !
!!                                                                          !
!! This file is part of libNEGF: a library for                              !
!! Non Equilibrium Green's Function calculation                             !
!!                                                                          !
!! libNEGF is free software: you can redistribute it and/or modify it       !
!! under the terms of the GNU Lesser General Public License as published    !
!! by the Free Software Foundation, either version 3 of the License, or     !
!! (at your option) any later version.                                      !
!!                                                                          !
!!  You should have received a copy of the GNU Lesser General Public        !
!!  License along with libNEGF.  If not, see                                !
!!  <http://www.gnu.org/licenses/>.                                         !
!!--------------------------------------------------------------------------!

!> Unit test of TMatrixCache%move_dp: the source is left empty, no extra memory is
!> used, the cached block is correct, and the memory tracker stays balanced.
program test_ln_cache_move

  use iso_c_binding, only : c_int64_t
  use ln_precision, only : dp
  use ln_allocation, only : getMem
  use mat_def, only : z_DNS, create, destroy
  use ln_cache

  implicit none

  real(dp), parameter :: tolerance = 1.0e-6_dp

  type(TMatrixCacheMem) :: cache
  type(TMatrixCacheDummy) :: dummyCache
  type(TMatLabel) :: label
  type(z_DNS) :: source
  integer(c_int64_t) :: memStart, memWithSource, mem

  call getMem(memStart)

  cache%tagname = 'test'
  call cache%init(1, 1, 2, 3, 1)
  label%energy_point = 1
  label%kpoint = 1
  label%spin = 1
  label%row_block = 2
  label%col_block = 1

  ! Move into an empty cache entry
  call make_block(source, 1.0_dp)
  call getMem(memWithSource)
  call cache%move_dp(source, label)
  call assert_empty(source, 'moved source')
  call getMem(mem)
  call assert(mem <= memWithSource, 'move must not allocate additional memory')
  call assert_cached(cache, label, 1.0_dp)

  ! Move onto an occupied entry: the old block is replaced and freed
  call make_block(source, 2.0_dp)
  call cache%move_dp(source, label)
  call assert_empty(source, 'source moved onto an occupied entry')
  call getMem(mem)
  call assert(mem <= memWithSource, 'the replaced block must be freed')
  call assert_cached(cache, label, 2.0_dp)

  ! A plain add copies and leaves the source untouched
  call make_block(source, 3.0_dp)
  call cache%add(source, label)
  call assert(allocated(source%val), 'add must keep the source')
  call assert_cached(cache, label, 3.0_dp)
  call destroy(source)

  call cache%destroy()
  call getMem(mem)
  call assert(mem == memStart, 'memory tracker not balanced after destroying the cache')

  ! Default implementation: copy to the cache, then free the source
  call make_block(source, 4.0_dp)
  call dummyCache%move_dp(source, label)
  call assert_empty(source, 'source moved into the dummy cache')
  call getMem(mem)
  call assert(mem == memStart, 'memory tracker not balanced with the default move')

  write(*,*) 'Done'

contains

  !> Rectangular block with values depending on factor
  subroutine make_block(block, factor)
    type(z_DNS), intent(inout) :: block
    real(dp), intent(in) :: factor

    integer :: ii, jj

    call create(block, 3, 2)
    do jj = 1, 2
      do ii = 1, 3
        block%val(ii,jj) = factor * cmplx(ii, jj, dp)
      end do
    end do

  end subroutine make_block


  subroutine assert_empty(block, name)
    type(z_DNS), intent(in) :: block
    character(len=*), intent(in) :: name

    call assert(.not.allocated(block%val), name//' is still allocated')
    call assert(block%nrow == 0 .and. block%ncol == 0, name//' keeps its dimensions')

  end subroutine assert_empty


  subroutine assert_cached(matrixCache, matrixLabel, factor)
    type(TMatrixCacheMem), intent(inout) :: matrixCache
    type(TMatLabel), intent(in) :: matrixLabel
    real(dp), intent(in) :: factor

    type(z_DNS) :: retrieved
    type(z_DNS) :: expected

    call matrixCache%retrieve(retrieved, matrixLabel)
    call make_block(expected, factor)
    call assert(retrieved%nrow == 3 .and. retrieved%ncol == 2, 'wrong shape of the cached block')
    call assert(maxval(abs(retrieved%val - expected%val)) < tolerance, &
        & 'wrong values in the cached block')
    call destroy(retrieved)
    call destroy(expected)

  end subroutine assert_cached


  subroutine assert(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not.condition) then
      write(*,*) 'ERROR: ', message
      error stop 1
    end if

  end subroutine assert

end program test_ln_cache_move
