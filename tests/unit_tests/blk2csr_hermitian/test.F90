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

!> Unit test of blk2csr with a Hermitian matrix of which only the diagonal and lower
!> blocks are stored: the upper part must be rebuilt by conjugation.
program test_blk2csr_hermitian

  use ln_precision, only : dp
  use mat_def, only : z_DNS, z_CSR, create, destroy
  use ln_structure, only : TStruct_Info
  use iterative, only : blk2csr

  implicit none

  integer, parameter :: nOrb = 5
  real(dp), parameter :: tolerance = 1.0e-13_dp

  type(z_DNS), allocatable :: G(:,:)
  type(z_CSR) :: pattern, csr
  type(TStruct_Info) :: struct
  complex(dp) :: full(nOrb,nOrb)
  integer :: ii, jj

  ! Hermitian matrix, two layers with 2 and 3 orbitals
  do jj = 1, nOrb
    do ii = 1, nOrb
      full(ii,jj) = cmplx(0.1_dp * ii * jj, 0.25_dp * (ii - jj), dp)
    end do
  end do

  struct%num_PLs = 2
  struct%mat_PL_start = [1, 3, 6]
  struct%mat_PL_end = [2, 5]

  allocate(G(2,2))
  call make_block(G(1,1), full(1:2,1:2))
  call make_block(G(2,2), full(3:5,3:5))
  call make_block(G(2,1), full(3:5,1:2))

  ! Dense sparsity pattern, sorted by columns in each row
  call create(pattern, nOrb, nOrb, nOrb*nOrb)
  pattern%rowpnt = [(1 + nOrb*(ii-1), ii = 1, nOrb+1)]
  pattern%colind = [((jj, jj = 1, nOrb), ii = 1, nOrb)]
  pattern%nzval = (0.0_dp, 0.0_dp)

  call blk2csr(G, struct, pattern, csr, hermitian=.true.)
  do ii = 1, nOrb
    do jj = 1, nOrb
      call assert(abs(element(csr, ii, jj) - full(ii,jj)) < tolerance, &
          & 'wrong element with hermitian=.true.')
    end do
  end do
  call destroy(csr)

  ! Without the flag the missing upper block is skipped
  call blk2csr(G, struct, pattern, csr)
  do jj = 3, nOrb
    do ii = 1, 2
      call assert(abs(element(csr, ii, jj)) < tolerance, 'upper block must be empty without the flag')
    end do
  end do
  do jj = 1, 2
    do ii = 3, nOrb
      call assert(abs(element(csr, ii, jj) - full(ii,jj)) < tolerance, 'wrong lower block')
    end do
  end do
  call destroy(csr)

  write(*,*) 'Done'

contains

  subroutine make_block(block, values)
    type(z_DNS), intent(inout) :: block
    complex(dp), intent(in) :: values(:,:)

    call create(block, size(values,1), size(values,2))
    block%val = values

  end subroutine make_block


  !> Element (row, col) of a matrix with the dense pattern built above
  function element(matrix, row, col) result(val)
    type(z_CSR), intent(in) :: matrix
    integer, intent(in) :: row, col
    complex(dp) :: val

    val = matrix%nzval(matrix%rowpnt(row) + col - 1)

  end function element


  subroutine assert(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not.condition) then
      write(*,*) 'ERROR: ', message
      error stop 1
    end if

  end subroutine assert

end program test_blk2csr_hermitian
