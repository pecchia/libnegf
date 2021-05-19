submodule (iterative) iterative_cpu
use inversions



contains

  !***********************************************************************
  !
  !  g_small right (gsmr) calculation - write on memory
  !
  !***********************************************************************

  module subroutine calculate_gsmr_blocks(ESH,sbl,ebl,keepall)

    !***********************************************************************
    !Input:
    !ESH: sparse matrices array ESH(nbl,nbl)
    !
    !nbl (number of layers)  global needed   (indblk not anymore)
    !
    !Output:
    !sparse matrices array global variable gsmr(nbl) is available in memory
    !single blocks are allocated internally, array Gr(nbl,nbl)
    !must be allocated externally
    !***********************************************************************

    implicit none

    !In/Out
    type(z_DNS), dimension(:,:), intent(in) :: ESH
    integer, intent(in) :: sbl,ebl                 ! start block, end block
    logical, intent(in), optional :: keepall

    !Work
    !type(z_DNS), dimension(:,:), allocatable :: INV
    type(z_DNS) :: work1, work2
    integer :: nrow, M, N
    integer :: i, nbl
    logical :: keep

    if (sbl.lt.ebl) return

    nbl = size(ESH,1)

    if (nbl.eq.1) return

    keep = .true.
    if (present(keepall)) then
      keep = keepall
    end if

    nrow=ESH(sbl,sbl)%nrow

    call create(gsmr(sbl),nrow,nrow)

    call compGreen(gsmr(sbl),ESH(sbl,sbl),nrow)


    do i=sbl-1,ebl,-1

      call prealloc_mult(ESH(i,i+1),gsmr(i+1),(-1.0_dp, 0.0_dp),work1)

      if (.not.keep) then
        call destroy(gsmr(i+1))
      end if

      call prealloc_mult(work1,ESH(i+1,i),work2)

      call destroy(work1)

      call prealloc_sum(ESH(i,i),work2,work1)

      call destroy(work2)

      call create(gsmr(i),work1%nrow,work1%nrow)

      call compGreen(gsmr(i),work1,work1%nrow)

      call destroy(work1)

    end do


  end subroutine calculate_gsmr_blocks

      
end submodule

