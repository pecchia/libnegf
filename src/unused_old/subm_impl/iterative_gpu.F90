submodule(iterative) iterative_gpu
#:if defined("GPU")
   use cublas_v2
  use cusolverDn
#:endif
contains
  
  module subroutine calculate_gsmr_blocks_sp(negf,ESH,sbl,ebl,keepall)

    !In/Out
    type(Tnegf), intent(in) :: negf
    type(c_DNS), dimension(:,:), intent(in) :: ESH
    integer, intent(in) :: sbl,ebl                 ! start block, end block
    logical, intent(in), optional :: keepall

    !Work
    type(CublasHandle) :: hh
    type(CusolverDnHandle) :: hhsol
    complex(sp), parameter :: one = (1.0, 0.0)
    complex(sp), parameter :: mone = (-1.0, 0.0)
    complex(sp), parameter :: zero = (0.0, 0.0)
    type(c_DNS) :: work1, work2
    integer :: nrow, ncol
    integer :: i, nbl, istat
    real(sp) :: sum1
    logical :: keep

    call copy_ESH_toGPU(ESH)

    if (sbl.lt.ebl) return

    nbl = size(ESH,1)
    if (nbl.eq.1) return

    keep = .true.
    if (present(keepall)) then
       keep = keepall
    end if

    nrow=ESH(sbl,sbl)%nrow

    hh = negf%hcublas
    hhsol = negf%hcusolver

    call createAll(gsmr(sbl),nrow,nrow)
    call inverse_gpu(hh, hhsol, ESH(sbl,sbl)%val, gsmr(sbl)%val, istat)

    do i=sbl-1,ebl,-1

       call createAll(work1, ESH(i,i+1)%nrow, gsmr(i+1)%ncol)
       call matmul_gpu(hh, one, ESH(i,i+1)%val, gsmr(i+1)%val, zero, work1%val)

       if (.not.keep) then
          call destroyAll(gsmr(i+1))
       end if

       call createAll(work2, work1%nrow, ESH(i+1,i)%ncol)
       call matmul_gpu(hh, one, work1%val, ESH(i+1,i)%val, zero, work2%val)
       call destroyAll(work1)

       call createAll(work1,  ESH(i,i)%nrow, ESH(i,i)%ncol)
       call add_cublas(hh, work1%val, ESH(i,i)%val, mone, work2%val)
       call destroyAll(work2)

       call createAll(gsmr(i),work1%nrow,work1%ncol)
       call inverse_gpu(hh, hhsol, work1%val, gsmr(i)%val, istat)    
       call destroyAll(work1)

    end do
  end subroutine calculate_gsmr_blocks_sp

  module subroutine calculate_gsmr_blocks_dp(negf,ESH,sbl,ebl,keepall)

    !In/Out
    type(Tnegf), intent(in) :: negf
    type(z_DNS), dimension(:,:), intent(in) :: ESH
    integer, intent(in) :: sbl,ebl                 ! start block, end block
    logical, intent(in), optional :: keepall

    !Work
    type(CublasHandle) :: hh
    type(CusolverDnHandle) :: hhsol
    complex(dp), parameter :: one = (1.0, 0.0)
    complex(dp), parameter :: mone = (-1.0, 0.0)
    complex(dp), parameter :: zero = (0.0, 0.0)
    type(z_DNS) :: work1, work2
    integer :: nrow, ncol
    integer :: i, nbl, istat
    real(dp) :: sum1
    logical :: keep

    call copy_ESH_toGPU(ESH)
    
    if (sbl.lt.ebl) return

    nbl = size(ESH,1)
    if (nbl.eq.1) return
    
    keep = .true.
    if (present(keepall)) then
       keep = keepall
    end if


    nrow=ESH(sbl,sbl)%nrow

    hh = negf%hcublas
    hhsol = negf%hcusolver

    call createAll(gsmr(sbl),nrow,nrow)
    call inverse_gpu(hh, hhsol, ESH(sbl,sbl)%val, gsmr(sbl)%val, istat)

    do i=sbl-1,ebl,-1

       call createAll(work1, ESH(i,i+1)%nrow, gsmr(i+1)%ncol)    
       call matmul_gpu(hh, one, ESH(i,i+1)%val, gsmr(i+1)%val, zero, work1%val)

       if (.not.keep) then
          call destroyAll(gsmr(i+1))
       end if

       call createAll(work2, work1%nrow, ESH(i+1,i)%ncol)
       call matmul_gpu(hh, one, work1%val, ESH(i+1,i)%val, zero, work2%val)
       call destroyAll(work1)

       call createAll(work1,  ESH(i,i)%nrow, ESH(i,i)%ncol)
       call add_cublas(hh, work1%val, ESH(i,i)%val, mone, work2%val)
       call destroyAll(work2)

       call createAll(gsmr(i), work1%nrow, work1%ncol)
       call inverse_gpu(hh, hhsol, work1%val, gsmr(i)%val, istat)    
       call destroyAll(work1)

    end do

  end subroutine calculate_gsmr_blocks_dp

  module subroutine calculate_gsml_blocks_sp(negf,ESH,sbl,ebl)

    !In/Out
    type(Tnegf), intent(in) :: negf
    type(c_DNS), dimension(:,:), intent(in) :: ESH
    integer, intent(in) :: sbl,ebl                 ! start block, end block

    !Work
    type(CublasHandle) :: hh
    type(CusolverDnHandle) :: hhsol
    complex(sp), parameter :: one = (1.0, 0.0)
    complex(sp), parameter :: mone = (-1.0, 0.0)
    complex(sp), parameter :: zero = (0.0, 0.0)
    type(c_DNS) :: work1, work2
    integer :: nrow, ncol
    integer :: i, nbl, istat
    real(sp) :: sum1

    call copy_ESH_toGPU(ESH)

    if (sbl.gt.ebl) return

    nbl = size(ESH,1)
    if (nbl.eq.1) return

    nrow=ESH(sbl,sbl)%nrow

    hh = negf%hcublas
    hhsol = negf%hcusolver

    call createAll(gsml(sbl),nrow,nrow)
    call inverse_gpu(hh, hhsol, ESH(sbl,sbl)%val, gsml(sbl)%val, istat)

    do i=sbl+1,ebl

       call createAll(work1, ESH(i,i-1)%nrow, gsml(i-1)%ncol)    
       call matmul_gpu(hh, one, ESH(i,i-1)%val, gsml(i-1)%val, zero, work1%val)

       call createAll(work2, work1%nrow, ESH(i-1,i)%ncol)
       call matmul_gpu(hh, one, work1%val, ESH(i-1,i)%val, zero, work2%val)
       call destroyAll(work1)

       call createAll(work1,  ESH(i,i)%nrow, ESH(i,i)%ncol)
       call add_cublas(hh, work1%val, ESH(i,i)%val, mone, work2%val)
       call destroyAll(work2)

       call createAll(gsml(i), work1%nrow, work1%ncol)
       call inverse_gpu(hh, hhsol, work1%val, gsml(i)%val, istat)    
       call destroyAll(work1)

    end do

  end subroutine calculate_gsml_blocks_sp

  module subroutine calculate_gsml_blocks_dp(negf,ESH,sbl,ebl)
    !In/Out
    type(Tnegf), intent(in) :: negf
    type(z_DNS), dimension(:,:), intent(in) :: ESH
    integer, intent(in) :: sbl,ebl                 ! start block, end block

    !Work
    type(CublasHandle) :: hh
    type(CusolverDnHandle) :: hhsol
    complex(dp), parameter :: one = (1.0, 0.0)
    complex(dp), parameter :: mone = (-1.0, 0.0)
    complex(dp), parameter :: zero = (0.0, 0.0)
    type(z_DNS) :: work1, work2
    integer :: nrow, ncol
    integer :: i, nbl, istat
    real(dp) :: sum1

    call copy_ESH_toGPU(ESH)

    if (sbl.gt.ebl) return

    nbl = size(ESH,1)
    if (nbl.eq.1) return

    nrow=ESH(sbl,sbl)%nrow

    hh = negf%hcublas
    hhsol = negf%hcusolver

    call createAll(gsml(sbl),nrow,nrow)
    call inverse_gpu(hh, hhsol, ESH(sbl,sbl)%val, gsml(sbl)%val, istat)

    do i=sbl+1,ebl

       call createAll(work1, ESH(i,i-1)%nrow, gsml(i-1)%ncol)    
       call matmul_gpu(hh, one, ESH(i,i-1)%val, gsml(i-1)%val, zero, work1%val)

       call createAll(work2, work1%nrow, ESH(i-1,i)%ncol)
       call matmul_gpu(hh, one, work1%val, ESH(i-1,i)%val, zero, work2%val)
       call destroyAll(work1)

       call createAll(work1,  ESH(i,i)%nrow, ESH(i,i)%ncol)
       call add_cublas(hh, work1%val, ESH(i,i)%val, mone, work2%val)
       call destroyAll(work2)

       call createAll(gsml(i), work1%nrow, work1%ncol)
       call inverse_gpu(hh, hhsol, work1%val, gsml(i)%val, istat)    
       call destroyAll(work1)

    end do

  end subroutine calculate_gsml_blocks_dp

  module subroutine calculate_Gr_tridiag_blocks_sp(negf,ESH,sbl,ebl)
    !In/Out
    type(Tnegf), intent(in) :: negf
    type(c_DNS), dimension(:,:), intent(in) :: ESH
    integer, intent(in) :: sbl
    integer, intent(in), optional :: ebl

    !Work
    type(CublasHandle) :: hh
    type(CusolverDnHandle) ::hhsol
    integer :: i,nrow,nbl
    type(c_DNS), target :: work1, work2, work3
    complex(sp), parameter :: one = (1.0, 0.0)
    complex(sp), parameter :: mone = (-1.0, 0.0)
    complex(sp), parameter :: zero = (0.0, 0.0)
    integer :: istat

    call copy_ESH_toGPU(ESH)

    nbl = size(ESH,1)
    hh = negf%hcublas
    hhsol = negf%hcusolver

    if (sbl.gt.nbl) return
    if (sbl.lt.1) return

    if (.not.present(ebl)) then
       if (nbl.eq.1) then
          call createAll(Gr(sbl,sbl), ESH(sbl,sbl)%nrow, ESH(sbl,sbl)%ncol)
          call inverse_gpu(hh, hhsol, ESH(sbl,sbl)%val, Gr(sbl,sbl)%val, istat)
       else
          call createAll(work1, ESH(sbl,sbl)%nrow, ESH(sbl,sbl)%ncol)
          call copy_mat_gpu(hh, ESH(sbl,sbl)%val, work1%val)

          if (sbl+1.le.nbl) then
             call createAll(work2, ESH(sbl,sbl+1)%nrow, gsmr(sbl+1)%ncol)
             call matmul_gpu(hh, one, ESH(sbl,sbl+1)%val, gsmr(sbl+1)%val, zero, work2%val)

             call createAll(work3, work2%nrow, ESH(sbl+1,sbl)%ncol)
             call matmul_gpu(hh, one, work2%val, ESH(sbl+1,sbl)%val, zero, work3%val)
             call add_cublas(hh, work1%val, ESH(sbl,sbl)%val, mone, work3%val)

             call destroyAll(work2)
             call destroyAll(work3)
          end if
          if (sbl-1.ge.1) then
             call createAll(work2, ESH(sbl,sbl-1)%nrow, gsml(sbl-1)%ncol)
             call matmul_gpu(hh, one, ESH(sbl,sbl-1)%val, gsml(sbl-1)%val, zero, work2%val)

             call createAll(work3, work2%nrow, ESH(sbl-1,sbl)%ncol)
             call matmul_gpu(hh, one, work2%val, ESH(sbl-1,sbl)%val, zero, work3%val)
             call add_cublas(hh, work1%val, ESH(sbl,sbl)%val, mone, work3%val)

             call destroyAll(work2)
             call destroyAll(work3)
          end if

          call createAll(Gr(sbl,sbl), work1%nrow, work1%ncol)
          call inverse_gpu(hh, hhsol, work1%val, Gr(sbl,sbl)%val, istat)
          call destroyAll(work1)
       endif
       return
    endif
    !***
    !Diagonal, Subdiagonal and Superdiagonal blocks
    !***
    if ((ebl.ge.sbl).and.(ebl.gt.1).and.(sbl.gt.1)) THEN
       do i=sbl,ebl,1
          call createAll(work1, gsmr(i)%nrow, ESH(i,i-1)%ncol)
          call createAll(Gr(i,i-1), work1%nrow, Gr(i-1,i-1)%ncol)
          call matmul_gpu(hh, one, gsmr(i)%val, ESH(i,i-1)%val, zero, work1%val)
          call matmul_gpu(hh, mone, work1%val, Gr(i-1,i-1)%val, zero, Gr(i,i-1)%val)

          call destroyAll(work1)
          call createAll(work2, ESH(i-1,i)%nrow, gsmr(i)%ncol)
          call createAll(Gr(i-1,i), Gr(i-1,i-1)%nrow, work2%ncol)
          call matmul_gpu(hh, one, ESH(i-1,i)%val, gsmr(i)%val, zero, work2%val)
          call matmul_gpu(hh, mone, Gr(i-1,i-1)%val, work2%val, zero, Gr(i-1,i)%val)

          call createAll(work1, Gr(i,i-1)%nrow, work2%ncol)
          call matmul_gpu(hh, mone, Gr(i,i-1)%val, work2%val, zero, work1%val)

          call destroyAll(work2)
          call createAll(Gr(i,i), gsmr(i)%nrow, gsmr(i)%ncol)
          call add_cublas(hh, Gr(i,i)%val, gsmr(i)%val, one, work1%val)
          call destroyAll(work1)
       end do
    else
       do i=sbl,ebl,-1
          call createAll(work1, gsml(i)%nrow, ESH(i,i+1)%ncol)
          call createAll(Gr(i,i+1),work1%nrow,Gr(i+1,i+1)%ncol)
          call matmul_gpu(hh, one, gsml(i)%val, ESH(i,i+1)%val, zero, work1%val)
          call matmul_gpu(hh, mone, work1%val, Gr(i+1,i+1)%val, zero, Gr(i,i+1)%val)

          call destroyAll(work1)
          call createAll(work2, ESH(i+1,i)%nrow, gsml(i)%ncol)
          call createAll(Gr(i+1,i), Gr(i+1,i+1)%nrow, work2%ncol)
          call matmul_gpu(hh, one, ESH(i+1,i)%val, gsml(i)%val, zero, work2%val)
          call matmul_gpu(hh, mone, Gr(i+1,i+1)%val, work2%val, zero, Gr(i+1,i)%val)

          call createAll(work1, Gr(i,i+1)%nrow, work2%ncol)
          call matmul_gpu(hh, mone, Gr(i,i+1)%val, work2%val, zero, work1%val)

          call destroyAll(work2)
          call createAll(Gr(i,i), gsml(i)%nrow, gsml(i)%ncol)
          call add_cublas(hh, Gr(i,i)%val, gsml(i)%val, one, work1%val)
          call destroyAll(work1)
       end do
    endif

  end subroutine calculate_Gr_tridiag_blocks_sp

  module subroutine calculate_Gr_tridiag_blocks_dp(negf,ESH,sbl,ebl)
    !In/Out
    type(Tnegf), intent(in) :: negf
    type(z_DNS), dimension(:,:), intent(in) :: ESH
    integer, intent(in) :: sbl
    integer, intent(in), optional :: ebl

    !Work
    type(CublasHandle) :: hh
    type(CusolverDnHandle) :: hhsol
    integer :: i, nrow, nbl, istat
    type(z_DNS), target :: work1, work2, work3
    complex(dp), parameter :: one = (1.0, 0.0)
    complex(dp), parameter :: mone = (-1.0, 0.0)
    complex(dp), parameter :: zero = (0.0, 0.0)

    call copy_ESH_toGPU(ESH)

    nbl = size(ESH,1)
    hh = negf%hcublas
    hhsol = negf%hcusolver

    if (sbl.gt.nbl) return
    if (sbl.lt.1) return

    if (.not.present(ebl)) then
       if (nbl.eq.1) then
          call createAll(Gr(sbl,sbl), ESH(sbl,sbl)%nrow, ESH(sbl,sbl)%ncol)
          call inverse_gpu(hh, hhsol, ESH(sbl,sbl)%val, Gr(sbl,sbl)%val, istat)
       else
          call createAll(work1, ESH(sbl,sbl)%nrow, ESH(sbl,sbl)%ncol)
          call copy_mat_gpu(hh, ESH(sbl,sbl)%val, work1%val)

          if (sbl+1.le.nbl) then
             call createAll(work2, ESH(sbl,sbl+1)%nrow, gsmr(sbl+1)%ncol)
             call matmul_gpu(hh, one, ESH(sbl,sbl+1)%val, gsmr(sbl+1)%val, zero, work2%val)

             call createAll(work3, work2%nrow, ESH(sbl+1,sbl)%ncol)
             call matmul_gpu(hh, one, work2%val, ESH(sbl+1,sbl)%val, zero, work3%val)
             call add_cublas(hh, work1%val, ESH(sbl,sbl)%val, mone, work3%val)

             call destroyAll(work2)
             call destroyAll(work3)
          end if
          if (sbl-1.ge.1) then
             call createAll(work2, ESH(sbl,sbl-1)%nrow, gsml(sbl-1)%ncol)
             call matmul_gpu(hh, one, ESH(sbl,sbl-1)%val, gsml(sbl-1)%val, zero, work2%val)

             call createAll(work3, work2%nrow, ESH(sbl-1,sbl)%ncol)
             call matmul_gpu(hh, one, work2%val, ESH(sbl-1,sbl)%val, zero, work3%val)
             call add_cublas(hh, work1%val, ESH(sbl,sbl)%val, mone, work3%val)

             call destroyAll(work2)
             call destroyAll(work3)
          end if

          call createAll(Gr(sbl,sbl), work1%nrow, work1%ncol)
          call inverse_gpu(hh, hhsol, work1%val, Gr(sbl,sbl)%val, istat)
          call destroyAll(work1)
       endif
       return
    endif
    !***
    !Diagonal, Subdiagonal and Superdiagonal blocks
    !***
    if ((ebl.ge.sbl).and.(ebl.gt.1).and.(sbl.gt.1)) THEN
       do i=sbl,ebl,1
          call createAll(work1, gsmr(i)%nrow, ESH(i,i-1)%ncol)
          call createAll(Gr(i,i-1), work1%nrow, Gr(i-1,i-1)%ncol)
          call matmul_gpu(hh, one, gsmr(i)%val, ESH(i,i-1)%val, zero, work1%val)
          call matmul_gpu(hh, mone, work1%val, Gr(i-1,i-1)%val, zero, Gr(i,i-1)%val)

          call destroyAll(work1)
          call createAll(work2, ESH(i-1,i)%nrow, gsmr(i)%ncol)
          call createAll(Gr(i-1,i), Gr(i-1,i-1)%nrow, work2%ncol)
          call matmul_gpu(hh, one, ESH(i-1,i)%val, gsmr(i)%val, zero, work2%val)
          call matmul_gpu(hh, mone, Gr(i-1,i-1)%val, work2%val, zero, Gr(i-1,i)%val)

          call createAll(work1, Gr(i,i-1)%nrow, work2%ncol)
          call matmul_gpu(hh, mone, Gr(i,i-1)%val, work2%val, zero, work1%val)

          call destroyAll(work2)
          call createAll(Gr(i,i), gsmr(i)%nrow, gsmr(i)%ncol)
          call add_cublas(hh, Gr(i,i)%val, gsmr(i)%val, one, work1%val)
          call destroyAll(work1)
       end do
    else
       do i=sbl,ebl,-1
          call createAll(work1, gsml(i)%nrow, ESH(i,i+1)%ncol)
          call createAll(Gr(i,i+1),work1%nrow,Gr(i+1,i+1)%ncol)
          call matmul_gpu(hh, one, gsml(i)%val, ESH(i,i+1)%val, zero, work1%val)
          call matmul_gpu(hh, mone, work1%val, Gr(i+1,i+1)%val, zero, Gr(i,i+1)%val)

          call destroyAll(work1)
          call createAll(work2, ESH(i+1,i)%nrow, gsml(i)%ncol)
          call createAll(Gr(i+1,i), Gr(i+1,i+1)%nrow, work2%ncol)
          call matmul_gpu(hh, one, ESH(i+1,i)%val, gsml(i)%val, zero, work2%val)
          call matmul_gpu(hh, mone, Gr(i+1,i+1)%val, work2%val, zero, Gr(i+1,i)%val)

          call createAll(work1, Gr(i,i+1)%nrow, work2%ncol)
          call matmul_gpu(hh, mone, Gr(i,i+1)%val, work2%val, zero, work1%val)

          call destroyAll(work2)
          call createAll(Gr(i,i), gsml(i)%nrow, gsml(i)%ncol)
          call add_cublas(hh, Gr(i,i)%val, gsml(i)%val, one, work1%val)
          call destroyAll(work1)
       end do
    endif

  end subroutine calculate_Gr_tridiag_blocks_dp

  module subroutine calculate_Gn_tridiag_blocks_sp(negf, ESH, SelfEneR, frm, ref, struct, Gn) 
    !In/Out
    type(Tnegf), intent(in) :: negf
    type(c_DNS), dimension(:,:), intent(in) :: ESH
    type(c_DNS), dimension(:,:), intent(inout) :: Gn
    type(c_DNS), dimension(:), intent(in) :: SelfEneR
    type(Tstruct_info), intent(in) :: struct
    real(sp), dimension(:), intent(in) :: frm
    integer, intent(in) :: ref

    !Work
    type(CublasHandle) :: hh
    Type(c_DNS) :: Gam
    type(c_DNS) :: work1, work2, work3
    integer :: i,j
    integer :: cb, nbl, ncont, istat
    complex(sp) :: frmdiff
    complex(sp), parameter :: one = (1.0, 0.0)
    complex(sp), parameter :: mone = (-1.0, 0.0)
    complex(sp), parameter :: zero = (0.0, 0.0)
    real(sp) :: somma

    call copy_ESH_toGPU(ESH)

    ncont = struct%num_conts
    nbl = struct%num_PLs
    hh = negf%hcublas


    do j=1,ncont

       if (j.NE.ref .AND. ABS(frm(j)-frm(ref)).GT.EPS) THEN

          cb=struct%cblk(j) ! block corresponding to contact j    

          !Creating Gam = i (SE -  SE^+)
          call createAll(Gam, SelfEneR(j)%nrow, SelfEneR(j)%ncol)
          call zspectral_gpu(SelfEneR(j),Gam)

          frmdiff = cmplx(frm(j)-frm(ref),0.0_sp,sp)

          ! Computation of Gn(cb,cb) = Gr(cb,cb) Gam(cb) Gr(cb,cb)^+
          call createAll(work1, Gam%nrow, Gr(cb,cb)%nrow) !it's Gr^+ --> I take nrow
          call matmul_gpu_dag(hh, frmdiff, Gam%val, Gr(cb,cb)%val, zero, work1%val)
          call matmul_gpu(hh, one, Gr(cb,cb)%val, work1%val, zero, Gn(cb,cb)%val)
          call destroyAll(work1)


          !***************************************************************
          !*** Gr(cb-1, cb) already allocated: must insert if in loop ****
          !***************************************************************

          do i=cb-1, 1, -1
             !Gr(i,cb) = - gL(i) ESH(i,i+1) Gr(i+1,cb)
             if (i .eq. cb-1) then
                call createAll(work1, ESH(i,i+1)%nrow, Gr(i+1,cb)%ncol)
             else 
                call createAll(work1, ESH(i,i+1)%nrow, Gr(i+1,cb)%ncol)   
                call createAll(Gr(i,cb), gsml(i)%nrow, work1%ncol)
             end if
             call matmul_gpu(hh, frmdiff, ESH(i,i+1)%val, Gr(i+1,cb)%val,zero, work1%val)
             call matmul_gpu(hh, mone, gsml(i)%val, work1%val, zero, Gr(i,cb)%val)
             call destroyAll(work1)

             !Gn(i, i) = Gr(i, cb) Gam(cb) Gr(i, cb)^+
             call createAll(work2, Gam%nrow, Gr(i,cb)%nrow)
             call matmul_gpu_dag(hh, frmdiff, Gam%val, Gr(i,cb)%val, zero, work2%val)
             call matmul_gpu(hh, one, Gr(i,cb)%val, work2%val, zero, Gn(i,i)%val)

             !Gn(i+1,i)  = Gr(i+1, cb) Gam(cb) Gr(i, cb)^+
             call matmul_gpu(hh, one, Gr(i+1,cb)%val, work2%val, zero, Gn(i+1,i)%val)
             call destroyAll(work2)

             !Gn(i,i+1)  = Gr(i, cb) Gam(cb) Gr(i+1, cb)^+
             call createAll(work3, Gam%nrow, Gr(i+1,cb)%nrow)
             call matmul_gpu_dag(hh, frmdiff, Gam%val, Gr(i+1,cb)%val, zero, work3%val)
             call matmul_gpu(hh, one, Gr(i,cb)%val, work3%val, zero, Gn(i,i+1)%val)
             call destroyAll(work3)
          end do

          !***************************************************************
          !*** Gr(cb+1, cb) already allocated: must insert if in loop ****
          !***************************************************************

          do i=cb+1, nbl
             !Gr(i,cb) = - gR(i) ESH(i,i-1) Gr(i-1,cb)
             if (i .eq. cb+1) then
                call createAll(work1, ESH(i,i-1)%nrow, Gr(i-1,cb)%ncol)

             else 
                call createAll(work1, ESH(i,i-1)%nrow, Gr(i-1,cb)%ncol)
                call createAll(Gr(i,cb), gsmr(i)%nrow, work1%ncol)
             end if
             call matmul_gpu(hh, frmdiff, ESH(i,i-1)%val, Gr(i-1,cb)%val,zero, work1%val)
             call matmul_gpu(hh, mone, gsmr(i)%val, work1%val, zero, Gr(i,cb)%val)
             call destroyAll(work1)


             !Gn(i, i) = Gr(i, cb) Gam(cb) Gr(i, cb)^+
             call createAll(work2, Gam%nrow, Gr(i,cb)%nrow)
             call matmul_gpu_dag(hh, one, Gam%val, Gr(i,cb)%val, zero, work2%val)
             call matmul_gpu(hh, frmdiff, Gr(i,cb)%val, work2%val, zero, Gn(i,i)%val)

             !Gn(i-1,i)  = Gr(i-1, cb) Gam(cb) Gr(i, cb)^+
             call matmul_gpu(hh, one, Gr(i-1,cb)%val, work2%val, zero, Gn(i-1,i)%val)
             call destroyAll(work2)

             !Gn(i,i-1)  = Gr(i, cb) Gam(cb) Gr(i-1, cb)^+
             call createAll(work3, Gam%nrow, Gr(i-1,cb)%nrow)
             call matmul_gpu_dag(hh, one, Gam%val, Gr(i-1,cb)%val, zero, work3%val)
             call matmul_gpu(hh, frmdiff, Gr(i,cb)%val, work3%val, zero, Gn(i,i-1)%val)
             call destroyAll(work3)
          end do

          do i=cb-2, 1, -1
             call destroyAll(Gr(i,cb))
          end do

          do i=cb+2, nbl
             call destroyAll(Gr(i,cb))
          end do
          call destroy(Gam)
       endif
    end do

  end subroutine calculate_Gn_tridiag_blocks_sp

  module subroutine calculate_Gn_tridiag_blocks_dp(negf, ESH, SelfEneR, frm, ref, struct, Gn)
    !In/Out
    type(Tnegf) :: negf
    type(z_DNS), dimension(:,:), intent(in) :: ESH
    type(z_DNS), dimension(:,:), intent(inout) :: Gn
    type(z_DNS), dimension(:), intent(in) :: SelfEneR
    type(Tstruct_info), intent(in) :: struct
    real(dp), dimension(:), intent(in) :: frm
    integer, intent(in) :: ref

    !Work
    type(CublasHandle) :: hh
    Type(z_DNS) :: Gam
    type(z_DNS) :: work1, work2, work3
    integer :: i, j, istat
    integer :: cb, nbl, ncont
    complex(dp) :: frmdiff
    complex(dp), parameter :: one = (1.0, 0.0)
    complex(dp), parameter :: mone = (-1.0, 0.0)
    complex(dp), parameter :: zero = (0.0, 0.0)
    real(dp) :: somma

    call copy_ESH_toGPU(ESH)

    ncont = struct%num_conts
    nbl = struct%num_PLs
    hh = negf%hcublas


    do j=1,ncont

       if (j.NE.ref .AND. ABS(frm(j)-frm(ref)).GT.EPS) THEN

          cb=struct%cblk(j) ! block corresponding to contact j    

          !Creating Gam = i (SE -  SE^+)
          call createAll(Gam, SelfEneR(j)%nrow, SelfEneR(j)%ncol)
          call zspectral_gpu(SelfEneR(j),Gam)

          frmdiff = cmplx(frm(j)-frm(ref),0.0_sp,sp)

          ! Computation of Gn(cb,cb) = Gr(cb,cb) Gam(cb) Gr(cb,cb)^+
          call createAll(work1, Gam%nrow, Gr(cb,cb)%nrow) !it's Gr^+ --> I take nrow
          call matmul_gpu_dag(hh, frmdiff, Gam%val, Gr(cb,cb)%val, zero, work1%val)
          call matmul_gpu(hh, one, Gr(cb,cb)%val, work1%val, zero, Gn(cb,cb)%val)
          call destroyAll(work1)


          !***************************************************************
          !*** Gr(cb-1, cb) already allocated: must insert if in loop ****
          !***************************************************************

          do i=cb-1, 1, -1
             !Gr(i,cb) = - gL(i) ESH(i,i+1) Gr(i+1,cb)
             if (i .eq. cb-1) then
                call createAll(work1, ESH(i,i+1)%nrow, Gr(i+1,cb)%ncol)
             else 
                call createAll(work1, ESH(i,i+1)%nrow, Gr(i+1,cb)%ncol)   
                call createAll(Gr(i,cb), gsml(i)%nrow, work1%ncol)
             end if
             call matmul_gpu(hh, frmdiff, ESH(i,i+1)%val, Gr(i+1,cb)%val,zero, work1%val)
             call matmul_gpu(hh, mone, gsml(i)%val, work1%val, zero, Gr(i,cb)%val)
             call destroyAll(work1)

             !Gn(i, i) = Gr(i, cb) Gam(cb) Gr(i, cb)^+
             call createAll(work2, Gam%nrow, Gr(i,cb)%nrow)
             call matmul_gpu_dag(hh, frmdiff, Gam%val, Gr(i,cb)%val, zero, work2%val)
             call matmul_gpu(hh, one, Gr(i,cb)%val, work2%val, zero, Gn(i,i)%val)

             !Gn(i+1,i)  = Gr(i+1, cb) Gam(cb) Gr(i, cb)^+
             call matmul_gpu(hh, one, Gr(i+1,cb)%val, work2%val, zero, Gn(i+1,i)%val)
             call destroyAll(work2)

             !Gn(i,i+1)  = Gr(i, cb) Gam(cb) Gr(i+1, cb)^+
             call createAll(work3, Gam%nrow, Gr(i+1,cb)%nrow)
             call matmul_gpu_dag(hh, frmdiff, Gam%val, Gr(i+1,cb)%val, zero, work3%val)
             call matmul_gpu(hh, one, Gr(i,cb)%val, work3%val, zero, Gn(i,i+1)%val)
             call destroyAll(work3)
          end do

          !***************************************************************
          !*** Gr(cb+1, cb) already allocated: must insert if in loop ****
          !***************************************************************

          do i=cb+1, nbl
             !Gr(i,cb) = - gR(i) ESH(i,i-1) Gr(i-1,cb)
             if (i .eq. cb+1) then
                call createAll(work1, ESH(i,i-1)%nrow, Gr(i-1,cb)%ncol)

             else 
                call createAll(work1, ESH(i,i-1)%nrow, Gr(i-1,cb)%ncol)
                call createAll(Gr(i,cb), gsmr(i)%nrow, work1%ncol)
             end if
             call matmul_gpu(hh, frmdiff, ESH(i,i-1)%val, Gr(i-1,cb)%val,zero, work1%val)
             call matmul_gpu(hh, mone, gsmr(i)%val, work1%val, zero, Gr(i,cb)%val)
             call destroyAll(work1)


             !Gn(i, i) = Gr(i, cb) Gam(cb) Gr(i, cb)^+
             call createAll(work2, Gam%nrow, Gr(i,cb)%nrow)
             call matmul_gpu_dag(hh, one, Gam%val, Gr(i,cb)%val, zero, work2%val)
             call matmul_gpu(hh, frmdiff, Gr(i,cb)%val, work2%val, zero, Gn(i,i)%val)

             !Gn(i-1,i)  = Gr(i-1, cb) Gam(cb) Gr(i, cb)^+
             call matmul_gpu(hh, one, Gr(i-1,cb)%val, work2%val, zero, Gn(i-1,i)%val)
             call destroyAll(work2)

             !Gn(i,i-1)  = Gr(i, cb) Gam(cb) Gr(i-1, cb)^+
             call createAll(work3, Gam%nrow, Gr(i-1,cb)%nrow)
             call matmul_gpu_dag(hh, one, Gam%val, Gr(i-1,cb)%val, zero, work3%val)
             call matmul_gpu(hh, frmdiff, Gr(i,cb)%val, work3%val, zero, Gn(i,i-1)%val)
             call destroyAll(work3)
          end do

          do i=cb-2, 1, -1
             call destroyAll(Gr(i,cb))
          end do

          do i=cb+2, nbl
             call destroyAll(Gr(i,cb))
          end do
          call destroy(Gam)
       endif
    end do

  end subroutine calculate_Gn_tridiag_blocks_dp
  
end submodule iterative_gpu
