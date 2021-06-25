module iterative_cpu
  use ln_precision
  use ln_constants, only : pi
  use ln_allocation
  use mat_def
  use sparsekit_drv
  use inversions
  use elph
  use ln_structure, only : TStruct_Info
  use lib_param, only : MAXNCONT, Tnegf, intarray


  implicit none
  private
  public :: calculate_gsmr_blocks
  public :: calculate_gsml_blocks
  public :: calculate_Gr_tridiag_blocks
  public :: calculate_Gn_tridiag_blocks
  public :: calculate_single_transmission_2_contacts
  public :: calculate_single_transmission_N_contacts
  public :: check_convergence_trid

contains

  !***********************************************************************
  !
  !  g_small right (gsmr) calculation - write on memory
  !
  !***********************************************************************

  subroutine calculate_gsmr_blocks(negf,ESH,sbl,ebl,gsmr,keepall)

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
    type(z_DNS), dimension(:), intent(inout) :: gsmr
    type(Tnegf), intent(in) :: negf
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


  !***********************************************************************
  !
  !  g_small left (gsml) calculation - write on memory
  !
  !***********************************************************************

  subroutine calculate_gsml_blocks(negf,ESH,sbl,ebl,gsml)

    !***********************************************************************
    !Input:
    !ESH: sparse matrices array ESH(nbl,nbl)
    !
    !nbl (number of layers), indblk(nbl+1) global variables needed
    !
    !Output:
    !sparse matrices array global variable gsml(nbl) is available in memory
    !single blocks are allocated internally, array Gr(nbl,nbl)
    !must be allocated externally
    !***********************************************************************


    implicit none

    !In/Out
    type(z_DNS), dimension(:), intent(inout) :: gsml
    type(Tnegf), intent(in) :: negf
    type(z_DNS), dimension(:,:), intent(in) :: ESH
    integer, intent(in) :: sbl,ebl                       ! start block, end block

    !Work
    type(z_DNS) :: work1, work2
    integer :: nrow
    integer :: i, nbl
    !type(z_DNS) :: INV(sbl,sbl)

    if (sbl.gt.ebl) return

    nbl = size(ESH,1)

    if (nbl.eq.1) return

    nrow=ESH(sbl,sbl)%nrow

    call create(gsml(sbl),nrow,nrow)

    call compGreen(gsml(sbl),ESH(sbl,sbl),nrow)


    do i=sbl+1,ebl

       nrow=ESH(i,i)%nrow

       call prealloc_mult(ESH(i,i-1),gsml(i-1),(-1.0_dp, 0.0_dp),work1)

       call prealloc_mult(work1,ESH(i-1,i),work2)

       call destroy(work1)

       call prealloc_sum(ESH(i,i),work2,work1)

       call destroy(work2)

       call create(gsml(i),work1%nrow,work1%nrow)

       call compGreen(gsml(i),work1,work1%nrow)

       call destroy(work1)

    end do


  end subroutine calculate_gsml_blocks

  !***********************************************************************
  !
  !  Diagonal, Subdiagonal, Superdiagonal blocks of Green Retarded
  !  Gr(nbl,nbl) - writing on memory
  !
  !***********************************************************************

  subroutine calculate_Gr_tridiag_blocks(negf,ESH,gsml,gsmr,Gr,sbl,ebl)

    !***********************************************************************
    !Input:
    !ESH: dense matrices array ESH(nbl,nbl)
    !sbl, ebl : block indexes
    ! If only sbl is specified, it calculates Gr(sbl, sbl)
    ! If sbl > ebl, it calculates Gr(ebl:sbl, ebl:sbl), Gr(ebl:sbl + 1, ebl:sbl),
    !    Gr(ebl:sbl, ebl:sbl + 1) (need gsml)
    ! If sbl < ebl, it calculates Gr(sbl:ebl, sbl:ebl), Gr(sbl:ebl - 1, sbl:ebl),
    !    Gr(sbl:ebl, sbl:ebl - 1) (need gsmr)
    !
    !
    !Output:
    !sparse matrices array global variable Gr(nbl,nbl) is available in
    !memory - single blocks are allocated internally, array Gr(nbl,nbl)
    !must be allocated externally
    !***********************************************************************

    implicit none

    !In/Out
    type(z_DNS), dimension(:), intent(in) :: gsmr, gsml
    type(z_DNS), dimension(:,:), intent(inout) :: Gr
    type(Tnegf), intent(in) :: negf
    type(z_DNS), dimension(:,:), intent(in) :: ESH
    integer, intent(in) :: sbl
    integer, intent(in), optional :: ebl

    !Work
    integer :: i,nrow,nbl
    type(z_DNS) :: work1, work2, work3

    nbl = size(ESH,1)

    if (sbl.gt.nbl) return
    if (sbl.lt.1) return

    if (.not.present(ebl)) then
       if (nbl.eq.1) then
          nrow = ESH(sbl,sbl)%nrow
          call create(Gr(sbl,sbl),nrow,nrow)
          call compGreen(Gr(sbl,sbl),ESH(sbl,sbl),nrow)
       else
          nrow = ESH(sbl,sbl)%nrow
          call create(work1,nrow,nrow)
          work1%val = ESH(sbl,sbl)%val
          if (sbl+1.le.nbl) then
             call prealloc_mult(ESH(sbl,sbl+1),gsmr(sbl+1),work2)
             call prealloc_mult(work2,ESH(sbl+1,sbl),work3)
             call destroy(work2)
             call prealloc_sum(work1,work3,(-1.0_dp, 0.0_dp),work2)
             call destroy(work3)
             work1%val = work2%val
             call destroy(work2)
          endif
          if (sbl-1.ge.1) then
             call prealloc_mult(ESH(sbl,sbl-1),gsml(sbl-1),work2)
             call prealloc_mult(work2,ESH(sbl-1,sbl),work3)
             call destroy(work2)
             call prealloc_sum(work1,work3,(-1.0_dp, 0.0_dp),work2)
             call destroy(work3)
             work1%val = work2%val
             call destroy(work2)
          endif

          call create(Gr(sbl,sbl),nrow,nrow)
          call compGreen(Gr(sbl,sbl),work1,nrow)
          call destroy(work1)
       endif
       return
    endif


    !***
    !Diagonal, Subdiagonal and Superdiagonal blocks
    !***
    if ((ebl.ge.sbl).and.(ebl.gt.1).and.(sbl.gt.1)) THEN
       do i=sbl,ebl,1
          call prealloc_mult(gsmr(i),ESH(i,i-1),work1)
          call prealloc_mult(work1,Gr(i-1,i-1),(-1.0_dp,0.0_dp),Gr(i,i-1))
          call destroy(work1)

          call prealloc_mult(ESH(i-1,i),gsmr(i),work2)
          call prealloc_mult(Gr(i-1,i-1),work2,(-1.0_dp, 0.0_dp),Gr(i-1,i))

          call prealloc_mult(Gr(i,i-1),work2,(-1.0_dp,0.0_dp),work1)
          call destroy(work2)

          call prealloc_sum(gsmr(i),work1,Gr(i,i))
          call destroy(work1)
       end do
    ELSE
       do i=sbl,ebl,-1
          call prealloc_mult(gsml(i),ESH(i,i+1),work1)
          call prealloc_mult(work1,Gr(i+1,i+1),(-1.0_dp,0.0_dp),Gr(i,i+1))
          call destroy(work1)

          call prealloc_mult(ESH(i+1,i),gsml(i),work2)
          call prealloc_mult(Gr(i+1,i+1),work2,(-1.0_dp, 0.0_dp),Gr(i+1,i))

          call prealloc_mult(Gr(i,i+1),work2,(-1.0_dp,0.0_dp),work1)
          call destroy(work2)

          call prealloc_sum(gsml(i),work1,Gr(i,i))
          call destroy(work1)
       end do
    endif

  end subroutine calculate_Gr_tridiag_blocks

  !****************************************************************************
  !
  ! Calculate G_n contributions for all contacts (except reference)
  !
  !   Sum   [f_j(E)-f_r(E)] Gr Gam_j Ga
  !   j!=r
  !
  ! TRICK: in order to have the usual Gn set ref>ncont and make sure f_r = 0
  !
  !****************************************************************************

  subroutine calculate_Gn_tridiag_blocks(negf,ESH,SelfEneR,frm,ref,struct,gsml,gsmr,Gr,Gn)

    !******************************************************************************
    !Input:
    !ESH(nbl,nbl): sparse matrices array ES-H
    !SelfEneR(ncont): sparse matrices array with contacts Self Energy
    !frm(ncont): Fermi distribution value for each contact
    !ref:  reference contact
    !
    !global variables needed: nbl, indblk(nbl+1), ncont, cblk(ncont), Gr(:,:)
    !diagonal, subdiagonal, overdiagonal and in colums Gr(:,cb) where cb are
    !layers interacting with all contacts but collector
    !
    !Output:
    !Gl: sparse matrix containing G  contacts contribution
    !*******************************************************************************

    implicit none

    !In/Out
    type(z_DNS), dimension(:), intent(in) :: gsmr, gsml
    type(Tnegf), intent(in) :: negf
    type(z_DNS), dimension(:,:), intent(in) :: ESH, Gr
    type(z_DNS), dimension(:,:), intent(inout) :: Gn
    type(z_DNS), dimension(:), intent(in) :: SelfEneR
    type(Tstruct_info), intent(in) :: struct
    real(dp), dimension(:), intent(in) :: frm
    integer, intent(in) :: ref

    !Work
    Type(z_DNS) :: Gam
    type(z_DNS) :: work1,Ga
    integer :: i,j,cb
    integer :: ncont, nbl
    complex(dp) :: frmdiff

    ncont = struct%num_conts
    nbl = struct%num_PLs

    !*******************************************
    ! Contact Iteration
    !*******************************************
    do j=1,ncont

       ! NOTE: this soubroutine uses prealloc_mult that performs
       ! C = C + A*B
       if (j.NE.ref .AND. ABS(frm(j)-frm(ref)).GT.EPS) THEN

          cb=struct%cblk(j) ! block corresponding to contact j

          call zspectral(SelfEneR(j),SelfEneR(j),0,Gam)

          frmdiff = cmplx(frm(j)-frm(ref),0.0_dp,dp)
          ! Computation of Gl(1,1) = Gr(1,cb) Gam(cb) Ga(cb,1)
          if (allocated(Gr(1,cb)%val)) then
             call prealloc_mult(Gr(1,cb),Gam,work1)
             call zdagger(Gr(1,cb),Ga)
             call prealloc_mult(work1,Ga,frmdiff,Gn(1,1))
             call destroy(work1, Ga)
          else
             Gn(1,1)%val=(0.0_dp,0.0_dp)
          endif

          ! Computation of all tridiagonal blocks
          do i=2,nbl

             ! Computation of Gl(i,j) = Gr(i,cb) Gam(cb) Ga(cb,j)
             ! Both Gr(i,cb) and Gr(j,cb) must be non-zero
             if (Gr(i-1,cb)%nrow.gt.0 .and. Gr(i,cb)%nrow.gt.0) then
                call prealloc_mult(Gr(i-1,cb),Gam,work1)
                call zdagger(Gr(i,cb),Ga)
                call prealloc_mult(work1,Ga,frmdiff,Gn(i-1,i))
                call destroy(work1)

                call prealloc_mult(Gr(i,cb),Gam,work1)
                call prealloc_mult(work1,Ga,frmdiff,Gn(i,i))

                call destroy(work1, Ga)

                call prealloc_mult(Gr(i,cb),Gam,work1)
                call zdagger(Gr(i-1,cb),Ga)
                call prealloc_mult(work1,Ga,frmdiff,Gn(i,i-1))

                call destroy(work1, Ga)
             else
                Gn(i  ,i)%val=(0.0_dp,0.0_dp)
                Gn(i-1,i)%val=(0.0_dp,0.0_dp)
                Gn(i,i-1)%val=(0.0_dp,0.0_dp)
             endif

          end do

          call destroy(Gam)

       endif
    end do

  end subroutine calculate_Gn_tridiag_blocks
  
  subroutine check_convergence_trid(negf,T,nome,gpu)
    type(Tnegf), intent(in) :: negf      
    type(z_DNS), dimension(:,:), intent(in) :: T
    character(3), intent(in) :: nome
    logical, intent(in) :: gpu

    integer :: nbl, i
    real(dp) :: summ 

    nbl = size(T,1)
    
    if (gpu .eq. .true.) then
       write(*,*) '~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-'     
    else
       write(*,*) '~-~-~-~-',nome,' check convergence CPU: ~-~-~-~-'     
       write(*,*) '       ',nome,'(',1,1,')=', sum(ABS(T(1,1)%val))

       do i= 2,nbl
          write(*,*) '       ',nome,'(',i,i,')=', sum(ABS(T(i,i)%val))
          
          write(*,*) '       ',nome,'(',i,i-1,')=', sum(ABS(T(i,i-1)%val))
          
          write(*,*) '       ',nome,'(',i-1,i,')=', sum(ABS(T(i-1,i)%val))
       end do
       write(*,*) '~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-'     
    endif     
  end subroutine check_convergence_trid

end module iterative_cpu

